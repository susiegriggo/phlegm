#!/usr/bin/env python
"""
PHLEGM: predict homo-oligomeric states of proteins from their sequences.

Commands
  phlegm msa      Build one MSA per protein                 (CPU; local DB or ColabFold server)
  phlegm fold     Predict 2-mer..10-mer complexes           (GPU; needs only AlphaFold weights)
  phlegm predict  Call oligomeric states from mean ipTM     (anywhere; seconds)
  phlegm run      msa + fold + predict on one machine
  phlegm download-weights   Fetch the AlphaFold-Multimer weights (~4 GB) ahead of time,
                            e.g. on an HPC login node when GPU nodes have no internet

All commands share one output directory and skip work that is already done,
so any command can be re-run to resume.

Examples
  phlegm run -i proteins.fasta -o out
  phlegm msa -i proteins.fasta -o out --db /path/to/colabfold_db    # CPU node
  phlegm fold -o out --shard 3/20                                   # GPU job-array task
  phlegm predict -o out
"""

import argparse
import csv
import logging
import os
import pickle
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from predict_oligomeric_state import predict_from_directory  # noqa: E402

__version__ = "1.0.1"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYXBZUO")
REGISTRY = Path("input") / "proteins.tsv"

log = logging.getLogger("phlegm")


class PhlegmError(Exception):
    """An error with a message meant for the user (no traceback)."""


# --------------------------------------------------------------------------- #
# Input handling
# --------------------------------------------------------------------------- #

def read_fasta(path):
    """Return a list of (header, sequence) tuples."""
    records, header, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header, seq = line[1:].strip(), []
            else:
                if header is None:
                    raise PhlegmError(f"{path}: sequence data found before the first '>' header")
                seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    if not records:
        raise PhlegmError(f"{path}: no sequences found")
    return records


def safe_id(header):
    """First word of the header, restricted to characters safe for file names.

    '.' is replaced as well because PHLEGM uses '<id>.<subunits>' to name jobs.
    """
    first = header.split()[0] if header.split() else "protein"
    return re.sub(r"[^A-Za-z0-9_-]", "_", first)


def load_fasta_proteins(path):
    proteins, seen = [], {}
    for header, seq in read_fasta(path):
        pid = safe_id(header)
        seq = seq.upper().replace(" ", "").rstrip("*")
        bad = set(seq) - VALID_AA
        if not seq:
            raise PhlegmError(f"'{header}' has an empty sequence")
        if bad:
            raise PhlegmError(f"'{header}' contains invalid characters: {''.join(sorted(bad))}")
        if pid in seen:
            raise PhlegmError(
                f"Headers '{seen[pid]}' and '{header}' both map to the ID '{pid}'. "
                "Please make the first word of each FASTA header unique.")
        seen[pid] = header
        proteins.append({"id": pid, "header": header, "sequence": seq})
    return proteins


def save_registry(proteins, outdir):
    """Record the run's proteins so fold/predict need only the output directory."""
    path = outdir / REGISTRY
    if path.is_file():
        old = load_registry(outdir)
        if [(p["id"], p["sequence"]) for p in old] != [(p["id"], p["sequence"]) for p in proteins]:
            raise PhlegmError(
                f"{outdir} already contains a run for a different set of sequences. "
                "Use a new output directory.")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein_id", "sequence", "header"])
        for p in proteins:
            w.writerow([p["id"], p["sequence"], p["header"]])
    write_fasta(proteins, outdir / "input" / "queries.fasta")


def load_registry(outdir):
    path = outdir / REGISTRY
    if not path.is_file():
        raise PhlegmError(f"{path} not found. Run 'phlegm msa' (or 'phlegm run') "
                          "with this output directory first.")
    with open(path, newline="") as fh:
        return [{"id": r["protein_id"], "sequence": r["sequence"], "header": r["header"]}
                for r in csv.DictReader(fh, delimiter="\t")]


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def run_cmd(cmd, dry_run=False):
    log.info("Running: %s", " ".join(shlex.quote(str(c)) for c in cmd))
    if dry_run:
        return
    # Unbuffered, so ColabFold's progress shows up in SLURM/PBS logs as it happens.
    env = {**os.environ, "PYTHONUNBUFFERED": "1"}
    result = subprocess.run([str(c) for c in cmd], env=env)
    if result.returncode != 0:
        raise PhlegmError(f"Command failed with exit code {result.returncode}: {cmd[0]}")


def require(executable):
    if shutil.which(executable) is None:
        raise PhlegmError(f"'{executable}' was not found on PATH. Did you run 'bash install.sh' "
                          "and 'conda activate phlegm'?")


def write_fasta(proteins, path):
    with open(path, "w") as fh:
        for p in proteins:
            fh.write(f">{p['id']}\n{p['sequence']}\n")


def read_a3m(path):
    """Parse an A3M file into (header, sequence) pairs.

    Removes ColabFold '#len<TAB>cardinality' lines and NUL separators, joins
    wrapped sequence lines and collapses whitespace in headers, so the output
    can be written back with exactly one header line and one sequence line per
    entry (required by ColabFold's A3M parser for complexes).
    """
    entries, header, seq = [], None, []
    with open(path) as fh:
        text = fh.read().replace("\x00", "")
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq)))
            header, seq = re.sub(r"\s+", " ", line), []
        else:
            seq.append(line)
    if header is not None:
        entries.append((header, "".join(seq)))
    return entries


def parse_shard(text):
    m = re.fullmatch(r"(\d+)/(\d+)", text or "")
    if not m or not 1 <= int(m.group(1)) <= int(m.group(2)):
        raise argparse.ArgumentTypeError("expected I/N with 1 <= I <= N, e.g. 3/20")
    return int(m.group(1)), int(m.group(2))


def check_gpu(allow_cpu):
    """Fail early if JAX would silently run AlphaFold on the CPU."""
    probe = subprocess.run(
        [sys.executable, "-c", "import jax; print(jax.devices()[0].platform)"],
        capture_output=True, text=True)
    platform = probe.stdout.strip().splitlines()[-1] if probe.returncode == 0 and probe.stdout.strip() else None
    if platform is None:
        log.warning("Could not check for a GPU (JAX not importable in this Python); continuing")
    elif platform == "cpu":
        msg = ("JAX cannot see a GPU, so AlphaFold would run on the CPU (very slow). "
               "Run 'phlegm fold' on a GPU node (on Linux, 'bash install.sh' without --cpu installs GPU JAX).")
        if not allow_cpu:
            raise PhlegmError(msg + " Use --allow-cpu to run anyway.")
        log.warning(msg)
    else:
        log.info("JAX is using: %s", platform)


# --------------------------------------------------------------------------- #
# AlphaFold weights
# --------------------------------------------------------------------------- #

MODEL_TYPE = "alphafold2_multimer_v3"
WEIGHTS_URL = "https://storage.googleapis.com/alphafold/alphafold_params_colab_2022-12-06.tar"
WEIGHTS_MARKER = Path("params") / "download_complexes_multimer_v3_finished.txt"


def default_weights_dir():
    """Where weights live: $PHLEGM_WEIGHTS_DIR, else ColabFold's default cache."""
    if os.environ.get("PHLEGM_WEIGHTS_DIR"):
        return Path(os.environ["PHLEGM_WEIGHTS_DIR"])
    try:
        from colabfold.download import default_data_dir
        return Path(default_data_dir)
    except Exception:
        if sys.platform == "darwin":
            return Path.home() / "Library" / "Caches" / "colabfold"
        return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "colabfold"


def weights_present(weights_dir):
    return (Path(weights_dir) / WEIGHTS_MARKER).is_file()


def download_weights(weights_dir):
    """Download the AlphaFold-Multimer v3 weights with ColabFold's own downloader."""
    import urllib.request
    try:
        urllib.request.urlopen(urllib.request.Request(WEIGHTS_URL, method="HEAD"), timeout=15)
    except Exception as e:
        raise PhlegmError(
            f"Cannot reach {WEIGHTS_URL} ({e}). This machine seems to have no internet access. "
            "Run 'phlegm download-weights' on a machine that does (e.g. an HPC login node), "
            "using the same --weights-dir.")
    try:
        from colabfold.download import download_alphafold_params
    except ImportError:
        raise PhlegmError("ColabFold is not installed in this environment; see the README.")
    log.info("Downloading AlphaFold-Multimer v3 weights (~4 GB) to %s", weights_dir)
    download_alphafold_params(MODEL_TYPE, Path(weights_dir))
    if not weights_present(weights_dir):
        raise PhlegmError(f"Weight download to {weights_dir} did not complete")
    log.info("Weights ready in %s", weights_dir)


def ensure_weights(args):
    """Make sure the weights exist before any GPU work starts."""
    wdir = Path(args.weights_dir)
    if weights_present(wdir):
        log.info("AlphaFold weights: %s", wdir)
        return
    if getattr(args, "shard", None):
        # Parallel shards must not all download into the same directory at once.
        raise PhlegmError(
            f"AlphaFold weights not found in {wdir}. Before submitting sharded fold jobs, "
            "run 'phlegm download-weights' once (with the same --weights-dir, if you set one).")
    log.info("AlphaFold weights not found in %s; downloading them now (one-off)", wdir)
    download_weights(wdir)


# --------------------------------------------------------------------------- #
# Step 1: monomer MSAs
# --------------------------------------------------------------------------- #

def step_msa(proteins, args, outdir):
    msa_dir = outdir / "msas"
    msa_dir.mkdir(parents=True, exist_ok=True)

    todo = [p for p in proteins if not (msa_dir / f"{p['id']}.a3m").is_file()]
    if not todo:
        log.info("MSA: all %d MSAs already present", len(proteins))
        return
    log.info("MSA: generating MSAs for %d protein(s)", len(todo))

    if args.msa_dir:
        for p in todo:
            src = Path(args.msa_dir) / f"{p['id']}.a3m"
            if not src.is_file():
                raise PhlegmError(f"--msa-dir: expected {src}")
            shutil.copy(src, msa_dir / src.name)
        return

    work = msa_dir / "_tmp"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir()
    query = work / "queries.fasta"
    write_fasta(todo, query)

    if args.db:
        require("colabfold_search")
        cmd = ["colabfold_search", query, args.db, work, "--threads", args.threads]
        if args.db_load_mode is not None:
            cmd += ["--db-load-mode", args.db_load_mode]
        run_cmd(cmd, args.dry_run)
    else:
        require("colabfold_batch")
        log.info("Using the public ColabFold MMseqs2 server for MSAs "
                 "(use --db for a local database on large datasets)")
        run_cmd(["colabfold_batch", query, work, "--msa-only"], args.dry_run)

    if args.dry_run:
        return

    # colabfold_search names monomer MSAs by query index (0.a3m, 1.a3m, ...);
    # colabfold_batch --msa-only names them by query ID. Accept either.
    for i, p in enumerate(todo):
        candidates = [work / f"{p['id']}.a3m", work / f"{i}.a3m"]
        found = next((c for c in candidates if c.is_file()), None)
        if found is None:
            raise PhlegmError(f"No MSA was produced for {p['id']} in {work}")
        entries = read_a3m(found)
        if not entries or entries[0][1].replace("-", "").upper() != p["sequence"]:
            raise PhlegmError(f"First sequence of {found} does not match the query {p['id']}")
        shutil.move(str(found), msa_dir / f"{p['id']}.a3m")
    shutil.rmtree(work)


# --------------------------------------------------------------------------- #
# Step 2: structure prediction
# --------------------------------------------------------------------------- #

def write_oligomer_msa(protein, monomer_a3m, n_subunits, out_path):
    """Reuse a monomer MSA as a homo-oligomer input for colabfold_batch.

    The '#<length>\\t<n>' header tells ColabFold to model n copies of the query.
    """
    seq = protein["sequence"]
    tmp = out_path.with_suffix(".tmp")
    with open(tmp, "w") as fh:
        fh.write(f"#{len(seq)}\t{n_subunits}\n")
        fh.write(f">{protein['id']}\n{seq}\n")
        for header, s in read_a3m(monomer_a3m):
            fh.write(f"{header}\n{s}\n")
    tmp.replace(out_path)


def check_fold_complete(protein, args, outdir):
    """colabfold_batch exits successfully even when a job fails (e.g. out of GPU
    memory), so check that every stoichiometry left its '.done.txt' marker."""
    sdir = outdir / "structures" / protein["id"]
    missing = [n for n in range(args.min_oligomer, args.max_oligomer + 1)
               if not (sdir / f"{protein['id']}.{n}.done.txt").is_file()]
    if missing:
        raise PhlegmError(
            f"ColabFold did not finish the {', '.join(f'{n}-mer' for n in missing)} "
            f"(often because the GPU ran out of memory; see {sdir / 'log.txt'}). "
            "Re-running 'phlegm fold' retries only the missing stoichiometries, "
            "e.g. on a GPU with more memory.")


def step_fold(proteins, args, outdir):
    if getattr(args, "shard", None):
        i, n = args.shard
        proteins = proteins[i - 1::n]
        log.info("Fold: shard %d/%d has %d protein(s)", i, n, len(proteins))
        if not proteins:
            return []

    missing = [p["id"] for p in proteins if not (outdir / "msas" / f"{p['id']}.a3m").is_file()]
    if missing:
        raise PhlegmError(f"Missing MSAs for {len(missing)} protein(s) "
                          f"({', '.join(missing[:5])}{', ...' if len(missing) > 5 else ''}). "
                          "Run 'phlegm msa' first.")

    require("colabfold_batch")
    if not args.dry_run:
        check_gpu(args.allow_cpu)
        ensure_weights(args)

    log.info("Fold: %d-mer to %d-mer with AlphaFold-Multimer v3", args.min_oligomer, args.max_oligomer)
    failed = []
    for k, p in enumerate(proteins, 1):
        pdir = outdir / "oligomer_msas" / p["id"]
        pdir.mkdir(parents=True, exist_ok=True)
        for n in range(args.min_oligomer, args.max_oligomer + 1):
            out = pdir / f"{p['id']}.{n}.a3m"
            if not out.is_file():
                write_oligomer_msa(p, outdir / "msas" / f"{p['id']}.a3m", n, out)

        log.info("[%d/%d] %s (%d aa)", k, len(proteins), p["id"], len(p["sequence"]))
        cmd = ["colabfold_batch", pdir, outdir / "structures" / p["id"],
               "--model-type", "alphafold2_multimer_v3",
               "--num-models", args.num_models,
               "--num-recycle", args.num_recycle,
               "--data", args.weights_dir]
        cmd += shlex.split(args.colabfold_args or "")
        try:
            run_cmd(cmd, args.dry_run)
            if not args.dry_run:
                check_fold_complete(p, args, outdir)
        except PhlegmError as e:
            if not args.keep_going:
                raise
            log.error("%s: %s (continuing because of --keep-going)", p["id"], e)
            failed.append(p["id"])
    if failed:
        log.error("Fold failed for %d protein(s): %s", len(failed), ", ".join(failed))
    return failed


# --------------------------------------------------------------------------- #
# Step 3: oligomeric state prediction
# --------------------------------------------------------------------------- #

def step_predict(proteins, args, outdir):
    log.info("Predict: mean ipTM threshold %.2f", args.threshold)
    subunits = list(range(args.min_oligomer, args.max_oligomer + 1))
    fields = (["protein_id", "header", "length", "predicted_subunits", "predicted_state",
               "max_mean_iptm", "top_multimer_subunits"]
              + [f"mean_iptm_{n}mer" for n in subunits] + ["status"])
    all_results, rows = {}, []

    for p in proteins:
        row = {"protein_id": p["id"], "header": p["header"], "length": len(p["sequence"])}
        sdir = outdir / "structures" / p["id"]
        if not sdir.is_dir():
            row["status"] = "not folded"
            rows.append(row)
            continue
        try:
            _, res = predict_from_directory(str(sdir), p["id"], args.threshold, args.min_structures,
                                            subunits=subunits)
        except (ValueError, OSError) as e:
            log.warning("%s: could not predict state: %s", p["id"], e)
            row["status"] = f"failed: {e}"
            rows.append(row)
            continue
        all_results[p["id"]] = res
        means = res["per_subunits_mean_iptm"]
        n_scored = sum(1 for n in subunits if n in means)
        row.update({
            "predicted_subunits": res["predicted_state"] or 1,
            "predicted_state": res["predicted_state_bin"],
            "max_mean_iptm": f"{res['max_iptm']:.4f}",
            "top_multimer_subunits": res["top_multimer_subunits"] or "",
            "status": "ok" if n_scored == len(subunits) else
                      f"incomplete: {n_scored}/{len(subunits)} stoichiometries scored",
        })
        for n in subunits:
            row[f"mean_iptm_{n}mer"] = f"{means[n]:.4f}" if n in means else ""
        rows.append(row)
        log.info("%s: %s (max mean ipTM %.3f)", p["id"], res["predicted_state_bin"], res["max_iptm"])

    tsv = outdir / "phlegm_predictions.tsv"
    with open(tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", restval="")
        writer.writeheader()
        writer.writerows(rows)
    with open(outdir / "phlegm_predictions.pkl", "wb") as fh:
        pickle.dump(all_results, fh)

    if args.extended_metrics:
        mdir = outdir / "metrics"
        mdir.mkdir(exist_ok=True)
        for p in proteins:
            if p["id"] in all_results:
                run_cmd([sys.executable, SCRIPT_DIR / "extract_metrics.py",
                         "-i", outdir / "structures" / p["id"], "-p", p["id"],
                         "-o", mdir / f"{p['id']}_metrics"])  # script appends .pkl

    n_ok = sum(1 for r in rows if r.get("status") == "ok")
    log.info("Predictions written to %s (%d/%d protein(s) fully scored)", tsv, n_ok, len(rows))
    return n_ok == len(rows)


# --------------------------------------------------------------------------- #
# Commands
# --------------------------------------------------------------------------- #

def cmd_msa(args, outdir):
    proteins = load_fasta_proteins(args.input)
    save_registry(proteins, outdir)
    log.info("Loaded %d protein(s) from %s", len(proteins), args.input)
    step_msa(proteins, args, outdir)
    return 0


def cmd_fold(args, outdir):
    failed = step_fold(load_registry(outdir), args, outdir)
    return 1 if failed else 0


def cmd_predict(args, outdir):
    return 0 if step_predict(load_registry(outdir), args, outdir) else 2


def cmd_download_weights(args, outdir=None):
    wdir = Path(args.weights_dir)
    if weights_present(wdir):
        log.info("AlphaFold weights already present in %s", wdir)
        return 0
    download_weights(wdir)
    return 0


def cmd_run(args, outdir):
    proteins = load_fasta_proteins(args.input)
    save_registry(proteins, outdir)
    log.info("Loaded %d protein(s) from %s", len(proteins), args.input)
    step_msa(proteins, args, outdir)
    if args.dry_run:
        log.info("Dry run: fold commands depend on the MSAs and are not shown")
        return 0
    failed = step_fold(proteins, args, outdir)
    ok = step_predict(proteins, args, outdir)
    return 1 if failed else (0 if ok else 2)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def build_parser():
    fmt = argparse.ArgumentDefaultsHelpFormatter

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-o", "--outdir", required=True,
                        help="Output directory shared by all commands (re-use it to resume)")

    inp = argparse.ArgumentParser(add_help=False)
    inp.add_argument("-i", "--input", required=True, help="FASTA file with one or more protein sequences")

    dry = argparse.ArgumentParser(add_help=False)
    dry.add_argument("--dry-run", action="store_true", help="Print commands without running them")

    stoich = argparse.ArgumentParser(add_help=False)
    stoich.add_argument("--min-oligomer", type=int, default=2, help="Smallest number of subunits")
    stoich.add_argument("--max-oligomer", type=int, default=10, help="Largest number of subunits")

    msa = argparse.ArgumentParser(add_help=False)
    g = msa.add_argument_group("MSA options")
    g.add_argument("--db", help="Local ColabFold database directory for colabfold_search. "
                                "If omitted, the public ColabFold MMseqs2 server is used.")
    g.add_argument("--msa-dir", help="Directory of precomputed monomer MSAs named <id>.a3m")
    g.add_argument("--threads", type=int, default=os.cpu_count() or 4, help="Threads for colabfold_search")
    g.add_argument("--db-load-mode", type=int, choices=[0, 1, 2, 3],
                   help="MMseqs2 --db-load-mode for colabfold_search")

    fold = argparse.ArgumentParser(add_help=False)
    g = fold.add_argument_group("Fold options (defaults match the PHLEGM manuscript)")
    g.add_argument("--num-models", type=int, default=5, choices=range(1, 6), metavar="{1..5}",
                   help="AlphaFold-Multimer models per stoichiometry")
    g.add_argument("--num-recycle", type=int, default=3, help="Recycles per model")
    g.add_argument("--colabfold-args", default="",
                   help="Extra arguments for colabfold_batch, e.g. \"--amber --num-relax 1\"")
    g.add_argument("--keep-going", action="store_true",
                   help="Continue with other proteins if one fails")
    g.add_argument("--allow-cpu", action="store_true",
                   help="Run even if JAX cannot see a GPU (very slow)")

    weights = argparse.ArgumentParser(add_help=False)
    weights.add_argument("--weights-dir", default=str(default_weights_dir()),
                         help="AlphaFold weights directory (shared with ColabFold's cache; "
                              "can also be set with $PHLEGM_WEIGHTS_DIR)")

    pred = argparse.ArgumentParser(add_help=False)
    g = pred.add_argument_group("Predict options")
    g.add_argument("-t", "--threshold", type=float, default=0.65,
                   help="Mean ipTM threshold below which a protein is called a monomer")
    g.add_argument("-m", "--min-structures", type=int, default=2,
                   help="Minimum models with an ipTM score needed per stoichiometry")
    g.add_argument("--extended-metrics", action="store_true",
                   help="Also run extract_metrics.py (pDockQ2, interface contacts, PAE)")

    ap = argparse.ArgumentParser(prog="phlegm", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    sub = ap.add_subparsers(dest="command", required=True, metavar="COMMAND")

    sub.add_parser("msa", parents=[inp, common, msa, dry], formatter_class=fmt,
                   help="Build one MSA per protein (CPU)",
                   description="Build one MSA per protein, with a local ColabFold database (--db) "
                               "or the public ColabFold MMseqs2 server.")
    p = sub.add_parser("fold", parents=[common, stoich, fold, weights, dry], formatter_class=fmt,
                       help="Predict complexes for each stoichiometry (GPU)",
                       description="Reuse each protein's MSA to predict 2-mer..10-mer complexes "
                                   "with AlphaFold-Multimer. Needs only the AlphaFold weights.")
    p.add_argument("--shard", type=parse_shard, metavar="I/N",
                   help="Process only every N-th protein starting at I (for job arrays)")
    sub.add_parser("predict", parents=[common, stoich, pred], formatter_class=fmt,
                   help="Call oligomeric states from ipTM scores",
                   description="Call each protein's oligomeric state from the mean ipTM of each "
                               "stoichiometry and write phlegm_predictions.tsv.")
    sub.add_parser("run", parents=[inp, common, msa, stoich, fold, weights, pred, dry], formatter_class=fmt,
                   help="Run msa, fold and predict on one machine",
                   description="Run the whole pipeline (msa, fold, predict) on one machine.")
    sub.add_parser("download-weights", parents=[weights], formatter_class=fmt,
                   help="Download the AlphaFold-Multimer weights (~4 GB)",
                   description="Download the AlphaFold-Multimer v3 weights. 'phlegm fold' does this "
                               "automatically when needed, but on HPC run it once on a login node, "
                               "because GPU nodes often have no internet access.")
    return ap


def parse_args(argv=None):
    ap = build_parser()
    args = ap.parse_args(argv)
    if getattr(args, "db", None) and getattr(args, "msa_dir", None):
        ap.error("--db and --msa-dir are mutually exclusive")
    if getattr(args, "db", None) and not Path(args.db).is_dir():
        ap.error(f"--db directory not found: {args.db}")
    if hasattr(args, "min_oligomer") and not 2 <= args.min_oligomer <= args.max_oligomer:
        ap.error("require 2 <= --min-oligomer <= --max-oligomer")
    if hasattr(args, "num_models") and hasattr(args, "min_structures") \
            and args.min_structures > args.num_models:
        ap.error("--min-structures cannot exceed --num-models")
    return args


def setup_logging(outdir, args):
    name = f"phlegm_{args.command}"
    if getattr(args, "shard", None):
        name += f"_shard{args.shard[0]}of{args.shard[1]}"
    (outdir / "logs").mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%Y-%m-%d %H:%M:%S")
    log.setLevel(logging.INFO)
    for handler in (logging.StreamHandler(sys.stdout),
                    logging.FileHandler(outdir / "logs" / f"{name}.log")):
        handler.setFormatter(fmt)
        log.addHandler(handler)


def main(argv=None):
    args = parse_args(argv)
    if args.command == "download-weights":
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s",
                                               "%Y-%m-%d %H:%M:%S"))
        log.setLevel(logging.INFO)
        log.addHandler(handler)
        try:
            return cmd_download_weights(args)
        except (PhlegmError, OSError) as e:
            log.error("%s", e)
            return 1

    outdir = Path(args.outdir).resolve()
    if args.command in ("fold", "predict") and not outdir.is_dir():
        sys.exit(f"Error: output directory {outdir} not found. Run 'phlegm msa' first.")
    outdir.mkdir(parents=True, exist_ok=True)
    setup_logging(outdir, args)
    log.info("PHLEGM %s: %s", __version__, args.command)

    commands = {"msa": cmd_msa, "fold": cmd_fold, "predict": cmd_predict, "run": cmd_run}
    try:
        return commands[args.command](args, outdir)
    except (PhlegmError, OSError) as e:
        log.error("%s", e)
        return 1


if __name__ == "__main__":
    sys.exit(main())
