#!/usr/bin/env python
"""
PHLEGM: run the full homo-oligomeric state prediction pipeline from a FASTA file.

Steps
  1. Build one MSA per protein (ColabFold MMseqs2 server, or a local ColabFold
     database with colabfold_search).
  2. Reuse each MSA to build inputs for every stoichiometry (2-mer..10-mer).
  3. Predict complexes with colabfold_batch (AlphaFold-Multimer v3).
  4. Predict the oligomeric state from the mean ipTM of each stoichiometry.

Every step skips work that is already finished, so a stopped run can be resumed
by re-running the same command.

Example
  python phlegm.py -i proteins.fasta -o phlegm_out
  python phlegm.py -i proteins.fasta -o phlegm_out --db /path/to/colabfold_db
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

__version__ = "1.0.0"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYXBZUO")

log = logging.getLogger("phlegm")


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
                    raise ValueError(f"{path}: sequence data found before the first '>' header")
                seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    if not records:
        raise ValueError(f"{path}: no sequences found")
    return records


def safe_id(header):
    """First word of the header, restricted to characters safe for file names.

    '.' is replaced as well because PHLEGM uses '<id>.<subunits>' to name jobs.
    """
    first = header.split()[0] if header.split() else "protein"
    return re.sub(r"[^A-Za-z0-9_-]", "_", first)


def load_proteins(path):
    proteins, seen = [], {}
    for header, seq in read_fasta(path):
        pid = safe_id(header)
        seq = seq.upper().replace(" ", "").rstrip("*")
        bad = set(seq) - VALID_AA
        if not seq:
            raise ValueError(f"'{header}' has an empty sequence")
        if bad:
            raise ValueError(f"'{header}' contains invalid characters: {''.join(sorted(bad))}")
        if pid in seen:
            raise ValueError(
                f"Headers '{seen[pid]}' and '{header}' both map to the ID '{pid}'. "
                "Please make the first word of each FASTA header unique.")
        seen[pid] = header
        proteins.append({"id": pid, "header": header, "sequence": seq})
    return proteins


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def run(cmd, dry_run=False):
    log.info("Running: %s", " ".join(shlex.quote(str(c)) for c in cmd))
    if dry_run:
        return
    result = subprocess.run([str(c) for c in cmd])
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {cmd[0]}")


def require(executable):
    if shutil.which(executable) is None:
        sys.exit(f"Error: '{executable}' was not found on PATH. Did you run 'bash install.sh' "
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


# --------------------------------------------------------------------------- #
# Step 1: monomer MSAs
# --------------------------------------------------------------------------- #

def generate_msas(proteins, args, outdir):
    msa_dir = outdir / "msas"
    msa_dir.mkdir(parents=True, exist_ok=True)

    todo = [p for p in proteins if not (msa_dir / f"{p['id']}.a3m").is_file()]
    if not todo:
        log.info("Step 1/4: all MSAs already present, skipping")
        return
    log.info("Step 1/4: generating MSAs for %d protein(s)", len(todo))

    if args.msa_dir:
        for p in todo:
            src = Path(args.msa_dir) / f"{p['id']}.a3m"
            if not src.is_file():
                raise FileNotFoundError(f"--msa-dir: expected {src}")
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
        run(cmd, args.dry_run)
    else:
        require("colabfold_batch")
        log.info("Using the public ColabFold MMseqs2 server for MSAs "
                 "(use --db for a local database on large datasets)")
        run(["colabfold_batch", query, work, "--msa-only"], args.dry_run)

    if args.dry_run:
        return

    # colabfold_search names monomer MSAs by query index (0.a3m, 1.a3m, ...);
    # colabfold_batch --msa-only names them by query ID. Accept either.
    for i, p in enumerate(todo):
        candidates = [work / f"{p['id']}.a3m", work / f"{i}.a3m"]
        found = next((c for c in candidates if c.is_file()), None)
        if found is None:
            raise FileNotFoundError(f"No MSA was produced for {p['id']} in {work}")
        entries = read_a3m(found)
        if not entries or entries[0][1].replace("-", "").upper() != p["sequence"]:
            raise ValueError(f"First sequence of {found} does not match the query {p['id']}")
        shutil.move(str(found), msa_dir / f"{p['id']}.a3m")
    shutil.rmtree(work)


# --------------------------------------------------------------------------- #
# Step 2: stoichiometry-specific MSAs
# --------------------------------------------------------------------------- #

def write_oligomer_msa(protein, monomer_a3m, n_subunits, out_path):
    """Reuse a monomer MSA as a homo-oligomer input for colabfold_batch.

    The '#<length>\\t<n>' header tells ColabFold to model n copies of the query.
    """
    seq = protein["sequence"]
    with open(out_path, "w") as fh:
        fh.write(f"#{len(seq)}\t{n_subunits}\n")
        fh.write(f">{protein['id']}\n{seq}\n")
        for header, s in read_a3m(monomer_a3m):
            fh.write(f"{header}\n{s}\n")


def generate_oligomer_msas(proteins, args, outdir):
    log.info("Step 2/4: building %d-mer to %d-mer inputs", args.min_oligomer, args.max_oligomer)
    for p in proteins:
        monomer = outdir / "msas" / f"{p['id']}.a3m"
        pdir = outdir / "oligomer_msas" / p["id"]
        pdir.mkdir(parents=True, exist_ok=True)
        for n in range(args.min_oligomer, args.max_oligomer + 1):
            out = pdir / f"{p['id']}.{n}.a3m"
            if not out.is_file() and not args.dry_run:
                write_oligomer_msa(p, monomer, n, out)


# --------------------------------------------------------------------------- #
# Step 3: structure prediction
# --------------------------------------------------------------------------- #

def predict_structures(proteins, args, outdir):
    require("colabfold_batch")
    log.info("Step 3/4: predicting structures with colabfold_batch")
    failed = []
    for i, p in enumerate(proteins, 1):
        log.info("[%d/%d] %s (%d aa)", i, len(proteins), p["id"], len(p["sequence"]))
        cmd = ["colabfold_batch",
               outdir / "oligomer_msas" / p["id"],
               outdir / "structures" / p["id"],
               "--model-type", "alphafold2_multimer_v3",
               "--num-models", args.num_models,
               "--num-recycle", args.num_recycle]
        cmd += shlex.split(args.colabfold_args or "")
        try:
            run(cmd, args.dry_run)
        except RuntimeError as e:
            if not args.keep_going:
                raise
            log.error("%s: %s (continuing because of --keep-going)", p["id"], e)
            failed.append(p["id"])
    return failed


# --------------------------------------------------------------------------- #
# Step 4: oligomeric state prediction
# --------------------------------------------------------------------------- #

def predict_states(proteins, args, outdir):
    log.info("Step 4/4: predicting oligomeric states (mean ipTM threshold %.2f)", args.threshold)
    subunits = list(range(max(2, args.min_oligomer), args.max_oligomer + 1))
    fields = (["protein_id", "header", "length", "predicted_subunits", "predicted_state",
               "max_mean_iptm", "top_multimer_subunits"]
              + [f"mean_iptm_{n}mer" for n in subunits] + ["status"])
    all_results, rows = {}, []

    for p in proteins:
        row = {"protein_id": p["id"], "header": p["header"], "length": len(p["sequence"])}
        try:
            _, res = predict_from_directory(str(outdir / "structures" / p["id"]), p["id"],
                                            args.threshold, args.min_structures)
        except (ValueError, FileNotFoundError) as e:
            log.warning("%s: could not predict state: %s", p["id"], e)
            row["status"] = f"failed: {e}"
            rows.append(row)
            continue
        all_results[p["id"]] = res
        means = res["per_subunits_mean_iptm"]
        row.update({
            "predicted_subunits": res["predicted_state"] or 1,
            "predicted_state": res["predicted_state_bin"],
            "max_mean_iptm": f"{res['max_iptm']:.4f}",
            "top_multimer_subunits": res["top_multimer_subunits"] or "",
            "status": "ok" if len(means) == len(subunits) else
                      f"incomplete: {len(means)}/{len(subunits)} stoichiometries scored",
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
    log.info("Predictions written to %s", tsv)
    return rows


def extended_metrics(proteins, args, outdir):
    log.info("Computing extended metrics (pDockQ2, interface contacts, PAE)")
    mdir = outdir / "metrics"
    mdir.mkdir(exist_ok=True)
    for p in proteins:
        run([sys.executable, SCRIPT_DIR / "extract_metrics.py",
             "-i", outdir / "structures" / p["id"],
             "-p", p["id"],
             "-o", mdir / f"{p['id']}_metrics"], args.dry_run)  # script appends .pkl


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        prog="phlegm",
        description="Predict homo-oligomeric states from protein sequences (FASTA -> predictions).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="FASTA file with one or more protein sequences")
    ap.add_argument("-o", "--outdir", required=True, help="Output directory (re-use it to resume a run)")
    ap.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")

    msa = ap.add_argument_group("MSA generation")
    msa.add_argument("--db", help="Local ColabFold database directory for colabfold_search. "
                                  "If omitted, the public ColabFold MMseqs2 server is used.")
    msa.add_argument("--msa-dir", help="Directory of precomputed monomer MSAs named <id>.a3m "
                                       "(skips MSA generation)")
    msa.add_argument("--threads", type=int, default=os.cpu_count() or 4,
                     help="Threads for colabfold_search")
    msa.add_argument("--db-load-mode", type=int, choices=[0, 1, 2, 3],
                     help="MMseqs2 --db-load-mode for colabfold_search")

    pred = ap.add_argument_group("Structure prediction (defaults match the PHLEGM manuscript)")
    pred.add_argument("--min-oligomer", type=int, default=2, help="Smallest number of subunits to model")
    pred.add_argument("--max-oligomer", type=int, default=10, help="Largest number of subunits to model")
    pred.add_argument("--num-models", type=int, default=5, choices=range(1, 6), metavar="{1..5}",
                      help="AlphaFold-Multimer models per stoichiometry")
    pred.add_argument("--num-recycle", type=int, default=3, help="Recycles per model")
    pred.add_argument("--colabfold-args", default="",
                      help="Extra arguments passed to colabfold_batch, e.g. \"--amber --num-relax 1\"")

    score = ap.add_argument_group("Oligomeric state prediction")
    score.add_argument("-t", "--threshold", type=float, default=0.65,
                       help="Mean ipTM threshold below which a protein is called a monomer")
    score.add_argument("-m", "--min-structures", type=int, default=2,
                       help="Minimum models with an ipTM score needed per stoichiometry")
    score.add_argument("--extended-metrics", action="store_true",
                       help="Also run extract_metrics.py (pDockQ2, interface contacts, PAE)")

    run_g = ap.add_argument_group("Run control")
    run_g.add_argument("--stop-after", choices=["msa", "alignments", "structures"],
                       help="Stop after this step (e.g. build MSAs on a CPU node, then "
                            "re-run on a GPU node to continue)")
    run_g.add_argument("--keep-going", action="store_true",
                       help="Continue with other proteins if structure prediction fails for one")
    run_g.add_argument("--dry-run", action="store_true", help="Print commands without running them")

    args = ap.parse_args(argv)
    if args.db and args.msa_dir:
        ap.error("--db and --msa-dir are mutually exclusive")
    if not 2 <= args.min_oligomer <= args.max_oligomer:
        ap.error("require 2 <= --min-oligomer <= --max-oligomer")
    if args.min_structures > args.num_models:
        ap.error("--min-structures cannot exceed --num-models")
    if args.db and not Path(args.db).is_dir():
        ap.error(f"--db directory not found: {args.db}")
    return args


def setup_logging(outdir):
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%Y-%m-%d %H:%M:%S")
    log.setLevel(logging.INFO)
    for handler in (logging.StreamHandler(sys.stdout), logging.FileHandler(outdir / "phlegm.log")):
        handler.setFormatter(fmt)
        log.addHandler(handler)


def main(argv=None):
    args = parse_args(argv)
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    setup_logging(outdir)
    log.info("PHLEGM %s", __version__)

    try:
        proteins = load_proteins(args.input)
    except (OSError, ValueError) as e:
        log.error("Could not read %s: %s", args.input, e)
        return 1
    log.info("Loaded %d protein(s) from %s", len(proteins), args.input)
    (outdir / "input").mkdir(exist_ok=True)
    write_fasta(proteins, outdir / "input" / "queries.fasta")

    generate_msas(proteins, args, outdir)
    if args.stop_after == "msa":
        return 0
    if args.dry_run:
        log.info("Dry run: remaining steps depend on MSAs; showing structure commands only")
    generate_oligomer_msas(proteins, args, outdir)
    if args.stop_after == "alignments":
        return 0
    failed = predict_structures(proteins, args, outdir)
    if args.stop_after == "structures" or args.dry_run:
        return 0

    rows = predict_states(proteins, args, outdir)
    if args.extended_metrics:
        extended_metrics([p for p in proteins if p["id"] not in failed], args, outdir)

    n_ok = sum(1 for r in rows if r.get("status") == "ok")
    log.info("Done: %d/%d protein(s) fully scored", n_ok, len(rows))
    return 0 if n_ok == len(rows) else 2


if __name__ == "__main__":
    sys.exit(main())
