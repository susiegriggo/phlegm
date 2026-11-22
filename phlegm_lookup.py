"""
phlegm lookup: assign homo-oligomeric states to proteins from the precomputed
PHLEGM resource, using each protein's PHROG assignment (e.g. from pharokka or phold).

The resource (data/phlegm_phrog_predictions.csv.gz; Supplementary Table of the
PHLEGM manuscript) holds the PHLEGM prediction for the representative sequence of
22,554 PHROGs. A protein inherits the prediction of its PHROG's representative.
Proteins whose PHROG is not in the resource can be folded with 'phlegm run'.
"""

import csv
import gzip
import logging
import re
from pathlib import Path

log = logging.getLogger("phlegm")

DEFAULT_RESOURCE = Path(__file__).resolve().parent / "data" / "phlegm_phrog_predictions.csv.gz"
STRUCTURE_BASE = "https://projects.pawsey.org.au/phrogstructures-public"
MAX_MODELLED_LENGTH = 350   # PHROG representatives longer than this were not modelled

OUTPUT_COLUMNS = [
    "protein_id", "phrog", "status", "predicted_state", "predicted_subunits", "max_mean_iptm",
    "phrog_annotation", "phrog_category", "protein_length_aa", "representative_length_aa",
    "identity_to_representative", "coverage_of_protein", "structure_url", "note",
]


class PhlegmLookupError(Exception):
    pass


def parse_phrog(value):
    """'phrog_123', '123' or 'PHROG 123' -> 123; 'No_PHROG', '', 'NA', ... -> None."""
    if value is None:
        return None
    m = re.fullmatch(r"\s*(?:phrog[_ ]?)?0*(\d+)\s*", str(value), flags=re.IGNORECASE)
    return int(m.group(1)) if m else None


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="", encoding="utf-8-sig")
    return open(path, newline="", encoding="utf-8-sig")


def load_resource(path=DEFAULT_RESOURCE):
    """Read the per-PHROG PHLEGM table into {phrog_number: row}."""
    path = Path(path)
    if not path.is_file():
        raise PhlegmLookupError(f"PHLEGM resource table not found: {path}")
    with open_text(path) as fh:
        sample = fh.read(4096)
        fh.seek(0)
        dialect = "excel-tab" if sample.count("\t") > sample.count(",") else "excel"
        rows = list(csv.DictReader(fh, dialect=dialect))
    resource = {}
    for row in rows:
        key = parse_phrog((row.get("phrog") or row.get("﻿phrog") or ""))
        if key is not None:
            resource[key] = row
    if not resource:
        raise PhlegmLookupError(f"No PHROG rows could be read from {path}")
    return resource


def read_annotations(path):
    """
    Return [(protein_id, phrog_or_None)] from one of:
      - phold:    *_per_cds_predictions.tsv          (columns cds_id, phrog)
      - pharokka: *_cds_final_merged_output.tsv      (columns gene, phrog)
      - any TSV/CSV whose first column is the protein ID and second column the PHROG
        (a header line is optional)
    """
    with open_text(path) as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        rows = [r for r in csv.reader(fh, delimiter=delim) if r and not r[0].startswith("#")]
    if not rows:
        raise PhlegmLookupError(f"No rows in {path}")

    header = [h.strip().lstrip("﻿").lower() for h in rows[0]]
    for id_col in ("cds_id", "gene", "protein_id", "id"):
        if id_col in header and "phrog" in header:
            i, j = header.index(id_col), header.index("phrog")
            fmt = {"cds_id": "phold", "gene": "pharokka"}.get(id_col, "table")
            log.info("Reading %s annotations from %s", fmt, path)
            return [(r[i].strip(), parse_phrog(r[j])) for r in rows[1:] if len(r) > max(i, j)]

    if len(rows[0]) < 2:
        raise PhlegmLookupError(f"{path}: expected at least two columns (protein ID, PHROG)")
    # generic two-column table; skip a header if its PHROG cell is not a PHROG
    start = 1 if parse_phrog(rows[0][1]) is None and rows[0][1].strip().lower() in (
        "phrog", "phrog_id", "phrogs") else 0
    log.info("Reading two-column (protein ID, PHROG) table from %s", path)
    return [(r[0].strip(), parse_phrog(r[1])) for r in rows[start:] if len(r) >= 2]


def read_proteins(path):
    """FASTA -> {id: sequence}. IDs are the first word of the header; for headers like
    'contig:CDS_0001' (phold) the part after the last ':' is also accepted."""
    seqs, name = {}, None
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0] if len(line) > 1 else ""
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.replace("*", ""))
    out = {}
    for name, parts in seqs.items():
        seq = "".join(parts).upper()
        out[name] = seq
        if ":" in name:
            out.setdefault(name.rsplit(":", 1)[1], seq)
    return out


def make_aligner():
    """Local protein aligner (Smith-Waterman, BLOSUM62), or None without Biopython."""
    try:
        from Bio import Align
        from Bio.Align import substitution_matrices
    except ImportError:
        log.warning("Biopython not available: identity to the PHROG representative is not reported")
        return None
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


_NONSTANDARD = re.compile(r"[^ARNDCQEGHILKMFPSTWYVBZX]")


def identity_and_coverage(aligner, query, target):
    """Percent identity over aligned columns (ignoring masked 'X' residues), and the percent
    of the query covered by the local alignment."""
    query, target = _NONSTANDARD.sub("X", query), _NONSTANDARD.sub("X", target)
    aln = aligner.align(query, target)[0]
    blocks_q, blocks_t = aln.aligned
    if len(blocks_q) == 0:
        return 0.0, 0.0
    compared = identical = 0
    for (qs, qe), (ts, te) in zip(blocks_q, blocks_t):
        for a, b in zip(query[qs:qe], target[ts:te]):
            if a == "X" or b == "X":
                continue
            compared += 1
            identical += a == b
    covered = blocks_q[-1][1] - blocks_q[0][0]
    identity = 100.0 * identical / compared if compared else 0.0
    return identity, 100.0 * covered / len(query)


def structure_url(phrog, oligomer):
    """Stable URL of the PHLEGM resource structure (oligomer for multimer calls)."""
    group = f"phrog_{phrog // 100:03d}"
    kind = "oligomer_pdb" if oligomer else "pdb"
    return f"{STRUCTURE_BASE}/{kind}/{group}/phrog_{phrog:05d}.pdb.gz"


def run_lookup(annotations, outdir, proteins_fasta=None, resource_path=DEFAULT_RESOURCE):
    outdir = Path(outdir)
    resource = load_resource(resource_path)
    log.info("PHLEGM resource: %d PHROGs (%s)", len(resource), resource_path)
    entries = read_annotations(annotations)
    proteins = read_proteins(proteins_fasta) if proteins_fasta else {}
    aligner = make_aligner() if proteins else None
    if proteins:
        n_x = sum(s.count("X") for s in proteins.values())
        n_all = sum(len(s) for s in proteins.values()) or 1
        if n_x / n_all > 0.02:
            log.warning("%.0f%% of residues in %s are 'X' (masked); identity to the PHROG "
                        "representative is computed on unmasked residues only. Unmasked "
                        "sequences (e.g. pharokka's prodigal-gv.faa) give better estimates.",
                        100 * n_x / n_all, proteins_fasta)

    missing_seq = [pid for pid, _ in entries if proteins and pid not in proteins]
    if missing_seq:
        log.warning("%d protein(s) in %s have no sequence in %s (e.g. %s)",
                    len(missing_seq), annotations, proteins_fasta, missing_seq[0])

    results = []
    for pid, phrog in entries:
        seq = proteins.get(pid)
        row = {c: "" for c in OUTPUT_COLUMNS}
        row.update(protein_id=pid, phrog=f"phrog_{phrog}" if phrog is not None else "",
                   protein_length_aa=len(seq) if seq else "")
        hit = resource.get(phrog) if phrog is not None else None
        if phrog is None:
            row.update(status="no_phrog", note="no PHROG assigned; fold with 'phlegm run'")
        elif hit is None:
            row.update(status="not_in_resource")
            if seq and len(seq) > MAX_MODELLED_LENGTH:
                row["note"] = (f"PHROG not in resource; protein is >{MAX_MODELLED_LENGTH} aa "
                               "(beyond the modelled size range)")
            else:
                row["note"] = ("PHROG not in resource (representative >350 aa or low-confidence "
                               "monomer model); fold with 'phlegm run'")
        else:
            state = hit.get("phlegm_predicted_state", "")
            subunits = hit.get("phlegm_predicted_subunits", "")
            iptm = hit.get("phlegm_max_mean_iptm", "")
            row.update(
                status="predicted",
                predicted_state=state,
                predicted_subunits=subunits,
                max_mean_iptm=f"{float(iptm):.3f}" if iptm not in ("", None) else "",
                phrog_annotation=hit.get("functional_annotation", ""),
                phrog_category=hit.get("functional_category", ""),
                representative_length_aa=hit.get("sequence_length_aa", ""),
                structure_url=structure_url(phrog, oligomer=state not in ("", "Monomer")),
            )
            rep = (hit.get("representative_sequence") or "").upper()
            if aligner and seq and rep:
                ident, cov = identity_and_coverage(aligner, seq, rep)
                row.update(identity_to_representative=f"{ident:.1f}",
                           coverage_of_protein=f"{cov:.1f}")
                if ident < 30 or cov < 50:
                    row["note"] = ("low sequence similarity to the PHROG representative (PHROGs "
                                   "assigned by structure can be remote homologs); the prediction "
                                   "is less certain, consider 'phlegm run'")
        results.append(row)

    out_tsv = outdir / "phlegm_lookup.tsv"
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(results)

    to_fold = [r for r in results if r["status"] != "predicted"]
    if proteins:
        out_faa = outdir / "to_fold.faa"
        n_written = 0
        with open(out_faa, "w") as fh:
            for r in to_fold:
                seq = proteins.get(r["protein_id"])
                if seq:
                    fh.write(f">{r['protein_id']}\n{seq}\n")
                    n_written += 1
        n_long = sum(1 for r in to_fold if r["protein_length_aa"] and
                     int(r["protein_length_aa"]) > MAX_MODELLED_LENGTH)
        log.info("Wrote %d protein(s) without a lookup result to %s", n_written, out_faa)
        if n_long:
            log.info("  %d of these are >%d aa: AlphaFold-Multimer up to 10-mers needs a "
                     "large GPU for these", n_long, MAX_MODELLED_LENGTH)

    n = len(results)
    n_pred = sum(r["status"] == "predicted" for r in results)
    n_multi = sum(r["status"] == "predicted" and r["predicted_state"] not in ("", "Monomer")
                  for r in results)
    log.info("%d protein(s): %d predicted from the resource (%d multimer, %d monomer), "
             "%d not in resource, %d without a PHROG",
             n, n_pred, n_multi, n_pred - n_multi,
             sum(r["status"] == "not_in_resource" for r in results),
             sum(r["status"] == "no_phrog" for r in results))
    log.info("Results written to %s", out_tsv)
    if to_fold:
        hint = f"phlegm run -i {outdir / 'to_fold.faa'} -o <outdir>" if proteins else \
            "phlegm run -i <their sequences>.faa -o <outdir>"
        log.info("To predict the remaining %d protein(s): %s", len(to_fold), hint)
    return results
