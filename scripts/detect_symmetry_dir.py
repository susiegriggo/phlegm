#!/usr/bin/env python3
"""
detect_symmetry_dir.py
----------------------
Scan a directory of PDB/mmCIF files, estimate point-group symmetry (Cn or Dn)
for homooligomeric assemblies, and write a TSV summary.

Usage:
  python detect_symmetry_dir.py /path/to/dir --glob "*.pdb" --out symmetries.tsv
  python detect_symmetry_dir.py /path/to/dir --glob "*.cif" --recursive --out symmetries.tsv

Notes:
- Output is a TSV with one row per file.
- Focuses on Cn/Dn symmetries, which cover most homomeric enzymes.
- Requires: numpy, biopython
"""

import os
import sys
import csv
import argparse
import math
from collections import defaultdict
import numpy as np

from Bio.PDB import PDBParser, MMCIFParser, is_aa

# ------------------------- math utils -------------------------

def kabsch(P, Q):
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)
    P0 = P - Pc
    Q0 = Q - Qc
    H = P0.T @ Q0
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = Qc - R @ Pc
    P_aligned = (R @ P.T).T + t
    rmsd = np.sqrt(np.mean(np.sum((P_aligned - Q) ** 2, axis=1)))
    return R, t, rmsd

def unit(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    return v / n if n > 0 else v

def project_to_plane(v, n):
    n = unit(n)
    return v - np.dot(v, n) * n

def circular_std(angles):
    angles = np.asarray(angles, dtype=float)
    C = np.mean(np.cos(angles))
    S = np.mean(np.sin(angles))
    R = math.hypot(C, S)
    if R <= 1e-12:
        return math.pi / np.sqrt(3)
    return math.sqrt(-2 * math.log(R))

# ------------------------- structure utils -------------------------

def load_structure(path):
    fname = os.path.basename(path).lower()
    if fname.endswith(".cif") or fname.endswith(".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure(id="X", file=path)
    return structure

def get_chain_ca_coords(chain):
    coords = []
    for res in chain.get_residues():
        if not is_aa(res, standard=True):
            continue
        if "CA" in res:
            coords.append(res["CA"].get_coord())
    return np.array(coords, dtype=float)

def collect_chains(structure):
    model = next(structure.get_models())
    chains = {}
    for chain in model.get_chains():
        coords = get_chain_ca_coords(chain)
        if len(coords) < 20:
            continue
        com = coords.mean(axis=0)
        chains[chain.id] = {"coords": coords, "com": com}
    return chains

def cluster_by_length(chains, len_tol=0.05):
    items = [(cid, len(info["coords"])) for cid, info in chains.items()]
    if not items:
        return []
    items.sort(key=lambda x: x[1])
    clusters = []
    for cid, L in items:
        placed = False
        for cl in clusters:
            Lref = cl[0][1]
            if abs(L - Lref) / max(Lref, 1) <= len_tol:
                cl.append((cid, L))
                placed = True
                break
        if not placed:
            clusters.append([(cid, L)])
    clusters = [[cid for cid, _ in cl] for cl in clusters]
    clusters.sort(key=len, reverse=True)
    return clusters

# ------------------------- symmetry detection -------------------------

def estimate_cn_from_com(coms, tol_angle_deg=8.0):
    coms = np.asarray(coms, dtype=float)
    K = len(coms)
    if K < 2:
        return False, 1, np.array([0,0,1.0]), np.zeros(1), 0.0, 0.0, 0.0

    G = coms.mean(axis=0)
    X = coms - G
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    axis = unit(Vt[0])

    # Plane basis
    tmp = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = unit(np.cross(axis, tmp))
    e2 = unit(np.cross(axis, e1))

    XY = np.array([project_to_plane(v, axis) for v in X])
    xs = XY @ e1
    ys = XY @ e2
    radii = np.hypot(xs, ys)
    theta = np.arctan2(ys, xs)

    order = np.argsort(theta)
    theta = theta[order]
    radii = radii[order]
    zvals = (X @ axis)[order]

    target = 2 * math.pi / K
    # circular dispersion around uniform spacing
    diffs = (theta - (theta[0] + np.arange(K) * target))
    # wrap to (-pi, pi]
    diffs = (diffs + math.pi) % (2 * math.pi) - math.pi
    angle_std_deg = float(np.degrees(circular_std(diffs)))
    radius_cv = float(np.std(radii) / (np.mean(radii) + 1e-8))
    z_span = float(np.max(zvals) - np.min(zvals))

    is_uniform = angle_std_deg <= tol_angle_deg and radius_cv < 0.2
    return is_uniform, K, axis, theta, angle_std_deg, radius_cv, z_span

def test_dn_symmetry(coms, axis, theta, tol_match_deg=12.0, tol_z=3.0, tol_r_frac=0.25):
    coms = np.asarray(coms, dtype=float)
    K = len(coms)
    G = coms.mean(axis=0)
    X = coms - G

    tmp = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = unit(np.cross(axis, tmp))
    e2 = unit(np.cross(axis, e1))

    XY = np.array([project_to_plane(v, axis) for v in X])
    xs = XY @ e1
    ys = XY @ e2
    radii = np.hypot(xs, ys)
    th = np.arctan2(ys, xs)
    z = X @ axis

    used = np.zeros(K, dtype=bool)
    tol_angle = math.radians(tol_match_deg)

    for i in range(K):
        target_theta = (-th[i]) % (2*math.pi)
        r_i, z_i = radii[i], z[i]
        best_j = -1
        for j in range(K):
            if used[j] or j == i:
                continue
            # angular distance modulo 2pi
            d = abs(((th[j] - target_theta + math.pi) % (2*math.pi)) - math.pi)
            r_ok = abs(radii[j] - r_i) <= tol_r_frac * max(1.0, r_i)
            z_ok = abs(z[j] + z_i) <= tol_z
            if d <= tol_angle and r_ok and z_ok:
                best_j = j
                break
        if best_j < 0:
            return False
        used[i] = True
        used[best_j] = True
    return True

def superposition_rmsds(chains, ref_id):
    ref = chains[ref_id]["coords"]
    rmsds = {}
    for cid, info in chains.items():
        P = info["coords"]
        L = min(len(P), len(ref))
        _, _, rmsd = kabsch(P[:L], ref[:L])
        rmsds[cid] = rmsd
    arr = np.array(list(rmsds.values()))
    return rmsds, float(np.mean(arr)), float(np.max(arr))

def analyze_pdb(path):
    res = {
        "file": path,
        "basename": os.path.basename(path),
        "n_chains_total": 0,
        "cluster_size": 0,
        "symmetry": "None",
        "n": 1,
        "angle_std_deg": np.nan,
        "radius_cv": np.nan,
        "z_span": np.nan,
        "rmsd_mean": np.nan,
        "rmsd_max": np.nan,
        "status": "OK"
    }
    try:
        structure = load_structure(path)
        chains = collect_chains(structure)
        res["n_chains_total"] = len(chains)
        if not chains:
            res["status"] = "NO_VALID_CHAINS"
            return res

        clusters = cluster_by_length(chains)
        if not clusters:
            res["status"] = "NO_CLUSTERS"
            return res

        main = clusters[0]
        res["cluster_size"] = len(main)

        if len(main) == 1:
            res["symmetry"] = "C1"
            res["n"] = 1
            rmsds, mean_rmsd, max_rmsd = superposition_rmsds({main[0]: chains[main[0]]}, main[0])
            res["rmsd_mean"] = mean_rmsd
            res["rmsd_max"] = max_rmsd
            return res

        coms = np.array([chains[cid]["com"] for cid in main])
        is_cn, n, axis, theta, angle_std_deg, radius_cv, z_span = estimate_cn_from_com(coms)
        res["angle_std_deg"] = float(angle_std_deg)
        res["radius_cv"] = float(radius_cv)
        res["z_span"] = float(z_span)

        rmsds, mean_rmsd, max_rmsd = superposition_rmsds({cid: chains[cid] for cid in main}, main[0])
        res["rmsd_mean"] = mean_rmsd
        res["rmsd_max"] = max_rmsd

        if not is_cn:
            res["symmetry"] = "None"
            res["n"] = len(main)
            return res

        is_dn = test_dn_symmetry(coms, axis, theta)
        if is_dn and len(main) % 2 == 0:
            res["symmetry"] = f"D{len(main)//2}"
            res["n"] = len(main) // 2
        else:
            res["symmetry"] = f"C{len(main)}"
            res["n"] = len(main)

        return res

    except Exception as e:
        res["status"] = f"ERROR: {e}"
        return res

def collect_paths(root, pattern, recursive=False):
    import glob
    if recursive:
        paths = glob.glob(os.path.join(root, "**", pattern), recursive=True)
    else:
        paths = glob.glob(os.path.join(root, pattern))
    # Filter typical structure extensions if pattern is generic
    return sorted(paths)

def main():
    ap = argparse.ArgumentParser(description="Batch-detect Cn/Dn symmetries from PDB/mmCIF files and write a TSV.")
    ap.add_argument("directory", help="Directory containing PDB/mmCIF files")
    ap.add_argument("--glob", default="*.pdb", help="Glob pattern (default: *.pdb). Try '*.cif' for mmCIF.")
    ap.add_argument("--recursive", action="store_true", help="Recurse into subdirectories")
    ap.add_argument("--out", default="symmetries.tsv", help="Output TSV path (default: symmetries.tsv)")
    args = ap.parse_args()

    if not os.path.isdir(args.directory):
        print("Error: input must be a directory", file=sys.stderr)
        sys.exit(2)

    paths = collect_paths(args.directory, args.glob, recursive=args.recursive)
    if not paths:
        print("No files matched the pattern.", file=sys.stderr)
        sys.exit(1)

    fields = [
        "file","basename","symmetry","n","cluster_size","n_chains_total",
        "rmsd_mean","rmsd_max","angle_std_deg","radius_cv","z_span","status"
    ]

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for p in paths:
            r = analyze_pdb(p)
            w.writerow({k: r.get(k, "") for k in fields})

    print(f"Wrote {len(paths)} rows to {args.out}")

if __name__ == "__main__":
    main()
