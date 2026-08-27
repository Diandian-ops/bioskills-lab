"""
Faithful reproduction of bioSkills `alignment/structural-alignment` SKILL.md core patterns.

Covers:
  P1. Bio.PDB.Superimposer  — known-correspondence superposition + RMSD
  P2. TMalign (compiled binary) — pairwise TM-score / RMSD on real PDBs
  P3. Foldseek easy-search / easy-cluster — structural homolog search & clustering

All inputs are real downloaded PDBs from RCSB:
  1ubq.pdb  (ubiquitin, chain A)  — reference
  1ubi.pdb  (ubiquitin, chain A)  — close homolog of 1ubq
  1fmb.pdb  (different fold)      — negative / cross-fold control
"""
import json
import re
import subprocess
from pathlib import Path

from Bio.PDB import PDBParser, Superimposer

BASE = Path("/Users/zhangdiandian/RedBook/content/素材/010-structural-alignment")
PD = BASE / "pdbs"
TMALIGN = BASE / "tools" / "TMalign"
FOLDSEEK = "/Applications/anaconda3/envs/foldseek/bin/foldseek"

PDBS = {"1ubq": "1ubq.pdb", "1ubi": "1ubi.pdb", "1fmb": "1fmb.pdb"}


def section(title):
    print("\n" + "=" * 70 + f"\n{title}\n" + "=" * 70)


# ---------------------------------------------------------------------------
# P1. Bio.PDB.Superimposer — known correspondence (same-family ubiquitin)
# ---------------------------------------------------------------------------
def superimposer(ref_path, mob_path):
    p = PDBParser(QUIET=True)
    ref = p.get_structure("ref", str(ref_path))
    mob = p.get_structure("mob", str(mob_path))
    ca_ref = [a for a in ref.get_atoms() if a.get_id() == "CA"]
    ca_mob = [a for a in mob.get_atoms() if a.get_id() == "CA"]
    n = min(len(ca_ref), len(ca_mob))
    sup = Superimposer()
    sup.set_atoms(ca_ref[:n], ca_mob[:n])
    sup.apply(mob.get_atoms())
    return sup.rms, n, len(ca_ref), len(ca_mob)


# ---------------------------------------------------------------------------
# P2. TMalign binary parsing
# ---------------------------------------------------------------------------
def run_tmalign(a, b):
    out = subprocess.run(
        [str(TMALIGN), str(a), str(b), "-outfmt", "2"],
        capture_output=True, text=True,
    )
    text = out.stdout + out.stderr
    # -outfmt 2 -> tabular: PDBchain1 PDBchain2 TM1 TM2 RMSD ID1 ID2 IDali L1 L2 Lali
    tm = None
    rmsd = None
    aln_len = None
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) >= 5 and "pdb" in cols[0].lower():
            try:
                tm1, tm2 = float(cols[2]), float(cols[3])
                rmsd = float(cols[4])
                aln_len = int(cols[10]) if len(cols) > 10 else None
                tm = max(tm1, tm2)  # larger = normalisation by shorter chain
            except (ValueError, IndexError):
                pass
    return {
        "tm_score": tm,
        "tm1": None, "tm2": None,
        "rmsd": rmsd,
        "aligned_len": aln_len,
        "rc": out.returncode,
        "stderr_first": (out.stderr.strip().splitlines()[:1] or [""])[0],
    }


# ---------------------------------------------------------------------------
# P3. Foldseek
# ---------------------------------------------------------------------------
def run_foldseek_search(query, target_dir, out_prefix):
    tmp = str(BASE / "foldseek_tmp")
    out = str(BASE / f"{out_prefix}.m8")
    r = subprocess.run(
        [FOLDSEEK, "easy-search", str(query), str(target_dir), out, tmp,
         "--format-output", "query,target,alntmscore,qtmscore,ttmscore,lddt,bits,rmsd"],
        capture_output=True, text=True,
    )
    rows = []
    if Path(out).exists():
        for line in Path(out).read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) >= 7:
                rows.append({
                    "query": cols[0], "target": cols[1],
                    "alntmscore": cols[2], "qtmscore": cols[3],
                    "ttmscore": cols[4], "lddt": cols[5],
                    "bits": cols[6], "rmsd": cols[7] if len(cols) > 7 else None,
                })
    return {"rc": r.returncode, "rows": rows,
            "stderr_first": (r.stderr.strip().splitlines()[:1] or [""])[0]}


def run_foldseek_cluster(pdb_dir, out_prefix):
    tmp = str(BASE / "foldseek_tmp_cl")
    out = str(BASE / f"{out_prefix}_cluster")
    r = subprocess.run(
        [FOLDSEEK, "easy-cluster", str(pdb_dir) + "/", out, tmp, "--tmscore-threshold", "0.5"],
        capture_output=True, text=True,
    )
    # clustering writes <out>_cluster_cluster.tsv (member<TAB>representative)
    clustered = str(BASE / f"{out_prefix}_cluster_cluster.tsv")
    rows = []
    if Path(clustered).exists():
        for line in Path(clustered).read_text().splitlines():
            if not line.strip():
                continue
            cols = line.split("\t")
            member = cols[0]
            rep = cols[1] if len(cols) > 1 else cols[0]
            rows.append({"member": member, "rep": rep})
    return {"rc": r.returncode, "rows": rows,
            "stderr_first": (r.stderr.strip().splitlines()[:1] or [""])[0]}


def main():
    result = {}

    # P1
    section("P1 — Bio.PDB.Superimposer (1ubq -> 1ubi, known correspondence)")
    rms, n, nref, nmob = superimposer(str(PD / PDBS["1ubq"]), str(PD / PDBS["1ubi"]))
    print(f"Reference 1ubq CA atoms = {nref}")
    print(f"Mobile   1ubi CA atoms = {nmob}")
    print(f"Aligned (min) = {n}")
    print(f"RMSD = {rms:.4f} A over {n} CA atoms")
    result["superimposer"] = {"rmsd": round(rms, 4), "n_aligned": n,
                              "n_ref": nref, "n_mob": nmob}

    # P2
    section("P2 — TMalign pairwise (real binary)")
    pairs = [("1ubq", "1ubi"), ("1ubq", "1fmb"), ("1ubi", "1fmb")]
    tm_matrix = {}
    for a, b in pairs:
        r = run_tmalign(str(PD / PDBS[a]), str(PD / PDBS[b]))
        tm_matrix[f"{a}|{b}"] = r
        print(f"{a} vs {b}: TM-score={r['tm_score']}  RMSD={r['rmsd']}  "
              f"aligned_len={r['aligned_len']}  rc={r['rc']}")
        if r["stderr_first"]:
            print(f"   [note] {r['stderr_first']}")
    result["tmalign"] = tm_matrix

    # P3
    section("P3 — Foldseek easy-search (1ubq query vs pdbs folder)")
    fs = run_foldseek_search(str(PD / PDBS["1ubq"]), str(PD), "fs_search")
    print(f"easy-search rc={fs['rc']}")
    for row in fs["rows"]:
        print(f"  {row['query']} -> {row['target']}  alnTM={row['alntmscore']}  "
              f"qTM={row['qtmscore']}  tTM={row['ttmscore']}  LDDT={row['lddt']}")
    if fs["stderr_first"]:
        print(f"   [note] {fs['stderr_first']}")
    result["foldseek_search"] = fs

    section("P3 — Foldseek easy-cluster (all-vs-all, TM>0.5)")
    cl = run_foldseek_cluster(str(PD), "fs")
    print(f"easy-cluster rc={cl['rc']}")
    for row in cl["rows"]:
        print(f"  rep={row['rep']}  member={row['member']}")
    if cl["stderr_first"]:
        print(f"   [note] {cl['stderr_first']}")
    result["foldseek_cluster"] = cl

    Path(BASE / "structural_results.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))
    print("\nWrote structural_results.json")
    return result


if __name__ == "__main__":
    main()
