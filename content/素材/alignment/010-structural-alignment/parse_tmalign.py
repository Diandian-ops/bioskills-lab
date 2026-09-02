import re, json, glob, os
out = {"tool": "TM-align 20240303", "reference_metric": "max of the two length-normalised TM-scores (normalised by shorter chain)"}
pairs = []
for f in sorted(glob.glob("tmalign_*.txt")):
    if f == "tmalign__.txt":
        continue
    m = re.search(r"tmalign_(.+?)_(.+?)\.txt", f)
    if not m:
        continue
    a, b = m.group(1), m.group(2)
    txt = open(f).read()
    al = re.search(r"Aligned length=\s*(\d+)", txt)
    rmsd = re.search(r"RMSD=\s*([\d.]+)", txt)
    tms = re.findall(r"TM-score=\s*([\d.]+)", txt)
    l1 = re.search(r"Length of Structure_1:\s*(\d+)", txt)
    l2 = re.search(r"Length of Structure_2:\s*(\d+)", txt)
    rec = {
        "pair": f"{a} vs {b}",
        "query": a, "target": b,
        "len_query": int(l1.group(1)) if l1 else None,
        "len_target": int(l2.group(1)) if l2 else None,
        "aligned_length": int(al.group(1)) if al else None,
        "rmsd": float(rmsd.group(1)) if rmsd else None,
        "tm_score_norm_query": float(tms[0]) if len(tms) > 0 else None,
        "tm_score_norm_target": float(tms[1]) if len(tms) > 1 else None,
        "tm_score_fold_sim": max(float(tms[0]), float(tms[1])) if len(tms) >= 2 else None,
    }
    rec["same_fold_TM_gt_0.5"] = bool(rec["tm_score_fold_sim"] is not None and rec["tm_score_fold_sim"] > 0.5)
    pairs.append(rec)
out["pairs"] = pairs
json.dump(out, open("structural_results.json", "w"), indent=2)
print("WROTE structural_results.json with", len(pairs), "pairs")
for p in pairs:
    print(p["pair"], "TM=%.4f" % p["tm_score_fold_sim"], "RMSD=%.2f" % p["rmsd"], "aligned=%d" % p["aligned_length"], "same_fold=%s" % p["same_fold_TM_gt_0.5"])
