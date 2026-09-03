#!/usr/bin/env python3
# audit.py -- parse all produced files (skill口径: 9-column walk, attribute split),
# check conversions / phase / namespaces / FASTA consistency, write
# audit_report.txt + audit.json (figure data).
import json
import os
import re
from statistics import median

OUT = os.path.dirname(os.path.abspath(__file__))
truth = json.load(open(os.path.join(OUT, "truth.json")))
VER = re.compile(r"\.\d+$")


def attrs_gtf(s):
    d = {}
    for part in s.strip().split(";"):
        part = part.strip()
        if part:
            k, _, v = part.partition(" ")
            d[k] = v.strip().strip('"')
    return d


def attrs_gff3(s):
    d = {}
    for part in s.strip().split(";"):
        if part:
            k, _, v = part.partition("=")
            d[k] = v
    return d


def read_gxf(path):
    """returns rows plus per-transcript exon/cds coordinate lists"""
    feats, seqids, keys, idver = {}, set(), set(), 0
    tx_ex, tx_cds, phase_bad = {}, {}, 0
    col_ok = True
    gff3 = path.endswith(".gff3")
    for ln in open(path):
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.rstrip("\n").split("\t")
        if len(p) != 9:
            col_ok = False
            continue
        seqid, feat, s, e, strand, phase, attr = p[0], p[2], int(p[3]), int(p[4]), p[6], p[7], p[8]
        feats[feat] = feats.get(feat, 0) + 1
        seqids.add(seqid)
        a = attrs_gff3(attr) if gff3 else attrs_gtf(attr)
        keys.update(a.keys())
        tid = a.get("transcript_id") or a.get("Parent", "").replace("transcript:", "")
        if tid and tid not in tx_ex:
            tx_ex[tid] = []
            tx_cds[tid] = []
        if feat == "exon":
            tx_ex[tid].append((s, e))
        elif feat == "CDS":
            tx_cds[tid].append((s, e, int(phase), strand))
        for k in ("gene_id", "transcript_id", "ID"):
            if k in a and VER.search(a[k]):
                idver += 1
                break
    # phase audit: recompute expected phase from cumulative coding length
    for tid, segs in tx_cds.items():
        if not segs:
            continue
        strand = segs[0][3]
        cum = 0
        for s, e, ph, _ in sorted(segs, key=lambda x: x[0], reverse=(strand == "-")):
            exp = (3 - (cum % 3)) % 3
            if exp != ph:
                phase_bad += 1
            cum += e - s + 1
    return {"feats": feats, "seqids": seqids, "keys": sorted(keys), "idver": idver,
            "tx_ex": tx_ex, "tx_cds": tx_cds, "phase_bad": phase_bad, "col_ok": col_ok}


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for ln in open(path):
        ln = ln.strip()
        if ln.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = ln[1:].split()[0], []
        elif ln:
            buf.append(ln)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def read_fasta_list(path):
    """ordered list of (header, seq)"""
    out, name, buf = [], None, []
    for ln in open(path):
        ln = ln.strip()
        if ln.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = ln[1:].split()[0], []
        elif ln:
            buf.append(ln)
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def concat_by_tx(entries):
    """bedtools -name headers look like 'txid::chr1:100-200(+)'; group in file order"""
    out = {}
    for h, s in entries:
        txid = h.split("::")[0]
        out.setdefault(txid, "")
        out[txid] += s
    return out


def norm(x):
    x = VER.sub("", x)
    return re.sub(r"^(transcript|gene):", "", x)


R = {}
F = {n: read_gxf(os.path.join(OUT, n)) for n in
     ["annotation.gtf", "annotation.gff3", "annotation_chr1.gff3",
      "conv_from_gff3.gtf", "conv_from_gtf.gff3"]}

# -- structure: column count, seqid namespaces --
R["col9_all_rows"] = {n: f["col_ok"] for n, f in F.items()}
R["seqids"] = {n: sorted(f["seqids"]) for n, f in F.items()}
R["seqid_intersection_gtf_vs_gff3"] = sorted(
    F["annotation.gtf"]["seqids"] & F["annotation.gff3"]["seqids"])
R["naive_shared_features"] = sum(
    1 for ln in open(os.path.join(OUT, "annotation.gtf"))
    if not ln.startswith("#") and ln.split("\t")[0] in F["annotation.gff3"]["seqids"])

# -- attribute key inventory & ID versioning --
R["attr_keys"] = {n: f["keys"] for n, f in F.items()}
R["versioned_id_rows"] = {n: f["idver"] for n, f in F.items()}

# -- feature counts --
R["feature_counts"] = {n: f["feats"] for n, f in F.items()}

# -- phase audit --
R["phase_mismatches"] = {n: f["phase_bad"] for n, f in F.items()}

# -- conversion consistency: annotation.gtf vs conv_from_gff3.gtf --
a, b = F["annotation.gtf"], F["conv_from_gff3.gtf"]
ex_match = ex_total = cds_delta = 0
deltas = {}
for tid, exons in a["tx_ex"].items():
    t2 = norm(tid)
    hit = None
    for cand in b["tx_ex"]:
        if norm(cand) == t2:
            hit = cand
            break
    if hit is None:
        continue
    ex_total += 1
    if set(exons) == set(b["tx_ex"][hit]):
        ex_match += 1
    l1 = sum(e - s + 1 for s, e, _, _ in a["tx_cds"][tid]) if a["tx_cds"][tid] else 0
    l2 = sum(e - s + 1 for s, e, _, _ in b["tx_cds"][hit]) if b["tx_cds"][hit] else 0
    deltas[tid] = l2 - l1
R["exon_set_match"] = {"matched_tx": ex_total, "identical_exon_sets": ex_match}
R["cds_len_delta_gff3_minus_gtf"] = deltas
R["cds_delta_all_plus3"] = all(v == 3 for v in deltas.values()) and len(deltas) > 0

# roundtrip GTF->GFF3
c = F["conv_from_gtf.gff3"]
rt_tx = sum(1 for tid in a["tx_ex"] if any(norm(k) == norm(tid) for k in c["tx_ex"]))
rt_ex = sum(1 for tid in a["tx_ex"]
            if any(norm(k) == norm(tid) and set(v) == set(a["tx_ex"][tid])
                   for k, v in c["tx_ex"].items()))
R["roundtrip_gtf_to_gff3"] = {"tx_matched": rt_tx, "tx_exons_identical": rt_ex}

# -- FASTA comparisons --
tx_gtf = read_fasta(os.path.join(OUT, "tx_from_gtf.fa"))
tx_gff = read_fasta(os.path.join(OUT, "tx_from_gff3.fa"))
cds_gtf = read_fasta(os.path.join(OUT, "cds_from_gtf.fa"))
cds_gff = read_fasta(os.path.join(OUT, "cds_from_gff3.fa"))
prot_gtf = read_fasta(os.path.join(OUT, "prot_from_gtf.fa"))
prot_gff = read_fasta(os.path.join(OUT, "prot_from_gff3.fa"))
tx_asis = read_fasta(os.path.join(OUT, "tx_from_gff3_asis.fa"))
bedseqs = concat_by_tx(read_fasta_list(os.path.join(OUT, "exons_bed_sense.fa")))
overseqs = concat_by_tx(read_fasta_list(os.path.join(OUT, "exons_bed_overcor.fa")))

R["fasta_counts"] = {
    "tx_gtf": len(tx_gtf), "tx_gff3": len(tx_gff), "cds_gtf": len(cds_gtf),
    "cds_gff3": len(cds_gff), "prot_gtf": len(prot_gtf), "prot_gff3": len(prot_gff),
    "tx_gff3_asis": len(tx_asis), "bed_concat": len(bedseqs),
    "bed_concat_overcorrected": len(overseqs)}

# header namespace: BED names vs versioned GTF fasta headers
bed_names = set(bedseqs)
gtf_hdr = set(tx_gtf)
R["header_join"] = {
    "naive_intersection": len(bed_names & gtf_hdr),
    "after_strip_version": len({norm(x) for x in gtf_hdr} & bed_names)}

# per-transcript exact sequence checks
same_len = same_seq = 0
per_tx_len = {}
for name, s in bedseqs.items():
    hit = tx_gff.get(name) or tx_gff.get(name + ".1")
    for k in tx_gff:
        if norm(k) == name:
            hit = tx_gff[k]
    per_tx_len[name] = len(s)
    if hit is not None:
        same_len += (len(hit) == len(s))
        same_seq += (hit == s)
R["bed_vs_gffread_tx"] = {"n": len(bedseqs), "same_len": same_len, "same_seq_exact": same_seq}

over_delta = {}
for name, s in overseqs.items():
    over_delta[name] = per_tx_len[name] - len(s)
R["overcorrected_deficit"] = over_delta
R["overcorrected_total_bp_lost"] = sum(over_delta.values())

# CDS length delta in extracted FASTAs
cds_deltas = {}
prot_trailing_star = {"gtf": 0, "gff3": 0}
prot_identical = 0
for k, s in cds_gtf.items():
    hit = None
    for j in cds_gff:
        if norm(j) == norm(k):
            hit = cds_gff[j]
    if hit is not None:
        cds_deltas[norm(k)] = len(hit) - len(s)
    p1 = prot_gtf.get(k) or next((prot_gtf[j] for j in prot_gtf if norm(j) == norm(k)), None)
    p2 = prot_gff.get(k) or next((prot_gff[j] for j in prot_gff if norm(j) == norm(k)), None)
    if p1 is not None and p2 is not None:
        prot_identical += (p1 == p2)
        prot_trailing_star["gtf"] += p1.endswith("*")
        prot_trailing_star["gff3"] += p2.endswith("*")
R["cds_fasta_delta"] = cds_deltas
R["protein"] = {"identical": prot_identical, "trailing_star_gtf": prot_trailing_star["gtf"],
                "trailing_star_gff3": prot_trailing_star["gff3"], "n": len(prot_gtf)}
R["tx_len_stats"] = {
    "n_tx": len(tx_gtf),
    "median_tx_len": int(median(len(s) for s in tx_gtf.values())),
    "min_tx_len": min(len(s) for s in tx_gtf.values()),
    "max_tx_len": max(len(s) for s in tx_gtf.values()),
    "median_cds_len": int(median(len(s) for s in cds_gtf.values())),
    "median_prot_len": int(median(len(s) for s in prot_gtf.values()))}

# per-transcript table for figures
per_tx = {}
for k, s in tx_gtf.items():
    n = norm(k)
    c1 = next((cds_gtf[j] for j in cds_gtf if norm(j) == n), "")
    c2 = next((cds_gff[j] for j in cds_gff if norm(j) == n), "")
    src = next((tid for tid in a["tx_ex"] if norm(tid) == n), None)
    per_tx[n] = {
        "tx_len": len(s),
        "cds_len_gtf": len(c1), "cds_len_gff3": len(c2),
        "prot_len": len(prot_gtf.get(k, "")),
        "n_exons": len(a["tx_ex"][src]) if src else 0,
        "bed_deficit": over_delta.get(n, 0)}
R["per_tx"] = per_tx

json.dump(R, open(os.path.join(OUT, "audit.json"), "w"), indent=1)

L = []
A = L.append
A("== 049 audit report ==")
A("9 columns in every row: %s" % R["col9_all_rows"])
A("seqids: %s" % R["seqids"])
A("gtf-vs-gff3 seqid intersection: %s -> naive shared features = %d"
  % (R["seqid_intersection_gtf_vs_gff3"], R["naive_shared_features"]))
A("feature counts: %s" % json.dumps(R["feature_counts"]))
A("attr keys gtf: %s" % R["attr_keys"]["annotation.gtf"])
A("attr keys gff3: %s" % R["attr_keys"]["annotation.gff3"])
A("versioned-ID rows: %s" % R["versioned_id_rows"])
A("phase mismatches vs recomputed: %s" % R["phase_mismatches"])
A("exon sets gtf vs gff3->gtf: %d/%d identical" % (ex_match, ex_total))
A("CDS len delta (gff3 - gtf) all +3: %s; deltas: %s"
  % (R["cds_delta_all_plus3"], sorted(set(deltas.values()))))
A("roundtrip gtf->gff3: %s" % R["roundtrip_gtf_to_gff3"])
A("fasta counts: %s" % json.dumps(R["fasta_counts"]))
A("header join bed vs versioned gtf fasta: %s" % R["header_join"])
A("bed(start-1) concat vs gffread -w: %s" % R["bed_vs_gffread_tx"])
A("over-corrected bed deficit per tx (bp): total=%d across %d tx"
  % (R["overcorrected_total_bp_lost"], len(over_delta)))
A("CDS fasta delta per tx: %s" % sorted(set(cds_deltas.values())))
A("proteins: identical=%d/%d, trailing '*': gtf=%d gff3=%d"
  % (prot_identical, R["protein"]["n"], prot_trailing_star["gtf"], prot_trailing_star["gff3"]))
A("tx/cds/prot lengths: %s" % json.dumps(R["tx_len_stats"]))
open(os.path.join(OUT, "audit_report.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
