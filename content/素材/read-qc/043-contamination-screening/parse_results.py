#!/usr/bin/env python3
"""043 contamination-screening: parse kraken2 reports and FastQ Screen reports
into results.json for figure generation."""
import json
import os

BASE = os.path.dirname(os.path.abspath(__file__))
SAMPLES = ["S1", "S2", "S3"]
CONFS = ["0.0", "0.05", "0.1", "0.2"]


def parse_design():
    design = {s: {} for s in SAMPLES}
    with open(os.path.join(BASE, "design.tsv")) as fh:
        next(fh)
        for line in fh:
            s, taxon, pairs, pct = line.strip().split("\t")
            design[s][taxon] = float(pct)
    return design


def parse_kreport(path):
    out = {"unclassified_pct": 0.0, "species": {}, "n_species": 0}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            pct, rank, name = float(f[0]), f[3], f[5].strip()
            if rank == "U":
                out["unclassified_pct"] = pct
            elif rank == "S" and pct > 0:
                out["species"][name] = pct
                out["n_species"] += 1
    return out


def parse_screen_txt(path):
    """Parse FastQ Screen *_screen.txt; adapt to actual header layout."""
    rows, header = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if header is None and "Genome" in f[0] and any("hit" in c.lower() for c in f[1:]):
                header = [c.strip() for c in f]
                continue
            if header is None:
                continue
            name = f[0].strip()
            rows[name] = {header[i]: float(f[i]) for i in range(1, len(f))}
    return rows


def main():
    results = {"design": parse_design(), "kraken": {}, "screen": {}, "conf_sweep": {}}

    for s in SAMPLES:
        results["kraken"][s] = parse_kreport(os.path.join(BASE, "kraken2", f"{s}.kreport"))

    for s in SAMPLES:
        # find this sample's screen txt (fastq_screen writes <base>_screen.txt)
        path = os.path.join(BASE, "screen", f"{s}_1_screen.txt")
        if not os.path.exists(path):
            cands = [x for x in os.listdir(os.path.join(BASE, "screen")) if x.startswith(s)]
            path = os.path.join(BASE, "screen", sorted(cands)[0])
        results["screen"][s] = parse_screen_txt(path)

    for s in ["S1", "S2"]:
        results["conf_sweep"][s] = {}
        for c in CONFS:
            results["conf_sweep"][s][c] = parse_kreport(
                os.path.join(BASE, "kraken2", f"{s}.conf{c}.kreport"))

    with open(os.path.join(BASE, "results.json"), "w") as fh:
        json.dump(results, fh, indent=2)

    # console summary
    print("== design vs kraken2 (species-level % of read pairs)")
    for s in SAMPLES:
        k = results["kraken"][s]
        row = {"unclassified": k["unclassified_pct"]}
        row.update({n.split(",")[0]: p for n, p in k["species"].items()})
        print(f"{s}: design={results['design'][s]} kraken={row}")
    print("== FastQ Screen (% one_hit_one_genome per genome)")
    for s in SAMPLES:
        print(s, {k: v for k, v in results["screen"][s].items()})
    print("== confidence sweep")
    for s in ["S1", "S2"]:
        for c in CONFS:
            r = results["conf_sweep"][s][c]
            print(f"{s} conf={c}: unclassified={r['unclassified_pct']}% "
                  f"n_species={r['n_species']} species={r['species']}")
    print("results.json written")


if __name__ == "__main__":
    main()
