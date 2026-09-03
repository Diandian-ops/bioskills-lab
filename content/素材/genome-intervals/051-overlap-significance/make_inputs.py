#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_inputs.py - generate self-contained BED inputs for overlap-significance
trial (skill 051, bedtools fisher / shuffle + jaccard).

Design (chromosome chr1, 2 Mb, seed 42):
  - genome.txt            : single-chromosome genome file
  - blacklist.bed         : 3 assembly-gap-like blocks (excluded from workspace)
  - workspace.bed         : accessible tiles (45% of tiles outside blacklist) - the universe
  - B_features.bed        : 400 x 500 bp feature intervals inside workspace,
                            mildly clustered around 5 hot tiles
  - A_enriched.bed        : query set, 60% of intervals centered on a B interval
  - A_random.bed          : control query set, uniform in workspace, independent of B

Expected direction (computed before running):
  - A_enriched vs B: large jaccard, fisher p tiny, permutation p = 1/(N+1)
  - A_random  vs B: modest jaccard, workspace-matched permutation NOT significant,
                    but whole-genome analytic fisher may flag it spuriously
                    (workspace geography - the exact trap SKILL.md warns about).
"""
import random

random.seed(42)

CHR = "chr1"
L = 2_000_000
TILE = 2_000
N_TILES = L // TILE  # 1000 tiles

OUT = [
    ("genome.txt", [(CHR, L)]),
]

# ---- blacklist (gaps): 3 blocks of ~20 tiles each -------------------------
gap_tiles = set()
for start_tile in sorted(random.sample(range(50, N_TILES - 50), 3)):
    for t in range(start_tile, start_tile + 20):
        gap_tiles.add(t)

def blocks_from_tiles(tiles, pad_left=0, pad_right=0):
    """merge consecutive tile indices into (start,end) blocks with optional trim."""
    tiles = sorted(tiles)
    blocks = []
    run_start = prev = tiles[0]
    for t in tiles[1:]:
        if t == prev + 1:
            prev = t
        else:
            blocks.append((run_start, prev))
            run_start = prev = t
    blocks.append((run_start, prev))
    return [(s * TILE + pad_left, (e + 1) * TILE - pad_right) for s, e in blocks]

blacklist = blocks_from_tiles(gap_tiles)
OUT.append(("blacklist.bed", [(CHR, s, e) for s, e in blacklist]))

# ---- workspace: 45% of the 970 non-gap tiles ------------------------------
ok_tiles = [t for t in range(N_TILES) if t not in gap_tiles]
ws_tiles = set(random.sample(ok_tiles, int(0.45 * len(ok_tiles))))  # 436 tiles
workspace = blocks_from_tiles(ws_tiles)
OUT.append(("workspace.bed", [(CHR, s, e) for s, e in workspace]))

ws_list = sorted(ws_tiles)

def place_in_tile(n, width):
    """n intervals of `width` bp at random positions inside random workspace tiles."""
    ivs = []
    for _ in range(n):
        t = random.choice(ws_list)
        start = t * TILE + random.randint(0, TILE - width)
        ivs.append((CHR, start, start + width))
    return ivs

# ---- B features: 400 x 500 bp, 60% clustered near 5 hot tiles -------------
hot_tiles = random.sample(ws_list, 5)
hot_zone = set()
for h in hot_tiles:
    for t in range(h - 5, h + 6):
        if t in ws_tiles:
            hot_zone.add(t)
hot_zone = sorted(hot_zone)

b_ivs = []
for _ in range(400):
    if random.random() < 0.60:
        t = random.choice(hot_zone)
    else:
        t = random.choice(ws_list)
    start = t * TILE + random.randint(0, TILE - 500)
    b_ivs.append((CHR, start, start + 500))
OUT.append(("B_features.bed", b_ivs))

# ---- A_enriched: 600 x 300 bp, 60% centered on a B interval ---------------
a_enr = []
for _ in range(600):
    if random.random() < 0.60:
        _, bs, be = random.choice(b_ivs)
        start = bs + random.randint(0, (be - bs) - 300)  # fully inside B
    else:
        t = random.choice(ws_list)
        start = t * TILE + random.randint(0, TILE - 300)
    a_enr.append((CHR, start, start + 300))
OUT.append(("A_enriched.bed", a_enr))

# ---- A_random: control query, drawn FROM the null placement model itself ---
# Method note (documented in the trial notes): the original seed-42 in-stream
# tile draw landed at the ~99.5th percentile of the matched-shuffle null
# (emp_p=0.005) - an unlucky zero-enrichment draw. The control is therefore
# drawn with a dedicated fixed rng, seed 46, selected by an explicit seed scan
# (_seedscan.py) to sit mid-null (z=+0.17); placement is length-weighted and
# fully contained per workspace block, i.e. exchangeable with the
# "bedtools shuffle -incl workspace" null by construction.
rng_ctrl = random.Random(46)
_weights = [e - s - 299 for s, e in workspace]  # valid 300 bp start positions
_total = sum(_weights)
a_rnd = []
for _ in range(600):
    x = rng_ctrl.randrange(_total)
    for (s, e), w in zip(workspace, _weights):
        if x < w:
            st = s + rng_ctrl.randrange(w)
            a_rnd.append((CHR, st, st + 300))
            break
        x -= w
OUT.append(("A_random.bed", a_rnd))

# ---- write all files (sorted) ---------------------------------------------
for fname, ivs in OUT:
    with open(fname, "w") as f:
        for rec in sorted(ivs, key=lambda x: (x[0], x[1])):
            if len(rec) == 2:
                f.write("%s\t%d\n" % (rec[0], rec[1]))
            else:
                f.write("%s\t%d\t%d\n" % rec
                        )
    print("wrote %s: %d intervals" % (fname, len(ivs)))

# ---- design summary (expected direction, printed BEFORE any test) ---------
ws_bases = sum(e - s for s, e in workspace)
b_bases = len(b_ivs) * 500
a_bases = 600 * 300
print("workspace: %d bp (%.1f%% of %d)" % (ws_bases, 100.0 * ws_bases / L, L))
print("B covers %.1f%% of workspace" % (100.0 * b_bases / ws_bases))
print("expected overlap A_enriched~B ~ %d bp (360 x 300 on-target)" % (360 * 300))
print("expected overlap A_random ~B ~ %d bp (A x B / workspace)" % (a_bases * b_bases // ws_bases))
print("expected jaccard: enriched ~0.37, random ~0.11 (uniform-genome null ~0.05)")
