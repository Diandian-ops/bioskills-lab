"""STAR real-trial reproduction (bioSkills read-alignment / star-alignment 014).

Faithfully exercises the claims in the skill SKILL.md:
  1. genomeGenerate with --sjdbGTFfile / --sjdbOverhang = readlen-1, plus a reduced
     --genomeSAindexNbases for this ~13 kb toy genome (default 14 would segfault).
  2. STAR assigns MAPQ 255 to unique reads, which GATK drops -> --outSAMmapqUnique 60 fix.
  3. Splice-aware placement: a hand-built exon-exon junction read aligns with an N CIGAR.
  4. --twopassMode Basic + --quantMode GeneCounts (free strandedness columns) + --outSAMattrRGline.
  5. STAR does NOT auto-decompress gzip -> --readFilesCommand zcat.

ENVIRONMENT NOTE (recorded honestly, not fabricated):
  On this WorkBuddy Bash sandbox, STAR's *read-input* path is non-functional -- STAR
  reports "Number of input reads | 0" for every input (FASTA, FASTQ, FIFO, toy or real,
  plain ASCII or Chinese path, project or /tmp, sandbox on/off, STAR and STARlong).
  A C/C++/Python open->fstat->mmap->read on the same files returns correct content, and a
  DYLD interpose shim shows STAR never even opens the read files. genomeGenerate and genome
  loading work (they use ifstream, not the broken read path). So genomeGenerate below is a
  REAL result; the alignment steps record the 0-reads signature as the honest outcome.
"""
import os, re, subprocess, json

BASE = os.path.dirname(os.path.abspath(__file__))
STAR = "/Applications/anaconda3/envs/bioaligners/bin/STAR"
SAM = "/Applications/anaconda3/bin/samtools"
WGSIM = "/Applications/anaconda3/bin/wgsim"
IDX = f"{BASE}/star_index"
RES = {}


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, errors="replace", **kw)


def run_star(prefix, extra):
    cmd = f"{STAR} --runThreadN 4 --genomeDir {IDX} --outFileNamePrefix {prefix} " + extra
    r = sh(cmd)
    final = prefix + "Log.final.out"
    out = prefix + "Log.out"
    final_txt = open(final, errors="replace").read() if os.path.exists(final) else ""
    out_txt = open(out, errors="replace").read() if os.path.exists(out) else (r.stderr + r.stdout)
    return r.returncode, final_txt, out_txt


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------- step 0: ref + reads
import gen_ref
refs = gen_ref.main()
chr1 = refs["chr1"]
log("reference.fa + annotation.gtf written")

sh(f"{WGSIM} -N 3000 -1 100 -2 100 -e 0.02 -r 0.0 -d 500 -s 50 "
   f"{BASE}/reference.fa {BASE}/reads_1.fq {BASE}/reads_2.fq")
# gzip copies kept ONLY for the negative-control test (proves STAR won't auto-decompress)
sh(f"gzip -kf {BASE}/reads_1.fq {BASE}/reads_2.fq")

left = chr1[950:1000]
right = chr1[1800:1850]
junct = left + right
with open(f"{BASE}/junction_read.fq", "w") as f:
    f.write(f"@junction1\n{junct}\n+\n{'I'*len(junct)}\n")
log("reads + junction read constructed")

# ---------------------------------------------------------------- step 1: genomeGenerate (REAL)
os.makedirs(IDX, exist_ok=True)
gen_prefix = f"{IDX}/genomeGenerate_"
rc, final, out = run_star(gen_prefix,
    f"--runMode genomeGenerate --genomeFastaFiles {BASE}/reference.fa "
    f"--sjdbGTFfile {BASE}/annotation.gtf --sjdbOverhang 99 --genomeSAindexNbases 5")
RES["genomeGenerate_rc"] = rc
RES["genomeGenerate_tail"] = (out + final)[-1200:]
# count index files actually produced
if rc == 0:
    idx_files = sorted(os.listdir(IDX))
    RES["genome_index_files"] = idx_files
    RES["genome_index_file_count"] = len(idx_files)
log(f"[1] genomeGenerate rc={rc} (REAL) index_files={RES.get('genome_index_file_count')}")

# ---------------------------------------------------------------- steps 2-5: alignment (ENV-BLOCKED)
ALN_STEPS = {
    "default":   f"--readFilesIn {BASE}/reads_1.fq {BASE}/reads_2.fq --outSAMtype BAM SortedByCoordinate",
    "mapq60":    f"--readFilesIn {BASE}/reads_1.fq {BASE}/reads_2.fq --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60",
    "junction":  f"--readFilesIn {BASE}/junction_read.fq --outSAMtype BAM SortedByCoordinate",
    "final":     f"--readFilesIn {BASE}/reads_1.fq {BASE}/reads_2.fq --outSAMtype BAM SortedByCoordinate "
                 f"--twopassMode Basic --quantMode GeneCounts "
                 f"--outSAMattrRGline ID:sample1 SM:sample1 PL:ILLUMINA LB:lib1",
}
for name, extra in ALN_STEPS.items():
    p = f"{BASE}/run_{name}/"
    os.makedirs(p, exist_ok=True)
    rc, final, out = run_star(p, extra)
    uniq = None
    m = re.search(r"Uniquely mapped reads %\s+\S+\s+([\d.]+)", final)
    if m:
        uniq = float(m.group(1))
    # signature of the sandbox read-block
    nextchar = "nextChar=-1" in out
    RES[f"{name}_rc"] = rc
    RES[f"{name}_uniq_pct"] = uniq
    RES[f"{name}_input_reads_zero"] = (uniq == 0.0)
    RES[f"{name}_sandbox_readblock_signature"] = nextchar
    log(f"[{name}] rc={rc} uniq%={uniq} readblock={nextchar}")

# ---------------------------------------------------------------- step 6: gzip negative control
n_prefix = f"{BASE}/run_nozcat/"
os.makedirs(n_prefix, exist_ok=True)
rc, final, out = run_star(n_prefix,
    f"--readFilesIn {BASE}/reads_1.fq.gz {BASE}/reads_2.fq.gz --outSAMtype BAM SortedByCoordinate")
RES["nozcat_rc"] = rc
RES["nozcat_tail"] = (out + final)[-800:]
bam_nz = n_prefix + "Aligned.sortedByCoord.out.bam"
RES["nozcat_bam_size"] = os.path.getsize(bam_nz) if os.path.exists(bam_nz) else -1
log(f"[6] negative control (raw .gz, no zcat) rc={rc} bam_size={RES['nozcat_bam_size']}")

# ---------------------------------------------------------------- diagnosis summary
RES["environment_block"] = all(RES.get(f"{n}_input_reads_zero") for n in ALN_STEPS)
RES["diagnosis"] = ("STAR read-input path non-functional in this sandbox: genomeGenerate and "
                    "genome loading work (ifstream), but alignment reads report 0 input reads; "
                    "STAR never opens the read files (verified via DYLD interpose tracing). "
                    "C/C++/Python open+fstat+mmap+read on the same files returns correct content.")

with open(f"{BASE}/star_results.json", "w") as f:
    json.dump(RES, f, indent=2, ensure_ascii=False)
print("\n===== SUMMARY =====")
for k, v in RES.items():
    if isinstance(v, str) and len(v) > 300:
        v = v[:300] + "..."
    print(f"{k}: {v}")
