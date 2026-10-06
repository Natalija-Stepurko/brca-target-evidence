"""Stage `wes`: the SEQC2 HCC1395 / HCC1395BL exome pair from raw reads to scored somatic calls (DESIGN §13).

Steps, each skipped when its output exists:
  1. verify the FASTQ md5 sums ENA publishes; record every input in the manifest
  2. unpack the SEQC2 GRCh38.d1.vd1 reference and its BWA index (the truth set's own reference)
  3. pilot: align the first 1,000,000 read pairs of the tumour and time it
  4. align both samples (bwa mem, 8 threads), sort, mark duplicates (samtools markdup)
  5. Mutect2 tumour vs normal on the exome target (padded 100 bp), FilterMutectCalls
  6. score PASS calls against the SEQC2 high-confidence truth set inside target ∩ high-confidence regions
  7. map PASS calls to genes through the target BED's gene names; join HCC1395's own DepMap gene effects
     and the basal-like candidate list

Tools come from a micromamba environment with no root access: bwa 0.7.19, samtools 1.24, GATK 4.6.2,
bcftools, bedtools. The run is CPU-bound on 4 physical cores; timings are written to results/wes/timings.json.
BQSR, the germline resource and the contamination estimate are omitted and stated as omitted.
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from bte import config as C
from bte import provenance as P

OUT = C.RESULTS / "wes"
DATA = C.DATA / "wes"
WORK = Path(os.environ.get("BTE_WES_WORK", "/scratch/bte-wes/work"))
REF_DIR = Path(os.environ.get("BTE_WES_REF", "/scratch/bte-wes/ref"))
MAMBA = Path(os.environ.get("BTE_WES_MAMBA", "/scratch/bte-wes/bin/micromamba"))
MAMBA_ROOT = os.environ.get("MAMBA_ROOT_PREFIX", "/scratch/bte-wes/mamba")
THREADS = int(os.environ.get("BTE_WES_THREADS", "8"))
PAD = 100
PILOT_PAIRS = 1_000_000

SAMPLES = {"tumour": ("SRR7890850", "HCC1395"), "normal": ("SRR7890851", "HCC1395BL")}
ENA_MD5 = {"SRR7890850_1.fastq.gz": "40ad1990b670a07ac855236c63782585",
           "SRR7890850_2.fastq.gz": "d849aae1df29aa8d78038ebd5ee07792",
           "SRR7890851_1.fastq.gz": "872a47b05e8ebdc7cb7d2a6de4adedcd",
           "SRR7890851_2.fastq.gz": "1240d5e5d402fbd192437cec09388c83"}
SEQC2 = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG"
INPUTS = {
    "GRCh38.d1.vd1_BWA.tar.gz": f"{SEQC2}/technical/reference_genome/GRCh38/GRCh38.d1.vd1_BWA.tar.gz",
    "GRCh38.d1.vd1.fa.tar.gz": f"{SEQC2}/technical/reference_genome/GRCh38/GRCh38.d1.vd1.fa.tar.gz",
    "GRCh38.d1.vd1.fa.fai": f"{SEQC2}/technical/reference_genome/GRCh38/GRCh38.d1.vd1.fa.fai",
    "GRCh38.d1.vd1.dict": f"{SEQC2}/technical/reference_genome/GRCh38/GRCh38.d1.vd1.dict",
    "S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed6.gz":
        f"{SEQC2}/technical/reference_genome/Exome_Target_bed/S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed6.gz",
    "high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz": f"{SEQC2}/release/latest/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz",
    "high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz": f"{SEQC2}/release/latest/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz",
    "High-Confidence_Regions_v1.2.bed": f"{SEQC2}/release/latest/High-Confidence_Regions_v1.2.bed",
}
for _run in ("SRR7890850", "SRR7890851"):
    for _i in (1, 2):
        INPUTS[f"{_run}_{_i}.fastq.gz"] = (f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR789/00{_run[-1]}/{_run}/"
                                           f"{_run}_{_i}.fastq.gz")

TIMES = {}


def sh(cmd: str, log: str | None = None, env_tools=True):
    """Run a shell command inside the micromamba env; raise on failure; time it."""
    full = f"{MAMBA} run -n wes bash -o pipefail -c {json.dumps(cmd)}" if env_tools else cmd
    t0 = time.time()
    r = subprocess.run(full, shell=True, env={**os.environ, "MAMBA_ROOT_PREFIX": MAMBA_ROOT},
                       capture_output=True, text=True)
    if log:
        TIMES[log] = round(time.time() - t0, 1)
    if r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {cmd}\n{r.stderr[-3000:]}")
    return r.stdout


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while b := f.read(1 << 22):
            h.update(b)
    return h.hexdigest()


def step(name, outputs, fn):
    outputs = [Path(o) for o in outputs]
    if all(o.exists() for o in outputs):
        print(f"  skip    {name} (outputs exist)"); return
    print(f"  run     {name}", flush=True)
    t0 = time.time()
    fn()
    TIMES[name] = round(time.time() - t0, 1)
    print(f"          {TIMES[name]:,.0f}s", flush=True)


def run():
    t_all = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    REF_DIR.mkdir(parents=True, exist_ok=True)

    # 1. inputs
    for name, url in INPUTS.items():
        p = DATA / name
        assert p.exists(), f"missing input {p}; download it first (see data/wes/urls.txt)"
        if name in ENA_MD5:
            got = md5(p)
            assert got == ENA_MD5[name], f"md5 mismatch for {name}: {got}"
        P.record_file(p, url, "SEQC2 / ENA PRJNA489865")
    print(f"  inputs  {len(INPUTS)} files verified and recorded")

    # 2. reference
    fa = REF_DIR / "GRCh38.d1.vd1.fa"

    def unpack():
        sh(f"tar -xzf {DATA / 'GRCh38.d1.vd1.fa.tar.gz'} -C {REF_DIR}", env_tools=False)
        sh(f"tar -xzf {DATA / 'GRCh38.d1.vd1_BWA.tar.gz'} -C {REF_DIR}", env_tools=False)
        found = list(REF_DIR.rglob("GRCh38.d1.vd1.fa*"))
        for f in found:                                   # flatten whatever directory the tarballs used
            if f.parent != REF_DIR:
                shutil.move(str(f), REF_DIR / f.name)
        shutil.copy(DATA / "GRCh38.d1.vd1.fa.fai", REF_DIR / "GRCh38.d1.vd1.fa.fai")
        shutil.copy(DATA / "GRCh38.d1.vd1.dict", REF_DIR / "GRCh38.d1.vd1.dict")
    step("unpack reference", [fa, REF_DIR / "GRCh38.d1.vd1.fa.bwt"], unpack)

    # target and high-confidence regions
    target = WORK / "target.bed"
    target_pad = WORK / "target.pad100.bed"
    hc = WORK / "hc.bed"
    eval_bed = WORK / "eval.bed"

    def regions():
        # the lifted-over target carries alt-contig intervals; GRCh38.d1.vd1 is a no-alt build, so keep only
        # intervals on contigs the reference has
        sh(f"zcat {DATA / 'S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed6.gz'} | cut -f1-4 "
           f"| awk 'NR==FNR {{ok[$1]=1; next}} ok[$1]' {REF_DIR / 'GRCh38.d1.vd1.fa.fai'} - "
           f"| sort -k1,1 -k2,2n > {target}")
        sh(f"bedtools slop -i {target} -g {REF_DIR / 'GRCh38.d1.vd1.fa.fai'} -b {PAD} | bedtools merge > {target_pad}")
        sh(f"sort -k1,1 -k2,2n {DATA / 'High-Confidence_Regions_v1.2.bed'} | cut -f1-3 > {hc}")
        sh(f"bedtools intersect -a {target} -b {hc} | cut -f1-3 | bedtools merge > {eval_bed}")
    step("regions", [target, target_pad, hc, eval_bed], regions)

    # 3. pilot
    pilot = WORK / "pilot.bam"

    def do_pilot():
        r1, r2 = DATA / "SRR7890850_1.fastq.gz", DATA / "SRR7890850_2.fastq.gz"
        n = PILOT_PAIRS * 4
        t0 = time.time()
        sh(f"bwa mem -t {THREADS} -R '@RG\\tID:pilot\\tSM:pilot\\tPL:ILLUMINA' {fa} "
           f"<(zcat {r1} | head -n {n}) <(zcat {r2} | head -n {n}) | samtools sort -@ 2 -o {pilot} -")
        TIMES["pilot_seconds_per_million_pairs"] = round(time.time() - t0, 1)
    step("pilot 1M pairs", [pilot], do_pilot)
    if "pilot_seconds_per_million_pairs" in TIMES:
        pairs = {"tumour": 8.08e9 / 2 / 151, "normal": 10.54e9 / 2 / 151}
        est = sum(pairs.values()) / 1e6 * TIMES["pilot_seconds_per_million_pairs"] / 3600
        print(f"          estimated alignment for both samples: {est:.1f} h")

    # 4. align, sort, mark duplicates
    bams = {}
    for role, (run_id, sm) in SAMPLES.items():
        bam = WORK / f"{sm}.md.bam"
        bams[role] = bam

        def align(run_id=run_id, sm=sm, bam=bam):
            r1, r2 = DATA / f"{run_id}_1.fastq.gz", DATA / f"{run_id}_2.fastq.gz"
            raw = WORK / f"{sm}.namesorted.bam"
            sh(f"bwa mem -t {THREADS} -R '@RG\\tID:{run_id}\\tSM:{sm}\\tPL:ILLUMINA\\tLB:{sm}' {fa} {r1} {r2} "
               f"| samtools fixmate -@ 2 -m - - | samtools sort -@ 4 -m 2G -o {raw} -", log=f"align {sm}")
            sh(f"samtools markdup -@ {THREADS} -s {raw} {bam} 2> {WORK / f'{sm}.markdup.txt'} && samtools index {bam}",
               log=f"markdup {sm}")
            raw.unlink(missing_ok=True)
        step(f"align {sm}", [bam, Path(str(bam) + ".bai")], align)

    # coverage on target
    cov = OUT / "coverage.json"

    def coverage():
        out = {}
        for role, (_run_id, sm) in SAMPLES.items():
            txt = sh(f"samtools depth -a -b {target} {bams[role]} | awk '{{s+=$3; n++; if($3>=20) c++}} "
                     f"END {{printf \"%.1f %.4f\", s/n, c/n}}'")
            mean, frac20 = txt.split()
            flag = sh(f"samtools flagstat -@ 4 {bams[role]}")
            out[sm] = {"mean_target_depth": float(mean), "fraction_target_ge_20x": float(frac20),
                       "flagstat": flag.strip().splitlines()[:6]}
        cov.write_text(json.dumps(out, indent=1) + "\n")
    step("coverage", [cov], coverage)

    # 5. Mutect2
    raw_vcf, filt_vcf = WORK / "mutect2.unfiltered.vcf.gz", WORK / "mutect2.filtered.vcf.gz"

    def mutect():
        sh(f"gatk --java-options '-Xmx24g' Mutect2 -R {fa} -I {bams['tumour']} -I {bams['normal']} "
           f"-normal HCC1395BL -L {target_pad} --native-pair-hmm-threads {THREADS} "
           f"-O {raw_vcf} --f1r2-tar-gz {WORK / 'f1r2.tar.gz'} > {WORK / 'mutect2.log'} 2>&1", log="mutect2")
        sh(f"gatk LearnReadOrientationModel -I {WORK / 'f1r2.tar.gz'} -O {WORK / 'rom.tar.gz'}", log="orientation")
        sh(f"gatk FilterMutectCalls -R {fa} -V {raw_vcf} --ob-priors {WORK / 'rom.tar.gz'} -O {filt_vcf} "
           f"> {WORK / 'filter.log'} 2>&1", log="filter")
    step("mutect2", [filt_vcf], mutect)

    # 6. benchmark against the truth set inside target ∩ high-confidence regions
    bench = OUT / "benchmark.json"

    def benchmark():
        calls = WORK / "calls.pass.norm.vcf.gz"
        sh(f"bcftools view -f PASS -R {eval_bed} {filt_vcf} | bcftools norm -f {fa} -m -both -Oz -o {calls} "
           f"&& bcftools index -f {calls}")
        res = {}
        for kind, truth_name in (("snv", "high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"),
                                 ("indel", "high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz")):
            truth = WORK / f"truth.{kind}.vcf.gz"
            sh(f"bcftools view -R {eval_bed} {DATA / truth_name} | bcftools norm -f {fa} -m -both -Oz -o {truth} "
               f"&& bcftools index -f {truth}")
            typ = "snps" if kind == "snv" else "indels"
            d = WORK / f"isec.{kind}"
            sh(f"bcftools isec -v {typ} -p {d} {calls} {truth}")
            n = {k: int(sh(f"grep -vc '^#' {d}/{f}.vcf || true").strip() or 0)
                 for k, f in (("calls_only", "0000"), ("truth_only", "0001"), ("shared", "0002"))}
            tp, fp, fn = n["shared"], n["calls_only"], n["truth_only"]
            res[kind] = {**n, "precision": tp / max(1, tp + fp), "recall": tp / max(1, tp + fn),
                         "f1": 2 * tp / max(1, 2 * tp + fp + fn)}
        res["eval_region_bp"] = int(sh(f"awk '{{s+=$3-$2}} END {{print s}}' {eval_bed}").strip())
        res["note"] = ("PASS Mutect2 calls vs SEQC2 v1.2.1 high-confidence somatic calls, both restricted to the "
                       "exome target intersected with the high-confidence regions; BQSR, germline resource and "
                       "contamination estimate omitted")
        bench.write_text(json.dumps(res, indent=1) + "\n")
        print("          " + json.dumps({k: {kk: round(vv, 3) if isinstance(vv, float) else vv
                                              for kk, vv in v.items()} for k, v in res.items() if k in ("snv", "indel")}))
    step("benchmark", [bench], benchmark)

    # 7. genes with PASS somatic calls -> HCC1395's own dependencies and the basal-like candidates
    genes_out = OUT / "mutated_genes.csv"

    def genes():
        calls = WORK / "calls.pass.norm.vcf.gz"
        txt = sh(f"bcftools query -f '%CHROM\\t%POS0\\t%END\\t%REF\\t%ALT\\n' {calls} "
                 f"| bedtools intersect -a - -b {target} -wa -wb")
        rows = []
        for line in txt.strip().splitlines():
            f = line.split("\t")
            names = f[8] if len(f) > 8 else ""
            syms = sorted({x.split("|")[1] for x in names.split(",") if x.startswith("ref|")})
            for g in syms:
                rows.append({"chrom": f[0], "pos": int(f[1]) + 1, "ref": f[3], "alt": f[4], "gene": g})
        v = pd.DataFrame(rows).drop_duplicates()
        ge = pd.read_csv(C.DATA / "raw" / "depmap_24Q4_CRISPRGeneEffect.csv", index_col=0)
        ge = ge.loc[C.HCC1395_DEPMAP_ID]
        ge.index = [c.split(" (")[0] for c in ge.index]
        per_gene = v.groupby("gene").size().rename("n_variants").to_frame()
        per_gene["hcc1395_gene_effect"] = ge.reindex(per_gene.index)
        lists = pd.read_csv(C.RESULTS / "nominate" / "lists_discovery.csv")
        basal = lists[(lists.k == C.K) & (lists.subtype == "basal")]
        for arm, d in basal.groupby("arm"):
            per_gene[f"in_basal_top50_{arm}"] = per_gene.index.isin(d.gene)
        t1 = pd.read_csv(C.RESULTS / "truth" / "t1_gene_effect.csv", index_col=0)
        per_gene["mean_effect_basal_lines"] = t1.mean_basal.reindex(per_gene.index)
        per_gene["common_essential"] = t1.common_essential.reindex(per_gene.index)
        per_gene.sort_values("hcc1395_gene_effect").to_csv(genes_out)
        v.to_csv(OUT / "pass_variants_in_target.csv", index=False)
        dep = per_gene[per_gene.hcc1395_gene_effect <= C.HIT_GENE_EFFECT]
        print(f"          {len(v)} PASS variants in {len(per_gene)} target genes; "
              f"{len(dep)} of those genes are HCC1395 dependencies (effect <= {C.HIT_GENE_EFFECT})")
    step("genes", [genes_out], genes)

    (OUT / "timings.json").write_text(json.dumps({**TIMES, "threads": THREADS}, indent=1) + "\n")
    P.log_run("wes", {"threads": THREADS, "pad": PAD, "pilot_pairs": PILOT_PAIRS,
                      "tools": "bwa 0.7.19, samtools 1.24, GATK 4.6.2.0, bcftools, bedtools (micromamba, no root)"},
              [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*"))], time.time() - t_all)
    print(f"done in {(time.time() - t_all) / 3600:.1f} h")


if __name__ == "__main__":
    sys.exit(run())
