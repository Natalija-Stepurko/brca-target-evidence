# Does the protein layer pick better drug targets?

A pre-registered evaluation of multi-omics target nomination in breast cancer. The second study in a
series; the first, [single-cell cell states and survival](https://github.com/Natalija-Stepurko/single-cell-fm-probing),
asked whether foundation-model cell states predict patient survival.

Study 3, [brca-spatial-evidence](https://github.com/Natalija-Stepurko/brca-spatial-evidence), takes this
study's candidates to spatial transcriptomics sections and asks where they are expressed and whether
histology can see them.

**Status (2026-10-06):** design pre-registered in [`docs/DESIGN.md`](docs/DESIGN.md) before any data were
downloaded. No results yet. Project page: <https://natalija-stepurko.github.io/brca-target-evidence/>.

## The question

Three proteogenomic studies have nominated breast-cancer drug targets by combining DNA, RNA and protein
measurements from the same tumours (Mertins 2016, Krug 2020, Savage 2024). Whether the extra layers
improve the shortlist has never been measured: no study compares integrated with RNA-only nomination
at equal list size, scores both against evidence that played no part in the nomination, and controls
for how well studied each gene is. The two CPTAC breast cohorts do not even agree with each other.

> In breast cancer, does adding protein (and DNA) evidence to RNA nominate drug targets that are more
> often confirmed by held-out evidence than RNA alone — beyond what random genes matched for popularity
> achieve — and do the nominations replicate in an independent cohort?

A second question follows: when an AI agent writes the rationale for a nominated target, how many of its
claims survive a check against the held-out data, and does it prefer famous genes?

## What will be done

1. **Discovery cohort.** The 77 TCGA-BRCA tumours that CPTAC profiled (Mertins 2016): somatic mutations,
   copy number, RNA-seq and the mass-spectrometry proteome on the same tumours, with PAM50 subtypes.
2. **Six nomination arms, one gene universe, K = 50 per subtype.** RNA alone (R); RNA + protein (R+P);
   RNA + DNA (R+D); all three (R+P+D); MOFA+ on the three views; and RNA alone on all ~1,000 TCGA-BRCA
   tumours (R-all), the baseline for the practical question of what a programme can buy.
3. **Held-out truth no arm sees.** Subtype-matched CRISPR dependency in the 53 breast cell lines of
   DepMap Public 24Q4 (primary; continuous endpoint), and clinical precedence from Open Targets 26.09
   (secondary).
4. **The ladder.** A floor of 1,000 random lists matched on PubMed count, expression and protein
   detectability; a subtype-label permutation null; paired bootstrap intervals on arm differences;
   positive controls (ERBB2, ESR1); replication of the frozen pipeline in the 122-tumour Krug 2020 cohort.
5. **Predictions P1–P8**, with a smallest effect of interest, fixed in `docs/DESIGN.md` §8.
6. **Candidate dossiers.** The top genes of the best arm per subtype, each with every rung: dependency
   and which lines, tractability, normal-tissue expression, tumour-cell vs stroma expression in the
   50,002-cell atlas from study 1, survival association, replication status.
7. **Agent audit.** An agent writes dossiers for 30 candidates and 30 matched decoys without knowing
   which is which; every typed claim is checked against the frozen tables and labelled supported,
   contradicted or unsupported. Reported: faithfulness, decoy discrimination, correlation with PubMed count.
8. **Exome module.** The SEQC2 HCC1395 / HCC1395BL exome pair from raw reads to somatic calls, scored
   against the published truth set, then traced into HCC1395's own DepMap dependencies.

## Pipeline

```bash
export UV_PROJECT_ENVIRONMENT=.venv      # any path; the project is a uv package
uv sync --extra mofa --group dev         # core dependencies, MOFA+, test tools
uv run bte stages                        # the ten stages
uv run bte run data                      # downloads ~1.5 GB with sha256 provenance; builds the gene universe
uv run bte run scores && uv run bte run nominate && uv run bte run truth
uv run bte run ladder                    # ~3 h on 4 cores: 1,000 matched lists, 1,000 permutations, 1,000 bootstraps
uv run bte run replicate && uv run bte run dossier && uv run bte run report
uv run python docs/site/build.py         # the page, every number read from results/
```

| Stage | Module | Writes | Does |
|---|---|---|---|
| `data` | `stages/data.py` | `results/data/`, `results/MANIFEST.sha256` | downloads every input (`sources.py`), records sha256, builds the 8,203-gene universe |
| `scores` | `stages/scores.py` | `results/scores/` | per-layer Hedges' *g* per subtype for the discovery, TCGA-all and Krug cohorts; the tumour-vs-normal gate |
| `nominate` | `stages/nominate.py`, `arms.py` | `results/nominate/` | the six arms' top-K lists; MOFA+ fit; positive controls |
| `truth` | `stages/truth.py` | `results/truth/` | DepMap 24Q4 subtype-matched dependency; Open Targets 26.09 clinical precedence; tractability; GTEx |
| `ladder` | `stages/ladder.py` | `results/ladder/` | matched-random floor, label-permutation null, paired bootstrap, arm differences |
| `replicate` | `stages/replicate.py` | `results/replicate/` | the frozen pipeline on Krug 2020; overlap against a permutation null |
| `dossier` | `stages/dossier.py` | `results/dossier/` | every rung for the best arm's candidates, including the study-1 atlas and TCGA survival |
| `wes` | `stages/wes.py` | `results/wes/` | SEQC2 exome pair from raw reads to scored somatic calls (needs the micromamba environment below) |
| `audit` | — | `results/audit/` | agent-written dossiers checked claim by claim (to come) |
| `report` | `stages/report.py` | `results/report/` | verdicts on P1–P6 from the tables, and the figures |

`tests/test_design_consistency.py` checks every number in `src/bte/config.py` against the text of
`docs/DESIGN.md`, so the code cannot drift from the design unnoticed. CI runs ruff, pytest and the page build.

**Exome module tools, with no root access.** `bwa`, `samtools`, `gatk4`, `bcftools` and `bedtools` from
bioconda in a micromamba environment:

```bash
curl -sL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /scratch/bte-wes bin/micromamba
export MAMBA_ROOT_PREFIX=/scratch/bte-wes/mamba
/scratch/bte-wes/bin/micromamba create -n wes -c conda-forge -c bioconda bwa samtools gatk4 bcftools bedtools
# inputs: docs/wes_inputs.txt lists the SEQC2 reference, target BED, truth set and ENA FASTQ files (~15 GB);
# download them into data/wes/
uv run bte run wes                       # ~6 h on 4 cores; resumable step by step
```

## Repository

| Path | What |
|---|---|
| [`docs/DESIGN.md`](docs/DESIGN.md) | the pre-registered design; §12 logs every change made after first contact with the data |
| [`research/literature.md`](research/literature.md) | the survey behind the design, with a strength tag on every entry |
| [`docs/site/build.py`](docs/site/build.py) → `docs/index.html` | the project page; every number is read from `results/` when it is built |
| `src/bte/` | the package: `config.py` (pre-registered parameters), `sources.py` (inputs), `cohorts.py`, `arms.py`, `stages/` |
| `results/` | small tables and figures from every stage (tracked); downloaded data and large intermediates are not |

## Data and licences

| Dataset | Route | Licence |
|---|---|---|
| TCGA-BRCA mutations, copy number, RNA-seq, PAM50, survival | UCSC Xena | open; the results here are based upon data generated by the TCGA Research Network |
| CPTAC breast proteome, Mertins 2016 (77 tumours, TCGA barcodes) | PDC000173 / cBioPortal | CC BY |
| CPTAC breast cohort, Krug 2020 (122 tumours) | LinkedOmics / cBioPortal | CC BY |
| DepMap Public 24Q4 CRISPR gene effect and model annotations | figshare 27993248 | CC BY 4.0 — the last release under it; later releases carry non-commercial terms, so 24Q4 is pinned |
| Open Targets Platform 26.09 clinical datasets, tractability | EBI FTP (parquet) | CC0; ChEMBL-derived fields CC BY-SA 3.0 |
| GTEx v10 median gene TPM | GTEx portal | open |
| NCBI gene2pubmed | NCBI FTP | open |
| SEQC2 HCC1395 / HCC1395BL exomes; somatic truth set v1.2.1 | ENA PRJNA489865; NCBI FTP | open |

Every downloaded file is recorded with URL, release and sha256 in `results/MANIFEST.sha256`. No data
from Sanger Project Score or Cell Model Passports are used.

## Licence

MIT — see [`LICENSE`](LICENSE).
