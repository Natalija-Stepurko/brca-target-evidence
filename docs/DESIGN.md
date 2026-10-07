# Design — does the protein layer pick better drug targets?

*This document is the experimental design: what is measured, against what, and in what order. The
related work is in [`../research/literature.md`](../research/literature.md). Sections 1–9 are
committed before any dataset is downloaded for this study (the commit that merges this file into
`main` is what "pre-specified" refers to; there is no external registry). Any later change to §1–§9 is
listed in §12 with the date, the commit and the reason; nothing in them is edited silently. §11
("what was known before this commit") records the feasibility checks that preceded the design.*

## 1. The question

Target identification from patient tumours now routinely combines DNA, RNA and protein measurements
on the same samples. Three proteogenomic studies have nominated breast-cancer targets this way
(Mertins 2016, Krug 2020, Savage 2024; `literature.md` §1), and the field's working assumption is that
each added layer improves the shortlist. That assumption has not been tested. No published study
compares an integrated nomination with an RNA-only nomination at the same list size, scores both
against evidence that played no part in the nomination, and controls for how well studied each gene
is. The two CPTAC breast cohorts do not even agree with each other: Krug 2020 reports that its
candidate drivers do not overlap with Mertins 2016's, and nobody has measured how reproducible such
nominations are.

This study asks, in breast cancer:

> **Does adding protein (and DNA) evidence to RNA nominate drug targets that are more often confirmed
> by held-out evidence than RNA alone — beyond what random genes matched for popularity achieve — and
> do the nominations replicate in an independent cohort?**

It also asks a second question that the literature leaves open (`literature.md` §7): when an AI agent
writes the rationale for a nominated target, **how many of its claims survive a check against the
held-out data, and does it prefer famous genes?**

**Indication: breast cancer.** It has the only two public tumour cohorts with matched DNA, RNA and
mass-spectrometry proteomics on the same samples (both CPTAC), a large RNA-only cohort to measure what
RNA alone can do with more patients (TCGA-BRCA, ~1,000 tumours), 53 CRISPR-screened cell lines in
DepMap for held-out dependency evidence, and well-known approved targets (ERBB2, ESR1, CDK4/6) that
serve as positive controls.

## 2. Why this is not another target list

A ranked list of genes means nothing on its own. Every list produced here is read against a ladder:

| Rung | What it is | What it establishes |
|---|---|---|
| **Matched-random floor** | 1,000 random gene lists of the same size, each gene matched to a nominated gene on publication count, expression level and protein detectability | what any list of that shape scores — the honest zero, with the popularity bias removed |
| **Label-permutation null** | subtype labels shuffled across tumours, every arm re-run | what the whole pipeline produces when there is no subtype signal |
| **Baseline arm** | RNA alone, same tumours, same gene universe, same list size | the thing every other arm has to beat |
| **Positive controls** | ERBB2 in the HER2-enriched list, ESR1 in the luminal list | if the pipeline cannot find these, it is wrong |
| **Replication** | the frozen pipeline run on the second CPTAC cohort | whether the nominations are a property of breast cancer or of one cohort |

The expected answer, stated here before any data: **little or no gain from the protein layer at the list
level**, because the genes that are druggable tend to have the higher mRNA–protein correlation
(Savage 2024), with gene-specific exceptions. A clean null is informative: it tells a target-discovery
team when proteomics is worth paying for. Either answer is reported.

## 3. Cohorts

| Role | Cohort | n | Layers | Source |
|---|---|---|---|---|
| **Discovery** | TCGA-BRCA tumours profiled by CPTAC (Mertins 2016) | 105 profiled, **77 passed CPTAC's quality control**; the 77 are used | DNA (masked somatic mutations, GISTIC2 copy number), RNA-seq, proteome, PAM50 | TCGA via UCSC Xena; proteome via PDC000173 / cBioPortal |
| **RNA-only reference** | all TCGA-BRCA primary tumours with RNA-seq and PAM50 | ~1,000 | RNA-seq | UCSC Xena |
| **Replication** | prospective CPTAC breast cohort (Krug 2020) | 122 | WES, RNA-seq, proteome, PAM50 | LinkedOmics / cBioPortal `brca_cptac_2020` |
| **Normal reference** | TCGA-BRCA adjacent normals; GTEx v10 breast and other tissues | ~110; 54 tissues | RNA-seq | Xena; GTEx portal |

The discovery and replication cohorts share no patients (Krug 2020 uses CPTAC identifiers, the
Mertins tumours carry TCGA barcodes). CPTAC breast has **no normal-tissue proteomics**, so the
tumour-vs-normal gate (§5.2) is RNA-based in every arm; this asymmetry is a limitation (§10), not a
choice to be hidden.

**Subtypes.** PAM50 **basal-like**, **HER2-enriched** and **luminal** (luminal A and B pooled).
Normal-like tumours are excluded. Three groups because the held-out cell lines support no finer split
(§6.1).

## 4. Gene universe

One universe for every arm, fixed before scoring: protein-coding genes (HGNC) that are (i) quantified
in the discovery proteome in at least 50% of the 77 tumours, (ii) quantified in RNA-seq, and
(iii) present in the DepMap 24Q4 CRISPR gene-effect matrix. Without a shared universe an arm could win
or lose on coverage alone. The universe size is recorded in `results/MANIFEST` when stage `data` runs
(expected: 7,000–9,000 genes).

## 5. Nomination arms

### 5.1 Per-layer evidence for gene *g* in subtype *s*

| Layer | Statistic |
|---|---|
| **RNA** | standardised mean difference (Hedges' *g*) of log2 expression, tumours of subtype *s* vs the cohort's other tumours |
| **Protein** | the same statistic on the CPTAC proteome log-ratios |
| **DNA** | the larger of two rank percentiles: (a) Hedges' *g* of the GISTIC2 thresholded copy-number score, *s* vs others, kept only if copy number and expression correlate in cis (Spearman ≥ 0.3 across the cohort); (b) the fraction of *s* tumours with a non-silent somatic mutation in *g*, minus the fraction in other tumours |

Each statistic is converted to a rank percentile within the universe, per subtype, so layers are
combined on a common scale.

### 5.2 Gate

A gene enters any list only if it is **up in tumour vs normal**: Hedges' *g* ≥ 0.5 for TCGA-BRCA
tumours vs TCGA adjacent normals in RNA-seq. The same gate, from the same data, applies to every arm.

### 5.3 Arms

Each arm produces a ranked list of **K = 50** genes per subtype (K = 25 and 100 as sensitivity analyses).

| Arm | Evidence | Combination | Tumours |
|---|---|---|---|
| **R** | RNA | — | the 77 discovery tumours |
| **R+P** | RNA, protein | mean of the two rank percentiles | 77 |
| **R+D** | RNA, DNA | mean of the two | 77 |
| **R+P+D** | all three | mean of the three | 77 |
| **MOFA+** | RNA, protein, copy number as three views | §5.4 | 77 |
| **R-all** | RNA | — | all ~1,000 TCGA-BRCA tumours |

R is the baseline for the information question (does the protein layer carry extra information about
the same tumours?). R-all is the baseline for the practical question (does proteomics on 77 tumours beat
RNA on 1,000, which is what a programme choosing between the two can actually buy?).

### 5.4 MOFA+ arm

MOFA+ is fitted on the 77 tumours with three views (RNA, protein, GISTIC2 scores), 10 factors, default
priors, seed 0. For each subtype the factor with the largest absolute point-biserial correlation with
membership of *s* is chosen; genes are ranked by the sum over views of their absolute standardised
weights on that factor, with the sign required to agree with "up in *s*". The literature predicts no
gain over rank stacking (`literature.md` §2); the arm is included because "integration method" is part
of the claim under test.

## 6. Held-out truth

Nothing in this section is used by any arm.

### 6.1 T1 — dependency in subtype-matched cell lines (primary)

DepMap Public **24Q4** CRISPR gene effect (Chronos), the last release under CC BY 4.0, pinned by
figshare article and sha256. Breast lines with a CRISPR screen (53 in 24Q4). Each line is assigned to
basal-like (triple-negative), HER2-enriched (ERBB2-amplified) or luminal (ER-positive, HER2-negative)
from DepMap's `ModelSubtypeFeatures` plus the receptor-status table of Dai et al. 2017; the assignment
table is committed before stage `truth` runs and its source cited per line. Expected counts are roughly
25 / 19 / 6–11, and the luminal arm is therefore the weakest.

**Primary endpoint (continuous):** the **mean Chronos gene effect of the top-K list in the *s*-matched
lines** (more negative = the list is enriched for genes those lines depend on), reported as its rank
within the matched-random floor distribution and as a difference between arms.

**Secondary endpoint (binary):** precision@K against the hit set H(*s*) = genes with mean gene effect
≤ −0.5 in *s*-lines, at least 0.25 more negative than in the other breast lines, and not in DepMap's
`CRISPRInferredCommonEssentials`. Pan-essential genes make poor targets (Chang 2021) and are excluded
from both endpoints.

### 6.2 T2 — clinical precedence (secondary)

Open Targets Platform **26.09** `clinical_target` and `clinical_indication` datasets, frozen by sha256:
genes that are the target of a drug that reached any clinical phase for a breast-cancer indication (EFO
`breast carcinoma` and descendants). Endpoint: count in the top-K list, against the matched-random
floor. This set is small and biased towards well-studied genes; the matched floor is what makes it
usable at all. Open Targets' *association score* is never used anywhere: it contains Project Score
dependencies and literature mining, which would make T1 circular and T2 popularity-driven.

### 6.3 Positive controls

Before any arm is scored against T1 or T2, the RNA arm must rank **ERBB2** in the HER2-enriched top 50
and **ESR1** in the luminal top 50. If either fails the pipeline is debugged and the failure recorded
in §12. ESR1- and ERBB2-type expression addictions are additionally flagged in every list, because they
do not reproduce in organoid screens (Neiswender 2026).

## 7. Statistics

- **Matched-random floor.** Each nominated gene is replaced by a universe gene from the same joint
  stratum of PubMed count (quintile, from NCBI gene2pubmed), mean RNA expression in the cohort
  (quintile) and protein detection fraction (tercile); strata with fewer than 20 genes are merged with
  the neighbouring expression stratum. 1,000 lists per arm and subtype. The arm's endpoint is reported
  as a percentile within that distribution.
- **Label-permutation null.** 1,000 permutations of subtype labels across the discovery tumours;
  every arm is re-nominated and re-scored.
- **Arm differences.** Paired bootstrap: 1,000 resamples of the 77 tumours with replacement, stratified
  by subtype, re-nominating every arm on each resample and scoring it; a percentile 95% interval on the
  difference in the primary endpoint between arms. Every interval is the spread of the statistic across
  resamples, never a standard error of a mean.
- **Smallest effect of interest** for P2–P4 and P6: **0.05** in mean Chronos gene effect (about one
  tenth of the dependency threshold) or **5 genes** in precision@50.
- **Replication.** Overlap of each arm's top-50 between the discovery and the Krug cohort, against the
  distribution of overlaps under the label-permutation null; T1 and T2 re-computed on the replication
  lists.
- Multiplicity: three subtypes × five arm comparisons are reported in full; no claim rests on a single
  cell of that table.

## 8. Predictions, fixed before any data are downloaded

| | Prediction | If it fails |
|---|---|---|
| **P1** | Every arm's top-50 beats its matched-random floor on the primary endpoint (percentile ≤ 5%) in at least two of the three subtypes | the truth is uninformative at this sample size; P2–P6 are reported but not interpreted |
| **P2** (primary) | R+P does not beat R by the smallest effect of interest; the interval on the difference is reported whichever way it falls | the protein layer carries list-level information about targets that RNA misses |
| **P3** | R+D adds no more over R than R+P does | DNA evidence is the more useful second layer |
| **P4** | MOFA+ does not beat R+P+D rank stacking | the integration method matters, not only the layers |
| **P5** | Every arm's discovery–replication overlap exceeds the permutation null, and the protein arm replicates no worse than RNA | nominations are cohort-specific; the layer that replicates worse is the less trustworthy one |
| **P6** | R-all (RNA on ~1,000 tumours) scores at least as well as R+P (77 tumours) | proteomics on fewer patients beats RNA on many — a result a programme would act on |
| **P7** | Fewer than 80% of the agent's checkable claims are supported by the held-out data (§9) | agent-written rationales are more faithful than the literature on hallucination suggests |
| **P8** | The agent's target scores correlate with PubMed count at Spearman ρ ≥ 0.5, and separate true candidates from matched decoys with AUC < 0.7 | the agent judges on evidence more than on fame |

## 9. The agent audit

Agentic systems now write target rationales; none has been scored claim by claim against held-out
evidence with popularity controlled (`literature.md` §7). This stage runs after everything above is
frozen.

- **Items.** The top-10 genes per subtype from the best-performing arm (30 candidates) and **30 decoys**:
  genes drawn from the same matched strata (§7) that no arm nominated. The agent is not told which is
  which.
- **Agent.** One named frontier language model, stated with its version and knowledge-cutoff date, in
  two conditions: **closed-book** (no tools) and **open-book** (read access to this repository's frozen
  evidence tables, nothing else). It is asked for a dossier per gene in a fixed schema: an overall
  target score (0–10) and a list of claims, each typed as one of *overexpressed in the subtype*,
  *up in tumour vs normal*, *dependency in subtype-matched lines*, *clinical precedence*,
  *restricted in normal tissue*, with a direction and a confidence.
- **Verification.** Every typed claim is checked automatically against the frozen tables (the cohort
  statistics, DepMap 24Q4, Open Targets 26.09, GTEx v10) and labelled **supported**, **contradicted**
  or **unsupported** (no evidence either way). Untyped free-text claims are counted but not scored.
- **Metrics.** Faithfulness = supported / (supported + contradicted), per condition and claim type;
  decoy discrimination = AUC of the agent's score for candidates vs decoys; popularity = Spearman ρ
  between the agent's score and PubMed count, and the partial correlation after adjusting for the
  primary endpoint.
- **Knowledge-cutoff caveat.** DepMap 24Q4 and the CPTAC cohorts pre-date current models' training
  data, so a closed-book model may remember some of the truth. The decoy test and the claim-level
  checks do not depend on the cutoff; as a secondary time-split check, T2 entries whose first clinical
  record post-dates the model's stated cutoff are scored separately.

## 10. Known limitations, stated in advance

- **Small cohorts.** 77 discovery tumours with protein, 122 for replication; 53 breast lines, with
  6–11 luminal. Coarse subtypes and intervals everywhere; the luminal arm may be uninformative.
- **Weak and partly confounded truth.** Dependency tracks expression (hence the matched floor); cell
  lines drift from tumours (Celligner); clinical precedence is small and popularity-biased.
- **Asymmetric layers.** No breast protein normals: the tumour-vs-normal gate is RNA in every arm.
- **Equal-weight rank stacking** is one way to integrate; MOFA+ is a second. Neither is tuned on the
  truth, which keeps the comparison honest and leaves better integration methods untested.
- **Likely null.** The design is built so that a null is a measurement, with a stated smallest effect
  of interest, and not an absence of a result.

## 11. What was known before this commit

A feasibility check on 2026-10-06, before this design was written, downloaded DepMap 24Q4
`Model.csv`, `CRISPRScreenMap.csv` and the gene-effect matrix to count breast lines with a screen (53)
and read their subtype annotations, and queried GDC, cBioPortal, PDC, ENA and figshare metadata for
cohort sizes and licences. No per-gene statistic was computed against any nomination, no tumour data
were downloaded, and none of those files are part of this repository; stage `data` re-downloads
everything with recorded sha256 sums. The literature survey behind §1 and `literature.md` was run the
same day.

## 12. Deviations from §1–§9

Each entry carries the date, what changed and why. No prediction, endpoint or threshold has changed.

1. **2026-10-06 — positive controls read on the pre-gate ranking (§6.3).** ESR1 ranks 1st in the luminal
   RNA ranking and ERBB2 4th in the HER2-enriched one, so the ranking passes both controls. The
   tumour-vs-normal gate (§5.2) removes ESR1 from every luminal list: normal breast epithelium
   expresses it (Hedges' *g* tumour vs normal = 0.19, below the 0.5 gate). The gate is doing what §5.2
   specifies, and the control tests the ranking, so the controls are evaluated before the gate. A
   consequence worth stating: a tumour-vs-normal gate excludes lineage targets, ESR1 among them.
2. **2026-10-06 — breast-carcinoma identifier (§6.2).** Open Targets 26.09 indexes "breast carcinoma" as
   `MONDO_0004989` (61 descendants); `EFO_0000305` is not a term in that release's disease index. The
   definition (the term and its descendants) is unchanged.
3. **2026-10-06 — receptor status of the cell lines (§6.1).** Lines are assigned from DepMap's own
   `ModelSubtypeFeatures` annotation where it names HER2+, ER+ or TNBC, and for the rest from DepMap
   expression: ERBB2 log2(TPM+1) ≥ 8 → HER2-enriched, ESR1 ≥ 3 → luminal, otherwise basal. The Dai et al.
   2017 table named in §6.1 was not used; the DepMap data in hand answer the same question and the basis
   is recorded per line in `results/truth/depmap_breast_lines.csv`. Result: 26 basal, 19 HER2-enriched,
   6 luminal, 2 lines with neither annotation nor expression left out.
4. **2026-10-06 — MOFA+ bootstrap (§7).** A MOFA+ fit takes about 14 s on this machine, so the MOFA+ arm
   is refitted on the first 200 of the 1,000 bootstrap resamples; the rank-stacking arms use all 1,000.
   The MOFA+ intervals are therefore wider than the others. The permutation null is unaffected: the MOFA+
   factors do not depend on the labels, only the factor choice does, which is re-done on every permutation.
5. **2026-10-06 — one source added (§3, §14).** Mertins et al. 2016 Supplementary Table 1 (Nature) is the
   source of the per-tumour QC verdict (77 pass, 28 fail) and PAM50 call for the discovery cohort; it is
   downloaded and recorded like every other input.
6. **2026-10-06 — cohort sizes as found.** TCGA-BRCA primary tumours with a PAM50 call other than
   normal-like: 821 (§3 said "~1,000"); adjacent normals with RNA-seq: 114. Krug 2020: 117 of 122 tumours
   have a non-normal-like PAM50 call. Discovery: 47 luminal, 18 basal-like, 12 HER2-enriched.
7. **2026-10-06 — matching covariates for a gene absent from a layer (§7).** In the replication cohort a
   gene missing from Krug's RNA table takes the median expression covariate when strata are formed.

## 13. Exome module

The tumour mutation calls in §3 arrive pre-processed. To show the step from raw reads to variants, the
SEQC2 reference pair **HCC1395 / HCC1395BL** (a triple-negative breast-cancer line and its matched
normal; ENA project PRJNA489865, exome pair WES_LL, ~89× on target) is aligned (BWA-MEM2) and called
(GATK Mutect2) and the calls are scored against the SEQC2 high-confidence somatic truth set v1.2.1
inside the exome target: precision and recall for SNVs (1,087 high-confidence SNVs in the target);
indel figures are indicative only (45). HCC1395 is DepMap model ACH-000699, so its called mutations are
then traced into its own CRISPR dependencies and into the basal-like candidate list. This module makes
no claim about §8; it is a methods demonstration with a published answer key.

## 14. Compute and provenance

- CPU only: 4 physical cores (8 threads), 165 GB RAM, no GPU. Every stage except the exome module runs
  in minutes on tables of a few hundred tumours and 53 cell lines; the exome pair is estimated at 6–10
  hours and is piloted on 1 million reads first.
- Every downloaded file is recorded in `results/MANIFEST.sha256` with its URL and release; every stage
  writes its parameters and git commit to `results/run_log.json`.
- Licences: TCGA (no publication restriction; the TCGA Research Network is acknowledged), CPTAC (CC BY),
  DepMap 24Q4 (CC BY 4.0), Open Targets (CC0; ChEMBL-derived fields CC BY-SA 3.0, attributed), GTEx and
  gene2pubmed (open), SEQC2 (open). No data from Sanger Project Score or Cell Model Passports are used.
