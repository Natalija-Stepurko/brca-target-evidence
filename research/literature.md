# Related work and data access

Survey run 2026-10-06 as three parallel searches (what has been done; how targets are validated and how
prioritisation is judged; data access and licences), then cross-checked. The design that follows from
it is [`../docs/DESIGN.md`](../docs/DESIGN.md); the section numbers cited there (§1, §2, §7) are this
file's.

Strength tags: **[F]** read in the full text · **[Ab]** abstract or record only · **[S]** secondary
source (press, blog, snippet) · **[U]** unverified. Entries dated 2026 post-date the model's training
data and rest entirely on the search; re-check before citing anything in public.

---

## 1. Multi-omics target nomination already published

- **Mertins et al. 2016, Nature** — 105 TCGA breast tumours, 77 passed QC; proteome and phosphoproteome.
  Nominated amplicon-linked hyper-phosphorylated kinases (CDK12, PAK1, PTK2, RIPK2, TLK2) and a GPCR
  phospho-cluster not seen at mRNA level. Validation: connectivity to LINCS shRNA signatures with a
  permutation FDR. No CRISPR test, no RNA-only arm. [F]
- **Krug et al. 2020, Cell** — 122 prospective tumours with phospho- and acetyl-proteomes. Nominations
  include Rb protein as a CDK4/6-inhibitor marker (checked against palbociclib response in GDSC breast lines,
  rho −0.61) and immune-checkpoint candidates. States that its candidate drivers **do not overlap**
  with Mertins 2016's (attributed to platform and CMap changes). [F]
- **Savage et al. 2024, Cell** — 1,043 patients, 10 cancer types. A target = protein (or activating
  phosphosite) up in tumour vs normal **and** a lineage-matched DepMap dependency. Against PRISM drug
  response (base rate 15%): CRISPR alone 22%, protein alone 29%, combined 39%. No mRNA arm; unmatched
  null; **breast excluded** from the tumour-vs-normal and PRISM analyses (no breast normals). [F]
- **Lin, Dai & Pusztai 2025, Breast Cancer Res Treat** — DepMap dependency in 48 breast lines, minus
  pan-essentials and unexpressed genes, ranked by DGIdb druggability: 66 ER+, 53 HER2+, 29 TNBC
  targets. No proteomics, no null. The closest match to the validation arm. [Ab]
- **Pacini et al. 2024, Cancer Cell** — 930 lines, 370 priority targets across 27 cancer types; most
  addictions are gain-of-function, not synthetic-lethal. [Ab]
- **Behan et al. 2019, Nature** — CRISPR priority-scoring framework (324 lines; WRN in MSI). [Ab]
- Also: Li et al. 2023 Cell (drivers → functional states) [Ab]; Anurag et al. 2022 Cancer Discov (TNBC
  chemo response) [Ab]; Neill et al. 2025 Cancer Discov (PTPN12 loss in TNBC) [Ab]; Ding et al. 2025
  Cancer Biol Med (dependency markers, 47 breast lines) [Ab]; Deng et al. 2024 Cell Rep Methods
  (Multiomics2Targets on CPTAC3) [Ab].

**Not found anywhere:** integrated vs single-omics nomination at equal list size; popularity- or
expression-matched random genes; held-out ground truth; cross-cohort reproducibility of nominations.

## 2. Does integration beat single-omics?

- **Rappoport & Shamir 2018, NAR** — multi-omics "does not consistently provide better prognostic value
  and clinical significance" than single-omics; expression alone gave the best clinical enrichment. [F]
- **Herrmann et al. 2021, Brief Bioinform** — 18 TCGA datasets, 11 methods; only block forest beat a
  clinical-only Cox model, "and only slightly". [F]
- **Cantini et al. 2021, Nat Commun** — 9 joint dimensionality-reduction methods (incl. MOFA v1) on 10
  TCGA cancers; survival-linked factors depend more on cancer type than method. No single-omics arm. [F]
- **Wissel et al. 2023, Cell Rep Methods** — adding modalities reduces noise robustness; use only
  modalities with known value. [Ab]
- **Bernett et al. 2026, Nat Commun (DrEval)** — expression and drug features carry the signal; CNV,
  methylation, mutation add little; proteomics ≈ expression; deep models barely beat a mean-effect
  baseline. [F]
- **Dempster et al. 2020, bioRxiv** — expression beats DNA features for predicting vulnerabilities. [Ab]

## 3. What protein adds beyond mRNA

- Tumour mRNA–protein median gene-wise correlation: 0.39 (Mertins 2016), 0.41 (Krug 2020), 0.34–0.61
  by cohort (Savage 2024, overall 0.48). Druggable genes correlate better; 19 druggable genes (mostly
  secreted, plus HDAC3, CDK9) show no positive correlation in any cohort. [F]
- **Nusinow et al. 2020, Cell** — 375 lines, mean correlation 0.48; protein-complex co-variation is
  largely absent from RNA. [F]
- **Gonçalves et al. 2022, Cancer Cell** — 949 lines; protein predicts drug response and CRISPR
  essentiality "very similar" to the transcriptome; 1,500 random proteins keep 88% of drug-prediction
  power; thousands of protein–vulnerability links not significant at transcript level. [F]
- **Chen et al. 2020, JCO CCI** — among RPPA genes, protein is the best same-gene predictor of CRISPR
  dependency for 38.5% of genes vs 20.5% for mRNA. [F]

## 4. DepMap as validation evidence

- DepMap 24Q4: 1,178 CRISPR-screened models, **53 breast** (of 96 breast models); 26Q1 still 53.
  Subtype annotation thin and free-text: HER2+ ~19, TNBC-labelled ~25, **ER+/HER2− only ~6**. 20 lines
  screened in both Avana and Sanger KY libraries (a screen-replication rung). Computed from the
  release files. [F]
- **Dempster et al. 2021** — Chronos, the current gene-effect model. [Ab]
- **Warren et al. 2021, Nat Commun (Celligner)** — breast lines align comparatively well with tumours;
  a mesenchymal/undifferentiated cluster of lines (incl. some basal-B breast) has distinct
  dependencies. [F]
- **Neiswender et al. 2026, Nature** — in 28 breast organoids, ESR1 expression-addiction is absent
  (r = 0.033 vs −0.637 in traditional lines) and ERBB2 addiction nearly absent. [F]
- **Lin et al. 2017 eLife; Lin et al. 2019 Sci Transl Med** — RNAi-era targets (MELK in TNBC; 10
  clinical-stage targets) not essential by CRISPR. [Ab]
- **Chang et al. 2021, Cancer Cell** — pan-essential genes make poor targets. [Ab]
- **Chiu et al. 2021 (DeepDEP)** — predicted tumour dependency maps; predictions, not evidence. [Ab]

## 5. What predicts clinical success

- **Minikel et al. 2024, Nature** — genetic support: relative success 2.6 overall; somatic (IntOGen)
  evidence RS 2.3 in oncology; 44% of genetically supported oncology-kinase pairs already reached
  phase I. [F]
- **Nelson et al. 2015, Nat Genet; King et al. 2019** — genetic support roughly doubles success. [Ab]
- **Wong et al. 2019** — oncology probability of success 3.4%; biomarker-stratified trials 10.3% vs
  5.5%. [F]
- **Czech et al. 2024 (preprint)** — time-aware evidence model; top 2% of pairs 4–5× likelier to pass
  phase 2. [Ab]
- **Gap:** no peer-reviewed estimate of whether CRISPR-dependency evidence predicts clinical
  advancement. [F, by absence]

## 6. Benchmarks, ground truth and popularity bias

- **Open Targets (Buniello et al. 2025, NAR)** — association score includes Project Score dependencies
  and Europe PMC text mining; using it as truth for a DepMap-informed method is circular. From 26.03
  drug data come from new `clinical_*` datasets. [F/Ab]
- **Leung et al. 2026, Sci Rep (TargetBench 1.0)** — clinical-stage positives, positives-unlabelled
  negatives, nested CV with no time split, Open Targets and 7 LLMs as baselines; 38 diseases, no breast;
  warns of temporal leakage in text features. [F]
- **Stoeger et al. 2018, PLoS Biol; Haynes et al. 2018, Sci Rep; Richardson et al. 2024, eLife** —
  research concentrates on ~2,000 genes; annotation inequality steers hit selection. [Ab]
- **Khan et al. 2025, Bioinformatics** — naive LLM gene scores vs PubMed frequency ρ = 0.795. [F]
- **Neeley et al. 2025, Bioinform Adv** — LLM bias to well-studied genes, input-order sensitivity. [Ab]
- **Bonner et al. 2022, Brief Bioinform; Guney et al. 2016, Nat Commun** — node-degree bias in graph
  scorers; degree-preserving randomisation as control. [Ab/S]
- **Not found:** random controls matched jointly on publication count, expression and essentiality. [F, by absence]

## 7. Agentic and LLM systems for target work

- **Biomni (Huang et al. 2026, Science)** — 433 instances, 10 tasks; gene retrieval against published
  answers, not target nomination against outcomes. [F/Ab]
- **AI co-scientist (Gottweis et al. 2026, Nature)** — AML repurposing (5 tested, 3 active) and liver
  fibrosis (2 of 3 active in organoids); novelty judge unreliable without search. [F]
- **Robin (Ghareeb et al. 2026, Nature)** — one disease (dry AMD); 44.5% hallucinated references when
  literature agents are swapped for a general LLM. [F]
- **Virtual Lab (Swanson et al. 2025, Nature)** — nanobody design, not target ID. [Ab]
- **TxAgent (Gao et al. 2025, arXiv)** — therapeutic reasoning, no target-ID evaluation. [F]
- **Gap:** no study scores agent-written target rationales claim by claim against time-split held-out
  evidence with popularity control. [F, by absence]

## 8. Data access facts (checked 2026-10-06)

- TCGA-BRCA open tier: 778 patients with mutation + CNV + RNA + RPPA from a primary tumour (GDC DR46);
  852 samples in the PanCanAtlas cBioPortal study.
- Mertins 2016 tumours carry TCGA barcodes, all with PanCanAtlas DNA/RNA/RPPA. Krug 2020 uses CPTAC IDs,
  disjoint from TCGA.
- CPTAC data CC BY (PDC). The `cptac` package (1.5.14, June 2024) still wires BRCA but its Zenodo
  downloads fail (issue #80).
- DepMap: 24Q4 is the last CC BY 4.0 release on figshare; later terms bar commercial and AI-product use
  and require reposting the terms with any rehosted data; the portal is captcha-gated.
- Sanger Project Score / Cell Model Passports: non-commercial, no rehosting.
- Open Targets 26.09 parquet on EBI FTP, CC0 (ChEMBL parts CC BY-SA 3.0).
- SEQC2 HCC1395 exomes: PRJNA489865 on ENA, public; truth set v1.2.1 (GRCh38): 1,087 high-confidence
  SNVs and 45 indels inside the exome target. HCC1395 = ACH-000699 (CRISPR-screened); HCC1395BL has
  no screen.

## References

1. Mertins P et al. 2016. Nature 534:55–62. doi:10.1038/nature18003
2. Krug K et al. 2020. Cell 183:1436–1456.e31. doi:10.1016/j.cell.2020.10.036
3. Savage SR et al. 2024. Cell 187:4389–4407.e15. doi:10.1016/j.cell.2024.05.039
4. Lin HK, Dai J, Pusztai L. 2025. Breast Cancer Res Treat 214:319–327. doi:10.1007/s10549-025-07817-0
5. Pacini C et al. 2024. Cancer Cell 42:301–316.e9. doi:10.1016/j.ccell.2023.12.016
6. Behan FM et al. 2019. Nature 568:511–516. doi:10.1038/s41586-019-1103-9
7. Li Y et al. 2023. Cell 186:3921–3944.e25. doi:10.1016/j.cell.2023.07.014
8. Anurag M et al. 2022. Cancer Discov 12:2586–2605. doi:10.1158/2159-8290.CD-22-0200
9. Neill NJ et al. 2025. Cancer Discov 15:2326–2343. doi:10.1158/2159-8290.CD-23-1173
10. Ding R, Shao Z, Yu T. 2025. Cancer Biol Med 22. doi:10.20892/j.issn.2095-3941.2025.0290
11. Deng EZ et al. 2024. Cell Rep Methods 4:100839. doi:10.1016/j.crmeth.2024.100839
12. Rappoport N, Shamir R. 2018. Nucleic Acids Res 46:10546–10562. doi:10.1093/nar/gky889
13. Herrmann M et al. 2021. Brief Bioinform 22:bbaa167. doi:10.1093/bib/bbaa167
14. Cantini L et al. 2021. Nat Commun 12:124. doi:10.1038/s41467-020-20430-7
15. Wissel D, Rowson D, Boeva V. 2023. Cell Rep Methods 3:100461. doi:10.1016/j.crmeth.2023.100461
16. Bernett J et al. 2026. Nat Commun 17:4238. doi:10.1038/s41467-026-72903-w
17. Dempster JM et al. 2020. bioRxiv. doi:10.1101/2020.02.21.959627
18. Nusinow DP et al. 2020. Cell 180:387–402.e16. doi:10.1016/j.cell.2019.12.023
19. Gonçalves E et al. 2022. Cancer Cell 40:835–849.e8. doi:10.1016/j.ccell.2022.06.010
20. Chen MM et al. 2020. JCO Clin Cancer Inform 4:357–366. doi:10.1200/CCI.19.00144
21. Haynes WA, Tomczak A, Khatri P. 2018. Sci Rep 8:1362. doi:10.1038/s41598-018-19333-x
22. Minikel EV et al. 2024. Nature 629:624–629. doi:10.1038/s41586-024-07316-0
23. Neiswender et al. 2026. Nature 657:789–799 (DOI not recorded by the search)
24. Warren A et al. 2021. Nat Commun (Celligner) (DOI not recorded)
25. Leung et al. 2026. Sci Rep 16:17018 (TargetBench 1.0)
26. Khan et al. 2025. Bioinformatics 41:btaf541
27. Stoeger T et al. 2018. PLoS Biol 16:e2006643
28. Richardson et al. 2024. eLife 12:RP93429
29. Buniello A et al. 2025. Nucleic Acids Res 53:D1467
30. Huang K et al. 2026. Science 393(6813). doi:10.1126/science.adz4351
31. Gottweis J et al. 2026. Nature 655:487–496
32. Ghareeb AE et al. 2026. Nature 655:497–505
33. Swanson K et al. 2025. Nature 646:716–723
34. Fang LT et al. 2021. Nat Biotechnol. doi:10.1038/s41587-021-00993-6
35. Xiao W et al. 2021. Nat Biotechnol. doi:10.1038/s41587-021-00994-5
36. Liu J et al. 2018. Cell (TCGA-CDR). doi:10.1016/j.cell.2018.02.052
