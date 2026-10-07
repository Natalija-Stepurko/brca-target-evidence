"""Every input file: where it comes from, which release, and the local name it is saved under.

Nothing here is downloaded until `bte run data` is called, and that stage records the sha256 of every
file in results/MANIFEST.sha256. Sizes are what the servers reported on 2026-10-06 (HEAD requests) and
are used only to show progress.
"""
from dataclasses import dataclass

XENA = "https://tcga.xenahubs.net/download"
LINKED_TCGA = "https://linkedomics.org/data_download/TCGA-BRCA"
LINKED_CPTAC = "https://www.linkedomics.org/data_download/CPTAC-BRCA"
FIGSHARE = "https://ndownloader.figshare.com/files"
OT = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output"
NCBI = "https://ftp.ncbi.nlm.nih.gov/gene/DATA"


@dataclass(frozen=True)
class Source:
    key: str            # short name used by the stages
    url: str
    local: str          # file name under data/raw/
    release: str
    licence: str
    size_mb: float      # as reported by the server, for progress only
    note: str = ""


SOURCES = [
    # --- TCGA-BRCA (UCSC Xena, TCGA hub) -------------------------------------------------------
    Source("tcga_rna", f"{XENA}/TCGA.BRCA.sampleMap/HiSeqV2.gz", "tcga_brca_HiSeqV2.gz",
           "Xena TCGA hub, 2017-10 build", "TCGA: no publication restriction", 64.3,
           "log2(RSEM+1), 1,218 samples including ~114 adjacent normals (sample code -11)"),
    Source("tcga_cnv", f"{XENA}/TCGA.BRCA.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz",
           "tcga_brca_gistic2_thresholded.gz", "Xena TCGA hub", "TCGA", 2.4,
           "GISTIC2 thresholded calls: -2, -1, 0, 1, 2 per gene and sample"),
    Source("tcga_mut", f"{XENA}/mc3/BRCA_mc3.txt.gz", "tcga_brca_mc3.txt.gz",
           "MC3 public MAF, BRCA subset", "TCGA", 2.1, "one row per somatic variant"),
    Source("tcga_clin", f"{XENA}/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix", "tcga_brca_clinicalMatrix.tsv",
           "Xena TCGA hub", "TCGA", 1.7, "PAM50Call_RNAseq column gives the subtype"),
    Source("tcga_surv", "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/"
           "Survival_SupplementalTable_S1_20171025_xena_sp", "tcga_cdr_survival.tsv",
           "TCGA-CDR (Liu 2018)", "TCGA", 1.9),
    # --- CPTAC proteome on the TCGA tumours (Mertins 2016), via LinkedOmics -----------------------
    Source("cptac16_prot", f"{LINKED_TCGA}/Human__TCGA_BRCA__BI__Proteome__QExact__01_28_2016__BI__Gene__"
           "CDAP_iTRAQ_UnsharedLogRatio_r2.cct.gz", "cptac2016_proteome_gene_itraq_logratio.cct.gz",
           "CPTAC CDAP r2, 2016-01-28", "CPTAC: CC BY", 2.5,
           "gene-level iTRAQ log-ratios, TCGA sample barcodes as columns"),
    Source("mertins_tables", "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature18003/"
           "MediaObjects/41586_2016_BFnature18003_MOESM111_ESM.zip", "mertins2016_supplementary_tables.zip",
           "Nature 534:55 (2016), Supplementary Tables 1-19", "Springer Nature supplementary material", 139.0,
           "Table 1 carries the per-sample QC status (77 pass, 28 fail) and PAM50; extracted by the data stage"),
    # --- CPTAC prospective breast cohort (Krug 2020), via LinkedOmics ----------------------------
    Source("krug_prot", f"{LINKED_CPTAC}/HS_CPTAC_BRCA_2018_Proteome_Ratio_Norm_gene_Median.cct",
           "krug2020_proteome_gene_median.cct", "LinkedOmics CPTAC-BRCA 2018", "CPTAC: CC BY", 12.0),
    Source("krug_rna", f"{LINKED_CPTAC}/HS_CPTAC_BRCA_2018_RNA_GENE.cct", "krug2020_rna_gene.cct",
           "LinkedOmics CPTAC-BRCA 2018", "CPTAC: CC BY", 20.0),
    Source("krug_cnv", f"{LINKED_CPTAC}/HS_CPTAC_BRCA_2018_CNA.cct", "krug2020_cna_gene.cct",
           "LinkedOmics CPTAC-BRCA 2018", "CPTAC: CC BY", 15.0),
    Source("krug_mut", f"{LINKED_CPTAC}/HS_CPTAC_BRCA_2018_MUT_GENE.cbt", "krug2020_mutation_gene.cbt",
           "LinkedOmics CPTAC-BRCA 2018", "CPTAC: CC BY", 3.0, "binary gene x sample mutation matrix"),
    Source("krug_clin", f"{LINKED_CPTAC}/HS_CPTAC_BRCA_2018_CLI.tsi", "krug2020_clinical.tsi",
           "LinkedOmics CPTAC-BRCA 2018", "CPTAC: CC BY", 0.1, "includes PAM50"),
    # --- DepMap Public 24Q4 (figshare 27993248): the last CC BY 4.0 release --------------------
    Source("depmap_effect", f"{FIGSHARE}/51064667", "depmap_24Q4_CRISPRGeneEffect.csv",
           "DepMap Public 24Q4", "CC BY 4.0", 428.0, "Chronos gene effect, models x genes"),
    Source("depmap_model", f"{FIGSHARE}/51065297", "depmap_24Q4_Model.csv", "DepMap Public 24Q4", "CC BY 4.0", 0.6),
    Source("depmap_common_essential", f"{FIGSHARE}/51064916", "depmap_24Q4_CRISPRInferredCommonEssentials.csv",
           "DepMap Public 24Q4", "CC BY 4.0", 0.1),
    Source("depmap_screenmap", f"{FIGSHARE}/51065159", "depmap_24Q4_CRISPRScreenMap.csv",
           "DepMap Public 24Q4", "CC BY 4.0", 0.1, "which library screened which model"),
    Source("depmap_expr", f"{FIGSHARE}/51065489", "depmap_24Q4_OmicsExpressionProteinCodingGenesTPMLogp1.csv",
           "DepMap Public 24Q4", "CC BY 4.0", 506.0, "for ESR1 / ERBB2 receptor status of the lines"),
    Source("depmap_readme", f"{FIGSHARE}/51065795", "depmap_24Q4_README.txt", "DepMap Public 24Q4", "CC BY 4.0", 0.1),
    # --- Open Targets Platform 26.09 ---------------------------------------------------------
    Source("ot_clinical_target", f"{OT}/clinical_target/00000000.parquet", "ot_26.09_clinical_target.parquet",
           "Open Targets 26.09", "CC0; ChEMBL-derived fields CC BY-SA 3.0", 2.7),
    Source("ot_clinical_indication", f"{OT}/clinical_indication/00000000.parquet",
           "ot_26.09_clinical_indication.parquet", "Open Targets 26.09", "CC0", 5.0),
    Source("ot_disease", f"{OT}/disease/00000000.parquet", "ot_26.09_disease.parquet", "Open Targets 26.09",
           "CC0", 30.0, "EFO ontology with descendants, to define breast-carcinoma indications"),
    Source("ot_tractability_0", f"{OT}/target_tractability/part-00000-e2201205-916c-4b28-971d-6511ae9ea691-c000.zstd.parquet",
           "ot_26.09_target_tractability_0.parquet", "Open Targets 26.09", "CC0", 0.1),
    Source("ot_tractability_1", f"{OT}/target_tractability/part-00001-e2201205-916c-4b28-971d-6511ae9ea691-c000.zstd.parquet",
           "ot_26.09_target_tractability_1.parquet", "Open Targets 26.09", "CC0", 0.1),
    # --- normal tissue, gene identifiers, popularity ---------------------------------------------
    Source("gtex_median", "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/"
           "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz", "gtex_v10_gene_median_tpm.gct.gz",
           "GTEx v10", "open", 8.8),
    Source("hgnc", "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
           "hgnc_complete_set.txt", "HGNC, current", "CC0", 17.0, "protein-coding set and symbol history"),
    Source("gene2pubmed", f"{NCBI}/gene2pubmed.gz", "ncbi_gene2pubmed.gz", "NCBI Gene, current", "open", 288.0,
           "all species; filtered to tax_id 9606 on load"),
    Source("gene_info", f"{NCBI}/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz", "ncbi_Homo_sapiens.gene_info.gz",
           "NCBI Gene, current", "open", 5.2, "Entrez id to symbol"),
]

BY_KEY = {s.key: s for s in SOURCES}
