"""Fixed study parameters. Every value here is stated in docs/DESIGN.md; tests check that they agree."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RESULTS = ROOT / "results"

# §3 cohorts
SUBTYPES = ("basal", "her2", "luminal")           # PAM50 basal-like, HER2-enriched, luminal A+B pooled

# §5 nomination
K = 50                                            # list size per subtype
K_SENSITIVITY = (25, 100)
GATE_HEDGES_G = 0.5                               # tumour vs adjacent normal, RNA, every arm
CIS_SPEARMAN_MIN = 0.3                            # copy number must track expression to count
MOFA_FACTORS = 10
MOFA_SEED = 0
ARMS = ("R", "R+P", "R+D", "R+P+D", "MOFA+", "R-all")

# §6 held-out truth
DEPMAP_RELEASE = "24Q4"                           # last CC BY 4.0 release; pinned
DEPMAP_FIGSHARE_ARTICLE = 27993248
HIT_GENE_EFFECT = -0.5                            # mean Chronos effect in subtype-matched lines
HIT_SELECTIVITY = 0.25                            # more negative than the other breast lines by at least this
OPEN_TARGETS_RELEASE = "26.09"

# §7 statistics
N_MATCHED_RANDOM = 1000
N_PERMUTATIONS = 1000
N_BOOTSTRAP = 1000
SEOI_GENE_EFFECT = 0.05                           # smallest effect of interest, mean Chronos units
SEOI_PRECISION = 5                                # genes in 50
MATCH_STRATA = {"pubmed": 5, "expression": 5, "protein_detection": 3}
MIN_STRATUM = 20

# §9 agent audit
AUDIT_CANDIDATES_PER_SUBTYPE = 10
AUDIT_DECOYS = 30
CLAIM_TYPES = ("overexpressed_in_subtype", "tumour_vs_normal", "dependency_in_subtype_lines",
               "clinical_precedence", "normal_tissue_restricted")

# §13 exome module
SEQC2_PROJECT = "PRJNA489865"
SEQC2_WES_PAIR = {"tumour": "SRR7890850", "normal": "SRR7890851"}   # WES_LL, ~89x on target
SEQC2_TRUTH_VERSION = "v1.2.1"
HCC1395_DEPMAP_ID = "ACH-000699"
