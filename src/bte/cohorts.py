"""Load each cohort into one shape: a dict of gene x sample matrices plus a per-sample subtype label.

Discovery  = the 77 QC-pass TCGA-BRCA tumours with CPTAC proteomics (Mertins 2016).
TCGA-all   = every TCGA-BRCA primary tumour with a PAM50 call (the R-all arm) plus adjacent normals.
Krug       = the 122-tumour prospective CPTAC cohort (replication).
Subtypes are collapsed to DESIGN §3: basal (PAM50 Basal), her2 (Her2), luminal (LumA + LumB);
normal-like tumours are dropped.
"""
import gzip
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from bte import config as C

RAW = C.DATA / "raw"
PAM50_TO_SUBTYPE = {"Basal": "basal", "Her2": "her2", "LumA": "luminal", "LumB": "luminal"}
NON_SILENT = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
              "In_Frame_Ins", "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"}


@dataclass
class Cohort:
    name: str
    subtype: pd.Series                      # sample -> basal / her2 / luminal
    rna: pd.DataFrame                       # gene x sample, log2 scale
    protein: pd.DataFrame | None = None     # gene x sample, log ratio
    cnv: pd.DataFrame | None = None         # gene x sample, GISTIC thresholds or log2 ratio
    mutation: pd.DataFrame | None = None    # gene x sample, 0/1 non-silent
    normals: pd.DataFrame | None = None     # gene x sample RNA of adjacent normals (TCGA only)
    notes: dict = field(default_factory=dict)

    def layers(self):
        return {k: v for k, v in (("rna", self.rna), ("protein", self.protein), ("cnv", self.cnv),
                                  ("mutation", self.mutation)) if v is not None}


def _read_tsv(path, **kw):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        return pd.read_csv(f, sep="\t", index_col=0, **kw)


def universe() -> pd.DataFrame:
    return pd.read_csv(C.RESULTS / "data" / "universe.csv")


def restrict(m: pd.DataFrame, genes) -> pd.DataFrame:
    """Keep universe genes, in universe order; genes absent from a layer become NaN rows."""
    m = m[~m.index.duplicated(keep="first")]
    return m.reindex(genes)


# ----------------------------------------------------------------------------- TCGA
def tcga_rna() -> pd.DataFrame:
    return _read_tsv(RAW / "tcga_brca_HiSeqV2.gz")


def tcga_subtypes() -> pd.Series:
    c = pd.read_csv(RAW / "tcga_brca_clinicalMatrix.tsv", sep="\t", index_col=0,
                    usecols=["sampleID", "PAM50Call_RNAseq", "sample_type"])
    c = c[c.sample_type == "Primary Tumor"]
    return c.PAM50Call_RNAseq.map(PAM50_TO_SUBTYPE).dropna()


def tcga_mutations(samples) -> pd.DataFrame:
    m = pd.read_csv(RAW / "tcga_brca_mc3.txt.gz", sep="\t", usecols=["sample", "gene", "effect"])
    m = m[m.effect.isin(NON_SILENT) & m["sample"].isin(samples)]
    return pd.crosstab(m.gene, m["sample"]).clip(upper=1).reindex(columns=list(samples), fill_value=0)


def discovery_table() -> pd.DataFrame:
    return pd.read_csv(C.RESULTS / "data" / "discovery_samples.csv")


def load_discovery(genes) -> Cohort:
    d = discovery_table()
    d = d[d.qc == "pass"]
    sample_ids = [f"{t}-01" for t in d.tcga_id]              # Xena uses the -01 primary-tumour suffix
    sub = pd.Series(d.pam50.map(PAM50_TO_SUBTYPE).values, index=sample_ids).dropna()
    rna_all = tcga_rna()
    missing = [s for s in sample_ids if s not in rna_all.columns]
    assert not missing, f"discovery tumours absent from TCGA RNA: {missing}"
    prot = _read_tsv(RAW / "cptac2016_proteome_gene_itraq_logratio.cct.gz")
    prot.columns = prot.columns.str.replace(".", "-", regex=False) + "-01"
    cnv = _read_tsv(RAW / "tcga_brca_gistic2_thresholded.gz")
    normals = rna_all[[c for c in rna_all.columns if c.endswith("-11")]]
    return Cohort("discovery", sub[sub.index], restrict(rna_all[sub.index], genes),
                  restrict(prot[sub.index], genes), restrict(cnv.reindex(columns=sub.index), genes),
                  restrict(tcga_mutations(sub.index), genes).fillna(0), restrict(normals, genes),
                  notes={"n_qc_pass": len(d), "n_with_subtype": len(sub)})


def load_tcga_all(genes) -> Cohort:
    rna_all = tcga_rna()
    sub = tcga_subtypes()
    sub = sub[sub.index.isin(rna_all.columns)]
    normals = rna_all[[c for c in rna_all.columns if c.endswith("-11")]]
    return Cohort("tcga_all", sub, restrict(rna_all[sub.index], genes), normals=restrict(normals, genes),
                  notes={"n_with_subtype": len(sub), "n_normals": normals.shape[1]})


# ----------------------------------------------------------------------------- Krug 2020
def load_krug(genes) -> Cohort:
    cli = _read_tsv(RAW / "krug2020_clinical.tsi").drop(index="IDX", errors="ignore")
    sub = cli.PAM50.map(PAM50_TO_SUBTYPE).dropna()
    rna = _read_tsv(RAW / "krug2020_rna_gene.cct")
    prot = _read_tsv(RAW / "krug2020_proteome_gene_median.cct")
    cna = _read_tsv(RAW / "krug2020_cna_gene.cct")
    mut = _read_tsv(RAW / "krug2020_mutation_gene.cbt")
    common = [s for s in sub.index if s in rna.columns and s in prot.columns and s in cna.columns]
    sub = sub[common]
    mut = mut.reindex(columns=common).fillna(0)
    return Cohort("krug", sub, restrict(rna[common].astype(float), genes),
                  restrict(prot[common].astype(float), genes), restrict(cna[common].astype(float), genes),
                  restrict(mut, genes).fillna(0),
                  notes={"n_with_subtype": len(sub), "pam50": cli.PAM50.value_counts().to_dict()})


# ----------------------------------------------------------------------------- statistics
def hedges_g(x: pd.DataFrame, mask: np.ndarray) -> pd.Series:
    """Standardised mean difference, columns in `mask` vs the rest, per row; NaN-aware."""
    a, b = x.loc[:, mask], x.loc[:, ~mask]
    na, nb = a.notna().sum(axis=1), b.notna().sum(axis=1)
    ma, mb = a.mean(axis=1), b.mean(axis=1)
    va, vb = a.var(axis=1, ddof=1), b.var(axis=1, ddof=1)
    sp = np.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    g = (ma - mb) / sp.replace(0, np.nan)
    j = 1 - 3 / (4 * (na + nb) - 9)                                     # small-sample correction
    return (g * j).where((na >= 3) & (nb >= 3))
