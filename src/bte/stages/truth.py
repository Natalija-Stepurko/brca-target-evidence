"""Stage `truth`: the held-out evidence (DESIGN §6). Nothing here is seen by any arm.

T1  DepMap Public 24Q4 CRISPR gene effect (Chronos) in the 53 breast lines with a screen, each assigned
    to basal / her2 / luminal; mean effect per subtype and the binary hit set H(s).
T2  Open Targets 26.09 clinical precedence: targets of drugs that reached any clinical phase for a
    breast-carcinoma indication (EFO_0000305 and descendants).
Also saved for the dossiers: tractability per target, and GTEx tissue expression.
"""
import re
import sys
import time

import pandas as pd

from bte import cohorts as K
from bte import config as C
from bte import provenance as P

RAW = C.DATA / "raw"
OUT = C.RESULTS / "truth"
# "breast carcinoma" in Open Targets 26.09 is MONDO_0004989 (the EFO id named in DESIGN §6.2 is no longer a
# term in the 26.09 disease index; see DESIGN §12)
BREAST_CARCINOMA_EFO = "MONDO_0004989"
# Open Targets 26.09 clinical stage labels as an ordinal: 1-3 trial phases, 4 = approved or post-marketing
STAGE_ORDER = {"PHASE_1": 1, "EARLY_PHASE_1": 1, "PHASE_1_2": 2, "PHASE_2": 2, "PHASE_2_3": 3, "PHASE_3": 3,
               "PHASE_4": 4, "APPROVAL": 4, "WITHDRAWAL": 4}
# receptor-status rules for lines whose DepMap annotation leaves the subtype open (log2(TPM+1))
ERBB2_HIGH, ESR1_HIGH = 8.0, 3.0


# ----------------------------------------------------------------------------- DepMap
def depmap_subtypes() -> pd.DataFrame:
    m = pd.read_csv(RAW / "depmap_24Q4_Model.csv")
    screened = set(pd.read_csv(RAW / "depmap_24Q4_CRISPRGeneEffect.csv", usecols=[0]).iloc[:, 0])
    b = m[(m.OncotreeLineage == "Breast") & m.ModelID.isin(screened)].copy()
    ex = pd.read_csv(RAW / "depmap_24Q4_OmicsExpressionProteinCodingGenesTPMLogp1.csv", index_col=0)
    col = {c.split(" (")[0]: c for c in ex.columns if c.split(" (")[0] in ("ESR1", "ERBB2")}
    e = ex.reindex(b.ModelID)[[col["ESR1"], col["ERBB2"]]]
    e.columns = ["esr1_log2tpm", "erbb2_log2tpm"]
    b = b.set_index("ModelID").join(e)

    def assign(r):
        ann = str(r.ModelSubtypeFeatures) if pd.notna(r.ModelSubtypeFeatures) else ""
        if "HER2+" in ann:
            return "her2", "annotation HER2+"
        if re.search(r"ER\+|ER,", ann):
            return "luminal", "annotation ER+"
        if "TNBC" in ann:
            return "basal", "annotation TNBC"
        if pd.isna(r.erbb2_log2tpm):
            return None, "no annotation, no expression"
        if r.erbb2_log2tpm >= ERBB2_HIGH:
            return "her2", f"ERBB2 log2(TPM+1) >= {ERBB2_HIGH}"
        if r.esr1_log2tpm >= ESR1_HIGH:
            return "luminal", f"ESR1 log2(TPM+1) >= {ESR1_HIGH}"
        return "basal", "ESR1 and ERBB2 low"

    out = b.apply(assign, axis=1, result_type="expand")
    b["subtype"], b["basis"] = out[0], out[1]
    keep = ["StrippedCellLineName", "ModelSubtypeFeatures", "esr1_log2tpm", "erbb2_log2tpm",
            "subtype", "basis"]
    return b[keep]


def gene_effect(models) -> pd.DataFrame:
    ge = pd.read_csv(RAW / "depmap_24Q4_CRISPRGeneEffect.csv", index_col=0)
    ge = ge.loc[ge.index.isin(models)]
    ge.columns = [c.split(" (")[0] for c in ge.columns]
    return ge.loc[:, ~ge.columns.duplicated()]


def common_essentials() -> set:
    ce = pd.read_csv(RAW / "depmap_24Q4_CRISPRInferredCommonEssentials.csv").iloc[:, 0]
    return set(ce.str.split(" (", regex=False).str[0])


def t1(universe_genes) -> tuple[pd.DataFrame, pd.DataFrame]:
    lines = depmap_subtypes()
    used = lines.dropna(subset=["subtype"])
    ge = gene_effect(used.index).reindex(columns=universe_genes)
    ess = common_essentials()
    t = pd.DataFrame(index=universe_genes)
    t["mean_all_breast"] = ge.mean(axis=0)
    t["n_lines_all"] = ge.notna().sum(axis=0)
    for s in C.SUBTYPES:
        ids = used.index[used.subtype == s]
        others = used.index[used.subtype != s]
        t[f"mean_{s}"] = ge.loc[ids].mean(axis=0)
        t[f"n_dependent_{s}"] = (ge.loc[ids] <= C.HIT_GENE_EFFECT).sum(axis=0)
        t[f"selectivity_{s}"] = t[f"mean_{s}"] - ge.loc[others].mean(axis=0)
    t["common_essential"] = t.index.isin(ess)
    for s in C.SUBTYPES:
        t[f"hit_{s}"] = ((t[f"mean_{s}"] <= C.HIT_GENE_EFFECT) & (t[f"selectivity_{s}"] <= -C.HIT_SELECTIVITY)
                         & ~t.common_essential)
    t["hit_any_breast"] = (t.mean_all_breast <= C.HIT_GENE_EFFECT) & ~t.common_essential
    return lines, t.rename_axis("gene")


# ----------------------------------------------------------------------------- Open Targets
def ensembl_to_symbol() -> pd.Series:
    h = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", low_memory=False,
                    usecols=["symbol", "ensembl_gene_id"])
    h = h.dropna().drop_duplicates("ensembl_gene_id")
    return h.set_index("ensembl_gene_id").symbol


def breast_indications() -> set:
    d = pd.read_parquet(RAW / "ot_26.09_disease.parquet", columns=["id", "name", "ancestors"])
    desc = d[d.ancestors.apply(lambda a: BREAST_CARCINOMA_EFO in list(a) if a is not None else False)]
    return set(desc.id) | {BREAST_CARCINOMA_EFO}


def t2(universe_genes) -> pd.DataFrame:
    sym = ensembl_to_symbol()
    ind = breast_indications()
    ct = pd.read_parquet(RAW / "ot_26.09_clinical_target.parquet")
    rows = []
    for r in ct.itertuples(index=False):
        ds = {d["diseaseId"] for d in (r.diseases if r.diseases is not None else []) if d}
        if ds & ind:
            rows.append({"targetId": r.targetId, "drugId": r.drugId, "stage": r.maxClinicalStage,
                         "max_stage": STAGE_ORDER.get(r.maxClinicalStage, 0),
                         "n_breast_indications": len(ds & ind)})
    t = pd.DataFrame(rows)
    t["gene"] = t.targetId.map(sym)
    t = t[t.max_stage >= 1].dropna(subset=["gene"])     # a clinical phase, not IND / preclinical / unknown
    agg = (t.groupby("gene").agg(max_stage=("max_stage", "max"), n_drugs=("drugId", "nunique"))
           .reset_index())
    agg["in_universe"] = agg.gene.isin(universe_genes)
    return agg.sort_values(["max_stage", "n_drugs"], ascending=False)


def tractability() -> pd.DataFrame:
    sym = ensembl_to_symbol()
    parts = [pd.read_parquet(RAW / f"ot_26.09_target_tractability_{i}.parquet") for i in (0, 1)]
    tr = pd.concat(parts)
    tr["gene"] = tr.targetId.map(sym)
    tr = tr.dropna(subset=["gene"])
    tr = tr[tr.value.astype(bool)]
    return tr.groupby(["gene", "modality"]).category.apply(lambda s: "|".join(sorted(set(s)))).unstack()


def gtex() -> pd.DataFrame:
    g = pd.read_csv(RAW / "gtex_v10_gene_median_tpm.gct.gz", sep="\t", skiprows=2).drop(columns=["Name"])
    g = g.groupby("Description").max()
    out = pd.DataFrame({"gtex_breast_tpm": g["Breast_Mammary_Tissue"],
                        "gtex_max_other_tpm": g.drop(columns=["Breast_Mammary_Tissue"]).max(axis=1),
                        "gtex_max_other_tissue": g.drop(columns=["Breast_Mammary_Tissue"]).idxmax(axis=1)})
    return out.rename_axis("gene")


def run():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    genes = list(K.universe().symbol)

    lines, t = t1(genes)
    lines.to_csv(OUT / "depmap_breast_lines.csv")
    t.to_csv(OUT / "t1_gene_effect.csv")
    counts = lines.subtype.value_counts(dropna=False).to_dict()
    print(f"  DepMap breast lines with a screen: {len(lines)}; assigned {counts}")
    for s in C.SUBTYPES:
        print(f"  H({s}): {int(t[f'hit_{s}'].sum()):>4} hits; mean effect over the universe "
              f"{t[f'mean_{s}'].mean():.3f}")
    print(f"  common-essential genes in universe: {int(t.common_essential.sum()):,}")

    prec = t2(genes)
    prec.to_csv(OUT / "t2_clinical_precedence.csv", index=False)
    print(f"  T2: {len(prec)} genes with a clinical-stage drug for a breast indication; "
          f"{int(prec.in_universe.sum())} in the universe; approved (stage 4): "
          f"{int((prec.max_stage >= 4).sum())}")

    tractability().to_csv(OUT / "tractability.csv")
    gtex().to_csv(OUT / "gtex.csv")

    summary = {"n_lines": int(len(lines)), "lines_per_subtype": {k if k == k else "unassigned": int(v)
               for k, v in counts.items()}, "hits": {s: int(t[f"hit_{s}"].sum()) for s in C.SUBTYPES},
               "hit_any_breast": int(t.hit_any_breast.sum()),
               "common_essential_in_universe": int(t.common_essential.sum()),
               "t2_genes": int(len(prec)), "t2_in_universe": int(prec.in_universe.sum()),
               "erbb2_high": ERBB2_HIGH, "esr1_high": ESR1_HIGH}
    pd.Series(summary).to_json(OUT / "summary.json", indent=1)
    P.log_run("truth", {"hit_gene_effect": C.HIT_GENE_EFFECT, "hit_selectivity": C.HIT_SELECTIVITY,
                        "depmap": C.DEPMAP_RELEASE, "open_targets": C.OPEN_TARGETS_RELEASE},
              [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*"))], time.time() - t0)
    print(f"done in {time.time() - t0:,.0f}s")


if __name__ == "__main__":
    sys.exit(run())
