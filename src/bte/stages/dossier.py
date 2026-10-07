"""Stage `dossier`: every rung for the nominated candidates (DESIGN §3.5 of the plan; README step 6).

For each subtype, the best arm is the one whose top-50 has the most negative mean gene effect in the
matched lines (results/ladder/observed.csv). Its 50 genes get: dependency (mean effect, lines dependent,
selectivity, hit flags), clinical precedence and tractability, GTEx breast vs other tissues, malignant
vs non-malignant expression in the 50,002-cell atlas of study 1, survival association in TCGA-BRCA
(Cox with age and stage), and replication status in the Krug list of the same arm.
"""
import json
import sys
import time

import numpy as np
import pandas as pd

from bte import cohorts as K
from bte import config as C
from bte import provenance as P

OUT = C.RESULTS / "dossier"
RAW = C.DATA / "raw"
ATLAS_PATH = "/data/scfm/data/atlas.h5ad"
TOP = 10


# ----------------------------------------------------------------------------- atlas
def atlas_specificity(genes) -> pd.DataFrame:
    """Mean log1p expression and fraction of cells expressing, malignant vs every other cell."""
    import anndata as ad
    a = ad.read_h5ad(ATLAS_PATH)
    keep = [g for g in genes if g in a.var_names]
    x = a[:, keep].X.tocsc()
    mal = (a.obs.cell_type == "malignant cell").values
    out = pd.DataFrame(index=keep)
    for lab, m in (("malignant", mal), ("non_malignant", ~mal)):
        sub = x[m]
        out[f"atlas_mean_{lab}"] = np.asarray(sub.mean(axis=0)).ravel()
        out[f"atlas_frac_{lab}"] = np.asarray((sub > 0).mean(axis=0)).ravel()
    out["atlas_malignant_ratio"] = (out.atlas_mean_malignant + 1e-3) / (out.atlas_mean_non_malignant + 1e-3)
    out["atlas_n_malignant"], out["atlas_n_other"] = int(mal.sum()), int((~mal).sum())
    return out.rename_axis("gene")


# ----------------------------------------------------------------------------- survival
def tcga_survival(genes) -> pd.DataFrame:
    """Cox proportional hazards per gene (z-scored log2 expression) with age and stage, overall survival,
    TCGA-BRCA primary tumours; hazard ratio per SD and its p-value; univariate concordance."""
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    rna = K.tcga_rna()
    cdr = pd.read_csv(RAW / "tcga_cdr_survival.tsv", sep="\t", index_col=0)
    cdr = cdr[cdr["cancer type abbreviation"] == "BRCA"]
    stage = cdr["ajcc_pathologic_tumor_stage"].astype(str).str.extract(r"Stage (IV|III|II|I)")[0]
    df = pd.DataFrame({"os": cdr["OS"], "time": cdr["OS.time"],
                       "age": cdr["age_at_initial_pathologic_diagnosis"],
                       "stage": stage.map({"I": 1, "II": 2, "III": 3, "IV": 4})}).dropna()
    df = df[df.index.str.endswith("-01") & df.index.isin(rna.columns) & (df.time > 0)]
    rows = []
    cph = CoxPHFitter()
    for g in genes:
        if g not in rna.index:
            continue
        x = rna.loc[g, df.index]
        if x.std() == 0:
            continue
        d = df.assign(gene=(x - x.mean()) / x.std())
        try:
            cph.fit(d, duration_col="time", event_col="os")
            hr, p = float(np.exp(cph.params_["gene"])), float(cph.summary.loc["gene", "p"])
        except Exception:
            hr, p = np.nan, np.nan
        c = concordance_index(d.time, -d.gene, d.os)
        rows.append({"gene": g, "cox_hr_per_sd": hr, "cox_p": p, "cindex_high_is_risk": c})
    out = pd.DataFrame(rows).set_index("gene")
    out["survival_n"], out["survival_events"] = int(len(df)), int(df.os.sum())
    return out


# ----------------------------------------------------------------------------- assemble
def run():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    obs = pd.read_csv(C.RESULTS / "ladder" / "observed.csv")
    best = obs.loc[obs.groupby("subtype").mean_effect.idxmin(), ["subtype", "arm", "mean_effect"]]
    best = best.set_index("subtype")
    print("  best arm per subtype (most negative mean gene effect):\n" + best.to_string())

    lists = pd.read_csv(C.RESULTS / "nominate" / "lists_discovery.csv")
    lists = lists[lists.k == C.K]
    cand = pd.concat([lists[(lists.subtype == s) & (lists.arm == best.loc[s, "arm"])] for s in C.SUBTYPES])
    genes = sorted(set(cand.gene))

    t1 = pd.read_csv(C.RESULTS / "truth" / "t1_gene_effect.csv", index_col=0)
    t2 = pd.read_csv(C.RESULTS / "truth" / "t2_clinical_precedence.csv").set_index("gene")
    tract = pd.read_csv(C.RESULTS / "truth" / "tractability.csv", index_col=0)
    gtex = pd.read_csv(C.RESULTS / "truth" / "gtex.csv", index_col=0)
    gate = pd.read_csv(C.RESULTS / "scores" / "gate.csv", index_col=0)
    rep_path = C.RESULTS / "replicate" / "lists_krug.csv"
    rep = pd.read_csv(rep_path) if rep_path.exists() else None
    print("  atlas specificity …", flush=True)
    atlas = atlas_specificity(genes)
    print("  survival …", flush=True)
    surv = tcga_survival(genes)

    rows = []
    for r in cand.itertuples(index=False):
        g, s = r.gene, r.subtype
        row = {"subtype": s, "arm": r.arm, "rank": r.rank, "gene": g, "nomination_score": r.score,
               "mean_effect_subtype": t1.loc[g, f"mean_{s}"],
               "n_dependent_lines": t1.loc[g, f"n_dependent_{s}"],
               "selectivity": t1.loc[g, f"selectivity_{s}"], "hit_subtype": bool(t1.loc[g, f"hit_{s}"]),
               "common_essential": bool(t1.loc[g, "common_essential"]),
               "clinical_stage": t2.max_stage.get(g, 0), "clinical_drugs": t2.n_drugs.get(g, 0),
               "tumour_vs_normal_g": gate.loc[g, "tumour_vs_normal_g"]}
        for col in tract.columns:
            row[f"tractability_{col.lower().replace(' ', '_')}"] = tract[col].get(g, "")
        for col in gtex.columns:
            row[col] = gtex[col].get(g, np.nan)
        for col in atlas.columns:
            row[col] = atlas[col].get(g, np.nan)
        for col in surv.columns:
            row[col] = surv[col].get(g, np.nan)
        if rep is not None:
            row["in_krug_list_same_arm"] = g in set(rep[(rep.arm == r.arm) & (rep.subtype == s)].gene)
        rows.append(row)
    d = pd.DataFrame(rows).sort_values(["subtype", "rank"])
    d.to_csv(OUT / "candidates.csv", index=False)
    top = d[d["rank"] <= TOP]
    top.to_csv(OUT / "top10.csv", index=False)

    cols = ["gene", "mean_effect_subtype", "n_dependent_lines", "hit_subtype", "clinical_stage",
            "atlas_malignant_ratio", "cox_hr_per_sd", "cox_p"]
    if rep is not None:
        cols.append("in_krug_list_same_arm")
    for s in C.SUBTYPES:
        print(f"  {s} (arm {best.loc[s, 'arm']}):\n"
              + top[top.subtype == s][cols].round(3).to_string(index=False))

    summary = {"best_arm": best.arm.to_dict(), "n_candidates": int(len(d)),
               "survival_n": int(surv.survival_n.iloc[0]),
               "survival_events": int(surv.survival_events.iloc[0]),
               "atlas_n_malignant": int(atlas.atlas_n_malignant.iloc[0]),
               "atlas_n_other": int(atlas.atlas_n_other.iloc[0])}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    P.log_run("dossier", {"top": TOP, "atlas": ATLAS_PATH},
              [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*"))], time.time() - t0)
    print(f"done in {time.time() - t0:,.0f}s")


if __name__ == "__main__":
    sys.exit(run())
