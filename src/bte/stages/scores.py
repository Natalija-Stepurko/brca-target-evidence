"""Stage `scores`: per-layer evidence for every universe gene and subtype (DESIGN §5.1) and the
tumour-vs-normal gate (§5.2), for the discovery cohort, the TCGA-all cohort and the Krug cohort.

The functions here are reused by `nominate`, `ladder` and `replicate` with permuted labels or
resampled tumours, so they take a cohort and a label vector and return a table.
"""
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from bte import cohorts as K
from bte import config as C
from bte import provenance as P

OUT = C.RESULTS / "scores"


def cis_correlation(cnv: pd.DataFrame, rna: pd.DataFrame) -> pd.Series:
    """Spearman correlation of copy number and expression across tumours, per gene."""
    rho = pd.Series(np.nan, index=cnv.index)
    both = cnv.notna().sum(axis=1) >= 10
    for g in cnv.index[both]:
        c, r = cnv.loc[g], rna.loc[g]
        ok = c.notna() & r.notna()
        if ok.sum() >= 10 and c[ok].nunique() > 1 and r[ok].nunique() > 1:
            rho[g] = spearmanr(c[ok], r[ok]).statistic
    return rho


def gate(tcga_all: K.Cohort) -> pd.Series:
    """Up in tumour vs adjacent normal in TCGA RNA-seq: Hedges' g >= GATE_HEDGES_G. One gate for every arm."""
    both = pd.concat([tcga_all.rna, tcga_all.normals], axis=1)
    mask = np.array([c in tcga_all.rna.columns for c in both.columns])
    g = K.hedges_g(both, mask)
    return (g >= C.GATE_HEDGES_G).rename("gate"), g.rename("tumour_vs_normal_g")


def score_cohort(co: K.Cohort, labels: pd.Series, cis: pd.Series | None = None) -> pd.DataFrame:
    """Long table: gene, subtype, and one column per layer statistic."""
    labels = labels.reindex(co.rna.columns)
    rows = []
    for s in C.SUBTYPES:
        mask = (labels == s).values
        if mask.sum() < 3:
            continue
        d = pd.DataFrame({"gene": co.rna.index, "subtype": s})
        d["rna_g"] = K.hedges_g(co.rna, mask).values
        if co.protein is not None:
            d["protein_g"] = K.hedges_g(co.protein.reindex(columns=co.rna.columns), mask).values
        if co.cnv is not None:
            d["cnv_g"] = K.hedges_g(co.cnv.reindex(columns=co.rna.columns), mask).values
            if cis is not None:
                d["cis_rho"] = cis.reindex(co.rna.index).values
                d.loc[~(d.cis_rho >= C.CIS_SPEARMAN_MIN), "cnv_g"] = np.nan
        if co.mutation is not None:
            m = co.mutation.reindex(columns=co.rna.columns)
            d["mut_diff"] = (m.loc[:, mask].mean(axis=1) - m.loc[:, ~mask].mean(axis=1)).values
            d["mut_freq"] = m.loc[:, mask].mean(axis=1).values
        rows.append(d)
    return pd.concat(rows, ignore_index=True)


def run():
    t0 = time.time()
    genes = list(K.universe().symbol)
    OUT.mkdir(parents=True, exist_ok=True)

    tcga_all = K.load_tcga_all(genes)
    gt, g_tn = gate(tcga_all)
    pd.concat([g_tn, gt], axis=1).rename_axis("gene").to_csv(OUT / "gate.csv")
    print(f"  gate: {int(gt.sum()):,} of {len(gt):,} genes up in tumour vs normal "
          f"({tcga_all.rna.shape[1]} tumours, {tcga_all.normals.shape[1]} normals)")

    summary = {"gate_pass": int(gt.sum()), "gate_total": int(len(gt))}
    for co in (K.load_discovery(genes), tcga_all, K.load_krug(genes)):
        cis = cis_correlation(co.cnv, co.rna) if co.cnv is not None else None
        tab = score_cohort(co, co.subtype, cis)
        tab.to_csv(OUT / f"{co.name}.csv", index=False)
        counts = co.subtype.value_counts().to_dict()
        summary[co.name] = {"n": int(len(co.subtype)), "subtypes": counts, **co.notes,
                            "layers": list(co.layers())}
        if cis is not None:
            summary[co.name]["cis_rho_ge_min"] = int((cis >= C.CIS_SPEARMAN_MIN).sum())
        print(f"  {co.name:<10} n={len(co.subtype):>4}  {counts}  layers={list(co.layers())}")
    pd.Series(summary).to_json(OUT / "summary.json", indent=1)
    P.log_run("scores", {"gate_hedges_g": C.GATE_HEDGES_G, "cis_spearman_min": C.CIS_SPEARMAN_MIN},
              [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*"))], time.time() - t0)
    print(f"done in {time.time() - t0:,.0f}s")


if __name__ == "__main__":
    sys.exit(run())
