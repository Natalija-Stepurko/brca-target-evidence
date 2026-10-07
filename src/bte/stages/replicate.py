"""Stage `replicate`: the frozen pipeline on the Krug 2020 cohort (DESIGN §7, P5).

Every arm is re-nominated on the 117 Krug tumours with a subtype (MOFA+ fitted afresh on that cohort,
labels unseen). For each arm and subtype: overlap of the Krug top-50 with the discovery top-50,
against the distribution of overlaps when Krug's subtype labels are permuted; and the held-out
endpoints of the Krug lists with their matched-random floor.
"""
import os
import pickle
import sys
import time

import numpy as np
import pandas as pd

from bte import arms as A
from bte import cohorts as K
from bte import config as C
from bte import provenance as P
from bte.stages import ladder as L
from bte.stages import scores as S

OUT = C.RESULTS / "replicate"
SMOKE = os.environ.get("BTE_SMOKE") == "1"
N_PERM = 10 if SMOKE else C.N_PERMUTATIONS
if SMOKE:
    OUT = C.RESULTS / "smoke" / "replicate"


def overlaps(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (arm, s), d in a.groupby(["arm", "subtype"], sort=False):
        other = set(b[(b.arm == arm) & (b.subtype == s)].gene)
        rows.append({"arm": arm, "subtype": s, "overlap": len(set(d.gene) & other)})
    return pd.DataFrame(rows)


def run():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1)
    genes = list(K.universe().symbol)
    truth = L.Truth()
    gate = pd.read_csv(C.RESULTS / "scores" / "gate.csv", index_col=0)["gate"]
    disc = pd.read_csv(C.RESULTS / "nominate" / "lists_discovery.csv")
    disc = disc[(disc.k == C.K) & (disc.arm != "R-all")]

    co = K.load_krug(genes)
    cis = S.cis_correlation(co.cnv, co.rna)
    print(f"  Krug cohort: {len(co.subtype)} tumours {co.subtype.value_counts().to_dict()}; fitting MOFA+ …",
          flush=True)
    mofa = A.Mofa({"rna": co.rna, "protein": co.protein, "cnv": co.cnv})
    with open(OUT / "mofa_krug.pkl", "wb") as f:
        pickle.dump(mofa, f)
    rep = L.nominate_all(co, co.subtype, cis, gate, mofa)
    rep.to_csv(OUT / "lists_krug.csv", index=False)

    # observed overlap with the discovery lists
    obs = overlaps(rep, disc)

    # permutation null: shuffle Krug's labels, re-nominate, overlap with the fixed discovery lists
    null = []
    for i in range(N_PERM):
        lab = pd.Series(rng.permutation(co.subtype.values), index=co.subtype.index)
        null.append(overlaps(L.nominate_all(co, lab, cis, gate, mofa), disc).assign(perm=i))
        if (i + 1) % 100 == 0:
            print(f"    permutation {i + 1}/{N_PERM}  {time.time() - t0:,.0f}s", flush=True)
    null = pd.concat(null, ignore_index=True)
    null.to_csv(OUT / "overlap_null.csv", index=False)
    summ = (null.groupby(["arm", "subtype"]).overlap
            .agg(null_median="median", null_p95=lambda x: np.percentile(x, 95)).reset_index().merge(obs))
    summ["null_pct"] = [float(np.mean(null[(null.arm == a) & (null.subtype == s)].overlap >= o))
                        for a, s, o in zip(summ.arm, summ.subtype, summ.overlap, strict=True)]
    summ["above_null"] = summ.overlap > summ.null_p95
    summ.to_csv(OUT / "overlap.csv", index=False)
    txt = summ.overlap.astype(str) + " (" + summ.null_p95.round(0).astype(int).astype(str) + ")"
    print("  overlap of Krug and discovery top-50 (null p95 in brackets):\n"
          + summ.assign(txt=txt).pivot(index="arm", columns="subtype", values="txt").to_string())

    # held-out endpoints of the Krug lists, with their matched floor
    end = L.endpoint_table(rep, truth)
    strat = L.strata(K.universe(), co.rna.mean(axis=1))
    rows = []
    for r in end.itertuples(index=False):
        g = rep[(rep.arm == r.arm) & (rep.subtype == r.subtype)].gene
        dist = L.matched_floor(g, strat, gate, truth, r.subtype, rng, n=20 if SMOKE else C.N_MATCHED_RANDOM)
        rows.append({"arm": r.arm, "subtype": r.subtype, "mean_effect": r.mean_effect,
                     "precision": r.precision, "n_t2": r.n_t2,
                     "floor_pct_mean_effect": float(np.mean(dist["mean_effect"] <= r.mean_effect))})
    end = pd.DataFrame(rows)
    end.to_csv(OUT / "endpoints_krug.csv", index=False)
    print("  Krug lists, mean gene effect in matched lines:\n"
          + end.pivot(index="arm", columns="subtype", values="mean_effect").round(3).to_string())

    if SMOKE:
        print(f"smoke run finished in {time.time() - t0:,.0f}s"); return
    P.log_run("replicate", {"n_permutations": N_PERM, "seed": 1, "n_krug": int(len(co.subtype))},
              [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*.csv"))], time.time() - t0)
    print(f"done in {time.time() - t0:,.0f}s")


if __name__ == "__main__":
    sys.exit(run())
