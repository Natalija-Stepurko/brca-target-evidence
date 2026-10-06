"""Stage `ladder`: every list against the ladder (DESIGN §2, §7) and the verdicts on P1-P4 and P6.

Endpoints per list (arm x subtype):
  primary    mean Chronos gene effect of the list in the subtype-matched lines (more negative = better)
  secondary  precision@K against the hit set H(s); count of T2 clinical-precedence genes
Rungs:
  matched-random floor   1,000 lists, each gene swapped for one from the same (PubMed, expression,
                         protein-detection) stratum
  label-permutation null 1,000 shuffles of the subtype labels; every arm re-nominated (MOFA+ re-selects
                         its factor from the fitted model, which never saw the labels)
  paired bootstrap       1,000 resamples of tumours, stratified by subtype; every arm re-nominated;
                         percentile interval on the difference between arms
"""
import gzip
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
from bte.stages import scores as S

OUT = C.RESULTS / "ladder"
RAW = C.DATA / "raw"
MOFA_BOOT = 200          # MOFA+ refits per bootstrap resample; see DESIGN §12
SMOKE = os.environ.get("BTE_SMOKE") == "1"          # tiny counts, separate output dir: checks the code path
N_RANDOM, N_PERM, N_BOOT, N_MOFA = (20, 10, 10, 2) if SMOKE else (C.N_MATCHED_RANDOM, C.N_PERMUTATIONS,
                                                                 C.N_BOOTSTRAP, MOFA_BOOT)
if SMOKE:
    OUT = C.RESULTS / "smoke" / "ladder"


# ----------------------------------------------------------------------------- truth lookups
class Truth:
    def __init__(self):
        t = pd.read_csv(C.RESULTS / "truth" / "t1_gene_effect.csv", index_col=0)
        self.effect = {s: t[f"mean_{s}"] for s in C.SUBTYPES}
        self.hits = {s: set(t.index[t[f"hit_{s}"]]) for s in C.SUBTYPES}
        p = pd.read_csv(C.RESULTS / "truth" / "t2_clinical_precedence.csv")
        self.t2 = set(p.gene)

    def score(self, genes, s) -> dict:
        genes = list(genes)
        eff = self.effect[s].reindex(genes)
        return {"mean_effect": float(eff.mean()), "n_scored": int(eff.notna().sum()),
                "precision": float(np.mean([g in self.hits[s] for g in genes])) if genes else np.nan,
                "n_t2": int(sum(g in self.t2 for g in genes))}


# ----------------------------------------------------------------------------- matched strata
def pubmed_counts(genes) -> pd.Series:
    info = pd.read_csv(RAW / "ncbi_Homo_sapiens.gene_info.gz", sep="\t", usecols=["GeneID", "Symbol"])
    with gzip.open(RAW / "ncbi_gene2pubmed.gz", "rt") as f:
        g2p = pd.read_csv(f, sep="\t", usecols=["#tax_id", "GeneID"])
    g2p = g2p[g2p["#tax_id"] == 9606]
    n = g2p.groupby("GeneID").size()
    sym = info.set_index("GeneID").Symbol
    counts = n.rename(index=sym).groupby(level=0).sum()
    return counts.reindex(genes).fillna(0).astype(int)


def strata(universe: pd.DataFrame, rna_mean: pd.Series) -> pd.DataFrame:
    """Joint stratum per gene: PubMed quintile x expression quintile x protein-detection tercile; strata
    with fewer than MIN_STRATUM genes are merged with the neighbouring expression stratum."""
    u = universe.set_index("symbol")
    d = pd.DataFrame({"pubmed": pubmed_counts(u.index), "expr": rna_mean.reindex(u.index),
                      "det": u.protein_detection})
    d = d.fillna(d.median(numeric_only=True))               # a gene missing from a layer takes the median
    d["q_pub"] = pd.qcut(d.pubmed.rank(method="first"), C.MATCH_STRATA["pubmed"], labels=False)
    d["q_expr"] = pd.qcut(d.expr.rank(method="first"), C.MATCH_STRATA["expression"], labels=False)
    d["q_det"] = pd.qcut(d.det.rank(method="first"), C.MATCH_STRATA["protein_detection"], labels=False)
    d["stratum"] = d.q_pub.astype(str) + "-" + d.q_expr.astype(str) + "-" + d.q_det.astype(str)
    sizes = d.stratum.value_counts()
    for st in sizes.index[sizes < C.MIN_STRATUM]:
        p, e, t = map(int, st.split("-"))
        e2 = e + 1 if e + 1 < C.MATCH_STRATA["expression"] else e - 1
        d.loc[d.stratum == st, "stratum"] = f"{p}-{e2}-{t}"
    return d


def matched_floor(genes, strat: pd.DataFrame, gate: pd.Series, truth: Truth, s: str, rng, n=None):
    """n random lists: each nominated gene replaced by a random gene from the same stratum (gate-passing,
    like the nominated gene itself); returns the endpoint distributions."""
    n = n or N_RANDOM
    pool = {st: np.array(list(idx)) for st, idx in strat[strat.index.map(gate).fillna(False).astype(bool)]
            .groupby("stratum").groups.items()}
    sts = strat.stratum.reindex(genes).values
    out = {"mean_effect": [], "precision": [], "n_t2": []}
    for _ in range(n):
        draw = [rng.choice(pool[st]) for st in sts]
        sc = truth.score(draw, s)
        for k in out:
            out[k].append(sc[k])
    return {k: np.array(v) for k, v in out.items()}


# ----------------------------------------------------------------------------- re-nomination
def nominate_all(co: K.Cohort, labels: pd.Series, cis: pd.Series, gate: pd.Series, mofa: A.Mofa | None,
                 k=C.K) -> pd.DataFrame:
    tab = S.score_cohort(co, labels, cis)
    lists = A.stack_lists(A.percentiles(tab), gate, k)
    if mofa is not None:
        lists = pd.concat([lists, mofa.lists(labels, gate, k)], ignore_index=True)
    return lists


def endpoint_table(lists: pd.DataFrame, truth: Truth) -> pd.DataFrame:
    rows = []
    for (arm, s), d in lists.groupby(["arm", "subtype"], sort=False):
        rows.append({"arm": arm, "subtype": s, "n": len(d), **truth.score(d.gene, s)})
    return pd.DataFrame(rows)


def run():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    genes = list(K.universe().symbol)
    universe = K.universe()
    truth = Truth()
    gate = pd.read_csv(C.RESULTS / "scores" / "gate.csv", index_col=0)["gate"]
    co = K.load_discovery(genes)
    cis = S.cis_correlation(co.cnv, co.rna)
    with open(C.RESULTS / "nominate" / "mofa_discovery.pkl", "rb") as f:
        mofa = pickle.load(f)
    lists = pd.read_csv(C.RESULTS / "nominate" / "lists_discovery.csv")
    lists = lists[lists.k == C.K]

    # 1. observed endpoints
    obs = endpoint_table(lists, truth)
    obs.to_csv(OUT / "observed.csv", index=False)
    print("  observed primary endpoint (mean gene effect in matched lines):\n"
          + obs.pivot(index="arm", columns="subtype", values="mean_effect").round(3).to_string())

    # 2. matched-random floor
    strat = strata(universe, co.rna.mean(axis=1))
    strat.to_csv(OUT / "strata.csv")
    floor_rows, floor_dist = [], {}
    for r in obs.itertuples(index=False):
        g = lists[(lists.arm == r.arm) & (lists.subtype == r.subtype)].gene
        dist = matched_floor(g, strat, gate, truth, r.subtype, rng)
        floor_dist[(r.arm, r.subtype)] = dist["mean_effect"]
        floor_rows.append({"arm": r.arm, "subtype": r.subtype,
                           "floor_pct_mean_effect": float(np.mean(dist["mean_effect"] <= r.mean_effect)),
                           "floor_mean_effect_median": float(np.median(dist["mean_effect"])),
                           "floor_mean_effect_p05": float(np.percentile(dist["mean_effect"], 5)),
                           "floor_pct_precision": float(np.mean(dist["precision"] >= r.precision)),
                           "floor_precision_median": float(np.median(dist["precision"])),
                           "floor_pct_n_t2": float(np.mean(dist["n_t2"] >= r.n_t2)),
                           "floor_n_t2_median": float(np.median(dist["n_t2"]))})
    floor = pd.DataFrame(floor_rows)
    floor.to_csv(OUT / "matched_floor.csv", index=False)
    np.savez_compressed(OUT / "floor_dist.npz", **{f"{a}|{s}": v for (a, s), v in floor_dist.items()})
    print("  floor percentile of the observed mean effect (lower = better than more random lists):\n"
          + floor.pivot(index="arm", columns="subtype", values="floor_pct_mean_effect").round(3).to_string())

    # 3. label-permutation null
    perm = []
    for i in range(N_PERM):
        lab = pd.Series(rng.permutation(co.subtype.values), index=co.subtype.index)
        e = endpoint_table(nominate_all(co, lab, cis, gate, mofa), truth).assign(perm=i)
        perm.append(e)
        if (i + 1) % 100 == 0:
            print(f"    permutation {i + 1}/{N_PERM}  {time.time() - t0:,.0f}s", flush=True)
    perm = pd.concat(perm, ignore_index=True)
    perm.to_csv(OUT / "permutation_null.csv", index=False)
    null = (perm.groupby(["arm", "subtype"]).mean_effect
            .agg(null_median="median", null_p05=lambda x: np.percentile(x, 5)).reset_index())
    null = null.merge(obs[["arm", "subtype", "mean_effect"]])
    null["null_pct_mean_effect"] = [
        float(np.mean(perm[(perm.arm == a) & (perm.subtype == s)].mean_effect <= m))
        for a, s, m in zip(null.arm, null.subtype, null.mean_effect, strict=True)]
    null.to_csv(OUT / "permutation_summary.csv", index=False)

    # 4. paired bootstrap over tumours (stack arms: N_BOOTSTRAP; MOFA+ refits: MOFA_BOOT)
    boot = []
    by_sub = {s: co.subtype.index[co.subtype == s] for s in C.SUBTYPES}
    for i in range(N_BOOT):
        idx = np.concatenate([rng.choice(ids, len(ids), replace=True) for ids in by_sub.values()])
        sub = co.subtype.loc[idx]
        sub.index = [f"{x}#{j}" for j, x in enumerate(idx)]            # duplicates need unique names
        rs = K.Cohort("boot", sub, co.rna[idx].set_axis(sub.index, axis=1),
                      co.protein[idx].set_axis(sub.index, axis=1), co.cnv[idx].set_axis(sub.index, axis=1),
                      co.mutation[idx].set_axis(sub.index, axis=1))
        m = None
        if i < N_MOFA:
            m = A.Mofa({"rna": rs.rna, "protein": rs.protein, "cnv": rs.cnv})
        e = endpoint_table(nominate_all(rs, sub, cis, gate, m), truth).assign(boot=i)
        boot.append(e)
        if (i + 1) % 50 == 0:
            print(f"    bootstrap {i + 1}/{N_BOOT}  {time.time() - t0:,.0f}s", flush=True)
    boot = pd.concat(boot, ignore_index=True)
    boot.to_csv(OUT / "bootstrap.csv", index=False)

    # 5. arm differences with intervals, and R-all vs R+P (P6) from the observed R-all list
    diffs = []
    wide = boot.pivot_table(index=["boot", "subtype"], columns="arm", values="mean_effect")
    for s in C.SUBTYPES:
        w = wide.xs(s, level="subtype")
        for a, b in (("R+P", "R"), ("R+D", "R"), ("R+P+D", "R"), ("MOFA+", "R+P+D"), ("R+P+D", "R+P")):
            if a not in w or b not in w:
                continue
            d = (w[a] - w[b]).dropna()
            o = obs.set_index(["arm", "subtype"]).mean_effect
            diffs.append({"subtype": s, "arm": a, "vs": b, "observed_diff": float(o[(a, s)] - o[(b, s)]),
                          "boot_lo": float(np.percentile(d, 2.5)), "boot_hi": float(np.percentile(d, 97.5)),
                          "n_boot": int(len(d))})
    diffs = pd.DataFrame(diffs)
    diffs.to_csv(OUT / "arm_differences.csv", index=False)
    print("  arm differences in mean gene effect (negative favours the first arm):\n"
          + diffs.round(3).to_string(index=False))

    if SMOKE:
        print(f"smoke run finished in {time.time() - t0:,.0f}s; nothing logged"); return
    P.log_run("ladder", {"n_matched_random": N_RANDOM, "n_permutations": N_PERM,
                         "n_bootstrap": N_BOOT, "mofa_bootstrap_refits": N_MOFA, "seed": 0},
              [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*.csv"))], time.time() - t0)
    print(f"done in {time.time() - t0:,.0f}s")


if __name__ == "__main__":
    sys.exit(run())
