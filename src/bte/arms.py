"""The nomination arms (DESIGN §5.3-§5.4): from a per-layer score table to a ranked gene list per subtype.

Rank stacking: each layer statistic becomes a rank percentile within the universe (per subtype), and an
arm's score is the mean of its layers' percentiles, ignoring layers where the gene has no value. The
gate is applied after ranking. MOFA+ is fitted once per cohort (it does not see the labels); only the
factor choice depends on labels, so permutations and resamples re-select the factor cheaply.
"""
import warnings

import numpy as np
import pandas as pd

from bte import config as C

LAYERS = {"R": ["rna"], "R+P": ["rna", "protein"], "R+D": ["rna", "dna"], "R+P+D": ["rna", "protein", "dna"]}
STACK_ARMS = tuple(LAYERS)


def percentiles(scores: pd.DataFrame) -> pd.DataFrame:
    """Per subtype, rank percentile (0-1, higher = more evidence) of each layer statistic."""
    out = scores[["gene", "subtype"]].copy()
    grp = scores.groupby("subtype", sort=False)
    out["rna"] = grp["rna_g"].rank(pct=True)
    if "protein_g" in scores:
        out["protein"] = grp["protein_g"].rank(pct=True)
    if "cnv_g" in scores or "mut_diff" in scores:
        parts = []
        if "cnv_g" in scores:
            parts.append(grp["cnv_g"].rank(pct=True))
        if "mut_diff" in scores:
            parts.append(grp["mut_diff"].rank(pct=True))
        out["dna"] = pd.concat(parts, axis=1).max(axis=1)        # the larger of the two DNA percentiles
    return out


def stack_lists(pct: pd.DataFrame, gate: pd.Series, k: int = C.K, arms=STACK_ARMS) -> pd.DataFrame:
    """Top-k per arm and subtype after the gate. Returns long table: arm, subtype, rank, gene, score."""
    rows = []
    for arm in arms:
        cols = [c for c in LAYERS[arm] if c in pct]
        if cols != LAYERS[arm]:
            continue                                                # a layer this cohort lacks
        score = pct[cols].mean(axis=1)
        d = pct[["gene", "subtype"]].assign(score=score)
        d = d[d.gene.map(gate).fillna(False).astype(bool) & d.score.notna()]
        for s, sd in d.groupby("subtype", sort=False):
            top = sd.nlargest(k, "score").reset_index(drop=True)
            rows.append(pd.DataFrame({"arm": arm, "subtype": s, "rank": np.arange(1, len(top) + 1),
                                      "gene": top.gene.values, "score": top.score.values}))
    return pd.concat(rows, ignore_index=True)


def ungated_rank(pct: pd.DataFrame, arm: str, subtype: str, gene: str) -> int:
    """Rank of a gene in an arm's ranking before the gate; for the positive controls (§6.3)."""
    cols = LAYERS[arm]
    d = pct[pct.subtype == subtype]
    score = d[cols].mean(axis=1)
    return int((score > score[d.gene == gene].iloc[0]).sum() + 1)


# ----------------------------------------------------------------------------- MOFA+
class Mofa:
    """MOFA+ on RNA, protein and copy number: weights per view (genes x factors), factor scores per tumour."""

    def __init__(self, views: dict[str, pd.DataFrame], n_factors=C.MOFA_FACTORS, seed=C.MOFA_SEED):
        from mofapy2.run.entry_point import entry_point
        samples = list(next(iter(views.values())).columns)
        data = []
        self.genes = {}
        for name, m in views.items():
            m = m.reindex(columns=samples)
            m = m[m.notna().mean(axis=1) >= 0.5]                      # features seen in half the tumours
            m = m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0)
            self.genes[name] = list(m.index)
            data.append([m.T.values])                                  # one group, samples x features
        ent = entry_point()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ent.set_data_options(scale_views=True)
            ent.set_data_matrix(data, views_names=list(views), groups_names=["all"], samples_names=[samples],
                                features_names=[[f"{v}:{g}" for g in self.genes[v]] for v in views])
            ent.set_model_options(factors=n_factors, spikeslab_weights=True, ard_weights=True)
            ent.set_train_options(iter=1000, convergence_mode="fast", seed=seed, verbose=False, quiet=True)
            ent.build()
            ent.run()
        model = ent.model
        exp = model.getExpectations()
        self.factors = pd.DataFrame(exp["Z"]["E"], index=samples)
        self.weights = {v: pd.DataFrame(exp["W"][i]["E"], index=self.genes[v]) for i, v in enumerate(views)}
        self.samples = samples

    def rank_genes(self, labels: pd.Series, subtype: str) -> pd.Series:
        """Choose the factor most correlated with membership of `subtype`; rank genes by the sum over views
        of |standardised weight|, sign agreeing with 'up in subtype'."""
        y = (labels.reindex(self.samples) == subtype).astype(float).values
        r = np.array([np.corrcoef(self.factors[f].values, y)[0, 1] for f in self.factors.columns])
        f = int(np.nanargmax(np.abs(r)))
        sign = np.sign(r[f]) or 1.0
        total = None
        for w in self.weights.values():
            wf = w[f] * sign
            z = wf / wf.std()
            total = z if total is None else total.add(z, fill_value=0)
        up = total.where(total > 0, 0)                                   # genes pointing the subtype's way
        return up.sort_values(ascending=False)

    def lists(self, labels: pd.Series, gate: pd.Series, k: int = C.K) -> pd.DataFrame:
        rows = []
        for s in C.SUBTYPES:
            if (labels == s).sum() < 3:
                continue
            r = self.rank_genes(labels, s)
            r = r[r.index.map(gate).fillna(False).astype(bool) & (r > 0)]
            top = r.head(k)
            rows.append(pd.DataFrame({"arm": "MOFA+", "subtype": s, "rank": np.arange(1, len(top) + 1),
                                      "gene": top.index.values, "score": top.values}))
        return pd.concat(rows, ignore_index=True)
