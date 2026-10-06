"""Stage `nominate`: the six arms' top-K lists per subtype on the discovery cohort (and R-all on TCGA),
the positive controls (DESIGN §6.3), and the fitted MOFA+ model saved for the ladder stage.
"""
import pickle
import sys
import time

import pandas as pd

from bte import arms as A
from bte import cohorts as K
from bte import config as C
from bte import provenance as P

OUT = C.RESULTS / "nominate"
SCORES = C.RESULTS / "scores"
CONTROLS = (("her2", "ERBB2"), ("luminal", "ESR1"))


def load_gate() -> pd.Series:
    return pd.read_csv(SCORES / "gate.csv", index_col=0)["gate"]


def run():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    gate = load_gate()
    genes = list(K.universe().symbol)

    disc = pd.read_csv(SCORES / "discovery.csv")
    pct = A.percentiles(disc)
    lists = [A.stack_lists(pct, gate, k) .assign(k=k) for k in (C.K, *C.K_SENSITIVITY)]

    # R-all: RNA alone on every TCGA-BRCA tumour with a subtype
    tall = pd.read_csv(SCORES / "tcga_all.csv")
    pct_all = A.percentiles(tall)
    for k in (C.K, *C.K_SENSITIVITY):
        lists.append(A.stack_lists(pct_all, gate, k, arms=("R",)).assign(arm="R-all", k=k))

    # MOFA+: fitted once on the discovery tumours (labels unseen); factor chosen per subtype
    co = K.load_discovery(genes)
    print("  fitting MOFA+ …", flush=True)
    t1 = time.time()
    mofa = A.Mofa({"rna": co.rna, "protein": co.protein, "cnv": co.cnv})
    print(f"  MOFA+ fitted in {time.time() - t1:,.0f}s on {len(mofa.samples)} tumours; features per view: "
          f"{ {v: len(g) for v, g in mofa.genes.items()} }")
    with open(OUT / "mofa_discovery.pkl", "wb") as f:
        pickle.dump(mofa, f)
    for k in (C.K, *C.K_SENSITIVITY):
        lists.append(mofa.lists(co.subtype, gate, k).assign(k=k))

    allk = pd.concat(lists, ignore_index=True)
    allk.to_csv(OUT / "lists_discovery.csv", index=False)
    main = allk[allk.k == C.K]
    sizes = main.groupby(["arm", "subtype"]).size().unstack()
    print("  list sizes at K=50:\n" + sizes.to_string())

    # positive controls on the RNA ranking before the gate (§6.3; see DESIGN §12 for why pre-gate)
    ctrl = []
    for s, gene in CONTROLS:
        r_pre = A.ungated_rank(pct, "R", s, gene)
        in_list = gene in set(main[(main.arm == "R") & (main.subtype == s)].gene)
        ctrl.append({"subtype": s, "gene": gene, "rank_pre_gate": r_pre, "in_gated_top50": in_list,
                     "passes_gate": bool(gate.get(gene, False))})
        print(f"  control {gene:<6} {s:<8} RNA rank before gate = {r_pre:>3}; gate = {gate.get(gene)}; "
              f"in gated top-50 = {in_list}")
    ctrl = pd.DataFrame(ctrl)
    ctrl.to_csv(OUT / "positive_controls.csv", index=False)
    assert (ctrl.rank_pre_gate <= C.K).all(), "a positive control is outside the pre-gate top-K"

    P.log_run("nominate", {"K": C.K, "K_sensitivity": list(C.K_SENSITIVITY), "mofa_factors": C.MOFA_FACTORS},
              [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*.csv"))], time.time() - t0)
    print(f"done in {time.time() - t0:,.0f}s")


if __name__ == "__main__":
    sys.exit(run())
