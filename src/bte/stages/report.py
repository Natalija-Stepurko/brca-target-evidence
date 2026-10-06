"""Stage `report`: verdicts on P1-P6 from results/ladder and results/replicate, and the figures.

Every verdict is a boolean derived here from the tables, with the rule written next to it; the page
reads results/report/summary.json and never types a verdict in.
"""
import json
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from bte import config as C  # noqa: E402
from bte import provenance as P  # noqa: E402

OUT = C.RESULTS / "report"
LAD, REP, DOS = C.RESULTS / "ladder", C.RESULTS / "replicate", C.RESULTS / "dossier"
INK, MUTED, RNA, PROT, DNA, TRUTH = "#16191D", "#5B646E", "#2D5BD1", "#C06014", "#6E7880", "#0E7C7B"
ARM_COL = {"R": RNA, "R+P": PROT, "R+D": DNA, "R+P+D": "#7A3E9D", "MOFA+": "#B8860B", "R-all": "#5B8DD6"}
ARMS = ["R", "R+P", "R+D", "R+P+D", "MOFA+", "R-all"]
SUB_LABEL = {"basal": "basal-like", "her2": "HER2-enriched", "luminal": "luminal"}
plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": "#B9C0C6", "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "figure.dpi": 150})


# ----------------------------------------------------------------------------- verdicts
def verdicts(obs, floor, diffs, rep_overlap, rep_end) -> dict:
    o = obs.set_index(["arm", "subtype"]).mean_effect
    f = floor.set_index(["arm", "subtype"]).floor_pct_mean_effect
    d = diffs.set_index(["subtype", "arm", "vs"])
    seoi = C.SEOI_GENE_EFFECT
    v = {}

    # P1: every arm beats its floor (percentile <= 5%) in at least two of three subtypes
    per_arm = {a: int(sum(f[(a, s)] <= 0.05 for s in C.SUBTYPES)) for a in ARMS}
    v["P1"] = {"rule": "each arm: floor percentile of the mean gene effect <= 0.05 in >= 2 of 3 subtypes",
               "subtypes_beating_floor": per_arm, "holds": all(n >= 2 for n in per_arm.values()),
               "arms_holding": [a for a, n in per_arm.items() if n >= 2]}

    def arm_diff(a, b, s):
        r = d.loc[(s, a, b)]
        return {"observed": float(r.observed_diff), "lo": float(r.boot_lo), "hi": float(r.boot_hi)}

    # P2: R+P does not beat R by the SEOI (negative difference favours R+P)
    p2 = {s: arm_diff("R+P", "R", s) for s in C.SUBTYPES}
    v["P2"] = {"rule": f"(R+P minus R) in mean gene effect > -{seoi} in every subtype; interval reported",
               "per_subtype": p2, "holds": all(x["observed"] > -seoi for x in p2.values()),
               "protein_hurts": {s: bool(x["lo"] > 0) for s, x in p2.items()}}
    # P3: R+D adds no more over R than R+P does
    p3 = {s: {"R+D_minus_R": arm_diff("R+D", "R", s)["observed"], "R+P_minus_R": p2[s]["observed"]}
          for s in C.SUBTYPES}
    v["P3"] = {"rule": f"(R+D minus R) >= (R+P minus R) - {seoi} in every subtype", "per_subtype": p3,
               "holds": all(x["R+D_minus_R"] >= x["R+P_minus_R"] - seoi for x in p3.values())}
    # P4: MOFA+ does not beat R+P+D
    p4 = {s: arm_diff("MOFA+", "R+P+D", s) for s in C.SUBTYPES}
    v["P4"] = {"rule": f"(MOFA+ minus R+P+D) > -{seoi} in every subtype", "per_subtype": p4,
               "holds": all(x["observed"] > -seoi for x in p4.values())}
    # P6: R-all scores at least as well as R+P
    p6 = {s: {"R-all": float(o[("R-all", s)]), "R+P": float(o[("R+P", s)])} for s in C.SUBTYPES}
    v["P6"] = {"rule": f"mean gene effect of R-all <= that of R+P + {seoi} in every subtype",
               "per_subtype": p6, "holds": all(x["R-all"] <= x["R+P"] + seoi for x in p6.values())}
    # P5: replication
    if rep_overlap is not None:
        ro = rep_overlap.set_index(["arm", "subtype"])
        above = {a: int(sum(bool(ro.above_null[(a, s)]) for s in C.SUBTYPES)) for a in ARMS if a != "R-all"}
        prot_ok = all(ro.overlap[("R+P", s)] >= ro.overlap[("R", s)] - C.SEOI_PRECISION for s in C.SUBTYPES)
        v["P5"] = {"rule": "each arm's discovery-Krug overlap above the permutation null p95 in every "
                           f"subtype; overlap(R+P) >= overlap(R) - {C.SEOI_PRECISION}",
                   "subtypes_above_null": above, "protein_replicates_no_worse": prot_ok,
                   "holds": all(n == 3 for n in above.values()) and prot_ok}
    return v


# ----------------------------------------------------------------------------- figures
def fig_floor(obs, floor, dist):
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, s in zip(axes, C.SUBTYPES, strict=True):
        for i, a in enumerate(ARMS):
            key = f"{a}|{s}"
            if key in dist:
                v = dist[key]
                ax.plot([i - .28, i + .28], [np.median(v)] * 2, color="#B9C0C6", lw=1.2)
                ax.fill_between([i - .28, i + .28], np.percentile(v, 5), np.percentile(v, 95),
                                color="#E6EAEC", zorder=0)
            m = obs[(obs.arm == a) & (obs.subtype == s)].mean_effect.iloc[0]
            ax.scatter([i], [m], color=ARM_COL[a], s=46, zorder=3)
        ax.set_xticks(range(len(ARMS)), ARMS, rotation=35, ha="right")
        ax.set_title(SUB_LABEL[s], fontsize=10, color=INK)
        ax.axhline(0, color="#DDE1E4", lw=.8)
    axes[0].set_ylabel("mean gene effect of the top-50\nin subtype-matched lines")
    fig.text(0.5, -0.02, "dot = observed; grey band = 5th-95th percentile of 1,000 matched random lists, "
             "line = their median. More negative = more dependency.", ha="center", color=MUTED, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(OUT / "fig_floor.png", bbox_inches="tight"); plt.close(fig)


def fig_diffs(diffs):
    pairs = [("R+P", "R"), ("R+D", "R"), ("R+P+D", "R"), ("MOFA+", "R+P+D")]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.2), sharex=True)
    for ax, s in zip(axes, C.SUBTYPES, strict=True):
        for j, (a, b) in enumerate(pairs):
            r = diffs[(diffs.subtype == s) & (diffs.arm == a) & (diffs.vs == b)]
            if r.empty:
                continue
            r = r.iloc[0]
            ax.plot([r.boot_lo, r.boot_hi], [j, j], color=ARM_COL[a], lw=2)
            ax.scatter([r.observed_diff], [j], color=ARM_COL[a], s=40, zorder=3)
        ax.axvline(0, color=INK, lw=.8)
        ax.axvspan(-C.SEOI_GENE_EFFECT, C.SEOI_GENE_EFFECT, color="#EEF1F2", zorder=0)
        ax.set_yticks(range(len(pairs)), [f"{a} − {b}" for a, b in pairs])
        ax.invert_yaxis()
        ax.set_title(SUB_LABEL[s], fontsize=10, color=INK)
    axes[1].set_xlabel("difference in mean gene effect (negative favours the first arm); "
                       "95% paired bootstrap interval; shaded = smallest effect of interest")
    fig.tight_layout()
    fig.savefig(OUT / "fig_diffs.png", bbox_inches="tight"); plt.close(fig)


def fig_replication(ro):
    arms = [a for a in ARMS if a != "R-all"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.2), sharey=True)
    for ax, s in zip(axes, C.SUBTYPES, strict=True):
        d = ro[ro.subtype == s].set_index("arm").reindex(arms)
        ax.bar(range(len(arms)), d.overlap, color=[ARM_COL[a] for a in arms], width=.6)
        ax.plot(range(len(arms)), d.null_p95, "_", color=INK, markersize=18, mew=1.5)
        ax.set_xticks(range(len(arms)), arms, rotation=35, ha="right")
        ax.set_title(SUB_LABEL[s], fontsize=10, color=INK)
    axes[0].set_ylabel("genes shared by the discovery\nand Krug top-50 lists")
    fig.text(0.5, -0.02, "bar = observed overlap; tick = 95th percentile of the overlap when Krug's subtype "
             "labels are shuffled 1,000 times", ha="center", color=MUTED, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(OUT / "fig_replication.png", bbox_inches="tight"); plt.close(fig)


def run():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    obs = pd.read_csv(LAD / "observed.csv")
    floor = pd.read_csv(LAD / "matched_floor.csv")
    diffs = pd.read_csv(LAD / "arm_differences.csv")
    dist = dict(np.load(LAD / "floor_dist.npz"))
    ro = pd.read_csv(REP / "overlap.csv") if (REP / "overlap.csv").exists() else None
    re_ = pd.read_csv(REP / "endpoints_krug.csv") if (REP / "endpoints_krug.csv").exists() else None

    v = verdicts(obs, floor, diffs, ro, re_)
    for k, x in v.items():
        print(f"  {k}: {'holds' if x['holds'] else 'fails'}   ({x['rule']})")
    fig_floor(obs, floor, dist)
    fig_diffs(diffs)
    if ro is not None:
        fig_replication(ro)

    summary = {"verdicts": v,
               "observed": obs.to_dict("records"), "floor": floor.to_dict("records"),
               "differences": diffs.to_dict("records"),
               "replication": ro.to_dict("records") if ro is not None else None,
               "replication_endpoints": re_.to_dict("records") if re_ is not None else None,
               "dossier": (json.loads((DOS / "summary.json").read_text())
                           if (DOS / "summary.json").exists() else None),
               "n_matched_random": C.N_MATCHED_RANDOM, "n_permutations": C.N_PERMUTATIONS,
               "n_bootstrap": C.N_BOOTSTRAP, "seoi_gene_effect": C.SEOI_GENE_EFFECT}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    P.log_run("report", {}, [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*"))], time.time() - t0)
    print(f"done in {time.time() - t0:,.0f}s")


if __name__ == "__main__":
    sys.exit(run())
