"""Stage `audit`: agent-written target dossiers checked claim by claim (DESIGN §9; P7, P8).

Items: the top-10 genes per subtype from the best arm (30 candidates) and 30 decoys drawn from the same
matched strata that no arm nominated. For each gene and subtype the agent returns, in a fixed JSON
schema, a target score (0-10) and typed claims with a direction and a confidence. Two conditions:
closed-book (the model alone) and open-book (the model is given this repository's frozen evidence
table for that gene). Every typed claim is checked against the frozen tables and labelled supported,
contradicted or unsupported.

Needs ANTHROPIC_API_KEY in the environment and the `audit` extra (`uv sync --extra audit`). The model
and its knowledge cutoff are recorded in results/audit/summary.json.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

from bte import config as C
from bte import provenance as P

OUT = C.RESULTS / "audit"
MODEL = os.environ.get("BTE_AUDIT_MODEL", "claude-sonnet-5-5")
CUTOFF = os.environ.get("BTE_AUDIT_CUTOFF", "stated by the provider for the model above")
SUB_LABEL = {"basal": "basal-like", "her2": "HER2-enriched", "luminal": "luminal (ER-positive)"}
# thresholds that turn a frozen statistic into a yes/no the claim can be checked against
THRESH = {"overexpressed_in_subtype": ("rna_g", 0.5),          # Hedges' g, subtype vs other tumours
          "tumour_vs_normal": ("tumour_vs_normal_g", 0.5),
          "dependency_in_subtype_lines": ("mean_effect_subtype", -0.5),
          "clinical_precedence": ("clinical_stage", 1),
          "normal_tissue_restricted": ("gtex_ratio", 2.0)}   # breast or tumour-relevant: max other tissue low

SCHEMA = """Return only JSON of this shape:
{"gene": "...", "subtype": "...", "target_score": 0-10,
 "claims": [{"type": "overexpressed_in_subtype" | "tumour_vs_normal" | "dependency_in_subtype_lines" |
             "clinical_precedence" | "normal_tissue_restricted" | "other",
             "direction": "yes" | "no", "confidence": 0.0-1.0, "text": "one sentence"}],
 "rationale": "three sentences at most"}
Claim types mean: overexpressed_in_subtype = expression is higher in this subtype than in other breast
tumours; tumour_vs_normal = higher in breast tumours than in normal breast; dependency_in_subtype_lines =
CRISPR knockout reduces growth of breast cancer cell lines of this subtype (DepMap gene effect at or below
-0.5); clinical_precedence = a drug against this gene has reached a clinical trial for breast cancer;
normal_tissue_restricted = expression in other normal tissues is low (at most half the breast level)."""


def prompt(gene, subtype, evidence=None):
    p = (f"You are assessing {gene} as a drug target for {SUB_LABEL[subtype]} breast cancer. Give a target "
         f"score from 0 (not a target) to 10 (as strong as ERBB2 in HER2-enriched disease) and list the "
         f"claims your assessment rests on, each typed, with a direction and your confidence. ")
    if evidence is not None:
        p += ("You may use this evidence table, computed from TCGA, CPTAC, DepMap and Open Targets; it is "
              f"the only source you may rely on:\n{json.dumps(evidence, indent=0)}\n")
    else:
        p += "Use only what you already know; you have no tools and no documents.\n"
    return p + SCHEMA


def call(client, text):
    r = client.messages.create(model=MODEL, max_tokens=900, temperature=0,
                               messages=[{"role": "user", "content": text}])
    out = r.content[0].text
    start, end = out.find("{"), out.rfind("}")
    return json.loads(out[start:end + 1])


def evidence_table(cands: pd.DataFrame, decoys: pd.DataFrame) -> pd.DataFrame:
    """One frozen row per (gene, subtype) for both candidates and decoys."""
    scores = pd.read_csv(C.RESULTS / "scores" / "discovery.csv")
    gate = pd.read_csv(C.RESULTS / "scores" / "gate.csv", index_col=0)
    t1 = pd.read_csv(C.RESULTS / "truth" / "t1_gene_effect.csv", index_col=0)
    t2 = pd.read_csv(C.RESULTS / "truth" / "t2_clinical_precedence.csv").set_index("gene")
    gtex = pd.read_csv(C.RESULTS / "truth" / "gtex.csv", index_col=0)
    items = pd.concat([cands.assign(role="candidate"), decoys.assign(role="decoy")], ignore_index=True)
    rows = []
    for r in items.itertuples(index=False):
        sc = scores[(scores.gene == r.gene) & (scores.subtype == r.subtype)].iloc[0]
        rows.append({"gene": r.gene, "subtype": r.subtype, "role": r.role,
                     "rna_g": float(sc.rna_g), "protein_g": float(sc.protein_g),
                     "tumour_vs_normal_g": float(gate.loc[r.gene, "tumour_vs_normal_g"]),
                     "mean_effect_subtype": float(t1.loc[r.gene, f"mean_{r.subtype}"]),
                     "n_dependent_lines": int(t1.loc[r.gene, f"n_dependent_{r.subtype}"]),
                     "clinical_stage": int(t2.max_stage.get(r.gene, 0)),
                     "gtex_breast_tpm": float(gtex.gtex_breast_tpm.get(r.gene, np.nan)),
                     "gtex_max_other_tpm": float(gtex.gtex_max_other_tpm.get(r.gene, np.nan))})
    e = pd.DataFrame(rows)
    e["gtex_ratio"] = e.gtex_breast_tpm / e.gtex_max_other_tpm.replace(0, np.nan)
    return e


def check(claim, ev) -> str:
    col, thr = THRESH.get(claim["type"], (None, None))
    if col is None or pd.isna(ev.get(col)):
        return "unsupported"
    truth_yes = (ev[col] <= thr) if col == "mean_effect_subtype" else (ev[col] >= thr)
    said_yes = claim.get("direction") == "yes"
    return "supported" if truth_yes == said_yes else "contradicted"


def pubmed_counts(genes):
    strat = pd.read_csv(C.RESULTS / "ladder" / "strata.csv", index_col=0)
    return strat.pubmed.reindex(genes)


def run():
    t0 = time.time()
    import anthropic
    from scipy.stats import spearmanr
    from sklearn.metrics import roc_auc_score
    OUT.mkdir(parents=True, exist_ok=True)
    client = anthropic.Anthropic()

    top = pd.read_csv(C.RESULTS / "dossier" / "top10.csv")
    cands = top[["gene", "subtype", "arm"]]
    # decoys: same strata as the candidates, nominated by no arm, gate-passing
    strat = pd.read_csv(C.RESULTS / "ladder" / "strata.csv", index_col=0)
    lists = pd.read_csv(C.RESULTS / "nominate" / "lists_discovery.csv")
    gate = pd.read_csv(C.RESULTS / "scores" / "gate.csv", index_col=0)["gate"]
    nominated = set(lists.gene)
    rng = np.random.default_rng(7)
    drows = []
    for r in cands.itertuples(index=False):
        pool = strat[(strat.stratum == strat.loc[r.gene, "stratum"]) & ~strat.index.isin(nominated)
                     & strat.index.map(gate).fillna(False).astype(bool)].index
        pool = [g for g in pool if g not in {d["gene"] for d in drows}]
        drows.append({"gene": rng.choice(pool), "subtype": r.subtype, "arm": "decoy"})
    decoys = pd.DataFrame(drows)
    ev = evidence_table(cands, decoys)
    ev.to_csv(OUT / "items.csv", index=False)

    results = []
    for cond in ("closed", "open"):
        for r in ev.itertuples(index=False):
            evid = None if cond == "closed" else {k: (round(v, 3) if isinstance(v, float) else v)
                                                  for k, v in r._asdict().items() if k not in ("role",)}
            try:
                d = call(client, prompt(r.gene, r.subtype, evid))
            except Exception as e:                                   # one bad response must not stop the run
                results.append({"condition": cond, "gene": r.gene, "subtype": r.subtype, "role": r.role,
                                "error": str(e)[:200]})
                continue
            for cl in d.get("claims", []):
                results.append({"condition": cond, "gene": r.gene, "subtype": r.subtype, "role": r.role,
                                "target_score": d.get("target_score"), "claim_type": cl.get("type"),
                                "direction": cl.get("direction"), "confidence": cl.get("confidence"),
                                "text": cl.get("text"), "label": check(cl, r._asdict())})
            if not d.get("claims"):
                results.append({"condition": cond, "gene": r.gene, "subtype": r.subtype, "role": r.role,
                                "target_score": d.get("target_score"), "claim_type": None})
        print(f"  {cond}-book: {len(ev)} dossiers", flush=True)
    res = pd.DataFrame(results)
    res.to_csv(OUT / "claims.csv", index=False)

    summary = {"model": MODEL, "knowledge_cutoff": CUTOFF, "n_candidates": int(len(cands)),
               "n_decoys": int(len(decoys)), "conditions": {}}
    for cond, d in res.groupby("condition"):
        typed = d[d.label.isin(["supported", "contradicted"])]
        faith = typed.groupby("claim_type").label.apply(lambda s: float((s == "supported").mean())).to_dict()
        per_gene = (d.dropna(subset=["target_score"]).groupby(["gene", "role"]).target_score.first()
                    .reset_index())
        auc = float(roc_auc_score((per_gene.role == "candidate").astype(int), per_gene.target_score)) \
            if per_gene.role.nunique() == 2 else None
        rho = float(spearmanr(per_gene.target_score, pubmed_counts(per_gene.gene).values).statistic)
        summary["conditions"][cond] = {
            "faithfulness_overall": float((typed.label == "supported").mean()) if len(typed) else None,
            "faithfulness_by_claim_type": faith, "n_claims_checked": int(len(typed)),
            "n_claims_unsupported": int((d.label == "unsupported").sum()),
            "decoy_auc": auc, "pubmed_spearman": rho,
            "P7_holds": bool((typed.label == "supported").mean() < 0.8) if len(typed) else None,
            "P8_holds": bool(rho >= 0.5 and (auc is not None and auc < 0.7))}
        print(f"  {cond}: faithfulness {summary['conditions'][cond]['faithfulness_overall']:.2f}, "
              f"decoy AUC {auc}, PubMed rho {rho:.2f}")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    P.log_run("audit", {"model": MODEL, "seed": 7},
              [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*"))], time.time() - t0)
    print(f"done in {time.time() - t0:,.0f}s")


if __name__ == "__main__":
    sys.exit(run())
