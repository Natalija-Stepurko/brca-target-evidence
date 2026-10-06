"""Build the project page as one self-contained HTML file: docs/index.html (served by GitHub Pages).

One scrollable page: the question, why the study exists, the design, the pre-registered predictions,
the results (slots that fill from results/ once the stages have run), the agent audit, the exome
module, data and licences, literature, limits. Design facts (cohort sizes, K, thresholds) are stated
once in DESIGN and repeated here; the page carries no number that a stage produces until that stage
has written it to results/.

    python docs/site/build.py        -> docs/index.html
"""
import html
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "index.html"
RES = ROOT / "results"
REPO = "https://github.com/Natalija-Stepurko/brca-target-evidence"
DESIGN = f"{REPO}/blob/main/docs/DESIGN.md"
LIT = f"{REPO}/blob/main/research/literature.md"
SISTER = "https://natalija-stepurko.github.io/single-cell-fm-probing/"

INK, MUTED, RULE, PANEL = "#16191D", "#5B646E", "#DDE1E4", "#FFFFFF"
RNA, PROT, DNA, TRUTH = "#2D5BD1", "#C06014", "#6E7880", "#0E7C7B"
SANS = "ui-sans-serif,system-ui,-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif"
MONO = "ui-monospace,SFMono-Regular,Menlo,Consolas,monospace"


def esc(s):
    return html.escape(str(s), quote=False)


# ----------------------------------------------------------------------------- status
def stage_status():
    """Which stages have written results. Everything is pending until results/run_log.json exists."""
    log = RES / "run_log.json"
    if not log.exists():
        return {}
    return {r["stage"]: r for r in json.loads(log.read_text()).get("runs", [])}


STATUS = stage_status()
PRE_REGISTERED = "2026-10-06"


def pending(stage, what):
    """A result slot. Fills from results/ once `stage` has run; a dashed box until then."""
    if stage in STATUS:
        raise NotImplementedError(f"stage {stage} has results but the page has no renderer for it yet")
    return (f'<div class="pending"><span class="chip">pending</span> {esc(what)} — fills from '
            f'<code>results/</code> once stage <code>{esc(stage)}</code> has run.</div>')


# ----------------------------------------------------------------------------- figures
def fig_flow():
    """Schematic: layers -> arms -> lists -> held-out truth."""
    w, h = 980, 300
    o = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="Study flow: three omics layers from 77 '
         f'tumours feed six nomination arms; each produces a top-50 list per subtype, scored against '
         f'held-out dependency and clinical evidence that no arm sees.">']

    def box(x, y, bw, bh, label, sub="", fill=PANEL, stroke=RULE, col=INK, dash=""):
        o.append(f'<rect x="{x}" y="{y}" width="{bw}" height="{bh}" rx="4" fill="{fill}" stroke="{stroke}" '
                 f'stroke-width="1.2"{" stroke-dasharray=" + chr(34) + dash + chr(34) if dash else ""}/>')
        ty = y + bh / 2 + (4 if not sub else -2)
        o.append(f'<text x="{x + bw / 2}" y="{ty:.0f}" text-anchor="middle" font-size="13" '
                 f'font-weight="600" fill="{col}" font-family="{SANS}">{esc(label)}</text>')
        if sub:
            o.append(f'<text x="{x + bw / 2}" y="{ty + 16:.0f}" text-anchor="middle" font-size="11" '
                     f'fill="{MUTED}" font-family="{MONO}">{esc(sub)}</text>')

    def arrow(x1, y1, x2, y2):
        o.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{MUTED}" stroke-width="1.2" '
                 f'marker-end="url(#ah)"/>')

    o.append(f'<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
             f'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{MUTED}"/></marker></defs>')
    # column titles
    for x, t in ((20, "77 discovery tumours"), (300, "nomination arms"), (560, "top-50 per subtype"),
                 (790, "held-out truth")):
        o.append(f'<text x="{x}" y="22" font-size="11" fill="{MUTED}" font-family="{MONO}" '
                 f'letter-spacing=".06em">{esc(t.upper())}</text>')
    # layers
    box(20, 44, 200, 44, "RNA-seq", "TCGA", col=RNA)
    box(20, 104, 200, 44, "Proteome", "CPTAC, Mertins 2016", col=PROT)
    box(20, 164, 200, 44, "DNA", "mutations · copy number", col=DNA)
    box(20, 230, 200, 44, "RNA-seq, ~1,000 tumours", "all TCGA-BRCA", col=RNA, dash="4 3")
    # arms
    arms = ["R", "R+P", "R+D", "R+P+D", "MOFA+", "R-all"]
    for i, a in enumerate(arms):
        y = 44 + i * 40
        box(300, y, 160, 30, a, col=INK if a != "R" else TRUTH)
        for ly, want in ((66, "R" in a), (126, "P" in a or a == "MOFA+"), (186, "D" in a or a == "MOFA+")):
            if want and a != "R-all":
                arrow(220, ly, 300, y + 15)
        if a == "R-all":
            arrow(220, 252, 300, y + 15)
        arrow(460, y + 15, 560, y + 15)
        box(560, y, 150, 30, "basal · HER2 · luminal", col=MUTED)
        arrow(710, y + 15, 790, 150)
    # truth
    box(790, 90, 170, 50, "Dependency", "DepMap 24Q4, 53 lines", fill="#EAF4F4", stroke=TRUTH, col=TRUTH)
    box(790, 160, 170, 50, "Clinical precedence", "Open Targets 26.09", fill="#EAF4F4", stroke=TRUTH, col=TRUTH)
    o.append(f'<text x="875" y="250" text-anchor="middle" font-size="11" fill="{MUTED}" '
             f'font-family="{MONO}">never used in nomination</text>')
    o.append("</svg>")
    return "".join(o)


def fig_ladder():
    """Schematic of the ladder a list is read against."""
    rungs = [("Replication", "frozen pipeline on 122 Krug 2020 tumours", TRUTH),
             ("Positive controls", "ERBB2 in HER2-enriched, ESR1 in luminal", INK),
             ("Baseline arm", "RNA alone, same tumours, same universe, same K", RNA),
             ("Label-permutation null", "subtype labels shuffled, every arm re-run", MUTED),
             ("Matched-random floor", "1,000 lists matched on PubMed count, expression, protein detection", MUTED)]
    w, rh = 980, 46
    h = 30 + rh * len(rungs)
    o = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="The ladder every list is read against, '
         f'from the matched-random floor at the bottom to replication at the top.">']
    for i, (name, what, col) in enumerate(rungs):
        y = 20 + i * rh
        o.append(f'<line x1="40" y1="{y + 22}" x2="{w - 40}" y2="{y + 22}" stroke="{RULE}"/>')
        o.append(f'<circle cx="40" cy="{y + 22}" r="6" fill="{col}"/>')
        o.append(f'<text x="60" y="{y + 18}" font-size="13" font-weight="600" fill="{INK}" '
                 f'font-family="{SANS}">{esc(name)}</text>')
        o.append(f'<text x="60" y="{y + 36}" font-size="11.5" fill="{MUTED}" font-family="{MONO}">'
                 f'{esc(what)}</text>')
    o.append("</svg>")
    return "".join(o)


# ----------------------------------------------------------------------------- sections
SECTIONS = [("question", "Question"), ("why", "Why"), ("design", "Design"), ("predictions", "Predictions"),
            ("results", "Results"), ("audit", "Agent audit"), ("exome", "Exome"), ("data", "Data"),
            ("literature", "Literature"), ("limits", "Limits")]


def nav():
    links = "".join(f'<a href="#{i}">{esc(t)}</a>' for i, t in SECTIONS)
    return (f'<nav class="topnav" aria-label="Sections"><a class="brand" href="#top">Target evidence</a>'
            f'<div class="navlinks">{links}</div>'
            f'<div class="navext"><a href="{REPO}">Code</a><a href="{DESIGN}">Design</a></div></nav>')


def header():
    return f"""<header class="page" id="top">
  <p class="eyebrow">Breast cancer · multi-omics · target identification · study 2 of a series</p>
  <h1>Does the protein layer pick better drug targets?</h1>
  <p class="lede">Three proteogenomic studies have nominated breast-cancer drug targets by combining DNA,
  RNA and protein measurements from the same tumours. Whether the extra layers improve the shortlist has
  never been measured. This study measures it: the same tumours, the same gene universe, the same list
  size, with and without each layer, scored against dependency and clinical evidence that played no part
  in the nomination, against random genes matched for how well studied they are.</p>
  <div class="status"><span class="chip live">pre-registered {PRE_REGISTERED}</span>
  <span class="chip">no data downloaded yet</span><span class="chip">results pending</span></div>
  <p class="links"><a href="{DESIGN}">Design (pre-registered)</a> · <a href="{LIT}">Literature</a> ·
  <a href="{REPO}">Repository</a> · <a href="{SISTER}">Study 1: single-cell cell states and survival</a></p>
</header>"""


def s_question():
    return f"""<section class="sec" id="question">
  <h2>The question</h2>
  <p>In breast cancer, does adding protein (and DNA) evidence to RNA nominate drug targets that are more
  often confirmed by held-out evidence than RNA alone — beyond what random genes matched for popularity
  achieve — and do the nominations replicate in an independent cohort?</p>
  <p>And a second question: when an AI agent writes the rationale for a nominated target, how many of its
  claims survive a check against the held-out data, and does it prefer famous genes?</p>
  <p class="note">Expected answer, stated before any data: little or no list-level gain from the protein
  layer, because the genes that are druggable tend to be the ones whose mRNA and protein levels agree. A
  clean null is a measurement here: it tells a target-discovery team when proteomics is worth paying for.
  Either answer is reported.</p>
</section>"""


def s_why():
    return f"""<section class="sec" id="why">
  <h2>Why this study</h2>
  <div class="cards">
    <div class="card"><h3>Already done</h3><p>Nominating breast-cancer targets from integrated tumour
    omics: Mertins 2016 (77 TCGA tumours with proteomics), Krug 2020 (122 prospective tumours), Savage
    2024 (1,043 patients, 10 cancer types), Lin 2025 (DepMap-filtered breast targets), Pacini 2024 (370
    priority targets across 27 cancer types). A sixth list adds little.</p></div>
    <div class="card"><h3>Never measured</h3><p>Whether the extra layers help. No study compares
    integrated with RNA-only nomination at equal list size, scores both against held-out evidence, and
    controls for gene popularity. Savage 2024 never ran an RNA arm and excluded breast. Krug 2020 reports
    zero overlap of its drivers with Mertins 2016's; nobody has measured how reproducible nominations
    are.</p></div>
    <div class="card"><h3>What the evidence predicts</h3><p>For prediction tasks, RNA carries most of the
    signal and protein is roughly at parity (Rappoport 2018, Herrmann 2021, Gonçalves 2022, DrEval
    2026). Druggable genes have the higher mRNA–protein correlation (Savage 2024). None of that was
    tested on target nomination.</p></div>
  </div>
  <p>Study 1 of this series asked a related question of single-cell foundation models — do the cell
  states they find in breast tumours predict survival better than a linear baseline? — and read the
  answer against the same kind of ladder. <a href="{SISTER}">Its result was a clean null.</a></p>
</section>"""


def s_design():
    return f"""<section class="sec" id="design">
  <h2>Design</h2>
  <div class="figwrap">{fig_flow()}</div>
  <p class="figcap">Three layers from the 77 TCGA-BRCA tumours that CPTAC profiled (Mertins 2016) feed
  six nomination arms. Each arm ranks one shared gene universe and returns its top 50 per PAM50 subtype
  (basal-like, HER2-enriched, luminal). The lists are scored against evidence no arm sees.</p>
  <h3>The arms</h3>
  <div class="scroll"><table>
    <thead><tr><th>Arm</th><th>Evidence</th><th>Combination</th><th>Tumours</th><th>Role</th></tr></thead>
    <tbody>
    <tr><td><b>R</b></td><td>RNA</td><td>—</td><td>77</td><td>baseline for the information question</td></tr>
    <tr><td><b>R+P</b></td><td>RNA, protein</td><td>mean of rank percentiles</td><td>77</td><td>the primary comparison, P2</td></tr>
    <tr><td><b>R+D</b></td><td>RNA, DNA</td><td>mean of rank percentiles</td><td>77</td><td>is DNA the more useful second layer? P3</td></tr>
    <tr><td><b>R+P+D</b></td><td>all three</td><td>mean of rank percentiles</td><td>77</td><td>everything a programme would have</td></tr>
    <tr><td><b>MOFA+</b></td><td>RNA, protein, copy number</td><td>factor weights</td><td>77</td><td>does the integration method matter? P4</td></tr>
    <tr><td><b>R-all</b></td><td>RNA</td><td>—</td><td>~1,000</td><td>baseline for the practical question, P6</td></tr>
    </tbody></table></div>
  <p>Per layer, a gene's evidence is a standardised mean difference between the subtype's tumours and the
  cohort's others (for DNA: copy-number difference with a cis-expression check, or mutation recurrence).
  A gene enters any list only if it is up in tumour against adjacent normal tissue; the same gate, from
  the same RNA data, applies to every arm.</p>
  <h3>Held-out truth</h3>
  <p><b>Dependency (primary).</b> DepMap 24Q4 CRISPR gene effect in the 53 breast cell lines with a
  screen, each assigned to a subtype. The endpoint is the mean gene effect of a list in the matched lines
  (more negative means the list is enriched for genes those cells depend on); precision against a binary
  hit set is secondary. Pan-essential genes are excluded. <b>Clinical precedence (secondary).</b> Targets
  of drugs that reached any clinical phase for a breast indication, Open Targets 26.09, frozen.</p>
  <h3>The ladder</h3>
  <div class="figwrap">{fig_ladder()}</div>
  <p class="figcap">Every list is read against all five rungs. Arm differences carry a paired bootstrap
  interval over tumours; the smallest effect of interest is 0.05 in mean gene effect or 5 genes in 50.</p>
</section>"""


PREDICTIONS = [
    ("P1", "Every arm's top-50 beats its matched-random floor (percentile ≤ 5%) in at least two of three subtypes",
     "if not, the truth is uninformative at this sample size and P2–P6 are reported but not interpreted"),
    ("P2", "R+P does not beat R by the smallest effect of interest; the interval is reported either way",
     "primary prediction"),
    ("P3", "R+D adds no more over R than R+P does", ""),
    ("P4", "MOFA+ does not beat R+P+D rank stacking", ""),
    ("P5", "Each arm's discovery–replication overlap exceeds the permutation null; the protein arm replicates no worse than RNA", ""),
    ("P6", "RNA on ~1,000 tumours scores at least as well as RNA + protein on 77", "the practical question"),
    ("P7", "Fewer than 80% of the agent's checkable claims are supported by the held-out data", "agent audit"),
    ("P8", "The agent's scores track PubMed count (ρ ≥ 0.5) and separate candidates from matched decoys at AUC < 0.7", "agent audit"),
]


def s_predictions():
    rows = "".join(f'<tr><td><b>{p}</b></td><td>{esc(t)}</td><td class="muted">{esc(n)}</td>'
                   f'<td><span class="chip">pending</span></td></tr>' for p, t, n in PREDICTIONS)
    return f"""<section class="sec" id="predictions">
  <h2>Predictions, fixed before any data</h2>
  <div class="scroll"><table>
    <thead><tr><th></th><th>Prediction</th><th>Note</th><th>Outcome</th></tr></thead><tbody>{rows}</tbody>
  </table></div>
  <p class="note">Committed in <a href="{DESIGN}">DESIGN.md §8</a> on {PRE_REGISTERED}. Any later change to
  the design is logged in its §12 with the commit and the reason.</p>
</section>"""


def s_results():
    return f"""<section class="sec" id="results">
  <h2>Results</h2>
  {pending("ladder", "Each arm against its matched-random floor and the permutation null, per subtype (P1)")}
  {pending("ladder", "Paired differences between arms with bootstrap intervals (P2–P4, P6)")}
  {pending("replicate", "Discovery–replication overlap per arm against the permutation null (P5)")}
  {pending("dossier", "The candidate dossiers: the top genes of the best arm per subtype, each with every rung")}
</section>"""


def s_audit():
    return """<section class="sec" id="audit">
  <h2>The agent audit</h2>
  <p>Agentic systems now write target rationales. Published evaluations ask whether the system found the
  right gene or whether a wet-lab experiment worked; none has scored the rationale itself, claim by
  claim, against held-out evidence while controlling for how famous each gene is.</p>
  <ol>
    <li><b>Items.</b> The top-10 genes per subtype from the best arm (30 candidates) and 30 decoys drawn
    from the same popularity, expression and detection strata that no arm nominated. The agent is not
    told which is which.</li>
    <li><b>Agent.</b> One named model, version and knowledge cutoff stated, closed-book and open-book
    (read access to this repository's frozen tables, nothing else). It returns a target score and a list
    of typed claims: over-expressed in the subtype, up in tumour vs normal, dependency in matched lines,
    clinical precedence, restricted in normal tissue.</li>
    <li><b>Check.</b> Every typed claim is verified automatically against the frozen tables and labelled
    supported, contradicted or unsupported.</li>
    <li><b>Report.</b> Faithfulness per claim type and condition; whether the agent's scores separate
    candidates from decoys; how strongly its scores track PubMed count.</li>
  </ol>
  """ + pending("audit", "Faithfulness, decoy discrimination and popularity correlation (P7, P8)") + """
</section>"""


def s_exome():
    return """<section class="sec" id="exome">
  <h2>From raw reads to a dependency</h2>
  <p>The tumour mutation calls above arrive pre-processed. To show the step from sequencer output to
  variants, the SEQC2 reference pair HCC1395 / HCC1395BL — a triple-negative breast-cancer line and its
  matched normal, with a published high-confidence somatic truth set — is aligned and called from raw
  exome reads (BWA-MEM2, GATK Mutect2) and the calls are scored against that truth set. HCC1395 is also a
  CRISPR-screened DepMap line, so its called mutations are traced into its own dependencies and into the
  basal-like candidate list.</p>
  """ + pending("wes", "SNV precision and recall against the SEQC2 truth set; mutations traced into dependencies") + """
</section>"""


DATA = [
    ("TCGA-BRCA: mutations, copy number, RNA-seq, PAM50, survival", "UCSC Xena", "open; TCGA Research Network acknowledged"),
    ("CPTAC breast proteome, Mertins 2016 (77 tumours, TCGA barcodes)", "PDC000173 / cBioPortal", "CC BY"),
    ("CPTAC breast cohort, Krug 2020 (122 tumours)", "LinkedOmics / cBioPortal", "CC BY"),
    ("DepMap Public 24Q4 CRISPR gene effect, model annotations", "figshare 27993248", "CC BY 4.0 (the last release under it; pinned)"),
    ("Open Targets Platform 26.09 clinical datasets, tractability", "EBI FTP (parquet)", "CC0; ChEMBL-derived fields CC BY-SA 3.0"),
    ("GTEx v10 median gene TPM by tissue", "GTEx portal", "open"),
    ("NCBI gene2pubmed (publication counts)", "NCBI FTP", "open"),
    ("SEQC2 HCC1395 / HCC1395BL exomes and somatic truth set v1.2.1", "ENA PRJNA489865; NCBI FTP", "open"),
    ("50,002-cell breast tumour atlas (study 1)", "CELLxGENE, built in study 1", "CC BY 4.0"),
]


def s_data():
    rows = "".join(f"<tr><td>{esc(a)}</td><td>{esc(b)}</td><td>{esc(c)}</td></tr>" for a, b, c in DATA)
    return f"""<section class="sec" id="data">
  <h2>Data and licences</h2>
  <div class="scroll"><table><thead><tr><th>Dataset</th><th>Route</th><th>Licence</th></tr></thead>
  <tbody>{rows}</tbody></table></div>
  <p class="note">Every file is recorded with its URL, release and sha256 in <code>results/MANIFEST.sha256</code>
  when stage <code>data</code> runs. No data from Sanger Project Score or Cell Model Passports are used
  (non-commercial terms). Open Targets' association score is never an input: it contains dependency and
  literature evidence that would make the truth circular.</p>
</section>"""


REFS = [
    ("Mertins P et al.", "Proteogenomics connects somatic mutations to signalling in breast cancer", "Nature 2016", "10.1038/nature18003"),
    ("Krug K et al.", "Proteogenomic landscape of breast cancer tumorigenesis and targeted therapy", "Cell 2020", "10.1016/j.cell.2020.10.036"),
    ("Savage SR et al.", "Pan-cancer proteogenomics expands the landscape of therapeutic targets", "Cell 2024", "10.1016/j.cell.2024.05.039"),
    ("Lin HK, Dai J, Pusztai L", "Breast-cancer targets from CRISPR dependency screens, by subtype", "Breast Cancer Res Treat 2025", "10.1007/s10549-025-07817-0"),
    ("Pacini C et al.", "A comprehensive clinically informed map of dependencies in cancer cells and framework for target prioritisation", "Cancer Cell 2024", "10.1016/j.ccell.2023.12.016"),
    ("Rappoport N, Shamir R", "Multi-omic and multi-view clustering algorithms: review and cancer benchmark", "Nucleic Acids Res 2018", "10.1093/nar/gky889"),
    ("Herrmann M et al.", "Large-scale benchmark study of survival prediction methods using multi-omics data", "Brief Bioinform 2021", "10.1093/bib/bbaa167"),
    ("Gonçalves E et al.", "Pan-cancer proteomic map of 949 human cell lines", "Cancer Cell 2022", "10.1016/j.ccell.2022.06.010"),
    ("Minikel EV et al.", "Refining the impact of genetic evidence on clinical success", "Nature 2024", "10.1038/s41586-024-07316-0"),
    ("Haynes WA, Tomczak A, Khatri P", "Gene annotation bias impedes biomedical research", "Sci Rep 2018", "10.1038/s41598-018-19333-x"),
    ("Chang L et al.", "Systematic profiling of conditional pathway activation identifies context-dependent synthetic lethalities", "Cancer Cell 2021", ""),
    ("Fang LT et al.", "Establishing community reference samples, data and call sets for benchmarking cancer mutation detection (SEQC2)", "Nat Biotechnol 2021", "10.1038/s41587-021-00993-6"),
]


def s_literature():
    items = []
    for who, title, venue, doi in REFS:
        link = f' <a href="https://doi.org/{doi}">doi</a>' if doi else ""
        items.append(f"<li>{esc(who)}. <i>{esc(title)}</i>. {esc(venue)}.{link}</li>")
    return f"""<section class="sec" id="literature">
  <h2>Literature</h2>
  <p>The survey behind the design, with a strength tag on every entry, is in
  <a href="{LIT}">research/literature.md</a>. The references the design rests on most:</p>
  <ol class="refs">{"".join(items)}</ol>
</section>"""


def s_limits():
    return """<section class="sec" id="limits">
  <h2>Limits, stated in advance</h2>
  <ul>
    <li><b>Small cohorts.</b> 77 discovery tumours with protein, 122 for replication; 53 breast cell
    lines, of which only 6–11 are luminal. The luminal arm may be uninformative.</li>
    <li><b>Weak, partly confounded truth.</b> Dependency tracks expression (hence the matched floor);
    cell lines drift from tumours; ESR1- and ERBB2-type expression addictions do not reproduce in
    organoid screens and are flagged. Clinical precedence is small and popularity-biased.</li>
    <li><b>Asymmetric layers.</b> CPTAC breast has no normal-tissue proteomics, so the tumour-vs-normal
    gate is RNA-based in every arm.</li>
    <li><b>Two integration methods only.</b> Equal-weight rank stacking and MOFA+, neither tuned on the
    truth. Better integration methods remain untested here.</li>
    <li><b>Knowledge cutoff.</b> The held-out data pre-date current language models' training data, so
    the agent audit leans on decoys and claim-level checks, with a time-split secondary check.</li>
  </ul>
</section>"""


def footer():
    return f"""<footer><p>Built {date.today().isoformat()} from <code>docs/site/build.py</code>. Every number
  that a pipeline stage produces is read from <code>results/</code> when this page is built; none is typed
  in. <a href="{REPO}">{REPO.split("//")[1]}</a></p></footer>"""


# ----------------------------------------------------------------------------- page
CSS = f"""
:root{{--paper:#F7F8F9;--panel:{PANEL};--ink:{INK};--muted:{MUTED};--rule:{RULE};--rna:{RNA};--prot:{PROT};
  --truth:{TRUTH};--sans:{SANS};--mono:{MONO}}}
*{{box-sizing:border-box}}
html{{scroll-behavior:smooth;scroll-padding-top:58px}}
body{{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);line-height:1.6;
  -webkit-font-smoothing:antialiased;padding-inline:16px}}
.wrap{{max-width:1040px;margin:0 auto;padding:24px 0 80px}}
a{{color:var(--ink)}}
h1{{font-size:clamp(1.6rem,4vw,2.4rem);line-height:1.15;margin:.3em 0 .5em;letter-spacing:-.01em}}
h2{{font-size:1.45rem;margin:2.2em 0 .6em;letter-spacing:-.01em}}
h3{{font-size:1.05rem;margin:1.6em 0 .4em}}
.eyebrow{{font-family:var(--mono);font-size:.72rem;color:var(--muted);letter-spacing:.06em;text-transform:uppercase;margin:0}}
.lede{{font-size:1.08rem;max-width:62ch}}
.note{{color:var(--muted);font-size:.95rem}}
.muted{{color:var(--muted)}}
.status{{display:flex;flex-wrap:wrap;gap:8px;margin:12px 0}}
.chip{{display:inline-block;font-family:var(--mono);font-size:.7rem;padding:3px 9px;border-radius:999px;
  border:1px solid var(--rule);background:var(--panel);color:var(--muted);white-space:nowrap}}
.chip.live{{border-color:var(--truth);color:var(--truth)}}
.links{{font-size:.95rem}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px;margin:1em 0}}
.card{{background:var(--panel);border:1px solid var(--rule);border-radius:8px;padding:14px 16px}}
.card h3{{margin:0 0 .4em;font-size:.95rem}}
.card p{{margin:0;font-size:.95rem}}
.figwrap{{background:var(--panel);border:1px solid var(--rule);border-radius:8px;padding:12px;overflow-x:auto}}
.figwrap svg{{display:block;width:100%;min-width:640px;height:auto}}
.figcap{{color:var(--muted);font-size:.9rem;margin:.5em 0 1.4em}}
.scroll{{overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:.93rem;background:var(--panel)}}
th,td{{text-align:left;padding:8px 10px;border-bottom:1px solid var(--rule);vertical-align:top}}
th{{font-family:var(--mono);font-size:.72rem;letter-spacing:.05em;text-transform:uppercase;color:var(--muted)}}
.pending{{border:1px dashed var(--rule);border-radius:8px;padding:14px 16px;margin:12px 0;color:var(--muted);font-size:.95rem}}
.pending .chip{{margin-right:6px}}
code{{font-family:var(--mono);font-size:.88em;background:#EEF1F2;padding:1px 5px;border-radius:4px}}
ol.refs{{font-size:.92rem;padding-left:1.3em}}
ol.refs li{{margin:.35em 0}}
footer{{margin-top:60px;border-top:1px solid var(--rule);padding-top:14px;color:var(--muted);font-size:.85rem}}
/* sticky section navbar */
.topnav{{position:sticky;top:0;z-index:50;display:flex;align-items:center;gap:18px;margin:0 -16px;
  padding:0 max(16px,calc((100% - 1040px)/2 + 16px));height:48px;background:rgba(247,248,249,.94);
  backdrop-filter:blur(6px);border-bottom:1px solid var(--rule);font-family:var(--mono);font-size:.7rem}}
.topnav a{{color:var(--muted);text-decoration:none;white-space:nowrap}}
.topnav a:hover{{color:var(--truth)}}
.topnav .brand{{color:var(--ink);font-weight:600;flex:none}}
.navlinks{{display:flex;gap:16px;overflow-x:auto;scrollbar-width:none;flex:1 1 auto;min-width:0;padding:0 12px}}
.navlinks::-webkit-scrollbar{{display:none}}
.navlinks a{{padding:14px 0 12px;border-bottom:2px solid transparent}}
.navlinks a.on{{color:var(--ink);border-bottom-color:var(--truth)}}
.navext{{display:flex;gap:14px;flex:none}}
.navext a{{color:var(--truth);font-weight:600}}
@media(max-width:640px){{.topnav .brand{{display:none}}.topnav{{gap:10px}}}}
"""

JS = """<script>
(function () {
  const links = [...document.querySelectorAll('.topnav .navlinks a')];
  const secs = links.map(a => document.getElementById(a.getAttribute('href').slice(1))).filter(Boolean);
  let ticking = false;
  function update() {
    ticking = false;
    const y = window.scrollY + 70;
    let cur = secs[0];
    for (const s of secs) if (s.offsetTop <= y) cur = s;
    if (window.innerHeight + window.scrollY >= document.body.offsetHeight - 2) cur = secs[secs.length - 1];
    links.forEach(a => a.classList.toggle('on', a.getAttribute('href') === '#' + cur.id));
  }
  addEventListener('scroll', () => { if (!ticking) { ticking = true; requestAnimationFrame(update); } });
  update();
})();
</script>"""

FAVICON = ("<link rel=\"icon\" href=\"data:image/svg+xml,"
           "%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'%3E"
           "%3Crect x='4' y='18' width='6' height='10' fill='%232D5BD1'/%3E"
           "%3Crect x='13' y='11' width='6' height='17' fill='%23C06014'/%3E"
           "%3Crect x='22' y='4' width='6' height='24' fill='%230E7C7B'/%3E%3C/svg%3E\">")


def build():
    body = "\n\n".join([header(), s_question(), s_why(), s_design(), s_predictions(), s_results(),
                        s_audit(), s_exome(), s_data(), s_literature(), s_limits(), footer()])
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Does the protein layer pick better drug targets?</title>
<meta name="description" content="A pre-registered evaluation of multi-omics target nomination in breast
cancer: RNA alone against RNA plus protein and DNA, scored against held-out dependency and clinical
evidence with popularity-matched nulls, replicated in an independent cohort.">
{FAVICON}
<style>{CSS}</style>
</head>
<body>
{nav()}
<div class="wrap">
{body}
</div>
{JS}
</body>
</html>
"""
    for bad in ("rather than", "instead of", "None<"):
        assert bad not in page, f"banned text on the page: {bad!r}"
    OUT.write_text(page)
    print(f"wrote {OUT.relative_to(ROOT)} ({len(page):,} bytes); stages with results: {sorted(STATUS) or 'none'}")


if __name__ == "__main__":
    build()
