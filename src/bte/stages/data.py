"""Stage `data`: download every input with sha256 provenance, then build the shared gene universe.

Downloads are resumable and skipped when the file is already present with the recorded sha256. The
universe (DESIGN §4) is protein-coding genes quantified in the discovery proteome in at least half of
the tumours, quantified in TCGA RNA-seq, and present in the DepMap gene-effect matrix.
"""
import gzip
import sys
import time
from pathlib import Path

import pandas as pd
import requests

from bte import config as C
from bte import provenance as P
from bte.sources import SOURCES

RAW = C.DATA / "raw"
OUT = C.RESULTS / "data"
MIN_PROTEIN_DETECTION = 0.5


# ----------------------------------------------------------------------------- download
def fetch(url: str, dest: Path, size_mb: float, retries=4):
    """Stream to dest.part with resume, then rename. Prints progress every ~50 MB."""
    part = dest.with_suffix(dest.suffix + ".part")
    have = part.stat().st_size if part.exists() else 0
    for attempt in range(retries):
        try:
            headers = {"Range": f"bytes={have}-"} if have else {}
            with requests.get(url, stream=True, timeout=120, headers=headers, allow_redirects=True) as r:
                if r.status_code == 416:            # already complete
                    break
                r.raise_for_status()
                mode = "ab" if r.status_code == 206 else "wb"
                if mode == "wb": have = 0
                with part.open(mode) as f:
                    last = have
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk); have += len(chunk)
                        if have - last > 50 << 20:
                            last = have
                            print(f"    {dest.name}: {have / 1e6:,.0f} / ~{size_mb:,.0f} MB", flush=True)
            break
        except (requests.RequestException, OSError) as e:
            if attempt == retries - 1: raise
            print(f"    retry {attempt + 1} after {type(e).__name__}", flush=True)
            time.sleep(10 * (attempt + 1))
    part.rename(dest)


def download_all():
    RAW.mkdir(parents=True, exist_ok=True)
    recorded = {}
    if P.MANIFEST.exists():
        for line in P.MANIFEST.read_text().splitlines():
            digest, rel, *_ = line.split("  ")
            recorded[rel] = digest
    for s in SOURCES:
        dest = RAW / s.local
        rel = str(dest.relative_to(C.ROOT))
        if dest.exists() and recorded.get(rel) == P.sha256(dest):
            print(f"  ok      {s.local}"); continue
        print(f"  fetch   {s.local}  ({s.release}, ~{s.size_mb:,.0f} MB)", flush=True)
        fetch(s.url, dest, s.size_mb)
        digest = P.record_file(dest, s.url, s.release)
        print(f"  sha256  {digest[:16]}…  {s.local}")


# ----------------------------------------------------------------------------- universe
def read_cct(path: Path) -> pd.DataFrame:
    """LinkedOmics .cct: tab-separated, genes as rows, samples as columns."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return pd.read_csv(f, sep="\t", index_col=0)


def depmap_genes(path: Path) -> pd.Series:
    """DepMap columns are 'SYMBOL (EntrezID)'; return symbol -> entrez."""
    cols = pd.read_csv(path, nrows=0).columns[1:]
    sym = cols.str.extract(r"^(\S+) \((\d+)\)$")
    return pd.Series(sym[1].astype(int).values, index=sym[0].values, name="entrez")


def hgnc_protein_coding(path: Path) -> pd.DataFrame:
    h = pd.read_csv(path, sep="\t", low_memory=False, usecols=["symbol", "entrez_id", "locus_group",
                                                               "prev_symbol", "alias_symbol"])
    return h[h.locus_group == "protein-coding gene"].dropna(subset=["entrez_id"]).astype({"entrez_id": int})


def discovery_samples() -> pd.DataFrame:
    """Mertins 2016 Supplementary Table 1: TCGA id, PAM50, receptor status and CPTAC's QC verdict."""
    import zipfile
    z = zipfile.ZipFile(RAW / "mertins2016_supplementary_tables.zip")
    name = "nature18003-s2/CPTAC_BC_SupplementaryTable01.xlsx"
    xlsx = RAW / "mertins2016_supplementary_table01.xlsx"
    if not xlsx.exists():
        xlsx.write_bytes(z.read(name))
    d = pd.read_excel(xlsx, sheet_name=0)
    d = d.rename(columns={"TCGA ID": "tcga_id", "PAM50": "pam50", "QC Status": "qc", "ER Status": "er",
                          "PR Status": "pr", "HER2 Status": "her2"})
    return d[["tcga_id", "pam50", "er", "pr", "her2", "qc"]]


def build_universe():
    OUT.mkdir(parents=True, exist_ok=True)
    hg = hgnc_protein_coding(RAW / "hgnc_complete_set.txt")
    samples = discovery_samples()
    passed = samples[samples.qc == "pass"]
    samples.to_csv(OUT / "discovery_samples.csv", index=False)
    prot = read_cct(RAW / "cptac2016_proteome_gene_itraq_logratio.cct.gz")
    prot.columns = prot.columns.str.replace(".", "-", regex=False)      # TCGA.A2.A0D0 -> TCGA-A2-A0D0
    missing = set(passed.tcga_id) - set(prot.columns)
    assert not missing, f"QC-pass tumours absent from the proteome file: {sorted(missing)}"
    prot = prot[list(passed.tcga_id)]                                   # the 77 that passed CPTAC QC
    rna_cols = pd.read_csv(RAW / "tcga_brca_HiSeqV2.gz", sep="\t", nrows=0).columns
    rna_genes = pd.read_csv(RAW / "tcga_brca_HiSeqV2.gz", sep="\t", usecols=[0]).iloc[:, 0]
    dm = depmap_genes(RAW / "depmap_24Q4_CRISPRGeneEffect.csv")

    tumour_cols = list(prot.columns)
    detection = prot.notna().mean(axis=1)
    in_prot = set(detection[detection >= MIN_PROTEIN_DETECTION].index)
    universe = hg[hg.symbol.isin(in_prot) & hg.symbol.isin(set(rna_genes)) & hg.symbol.isin(set(dm.index))]
    universe = universe.assign(protein_detection=universe.symbol.map(detection)).sort_values("symbol")
    universe[["symbol", "entrez_id", "protein_detection"]].to_csv(OUT / "universe.csv", index=False)

    counts = {
        "hgnc_protein_coding": int(len(hg)),
        "proteome_samples_profiled": int(len(samples)), "proteome_samples_qc_pass": int(len(tumour_cols)),
        "proteome_genes": int(len(prot)),
        "proteome_genes_detected_50pct": int(len(in_prot)),
        "tcga_rna_genes": int(len(rna_genes)), "tcga_rna_samples": int(len(rna_cols) - 1),
        "depmap_genes": int(len(dm)),
        "universe": int(len(universe)),
    }
    pd.Series(counts, name="n").to_csv(OUT / "counts.csv")
    for k, v in counts.items():
        print(f"  {k:<32} {v:>8,}")
    return counts


def run():
    t0 = time.time()
    print("downloading", flush=True)
    download_all()
    print("building the gene universe", flush=True)
    counts = build_universe()
    P.log_run("data", {"min_protein_detection": MIN_PROTEIN_DETECTION, "n_sources": len(SOURCES)},
              ["results/data/universe.csv", "results/data/counts.csv", "results/data/discovery_samples.csv",
               "results/MANIFEST.sha256"],
              time.time() - t0)
    print(f"done in {time.time() - t0:,.0f}s; universe = {counts['universe']:,} genes")


if __name__ == "__main__":
    sys.exit(run())
