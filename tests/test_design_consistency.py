"""The numbers in src/bte/config.py must be the ones docs/DESIGN.md pre-registers."""
import re
from pathlib import Path

from bte import config as C
from bte.cli import STAGES

DESIGN = (Path(__file__).resolve().parents[1] / "docs" / "DESIGN.md").read_text()


def test_eight_predictions_are_pre_registered():
    assert re.findall(r"^\| \*\*P(\d)\*\*", DESIGN, flags=re.M) == [str(i) for i in range(1, 9)]


def test_list_size_and_sensitivity():
    assert f"K = {C.K}" in DESIGN
    assert f"K = {C.K_SENSITIVITY[0]} and {C.K_SENSITIVITY[1]} as sensitivity" in DESIGN


def test_thresholds_match():
    assert f"≤ {C.HIT_GENE_EFFECT}".replace("-", "−") in DESIGN   # typographic minus in the text
    assert f"{C.HIT_SELECTIVITY} more negative" in DESIGN
    assert f"Spearman ≥ {C.CIS_SPEARMAN_MIN}" in DESIGN
    assert f"*g* ≥ {C.GATE_HEDGES_G}" in DESIGN
    assert f"**{C.SEOI_GENE_EFFECT}** in mean Chronos" in DESIGN
    assert f"**{C.SEOI_PRECISION} genes**" in DESIGN


def test_releases_pinned():
    assert f"Public **{C.DEPMAP_RELEASE}**" in DESIGN
    assert f"Platform **{C.OPEN_TARGETS_RELEASE}**" in DESIGN
    assert C.SEQC2_PROJECT in DESIGN and C.HCC1395_DEPMAP_ID in DESIGN


def test_resample_counts():
    assert f"{C.N_MATCHED_RANDOM:,} lists" in DESIGN
    assert f"{C.N_PERMUTATIONS:,} permutations" in DESIGN
    assert f"{C.N_BOOTSTRAP:,} resamples" in DESIGN


def test_arms_and_stages_named_in_design():
    for arm in C.ARMS:
        assert f"**{arm}**" in DESIGN
    for stage, _ in STAGES:
        if stage not in ("scores", "nominate", "report"):
            assert f"`{stage}`" in DESIGN or stage in DESIGN.lower()


def test_no_banned_phrases():
    for path in (Path("docs/DESIGN.md"), Path("README.md"), Path("research/literature.md")):
        p = Path(__file__).resolve().parents[1] / path
        if p.exists():
            assert not re.search(r"rather than|instead of", p.read_text(), flags=re.I), path
