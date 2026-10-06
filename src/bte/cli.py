"""`bte run <stage>`: the pipeline, one stage per module. Stages arrive one pull request at a time;
until a stage exists, asking for it says so and exits non-zero."""
import argparse
import sys

STAGES = [
    ("data", "download every input with sha256 provenance; build the gene universe"),
    ("scores", "per-layer evidence per gene and subtype; the tumour-vs-normal gate"),
    ("nominate", "the six arms, top-K per subtype; positive controls"),
    ("truth", "subtype-matched DepMap dependency and Open Targets clinical precedence"),
    ("ladder", "matched-random floor, label-permutation null, paired bootstrap; P1-P4, P6"),
    ("replicate", "the frozen pipeline on the Krug 2020 cohort; P5"),
    ("dossier", "candidate dossiers with every rung"),
    ("wes", "SEQC2 exome pair from raw reads to scored somatic calls"),
    ("audit", "agent-written dossiers for candidates and decoys, checked claim by claim; P7, P8"),
    ("report", "figures and tables for the project page"),
]


def main(argv=None):
    ap = argparse.ArgumentParser(prog="bte", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    run = sub.add_parser("run", help="run one stage")
    run.add_argument("stage", choices=[s for s, _ in STAGES])
    sub.add_parser("stages", help="list the stages")
    a = ap.parse_args(argv)
    if a.cmd == "stages":
        for s, d in STAGES:
            print(f"{s:<10} {d}")
        return 0
    print(f"stage {a.stage!r} is not implemented yet; see the open pull requests", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
