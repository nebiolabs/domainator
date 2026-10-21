"""Time a few common Domainator steps, to catch performance regressions.

    python scripts/bench/bench.py                       # print the timings
    python scripts/bench/bench.py --json_out base.json  # save them
    python scripts/bench/bench.py --baseline base.json  # compare with saved timings

Each step runs --repeat times and reports the fastest run, on the fixtures in test/data
(Staph_phages.gb is one 1.9 Mb genome with 1,865 CDSs). Run it on a quiet machine;
differences under about 10% are usually noise.
"""
import argparse
import io
import json
import os
import sys
import tempfile
import time
from pathlib import Path

DATA = Path(__file__).resolve().parents[2] / "test" / "data"
GENOMES = str(DATA / "Staph_phages.gb")
HMMS = [str(DATA / "pdonr_hmms.hmm"), str(DATA / "SPR.hmm")]


def _parse(parser):
    from domainator.seqio import parse_seqfiles
    def run(_):
        for _record in parse_seqfiles([GENOMES], genbank_parser=parser):
            pass
    return run


def _write(_):
    from domainator.seqio import parse_seqfiles, write_genbank
    records = list(parse_seqfiles([GENOMES]))
    start = time.perf_counter()
    write_genbank(records, io.StringIO())
    return time.perf_counter() - start


def _domainate(cpu):
    def run(tmp):
        from domainator.domainate import main
        main(["-i", GENOMES, "-r", *HMMS, "-o", os.path.join(tmp, "domainate.gb"), "--cpu", str(cpu), "-Z", "1000"])
    return run


def _domain_search(cpu):
    def run(tmp):
        from domainator.domain_search import main
        main(["-i", GENOMES, "-r", str(DATA / "SPR.hmm"), "-o", os.path.join(tmp, "domain_search.gb"),
              "--cpu", str(cpu), "-Z", "1000", "--cds_range", "5"])
    return run


def _extract_domains(tmp):
    from domainator.domainate import main as domainate
    from domainator.extract_domains import main as extract_domains
    annotated = os.path.join(tmp, "annotated.gb")
    if not os.path.exists(annotated):
        domainate(["-i", GENOMES, "-r", *HMMS, "-o", annotated, "--cpu", "0", "-Z", "1000"])
    start = time.perf_counter()
    extract_domains(["-i", annotated, "-o", os.path.join(tmp, "domains.gb")])
    return time.perf_counter() - start


def steps(cpu):
    return {
        "parse (biopython)": _parse("biopython"),
        "parse (lean)": _parse("lean"),
        "write_genbank": _write,
        "domainate": _domainate(cpu),
        "domain_search": _domain_search(cpu),
        "extract_domains": _extract_domains,
    }


def time_step(function, tmp, repeat):
    best = float("inf")
    for _ in range(repeat):
        start = time.perf_counter()
        measured = function(tmp)  # a step that needs setup times itself and returns its own time
        best = min(best, measured if measured is not None else time.perf_counter() - start)
    return best


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repeat", type=int, default=3, help="Runs per step; the fastest is reported.")
    parser.add_argument("--cpu", type=int, default=1, help="--cpu for the search steps (1 makes them comparable across machines).")
    parser.add_argument("--only", nargs="+", default=None, help="Run only the steps whose names contain one of these.")
    parser.add_argument("--json_out", default=None, help="Save the timings to this file.")
    parser.add_argument("--baseline", default=None, help="Compare with timings saved by --json_out.")
    params = parser.parse_args(argv)

    import warnings
    warnings.simplefilter("ignore")
    baseline = json.loads(Path(params.baseline).read_text()) if params.baseline else {}
    results = {}
    with tempfile.TemporaryDirectory() as tmp:
        for name, function in steps(params.cpu).items():
            if params.only and not any(word in name for word in params.only):
                continue
            seconds = time_step(function, tmp, params.repeat)
            results[name] = seconds
            line = f"{name:20s} {seconds:8.3f} s"
            if name in baseline:
                line += f"   {seconds / baseline[name]:5.2f}x baseline ({baseline[name]:.3f} s)"
            print(line, flush=True)
    if params.json_out:
        Path(params.json_out).write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    sys.exit(main())
