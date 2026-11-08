"""Reproduce the finite checks in On Fourier and Hadamard matrices of prime order."""
import argparse
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
GROUPS = ("w19", "order11", "minors", "orbits", "corners")


def python(script, *args):
    subprocess.run([sys.executable, "-B", script, *args], check=True)


def compile_cpp(source, work):
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    binary = work / source.stem
    subprocess.run([*compiler, "-std=c++17", "-O3", source, "-o", binary], check=True)
    return binary


def read_json(path):
    # Wall-clock time is the only field that varies between machines.
    return json.loads(path.read_text(), object_hook=lambda obj: {
        key: value for key, value in obj.items() if key != "elapsed_seconds"
    })


def order11(work):
    original = ROOT / "order11"
    branch, square = work / "a_branch", work / "a_square"
    for directory in (branch, square):
        directory.mkdir()
        for source in (original / directory.name).iterdir():
            if source.suffix in (".py", ".cpp"):
                shutil.copyfile(source, directory / source.name)
    for script in ("compute.py", "enumerate_labels.py", "verify.py"):
        python(branch / script)
    brute = compile_cpp(branch / "brute_labels.cpp", work)
    with (branch / "brute_input.txt").open() as inputs:
        report = subprocess.check_output([brute], stdin=inputs, text=True)
    print(report, end="")
    records = re.findall(r": (\d+) normalized labels, (\d+) assignments; expected (\d+)", report)
    expected = read_json(original / branch.name / "independent_brute_counts.json")["counts"]
    assert len(records) == len(expected) == 6
    for (labels, assignments, agreed), ref in zip(records, expected):
        assert int(labels) == int(agreed) == ref["normalized_labels"]
        assert int(assignments) == ref["full_assignments"] and ref["agrees_with_python"]

    python(square / "generate_rows.py")
    clique = compile_cpp(square / "complete.cpp", work)
    subprocess.run([clique, square / "rows.txt", square / "completions.json"], check=True)
    python(square / "verify.py")
    for directory in (branch, square):
        for reference in (original / directory.name).iterdir():
            output = directory / reference.name
            if reference.name == "independent_brute_counts.json":
                continue  # Independently checked from the C++ output above.
            if reference.suffix == ".json":
                assert read_json(output) == read_json(reference), reference.name
            if reference.suffix == ".txt":
                assert output.read_bytes() == reference.read_bytes(), reference.name
    print("Order 11: all archived counts, witnesses, and exclusion certificates agree.")


def orbits(work):
    src = ROOT / "orbits"
    sums = compile_cpp(src / "pair_sums.cpp", work)
    diffs = compile_cpp(src / "differences.cpp", work)
    python(src / "small.py", sums, diffs)
    for name in ("odd.py", "even.py", "parity.py", "prime109.py"):
        python(src / name)


def corners(work):
    src = ROOT / "corners"
    program = compile_cpp(src / "deletions.cpp", work)
    report = subprocess.check_output([program], text=True)
    # An extra nonlinear leaf is printed separately by the enumerator.
    rows = [tuple(map(int, line.split())) for line in report.splitlines()]
    primes = [19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103]
    assert [p for p, _, _ in rows] == primes
    assert rows[-1] == (103, 38060, 1)
    assert all(count == 1 for _, _, count in rows)
    for low, high, expected in [(19,31,737),(37,59,9256),(61,79,28336),(83,103,125984)]:
        nodes = sum(n for p,n,_ in rows if low <= p <= high)
        assert nodes == expected
        print(f"p = {low}..{high}: {nodes} visited assignments")
    for name in ("anchors.py", "strips.py", "inverse.py"):
        python(src / name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("group", nargs="?", default="all", choices=("all", *GROUPS))
    group = parser.parse_args().group
    if not __debug__ or os.environ.get("PYTHONOPTIMIZE"):
        parser.error("Run with Python assertions enabled (unset PYTHONOPTIMIZE).")
    for name in GROUPS if group == "all" else (group,):
        print(f"\n=== {name} ===", flush=True)
        with tempfile.TemporaryDirectory(prefix="code003-") as temporary:
            work = Path(temporary)
            if name == "w19":
                python(ROOT / name / "completions.py")
            elif name == "minors":
                python(ROOT / name / "check.py")
            else:
                {"order11": order11, "orbits": orbits, "corners": corners}[name](work)
        print(f"{name}: PASS", flush=True)


if __name__ == "__main__":
    main()
