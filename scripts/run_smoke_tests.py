"""Run data-independent numerical checks; never launch measured experiments."""
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
TESTS = (
    "test_pearson_scope.py",
    "test_signed_waveform.py",
    "test_tidal_waveform.py",
)


def main():
    # Runners expect this local, git-ignored directory on fresh clones.
    (ROOT / "results").mkdir(exist_ok=True)
    for name in TESTS:
        print(f"Running {name}", flush=True)
        subprocess.run(
            [sys.executable, str(ROOT / "tests" / name)],
            cwd=ROOT,
            check=True,
        )
    print(f"PASS: {len(TESTS)} smoke scripts (numerical checks only)", flush=True)


if __name__ == "__main__":
    main()
