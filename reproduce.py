"""Build the generator, regenerate the proof files, and verify them in order."""
import argparse
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from time import monotonic_ns


ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--through", type=int, choices=(10, 11), default=11)
    parser.add_argument("--out", type=Path, default=ROOT / "regenerated")
    parser.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    args = parser.parse_args()
    dest = args.out.resolve()
    dest.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="rg-build-", dir=dest) as directory:
        temp = Path(directory)
        binary = temp / ("search.exe" if os.name == "nt" else "search")
        subprocess.run([args.cxx, "-std=c++17", "-O3", str(ROOT / "search.cpp"),
                        "-o", str(binary)], check=True)
        for n in range(5, args.through + 1):
            raw = temp / f"n{n}.rgc"
            start = monotonic_ns()
            subprocess.run([str(binary), str(n), str(raw)], check=True)
            milliseconds = (monotonic_ns() - start) // 1_000_000
            if n == 10 and milliseconds > 90_000:
                print(f"10-mark search took {milliseconds} ms. The 11-mark tree may be much larger.",
                      flush=True)
            packed = temp / (raw.name + ".gz")
            with raw.open("rb") as inp, packed.open("wb") as sink:
                with gzip.GzipFile(filename="", mode="wb", fileobj=sink,
                                   mtime=0, compresslevel=6) as archive:
                    shutil.copyfileobj(inp, archive)
            raw.unlink()
            if n >= 10:
                subprocess.run([sys.executable, "-I", str(ROOT / "check.py"),
                                "--through", str(n), "--dir", str(temp)], check=True)
        for n in range(5, args.through + 1):
            name = f"n{n}.rgc.gz"
            (temp / name).replace(dest / name)


if __name__ == "__main__":
    main()
