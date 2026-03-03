"""Build a Pyodide-targeted wheel from an existing Cantera Python sdist."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tarfile
import tempfile


def _find_latest_sdist(dist_dir: Path) -> Path:
    candidates = sorted(
        dist_dir.glob("cantera-*.tar.gz", case_sensitive=False),
        key=lambda candidate: candidate.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(f"No Cantera sdist found in '{dist_dir}'.")
    return candidates[0]


def _run(command: list[str], cwd: Path, env: dict[str, str] | None = None) -> None:
    print(f"+ {' '.join(command)}")
    subprocess.run(command, cwd=cwd, env=env, check=True)


def _extract_sdist(sdist: Path) -> Path:
    temp_root = Path(tempfile.mkdtemp(prefix="cantera-pyodide-src-"))
    with tarfile.open(sdist, mode="r:gz") as archive:
        archive.extractall(path=temp_root)
    candidates = [path for path in temp_root.iterdir() if path.is_dir()]
    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected exactly one source directory in extracted sdist at '{temp_root}'."
        )
    return candidates[0]


def main() -> int:
    repo_root = Path(subprocess.getoutput("git rev-parse --show-toplevel").strip())
    dist_dir = repo_root / "build" / "python_sdist" / "dist"
    sdist = _find_latest_sdist(dist_dir)
    if not sdist.is_file():
        raise FileNotFoundError(f"Cantera sdist not found: '{sdist}'.")

    outdir = repo_root / "build" / "pyodide_dist"
    outdir.mkdir(parents=True, exist_ok=True)

    env = dict(os.environ)
    env["CANTERA_PYODIDE"] = "ON"
    cmake_args = env.get("CMAKE_ARGS", "").strip()
    pyodide_cmake_args = "-DCMAKE_CXX_SCAN_FOR_MODULES=OFF"
    env["CMAKE_ARGS"] = f"{cmake_args} {pyodide_cmake_args}".strip()
    if "Boost_INCLUDE_DIRS" not in env and (conda_prefix := env.get("CONDA_PREFIX")):
        boost_include = Path(conda_prefix) / "include"
        if boost_include.is_dir():
            env["Boost_INCLUDE_DIRS"] = str(boost_include)

    print(f"Using sdist: {sdist}")
    print(f"Output dir: {outdir}")
    source_dir = _extract_sdist(sdist)
    print(f"Extracted source: {source_dir}")
    _run(["pyodide", "build", str(source_dir), "--outdir", str(outdir)],
         cwd=repo_root, env=env)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
