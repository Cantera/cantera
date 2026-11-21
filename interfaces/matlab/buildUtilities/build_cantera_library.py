#!/usr/bin/env python3
"""Build a self-contained Cantera shared library for the MATLAB toolbox.

This script produces a library by building Cantera from source with all
third-party dependencies vendored and statically linked, and verifies the result
against a per-platform allowlist of system libraries (see audit_library.py).

It is deliberately kept out of MATLAB: MATLAB injects its own library paths into
child processes, which breaks compilers and linkers.

Usage:
    python build_cantera_library.py --cantera-root <src> --prefix <dir>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

STAMP_NAME = ".build-stamp"

# Bump when the recipe changes in a way the SCons flag list does not capture.
# 2: normalize_layout also stages the Windows import library into lib/.
RECIPE_VERSION = 2

# Submodules the vendored build compiles directly. HighFive is excluded because
# hdf_support=n means it is never used.
REQUIRED_SUBMODULES = ("ext/fmt", "ext/sundials", "ext/eigen", "ext/yaml-cpp")


def check_submodules(cantera_root: Path) -> list[str]:
    """Check the vendored dependencies are at the commits Cantera records.

    A stale submodule looks healthy but fails later with missing headers.
    """
    try:
        out = subprocess.run(
            ["git", "-C", str(cantera_root), "submodule", "status"],
            capture_output=True, text=True, check=True,
        ).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        # Not a git checkout; fall back to checking the directories exist.
        empty = [p for p in REQUIRED_SUBMODULES
                 if not any((cantera_root / p).glob("*"))]
        if empty:
            return [f"Vendored dependencies are missing: {', '.join(empty)}"]
        return []

    # "<marker><sha> <path> (<describe>)"; '-' uninitialized, '+' wrong commit.
    stale = []
    for line in out.splitlines():
        if not line.strip():
            continue
        parts = line[1:].split()
        if len(parts) < 2:
            continue
        path = parts[1].replace("\\", "/")
        if path not in REQUIRED_SUBMODULES:
            continue
        if line[0] == "-":
            stale.append(f"{path}: not initialized")
        elif line[0] == "+":
            described = parts[2].strip("()") if len(parts) > 2 else "unknown"
            stale.append(f"{path}: at {described}, not the recorded commit")

    if stale:
        return ["Cantera submodules are not at the recorded commits:\n"
                + "\n".join(f"      {s}" for s in stale)
                + f"\n    git -C {cantera_root} submodule update --init "
                  f"--recursive"]

    return []


def check_prerequisites(boost_inc_dir: str | None = None) -> list[str]:
    """Report missing build tools before starting a long build."""
    problems = []

    # cantera_clib is generated at build time by sourcegen, which reads a
    # Doxygen tag file, so doxygen is a build requirement and not a docs-only
    # one. Cantera's Doxyfile sets CITE_BIB_FILES, so doxygen also needs perl.
    if shutil.which("doxygen") is None:
        problems.append(
            "doxygen is not on PATH; it generates the tag file sourcegen "
            "reads to produce cantera_clib.\n"
            "    conda install -c conda-forge doxygen"
        )

    if shutil.which("perl") is None:
        problems.append(
            "perl is not on PATH; doxygen needs it to run bibtex.\n"
            "    Windows: run from Git Bash, or install Strawberry Perl."
        )

    if boost_inc_dir:
        boost = Path(boost_inc_dir)
        if not boost.is_dir():
            problems.append(f"BOOST_INC_DIR does not exist: {boost}")
        elif not (boost / "boost" / "version.hpp").is_file():
            problems.append(
                f"BOOST_INC_DIR has no boost/version.hpp: {boost}\n"
                f"    It must contain boost/, not be boost/ itself."
            )

    # packaging is imported by SConstruct itself; the other two by sourcegen.
    missing = []
    for module in ("packaging", "jinja2", "ruamel.yaml"):
        try:
            __import__(module)
        except ImportError:
            missing.append(module)

    if missing:
        problems.append(
            f"The Cantera build requires {', '.join(missing)}, not importable "
            f"from {sys.executable}.\n"
            f"    python -m pip install {' '.join(missing)}"
        )

    return problems


def warn_about_inherited_options(cantera_root: Path, flags: list[str]) -> None:
    """Warn about cantera.conf options this script does not set.

    SCons rewrites cantera.conf every run, so only options absent from `flags`
    are actually inherited.
    """
    conf = cantera_root / "cantera.conf"
    if not conf.is_file():
        return

    ours = {flag.split("=", 1)[0] for flag in flags}
    try:
        theirs = set(re.findall(r"^(\w+)\s*=", conf.read_text(), re.M))
    except OSError:
        return

    inherited = sorted(theirs - ours)
    if inherited:
        print(f"warning: {conf} sets options this script does not control, "
              f"which will be\n         inherited by the build: "
              f"{', '.join(inherited)}.\n"
              f"         Rename the file for a build that matches CI.",
              file=sys.stderr)


def boost_flags(boost_inc_dir: str | None) -> list[str]:
    """Point SCons at Boost headers when not on the default include path.

    Boost (>= 1.83) is the one dependency Cantera does not vendor. It is
    header-only here, so it adds no runtime dependency.
    """
    if not boost_inc_dir:
        return []
    return [f"boost_inc_dir={Path(boost_inc_dir).as_posix()}"]


def scons_flags(prefix: Path) -> list[str]:
    """SCons options shared by every platform.

    Each ``system_*=n`` statically links the copy vendored under ``ext/``,
    leaving no runtime dependency behind.
    """
    return [
        f"prefix={prefix.as_posix()}",
        "layout=compact",
        "libdirname=lib",
        "system_eigen=n",
        "system_fmt=n",
        "system_yamlcpp=n",
        "system_sundials=n",
        "system_highfive=n",
        "system_blas_lapack=n",
        "hdf_support=n",
        "python_package=n",
        "f90_interface=n",
        "googletest=none",
        "package_build=y",
        "optimize=y",
        "debug=n",
        "renamed_shared_libraries=y",
        "versioned_shared_library=y",
        "use_rpath_linkage=n",
    ]


def platform_flags() -> list[str]:
    """Per-platform flags that make the library relocatable.

    ``no_debug_linker_flags`` reaches LINKFLAGS because every build here is
    ``debug=n``.
    """
    system = platform.system()

    if system == "Linux":
        link = [
            # MATLAB ships its own, usually older, libstdc++.
            "-static-libstdc++",
            "-static-libgcc",
            # Hide symbols from the vendored static archives so they cannot
            # collide with copies already loaded into MATLAB. Do NOT add
            # -fvisibility=hidden: it would hide Cantera's own exports too.
            "-Wl,--exclude-libs,ALL",
            "-Wl,-z,origin",
            # SCons substitutes '$$' to a literal '$'.
            "-Wl,-rpath,$$ORIGIN",
        ]
        return [f"no_debug_linker_flags={' '.join(link)}"]

    if system == "Darwin":
        # libc++ is part of the OS, so there is no C++ runtime to bundle.
        return ["no_debug_linker_flags=-Wl,-rpath,@loader_path"]

    if system == "Windows":
        # ext/sundials_export.h defines SUNDIALS_DEPRECATED with GCC's
        # __attribute__ syntax unconditionally, which MSVC cannot parse. The
        # header guards it with #ifndef, so pre-defining it empty sidesteps
        # the problem. Remove once fixed upstream.
        #
        # cc_flags REPLACES SConstruct's default, so the 'cl' default is
        # repeated here and must be kept in sync. Dropping /MD would switch
        # cantera_shared to the static CRT and give it a private heap.
        msvc_defaults = ("/MD /nologo /D_SCL_SECURE_NO_WARNINGS "
                         "/D_CRT_SECURE_NO_WARNINGS")
        return [
            # SConstruct guesses mingw whenever g++ is on PATH and cl.exe is
            # not -- true on GitHub's Windows images -- and the MSVC flags
            # below would then go to g++.
            "toolchain=msvc",
            f"cc_flags={msvc_defaults} /DSUNDIALS_DEPRECATED=",
        ]

    raise RuntimeError(f"Unsupported platform: {system}")


def environment() -> dict[str, str]:
    """Environment for the SCons subprocess."""
    env = os.environ.copy()

    if platform.system() == "Darwin":
        # Load on older macOS than the build machine.
        env.setdefault("MACOSX_DEPLOYMENT_TARGET", "11.0")

    return env


def stamp_value(cantera_root: Path, flags: list[str]) -> dict:
    """Identify this build by Cantera commit and recipe.

    Hashing the source tree would be accurate but far too slow.
    """
    try:
        sha = subprocess.run(
            ["git", "-C", str(cantera_root), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        sha = "unknown"

    digest = hashlib.sha256("\n".join(sorted(flags)).encode()).hexdigest()[:16]

    return {
        "cantera_commit": sha,
        "flags_sha256": digest,
        "recipe_version": RECIPE_VERSION,
        "platform": f"{platform.system()}-{platform.machine()}",
    }


def is_current(prefix: Path, want: dict) -> bool:
    """Report whether the existing install matches this build's stamp."""
    stamp = prefix / STAMP_NAME
    if not stamp.is_file():
        return False
    try:
        return json.loads(stamp.read_text()) == want
    except (json.JSONDecodeError, OSError):
        return False


def scons_command() -> list[str]:
    """Prefer the SCons module in the active interpreter over a bare script."""
    try:
        subprocess.run([sys.executable, "-m", "SCons", "--version"],
                       capture_output=True, check=True)
        return [sys.executable, "-m", "SCons"]
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    scons = shutil.which("scons")
    if scons is None:
        raise RuntimeError(
            "SCons not found. Install it with 'pip install scons' or make the "
            "'scons' executable available on PATH."
        )
    return [scons]


def run_scons(cantera_root: Path, targets: list[str], flags: list[str],
              jobs: int, verbose: bool) -> None:
    """Run SCons in the Cantera source tree."""
    cmd = scons_command() + targets + [f"-j{jobs}"] + flags
    print(f"[build_cantera_library] {' '.join(cmd)}", flush=True)

    result = subprocess.run(
        cmd, cwd=cantera_root, env=environment(),
        stdout=None if verbose else subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True,
    )
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stderr)
        raise RuntimeError(f"SCons failed with exit code {result.returncode}")


def normalize_layout(prefix: Path) -> None:
    """Adjust the install so ct.buildInterface can consume it directly.

    The compact layout already provides ``<include>/cantera_clib``.
    """
    if platform.system() == "Windows":
        # SCons installs both the DLL and its import library to bin/, while
        # lib/ receives the static library. ctLib.m looks for the DLL in one
        # directory and MSVC needs the import library to link, so copy both
        # into lib/. The .exp is a link byproduct and is not needed.
        bindir, libdir = prefix / "bin", prefix / "lib"
        libdir.mkdir(parents=True, exist_ok=True)
        for pattern in ("*.dll", "*_shared.lib"):
            for artifact in bindir.glob(pattern):
                shutil.copy2(artifact, libdir / artifact.name)
                print(f"[build_cantera_library] staged {artifact.name} "
                      f"into lib/")

    if platform.system() == "Darwin":
        # Locate the library relative to wherever it is unpacked.
        for dylib in (prefix / "lib").glob("libcantera_shared*.dylib"):
            if dylib.is_symlink():
                continue
            subprocess.run(
                ["install_name_tool", "-id", f"@rpath/{dylib.name}", str(dylib)],
                check=True,
            )
            print(f"[build_cantera_library] set install_name for {dylib.name}")


def audit(prefix: Path) -> None:
    """Verify the built library has no non-system dependencies."""
    script = Path(__file__).with_name("audit_library.py")
    subprocess.run([sys.executable, str(script), "--prefix", str(prefix)],
                   check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cantera-root", default=os.environ.get("CANTERA_ROOT"),
                        help="Cantera source checkout (default: $CANTERA_ROOT)")
    parser.add_argument("--prefix", required=True,
                        help="Install prefix for the built library")
    parser.add_argument("--boost-inc-dir", default=os.environ.get("BOOST_INC_DIR"),
                        help="Boost header directory (default: $BOOST_INC_DIR). "
                             "Only needed if Boost >= 1.83 is not on the "
                             "compiler's default include path.")
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    parser.add_argument("--force", action="store_true",
                        help="Rebuild even if the stamp is current")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if not args.cantera_root:
        parser.error("--cantera-root not given and CANTERA_ROOT is not set")

    cantera_root = Path(args.cantera_root).resolve()
    if not (cantera_root / "SConstruct").is_file():
        parser.error(f"No SConstruct in {cantera_root}; not a Cantera checkout")

    problems = (check_submodules(cantera_root)
                + check_prerequisites(args.boost_inc_dir))
    if problems:
        print("error: missing build prerequisites:\n\n"
              + "\n\n".join(f"  - {p}" for p in problems),
              file=sys.stderr)
        return 1

    prefix = Path(args.prefix).resolve()
    flags = (scons_flags(prefix) + platform_flags()
             + boost_flags(args.boost_inc_dir))
    warn_about_inherited_options(cantera_root, flags)
    want = stamp_value(cantera_root, flags)

    if not args.force and is_current(prefix, want):
        print(f"[build_cantera_library] up to date at {prefix}; skipping "
              f"(--force to rebuild)")
        return 0

    if prefix.exists():
        shutil.rmtree(prefix)

    run_scons(cantera_root, ["build"], flags, args.jobs, args.verbose)
    run_scons(cantera_root, ["install"], flags, args.jobs, args.verbose)

    normalize_layout(prefix)
    audit(prefix)

    (prefix / STAMP_NAME).write_text(json.dumps(want, indent=2))
    print(f"[build_cantera_library] self-contained Cantera installed to {prefix}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
