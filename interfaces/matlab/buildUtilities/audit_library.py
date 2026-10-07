#!/usr/bin/env python3
"""Verify that a built Cantera library depends only on system libraries.

Inspects the dynamic dependencies of a built cantera_shared and fails if
anything outside a per-platform allowlist appears.

Run standalone against any prefix:
    python audit_library.py --prefix build/canteraLib
"""

from __future__ import annotations

import argparse
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Libraries that are part of the OS (or of MATLAB's guaranteed runtime).
# libstdc++ and libgcc_s are deliberately absent on Linux: the build links them
# statically, so their appearance means -static-libstdc++ did not take effect.
ALLOWLISTS = {
    "Linux": [
        r"^linux-vdso\.so\.1$",
        r"^libc\.so\.6$",
        r"^libm\.so\.6$",
        r"^libdl\.so\.2$",
        r"^libpthread\.so\.0$",
        r"^librt\.so\.1$",
        r"^ld-linux-.*\.so\.\d+$",
    ],
    "Darwin": [
        r"^/usr/lib/libSystem\.B\.dylib$",
        r"^/usr/lib/libc\+\+\.1\.dylib$",
        r"^/usr/lib/libobjc\.A\.dylib$",
        r"^/System/Library/Frameworks/Accelerate\.framework/",
        # The library's own id and its co-located siblings.
        r"^@rpath/libcantera_shared",
        r"^@loader_path/",
    ],
    "Windows": [
        r"^KERNEL32\.dll$",
        r"^ADVAPI32\.dll$",
        r"^USER32\.dll$",
        r"^SHELL32\.dll$",
        r"^ole32\.dll$",
        r"^OLEAUT32\.dll$",
        r"^bcrypt\.dll$",
        r"^dbgeng\.dll$",
        r"^api-ms-win-.*\.dll$",
        # The MSVC redistributable, which MATLAB ships.
        r"^MSVCP140.*\.dll$",
        r"^VCRUNTIME140.*\.dll$",
    ],
}


def find_library(prefix: Path) -> Path:
    """Locate the real (non-symlink) shared library under the prefix."""
    patterns = {
        "Linux": "libcantera_shared.so*",
        "Darwin": "libcantera_shared*.dylib",
        "Windows": "cantera_shared.dll",
    }[platform.system()]

    candidates = [p for p in (prefix / "lib").glob(patterns) if not p.is_symlink()]
    if not candidates:
        raise FileNotFoundError(
            f"No cantera_shared library found in {prefix / 'lib'}"
        )
    # Prefer the fully versioned file, which is the real object on Linux/macOS.
    return max(candidates, key=lambda p: len(p.name))


def dependencies(library: Path) -> list[str]:
    """Return the library's direct dynamic dependencies, as recorded names."""
    system = platform.system()

    if system == "Linux":
        # readelf reports DT_NEEDED. ldd would resolve the transitive graph
        # through this machine's paths and hide a missing dependency behind a
        # locally installed copy.
        out = _run(["readelf", "-d", str(library)])
        return re.findall(r"Shared library: \[([^\]]+)\]", out)

    if system == "Darwin":
        # First line is the file name; the rest are deps, starting with the
        # library's own install id.
        out = _run(["otool", "-L", str(library)])
        return [line.strip().split(" (")[0]
                for line in out.splitlines()[1:] if line.startswith("\t")]

    if system == "Windows":
        # Parsed in-process; dumpbin only exists in a Visual Studio developer
        # shell, and CI runs these steps under plain bash.
        return pe_imports(library)

    raise RuntimeError(f"Unsupported platform: {system}")


def pe_imports(library: Path) -> list[str]:
    """Read the DLL names from a PE file's import directory."""
    data = library.read_bytes()

    def u16(off: int) -> int:
        return int.from_bytes(data[off:off + 2], "little")

    def u32(off: int) -> int:
        return int.from_bytes(data[off:off + 4], "little")

    if data[:2] != b"MZ":
        raise ValueError(f"{library} is not a PE image")

    pe = u32(0x3C)
    if data[pe:pe + 4] != b"PE\0\0":
        raise ValueError(f"{library} has no PE signature")

    coff = pe + 4
    n_sections = u16(coff + 2)
    opt_size = u16(coff + 16)
    opt = coff + 20

    # Data directories sit at offset 96 in a PE32 optional header, 112 in
    # PE32+, where several fields widen to 8 bytes.
    magic = u16(opt)
    dir_off = opt + (96 if magic == 0x10B else 112)
    import_rva = u32(dir_off + 8 * 1)  # data directory 1 = imports
    if import_rva == 0:
        return []

    sections = []
    sec_table = opt + opt_size
    for i in range(n_sections):
        s = sec_table + 40 * i
        sections.append((u32(s + 12), u32(s + 8), u32(s + 20)))  # va, vsize, raw

    def to_offset(rva: int) -> int:
        for va, vsize, raw in sections:
            if va <= rva < va + max(vsize, 1):
                return raw + (rva - va)
        raise ValueError(f"RVA {rva:#x} is outside every section")

    def cstring(rva: int) -> str:
        start = to_offset(rva)
        end = data.index(b"\0", start)
        return data[start:end].decode("ascii", "replace")

    names, entry = [], to_offset(import_rva)
    while True:
        # IMAGE_IMPORT_DESCRIPTOR is 20 bytes with the name RVA at offset 12;
        # an all-zero descriptor terminates the table.
        descriptor = data[entry:entry + 20]
        if len(descriptor) < 20 or descriptor == b"\0" * 20:
            break
        names.append(cstring(u32(entry + 12)))
        entry += 20

    return names


def _run(cmd: list[str]) -> str:
    """Run a required inspection tool and return its stdout."""
    if shutil.which(cmd[0]) is None:
        raise RuntimeError(
            f"'{cmd[0]}' not found; it is required to audit the built library."
        )
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--prefix", help="Install prefix containing lib/")
    source.add_argument("--library",
                        help="Audit this library file directly, e.g. the copy "
                             "clibgen placed in the staged ctMatlab folder.")
    args = parser.parse_args()

    if args.library:
        library = Path(args.library).resolve()
        if not library.is_file():
            print(f"error: no such file: {library}", file=sys.stderr)
            return 1
    else:
        library = find_library(Path(args.prefix).resolve())
    allowed = [re.compile(p) for p in ALLOWLISTS[platform.system()]]

    deps = dependencies(library)
    unexpected = [d for d in deps if not any(p.search(d) for p in allowed)]

    print(f"[audit_library] {library.name} depends on {len(deps)} libraries:")
    for dep in deps:
        mark = "  " if dep not in unexpected else "->"
        print(f"  {mark} {dep}")

    if unexpected:
        sys.stdout.flush()  # keep the listing above the error in CI logs
        print(
            f"\nerror: {library.name} has {len(unexpected)} dependency/ies "
            f"outside the {platform.system()} allowlist:\n"
            + "\n".join(f"  {d}" for d in unexpected)
            + "\n\nEither vendor them into the build or, if genuinely part of "
              "the OS, add them to ALLOWLISTS in this file.",
            file=sys.stderr,
        )
        return 1

    print(f"\n[audit_library] OK: no dependencies outside the "
          f"{platform.system()} system allowlist.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
