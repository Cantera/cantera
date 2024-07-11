"""SDist of the Python Module"""

from pathlib import Path
import re
from textwrap import indent
import argparse
import shutil
import sys


sys.path.append("site_scons/site_tools")
from UnitsInterfaceBuilder import get_property_definition_strings


def replace_git_hash(source: Path, target: Path, git_commit: str) -> None:
    # Avoid having to set a C preprocessor define at compile time, since
    # the git commit is unknown from the sdist
    git_commit_replaced = re.sub(
        "#ifdef GIT_COMMIT.*?#endif",
        f"""    return "{git_commit}";""",
        source.read_text(),
        flags=re.DOTALL,
    )
    target.write_text(git_commit_replaced)


def _substitute_lines(lines: list[str], substitutions: dict[str, str]) -> list[str]:
    output_lines: list[str] = []
    for line in lines:
        for key, value in substitutions.items():
            if key in line:
                line = line.replace(key, value)
                break

        output_lines.append(line)
    return output_lines


def do_configure_substitution(
    configure_source: Path, cantera_version: str, cantera_short_version: str
) -> str:
    configure_template = configure_source.read_text().splitlines()
    configure_subst = {
        "@cantera_version@": cantera_version,
        "@cantera_short_version@": cantera_short_version,
    }
    configure_output = _substitute_lines(configure_template, configure_subst)
    return "\n".join(configure_output)


def do_pyproject_substitution(
    pyproject_toml_source: Path, py_requires_ver_str: str, cantera_version: str
) -> str:
    pyproject_toml_template = pyproject_toml_source.read_text().splitlines()
    pyproject_subst = {
        "@py_requires_ver_str@": py_requires_ver_str,
        "@cantera_version@": cantera_version,
    }
    pyproject_toml_output = _substitute_lines(pyproject_toml_template, pyproject_subst)
    return "\n".join(pyproject_toml_output)


def do_units_substitution(with_units_source: Path) -> str:
    common, thermophase, purefluid = get_property_definition_strings()
    substitution_dict = {
        "@common_properties@": indent(common, " " * 4),
        "@thermophase_properties@": indent(thermophase, " " * 4),
        "@purefluid_properties@": indent(purefluid, " " * 4),
    }
    with_units_template = with_units_source.read_text().splitlines()
    with_units_output = _substitute_lines(with_units_template, substitution_dict)
    return "\n".join(with_units_output)


def main(
    source_directory: Path,
    target_directory: Path,
    git_commit: str,
    py_requires_ver_str: str,
    cantera_version: str,
    cantera_short_version: str,
):
    src_source = source_directory / "src"
    src_target = target_directory / "src"
    shutil.copytree(
        src_source,
        src_target,
        ignore=shutil.ignore_patterns(
            "fortran",
            "clib",
            "pch",
            "global.cpp",
            "SCons*",
            "canteraStatic.cpp",
        ),
        dirs_exist_ok=True,
    )
    base_global_cpp = Path("base", "global.cpp")
    replace_git_hash(
        src_source / base_global_cpp, src_target / base_global_cpp, git_commit
    )

    include_source = source_directory / "include"
    include_target = target_directory / "include"
    # config.h needs to be filled on the user's machine by CMake. utils_utils.h is
    # handled by the replace_git_hash command.
    shutil.copytree(
        include_source,
        include_target,
        ignore=shutil.ignore_patterns(
            "clib",
            "ext",
            "config.h",
            "utils_utils.h",
        ),
        dirs_exist_ok=True,
    )
    utils_utils_h = Path("cantera", "cython", "utils_utils.h")
    replace_git_hash(
        include_source / utils_utils_h,
        include_target / utils_utils_h,
        git_commit,
    )

    cantera_source = source_directory / "interfaces" / "cython" / "cantera"
    cantera_target = target_directory / "cantera"
    shutil.copytree(
        cantera_source,
        cantera_target,
        ignore=shutil.ignore_patterns("__pycache__", "solution.py.in"),
        dirs_exist_ok=True,
    )

    with_units = Path("with_units", "solution.py")
    with_units_source = cantera_source / with_units.with_suffix(".py.in")
    with_units_target = cantera_target / with_units
    with_units_target.write_text(do_units_substitution(with_units_source))

    data_source = source_directory / "data"
    data_target = target_directory / "cantera" / "data"
    shutil.copytree(
        data_source,
        data_target,
        ignore=shutil.ignore_patterns("example_data"),
        dirs_exist_ok=True,
    )
    shutil.copytree(
        data_source / "example_data",
        data_target / "example_data",
        ignore=shutil.ignore_patterns(".git"),
        dirs_exist_ok=True,
    )

    HERE = Path(__file__).parent
    pyproject_toml_source = HERE / "pyproject.toml.in"
    pyproject_toml_target = target_directory / "pyproject.toml"
    pyproject_toml_target.write_text(
        do_pyproject_substitution(
            pyproject_toml_source, py_requires_ver_str, cantera_version
        )
    )

    configure_source = HERE / "src" / "configure.py.in"
    configure_target = target_directory / "src" / "configure.py"
    configure_target.write_text(
        do_configure_substitution(
            configure_source, cantera_version, cantera_short_version
        )
    )

    for cmakelist in ("", "cantera", "src"):
        location = Path(cmakelist, "CMakeLists.txt")
        shutil.copy2(HERE / location, target_directory / location)
    shutil.copy2(source_directory / "README.rst", target_directory / "README.rst")
    shutil.copy2(source_directory / "License.txt", target_directory / "License.txt")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_directory")
    parser.add_argument("target_directory")
    parser.add_argument("git_commit")
    parser.add_argument("python_version_req")
    parser.add_argument("cantera_version")
    parser.add_argument("cantera_short_version")
    args = parser.parse_args(argv)
    return args


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(
        Path(args.source_directory),
        Path(args.target_directory),
        args.git_commit,
        args.python_version_req,
        args.cantera_version,
        args.cantera_short_version,
    )
