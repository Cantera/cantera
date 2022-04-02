#!/usr/bin/env python3

from pathlib import Path
import re
import logging
from collections import Counter

def get_python_keywords(filename):
    """
    Get keywords defined in the file docstring of a Python example. Keywords
    start after a line starting with "Keywords:" and end with either the end of
    the docstring or with a blank line. Individual keywords are comma separated.
    """
    text = Path(filename).read_text()
    match = re.search(r"(\"\"\"|\'\'\')(.*?)\1", text, re.DOTALL | re.MULTILINE)
    if not match:
        logging.error(f"Couldn't parse docstring for {filename}")
        return False
    docstring = match.group(2) + "\n\n"

    match = re.search(r"\s*Keywords:(.*?)\n\n", docstring, re.DOTALL | re.MULTILINE)
    if not match:
        logging.warning(f"No keywords found in {filename}")
    keywords = set(kw.strip() for kw in match[1].split(","))
    return keywords


def get_matlab_keywords(filename):
    """
    Get keywords from the first comment block of a Matlab example. Keywords start after
    a line starting with "Keywords:" and end with either the end of  the comment block
    or with a line only containing the comment character (%). Individual keywords are
    comma separated.
    """
    text = Path(filename).read_text()
    match = re.search(r"(?:%.*?\n)+", text, re.DOTALL | re.MULTILINE)
    if not match:
        logging.error(f"Couldn't parse docstring for {filename}")
        return False
    docstring = match.group(0) + "\n\n"
    docstring = "\n".join(line.lstrip("% ") for line in docstring.split("\n"))

    match = re.search(r"\s*Keywords:(.*?)\n\n", docstring, re.DOTALL | re.MULTILINE)
    if not match:
        logging.warning(f"No keywords found in {filename}")
        return False
    keywords = set(kw.strip() for kw in match[1].split(","))
    return keywords


def get_cxx_keywords(filename):
    text = Path(filename).read_text()
    match = re.search(r"\/\*[!\*](.*?)\*\/", text, re.DOTALL | re.MULTILINE)
    if not match:
        logging.error(f"Couldn't parse docstring for {filename}")
        return False
    docstring = match.group(1) + "\n\n"
    docstring = "\n".join(line.lstrip("* ") for line in docstring.split("\n"))
    match = re.search(r"\s*Keywords:(.*?)\n\n", docstring, re.DOTALL | re.MULTILINE)
    if not match:
        logging.warning(f"No keywords found in {filename}")
        return False
    keywords = set(kw.strip() for kw in match[1].split(","))
    return keywords


def get_fortran_keywords(filename, comment_char):
    text = Path(filename).read_text()
    match = re.search(fr"(?:{comment_char}.*?\n)+", text,
                      re.DOTALL | re.MULTILINE | re.IGNORECASE)
    if not match:
        logging.error(f"Couldn't parse docstring for {filename}")
        return False
    docstring = match.group(0) + "\n\n"
    docstring = "\n".join(line.lstrip(f"{comment_char} ")
                          for line in docstring.split("\n"))

    match = re.search(r"\s*Keywords:(.*?)\n\n", docstring, re.DOTALL | re.MULTILINE)
    if not match:
        logging.warning(f"No keywords found in {filename}")
        return False
    keywords = set(kw.strip() for kw in match[1].split(","))
    return keywords


def get_all_keywords():
    """
    Read keywords from all Cantera examples and print out a summary list
    """
    all_keywords = Counter()
    for d in Path("interfaces/cython/cantera/examples").glob("**"):
        if d.is_dir():
            for f in d.glob("*.py"):
                if str(f.name).startswith("_"):
                    continue
                all_keywords.update(get_python_keywords(f))

    for f in Path("samples/matlab").glob("*.m"):
        kw = get_matlab_keywords(f)
        if kw:
            all_keywords.update(kw)

    for d in Path("samples/cxx").glob("**"):
        if d.is_dir():
            for f in d.glob("*.cpp"):
                all_keywords.update(get_cxx_keywords(f))

    for f in Path("samples/f77").glob("*.f"):
        kw = get_fortran_keywords(f, "c")
        if kw:
            all_keywords.update(kw)

    for f in Path("samples/f90").glob("*.f90"):
        kw = get_fortran_keywords(f, "!")
        if kw:
            all_keywords.update(kw)

    for kw, count in sorted(all_keywords.items()):
        print(f"{kw} ({count})")


if __name__ == "__main__":
    get_all_keywords()
