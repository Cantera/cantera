#!/usr/bin/env python3

"""
example-keywords.py

Parse Cantera examples for "Keywords" declarations to ensure that all examples have
keyword definitions and to help maintain consistency in the keywords chosen.

Usage:

    example-keywords.py print
        Print a list of keywords found and the number of occurrences of each

    example-keywords.py compare
        Compare the keywords appearing in the examples with the list of known keywords
        from example-keywords.txt. List any that only occur in the examples. Exits with
        an error status code if any keywords are not listed in the known keywords list.
    example-keywords.py save
       Save an updated list of known keywords list, example-keywords.txt
"""

from pathlib import Path
import re
import logging
from collections import Counter
import sys

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
        return False
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
                kw = get_python_keywords(f)
                if kw:
                    all_keywords.update(kw)

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

    return all_keywords


def compare():
    """
    Print a list of keywords appearing in examples that are not in the known keywords
    list. Return True if there are any such items.
    """
    text = (Path(__file__).parent / "example-keywords.txt").read_text()
    known = set(text.split("\n"))
    current = set(get_all_keywords())
    delta = current - known
    for kw in delta:
        logging.warning(f"Keyword {kw!r} not in known keywords list")
    return len(delta) > 0


def save_keywords():
    """
    Save an updated version of the known keywords list based on keywords appearing in
    any of the examples.
    """
    found_kw = get_all_keywords()
    with open(Path(__file__).parent / "example-keywords.txt", "w") as known:
        known.writelines(kw + "\n" for kw in sorted(found_kw))


def print_keywords():
    found_kw = get_all_keywords()
    for kw, count in found_kw.most_common():
        print(f"{kw} ({count})")


if __name__ == "__main__":
    if "compare" in sys.argv:
        delta = compare()
        if delta:
            sys.exit(1)
    elif "save" in sys.argv:
        save_keywords()
    elif "print" in sys.argv:
        print_keywords()
    else:
        print("Valid options are 'print', 'save', or 'compare'")
