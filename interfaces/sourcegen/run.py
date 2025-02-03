# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
import argparse
import textwrap

import sourcegen

def main(argv=None):
    parser = create_argparser()
    if argv is None and len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args(argv)
    lang = args.api
    output = args.output
    verbose = args.verbose
    sourcegen.generate_source(lang, output, verbose=verbose)

def create_argparser():
    parser = argparse.ArgumentParser(
        description=(
            "Experimental source generator for creating Cantera interface code."),
        epilog=textwrap.dedent(
            """
            The **sourcegen** utility is invoked as follows::

                python path/to/sourcegen/run.py --api=csharp --output=.

            where the relative path has to be provided as the utility is not installed.
            Currently supported API options are: 'csharp', 'clib' and 'yaml'.
            """),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="show additional logging output")
    parser.add_argument(
        "--api", default="",
        help="language of generated Cantera API code")
    parser.add_argument(
        "--output", default="",
        help="specifies the OUTPUT folder name")

    return parser

if __name__ == "__main__":
    main()
