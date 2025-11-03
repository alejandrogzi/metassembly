#!/usr/bin/env python

import argparse

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.1"


def run(args: argparse.Namespace) -> None:
    """
    Run the renaming process.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    with open(args.gtf) as f:
        with open(args.output, "w") as out:
            for line in f:
                if line.startswith("#"):
                    out.write(line)
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 8:
                    continue

                attributes = {
                    attr.strip().split(" ")[0]: attr.strip().split(" ")[1]
                    for attr in fields[8].split(";")
                    if len(attr.strip().split(" ")) > 1
                }

                attributes["gene_id"] = (
                    f"\"{args.prefix}_{attributes['gene_id'].strip('"')}\""
                )
                attributes["transcript_id"] = (
                    f"\"{args.prefix}_{attributes['transcript_id'].strip('"')}\""
                )

                fields[8] = ";".join(
                    [
                        f" {k} {v}" if k != "gene_id" else f"{k} {v}"
                        for k, v in attributes.items()
                    ]
                )
                out.write("\t".join(fields))
                out.write("\n")

    print(f"INFO: Renamed GTF file to {args.gtf}")
    return


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments

    Examples
    --------
    >>> args = parse_arguments()
    """
    parser = argparse.ArgumentParser(
        description="Rename gene_id and transcript_id in a GTF file using a prefix",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-g",
        "--gtf",
        required=True,
        metavar="GTF",
        help="GTF file",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=True,
        metavar="STR",
        help="Prefix for the every line of the GTF file",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        metavar="GTF",
        help="Path to output GTF file",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run(args)


if __name__ == "__main__":
    main()
