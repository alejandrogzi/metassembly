#!/usr/bin/env python3

"""
DEACON Pan-genome Index Builder

A command-line tool for building pan-genome indexes using DEACON with optional
background datasets (FDA-ARGOS and RefSeq viral) for differential analysis.
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Final, Sequence

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.1"

# Constants
FDA_ARGOS_URL: Final[str] = (
    "https://zenodo.org/records/15424142/files/argos988.fa.zst?download=1"
)
REFSEQ_VIRAL_URL: Final[str] = (
    "https://zenodo.org/records/15411280/files/rsviruses17900.fa.gz?download=1"
)
DEFAULT_KMER_LENGTH: Final[int] = 31
DEFAULT_WINDOW_SIZE: Final[int] = 15
DEFAULT_ENTROPY: Final[float] = 0.5
DEFAULT_THREADS: Final[int] = 8


class DeaconIndexBuilder:
    """Handles building and managing DEACON indexes for pan-genome analysis."""

    def __init__(
        self,
        outdir: Path,
        kmer_length: int = DEFAULT_KMER_LENGTH,
        window_size: int = DEFAULT_WINDOW_SIZE,
        entropy: float = DEFAULT_ENTROPY,
        threads: int = DEFAULT_THREADS,
    ) -> None:
        """
        Initialize the DEACON index builder.

        Parameters
        ----------
        outdir : Path
            Output directory for indexes and downloaded files
        kmer_length : int, optional
            K-mer length for index building (default: 31)
        window_size : int, optional
            Window size for index building (default: 15)
        entropy : float, optional
            Entropy threshold for filtering (default: 0.5)
        threads : int, optional
            Number of threads for parallel processing (default: 8)
        """
        self.outdir = outdir
        self.kmer_length = kmer_length
        self.window_size = window_size
        self.entropy = entropy
        self.threads = threads
        self.outdir.mkdir(exist_ok=True, parents=True)

    def build_genome_index(self, genome_path: Path) -> Path:
        """
        Build a DEACON index for a single genome.

        Parameters
        ----------
        genome_path : Path
            Path to the genome FASTA file

        Returns
        -------
        Path
            Path to the generated index file
        """
        output = self.outdir / f"{genome_path.name}.idx"
        cmd = (
            f"deacon index build "
            f"-k {self.kmer_length} "
            f"-w {self.window_size} "
            f"-e {self.entropy} "
            f"-t {self.threads} "
            f"{genome_path} > {output}"
        )
        self._run_command(cmd)
        return output

    def create_union_index(self, indexes: Sequence[Path], output_name: str) -> Path:
        """
        Create a union of multiple DEACON indexes.

        Parameters
        ----------
        indexes : Sequence[Path]
            List of index paths to union
        output_name : str
            Name for the output union index

        Returns
        -------
        Path
            Path to the generated union index
        """
        output = self.outdir / output_name
        index_paths = " ".join(str(idx) for idx in indexes)
        cmd = f"deacon index union {index_paths} --output {output}"
        self._run_command(cmd)
        return output

    def create_diff_index(self, main_index: Path, background_index: Path) -> Path:
        """
        Create a differential index by subtracting background from main index.

        Parameters
        ----------
        main_index : Path
            Path to the main genome index
        background_index : Path
            Path to the background index

        Returns
        -------
        Path
            Path to the generated differential index
        """
        output = self.outdir / "final.idx"
        cmd = f"deacon index diff {main_index} {background_index} --output {output}"
        self._run_command(cmd)
        return output

    def download_and_index_fda_argos(self) -> Path:
        """
        Download FDA-ARGOS dataset and build its index.

        Returns
        -------
        Path
            Path to the FDA-ARGOS index
        """
        print("INFO: Downloading FDA-ARGOS dataset")
        zst_output = self.outdir / "argos988.fa.zst"
        gz_output = self.outdir / "argos988.fa.gz"

        download_cmd = (
            f"wget -O {zst_output} {FDA_ARGOS_URL} && "
            f"zstdcat {zst_output} | pigz > {gz_output}"
        )
        self._run_command(download_cmd)

        print("INFO: Building index for FDA-ARGOS dataset")
        argos_idx_output = self.outdir / "argos988.idx"
        cmd = (
            f"deacon index build "
            f"-k {DEFAULT_KMER_LENGTH} "
            f"-w {DEFAULT_WINDOW_SIZE} "
            f"{gz_output} --output {argos_idx_output}"
        )
        self._run_command(cmd)
        return argos_idx_output

    def download_and_index_refseq_viral(self) -> Path:
        """
        Download RefSeq viral dataset and build its index.

        Returns
        -------
        Path
            Path to the RefSeq viral index
        """
        print("INFO: Downloading RefSeq viral dataset")
        output = self.outdir / "rsviruses17900.fa.gz"
        download_cmd = f"wget -O {output} {REFSEQ_VIRAL_URL}"
        self._run_command(download_cmd)

        print("INFO: Building index for RefSeq viral dataset")
        viral_idx_output = self.outdir / "rsviruses17900.idx"
        cmd = (
            f"deacon index build "
            f"-k {DEFAULT_KMER_LENGTH} "
            f"-w {DEFAULT_WINDOW_SIZE} "
            f"{output} --output {viral_idx_output}"
        )
        self._run_command(cmd)
        return viral_idx_output

    @staticmethod
    def _run_command(cmd: str) -> None:
        """
        Execute a shell command and print it for logging.

        Parameters
        ----------
        cmd : str
            Shell command to execute
        """
        print(f"INFO: Running -> {cmd}")
        subprocess.run(cmd, shell=True, check=True)


def validate_genome_paths(genome_paths: Sequence[Path]) -> None:
    """
    Validate that all genome files exist.

    Parameters
    ----------
    genome_paths : Sequence[Path]
        List of genome file paths to validate

    Raises
    ------
    SystemExit
        If any genome file does not exist
    """
    for genome in genome_paths:
        if not genome.exists():
            print(f"ERROR: Genome {genome} does not exist", file=sys.stderr)
            sys.exit(1)


def build_main_indexes(
    builder: DeaconIndexBuilder, genome_paths: Sequence[Path]
) -> Path:
    """
    Build indexes for all provided genomes and create union if multiple.

    Parameters
    ----------
    builder : DeaconIndexBuilder
        Index builder instance
    genome_paths : Sequence[Path]
        List of genome file paths

    Returns
    -------
    Path
        Path to the main index (union if multiple, single if one)

    Raises
    ------
    SystemExit
        If no indexes were created
    """
    indexes = [builder.build_genome_index(genome) for genome in genome_paths]

    if len(indexes) > 1:
        print(f"INFO: Found {len(indexes)} indexes. Creating multindex")
        main_index = builder.create_union_index(indexes, "multindex.idx")
        print(f"INFO: Multindex created at {main_index}")
    elif len(indexes) == 1:
        main_index = indexes[0]
        print(f"INFO: Index created at {main_index}")
    else:
        print("ERROR: No indexes created!", file=sys.stderr)
        sys.exit(1)

    return main_index


def build_background_index(
    builder: DeaconIndexBuilder, use_fda_argos: bool, use_refseq_viral: bool
) -> Path | None:
    """
    Build background index from optional datasets.

    Parameters
    ----------
    builder : DeaconIndexBuilder
        Index builder instance
    use_fda_argos : bool
        Whether to include FDA-ARGOS dataset
    use_refseq_viral : bool
        Whether to include RefSeq viral dataset

    Returns
    -------
    Path | None
        Path to background index if any datasets were used, None otherwise
    """
    background_indexes = []

    if use_fda_argos:
        background_indexes.append(builder.download_and_index_fda_argos())

    if use_refseq_viral:
        background_indexes.append(builder.download_and_index_refseq_viral())

    if len(background_indexes) > 1:
        print("INFO: Creating background multindex from multiple datasets")
        background = builder.create_union_index(background_indexes, "background.idx")
    elif len(background_indexes) == 1:
        background = background_indexes[0]
    else:
        return None

    print(f"INFO: Background created at {background}")
    return background


def run(args: argparse.Namespace) -> None:
    """
    Execute the pan-genome index building pipeline.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments
    """
    genome_paths = [Path(g) for g in args.genomes]
    validate_genome_paths(genome_paths)

    outdir = Path(args.outdir)
    print(f"Output directory: {outdir}")

    builder = DeaconIndexBuilder(
        outdir=outdir,
        kmer_length=args.kmer_length,
        window_size=args.window_size,
        entropy=args.entropy,
        threads=args.threads,
    )

    # Build main indexes
    main_index = build_main_indexes(builder, genome_paths)

    # Build background if requested
    background = build_background_index(
        builder, args.use_fda_argos, args.use_refseq_viral
    )

    # Create differential index if background exists
    if background:
        pan_output = builder.create_diff_index(main_index, background)
        print(f"INFO: Pan-genome index created at {pan_output}")
    else:
        print(f"INFO: No background. Final index at {main_index}")


class SplitArgsAction(argparse.Action):
    """Custom argparse action to split comma-separated values."""

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: str,
        option_string: str | None = None,
    ) -> None:
        """Split comma-separated string into list."""
        setattr(namespace, self.dest, values.split(","))


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments

    Examples
    --------
    >>> args = parse_arguments()
    >>> print(args.genomes)
    ['genome1.fa', 'genome2.fa']
    """
    parser = argparse.ArgumentParser(
        description="Build pan-genome indexes using DEACON with optional background datasets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -g genome1.fa,genome2.fa -o output/
  %(prog)s -g genome.fa -F -R -k 31 -w 15 -t 16
        """,
    )

    parser.add_argument(
        "-g",
        "--genomes",
        required=True,
        action=SplitArgsAction,
        metavar="GENOME1,GENOME2,...",
        help="Comma-separated paths to genome FASTA files",
    )
    parser.add_argument(
        "-e",
        "--entropy",
        type=float,
        default=DEFAULT_ENTROPY,
        metavar="FLOAT",
        help=f"Entropy threshold for filtering (default: {DEFAULT_ENTROPY})",
    )
    parser.add_argument(
        "-k",
        "--kmer-length",
        type=int,
        default=DEFAULT_KMER_LENGTH,
        metavar="INT",
        help=f"K-mer length for index building (default: {DEFAULT_KMER_LENGTH})",
    )
    parser.add_argument(
        "-w",
        "--window-size",
        type=int,
        default=DEFAULT_WINDOW_SIZE,
        metavar="INT",
        help=f"Window size for index building (default: {DEFAULT_WINDOW_SIZE})",
    )
    parser.add_argument(
        "-F",
        "--use-fda-argos",
        action="store_true",
        help="Include FDA-ARGOS dataset as background for differential analysis",
    )
    parser.add_argument(
        "-R",
        "--use-refseq-viral",
        action="store_true",
        help="Include RefSeq viral dataset as background for differential analysis",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        metavar="INT",
        help=f"Number of threads for parallel processing (default: {DEFAULT_THREADS})",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        default=".",
        metavar="DIR",
        help="Output directory for indexes and downloads (default: current directory)",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser.parse_args()


def main() -> None:
    """Main entry point for the application."""
    args = parse_arguments()
    run(args)


if __name__ == "__main__":
    main()
