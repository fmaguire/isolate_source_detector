"""Console script for isolate_source_detector."""
import sys
import click
from pathlib import Path
from isolate_source_detector import isolate_source_detector

@click.command()
@click.argument("isolates", nargs=-1)
@click.option("--metadata", "-m", required=True,
              help="Nextstrain metadata file", type=click.Path(exists=True))
@click.option("--fasta", "-f", required=True,
              help="Nextstrain fasta file", type=click.Path(exists=True))
@click.option("--tree", "-t", required=True,
              help="Augur refined phylogeny or global phylogeny",
              type=click.Path(exists=True))
@click.option("--tree_label_type", "-l", default='strain',
              help="Metadata field used to name the supplied tree e.g. "
                   "'strain' for augur trees or 'gisaid_epi_isl' for "
                   " Rob Lanfear's global phylogeny",
              type=click.Choice(['strain', 'gisaid_epi_isl']))
@click.option("--output_dir", "-o", default=Path("out"), type=click.Path(),
              help="Output directory")
@click.option("--debug", is_flag=True, help="Log debugging information")
@click.option("--use_existing", is_flag=True,
              help="Continue with files from previous runs")
def main(isolates, metadata, fasta, tree, tree_label_type, output_dir,
         debug, use_existing):
    """Console script for isolate_source_detector"""

    isolate_source_detector.isd(isolates, metadata, fasta, tree,
                                tree_label_type, output_dir, debug,
                                use_existing)
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
