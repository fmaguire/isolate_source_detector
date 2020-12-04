"""Console script for isolate_source_detector."""
import sys
import click
from isolate_source_detector import isolate_source_detector

@click.command()
@click.argument("isolates", nargs=-1)
@click.option("--metadata", "-m", required=True,
              help="Nextstrain metadata file", type=click.Path(exists=True))
@click.option("--fasta", "-f", required=True,
              help="Nextstrain fasta file", type=click.Path(exists=True))
@click.option("--augur_tree", "-t", required=True,
              help="Augur refined phylogeny", type=click.Path(exists=True))
@click.option("--output_dir", "-o", default="out",
              help="Output directory")
@click.option("--overwrite", is_flag=True, help="Overwrite previous output")
def main(isolates, fasta, metadata, augur_tree, output_dir, overwrite):
    """Console script for isolate_source_detector"""
    isolate_source_detector.isolate_source_detector(isolates, fasta, metadata,
                                                    augur_tree, output_dir,
                                                    overwrite)

    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
