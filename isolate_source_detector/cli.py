"""Console script for isolate_source_detector."""
import sys
import click
from pathlib import Path
from isolate_source_detector import __version__
from isolate_source_detector import isolate_source_detector

def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"isolate_source_detector version {__version__}")
    ctx.exit()

@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option("--isolates", "-i", required=True, type=click.Path(exists=True),
              help="File containing list of isolate strain names (one per line)")
@click.option("--metadata", "-m", required=True,
              help="Nextstrain metadata file", type=click.Path(exists=True))
@click.option("--fasta", "-f", required=True,
              help="Nextstrain fasta file", type=click.Path(exists=True))
@click.option("--tree", "-t", required=True,
              help="Augur refined time-corrected phylogeny",
              type=click.Path(exists=True))
#@click.option("--tree_label_type", "-l", default='strain',
#              help="Metadata field used to name the supplied tree e.g. "
#                   "'strain' for augur trees or 'gisaid_epi_isl' for "
#                   " Rob Lanfear's global phylogeny",
#              type=click.Choice(['strain', 'gisaid_epi_isl']))
@click.option("--traits", "-j",
              help="Augur inferred traits json e.g. traits.json",
              required=True,
              type=click.Path(exists=True))
@click.option("--output_dir", "-o", default=Path("out"), type=click.Path(),
              help="Output directory")
@click.option("--debug", is_flag=True, help="Log debugging information")
#@click.option("--use_existing", is_flag=True,
#              help="Continue with files from previous runs")
@click.option('--version', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
def main(isolates, metadata, fasta, tree, traits, output_dir,
         debug):
    """Console script for isolate_source_detector"""

    isolate_source_detector.isd(isolates, metadata, fasta, tree,
                                traits, output_dir, debug)
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
