
import subprocess
import logging
import sys
import ete3
import time
import tqdm
from pathlib import Path
import pandas as pd
from Bio import SeqIO

def check_dependencies():
    """
    Check all dependencies exist and work
    """
    missing=False
    for program in ['mash']:
        try:
            output = subprocess.run(program, shell=True, check=True,
                    stdout=subprocess.PIPE, encoding='utf-8')
            version = output.stdout.split('\n')[1]
            logging.debug(f"Tool {program} is installed: {version}")
        except:
            logging.error(f"Tool {program} is not installed")
            missing = True

    if missing:
        logging.error("One or more dependencies are missing please install")
        sys.exit(1)
    else:
        logging.debug("All dependencies found")


def set_up_logging(output_dir, debug):
    """
    Handle logging set-up
    """
    run_name = output_dir.name + str(int(time.time()))

    if debug:
        level = logging.DEBUG
    else:
        level = logging.WARNING

    logging.basicConfig(format='%(levelname)s:%(message)s',
                        level=level,
                        handlers=[logging.FileHandler(f"{run_name}.log"),
                                  logging.StreamHandler()])


def set_up_output_dir(output_dir, use_existing):
    """
    Configure output path
    """
    # set-up output dir
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir()
    elif output_dir.exists() and not use_existing:
        logging.error(f"{output_dir} already exists, either remove, change "
                       "output directory, or --use_existing to continue using "
                      f" results from the original {output_dir}")
        sys.exit(1)


def parse_input_files(isolates, metadata_fp, fasta_fp, tree_fp,
                      tree_label_type, output_dir, use_existing):
    """
    Parse input files and create the appropriate data structures
    and files for tree analysis and minimap analysis of closest
    older relatives to query strain
    """

    # nextstrain metadata
    metadata = pd.read_csv(metadata_fp, sep='\t')
    metadata['date'] = pd.to_datetime(metadata['date'].str.replace('-XX', ''))
    metadata = metadata.set_index('strain')

    # refined augur tree
    tree = ete3.Tree(tree_fp, format=1)

    # rename tree to use strain names if using Rob Lanfear's full tree
    if tree_label_type == 'gisaid_epi_isl':
        logging.info(f"Renaming {tree_fp} from gisaid_epi_isl names to strain")
        gisaid_to_strain = metadata.reset_index().set_index('gisaid_epi_isl')['strain'].to_dict()
        unknown_ix = 0
        for leaf in tree.iter_leaves():
            try:
                leaf.name = gisaid_to_strain[leaf.name]
            except KeyError:
                logging.warning(f"{leaf.name} missing in metadata, renaming "
                                f"to unknown_{unknown_ix}")
                leaf.name = f"unknown_{unknown_ix}"
                unknown_ix += 1

    # extract isolate query fasta from nextfasta if it doesn't already exist
    isolate_fasta_fp = output_dir / "isolate_genomes.fasta"
    if isolate_fasta_fp.exists() and use_existing:
        logging.warning(f"Using previously created {isolate_fasta_fp}")
    elif isolate_fasta_fp.exists() and not use_existing:
        logging.error(f"{isolate_fasta_fp} exists but shouldn't as use_existing is"
                       " False")
        sys.exit(1)
    elif not isolate_fasta_fp.exists():
        isolate_seqs = []
        for record in SeqIO.parse(fasta_fp, "fasta"):
            if record.id in isolates:
                isolate_seqs.append(record)

        ## create isolate_fasta
        with open(isolate_fasta_fp, 'w') as isolate_fasta_fh:
            SeqIO.write(isolate_seqs, isolate_fasta_fh, "fasta")

    # check isolates are in metadata and tree
    leaf_names = set([label for label in tree.iter_leaf_names() \
                            if not label.startswith('unknown')])

    for isolate in isolates:
        # checks full metadata table to see where the missing isolate
        # has gone missing
        if isolate not in metadata.index:
            logging.error(f"{isolate} could not be found in {metadata_fp} "
                          f"check strain name provided is correct")
            sys.exit(1)
        # check for augur subsampling
        if isolate not in leaf_names:
            logging.error(f"{isolate} could not be found in tree "
                          f"{tree_fp} check it wasn't removed in "
                           "subsampling by augur modify build or edit "
                           "include file")
            sys.exit(1)

    # sketch input _files
    sketch = output_dir / Path(fasta_fp).with_suffix('.msh').name
    num_threads = 8
    if sketch.exists() and use_existing:
        logging.warning(f"Using previously created sketch: {sketch}")
    elif isolate_fasta_fp.exists() and not use_existing:
        logging.error(f"{sketch} exists but shouldn't as use_existing is"
                       " False")
        sys.exit(1)
    elif not sketch.exists():
        logging.info(f"Sketching {fasta_fp} to {sketch}")
        subprocess.run(f"mash sketch -p {num_threads} -i {fasta_fp} "
                       f"-o {sketch}".split(), check=True)

    return isolate_fasta_fp, sketch, metadata, tree
