
import subprocess
import logging
import sys
import ete3
import time
import json
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


def set_up_output_dir(output_dir):
    """
    Configure output path
    """
    # set-up output dir
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir()
    elif output_dir.exists():
        logging.error(f"{output_dir} already exists, either remove or change "
                       "output directory")
        sys.exit(1)


def parse_metadata(metadata_fp):
    # nextstrain metadata
    metadata = pd.read_csv(metadata_fp, sep='\t')
    metadata['date'] = pd.to_datetime(metadata['date'].str.replace('-XX', ''))
    metadata = metadata.set_index('strain')
    return metadata


def parse_tree(tree_fp):
    """
    Parse newick tree (from augur or a global phylogeny)
    and rename leaves if needed or data is missing in metadata
    """
    tree = ete3.Tree(tree_fp, format=1)

    # rename tree to use strain names if using Rob Lanfear's full tree
    # disable this as it isn't time corrected so nearest neighbour by
    # phylogenetic distance can get screwy
    #if tree_label_type == 'gisaid_epi_isl':
    #    logging.info(f"Renaming {tree_fp} from gisaid_epi_isl names to strain")
    #    gisaid_to_strain = metadata.reset_index().set_index('gisaid_epi_isl')['strain'].to_dict()
    #    unknown_ix = 0
    #    for leaf in tree.iter_leaves():
    #        try:
    #            leaf.name = gisaid_to_strain[leaf.name]
    #        except KeyError:
    #            logging.warning(f"{leaf.name} missing in metadata, renaming "
    #                            f"to unknown_{unknown_ix}")
    #            leaf.name = f"unknown_{unknown_ix}"
    #            unknown_ix += 1
    return tree


def parse_fasta(fasta_fp, isolates, output_dir):
    """
    extract isolate query fasta from nextfasta if it doesn't already exist
    """
    isolate_fasta_fp = output_dir / "isolate_genomes.fasta"

    fasta_strains = []
    isolate_seqs = []
    for record in SeqIO.parse(fasta_fp, "fasta"):
        fasta_strains.append(record.id)
        if record.id in isolates:
            isolate_seqs.append(record)

    fasta_strains = set(fasta_strains)

    with open(isolate_fasta_fp, 'w') as isolate_fasta_fh:
        SeqIO.write(isolate_seqs, isolate_fasta_fh, "fasta")

    return fasta_strains, isolate_fasta_fp


def check_input_consistency(isolates, metadata, traits, tree, fasta_strains):
    """
    Checks for consistency in input between the tree,
    sequences, and metadata.

    Everything in the tree needs to be in the sequences and metadata
    but not everything in the metadata needs to be in

    Then returns the subset of isolates present in all 3 files
    """

    tree_strains = set([label for label in tree.iter_leaf_names() \
                            if not label.startswith('unknown')])
    tree_root = tree.get_tree_root()
    tree_nodes = set([node.name for node in tree_root.traverse()])

    metadata_strains = set(metadata.index)

    trait_strains = set(traits.keys())

    isolates = set(isolates)

    # check every leaf in tree has metadata (not all metadata has to be
    # in tree though)
    missing_tree_in_metadata = tree_strains - metadata_strains
    if len(missing_tree_in_metadata) > 0:
        logging.error(f"{len(missing_tree_in_metadata)} strains are in tree "
                       "but not in metadata. All tree strains must be in "
                       "metadata (although metadata can contain strains not in "
                       "tree).")
        logging.debug(f"{missing_tree_in_metadata}")
        sys.exit(1)

    # check every node in tree has an entry for traits
    missing_tree_nodes_in_traits = tree_nodes - trait_strains
    if len(missing_tree_nodes_in_traits) > 0:
        logging.error(f"{len(missing_tree_nodes_in_traits)} nodes in tree "
                       "don't have a corresponding entry in traits file "
                       "ensure you are using augur traits defined file for the "
                       "tree you are actually using")
        logging.debug(f"{missing_tree_nodes_in_traits}")
        sys.exit(1)

    # check every leaf in the tree has a fasta seq (but again not all fasta
    # seqs have to be in the tree)
    missing_tree_in_fasta = tree_strains - fasta_strains
    if len(missing_tree_in_fasta) > 0:
        logging.error(f"{len(missing_tree_in_fasta)} strains are in tree "
                       "but not in fasta file. All tree strains must be in "
                       "fasta (although fasta can contain strains not in "
                       "tree).")
        logging.debug(f"{missing_tree_in_fasta}")
        sys.exit(1)

    # check all fasta seqs have metadata
    missing_fasta_metadata = fasta_strains - metadata_strains
    if len(missing_fasta_metadata) > 0:
        logging.error(f"{len(missing_fasta_metadata)} strains are in fasta "
                       "but not in metadata. All fasta strains must have "
                       "an entry in the metadata (although metadata can "
                       "contain strains not in fasta).")
        logging.debug(f"{missing_fasta_metadata}")
        sys.exit(1)


    # check all isolates are in fasta, tree, and metadata
    missing_isolate_in_tree = isolates - tree_strains
    if len(missing_isolate_in_tree) > 0:
        logging.warning(f"{len(missing_isolate_in_tree)} isolates are missing "
                         "in the tree, will only use isolates present in all "
                         "input files.")
        logging.debug("f{missing_isolate_in_tree}")

    missing_isolate_in_metadata = isolates - metadata_strains
    if len(missing_isolate_in_metadata) > 0:
        logging.warning(f"{len(missing_isolate_in_metadata)} isolates are missing "
                         "in the metadata, will only use isolates present in "
                         "all input files.")
        logging.debug("f{missing_isolate_in_metadata}")

    missing_isolate_in_fasta = isolates - fasta_strains
    if len(missing_isolate_in_fasta) > 0:
        logging.warning(f"{len(missing_isolate_in_fasta)} isolates are missing "
                         "in the fasta, will only use isolates present in "
                         "all input files.")
        logging.debug(f"{missing_isolate_in_fasta}")


    # present isolates
    present_isolates = isolates.intersection(tree_strains, metadata_strains,
                                             fasta_strains, trait_strains)
    if len(present_isolates) != len(isolates):
        logging.warning(f"Querying {len(present_isolates)/len(isolates)*100}% due "
                        "to following isolates missing in input data: "
                        f"{isolates-present_isolates}")
        logging.debug(f"Included isolates: {present_isolates}")

    return present_isolates


def mash_sketch_input_fasta(fasta_fp, output_dir):
    # sketch input _files
    sketch = output_dir / Path(fasta_fp).with_suffix('.msh').name
    num_threads = 8
    if not sketch.exists():
        logging.info(f"Sketching {fasta_fp} to {sketch}")
        subprocess.run(f"mash sketch -p {num_threads} -i {fasta_fp} "
                       f"-o {sketch}".split(), check=True)
    elif sketch.exists() and use_existing:
        logging.warning(f"Using previously created sketch: {sketch}")
    elif isolate_fasta_fp.exists() and not use_existing:
        logging.error(f"{sketch} exists but shouldn't as use_existing is"
                       " False")
        sys.exit(1)
    return sketch


def parse_traits(trait_fp):
    """
    Read the augur inferred exposure traits
    """
    with open(trait_fp) as fh:
        traits = json.load(fh)
    traits = traits['nodes']
    return traits


def add_geo_location_from_metadata(closest_df, metadata):
    metadata_location = metadata[['region', 'country', 'division']]
    closest_df = pd.merge(closest_df, metadata_location, left_on='closest_ancestor',
                          right_index=True)
    return closest_df


