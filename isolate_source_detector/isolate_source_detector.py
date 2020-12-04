"""Main module."""

from Bio import SeqIO
import pandas as pd
import ete3
import sys
import shutil
import tqdm
import logging
import time
from pathlib import Path


def parse_input_files(nextfasta, nextmeta, augur_tree_fp, isolates, output_dir):
    """
    Parse input files and create the appropriate data structures
    and files for tree analysis and minimap analysis of closest
    older relatives to query strain
    """

    # nextstrain metadata
    metadata = pd.read_csv(nextmeta, sep='\t')
    metadata['date'] = pd.to_datetime(metadata['date'].str.replace('-XX', ''))
    metadata = metadata.set_index('strain')

    # refined augur tree
    augur_tree = ete3.Tree(augur_tree_fp, format=1)

    # extract query fasta from nextfasta
    #isolate_seqs = []
    #for record in SeqIO.parse(nextfasta, "fasta"):
    #    if record.id in isolates:
    #        isolate_seqs.append(record)

    ## create isolate_fasta
    isolate_fasta = output_dir / "isolate_genomes.fasta"
    #with open(isolate_fasta, 'w') as isolate_fasta_fh:
    #    SeqIO.write(isolate_seqs, isolate_fasta_fh, "fasta")

    # check isolates are in metadata and tree
    augur_leaf_names = set([strain for strain in augur_tree.iter_leaf_names()])
    for isolate in isolates:
        # checks full metadata table to see where the missing isolate
        # has gone missing
        if isolate not in metadata.index:
            logging.error(f"{isolate} could not be found in {nextmeta} "
                          f"check strain name provided is correct")
            sys.exit(1)
        # check for augur subsampling
        if isolate not in augur_leaf_names:
            logging.error(f"{isolate} could not be found in tree "
                          f"{augur_tree_fp} check it wasn't removed in "
                           "subsampling by augur modify build or edit "
                           "include file")
            sys.exit(1)

    return nextfasta, isolate_fasta, metadata, augur_tree


def get_closest_older_leaves_in_tree(augur_tree, metadata, isolates):
    """
    For each isolate find all nearest sequences in the tree with older
    collection dates than the input sample
    """

    # get just metadata for strains in tree
    augur_leaf_names = set([strain for strain in augur_tree.iter_leaf_names()])
    augur_metadata = metadata.loc[augur_leaf_names]
    augur_tree_distances = []

    isolate_pbar = tqdm.tqdm(isolates)
    for isolate in isolate_pbar:
        isolate_pbar.set_description(f"Searching tree for {isolate}")
        # checks full metadata table to see where the missing isolate
        # has gone missing
        if isolate not in metadata.index:
            logging.error(f"{isolate} could not be found in metadata "
                          f"check strain name provided is correct")
            sys.exit(1)
        # check for augur subsampling
        if isolate not in augur_leaf_names:
            logging.error(f"{isolate} could not be found in tree "
                          "(augur_tree) check it wasn't removed in "
                           "subsampling by augur modify build or edit "
                           "include file")
            sys.exit(1)

        #get metadata for all strains in tree older than query isolate
        strain_date = augur_metadata.loc[isolate, 'date']
        metadata_older = augur_metadata[augur_metadata['date'] < strain_date]
        older_strains = set(metadata_older.index)

        # for each leaf in the tree if its an older node get the distance
        # to the isolate node
        isolate_node = augur_tree.search_nodes(name=isolate)[0]
        isolate_tree_distances = {'query_strain': [],
                                  'isolate': [], 'distance': []}

        for leaf in augur_tree.iter_leaves():
            if leaf.name in older_strains and leaf.name != isolate:
                isolate_tree_distances['query_strain'].append(isolate)
                isolate_tree_distances['isolate'].append(leaf.name)
                distance = augur_tree.get_distance(isolate_node, leaf)
                isolate_tree_distances['distance'].append(distance)

        if len(isolate_tree_distances['distance']) == 0:
            logging.error(f"Could not find older strain in tree to {isolate}")
            sys.exit(1)
        else:
            augur_tree_distances.append(pd.DataFrame(isolate_tree_distances))

    augur_tree_distances = pd.concat(augur_tree_distances)

    # get closest relatives for each isolate
    closest_older_strains = augur_tree_distances.groupby('query_strain')['distance'].idxmin()
    closest_older_strains = augur_tree_distances.loc[closest_older_strains]

    #print(closest_older_strains)

    return closest_older_strains


def isolate_source_detector(isolates, nextfasta, nextmeta, augur_tree_fp,
                            output_dir, overwrite):
    """
    Main runner for ISD: finds closest relatives of isolates
    in both refined augur phylogeny and using minimap
    """

    run_name = output_dir + str(int(time.time()))

    logging.basicConfig(format='%(levelname)s:%(message)s',
                        level=logging.DEBUG,
                        handlers=[logging.FileHandler(f"{run_name}.log"),
                                  logging.StreamHandler()])

    # set-up output dir
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir()
    if output_dir.exists() and overwrite:
        shutil.rmtree(output_dir)
        output_dir.mkdir()
    else:
        logging.error(f"{output_dir} already exists, use alternative "
                       " or --overwrite")
        sys.exit(1)

    # get all metadata and data structure set up
    logging.info(f"Parsing input {nextfasta}, {nextmeta}, {augur_tree_fp}")
    nextfasta, isolate_fasta, \
            metadata, augur_tree = parse_input_files(nextfasta, nextmeta,
                                                     augur_tree_fp, isolates,
                                                     output_dir)


    # for each isolate find all nearest sequences in the tree with older
    # collection dates than the input sample
    logging.info("Searching for older closest relatives of isolates in "
                 f"{augur_tree_fp} ")

    closest_older_strains = get_closest_older_leaves_in_tree(augur_tree,
                                                             metadata,
                                                             isolates)

    # use minimap2 to find the closest sequence that is older to the input
    # grab strains from fasta
    # minimap2 of strains vs

    ## get all vs all node distances matrix from tree
    ## inefficient but simplifies later code greatly
    #augur_tree_distances = {'strain1': [], 'strain2': [], 'distance': []}
    #for leaf1 in augur_tree.iter_leaves():
    #    for leaf2 in augur_tree.iter_leaves():
    #        if leaf1 != leaf2:
    #            augur_tree_distances['strain1'].append(leaf1.name)
    #            augur_tree_distances['strain2'].append(leaf2.name)
    #            distance = augur_tree.get_distance(leaf1, leaf2)
    #            augur_tree_distances['distance'].append(distance)


