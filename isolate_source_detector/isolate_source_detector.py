"""Main module."""

import sys
import logging
from pathlib import Path

from isolate_source_detector import utils, closest

def isd(isolates_fp, metadata_fp, fasta_fp, tree_fp, traits_fp,
        output_dir, debug):
    """
    Main runner for ISD: finds closest relatives of isolates
    in both refined augur phylogeny and using minimap
    """
    output_dir = Path(output_dir)
    utils.set_up_logging(output_dir, debug)
    utils.check_dependencies()
    utils.set_up_output_dir(output_dir)

    # get all metadata and data structure set up
    logging.info(f"Parsing input isolates: {isolates_fp}")
    with open(isolates_fp) as fh:
        isolates = [line.strip() for line in fh]

    #
    #    Parse input files and create the appropriate data structures
    #    and files for tree analysis and minimap analysis of closest
    #    older relatives to query strain
    #
    logging.info(f"Parsing input metadata: {metadata_fp}")
    metadata = utils.parse_metadata(metadata_fp)

    logging.info(f"Parsing input traits: {traits_fp}")
    traits = utils.parse_traits(traits_fp)

    logging.info(f"Parsing input tree: {tree_fp}")
    tree = utils.parse_tree(tree_fp)

    logging.info(f"Parsing input fasta: {fasta_fp}")
    fasta_strains, isolate_fasta_fp = utils.parse_fasta(fasta_fp, isolates,
                                                        output_dir)

    logging.info("Checking input files for consistency")
    present_isolates = utils.check_input_consistency(isolates,
                                                     metadata,
                                                     traits,
                                                     tree,
                                                     fasta_strains)

    logging.info(f"Sketching input fasta: {fasta_fp}")
    sketch = utils.mash_sketch_input_fasta(fasta_fp, output_dir)

    # for each isolate find all nearest sequences in the tree with older
    # collection dates than the input sample using a phylogenetic distances
    # lens
    logging.info("Searching for older closest relatives of isolates in "
                 f"{tree_fp}")
    closest_older_in_tree = closest.get_closest_older_leaves_in_tree(present_isolates,
                                                                     tree,
                                                                     metadata)
    #closest_older_in_tree = utils.add_geo_location_from_metadata(closest_older_in_tree,
    #                                                             metadata)

    closest_older_in_tree.to_csv(output_dir / "phylo_distance.tsv", sep='\t')

    # for each isolate get parent node's inferred exposure information
    # inferred_ancestor_exposure = closest.get_aancestor_exposure(tree,


    # for each isolate find closest older genomes using mash
    logging.info(f"Searching {fasta_fp} using mash sketch {sketch} "
                 "for closest relatives to isolates")
    closest_older_in_mash = closest.get_closest_older_genomes(present_isolates,
                                                              isolate_fasta_fp,
                                                              sketch,
                                                              metadata,
                                                              output_dir)

    #closest_older_in_mash = utils.add_geo_location_from_metadata(closest_older_in_mash,
    #                                                             metadata)

    closest_older_in_mash.to_csv(output_dir / "mash.tsv", sep='\t')


    # for each isolate find the inferred trait inference for ancestor node
    # in the tree
    #logging.info(f"Searching for ancestor inferred trait in {tree_fp} using "
    #             f"{traits_fp} ")
    #ancestor_traits = closest.get_ancestor_traits(present_isolates,
    #                                              tree,
    #                                              metadata)

