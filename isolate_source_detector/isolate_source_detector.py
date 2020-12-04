"""Main module."""

import sys
import logging
from pathlib import Path

from isolate_source_detector import utils, closest

def isd(isolates, metadata_fp, fasta_fp, tree_fp, tree_label_type, output_dir,
        debug, use_existing):
    """
    Main runner for ISD: finds closest relatives of isolates
    in both refined augur phylogeny and using minimap
    """
    utils.set_up_logging(output_dir, debug)
    utils.check_dependencies()
    utils.set_up_output_dir(output_dir, use_existing)

    # get all metadata and data structure set up
    logging.info(f"Parsing input {metadata_fp}, {fasta_fp}, {tree_fp}")

    isolate_fasta_fp, sketch, metadata, tree = utils.parse_input_files(isolates,
                                                                metadata_fp,
                                                                fasta_fp,
                                                                tree_fp,
                                                                tree_label_type,
                                                                output_dir,
                                                                use_existing)

    # for each isolate find all nearest sequences in the tree with older
    # collection dates than the input sample
    #logging.info("Searching for older closest relatives of isolates in "
    #             f"{tree_fp}")
    #closest_older_in_tree = closest.get_closest_older_leaves_in_tree(tree,
    #                                                                 metadata,
    #                                                                 isolates)

    #print(closest_older_in_tree)
    # for each isolate find closest older genomes using mash
    logging.info(f"Searching {fasta_fp} using mash for closest relatives"
                  " to isolates")
    closest_older_in_mash = closest.get_closest_older_genomes(sketch,
                                                              isolate_fasta_fp,
                                                              metadata,
                                                              output_dir)

    print(closest_older_in_mash)
