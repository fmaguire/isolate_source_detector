import logging
import tqdm
import subprocess
import sys
import pandas as pd


def get_closest_older_leaves_in_tree(tree, metadata, isolates):
    """
    For each isolate find all nearest sequences in the tree with older
    collection dates than the input sample
    """
    # get just metadata for strains in tree
    leaf_names = set([strain for strain in tree.iter_leaf_names() \
                        if not strain.startswith('unknown')])

    # i.e. for augur tree "strain"
    tree_metadata = metadata.loc[leaf_names]
    tree_distances = []

    isolate_pbar = tqdm.tqdm(isolates)
    for isolate in isolate_pbar:
        isolate_pbar.set_description(f"Searching tree for {isolate}")
        # checks full metadata table to see where the missing isolate
        # has gone missing
        if isolate not in tree_metadata.index:
            logging.error(f"{isolate} could not be found in metadata "
                          f"check strain name provided is correct")
            sys.exit(1)
        # check for subsampling
        if isolate not in leaf_names:
            logging.error(f"{isolate} could not be found in tree "
                          "(tree) check it wasn't removed in "
                           "subsampling by augur modify build or edit "
                           "include file")
            sys.exit(1)

        #get metadata for all strains in tree older than query isolate
        strain_date = tree_metadata.loc[isolate, 'date']
        metadata_older = tree_metadata[tree_metadata['date'] < strain_date]
        older_strains = set(metadata_older.index)

        # for each leaf in the tree if its an older node get the distance
        # to the isolate node
        isolate_node = tree.search_nodes(name=isolate)[0]
        isolate_tree_distances = {'isolate': [],
                                  'tree_genome': [], 'phylo_distance': []}

        for leaf in tree.iter_leaves():
            if leaf.name in older_strains and leaf.name != isolate:
                isolate_tree_distances['isolate'].append(isolate)
                isolate_tree_distances['tree_genome'].append(leaf.name)
                distance = tree.get_distance(isolate_node, leaf)
                isolate_tree_distances['phylo_distance'].append(distance)

        if len(isolate_tree_distances['phylo_distance']) == 0:
            logging.error(f"Could not find older strain in tree to {isolate}")
            sys.exit(1)
        else:
            tree_distances.append(pd.DataFrame(isolate_tree_distances))

    tree_distances = pd.concat(tree_distances)

    # get closest relatives for each isolate
    # can't use idxmin because want to keep minimum ties

    min_ix = tree_distances.groupby('isolate')['phylo_distance']\
                    .nsmallest(1, keep='all').reset_index()['level_1'].values
    closest_older_strains = tree_distances.loc[min_ix]

    return closest_older_strains


def get_closest_older_genomes(sketch, isolate_fasta_fp, isolates, metadata,
                              output_dir):
    """
    Use mash to get closest older genomes to isolate genomes via minimap2
    """

    mash_hits = output_dir / "isolate_mash_hits.tsv"
    subprocess.run(f"mash dist -i {isolate_fasta_fp} {sketch} > {mash_hits} "
                    " 2> /dev/null",
                    shell=True,
                    check=True)

    mash_hits = pd.read_csv(mash_hits, sep='\t', names=["isolate",
                                                        "mash_hit_genome",
                                                        "mash_distance",
                                                        "p-value",
                                                        "shared-hashes"])

    # filter out self-hits
    mash_hits = mash_hits[mash_hits['isolate'] != mash_hits['mash_hit_genome']]

    # get older strain names
    closest_older_mash_hits = []
    for isolate in isolates:
        isolate_date = metadata.loc[isolate, 'date']
        older_strains = set(metadata[metadata['date'] < isolate_date].index)
        isolate_hits = mash_hits.loc[mash_hits['isolate'] == isolate]
        older_isolate_hits = isolate_hits[isolate_hits['mash_hit_genome'].isin(older_strains)]
        closest_older_mash_hits.append(older_isolate_hits)

    closest_older_mash_hits = pd.concat(closest_older_mash_hits)

    # can't use idxmin because want to keep minimum ties
    min_ix = closest_older_mash_hits.groupby('isolate')['mash_distance']\
                    .nsmallest(1, keep='all').reset_index()['level_1'].values

    closest_older_mash_hits = closest_older_mash_hits.loc[min_ix,
                                                          ['isolate',
                                                           'mash_hit_genome',
                                                           'mash_distance']]

    return closest_older_mash_hits

