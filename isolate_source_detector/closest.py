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
        isolate_tree_distances = {'query_strain': [],
                                  'isolate': [], 'distance': []}

        for leaf in tree.iter_leaves():
            if leaf.name in older_strains and leaf.name != isolate:
                isolate_tree_distances['query_strain'].append(isolate)
                isolate_tree_distances['isolate'].append(leaf.name)
                distance = tree.get_distance(isolate_node, leaf)
                isolate_tree_distances['distance'].append(distance)

        if len(isolate_tree_distances['distance']) == 0:
            logging.error(f"Could not find older strain in tree to {isolate}")
            sys.exit(1)
        else:
            tree_distances.append(pd.DataFrame(isolate_tree_distances))

    tree_distances = pd.concat(tree_distances)

    # get closest relatives for each isolate
    closest_older_strains = tree_distances.groupby('query_strain')['distance'].idxmin()
    closest_older_strains = tree_distances.loc[closest_older_strains]

    return closest_older_strains


def get_closest_older_genomes(sketch, isolate_fasta_fp, metadata,
                              output_dir):
    """
    Use mash to get closest older genomes to isolate genomes via minimap2
    """

    mash_hits = output_dir / "isolate_mash_hits.tsv"
    subprocess.run(f"mash dist -i {isolate_fasta_fp} {sketch} > {mash_hits}",
                    shell=True,
                    check=True)

    mash_hits = pd.read_csv(mash_hits, sep='\t', names=["mash_hit_genome",
                                                        "isolate",
                                                        "mash_distance",
                                                        "p-value",
                                                        "shared-hashes"])

    mash_hits = pd.melt(mash_hits, id_vars='#query',
                        var_name='mash_hit_genome', value_name='mash_dist')

    # filter out self-hits
    mash_hits = mash_hits[mash_hits['mash_hit_genome'] != mash_hits['isolate']]

    # get older strain names
    closest_older_mash_hits = []
    for isolate in isolates:
        isolate_date = metadata.loc[isolate, 'date']
        older_strains = metadata.loc[metadata['date'] < isolate_date, 'strain']

        isolate_hits = mash_hits.loc[mash_hits['isolate'] == isolate]
        older_isolate_hits = isolate_hits.loc[isolate_hits['mash_hit_genome'].isin(older_strains)]
        closest_older_mash_hits.append(older_isolate_hits)

    closest_older_mash_hits = pd.concat(closest_older_mash_hits)

    closest_older_mash_hits = closest_older_mash_hits.groupby('query_strain')['distance'].idxmin()

    return closest_older_mash_hits

