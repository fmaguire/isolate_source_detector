import logging
import tqdm
import subprocess
import sys
import pandas as pd
import pickle
from multiprocessing import Pool

def older_leaves_for_isolate(input_data):
    """
    Inner loop worker for parallelising extract of older leaves for a specific
    isolate
    """
    # extract input data
    isolate, tree_metadata, tree = input_data

    #get metadata for all strains in tree older than query isolate
    strain_date = tree_metadata.loc[isolate, 'date']
    metadata_older = tree_metadata[tree_metadata['date'] < strain_date]
    older_strains = set(metadata_older.index)

    # for each leaf in the tree if its an older node get the distance
    # to the isolate node
    isolate_node = tree.search_nodes(name=isolate)[0]
    isolate_tree_distances = {'isolate': [], 'closest_ancestor': [],
                              'metric': [], 'distance': []}
    for leaf in tree.iter_leaves():
        if leaf.name in older_strains and leaf.name != isolate:
            isolate_tree_distances['isolate'].append(isolate)
            isolate_tree_distances['closest_ancestor'].append(leaf.name)
            distance = tree.get_distance(isolate_node, leaf)
            isolate_tree_distances['metric'].append('phylo_distance')
            isolate_tree_distances['distance'].append(distance)

    if len(isolate_tree_distances['distance']) == 0:
        logging.error(f"Could not find older strain in tree to {isolate}")
        sys.exit(1)
    else:
        return pd.DataFrame(isolate_tree_distances)


def get_closest_older_leaves_in_tree(isolates, tree, metadata, num_processes):
    """
    For each isolate find all nearest sequences in the tree with older
    collection dates than the input sample
    """
    # get just metadata for strains in tree
    leaf_names = set([strain for strain in tree.iter_leaf_names() \
                        if not strain.startswith('unknown')])

    # i.e. for augur tree "strain"
    tree_metadata = metadata.loc[leaf_names]

    # parallelise
    pool = Pool(num_processes)
    parallel_input = [(isolate, tree_metadata, tree) for isolate in isolates]
    tree_distances = pool.map(older_leaves_for_isolate, parallel_input)
    pool.close()

    tree_distances = pd.concat(tree_distances)

    # get closest relatives for each isolate
    # can't use idxmin because want to keep minimum ties
    min_ix = tree_distances.groupby('isolate')['distance']\
                    .nsmallest(1, keep='all').reset_index()['level_1'].values
    closest_older_strains = tree_distances.loc[min_ix]

    return closest_older_strains


def get_closest_older_genomes(isolates, isolate_fasta_fp, sketch, metadata,
                              output_dir, num_processes):
    """
    Use mash to get closest older genomes to isolate genomes via minimap2
    """
    mash_hits = output_dir / "isolate_mash_hits.tsv"
    subprocess.run(f"mash dist -p {num_processes} -i {isolate_fasta_fp} {sketch} > {mash_hits} "
                    " 2> /dev/null",
                    shell=True,
                    check=True)

    mash_hits = pd.read_csv(mash_hits, sep='\t', names=["isolate",
                                                        "closest_ancestor",
                                                        "distance",
                                                        "p-value",
                                                        "shared-hashes"])

    # add metric type
    mash_hits['metric'] = 'mash'

    # filter out self-hits
    mash_hits = mash_hits[mash_hits['isolate'] != mash_hits['closest_ancestor']]

    # get older strain names
    closest_older_mash_hits = []
    for isolate in isolates:
        isolate_date = metadata.loc[isolate, 'date']
        older_strains = set(metadata[metadata['date'] < isolate_date].index)
        isolate_hits = mash_hits.loc[mash_hits['isolate'] == isolate]
        older_isolate_hits = isolate_hits[isolate_hits['closest_ancestor'].isin(older_strains)]
        closest_older_mash_hits.append(older_isolate_hits)

    closest_older_mash_hits = pd.concat(closest_older_mash_hits)

    # can't use idxmin because want to keep minimum ties
    min_ix = closest_older_mash_hits.groupby('isolate')['distance']\
                    .nsmallest(1, keep='all').reset_index()['level_1'].values

    closest_older_mash_hits = closest_older_mash_hits.loc[min_ix,
                                                          ['isolate',
                                                           'closest_ancestor',
                                                           'metric',
                                                           'distance']]
    return closest_older_mash_hits


def get_ancestor_traits(present_isolates, tree, traits):
    """
    Extract ancestor traits from augur's inferred traits json dictionary
    """
    data = {'isolate': [], 'closest_ancestor': [],
            'geographic_scale': [], 'geographic_loc': [],
            'inference_confidence': []}

    for isolate in present_isolates:
        node = tree.search_nodes(name=isolate)[0]
        node = node.up
        ancestor_traits = traits[node.name]
        for geo_scale in ['region_confidence', 'country_confidence', 'division_confidence']:
            for location, confidence in ancestor_traits[geo_scale].items():
                data['isolate'].append(isolate)
                data['closest_ancestor'].append(node.name)
                data['geographic_scale'].append(geo_scale.split('_')[0])
                data['geographic_loc'].append(location)
                data['inference_confidence'].append(confidence)

    return pd.DataFrame(data)
