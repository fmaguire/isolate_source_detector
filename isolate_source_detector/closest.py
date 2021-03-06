import logging
import tqdm
import subprocess
import sys
import pandas as pd
import pickle
from io import StringIO
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


def older_mash_for_isolate(input_data):
    """
    Inner loop for parallelising mash across isolates without ending up with
    a 60GB mash dist file containing all hits
    """
    isolate, isolate_fasta, metadata, ref_sketch = input_data

    result = subprocess.run(f"mash dist -p 1 {ref_sketch} {isolate_fasta} 2> /dev/null",
                             stdout=subprocess.PIPE,
                             shell=True,
                             check=True)

    isolate_hits = pd.read_csv(StringIO(result.stdout.decode('utf8')),
                              sep='\t', names=["closest_ancestor",
                                               "isolate",
                                               "distance",
                                               "p-value",
                                               "shared-hashes"])
    # filter out self-hits
    isolate_hits = isolate_hits[isolate_hits['isolate'] != isolate_hits['closest_ancestor']]
    isolate_hits['metric'] = 'mash'
    # overwrite isolate name to replace fasta path with just isolate name
    isolate_hits['isolate'] = isolate

    isolate_date = metadata.loc[isolate, 'date']
    older_strains = set(metadata[metadata['date'] < isolate_date].index)
    older_isolate_hits = isolate_hits[isolate_hits['closest_ancestor'].isin(older_strains)]

    top_hit_ix = older_isolate_hits['distance'].nsmallest(1, keep='all').index
    top_hits = older_isolate_hits.loc[top_hit_ix, ['isolate', 'closest_ancestor',
                                                   'metric', 'distance']]
    return top_hits


def get_closest_older_genomes(isolate_fasta_index, ref_sketch, metadata,
                              output_dir, num_processes):
    """
    Use mash to get closest older genomes to isolate genomes via minimap2
    """
    # parallelise
    pool = Pool(num_processes)
    parallel_input = [(isolate, isolate_path, metadata, ref_sketch) \
                        for isolate, isolate_path in isolate_fasta_index.items()]
    mash_distances = pool.map(older_mash_for_isolate, parallel_input)
    pool.close()

    closest_older_mash_hits = pd.concat(mash_distances)

    return closest_older_mash_hits


def get_ancestor_traits(present_isolates, tree, traits):
    """
    Extract ancestor traits from augur's inferred traits json dictionary
    """
    data = {'isolate': [], 'trait_type': [],
            'trait': []}
    for isolate in present_isolates:
        node = tree.search_nodes(name=isolate)[0]
        node = node.up
        ancestor_traits = traits[node.name]

	for trait in ['region', 'country', 'division']:
	    data['isolate'].append(isolate)
	    data['trait_type'].append(trait)
	    data['trait_value'].append(ancestor_traits[trait])
	    # add exposure trait too
	    data['isolate'].append(isolate)
	    data['trait_type'].append(trait + "_exposure")
	    data['trait_value'].append(ancestor_traits[trait + "_exposure"])

    return pd.DataFrame(data)
