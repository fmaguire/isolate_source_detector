# Isolate Source Detector

Script to produce summaries of the likely source of a given isolate.
This is primarily based on augur (treetime)-corrected phylogenies trait inference
with back-up checks of localation metadata for closest older relatives using mash and within the phylogeny.

The mash-based closest genome check should ensure that subsampling in the augur phylogeny hasn't accidentally removed any key taxa.


Specifically, for a set of supplied isolate names this script:

- Extracts the location of the nearest older sequence in the tree 
  (by phylogenetic distance). As this is a time-corrected tree we 
  don't have to deal with complications like a long branching sister leaf 

- Extracts the inferred ancestor location trait from the results of 
  `augur trait`.

- Extracts the location of the nearest older genome to the isolate as inferred
  by mash distance.

## Installation

As an external dependency this only requires [mash](https://github.com/marbl/Mash) (although it does require the tree and traits file from an ncov/augur run with appropriate configuration e.g., [here](https://github.com/fmaguire/ncov).

This can be installed using conda (with bioconda channels appropriately configured):

        conda create -n isd mash pip

Other python dependencies should be installed automatically:

        git clone https://github.com/fmaguire/isolate_source_detector
        pip install isolate_source_detector

## Usage

        Usage: isd [OPTIONS]
        
          Console script for isolate_source_detector
        
        Options:
          -i, --isolates PATH          File containing list of isolate strain names
                                       (one per line)  [required]
        
          -m, --metadata PATH          Nextstrain metadata file  [required]
          -f, --fasta PATH             Nextstrain fasta file  [required]
          -t, --tree PATH              Augur refined time-corrected phylogeny
                                       [required]
        
          -j, --traits PATH            Augur inferred traits json e.g. traits.json
                                       [required]
        
          -o, --output_dir PATH        Output directory
          --debug                      Log debugging information
          --version
          -n, --num_processes INTEGER  Number of processes to execute ISD using.
          -h, --help                   Show this message and exit.


 e.g. using an augur run (with inference of region, country, and division enabled):
        
        isd --isolates list_of_isolates.txt --metadata ../ncov/data/metadata_2020-12-14_11-32.tsv \
                --fasta ../ncov/data/sequences_2020-12-14_08-02.fasta --tree ../ncov/results/north-america_canada/tree.nwk \
                --traits ../ncov/results/north-america_canada/traits.json 

## Output

Currently outputs 3 separate tsvs (plan is to tidy this up and create a bit more of a report).

Table containing the closest genome and its metadata location to

	out/mash.tsv - closest older genome by mash distance with location metadata 
	out/phylo_distance.tsv - closest older isolate by phylogenetic distance with location metadata
	out/inferred_traits.tsv - inferred location data for ancestral node for each isolate
