=======================
isolate_source_detector
=======================

Script to extract likely isolate sources from augur time-corrected phylogenies, 
augur trait inference, genomic fasta sequences and nextstrain metadata.

For a set of supplied isolate names:

- Extracts the location of the nearest older sequence in the tree 
  (by phylogenetic distance). As this is a time-corrected tree we 
  don't have to deal with complications like a long branching sister leaf 

- Extracts the inferred ancestor location trait from the results of 
  `augur trait`.

- Extracts the location of the nearest older genome to the isolate as inferred
  by mash distance.

