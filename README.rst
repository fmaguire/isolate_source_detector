=======================
isolate_source_detector
=======================

Script to extract likely isolate sources from phylogenies, augur trait inference,
and nextstrain metadata.

For a set of supplied isolate names:

- Extracts the location of the nearest older sequence in the tree 
  (by phylogenetic distance). As ideally this is a time-corrected tree we 
  don't have to deal with complications like a long branching sister leaf 
  (note: Rob Lanfear's global phylogeny is not time-corrected)

- Extracts the inferred ancestor location trait from the results of 
  `augur trait`.

- Extracts the location of the nearest older genome to the isolate as inferred
  by mash distance.


