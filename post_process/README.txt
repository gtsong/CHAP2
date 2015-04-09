
expand_branches:
  This pipeline component expands the tree_branch field in the
  all.gc file to a comma-separated list of tree edges (sub-lineage).
  The original all.gc file has only a single, upper-bound edge,
  but this program examines the detection patterns in the primary
  and outgroup species to identify lower bounds as well.  If it
  finds conflicting information, it prints a question mark '?'
  for the lower bound.

reformat_nonredundant:
  This pipeline component converts the all.remove_redundancy.gc
  file to the same format as all.gc, by extracting entire lines from
  all.gc that correspond to the ones in all.remove_redundancy.gc.

gc-info:
  This utility is run by the user as needed, rather than by the
  pipeline.  It prints a summary of the detected conversions
  in each species.

