#! /bin/bash

# Run CD-HIT-OTU Command to create count table

clstr_sample_count_matrix.pl _ ./Cluster_results/OTU.nr2nd.clstr

# Copy the OTU count table to the work folder.
cat ./Cluster_results/OTU.nr2nd.clstr.otu.txt | awk '{ printf $1; for(i=2; i<NF; i++) printf "\t"$i; printf "\n"; }' > ./OTU_counts.table.txt


