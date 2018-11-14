#! /bin/bash

# Run CD-Hit-OTU to pick OTUs
cd-hit-otu-all.pl  -i input_sequences.fas -c 0.97 -o Cluster_results | tee ./Cluster_results/clustering_output.txt

# Rename representative sequences to OTU number
cut -d' ' -f1  ./Cluster_results/OTU | sed '/>/s/$/</' | tr -d '\n' |tr '>' '\n' | grep . | awk '{ n+=1; print ">OTU"n" "$0 }' | tr '<' '\n' | cut -d' ' -f1 > OTU.fas


