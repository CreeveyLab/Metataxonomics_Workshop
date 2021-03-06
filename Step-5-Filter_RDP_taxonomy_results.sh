#! /bin/bash

# Set cut off for filtering of taxonomy assignments from RDP
PID_CUTOFF=80;
TAXONOMY_FILE="fixrank_OTU.fas_classified.txt"

# Write Header of column names to output file first.
echo "#OTUID Domain Phylum Class Order Family Genus" > filtered.$TAXONOMY_FILE

# Run filter of taxonomy results from RDP
grep "^OTU[0-9]*;" $TAXONOMY_FILE | tr ';' ' ' | tr -d '{%"}' | awk -v PID_CUTOFF=$PID_CUTOFF '{ printf $1" "; taxa=3; PID=4;  while( taxa < NF ) { if( $PID >= PID_CUTOFF ) { printf $taxa" "} taxa=taxa+2; PID=PID+2 } printf "\n" }' >> filtered.$TAXONOMY_FILE

# This reults in a filtered taxonomy file with the word "filtered." appended to the front of the input taxonomy file

