#! /bin/bash

# Take the files from the folder "raw_seqs" and rename them by adding the filename to the name of ervey sequence and then combine all renamed sequences into a file called "input_sequences.fas" in the previous directory.
rm -rf input_sequences.fas
cd ./raw_seqs;
for i in *.fas; do
    name=`echo $i | sed 's/.fas//g'`;
    sed "s/^>/>${name}_/g" < $i >> ../input_sequences.fas;
done


