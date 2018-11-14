#! /bin/bash

# Define input files
FILTERED_TAXONOMY="filtered.fixrank_OTU.fas_classified.txt"
OTU_COUNTS_TABLE="OTU_counts.table.txt"

# Define Taxonomic Level to summarise counts at
# Choices are: Domain Phylum Class Order Family Genus - NOTE: make sure to mactch spelling exactly!
TAXONOMIC_LEVEL="Genus"




# Find the Column that has the Taxonomic level defined
COLUMN_ID=`head -1 $FILTERED_TAXONOMY | awk -v TAXONOMIC_LEVEL=$TAXONOMIC_LEVEL '{ for (i=1; i<=NF; i++ ) { if ( $i == TAXONOMIC_LEVEL ) print i }}'`

# Check to see that the taxonomic ID specified is valid
if [ "$COLUMN_ID" == "" ]; then 
	echo "Error: The taxonomic level defined is not valid, it should be one of:"
	echo "Domain Phylum Class Order Family Genus"
	echo "You entered \"$TAXONOMIC_LEVEL\" - please check the definition on line 9 of this script"
else
	echo "Found $TAXONOMIC_LEVEL in column $COLUMN_ID in the file $FILTERED_TAXONOMY"


	# filter the OTU assignments to those OTUs that have an assignment at the taxonomic level chosen
	rm -rf dropped_OTUs.txt;
	cut -d' ' -f1,$COLUMN_ID < filtered.fixrank_OTU.fas_classified.txt | grep -v "#" | sort -k2,2 | awk '{ if( NF < 2 ) { print $1 >> "dropped_OTUs.txt" } else { print $0 }}' > temp_taxfile.txt




	# Summarise the OTU counts for each sample in the file OTU_COUNTS_TABLE (see line 2) at the taxonomic level defined 
	head -1 $OTU_COUNTS_TABLE > $TAXONOMIC_LEVEL.$OTU_COUNTS_TABLE
	
	sed '1d' $OTU_COUNTS_TABLE | awk 'BEGIN{
		while(( getline line < "temp_taxfile.txt" ) > 0 ) { 
		   split(line, a, " "); 
           assignment[a[1]]=a[2]; 
           }
	};{
		if( $1 in assignment )
			{
			printf assignment[$1];
			for(i=2; i<=NF; i++) printf "\t"$i;
			printf "\n";
			}
	}' |  sort -k1,1 | awk 'BEGIN{ PREV="";
	};{
		if($1 != PREV)
			{
			if(PREV != "")
				{
				printf PREV;
				for(j=2; j<=NF; j++) printf "\t"sum[j];
				printf "\n";
				}
			PREV=$1
			for(j=2; j<=NF; j++) sum[j]=$j;
			}
		else
			{
			for(j=2; j<=NF; j++) sum[j]=sum[j]+$j;
			}
	}END{
		printf PREV;
		for(j=2; j<=NF; j++) printf "\t"sum[j];
		printf "\n";
	}' >> $TAXONOMIC_LEVEL.$OTU_COUNTS_TABLE
	echo "Summarised count file at $TAXONOMIC_LEVEL level has been written to file $TAXONOMIC_LEVEL.$OTU_COUNTS_TABLE"
fi

rm -rf temp_taxfile.txt

