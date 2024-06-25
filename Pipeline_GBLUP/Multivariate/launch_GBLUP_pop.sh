#!/bin/bash


## Script to launch GBLUP on severale population
## First argument tsv file : first column path of vcf file, second column name of the output

# Check if the file exists
if [ ! -f $1 ]; then
	echo "File not found."
	exit 1
fi

# Read the file line by line
while IFS= read -r line; do
	echo "$line"
	bash launch_GBLUP.sh $line &
	 # Process each line as needed
done < $1

wait
