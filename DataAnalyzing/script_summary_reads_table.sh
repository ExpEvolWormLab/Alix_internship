#!/bin/bash
## Script to summarize information about read for lines in a specific folder
## folder has to contains fastp_reports folder and a filtred_bam folder
## File names as to be specified in the fastp_reports folder as Name.*.json
## File names as to be specified in the filtered_bam folder as Name_*.bam
#Path is given as $1

echo -e "Line\tReads\tReads_after_fastp\tReads_final\tCoverage" > summary_reads.tsv

#In report_fastp repertory
for file in $(ls $1/fastp_reports/*json)
do
	#Keep the name of the line
	line=${file%.*}
	line=${line##*/}
	#Convert the name
	if grep -rq "\b$line\b" correspondances.txt
	then
		line=$(grep -r "\b$line\b" correspondances.txt|awk '{print $2}')
	fi
	before=$(jq ".summary.before_filtering.total_reads" $file)
	after=$(jq ".summary.after_filtering.total_reads" $file)
	echo -e $line"\t"$before"\t"$after >> summary_reads.tsv
	#echo -e $line"\t"$before"\t"$after
done

#In filtered_bam repertory
for file in $(ls $1/filtered_bam/*bam)
do
	#Keep the name of the line
        line=${file%_*}
        line=${line##*/}
	reads=$(samtools view -c $file)
	coverage=$(samtools depth $file | awk '{sum+=$3} END { print sum/NR}')
	#Add at the good line
	sed -i  "/\b$line\b/s/$/\t$reads\t$coverage/" summary_reads.tsv
done
