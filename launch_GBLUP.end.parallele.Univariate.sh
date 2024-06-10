#!bin/bash

# Script which launch GBLUP steps for each kinship matrix
# First argument kinship matrix
# Second argument : output name (extension of file not a new directory)
# Third argument : M matrix
# Fourth argument : snp_position

# Record the start time
start_time=$(date +%s)




for i in "NGM"
do
	for trait in "SF" "SB" "FS" "FB" "BS" "BF"
	do
	        Rscript ../4-Diagnostic.Univariate.R "$2_${i}_${trait}_MCMCmodel_Sol.csv" "$2_${i}_${trait}_MCMCmodel_VCV.csv" 'session_id' $trait $2_${i}_$trait
        	echo "$2_${i}: Diagnostic executed"

		Rscript ../5-Backsolving.Univariate.R "$2_${i}_${trait}_MCMCmodel_Sol.csv" $3 $4 $trait $2_${i}_$trait
		echo "$2_${i}: Backsolving executed"
	done
done &

for i in "NaCl" 
do
	for trait in "SF" "SB" "FS" "FB" "BS" "BF"
        do
                Rscript ../4-Diagnostic.Univariate.R "$2_${i}_${trait}_MCMCmodel_Sol.csv" "$2_${i}_${trait}_MCMCmodel_VCV.csv" 'session_id' $trait $2_${i}_$trait
                echo "$2_${i}: Diagnostic executed"

                Rscript ../5-Backsolving.Univariate.R "$2_${i}_${trait}_MCMCmodel_Sol.csv" $3 $4 $trait $2_${i}_$trait
                echo "$2_${i}: Backsolving executed"
        done
done &

wait


# Calculate the elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "GBLUP $2 executed successfully in $elapsed_time seconds."
