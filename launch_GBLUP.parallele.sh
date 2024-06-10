#!bin/bash

# Script which launch GBLUP steps for each kinship matrix
# First argument kinship matrix
# Second argument : output name (extension of file not a new directory)
# Third argument : M matrix
# Fourth argument : snp_position

# Record the start time
start_time=$(date +%s)



Rscript 2-Invert_GRM.R $1 $2
echo "$2: 2-Invert_GRM executed successfully."

for i in "NGM"
do
	Rscript 3-MCMCglmm_model.R "Inverted_kinship_matrix_$2.csv" $3 $i $2_$i
	echo "$2_${i}: MCMCglmm_model executed"

	if [ $? -eq 0 ]; then #Different sessions
		Rscript 4-Diagnostic.R "$2_${i}_MCMCmodel_Sol.csv" "$2_${i}_MCMCmodel_VCV.csv" 'session_id' $2_$i
		echo "$2_${i}: Diagnostic executed"
	else
		Rscript 4-Diagnostic.R "$2_${i}_MCMCmodel_Sol.csv" "$2_${i}_MCMCmodel_VCV.csv" '' $2_$i
                echo "$2_${i}: Diagnostic executed"
	fi

	Rscript 5-Backsolving.R "$2_${i}_MCMCmodel_Sol.csv" $3 $4 $2_$i
	echo "$2_${i}: Backsolving executed"
done &

for i in "NaCl" 
do
        Rscript 3-MCMCglmm_model.R "Inverted_kinship_matrix_$2.csv" $3 $i $2_$i
        echo "$2_${i}: MCMCglmm_model executed"

        if [ $? -eq 0 ]; then #Different sessions
                Rscript 4-Diagnostic.R "$2_${i}_MCMCmodel_Sol.csv" "$2_${i}_MCMCmodel_VCV.csv" 'session_id' $2_$i
                echo "$2_${i}: Diagnostic executed"
        else
                Rscript 4-Diagnostic.R "$2_${i}_MCMCmodel_Sol.csv" "$2_${i}_MCMCmodel_VCV.csv" '' $2_$i
                echo "$2_${i}: Diagnostic executed"
        fi

        Rscript 5-Backsolving.R "$2_${i}_MCMCmodel_Sol.csv" $3 $4 $2_$i
        echo "$2_${i}: Backsolving executed"
done &

wait


# Calculate the elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "GBLUP $2 executed successfully in $elapsed_time seconds."
