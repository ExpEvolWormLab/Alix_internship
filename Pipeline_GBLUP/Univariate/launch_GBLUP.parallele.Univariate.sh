#!bin/bash

# Script which launch GBLUP steps for each kinship matrix
# First argument kinship matrix
# Second argument : output name (extension of file not a new directory)
# Third argument : M matrix
# Fourth argument : snp_position

# Record the start time
start_time=$(date +%s)



Rscript 2-Invert_GRM.Univariate.R $1 $2
echo "$2: 2-Invert_GRM executed successfully."

for i in "NGM"
do
	for trait in "SF" "SB" "FS" "FB" "BS" "BF"
	do
		Rscript 3-MCMCglmm_model.Univariate.R "Inverted_kinship_matrix_$2.csv" $3 $i $trait $2_${i}_$trait
		echo "$2_${i}: MCMCglmm_model executed"

		if [ $? -eq 0 ]; then #Different sessions
	                Rscript 4-Diagnostic.Univariate.R "$2_${i}_${trait}_MCMCmodel_Sol.csv" "$2_${i}_${trait}_MCMCmodel_VCV.csv" 'session_id' $trait $2_${i}_$trait
        	        echo "$2_${i}: Diagnostic executed"
        	else
                	Rscript 4-Diagnostic.Univariate.R "$2_${i}_${trait}_MCMCmodel_Sol.csv" "$2_${i}_${trait}_MCMCmodel_VCV.csv" '' $trait $2_${i}_$trait
	                echo "$2_${i}: Diagnostic executed"
       		fi

		Rscript 5-Backsolving.Univariate.R "$2_${i}_${trait}_MCMCmodel_Sol.csv" $3 $4 $trait $2_${i}_$trait
		echo "$2_${i}: Backsolving executed"
	done
done &

for i in "NaCl" 
do
	for trait in "SF" "SB" "FS" "FB" "BS" "BF"
        do
                Rscript 3-MCMCglmm_model.Univariate.R "Inverted_kinship_matrix_$2.csv" $3 $i $trait $2_${i}_$trait
                echo "$2_${i}: MCMCglmm_model executed"

		if [ $? -eq 0 ]; then #Different sessions
                        Rscript 4-Diagnostic.Univariate.R "$2_${i}_${trait}_MCMCmodel_Sol.csv" "$2_${i}_${trait}_MCMCmodel_VCV.csv" 'session_id' $trait $2_${i}_$trait
                        echo "$2_${i}: Diagnostic executed"
                else
                        Rscript 4-Diagnostic.Univariate.R "$2_${i}_${trait}_MCMCmodel_Sol.csv" "$2_${i}_${trait}_MCMCmodel_VCV.csv" '' $trait $2_${i}_$trait
                        echo "$2_${i}: Diagnostic executed"
                fi

                Rscript 5-Backsolving.Univariate.R "$2_${i}_${trait}_MCMCmodel_Sol.csv" $3 $4 $trait $2_${i}_$trait
                echo "$2_${i}: Backsolving executed"
        done
done &

wait


# Calculate the elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "GBLUP $2 executed successfully in $elapsed_time seconds."
