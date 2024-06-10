for i in "VanRaden_" "noia_" "Luke_"
do
	bash ../launch_GBLUP.end.parallele.Univariate.sh "Kinship_matrix_${i}$2.csv" "${i}$2" "$2_t_convert_genotype.csv" "$2_snp_positions.csv"
done
wait
