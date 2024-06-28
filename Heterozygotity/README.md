## Set of script to get heterozygotity along chromosomes and to plot it

### 1-Get_hetero.R :
  Rscript 1-Get_hetero.R het_flag.vcf
Parameters: het_flag.vcf - vcf file when heterozygotes loci where flag as is_het (combined.annotated.filtred.vcf.gz)
            !! Some parameters as chromosomes to plot, population to plot or number of SNPs in a bin can be changed directly on the script
**Output** : Heterozygosity*csv - Mean of lines' mean heterozygotity for the windows of SNPs definite by the user along chromosomes

### 2-Plotting_hetero.R :
  To run in local
Hyper divergent file is necessary 
**Outputs** : Heterozygotity_HDregion - pdf heterozygotity along chromosomes with hyper divergent regions plotted as black dots
          Heterozygotity_HDregion_Zoom_V - pdf heterozygotity for the remarkable regions of chromosome V with hyper divergent regions plotted as black dots
          Randomisation - pdf comparison of observed proportion of heterozygotes associated with HDR for all choromosomes compared to the distribution obtained by chance 
          Randomisation_Chrom5 - pdf comparison of observed proportion of heterozygotes associated with HDR for the remarkable regions of chromosome V compared to the distribution obtained by chance 

### Plotting_Chrom5 :
  To run in local
Hyper divergent file is necessary
All pop should be in the working directory name as Heterozygosity_pop.csv 
**Output** : Chromosome_5_Allpop.pdf - pdf heterozygotity hromosomes with hyper divergent regions plotted as black dots for all pop for the remarkable region of chromosome V