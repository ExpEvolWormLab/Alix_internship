# Scripts to Perform Imputation Per Population, Filter It, and Prune It

**Note:**
- Don't forget to modify `path_gatk`.
- `SNP_2exclude.txt` must be in the same directory.

**Packages:**
- Beagle5
- bcftools
- GATK v4.5.0
- plink
- vcftools
- bgzip

## Usage
Run the following command:
```bash
bash launch_imputation_pruning.sh final.vcf.gz output_directory
```

### `launch_imputation_pruning.sh`
This script performs:
1. **Imputation Step**:
   - Tools: `bcftools` (query, view), `beagle.22Jul22.46e.jar`
2. **Filtration Step**:
   - Tools: `bcftools` (view, filter), GATK (VariantFiltration, SelectVariants), `vcftools`
3. **Pruning Step**:
   - Tools: `plink`

## Parameters
- `final.vcf.gz`: Final VCF with lines filtered for high heterozygosity and isotype groups.
- `output_directory`: Directory where output files will be stored (does not need to already exist).

## Outputs
- `output_directory/.sample.name`: Names of all lines belonging to the population.
- `output_directory/_final.vcf.gz`: VCF files containing data for the specific population.
- `output_directory/imputed_${chrom}_${pop}.vcf.gz`: Imputed VCF files for each chromosome for each population.
- `output_directory/imputed_All_*.vcf.gz`: Imputed VCF files for all chromosomes for each population.
- `output_directory/pop/Filtered/imputed.SNP.vcf.gz`: VCF with only SNPs.
- `output_directory/pop/Filtered/imputed.SNP.het.vcf.gz`: Soft filter VCF for heterozygosity.
- `output_directory/pop/Filtered/imputed.SNP.het1.vcf.gz`: Soft filter VCF for high heterozygosity.
- `output_directory/pop/Filtered/imputed.SNP.without_hetero.vcf.gz`: VCF with heterozygote loci set as missing values.
- `output_directory/pop/Filtered/imputed.SNP.filtered.vcf.gz`: Soft filter VCF for high missing data.
- `output_directory/pop/Filtered/imputed.SNP.hard_filtered.vcf.gz`: Hard filter VCF, removing non-pass, multiallelic, monomorphic sites.
- `output_directory/pop/Filtered/imputed.SNP.filtered.final.vcf.gz`: Hard filter VCF, removing MAF < 0.05.
- `output_directory/pop/Pruned/prune.in`: SNPs not in linkage disequilibrium.
- `output_directory/pop/Pruned/prune.out`: SNPs in linkage disequilibrium.
- `output_directory/pop/Pruned/*.vcf.gz`: VCF files with only SNPs not in linkage disequilibrium.
