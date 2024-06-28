# Pipeline to Launch GATK Pipeline to Call Variants Using Tools from GATK v.4.5.0.0


This document provides an overview of how to launch the GATK pipeline to call variants using GATK v.4.5.0.0 tools. The steps include pre-processing, variant calling, genotyping, and filtration, with detailed descriptions of the commands, scripts, parameters, and outputs for each stage.

## 1. Pre-processing and Variant Calling

### Command:
```bash
bash launch_GATK.part1.sh folder output_directory
```

### Script: `launch_GATK.part1.sh`
#### Description:
This script performs pre-processing and variant calling steps using GATK tools.

#### Steps:
- **Pre-processing:**
  - `MarkDuplicates`
  - `AddOrReplaceReadGroups`
  - `BaseRecalibrator`
  - `ApplyBQSR`
- **Variant Calling:**
  - `HaplotypeCaller`

#### Parameters:
- **folder**: Folder where `.bam` files are stored.
- **output_directory**: Name of the output directory (does not need to exist).

#### Outputs:
- `output_directory/dedub_files/*.sorted.dedub.bam`: Files with duplicated reads flagged.
- `output_directory/dedub_files/*.sorted.dedub.metrics.txt`: Duplication metrics files.
- `output_directory/dedub_files/*.sorted.dedub.RG.bam`: Files with read groups added.
- `output_directory/BQSR_files/*_recal_data.table`: Recalibration tables.
- `output_directory/BQSR_files/*.sorted.dedub.corrected.bam`: BAM files with base quality scores corrected.
- `output_directory/HC_files/*.g.vcf.gz`: GVCF files containing all variants called.

## 2. Create a List of GVCF Files

### Command:
```bash
find -name *reheader.g.vcf.gz > .list | sed "s/\.\//$pwd/g"
```

#### Description:
This command finds all GVCF files and creates a list of their paths in `.list`. This list is required before running `launch_GATK.part2.sh`.

#### Output:
- `.list`: Contains paths of all GVCF files to be combined.

## 3. Genotyping and Filtration

### Command:
```bash
bash launch_GATK.part2.sh .list output_directory
```

### Script: `launch_GATK.part2.sh`
#### Description:
This script performs genotyping and filtration steps using GATK and bcftools.

#### Steps:
- **Genotyping:**
  - `CombineGVCFs`
  - `GenotypeGVCFs`
- **Filtration:**
  - `SelectVariants`
  - `VariantFiltration`
  - `bcftools` (view, filter)

#### Parameters:
- **.list**: Contains paths of all GVCF files to be combined.
- **output_directory**: Name of the output directory (does not need to exist).

#### Outputs:
- `output_directory/All_combined.${chrom}.g.vcf.gz`: Combined GVCF for each chromosome.
- `output_directory/All_combined.final.${chrom}.vcf.gz`: Genotyped combined GVCF for each chromosome.
- `output_directory/All_Lines_Chrom.final.vcf.gz`: VCF file containing all chromosomes.
- `output_directory/Filtred/combined.SNP.annotated.vcf.gz`: VCF file containing only SNPs.
- `output_directory/Filtred/combined.SNP.annotated.1.vcf.gz`: VCF file with monomorphic SNPs removed.
- `output_directory/Filtred/combined.annotated.filtred.vcf.gz`: Soft-filtered VCF with caeNDR parameters.
- `output_directory/Filtred/combined.annotated.filtred.1.vcf.gz`: Soft-filtered VCF with additional filter for high heterozygosity (>85%).
- `output_directory/Filtred/combined.without_hetero.vcf.gz`: VCF with heterozygote loci set as missing.
- `output_directory/Filtred/combined.final.vcf.gz`: Soft-filtered VCF with additional filter for high missing values (≥ 95%).
- `output_directory/Filtred/SNP_2exclude.txt`: File containing SNPs flagged as NON PASS, singletons, multi-allelic (to exclude from imputation).
- `output_directory/Filtred/summary.txt`: Summary of how many SNPs need to be excluded in each category (a SNP can be in several categories).
- `output_directory/Filtred/combined.final.hard_filter.vcf.gz`: Hard-filtered VCF with all SNPs in `SNP_2exclude.txt` removed.

