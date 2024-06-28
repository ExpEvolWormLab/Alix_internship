# Set of Scripts to Analyze Sequencing and Mapping Step

This document provides an overview of how to perform the analysis of sequencing and mapping steps using the provided scripts and tools. The steps include running the initial summarization script, followed by analysis using the interactive R script to generate informative plots.

## First Step
Run the following command:
```bash
bash script_summary_reads_table.sh path
```
- **path**: Path to a folder which contains `fastp_reports` and `filtered_bam` directories.
  - **Note:**
    - In `fastp_reports`, files must be named as `Line.*json`.
    - In `filtered_bam`, files must be named as `Line_*bam`.

### Output
- `summary_reads.tsv` - Contains `Line`, `Reads`, `Reads_after_fastp`, `Reads_final`, `Coverage`.

This command has been run for each different batch (NYC, CeMee, Tom, CeMee_2):
- **path NYC**: `/mnt/data5/mallard/2024_nyuline (GEVPC05)` (_script needs to be adapted for this one, `fastp_reports[_01/_02]`_)
- **path CeMee**: `/mnt/data2/mallard/2023_Mapping_RILS/ (GEVPC10)`
- **path CeMee_2**: `/mnt/usb_FM/missing_SRR_NCBI (GEVPC10)`
- **path Tom**: `/mnt/data4/mallard/2024_RILs_Tom/ (GEVPC05)`

Each summary table computed for these paths can be found in `Alix_internship/DataAnalyzing/Outputs`.

## Analyzing `summary_reads` Tables
### `Plot_summary_reads_tables.R`
Interactive script

### Outputs
- `Lines_Percentage_final_read_X_Percentage_read_kept_after_fastp.pdf`
- `Lines_Percentage_final_read_X_Coverage.pdf`
- `Zoom_Lines_Percentage_final_read_X_Percentage_read_kept_after_fastp.pdf`
- `Zoom_Lines_Percentage_final_read_X_Coverage.pdf`
