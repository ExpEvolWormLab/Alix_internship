## Set of script to analyze sequencing and mapping step

### First, run the command : 
      bash script_summary_reads_table.sh path
  path - Path to a folder which contains a fastp_reports and a filtered_bam repositories
        !! In fastp_reports : files have to be name as Line.*json
        !! in filtered_bam : files have to be name as Line_*bam
  ouput :
        summary_reads.tsv - Contains Line Reads Reads_after_fastp Reads_final Coverage

This command has been run for each different batchs (NYC, CeMee, Tom, CeMee_2)
path NYC - /mnt/data5/mallard/2024_nyuline (GEVPC05) _script need to be adapted for this one, fastp_reports[_01/_02]_
path CeMee - /mnt/data2/mallard/2023_Mapping_RILS/ (GEVPC10)
path CeMee_2 - /mnt/usb_FM/missing_SRR_NCBI (GEVPC10)
path Tom - /mnt/data4/mallard/2024_RILs_Tom/ (GEVPC05)

Each summary_table computed for theses paths can be found in Alix_internship/DataAnalyzing/Outputs

### Analyzing summary_reads tables
#### Plot_summary_reads_tables.R
Interactif script
Output :Lines_Percentage_final_read_X_Percentage_read_kept_after_fastp.pdf
        Lines_Percentage_final_read_X_Coverage.pdf
        Zoom_Lines_Percentage_final_read_X_Percentage_read_kept_after_fastp.pdf
        Zoom_Lines_Percentage_final_read_X_Coverage.pdf
