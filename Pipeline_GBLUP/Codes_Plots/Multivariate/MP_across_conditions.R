### Script to compare credible SNP between conditions
## Output : Manhattan Plot - effet median a postefiori of SNP along chromosome colored by association in the two environments

setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Pruned")
output <- 'VanRaden_A6'
pruning <- '0.99'
# Charger les fichiers CSV
motif_NGM <- gsub('PPP',pruning,gsub('XXX',output,"Summary_[A-Z]*_XXX_NGM_PPP.csv"))
motif1_NGM <- gsub('PPP',pruning,gsub('XXX',output,"Summary_(.*)_XXX_NGM_PPP.csv"))
files <- list.files(pattern = motif_NGM)
file_outputs_MT <- sub(motif1_NGM, "\\1", files)
dfs_NGM <- lapply(files, read.csv)


motif_NaCl <- gsub('PPP',pruning,gsub('XXX',output,"Summary_[A-Z]*_XXX_NaCl_PPP.csv"))
motif1_NaCl <- gsub('PPP',pruning,gsub('XXX',output,"Summary_(.*)_XXX_NaCl_PPP.csv"))
files <- list.files(pattern = motif_NaCl)
file_outputs_ST <- sub(motif1_NaCl, "\\1", files)
dfs_NaCl <- lapply(files, read.csv)


# Combine the NGM and NaCl data frames
combined_dfs <- list()
for (i in 1:length(dfs_NGM)) {
  df_MT <- dfs_NGM[[i]]
  df_ST <- dfs_NaCl[[i]]
  
  df_MT$source <- "NGM"
  df_ST$source <- "NaCl"
  
  combined_df <- rbind(df_MT, df_ST)
  combined_dfs[[i]] <- combined_df
}


# Plotting
plots <- list()
plots_mp <- list()

for (i in 1:length(combined_dfs)) {
  df <- combined_dfs[[i]]
  df <- separate(df, X, into = c("CHROM", "POS"), sep = "_", remove = FALSE)
  df$POS <- as.numeric(df$POS)
  df$CHROM <- gsub('23','X',df$CHROM)
  
  trait <- file_outputs_MT[[i]]
  
  # Add new credible status for coloring
  df$credible_status <- "Not Credible"
  df$credible_status[df$credible == 'Credible' & df$source == 'NGM'] <- "Credible NGM"
  df$credible_status[df$credible == 'Credible' & df$source == 'NaCl'] <- "Credible NaCl"
  df$credible_status <- factor(df$credible_status, levels = c("Not Credible", "Credible NGM", "Credible NaCl"))
  
  # Manhattan plot for combined data
  manhattan_plot <- ggplot(df, aes(x = POS, y = median, color = credible_status, alpha = credible_status)) +
    geom_point(size = 0.5) + 
    scale_color_manual(values = c("Not Credible" = 'grey', "Credible NGM" = 'darkorange', "Credible NaCl" = 'darkgreen')) +
    scale_alpha_manual(values = c("Not Credible" = 0.2, "Credible NGM" = 1, "Credible NaCl" = 1)) +
    theme_minimal() +
    ggtitle(paste("Manhattan Plot Bayésien", trait)) +
    labs(x = "SNPs",
         y = "Effet médian à posteriori") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    facet_wrap(~ CHROM, ncol = 2, scales = "free_x")
  
  ggsave(paste0('Manhattan_plots_', output, '_',pruning,'_',trait,'.pdf'))
  plots_mp[[i]] <- manhattan_plot
  
}

# Function to extract the legend from a ggplot
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Get the legend from the first plot
legend_mp <- get_legend(plots_mp[[1]])

# Remove the legend from all plots
plots_mp <- lapply(plots_mp, function(plot) plot + theme(legend.position = "none"))

# Arrange plots with the legend
pdf(paste0('Manhattan_plots_', output, '_',pruning,'.pdf'), width = 14, height = 10)
grid.arrange(grobs = plots_mp, ncol = 2, top = textGrob("Manhattan Plots", gp = gpar(fontsize = 20, fontface = "bold")))
dev.off()

pdf('Legend_MP.pdf')
grid.arrange(legend_mp)
dev.off()



