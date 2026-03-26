# ============================================================================
# Microbiota Composition Analysis: AD vs Controls across timepoints
# Shannon Diversity, PCoA, and PERMANOVA Analysis
# ============================================================================

# 1) Clear environment
# ============================================================================
rm(list = ls())
gc()

# Load required libraries
library(vegan)        
library(ggplot2)      
library(dplyr)      
library(tidyr)       
library(ape)       
library(gridExtra)   
library(grid)         
library(readxl)       
library(writexl)      
library(RColorBrewer) 

# ============================================================================
# 2) Import data - CSV file
# ============================================================================
data <- read.csv("Normalized_abundance_Species.csv")

# ============================================================================
# 3) Data preparation and define timepoints
# ============================================================================
timepoints_of_interest <- c("PGV2", "PGV5", "W3", "M3", "M6")

data_filtered <- data %>%
  filter(Timepoints %in% timepoints_of_interest) %>%
  filter(!is.na(Eczema_M18_2Groups))

data_filtered$AD_status <- factor(data_filtered$Eczema_M18_2Groups,
                                  levels = c(0, 1),
                                  labels = c("Controls", "AD"))

data_filtered$Timepoints <- factor(data_filtered$Timepoints,
                                   levels = timepoints_of_interest)

metadata_cols <- c("Subject_TP", "Timepoints", "Eczema_M18_2Groups", "AD_status")
species_cols <- setdiff(names(data_filtered), metadata_cols)

species_abundance <- data_filtered[, species_cols]
species_abundance[is.na(species_abundance)] <- 0  

# ============================================================================
# 4) Shannon Diversity Analysis
# ============================================================================
data_filtered$Shannon_Diversity <- diversity(species_for_analysis, index = "shannon")

shannon_results_list <- list()

for (tp in timepoints_of_interest) {
  # Subset data for current timepoint
  tp_data <- data_filtered %>% filter(Timepoints == tp)
  
  if (nrow(tp_data) > 0 && length(unique(tp_data$AD_status)) == 2) {
    # Wilcoxon test (Mann-Whitney U test)
    test_result <- wilcox.test(Shannon_Diversity ~ AD_status, data = tp_data)
    
    # Calculate summary statistics
    summary_stats <- tp_data %>%
      group_by(AD_status) %>%
      summarise(
        n = n(),
        mean = mean(Shannon_Diversity, na.rm = TRUE),
        sd = sd(Shannon_Diversity, na.rm = TRUE),
        median = median(Shannon_Diversity, na.rm = TRUE),
        IQR = IQR(Shannon_Diversity, na.rm = TRUE)
      )
    
    # Store results
    shannon_results_list[[tp]] <- data.frame(
      Timepoint = tp,
      summary_stats,
      p_value = test_result$p.value,
      test_method = "Wilcoxon rank-sum test"
    )
  }
}

# Combine results
shannon_results <- bind_rows(shannon_results_list)

# ============================================================================
# 5) Beta Diversity Analysis - Bray-Curtis Dissimilarity and PCoA
# ============================================================================
# Calculate Bray-Curtis dissimilarity matrix
bray_dist <- vegdist(species_transformed, method = "bray")

# Perform PCoA (Principal Coordinates Analysis)
pcoa_result <- pcoa(bray_dist)

# Extract coordinates for plotting (first 2 axes)
pcoa_coords <- as.data.frame(pcoa_result$vectors[, 1:2])
colnames(pcoa_coords) <- c("PCoA1", "PCoA2")

# Add metadata to PCoA coordinates
pcoa_data <- cbind(pcoa_coords, data_filtered[, c("Timepoints", "AD_status", "Subject_TP")])

# Calculate variance explained by each axis
variance_explained <- pcoa_result$values$Relative_eig[1:2] * 100

# ============================================================================
# 6) PERMANOVA Analysis (Beta Diversity Testing)
# ============================================================================
permanova_results_list <- list()

for (tp in timepoints_of_interest) {
  # Get indices for current timepoint
  tp_indices <- which(data_filtered$Timepoints == tp)
  
  if (length(tp_indices) > 0 && length(unique(data_filtered$AD_status[tp_indices])) == 2) {
    # Subset distance matrix and metadata
    tp_dist <- as.dist(as.matrix(bray_dist)[tp_indices, tp_indices])
    tp_metadata <- data_filtered[tp_indices, ]
    
    # Perform PERMANOVA
    set.seed(123) 
    permanova_result <- adonis2(tp_dist ~ AD_status, 
                                data = tp_metadata, 
                                permutations = 999)
    
    # Test for homogeneity of dispersions (betadisper)
    set.seed(123)
    disp <- betadisper(tp_dist, tp_metadata$AD_status)
    disp_test <- permutest(disp, permutations = 999)
    
    # Extract results
    permanova_results_list[[tp]] <- data.frame(
      Timepoint = tp,
      n_samples = length(tp_indices),
      n_AD = sum(tp_metadata$AD_status == "AD"),
      n_Controls = sum(tp_metadata$AD_status == "Controls"),
      R2 = permanova_result$R2[1],
      F_statistic = permanova_result$F[1],
      PERMANOVA_p_value = permanova_result$`Pr(>F)`[1],
      Dispersion_F = disp_test$tab$F[1],
      Dispersion_p_value = disp_test$tab$`Pr(>F)`[1],
      Interpretation = ifelse(permanova_result$`Pr(>F)`[1] < 0.05 & 
                                disp_test$tab$`Pr(>F)`[1] >= 0.05,
                              "Significant difference in composition",
                              ifelse(permanova_result$`Pr(>F)`[1] < 0.05 & 
                                       disp_test$tab$`Pr(>F)`[1] < 0.05,
                                     "Significant - but check dispersion",
                                     "No significant difference"))
    )
  }
}

# Combine PERMANOVA results
permanova_results <- bind_rows(permanova_results_list)

# ============================================================================
# 7) PLOTTING
# ============================================================================

# Plot 1: Shannon Diversity Boxplots by Timepoint with Significance
p1 <- ggplot(data_filtered, aes(x = Timepoints, y = Shannon_Diversity, fill = AD_status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2, outlier.alpha = 0.6) +
  scale_fill_manual(values = c("Controls" = "#4DAF4A", "AD" = "#E41A1C")) +
  labs(title = "Shannon Diversity: AD vs Controls across Timepoints",
       x = "Timepoints",
       y = "Shannon Diversity Index",
       fill = "Group") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold")
  )

if (nrow(shannon_results) > 0) {
  y_max <- max(data_filtered$Shannon_Diversity, na.rm = TRUE)
  y_range <- diff(range(data_filtered$Shannon_Diversity, na.rm = TRUE))
  y_position <- y_max + 0.1 * y_range
  
  sig_labels <- shannon_results %>%
    select(Timepoint, p_value) %>%
    distinct(Timepoint, .keep_all = TRUE) %>%
    mutate(
      sig_label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      display_label = ifelse(sig_label == "ns", "", sig_label)
    )
  
  for (i in 1:nrow(sig_labels)) {
    if (sig_labels$display_label[i] != "") {
      p1 <- p1 + 
        annotate("text", 
                 x = which(timepoints_of_interest == sig_labels$Timepoint[i]),
                 y = y_position,
                 label = sig_labels$display_label[i],
                 size = 6,
                 fontface = "bold",
                 color = "black")
    }
  }
  
  p1 <- p1 +
    annotate("text",
             x = match(sig_labels$Timepoint, timepoints_of_interest),
             y = y_position - 0.05 * y_range,
             label = sprintf("p=%.3f", sig_labels$p_value),
             size = 3,
             color = "gray30",
             fontface = "italic")
}

ggsave("Shannon_Diversity_Boxplot.png", plot = p1, width = 10, height = 6, dpi = 300)
print(p1)

# ============================================================================
# Individual PCoA Plots 
# ============================================================================

pcoa_plots_enhanced <- list()

for (tp in timepoints_of_interest) {
  # Get indices for current timepoint
  tp_indices <- which(data_filtered$Timepoints == tp)
  
  if (length(tp_indices) > 0 && length(unique(data_filtered$AD_status[tp_indices])) == 2) {
    
    # Subset distance matrix for this timepoint
    tp_dist <- as.dist(as.matrix(bray_dist)[tp_indices, tp_indices])
    
    # Recalculate PCoA for this timepoint subset
    tp_pcoa_result <- pcoa(tp_dist)
    
    # Extract coordinates
    tp_pcoa_coords <- as.data.frame(tp_pcoa_result$vectors[, 1:2])
    colnames(tp_pcoa_coords) <- c("PCoA1", "PCoA2")
    
    # Calculate variance 
    tp_variance_explained <- tp_pcoa_result$values$Relative_eig[1:2] * 100
    
    # Add metadata
    tp_pcoa_data <- cbind(tp_pcoa_coords, data_filtered[tp_indices, c("Timepoints", "AD_status", "Subject_TP")])
    
    # Get PERMANOVA results for this timepoint
    perm_result <- permanova_results %>% filter(Timepoint == tp)
    
    if (nrow(perm_result) > 0) {
      p_val <- perm_result$PERMANOVA_p_value
      r2_val <- perm_result$R2
      f_val <- perm_result$F_statistic
      disp_p <- perm_result$Dispersion_p_value
      n_ad <- sum(tp_pcoa_data$AD_status == "AD")
      n_controls <- sum(tp_pcoa_data$AD_status == "Controls")
      
      # Determine significance
      if (p_val < 0.001) {
        sig_label <- "p < 0.001***"
        sig_status <- "Significant"
      } else if (p_val < 0.01) {
        sig_label <- sprintf("p = %.3f**", p_val)
        sig_status <- "Significant"
      } else if (p_val < 0.05) {
        sig_label <- sprintf("p = %.3f*", p_val)
        sig_status <- "Significant"
      } else {
        sig_label <- sprintf("p = %.3f (ns)", p_val)
        sig_status <- "Not Significant"
      }
      
      # Create custom legend labels
      legend_labels <- c(
        sprintf("Controls (n=%d)", n_controls),
        sprintf("AD (n=%d)", n_ad)
      )
      
      p <- ggplot(tp_pcoa_data, aes(x = PCoA1, y = PCoA2, color = AD_status, fill = AD_status)) +
        geom_point(size = 4, alpha = 0.8, shape = 16) +
        stat_ellipse(aes(group = AD_status), type = "t", level = 0.95, 
                     linetype = 2, linewidth = 1, alpha = 0.2) +
        scale_color_manual(
          name = NULL,
          values = c("Controls" = "#4DAF4A", "AD" = "#E41A1C"),
          labels = legend_labels
        ) +
        scale_fill_manual(
          name = NULL,
          values = c("Controls" = "#4DAF4A", "AD" = "#E41A1C"),
          labels = legend_labels,
          guide = "none"
        ) +
        labs(
          title = sprintf("Timepoint: %s", tp),
          subtitle = sprintf("PERMANOVA: %s (R² = %.3f)", sig_status, r2_val),
          x = paste0("PCoA1 (", round(tp_variance_explained[1], 2), "%)"),
          y = paste0("PCoA2 (", round(tp_variance_explained[2], 2), "%)")
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11, 
                                       color = ifelse(p_val < 0.05, "#E41A1C", "#666666"),
                                       face = ifelse(p_val < 0.05, "bold", "plain")),
          legend.position = "bottom",
          legend.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
          legend.box.background = element_rect(fill = "grey95", color = "black", linewidth = 1),
          legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
          legend.text = element_text(size = 11),
          legend.key.size = unit(1.2, "lines"),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10)
        ) +
        guides(color = guide_legend(
          title = sprintf("PERMANOVA: %s\nF = %.2f, %s\nBetadisper p = %.3f",
                          sig_status, f_val, sig_label, disp_p),
          title.theme = element_text(size = 10, face = "bold.italic", lineheight = 1.3),
          override.aes = list(size = 5, alpha = 1)
        ))
      
      # Store and save plot
      pcoa_plots_enhanced[[tp]] <- p
      filename <- sprintf("PCoA_%s_Enhanced.png", tp)
      ggsave(filename, plot = p, width = 9, height = 8, dpi = 300)
      
      print(p)

    }
  }
}

# Create combined grid plot
combined_plot <- grid.arrange(
  grobs = pcoa_plots_enhanced,
  ncol = 3,
  top = grid::textGrob("PCoA Analysis by Timepoint: AD vs Controls\n(Bray-Curtis Dissimilarity)", 
                       gp = grid::gpar(fontface = "bold", fontsize = 16))
)

ggsave("PCoA_All_Timepoints_Individual.png", 
       plot = combined_plot, 
       width = 18, height = 12, dpi = 300)

