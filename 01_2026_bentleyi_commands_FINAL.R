
# ==============================================================================
# Population Genetics Analysis for Corallorhiza bentleyi & C. involuta
# ==============================================================================

# 0. Setup & Package Loading --------------------------------------------------
required_packages <- c(
  "vcfR", "adegenet", "poppr", "hierfstat", "pegas", 
  "LEA", "ggplot2", "ggsci", "viridis", "gridExtra", 
  "scatterpie", "ggrepel", "dplyr", "tidyr", "geosphere", "vegan"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
  library(poppr)
  library(hierfstat)
  library(pegas)
  library(LEA)
  library(ggplot2)
  library(ggsci)
  library(viridis)
  library(gridExtra)
  library(scatterpie)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(geosphere)
  library(vegan)
  library(reshape2)
})


# vcf files
# bentleyi_involuta_filtered_2026_09_14.vcf		2,347 variants		# C. bentleyi, C. involuta, missing data >20% filtered, LD-thinned
# bentleyi_only_tassel_filtered.vcf				52,718 variants		# C. bentleyi only, non-filtered
# filtered_variants.vcf							15,958 variants		# C. bentleyi, C. involuta, missing data >20% filtered	
# coords.cvs														# metadata file: sample	species	pop	latitude	longitude


# Create output directory
if(!dir.exists("figures")) dir.create("figures")

# Load VCF and Metadata
vcf_file  <- "bentleyi_involuta_filtered_2026_09_14.vcf"
meta_file <- "coords.csv"

vcf <- read.vcfR(vcf_file, verbose = FALSE)
coords <- read.csv(meta_file, stringsAsFactors = FALSE)

# Convert to genind and genclone objects
gid <- vcfR2genind(vcf)
pop_info <- coords$pop[match(indNames(gid), coords$sample)]
pop(gid) <- pop_info
strata(gid) <- coords[match(indNames(gid), coords$sample), c("species", "pop")]
gclone <- as.genclone(gid)

# Palette function
igv_pal <- pal_igv()(length(unique(pop_info)))


# ------------------------------------------------------------------------------
# 1. Diversity Metrics (Ho, He, Ho/He, Fis), Stat Tests & plots
# ------------------------------------------------------------------------------

div_summary <- do.call(rbind, lapply(names(pop_list), function(p) {
  sub_gid <- pop_list[[p]]
  stats   <- summary(sub_gid)
  
  # Calculate population means across loci
  mean_Ho  <- mean(stats$Hobs, na.rm = TRUE)
  mean_He  <- mean(stats$Hexp, na.rm = TRUE)
  ho_he    <- mean_Ho / mean_He
  mean_Fis <- 1 - ho_he
  
  data.frame(
    Pop     = p,
    N_Ind   = nInd(sub_gid),
    Mean_Ho = mean_Ho,
    Mean_He = mean_He,
    Ho_He   = ho_he,
    Fis     = mean_Fis
  )
}))

print(div_summary)
write.csv(div_summary, "figures/Table1_Diversity_Metrics.csv", row.names = FALSE)

# 1. Build the per-locus dataset (div_stats) across populations
div_stats <- do.call(rbind, lapply(names(pop_list), function(p) {
  sub_gid <- pop_list[[p]]
  stats   <- summary(sub_gid)
  
  data.frame(
    Locus = names(stats$Hexp),
    Pop   = p,
    Ho    = stats$Hobs,
    He    = stats$Hexp,
    Fis   = 1 - (stats$Hobs / stats$Hexp),
    stringsAsFactors = FALSE
  )
}))

# 2. Perform Kruskal-Wallis Test on Expected Heterozygosity (He)
kw_res <- kruskal.test(He ~ Pop, data = div_stats)
print(kw_res)

# 3. Post-hoc Dunn's Test (requires 'FSA' or 'dunn.test' package)
if (!requireNamespace("FSA", quietly = TRUE)) install.packages("FSA")
library(FSA)

dunn_res <- dunnTest(He ~ Pop, data = div_stats, method = "bh")
print(dunn_res)


# 1. Add the missing Ho_He column to div_stats
div_stats$Ho_He <- div_stats$Ho / div_stats$He

# 2. Build the 4-panel boxplot
p1_1 <- ggplot(div_stats, aes(x = Pop, y = Ho, fill = Pop)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_fill_igv() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + 
  labs(title = "Observed Heterozygosity (Ho)", x = "")

p1_2 <- ggplot(div_stats, aes(x = Pop, y = He, fill = Pop)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_fill_igv() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + 
  labs(title = "Expected Heterozygosity (He)", x = "")

p1_3 <- ggplot(div_stats, aes(x = Pop, y = Ho_He, fill = Pop)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_fill_igv() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + 
  labs(title = "Ho / He Ratio", x = "")

p1_4 <- ggplot(div_stats, aes(x = Pop, y = Fis, fill = Pop)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_fill_igv() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + 
  labs(title = "Inbreeding Coefficient (Fis)", x = "")

# 3. Combine and Save Figure
p1_all <- grid.arrange(p1_1, p1_2, p1_3, p1_4, ncol = 2)
ggsave("figures/Fig1_Diversity_Boxplots.pdf", p1_all, width = 10, height = 8)



# ------------------------------------------------------------------------------
# 2. Isolation by Distance (IBD)
# ------------------------------------------------------------------------------

# Filter to C. bentleyi
bent_samples <- coords$sample[coords$species == "bentleyi"]
gid_bent <- gid[bent_samples, ]
coords_bent <- coords[coords$sample %in% bent_samples, ]

gen_dist <- dist.genpop(genind2genpop(gid_bent, quiet=TRUE))
pop_coords <- coords_bent %>% group_by(pop) %>% summarise(lat = mean(latitude), lon = mean(longitude))
geo_dist <- distm(pop_coords[, c("lon", "lat")], fun = distHaversine) / 1000 # km
geo_dist <- as.dist(geo_dist)

mantel_res <- mantel(gen_dist, geo_dist, permutations = 999)

ibd_df <- data.frame(GenDist = as.vector(gen_dist), GeoDist = as.vector(geo_dist))
p2 <- ggplot(ibd_df, aes(x = GeoDist, y = GenDist)) +
  geom_point(color = "#1F77B4", size = 3) +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  scale_color_igv() +
  theme_bw() +
  labs(title = paste0("Isolation by Distance (Mantel r = ", round(mantel_res$statistic, 3), ", p = ", mantel_res$signif, ")"),
       x = "Geographic Distance (km)", y = "Genetic Distance")
print(p2)
ggsave("figures/Fig2_Isolation_By_Distance.pdf", p2, width = 7, height = 5)




# ------------------------------------------------------------------------------
# 3. Linkage Disequilibrium (IA and rbar-D)
# ------------------------------------------------------------------------------

ld_species <- ia(gclone, sample = 199, quiet = TRUE)

ld_pops <- do.call(rbind, lapply(popNames(gclone), function(p) {
  sub_gc <- popsub(gclone, sublist = p)
  res <- ia(sub_gc, sample = 199, quiet = TRUE)
  data.frame(Pop = p, Ia = res["Ia"], p_Ia = res["p.Ia"], rbarD = res["rbarD"], p_rbarD = res["p.rD"])
}))
ld_table <- rbind(data.frame(Pop = "Overall_Species", Ia = ld_species["Ia"], p_Ia = ld_species["p.Ia"], rbarD = ld_species["rbarD"], p_rbarD = ld_species["p.rD"]), ld_pops)
write.csv(ld_table, "figures/Table2_Linkage_Disequilibrium.csv", row.names = FALSE)




# ------------------------------------------------------------------------------
# 4. Selfing Rate Prediction
# ------------------------------------------------------------------------------

div_stats_ind <- do.call(rbind, lapply(names(pop_list), function(p) {
  sub_gid <- pop_list[[p]]
  stats   <- summary(sub_gid)
  
  # Individual Ho across loci
  ind_Ho  <- rowMeans(tab(sub_gid, NA.method = "zero") == 1, na.rm = TRUE)
  mean_He <- mean(stats$Hexp, na.rm = TRUE)
  ind_Fis <- 1 - (ind_Ho / mean_He)
  
  data.frame(
    Sample = names(ind_Ho),
    Pop    = p,
    Fis    = ind_Fis,
    stringsAsFactors = FALSE
  )
}))

# Calculate s and bound between 0 and 1
div_stats_ind$Selfing_Rate <- (2 * div_stats_ind$Fis) / (1 + div_stats_ind$Fis)
div_stats_ind$Selfing_Rate[div_stats_ind$Selfing_Rate < 0] <- 0
div_stats_ind$Selfing_Rate[div_stats_ind$Selfing_Rate > 1] <- 1

write.csv(div_stats_ind, "figures/Table3_Selfing_Rates_Individual.csv", row.names = FALSE)




# ------------------------------------------------------------------------------
# 5. Private Alleles
# ------------------------------------------------------------------------------


# 1. Private alleles by Population
priv_pop <- private_alleles(gclone)
priv_pop_counts <- rowSums(priv_pop)

# 2. Private alleles by Species
gclone_sp <- gclone
setPop(gclone_sp) <- ~species
priv_sp <- private_alleles(gclone_sp)
priv_sp_counts <- rowSums(priv_sp)

# 3. Combine into Summary Data Frame
priv_df <- data.frame(
  Group = c(names(priv_pop_counts), names(priv_sp_counts)),
  Level = c(rep("Population", length(priv_pop_counts)), rep("Species", length(priv_sp_counts))),
  Private_Alleles = c(priv_pop_counts, priv_sp_counts)
)

print(priv_df)
write.csv(priv_df, "figures/Table4_Private_Alleles.csv", row.names = FALSE)

# 4. Visualization
p5 <- ggplot(priv_df[priv_df$Level == "Population", ], aes(x = Group, y = Private_Alleles, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_igv() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Number of Private Alleles per Population", x = "Population", y = "Private Alleles")

print(p5)
ggsave("figures/Fig5_Private_Alleles.pdf", p5, width = 7, height = 5)




# ------------------------------------------------------------------------------
# 6. AMOVA Analysis
# ------------------------------------------------------------------------------

# --- Task 6: AMOVA Analysis & Significance Testing ---

library(ade4)

# 1. Run AMOVA with within = FALSE to avoid ploidy/dosage warnings
amova_res <- poppr.amova(gclone, ~species/pop, within = FALSE, quiet = TRUE)

# Print Summary Table
print(amova_res)

# Save Components of Variation to CSV
amova_table <- amova_res$componentsofavariation
write.csv(amova_table, "figures/Table5_AMOVA_Results.csv")

# 2. Permutation Significance Test (randtest)
# nrepet = 999 permutations for robust p-values
set.seed(123)
amova_test <- randtest(amova_res, nrepet = 999)

# Display Permutation Test Results
print(amova_test)

# Save Permutation Test Results
randtest_df <- data.frame(
  Test = amova_test$names,
  Observed = amova_test$obs,
  Std_Stat = amova_test$alter,
  P_Value = amova_test$pvalue
)
write.csv(randtest_df, "figures/Table5_AMOVA_Significance.csv", row.names = FALSE)

# 3. Plot Permutation Distribution
pdf("figures/Fig6_AMOVA_Significance_Test.pdf", width = 8, height = 6)
plot(amova_test)
dev.off()




# ------------------------------------------------------------------------------
# 7. Global & Pairwise Fst Matrix Heatmap
# ------------------------------------------------------------------------------


library(reshape2)
hf_data <- genind2hierfstat(gid)
pairwise_fst <- pairwise.WCfst(hf_data)

write.csv(pairwise_fst, "figures/Table6_Pairwise_Fst.csv")

fst_melt <- melt(as.matrix(pairwise_fst))
p7 <- ggplot(fst_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", value)), color = "white", size = 4) +
  scale_fill_viridis_c(option = "magma", na.value = "grey50") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pairwise Fst Matrix", x = "", y = "", fill = "Fst")
print(p7)
ggsave("figures/Fig7_Pairwise_Fst_Heatmap.pdf", p7, width = 8, height = 6)



# ------------------------------------------------------------------------------
# 8. PCA (All samples vs C. bentleyi only)
# ------------------------------------------------------------------------------


# --- 1. Compute Percent Variance ---
# Calculate variance percentage for all samples
var_all <- (pca_all$eig / sum(pca_all$eig)) * 100
# Calculate variance percentage for C. bentleyi samples
var_bent <- (pca_bent$eig / sum(pca_bent$eig)) * 100

# --- 2. Create Named Color Palette ---
# Get all unique populations across both datasets
all_pops <- unique(c(as.character(pop(gid)), as.character(pop(gid_bentleyi))))

# Map each unique population to a fixed color from ggsci::pal_igv()
library(ggsci)
pop_colors <- setNames(pal_igv()(length(all_pops)), all_pops)

# --- 3. Build Plots ---
p8_1 <- ggplot(df_pca_all, aes(x=Axis1, y=Axis2, color=Pop)) +
  geom_point(size=3) +
  scale_color_manual(values = pop_colors) +
  theme_bw() +
  labs(
    title = "All Samples: PC1 vs PC2",
    x = sprintf("Axis 1 (%.2f%%)", var_all[1]),
    y = sprintf("Axis 2 (%.2f%%)", var_all[2])
  )

p8_2 <- ggplot(df_pca_all, aes(x=Axis3, y=Axis4, color=Pop)) +
  geom_point(size=3) +
  scale_color_manual(values = pop_colors) +
  theme_bw() +
  labs(
    title = "All Samples: PC3 vs PC4",
    x = sprintf("Axis 3 (%.2f%%)", var_all[3]),
    y = sprintf("Axis 4 (%.2f%%)", var_all[4])
  )

p8_3 <- ggplot(df_pca_bent, aes(x=Axis1, y=Axis2, color=Pop)) +
  geom_point(size=3) +
  scale_color_manual(values = pop_colors) +
  theme_bw() +
  labs(
    title = "C. bentleyi: PC1 vs PC2",
    x = sprintf("Axis 1 (%.2f%%)", var_bent[1]),
    y = sprintf("Axis 2 (%.2f%%)", var_bent[2])
  )

p8_4 <- ggplot(df_pca_bent, aes(x=Axis3, y=Axis4, color=Pop)) +
  geom_point(size=3) +
  scale_color_manual(values = pop_colors) +
  theme_bw() +
  labs(
    title = "C. bentleyi: PC3 vs PC4",
    x = sprintf("Axis 3 (%.2f%%)", var_bent[3]),
    y = sprintf("Axis 4 (%.2f%%)", var_bent[4])
  )

p8_all <- grid.arrange(p8_1, p8_2, p8_3, p8_4, ncol=2)
p8_all






# ------------------------------------------------------------------------------
# 9. DAPC Analysis
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 9: 4-Panel DAPC Analysis (With/Without Involuta & Prior vs Find.Clusters)
# ------------------------------------------------------------------------------


library(adegenet)
library(ggplot2)
library(ggsci)
library(gridExtra)

# ------------------------------------------------------------------------------
# Helper Function: Run DAPC with optimal PC selection via Cross-Validation
# ------------------------------------------------------------------------------
run_dapc_auto <- function(gen_obj, pop_grp) {
  # Perform cross-validation to select optimal number of PCs
  set.seed(123)
  xval <- xvalDapc(tab(gen_obj, NA.method = "mean"), pop_grp,
                   n.pca.max = min(50, nInd(gen_obj) - 1),
                   training.set = 0.9, result = "groupMean",
                   center = TRUE, scale = FALSE,
                   n.pca = NULL, n.rep = 30, xval.plot = FALSE)
  
  opt_pca <- xval$DAPC$n.pca
  cat("  -> Optimal PCs selected via xval:", opt_pca, "\n")
  
  # Run DAPC with optimal PCs
  dapc_obj <- dapc(gen_obj, pop_grp, n.pca = opt_pca, n.da = 3)
  return(dapc_obj)
}

# ------------------------------------------------------------------------------
# Dataset Setup
# ------------------------------------------------------------------------------
# Full Dataset (All samples)
gid_full <- gid

# Subset Dataset (Exclude C. involuta: Mex-MX and Mor-MX)
# Correct syntax for popsub in poppr:
bentleyi_pops <- setdiff(popNames(gid_full), c("Mex-MX", "Mor-MX"))
gid_bentleyi  <- popsub(gid_full, sublist = bentleyi_pops)

# ------------------------------------------------------------------------------
# Panel 1: Full Dataset + Population Prior
# ------------------------------------------------------------------------------

dapc1 <- run_dapc_auto(gid_full, pop(gid_full))
df1 <- data.frame(dapc1$ind.coord, Group = pop(gid_full))

# Handle 1D vs 2D coordinate plotting
if (!"LD2" %in% colnames(df1)) df1$LD2 <- 0

p1 <- ggplot(df1, aes(x = LD1, y = LD2, color = Group, fill = Group)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", alpha = 0.2, na.rm = TRUE) +
  scale_color_igv() + scale_fill_igv() + theme_bw() +
  labs(title = "A) All Samples (Pop Prior)", x = "LD1", y = "LD2") +
  theme(legend.position = "right")

# ------------------------------------------------------------------------------
# Panel 2: Full Dataset + No Prior (find.clusters)
# ------------------------------------------------------------------------------

set.seed(123)
grp_full <- find.clusters(gid_full, max.n.clust = 10, n.pca = 30, choose.n.clust = FALSE)
cat("  -> Optimal K identified:", max(as.numeric(grp_full$grp)), "\n")

dapc2 <- run_dapc_auto(gid_full, grp_full$grp)
df2 <- data.frame(dapc2$ind.coord, Group = grp_full$grp, Pop = pop(gid_full))
if (!"LD2" %in% colnames(df2)) df2$LD2 <- 0

# Correct function name is scale_color_d3()
p2 <- ggplot(df2, aes(x = LD1, y = LD2, color = Group, fill = Group, shape = Pop)) + 
  geom_point(size = 3) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, na.rm = TRUE) + 
  scale_color_d3() + scale_fill_d3() + theme_bw() + 
  labs(title = paste0("B) All Samples (No Prior, K=", max(as.numeric(grp_full$grp)), ")"), x = "LD1", y = "LD2") + 
  theme(legend.position = "right")


# ------------------------------------------------------------------------------
# Panel 3: C. bentleyi Only + Population Prior
# ------------------------------------------------------------------------------

dapc3 <- run_dapc_auto(gid_bentleyi, pop(gid_bentleyi))
df3 <- data.frame(dapc3$ind.coord, Group = pop(gid_bentleyi))
if (!"LD2" %in% colnames(df3)) df3$LD2 <- 0

p3 <- ggplot(df3, aes(x = LD1, y = LD2, color = Group, fill = Group)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", alpha = 0.2, na.rm = TRUE) +
  scale_color_igv() + scale_fill_igv() + theme_bw() +
  labs(title = "C) C. bentleyi Only (Pop Prior)", x = "LD1", y = "LD2") +
  theme(legend.position = "right")

# ------------------------------------------------------------------------------
# Panel 4: C. bentleyi Only + No Prior (find.clusters)
# ------------------------------------------------------------------------------

set.seed(123)
grp_bentleyi <- find.clusters(gid_bentleyi, max.n.clust = 8, n.pca = 30, choose.n.clust = FALSE)
cat("  -> Optimal K identified:", max(as.numeric(grp_bentleyi$grp)), "\n")

dapc4 <- run_dapc_auto(gid_bentleyi, grp_bentleyi$grp)
df4 <- data.frame(dapc4$ind.coord, Group = grp_bentleyi$grp, Pop = pop(gid_bentleyi))
if (!"LD2" %in% colnames(df4)) df4$LD4 <- 0

# Apply the same fix for Panel 4:
p4 <- ggplot(df4, aes(x = LD1, y = LD2, color = Group, fill = Group, shape = Pop)) + 
  geom_point(size = 3) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, na.rm = TRUE) + 
  scale_color_d3() + scale_fill_d3() + theme_bw() + 
  labs(title = paste0("D) C. bentleyi Only (No Prior, K=", max(as.numeric(grp_bentleyi$grp)), ")"), x = "LD1", y = "LD2") + 
  theme(legend.position = "right")

# ------------------------------------------------------------------------------
# Combine and Save 4-Panel Grid
# ------------------------------------------------------------------------------
p_grid <- grid.arrange(p1, p2, p3, p4, ncol = 2)

ggsave("figures/Fig9_DAPC_4Panel_Comparison.pdf", p_grid, width = 14, height = 11)
ggsave("figures/Fig9_DAPC_4Panel_Comparison.png", p_grid, width = 14, height = 11, dpi = 300)




# ------------------------------------------------------------------------------
# Plot DAPC ancestry components
# ------------------------------------------------------------------------------


library(adegenet)
library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyr)
library(gridExtra)

# ------------------------------------------------------------------------------
# 1. Map Counties from Metadata or Population Names
# ------------------------------------------------------------------------------
# Map population names to counties based on your population string patterns
# (e.g., "PM-Monroe-WV" -> "Monroe", "Mex-MX" -> "Mexico", etc.)
get_county <- function(pop_vec) {
  sapply(as.character(pop_vec), function(p) {
    if (grepl("Monroe", p))  return("Monroe (WV)")
    if (grepl("Giles", p))   return("Giles (VA)")
    if (grepl("Allegany", p)) return("Allegany (VA)")
    if (grepl("Bath", p))     return("Bath (VA)")
    if (grepl("-MX", p))      return("Mexico (Outgroup)")
    return("Other")
  })
}

# Add County metadata into genind objects
gid_full$other$county     <- get_county(pop(gid_full))
gid_bentleyi$other$county <- get_county(pop(gid_bentleyi))

# ------------------------------------------------------------------------------
# 2. Helper Function: Prepare Data & Order Individuals by Ancestry
# ------------------------------------------------------------------------------
prep_ordered_ancestry <- function(dapc_obj, gen_obj, panel_title) {
  post_mat <- as.data.frame(dapc_obj$posterior)
  
  # Ensure standard column naming (Cluster1, Cluster2, ...)
  colnames(post_mat) <- paste0("Cluster_", seq_len(ncol(post_mat)))
  
  df_wide <- post_mat %>%
    mutate(
      Individual = rownames(post_mat),
      Population = as.character(pop(gen_obj)),
      County     = gen_obj$other$county
    )
  
  # Identify primary cluster to order samples within each county
  primary_col <- "Cluster_1"
  
  # Sort individuals within each county by Cluster_1 proportion
  df_ordered <- df_wide %>%
    group_by(County) %>%
    arrange(desc(!!sym(primary_col)), .by_group = TRUE) %>%
    ungroup() %>%
    mutate(Individual = factor(Individual, levels = unique(Individual)))
  
  # Pivot to long format for ggplot
  df_long <- df_ordered %>%
    pivot_longer(
      cols = starts_with("Cluster_"),
      names_to = "Cluster",
      values_to = "Ancestry"
    ) %>%
    mutate(Panel = panel_title)
  
  return(df_long)
}

# ------------------------------------------------------------------------------
# 3. Process Ancestry Data for All 4 DAPC Models
# ------------------------------------------------------------------------------
anc1 <- prep_ordered_ancestry(dapc1, gid_full,     "A) Full Data: Pop Prior")
anc2 <- prep_ordered_ancestry(dapc2, gid_full,     paste0("B) Full Data: No Prior (K=", max(as.numeric(grp_full$grp)), ")"))
anc3 <- prep_ordered_ancestry(dapc3, gid_bentleyi, "C) C. bentleyi Only: Pop Prior")
anc4 <- prep_ordered_ancestry(dapc4, gid_bentleyi, paste0("D) C. bentleyi Only: No Prior (K=", max(as.numeric(grp_bentleyi$grp)), ")"))

# ------------------------------------------------------------------------------
# 4. Plotting Helper Function
# ------------------------------------------------------------------------------
plot_ancestry_panel <- function(df_anc, use_igv = TRUE) {
  p <- ggplot(df_anc, aes(x = Individual, y = Ancestry, fill = Cluster)) +
    geom_col(width = 1, color = "black", linewidth = 0.1) +
    facet_grid(~ County, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(
      axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      axis.ticks.x     = element_blank(),
      panel.spacing    = unit(0.15, "lines"),
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(size = 8, face = "bold"),
      legend.position  = "bottom",
      plot.title       = element_text(size = 11, face = "bold")
    ) +
    labs(
      title = unique(df_anc$Panel),
      x = NULL,
      y = "Posterior Prob.",
      fill = "Cluster"
    )
  
  if (use_igv) {
    p <- p + scale_fill_igv()
  } else {
    p <- p + scale_fill_d3()
  }
  
  return(p)
}

# Generate Individual Bar Plots
bar1 <- plot_ancestry_panel(anc1, use_igv = TRUE)
bar2 <- plot_ancestry_panel(anc2, use_igv = TRUE)
bar3 <- plot_ancestry_panel(anc3, use_igv = TRUE)
bar4 <- plot_ancestry_panel(anc4, use_igv = TRUE)

# ------------------------------------------------------------------------------
# 5. Display 4-Panel Grid
# ------------------------------------------------------------------------------
grid_bar <- grid.arrange(bar1, bar2, bar3, bar4, ncol = 1)

# Save to file
ggsave("figures/DAPC_Ancestry_Barplots_County.pdf", grid_bar, width = 12, height = 14)







# ------------------------------------------------------------------------------
# 10. LEA sNMF Ancestry Barplots & Scatterpie Map
# ------------------------------------------------------------------------------

# Write geno format for LEA
geno_file <- vcf2geno(vcf_file, force = TRUE)

obj_snmf <- snmf(geno_file, K = 2:6, entropy = TRUE, repetitions = 5, project = "new")
best_k <- 3

Q_mat <- Q(obj_snmf, K = best_k)
df_q <- cbind(coords, Q_mat)

# Ancestry Barplot
df_q_long <- pivot_longer(df_q, cols = starts_with("V"), names_to = "Cluster", values_to = "Ancestry")
p10_bar <- ggplot(df_q_long, aes(x = sample, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity") +
  facet_wrap(~pop, scales = "free_x") +
  scale_fill_igv() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = paste0("sNMF Ancestry (K = ", best_k, ")"))
print(p10_bar)
ggsave("figures/Fig10_sNMF_Barplot.pdf", p10_bar, width = 11, height = 6)

# Scatterpie Plot (C. bentleyi mapped geographically)
df_q_bent <- df_q %>% filter(species == "bentleyi") %>%
  group_by(pop) %>%
  summarise(lat = mean(latitude), lon = mean(longitude),
            V1 = mean(V1), V2 = mean(V2), V3 = mean(V3))

p10_map <- ggplot() +
  geom_scatterpie(data = df_q_bent, aes(x = lon, y = lat, r = 0.015),
                  cols = c("V1", "V2", "V3")) +
  scale_fill_igv() +
  coord_fixed() +
  theme_bw() +
  labs(title = "C. bentleyi Ancestry Map", x = "Longitude", y = "Latitude")
print(p10_map)
ggsave("figures/Fig10_sNMF_Spatial_Map.pdf", p10_map, width = 8, height = 6)





##################################
##################################
##################################

### Left off here 8/14/26

# Save entire session
save.image(file = "full_session_workspace.RData")

# Restores all objects back into your environment
load("full_session_workspace.RData")


##################################
##################################
##################################


# ------------------------------------------------------------------------------
# 11. Further investigation of excess Ho and negative Fis 
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Test 2A: Fis vs. Minor Allele Frequency (MAF)
# ------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)

# Extract MAF and Fis per locus
maf_vec <- minorAllele(gid_bentleyi)
fis_vec <- rowMeans(pop_stats$Fis, na.rm = TRUE)

df_maf_fis <- data.frame(MAF = maf_vec, Fis = fis_vec) %>% filter(!is.na(Fis))

# Fit a linear or non-linear model
lm_fit <- lm(Fis ~ MAF, data = df_maf_fis)
summary(lm_fit)

# Plot Fis vs MAF
fisplot <- ggplot(df_maf_fis, aes(x = MAF, y = Fis)) +
  geom_point(alpha = 0.3, color = "darkblue") +
  geom_smooth(method = "lm", color = "red") +
  theme_bw() +
  labs(title = "Locus-Level Fis vs. Minor Allele Frequency",
       subtitle = "Negative correlation indicates overdominance/purging",
       x = "Minor Allele Frequency (MAF)", y = "Fis")
fisplot	   


# Call:
# lm(formula = Fis ~ MAF, data = df_maf_fis)
# 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -1.01535 -0.07756 -0.04382 -0.00181  1.25519 
# 
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.106136   0.005865   18.10   <2e-16 ***
# MAF         -1.343663   0.039087  -34.38   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1672 on 2247 degrees of freedom
# Multiple R-squared:  0.3447,	Adjusted R-squared:  0.3444 
# F-statistic:  1182 on 1 and 2247 DF,  p-value: < 2.2e-16

# ------------------------------------------------------------------------------
# Test 2B: Fis as a function of Minor Allele Frequency (MAF), across all populations
# ------------------------------------------------------------------------------

library(hierfstat)
library(dplyr)

# Extract locus-level basic stats
pop_stats <- basic.stats(gid_bentleyi)

# Get sample sizes per population
pop_sizes <- table(pop(gid_bentleyi))

# Iterate through populations and calculate expected vs observed Fis
results_list <- list()

for (p in names(pop_sizes)) {
  n_ind <- pop_sizes[p]
  
  # Skip populations with N < 3 as variance is too high/undefined
  if (n_ind < 3) next 
  
  exp_fis <- -1 / (2 * n_ind - 1)
  
  # Extract locus Fis vector
  fis_vec <- pop_stats$Fis[, p]
  fis_vec <- fis_vec[!is.nan(fis_vec) & !is.na(fis_vec)]
  
  # Perform t-test comparing observed mean against small N expectation
  ttest <- t.test(fis_vec, mu = exp_fis, alternative = "less")
  
  results_list[[p]] <- data.frame(
    Population = p,
    N_Ind = n_ind,
    Observed_Fis = mean(fis_vec),
    Expected_Fis_SmallN = exp_fis,
    t_stat = ttest$statistic,
    p_value = ttest$p.value
  )
}

fis_summary_table <- bind_rows(results_list)
print(fis_summary_table)


#                    Population N_Ind Observed_Fis Expected_Fis_SmallN     t_stat      p_value
# PM-Monroe-WV     PM-Monroe-WV     5  -0.14032680         -0.11111111 -3.9308784 4.425192e-05
# OR-Giles-VA       OR-Giles-VA     8  -0.11494206         -0.06666667 -7.6495779 1.623988e-14
# WR-Giles-VA       WR-Giles-VA     5  -0.02429332         -0.11111111 10.1163131 1.000000e+00
# RG-Allegany-VA RG-Allegany-VA     5  -0.09717424         -0.11111111  1.7765558 9.620637e-01
# CG-Bath-VA         CG-Bath-VA    14  -0.04184985         -0.03703704 -0.7964541 2.129325e-01
# 
# 


# ------------------------------------------------------------------------------
# Test 2C: Two-tailed test for populations that are either more negative or less negative than the small-N expectation
# ------------------------------------------------------------------------------


library(hierfstat)
library(dplyr)

pop_stats <- basic.stats(gid_bentleyi)
pop_sizes <- table(pop(gid_bentleyi))

results_list <- list()

for (p in names(pop_sizes)) {
  n_ind <- pop_sizes[p]
  if (n_ind < 3) next 
  
  exp_fis <- -1 / (2 * n_ind - 1)  # The Robertson (1965) expectation under small-N sampling
  fis_vec <- pop_stats$Fis[, p]
  fis_vec <- fis_vec[!is.nan(fis_vec) & !is.na(fis_vec)]
  
  # Two-tailed test
  ttest <- t.test(fis_vec, mu = exp_fis, alternative = "two.sided")
  
  results_list[[p]] <- data.frame(
    Population = p,
    N = n_ind,
    Obs_Fis = round(mean(fis_vec), 4),
    Exp_Fis_SmallN = round(exp_fis, 4),
    t_stat = round(ttest$statistic, 3),
    p_value = format.pval(ttest$p.value, digits = 3),
    Inference = case_when(
      ttest$p.value > 0.05 ~ "Small-N Sampling Bias",
      mean(fis_vec) < exp_fis ~ "Overdominance / Purging",
      mean(fis_vec) > exp_fis ~ "Inbreeding / Selfing Shift"
    )
  )
}

final_table <- bind_rows(results_list)
print(final_table)



# ------------------------------------------------------------------------------
# Test 2C: Locus-Level Heterozygosity Excess (Ho-He)
# ------------------------------------------------------------------------------


# Load Unpruned VCF and Create `gid_bent_unpruned`

library(vcfR)
library(adegenet)
library(poppr)
library(hierfstat)
library(ggplot2)
library(dplyr)

vcf_unpruned <- read.vcfR("bentleyi_only_tassel_filtered.vcf", verbose = FALSE)
gid_unpruned <- vcfR2genind(vcf_unpruned)

# Assign populations to samples using metadata
pop(gid_unpruned) <- coords$pop[match(indNames(gid_unpruned), coords$sample)]

# Subset to C. bentleyi only (excluding C. involuta outgroups if present, just as a check)
bentleyi_pops <- setdiff(popNames(gid_unpruned), c("Mex-MX", "Mor-MX"))
gid_bent_unpruned <- popsub(gid_unpruned, sublist = bentleyi_pops)

# get basic statts
basic_unpruned <- basic.stats(gid_bent_unpruned)

# Delta H = Ho-Hs
ho_vec <- rowMeans(basic_unpruned$Ho, na.rm = TRUE)
he_vec <- rowMeans(basic_unpruned$Hs, na.rm = TRUE)
delta_h <- ho_vec - he_vec

library(ggplot2)

#  Plot
p_dh <- ggplot(df_dh, aes(x = Delta_H)) +
  geom_histogram(bins = 50, fill = "midnightblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  theme_bw() +
  labs(
    title = "Locus-Level Heterozygosity Excess (Ho - He) in C. bentleyi",
    subtitle = paste0("Mean Delta H = +", round(mean(df_dh$Delta_H), 4), 
                      " (t = ", round(t_test_dh$statistic, 2), ", p < 2.2e-16)"),
    x = "Observed - Expected Heterozygosity (Ho - He)",
    y = "Locus Count"
  )

print(p_dh)


# ------------------------------------------------------------------------------
# 3. Test 2D: Individual Multilocus Heterozygosity (MLH)
# ------------------------------------------------------------------------------

library(inbreedR)

# Convert genind object to 0/1/2 matrix
gen_mat_unpruned <- tab(gid_bent_unpruned, NA.method = "mean")
mlh_inds <- MLH(gen_mat_unpruned)

df_mlh <- data.frame(
  Sample = names(mlh_inds),
  MLH = mlh_inds,
  Population = pop(gid_bent_unpruned)
)

print(head(df_mlh))

# Plot MLH across populations
p_mlh <- ggplot(df_mlh, aes(x = Population, y = MLH, fill = Population)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2.5, color = "black") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "Individual Multilocus Heterozygosity (MLH) Across Populations",
    x = "Population",
    y = "Proportion of Heterozygous Loci (MLH)"
  )



# 	One Sample t-test
# 
# data:  df_dh$Delta_H
# t = 98.778, df = 52717, p-value < 2.2e-16
# alternative hypothesis: true mean is greater than 0
# 95 percent confidence interval:
#  0.05893114        Inf
# sample estimates:
#  mean of x 
# 0.05992909 
# 
# 
# 
# print(p_mlh)


# df_mlh


# ------------------------------------------------------------------------------
# 3. Test 2E:  Individual Heterozygosity Concordance
# ------------------------------------------------------------------------------


library(vcfR)

# Extract GT matrix on unpruned data
vcf_unpruned <- read.vcfR("bentleyi_only_tassel_filtered.vcf", verbose = FALSE)
gt_mat <- extract.gt(vcf_unpruned, element = "GT")

# Calculate proportion of samples that are heterozygous per locus
is_het <- gt_mat == "0/1" | gt_mat == "1/0" | gt_mat == "0|1" | gt_mat == "1|0"
het_rate_per_locus <- rowMeans(is_het, na.rm = TRUE)

# Count how many loci are heterozygous in >80% or >90% of ALL individuals
paralog_candidates <- sum(het_rate_per_locus >= 0.85, na.rm = TRUE)
total_loci <- length(het_rate_per_locus)

cat(sprintf("Loci Heterozygous in >= 85%% of Individuals: %d / %d (%.2f%%)\n",
            paralog_candidates, total_loci, (paralog_candidates / total_loci) * 100))
			
			
### Loci Heterozygous in >= 85% of Individuals: 2942 / 52718 (5.58%)

# Filter Out paralog candidates

# Identify true single-copy loci (het rate < 0.85)
keep_loci <- names(het_rate_per_locus[het_rate_per_locus < 0.85 & !is.na(het_rate_per_locus)])

# Subset genind object
gid_bent_cleaned <- gid_bent_unpruned[, loc = keep_loci]

cat(sprintf("Retained %d single-copy loci after removing paralog candidates.\n", nLoc(gid_bent_cleaned)))

# Recalculate Delta H on filtered dataset
basic_clean <- basic.stats(gid_bent_cleaned)
ho_clean    <- rowMeans(basic_clean$Ho, na.rm = TRUE)
he_clean    <- rowMeans(basic_clean$Hs, na.rm = TRUE)
delta_h_clean <- ho_clean - he_clean

t_clean <- t.test(delta_h_clean, mu = 0, alternative = "greater")
cat(sprintf("Filtered Mean Delta H = +%.4f (t = %.2f, p = %.4e)\n", 
            mean(delta_h_clean, na.rm = TRUE), t_clean$statistic, t_clean$p.value))

# Filtered Mean Delta H = +0.0382 (t = 77.95, p = 0.0000e+00)



