# ------------------------------------------------------------------------------
# 11. Further investigation of excess Ho and negative Fis
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Test 2A: Fis as a function of Minor Allele Frequency (MAF)
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
# Test 2B: Fis as a function of Minor Allele Frequency (MAF)
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

# ------------------------------------------------------------------------------
# Test 2B: Two-tailed test for populations that are either more negative or less negative than the small-$N$ expectation
# ------------------------------------------------------------------------------


library(hierfstat)
library(dplyr)

pop_stats <- basic.stats(gid_bentleyi)
pop_sizes <- table(pop(gid_bentleyi))

results_list <- list()

for (p in names(pop_sizes)) {
  n_ind <- pop_sizes[p]
  if (n_ind < 3) next 
  
  exp_fis <- -1 / (2 * n_ind - 1)
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
# 2. Test 2C: Locus-Level Heterozygosity Excess (Ho - He)
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

# Assign populations to samples using your coords data
pop(gid_unpruned) <- coords$pop[match(indNames(gid_unpruned), coords$sample)]

# Subset to C. bentleyi only (excluding C. involuta outgroups if present)
bentleyi_pops <- setdiff(popNames(gid_unpruned), c("Mex-MX", "Mor-MX"))
gid_bent_unpruned <- popsub(gid_unpruned, sublist = bentleyi_pops)

basic_unpruned <- basic.stats(gid_bent_unpruned)

# Delta H = Ho - Hs (Hs is gene diversity / expected heterozygosity He)
ho_vec <- rowMeans(basic_unpruned$Ho, na.rm = TRUE)
he_vec <- rowMeans(basic_unpruned$Hs, na.rm = TRUE)
delta_h <- ho_vec - he_vec

library(ggplot2)

# Corrected ggplot script using standard color names and 'linewidth'
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


> df_mlh
#        Sample       MLH     Population
# 2250     2250 0.5295994   PM-Monroe-WV
# 257A     257A 0.4767843    OR-Giles-VA
# 257D     257D 0.3814328    OR-Giles-VA
# 257G     257G 0.5188867    OR-Giles-VA
# 257K     257K 0.4872757    OR-Giles-VA
# 257L     257L 0.4073619    OR-Giles-VA
# 423a     423a 0.3063187    WR-Giles-VA
# 423c     423c 0.4330313    WR-Giles-VA
# 424a     424a 0.5117286   PM-Monroe-WV
# 424b     424b 0.4496961   PM-Monroe-WV
# 424c     424c 0.4703158   PM-Monroe-WV
# 424d     424d 0.4992015   PM-Monroe-WV
# 426a     426a 0.3707822 RG-Allegany-VA
# 426c     426c 0.4911656 RG-Allegany-VA
# 426d     426d 0.3481599 RG-Allegany-VA
# 64a       64a 0.4838121    OR-Giles-VA
# 64b       64b 0.5173461    OR-Giles-VA
# 65a       65a 0.5079523    OR-Giles-VA
# Alleg1 Alleg1 0.5002194 RG-Allegany-VA
# Alleg2 Alleg2 0.5255550 RG-Allegany-VA
# B10       B10 0.4296720     CG-Bath-VA
# B11       B11 0.4504450     CG-Bath-VA
# B12       B12 0.4501444     CG-Bath-VA
# B14       B14 0.4383589     CG-Bath-VA
# B15       B15 0.5074889     CG-Bath-VA
# B16       B16 0.4432704     CG-Bath-VA
# B18       B18 0.4570557     CG-Bath-VA
# B2         B2 0.4777550     CG-Bath-VA
# B20       B20 0.4795792     CG-Bath-VA
# B21       B21 0.3723921     CG-Bath-VA
# B3         B3 0.5191054     CG-Bath-VA
# B5         B5 0.4774355     CG-Bath-VA
# B8         B8 0.5542773     CG-Bath-VA
# B9         B9 0.4694312     CG-Bath-VA
# WR2       WR2 0.4354556    WR-Giles-VA
# WR3       WR3 0.4544836    WR-Giles-VA
# WR4       WR4 0.5365425    WR-Giles-VA



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
			
			
# Loci Heterozygous in >= 85% of Individuals: 2942 / 52718 (5.58%)


# ------------------------------------------------------------------------------
# 3. Test 2F:  Filter out paralog candidates
# ------------------------------------------------------------------------------

# 1. Identify true single-copy loci (het rate < 0.85)
keep_loci <- names(het_rate_per_locus[het_rate_per_locus < 0.85 & !is.na(het_rate_per_locus)])

# 2. Subset genind object
gid_bent_cleaned <- gid_bent_unpruned[, loc = keep_loci]

# 3. Recalculate Delta H on filtered dataset
basic_clean <- basic.stats(gid_bent_cleaned)
ho_clean    <- rowMeans(basic_clean$Ho, na.rm = TRUE)
he_clean    <- rowMeans(basic_clean$Hs, na.rm = TRUE)
delta_h_clean <- ho_clean - he_clean

t_clean <- t.test(delta_h_clean, mu = 0, alternative = "greater")

# Filtered Mean Delta H = +0.0382 (t = 77.95, p = 0.0000e+00)



######################



######################



########################



library(tidyverse)
library(ggsci)
library(cowplot)

out_dir <- "bentleyi_het_tests"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Color Palette
all_pops <- c("CG-Bath-VA", "OR-Giles-VA", "PM-Monroe-WV", "RG-Allegany-VA", "WR-Giles-VA")
pop_colors <- setNames(pal_igv()(length(all_pops)), all_pops)

# ------------------------------------------------------------------------------
# Panel A: Locus Ho Distribution (Rule out Clonality/Paralogy Spike)
# ------------------------------------------------------------------------------
df_ho <- data.frame(Ho = rowMeans(basic_unpruned$Ho, na.rm = TRUE)) %>% filter(!is.na(Ho))

panel_a <- ggplot(df_ho, aes(x = Ho)) +
  geom_histogram(binwidth = 0.04, fill = "#000080", color = "white", boundary = 0) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(-0.02, 1.03)) +
  theme_bw(base_size = 11) +
  labs(
    title = "A. Locus-Level Ho Distribution",
    subtitle = "Rules out clonality/paralogs (no spike at 1.0)",
    x = "Observed Heterozygosity (Ho)",
    y = "Locus Count"
  ) +
  theme(plot.title = element_text(face = "bold"))

# ------------------------------------------------------------------------------
# Panel B: Delta H (Fixed Subtitle)
# ------------------------------------------------------------------------------
panel_b <- ggplot(df_dh, aes(x = Delta_H)) +
  geom_histogram(bins = 45, fill = "#1E3A8A", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
  theme_bw(base_size = 11) +
  labs(
    title = expression(paste("B. Heterozygosity Excess (", H[o] - H[e], ")")),
    subtitle = expression(paste("Significant shift (", Delta, "H = +0.038, ", p < 10^-15, ")")),
    x = expression(paste(H[o] - H[e])),
    y = "Locus Count"
  ) +
  theme(plot.title = element_text(face = "bold"))
  
# ------------------------------------------------------------------------------
# Panel C: Fis vs. MAF (Purging Mechanism)
# ------------------------------------------------------------------------------
# Make sure hexbin is installed: install.packages("hexbin")

panel_c <- ggplot(df_maf_fis, aes(x = MAF, y = Fis)) +
    # 2D Hexagonal Binning reduces >50,000 points down to a tiny grid
    stat_bin_hex(bins = 60, color = NA) +
    scale_fill_viridis_c(
        option = "plasma", 
        trans = "log10", 
        name = "Locus\nCount"
    ) +
    geom_smooth(method = "lm", color = "black", linewidth = 2) +
    theme_bw(base_size = 11) +
    labs(
        title = expression(paste("C. ", F[IS], " vs. Minor Allele Frequency")),
        subtitle = expression(paste("Purging signal (", beta, " = -1.34, ", R^2, " = 0.34)")),
        x = "Minor Allele Frequency (MAF)",
        y = expression(F[IS])
    ) +
    theme(
        plot.title = element_text(face = "bold"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)
    )

# ------------------------------------------------------------------------------
# Panel D: Expanded Axis Range for Highly Negative Fis (PM-Monroe Fix)
# ------------------------------------------------------------------------------
plot_data_d <- df_small_n %>%
  mutate(
    Pop_Label = sprintf("%s (N=%d)", Population, N),
    Sig_Label = case_when(
      Inference == "Overdominance / Purging" ~ "***",
      Inference == "Inbreeding / Selfing Shift" ~ "***",
      TRUE ~ "ns"
    )
  )

panel_d <- ggplot(plot_data_d, aes(x = reorder(Pop_Label, N))) +
  # Segment connecting expected to observed
  geom_segment(
    aes(xend = Pop_Label, y = Exp_Fis_SmallN, yend = Obs_Fis), 
    color = "gray50", linewidth = 0.7, linetype = "dashed"
  ) +
  # White diamond for Expected Small-N Fis
  geom_point(
    aes(y = Exp_Fis_SmallN), 
    shape = 23, size = 3.5, fill = "white", color = "black", stroke = 1
  ) +
  # Colored point for Observed Fis
  geom_point(
    aes(y = Obs_Fis, color = Population), 
    size = 4.5
  ) +
  # Significance text label positioned cleanly to the left of the lowest value
  geom_text(
    aes(y = pmin(Obs_Fis, Exp_Fis_SmallN) - 0.012, label = Sig_Label), 
    color = "darkred", fontface = "bold", size = 4
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
  scale_color_manual(values = pop_colors, guide = "none") +
  # Expanded lower limit to -0.35 handles extreme PM-Monroe negative values easily
  coord_flip(ylim = c(-0.35, 0.02), clip = "off") +
  scale_y_continuous(breaks = seq(-0.35, 0.0, 0.05)) +
  labs(
    title = expression(paste("D. Observed vs. Expected Small-N ", F[IS])),
    subtitle = expression(paste(diamond, " = Small-N Model (", -1/(2*N - 1), ")")),
    x = "Population",
    y = expression(F[IS])
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 10)
  )
  

# ------------------------------------------------------------------------------
# Panel E: Individual MLH Across Populations
# ------------------------------------------------------------------------------
panel_e <- ggplot(df_mlh, aes(x = Population, y = MLH, fill = Population)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2, shape = 21, color = "black") +
  scale_fill_manual(values = pop_colors, guide = "none") +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1)
  ) +
  labs(
    title = "E. Individual Multilocus Heterozygosity",
    subtitle = "High baseline (0.45-0.50) across all populations",
    x = "Population",
    y = "Individual MLH"
  )

# ------------------------------------------------------------------------------
# Combine Panels into Master Composite Figure
# ------------------------------------------------------------------------------
top_row <- plot_grid(panel_a, panel_b, panel_c, ncol = 3)
bottom_row <- plot_grid(panel_d, panel_e, ncol = 2, rel_widths = c(1.1, 1))

master_figure <- plot_grid(top_row, bottom_row, nrow = 2, rel_heights = c(1, 1.1))

print(master_figure)

# Save high-resolution PDF for publication
ggsave(
  filename = file.path(out_dir, "Figure_Master_Heterozygosity_Synthesis.pdf"),
  plot = master_figure,
  width = 13,
  height = 8.5,
  useDingbats = FALSE
)




# ==============================================================================
# Dual-Dataset Bottleneck & Overdominance Pipeline
# Uses gid_full for Bottlenecks & gid_bent_unpruned for Overdominance
# ==============================================================================

library(adegenet)
library(pegas)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# 1. PRE-PROCESSING
# ------------------------------------------------------------------------------

gid_full_split <- seppop(gid_full)                # LD-thinned for bottlenecks
gid_unpruned_split <- seppop(gid_bent_unpruned)   # Unpruned for overdominance

# ------------------------------------------------------------------------------
# 2. BOTTLENECK TESTS (Run on LD-Thinned Dataset: gid_full)
# ------------------------------------------------------------------------------

run_bottleneck_analysis <- function(gid_list) {
  results <- lapply(names(gid_list), function(sp_name) {
    sub_gid <- gid_list[[sp_name]]
    
    pop_summary <- summary(sub_gid)
    n_ind <- nInd(sub_gid)
    n_loci <- nLoc(sub_gid)
    
    ho <- pop_summary$Hobs
    he <- pop_summary$Hexp
    
    # Wilcoxon signed-rank test for Heterozygosity Excess
    wilcox_test <- wilcox.test(ho, he, paired = TRUE, alternative = "greater")
    
    k_per_locus <- pop_summary$loc.n.all
    mean_k <- mean(k_per_locus, na.rm = TRUE)
    
    data.frame(
      Taxon = sp_name,
      N_Ind = n_ind,
      N_Loci = n_loci,
      Mean_Ho = round(mean(ho, na.rm = TRUE), 4),
      Mean_He = round(mean(he, na.rm = TRUE), 4),
      Mean_Alleles_Per_Locus = round(mean_k, 2),
      Wilcox_V = wilcox_test$statistic,
      Bottleneck_p_val = format.pval(wilcox_test$p.value, digits = 4),
      Bottleneck_Signal = ifelse(wilcox_test$p.value < 0.05, "Significant Excess", "Neutral/Deficit")
    )
  })
  
  return(do.call(rbind, results))
}

bottleneck_results <- run_bottleneck_analysis(gid_full_split)

print(bottleneck_results)

# ------------------------------------------------------------------------------
# 3. SELECTION TESTS (Run on Unpruned Dataset: gid_bent_unpruned)
# ------------------------------------------------------------------------------
# ==============================================================================
# Unified inbreedR-Native HWE Permutation Test & Plotting
# ==============================================================================

library(adegenet)
library(inbreedR)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# 1. NATIVE inbreedR PERMUTATION TEST FUNCTION
# ------------------------------------------------------------------------------
mlh_hwe_inbreedr_test <- function(gid, n_perm = 1000) {
  # Convert genind object directly to inbreedR raw 2-column format per locus
  # Extracts genotypes without altering missing data
  raw_mat <- tab(gid, freq = FALSE)
  
  # Calculate exact observed MLH per individual via inbreedR
  # inbreedR::MLH expects 0/1/2 dosage or two-column genind matrix
  ind_obs_mlh <- inbreedR::MLH(raw_mat)
  obs_mean_mlh <- mean(ind_obs_mlh, na.rm = TRUE)
  
  # Extract locus dosages (0, 1, 2) to calculate accurate allele frequencies
  locus_cols <- seq(1, ncol(raw_mat), by = 2)
  dosage_mat <- raw_mat[, locus_cols]
  
  # Calculate minor allele frequency p per locus
  p_vec <- colMeans(dosage_mat, na.rm = TRUE) / 2
  valid_p <- p_vec[!is.na(p_vec)]
  
  n_ind <- nInd(gid)
  
  # Simulate HWE null genotypes (0, 1, 2) directly from locus p-frequencies
  null_mlhs <- replicate(n_perm, {
    sim_genotypes <- sapply(valid_p, function(p) rbinom(n_ind, size = 2, prob = p))
    # Dosage of 1 corresponds to heterozygous locus
    mean(rowMeans(sim_genotypes == 1, na.rm = TRUE))
  })
  
  p_val <- sum(null_mlhs >= obs_mean_mlh) / n_perm
  
  return(list(
    summary = data.frame(
      Observed_MLH = round(obs_mean_mlh, 4),
      Expected_HWE_MLH = round(mean(null_mlhs), 4),
      p_value = round(p_val, 4)
    ),
    ind_mlh = ind_obs_mlh
  ))
}

# ------------------------------------------------------------------------------
# 2. RUN PIPELINE & GENERATE SUMMARY TABLE
# ------------------------------------------------------------------------------

mlh_results_list <- list()
mlh_plot_df <- data.frame()

for (sp in names(gid_unpruned_split)) {
  res <- mlh_hwe_inbreedr_test(gid_unpruned_split[[sp]], n_perm = 1000)
  
  # Store summary row
  sum_row <- res$summary
  sum_row$Taxon <- sp
  mlh_results_list[[sp]] <- sum_row[, c("Taxon", "Observed_MLH", "Expected_HWE_MLH", "p_value")]
  
  # Store individual MLH values for plotting
  sp_df <- data.frame(
    Taxon = sp,
    Individual = names(res$ind_mlh),
    MLH = as.numeric(res$ind_mlh)
  )
  mlh_plot_df <- rbind(mlh_plot_df, sp_df)
}

unified_mlh_table <- do.call(rbind, mlh_results_list)
print("--- UNIFIED inbreedR OVERDOMINANCE RESULTS ---")
print(unified_mlh_table)

# ------------------------------------------------------------------------------
# 3. PLOT INDIVIDUAL MLH DISTRIBUTIONS (inbreedR)
# ------------------------------------------------------------------------------
p_mlh <- ggplot(mlh_plot_df, aes(x = MLH, fill = Taxon)) +
  geom_histogram(bins = 15, color = "black", alpha = 0.75, boundary = 0) +
  facet_wrap(~Taxon, scales = "free_y") +
  scale_fill_viridis_d(option = "mako") +
  labs(
    title = "Individual Multilocus Heterozygosity (inbreedR Framework)",
    subtitle = "Calculated across 52,718 unpruned loci; right shift indicates overdominance",
    x = "Individual Multilocus Heterozygosity (MLH)",
    y = "Number of Individuals"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold")
  )

print(p_mlh)

