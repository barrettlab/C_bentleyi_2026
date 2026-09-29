#########################################################

# Set up simulations to test empirical data vs OD and NFDS

# Selection mode: OD, NFDS
# OD fintesses: AA, Aa, aa (1,1,1 neutral; 0.85,1,0.85 weak; 0.5,1,0.5 moderate; 0.15.1,0.15 strong; 0,1,0 lethal)
# NFDS fitnesses: 0.3, 0.8
# N: 10, 100
# Selfing: 90, 100%


#########################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(transport)
library(furrr)

# 1. Define # cored
plan(multisession, workers = availableCores() - 1)

# 2. Set up for running in parallel
simulate_factorial_fis <- function(n_loci = 500, N = 100, selfing_rate = 0.90, 
                                   generations = 50, mechanism = "OD", 
                                   selection_level = "weak") {
  fis_results <- numeric(n_loci)
  
  for (l in 1:n_loci) {
    # Define genotypes: 1 = AA, 2 = Aa, 3 = aa
    pop <- sample(c(1, 2, 3), size = N, replace = TRUE, prob = c(0.25, 0.5, 0.25))
    
    for (gen in 1:generations) {
      n_AA <- sum(pop == 1); n_Aa <- sum(pop == 2); n_aa <- sum(pop == 3)
      p <- (2 * n_AA + n_Aa) / (2 * N)
      if (p == 0 || p == 1) break
      
      # Assign Fitnesses
      if (mechanism == "OD") {
        if (selection_level == "neutral")    w <- c(1.00, 1.0, 1.00) # True Neutral Baseline
        if (selection_level == "weak")       w <- c(0.85, 1.0, 0.85)
        if (selection_level == "moderate")   w <- c(0.50, 1.0, 0.50)
        if (selection_level == "strong")     w <- c(0.15, 1.0, 0.15)
        if (selection_level == "lethal")     w <- c(0.00, 1.0, 0.00) # Balanced Lethal
      } else if (mechanism == "NFDS") {
        freqs <- c(n_AA, n_Aa, n_aa) / N
        s_coeff <- ifelse(selection_level == "weak", 0.3, 0.8)
        w <- 1 - (s_coeff * freqs)
        w[w < 0] <- 0.01
      }
      
      # Probability weight for parent sampling
      prob_weights <- w * c(n_AA, n_Aa, n_aa)
      
      # Safety check: If all individuals have 0 fitness (e.g. loss of heterozygotes), break early
      if (sum(prob_weights) == 0) break
      
      next_pop <- numeric(N)
      for (i in 1:N) {
        if (runif(1) < selfing_rate) {
          # Selfing branch
          parent <- sample(c(1, 2, 3), 1, prob = prob_weights)
          if (parent == 1) next_pop[i] <- 1
          else if (parent == 3) next_pop[i] <- 3
          else next_pop[i] <- sample(c(1, 2, 3), 1, prob = c(0.25, 0.5, 0.25))
        } else {
          # Outcrossing branch
          p_fit <- (2 * n_AA * w[1] + n_Aa * w[2]) / (2 * sum(prob_weights))
          g1 <- runif(1) < p_fit; g2 <- runif(1) < p_fit
          if (g1 && g2) next_pop[i] <- 1
          else if (!g1 && !g2) next_pop[i] <- 3
          else next_pop[i] <- 2
        }
      }
      pop <- next_pop
    }
    
    # Calculate final locus FIS
    n_AA <- sum(pop == 1); n_Aa <- sum(pop == 2); n_aa <- sum(pop == 3)
    p <- (2 * n_AA + n_Aa) / (2 * N)
    He <- 2 * p * (1 - p)
    Ho <- n_Aa / N
    if (He > 0) fis_results[l] <- (He - Ho) / He else fis_results[l] <- NA
  }
  
  return(data.frame(
    FIS = fis_results[!is.na(fis_results)],
    Selfing = paste0(selfing_rate * 100, "%"),
    PopSize = paste0("N = ", N),
    Mechanism = mechanism,
    Level = selection_level
  ))
}

# 3. Create Parameter Set (7 selection sets x 2 selfing rates x 2 pop sizes = 28 combinations)
param_grid <- expand.grid(
  selfing_rate = c(0.90, 1.0),
  N = c(10, 100),
  mechanism_level = c("OD_neutral", "OD_weak", "OD_moderate", "OD_strong", "OD_lethal", 
                      "NFDS_weak", "NFDS_strong"),
  stringsAsFactors = FALSE
) %>%
  separate(mechanism_level, into = c("mechanism", "selection_level"), sep = "_")

# 4. Run Parallel Simulations

set.seed(123)

df_all_sims_popsize <- future_pmap_dfr(
  list(
    selfing_rate = param_grid$selfing_rate,
    N = param_grid$N,
    mechanism = param_grid$mechanism,
    selection_level = param_grid$selection_level
  ),
  function(selfing_rate, N, mechanism, selection_level) {
    simulate_factorial_fis(
      n_loci = 500, 
      N = N, 
      selfing_rate = selfing_rate, 
      generations = 50, 
      mechanism = mechanism, 
      selection_level = selection_level
    )
  },
  .options = furrr_options(seed = TRUE)
) %>%
  mutate(Model_ID = paste(Mechanism, Level, sep = "_"))

plan(sequential) # Close background workers

# 5. Fit Model Rankings (Wasserstein Distance & KS Test)
fit_rankings_popsize <- df_all_sims_popsize %>%
  group_by(PopSize, Selfing, Mechanism, Level, Model_ID) %>%
  summarise(
    KS_Stat = ks.test(FIS, emp_fis)$statistic,
    Wasserstein_Dist = wasserstein1d(FIS, emp_fis),
    Mean_FIS = mean(FIS),
    .groups = "drop"
  ) %>%
  arrange(Wasserstein_Dist)

print(head(as.data.frame(fit_rankings_popsize), 10))


# 6. Summarize model weights

library(dplyr)

# Calculate  Model Weights across the 28  combinations
model_weights_kernel <- fit_rankings_popsize %>%
  mutate(
    # Set bandwidth based on the best-fitting model distance
    epsilon = min(Wasserstein_Dist) * 1.5,
    
    # Unnormalized Gaussian Likelihood
    Raw_Likelihood = exp(-0.5 * (Wasserstein_Dist / epsilon)^2)
  ) %>%
  mutate(
    # Normalize to sum to 1.0 (Model Weights)
    Model_Weight = Raw_Likelihood / sum(Raw_Likelihood)
  ) %>%
  arrange(desc(Model_Weight)) %>%
  select(PopSize, Selfing, Mechanism, Level, Model_ID, Wasserstein_Dist, Model_Weight)

# Display top models
print(as.data.frame(head(model_weights_kernel, 10)))


################################################################

# Plotting

################################################################

library(ggplot2)
library(dplyr)

# Panel 1: Updated Density Plot with Reversed Aesthetics
p_density <- ggplot() +
  # Empirical Data Layer (Purple Fill)
  geom_density(
    data = data.frame(FIS = emp_fis), 
    aes(x = FIS, fill = "Empirical Data"), 
    color = "purple4", linewidth = 1.1, alpha = 0.35
  ) +
  # All Simulated Models
  # Linetype = Mechanism (OD = solid, NFDS = dotted)
  # Color = Selection Intensity Level
  geom_density(
    data = df_all_sims_popsize, 
    aes(x = FIS, linetype = Mechanism, color = Level, group = Model_ID), 
    linewidth = 0.9, alpha = 0.85
  ) +
  facet_grid(PopSize ~ Selfing, labeller = label_both) +
  coord_cartesian(ylim = c(0, 6), xlim = c(-1.0, 1.0)) +
  
  # Legend Customization
  scale_fill_manual(name = "", values = c("Empirical Data" = "mediumorchid")) +
  
  # Mechanism Aesthetics: OD (solid) vs NFDS (dotted)
  scale_linetype_manual(
    name = "Mechanism", 
    values = c("OD" = "solid", "NFDS" = "dotted")
  ) +
  
  # Selection Intensity Palette (Distinct Color per Tier)
  scale_color_manual(
    name = "Selection Tier", 
    values = c(
      "lethal"   = "#00441B", # Dark Green (Top Performer)
      "strong"   = "#1B9E77", # Teal Green
      "moderate" = "#D95F02", # Orange
      "weak"     = "#E7298A", # Magenta/Pink
      "neutral"  = "#7570B3"  # Purple/Grey
    ),
    breaks = c("lethal", "strong", "moderate", "weak", "neutral")
  ) +
  
  theme_bw(base_size = 11) +
  labs(
    title = "Model Fitting: Population Size (N) vs. Selfing Rate",
    x = "Locus-Specific Inbreeding Coefficient (FIS)",
    y = "Density"
  ) +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold")
  )

print(p_density)

# Panel 2: Heatmap Plot
p_heatmap <- ggplot(df_heatmap, aes(x = Model_Label, y = Demographic_Scenario, fill = Model_Weight)) +
    geom_tile(color = "white", linewidth = 0.5) +
    # Overlay percentages inside cells
    geom_text(
        aes(label = ifelse(Model_Weight > 0.005, sprintf("%.1f%%", Model_Weight * 100), "<0.5%")), 
        size = 3, 
        color = ifelse(df_heatmap$Model_Weight > 0.10, "white", "black")
    ) +
    scale_fill_gradient(
        low = "white", 
        high = "black", 
        labels = scales::percent, 
        name = "Weight"
    ) +
    theme_bw(base_size = 11) +
    labs(
        title = "B) Relative Model Weights Grid",
        x = "Selection Model",
        y = "Demographic Context"
    ) +
    theme(
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
        panel.grid = element_blank(),
        legend.position = "right"
    )

