## R code for plotting C bentleyi floral form data (chasmogamous vs cleistogamous)

# Load required libraries
library(tidyverse)
library(patchwork)

# Load data
df <- read.csv("2022_bentleyi_cleistogamous.csv")

# Set consistent color palette across plots
pop_colors <- c("CG" = "#2b5c8f", "DK" = "#e07a5f", "PC" = "#588b76", "PM" = "#d4a373")

# Custom theme for publication consistency
theme_pub <- function() {
  theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 9, color = "black"),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.background = element_rect(fill = alpha("white", 0.8), color = "gray80"),
      panel.grid.major.y = element_line(color = "gray92", linewidth = 0.3)
    )
}

# ------------------------------------------------------------------------------
# Panel A: Proportion Open (Chasmogamous) Flowers by Site
# ------------------------------------------------------------------------------
p_a <- ggplot(df, aes(x = pop, y = prop_chasm, fill = pop)) +
  geom_boxplot(alpha = 0.75, width = 0.45, outlier.shape = NA) +
  geom_jitter(color = "black", alpha = 0.6, width = 0.15, size = 1.8) +
  scale_fill_manual(values = pop_colors) +
  scale_y_continuous(limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "A. Chasmogamous Proportion by Site",
    x = "Population",
    y = "Proportion Open (Chasmogamous)"
  ) +
  theme_pub() +
  theme(legend.position = "none")

# ------------------------------------------------------------------------------
# Panel B: Interannual Plasticity (PC and PM across years)
# ------------------------------------------------------------------------------
df_multi <- df %>% filter(pop %in% c("PC", "PM"))

p_b <- ggplot(df_multi, aes(x = factor(year), y = prop_chasm, fill = pop)) +
  geom_boxplot(alpha = 0.75, width = 0.55, outlier.shape = NA, position = position_dodge(0.7)) +
  geom_point(aes(group = pop), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), 
             color = "black", alpha = 0.6, size = 1.5) +
  scale_fill_manual(values = pop_colors, name = "Site") +
  scale_y_continuous(limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "B. Interannual Plasticity (PC & PM)",
    x = "Year",
    y = "Proportion Open (Chasmogamous)"
  ) +
  theme_pub() +
  theme(
    legend.position = c(0.82, 0.82),
    legend.margin = margin(3, 5, 3, 5)
  )

# ------------------------------------------------------------------------------
# Panel C: Total Floral Output vs Strategy Trajectory
# ------------------------------------------------------------------------------
p_c <- ggplot(df, aes(x = total_flr, y = prop_chasm)) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40", linetype = "dashed", fill = "gray85") +
  geom_point(aes(color = pop, shape = pop), size = 2.5, alpha = 0.85) +
  scale_color_manual(values = pop_colors, name = "Site") +
  scale_shape_manual(values = c("CG" = 17, "DK" = 15, "PC" = 16, "PM" = 18), name = "Site") +
  scale_y_continuous(limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "C. Floral Output vs Display Strategy",
    x = "Total Flowers per Individual",
    y = "Proportion Open (Chasmogamous)"
  ) +
  theme_pub() +
  theme(
    legend.position = c(0.82, 0.78),
    legend.margin = margin(3, 5, 3, 5)
  )

# ------------------------------------------------------------------------------
# Panel D: Stacked Bar Chart of Absolute Flower Counts
# ------------------------------------------------------------------------------
df_counts <- df %>%
  group_by(pop) %>%
  summarise(
    Cleistogamous = sum(cleist),
    Chasmogamous = sum(chasm),
    Total = sum(total_flr),
    Prop_Open = Chasmogamous / Total
  ) %>%
  pivot_longer(cols = c("Cleistogamous", "Chasmogamous"), names_to = "Type", values_to = "Count") %>%
  mutate(
    pop = factor(pop, levels = c("CG", "DK", "PC", "PM")),
    Type = factor(Type, levels = c("Chasmogamous", "Cleistogamous"))
  )

# Extract summary labels
df_labels <- df_counts %>%
  group_by(pop) %>%
  summarise(
    Total = max(Total),
    Prop_Open = max(Prop_Open)
  ) %>%
  mutate(Label = sprintf("N=%d\n(%.1f%% Open)", Total, Prop_Open * 100))

p_d <- ggplot(df_counts, aes(x = pop, y = Count, fill = Type)) +
  geom_bar(stat = "identity", width = 0.45, alpha = 0.9) +
  geom_text(data = df_labels, aes(x = pop, y = Total + 18, label = Label), 
            inherit.aes = FALSE, size = 2.8, vjust = 0, lineheight = 0.9) +
  scale_fill_manual(values = c("Cleistogamous" = "#3d405b", "Chasmogamous" = "#e9c46a"), name = NULL) +
  scale_y_continuous(limits = c(0, 560), expand = c(0, 0)) +
  labs(
    title = "D. Cumulative Flower Production",
    x = "Population",
    y = "Total Flower Count Observed"
  ) +
  theme_pub() +
  theme(
    legend.position = c(0.32, 0.85),
    legend.margin = margin(2, 4, 2, 4)
  )

# ------------------------------------------------------------------------------
# Combine all 4 panels into a single multi-pane figure
# ------------------------------------------------------------------------------
final_figure <- (p_a | p_b) / (p_c | p_d)

# Save high-resolution PNG for manuscript submission
ggsave("cleistogamy_summary_figure_R.png", plot = final_figure, width = 10, height = 8, dpi = 300)