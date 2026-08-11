library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)
library(readr)

# ============================================================
# 1. Load Tier 1 and Tier 2 data only
# ============================================================

tier1 <- read_rds("~/Phase 2/tier1.rds")
tier2 <- read_rds("~/Phase 2/tier2.rds")

# IMPORTANT:
# Create a new tier_label so it definitely says Tier 1 / Tier 2
# Do not use all_results here
tier12_df <- bind_rows(
  tier1 %>% mutate(tier_label = "Tier 1"),
  tier2 %>% mutate(tier_label = "Tier 2")
) %>%
  mutate(
    pair_plot_label = paste0(
      exposure_name,
      " -> ",
      outcome_name,
      " (n = ", n_snps,
      ", ", tier_label, ")"
    )
  )

# Check number of unique pairs
tier12_df %>%
  distinct(ancestry, exposure_name, outcome_name, tier_label) %>%
  count(ancestry, tier_label)

# ============================================================
# 2. Convert IVW, MR-Egger, Weighted median into long format
# ============================================================

plot_df <- tier12_df %>%
  select(
    ancestry,
    exposure_name,
    outcome_name,
    pair_plot_label,
    tier_label,
    n_snps,
    ivw_beta, ivw_se, ivw_pval,
    egger_beta, egger_se, egger_pval,
    wm_beta, wm_se, wm_pval
  ) %>%
  pivot_longer(
    cols = c(
      ivw_beta, ivw_se, ivw_pval,
      egger_beta, egger_se, egger_pval,
      wm_beta, wm_se, wm_pval
    ),
    names_to = c("method", ".value"),
    names_pattern = "(ivw|egger|wm)_(beta|se|pval)"
  ) %>%
  mutate(
    method = recode(
      method,
      "ivw" = "IVW",
      "egger" = "MR-Egger",
      "wm" = "Weighted median"
    ),
    
    # Reversed for plot dodge order
    # Top to bottom visually will be: IVW, MR-Egger, Weighted median
    method = factor(
      method,
      levels = c("Weighted median", "MR-Egger", "IVW")
    ),
    
    OR = exp(beta),
    OR_lower = exp(beta - 1.96 * se),
    OR_upper = exp(beta + 1.96 * se)
  ) %>%
  filter(
    !is.na(beta),
    !is.na(se),
    !is.na(pval),
    !is.na(ancestry)
  )

# ============================================================
# 3. Forest plot function
# ============================================================

make_tier_forest_plot <- function(data, ancestry_name, x_limits, x_breaks) {
  
  plot_data <- data %>%
    filter(ancestry == ancestry_name) %>%
    mutate(
      tier_label = factor(tier_label, levels = c("Tier 1", "Tier 2"))
    ) %>%
    arrange(tier_label, exposure_name, outcome_name) %>%
    mutate(
      pair_plot_label = factor(
        pair_plot_label,
        levels = rev(unique(pair_plot_label))
      )
    )
  
  ggplot(
    plot_data,
    aes(
      x = OR,
      y = pair_plot_label,
      xmin = OR_lower,
      xmax = OR_upper,
      colour = method,
      shape = method,
      group = method
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.5
    ) +
    geom_errorbarh(
      position = position_dodge(width = 0.65),
      height = 0.25,
      linewidth = 0.6
    ) +
    geom_point(
      position = position_dodge(width = 0.65),
      size = 2.5
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_breaks
    ) +
    coord_cartesian(
      xlim = x_limits
    ) +
    scale_colour_manual(
      breaks = c("IVW", "MR-Egger", "Weighted median"),
      values = c(
        "IVW" = "#d73027",
        "MR-Egger" = "#2b6cb0",
        "Weighted median" = "#1a9850"
      )
    ) +
    scale_shape_manual(
      breaks = c("IVW", "MR-Egger", "Weighted median"),
      values = c(
        "IVW" = 15,
        "MR-Egger" = 16,
        "Weighted median" = 17
      )
    ) +
    labs(
      title = paste0("Tier 1 and Tier 2 MR results: ", ancestry_name),
      subtitle = "Tier 1 pairs are shown first, followed by Tier 2 pairs",
      x = "Odds Ratio (95% CI)",
      y = NULL,
      colour = NULL,
      shape = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey30"),
      legend.position = "bottom",
      legend.text = element_text(size = 9)
    )
}

# ============================================================
# 4. Generate EUR and EAS plots
# ============================================================

eur_tier_plot <- make_tier_forest_plot(
  plot_df,
  ancestry_name = "EUR",
  x_limits = c(0, 5),
  x_breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4,4.5, 5)
)

eas_tier_plot <- make_tier_forest_plot(
  plot_df,
  ancestry_name = "EAS",
  x_limits = c(0, 5),
  x_breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4,4.5, 5)
)

eur_tier_plot
eas_tier_plot

# ============================================================
# 5. Save plots
# ============================================================

ggsave(
  "~/Phase 1/Results/Tier1_Tier2_forest_EUR_IVW_Egger_WM.png",
  eur_tier_plot,
  width = 13,
  height = 10,
  units = "in",
  dpi = 300
)

ggsave(
  "~/Phase 1/Results/Tier1_Tier2_forest_EAS_IVW_Egger_WM.png",
  eas_tier_plot,
  width = 13,
  height = 10,
  units = "in",
  dpi = 300
)