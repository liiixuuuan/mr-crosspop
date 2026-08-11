library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)

# ============================================================
# 1. Prepare IVW-only plotting data
# ============================================================

all_results = read_rds('~/Phase 2/all_results.rds')
n_tests <- 182
bonf_threshold <- 0.05 / n_tests

plot_df <- all_results %>%
  mutate(
    beta = ivw_beta,
    se = ivw_se,
    pval = ivw_pval,
    nsnp = n_snps,
    
    OR = exp(beta),
    OR_lower = exp(beta - 1.96 * se),
    OR_upper = exp(beta + 1.96 * se),
    
    sig_group = case_when(
      ivw_pval < bonf_threshold ~ "Significant",
      ivw_pval >= bonf_threshold & fdr_qval < 0.05 ~ "Nominal",
      TRUE ~ "Not significant"
    ),
    
    sig_group = factor(
      sig_group,
      levels = c("Not significant", "Nominal", "Significant")
    ),
    
    outcome_label = paste0(outcome_name, "  (n=", nsnp, ")"),
    
    # Use only trait name, remove dataset ID
    exposure_label = str_trim(str_remove(exposure_name, "\\s*\\|\\|.*$"))
  ) %>%
  filter(
    !is.na(beta),
    !is.na(se),
    !is.na(pval),
    !is.na(ancestry)
  )


make_forest_plot <- function(data, ancestry_name, x_limits, x_breaks) {
  
  plot_data <- data %>%
    filter(ancestry == ancestry_name) %>%
    arrange(exposure_label, outcome_name) %>%
    group_by(exposure_label) %>%
    mutate(outcome_label = fct_inorder(outcome_label)) %>%
    ungroup()
  
  ggplot(
    plot_data,
    aes(
      x = OR,
      y = outcome_label,
      xmin = OR_lower,
      xmax = OR_upper,
      colour = sig_group,
      shape = sig_group
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.5
    ) +
    geom_errorbarh(
      height = 0.18,
      linewidth = 0.6
    ) +
    geom_point(size = 2.4) +
    
    # Log-scale OR axis to show small OR values better
    scale_x_log10(
      breaks = x_breaks,
      labels = x_breaks
    ) +
    coord_cartesian(
      xlim = x_limits
    ) +
    
    facet_grid(
      exposure_label ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    scale_colour_manual(
      values = c(
        "Not significant" = "grey75",
        "Nominal" = "#2b6cb0",
        "Significant" = "#d73027"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Not significant" = 1,
        "Nominal" = 16,
        "Significant" = 15
      )
    ) +
    labs(
      title = paste0("Overall UVMR landscape: ", ancestry_name),
      subtitle = paste0(
        "IVW results only; Bonferroni p < ",
        signif(bonf_threshold, 3),
        "; FDR-BH q < 0.05."
      ),
      x = "Odds Ratio (95% CI, log scale)",
      y = NULL,
      colour = NULL,
      shape = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.placement = "outside",
      strip.background.y = element_rect(fill = "grey90", colour = NA),
      strip.text.y.left = element_text(
        angle = 0,
        face = "bold",
        size = 11,
        hjust = 0.5
      ),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey30"),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      panel.spacing.y = unit(0.25, "lines")
    )
}

eur_plot <- make_forest_plot(
  plot_df,
  ancestry_name = "EUR",
  x_limits = c(0.4, 15),
  x_breaks = c(0.5, 1, 2, 5, 10, 15)
)
eas_plot <- make_forest_plot(
  plot_df,
  ancestry_name = "EAS",
  x_limits = c(0.4, 30),
  x_breaks = c(0.5, 1, 2, 5, 10, 15, 30)
)

eur_plot
eas_plot

ggsave(
  "UVMR_forest_EUR_IVW_all_pairs.png",
  eur_plot,
  width = 11,
  height = 18,
  units = "in"
)

ggsave(
  "UVMR_forest_EAS_IVW_all_pairs.png",
  eas_plot,
  width = 11,
  height = 18,
  units = "in"
)