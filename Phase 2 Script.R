library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(purrr)

## ---- 1. Combine significant pairs ----
tier1 <- readRDS("~/Phase 2/tier1.rds")
tier2 <- readRDS("~/Phase 2/tier2.rds")
sig_pairs <- bind_rows(
  tier1 %>% mutate(tier = "Tier 1"),
  tier2 %>% mutate(tier = "Tier 2")
) %>%
  select(exposure, outcome, ancestry, tier, ivw_beta, ivw_pval, fdr_qval) %>%
  distinct()

## ---- 2. Build graph + metrics + cascades for one ancestry ----
build_mr_network <- function(data, ancestry_label) {
  
  edges <- data %>% filter(ancestry == ancestry_label)
  
  if (nrow(edges) == 0) {
    message("No edges for ", ancestry_label)
    return(NULL)
  }
  
  g <- graph_from_data_frame(edges, directed = TRUE)
  E(g)$sign  <- ifelse(E(g)$ivw_beta > 0, "positive", "negative")
  E(g)$width <- -log10(E(g)$ivw_pval)
  
  ## topology metrics
  node_metrics <- tibble(
    node        = V(g)$name,
    in_degree   = degree(g, mode = "in"),
    out_degree  = degree(g, mode = "out"),
    betweenness = betweenness(g, directed = TRUE, normalized = TRUE)
  ) %>% arrange(desc(betweenness))
  
  ## cascades: directed paths length >= 2 edges
  node_names <- V(g)$name
  cascades <- list()
  for (src in node_names) {
    for (tgt in node_names) {
      if (src == tgt) next
      paths <- all_simple_paths(g, from = src, to = tgt, mode = "out")
      long_paths <- Filter(function(p) length(p) >= 3, paths)
      if (length(long_paths) > 0) cascades[[paste(src, tgt, sep = "->")]] <- long_paths
    }
  }
  cascade_df <- map_dfr(names(cascades), function(k) {
    map_dfr(cascades[[k]], function(p) {
      tibble(path = paste(V(g)$name[p], collapse = " -> "), n_steps = length(p) - 1)
    })
  })
  
  ## --- fix: guarantee columns exist even with zero cascades ---
  if (nrow(cascade_df) == 0) {
    cascade_df <- tibble(path = character(0), n_steps = numeric(0))
  }
  cascade_df <- distinct(cascade_df)
  
  list(graph = g, node_metrics = node_metrics, cascades = cascade_df, ancestry = ancestry_label)
}

eur_net <- build_mr_network(sig_pairs, "EUR")
eas_net <- build_mr_network(sig_pairs, "EAS")

## quick look
eur_net$node_metrics
eas_net$node_metrics
eur_net$cascades %>% arrange(desc(n_steps))
eas_net$cascades %>% arrange(desc(n_steps))

## ---- 3. Plotting function ----
plot_mr_network <- function(net) {
  set.seed(123)
  if (is.null(net)) return(NULL)
  
  g <- net$graph
  
  ggraph(g, layout = "fr") +
    geom_edge_link(
      aes(color = sign, width = width),
      arrow = arrow(length = unit(3, "mm"), type = "closed"),
      end_cap = circle(3, "mm"),
      alpha = 0.7
    ) +
    scale_edge_color_manual(
      values = c(
        positive = "#d73027",
        negative = "#4575b4"
      )
    ) +
    scale_edge_width(range = c(0.4, 2.5)) +
    
    geom_node_point(
      aes(size = degree(g, mode = "all")),
      color = "steelblue"
    ) +
    
    # Bigger node labels
    geom_node_text(
      aes(label = name),
      repel = TRUE,
      size = 5,
      fontface = "bold"
    ) +
    
    theme_graph(base_family = "sans") +
    labs(
      title = paste("MR Causal Network -", net$ancestry),
      edge_color = "Effect direction",
      edge_width = "-log10(p)"
    ) +
    theme(
      plot.title = element_text(size = 18, face = "bold")
    )
}
p_eur <- plot_mr_network(eur_net)
p_eas <- plot_mr_network(eas_net)

p_eur
p_eas

ggsave("~/Phase 2/MR_network_EUR.png", p_eur, width = 12, height = 9)
ggsave("~/Phase 2/MR_network_EAS.png", p_eas, width = 12, height = 9)