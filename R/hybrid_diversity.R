### Explore relationship between hybrid method performance and method diversity
# hybrid_diversity.R: explore relationship between hybrid method performance and method diversity
library(stringr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

methods <- c("mv_lognc_scran(1mv)*",
             "logmv_lognc(2lmv)",
             "poisson_ct_scran(3pos)*",
             "disp_nc_seuratv1(4dis)*",
             "mean_max_nc(5max)",
             "1mv2lmv",
             "1mv3pos",
             "1mv4dis",
             "1mv5max",
             "2lmv3pos",
             "2lmv4dis",
             "2lmv5max",
             "3pos4dis",
             "3pos5max",
             "4dis5max",
             "1mv2lmv3pos",
             "1mv2lmv4dis",
             "1mv2lmv5max",
             "1mv3pos4dis",
             "1mv3pos5max",
             "1mv4dis5max",
             "2lmv3pos4dis",
             "2lmv3pos5max",
             "2lmv4dis5max",
             "3pos4dis5max",
             "1mv2lmv3pos4dis",
             "1mv2lmv3pos5max",
             "1mv2lmv4dis5max",
             "1mv3pos4dis5max",
             "2lmv3pos4dis5max",
             "1mv2lmv3pos4dis5max")

hybrid_methods <- methods[6:length(methods)]

base_methods <- methods[1:5]

baseline_hybrid <- c("mv_lognc_scran(1mv)*",
                     "logmv_lognc(2lmv)",
                     "poisson_ct_scran(3pos)*",
                     "disp_nc_seuratv1(4dis)*",
                     "mean_max_nc(5max)")

names(baseline_hybrid) <- c("1mv", "2lmv", "3pos", "4dis", "5max")

# Read HVGs for each dataset and method
hvg_list <- readRDS("/dcs05/hongkai/data/yli6/rzhao/mixhvg/full_hvg_list.rds")

# Compute geneset diversity for hybrid methods
get_hybrid_diversity <- function(hvg_list, hybrid_method, score_method = "entropy") {
  # Get baseline methods
  base_methods <- baseline_hybrid[unlist(str_extract_all(hybrid_method, "\\d[a-z]+"))]
  diversity <- sapply(names(hvg_list), function(dataset) {
    if (score_method == "pool") {
      # diversity = number of distinct HVGs in all HVGs identified by constituent baseline methods
      return(n_distinct(unlist(hvg_list[[dataset]][base_methods])))
    } else if (score_method == "entropy") {
      # diversity = Shannon's entropy of baseline methods HVGs composition of method-specific hybrid HVGs
      hybrid_hvgs <- hvg_list[[dataset]][[hybrid_method]]
      hvg_membership <- sapply(base_methods, function(method) {
        hybrid_hvgs %in% hvg_list[[dataset]][[method]]
      })
      hvg_membership <- hvg_membership[rowSums(hvg_membership) == 1, ]
      method_sp_distr <- colSums(hvg_membership) / sum(hvg_membership)
      method_sp_distr <- method_sp_distr[method_sp_distr > 0]
      - sum(method_sp_distr * log2(method_sp_distr))
    }
  })
  return(mean(diversity))
}

hybrid_entropy <- sapply(hybrid_methods, function(hybrid_method) {
  get_hybrid_diversity(hvg_list, hybrid_method, score_method = "entropy")
})

# Inverse of average ranks (larger value indicates better performance)
ranks <- readRDS("/dcs05/hongkai/data/yli6/rzhao/mixhvg/inv_avg_ranks.rds")

# Compute average baseline methods ranks for each hybrid method
avg_base_ranks <- sapply(hybrid_methods, function(hybrid_method) {
  base_methods <- baseline_hybrid[unlist(str_extract_all(hybrid_method, "\\d[a-z]+"))]
  mean(ranks[base_methods])
})

# The number of constituent baseline methods in each hybrid method
nhybrid <- sapply(hybrid_methods, function(hybrid_method) {
  length(unlist(str_extract_all(hybrid_method, "\\d[a-z]+")))
})

### Explore effect of diversity and baseline method performance on hybrid method rank

regr_df <- data.frame(
  method = methods, 
  ranks = ranks[methods], 
  diversity = c(rep(1, 5), nhybrid), 
  entropy = c(rep(0, 5), hybrid_entropy), 
  avg_base_ranks = c(ranks[methods[1:5]], avg_base_ranks), 
  method_type = factor(c(rep("Baseline Method", 5), paste0(nhybrid, "Mix Hybrid")), 
                       levels = c("Baseline Method", paste0(2:5, "Mix Hybrid")))
)

# Plot hybrid method rank against diversity score (include baseline methods)
g1 <- ggplot(regr_df, aes(diversity, ranks)) + 
  geom_point(aes(color = method_type)) + 
  geom_smooth(method = "lm") +
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(arrow = arrow(type = "closed", length = unit(10, "pt"))), # Add closed arrows to both axes
    axis.title.x = element_text(face = "bold"), 
    axis.title.y = element_text(face = "bold"), 
    legend.position = "none"
  ) +
  labs(x = "Number of Baseline Methods Included", 
       y = "Overall Method Rank") +
  scale_color_manual(
    values = c("grey", "#FDC086", "#D95F02", "#7570B3", "#66A61E")
  )

# Get hybrid method pairs with and without 5max
meanmax_hybrid <- methods[grepl("*5max$", methods)]
hybrid_5max_df <- data.frame(
  method_wo_5max = str_replace(meanmax_hybrid[5:length(meanmax_hybrid)], "5max", ""), 
  method_w_5max = meanmax_hybrid[5:length(meanmax_hybrid)]
)
hybrid_5max_df <- hybrid_5max_df %>% 
  mutate(
    rank_wo_5max = ranks[method_wo_5max], 
    rank_w_5max = ranks[method_w_5max]
  )

# Plot 5max hybrid vs hybrid without 5max ranks
g2 <- ggplot(hybrid_5max_df, aes(rank_wo_5max, rank_w_5max)) +
  geom_abline(slope = 1, intercept = 0, color = "red", size = 1.5, linetype = "dashed") + 
  geom_point() + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(arrow = arrow(type = "closed", length = unit(10, "pt"))), # Add closed arrows to both axes
    axis.title.x = element_text(face = "bold"), 
    axis.title.y = element_text(face = "bold"), 
    legend.position = "none"
  ) +
  xlim(c(1.6, 2.5)) + 
  ylim(c(1.6, 2.5)) +
  labs(x = "Rank of Hybrid Method without 5max", 
       y = "Rank of Hybrid Method with 5max")

# Plot hybrid method rank against average base ranks
g3 <- ggplot(regr_df, aes(avg_base_ranks, ranks)) + 
  geom_point(aes(color = method_type)) + 
  geom_smooth(method = "lm") +
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(arrow = arrow(type = "closed", length = unit(10, "pt"))), # Add closed arrows to both axes
    axis.title.x = element_text(face = "bold"), 
    axis.title.y = element_text(face = "bold"), 
    legend.position = "none"
  ) + 
  labs(x = "Average Baseline Methods Rank", 
       y = "Overall Method Rank") +
  scale_color_manual(
    values = c("grey", "#FDC086", "#D95F02", "#7570B3", "#66A61E")
  )

# Plot hybrid method rank against entropy (include baseline methods)
g4 <- ggplot(regr_df, aes(entropy, ranks)) + 
  geom_point(aes(color = method_type)) + 
  geom_smooth(method = "lm") +
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(arrow = arrow(type = "closed", length = unit(10, "pt"))), # Add closed arrows to both axes
    axis.title.x = element_text(face = "bold"), 
    axis.title.y = element_text(face = "bold"), 
    legend.title = element_text(face = "bold"), 
    legend.text = element_text(face = "bold")
  ) +
  labs(color = "Method Type", 
       x = "Shannon's Entropy", 
       y = "Overall Method Rank") +
  scale_color_manual(
    values = c("grey", "#FDC086", "#D95F02", "#7570B3", "#66A61E")
  )

top_row <- plot_grid(g1, g2, NULL, ncol = 3, 
                           scale = 0.9, 
                           rel_widths = c(1.5, 1.5, 0.6), 
                           labels = c("a", "b", ""))
bot_row <- plot_grid(g3, g4, ncol = 2, 
                           scale = 0.9, 
                           rel_widths = c(1.5, 2.1), 
                           labels = c("c", "d"))
combined_plot <- plot_grid(top_row, bot_row, ncol = 1)
ggsave("rank_diversity_plot.pdf", plot = combined_plot, width = 250, height = 200, unit = "mm")

# Get correlations between overall rank and baseline ranks and entropy
cor.test(regr_df$ranks, regr_df$avg_base_ranks, method = "pearson")
cor.test(regr_df$ranks, regr_df$entropy, method = "pearson")

# Linear regression on overall ranks using HVG diversity
lm_res <- lm(ranks ~ diversity, data = regr_df)
summary(lm_res)

# Linear regression on overall ranks using HVG diversity and average baseline ranks
lm_res <- lm(ranks ~ diversity + avg_base_ranks, data = regr_df)
summary(lm_res)
confint(lm_res)

# Significant difference between 5max hybrids and non-5max hybrids
hybrid_ranks <- ranks[hybrid_methods]
meanmax_methods_ind <- grepl("5max", hybrid_methods)
t.test(hybrid_ranks[meanmax_methods_ind], hybrid_ranks[!meanmax_methods_ind], alternative = "less")

regr_df <- regr_df %>% mutate(contains_5max = ifelse(grepl("5max", method), "5max", "no 5max"))

lm_res <- lm(ranks ~ diversity, data = regr_df %>% filter(contains_5max == "5max"))
summary(lm_res)

lm_res <- lm(ranks ~ diversity + avg_base_ranks, data = regr_df %>% filter(contains_5max == "5max"))
summary(lm_res)

lm_res <- lm(ranks ~ diversity, data = regr_df %>% filter(contains_5max == "no 5max"))
summary(lm_res)

lm_res <- lm(ranks ~ diversity + avg_base_ranks, data = regr_df %>% filter(contains_5max == "no 5max"))
summary(lm_res)
confint(lm_res)

# Linear regression using Shannon's entropy
lm_res <- lm(ranks ~ entropy + avg_base_ranks, data = regr_df)
summary(lm_res)
confint(lm_res)
lm_res <- lm(ranks ~ entropy + avg_base_ranks, data = regr_df %>% filter(contains_5max == "no 5max"))
summary(lm_res)
confint(lm_res)
