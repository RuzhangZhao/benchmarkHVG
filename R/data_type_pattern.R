### Explore mean-max methods good performance on Cell Sorting datasets
# data_type_pattern.R: explore baseline methods performance difference on 3 data types
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)
library(cowplot)

marker_dir <- "/dcs05/hongkai/data/yli6/rzhao/mixhvg/marker"

dataset_name_range = c("duo4",
                       "duo4un",
                       "duo8",
                       "GBM_sd",
                       "human",
                       "mouse",
                       "zheng", 
                       "cbmc_pbmc", #1
                       "cbmc8k",#2
                       "fetalBM",#3
                       "CD34",#4
                       "bmcite",#5
                       "seurat_cite",#6
                       "Sucovid", 
                       "pbmc3k",
                       "human_brain_3k",
                       "mouse_brain_fresh_5k",
                       "pbmc10k",
                       "lymphoma_14k")

cell_sorting_ds <- c("duo4",
                     "duo4un",
                     "duo8",
                     "GBM_sd",
                     "human",
                     "mouse",
                     "zheng")

cite_ds = c("cbmc_pbmc", #1
            "cbmc8k",#2
            "fetalBM",#3
            "CD34",#4
            "bmcite",#5
            "seurat_cite",#6
            "Sucovid")  #7

mult_ds <- c("pbmc3k",
             "human_brain_3k",
             "mouse_brain_fresh_5k",
             "pbmc10k",
             "lymphoma_14k")

dataset_info <- data.frame(dataset = dataset_name_range)
dataset_info <- dataset_info %>% 
  mutate(data_type = case_when(dataset %in% cell_sorting_ds ~ "Cell Sorting",
                               dataset %in% cite_ds ~ "CITE-seq",
                               dataset %in% mult_ds ~ "MultiomeATAC"))

methods <- c("random",
             "mv_ct",
             "mv_nc",
             "mv_lognc_scran(1mv)*",
             "mv_PFlogPF",
             "logmv_ct_seuratv3*",
             "logmv_nc",
             "logmv_lognc(2lmv)",
             "logmv_PFlogPF",
             "poisson_ct_scran(3pos)*",
             "disp_ct",
             "disp_nc_seuratv1(4dis)*",
             "mvp_nc_seuratv2*",
             "disp_lognc",
             "disp_PFlogPF",
             "scanpy_cell_ranger*",
             "mean_max_ct",
             "mean_max_nc(5max)",
             "mean_max_lognc",
             "mean_max_PFlogPF",
             "SCT*")

method_info <- data.frame(
  method = methods, 
  method_type = factor(c("random", rep("var~mean(mv)", 4), 
                         rep("logvar~logmean", 4), "poisson_scran(pos)", 
                         rep("dispersion(disp)", 6), 
                         rep("mean_max(max)", 4), "SCT"), 
                       levels = c("random", "var~mean(mv)", "logvar~logmean", 
                                  "poisson_scran(pos)", "dispersion(disp)", 
                                  "mean_max(max)", "SCT"))
)

avg_expr_list <- readRDS("/dcs05/hongkai/data/yli6/rzhao/mixhvg/avg_lognc.rds") # Average log-normalized counts
hvg_list <- readRDS("/dcs05/hongkai/data/yli6/rzhao/mixhvg/full_hvg_list.rds")
marker_list <- readRDS(file.path(marker_dir, "pct0_logfc0.25.rds"))

## Make marker percentage in top expression genes subplot

# Function to get average gene expression across cells
get_gene_avg_expr <- function(gene_avgs, marker_genes, sample_name) {
  # Identify markers and non-markers
  marker_genes <- intersect(marker_genes, names(gene_avgs))
  gene_type <- ifelse(names(gene_avgs) %in% marker_genes, "Marker", "NonMarker")
  
  tibble(
    Gene = names(gene_avgs),
    AvgExpression = gene_avgs,
    GeneSet = gene_type,
    Dataset = sample_name, 
    DatasetType = case_when(
      sample_name %in% cell_sorting_ds ~ "Cell Sorting", 
      sample_name %in% cite_ds ~ "CITE-seq", 
      sample_name %in% mult_ds ~ "MultiomeATAC"
    )
  )
}

# Apply to each dataset
plot_df <- imap_dfr(avg_expr_list, function(gene_avgs, sample_name) {
  markers <- marker_list[[sample_name]]
  get_gene_avg_expr(gene_avgs, markers, sample_name)
}) %>%
  group_by(Dataset) %>%
  mutate(ExpressionPercentile = percent_rank(AvgExpression)) %>%
  ungroup()

# Compute percentage of marker genes above certain k cutoff
pct_marker_cutoff <- plot_df %>% 
  group_by(Dataset, DatasetType) %>% 
  summarise(MarkerPercentage = sum(rank(1 - ExpressionPercentile) <= 2000 & GeneSet == "Marker") / sum(GeneSet == "Marker"))
g1 <- ggplot(pct_marker_cutoff, 
       aes(x = DatasetType, y = MarkerPercentage)) +
  geom_boxplot(aes(fill = DatasetType)) + 
  geom_jitter() +
  scale_fill_manual(values = c("Cell Sorting" = "#619CFF", "CITE-seq" = "#00BA38", "MultiomeATAC" = "#F8766D")) +
  labs(
    x = "Data Type", 
    y = "Percentage of Marker Genes", 
    fill = "Dataset Type"
  ) +
  ylim(c(0.05, 0.5)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold"), 
    axis.title.y = element_text(face = "bold"), 
    legend.title = element_text(face = "bold"), 
    legend.text = element_text(face = "bold"), 
    strip.text = element_text(face = "bold"), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

## Make overlap curves subplot

# Get marker gene overlap for top-k HVGs
get_overlap_pct <- function(dataset, method) {
  marker_genes <- marker_list[[dataset]]
  hvgs <- hvg_list[[dataset]][[method]]
  overlap_pct <- cumsum(hvgs %in% marker_genes) / length(marker_genes)
  if (length(overlap_pct) < 2000) {
    overlap_pct <- c(overlap_pct, 
                     rep(max(overlap_pct), 2000 - length(overlap_pct)))
  }
  overlap_pct[1:2000]
}

# Flatten into long format
results <- expand.grid(
  method = methods,
  dataset = dataset_name_range,
  stringsAsFactors = FALSE
) %>%
  left_join(dataset_info, by = "dataset") %>%
  rowwise() %>%
  mutate(k = list(1:2000), 
         overlap_pct = list(get_overlap_pct(dataset, method))) %>%
  unnest_longer(c(k, overlap_pct)) %>%
  left_join(method_info, by = "method")

# Aggregate: mean overlap across datasets per k, category, method
summary_df <- results %>%
  group_by(k, method, data_type, method_type) %>%
  summarise(mean_overlap_pct = mean(overlap_pct), .groups = 'drop')

# Plot overlap for each data type
g2 <- ggplot(summary_df, aes(x = k, y = mean_overlap_pct, color = method_type, group = method)) +
  geom_line(data = summary_df %>% subset(method_type != "mean_max(max)"), alpha = 0.8) +
  geom_line(data = summary_df %>% subset(method_type == "mean_max(max)"), alpha = 1) +
  facet_wrap(~data_type) +
  labs(x = "Top-k Selected Genes", y = "Avg. Marker Gene Overlap Percentage",
       color = "Method Type") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(), # Add lines to both axes
    axis.title.x = element_text(face = "bold"), 
    axis.title.y = element_text(face = "bold"), 
    legend.title = element_text(face = "bold"), 
    legend.text = element_text(face = "bold"), 
    strip.text = element_text(face = "bold")
  ) +
  scale_color_manual(
    values = c("#808080", "#00AFBB", "#1B9E77", "#FC074E", "#4682B4", 
               "#BEAED4", "#E7B800")
  )

combined_plot <- plot_grid(g1, g2, ncol = 2, 
                           scale = 0.9, 
                           rel_widths = c(1, 2), 
                           labels = "auto")
ggsave("mean_max_sorting_plot.pdf", plot = combined_plot, width = 300, height = 100, unit = "mm")
