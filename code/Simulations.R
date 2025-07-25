# Super simple
library(tidyverse)
library(phyloseq)
library(ggtree)
library(patchwork)
library(vegan)
library(GUniFrac)

source("code/plotting.R")

nonsym_tree <- ape::read.tree("data/nonsym_tree.tree")
sym_tree <- ape::read.tree("data/sym_tree.tree")

sym_tree_plot <- ggtree(sym_tree) + 
  geom_tiplab(label = paste("ASV_", c(1:4))) + 
  geom_treescale(x = 0, width = 0.25) + 
  scale_x_continuous(expand = expansion(c(0,.5))) + 
  theme(plot.margin = unit(c(0,0,0, 0), "mm"))

smol_tree_plot <- ggtree(nonsym_tree) + 
  geom_tiplab(label = paste("ASV_", c(1:4))) + 
  geom_treescale(x = 0, width = 0.25) + 
  scale_x_continuous(expand = expansion(c(0,.5))) + 
  theme(plot.margin = unit(c(0,0,0, 0), "mm"))

flipped_tree_plot <- ggtree(nonsym_tree) + 
  geom_tiplab(label = paste("ASV_", c(1:4)), hjust = 1) + 
  geom_treescale(x = 0, width = 0.25) + 
  scale_x_reverse(expand = expansion(c(.5,0))) + 
  theme(plot.margin = unit(c(0,0,0, 0), "mm"))


bf <- ggtree(nonsym_tree) + 
  theme(plot.margin = unit(c(0,5,0, 0), "mm"),
        panel.background = element_rect(fill = "transparent",
                                        color = "transparent"),
        plot.background = element_rect(fill = "transparent",
                                       color = "transparent"))
br <- ggtree(nonsym_tree) + 
  scale_x_reverse() +
  theme(plot.margin = unit(c(0,0,0, 5), "mm"),
        panel.background = element_rect(fill = "transparent",
                                        color = "transparent"),
        plot.background = element_rect(fill = "transparent",
                                       color = "transparent"))

facing <- bf + br 
ggsave(facing, filename = "figures/facing.pdf", width = 2, height = 1)
ggsave(smol_tree_plot, filename = "figures/smol_tree_plot.pdf", width = 2, height = 2)
ggsave(sym_tree_plot, filename = "figures/sym_tree_plot.pdf", width = 2, height = 2)
ggsave(flipped_tree_plot, filename = "figures/flipped_tree_plot.png", width = 2, height = 2)

possible_abunds <- expand.grid(rep(list(c(1,10,100)),4)) %>%
  filter(rowSums(.) != 0)


abs_abunds <- possible_abunds %>%
  as.data.frame()

colnames(abs_abunds) <- paste("ASV", c(1:4), sep = "_")
rownames(abs_abunds) <-  paste("S", c(1:nrow(abs_abunds)), sep = "_")

rel_abunds <- abs_abunds %>%
  mutate(Sum = rowSums(across(everything()))) %>%
  mutate(across(everything(),\(x)x/Sum)) %>%
  select(-Sum)

bray_rel <- vegan::vegdist(rel_abunds, method = "bray", binary = FALSE)
bray_abs <- vegan::vegdist(abs_abunds, method = "bray", binary = FALSE)

uni_rel <- as.dist(GUniFrac(otu.tab = abs_abunds, tree = nonsym_tree, alpha = 1, normalize_counts = TRUE)$unifracs[, , "d_1"])

uni_abs <- as.dist(GUniFrac(otu.tab = abs_abunds, tree = nonsym_tree, alpha = 1, normalize_counts = FALSE)$unifracs[, , "d_1"])

dists <- list(Bray_Relative = bray_rel,
              Bray_Absolute = bray_abs,
              Unifrac_Relative = uni_rel,
              Unifrac_Absolute = uni_abs)

map(dists, \(x)broom::tidy(x)) %>%
  bind_rows(.id = "Distance") %>%
  mutate(Comp_Key = paste(item1, item2, sep = ":")) %>%
  select(-item1, -item2) %>%
  pivot_wider(names_from = Distance, values_from = distance) %>%
  select(-Comp_Key) %>%
  corrr::correlate() %>%
  corrr::focus(Unifrac_Absolute)

map(dists, \(x)broom::tidy(x)) %>%
  bind_rows(.id = "Distance") %>%
  mutate(Comp_Key = paste(item1, item2, sep = ":")) %>%
  select(-item1, -item2) %>%
  pivot_wider(names_from = Distance, values_from = distance) %>%
  summarize(cors = list(cor.test(Bray_Absolute, Unifrac_Absolute))) %>%
  pull(cors)

map(dists, \(x)broom::tidy(x)) %>%
  bind_rows(.id = "Distance") %>%
  mutate(Comp_Key = paste(item1, item2, sep = ":")) %>%
  select(-item1, -item2) %>%
  pivot_wider(names_from = Distance, values_from = distance) %>%
  rowwise() %>%
  mutate(BR = Bray_Relative - Unifrac_Absolute,
         BA = Bray_Absolute - Unifrac_Absolute,
         UR = Unifrac_Relative - Unifrac_Absolute) %>%
  pivot_longer(BR:UR, names_to = "Comparison", values_to = "Difference") %>%
  group_by(Comparison) %>%
  slice_max(n = 1, order_by = Difference, with_ties = TRUE)

map(dists, \(x)broom::tidy(x)) %>%
  bind_rows(.id = "Distance") %>%
  mutate(Comp_Key = paste(item1, item2, sep = ":")) %>%
  select(-item1, -item2) %>%
  pivot_wider(names_from = Distance, values_from = distance) %>%
  rowwise() %>%
  mutate(BR = Bray_Relative - Unifrac_Absolute,
         BA = Bray_Absolute - Unifrac_Absolute,
         UR = Unifrac_Relative - Unifrac_Absolute) %>%
  pivot_longer(BR:UR, names_to = "Comparison", values_to = "Difference") %>%
  group_by(Comparison) %>%
  slice_min(n = 1, order_by = Difference, with_ties = TRUE)


abs_abunds[rownames(abs_abunds) %in% c("S_1","S_3", "S_6", "S_7", "S_8", "S_46", "S_54","S_64","S_72", "S_78","S_80", "S_81"),]


map(dists, \(x)broom::tidy(x)) %>%
  bind_rows(.id = "Distance") %>%
  mutate(Comp_Key = paste(item1, item2, sep = ":")) %>%
  select(-item1, -item2) %>%
  pivot_wider(names_from = Distance, values_from = distance) %>%
  rowwise() %>%
  mutate(BR = Bray_Relative - Unifrac_Absolute,
         BA = Bray_Absolute - Unifrac_Absolute,
         UR = Unifrac_Relative - Unifrac_Absolute) %>%
  pivot_longer(BR:UR, names_to = "Comparison", values_to = "Difference") %>%
  mutate(Comparison = factor(Comparison, levels = c("BR","UR", "BA"),
                             labels = c(expression(italic(BC^{R}-U^{A})),
                                        expression(italic(U^{R}-U^{A})),
                                        expression(italic(BC^{A}-U^{A}))))) %>%
  ggplot(aes(x = Difference)) + 
  geom_vline(xintercept = 0, linetype = 2, linewidth = .2) + 
  geom_histogram(bins = 50, alpha = 0.7, fill = "#6A3D9A") + 
  facet_wrap(~Comparison, labeller = label_parsed) + 
  theme_classic() + 
  scale_y_continuous(expand = expansion(mult = 0)) + 
  scale_x_continuous(limits = c(-1,1)) + 
  labs(y = "Count", x = "Difference in Dissimilarity") + 
  theme(strip.background = element_rect(color = "white"),
        strip.text.x = element_text(size = 10))

ggsave(filename = "figures/diff_hists.pdf", width = 6, height = 2)

# Supplementary stuff:
uni_abs_sym <- as.dist(GUniFrac(otu.tab = abs_abunds, tree = sym_tree, alpha = 1, normalize_counts = FALSE)$unifracs[, , "d_1"])

dists <- list(Bray_Absolute = bray_abs,
              Unifrac_Absolute = uni_abs_sym)

map(dists, \(x)broom::tidy(x)) %>%
  bind_rows(.id = "Distance") %>%
  mutate(Comp_Key = paste(item1, item2, sep = ":")) %>%
  select(-item1, -item2) %>%
  pivot_wider(names_from = Distance, values_from = distance) %>%
  rowwise() %>%
  mutate( BA = Bray_Absolute - Unifrac_Absolute) %>%
  pivot_longer(BA, names_to = "Comparison", values_to = "Difference") %>%
  mutate(Comparison = factor(Comparison, levels = c("BA"),
                             labels = c(expression(italic(BC^{A}-U^{A}))))) %>%
  ggplot(aes(x = Difference)) + 
  geom_vline(xintercept = 0, linetype = 2, linewidth = .2) + 
  geom_histogram(bins = 50, alpha = 0.7, fill = "#6A3D9A") + 
  facet_wrap(~Comparison, labeller = label_parsed) + 
  theme_classic() + 
  scale_y_continuous(expand = expansion(mult = 0)) + 
  scale_x_continuous(limits = c(-1,1)) + 
  labs(y = "Count", x = "Difference in Dissimilarity") + 
  theme(strip.background = element_rect(color = "white"),
        strip.text.x = element_text(size = 10))

ggsave(filename = "figures/sym_hists.pdf", width = 3, height = 3)

# Using real data

load("data/full_abs_physeq.RData")
load("data/ont_absolute_wunifrac.RData")
load("data/ont_relative_unifrac.RData")

# Maybe remove upwelling stations, cause that's a WHOLE other story
full_abs_physeq <- full_abs_physeq %>%
  subset_samples(!(Upwelling == "Upwelling"&month == "September"))

# ont_relative_wunifrac <- GUniFrac(otu.tab = as.matrix(otu_table(full_abs_physeq)), tree = phy_tree(full_abs_physeq), alpha = c(0,0.5,1), normalize_counts = TRUE)
# save(ont_relative_wunifrac, file = "ont_relative_unifrac.RData")
# 
# ont_absolute_wunifrac <- GUniFrac(otu.tab = as.matrix(otu_table(full_abs_physeq)), tree = phy_tree(full_abs_physeq), alpha = c(0,0.5,1), normalize_counts = FALSE)
# save(ont_absolute_wunifrac, file = "ont_absolute_wunifrac.RData")



abs_unifrac_dist <- as.dist(ont_absolute_wunifrac$unifracs[, , "d_1"])
abs_unifrac_dist_0.5 <- as.dist(ont_absolute_wunifrac$unifracs[, , "d_0.5"])
abs_unifrac_dist_0 <- as.dist(ont_absolute_wunifrac$unifracs[, , "d_0"])

rel_unifrac_dist <- as.dist(ont_relative_wunifrac$unifracs[, , "d_1"])
rel_unifrac_dist_0.5 <- as.dist(ont_relative_wunifrac$unifracs[, , "d_0.5"])
rel_unifrac_dist_0 <- as.dist(ont_relative_wunifrac$unifracs[, , "d_0"])

ordinate_and_plot <- function(dist, physeq, title){
  plot_ordination(physeq = physeq,
                  ordination = ordinate(physeq, method = "PCoA", distance = dist),
                  color = "Comp_Group_Hier") + 
    labs(title = title,
         color = "Group") + 
    scale_color_manual(values = comp_group_colors_three)+ 
    geom_point(alpha = 0.5, size = 4) + 
    theme(axis.title = element_text(size = 20),
          axis.title.x = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 25))
}

ont_dists <- list(Rel_0 = rel_unifrac_dist_0,
                  Rel_5 = rel_unifrac_dist_0.5,
                  Rel_1 = rel_unifrac_dist,
                  Abs_0 = abs_unifrac_dist_0,
                  Abs_5 = abs_unifrac_dist_0.5,
                  Abs_1 = abs_unifrac_dist)

plots <- map2(ont_dists, c("\u03b1 = 0.0",
                           "\u03b1 = 0.5",
                           "\u03b1 = 1.0",
                           "\u03b1 = 0.0",
                           "\u03b1 = 0.5",
                           "\u03b1 = 1.0"),
              \(x,y)ordinate_and_plot(dist = x, physeq = full_abs_physeq, title = y))



wrap_plots(list(plots[[4]] + scale_x_reverse(), plots[[5]], plots[[6]])) + 
  plot_layout(ncol = 3, guides = "collect")

ggsave(filename = "figures/ont_ordinations.png",width = 14, height = 5, dpi = 600)


wrap_plots(plots[1:3]) + 
  plot_layout(ncol = 3, guides = "collect")

ggsave(filename = "figures/ont_ordinations_rel.png",width = 14, height = 5, dpi = 600)

cell_dist <- sample_data(full_abs_physeq) %>%
  data.frame() %>%
  select(avg_cells_per_ml) %>% 
  dist()

ont_dists_bray <- list(Bray_Rel = distance(transform_sample_counts(full_abs_physeq, \(x)x/sum(x)), method = "bray"),
                       Bray_Abs = distance(full_abs_physeq, method = "bray"),
                       Rel_0 = rel_unifrac_dist_0,
                       Rel_0.5 = rel_unifrac_dist_0.5,
                       Rel_1 = rel_unifrac_dist,
                       Abs_0 = abs_unifrac_dist_0,
                       Abs_0.5 = abs_unifrac_dist_0.5,
                       Abs_1 = abs_unifrac_dist)

ont_mantels <- map(ont_dists_bray, \(x)vegan::mantel(cell_dist, x, permutations = 999)$statistic)

ont_mantels %>%
  bind_rows() %>%
  pivot_longer(everything(), names_to = "Metric", values_to = "Correlation") %>%
  filter(!str_detect(Metric, "Bray")) %>%
  separate_wider_delim(Metric, delim = "_", names = c("Normalization", "Alpha")) %>%
  mutate(Alpha = as.numeric(Alpha)) %>%
  ggplot(aes(x = Alpha, y = Correlation, color = Normalization)) + 
  annotate(geom = "segment", y = ont_mantels$Bray_Rel, yend = ont_mantels$Bray_Rel,
           x = 0, xend = 1, linetype = 2, color = "#196689", linewidth = 2) + 
  annotate(geom = "segment", y = ont_mantels$Bray_Abs, yend = ont_mantels$Bray_Abs,
           x = 0, xend = 1, linetype = 2, color = "#B6AA0D", linewidth = 2) +  
  labs(x = "\u03b1", y = expression(Mantel~italic(r)~"(dist. vs. cells)")) + 
  geom_line(linewidth = 2) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+ 
  scale_color_manual(values = c("#F2C621", "#008FF8"))+ 
  theme(legend.position = "none",
        axis.title = element_text(size = 20))


ggsave(filename = "figures/ont_mantel.png", width = 4, height = 3.5, dpi=600)

# Permanovas
sam_data_for_adonis <- data.frame(sample_data(full_abs_physeq))


ont_perms <- map(ont_dists_bray, \(x){
  adonis2(x ~ Comp_Group_Hier, data = sam_data_for_adonis, by = "margin") %>%
    broom::tidy()
})

ont_perms %>%
  bind_rows(.id = "Metric") %>%
  filter(!str_detect(Metric, "Bray")) %>%
  separate_wider_delim(Metric, delim = "_", names = c("Normalization", "Alpha")) %>%
  mutate(Alpha = as.numeric(Alpha)) %>%
  filter(term == "Comp_Group_Hier") %>%
  ggplot(aes(x = Alpha, y = R2, color = Normalization)) + 
  annotate(geom = "segment", y = ont_perms$Bray_Rel[[1,4]], yend = ont_perms$Bray_Rel[[1,4]],
           x = 0, xend = 1, linetype = 2, color = "#196689", linewidth = 2) + 
  annotate(geom = "segment", y = ont_perms$Bray_Abs[[1,4]], yend = ont_perms$Bray_Abs[[1,4]],
           x = 0, xend =1, linetype = 2, color = "#B6AA0D", linewidth = 2) +
  geom_line(linewidth = 2) + 
  labs(x = "\u03b1", y = expression(italic(R^2))) +
  scale_x_continuous(expand = expansion(mult = c(0,0.05)))+ 
  scale_color_manual(values = c("#F2C621", "#008FF8"))+ 
  theme(legend.position = "none",
        axis.title = element_text(size = 20))

ggsave(filename = "figures/ont_rsqr.png", width = 4, height = 3.5, dpi=600)

ont_perms %>%
  bind_rows(.id = "Metric") %>%
  filter(!str_detect(Metric, "Bray")) %>%
  separate_wider_delim(Metric, delim = "_", names = c("Normalization", "Alpha")) %>%
  mutate(Alpha = as.numeric(Alpha)) %>%
  filter(term == "Comp_Group_Hier") %>%
  ggplot(aes(x = Alpha, y = statistic, color = Normalization)) + 
  annotate(geom = "segment", y = ont_perms$Bray_Rel[[1,5]], yend = ont_perms$Bray_Rel[[1,5]],
           x = 0, xend = 1, linetype = 2, color = "#196689", linewidth = 2) + 
  annotate(geom = "segment", y = ont_perms$Bray_Abs[[1,5]], yend = ont_perms$Bray_Abs[[1,5]],
           x = 0, xend =1, linetype = 2, color = "#B6AA0D", linewidth = 2) +
  geom_line(linewidth = 2) + 
  labs(x = "\u03b1", y = expression(italic(F)-statistic)) +
  scale_x_continuous(expand = expansion(mult = c(0,0.05)))+ 
  scale_color_manual(values = c("#F2C621", "#008FF8"))+ 
  theme(legend.position = "none",
        axis.title = element_text(size = 20))

ggsave(filename = "figures/ont_fstat.png", width = 4, height = 3.5, dpi=600)

## Correlation of Axis1 and Cell counts

pcoa_res <- ordinate(full_abs_physeq, method = "PCoA", distance = abs_unifrac_dist)$vectors %>%as.data.frame()

cell_df <- sample_data(full_abs_physeq) %>%
  data.frame() %>%
  select(avg_cells_per_ml) %>%
  mutate(sample = row.names(.))

pcoa_res %>%
  select(Axis.1) %>%
  mutate(sample = row.names(.)) %>%
  left_join(cell_df) %>%
  summarize(cor = list(cor.test(Axis.1, avg_cells_per_ml, method = "spearman"))) %>%
  pull(cor)

pcoa_res <- ordinate(full_abs_physeq, method = "PCoA", distance = abs_unifrac_dist_0.5)$vectors %>%as.data.frame()

cell_df <- sample_data(full_abs_physeq) %>%
  data.frame() %>%
  select(avg_cells_per_ml) %>%
  mutate(sample = row.names(.))

pcoa_res %>%
  select(Axis.1) %>%
  mutate(sample = row.names(.)) %>%
  left_join(cell_df) %>%
  summarize(cor = list(cor.test(Axis.1, avg_cells_per_ml, method = "spearman"))) %>%
  pull(cor)

pcoa_res <- ordinate(full_abs_physeq, method = "PCoA", distance = abs_unifrac_dist_0)$vectors %>%as.data.frame()

cell_df <- sample_data(full_abs_physeq) %>%
  data.frame() %>%
  select(avg_cells_per_ml) %>%
  mutate(sample = row.names(.))

pcoa_res %>%
  select(Axis.1) %>%
  mutate(sample = row.names(.)) %>%
  left_join(cell_df) %>%
  summarize(cor = list(cor.test(Axis.1, avg_cells_per_ml, method = "spearman"))) %>%
  pull(cor)
