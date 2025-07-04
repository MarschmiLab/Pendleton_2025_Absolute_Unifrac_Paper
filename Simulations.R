# Super simple
library(tidyverse)
library(phyloseq)
library(ggtree)

load("~/Local_Desktop/AAH_ONT_2023/Whole_Fraction_16S/data/08_compositional_exports/full_abs_physeq.RData")
melted <- full_abs_physeq %>%
  psmelt()

keep_ASVS <- melted %>%
  group_by(ASV, Class) %>%
  summarize(total_Abund = sum(Abundance)) %>%
  filter(Class %in% c("Gammaproteobacteria","Actinobacteria")) %>%
  ungroup() %>%
  group_by(Class) %>%
  slice_max(n = 2, order_by = total_Abund) %>%
  pull(ASV)


smol_p <- full_abs_physeq %>%
  subset_taxa(ASV %in% keep_ASVS)

taxa_names(smol_p) <- paste("ASV", c(1:4), sep = "_")

smol_tree <- phy_tree(smol_p)
plot(smol_tree)


abs_table <- data.frame(row.names = paste("S", c(1:5), sep = "_"),
           ASV_1 = c(10,0,0,20,8),
           ASV_2 = c(0,10,0,0,0),
           ASV_3 = c(5,5,5,10,7),
           ASV_4 = c(0,0,10,0,0))

rel_table <- abs_table %>%
  mutate(Sum = rowSums(across(everything()))) %>%
  mutate(across(everything(),\(x)x/Sum)) %>%
  select(-Sum)

abs_bray <- vegan::vegdist(abs_table, method = "bray", binary = FALSE)
rel_bray <- vegan::vegdist(rel_table, method = "bray", binary = FALSE)

library(GUniFrac)

rel_uni_raw <- GUniFrac(otu.tab = abs_table, tree = smol_tree, alpha = 1, normalize_counts = TRUE)

abs_uni_raw <- GUniFrac(otu.tab = abs_table, tree = smol_tree, alpha = 1, normalize_counts = FALSE)

rel_uni <- as.dist(rel_uni_raw$unifracs[, , "d_1"])
abs_uni <- as.dist(abs_uni_raw$unifracs[, , "d_1"])

dists <- list(Bray_Relative = rel_bray,
     Bray_Absolute = abs_bray,
     Unifrac_Relative = rel_uni,
     Unifrac_Absolute = abs_uni)


map(dists, \(x)broom::tidy(x) %>% filter(item1 == "S_1")) %>%
  bind_rows(.id = "Distance") %>%
  ggplot(aes(y = Distance, x = distance)) + 
  geom_point() +
  facet_wrap(~item2, scales= "free")


expand.grid(rep(list(c(1,10,50)),4)) %>%
  filter(rowSums(.) != 0)

possible_abunds <- expand.grid(rep(list(c(1,10,50)),4)) %>%
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

uni_rel <- as.dist(GUniFrac(otu.tab = abs_abunds, tree = smol_tree, alpha = 1, normalize_counts = TRUE)$unifracs[, , "d_1"])

uni_abs <- as.dist(GUniFrac(otu.tab = abs_abunds, tree = smol_tree, alpha = 1, normalize_counts = FALSE)$unifracs[, , "d_1"])

dists <- list(Bray_Relative = bray_rel,
              Bray_Absolute = bray_abs,
              Unifrac_Relative = uni_rel,
              Unifrac_Absolute = uni_abs)

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
  slice_max(n = 1, order_by = Difference)

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
  slice_min(n = 1, order_by = Difference)

abs_abunds[rownames(abs_abunds) %in% c("S_1","S_3", "S_7", "S_8", "S_55", "S_63", "S_80", "S_81"),]


map(dists, \(x)broom::tidy(x)) %>%
  bind_rows(.id = "Distance") %>%
  mutate(Comp_Key = paste(item1, item2, sep = ":")) %>%
  select(-item1, -item2) %>%
  pivot_wider(names_from = Distance, values_from = distance) %>%
  rowwise() %>%
  mutate(BR = Bray_Relative - Bray_Absolute,
         BA = Bray_Absolute - Unifrac_Absolute,
         UR = Unifrac_Relative - Unifrac_Absolute) %>%
  pivot_longer(BR:UR, names_to = "Comparison", values_to = "Difference") %>%
  mutate(Comparison = factor(Comparison, levels = c("BR","BA","UR"),
                             labels = c("Rel.Bray - Abs.Bray",
                                        "Abs.Bray - Abs.Uni",
                                        "Rel.Uni - Abs.Uni"))) %>%
  ggplot(aes(x = Difference)) + 
  geom_density() + 
  facet_wrap(~Comparison) + 
  geom_vline(xintercept = 0) + 
  theme_classic()

# Using real data

load("full_abs_physeq.RData")
load("absolute_wunifrac.RData")
load("relative_unifrac.RData")

# Maybe remove upwelling stations, cause that's a WHOLE other story
full_abs_physeq <- full_abs_physeq %>%
  subset_samples(!(Upwelling == "Upwelling"&month == "September"))

#relative_wunifrac <- GUniFrac(otu.tab = as.matrix(otu_table(full_abs_physeq)), tree = phy_tree(full_abs_physeq), alpha = c(0,0.5,1), normalize_counts = TRUE)
#save(relative_wunifrac, file = "relative_unifrac.RData")

#absolute_wunifrac <- GUniFrac(otu.tab = as.matrix(otu_table(full_abs_physeq)), tree = phy_tree(full_abs_physeq), alpha = c(0,0.5,1), normalize_counts = FALSE)
#save(absolute_wunifrac, file = "absolute_wunifrac.RData")



abs_unifrac_dist <- as.dist(absolute_wunifrac$unifracs[, , "d_1"])
abs_unifrac_dist_0.5 <- as.dist(absolute_wunifrac$unifracs[, , "d_0.5"])
abs_unifrac_dist_0 <- as.dist(absolute_wunifrac$unifracs[, , "d_0"])

rel_unifrac_dist <- as.dist(relative_wunifrac$unifracs[, , "d_1"])
rel_unifrac_dist_0.5 <- as.dist(relative_wunifrac$unifracs[, , "d_0.5"])
rel_unifrac_dist_0 <- as.dist(relative_wunifrac$unifracs[, , "d_0"])

ordinate_and_plot <- function(dist, physeq, title){
  plot_ordination(physeq = physeq,
                  ordination = ordinate(physeq, method = "PCoA", distance = dist),
                  color = "Comp_Group_Hier_Colors") + 
    labs(title = title,
         color = "Group") + 
    scale_color_manual(values = comp_group_colors_hier)+ 
    geom_point(alpha = 0.5, size = 4)
}

ont_dists <- list(Rel_0 = rel_unifrac_dist_0,
                  Rel_5 = rel_unifrac_dist_0.5,
                  Rel_1 = rel_unifrac_dist,
                  Abs_0 = abs_unifrac_dist_0,
                  Abs_5 = abs_unifrac_dist_0.5,
                  Abs_1 = abs_unifrac_dist)

plots <- map2(ont_dists, c("WUnifac, Relative Abundance, 0.0",
                           "WUnifac, Relative Abundance, 0.5",
                           "WUnifac, Relative Abundance, 1.0",
                           "WUnifac, Absolute Abundance, 0.0",
                           "WUnifac, Absolute Abundance, 0.5",
                           "WUnifac, Absolute Abundance, 1.0"),
              \(x,y)ordinate_and_plot(dist = x, physeq = full_abs_physeq, title = y))

library(patchwork)

wrap_plots(plots) + 
  plot_layout(ncol = 3, guides = "collect")

cell_dist <- sample_data(full_abs_physeq) %>%
  data.frame() %>%
  select(avg_cells_per_ml) %>% 
  dist()

map(ont_dists, \(x)vegan::mantel(cell_dist, x, permutations = 999))

# Lake Erie

load("erie_abs_physeq.RData")
load("erie_absolute_wunifrac.RData")
load("erie_relative_unifrac.RData")

# This one sample is a confusing outlier, that I haven't really decided on yet
erie_abs_physeq <- erie_abs_physeq %>%
  subset_samples(DNA_ID != "KH_D076")

# erie_relative_wunifrac <- GUniFrac(otu.tab = t(as.matrix(otu_table(erie_abs_physeq))), tree = phy_tree(erie_abs_physeq), alpha = c(0,0.5,1), normalize_counts = TRUE)
# save(erie_relative_wunifrac, file = "relative_unifrac.RData")
# 
# erie_absolute_wunifrac <- GUniFrac(otu.tab = t(as.matrix(otu_table(erie_abs_physeq))), tree = phy_tree(erie_abs_physeq), alpha = c(0,0.5,1), normalize_counts = FALSE)
# save(erie_absolute_wunifrac, file = "erie_absolute_wunifrac.RData")



abs_unifrac_dist <- as.dist(erie_absolute_wunifrac$unifracs[, , "d_1"])
abs_unifrac_dist_0.5 <- as.dist(erie_absolute_wunifrac$unifracs[, , "d_0.5"])
abs_unifrac_dist_0 <- as.dist(erie_absolute_wunifrac$unifracs[, , "d_0"])

rel_unifrac_dist <- as.dist(erie_relative_wunifrac$unifracs[, , "d_1"])
rel_unifrac_dist_0.5 <- as.dist(erie_relative_wunifrac$unifracs[, , "d_0.5"])
rel_unifrac_dist_0 <- as.dist(erie_relative_wunifrac$unifracs[, , "d_0"])

ordinate_and_plot <- function(dist, physeq, title){
  plot_ordination(physeq = physeq,
                  ordination = ordinate(physeq, method = "PCoA", distance = dist),
                  color = "Month") + 
    labs(title = title,
         color = "Group") + 
    scale_color_manual(values = erie_month_colors)+ 
    geom_point(alpha = 0.5, size = 4)
}

ont_dists <- list(Rel_0 = rel_unifrac_dist_0,
                  Rel_5 = rel_unifrac_dist_0.5,
                  Rel_1 = rel_unifrac_dist,
                  Abs_0 = abs_unifrac_dist_0,
                  Abs_5 = abs_unifrac_dist_0.5,
                  Abs_1 = abs_unifrac_dist)

plots <- map2(ont_dists, c("WUnifac, Relative Abundance, 0.0",
                           "WUnifac, Relative Abundance, 0.5",
                           "WUnifac, Relative Abundance, 1.0",
                           "WUnifac, Absolute Abundance, 0.0",
                           "WUnifac, Absolute Abundance, 0.5",
                           "WUnifac, Absolute Abundance, 1.0"),
              \(x,y)ordinate_and_plot(dist = x, physeq = erie_abs_physeq, title = y))

wrap_plots(plots) + 
  plot_layout(ncol = 3, guides = "collect")

cell_dist <- sample_data(erie_abs_physeq) %>%
  data.frame() %>%
  select(avg_cells_per_ml) %>% 
  dist()

map(ont_dists, \(x)vegan::mantel(cell_dist, x, permutations = 999))
