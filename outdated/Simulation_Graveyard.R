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

wrap_plots(plots[4:6]) + 
  plot_layout(ncol = 3, guides = "collect")

ggsave(filename = "figures/erie_ordinations.png",width = 14, height = 5)

cell_dist <- sample_data(erie_abs_physeq) %>%
  data.frame() %>%
  select(avg_cells_per_ml) %>% 
  dist()

erie_dists_bray <- list(Bray_Rel = distance(transform_sample_counts(erie_abs_physeq, \(x)x/sum(x)), method = "bray"),
                        Bray_Abs = distance(erie_abs_physeq, method = "bray"),
                        Rel_0 = rel_unifrac_dist_0,
                        Rel_0.5 = rel_unifrac_dist_0.5,
                        Rel_1 = rel_unifrac_dist,
                        Abs_0 = abs_unifrac_dist_0,
                        Abs_0.5 = abs_unifrac_dist_0.5,
                        Abs_1 = abs_unifrac_dist)

erie_mantels <- map(erie_dists_bray, \(x)vegan::mantel(cell_dist, x, permutations = 999)$statistic)

erie_mantels %>%
  bind_rows() %>%
  pivot_longer(everything(), names_to = "Metric", values_to = "Correlation") %>%
  filter(!str_detect(Metric, "Bray")) %>%
  separate_wider_delim(Metric, delim = "_", names = c("Normalization", "Alpha")) %>%
  mutate(Alpha = as.numeric(Alpha)) %>%
  ggplot(aes(x = Alpha, y = Correlation, color = Normalization)) + 
  geom_line() + 
  annotate(geom = "segment", y = erie_mantels$Bray_Rel, yend = erie_mantels$Bray_Rel,
           x = 0, xend = 0.7, linetype = 2) + 
  annotate(geom = "segment", y = erie_mantels$Bray_Abs, yend = erie_mantels$Bray_Abs,
           x = 0, xend = 0.7, linetype = 2) +  
  annotate(geom = "text", x = 0.73, y = c(erie_mantels$Bray_Rel, erie_mantels$Bray_Abs),
           label = c("Rel. Bray", "Abs. Bray"), hjust = 0) + 
  labs(x = "Generalized alpha value", y = "Correlation with cell counts") + 
  scale_x_continuous(expand = expansion(mult = 0)) + 
  scale_color_manual(labels = c("Absolute Unifrac", "Relative Unifrac"),
                     name = "",
                     values = c("#B6AA0D", "#008FF8")) + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.17,0.85),
        legend.background = element_blank(),
        legend.box.background = element_blank())

ggsave(filename = "figures/erie_mantel.png", width = 5, height = 3.5)
