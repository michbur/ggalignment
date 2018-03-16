library(ggplot2)
library(dplyr)
library(reshape2)
library(xtable)

alph_color_list <- list(list(col = "#8dd3c7", aa = c("G", "A", "S", "T")),
                     list(col = "#ffffb3", aa = c("C", "V", "I", "L", "P", "F", "Y", "M", "W")),
                     list(col = "#bebada", aa = c("N", "H", "Q")),
                     list(col = "#fb8072", aa = c("D", "E")),
                     list(col = "#80b1d3", aa = c("K", "R")),
                     list(col = "white", aa = c("-")))

reduced_alph <- lapply(alph_color_list, function(i) i[["aa"]])

names(reduced_alph) <- c("Small nonpolar",
                         "Hydrophobic",
                         "Polar",
                         "Negatively charged",
                         "Positively charged",
                         "Gap")

cols <- sapply(alph_color_list, function(i) i[["col"]])


all_lines <- readLines("csgA.aln")

prot_id <- cumsum(grepl("^>", all_lines))

all_prots <- split(all_lines, prot_id)

aln_dat <- lapply(all_prots, function(ith_prot) {
  data.frame(species = sub("^>", "", ith_prot[1]),
    aa = strsplit(paste0(ith_prot[-1], collapse = ""), "")[[1]],
             stringsAsFactors = FALSE)
}) %>% bind_rows() %>% 
  group_by(species) %>% 
  mutate(pos = 1L:length(species),
         pos_disc = cut(pos, quantile(pos, probs = c(0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.84, 1)), include.lowest = TRUE),
         aa_groups = factor(biogram::degenerate(aa, reduced_alph),
                            levels = names(reduced_alph)),
         prone = FALSE)

cairo_ps("aln.eps", height = 4.8, width = 6.5)
ggplot(aln_dat, aes(x = pos, y = species, label = aa, fill = aa_groups)) +
  geom_tile(color = NA) +
  geom_text(size = 2.1, family = "Courier", fontface="bold") +
  facet_wrap(~ pos_disc, ncol = 1, scales = "free_x") +
  scale_x_continuous("Position", expand = c(0, 0)) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(color = "black"))) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.7, "lines"),
        strip.background = element_blank(),
        strip.text.x = element_blank())
dev.off()
