library(ggplot2)
library(dplyr)

all_lines <- readLines("nucleotides.fasta")

seq_id <- cumsum(grepl("^>", all_lines))

all_seqs <- split(all_lines, seq_id)

dat <- lapply(all_seqs, function(ith_seq) {
  data.frame(name = sub("^>", "", ith_seq[1]),
             unigram = strsplit(paste0(ith_seq[-1], collapse = ""), "")[[1]],
             stringsAsFactors = FALSE)
}) %>% bind_rows() %>% 
  mutate(pos = 1L:length(name),
         pos_disc = cut(pos, quantile(pos, probs = c(0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.84, 1)), include.lowest = TRUE),
         special = pos %in% 5L:9 | pos %in% 150L:170,
         special = ifelse(special, "HTH domain", "Other"))

cairo_pdf("nuc.pdf", height = 2.8, width = 6.5)
ggplot(dat, aes(x = pos, y = name, label = unigram, fill = special)) +
  geom_tile(color = "black") +
  facet_wrap(~ pos_disc, ncol = 1, scales = "free_x") +
  geom_text(size = 2.1) +
  scale_x_continuous("Position", expand = c(0, 0)) +
  scale_y_discrete("") +
  scale_fill_manual("", values = c("lightblue", "white")) +
  theme_bw(base_size = 7) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.7, "lines"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
dev.off()
