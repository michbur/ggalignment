library(dplyr)
library(ggplot2)
library(reshape2)

compute_pid <- function(file) {
  all_lines <- readLines(file)
  
  all_prots <- split(all_lines, cumsum(grepl("^>", all_lines)))
  
  all_prots <- all_prots[order(sapply(all_prots, first))]
  
  splitted_seqs <- sapply(all_prots, function(i) paste0(i[-1], collapse = "")) %>% 
    strsplit("")
  
  lapply(combn(1L:length(all_prots), 2, simplify = FALSE), function(i) {
    pid <- data.frame(S1 = splitted_seqs[[i[1]]], S2 = splitted_seqs[[i[2]]], stringsAsFactors = FALSE) %>% 
      filter(S1 != "-", S2 != "-") %>% 
      mutate(ident = S1 == S2) %>% 
      summarise(ident = mean(ident)) %>% 
      unlist %>% 
      unname()
    
    data.frame(S1 = all_prots[[i[1]]][[1]],
               S2 = all_prots[[i[2]]][[1]],
               pid = pid)
  }) %>% 
    do.call(rbind, .)
}

pid_data <- rbind(compute_pid("csga.syn.aln") %>% 
        mutate(protein = "CsgA"),
      compute_pid("CsgC_sekwencje.fasta.aln") %>% 
        mutate(protein = "CsgC")) %>% 
  mutate(S1 = sub(">", "", S1),
         S2 = sub(">", "", S2))

ggplot(pid_data, aes(x = S1, y = S2, fill = pid, label = round(pid, 2))) +
  geom_tile(color = "black") +
  geom_text(color = "red") +
  facet_wrap(~ protein) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


pid_data %>% 
  filter(S1 != "Alpha-synuclein") %>%
  dcast(S1 + S2 ~ protein, value.var = "pid") %>%  
  mutate(pid = CsgA - CsgC) %>% 
  ggplot(aes(x = S1, y = S2, fill = pid, label = round(pid, 2))) +
  geom_tile(color = "black") +
  geom_text(color = "red") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
