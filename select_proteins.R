library(dplyr)

protein_type = "CsgA"

all_lines <- readLines("example.fasta")

prot_names <- all_lines[grepl("^>", all_lines)] 

all_prots <- split(all_lines, cumsum(grepl("^>", all_lines)))
all_prots[sapply(all_prots, function(i)
  grepl(protein_type, i[1], ignore.case = TRUE))] %>% 
  unlist %>% 
  writeLines(con = "selected.fasta")
