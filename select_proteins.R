library(dplyr)

protein_type <- "CsgA"
input_file <- "example.fasta"

all_lines <- readLines(input_file)

all_prots <- split(all_lines, cumsum(grepl("^>", all_lines)))
all_prots[sapply(all_prots, function(i)
  grepl(protein_type, i[1], ignore.case = TRUE))] %>% 
  unlist %>% 
  writeLines(con = "selected.fasta")
