library(dplyr)

all_lines <- readLines("./change_trees/MALEDRZEWKO.aln")

prot_id <- cumsum(grepl("^>", all_lines))

all_prots <- split(all_lines, prot_id)

ids <- sapply(all_prots, first) %>% 
  strsplit(" ") %>% 
  sapply(first) %>% 
  sub(">", "", x = .) %>% 
  sub("/", "_", x = .)

ids_only <- sapply(all_prots, first) %>% 
  strsplit("[/ ]") %>% 
  sapply(first) %>% 
  sub(">", "", x = .)

org_names <- sapply(all_prots, first) %>% 
  strsplit("[", fixed = TRUE) %>% 
  sapply(last) %>% 
  sub("]", "", x = ., fixed = TRUE)

name_df <- data.frame(id = ids, 
                      full_name = paste(ids_only, org_names), 
                      stringsAsFactors = FALSE) %>% 
  mutate(full_name = gsub(" ", "_", full_name))

tree_lines <- readLines("./change_trees/MALEDRZEWKO.aln.treefile")

for(i in 1L:nrow(name_df)) {
  tree_lines <- sub(name_df[i, "id"], name_df[i, "full_name"], tree_lines)
}

cat(tree_lines, file = "./change_trees/MALEDRZEWKO_nice_names.aln.treefile")
