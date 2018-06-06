library(ggtree)

tr <- read.tree("./change_trees/MALEDRZEWKO_nice_names.aln.treefile")

cairo_pdf("tree.pdf", height = 8.6, width = 6.5)
ggtree(tr) + 
  geom_treescale() + 
  geom_tiplab(size = 2) +
  scale_x_continuous(limits = c(0, 4))
dev.off()

