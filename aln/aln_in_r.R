library(Biostrings)

seqs <- readAAStringSet("./aln/csgA_fapC.fasta")

aln <- pairwiseAlignment(seqs[[1]], seqs[[2]],
                         substitutionMatrix = "BLOSUM62",
                         gapOpening = 10, gapExtension = 0.5, type = "local")

aln_glob <- pairwiseAlignment(seqs[[1]], seqs[[2]],
                         substitutionMatrix = "BLOSUM62",
                         gapOpening = 10, gapExtension = 0.5, type = "global")

cat(">CsgA",
    as.character(aln@pattern),
    ">FapC",
    as.character(aln@subject), file = "./aln/csgA_fapC_local.aln", sep = "\n")

cat(">CsgA",
    as.character(aln_glob@pattern),
    ">FapC",
    as.character(aln_glob@subject), file = "./aln/csgA_fapC_global.aln", sep = "\n")

pid(aln_glob, type = "PID3")
