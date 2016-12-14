library(biomaRt)
# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Test
# query biomart
results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id"),
                 filters = "ensembl_peptide_id", values = "ENSP00000053867",
                 mart = mart)
results


diseases <- read.csv('D:\\Projects\\MultiBioNet\\data\\alz_DISEASES', sep = '\t', header = T, stringsAsFactors = F)
ensp <- c(diseases[1])

results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id"),
                 filters = "ensembl_peptide_id", values = ensp,
                 mart = mart)

results[1]

write.csv(results[1], 'D:\\Projects\\MultiBioNet\\data\\alz_DISEASES_ENSG')
