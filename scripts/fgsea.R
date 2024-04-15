setwd("~/Studijos/Transcriptomics/")

library(fgsea)
library(data.table)
library(reactome.db)
library(org.Hs.eg.db)

arg <- commandArgs(trailingOnly = TRUE)[1]

cts <- read.csv(paste("deseq/", arg, ".csv", sep = ""), header = TRUE, row.names = 1)
genes <- rownames(cts)
genes <- sub("\\..*$", "", genes)
rownames(cts) <- genes
genes <- mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "ENSEMBL")
valid_genes <- genes[!is.na(genes)]
cts <- cts[rownames(cts) %in% names(valid_genes), ]
row.names(cts) <- valid_genes
pathways <- reactomePathways(genes = valid_genes)
ranks <- setNames(cts$log2FoldChange, rownames(cts))

fgseaRes <- fgsea(pathways, ranks)
fgseaRes <- fgseaRes[order(fgseaRes$pval), ]

plotEnrichment(pathways[[fgseaRes$pathway[1]]], ranks)

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)

fgseaRes <- fgseaRes[, leadingEdge := mapIdsList(x=org.Hs.eg.db, keys=leadingEdge, keytype="ENTREZID", column="SYMBOL")]
fwrite(fgseaRes, file = paste("fgsea/", arg,"_pathways.tsv", sep = ""), sep = "\t", sep2 = c("", " ",""))

# fgseaRes <- head(fgseaRes)
# 
# round_numeric <- function(x, digits) {
#   if(is.numeric(x)) {
#     return(round(x, digits = digits))
#   } else {
#     return(x)
#   }
# }

#df_rounded <- lapply(fgseaRes, FUN = round_numeric, digits = 5)
#genes <- df_rounded[8]
#df_rounded <- data.frame(df_rounded[1:6])
#genes <- lapply(genes, function(inner) {lapply(inner, unname)})
#genes <- lapply(genes$leadingEdge, paste, collapse = ", ")
#df_rounded <- merge(df_rounded, data.frame(genes = unlist(genes)), by = "row.names")
#df_rounded <- df_rounded[, -1]
#fwrite(df_rounded, file = paste("fgsea/", arg,"_pathways_top.csv", sep = ""), sep =',')

