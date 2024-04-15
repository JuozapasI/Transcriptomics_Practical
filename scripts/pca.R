setwd("~/Studijos/Transcriptomics/")
library(ggfortify)
library(ggrepel)

arg <- commandArgs(trailingOnly = TRUE)

pattern <- ".*\\.(HBR|UHRR)\\..*S(\\d+)_.*"

cts_collibri <- read.csv(arg[1], sep = '\t', row.names = "Geneid", comment.char = '#')
cts_kapa <- read.csv(arg[2], sep = '\t', row.names = "Geneid", comment.char = '#')

colnames(cts_collibri) <- sub(pattern, "\\1_S\\2", colnames(cts_collibri))
colnames(cts_kapa) <- sub(pattern, "\\1_S\\2", colnames(cts_kapa))
cts_collibri <- cts_collibri[, 6:9]
cts_kapa <- cts_kapa[, 6:9]

cts = merge(cts_collibri, cts_kapa, by = "row.names", all = FALSE)
rownames(cts) = cts$Row.names
cts <- cts[, -1]

de_collibri <- read.csv(arg[3], header = TRUE, row.names = 1)
de_kapa <- read.csv(arg[4], header = TRUE, row.names = 1)

de_genes <- intersect(rownames(de_collibri), rownames(de_kapa))

cts <- cts[rownames(cts) %in% de_genes, ]

pca_cts <- prcomp(t(cts))

png("pca/pca.png")
autoplot(pca_cts) +
geom_text_repel(aes(label = colnames(cts)), size = 3)
dev.off()

