setwd("~/Studijos/Transcriptomics/")
library(DESeq2)
library(ggplot2)

arg <- commandArgs(trailingOnly = TRUE)[1]

pattern <- ".*\\.(HBR|UHRR)\\..*S(\\d+)_.*"

cts <- read.csv(paste("counts/", arg,".txt", sep=""),
                sep = '\t',
                row.names = "Geneid",
                comment.char = '#')

colnames(cts) <- sub(pattern, "\\1_S\\2", colnames(cts))

cts <- cts[, 6:9]

coldata <- data.frame(row.names = colnames(cts))
coldata$condition <- ifelse(grepl("^HBR", rownames(coldata)), "healthy", "cancer")

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds$condition <- factor(dds$condition, levels = c("healthy", "cancer"))
dds <-DESeq(dds)
res <- results(dds)
resLFC <- lfcShrink(dds, coef="condition_cancer_vs_healthy", type="apeglm")
resLFC <- resLFC[order(resLFC$pvalue), ]

res <- res[order(res$pvalue), ]
write.csv(as.data.frame(resLFC), 
          file=paste("deseq/",arg,".csv", sep =""))

#plotMA(res, ylim=c(-10,10))
plotMA(resLFC, ylim=c(-10,10))

png(file=paste("deseq/",arg,".png", sep =""))
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
ggplot(data.frame(res)[,c("log2FoldChange", "padj")], aes(x = log2FoldChange, y = -log10(padj))) +
     geom_point(alpha = 0.5, size = 1) +
     geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
     geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
     labs(x = "Log2 Fold Change", y = "-log10(p-value)", title = "Volcano Plot")
dev.off()


