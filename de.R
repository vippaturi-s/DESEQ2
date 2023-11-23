#!/usr/bin/env Rscript
# de.R

library(tximport)
library(readr)
library(DESeq2)
library(knitr)
tx2gene <- read.csv("tx2gene.csv")
head(tx2gene)
samples <- read.csv("Samples.csv", header=TRUE)
head(samples)
files <- file.path("quant", samples$Sample, "quant.sf")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Menthol + Vibrio)
dds$Vibrio <- relevel(dds$Vibrio, ref = "Control")
dds$Menthol <- relevel(dds$Menthol, ref = "Control")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
padj <- .05
minLog2FoldChange <- .5
dfAll <- data.frame()
# Get all DE results except Intercept, and "flatten" into a single file.
for (result in resultsNames(dds)){
    if(result != 'Intercept'){
        res <- results(dds, alpha=.05, name=result)
        dfRes <- as.data.frame(res)
        dfRes <- subset(subset(dfRes, select=c(log2FoldChange, padj)))
        dfRes$Factor <- result
        dfAll <- rbind(dfAll, dfRes)
    }
}

dfAll_Filter <- subset(dfAll, padj<0.05)
dfAll_Filter1 <- data.frame(dfAll_Filter, row.names=NULL)
dfAll_bind <- cbind(rownames(dfAll_Filter),dfAll_Filter1)
colnames(dfAll_bind)[1] <- "ko"
kable(dfAll_bind) 

P_File <- "/home/vippaturi.s/BINF6309/module04-vippaturi-s-1/path.txt"
Path_ways <- read.table(P_File, sep="\t", header=FALSE)
colnames(Path_ways) <- c("ko", "Pathways")
ko_File <- "/home/vippaturi.s/BINF6309/module04-vippaturi-s-1/ko"
pathway <- read.table(ko_File, sep="\t", header=FALSE)
colnames(pathway) <- c("Pathways", "Description")
df_Pathway <- merge(Path_ways, pathway)
deAnnotated <- merge(df_Pathway, dfAll_bind)
kable(deAnnotated)
write.csv(deAnnotated, file="deAnnotated.csv")
