#!/usr/bin/env Rscript
# mergeKo.R

library(knitr)

blastFile <- "/scratch/SampleDataFiles/Annotation/transcriptBlast.txt"
keggFile <- "/scratch/SampleDataFiles/Annotation/kegg.txt"
koFile <- "/scratch/SampleDataFiles/Annotation/ko.txt"

blast <- read.table(blastFile, sep="\t", header=FALSE)

colnames(blast) <- c("trans", "sp", "qlen", "slen", "bitscore", 
                     "length", "nident", "pident", "evalue", "ppos")
blast$cov <- blast$nident/blast$slen
blast <- subset(blast, cov > .5)
kable(head(blast))

kegg <- read.table(keggFile, sep="\t", header=FALSE)
colnames(kegg) <- c("sp", "kegg")

kegg$sp <- gsub("up:", "", kegg$sp)

kable(head(kegg))

blastKegg <- merge(blast, kegg)
kable(head(blastKegg))
ko <- read.table(koFile, sep="\t", header=FALSE)
colnames(ko) <- c("kegg", "ko")
kable(head(ko))
blastKo <- merge(blastKegg, ko)
kable(head(blastKo))
tx2gene <- unique(subset(blastKo, select=c(trans, ko)))
kable(head(tx2gene))
write.csv(tx2gene, file="tx2gene.csv", row.names=FALSE)
