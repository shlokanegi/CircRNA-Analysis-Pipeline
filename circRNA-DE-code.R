#---- circRNA Expression Analysis ----

# Set the working directory
setwd("~/data/analyses")

# Install and Load Required Packages
library("BiocManager")
library("dplyr")
library("tidyverse")
library("plyr")
library("data.table")
library("DESeq2")

# Collect circRNA detection files for all samples.
Samples <- c("SRR6674618", "SRR6674619", "SRR6674620", 
             "SRR6674622", "SRR6674623", "SRR6674624")

# Create an input list of CIRI outputs
inputs <- list()
for (i in 1:length(Samples)) {
  inputs[[i]] <- paste0("../ciri_out/", Samples[i], "_ciri.out")
}
names(inputs) <- Samples
head(inputs)
View(inputs)

# Form a combined CIRI output dataframe from all samples
combined.df <- ldply(inputs, function(x) {
  a <- read.table(file=x, sep="\t", header=T, comment.char="", 
                  stringsAsFactors=F) [,c(1:11)]; a})
colnames(combined.df)[1] <- "sampleID"
head(combined.df)
View(combined.df) # 36398 rows, 12 columns

coldata <- data.frame(sample_id=Samples, condition=c("B","B","B","M","M","M"))
rownames(coldata) <- coldata$sample_id
coldata$condition <- factor(coldata$condition)

# Making rmat -> count data matrix for all circRNAs across all samples
unicID = sort(unique(combined.df$circRNA_ID))
rmat = matrix(NA, nrow = length(unicID), ncol = 6, dimnames = list(unicID, Samples))

# Checking if all circRNA_IDs are unique in the ciri.out output files
for(s in Samples) {
  cID = combined.df[combined.df$sampleID==s,"circRNA_ID"] 
  print(paste(s,length(cID),length(unique(cID))))
}

# Filling rmat with count data (junction reads extracted from combined.df)
for(s in Samples) {
  m <- match(rownames(rmat), combined.df[combined.df$sampleID==s,"circRNA_ID"])
  rmat[,s] = combined.df[m, "X.junction_reads"]
}
View(rmat)

# Removing NA values (rows containing NA values)
counts <- na.omit(rmat)
View(counts)

# Implementing DESeq2 --> Column names of the count matrix must match the order of the coldata rownames.
#------------------------------------------------------------------
# Replace NA values with 0s
counts <- data.frame(dcast(circ.circ.xpr, formula = circ_id ~ sample_id, value.var="circ.reads", fill=0), row.names="circ_id") 

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
keep <- rowSums(counts(dds)) >= 5 #Filtering out circRNAs with very low counts
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)
res

# Visualizing results
#----------------------
summary(res) # 133 circRNAs are upregulated
             # 821 circRNAs are downregulated
View(res)
res_pval <- data.frame("circIDs" = res@rownames, "pvalue" = res$pvalue, "padj" = res$padj, "Log2FC" = res$log2FoldChange)
res_pval <- res_pval[order(res_pval$padj), ]

length(which(res_pval$padj < 0.05))
sum(res$padj < 0.05, na.rm=TRUE) # 622 circRNAs are significantly differentially expressed between B and M conditions.

#BiocManager::install("apeglm")
resLFC <- lfcShrink(dds, coef="condition_M_vs_B", type="apeglm")

plotMA(res)
plotMA(resLFC)

# Compare counts between B-Cell and Monocyte groups for our top 6 circRNAs
par(mfrow=c(2,3))
plotCounts(dds, gene = "4:143543509|143543972", intgroup="condition")
plotCounts(dds, gene = "5:123545417|123557564", intgroup="condition")
plotCounts(dds, gene = "5:97101759|97103094", intgroup="condition")
plotCounts(dds, gene = "8:51860845|51861246", intgroup="condition")
plotCounts(dds, gene = "5:137985257|137988315", intgroup="condition")
plotCounts(dds, gene = "4:152411303|152412529", intgroup="condition")


## VOLCANO PLOTS
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim = c(-13,13)))
# Add colored points: blue if padj<0.05, red if log2FC>7 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>7), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

### Volcano Plot using ggplot
Link - https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
# add a column of NAs
res_pval$diffexpressed <- "NO"
# if log2Foldchange > 7 and padj < 0.05, set as "UP" 
res_pval$diffexpressed[res_pval$Log2FC > 7 & res_pval$padj < 0.05] <- "UP"
# if log2Foldchange < -7 and pvalue < 0.05, set as "DOWN"
res_pval$diffexpressed[res_pval$Log2FC < -7 & res_pval$padj < 0.05] <- "DOWN"

p <- ggplot(data=res_pval, aes(x=Log2FC, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-7, 7), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2

res_pval$delabel <- NA
res_pval$delabel[res_pval$diffexpressed != "NO"] <- res_pval$circIDs[res_pval$diffexpressed != "NO"]
ggplot(data=res_pval, aes(x=Log2FC, y=-log10(pvalue), col=diffexpressed, label=delabel)) + geom_point() + theme_minimal() + geom_text()

library(ggrepel)
p3 <- ggplot(data=res_pval, aes(x=Log2FC, y=-log10(pvalue), col=diffexpressed, label=delabel)) + geom_point() + theme_minimal() + geom_text_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 10)) + scale_color_manual(values=c("blue", "black", "red")) + geom_vline(xintercept=c(-7, 7), col="red") + geom_hline(yintercept=-log10(0.05), col="red")
p3 

length(which(res_pval$diffexpressed == "DOWN")) #153 circRNAs downregulated in M, as compared to B
length(which(res_pval$diffexpressed == "UP")) #2 circRNAs upregulated in M, as compared to B
length(which(res_pval$diffexpressed == "NO")) #9649 circRNAs not significant

#### PCA
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="condition") #PC1 and PC2 explaining 90% variance
#look at how our samples group by conditions


#--------------------EXTRAS/NOT NEEDED--------------------------
'BiocManager::install("Rsamtools")
library("Rsamtools")
BiocManager::install("GenomicFeatures")
library("GenomicFeatures")
BiocManager::install("GenomicAlignments")
library("GenomicAlignments")
library(BiocParallel)

filenames <- file.path(paste0("../bwa_out/", Samples, "_bwa.bam"))
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])

gtffile <- file.path("../Homo_sapiens.GRCh38.100.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf")

ebg <- exonsBy(txdb, by="gene")
#ebg1 <- transcriptsBy(txdb, by="gene")

detectCores(all.tests = FALSE)
register(MulticoreParam(3))

se <- summarizeOverlaps(features=ebg, 
                        reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE)
class(assay(se))
combined.df <- as.data.frame(assay(se), stringsAsFactors=F)

#Extras
length(which(temp$circRNA_type!="exon"))
tempgene = separate_rows(temp, gene_id, sep = ",", convert = FALSE)
View(tempgene)
length(unique(tempgene$gene_id))
SRR18 <- unique(tempgene$gene_id)
head(SRR18)
which(SRR18 == "" | SRR18 == "n/a")
ixrm = which(SRR18 == "" | SRR18 == "n/a")
if (length(ixrm) > 0){
  SRR18 <- SRR18[-ixrm]
}

#---------------------------------------------------------------------------
SRR20_counts <- read.delim("~/data/SRR6674620_bks_linear_counts.tab", header=FALSE)
View(SRR20_counts)
temp20 <- read.delim("~/data/temp20", header=FALSE)
View(temp20)
#---------------------------------------------------------------------------'''
