# CircRNA-Analysis-Pipeline

# Project Overview

## Implementation and Optimization of Circular RNA Analysis Bioinformatics Pipeline: Preprocessing, Alignment, Detection, Differential Expression and Functional Annotation

### Overview
Circular RNAs (CircRNAs) are a subclass of lncRNAs, which are transcribed and backspliced using the same cellular machinery as that by linear RNAs. In the past few years, circRNAs have been studied with the aim of identifying new and possibly targetable disease mechanisms, including cancer. The suitability of circRNA as a diagnostic and prognostic marker needs to be assessed. 

Using bioinformatics analysis, we have tried to study circRNAs through backsplice region identification. In this project, we implement and optimize a computational protocol for circRNA detection and quantification from RNA-Seq data. Next, we test for circRNA differential expression and study the relationship between linear and circular transcript expression. Finally, we reconstruct the sequence of circRNAs and carry out in-silico circRNA functional characterization.

![Picture1](https://user-images.githubusercontent.com/66521525/135354639-799dfbf8-b60e-4da1-8d2c-ff5159f855e3.png)


## IMPLEMENTATION is divided into the following parts :-
#### 1) CircRNA Detection
#### 2) CircRNA Expression Analysis
#### 3) CircRNA Characterization

## RESULTS and CONCLUSION
### 1. Detection of CircRNAs with backsplice junction coordinates
CIRI2 tool was used to detect the backspliced junctions in the samples. We provided the bwa chimeric alignment SAM file to the tool. 10735 circRNAs were identified in sample SRR6674618, 9247 circRNAs were identified in sample SRR6674619, 4351 circRNAs were identified in sample SRR6674620, 4805 circRNAs were identified in sample SRR6674622, 4120 circRNAs were identified in sample SRR6674623 and 3140 circRNAs were identified in sample SRR6674624. Table 1 represents the head of the SRR6674618_ciri.out table.

### 2. Detection of Collinear Read Counts on Backspliced Junctions
bks_linear_counts.tab file was created as an output of bedtools coverage tool. We provided the hisat2 collinear alignment BAM file to the tool. The output file consisted of the counts of linear transcripts which flanked the identified backspliced junction of left and right sides. Total linear transcripts can be calculated by adding the left and right flanked counts.

### 3. 662 CircRNAs were found to be significantly differentially expressed between B-Cell and Monocyte conditions
We identified 133 circRNAs showing upregulated expression in Monocyte samples as compared to B-cell samples, whereas 821 circRNAs were downregulated in Monocyte samples as compared to B-cell samples. Setting the Padj threshold to 0.05, we identified 662 circRNAs which were differentially expressed between Monocyte and B-Cell conditions. circRNA 4:143543509|143543972 was found to be most significant (Padj = 3.14E-11) amongst all circRNAs, and had a log2FC value of -12.25.

![Picture1](https://user-images.githubusercontent.com/66521525/135355460-5f0e857f-72cb-491c-9496-cadcf73f29bd.png)

### 4. Two CircRNAs were Upregulated and 153 CircRNAs were Downregulated in Monocyte conditions as compared to B- Cell conditions
Violin plot represents circRNAs which were showed a significant differential expression between conditions and had a |log2 fold change| > 7. As shown, only 2 circRNAs were significantly upregulated, one of them belonging to the sex- linked chromosome. On the contrary, 153 circRNAs were significantly downregulated in monocyte conditions as compared to B-cell conditions.

![Picture1](https://user-images.githubusercontent.com/66521525/135355507-7fac8ff3-e14f-4260-80b6-d41ffe59626a.png)

### 5. PC1 and PC2 covered 90% of the variation in the data
PCA plot shows how well the samples group by conditions. Here, PC1 covers 80% or the variation in the data, whereas PC2 covers 10% of the variation in the data. Together, 900% variance is being explained.

![Picture1](https://user-images.githubusercontent.com/66521525/135355568-a5aaa00f-44fb-41ce-9fd0-e7e34ca32ffe.png)

### 6. 43 CircRNAs were found to be represent significant diversion in the expression of their circular and linear transcripts. 
To discern which circular and linear RNA isoform expressions diverge, counts of the reads that are linearly spliced on the backsplice junctions are used as an estimate of the parental gene expression (from the bks_linear_counts.tab.gz file, which had combined linear counts of all samples. 15:85113873|85115487 was identified as the most significant circRNA (p- value = 0.0089). The top 6 significant circRNAs and the abundance of host- gene and circular expression. In all cases, the overall circular reads, mapping to a given backsplice junction are lower than the linear reads.

![Picture2](https://user-images.githubusercontent.com/66521525/135355605-760dd463-d8e8-4924-8979-39abd7b308a5.png)
