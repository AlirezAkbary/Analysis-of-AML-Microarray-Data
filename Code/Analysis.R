####setting working directory, Please change the directory for running in your computer
setwd("/Users/alireza/Desktop/University/97-98\ Spring/Bioinformatics/Project/R")

####Loading needed libraries
library(Biobase)
library(GEOquery)
library(limma)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
library(plotly)
library(magrittr)
library(pheatmap)

####load data
series <- "GSE48558"
platform <- "GPL6244"
gset <- getGEO(series, GSEMatrix = T, AnnotGPL = T, destdir = "Data/")
if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- exprs(gset) #extract expression matrix 


####making grouping of samples according to dataset

##gr2 have AMLP(based on source name), Normal(based on phenotype), other(based on source name)
gr2 <- c(rep('AMLP', 13), rep("B_ALL", 4), rep("T_ALL", 2), "B_ALL", "T_ALL","B_ALL", "B_ALL", "T_ALL", "B_ALL", "B_ALL", rep("T_ALL", 2), rep("B_ALL", 5), "T_ALL", "T_ALL", "B_ALL", "T_ALL", "BP", "T_ALL","AMLCL", "Normal", "BP", "T_ALL", "AMLCL", "Normal", "BP", "T_ALL", "AMLCL", "BP", rep("AMLCL", 2), "BP",rep("AMLCL", 2), "BP", rep("AMLCL", 2), "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "Normal", "B_ALL", "Normal", "AMLCL", "BP", "B_ALL", "Normal", "B_ALL", "Normal", "Normal", "Normal", "Normal", "B_ALL", "Normal", "AMLCL", "B_ALL", rep("Normal", 2), "AMLCL", "B_ALL", rep("Normal", 2), "AMLCL", "Normal", "B_ALL", "Normal", "AMLCL", "Normal", "B_ALL", "Normal", "AMLCL", "Normal", "T_ALL", "TP","AMLCL", "Normal", "T_ALL", "TP","AMLCL", "Normal", "T_ALL", "TP", "AMLCL", "TP", rep("BP", 9), "TP", rep("BP", 7), rep("TP", 8), rep("Normal", 7), "AMLP", "AMLP", "Normal", rep("AMLP", 3), rep("Normal", 7), "Normal", rep("Normal", 4), "Normal", rep("Normal", 7))

##gr3 have AMLP(based on source name), Normals that divided into subgroup based on source name, other(based on source name)
gr3 <- c(rep('AMLP', 13), rep("B_ALL", 4), rep("T_ALL", 2), "B_ALL", "T_ALL","B_ALL", "B_ALL", "T_ALL", "B_ALL", "B_ALL", rep("T_ALL", 2), rep("B_ALL", 5), "T_ALL", "T_ALL", "B_ALL", "T_ALL", "BP", "T_ALL","AMLCL", "Granul_Norm", "BP", "T_ALL", "AMLCL", "Granul_Norm", "BP", "T_ALL", "AMLCL", "BP", rep("AMLCL", 2), "BP",rep("AMLCL", 2), "BP", rep("AMLCL", 2), "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "B_Cell_Norm", "B_ALL", "T_Cell_Norm", "AMLCL", "BP", "B_ALL", "Granul_Norm", "B_ALL", "Granul_Norm", "Mono_Norm", "Mono_Norm", "B_Cell_Norm", "B_ALL", "T_Cell_Norm", "AMLCL", "B_ALL", rep("T_Cell_Norm", 2), "AMLCL", "B_ALL", rep("T_Cell_Norm", 2), "AMLCL", "B_Cell_Norm", "B_ALL", "T_Cell_Norm", "AMLCL", "B_Cell_Norm", "B_ALL", "T_Cell_Norm", "AMLCL", "CD34_Norm", "T_ALL", "TP","AMLCL", "CD34_Norm", "T_ALL", "TP","AMLCL", "CD34_Norm", "T_ALL", "TP", "AMLCL", "TP", rep("BP", 9), "TP", rep("BP", 7), rep("TP", 8), rep("Granul_Norm", 7), "AMLP", "AMLP", "T_Cell_Norm", rep("AMLP", 3), rep("B_Cell_Norm", 7), "T_Cell_Norm", rep("Mono_Norm", 4), "Granul_Norm", rep("T_Cell_Norm", 7))




###these are used when we use "just" AMLP type and Normal(divided by their source name) samples.
df.ex <- as.data.frame(ex)
df.ex <- df.ex[,(gr3 == 'AMLP') | (grepl('Norm', gr3))]
gr_reduced2 <- gr3[(gr3 == "AMLP") | (grepl('Norm', gr3))]

################################################################################
#####quality control

##Checking Log Scale
dim(ex)
max(ex) #we see that max is 13, and it is not too high.So it does not need any change


##Checking Normalization of Sample
pdf("Results/boxplot.pdf", width = 15)
boxplot(ex)
dev.off()



##Correlation Heatmap
pdf("Results/CorHeatmap.pdf", width = 20, height = 20)
pheatmap(cor(ex), labels_row = gr2, labels_col = gr2,  color=bluered(256), border_color = NA)
pheatmap(cor(df.ex), labels_row = gr_reduced2, labels_col = gr_reduced2 ,color=bluered(256), border_color = NA)
dev.off()



##PCA(on genes)
pc <- prcomp(ex)

pdf("Results/PCA.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

##PCA applied on difference of each sample from its average(on genes)

ex.scale <- t(scale(t(ex), scale = F))
pc <- prcomp(ex.scale)

df.ex.scale <- t(scale(t(df.ex), scale = F))#df.ex
df.pc <- prcomp(df.ex.scale)#df.ex

pdf("Results/PCA_scaled.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

##PCA on samples

pcr <- data.frame(df.pc$rotation[, 1:3], Group = gr_reduced2)

pdf("Results/PCA_samples.pdf")
ggplot(pcr, aes(PC1, PC2, color= Group)) + geom_point(size=3) + theme_bw()
dev.off()

#ploting 3d
plot_ly(x= pcr$PC1, y= pcr$PC2, z = pcr$PC3, mode="markers", color = pcr$Group)%>%
  layout(
    title = "AML Patient and Normal samples on PC1, PC2, PC3",
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    ))

#end of quality control, if grouping of sample were good, we can rely on data
################################################################################


####Differntial Expression Analysis

#gr4 is groupping of samples and it is like gr2, which I descipe at top.
gr4 <- c(rep('AMLP', 13), rep("B_ALL", 4), rep("T_ALL", 2), "B_ALL", "T_ALL","B_ALL", "B_ALL", "T_ALL", "B_ALL", "B_ALL", rep("T_ALL", 2), rep("B_ALL", 5), "T_ALL", "T_ALL", "B_ALL", "T_ALL", "BP", "T_ALL","AMLCL", "Normal", "BP", "T_ALL", "AMLCL", "Normal", "BP", "T_ALL", "AMLCL", "BP", rep("AMLCL", 2), "BP",rep("AMLCL", 2), "BP", rep("AMLCL", 2), "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "Normal", "B_ALL", "Normal", "AMLCL", "BP", "B_ALL", "Normal", "B_ALL", "Normal", "Normal", "Normal", "Normal", "B_ALL", "Normal", "AMLCL", "B_ALL", rep("Normal", 2), "AMLCL", "B_ALL", rep("Normal", 2), "AMLCL", "Normal", "B_ALL", "Normal", "AMLCL", "Normal", "B_ALL", "Normal", "AMLCL", "Normal", "T_ALL", "TP","AMLCL", "Normal", "T_ALL", "TP","AMLCL", "Normal", "T_ALL", "TP", "AMLCL", "TP", rep("BP", 9), "TP", rep("BP", 7), rep("TP", 8), rep("Normal", 7), "AMLP", "AMLP", "Normal", rep("AMLP", 3), rep("Normal", 7), "Normal", rep("Normal", 4), "Normal", rep("Normal", 7))

gr4 <- factor(gr4)
gset$description <- gr4
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(gr4)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(AMLP - Normal, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
colnames(tT)
tT <- subset(tT, select=c("Gene.symbol", "Gene.ID","adj.P.Val", "logFC"))
write.table(tT, "Results/AMLP_Normal.txt", row.names=F, sep="\t", quote = F)


####Pathway & gene anthology Analysis

##Making Gene wiht increase of expression
aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(aml.up$Gene.symbol)
aml.up.genes <- unique(as.character(strsplit2(aml.up.genes, "///")))
write.table(aml.up.genes, "Results/AMLP_Normal_up.txt",quote = F, row.names = F, col.names = F)

##Making Gene wiht decrease of expression
aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(aml.down$Gene.symbol)
aml.down.genes <- unique(as.character(strsplit2(aml.down.genes, "///")))
write.table(aml.down.genes, "Results/AMLP_Normal_down.txt", row.names = F, col.names = F, quote = F)

################################################################################


#####Gene Set enrichment analysis
ex.df <- as.data.frame(ex)
ex.df <- cbind(rownames(ex), ex.df)
colnames(ex.df) <- c("Name", as.character(gr2)) 
rownames(ex.df) <- NULL
write.table(ex.df, "Results/gsea.txt", row.names = F, col.names = T, sep = "\t", quote = F)
