#Pavana Anur
#With contributions from Thuy Ngo and Breeshey Rosksams-Heiter
# Script to generate figure 3 and 4 for cfRNA manuscript
library(reshape2)
library(tsne)
library(ggplot2)
library(MASS)
library(cowplot)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggpubr)

# Set the working directory
setwd("../cfRNA/manuscript/figures/figure3/")
#############**** Import Data ******######################################################
data <- read.csv("../../tables/RPM_without_dates.csv",header=TRUE, sep=",")
metadata <- read.csv("../../tables/PP_metadata_keep_FINAL.csv",sep = ",", stringsAsFactors = FALSE)
colors <- read.delim("../../tables/colors.txt", stringsAsFactors = FALSE)

#########********** Format data *********************************#########################
data$gene=as.character(data$gene)

#Consolidate duplicate gene names
datam=melt(data,id="gene")
datac=dcast(datam,gene~variable,fun.aggregate = mean)
# Assign gene as rownames 
rownames(datac)=datac[,1]
#Remove first column, gene names
datac=datac[,-1]

######### PCA of genes selected by variance ################################################

# PCA function
GeneratePCA <- function(Data, Metadata, Baseline, Target1,Target2, BaselineName, TargetName1, TargetName2, Colors) {
  Metadata=Metadata[match(colnames(Data),Metadata$PP_ID),]
  sample_sel <- Metadata[Metadata$Status==Baseline|Metadata$Status==Target1|Metadata$Status==Target2,]$PP_ID
  datac_sel <- Data[,colnames(Data) %in% sample_sel]
  # Make genename become legal name in R
  names(datac_sel) <- make.names(names(datac_sel))
  # Transform data to log scale
  genecount <- log(datac_sel+1,2)
  # Which values in genecount are greater than 0.5?
  thresh <- genecount  > 1.5
  # This produces a logical matrix with TRUEs and FALSEs
  head(thresh)
  # Summary of how many TRUEs there are in each row
  table(rowSums(thresh))
  # we would like to keep genes that have at least 2 TRUES in each row of thresh
  keep <- rowSums(thresh) >= 5
  # Subset the rows of countdata to keep the more highly expressed genes
  counts.keep <- genecount[keep,]
  summary(keep)
  dim(counts.keep)
  # We estimate the variance for each row in the logcounts matrix
  var_genes <- apply(genecount, 1, var)
  head(var_genes)
  # Get the gene names for the top 500 most variable genes
  select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
  head(select_var)
  # Subset logcounts matrix
  highly_variable_lcpm <- genecount[select_var,]
  dim(highly_variable_lcpm)
  head(highly_variable_lcpm)
  Genes_wo_hit <- apply(highly_variable_lcpm,1,function(x) all(x==0))
  genecount_nozero <- highly_variable_lcpm[!Genes_wo_hit,]
  genecount_PCA <- t(genecount_nozero)
  
  # Run PCA
  PCA_out <- prcomp(log(genecount_PCA+1,2), scale. = TRUE)
  metadata_sel <- Metadata[Metadata$PP_ID %in% sample_sel,]
  rownames(metadata_sel)=metadata_sel$PP_ID
  metadata_sel=metadata_sel[,c(4,5)]
  Groups <- metadata_sel$Status
  # Generate the summary stats for prcomp
  eigs <- PCA_out$sdev^2
  summary <- rbind(
    SD = sqrt(eigs),
    Proportion = eigs/sum(eigs),
    Cumulative = cumsum(eigs)/sum(eigs))
  # Find the proprtion of PC1 and PC2
  propPC1 <- round(as.numeric(summary[2,1]) * 100)
  propPC2 <- round(as.numeric(summary[2,2]) * 100)
  # Find colours for baseline and target
  colourBaseline <- Colors[Colors$Status==Baseline,]$Colour
  colourTarget1 <- Colors[Colors$Status==Target1,]$Colour
  colourTarget2 <- Colors[Colors$Status==Target2,]$Colour
  
  # Plot PCA
  my.pca <- data.frame(PCA_out$x[,1:3], Groups)
  my.pca$Groups <- factor(my.pca$Groups, levels = c(Baseline, Target1,Target2))
  p <- ggplot(data = my.pca, mapping = aes(x = PC1,y = PC2, colour = Groups, size = 1)) +
    geom_point(size = 5) +
    theme_bw() +
    xlab(paste0("PC-1 (", propPC1, "%)")) +
    ylab(paste0("PC-2 (", propPC2, "%)")) +
    scale_color_manual(breaks=c(Baseline, Target1,Target2),
                       values=c(colourBaseline,colourTarget1,colourTarget2),
                       labels = c(BaselineName, TargetName1, TargetName2)) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          aspect.ratio=1,
          legend.title = element_blank()) 
  
  # Get legend
  legend <- cowplot::get_legend(p)
  
  pdf(paste0("Legend-", Baseline, "-vs-",Target1,"-vs-", Target2, ".pdf"), 5, 2)
  grid.newpage()
  grid.draw(legend)
  dev.off()
  
  pdf(paste0("PCA-top500-variance-", Baseline, "-vs-", Target1,"-vs-", Target2, ".pdf"))
  p <- p + theme(legend.position = "none")
  print(p)
  dev.off()
  
}

GeneratePCA(Data = datac, Metadata = metadata, Baseline = "HD", Target1 = "LCir", Target2 = "LCx",
            BaselineName = "Health Donor", TargetName1 = "Liver Cirrhosis",TargetName2 = "Liver Cancer", Colors = colors)
GeneratePCA(Data = datac, Metadata = metadata, Baseline = "HD", Target1 = "MG", Target2 = "MM",
            BaselineName = "Healthy Donor", TargetName1 = "MGUS",TargetName2 = "Multiple Myeloma ", Colors = colors)

################################################################################################
#boxplot for MG, MM and HD

#Select data
metadata=metadata[match(colnames(datac),metadata$PP_ID),]
sample_sel <- metadata[metadata$Status=="HD"|metadata$Status=="MG"|metadata$Status=="MM",]$PP_ID
metadata_sub=metadata[metadata$Status=="HD"|metadata$Status=="MG"|metadata$Status=="MM",]
genes=c("AIDA","CENPE","CPOX","ELL2","ELOVL6","GAB1","HBG1","MKI67","NEK2","NUSAP1")
datac_sel <- datac[rownames(datac) %in% genes,colnames(datac) %in% sample_sel]

# Make genename become legal name in R
names(datac_sel) <- make.names(names(datac_sel))

# Transform data to log scale
genecount <- log(datac_sel+1,2)

genecount$gene=rownames(genecount)
genecountm=melt(genecount,id="gene")
colnames(genecountm)=c("gene","PP_ID","counts")
metadata_group=metadata[,c(1,4)]
genecountm_group=merge(genecountm,metadata_group,by="PP_ID")

my_comparisons <- list(c("HD","MG"),c("MG","MM"),c("HD","MM"))

p=ggplot(genecountm_group, aes(x=Status, y=counts,color=Status)) + 
  geom_boxplot()+facet_wrap(~gene, scale="free",nrow=2)+geom_jitter(width = 0.1)+
  xlab("Status") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  ylab("Counts (log2(RPM+1))") +
  scale_colour_manual(values=c("#000000","#F375FF","#9000FF"))+
  stat_compare_means(method="t.test",comparisons = my_comparisons, label = "p.signif", size = 3)

fname <- "MM_boxplots.pdf"

pdf(fname, 10.5, 6.5)
print(p)
dev.off()


################################################################################################
#boxplots for LCir, LCx and HD

metadata=metadata[match(colnames(datac),metadata$PP_ID),]
sample_sel <- metadata[metadata$Status=="HD"|metadata$Status=="LCir"|metadata$Status=="LCx",]$PP_ID
metadata_sub=metadata[metadata$Status=="HD"|metadata$Status=="LCir"|metadata$Status=="LCx",]

genes=c("APOE","C3","CP","FGA","FGB","FGG","HRG","HULC","IFITM3","ORM1")
datac_sel <- datac[rownames(datac) %in% genes,colnames(datac) %in% sample_sel]

# Make genename become legal name in R
names(datac_sel) <- make.names(names(datac_sel))
# Transform data to log scale
genecount <- log(datac_sel+1,2)

genecount$gene=rownames(genecount)
genecountm=melt(genecount,id="gene")
colnames(genecountm)=c("gene","PP_ID","counts")
metadata_group=metadata[,c(1,4)]
genecountm_group=merge(genecountm,metadata_group,by="PP_ID")
my_comparisons <- list(c("HD","LCir"),c("LCir","LCx"),c("HD","LCx"))

p=ggplot(genecountm_group, aes(x=Status, y=counts,color=Status)) + 
  geom_boxplot()+facet_wrap(~gene, nrow=2,scale="free")+geom_jitter(width = 0.1)+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  ylab("Counts (log2(RPM+1))") +
  scale_colour_manual(values=c("#000000","#FFAE5C","#FF010A"))+
  stat_compare_means(method="t.test",comparisons = my_comparisons, label = "p.signif", size = 3)

fname <- "liver_boxplots.pdf"

pdf(fname, 10.5, 6.5)
print(p)
dev.off()


