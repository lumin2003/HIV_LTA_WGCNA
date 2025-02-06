library(WGCNA)
library(DESeq2)
library(tidyverse)
library(gridExtra)
library(stats)
library(ggplot2)
library(CorLevelPlot)
library(DESeq2)
library(dplyr)
library(ggthemes)
library(tidyr)
library(tidyverse)
library(dplyr)
library(stringr)
library(msigdb)
library(GSEABase)
library(tidyHeatmap)
library(data.table)



rawCountTable <- read.table("HIV1_HSC_htseq.txt", header = TRUE, sep = "\t", row.names = 1)

colnames(rawCountTable) <- c("CT1","CT2","CT3","gp1","gp2","gp3","Ba1","Ba2","Ba3","VT1","VT2","VT3")
spname <- colnames(rawCountTable[,-13])
batch <- substr(spname, start =3, stop =3)
treatment <- substr(spname,start = 1, stop=2)
group <-  substr(spname, start = 1, stop =2)
metadata <- data.frame(spname,treatment,group,batch)
rawCountTable <- rawCountTable[,-13]



# set up function to generate a DESeqDataSet
# using our htseq raw counts, sample annotations, and our model (which includes stage and Type)
library(DESeq2)
library(dplyr)
dds <- DESeqDataSetFromMatrix(countData = rawCountTable,
                              colData = metadata,
                              design= ~ treatment + batch)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds_norm <- vst(dds)
normalized_counts <- assay(dds_norm)

data <- data.matrix(normalized_counts)
length(which(grepl("CXCL8", rownames(data))))
##remove outlier
library(WGCNA)
gsg <- goodSamplesGenes(t(data))
data <- data[gsg$goodGenes==TRUE,]


htree <- hclust(dist(t(data)), method="average")
plot(htree)


dds <- DESeq(dds)
resultsNames(dds)



res <- results(dds, c("treatment","Ba","CT"))
resSig <- subset(res, pvalue < 0.05)
length(which(grepl("CXCL8", rownames(resSig))))


####Network Construction


data2 <- data.matrix(data)
data_norm <- t(data2)

library(doParallel)
# Detect the number of cores available
num_cores <- detectCores() - 1

# Create a cluster
cl <- makeCluster(num_cores)

# Register the parallel backend
registerDoParallel(cl)
library(WGCNA)
library(doParallel)
# Define power vector
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Run pickSoftThreshold with parallel processing enabled
sft <- pickSoftThreshold(data_norm, powerVector = power, networkType = "signed", verbose = 5)

# Stop the cluster after use
stopCluster(cl)

sft.data <- sft$fitIndices

library(gridExtra)
library(ggplot2)
a1 <- ggplot(sft.data, aes(power, SFT.R.sq, label= Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.8, color = "red") +
  labs(x="Power", y="Scale free topology model fit, signedR ^2")+
  theme_classic()


a2 <- ggplot(sft.data, aes(power, mean.k., label= Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.1) +
  labs(x="Power", y="Mean Connectivity")+
  theme_classic()

grid.arrange(a1,a2, nrow=2)


###momoery estimate w.r.r blocksize
# Choose a power based on the analysis above
norm.counts <- sapply(data_norm, as.numeric)
power = 22
temp_cor <- cor
cor <- WGCNA::cor

# Construct the network
bwnet = blockwiseModules(data_norm,
                         power = 22,
                       TOMType = "signed", 
                       masBlockSize = 14000,
                       mergeCutHeight = 0.25,
                       numericLabels = FALSE, 
                       randomSeed = 1234,
                       verbose = 3)
cor <- temp_cor
#Relabel blockwise modules

#2.B.2 Co-expression similarity and adjacency
softPower = 22;
adjacency = adjacency(data_norm, power = softPower);
#2.B.3 Adjacenct to Topological Overlap Matrix (TOM)
#Done to minimize effects of noise and spurious associations
TOM_tree = TOMsimilarity(adjacency);
dissTOM = 1-TOM_tree
#2.B.4 Clustering using TOM: produces dendrogram of genes
geneTree = hclust(as.dist(dissTOM), method = "average");
#Plot the dendogram
sizeGrWindow(12,9) 
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04);
#Set module size relatively high
minModuleSize = 30;


# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)
MEDissThres =0.25
merge = mergeCloseModules(data_norm, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
moduleLabels = bwnet$colors
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;



## module eigenenes
module_eigengenes <- bwnet$MEs

nSmaples <- nrow(data_norm)
nGenes <- ncol(data_norm)


gourp.out <- binarizeCategoricalColumns(metadata$treatment, includePairwise = TRUE, includeLevelVsAll = TRUE,minCount = 1)
traits <- gourp.out
library(dplyr)
traits <- metadata %>% mutate(viral_state_bin = ifelse(grepl('Ba',treatment),1,0)) %>% dplyr::select(5)
traits <- cbind(traits,gourp.out)
row.names(traits) <- metadata$spname

library(tidyr)
library(tidyverse)
module_trait_corr <- cor(module_eigengenes, traits, use ="p")
module_trait_corr_pvals <- corPvalueStudent(module_trait_corr, nSmaples)
heatmap.data <- merge(module_eigengenes, traits, by = "row.names")                              
heatmap.data <- heatmap.data %>% column_to_rownames(var='Row.names')
names(heatmap.data)
library(CorLevelPlot)


#show correlation value
CorLevelPlot(heatmap.data,
            x=names(heatmap.data)[22:30],
            y=names(heatmap.data)[1:21],
            col=c("blue1","skyblue","white","pink","red"))

####enrich p value sig gene; Sig_nodule for output
module.gene.mapping <- as.data.frame(bwnet$colors)

Sig_nodule1 <-  module.gene.mapping %>% dplyr::filter(`bwnet$colors`=="salmon"  ) 
Sig_nodule2 <-  module.gene.mapping %>% dplyr::filter(`bwnet$colors`=="green"  ) 
Sig_nodule3 <-  module.gene.mapping %>% dplyr::filter(`bwnet$colors`=="magenta"  ) 
Sig_nodule4 <-  module.gene.mapping %>% dplyr::filter(`bwnet$colors`=="pink" )
Sig_nodule5 <-  module.gene.mapping %>% dplyr::filter(`bwnet$colors`=="purple" )
Sig_nodule6 <-  module.gene.mapping %>% dplyr::filter(`bwnet$colors`=="cyan" )
Sig_nodule7 <-  module.gene.mapping %>% dplyr::filter(`bwnet$colors`=="lightgreen" )
Sig_nodule8 <-  module.gene.mapping %>% dplyr::filter(`bwnet$colors`=="grey60" )
Sig_nodule9 <-  module.gene.mapping %>% dplyr::filter(`bwnet$colors`=="yellow" )
Sig_nodule10 <-  module.gene.mapping %>% dplyr::filter(`bwnet$colors`=="greenyellow" )

#############for VT Vs Ba 
Sig_nodule <- bind_rows(Sig_nodule1,Sig_nodule2,Sig_nodule3,Sig_nodule4,Sig_nodule5,Sig_nodule6,Sig_nodule7,Sig_nodule8,Sig_nodule9,Sig_nodule10,IL8_nodule ) 



#### module Vs Gene p value (the correlation between genes and module)
module.membership.measure <- cor(module_eigengenes, data_norm, use='p')
module.membership.measure.pvalue <- corPvalueStudent(module.membership.measure, nSmaples)



### gene Vs treatment p value (The correlation of Genes with treatment)
 gene.signf.corr <- cor(data_norm, traits$data.VT.vs.Ba, use = 'p')
 gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSmaples)
 gene.signf.corr.pvals %>% as.data.frame() %>% arrange(V1) %>% head()

###gene.signf.corr.pvals of CXCL8 is 0.03982587, meaning VT increased IL8 comparedd to BaL infection
gene.signf.corr.pvals <- as.data.frame(gene.signf.corr.pvals)
gene_filter <- dplyr::filter(gene.signf.corr.pvals, V1 < 0.05)

length(which(grepl("HDAC1", gene_filter )))

###### select interesed module and export edge_list for Cytoscape, VisANT
expr_normalized <- data2
shared_genes <- intersect(rownames(gene_filter),rownames(Sig_nodule))
expr_of_interest <- expr_normalized[shared_genes,]

TOM3 = TOMsimilarityFromExpr(t(expr_of_interest))
row.names(TOM3) = row.names(expr_of_interest)
colnames(TOM3) = row.names(expr_of_interest)




# Select modules
modules = c("green", "pink","salmon","magenta","grey60","greenyellow","yellow","lightgreen","cyan","purple","grey");

# Export the network into edge and node list files for Cytoscape

cyt = exportNetworkToCytoscape(TOM3,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = row.names(expr_of_interest),
                               altNodeNames = NULL,
                               nodeAttr = NULL,
                               includeColNames = TRUE);




