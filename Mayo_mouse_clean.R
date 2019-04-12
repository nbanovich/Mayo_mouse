#### Import Libs ####
library(Seurat)
library(edgeR)
library(gplots)
library(dplyr)
library(reshape2)
library(cellrangerRkit)
library(ggplot2)
library(scde)
library(monocle)
#library(flashClust)
library(WGCNA)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(topGO)
library(scatterplot3d)
library(plotly)
library(AUCell)
library(RcisTarget)
library(SCENIC)
library(GENIE3)
library(doRNG)
options(stringsAsFactors=FALSE)

#### Data Prep and plots ####

for(i in 1:16){
  tmp =  paste("data", i, sep = "_")
  id1 =  0+i
  id2 = 16623 +i
  assign(tmp,  read.table(gzfile(paste("Data/Mayo_Mouse/MayoMouse_JJK_p0_CL_Whole_T",id1,"_X3SCR_K",id2,"_raw.txt.gz", sep = ""))))
}



for(i in 1:4){
  obj = 16 + i
  tmp =  paste("data", obj, sep = "_")
  id1 =  0+i
  id2 = 17077 +i
  assign(tmp,  read.table(gzfile(paste("Data/Mayo_Mouse/MayoMouse_JJK_p0_CL_Whole_T",id1,"_XCSCR_K",id2,"_raw.txt.gz", sep = ""))))
  }

all_list = list(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13,data_14, data_15, data_16, data_17, data_18, data_19, data_20)

for(i in 1:20){
  tmp = paste("data_object",i, sep = "_")
  assign(tmp, CreateSeuratObject(all_list[[i]], min.cells = 0, min.genes = 0, project = "MM_mouse"))
}


# Combine all data into a list

all_list = list(data_object_1, data_object_2, data_object_3, data_object_4, data_object_5, data_object_6, data_object_7, data_object_8, data_object_9, data_object_10, data_object_11, data_object_12, data_object_13,data_object_14, data_object_15, data_object_16, data_object_17, data_object_18, data_object_19, data_object_20)

# Remove individual data objects
rm(data_object_1, data_object_2, data_object_3, data_object_4, data_object_5, data_object_6, data_object_7, data_object_8, data_object_9, data_object_10, data_object_11, data_object_12, data_object_13,data_object_14, data_object_15, data_object_16, data_object_17, data_object_18, data_object_19, data_object_20)

rm(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13,data_14, data_15, data_16, data_17, data_18, data_19, data_20)

# Read in meta.data
metadata = read.table("Data/Mayo_Mouse/Mayo_metadata.txt", header = T, sep = "\t")

# Add treatment to the meta.data matrix


for(i in 1:20){
  all_list[[i]]@meta.data$treatment <- metadata[i,3]
}

# Change the orig.ident to the proxy number of the individual
for(i in 1:20){
  all_list[[i]]@meta.data$orig.ident <- i
}


## Plot UMI of raw data per cell by treatment

# Get nUMI for all samples
tmp = lapply(all_list, function(xx) return(xx@meta.data$nUMI))

# Make melted df
tmp_df = melt(tmp)

# Get treatment id from each mouse
for(i in 1:20){
  tmp_df[tmp_df$L1 == i,]$L1 <- as.character(metadata[i,]$Treatment)
}

# Generate violin plot from nUmi data
raw_violin = ggplot(tmp_df) +  geom_violin(aes(x = as.factor(L1), y = value, fill = as.factor(L1)), scale = "width", trim = F ) + geom_boxplot(aes(x = as.factor(L1), y = value), col = "black", width=.15, outlier.colour = NA) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ggtitle("nUMI")


#### Make all treatments ####

# Make a large object containing samples from all treatments - everything but last two controls
tmp = MergeSeurat(all_list[[1]], all_list[[2]], do.normalize = F, scale.factor = 1e4, add.cell.id1 = "cont1", add.cell.id2 = "cont2", do.scale = F, do.center = F)
#tmp = MergeSeurat(all_list[[1]], all_list[[4]], do.normalize = F, scale.factor = 1e4, add.cell.id1 = "cont1", add.cell.id2 = "untreat1", do.scale = F, do.center = F)

#tmp = MergeSeurat(tmp, all_list[[3]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "cont3", do.scale = F, do.center = F)

#tmp = MergeSeurat(tmp, all_list[[4]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "untreat1", do.scale = F, do.center = F)

tmp = MergeSeurat(tmp, all_list[[5]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "untreat1", do.scale = F, do.center = F)

tmp = MergeSeurat(tmp, all_list[[6]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "untreat2", do.scale = F, do.center = F)

tmp = MergeSeurat(tmp, all_list[[7]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "lcl1", do.scale = F, do.center = F)

#tmp = MergeSeurat(tmp, all_list[[8]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "lcl2", do.scale = F, do.center = F)

tmp = MergeSeurat(tmp, all_list[[9]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "lcl3", do.scale = F, do.center = F)

tmp = MergeSeurat(tmp, all_list[[10]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "cyc1", do.scale = F, do.center = F)

tmp = MergeSeurat(tmp, all_list[[11]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "cyc2", do.scale = F, do.center = F)

tmp = MergeSeurat(tmp, all_list[[12]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "bor1", do.scale = F, do.center = F)

tmp = MergeSeurat(tmp, all_list[[13]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "bor2", do.scale = F, do.center = F)

#tmp = MergeSeurat(tmp, all_list[[14]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "cont4", do.scale = F, do.center = F)

#tmp = MergeSeurat(tmp, all_list[[15]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "cont5", do.scale = F, do.center = F)

#tmp = MergeSeurat(tmp, all_list[[16]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "cont6", do.scale = F, do.center = F)

tmp = MergeSeurat(tmp, all_list[[17]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "lcl4", do.scale = F, do.center = F)

tmp = MergeSeurat(tmp, all_list[[18]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "untreat3", do.scale = F, do.center = F)

tmp = MergeSeurat(tmp, all_list[[19]], do.normalize = F, scale.factor = 1e4, add.cell.id2 = "lcl_cont", do.scale = F, do.center = F)

all_treatments = MergeSeurat(tmp, all_list[[20]], min.cells = 10, min.genes = 200, do.normalize = T, scale.factor = 1e4, add.cell.id2 = "cont3", do.scale = F, do.center = F)

#all_treatments = MergeSeurat(tmp, all_list[[14]], min.cells = 3, min.genes = 100, do.logNormalize = T, total.expr = 1e4, add.cell.id2 = "bor", do.scale = T, do.center = T)

rm(all_list)


#### Data clean up and initial analysis ####


# Generate MeanVarPlot data (identifies variable genes which can be used later)
all_treatments <- FindVariableGenes(all_treatments, x.low.cutoff = 0.0125, y.cutoff = 1.5, do.contour = F, do.plot = T)

# Compute PCs from the data - NOTE must run this way because the object is too large to scale data on my computer
all_treatments <- ScaleData(all_treatments, genes.use = all_treatments@var.genes)

all_treatments <- RunPCA(all_treatments, pc.genes = all_treatments@var.genes)

# Find optimal # of PCs to use - 15
PCElbowPlot(all_treatments)

# Build SSN for cluster identification
all_treatments <- BuildSNN(all_treatments, dims.use  = 1:20)

# Find clusters bases on PCs
all_treatments <- FindClusters(all_treatments, resolution = .03, print.output = 0, reuse.SNN = T)

# Run TSNE on PC data
all_treatments <- RunTSNE(all_treatments, dims.use = 1:20, do.fast = T, dim_embed = 2)


# Find markers for each cluster
cluster_marks = FindAllMarkers(all_treatments)

# Cluster 3 appears to be erythroid precursors. Remove these cells and also filter other low quality cells.

# Find cells with high absolute number of mitochondrial and HBB genes
mito_genes = grep(pattern = "^mt-", x = rownames(x = all_treatments@data), value = TRUE)
HBB_genes = grep(pattern = "^Hbb", x = rownames(x = all_treatments@data), value = TRUE)
percent_mito = colSums(all_treatments@raw.data[mito_genes, ]) / colSums(all_treatments@raw.data)
percent_HBB = colSums(all_treatments@raw.data[HBB_genes, ]) / colSums(all_treatments@raw.data)

all_treatments@meta.data$percent_HBB = percent_HBB
all_treatments@meta.data$percent_mito = percent_mito

VlnPlot(all_treatments, c("percent_mito", "percent_HBB"))

# Filter cells with HBB higher than 15% of reads from mitochondiral and/or 10% of reads from HBB 
all_treatments = FilterCells(object = all_treatments, subset.names = c("percent_HBB", "percent_mito"), low.thresholds = c(-Inf, -Inf), high.thresholds = c(0.1, 0.15))

# Remove remaning cells from 
tmp = row.names(all_treatments@meta.data[all_treatments@meta.data$res.0.03 != 3,])

all_treatments <- SubsetData(all_treatments, cells.use = tmp)

# Rerun analysis

all_treatments <- FindVariableGenes(all_treatments, x.low.cutoff = 0.0125, y.cutoff = 1.5, do.contour = F, do.plot = T)

# Compute PCs from the data - NOTE must run this way because the object is too large to scale data on my computer
all_treatments <- ScaleData(all_treatments, genes.use = all_treatments@var.genes, vars.to.regress = c("percent_mito","nUMI"))

all_treatments <- RunPCA(all_treatments, pc.genes = all_treatments@var.genes)

# Find optimal # of PCs to use - 15
PCElbowPlot(all_treatments)

# Build SSN for cluster identification
all_treatments <- BuildSNN(all_treatments, dims.use  = 1:20, force.recalc = T)

# Remove old clusterings from metadata
all_treatments@meta.data <- all_treatments@meta.data[,-5]

# Find clusters bases on PCs
all_treatments <- FindClusters(all_treatments, resolution = .02, print.output = 0, reuse.SNN = T)
all_treatments <- FindClusters(all_treatments, resolution = .05, print.output = 0, reuse.SNN = T)

# Run TSNE on PC data
all_treatments <- RunTSNE(all_treatments, dims.use = 1:20, do.fast = T, dim_embed = 2)

# Add mayo sample IDs to the object
current.cluster.ids = row.names(metadata)
new.cluster.ids = metadata$Mayo.ID
all_treatments@meta.data$Mayo_ID = plyr::mapvalues(x = all_treatments@meta.data$orig.ident, from = current.cluster.ids, to = new.cluster.ids)

TSNEPlot(all_treatments)

saveRDS(all_treatments, file = "Data/Mayo_Mouse/all_treatments_new.rds")


# MM cells
FeaturePlot(all_treatments, c("Sdc1"))

# T cells
FeaturePlot(all_treatments,  c("Il7r", "Cd3g", "Cd3d", "Cd3e", "Cd4", "Cd8a"), pt.size = .3)

# B cells
FeaturePlot(all_treatments,  "Ms4a1")

# NK cells
FeaturePlot(all_treatments,  "Nkg7")

# CD14 Monocytes
FeaturePlot(all_treatments,  c("Cd14","Lyz2"))

#FCGR3A Monocytes
FeaturePlot(all_treatments,  c("Fcgr3", "Ms4a7"))

# Dendritic cells
FeaturePlot(all_treatments,  c("Cst3", "Fcer1a"))


tmp = all_treatments@meta.data$res.0.05

tmp[tmp == 0 | tmp == 7 | tmp == 8 | tmp == 9 | tmp == 10] <- "Myeloma cells"
tmp[tmp == 1] <- "T cells"
tmp[tmp == 6] <- "NK cells"
tmp[tmp == 2] <- "B cells"
tmp[tmp == 3 | tmp == 4 | tmp == 5 | tmp == 11]  <- "Myeloid cells"

all_treatments@meta.data$cell_type <- tmp 

# saveRDS(all_treatments, file = "~/Data/Mayo_Mouse/all_treatments_new.rds")
# all_treatments <- readRDS("~/Data/Mayo_Mouse/all_treatments_new.rds")



#### Tumor ####

#tumor = readRDS(file = "Data/Mayo_Mouse/tumor_subset_new.rds")

temp = as.vector( rownames(all_treatments@meta.data[all_treatments@meta.data$cell_type == "Myeloma cells",]))

subs = SubsetData(all_treatments, cells.use = temp)
subs <- FindVariableGenes(subs)

subs <- ScaleData(subs, genes.use = subs@var.genes)

subs <- RunPCA(subs, pc.genes = subs@var.genes)

subs <- BuildSNN(subs, dims.use = 1:15)

subs <- FindClusters(subs, resolution = seq(.1,1.5,.20), print.output = 0, reuse.SNN = T)
subs <- FindClusters(subs, resolution = .02, print.output = 0, reuse.SNN = T)

subs <- RunTSNE(subs, dims.use = 1:15, do.fast = T)

subs <- SetIdent(subs, ident.use = subs@meta.data$treatment)

tumor = subs

saveRDS(file = "Data/Mayo_Mouse/tumor_subset_new.rds", tumor)

tumor <- SetIdent(tumor, ident.use = tumor@meta.data$treatment)

tumor_marks = FindAllMarkers(tumor, only.pos = T, do.print = T)

write.table(tumor_marks, file = "Data/Mayo_Mouse/tumor_marks.txt", quote = F, sep = "\t")

LCL161_upreg = tumor_marks[tumor_marks$cluster == "LCL161",]

# secreted_genes = read.table("Data/protein_class_Predicted.txt", header = T, sep = "\t", quote = "$$")
# 
# simpleCap <- function(x) {
#   x <- tolower(x)
#   s <- strsplit(x, " ")[[1]]
#   paste(toupper(substring(s, 1,1)), substring(s, 2),
#         sep="", collapse=" ")
# }
# 
# secreted_genes$Gene <- apply(as.data.frame(secreted_genes$Gene), 1, simpleCap)
# 


#### B-cell ####
#B_cell = readRDS(file = "Data/Mayo_Mouse/B_cell_subset_new.rds")

temp = as.vector( rownames(all_treatments@meta.data[all_treatments@meta.data$cell_type== "B cells" ,]))


subs = SubsetData(all_treatments, cells.use = temp)
subs <- FindVariableGenes(subs)

subs <- ScaleData(subs, genes.use = subs@var.genes)

subs <- RunPCA(subs, pc.genes = subs@var.genes)

subs <- BuildSNN(subs, dims.use = 1:15)

subs <- FindClusters(subs, resolution = seq(.1,1.5,.20), print.output = 0, reuse.SNN = T)
subs <- FindClusters(subs, resolution = .02, print.output = 0, reuse.SNN = T)

subs <- RunTSNE(subs, dims.use = 1:15, do.fast = T)

saveRDS(subs, file = "Data/Mayo_Mouse/B_cell_subset_new.rds")

B_cell <- subs

B_cell <- SetIdent(B_cell, ident.use = B_cell@meta.data$treatment)

B_cell_marks = FindAllMarkers(B_cell, only.pos = T, do.print = T)

write.table(B_cell_marks, file = "Data/Mayo_Mouse/B_cell_marks.txt", quote = F, sep = "\t")


 
#### T-Cell ####
#T_cell = readRDS(file = "Data/Mayo_Mouse/T_cell_subset_new.rds")
temp = as.vector( rownames(all_treatments@meta.data[all_treatments@meta.data$cell_type == "T cells",]))


subs = SubsetData(all_treatments, cells.use = temp)
subs <- FindVariableGenes(subs)

subs <- ScaleData(subs, genes.use = subs@var.genes)

subs <- RunPCA(subs, pc.genes = subs@var.genes)

subs <- BuildSNN(subs, dims.use = 1:15)

subs <- FindClusters(subs, resolution = seq(.1,1.5,.20), print.output = 0, reuse.SNN = T)
subs <- FindClusters(subs, resolution = .02, print.output = 0, reuse.SNN = T)

subs <- RunTSNE(subs, dims.use = 1:15, do.fast = T)

T_cell = subs

saveRDS(T_cell, file = "Data/Mayo_Mouse/T_cell_subset_new.rds")

T_cell <- SetIdent(T_cell, ident.use = T_cell@meta.data$treatment)

T_cell_marks = FindAllMarkers(T_cell, only.pos = T, do.print = T)

write.table(T_cell_marks, file = "Data/Mayo_Mouse/T_cell_marks.txt", quote = F, sep = "\t")



#### myeloid ####
# myeloid = readRDS(file = "~/Data/Mayo_Mouse/myeloid_subset.rds")

temp = as.vector( rownames(all_treatments@meta.data[all_treatments@meta.data$cell_type == "Myeloid cells" ,]))


subs = SubsetData(all_treatments, cells.use = temp)
subs <- FindVariableGenes(subs)

subs <- ScaleData(subs, genes.use = subs@var.genes)

subs <- RunPCA(subs, pc.genes = subs@var.genes)

subs <- BuildSNN(subs, dims.use = 1:15)

subs <- FindClusters(subs, resolution = seq(.1,1.5,.20), print.output = 0, reuse.SNN = T)
subs <- FindClusters(subs, resolution = .02, print.output = 0, reuse.SNN = T)
subs <- FindClusters(subs, resolution = .05, print.output = 0, reuse.SNN = T)

subs <- RunTSNE(subs, dims.use = 1:15, do.fast = T)

myeloid <- subs

saveRDS(myeloid, file = "Data/Mayo_Mouse/myeloid_subset.rds")

myeloid <- SetIdent(myeloid, ident.use = myeloid@meta.data$treatment)

myeloid_marks = FindAllMarkers(myeloid, only.pos = T, do.print = T)

write.table(myeloid_marks, file = "Data/Mayo_Mouse/myeloid_marks.txt", quote = F, sep = "\t")

myeloid_LCL161_marks = FindMarkers(myeloid, ident.1 = "LCL161", ident.2 = "LCL161_cont", only.pos = T, do.print = T, min.pct = 0)

myeloid_LCL161_marks <- myeloid_LCL161_marks[order(myeloid_LCL161_marks$avg_logFC),]


#### Analysis ####
tmp = as.matrix(myeloid@data[c("Cstb", "Ccl4"),myeloid@meta.data$treatment == "LCL161_cont"])

tmp = t(as.data.frame(tmp))

tmp[tmp == 0] <- NA

tmp = tmp[!(is.na(rowMeans(tmp))),]

cor(tmp)

plot(tmp)


# 
# # LCL specific signals in macrophages
# temp = as.vector( rownames(myeloid@meta.data[myeloid@meta.data$res.0.1 == 0,]))
# macros = subsetData(myeloid, cells.use = temp, do.scale = T, do.center = T)
# macros <- SetIdent(macros, ident.use = macros@meta.data$treatment)
# Macrophage_treatment_marks = FindMarkers(macros, ident.1 = "LCL161", test.use = "negbinom", only.pos = T)
# write.table(file = "Data/Mayo_Mouse/macrophage_LCL161_marks.txt", Macrophage_treatment_marks, sep = "\t", quote = F)
# 
# ids <- bitr(rownames(Macrophage_treatment_marks), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
# ids2 <- bitr(unique(Myeloid_marks$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
# ids3 <- bitr(rownames(Myeloid_marks[Myeloid_marks$cluster == 0,]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
# ids4 <- bitr(rownames(all_treatments@data), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
# 
# x1 <- enrichPathway(gene=ids[,2],pvalueCutoff=0.05, readable=T, organism = "mouse", universe = unique(ids4[,2]))
# x2 <- enrichPathway(gene=ids[,2],pvalueCutoff=0.05, readable=T, organism = "mouse", universe = unique(ids2[,2]))
# 
# 
# temp = as.vector( rownames(macros@meta.data[macros@meta.data$treatment == "LCL161",]))
# macros_lcl = myeloidetData(macros, cells.use = temp, use.imputed = TRUE)
# tmp = prcomp(macros_lcl@data[rownames(macros_lcl@data) %in% macros@var.genes[moduleColors == "turquoise"], ], scale = T, center = T)
# macros_lcl@pca.x <- as.data.frame(tmp$x)
# macros_lcl@pca.rot <- as.data.frame(tmp$rot)
# macros_lcl@pca.obj <- list(tmp)
# macros_lcl <- BuildSNN(macros_lcl, pc.use = 1:10, do.sparse = T)
# 
# macros_lcl <- FindClusters(macros_lcl, resolution = seq(.1,1.5,.20), print.output = 0, reuse.SNN = T)
# macros_lcl <- FindClusters(macros_lcl, resolution = .02, print.output = 0, reuse.SNN = T)
# 
# # Run TSNE on PC data
# macros_lcl <- RunTSNE(macros_lcl, dims.use = 1:10, do.fast = T)
# 
# macros_lcl <- SetIdent(macros_lcl, ident.use = macros_lcl@meta.data$res.0.1)
# 
# 

#### WGCNA ####

# Make data frame for the analysis

datExpr = as.data.frame(t(as.matrix(myeloid@data)))
#run this to check if there are gene outliers
gsg=goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

datExpr <- datExpr[,gsg$goodGenes]

## Check the scale free toplogy 
powers = c(c(1:13), seq(from = 14, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## NOTE: Pick the soft power threshold that maximise the SFT index while minimizing the soft power - want to be > 0.8

#softPower =  min(sft$fitIndices[sft$fitIndices[,2] > .89 ,1])
softPower = 2

# Run blockwise WGCNA on using the standard parameters

net = blockwiseModules(datExpr, power = softPower, detectCutHeight = .9998,
                       TOMType = "unsigned", minModuleSize = 5,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "~/Data/Mayo_Mouse/WGCNA.test",
                       verbose = 3, maxBlockSize = 20000, nThreads = 6)


# net = blockwiseModules(datExpr, power = 3, detectCutHeight = .995,
#                        TOMType = "unsigned", minModuleSize = 5,
#                        reassignThreshold = 1e-6, mergeCutHeight = 0.10,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE, impute = T,
#                        saveTOM =  FALSE,
#                        saveTOMFileBase = "~/Data/Mayo_Mouse/WGCNA.test",
#                        verbose = 3, maxBlockSize = 20000)


### Plot the dendrogram and the module colors underneath ###

net = readRDS("~/Data/Mayo_Mouse/WGCNA_net_object.rds")

mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



# Get module information and relable
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors, excludeGrey = T)$eigengenes
MEs = orderMEs(MEs0)

modNames = substring(names(MEs), 3)

# Plot Eigengenes

ME2 = MEs[, colnames(MEs) %in% paste("ME",names(table(moduleColors)[table(moduleColors) > 20]), sep = "")]

plotEigengeneNetworks(ME2, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

plotEigengeneNetworks(ME2, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)



TOM = TOMsimilarityFromExpr(datExpr, power = 2, TOMType = "signed")
# Select modules
modules = c("salmon")
# Select module probes
genes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modGenes = genes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modGenes, modGenes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("~/Data/Mayo_Mouse/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               #nodeFile = paste("~/Data/Mayo_Mouse/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.01,
                               nodeNames = modGenes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);



geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste(modNames);
names(MMPvalue) = paste("p.", modNames, sep="")

syms = rownames(geneModuleMembership)
ids <- bitr(syms, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

tmp = cbind(syms, moduleColors)

tmp2 = merge(ids, tmp)

rowMeans(as.matrix(myeloid@data[names(moduleColors[moduleColors == "blue"]), myeloid@meta.data$treatment == "LCL161_cont" ]))


## WGCNA LCL161 specific ##


datExpr = as.data.frame(t(as.matrix(myeloid@data)))
datExpr = datExpr[row.names(myeloid@meta.data[myeloid@meta.data$treatment == "LCL161",]),]

#run this to check if there are gene outliers
gsg=goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

datExpr <- datExpr[,gsg$goodGenes]

## Check the scale free toplogy 
powers = c(c(1:13), seq(from = 14, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## NOTE: Pick the soft power threshold that maximise the SFT index while minimizing the soft power - want to be > 0.8

#softPower =  min(sft$fitIndices[sft$fitIndices[,2] > .89 ,1])
softPower = 2

# Run blockwise WGCNA on using the standard parameters

net_LCL161 = blockwiseModules(datExpr, power = softPower, detectCutHeight = .98,
                       TOMType = "unsigned", minModuleSize = 5,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "~/Data/Mayo_Mouse/WGCNA.test",
                       verbose = 3, maxBlockSize = 20000, nThreads = 6)

#net = readRDS("~/Data/Mayo_Mouse/WGCNA_net_object.rds")

mergedColors = labels2colors(net_LCL161$colors)
plotDendroAndColors(net_LCL161$dendrograms[[1]], mergedColors[net_LCL161$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



# Get module information and relable
moduleLabels = net_LCL161$colors
moduleColors = labels2colors(net_LCL161$colors)
names(moduleColors) <- colnames(datExpr)
MEs = net_LCL161$MEs;
genet_LCL161ree = net_LCL161$dendrograms[[1]]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors, excludeGrey = T)$eigengenes
MEs = orderMEs(MEs0)

modNames = substring(names(MEs), 3)

# Plot Eigengenes

ME2 = MEs[, colnames(MEs) %in% paste("ME",names(table(moduleColors)[table(moduleColors) > 20]), sep = "")]

plotEigengeneNetworks(ME2, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

plotEigengeneNetworks(ME2, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


table(moduleColors[row.names(tail(myeloid_LCL161_marks[myeloid_LCL161_marks$pct.2 <= .1,], 50))])

TOM = TOMsimilarityFromExpr(datExpr, power = 2, TOMType = "signed")
# Select modules
modules = c("brown")
# Select module probes
genes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modGenes = genes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modGenes, modGenes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("~/Data/Mayo_Mouse/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               #nodeFile = paste("~/Data/Mayo_Mouse/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.12,
                               nodeNames = modGenes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);




tmp = datExpr[,names(moduleColors[moduleColors == "tan"])]

tmp2 = apply(tmp, 2, function(xx){length(xx[xx > 0])})

tmp = tmp[,order(tmp2)]

tmp = tmp[order(rowSums(tmp)),]

heatmap.2(as.matrix(tmp),dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')


## WGCNA LCL161_cont ##


datExpr = as.data.frame(t(as.matrix(myeloid@data)))
datExpr = datExpr[row.names(myeloid@meta.data[myeloid@meta.data$treatment == "LCL161_cont",]),]

#run this to check if there are gene outliers
gsg=goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

datExpr <- datExpr[,gsg$goodGenes]

## Check the scale free toplogy 
powers = c(c(1:13), seq(from = 14, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## NOTE: Pick the soft power threshold that maximise the SFT index while minimizing the soft power - want to be > 0.8

#softPower =  min(sft$fitIndices[sft$fitIndices[,2] > .89 ,1])
softPower = 2

# Run blockwise WGCNA on using the standard parameters

net_LCL161_cont = blockwiseModules(datExpr, power = softPower, detectCutHeight = .98,
                              TOMType = "unsigned", minModuleSize = 5,
                              reassignThreshold = 0, mergeCutHeight = 0.25,
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              saveTOMs = TRUE,
                              saveTOMFileBase = "~/Data/Mayo_Mouse/WGCNA.test",
                              verbose = 3, maxBlockSize = 20000, nThreads = 6)

#net = readRDS("~/Data/Mayo_Mouse/WGCNA_net_object.rds")

mergedColors = labels2colors(net_LCL161_cont$colors)
plotDendroAndColors(net_LCL161_cont$dendrograms[[1]], mergedColors[net_LCL161_cont$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



# Get module information and relable
moduleLabels = net_LCL161_cont$colors
moduleColors = labels2colors(net_LCL161_cont$colors)
names(moduleColors) <- colnames(datExpr)
MEs = net_LCL161_cont$MEs;
genet_LCL161_contree = net_LCL161_cont$dendrograms[[1]]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors, excludeGrey = T)$eigengenes
MEs = orderMEs(MEs0)

modNames = substring(names(MEs), 3)

# Plot Eigengenes

ME2 = MEs[, colnames(MEs) %in% paste("ME",names(table(moduleColors)[table(moduleColors) > 20]), sep = "")]

plotEigengeneNetworks(ME2, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

plotEigengeneNetworks(ME2, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)




TOM = TOMsimilarityFromExpr(datExpr, power = 2, TOMType = "signed")
# Select modules
modules = c("brown")
# Select module probes
genes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modGenes = genes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modGenes, modGenes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("~/Data/Mayo_Mouse/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               #nodeFile = paste("~/Data/Mayo_Mouse/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.12,
                               nodeNames = modGenes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);






#### TopGO ####

# DE_genes = FindMarkers(subs, ident.1 = 8, only.pos = T)
# 
# onion = subs@data
# onion = as.matrix(onion)
# 
# tmp = apply(onion, 1, function(xx){length(xx[xx >0])/length(xx)})
# 
# tmp2 = tmp[tmp > .01]
# 
# tmp2[names(tmp2) %in% rownames(DE_genes)] <-1
# tmp2[!(names(tmp2) %in% rownames(DE_genes))] <-0



get.sig.genes = function(x){x == 1}
genes = (rep(0, length(moduleColors)))
names(genes) <- names(moduleColors)
genes[moduleColors == "tan"] <-1


GOdata <- new("topGOdata",
              description = "Enrichment", ontology = "BP",
              allGenes = genes, 
              geneSel = get.sig.genes,
              nodeSize = 10,
              annot = annFUN.org, mapping="org.Mm.eg.db", ID="Symbol")

allGO = genesInTerm(object = GOdata) 

resultfisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

GO.re= GenTable(GOdata, classicfisher = resultfisher, topNodes = length(allGO))

SAM_ANOTATION = lapply(allGO,function(x) x[x %in% names(moduleColors[moduleColors == "salmon"])] )

SAM_ANOTATION["GO:1903555"]

