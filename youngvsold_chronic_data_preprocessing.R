####################################
# sc.mm.gam.gender                 #
# data preprocessing and analysis  #
####################################

# Load packages
library(parallel)
library(future)
library(Seurat)
library(Matrix)
options(future.globals.maxSize=256*1024^3)
plan(multiprocess)
library(ggplot2)
library(cowplot)
library(magrittr)
library(dplyr)
library(tidyverse)
library(purrr)

# Analysis parameters
anchor_dims <- 30            # number of anchor dimensions used for biological replicates integration
pca_dims <- 30               # number of PCA dimensions to compute and use in tSNE, UMAP                             
umap_n_neighbors <- 30       # UMAP parameter,                                                                    
clustering_resolution <- 0.6 # Resolution parameter for Seurat clustering
n_features <- 2000

# Read raw data (gene/cell count matrix from cellranger, filtered: use only detected cellular barcodes)
Aged_TBI.data <- Read10X(data.dir = "/Users/zhangyingchen/Desktop/SCR Chronic Aged v Young/samples/Aged/Aged_TBI/filtered_feature_bc_matrix/")
Aged_TBI <- CreateSeuratObject(counts = Aged_TBI.data, project = "Aged_TBI")

Aged_sham.data <- Read10X(data.dir = "/Users/zhangyingchen/Desktop/SCR Chronic Aged v Young/samples/Aged/Aged_sham/filtered_feature_bc_matrix/")
Aged_sham <- CreateSeuratObject(counts = Aged_sham.data, project = "Aged_sham", min.cells = 3, min.features = 200)

Young_TBI.data <- Read10X(data.dir = "/Users/zhangyingchen/Desktop/SCR Chronic Aged v Young/samples/Young/Young_TBI/filtered_feature_bc_matrix/")
Young_TBI <- CreateSeuratObject(counts = Young_TBI.data, project = "Young_TBI", min.cells = 3, min.features = 200)

Young_sham.data <- Read10X(data.dir = "/Users/zhangyingchen/Desktop/SCR Chronic Aged v Young/samples/Young/Young_sham/filtered_feature_bc_matrix/")
Young_sham <- CreateSeuratObject(counts = Young_sham.data, project = "Young_sham", min.cells = 3, min.features = 200)

# Filter cells (based on percent_mito and nFeature_RNA), normalize (log-normalize in columns) 
# and find n_features most variable features/genes (using Variance-stabilizing transformation)
# to the list of top variable genes add markers of microglia, macrophages and cell cycle                                              nfeatures = n_features, verbose = FALSE)

Young_sham$mitoRatio <- PercentageFeatureSet(object = Young_sham, pattern = "^mt-")
Young_sham$mitoRatio <- Young_sham@meta.data$mitoRatio / 100
Young_TBI$mitoRatio <- PercentageFeatureSet(object = Young_TBI, pattern = "^mt-")
Young_TBI$mitoRatio <- Young_TBI@meta.data$mitoRatio / 100
Aged_sham$mitoRatio <- PercentageFeatureSet(object = Aged_sham, pattern = "^mt-")
Aged_sham$mitoRatio <- Aged_sham@meta.data$mitoRatio / 100
Aged_TBI$mitoRatio <- PercentageFeatureSet(object = Aged_TBI, pattern = "^mt-")
Aged_TBI$mitoRatio <- Aged_TBI@meta.data$mitoRatio / 100

Young_sham <- subset(Young_sham,nCount_RNA < 10000 & nFeature_RNA < 5000 & mitoRatio < 0.2)
Young_TBI <- subset(Young_TBI, nCount_RNA < 10000  & nFeature_RNA < 5000 & mitoRatio < 0.2)
Aged_sham <- subset(Aged_sham,nCount_RNA < 10000 & nFeature_RNA < 5000 & mitoRatio < 0.2)
Aged_TBI <- subset(Aged_TBI, nCount_RNA < 10000 & nFeature_RNA < 5000 & mitoRatio < 0.2)


#since two objects are from the same batch, merge them
young_merge <- merge(Young_sham, Young_TBI, 
                      add.cell.ids = c("Young_sham","Young_TBI"), project = "young_merge")
young_merge <- FindVariableFeatures(object = young_merge,
                                    selection.method = "vst", 
                                    nfeatures = 3000, verbose = FALSE)
young_merge <- SCTransform(young_merge, vars.to.regress = c("mitoRatio"))
aged_merge <- merge(Aged_sham, Aged_TBI, 
                     add.cell.ids = c("Aged_sham","Aged_TBI"), project = "aged_merge")
aged_merge <- FindVariableFeatures(object = aged_merge,
                                    selection.method = "vst", 
                                    nfeatures = 3000, verbose = FALSE)
aged_merge <- SCTransform(aged_merge, vars.to.regress = c("mitoRatio"))

# Add meta data
#samples_objects$Young_sham$condition <- "sham"
#samples_objects$Young_TBI$condition <- "TBI"   

# Add selected (previously reported) genes to var.features

#samples_objects <- mclapply(seq_along(samples_objects), function(i) {
#  samples_objects[[i]]@assays$RNA@var.features <- unique(c(samples_objects[[i]]@assays$RNA@var.features, 
#                                                           s_genes, g2m_genes, microglia_markers, macrophages_markers))
#  samples_objects[[i]]})

#names(samples_objects) <- names(samples_raw_data)
#young_merge@assays$RNA@var.features <- unique(c(young_merge@assays$RNA@var.features, 
#                                    s_genes, g2m_genes, microglia_markers, macrophages_markers))
#aged_merge@assays$RNA@var.features <- unique(c(aged_merge@assays$RNA@var.features, 
#                                               s_genes, g2m_genes, microglia_markers, macrophages_markers))

# Integrate replicates within conditions, scale, regress out unwanted sources of variation, 
# calculate PCA, t-SNE and find cell clusters
#integrate for datasets normalized with the sctransform worflow
features <- SelectIntegrationFeatures(object.list = c(young_merge, aged_merge), nfeatures = 3000)
merge.list <- PrepSCTIntegration(object.list = c(young_merge, aged_merge), anchor.features = features)
samples_anchors <- FindIntegrationAnchors(object.list = merge.list, 
                                          dims = 1:anchor_dims,
                                          anchor.features = features,
                                          normalization.method="SCT")
samples_integrated <- IntegrateData(anchorset = samples_anchors, dims = 1:30, normalization.method = "SCT")
DefaultAssay(object=samples_integrated) <- "integrated"


#samples_integrated <- CellCycleScoring(object = samples_integrated, s.features = s_genes, 
#                                       g2m.features = g2m_genes, set.ident = TRUE)
#samples_integrated$CC_Difference <- samples_integrated$S.Score - samples_integrated$G2M.Score
#samples_integrated <- ScaleData(object = samples_integrated, verbose = FALSE, 
#                                vars.to.regress = c("nCount_RNA", 
#                                                    "percent_mito", 
#                                                    "CC_Difference"))

samples_integrated <- RunPCA(object = samples_integrated, verbose = FALSE)
samples_integrated <- RunUMAP(samples_integrated, dims = 1:30, verbose = FALSE)
samples_integrated <- FindNeighbors(object = samples_integrated, dims = 1:pca_dims)
samples_integrated <- FindClusters(object = samples_integrated, resolution = clustering_resolution) 
samples_integrated
DimPlot(samples_integrated, split.by = "orig.ident", ncol = 2,label = TRUE)


#Main figures
#microglia cluster 0,1,3,14,17
FeaturePlot(object = samples_integrated, features =  c("P2ry12","Tmem119","Hexb","C1qa","C1qb","Ctss"), order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)

#T cell cluster 2,6,8
T_cell <- c("Cd8a","Il7r","Cd3d","Cd4")
FeaturePlot(object = samples_integrated, features = T_cell, order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)

#CD4+ T cell cluster 
FeaturePlot(object = samples_integrated, features = c("Cd4","Ccr7"), sort.cell = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE, blend = TRUE)

#B cell cluster 5
B_cell <- c("Cd79a","Cd79b","Cd74","Ly6d","Ms4a1")
FeaturePlot(object = samples_integrated, features = B_cell, order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)

#NK cluster cluster 12
NK_cell <- c("Ncr1","Klre1") 
FeaturePlot(object = samples_integrated, features = c("Ncr1","Klre1"), order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)

#to find subclusters
#samples_integrated <-FindSubCluster(
#  samples_integrated,
#  8,
#  "integrated_snn",
#  subcluster.name = "NK",
#  resolution = 0.5,
#  algorithm = 1
#)
#DimPlot(samples_integrated, group.by = "NK", label = TRUE, label.size = 6)
#Idents(samples_integrated)=samples_integrated$NK

Idents(samples_integrated)=samples_integrated$seurat_clusters
samples_integrated_subset <- subset(samples_integrated, idents = c((0:3),(5:6),(8:9),(12:18)))
DimPlot(samples_integrated_subset,label = TRUE, label.size = 6)

#CAM cluster part of 9_1
Idents(samples_integrated_subset)=samples_integrated_subset$seurat_clusters

CAM <- c("Lyve1","Cd163","Siglec1")
FeaturePlot(object = samples_integrated_subset, features = c("Lyve1","Cd163","Siglec1"), order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)
VlnPlot(samples_integrated_subset,features = c("Lyve1","Cd163"))

#to find subclusters
samples_integrated_subset <-FindSubCluster(
  samples_integrated_subset,
  9,
  graph.name = "seurat_clusters",
  subcluster.name = "sub2",
  resolution = 0.1,
  algorithm = 1
)

DimPlot(samples_integrated_subset, group.by = "sub", label = TRUE)

#classical monocytes/macrophages cluster 
Mo_MΦ <- c("Ly6c1","Ly6c2", "Ccr2","Itga4","Tgfbi","Ifitm2","Ifitm3","S100a6")
FeaturePlot(object = samples_integrated, features = Mo_MΦ,  order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)
VlnPlot(Aged_M_TBI_chronic_SCT,features = c("Ncr1","Klre1"))

#UK cluster
UK <- c("Meg3","Clu")
FeaturePlot(object = samples_integrated, features = annot[match(UK, annot$Gene.name), "Gene.stable.ID"],  order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)

#neutrophil cluster 13 
NP <- c("Cxcr2", "Csf3r","Ly6g","S100a9","S100a8")
FeaturePlot(object = samples_integrated, features = c("Cxcr2", "Csf3r","Ly6g","S100a9","S100a8"), order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)

#Nonclassical monocytes cluster
NonMo <- c("Spn", "Ace")
FeaturePlot(object = samples_integrated, features = c("Spn", "Ace"),  order = TRUE, reduction = "umap", ncol=1, min.cutoff = 'q10', label = TRUE)
VlnPlot(samples_integrated,features = c("Spn", "Ace"))

DefaultAssay(samples_integrated) <- "RNA"

# Create function to get conserved markers in all conditions -useful when you are unsure of the identify for a cluster
#get_conserved <- function(cluster){
#  FindConservedMarkers(samples_integrated,
#                       ident.1 = cluster,
#                       grouping.var = "orig.ident",
#                       only.pos = TRUE) %>%
#    rownames_to_column(var = "gene") %>%
#    cbind(cluster_id = cluster, .)
#}
#conserved_markers <- map_dfr(c(0,18), get_conserved)

cluster.4.markers <- FindConservedMarkers(samples_integrated, ident.1 = 4, grouping.var = "orig.ident", verbose = FALSE, only.pos = TRUE)
cluster.7.markers <- FindConservedMarkers(samples_integrated, ident.1 = 7, grouping.var = "orig.ident", verbose = FALSE, only.pos = TRUE)
cluster.10.markers <- FindConservedMarkers(samples_integrated, ident.1 = 10, grouping.var = "orig.ident", verbose = FALSE, only.pos = TRUE)
cluster.11.markers <- FindConservedMarkers(samples_integrated, ident.1 = 11, grouping.var = "orig.ident", verbose = FALSE, only.pos = TRUE)
cluster.15.markers <- FindConservedMarkers(samples_integrated, ident.1 = 15, grouping.var = "orig.ident", verbose = FALSE, only.pos = TRUE)
cluster.16.markers <- FindConservedMarkers(samples_integrated, ident.1 = 16, grouping.var = "orig.ident", verbose = FALSE, only.pos = TRUE)
cluster.18.markers <- FindConservedMarkers(samples_integrated, ident.1 = 18, grouping.var = "orig.ident", verbose = FALSE, only.pos = TRUE)


#find all markers for each cluster
all_markers <-FindAllMarkers(samples_integrated, only.pos = TRUE, min.pct =  0.25,min.diff.pct = 0.25)
allcells_referece <-FindAllMarkers(samples_integrated, only.pos = TRUE, min.pct =  0.25,min.diff.pct = 0.25, logfc.threshold = 0)
write.csv(allcells_referece, 
          file = "/Users/zhangyingchen/Desktop/allcells_referece.csv", 
          quote = FALSE, 
          row.names = FALSE)

top15_comb <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 15, avg_log2FC)
# Save markers to file
write.csv(top15_comb, 
          file = "/Users/zhangyingchen/Desktop/top15_comb.csv", 
          quote = FALSE, 
          row.names = FALSE)

# Figure 

#"Spn", "Ace""Cxcr2", "Csf3r","Ly6g","S100a9","S100a8""Lyve1","Cd163","Siglec1""Ncr1","Klre1""Cd79a","Cd79b","Cd74","Ly6d","Ms4a1""Cd8a","Il7r","Cd3d","Cd4"
gene_panel <- c("Tmem119", "P2ry12","Ctss","C1qa","Pf4","Tgfbi", "Ifitm2", "Ifitm3", "Ly6c2",
                      "Ccr2", "Mrc1", "Cd163", "Cd24a", "Ncam1", "Klrk1", "Ncr1", "Cd2", "Cd3d",
                      "Cd4", "Cd8a", "Spn","Ace","Cxcr2","Lyve1","Cd79a")


DotPlot(samples_integrated_subset, features = gene_panel) + 
  RotatedAxis() +
  theme(legend.position="top")

#neutrophil cluster 13 
#CAM&MΦ cluster part of 9
#NK cluster cluster 12
#B cell cluster 5
#T cell cluster 2,6,8
#microglia cluster 0,1,3,14,17
#Mo/MΦ cluster 16, 18
#NCMo cluster 15 (comprising some nonclassical monocytes labeled by Spn, Ace)
#Y.UKs 4, 7
#others 10, 11

samples_integrated_subset <- subset(samples_integrated, idents = c((0:9),(12:18)))
VlnPlot(samples_integrated_subset,features = c("Pf4", "F13a1","Lyz2","Ifit3","Dab2","Cxcr2","Scd2","Csf1","Crybb1","Fcrls","Selplg","Mafb","Pmepa1","Cd14"), pt.size = 0, assay = "RNA")

samples_integrated_subset$cell_type_8_groups <- plyr::mapvalues(Idents(samples_integrated_subset), 
                                                        from=c(0:9, 12:18), 
                                                         to=c("MG", "MG", "T", "MG", "Y.UK", "B","T",
                                                              "Y.UK","T","CAM&MΦ","NK","NP",
                                                              "MG", "NCMo", "MoMΦ", "MG", "MoMΦ"))

Idents(samples_integrated_subset)=samples_integrated_subset$cell_type_8_groups
DimPlot(samples_integrated_subset, split.by = "orig.ident", ncol = 2,cols = c("#0DD1AD", "#C2B4FC","grey", "#EC5CA5","#8c42a3","#2DA7C8","#E98934","#FABF00","#7f9c00"))
#0DD1AD,C2B4FC, 2DA7C8,grey,EC5CA5,8c42a3,E98934,FABF00,7f9c00
pt <- table(Idents(samples_integrated_subset), samples_integrated_subset$orig.ident)
pt <- as.data.frame(pt)
pt_aged_sham = pt[c(1:8),]
pt_aged_TBI = pt[c(9:16),]
pt_young_sham = pt[c(17:24),]
pt_young_TBI = pt[c(25:32),]

pt_aged_sham$Var1 <- as.character(pt_aged_sham$Var1)
pt_aged_TBI$Var1 <- as.character(pt_aged_TBI$Var1)
pt_young_sham$Var1 <- as.character(pt_young_sham$Var1)
pt_young_TBI$Var1 <- as.character(pt_young_TBI$Var1)

library(ggplot2)
library(ggrepel)
library(ggpubr)
MG <- "#EC5CA5"
T <- "#8c42a3"
B <- "#EC5CA5"
CAMMPC <- "#2DA7C8"
MoMΦ<-"#FABF00"
NP<- "#E98934"
NK<-"#C2B4FC" 
UK<-"grey"
NCM<-"#7f9c00"

ggplot(pt_aged_sham, aes(x=" ", y=Freq, fill=Var1))+
  geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("#EC5CA5", "#8c42a3","#0DD1AD", "#FABF00","#2DA7C8","#E98934","#C2B4FC","grey"))+
  theme_light()+
  geom_label_repel(aes(label = scales::percent((Freq)/sum(Freq),0.5)), size=4, show.legend = F, nudge_x = 1)

ggplot(pt_aged_TBI, aes(x=" ", y=Freq, fill=Var1))+
  geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("#EC5CA5", "#8c42a3","#0DD1AD", "#FABF00","#2DA7C8","#E98934","#C2B4FC","grey"))+
  theme_light()+
  geom_label_repel(aes(label = scales::percent((Freq)/sum(Freq),0.5)), size=4, show.legend = F, nudge_x = 1)

ggplot(pt_young_sham, aes(x=" ", y=Freq, fill=Var1))+
  geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("#EC5CA5", "#8c42a3","#0DD1AD", "#FABF00","#2DA7C8","#E98934","#C2B4FC","grey"))+
  theme_light()+
  geom_label_repel(aes(label = scales::percent((Freq)/sum(Freq),0.5)), size=4, show.legend = F, nudge_x = 1)

ggplot(pt_young_TBI, aes(x=" ", y=Freq, fill=Var1))+
  geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("#EC5CA5", "#8c42a3","#0DD1AD", "#FABF00","#2DA7C8","#E98934","#C2B4FC","grey"))+
  theme_light()+
  geom_label_repel(aes(label = scales::percent((Freq)/sum(Freq),0.5)), size=4, show.legend = F, nudge_x = 1)

ggarrange(p1, p2, p3, p4, ncol=2)
