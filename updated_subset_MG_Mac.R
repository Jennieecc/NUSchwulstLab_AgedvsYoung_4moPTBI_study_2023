# Subset microglia and macrophage clusters
MG.cells <- subset(samples_integrated_subset,idents = c("MG"))

#NormalizeData,Scale,FindVariableFeatures replaced by sctransform, so go to run PCR, UMAP, Find Neighbors, Find Clusters
MG.cells <- SCTransform(MG.cells, vars.to.regress = c("mitoRatio"), verbose = FALSE)
#MG.cells <- FindVariableFeatures(object = MG.cells,
#                                selection.method = "vst", 
#                                nfeatures = 3000, verbose = FALSE)
MG.cells <- RunPCA(object = MG.cells, verbose = FALSE, assay = "SCT")
MG.cells <- RunUMAP(MG.cells, dims = 1:30, verbose = FALSE)
MG.cells <- FindNeighbors(object = MG.cells, dims = 1:30)
MG.cells_1 <- FindClusters(object = MG.cells,resolution = 1.0)

DimPlot(MG.cells_1, label = TRUE)
MG.cells_1<-subset(MG.cells_1, idents = c(0:13))

DefaultAssay(object=MG.cells_1) <- "RNA"
marker_MG.cells_1 <- FindAllMarkers(object = MG.cells_1,logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)
marker_MG.cells_1 <- marker_MG.cells_1[!grepl("Rp", rownames(marker_MG.cells_1)),]
marker_MG.cells_1 <- marker_MG.cells_1[!grepl("Mt", rownames(marker_MG.cells_1)),]
#reference
reference_MG.cells_1 <- FindAllMarkers(object = MG.cells_1,logfc.threshold = 0, only.pos = TRUE)


top15_MG <- marker_MG.cells_1 %>%
  group_by(cluster) %>%
  top_n(n = 15, avg_log2FC)
write.csv(reference_MG.cells_1, file="/Users/zhangyingchen/Desktop/reference_MG.cells_1.csv")

#annotation
#microglia subsets based on functionality
VlnPlot(subset_Mic,features=c("P2ry12","Tmem119","Lpl","Cst7","Apoe","Ccl3","Ccl4","Cst7","Lpl","Mafb"),stack = T,flip = T)

FeaturePlot(object = subset_Mic, features = c("Cd68","Ccl3","Ccl4","Cst7","Ifitm3","Rtp4","Spp1","Cd74","H2-Aa",
                                              "Ly86","Ccl2","Cxcl10","Mki67","Ccl5"), order = TRUE,reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
#Homeostatic less in 4 
FeaturePlot(MG.cells, features = c("P2ry12","Tmem119","Hexb"), order = TRUE,reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)

#Dam cluster 4
FeaturePlot(object = MG.cells, features = c("Lpl","Cst7","Apoe","Lyz2"), order = TRUE,reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)

#Inflammatory cluster 0, part of 4 and 6
FeaturePlot(object = MG.cells, features = c("Ccl3","Ccl4","Cst7"), order = TRUE,reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)

#T-cell recruitment 
FeaturePlot(object = MG.cells, features = c("Cxcl10","Ccl5"), order = TRUE,reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)


VlnPlot(MG.cells,features=c("Cd74","Ccl5","Cxcl10","Cd68","Ifitm3","Rtp4","Spp1","Cd74","H2-Aa","Ly86","Ccl2","Ndufb7","Mif"), stack = T, flip = T)

#microglia
FeaturePlot(object = MG.cells, features = c("Tmem119","Cx3cr1","P2ry12","P2ry13","Gpr34","Olfml3","Selplg","Sparc","Fcrls","Siglech","Slc2a5"), order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)
#microglia progenitors and early microglia
FeaturePlot(object = MG.cells, features = c("Pf4","F13a1","Lyz2","Ifit3","Mcm5","Dab2"),order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
#pre-microglia
FeaturePlot(object = MG.cells, features = c("Cxcr2","Scd2","Psat1","Csf1","Crybb1","Fcrls"), order = TRUE,reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
#adult microglia
FeaturePlot(object = MG.cells, features = c("Selplg","Mafb","Pmepa1","Cd14"), order = TRUE,reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
#Dam cluster 4
FeaturePlot(object = MG.cells_1, features = c("Lpl","Cst7","Apoe"), order = TRUE,reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)

FeaturePlot(object = subset_MG_MAC, features = c("Itga4","Tgfbi","Ifitm2","Ifitm3","Tagln2","F13a1",
                                                 "Fpr3","Kynu","S100a11","S100a6"), order = TRUE,reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)

FeaturePlot(object = subset_MG_MAC, features = c("Cd163","Mrc1","Lyve1","Siglec1"), order = TRUE,reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
FeaturePlot(object = subset_MG_MAC, features = c("Ly6c1","Ly6c2","Ccr2",
                                                 "Tgfbi","Spn"), order = TRUE,reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
#aging associated - more activation - Ifi204, Lilrb4, Arhgap, Oas1a, Cd244, Ildr2
FeaturePlot(object = MG.cells, features = c("Ifi204","Lilrb4","Arhgap","Oas1a","Cd244","Ildr2","Cd206","Cd36"), order = TRUE,reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)

pt <- table(MG.cells_1$orig.ident,MG.cells_1$MG_5_groups)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.8) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values=refCols)


#DAM cluster 4
#AgedMG cluster 1  
#YoungMG clustter 5, 6, 
#M1 cluster 0 (enriched by young) M2 cluster 2 (enriched by young) M3 cluster 3 (enriched by aged) 
#M4 cluster 7 (equal), M5 cluster 8, M 6 cluster 9

#rename
MG.cells_1$MG_5_groups <- plyr::mapvalues(Idents(MG.cells_1), 
                                                    from=c(0:13), 
                                                    to = c("A.MG2",  "MG1",  "Y.MG1",  "Y.MG2",  
                                                           "Y.MG1", "MG1","A.MG1","MG1","MG1","A.MG2",
                                                           "MG1","A.MG1","A.MG2","MG1"))
Idents(MG.cells_1)=MG.cells_1$MG_5_groups
DimPlot(MG.cells_1, label = TRUE, label.size = 5)

DefaultAssay(object=MG.cells) <- "RNA"
marker_MG.cells_annotated <- FindAllMarkers(object = MG.cells, only.pos = TRUE, logfc.threshold = 0.25 ,min.pct = 0.25)
marker_MG.cells_annotated <- marker_MG.cells_annotated[!grepl("Rp", rownames(marker_MG.cells_annotated)),]
marker_MG.cells_annotated <- marker_MG.cells_annotated[!grepl("mt", rownames(marker_MG.cells_annotated)),]
marker_MG.cells_annotated <- marker_MG.cells_annotated[!grepl("Mt", rownames(marker_MG.cells_annotated)),]

top15_MG_annotated <- marker_MG.cells_annotated %>%
  group_by(cluster) %>%
  top_n(n = 15, avg_log2FC)
write.csv(MG.cells, file="/Users/zhangyingchen/Desktop/MG.cells.csv")

FeaturePlot(object = MG.cells, features = c("Tmem119","Cd8a"),order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
MG.cells <- subset(MG.cells, idents=c("MG1",  "A.MG",  "MG2",  "MG3",  
                                      "DAM", "Y.MG1",  "Y.MG2", "MG4","MG5"))


#DE aged_sham vs aged_TBI
Idents(object = MG.cells_1) <- "orig.ident"
markers_agedTBI_vs_agedsham <- FindMarkers(object = MG.cells_1, ident.1 = "Aged_TBI", 
                                           ident.2 = "Aged_sham", only.pos = TRUE, min.pct = 0.25, 
                                           logfc.threshold = 0.25)
markers_agedTBI_vs_agedsham <- markers_agedTBI_vs_agedsham[!grepl("Rp", rownames(markers_agedTBI_vs_agedsham)),]
markers_agedTBI_vs_agedsham <- markers_agedTBI_vs_agedsham[!grepl("mt", rownames(markers_agedTBI_vs_agedsham)),]

markers_youngTBI_vs_youngsham <- FindMarkers(object = MG.cells_1, ident.1 = "Young_TBI", 
                                           ident.2 = "Young_sham", only.pos = TRUE, min.pct = 0.25, 
                                           logfc.threshold = 0.25)
markers_youngTBI_vs_youngsham <- markers_youngTBI_vs_youngsham[!grepl("Rp", rownames(markers_youngTBI_vs_youngsham)),]
markers_youngTBI_vs_youngsham <- markers_youngTBI_vs_youngsham[!grepl("Mt", rownames(markers_youngTBI_vs_youngsham)),]

genes <- unique(c(rownames(markers_agedTBI_vs_agedsham), rownames(markers_youngTBI_vs_youngsham)))
common_genes <- intersect(rownames(markers_agedTBI_vs_agedsham), rownames(markers_youngTBI_vs_youngsham))
markers_agedTBI_only <- rownames(markers_agedTBI_vs_agedsham)[!(rownames(markers_agedTBI_vs_agedsham) %in% common_genes)]
markers_youngTBI_only <- rownames(markers_youngTBI_vs_youngsham)[!(rownames(markers_youngTBI_vs_youngsham) %in% common_genes)]

sample_selected <- c("Aged_TBI", "Aged_sham", "Young_TBI","Young_sham")
names(sample_selected) <- c("Aged_TBI", "Aged_sham", "Young_TBI","Young_sham")

genes_mean_expr <- sapply(sample_selected, function(sample1) {
  Matrix::rowMeans(MG.cells_1@assays$RNA@data[genes, 
                                                   colnames(MG.cells_1)[MG.cells_1$orig.ident == sample1]],
                   na.rm = T)
})

genes_mean_expr <- as.data.frame(genes_mean_expr)
genes_mean_expr$gene <- rownames(genes_mean_expr)

#load genes to labeled 
labels<-c("S100a9","S100a8","Ly86","Cd9","Il1a","Ifitm3",
          "Ccr5","Rock1","Lrp1","Pdcd10","Itgam","Ifngr1","Pik3r1", "Ccl5")

genes_mean_expr$color <- 0
genes_mean_expr[genes_mean_expr$gene %in% markers_agedTBI_only, "color"]<-"Aged_TBI"
genes_mean_expr[genes_mean_expr$gene %in% markers_youngTBI_only, "color"]<-"Young_TBI"
genes_mean_expr[genes_mean_expr$gene %in% common_genes, "color"]<-"common"
genes_mean_expr$color<-factor(genes_mean_expr$color, levels=c("Aged_TBI", "Young_TBI", "common"))
genes_mean_expr_labeled <- genes_mean_expr[labels, ]

col_Aged_TBI<-"#9e0303"
col_Young_TBI<- "#ffe5e0"

library(ggrepel)
ggplot(genes_mean_expr, aes(x=Aged_TBI, y=Young_TBI))+
  geom_jitter(aes(fill=color),shape=21, color="grey", alpha=0.7, size=4)+
  geom_text_repel(data=genes_mean_expr_labeled, aes(label=gene), nudge_y=0.2, size=5,
                  direction="both",max.overlaps = Inf)+
  geom_abline(intercept=0, slope=1)+
  scale_fill_manual(values=c(col_Aged_TBI, col_Young_TBI, "black"))+ 
  xlim(0,10)+
  ylim(0,10)+
  xlab("expression level in Aged_TBI")+
  ylab("expression level Young_TBI")+
  coord_fixed()+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())


GO_Aged_TBI_Only_MG <- read.table(file = "/Users/zhangyingchen/Desktop/GOPROCESS_Aged_TBI_only_MG.xls", 
                     sep = "\t", header=TRUE)
ggplot(GO_Aged_TBI_Only_MG,aes(reorder(Description, Enrichment), Enrichment, fill=FDR.q.value, split='GO.Term')) + geom_col(width = 0.7) + 
  theme(axis.text.x=element_text(angle=-40, hjust=0, size = 18), axis.text.y = element_text(size = 16))+
  RotatedAxis()+coord_flip(xlim=c(length(unique(GO_Aged_TBI_Only_MG$GO.Term))-9,length(unique(GO_Aged_TBI_Only_MG$GO.Term))))+scale_fill_distiller(palette = "Reds")+
  theme(axis.title.x = element_text(size = 23), axis.title.y = element_text(size = 23))+
  xlab("Description")

# plot by microglia functions
#Homeostatic Siglech, Gpr34, P2ry12, phagocytosis Tyrobp, Trem2, Clec4a1,Axl, 
#Lipid metab Cst7, Lpl, Apoe, Inflamm Tnf, Irf1, Itgax, Lysosome Cln3, Tpp1, Hexa 
VlnPlot(MG.cells_1,features=c("P2ry12","Tmem119","Tyrobp","Trem2","Clec4a1","Axl","Lpl","Cst7","Apoe","Ccl3","Ccl4","Tnf","Irf1","Itgax","Cln3","Tpp1","Mafb"),stack =TRUE, sort = T,flip = TRUE, assay="SCT", pt.size = 0.2)
VlnPlot(MG.cells,features=c("P2ry12","Tmem119","Siglech","Gpr34"),stack =TRUE, sort = T,flip = TRUE, assay="SCT", pt.size = 0.2)
VlnPlot(MG.cells,features=c("Tyrobp","Trem2","Clec4a1","Clec7a","Axl"),stack =TRUE, sort = T,flip = TRUE, assay="SCT", pt.size = 0.2)
VlnPlot(MG.cells,features=c("Cst7","Lpl","Apoe"),stack =TRUE, sort = T,flip = TRUE, assay="SCT", pt.size = 0.2)
VlnPlot(MG.cells,features=c("Tnf","Irf1","Itgax","Ccl3","Ccl4"),stack =TRUE, sort = T,flip = TRUE, assay="SCT", pt.size = 0.2)
VlnPlot(MG.cells,features=c("Cln3","Tpp1","Hexa"),stack =TRUE, sort = T,flip = TRUE, assay="SCT", pt.size = 0.2)


#
DoHeatmap(avgexp, features = c("P2ry12","Tmem119","Tyrobp","Trem2","Clec4a1","Axl","Lpl","Cst7","Apoe","Ccl3","Ccl4","Tnf","Irf1","Itgax","Cln3","Tpp1","Mafb"),
          group.by = 'MG_8_groups', size = 3)+ scale_fill_gradient2(low = rev(c('#e5f5f9','99d8c9','#2ca25f')), 
                                                              mid = "white", high = rev(c('#c51b8a','#fa9fb5','#fde0dd')),
                                                              midpoint = 0, guide = "colourbar", aesthetics = "fill")

#

FeaturePlot(object = MG.cells, features = c("Itgb5","Sall1","Hexb","Ccl4","Ccl5","Ctss","Aif1","Mertk"),order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
FeaturePlot(object = subset_MG_MAC, features = c("Socs3","Jun","Txnip","Zfp36","Ctss"),order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
FeaturePlot(object = subset_MG_MAC, features = c("Rplp1","Cd209g","Ly6g","Fau","Rps29"),order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)


