# ================================================================
# SingleR-based macrophage annotation pipeline (mouse liver)
# ================================================================
# Author: Jingshuo Li
# Version: 1.0
# Description: End-to-end scRNA-seq workflow for mouse liver with
# DAG-QC, clustering, SingleR annotation (MCA reference), macrophage
# reclustering, marker/heatmap/violin/featureplots, and exports.
# ================================================================

# -------------------- Step 0: Parameters (EDIT HERE) --------------------
project_name <- "MouseLiver_SingleR"
data_dir     <- "data/mouse/liver"                     # raw 10X input
dag_path     <- "data/DAG_136_MouseGenes.csv"          # DAG list (first column = gene symbol)
ref_path     <- "refs/mouse/mca_liver_ref.rds"         # MCA liver SingleR reference (RDS)
results_dir  <- "results/mouse/liver/singler_pipeline" # output root
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- Step 1: Libraries --------------------
suppressPackageStartupMessages({
  library(GEOquery)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(openxlsx)
  library(cluster)
  library(SingleR)
  library(SingleCellExperiment)
  library(patchwork)
  library(plyr)
})
set.seed(111)

# -------------------- Step 2: GEO metadata (optional) --------------------
# Note: Will fail silently without internet
gsm_mouse_liver <- tryCatch({
  getGEO("GSM5687567", GSEMatrix = FALSE)
}, error = function(e) { message("GEO fetch skipped."); NULL })
if (!is.null(gsm_mouse_liver)) utils::str(gsm_mouse_liver)

# -------------------- Step 3: Load raw 10X gene expression data --------------------
mouse_liver <- Read10X(data.dir = data_dir) %>%
  CreateSeuratObject(project = project_name, min.cells = 3, min.features = 200)

# -------------------- Step 4: Load mouse DAG gene list --------------------
dag_df       <- read.csv(dag_path, stringsAsFactors = FALSE)
dag_genes    <- unique(dag_df[[1]])
overlap_dag  <- intersect(rownames(mouse_liver), dag_genes)

# -------------------- Step 5: Add QC metrics --------------------
mouse_liver[["percent.mito"]] <- PercentageFeatureSet(mouse_liver, pattern = "^mt-")
mouse_liver[["percent.dag"]]  <- PercentageFeatureSet(mouse_liver, features = overlap_dag)

print(summary(mouse_liver$percent.mito))
print(summary(mouse_liver$percent.dag))

qc_hist_dir <- file.path(results_dir, "QC_Histograms")
dir.create(qc_hist_dir, showWarnings = FALSE, recursive = TRUE)
png(file.path(qc_hist_dir, "Histogram_percent.mito.png"), width = 800, height = 600)
hist(mouse_liver$percent.mito, breaks = 50, col = "lightblue", main = "Mito %", xlab = "percent.mito"); dev.off()
png(file.path(qc_hist_dir, "Histogram_percent.dag.png"), width = 800, height = 600)
hist(mouse_liver$percent.dag, breaks = 50, col = "pink", main = "DAG %", xlab = "percent.dag"); dev.off()

# -------------------- Step 6: QC scatter plots --------------------
output_dir_qc <- file.path(results_dir, "QC_Plots")
dir.create(output_dir_qc, recursive = TRUE, showWarnings = FALSE)

plot1 <- FeatureScatter(mouse_liver, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(mouse_liver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(mouse_liver, feature1 = "nCount_RNA", feature2 = "percent.dag")
plot4 <- FeatureScatter(mouse_liver, feature1 = "percent.dag", feature2 = "percent.mito")

ggsave(file.path(output_dir_qc, "Scatter_nCount_vs_percent.mito.png"), plot1, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir_qc, "Scatter_nCount_vs_nFeature.png"),  plot2, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir_qc, "Scatter_nCount_vs_percent.dag.png"),plot3, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir_qc, "Scatter_percent.dag_vs_percent.mito.png"), plot4, width = 8, height = 8, dpi = 300, bg = "white")
rm(plot1, plot2, plot3, plot4); gc()

# -------------------- Step 7: Final QC filtering --------------------
mouse_liver_filtered <- subset(mouse_liver,
                               subset = nFeature_RNA > 500 &
                                 nFeature_RNA < 4500 &
                                 percent.mito   < 15 &
                                 percent.dag    < 5)
saveRDS(mouse_liver_filtered, file = file.path(results_dir, "mouse_liver_filtered.rds"))

# -------------------- Step 8: Normalization --------------------
mouse_liver_filtered <- NormalizeData(mouse_liver_filtered)

# -------------------- Step 9: Join layers (Seurat v5) --------------------
mouse_liver_filtered <- JoinLayers(mouse_liver_filtered)

# -------------------- Step 10: HVGs (mean.var.plot) --------------------
mouse_liver_meanvar <- FindVariableFeatures(
  mouse_liver_filtered,
  selection.method = "mean.var.plot",
  mean.cutoff = c(0.04, 3),
  dispersion.cutoff = c(0.4, Inf)
)
message("Variable features: ", length(VariableFeatures(mouse_liver_meanvar)))

output_dir_meanvar <- file.path(results_dir, "Method_meanvar")
dir.create(output_dir_meanvar, recursive = TRUE, showWarnings = FALSE)

# -------------------- Step 11: Scale + PCA --------------------
mouse_liver_meanvar <- ScaleData(mouse_liver_meanvar,
                                 vars.to.regress = c("nCount_RNA", "percent.mito"))
mouse_liver_meanvar <- RunPCA(mouse_liver_meanvar,
                              features = VariableFeatures(mouse_liver_meanvar))

# -------------------- Step 12: Elbow plot --------------------
p_elbow_meanvar <- ElbowPlot(mouse_liver_meanvar) + theme_classic()
ggsave(filename = file.path(output_dir_meanvar, "ElbowPlot_meanvar.png"),
       plot = p_elbow_meanvar, width = 6, height = 5, dpi = 300)
rm(p_elbow_meanvar)

# -------------------- Step 13: Neighbors + Clusters --------------------
mouse_liver_meanvar <- FindNeighbors(mouse_liver_meanvar, dims = 1:15)
mouse_liver_meanvar <- FindClusters(mouse_liver_meanvar, resolution = 1.0)

Idents(mouse_liver_meanvar) <- mouse_liver_meanvar$seurat_clusters
cluster_meanvar_ids <- Idents(mouse_liver_meanvar)
levels(cluster_meanvar_ids) <- as.character(as.numeric(levels(cluster_meanvar_ids)) + 1)
mouse_liver_meanvar$Cluster1based <- cluster_meanvar_ids
Idents(mouse_liver_meanvar) <- mouse_liver_meanvar$Cluster1based

# -------------------- Step 14: UMAP --------------------
mouse_liver_meanvar <- RunUMAP(mouse_liver_meanvar, reduction = "pca", dims = 1:15)

# -------------------- Step 15: tSNE --------------------
mouse_liver_meanvar <- RunTSNE(mouse_liver_meanvar, reduction = "pca", dims = 1:15)

# -------------------- Step 16: Save UMAP/tSNE (cluster labels) --------------------
p_umap_meanvar <- DimPlot(mouse_liver_meanvar, reduction = "umap", label = TRUE) + theme_classic()
ggsave(filename = file.path(output_dir_meanvar, "UMAP_FirstClusters_meanvar.png"),
       plot = p_umap_meanvar, width = 8, height = 6, dpi = 300)
rm(p_umap_meanvar)

p_tsne_meanvar <- DimPlot(mouse_liver_meanvar, reduction = "tsne", label = TRUE) + theme_classic()
ggsave(filename = file.path(output_dir_meanvar, "TSNE_FirstClusters_meanvar.png"),
       plot = p_tsne_meanvar, width = 8, height = 6, dpi = 300)
rm(p_tsne_meanvar)

# -------------------- Step 17: Silhouette score --------------------
emb_pca <- Embeddings(mouse_liver_meanvar, reduction = "pca")[, 1:15]
cluster_labels <- as.factor(mouse_liver_meanvar$Cluster1based)
dists <- dist(emb_pca)
sil <- silhouette(as.integer(cluster_labels), dists)

sil_df <- as.data.frame(sil)
sil_df$Cluster <- cluster_labels[as.numeric(rownames(sil_df))]
sil_summary <- sil_df %>% group_by(Cluster) %>% summarise(mean_sil = mean(sil_width)) %>% arrange(desc(mean_sil))
print(sil_summary)

output_dir_silhouette <- file.path(results_dir, "SilhouettePlot")
dir.create(output_dir_silhouette, showWarnings = FALSE, recursive = TRUE)
p_sil <- ggplot(sil_df, aes(x = Cluster, y = sil_width, fill = Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.4) + geom_jitter(width = 0.2, size = 0.3, alpha = 0.5) +
  theme_classic() + labs(y = "Silhouette width", title = "Silhouette Score per Cluster")
ggsave(filename = file.path(output_dir_silhouette, "Silhouette_PerCluster.png"), plot = p_sil, width = 8, height = 6, dpi = 300)

# ================================= SingleR ==================================
# -------------------- Step 18: SingleR annotation (MCA liver ref) --------------------
output_dir_singler <- file.path(results_dir, "AutoAnnotation_SingleR")
dir.create(output_dir_singler, showWarnings = FALSE, recursive = TRUE)

sce <- as.SingleCellExperiment(mouse_liver_meanvar)
ref_mca_liver <- readRDS(ref_path)
pred_mca <- SingleR(test = sce, ref = ref_mca_liver, labels = ref_mca_liver$label.main)
mouse_liver_meanvar$MCA_liver_label <- pred_mca$labels

# -------------------- Step 19: Save SingleR outputs & plots --------------------
pred_df <- as.data.frame(pred_mca)
write.csv(pred_df, file.path(output_dir_singler, "SingleR_MCA_liver_Result.csv"), row.names = TRUE)

Idents(mouse_liver_meanvar) <- mouse_liver_meanvar$MCA_liver_label
p_umap_mca <- DimPlot(mouse_liver_meanvar, reduction = "umap", label = FALSE, repel = TRUE) + ggtitle("SingleR UMAP – MCA liver") + theme_classic()
p_tsne_mca <- DimPlot(mouse_liver_meanvar, reduction = "tsne", label = FALSE, repel = TRUE) + ggtitle("SingleR tSNE – MCA liver") + theme_classic()
ggsave(file.path(output_dir_singler, "UMAP_SingleR_MCA_liver.png"), p_umap_mca, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir_singler, "tSNE_SingleR_MCA_liver.png"), p_tsne_mca, width = 8, height = 6, dpi = 300)
rm(p_umap_mca, p_tsne_mca)

saveRDS(mouse_liver_meanvar, file = file.path(output_dir_singler, "mouse_liver_annotated_MCA_liver.rds"))
utils::write.csv(as.data.frame(table(mouse_liver_meanvar$MCA_liver_label)),
                 file = file.path(output_dir_singler, "Counts_by_SingleR_Label.csv"), row.names = FALSE)

# -------------------- Step 20: DEGs per SingleR label --------------------
Idents(mouse_liver_meanvar) <- mouse_liver_meanvar$MCA_liver_label
deg_results <- FindAllMarkers(object = mouse_liver_meanvar, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(deg_results, file = file.path(output_dir_singler, "All_DEGs_SingleR_Celltypes.csv"), row.names = FALSE)

top30_deg <- deg_results %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 30)
write.csv(top30_deg, file = file.path(output_dir_singler, "Top30_DEGs_SingleR_Celltypes.csv"), row.names = FALSE)

# -------------------- Step 21: Violin (Kupffer cell) --------------------
vln_dir1 <- file.path(output_dir_singler, "Violin_Markers_KupfferCheck")
dir.create(vln_dir1, showWarnings = FALSE, recursive = TRUE)

subset_kupffer <- subset(mouse_liver_meanvar, subset = MCA_liver_label == "Kupffer cell")
subset_kupffer$MCA_liver_label <- droplevels(factor(subset_kupffer$MCA_liver_label))

marker_genes1 <- c("Csf1r","Cd68","Adgre1","Fcgr1","Mrc1","Timd4","Folr2","Lyve1","Ccr2","Plac8","S100a8","S100a9","H2-Aa","H2-Ab1")
plot_list <- lapply(marker_genes1, function(gene) {
  VlnPlot(subset_kupffer, features = gene, group.by = "MCA_liver_label", pt.size = 0) +
    ggtitle(gene) + theme_classic(base_size = 12) +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
})
p_combined <- wrap_plots(plotlist = plot_list, ncol = 5)
ggsave(filename = file.path(vln_dir1, "Violin_AllMarkers_Kupffer_Patched.png"), plot = p_combined, width = 30, height = 20, dpi = 600)
rm(p_combined); gc()

# -------------------- Step 22: Violin (four macrophage-related populations) --------------------
vln_dir2 <- file.path(output_dir_singler, "Violin_Markers_FourCell")
dir.create(vln_dir2, showWarnings = FALSE, recursive = TRUE)

four_cell_for_violin <- subset(mouse_liver_meanvar, idents = c("Kupffer cell","Macrophage_Chil3 high","Proliferating Macrophage","Monocyte"))
four_cell_for_violin$MCA_liver_label <- droplevels(factor(four_cell_for_violin$MCA_liver_label))
Idents(four_cell_for_violin) <- four_cell_for_violin$MCA_liver_label

marker_genes2 <- marker_genes1
for (gene in marker_genes2) {
  p <- VlnPlot(four_cell_for_violin, features = gene, group.by = "MCA_liver_label", pt.size = 0) +
    ggtitle(gene) + theme_classic(base_size = 12) +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks.x = element_blank())
  ggsave(filename = file.path(vln_dir2, paste0("Violin_", gene, ".png")), plot = p, width = 8, height = 8, dpi = 300)
}

# ============================ Macrophage Focus =============================
# -------------------- Step 23: Macrophage reclustering --------------------
output_dir_recluster <- file.path(results_dir, "recluster")
dir.create(output_dir_recluster, recursive = TRUE, showWarnings = FALSE)

Idents(mouse_liver_meanvar) <- mouse_liver_meanvar$MCA_liver_label
macro_subset <- subset(mouse_liver_meanvar, idents = c("Kupffer cell","Macrophage_Chil3 high"))
saveRDS(macro_subset, file = file.path(output_dir_recluster, "MouseLiver_MCA_Macrophages.rds"))

macro_subset <- FindVariableFeatures(macro_subset, selection.method = "mean.var.plot", mean.cutoff = c(0.04, 3), dispersion.cutoff = c(0.4, Inf))
message("Macrophage variable features: ", length(VariableFeatures(macro_subset)))

p_elbow_macro <- ElbowPlot(macro_subset) + theme_classic()
ggsave(filename = file.path(output_dir_recluster, "ElbowPlot_macro_subset.png"), plot = p_elbow_macro, width = 6, height = 5, dpi = 300)
rm(p_elbow_macro); gc()

macro_subset <- ScaleData(macro_subset)
macro_subset <- RunPCA(macro_subset, features = VariableFeatures(macro_subset))
macro_subset <- FindNeighbors(macro_subset, dims = 1:15)
macro_subset <- FindClusters(macro_subset, resolution = 0.25)

Idents(macro_subset) <- macro_subset$seurat_clusters
levels(Idents(macro_subset)) <- as.character(as.numeric(levels(Idents(macro_subset))) + 1)
macro_subset$MF_subtype <- Idents(macro_subset)
Idents(macro_subset) <- macro_subset$MF_subtype

macro_subset <- RunUMAP(macro_subset, dims = 1:15)
macro_subset <- RunTSNE(macro_subset, dims = 1:15)

ggsave(file.path(output_dir_recluster, "UMAP_MF_Subtype_Annotated.png"),
       DimPlot(macro_subset, reduction = "umap", group.by = "MF_subtype", label = TRUE, repel = TRUE) + theme_classic(),
       width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir_recluster, "tSNE_MF_Subtype_Annotated.png"),
       DimPlot(macro_subset, reduction = "tsne", group.by = "MF_subtype", label = TRUE, repel = TRUE) + theme_classic(),
       width = 8, height = 6, dpi = 300)

# -------------------- Step 24: MF subtype marker detection --------------------
mf_markers <- data.frame()
for (clust in levels(Idents(macro_subset))) {
  markers <- FindMarkers(macro_subset, ident.1 = clust, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2)
  if (nrow(markers) > 0) { markers$cluster <- clust; markers$gene <- rownames(markers); mf_markers <- rbind(mf_markers, markers) }
}
mf_sig <- mf_markers %>% filter(p_val_adj < 0.05) %>% filter(!grepl("^mt-", gene, ignore.case = TRUE)) %>% filter(!grepl("^Rp[sl]", gene, ignore.case = TRUE))
top30_mf <- mf_sig %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 30)
write.xlsx(mf_markers, file.path(output_dir_recluster, "DEGs_MFSubtype_Markers.xlsx"))
write.xlsx(top30_mf,  file.path(output_dir_recluster, "Top30_MFSubtype_Markers.xlsx"))

# -------------------- Step 25: Violin plots (MF subtype markers) --------------------
vln_dir3 <- file.path(output_dir_recluster, "Violin_MFSubtype_Markers")
dir.create(vln_dir3, showWarnings = FALSE, recursive = TRUE)
genes_to_plot <- c("Csf1r","Cd68","Adgre1","Fcgr1","Mrc1","Timd4","Folr2","Lyve1","Ccr2","Plac8","S100a8","H2-Aa","H2-Ab1","Ccnb2","Mki67","Clec4f","Vsig4","Marco","Cd5l","Cd3d")
for (gene in genes_to_plot) {
  p_vln <- VlnPlot(macro_subset, features = gene, group.by = "MF_subtype", pt.size = 0) + ggtitle(gene) + theme_classic()
  ggsave(filename = file.path(vln_dir3, paste0("Violin_", gene, ".png")), plot = p_vln, width = 5, height = 4, dpi = 300)
}

# -------------------- Step 26: Heatmap (Top30 per MF subtype) --------------------
top30_path <- file.path(output_dir_recluster, "Top30_MFSubtype_Markers.xlsx")
df_top30 <- read.xlsx(top30_path)
df_top30$cluster <- as.character(df_top30$cluster)
top30_genes <- unique(df_top30$gene)

set.seed(111)
subtypes <- unique(df_top30$cluster)
cells_use <- c()
for (stype in subtypes) {
  cells <- WhichCells(macro_subset, expression = MF_subtype == stype)
  if (length(cells) > 0) cells_use <- c(cells_use, sample(cells, size = min(50, length(cells))))
}
subset_obj <- subset(macro_subset, cells = cells_use)
DefaultAssay(subset_obj) <- "RNA"
subset_obj <- ScaleData(subset_obj, features = top30_genes, verbose = FALSE)
subset_obj$MF_subtype <- factor(subset_obj$MF_subtype, levels = subtypes)
Idents(subset_obj) <- subset_obj$MF_subtype

p_heat <- DoHeatmap(object = subset_obj, features = top30_genes) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 5), strip.text.y = element_text(size = 9, face = "bold")) +
  ggtitle("Top 30 Marker Genes per MF Subtype (Liver)")
ggsave(filename = file.path(output_dir_recluster, "Heatmap_Top30_MFSubtype_50cells.png"), plot = p_heat, width = 13, height = 14, dpi = 600)

# -------------------- Step 27: Assign MF subtype final labels (manual) --------------------
mf_labels <- c("1" = "Kupffer Macrophage", "2" = "Kupffer Macrophage", "3" = "Proliferative Macrophage", "4" = "TLF⁺ Macrophage")
macro_subset$MF_subtype_label <- plyr::mapvalues(x = Idents(macro_subset), from = names(mf_labels), to = mf_labels)
Idents(macro_subset) <- macro_subset$MF_subtype_label

ggsave(filename = file.path(output_dir_recluster, "UMAP_MF_Subtype_ManualLabel.png"),
       plot = DimPlot(macro_subset, group.by = "MF_subtype_label", label = TRUE, repel = TRUE) + theme_classic(), width = 8, height = 6, dpi = 300)
ggsave(filename = file.path(output_dir_recluster, "tSNE_MF_Subtype_ManualLabel.png"),
       plot = DimPlot(macro_subset, reduction = "tsne", group.by = "MF_subtype_label", label = TRUE, repel = TRUE) + theme_classic(), width = 8, height = 6, dpi = 300)

# -------------------- Step 28: FeaturePlots for key markers --------------------
key_markers <- c("Timd4","Folr2","Lyve1","Ccr2","Plac8","S100a8","H2-Aa","H2-Ab1","Ccnb2","Mki67","Clec4f","Vsig4","Marco","Cd5l")
p_feat_umap <- FeaturePlot(macro_subset, features = key_markers, reduction = "umap", order = TRUE, cols = c("lightgrey","red")) & theme_classic()
ggsave(filename = file.path(output_dir_recluster, "FeaturePlot_MFSubtypeMarkers_UMAP.png"), plot = p_feat_umap, width = 20, height = 20, dpi = 300)
p_feat_tsne <- FeaturePlot(macro_subset, features = key_markers, reduction = "tsne", order = TRUE, cols = c("lightgrey","red")) & theme_classic()
ggsave(filename = file.path(output_dir_recluster, "FeaturePlot_MFSubtypeMarkers_tSNE.png"), plot = p_feat_tsne, width = 20, height = 20, dpi = 300)

# -------------------- Step 29: DEGs & Heatmap by final labels + Save object --------------------
Idents(macro_subset) <- macro_subset$MF_subtype_label
deg_mf <- FindAllMarkers(object = macro_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top30_mf_final <- deg_mf %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 30)
write.xlsx(deg_mf,        file.path(output_dir_recluster, "DEGs_Reassigned_MFSubtype_All.xlsx"))
write.xlsx(top30_mf_final, file.path(output_dir_recluster, "Top30_Reassigned_MFSubtype.xlsx"))

top30_path_final <- file.path(output_dir_recluster, "Top30_Reassigned_MFSubtype.xlsx")
df_top30_final   <- read.xlsx(top30_path_final)
df_top30_final$cluster <- as.character(df_top30_final$cluster)
top30_genes_final <- unique(df_top30_final$gene)

set.seed(111)
subtypes_final <- unique(df_top30_final$cluster)
cells_use2 <- c()
for (stype in subtypes_final) {
  cells <- WhichCells(macro_subset, expression = MF_subtype_label == stype)
  if (length(cells) > 0) cells_use2 <- c(cells_use2, sample(cells, size = min(50, length(cells))))
}
subset_obj2 <- subset(macro_subset, cells = cells_use2)
DefaultAssay(subset_obj2) <- "RNA"
subset_obj2 <- ScaleData(subset_obj2, features = top30_genes_final, verbose = FALSE)
subset_obj2$MF_subtype_label <- factor(subset_obj2$MF_subtype_label, levels = subtypes_final)
Idents(subset_obj2) <- subset_obj2$MF_subtype_label

p_heat2 <- DoHeatmap(object = subset_obj2, features = top30_genes_final) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 5), strip.text.y = element_text(size = 10, face = "bold"))
ggsave(filename = file.path(output_dir_recluster, "Heatmap_Top30_Reassigned_MFSubtype_50cells.png"), plot = p_heat2, width = 13, height = 14, dpi = 600)

saveRDS(macro_subset, file = file.path(output_dir_recluster, "MouseLiver_MCA_MFSubtype_Final.rds"))
# ============================ END ================================
