# ================================================================
# SingleR-based macrophage annotation pipeline (mouse lung)
# ================================================================
# Author: Jingshuo Li
# Version: 1.0
# Description: End-to-end scRNA-seq workflow for mouse lung with
# DAG-QC, clustering, SingleR annotation (MCA reference), macrophage
# reclustering, marker/heatmap/violin/featureplots, and exports.
# ================================================================

# -------------------- Step 0: Parameters (EDIT HERE) --------------------
project_name <- "MouseLung_SingleR"
data_dir     <- "data/mouse/lung"                        # raw 10X input
dag_path     <- "data/DAG_136_MouseGenes.csv"            # DAG list (first column = gene symbol)
ref_path     <- "refs/mouse/mca_lung_ref.rds"            # MCA lung SingleR reference (RDS)
results_dir  <- "results/mouse/lung/singler_pipeline"    # output root

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
gsm_mouse_lung <- tryCatch({
  getGEO("GSM5687568", GSEMatrix = FALSE)
}, error = function(e) { message("GEO fetch skipped."); NULL })
if (!is.null(gsm_mouse_lung)) utils::str(gsm_mouse_lung)

# -------------------- Step 3: Load raw 10X gene expression data --------------------
mouse_lung <- Read10X(data.dir = data_dir) %>%
  CreateSeuratObject(project = project_name, min.cells = 3, min.features = 200)

# -------------------- Step 4: Load mouse DAG gene list --------------------
dag_df    <- read.csv(dag_path, stringsAsFactors = FALSE)
dag_genes <- unique(dag_df[[1]])
overlap_dag <- intersect(rownames(mouse_lung), dag_genes)

# -------------------- Step 5: Add QC metrics --------------------
# mouse mito genes use "^mt-"
mouse_lung[["percent.mito"]] <- PercentageFeatureSet(mouse_lung, pattern = "^mt-")
mouse_lung[["percent.dag"]]  <- PercentageFeatureSet(mouse_lung, features = overlap_dag)

# Display distributions
print(summary(mouse_lung$percent.mito))
print(summary(mouse_lung$percent.dag))

# QC histograms
qc_hist_dir <- file.path(results_dir, "QC_Histograms")
dir.create(qc_hist_dir, showWarnings = FALSE, recursive = TRUE)

png(file.path(qc_hist_dir, "Histogram_percent.mito.png"), width = 800, height = 600)
hist(mouse_lung$percent.mito, breaks = 50, col = "lightblue", main = "Mito %", xlab = "percent.mito")
dev.off()

png(file.path(qc_hist_dir, "Histogram_percent.dag.png"), width = 800, height = 600)
hist(mouse_lung$percent.dag, breaks = 50, col = "pink", main = "DAG %", xlab = "percent.dag")
dev.off()

# -------------------- Step 6: QC scatter plots --------------------
output_dir_qc <- file.path(results_dir, "QC_Plots")
dir.create(output_dir_qc, recursive = TRUE, showWarnings = FALSE)

plot1 <- FeatureScatter(mouse_lung, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(mouse_lung, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(mouse_lung, feature1 = "nCount_RNA", feature2 = "percent.dag")
plot4 <- FeatureScatter(mouse_lung, feature1 = "percent.dag", feature2 = "percent.mito")

ggsave(file.path(output_dir_qc, "Scatter_nCount_vs_percent.mito.png"), plot1, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir_qc, "Scatter_nCount_vs_nFeature.png"),  plot2, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir_qc, "Scatter_nCount_vs_percent.dag.png"),plot3, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir_qc, "Scatter_percent.dag_vs_percent.mito.png"), plot4, width = 8, height = 8, dpi = 300, bg = "white")

rm(plot1, plot2, plot3, plot4); gc()

# -------------------- Step 7: Final QC filtering --------------------
mouse_lung_filtered <- subset(mouse_lung,
                              subset = nFeature_RNA > 200 &
                                nFeature_RNA < 6000 &
                                percent.mito   < 15 &
                                percent.dag    < 5)

saveRDS(mouse_lung_filtered, file = file.path(results_dir, "mouse_lung_filtered.rds"))

# -------------------- Step 8: Normalization --------------------
mouse_lung_filtered <- NormalizeData(mouse_lung_filtered)

# -------------------- Step 9: Join layers (Seurat v5) --------------------
mouse_lung_filtered <- JoinLayers(mouse_lung_filtered)

# -------------------- Step 10: HVGs (mean.var.plot) --------------------
mouse_lung_meanvar <- FindVariableFeatures(
  mouse_lung_filtered,
  selection.method = "mean.var.plot",
  mean.cutoff = c(0.005, 4),
  dispersion.cutoff = c(0.25, Inf)
)
message("Variable features: ", length(VariableFeatures(mouse_lung_meanvar)))

output_dir_meanvar <- file.path(results_dir, "Method_meanvar")
dir.create(output_dir_meanvar, recursive = TRUE, showWarnings = FALSE)

# -------------------- Step 11: Scale + PCA --------------------
mouse_lung_meanvar <- ScaleData(mouse_lung_meanvar,
                                vars.to.regress = c("nCount_RNA", "percent.mito"))
mouse_lung_meanvar <- RunPCA(mouse_lung_meanvar,
                             features = VariableFeatures(mouse_lung_meanvar))

# -------------------- Step 12: Elbow plot --------------------
p_elbow_meanvar <- ElbowPlot(mouse_lung_meanvar) + theme_classic()
ggsave(filename = file.path(output_dir_meanvar, "ElbowPlot_meanvar.png"),
       plot = p_elbow_meanvar, width = 6, height = 5, dpi = 300)
rm(p_elbow_meanvar)

# -------------------- Step 13: Neighbors + Clusters --------------------
mouse_lung_meanvar <- FindNeighbors(mouse_lung_meanvar, dims = 1:20)
mouse_lung_meanvar <- FindClusters(mouse_lung_meanvar, resolution = 1.0)

# Use seurat clusters, then rename from 1
Idents(mouse_lung_meanvar) <- mouse_lung_meanvar$seurat_clusters
cluster_meanvar_ids <- Idents(mouse_lung_meanvar)
levels(cluster_meanvar_ids) <- as.character(as.numeric(levels(cluster_meanvar_ids)) + 1)
mouse_lung_meanvar$Cluster1based <- cluster_meanvar_ids
Idents(mouse_lung_meanvar) <- mouse_lung_meanvar$Cluster1based

# -------------------- Step 14: UMAP --------------------
mouse_lung_meanvar <- RunUMAP(mouse_lung_meanvar, reduction = "pca", dims = 1:20)

# -------------------- Step 15: tSNE --------------------
mouse_lung_meanvar <- RunTSNE(mouse_lung_meanvar, reduction = "pca", dims = 1:20)

# -------------------- Step 16: Save UMAP/tSNE (cluster labels) --------------------
p_umap_meanvar <- DimPlot(mouse_lung_meanvar, reduction = "umap", label = TRUE) + theme_classic()
ggsave(filename = file.path(output_dir_meanvar, "UMAP_FirstClusters_meanvar.png"),
       plot = p_umap_meanvar, width = 8, height = 6, dpi = 300)
rm(p_umap_meanvar)

p_tsne_meanvar <- DimPlot(mouse_lung_meanvar, reduction = "tsne", label = TRUE) + theme_classic()
ggsave(filename = file.path(output_dir_meanvar, "TSNE_FirstClusters_meanvar.png"),
       plot = p_tsne_meanvar, width = 8, height = 6, dpi = 300)
rm(p_tsne_meanvar)

# -------------------- Step 17: Silhouette score --------------------
emb_pca <- Embeddings(mouse_lung_meanvar, reduction = "pca")[, 1:20]
cluster_labels <- as.factor(mouse_lung_meanvar$Cluster1based)
dists <- dist(emb_pca)
sil <- silhouette(as.integer(cluster_labels), dists)

sil_df <- as.data.frame(sil)
sil_df$Cluster <- cluster_labels[as.numeric(rownames(sil_df))]
sil_summary <- sil_df %>%
  group_by(Cluster) %>%
  summarise(mean_sil = mean(sil_width)) %>%
  arrange(desc(mean_sil))
print(sil_summary)

output_dir_silhouette <- file.path(results_dir, "SilhouettePlot")
dir.create(output_dir_silhouette, showWarnings = FALSE, recursive = TRUE)

p_sil <- ggplot(sil_df, aes(x = Cluster, y = sil_width, fill = Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.5) +
  theme_classic() +
  labs(y = "Silhouette width", title = "Silhouette Score per Cluster")

ggsave(filename = file.path(output_dir_silhouette, "Silhouette_PerCluster.png"),
       plot = p_sil, width = 8, height = 6, dpi = 300)

# ================================= SingleR ==================================
# -------------------- Step 18: SingleR annotation (MCA lung ref) --------------------
output_dir_singler <- file.path(results_dir, "AutoAnnotation_SingleR")
dir.create(output_dir_singler, showWarnings = FALSE, recursive = TRUE)

# Convert to SCE
sce <- as.SingleCellExperiment(mouse_lung_meanvar)

# Load reference and annotate
ref_mca_lung <- readRDS(ref_path)
pred_mca <- SingleR(test = sce, ref = ref_mca_lung, labels = ref_mca_lung$label.main)
mouse_lung_meanvar$MCA_Lung_label <- pred_mca$labels

# -------------------- Step 19: Save SingleR outputs & plots --------------------
pred_df <- as.data.frame(pred_mca)
write.csv(pred_df, file.path(output_dir_singler, "SingleR_MCA_Lung_Result.csv"), row.names = TRUE)

Idents(mouse_lung_meanvar) <- mouse_lung_meanvar$MCA_Lung_label
p_umap_mca <- DimPlot(mouse_lung_meanvar, reduction = "umap", label = FALSE, repel = TRUE) +
  ggtitle("SingleR UMAP – MCA Lung") + theme_classic()
p_tsne_mca <- DimPlot(mouse_lung_meanvar, reduction = "tsne", label = FALSE, repel = TRUE) +
  ggtitle("SingleR tSNE – MCA Lung") + theme_classic()
ggsave(file.path(output_dir_singler, "UMAP_SingleR_MCA_Lung.png"), p_umap_mca, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir_singler, "tSNE_SingleR_MCA_Lung.png"), p_tsne_mca, width = 8, height = 6, dpi = 300)

saveRDS(mouse_lung_meanvar, file = file.path(output_dir_singler, "mouse_lung_annotated_MCA_Lung.rds"))
utils::write.csv(as.data.frame(table(mouse_lung_meanvar$MCA_Lung_label)),
                 file = file.path(output_dir_singler, "Counts_by_SingleR_Label.csv"), row.names = FALSE)

# -------------------- Step 20: DEGs per SingleR label --------------------
Idents(mouse_lung_meanvar) <- mouse_lung_meanvar$MCA_Lung_label
deg_results <- FindAllMarkers(
  object = mouse_lung_meanvar,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(deg_results, file = file.path(output_dir_singler, "All_DEGs_SingleR_Celltypes.csv"), row.names = FALSE)

top30_deg <- deg_results %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30)
write.csv(top30_deg, file = file.path(output_dir_singler, "Top30_DEGs_SingleR_Celltypes.csv"), row.names = FALSE)

# -------------------- Step 21: Violin (Myeloid cell_Elane high) --------------------
vln_dir1 <- file.path(output_dir_singler, "Violin_Markers_MyeloidCheck")
dir.create(vln_dir1, showWarnings = FALSE, recursive = TRUE)

subset_myeloid <- subset(mouse_lung_meanvar, subset = MCA_Lung_label == "Myeloid cell_Elane high")
subset_myeloid$MCA_Lung_label <- droplevels(factor(subset_myeloid$MCA_Lung_label))

marker_genes1 <- c(
  "Csf1r", "Cd68", "Adgre1", "Fcgr1", "Ear1",
  "Ccr2", "Ly6c2", "Plac8", "Il1b",
  "S100a8", "S100a9",
  "H2-Aa", "H2-Ab1"
)

plot_list <- lapply(marker_genes1, function(gene) {
  VlnPlot(
    subset_myeloid,
    features = gene,
    group.by = "MCA_Lung_label",
    pt.size = 0
  ) +
    ggtitle(gene) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
})
p_combined <- wrap_plots(plotlist = plot_list, ncol = 5)
ggsave(
  filename = file.path(vln_dir1, "Violin_AllMarkers_ElaneHigh_Patched.png"),
  plot = p_combined, width = 30, height = 20, dpi = 600
)
rm(p_combined); gc()

# -------------------- Step 22: Violin (four macrophage-related populations) --------------------
vln_dir2 <- file.path(output_dir_singler, "Violin_Markers_FourCell")
dir.create(vln_dir2, showWarnings = FALSE, recursive = TRUE)

four_cell_for_violin <- subset(mouse_lung_meanvar, idents = c(
  "Alveolar macrophage_Ear2 high",
  "Alveolar macrophage_Pclaf high",
  "Myeloid cell_Elane high",
  "Interstitial macrophage"
))
four_cell_for_violin$MCA_Lung_label <- droplevels(factor(four_cell_for_violin$MCA_Lung_label))
Idents(four_cell_for_violin) <- four_cell_for_violin$MCA_Lung_label

marker_genes2 <- c(
  "Csf1r", "Cd68", "Adgre1", "Fcgr1", "Ear1",
  "Ccr2", "Ly6c2", "Plac8", "Il1b",
  "S100a8", "S100a9",
  "H2-Aa", "H2-Ab1"
)

for (gene in marker_genes2) {
  p <- VlnPlot(
    four_cell_for_violin,
    features = gene,
    group.by = "MCA_Lung_label",
    pt.size = 0
  ) +
    ggtitle(gene) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks.x = element_blank()
    )
  ggsave(
    filename = file.path(vln_dir2, paste0("Violin_", gene, ".png")),
    plot = p, width = 8, height = 8, dpi = 300
  )
}

# ============================ Macrophage Focus =============================
# -------------------- Step 23: Macrophage reclustering --------------------
output_dir_recluster <- file.path(results_dir, "recluster")
dir.create(output_dir_recluster, recursive = TRUE, showWarnings = FALSE)

# Identity: SingleR labels
Idents(mouse_lung_meanvar) <- mouse_lung_meanvar$MCA_Lung_label

# NOTE: this selects only Interstitial macrophage (follow your code)
macro_subset <- subset(mouse_lung_meanvar, idents = c("Interstitial macrophage"))
saveRDS(macro_subset, file = file.path(output_dir_recluster, "Mouselung_MCA_Macrophages.rds"))

# HVGs, Scale, PCA, neighbors/clusters, UMAP/TSNE
macro_subset <- FindVariableFeatures(
  macro_subset,
  selection.method = "mean.var.plot",
  mean.cutoff = c(0.02, 4),
  dispersion.cutoff = c(0.3, Inf)
)
message("Macrophage variable features: ", length(VariableFeatures(macro_subset)))

p_elbow_macro <- ElbowPlot(macro_subset) + theme_classic()
ggsave(
  filename = file.path(output_dir_recluster, "ElbowPlot_macro_subset.png"),
  plot = p_elbow_macro, width = 6, height = 5, dpi = 300
)
rm(p_elbow_macro); gc()

macro_subset <- ScaleData(macro_subset)
macro_subset <- RunPCA(macro_subset, features = VariableFeatures(macro_subset))
macro_subset <- FindNeighbors(macro_subset, dims = 1:15)
macro_subset <- FindClusters(macro_subset, resolution = 0.4)

Idents(macro_subset) <- macro_subset$seurat_clusters
levels(Idents(macro_subset)) <- as.character(as.numeric(levels(Idents(macro_subset))) + 1)
macro_subset$MF_subtype <- Idents(macro_subset)
Idents(macro_subset) <- macro_subset$MF_subtype

macro_subset <- RunUMAP(macro_subset, dims = 1:15)
macro_subset <- RunTSNE(macro_subset, dims = 1:15)

ggsave(file.path(output_dir_recluster, "UMAP_MF_Subtype_Annotated.png"),
       DimPlot(macro_subset, reduction = "umap", group.by = "MF_subtype", label = TRUE, repel = TRUE),
       width = 8, height = 6, dpi = 300)

ggsave(file.path(output_dir_recluster, "tSNE_MF_Subtype_Annotated.png"),
       DimPlot(macro_subset, reduction = "tsne", group.by = "MF_subtype", label = TRUE, repel = TRUE),
       width = 8, height = 6, dpi = 300)

# -------------------- Step 24: MF subtype marker detection --------------------
mf_markers <- data.frame()
for (clust in levels(Idents(macro_subset))) {
  cat("Finding MF markers for cluster", clust, "...\n")
  markers <- FindMarkers(macro_subset,
                         ident.1 = clust,
                         only.pos = TRUE,
                         min.pct = 0.25,
                         logfc.threshold = 0.2)
  if (nrow(markers) > 0) {
    markers$cluster <- clust
    markers$gene <- rownames(markers)
    mf_markers <- rbind(mf_markers, markers)
  }
}

mf_sig <- mf_markers %>%
  filter(p_val_adj < 0.05) %>%
  filter(!grepl("^mt-", gene, ignore.case = TRUE)) %>%
  filter(!grepl("^Rp[sl]", gene, ignore.case = TRUE))

top30_mf <- mf_sig %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30)

write.xlsx(mf_markers, file.path(output_dir_recluster, "DEGs_MFSubtype_Markers.xlsx"))
write.xlsx(top30_mf,  file.path(output_dir_recluster, "Top30_MFSubtype_Markers.xlsx"))

# -------------------- Step 25: Violin plots (MF subtype markers) --------------------
vln_dir3 <- file.path(output_dir_recluster, "Violin_MFSubtype_Markers")
dir.create(vln_dir3, showWarnings = FALSE, recursive = TRUE)

genes_to_plot <- c(
  "Csf1r", "Cd68", "Adgre1", "Fcgr1", "Mrc1",
  "Timd4", "Folr2", "Lyve1",
  "Ccr2", "Plac8", "S100a8",
  "H2-Aa", "H2-Ab1",
  "Cd3d"
)

for (gene in genes_to_plot) {
  p_vln <- VlnPlot(
    macro_subset,
    features = gene,
    group.by = "MF_subtype",
    pt.size = 0
  ) + ggtitle(gene) + theme_classic()
  ggsave(filename = file.path(vln_dir3, paste0("Violin_", gene, ".png")),
         plot = p_vln, width = 5, height = 4, dpi = 300)
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
  if (length(cells) > 0) {
    sampled <- sample(cells, size = min(30, length(cells)))
    cells_use <- c(cells_use, sampled)
  }
}

subset_obj <- subset(macro_subset, cells = cells_use)
DefaultAssay(subset_obj) <- "RNA"
subset_obj <- ScaleData(subset_obj, features = top30_genes, verbose = FALSE)
subset_obj$MF_subtype <- factor(subset_obj$MF_subtype, levels = subtypes)
Idents(subset_obj) <- subset_obj$MF_subtype

p_heat <- DoHeatmap(
  object = subset_obj,
  features = top30_genes
) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 5),
    strip.text.y = element_text(size = 9, face = "bold")
  ) +
  ggtitle("Top 30 Marker Genes per MF Subtype (Lung)")

ggsave(
  filename = file.path(output_dir_recluster, "Heatmap_Top30_MFSubtype_50cells.png"),
  plot = p_heat, width = 13, height = 14, dpi = 600
)

# -------------------- Step 27: Assign MF subtype final labels (manual) --------------------
mf_labels <- c(
  "1" = "TLF⁺ Macrophage",
  "2" = "CCR2⁺ Macrophage"
)

macro_subset$MF_subtype_label <- plyr::mapvalues(
  x = Idents(macro_subset),
  from = names(mf_labels),
  to = mf_labels
)
Idents(macro_subset) <- macro_subset$MF_subtype_label

# Visualize final labels
ggsave(
  filename = file.path(output_dir_recluster, "UMAP_MF_Subtype_ManualLabel.png"),
  plot = DimPlot(macro_subset, group.by = "MF_subtype_label", label = TRUE, repel = TRUE) + theme_classic(),
  width = 8, height = 6, dpi = 300
)
ggsave(
  filename = file.path(output_dir_recluster, "tSNE_MF_Subtype_ManualLabel.png"),
  plot = DimPlot(macro_subset, reduction = "tsne", group.by = "MF_subtype_label", label = TRUE, repel = TRUE) + theme_classic(),
  width = 8, height = 6, dpi = 300
)

# -------------------- Step 28: FeaturePlots for key markers --------------------
key_markers <- c("Timd4", "Lyve1", "Folr2", "Ccr2", "Plac8", "S100a8")

p_feat_umap <- FeaturePlot(
  macro_subset,
  features = key_markers,
  reduction = "umap",
  order = TRUE,
  cols = c("lightgrey", "red")
) & theme_classic()
ggsave(
  filename = file.path(output_dir_recluster, "FeaturePlot_MFSubtypeMarkers_UMAP.png"),
  plot = p_feat_umap, width = 20, height = 16, dpi = 300
)

p_feat_tsne <- FeaturePlot(
  macro_subset,
  features = key_markers,
  reduction = "tsne",
  order = TRUE,
  cols = c("lightgrey", "red")
) & theme_classic()
ggsave(
  filename = file.path(output_dir_recluster, "FeaturePlot_MFSubtypeMarkers_tSNE.png"),
  plot = p_feat_tsne, width = 20, height = 16, dpi = 300
)

# -------------------- Step 29: DEGs & Heatmap by final labels + Save object --------------------
# DEGs by final labels
Idents(macro_subset) <- macro_subset$MF_subtype_label
deg_mf <- FindAllMarkers(
  object = macro_subset,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
top30_mf_final <- deg_mf %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30)

write.xlsx(deg_mf,        file.path(output_dir_recluster, "DEGs_Reassigned_MFSubtype_All.xlsx"))
write.xlsx(top30_mf_final, file.path(output_dir_recluster, "Top30_Reassigned_MFSubtype.xlsx"))

# Heatmap (final labels)
top30_path_final <- file.path(output_dir_recluster, "Top30_Reassigned_MFSubtype.xlsx")
df_top30_final   <- read.xlsx(top30_path_final)
df_top30_final$cluster <- as.character(df_top30_final$cluster)
top30_genes_final <- unique(df_top30_final$gene)

set.seed(111)
subtypes_final <- unique(df_top30_final$cluster)
cells_use2 <- c()
for (stype in subtypes_final) {
  cells <- WhichCells(macro_subset, expression = MF_subtype_label == stype)
  if (length(cells) > 0) {
    sampled <- sample(cells, size = min(30, length(cells)))
    cells_use2 <- c(cells_use2, sampled)
  }
}

subset_obj2 <- subset(macro_subset, cells = cells_use2)
DefaultAssay(subset_obj2) <- "RNA"
subset_obj2 <- ScaleData(subset_obj2, features = top30_genes_final, verbose = FALSE)
subset_obj2$MF_subtype_label <- factor(subset_obj2$MF_subtype_label, levels = subtypes_final)
Idents(subset_obj2) <- subset_obj2$MF_subtype_label

p_heat2 <- DoHeatmap(
  object = subset_obj2,
  features = top30_genes_final
) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 5),
    strip.text.y = element_text(size = 10, face = "bold")
  )
ggsave(
  filename = file.path(output_dir_recluster, "Heatmap_Top30_Reassigned_MFSubtype_50cells.png"),
  plot = p_heat2, width = 13, height = 14, dpi = 600
)

# Save final object
saveRDS(
  macro_subset,
  file = file.path(output_dir_recluster, "MouseLung_MCA_MFSubtype_Final.rds")
)

# ============================ END ================================
