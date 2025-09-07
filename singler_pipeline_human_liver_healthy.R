# ================================================================
# SingleR-based macrophage annotation pipeline (human liver – healthy)
# ================================================================
# Author: Jingshuo Li
# Version: 1.0
# Data source: GEO GSE136103 (healthy human liver; samples healthy1–healthy5)
# Description: End-to-end scRNA-seq workflow for HUMAN LIVER (HEALTHY)
# with 512-gene DAG-QC, scDblFinder doublet removal, clustering,
# SingleR (human liver reference), macrophage reclustering,
# markers/heatmaps/violin/featureplots, and exports.
# ================================================================

# -------------------- Step 0: Parameters (EDIT HERE) --------------------
set.seed(111)

organ     <- "liver"
condition <- "healthy"

# relative paths for GitHub
base_data_dir <- file.path("data", "human", organ, condition)
sample_ids    <- c("healthy1","healthy2","healthy3","healthy4","healthy5")
sample_paths  <- setNames(file.path(base_data_dir, sample_ids), sample_ids)

dag_path <- file.path("data", "DAG_512_HumanGenes.csv")  # 512-gene human DAG list (col1=gene symbol)
ref_path <- file.path("refs", "human", "liver_ref.rds")  # your liver reference RDS (e.g., IFALD-based)

results_dir <- file.path("results", "human", organ, condition, "singler_pipeline")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- Step 1: Load packages --------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(SingleR)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(openxlsx)
  library(cluster)
  library(patchwork)
  library(plyr)
})

# -------------------- Step 2: Define sample paths --------------------
stopifnot(file.exists(dag_path), file.exists(ref_path))
stopifnot(all(dir.exists(unname(sample_paths))))

# -------------------- Step 3: Read & create individual Seurat objects and merge --------------------
seurat_list <- lapply(names(sample_paths), function(sample) {
  m  <- Read10X(data.dir = sample_paths[[sample]])
  so <- CreateSeuratObject(counts = m, project = sample, min.cells = 3, min.features = 200)
  so$sample <- sample
  so$group  <- "Healthy"
  so
})
names(seurat_list) <- names(sample_paths)

liver_healthy <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "Liver_Healthy"
)

# -------------------- Step 4: Join layers (Seurat v5 compatibility) --------------------
liver_healthy <- JoinLayers(liver_healthy)

# -------------------- Step 5: Load human DAG gene list --------------------
dag_df      <- read.csv(dag_path, stringsAsFactors = FALSE)
dag_genes   <- unique(dag_df[[1]])
overlap_dag <- intersect(rownames(liver_healthy), dag_genes)

# -------------------- Step 6: Add QC metrics --------------------
# human mitochondrial genes: "^MT-"
liver_healthy[["percent.mito"]] <- PercentageFeatureSet(liver_healthy, pattern = "^MT-")
liver_healthy[["percent.dag"]]  <- PercentageFeatureSet(liver_healthy, features = overlap_dag)

print(summary(liver_healthy$percent.mito))
print(summary(liver_healthy$percent.dag))
print(summary(liver_healthy$nFeature_RNA))

# -------------------- Step 7: QC histograms --------------------
hist_dir <- file.path(results_dir, "QC_Histograms")
dir.create(hist_dir, showWarnings = FALSE, recursive = TRUE)

png(file.path(hist_dir, "Histogram_percent.mito.png"), width = 800, height = 600)
hist(liver_healthy$percent.mito, breaks = 50, col = "lightblue", main = "Mito %", xlab = "percent.mito")
dev.off()

png(file.path(hist_dir, "Histogram_percent.dag.png"), width = 800, height = 600)
hist(liver_healthy$percent.dag, breaks = 50, col = "pink", main = "DAG %", xlab = "percent.dag")
dev.off()

# -------------------- Step 8: QC scatter plots --------------------
qc_plot_dir <- file.path(results_dir, "QC_Plots")
dir.create(qc_plot_dir, showWarnings = FALSE, recursive = TRUE)

p1 <- FeatureScatter(liver_healthy, feature1 = "nCount_RNA", feature2 = "percent.mito")
p2 <- FeatureScatter(liver_healthy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p3 <- FeatureScatter(liver_healthy, feature1 = "nCount_RNA", feature2 = "percent.dag")
p4 <- FeatureScatter(liver_healthy, feature1 = "percent.dag", feature2 = "percent.mito")

ggsave(file.path(qc_plot_dir, "Scatter_nCount_vs_percent.mito.png"), p1, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(qc_plot_dir, "Scatter_nCount_vs_nFeature.png"),  p2, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(qc_plot_dir, "Scatter_nCount_vs_percent.dag.png"), p3, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(qc_plot_dir, "Scatter_percent.dag_vs_percent.mito.png"), p4, width = 8, height = 8, dpi = 300, bg = "white")
rm(p1,p2,p3,p4); gc()

# -------------------- Step 9: Final QC filtering + Doublet removal (scDblFinder) --------------------
liver_healthy_filtered <- subset(
  liver_healthy,
  subset = nFeature_RNA > 300 &
    nFeature_RNA < 8000 &
    percent.mito   < 15 &
    percent.dag    < 8
)
message("Raw cells: ", ncol(liver_healthy))
message("After QC filtering: ", ncol(liver_healthy_filtered))
message("Retention ratio: ", round(ncol(liver_healthy_filtered) / ncol(liver_healthy), 3))

# scDblFinder (per-sample doublet calling)
sce <- as.SingleCellExperiment(liver_healthy_filtered)
sce <- scDblFinder(sce, samples = colData(sce)$sample)
liver_healthy_filtered$doublet_class <- colData(sce)$scDblFinder.class
liver_healthy_filtered$doublet_score <- colData(sce)$scDblFinder.score

dbl_dir <- file.path(results_dir, "Doublets")
dir.create(dbl_dir, showWarnings = FALSE, recursive = TRUE)
dbl_tab <- as.data.frame(table(liver_healthy_filtered$doublet_class))
dbl_tab$Proportion <- round(dbl_tab$Freq / sum(dbl_tab$Freq), 4)
write.csv(dbl_tab, file = file.path(dbl_dir, "doublet_summary.csv"), row.names = FALSE)

liver_healthy_clean <- subset(liver_healthy_filtered, subset = doublet_class == "singlet")
saveRDS(liver_healthy_clean, file = file.path(results_dir, "liver_healthy_filtered_singlets.rds"))

# -------------------- Step 10: Normalize data & JoinLayers --------------------
liver_healthy_clean <- NormalizeData(liver_healthy_clean)
liver_healthy_clean <- JoinLayers(liver_healthy_clean)

# -------------------- Step 11: FindVariableFeatures (mean.var.plot) --------------------
liver_healthy_meanvar <- FindVariableFeatures(
  liver_healthy_clean,
  selection.method = "mean.var.plot",
  mean.cutoff = c(0.01, 5),
  dispersion.cutoff = c(0.5, Inf),
  nfeatures = 5000
)
message("Variable features: ", length(VariableFeatures(liver_healthy_meanvar)))

output_dir_meanvar <- file.path(results_dir, "Method_meanvar")
dir.create(output_dir_meanvar, recursive = TRUE, showWarnings = FALSE)

# -------------------- Step 12: Scale + PCA --------------------
liver_healthy_meanvar <- ScaleData(
  liver_healthy_meanvar,
  vars.to.regress = c("nCount_RNA", "percent.mito")
)
liver_healthy_meanvar <- RunPCA(
  liver_healthy_meanvar,
  features = VariableFeatures(liver_healthy_meanvar)
)

# -------------------- Step 13: Elbow plot --------------------
p_elbow <- ElbowPlot(liver_healthy_meanvar) + theme_classic()
ggsave(filename = file.path(output_dir_meanvar, "ElbowPlot_meanvar.png"),
       plot = p_elbow, width = 6, height = 5, dpi = 300)
rm(p_elbow)

# -------------------- Step 14: Neighbors + Clusters --------------------
liver_healthy_meanvar <- FindNeighbors(liver_healthy_meanvar, dims = 1:15)
liver_healthy_meanvar <- FindClusters(liver_healthy_meanvar, resolution = 1.0)

Idents(liver_healthy_meanvar) <- liver_healthy_meanvar$seurat_clusters
lvl <- Idents(liver_healthy_meanvar)
levels(lvl) <- as.character(as.numeric(levels(lvl)) + 1)
liver_healthy_meanvar$Cluster1based <- lvl
Idents(liver_healthy_meanvar) <- liver_healthy_meanvar$Cluster1based

# -------------------- Step 15: UMAP --------------------
liver_healthy_meanvar <- RunUMAP(liver_healthy_meanvar, reduction = "pca", dims = 1:15)

# -------------------- Step 16: tSNE --------------------
liver_healthy_meanvar <- RunTSNE(liver_healthy_meanvar, reduction = "pca", dims = 1:15)

# -------------------- Step 17: Save UMAP & tSNE (cluster labels) --------------------
p_umap <- DimPlot(liver_healthy_meanvar, reduction = "umap", label = TRUE) + theme_classic()
ggsave(file.path(output_dir_meanvar, "UMAP_FirstClusters_meanvar.png"),
       p_umap, width = 8, height = 6, dpi = 300)
rm(p_umap)

p_tsne <- DimPlot(liver_healthy_meanvar, reduction = "tsne", label = TRUE) + theme_classic()
ggsave(file.path(output_dir_meanvar, "TSNE_FirstClusters_meanvar.png"),
       p_tsne, width = 8, height = 6, dpi = 300)
rm(p_tsne)

# -------------------- Step 18: Silhouette score --------------------
emb_pca <- Embeddings(liver_healthy_meanvar, reduction = "pca")[, 1:15]
cluster_labels <- as.factor(liver_healthy_meanvar$Cluster1based)
dists <- dist(emb_pca)
sil   <- silhouette(as.integer(cluster_labels), dists)

sil_df <- as.data.frame(sil)
sil_df$Cluster <- cluster_labels[as.numeric(rownames(sil_df))]
sil_summary <- sil_df %>%
  group_by(Cluster) %>%
  summarise(mean_sil = mean(sil_width)) %>%
  arrange(desc(mean_sil))
print(sil_summary, n = 30)

sil_dir <- file.path(results_dir, "SilhouettePlot")
dir.create(sil_dir, showWarnings = FALSE, recursive = TRUE)
p_sil <- ggplot(sil_df, aes(x = Cluster, y = sil_width, fill = Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.5) +
  theme_classic() +
  labs(y = "Silhouette width", title = "Silhouette Score per Cluster")
ggsave(filename = file.path(sil_dir, "Silhouette_PerCluster.png"),
       plot = p_sil, width = 8, height = 6, dpi = 300)
rm(p_sil)

# ============================ SingleR ==================================
# -------------------- Step 19: SingleR annotation (human liver reference) --------------------
output_dir_singler <- file.path(results_dir, "AutoAnnotation_SingleR")
dir.create(output_dir_singler, showWarnings = FALSE, recursive = TRUE)

sce_singler <- as.SingleCellExperiment(liver_healthy_meanvar)
ref_liver   <- readRDS(ref_path)

# Optional parallel (commented for portability)
# library(BiocParallel); param <- MulticoreParam(workers = 4)
# pred_liver <- SingleR(test = sce_singler, ref = ref_liver, labels = ref_liver$label.main, BPPARAM = param)

pred_liver <- SingleR(test = sce_singler, ref = ref_liver, labels = ref_liver$label.main)
liver_healthy_meanvar$HCA_liver_label <- pred_liver$labels

# -------------------- Step 20: Save SingleR outputs & plots --------------------
write.csv(as.data.frame(pred_liver),
          file.path(output_dir_singler, "SingleR_HCA_liver_Result.csv"), row.names = TRUE)

Idents(liver_healthy_meanvar) <- liver_healthy_meanvar$HCA_liver_label
p_umap2 <- DimPlot(liver_healthy_meanvar, reduction = "umap", label = FALSE, repel = TRUE) +
  ggtitle("SingleR UMAP – Human liver") + theme_classic() + coord_fixed()
p_tsne2 <- DimPlot(liver_healthy_meanvar, reduction = "tsne", label = FALSE, repel = TRUE) +
  ggtitle("SingleR tSNE – Human liver") + theme_classic() + coord_fixed()

ggsave(file.path(output_dir_singler, "UMAP_SingleR_HCA_liver.png"),
       p_umap2, width = 10, height = 7, dpi = 300)
ggsave(file.path(output_dir_singler, "tSNE_SingleR_HCA_liver.png"),
       p_tsne2, width = 10, height = 7, dpi = 300)
rm(p_umap2, p_tsne2)

saveRDS(liver_healthy_meanvar, file = file.path(output_dir_singler, "liver_healthy_annotated_HCA_liver.rds"))
write.csv(as.data.frame(table(liver_healthy_meanvar$HCA_liver_label)),
          file = file.path(output_dir_singler, "SingleR_Label_Summary.csv"), row.names = FALSE)

# (Optional) Myeloid marker violin check across selected SingleR labels
vln_dir_macrocheck <- file.path(output_dir_singler, "Violin_MacrophageMarkers")
dir.create(vln_dir_macrocheck, recursive = TRUE, showWarnings = FALSE)
seurat_vln_subset <- subset(
  liver_healthy_meanvar,
  subset = HCA_liver_label %in% c("macrophage","Kupffer cell","mononuclear phagocyte","cycling myeloid cell")
)
if (ncol(seurat_vln_subset) > 0) {
  seurat_vln_subset$HCA_liver_label <- droplevels(factor(seurat_vln_subset$HCA_liver_label))
  Idents(seurat_vln_subset) <- seurat_vln_subset$HCA_liver_label
  marker_genes <- c("CSF1R","CD68","FCGR1A","CD163","MRC1","MARCO","VSIG4","CCR2","IL1B","PLAC8",
                    "MKI67","TOP2A","NUSAP1","HLA-DRA","HLA-DRB1","S100A8","S100A9",
                    "CD14","FCN1","CD1C","CLEC9A","LYZ")
  for (gene in marker_genes) {
    if (gene %in% rownames(seurat_vln_subset)) {
      p <- VlnPlot(seurat_vln_subset, features = gene, group.by = "HCA_liver_label", pt.size = 0) +
        ggtitle(gene) + theme_classic(base_size = 12) +
        theme(legend.position = "none", axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks.x = element_blank())
      ggsave(file.path(vln_dir_macrocheck, paste0("Violin_", gene, ".png")), p, width = 5, height = 5, dpi = 300)
    } else message("Gene ", gene, " not found; skipped.")
  }
}

# -------------------- Step 21: DEGs per SingleR label --------------------
Idents(liver_healthy_meanvar) <- liver_healthy_meanvar$HCA_liver_label
deg_results <- FindAllMarkers(
  object = liver_healthy_meanvar,
  only.pos = TRUE,
  min.pct = 0.3,
  logfc.threshold = 0.25
)
write.csv(deg_results, file = file.path(output_dir_singler, "All_DEGs_SingleR_Celltypes.csv"), row.names = FALSE)

top30_deg <- deg_results %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 30)
write.csv(top30_deg, file = file.path(output_dir_singler, "Top30_DEGs_SingleR_Celltypes.csv"), row.names = FALSE)

# ============================ Macrophage Focus =============================
# -------------------- Step 22: Macrophage reclustering --------------------
output_dir_recluster <- file.path(results_dir, "recluster")
dir.create(output_dir_recluster, recursive = TRUE, showWarnings = FALSE)

Idents(liver_healthy_meanvar) <- liver_healthy_meanvar$HCA_liver_label
macro_cells <- colnames(liver_healthy_meanvar)[grepl("macroph|kupffer|mononuclear phagocyte", 
                                                     liver_healthy_meanvar$HCA_liver_label, ignore.case = TRUE)]
liver_healthy_macro_subset <- subset(liver_healthy_meanvar, cells = macro_cells)
saveRDS(liver_healthy_macro_subset, file = file.path(output_dir_recluster, "HumanLiver_HCA_Macrophages.rds"))

liver_healthy_macro_subset <- FindVariableFeatures(
  liver_healthy_macro_subset,
  selection.method = "mean.var.plot",
  mean.cutoff = c(0.1, 5),
  dispersion.cutoff = c(0.5, Inf),
  nfeatures = 5000
)
p_elbow_macro <- ElbowPlot(liver_healthy_macro_subset) + theme_classic()
ggsave(file.path(output_dir_recluster, "ElbowPlot_macro_subset.png"),
       p_elbow_macro, width = 6, height = 5, dpi = 300)
rm(p_elbow_macro); gc()

liver_healthy_macro_subset <- ScaleData(liver_healthy_macro_subset)
liver_healthy_macro_subset <- RunPCA(liver_healthy_macro_subset, features = VariableFeatures(liver_healthy_macro_subset))
liver_healthy_macro_subset <- FindNeighbors(liver_healthy_macro_subset, dims = 1:15)
liver_healthy_macro_subset <- FindClusters(liver_healthy_macro_subset, resolution = 0.2)

Idents(liver_healthy_macro_subset) <- liver_healthy_macro_subset$seurat_clusters
levels(Idents(liver_healthy_macro_subset)) <- as.character(as.numeric(levels(Idents(liver_healthy_macro_subset))) + 1)
liver_healthy_macro_subset$MF_subtype <- Idents(liver_healthy_macro_subset)
Idents(liver_healthy_macro_subset) <- liver_healthy_macro_subset$MF_subtype

liver_healthy_macro_subset <- RunUMAP(liver_healthy_macro_subset, dims = 1:15)
liver_healthy_macro_subset <- RunTSNE(liver_healthy_macro_subset, dims = 1:15)

ggsave(file.path(output_dir_recluster, "UMAP_MF_Subtype_Annotated.png"),
       DimPlot(liver_healthy_macro_subset, reduction = "umap", group.by = "MF_subtype", label = TRUE, repel = TRUE),
       width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir_recluster, "tSNE_MF_Subtype_Annotated.png"),
       DimPlot(liver_healthy_macro_subset, reduction = "tsne", group.by = "MF_subtype", label = TRUE, repel = TRUE),
       width = 8, height = 6, dpi = 300)

# -------------------- Step 23: MF subtype marker detection --------------------
mf_markers <- data.frame()
for (clust in levels(Idents(liver_healthy_macro_subset))) {
  cat("Finding MF markers for cluster", clust, "...\n")
  markers <- FindMarkers(liver_healthy_macro_subset,
                         ident.1 = clust,
                         only.pos = TRUE,
                         min.pct = 0.3,
                         logfc.threshold = 0.25)
  if (nrow(markers) > 0) {
    markers$cluster <- clust
    markers$gene    <- rownames(markers)
    mf_markers      <- rbind(mf_markers, markers)
  }
}
mf_sig <- mf_markers %>%
  filter(p_val_adj < 0.05) %>%
  filter(!grepl("^MT-", gene, ignore.case = TRUE)) %>%
  filter(!grepl("^RP[SL]", gene, ignore.case = TRUE))
top30_mf <- mf_sig %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 30)

write.xlsx(mf_markers, file.path(output_dir_recluster, "DEGs_MFSubtype_Markers.xlsx"))
write.xlsx(top30_mf,  file.path(output_dir_recluster, "Top30_MFSubtype_Markers.xlsx"))

# -------------------- Step 24: Violin plots (MF subtype markers) --------------------
vln_dir <- file.path(output_dir_recluster, "Violin_MFSubtype_Markers")
dir.create(vln_dir, showWarnings = FALSE, recursive = TRUE)

genes_to_plot <- c(
  # General macrophage markers
  "CD68","CSF1R","AIF1","LILRB4",
  # TLF⁺ resident-like (LYVE1/FOLR2/TIMD4)
  "LYVE1","FOLR2","FABP5","GPNMB","RNASE1","TIMD4",
  # Kupffer-specific
  "CD5L","VCAM1","MARCO","VSIG4","CD163",
  # Inflammatory
  "IL1B","S100A8","CCL4","CCL3L3","CXCL8","CCR2","NFKBID","VCAN","NEAT1","BID","FCN1",
  # Proliferation
  "MKI67","TOP2A","NUSAP1",
  # LAM / lipid & remodeling
  "TREM2","FABP4","MMP9","APOE","CD9","SPP1","PLP2",
  # MHC-II / IFN
  "HLA-DRA","CD74","B2M","IFI6","STAT1"
)
for (gene in genes_to_plot) {
  p_vln <- VlnPlot(liver_healthy_macro_subset, features = gene, group.by = "MF_subtype", pt.size = 0) +
    ggtitle(gene) + theme_classic()
  ggsave(file.path(vln_dir, paste0("Violin_", gene, ".png")), p_vln, width = 5, height = 4, dpi = 300)
}
rm(p_vln)

# -------------------- Step 25: Heatmap (Top30 per MF subtype) --------------------
top30_path <- file.path(output_dir_recluster, "Top30_MFSubtype_Markers.xlsx")
df_top30   <- read.xlsx(top30_path)
df_top30$cluster <- as.character(df_top30$cluster)
top30_genes <- unique(df_top30$gene)

set.seed(111)
subtypes <- unique(df_top30$cluster)
cells_use <- c()
for (stype in subtypes) {
  cells <- WhichCells(liver_healthy_macro_subset, idents = stype)
  if (length(cells) > 0) cells_use <- c(cells_use, sample(cells, size = min(50, length(cells))))
}

subset_obj <- subset(liver_healthy_macro_subset, cells = cells_use)
DefaultAssay(subset_obj) <- "RNA"
subset_obj <- ScaleData(subset_obj, features = top30_genes, verbose = FALSE)
subset_obj$MF_subtype <- factor(subset_obj$MF_subtype, levels = subtypes)
Idents(subset_obj) <- subset_obj$MF_subtype

p_heat <- DoHeatmap(subset_obj, features = top30_genes) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6),
        strip.text.y = element_text(size = 9, face = "bold")) +
  ggtitle("Top 30 Marker Genes per MF Subtype (Liver)")
ggsave(file.path(output_dir_recluster, "Heatmap_Top30_MFSubtype_50cells.png"),
       p_heat, width = 13, height = 14, dpi = 600)
rm(p_heat)

# -------------------- Step 26: Assign MF subtype final labels (manual) --------------------
mf_colors <- c(
  "TLF-A Macrophage (Homeostatic)"                 = "#3B4992",
  "Inflammatory Macrophage"                        = "#EE7733",
  "Stress-Responsive Macrophage"                   = "#009988",
  "Plasma-like Macrophage (Ig⁺ / Stress-Adaptive)" = "#F2CC0C",
  "TLF⁺ Macrophage (Kupffer-like)"                 = "#CC79A7",
  "TLF-B Macrophage (Inflammatory-Primed)"         = "#D62828",
  "TLF-C Macrophage (Proliferating)"               = "#00B0F6",
  "ECM Remodeling & Lipid Handling Macrophage"     = "#61D04F",
  "Monocyte-Derived Macrophage (IFN-Responsive)"   = "#DC267F",
  "Metabolically-Active Macrophage"                = "#8B5E3C",
  "Antigen-Presenting Macrophage"                  = "#6A5ACD"
)

mf_labels <- c(
  "1" = "Inflammatory Macrophage",
  "2" = "Stress-Responsive Macrophage",
  "3" = "Antigen-Presenting Macrophage",
  "4" = "TLF⁺ Macrophage (Kupffer-like)"
)

liver_healthy_macro_subset$MF_subtype_label <- plyr::mapvalues(
  x = Idents(liver_healthy_macro_subset),
  from = names(mf_labels),
  to   = mf_labels
)
Idents(liver_healthy_macro_subset) <- liver_healthy_macro_subset$MF_subtype_label

ggsave(file.path(output_dir_recluster, "UMAP_MF_Subtype_ManualLabel.png"),
       DimPlot(liver_healthy_macro_subset, reduction = "umap",
               group.by = "MF_subtype_label", cols = mf_colors, label = FALSE, repel = TRUE) +
         ggtitle("UMAP of MF Subtypes") + theme_classic() + coord_fixed(),
       width = 14, height = 8, dpi = 300)

ggsave(file.path(output_dir_recluster, "tSNE_MF_Subtype_ManualLabel.png"),
       DimPlot(liver_healthy_macro_subset, reduction = "tsne",
               group.by = "MF_subtype_label", cols = mf_colors, label = FALSE, repel = TRUE) +
         ggtitle("tSNE of MF Subtypes") + theme_classic() + coord_fixed(),
       width = 14, height = 8, dpi = 300)

# -------------------- Step 27: FeaturePlots for key markers --------------------
key_markers <- c(
  # TLF⁺ (Kupffer-like)
  "FOLR2","LYVE1","CD163","VSIG4",
  # Inflammatory
  "IL1B","S100A8","NFKBID",
  # Stress-Responsive
  "NEAT1","STAT1","BID",
  # Antigen-Presenting
  "HLA-DOA","HLA-DPA1","HLA-DPB1"
)
p_feat_umap <- FeaturePlot(liver_healthy_macro_subset, features = key_markers,
                           reduction = "umap", order = TRUE, cols = c("lightgrey", "red")) & theme_classic()
ggsave(file.path(output_dir_recluster, "FeaturePlot_MFSubtypeMarkers_UMAP.png"),
       p_feat_umap, width = 30, height = 20, dpi = 300)

p_feat_tsne <- FeaturePlot(liver_healthy_macro_subset, features = key_markers,
                           reduction = "tsne", order = TRUE, cols = c("lightgrey", "red")) & theme_classic()
ggsave(file.path(output_dir_recluster, "FeaturePlot_MFSubtypeMarkers_tSNE.png"),
       p_feat_tsne, width = 30, height = 20, dpi = 300)

# -------------------- Step 28: DEGs & Heatmap by final labels --------------------
Idents(liver_healthy_macro_subset) <- liver_healthy_macro_subset$MF_subtype_label
deg_mf_final <- FindAllMarkers(
  object = liver_healthy_macro_subset,
  only.pos = TRUE,
  min.pct = 0.3,
  logfc.threshold = 0.25
)
top30_mf_final <- deg_mf_final %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 30)

write.xlsx(deg_mf_final,  file.path(output_dir_recluster, "DEGs_Reassigned_MFSubtype_All.xlsx"))
write.xlsx(top30_mf_final, file.path(output_dir_recluster, "Top30_Reassigned_MFSubtype.xlsx"))

df_top30_final   <- read.xlsx(file.path(output_dir_recluster, "Top30_Reassigned_MFSubtype.xlsx"))
df_top30_final$cluster <- as.character(df_top30_final$cluster)
top30_genes_final <- unique(df_top30_final$gene)

set.seed(111)
subtypes_final <- unique(df_top30_final$cluster)
cells_use2 <- c()
for (stype in subtypes_final) {
  cells <- WhichCells(liver_healthy_macro_subset, idents = stype)
  if (length(cells) > 0) cells_use2 <- c(cells_use2, sample(cells, size = min(50, length(cells))))
}

subset_obj2 <- subset(liver_healthy_macro_subset, cells = cells_use2)
DefaultAssay(subset_obj2) <- "RNA"
subset_obj2 <- ScaleData(subset_obj2, features = top30_genes_final, verbose = FALSE)
subset_obj2$MF_subtype_label <- factor(subset_obj2$MF_subtype_label, levels = subtypes_final)
Idents(subset_obj2) <- subset_obj2$MF_subtype_label

p_heat2 <- DoHeatmap(subset_obj2, features = top30_genes_final) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6),
        strip.text.y = element_text(size = 5, face = "bold"))
ggsave(file.path(output_dir_recluster, "Heatmap_Top30_Reassigned_MFSubtype_50cells.png"),
       p_heat2, width = 14, height = 15, dpi = 600)
rm(p_heat2)

# -------------------- Step 29: Save final Seurat object --------------------
saveRDS(liver_healthy_macro_subset,
        file = file.path(output_dir_recluster, "HumanLiver_HCA_MFSubtype_Final.rds"))

# -------------------- Step 30: Barplot of MF subtype counts (UMAP-matched colors) --------------------
mf_counts <- as.data.frame(table(liver_healthy_macro_subset$MF_subtype_label))
colnames(mf_counts) <- c("MF_Subtype", "Cell_Count")
mf_counts$Percentage <- round(mf_counts$Cell_Count / sum(mf_counts$Cell_Count) * 100, 1)
mf_counts <- mf_counts[order(-mf_counts$Cell_Count), ]
mf_counts$MF_Subtype <- factor(mf_counts$MF_Subtype, levels = mf_counts$MF_Subtype)

p_bar <- ggplot(mf_counts, aes(x = MF_Subtype, y = Cell_Count, fill = MF_Subtype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Cell_Count, "\n", Percentage, "%")),
            vjust = -0.5, size = 4.2) +
  scale_fill_manual(values = mf_colors[names(mf_colors) %in% mf_counts$MF_Subtype]) +
  expand_limits(y = max(mf_counts$Cell_Count) * 1.15) +
  theme_classic() +
  xlab(NULL) + ylab("Cell Count") +
  ggtitle("Cell Number and Percentage per Macrophage Subtype (Liver Healthy)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text  = element_text(size = 10))
ggsave(file.path(output_dir_recluster, "Barplot_MFSubtype_LiverHealthy_UMAPcolors.png"),
       p_bar, width = 9, height = 6, dpi = 300)
rm(p_bar)
