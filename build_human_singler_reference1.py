#!/usr/bin/env python3.10.18
# -*- coding: utf-8 -*-

# ================================================================
# Build a SingleR-ready reference (genes × cells + labels) from .h5ad
# ================================================================
# Dataset sources (CellxGene, for citation):
# - Heart (Tabula Sapiens – heart subset):
#   https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5
# - Lung (An Integrated Cell Atlas of the Human Lung in Health and Disease):
#   https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293
# - Liver (A Single-Cell Atlas of Human Pediatric Liver Reveals Age-Related Hepatic Gene Signatures):
#   https://cellxgene.cziscience.com/collections/ff69f0ee-fef6-4895-9f48-6c64a68c8289
# ================================================================

# -------------------- Step 0: Parameters (EDIT HERE) --------------------
ORGAN       = "lung"   # "lung" | "heart" | "liver"
INPUT_H5AD  = "data/human/lung/reference/hlca_core.h5ad"             # <- path to .h5ad
OUTPUT_DIR  = "refs/human/lung/HLCA_core"                            # <- output folder
LABEL_MAIN  = "cell_type"                                            # <- obs column for main labels
LABEL_FINE  = "cell_subtype"                                         # <- obs column for fine labels (fallback to main)
USE_RAW     = True                                                   # use adata.raw.X if available
WRITE_CSV   = True                                                   # write genes×cells csv.gz (may be large)
WRITE_MTX   = False                                                  # write 10X mtx (set True for huge datasets)

# -------------------- Step 1: Imports --------------------
import os, gzip
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

# -------------------- Step 2: Load AnnData --------------------
os.makedirs(OUTPUT_DIR, exist_ok=True)
adata = ad.read_h5ad(INPUT_H5AD)

X   = adata.raw.X if (USE_RAW and adata.raw is not None) else adata.X
VAR = adata.raw.var if (USE_RAW and adata.raw is not None) else adata.var
CELLS = adata.obs_names.rename("cell")

# -------------------- Step 3: Pick gene symbols & de-duplicate --------------------
CAND_GENE_COLS = ["gene_symbol","gene_symbols","GeneSymbol","SYMBOL","symbol","hgnc_symbol"]
if any(c in VAR.columns for c in CAND_GENE_COLS):
    for c in CAND_GENE_COLS:
        if c in VAR.columns:
            genes = VAR[c].astype(str)
            if genes.notna().sum() and (genes != "").any():
                break
else:
    genes = pd.Index(VAR.index.astype(str), name="gene")

# make unique: APPEND .2/.3/... to duplicates
if pd.Index(genes).has_duplicates:
    counts = {}
    uniq = []
    for g in genes:
        if g not in counts:
            counts[g] = 1; uniq.append(g)
        else:
            counts[g] += 1; uniq.append(f"{g}.{counts[g]}")
    genes = pd.Index(uniq, name="gene")
else:
    genes = pd.Index(genes, name="gene")

# -------------------- Step 4: Build labels (with simple fallback) --------------------
obs = adata.obs
label_main = obs.get(LABEL_MAIN, pd.Series(["unknown"]*obs.shape[0], index=obs.index)).astype(str)
label_fine = obs.get(LABEL_FINE, label_main).astype(str)
labels = pd.DataFrame({"label.main": label_main, "label.fine": label_fine}, index=CELLS)
labels.index.name = "cell"
labels.to_csv(os.path.join(OUTPUT_DIR, "labels.csv"))

# -------------------- Step 5: Write expression (genes × cells) and/or 10X mtx --------------------
if WRITE_CSV:
    # NOTE: CSV requires dense; for huge data, set WRITE_CSV=False and WRITE_MTX=True
    Xdense = X.toarray() if sparse.issparse(X) else np.asarray(X)
    expr_df = pd.DataFrame(Xdense, index=CELLS, columns=genes, copy=False).T  # genes × cells
    expr_df.to_csv(os.path.join(OUTPUT_DIR, "expression.csv.gz"), compression="gzip")

if WRITE_MTX:
    # matrix.mtx.gz + features.tsv.gz + barcodes.tsv.gz
    Xmtx = sparse.csc_matrix(X) if not sparse.issparse(X) else X.tocsc()
    with gzip.open(os.path.join(OUTPUT_DIR,"matrix.mtx.gz"), "wb") as f:
        mmwrite(f, Xmtx.tocoo())
    with gzip.open(os.path.join(OUTPUT_DIR,"features.tsv.gz"), "wt", encoding="utf-8") as f:
        for g in genes: f.write(f"{g}\t{g}\tGene Expression\n")
    with gzip.open(os.path.join(OUTPUT_DIR,"barcodes.tsv.gz"), "wt", encoding="utf-8") as f:
        for c in CELLS: f.write(f"{c}\n")

# -------------------- Step 6: Log --------------------
print("[OK] SingleR reference exported")
print(f" Organ: {ORGAN}")
print(f" Input: {INPUT_H5AD}")
print(f" Output dir: {OUTPUT_DIR}")
print(f" Cells: {len(CELLS):,} | Genes: {len(genes):,}")
print(f" Wrote CSV: {WRITE_CSV} | Wrote 10X mtx: {WRITE_MTX}")
