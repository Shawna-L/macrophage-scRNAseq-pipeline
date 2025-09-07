# macrophage-scRNAseq-pipeline
Analysis pipeline for macrophage heterogeneity in fibrotic versus healthy heart, liver and lung using scRNA-seq, combining SingleR reference libraries and marker-based workflows.
This repository provides analysis pipelines to study macrophage heterogeneity in **human** (heart, lung, liver) and **mouse** tissues under healthy and fibrotic conditions using single-cell RNA sequencing (scRNA-seq).  
Workflows integrate **marker-based annotation** and **SingleR reference-based annotation**, ensuring reproducible identification of macrophage subsets across organs and species.

---

## 📂 Contents

- `build_human_singler_reference1.py` – Python script to generate human SingleR reference  
- `build_human_singler_reference2.R` – R script to generate human SingleR reference  
- `dick2022_marker_based_mouse_lung.R` – Marker-based mouse lung pipeline (from Dick et al. 2022)  
- `mca3_mouse_ref_builder.R` – Builds mouse reference using MCA3 data  
- `singler_pipeline_human_*` – Human pipelines: heart, lung, liver (healthy vs fibrotic)  
- `singler_pipeline_mouse_*` – Mouse pipelines: heart, lung, liver  

---

## 🔧 Requirements

- **R (version 4.5.1)** with packages:  
  `Seurat`, `SingleR`, `Matrix`, `scDblFinder`, `org.Hs.eg.db`, `AnnotationDbi`, `ggplot2`, `dplyr`, `openxlsx`

- **Python (version 3.10.18)** with packages:  
  `scanpy`, `pandas`

---

## 🚀 Usage

1. Clone this repository:
   ```bash
   git clone https://github.com/Shawna-L/macrophage-scRNAseq-pipeline.git
   cd macrophage-scRNAseq-pipeline

## Resources

Reference dissociation-associated gene (DAG) lists are provided for reproducibility:

512_Human_Dissociation_genes.csv – Human DAG list (O’Flanagan et al., 2019)

DAG_136_MouseGenes_From_van den Brink et al. 2017_SuppTable5.csv – Mouse DAG list (van den Brink et al., 2017)


## References

Ramachandran P, Dobie R, Wilson-Kanamori JR, Dora EF, Henderson BEP, Luu NT, et al. Resolving the fibrotic niche of human liver cirrhosis at single-cell level. Nature. 2019;575(7783):512–8. doi:10.1038/s41586-019-1631-3.

Vieira Braga FA, Kar G, Berg M, Carpaij OA, Polanski K, Simon LM, et al. A cellular census of human lungs identifies novel cell states in health and in asthma. Nat Med. 2019;25(7):1153–63. doi:10.1038/s41591-019-0468-5.

Tucker NR, Chaffin M, Fleming SJ, Hall AW, Parsons VA, Bedi KC Jr, et al. Transcriptional and cellular diversity of the human heart. Nature. 2020;588(7837):466–72. doi:10.1038/s41586-020-2797-4.

Dick SA, Wong A, Hamidzada H, Nejat S, Nechanitzky R, Vohra S, et al. Three tissue resident macrophage subsets coexist across organs with conserved origins and life cycles. Sci Immunol. 2022;7(74):eabf7777. doi:10.1126/sciimmunol.abf7777.

Han X, Zhou Z, Fei L, Sun H, Wang R, Chen Y, et al. Construction of a human cell landscape at single-cell level. Nucleic Acids Res. 2023;51(2):501–16. doi:10.1093/nar/gkac1014.

Tabula Sapiens Consortium, Jones RC, Karkanias J, Krasnow MA, Pisco AO, Quake SR, et al. The Tabula Sapiens: a multiple-organ single-cell transcriptomic atlas of humans. Science. 2022;376(6594):eabl4896. doi:10.1126/science.abl4896.

Sikkema L, Ramírez-Suástegui C, Strobl DC, Genga RMJ, Marconato M, Pratapa A, et al. An integrated cell atlas of the human lung in health and disease. Nat Med. 2023;29(7):1563–77. doi:10.1038/s41591-023-02327-2.

Gao J, Yang X, Zhang Y, Liu Y, Sun H, Zhou X, et al. Spatially restricted and ontogenically distinct hepatic macrophages are required for tissue repair. bioRxiv. 2025; doi:10.1101/2025.04.16.649149.

Aran D, Looney AP, Liu L, Wu E, Fong V, Hsu A, et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nat Immunol. 2019;20(2):163–72. doi:10.1038/s41590-018-0276-y.

Aran D, Butte AJ. SingleR: annotation of single-cell RNA-seq data. [Internet]. South Wales Bioinformatics; 2023 [cited 2025 Sep 7]. Available from: https://swbioinf.github.io/scRNAseqInR_Doco/singler.html

Heumos L, Schaar AC, Lance C, Ji Y, Tarashansky AJ, Zhang J, et al. Best practices for single-cell analysis across modalities. Nat Rev Genet. 2023;24(9):550–72. doi:10.1038/s41576-023-00586-w.

van den Brink SC, Sage F, Vértesy Á, Spanjaard B, Peterson-Maduro J, Baron CS, et al. Single-cell sequencing reveals dissociation-induced gene expression in tissue subpopulations. Nat Methods. 2017;14(10):935–6. doi:10.1038/nmeth.4437.

O’Flanagan CH, Campbell KR, Zhang AW, Kabeer F, Lim JLP, Biele J, et al. Dissociation of solid tumor tissues with cold active protease for single-cell RNA-seq minimizes conserved collagenase-associated stress responses. Genome Biol. 2019;20(1):210. doi:10.1186/s13059-019-1830-0.
