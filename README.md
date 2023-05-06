# Single-cell-RNA-Sequencing_SCLC_chemotherapy_Pipeline

This pipeline was used for visualization of clinical information of enrolled small cell lung cancer (SCLC) patients, generation of gene expression matrices, scRNA-seq profiling of cellular heterogeneity in SCLC patients, transcriptional characters and alterations of NE (neuroendocrine) cells in SCLC patients following chemotherapy, as well as chemotherapy-induced transcriptional changes in endothelial cells, fibroblasts, B cells and T cells in SCLC patients based on the single-cell RNA sequencing data.



This pipeline was described in the publication: XXX



The raw single-cell RNA-seq data generated in this study is available in the National Genomics Data Center (NGDC) Genome Sequence Archive (GSA) database (https://bigd.big.ac.cn/gsa485 human/) under accession number HRA003150.




## Structure of this repository

Directory: src



Part1_clinicl_information.R

Clinical information of enrolled patients.



Part1_count.sh

Generation of gene expression matrices. 



Part1_cell-assign.R

UMAP representation of all cells analyzed in tissues from treated and treatment-naïve patients.  The expression of marker genes in the indicated cell types. Fractions of immune and non-immune cells from libraries prepared from samples depleted or not of CD45+ cells. The proportion of cell types in tumor and paired normal tissue in treated and untreated patients.



Part2_epithelial-cell.R

UMAP representation of subclustered epithelial cells.



Part2_NE-cell_treated-vs-untreated.R

Heatmap of genes differentially expressed between NE cells and other types of epithelial cells. GO and KEGG pathway enrichment of genes in NE cells differentially expressed between treated and treatment-naïve patients. Violin plots showing differential expression of genes associated with antigen processing and presentation as well as cellular senescence between treated and treatment-naïve patients.  Violin plots showing differences in gene set activity scores for genes associated with DNA repair and therapy targets between treated and treatment-naïve patients. 



Part3_endothelial-cell_treated-vs-untreated.R

UMAP visualization of subclustered endothelial cells. Heatmap showing correlation of endothelial cell types with expression of their established markers.  Heatmap of genes in stalk-like endothelial cells differentially expressed between treated and treatment-naïve patients.  GO and KEGG pathway enrichment analyses of genes in stalk-like endothelial cells differentially expressed between treated and treatment-naïve patients. Violin plots showing the difference in expression of selected genes associated with angiogenesis and fluid shear stress in stalk-like endothelial cells between treated and treatment-naïve patients. Violin plots showing the difference in expression of genes involved in VEGFR signaling and transendothelial migration of lymphatic endothelial cells between treated and treatment-naïve patients.



Part4_fibroblast_treated-vs-untreated.R

UMAP plot of subclustered fibroblasts. Heatmap showing expression level of classic gene markers of fibroblast subtypes. The GSVA score for most significant pathways enriched in plasma cells from treated patients. Unsupervised clustering of GSVA scores of fibroblast subtypes in treated and treatment-naïve patients. GO and KEGG pathway enrichment analyses of genes in fibroblasts differentially expressed between treated and treatment-naïve patients. Ridgeline plot of collagen gene expression in fibroblasts.


Part5_B-cell_treated-vs-untreated.R

Immune cell type composition of each patient. UMAP plot of subclustered B cells in patients. Heatmap showing expression of canonical markers of B cell subtypes. Expression heatmap of top 100 genes in plasma cells differentially expressed between treated and treatment-naïve patients. GO enrichment of genes in plasma cells differentially expressed between treated and treatment-naïve patients. Bar chart illustrating the GSVA score for most significant pathways enriched in treatment patients compared to treatment-naïve ones in plasma cells. Unsupervised clustering of GSVA scores for plasma cells from treated and treatment-naïve patients. Potential developmental trajectory of CD8+ T cells inferred by Monocle. 



utilities.R

Dependencies, variables and functions used in R script.
