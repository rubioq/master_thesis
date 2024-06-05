
This repository contains the code used in the work "Stress granule protein Tudor Staphylococcal Nuclease mediates heat response at the multi-omic level in *Arabidopsis thaliana*". This study represents a Master's thesis (Master's degree in Omic Data Analysis and Systems Biology).  

## Scripts

- **functions.R** - Used to define various functions helpful for the processing of phosphoproteomic, proteomic and transcriptomic data, as well as downstream analyses specific to our research field. 

- **data_preprocessing_protein.R** - Script for the preprocessing of proteomic data, with steps such as filtering and quantile normalization. It also includes Principal Component Analysis (PCA).

- **data_preprocessing_phospho.R** - Analogous to the previous script, but for phosphoproteomic data, where similar steps are carried out. In this case, we also include a correction with proteomic data step and the generation of pie charts summarizing phosphoproteomic identification results.

- **phospho_diff_analysis.R** - Differential phosphorylation analysis using *limma* and GO term Over Representation Analysis (ORA). Using a custom function, results can be exported to a plain text file.

- **iupred_analysis.ipynb** - Python script for the *in silico* prediction of intrinsic disorder using IUPred2A library (https://iupred2a.elte.hu/). All *Arabidopsis* sequences from the TAIR10 database were analyzed in our work.

- **phospho_downstream_analysis.R** - Includes additional analyses performed on differentially phosphorylated sites (DPSs) identified in phospho_diff_analysis.R. We mainly studied cellular localization of proteins, their physicochemical characteristics (in relation to Liquid-Liquid Phase Separation and Intrinsic Disorder) and carried out a clustering analysis.

- **protein_diff_analysis.R** and **transc_diff_analysis.R** - Used for the identification of differentially expressed proteins (DEPs) and genes (DEGs), respectively, as well as the corresponding ORAs.

- **multi_omic_venn_diagrams.R** - Generates Venn diagrams showing the overlap between differentially phosphorylated proteins, DEPs and DEGs in the three contrasts of interest in our work.