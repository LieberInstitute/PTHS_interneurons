# PTHS_interneurons
Genomics analyses of PTHS mouse models investigating interneuron populations.
by BaDoi Phan (@ badoi dot phan at pitt dot edu)

## Methods
The [**Tasic_2018_mouse_single-cell**](Tasic_2018_mouse_single-cell) folder contains scripts to perform 
- Cell type-specific expression analyses
- Cell type deconvolution analyses

### scripts
0a. [Download script to grab Tasic _et al_ 2018, Nature](Tasic_2018_mouse_single-cell/download.sh)

0b. [Script to tidy Tasic _et al_ data into R and Bioconductor](Tasic_2018_mouse_single-cell/step0_make_tasic_rse.R)

1. [Script to find marker genes for Tasic _et al_ cell types](Tasic_2018_mouse_single-cell/step1_findMarkerGenes.R)

3. [Script for Cell type-specific Expression Analyses](Tasic_2018_mouse_single-cell/step2_enrichMarkers_PTHS_mouse.R)

3a. [Script to prepare bulk and scRNA data for deconvolution analyses](Tasic_2018_mouse_single-cell/step3a_prepare_deconvolution.R)

3a. [Script to perform deconvolution of bulk RNA-seq samples with scRNA-seq markers](Tasic_2018_mouse_single-cell/step3b_deconvolve_PTHS_mouse_RNA-seq.py)

4. [Script to visualize deconvolution cell types from PTHS samples](Tasic_2018_mouse_single-cell/step4_analyze_allenCellular_deconvolution.R)

### outputs
The [**Tasic_2018_mouse_single-cell/plots**](Tasic_2018_mouse_single-cell/plots) folder contains source plots for CSEA and deconvolution figures.
The [**Tasic_2018_mouse_single-cell/tables**](Tasic_2018_mouse_single-cell/tables) folder contains numerical and tables underlying the figures.
