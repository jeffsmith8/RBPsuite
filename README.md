![alt text](https://github.com/jeffsmith8/RBPsuite/blob/master/RBPsuite%20logo.png?raw=true)
  
__RBPsuite v1.0: A Proteomics Pipeline for RBP Affinity Capture__  
_Written by Jeff Smith, The Walter and Eliza Hall Institute, 2020_  
Platform: Linux, Windows, MACOSX
  
_RBPsuite_ is a group of python tools useful for processing MS data. It is intended to provide a complete pipeline for the qualitative identification of RNA-binding proteins captured via the accompanying RNP purification method. For comparative purposes tools are also provided for a conventional affinity analysis as would be conducted for RNA-Interactome Capture (Castello/Baltz 2012) or for standard MS identification of binding partners by immunoprecipitation.  
  
A fully realised object-oriented programming (OOP) structure is not intended for this package. It is intended to be usable by researchers with beginner-intermediate python skills, and to do so with minimal notebook code to improve overall readability of an experimental analyses.  
  
A brief overview of _some_ of the features gives users the following options:  
  
__Routine Processing__  
_These tools clean up MaxQuant outputs_  
* Manipulation, renaming of samples, and rewriting of Max Quant .txt outputs  
* Re-annotation of protein and gene names with any user-selected convention utilising a stable picking algorithm  
* Limits the impact of peptide false transfer rates (due to Match Between Runs) on the identification of proteins by Intensity or iBAQ intensity  
* Removes gene duplication by disqualifying protein isoforms from generating multiple protein groups for disparate parents  
  
__Experimental QC__  
_All tools can be applied on a group-wise or individual sample level_  
* Extraction and review of sample contaminants  
* Assessment of digestion efficiencies  
* Sample clustering for routine experimental QC  
* Counting of unique genes and proteins  
* Counting of peptides  
* Analysis of normalisation impacts where MaxQuant LFQ algorithms are applied  
* Analysis of peptide Intensity, LFQ intensity or iBAQ intensity distribution  
* Review sequence coverage  
* Review replicate correlations and pearson-r statistics permuted between samples of any group  
* Review missing values and imputation depth  
  
__Affinity Analysis__  
_Tools for identifying enriched or purified proteins/genes_  
* Quantitative T-Testing, Imputation, FDR correction (BH, q-val), Fold change calc  
* Qualitative Protein/Gene Classifier for Affinity Captures by Purification  
  
__Molecular Analysis__  
_Tools for Characterising Protein Biochemistry and Exploring Reported Gene Functions_  
* Isoelectric Point Calculations  
* Real-time Quick-Go query and record extraction  
* Quickgo record manipulation and filtering tools  
  
__Visualisation Aids__  
_Prerolled functions to rapidly plot analyses_  
* Static and interactive volcano plots  
* Heatmaps  
* Grouped/individual/stacked barplots and histograms  
* Scatterplots and correlation reporting  
* Dendrogram  
* Split sample violin plots  
* Box and whisker plots  
* KDE/Density plots  
* Upsetplots  
  
__Utility Tools__  
_Tools to make life easier_  
* DataFrame manipulation for plotting, i.e. table reorganisation by metadata classes, long/wide form tabulation  
* Calculation and extraction of descriptive statistics by metadata classes and intensity types  
* Conversion of flat dataframes to multi-indexed dataframes  
* Reading in mixed files type from folders _en masse_  
* Dataframe annotation by frequency counts and other statistics   
