# Purpose
This repo has been created to share the code and details of the modes of gene expression regulation project for the iModEst paper and tool development. 

The modes of regulation project seeks to identify the relationship between different gene expression regulators and gene expression itself. This has been done by fitting elastic net models associating each of a set of regulators with the genes they are proposed to target. Details on how the regulators were chosen and an in-depth description of the methodology can be found in the [iModEst manuscript draft](https://docs.google.com/document/d/1Cr4qanLEvsm4gJR2_z3cMIneZFDzSYjcFGg4n0TVO1s/edit?usp=sharing). 

# Organization of this repo 
This repo contains the scripts used to create the figures and conduct the analyses detailed in the [iModEst manuscript draft](https://docs.google.com/document/d/1Cr4qanLEvsm4gJR2_z3cMIneZFDzSYjcFGg4n0TVO1s/edit?usp=sharing). \
   (1) **generate_model_results**: Scripts to generate the PRESS R&#x00B2; and model coefficients. These scripts reproduce the fundamental analyses of our project and can be used for reference but don't need to be run (they take a long time since they implementent a leave-one-out procedure on all of our data). \
   (2) **process_model_results**: Scripts to generate the added the PRESS R&#x00B2; for a given regulator group and gene as well as the average and weighted average coefficients for a given individual regulator and gene.\
   (3) **graphs**: Scripts to reproduce the graphs in the paper with the same color. \
   (4) **annotate_regs**: Scripts to annotate regulators. \
   (5) To potentially be added: preprocessing scripts for reference

# Data
The data for this project comes in several forms. You only need to worry about the PRESS R&#x00B2; and model coefficient data: \
   (1) **Raw data**: raw expression, methylation, CNV, and SNP data

   (2) **Preprocessed data**: expression, methylation, CNV, and SNP data which has been pre-processed for quality and corrected for confounders. Preprocessing steps detailed in the [iModEst manuscript draft](https://docs.google.com/document/d/1Cr4qanLEvsm4gJR2_z3cMIneZFDzSYjcFGg4n0TVO1s/edit?usp=sharing).

   (3) **PRESS R&#x00B2; data**: The results from our analyses represent 64 [PRESS R&#x00B2;](https://en.wikipedia.org/wiki/PRESS_statistic) statistics across over 10,000 genes in each of 21 cancers. These data are modeled unto themselves so that users can gain useful insights from broader trends in the genomic data. These data are stored in ".rds" files and can be read in to R using the readRDS function in base R. The PRESS R&#x00B2; data come in two forms: \
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(a) Raw PRESS R&#x00B2; data in which each .rds file (named CANCER_0_CANCER_0_r_square_1.rds) is organized in a matrix format with models in the rows and genes (labelled by their ensemble IDs) in the columns. The models are labelled with "glm_<model>" where <model> contains some combination of the numbers 1-6 (e.g. glm_1456). The numbers are always in sequential order and represent which regulator groups were included in that given model. The names are mapped from numbers to regulators in the following way: 1 = miRNA, 2 = TF, 3 = lncRNA, 4 = methylation, 5 = CNV, 6 = SNPs. \
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(b) Processed PRESS R&#x00B2; data in which each .rds file (named CANCER-var-added-list.rds) is organized into a list with elements named by gene symbol. Within each element of the list (i.e. for each gene), there is a dataframe with the following columns: "var.added","based.model","r2.growth","mod.r2","var.added.f","gene". Where "var.added" = the variable (1-6) that has been added to the model, "base.model" = the model that the variable was added to, "r2.growth" = the increase (therefore negative values represent a decrese) in PRESS R&#x00B2; when "var.added" was added to the "base.model", "mod.r2" = PRESS R&#x00B2; value of the model with the "var.added" *included*, "var.added.f" = same as "var.added" but the regulator type is a factor label (e.g. "var.added" = 1, "var.added.f" = "miRNA"), "gene" = ensemble ID of the gene. 
   So, an example element of this file would be a list element named "LYSMD2". One row of the dataframe corresponding to this gene would have "var.added" = 1, "base.model" = glm_23456 (therefore, you're evaluating the increase in PRESS R&#x00B2; from glm_23456 to glm_123456), "r2.growth" = -0.05539 (meaning the PRESS R&#x00B2; went down when miRNA were added to the glm_23456 model), "mod.r2" = 0.732 (meaning the PRESS R&#x00B2; for glm_123456 is 0.732), "var.added.f" = miRNA (since "var.added" = 1), and "gene" = "ENSG00000000457.9" (since this is the ensemble ID for LYSMD2). 

   (4) **Model coefficient data**: Model coefficient data comes in two forms \
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(a) Raw model coeffiecients direct from each of the models that were fit (note that the model coefficients here are coefficients computed with the full dataset, not a leave-one-out dataset). There are separate directories for each cancer which contain an .rds file for each model fit. These .rds files contain a list object with list elements named by gene. Each element of the list is a vector of the coefficients for that gene for the given model. 

   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(b) Weighted (and unweighted) average model coefficients in which the average coefficient weighted by the model's PRESS R&#x00B2; has already been computed. These data are stored in .rds files which contain a list named by gene ensemble IDs. Each list element has 2 sublist elements named "nowt" and "wt" where "nowt" is a vector of the unweighted coefficient averages and "wt" is a vector of the weighted coefficient averages. These can be found in:



# Data locations on HPF
(1) Press R&#x00B2; statistics for each cancer:
```
/hpf/largeprojects/agoldenb/mingjie/temp/predict_on_same_tissue_cohort_Sep19_cached/CANCER_0_tumor_r_square_1.rds \
```
(2) Summary files for each cancer:
```
/hpf/largeprojects/agoldenb/mingjie/temp/data_summary_Sep16_cached 
```
(3) Press R&#x00B2; statistics when models are applied to tumour adjacent tissue for each cancer \
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(a) Raw PRESS R&#x00B2; data
```
/hpf/largeprojects/agoldenb/mingjie/temp/predict_on_same_tissue_cohort_Sep19_cached/CANCER_0_CANCER_0_r_square_1.rds \
```
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(b) Processed PRESS R&#x00B2; data
```
/hpf/largeprojects/agoldenb/lauren/Modes/gene-clustering-out/CANCER-var-added-list.rds
```
(4) Model coefficients for each cancer \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(a) Raw model coeffiecients
```
/hpf/largeprojects/agoldenb/mingjie/temp/model_coeff_Sep19_cached
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(b) Weighted (and unweighted) average model coefficients
```
/hpf/largeprojects/agoldenb/lauren/Modes/model_coef_summary
```
Note CANCER follows the [TCGA abbreviations](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations) and can take on one of the following values: BRCA, BLCA, CESC, ESCA, HNSC, KIRC, KIRP, LGG, LIHC, LUAD, LUSC, PAAD, PCPG, PRAD, READ-COAD, SARC, SKCM,  STAD, TGCT, THCA, UCEC


