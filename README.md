# Neonatal-HSV-Proteomics

The R scripts conducts data analysis for OLINK proteomics data used in article the "Proteomic profiling of neonatal herpes simplex virus infection on dried blood spots: exploring screening perspectives" published in Communications Medicine.
The analysis involving data quality control, principal component analysis (PCA) for outlier detection, group-based association testing, ANOVA, and permutation tests. Additionally, it generates several visualizations, including PCA plots, forest plots, heatmaps, and boxplots used in the publication. 


Usage:
Data Preparation
Update the file paths: covariates and olink_data.
covariates : Contains covariate data: SampleID, Bday(Birthday), GA(gestational age, weeks), BW(Birthweigth, grams), Sex(M/F), Severity(Control, SEM, CNS, DIS), caco (case/control)
olink_data: OLINK assay data, including variables like SampleID, Assay, QC_Warning, NPX etc..

Requirements:
The script installs any missing packages automatically.

Script Structure and Functions:
Import Packages and Data
Data Merging and Initial Quality Control
PCA Outlier Detection:
  - This section performs Principal Component Analysis (PCA) for each assay panel to detect and mark outliers based on the 3-standard-deviation rule.
REDCAP(group) Quality Control
  - This section checks are conducted to ensure that redcap groups match in sex, GA (gestational age), and BW (body weight), and any mismatched entries are removed.
REDCAP(group)-pairs Quality Control
  - This section ensures that each assay group has both case and control samples
Covariate Associations:
  -  This section checks for associations between NPX values and covariates, including SEX, BW, GA, and YEAR
ANOVA Testing and posthoc test:
  - ANOVA is performed to assess the association between protein expression (NPX) and severity categories. Significant results undergo a posthoc test
Permutation Testing 
  - Very computationally intensive section that runs permutation tests on the ANOVA results to ensure robustness

Plots:
PCA outlier detection
Forest
Heatmap
Boxplot
    
