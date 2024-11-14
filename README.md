# **Neonatal-HSV-Proteomics**

This repository contains R scripts that conduct data analysis for OLINK proteomics data used in the article, *Proteomic profiling of neonatal herpes simplex virus infection on dried blood spots: exploring screening perspectives,* published in *Communications Medicine.* The analysis includes data quality control, principal component analysis (PCA) for outlier detection, group-based association testing, ANOVA, and permutation tests. Additionally, it generates several visualizations, including PCA plots, forest plots, heatmaps, and boxplots, used in the publication.

---

## **Usage**

### **Data Preparation**

1. **Update the file paths**: Prepare your covariate and OLINK assay data files:
   * **covariates**: Contains data such as `SampleID`, `Bday` (Birthday), `GA` (Gestational Age in weeks), `BW` (Birthweight in grams), `Sex` (M/F), `Severity` (Control, SEM, CNS, DIS), and `caco` (case/control).
   * **olink_data**: Contains OLINK assay data, including variables such as `SampleID`, `Assay`, `QC_Warning`, `NPX`, etc.

### **Requirements**

* The script automatically installs any missing packages required for the analysis.

---

## **Script Structure and Functions**

### **1. Import Packages and Data**
* Checks for required packages and loads them. Automatically installs any missing packages.

### **2. Data Merging and Initial Quality Control**
* Merges covariate and assay data.
* Filters out "CTRL" assays for further analysis.

### **3. PCA Outlier Detection**
   * Performs Principal Component Analysis (PCA) for each assay panel to detect and flag outliers based on a 3-standard-deviation rule.

### **4. REDCAP (Group) Quality Control**
   * Ensures redcap groups match in `sex`, `GA` (gestational age), and `BW` (body weight).
   * Removes any entries with mismatches in these criteria.

### **5. REDCAP (Group)-Pairs Quality Control**
   * Ensures that each assay group contains both case and control samples.
   * Excludes unpaired assays from analysis.

### **6. Covariate Associations**
   * Checks for associations between NPX values and covariates, including `SEX`, `BW`, `GA`, and `YEAR`.

### **7. ANOVA Testing and Posthoc Test**
   * Conducts ANOVA to assess the association between protein expression (NPX) and severity categories.
   * Significant results undergo a posthoc test for further insights.

### **8. Permutation Testing**
   * Computationally intensive section that runs permutation tests on ANOVA results to ensure robustness.

---

## **Plots**

### **PCA Outlier Detection**
* Visualizes PCA outliers across assay panels.

### **Forest Plot**
* Shows significant differences between control and disease groups.

### **Heatmap**
* Visualizes NPX values for proteins with posthoc significance across severity groups.

### **Boxplot**
* Compares NPX levels across assays and phenotypes for selected significant assays.

---

This README provides a clear structure for users and researchers to understand and use the repository effectively.
