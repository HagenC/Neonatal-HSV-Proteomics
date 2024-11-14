#Packages: ----
packages <- c("data.table", 
              "tidyr",
              "dplyr",
              "forcats",
              "ggplot2", 
              "tibble", 
              "OlinkAnalyze",
              "broom.mixed",
              "lmerTest",
              "doParallel",
              "pheatmap")

# check if installed and load 
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}


#Data paths:
covariates <- "<filepath>"
olink_data <- "<filepath>"
#Importing data ----
COVAR <- fread(covariates, na.strings = "") %>%   filter(!if_all(everything(), is.na))
#Variables: #"bxpID","Bday","SampleID","redcap","caco","GA" ,"BW","Sex" ,"Severity"  
OLINK  <- fread(olink_dat, na.strings = "")
#Variables "SampleID", "UniProt","Assay","MissingFreq","Panel","Panel_Lot_Nr","PlateID","QC_Warning","LOD","NPX","Normalization","Assay_Warning"

#Merging OLINK data and covariate data and removing CTRL assays
OLINK_COVAR <- left_join(COVAR,OLINK, by = c("bxpID" = "SampleID")) %>% filter(Assay != "CTRL")

#Extracting PASS (passed) QC and Assay-warnings:
OLINK_COVAR_PASS <- OLINK_COVAR  %>% filter(Assay_Warning == "PASS" & QC_Warning == "PASS") 

#PCA outlier detection: ----
#NOTE: apply serveral iterations for more homogenious population.
#Removing incomplete Assays:
incomplete.assays <- OLINK_COVAR %>%  group_by(Panel, Assay, QC_Warning) %>% count(Assay_Warning) %>% filter(QC_Warning != "PASS" | Assay_Warning != "PASS") %>% unite(SELECTION, c(Assay, Panel), sep = "/")
complete.assay <- OLINK_COVAR %>%  unite(SELECTION, c(Assay, Panel), sep = "/") %>% filter(!SELECTION %in% incomplete.assays$SELECTION) %>% count(SELECTION)
#Listing panels:
panel.list <- complete.assay  %>% separate(SELECTION, c("Assay", "Panel"), sep = "/") %>% count(Panel) %>% pull(Panel)

#Plotting and marking PCA outliers:
outliers <- list()
C = 1
for(i in panel.list){
  PCA <- OLINK_COVAR %>%  filter(Panel == i) %>% unite(SELECTION, c(Assay, Panel), sep = "/") %>%  filter(SELECTION %in% complete.assay$SELECTION) %>% 
    separate(SELECTION, c("Assay", "Panel"), sep = "/") %>%  select(SampleID, Assay, NPX) %>% 
    spread(Assay, NPX) %>% column_to_rownames(., "SampleID") %>% prcomp(., scale = TRUE) 
  PCs <- data.frame(-1*PCA$x)
  SD3.PC1 <- PCs %>% summarise(PC1.m = mean(PC1), PC1.sd = sd(PC1)) %>% mutate(SD3 = PC1.m + 3*PC1.sd)
  SD3.PC2 <- PCs %>% summarise(PC2.m = mean(PC2), PC2.sd = sd(PC2)) %>% mutate(SD3 = PC2.m + 3*PC2.sd)
  explained <- data.frame(EXP = PCA$sdev^2 / sum(PCA$sdev^2)) %>% mutate(EXP = format(EXP, scientific = F)) %>% mutate(EXP = round(as.numeric(EXP)*100, digits = 1))
  PCs.color <- PCs %>% mutate(outlier = ifelse(PC1 < -SD3.PC1$SD3 | PC1 > SD3.PC1$SD3 |PC2 < -SD3.PC2$SD3 | PC2 > SD3.PC2$SD3, "outlier", "" ))
  outliers[[C]] <- PCs.color %>% filter(outlier == "outlier")  %>% mutate(SampleID = row.names(.), Panel = i) %>% select(SampleID, Panel) 
  print(i)
  C <- C+1
  fig.out <- PCs.color %>% ggplot(., aes(PC1, PC2, color = outlier)) + geom_point() + geom_vline(xintercept = SD3.PC1$SD3, linetype = "dotted", color = "red") + geom_vline(xintercept = -SD3.PC1$SD3, linetype = "dotted", color = "red" )  +
    geom_hline(yintercept = SD3.PC2$SD3, linetype = "dotted", color = "red") + geom_hline(yintercept = -SD3.PC2$SD3, linetype = "dotted", color = "red") +
    scale_color_manual(values=c("#99d594", "#d53e4f")) + theme(legend.position = "none") + ggtitle(i) + xlab(paste0("PC1 ", explained$EXP[1], "%")) +
    ylab(paste0("PC2 ", explained$EXP[2], "%"))
  print(fig.out)
}


outlier_samples <- do.call("rbind", outliers) %>% left_join(., COVAR, by = c("SampleID" = "SampleID")) 
outlier_samples_removelist <- outlier_samples %>% select(SampleID, Panel) %>% mutate(Outlier = "Outlier") 

#Removing outliers:
OLINK_COVAR_PCA <- OLINK_COVAR_PASS %>% left_join(., outlier_samples_removelist, by = c("SampleID", "Panel")) %>% filter(is.na(Outlier))


#REDCAP-sex: Checking for grouping (redcap) sex match ----
redcap_sex_check <- COVAR %>% group_by(redcap) %>% count(Sex) %>% spread(Sex, n)  %>% filter(!is.na(M) & !is.na(F)) %>% left_join(., COVAR, by = "redcap") 
  #Removing mismatches:
 redcap_sex_check
 OLINK_COVAR_PCA_SEX <- OLINK_COVAR_PCA  %>% filter(!bxpID %in% c("BX0345_051"))

#REDCAP-GA: hecking for grouping (redcap) GA mismatch setting limits to 1 (week) ----
redcap_GA_check <- COVAR %>% group_by(redcap) %>% select(redcap, GA) %>% group_by(redcap) %>% summarise(mean.GA = mean(GA), max.GA = max(GA), min.GA = min(GA), SD.GA = sd(GA)) %>% filter(max.GA - min.GA > 1) 
redcap_GA_check
#REDCAP-BW : Checking for grouping (redcap) BW where standard deviation of body weight (BW) is is unusually high compared to other groups (> Q3 + 1.5 IQR).----
redcap_BW_check <- COVAR %>% filter(!is.na(BW)) %>% select(redcap, BW) %>% group_by(redcap) %>% 
  summarise( SD.BW = sd(BW)) %>% ungroup() %>%  
  mutate(Q3 = quantile(SD.BW, 0.75), IQR = IQR(SD.BW))  %>% mutate(H.limit = Q3 + 1.5*IQR) %>% 
  filter(H.limit < SD.BW) %>% left_join(., COVAR, by = "redcap") 
  #Removing outliers:
  redcap_BW_check
  OLINK_COVAR_PCA_SEX_BW <- OLINK_COVAR_PCA_SEX  %>% filter(!bxpID %in% c("BX0345_004"))

#REDCAP-pairs: Removing unpaired redcap Assays i.e Assayes where either case or controls is missing. ---- 
caco_complete_list <- OLINK_COVAR_PCA_SEX_BW %>% select(Assay, Panel, redcap, caco) %>%  group_by(Assay, Panel, redcap) %>% count(caco) %>% spread(caco, n) %>% 
    filter(!is.na(Case) & !is.na(Control)) %>% unite(Assay_Panel_redcap, c("Assay", "Panel", "redcap"), sep = "-")
OLINK_COVAR_PCA_SEX_BW_CACO <- OLINK_COVAR_PCA_SEX_BW %>%  unite(Assay_Panel_redcap, c("Assay", "Panel", "redcap"), sep = "-", remove = FALSE) %>% 
  filter(Assay_Panel_redcap %in% caco_complete_list$Assay_Panel_redcap) 

#Renaming for ease:
MODEL_data <- OLINK_COVAR_PCA_SEX_BW_CACO %>% unite(SELECTION, c(Assay, Panel), sep = "-", remove = FALSE) 

#cleaning
       rm(list = setdiff(ls(), "MODEL_data"))
       gc()
#Covar-assocation ----
#Checking for NPX association with covariates: SEX, BW, GA and YEAR
#Year with linear model:
lm_YEAR <-MODEL_data %>% mutate(Year = as.integer(substr(Bday, 7,10))) %>%  select(SELECTION, Year, NPX) %>%
  group_by(SELECTION ) %>%  dplyr::group_modify(~generics::tidy(stats::lm(Year ~NPX, data = .)))
lm_YEAR.term <- lm_YEAR %>% filter(term == "NPX")  %>% arrange(p.value) %>% 
  mutate(FDR.p.value = p.adjust(p.value, "fdr", nrow(.)))
nrow(lm_YEAR.term %>% filter(FDR.p.value < 0.05))/nrow(lm_YEAR.term)*100
#

# GA mixed model with redcap as a random effect using lmerTest
lm_GA_mixed <- MODEL_data %>%
  filter(caco == "Control") %>%
  select(SELECTION, GA, NPX, redcap) %>%
  group_by(SELECTION) %>%
  dplyr::group_modify(~ broom.mixed::tidy(lmerTest::lmer(GA ~ NPX + (1 | redcap), data = .), effects = "fixed"))
lm_GA_mixed.term <- lm_GA_mixed %>%
  filter(term == "NPX") %>% ungroup() %>% 
  mutate(FDR.p.value = p.adjust(p.value, method = "fdr")) 
nrow(lm_GA_mixed.term %>% filter(FDR.p.value < 0.05))/nrow(lm_GA_mixed.term)*100

# BW mixed model with redcap as a random effect using lmerTest
lm_BW_mixed <- MODEL_data %>%
  filter(caco == "Control") %>%
  select(SELECTION, BW, NPX, redcap) %>%
  group_by(SELECTION) %>%
  dplyr::group_modify(~ broom.mixed::tidy(lmerTest::lmer(BW ~ NPX + (1 | redcap), data = .), effects = "fixed"))
lm_BW_mixed.term <- lm_BW_mixed %>%
  filter(term == "NPX") %>% ungroup() %>% 
  mutate(FDR.p.value = p.adjust(p.value, method = "fdr")) %>% filter(FDR.p.value < 0.05) 

# Sex logistic regression 
lm_sex <- MODEL_data %>% mutate(sex.binary = ifelse(Sex == "M",0, 1)) %>% 
  select(SELECTION, sex.binary, NPX) %>% 
  group_by(SELECTION) %>%
  dplyr::group_modify(~ broom::tidy(glm(sex.binary ~ NPX, data = ., family = binomial)))
lm_sex.term <- lm_sex %>%
  filter(term == "NPX") %>% ungroup() %>% 
  mutate(FDR.p.value = p.adjust(p.value, method = "fdr")) %>% filter(FDR.p.value < 0.05) 

#ANOVA ----
MODEL_data.ANOVA <- MODEL_data %>%  mutate(Severity=fct_relevel(Severity,c("Control","SEM","CNS", "DIS")))
anova_results <- olink_anova(df = MODEL_data.ANOVA, variable = "Severity") %>% unite(SELECTION, c(Assay, Panel), sep = "-", remove = FALSE)
anova_results_sub <- anova_results %>% select(SELECTION,  p.value, Adjusted_pval, Threshold) %>% rename(anova.p.value = p.value, anova.Adjusted_pval = Adjusted_pval, anova.Threshold = Threshold)# %>% filter(Threshold == "Significant")

sig_anova_results <- anova_results %>% distinct() %>% pull(OlinkID) 
anova_posthoc_results <- olink_anova_posthoc(MODEL_data.ANOVA,
                                                 variable=c("Severity"),
                                                 covariates=NULL,
                                                 olinkid_list = sig_anova_results,
                                                 effect = "Severity")
results_DB <- anova_posthoc_results %>%   unite(SELECTION, c(Assay, Panel), sep = "-", remove = FALSE) %>% left_join(., anova_results_sub, by = "SELECTION")

#Plots
#Severity combinations:
contrast_combinations  <- results_DB %>% count(contrast) %>% separate(contrast, c("A", "B"), sep = " - ")
#mean difference in NPX grouped by Severity
mean.diff <- MODEL_data.ANOVA %>% select(Severity, SELECTION, NPX) %>% group_by(SELECTION,Severity ) %>% summarise(Mean = mean(NPX)) %>% spread(Severity, Mean) 

collect_diff <- list()
for(i in 1:nrow(contrast_combinations)){
  collect_diff[[i]] <- mean.diff %>% select(SELECTION, contrast_combinations$A[i], contrast_combinations$B[i]) %>% 
    mutate(contrast = paste0(colnames(.)[2], " - ", colnames(.)[3])) %>% 
    rename(A = 2, B = 3) %>% 
    mutate(NPX.diff = ifelse(  A < 0 & B < 0 & A < B, abs(A - B), NA)) %>% 
    mutate(NPX.diff = ifelse( is.na(NPX.diff) & A < 0 & B < 0 & A > B,  A - B, NPX.diff)) %>% 
    mutate(NPX.diff = ifelse( is.na(NPX.diff) & A < 0 & B > 0 ,  B - A, NPX.diff)) %>% 
    mutate(NPX.diff = ifelse( is.na(NPX.diff) & A > 0 & B < 0 ,  -(A - B), NPX.diff)) %>% 
    mutate(NPX.diff = ifelse( is.na(NPX.diff) & A > 0 & B > 0 & A > B , B - A , NPX.diff)) %>% 
    mutate(NPX.diff = ifelse( is.na(NPX.diff) & A > 0 & B > 0 & A < B , B - A , NPX.diff)) 
}

collect_diff_df <- do.call("rbind",collect_diff ) 

results_DBs_NPXdiff <- results_DB %>% left_join(., collect_diff_df, by = c("SELECTION", "contrast")) %>% rename(NPX.mean.difference = NPX.diff) %>% 
  mutate(Significans = ifelse(Threshold == "Significant" & anova.Threshold == "Significant" , "Post hoc significant", "Non-significant"))

#Plot:
#Lines:
contrast_list <- results_DBs_NPXdiff %>% filter(!duplicated(contrast)) %>% pull(contrast)
for(i in contrast_list){
max_postsignificant  <- results_DBs_NPXdiff %>% filter(contrast == i) %>% filter(Adjusted_pval < 0.05 & Threshold == "Significant") %>% filter(Adjusted_pval == max(Adjusted_pval)) %>% pull(Adjusted_pval)
sig_control  <- results_DBs_NPXdiff %>% filter(contrast == i) %>% filter(Adjusted_pval < 0.05 & NPX.mean.difference < 0 & Threshold == "Significant") %>% filter(NPX.mean.difference == max(NPX.mean.difference)) %>% pull(NPX.mean.difference)
sig_cases  <-  results_DBs_NPXdiff %>% filter(contrast == i) %>% filter(Adjusted_pval < 0.05 & NPX.mean.difference > 0  & Threshold == "Significant" ) %>% filter(NPX.mean.difference == min(NPX.mean.difference)) %>% pull(NPX.mean.difference)
P <- results_DBs_NPXdiff %>% filter(contrast == i) %>% ggplot(., aes(NPX.mean.difference , -log10(Adjusted_pval), color = Threshold)) + geom_point() +
  geom_text(data = results_DBs_NPXdiff %>%  
              mutate(color.group = ifelse(Significans == "Post hoc significant", "POST.SIG", "GRAY" )) %>% filter(contrast == i & Significans == "Post hoc significant"),
            aes(NPX.mean.difference, -log10(Adjusted_pval), label = SELECTION), size = 3, vjust = 1, position = position_jitter(height =  0.2)) + 
            geom_hline(yintercept = -log10(max_postsignificant ), lty = 2) + 
            geom_vline(xintercept = sig_control, lty = 2) + 
            geom_vline(xintercept = sig_cases, lty = 2) + 
            scale_color_manual(values = c("#bdbdbd","#f03b20"))  + ggtitle(paste0("ANOVA (no-selection): ",i), subtitle = "Post-hoc significant assays") + xlab("NPX mean difference: ") + theme(legend.position = "none")
 print(P)
}

#Permutations: (Warning - may take days depending on aviable cores!) ----
registerDoParallel(12) #adjust for speed

SELECTION.list <- MODEL_data.ANOVA %>% filter(!duplicated(SELECTION)) %>% pull(SELECTION)

#Anova function
anova_F_function <- function(df) {
  result.ANOVA <- summary(aov(NPX ~ Severity, data = df))
  f.stat <- result.ANOVA[[1]]$`F value`[1]
  return(f.stat)
}


k <- 1
#SELECTION wise permuations:
foreach(k = 1:length(SELECTION.list)) %dopar% {
  data.SELECTION <- MODEL_data.ANOVA %>% filter(SELECTION == SELECTION.list[k])
  obs.f.stat <- anova_F_function(data.SELECTION)
  obs.p.value <- summary(aov(NPX ~ Severity, data = data.SELECTION))[[1]]$`Pr(>F)`[1]
  perm.F.stats <- list()
  N <- 0
  above.count <- 0
  while(above.count < 1000 & N <= 20000000){
    N <- N + 1
    mix.group <- sample(data.SELECTION$Severity)
    perm.df <- data.frame(NPX = data.SELECTION$NPX, Severity = mix.group)
    perm.F.stats[N] <- anova_F_function(perm.df)
    if((anova_F_function(perm.df) >= obs.f.stat) == TRUE){
      above.count <- above.count + 1
    }
  }
  RETURN <- data.frame(perm.P.value = sum(perm.F.stats >= obs.f.stat) / N, org.P.value = obs.p.value, N = N)
  fwrite(RETURN, paste0(SELECTION.list[k],".stats"), sep = "\t")
}

#Import (When finished )
perm.list <- list.files(".", "*.stats$")
collect.perm <- list()
for(i in 1:length(perm.list)){
  collect.perm[[i]] <- fread(paste0(perm.list[i])) %>% mutate(ID = perm.list[i] )
}
permutation_df <- do.call("rbind", collect.perm)
permutation_df_annotated <- permutation_df %>% mutate(ID = gsub(".stats", "", ID)) %>% mutate(P.org.adjusted = p.adjust(org.P.value, "fdr", n())) %>% mutate(P.perm.adjusted = p.adjust(perm.P.value, "fdr", n())) %>% 
  mutate(LABEL.org = ifelse(P.org.adjusted < 0.05, ID, NA)) %>% separate(LABEL.org, c("Assay.org", "Panel"), sep = "-") %>% mutate(LABEL.perm = ifelse(P.perm.adjusted < 0.05, ID, NA)) %>% 
  separate(LABEL.perm, c("Assay.perm", "Panel"), sep = "-") %>% arrange(perm.P.value, org.P.value)
#Correlation:
summary(lm(perm.P.value ~ org.P.value, data = permutation_df_annotated)) 



#Forest plot: Controls -DIS ----
ControlDis <- results_DBs_NPXdiff %>% filter(contrast == "Control - DIS" & Significans == "Post hoc significant") %>% select(Assay, estimate, conf.low, conf.high) %>% 
  mutate(Assay = factor(Assay, levels = unique(.$Assay))) %>%  mutate(across(c(estimate, conf.low, conf.high), ~ . * -1))

# Forest plot
ggplot(ControlDis, aes(x = Assay, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange(color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "orange") +
  coord_flip() +
  labs(title = "DIS vs controls - posthoc proteins",
       x = "",
       y = "npx + 95% CI") +
  theme_minimal()


#Heatmap: ----
sig_proteins <- results_DBs_NPXdiff %>% filter(Significans == "Post hoc significant") %>% filter(!duplicated(SELECTION)) %>% filter(!SELECTION %in% c("IL6-Cardiometabolic", "IL6-Oncology", "IL6-Neurology")) %>% pull(SELECTION)
prep_heat <- MODEL_data.ANOVA %>% filter(SELECTION %in% sig_proteins) %>% select(Severity, SELECTION, NPX) %>% group_by(SELECTION,Severity ) %>% 
  summarise(Mean = mean(NPX)) %>% mutate(Mean = as.numeric(Mean)) %>% spread(SELECTION, Mean) %>%   mutate(across(-Severity, ~ (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))))
prep_heat.matrix <- as.matrix(prep_heat[,2:length(prep_heat)]) 
rownames(prep_heat.matrix) <- prep_heat$Severity
pheatmap(prep_heat.matrix)

#Boxplot: ----
sig_Assay <- results_DBs_NPXdiff %>% filter(Significans == "Post hoc significant") %>% filter(!duplicated(Assay))  %>% pull(Assay)
boxplot_selection <- MODEL_data.ANOVA %>% filter(Assay %in% sig_Assay)

boxplot_selection %>% rename(Phenotype = Severity) %>%
  ggplot(aes(Assay, NPX, fill = Phenotype, color = Phenotype)) + 
  geom_boxplot(outlier.size = 2, alpha = 0.6, width = 1, position = position_dodge(width = 0.8)) + 
  geom_jitter(position = position_dodge(width = 0.8), alpha = 0.6, size = 1.5) +  
  scale_fill_manual(values = c("#FF7F50", "#4682B4", "#6A5ACD", "#B8860B")) + 
  scale_color_manual(values = c("#8B3E2F", "#1C3D6E", "#4B0082", "#B8860B")) + 
  labs(x = "Assay", y = "NPX") +
  theme_minimal(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()) 
