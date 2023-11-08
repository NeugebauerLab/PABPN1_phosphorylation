#load bakR and tidyverse
library(bakR)
library(tidyverse)
library(pheatmap)
library(dplyr) 
library(magrittr) 
library(ggplot2) 
library(stats)

#create bakR object from the cB file
file <- "/home/jmg292/project/Exp31/bakR/cB-230512.csv"
bakr_file <- tibble(read_csv(file))
head(bakr_file)

#load the metadf file that contains sample information
metadf <- readRDS("/home/jmg292/project/Exp31/bakR/IWV_metadf.rds")

#4SA vs 4SD
#Filter out bakr_file and meta_df to compares 4SD vs 4SA
subset_list = c("Sample_CB230401_10", "Sample_CB230401_11", "Sample_CB230401_12", "Sample_CB230401_13", "Sample_CB230401_14", "Sample_CB230401_15", "Sample_CB230401_16", "Sample_CB230401_9")
tl <- c(2,2,0,2,2,2,0,2)
Exp_ID <- c(2,2,2,1,1,1,1,2)
metadf_2 <- data.frame(tl, Exp_ID)
row.names(metadf_2) <- subset_list

bakr_sub <- bakr_file %>% filter(sample %in% subset_list)

print(rownames(metadf_2))
print(unique(bakr_sub$sample))


#create bakR object
bakRData <- bakRData(bakr_sub, metadf_2)

#Fit the bakR data using the fast model
Fit <- bakRFit(bakRData)

#PCA for bakR data
FnPCA2(Fit)

#Volcano plot
plotVolcano(Fit$Fast_Fit)

#Cumulative distrubution plot
ggplot(Fit$Fast_Fit$Fn_Estimates, aes(kdeg, colour = factor(Exp_ID))) +
  stat_ecdf()

write_delim(Fit$Fast_Fit$Effects_df, "//home/jmg292/palmer_scratch/Exp31/results/bakR/Fit_effects_DF_4SA_4SD.csv")
write_delim(Fit$Fast_Fit$Fn_Estimates, "//home/jmg292/palmer_scratch/Exp31/results/bakR/Fit_Fn_Estimates_4SA_4SD.csv")
write_delim(Fit$Fast_Fit$Regularized_ests, "//home/jmg292/palmer_scratch/Exp31/results/bakR/Fit_Regularized_ests_4SA_4SD.csv")

#Run DESeq2 for 4SA vs 4SD.
library(DESeq2)

Counts <- Fit$Data_lists$Count_Matrix
conditions <- as.factor(c("4SA", "4SA", "4SA", "4SD", "4SD", "4SD", "4SD", "4SA"))

# Make the colData input for DESeq2
colData <- data.frame(conditions = conditions)
rownames(colData) <- colnames(Counts)

dds <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = colData,
                              design = ~conditions)


ddso <- DESeq(dds)

reso <- results(ddso, contrast = c("conditions", "4SD", "4SA"))

DE_df <- data.frame(XF = row.names(reso),
                    L2FC_RNA = reso$log2FoldChange,
                    DE_score = reso$stat,
                    DE_se = reso$lfcSE,
                    DE_pval = reso$pvalue,
                    DE_padj = reso$padj)

ksyn_df <- data.frame(L2FC = reso$log2FoldChange + Fit$Fast_Fit$Effects_df$L2FC_kdeg,
                      Gene = Fit$Fast_Fit$Effects_df$XF)

write_delim(ksyn_df, "//home/jmg292/palmer_scratch/Exp31/results/bakR/ksyn_df_4SA_4SD.csv")


ksyn_df$se <- sqrt(reso$lfcSE^2 + (Fit$Fast_Fit$Effects_df$se*log2(exp(1)))^2 ) 

# Calculate p-value using asymptotic Wald test
ksyn_df <- ksyn_df %>%
  mutate(pval = 2*pnorm(-abs(L2FC/se)),
         padj = p.adjust(pval, method = "BH"))

# Add conclusion at 0.01 FDR control
ksyn_df <- ksyn_df %>%
  mutate(conclusion = as.factor(ifelse(padj < 0.01, 
                                       ifelse(L2FC < 0, "Decreased txn", "Increased txn"),
                                       "Not sig.")))

# Make volcano plot
ksyn_volc <- ggplot(ksyn_df, aes(x = L2FC, y = -log10(padj), color = conclusion)) + 
  theme_classic() + 
  geom_point(size = 1) + 
  xlab("L2FC(ksyn)") + 
  ylab("-log10(padj)") + 
  scale_color_manual(values = c("blue", "orange", "gray"))

# Observe volcano plot
ksyn_volc










#write_delim(DE_df, "//home/jmg292/palmer_scratch/Exp31/results/bakR/DE_results_4SA_4SD.csv")

###Run bakR for the rest of the samples for comparison. Save dfs accordingly.
#create bakR object from the cB file
file <- "/home/jmg292/palmer_scratch/Exp31/results/bakR/cB-230512.csv"
bakr_file <- tibble(read_csv(file))
head(bakr_file)

#load the metadf file that contains sample information
metadf <- readRDS("/home/jmg292/palmer_scratch/Exp31/results/bakR/IWV_metadf.rds")

#Ensure that the columns in metadf match the samples from bakr file
print(rownames(metadf))
print(unique(bakr_file$sample))

#create bakR object
bakRData <- bakRData(bakr_file, metadf)

#Fit the bakR data using the fast model
Fit <- bakRFit(bakRData)
ggplot(Fit$Fast_Fit$Fn_Estimates, aes(kdeg, colour = factor(Exp_ID))) +
  stat_ecdf()

Fit_c <- CorrectDropout(Fit)

ggplot(Fit_c$Fast_Fit$Fn_Estimates, aes(kdeg, colour = factor(Exp_ID))) +
  stat_ecdf()

Counts <- Fit$Data_lists$Count_Matrix

conditions <- as.factor(c("4SA", "4SA", "4SA", "4SD", "4SD", "4SD", "4SD", "EV" ,"EV", "EV", "EV", "WT", "WT", "WT", "WT", "4SA"))

# Make the colData input for DESeq2
colData <- data.frame(conditions = conditions)
rownames(colData) <- colnames(Counts)

dds <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = colData,
                              design = ~conditions)


ddso <- DESeq(dds)

reso <- results(ddso, contrast = c("conditions", "4SD", "4SA"))

DE_df <- data.frame(XF = row.names(reso),
                    L2FC_RNA = reso$log2FoldChange,
                    DE_score = reso$stat,
                    DE_se = reso$lfcSE,
                    DE_pval = reso$pvalue,
                    DE_padj = reso$padj)



ksyn_df <- data.frame(L2FC = reso$log2FoldChange + Fit$Fast_Fit$Effects_df$L2FC_kdeg,
                      Gene = Fit$Fast_Fit$Effects_df$XF)

ksyn_df$se <- sqrt(reso$lfcSE^2 + (Fit$Fast_Fit$Effects_df$se*log2(exp(1)))^2 ) 

# Calculate p-value using asymptotic Wald test
ksyn_df <- ksyn_df %>%
  mutate(pval = 2*pnorm(-abs(L2FC/se)),
         padj = p.adjust(pval, method = "BH"))

# Add conclusion at 0.01 FDR control
ksyn_df <- ksyn_df %>%
  mutate(conclusion = as.factor(ifelse(padj < 0.01, 
                                       ifelse(L2FC < 0, "Decreased txn", "Increased txn"),
                                       "Not sig.")))

# Make volcano plot
ksyn_volc <- ggplot(ksyn_df, aes(x = L2FC, y = -log10(padj), color = conclusion)) + 
  theme_classic() + 
  geom_point(size = 1) + 
  xlab("L2FC(ksyn)") + 
  ylab("-log10(padj)") + 
  scale_color_manual(values = c("blue", "orange", "gray"))

# Observe volcano plot
ksyn_volc

#write_delim(Fit$Fast_Fit$Effects_df, "//home/jmg292/palmer_scratch/Exp31/results/bakR/Fit_effects_DF_all.csv")
#write_delim(Fit$Fast_Fit$Fn_Estimates, "//home/jmg292/palmer_scratch/Exp31/results/bakR/Fit_Fn_Estimates_all.csv")
#write_delim(Fit$Fast_Fit$Regularized_ests, "//home/jmg292/palmer_scratch/Exp31/results/bakR/Fit_Regularized_ests_all.csv")



################################################################
#WT vs 4SA and 4SD
#Filter out bakr_file and meta_df, set WT as control condition
subset_list = c("Sample_CB230401_10", "Sample_CB230401_11", "Sample_CB230401_12", "Sample_CB230401_13", "Sample_CB230401_14", "Sample_CB230401_15", "Sample_CB230401_16", "Sample_CB230401_1", "Sample_CB230401_2", "Sample_CB230401_3", "Sample_CB230401_4", "Sample_CB230401_5", "Sample_CB230401_6", "Sample_CB230401_7", "Sample_CB230401_8", "Sample_CB230401_9")
tl <- c(2,2,0,2,2,2,0,2,2,2,0,2,2,2,0,2)
Exp_ID <- c(2,2,2,3,3,3,3,4,4,4,4,1,1,1,1,2)
metadf_2 <- data.frame(tl, Exp_ID)
row.names(metadf_2) <- subset_list

bakr_sub <- bakr_file %>% filter(sample %in% subset_list)

print(rownames(metadf_2))
print(unique(bakr_file$sample))


#create bakR object
bakRData <- bakRData(bakr_sub, metadf_2)

#Fit the bakR data using the fast model
Fit <- bakRFit(bakRData)

#PCA for bakR data
FnPCA2(Fit)

#Volcano plot
plotVolcano(Fit$Fast_Fit, Exps = 4)

Fit_c <- CorrectDropout(Fit)

ggplot(Fit$Fast_Fit$Fn_Estimates, aes(kdeg, colour = factor(Exp_ID))) +
  stat_ecdf()

write_delim(Fit$Fast_Fit$Effects_df, "//home/jmg292/palmer_scratch/Exp31/results/bakR/Fit_effects_DF_WTref.csv")
write_delim(Fit$Fast_Fit$Fn_Estimates, "//home/jmg292/palmer_scratch/Exp31/results/bakR/Fit_Fn_Estimates_WTref.csv")
write_delim(Fit$Fast_Fit$Regularized_ests, "//home/jmg292/palmer_scratch/Exp31/results/bakR/Fit_Regularized_ests_WTref.csv")

