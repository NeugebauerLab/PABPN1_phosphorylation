#load bakR and tidyverse
library(bakR)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(dplyr) 
library(magrittr) 
library(ggplot2) 
library(stats)

#create bakR object from the cB file
file <- "../results/cB/cB.csv.gz"
#file <- "/home/jmg292/project/Exp31/bakR/cB-230512.csv"
bakr_file <- tibble(read_csv(file))
head(bakr_file)

sample_names <- c("4SA_1", "4SA_2", "4SA_3", "4SA_ctl", "4SD_1", "4SD_2", "4SD_3", "4SD_ctl", "EV_1", "EV_2", "EV_3", "EV_ctl", "WT_1", "WT_2", "WT_3", "WT_ctl")
tl <- c(120,120,120,0,120,120,120,0,120,120,120,0,120,120,120,0)
Exp_ID <- c(3,3,3,3,4,4,4,4,1,1,1,1,2,2,2,2)
#Exp_ID <- c(3,3,3,3,4,4,4,4,2,2,2,2,1,1,1,1)
metadf <- data.frame(tl, Exp_ID)
#metadf <- readRDS("/home/jmg292/project/Exp31/bakR/IWV_metadf.rds")
row.names(metadf) <- sample_names

print(rownames(metadf))
print(unique(bakr_file$sample))

#create bakR object
bakRData <- bakRData(bakr_file, metadf)

#Fit the bakR data using the fast model
Fit_s <- bakRFit(bakRData)

plotVolcano(Fit_s$Fast_Fit, Exps = 3)

#Visualize dropout of 4sU from samples
Vis_DO <- VisualizeDropout(Fit_s)
Vis_DO$ExpID_1_Rep_3

#Correct for dropout
Fit_c <- CorrectDropout(Fit_s)

#Volcano plot to visualize results
plotVolcano(Fit_c$Fast_Fit, Exps = 3)
FnPCA2(Fit_s)


#CDF plot
ggplot(Fit_s$Fast_Fit$Fn_Estimates, aes(kdeg, colour = factor(Exp_ID))) +
  stat_ecdf()

#Save dataframes for making figures in python
write_delim(Fit_s$Fast_Fit$Effects_df, "12142023_Fit_effects_DF_all.csv")
write_delim(Fit_s$Fast_Fit$Fn_Estimates, "12142023_Fit_Fn_Estimates_all.csv")
write_delim(Fit_s$Fast_Fit$Regularized_ests, "12142023_Fit_Regularized_ests_all.csv")


#########################################################
#Rerun the bakR analysis, subset 4SA and 4SD to directly compare and to run DESeq2
subset_list = c("4SA_1", "4SA_2", "4SA_3", "4SA_ctl", "4SD_1", "4SD_2", "4SD_3", "4SD_ctl")
tl <- c(2,2,2,0,2,2,2,0)
Exp_ID <- c(2,2,2,2,1,1,1,1)
metadf_2 <- data.frame(tl, Exp_ID)
row.names(metadf_2) <- subset_list

bakr_sub <- bakr_file %>% filter(sample %in% subset_list)

print(rownames(metadf_2))
print(unique(bakr_sub$sample))

#create bakR object
bakRData <- bakRData(bakr_sub, metadf_2)

#Fit the bakR data using the fast model
Fit_s <- bakRFit(bakRData)

#Correct for dropout
Fit_c <- CorrectDropout(Fit_s)

#Volcano plot to visualize results
plotVolcano(Fit_s$Fast_Fit, Exps = 2)

#Save dataframes for making figures in python
write_delim(Fit_s$Fast_Fit$Effects_df, "12142023_Fit_effects_DF_4SA_4SD.csv")
write_delim(Fit_s$Fast_Fit$Fn_Estimates, "12142023_Fit_Fn_Estimates_4SA_4SD.csv")
write_delim(Fit_s$Fast_Fit$Regularized_ests, "12142023_Fit_Regularized_ests_4SA_4SD.csv")



##########################################################################################
#Run DESeq2 for 4SA vs 4SD.
library(DESeq2)

Counts <- Fit_s$Data_lists$Count_Matrix
conditions <- as.factor(c("4SA", "4SA", "4SA", "4SA", "4SD", "4SD", "4SD", "4SD"))

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


# Add conclusion at 0.05 FDR control
DE_df <- DE_df %>%
  mutate(conclusion = as.factor(ifelse(DE_padj < 0.01, 
                                       ifelse(L2FC_RNA < 0, "Decreased txn", "Increased txn"),
                                       "Not sig.")))
DE_df
write_delim(DE_df, "12042023_DE_df_4SA_4SD.csv")

# Make volcano plot
ggplot(DE_df, aes(x = L2FC_RNA, y = -log10(DE_padj), color = conclusion)) + 
  theme_classic() + 
  geom_point(size = 1) + 
  xlab("L2FC(ksyn)") + 
  ylab("-log10(padj)") + 
  scale_color_manual(values = c("blue", "orange", "gray"))



########################################################################################################
ksyn_df <- data.frame(L2FC = reso$log2FoldChange + Fit_s$Fast_Fit$Effects_df$L2FC_kdeg,
                      Gene = Fit_s$Fast_Fit$Effects_df$XF)

ksyn_df$se <- sqrt(reso$lfcSE^2 + (Fit_s$Fast_Fit$Effects_df$se*log2(exp(1)))^2 ) 

# Calculate p-value using asymptotic Wald test
ksyn_df <- ksyn_df %>%
  mutate(pval = 2*pnorm(-abs(L2FC/se)),
         padj = p.adjust(pval, method = "BH"))

# Add conclusion at 0.01 FDR control
ksyn_df <- ksyn_df %>%
  mutate(conclusion = as.factor(ifelse(padj < 0.05, 
                                       ifelse(L2FC < 0, "Decreased txn", "Increased txn"),
                                       "Not sig.")))

write_delim(ksyn_df, "12042023_ksyn_df_4SA_4SD.csv")

# Make volcano plot
ggplot(ksyn_df, aes(x = L2FC, y = -log10(padj), color = conclusion)) + 
  theme_classic() + 
  geom_point(size = 1) + 
  xlab("L2FC(ksyn)") + 
  ylab("-log10(padj)") + 
  scale_color_manual(values = c("blue", "orange", "gray"))

