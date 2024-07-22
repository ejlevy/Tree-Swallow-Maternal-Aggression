### INCUBATION MS SCRIPT, SPRING 2024
# DESEQ2 on treatment, aggression scores, testosterone, capture latency, mass, and incubation day.
# I created dds files that tested the effect of each variable on its own (as opposed to a design with a few effects)



# Get data together, load libraries ----------------------------------------
library('DESeq2')
library(readxl)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = FALSE);

## Pull in raw counts
vSBNraw <- read_xlsx("rawCounts.vSBN.xlsx", sheet = 1) # 14782 x 21
vSBNraw <- vSBNraw[,-2]# remove F21-vSBN from vSBN dataset, she is a massive outlier (and the only outlier based on correlation and PCA analysis)

dSBNraw <- read_xlsx("rawCounts.dSBN.xlsx", sheet = 1) # 14782 x 21


## Get dataset in order - make rownames the transcript names, and load treatment data
vSBNraw2 <- vSBNraw[,-1]
rownames(vSBNraw2) <- vSBNraw$CONTIG
dim(vSBNraw2) # 14782    19

dSBNraw2 <- dSBNraw[,-1]
rownames(dSBNraw2) <- dSBNraw$CONTIG
dim(dSBNraw2) # 14782    20


## And versions with only the STI birds
# These are the STI birds (but if you are running you should double-check, in case things got shifted around!)
  # STI for vSBN: PATSF23	PATSF25	PATSF26	PATSF29	PATSF30	PATSF34	PATSF36	PATSF39	PATSF40
vSBNraw2_STIonly <- vSBNraw2[,c(2, 4, 5, 8, 9, 13, 15, 18, 19)] 
rownames(vSBNraw2_STIonly) <- vSBNraw$CONTIG

# These are the STI birds  (but if you are running you should double-check, in case things got shifted around!)
  # STI for dSBN: PATSF21	PATSF23	PATSF25	PATSF26	PATSF29	PATSF30	PATSF34	PATSF36	PATSF39	PATSF40
dSBNraw2_STIonly <- dSBNraw2[,c(1, 3, 5, 6, 9, 10, 14, 16, 19, 20)]
rownames(dSBNraw2_STIonly) <- dSBNraw$CONTIG




## Traits
# dSBN includes all birds; vSBN excludes F21
allTraits_dSBN <- read_xlsx("Full trait data.xlsx", sheet = 1) # 20 x 21
# Make a character column for treatment, just to make interpretation a little easier.
allTraits_dSBN$treatment <- ifelse(allTraits_dSBN$treatment_0con_1sti == 0, "Control", "STI")
# Remove PATSF21 from the vSBN version
allTraits_vSBN <- allTraits_dSBN[-1,]

# Make rownames the ID of the bird
allTraits_dSBN2 <- allTraits_dSBN[,-1]
rownames(allTraits_dSBN2) <- allTraits_dSBN$id

allTraits_vSBN2 <- allTraits_vSBN[,-1]
rownames(allTraits_vSBN2) <- allTraits_vSBN$id


## Trait data for STI_only birds
allTraits_dSBN_STI <- allTraits_dSBN2 %>% filter(treatment == "STI")
allTraits_vSBN_STI <- allTraits_vSBN2 %>% filter(treatment == "STI")



## Transcript-to-gene info, for creating helpful datasets later on.
info <- read_xlsx("Gene Info_Incubation TRES Aggression.xlsx", sheet = 1)

# Run Deseq2 --------------------------------------------------------------


#### Treatment deseq2 ####

# Create dds of Gene Expression ~ Treatment
ddsvSBN_treatment <- DESeqDataSetFromMatrix(countData = vSBNraw2 %>% 
                                    as.matrix(),
                                  colData = allTraits_vSBN2,
                                  design = ~ treatment)
ddsdSBN_treatment <- DESeqDataSetFromMatrix(countData = dSBNraw2 %>% 
                                   as.matrix(),
                                 colData = allTraits_dSBN2,
                                 design = ~ treatment)
# check data
ddsvSBN_treatment # 14782 19 
ddsdSBN_treatment  # 14782 20 

# Control should be the reference anyway, but let's make extra sure:
ddsvSBN_treatment$treatment<-relevel(ddsvSBN_treatment$treatment, ref = "Control")
ddsdSBN_treatment$treatment<-relevel(ddsdSBN_treatment$treatment, ref = "Control")

# Run Deseq
ddssvSBN_treatment<-DESeq(ddsvSBN_treatment)
ddssdSBN_treatment<-DESeq(ddsdSBN_treatment)

# Specify which results we want to look at
results_vSBN_treatment<-results(ddssvSBN_treatment, alpha = 0.1, contrast = c("treatment", "STI", "Control"))
results_dSBN_treatment<-results(ddssdSBN_treatment, alpha = 0.1, contrast = c("treatment", "STI", "Control"))

# Take a look at results
results_vSBN_treatment
results_dSBN_treatment

# And the summary
summary(results_vSBN_treatment)
summary(results_dSBN_treatment)

# This part is only helpful if there actually are DEGs...
subset(results_vSBN_treatment, padj < 0.1)
subset(results_dSBN_treatment, padj < 0.1)

# Plot
DESeq2::plotMA(results_vSBN_treatment, main = "Log2-Fold Change vs. Base Mean", ylim=c(-5,5))
DESeq2::plotMA(results_dSBN_treatment, main = "Log2-Fold Change vs. Base Mean", ylim=c(-5,5))

# Merge in the gene info for each transcript.
results_vSBN_treatment$TranscriptID <- rownames(results_vSBN_treatment) # Currently the row names are the transcript ID; making a new column for that
vSBN_treatment_merge <- merge(as.data.frame(results_vSBN_treatment), info, by = "TranscriptID", all.x = T)

results_dSBN_treatment$TranscriptID <- rownames(results_dSBN_treatment) # Currently the row names are the transcript ID; making a new column for that
dSBN_treatment_merge <- merge(as.data.frame(results_dSBN_treatment), info, by = "TranscriptID", all.x = T)

# Write data frame to csv
write.csv(vSBN_treatment_merge[order(vSBN_treatment_merge$pvalue),], file = "vSBN deseq by treatment.csv")
write.csv(dSBN_treatment_merge[order(dSBN_treatment_merge$pvalue),], file = "dSBN deseq by treatment.csv")




#### 5 min aggscore deseq2 ####
# Running a new design that uses STI-only birds and asks whether agg score predicts gene expression


## Create dds
ddsvSBNSTI5 <- DESeqDataSetFromMatrix(countData = vSBNraw2_STIonly %>% 
                                    as.matrix(),
                                  colData = allTraits_vSBN_STI,
                                  design = ~ agg5min)
ddsdSBNSTI5 <- DESeqDataSetFromMatrix(countData = dSBNraw2_STIonly %>% 
                                   as.matrix(),
                                 colData = allTraits_dSBN_STI,
                                 design = ~ agg5min)
# We get a warning when not scaling the variables; 
  # however, I ran it with scaled scores and compared the values; the p-values are essentially identical
  # r = 1.000000000 in the genes with adjusted p < 0.1. So I will continue with the raw values to make interpretation easier.
    # The warning (first part is not about the scaling): 
    # "the design formula contains one or more numeric variables that have mean or
    # standard deviation larger than 5 (an arbitrary threshold to trigger this message).
    # Including numeric variables with large mean can induce collinearity with the intercept.
    # Users should center and scale numeric variables in the design to improve GLM convergence."
  
ddssvSBNSTI5<-DESeq(ddsvSBNSTI5)
ddssdSBNSTI5<-DESeq(ddsdSBNSTI5)

results_vSBNSTI5<-results(ddssvSBNSTI5, alpha = 0.1, name = c("agg5min"))
results_dSBNSTI5<-results(ddssdSBNSTI5, alpha = 0.1, name = c("agg5min"))

results_vSBNSTI5
results_dSBNSTI5

summary(results_vSBNSTI5) # LFC up = 221 LFC down = 195. outliers [1]: 0, 0%; low counts [2]: 1699, 12%; (mean count < 6)
summary(results_dSBNSTI5) # LFC up = 18. LFC down = 0. outliers [1]  : 0, 0%; low counts [2]: 0, 0%; (mean count < 0)

subset(results_vSBNSTI5, padj < 0.1)
subset(results_dSBNSTI5, padj < 0.1)

DESeq2::plotMA(results_vSBNSTI5, main = "Log2-Fold Change vs. Base Mean", ylim=c(-.25,.25))
DESeq2::plotMA(results_dSBNSTI5, main = "Log2-Fold Change vs. Base Mean", ylim=c(-.25,.25)) # All the dSBN degs are higher in STI than control

results_vSBNSTI5$TranscriptID <- rownames(results_vSBNSTI5) # Currently the row names are the transcript ID; making a new column for that
vSBN_STI5_merge <- merge(as.data.frame(results_vSBNSTI5), info, by = "TranscriptID", all.x = T)

results_dSBNSTI5$TranscriptID <- rownames(results_dSBNSTI5) # Currently the row names are the transcript ID; making a new column for that
dSBN_STI5_merge <- merge(as.data.frame(results_dSBNSTI5), info, by = "TranscriptID", all.x = T)

write.csv(vSBN_STI5_merge[order(vSBN_STI5_merge$pvalue),], file = "vSBN deseq by 5min agg score.csv")
write.csv(dSBN_STI5_merge[order(dSBN_STI5_merge$pvalue),], file = "dSBN deseq by 5min agg score.csv")

# I also tested this out with a few dummy genes, it worked great.

# Pulling normcounts from vSBN to use for figure 4
normcounts_vSBN1 <- estimateSizeFactors(ddssvSBNSTI5)
normcounts_vSBN2 <- counts(normcounts_vSBN1, normalized = T)
write.csv(normcounts_vSBN2, "NormCounts vSBN 5min agg Deseq2.csv")


#### 30 min aggscore deseq2 ####


ddsvSBNSTI30 <- DESeqDataSetFromMatrix(countData = vSBNraw2_STIonly %>% 
                                       as.matrix(),
                                     colData = allTraits_vSBN_STI,
                                     design = ~ agg30min)
ddsdSBNSTI30 <- DESeqDataSetFromMatrix(countData = dSBNraw2_STIonly %>% 
                                      as.matrix(),
                                    colData = allTraits_dSBN_STI,
                                    design = ~ agg30min)

ddssvSBNSTI30<-DESeq(ddsvSBNSTI30)
ddssdSBNSTI30<-DESeq(ddsdSBNSTI30)

results_vSBNSTI30<-results(ddssvSBNSTI30, alpha = 0.1, name = c("agg30min"))
results_dSBNSTI30<-results(ddssdSBNSTI30, alpha = 0.1, name = c("agg30min"))

results_vSBNSTI30
results_dSBNSTI30

summary(results_vSBNSTI30) # LFC up = 214 LFC down = 32. outliers [1]: 0, 0%; low counts [2]: 2548, 17%; (mean count < 14)
summary(results_dSBNSTI30) # LFC up = 0 LFC down = 0.outliers [1]  : 0, 0%; low counts [2]: 0, 0%; (mean count < 0)

subset(results_vSBNSTI30, padj < 0.1)
subset(results_dSBNSTI30, padj < 0.1)

DESeq2::plotMA(results_vSBNSTI30, main = "Log2-Fold Change vs. Base Mean", ylim=c(-.25,.25))
DESeq2::plotMA(results_dSBNSTI30, main = "Log2-Fold Change vs. Base Mean", ylim=c(-.5,.5)) # All the dSBN degs are higher in STI than control

results_vSBNSTI30$TranscriptID <- rownames(results_vSBNSTI30) # Currently the row names are the transcript ID; making a new column for that
vSBN_STI30_merge <- merge(as.data.frame(results_vSBNSTI30), info, by = "TranscriptID", all.x = T)

results_dSBNSTI30$TranscriptID <- rownames(results_dSBNSTI30) # Currently the row names are the transcript ID; making a new column for that
dSBN_STI30_merge <- merge(as.data.frame(results_dSBNSTI30), info, by = "TranscriptID", all.x = T)

write.csv(vSBN_STI30_merge[order(vSBN_STI30_merge$pvalue),], file = "vSBN deseq by 30min agg score.csv")
write.csv(dSBN_STI30_merge[order(dSBN_STI30_merge$pvalue),], file = "dSBN deseq by 30min agg score.csv")



#### Testosterone deseq2 ####

## Create dds
ddsvSBNT <- DESeqDataSetFromMatrix(countData = vSBNraw2 %>% 
                                     as.matrix(),
                                   colData = allTraits_vSBN2,
                                   design = ~ testosterone)
ddsdSBNT <- DESeqDataSetFromMatrix(countData = dSBNraw2 %>% 
                                    as.matrix(),
                                  colData = allTraits_dSBN2,
                                  design = ~ testosterone)

ddssvSBNT<-DESeq(ddsvSBNT)
ddssdSBNT<-DESeq(ddsdSBNT)

results_vSBNT<-results(ddssvSBNT, alpha = 0.1, name = c("testosterone"))
results_dSBNT<-results(ddssdSBNT, alpha = 0.1, name = c("testosterone"))

results_vSBNT
results_dSBNT

summary(results_vSBNT) # LFC up = 0 LFC down = 1. outliers [1]: 10, 0.068%; low counts [2]: 0, 0%; (mean count < 0)
summary(results_dSBNT) # LFC up = 0 LFC down = 1. outliers [1]  : 1, 0.0068%; low counts [2]: 0, 0%; (mean count < 0)

subset(results_vSBNT, padj < 0.1)
subset(results_dSBNT, padj < 0.1)
# It's the same gene in both brain regions!! It's uncharacterized, no good hits in blast. XP_014129665.1

DESeq2::plotMA(results_vSBNT, main = "Log2-Fold Change vs. Base Mean", ylim=c(-0.025,0.025))
DESeq2::plotMA(results_dSBNT, main = "Log2-Fold Change vs. Base Mean", ylim=c(-0.025,0.025)) # All the dSBN degs are higher in STI than control

results_vSBNT$TranscriptID <- rownames(results_vSBNT) # Currently the row names are the transcript ID; making a new column for that
vSBN_T_merge <- merge(as.data.frame(results_vSBNT), info, by = "TranscriptID", all.x = T)

results_dSBNT$TranscriptID <- rownames(results_dSBNT) # Currently the row names are the transcript ID; making a new column for that
dSBN_T_merge <- merge(as.data.frame(results_dSBNT), info, by = "TranscriptID", all.x = T)

write.csv(vSBN_T_merge[order(vSBN_T_merge$pvalue),], file = "vSBN deseq by T.csv")
write.csv(dSBN_T_merge[order(dSBN_T_merge$pvalue),], file = "dSBN deseq by T.csv")




#### Capture latency deseq2 ####

## Create dds
ddsvSBNLat <- DESeqDataSetFromMatrix(countData = vSBNraw2_STIonly %>% 
                                       as.matrix(),
                                     colData = allTraits_vSBN_STI,
                                     design = ~ sti_start_to_cap_min)
ddsdSBNLat <- DESeqDataSetFromMatrix(countData = dSBNraw2_STIonly %>% 
                                      as.matrix(),
                                    colData = allTraits_dSBN_STI,
                                    design = ~ sti_start_to_cap_min)

ddssvSBNLat<-DESeq(ddsvSBNLat)
ddssdSBNLat<-DESeq(ddsdSBNLat)

results_vSBNLat<-results(ddssvSBNLat, alpha = 0.1, name = c("sti_start_to_cap_min"))
results_dSBNLat<-results(ddssdSBNLat, alpha = 0.1, name = c("sti_start_to_cap_min"))

results_vSBNLat
results_dSBNLat

summary(results_vSBNLat) # LFC up = 0 LFC down = 2. outliers [1]: 0, 0%; low counts [2]: 0, 0%; (mean count < 0)
summary(results_dSBNLat) # LFC up = 1 LFC down = 0. outliers [1]  : 0, 0%; low counts [2]: 0, 0%; (mean count < 0)

subset(results_vSBNLat, padj < 0.1)
subset(results_dSBNLat, padj < 0.1)

DESeq2::plotMA(results_vSBNLat, main = "Log2-Fold Change vs. Base Mean", ylim=c(-.25,.25))
DESeq2::plotMA(results_dSBNLat, main = "Log2-Fold Change vs. Base Mean", ylim=c(-.25,.25)) # All the dSBN degs are higher in STI than control

results_vSBNLat$TranscriptID <- rownames(results_vSBNLat) # Currently the row names are the transcript ID; making a new column for that
vSBN_Lat_merge <- merge(as.data.frame(results_vSBNLat), info, by = "TranscriptID", all.x = T)

results_dSBNLat$TranscriptID <- rownames(results_dSBNLat) # Currently the row names are the transcript ID; making a new column for that
dSBN_Lat_merge <- merge(as.data.frame(results_dSBNLat), info, by = "TranscriptID", all.x = T)

write.csv(vSBN_Lat_merge[order(vSBN_Lat_merge$pvalue),], file = "vSBN deseq by Sti-to-cap latency.csv")
write.csv(dSBN_Lat_merge[order(dSBN_Lat_merge$pvalue),], file = "dSBN deseq by Sti-to-cap latency.csv")



#### Mass deseq2 ####

## Create dds
ddsvSBNMass <- DESeqDataSetFromMatrix(countData = vSBNraw2 %>% 
                                       as.matrix(),
                                     colData = allTraits_vSBN2,
                                     design = ~ mass)
ddsdSBNMass <- DESeqDataSetFromMatrix(countData = dSBNraw2 %>% 
                                      as.matrix(),
                                    colData = allTraits_dSBN2,
                                    design = ~ mass)

ddssvSBNMass<-DESeq(ddsvSBNMass)
ddssdSBNMass<-DESeq(ddsdSBNMass)

results_vSBNMass<-results(ddssvSBNMass, alpha = 0.1, name = c("mass"))
results_dSBNMass<-results(ddssdSBNMass, alpha = 0.1, name = c("mass"))

results_vSBNMass
results_dSBNMass

summary(results_vSBNMass) # LFC up = 1 LFC down = 1. outliers [1]: 0, 0%; low counts [2]: 0, 0%; (mean count < 0)
summary(results_dSBNMass) # LFC up = 0 LFC down = 0. outliers [1]  : 0, 0%; low counts [2]: 0, 0%; (mean count < 0)

subset(results_vSBNMass, padj < 0.1)
subset(results_dSBNMass, padj < 0.1)

DESeq2::plotMA(results_vSBNMass, main = "Log2-Fold Change vs. Base Mean", ylim=c(-2,2))
DESeq2::plotMA(results_dSBNMass, main = "Log2-Fold Change vs. Base Mean", ylim=c(-2,2)) # All the dSBN degs are higher in STI than control

results_vSBNMass$TranscriptID <- rownames(results_vSBNMass) # Currently the row names are the transcript ID; making a new column for that
vSBN_Mass_merge <- merge(as.data.frame(results_vSBNMass), info, by = "TranscriptID", all.x = T)

results_dSBNMass$TranscriptID <- rownames(results_dSBNMass) # Currently the row names are the transcript ID; making a new column for that
dSBN_Mass_merge <- merge(as.data.frame(results_dSBNMass), info, by = "TranscriptID", all.x = T)

write.csv(vSBN_Mass_merge[order(vSBN_Mass_merge$pvalue),], file = "vSBN deseq by Mass.csv")
write.csv(dSBN_Mass_merge[order(dSBN_Mass_merge$pvalue),], file = "dSBN deseq by Mass.csv")




#### Incubation day deseq2 ####

ddsvSBNInc <- DESeqDataSetFromMatrix(countData = vSBNraw2 %>% 
                                     as.matrix(),
                                   colData = allTraits_vSBN2,
                                   design = ~ incubation_day_EJL_estimates)
ddsdSBNInc <- DESeqDataSetFromMatrix(countData = dSBNraw2 %>% 
                                    as.matrix(),
                                  colData = allTraits_dSBN2,
                                  design = ~ incubation_day_EJL_estimates)

ddssvSBNInc<-DESeq(ddsvSBNInc)
ddssdSBNInc<-DESeq(ddsdSBNInc)

results_vSBNInc<-results(ddssvSBNInc, alpha = 0.1, name = c("incubation_day_EJL_estimates"))
results_dSBNInc<-results(ddssdSBNInc, alpha = 0.1, name = c("incubation_day_EJL_estimates"))

results_vSBNInc
results_dSBNInc

summary(results_vSBNInc) # LFC up = 1 LFC down = 0. outliers [1]: 42, 0.29%; low counts [2]: 0, 0%; (mean count < 0)
summary(results_dSBNInc) # LFC up = 0 LFC down = 0. outliers [1]  : 46, 0.31%; low counts [2]: 0, 0%; (mean count < 0)

subset(results_vSBNInc, padj < 0.1)
subset(results_dSBNInc, padj < 0.1)

DESeq2::plotMA(results_vSBNInc, main = "Log2-Fold Change vs. Base Mean", ylim=c(-2,2))
DESeq2::plotMA(results_dSBNInc, main = "Log2-Fold Change vs. Base Mean", ylim=c(-2,2)) 

results_vSBNInc$TranscriptID <- rownames(results_vSBNInc) # Currently the row names are the transcript ID; making a new column for that
vSBN_Inc_merge <- merge(as.data.frame(results_vSBNInc), info, by = "TranscriptID", all.x = T)

results_dSBNInc$TranscriptID <- rownames(results_dSBNInc) # Currently the row names are the transcript ID; making a new column for that
dSBN_Inc_merge <- merge(as.data.frame(results_dSBNInc), info, by = "TranscriptID", all.x = T)

write.csv(vSBN_Inc_merge[order(vSBN_Inc_merge$pvalue),], file = "vSBN deseq by Incubation day.csv")
write.csv(dSBN_Inc_merge[order(dSBN_Inc_merge$pvalue),], file = "dSBN deseq by Incubation day.csv")