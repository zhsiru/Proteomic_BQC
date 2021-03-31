setwd("C:\\Users\\sirui.zhou\\work\\Works-home\\Covid\\SomaLogic")
setwd("C:\\Users\\Sirui\\Desktop\\WORKS\\OAS1.NM.revision\\McGill-Richards-C-19 - SomaScan Data")
install.packages(c("readr", "stringr"))
install.packages("SomaDataIO_3.1.0.tar.gz", repos = NULL, type = "source")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install()
BiocManager::install(c("Biobase"))
BiocManager::install(c("limma"))
BiocManager::install(c("pvca"))
install.packages("factoextra")
install.packages("gbm")
install.packages("GGally")
install.packages("ggdendro")
install.packages("glmnet")
install.packages("pROC")
install.packages("randomForest")
install.packages("rmarkdown")
install.packages("rpart.plot")

library(Biobase)
library(limma)
library(pvca)
library(dplyr)      # data manipulation and pipe operation
library(factoextra) # extract and visualize results of multivariate data analysis
library(forestplot) # forest plot
library(gbm)        # generalized boosted regression models
library(GGally)     # extension to 'ggplot2'
library(ggdendro)   # dendrogram for data clustering
library(ggfortify)  # data visualization tools for statistical analysis results
library(ggplot2)    # create graphics based on "The Grammer of Graphics"
library(ggrepel)    # tidy text display in ggplot
library(glmnet)     # cross validation plot for glmnet
library(grDevices)  # R graphics devices and support for colors and fonts
library(gridExtra)  # Grid graphics
library(knitr)      # dynamic report generation
library(methods)    # formal methods and classes
library(pROC)       # display and analyze ROC curves
library(randomForest) # Random forest variable importance
library(reshape2)   # flexibly reshape data
library(rmarkdown)  # dynamic documents for R
library(rpart.plot) # plots for recursive partitioning for classification, regression and survival trees
library(tibble)     # simple data frames
library(stats)      # basic statistical functions

install.packages("statVisual")
library(statVisual)

library(SomaDataIO)
my_adat <- read.adat("SS-200150_v4_ACDPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.medNormRefSMP.adat")
dat <- read.adat("SS-200150_v4_ACDPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.adat")
is.soma_adat(my_adat)
my_adat
dim(my_adat)
print(my_adat, show_header = TRUE)
OAS1  <- grep("^OAS1", getFeatures(my_adat), value = TRUE)
#indicator_prptein  <- grep("^MAX | ^CETN2 | ^SOST | ^C43BP | ^CES1 | ^K6PF | ^NAA10 | ^CTF1 | ^LEG9 | ^SMAD1", getFeatures(my_adat), value = TRUE)
#VDR  <- grep("^VDR", getFeatures(my_adat), value = TRUE)
OAS1 <- dplyr::filter(getFeatureData(my_adat), AptName %in% OAS1)
library(dplyr)
OAS1.dat <- my_adat %>% select(PlateId, RowCheck, NormScale_20, NormScale_0_005, NormScale_0_5, SubjectID, OAS1.10361.25, matches('MAX|CETN2|SOST.13101.60|COL4A3BP|CES1|PFKM|NAA10|CTF1|LGALS9|SMAD1|USP8|SOD1|STIP1|PRDX3|APEX1|CRKL')) 
OAS1.dat.2 <- dat %>% select(PlateId, RowCheck, NormScale_20, NormScale_0_005, NormScale_0_5, SubjectID, OAS1.10361.25, matches('MAX|CETN2|SOST.13101.60|COL4A3BP|CES1|PFKM|NAA10|CTF1|LGALS9|SMAD1|USP8|SOD1|STIP1|PRDX3|APEX1|CRKL'))
#VDR.dat <- my_adat %>% select(SubjectID, VDR.10023.32)
#OAS1.dat.all <- OAS1.dat[!is.na(OAS1.dat$OAS1.10361.25), ] 
OAS1.dat.all.medNormRefSMP <- OAS1.dat[grepl("VAP|CHUM", OAS1.dat$SubjectID), , drop = FALSE]
OAS1.dat.all <- OAS1.dat.2[grepl("VAP|CHUM", OAS1.dat$SubjectID), , drop = FALSE]
#VDR.dat.JGH <- VDR.dat[grepl("VAP", VDR.dat$SubjectID), , drop = FALSE]
All <- read.table("Sample_full.txt", header=T)
All_OAS1_merge_medNormRefSMP <- merge(OAS1.dat.all.medNormRefSMP , All, by="SubjectID")
All_OAS1_merge <- merge(OAS1.dat.all, All, by="SubjectID")

Sample <- readRDS("../basic_JGH_CHUM.rds")
#Sample$covid19_test_date <- as.factor(Sample$covid19_test_date)
#Sample$covid19_first_symptoms_date <- as.factor(Sample$covid19_first_symptoms_date)

All_OAS1_merge_medNormRefSMP_sampleINFO <- merge(All_OAS1_merge_medNormRefSMP , Sample, by="anonymized_patient_id")
All_OAS1_merge_sampleINFO <- merge(All_OAS1_merge , Sample, by="anonymized_patient_id") 

write.table(All_OAS1_merge_medNormRefSMP_sampleINFO, file="All_OAS1_merge_medNormRefSMP_sampleINFO.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(All_OAS1_merge_sampleINFO, file="All_OAS1_merge_sampleINFO.txt",col.names=T,row.names=F,quote=F,sep="\t")

indicator.EDTA.pca <- prcomp(All_OAS1_merge_medNormRefSMP_sampleINFO[,c(10,12:14,16:18,20,21,24,25)], center = TRUE,scale. = TRUE)
indicator.citrated.pca <- prcomp(All_OAS1_merge_medNormRefSMP_sampleINFO[,c(9:11,13:17,19,22,23)], center = TRUE,scale. = TRUE)

indicator.EDTA.pca2 <- prcomp(All_OAS1_merge_sampleINFO[,c(10,12:14,16:18,20,21,24,25)], center = TRUE,scale. = TRUE)
indicator.citrated.pca2 <- prcomp(All_OAS1_merge_sampleINFO[,c(9:11,13:17,19,22,23)], center = TRUE,scale. = TRUE)

EDTA <- as.data.frame(indicator.EDTA.pca$x)
citrated <- as.data.frame(indicator.citrated.pca$x)

All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_EDTA_pca <- cbind(All_OAS1_merge_medNormRefSMP_sampleINFO, EDTA)
names(All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_EDTA_pca)[names(All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_EDTA_pca) == "PC1"] <- "EDTA.PC1"
names(All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_EDTA_pca)[names(All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_EDTA_pca) == "PC2"] <- "EDTA.PC2"
All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_EDTA_pca <- All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_EDTA_pca[1:(length(All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_EDTA_pca)-9)]
All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_pca <- cbind(All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_EDTA_pca, citrated)
names(All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_pca)[names(All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_pca) == "PC1"] <- "citrated.PC1"
names(All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_pca)[names(All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_pca) == "PC2"] <- "citrated.PC2"
All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_pca <- All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_pca[1:(length(All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_pca)-9)]

Control_viral_infection <- read.table("Control_viral_infection_505.txt",header=T)
All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_pca <- merge(All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_pca, Control_viral_infection, by="SubjectID")
All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_pca <- All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_pca[-c(9:25)]
OAS1_dat <- All_OAS1_merge_medNormRefSMP_sampleINFO_indicator_pca
write.table(OAS1_dat, file="All_OAS1_merge_medNormRefSMP_sampleINFO.txt",col.names=T,row.names=F,quote=F,sep="\t")

OAS1_dat <- read.table("All_OAS1_merge_medNormRefSMP_sampleINFO.txt", header=T)
sym.update <- read.table("symp.update.505.txt",header=T)
Pro_hours <- read.table("process_time.txt", header = T)

OAS1_dat <- merge(OAS1_dat, sym.update, by = "SubjectID")
OAS1_dat <- merge(OAS1_dat, Pro_hours, by = "SubjectID")
OAS1_dat <- OAS1_dat[!is.na(OAS1_dat$OAS1.10361.25),]
OAS1_dat$LOG_OAS1 <- log(OAS1_dat$OAS1.10361.25)

write.table(OAS1_dat, file="All_OAS1_merge_medNormRefSMP_sampleINFO_update_505.txt",col.names=T,row.names=F,quote=F,sep="\t")
OAS1_dat <- read.table("All_OAS1_merge_medNormRefSMP_sampleINFO_update.txt", header=T)

OAS1_no_norm <- dat %>% select(SubjectID, OAS1.10361.25)

OAS1_dat_withnonorm <- merge(OAS1_dat, OAS1_no_norm, by="SubjectID")

library(ggpubr)
par(mar=c(5,8,8,8))
ggscatter(OAS1_dat, x = "citrated.PC1", y = "citrated.PC2", 
          color = "SampleGroup", size = 3, alpha = 0.6, 
          label.select = list(criteria = "`y` > 10 | `x` > 11"), 
          label="SubjectID", 
          font.label = c(9, "bold", "#756bb1"), 
          repel = TRUE,
          title = "citrated indicator protein PCs for different sites")

D0 <- OAS1_dat[ which(OAS1_dat$Draw == "D0"), ]
fit = glm(C1 ~  LOG_OAS1 + age_at_diagnosis + sex + SampleGroup + citrated.PC1, family=binomial(), data = D0)
summary(fit)
nrow(D0[which(D0$C1 == "0"), ])
nrow(D0[which(D0$C1 == "1"), ])
D0$SampleGroup <- gsub('JGH', '1', D0$SampleGroup)
D0$SampleGroup <- gsub('CHUM', '0', D0$SampleGroup)
fit <- glm(as.numeric(SampleGroup) ~  LOG_OAS1 + age_at_diagnosis, data = D0)
summary(fit)

D30_nov <- OAS1_dat[ which(OAS1_dat$Draw == "D30" | OAS1_dat$Controls_without_viral_infection_stringent == "1"), ]
fit = glm(A1 ~  LOG_OAS1 + age_at_diagnosis + sex , family=binomial(), data = D30)
summary(fit)
nrow(D30[which(D30$C1 == "0"), ])
nrow(D30[which(D30$C1 == "1"), ])

D30 <- OAS1_dat[ which(OAS1_dat$Draw == "D30" | OAS1_dat$Draw == "D0" & OAS1_dat$CaseUpdate == "Control"), ]
fit = glm(C1 ~  LOG_OAS1 + age_at_diagnosis + sex + SampleGroup + citrated.PC1, family=binomial(), data = D30)
summary(fit)
nrow(D30[which(D30$C1 == "0"), ])
nrow(D30[which(D30$C1 == "1"), ])



###test Celia's model
install.packages("mgcv")
library(mgcv)
library(dplyr)
OAS1_dat <- read.table("All_OAS1_merge_medNormRefSMP_sampleINFO_update.txt", header=T)
OAS1_dat <- OAS1_dat[order(OAS1_dat$Days_symptom_update3),]
OAS1_dat$T1 = OAS1_dat$Days_symptom_update3
OAS1_dat <- within(OAS1_dat, T1[CaseUpdate == 'Control'] <- '0')
OAS1_dat_1st <- OAS1_dat %>% 
  group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == first(OAS1.10361.25))
OAS1_dat_1st$LOG_OAS <- log(OAS1_dat_1st$OAS1.10361.25)
OAS1_dat_1st$OAS_norm <- scale(residuals(glm(LOG_OAS ~ ProcessTime.hrs., data=OAS1_dat_1st)))
OAS1_dat_1st$ix1 <- as.numeric(OAS1_dat_1st$T1)*OAS1_dat_1st$OAS_norm
OAS1_dat_1st$ix2 <- as.numeric(OAS1_dat_1st$T1)*OAS1_dat_1st$LOG_OAS
fit2 <- gam(A2 ~ sex + s(age_at_diagnosis) + PlateId + s(LOG_OAS) + s(T1, by = LOG_OAS) , family=binomial(), data = OAS1_dat_1st[which(OAS1_dat_1st$ix2 <200),])
summary(fit2)
plot(fit2,select=3,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=5,xlim = c(0,150), ylim = c(-10, 10))
abline(h=0)

test <- OAS1_dat_1st[which(OAS1_dat_1st$ix2 <200),]
test$T1 <- as.numeric(test$T1)
OAS1_dat_1st$T1 <- as.numeric(OAS1_dat_1st$T1)
OAS1_dat_last$T1 <- as.numeric(OAS1_dat_last$T1)

fit2 <- gam(A2 ~ sex + s(age_at_diagnosis) + PlateId + s(LOG_OAS) + s(ix2) , family=binomial(), data = test)
fit2 <- gam(A2 ~ sex + s(age_at_diagnosis) + PlateId + s(T1, by = LOG_OAS) , family=binomial(), data = OAS1_dat_last)
fit2 <- gam(A2 ~ sex + s(age_at_diagnosis) + PlateId + s(ix2) , family=binomial(), data = OAS1_dat_1st)
summary(fit2)
plot(fit2,select=2,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=5, ylim = c(-1, 1))
abline(h=0)






OAS1_dat_last <- OAS1_dat %>% 
  group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == last(OAS1.10361.25))
OAS1_dat_last$LOG_OAS <- log(OAS1_dat_last$OAS1.10361.25)
OAS1_dat_last$OAS_norm <- scale(residuals(glm(LOG_OAS ~ ProcessTime.hrs., data=OAS1_dat_last)))
OAS1_dat_last$ix1 <- as.numeric(OAS1_dat_last$T1)*OAS1_dat_last$OAS_norm
fit2 <- gam(A2 ~ sex + s(age_at_diagnosis) + PlateId + s(OAS_norm) + s(ix1) , family=binomial(), data = OAS1_dat_last)
summary(fit2)



OAS1_dat$LOG_OAS <- log(OAS1_dat$OAS1.10361.25)
OAS1_dat$OAS_norm <- scale(residuals(glm(LOG_OAS ~ ProcessTime.hrs., data=OAS1_dat)))
OAS1_dat$ix1 <- as.numeric(OAS1_dat$T1)*OAS1_dat$OAS_norm
fit3 <- gam(A2 ~ sex + s(age_at_diagnosis) + PlateId + s(OAS_norm) + s(ix1) , family=binomial(), data = OAS1_dat)
summary(fit3)
plot(fit3,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=5)





###Covid infectious: 1-14 days of sym onset, >31 days non-infectious

OAS1_dat <- read.table("All_OAS1_merge_medNormRefSMP_sampleINFO_update.txt", header=T)
OAS1_dat$age2 <- OAS1_dat$age_at_diagnosis^2
OAS1_dat <- OAS1_dat[order(OAS1_dat$Days_symptom_update3),] 

INFE <- OAS1_dat[ which(OAS1_dat$Days_symptom_update3 < 15 | OAS1_dat$CaseUpdate == "Control"), ]
Non_INFE_31 <- OAS1_dat[ which(OAS1_dat$Days_symptom_update3 > 31 | OAS1_dat$CaseUpdate == "Control"), ]
Non_INFE_31_rm <- OAS1_dat[ which(OAS1_dat$Days_symptom_update3 > 31 | OAS1_dat$CaseUpdate == "Control" & OAS1_dat$Controls_without_viral_infection_stringent == "0"), ]


#Non_INFE <- Non_INFE %>% group_by(anonymized_patient_id) %>% mutate(OAS_mean = mean(OAS1.10361.25, na.rm = TRUE))
##mean, first, latest
INFE <- INFE %>% 
  group_by(anonymized_patient_id) %>% 
  mutate(OAS_mean = mean(OAS1.10361.25, na.rm = TRUE),
         ProcessTime_mean = mean(ProcessTime.hrs., na.rm = TRUE)) %>% 
  add_count(anonymized_patient_id)

Non_INFE_31 <- Non_INFE_31 %>% 
  group_by(anonymized_patient_id) %>% 
  mutate(OAS_mean = mean(OAS1.10361.25, na.rm = TRUE),
         ProcessTime_mean = mean(ProcessTime.hrs., na.rm = TRUE)) %>% 
  add_count(anonymized_patient_id)


Non_INFE_31_1st <- Non_INFE_31 %>% 
  group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == first(OAS1.10361.25))

Non_INFE_31_last <- Non_INFE_31 %>% 
  group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == last(OAS1.10361.25))


INFE_1st <- INFE %>% 
  group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == first(OAS1.10361.25))

INFE_last <- INFE %>% 
  group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == last(OAS1.10361.25))

INFE_1st$LOG_OAS <- log(INFE_1st$OAS1.10361.25)
INFE_1st$OAS_norm <- scale(residuals(glm(LOG_OAS ~ ProcessTime.hrs., data=INFE_1st)))

INFE_last$LOG_OAS <- log(INFE_last$OAS1.10361.25)
INFE_last$OAS_norm <- scale(residuals(glm(LOG_OAS ~ ProcessTime.hrs., data=INFE_last)))


Non_INFE_31_1st$LOG_OAS <- log(Non_INFE_31_1st$OAS1.10361.25)
Non_INFE_31_1st$OAS_norm <- scale(residuals(glm(LOG_OAS ~ ProcessTime.hrs., data=Non_INFE_31_1st)))

Non_INFE_31_last$LOG_OAS <- log(Non_INFE_31_last$OAS1.10361.25)
Non_INFE_31_last$OAS_norm <- scale(residuals(glm(LOG_OAS ~ ProcessTime.hrs., data=Non_INFE_31_last)))


fit = glm(C1 ~  OAS_norm + age2 + sex + age_at_diagnosis + SampleGroup + PlateId, family=binomial(), data = Non_INFE_31_1st)
summary(fit)

fit = glm(C1 ~  OAS_norm + age2 + sex + age_at_diagnosis + SampleGroup + PlateId, family=binomial(), data = Non_INFE_31_last)
summary(fit)


fit = glm(C1 ~  OAS_norm + age2 + sex + age_at_diagnosis + SampleGroup + PlateId, family=binomial(), data = INFE_1st)
summary(fit)

fit = glm(C1 ~  OAS_norm + age2 + sex + age_at_diagnosis + SampleGroup + PlateId, family=binomial(), data = INFE_last)
summary(fit)


write.table(Non_INFE_31_last, file="Non_INFE_31_last_505.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(INFE_1st, file="INFE_1st_505.txt",col.names=T,row.names=F,quote=F,sep="\t")


A <- read.table("INFE_1st.txt", header=T)
B <- read.table("Non_INFE_31_last.txt", header=T)
A_nonorm <- merge(A, OAS1_no_norm, by="SubjectID")
B_nonorm <- merge(B, OAS1_no_norm, by="SubjectID")

A_nonorm$LOG_OAS_nonorm <- log(A_nonorm$OAS1.10361.25.y)
B_nonorm$LOG_OAS_nonorm <- log(B_nonorm$OAS1.10361.25.y)

hist(B_nonorm$LOG_OAS_nonorm, breaks = 100)
hist(B_nonorm$LOG_OAS1, breaks = 100)

A_QC1 <- A[which(A$ProcessTime.hrs. < 50 & A$LOG_OAS < 8), ]
B_QC1 <- B[which(B$ProcessTime.hrs. < 50 & B$LOG_OAS < 8), ]


A_nonorm_QC1 <- A_nonorm[which(A_nonorm$ProcessTime.hrs. < 50 & A_nonorm$LOG_OAS_nonorm < 7.5 & A_nonorm$LOG_OAS_nonorm > 5), ]
B_nonorm_QC1 <- B_nonorm[which(B_nonorm$ProcessTime.hrs. < 50 & B_nonorm$LOG_OAS_nonorm < 7 & B_nonorm$LOG_OAS_nonorm > 5), ]

write.table(A_QC1, file="Non_INFE_31_last.QC1.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(B_QC1, file="INFE_1st.QC1.txt",col.names=T,row.names=F,quote=F,sep="\t")


write.table(A_nonorm_QC1, file="Non_INFE_31_last.QC1.nonorm.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(B_nonorm_QC1, file="INFE_1st.QC1.nonorm.txt",col.names=T,row.names=F,quote=F,sep="\t")

A_QC1$LOG_OAS2 <- log(A_QC1$OAS1.10361.25)
A_QC1$OAS_norm2 <- scale(residuals(glm(LOG_OAS2 ~ ProcessTime.hrs., data=A_QC1)))
B_QC1$LOG_OAS2 <- log(B_QC1$OAS1.10361.25)
B_QC1$OAS_norm2 <- scale(residuals(glm(LOG_OAS2 ~ ProcessTime.hrs., data=B_QC1)))

A_nonorm_QC1$OAS_nomednorm_norm <- scale(residuals(glm(LOG_OAS_nonorm ~ ProcessTime.hrs., data=A_nonorm_QC1)))
B_nonorm_QC1$OAS_nomednorm_norm <- scale(residuals(glm(LOG_OAS_nonorm ~ ProcessTime.hrs., data=B_nonorm_QC1)))


fit = glm(A2 ~  OAS_norm2 + age2 + sex + age_at_diagnosis + SampleGroup + PlateId, family=binomial(), data = A_QC1)
summary(fit)

fit = glm(A2 ~  OAS_nomednorm_norm + age2 + sex + age_at_diagnosis + SampleGroup + PlateId, family=binomial(), data = A_nonorm_QC1)
summary(fit)
fit = glm(A2 ~  OAS_nomednorm_norm + age2 + sex + age_at_diagnosis + SampleGroup + PlateId, family=binomial(), data = B_nonorm_QC1)
summary(fit)



PCS <- read.table("BQC19.PC.eigenvec", header = T)
A <- merge(A, PCS, by = "anonymized_patient_id")
B <- merge(B, PCS, by = "anonymized_patient_id")



A_QC1$GROUP <- c("During Infection")
B_QC1$GROUP <- c("Non-infectious State")
BQC <- rbind(A_QC1, B_QC1)




t.test(OAS1.10361.25 ~ GROUP, data = BQC, alternative = c("two.sided"))


###new violin###

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(hrbrthemes)
library(viridis)

BQC_A2 <- BQC %>% 
  mutate(A2 = recode(A2, '0' = 'Control', '1' = 'Case')) %>%
  mutate(ICDA = "Very Severe COVID-19") %>%
  select(c(A2,ICDA,GROUP,LOG_OAS)) %>%
  rename(CC = A2)

BQC_B2 <- BQC %>% 
  mutate(B2 = recode(B2, '0' = 'Control', '1' = 'Case')) %>%
  mutate(ICDA = "Hospitalization") %>%
  select(c(B2,ICDA,GROUP,LOG_OAS)) %>%
  rename(CC = B2)

BQC_C1 <- BQC %>% 
  mutate(C1 = recode(C1, '0' = 'Control', '1' = 'Case')) %>%
  mutate(ICDA = "Susceptibility") %>%
  select(c(C1,ICDA,GROUP,LOG_OAS)) %>%
  rename(CC = C1)

BQC2 <- rbind(BQC_A2, BQC_B2, BQC_C1)

BQC2$ICDA = factor(BQC2$ICDA, levels=c("Very Severe COVID-19", "Hospitalization", "Susceptibility"))
p <- ggplot(aes(y = LOG_OAS,
           x = `CC`,
           fill = CC),
       data = BQC2) +
  #geom_violin(position=position_dodge(),draw_quantiles=c(0.5)) +
  geom_violin(width=0.7, position="dodge") +
  geom_boxplot(width=0.3,color="white", alpha=0.2, position = position_dodge(width =0.9))+
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum() +
  labs(title="",
       x="", 
       y="LOG 0AS1 Level") +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size=11)
  ) +
  scale_fill_brewer(palette="Set1") +
  facet_grid(`ICDA`~`GROUP`)

p + labs(fill = "Outcome")


###density



###After sample selection and QC
BQC$T1 = BQC$Days_symptom_update3
BQC <- within(BQC, T1[CaseUpdate == 'Control'] <- '0')
BQC$ProcessTime.hrs. <- as.numeric(BQC$ProcessTime.hrs.)
BQC$T1 <- as.numeric(BQC$T1)


library(ggExtra)

p <- ggplot(BQC, aes(x=T1, y=ProcessTime.hrs., color=LOG_OAS)) +
  geom_point(size = 3, alpha=0.8) +
  theme(legend.position=c(1.1, 1.1)) +
  labs(title="",
       x="Days Since Symptom Onset", 
       y="Sample Processing Time")

p2 <- ggMarginal(p, type="histogram", fill = "#3182bd", xparams = list(  bins=50))


p2


#####All samples

OAS1_dat <- read.table("All_OAS1_merge_medNormRefSMP_sampleINFO_update.txt", header=T)
OAS1_dat$age2 <- OAS1_dat$age_at_diagnosis^2
OAS1_dat <- OAS1_dat[order(OAS1_dat$Days_symptom_update3),] 

OAS1_dat$T1 = OAS1_dat$Days_symptom_update3
OAS1_dat <- within(OAS1_dat, T1[CaseUpdate == 'Control'] <- '0')
OAS1_dat$ProcessTime.hrs. <- as.numeric(OAS1_dat$ProcessTime.hrs.)
OAS1_dat$T1 <- as.numeric(OAS1_dat$T1)


p <- ggplot(OAS1_dat, aes(x=T1, y=ProcessTime.hrs., color=LOG_OAS1)) +
  geom_point(size = 3, alpha=0.8) +
  theme(legend.position=c(1.1, 1.1)) +
  labs(title="",
       x="Days Since Symptom Onset", 
       y="Sample Processing Time")

p2 <- ggMarginal(p, type="histogram", fill = "#3182bd", xparams = list(  bins=50))


p2


###simple density for all samples

library(ggplot2)
library(dplyr)
library(gridExtra)

p1 <- OAS1_dat %>%
  filter( CaseUpdate != 'Control' ) %>%
  ggplot( aes(x=T1)) +
  geom_density(fill="#6baed6", alpha=0.8) +
  labs(title="",
       x="Days Since Symptom Onset", 
       y="") +
  ggtitle("A. Days Since Symptom Onset for COVID-19 Patients")


p2 <- OAS1_dat %>%
  #filter( CaseUpdate != 'Control' ) %>%
  ggplot( aes(x=ProcessTime.hrs.)) +
  geom_density(fill="#74c476", alpha=0.8) +
  labs(title="",
       x="Sample Processing Time (hrs)", 
       y="") +
  ggtitle("B. Sample Processing Time for BQC19 Samples")


###density of samples remove outliers###
BQC$T1 = BQC$Days_symptom_update3
library(ggplot2)
library(dplyr)
p3 <- BQC %>%
  filter( CaseUpdate != 'Control' ) %>%
  ggplot( aes(x=T1)) +
  geom_density(fill="#6baed6", alpha=0.8) +
  labs(title="",
       x="Days Since Symptom Onset", 
       y="") +
  ggtitle("C. Days Since Symptom Onset for COVID-19 Patients")


p4 <- BQC %>%
  #filter( CaseUpdate != 'Control' ) %>%
  ggplot( aes(x=ProcessTime.hrs.)) +
  geom_density(fill="#74c476", alpha=0.8) +
  labs(title="",
       x="Sample Processing Time (hrs)", 
       y="") +
  ggtitle("D. Sample Processing Time for BQC19 Samples")

grid.arrange(p1, p2, p3, p4, nrow = 2)

####old violin##


BQC3 <- BQC2 %>%
  mutate(ICDA_CC = fct_reorder(ICDA_CC, LOG_OAS)) %>%
  mutate(ICDA_CC = factor(ICDA_CC, levels=c("A2 Cases", "A2 Controls", "B2 Cases", "B2 Controls", "C1 Cases", "C1 Controls")))
p <- ggplot(aes(fill=GROUP, y=LOG_OAS, x=ICDA_CC), data = BQC3) + 
  geom_violin(width=1, position="dodge") +
  geom_boxplot(width=1, color="dark grey", alpha=0.2) + 
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum() +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
    plot.title = element_text(size=11)
  ) +
  scale_fill_brewer(palette = "Set1") +
  xlab("ICDA Case and Controls") +
  ylab("OAS1 raw level") +
  ylim(5,7.5)

p


OAS1_dat <- read.table("All_OAS1_merge_medNormRefSMP_sampleINFO_update.txt", header=T)
OAS1_dat <- OAS1_dat[order(OAS1_dat$Days_symptom_update3),] 
OAS1_dat_QC1 <- OAS1_dat[which(OAS1_dat$ProcessTime.hrs. < 50 & OAS1_dat$LOG_OAS1 < 8), ]
OAS1_dat_QC1 <- dplyr::add_count(OAS1_dat_QC1, anonymized_patient_id)
OAS1_dat_4 <- OAS1_dat_QC1[ which(OAS1_dat_QC1$n == "4"), ]
OAS1_dat_4 <- OAS1_dat_4 %>% 
  mutate(A2 = recode(A2, '0' = 'Non severe COVID-19 patients', '1' = 'severe COVID-19 patients'))


#OAS1_dat_4_case <- OAS1_dat_4[ which(OAS1_dat_4$C1 == "1"), ]
p=statVisual(type = "LinePlot",
             data = OAS1_dat_4,
             x = 'Days_symptom_update3',
             y = 'OAS1.10361.25',
             sid = 'anonymized_patient_id',
             xlab = "Days since symptom onset",
             ylab = "OAS1 raw level", group = "A2") + labs(title = "OAS1 level trajectory with COVID19 severity") + labs(color='COVID-19 Outcomes') + theme_classic()
p + scale_x_discrete(breaks = seq(0, 140, by = 10)) + scale_color_manual(values=c('#4292c6','#d94801'))



library(statVisual)
OAS1_dat <- dplyr::add_count(OAS1_dat, anonymized_patient_id)
OAS1_dat_4 <- OAS1_dat[ which(OAS1_dat$n == "4"), ]
#OAS1_dat_4_case <- OAS1_dat_4[ which(OAS1_dat_4$C1 == "1"), ]
p=statVisual(type = "LinePlot",
             data = OAS1_dat_4,
             x = 'Days_symptom_update3',
             y = 'LOG_OAS1',
             sid = 'anonymized_patient_id',
             xlab = "Days since symptom onset",
             ylab = "Log OAS1 level", group = "SampleGroup") + labs(title = "OAS1 level trajectory in COVID19 patients") + labs(color='COVID-19 Outcomes') + theme_classic()
p + scale_x_discrete(breaks = seq(0, 140, by = 10))


library(ggpubr)
B_rev <- B %>% 
  mutate(C1 = recode(C1, '0' = 'PCR negative', '1' = 'PCR positive'))
ggscatter(B_rev, x = "age_at_diagnosis", y = "LOG_OAS", 
          add = "reg.line",
          conf.int = TRUE, color = "C1", size = 3, alpha = 0.6, 
          cor.coef = T, cor.method = "pearson",
          title = "OAS1 change with Age in Late Stage",
          xlab = "Age", ylab = "Log OAS1 level", repel = TRUE)

ggscatter(B_QC1, x = "age_at_diagnosis", y = "LOG_OAS", 
          add = "reg.line",
          conf.int = TRUE, color = "sex", size = 3, alpha = 0.5, 
          palette = c("#e6550d", "#43a2ca"), 
          cor.coef = T, cor.method = "pearson",
          title = "OAS1 change with age and sex in samples of non-infectious state",
          xlab = "Age", ylab = "Log OAS1 Level", repel = TRUE)

ggscatter(BQC, x = "ProcessTime.hrs.", y = "LOG_OAS", 
          add = "reg.line",
          conf.int = TRUE, 
          color = "GROUP", palette = c("#d7301f", "#78c679"), size = 3, alpha = 0.5, 
          cor.coef = T, cor.method = "pearson",
          title = "OAS1 v Processing time",
          xlab = "Processing Time", ylab = "Log OAS1 Level", repel = TRUE)

ggscatter(BQC, x = "citrated.PC1", y = "LOG_OAS", 
          add = "reg.line",
          conf.int = TRUE, 
          color = "GROUP", size = 3, alpha = 0.6, 
          cor.coef = T, cor.method = "pearson",
          title = "OAS1 v indicator protein",
          xlab = "Indicator protein PC1", ylab = "Log OAS1 level", repel = TRUE)



A$GROUP <- c("During Infection")
B$GROUP <- c("Non-infectious State")
BQC <- rbind(A, B)

install.packages("ggridges")

library(ggridges)
library(gridExtra)
library(ggplot2)
library(viridis)
library(hrbrthemes)

ggplot(BQC, aes(x = `LOG_OAS`, y = `GROUP`)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Days_symptom_update3", option = "E") +
  #labs(title = 'OAS1 Level Distribution in Two Groups Before QC') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )


ggplot(BQC, aes(x=LOG_OAS, y=GROUP, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 10, quantile_lines = TRUE
  ) +
  theme_ipsum() +
  scale_fill_viridis_d(name = "Quantiles", option = "magma")


fit = glm(B2 ~  LOG_OAS + age2 + sex + age_at_diagnosis , family=binomial(), data = A)
summary(fit)



#Non_INFE_30 <- Non_INFE_30 %>% group_by(anonymized_patient_id) %>% mutate(OAS_mean = mean(OAS1.10361.25, na.rm = TRUE))
###take last measurement for 31


#Non_INFE_31_rm <- dplyr::add_count(Non_INFE_31_rm, anonymized_patient_id)
Non_INFE_31_rm <- Non_INFE_31_rm[!is.na(Non_INFE_31_rm$OAS_mean),]
Non_INFE_31_rm$LOG_OAS_mean <- log(Non_INFE_31_rm$OAS_mean)
Non_INFE_31_rm_2 <- unique(Non_INFE_31_rm[, c("age_at_diagnosis", "sex", "SampleGroup", "LOG_OAS_mean", "A1", "A2", "B1", "B2", "C1", "ProcessTime.hrs.")])

Non_INFE_31_rm_2$LOG_OAS_mean_res_scale <- scale(residuals(glm(LOG_OAS_mean ~ ProcessTime.hrs., data=Non_INFE_31_rm_2)))
fit = glm(A2 ~  LOG_OAS_mean_res_scale + age2 + sex + age_at_diagnosis + sex*age_at_diagnosis, family=binomial(), data = Non_INFE_31_rm_2)
summary(fit)



Non_INFE <- dplyr::add_count(Non_INFE, anonymized_patient_id)
INFE <- dplyr::add_count(INFE, anonymized_patient_id)
Non_INFE_30 <- dplyr::add_count(Non_INFE_30, anonymized_patient_id)
Non_INFE_31 <- dplyr::add_count(Non_INFE_31, anonymized_patient_id)

Non_INFE <- Non_INFE[!is.na(Non_INFE$OAS_mean),]
Non_INFE$LOG_OAS_mean <- log(Non_INFE$OAS_mean)

INFE <- INFE[!is.na(INFE$OAS_mean),]
INFE$LOG_OAS_mean <- log(INFE$OAS_mean)

Non_INFE_30 <- Non_INFE_30[!is.na(Non_INFE_30$OAS_mean),]
Non_INFE_30$LOG_OAS_mean <- log(Non_INFE_30$OAS_mean)

Non_INFE_31 <- Non_INFE_31[!is.na(Non_INFE_31$OAS_mean),]
Non_INFE_31$LOG_OAS_mean <- log(Non_INFE_31$OAS_mean)

write.table(Non_INFE, file="All_OAS1_merge_medNormRefSMP_sampleINFO.Non_INFE.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(INFE, file="All_OAS1_merge_medNormRefSMP_sampleINFO.INFE.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(Non_INFE_30, file="All_OAS1_merge_medNormRefSMP_sampleINFO.Non_INFE_30.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(Non_INFE_31, file="All_OAS1_merge_medNormRefSMP_sampleINFO.Non_INFE_31.txt",col.names=T,row.names=F,quote=F,sep="\t")

INFE_2 <- unique(INFE[, c("age_at_diagnosis", "sex", "SampleGroup", "LOG_OAS_mean", "A1", "A2", "B1", "B2", "C1", "ProcessTime.hrs.")])

Non_INFE_30_2 <- unique(Non_INFE_30[, c("age_at_diagnosis", "sex", "SampleGroup", "LOG_OAS_mean", "A1", "A2", "B1", "B2", "C1")])
Non_INFE_31_2 <- unique(Non_INFE_31[, c("age_at_diagnosis", "sex", "SampleGroup", "LOG_OAS_mean", "A1", "A2", "B1", "B2", "C1", "ProcessTime.hrs.")])
Non_INFE_31_2$age2 <- Non_INFE_31_2$age_at_diagnosis^2
Non_INFE_31_2$LOG_OAS_mean_res_scale <- scale(residuals(glm(LOG_OAS_mean ~ ProcessTime.hrs., data=Non_INFE_31_2)))
fit = glm(A2 ~  LOG_OAS_mean_res_scale + age2 + sex + age_at_diagnosis + sex*age_at_diagnosis, family=binomial(), data = Non_INFE_31_2)
summary(fit)


nrow(INFE_2[which(INFE_2$A1 == "1"), ])

nrow(INFE[which(Non_INFE_35_2$A1 == "1"), ])

fit = glm(A1 ~  LOG_OAS_mean + age_at_diagnosis + sex , family=binomial(), data = Non_INFE_31_2)
summary(fit)


fit = glm(C1 ~  LOG_OAS_mean + age_at_diagnosis + sex + SampleGroup, family=binomial(), data = INFE_2)
summary(fit)



library(ggplot2)
#test$Protein = factor(test$Protein, levels=c("OAS1","IL10RB","ABO","LBP"))

test$Group = factor(test$Group, levels=c("Severe", "Hospitalized", "Susceptible"))
test=data.frame(read.table("OAS1.res.txt",header=T,sep = "\t"))
p = ggplot(data=test,
           aes(x = COVID, y = OR, ymin = CIL, ymax = CIU ), width=0.1, cex=1.2)+
  geom_pointrange(aes(col=COVID))+
  geom_hline(aes(fill=COVID), yintercept =1, linetype=2)+
  xlab('Groups based on active infection state')+ ylab("Odds Ratio (95% Confidence Interval) per 1SD increase of log OAS1 level") +
  labs(color='COVID-19 Outcomes') + 
  geom_errorbar(aes(ymin=CIL, ymax=CIU,col=COVID),width=0.1,cex=1.2)+ 
  facet_wrap(~Group,strip.position="left",nrow=4,scales = "free_y") +
  theme_classic() + 
  theme(plot.title=element_text(size=14,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"), strip.text.y = element_text(size=10, face="bold")) +
  scale_y_continuous(breaks=seq(0,10,0.5))+
  coord_flip()
p


########

library(ggplot2)
library(grid)
library(gridExtra)
require(scales)
test=data.frame(read.table("OAS1_MR2.txt",header=T,sep = "\t"))

test$outcome <- factor(test$outcome, levels=c("Very Severe", "Hospitalization", "Susceptibility"))
test$exposure <- factor(test$exposure, levels=c("OAS1 pQTL", "OAS1 eQTL", "ABO pQTL", "OAS1 sQTL", "IL10RB pQTL"))

p <- ggplot(data=test, 
	aes(x = outcome, y = OR, ymin = LL, ymax = UL), width=0.1, cex=1.2) +
  geom_pointrange(aes(col=outcome))+
  geom_hline(aes(fill=outcome), yintercept =1, linetype=2)+
  xlab('Exposures')+ ylab("Odds Ratio (95% Confidence Interval) per unit* increase of log exposure level") +
  labs(color='COVID-19 Outcomes') + 
  geom_errorbar(aes(ymin=LL, ymax=UL, col=outcome), width=0.1,cex=1.2)+ 
  facet_wrap(~exposure,strip.position="left",nrow=4,scales = "free_y") +
  theme_classic() + 
  theme(plot.title=element_text(size=14,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"), strip.text.y = element_text(size=10, face="bold")) +
  scale_y_continuous(breaks=seq(0,1.3,0.1)) +
  coord_flip()
  
 p + scale_color_manual(values=c('#f16913','#9e9ac8','#6baed6')) + theme(legend.position = c(0.8, 0.1)) + guides(colour = guide_legend(reverse=TRUE)) 
  
#names(test) <- c("labels", "eventnum", "arr")
test$Title <- factor(test$Title, rev(levels(test$Title)))
test$colour <- rep(c("white"), 9)


data_table <- ggplot(data = test, aes(y = Title)) +
  geom_hline(aes(yintercept = Title, colour = colour), size = 3) +
  geom_text(aes(x = 0, label = Title), hjust = 0) +
  geom_text(aes(x = 10, label = OR)) +
  geom_text(aes(x = 15, label = pval), hjust = 1) +
  scale_colour_identity() +
  theme_void() + 
  theme(plot.margin = margin(5, 0, 35, 0))

grid.arrange(data_table,p, ncol = 2)


###

test=data.frame(read.table("OAS1.res.txt",header=T,sep = "\t"))
test$outcome <- factor(test$outcome, levels=c("Very Severe", "Hospitalization", "Susceptibility"))
p = ggplot(data=test,
           aes(x = outcome, y = OR, ymin = CIL, ymax = CIU ), width=0.1, cex=1.2)+
  geom_pointrange(aes(col=outcome))+
  geom_hline(aes(fill=outcome), yintercept =1, linetype=2)+
  xlab('')+ ylab("Odds Ratio (95% Confidence Interval) per 1SD increase of log OAS1 level") +
  labs(color='COVID-19 Outcomes') + 
  geom_errorbar(aes(ymin=CIL, ymax=CIU,col=outcome),width=0.1,cex=1.2)+ 
  facet_wrap(~Group,strip.position="left",nrow=4,scales = "free_y") +
  theme_classic() + 
  theme(plot.title=element_text(size=14,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"), strip.text.y = element_text(size=10, face="bold")) +
  scale_y_continuous(breaks=seq(0,7,0.5))+
  #scale_y_continuous(trans = log10_trans(),
   #                  breaks = trans_breaks("log10", function(x) 10^x, n=8),
    #                 labels = trans_format("log10", math_format(10^.x))) +
  coord_flip()
 p + scale_color_manual(values=c('#f16913','#9e9ac8','#6baed6')) + theme(legend.position = c(0.8, 0.3)) + guides(colour = guide_legend(reverse=TRUE))

###some plots###


library(ggpubr)

ggscatter(OAS1_dat, x = "citrated.PC1", y = "LOG_OAS1", 
         # add = "reg.line",
          conf.int = TRUE, color = "SampleGroup", size = 3, alpha = 0.6, 
          cor.coef = T, cor.method = "pearson",
          title = "Indicator protein v. log OAS1 level",
          xlab = "Indicator protein PC1", ylab = "log OAS1 level", repel = TRUE)


ggboxplot(D0_30_merge_B2, x = "Draw", y = "invnorm_OAS1",
          color = "B2", fill = "B2", palette = "jco",
          alpha = 0.5, ggtheme = theme_bw(), xlab = "Sampling date", title = "OAS1 level to B2 in D0 and D30")
par(mar=c(5,11,2,3))
par(cex.axis=0.8)
boxplot(D0[ which(D0$B1 == "1"), ]$LOG_OAS1, D0[ which(D0$B1 == "0"), ]$LOG_OAS1, 
        D30[ which(D30$B1 == "1"), ]$LOG_OAS1, D30[ which(D30$B1 == "0"), ]$LOG_OAS1,
        D30_nov[ which(D30_nov$B1 == "1"), ]$LOG_OAS1, D30_nov[ which(D30_nov$B1 == "0"), ]$LOG_OAS1,
        at = c(1,2,4,5,7,8),
        names = c("Case at D0", "Control at D0", "Case at D30/D0-", "Control at D30/D0-", 
                  "Case at D30/D0- no infection", "Control at D30/D0- no infection"),
        las = 2,
        xlab="Log OAS1 level",
        main="Hospitalized covid vs. not hospitalized covid (B1)",
        col = c("#3182bd","#31a354"),
        border = "#756bb1",
        horizontal = TRUE,
        notch = TRUE
)

#C1D0 <- OAS1_dat[ which(OAS1_dat$C1 == "1" & OAS1_dat$Draw == "D0"), ]


####check indicator proteins

fit = glm(C1 ~  EDTA.PC1, data = D0)
summary(fit)
fit = glm(C1 ~  EDTA.PC1, data = D30)
summary(fit)
fit = glm(C1 ~  citrated.PC1, data = D0)
summary(fit)
fit = glm(C1 ~  citrated.PC1, data = D30)
summary(fit)



fit = glm(OAS1.10361.25 ~  EDTA.PC1, data = D0)
summary(fit)
fit = glm(OAS1.10361.25 ~  citrated.PC1, data = D0)
summary(fit)
fit = glm(OAS1.10361.25 ~  EDTA.PC1, data = D30)
summary(fit)
fit = glm(OAS1.10361.25 ~  citrated.PC1, data = D30)
summary(fit)

library(statVisual)
p=statVisual(type = "LinePlot",
           data = OAS1_dat,
           x = 'Days_symptom',
           y = 'LOG_OAS1',
           sid = 'anonymized_patient_id',
           group = 'B2')
p + scale_x_discrete(breaks = seq(0, 140, by = 10))

#D30_ID <- OAS1_dat[ which(OAS1_dat$Draw == "D30"), 1:2]
#D30_ID_OAS1 <- merge(OAS1_dat, D30_ID, by = "anonymized_patient_id")
#D30_ID_OAS1$B2 <- gsub('1', 'Case', D30_ID_OAS1$B2)
#D30_ID_OAS1$B2 <- gsub('0', 'Control', D30_ID_OAS1$B2)
#D30_ID_OAS1$A2 <- gsub('1', 'Case', D30_ID_OAS1$A2)
#D30_ID_OAS1$A2 <- gsub('0', 'Control', D30_ID_OAS1$A2)
#write.table(D30_ID_OAS1, file="OAS1_4_timepoints.txt",col.names=T,row.names=F,quote=F,sep="\t")

OAS1_dat_count_sort <- dplyr::add_count(OAS1_dat, anonymized_patient_id)

OAS1_dat_count_sort <- OAS1_dat_count_sort[order(OAS1_dat_count_sort$Days_symptom), ]
OAS1_dat_count_sort$B2 <- gsub('1', 'Case', OAS1_dat_count_sort$B2)
OAS1_dat_count_sort$B2 <- gsub('0', 'Control', OAS1_dat_count_sort$B2)
OAS1_dat_count_sort$A2 <- gsub('1', 'Case', OAS1_dat_count_sort$A2)
OAS1_dat_count_sort$A2 <- gsub('0', 'Control', OAS1_dat_count_sort$A2)
OAS1_dat_count_sort_4timepoints <- OAS1_dat_count_sort[ which(OAS1_dat_count_sort$n == "4"), ]




D0 <- read.table("D0.pheno1.txt",header=T)
D0_merge <- merge(OAS1.dat.all , D0, by="SubjectID")
D30 <- read.table("D30.pheno1.txt",header=T)
D30_merge <- merge(OAS1.dat.all , D30, by="SubjectID")
CTRL <- read.table("D30_CTRL.txt",header=T)
CTRL_merge <- merge(OAS1.dat.all , CTRL, by="SubjectID")
D0_30 <- read.table("D0_D30.pheno1.txt",header=T)
D0_30_merge <- merge(OAS1.dat.all , D0_30, by="SubjectID")

D0_30_merge<-D0_30_merge[!is.na(D0_30_merge$OAS1.10361.25),]
OAS1<-log(D0_30_merge$OAS1.10361.25)
Model<-lm(OAS1 ~ PlateId + SampleGroup, data=D0_30_merge)
Residuals<-residuals(Model)
InvNorm_var<-qnorm((rank(Residuals,na.last="keep")-0.5)/sum(!is.na(Residuals)))
D0_30_merge$invnorm_OAS1<-InvNorm_var

D0_merge<-D0_merge[!is.na(D0_merge$OAS1.10361.25),]
OAS1<-log(D0_merge$OAS1.10361.25)
Model<-lm(OAS1 ~ PlateId + SampleGroup, data=D0_merge)
Residuals<-residuals(Model)
InvNorm_var<-qnorm((rank(Residuals,na.last="keep")-0.5)/sum(!is.na(Residuals)))
D0_merge$invnorm_OAS1<-InvNorm_var

CTRL_merge<-CTRL_merge[!is.na(CTRL_merge$OAS1.10361.25),]
OAS1<-log(CTRL_merge$OAS1.10361.25)
Model<-lm(OAS1 ~ PlateId + SampleGroup, data=CTRL_merge)
Residuals<-residuals(Model)
InvNorm_var<-qnorm((rank(Residuals,na.last="keep")-0.5)/sum(!is.na(Residuals)))
CTRL_merge$invnorm_OAS1<-InvNorm_var

hist(CTRL_merge$OAS1.10361.25)
hist(CTRL_merge$invnorm_OAS1)

fit = lm(invnorm_OAS1 ~  Age + Sex, data = D0_merge)
summary(fit)
fit = lm(invnorm_OAS1 ~  Age + Sex, data = CTRL_merge)
summary(fit)



D0_30_merge_B1<-D0_30_merge[!is.na(D0_30_merge$B1),]
#D0_30_merge_B1$B1 <- gsub('0', 'Control', D0_30_merge_B1$B1)
#D0_30_merge_B1$B1 <- gsub('1', 'Case', D0_30_merge_B1$B1)

D0_merge_B1 <- D0_30_merge_B1[ which(D0_30_merge_B1$Draw == "D0"), ]
D30_merge_B1 <- D0_30_merge_B1[ which(D0_30_merge_B1$Draw == "D30"), ]


D0_30_merge_B2<-D0_30_merge[!is.na(D0_30_merge$B2),]
#D0_30_merge_B2$B2 <- gsub('0', 'Control', D0_30_merge_B2$B2)
#D0_30_merge_B2$B2 <- gsub('1', 'Case', D0_30_merge_B2$B2)

D0_30_merge_B2$B2 <- gsub('Control', '0', D0_30_merge_B2$B2)
D0_30_merge_B2$B2 <- gsub('Case', '1', D0_30_merge_B2$B2)

D0_merge_B2 <- D0_30_merge_B2[ which(D0_30_merge_B2$Draw == "D0"), ]
D30_merge_B2 <- D0_30_merge_B2[ which(D0_30_merge_B2$Draw == "D30"), ]

fit = glm(B2 ~  invnorm_OAS1 + Age + Sex, data = D30_merge_B2)
summary(fit)
fit = glm(B2 ~  invnorm_OAS1 + Age + Sex, data = D0_merge_B2)
summary(fit)

fit = lm(B1 ~  invnorm_OAS1 + Age + Sex, data = D30_merge_B1)
summary(fit)

D0_merge_B1<-D0_merge[!is.na(D0_merge$B1),]

ggboxplot(D0_30_merge_B2, x = "Draw", y = "invnorm_OAS1",
          color = "B2", fill = "B2", palette = "jco",
          alpha = 0.5, ggtheme = theme_bw(), xlab = "Sampling date", title = "OAS1 level to B2 in D0 and D30")

ggboxplot(D0_merge, x = "Sex", y = "invnorm_OAS1",
          color = "CaseUpdate", fill = "CaseUpdate", palette = "jco",
          alpha = 0.5, ggtheme = theme_bw(), xlab = "Sex")

#B2D0_VDR <- merge(VDR.dat.JGH , D0, by="SubjectID")
library(ggpubr)
p=ggplot(B2D0, aes(x = OAS1.10361.25, fill=B2)) + geom_density(alpha = 0.7) + xlab("OAS1 level") + labs(fill = "COVID19")
p+scale_fill_manual(values=c("#fbb4ae", "#b3cde3"))
dim(B2D0)

Case_OAS1 <- B2D0[ which(B2D0$B2 == 1 & B2D0$OAS1.10361.25<1500), 2]
Control_OAS1 <- B2D0[ which(B2D0$B2 == 0 & B2D0$OAS1.10361.25<1500), 2]

Case_VDR <- B2D0_VDR[ which(B2D0$B2 == 1), 2]
Control_VDR <- B2D0_VDR[ which(B2D0$B2 == 0), 2]

Case_OAS1_D0 <- B2D0[ which(B2D0$B2 == 1), 2]
Control_OAS1_D0 <- B2D0[ which(B2D0$B2 == 0), 2]
Case_OAS1_D30 <- B2D30[ which(B2D30$B2 == 1), 2]
Control_OAS1_D30 <- B2D30[ which(B2D30$B2 == 0), 2]
par(mar=c(4,9,2,3))
boxplot(log(Case_OAS1_D0), log(Control_OAS1_D0), log(Case_OAS1_D30), log(Control_OAS1_D30),
        at = c(1,2,4,5),
        names = c("CASE_B2D0", "CONTROL_B2D0", "CASE_B2D30", "CONTROL_B2D30"),
        las = 2,
        xlab="Log raw OAS1 level",
        col = c("orange","red"),
        border = "brown",
        horizontal = TRUE,
        notch = TRUE
)

B2_a_s <- read.table("B2.age.sex.txt",header=T)
B2D0_VDR <- merge(B2D0_VDR , B2_a_s, by="SubjectID")

fit = lm(log(B2D0_VDR$VDR.10023.32)~B2D0$B2 + B2D0$Age + B2D0$Sex)
summary(fit)

B2D30 <- read.table("B2D30.txt",header=T)
B2D30 <- merge(OAS1.dat.JGH , B2D30, by="SubjectID")

boxplot(log(B2D30$OAS1.10361.25)~B2D30$B2,
        at = c(1,2),
        #names = c("CASE_B2", "CONTROL_B2"),
        las = 2,
        ylim = c(5.5,8),
        col = c("orange","red"),
        border = "brown",
        horizontal = TRUE,
        notch = TRUE
)

fit = lm(log(B2D30$OAS1.10361.25)~B2D30$B2 + B2D30$Age + B2D30$Sex)
summary(fit)




############add genotype#####################

OAS1 <- read.table("All_OAS1_merge_medNormRefSMP_sampleINFO_update_505.txt",header=T)
Gen <- read.table("SOMA_OAS1_Genotype.txt",header=T)
OAS1 <- merge(OAS1, Gen, by = "anonymized_patient_id")

OAS1 <- OAS1[order(OAS1$Days_symptom_update3),] 

baseline <- OAS1 %>% 
  group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == first(OAS1.10361.25))

baseline_QC <- baseline[ which(baseline$LOG_OAS1 < 8), ]


OAS1_no_norm <- dat %>% select(SubjectID, OAS1.10361.25)

baseline <- merge(baseline, OAS1_no_norm, by ="SubjectID")


baseline$LOG_OAS1_no_norm <- log(baseline$OAS1.10361.25.y) 




####violin+box####

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(hrbrthemes)
library(viridis)




p <- ggplot(aes(y = LOG_OAS1,
                x = sQTL),
                #fill = Ancestry),
            data = baseline_QC) +
  #geom_violin(position=position_dodge(),draw_quantiles=c(0.5)) +
  geom_violin(width=0.7, position="dodge") +
  geom_boxplot(width=0.4,color="black", alpha=0.5, position = position_dodge(width =0.7))+
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..)),
               position=position_nudge(x=0.33), size=3.5) +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum() +
  labs(title="",
       x="", 
       y="0AS1 Level") +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size=11)
  ) +
 # facet_grid(`CaseUpdate`~`Ancestry`) + 
  scale_fill_brewer(palette="Accent")

p


res.aov <- aov(LOG_OAS1 ~ sQTL, data = baseline_QC)
summary(res.aov)


viral <- read.table("Viral.txt",header=T)

fit = lm(viral$oas1 ~ viral$E + viral$Age + viral$Sex)
summary(fit)


###select 100 JGH samples for ELISA####
OAS1 <- read.table("All_OAS1_merge_medNormRefSMP_sampleINFO_update_505.txt",header=T)
ELISA <- read.table("ELISA_sample.txt",header=T)
OAS1_ELISA <- merge(OAS1, ELISA, by="SubjectID")



OAS1_ELISA <- OAS1_ELISA[order(OAS1_ELISA$LOG_OAS1),]


library(Hmisc)

OAS1_ELISA$Bins <- as.numeric(cut_number(OAS1_ELISA$LOG_OAS1,10))

library(dplyr)
OAS1_ELISA_100 <- OAS1_ELISA %>% group_by(Bins) %>% sample_n(10)

write.table(OAS1_ELISA_100, file="OAS1_ELISA_100.txt",col.names=T,row.names=F,quote=F,sep="\t")
