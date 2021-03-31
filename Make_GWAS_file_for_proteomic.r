setwd("C:\\Users\\sirui.zhou\\work\\BQC_protein_GWAS")
setwd("C:\\Users\\Sirui\\Desktop\\WORKS\\BQC_protein_GWAS")
library(dplyr)  
library(data.table)
library(ggplot2)
library(reshape2)
library(SomaDataIO)

my_adat <- read.adat("SS-200150_v4_ACDPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.medNormRefSMP.adat")
measurement <- my_adat[c(23, 34:5317)]


Soma <- read.table("Somalogic_503_new_A2.txt",header=T)
Soma <- Soma[order(Soma$Days_symptom_update3),] 
Baseline <- Soma %>% 
  group_by(anonymized_patient_id) %>% 
  filter(SubjectID == first(SubjectID))

write.table(Baseline, file="Baseline_503.txt",col.names=T,row.names=F,quote=F,sep="\t")

Soma_inf <- read.table("Somalogic_Inf_417_new_A2.txt",header=T)

Soma_inf_case <- merge(Soma_inf[which(Soma_inf$CaseUpdate == "Positive"), ], measurement, by = "SubjectID")
Soma_inf_case_log <- Soma_inf_case
Soma_inf_case_log[c(35:5318)] <- log(Soma_inf_case[c(35:5318)])


EUR_Soma_PC <- fread("Soma_inf_covidpos_EUR_PC.eigenvec", header=F, select = c(2:7))
colnames(EUR_Soma_PC) <- c("GWAS_ID","PC1","PC2","PC3","PC4","PC5")


Soma_inf_case_EUR <- merge(Soma_inf_case_log, EUR_Soma_PC, by = "GWAS_ID")
Soma_inf_case_EUR_pheno <- Soma_inf_case_EUR[c(1,35:5318)]
write.csv(Soma_inf_case_EUR_pheno, file="Soma_inf_case_EUR_pheno.csv",row.names=F)

Soma_inf_case_EUR_qcov <- Soma_inf_case_EUR[c(1,1,11,5319:5323)]
Soma_inf_case_EUR_cov <- Soma_inf_case_EUR[c(1,1,12)]
write.table(Soma_inf_case_EUR_qcov, file="Soma_inf_case_EUR_qcov.txt",col.names=F,row.names=F,quote=F,sep=" ")
write.table(Soma_inf_case_EUR_cov, file="Soma_inf_case_EUR_cov.txt",col.names=F,row.names=F,quote=F,sep=" ")



Lum <- read.table("Luminex_adj_rm_320_QC1.txt",header=T)

Lum <- Lum[order(Lum$Date_sampling), ]
Lum$DaysSinceSymptomOnset <- as.numeric(Lum$DaysSinceSymptomOnset)
Lum_1st <- Lum %>% 
  group_by(GWAS_ID) %>% 
  filter(CCL2 == first(CCL2))

Lum_case <- Lum[ which(Lum$GWAS_ID != "NA" & Lum$covidstatus == "POSITIF"), ]

Lum_case_max <- Lum_case %>%
  group_by(GWAS_ID) %>%  
  summarise_at(vars(2:7,22:42), max)
Lum_case_max_norm <- Lum_case_max
Lum_case_max_norm[c(2:28)] <- scale(log(Lum_case_max[c(2:28)]))

Lum_case_mean <- Lum_case %>%
  group_by(GWAS_ID) %>%  
  summarise_at(vars(2:7,22:42), mean)
Lum_case_mean_norm <- Lum_case_mean
Lum_case_mean_norm[c(2:28)] <- scale(log(Lum_case_mean[c(2:28)]))

Lum_case_inf <- Lum_case[ which(Lum_case$DaysSinceSymptomOnset < 14), ]

Lum_case_inf_max <- Lum_case_inf %>%
  group_by(GWAS_ID) %>%  
  summarise_at(vars(2:7,22:42), max)
Lum_case_inf_max_norm <- Lum_case_inf_max
Lum_case_inf_max_norm[c(2:28)] <- scale(log(Lum_case_inf_max[c(2:28)]))

Lum_case_inf_mean <- Lum_case_inf %>%
  group_by(GWAS_ID) %>%  
  summarise_at(vars(2:7,22:42), mean)
Lum_case_inf_mean_norm <- Lum_case_inf_mean
Lum_case_inf_mean_norm[c(2:28)] <- scale(log(Lum_case_inf_mean[c(2:28)]))



ggplot(melt(Lum_case_mean_norm),aes(x=value)) + geom_histogram() + facet_wrap(~variable)

write.csv(Lum_case_max_norm, file="Lum_case_max_norm.csv",row.names=F)
write.csv(Lum_case_mean_norm, file="Lum_case_mean_norm.csv",row.names=F)
write.csv(Lum_case_inf_max_norm, file="Lum_case_inf_max_norm.csv",row.names=F)
write.csv(Lum_case_inf_mean_norm, file="Lum_case_inf_mean_norm.csv",row.names=F)



EUR_Lum_PC <- fread("Lum_case_EUR_PC.eigenvec", header=F, select = c(2:7))
colnames(EUR_Lum_PC) <- c("GWAS_ID","PC1","PC2","PC3","PC4","PC5")
EUR_Lum_PC$POP <- "EUR"
AFR_Lum_PC <- fread("Lum_case_AFR_PC.eigenvec", header=F, select = c(2:7))
colnames(AFR_Lum_PC) <- c("GWAS_ID","PC1","PC2","PC3","PC4","PC5")
AFR_Lum_PC$POP <- "AFR"
EAS_Lum_PC <- fread("Lum_case_EAS_PC.eigenvec", header=F, select = c(2:7))
colnames(EAS_Lum_PC) <- c("GWAS_ID","PC1","PC2","PC3","PC4","PC5")
EAS_Lum_PC$POP <- "EAS"
SAS_Lum_PC <- fread("Lum_case_SAS_PC.eigenvec", header=F, select = c(2:7))
colnames(SAS_Lum_PC) <- c("GWAS_ID","PC1","PC2","PC3","PC4","PC5")
SAS_Lum_PC$POP <- "SAS"
PC <- unique(rbind(EUR_Lum_PC, AFR_Lum_PC, EAS_Lum_PC, SAS_Lum_PC))

Lum_case_max_EUR <- merge(Lum_case_max_norm, EUR_Lum_PC, by = "GWAS_ID")[c(1,8:28)]
write.csv(Lum_case_max_EUR, file="Lum_case_max_EUR.csv",row.names=F)
Lum_case_mean_EUR <- merge(Lum_case_mean_norm, EUR_Lum_PC, by = "GWAS_ID")[c(1,8:28)]
write.csv(Lum_case_mean_EUR, file="Lum_case_mean_EUR.csv",row.names=F)

Lum_case_inf_max_EUR <- merge(Lum_case_inf_max_norm, EUR_Lum_PC, by = "GWAS_ID")[c(1,8:28)]
write.csv(Lum_case_inf_max_EUR, file="Lum_case_inf_max_EUR.csv",row.names=F)
Lum_case_inf_mean_EUR <- merge(Lum_case_inf_mean_norm, EUR_Lum_PC, by = "GWAS_ID")[c(1,8:28)]
write.csv(Lum_case_inf_mean_EUR, file="Lum_case_inf_mean_EUR.csv",row.names=F)



Lum_case_cov <- unique(Lum_case[c(12,10,17)])
Lum_case_cov_EUR <- merge(Lum_case_cov, EUR_Lum_PC, by = "GWAS_ID")[c(1,1,2,3)]
Lum_case_qcov <- unique(Lum_case[c(12,16)])
Lum_case_qcov_EUR <- merge(Lum_case_qcov, EUR_Lum_PC, by = "GWAS_ID")[c(1,1:7)]
write.table(Lum_case_cov_EUR, file="Lum_case_cov_EUR.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(Lum_case_qcov_EUR, file="Lum_case_qcov_EUR.txt",col.names=T,row.names=F,quote=F,sep="\t")


Lum_case_inf_cov <- unique(Lum_case_inf[c(12,10,17)])
Lum_case_inf_cov_EUR <- merge(Lum_case_inf_cov, EUR_Lum_PC, by = "GWAS_ID")[c(1,1,2,3)]
Lum_case_inf_qcov <- unique(Lum_case_inf[c(12,16)])
Lum_case_inf_qcov_EUR <- merge(Lum_case_inf_qcov, EUR_Lum_PC, by = "GWAS_ID")[c(1,1:7)]
write.table(Lum_case_inf_cov_EUR, file="Lum_case_inf_cov_EUR.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(Lum_case_inf_qcov_EUR, file="Lum_case_inf_qcov_EUR.txt",col.names=T,row.names=F,quote=F,sep="\t")


Lum_case_cov_AFR <- merge(Lum_case_cov, AFR_Lum_PC, by = "GWAS_ID")[c(1,2,3)]
Lum_case_qcov_AFR <- merge(Lum_case_qcov, AFR_Lum_PC, by = "GWAS_ID")[c(1:7)]


Lum_case_inf_cov_AFR <- merge(Lum_case_inf_cov, AFR_Lum_PC, by = "GWAS_ID")[c(1,2,3)]
Lum_case_inf_qcov_AFR <- merge(Lum_case_inf_qcov, AFR_Lum_PC, by = "GWAS_ID")[c(1:7)]

#Lum_1st_case_norm <- Lum_1st_case %>% mutate(across(23:43, scale))

write.csv(Lum_inf_1st_case_norm, file="Lum_inf_1st_case_norm.csv",row.names=F)
write.csv(Lum_1st_case_norm, file="Lum_1st_case_norm.csv",row.names=F)


###plots
library(CMplot) 
SPD <- read.table("SPD_inf_max_EUR_lr.fastGWA.0.05",header=T)

CMplot(SPD[,c(2,1,3,10)], plot.type="m",LOG10=TRUE,col=c("#bdc9e1","#045a8d"),threshold=5e-8,threshold.lty=2,file="jpg",memo="SPD",   
       amplify=FALSE,file.output=T,verbose=TRUE,width=14,height=6)

