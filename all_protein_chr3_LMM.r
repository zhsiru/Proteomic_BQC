setwd("C:\\Users\\sirui.zhou\\work\\Works-home\\Covid\\SomaLogic")
setwd("C:\\Users\\Sirui\\Desktop\\WORKS\\OAS1.NM.revision\\McGill-Richards-C-19 - SomaScan Data")

library(SomaDataIO)
Norm <- read.adat("SS-200150_v4_ACDPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.medNormRefSMP.adat")
#Unnorm <- read.adat("SS-200150_v4_ACDPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.adat")
library(dplyr)
Norm2 <- Norm %>% select(-(1:22),-(24:33))

Non_InF <- read.table("Non_INFE_31_last_505_short.txt",header=T)
InF <- read.table("INFE_1st_505_short.txt",header=T)
All <- read.table("All_OAS1_merge_medNormRefSMP_sampleINFO_update_505_short.txt",header=T)
All_QC <- All[ which(All$anonymized_patient_id != "JGH_366"), ]
InF_QC <- InF[ which(InF$anonymized_patient_id != "JGH_366"), ]
Non_InF_QC <- Non_InF[ which(Non_InF$anonymized_patient_id != "JGH_366"), ]


chr3 <- read.table("chr3_BQC.txt",header=T)
EUR <- read.table("EUR.sample.PCs", header=T)

chr3_eur_pc <- merge(chr3, EUR, by="IID")

Non_InF_chr3_eur <- merge(Non_InF, chr3_eur_pc, by="anonymized_patient_id")

InF_chr3_eur <- merge(InF, chr3_eur_pc, by="anonymized_patient_id")

Non_InF_chr3_eur_QC <- Non_InF_chr3_eur[ which(Non_InF_chr3_eur$NormScale_20 < 4 & Non_InF_chr3_eur$NormScale_20 > 0.25), ]
Non_InF_chr3_eur_QC <- Non_InF_chr3_eur_QC[ which(Non_InF_chr3_eur_QC$NormScale_0_005 < 4 & Non_InF_chr3_eur_QC$NormScale_0_005 > 0.25), ]
Non_InF_chr3_eur_QC <- Non_InF_chr3_eur_QC[ which(Non_InF_chr3_eur_QC$NormScale_0_5 < 4 & Non_InF_chr3_eur_QC$NormScale_0_5 > 0.25), ]


InF_chr3_eur_QC <- InF_chr3_eur[ which(InF_chr3_eur$NormScale_20 < 4 & InF_chr3_eur$NormScale_20 > 0.25), ]
InF_chr3_eur_QC <- InF_chr3_eur_QC[ which(InF_chr3_eur_QC$NormScale_0_005 < 4 & InF_chr3_eur_QC$NormScale_0_005 > 0.25), ]
InF_chr3_eur_QC <- InF_chr3_eur_QC[ which(InF_chr3_eur_QC$NormScale_0_5 < 4 & InF_chr3_eur_QC$NormScale_0_5 > 0.25), ]

Non_InF_chr3_eur_QC_all_p <- merge(Non_InF_chr3_eur_QC, Norm2, by="SubjectID")
InF_chr3_eur_QC_all_p <- merge(InF_chr3_eur_QC, Norm2, by="SubjectID")


#####define outcome exposure####

# outcome
out_start= 47
out_end= 5330
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

# exposure
exp_start=34
exp_end=34
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)

number=1

inp=Non_InF_chr3_eur_QC_all_p

####LMM######

library(lme4)
library(tidyverse)
for (i in out_start:out_end){
  outcome = colnames(inp)[i]
  outcome2 = scale(residuals(glm(log(get(outcome)) ~ ProcessTime , data = inp)))
  for (j in exp_start:exp_end){
    exposure = colnames(inp)[j]
    model <- lmer(outcome2 ~ get(exposure) + age_at_diagnosis + sex + (1|PlateId) + (1|SampleGroup) + PC1 + PC2 + PC3 + PC4 + PC5,
                  na.action = na.exclude,
                  data=inp)
    Vcov <- vcov(model, useScale = FALSE)
    beta <- fixef(model)
    se <- sqrt(diag(Vcov))
    zval <- beta / se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    
    out_beta[number] = as.numeric(beta[2])
    out_se[number] = as.numeric(se[2])
    out_pvalue[number] = as.numeric(pval[2])
    out_variable[number] = outcome
    number = number + 1
    
    exp_beta[number] = as.numeric(beta[2])
    exp_se[number] = as.numeric(se[2])
    exp_pvalue[number] = as.numeric(pval[2])
    exp_variable[number] = exposure
    number = number + 1
  }
}


###GLM####

for (i in out_start:out_end){
  outcome = colnames(inp)[i]
  outcome2 = scale(residuals(glm(log(get(outcome)) ~ ProcessTime , data = inp)))
  for (j in exp_start:exp_end){
    exposure = colnames(inp)[j]
    model <- glm(outcome2 ~ get(exposure) + age_at_diagnosis + sex + PlateId + SampleGroup + PC1 + PC2 + PC3 + PC4 + PC5,
                  data=inp)
    beta <- coef(summary(model))[,1]
    se <- coef(summary(model))[,2]
    pval <- coef(summary(model))[,4]
    
    out_beta[number] = as.numeric(beta[2])
    out_se[number] = as.numeric(se[2])
    out_pvalue[number] = as.numeric(pval[2])
    out_variable[number] = outcome
    number = number + 1
    
    exp_beta[number] = as.numeric(beta[2])
    exp_se[number] = as.numeric(se[2])
    exp_pvalue[number] = as.numeric(pval[2])
    exp_variable[number] = exposure
    number = number + 1
  }
}

#####Output transformation####


outcome = na.omit(data.frame(out_variable, out_beta, out_se, out_pvalue))

exposure = na.omit(data.frame(exp_variable, exp_beta, exp_se, exp_pvalue))




outcome = outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    pvalue = out_pvalue,
  )

outcome$variable <- sub("^", "Somamer_", outcome$variable )

exposure = exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    pvalue = exp_pvalue,
  )
all = rbind(outcome, exposure)


data = all %>% 
  mutate(
    type = substr(variable, 1, 2)
  ) %>% 
  spread(type, variable) %>% 
rename(
   Somamer  = So,
   rs10490770 = a1
  ) %>% 
  mutate (
    beta = round(beta, 5),
    se = round(se, 5),
    pvalue = round(pvalue, 5)
  ) %>% 
select(Somamer, rs10490770, beta, se, pvalue)

write.table(data, file="NONINFE_QC_EUR_All_protein_rs10490770_LMM_txt",col.names=T,row.names=F,quote=F,sep="\t")


####test###
FLT1_stand = scale(residuals(glm(log(FLT1.8231.122) ~ ProcessTime , data = InF_chr3_eur_QC_all_p)))
model <- lmer(FLT1_stand ~ a1 + age_at_diagnosis + sex + (1|PlateId) + (1|SampleGroup) + PC1 + PC2 + PC3 + PC4 + PC5,
              na.action = na.exclude,
              data=InF_chr3_eur_QC_all_p)
summary(model)


####make new A2####

All_QC$A2_new <- ifelse(All_QC$resp_severe == "1" | All_QC$resp_mild == "1" | All_QC$death == "1" , '1', '0')
nrow(All_QC[ which(All_QC$A2_new == "1"), ])

InF_QC$A3 <- ifelse(InF_QC$resp_severe == "1" | InF_QC$resp_mild == "1" | InF_QC$death == "1" , '1', '0')
nrow(InF_QC[ which(InF_QC$A2_new == "1"), ])  

Non_InF_QC$A3 <- ifelse(Non_InF_QC$resp_severe == "1" | Non_InF_QC$resp_mild == "1" | Non_InF_QC$death == "1" , '1', '0')


write.table(Non_InF_QC, file="Somalogic_NonInf_220_A3.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(InF_QC, file="Somalogic_Inf_418_A3.txt",col.names=T,row.names=F,quote=F,sep="\t")
