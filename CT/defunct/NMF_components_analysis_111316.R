###################################
# some analyses for CT-DD project #

require(mgcv)
library(visreg)
library(ggplot2)
library(reshape2)
library(corrplot)

# add R functions 
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/regDataGam.R") # regional analysis
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/ROIanalysis.R")
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/ROIanalysis_flexible.R")
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/ROIanalysis_flexible2.R")
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/ROIanalysis_CtDep.R")
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/regDataGamMultOut.R")

# read in data
data<-read.table('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n452_pnc_itc_whole_sample_20161226.csv',header=T,sep=',')

# medications data 
meds<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/go1_go2_go3_medications_92816_ts_st_tsFinal.csv',header=T,sep=',')
names(meds)[names(meds)=="study"]<-"timepoint"
meds$medsComb<-meds$psychoactiveMedPsych_20160930+meds$psychoactiveMedMedical_20160930

# NMF names 
nmf20Names <-names(data)[1467:1486] 
nmf19Names <-names(data)[c(1467:1482,1484:1486)]

# subset data for models
data_k<-data[which(data$jlfCTInclude_k==1),]
# adding squared term for logk
data_k$ageSq<-(data_k$ageAtScan-mean(data_k$ageAtScan))^2
# logk quartiles
data_k$logkQ<-factor(cut(data_k$logk,breaks=quantile(data_k$logk,probs=seq(0,1, by=0.25),na.rm=TRUE),include.lowest=TRUE),labels=c("Q1","Q2","Q3","Q4"))
# ordered logk variable, need for GAM
data_k$logkQ.o <- ordered(data_k$logkQ)

#save(data_k,file='/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/DataK_n427_030717.Rda')

# merge 
comb<-merge(data_k, meds[,c("bblid","scanid","timepoint","medsComb")], by=c("bblid","timepoint"), all.x=T)
data_k_med <-comb[which(comb$medsComb==0),] # subsetting data with subjects NOT using medications

### saving data for matlab
#data_sub<-data_k[,c(1,171,203:207,1253,1255:1257,1467:1486,1659:1662,1266,1749,1752:1755,160,158)]
#write.csv(data_sub, file = "/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n427_forSVR_matlab_20170117.csv",row.names=F)

# descriptive statistics for k
mean(data_k$k)
sd(data_k$k)

#### demographic analyses for manuscript
# association with age
cor.test(data_k$logk,data_k$ageAtScan,method="pearson")
# association with sex
mean(data_k$logk[which(data_k$sex.o==1)])
sd(data_k$logk[which(data_k$sex.o==1)])

mean(data_k$logk[which(data_k$sex.o==2)])
sd(data_k$logk[which(data_k$sex.o==2)])
t.test(data_k$logk[which(data_k$sex.o==1)],data_k$logk[which(data_k$sex.o==2)])

wilcox.test(data_k$logk~data_k$sex.o)

# associations with age and sex
gmAgeSex1 <- gam(logk ~ s(ageAtScan) + sex.o, data=data_k)
gmAgeSex2 <- gam(logk ~ s(ageAtScan) + sex.o + s(ageAtScan, by=sex.o), data=data_k)

# create vector of model statements
ms<-paste0(nmf20Names[-17],"~logk+sex.o+s(ageAtScan)")
ms_cog<-paste0(nmf20Names[-17],"~F1_Exec_Comp_Res_Accuracy+sex.o+s(ageAtScan)")
ms_medu<-paste0(nmf20Names[-17],"~meduCnbGo1+sex.o+s(ageAtScan)")

ms1<-paste0("logk~",nmf20Names[-17],"+sex.o+s(ageAtScan)")
ms2<-paste0("logk~",nmf20Names[-17],"+sex.o+s(ageAtScan)")

# correlations between logk and CT components
cors_s<-cor(data_k$logk,data_k[,nmf20Names[-17]],method="spearman")
cors_p<-cor(data_k$logk,data_k[,nmf20Names[-17]],method="pearson")

# running basic models for association between NMF components and logk
modOut<-regDataGam(ms,data_k,nmf20Names[-17],1) # w/o the noise component
sigComp<-rownames(modOut[which(modOut$fdr.p<.05),])

# ECRA
modCog<-regDataGam(ms_cog,data_k,nmf20Names[-17],1)
sigComp_CR<-rownames(modCog[which(modCog$fdr.p<.05),])

modOut_ctDep<-regDataGam(ms,data,nmf20Names[-17],1) 

modOut_Cor<-cbind(modOut,t(cors_p),t(cors_s))
names(modOut_Cor)[8:9]<-c("pearsCor","spearCor")

modOut_ctDep_Cor<-cbind(modOut_ctDep,t(cors_p),t(cors_s))
names(modOut_ctDep_Cor)[8:9]<-c("pearsCor","spearCor")

# Medu
modMedu<-regDataGam(ms_medu,data_k,nmf20Names[-17],1)
sigComp_medu<-rownames(modMedu[which(modMedu$fdr.p<.05),])
modMedu[which(modMedu$fdr.p<.05),]

# saving results
saveRDS(modOut,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/logkCtNmf20Results_19comp.rds")
write.table(modOut_Cor,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/logkCtNmf20Results19compWithEffectSize.csv",row.names=T,col.names=T,sep=",")
write.table(modOut_ctDep_Cor,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/logkCtNmf20Results19compWithEffectSize_ctDep.csv",row.names=T,col.names=T,sep=",")

# MEDIATION MODELS, only for regions FDR-sig for logk

# cognitive reasoning
medResults<-ROIanalysis(sigComp,data)
saveRDS(medResults,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/mediationCrNmfCt20_just11sigComp.rds")

# BAS drive
medBASdr<-ROIanalysis_flexible(sigComp,data,"BAS_dr","sex.o + s(ageAtScan)")
write.table(medBASdr,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/mediationBASdrNmfCt20_just11sigComp.csv",row.names=F,col.names=T,sep=",")

# BAS reward responsiveness
medBASrr<-ROIanalysis_flexible(sigComp,data,"BAS_rr","sex.o + s(ageAtScan)")
write.table(medBASrr,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/mediationBASrrNmfCt20_just11sigComp.csv",row.names=F,col.names=T,sep=",")

# mother's education (no race)
medMeduNoRace<-ROIanalysis_flexible(sigComp,data,"meduCnbGo1","sex.o + s(ageAtScan)")
write.table(medMeduNoRace,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/mediationMedu1NoRaceNmfCt20_just11sigComp.csv",row.names=F,col.names=T,sep=",")

gm.test <- gam(Nmf20C14 ~ logk + meduCnbGo1 + F1_Exec_Comp_Res_Accuracy+  s(ageAtScan) + sex.o, data=data)
gm.test <- gam(Nmf20C15 ~ logk + meduCnbGo1 + F1_Exec_Comp_Res_Accuracy+  s(ageAtScan) + sex.o, data=data)
gm.test <- gam(Nmf20C20 ~ logk + meduCnbGo1 + F1_Exec_Comp_Res_Accuracy+  s(ageAtScan) + sex.o, data=data)

# mother's education (controlling for race)
medMedu <-ROIanalysis_flexible2(sigComp ,data, "meduCnbGo1")
write.table(medMedu,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/mediationMedu1NmfCt20_just11sigComp.csv",row.names=F,col.names=T,sep=",")

# environmental score (adjusting for race)
medEnv <-ROIanalysis_flexible2(sigComp ,data, "goassessEnvSes")
write.table(medEnv,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/mediationEnvNmfCt20_just11sigComp.csv",row.names=F,col.names=T,sep=",")

# environmental score (no race)
medEnvNoRace <-ROIanalysis_flexible(sigComp,data, "goassessEnvSes","sex.o + s(ageAtScan)")
write.table(medEnv,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/mediationEnvNmfCt20_just11sigComp.csv",row.names=F,col.names=T,sep=",")

# visualizations of CT-logk effects
gm.comp14 <- gam(Nmf20C14 ~ logk + s(ageAtScan) + sex.o, data=data_k)      
gm.comp15 <- gam(Nmf20C15 ~ logk + s(ageAtScan) + sex.o, data=data_k)
gm.comp18 <- gam(Nmf20C18 ~ logk + s(ageAtScan) + sex.o, data=data_k)

setwd("/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/pdfs")

# Comp14  without labels
pdf("CT-Comp14-logkEffect_noLabels_0417.pdf", width = 8, height = 8)
visreg(gm.comp14, 'logk',xlab="",ylab="", overlay=T,cex.axis=2,cex.main=2,points=list(cex=1.03))
dev.off()
# Comp 14 with labels
pdf("CT-Comp14-logkEffect_0417.pdf", width = 8, height = 8)
visreg(gm.comp14, 'logk', overlay=T,ylab="CT Scores in Network 14",xlab="Discount Rate (log k)",cex.axis=1.5,cex.lab=2,cex.main=2,points=list(cex=1.01))
dev.off()

# Comp 15 without labels
pdf("CT-Comp15-logkEffect_noLabels_0417.pdf", width = 8, height = 8)
visreg(gm.comp15, 'logk',xlab="",ylab="", overlay=T,cex.axis=2,cex.main=2,points=list(cex=1.03))
dev.off()
# Comp 15 with labels
pdf("CT-Comp15-logkEffect_0417.pdf", width = 8, height = 8)
visreg(gm.comp15, 'logk', overlay=T,ylab="CT Scores in Network 15",cex.axis=1.5,cex.lab=1.5,cex.main=2,points=list(cex=1.01))
dev.off()

# association with cognitive variables
cor.test(data_k$logk, data_k$Overall_Accuracy, method="pearson")
cor.test(data_k$logk, data_k$F1_Exec_Comp_Res_Accuracy, method="pearson")
cor.test(data_k$logk, data_k$F2_Social_Cog_Accuracy, method="pearson")        
cor.test(data_k$logk, data_k$F3_Memory_Accuracy, method="pearson")

###############################################
# multivariate prediction with NMF components #
###############################################

# baseline covariates
covs1<-"sex.o+ageAtScan"
#covs<-"sex.o+ageAtScan+F1_Exec_Comp_Res_Accuracy"
covs2<-"sex.o+ageAtScan+F1_Exec_Comp_Res_Accuracy+F3_Memory_Accuracy+F2_Social_Cog_Accuracy"
covs2<-"sex.o+ageAtScan+F1_Exec_Comp_Res_Accuracy+F3_Memory_Accuracy"

gamCovs1<-"sex.o+s(ageAtScan)"
#covs<-"sex.o+ageAtScan+F1_Exec_Comp_Res_Accuracy"
gamCovs2<-"sex.o+s(ageAtScan)+F1_Exec_Comp_Res_Accuracy+F3_Memory_Accuracy"

nmfcovs<-"Nmf20C1+Nmf20C2+Nmf20C3+Nmf20C4+Nmf20C5+Nmf20C6+Nmf20C7+Nmf20C8+Nmf20C9+Nmf20C10+Nmf20C11+Nmf20C12+Nmf20C13+Nmf20C14+Nmf20C15+Nmf20C16+Nmf20C18+Nmf20C19+Nmf20C20+"

#nmfcovs<-"Nmf20C3+Nmf20C7+Nmf20C8+Nmf20C9+Nmf20C11+Nmf20C12+Nmf20C13+Nmf20C14+Nmf20C15+Nmf20C18+Nmf20C20+"
onlyNmfCovs<-"Nmf20C1+Nmf20C2+Nmf20C3+Nmf20C4+Nmf20C5+Nmf20C6+Nmf20C7+Nmf20C8+Nmf20C9+Nmf20C10+Nmf20C11+Nmf20C12+Nmf20C13+Nmf20C14+Nmf20C15+Nmf20C16+Nmf20C18+Nmf20C19+Nmf20C20"

# linear models
lmBase1<-lm(as.formula(paste0("logk~",covs1)),data=data_k)
lmBase2<-lm(as.formula(paste0("logk~",covs2)),data=data_k)
lmBaseNMF<-lm(as.formula(paste0("logk~",onlyNmfCovs)),data=data_k)
lmCt1<-lm(as.formula(paste0("logk~",nmfcovs,covs1)),data=data_k)
lmCt2<-lm(as.formula(paste0("logk~",nmfcovs,covs2)),data=data_k)

lmCt<-lm(logk~Nmf20C1+Nmf20C2+Nmf20C3+Nmf20C4+Nmf20C5+Nmf20C6+Nmf20C7+Nmf20C8+Nmf20C9+Nmf20C10+Nmf20C11+Nmf20C12+Nmf20C13+Nmf20C14+Nmf20C15+Nmf20C16+Nmf20C18+Nmf20C19+Nmf20C20+sex.o+ageAtScan+ageSq+F1_Exec_Comp_Res_Accuracy+F3_Memory_Accuracy,data=data_k)

gamCt<-gam(logk~Nmf20C1+Nmf20C2+Nmf20C3+Nmf20C4+Nmf20C5+Nmf20C6+Nmf20C7+Nmf20C8+Nmf20C9+Nmf20C10+Nmf20C11+Nmf20C12+Nmf20C13+Nmf20C14+Nmf20C15+Nmf20C16+Nmf20C18+Nmf20C19+Nmf20C20+sex.o+s(ageAtScan)+F1_Exec_Comp_Res_Accuracy+F3_Memory_Accuracy,data=data_k)

gamBase1<-gam(logk~sex.o+s(ageAtScan)+F1_Exec_Comp_Res_Accuracy+F3_Memory_Accuracy,data=data_k)

# GAMs
gam(Nmf20C14 ~ logk + meduCnbGo1 + F1_Exec_Comp_Res_Accuracy+  s(ageAtScan) + sex.o, data=data)

# model with just demographics + cognition
cor.test(predict(lmBase2),data_k$logk) 

# CT model against baseline of just age and sex
anova(lmCt1,lmBase1)
cor.test(predict(lmCt1),data_k$logk)

# CT model against baseline of just age and sex + cognition
anova(lmCt2,lmBase2)
cor.test(predict(lmCt2),data_k$logk)

# just baseline CT model, no covariates
cor.test(predict(lmBaseNMF),data_k$logk)

par("las"=1)
pdf("NmfCompMultivarPredNoCovariates.pdf", width = 8, height = 8)
par("las"=1)
plot(predict(lmBaseNMF),data_k$logk,main=expression(paste('CT multivariate logk prediction, ',italic(r),' = 0.31')),ylab="Actual logk",xlab="Predicted logk",cex.axis=1.5,cex.lab=1.5,cex.main=1.7,cex=1.01,pch=20,col="blue")
abline(lm(data_k$logk~predict(lmCt1)), col="red")
dev.off()

visreg(lm(predict(lmBaseNMF)~data_k$logk), 'logk', overlay=T, main='Association between CT-Comp18 and logk',ylab="CT Scores in Comp18",cex.axis=1.5,cex.lab=1.5,cex.main=2,points=list(cex=1.01))


visreg(gm.comp18, 'logk', overlay=T, main='Association between CT-Comp18 and logk',ylab="CT Scores in Comp18",cex.axis=1.5,cex.lab=1.5,cex.main=2,points=list(cex=1.01))
dev.off()

### exploring covariates
nmfSigComp<-c("Nmf20C3","Nmf20C7","Nmf20C8","Nmf20C9","Nmf20C11","Nmf20C12","Nmf20C13","Nmf20C14","Nmf20C15","Nmf20C18","Nmf20C20")
gm.mod <- gam(Nmf20C1 ~ logk + s(ageAtScan) + sex.o, data=data_k)
gm.mod1 <- gam(Nmf20C1 ~ logk + s(ageAtScan) + sex.o + averageRating, data=data_k)

# adding T1 rating as covariate
ms5 <-paste0(nmf20Names[-17],"~logk+averageRating+sex.o+s(ageAtScan)")
ms5_model<-regDataGamMultOut(ms5,data_k,nmf20Names[-17],1,2)
write.table(ms5_model,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMFlogkT1CovariateModels.csv",row.names=T,col.names
=T,sep=",")
justSig<-c(3,7:9,11:15,17,19)

# adding Medu as a covariate 
ms6 <-paste0(nmf20Names[-17],"~logk+meduCnbGo1+sex.o+s(ageAtScan)")
ms6_model<-regDataGamMultOut(ms6,data_k,nmf20Names[-17],1,2)
write.table(ms6_model,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMFlogkMeduCovariateModels.csv",row.names=T,col.names
=T,sep=",")

# adding TBV rating as covariate
ms7 <-paste0(nmf20Names[-17],"~logk+tbv+sex.o+s(ageAtScan)")
ms7_model<-regDataGamMultOut(ms7,data_k,nmf20Names[-17],1,2)
write.table(ms7_model,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMFlogkTBVCovariateModels.csv",row.names=T,col.names
=T,sep=",")

# getting FDR-p values just for significant components
ms5_sig <-paste0(nmfSigComp,"~logk+averageRating+sex.o+s(ageAtScan)")
ms5_model_sig<-regDataGamMultOut(ms5_sig,data_k,nmfSigComp,1,2)

# only restricting analysis to T1 imgages with a score of 2
data_T12 <- data_k[which(data_k$averageRating==2),]
msT12<-paste0(nmf20Names[-17],"~logk+sex.o+s(ageAtScan)")
modT12<-regDataGam(msT12,data_T12,nmf20Names[-17],1) # w/o the noise component

# only testing in significant components
msT12_sig<-paste0(nmfSigComp,"~logk+sex.o+s(ageAtScan)")
modT12sig<-regDataGam(msT12_sig,data_T12,nmfSigComp,1)

# sensitivity analysis on subjects not using medication
modOutNoMeds<-regDataGam(ms,data_k_med,nmf20Names[-17],1) # w/o the noise component
sigComp<-rownames(modOutNoMeds[which(modOutNoMeds$fdr.p<.05),])
write.table(modOutNoMeds,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NmfMedicationSensitivity.csv",row.names=T,col.names=T,sep=",")

