#####################################
# most analyses for CT-DD paper, MP #
#####################################

require(mgcv)
library(visreg)
library(ggplot2)
library(reshape2)
library(corrplot)
library(voxel)

# add R functions 
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/regDataGam.R") # regional analysis
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/regDataGamMultOut.R")

# load data
load("/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/DataK_n427_030717.Rda")
# additional age vars
data_k$ageScanY=data_k$ageAtScan/12 # age in years
# data with top and bottom logk quartiles
dataQ14<-data_k[data_k$logkQ.o %in% c('Q1','Q4'),]

# NMF names 
nmf20Names <-names(data_k)[1467:1486] 

### demographic + cognitive analyses for manuscript
# descriptive statistics for k
mean(data_k$k)
sd(data_k$k)

# associations with age and sex in GAM models
gmAgeSex1 <- gam(logk ~ s(ageAtScan) + sex.o, data=data_k)
gmAgeSex2 <- gam(logk ~ s(ageAtScan) + sex.o + s(ageAtScan, by=sex.o), data=data_k)

# association with cognitive variables (we ended up using partial cor. analyses in matlab)
# code for partial correlations in NMF_matlab_calculations_011717.m
cor.test(data_k$logk, data_k$Overall_Accuracy, method="pearson")
cor.test(data_k$logk, data_k$F1_Exec_Comp_Res_Accuracy, method="pearson")
cor.test(data_k$logk, data_k$F2_Social_Cog_Accuracy, method="pearson")
cor.test(data_k$logk, data_k$F3_Memory_Accuracy, method="pearson")

### main analyses of NMF components
# create vector of model statements
ms<-paste0(nmf20Names[-17],"~logk+sex.o+s(ageAtScan)")
# running basic models for association between NMF components and logk
modOut<-regDataGam(ms,data_k,nmf20Names[-17],1) # w/o the noise component
sigComp<-rownames(modOut[which(modOut$fdr.p<.05),])

# correlations between logk and CT components (also ended up using partial correlations in matlab script)
## code for partial correlations in NMF_matlab_calculations_011717.m
cors_p<-cor(data_k$logk,data_k[,nmf20Names[-17]],method="pearson")

# saving results
#saveRDS(modOut,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/logkCtNmf20Results_19comp.rds")

### visualizations of CT-logk effects
gm.comp14 <- gam(Nmf20C14 ~ logk + s(ageAtScan) + sex.o, method="REML", data=data_k)      
gm.comp15 <- gam(Nmf20C15 ~ logk + s(ageAtScan) + sex.o, method="REML", data=data_k)

setwd("/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/pdfs")

# Comp 14 with labels
pdf("CT-Comp14-logkEffect_0417.pdf", width = 8, height = 8)
visreg(gm.comp14, 'logk', overlay=T,ylab="CT Scores in Network 14",xlab="Discount Rate (log k)",cex.axis=1.5,cex.lab=2,cex.main=2,points=list(cex=1.01))
dev.off()
# Comp 15 with labels
pdf("CT-Comp15-logkEffect_0417.pdf", width = 8, height = 8)
visreg(gm.comp15, 'logk', overlay=T,ylab="CT Scores in Network 15",xlab="Discount Rate (log k)",cex.axis=1.5,cex.lab=1.5,cex.main=2,points=list(cex=1.01))
dev.off()

### age effects 
# interaction models in signficant components
mods<-paste0(sigComp,"~ ageScanY*logk + logk + sex.o + ageScanY + ageSq")
pvals = NA
for (i in 1:length(mods)){
        foo<-lm(as.formula(mods[i]),data=data_k)
        tab<-summary(foo)
        pvals[i]<-tab$coefficients[6,4]
}
median(pvals)
min(pvals)
max(pvals)

### Interaction plots (Figure 5)
# network 14
gm14 <- gam(Nmf20C14 ~ s(ageScanY) + logkQ.o, data=dataQ14, method="REML")
plot14<-plotGAM(gm14, "ageScanY", "logkQ.o", orderedAsFactor = T, rawOrFitted = "raw", plotCI = T)

pdf("LogkQuartiles-Comp14_0517_2.pdf", width = 8, height = 8)
plot14 + theme_bw() + ylab("CT Scores in Network 14") + xlab("Age") +
       ggtitle("") + scale_size(range=c(8,20))+
       theme(plot.title = element_text(hjust = 0.5,size=22,face="bold"), panel.background = element_blank(), panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),axis.title=element_text(size=24),legend.text = element_text(size = 16, hjust = 3, vjust = 3),
       legend.title=element_text(size=16),axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)) +
       scale_colour_discrete(name="DD")+scale_x_continuous(breaks=seq(10,24,2))
dev.off()

# network 15 
gm15 <- gam(Nmf20C15 ~ s(ageScanY) + logkQ.o, data=dataQ14, method="REML")
plot15<-plotGAM(gm15, "ageScanY", "logkQ.o", orderedAsFactor = T, rawOrFitted = "raw", plotCI = T)

pdf("LogkQuartiles-Comp15_0517.pdf", width = 8, height = 8)
plot15 + theme_bw() + ylab("CT Scores in Network 15") + xlab("Age") +
       ggtitle("") + scale_size(range=c(8,20))+
       theme(plot.title = element_text(hjust = 0.5,size=22,face="bold"), panel.background = element_blank(), panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),axis.title=element_text(size=24),legend.text = element_text(size = 16, hjust = 3, vjust = 3),
       legend.title=element_text(size=16),axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)) +
       scale_colour_discrete(name="DD")+scale_x_continuous(breaks=seq(10,24,2))
dev.off()

### multivariate prediction with NMF components 
# baseline covariates
covs1<-"sex.o+ageAtScan"
covs2<-"sex.o+ageAtScan+F1_Exec_Comp_Res_Accuracy+F3_Memory_Accuracy"
nmfcovs<-"Nmf20C1+Nmf20C2+Nmf20C3+Nmf20C4+Nmf20C5+Nmf20C6+Nmf20C7+Nmf20C8+Nmf20C9+Nmf20C10+Nmf20C11+Nmf20C12+Nmf20C13+Nmf20C14+Nmf20C15+Nmf20C16+Nmf20C18+Nmf20C19+Nmf20C20+"

# linear models
lmBase1<-lm(as.formula(paste0("logk~",covs1)),data=data_k)
lmBase2<-lm(as.formula(paste0("logk~",covs2)),data=data_k)
lmCt1<-lm(as.formula(paste0("logk~",nmfcovs,covs1)),data=data_k)
lmCt2<-lm(as.formula(paste0("logk~",nmfcovs,covs2)),data=data_k)

# model with just demographics
cor.test(predict(lmBase1), data_k$logk)
# model with just demographics + cognition
cor.test(predict(lmBase2),data_k$logk) 
# CT model against baseline of just age and sex
anova(lmCt1,lmBase1)
cor.test(predict(lmCt1),data_k$logk)
# CT model against baseline of just age and sex + cognition
anova(lmCt2,lmBase2)
cor.test(predict(lmCt2),data_k$logk)

### plotting multivariate prediction (Figure 6)
# plotting values predicted from model with CT + demographics against actual log k values
gamCovs1<-"sex.o+s(ageAtScan)"
nmfcovs<-"Nmf20C1+Nmf20C2+Nmf20C3+Nmf20C4+Nmf20C5+Nmf20C6+Nmf20C7+Nmf20C8+Nmf20C9+Nmf20C10+Nmf20C11+Nmf20C12+Nmf20C13+Nmf20C14+Nmf20C15+Nmf20C16+Nmf20C18+Nmf20C19+Nmf20C20+"
nmfModel2<-gam(as.formula(paste0("logk~",nmfcovs,gamCovs1)),data=data_k,method="REML")
data_k$nmfPred<-predict(nmfModel2)

plotNmfDemoModel <- lm(logk~nmfPred,data=data_k)
cor.test(data_k$logk,data_k$nmfPred,method="pearson")

# plot + saving plot (Note: in paper draft, dimensions might be slightly different 
pdf("multivariatePrediction_0717.pdf", width = 8, height = 8)
par("las"=1)
visreg(plotNmfDemoModel, 'nmfPred',xlab="Predicted Discount Rate (Log k)",ylab="Actual Discount Rate (Log k)",cex.axis=1.3,cex.lab=1.5,cex.main=1.8,points=list(cex=0.7),line.par=list(col="blue"))
legend(x='bottomright', legend='r = 0.33, p < 0.0001',bty = "n",cex=1.235)
dev.off()

### sensitivity analyses
## adding Medu as a covariate 
ms6 <-paste0(nmf20Names[-17],"~logk+meduCnbGo1+sex.o+s(ageAtScan)")
ms6_model<-regDataGamMultOut(ms6,data_k,nmf20Names[-17],1,2)
write.table(ms6_model,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMFlogkMeduCovariateModels.csv",row.names=T,col.names
=T,sep=",")

## adding TBV rating as covariate
ms7 <-paste0(nmf20Names[-17],"~logk+tbv+sex.o+s(ageAtScan)")
ms7_model<-regDataGamMultOut(ms7,data_k,nmf20Names[-17],1,2)
write.table(ms7_model,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMFlogkTBVCovariateModels.csv",row.names=T,col.names
=T,sep=",")

## sensitivity analysis on subjects not using medication
# medications data 
meds<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/go1_go2_go3_medications_92816_ts_st_tsFinal.csv',header=T,sep=',')
names(meds)[names(meds)=="study"]<-"timepoint"
meds$medsComb<-meds$psychoactiveMedPsych_20160930+meds$psychoactiveMedMedical_20160930

# merge with meds data
comb<-merge(data_k, meds[,c("bblid","scanid","timepoint","medsComb")], by=c("bblid","timepoint"), all.x=T)
data_k_med <-comb[which(comb$medsComb==0),] # subsetting data with subjects NOT using medications

modOutNoMeds<-regDataGam(ms,data_k_med,nmf20Names[-17],1) # w/o the noise component
medSigComp<-rownames(modOutNoMeds[which(modOutNoMeds$fdr.p<.05),])
write.table(modOutNoMeds,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NmfMedicationSensitivity.csv",row.names=T,col.names=T,sep=",")


if(FALSE) {
# adding T1 rating as covariate, not currently in manuscript
ms5 <-paste0(nmf20Names[-17],"~logk+averageRating+sex.o+s(ageAtScan)")
ms5_model<-regDataGamMultOut(ms5,data_k,nmf20Names[-17],1,2)
write.table(ms5_model,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMFlogkT1CovariateModels.csv",row.names=T,col.names
=T,sep=",")
justSig<-c(3,7:9,11:15,17,19)

# getting FDR-p values just for significant components
ms5_sig <-paste0(sigComp,"~logk+averageRating+sex.o+s(ageAtScan)")
ms5_model_sig<-regDataGamMultOut(ms5_sig,data_k,sigComp,1,2)

# only restricting analysis to T1 imgages with a score of 2
data_T12 <- data_k[which(data_k$averageRating==2),]
msT12<-paste0(nmf20Names[-17],"~logk+sex.o+s(ageAtScan)")
modT12<-regDataGam(msT12,data_T12,nmf20Names[-17],1) # w/o the noise component
write.table(modT12,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMFlogkOnly2Scores.csv",row.names=T,col.names
=T,sep=",")

# only testing in significant components
msT12_sig<-paste0(sigComp,"~logk+sex.o+s(ageAtScan)")
modT12sig<-regDataGam(msT12_sig,data_T12,nmfSigComp,1)
}


