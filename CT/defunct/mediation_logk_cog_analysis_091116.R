# script to run GAM models with JLF regional CT data
require(mgcv)
library(visreg)

# read in data
data<-readRDS('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n452_pnc_itc_whole_sample_20160905.rds')

# extract jlf CT names
jlfCtNames <- names(data)[grep("mprage_jlf_ct",names(data))]

# exclude some regions including subcortical structures
rem<-unique(c(grep("WM",jlfCtNames),grep("Vessel",jlfCtNames), grep("Vent",jlfCtNames),grep("White",jlfCtNames), grep("corpus",jlfCtNames),grep("Brain_Stem",jlfCtNames), grep("CSF",jlfCtNames),grep("fornix",jlfCtNames), grep("Cer",jlfCtNames),grep("DC",jlfCtNames), grep("OpticChiasm", jlfCtNames), grep("Putamen",jlfCtNames),grep("Amygdala",jlfCtNames),grep("Accumbens",jlfCtNames),grep("Caudate",jlfCtNames),grep("Thal",jlfCtNames),grep("Hippocampus",jlfCtNames), grep("Pallidum",jlfCtNames), grep("BasForebr",jlfCtNames)))

jlfCtNames1<- jlfCtNames[-rem]
# number of ROIs
n<-length(jlfCtNames1)

# subset data for models
data_k<-data[which(data$jlfCTInclude_k==1),]
data_a<-data[which(data$jlfCTInclude_alpha==1),]
# adding squared term for logk
data_k$logkSq<-(data_k$logk-mean(data_k$logk))^2
# logk quartiles
data_k$logkQ<-factor(cut(data_k$logk,breaks=quantile(data_k$logk,probs=seq(0,1, by=0.25),na.rm=TRUE),include.lowest=TRUE),labels=c("Q1","Q2","Q3","Q4"))
data_k$logkD<-factor(cut(data_k$logk,breaks=quantile(data_k$logk,probs=seq(0,1, by=0.1),na.rm=TRUE),include.lowest=TRUE),labels=c("D1","D2","D3","D4","D5","D6","D7","D8","D9","D10"))

predK<-readRDS("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/logk98pcs.rds")
data_k$pcaK<-as.vector(predK)
# add R function for regional analysis
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/regDataGam.R")
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/ROIanalysis.R")

#####################################################################
# vectors with model statements, variable of interest is listed first
ms_f1_exra <-paste0(jlfCtNames1,"~F1_Exec_Comp_Res_Accuracy+sex.o+s(ageAtScan)")
gamExra<-regDataGam(ms_f1_exra,data_k,jlfCtNames1,1)
ms_k_exra <-paste0(jlfCtNames1,"~logk+F1_Exec_Comp_Res_Accuracy+sex.o+s(ageAtScan)")
gamKexra<-regDataGam(ms_k_exra,data_k,jlfCtNames1,1)

ms_k<-paste0(jlfCtNames1,"~logk+sex.o+s(ageAtScan)")
ms_ksq<-paste0(jlfCtNames1,"~logk+logkSq+sex.o+s(ageAtScan)")
ms_ksq_exra<-paste0(jlfCtNames1,"~logk+logkSq+F1_Exec_Comp_Res_Accuracy+sex.o+s(ageAtScan)")
ms_k_med1 <-paste0(jlfCtNames1,"~logk+meduCnbGo1+sex.o+s(ageAtScan)")
ms_k_med2 <-paste0(jlfCtNames1,"~logk+meduCnbGo2+sex.o+s(ageAtScan)")

#ms_mEd <-paste0(jlfCtNames1,"~meduCnbGo1+sex.o+s(ageAtScan)")
#ms_smEd <-paste0(jlfCtNames1,"~s(meduCnbGo1)+sex.o+s(ageAtScan)")

# MODELS NEEDED FOR MEDIATION ANALYSES
ms_kI <-paste0("logk~",jlfCtNames1,"+sex.o+s(ageAtScan)")
ms_ecraII <-paste0("F1_Exec_Comp_Res_Accuracy~",jlfCtNames1,"+sex.o+s(ageAtScan)")
ms_III <-paste0("logk~F1_Exec_Comp_Res_Accuracy+",jlfCtNames1,"+sex.o+s(ageAtScan)")

# running models
gam_kI<-regDataGam(ms_kI,data_k,jlfCtNames1,1)
sigI<-row.names(gam_kI[which(gam_kI$fdr.p<0.05),])
saveRDS(gam_kI,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_ecraMediationModel1_Ct_%drois.rds",dim(data_k)[1],n))

gam_ecraII<-regDataGam(ms_ecraII,data_k,jlfCtNames1,1)
sigII<-row.names(gam_ecraII[which(gam_ecraII$fdr.p<0.05),])
saveRDS(gam_ecraII,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_ecraMediationModel2_Ct_%drois.rds",dim(data_k)[1],n))

gam_III<-regDataGam(ms_III,data_k,jlfCtNames1,1)
sigIII<-row.names(gam_III[which(gam_III$fdr.p<0.05),])
saveRDS(gam_III,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_ecraMediationModel3_Ct_%drois.rds",dim(data_k)[1],n))

# L MSFG (reg # 153) -- one region where there is no significant mediation 
summary(gam(as.formula(ms_kI[48]), data=data_k, method="REML"))
summary(gam(as.formula(ms_ecraII[48]), data=data_k, method="REML"))
summary(gam(as.formula(ms_III[48]), data=data_k, method="REML"))

# mediation results
regs<-intersect(sigI,sigII) # regions to look at for mediation
write.csv(regs,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/CT_logk_sigMed.csv")

notMed<-setdiff(sigI,sigII)

# mediation on all regions significant for logk
medResults <-ROIanalysis(sigI ,data_k)

# mediation on intersection of logk and ct sig  
medResultsI <-ROIanalysis(regs,data_k)
saveRDS(medResultsI,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_logkEcraMediation_Ct_SobelZ.rds",dim(data_k)[1]))

# areas of significant mediation by sobel's
sigMed<-jlfCtNames1[which(medResults$sobelFDR<0.05)]

medResults[which(jlfCtNames1=="mprage_jlf_ct_L_MFC"),]

# see which ROIs are significant
gamOut1<-gam_ecraII
gamOut1[which(gamOut1$fdr.p<0.05),]
           
###### random code in preparation for committee meeting
par(las=1)
plot(data_k$logk, data_k$F1_Exec_Comp_Res_Accuracy,pch=16,xlab="logk",ylab="Complex Reasoning Accuracy",col="blue",cex.lab = 1.2,cex.axis = 1.2)
boxplot(F1_Exec_Comp_Res_Accuracy~logkD,data=data_k)
# testing non-linearity
fit1<-lm(F1_Exec_Comp_Res_Accuracy~logk,data=data_k)
fit2<-lm(F1_Exec_Comp_Res_Accuracy~logk+logkSq,data=data_k)
anova(fit1,fit2)

# visreg
# R MORG 146 mprage_jlf_ct_R_MOrG 
# also try mprage_jlf_ct_L_FRP 121

pdf('/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/pdfs/logkL_FRP.pdf')
gm.test <- gam(mprage_jlf_ct_L_FRP ~ logk + s(ageAtScan) + sex.o, data=data_k)
visreg(gm.test, 'logk', overlay=T, main='CT in Left Frontal Pole',ylab="CT")
dev.off()

pdf('/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/pdfs/logkL_FRPs.pdf')
gm.testS <- gam(mprage_jlf_ct_L_FRP ~ s(logk) + s(ageAtScan) + sex.o, data=data_k)
visreg(gm.testS, 'logk', overlay=T, main='CT in Left Frontal Pole',ylab="CT")
dev.off()

pdf('/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/pdfs/logkL_FuG.pdf')
gm.testF <- gam(mprage_jlf_ct_L_FuG ~ logk + s(ageAtScan) + sex.o, data=data_k)
visreg(gm.testF, 'logk', overlay=T, main='CT in Left Fusiform Gyrus',ylab="CT")
dev.off()
 
pdf('/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/pdfs/logkL_FuGS.pdf')
gm.testFs <- gam(mprage_jlf_ct_R_FuG ~ s(logk) + s(ageAtScan) + sex.o, data=data_k)
visreg(gm.testFs, 'logk', overlay=T, main='CT in Left Fusiform Gyrus',ylab="CT")
dev.off()

# COGNITION
pdf('/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/pdfs/ecraL_FRP.pdf')
gm.test <- gam(mprage_jlf_ct_L_FRP ~ F1_Exec_Comp_Res_Accuracy + s(ageAtScan) + sex.o, data=data_k)
visreg(gm.test, 'F1_Exec_Comp_Res_Accuracy', overlay=T, main='CT in Left Frontal Pole',ylab="CT",xlab="Cognitive Reasoning Accuracy")
dev.off()

pdf('/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/pdfs/ecraL_FRPs.pdf')
gm.testS <- gam(mprage_jlf_ct_L_FRP ~ s(F1_Exec_Comp_Res_Accuracy) + s(ageAtScan) + sex.o, data=data_k)
visreg(gm.testS, 'F1_Exec_Comp_Res_Accuracy', overlay=T, main='CT in Left Frontal Pole',ylab="CT",xlab="Cognitive Reasoning Accuracy")
dev.off()

pdf('/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/pdfs/ecraL_MORG.pdf')
gm.test1 <- gam(mprage_jlf_ct_L_MOrG ~ F1_Exec_Comp_Res_Accuracy + s(ageAtScan) + sex.o, data=data_k)
visreg(gm.test1, 'F1_Exec_Comp_Res_Accuracy', overlay=T, main='CT in Medial Orbital Gyrus',ylab="CT",xlab="Cognitive Reasoning Accuracy")
dev.off()

# LOGK-COGNITION
pdf('/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/pdfs/ecraL_logk.pdf')
gm.test <- gam(F1_Exec_Comp_Res_Accuracy ~ logk + s(ageAtScan) + sex.o, data=data_k)
visreg(gm.test, 'logk', overlay=T, main='Discounting-Cognitive Reasoning Accuracy',ylab="CRA",xlab="logk")
dev.off()

pdf('/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/pdfs/ecraL_logkS.pdf')
gm.testS <- gam(F1_Exec_Comp_Res_Accuracy ~ s(logk) + s(ageAtScan) + sex.o, data=data_k)
visreg(gm.testS, 'logk', overlay=T, main='Discounting-Cognitive Reasoning Accuracy',ylab="CRA",xlab="logk")
dev.off()

# mother's education explorations
dataEdu<-data_k[which(data_k$meduCnbGo2!=0),] # 18 NAs, 29 0s
dataEdu$mEdu<-(dataEdu$meduCnbGo1+dataEdu$meduCnbGo2)/2

ms_k <-paste0(jlfCtNames1,"~logk+logSq+sex.o+s(ageAtScan)")
gam1<-regDataGam(ms_k,data_k,jlfCtNames1,1)

ms_med <-paste0(jlfCtNames1,"~mEdu+sex.o+s(ageAtScan)")
gam_med<-regDataGam(ms_med,dataEdu,jlfCtNames1,1)
saveRDS(gam_med,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_Medu_Ct.rds",dim(dataEdu)[1]))

ms_k_med <-paste0(jlfCtNames1,"~logk+mEdu+sex.o+s(ageAtScan)")
ms_med_race <-paste0(jlfCtNames1,"~mEdu+as.factor(race2)+sex.o+s(ageAtScan)")
gam_med_race<-regDataGam(ms_med_race,dataEdu,jlfCtNames1,1)

# MODELS NEEDED FOR MEDIATION ANALYSES
ms_kI <-paste0("logk~",jlfCtNames1,"+sex.o+s(ageAtScan)")
ms_meduII <-paste0("mEdu~",jlfCtNames1,"+sex.o+s(ageAtScan)")
ms_III <-paste0("logk~mEdu+",jlfCtNames1,"+sex.o+s(ageAtScan)")

gam_kI<-regDataGam(ms_kI,dataEdu,jlfCtNames1,1)
sigI<-row.names(gam_kI[which(gam_kI$fdr.p<0.05),])

gam_eduII<-regDataGam(ms_meduII,dataEdu,jlfCtNames1,1)
sigII<-row.names(gam_eduII[which(gam_eduII$fdr.p<0.05),])

gam_III<-regDataGam(ms_III,dataEdu,jlfCtNames1,1)
sigIII<-row.names(gam_III[which(gam_III$fdr.p<0.05),])

medResults <-ROIanalysis(sigI ,dataEdu)


gam_k<-regDataGam(ms_k,dataEdu,jlfCtNames1,1)
gam_med<-regDataGam(ms_med,dataEdu,jlfCtNames1,1)
gam_k_med<-regDataGam(ms_k_med,dataEdu,jlfCtNames1,1)

saveRDS(gamSlogk,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_SlogkJlfCtResults_%drois.rds",dim(data_k)[1],n))

# correlation matrix for cognitive variables
library(corrplot)
M<-cor(data[,1252:1264])
corrplot(M, method="number",type="lower",tl.cex = .6)
corrplot.mixed(M,tl.cex=.5)

# look at correlation matrix with TBV
ctData<-data_k[,jlfCtNames1] # all good CT regions
ctDataSig<- data_k[,row.names(gamLogk[which(gamLogk$fdr.p<0.05),])] # just FDR-significant CT regions
cors<-cor(ctData,data_k$tbv)

barplot(cors[1:98],names.arg=gsub("mprage_jlf_ct_","",row.names(cors)[1:98]))

par(las=2)
barplot(cors[,1],names.arg=gsub("mprage_jlf_ct_","",row.names(cors)),main="CT correlations with TBV")
barplot(cors[,1],names.arg=gsub("mprage_jlf_ct_","",row.names(cors)))





