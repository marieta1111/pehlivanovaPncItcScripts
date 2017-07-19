# script to run GAM models with JLF regional CT data
require(mgcv)
library(visreg)

# read in data
data<-read.table('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n452_pnc_itc_whole_sample_20160905_wSubsData.csv',header=T,sep=',')

# read in jlf CT names for 98 "good" regions
jlfCtNames1 <-as.matrix(read.csv("/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/jlf/jlfCtGoodRegions.csv",header=F))

# read in CT PCs
pcsData<-read.csv("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/ctPCs.csv",header=T,sep=",")

# number of ROIs
n<-length(jlfCtNames1)

# subset data for models
data_k<-data[which(data$jlfCTInclude_k==1),]
# adding squared term for logk
data_k$logkSq<-(data_k$logk-mean(data_k$logk))^2
# logk quartiles
data_k$logkQ<-factor(cut(data_k$logk,breaks=quantile(data_k$logk,probs=seq(0,1, by=0.25),na.rm=TRUE),include.lowest=TRUE),labels=c("Q1","Q2","Q3","Q4"))
data_k$logkD<-factor(cut(data_k$logk,breaks=quantile(data_k$logk,probs=seq(0,1, by=0.1),na.rm=TRUE),include.lowest=TRUE),labels=c("D1","D2","D3","D4","D5","D6","D7","D8","D9","D10"))

# merge datasets with pcs
data<-merge(data_k,pcsData,by=c("bblid"), all.x=T)

#predK<-readRDS("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/logk98pcs.rds")
#data_k$pcaK<-as.vector(predK)
# add R function for regional analysis
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/regDataGam.R")
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/ROIanalysis.R")
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/ROIanalysis_flexible.R")
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/ROIanalysis_flexible2.R")
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/ROIanalysis_CtDep.R")

#####################################################################
# MODELS NEEDED FOR MEDIATION ANALYSES

medVar<-"meduCnbGo1" # variable for mediation

ms_kI <-paste0("logk~",jlfCtNames1,"+sex.o+s(ageAtScan)")
ms_II <-paste0(medVar,"~",jlfCtNames1,"+sex.o+s(ageAtScan)")
ms_III <-paste0("logk~",medVar,"+",jlfCtNames1,"+sex.o+s(ageAtScan)")

# running models
gam_kI<-regDataGam(ms_kI,data_k,jlfCtNames1,1)
sigI<-row.names(gam_kI[which(gam_kI$fdr.p<0.05),])

gam_II<-regDataGam(ms_II,data_k,jlfCtNames1,1)
sigII<-row.names(gam_II[which(gam_II$fdr.p<0.05),])

gam_III<-regDataGam(ms_III,data_k,jlfCtNames1,1)
sigIII<-row.names(gam_III[which(gam_III$fdr.p<0.05),])

# L MSFG (reg # 153) -- one region where there is no significant mediation 
summary(gam(as.formula(logk~BAS_dr+mprage_jlf_ct_L_TTG+sex.o+s(ageAtScan)), data=data_k, method="REML"))
summary(gam(as.formula(ms_ecraII[48]), data=data_k, method="REML"))
summary(gam(as.formula(ms_III[48]), data=data_k, method="REML"))

# mediation results
regs<-intersect(sigI,sigII) # regions to look at for mediation
write.csv(regs,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/CT_logk_sigMed.csv")

# mediation on all regions significant for logk
medResults<-ROIanalysis(sigI,data_k)
medResults <-ROIanalysis_flexible2(sigI ,data_k, "meduCnbGo1")
medResults <-ROIanalysis_flexible2(sigI ,data_k, "goassessEnvSes")
medResults2 <-ROIanalysis_flexible(sigI ,data_k, "meduCnbGo1")
medResults2 <-ROIanalysis_flexible(sigI ,data_k, "meduCnbGo1","as.factor(race2)+ sex.o + s(ageAtScan)")
saveRDS(medResults,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_logkBASdrMediation_Ct_SobelZ.rds",dim(data_k)[1]))
# areas of significant mediation by sobel's
sigMed<-jlfCtNames1[which(medResults$sobelFDR<0.05)]
           
##############################################################################
# Exploring effect of race: associations between variables of interest and CT, 
# with and without controlling for race 

resLogk<-ROIanalysis_CtDep(sigI,data_k,"logk")
resCR<-ROIanalysis_CtDep(sigI,data_k,"F1_Exec_Comp_Res_Accuracy")
resMedu<-ROIanalysis_CtDep(sigI,data_k,"meduCnbGo1")
resBasDr<-ROIanalysis_CtDep(sigI,data_k,"BAS_dr")

##############################################################################
# multivariate mediation with PCs
covs<-"sex.o+ageAtScan+meduCnbGo1+F1_Exec_Comp_Res_Accuracy"
covs<-"sex.o+ageAtScan+meduCnbGo1"
pccovs<-"PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+"
lmBase<-lm(as.formula(paste0("logk~",covs)),data=data)
lmCt<-lm(as.formula(paste0("logk~",pccovs,covs)),data=data)
anova(lmBase,lmCt)


mdl2<-summary(gam(as.formula(mprage_jlf_ct_L_TMP ~ logk + sex.o + s(ageAtScan) + as.factor(race2)), data=data_k, method="REML"))

# correlation matrix for cognitive variables
library(corrplot)
M<-cor(data[,1252:1264])
corrplot(M, method="number",type="lower",tl.cex = .6)
corrplot.mixed(M,tl.cex=.5)






