# script to run GAM models with JLF regional CT data
require(mgcv)

# read in data
data<-readRDS('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n452_pnc_itc_whole_sample_20160905.rds')

# extract jlf CT names
jlfCtNames <- names(data)[grep("mprage_jlf_ct",names(data))]

# exclude some regions including subcortical structures
rem<-unique(c(grep("WM",jlfCtNames),grep("Vessel",jlfCtNames), grep("Vent",jlfCtNames),grep("White",jlfCtNames), grep("corpus",jlfCtNames),grep("Brain_Stem",jlfCtNames), grep("CSF",jlfCtNames),grep("fornix",jlfCtNames), grep("Cer",jlfCtNames),grep("DC",jlfCtNames), grep("OpticChiasm", jlfCtNames), grep("Putamen",jlfCtNames),grep("Amygdala",jlfCtNames),grep("Accumbens",jlfCtNames),grep("Caudate",jlfCtNames),grep("Thal",jlfCtNames),grep("Hippocampus",jlfCtNames), grep("Pallidum",jlfCtNames), grep("BasForebr",jlfCtNames)))

jlfCtNames1<- jlfCtNames[-rem]

# from Ted:
#grep("DC",museNames[,1]), grep("InC",museNames[,1]), grep("fornix",museNames[,1]))))] #119 as in atlas

n<-length(jlfCtNames1)

# subset data for models
data_k<-data[which(data$jlfCTInclude_k==1),]
data_a<-data[which(data$jlfCTInclude_alpha==1),]
# adding squared term for logk
data_k$logkSq<-(data_k$logk-mean(data_k$logk))^2
# logk quartiles
data_k$logkQ<-factor(cut(data_k$logk,breaks=quantile(data_k$logk,probs=seq(0,1, by=0.25),na.rm=TRUE),include.lowest=TRUE),labels=c("Q1","Q2","Q3","Q4"))
data_k$logkD<-factor(cut(data_k$logk,breaks=quantile(data_k$logk,probs=seq(0,1, by=0.1),na.rm=TRUE),include.lowest=TRUE),labels=c("D1","D2","D3","D4","D5","D6","D7","D8","D9","D10"))

#####################################################################
# vectors with model statements, variable of interest is listed first
ms_oe     <-paste0(jlfCtNames1,"~Overall_Efficiency+sex.o+s(ageAtScan)")
ms_oa     <-paste0(jlfCtNames1,"~Overall_Accuracy+sex.o+s(ageAtScan)")
ms_os     <-paste0(jlfCtNames1,"~Overall_Speed+sex.o+s(ageAtScan)")
ms_f1_exra <-paste0(jlfCtNames1,"~F1_Exec_Comp_Res_Accuracy+sex.o+s(ageAtScan)")
ms_f2_sca <-paste0(jlfCtNames1,"~F2_Social_Cog_Accuracy+sex.o+s(ageAtScan)")
ms_f3_ma <-paste0(jlfCtNames1,"~F3_Memory_Accuracy+sex.o+s(ageAtScan)")
ms_f1_cre <-paste0(jlfCtNames1,"~F1_Complex_Reasoning_Efficiency+sex.o+s(ageAtScan)")
ms_f2_me <-paste0(jlfCtNames1,"~F2_Memory.Efficiency+sex.o+s(ageAtScan)")
ms_f3_ee <-paste0(jlfCtNames1,"~F3_Executive_Efficiency+sex.o+s(ageAtScan)")
ms_f4_sce <-paste0(jlfCtNames1,"~F4_Social_Cognition_Efficiency+sex.o+s(ageAtScan)")
ms_f1_ss <-paste0(jlfCtNames1,"~F1_Slow_Speed+sex.o+s(ageAtScan)")
ms_f2_fs <-paste0(jlfCtNames1,"~F2_Fast_Speed+sex.o+s(ageAtScan)")
ms_f3_ms <-paste0(jlfCtNames1,"~F3_Memory_Speed+sex.o+s(ageAtScan)")
ms_mEd <-paste0(jlfCtNames1,"~meduCnbGo1+sex.o+s(ageAtScan)")
ms_smEd <-paste0(jlfCtNames1,"~s(meduCnbGo1)+sex.o+s(ageAtScan)")

ms_f1_exra <-paste0(jlfCtNames1,"~F1_Exec_Comp_Res_Accuracy+sex.o+s(ageAtScan)")
ms_k_exra <-paste0(jlfCtNames1,"~logk+F1_Exec_Comp_Res_Accuracy+sex.o+s(ageAtScan)")
ms_k<-paste0(jlfCtNames1,"~logk+sex.o+s(ageAtScan)")
ms_ksq<-paste0(jlfCtNames1,"~logk+logkSq+sex.o+s(ageAtScan)")
ms_ksq_exra<-paste0(jlfCtNames1,"~logk+logkSq+F1_Exec_Comp_Res_Accuracy+sex.o+s(ageAtScan)")
ms_k_med1 <-paste0(jlfCtNames1,"~logk+meduCnbGo1+sex.o+s(ageAtScan)")
ms_k_med2 <-paste0(jlfCtNames1,"~logk+meduCnbGo2+sex.o+s(ageAtScan)")

# add R function for regional analysis
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/regDataGam.R")

# running models
gamOE<-regDataGam(ms_oe,data_k,jlfCtNames1,1)
gamOA<-regDataGam(ms_oa,data_k,jlfCtNames1,1)
gamOS<-regDataGam(ms_os,data_k,jlfCtNames1,1)
gamf1_exra<-regDataGam(ms_f1_exra,data_k,jlfCtNames1,1)
gamf2_sca<-regDataGam(ms_f2_sca,data_k,jlfCtNames1,1)
gamf3_ma<-regDataGam(ms_f3_ma,data_k,jlfCtNames1,1)
gamf1_cre<-regDataGam(ms_f1_cre,data_k,jlfCtNames1,1)
gamf2_me<-regDataGam(ms_f2_me,data_k,jlfCtNames1,1)
gamf3_ee<-regDataGam(ms_f3_ee,data_k,jlfCtNames1,1)
gamf4_sce<-regDataGam(ms_f4_sce,data_k,jlfCtNames1,1)
gamf1_ss<-regDataGam(ms_f1_ss,data_k,jlfCtNames1,1)
gamf2_fs<-regDataGam(ms_f2_fs,data_k,jlfCtNames1,1)
gamf3_ms<-regDataGam(ms_f3_ms,data_k,jlfCtNames1,1)
gam_mEd<-regDataGam(ms_mEd,data_k,jlfCtNames1,1)
gam_smEd<-regDataGam(ms_smEd,data_k,jlfCtNames1,2)

gamf1_exra<-regDataGam(ms_f1_exra,data_k,jlfCtNames1,1)
gamK_exra<-regDataGam(ms_k_exra,data_k,jlfCtNames1,1)
gamK<-regDataGam(ms_k,data_k,jlfCtNames1,1)
gamKsq<-regDataGam(ms_ksq,data_k,jlfCtNames1,1)
gamKsq_exra<-regDataGam(ms_ksq_exra,data_k,jlfCtNames1,1)
gam_k_med1<-regDataGam(ms_k_med1,data_k,jlfCtNames1,1)
gam_k_med2<-regDataGam(ms_k_med2,data_k,jlfCtNames1,1)

# see which ROIs are significant
gamOut1<-gam_k_med1
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
gm.test <- gam(mprage_jlf_ct_L_FuG ~ logk + s(ageAtScan) + sex.o, data=data_k)
gm.test1 <- gam(mprage_jlf_ct_L_FuG ~ F1_Exec_Comp_Res_Accuracy + s(ageAtScan) + sex.o, data=data_k)
visreg(gm.test, 'logk', overlay=T, main='CT in L fusiform gyrus') 
visreg(gm.test1, 'F1_Exec_Comp_Res_Accuracy', overlay=T, main='CT in L fusiform gyrus')

# mother's education explorations
dataEdu<-data_k[which(data_k$meduCnbGo2!=0),] # 18 NAs, 29 0s
dataEdu$mEdu<-(dataEdu$meduCnbGo1+dataEdu$meduCnbGo2)/2

ms_k <-paste0(jlfCtNames1,"~logk+sex.o+s(ageAtScan)")
ms_med <-paste0(jlfCtNames1,"~mEdu+sex.o+s(ageAtScan)")
ms_k_med <-paste0(jlfCtNames1,"~logk+mEdu+sex.o+s(ageAtScan)")

gam_k<-regDataGam(ms_k,dataEdu,jlfCtNames1,1)
gam_med<-regDataGam(ms_med,dataEdu,jlfCtNames1,1)
gam_k_med<-regDataGam(ms_k_med,dataEdu,jlfCtNames1,1)

# saving results
saveRDS(gamOE,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_overEffCtResults_%drois.rds",dim(data_k)[1],n))
saveRDS(gamOA,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_overAccCtResults_%drois.rds",dim(data_k)[1],n))
saveRDS(gamOS,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_overSpeedCtResults_%drois.rds",dim(data_k)[1],n))
saveRDS(gamf1_exra,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_ecraCtResults_%drois.rds",dim(data_k)[1],n))
saveRDS(gamf2_sca,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_scaCtResults_%drois.rds",dim(data_k)[1],n))
saveRDS(gamf3_ma,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_maCtResults_%drois.rds",dim(data_k)[1],n))
saveRDS(gamf1_cre,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_creCtResults_%drois.rds",dim(data_k)[1],n))
saveRDS(gamf2_me,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_meCtResults_%drois.rds",dim(data_k)[1],n))
saveRDS(gamf3_ee,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_eeCtResults_%drois.rds",dim(data_k)[1],n))
saveRDS(gam_mEd,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_mEd1CtResults_%drois.rds",dim(data_k)[1],n))


saveRDS(gamLogk,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_logkJlfCtResults_%drois.rds",dim(data_k)[1],n))
saveRDS(gamLogaTbv,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_logAlphaJlfCtResults_%drois_wTBV.rds",dim(data_a)[1],n-1))
saveRDS(gamLoga,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_logAlphaJlfCtResults_%drois.rds",dim(data_a)[1],n-1))
saveRDS(gamSlogkTbv,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_SlogkJlfCtResults_%drois_wTBV.rds",dim(data_k)[1],n))
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




if(FALSE){
# USING A LOOP
# preallocate output
gamOut<-matrix(NA,nrow=n,ncol=5)
dimnames(gamOut) <- list(jlfCtNames1, c("coef","se","tval","p","fdr.p"))

for (i in c(1:9,11:length(jlfCtNames1))){
        #print(as.character(jlfCtNames1))
        print(i)
	x<-paste0(jlfCtNames1[i],"~logAlpha+sex.o+s(ageAtScan)+tbv")
        fit<-summary(gam(as.formula(x), data=data_a, method="REML"))
        gamOut[i,1:4]<-fit$p.table[2,]
}
model.results <- lapply(ms, function(x) {
    foo <- summary(gam(as.formula(x), data=dataIncl, method="REML"))
    return(foo$p.table[2,])
})

# convert results to data frame
gamOut1 <- data.frame(matrix(unlist(model.results), nrow=n, byrow=T))
names(gamOut1) <-c("coef","se","tval","p")
row.names(gamOut1) <- museNames1

# multiple comparison testing
gamOut1[,5]<-p.adjust(gamOut1[,4], method = "fdr", n = n)
names(gamOut1)[5]<-"fdr.p"

# saving results
saveRDS(gamOut1,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/regionalVolumes/n%d_logkMuseResults_119rois.rds",dim(dataIncl)[1]))

# see which ROIs are significant
gamOut1[which(gamOut1$fdr.p<0.05),]

}}
