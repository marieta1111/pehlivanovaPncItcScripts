# script to run GAM models with JLF regional CT data
require(mgcv)

# read in data
data<-readRDS('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n452_pnc_itc_whole_sample_20160905.rds')

# extract jlf CT names
jlfCtNames <- names(data)[grep("mprage_jlf_ct",names(data))]

# exclude some regions including subcortical structures
rem<-unique(c(grep("WM",jlfCtNames),grep("Vessel",jlfCtNames), grep("Vent",jlfCtNames),grep("White",jlfCtNames), grep("corpus",jlfCtNames),grep("Brain_Stem",jlfCtNames), grep("CSF",jlfCtNames),grep("fornix",jlfCtNames), grep("Cer",jlfCtNames),grep("DC",jlfCtNames), grep("OpticChiasm", jlfCtNames), grep("Putamen",jlfCtNames),grep("Amygdala",jlfCtNames),grep("Accumbens",jlfCtNames),grep("Caudate",jlfCtNames),grep("Thal",jlfCtNames),grep("Hippocampus",jlfCtNames), grep("Pallidum",jlfCtNames), grep("BasForebr",jlfCtNames)))

jlfCtNames1<- jlfCtNames[-rem]

n<-length(jlfCtNames1)

# subset data for models
data_k<-data[which(data$jlfCTInclude_k==1),]
# adding squared term for logk
data_k$logkSq<-(data_k$logk-mean(data_k$logk))^2
# logk quartiles
data_k$logkQ<-factor(cut(data_k$logk,breaks=quantile(data_k$logk,probs=seq(0,1, by=0.25),na.rm=TRUE),include.lowest=TRUE),labels=c("Q1","Q2","Q3","Q4"))
data_k$logkD<-factor(cut(data_k$logk,breaks=quantile(data_k$logk,probs=seq(0,1, by=0.1),na.rm=TRUE),include.lowest=TRUE),labels=c("D1","D2","D3","D4","D5","D6","D7","D8","D9","D10"))

########################################################################################
# AvgT1Rating (data quality)
# Race
# Race + Medu (Can use mEdu score from TP1 only)
# Would know what relationship of all 3 variables (race, quality, medu) are w/ LogK
# Also would ask whether mEdu has a relationship w/ CT over and above race (which we know has a relationship w/ CT)
# Environment score from timepoint 1-- called goassessEnvSes in data release.  Note that for this score race MUST be in model. (like mEdu).

x<-"mprage_jlf_ct_L_STG~logk+meduCnbGo1+as.factor(race2)+sex.o+s(ageAtScan)"
test<-gam(as.formula(x), data=data_k, method="REML")
foo<-summary(test)

# add R function for regional analysis
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/regDataGam.R")
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/regDataGamMultOut.R")

#####################################################################
# vectors with model statements, variable of interest is listed first
# then running model with custom function

# basic model with just logk
ms0 <-paste0(jlfCtNames1,"~logk+sex.o+s(ageAtScan)")
ms0_model<-regDataGamMultOut(ms0,data_k,jlfCtNames1,1,1)
# regions significant for logk
sigReg<-row.names(ms0_model[which(ms0_model$fdr.p<0.05),])
sigRegN<-which(jlfCtNames1 %in% sigReg)

#### just mother's education
ms1 <-paste0(jlfCtNames1,"~logk+meduCnbGo1+sex.o+s(ageAtScan)")
ms1_model<-regDataGamMultOut(ms1,data_k,jlfCtNames1,1,2)
# get standardized coefficients
ttt<-paste0("data_k$",jlfCtNames1,"[which(!is.na(data_k$meduCnbGo1))]")
sds<-unname(sapply(ttt,function(x) sd(eval(parse(text=x))))) 
ms1_model$logk_st <- ms1_model$coef*sd(data_k$logk[which(!is.na(data_k$meduCnbGo1))])/sds
ms1_model$medu_st <- ms1_model$coef.1*sd(data_k$meduCnbGo1[which(!is.na(data_k$meduCnbGo1))])/sds

ms1_sig<-ms1_model[sigRegN,] # in all models, logk is still significant, but weaker association
ms1_sig[which(ms1_sig[,9]<0.05),] # in these regions medu is sig. associated with CT

#### mother's education + race
ms2 <-paste0(jlfCtNames1,"~logk+meduCnbGo1+as.factor(race2)+sex.o+s(ageAtScan)")
ms2_model<-regDataGamMultOut(ms2,data_k,jlfCtNames1,1,4)

ms2_sig<-ms2_model[sigRegN,] 

#### race
ms3 <-paste0(jlfCtNames1,"~logk+as.factor(race2)+sex.o+s(ageAtScan)")
ms3_model<-regDataGamMultOut(ms3,data_k,jlfCtNames1,1,3)

ms3_sig<-ms3_model[sigRegN,]

#### Environment score + race
ms4 <-paste0(jlfCtNames1,"~logk+goassessEnvSes+as.factor(race2)+sex.o+s(ageAtScan)")
ms4_model<-regDataGamMultOut(ms4,data_k,jlfCtNames1,1,4)

ms4_sig<-ms4_model[sigRegN,]

#### average T1 rating
ms5 <-paste0(jlfCtNames1,"~logk+averageRating+sex.o+s(ageAtScan)")
ms5_model<-regDataGamMultOut(ms5,data_k,jlfCtNames1,1,2)
ms5_sig<-ms5_model[sigRegN,]

x<-"mprage_jlf_ct_R_CO~logk+meduCnbGo1+as.factor(race2)+sex.o+s(ageAtScan)"
test<-gam(as.formula(x), data=data_k, method="REML")
foo<-summary(test)

#### Mother's education + Race, w/o logk
ms6<-paste0(jlfCtNames1,"~meduCnbGo1+as.factor(race2)+sex.o+s(ageAtScan)")
ms6_model<-regDataGamMultOut(ms6,data_k,jlfCtNames1,1,3)

# regions where medu is significant, only 7
row.names(ms6_model[ms6_model[,5]<0.05,])
# get standardized coefficients
ttt<-paste0("data_k$",jlfCtNames1,"[which(!is.na(data_k$meduCnbGo1))]")
sds<-unname(sapply(ttt,function(x) sd(eval(parse(text=x)))))
ms6_model$medu_st <- ms6_model$coef*sd(data_k$meduCnbGo1[which(!is.na(data_k$meduCnbGo1))])/sds
ms6_model$aaVw_st <- ms6_model$coef.1*sd(data_k$meduCnbGo1[which(!is.na(data_k$meduCnbGo1))])/sds




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
