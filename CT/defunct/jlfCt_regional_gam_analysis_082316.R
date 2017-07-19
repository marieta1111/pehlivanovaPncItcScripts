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

#####################################################################
# vectors with model statements, variable of interest is listed first
ms_k_tbv <-paste0(jlfCtNames1,"~logk+sex.o+s(ageAtScan)+tbv")
ms_k     <-paste0(jlfCtNames1,"~logk+sex.o+s(ageAtScan)")
ms_a_tbv <-paste0(jlfCtNames1,"~logAlpha+sex.o+s(ageAtScan)+tbv")
ms_a     <-paste0(jlfCtNames1,"~logAlpha+sex.o+s(ageAtScan)")
ms_sk_tbv<-paste0(jlfCtNames1,"~s(logk)+sex.o+s(ageAtScan)+tbv")
ms_sk    <-paste0(jlfCtNames1,"~s(logk)+sex.o+s(ageAtScan)")
ms_sa_tbv<-paste0(jlfCtNames1,"~s(logAlpha)+sex.o+s(ageAtScan)+tbv")
ms_sa    <-paste0(jlfCtNames1,"~s(logAlpha)+sex.o+s(ageAtScan)")

ms_sk_ecra    <-paste0(jlfCtNames1,"~s(logk)+F1_Exec_Comp_Res_Accuracy+sex.o+s(ageAtScan)")
gamSlogkEcra   <-regDataGam(ms_sk_ecra,data_k,jlfCtNames1,2)

ms_k_ecra     <-paste0(jlfCtNames1,"~logk+F1_Exec_Comp_Res_Accuracy+sex.o+s(ageAtScan)")
gamLogkEcra   <-regDataGam(ms_k_ecra,data_k,jlfCtNames1,1)

# GAM model function to run models and extract stats for logk or logAlpha
gamFunc <- function(ms,dataset,varNames,s) {
	# s is a splines variable, if 1, no splines on logk
        # if 2, splines on logk
	n<-length(varNames)
	model.results<-lapply(ms, function(x) {
    	foo<-summary(gam(as.formula(x), data=dataset, method="REML"))
    	# takes stats for first covariate, specifically
        if (s==1) {
                return(c(foo$n,foo$p.table[2,]))
        } else if (s==2) {
                return(c(foo$n,foo$s.table[1,]))
        }
        })
	# convert results to data frame
	gamOut1<-data.frame(matrix(unlist(model.results), nrow=n, byrow=T))
	if (s==1) {
                 names(gamOut1) <-c("n","coef","se","tval","p")
        } else if (s==2) {
                 names(gamOut1) <-c("n","edf","rfDF","F","p")
        }
	row.names(gamOut1)<-varNames
	# multiple comparison testing
	gamOut1[,6]<-p.adjust(gamOut1[,5], method = "fdr", n = n)
	gamOut1[,7]<-ms
        names(gamOut1)[6:7]<-c("fdr.p","model")
	return(gamOut1)
}

# run GAM models
gamLogkTbv <-gamFunc(ms_k_tbv,data_k,jlfCtNames1,1)
gamLogk    <-gamFunc(ms_k,data_k,jlfCtNames1,1)
gamSlogkTbv<-gamFunc(ms_sk_tbv,data_k,jlfCtNames1,2)
gamSlogk   <-gamFunc(ms_sk,data_k,jlfCtNames1,2)
# THERE IS AN OVERFITTING PROBLEM WITH ALPHA AND REGION # 10 mprage_jlf_ct_L_Pallidum
# Can't use the gamFunc for that reason, will exclude this region
gamLogaTbv<-gamFunc(ms_a_tbv[-10],data_a,jlfCtNames1[-10],1)
gamLoga	  <-gamFunc(ms_a[-10],data_a,jlfCtNames1[-10],1)
gamSlogaTbv<-gamFunc(ms_sa_tbv[-10],data_a,jlfCtNames1[-10],2)
gamSloga   <-gamFunc(ms_sa[-10],data_a,jlfCtNames1[-10],2)

# saving results
saveRDS(gamLogkTbv,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_logkJlfCtResults_%drois_wTBV.rds",dim(data_k)[1],n))
saveRDS(gamLogk,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_logkJlfCtResults_%drois.rds",dim(data_k)[1],n))
saveRDS(gamLogaTbv,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_logAlphaJlfCtResults_%drois_wTBV.rds",dim(data_a)[1],n-1))
saveRDS(gamLoga,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_logAlphaJlfCtResults_%drois.rds",dim(data_a)[1],n-1))
saveRDS(gamSlogkTbv,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_SlogkJlfCtResults_%drois_wTBV.rds",dim(data_k)[1],n))
saveRDS(gamSlogk,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n%d_SlogkJlfCtResults_%drois.rds",dim(data_k)[1],n))


# see which ROIs are significant
gamOut1<-gamLogk
gamOut1[which(gamOut1$fdr.p<0.05),]

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
