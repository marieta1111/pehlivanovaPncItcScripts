# code to run test association between various predcitors of CT, with or w/o race  
ROIanalysis_CtDep<- function(roiNames,data,var1){
require(mgcv)

# roiNames:  string vector of names of ROIs for which to do mediation
# data:      dataset
# var1:      variable to be tested for association with CT 

nn<-length(roiNames) # number of ROIs to look at

# create data frame to store output
ROIresults<-data.frame(roiName=double(),Mod1_coef=double(),Mod1_SE=double(),Mod1_p=double(),Mod2_coef=double(),Mod2_SE=double(),Mod2_p=double(),sobelFDR=double())

# do that for these clusters (all or subset)
for (i in 1:nn){

	# assign brain region 
	y<-as.name(eval(roiNames[i]))

	# variable of interest
	var<-as.name(eval(var1))

	# ROI
	ROIresults[i,1]<-i    

	# Model 1: CT = [var] + age + sex
	mdl1<-gam(as.formula(paste0("eval(y) ~ eval(var) + sex.o + s(ageAtScan)")),data=data, method="REML") 
	coefs1<-summary(mdl1)$p.table[2,1:4]
	ROIresults[i,2]<-coefs1[1] # coef
	ROIresults[i,3]<-coefs1[2] # SE
	ROIresults[i,4]<-coefs1[4] # p-value

	# Model 2: CT = [var] + age + sex + race
	mdl2<-gam(as.formula(paste0("eval(y) ~ eval(var) + sex.o + s(ageAtScan) + as.factor(race2)")), data=data, method="REML")
	coefs2<-summary(mdl2)$p.table[2,1:4]
	ROIresults[i,5]<-coefs2[1] # coef
	ROIresults[i,6]<-coefs2[2] # SE
	ROIresults[i,7]<-coefs2[4] # p-value

} # for loop
 
ROIresults[,8]<-p.adjust(ROIresults[,7], method = "fdr", n=nn)

return(ROIresults)
} # end of function

