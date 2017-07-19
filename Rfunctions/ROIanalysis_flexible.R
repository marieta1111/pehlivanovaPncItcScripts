# code to run ROI mediation analysis  for specific clusters/regions 
# INPUT
# roiNames:  string vector of names of ROIs for which to do mediation
# data:      dataset
# medVar:    mediating variable to be tested, as character
# covs:      covariates in the model, of the type: "sex.o + s(ageAtScan)"

ROIanalysis_flexible<- function(roiNames,data,medVar,covs){

require(mgcv)

covs<-"sex.o + s(ageAtScan)"

nn<-length(roiNames) # number of ROIs to look at

# create data frame to store output
ROIresults<-data.frame(roiName=double(),Mod1_coef=double(),Mod1_SE=double(),Mod1_p=double(),Mod2_coef=double(),Mod2_SE=double(),Mod2_p=double(),Mod3_medCoef=double(),Mod3_medSE=double(),Mod3_medP=double(),Mod3_ivCoef=double(),Mod3_ivSE=double(),Mod3_ivP=double(),sobelZ=double(),sobelP=double(),crit1=logical(),crit2=logical(),crit3=logical(),crit4=logical(),allCrit=double(),sobelFDR=double())

# do that for these clusters (all or subset)
for (i in 1:nn){

	# assign brain region 
	y<-as.name(eval(roiNames[i]))

	# mediator
	med<-as.name(eval(medVar))

	# ROI
	ROIresults[i,1]<-i    

	#------------------ code for mediation analysis and sobel's test
	# Model 1: logk = CT + sex + age 
	mdl1<-gam(as.formula(paste0("logk ~ eval(y) +", covs)),data=data, method="REML") 
	coefs1<-summary(mdl1)$p.table[2,1:4]
	ROIresults[i,2]<-coefs1[1] # coef
	ROIresults[i,3]<-coefs1[2] # SE
	ROIresults[i,4]<-coefs1[4] # p-value

	# Model 2: Mediator = CT + sex + age
	mdl2<-gam(as.formula(paste0("eval(med) ~ eval(y) +",covs)), data=data, method="REML")
	coefs2<-summary(mdl2)$p.table[2,1:4]
	a<-coefs2[1]
	Sa<-coefs2[2]
	ROIresults[i,5]<-coefs2[1] # coef
	ROIresults[i,6]<-coefs2[2] # SE
	ROIresults[i,7]<-coefs2[4] # p-value

	# Model 3: logk = Mediator + CT + sex + age 
	mdl3<-gam(as.formula(paste0("logk ~ eval(med) + eval(y) +",covs)), data=data, method="REML")
	coefs3<-c(summary(mdl3)$p.table[2,c(1,2,4)],summary(mdl3)$p.table[3,c(1,2,4)])
	b<-coefs3[1] # Mediator coefficient in predicting logk, controlling for CT
	Sb<-coefs3[2] # SE of the b coefficient
	ROIresults[i,8]<-coefs3[1] # coef for mediator
	ROIresults[i,9]<-coefs3[2] # SE for mediator 
	ROIresults[i,10]<-coefs3[3] # p-value for mediator
	ROIresults[i,11]<-coefs3[4] # coef for brain
	ROIresults[i,12]<-coefs3[5] # SE for brain
	ROIresults[i,13]<-coefs3[6] # p-value for brain

	sobelZ<-a*b/sqrt((b^2)*(Sa^2) + (a^2)*(Sb^2) + (Sa^2)*(Sb^2))
	ROIresults[i,14]<-sobelZ
	ROIresults[i,15]<-2*pnorm(-abs(sobelZ))
	# mediation criteria
	ROIresults[i,16]<-ROIresults[i,4]<0.05 # criterion 1
	ROIresults[i,17]<-ROIresults[i,7]<0.05 # criterion 2
	ROIresults[i,18]<-ROIresults[i,10]<0.05 # criterion 3
	ROIresults[i,19]<-abs(ROIresults[i,2])>abs(ROIresults[i,11]) # criterion 4
	ROIresults[i,20]<-sum(ROIresults[i,16],ROIresults[i,17],ROIresults[i,18],ROIresults[i,19]) # all criteria

} # for loop 

ROIresults[,21]<-p.adjust(ROIresults[,15], method = "fdr", n =nn)

return(ROIresults)
} # end of function

