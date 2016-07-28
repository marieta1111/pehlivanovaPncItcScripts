# script to run GAM model with FA DTI data
require(mgcv)

# read in data
data<-readRDS('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n448_pnc_itc_whole_sample_20160711.rds')
# read in ROI names
faNames<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/dti_fa_names.csv', header=F)

n<-dim(faNames)[1]

# preallocate output
gamOut<-matrix(NA,nrow=n,ncol=5)
dimnames(gamOut) <- list(faNames[,1], c("coef","se","tval","p","fdr.p"))

dataIncl<-data[which(data$faDtiInclude==1),]

for (i in 1:n){
        print(faNames[i,])
        x<-paste(faNames[i,],"~logk+sex.o+s(ageAtScan)", sep="")
        fit<-summary(gam(as.formula(x), data=dataIncl, method="REML"))
        gamOut[i,1:4]<-fit$p.table[2,]
}

# mult. comparison testing
gamOut[,5]<-p.adjust(gamOut[,4], method = "fdr", n = n)

# saving results
saveRDS(gamOut,"/data/joy/BBL/projects/pehlivanovaPncItc/diffusion/dtiFaTractsResults_noTbv.rds" )
