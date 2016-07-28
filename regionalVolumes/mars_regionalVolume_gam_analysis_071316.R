# script to run GAM models with Mars regional volume data
require(mgcv)

# read in data
data<-readRDS('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n448_pnc_itc_whole_sample_20160711.rds')
# read in ROI names
marsNames<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/volume_mars_names.csv', header=F)

# remove first row: ICV
marsNames1<- marsNames[-1,]
n<-length(marsNames1)

# subset data for models
dataIncl<-data[which(data$marsRegInclude==1),]

#########################
# vector with model specifications
ms<-paste0(marsNames1,"~logk+sex.o+s(ageAtScan)+tbv")

# run gam models and extracts stats for logk
model.results <- lapply(ms, function(x) {
    foo <- summary(gam(as.formula(x), data=dataIncl, method="REML"))
    return(foo$p.table[2,])
})

# convert results to data frame
gamOut1 <- data.frame(matrix(unlist(model.results), nrow=n, byrow=T))
names(gamOut1) <-c("coef","se","tval","p")
row.names(gamOut1) <- marsNames1

# multiple comparison testing
gamOut1[,5]<-p.adjust(gamOut1[,4], method = "fdr", n = n)
names(gamOut1)[5]<-"fdr.p"

# saving results
saveRDS(gamOut1, sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/regionalVolumes/n%d_logkMarsResults.rds",dim(dataIncl)[1]))


# USING A LOOP

# preallocate output
#gamOut<-matrix(NA,nrow=n,ncol=5)
#dimnames(gamOut) <- list(marsNames1, c("coef","se","tval","p","fdr.p"))

#for (i in 1:n){
 #       print(as.character(marsNames1[i]))
 #       x<-paste0(marsNames1[i],"~logk+sex.o+s(ageAtScan)+tbv")
 #       fit<-summary(gam(as.formula(x), data=dataIncl, method="REML"))
 #       gamOut[i,1:4]<-fit$p.table[2,]
#}

