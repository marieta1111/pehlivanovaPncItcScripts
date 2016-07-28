# script to run GAM models with Muse regional volume data
require(mgcv)

# read in data
data<-readRDS('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n448_pnc_itc_whole_sample_20160711.rds')
# read in ROI names
museNames<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/volume_muse_names.csv', header=F)

rem<-unique(c(1:116,grep("WM",museNames[,1]),grep("Vent",museNames[,1]),grep("White",museNames[,1]),grep("corpus",museNames[,1]),grep("Brain_Stem",museNames[,1]), grep("CSF",museNames[,1]),grep("fornix",museNames[,1]), grep("Cer",museNames[,1]) ))

museNames1<-museNames[-rem,]

# what about these two?  vol_R_vol_PLIC_CerPed_R vol_L_vol_PLIC_CerPed_R
# from Ted:
#grep("DC",museNames[,1]), grep("InC",museNames[,1]), grep("fornix",museNames[,1]))))] #119 as in atlas

n<-length(museNames1)

# subset data for models
dataIncl<-data[which(data$museRegInclude==1),]

#########################
# vector with model specifications
ms<-paste0(museNames1,"~logk+sex.o+s(ageAtScan)+tbv")

# run gam models and extracts stats for logk
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
