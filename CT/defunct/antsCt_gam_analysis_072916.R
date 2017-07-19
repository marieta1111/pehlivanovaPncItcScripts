# script to run GAM models with Muse regional volume data
require(mgcv)

# read in data
data<-readRDS('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n448_pnc_itc_whole_sample_20160730.rds')

# CT variables
ctVars <- names(data)[grep("CT_", colnames(data))]

# remove some ROIs that we don't need, such as WM, ventricles, csf, etc.
rem<-unique(c(grep("WM",ctVars),grep("Vent",ctVars),grep("White",ctVars),grep("capsule",ctVars),grep("corpus",ctVars),grep("BrainStem",ctVars), grep("CSF",ctVars),grep("fornix",ctVars), grep("Cer",ctVars),grep("DC",ctVars) ))

ctVars1<-ctVars[-rem]

# what about these?
# pallidum? 
# "CT_LeftBasalForebrain" <-- includes nucleus accumbens                                    
# "CT_RightBasalForebrain"  
#grep(grep("InC",museNames[,1])))] #119 as in atlas

n<-length(ctVars1)

# subset data for models
dataIncl<-data[which(data$scanIncludeAll==1),]

#########################
# vector with model specifications
ms<-paste0(ctVars1,"~logk+sex.o+s(ageAtScan)")

# run GAM models and extract stats *for logk* and n for individual models
model.results <- lapply(ms, function(x) {
    foo <- summary(gam(as.formula(x), data=dataIncl, method="REML"))
    return(c(foo$n, foo$p.table[2,]))
    #return(foo$p.table[2,])
})

# convert results to data frame
gamOut1 <- data.frame(matrix(unlist(model.results), nrow=n, byrow=T))
names(gamOut1) <-c("n","coef","se","tval","p")
row.names(gamOut0) <- ctVars1

# multiple comparison testing
gamOut1[,6]<-p.adjust(gamOut1[,5], method = "fdr", n = n)
names(gamOut1)[6]<-"fdr.p"

# saving results
#saveRDS(gamOut1,sprintf("/data/joy/BBL/projects/pehlivanovaPncItc/antsCT/logkAntsCtResults_120rois.rds",dim(dataIncl)[1]))
saveRDS(gamOut1,"/data/joy/BBL/projects/pehlivanovaPncItc/antsCT/logkAntsCtResults_116rois_wo_tbv.rds")

# see which ROIs are significant
gamOut1[which(gamOut1$fdr.p<0.05),]
