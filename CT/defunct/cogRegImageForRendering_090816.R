library(ANTsR)

# read in jlf ROI numbers
labelsMP<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/jlfCtRoiNumbers.csv',header=T) 
# read in data from model of interest
modData<-readRDS('/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n380_Medu_Ct.rds')
# file names to save nifti output
fileName<-"/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/nifti/Ct_med.nii.gz"

# just subset to significant results
sigRes<-modData[which(modData$fdr.p<0.05),]
names<-row.names(sigRes)

# only for mediation results
sigRes<-modData
#theNames<-read.csv("/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/CT_logk_sigMed.csv",header=F)
names<-c("mprage_jlf_ct_R_CO", "mprage_jlf_ct_L_CO", "mprage_jlf_ct_R_Cun", "mprage_jlf_ct_R_FRP","mprage_jlf_ct_L_FRP", "mprage_jlf_ct_L_FuG","mprage_jlf_ct_L_Gre","mprage_jlf_ct_R_IOG","mprage_jlf_ct_L_MFC","mprage_jlf_ct_R_MOrG","mprage_jlf_ct_L_MOrG","mprage_jlf_ct_R_MPrG","mprage_jlf_ct_R_MTG","mprage_jlf_ct_L_OrIFG","mprage_jlf_ct_R_Pins","mprage_jlf_ct_L_POrG","mprage_jlf_ct_L_PP","mprage_jlf_ct_L_PrG","mprage_jlf_ct_R_PT","mprage_jlf_ct_R_STG","mprage_jlf_ct_L_STG","mprage_jlf_ct_R_TMP","mprage_jlf_ct_L_TMP")

# identify region numbers of significant areas
sigRes$roiNum <- sapply(names, function(x) {
        ids<-grep(x,as.character(labelsMP$region),ignore.case=F)
        l<-labelsMP$regNumber[ids[1]]
        return(l)
})

# for visual comparison with OASIS30Labels.csv
cbind(names,sigRes$roiNum)

numSig<-sigRes$roiNum

# to create a function, would need: 
# 1. names of regions that overlap
# 2. what value to assign to each region

# read in template
temp<-antsImageRead('/data/joy/BBL/templates/MNI/jlf/Younger24/mniTemplateJLF_Labels.nii.gz',3)
test<-antsImageClone(temp)

# count numbers of voxels in mask per region
boo <- data.frame(nums=integer(0), numVox=integer(0))
for (j in 1:207){
	print(j)
	boo[j,1]<-j
	boo[j,2]<-length(test[test==as.numeric(j)])
}

# identify areas which exist in mask, i.e. are not 0 voxels
numExist<-boo$nums[which(boo$numVox!=0)]

# identify voxels which are in mask, but are not sig,
# i.e. should be set to 0
numSet0<-setdiff(numExist, numSig)
 
# set the Significant regions to t-value
for (j in numSig){
	print(j)
#	t-value for associations
	tval<-sigRes$tval[which(numSig==j)]
#  	sobel's z for mediation
#	tval<-sigRes$sobelZ[which(numSig==j)]
	print(tval)
	print(test[test==as.numeric(j)])
	test[test==as.numeric(j)]<-tval
}

# set all other existing regions to 0
for (j in numSet0){
	print(j)
	print(test[test==as.numeric(j)])
	test[test==as.numeric(j)]<-0
}

antsImageWrite(test,fileName)





if (FALSE){
# read in ROI labels
labels<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/regionalVolumes/OASIS30Labels.csv')

# identify region numbers
roiNums <- sapply(row.names(sigRes), function(x) {
        ids<-grep(substr(x,17,nchar(x)+1),labels1$Label.Name,ignore.case=TRUE)
        # determine if left or right
        if (length(grep("_L",x))==1){
        l<-labels1$Label.Number[ids[2]]
        } else if (length(grep("_L",x))==0) {
        l<-labels1$Label.Number[ids[1]]
        }
        return(l)
})
}

