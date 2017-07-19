library(ANTsR)

# read in jlf ROI numbers
labelsMP<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/jlfCtRoiNumbers.csv',header=T) 

# read in data from model of interest
modData<-readRDS('/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/n427_logkJlfCtResults_98rois.rds')
# just subset to significant results
sigRes<-modData[which(modData$fdr.p<0.05),]
row.names(sigRes)

# identify region numbers of significant areas
sigRes$roiNum <- sapply(row.names(sigRes), function(x) {
        ids<-grep(x,as.character(labelsMP$region),ignore.case=TRUE)
        l<-labelsMP$regNumber[ids]
        return(l)
})

numSig<-sigRes$roiNum

# read in template
temp<-antsImageRead('/data/joy/BBL/templates/MNI/jlf/Younger24/mniTemplateJLF_Labels.nii.gz',3)
test<-antsImageClone(temp)

# count numbers of voxels in mask per region
boo <- data.frame(nums=integer(0), numVox= integer(0))
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
	tval<-sigRes$tval[which(numSig==j)]
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

antsImageWrite(test,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/jlf/nifti/logkCtImageSigAreas_tvals.nii.gz")

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

