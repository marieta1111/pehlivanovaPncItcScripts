#create nifti image for the association between logk and the NMF components

library(ANTsR)

############################
# for CR mediation analysis
img<-antsImageRead('/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/all20CompsNMFWarped.nii.gz',4) #this 4d set of merged images for each comp.
mask<-antsImageRead('/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/n1359_antsCT_mask_warped_bin.nii.gz',3)

# create matrix: rows are components, columns are voxels
seed.mat<-timeseries2matrix(img, mask)

# find, for each voxel, which component has highest loading
whichCompStr <-apply(seed.mat,2,which.max) # however, some of those are all zeros, need to remove
foo <-apply(seed.mat,2,sum) 		   # this is sum of loadings across column; if 0, entire column is 0
whichCompStr[which(foo==0)]<-0 		   # assign 0-columns to 0

# writing that to an image where every voxel is assigned to one component
newImg<-antsImageClone(mask)               # prep for writing out image
newImg[mask==1]<-as.matrix(whichCompStr)   # put assigned values back in the image	
antsImageWrite(newImg,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/NmfHardPartition_fromRawThr004.nii.gz")

# assign t-values for components where logk is FDR-sig. associated with logk
# these are t-values but I've mistakenly named it "pval"
pvalMap<-antsImageClone(mask)
pvalVoxVector<-whichCompStr						       # a vector that contains # of Comp for each voxel
pvalVoxVector[which(pvalVoxVector %in% c(1,2,4,5,6,7,9,10,13,16,17,18,19))]<-0 # assign all components that you don't want to 0
pvalVoxVector[which(pvalVoxVector==3)]<--2.89				       # for the remaining components, assign p- or t-value
pvalVoxVector[which(pvalVoxVector==8)]<--2.34
pvalVoxVector[which(pvalVoxVector==11)]<--3.03
pvalVoxVector[which(pvalVoxVector==12)]<--2.24
pvalVoxVector[which(pvalVoxVector==14)]<--4.09
pvalVoxVector[which(pvalVoxVector==15)]<--4.21
pvalVoxVector[which(pvalVoxVector==20)]<--3.06

pvalMap[mask==1]<-as.matrix(pvalVoxVector)				       # essentially this is a replaced image with t-values or 0s instead of comp numbers
antsImageWrite(pvalMap,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/pvalsCrMediation.nii.gz")

