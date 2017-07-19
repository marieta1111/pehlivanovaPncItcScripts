# create nifti image for the association between logk and the NMF components

library(ANTsR)

####################
# logk-NMF ct associations
# combine warped images first
img<-antsImageRead('/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/all20CompsNMFWarped.nii.gz',4) #this 4d set of merged images
mask<-antsImageRead('/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/n1359_antsCT_mask_warped_bin.nii.gz',3)

# create matrix, rows are components, columns are voxels
seed.mat<-timeseries2matrix(img, mask)

# find for each voxels which component has highest loading
whichCompStr <-apply(seed.mat,2,which.max) # however, some of those well all zeros, need to remove
foo <-apply(seed.mat,2,sum) # this is sum across column, if 0, entire column is 0
whichCompStr[which(foo==0)]<-0 # assign 0-columns to 0

# writing that to an image
newImg<-antsImageClone(mask)  #prep for writing out image
newImg[mask==1]<-as.matrix(whichCompStr)
antsImageWrite(newImg,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/NmfHardPartition_fromWarpedThrImages.nii.gz")

# assign t-values for components where logk is Fdr-sig. associated with logk
pvalMap<-antsImageClone(mask)
pvalVoxVector<-whichCompStr
pvalVoxVector[which(pvalVoxVector %in% c(1,2,4,5,6,10,16,17,19))]<-0
pvalVoxVector[which(pvalVoxVector==3)]<--2.83
pvalVoxVector[which(pvalVoxVector==7)]<--2.62
pvalVoxVector[which(pvalVoxVector==8)]<--2.54
pvalVoxVector[which(pvalVoxVector==9)]<--2.21
pvalVoxVector[which(pvalVoxVector==11)]<--2.47
pvalVoxVector[which(pvalVoxVector==12)]<--2.57
pvalVoxVector[which(pvalVoxVector==13)]<--2.74
pvalVoxVector[which(pvalVoxVector==14)]<--4.76
pvalVoxVector[which(pvalVoxVector==15)]<--4.14
pvalVoxVector[which(pvalVoxVector==18)]<--2.91
pvalVoxVector[which(pvalVoxVector==20)]<--2.91

pvalMap[mask==1]<-as.matrix(pvalVoxVector)
antsImageWrite(pvalMap,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/pvalsIn11sigComps_fromwarpedImages.nii.gz")

if(FALSE) {
############################
# for CR mediation analysis
img<-antsImageRead('/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/all20CompsNMFWarped.nii.gz',4) #this 4d set of merged images
mask<-antsImageRead('/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/n1359_antsCT_mask_warped_bin.nii.gz',3)

# create matrix, rows are components, columns are voxels
seed.mat<-timeseries2matrix(img, mask)

# find for each voxels which component has highest loading
whichCompStr <-apply(seed.mat,2,which.max) # however, some of those well all zeros, need to remove
foo <-apply(seed.mat,2,sum) # this is sum across column, if 0, entire column is 0
whichCompStr[which(foo==0)]<-0 # assign 0-columns to 0

# writing that to an image
newImg<-antsImageClone(mask)  #prep for writing out image
newImg[mask==1]<-as.matrix(whichCompStr)
antsImageWrite(newImg,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/NmfHardPartition_fromRawThr004.nii.gz")

# assign t-values for components where logk is Fdr-sig. associated with logk
pvalMap<-antsImageClone(mask)
pvalVoxVector<-whichCompStr
pvalVoxVector[which(pvalVoxVector %in% c(1,2,4,5,6,7,9,10,13,16,17,18,19))]<-0
pvalVoxVector[which(pvalVoxVector==3)]<--2.89
pvalVoxVector[which(pvalVoxVector==8)]<--2.34
pvalVoxVector[which(pvalVoxVector==11)]<--3.03
pvalVoxVector[which(pvalVoxVector==12)]<--2.24
pvalVoxVector[which(pvalVoxVector==14)]<--4.09
pvalVoxVector[which(pvalVoxVector==15)]<--4.21
pvalVoxVector[which(pvalVoxVector==20)]<--3.06

pvalMap[mask==1]<-as.matrix(pvalVoxVector)
antsImageWrite(pvalMap,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/pvalsCrMediation.nii.gz")

#######################################
# for Mother's education mediation analysis
img<-antsImageRead('/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/all20CompsNMFWarped.nii.gz',4) #this 4d set of merged images
mask<-antsImageRead('/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/n1359_antsCT_mask_warped_bin.nii.gz',3)

# create matrix, rows are components, columns are voxels
seed.mat<-timeseries2matrix(img, mask)

# find for each voxels which component has highest loading
whichCompStr <-apply(seed.mat,2,which.max) # however, some of those well all zeros, need to remove
foo <-apply(seed.mat,2,sum) # this is sum across column, if 0, entire column is 0
whichCompStr[which(foo==0)]<-0 # assign 0-columns to 0

# assign t-values for components where Medu mediation is sig (but w/o FDR correction)
pvalMap<-antsImageClone(mask)
pvalVoxVector<-whichCompStr
pvalVoxVector[which(pvalVoxVector %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,16,17,18,19))]<-0
pvalVoxVector[which(pvalVoxVector==14)]<--2.12
pvalVoxVector[which(pvalVoxVector==15)]<--2.19
pvalVoxVector[which(pvalVoxVector==20)]<--1.97

pvalMap[mask==1]<-as.matrix(pvalVoxVector)
antsImageWrite(pvalMap,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/tvalsMeduMediationNonCorrected.nii.gz")

####################
# CR-NMF ct associations
# combine warped images first
img<-antsImageRead('/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/all20CompsNMFWarped.nii.gz',4) #this 4d set of merged images
mask<-antsImageRead('/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/n1359_antsCT_mask_warped_bin.nii.gz',3)

# create matrix, rows are components, columns are voxels
seed.mat<-timeseries2matrix(img, mask)

# find for each voxels which component has highest loading
whichCompStr <-apply(seed.mat,2,which.max) # however, some of those well all zeros, need to remove
foo <-apply(seed.mat,2,sum) # this is sum across column, if 0, entire column is 0
whichCompStr[which(foo==0)]<-0 # assign 0-columns to 0

# writing that to an image
#newImg<-antsImageClone(mask)  #prep for writing out image
#newImg[mask==1]<-as.matrix(whichCompStr)
#antsImageWrite(newImg,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/NmfHardPartition_fromWarpedThrImages.nii.gz")

# assign t-values for components where logk is Fdr-sig. associated with logk
pvalMap<-antsImageClone(mask)
pvalVoxVector<-whichCompStr
pvalVoxVector[which(pvalVoxVector %in% c(1,4,6,7,9,13,16,17,18))]<-0
pvalVoxVector[which(pvalVoxVector==2)]<-2.33
pvalVoxVector[which(pvalVoxVector==3)]<-3.34
pvalVoxVector[which(pvalVoxVector==5)]<-3.12
pvalVoxVector[which(pvalVoxVector==8)]<-2.55
pvalVoxVector[which(pvalVoxVector==10)]<-2.38
pvalVoxVector[which(pvalVoxVector==11)]<-3.53
pvalVoxVector[which(pvalVoxVector==12)]<-2.46
pvalVoxVector[which(pvalVoxVector==14)]<-7.08
pvalVoxVector[which(pvalVoxVector==15)]<-7.19
pvalVoxVector[which(pvalVoxVector==19)]<-3.61
pvalVoxVector[which(pvalVoxVector==20)]<-3.58

pvalMap[mask==1]<-as.matrix(pvalVoxVector)
antsImageWrite(pvalMap,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/pvalsCR11sigComps_fromwarpedImages.nii.gz")

#########################
# Medu-NMF CT associations
# combine warped images first
img<-antsImageRead('/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/all20CompsNMFWarped.nii.gz',4) #this 4d set of merged images
mask<-antsImageRead('/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/n1359_antsCT_mask_warped_bin.nii.gz',3)

# create matrix, rows are components, columns are voxels
seed.mat<-timeseries2matrix(img, mask)

# find for each voxels which component has highest loading
whichCompStr <-apply(seed.mat,2,which.max) # however, some of those well all zeros, need to remove
foo <-apply(seed.mat,2,sum) # this is sum across column, if 0, entire column is 0
whichCompStr[which(foo==0)]<-0 # assign 0-columns to 0

# assign t-values for components where logk is Fdr-sig. associated with logk
pvalMap<-antsImageClone(mask)
pvalVoxVector<-whichCompStr
pvalVoxVector[which(pvalVoxVector %in% c(1,4,6,7,9,13,16,17,18))]<-0
pvalVoxVector[which(pvalVoxVector==2)]<-2.33
pvalVoxVector[which(pvalVoxVector==3)]<-3.34
pvalVoxVector[which(pvalVoxVector==5)]<-3.12
pvalVoxVector[which(pvalVoxVector==8)]<-2.55
pvalVoxVector[which(pvalVoxVector==10)]<-2.38
pvalVoxVector[which(pvalVoxVector==11)]<-3.53
pvalVoxVector[which(pvalVoxVector==12)]<-2.46
pvalVoxVector[which(pvalVoxVector==14)]<-7.08
pvalVoxVector[which(pvalVoxVector==15)]<-7.19
pvalVoxVector[which(pvalVoxVector==19)]<-3.61
pvalVoxVector[which(pvalVoxVector==20)]<-3.58

pvalMap[mask==1]<-as.matrix(pvalVoxVector)
antsImageWrite(pvalMap,"/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20images/pvalsCR11sigComps_fromwarpedImages.nii.gz")
}
