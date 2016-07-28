#####################################################
# 7/20/16: script to subset data for MIDAS testing  #
#####################################################

# read in full sample data
data <- read.csv("/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n448_pnc_itc_whole_sample_20160711.csv", header=T, sep=",")

# make a new sex indicator variable 1 vs. 0
data$sexInd <- data$sex-1 

# subset only by rows needed for MIDAS-RAVENS analysis
subData <- data[which(data$scanIncludeAll==1),c("ravensPath","logk","sexInd","ageAtScan","tbv")]

# save .csv to use for MIDAS
write.csv(subData, '/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/n425_ravens_for_midas_072016.csv', row.names=F, quote=F)

 
