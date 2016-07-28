#----------------------------------------------------------------------------------------------------------------------------
# JUNE16: this script combines and formats various datasets relevant to the
#
# PNC-ITC project; dataset were compiled by Kosha on 6/7/16
#
# It produces several datasets for further analysis or processing:
#
# 1. dataset with n=453 subjects who have full ITC data to 1) pull out scan dates and 2) to caclulate k-values from ITC data
# 	/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n453_ITC_unique_bblids.csv
#	
#	Currently commented out
#
# 2. dataset with n=448 subjects with processed ITC (*only* usable k-values), paths for imaging analysis, and imaging inclusion variables
#	/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/nXXX_pnc_itc_whole_sample_YYYYMMDD.csv
#	/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/nXXX_pnc_itc_whole_sample_YYYYMMDD.rds
#
#  	Saves data with current date when running code
#------------------------------------------------------------------------------------------------------------------------------

##############################################
#           READING IN ORIGINAL DATA         #
##############################################

# general demographics, includes age at CNB and scans
demo <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n9498_demographics_go1_go2_113015.csv',head=T,sep=',')
# ITC: this has n=653 and not 654 as the file name states
itc <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n656_scales_mfq_asrm_panas_spq_ddisc_rdisc_wolfq.csv',head=T,sep=',')
# CNB measures
cnb <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/GO1_GO2_CNB_Factor_Scores_Non-Age-Regressed.csv',head=T,sep=',')
# data release
data <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n1601_go1_datarel_020716.csv',head=T,sep=',')
# imaging QA measures
qa_go1 <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n1601_t1RawManualQA.csv',header=T,sep=',')
qa_go2 <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n368_t1RawManualQA_GO2.csv',header=T,sep=',')

#---- MIGHT ADD THESE LATER WHEN LOOKING AT THAT
# substance use, n=9700, 9267 unique
#subs <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/pnc_substance_cnb_go1_go2_10215.csv', head=T,sep=',')

##############################################
#             FORMATTING ITC DATA            #
##############################################

# order itc data by bblid and timepoint
itc_o <- itc[with(itc, order(bblid, timepoint)), ]
# just sample with observations with full ITC data
itc_full_o <- itc_o[which(!is.na(itc_o$kddisc_q_01)),c(1:4,207:275)]
# keep just one record per person, FIRST observation/scan
# 32 subjects have duplicate ITC data
itc_o_single <- itc_full_o[!duplicated(itc_full_o$bblid),-3]
# save data to pull out scan dates
#write.csv(itc_o_single,'/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n453_ITC_unique_bblids.csv',row.names=FALSE)

##############################################
#          FORMATTING SCAN DATES DATA        #
##############################################

#************************************************************************************************************
#                                  HERE RUN TO OBTAIN DATES FROM CFN PATHS                                  *
# /data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/sampleAssembly/ravens_get_scan_dates.sh *
#************************************************************************************************************

# scan dates from CFN cluster, 3 subjects out of 453 unique ITC subjects don't have RAVENS folders: 88677,89115,109980
dates_cfn<-read.table('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/list_of_ravens_scan_dates.txt',header=T,sep=',')
# dates of CNB, scans, from Kosha: n1601_dates_cnb_imaging.csv
dates_k <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n1601_dates_cnb_imaging.csv',head=T,sep=',')

#------ dates from CFN paths --------------------
# scan dates from CNF imaging paths
dates_cfn$scan1_d_cfn <- as.Date(substr(dates_cfn$dates,1,8), "%Y%m%d")
dates_cfn$scan2_d_cfn <- as.Date(substr(dates_cfn$dates,15,22), "%Y%m%d") # did checks to make sure it's right
# 344 have both scans
sum(dates_cfn$scan2_d_cfn[!is.na(dates_cfn$scan2_d_cfn)]>dates_cfn$scan1_d_cfn[!is.na(dates_cfn$scan2_d_cfn)],na.rm=T) # scan 2 date always later than scan 1, ok!

#------ dates dataset from Kosha -----------------
# CNB and scan dates from Kosha: n1601_dates_cnb_imaging.csv 
dates_k$cnb1_d_k <- as.Date(dates_k$go1_date_of_cnb,"%m/%d/%y")
dates_k$cnb2_d_k <- as.Date(dates_k$go2_date_of_cnb,"%m/%d/%y")
dates_k$scan1_d_k <- as.Date(dates_k$go1_date_of_scan, "%m/%d/%y")
dates_k$scan2_d_k <- as.Date(dates_k$go2_date_of_scan, "%m/%d/%y") # scan2 is always later than scan 1, good.

# modifying this because it's the only discrepancy between kosha's and cfn path dates
dates_k$scan2_d_k[which(dates_k$bblid==118513)]<-"2012-07-21" # wrong date was 2012-07-28

##############################################
#        MERGING DATA AND FORMATTING         #
##############################################
require(plyr)
itc_all<- join_all(list(itc_o_single, demo, dates_k, dates_cfn, data[,c("bblid","healthExclude")]), by="bblid")

# CNB2 is date of ITC
itc_all$difCnb2Scan2 <- abs(as.numeric(difftime(itc_all$scan2_d_cfn, itc_all$cnb2_d_k, unit="days"))/30)
itc_all$difCnb2Scan1 <- abs(as.numeric(difftime(itc_all$scan1_d_cfn, itc_all$cnb2_d_k, unit="days"))/30)
# pairwise minimium diff <-- this is the last and relevant time difference to look at, 6/15/2016
itc_all$min_diff_with_scan <-pmin(itc_all$difCnb2Scan1,itc_all$difCnb2Scan2,na.rm=T) # time difference from closest scan
scan1_ids <- which(itc_all$min_diff_with_scan == itc_all$difCnb2Scan1) # ids of subjects with scan 1 as closer scan
scan2_ids <- which(itc_all$min_diff_with_scan == itc_all$difCnb2Scan2) # ids of subjects with scan 2 as closer scan
# note that this is not exactly the same as GO1 and GO2 because if Go1 data was bad, scan1 date could be from Go2 scan

itc_all$date_closest_scan[scan1_ids] <- substr(itc_all$dates,1,8)[scan1_ids] # getting exact date of scan for path, file extraction
itc_all$date_closest_scan[scan2_ids] <- substr(itc_all$dates,15,22)[scan2_ids]

itc_all$scanid_closest_scan[scan1_ids] <- as.numeric(substr(itc_all$dates,10,13)[scan1_ids]) # getting exact scanid of closest scan for file extraction
itc_all$scanid_closest_scan[scan2_ids] <- as.numeric(substr(itc_all$dates,24,27)[scan2_ids])

# figuring out which age variable to use
# we're using Kosha's scan1 and scan2 ids; can't use the cfn paths dates, because sometimes scan1 is actually scan 2 if the first scan was bad
from_scan1 <- which(itc_all$scanid_closest_scan==itc_all$go1_scan_id)
from_scan2 <- which(itc_all$scanid_closest_scan==itc_all$go2_scan_id)

length(from_scan1) # 109 scans from Go1
length(from_scan2) # 341 scans from Go2

itc_all$timepoint[from_scan1] <- "go1" # whether scan closest to ITC was GO1 or GO2
itc_all$timepoint[from_scan2] <- "go2"

itc_all$ageAtScan[from_scan1] <- itc_all$ageAtGo1Scan[from_scan1] # age at scan comes from scan1
itc_all$ageAtScan[from_scan2] <- itc_all$ageAtGo2Scan[from_scan2] # age at scan comes from scan2

# ordered sex variable, need for GAM
itc_all$sex.o <- ordered(itc_all$sex)

# drop original scanid variable and replace with "scanid_closest_scan"
comb2<-itc_all[,-c(2,3)]
names(comb2)[which(names(comb2)=="scanid_closest_scan")]<-"scanid"

# add CNB data matched by timepoint 
# there's no scanid in CNB data, but presumably timepoint 1=Go1, timepoint 2=Go2
# ADD CNB data LATER!!!

qa_data<- rbind(qa_go1[,c("bblid","scanid","averageRating")], qa_go2[,c("bblid","scanid","averageRating")])
comb3<-merge(comb2, qa_data, by=c("bblid","scanid"), all.x=T) # will complete with NAs for subjects who don't match by bblid and scanid

##############################################
#               ADDING K-VALUES              #
##############################################

#----------- add k-values from matlab script based on n453_ITC_unique_bblids.csv
ks <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/itcRiskData/n453_Kvalues_06_2016.csv', header=T, sep=',')
comb4 <- merge(comb3,ks[-which(ks$mrnR2<.30),],by="bblid") # discard 5 observations with bad k-values

#--------- save dataset with everyone who has ITC data with processed k values and demo ---------------------
#write.csv(comb4,'/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n453_k_demo_June2016.csv',row.names=FALSE)

##############################################
#	      ADDING RAVENS DATA             #
##############################################

#*******************************************************************************************************************
#                                    HERE RUN TO OBTAIN CNF RAVENS PATHS:
# /data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/sampleAssembly/ravens_get_scan_paths_061616.sh *
#*******************************************************************************************************************

#*******************************************************************************************************************
#                                    HERE RUN TO OBTAIN TOTAL BRAIN VOLUME:
# /data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/ravens/ravens_get_tbv_062216.sh		   *
#*******************************************************************************************************************

# image paths for GAM wrapper
ravensPaths<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/list_of_ravens_image_paths_n449.txt',header=T,sep=',')
# total brain volume
tbv <- read.table('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/tbv_values.txt', header=T, sep=',')

ravensData<-merge(tbv, ravensPaths, by=c("bblid","scanid"), all.x=T)
comb5<-merge(comb4, ravensData, by=c("bblid","scanid"), all.x=T)

# ravens analysis exclusion
# exclude 1. 3 subjects who don't have a scan 
#	  2. 1 subject whose scan is at 23 months of difference
#         3. 18 subjects for health reasons, from healthExclude
#	  4. 1 subject with a T1 QA rating of 0
#     Total: 23

comb5$scanIncludeAll <- 1

# note that this variable should be the minimal exclusion for ANY imaging modality
excl <- unique(c(which(is.na(comb5$scanid)),which(comb5$min_diff_with_scan>23),which(comb5$healthExclude==1),which(comb5$averageRating==0)))

comb5$scanIncludeAll[excl]<-0

##############################################
#                 DTI DATA                   #
##############################################

# FA ROI 
faDti <-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/go1_go2_longitudinal_dti_fa_JHUTracts.csv',header=T,sep=',') 
comb6<-merge(comb5, faDti, by=c("bblid","scanid"), all.x=T)
comb6$faDtiInclude <- 1

# exclude NA fa data
comb6$faDtiInclude[which(is.na(comb6$ATR_L_fa))]<-0 
comb6$faDtiInclude[which(comb5$scanIncludeAll==0)]<-0

##############################################
#            REGIONAL VOLUME DATA            #
##############################################

# MARS data
marsData<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/go1_go2_longitudinal_volume_mars.csv',header=T,sep=',')
comb7<-merge(comb6, marsData[,-3], by=c("bblid","scanid"), all.x=T)
# mars regional data exclusions
comb7$marsRegInclude<-1 
comb7$marsRegInclude[which(is.na(comb7$mprage_mars_vol_L_TTG))]<-0
comb7$marsRegInclude[which(comb5$scanIncludeAll==0)]<-0

# mprageMassICV.y from mars data is THE SAME as the TBV I caclulated = great.

# Muse data 
museData<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/go1_go2_longitudinal_volume_muse.csv',header=T,sep=',')
comb8<-merge(comb7, museData[,-269], by=c("bblid","scanid"), all.x=T)
# muse regional data exclusions
comb8$museRegInclude<-1
comb8$museRegInclude[which(is.na(comb8$vol_L_TTG))]<-0
comb8$museRegInclude[which(comb5$scanIncludeAll==0)]<-0

# NB: ICV from muse data is highly correlated with TBV (r=.99), but NOT the same, several large diffs

##############################################
#                  ANTS CT                   #
##############################################

#**********************************************************************************************************************
#                                    HERE RUN TO OBTAIN ANTS CT IMAGE PATHS:
# /data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/sampleAssembly/antsCT_get_image_paths_070516.sh   *
#**********************************************************************************************************************

antsCtData<-read.table('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/list_of_antsCT_image_paths.txt',header=T,sep=',')
comb9<-merge(comb8, antsCtData, by=c("bblid","scanid"), all.x=T)

# use the scanIncludeAll for exlusions in antsCT analyses

##############################################
#            	 ASL, CBF                    #
##############################################

#*******************************************************************************************************************
#                                    HERE RUN TO OBTAIN CBF IMAGE PATHS:
# /data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/sampleAssembly/cbf_get_image_paths_070516.sh   *
#*******************************************************************************************************************

# as of 7/22/16, CBF data is supposedly not fully processed, missing lots of it

cbfData <- read.table('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/asl/list_of_cbf_image_paths.txt',header=T,sep=',')
comb10<-merge(comb9, cbfData, by=c("bblid","scanid"), all.x=T)

comb10$cbfRegInclude<-1
comb10$cbfRegInclude[which(is.na(comb10$cbfPath))]<-0
comb10$cbfRegInclude[which(comb10$scanIncludeAll==0)]<-0


##############################################
#                SAVING DATA                 #
##############################################
#whichData<-sort(ls(pattern= "comb"),decreasing=T)[1] # get the highest numbered comb data frame, i.e. more recent
whichData <- comb10		 # define data to save
theDate<-gsub("-","",Sys.Date()) # add current date to save with most recent date
theN<-dim(whichData)[1] 	 # number of observations
write.csv(whichData, sprintf('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n%d_pnc_itc_whole_sample_%s.csv',theN,theDate),row.names=F)
saveRDS(whichData, sprintf('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n%d_pnc_itc_whole_sample_%s.rds',theN,theDate))

