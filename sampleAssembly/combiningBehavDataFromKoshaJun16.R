#----------------------------------------------------------------------------------------------------------------------------
# JUNE16: this script combines and formats various datasets relevant to the
#
# PNC-ITC project; most individual datasets were compiled by Kosha 
#
# It produces several datasets for further analysis or processing:
#
# 1. dataset with n=453 subjects who have full ITC data to 1) pull out scan dates and 2) to caclulate k-values from ITC data
# 	/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n453_ITC_unique_bblids.csv
#	
#	NB: Currently commented out
#
# 2. dataset with n=427 subjects with processed ITC (*only* usable k-values), paths for imaging analysis, and imaging inclusion variables
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
# ITC: this has n=655 and not 656 observations as the file name states
itc <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n656_scales_mfq_asrm_panas_spq_ddisc_rdisc_wolfq.csv',head=T,sep=',')
# CNB measures
cnb <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/GO1_GO2_CNB_Factor_Scores_Non-Age-Regressed.csv',head=T,sep=',')
# data release
data <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n1601_go1_datarel_020716.csv',head=T,sep=',')
# imaging QA measures
qa_go1 <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n1601_t1RawManualQA.csv',header=T,sep=',')
qa_go2 <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n368_t1RawManualQA_GO2.csv',header=T,sep=',')
# substance use
subs <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/pnc_substance_cnb_go1_go2_10215.csv', head=T,sep=',')

########################################################
#             FORMATTING ITC AND RISK  DATA            #
#######################################################

# order itc data by bblid and timepoint
itc_risk <- itc[with(itc, order(bblid, timepoint)), ]
# 			ITC
# just sample with observations with full ITC data
itc_full_o <- itc_risk[which(!is.na(itc_risk$kddisc_q_01)),c(1:2,4,207:275)]
# keep just one record per person, FIRST observation/scan
# 32 subjects have duplicate ITC data
itc_o_single <- itc_full_o[!duplicated(itc_full_o$bblid),]
# save data to pull out scan dates
#write.csv(itc_o_single,'/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n453_ITC_unique_bblids.csv',row.names=FALSE)

# 			Risk
risk_full_o <- itc_risk[which(!is.na(itc_risk$krdisc_q_01)),c(1,4,276:358)]
# 32 subjects have duplicate Risk data, remove those
# bblid 102845 has risk data, but not ITC data
risk_o_single <- risk_full_o[!duplicated(risk_full_o$bblid),]
# save data to run risk model
#write.csv(risk_o_single,'/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n454_Risk_data_unique_bblids.csv',row.names=FALSE)

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
itc_all <- join_all(list(risk_o_single[,-2],itc_o_single[,-c(2,3)], demo, dates_k, dates_cfn, data[,c("bblid","healthExclude")]), by="bblid")

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
# we're using Kosha's scan1 and scan2 ids; can't use the CFN paths dates, because sometimes scan1 is actually scan2 if the first scan was bad
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
comb2<-itc_all
names(comb2)[which(names(comb2)=="scanid_closest_scan")]<-"scanid"

# add CNB data matched by timepoint 
# there's no scanid in CNB data, but presumably timepoint 1=Go1, timepoint 2=Go2

qa_data<- rbind(qa_go1[,c("bblid","scanid","averageRating")], qa_go2[,c("bblid","scanid","averageRating")])
comb3<-merge(comb2, qa_data, by=c("bblid","scanid"), all.x=T) # will complete with NAs for subjects who don't match by bblid and scanid

# 08/15/16: save data for extraction of scan paths
#write.csv(comb3[,c("bblid","date_closest_scan","scanid","min_diff_with_scan")],'/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/n454_ITC_Risk_ids_for_scans_Aug2016.csv',row.names=FALSE)

##############################################
#         ADDING K- AND ALPHA-VALUES         #
##############################################

#----------- add k-values from matlab script based on n453_ITC_unique_bblids.csv
#----------- and alpha-values from from matlab script based on n454_Risk_data_unique_bblids.csv

ks <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/itcRiskData/n453_Kvalues_08162016.csv', header=T, sep=',')
alphas <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/itcRiskData/n454_Alphas_08162016.csv', header=T, sep=',')
# discard observations with k- and alpha-values that don't pass QC
comb_k_alpha <- merge(ks[-which(ks$kTjur<.20),],alphas[-which(alphas$alphaTjur<.20),], by="bblid",all=T) 

comb4 <- merge(comb3, comb_k_alpha, by="bblid") 

#--------- save dataset with everyone who has ITC data with processed k values and demo ---------------------
#write.csv(comb4,'/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n453_k_demo_June2016.csv',row.names=FALSE)

##############################################
#	      ADDING RAVENS DATA             #
##############################################

#************************************************************************************************************************
#                                    HERE RUN TO OBTAIN CNF RAVENS PATHS:
# /data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/sampleAssembly/ravens_get_scan_paths_081516.sh
# old: /data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/sampleAssembly/ravens_get_scan_paths_061616.sh *
#************************************************************************************************************************

#**************************************************************************************************************************
#                                    HERE RUN TO OBTAIN TOTAL BRAIN VOLUME:
# /data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/ravens/ravens_get_tbv_081516.sh
# old: /data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/ravens/ravens_get_tbv_062216.sh	          *
#**************************************************************************************************************************

# image paths for GAM wrapper
ravensPaths<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/list_of_ravens_image_paths_n454_Aug16.txt',header=T,sep=',')
# total brain volume
tbv <- read.table('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/tbv_values_Aug16.txt', header=T, sep=',')

ravensData<-merge(tbv, ravensPaths, by=c("bblid","scanid"), all.x=T)
comb5<-merge(comb4, ravensData, by=c("bblid","scanid"), all.x=T)

# ravens analysis exclusion
# exclude 1. 4 subjects who don't have a scan 
#	  2. 1 subject whose scan is at 23 months of difference
#         3. 19 subjects for health reasons, from healthExclude
#	  4. 1 subject with a T1 QA rating of 0
#     Total: 25
#     Plus exlcusions for missing behavioral data

# JUST FOR LOGK out of 451 subjects with good ITC data: 
#	 3 subjects who don't have a scan;  
#        1 subject whose scan is at 23 months of difference
#        19 subjects for health reasons, from healthExclude
#        1 subject with a T1 QA rating of 0


# general inclusion variable
comb5$scanIncludeAll <- 1 
# inclusion for delay discounting analyses
comb5$scanIncludeAll_k <- 1
# inclusion for risk analyses 
comb5$scanIncludeAll_alpha <- 1

# note that this variable should be the minimal exclusion for ANY imaging modality
excl <- unique(c(which(is.na(comb5$scanid)),which(comb5$min_diff_with_scan>23),which(comb5$healthExclude==1),which(comb5$averageRating==0)))
k_nas <- which(is.na(comb5$k))
a_nas <- which(is.na(comb5$alpha)) 
excl_k <- unique(c(excl,k_nas))
excl_a <- unique(c(excl,a_nas))

comb5$scanIncludeAll[excl]<-0
comb5$scanIncludeAll_k[excl_k] <-0
comb5$scanIncludeAll_alpha[excl_a] <-0

##############################################
#                 DTI DATA                   #
##############################################

# FA ROI 
faDti <-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/go1_go2_longitudinal_dti_fa_JHUTracts.csv',header=T,sep=',') 
comb6<-merge(comb5, faDti[,-2], by=c("bblid","scanid"), all.x=T)
comb6$faDtiInclude <- 1
comb6$faDtiInclude_k <- 1
comb6$faDtiInclude_alpha <- 1

# exclusions:
e<-unique(c(which(is.na(comb6$ATR_L_fa)), which(comb5$scanIncludeAll==0))) #NA fa data and general scan and health exclusions
e_k <- unique(c(e,k_nas))
e_a <- unique(c(e,a_nas))

comb6$faDtiInclude[e]<-0 
comb6$faDtiInclude_k[e_k]<-0
comb6$faDtiInclude_alpha[e_a]<-0

##############################################
#            REGIONAL VOLUME DATA            #
##############################################

# MARS data
marsData<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/go1_go2_longitudinal_volume_mars.csv',header=T,sep=',')
comb7<-merge(comb6, marsData[,-3], by=c("bblid","scanid"), all.x=T)
# mars regional data exclusions
comb7$marsRegInclude<-1 
comb7$marsRegInclude_k<-1
comb7$marsRegInclude_alpha<-1

# exclusions
e_mars <- unique(c(which(is.na(comb7$mprage_mars_vol_L_TTG)),which(comb7$scanIncludeAll==0)))
e_mars_k <- unique(c(e_mars,k_nas))
e_mars_a <- unique(c(e_mars,a_nas))

comb7$marsRegInclude[e_mars]<-0
comb7$marsRegInclude_k[e_mars_k]<-0
comb7$marsRegInclude_alpha[e_mars_a]<-0

# mprageMassICV.y from mars data is THE SAME as the TBV I caclulated = great.

# Muse data 
museData<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/go1_go2_longitudinal_volume_muse.csv',header=T,sep=',')
comb8<-merge(comb7, museData[,-c(269,270)], by=c("bblid","scanid"), all.x=T)
# muse regional data exclusions
comb8$museRegInclude<-1
comb8$museRegInclude_k<-1
comb8$museRegInclude_alpha<-1

# exclusions
e_muse <- unique(c(which(is.na(comb8$vol_L_TTG)),which(comb8$scanIncludeAll==0)))
e_muse_k <- unique(c(e_muse,k_nas))
e_muse_a <- unique(c(e_muse,a_nas))

comb8$museRegInclude[e_muse]<-0
comb8$museRegInclude_k[e_muse_k]<-0
comb8$museRegInclude_alpha[e_muse_a]<-0

# NB: ICV from muse data is highly correlated with TBV (r=.99), but NOT the same, several large diffs

##############################################
#                  ANTS CT                   #
##############################################

antsCtData<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n1601_Corticalthickness_Longitudinal.csv',header=T,sep=',')

# names are reading in with dots and ellipses (were blanks in original csv); want to fix that
# removing dots and add CT in front of CT variables 
antsCtData_goodNames<-antsCtData
names(antsCtData_goodNames)<-c("X","bblid","scanid",paste0("CT_",gsub("[.]","",names(antsCtData)[4:151])))

comb9<-merge(comb8, antsCtData_goodNames[,-1], by=c("bblid","scanid"), all.x=T)

# use the scanIncludeAll for exlusions in antsCT analyses

##############################################
#            	 ASL	                     #
##############################################

#************************************************************************************************************************
#                                    HERE RUN TO OBTAIN CBF IMAGE PATHS:
#      /data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/sampleAssembly/asl_get_image_paths_082216.sh   *
# old: /data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/sampleAssembly/asl_get_image_paths_070516.sh   *
#************************************************************************************************************************

cbfData <- read.table('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/asl/list_of_asl_image_paths_Aug16.txt',header=T,sep=',')

# ASL QA from Adon
qa_asl <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n1895_FlagStatus.csv',header=T,sep=',')
names(qa_asl)[1:2]=c("bblid","dateXscanid")
qa_asl$scanid <- as.numeric(substr(qa_asl$dateXscanid,10,13))
# combine all ASL data
asl_data <- merge(cbfData,qa_asl[,-2],by=c("bblid","scanid"), all.x=T)

# combine with other imaging and behavioral dataa 
comb10<-merge(comb9, asl_data, by=c("bblid","scanid"), all.x=T)

# ASL exlcusions, missing data, scan + health exclusions + QA image exclusion
comb10$aslRegInclude<-1
comb10$aslRegInclude[which(is.na(comb10$aslPath))]<-0
comb10$aslRegInclude[which(comb10$scanIncludeAll==0)]<-0
comb10$aslRegInclude[which(comb10$pcaslExclude==1)]<-0

# k analyses exclusions
comb10$aslRegInclude_k <-comb10$aslRegInclude
comb10$aslRegInclude_k[k_nas]<-0

# alpha analyses exclusions
comb10$aslRegInclude_alpha <-comb10$aslRegInclude
comb10$aslRegInclude_alpha[a_nas]<-0

##############################################
#              JLF Regional Data             #
##############################################

#------------- JLF Volume -------------------#
#jlfVol<-read.csv('/data/joy/BBL/projects/pncReproc2015/jlf/volumeValues/n1960_jlfVol_antsCTVol_T1QA_mm3.csv',header=T,sep=',')
jlfVol<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/allJlfVolValuesmm3.csv', header=T,sep=',')
comb11<-merge(comb10, jlfVol[,-3], by=c("bblid","scanid"), all.x=T)

# exclusions (there's no missing JLF Vol data)
comb11$jlfVolInclude<-1
comb11$jlfVolInclude[which(comb11$scanIncludeAll==0)]<-0
# k analyses exclusions
comb11$jlfVolInclude_k <- comb11$jlfVolInclude
comb11$jlfVolInclude_k[k_nas]<-0
# alpha analyses exclusions
comb11$jlfVolInclude_alpha <- comb11$jlfVolInclude
comb11$jlfVolInclude_alpha[a_nas]<-0

#------------- JLF CT ------------------------#
jlfCT <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n2187_jlfCT.csv', header=T, sep=',')
comb12<- merge(comb11, jlfCT[,-3], by=c("bblid","scanid"), all.x=T)

# exclusions (there's no missing JLF CT data)
comb12$jlfCTInclude<-1
comb12$jlfCTInclude[which(comb12$scanIncludeAll==0)]<-0
# k analyses exclusions
comb12$jlfCTInclude_k <- comb12$jlfCTInclude
comb12$jlfCTInclude_k[k_nas]<-0
# alpha analyses exclusions
comb12$jlfCTInclude_alpha <- comb12$jlfCTInclude
comb12$jlfCTInclude_alpha[a_nas]<-0

#------------- JLF ASL ------------------------#
jlfSST1<- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/pcasl_JLF_ssT1-correctHeaders.csv', header=T, sep=',')
comb13<- merge(comb12, jlfSST1, by=c("bblid","scanid"), all.x=T)

# exclusions (there's no missing JLF CT data)
comb13$jlfCbfInclude<-1
comb13$jlfCbfInclude[which(comb13$scanIncludeAll==0)]<-0
comb13$jlfCbfInclude[which(is.na(comb13$pcasl_jlf_cbf_L_Cun))]<-0
# k analyses exclusions
comb13$jlfCbfInclude_k <- comb13$jlfCbfInclude
comb13$jlfCbfInclude_k[k_nas]<-0
# alpha analyses exclusions
comb13$jlfCbfInclude_alpha <- comb13$jlfCbfInclude
comb13$jlfCbfInclude_alpha[a_nas]<-0

##############################################
#                   CNB DATA                 #
##############################################

# just keep timepoint 2, CNB 2 because ITC data was collected as part of CNB2
cnb_t2<-cnb[which(cnb$timepoint==2),]
comb14 <- merge(comb13, cnb_t2[,-2], by=c("bblid"), all.x=T)

# take data from bigger dataset
# NOTE, THIS IS ONLY for Go1!!! does not seem appropriate for time sensitive data 
data_sub<-data[,c(1:2,7,2275:2335,2337:2475)]
comb15 <-merge(comb14, data_sub[,-2], by=c("bblid"), all.x=T)

##############################################
#             Voxel-wise CT maps             #
##############################################

#**************************************************************************************************************************
#                                    HERE RUN TO OBTAIN CT IMAGE PATHS:							  *
#      /data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/sampleAssembly/antsCT_get_image_paths_092516.sh  *
#      then ran command to get a list of all the files ls *ToTemplate.nii.gz >> listOfCtImages.txt	 		  *
#**************************************************************************************************************************

ctPaths <- read.table('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages/listOfCtImages.txt',header=F)
names(ctPaths)<-"ctImagePath"
ctPaths1 <- ctPaths
ctPaths1$ctImagePath <- paste0("/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages/",ctPaths$ctImagePath)
# extract bblid to match with other data
ctPaths1$bblid <- as.numeric(sapply(strsplit(as.character(ctPaths[,1]),"_"), "[", 1))
ctPaths1$scanid <- as.numeric(sapply(strsplit(as.character(ctPaths[,1]),"_"), "[", 2))

comb16 <-merge(comb15, ctPaths1, by=c("bblid","scanid"), all.x=T)

# can use these inclusion variables for voxel-wise CT data: scanIncludeAll, scanIncludeAll_k, scanIncludeAll_alpha

##############################################
#             NMF COMPONENT RESULTS          #
##############################################

nmf20<-read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/NMF/NmfResults20Bases.csv',header=T,sep=',')
orderBblid<-read.table('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages_2mm/listOf427Ct2mmImages.txt',header=F,sep='')
bblid<-as.numeric(sapply(strsplit(as.character(orderBblid[,1]),"_"),"[", 1))
nmf20ord<-cbind(bblid,nmf20)

comb17 <-merge(comb16, nmf20ord, by=c("bblid"), all.x=T)

##############################################
#             SUBSTANCE USE DATA             #
##############################################

# data from Dan Wolf
dw1 <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/itc_substance_dw4marieta_102716_sheet1.csv', head=T,sep=',')
dw2 <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/itc_substance_dw4marieta_102716_sheet2.csv', head=T,sep=',')

# from general substance data
subs$timepoint1 <- subs$timepoint
subs$timepoint[which(subs$timepoint==1)]<-"go1"
subs$timepoint[which(subs$timepoint==2)]<-"go2"

# these bblids of interest have duplicate go2 records: 81865 83253 83835 85369 92896
# will remove the empty ones, or the ones will less data
subs <- subs[-which(subs$cnb_datasetid %in% c(49516, 60870, 49117, 49183, 50867)),]

# calculate BIS BAS variables
itc$BAS_dr <- itc$bisbas_c_q_02 + itc$bisbas_c_q_07 + itc$bisbas_c_q_09 + itc$bisbas_c_q_17
itc$BAS_fs <- itc$bisbas_c_q_04 + itc$bisbas_c_q_08 + itc$bisbas_c_q_12 + itc$bisbas_c_q_16
itc$BAS_rr <- itc$bisbas_c_q_03 + itc$bisbas_c_q_05 + itc$bisbas_c_q_11 + itc$bisbas_c_q_14 + itc$bisbas_c_q_19
itc$BIS_   <- (3-itc$bisbas_c_q_01) + itc$bisbas_c_q_06 + itc$bisbas_c_q_10 + itc$bisbas_c_q_13 + itc$bisbas_c_q_15 + itc$bisbas_c_q_18 + itc$bisbas_c_q_20

comb18<-merge(comb17, dw1[,-c(2:9,15:17,23,95:128,216:230)], by=c("bblid"), all.x=T) # remove some variables that repeat with my combined dataset 
comb19<-merge(comb18, dw2, by=c("bblid"), all.x=T)
comb20<-merge(comb19, itc[!duplicated(itc$bblid),c(1,359:362)], by=c("bblid"), all.x=T)
comb21<-merge(comb20, subs, by=c("bblid","timepoint"),all.x=T)

# testing dan's substance use variables
# flip direction of freq variables so now more means more frequent use
# set NAs to 0, no use
comb21$substance_tob_freq_rev <- 7-comb21$cnb_substance_tob_060
comb21$substance_tob_freq_rev[is.na(comb21$substance_tob_freq_rev)]<-0
comb21$substance_alc_freq_rev <- 7-comb21$cnb_substance_alc_050
comb21$substance_alc_freq_rev[is.na(comb21$substance_alc_freq_rev)]<-0
comb21$substance_mar_freq_rev <- 7-comb21$cnb_substance_mar_040
comb21$substance_mar_freq_rev[is.na(comb21$substance_mar_freq_rev)]<-0
comb21$substance_all_avg_freq<- (comb21$substance_tob_freq_rev + comb21$substance_alc_freq_rev + comb21$substance_mar_freq_rev)/3

##############################################
#                SAVING DATA                 #
##############################################

#whichData<-sort(ls(pattern= "comb"),decreasing=T)[1] # get the highest numbered comb data frame, i.e. more recent
whichData <- comb21		 # which data to save
theDate<-gsub("-","",Sys.Date()) # add current date to save with most recent date
theN<-dim(whichData)[1] 	 # number of observations
write.csv(whichData, sprintf('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n%d_pnc_itc_whole_sample_%s.csv',theN,theDate),row.names=F)
saveRDS(whichData, sprintf('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n%d_pnc_itc_whole_sample_%s.rds',theN,theDate))

