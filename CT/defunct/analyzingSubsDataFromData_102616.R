# read in data
itc <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n656_scales_mfq_asrm_panas_spq_ddisc_rdisc_wolfq.csv',head=T,sep=',')

# substance data from dan wolf
dw1 <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/itc_substance_dw4marieta_102716_sheet1.csv', head=T,sep=',')
dw2 <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/itc_substance_dw4marieta_102716_sheet2.csv', head=T,sep=',')

# all substance data
subs <- read.csv('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/pnc_substance_cnb_go1_go2_10215.csv', head=T,sep=',')
subs$timepoint1 <- subs$timepoint
subs$timepoint[which(subs$timepoint==1)]<-"go1"
subs$timepoint[which(subs$timepoint==2)]<-"go2"
# these bblids of interest have duplicate go2 records: 81865 83253 83835 85369 92896
# will remove the empty ones, or the ones will less data
subs <- subs[-which(subs$cnb_datasetid %in% c(49516, 60870, 49117, 49183, 50867)),]

# all my data
data<-readRDS('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n452_pnc_itc_whole_sample_20160905.rds')

# calculate BIS BAS variables
itc$BAS_dr <- itc$bisbas_c_q_02 + itc$bisbas_c_q_07 + itc$bisbas_c_q_09 + itc$bisbas_c_q_17
itc$BAS_fs <- itc$bisbas_c_q_04 + itc$bisbas_c_q_08 + itc$bisbas_c_q_12 + itc$bisbas_c_q_16
itc$BAS_rr <- itc$bisbas_c_q_03 + itc$bisbas_c_q_05 + itc$bisbas_c_q_11 + itc$bisbas_c_q_14 + itc$bisbas_c_q_19
itc$BIS_   <- (3-itc$bisbas_c_q_01) + itc$bisbas_c_q_06 + itc$bisbas_c_q_10 + itc$bisbas_c_q_13 + itc$bisbas_c_q_15 + itc$bisbas_c_q_18 + itc$bisbas_c_q_20

# combine datasets
comb<-merge(data, dw1[,-c(2:9,15:17,95:128,216:230)], by=c("bblid"), all.x=T) # remove some variables that repeat with my combined dataset 
comb2<-merge(comb, dw2, by=c("bblid"), all.x=T)
comb3<-merge(comb2, itc[!duplicated(itc$bblid),c(1,359:362)], by=c("bblid"), all.x=T)
comb4<-merge(comb3,subs,by=c("bblid","timepoint"),all.x=T)

# testing dan's substance use variables
# flip direction of freq variables so now more means more frequent use
# set NAs to 0, no use
comb4$substance_tob_freq_rev <- 7-comb4$cnb_substance_tob_060
comb4$substance_tob_freq_rev[is.na(comb4$substance_tob_freq_rev)]<-0
comb4$substance_alc_freq_rev <- 7-comb4$cnb_substance_alc_050
comb4$substance_alc_freq_rev[is.na(comb4$substance_alc_freq_rev)]<-0
comb4$substance_mar_freq_rev <- 7-comb4$cnb_substance_mar_040
comb4$substance_mar_freq_rev[is.na(comb4$substance_mar_freq_rev)]<-0

comb4$substance_all_avg_freq<- (comb4$substance_tob_freq_rev + comb4$substance_alc_freq_rev + comb4$substance_mar_freq_rev)/3

# all 0s in the summary measures are also 0s in the edited frequency variables
sum(comb4$substance_tob_freq_rev[which(comb4$subs_smry_sub_tob==0)])
sum(comb4$substance_alc_freq_rev[which(comb4$subs_smry_sub_alc==0)])
sum(comb4$substance_mar_freq_rev[which(comb4$subs_smry_sub_mar==0)])

ss<-which(comb4$scanIncludeAll==1) # subjects who have good scans, n=427
a15<-which(comb4$ageAtScan>=180) # subjects of age 15+
a16<-which(comb4$ageAtScan>=192) # subjects of age 16+
a17<-which(comb4$ageAtScan>=204) # subjects of age 17+

cor.test(comb4$substance_tob_freq_rev,comb4$substance_alc_freq_rev,method="spearman")

# relationship of substance variables with logk
cor.test(comb4$substance_tob_freq_rev[],comb4$logk[],method="spearman") # sig with pearson, but very small)
cor.test(comb4$substance_alc_freq_rev[],comb4$logk[],method="spearman")
cor.test(comb4$substance_mar_freq_rev[],comb4$logk[],method="spearman") # sig with pearson, but very small
cor.test(comb4$substance_all_avg_freq[],comb4$logk[],method="spearman") # sig with pearson, but very small 
cor.test(comb4$subs_smry_sub_tot[],comb4$logk[],method="spearman")

# personality and logk relationship
cor.test(comb4$BIS_[],comb4$logk[],method="spearman")
cor.test(comb4$BAS_dr[],comb4$logk[],method="spearman") # 0.176
cor.test(comb4$BAS_fs[],comb4$logk[],method="spearman")
cor.test(comb4$BAS_rr[],comb4$logk[],method="spearman") # 0.10

# BMI
cor.test(comb4$bmi[],comb4$logk[],method="spearman")

#2 find BMI data matched by timepoint
#3 conduct formal mediation for any variable that's significantly associated with logk including medu

# "bmi"                                       
# "bmiPercentile"                             
# "bmiAge"                                    
# "bmiHeight"                                 
# "bmiWeight" 

