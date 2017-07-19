# script to create path variable for RAVENS images to run GAM models
# already have list of all relevant bblids and scanids
# nothing should be missing because list of ids was created based on which scans are available 

input_list=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/n454_ITC_Risk_ids_for_scans_Aug2016.csv

logdir=/data/joy/BBL/projects/pehlivanovaPncItc/
image_log=$logdir/subjectData/ravens/list_of_ravens_image_paths_n454_Aug16.txt
missing_log=$logdir/logs/list_of_missing_ravens_paths_Aug16.txt

rm -f $image_log
rm -f $missing_log

echo "bblid","scanid","ravensPath" >> $image_log

for i in $(cat $input_list); do
        echo ""
        echo "--- bblid & scanid ---"

        bblid=$(echo $i|cut -d, -f1)
        scanid=$(echo $i|cut -d, -f3)
	echo $bblid , $scanid
	
	image_file=$(ls /data/joy/BBL/studies/pnc/processedData/structural/ravens/${bblid}/*${scanid}/*${bblid}*${scanid}*_RAVENS_150_2mm.nii.gz 2> /dev/null)
	
	# if downsampled image is missing 	
	if [ ! -e "$image_file" ]; then
		echo "downsampled image not present"
		echo "$image_file" >> $missing_log

		#look for image that is not downsampled
		highRes=$(ls /data/joy/BBL/studies/pnc/processedData/structural/ravens/${bblid}/*${scanid}/*${bblid}*${scanid}*_RAVENS_150.nii.gz 2> /dev/null)
		
		# high resolution image missing as well
		if [ ! -e "$highRes" ]; then
			echo "******HIGH RES IMAGE NOT PRESENT!!! SKIPPING THIS SUBJECT AND LOGGING!!!!*****"
			echo "$bblid,$scanid" >> $missing_log
			continue
		fi
	
		imageName=$(echo $highRes | cut -d. -f1)
		echo "downsampling $imageName"

		# how does the path work here??
		# in any case, all images are there

		fslmaths $highRes -subsamp2 ${imageName}_2mm

	fi
	
	 image_file_final=$(ls /data/joy/BBL/studies/pnc/processedData/structural/ravens/${bblid}/*${scanid}/*${bblid}*${scanid}*_RAVENS_150_2mm.nii.gz) # 2> /dev/null)
	
	if [ ! -e "$image_file_final" ]; then
		echo "*****IMAGE NOT CREATED AS EXPECTED???*******"
	        echo "$bblid,$scanid" >> $missing_log
	fi
	echo "$bblid, $scanid, $image_file_final" >> $image_log

done
