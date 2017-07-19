# script to create path variable for antsCT images to run GAM models
# start out with a list of all relevant bblids and scanids

input_list=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/n454_ITC_Risk_ids_for_scans_Aug2016.csv

dir1=/data/joy/BBL/projects/pehlivanovaPncItc/
image_log=$dir1/subjectData/antsCT/list_of_antsCT_image_paths_092516.txt
missing_log=$dir1/logs/list_of_missing_CT_paths_092516.txt
copy_dir=$dir1/subjectData/antsCT/ctImages

rm -f $image_log
rm -f $missing_log

echo "bblid","scanid","antsCtPath" >> $image_log

for i in $(cat $input_list); do
        echo ""
        echo "--- bblid & scanid ---"

        bblid=$(echo $i|cut -d, -f1)
        scanid=$(echo $i|cut -d, -f3)
	echo $bblid , $scanid

	image_file=$(ls /data/joy/BBL/studies/pnc/processedData/structural/antsCorticalThickness/${bblid}/*${scanid}/CorticalThicknessNormalizedToTemplate.nii.gz 2> /dev/null)

	echo $image_file
	
	# if image is missing 	
	if [ ! -e "$image_file" ]; then

     		echo "image not present"
                echo "$bblid, $scanid" >> $missing_log

	# output image if present
	elif [ -e "$image_file" ]; then
	 	echo "$bblid, $scanid, $image_file" >> $image_log
		cp $image_file $copy_dir/${bblid}_${scanid}_CorticalThicknessNormalizedToTemplate.nii.gz
	fi
	
done
