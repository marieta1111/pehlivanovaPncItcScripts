# script to create path variable for cbf images to run GAM models
# start out with a list of all relevant bblids and scanids

input_list=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/n449_ITC_ids_for_scans_June2016.csv

logdir=/data/joy/BBL/projects/pehlivanovaPncItc/
image_log=$logdir/subjectData/asl/list_of_asl_image_paths.txt
missing_log=$logdir/logs/list_of_missing_asl_paths.txt

rm -f $image_log
rm -f $missing_log

echo "bblid","scanid","aslPath" >> $image_log

for i in $(cat $input_list); do
        echo ""
        echo "--- bblid & scanid ---"

        bblid=$(echo $i|cut -d, -f1)
        scanid=$(echo $i|cut -d, -f3)
	echo $bblid , $scanid

	image_file=$(ls /data/joy/BBL/studies/pnc/processedData/pcasl/pcasl_201607291423/${bblid}/*${scanid}/norm/${bblid}_*${scanid}_asl_quant_ssT1Std.nii.gz 2> /dev/null)

	#/data/joy/BBL/studies/pnc/processedData/pcasl/pcasl_201606231423/100031/20100918x3818/norm/100031_20100918x3818_asl_quant_ssT1Std.nii.gz

	echo $image_file
	
	# if image is missing 	
	if [ ! -e "$image_file" ]; then
	        echo "image not present"
                echo "$bblid, $scanid" >> $missing_log

	# output image if present
	elif [ -e "$image_file" ]; then
	 	echo "$bblid, $scanid, $image_file" >> $image_log
	fi
	
done
