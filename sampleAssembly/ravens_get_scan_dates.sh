input_list=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n453_ITC_unique_bblids.csv

logdir=/data/joy/BBL/projects/pehlivanovaPncItc
image_log=$logdir/subjectData/ravens/list_of_ravens_scan_dates.txt
missing_log=$logdir/logs/list_of_missing_ravens_images.txt

rm -f $image_log
rm -f $missing_log

echo "bblid","dates" >> $image_log
for i in $(cat $input_list); do
        echo ""
        echo "--- bblid ---"

        bblid=$(echo $i|cut -d, -f1)
	
	echo $bblid 
	
	# list folders for specific bblid

	folders=$(ls /data/joy/BBL/studies/pnc/processedData/structural/ravens/${bblid} 2> /dev/null)
	
	# print said folders	
	echo $folders
	folder1=$(echo $folders|cut -d, -f1)
	
	# if folders exist, add to text file next to bblid
	if [ ! -z "$folders" ]; then
		
		echo "$bblid,$folder1" >> $image_log
	# if folders var is empty, i.e. no folder, log in the missing scan file
	elif [ -z "$folders" ]; then
		echo "****** NO RAVENS DATA *****"
		echo "$bblid" >> $missing_log
	fi
done
