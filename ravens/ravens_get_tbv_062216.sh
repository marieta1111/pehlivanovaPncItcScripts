# 6/22/16, script to calculate total brain volume (TBV) from ravens images
# by MP
# example: fslstats /data/joy/BBL/studies/pnc/processedData/structural/ravens/91025/20130629x8181/091025_008181_mprage_brain_bc.nii.gz -V

input_list=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/n449_ITC_ids_for_scans_June2016.csv

missingImages=/data/joy/BBL/projects/pehlivanovaPncItc/logs/missingImagesForTBV.txt
# file that will contain the TBV values and other relevant info
tbvFile=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/tbv_values.txt

rm -f $missingImages
rm -f $tbvFile

# which image is whether the _mprage_brain_bc.nii.gz or _mprage_brain.nii.gz was used
echo "bblid","scanid","whichImage","tbv" >> $tbvFile
for i in $(cat $input_list); do

        echo ""
        echo "--- bblid & scanid ---"

        bblid=$(echo $i|cut -d, -f1)
        scanid=$(echo $i|cut -d, -f3)
        echo $bblid , $scanid

	# define image path for the most common path template
        image1=$(ls /data/joy/BBL/studies/pnc/processedData/structural/ravens/${bblid}/*${scanid}/*${bblid}_*${scanid}_mprage_brain_bc.nii.gz)

	if [ -e "$image1" ]; then
		
		# if the image exits, calculate TBV and input into txt file
		tbvImage=$(ls /data/joy/BBL/studies/pnc/processedData/structural/ravens/${bblid}/*${scanid}/*${bblid}_*${scanid}_mprage_brain_bc.nii.gz)
		tbv=$(fslstats $tbvImage -V |cut -d" " -f2)
		echo "$bblid, $scanid, brain_bc, $tbv" >> $tbvFile
	
	elif [ ! -e "$image1" ]; then 
 		# if image1 does not exist, check if subject has other data format instead
		image2=$(ls /data/joy/BBL/studies/pnc/processedData/structural/ravens/${bblid}/*${scanid}/*${bblid}_*${scanid}_mprage_brain.nii.gz)
  
                if [ -e "$image2" ]; then
			# NOTE: if I put this at the end after devining one "final" image (either one), it adds some weird output ot the text file...	
			tbvImage=$(ls /data/joy/BBL/studies/pnc/processedData/structural/ravens/${bblid}/*${scanid}/*${bblid}_*${scanid}_mprage_brain.nii.gz)
                	tbv=$(fslstats $tbvImage -V |cut -d" " -f2)
	                echo "$bblid, $scanid, brain , $tbv" >> $tbvFile
		
		# if none of the images are there, log the subject data as missing
		# NOTE: do not understand why this file gets created even though all 449 subjects have a good image
		elif [ !  -e "$image2" ]; then
			echo "ravens bc image is missing"
			echo "$bblid","$scanid" >> $missingImages
		fi

	fi

done































	
