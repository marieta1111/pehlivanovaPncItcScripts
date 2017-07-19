# script to threshold CT images

input_list=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages_01thr/listOfThreshImages.txt

outDir=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages_01thr

cp $outDir/thr01_100050_8309_CorticalThicknessNormalizedToTemplate.nii.gz $outDir/sumGoodVoxN427.nii.gz

for i in $(cat $input_list | tail -426); do
        echo ""
	
       	image=$(echo $i|cut -d, -f1)
	echo $image 

	fslmaths $outDir/sumGoodVoxN427.nii.gz -add $outDir/$image $outDir/sumGoodVoxN427.nii.gz	
	
done
