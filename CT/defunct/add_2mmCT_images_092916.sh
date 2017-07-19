# script to threshold CT images

input_list=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages_2mm/listOf2mmCtImages.txt

outDir=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages_2mm

cp $outDir/100050_8309_CorticalThicknessNormalizedToTemplate_2mm.nii.gz $outDir/sum2mmCtImagesN427.nii.gz

for i in $(cat $input_list | tail -426); do
        echo ""
	
       	image=$(echo $i|cut -d, -f1)
	echo $image 

	fslmaths $outDir/sum2mmCtImagesN427.nii.gz -add $outDir/$image $outDir/sum2mmCtImagesN427.nii.gz	
	
done
