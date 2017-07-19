# script to threshold CT images

input_list=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages_2mm/listOf2mmCtImages.txt

inDir=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages_2mm
outDir=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages_2mm_s1mm

for i in $(cat $input_list); do
        echo ""

       	image=$(echo $i|cut -d, -f1)
	imageName=$(basename $image | cut -d. -f1)
	echo $image , $imageName

	fslmaths $inDir/$image -s 0.4246609 $outDir/${imageName}_s1mm.nii.gz	
	
done
