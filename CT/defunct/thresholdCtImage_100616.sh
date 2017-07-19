input_list=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages/listOf427CtImages.txt

inDir=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages
outDir=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages_thr_ants

for i in $(cat $input_list); do
        echo ""

        image=$(echo $i|cut -d, -f1)
        imageName=$(basename $image | cut -d. -f1)
        echo $image	
	ThresholdImage 3 $inDir/$image $outDir/${imageName}_01mask.nii.gz 0.1 Inf

done

AverageImages 3 ctThicknessCoverage.nii.gz 0  $outDir/*mask.nii.gz
