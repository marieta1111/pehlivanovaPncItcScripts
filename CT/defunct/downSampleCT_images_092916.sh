# script to downsample CT images

input_list=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages/listOf427CtImages.txt

inDir=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages
outDir=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages_2mm

for i in $(cat $input_list); do
        echo ""

       	image=$(echo $i|cut -d, -f1)
	imageName=$(basename $image | cut -d. -f1)
	echo $image

	fslmaths $inDir/$image -subsamp2 $outDir/${imageName}_2mm	
	
done
