# script to threshold CT images

input_list=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages/listOf427CtImages.txt

inDir=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages
outDir=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages_01thr

for i in $(cat $input_list); do
        echo ""

       	image=$(echo $i|cut -d, -f1)
	echo $image 

	fslmaths $inDir/$image -thr 0.1 -bin $outDir/thr01_${image}	
	
done
