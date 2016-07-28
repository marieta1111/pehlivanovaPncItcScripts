input_list=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/n425_ravens_for_midas_072016.csv

mask=/data/joy/BBL/projects/pehlivanovaPncItc/RAVENS/mniGreyCombEons2mmEro.nii.gz

outdir=/data/joy/BBL/projects/pehlivanovaPncItc/RAVENS/smoothedMaskedImages/

for i in $(cat $input_list); do

        echo "--- image path ---"

        imagePath=$(echo $i|cut -d, -f1)
	imageName=$(basename $i | cut -d. -f1)
        echo  $imagePath, ${imageName}

	fslmaths $imagePath -s 0.85 -mas $mask $outdir/${imageName}_s2mm_masked.nii.gz
done

