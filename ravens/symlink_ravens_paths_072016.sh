# 7/20/16, script to symlink ravens images to copy over to cbica
# by MP

input_list=/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/ravens/n425_ravens_for_midas_072016.csv

#rm -f $tbvFile

for i in $(cat $input_list); do

        echo "--- image path ---"

        imagePath=$(echo $i|cut -d, -f1)
        echo $imagePath

	if [ ! -e "$imagePath" ]; then
		echo "--- IMAGE IS MISSING ---"

	elif [ -e "$imagePath" ]; then
	
		ln -sf $imagePath /data/joy/BBL/projects/pehlivanovaPncItc/RAVENS/tmpRavenSymLink/

	fi 

done































	
