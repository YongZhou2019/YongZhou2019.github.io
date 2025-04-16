#!/bin/bash

#for fname in `ls /home/thimmamp/MAGIC16/indels_results/*_Core_Noduplication_genic_InDelsStats_Details_testing.txt`
#for fname in `ls /home/thimmamp/MAGIC16/indels_results/genome*_Dispensable_with_unique__orthogenes_genic_InDelsStats_Details_testing.txt`
#for fname in `ls /home/thimmamp/MAGIC16/indels_results/genome*Dispensable_with_all__orthogenes_genic_InDelsStats_Details_testing.txt`
for fname in `ls /home/thimmamp/MAGIC16/indels_results/genome*_Core_Allwithduplication_orthogenes_genic_InDelsStats_Details_testing.txt`
do
	#echo "Reading file: " $fname
	INS=`awk -F'\t' 'BEGIN{total=0}{total+=$3}END{print total}' $fname`
	DELS=`awk -F'\t' 'BEGIN{total=0}{total+=$4}END{print total}' $fname`
	echo $fname $INS $DELS
done
