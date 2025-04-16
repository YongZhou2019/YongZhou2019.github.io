#!/bin/bash
#for h in "MH63RS3" "ZS97RS3" "117425" "125619" "125827" "127518" "127564" "127652" "127742" "128077" "132278" "132424" "Azucena" "IR64RS2" "Nipponbare" "117534":
#for h in "ZS97RS3" "117425" "125619" "125827" "127518" "127564" "127652" "127742" "128077" "132278" "132424" "Azucena" "IR64RS2" "117534" ;
for h in "MH63RS3" "Nipponbare" ;
do
	#echo "Reading column header: " ${h}
	#python GetAbsentGenesForeachGenome.py ${h} "unique"
	#python GetAbsentGenesForeachGenome.py ${h}
	python AnnotateMagic16Genomes_indels_dev_testingdispensable_parallel.py ${h}
done
