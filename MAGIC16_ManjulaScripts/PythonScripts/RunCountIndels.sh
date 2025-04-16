#!/bin/bash
#for h in "MH63RS3" "ZS97RS3" "117425" "125619" "125827" "127518" "127564" "127652" "127742" "128077" "132278" "132424" "Azucena" "IR64RS2" "Nipponbare" "117534":
#for h in "125619" "125827" "127518" "127564" "127652" "127742" "128077" "132278" "132424" "Azucena" "IR64RS2" "Nipponbare" "117534":
for h in "117534" ;
do
	echo "Reading column header: " ${h}
	python AnnotateMagic16Genomes_indels_dev_testingcore_parallel.py ${h} ~/MAGIC16/gramene_genome33_Only_YongUpdatedFromZhiChao/1.IR64RS2.groupset.sorted.core.16.homologs
done
