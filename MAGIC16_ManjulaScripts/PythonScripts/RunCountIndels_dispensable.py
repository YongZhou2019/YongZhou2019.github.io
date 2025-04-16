#!/bin/bash
## Core header
for h in "MH63RS3" "ZS97RS3" "117425" "125619" "125827" "127518" "127564" "127652" "127742" "128077" "132278" "132424" "Azucena" "IR64RS2" "Nipponbare" "117534":
do
	echo "Reading column header: " ${h}
	python AnnotateMagic16Genomes_indels_dev_testingdispensable_parallel.py ${h} ~/MAGIC16/gramene_genome33_Only_YongUpdatedFromZhiChao/2.IR64RS2.groupset.sorted.dispensable.16.homologs "all"
        python AnnotateMagic16Genomes_indels_dev_testingdispensable_parallel.py ${h} ~/MAGIC16/gramene_genome33_Only_YongUpdatedFromZhiChao/2.IR64RS2.groupset.sorted.dispensable.16.homologs "unique"
done
