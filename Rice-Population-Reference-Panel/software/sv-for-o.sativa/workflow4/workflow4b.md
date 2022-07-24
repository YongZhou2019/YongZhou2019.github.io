**The header lines are used for submit a jog on ibex by sbatch**

```
#!/bin/bash
#SBATCH --job-name=workflow4 #Your Job Name
#SBATCH --nodes=1 #Number of Nodes desired e.g 1 node
#SBATCH --time=50:00:00 #Walltime: Duration for the Job to run HH:MM:SS
#SBATCH --mail-user=yong.zhou@kaust.edu.sa #Your Email address assigned for your job
#SBATCH --mail-type=ALL #Receive an email for ALL Job Statuses
#SBATCH --error=JobName.%J.log #The .log file name
#SBATCH --output=JobName.%J.out #The .output file name
#SBATCH --cpus-per-task=32 #Number of cores assigned per task
#SBATCH --mem=20GB #Memory requested for e.g 2GB
```

** scripts;**

the output will be in $INV/"$Ref"_"$Qry".INV.summary

```
###============generate inversion vcf files================
CurPath=`pwd` ### please keep the genomes in this path;

Ref=genome1 ### could be edit as input parameters, I personally do it here and also could be edit as loops as well;
Qry=genome2

INV=$CurPath/01_genome-wide/"$Ref"."$Qry"/"$Ref"."$Qry".m_i80_l100.delta.SV

cat $INV/"$Ref"."$Qry"syri.vcf | awk 'NR==31{print}' > $INV/"$Ref"_"$Qry"syri.vcf.tmp.head
###two levels are observed；
###Level1
cat $INV/"$Ref"."$Qry"syri.vcf | grep 'PASS' | grep '<INV>' > $INV/"$Ref"_"$Qry"syri.invertion_Level1.tmp1
###Level2
cat $INV/"$Ref"."$Qry"syri.vcf | grep 'PASS' | grep '<INVAL>' > $INV/"$Ref"_"$Qry"syri.invertion_Level2.tmp1
###2level together
cat $INV/"$Ref"_"$Qry"syri.invertion_Level1.tmp1 $INV/"$Ref"_"$Qry"syri.invertion_Level2.tmp1 | sort -k1,1 -k2,2n > $INV/"$Ref"_"$Qry"syri.invertion_2Level.tmp1
cat $INV/"$Ref"_"$Qry"syri.vcf.tmp.head $INV/"$Ref"_"$Qry"syri.invertion_2Level.tmp1 > $INV/"$Ref"_"$Qry"syri.invertion_2Level.tmp2
###change the format；
cat $INV/"$Ref"_"$Qry"syri.invertion_2Level.tmp2 | sed 's/END=/SVTYPE=INV;END=/g' > $INV/"$Ref"_"$Qry"syri.invertion_2Level.tmp3
cat $INV/"$Ref"_"$Qry"syri.invertion_2Level.tmp3 | awk -F '[\t=;]' 'NR!=1 {print $1"\t"$2"\t"$13"|"$15"-"$17"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"="$9";"$10"="$11";"$12"="$13";"$14"="$15";"$16"="$17";"$18"="$19";"$20"="$21";"$21"="$22}' > $INV/"$Ref"_"$Qry"syri.invertion_2Level.tmp4
cat $INV/"$Ref"_"$Qry"syri.vcf.tmp.head $INV/"$Ref"_"$Qry"syri.invertion_2Level.tmp4 > $INV/"$Ref"_"$Qry"syri.invertion_2Level.tmp5
mv $INV/"$Ref"_"$Qry"syri.invertion_2Level.tmp5 $INV/"$Ref"_"$Qry"syri.invertion_2Level.vcf
###Level1
cat $INV/"$Ref"_"$Qry"syri.invertion_2Level.vcf | grep '<INV>' > $INV/"$Ref"_"$Qry"syri.invertion_Level1.tmp6
cat $INV/"$Ref"_"$Qry"syri.vcf.tmp.head $INV/"$Ref"_"$Qry"syri.invertion_Level1.tmp6 > $INV/"$Ref"_"$Qry"syri.invertion_Level1.vcf
cat $INV/"$Ref"_"$Qry"syri.invertion_Level1.vcf | awk -F '[\t=;|-]' 'NR!=1 {print $1"\t"$2"\t"$13"\t"$3"\t"$4"\t"$5"\t"$13-$2"\t"$5-$4}' > $INV/"$Ref"_"$Qry"syri.invertion_Level1.tmp8
echo Ref_Chr"\t"Ref_Start"\t"Ref_End"\t"Qry_Chr"\t"Qry_Start"\t"Qry_End"\t"Ref_Length"\t"Qry_Length > $INV/"$Ref"_"$Qry"syri.invertion_Level1.tmp9
cat $INV/"$Ref"_"$Qry"syri.invertion_Level1.tmp9 $INV/"$Ref"_"$Qry"syri.invertion_Level1.tmp8 > $INV/"$Ref"_"$Qry"syri.invertion_Level1.bed
###Level2
cat $INV/"$Ref"_"$Qry"syri.invertion_2Level.vcf | grep 'INVAL' > $INV/"$Ref"_"$Qry"syri.invertion_Level1.tmp7
cat $INV/"$Ref"_"$Qry"syri.vcf.tmp.head $INV/"$Ref"_"$Qry"syri.invertion_Level1.tmp7 > $INV/"$Ref"_"$Qry"syri.invertion_Level2.vcf
cat $INV/"$Ref"_"$Qry"syri.invertion_Level2.vcf | awk -F '[\t=;|-]' 'NR!=1 {print $1"\t"$2"\t"$13"\t"$3"\t"$4"\t"$5"\t"$13-$2"\t"$5-$4}' > $INV/"$Ref"_"$Qry"syri.invertion_Level2.tmp8
echo Ref_Chr"\t"Ref_Start"\t"Ref_End"\t"Qry_Chr"\t"Qry_Start"\t"Qry_End"\t"Ref_Length"\t"Qry_Length > $INV/"$Ref"_"$Qry"syri.invertion_Level2.tmp9
cat $INV/"$Ref"_"$Qry"syri.invertion_Level2.tmp9 $INV/"$Ref"_"$Qry"syri.invertion_Level2.tmp8 > $INV/"$Ref"_"$Qry"syri.invertion_Level2.bed


rm $INV/*tmp*
### count and summary the INV; 
### the output will be in $INV/"$Ref"_"$Qry".INV.summary
a1=`cat $INV/"$Ref"_"$Qry"syri.invertion_Level1.bed | wc -l`
a2=`expr $a1 - 1`
a3=`cat $INV/"$Ref"_"$Qry"syri.invertion_Level1.bed | awk '{sum += $7};END {print sum}'`
a4=`cat $INV/"$Ref"_"$Qry"syri.invertion_Level1.bed | awk '{sum += $8};END {print sum}'`

echo "The number of Inversions between" $Ref.$Qry  $a2 
echo "The length of Insertion in Reference" $Ref.$Qry  $a3 
echo "The length of Insertion in Query" $Ref.$Qry  $a4 

echo "The number of Inversions between" $Ref.$Qry  $a2 > $INV/"$Ref"_"$Qry".INV.summary
echo "The length of Insertion in Reference" $Ref.$Qry  $a3 >> $INV/"$Ref"_"$Qry".INV.summary
echo "The length of Insertion in Query" $Ref.$Qry  $a4 >> $INV/"$Ref"_"$Qry".INV.summary

echo "DONE"
```

