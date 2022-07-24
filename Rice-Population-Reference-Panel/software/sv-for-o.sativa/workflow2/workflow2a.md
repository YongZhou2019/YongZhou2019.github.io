### workflow 2: call SV, by NGMLR + sniffles, i.e. INV and DEL (PAV);
### workflow 2a: call SV by each/mutiple genome(s)
### use 12 chromsoomes sequence as genome input;


```
#!/bin/bash
#SBATCH --job-name=workflow2 #Your Job Name
#SBATCH -N 1
#SBATCH --cpus-per-task=30
#SBATCH --partition=batch
#SBATCH --error=JobName.%J.log #The .log file name
#SBATCH --output=JobName.%J.out #The .output file name
#SBATCH --mail-user=yong.zhou@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH --time=20:00:00
#SBATCH --mem=50GB #Memory requested for e.g 2GB
```

### load modules/tools
```
module load emboss
module load samtools
module load ngmlr/0.2.7/anaconda3env
module load sniffles
```


### set the path and input, we have 4 parameters need to be changed;
```
genome=/ibex/scratch/projects/c2016/SV_genomes/workflow2/Genome   ### modify 1
mkdir -p /ibex/scratch/projects/c2016/SV_genomes/workflow2/SVcall   ### modify 2
SVcall=/ibex/scratch/projects/c2016/SV_genomes/workflow2/SVcall  ### modify 2
```

```
### could be edit as input parameters, I personally do it here and also could be edit as loops as well;
Ref=genome1 ### modify 3 
### could be edit as a loop, then can be easily with couple of queries;
Qry=genome2 ### modify 4
```


### run workflow2a;
```
### Here we used genome1.chr as the reference and genome2.chr as the query (both with only 2 chromsomes);
### split the query genome by window size:50Kb and step window size: 5Kb; this is the best option for balance the accurate and time comsum
splitter -sequence $genome/$Qry.chr -outseq $genome/$Qry.fasta -size 50000 -overlap 45000 
### run mapping 
ngmlr -t 26 -r $genome/$Ref.chr -q $genome/$Qry.fasta --skip-write | samtools view -bS | samtools sort -o $SVcall/$Ref.$Qry.sorted.bam
### run SV calling
sniffles -m $SVcall/$Ref.$Qry.sorted.bam -v $sniSV/$Ref.$Qry.sniffles.vcf -s 6 -t 4 -l 50 #-b $Ref.$Qry.bed
```


```
### rm the tmp files
# rm $SVcall/$Ref.$Qry.sorted.bam ###  but could be keep for a while for mapping check;
rm $genome/$Qry.fasta
rm -rf $SVcall/$Ref.$Qry
```

```
## sort out insertion; please note that the results have the only coordinates based on Reference
cat $SVcall/$Ref.$Qry.sniffles.vcf | grep 'CHROM' > $SVcall/$Ref.$Qry.Insertion.vcf
cat $SVcall/$Ref.$Qry.sniffles.vcf | grep 'INS' | grep 'PASS' >> $SVcall/$Ref.$Qry.Insertion.vcf
### sort out Deletion; please note that the results have the only coordinates based on Reference
cat $SVcall/$Ref.$Qry.sniffles.vcf | grep 'CHROM' > $SVcall/$Ref.$Qry.Deletion.vcf
cat $SVcall/$Ref.$Qry.sniffles.vcf | grep '<DEL>' | grep 'PASS' >> $SVcall/$Ref.$Qry.Deletion.vcf
### the output is stored in $SVcall/$Ref.$Qry.Insertion.vcf and $SVcall/$Ref.$Qry.Deletion.vcf
```


```
### convert to bed files
cat $SVcall/$Ref.$Qry.Insertion.vcf | awk -F '[\t=;]' '{print $1,$2,$14,sqrt($28*$28),$3}' | sed 's/ /    /g' > $SVcall/$Ref.$Qry.Insertion.bed
cat $SVcall/$Ref.$Qry.Deletion.vcf | awk -F '[\t=;]' '{print $1,$2,$14,sqrt($28*$28),$3}' | sed 's/ /   /g' > $SVcall/$Ref.$Qry.Deletion.bed

### count and summary the INS and DELs; 
### the output will be in $SVcall/$Ref.$Qry.SV.summary
a1=`cat $SVcall/$Ref.$Qry.Insertion.vcf | wc -l`
a2=`expr $a1 - 1`
a3=`cat $SVcall/$Ref.$Qry.Insertion.bed | awk '{sum += $4};END {print sum}'`
echo "The number of Insertion between" $Ref.$Qry  $a2 > $SVcall/$Ref.$Qry.SV.summary
echo "The length of Insertion between" $Ref.$Qry  $a3 >> $SVcall/$Ref.$Qry.SV.summary

b1=`cat $SVcall/$Ref.$Qry.Deletion.vcf | wc -l`
b2=`expr $b1 - 1`
b3=`cat $SVcall/$Ref.$Qry.Deletion.bed | awk '{sum += $4};END {print sum}'`

echo "The number of Deletion of" $Ref.$Qry  $b2 >> $SVcall/$Ref.$Qry.SV.summary
echo "The length of Deletion of" $Ref.$Qry  $b3 >> $SVcall/$Ref.$Qry.SV.summary
```


```
date
echo "workflow2a DONE"
```


