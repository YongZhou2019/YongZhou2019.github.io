### The header lines are used for submit a jog on ibex by sbatch

```
#!/bin/bash
#SBATCH --job-name=workflow3 #Your Job Name
#SBATCH --nodes=1 #Number of Nodes desired e.g 1 node
#SBATCH --time=50:00:00 #Walltime: Duration for the Job to run HH:MM:SS
#SBATCH --mail-user=yong.zhou@kaust.edu.sa #Your Email address assigned for your job
#SBATCH --mail-type=ALL #Receive an email for ALL Job Statuses
#SBATCH --error=JobName.%J.log #The .log file name
#SBATCH --output=JobName.%J.out #The .output file name
#SBATCH --cpus-per-task=32 #Number of cores assigned per task
#SBATCH --mem=20GB #Memory requested for e.g 2GB
```
### Load the modules and set the tool paths that we need in this workflow3
```
module load minimap2
module load syri/1.4
module load samtools
### set the path and input, we have 4 parameters need to be changed;
genome=/ibex/scratch/projects/c2016/SV_genomes/workflow3/Genome   ### modify 1
mkdir -p /ibex/scratch/projects/c2016/SV_genomes/workflow3/SVcall   ### modify 2
SVcall=/ibex/scratch/projects/c2016/SV_genomes/workflow3/SVcall  ### modify 2
Script=syri -c $SVcall/$Ref.$Qry.bam -r $genome/$Ref.chr -q $genome/$Qry.chr -k -F B
```

### could be edit as input parameters, I personally do it here and also could be edit as loops as well;
```
Ref=genome1 ### modify 3 
### could be edit as a loop, then can be easily with couple of queries;
Qry=genome2 ### modify 4
```

### this part is used for run minimap2 alignment, which will be used for SV calling; Also here will generate a quick dotplot as a generall comparasion
```
minimap2 -ax asm5 --eqx $genome/$Ref.chr $genome/$Qry.chr > $SVcall/$Ref.$Qry.sam
samtools view -b $SVcall/$Ref.$Qry.sam > $SVcall/$Ref.$Qry.bam
samtools index $SVcall/$Ref.$Qry.bam
```
### this part is used for run nucmer alignment, which will be used for SV calling; Also here will generate a quick dotplot as a generall comparasion
```
syri -c $SVcall/$Ref.$Qry.bam -r $genome/$Ref.chr -q $genome/$Qry.chr -k -F B --dir $SVcall --prefix $Ref.$Qry
plotsr $SVcall/$Ref.$Qry/yri.out $genome/$Ref.chr $genome/$Qry.chr -H 8 -W 5
mv syri.pdf $SVcall/$Ref.$Qry.syri.pdf
```

```
date
echo "workflow3a DONE"
```



