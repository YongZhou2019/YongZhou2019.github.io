### workflow 2a: call SV by NGMLR + sniffles for each/mutiple genome(s);
### workflow 2b: merge mutiple genomes' SVs
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
#SBATCH --time=1:00:00
#SBATCH --mem=50GB #Memory requested for e.g 2GB
```

### installed the tool/paths in the following path
```
SURVIVOR=/ibex/scratch/projects/c2016/SV_genomes/workflow4/SURVIVOR/Debug
SVcall=/ibex/scratch/projects/c2016/SV_genomes/workflow2/SVcall  ### the dir of output from workflow2a
Ref=genome1
```



### merge all Insersion;
```
ls $SVcall/$Ref.*Insertion.vcf > $SVcall/$Ref.Insertion.tmp.txt
$SURVIVOR merge $SVcall/$Ref.Insertion.tmp.txt  50 1 1 -1 -1 -1 $SVcall/magic16.$Ref.Insertion.merge.tmp.vcf ### 50bp means we allow 50bp differnece at the edge/breakpoint of the SVs
rm $SVcall/$Ref.Insertion.tmp.txt
```


### merge all Deletion;
```
ls $SVcall/$Ref.*Deletion.vcf > $SVcall/$Ref.Deletion.tmp.txt
$SURVIVOR merge $SVcall/$Ref.Deletion.tmp.txt  50 1 1 -1 -1 -1 $SVcall/magic16.$Ref.Deletion.merge.tmp.vcf ### 50bp means we allow 50bp differnece at the edge/breakpoint of the SVs
rm $SVcall/$Ref.Deletion.tmp.txt
```


### more parameters about SURVIVOR please visit, https://github.com/fritzsedlazeck/SURVIVOR


```
date
echo "workflow2b DONE"
```
