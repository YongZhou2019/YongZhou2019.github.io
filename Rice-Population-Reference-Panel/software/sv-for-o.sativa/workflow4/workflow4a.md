**The header lines are used for submit a jog on ibex by sbatch**

```
#!/bin/bash
#SBATCH --job-name=workflow4 #Your Job Name
#SBATCH --nodes=1 #Number of Nodes desired e.g 1 node
#SBATCH --time=50:00:00 #Walltime: Duration for the Job to run HH:MM:SS
#SBATCH --mail-user=user@kaust.edu.sa #Your Email address assigned for your job
#SBATCH --mail-type=ALL #Receive an email for ALL Job Statuses
#SBATCH --error=JobName.%J.log #The .log file name
#SBATCH --output=JobName.%J.out #The .output file name
#SBATCH --cpus-per-task=32 #Number of cores assigned per task
#SBATCH --mem=20GB #Memory requested for e.g 2GB
```

**Load the modules and set the tool paths that we need in this workflow4**

```
module load mummer/4.0.0.2/gnu-6.4.0
module load gnuplot/5.0.0
module load syri/1.4
```

**could be edit as input parameters, I personally do it here and also could be edit as loops as well;**

```
mummerplot=/ibex/scratch/projects/c2016/SV_genomes/workflow4/mummerplot
### 
Ref=genome1 ### could be edit as input parameters, I personally do it here and also could be edit as loops as well;
Qry=genome2 ### change your refernce and queries;
```

**this part is used for run nucmer alignment, which will be used for SV calling; Also here will generate a quick dotplot as a generall comparasion**
```
CurPath=`pwd` ### please keep the genomes in this path;
### generate the working path for Ref and Qry
mkdir -p $CurPath/01_genome-wide/$Ref.$Qry  
WorkingPath=$CurPath/01_genome-wide/$Ref.$Qry
echo "========nucmer begin======="
### compare the Qry to Ref
nucmer -t 32 -p $WorkingPath/$Ref.$Qry $CurPath/$Ref.chr $CurPath/$Qry.chr
### show the quick comparasion with dotplot
### comparasion filter
delta-filter -i 95 -l 15000 -o 95 $WorkingPath/$Ref.$Qry.delta -1 > $WorkingPath/$Ref.$Qry.best.delta
### show coords
show-snps -Clr $WorkingPath/$Ref.$Qry.best.delta > $WorkingPath/$Ref.$Qry.best.snps
$mummerplot -p $WorkingPath/$Ref.$Qry.best $WorkingPath/$Ref.$Qry.best.delta -t postscript
ps2pdf $WorkingPath/$Ref.$Qry.best.ps $WorkingPath/$Ref.$Qry.best.pdf
convert -density 300 $WorkingPath/$Ref.$Qry.best.pdf $WorkingPath/$Ref.$Qry.best.png
echo "========nucmer done======="
```
**syri call Structure varations, there is a summary result there**
```
echo "========syri call SV ======="
### nucmer result filter again for INV calling; modify here if necessary, i.e. iden=80, length=100bp
delta-filter -m -i 80 -l 100 $WorkingPath/$Ref.$Qry.delta > $WorkingPath/$Ref.$Qry.m_i80_l100.delta
### show coords, human read
show-coords -THrd $WorkingPath/$Ref.$Qry.m_i80_l100.delta > $WorkingPath/$Ref.$Qry.m_i80_l100.coords

### run syri call , there is a summary result there;
mkdir -p $WorkingPath/$Ref.$Qry.m_i80_l100.delta.SV
SVPath=$WorkingPath/$Ref.$Qry.m_i80_l100.delta.SV
syri -c $WorkingPath/$Ref.$Qry.m_i80_l100.coords -d $WorkingPath/$Ref.$Qry.m_i80_l100.delta -r $Ref.chr -q $Qry.chr --dir $SVPath --prefix $Ref.$Qry --nc 20

echo "syri finished on"
date
```
