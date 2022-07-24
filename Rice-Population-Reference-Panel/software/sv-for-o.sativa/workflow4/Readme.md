workflow 4: call INV, by numcer and syri;
Uploaded by  Yong, Andrea, Noor and Saule, on 2021-05-25

Usage:
1. chagne your email address;
2. edit the reference genome and query genome only, in both the workflow4a and workflow4b;
3. run workflow4a;
4. run workflow4b; when workflow4a is fished;

Examples:

Phase1: **sbatch workflow4a**

Phase2: **sbatch workflow4b**

There will be 2 summary results in the output dir, i.e. 
genome1.genome2syri.summary: all the SV results; 
genome1_genome2.INV.summary: the inversions only;
