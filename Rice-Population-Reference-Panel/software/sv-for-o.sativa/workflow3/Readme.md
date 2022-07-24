workflow 3: call SVs, including INV, by minimap2 and syri;
Uploaded by  Yong, Andrea, Noor and Saule, on 2022-01-02
Usage:

chagne your email address;
edit the reference genome and query genome only, in both the workflow3a and workflow3b;
run workflow3a;
run workflow3b; when workflow3a is fished;

Examples:
Phase1: sbatch workflow4a
Phase2: sbatch workflow4b
There will be 2 summary results in the output dir, i.e.
genome1.genome2syri.summary: all the SV results;
genome1_genome2.INV.summary: the inversions only;
