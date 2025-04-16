import sys, os
from matplotlib import pyplot as plt
from matplotlib_venn import venn3
from matplotlib_venn import venn2
import itertools

## usage : thimmamp@kw61043:~/PythonPrograms/utilities$ python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome27_SNP_7K /data_b/IRGSP_SNP_InDels_Venn/genome27_SNP_3K
## usage: python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome27_SNP_3K /data_b/IRGSP_SNP_InDels_Venn/genome27_SNP_7K /data_b/IRGSP_SNP_InDels_Venn/genome27_SNP_10K 
#Section 3K: 2066255 entries
#Section 7K: 6377315 entries
#Section 10K: 11807359 entries
#Section 3K_7K: 2077052 entries
#Section 3K_10K: 3067827 entries
#Section 7K_10K: 10139758 entries
#Section 3K_7K_10K: 21223929 entries
## usage : thimmamp@kw61043:~/PythonPrograms/utilities$ python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome1.3K.biallelic.base.genomewide /data_b/IRGSP_SNP_InDels_Venn/genome1.7K.biallelic.base.genomewide /data_b/IRGSP_SNP_InDels_Venn/genome1.10K.biallelic.base.genomewide 
## usage : python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome1.3K.biallelic.genomewide.SNPs  /data_b/IRGSP_SNP_InDels_Venn/genome1.7K.biallelic.genomewide.SNP /data_b/IRGSP_SNP_InDels_Venn/genome1.10K.biallelic.genomewide.SNPs
## /data_b/IRGSP_SNP_InDels_Venn/genome1.3K.genomewide.SNPs.withID /data_b/IRGSP_SNP_InDels_Venn/genome1.7K.genomewide.SNPs.withID /data_b/IRGSP_SNP_InDels_Venn/genome1.10K.genomewide.SNPs.withID
## python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome1.3K.biallelic.base.missing0 /data_b/IRGSP_SNP_InDels_Venn/genome1.7K.biallelic.base.missing0 /data_b/IRGSP_SNP_InDels_Venn/genome1.3K.biallelic.base.missing0
##python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome1.3K.biallelic.base.missing0.prune2 /data_b/IRGSP_SNP_InDels_Venn/genome1.7K.biallelic.base.missing0.prune2 /data_b/IRGSP_SNP_InDels_Venn/genome1.10K.biallelic.base.missing0.prune2

##  python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome27.3K.biallelic.base.genomewide /data_b/IRGSP_SNP_InDels_Venn/genome27.7K.biallelic.base.genomewide /data_b/IRGSP_SNP_InDels_Venn/genome27.10K.biallelic.base.genomewide
## python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome27.3K.biallelic.genomewide.SNPs /data_b/IRGSP_SNP_InDels_Venn/genome27.7K.biallelic.genomewide.SNPs /data_b/IRGSP_SNP_InDels_Venn/genome27.10K.biallelic.genomewide.SNPs
## python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome27.3K.genomewide.SNPs.withID //data_b/IRGSP_SNP_InDels_Venn/genome27.7K.genomewide.SNPs.withID /data_b/IRGSP_SNP_InDels_Venn/genome27.10K.genomewide.SNPs.withID
## python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome27.3K.biallelic.base.missing0 /data_b/IRGSP_SNP_InDels_Venn/genome27.7K.biallelic.base.missing0 /data_b/IRGSP_SNP_InDels_Venn/genome27.10K.biallelic.base.missing0
## python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome27.3K.biallelic.base.missing0.prune2 /data_b/IRGSP_SNP_InDels_Venn/genome27.7K.biallelic.base.missing0.prune2 /data_b/IRGSP_SNP_InDels_Venn/genome27.10K.biallelic.base.missing0.prune2


##  python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome33.3K.biallelic.base.genomewide /data_b/IRGSP_SNP_InDels_Venn/genome33.7K.biallelic.base.genomewide /data_b/IRGSP_SNP_InDels_Venn/genome33.10K.biallelic.base.genomewide
## python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome33.3K.biallelic.genomewide.SNPs /data_b/IRGSP_SNP_InDels_Venn/genome33.7K.biallelic.genomewide.SNPs /data_b/IRGSP_SNP_InDels_Venn/genome33.10K.biallelic.genomewide.SNPs
## python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome33.3K.genomewide.SNPs.withID //data_b/IRGSP_SNP_InDels_Venn/genome33.7K.genomewide.SNPs.withID /data_b/IRGSP_SNP_InDels_Venn/genome33.10K.genomewide.SNPs.withID
## python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome33.3K.biallelic.base.missing0 /data_b/IRGSP_SNP_InDels_Venn/genome33.7K.biallelic.base.missing0 /data_b/IRGSP_SNP_InDels_Venn/genome33.10K.biallelic.base.missing0
## python CreateVennDiagram.py /data_b/IRGSP_SNP_InDels_Venn/genome33.3K.biallelic.base.missing0.prune2 /data_b/IRGSP_SNP_InDels_Venn/genome33.7K.biallelic.base.missing0.prune2 /data_b/IRGSP_SNP_InDels_Venn/genome33.10K.biallelic.base.missing0.prune2

## Venn diagram for new datasets from plink2 and plink Nov'24
#python CreateVennDiagram.py /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.3K_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.7K_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.10K_forvenn 

#python CreateVennDiagram.py /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.3K.biallelic_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.7K.biallelic_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.10K.biallelic_forvenn


#python CreateVennDiagram.py /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.3K.biallelic.base_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.7K.biallelic.base_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.10K.biallelic.base_forvenn

#python CreateVennDiagram.py /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.3K.biallelic.base.missing0.2.maf0.01_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.7K.biallelic.base.missing0.2.maf0.01_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.10K.biallelic.base.missing0.2.maf0.01_forvenn

#python CreateVennDiagram.py /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.3K.biallelic.base.missing0.2_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.7K.biallelic.base.missing0.2_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.10K.biallelic.base.missing0.2_forvenn

#python CreateVennDiagram.py /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.3K.biallelic.base.maf0.01_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.7K.biallelic.base.maf0.01_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.10K.biallelic.base.maf0.01_forvenn

#python CreateVennDiagram.py /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.3K.biallelic.base.missing0.2.maf0.01.prune2_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.7K.biallelic.base.missing0.2.maf0.01.prune2_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1.10K.biallelic.base.missing0.2.maf0.01.prune2_forvenn
#-------


## venn for first set
#python CreateVennDiagram.py /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1_phase5_plink2_3K_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1_phase5_plink2_7K_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1_phase5_plink2_10K_forvenn

#python CreateVennDiagram.py /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1_biallelic_phase5_plink2_3K_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1_biallelic_phase5_plink2_7K_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1_biallelic_phase5_plink2_10K_forvenn

#python CreateVennDiagram.py /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1_biallelic_fc5phase5_plink2_3K_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1_biallelic_fc5phase5_plink2_7K_forvenn /data-A/IRGSP_SNP_InDels_Venn/DataForVennPlots/genome1_biallelic_fc5phase5_plink2_10K_forvenn


location = "/data-A/IRGSP_SNP_InDels_Venn/VennResults/"
#datatype = "genome1_phase5_plink2"
datatype = "genome1_biallelic_phase5_plink2"
#datatype = "genome1_biallelic_phase5_filterCount5_plink2"

# Read the files and extract unique strings
file1_data = set()
file2_data = set()
file3_data = set()

# Read and process file 1
with open(sys.argv[1], "r") as file1:
    for line in file1:
        file1_data.add(line.strip())

# Read and process file 2
with open(sys.argv[2], "r") as file2:
    for line in file2:
        file2_data.add(line.strip())

# Read and process file 3
with open(sys.argv[3], "r") as file3:
    for line in file3:
        file3_data.add(line.strip())

# Generate the Venn diagram
venn = venn3([file1_data, file2_data, file3_data], ('genome1_3K', 'genome1_7K', 'genome1_10K'))
#venn = venn2([file1_data, file2_data], ('genome1_7K_SNP', 'genome1_3K_SNP'))

# Get the entries of each section
#venn_entries = {
#    "100": file1_data - file2_data - file3_data,
#    "010": file2_data - file1_data - file3_data,
#    "001": file3_data - file1_data - file2_data,
#    "110": file1_data & file2_data - file3_data,
#    "101": file1_data & file3_data - file2_data,
#    "011": file2_data & file3_data - file1_data,
#    "111": file1_data & file2_data & file3_data
#}

#venn_entries = {
#    "genome1_10": file1_data - file2_data,
#    "genome1_1": file2_data - file1_data,
#    "genome1_10_genome1_1": file1_data & file2_data
#    }

venn_entries = {
    "3K": file1_data - file2_data - file3_data,
    "7K": file2_data - file1_data - file3_data,
    "10K": file3_data - file1_data - file2_data,
    "3K_7K": file1_data & file2_data - file3_data,
    "3K_10K": file1_data & file3_data - file2_data,
    "7K_10K": file2_data & file3_data - file1_data,
    "3K_7K_10K": file1_data & file2_data & file3_data
}

# Print the entries of each section

##write output
for key, value in venn_entries.items():
    print(f"Section {key}: {len(value)} entries")
    fname = location+datatype+"_"+key+"_entries.txt"
    with open(fname, 'w') as f:
        #f.write("Dataset\tNumber\tEntries\n")
        #f.write(key+'\t'+str(len(value))+'\t'+' '.join(list(venn_entries[key]))+'\n')
        f.write("\n".join(list(venn_entries[key]))+'\n')

# Save the plot to a PDF file
filename=datatype+"3K_7K_10K_Venndiagram.pdf"
plt.savefig(location+filename)
#genome33_3K_7K_10K_biallelic_base_missing0_prune2_Venndiagram.pdf")

# Display the plot
#plt.show()
