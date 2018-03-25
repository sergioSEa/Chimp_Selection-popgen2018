from Bio import SeqIO
from Bio import Entrez
import numpy as np
import re

### GET LIST OF GENES

pattern = 'NC_\d+'
pattern = re.compile(pattern)

snps=[]
f = open('../ChimpExome/PBSrelevant.txt', 'r')
for line in f:
    line=line.split(sep=" ")
    snps.append([float(line[0]),float(line[1]),float(line[2].rstrip())])
snps = np.array(snps)

zeros_c = np.zeros((50,1))
snps = np.append(snps,zeros_c,axis=1)
print(snps)
# Separated list of genes
gene_list = []


new_pbs_gene = []

for seq_record in SeqIO.parse("../human_gb/GCF_000001405.37_GRCh38.p11_genomic.gbff",
                        "genbank"):
    if re.match(pattern, seq_record.name):
        name = seq_record.name
        name = name.split("_")[1]
        chromosome = int(name)
        snpsChr = snps[snps[:,1]==chromosome]
        for feat in seq_record.features:
            for snp in snpsChr:
                if snp[2] > feat.location.start and snp[2] < feat.location.end:
                    # Gene found --> get ID
                    gene = feat.qualifiers.get('db_xref')
                    for id_gene in gene:
                        id_gene = id_gene.split(':')
                        if id_gene[0] == 'GeneID':
                            gene_list.append(id_gene[1])
                            # snps[snps[:,2]==snp[2],-1] = id_gene[1]
                            item = [str(float(snp[0])), str(int(snp[1])), str(int(snp[2])),str(int(id_gene[1]))]
                            if item not in new_pbs_gene:
                                new_pbs_gene.append(item)



gene_list = set(gene_list)
print(gene_list)
print(new_pbs_gene)
fw = open('../human_gb/gene_PBS.txt', 'w')
for line in new_pbs_gene:
    line = ','.join(line)
    fw.write(line+'\n')
fw.close()

# ['GeneID:143187', 'HGNC:HGNC:17792', 'MIM:614316']
