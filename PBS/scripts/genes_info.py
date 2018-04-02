from Bio import SeqIO
from Bio import Entrez
import numpy as np
import re
import sys

Entrez.email = "lvz544@alumni.ku.dk"

def retrieve_annotation(id_list):

    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""

    request = Entrez.epost("gene", rettype="", retmode="xml", id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print("An error occurred while retrieving the annotations.")
        print ("The error returned was %s" % e)
        sys.exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)

    print("Retrieved %d annotations for %d genes" % (len(annotations), len(id_list)))
    return annotations

def fetch_genes(id_list):
    """Fetch Entrez Gene records using Bio.Entrez, in particular epost
    (to submit the data to NCBI) and efetch to retrieve the
    information, then use Entrez.read to parse the data.
    Returns a list of parsed gene records.
    """

    request = Entrez.epost("gene",id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print("An error occurred while retrieving the annotations.")
        print("The error returned was %s" % e)
        sys.exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    efetch_result = Entrez.efetch(db="gene", webenv=webEnv, query_key = queryKey, retmode="xml")
    genes = Entrez.read(efetch_result)
    #print "Retrieved %d records for %d genes" % (len(genes),len(id_list))
    return genes

def print_data(annotation):
    for gene_data in annotation:
        gene_id = gene_data["Id"]
        gene_symbol = gene_data["NomenclatureSymbol"]
        gene_name = gene_data["Description"]
        print("ID: %s - Gene Symbol: %s - Gene Name: %s" % (gene_id, gene_symbol, gene_name))

data = np.loadtxt('../human_gb/gene_PBS.txt')

gene_list = []
for snp in data:
    gene_id = int(snp[3])
    if gene_id == 0:
        continue
    gene_list.append(str(gene_id))
print(gene_list)

gene_list = ['9927']


# Genes
# genes = fetch_genes(gene_list)
# print(genes)
# for g in genes:
#     print(g)
#     print(g['Entrezgene_gene'])
#     print("\n\n")
    #print("%s\t%s\t%s\t%s\t%s" % (g["taxname"],g["entrez_id"],g["official_symbol"],",".join(g["refseq_ids"]),g["official_full_name"]))


#
# print("\n\n\n\n")
# Annotations
annotations = retrieve_annotation(gene_list)
for annot in annotations['DocumentSummarySet']['DocumentSummary']:
    print(annot)
    print('\n---------\n')



# {'Name': 'MFN2', 'Description': 'mitofusin 2', 'Status': '0', 'CurrentID': '0',
# 'Chromosome': '1', 'GeneticSource': 'genomic', 'MapLocation': '1p36.22',
# 'OtherAliases': 'CMT2A, CMT2A2, CMT2A2A, CMT2A2B, CPRP1, HMSN6A, HSG, MARF',
# 'OtherDesignations': 'mitofusin-2|hyperplasia suppressor|mitochondrial assembly regulatory factor|transmembrane GTPase MFN2',
# 'NomenclatureSymbol': 'MFN2', 'NomenclatureName': 'mitofusin 2',
# 'NomenclatureStatus': 'Official', 'Mim': ['608507'],
# 'GenomicInfo': [{'ChrLoc': '1', 'ChrAccVer': 'NC_000001.11', 'ChrStart': '11980180', 'ChrStop': '12013514', 'ExonCount': '20'}],
# 'GeneWeight': '22397',
# 'Summary': 'This gene encodes a mitochondrial membrane protein that participates in mitochondrial fusion and contributes to the maintenance and operation of the mitochondrial network. This protein is involved in the regulation of vascular smooth muscle cell proliferation, and it may play a role in the pathophysiology of obesity. Mutations in this gene cause Charcot-Marie-Tooth disease type 2A2, and hereditary motor and sensory neuropathy VI, which are both disorders of the peripheral nervous system. Defects in this gene have also been associated with early-onset stroke. Two transcript variants encoding the same protein have been identified. [provided by RefSeq, Jul 2008]',
# 'ChrSort': '01', 'ChrStart': '11980180',
# 'Organism': {'ScientificName': 'Homo sapiens', 'CommonName': 'human', 'TaxID': '9606'}
#
