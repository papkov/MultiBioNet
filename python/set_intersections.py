from matplotlib_venn import venn3
from matplotlib import pyplot as plt, pylab


def set_of_ensg(file_path, flag):
    s = set()
    with open(file_path, "r") as f:
        if flag == "ppi":
            lines = f.readlines()[1:]
            for ppi in lines:
                ppi = ppi.split()
                s.add(ppi[0])
                s.add(ppi[1])

        elif flag == "GWAS":
            lines = str(f.readlines()[0])
            lines = lines.split('\r')
            for genes in lines:
                genes = genes.strip().split()
                s.add(genes[1])

        elif flag == "pathway":
            lines = str(f.readlines()[0])
            lines = lines.split('\r')
            for genes in lines:
                genes = genes.strip().split()
                if genes:
                    s.add(genes[0])
    return s


alz_intact_set = set_of_ensg("../data/alz_intact_int_PPI.txt", "ppi")
alz_coexpr_set = set_of_ensg("../data/alzcoexp_int_00001_10102016.txt", "ppi")
all_intact_set = set_of_ensg("../data/intact_int.txt", "ppi")
parkinson_intact_set = set_of_ensg("../data/parkinson_intact_int_PPI.txt", "ppi")
synapse_intact_set = set_of_ensg("../data/synapse_intact_int.txt", "ppi")

alz_gwas_set = set_of_ensg("../data/GWAS_genes.txt", "GWAS")
alz_pathway_set = set_of_ensg("../data/alzheimer_pathway_genes.txt", "pathway")


plt.figure(figsize=(10, 10))
venn3([all_intact_set, alz_gwas_set, alz_pathway_set], ('all_intact', 'alz_gwas', 'alz_pathway'))
pylab.savefig('all+alz-genes_datasets.png')
