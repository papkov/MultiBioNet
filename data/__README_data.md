Dataset consists of 5 files.
First 4 files describe protein-protein interactions(PPI) in homo sapiens. They are obtained from IntAct database(http://www.ebi.ac.uk/intact/).

Selected interaction are filtered based on MI score >= 0.45 from IntAct and considered to be highly confident.
Protein names are translated to unique Ensemble(http://www.ensembl.org/index.html) gene identifiers (ENSG). 

1. Dataset containing expert curated interactions related to Parkinson's disease parkinson_intact_int_PPI.txt
2. All highly confident interactions in hs available in IntAct database intact_int.txt
3. Aumatically curated interactions related to synaptic activity synapse_intact_int.txt
4. Dataset containing expert curated interactions related to Alzheimer's diseasealz_intact_int_PPI.txt

The 5th file describes genes, that are coexpressed  in the Alzheimer's patients and healthy individuals.
Correlation profiles are calculated using the set of microarray data obtained from ArrayExpress database (https://www.ebi.ac.uk/arrayexpress/). Final coexpression score is calculated by applying Robust Rank Aggregation method to aggregate correlation scores in individual microarray datasets.

Coexpression scores are filtered based on pvalue <=0.00001.
Gene names are translated to unique Ensemble(http://www.ensembl.org/index.html) gene identifiers (ENSG). 

alzcoexp_int_31052016.txt