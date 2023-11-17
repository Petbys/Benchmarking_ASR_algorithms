# Benchmarking_ASR_algorithms
Benchmarking Ancestral Sequence Reconstruction algorithms with Ancient DNA

# Data preparation
## protein coding sequences
To benchmark the different algorithms over different datases, the data has to be in the same sequence space. As the main utility for ancestral sequence reconstruction lies in protein coding space the utilized data sets will be converted into protein coding space. 

The _y. Sernia_ data sets consist of aligned snps between 250 different strains. To translate the snps into protein coding space, all corresponding codons to the snps are found using $snp_to_codon.py$. in this program the corresponding gene to each snp is found and based on its start, stop location and orientation the remaining codon bases are extracted from the reference genome to construct an codon alignment fasta file.

Intital files:

|Yrsernia Pesti|full alignment|snp alignment|Information file|
|--|--|--|--|
|2022 spyrou|.fasta, 252 species|.fasta, 252 species|.xlsx, snp locations, decent|
|2018 spyrou||.fasta||





## ASR algorithms and software to be investigate

|Name|Function|Availability|Language|
|--|--|--|--|
|PAML||||
|FastML||||
|ARPIP|ASR using Poisson Indel Process, reconstructing ancestral sequences based on indels |https://github.com/acg-team/bpp-ARPI|C++|

## Data Availability

### Ancient Mammals
#### Mammuth
|Species name|Amount|Format|Availability|
|--|--|--|--|
|Mammuth|16 newly sequenced,23 in total|FastQ|Newly sequenced: https://www.ebi.ac.uk/ena/browser/view/PRJEB59491, others belong to other IDs|
|Elephant|28 individuals|FastQ|28 individual entries|

#### Human
|Species name|Amount|Format|Availability|
|--|--|--|--|
|Human|>1.000.000 SNPs|.anno?|[AADR](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW)|

### Ancient Shotgun metagenomics 
|Type|Function|
|--|--|
|AncientMetagenomeDir| Curated collection of all published data using shotgun metagenomics or microbial genomes|

### Bacteria
