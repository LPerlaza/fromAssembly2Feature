# From Assemblies to Feature (Prokaryotes)

These scripts detect features of interest in prokaryotic assemblies. Each perl script uses genome assembly and the sequence of your feature of interest and outputs a descriptive table, the alignment of your feature against the feature in the assembly and the peptide sequences (when applies). These scripts have been designed to emulate manually curated alignmets, they are also execute in parallel. These characteristics allow you to go throught a large number of interesting features without the time comsuming task of exploring each sequence manually. Please reach out if you find inconsistences.

## Requires
 
- The required inputs are assemblies and genes in fasta formats. 
- The genomes must be located in a folder together with not additional files.  
- The genes of interest must be located in a folder together with not additional files. 
- The genes must have the extension ".fasta" as it is used as a tag for handling file in the program.  
- The assemblies can have any kind of extension.  

## Dependencies 
Please make sure you have the following dependencies in your computer:											

- Programs in $PATH: Blast, prodigal and samtools. 
- R packages: [msa](https://bioconductor.org/packages/release/bioc/html/msa.html), [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html), [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html)

optional (only if you want to predict plasmids)

- plasmid prediction: https://github.com/Shamir-Lab/PlasClass
 
## Conda env

You can create a conda environment with all the dependencies:

```
conda create --name fromAssembly2feature
conda activate fromAssembly2feature

conda install -c bioconda prodigal
conda install -c bioconda samtools
conda install -c bioconda blast
conda install -c bioconda plasclass
conda install -c r r-essentials

#open R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("msa")
BiocManager::install("Biostrings")

install.packages("seqinr")
install.packages("reshape2")

```

## INSTALL 
```
git clone git@github.com:LPerlaza/fromAssembly2Feature.git
cd Assembly2Gene
chmod a+x fromAssembly2gene.pl
```	 		

# From Assembly to Gene (fromAssembly2gene.pl)
Detect coding genes in an assembly, and get their alignments and description using fromAssembly2gene. fromAssembly2gene is a perl script run in command line that uses several available programs and R packages to identify genes of your interests in an assembled genome and outputs a descriptive table, the alignment of your gene against the gene in the assembly and the predicted peptide.


It runs several steps:

1. Predicts genes using prodigal https://github.com/hyattpd/Prodigal
2. Finds matches of the genes of interest in the predicted genes in the assembly using local blast.
3. Refines the alignments using an R scripts using the "msa", "reshape2", "Biostrings", "seqinr" packages
4. Using the alignments prints out a table that describes the findings. It generates descriptive tables of presence, absence and truncated genes, curated alignments and peptides predictions

The alternative run plasmids predictions:

```--plasmid -p``` options. This setting works similarly with the exception that before the prediction of genes, it uses "PlasClass" to identify plasmids (https://github.com/Shamir-Lab/PlasClass)
it separates chromosome and plasmids from assembled genomes and find the genes of your interest identifying if they are in the plasmid or in the chromosome.	


## RUN
																						
 From Assemblies to genes (Version  22th June 2021)																																															
 																																	
Detect genes in an assembly, and get their detailed alignments and description using fromAssembly2gene. 		
It generates descriptive tables of presence, absence and truncated  genes, curated alignments and peptides predictions.													
Dependencies: Blast, prodigal and samtools. R packages "msa", "reshape2", "Biostrings", "seqinr".
the program PlasClass is used to predict pasmids.
																																	
run like: 
```
./fromAssembly2gene.pl -g gene/*fasta -a genome/* -o GenesInterest -c 10 												


OPTIONS
	--assemblies  -a put all your files (assemblies) in a folder and write here the path with * at the end.							
	--genes 	  -g put all your files (genes in nucleotide sequences) in a folder and write here the path with *fasta at the end. extension ".fasta".	
	--out 		  -o prefix for output folders																						
	--cores		  -c number of cores to use (4 default).
	--percentage  -i number that indicates the minimum percentage identity considered. defaults is 80%. if your genes are not not highly conserved consider lowering this value
	--plasmid	  -p prediction of plasmids. (if not specify it won't predict plasmids)
	--genes_type  -t type of sequence of the genes input (nucleotypes or aminoacid sequences). default is nucl. if you have protein type: protein 
    --help        -h print this help
    --version     -v version

NOTE 1: When you write * at the end of the path the program will take every file in the folder. 
NOTE 2: The genes must have the extension ".fasta" as it is used as a tag for handling file in the program. 
NOTE 3: The blast indexes are going to be created inside the genes folder. Have this in mind when runing several times, use *fasta to avoid using the index files as inputs.	
NOTE 4: Make sure your input genes have a starting and c codon, and are not in reverse. Unexpected mistakes on the identification of the gene can happen when these are not met. We could account for this but the classification of possibly truncated matches is based on the movement of the stop codon in comparison with the reference.
NOTE 5: If you failed to note 4, the program will add a ATG at the start of your sequence and a TAG at the end.
NOTE 6: If you are completely SURE you dont have plasmids in your assemblies you don't need to set the option and the program will assumed all your assembly is chromosome. 	

```

## OUTPUTS

folder **test_output**:
	This folder contains all the intermedia files used for the program. These files will help you to check in detail where your alignment come from. In case you are puzzle by your final table. Each folder contains:

  **GenomeName.chr.genes.faa**: all predicted genes in AA  
  **GenomeName.chr.genes.fasta**: all predicted genes in nt  
  **GenomeName.chr.genes.gff**: all predicted genes in gff format with the fasta file at the end  
  **GenomeName.chromosome.GeneName.Blast.txt**: blast results of the gene against the genome  
  **GenomeName.chromosome.GeneName.fasta.plusRef.fasta**: fasta of gene reference and gene in the genome  

folder *test_results*  
	This folder contains all final results files and a folder with the predicted peptides that match with the genes of interest  
  **GeneName.fasta.nt_alignment.fasta**: The alignment of each gene of interest for all the genomes analysed  
  **test.alignments.description.txt**: table with the descriptive information of the alignments, stop codons, gaps, insertions, SNPs, N.copies (numbers of copies)  
  **Peptides** (folder): predicted peptides that match with the genes of interest  

Additional files when running any of the --Kleb, --Esch or --Ent options

folder **test_output**:  
	Each assembly has two folders one for the chromosome and one for the plasmid. Examples here are about the chromosome, *.plasmid_* are for plasmids   
   **GenomeName.chromosome_contigslist.txt**:list of contigs in the assembly that are chromosomal  
   **GenomeName.chromosome.fasta**: fasta file of chromosomal contigs  
   **GenomeName.fsa_nt.chromosomesummary.txt**: summary results from chromosome prediction from mlplasmid  
	
