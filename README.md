# Neighbourhood_extractor

Author: Murat Buyukyoruk

        Neighbourhood_extractor help:

This script is developed to fetch flanking ORFs from a list of ORFs to perform neighbourhood analysis. The extracted ORFs can be submited to EggNOG-mapper to get COG gene IDs and parsed with eggnogg_mapper.py to use it for Phylo_2_genemap.py script. 

SeqIO package from Bio is required to fetch sequences. IMPORTANT!!! The SeqIO update is required to hangel PRODIGAL ORFs! See BioSeqIO_Update repositories to setup any related updates.

Additionally, tqdm is required to provide a progress bar since some multifasta files can contain long and many sequences.
        
### Example list file format:        

Following bash script can be used to generate following list from an PRODIGAL ORF .fasta file:

bash function:

    getfastalist ()
    { 
        echo -e "Accession \t ORF \t start \t end \t strand" > "${1%.*}"_header_list.txt;
        grep --color=auto ">" "$1" | cut -d">" -f2 | cut -d' ' -f1 | rev | cut -d'_' -f2- | rev > file_1.txt;
        grep --color=auto ">" "$1" | cut -d">" -f2 | cut -d' ' -f1,3,5,7 --output-delimiter '	' > file_2.txt;
        paste file_1.txt file_2.txt >> "${1%.*}"_header_list.txt;
        rm file_1.txt file_2.txt
    }

Execute function:

    getfastalist PRODIGAL_ORF.fasta

Accession file output:

    Accession   Protein_accession/array_no  Start   Stop    Strand
    NZ_CP006019 NZ_CP006019_1756            1875203 1877050 -1
    CP000472.1  CP000472.1_235              123     975     1
        
Syntax:

    python Neighbourhood_extractor.py -i demo.fasta -l demo_acc_list.txt -o demo_flanks_3kb.fasta -f 3000

OR
        
    python Neighbourhood_extractor.py -i demo.fasta -l demo_acc_list.txt -o demo_flanks_all.fasta -f all

Neighbourhood_extractor dependencies:

Bio module and SeqIO available in this package      refer to https://biopython.org/wiki/Download

Specific Bio update is required                     refer to https://github.com/mbuyukyoruk?tab=repositories

tqdm                                                refer to https://pypi.org/project/tqdm/
	
Input Paramaters (REQUIRED):
----------------------------
	-i/--input		FASTA			Specify a PRODIGAL fasta file.)

	-l/--list		List			Specify a list of PRODIGAL ORFs and their position specific intormations.

	-o/--output		output file             Specify a output file name that should contain fetched sequences.

Parameters [optional]:
----------------------
	-f/--flank		5000			This is the default length of flanks that is fetched.

Basic Options:
--------------
	-h/--help		HELP			Shows this help text and exits the run.
