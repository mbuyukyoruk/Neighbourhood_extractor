import argparse
import sys
import os
import subprocess
import re
import textwrap
import uuid
import time

timestr = time.strftime("%Y%m%d_%H%M%S")

try:
    from Bio import SeqIO
except:
    print("SeqIO module is not installed! Please install SeqIO and try again.")
    sys.exit()

try:
    import tqdm
except:
    print("tqdm module is not installed! Please install tqdm and try again.")
    sys.exit()

parser = argparse.ArgumentParser(prog='python Neighbourhood_extractor.py',
      formatter_class=argparse.RawDescriptionHelpFormatter,
      epilog=textwrap.dedent('''\

# Neighbourhood_extractor

Author: Murat Buyukyoruk

        Neighbourhood_extractor help:

This script is developed to fetch flanking ORFs from a list of ORFs to perform neighbourhood analysis. The extracted ORFs can be submited to EggNOG-mapper to get COG gene IDs and parsed with eggnogg_mapper.py to use it for Phylo_2_genemap.py script. 

SeqIO package from Bio is required to fetch sequences. IMPORTANT!!! The SeqIO update is required to hangel PRODIGAL ORFs! See BioSeqIO_Update repositories to setup any related updates.

Additionally, tqdm is required to provide a progress bar since some multifasta files can contain long and many sequences.
        
### Example list file format:        

Following bash script can be used to generate following list from an PRODIGAL ORF .fasta file:

getfastalist ()
{ 
    echo -e "Accession \t ORF \t start \t end \t strand" > "${1%.*}"_header_list.txt;
    grep --color=auto ">" "$1" | cut -d">" -f2 | cut -d' ' -f1 | rev | cut -d'_' -f2- | rev > file_1.txt;
    grep --color=auto ">" "$1" | cut -d">" -f2 | cut -d' ' -f1,3,5,7 --output-delimiter '	' > file_2.txt;
    paste file_1.txt file_2.txt >> "${1%.*}"_header_list.txt;
    rm file_1.txt file_2.txt
}

getfastalist PRODIGAL_ORF.fasta

File:
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

	-o/--output		output file	    Specify a output file name that should contain fetched sequences.

Parameters [optional]:
----------------------
	-f/--flank		5000			This is the default length of flanks that is fetched.

Basic Options:
--------------
	-h/--help		HELP			Shows this help text and exits the run.

      	'''))
parser.add_argument('-i', '--input', required=True, type=str, dest='filename',
                    help='Specify a original fasta file.\n')
parser.add_argument('-l', '--list', required=False, type=str, dest='list_file', default= None,
                    help='Specify a list of accession numbers to fetch.\n')
parser.add_argument('-o', '--output', required=False, dest='out', default = "out_" + timestr + ".fasta",
                    help='Specify a output fasta file name.\n')
parser.add_argument('-f', '--flank', required=False, type=str, dest='flank', default=5000,
                    help='Specify flank i.e., 5000.\n')

results = parser.parse_args()
filename = results.filename
list_file = results.list_file
out = results.out
flank = results.flank

dict = {}
dict_ORF = {}
all = None

try:
    flank = int(flank)
except:
    pass

if isinstance(flank, str):
    if flank.lower() in ["all","a"]:
        all =True
        flank_decide = 0
    else:
        print("Invalid option for flank. Use number (i.e., 1000) or 'all','All','a','A' to grab all ORFs from rlated genomes.")
        sys.exit()
elif isinstance(flank, int):
    all = False
    flank_decide = flank
else:
    print("Invalid option for flank. Use number (i.e., 1000) or 'all','All','a','A' to grab all ORFs from rlated genomes.")
    sys.exit()

flank = flank_decide

print("Locating ORF index file! An index file will be generated if the fasta file is provided for the first time to generate DB.")

index = filename + '.ORF.index'
record_dict = SeqIO.ORF_to_genome_index_db(index,filename,"fasta",key_function=lambda rec : rec.rsplit("_",1)[0])

if list_file != None:
    proc = subprocess.Popen("wc -l < " + list_file, shell=True, stdout=subprocess.PIPE, text=True)
    length = int(proc.communicate()[0].split('\n')[0])
    uniq_id = uuid.uuid4()
    tmp = "/tmp/" + str(uniq_id)
    f = open(tmp, 'a')
    sys.stdout = f
    with tqdm.tqdm(range(length)) as pbar:
        pbar.set_description('Reading accession list...')
        with open(list_file, 'r') as file:
            for line in file:
                pbar.update()
                if "Accession" not in line:
                    if len(line.split()) != 0:
                        arr = line.replace("\n","").split("\t")
                        genome = arr[0]
                        ORF_target = arr[1]
                        start = int(arr[2])
                        stop = int(arr[3])
                        strand = arr[4]
                        dict_ORF[ORF_target] = range(start - flank - 1, stop + flank + 1)
                        if genome not in dict:
                            dict[genome] = [ORF_target]
                            record_ask = record_dict[genome]
                            # record_ask = record_dict[ORF_target]
                            if isinstance(record_ask, list):
                                for i in range(len(record_ask)):
                                    print(record_ask[i].format("fasta"))
                            else:
                                print(record_ask.format("fasta"))
                        else:
                            dict[genome].append(ORF_target)

    if all == True:
        os.system("mv " + tmp + " " + out)
        sys.exit()
    elif all == False:
        os.system('> ' + out)
        f = open(out, 'a')
        sys.stdout = f
        proc = subprocess.Popen("grep -c '>' " + tmp, shell=True, stdout=subprocess.PIPE, text=True)
        length = int(proc.communicate()[0].split('\n')[0])
        with tqdm.tqdm(range(length)) as pbar:
            pbar.set_description('Grabbing...')
            for record in SeqIO.parse(tmp, "fasta"):
                pbar.update()
                genome_ask = record.id.rsplit('_', 1)[0]
                if genome_ask in dict:
                    start_ask = int(record.description.split(" # ")[1])
                    stop_ask = int(record.description.split(" # ")[2])

                    for i in range(len(dict[genome_ask])):
                        if set(range(start_ask, stop_ask)).issubset(dict_ORF[dict[genome_ask][i]]):
                            print(">" + record.id + "|" + dict[genome_ask][i] + record.description.replace(record.id, ""))
                            print(re.sub("(.{60})", "\\1\n", str(record.seq), 0, re.DOTALL))

        time.sleep(0.5)
        os.system('rm ' + tmp)
    else:
        print("Something is wrong! Exiting")
        sys.exit()

