#importing modules
import os
import requests
from bs4 import BeautifulSoup
import pandas as pd
import urllib.request
import argparse
import gzip
import shutil
import subprocess

##PARSER##

parser = argparse.ArgumentParser(description="Downloading all the available species from the WormBase database, indexed for using the download function on WormQuest.")

parser.add_argument('species', nargs="?", default = None, help = "Involved with the available WormBase species. Defaults to help message.")
parser.add_argument('-m', '--mine', action='store_true', help="Allows to mine the WormBase database for species names and their respective Bioproject code.")
parser.add_argument('-dl', '--download', nargs="+", type=int, help="Downloads the specified species from the Wormbase database using indexes.")

parser.add_argument('pfam', nargs="?", help="Run jobs related to PFam files.")
parser.add_argument('-dl_hmm', '--download_hmm', action='store_true', help="Download the associated HMM annotation for the Pfam identifier.")

parser.add_argument('hpc', nargs="?", help="Commands with ALICE HPC System.")
parser.add_argument('-tr', '--transfer', action='store_true', help="Transfers generated directoriies to ALICE.")
parser.add_argument('-hmms', '--hmmsearch', action='store_true', help="Generates a HMMer hmmsearch script to run on ALICE.")

parser.add_argument('analysis', nargs="?", help="WormQuest analysis tools.")
parser.add_argument('-d', '--demo', action='store_true', help="Runs through a WormQuest HMMer analysis demo.")

args = parser.parse_args()

##FUNCTION 1: DOWNLOADS AND UNPACKING##

##SECTION A: RETRIEVE FASTA.GZ FILES##

#wormquest.py species -m = allows the user to get the list of species to help choose.

if args.species is None:
    print('''
          WormQuest 1.0.0
          'species' allows you to interact with Wormbase species data in different ways:
          >>> -m or --mine: Mines the WormBase database for all the available species and their respective BioRender. 
                            It will generate a .csv file needed to download files from WormBase.
          >>> -dl or --download "species":
                            Downloads the specified species from the WormBase database. 
                            It requires the user to use species indexes from the .csv file.
          ''')
    exit()
if args.mine:    
    #Using BeautifulSoup, the script scrapes the html content of the wormbase main page for species, and adds them to a .csv for user reference:

    #mining database
    url = "https://parasite.wormbase.org/ftp.html"
    response = requests.get(url)
    soup = BeautifulSoup(response.content, 'html.parser')
    species_list = []

    #populating dataframe
    for row in soup.find_all('tr')[10:]: #skips all the detail before the species table
        columns = row.find_all('td')
        if len(columns) > 0:
            species = columns[0].text.strip()
            bioproject = columns[1].text.strip()
            species_list.append({"Species": species, "BioProject": bioproject})

    species_data = pd.DataFrame(species_list)
    species_data.to_csv('species_data.csv', index = True)
    print(f"'Species_data.csv' has been created in your working directory.")
    
    exit()

#wormquest.py species -dl [speciesindex1] speciesindex2] [speciesindex3] = downloads 3 species specified by the user

#Reading csv
wb_dir = "wormbase"

if args.download:
    try:
        species_data = pd.read_csv('species_data.csv', index_col=0)
    except FileNotFoundError:
        print("Error: 'species_data.csv' not found in your working directory.\nPlease run -m or --mine first to retrieve the data.")
        exit()

    ##Create a directory in the current working directory for FASTA files

    os.makedirs(wb_dir, exist_ok=True)
    os.chdir(wb_dir)

#formatting function

    def species_name_formatting(species):
        species= species.lower()
        species = species.replace("_","")
        species= species.replace(" sp.", "") 
        species= species.replace(" ", "_")
        species= species.replace(".", "") 
        return species
    

#download function

    def wormbase_download_file(species_name,bioproject,base_url,max_retries=2):
        for attempt in range(max_retries + 1):
            prefix = f"TD{attempt + 1}_" if attempt > 0 else ""
            url = f"{base_url}/{species_name}/{prefix}{bioproject}/{species_name}.{prefix}{bioproject}.WBPS19.protein.fa.gz"
            filename = f"{species_name}.{prefix}{bioproject}.WBPS19.protein.fa.gz"

            try:
                print(f"Trying URL: {url}")
                urllib.request.urlretrieve(url, filename) #Code from https://www.squash.io/how-to-download-a-file-over-http-in-python/
                print(f"Successfully downloaded: {filename}")
                return True  # Successful download
            except Exception as e:
                print(f"Failed to download with URL: {url} - {e}")    

    wb_base_url = "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species"

    #Retrieving data from csv
    for index in args.download:
        row = species_data.loc[index]
        raw_species_name = row["Species"]
        species_name = species_name_formatting(raw_species_name)
        bioproject = row["BioProject"]

#Download FASTA from database
    #Generating download link
        wb_download = wormbase_download_file(species_name, bioproject, wb_base_url)
        if not wb_download:
            print(f"Could not download {raw_species_name} with BioProject {bioproject}.")

    print(f"Your files have been successfully downloaded.\nAvailable in the {wb_dir} directory.")

##SECTION B: GATHER AND UNZIP HMM FILES##

#wormquest.py pfam -dl_hmm = extracts all PFam identifiers from .tsv file given, and downloads them into a new directory.
pf_dir = "PFam_HMM"

if args.download_hmm:
    #reading tsv
    SD_matches = pd.read_csv('SearchResults-succinatedehydrogenase.tsv', sep="\t", header = 0)

    #Find PFam identifiers
    PF_identifiers = SD_matches[SD_matches['Accession'].str.startswith('PF', na=False)]['Accession'].tolist()

    #Downloading HMMs from PFam

    ##Create a directory in the current working directory for HMM files
    os.makedirs(pf_dir, exist_ok=True)
    os.chdir(pf_dir)

    #defining download function

    def pfam_download_file(identifier,base_url):
        url = f"{base_url}{identifier}?annotation=hmm"
        filename = f"{identifier}.gz"

        try:
            print(f"Trying URL: {url}")
            urllib.request.urlretrieve(url, filename)
            print(f"Successfully downloaded: {filename}")
            return True  # Successful download
        except Exception as e:
            print(f"Failed to download with URL: {url} - {e}")  

    pfam_base_url = f"https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/"

#Downloads the HMM files associated to mined PFam identifiers
#Also gunzips all files; sourced from: https://stackoverflow.com/questions/31028815/how-to-unzip-gz-file-using-python

    for pf_id in PF_identifiers:
        hmm_download = pfam_download_file(pf_id, pfam_base_url)
        with gzip.open(f'{pf_id}.gz', 'rb') as f_in:
            with open(f'{pf_id}.txt', 'wb') as f_out:
                shutil.copyfileobj(f_in,f_out)
        print("The HMM annotations have been successfully unzipped.")  
        if not hmm_download:
            print(f"Could not download {pf_id} identifier from PFam database.")
    print(f"Total number of files = {len(PF_identifiers)}.\nAvailable in the {pf_dir} directory.")

    ##FUNCTION 2: GENERATE SLURP SUBMISSION SCRIPT FOR HMMSEARCH##
def transfer_directory(dir, username):
        try:
            # Construct the SCP command
            command = ["scp", "-r", dir, f"{username}@alice.le.ac.uk:~"]
            
            # Execute the command
            subprocess.run(command, check=True)
            print(f"Successfully transferred {dir} to {username}@alice.le.ac.uk:~")
        except subprocess.CalledProcessError as e:
            print(f"Error during transfer: {e}")

alice_username = ()
           
if args.transfer:


    def transfer_directory(dir, username):
        try:
            # Construct the SCP command
            command = ["scp", "-r", dir, f"{username}@alice.le.ac.uk:~"]
            
            # Execute the command
            subprocess.run(command, check=True)
            print(f"Successfully transferred {dir} to {username}@alice.le.ac.uk:~")
        except subprocess.CalledProcessError as e:
            print(f"Error during transfer: {e}")
   
    alice_username = input("Input your ALICE username (ex. 'ab123'): ").lower()
    transfer_directory(wb_dir, alice_username)
    transfer_directory(pf_dir, alice_username)

if args.hmmsearch:

    # Prompt the user for input
    username = input("Enter your username: ")
    user_mail = input("Enter your email address: ")

    # Create the username_initial (first letter of the username)
    username_initial = username[0].lower()

    # Define directories
    pfam_hmm_dir = './PFam_HMM'  # Ensure this is the correct path
    wormbase_dir = './wormbase'  # Ensure this is the correct path

    # Get all *.txt files in PFam_HMM directory
    try:
        txt_files = [os.path.splitext(f)[0] for f in os.listdir(pfam_hmm_dir) if f.endswith('.txt')]
    except FileNotFoundError:
        print(f"Error: The directory {pfam_hmm_dir} does not exist.")
        txt_files = []

    # Get all files in the wormbase directory
    try:
        dbase_files = [f for f in os.listdir(wormbase_dir) if f.endswith('.fa.gz')]
    except FileNotFoundError:
        print(f"Error: The directory {wormbase_dir} does not exist.")
        dbase_files = []

    # SLURM script template
    slurm_script = f"""#!/usr/bin/bash
#SBATCH --job-name=hmmsearch_{username}
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=8gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user={user_mail}
#SBATCH --export=NONE

# setting files/directories
out=/home/{username_initial}/{username}/
hmm=/home/{username_initial}/{username}/PFam_HMM
dbase=/home/{username_initial}/{username}/wormbase

# executable from ALICE
hmmsearch=/cm/shared/spack/opt/spack/linux-rocky9-x86_64_v3/gcc-12.3.0/hmmer-3.3.2-6xbbubuci3m4xefi3yqu6thg5loc4jju/bin/hmmsearch

# module needed for using HPC installed software
module load gcc/12.3.0-yxgv2bl
module load openmpi/4.1.5-fzc7xdf
module load hmmer/3.3.2-6xbbubu

# Execute the job code
"""

    # Check that both txt_files and dbase_files have content
    if txt_files and dbase_files:
        # Generate the commands for each .txt and .dbase file pair
        for txt_file in txt_files:
            for dbase_file in dbase_files:
                dbase_name = dbase_file.split('.')[0]
                slurm_script += f"$hmmsearch --tblout $out/{txt_file}_{dbase_name}.out -E 0.1 --noali $hmm/{txt_file}.txt $dbase/{dbase_file}\n"
    else:
        print("No files found in the required directories.")
        exit(1)

    # Save the SLURM script to a file
    slurm_filename = f"hmmsearch_{username}.sh"
    with open(slurm_filename, "w") as slurm_file:
        slurm_file.write(slurm_script)

    print(f"SLURM script for {username} has been generated and saved as {slurm_filename}")
    alice_transfer = input("Would you like to transfer the script to your ALICE working directory? [y/n]: ")
    if alice_transfer.lower() == 'y':

        def transfer_file(file, username):
            try:
                # Construct the SCP command
                command = ["scp",file, f"{username}@alice.le.ac.uk:~"]
                
                # Execute the command
                subprocess.run(command, check=True)
                print(f"Successfully transferred {file} to {username}@alice.le.ac.uk:~")
            except subprocess.CalledProcessError as e:
                print(f"Error during transfer: {e}")
        
        alice_username = input("Input your ALICE username (ex. 'ab123'): ").lower()
        transfer_file(slurm_filename, alice_username)

    else:
        exit()

###FUNCTION 3: PARSE AND ANALYSE HMMER FILES###

#wormquest.py analysis -d = demos the parsing and analysis of HMMER files

if args.demo:

    print(f'''Welcome to the WormQuest demo!

This demo will show how WormQuest can parse and analysise HMMsearch output files.

WormQuest can create scripts for HMMer searches using the HPC ALICE system (refer to [hpc -hmms] tool).
Within the WormQuest GitHub Repo, you will be able to access the 'demo' directory. Please ensure you have acquired this in your working directory before proceeding.
    ''')

    dir_check = input("Do you have the 'demo' directory within your current working directory? [y/n]: ")
    if dir_check.lower() == "y":
        print(f'''The 'demo' directory contains 18 HMMsearch output files using six HMM PFAM files and three nematode species from the WormBase database.
WormQuest is able to parse through them, and generate a variety of tables, figures and summaries.
        ''')
        print('''Please select an option:
[a] Graphs
[b] Tables
[c] Summaries
[d] Exit''')




    #Parse through the files (function?)

        exit()