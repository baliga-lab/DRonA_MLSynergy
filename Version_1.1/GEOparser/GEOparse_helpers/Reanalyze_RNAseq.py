# -*- coding: utf-8 -*-
"""@Author: Vivek Srinivas - Baliga Lab, ISB
    - A helper function to reanalyze RNA GSM samples present on GEO 
    - Requires GSM id from GEO dataset
    """

import os, sys, glob, time, datetime, paramiko, getpass, os.path
import subprocess as sp
from os import path
import pandas as pd

# Create temporary data folder 
if "temp_data" not in os.listdir(os.getcwd()):
    temp_path_for_metadata_extraction = os.makedirs("temp_data")
    
# Define path variable
edirect = "/Users/vivek/edirect/"
fastq_dump_path = "/usr/local/Cellar/sratoolkit/2.9.2/bin/"
aspera_client_path = "/Users/vivek/Applications/Aspera_Connect.app/Contents/Resources/"
aspera_private_key_path = aspera_client_path+"asperaweb_id_dsa.openssh "
kallisto_path = ""

# functions to get read info from sra db


def run_get_info(GSM_ID):
    t0 = time.time()
    run_command = edirect+"esearch -db SRA -query "+GSM_ID+" | "+edirect+"efetch -format runinfo >temp_data/"+GSM_ID+"_SRA_data.txt"
    print(run_command)
    os.system(run_command)
    text = open("temp_data/%s_SRA_data.txt"%GSM_ID,"r").readlines()
    os.system("rm temp_data/%s_SRA_data.txt"%GSM_ID)
    return {k:v for k,v in zip(text[0].split(","),text[1].split(","))}

lastrun = [time.time()]
def get_info(GSM_ID):
    t0 = time.time()
    if (t0 - lastrun[0]) < 60:
        time.sleep(61)
        r = run_get_info(GSM_ID)
    else:
        r = run_get_info(GSM_ID)
    if r:
        pass
    else:
        print("Run info for %s was not collected, check GSM ID or internet connection"%GSM_ID)
    lastrun[0] = time.time()
    return r

# functions to download fastq files of the sra

def run_aspera(srr_number):
    command = aspera_client_path+\
              "ascp -T -k 1 -i "+\
              aspera_private_key_path+\
              "anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/"+\
              "%s/%s/%s/%s.sra "%(srr_number[:3],srr_number[:6],srr_number,srr_number)+\
              " temp_data/%s.sra"%(srr_number)
    print(command)
    os.system(command)

def dump_fastq(srr_number):
    command = fastq_dump_path+"fastq-dump temp_data/%s.sra -O temp_data/%s --split-files"%(srr_number,srr_number)
    print(command)
    os.system(command)
    os.system("rm temp_data/%s.sra"%srr_number)

def download_fastq(srr_number):
    if path.isdir("temp_data/%s"%srr_number):
        if len([i for i in os.listdir("temp_data/%s"%srr_number) if i.split(".")[-1] == "fastq"]) > 0:
            pass
    else:
        t0 = time.time()
        run_aspera(srr_number) #step 1
        dump_fastq(srr_number)
        tn = time.time()
        print("FastQ downloaded in %0.2f secs"%(tn-t0))
    if os.path.exists("temp_data/%s"%srr_number):
        fastq_existence = True
    else:
        fastq_existence = False
    return fastq_existence

# Functions to reanalyze RNAseq
def create_kallisto_index(fasta_file):
    parent_dir = "/".join(fasta_file.split("/")[:-1])
    fasta = fasta_file.split("/")[-1]
    index = parent_dir+"/"+fasta.split(".")[0]+"_index"
    if path.isfile(index):
        pass
    else:
        command = "kallisto index -i %s %s"%(index,fasta_file)
        os.system(command)
    return index
    
def reanalyze_with_kallisto(srr_number,srr_info,indexed_database):
    t0 = time.time()
    if path.isfile("temp_data/%s/abundance.tsv"%srr_number) is False:
        if srr_info["LibraryLayout"].lower() == "single":
            average = srr_info["avgLength"]
            deviation = int(srr_info["InsertDev"])+1
            fastq_files = " ".join(["temp_data/%s/%s"%(srr_number,i) for i in os.listdir("temp_data/%s"%srr_number) if i.split(".")[-1]=="fastq"])
            command = "kallisto quant -i %s -o temp_data/%s --single -l %s -s %s %s"%(indexed_database,srr_number,average,deviation,fastq_files)
        elif srr_info["LibraryLayout"].lower() == "paired":
            fastq_files = " ".join(["temp_data/%s/%s"%(srr_number,i) for i in os.listdir("temp_data/%s"%srr_number) if i.split(".")[-1]=="fastq"])
            command = "kallisto quant -i %s -o temp_data/%s %s"%(indexed_database,srr_number,fastq_files)
        os.system(command)
        tn = time.time()
        print("The sample was reanalyzed in %0.2f secs"%(tn-t0))
        if path.isfile("temp_data/%s/abundance.tsv"%srr_number):
            return_dict = pd.read_csv("temp_data/%s/abundance.tsv"%srr_number, sep = "\t", index_col = 0)
        else:
            return_dict = pd.DataFrame()
    else:
        return_dict = pd.read_csv("temp_data/%s/abundance.tsv"%srr_number, sep = "\t", index_col = 0)
    os.system("rm -r temp_data/%s"%srr_number)
    return return_dict


# Compiler function
def get_reanalyzed_sequencing_data(GSM_ID, gene_seq_database):
    run_info = get_info(GSM_ID)
    srr_number = run_info["Run"]
    fastq_existence = download_fastq(srr_number)
    if fastq_existence:
        indexed_database = create_kallisto_index(gene_seq_database)
        ret = reanalyze_with_kallisto(srr_number,run_info,indexed_database)
        ids = [i.capitalize() for i in ret.index]
        try:
            ret_data = pd.series(data = ret.est_counts.values, index = ids)
        except:
            pass
    return ret_data


