#!/usr/bin/python2.7
"""@Author: Vivek Srinivas - Baliga Lab, ISB
    - A helper function to obtain metadata
    - Requires GSM id from GEO dataset
    """

import os, re, time

# Define path variable

edirect = "/Users/vivek/edirect/"

# create a temproary folder if it doesnt exist

if "temp_data" not in os.listdir(os.getcwd()):
    temp_path_for_metadata_extraction = os.makedirs("temp_data")

# get metadata

def parse_metadata(GSMID):
    run_command = edirect+"esearch -db gds -query "+GSMID+" | "+edirect+"efetch > temp_data/"+GSMID+".txt"
    #print(run_command)
    os.system(run_command)
    text = open("temp_data/%s.txt"%GSMID,"r").read()
    splits = re.split("\n\d\.",text)
    exp_type = re.findall("Type\:\W+([\w\s]+)\s",text)[0]#[0])[0]
    exp_id = re.findall("Accession\:\s?([\w]+)\s?\t?\n?",text)[0]#splits[0])[0]
    exp_ftp = re.findall("ftp\:\/+([\w\.\/]+)GSE",text)[0]+exp_id#splits[0])[0]
    platform_id = re.findall("Accession\:\s?(GPL[\w]+)\s?\t?\n?",text)[0]#splits[-2])[0]
    platform_ftp = re.findall("ftp\:\/+([\w\.\/]+)GPL",text)[0]+platform_id#splits[0])[0]
    sample_id = re.findall("Accession\:\s?(GSM[\w]+)\s?\t?\n?",text)[0]#splits[-1])[0]
    sample_org = re.findall("Organism\:\\t([\W\w]+?)\\n?Type?",text)[0]#splits[-1])[0]
    try:
        exp_desc = re.findall("^\n\d?\.?([\w\W]+?)[\n]?Organism",text)[0]#splits[0])[0]
    except:
        exp_desc = "dnf"
    try:
        sample_desc = re.findall("^^([\w\W]+?)[\n]?Organism",splits[-1])[0]
    except:
        sample_desc = "dnf"
    try:
        sample_ftp = re.findall("ftp\:\/+([\w\.\/]+)GSM",text)[0]+sample_id#splits[-1])[0]
    except:
        sample_ftp = "dnf"
    try:
        RNAseq_acc = re.findall("acc=([\w]+)[\n]?Sample",text)[0]#splits[-1])[0]
    except:
        RNAseq_acc = "dnf"
    try:
        source = re.findall("Source\sname\:\t?([\W\w]+?)\n?Platform",splits[-1])[0]
    except:
        source = "dnf"
    collected = "Y"
    os.system("rm temp_data/%s.txt"%GSMID)
    return {k:v for k,v in zip(["Collected","Exp_type","Exp_ID","Exp_FTP","Pltf_ID","Pltf_FTP","Sample_desc","Exp_desc","Sample_ID","Sample_FTP","Sample_ORG","Sample_acc","Source"],\
                               [collected,exp_type,exp_id,exp_ftp,platform_id,platform_ftp,sample_desc,exp_desc,sample_id,sample_ftp,sample_org,RNAseq_acc,source])}





