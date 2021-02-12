# -*- coding: utf-8 -*-
"""@Author: Vivek Srinivas - Baliga Lab, ISB
- This is program to get microarray based gene expression for GEO samples
- Requires GSM id, GPL id and FTP of GPL
"""


import os,re, shutil, csv, pathlib
import urllib2 as urllib3
import pandas as pd
import numpy as np
from contextlib import closing
from ftplib import FTP
import xml.etree.ElementTree as ET
from collections import OrderedDict as OD

    
# Create temporary data folder 
if "temp_data" not in os.listdir(os.getcwd()):
    temp_path_for_metadata_extraction = os.makedirs("temp_data")
else:
    pass

# functions for downloading microarray gene expression data

def download_ex_genex_array_for_platforms(GSP_FTP, GSP):
    with closing(urllib3.urlopen("ftp://"+GSP_FTP+'miniml/%s_family.xml.tgz'%GSP)) as r:
        with open('temp_data/%s.tgz'%GSP, 'wb') as f:
            shutil.copyfileobj(r, f)
    if GSP not in os.listdir("temp_data/"):
        os.mkdir("temp_data/%s"%GSP)
    os.system("tar -xvzf %s -C %s"%('temp_data/%s.tgz'%GSP,'temp_data/%s'%GSP))
    os.system("rm %s"%('temp_data/%s.tgz'%GSP))


def get_platform_table_values(platform_xml):
    tree = ET.parse(platform_xml)
    root = tree.getroot()
    ret_dict = {}
    for elements in root:
        if elements.tag.split("}")[-1]=="Platform":
            for datatable in elements.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Data-Table"):
                for cols in datatable.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Column"):
                    ret_dict.update({cols.find("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Name").text.capitalize():\
                                     cols.attrib["position"]})
    return ret_dict

def get_sample_table_values(platform_xml):
    tree = ET.parse(platform_xml)
    root = tree.getroot()
    ret_dict = {}
    for elements in root:
        if elements.tag.split("}")[-1]=="Sample":
            for datatable in elements.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Data-Table"):
                for cols in datatable.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Column"):
                    ret_dict.update({cols.find("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Name").text.capitalize():\
                                     cols.attrib["position"]})
            break
    return ret_dict

def create_dict_for_platforms(gpl_id):
    gpl_xml = "temp_data/%s/%s_family.xml"%(gpl_id,gpl_id)
    gpl_table = "temp_data/%s/%s-tbl-1.txt"%(gpl_id,gpl_id)
    print(gpl_xml)
    print(gpl_table)
    platform_details = get_platform_table_values(gpl_xml)
    Id = int(platform_details["Id"])-1
    for i in ["Orf","Gene_symbol","Accession_string","Gene"]:
        if i in platform_details.keys():
            Gene_names = int(platform_details[i])-1
            break
    #Gene_names = 3 # Overriding options
    print ("%s-%s"%(Gene_names, Id))
    ret_dict = pd.DataFrame({},columns = ["Positions"])
    for line in open(gpl_table,"r").readlines():
        split = line.split("\t")
        gene_id = split[Gene_names].lower()
        #gene_id = split[Gene_names].lower().split("=")[-1] # Overriding options
        if gene_id[:2] in ["mt","rv"]:
        #if split[Gene_names].lower()[:2] in ["mt","rv"]:
            if gene_id not in ret_dict.index:
                ret_dict.at[gene_id,"Positions"]=[split[Id]]
            else:
                ret_dict.loc[gene_id].Positions.append(split[Id])
        else:
            pass
    if ret_dict.shape[0]> 100:
        ret_dict.to_pickle("temp_data/Pickles/%s.pkl"%gpl_id)
        print("The indexes of the recent written %s platform dictionary are %s"%(gpl_id,str(ret_dict.index[:5])))
        status="s"
    else:
        print("Dictionary was not made")
        status="f"
        

def create_dict_for_GSM(fname,spot,value):
    ret_dict = OD()
    with open(fname,"r") as f:
        csvreader = csv.reader(f,delimiter = "\t")
        for i in csvreader:
            try:
                ret_dict[i[0]] = float(i[1])
            except:
                pass
    return ret_dict

def match_GPL_to_GSM(GSM,GPL,GPL_dict):
    gpl_xml = "temp_data/%s/%s_family.xml"%(GPL,GPL)
    gpl_table = "temp_data/%s/%s-tbl-1.txt"%(GPL,GPL)
    sample_details = get_sample_table_values(gpl_xml)
    spot = int(sample_details["Id_ref"])-1
    for i in ["Value"]:
        if i in sample_details.keys():
            value = int(sample_details[i])-1
            break
    gsm = create_dict_for_GSM("temp_data/%s/%s-tbl-1.txt"%(GPL,GSM),spot,value)
    ret_series = pd.Series()
    MTB_id = re.compile("[RrVvMmTt]{2}\d+[\w\.\d]?[^\,\:]$")
    for k,v in GPL_dict.iterrows():
        try:
            ret_series.loc[k.capitalize()] = np.median([gsm[i] for i in v.Positions])
        except:
            pass
    return ret_series
   
def get_expression_arrays(gpl_id,gpl_ftp,gsm_id):
    """Requires gpl_id, gpl_ftp and gsm_ids.
    The function will create a csv file in project name and write data to it.
    Provide same project name if the data should be written in one file"""
    if gpl_id not in os.listdir("temp_data/"):
            download_ex_genex_array_for_platforms(gpl_ftp,gpl_id)
    else:
        pass
    if "%s-tbl-1.txt"%gpl_id in os.listdir("temp_data/%s/"%gpl_id):
        if "%s.pkl"%gpl_id not in os.listdir("temp_data/Pickles/"):
            print("Preparing platform dictionaries")
            create_dict_for_platforms(gpl_id)
        else:
            pass
    if "%s.pkl"%gpl_id in os.listdir("temp_data/Pickles/"):
        platform_dict = pd.read_pickle("temp_data/Pickles/%s.pkl"%gpl_id)
        try:
                matched_values = match_GPL_to_GSM(gsm_id,gpl_id,platform_dict)
        except:
                pass
    return matched_values





