#!/usr/bin/python3.7
"""@Author: Vivek Srinivas - Baliga Lab, ISB
- This is program to get gene expression for GEO samples
- Requires GSM id from GEO dataset
"""

from GEOparse_helpers import *
import time, os, sys, getopt, re
import pandas as pd
from Bio import SeqIO


last_run_gmd = [time.time()]

def get_metadata(gsm_id):
    if time.time() - last_run_gmd[0] < 30:
        time.sleep(35)
        try:
            gsm_metadata = Obtain_metadata.parse_metadata(gsm_id)
        except:
            pass
            print("Metadata was not obtained, check GSM ID for changes in format"%gsm_id)
    else:
        try:
            gsm_metadata = Obtain_metadata.parse_metadata(gsm_id)
        except:
            pass            
            print("Metadata was not obtained, check GSM ID for changes in format"%gsm_id)
    last_run_gmd[0] = time.time()
    return gsm_metadata
        
def get_expression_data(project_name,collection_directory,gene_seq_database, gsm_list, **kwargs):
    """ project_name : Str unique to project
        gene_seq_database : Link to fasta file having transcript sequences and their RefSeq IDs
        gsm_list : list of GSM ids
        ORG_filter [optional] : Boolean, if True will only get sequencing data from samples belonging to mentioned ORG
        ORG [required if ORG_filter = True] : string of organism name. 
        Build_metadata [optional] : Boolean, if False will not build metadata table
        Build_expression_data [optional] : Boolean, if False will not build expression table"""
    ##1 Generate IDs to collect
    ids_to_collect_data = [i.id.capitalize() for i in SeqIO.parse(gene_seq_database,"fasta")]
    ##2 Create gene expression collection directory
    expression_file = "%s/%s_gene_expression.csv"%(collection_directory,project_name)
    if kwargs["Build_expression_data"] == "Y":
        if os.path.exists(expression_file):
            gene_expression_collected_stage1 = pd.read_csv(expression_file,index_col=0)
        else:
            gene_expression_collected_stage1 = pd.DataFrame(index=ids_to_collect_data)
            gene_expression_collected_stage1.to_csv(expression_file)
    elif kwargs["Build_expression_data"] == "N":
        pass
    ##3 Create metadata collection directory
    metadata_file = "%s/%s_metadata.csv"%(collection_directory,project_name)
    if kwargs["Build_metadata"] == "Y":
        if os.path.exists(metadata_file):
            metadata_collected_stage1 = pd.read_csv(metadata_file,index_col=0)
        else:
            md_columns = ["Collected","Exp_type","Exp_ID","Exp_FTP","Pltf_ID","Pltf_FTP","Sample_desc","Exp_desc","Sample_ID","Sample_FTP","Sample_ORG","Sample_acc","Source"]
            metadata_collected_stage1 = pd.DataFrame(columns = md_columns)
            metadata_collected_stage1.to_csv(metadata_file)
    elif kwargs["Build_metadata"] == "N":
        pass
    ##4 Iterate over GSM
    for gsm_id in gsm_list:
        gsm_id = gsm_id.strip()
        print("Initiated data collection for %s"%gsm_id)
        ##4.1 Check if GSM exist in GE collection directory
        if kwargs["Build_expression_data"] == "Y":
            if gsm_id in gene_expression_collected_stage1.columns:
                gsm_cols = gene_expression_collected_stage1[gsm_id].dropna()
                if gsm_cols.size>1000:
                    check_for_GE_collection = False
                    print("GE data for %s exists in catalogue"%gsm_id)
                else:
                    check_for_GE_collection = True
            else:
                check_for_GE_collection = True
        else:
            check_for_GE_collection = False
        ##4.3 Check if GSM exist in MD collection directory
        if kwargs["Build_metadata"] == "Y":
            if gsm_id in metadata_collected_stage1.index:
                check_for_MD_collection = False
                gsm_in_MD_collected = True
                print("Metadata for %s exists in catalogue"%gsm_id)
            else:
                check_for_MD_collection = True
                gsm_in_MD_collected = False
        else:
            check_for_MD_collection = False
            gsm_in_MD_collected = False
        ## 4.4 Collect metadata
        if gsm_in_MD_collected:
            gsm_metadata = metadata_collected_stage1.loc[gsm_id]
            check_for_metadata = True
        else:
            if check_for_GE_collection or check_for_MD_collection:
                try:
                    gsm_metadata = get_metadata(gsm_id)
                    check_for_metadata = True
                except:
                    check_for_metadata = False
            else:
                print("%s will be skipped as GE and MD cataloging is not required"%gsm_id)
                check_for_metadata = False
        ## 4.5 Catalog MD and GE
        if check_for_metadata:
            ## 4.5.1 Catalogue MD
            if check_for_metadata:
                metadata_collected_stagen = pd.read_csv(metadata_file,index_col=0)
                metadata_collected_stagen.loc[gsm_id] = gsm_metadata
                metadata_collected_stagen.to_csv(metadata_file)
            else:
                print("Something went wrong in MD collection\ Metadata for %s was not cataloged"%gsm_id)
                pass
            ## 4.5.2 Catalog GE data
            if check_for_GE_collection:
                gene_expression_collected_stagen = pd.read_csv(expression_file,index_col=0)
                reindex_index = gene_expression_collected_stagen.index
                if gsm_metadata["Exp_type"] == 'Expression profiling by array':
                    pltfid,pltf_ftp,sample_id = gsm_metadata["Pltf_ID"],\
                                                gsm_metadata["Pltf_FTP"],\
                                                gsm_metadata["Sample_ID"]
                    try:
                        expression_data = get_microarray_data.get_expression_arrays(pltfid,pltf_ftp,sample_id)## output has to be pd.series
                        expression_data_reindexed = expression_data.reindex(index=reindex_index)
                        gene_expression_collected_stagen[gsm_id] = expression_data_reindexed
                        gene_expression_collected_stagen.to_csv(expression_file)
                    except:
                        print("Something went wrong in GE data collection\nExpression data for %s was not cataloged"%gsm_id)
                        pass
                elif gsm_metadata["Exp_type"] == 'Expression profiling by high throughput sequencing':
                    try:
                        expression_data = Reanalyze_RNAseq.get_reanalyzed_sequencing_data(gsm_id,gene_seq_database)## output has to be pd.series
                        expression_data_reindexed = expression_data.reindex(index=ids_to_collect_data)
                        gene_expression_collected_stagen[gsm_id] = expression_data_reindexed
                        gene_expression_collected_stagen.to_csv(expression_file)
                    except:
                        print("Something went wrong in GE data collection\nExpression data for %s was not cataloged"%gsm_id)
                        pass
                else:
                    print("%s is not a microarray or a RNAseq sample"%gsm_id)
                    pass
        else:
            print("Matadata for %s was not obtained and GE or MD was not cataloged"%gsm_id)
            pass
        print("-----------######----------------")
            
                
        


def main(argv):
    help_string = "GEOparser.py --help --prjname <> --coldir <> --refseq <> -gsmidfile <> [OPTIONS]"
    help_string2 = "--prjname or -p : <str> Project name under which metadata and gene expression files will be built\n\
--coldir or -c : <str> Collection directory in which metadata and expression files will be stored\n\
--refseq or -ref : <str> Location of reference gene sequences for the organism in fasta format\n\
--gsmidfile or -i : <str> Location of the file (.txt) having GSMIDs required to be collected\n\
OPTIONS:\n\
--Build_metadata or -m: <bool> True or False for metadata collection required; default: True\n\
--Build_expression_data or -e: <bool> True or False for gene expression data collection required; default: True"
    
    project_name,collection_directory,gene_seq_database, gsm_list, Build_meta, Build_exp = "","","",[],"",""
    try:
        opts, args = getopt.getopt(argv,"h:p:c:ref:i:m:e:",["help","prjname=","coldir=","refseq=","gsmidfile=","Build_metadata=","Build_expression_data="])
    except getopt.GetoptError:
        print(help_string+"\n"+help_string2)
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            print(help_string+"\n"+help_string2)
            sys.exit(2)
        elif opt in ("-p","--prjname"):
            project_name = arg
        elif opt in ("-c","--coldir"):
            collection_directory = arg
        elif opt in ("-ref","--refseq"):
            gene_seq_database = arg
        elif opt in ("-i","--gsmidfile"):
            gsm_list = list(set(open(arg,"r").read().split(",")))
        elif opt in ("-m","--Build_metdata"):
            Build_meta = arg
        elif opt in ("-e","--Build_expression_data"):
            Build_exp = arg
    if project_name and collection_directory and gene_seq_database and gsm_list and Build_meta and Build_exp:
        get_expression_data(project_name,collection_directory,gene_seq_database,gsm_list, Build_metadata = Build_meta, Build_expression_data = Build_exp)
    else:
        print("Insufficient arguements\n%s"%help_string)
        sys.exit(2)
        
if __name__ == "__main__":
    main(sys.argv[1:])


