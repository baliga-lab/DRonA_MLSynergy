#!/usr/bin/python3.0
"""@Author: Vivek Srinivas - Baliga Lab, ISB
- Data preparation support for DRonA and MLSynergy
"""

# Import functions

import pandas as pd
import numpy as np
from scipy.stats import rankdata
import sys

# Functions

def reindex_data(data,index_file):
    index = [i.strip() for i in open(index_file,"r").readlines()]
    if len(index) < 1:
        print("Provide location of ID file in .txt format and IDs in new line")
        sys.exit(2)
    return data.reindex(index)

def clean_data(data, threshold):
    """The function takes two arguements:
        data: pandas dataframe with genes in rows and GSM samples in columns
        threshold: float (0 to 1)"""
    copy = data.copy()
    for columns in data.columns:
        useless_vals = data[columns].isin([np.Inf,-np.Inf,np.NAN]).sum()
        if (useless_vals/float(data[columns].size)) > threshold:
            copy = copy.drop(columns = columns)
        else:
            pass
    print("Data of %s was reduced to data of %s"%(str(data.shape),str(copy.shape)))
    return copy

def rank_normalize_across_genes(data, method):
    normalized_data = data.copy()
    for columns in data.columns:
        normalized_data[columns]= rankdata(data[columns].values,method=method)
    return normalized_data

