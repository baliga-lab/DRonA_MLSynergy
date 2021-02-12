#!/usr/bin/python3.7
"""@Author: Vivek Srinivas - Baliga Lab, ISB, Seattle, WA, USA
- MLSynergy
"""


#1. Imports
import os, sys, pickle, getopt, Data_prep, itertools
import pandas as pd
import numpy as np
from sklearn.svm import OneClassSVM as OCSVM
from scipy import stats
## Chord diagram
import holoviews as hv
from holoviews import opts, dim
from bokeh.sampledata.les_mis import data
from bokeh.plotting import show, output_file

#2. Functions

def give_combinations(listi,n,joinstr,**kwargs):
    ret_list = []
    if "drug_repeat" in kwargs.keys():
        val = kwargs["drug_repeat"]
    else:
        val = "F"
    for i in itertools.combinations(listi,n):
        len_tr = len(list(set(z.split("_")[2] for z in i)))
        if len_tr ==1:
            len_d = len(list(set(z.split("_")[0] for z in i)))
            if val == "F":
                if len_d == n:
                    jj = list(i)
                    jj.sort()
                    ret_list.append(joinstr.join(jj))
            elif val == "T":
                jj = list(i)
                jj.sort()
                ret_list.append(joinstr.join(jj))   
    return ret_list

def prepare_combinations_file(metadata_loc,N,tr_times):
    metadata = pd.read_csv(metadata_loc,index_col=0)
    d_c_t_ps = ["%s_%0.2f_%0.0f"%names for names, group in \
             metadata.groupby(["Drug","Concentration","Treatment time (hrs)"])\
                if float(names[2]) >= min(tr_times) and float(names[2])<=max(tr_times) and names[0] not in ["No drug","Untreated"]]
    d_c_t = d_c_t_ps
    drug_comb = []
    for nbd in N:
        drug_comb.extend(give_combinations(d_c_t,nbd,"#"))
    ret_df = pd.DataFrame()
    for combs in drug_comb:
        ds = []
        for nz, dct in enumerate(combs.split("#")):
            nn = nz+1.
            d,c,t = dct.split("_")
            ret_df.at[combs,"Drug%0.0f"%nn] = d
            ret_df.at[combs,"D%0.0f_conc"%nn] = c
            ds.append(d)
        ds.sort()
        ret_df.at[combs,"Treatment time"] = t
        ret_df.at[combs,"Combination"] = "_".join(ds)
        ret_df.at[combs,"NumbDrugs"] = len(ds)
    return ret_df.sort_index(axis=1)


def make_sample_dict(metadata):
    ret_dict = {}
    for name, group in metadata.groupby(["Drug","Concentration","Treatment time (hrs)"]):
        key = "%s_%0.2f_%0.0f"%name
        ret_dict[key]=list(group.index)
    return ret_dict


def triangulate(data,sample_dict,combination,Tr):
    cs = combination.split("#")
    trans_pred = []
    UT = "No drug_0.00_%0.0f"%Tr
    samples_for_triangulation = sample_dict[UT]
    for c in cs:
        trans_pred.append(stats.mstats.gmean(data[sample_dict[c]],axis = 1))
        samples_for_triangulation.extend(sample_dict[c])
    trans_D1_D2_half = stats.mstats.gmean(data[samples_for_triangulation],axis = 1)
    trans_pred.append(trans_D1_D2_half)
    return trans_pred


def mlsyn(trained_classifier,data,sd,combination,Tr,ND):
    ref_CVS = 2.32497701e+12
    classifier = pickle.load(open(trained_classifier,"rb"))
    pred_trans = triangulate(data,sd,combination,Tr)
    CVS_di = [ref_CVS - classifier.score_samples([di]) for di in pred_trans[:-1]]
    CVS_dcomb_e = np.mean(CVS_di)
    CVS_dcomb_o = ref_CVS - classifier.score_samples([pred_trans[-1]])
    return CVS_di+[CVS_dcomb_e,CVS_dcomb_o]


def get_mlsyn_score(trained_classifier,data,metadata,combination_metadata):
    metabl = pd.read_csv(metadata,index_col=0)
    sd = make_sample_dict(metabl)
    metacomb = pd.read_csv(combination_metadata,index_col=0)
    for ni, (n, row) in enumerate(metacomb.iterrows()):
        trt = row["Treatment time"]
        if row.NumbDrugs in [1,2]:
            d1,d2,d12_e,d12_o = mlsyn(trained_classifier,data,sd,n,trt,2)
            mlsyn_sc = d12_e/d12_o
            expected, observed = d12_e, d12_o
        elif row.NumbDrugs == 3:
            d1,d2,d3,d123_e,d123_o = mlsyn(trained_classifier,data,sd,n,trt,3)
            mlsyn_sc = d123_e/d123_o
            expected, observed = d123_e,d123_o
        mlsyn_score = np.log2(mlsyn_sc + 0.3329)*10#Score corrected with the avergage ratio of additive combinations
        metacomb.at[n,"MLSynergy_score"] = mlsyn_score
        if mlsyn_score < 0:
            metacomb.at[n,"MLSynergy_label"] = "Synergy"
        elif mlsyn_score >= 0:
            metacomb.at[n,"MLSynergy_label"] = "Antagony"
    return metacomb


def plot_chord(predictions,filename):
    to_use = predictions.copy()
    for n,row in predictions.iterrows():
        if row.MLSynergy_score < 0:
            to_use.at[n,"MLSynergy_score"] = row.MLSynergy_score * -1
            to_use.at[n,"Interaction"] = "Synergy"
        else:
            to_use.at[n,"Interaction"] = "Antagony"
    hv.extension('bokeh')
    hv.output(size=200)
    to_use2 = to_use[to_use.NumbDrugs == 2]
    links = to_use2[["Drug1","Drug2","MLSynergy_score","Interaction"]]
    drugs = list(links["Drug1"].unique())+list(links["Drug2"].unique())
    nodes = hv.Dataset(drugs, 'Drug')
    chord = hv.Chord((links, nodes)).select(value=(1, None))
    chord.opts(opts.Chord(cmap='Rainbow', edge_cmap='Rainbow',\
                          edge_color=dim('Interaction').str(), labels='Drug',\
                          node_color=dim('Drug').str()))
    output_file(filename)
    show(hv.render(chord))
    return to_use2
        

def main(argv):
    #### help strings
    help_str1 = "MLSynergy.py --command <> --[conditional]classifier <> --[conditional]metadata <> --[conditional]combinaitons --[conditional]data <> --[conditional]ref_ids <> --output [OPTIONS]"
    options_str2 = "--help or -h: Display description of options\n\
--command or -c: <str> Make_combinations or Score_combinations\n\
--[conditional]classifier or -f: <str> If command = Score_combinations; location of the trained drona model (folder) \n\
--metadata or -m: <str> metadata file of single drug treatments (csv format)\n\
--[conditional]combinations or -p: <str> If command = Score_combinations; combinations file generated by Make_combination command (csv format)\n\
--[conditional]data or -d: <str> If command = Score_combinations; datafile of single drug treatment with gene expession data (csv format)\n\
--[conditional]output or -o: <str> If command = Make_combinations; file to write the combinations as (csv format) in\n\
OPTIONS:\n\
--plot_chord or -u: <bool> Yes or No for plotting predictions of two drug interaction; default: False\n\
--comb_number or -n: <comma sperated> If command = Make_combinations; Number of drug to make the combinations with; default = 2,3\n\
--treatment_durations or -t: <comma sperated> If command = Make_combinations; Treatment duration (in hours) to concider to make the combinations with; default = 72\n\
--rank_normalize or -z: <bool> True or False for rank normalisation; default: True\n\
--clean_data or -l: <bool> True or False for removing transcritomes with large number of missing data; default: False\n\
--[conditional]norm_method or -a: <str> If rank_normalize is True, select from [average,min,max,dense,ordinal]; default: min"
    #### define classical python input
    command,classifier,metadata,combinations,data,ref_ids,output = "","","","","","",""
    clean_data,normalize,norm_method,N,td,chrd = "False","True","min",[2,3],[72],"No"
    
    try:
        opts, args = getopt.getopt(argv,"h:c:f:m:p:d:o:z:l:a:n:t:u:",\
                                   ["help","command=","classifier=","metadata=","combinations=","data=","output=",\
                                    "rank_normalize=","clean_data=","norm_method=","comb_number=","treatment_durations=","plot_chord="])
    except getopt.GetoptError:
        print(help_str1+"\n"+options_str2)
        sys.exit(2)
        
    #### assign variable
    for opt, arg in opts:
        if opt in ("-h","--help"):
            print(help_str1+"\n"+options_str2)
            sys.exit()
        elif opt in ("-c","--command"):
            command = arg
        elif opt in ("-f","--classifier"):
            model_name = arg.split("/")[-1]
            classifier = "%s/%s_DRonA_model"%(arg,model_name)
            ref_ids = "%s/%s_ref_id"%(arg,model_name)
        elif opt in ("-m","--metadata"):
            metadata = arg
        elif opt in ("-p","--combinations"):
            combinations = arg
        elif opt in ("-d","--data"):
            data = arg
        elif opt in ("-o","--output"):
            output = arg
        elif opt in ("-z","--rank_normalize"):
            normalize = arg
        elif opt in ("-l","--clean_data"):
            clean_data = arg   
        elif opt in ("-a","--norm_method"):
            norm_method = arg
        elif opt in ("-n","--comb_number"):
            N = [int(i) for i in arg.split(",")]
        elif opt in ("-t","--treatment_durations"):
            td = [int(i) for i in arg.split(",")]
        elif opt in ("-u","--plot_chord"):
            chrd = arg
        else:
            pass

    #### run commands
    if command:
        ### If train command
        if command == "Make_combinations":
            if metadata:
                combinations = prepare_combinations_file(metadata,N,td)
                combinations.to_csv(output)
            else:
                print(help_str1+"\n"+options_str2)
                sys.exit(2)
        ### If score command
        elif command == "Score_combinations":
            if classifier and data and combinations and metadata:
                ### Read data
                exp_data = pd.read_csv(data,index_col=0)
                ### Reindex data
                reInd_exp_data = Data_prep.reindex_data(exp_data,ref_ids)
                ### Clean data
                if clean_data == "True":
                    cleaned_reInd_exp_data = Data_prep.clean_data(reInd_exp_data,0.70)
                else:## this is default
                    cleaned_reInd_exp_data = reInd_exp_data.copy()
                ### Rank normalized
                if normalize == "True":## this is default
                    normed_cleaned_reInd_exp_data = Data_prep.rank_normalize_across_genes(cleaned_reInd_exp_data,norm_method)
                else:
                    normed_cleaned_reInd_exp_data = cleaned_reInd_exp_data.copy()
                ### Get MLSynergy score
                mlsyn_score = get_mlsyn_score(classifier,normed_cleaned_reInd_exp_data,metadata,combinations)
                csvfilename = "%s_predicted.csv"%combinations.split(".")[0]
                mlsyn_score.to_csv(csvfilename)
                if chrd == "Yes":
                    htmlfilename = "%s_chordplot.html"%combinations.split(".")[0]
                    plot_chord(mlsyn_score,htmlfilename)
            else:
                print(help_str1+"\n"+options_str2)
                sys.exit(2)
        else:
            print(help_str1+"\n"+options_str2)
            sys.exit(2)
    else:
        print(help_str1+"\n"+options_str2)
        sys.exit(2)


if __name__ == "__main__":
    main(sys.argv[1:])

