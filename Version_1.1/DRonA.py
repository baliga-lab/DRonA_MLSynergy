#!/usr/bin/python3.7
"""@Author: Vivek Srinivas - Baliga Lab, ISB
- DRonA
"""


#1. Imports
import os, sys, pickle, getopt, Data_prep
import pandas as pd
import numpy as np
from sklearn.svm import OneClassSVM as OCSVM
from scipy import stats

#2. Functions

def create_gsm_dict(values, annotation_file):
    data = pd.read_csv(annotation_file, index_col = 0)
    ret_dict = {}
    for i,v in data.iterrows():
        ret_dict[i] = v[values]
    return ret_dict

def DRonA_trainer(metadata, data,**kwargs):
    """Transcript data of optimally grown cultures which are
    required for training classifiers to predict drug effect
    are not annotated and it is hard to sift through the dat-
    -asets to manually annotate them, therefore we would use
    this function to populate untreated samples by starting
    from a conservative set of manually annotated samples """
    #### precautionary filling of missing values
    data = data.fillna(0)
    #### define hue
    hue = kwargs["classify_by"]#kwargs
    #### Data curation
    gsm_labels = create_gsm_dict(hue,metadata)
    gsms_to_include = [i for i in data.columns if i in gsm_labels.keys()]
    data = data[gsms_to_include]
    print("Shape of data provided is %s,%s"%data.shape)
    #### Labels for classifying GSMs
    labels = kwargs["labels"]#kwargs
    #### Prepare initial starting training set
    untreated_GSMs = [gsm for gsm, label in gsm_labels.items() if label in [labels["seed_positives"]] and gsm in data.columns]
    treated_GSMs = [gsm for gsm, label in gsm_labels.items() if label in [labels["seed_negatives"]] and gsm in data.columns]
    null_GSMs = [gsm for gsm, label in gsm_labels.items() if label not in [labels["seed_negatives"],labels["seed_positives"]]\
                 and gsm in data.columns]
    initial_untreated_data = data[untreated_GSMs]
    treated_data_used_for_accuracy = data[treated_GSMs]
    null_data_used_for_iterative_training = data[null_GSMs]
    print("Shape of initial dataset is - %s,%s"%initial_untreated_data.shape)
    print("Lenght of initial test data - %s,%s"%treated_data_used_for_accuracy.shape)
    print("Lenght of initial null data used for inclusion - %s,%s"%null_data_used_for_iterative_training.shape)
    
    #### Make arrays used for training and testing
    initial_training_array = initial_untreated_data.transpose().values
    treated_data_array = treated_data_used_for_accuracy.transpose().values
    null_data_array = null_data_used_for_iterative_training.transpose().values

    
    #### Classifier parameters
    kernel_p = kwargs["kernel"]#kwargs
    degree_p = kwargs["degree"]#kwargs
    gamma_p = kwargs["gamma"]#kwargs
    gain = kwargs["gain"]#kwargs
    tolerance = kwargs["tolerance"]#kwargs
    threshold = kwargs["threshold"]#kwargs
    #### Define initial classifier   
    classifierID = OCSVM(kernel = kernel_p,degree = degree_p, gamma = gamma_p)
    
    #### Train initial classifier
    classifierID.fit(initial_training_array) ## Train
    pos_score = classifierID.score_samples(initial_training_array)## get pos scores
    neg_score = classifierID.score_samples(treated_data_array)## get neg scores
    pos_mean, pos_sigma = np.mean(pos_score), np.std(pos_score)
    initial_training_min_scores, initial_training_min_scores = stats.norm.interval(gain, loc=pos_mean, scale=pos_sigma)
    neg_mean, neg_sigma = np.mean(neg_score), np.std(neg_score)
    initial_test_min_scores, initial_test_max_scores = stats.norm.interval(tolerance, loc=neg_mean, scale=neg_sigma)
    print("Range of initial viable class is %s-%s"%(initial_training_min_scores,initial_training_min_scores))
    print("Range of non-viable class is %s-%s"%(initial_test_min_scores,initial_test_max_scores))
    false_positive = [i for i in neg_score if i > initial_training_min_scores]
    false_negative = [i for i in pos_score if i < initial_test_max_scores]
    tp,tn = len(pos_score)-len(false_negative), len(neg_score)-len(false_negative)
    fp,fn = float(len(false_positive)),len(false_negative)
    print([tp,fp,tn,fn])
    initial_accuracy = (tp+tn)/(tp+tn+fp+fn)
    print("Initial accuracy is %0.3f"%initial_accuracy)
    
    #### Iterative training
    run_details = pd.DataFrame()
    val1, val2 = 999, null_data_used_for_iterative_training.shape[1]
    new_accuracy = 100
    runN = 0
    classified_GSMs = list(initial_untreated_data.columns)
    classifiers ={"Classifier_0":classifierID}
    maxs = {"Classifier_0":initial_test_max_scores}
    run_details.at[runN,"Classifier"] = "Classifier_%s"%runN
    run_details.at[runN,"GSM_included"] = ",".join(initial_untreated_data.columns)
    run_details.at[runN,"Total_GSMs"] = len(classified_GSMs)
    run_details.at[runN,"Accuracy"] = initial_accuracy
    run_details.at[runN,"Run_number"] = runN
    run_details.at[runN,"Number_of_GSMs_included"] = len(classified_GSMs)
    for zzz in range(val2):
        classifierO = classifiers["Classifier_%s"%runN]
        max0 = maxs["Classifier_%s"%runN]
        ### Start run count
        runN += 1
        print("Run number:%s"%runN)
        scores = classifierO.score_samples(null_data_array)
        r1_GSMs_to_include = [GSM for GSM, score in zip(null_data_used_for_iterative_training.columns,scores)\
                              if score > max0 and GSM not in classified_GSMs]
        val1 = len(r1_GSMs_to_include)
        print("Number of new GSMs included = %s"%val1)
        if val1 > 0:
            classified_GSMs.extend(r1_GSMs_to_include)
            new_untreated_data = data[classified_GSMs]
            new_untreated_data_array = new_untreated_data.transpose().values
            print("s")
            classifierN = OCSVM(kernel = kernel_p,degree = degree_p, gamma = gamma_p)
            classifierN.fit(new_untreated_data_array)
            print(new_untreated_data_array.shape[0])
            print("e")
            new_pos_score = classifierN.score_samples(new_untreated_data_array)
            new_neg_score = classifierN.score_samples(treated_data_array)
            new_pos_mean, new_pos_sigma = np.mean(new_pos_score), np.std(new_pos_score)
            new_training_min_scores, new_training_min_scores = stats.norm.interval(0.9, loc=new_pos_mean, scale=new_pos_sigma)
            new_neg_mean, new_neg_sigma = np.mean(new_neg_score), np.std(new_neg_score)
            new_test_min_scores, new_test_max_scores = stats.norm.interval(0.7, loc=new_neg_mean, scale=new_neg_sigma)
            maxs["Classifier_%s"%runN] = new_test_max_scores
            print("Range of training is %s-%s"%(new_training_min_scores,new_training_min_scores))
            print("Range of testing is %s-%s"%(new_test_min_scores,new_test_max_scores))
            new_false_positive = [i for i in new_neg_score if i > new_training_min_scores]
            new_false_negative = [i for i in new_pos_score if i < new_test_max_scores]
            new_tp,new_tn = len(new_pos_score)-len(new_false_negative), len(new_neg_score)-len(new_false_negative)
            new_fp,new_fn = float(len(new_false_positive)),len(new_false_negative)
            print([new_tp,new_fp,new_tn,new_fn])
            new_accuracy = (new_tp+new_tn)/(new_tp+new_tn+new_fp+new_fn)
            print("\tAccuracy in identifying known drug treated sets = %0.3f"%new_accuracy)
            #--------------------
            run_details.at[runN,"Total_GSMs"] = new_untreated_data_array.shape[0]
            run_details.at[runN,"Accuracy"] = new_accuracy
            run_details.at[runN,"Run_number"] = runN
            run_details.at[runN,"Number_of_GSMs_included"] = len(r1_GSMs_to_include)
            if new_accuracy >=threshold:
                classifiers["Classifier_%s"%runN] = classifierN
                run_details.at[runN,"GSM_included"] = ",".join(r1_GSMs_to_include)
                run_details.at[runN,"Classifier"] = "Classifier_%s"%runN
            else:
                run_details.at[runN,"GSM_included"] = "NA"
                run_details.at[runN,"Classifier"] = "NA"
                print("\tIteration broken at run %s as accuracy of the new classifier was below %s"%(runN,accuracy_tresh))
                final_classifier = classifierN
                break
        else:
            run_details.at[runN,"Total_GSMs"] = new_untreated_data_array.shape[0]
            run_details.at[runN,"Accuracy"] = new_accuracy
            run_details.at[runN,"Run_number"] = runN
            run_details.at[runN,"GSM_included"] = "NA"
            run_details.at[runN,"GSM_included"] = "NA"
            run_details.at[runN,"Number_of_GSMs_included"] = 0
            run_details.at[runN,"Classifier"] = "NA"
            print("\tIteration broken at run %s as no new GSMs were added to the viable set training data"%runN)
            final_classifier = classifierN
            break
        print("Overall shape of the data at the end of run %s = %s"%(runN,new_untreated_data.shape))
    #### returns    
    return run_details, final_classifier

def return_scores2(classifier,data,metadata):
    meta_data = metadata.copy()
    scores = {gsm:score for gsm,score in zip(data.columns,classifier.score_samples(data.transpose().values))}
    ref = max(scores.values())
    for n, row in meta_data.iterrows():
        try:
            meta_data.at[n,"CVS"] = scores[n] - ref
        except:
            pass
    return meta_data

def add_CVS(classifier_loc,data, meta_data_loc):
    classifier = pickle.load(open(classifier_loc,"rb"))
    meta_data = pd.read_csv(meta_data_loc,index_col=0)
    expt_groups = meta_data.groupby("Experiment")
    scores_list = []
    for name,group in expt_groups:
            samples = [i for i in group.index if i in data.columns]
            sub_data = data[samples]
            scores_list.append(return_scores2(classifier,sub_data, group))
    ret_data = pd.concat(scores_list)
    return ret_data



def main(argv):
    #### help strings
    help_str1 = "DRonA.py --command <> --classifier[conditional] <> --metadata <> --data <> --ref_ids <> --output <> [OPTIONS]"
    options_str2 = "--help or -h: <str> Lists inputs and their descriptions\n\
--command or -c: <str> Train or Score\n\
--[conditional]classifier or -f: <str> If command = Score_combinations; location of the trained drona model (folder) \n\
--metadata or -m: <str> Metadata file (csv format)\n\
--data or -d: <str> Datafile (csv format) with gene expession data\n\
OPTIONS:\n\
--output or -o: If command = Train; Provide name for the trained model; default will be Unnamed\n\
--rank_normalize or -z: <bool> True or False for rank normalisation; default: True\n\
--clean_data or -l: <bool> True or False for removing transcritomes with large number of missing data; default: False\n\
--cb or -q: <str> Label of row in metadata file having data labels for training and testing; default: Sample_class\n\
--sp or -p: <str> Label for iniital seed data; default: Manually IDd viable\n\
--sn or -n: <str> Label for test data; default: Manually IDd non-viable\n\
--k or -k: <str> Kernel [linear,poly,rbf,sigmoid,precomputed]; default: linear\n\
--[conditional]dg or -t: <int: 1 to inf> if kernel is poly; default: 1\n\
--ga or -g: <str> Gamma [scale,auto]; default: auto\n\
--gn or -w: <float: 0.0 to 1.0> Gain; default: 0.9\n\
--to or -x: <float: 0.0 to 1.0> Tolerance; default: 0.9\n\
--th or -y: <float: 0.0 to 1.0> Threshold; default: 0.85\n\
--[conditional]norm_method or -a: <str> If rank_normalize is True, select from[average,min,max,dense,ordinal]; default: min"
    #### define classical python input
    command,classifier,metadata,data,ref_ids = "","","","",""
    output,clean_data,normalize,norm_method,classify_by,seed_positives,seed_negatives = "Unnamed","False","False","min","Sample_class","Manually IDd viable","Manually IDd non-viable"
    kernel,degree,gamma,gain,tolerance,threshold = "linear",1,"auto",0.9,0.9,0.85
    try:
        opts, args = getopt.getopt(argv,"h:c:f:m:d:r:o:z:l:q:p:n:k:t:g:w:x:y:a:",\
                                   ["help","command=","classifier=","metadata=","data=","ref_ids=","output=",\
                                    "rank_normalize=","clean_data=","cb=","sp=","sn=","k=","dg=",\
                                    "ga=","gn=","to=","th=","norm_method="])
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
            ref_ids = "%s/%s_ref_ID"%(arg,model_name)
        elif opt in ("-m","--metadata"):
            metadata = arg
        elif opt in ("-d","--data"):
            data = arg
        elif opt in ("-o","--output"):
            output = arg
        elif opt in ("-z","--rank_normalize"):
            normalize = arg
        elif opt in ("-l","--clean_data"):
            clean_data = arg   
        elif opt in ("-q","--cb"):
            classify_by = arg
        elif opt in ("-p","--sp"):
            seed_positives = arg
        elif opt in ("-n","--sn"):
            seed_negatives = arg
        elif opt in ("-k","--k"):
            kernel = arg
        elif opt in ("-t","--dg"):
            degree = int(arg)
        elif opt in ("-g","--ga"):
            gamma = arg
        elif opt in ("-w","--gn"):
            gain = float(arg)
        elif opt in ("-x","--to"):
            tolerance = float(arg)
        elif opt in ("-y","--th"):
            threshold = float(arg)
        elif opt in ("-a","--norm_method"):
            norm_method = arg
            
    labels = {"seed_positives":seed_positives, "seed_negatives":seed_negatives}
        
    #### run commands
    if command and metadata and data:
        ### Read data
        exp_data = pd.read_csv(data,index_col=0)
        ### Reindex data
        if command == "Score":
            print("Re-indexing the gene expression data")
            reInd_exp_data = Data_prep.reindex_data(exp_data,ref_ids)
        else:
            reInd_exp_data = exp_data.copy()
        ### Clean data
        if clean_data == "True":
            print("Cleaning the gene expression data; Removing samples with missing values")
            cleaned_reInd_exp_data = Data_prep.clean_data(reInd_exp_data,0.70)
        else:## this is default
            cleaned_reInd_exp_data = reInd_exp_data.copy()
        ### Rank normalized
        if normalize == "True":## this is default
            print("Normalizing the gene expression data with %s method"%norm_method)
            normed_cleaned_reInd_exp_data = Data_prep.rank_normalize_across_genes(cleaned_reInd_exp_data,norm_method)
        else:
            normed_cleaned_reInd_exp_data = cleaned_reInd_exp_data.copy()
        ### If train command
        if command == "Train":
            if output:
                pass
            else:
                print("Trained model will be saved as Unnamed")
                output = "Unnamed"
            if os.path.exists(output):
                pass
            else:
                os.mkdir(output)
            rd,cls = DRonA_trainer(metadata,normed_cleaned_reInd_exp_data,\
                                   classify_by=classify_by,labels=labels,kernel=kernel,degree=degree,\
                                   gamma=gamma,gain=gain,tolerance=tolerance,threshold=threshold)
            rd.to_csv("%s/%s_run_details.csv"%(output,output))
            with open("%s/%s_ref_ID"%(output,output),"w") as f:
                f.writelines("\n".join(normed_cleaned_reInd_exp_data.index))
            pickle.dump(cls,open("%s/%s_DRonA_model"%(output,output),"wb"))
            os.system("cp %s %s/%s_ref_ID"%(ref_ids,output,output))
            
        ### If score command
        elif command == "Score":
            if classifier:
                print("Calculating scores")
                scores = add_CVS(classifier,normed_cleaned_reInd_exp_data,metadata)
                scores.to_csv("%s_scored.csv"%".".join(metadata.split(".")[:-1]))
                print("Write the scores to %s file"%"%s_scored.csv"%".".join(metadata.split(".")[:-1]))
            else:
                print(help_str1+"\n"+"Provide location of trained classifier")
        else:
            print(help_str1+"\n"+"Provide location of trained classifier")
            sys.exit(2)
    else:
        print(help_str1+"\n"+options_str2)
        sys.exit(2)
        
if __name__ == "__main__":
    main(sys.argv[1:])
