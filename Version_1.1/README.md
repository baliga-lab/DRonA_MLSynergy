# DRonA and MLSynergy

_This is a repository for GEOparser, DRonA and MLSynergy algorithms, developed in [__Baliga Lab,ISB__](https://baliga.systemsbiology.net/) for the study titled - [__"Transcriptome signature of cell viability predicts drug response and drug interaction for tuberculosis"__](https://www.biorxiv.org/content/10.1101/2021.02.09.430468v2)_
##
### Requirements -
- Python 2.7
- Python 3.0 or more
## 

#### Downloading DRonA and MLSynergy from GitHub:

```
git clone https://github.com/viveks-baligalab/DRonA_MLSynergy.git
```  

#### Installing the python packages:

```
pip2 install -r GEOparser_requirements.txt
```
```
pip3 install -r DRonA_MLSynergy_requirements.txt
```

## GEOparser

*_Before running GEOparser, download and install (if required)_*

-  [kallisto](https://pachterlab.github.io/kallisto/)

-  [Entrez direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

-  [SRA toolkit](https://ncbi.github.io/sra-tools/)

-  [Aspera connect](https://www.ibm.com/aspera/connect/)

-  [Aspera key path](https://www.ncbi.nlm.nih.gov/sra/docs/aspera-key-pairs/)

*and update the path names in paths.txt file in GEOparse_helpers folder*

*In order to download transcriptome data from GEO with GEOparser, create a text file with required GSM ids (seperated by comma)*. Eg. -
```
GSM3230595,GSM3230594,GSM3230593,GSM3229928
```
#### Getting help - 

```
python2 GEOparser.py --help --prjname <> --coldir <> --refseq <> -gsmidfile <> [OPTIONS]
```

#### Running GEOparser -


```
python2 GEOparser.py --prjname GEO_data --coldir Training_data --refseq Reference_genomes/MTB_H37Rv_transcripts.fa -gsmidfile Training_data/MTB_GSM_sample_list.txt
```

## DRonA

#### Data preperation for DRonA -
+ *__Metadata file__ associated with the transcriptome data, used for __training__ of DRonA has to be formatted in this manner -*

|         | Sample_class           | 
| ------------- |:-------------|
| GSM3230595    | Manually IDd non-viable | 
| GSM3230594    | Manually IDd viable      |   
| GSM3230593 |       | 
| ... | ...      |   

If the name of the column (__Sample_class__) or labels (__Manually IDd non-viable, Manually IDd viable__) of the transcriptomes are changed, then use options (__cb, sp__ and __sn__) in DRonA to specify the changes.

+ *__Metadata file__ associated with the transcriptome data, used for __prediction of drug response__ has to be formatted in this manner -*

|     Samples    | Experiment           | 
| ------------- |:-------------|
| Drug1_treated_E1    | Exp1 | 
| Drug1_treated_E2    | Exp2      |   
| Drug3_treated_E1 |     Exp1  | 
| ... | ...      |   

 Additional columns describing the samples can be added to the metadata file. 

+ *__Trancriptome data__ used for __training of DRonA__ or __prediction of drug response__ has to be formatted in this manner -*

|      Locus_tag   | GSM3230595           | GSM28267|...|
| ------------- |:-------------|:-------------|:-------------|
|   Rv0001  | 3178 |1234 |...|
|  Rv0002   | 1235      |  2453 |...
| Rv0003 |   2345    | 1111|...|
| ... | ...      | ...|...|

 __or__
 
|      Locus_tag   | Drug1_treated_E1           | Drug1_treated_E2|...|
| ------------- |:-------------|:-------------|:-------------|
|   Rv0001  | 3178 |1234 |...|
|  Rv0002   | 1235      |  2453 |...
| Rv0003 |   2345    | 1111|...|
| ... | ...      | ...|...|

#### Getting help -

```
python3 DRonA.py -h
```

#### Training a new DRonA model -

```
python3 DRonA.py -c Train -m Training_data/GEO_metadata_for_DRonA.csv -d Training_data/GEO_data_for_DRonA.csv -o New_model
```

#### Using trained DRonA model to predict drug response - 

```
python3 DRonA.py -c Score -f DRonA_trained_models/MTB_2020 -m Single_drug_treatments/Single_drug_treatments_metadata.csv -d Single_drug_treatments/Single_drug_treatments_data.csv
```

## MLSynergy -

#### Data preperation for MLSynergy -

+ *__Metadata file__ associated with the transcriptome data, used for __generation of the drug combinations file__, and  __predicition of drug interactions__ has to be formatted in this manner -*

|     Samples    | Drug           | Concentration| Treatment time (hrs)|
| ------------- |:-------------|:-------------|:-------------|
| Drug1_treated_E1    | Drug1 | 1.0|72|
| Drug1_treated_E2    | Drug1      |  2.0| 72 |
| Drug3_treated_E1 |     Drug3  | 1.0 | 24|
| No_drug1_E1 |     No drug  | 0.0 | 24 |
| ... | ...      |   ...|...|

Additional columns describing the samples can be added to the metadata file. 

+ *__Single drug treatments table__ used for __generation of the drug combinations file__ has to be formatted in this manner -*

| Drug           | Concentration| Treatment time (hrs)|
| -------------|:-------------|:-------------|
| Drug1 | 1.0|72|
| Drug2 |  2.0| 72 |
| Drug3  | 1.0 | 72|
| ...      |   ...|...|

+ *__Trancriptome data__ used for __predicition of drug interactions__ has to be formatted in this manner -*

|      Locus_tag   | Drug1_treated_E1           | Drug1_treated_E2|...|
| ------------- |:-------------|:-------------|:-------------|
|   Rv0001  | 3178 |1234 |...|
|  Rv0002   | 1235      |  2453 |...
| Rv0003 |   2345    | 1111|...|
| ... | ...      | ...|...|

#### Getting help -

```
python3 DRonA.py -h
```

#### Creating the drug combination metadata file - 

```
python3 MLSynergy.py -c Make_combinations -m Drug_synergy_prediction/SDTs_used_for_predicting_drug_synergy.csv -o Drug_combinations_2020.csv
```

#### Predicting drug interactions with the trained DRonA model*

```
python3 MLSynergy.py -c Score_combinations -f DRonA_trained_models/MTB_2020 -m Single_drug_treatments/Single_drug_treatments_metadata.csv -d Single_drug_treatments/Single_drug_treatments_data.csv -p Drug_synergy_prediction/Drug_combinations_2020.csv
```
