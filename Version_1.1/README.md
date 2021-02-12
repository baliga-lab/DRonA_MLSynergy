# DRonA and MLSynergy
**Repository for GEOparser, DRonA and MLSynergy algorithms**

##### Requires Python 2.7 for GEOparser and Python 3.0 for DRonA and MLSynergy

## Downloading git contents:

```
git clone https://github.com/viveks-baligalab/DRonA_MLSynergy.git
```

## Installing packages required to run algorithms:
```
pip2 install -r GEOparser_requirements.txt
pip3 install -r DRonA_MLSynergy_requirements.txt
```

## Running instructions:
**GEOparser -**
```
python2 GEOparser.py --help --prjname <> --coldir <> --refseq <> -gsmidfile <> [OPTIONS]
```
**DRonA -**

*Get help*
```
python3 DRonA.py -h

```
*Train DRonA model*
```
python3 DRonA.py -c Train -m Training_data/GEO_metadata_for_DRonA.csv -d Training_data/GEO_data_for_DRonA.csv -o New_model
```
*Use trained DRonA model to score drug treatments*
```
python3 DRonA.py -c Score -f DRonA_trained_models/MTB_2020 -m Single_drug_treatments/Single_drug_treatments_metadata.csv -d Single_drug_treatments/Single_drug_treatments_data.csv
```

**MLSynergy -**

*Get help*
```
python3 DRonA.py -h

```
*Use MLSynergy to make the drug combination metadata file*
```
python3 MLSynergy.py -c Make_combinations -m Drug_synergy_prediction/SDTs_used_for_predicting_drug_synergy.csv -o Drug_combinations_2020.csv

```
*Use MLSynergy to score the drug combinations*
```
python3 MLSynergy.py -c Score_combinations -f DRonA_trained_models/MTB_2020 -m Single_drug_treatments/Single_drug_treatments_metadata.csv -d Single_drug_treatments/Single_drug_treatments_data.csv -p Drug_synergy_prediction/Drug_combinations_2020.csv

```
