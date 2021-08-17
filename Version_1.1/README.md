# DRonA and MLSynergy
**Repository for GEOparser, DRonA and MLSynergy algorithms**

*DRonA: Drug responce analyzer, an algorithm for predicting drug responce in bacteria with transcriptomes from drug treated cultures*
*MLSynergy: An algorithm for predicting drug synergism and antagonism in bacteria with  transcriptomes from single drug treated cultures*
*GEOParser: An algorithm for collecting and reanalyzing microarray and RNAseq data from GEO*

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

*Before running GEOparser, download and install (if required)*
*-[kallisto](https://pachterlab.github.io/kallisto/)*
*-[Entrez direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)*
*-[SRA toolkit](https://ncbi.github.io/sra-tools/)*
*-[Aspera connect](https://www.ibm.com/aspera/connect/)*
*-[Aspera key path](https://www.ncbi.nlm.nih.gov/sra/docs/aspera-key-pairs/)*
*and update the path names in GEOparse_helpers/paths.txt file*

*Get help*
```
python2 GEOparser.py --help --prjname <> --coldir <> --refseq <> -gsmidfile <> [OPTIONS]
```
*Running GEOparser*
```
python2 GEOparser.py --prjname GEO_data --coldir Training_data --refseq Reference_genomes/MTB_H37Rv_transcripts.fa -gsmidfile Training_data/MTB_GSM_sample_list.txt
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
