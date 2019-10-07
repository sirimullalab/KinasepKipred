# KinasepKipred
This repository contains codes for three different models: Random Forest, XGBoost, and Artificial Neural Network. This work is specially done for Kinases to predict binding affinity interms of pKI (where pKI is decadic logarithm of Ki). We used only Ki values to train and test the models.

## Requirements:
Please, install the following in order to run the codes.
* python==2.7.16
* pydpi==1.0
* numpy==1.16.5
* pandas==0.24.2
* tqdm==4.36.1
* scipy==1.2.2
* scikit-learn==0.20.4
* rdkit==2018.03

## For the protein features:
Download pydpi-1.0 (by Dongsheng Cao and Yizeng Liang) which is available in this same repository or pip install pydpi
## For Ligand Features:
Please download MayaChemTools from http://www.mayachemtools.org/Download.html
We have used TPATPF fingerprints. 'generate_ligand_fp.py' can be used to generate fp.

## For XGboost:
* pip install xgboost

## For Artificial Neural Network:

* pip install keras
* pip install tensorflow

## For plots
* pip install matplotlib
* pip install docopt

## Usage:
How to predict the ligand kinase pKi:
```
python get_kinase_pki.py protein_kinase_sequence "SMILE_pattern"
```
(You can use 'test_protein_sequence.txt' and 'smiles.txt' as test files)
