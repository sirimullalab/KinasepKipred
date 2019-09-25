# kinasepkipred
This repository contains codes for three different models: Random Forest, XGBoost, and Artificial Neural Network. This work is specially done for Kinases to predict binding affinity interms of pKI (where pKI is decadic logarithm of Ki). We used only Ki values to train and test the models.

## Requirements:
Please, install the following in order to run the codes.
* python>= 3.6
* numpy==1.16.4
* pandas==0.24.2
* tqdm==4.32.1
* scipy==1.2.2
* scikit-learn==0.21.3
* rdkit==2019.03.2

## For the protein features:
Use pydpi (by Dongsheng Cao and Yizeng Liang) which is available in this same repository.

## For Ligand Features:
Download MayaChemTools from http://www.mayachemtools.org/Download.html. 
We have used TPATPF fingerprints. See 'utility.py' file for details.

## For XGboost:
* pip install xgboost

## For Artificial Neural Network:

* pip install keras
* pip install tensorflow



## Usage:
How to use to predict the ligand kinase binding affinity interms of pKI:
```
python get_kinase_pki.py protein_kinase_sequence "smiles"
```
(You can use 'test_protein_sequence.txt' and 'smiles.txt' as test files)
