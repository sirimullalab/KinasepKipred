# KinasepKipred
This repository contains codes for three different models: Random Forest, XGBoost, and Artificial Neural Network. This work is specifically done for Kinases to predict Inhibitor constant interms of pKI (where pKi is decadic logarithm of Ki). We used the data points that were specifically estimated for Ki values to train and test the models.

## Requirements:
Install the following in order to run the codes.
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
Download MayaChemTools from http://www.mayachemtools.org/Download.html
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
Run `./test.sh`. It basically runs the python script `get_kinase_pki.py`

## Docker 
* Download the docker image `docker pull sirimullalab/kinasepkipred:py2`
* Run the container `docker run --rm sirimullalab/kinasepkipred:py2 <protein_sequence> <compound_smiles>`. To run with a built-in sample, do  `docker run --rm sirimullalab/kinasepkipred:py2`
