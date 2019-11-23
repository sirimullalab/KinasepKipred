# KinasepKipred
Kinases to predict Inhibitor constant in terms of pKI (where pKi is decadic logarithm of Ki). We used the data points that were specifically represent Ki values to train and test the models.

## Requirements:
* python==2.7.16
* pydpi==1.0
* numpy==1.16.5
* pandas==0.24.2
* tqdm==4.36.1
* scipy==1.2.2
* scikit-learn==0.20.4
* rdkit==2018.03

## Setup
* Download the model file 
```bash
./download_model.sh
```
* Set up the conda environment and activate it
```bash
conda env create -f environment.yml
conda activate kinasepki
```

## Usage:
Run `./test.sh` to get the prediction for an example pair of protein sequence and a ligand SMILES. For any other inputs, run the following

```bash
python2 get_kinase_pki.py "<PROTEIN_SEQUENCE>" "<LIGAND_SMILES>"
```

## Docker 
* Download the docker image `docker pull sirimullalab/kinasepkipred:py2`
* Run the container `docker run --rm sirimullalab/kinasepkipred:py2 <protein_sequence> <compound_smiles>`. To run with a built-in sample, do  `docker run --rm sirimullalab/kinasepkipred:py2`
