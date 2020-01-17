# Scripts for combining protein and ligand features
# *************************************************************
# Author: Govinda KC                                          #
# UTEP, Computational Science                                 #
# Last modified: 1/17/2020                                    #
# *************************************************************

import pickle
import numpy as np
import pandas as pd
import glob
import os
from tqdm import *
import random
import json

DATA_DIR = '../data/'

print "Reading gene dictionary"
with open(os.path.join(DATA_DIR, 'gene_dict.json'), 'r') as f:
    gene_dict = json.load(f)

print "Reading tpatf features"
data_files = glob.glob(os.path.join(DATA_DIR, 'ligand_npy/*.npy'))
data_genes = [os.path.split(i)[1][:-4] for i in data_files]

# Filter out the genes that we don't have a sequence for
data_genes = [i for i in data_genes if i in gene_dict.keys()]
print "Number of genes: " + str(len(data_genes))
print "Extracting features..."
# First sample
indices = {} # records the indices of the gene samples
data = np.load(os.path.join(DATA_DIR, 'ligand_npy/' + data_genes[0] + '.npy'))
indices[data_genes[0]] = [i for i in range(data.shape[0])]
protein_feature = gene_dict[data_genes[0]]['features']
features = np.hstack((np.tile(protein_feature, (data.shape[0], 1)), data)).astype(np.float32)

for i in tqdm(range(1, len(data_genes))):
    _data = np.load(os.path.join(DATA_DIR, 'ligand_npy/' + data_genes[i] + '.npy'))
    indices[data_genes[i]]= [k for k in range(features.shape[0], features.shape[0] + _data.shape[0])]
    
    _protein_feature = gene_dict[data_genes[i]]['features']
    _features = np.hstack((np.tile(_protein_feature, (_data.shape[0], 1)), _data)).astype(np.float32)
    features = np.vstack((features, _features))
np.save(DATA_DIR+'kinase_and_lignad_features.npy',features)
print "File is saved"
