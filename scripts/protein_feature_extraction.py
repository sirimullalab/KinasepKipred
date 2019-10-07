# Script for extracting protein features
# Author: Md Mahmudulla Hassan
# Last Modified: 03/15/2019
# Works only with python 2.7 because of PyPro's dependency on this python version

from pydpi.pypro import PyPro
import pickle
import pandas as pd
import json
import multiprocessing as mp
import os

DATA_DIR = './data'
# Read the gene seqeuences from the csv file

print("Reading gene sequences...")
gene_sequences = pd.read_csv(os.path.join(DATA_DIR, "gene_sequences_test.csv"), index_col=0)
sequence_list = gene_sequences.sequence.tolist()

def get_features(seq):
    p = PyPro()
    try:
      p.ReadProteinSequence(seq)
      features = list(p.GetALL().values())
      return features
    except:
      return ''

#print "Extracting protein features using all " + str(mp.cpu_count()) + " threads..."
pool = mp.Pool(mp.cpu_count())
result = pool.map(get_features, sequence_list)
gene_list = gene_sequences.index.tolist()
gene_dict = {}
for g, s, f in zip(gene_list, sequence_list, result):
    gene_dict[g] = {'sequence': s, 'features': f}
#with open('gene_dict_2.pickle', 'wb') as f:
#    pickle.dump(gene_dict, f)

#print "Features have been extracted for " + str(len(gene_dict)) + " proteins"
#print "Writing protein features into data/gene_dict.json ..."
with open(os.path.join(DATA_DIR, "gene_dict_test.json"), 'w') as f:
    json.dump(gene_dict, f)
pool.close()

print("Done")
