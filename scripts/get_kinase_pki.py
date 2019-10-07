# Script to predict (or test) the model using protein (kinase) sequence and SMILE pattern of a compound.
# Usage: python2 get_kinase_pki.py protein_sequence "SMILE_Pattern"

import numpy as np
from pydpi.pypro import PyPro
import pandas as pd
import json
import multiprocessing as mp
from sklearn.externals import joblib
import os
import sys
import numpy as np
from utility import FeatureGenerator
#from keras.models import load_model
import pickle
def get_smi_features(smiles):
    try:
        feat_gen = FeatureGenerator(smiles)
        features = feat_gen.toTPATF()
        return features
    except:
        return None


def get_features(seq, smi):
    p = PyPro()
    try:
      p.ReadProteinSequence(seq)
      features = list(p.GetALL().values())
      
      smi_features = get_smi_features(smi)
      smi_features2 = list(np.array([f for f in smi_features], dtype=np.float32))
      
      total_features = np.array(features+smi_features2)[np.newaxis, :]
#      total_features = np.array(smi_features2+features)[np.newaxis, :] # does not work...!
      return total_features
    except Exception as e:
      print str(e)
      return None

if __name__=='__main__':
    protein_feature = get_features(sys.argv[1], sys.argv[2])
    load_model = joblib.load('Random_forest_gridsearch_py27.mdl')
    y_pred = load_model.predict(protein_feature)
    
    print y_pred
