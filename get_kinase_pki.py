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
      total_features = np.array(smi_features2+features)[np.newaxis, :]
      return total_features
    except Exception as e:
      print((str(e)))
      return None

protein_feature = get_features(sys.argv[1], sys.argv[2])
load_model = joblib.load('models/random_forest_gridsearch_pki.mdl')
y_pred = load_model.predict(protein_feature)
print(y_pred)
