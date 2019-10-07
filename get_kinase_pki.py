from __future__ import print_function, absolute_import 

# Script to predict (or test) the model using protein (kinase) sequence and SMILE pattern of a compound.
# Usage: python2 get_kinase_pki.py protein_sequence "SMILE_Pattern"

import numpy as np
from pydpi.pypro import PyPro
import pandas as pd
import json
import multiprocessing as mp
import os
import sys
import numpy as np
from sklearn.externals import joblib
from scripts import FeatureGenerator
#from keras.models import load_model
import pickle

class pKiPred(object):
    
    def __init__(self):
        self.model = joblib.load('models/Random_forest_gridsearch_py27.mdl')

    def get_smi_features(self, smiles):
        try:
            feat_gen = FeatureGenerator(smiles)
            features = feat_gen.toTPATF()
            return features
        except:
            return None

    def get_features(self, seq, smi):
        p = PyPro()
        try:
          p.ReadProteinSequence(seq)
          features = list(p.GetALL().values())
      
          smi_features = self.get_smi_features(smi)
          smi_features2 = list(np.array([f for f in smi_features], dtype=np.float32))
      
          total_features = np.array(features+smi_features2)[np.newaxis, :]
    #      total_features = np.array(smi_features2+features)[np.newaxis, :] # does not work...!
          return total_features
        except Exception as e:
          print(str(e))
        return None

    def predict(self, seq, smi):
        protein_feature = self.get_features(seq, smile)
        return self.model.predict(protein_feature)



if __name__=='__main__':
    seq = "MGCGCSSHPEDDWMENIDVCENCHYPIVPLDGKGTLLIRNGSEVRDPLVTYEGSNPPASPLQDNLVIALHSYEPSHDGDLGFEKGEQLRILEQSGEWWKAQSLTTGQEGFIPFNFVAKANSLEPEPWFFKNLSRKDAERQLLAPGNTHGSFLIRESESTAGSFSLSVRDFDQNQGEVVKHYKIRNLDNGGFYISPRITFPGLHELVRHYTNASDGLCTRLSRPCQTQKPQKPWWEDEWEVPRETLKLVERLGAGQFGEVWMGYYNGHTKVAVKSLKQGSMSPDAFLAEANLMKQLQHQRLVRLYAVVTQEPIYIITEYMENGSLVDFLKTPSGIKLTINKLLDMAAQIAEGMAFIEERNYIHRDLRAANILVSDTLSCKIADFGLARLIEDNEYTAREGAKFPIKWTAPEAINYGTFTIKSDVWSFGILLTEIVTHGRIPYPGMTNPEVIQNLERGYRMVRPDNCPEELYQLMRLCWKERPEDRPTFDYLRSVLEDFFTATEGQYQPQP"
    smile = "CC(C)Oc1ccc(cc1Cl)c2noc(n2)c3ccc(N[C@H]4CC[C@H](C4)C(=O)O)cc3"
    pkipred = pKiPred()

    if len(sys.argv) == 1:
        print(pkipred.predict(seq, smile))
    else:
        print(pkipred.predict(sys.argv[1], sys.argv[2]))
