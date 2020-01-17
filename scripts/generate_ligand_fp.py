# Script to generate TPATFP fingerprints based MayaChemTools
# *************************************************************
# Author: Govinda KC                                          #
# UTEP, Computational Science                                 #
# Last modified: 1/17/2020                                    #
# *************************************************************

from rdkit.Chem import AllChem, Descriptors
from rdkit import Chem
from tqdm import tqdm
from sklearn.preprocessing import LabelEncoder
import sys
import pandas as pd
import numpy as np
from rdkit.Chem import MACCSkeys
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import RDKFingerprint
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
import os, glob
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.preprocessing import MinMaxScaler
import tempfile
import shutil

class Ligand_Features:
    def __init__(self, csv_path):
        self.csv_path = csv_path
        self.temp_dir = tempfile.mkdtemp()
        
    def toSDF(self, smiles):
        """
        Converts smiles into sdf format which is used to generate fingerprints in TPATF, TPAPF, and PHYC
        """
        # Get mol format of smiles
        mol = Chem.MolFromSmiles(smiles)
        
        # Compute 2D coordinates
        AllChem.Compute2DCoords(mol)
        mol.SetProp("smiles", smiles)
        
        w = Chem.SDWriter(os.path.join(self.temp_dir, "temp.sdf"))
        w.write(mol)
        w.flush()

    def _cleanup(self):
        """
        Removes the temporary files temp.sdf and temp.csv files
        """
        shutil.rmtree(self.temp_dir)

    def toTPATF(self):
        """
        Calculates the topological pharmacophore atomic triplets fingerprints
        Parameters
        ----------
        input : sdf file
            Sdf file is created using 'to_SDF()' 
        
        return : list
            returns the features in list form
        """
        features = []
        script_path = "/Users/gvin/Downloads/mayachemtools/bin/TopologicalPharmacophoreAtomTripletsFingerprints.pl"
          
        # Now generate the TPATF features
        # Check if the sdf file exists
        
        if not os.path.isfile(os.path.join(self.temp_dir, "temp.sdf")): return None
        command = "perl " + script_path + " -r " + os.path.join(self.temp_dir, "temp") + " --AtomTripletsSetSizeToUse FixedSize -v ValuesString -o " + os.path.join(self.temp_dir, "temp.sdf")
        os.system(command)
        
        with open(os.path.join(self.temp_dir, "temp.csv"), 'r') as f:
            for line in f.readlines():
                if "Cmpd" in line:
                    line = line.split(';')[5].replace('"','')
                    features = [int(i) for i in line.split(" ")]

        return features
        
        print('Final Numpy array shape: {}'.format(final_array_unique.shape))
        print('Type of final array: {}'.format(type(final_array_unique)))
        final_numpy_array = np.asarray((final_array_unique), dtype = np.float32)
        return final_numpy_array

    def tpatf_fp(self):
        """
        Receives the csv file which is used to generate TPATF fingerprints (2692) and saves as numpy file
        
        Parameter
        ---------
        
        input smiles : str
            Compouds in the form of smiles are used
    
        return : np.array
            Features are saved in the form of numpy files
        
        """
        df= pd.read_csv(self.csv_path)
        smiles_list = df['smiles'].tolist()
        
        fingerprints = []
        not_found = []

        for i in tqdm(range(len(smiles_list))):
            try:
                self.toSDF(smiles_list[i])  
                
                features = fg.toTPATF()
                
                fingerprints.append(features)
            except:
                
                fingerprints.append(np.nan)
                not_found.append(i)
                pass

        # Clean up the temporary files
        self._cleanup()
        
        # Drops rowns if features are not found
        df.drop(not_found, axis=0,inplace=True)
        print('Number of FPs not found: {}'.format(len(not_found)))
        
        df.reset_index(drop=True, inplace=True)
        
        # Encoding categorical data
        labelencoder = LabelEncoder()                       
        Y = labelencoder.fit_transform(df['pki_value'].values)
        Y = Y.reshape(Y.shape[0],1)
        
        print('Output shape: {}'.format(Y.shape))
        
        fp_array = ( np.asarray((fingerprints), dtype=object) )
        
        # Drop rows from array where FP not generated
        X = np.delete(fp_array, not_found, axis=0)
        
        X = np.vstack(X).astype(np.float32)                 
        
        print('Input shape: {}'.format(X.shape))
        
        #Concatenating input and output array
        final_array = np.concatenate((X, Y), axis=1)
        
        # Removing rows, from final_array, where duplicate FPs are present
        final_array_slice = final_array[:, 0:(final_array.shape[1]-1)]
        _, unq_row_indices = np.unique(final_array_slice,return_index=True,axis=0)
        final_array_unique = final_array[unq_row_indices]
        
        final_numpy_array = np.asarray((final_array_unique), dtype = np.float32)
        print('Number of Duplicate FPs: {}'.format(final_array.shape[0] - final_array_unique.shape[0]))
    
        print('Final Numpy array shape: {}'.format(final_array_unique.shape))
        return final_numpy_array  

def main():
    numpy_file = fg.tpatf_fp()
    np.save(numpy_folder+'/'+filename+'.npy',numpy_file)

if __name__ == "__main__":
    
    numpy_folder = '../data/ligand_numpy'
    
    # Creates folder if it does not exists
    if not os.path.isdir(numpy_folder): os.mkdir(numpy_folder)
    
    files = glob.glob('../data/test_data/*.csv')
    for f in files:
        filename = os.path.split(f)[1][:-4]
        fg = Ligand_Features(f)
        main()
