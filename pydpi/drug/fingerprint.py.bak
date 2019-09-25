# -*- coding: utf-8 -*-
"""
##############################################################################
This module is to compute the various fingerprints  based on the provided 

fingerprint system. If you have any question please contact me via email.

My email adress is : orientalcds@gmail.com

2012.09.25

@author: Dongsheng Cao
##############################################################################
"""
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from rdkit import Chem
from estate import CalculateEstateFingerprint as EstateFingerprint
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit import DataStructs


Version=1.0
similaritymeasure=[i[0] for i in DataStructs.similarityFunctions]
################################################################

def CalculateDaylightFingerprint(mol):
    """
    #################################################################
    Calculate Daylight-like fingerprint or topological fingerprint
    
    (2048 bits).
    
    Usage:
        
        result=CalculateDaylightFingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res={}
    NumFinger=2048
    bv=FingerprintMols.FingerprintMol(mol)
    temp=tuple(bv.GetOnBits())
    for i in temp:
        res.update({i:1})

    return NumFinger,res,bv



def CalculateMACCSFingerprint(mol):
    """
    #################################################################
    Calculate MACCS keys (166 bits).
    
    Usage:
        
        result=CalculateMACCSFingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res={}
    NumFinger=166
    bv=MACCSkeys.GenMACCSKeys(mol)
    temp=tuple(bv.GetOnBits())
    for i in temp:
        res.update({i:1})

    return NumFinger,res,bv
            

def CalculateFP4Fingerprint(mol):
    """
    #################################################################
    Calculate FP4 fingerprints (307 bits).
    
    Usage:
        
        result=CalculateFP4Fingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    import pybel
    res={}
    NumFinger=307
    m=pybel.readstring('smi',mol)
    temp=m.calcfp('FP4').bits
    for i in temp:
        res.update({i:1})

    return NumFinger,res
    


def CalculateEstateFingerprint(mol):
    """
    #################################################################
    Calculate E-state fingerprints (79 bits).
    
    Usage:
        
        result=CalculateEstateFingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    NumFinger=79
    res={}
    temp=EstateFingerprint(mol)
    for i in temp:
        if temp[i]>0:
            res[i[7:]]=1
       
    return NumFinger,res,temp


def CalculateAtomPairsFingerprint(mol):
    """
    #################################################################
    Calculate atom pairs fingerprints
    
    Usage:
        
        result=CalculateAtomPairsFingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    
    res=Pairs.GetAtomPairFingerprint(mol)
    
    return res.GetLength(),res.GetNonzeroElements(),res



def CalculateTopologicalTorsionFingerprint(mol):
    """
    #################################################################
    Calculate Topological Torsion Fingerprints
    
    Usage:
        
        result=CalculateTopologicalTorsionFingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res=Torsions.GetTopologicalTorsionFingerprint(mol)

    return res.GetLength(),res.GetNonzeroElements(),res



def CalculateMorganFingerprint(mol,radius=2):
    """
    #################################################################
    Calculate Morgan
    
    Usage:
        
        result=CalculateMorganFingerprint(mol)
        
        Input: mol is a molecule object.
        
        radius is a radius.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res=AllChem.GetMorganFingerprint(mol,radius)
    
    return res.GetLength(),res.GetNonzeroElements(),res


def CalculateSimilarity(fp1,fp2,similarity="Tanimoto"):
    """
    #################################################################
    Calculate similarity between two molecules.
    
    Usage:
        
        result=CalculateSimilarity(fp1,fp2)
        
        Input: fp1 and fp2 are two DataStructs.
        
        Output: result is a similarity value.
    #################################################################
    """
    temp=DataStructs.similarityFunctions
    for i in temp:
        if similarity in i[0]:
            similarityfunction=i[1]
        else:
            similarityfunction=temp[0][1]
            
    res=similarityfunction(fp1,fp2)
    return round(res,3)


_FingerprintFuncs={'topological':CalculateDaylightFingerprint,
                 'Estate':CalculateEstateFingerprint,
                 'FP4':CalculateFP4Fingerprint,
                 'atompairs':CalculateAtomPairsFingerprint,
                 'torsions':CalculateTopologicalTorsionFingerprint,
                 'morgan':CalculateMorganFingerprint,
                 'MACCS':CalculateMACCSFingerprint}
################################################################


if __name__=="__main__":
    
    print "fingerprint......"
    
    ms = [Chem.MolFromSmiles('CCOC=N'), Chem.MolFromSmiles('CCO')]
    res1=CalculateTopologicalTorsionFingerprint(ms[0])
    print res1
    res2=CalculateTopologicalTorsionFingerprint(ms[1])
    print res2
    print CalculateSimilarity(res1[2],res2[2])
    