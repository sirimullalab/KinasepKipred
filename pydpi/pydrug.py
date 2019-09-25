# -*- coding: utf-8 -*-
"""
##############################################################################

A class used for computing different types of drug descriptors! 

You can freely use and distribute it. If you have any problem, 

you could contact with us timely.

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2012.09.24

Email: oriental-cds@163.com

##############################################################################
"""

from .drug import getmol
from .drug import kappa
from .drug import charge
from .drug import connectivity 
from .drug import constitution 
from .drug import estate
from .drug import geary
from .drug import moe
from .drug import molproperty 
from .drug import moran
from .drug import moreaubroto 
from .drug import topology
from .drug import fingerprint
from .drug import basak
from .drug import bcut

from rdkit import Chem
import string

Version=1.0
FingerprintName=['topological','Estate','FP4','atompairs','torsions',
                     'morgan','MACCS']
############################################################################## 

class PyDrug:
    

    """
    #################################################################
    A PyDrug class used for computing drug descriptors.
    #################################################################
    """
    def __init__(self):
        """
        #################################################################
        constructor of pydrug.
        #################################################################
        """
        pass
    
    
    def ReadMolFromMOL(self,filename=""):
        """
        #################################################################
        Read a molecule by SDF or MOL file format.
        
        Usage:
            
            res=ReadMolFromFile(filename)
            
            Input: filename is a file name.
            
            Output: res is a molecule object.
        #################################################################
        """
        self.mol=Chem.MolFromMolMOL(filename)
        return self.mol
    
    
    def ReadMolFromSmile(self,smi=""):
        """
        #################################################################
        Read a molecule by SMILES string.
        
        Usage:
            
            res=ReadMolFromSmile(smi)
            
            Input: smi is a SMILES string.
            
            Output: res is a molecule object.
        #################################################################
        """
        self.mol = Chem.MolFromSmiles(string.strip(smi))
        
        return self.mol
        
        
    def ReadMolFromInchi(self,inchi=""):
        """
        #################################################################
        Read a molecule by Inchi string.
        
        Usage:
            
            res=ReadMolFromInchi(inchi)
            
            Input: inchi is a InChi string.
            
            Output: res is a molecule object.
        #################################################################
        """
        import pybel
        temp=pybel.readstring("inchi",inchi)
        smi=temp.write("smi")
        self.mol = Chem.MolFromSmiles(string.strip(smi))
        
        return self.mol
 
       
    def ReadMolFromMol(self,filename=""):
        """
        #################################################################
        Read a molecule with mol file format.
        
        Usage:
            
            res=ReadMolFromMol(filename)
            
            Input: filename is a file name.
            
            Output: res is a molecule object.
        #################################################################
        """
        self.mol=Chem.MolFromMolFile(filename)
        return self.mol
 
  
   
    def GetMolFromNCBI(self,ID=""):
        """
        #################################################################
        Get a molecule by NCBI id (e.g., 2244).
        
        Usage:
            
            res=GetMolFromNCBI(ID)
            
            Input: ID is a compound ID (CID) in NCBI.
            
            Output: res is a SMILES string.
        #################################################################
        """
        res=getmol.GetMolFromNCBI(cid=ID)
        return res
 
 
   
    def GetMolFromEBI(self,ID=""):
        """
        #################################################################
        Get a molecule by EBI id.

        Usage:
            
            res=GetMolFromEBI(ID)
            
            Input: ID is a compound identifier in EBI.
            
            Output: res is a SMILES string.
        #################################################################
        """
        res=getmol.GetMolFromEBI(ID)
        return res
 
 
   
    def GetMolFromCAS(self,ID=""):
        """
        #################################################################
        Get a molecule by kegg id (e.g., 50-29-3).
        
        Usage:
            
            res=GetMolFromCAS(ID)
            
            Input: ID is a CAS identifier.
            
            Output: res is a SMILES string.
        #################################################################
        """
        res=getmol.GetMolFromCAS(casid=ID)
        return res

   
     
    def GetMolFromKegg(self,ID=""):
        """
        #################################################################
        Get a molecule by kegg id (e.g., D02176).
        
        Usage:
            
            res=GetMolFromKegg(ID)
            
            Input: ID is a compound identifier in KEGG.
            
            Output: res is a SMILES string.
        #################################################################
        """
        res=getmol.GetMolFromKegg(kid=ID)
        return res

 
 
    def GetMolFromDrugbank(self,ID=""):
        """
        #################################################################
        Get a molecule by drugbank id (e.g.,DB00133).
        
        Usage:
            
            res=GetMolFromDrugbank(ID)
            
            Input: ID is a compound identifier in Drugbank.
            
            Output: res is a SMILES string.
        #################################################################
        """
        res=getmol.GetMolFromDrugbank(dbid=ID)
        return res
  
    
    def GetKappa(self):
        """
        #################################################################
        Calculate all kappa descriptors (7).
        
        Usage:
            
            res=GetKappa()
            
            res is a dict form.
        #################################################################
        """
        res=kappa.GetKappa(self.mol)
        return res
    
    
    def GetCharge(self):
        """
        #################################################################
        Calculate all charge descriptors (25).
        
        Usage:
            
            res=GetCharge()
            
            res is a dict form.
        #################################################################
        """
        res=charge.GetCharge(self.mol)
        return res
    
    
    def GetConnectivity(self):
        """
        #################################################################
        Calculate all conenctivity descriptors (44).
        
        Usage:
            
            res=GetConnectivity()
            
            res is a dict form.
        #################################################################
        """
        res=connectivity.GetConnectivity(self.mol)
        return res
        
    
    
    def GetConstitution(self):
        """
        #################################################################
        Calculate all constitutional descriptors (30).
        
        Usage:
            
            res=GetConstitution()
            
            res is a dict form.
        #################################################################
        """
        res=constitution.GetConstitutional(self.mol)
        return res
        
        
    def GetBasak(self):
        """
        #################################################################
        Calculate all basak's information content descriptors (21).
        
        Usage:
            
            res=GetBasak()
            
            res is a dict form.
        #################################################################
        """
        res=basak.Getbasak(self.mol)
        return res    


    def GetBurden(self):
        
        """
        #################################################################
        Calculate all Burden descriptors (64).
        
        Usage:
            
            res=GetBurden()
            
            res is a dict form.
        #################################################################
        """
        res=bcut.GetBurden(self.mol)
        return res        
        
    
    def GetEstate(self):
        """
        #################################################################
        Calculate estate descriptors (316).
        
        Usage:
            
            res=GetEstate()
            
            res is a dict form.
        #################################################################
        """
        res=estate._GetEstate(self.mol)
        return res
    
    
    def GetGeary(self):
        """
        #################################################################
        Calculate all Geary autocorrelation descriptors (32).
        
        Usage:
            
            res=GetGeary()
            
            res is a dict form.
        #################################################################
        """
        res=geary.GetGearyAuto(self.mol)
        return res
    
    
    def GetMOE(self):
        """
        #################################################################
        Calculate all MOE-type descriptors (60).
        
        Usage:
            
            res=GetMOE()
            
            res is a dict form.
        #################################################################
        """
        res=moe.GetMOE(self.mol)
        return res
    
    
    def GetMolProperty(self):
        """
        #################################################################
        Calculate all molecular properties (6).
        
        Usage:
            
            res=GetMolProperty()
            
            res is a dict form.
        #################################################################
        """
        res=molproperty.GetMolecularProperty(self.mol)
        return res
    
    
    def GetMoran(self):
        """
        #################################################################
        Calculate all Moran autocorrealtion descriptors (32).
        
        Usage:
#            
            res=GetMoran()
            
            res is a dict form.
        #################################################################
        """
        res=moran.GetMoranAuto(self.mol)
        return res
    
    
    def GetMoreauBroto(self):
        """
        #################################################################
        Calculate all Moreau-Broto autocorrelation descriptors(32).
        
        Usage:
            
            res=GetMoreauBroto()
            
            res is a dict form.
        #################################################################
        """
        res=moreaubroto.GetMoreauBrotoAuto(self.mol)
        return res
    
    
    def GetTopology(self):
        """
        #################################################################
        Calculate all topological descriptors (25).
        
        Usage:
            
            res=GetTopology()
            
            res is a dict form.
        #################################################################
        """
        res=topology.GetTopology(self.mol)
        return res
    
    
    def GetFingerprint(self,FPName='topological'):
        """
        #################################################################
        Calculate all fingerprint descriptors.
        
        see the fingerprint type in FingerprintName
        
        Usage:
            
            res=GetFingerprint(FPName='topological')
            
            res is a tuple form.
        #################################################################
        """
        
        if FPName in FingerprintName:
            temp=fingerprint._FingerprintFuncs[FPName]
            res=temp(self.mol)
        else:
            res=fingerprint.CalculateDaylightFingerprint(self.mol)
        
        
        return res[0],res[1]
            
        
    def GetAllDescriptor(self):
        """
        #################################################################
        Calculate all descriptors (608).
        
        Usage:
            
            res=GetAllDescriptor()
            
            res is a dict form.
        #################################################################
        """
        res={}
        res.update(self.GetKappa())
        res.update(self.GetCharge())
        res.update(self.GetConnectivity())
        res.update(self.GetConstitution())
        res.update(self.GetEstate())
        res.update(self.GetGeary())
        res.update(self.GetMOE())
        res.update(self.GetMoran())
        res.update(self.GetMoreauBroto())
        res.update(self.GetTopology())
        res.update(self.GetMolProperty())
        res.update(self.GetBasak())
        res.update(self.GetBurden())
        
        return res
##############################################################################    

if __name__=="__main__":
    
    drugclass=PyDrug()
    drugclass.ReadMolFromSmile("CCC1(c2ccccc2)C(=O)N(C)C(=N1)O")
    print(drugclass.GetKappa())
    print(len(drugclass.GetKappa()))
    print(drugclass.GetTopology())
    print(len(drugclass.GetTopology()))
    print(drugclass.GetMoreauBroto())
    res=drugclass.GetAllDescriptor()
    #print drugclass.GetMolFromDrugbank(ID="DB00133")
    #res=drugclass.GetFingerprint(FPName='Estate')
    print(res)
    print(len(res))
