# -*- coding: utf-8 -*-
"""
##############################################################################

A class used for computing drug-target interaction, protein-protein interaction

features. The PyDPI class inherits the PyDrug and PyPro classes, So you can easily

use them. You can freely use and distribute it. If you hava  any problem, 

you could contact with us timely!

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2012.09.24

Email: oriental-cds@163.com

##############################################################################
"""

from pydrug import *
from pypro import *

Version=1.0
#############################################################################
class PyDPI(PyDrug,PyPro):
    """
    #################################################################
    
    A PyDPI class used for generating drug-target interaction features.
    
    #################################################################
    """
    def __init__(self):
        """
        #################################################################
        
        constructor of PyDPI.
        
        #################################################################
        """
        pass
    
    def GetDPIFeature1(self,ddict={},pdict={}):
        """
        #################################################################
        
        Calculate the drug-target interaction features by combining drug 
        
        features and protein features.(nd+np)
        
        Usage:
            
            res=GetDPIFeature1(ddict,pdict)
            
            Input: ddict is a dict form containing drug features.
            
                   pdict is a dict form containing protein features.
                   
            Output: res is a dict form containing drug-target interaction
            
            features.
            
        #################################################################
        """
        result={}
        result.update(ddict)
        result.update(pdict)
        return result

    def GetDPIFeature2(self,ddict={},pdict={}):
        """
        #################################################################
        Calculate the drug-target interaction features by  the tensor product.
        
        (nd*np)
        
        Usage:
            
            res=GetDPIFeature2(ddict,pdict)
            
            Input: ddict is a dict form containing drug features.
            
                   pdict is a dict form containing protein features.
                   
            Output: res is a dict form containing drug-target interaction
            
            features.
        #################################################################
        """
        res={}
        for i in ddict:
            for j in pdict:
                res[i+'*'+j]=round(ddict[i]*pdict[j],3)
        return res

######################################################################

class PyPPI(PyPro):
    """
    #################################################################
    
    A PyPPI class used for generating protein-protein interaction features.
    
    #################################################################
    """
    
    def __init__(self):
        """
        #################################################################
        
        constructor of PyPPI.
        
        #################################################################
        """
        pass  
 
    def GetPPIFeature1(self,pdict1={},pdict2={}):
        """
        #################################################################
        Calculate the protein-protein interaction features by 
        
        F=[Fa(i)+Fb(i)),Fa(i)*Fb(i)] (2n)

        Usage:
            
            res=GetPPIFeature1(pdict)
            
            Input: pdict1 and pdict2 are dict forms containing protein features.
                   
            Output: res is a dict form containing protein-protein interaction
            
            features.
        #################################################################
        """
        res={}
        for i in pdict1:
            res[i+'-'+i]=round(pdict1[i]*pdict2[i],3)
            res[i+'*'+i]=round(pdict1[i]*pdict2[i],3)
        return res
    
    def GetPPIFeature2(self,pdict1={},pdict2={}):
        """
        #################################################################
        Calculate the protein-protein interaction features by the tensor product.

        (n^2)

        Usage:
            
            res=GetPPIFeature2(pdict)
            
            Input: pdict1 and pdict2 are dict forms containing protein features.
                   
            Output: res is a dict form containing protein-protein interaction
            
            features.
        #################################################################
        """
        res={}
        for i in pdict1:
            for j in pdict2:
                res[i+'*'+j]=round(pdict1[i]*pdict2[j],3)
        return res
        
#############################################################################

if __name__=="__main__":

    print "calculating......"
    
    cds=PyDPI()
    smi=cds.GetMolFromNCBI("2244") ##connect internet
    cds.ReadMolFromSmile(smi)
    drugdes=cds.GetConnectivity()
    print drugdes
    
    ps=cds.GetProteinSequenceFromID("P48039") ##connect internet
    cds.ReadProteinSequence(ps)
    prodes=cds.GetAPAAC()
    print prodes
    
    dpi1=cds.GetDPIFeature1(drugdes,prodes)
    print len(dpi1)
    
    dpi2=cds.GetDPIFeature2(drugdes,prodes)
    print len(dpi2)
    
    temp=PyPPI()
    ps=temp.GetProteinSequenceFromID("P48039") ##connect internet
    temp.ReadProteinSequence(ps)
    prodes=temp.GetAPAAC()
    print prodes
    ppi1=temp.GetPPIFeature1(prodes)
    print len(ppi1)

    ppi2=temp.GetPPIFeature2(prodes)
    print len(ppi2)

