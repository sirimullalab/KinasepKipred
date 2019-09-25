# -*- coding: utf-8 -*-
"""
##############################################################################

This script is used for providing various application examples guiding how

the user use our provided functionality to calculate different types of 

molecular features. You can freely use and distribute it. If you hava 

any problem, you could contact with us timely!

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2012.09.24

Email: oriental-cds@163.com 

##############################################################################
"""

import pydrug, pypro, pydpi

def main():
    """
    An example
    You must connect the internet to run the example.
    """
    smi=pydrug.getmol.GetMolFromNCBI('2244')
    print smi
    drugobject=pydrug.PyDrug()
    print drugobject
    mol=drugobject.ReadMolFromSmile(smi)
    print mol
    print 'calculating constitutional descsriptors'
    print drugobject.GetConstitution()
    print 'calculating topology descriptors'
    print drugobject.GetTopology()
    print 'calculating connectivity indices'
    print drugobject.GetConnectivity()
    print 'calculating E-state descriptors'
    print drugobject.GetEstate() 
    print 'calculating Basak descriptors'
    print drugobject.GetBasak()
    print 'calculating Burden descriptors'
    print drugobject.GetBurden()
    print 'calculating kappa shape descriptors'
    print drugobject.GetKappa()
    print 'calculating MOE-type descriptors'
    print drugobject.GetMOE()
    print 'calculating Geary autocorrelation descriptors'
    print drugobject.GetGeary()
    print 'calculating Moran autocorrelation descriptors'
    print drugobject.GetMoran()
    print 'calculating Moreau-Broto autocorrelation descriptors'
    print drugobject.GetMoreauBroto()
    print 'calculating charge descriptors'
    print drugobject.GetCharge()
    print 'calculating molecular properties'
    print drugobject.GetMolProperty()
    print 'calculating all descriptors'
    print drugobject.GetAllDescriptor()
    
    proteinobject=pypro.PyPro()
    ps=proteinobject.GetProteinSequenceFromID('P48039')
    print ps
    proteinobject.ReadProteinSequence(ps)
    print 'calculating amino acid composition'
    print proteinobject.GetAAComp()
    print 'calculating CTD descriptors'
    print proteinobject.GetCTD()
    print 'calculating dipeptide composition'
    print proteinobject.GetDPComp()
    print 'calculating Moran autocorrelation descriptors'
    print proteinobject.GetMoranAuto()
    print 'calculating APPAC descriptors'
    print proteinobject.GetAPAAC()
    print 'calculating QSO descriptors'
    print proteinobject.GetQSO()
    print 'calculating conjoint triad features'
    print proteinobject.GetTriad()
    print 'calculating SOCN descriptors'
    print proteinobject.GetSOCN()
    print 'calculating tripeptide composition'
    print proteinobject.GetTPComp()
    print 'calcualing pseudo amino acid composition'
    print proteinobject.GetPAAC()
    
    dpi=pydpi.PyDPI()
    smi=dpi.GetMolFromNCBI('2244')
    print smi
    dpi.ReadMolFromSmile(smi)
    ddict=dpi.GetConnectivity()
    print ddict
    ddict.update(dpi.GetKappa())
    print ddict
    ps=dpi.GetProteinSequenceFromID('P48039')
    print ps
    dpi.ReadProteinSequence(ps)
    pdict=dpi.GetAAComp()
    print pdict
    pdict.update(dpi.GetAPAAC())
    print pdict
    print dpi.GetDPIFeature1(ddict,pdict)
    print dpi.GetDPIFeature2(ddict,pdict)
    
    
    ppi=pydpi.PyPPI()
    ps=ppi.GetProteinSequenceFromID('P48039')
    print ps
    ppi.ReadProteinSequence(ps)
    
    pdict1=ppi.GetPAAC(lamda=10,weight=0.05)
    print pdict1
    
    ppi=pydpi.PyPPI()
    ps=ppi.GetProteinSequenceFromID('P08172')
    print ps
    ppi.ReadProteinSequence(ps)
    
    pdict2=ppi.GetPAAC(lamda=10,weight=0.05)
    
    print pdict2
    
    print ppi.GetPPIFeature1(pdict1,pdict2)
    print ppi.GetPPIFeature2(pdict1,pdict2)



if __name__=="__main__":

    main()






