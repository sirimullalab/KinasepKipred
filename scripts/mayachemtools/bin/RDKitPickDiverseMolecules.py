#!/bin/env python
#
# File: RDKitPickDiverseMolecules.py
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2018 Manish Sud. All rights reserved.
#
# The functionality available in this script is implemented using RDKit, an
# open source toolkit for cheminformatics developed by Greg Landrum.
#
# This file is part of MayaChemTools.
#
# MayaChemTools is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# MayaChemTools is distributed in the hope that it will be useful, but without
# any warranty; without even the implied warranty of merchantability of fitness
# for a particular purpose.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MayaChemTools; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation Inc., 59 Temple Place, Suite 330,
# Boston, MA, 02111-1307, USA.
#

from __future__ import print_function

# Add local python path to the global path and import standard library modules...
import os
import sys;  sys.path.insert(0, os.path.join(os.path.dirname(sys.argv[0]), "..", "lib", "Python"))
import time
import re

# RDKit imports...
try:
    from rdkit import rdBase
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import DataStructs
    from rdkit.Chem.Fingerprints import FingerprintMols
    from rdkit.Chem import rdMolDescriptors
    from rdkit.SimDivFilters import rdSimDivPickers
    from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
    from rdkit.SimDivFilters.rdSimDivPickers import HierarchicalClusterPicker
except ImportError as ErrMsg:
    sys.stderr.write("\nFailed to import RDKit module/package: %s\n" % ErrMsg)
    sys.stderr.write("Check/update your RDKit environment and try again.\n\n")
    sys.exit(1)

# MayaChemTools imports...
try:
    from docopt import docopt
    import MiscUtil
    import RDKitUtil
except ImportError as ErrMsg:
    sys.stderr.write("\nFailed to import MayaChemTools module/package: %s\n" % ErrMsg)
    sys.stderr.write("Check/update your MayaChemTools environment and try again.\n\n")
    sys.exit(1)

ScriptName = os.path.basename(sys.argv[0])
Options = {}
OptionsInfo = {}

def main():
    """Start execution of the script"""
    
    MiscUtil.PrintInfo("\n%s (RDK v%s; %s): Starting...\n" % (ScriptName, rdBase.rdkitVersion, time.asctime()))
    
    (WallClockTime, ProcessorTime) = MiscUtil.GetWallClockAndProcessorTime()
    
    # Retrieve command line arguments and options...
    RetrieveOptions()
    
    # Process and validate command line arguments and options...
    ProcessOptions()
    
    # Perform actions required by the script...
    PickDiverseMolecules()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PickDiverseMolecules():
    """Pick diverse molecules."""

    Mols = RetrieveMolecules()
    MolsFingerprints = GenerateFingerprints(Mols)
    DiverseMols = SelectMolecules(Mols, MolsFingerprints)
    
    WriteMolecules(DiverseMols)

def SelectMolecules(Mols, MolsFingerprints):
    """Select diverse molecules."""

    if OptionsInfo["NumMols"] > len(Mols):
        MiscUtil.PrintError("The number of diverse molecules to pick, %d, specified using \"-n, --numMols\" must be less than total number of valid molecules, %d" % (OptionsInfo["NumMols"], len(Mols)))
    
    DiverseMols = []
    if re.match("^MaxMin$", OptionsInfo["Mode"], re.I):
        return SelectMoleculesUsingMaxMin(Mols, MolsFingerprints)
    elif re.match("^HierarchicalClustering$", OptionsInfo["Mode"], re.I):
        return SelectMoleculesUsingHierarchicalClustering(Mols, MolsFingerprints)
    else:
        MiscUtil.PrintError("The mode vaue, %s, is not a valid mode." % OptionsInfo["Mode"])
    
    return DiverseMols

def SelectMoleculesUsingMaxMin(Mols, MolsFingerprints):
    """Select diverse molecules using MaxMin methodology."""

    MiscUtil.PrintInfo("\nSelecting diverse molecules using MaxMin methodology and %s similarity metric..." % OptionsInfo["SimilarityMetric"])
    
    DiverseMols = []
    
    PoolSize = len(MolsFingerprints)
    PickSize = OptionsInfo["NumMols"]
    SimilarityFunction = OptionsInfo["SimilarityFunction"]

    Picker = MaxMinPicker()
    PairwiseDistance = lambda i, j: 1 - SimilarityFunction(MolsFingerprints[i], MolsFingerprints[j])

    MolIndices = Picker.LazyPick(PairwiseDistance, PoolSize, PickSize)
            
    for Index in list(MolIndices):
        DiverseMols.append(Mols[Index])
    
    return DiverseMols

def SelectMoleculesUsingHierarchicalClustering(Mols, MolsFingerprints):
    """Select diverse molecules using hierarchical clustering  methodology."""

    try:
        import numpy
    except ImportError:
        MiscUtil.PrintError("Failed to import numpy python module. This is required for picking diverse molecules using hierarchical for clustering.")
    
    MiscUtil.PrintInfo("\nSelecting diverse molecules using %s hierarchical clustering methodology..." % OptionsInfo["SpecifiedClusteringMethod"])
    
    DiverseMols = []
    
    PoolSize = len(MolsFingerprints)
    PickSize = OptionsInfo["NumMols"]
    DistanceMatrix = GenerateLowerTriangularDistanceMatrix(MolsFingerprints)
    
    ClusterPicker = HierarchicalClusterPicker(OptionsInfo["SpecifiedClusteringMethodID"])
    MolIndices = ClusterPicker.Pick(numpy.asarray(DistanceMatrix), PoolSize, PickSize)
    
    for Index in MolIndices:
        DiverseMols.append(Mols[Index])
    
    return DiverseMols

def RetrieveMolecules():
    """Retrieve molecules."""

    Infile = OptionsInfo["Infile"]
    
    # Read molecules...
    MiscUtil.PrintInfo("\nReading file %s..." % Infile)
    
    OptionsInfo["InfileParams"]["AllowEmptyMols"] = False
    ValidMols, MolCount, ValidMolCount  = RDKitUtil.ReadAndValidateMolecules(Infile, **OptionsInfo["InfileParams"])
    
    MiscUtil.PrintInfo("Total number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))

    return ValidMols

def GenerateFingerprints(Mols):
    """Generate fingerprints."""

    FingerprintsName = OptionsInfo["SpecifiedFingerprints"]
    
    MolsFingerprints = []
    if re.match("^AtomPairs$", FingerprintsName, re.I):
        return GenerateAtomPairsFingerprints(Mols)
    elif re.match("^MACCS166Keys$", FingerprintsName, re.I):
        return GenerateMACCS166KeysFingerprints(Mols)
    elif re.match("^Morgan$", FingerprintsName, re.I):
        return GenerateMorganFingerprints(Mols)
    elif re.match("^MorganFeatures$", FingerprintsName, re.I):
        return GenerateMorganFeaturesFingerprints(Mols)
    elif re.match("^PathLength$", FingerprintsName, re.I):
        return GeneratePathLengthFingerprints(Mols)
    elif re.match("^TopologicalTorsions$", FingerprintsName, re.I):
        return GenerateTopologicalTorsionsFingerprints(Mols)
    else:
        MiscUtil.PrintError("Fingerprints name, %s, is not a valid name" % FingerprintsName)
    
    return MolsFingerprints

def GenerateAtomPairsFingerprints(Mols):
    """Generate AtomPairs fingerprints."""

    MiscUtil.PrintInfo("\nGenerating AtomPairs fingerprints...")
    
    MinLength = OptionsInfo["FingerprintsParams"]["AtomPairs"]["MinLength"]
    MaxLength = OptionsInfo["FingerprintsParams"]["AtomPairs"]["MaxLength"]
    UseChirality = OptionsInfo["FingerprintsParams"]["AtomPairs"]["UseChirality"]

    if OptionsInfo["GenerateBitVectFingerints"]:
        # Generate ExplicitBitVect fingerprints...
        FPSize = 2048
        BitsPerHash = 4
        MolsFingerprints = [rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Mol, minLength = MinLength, maxLength = MaxLength, includeChirality = UseChirality, nBits = FPSize, nBitsPerEntry = BitsPerHash) for Mol in Mols]
    else:
        # Generate IntSparseIntVect fingerprints...
        MolsFingerprints = [rdMolDescriptors.GetAtomPairFingerprint(Mol, minLength = MinLength, maxLength = MaxLength, includeChirality = UseChirality) for Mol in Mols]

    return MolsFingerprints

def GenerateMACCS166KeysFingerprints(Mols):
    """Generate MACCS166Keys fingerprints."""

    MiscUtil.PrintInfo("\nGenerating MACCS166Keys fingerprints...")

    # Generate ExplicitBitVect fingerprints...
    MolsFingerprints = [rdMolDescriptors.GetMACCSKeysFingerprint(Mol) for Mol in Mols]

    return MolsFingerprints

def GenerateMorganFingerprints(Mols):
    """Generate Morgan fingerprints."""

    MiscUtil.PrintInfo("\nGenerating  Morgan fingerprints...")
    
    Radius = OptionsInfo["FingerprintsParams"]["Morgan"]["Radius"]
    UseChirality = OptionsInfo["FingerprintsParams"]["Morgan"]["UseChirality"]
    UseFeatures = False

    if OptionsInfo["GenerateBitVectFingerints"]:
        # Generate ExplicitBitVect fingerprints...
        FPSize = 2048
        MolsFingerprints = [rdMolDescriptors.GetMorganFingerprintAsBitVect(Mol, Radius, useFeatures = UseFeatures, useChirality = UseChirality, nBits = FPSize) for Mol in Mols]
    else:
        # Generate UIntSparseIntVect fingerprints...
        MolsFingerprints = [rdMolDescriptors.GetMorganFingerprint(Mol, Radius, useFeatures = UseFeatures, useChirality = UseChirality) for Mol in Mols]

    return MolsFingerprints

def GenerateMorganFeaturesFingerprints(Mols):
    """Generate MorganFeatures fingerprints."""

    MiscUtil.PrintInfo("\nGenerating  MorganFeatures fingerprints...")
    
    # Setup fingerprints parameters...
    Radius = OptionsInfo["FingerprintsParams"]["MorganFeatures"]["Radius"]
    UseChirality = OptionsInfo["FingerprintsParams"]["MorganFeatures"]["UseChirality"]
    UseFeatures = True
    
    if OptionsInfo["GenerateBitVectFingerints"]:
        # Generate ExplicitBitVect fingerprints...
        FPSize = 2048
        MolsFingerprints = [rdMolDescriptors.GetMorganFingerprintAsBitVect(Mol, Radius, useFeatures = UseFeatures, useChirality = UseChirality, nBits = FPSize) for Mol in Mols]
    else:
        # Generate UIntSparseIntVect fingerprints...
        MolsFingerprints = [rdMolDescriptors.GetMorganFingerprint(Mol, Radius, useFeatures = UseFeatures, useChirality = UseChirality) for Mol in Mols]

    return MolsFingerprints

def GeneratePathLengthFingerprints(Mols):
    """Generate PathLength fingerprints."""

    MiscUtil.PrintInfo("\nGenerating PathLength fingerprints ...")
    
    MinPath = OptionsInfo["FingerprintsParams"]["PathLength"]["MinPath"]
    MaxPath = OptionsInfo["FingerprintsParams"]["PathLength"]["MaxPath"]
    FPSize = OptionsInfo["FingerprintsParams"]["PathLength"]["FPSize"]
    BitsPerHash = OptionsInfo["FingerprintsParams"]["PathLength"]["BitsPerHash"]
    UseHs = False
    TargetDensity = 0.3
    MinSize = 54

    # Generate ExplicitBitVect fingerprints...
    MolsFingerprints = [FingerprintMols.FingerprintMol(Mol, minPath = MinPath, maxPath = MaxPath, fpSize = FPSize, bitsPerHash = BitsPerHash, useHs = UseHs, tgtDensity = TargetDensity, minSize = MinSize) for Mol in Mols]

    return MolsFingerprints

def GenerateTopologicalTorsionsFingerprints(Mols):
    """Generate TopologicalTorsions fingerprints."""

    MiscUtil.PrintInfo("\nGenerating TopologicalTorsions fingerprints...")
    
    UseChirality = OptionsInfo["FingerprintsParams"]["TopologicalTorsions"]["UseChirality"]

    if OptionsInfo["GenerateBitVectFingerints"]:
        FPSize = 2048
        BitsPerHash = 4
        MolsFingerprints = [rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(Mol,  includeChirality = UseChirality, nBits = FPSize, nBitsPerEntry = BitsPerHash) for Mol in Mols]
    else:
        # Generate LongSparseIntVect fingerprint...
        MolsFingerprints = [rdMolDescriptors.GetTopologicalTorsionFingerprint(Mol,  includeChirality = UseChirality) for Mol in Mols]

    return MolsFingerprints

def GenerateLowerTriangularDistanceMatrix(MolsFingerprints):
    """Generate a lower triangular distance matrix without the diagonal."""

    SimilarityFunction = OptionsInfo["SimilarityFunction"]

    DistanceMatrix = []
    NumFPs = len(MolsFingerprints)
    for Index1 in range(0, NumFPs):
        for Index2 in range(0, Index1):
            Distance =  1 - SimilarityFunction(MolsFingerprints[Index1], MolsFingerprints[Index2],)
            DistanceMatrix.append(Distance)

    return DistanceMatrix

def WriteMolecules(Mols):
    """Write out molecules."""

    Outfile = OptionsInfo["Outfile"]
    
    # Set up a molecule writer...
    Writer = None
    Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
    MiscUtil.PrintInfo("\nGenerating file %s...\n" % Outfile)

    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    
    # Write out molecules...
    for Mol in Mols:
        if Compute2DCoords:
            AllChem.Compute2DCoords(Mol)
        Writer.write(Mol)
        
    if Writer is not None:
        Writer.close()
        
    MiscUtil.PrintInfo("Total number of diverse molecules selected: %d" % (len(Mols)))

def ProcessFingerprintsParameters():
    """Set up and process fingerprints parameters."""

    SetupFingerprintsNamesAndParameters()
    ProcessSpecifiedFingerprintsName()
    ProcessSpecifiedFingerprintsParameters()

def SetupFingerprintsNamesAndParameters():
    """Set up fingerprints parameters."""
    
    OptionsInfo["FingerprintsNames"] = ["AtomPairs", "MACCS166Keys", "Morgan", "MorganFeatures", "PathLength", "TopologicalTorsions"]
    
    OptionsInfo["FingerprintsParams"] = {}
    OptionsInfo["FingerprintsParams"]["AtomPairs"] = {"MinLength": 1, "MaxLength": 30, "UseChirality": False}
    OptionsInfo["FingerprintsParams"]["MACCS166Keys"] = {}
    OptionsInfo["FingerprintsParams"]["Morgan"] = {"Radius": 2, "UseChirality": False}
    OptionsInfo["FingerprintsParams"]["MorganFeatures"] = {"Radius": 2, "UseChirality": False}
    OptionsInfo["FingerprintsParams"]["TopologicalTorsions"] = {"UseChirality": False}
    OptionsInfo["FingerprintsParams"]["PathLength"] = {"MinPath": 1, "MaxPath": 7, "FPSize": 2048, "BitsPerHash": 2}

def ProcessSpecifiedFingerprintsName():
    """Process specified fingerprints name."""

    #  Set up a canonical fingerprints name map...
    CanonicalFingerprintsNamesMap = {}
    for Name in OptionsInfo["FingerprintsNames"]:
        CanonicalName = Name.lower()
        CanonicalFingerprintsNamesMap[CanonicalName] = Name

    # Validate specified fingerprints name...
    CanonicalFingerprintsName = OptionsInfo["Fingerprints"].lower()
    if CanonicalFingerprintsName not in CanonicalFingerprintsNamesMap:
        MiscUtil.PrintError("The fingerprints name, %s, specified using \"-f, --fingerprints\" option is not a valid name." % (OptionsInfo["Fingerprints"]))
    
    OptionsInfo["SpecifiedFingerprints"] = CanonicalFingerprintsNamesMap[CanonicalFingerprintsName]

def ProcessSpecifiedFingerprintsParameters():
    """Process specified fingerprints parameters."""

    if re.match("^auto$", OptionsInfo["ParamsFingerprints"], re.I):
        # Nothing to process...
        return

    SpecifiedFingerprintsName = OptionsInfo["SpecifiedFingerprints"]
    
    # Parse specified fingerprints parameters...
    ParamsFingerprints = re.sub(" ", "", OptionsInfo["ParamsFingerprints"])
    if not ParamsFingerprints:
        MiscUtil.PrintError("No valid parameter name and value pairs specified using \"-p, --paramsFingerprints\" option corrresponding to fingerprints %s." % (SpecifiedFingerprintsName))

    ParamsFingerprintsWords = ParamsFingerprints.split(",")
    if len(ParamsFingerprintsWords) % 2:
        MiscUtil.PrintError("The number of comma delimited paramater names and values, %d, specified using \"-p, --paramsFingerprints\" option must be an even number." % (len(ParamsFingerprintsWords)))

    # Setup a canonical parameter names for specified fingerprints...
    ValidParamNames = []
    CanonicalParamNamesMap = {}
    for ParamName in sorted(OptionsInfo["FingerprintsParams"][SpecifiedFingerprintsName]):
        ValidParamNames.append(ParamName)
        CanonicalParamNamesMap[ParamName.lower()] = ParamName

    # Validate and set paramater names and value...
    for Index in range(0, len(ParamsFingerprintsWords), 2):
        Name = ParamsFingerprintsWords[Index]
        Value = ParamsFingerprintsWords[Index + 1]

        CanonicalName = Name.lower()
        if  not CanonicalName in CanonicalParamNamesMap:
            MiscUtil.PrintError("The parameter name, %s, specified using \"-p, --paramsFingerprints\" option for fingerprints, %s, is not a valid name. Supported parameter names: %s" % (Name, SpecifiedFingerprintsName, " ".join(ValidParamNames)))

        ParamName = CanonicalParamNamesMap[CanonicalName]
        if re.match("^UseChirality$", ParamName, re.I):
            if not re.match("^(Yes|No|True|False)$", Value, re.I):
                MiscUtil.PrintError("The parameter value, %s, specified using \"-p, --paramsFingerprints\" option for fingerprints, %s, is not a valid value. Supported values: Yes No True False" % (Value, SpecifiedFingerprintsName))
            ParamValue = False
            if re.match("^(Yes|True)$", Value, re.I):
                ParamValue = True
        else:
            ParamValue = int(Value)
            if ParamValue <= 0:
                MiscUtil.PrintError("The parameter value, %s, specified using \"-p, --paramsFingerprints\" option for fingerprints, %s, is not a valid value. Supported values: > 0" % (Value, SpecifiedFingerprintsName))
        
        # Set value...
        OptionsInfo["FingerprintsParams"][SpecifiedFingerprintsName][ParamName] = ParamValue

def ProcessSimilarityMetricParameter():
    """Process specified similarity metric value."""

    SimilarityInfoMap = {}
    CanonicalNameMap = {}
    
    for SimilarityFunctionInfo in DataStructs.similarityFunctions:
        Name = SimilarityFunctionInfo[0]
        Function = SimilarityFunctionInfo[1]
        
        SimilarityInfoMap[Name] = Function
        CanonicalName = Name.lower()
        CanonicalNameMap[CanonicalName] = Name
    
    SpecifiedCanonicalName = OptionsInfo["SimilarityMetric"].lower()
    SimilarityFunction = None
    if  SpecifiedCanonicalName in CanonicalNameMap:
        SimilarityName = CanonicalNameMap[SpecifiedCanonicalName]
        SimilarityFunction = SimilarityInfoMap[SimilarityName]
    else:
        MiscUtil.PrintError("Similarity metric name, %s, is not a valid name. " % OptionsInfo["SimilarityMetric"])
        
    OptionsInfo["SimilarityMetric"] = SimilarityName
    OptionsInfo["SimilarityFunction"] = SimilarityFunction

    # RDKit similarity functions, besides Dice and Tanimoto, are not able to handle int bit vectors...
    GenerateBitVectFingerints = False
    if not re.match("^(Tanimoto|Dice)$", SimilarityName, re.I):
        GenerateBitVectFingerints = True
    OptionsInfo["GenerateBitVectFingerints"] = GenerateBitVectFingerints
    
def ProcessClusteringMethodParameter():
    """Process specified clustering method parameter."""

    OptionsInfo["SpecifiedClusteringMethod"] = ""
    OptionsInfo["SpecifiedClusteringMethodID"] = ""
    
    if not re.match("^HierarchicalClustering$", OptionsInfo["Mode"], re.I):
        # Nothing to process...
        return

    # Setup a canonical cluster method name map..
    ClusteringMethodInfoMap = {}
    CanonicalClusteringMethodNameMap = {}
    for Name in sorted(rdSimDivPickers.ClusterMethod.names):
        NameID =  rdSimDivPickers.ClusterMethod.names[Name]
        ClusteringMethodInfoMap[Name] = NameID
        
        CanonicalName = Name.lower()
        CanonicalClusteringMethodNameMap[CanonicalName] = Name

    CanonicalName = OptionsInfo["ClusteringMethod"].lower()
    if not CanonicalName in CanonicalClusteringMethodNameMap:
        MiscUtil.PrintError("The cluster method, %s, specified using \"-c, --clusteringMethod\" option is not a valid name." % (OptionsInfo["ClusteringMethod"]))

    SpecifiedClusteringMethodName = CanonicalClusteringMethodNameMap[CanonicalName]
    OptionsInfo["SpecifiedClusteringMethod"] = SpecifiedClusteringMethodName
    OptionsInfo["SpecifiedClusteringMethodID"] = ClusteringMethodInfoMap[SpecifiedClusteringMethodName] 
    
def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["Mode"] = Options["--mode"]
    OptionsInfo["Fingerprints"] = Options["--fingerprints"]
    
    OptionsInfo["ClusteringMethod"] = Options["--clusteringMethod"]
    ProcessClusteringMethodParameter()

    OptionsInfo["NumMols"] = int(Options["--numMols"])
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"], Options["--infile"], Options["--outfile"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]

    OptionsInfo["SimilarityMetric"] = Options["--similarityMetric"]
    ProcessSimilarityMetricParameter()
    
    OptionsInfo["ParamsFingerprints"] = Options["--paramsFingerprints"]
    ProcessFingerprintsParameters()
    
def RetrieveOptions():
    """Retrieve command line arguments and options"""
    
    # Get options...
    global Options
    Options = docopt(_docoptUsage_)
    
    # Set current working directory to the specified directory...
    WorkingDir = Options["--workingdir"]
    if WorkingDir:
        os.chdir(WorkingDir)
    
    # Handle examples option...
    if "--examples" in Options and Options["--examples"]:
        MiscUtil.PrintInfo(MiscUtil.GetExamplesTextFromDocOptText(_docoptUsage_))
        sys.exit(0)

def ValidateOptions():
    """Validate option values"""
    
    MiscUtil.ValidateOptionTextValue("-c, --clusteringMethod", Options["--clusteringMethod"], "Centroid CLink Gower McQuitty SLink UPGMA Ward")
    MiscUtil.ValidateOptionTextValue("-f, --fingerprints", Options["--fingerprints"], "AtomPairs MACCS166Keys Morgan MorganFeatures PathLength TopologicalTorsions")
    
    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "MaxMin HierarchicalClustering")
    MiscUtil.ValidateOptionIntegerValue("-n, --numMols", Options["--numMols"], {">": 0})
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd smi txt csv tsv")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd smi")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
        
    MiscUtil.ValidateOptionTextValue("-s, --similarityMetric", Options["--similarityMetric"], "BraunBlanquet Cosine Dice Kulczynski RogotGoldberg Russel Sokal Tanimoto")
    
# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitPickDiverseMolecules.py - Pick a diverse subset of molecules

Usage:
    RDKitPickDiverseMolecules.py [--clusteringMethod <Centroid, CLink...>]
                                 [--fingerprints <MACCS166Keys, Morgan, PathLength...>]
                                 [--infileParams <Name,Value,...>] [--mode <MaxMin or HierarchicalClustering>]
                                 [--numMols <number>]  [--outfileParams <Name,Value,...>] 
                                 [--overwrite] [--paramsFingerprints <Name,Value,...>]
                                 [--similarityMetric <Dice, Tanimoto...>] [-w <dir>] -i <infile> -o <outfile> 
    RDKitPickDiverseMolecules.py -h | --help | -e | --examples

Description:
    Pick a subset of diverse molecules  based on a variety of 2D fingerprints using
    MaxMin [ Ref 135 ] or an available hierarchical clustering methodology and write
    them to a file.

    The default fingerprints types for various fingerprints are shown below:

        AtomPairs              IntSparseIntVect
        MACCS166Keys           ExplicitBitVect
        Morgan                 UIntSparseIntVect
        MorganFeatures         UIntSparseIntVect
        PathLength             ExplicitBitVect
        TopologicalTorsions    LongSparseIntVect
 
    The Dice and Tanimoto similarity functions available in RDKit are able to
    handle fingerprints corresponding to both IntVect and BitVect. All other
    similarity functions, however, expect BitVect fingerprints to calculate
    pairwise similarity. Consequently, ExplicitBitVect fingerprints are generated
    for AtomPairs, Morgan, MorganFeatures, and TopologicalTorsions for
    similarity calculations instead of default IntVect fingerprints.

    The supported input file formats are: SD (.sdf, .sd), SMILES (.smi, .csv, .tsv, .txt)

    The supported output file formats are: SD (.sdf, .sd), SMILES (.smi)

Options:
    -c, --clusteringMethod <Centroid, CLink...>  [default: Centroid]
        Clustering method to use for picking a subset of diverse molecules during
        hierarchical clustering. Supported values: Centroid, CLink, Gower,
        McQuitty, SLink, UPGMA, Ward. This option is ignored for 'MaxMin' value
        of '-m, --mode' option. The Clink and SLink corresponding to CompleteLink
        and SingleLink cluster method.
    -f, --fingerprints <MACCS166Keys, Morgan, PathLength...>  [default: Morgan]
        Fingerprints to use for calculating similarity/distance between molecules.
        Supported values: AtomPairs, MACCS166Keys, Morgan, MorganFeatures, PathLength,
        TopologicalTorsions. The PathLength fingerprints are Daylight like fingerprints.
        The Morgan and MorganFeature fingerprints are circular fingerprints, corresponding
        Scitegic's Extended Connectivity Fingerprints (ECFP) and Features Connectivity
        Fingerprints (FCFP). The values of default parameters for generating fingerprints
        can be modified using '-p, --paramsFingerprints' option.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -i, --infile <infile>
        Input file name.
    --infileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for reading
        molecules from files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD, MOL: removeHydrogens,yes,sanitize,yes,strictParsing,yes
            SMILES: smilesColumn,1,smilesNameColumn,2,smilesDelimiter,space,
                smilesTitleLine,auto,sanitize,yes
            
        Possible values for smilesDelimiter: space, comma or tab.
    -m, --mode <MaxMin or HierarchicalClustering>  [default: MaxMin]
        Pick a diverse subset of molecules using MaxMin or hierarchical clustering
        methodology.
    -n, --numMols <number>  [default: 25]
        Number of diverse molecules to pick.
    -o, --outfile <outfile>
        Output file name.
    --outfileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for writing
        molecules to files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD: compute2DCoords,auto,kekulize,no
            SMILES: kekulize,no,smilesDelimiter,space, smilesIsomeric,yes,
                smilesTitleLine,yes
            
        Default value for compute2DCoords: yes for SMILES input file; no for all other
        file types.
    --overwrite
        Overwrite existing files.
    -p, --paramsFingerprints <Name,Value,...>  [default: auto]
        Parameter values to use for generating fingerprints. The default values
        are dependent on the value of '-f, --fingerprints' option. In general, it is a
        comma delimited list of parameter name and value pairs for the name of
        the fingerprints specified using '-f, --fingerprints' option. The supported
        parameter names along with their default values for valid fingerprints
        names are shown below:
            
            AtomPairs: minLength,1 ,maxLength,30, useChirality,No
            Morgan:   radius,2, useChirality,No
            MorganFeatures:   radius,2, useChirality,No
            PathLength: minPath,1, maxPath,7, fpSize, 2048, bitsPerHash,2
            TopologicalTorsions: useChirality,No
            
    -s, --similarityMetric <Dice, Tanimoto...>  [default: Tanimoto]
        Similarity metric to use for calculating similarity/distance between molecules.
        Possible values: BraunBlanquet, Cosine, Dice, Kulczynski, RogotGoldberg,
        Russel, Sokal, Tanimoto.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To pick 25 diverse molecules using MaxMin methodology, Tanimoto similarity
    metric corresponding to Morgan fingerprints with radius of 2, and write
    out a SMILES file, type:

        % RDKitPickDiverseMolecules.py  -i Sample.smi -o SampleOut.smi

    To pick 50 diverse molecules using MaxMin methodology, Dice similarity metric
    corresponding to PathLength fingerprints with max path length of 6, and write
    out a SD file, type:

        % RDKitPickDiverseMolecules.py  -m MaxMin -f PathLength -s Dice -n 50
          -p 'maxPath,6' -i Sample.sdf -o SampleOut.sdf

    To pick 25 diverse molecules using Centroid hierarchical clustering methodology,
    Tanimoto similarity metric corresponding to Morgan fingerprints with radius of 2,
    and write out a SMILES file, type:

        % RDKitPickDiverseMolecules.py  -m HierarchicalClustering -i Sample.smi
          -o SampleOut.smi

    To pick 50 diverse molecules using Ward hierarchical methodology methodology,
    Dice similarity metric corresponding to MorganFeatures fingerprints with radius
    of 2 along with deploying chirality, and write out a SD file, type:

        % RDKitPickDiverseMolecules.py  -m HierarchicalClustering -c Ward -n 50
          -f MorganFeatures -p 'radius,2,useChirality,No' -i Sample.sdf -o
          SampleOut.sdf

    To pick 25 diverse molecules using MaxMin methodology, Tanimoto similarity
    metric corresponding to Morgan fingerprints with radius of 2 from a CSV SMIKES
    file , SMILES strings in column 1, name in olumn 2, and write out a SD file, type:

        % RDKitPickDiverseMolecules.py  --infileParams
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitClusterMolecules.py, RDKitConvertFileFormat.py, RDKitSearchFunctionalGroups.py,
    RDKitSearchSMARTS.py

Copyright:
    Copyright (C) 2018 Manish Sud. All rights reserved.

    The functionality available in this script is implemented using RDKit, an
    open source toolkit for cheminformatics developed by Greg Landrum.

    This file is part of MayaChemTools.

    MayaChemTools is free software; you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your option) any
    later version.

"""

if __name__ == "__main__":
    main()
