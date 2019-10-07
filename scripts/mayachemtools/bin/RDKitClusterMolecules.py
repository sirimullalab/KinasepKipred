#!/bin/env python
#
# File: RDKitClusterMolecules.py
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
    from rdkit.ML.Cluster import Butina
    from rdkit.SimDivFilters import rdSimDivPickers
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
    ClusterMolecules()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def ClusterMolecules():
    """Cluster molecules."""

    Mols = RetrieveMolecules()
    MolsFingerprints = GenerateFingerprints(Mols)
    MolsClusters = PerformClustering(Mols, MolsFingerprints)
    
    WriteMolecules(MolsClusters)

def PerformClustering(Mols, MolsFingerprints):
    """Perform clustering."""

    ClusteredMols = []
    if re.match("^Butina$", OptionsInfo["ClusteringMethod"], re.I):
        return PerformButinaClustering(Mols, MolsFingerprints)
    else:
        return PerformHierarchicalClustering(Mols, MolsFingerprints)
    
    return ClusteredMols

def PerformButinaClustering(Mols, MolsFingerprints):
    """Perform clustering using Butina methodology."""

    MiscUtil.PrintInfo("\nClustering molecules using Butina methodology and %s similarity metric..." % OptionsInfo["SimilarityMetric"])
    
    FingerprintsCount = len(MolsFingerprints)
    DistanceCutoff =  1 - OptionsInfo["ButinaSimilarityCutoff"]
    Reordering = OptionsInfo["ButinaReordering"]
    
    DistanceMatrix = GenerateLowerTriangularDistanceMatrix(MolsFingerprints)

    ClusteredMolIndices = Butina.ClusterData(DistanceMatrix, FingerprintsCount, DistanceCutoff, reordering = Reordering, isDistData = True)

    MolsClusters = []
    for Cluster in ClusteredMolIndices:
        MolsCluster = [Mols[MolIndex] for MolIndex in Cluster]
        MolsClusters.append(MolsCluster)

    return MolsClusters

def PerformHierarchicalClustering(Mols, MolsFingerprints):
    """Perform hierarchical clustering."""

    try:
        import numpy
    except ImportError:
        MiscUtil.PrintError("Failed to import numpy python module. This is required to cluster molecules using hierarchical clustering methodology.")
    
    if OptionsInfo["NumClusters"] > len(Mols):
        MiscUtil.PrintError("The number of clusters, %d, specified using \"-n, --numClusters\" must be less than total number of valid molecules, %d" % (OptionsInfo["NumClusters"], len(Mols)))
    
    MiscUtil.PrintInfo("\nCluster molecules using %s hierarchical clustering methodology and %s similarity metric......" % (OptionsInfo["SpecifiedHierarchicalClusteringMethod"], OptionsInfo["SimilarityMetric"]))
    
    NumFingerprints = len(MolsFingerprints)
    NumClusters = OptionsInfo["NumClusters"]
    DistanceMatrix = GenerateLowerTriangularDistanceMatrix(MolsFingerprints)
    
    ClusterPicker = HierarchicalClusterPicker(OptionsInfo["SpecifiedHierarchicalClusteringMethodID"])
    ClusteredMolIndices = ClusterPicker.Cluster(numpy.asarray(DistanceMatrix), NumFingerprints, NumClusters)

    MolsClusters = []
    for Cluster in ClusteredMolIndices:
        MolsCluster = [Mols[MolIndex] for MolIndex in Cluster]
        MolsClusters.append(MolsCluster)
    
    return MolsClusters

def WriteMolecules(MolsClusters):
    """Write out molecules for each cluster along with cluster numbers."""

    ClustersCount = len(MolsClusters)
    
    SingleOutFileMode = OptionsInfo["SingleOutFileMode"]
    TextOutFileMode = OptionsInfo["TextOutFileMode"]
    TextOutFileDelim = OptionsInfo["TextOutFileDelim"]

    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    
    SMILESIsomeric = OptionsInfo["OutfileParams"]["SMILESIsomeric"]
    SMILESKekulize = OptionsInfo["OutfileParams"]["Kekulize"]
    
    # Setup outfile names and writers...
    SetupClustersOutFilesNames(len(MolsClusters))
    SingleClusterWriter, ClustersOutfilesWriters = SetupMoleculeWriters(ClustersCount)

    MolCount = 0
    SingleMolClustersCount = 0
    
    if SingleOutFileMode:
        Writer = SingleClusterWriter
    
    for ClusterIndex in range(0, ClustersCount):
        MolsCluster = MolsClusters[ClusterIndex]
        ClusterNum = ClusterIndex + 1

        if len(MolsCluster) == 1:
            SingleMolClustersCount += 1
        
        if not SingleOutFileMode:
            Writer = ClustersOutfilesWriters[ClusterIndex]
            
        for Mol in MolsCluster:
            MolCount += 1

            if TextOutFileMode:
                # Write out text file including SMILES file...
                SMILES = Chem.MolToSmiles(Mol, isomericSmiles = SMILESIsomeric, kekuleSmiles = SMILESKekulize)
                MolName = RDKitUtil.GetMolName(Mol, MolCount)
                Line = TextOutFileDelim.join([SMILES, MolName, "%d" % ClusterNum])
                Writer.write("%s\n" % Line)
            else:
                # Write out SD file...
                Mol.SetProp("ClusterNumber", "%s" % ClusterNum)
                if Compute2DCoords:
                    AllChem.Compute2DCoords(Mol)
                Writer.write(Mol)
    
    if SingleClusterWriter is not None:
        SingleClusterWriter.close()
    for ClusterOutfileWriter in ClustersOutfilesWriters:
        ClusterOutfileWriter.close()

    MiscUtil.PrintInfo("\nTotal number of clusters: %d" % ClustersCount)

    if ClustersCount > 0:
        MiscUtil.PrintInfo("\nNumber of clusters containing only a single molecule: %d" % SingleMolClustersCount)
        MiscUtil.PrintInfo("Average number of molecules per cluster: %.1f" % (MolCount/ClustersCount))
    
        MiscUtil.PrintInfo("\nNumber of molecules in each cluster:")
        MiscUtil.PrintInfo("ClusterNumber,MolCount")
        ClusterNum = 0
        for MolsCluster in MolsClusters:
            ClusterNum += 1
            MiscUtil.PrintInfo("%d,%d" % (ClusterNum, len(MolsCluster)))

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

def SetupMoleculeWriters(ClustersCount):
    """Set up molecule writers for SD and text files."""
    
    Writer = None
    ClustersOutfilesWriters = []

    TextOutFileMode = OptionsInfo["TextOutFileMode"]
    TextOutFileDelim = OptionsInfo["TextOutFileDelim"]
    TextOutFileTitleLine = OptionsInfo["TextOutFileTitleLine"]
    
    if OptionsInfo["SingleOutFileMode"]:
        Outfile = OptionsInfo["Outfile"]
        if TextOutFileMode:
            Writer = open(Outfile, "w")
        else:
            Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
        if Writer is None:
            MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
        
        if TextOutFileMode:
            if TextOutFileTitleLine:
                WriteTextFileHeaderLine(Writer, TextOutFileDelim)
        
        MiscUtil.PrintInfo("Generating file %s..." % Outfile)
    else:
        for ClusterIndex in range(0, ClustersCount):
            Outfile = OptionsInfo["ClustersOutfiles"][ClusterIndex]
            if TextOutFileMode:
                ClusterWriter = open(Outfile, "w")
            else:
                ClusterWriter = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
            if ClusterWriter is None:
                MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
            
            if TextOutFileMode:
                if TextOutFileTitleLine:
                    WriteTextFileHeaderLine(ClusterWriter, TextOutFileDelim)
        
            ClustersOutfilesWriters.append(ClusterWriter)
        
        if ClustersCount > 4:
            MiscUtil.PrintInfo("Generating %d output files with the following file name format: %s_Cluster<Num>.%s" % (ClustersCount, OptionsInfo["OutfileBasename"], OptionsInfo["OutfileExt"]))
        else:
            Delmiter = ','
            OutfileNames = Delmiter.join(OptionsInfo["ClustersOutfiles"])
            MiscUtil.PrintInfo("Generating %d output files: %s..." % (ClustersCount, OutfileNames))
        
    return (Writer, ClustersOutfilesWriters)

def WriteTextFileHeaderLine(Writer, TextOutFileDelim):
    """Write out a header line for text files including SMILEs file."""
    
    Line = TextOutFileDelim.join(["SMILES", "Name", "ClusterNumber"])
    Writer.write("%s\n" % Line)

def SetupClustersOutFilesNames(ClustersCount):
    """Set up out file names for clusters."""

    OptionsInfo["ClustersOutfiles"] = []
    if OptionsInfo["SingleOutFileMode"] or ClustersCount == 0:
        # Nothing to do...
        return

    OutfileBasename = OptionsInfo["OutfileBasename"]
    OutfileExt = OptionsInfo["OutfileExt"]
    
    ClusterOutfiles = []
    for ClusterIndex in range(0, ClustersCount):
        ClusterNum = ClusterIndex + 1
        ClusterOutfile = "%s_Cluster%d.%s" % (OutfileBasename, ClusterNum, OutfileExt)
        ClusterOutfiles.append(ClusterOutfile)
    
    OptionsInfo["ClustersOutfiles"] = ClusterOutfiles
    
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

    # Setup canonical parameter names for specified fingerprints...
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

    OptionsInfo["SpecifiedHierarchicalClusteringMethod"] = ""
    OptionsInfo["SpecifiedHierarchicalClusteringMethodID"] = ""
    
    if re.match("^Butina$", OptionsInfo["ClusteringMethod"], re.I):
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
        MiscUtil.PrintError("The clustering method, %s, specified using \"-c, --clusteringMethod\" option is not a valid name." % (OptionsInfo["ClusteringMethod"]))

    SpecifiedHierarchicalClusteringMethodName = CanonicalClusteringMethodNameMap[CanonicalName]
    OptionsInfo["SpecifiedHierarchicalClusteringMethod"] = SpecifiedHierarchicalClusteringMethodName
    OptionsInfo["SpecifiedHierarchicalClusteringMethodID"] = ClusteringMethodInfoMap[SpecifiedHierarchicalClusteringMethodName] 
    
def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["ButinaSimilarityCutoff"] = float(Options["--butinaSimilarityCutoff"])
    OptionsInfo["ButinaReordering"] = False
    if re.match("^Yes$", Options["--butinaReordering"], re.I):
        OptionsInfo["ButinaReordering"] = True
    
    OptionsInfo["Fingerprints"] = Options["--fingerprints"]
    
    OptionsInfo["ClusteringMethod"] = Options["--clusteringMethod"]
    ProcessClusteringMethodParameter()

    OptionsInfo["NumClusters"] = int(Options["--numClusters"])
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"], Options["--infile"], Options["--outfile"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]

    OptionsInfo["OutFileMode"] = Options["--outfileMode"]
    SingleOutFileMode = True
    if not re.match("^SingleFile$", Options["--outfileMode"], re.I):
        SingleOutFileMode = False
    OptionsInfo["SingleOutFileMode"] = SingleOutFileMode
    
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(Options["--outfile"])
    OptionsInfo["OutfileBasename"] = FileName
    OptionsInfo["OutfileExt"] = FileExt

    TextOutFileMode = False
    TextOutFileDelim = ""
    TextOutFileTitleLine = True
    
    if MiscUtil.CheckFileExt(Options["--outfile"], "csv"):
        TextOutFileMode = True
        TextOutFileDelim = ","
    elif MiscUtil.CheckFileExt(Options["--outfile"], "tsv txt"):
        TextOutFileMode = True
        TextOutFileDelim = "\t"
    elif MiscUtil.CheckFileExt(Options["--outfile"], "smi"):
        TextOutFileMode = True
        TextOutFileDelim = OptionsInfo["OutfileParams"]["SMILESDelimiter"]
        TextOutFileTitleLine = OptionsInfo["OutfileParams"]["SMILESTitleLine"]
        
    OptionsInfo["TextOutFileMode"] = TextOutFileMode
    OptionsInfo["TextOutFileDelim"] = TextOutFileDelim
    OptionsInfo["TextOutFileTitleLine"] = TextOutFileTitleLine
    
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
    
    MiscUtil.ValidateOptionFloatValue("-b, --butinaSimilarityCutoff", Options["--butinaSimilarityCutoff"], {">": 0.0, "<=" : 1.0})
    MiscUtil.ValidateOptionTextValue("--butinaReordering", Options["--butinaReordering"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("-c, --clusteringMethod", Options["--clusteringMethod"], "Butina Centroid CLink Gower McQuitty SLink UPGMA Ward")
    MiscUtil.ValidateOptionTextValue("-f, --fingerprints", Options["--fingerprints"], "AtomPairs MACCS166Keys Morgan MorganFeatures PathLength TopologicalTorsions")
    
    MiscUtil.ValidateOptionIntegerValue("-n, --numClusters", Options["--numClusters"], {">": 0})
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi csv tsv txt")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd smi csv tsv txt")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
        
    MiscUtil.ValidateOptionTextValue("--outfileMode", Options["--outfileMode"], "SingleFile MultipleFiles")
    
    MiscUtil.ValidateOptionTextValue("-s, --similarityMetric", Options["--similarityMetric"], "BraunBlanquet Cosine Dice Kulczynski RogotGoldberg Russel Sokal Tanimoto")
    
# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitClusterMolecules.py - Cluster molecules using 2D fingerprints

Usage:
    RDKitClusterMolecules.py [--butinaSimilarityCutoff <number>]  [--butinaReordering <yes or no>]
                             [--clusteringMethod <Butina, Centroid, CLink...>]
                             [--fingerprints <MACCS166Keys, Morgan, PathLength...> ] [--infileParams <Name,Value,...>]
                             [--numClusters <number>] [--outfileMode <SingleFile or MultipleFiles>]
                             [ --outfileParams <Name,Value,...> ] [--overwrite] [--paramsFingerprints <Name,Value,...>]
                             [--similarityMetric <Dice, Tanimoto...>] [-w <dir>] -i <infile> -o <outfile> 
    RDKitClusterMolecules.py -h | --help | -e | --examples

Description:
    Cluster molecules based on a variety of 2D fingerprints using Butina [ Ref 136 ] or any
    other available hierarchical clustering methodology and write them to output file(s).

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
    pairwise similarity. Consequently, ExplicitBitVect fingerprints, instead of
    default IntVect fingerprints, are generated for AtomPairs, Morgan,
    MorganFeatures, and TopologicalTorsions to calculate similarity.

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd), SMILES (.smi,
    .txt, .csv, .tsv)

    The supported output file formats are: SD (.sdf, .sd), SMILES (.smi), CSV/TSV
    (.csv, .tsv, .txt)

Options:
    -b, --butinaSimilarityCutoff <number>  [default: 0.55]
        Similarity cutoff to use during Butina clustering. The molecule pairs with
        similarity value greater than specified value or distance less than '1 - specified
        value' are considered neighbors. This value is only used during 'Butina' value
        of '-c, --clusteringMethod' option and determines the number of clusters
        during the clustering of molecules. It is ignored for all other clustering methods.
    --butinaReordering <yes or no>  [default: no]
        Update number of neighbors for unassigned molecules after creating a new
        cluster in order to insure that the molecule with the largest number of
        unassigned neighbors is selected as the next cluster center.
    -c, --clusteringMethod <Butina, Centroid, CLink...>  [default: Butina]
        Clustering method to use for clustering molecules. Supported values:
        Butina, Centroid, CLink, Gower, McQuitty, SLink, UPGMA, Ward.
        Butina is an unsupervised database clustering method to automatically
        cluster small and large data sets. All other clustering methods correspond
        to hierarchical clustering and require a priori specification of number of
        clusters to be generated.
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
    -n, --numClusters <number>  [default: 10]
        Number of clusters to generate during hierarchical clustering. This option is
        ignored for 'Butina' value of '-c, --clusteringMethod' option.
    -o, --outfile <outfile>
        Output file name.
    --outfileMode <SingleFile or MultipleFiles>  [default: SingleFile]
        Write out a single file containing molecule clusters or generate an individual file
        for each cluster. Possible values: SingleFile or MultipleFiles. The molecules are
        grouped for each cluster before they are written to output file(s) along with
        appropriate cluster numbers. The cluster number is also appended to output
        file names during generation of multiple output files.
    --outfileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for writing
        molecules to files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD: compute2DCoords,auto,kekulize,no
            SMILES: kekulize,no,smilesDelimiter,space, smilesIsomeric,yes,
                smilesTitleLine,yes
            
        Default value for compute2DCoords: yes for SMILES input file; no for all other
        file types. The kekulize and smilesIsomeric parameters are also used during
        generation of SMILES strings for CSV/TSV files.
    --overwrite
        Overwrite existing files.
    -p, --paramsFingerprints <Name,Value,...>  [default: auto]
        Parameter values to use for generating fingerprints. The default values
        are dependent on the value of '-f, --fingerprints' option. In general, it is a
        comma delimited list of parameter name and value pairs for the name of
        fingerprints specified using '-f, --fingerprints' option. The supported
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
    To cluster molecules using Butina methodology at a similarity cutoff of 0.55
    with automatic determination of number of clusters, Tanimoto similarity
    metric corresponding to Morgan fingerprints with radius of 2, and write out
    a single SMILES file containing clustered molecules along with cluster number
    for each molecule, type:

        % RDKitClusterMolecules.py  -i Sample.smi -o SampleOut.smi

    To cluster molecules using Butina methodology at similarity cutoff of 0.45
    with automatic determination of number of clusters, Dice similarity metric
    corresponding to Morgan fingerprints with radius of 2, and write out multiple
    SD files containing clustered molecules for each cluster, type:

        % RDKitClusterMolecules.py  -b 0.45 --outfileMode MultipleFiles
          -i Sample.smi -o SampleOut.sdf

    To cluster molecules using Ward hierarchical methodology to generate 15
    clusters, Dice similarity metric corresponding to Pathlength fingerprints with 
    path length between 1 and 7,  and write out a single TSV file for clustered
    molecules along with cluster numner for each molecule, type:

        % RDKitClusterMolecules.py  -c Ward -f PathLength -n 15
          -p 'minPath,1, maxPath,7' -i Sample.sdf -o SampleOut.tsv

    To cluster molecules using Centroid hierarchical methodology to generate 5
    clusters, Dice similarity metric corresponding to MACCS166Keys fingerprints
    for molecules in a SMILES CSV file, SMILES strings in column 1, name in
    column 2, and write out a single SD file for clustered molecules along with
    cluster numner for each molecule, type:

        % RDKitClusterMolecules.py  -c Centroid -f MACCS166Keys --infileParams
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitPickDiverseMolecules.py, RDKitSearchFunctionalGroups.py,
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
