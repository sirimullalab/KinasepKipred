#!/bin/env python
#
# File: RDKitPerformRGroupDecomposition.py
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
    from rdkit.Chem import rdRGroupDecomposition as rgd
    from rdkit.Chem import rdFMCS
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
    PerformRGroupDecomposition()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformRGroupDecomposition():
    """Perform R group decomposition."""

    # Retrieve molecules...
    Mols = RetrieveMolecules()

    # Identify R groups and write them out...
    RGroups, UnmatchedMolIndices = RetrieveRGroups(Mols)
    WriteRGroups(Mols, RGroups, UnmatchedMolIndices)
    
def RetrieveRGroups(Mols):
    """Retrieve R groups"""
    
    CoreMols = SetupCoreScaffolds(Mols)
    DecompositionParams = SetupRGroupDecompositionParams()
    RGroupDecompositionObject = rgd.RGroupDecomposition(CoreMols, DecompositionParams)
    
    MiscUtil.PrintInfo("\nPerforming R group decomposition...")
    
    UnmatchedMolIndices = []
    for MolIndex, Mol in enumerate(Mols):
        Status = RGroupDecompositionObject.Add(Mol)
        if Status < 0:
            UnmatchedMolIndices.append(MolIndex)
    
    if not RGroupDecompositionObject.Process():
        MiscUtil.PrintWarning("R group decomposition failed to match any molecule to core scaffold(s)...")
    
    RGroups = RGroupDecompositionObject.GetRGroupsAsColumns(asSmiles=True)

    return (RGroups, UnmatchedMolIndices)
    
def SetupCoreScaffolds(Mols):
    """Setup core scaffold molecules(s)"""

    if re.match("^(BySMARTS|BySMILES)$", OptionsInfo["CoreScaffold"], re.I):
        return SetupCoreScaffoldsBySMARTSOrSMILES()
    elif re.match("^ByMCS$", OptionsInfo["CoreScaffold"], re.I):
        return SetupCoreScaffoldsByMCS(Mols)
    else:
        MiscUtil.PrintError("The  value, %s, specified for  \"-c, --coreScaffold\" option is not supported." % (OptionsInfo["CoreScaffold"]))
        
def SetupCoreScaffoldsBySMARTSOrSMILES():
    """Setup core scaffold molecules(s) using specified SMARTS or SMILES."""

    BySMARTS = True if re.match("^BySMARTS$", OptionsInfo["CoreScaffold"], re.I) else False
    CoreScaffoldList = OptionsInfo["SMARTSOrSMILESCoreScaffoldList"]

    if BySMARTS:
        MiscUtil.PrintInfo("\nSetting up core scaffold(s) using SMARTS...\nSMARTS core scaffold(s): %s" % " ".join(CoreScaffoldList))
    else:
        MiscUtil.PrintInfo("\nSetting up core scaffold(s) using SMILES...\nSMILES core scaffold(s): %s" % " ".join(CoreScaffoldList))
    
    CoreMols = []
    for Core in CoreScaffoldList:
        if BySMARTS:
            CoreMol = Chem.MolFromSmarts(Core)
        else:
            CoreMol = Chem.MolFromSmiles(Core)
        if CoreMol is None:
            MiscUtil.PrintError("Failed to generate mol for core scaffold: %s" % (Core))
        CoreMols.append(CoreMol)

    return CoreMols

def SetupCoreScaffoldsByMCS(Mols):
    """Setup core scaffold molecule using MCS."""

    MiscUtil.PrintInfo("\nSetting up core scaffold using MCS...")

    MCSParams = OptionsInfo["MCSParams"]
    
    CoreMols = []

    MCSResultObject = rdFMCS.FindMCS(Mols, maximizeBonds = MCSParams["MaximizeBonds"], threshold = MCSParams["Threshold"], timeout = MCSParams["TimeOut"], verbose = MCSParams["Verbose"], matchValences = MCSParams["MatchValences"], ringMatchesRingOnly = MCSParams["RingMatchesRingOnly"], completeRingsOnly = MCSParams["CompleteRingsOnly"], matchChiralTag = MCSParams["MatchChiralTag"], atomCompare = MCSParams["AtomCompare"], bondCompare = MCSParams["BondCompare"], seedSmarts = MCSParams["SeedSMARTS"]) 

    if MCSResultObject.canceled:
        MiscUtil.PrintError("MCS failed to identify a core scaffold. Specify a different set of parameters using \"-m, --mcsParams\" option and try again.")
    
    CoreNumAtoms = MCSResultObject.numAtoms
    CoreNumBonds = MCSResultObject.numBonds
    SMARTSCore = MCSResultObject.smartsString
    
    if not len(SMARTSCore):
        MiscUtil.PrintError("MCS failed to identify a core scaffold. Specify a different set of parameters using \"-m, --mcsParams\" option and try again.")
        
    MiscUtil.PrintInfo("SMARTS core scaffold: %s\nNumber of atoms in core scaffold: %s\nNumber of bonds in core scaffold: %s" % (SMARTSCore, CoreNumAtoms, CoreNumBonds))

    if CoreNumAtoms < MCSParams["MinNumAtoms"]:
        MiscUtil.PrintError("Number of atoms, %d, in core scaffold identified by MCS is less than, %d, as specified by \"minNumAtoms\" parameter in  \"-m, --mcsParams\" option." % (CoreNumAtoms, MCSParams["MinNumAtoms"]))
    
    if CoreNumBonds < MCSParams["MinNumBonds"]:
        MiscUtil.PrintError("Number of bonds, %d, in core scaffold identified by MCS is less than, %d, as specified by \"minNumBonds\" parameter in  \"-m, --mcsParams\" option." % (CoreNumBonds, MCSParams["MinNumBonds"]))
    
    CoreMol = Chem.MolFromSmarts(SMARTSCore)
    CoreMols.append(CoreMol)

    return CoreMols

def SetupRGroupDecompositionParams():
    """Setup R group decomposition parameters"""

    DecompositionParams = rgd.RGroupDecompositionParameters()

    DecompositionParams.alignment = OptionsInfo["DecompositionParams"]["RGroupCoreAlignment"]
    DecompositionParams.chunkSize = OptionsInfo["DecompositionParams"]["chunkSize"]
    DecompositionParams.matchingStrategy = OptionsInfo["DecompositionParams"]["RGroupMatching"]
    DecompositionParams.onlyMatchAtRGroups = OptionsInfo["DecompositionParams"]["matchOnlyAtRGroups"]
    DecompositionParams.removeAllHydrogenRGroups = OptionsInfo["DecompositionParams"]["removeHydrogenOnlyGroups"]
    DecompositionParams.removeHydrogensPostMatch = OptionsInfo["DecompositionParams"]["removeHydrogensPostMatch"]

    return DecompositionParams
    
def WriteRGroups(Mols, RGroups, UnmatchedMolIndices):
    """Write out R groups"""

    Outfile = OptionsInfo["Outfile"]
    UnmatchedOutfile = OptionsInfo["UnmatchedOutfile"]
    RemoveUnmatchedMode = OptionsInfo["RemoveUnmatchedMode"]
    
    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    
    TextOutFileMode = OptionsInfo["TextOutFileMode"]
    TextOutFileDelim = OptionsInfo["TextOutFileDelim"]
    Quote = OptionsInfo["TextOutQuote"]

    SMILESIsomeric = OptionsInfo["OutfileParams"]["SMILESIsomeric"]
    SMILESKekulize = OptionsInfo["OutfileParams"]["Kekulize"]
    
    # Setup writers...
    Writer = None
    UnmatchedWriter = None
    if TextOutFileMode:
        Writer = open(Outfile, "w")
        if RemoveUnmatchedMode:
            UnmatchedWriter = open(UnmatchedOutfile, "w")
    else:
        Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
        if RemoveUnmatchedMode:
            UnmatchedWriter = RDKitUtil.MoleculesWriter(UnmatchedOutfile, **OptionsInfo["OutfileParams"])
        
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
    if RemoveUnmatchedMode:
        if UnmatchedWriter is None:
            MiscUtil.PrintError("Failed to setup a writer for output fie %s " % UnmatchedOutfile)
    
    if RemoveUnmatchedMode:
        MiscUtil.PrintInfo("\nGenerating files: %s %s..." % (Outfile, UnmatchedOutfile))
    else:
        MiscUtil.PrintInfo("\nGenerating file %s..." % Outfile)
        
    # Set up data keys and labels for  core and R groups...
    RGroupsDataKeys = []
    RGroupsDataLabels = []
    CoreDataLabelPresent = False
    for DataLabel in sorted(RGroups):
        if re.match("^Core", DataLabel, re.I):
            CoreDataLabelPresent = True
            RGroupsDataKeys.append(DataLabel)
            RGroupsDataLabels.append("SMILES%s" % DataLabel)
        elif re.match("^R", DataLabel, re.I):
            RGroupsDataKeys.append(DataLabel)
            RGroupsDataLabels.append("SMILES%s" % DataLabel)
        else:
            MiscUtil.PrintWarning("Ignoring unknown R group data label, %s, found during R group decomposition..." % DataLabel)

    if CoreDataLabelPresent:
        RGroupsCategoriesCount = len(RGroupsDataLabels) - 1
    else:
        RGroupsCategoriesCount = len(RGroupsDataLabels)
    RGroupsMolUnmatchedCount = len(UnmatchedMolIndices)
    RGroupsMolMatchedCount = len(Mols) - RGroupsMolUnmatchedCount
    
    # Wite out headers for a text file...
    if TextOutFileMode:
        LineWords = ["SMILES","Name"]
        
        if RemoveUnmatchedMode:
            Line = MiscUtil.JoinWords(LineWords, TextOutFileDelim, Quote)
            UnmatchedWriter.write("%s\n" % Line)
        
        LineWords.extend(RGroupsDataLabels)
        Line = MiscUtil.JoinWords(LineWords, TextOutFileDelim, Quote)
        Writer.write("%s\n" % Line)

    MolCount = 0
    RGroupsResultIndex = -1
    
    for MolIndex, Mol in enumerate(Mols):
        MolCount += 1

        UnmatchedMol = False
        if MolIndex in UnmatchedMolIndices:
            UnmatchedMol = True
        
        if UnmatchedMol:
            RGroupsDataSMILES = [""] * len(RGroupsDataKeys)
        else:
            RGroupsResultIndex += 1
            RGroupsDataSMILES = [RGroups[RGroupsDataKey][RGroupsResultIndex] for RGroupsDataKey in RGroupsDataKeys]

        if TextOutFileMode:
            # Write out text file including SMILES file...
            MolSMILES = Chem.MolToSmiles(Mol, isomericSmiles = SMILESIsomeric, kekuleSmiles = SMILESKekulize)
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            LineWords = [MolSMILES, MolName]

            if UnmatchedMol and RemoveUnmatchedMode:
                Line = MiscUtil.JoinWords(LineWords, TextOutFileDelim, Quote)
                UnmatchedWriter.write("%s\n" % Line)
            else:
                LineWords.extend(RGroupsDataSMILES)
                Line = MiscUtil.JoinWords(LineWords, TextOutFileDelim, Quote)
                Writer.write("%s\n" % Line)
        else:
            # Write out SD file...
            if Compute2DCoords:
                AllChem.Compute2DCoords(Mol)
            
            if UnmatchedMol and RemoveUnmatchedMode:
                UnmatchedWriter.write(Mol)
            else:
                for (Name, Value) in zip(RGroupsDataLabels, RGroupsDataSMILES):
                    Mol.SetProp(Name, Value)
                Writer.write(Mol)
    
    if Writer is not None:
        Writer.close()
    if UnmatchedWriter is not None:
        UnmatchedWriter.close()

    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of  R group categories: %d" % RGroupsCategoriesCount)
    MiscUtil.PrintInfo("Number of  matched molecules containing core scaffold(s): %d" % RGroupsMolMatchedCount)
    MiscUtil.PrintInfo("Number of  unmatched molecules containing no core scaffold(s): %d" % RGroupsMolUnmatchedCount)
    
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

def ProcessDecompositionParameters():
    """Set up and process decomposition parameters."""

    SetupDecompositionParameters()
    ProcessSpecifiedDecompositionParameters()

def SetupDecompositionParameters():
    """Set up default decomposition parameters."""
    
    OptionsInfo["DecompositionParams"] = {"RGroupCoreAlignment": rgd.RGroupCoreAlignment.MCS, "RGroupMatching": rgd.RGroupMatching.GreedyChunks, "chunkSize": 5, "matchOnlyAtRGroups": False, "removeHydrogenOnlyGroups": True, "removeHydrogensPostMatch": False}
    
def ProcessSpecifiedDecompositionParameters():
    """Process specified decomposition parameters."""

    if re.match("^auto$", OptionsInfo["SpecifiedDecompositionParams"], re.I):
        # Nothing to process...
        return

    # Parse specified parameters...
    DecompositionParams = re.sub(" ", "", OptionsInfo["SpecifiedDecompositionParams"])
    if not DecompositionParams:
        MiscUtil.PrintError("No valid parameter name and value pairs specified using \"-d, --decompositionParams\" option.")

    DecompositionParamsWords = DecompositionParams.split(",")
    if len(DecompositionParamsWords) % 2:
        MiscUtil.PrintError("The number of comma delimited paramater names and values, %d, specified using \"-d, --decompositionParams\" option must be an even number." % (len(DecompositionParamsWords)))
    
    # Setup  canonical parameter names...
    ValidParamNames = []
    CanonicalParamNamesMap = {}
    for ParamName in sorted(OptionsInfo["DecompositionParams"]):
        ValidParamNames.append(ParamName)
        CanonicalParamNamesMap[ParamName.lower()] = ParamName

    # Validate and set paramater names and value...
    for Index in range(0, len(DecompositionParamsWords), 2):
        Name = DecompositionParamsWords[Index]
        Value = DecompositionParamsWords[Index + 1]

        CanonicalName = Name.lower()
        if  not CanonicalName in CanonicalParamNamesMap:
            MiscUtil.PrintError("The parameter name, %s, specified using \"-d, --decompositionParams\" option is not a valid name. Supported parameter names: %s" % (Name,  " ".join(ValidParamNames)))

        ParamName = CanonicalParamNamesMap[CanonicalName]
        if re.match("^RGroupCoreAlignment$", ParamName, re.I):
            if re.match("^MCS$", Value, re.I):
                ParamValue = rgd.RGroupCoreAlignment.MCS
            elif re.match("^None$", Value, re.I):
                ParamValue = rgd.RGroupCoreAlignment.names['None']
            else:
                MiscUtil.PrintError("The parameter value, %s, specified using \"-d, --decompositionParams\" option  for parameter, %s, is not a valid value. Supported values: MCS None" % (Value, Name))
        elif re.match("^RGroupMatching$", ParamName, re.I):
            if re.match("^Greedy$", Value, re.I):
                ParamValue = rgd.RGroupMatching.Greedy
            elif re.match("^GreedyChunks$", Value, re.I):
                ParamValue = rgd.RGroupMatching.GreedyChunks
            elif re.match("^Exhaustive$", Value, re.I):
                ParamValue = rgd.RGroupMatching.Exhaustive
            else:
                MiscUtil.PrintError("The parameter value, %s, specified using \"-d, --decompositionParams\" option  for parameter, %s, is not a valid value. Supported values: Greedy GreedyChunks Exhaustive" % (Value, Name))
        elif re.match("^chunkSize$", ParamName, re.I):
            Value = int(Value)
            if Value <= 0 :
                MiscUtil.PrintError("The parameter value, %s, specified using \"-d, --decompositionParams\" option  for parameter, %s, is not a valid value. Supported values: > 0" % (Value, Name))
            ParamValue = Value
        else:
            if not re.match("^(Yes|No|True|False)$", Value, re.I):
                MiscUtil.PrintError("The parameter value, %s, specified using \"-d, --decompositionParams\" option  for parameter, %s, is not a valid value. Supported values: Yes No True False" % (Value, Name))
            ParamValue = False
            if re.match("^(Yes|True)$", Value, re.I):
                ParamValue = True
        
        # Set value...
        OptionsInfo["DecompositionParams"][ParamName] = ParamValue

def ProcessMCSParameters():
    """Set up and process MCS parameters."""

    SetupMCSParameters()
    ProcessSpecifiedMCSParameters()

def SetupMCSParameters():
    """Set up default MCS parameters."""
    
    OptionsInfo["MCSParams"] = {"MaximizeBonds": True, "Threshold": 0.9, "TimeOut": 3600, "Verbose": False, "MatchValences": True, "MatchChiralTag": False, "RingMatchesRingOnly": True, "CompleteRingsOnly": True, "AtomCompare": rdFMCS.AtomCompare.CompareElements, "BondCompare": rdFMCS.BondCompare.CompareOrder, "SeedSMARTS": "", "MinNumAtoms": 1, "MinNumBonds": 0}
    
def ProcessSpecifiedMCSParameters():
    """Process specified MCS parameters."""

    if re.match("^auto$", OptionsInfo["SpecifiedMCSParams"], re.I):
        # Nothing to process...
        return
    
    # Parse specified parameters...
    MCSParams = re.sub(" ", "", OptionsInfo["SpecifiedMCSParams"])
    if not MCSParams:
        MiscUtil.PrintError("No valid parameter name and value pairs specified using \"-m, --mcsParams\" option.")

    MCSParamsWords = MCSParams.split(",")
    if len(MCSParamsWords) % 2:
        MiscUtil.PrintError("The number of comma delimited paramater names and values, %d, specified using \"-m, --mcsParams\" option must be an even number." % (len(MCSParamsWords)))
    
    # Setup  canonical parameter names...
    ValidParamNames = []
    CanonicalParamNamesMap = {}
    for ParamName in sorted(OptionsInfo["MCSParams"]):
        ValidParamNames.append(ParamName)
        CanonicalParamNamesMap[ParamName.lower()] = ParamName

    # Validate and set paramater names and value...
    for Index in range(0, len(MCSParamsWords), 2):
        Name = MCSParamsWords[Index]
        Value = MCSParamsWords[Index + 1]

        CanonicalName = Name.lower()
        if  not CanonicalName in CanonicalParamNamesMap:
            MiscUtil.PrintError("The parameter name, %s, specified using \"-m, --mcsParams\" option is not a valid name. Supported parameter names: %s" % (Name,  " ".join(ValidParamNames)))

        ParamName = CanonicalParamNamesMap[CanonicalName]
        if re.match("^Threshold$", ParamName, re.I):
            Value = float(Value)
            if Value <= 0.0 or Value > 1.0 :
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: > 0 and <= 1.0" % (Value, Name))
            ParamValue = Value
        elif re.match("^Timeout$", ParamName, re.I):
            Value = int(Value)
            if Value <= 0:
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: > 0" % (Value, Name))
            ParamValue = Value
        elif re.match("^MinNumAtoms$", ParamName, re.I):
            Value = int(Value)
            if Value < 1:
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: >= 1" % (Value, Name))
            ParamValue = Value
        elif re.match("^MinNumBonds$", ParamName, re.I):
            Value = int(Value)
            if Value < 0:
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: >=0 " % (Value, Name))
            ParamValue = Value
        elif re.match("^AtomCompare$", ParamName, re.I):
            if re.match("^CompareAny$", Value, re.I):
                ParamValue = rdFMCS.AtomCompare.CompareAny
            elif re.match("^CompareElements$", Value, re.I):
                ParamValue = Chem.rdFMCS.AtomCompare.CompareElements
            elif re.match("^CompareIsotopes$", Value, re.I):
                ParamValue = Chem.rdFMCS.AtomCompare.CompareIsotopes
            else:
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: CompareAny CompareElements CompareIsotopes" % (Value, Name))
        elif re.match("^BondCompare$", ParamName, re.I):
            if re.match("^CompareAny$", Value, re.I):
                ParamValue = Chem.rdFMCS.BondCompare.CompareAny
            elif re.match("^CompareOrder$", Value, re.I):
                ParamValue = rdFMCS.BondCompare.CompareOrder
            elif re.match("^CompareOrderExact$", Value, re.I):
                ParamValue = rdFMCS.BondCompare.CompareOrderExact
            else:
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: CompareAny CompareOrder CompareOrderExact" % (Value, Name))
        elif re.match("^SeedSMARTS$", ParamName, re.I):
            if not len(Value):
                MiscUtil.PrintError("The parameter value specified using \"-m, --mcsParams\" option  for parameter, %s, is empty. " % (Name))
            ParamValue = Value
        else:
            if not re.match("^(Yes|No|True|False)$", Value, re.I):
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: Yes No True False" % (Value, Name))
            ParamValue = False
            if re.match("^(Yes|True)$", Value, re.I):
                ParamValue = True
        
        # Set value...
        OptionsInfo["MCSParams"][ParamName] = ParamValue

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["CoreScaffold"] = Options["--coreScaffold"]
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"], Options["--infile"], Options["--outfile"])
    
    TextOutFileMode = False
    TextOutFileDelim = ""
    
    if MiscUtil.CheckFileExt(Options["--outfile"], "csv"):
        TextOutFileMode = True
        TextOutFileDelim = ","
    elif MiscUtil.CheckFileExt(Options["--outfile"], "tsv txt"):
        TextOutFileMode = True
        TextOutFileDelim = "\t"
    
    OptionsInfo["TextOutFileMode"] = TextOutFileMode
    OptionsInfo["TextOutFileDelim"] = TextOutFileDelim

    TextOutQuote = False
    if re.match("^auto$", Options["--quote"], re.I):
        if MiscUtil.CheckFileExt(Options["--outfile"], "csv"):
            TextOutQuote = True
    else:
        if re.match("^yes$", Options["--quote"], re.I):
            TextOutQuote = True
    OptionsInfo["TextOutQuote"] = TextOutQuote
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]

    RemoveUnmatchedMode = False
    UnmatchedOutfile = None
    if re.match("^yes$", Options["--removeUnmatched"], re.I):
        RemoveUnmatchedMode = True
        FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Outfile"])
        UnmatchedOutfile = "%sUnmatched.%s" % (FileName, FileExt)
    OptionsInfo["RemoveUnmatchedMode"] = RemoveUnmatchedMode
    OptionsInfo["UnmatchedOutfile"] = UnmatchedOutfile

    OptionsInfo["SpecifiedDecompositionParams"] = Options["--decompositionParams"]
    ProcessDecompositionParameters()
    
    OptionsInfo["SpecifiedMCSParams"] = Options["--mcsParams"]
    ProcessMCSParameters()
    
    SMARTSOrSMILESCoreScaffold = ""
    SMARTSOrSMILESCoreScaffoldList = []
    if not re.match("^none$", Options["--smartsOrSmilesCoreScaffold"], re.I) or len(Options["--smartsOrSmilesCoreScaffold"]):
        if re.match("^(BySMARTS|BySMILES)$", Options["--coreScaffold"], re.I):
            SMARTSOrSMILESCoreScaffold = re.sub(" ", "", Options["--smartsOrSmilesCoreScaffold"])
            if not SMARTSOrSMILESCoreScaffold:
                MiscUtil.PrintError("A non empty value must be specified for \"-s, --smartsOrSmilesCoreScaffold\" during %s value of \"-c, --coreScaffold\" option " % (Options["--coreScaffold"]))
            SMARTSOrSMILESCoreScaffoldList = SMARTSOrSMILESCoreScaffold.split(",")
    OptionsInfo["SMARTSOrSMILESCoreScaffold"] = SMARTSOrSMILESCoreScaffold
    OptionsInfo["SMARTSOrSMILESCoreScaffoldList"] = SMARTSOrSMILESCoreScaffoldList
    
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
    
    MiscUtil.ValidateOptionTextValue("-c, --coreScaffold", Options["--coreScaffold"], "ByMCS BySMARTS BySMILES")
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi csv tsv txt")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd csv tsv txt")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    
    if re.match("^none$", Options["--smartsOrSmilesCoreScaffold"], re.I) or not len(Options["--smartsOrSmilesCoreScaffold"]):
        if re.match("^(BySMARTS|BySMILES)$", Options["--coreScaffold"], re.I):
            MiscUtil.PrintError("A non empty value must be specified for \"-s, --smartsOrSmilesCoreScaffold\" during %s value of \"-c, --coreScaffold\" option " % (Options["--coreScaffold"]))
    else:
        if not re.match("^(BySMARTS|BySMILES)$", Options["--coreScaffold"], re.I):
            MiscUtil.PrintError("%s value of \"-s, --smartsOrSmilesCoreScaffold\" is not allowed during %s value of \"-c, --coreScaffold\" option " % (Options["--smartsOrSmilesCoreScaffold"], Options["--coreScaffold"]))
    
    MiscUtil.ValidateOptionTextValue("-q, --quote", Options["--quote"], "yes no auto")
    MiscUtil.ValidateOptionTextValue("-r, --removeUnmatched", Options["--removeUnmatched"], "yes no")
    
# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitPerformRGroupDecomposition.py - Perform R group decomposition analysis

Usage:
    RDKitPerformRGroupDecomposition.py [--coreScaffold <ByMCS, BySMARTS or BySMILES>]
                                       [--decompositionParams <Name,Value,...>]
                                       [--infileParams <Name,Value,...>] [--mcsParams <Name,Value,...>]
                                       [--outfileParams <Name,Value,...>] [--overwrite] [--quote <yes or no>]
                                       [--removeUnmatched <yes or no>] [--smartsOrSmilesCoreScaffold <text>]
                                       [-w <dir>] -i <infile> -o <outfile> 
    RDKitPerformRGroupDecomposition.py -h | --help | -e | --examples

Description:
    Perform R group decomposition for a set of molecules in a series containing
    a common core scaffold. The core scaffold is identified by a SMARTS string,
    SMILES string, or using maximum common substructure (MCS) search.
    Multiple core scaffolds may be specified using SMARTS or SMILES strings for
    set of molecules corresponding to multiple series.

    The core scaffolds along with appropriate R groups are written out as SMILES
    strings to a SD or text file. The unmatched molecules without any specified
    core scaffold are written to a different output file.

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd), SMILES (.smi,
    .txt, .csv, .tsv)

    The supported output file formats are: SD (.sdf, .sd), CSV/TSV (.csv, .tsv, .txt)

Options:
    -c, --coreScaffold <ByMCS, BySMARTS or BySMILES>  [default: ByMCS]
        Specify a core scaffold for a set of molecules in a series. The core scaffold
        is identified by an explicit SMARTS string, SMILES string, or using maximum
        common substructure (MCS) search. Multiple core scaffolds may be specified
        using SMARTS or SMILES strings for set of molecules corresponding to multiple
        series.
    -d, --decompositionParams <Name,Value,...>  [default: auto]
        Parameter values to use during R group decomposition for a series of molecules.
        In general, it is a comma delimited list of parameter name and value pairs. The
        supported parameter names along with their default values are shown below:
            
            RGroupCoreAlignment,MCS, RGroupMatching,GreedyChunks,chunkSize,5,
            matchOnlyAtRGroups,no,removeHydrogenOnlyGroups,yes,
            removeHydrogensPostMatch,no
            
        A brief description of each supported parameter taken from  RDKit documentation,
        along with their possible values, is as follows.
        
        RGroupCoreAlignment - Mapping of core labels:
        
            MCS - Map core labels to each other using MCS
            None - No mapping
        
        RGroupMatching: Greedy, GreedyChunks, Exhaustive
        
        matchOnlyAtRGroups - Allow R group decomposition only at specified R groups.
        Possible values: yes, no.
        
        removeHydrogenOnlyGroups - Remove all R groups that only have hydrogens.
        Possible values: yes, no.
        
        removeHydrogensPostMatch - Remove all hydrogens from the output molecules.
        Possible values: yes, no.
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
    -m, --mcsParams <Name,Value,...>  [default: auto]
        Parameter values to use for identifying a maximum common substructure
        (MCS) in a series of molecules. In general, it is a comma delimited list of
        parameter name and value pairs. The supported parameter names along with
        their default values are shown below:
            
            atomCompare,CompareElements,bondCompare,CompareOrder,
            maximizeBonds,yes,matchValences,yes,matchChiralTag,no,
            minNumAtoms,1,minNumBonds,0,ringMatchesRingOnly,yes,
            completeRingsOnly,yes,threshold,1.0,timeOut,3600,seedSMARTS,none
            
        Possible values for atomCompare: CompareAny, CompareElements,
        CompareIsotopes. Possible values for bondCompare: CompareAny,
        CompareOrder, CompareOrderExact.
        
        A brief description of MCS parameters taken from RDKit documentation is
        as follows:
            
            atomCompare - Controls match between two atoms
            bondCompare - Controls match between two bonds
            maximizeBonds - Maximize number of bonds instead of atoms
            matchValences - Include atom valences in the MCS match
            matchChiralTag - Include atom chirality in the MCS match
            minNumAtoms - Minimum number of atoms in the MCS match
            minNumBonds - Minimum number of bonds in the MCS match
            ringMatchesRingOnly - Ring bonds only match other ring bonds
            completeRingsOnly - Partial rings not allowed during the match
            threshold - Fraction of the dataset that must contain the MCS
            seedSMARTS - SMARTS string as the seed of the MCS
            timeout - Timeout for the MCS calculation in seconds
            
    -o, --outfile <outfile>
        Output file name.
    --outfileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for writing
        molecules to files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD: compute2DCoords,auto,kekulize,no
            SMILES: kekulize,no,smilesIsomeric,yes
            
        Default value for compute2DCoords: yes for SMILES input file; no for all other
        file types. The kekulize and smilesIsomeric parameters are also used during
        generation of SMILES strings for CSV/TSV files.
    --overwrite
        Overwrite existing files.
    -q, --quote <yes or no>  [default: auto]
        Quote SMILES strings and molecule names before writing them out to text
        files. Possible values: yes or no. Default: yes for CSV (.csv) text files; no for
        TSV (.tsv) and TXT (.txt) text files.
    -r, --removeUnmatched <yes or no>  [default: no]
        Remove unmatched molecules containing no specified core scaffold from the
        output file and write them to a different output file.
    -s, --smartsOrSmilesCoreScaffold <text>  [default: none]
        SMARTS or SMILES string to use for core scaffold during 'SMARTS' or 'SMILES'
        value of '-c, --coreScaffold' option. Multiple core scaffolds may be specified using a
        comma delimited set of SMARTS or SMILES strings.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To perform R group decomposition for a set of molecules in a series using MCS
    to identify a core scaffold and write out a CSV file containing R groups, type:

        % RDKitPerformRGroupDecomposition.py -i SampleSeriesD3R.smi
          -o SampleSeriesD3ROut.csv

    To perform R group decomposition for a set of molecules in a series using a
    specified core scaffold and write out a SD file containing R groups, type:

        % RDKitPerformRGroupDecomposition.py  -c BySMARTS
          -s "Nc1nccc(-c2cnc(CNCc3ccccc3)c2)n1" -i SampleSeriesD3R.smi
          -o SampleSeriesD3ROut.sdf

    To perform R group decomposition for a set of molecules in a series using MCS
    to identify a core scaffold and write out CSV files containing matched and
    unmatched molecules without quoting values, type:

        % RDKitPerformRGroupDecomposition.py -c ByMCS -r yes -q no
          -i SampleSeriesD3R.sdf -o SampleSeriesD3ROut.csv

    To perform R group decomposition for a set of molecules in multiple series using
    specified core scaffolds and write out a TSV file containing R groups, type:

        % RDKitPerformRGroupDecomposition.py  -c BySMARTS
          -s "Nc1nccc(-c2cnc(CNCc3ccccc3)c2)n1,[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:
          [#6]:1" -i SampleMultipleSeriesD3R.smi -o
          SampleMultipleSeriesD3ROut.tsv

    To perform R group decomposition for a set of molecules in a CSV SMILES file,
    SMILES strings in  olumn 1, name in column 2, and write out a CSV file containing
    R groups, type:
     
        % RDKitPerformRGroupDecomposition.py --infileParams 
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSeriesD3R.smi -o SampleSeriesD3ROut.csv

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitSearchFunctionalGroups.py, RDKitSearchSMARTS.py

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
