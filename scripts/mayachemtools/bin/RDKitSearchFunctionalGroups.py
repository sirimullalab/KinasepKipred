#!/bin/env python
#
# File: RDKitSearchFunctionalGroups.py
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
    from rdkit.Chem import FunctionalGroups
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

FunctionalGroupsMap = {}

def main():
    """Start execution of the script"""
    
    MiscUtil.PrintInfo("\n%s (RDK v%s; %s): Starting...\n" % (ScriptName, rdBase.rdkitVersion, time.asctime()))
    
    (WallClockTime, ProcessorTime) = MiscUtil.GetWallClockAndProcessorTime()
    
    # Retrieve command line arguments and options...
    RetrieveOptions()
    
    # Process and validate command line arguments and options...
    ProcessOptions()
    
    # Perform actions required by the script...
    PerformFunctionalGroupsSearch()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformFunctionalGroupsSearch():
    """Retrieve functional groups information and perform search."""

    ProcessFunctionalGroupsInfo()
    PerformSearch()

def PerformSearch():
    """Perform search using SMARTS pattern for specified functional groups."""

    Infile = OptionsInfo["Infile"]
    Outfile = OptionsInfo["Outfile"]

    Groups = OptionsInfo["SpecifiedFunctionalGroups"]
    GroupsNegateMatch = OptionsInfo["SpecifiedFunctionalGroupsNegateMatch"]
    GroupsPatternMols = OptionsInfo["SpecifiedFunctionalGroupsPatternMols"]
    
    GroupsOutfiles = OptionsInfo["SpecifiedFunctionalGroupsOutfiles"]
    GroupsCount = len(Groups)

    CombineMatchResults = OptionsInfo["CombineMatchResults"]
    AndCombineOperatorMode = OptionsInfo["AndCombineOperatorMode"]
    
    CountMode = OptionsInfo["CountMode"]
    UseChirality = OptionsInfo["UseChirality"]

    # Setup a molecule reader...
    MiscUtil.PrintInfo("\nProcessing file %s..." % Infile)
    Mols  = RDKitUtil.ReadMolecules(Infile, **OptionsInfo["InfileParams"])
    
    # Set up  molecule writers...
    Writer = None
    GroupOutfilesWriters = []
    if not CountMode:
        Writer,  GroupOutfilesWriters = SetupMoleculeWriters(CombineMatchResults, Outfile, GroupsOutfiles)
    
    # Initialize pattern mols match count and status...
    GroupsPatternMolsMatchCount = [0] * GroupsCount
    GroupsPatternMolsMatchStatus = [False] * GroupsCount

    # Process molecules...
    MolCount = 0
    ValidMolCount = 0
    MatchCount = 0

    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    
    for Mol in Mols:
        MolCount += 1
        
        if Mol is None:
            continue
        
        if RDKitUtil.IsMolEmpty(Mol):
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
            continue
        
        ValidMolCount += 1

        # Match pattern mols...
        for GroupIndex in range(0, GroupsCount):
            GroupsPatternMolsMatchStatus[GroupIndex] = DoesPatternMolMatch(GroupsPatternMols[GroupIndex], Mol, UseChirality, GroupsNegateMatch[GroupIndex])
            if GroupsPatternMolsMatchStatus[GroupIndex]:
                GroupsPatternMolsMatchCount[GroupIndex] += 1
        
        # Match mol against all specified criteria...
        MolMatched = DoesMolMeetSpecifiedMatchCriteria(GroupsPatternMolsMatchStatus, CombineMatchResults, AndCombineOperatorMode)
        if MolMatched:
            MatchCount += 1
            
        # Nothing to write...
        if CountMode or (not MolMatched):
            continue
        
        # Write out matched molecules...
        if Compute2DCoords:
            AllChem.Compute2DCoords(Mol)
        
        if CombineMatchResults:
            Writer.write(Mol)
        else:
            for GroupIndex in range(0, GroupsCount):
                if GroupsPatternMolsMatchStatus[GroupIndex]:
                    GroupOutfilesWriters[GroupIndex].write(Mol)
    
    if Writer is not None:
        Writer.close()
    for GroupOutfileWriter in GroupOutfilesWriters:
        GroupOutfileWriter.close()
        
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))

    MiscUtil.PrintInfo("\nTotal number of molecules matched against specified match criteria: %d" % MatchCount)
    
    MiscUtil.PrintInfo("\nNumber of molecuels matched against individual functional groups:")
    MiscUtil.PrintInfo("FunctionalGroupName,MatchCount")
    
    for GroupIndex in range(0, GroupsCount):
        GroupName = Groups[GroupIndex]
        if GroupsNegateMatch[GroupIndex]:
            GroupName = '!' + GroupName
        GroupMatchCount = GroupsPatternMolsMatchCount[GroupIndex]
        MiscUtil.PrintInfo("%s,%d" % (GroupName, GroupMatchCount))

def DoesMolMeetSpecifiedMatchCriteria(GroupsPatternMolsMatchStatus,  CombineMatchResults, AndCombineOperatorMode):
    """Match molecule using specified match criteia."""

    if CombineMatchResults and AndCombineOperatorMode:
        # Must match all specified SMARTS
        Status = True
        for MatchStatus in GroupsPatternMolsMatchStatus:
            if not MatchStatus:
                Status = False
                break
    else:
        # One match is enough...
        Status = False
        for MatchStatus in GroupsPatternMolsMatchStatus:
            if MatchStatus:
                Status = True
                break
    
    return Status
    
def SetupMoleculeWriters(CombineMatchResults, Outfile, GroupsOutfiles):
    """Set up molecule writers for output files."""

    Writer = None
    GroupOutfilesWriters = []
    
    if CombineMatchResults:
        Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
        if Writer is None:
            MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
        MiscUtil.PrintInfo("Generating file %s..." % Outfile)
    else:
        for GroupOutfile in GroupsOutfiles:
            GroupOutfileWriter = RDKitUtil.MoleculesWriter(GroupOutfile, **OptionsInfo["OutfileParams"])
            if GroupOutfileWriter is None:
                MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Writer)
            GroupOutfilesWriters.append(GroupOutfileWriter)
        
        GroupsCount = len(GroupsOutfiles)
        if GroupsCount > 4:
            MiscUtil.PrintInfo("Generating %d output files with the following file name format: %s<GroupName>.%s" % (GroupsCount, OptionsInfo["OutfileBasename"], OptionsInfo["OutfileExt"]))
        else:
            Delmiter = ', '
            OutfileNames = Delmiter.join(GroupsOutfiles)
            MiscUtil.PrintInfo("Generating %d output files: %s..." % (GroupsCount, OutfileNames))

    return (Writer, GroupOutfilesWriters)
    
def DoesPatternMolMatch(PatternMol, Mol, UseChirality, NegateMatch):
    """Perform a substructure match for the presence of pattern molecule in a molecule."""

    MolMatched = Mol.HasSubstructMatch(PatternMol, useChirality = UseChirality)
    if NegateMatch:
        if MolMatched:
            MolMatched = False
        else:
            MolMatched = True
    
    return MolMatched
    
def ProcessFunctionalGroupsInfo():
    """Process functional groups information."""
    
    RetrieveFunctionalGroupsInfo()
    ProcessSpecifiedFunctionalGroups()

    SetupFunctionalGroupsSMARTSPatterns()
    SetupFunctionalGroupsOutputFileNames()
    
def ProcessSpecifiedFunctionalGroups():
    """Process and validate specified functional groups"""
    
    OptionsInfo["SpecifiedFunctionalGroups"] = []
    OptionsInfo["SpecifiedFunctionalGroupsNegateMatch"] = []
    
    if re.match("^All$", OptionsInfo["FunctionalGroups"], re.I):
        OptionsInfo["SpecifiedFunctionalGroups"] = FunctionalGroupsMap['Names']
        OptionsInfo["SpecifiedFunctionalGroupsNegateMatch"] = [False] * len(OptionsInfo["SpecifiedFunctionalGroups"])
        return
    
    # Set up a map of valid group names for checking specified group names...
    CanonicalGroupNameMap = {}
    for GroupName in FunctionalGroupsMap['Names']:
        CanonicalGroupNameMap[GroupName.lower()] = GroupName
        
    # Parse and validate specified names...
    GroupNames = re.sub(" ", "", OptionsInfo["FunctionalGroups"])
    if not GroupNames:
        MiscUtil.PrintError("No functional group name specified for \"-f, --functionalGroups\" option")
    
    SpecifiedFunctionalGroups = []
    SpecifiedNegateMatchStatus = []
    
    for GroupName in GroupNames.split(","):
        CanonicalGroupName = GroupName.lower()
        NegateMatchStatus = False
        if re.match("^!", CanonicalGroupName, re.I):
            NegateMatchStatus = True
            CanonicalGroupName = re.sub("^!", "", CanonicalGroupName)
        if CanonicalGroupName in CanonicalGroupNameMap:
            SpecifiedFunctionalGroups.append(CanonicalGroupNameMap[CanonicalGroupName])
            SpecifiedNegateMatchStatus.append(NegateMatchStatus)
        else:
            MiscUtil.PrintWarning("The functional group name, %s, specified using \"-f, --functionalGroups\" option is not a valid name." % (GroupName))

    if not len(SpecifiedFunctionalGroups):
        MiscUtil.PrintError("No valid functional group names specified for \"-f, --functionalGroups\" option")
        
    OptionsInfo["SpecifiedFunctionalGroups"] = SpecifiedFunctionalGroups
    OptionsInfo["SpecifiedFunctionalGroupsNegateMatch"] = SpecifiedNegateMatchStatus

def SetupFunctionalGroupsSMARTSPatterns():
    """Setup SMARTS patterns for specified functional groups."""

    OptionsInfo["SpecifiedFunctionalGroupsSMARTSPatterns"] = []
    OptionsInfo["SpecifiedFunctionalGroupsPatternMols"] = []
    
    for Name in OptionsInfo["SpecifiedFunctionalGroups"]:
        SMARTSPattern = FunctionalGroupsMap['SMARTSPattern'][Name]
        PatternMol = Chem.MolFromSmarts(SMARTSPattern)
        if PatternMol is None:
            MiscUtil.PrintError("Failed to parse SMARTS pattern, %s, for function group, %s" % (SMARTSPattern, Name))
        
        OptionsInfo["SpecifiedFunctionalGroupsSMARTSPatterns"].append(SMARTSPattern)
        OptionsInfo["SpecifiedFunctionalGroupsPatternMols"].append(PatternMol)

def SetupFunctionalGroupsOutputFileNames():
    """Setup output file names for specified functional group names."""

    OptionsInfo["SpecifiedFunctionalGroupsOutfiles"] = []
    
    if OptionsInfo["CountMode"]:
        # No need of any output file...
        return
    
    if OptionsInfo["CombineMatchResults"]:
        # No need of output files for specified functional groups...
        return
    
    OutfileBasename = OptionsInfo["OutfileBasename"]
    OutfileExt = OptionsInfo["OutfileExt"]
    SpecifiedFunctionalGroupsOutfiles = []

    GroupsCount = len(OptionsInfo["SpecifiedFunctionalGroups"])
    for GroupIndex in range(0, GroupsCount):
        GroupName = OptionsInfo["SpecifiedFunctionalGroups"][GroupIndex]
        if OptionsInfo["SpecifiedFunctionalGroupsNegateMatch"][GroupIndex]:
            GroupName = "Not" + GroupName
        GroupName = re.sub("\.", "", GroupName)
        
        GroupOutfile = "%s%s.%s" % (OutfileBasename, GroupName, OutfileExt)
        SpecifiedFunctionalGroupsOutfiles.append(GroupOutfile)

    OptionsInfo["SpecifiedFunctionalGroupsOutfiles"] = SpecifiedFunctionalGroupsOutfiles
    
def RetrieveFunctionalGroupsInfo():
    """Retrieve functional groups information"""

    MiscUtil.PrintInfo("\nRetrieving data from default RDKit functional groups hierarchy file Functional_Group_Hierarchy.txt...")
    
    FunctionalGroupNamesFile = OptionsInfo["GroupNamesFile"]
    FunctionalGroupsNodes = FunctionalGroups.BuildFuncGroupHierarchy(FunctionalGroupNamesFile)

    FunctionalGroupsMap['Names'] = []
    FunctionalGroupsMap['SMARTSPattern'] = {}
    
    RetrieveDataFromFunctionalGroupsHierarchy(FunctionalGroupsNodes)

    if not len(FunctionalGroupsMap['Names']):
        MiscUtil.PrintError("Failed to retrieve any functional group names and SMARTS patterns...")
        
    MiscUtil.PrintInfo("Total number of functional groups present functional group hierarchy: %d" % (len(FunctionalGroupsMap['Names'])))
    
def RetrieveDataFromFunctionalGroupsHierarchy(FGNodes):
    """Retrieve functional groups data by recursively visiting functional group nodes."""
    
    for FGNode in FGNodes:
        Name = FGNode.label
        SMARTSPattern = FGNode.smarts

        if Name in FunctionalGroupsMap['SMARTSPattern']:
            MiscUtil.PrintWarning("Ignoring duplicate functional group name: %s..." % Name)
        else:
            FunctionalGroupsMap['Names'].append(Name)
            FunctionalGroupsMap['SMARTSPattern'][Name] = SMARTSPattern

        RetrieveDataFromFunctionalGroupsHierarchy(FGNode.children)

def ListFunctionalGroupsInfo():
    """List functional groups information"""

    MiscUtil.PrintInfo("\nListing available functional groups names and SMARTS patterns...")
    MiscUtil.PrintInfo("\nFunctionalGroupName\tSMARTSPattern")
    
    for Name in sorted(FunctionalGroupsMap['Names']):
        SMARTSPattern = FunctionalGroupsMap['SMARTSPattern'][Name]
        MiscUtil.PrintInfo("%s\t%s" % (Name, SMARTSPattern))
    
    MiscUtil.PrintInfo("")

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["CombineMatches"] = Options["--combineMatches"]
    
    OptionsInfo["CombineMatchResults"] = True
    if re.match("^No$", Options["--combineMatches"], re.I):
        OptionsInfo["CombineMatchResults"] = False
        if Options["--outfile"]:
            FileDir, FileName, FileExt = MiscUtil.ParseFileName(Options["--outfile"])
            OptionsInfo["OutfileBasename"] = FileName
            OptionsInfo["OutfileExt"] = FileExt
    
    OptionsInfo["CombineOperator"] = Options["--combineOperator"]
    OptionsInfo["AndCombineOperatorMode"] = True
    if re.match("^or$", Options["--combineOperator"], re.I):
        OptionsInfo["AndCombineOperatorMode"] = False
    
    OptionsInfo["GroupNamesFile"] = None
    if not re.match("^auto$", Options["--groupNamesFile"], re.I):
        OptionsInfo["GroupNamesFile"] = Options["--groupNamesFile"]
        
    OptionsInfo["FunctionalGroups"] = Options["--functionalGroups"]
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"], Options["--infile"], Options["--outfile"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    
    OptionsInfo["CountMode"] = False
    if re.match("^count$", Options["--mode"], re.I):
        OptionsInfo["CountMode"] = True
    
    OptionsInfo["UseChirality"] = False
    if re.match("^yes$", Options["--useChirality"], re.I):
        OptionsInfo["UseChirality"] = True

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
    
    # Handle listing of functional group information...
    if  Options and Options["--list"]:
        ProcessListFunctionalGroupsOption()
        sys.exit(0)

def ProcessListFunctionalGroupsOption():
    """Process list functional groups information."""

    # Validate and process dataFile option for listing functional groups information...
    OptionsInfo["GroupNamesFile"] = None
    if not re.match("^auto$", Options["--groupNamesFile"], re.I):
        MiscUtil.ValidateOptionFilePath("-g, --groupNamesFile", Options["--groupNamesFile"])
        OptionsInfo["GroupNamesFile"] = Options["--groupNamesFile"]
    
    RetrieveFunctionalGroupsInfo()
    ListFunctionalGroupsInfo()
    
def ValidateOptions():
    """Validate option values"""
    
    MiscUtil.ValidateOptionTextValue("-c, --combineMatches", Options["--combineMatches"], "yes no")
    MiscUtil.ValidateOptionTextValue("--combineOperator", Options["--combineOperator"], "and or")
    
    if not re.match("^auto$", Options["--groupNamesFile"], re.I):
        MiscUtil.ValidateOptionFilePath("-g, groupNamesFile", Options["--groupNamesFile"])
        
    if re.match("^none$", Options["--functionalGroups"], re.I):
        MiscUtil.PrintError("The name(s) of functional groups must be specified using \"-f, --functionalGroups\" option")
        
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd smi txt csv tsv")
    if Options["--outfile"]:
        MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd smi")
        MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
        MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    
    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "retrieve count")
    if re.match("^retrieve$", Options["--mode"], re.I):
        if not Options["--outfile"]:
            MiscUtil.PrintError("The outfile must be specified using \"-o, --outfile\" during \"retrieve\" value of \"-m, --mode\" option")
        
    MiscUtil.ValidateOptionTextValue("--useChirality", Options["--useChirality"], "yes no")

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitSearchFunctionalGroups.py - Search for functional groups using SMARTS patterns

Usage:
    RDKitSearchFunctionalGroups.py  [--combineMatches <yes or no>] [--combineOperator <and or or>]
                                           [--groupNamesFile <FileName or auto>] [--infileParams <Name,Value,...>]
                                           [--mode <retrieve or count>] [--negate <yes or no>]
                                           [--outfileParams <Name,Value,...>] [--overwrite] [--useChirality <yes or no>]
                                           [-w <dir>] [-o <outfile>] -i <infile> -f <Name1,Name2,Name3... or All>
    RDKitSearchFunctionalGroups.py [--groupNamesFile <FileName or auto>] -l | --list
    RDKitSearchFunctionalGroups.py -h | --help | -e | --examples

Description:
    Perform a substructure search in an input file using SMARTS patterns for functional
    groups and write out the matched molecules to an output file or simply count the
    number of matches.

    The SMARTS patterns for specified functional group(s) are retrieved from file,
    Functional_Group_Hierarchy.txt, available in RDKit data directory.

    The names of valid functional groups and hierarchies  are dynamically retrieved from the
    functional groups hierarchy file and are shown below:

        AcidChloride, AcidChloride.Aromatic, AcidChloride.Aliphatic
        Alcohol, Alcohol.Aromatic, Alcohol.Aliphatic
        Aldehyde, Aldehyde.Aromatic, Aldehyde.Aliphatic
        Amine, Amine.Primary, Amine.Primary.Aromatic, Amine.Primary.Aliphatic,
        Amine.Secondary, Amine.Secondary.Aromatic, Amine.Secondary.Aliphatic
        Amine.Tertiary, Amine.Tertiary.Aromatic, Amine.Tertiary.Aliphatic
        Amine.Aromatic, Amine.Aliphatic, Amine.Cyclic
        Azide, Azide.Aromatic, Azide.Aliphatic
        BoronicAcid, BoronicAcid.Aromatic, BoronicAcid.Aliphatic
        CarboxylicAcid, CarboxylicAcid.Aromatic, CarboxylicAcid.Aliphatic,
        CarboxylicAcid.AlphaAmino
        Halogen, Halogen.Aromatic, Halogen.Aliphatic
        Halogen.NotFluorine, Halogen.NotFluorine.Aliphatic,
        Halogen.NotFluorine.Aromatic
        Halogen.Bromine, Halogen.Bromine.Aliphatic, Halogen.Bromine.Aromatic,
        Halogen.Bromine.BromoKetone
        Isocyanate, Isocyanate.Aromatic, Isocyanate.Aliphatic
        Nitro, Nitro.Aromatic, Nitro.Aliphatic,
        SulfonylChloride, SulfonylChloride.Aromatic, SulfonylChloride.Aliphatic
        TerminalAlkyne

    The supported input file formats are: SD (.sdf, .sd), SMILES (.smi, .csv, .tsv, .txt)

    The supported output file formats are: SD (.sdf, .sd), SMILES (.smi)

Options:
    -c, --combineMatches <yes or no>  [default: yes]
        Combine search results for matching SMARTS patterns of specified functional groups
        against a molecule. Possible values: yes or no.
        
        The matched molecules are written to a single output file for "yes" value. Otherwise,
        multiple output files are generated, one for each functional group. The names of  
        these files correspond to a combination of the basename of the specified output file
        and the name of the functional group.
        
        No output files are generated during "count" value of "-m, --mode" option.
    --combineOperator <and or or>  [default: and]
        Logical operator to use for combining match results corresponding to specified
        functional group names before writing out a single file. This option is ignored
        during "No" value of  "-c, --combineMatches" option.
    -e, --examples
        Print examples.
    -g, --groupNamesFile <FileName or auto>  [default: auto]
        Specify a file name containing data for functional groups hierarchy or use functional
        group hierarchy file, Functional_Group_Hierarchy.txt, available in RDKit data directory.
        
        RDKit data format: Name<tab>Smarts<tab>Label<tab>RemovalReaction (optional)
        
        The format of data in local functional group hierarchy must match format of the
        data in functional group file available in RDKit data directory.
    -f, --functionalGroups <Name1,Name2,Name3... or All>  [default: none]
        Functional group names for performing substructure SMARTS search. Possible values:
        Comma delimited list of valid functional group names or All. The current set of valid
        functional group names are listed in the description section.
        
        The match results for multiple functional group names are combined using 'and'
        operator before writing them out to single file. No merging of match results takes
        place during generation of individual result files corresponding to fictional group
        names.
        
        The functional group name may be started with an exclamation mark to negate
        the match result for that fictional group.
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
    -l, --list
        List functional groups information without performing any search.
    -m, --mode <retrieve or count>  [default: retrieve]
        Specify whether to retrieve and write out matched molecules to an output
        file or simply count the number of matches.
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
    -u, --useChirality <yes or no>  [default: no]
        Use stereochemistry information for SMARTS search.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To list names of all available functional groups along with their SMARTS patterns, type:

        % RDKitSearchFunctionalGroups.py -l

    To retrieve molecule containing amine functional group and write out a SMILES file, type: 

        % RDKitSearchFunctionalGroups.py -f Amine -i Sample.smi -o SampleOut.smi

    To retrieve molecules containing amine functional group but not halogens and carboxylic
    acid functional groups and write out a SMILES file, type: 

        % RDKitSearchFunctionalGroups.py -f 'Amine,!Halogen,!CarboxylicAcid'
          -i Sample.smi -o SampleOut.smi

    To retrieve molecules containing amine, halogens or carboxylic  acid functional groups
    and write out a SMILES file, type: 

        % RDKitSearchFunctionalGroups.py -f 'Amine,Halogen,CarboxylicAcid'
          --combineOperator or -i Sample.smi -o SampleOut.smi

    To retrieve molecules containing amine and carboxylic acid functional groups defined in
    a local functional groups hierarchy file and write out individual SD files for each
    funcitonal group, type: 

        % RDKitSearchFunctionalGroups.py -f 'Amine,CarboxylicAcid' -i Sample.sdf 
          -g Custom_Functional_Group_Hierarchy.txt --combineMatches No -o SampleOut.sdf

    To count number of all functional groups in molecules without writing out an output
    files, type:

        % RDKitSearchFunctionalGroups.py -m count -f All --combineMatches no -i Sample.smi

    To retrieve molecule not containing aromatic alcohol and aromatic halogen functional
    group along with the use of chirality during substructure search and write out individual
    SMILES files for eeah functional group, type: 

        % RDKitSearchFunctionalGroups.py --combineMatches no -u yes
           -f '!Alcohol.Aromatic,!Halogen.Aromatic' -i Sample.smi -o SampleOut.smi

    To retrieve molecule containing amine functional group from a CSV SMILES file,
    SMILES strings in column 1, name in column 2, and write out a SD file, type: 

        % RDKitSearchFunctionalGroups.py -f Amine --infileParams
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitFilterPAINS.py, RDKitSearchSMARTS.py

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
