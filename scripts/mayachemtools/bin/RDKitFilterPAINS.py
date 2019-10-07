#!/bin/env python
#
# File: RDKitFilterPAINS.py
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
import csv

# RDKit imports...
try:
    from rdkit import rdBase
    from rdkit import Chem
    from rdkit.Chem import AllChem
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
    PerformFiltering()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformFiltering():
    """Filter molecules using SMARTS specified in PAINS filter file."""
    
    Infile = OptionsInfo["Infile"]
    Outfile = OptionsInfo["Outfile"]
    
    CountMode = OptionsInfo["CountMode"]
    NegateMatch = OptionsInfo["NegateMatch"]
    PAINSMode = OptionsInfo["PAINSMode"]

    # Setup PAINS patterns and pattern mols...
    PAINSPatterns = RetrievePAINSPatterns(PAINSMode)
    PAINSPatternMols = SetupPAINSPatternMols(PAINSPatterns)
    
    # Setup a molecule reader...
    MiscUtil.PrintInfo("\nProcessing file %s..." % Infile)
    Mols  = RDKitUtil.ReadMolecules(Infile, **OptionsInfo["InfileParams"])
    
    # Set up a molecule writer...
    Writer = None
    if not CountMode:
        Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
        if Writer is None:
            MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
        MiscUtil.PrintInfo("Generating file %s..." % Outfile)

    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    
    # Process molecules...
    MolCount = 0
    ValidMolCount = 0
    FilteredCount = 0

    for Mol in Mols:
        MolCount += 1
        
        if Mol is None:
            continue
        
        if RDKitUtil.IsMolEmpty(Mol):
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
            continue
        
        ValidMolCount += 1
        
        MolMatched = DoesMoleculeContainsPAINSPattern(Mol, PAINSPatternMols)
        if MolMatched == NegateMatch:
            FilteredCount += 1
            if not CountMode:
                if Compute2DCoords:
                    AllChem.Compute2DCoords(Mol)
                Writer.write(Mol)
    
    if Writer is not None:
        Writer.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))

    MiscUtil.PrintInfo("\nTotal number of filtered molecules: %d" % FilteredCount)

def DoesMoleculeContainsPAINSPattern(Mol, PAINSPatternMols):
    """Check presence of PAINS pattern in the molecule"""

    MolMatched = False
    
    for PatternMol in PAINSPatternMols:
        if Mol.HasSubstructMatch(PatternMol, useChirality = True):
            MolMatched = True
            break
        
    return MolMatched
    
def RetrievePAINSPatterns(PAINSFilterMode):
    """Retrieve PAINS patterns for specified PAINS mode"""

    MayaChemToolsDataDir = MiscUtil.GetMayaChemToolsLibDataPath()
    PAINSFiltersFilePath = os.path.join(MayaChemToolsDataDir, "PAINSFilters.csv")
    
    MiscUtil.PrintInfo("\nRetrieving PAINS SMARTS patterns for PAINS filter type, %s, from file %s" % (PAINSFilterMode, PAINSFiltersFilePath))

    if not os.path.exists(PAINSFiltersFilePath):
        MiscUtil.PrintError("The PAINS filters file, %s, doesn't exist.\n" % (PAINSFiltersFilePath))
        
    FilterFile = open(PAINSFiltersFilePath, "r")
    if FilterFile is None:
        MiscUtil.PrintError("Couldn't open PAINS filter file: %s.\n" % (PAINSFiltersFilePath))

    # Collect all PAINS filter lines...
    HeaderLine = True
    FiltersLines = []
    for Line in FilterFile:
        Line = Line.strip()
        # Ignore comments...
        if re.match("^#", Line, re.I):
            continue
        # Ignore header line...
        if HeaderLine:
            HeaderLine = False
            continue
        FiltersLines.append(Line)
        
    # Process PAINS filter lines using csv reader...
    SMARTSPatterns = []
    
    FiltersReader = csv.reader(FiltersLines, delimiter=',', quotechar='"')
    for LineWords in FiltersReader:
        FilterType = LineWords[0]
        ID = LineWords[1]
        SMARTS = LineWords[2]

        if re.match("^All$", PAINSFilterMode, re.I) or FilterType.lower() == PAINSFilterMode.lower():
            SMARTSPatterns.append(SMARTS)
        
    FilterFile.close()

    MiscUtil.PrintInfo("Total number of PAINS SMARTS patterns: %d" % (len(SMARTSPatterns)))
    
    return SMARTSPatterns

def SetupPAINSPatternMols(PAINSPatterns):
    """Set up PAINS pattern mols for substructure search"""

    MiscUtil.PrintInfo("\nSetting up PAINS pattern molecules for performins substructure search...")
    PatternMols = []
    for Pattern in PAINSPatterns:
        PatternMol = Chem.MolFromSmarts(Pattern)
        PatternMols.append(PatternMol)
        
    return PatternMols    

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"], Options["--infile"], Options["--outfile"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]

    OptionsInfo["CountMode"] = False
    if re.match("^count$", Options["--mode"], re.I):
        OptionsInfo["CountMode"] = True
    
    OptionsInfo["NegateMatch"] = False
    if re.match("^yes$", Options["--negate"], re.I):
        OptionsInfo["NegateMatch"] = True
        
    OptionsInfo["PAINSMode"] = Options["--painsMode"]
    
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
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd smi txt csv tsv")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd smi")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])

    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "filter count")
    if re.match("^filter$", Options["--mode"], re.I):
        if not Options["--outfile"]:
            MiscUtil.PrintError("The outfile must be specified using \"-o, --outfile\" during \"filter\" value of \"-m, --mode\" option")
        
    MiscUtil.ValidateOptionTextValue("-n, --negate", Options["--negate"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("-p, --painsMode", Options["--painsMode"], "All A B C")

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitFilterPAINS.py - Filter PAINS molecules

Usage:
    RDKitFilterPAINS.py  [--infileParams <Name,Value,...>] [--mode <filter or count>]
                         [ --outfileParams <Name,Value,...> ] [--painsMode <All, A, B or C>]
                         [--negate <yes or no>] [--overwrite] [-w <dir>] -i <infile> -o <outfile>
    RDKitFilterPAINS.py -h | --help | -e | --examples

Description:
    Filter Pan-assay Interference molecules (PAINS) [ Ref 130 - 131 ] from an input
    file by performing a substructure search using SMARTS pattern specified in
    MAYACHEMTOOLS/lib/data/PAINSFilter.csv file and write out appropriate
    molecules to an output file or simply count the number of filtered molecules.

    The supported input file formats are: SD (.sdf, .sd), SMILES (.smi, .csv,
     .tsv, .txt)

    The supported output file formats are: SD (.sdf, .sd), SMILES (.smi)

Options:
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
            
            SD: removeHydrogens,yes,sanitize,yes,strictParsing,yes
            SMILES: smilesColumn,1,smilesNameColumn,2,smilesDelimiter,space,
                smilesTitleLine,auto,sanitize,yes
            
        Possible values for smilesDelimiter: space, comma or tab.
    -m, --mode <filter or count>  [default: filter]
        Specify whether to filter the matched molecules and write out the rest of the 
        molecules to an outfile or simply count the number of matched molecules
        marked for filtering.
    -n, --negate <yes or no>  [default: no]
        Specify whether to filter molecules not matching the PAINS filters specified by
        SMARTS patterns.
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
    -p, --painsMode <All, A, B, or C>  [default: All]
        Specify PAINS filter family type to used for filtering molecules. 
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To count the number of molecules not containing any substructure corresponding to
    PAINS SMARTS patterns and write out a SMILES file, type: 

        % RDKitFilterPAINS.py -i Sample.smi -o SampleOut.smi

    To only count the number of molecules not containing any substructure corresponding
    to PAINS SMARTS patterns without writing out any file file, type: 

        % RDKitFilterPAINS.py -m count -i Sample.sdf -o SampleOut.smi

    To count the number of molecules containing any substructure corresponding to
    PAINS SMARTS patterns and write out a SD file with computed 2D coordinates,
    type: 

        % RDKitFilterPAINS.py -n yes -i Sample.smi -o SampleOut.sdf

    To count the number of molecules not containing any substructure corresponding to
    PAINS SMARTS patterns family of Type A in a CSV SMILS file and write out a SD file, type: 

        % RDKitFilterPAINS.py --painsMode A --infileParams
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitSearchSMARTS.py

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
