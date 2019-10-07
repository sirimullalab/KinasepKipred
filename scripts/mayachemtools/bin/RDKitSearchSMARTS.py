#!/bin/env python
#
# File: RDKitSearchSMARTS.py
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
    PerformSearch()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformSearch():
    """Perform search using specified SMARTS pattern."""
    
    Infile = OptionsInfo["Infile"]
    Outfile = OptionsInfo["Outfile"]
    
    CountMode = OptionsInfo["CountMode"]
    NegateMatch = OptionsInfo["NegateMatch"]
    UseChirality = OptionsInfo["UseChirality"]

    # Set up a pattern molecule...
    PatternMol = Chem.MolFromSmarts(OptionsInfo["Pattern"])
    
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

        MolMatched = Mol.HasSubstructMatch(PatternMol, useChirality = UseChirality)
        if MolMatched != NegateMatch:
            MatchCount += 1
            if not CountMode:
                if Compute2DCoords:
                    AllChem.Compute2DCoords(Mol)
                Writer.write(Mol)
    
    if Writer is not None:
        Writer.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))

    if NegateMatch:
        MiscUtil.PrintInfo("\nTotal number of not matched molecules: %d" % MatchCount)
    else:
        MiscUtil.PrintInfo("\nTotal number of matched molecules: %d" % MatchCount)

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
    
    OptionsInfo["Pattern"] = Options["--pattern"]
    
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

def ValidateOptions():
    """Validate option values"""
    
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
        
    MiscUtil.ValidateOptionTextValue("-n, --negate", Options["--negate"], "yes no")
    
    PatternMol = Chem.MolFromSmarts(Options["--pattern"])
    if PatternMol is None:
        MiscUtil.PrintError("The value specified, %s, using option \"-p, --pattern\" is not a valid SMARTS: Failed to create pattern molecule" % Options["--pattern"])
    
    MiscUtil.ValidateOptionTextValue("--useChirality", Options["--useChirality"], "yes no")

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitSearchSMARTS.py - Perform a substructure search using SMARTS pattern

Usage:
    RDKitSearchSMARTS.py  [--infileParams <Name,Value,...>] [--mode <retrieve or count>]
                          [--negate <yes or no>] [--outfileParams <Name,Value,...>] [--overwrite]
                          [--useChirality <yes or no>] [-w <dir>] [-o <outfile>] -p <SMARTS> -i <infile>
    RDKitSearchSMARTS.py -h | --help | -e | --examples

Description:
    Perform a substructure search in an input file using specified SMARTS pattern and
    write out the matched molecules to an output file or simply count the number
    of matches.

    The supported input file formats are: SD (.sdf, .sd), SMILES (.smi., csv, .tsv, .txt)

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
            
            SD, MOL: removeHydrogens,yes,sanitize,yes,strictParsing,yes
            SMILES: smilesColumn,1,smilesNameColumn,2,smilesDelimiter,space,
                smilesTitleLine,auto,sanitize,yes
            
        Possible values for smilesDelimiter: space, comma or tab.
    -m, --mode <retrieve or count>  [default: retrieve]
        Specify whether to retrieve and write out matched molecules to an output
        file or simply count the number of matches.
    -n, --negate <yes or no>  [default: no]
        Specify whether to find molecules not matching the specified SMARTS pattern.
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
    -p, --pattern <SMARTS>  [default: none]
        SMARTS pattern for performing search.
    -u, --useChirality <yes or no>  [default: no]
        Use stereochemistry information for SMARTS search.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To retrieve molecules containing the substructure corresponding to a specified
    SMARTS pattern and write out a SMILES file, type: 

        % RDKitSearchSMARTS.py -p 'c1ccccc1' -i Sample.smi -o SampleOut.smi

    To only count the number of molecules containing the substructure corresponding
    to a specified SMARTS pattern without writing out any file, type: 

        % RDKitSearchSMARTS.py -m count -p 'c1ccccc1' -i Sample.smi

    To count the number of molecules in a SD file not containing the substructure
    corresponding to a specified SMARTS pattern and write out a SD file, type: 

        % RDKitSearchSMARTS.py -n yes -p 'c1ccccc1' -i Sample.sdf -o SampleOut.sdf

    To retrieve molecules containing the substructure corresponding to a specified
    SMARTS pattern from a CSV SMILES file, SMILES strings in column 1, name in
    and write out a SD file, type: 

        % RDKitSearchSMARTS.py -p 'c1ccccc1' --infileParams
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitFilterPAINS.py, RDKitSearchFunctionalGroups.py 

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
