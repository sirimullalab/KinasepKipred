#!/bin/env python
#
# File: RDKitRemoveDuplicateMolecules.py
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
    RemoveDuplicates()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def RemoveDuplicates():
    """Identify and remove duplicate molecules based on canonical SMILES"""
    
    Infile = OptionsInfo["Infile"]
    Outfile = OptionsInfo["Outfile"]
    DuplicatesOutfile = OptionsInfo["DuplicatesOutfile"]
    
    CountMode = OptionsInfo["CountMode"]
    UseChirality = OptionsInfo["UseChirality"]
    
    # Setup a molecule reader...
    MiscUtil.PrintInfo("\nProcessing file %s..." % Infile)
    Mols  = RDKitUtil.ReadMolecules(Infile, **OptionsInfo["InfileParams"])
    
    # Set up a molecule writer...
    Writer = None
    DuplicatesWriter = None
    if not CountMode:
        Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
        if Writer is None:
            MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
        DuplicatesWriter = RDKitUtil.MoleculesWriter(DuplicatesOutfile, **OptionsInfo["OutfileParams"])
        if DuplicatesWriter is None:
            MiscUtil.PrintError("Failed to setup a writer for output fie %s " % DuplicatesOutfile)
        
        MiscUtil.PrintInfo("Generating files %s and %s..." % (Outfile, DuplicatesOutfile))

    # Process molecules...
    MolCount = 0
    ValidMolCount = 0
    
    UniqueMolCount = 0
    DuplicateMolCount = 0
    
    CanonicalSMILESMap = {}
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
        
        CanonicalSMILES = Chem.MolToSmiles(Mol, isomericSmiles = UseChirality, canonical = True)
        
        if Compute2DCoords:
            if not CountMode:
                AllChem.Compute2DCoords(Mol)
        
        if CanonicalSMILES in CanonicalSMILESMap:
            DuplicateMolCount += 1
            if not CountMode:
                DuplicatesWriter.write(Mol)
        else:
            UniqueMolCount += 1
            CanonicalSMILESMap[CanonicalSMILES] = CanonicalSMILES
            if not CountMode:
                Writer.write(Mol)
    
    if Writer is not None:
        Writer.close()
    
    if DuplicatesWriter is not None:
        DuplicatesWriter.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))

    MiscUtil.PrintInfo("\nTotal number of unique molecules: %d" % UniqueMolCount)
    MiscUtil.PrintInfo("Total number of duplicate molecules: %d" % DuplicateMolCount)

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
    
    # Setup outfile for writing out duplicates...
    OptionsInfo["DuplicatesOutfile"] = ""
    if not OptionsInfo["CountMode"] :
        FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Outfile"])
        OptionsInfo["DuplicatesOutfile"] = "%sDuplicates.%s" % (FileName, FileExt)

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

    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "remove count")
    if re.match("^remove$", Options["--mode"], re.I):
        if not Options["--outfile"]:
            MiscUtil.PrintError("The outfile must be specified using \"-o, --outfile\" during \"remove\" value of \"-m, --mode\" option")
        
    MiscUtil.ValidateOptionTextValue("--useChirality", Options["--useChirality"], "yes no")
    
# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitRemoveDuplicateMolecules.py - Remove duplicate molecules

Usage:
    RDKitRemoveDuplicateMolecules.py  [--infileParams <Name,Value,...>]
                              [--mode <remove or count>] [ --outfileParams <Name,Value,...> ] 
                              [--overwrite] [--useChirality <yes or no>] [-w <dir>] [-o <outfile>]  -i <infile>
    RDKitRemoveDuplicateMolecules.py -h | --help | -e | --examples

Description:
    Identify and remove duplicate molecules based on canonical SMILES strings or
    simply count the number of duplicate molecules.

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
            
            SD: removeHydrogens,yes,sanitize,yes,strictParsing,yes
            SMILES: smilesColumn,1,smilesNameColumn,2,smilesDelimiter,space,
                smilesTitleLine,auto,sanitize,yes
            
        Possible values for smilesDelimiter: space, comma or tab.
    -m, --mode <remove or count>  [default: remove]
        Specify whether to remove duplicate molecules and write out filtered molecules
        to output files or or simply count the number of duplicate molecules.
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
    -u, --useChirality <yes or no>  [default: yes]
        Use stereochemistry information for generation of canonical SMILES strings
        to identify duplicate molecules.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To remove duplicate molecules and generate output files containing unique and
    duplicate SMILES strings, type:

        % RDKitRemoveDuplicateMolecules.py -i Sample.smi -o SampleOut.smi

    To remove duplicate molecules without using stereochemistry information for
    generation of canonical SMILES and generate output files containing unique and
    duplicate SMILES strings, type:

        % RDKitRemoveDuplicateMolecules.py -u no -i Sample.sdf -o SampleOut.sdf

    To count number of unique and duplicate molecules without generating any
    output files, type:

        % RDKitRemoveDuplicateMolecules.py -m count -i Sample.sdf

    To remove duplicate molecules from a CSV SMILES file, SMILES strings in
    column 1, name in column 2, and generate output SD files containing unique and
    duplicate molecules, type:

        % RDKitRemoveDuplicateMolecules.py --infileParams 
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitRemoveSalts, RDKitSearchFunctionalGroups.py,
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
