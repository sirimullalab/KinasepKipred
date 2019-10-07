#!/bin/env python
#
# File: RDKitEnumerateStereoisomers.py
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
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
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
    PerformEnumeration()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformEnumeration():
    """Enumerate stereoisomers."""
    
    Infile = OptionsInfo["Infile"]
    Outfile = OptionsInfo["Outfile"]
    
    # Setup a molecule reader...
    MiscUtil.PrintInfo("\nProcessing file %s..." % Infile)
    Mols  = RDKitUtil.ReadMolecules(Infile, **OptionsInfo["InfileParams"])
    
    # Set up a molecule writer...
    Writer = None
    Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
    MiscUtil.PrintInfo("Generating file %s...\n" % Outfile)

    # Setup stereo enumeration options...
    StereoOptions = StereoEnumerationOptions(tryEmbedding = OptionsInfo["DiscardNonPhysical"], onlyUnassigned = OptionsInfo["UnassignedOnly"], maxIsomers = OptionsInfo["MaxIsomers"])
    
    # Process molecules...
    MolCount = 0
    ValidMolCount = 0

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
        
        MolName = RDKitUtil.GetMolName(Mol, MolCount)
        
        # Generate and process stereoisomers...
        StereoisomersMols = EnumerateStereoisomers(Mol, options = StereoOptions)
        IsomerCount = 0
        for IsomerMol in StereoisomersMols:
            IsomerCount += 1
            
            # Set isomer mol name...
            IsomerMolName = "%s_Isomer%d" % (MolName, IsomerCount)
            IsomerMol.SetProp("_Name", IsomerMolName)
            
            if Compute2DCoords:
                AllChem.Compute2DCoords(IsomerMol)
                
            Writer.write(IsomerMol)
        
        MiscUtil.PrintInfo("Number of stereoisomers written for %s: %d" % (MolName, IsomerCount))
            
    if Writer is not None:
        Writer.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["DiscardNonPhysical"] = True
    if re.match("^no$", Options["--discardNonPhysical"], re.I):
        OptionsInfo["DiscardNonPhysical"] = False
    
    OptionsInfo["Mode"] = Options["--mode"]
    UnassignedOnly = True
    if re.match("^All$", Options["--mode"], re.I):
        UnassignedOnly = False
    OptionsInfo["UnassignedOnly"] = UnassignedOnly
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"], Options["--infile"], Options["--outfile"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]

    OptionsInfo["MaxIsomers"] = int(Options["--maxIsomers"])

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
    
    MiscUtil.ValidateOptionTextValue("-d, --discardNonPhysical", Options["--discardNonPhysical"], "yes no")
    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "UnassignedOnly All")
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi txt csv tsv")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd smi")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
        
    MiscUtil.ValidateOptionIntegerValue("--maxIsomers", Options["--maxIsomers"], {">=": 0})

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitEnumerateStereoisomers.py - Enumerate stereoisomers of molecules

Usage:
    RDKitEnumerateStereoisomers.py [--discardNonPhysical <yes or no>]
                                [--infileParams <Name,Value,...>] [--mode <UnassignedOnly or All>]
                                [--maxIsomers <number>] [--outfileParams <Name,Value,...>] 
                                [--overwrite] [-w <dir>] -i <infile> -o <outfile> 
    RDKitEnumerateStereoisomers.py -h | --help | -e | --examples

Description:
    Perform a combinatorial enumeration of stereoisomers for molecules around all
    or unassigned chiral atoms and bonds.

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd), SMILES (.smi,
    .csv, .tsv, .txt)

    The supported output file format are: SD (.sdf, .sd), SMILES (.smi)

Options:
    -d, --discardNonPhysical <yes or no>  [default: yes]
        Discard stereoisomers with non-physical structures. Possible values: yes or no.
        The non-physical nature of a stereoisomer is determined by embedding the
        structure to generate a conformation for the stereoisomer using standard
        distance geometry methodology.
        
        A word to the wise from RDKit documentation: this is computationally expensive
        and uses a heuristic that could result in loss of stereoisomers.
    -e, --examples
        Print examples.
    -m, --mode <UnassignedOnly or All>  [default: UnassignedOnly]
        Enumerate unassigned or all chiral centers. The chiral atoms and bonds with
        defined stereochemistry are preserved.
    --maxIsomers <number>  [default: 50]
        Maximum number of stereoisomers to generate for each molecule. A  value of zero
        indicates generation of all possible steroisomers.
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
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To enumerate only unassigned atom and bond chiral centers along with discarding
    of non-physical structures, keeping a maximum of 50 stereoisomers for each molecule,
    and write out a SMILES file, type:

        % RDKitEnumerateStereoisomers.py  -i Sample.smi -o SampleOut.smi

    To enumerate only unassigned atom and bond chiral centers along with discarding
    any non-physical structures, keeping a maximum of 250 stereoisomers for a molecule,
    and write out a SD file, type:

        % RDKitEnumerateStereoisomers.py  --maxIsomers 0 -i Sample.smi
           --maxIsomers 250 -o SampleOut.sdf

    To enumerate all possible assigned and unassigned atom and bond chiral centers,
    without discarding any non-physical structures, keeping a maximum of 500 
    stereoisomers for a molecule, and write out a SD file, type:

        % RDKitEnumerateStereoisomers.py  -d no -m all --maxIsomers 500
          -i Sample.smi -o SampleOut.sdf

    To enumerate only unassigned atom and bond chiral centers along with discarding
    of non-physical structures, keeping a maximum of 50 stereoisomers for each molecule
    in a CSV SMILES file, SMILES strings in column 1, name in column 2, and write out a
    SD file with kekulization, type:

        % RDKitEnumerateStereoisomers.py  --infileParams 
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes,
          kekulize,yes" -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitEnumerateCompoundLibrary.py, RDKitGenerateConformers.py,
    RDKitGenerateMolecularFrameworks.py

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
