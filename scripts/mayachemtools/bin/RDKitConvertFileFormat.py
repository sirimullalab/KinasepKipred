#!/bin/env python
#
# File: RDKitConvertFileFormat.py
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
    ConvertFileFormat()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def ConvertFileFormat():
    """Convert between  file formats"""
    
    Infile = OptionsInfo["Infile"]
    Outfile = OptionsInfo["Outfile"]
    
    # Read molecules...
    MiscUtil.PrintInfo("\nReading file %s..." % Infile)
    Mols = RDKitUtil.ReadMolecules(Infile, **OptionsInfo["InfileParams"])
    MiscUtil.PrintInfo("Total number of molecules: %d" % len(Mols))
    
    # Write molecules...
    MiscUtil.PrintInfo("\nGenerating file %s..." % Outfile)
    MolCount, ProcessedMolCount = RDKitUtil.WriteMolecules(Outfile, Mols, **OptionsInfo["OutfileParams"])
    
    MiscUtil.PrintInfo("Total number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of molecules processed: %d" % ProcessedMolCount)
    MiscUtil.PrintInfo("Number of molecules ignored: %d" % (MolCount - ProcessedMolCount))

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    # Process and setup options for RDKit functions...
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"], Options["--infile"], Options["--outfile"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]

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
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi txt csv tsv mol2 pdb")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd mol smi pdb")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    
# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitConvertFileFormat.py - Convert between molecular file formats

Usage:
    RDKitConvertFileFormat.py [--infileParams <Name,Value,...>]
                              [ --outfileParams <Name,Value,...> ] [--overwrite]
                              [-w <dir>] -i <infile> -o <outfile>
    RDKitConvertFileFormat.py -h | --help | -e | --examples

Description:
    Convert between molecular file formats.

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd), SMILES (.smi,
    .txt, .csv, .tsv), MOL2 (.mol2), PDB (.pdb)

    The supported output file formats are: SD (.sdf, .sd), SMILES (.smi), PDB (.pdb)

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
            MOL2: removeHydrogens,yes,sanitize,yes
            SMILES: smilesColumn,1,smilesNameColumn,2,smilesDelimiter,space,
                smilesTitleLine,auto,sanitize,yes
            PDB: removeHydrogens,yes,sanitize,yes
            
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
    To convert a SD file  into a isomeric SMILES file, type:

        % RDKitConvertFileFormat.py -i Sample.sdf -o SampleOut.smi

    To convert a SD file into a non isomeric SMILES file, type

        % RDKitConvertFileFormat.py --outfileParams "smilesIsomeric,no"
          -i Sample.sdf -o SampleOut.smi

    To convert a SMILES file into a SD file along with calculation of 2D
    coordinates, type:

        % RDKitConvertFileFormat.py -i Sample.smi -o SampleOut.sdf

    To convert a MDL MOL file into a PDB file, type:

        % RDKitConvertFileFormat.py -i Sample.mol -o SampleOut.pdb

    To convert a CSV SMILES file  with column headers, SMILES strings
    in column 1, and name in column 2 into a SD file containing 2D coordinates, type:

        % RDKitConvertFileFormat.py --infileParams "smilesDelimiter,comma,
          smilesTitleLine,yes,smilesColumn,1,smilesNameColumn,2" -i Sample.csv
          -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitDrawMolecules.py, RDKitRemoveDuplicateMolecules.py, RDKitSearchFunctionalGroups.py,
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
