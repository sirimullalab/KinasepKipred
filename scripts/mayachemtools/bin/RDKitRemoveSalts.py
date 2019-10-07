#!/bin/env python
#
# File: RDKitRemoveSalts.py
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
    from rdkit.Chem.SaltRemover import SaltRemover
    from rdkit.Chem.SaltRemover import InputFormat
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
    RemoveSalts()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def RemoveSalts():
    """Identify and remove salts from molecules"""
    
    Infile = OptionsInfo["Infile"]
    Outfile = OptionsInfo["Outfile"]
    
    CountMode = OptionsInfo["CountMode"]
    
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
        
        MiscUtil.PrintInfo("Generating file %s..." % (Outfile))

    # Set up a salt remover...
    SaltsByComponentsMode = OptionsInfo["SaltsByComponentsMode"]
    Remover = None
    if not SaltsByComponentsMode:
        Remover = SaltRemover(defnFilename = OptionsInfo["SaltsFile"], defnData = OptionsInfo["SaltsSMARTS"], defnFormat = InputFormat.SMARTS)
    
    # Process molecules...
    MolCount = 0
    ValidMolCount = 0
    SaltsMolCount = 0
    
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
        
        UnsaltedMol, SaltyStatus = RemoveMolSalts(Mol, Remover, MolCount)
        
        if SaltyStatus:
            SaltsMolCount += 1

        if not CountMode:
            if Compute2DCoords:
                AllChem.Compute2DCoords(UnsaltedMol)
            Writer.write(UnsaltedMol)
    
    if Writer is not None:
        Writer.close()
    
    if DuplicatesWriter is not None:
        DuplicatesWriter.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))
    
    MiscUtil.PrintInfo("\nNumber of molecules coontaining salts: %d" % (SaltsMolCount))

def RemoveMolSalts(Mol, Remover, MolCount):
    """Remove salts from mol and return unsalted mol along with mol salty status."""

    UnsaltedMol = Mol
    SaltyStatus = False
    
    if Remover is not None:
        KeptMol, DeletedMols = Remover.StripMolWithDeleted(Mol, dontRemoveEverything = False)
        if len(DeletedMols) >= 1:
            SaltyStatus = True
        if RDKitUtil.IsMolEmpty(KeptMol):
            if len(DeletedMols) >= 1:
                # Take the larged fragment from DeletedMols
                UnsaltedMol = GetLargestMol(DeletedMols)
    else:
        # Use largest fragment as unsalted molecule...
        MolFrags = Chem.GetMolFrags(Mol, asMols = True)
        if len(MolFrags) > 1:
            # Keep the largest fragment as unsalted molecule...
            SaltyStatus = True
            UnsaltedMol = GetLargestMol(MolFrags)

    if SaltyStatus:
        Chem.SanitizeMol(UnsaltedMol)
        MolName = RDKitUtil.GetMolName(Mol, MolCount)
        if len(MolName):
            UnsaltedMol.SetProp("_Name", MolName)
    
    return (UnsaltedMol, SaltyStatus)

def GetLargestMol(Mols):
    """Get largest mol from list of mols"""

    LargestMol = None
    LargestMolSize = -1
    for Mol in Mols:
        Size = Mol.GetNumAtoms()
        if Size > LargestMolSize:
            LargestMol = Mol
            LargestMolSize = Size

    return LargestMol
    
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

    SaltsByComponentsMode = False
    SaltsBySMARTSFileMode = False
    SaltsBySMARTSMode = False
    if re.match("^ByComponent$", Options["--saltsMode"], re.I):
        SaltsByComponentsMode = True
    elif re.match("^BySMARTSFile$", Options["--saltsMode"], re.I):
        SaltsBySMARTSFileMode = False
    elif re.match("^BySMARTS$", Options["--saltsMode"], re.I):
        SaltsBySMARTSMode = True
    else:
        MiscUtil.PrintError("The salts mode specified, %s, using \"--saltsMode\" option is not valid." % Options["--saltsMode"])
    OptionsInfo["SaltsByComponentsMode"]  = SaltsByComponentsMode
    OptionsInfo["SaltsBySMARTSFileMode"]  = SaltsBySMARTSFileMode
    OptionsInfo["SaltsBySMARTSMode"]  = SaltsBySMARTSMode

    SaltsFile = None
    if re.match("^BySMARTSFile$", Options["--saltsMode"], re.I):
        if not re.match("^auto$", Options["--saltsFile"], re.I):
            SaltsFile = Options["--saltsFile"]
    OptionsInfo["SaltsFile"] = SaltsFile
    
    SaltsSMARTS = None
    if re.match("^BySMARTS$", Options["--saltsMode"], re.I):
        if not Options["--saltsSMARTS"]:
            MiscUtil.PrintError("No salts SMARTS pattern specified using \"--saltsSMARTS\" option during \"BySMARTS\" value of \"-s, --saltsMode\" option")
        SaltsSMARTS = Options["--saltsSMARTS"].strip(" ")
        if not len(SaltsSMARTS):
            MiscUtil.PrintError("Empty SMARTS pattern specified using \"--saltsSMARTS\" option during \"BySMARTS\" value of \"-s, --saltsMode\" option")
        if re.search(" ", SaltsSMARTS):
            SaltsSMARTS = re.sub('[ ]+', '\n', SaltsSMARTS)
        
    OptionsInfo["SaltsSMARTS"] = SaltsSMARTS
    
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
    
    MiscUtil.ValidateOptionTextValue("--saltsMode", Options["--saltsMode"], "ByComponent BySMARTSFile BySMARTS")
    
    if re.match("^BySMARTSFile$", Options["--saltsMode"], re.I):
        if not re.match("^auto$", Options["--saltsFile"], re.I):
            MiscUtil.ValidateOptionFilePath("--saltsFile", Options["--saltsFile"])

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitRemoveSalts.py - Remove salts

Usage:
    RDKitRemoveSalts.py  [--infileParams <Name,Value,...>] [--mode <remove or count>] [--outfileParams <Name,Value,...> ] [--overwrite]
                         [--saltsMode <ByComponent, BySMARTSFile, BySMARTS>] [--saltsFile <FileName or auto>] [--saltsSMARTS <SMARTS>]
                         [-w <dir>] [-o <outfile>]  -i <infile>
    RDKitRemoveSalts.py -h | --help | -e | --examples

Description:
    Remove salts from molecules or simply count the number of molecules containing
    salts. Salts are identified and removed based on either SMARTS strings or by selecting
    the largest disconnected components in molecules as non-salt portion of molecules.

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
        Specify whether to remove salts from molecules and write out molecules
        or or simply count the number of molecules containing salts.
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
    -s, --saltsMode <ByComponent, BySMARTSFile, BySMARTS>  [default: ByComponent]
        Specify whether to identify and remove salts based on SMARTS strings or
        by selecting the largest disconnected component as non-salt portion of a
        molecule. Possible values: ByComponent, BySMARTSFile or BySMARTS.
    --saltsFile <FileName or auto>  [default: auto]
        Specify a file name containing specification for SMARTS corresponding to salts or
        use default salts file, Salts.txt, available in RDKit data directory. This option is only
        used during 'BySMARTSFile' value of '-s, --saltsMode' option.
        
        RDKit data format: Smarts<tab>Name(optional)
        
        For example:
            
            [Cl,Br,I]
            [N](=O)(O)O
            [CH3]C(=O)O	  Acetic acid
            
    --saltsSMARTS <SMARTS text>
        Space delimited SMARTS specifications to use for salts identification instead
        their specifications in '--saltsFile'. This option is only used during 'BySMARTS'
        value of '-s, --saltsMode' option.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To remove salts from molecules in a SMILES file by keeping largest disconnected
    components as non-salt portion of molecules and write out a SMILES file, type:

        % RDKitRemoveSalts.py -i Sample.smi -o SampleOut.smi

    To count number of molecule containing salts from in a SD file, using largest
    components as non-salt portion of molecules, without generating any output
    file, type:

        % RDKitRemoveSalts.py -m count -i Sample.sdf

    To remove salts from molecules in a SMILES file using SMARTS strings in default
    Salts.txt distributed with RDKit to identify salts and write out a SMILES file, type:

        % RDKitRemoveSalts.py -m remove -s BySMARTSFile -i Sample.smi
          -o SampleOut.smi

    To remove salts from molecules in a SD file using SMARTS strings in a local
    CustomSalts.txt to identify salts and write out a SMILES file, type:

        % RDKitRemoveSalts.py -m remove -s BySMARTSFile --saltsFile
          CustomSalts.txt -i Sample.sdf -o SampleOut.smi

    To remove salts from molecules in a SD file using specified SMARTS to identify
    salts and write out a SD file, type:

        % RDKitRemoveSalts.py -m remove -s BySMARTS  --saltsSMARTS
          '[Cl,Br,I]  [N](=O)(O)O [N](=O)(O)O'
          -i Sample.sdf -o SampleOut.smi

    To remove salts form  molecules from a CSV SMILES file, SMILES strings in column 1,
    name in column 2, and generate output SD file, type:

        % RDKitRemoveSalts.py --infileParams 
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitRemoveDuplicateMolecules.py,
    RDKitSearchFunctionalGroups.py, RDKitSearchSMARTS.py

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
