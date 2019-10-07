#!/bin/env python
#
# File: RDKitCalculatePartialCharges.py
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
    from rdkit.Chem import rdPartialCharges
except ImportError as ErrMsg:
    sys.stderr.write("\nFailed to import RDKit module/package: %s\n" % ErrMsg)
    sys.stderr.write("Check/update your RDKit environment and try again.\n\n")
    sys.exit(1)

# RDKit dependency imports...
import numpy

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

DescriptorNamesMap = {}

def main():
    """Start execution of the script"""
    
    MiscUtil.PrintInfo("\n%s (RDK v%s; %s): Starting...\n" % (ScriptName, rdBase.rdkitVersion, time.asctime()))
    
    (WallClockTime, ProcessorTime) = MiscUtil.GetWallClockAndProcessorTime()
    
    # Retrieve command line arguments and options...
    RetrieveOptions()
    
    # Process and validate command line arguments and options...
    ProcessOptions()
    
    # Perform actions required by the script...
    CalculatePartialCharges()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def CalculatePartialCharges():
    """Calculate partial atomic charges."""

    Infile = OptionsInfo["Infile"]
    Outfile = OptionsInfo["Outfile"]

    MiscUtil.PrintInfo("Calculating partial atomic charges...")
    
    # Setup a molecule reader...
    MiscUtil.PrintInfo("\nProcessing file %s..." % Infile)
    Mols = RDKitUtil.ReadMolecules(Infile, **OptionsInfo["InfileParams"])
    
    # Setup a writer...
    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
        
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)

    MiscUtil.PrintInfo("Generating file %s..." % Outfile)

    # Process molecules...
    MolCount, ValidMolCount, CalcFailedCount = [0] * 3 

    for Mol in Mols:
        MolCount += 1
        if Mol is None:
            continue
        
        if RDKitUtil.IsMolEmpty(Mol):
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
            continue
        
        ValidMolCount += 1
        MolWithHs = Chem.AddHs(Mol)

        # Retrieve charges...
        PartialCharges = CalculateMolPartialCharges(MolWithHs, MolCount)
        if not len(PartialCharges):
            CalcFailedCount += 1
            continue
            
        # Write out charges...
        WriteMolPartialCharges(Writer, MolWithHs, PartialCharges, Compute2DCoords)
        
    if Writer is not None:
        Writer.close()
        
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of molecules failed during calculation of partial charges: %d" % CalcFailedCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount + CalcFailedCount))

def CalculateMolPartialCharges(Mol, MolCount):
    """Calculate partial atomic charges for a molecule."""
    
    PartialCharges = []
    if OptionsInfo["MMFFChargesMode"]:
        if  AllChem.MMFFHasAllMoleculeParams(Mol):
            MMFFProp = AllChem.MMFFGetMoleculeProperties(Mol)
            PartialCharges = [MMFFProp.GetMMFFPartialCharge(AtomIndex) for AtomIndex in range(Mol.GetNumAtoms())]
        else:
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            MiscUtil.PrintWarning("Failed to calculate MMFF partial charges for molecule, %s: Missing forcefield parameters" % MolName)
            return PartialCharges
    else:
        rdPartialCharges.ComputeGasteigerCharges(Mol, nIter = OptionsInfo["NumIters"], throwOnParamFailure = OptionsInfo["AllowParamFailure"])
        PartialCharges = [Atom.GetProp("_GasteigerCharge") for Atom in Mol.GetAtoms()]
    
    # Format charges...
    PartialCharges = ["%.*f" % (OptionsInfo["Precision"], float(Value)) for Value in PartialCharges]
    
    return PartialCharges
    
def WriteMolPartialCharges(Writer, Mol, PartialCharges, Compute2DCoords):
    """Write out partial atomic charges for a molecule."""
    
    if OptionsInfo["AtomAliasesFormatMode"]:
        for Atom, PartialCharge in zip(Mol.GetAtoms(), PartialCharges):
            Atom.SetProp('molFileAlias', PartialCharge)
    else:
        ChargesValues = "\n".join(PartialCharges)
        Mol.SetProp(OptionsInfo["DataFieldLabel"], ChargesValues)
    
    if Compute2DCoords:
        AllChem.Compute2DCoords(Mol)
    
    Writer.write(Mol)
    
def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()

    AllowParamFailure = True
    if re.match("^No", Options["--allowParamFailure"], re.I):
        AllowParamFailure = False
    OptionsInfo["AllowParamFailure"] = AllowParamFailure
    
    AtomAliasesFormatMode = True
    if re.match("^DataField", Options["--chargesSDFormat"], re.I):
        AtomAliasesFormatMode = False
    OptionsInfo["AtomAliasesFormatMode"] = AtomAliasesFormatMode

    OptionsInfo["DataFieldLabel"] = Options["--dataFieldLabel"]
    
    MMFFChargesMode = False
    if re.match("^MMFF", Options["--mode"], re.I):
        MMFFChargesMode = True
    OptionsInfo["Mode"] = Options["--mode"]
    OptionsInfo["MMFFChargesMode"] = MMFFChargesMode
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"], Options["--infile"], Options["--outfile"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    
    OptionsInfo["NumIters"] = int(Options["--numIters"])
    OptionsInfo["Precision"] = int(Options["--precision"])

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

    MiscUtil.ValidateOptionTextValue("-a, --allowParamFailure", Options["--allowParamFailure"], "yes no")
    MiscUtil.ValidateOptionTextValue("-c, --chargesSDFormat", Options["--chargesSDFormat"], "AtomAliases DataField")
    
    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "Gasteiger MMFF")
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi csv tsv txt")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])

    MiscUtil.ValidateOptionIntegerValue("-n, --numIters", Options["--numIters"], {">": 0})
    MiscUtil.ValidateOptionIntegerValue("-p, --precision", Options["--precision"], {">": 0})
        
# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitCalculatePartialCharges.py - Calculate partial atomic charges

Usage:
    RDKitCalculatePartialCharges.py [--allowParamFailure <yes or no>]
                                          [--chargesSDFormat <AtomAliases or DataField>]  [--dataFieldLabel <text>]
                                          [--infileParams <Name,Value,...>] [--mode <Gasteiger or MMFF>]  [--numIters <number>]
                                          [--outfileParams <Name,Value,...>] [--precision <number>] [--overwrite]
                                          [-w <dir>] -i <infile> -o <outfile> 
    RDKitCalculatePartialCharges.py -h | --help | -e | --examples

Description:
    Calculate partial charges for atoms in molecules and write them out to a SD file.
    The hydrogens are automatically added to molecules before calculating partial
    charges.

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd), SMILES (.smi,
    .txt, .csv, .tsv)

    The supported output file format are: SD File (.sdf, .sd)

Options:
    -a, --allowParamFailure <yes or no>  [default: yes]
        Allow calculation of Gasteiger partial charges to proceed for molecules
        containing atoms with unknown parameters. The atoms with unknown
        parameters are removed from the calculations by setting their values to
        zero.
    -c, --chargesSDFormat <AtomAliases or DataField>  [default: AtomAliases]
        Format for writing out partial atomic charges to SD file. Possible values:
        AtomAliases or DataField.
        
        The charges are stored as atom property named 'molFileAlias' for
        'AtomAliases' format and may be retrieved using the RDKit function
        'GetProp' for atoms: Aotm.GetProp('molFileAliases').
        
        The charges are stored under a data field label speciefied using
        '-d, --dataFieldLabel' for 'DataField' format and may be retrieved using the
        RDKit function 'GetProp' for molecues.
    -d, --dataFieldLabel <text>  [default: PartialCharges]
        Data field label to use for storing charged in SD file during 'DataField' value
        of '-c, --chargesSDFormat'.
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
    -m, --mode <Gasteiger or MMFF>  [default: Gasteiger]
        Type of partial atomic charges to calculate. Possible values: Gasteiger
        [ Ref 138 ] or Merk Molecular Mechanics Fore Field (MMFF) [ Ref 83-87 ].
    -n, --numIters <number>  [default: 12]
        Number of iterations to perform during calculation of Gasteiger charges.
    -o, --outfile <outfile>
        Output file name.
    --outfileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for writing
        molecules to files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD: compute2DCoords,auto,kekulize,no
            
        Default value for compute2DCoords: yes for SMILES input file; no for all other
        file types.
    -p, --precision <number>  [default: 3]
        Floating point precision for writing the calculated partial atomic charges.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To calculate Gasteiger partial atomic charges for molecules in a SMILES
    file and write them out to a SD file as atom aliases, type:

        % RDKitCalculatePartialCharges.py  -i Sample.smi -o SampleOut.sdf

    To calculate MMFF forcefield partial atomic charges for molecules in a SD
    file and write them out to a SD file under 'PartialCharges' data field, type:

        % RDKitCalculatePartialCharges.py  -m MMFF -c DataField -i Sample.sdf
          -o SampleOut.sdf

    To calculate Gasteiger partial atomic charges for molecules in a SMILES
    file and write them out to a SD file under a data field named 'GasteigerCharges',
    type:

        % RDKitCalculatePartialCharges.py  -m Gasteiger -c DataField
          -d GasteigerCharges -p 4 -i Sample.smi -o SampleOut.sdf

    To calculate Gasteiger partial atomic charges for molecules in a CSV SMILES
    file, SMILES strings in column 1, name in column 2, and write out a SD file
    containing charges as atom aliases, type:

        % RDKitCalculatePartialCharges.py --infileParams
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitCalculateMolecularDescriptors.py, RDKitCalculateRMSD.py,
    RDKitCompareMoleculeShapes.py, RDKitConvertFileFormat.py,

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
