#!/bin/env python
#
# File: RDKitCalculateRMSD.py
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
    from rdkit.Chem import AllChem, rdMolAlign
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
    CalculateRMSD()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def CalculateRMSD():
    """Calculate RMSD values."""
    
    if not re.match("^(OneToOne|AllToAll|FirstToAll)$", OptionsInfo["Mode"], re.I):
        MiscUtil.PrintError("RMSD couldn't be calculated: Specified mode, %s, is not supported" % OptionsInfo["Mode"])
        
    RefFile = OptionsInfo["RefFile"]
    ProbeFile = OptionsInfo["ProbeFile"]
    
    Outfile = OptionsInfo["Outfile"]
    OutDelim = OptionsInfo["OutDelim"]

    # Read reference and probe molecules...
    OptionsInfo["InfileParams"]["AllowEmptyMols"] = False
    
    MiscUtil.PrintInfo("\nProcessing file %s..." % (RefFile))
    ValidRefMols, RefMolCount, ValidRefMolCount  = RDKitUtil.ReadAndValidateMolecules(RefFile, **OptionsInfo["InfileParams"])
    
    MiscUtil.PrintInfo("Processing file %s..." % (ProbeFile))
    ValidProbeMols, ProbeMolCount, ValidProbeMolCount  = RDKitUtil.ReadAndValidateMolecules(ProbeFile, **OptionsInfo["InfileParams"])

    # Set up output file...
    MiscUtil.PrintInfo("\nGenerating file %s...\n" % Outfile)
    OutFH = open(Outfile, "w")
    if OutFH is None:
        MiscUtil.PrintError("Couldn't open output file: %s.\n" % (Outfile))
    
    Line = "RefMolID%sProbeMolID%sRMSD\n" % (OutDelim, OutDelim)
    OutFH.write(Line)
    
    if re.match("^OneToOne$", OptionsInfo["Mode"], re.I):
        CalculateOneToOneRMSDValues(ValidRefMols, ValidProbeMols, OutFH, OutDelim)
    elif re.match("^AllToAll$", OptionsInfo["Mode"], re.I):
        CalculateAllToAllRMSDValues(ValidRefMols, ValidProbeMols, OutFH, OutDelim)
    elif re.match("^FirstToAll$", OptionsInfo["Mode"], re.I):
        CalculateFirstToAllRMSDValues(ValidRefMols, ValidProbeMols, OutFH, OutDelim)
    else:
        MiscUtil.PrintError("RMSD couldn't be calculated: Specified mode, %s, is not supported" % OptionsInfo["Mode"])

    OutFH.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: Reference - %d; Probe - %d" % (RefMolCount, ProbeMolCount))
    MiscUtil.PrintInfo("Number of valid molecules: Reference - %d; Probe - %d" % (ValidRefMolCount, ValidProbeMolCount))
    MiscUtil.PrintInfo("Number of ignored molecules:  Reference - %d; Probe - %d" % ((RefMolCount - ValidRefMolCount), (ProbeMolCount - ValidProbeMolCount)))
    
def CalculateOneToOneRMSDValues(ValidRefMols, ValidProbeMols, OutFH, OutDelim):
    """Calculate pairwise RMSD values."""
    
    ValidRefMolCount = len(ValidRefMols)
    ValidProbeMolCount = len(ValidProbeMols)
    
    MolCount = ValidRefMolCount
    if ValidRefMolCount > ValidProbeMolCount:
        MolCount = ValidProbeMolCount
    
    if ValidRefMolCount != ValidProbeMolCount:
        MiscUtil.PrintWarning("Number of valid reference molecules, %d,  is not equal to number of valid probe molecules, %d .\n" % (ValidRefMolCount, ValidProbeMolCount))
        MiscUtil.PrintWarning("Pairwise RMSD will be calculated only for first %s molecules.\n" % (MolCount))

    # Process molecules...
    for MolIndex in range(0, MolCount):
        RefMol = ValidRefMols[MolIndex]
        ProbeMol = ValidProbeMols[MolIndex]

        RefMolName = RDKitUtil.GetMolName(RefMol, (MolIndex + 1))
        ProbeMolName = RDKitUtil.GetMolName(ProbeMol, (MolIndex + 1))

        RMSD = CalculateRMSDValue(RefMol, ProbeMol)
        
        Line = "%s%s%s%s%s\n" % (RefMolName, OutDelim, ProbeMolName, OutDelim, RMSD)
        OutFH.write(Line)
        
def CalculateAllToAllRMSDValues(ValidRefMols, ValidProbeMols, OutFH, OutDelim):
    """Calculate RMSD values between all pairs of molecules."""
    
    # Process molecules...
    RefMolCount = 0
    for RefMol in ValidRefMols:
        RefMolCount += 1
        RefMolName = RDKitUtil.GetMolName(RefMol, RefMolCount)

        ProbeMolCount = 0
        for ProbeMol in ValidProbeMols:
            ProbeMolCount += 1
            ProbeMolName = RDKitUtil.GetMolName(ProbeMol, ProbeMolCount)
            
            RMSD = CalculateRMSDValue(RefMol, ProbeMol)

            Line = "%s%s%s%s%s\n" % (RefMolName, OutDelim, ProbeMolName, OutDelim, RMSD)
            OutFH.write(Line)
        
def CalculateFirstToAllRMSDValues(ValidRefMols, ValidProbeMols, OutFH, OutDelim):
    """Calculate RMSD values between first reference molecues and all probe molecules. """
    
    # Process molecules...
    RefMol = ValidRefMols[0]
    RefMolCount = 1
    RefMolName = RDKitUtil.GetMolName(RefMol, RefMolCount)

    ProbeMolCount = 0
    for ProbeMol in ValidProbeMols:
        ProbeMolCount += 1
        ProbeMolName = RDKitUtil.GetMolName(ProbeMol, ProbeMolCount)
            
        RMSD = CalculateRMSDValue(RefMol, ProbeMol)

        Line = "%s%s%s%s%s\n" % (RefMolName, OutDelim, ProbeMolName, OutDelim, RMSD)
        OutFH.write(Line)
        
def CalculateRMSDValue(RefMol, ProbeMol):
    """Calculate RMSD value for a pair of molecules and return it as a string"""

    try:
        if OptionsInfo["UseBestRMSD"]:
            RMSD = AllChem.GetBestRMS(RefMol, ProbeMol)
        else:
            RMSD = rdMolAlign.AlignMol(ProbeMol, RefMol, maxIters = OptionsInfo["MaxIters"])
        RMSD = "%.2f" % RMSD
    except (RuntimeError, ValueError):
        RMSD = "None"
        
    return RMSD

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["CalcRMSD"] = Options["--calcRMSD"]
    OptionsInfo["UseBestRMSD"] = False
    if re.match("^BestRMSD$", OptionsInfo["CalcRMSD"], re.I):
        OptionsInfo["UseBestRMSD"] = True
    
    OptionsInfo["MaxIters"] = int(Options["--maxIters"])
    
    OptionsInfo["Mode"] = Options["--mode"]
    
    OptionsInfo["RefFile"] = Options["--reffile"]
    OptionsInfo["ProbeFile"] = Options["--probefile"]

    # No need for any RDKit specific --outfileParams....
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    
    OptionsInfo["OutDelim"] = " "
    if MiscUtil.CheckFileExt(OptionsInfo["Outfile"], "csv"):
        OptionsInfo["OutDelim"] = ","
    elif MiscUtil.CheckFileExt(OptionsInfo["Outfile"], "tsv txt"):
        OptionsInfo["OutDelim"] = "\t"
    else:
        MiscUtil.PrintError("The file name specified , %s, for option \"--outfile\" is not valid. Supported file formats: csv tsv txt\n" % (OptionsInfo["Outfile"]))

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
    
    MiscUtil.ValidateOptionTextValue("--calcRMSD", Options["--calcRMSD"], "RMSD BestRMSD")
    
    MiscUtil.ValidateOptionIntegerValue("--maxIters", Options["--maxIters"], {">": 0})
    MiscUtil.ValidateOptionTextValue("--mode", Options["--mode"], "OneToOne  AllToAll FirstToAll")
    
    MiscUtil.ValidateOptionFilePath("-r, --reffile", Options["--reffile"])
    MiscUtil.ValidateOptionFileExt("-r, --reffile", Options["--reffile"], "sdf sd mol")
    
    MiscUtil.ValidateOptionFilePath("-p, --probefile", Options["--probefile"])
    MiscUtil.ValidateOptionFileExt("-p, --probefile", Options["--probefile"], "sdf sd mol")
        
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "csv tsv txt")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-r, --reffile", Options["--reffile"], "-o, --outfile", Options["--outfile"])
    MiscUtil.ValidateOptionsDistinctFileNames("-p, --probefile", Options["--probefile"], "-o, --outfile", Options["--outfile"])
    
# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitCalculateRMSD.py - Calculate RMSD between molecules

Usage:
    RDKitCalculateRMSD.py [--calcRMSD <RMSD, BestRMSD>] [--infileParams <Name,Value,...>]
                          [--maxIters <number>] [--mode <OneToOne, AllToAll, FirstToAll>]
                          [--overwrite] [-w <dir>] -r <reffile> -p <probefile> -o <outfile> 
    RDKitCalculateRMSD.py -h | --help | -e | --examples

Description:
    Calculate Root Mean Square Distance (RMSD) between a set of similar molecules in
    reference and probe input files. The RDKit function fails to calculate RMSD values for
    dissimilar molecules. Consequently, a text string 'None' is written out as a RMSD value
    for dissimilar molecule pairs.

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd)

    The supported output file formats are:  CSV (.csv), TSV (.tsv, .txt)

Options:
    -c, --calcRMSD <RMSD, BestRMSD>  [default: RMSD]
        Methodology for calculating RMSD values. Possible values: RMSD, BestRMSD.
        During BestRMSMode mode, the RDKit 'function AllChem.GetBestRMS' is used to
        align and calculate RMSD. This function calculates optimal RMSD for aligning two
        molecules, taking symmetry into account. Otherwise, the RMSD value is calculated
        using 'AllChem.AlignMol function' without changing the atom order. A word to the
        wise from RDKit documentation: The AllChem.GetBestRMS function will attempt to
        align all permutations of matching atom orders in both molecules, for some molecules
        it will lead to 'combinatorial explosion'.
    --infileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for reading
        molecules from files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD, MOL: removeHydrogens,yes,sanitize,yes,strictParsing,yes
            
    --maxIters <number>  [default: 50]
        Maximum number of iterations to perform for each molecule pair during minimization
        of RMSD values. This option is ignored during BestRMSD mode.
    -m, --mode <OneToOne, AllToAll, FirstToAll>  [default: OneToOne]
        Specify how molecules are handled in reference and probe input files during
        calculation of RMSD between reference and probe molecules.  Possible values:
        OneToOne, AllToAll and AllToFirst. For OneToOne mode, the number of molecules
        in reference file must be equal to the number of molecules in probe file. The RMSD
        is calculated for each pair of molecules in the reference and probe file and written
        to the output file. For AllToAll mode, the RMSD is calculated for each reference
        molecule against all probe molecules. For FirstToAll mode, however, the RMSD
        is only calculated for the first reference molecule against all probe molecules.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -p, --probefile <probefile>
        Probe input file name.
    -r, --reffile <reffile>
        Reference input file name.
    -o, --outfile <outfile>
        Output file name for writing out RMSD values. Supported text file extensions: csv or tsv.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To calculate RMSD between pair of molecules in reference and probe molecules in
    3D SD files and write out a CSV file containing calculated RMSD values along with
    appropriate molecule IDs, type:

        % RDKitCalculateRMSD.py  -r Sample3DRef.sdf -p Sample3DProb.sdf
          -o SampleOut.csv

    To calculate RMSD between all molecules in reference and probe molecules in
    3D SD files and write out a CSV file containing calculated RMSD values along with
    appropriate molecule IDs, type:

        % RDKitCalculateRMSD.py  -m AllToAll -r Sample3DRef.sdf -p
          Sample3DProb.sdf -o SampleOut.csv

    To calculate best RMSD between first  molecule in reference all probe molecules
    in 3D SD files and write out a TSV file containing calculated RMSD values along with
    appropriate molecule IDs, type:

        % RDKitCalculateRMSD.py  -m FirstToAll --calcRMSD BestRMSD -r
          Sample3DRef.sdf -p Sample3DProb.sdf -o SampleOut.tsv

    To calculate RMSD between all molecules in reference and probe molecules in
    3D SD files without removing hydrogens and write out a TSV file containing
    calculated RMSD values along with appropriate molecule IDs, type:

        % RDKitCalculateRMSD.py  -m AllToAll --infileParams
          "removeHydrogens,no" -r Sample3DRef.sdf  -p Sample3DProb.sdf
          -o SampleOut.tsv

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitCalculateMolecularDescriptors.py, RDKitCompareMoleculeShapes.py, RDKitConvertFileFormat.py,
    RDKitGenerateConformers.py, RDKitPerformMinimization.py

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
