#!/bin/env python
#
# File: RDKitCompareMoleculeShapes.py
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
    from rdkit.Chem import AllChem, rdMolAlign, rdShapeHelpers
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
    CompareMoleculeShapes()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def CompareMoleculeShapes():
    """Compare shape of molecules."""
    
    if not re.match("^(OneToOne|AllToAll|FirstToAll)$", OptionsInfo["Mode"], re.I):
        MiscUtil.PrintError("Shape comparison couldn't be performed: Specified mode, %s, is not supported" % OptionsInfo["Mode"])
        
    if not re.match("^(Open3A|CrippenOpen3A)$", OptionsInfo["Alignment"], re.I):
        MiscUtil.PrintError("Shape couldn't be performed: Specified alignment mode, %s, is not supported" % OptionsInfo["Alignment"])
        
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
    MiscUtil.PrintInfo("Generating file %s...\n" % Outfile)
    OutFH = open(Outfile, "w")
    if OutFH is None:
        MiscUtil.PrintError("Couldn't open output file: %s.\n" % (Outfile))

    if OptionsInfo["UseCrippenOpen3A"]:
        Line = "RefMolID%sProbeMolID%sCrippenOpen3AScore" % (OutDelim, OutDelim)
    else:
        Line = "RefMolID%sProbeMolID%sOpen3AScore" % (OutDelim, OutDelim)
        
    if OptionsInfo["CalcTanimotoDistance"]:
        Line = "%s%sTanimotoDistance" % (Line, OutDelim)
    if OptionsInfo["CalcProtrudeDistance"]:
        Line = "%s%sProtrudeDistance" % (Line, OutDelim)
    OutFH.write("%s\n" % Line)
        
    if re.match("^OneToOne$", OptionsInfo["Mode"], re.I):
        PerformOneToOneShapeComparison(ValidRefMols, ValidProbeMols, OutFH, OutDelim)
    elif re.match("^AllToAll$", OptionsInfo["Mode"], re.I):
        PerformAllToAllShapeComparison(ValidRefMols, ValidProbeMols, OutFH, OutDelim)
    elif re.match("^FirstToAll$", OptionsInfo["Mode"], re.I):
        PerformFirstToAllShapeComparison(ValidRefMols, ValidProbeMols, OutFH, OutDelim)
    else:
        MiscUtil.PrintError("Shape comaprison couldn't be performed: Specified mode, %s, is not supported" % OptionsInfo["Mode"])

    OutFH.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: Reference - %d; Probe - %d" % (RefMolCount, ProbeMolCount))
    MiscUtil.PrintInfo("Number of valid molecules: Reference - %d; Probe - %d" % (ValidRefMolCount, ValidProbeMolCount))
    MiscUtil.PrintInfo("Number of ignored molecules:  Reference - %d; Probe - %d" % ((RefMolCount - ValidRefMolCount), (ProbeMolCount - ValidProbeMolCount)))
    
def PerformOneToOneShapeComparison(ValidRefMols, ValidProbeMols, OutFH, OutDelim):
    """Perform pairwise shape comparison"""
    
    ValidRefMolCount = len(ValidRefMols)
    ValidProbeMolCount = len(ValidProbeMols)
    
    MolCount = ValidRefMolCount
    if ValidRefMolCount > ValidProbeMolCount:
        MolCount = ValidProbeMolCount
    
    if ValidRefMolCount != ValidProbeMolCount:
        MiscUtil.PrintWarning("Number of valid reference molecules, %d,  is not equal to number of valid probe molecules, %d .\n" % (ValidRefMolCount, ValidProbeMolCount))
        MiscUtil.PrintWarning("Pairwise shape comparison will be performed only for first %s molecules.\n" % (MolCount))

    # Process molecules...
    for MolIndex in range(0, MolCount):
        RefMol = ValidRefMols[MolIndex]
        ProbeMol = ValidProbeMols[MolIndex]

        RefMolName = RDKitUtil.GetMolName(RefMol, (MolIndex + 1))
        ProbeMolName = RDKitUtil.GetMolName(ProbeMol, (MolIndex + 1))

        PerformShapeComparisonAndWrieOutput(RefMol, ProbeMol, RefMolName, ProbeMolName, OutFH, OutDelim)
        
def PerformAllToAllShapeComparison(ValidRefMols, ValidProbeMols, OutFH, OutDelim):
    """Perform shape comparison between all pairs of molecules."""
    
    # Process molecules...
    RefMolCount = 0
    for RefMol in ValidRefMols:
        RefMolCount += 1
        RefMolName = RDKitUtil.GetMolName(RefMol, RefMolCount)

        ProbeMolCount = 0
        for ProbeMol in ValidProbeMols:
            ProbeMolCount += 1
            ProbeMolName = RDKitUtil.GetMolName(ProbeMol, ProbeMolCount)
            
            PerformShapeComparisonAndWrieOutput(RefMol, ProbeMol, RefMolName, ProbeMolName, OutFH, OutDelim)
        
def PerformFirstToAllShapeComparison(ValidRefMols, ValidProbeMols, OutFH, OutDelim):
    """Perform shape comparison between first reference molecues and all probe molecules. """
    
    # Process molecules...
    RefMol = ValidRefMols[0]
    RefMolCount = 1
    RefMolName = RDKitUtil.GetMolName(RefMol, RefMolCount)

    ProbeMolCount = 0
    for ProbeMol in ValidProbeMols:
        ProbeMolCount += 1
        ProbeMolName = RDKitUtil.GetMolName(ProbeMol, ProbeMolCount)
        
        PerformShapeComparisonAndWrieOutput(RefMol, ProbeMol, RefMolName, ProbeMolName, OutFH, OutDelim)

def PerformShapeComparisonAndWrieOutput(RefMol, ProbeMol, RefMolName, ProbeMolName, OutFH, OutDelim):
    """Perform shape comparison and write to output file."""

    AlignmentScore = PerformShapeAlignment(RefMol, ProbeMol)
    AlignmentScore = "%.2f" % AlignmentScore
    
    LineWords = []
    LineWords.extend([RefMolName, ProbeMolName, AlignmentScore])
    
    if OptionsInfo["CalcTanimotoDistance"]:
        TanimotoDistance = CalculateTanimotoShapeDistance(RefMol, ProbeMol)
        LineWords.append(TanimotoDistance)
    
    if OptionsInfo["CalcProtrudeDistance"]:
        ProtrudeDistance = CalculateProtrudeShapeDistance(RefMol, ProbeMol)
        LineWords.append(ProtrudeDistance)
    
    Line = OutDelim.join(LineWords)
    OutFH.write("%s\n" % Line)
    
def PerformShapeAlignment(RefMol, ProbeMol):
    """Perform shape alignment and return alignment score."""
    
    if OptionsInfo["UseCrippenOpen3A"]:
        CrippenO3A = rdMolAlign.GetCrippenO3A(ProbeMol, RefMol)
        Score = CrippenO3A.Align()
    else:
        O3A = rdMolAlign.GetO3A(ProbeMol, RefMol)
        Score = O3A.Align()

    return Score
        
def CalculateTanimotoShapeDistance(RefMol, ProbeMol):
    """Calculate Tanimoto shape for a pair of already aligned molecules and return it as a string"""

    Distance = rdShapeHelpers.ShapeTanimotoDist(ProbeMol, RefMol)
    Distance = "%.2f" % Distance

    return Distance

def CalculateProtrudeShapeDistance(RefMol, ProbeMol):
    """Calculate protrude shape for a pair of already aligned molecules and return it as a string"""

    Distance = rdShapeHelpers.ShapeProtrudeDist(ProbeMol, RefMol)
    Distance = "%.2f" % Distance

    return Distance

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["Alignment"] = Options["--alignment"]
    OptionsInfo["UseCrippenOpen3A"] = False
    if re.match("^CrippenOpen3A$", OptionsInfo["Alignment"], re.I):
        OptionsInfo["UseCrippenOpen3A"] = True
    
    OptionsInfo["Distance"] = Options["--distance"]
    
    OptionsInfo["CalcTanimotoDistance"] = False
    OptionsInfo["CalcProtrudeDistance"] = False
    if re.match("^Tanimoto$", OptionsInfo["Distance"], re.I):
        OptionsInfo["CalcTanimotoDistance"] = True
    elif re.match("^Protrude$", OptionsInfo["Distance"], re.I):
        OptionsInfo["CalcProtrudeDistance"] = True
    else:
        OptionsInfo["CalcTanimotoDistance"] = True
        OptionsInfo["CalcProtrudeDistance"] = True
    
    OptionsInfo["MaxIters"] = int(Options["--maxIters"])
    OptionsInfo["Mode"] = Options["--mode"]
    
    OptionsInfo["RefFile"] = Options["--reffile"]
    OptionsInfo["ProbeFile"] = Options["--probefile"]
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    
    # No need for any RDKit specific --outfileParams....
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"])
    
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
    
    MiscUtil.ValidateOptionTextValue("Alignment", Options["--alignment"], "Open3A CrippenOpen3A")
    MiscUtil.ValidateOptionTextValue("Distance", Options["--distance"], "Tanimoto Protrude Both")
    
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
RDKitCompareMoleculeShapes.py - Compare shapes of molecules

Usage:
    RDKitCompareMoleculeShapes.py [--alignment <Open3A, CrippenOpen3A>]
                                  [--distance <Tanimoto, Protrude, Both>]  [--infileParams <Name,Value,...>]
                                  [--maxIters <number>] [--mode <OneToOne, AllToAll, FirstToAll>]
                                  [--overwrite] [-w <dir>] -r <reffile> -p <probefile> -o <outfile> 
    RDKitCompareMoleculeShapes.py -h | --help | -e | --examples

Description:
    Compare shapes of molecules between a set of molecules in reference and probe
    input files. The molecules are aligned using either Open 3DAlign or Crippen Open
    3DAlign before calculating shape Tanimoto and protrude distances.

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd)

    The supported output file formats are: CSV (.csv), TSV (.tsv, .txt)

Options:
    -a, --alignment <Open3A, CrippenOpen3A>  [default: Open3A]
        Alignment methodology to use for aligning molecules before calculating Tanimoto and
        protrude shape distances. Possible values: Open3A or CrippenOpen3A. Open 3DAlign
        (Open3A) [ Ref 132 ] overlays molecules based on MMFF atom types and charges.
        Crippen Open 3DAlign (CrippenOpen3A) uses Crippen logP contributions to overlay
        molecules.
    -d, --distance <Tanimoto, Protrude, Both>  [default: Both]
        Shape comparison distance to calculate for comparing shapes of molecules. Possible
        values: Tanimoto, Protrude, or Both. Shape Tanimoto distance takes the volume
        overlay into account during the calculation of distance. Shape protrude distance,
        however, focuses on the volume mismatch.
    --infileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for reading
        molecules from files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD, MOL: removeHydrogens,yes,sanitize,yes,strictParsing,yes
            
    --maxIters <number>  [default: 50]
        Maximum number of iterations to perform for each molecule pair during alignment.
    -m, --mode <OneToOne, AllToAll, FirstToAll>  [default: OneToOne]
        Specify how molecules are handled in reference and probe input files during
        comparison of shapes between reference and probe molecules.  Possible values:
        OneToOne, AllToAll and AllToFirst. For OneToOne mode, the molecule shapes are
        calculated for each pair of molecules in the reference and probe file and the shape
        distances are written to the output file. For AllToAll mode, the shape distances are
        calculated for each reference molecule against all probe molecules. For FirstToAll mode,
        however, the shape distances are only calculated for the first reference molecule
        against all probe molecules.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -p, --probefile <probefile>
        Probe input file name.
    -r, --reffile <reffile>
        Reference input file name.
    -o, --outfile <outfile>
        Output file name for writing out shape distances. Supported text file extensions: csv or tsv.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To perform shape alignment using Open3A methodology between pair of molecules in
    reference and probe molecules in 3D SD files, calculate both Tanimoto and protrude
    distances, and write out a CSV file containing calculated distance values along with
    appropriate molecule IDs, type:

        % RDKitCompareMoleculeShapes.py  -r Sample3DRef.sdf -p Sample3DProb.sdf
          -o SampleOut.csv

    To perform shape alignment using Crippen Open3A methodology between all molecules in
    reference and probe molecules in 3D SD files, calculate only Tanimoto distance, and write
    out a TSV file containing calculated distance value along with appropriate molecule IDs, type:

        % RDKitCompareMoleculeShapes.py  -m AllToAll -a CrippenOpen3A -d Tanimoto
          -r Sample3DRef.sdf -p Sample3DProb.sdf -o SampleOut.csv

    To perform shape alignment using Open3A methodology between first reference molecule
    against all probe molecules in 3D SD files without removing hydrogens , calculate both
    Tanimoto and protrude distances, and write out a CSV file containing calculated distance values along with
    appropriate molecule IDs, type:

        % RDKitCompareMoleculeShapes.py -m FirstToAll -a Open3A -d Both 
          --infileParams "removeHydrogens,no" -r Sample3DRef.sdf
          -p Sample3DProb.sdf -o SampleOut.csv

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitCalculateRMSD.py, RDKitCalculateMolecularDescriptors.py, RDKitConvertFileFormat.py,
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
