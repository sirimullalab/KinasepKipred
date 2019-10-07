#!/bin/env python
#
# File: RDKitGenerateMolecularFrameworks.py
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
    from rdkit.Chem.Scaffolds import MurckoScaffold
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
    GenerateMolecularFrameworks()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def GenerateMolecularFrameworks():
    """Generate Bemis Murcko molecular framworks."""
    
    Infile = OptionsInfo["Infile"]
    Outfile = OptionsInfo["Outfile"]

    UseChirality = OptionsInfo["UseChirality"]

    RemoveDuplicateFrameworks = OptionsInfo["RemoveDuplicateFrameworks"]
    UseGraphFrameworks = OptionsInfo["UseGraphFrameworks"]
    
    SortFrameworks = OptionsInfo["SortFrameworks"]
    if SortFrameworks:
        FrameworkMolIDs = []
        FrameworkMolIDToMolMap = {}
        FrameworkMolIDToAtomCountMap = {}
        
        DuplicateFrameworkMolIDs = []
        DuplicateFrameworkMolIDToMolMap = {}
        DuplicateFrameworkMolIDToAtomCountMap = {}
        
    DuplicatesOutfile = ""
    if RemoveDuplicateFrameworks:
        DuplicatesOutfile = OptionsInfo["DuplicatesOutfile"]

    # Setup a molecule reader...
    MiscUtil.PrintInfo("\nProcessing file %s..." % Infile)
    Mols  = RDKitUtil.ReadMolecules(Infile, **OptionsInfo["InfileParams"])
    
    # Set up a molecular framework  writer...
    Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
    
    # Set up a duplicate molecular framework writer...    
    if RemoveDuplicateFrameworks:
        DuplicatesWriter = RDKitUtil.MoleculesWriter(DuplicatesOutfile, **OptionsInfo["OutfileParams"])
        if Writer is None:
            MiscUtil.PrintError("Failed to setup a writer for duplicates output fie %s " % DuplicatesOutfile)
        
    if RemoveDuplicateFrameworks:
        MiscUtil.PrintInfo("Generating files: %s and %s..." % (Outfile, DuplicatesOutfile))
    else:
        MiscUtil.PrintInfo("Generating file %s..." % Outfile)

    # Process molecules...
    MolCount = 0
    ValidMolCount = 0
    
    FrameworksCount = 0
    UniqueFrameworksCount = 0
    DuplicateFrameworksCount = 0
    
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

        if UseGraphFrameworks:
            FrameworksMol = MurckoScaffold.MakeScaffoldGeneric(Mol)
        else:
            FrameworksMol = MurckoScaffold.GetScaffoldForMol(Mol)

        if Compute2DCoords:
            AllChem.Compute2DCoords(FrameworksMol)
            
        if SortFrameworks:
            HeavyAtomCount = FrameworksMol.GetNumHeavyAtoms()

        FrameworksCount += 1
        
        if RemoveDuplicateFrameworks:
            CanonicalSMILES = Chem.MolToSmiles(FrameworksMol, isomericSmiles = UseChirality, canonical = True)
            if CanonicalSMILES in CanonicalSMILESMap:
                DuplicateFrameworksCount += 1
                if SortFrameworks:
                    # Track duplicate frameworks...
                    DuplicateFrameworkMolIDs.append(DuplicateFrameworksCount)
                    DuplicateFrameworkMolIDToMolMap[DuplicateFrameworksCount] = FrameworksMol
                    DuplicateFrameworkMolIDToAtomCountMap[DuplicateFrameworksCount] = HeavyAtomCount
                else:
                    # Write it out...
                    DuplicatesWriter.write(FrameworksMol)
            else:
                UniqueFrameworksCount += 1
                CanonicalSMILESMap[CanonicalSMILES] = CanonicalSMILES
                if SortFrameworks:
                    # Track unique frameworks...
                    FrameworkMolIDs.append(UniqueFrameworksCount)
                    FrameworkMolIDToMolMap[UniqueFrameworksCount] = FrameworksMol
                    FrameworkMolIDToAtomCountMap[UniqueFrameworksCount] = HeavyAtomCount
                else:
                    # Write it out...
                    Writer.write(FrameworksMol)
        elif SortFrameworks:
            # Track for sorting...
            FrameworkMolIDs.append(FrameworksCount)
            FrameworkMolIDToMolMap[FrameworksCount] = FrameworksMol
            FrameworkMolIDToAtomCountMap[FrameworksCount] = HeavyAtomCount
        else:
            # Write it out...
            Writer.write(FrameworksMol)
            
    if SortFrameworks:
        ReverseOrder = OptionsInfo["DescendingSortOrder"]
        SortAndWriteFrameworks(Writer, FrameworkMolIDs, FrameworkMolIDToMolMap, FrameworkMolIDToAtomCountMap, ReverseOrder)
        if RemoveDuplicateFrameworks:
            SortAndWriteFrameworks(DuplicatesWriter, DuplicateFrameworkMolIDs, DuplicateFrameworkMolIDToMolMap, DuplicateFrameworkMolIDToAtomCountMap, ReverseOrder)
    
    Writer.close()
    if RemoveDuplicateFrameworks:
        DuplicatesWriter.close()

    MiscUtil.PrintInfo("\nTotal number of molecular frameworks: %d" % FrameworksCount)
    if RemoveDuplicateFrameworks:
        MiscUtil.PrintInfo("Number of unique molecular frameworks: %d" % UniqueFrameworksCount)
        MiscUtil.PrintInfo("Number of duplicate molecular frameworks: %d" % DuplicateFrameworksCount)
        
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))

def SortAndWriteFrameworks(MolWriter, MolIDs, MolIDToMolMap, MolIDToAtomCountMap, ReverseOrder):
    """Sort frameworks and write them out."""
    SortedMolIDs = sorted(MolIDs, key = lambda MolID: MolIDToAtomCountMap[MolID], reverse = ReverseOrder)
    for MolID in SortedMolIDs:
        FrameworksMol = MolIDToMolMap[MolID]
        MolWriter.write(FrameworksMol)

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

    OptionsInfo["Mode"] = Options["--mode"]
    OptionsInfo["UseGraphFrameworks"] = False
    if re.match("^GraphFrameworks$", OptionsInfo["Mode"], re.I):
        OptionsInfo["UseGraphFrameworks"] = True
    
    OptionsInfo["RemoveDuplicates"] = Options["--removeDuplicates"]
    OptionsInfo["RemoveDuplicateFrameworks"] = False
    if re.match("^Yes$", OptionsInfo["RemoveDuplicates"], re.I):
        OptionsInfo["RemoveDuplicateFrameworks"] = True
        
    # Setup outfile for writing out duplicates...
    OptionsInfo["DuplicatesOutfile"] = ""
    if OptionsInfo["RemoveDuplicateFrameworks"]:
        FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Outfile"])
        OptionsInfo["DuplicatesOutfile"] = "%sDuplicates.%s" % (FileName, FileExt)

    OptionsInfo["Sort"] = Options["--sort"]
    OptionsInfo["SortFrameworks"] = False
    if re.match("^Yes$", OptionsInfo["Sort"], re.I):
        OptionsInfo["SortFrameworks"] = True
    
    OptionsInfo["SortOrder"] = Options["--sortOrder"]
    OptionsInfo["DescendingSortOrder"] = False
    if re.match("^Descending$", OptionsInfo["SortOrder"], re.I):
        OptionsInfo["DescendingSortOrder"] = True
    
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
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd smi")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    
    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "GraphFrameworks AtomicFrameworks")
    
    MiscUtil.ValidateOptionTextValue("-r, --removeDuplicates", Options["--removeDuplicates"], "yes no")
    MiscUtil.ValidateOptionTextValue("-s, --sort", Options["--sort"], "yes no")
    MiscUtil.ValidateOptionTextValue("--sortOrder", Options["--sortOrder"], "ascending descending")

    MiscUtil.ValidateOptionTextValue("--useChirality", Options["--useChirality"], "yes no")
    
# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitGenerateMolecularFrameworks.py - Generate Bemis Murcko molecular frameworks

Usage:
    RDKitGenerateMolecularFrameworks.py [--infileParams <Name,Value,...>]
                                         [--mode <GraphFrameworks or AtomicFrameworks> ]
                                         [ --outfileParams <Name,Value,...> ]  [--overwrite] [--removeDuplicates <yes or no>]
                                         [--sort <yes or no>] [--sortOrder <ascending or descending>]
                                         [--useChirality <yes or no>] [-w <dir>] -i <infile> -o <outfile>
    RDKitGenerateMolecularFrameworks.py -h | --help | -e | --examples

Description:
    Generate Bemis Murcko [ Ref 133 ] molecular frameworks for molecules. Two types of molecular
    frameworks can be generated: Graph or atomic frameworks. The graph molecular framework
    is a generic framework. The atom type, hybridization, and bond order is ignore during its
    generation. All atoms are set to carbon atoms and all bonds are single bonds. The atom type,
    hybridization, and bond order is preserved during generation of atomic molecular frameworks.

    The supported input file formats are: SD (.sdf, .sd), SMILES (.smi, .csv, .tsv, .txt)

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
    -m, --mode <GraphFrameworks or AtomicFrameworks>  [default: GraphFrameworks]
        Type of molecular frameworks to generate for molecules. Possible values: GraphFrameworks
         or AtomicFrameworks. The graph molecular framework is a generic framework. The atom type,
         hybridization, and bond order is ignore during its generation. All atoms are set to carbon atoms
         and all bonds are single bonds. The atom type, hybridization, and bond order is preserved
         during the generation of atomic molecular frameworks.
    -o, --outfile <outfile>
        Output file name.
    --outfileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for writing
        molecules to files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD: compute2DCoords,auto
            SMILES: smilesDelimiter,space,smilesTitleLine,yes
            
        Default value for compute2DCoords: yes for SMILES input file; no for all other
        file types.
    --overwrite
        Overwrite existing files.
    -r, --removeDuplicates <yes or no>  [default: no]
        Remove duplicate molecular frameworks. Possible values: yes or no. The duplicate
        molecular franworks are identified using canonical SMILES. The removed frameworks
        are written to a separate output file.
    -s, --sort <yes or no>  [default: no]
        Sort molecular frameworks by heavy atom count. Possible values: yes or no.
    --sortOrder <ascending or descending>  [default: ascending]
        Sorting order for molecular frameworks. Possible values: ascending or descending.
    -u, --useChirality <yes or no>  [default: yes]
        Use stereochemistry for generation of canonical SMILES strings to identify
        duplicate molecular frameworks.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To generate graph molecular framworks for molecules and write out a SMILES file,
    type:

        % RDKitGenerateMolecularFrameworks.py -i Sample.smi -o SampleOut.smi

    To generate graph molecular framworks, remove duplicate frameworks for molecules
    and write out SD files for unique and duplicate frameworks, type:

        % RDKitGenerateMolecularFrameworks.py -m GraphFrameworks -r yes
          -i Sample.sdf -o SampleOut.sdf

    To generate atomic molecular framworks, remove duplicate frameworks, sort
    framworks by heavy atom count in ascending order, write out SMILES files for
    unique and duplicate frameworks, type:

        % RDKitGenerateMolecularFrameworks.py -m AtomicFrameworks -r yes
          -s yes -i Sample.smi -o SampleOut.smi

    To generate graph molecular framworks for molecules in a CSV SMILES file,
    SMILES strings in column 1, name in olumn 2, emove duplicate frameworks,
    sort framworks by heavy atom count in decending order and write out a SD
    file, type:

        % RDKitGenerateMolecularFrameworks.py -m AtomicFrameworks
          --removeDuplicates yes -s yes --sortOrder descending --infileParams
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitDrawMolecules.py, RDKitSearchFunctionalGroups.py,
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
