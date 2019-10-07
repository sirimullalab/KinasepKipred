#!/bin/env python
#
# File: RDKitPerformMinimization.py
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
    from rdkit.Chem import Descriptors
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
    PerformMinimization()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformMinimization():
    """Perform minimization."""
    
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

    # Process molecules...
    MolCount = 0
    ValidMolCount = 0
    MinimizationFailedCount = 0

    if OptionsInfo["SkipConformerGeneration"]:
        MiscUtil.PrintInfo("Performing minimization without generation of conformers...\n")
    else:
        MiscUtil.PrintInfo("Performing minimization with generation of conformers...\n")
        
    for Mol in Mols:
        MolCount += 1
        
        if Mol is None:
            continue
        
        if RDKitUtil.IsMolEmpty(Mol):
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
            continue
        
        ValidMolCount += 1

        if OptionsInfo["SkipConformerGeneration"]:
            Status = MinimizeMolecule(Mol, MolCount, Writer)
        else:
            Status = GenerateAndMinimizeConformers(Mol, MolCount, Writer)
        
        if not Status:
            MinimizationFailedCount += 1
            
    if Writer is not None:
        Writer.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of molecules failed during conformation generation or minimization: %d" % MinimizationFailedCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount + MinimizationFailedCount))

def MinimizeMolecule(Mol, MolCount, Writer):
    "Minimize moleculer and write it out"
    
    if  OptionsInfo["AddHydrogens"]:
        Mol = Chem.AddHs(Mol)

    Status = 0
    try:
        if OptionsInfo["UseUFF"]:
            Status = AllChem.UFFOptimizeMolecule(Mol, maxIters = OptionsInfo["MaxIters"])
        elif OptionsInfo["UseMMFF"]:
            Status = AllChem.MMFFOptimizeMolecule(Mol,  maxIters = OptionsInfo["MaxIters"])
        else:
            MiscUtil.PrintError("Minimization couldn't be performed: Specified forcefield, %s, is not supported" % OptionsInfo["ForceField"])
    except RuntimeError as ErrMsg:
        MolName = RDKitUtil.GetMolName(Mol, MolCount)
        MiscUtil.PrintWarning("Minimization couldn't be performed for molecule %s:\n%s\n" % (MolName, ErrMsg))
        return False
    
    if Status != 0:
        MolName = RDKitUtil.GetMolName(Mol, MolCount)
        MiscUtil.PrintWarning("Minimization failed to converge for molecule %s in %d steps. Try using higher value for \"--maxIters\" option...\n" % (MolName, OptionsInfo["MaxIters"]))
        
    if  OptionsInfo["RemoveHydrogens"]:
        Mol = Chem.RemoveHs(Mol)
    
    Writer.write(Mol)
    
    return True

def GenerateAndMinimizeConformers(Mol, MolCount, Writer):
    "Generate and mininize conformers for a molecule and write out the lowest energy conformer."
    
    ConfIDs = EmbedMolecule(Mol)
    if not len(ConfIDs):
        MolName = RDKitUtil.GetMolName(Mol, MolCount)
        MiscUtil.PrintWarning("Minimization couldn't be performed for molecule %s: Embedding failed...\n" % MolName)
        return False

    CalcEnergyMap = {}
    for ConfID in ConfIDs:
        try:
            if OptionsInfo["UseUFF"]:
                Status = AllChem.UFFOptimizeMolecule(Mol, confId = ConfID, maxIters = OptionsInfo["MaxIters"])
            elif OptionsInfo["UseMMFF"]:
                Status = AllChem.MMFFOptimizeMolecule(Mol, confId = ConfID, maxIters = OptionsInfo["MaxIters"])
            else:
                MiscUtil.PrintError("Minimization couldn't be performed: Specified forcefield, %s, is not supported" % OptionsInfo["ForceField"])
        except RuntimeError as ErrMsg:
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            MiscUtil.PrintWarning("Minimization couldn't be performed for molecule %s:\n%s\n" % (MolName, ErrMsg))
            return False
        
        EnergyStatus, Energy = GetConformerEnergy(Mol, ConfID)
        if not EnergyStatus:
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            MiscUtil.PrintWarning("Failed to retrieve calculated energy for conformation number %d of molecule %s. Try again after removing any salts or cleaing up the molecule...\n" % (ConfID, MolName))
            return False
        
        if Status != 0:
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            MiscUtil.PrintWarning("Minimization failed to converge for conformation number %d of molecule %s in %d steps. Try using higher value for \"--maxIters\" option...\n" % (ConfID, MolName, OptionsInfo["MaxIters"]))
            
        CalcEnergyMap[ConfID] = Energy

    SortedConfIDs = sorted(ConfIDs, key = lambda ConfID: CalcEnergyMap[ConfID])
    MinEnergyConfID = SortedConfIDs[0]
        
    if  OptionsInfo["RemoveHydrogens"]:
        Mol = Chem.RemoveHs(Mol)
    
    Writer.write(Mol, confId = MinEnergyConfID)
    
    return True
    
def GetConformerEnergy(Mol, ConfID):
    "Calculate conformer energy"

    Status = True
    Energy = 0.0
    
    if OptionsInfo["UseUFF"]:
        UFFMoleculeForcefield = AllChem.UFFGetMoleculeForceField(Mol, confId = ConfID)
        if UFFMoleculeForcefield is None:
            Status = False
        else:
            Energy = UFFMoleculeForcefield.CalcEnergy()
    elif OptionsInfo["UseMMFF"]:
        MMFFMoleculeProperties = AllChem.MMFFGetMoleculeProperties(Mol)
        MMFFMoleculeForcefield = AllChem.MMFFGetMoleculeForceField(Mol, MMFFMoleculeProperties, confId = ConfID)
        if MMFFMoleculeForcefield is None:
            Status = False
        else:
            Energy = MMFFMoleculeForcefield.CalcEnergy()
    else:
        MiscUtil.PrintError("Couldn't retrieve conformer energy: Specified forcefield, %s, is not supported" % OptionsInfo["ForceField"])
    
    return (Status, Energy)
    
def EmbedMolecule(Mol):
    "Embed conformations"
    
    ConfIDs = []
    
    MaxConfs = OptionsInfo["MaxConfs"]
    RandomSeed = OptionsInfo["RandomSeed"]
    EnforceChirality = OptionsInfo["EnforceChirality"]
    UseExpTorsionAnglePrefs = OptionsInfo["UseExpTorsionAnglePrefs"]
    UseBasicKnowledge = OptionsInfo["UseBasicKnowledge"]

    ConfIDs = AllChem.EmbedMultipleConfs(Mol, numConfs = MaxConfs, randomSeed = RandomSeed, enforceChirality = EnforceChirality, useExpTorsionAnglePrefs = UseExpTorsionAnglePrefs, useBasicKnowledge = UseBasicKnowledge)
    
    return ConfIDs

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]

    OptionsInfo["AddHydrogens"] = True
    if re.match("^no$", Options["--addHydrogens"], re.I):
        OptionsInfo["AddHydrogens"] = False

    if re.match("^ETDG$", Options["--conformerGenerator"], re.I):
        ConformerGenerator = "ETDG"
        UseExpTorsionAnglePrefs = True
        UseBasicKnowledge = False
        SkipConformerGeneration = False
    elif re.match("^KDG$", Options["--conformerGenerator"], re.I):
        ConformerGenerator = "KDG"
        UseExpTorsionAnglePrefs = False
        UseBasicKnowledge = True
        SkipConformerGeneration = False
    elif re.match("^ETKDG$", Options["--conformerGenerator"], re.I):
        ConformerGenerator = "ETKDG"
        UseExpTorsionAnglePrefs = True
        UseBasicKnowledge = True
        SkipConformerGeneration = False
    elif re.match("^SDG$", Options["--conformerGenerator"], re.I):
        ConformerGenerator = "SDG"
        UseExpTorsionAnglePrefs = False
        UseBasicKnowledge = False
        SkipConformerGeneration = False
    else:
        ConformerGenerator = "None"
        UseExpTorsionAnglePrefs = False
        UseBasicKnowledge = False
        SkipConformerGeneration = True
    
    OptionsInfo["SkipConformerGeneration"] = SkipConformerGeneration
    OptionsInfo["ConformerGenerator"] = ConformerGenerator
    OptionsInfo["UseExpTorsionAnglePrefs"] = UseExpTorsionAnglePrefs
    OptionsInfo["UseBasicKnowledge"] = UseBasicKnowledge

    if re.match("^UFF$", Options["--forceField"], re.I):
        ForceField = "UFF"
        UseUFF = True
        UseMMFF = False
    elif re.match("^MMFF$", Options["--forceField"], re.I):
        ForceField = "MMFF"
        UseUFF = False
        UseMMFF = True
    
    OptionsInfo["ForceField"] = ForceField
    OptionsInfo["UseMMFF"] = UseMMFF
    OptionsInfo["UseUFF"] = UseUFF
        
    OptionsInfo["EnforceChirality"] = True
    if re.match("^no$", Options["--enforceChirality"], re.I):
        OptionsInfo["EnforceChirality"] = False
    
    OptionsInfo["MaxIters"] = int(Options["--maxIters"])
    OptionsInfo["MaxConfs"] = int(Options["--maxConfs"])
    
    RandomSeed = -1
    if not re.match("^auto$", Options["--randomSeed"], re.I):
        RandomSeed = int(Options["--randomSeed"])
    OptionsInfo["RandomSeed"] = RandomSeed
    
    OptionsInfo["RemoveHydrogens"] = True
    if re.match("^no$", Options["--removeHydrogens"], re.I):
        OptionsInfo["RemoveHydrogens"] = False

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
    
    MiscUtil.ValidateOptionTextValue("-a, --addHydrogens", Options["--addHydrogens"], "yes no")
    MiscUtil.ValidateOptionTextValue("-c, --conformerGenerator", Options["--conformerGenerator"], "SDG ETDG KDG ETKDG None")
    MiscUtil.ValidateOptionTextValue("-f, --forceField", Options["--forceField"], "UFF MMFF")
    
    MiscUtil.ValidateOptionTextValue("--enforceChirality ", Options["--enforceChirality"], "yes no")
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi txt csv tsv")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
        
    MiscUtil.ValidateOptionIntegerValue("--maxConfs", Options["--maxConfs"], {">": 0})
    MiscUtil.ValidateOptionIntegerValue("--maxIters", Options["--maxIters"], {">": 0})
    
    if not re.match("^auto$", Options["--randomSeed"], re.I):
        MiscUtil.ValidateOptionIntegerValue("--randomSeed", Options["--randomSeed"], {})
        RandomSeed = int(Options["--randomSeed"])
    
    MiscUtil.ValidateOptionTextValue("-r, --removeHydrogens", Options["--removeHydrogens"], "yes no")

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitPerformMinimization.py - Perform structure minimization

Usage:
    RDKitPerformMinimization.py [--addHydrogens <yes or no>] [--conformerGenerator <SDG, ETDG, KDG, ETKDG, None> ]
                                [--forceField <UFF or MMFF>] [--enforceChirality <yes or no>] [--infileParams <Name,Value,...>]
                                [--maxConfs <number>] [--maxIters <number>] [ --outfileParams <Name,Value,...> ]
                                [--overwrite] [ --removeHydrogens <yes or no>] [--randomSeed <number>] [-w <dir>] -i <infile> -o <outfile> 
    RDKitPerformMinimization.py -h | --help | -e | --examples

Description:
    Generate 3D structures for molecules using combination of distance geometry
    and forcefield minimization or minimize existing 3D structures using a specified
    forcefield.

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd), SMILES (.smi)
    .csv, .tcsv .txt)

    The supported output file formats are: SD (.sdf, .sd)

Options:
    -a, --addHydrogens <yes or no>  [default: yes]
        Add hydrogens before minimization.
    -c, --conformerGenerator <SDG, ETDG, KDG, ETKDG, None>  [default: ETKDG]
        Conformation generation methodology for generating initial 3D coordinates.
        Possible values: Standard Distance Geometry, (SDG), Experimental Torsion-angle
        preference with Distance Geometry (ETDG), basic Knowledge-terms with Distance
        Geometry (KDG),  and Experimental Torsion-angle preference along with basic
        Knowledge-terms with Distance Geometry (ETKDG) [Ref 129] .
        
        The conformation generation step may be skipped by specifying 'None' value to
        perform only forcefield minimization of molecules with 3D structures in input
        file.  This doesn't work for molecules in SMILES file or molecules in SD/MOL files
        containing 2D structures.
    -f, --forceField <UFF or MMFF>  [default: MMFF]
        Forcefield method to use for energy minimization. Possible values: Universal Force
        Field (UFF) [ Ref 81 ] or Merck Molecular Mechanics Force Field (MMFF) [ Ref 83-87 ] .
    --enforceChirality <yes or no>  [default: Yes]
        Enforce chirality for defined chiral centers.
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
    --maxConfs <number>  [default: 250]
        Maximum number of conformations to generate for each molecule by conformation
        generation methodology for initial 3D coordinates. The conformations are minimized
        using the specified forcefield and the lowest energy conformation is written to the
        output file. This option is ignored during 'None' value of '-c --conformerGenerator'
        option.
    --maxIters <number>  [default: 500]
        Maximum number of iterations to perform for each molecule during forcefield
        minimization.
    -o, --outfile <outfile>
        Output file name.
    --outfileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for writing
        molecules to files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD: kekulize,no
            
    --overwrite
        Overwrite existing files.
    -r, --removeHydrogens <yes or no>  [default: Yes]
        Remove hydrogens after minimization.
    --randomSeed <number>  [default: auto]
        Seed for the random number generator for reproducing 3D coordinates.
        Default is to use a random seed.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To generate up to 250 conformations using ETKDG methodology followed by MMFF
    forcefield minimization for a maximum of 500 iterations for molecules in a SMILES file
    and write out a SD file containing minimum energy structurse corresponding to each
    molecule, type

        % RDKitPerformMinimization.py  -i Sample.smi -o SampleOut.sdf

    To generate up to 150 conformations using ETKDG methodology followed by MMFF
    forcefield minimization for a maximum of 250 iterations along with a specified random
    seed  for molecules in a SMILES file and write out a SD file containing minimum energy
    structures corresponding to each molecule, type

        % RDKitPerformMinimization.py  --maxConfs 150  --randomSeed 201780117 
          --maxIters 250  -i Sample.smi -o SampleOut.sdf

    To minimize structures in a 3D SD file using UFF forcefield for a maximum of 150
    iterations without generating any conformations and write out a SD file containing
    minimum energy structures corresponding to each molecule, type

        % RDKitPerformMinimization.py  -c None -f UFF --maxIters 150
          -i Sample3D.sdf -o SampleOut.sdf

    To generate up to 50 conformations using SDG methodology followed
    by UFF forcefield minimization for a maximum of 50 iterations for 
    molecules in a CSV SMILES file, SMILES strings in column 1, name in
    column 2, and write out a SD file, type:

        % RDKitPerformMinimization.py  --maxConfs 50  --maxIters 50 -c SDG
          -f UFF --infileParams "smilesDelimiter,comma,smilesTitleLine,yes,
          smilesColumn,1,smilesNameColumn,2"  -i SampleSMILES.csv
          -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitCalculateRMSD.py, RDKitCalculateMolecularDescriptors.py, RDKitCompareMoleculeShapes.py,
    RDKitConvertFileFormat.py, RDKitGenerateConformers.py

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
