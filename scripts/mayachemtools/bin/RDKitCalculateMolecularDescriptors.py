#!/bin/env python
#
# File: RDKitCalculateMolecularDescriptors.py
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
    from rdkit.Chem import rdMolDescriptors
    from rdkit.Chem import Descriptors
    from rdkit.Chem import Descriptors3D
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
    CalculateMolecularDescriptors()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def CalculateMolecularDescriptors():
    """Calculate molecular descriptors."""

    ProcessMolecularDescriptorsInfo()
    PerformCalculations()

def ProcessMolecularDescriptorsInfo():
    """Process descriptors information."""

    RetrieveMolecularDescriptorsInfo()
    ProcessSpecifiedDescriptorNames()

def PerformCalculations():
    """Calculate descriptors for a specified list of descriptors."""

    Infile = OptionsInfo["Infile"]
    Outfile = OptionsInfo["Outfile"]
    
    DescriptorNames = OptionsInfo["SpecifiedDescriptorNames"]
    DescriptorCount = len(DescriptorNames)
    
    TextOutFileMode = OptionsInfo["TextOutFileMode"]
    TextOutFileDelim = OptionsInfo["TextOutFileDelim"]
    SMILESOut = OptionsInfo["SMILESOut"]

    Precision = OptionsInfo["Precision"]
    
    MiscUtil.PrintInfo("Calculating %d molecular descriptor(s) for each molecule..." % DescriptorCount)
    
    # Setup a molecule reader...
    MiscUtil.PrintInfo("\nProcessing file %s..." % Infile)
    Mols = RDKitUtil.ReadMolecules(Infile, **OptionsInfo["InfileParams"])
    
    # Setup a writer...
    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    if TextOutFileMode:
        Writer = open(Outfile, "w")
    else:
        Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
        
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)

    MiscUtil.PrintInfo("Generating file %s..." % Outfile)

    # Wite out headers for a text file...
    if TextOutFileMode:
        LineWords = []
        if SMILESOut:
            LineWords.append("SMILES")
        LineWords.append("MolID")
        LineWords.extend(DescriptorNames)
        Line = TextOutFileDelim.join(LineWords)
        Writer.write("%s\n" % Line)
        
    # Process molecules...
    MolCount = 0
    ValidMolCount = 0
    MatchCount = 0

    for Mol in Mols:
        MolCount += 1
        if Mol is None:
            continue
        
        if RDKitUtil.IsMolEmpty(Mol):
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
            continue
        
        ValidMolCount += 1

        # Calculate descriptors...
        CalculatedValues = []
        for Index in range(0, DescriptorCount):
            Name = DescriptorNames[Index]
            ComputeFunction = DescriptorNamesMap["ComputeFunction"][Name]
            Value = FormatCalculatedValue(ComputeFunction(Mol), Precision)
            CalculatedValues.append(Value)

        # Write out calculated values...
        if TextOutFileMode:
            LineWords = []
            if SMILESOut:
                SMILES = Chem.MolToSmiles(Mol, isomericSmiles = True, canonical = True)
                LineWords.append(SMILES)
            
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            LineWords.append(MolName)
            LineWords.extend(CalculatedValues)
            Line = TextOutFileDelim.join(LineWords)
            Writer.write("%s\n" % Line)
        else:
            for Index in range(0, DescriptorCount):
                Name = DescriptorNames[Index]
                Value = CalculatedValues[Index]
                Mol.SetProp(Name, Value)
                
            if Compute2DCoords:
                AllChem.Compute2DCoords(Mol)
                
            Writer.write(Mol)
    
    if Writer is not None:
        Writer.close()
        
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))

def FormatCalculatedValue(Value, Precision):
    """Format calculated value of descriptor based on its type."""

    if (type(Value) is float) or (type(Value) is numpy.float64):
        FormattedValue = "%.*f" % (Precision, Value)
        if not re.search("[1-9]", FormattedValue):
            FormattedValue = "0.0"
    elif type(Value) is list:
        FormattedValue = "%s" % Value
        FormattedValue = re.sub('\[|\]|,', '', FormattedValue)
    else:
        FormattedValue = "%s" % Value

    return FormattedValue

def ProcessSpecifiedDescriptorNames():
    """Process and validate specified decriptor names."""

    OptionsInfo["SpecifiedDescriptorNames"] = []

    if not re.match("^(2D|3D|All|FragmentCountOnly|Specify)$", OptionsInfo["Mode"], re.I):
        MiscUtil.PrintError("Mode value, %s, using \"-m, --mode\" option is not a valid value." % OptionsInfo["Mode"])
    
    if re.match("^2D$", OptionsInfo["Mode"], re.I):
        OptionsInfo["SpecifiedDescriptorNames"] = DescriptorNamesMap["2D"]["Names"]
        if OptionsInfo["FragmentCount"]:
            OptionsInfo["SpecifiedDescriptorNames"].extend(DescriptorNamesMap["FragmentCount"]["Names"])
        return
    elif re.match("^3D$", OptionsInfo["Mode"], re.I):
        OptionsInfo["SpecifiedDescriptorNames"] = DescriptorNamesMap["3D"]["Names"]
        return
    elif re.match("^All$", OptionsInfo["Mode"], re.I):
        OptionsInfo["SpecifiedDescriptorNames"] = DescriptorNamesMap["2D"]["Names"]
        if OptionsInfo["FragmentCount"]:
            OptionsInfo["SpecifiedDescriptorNames"].extend(DescriptorNamesMap["FragmentCount"]["Names"])
        OptionsInfo["SpecifiedDescriptorNames"].extend(DescriptorNamesMap["3D"]["Names"])
        return
    elif re.match("^FragmentCountOnly$", OptionsInfo["Mode"], re.I):
        OptionsInfo["SpecifiedDescriptorNames"] = DescriptorNamesMap["FragmentCount"]["Names"]
        return

    # Set up a canonical descriptor names map for checking specified names...
    CanonicalNameMap = {}
    for Name in  DescriptorNamesMap["ComputeFunction"]:
        CanonicalNameMap[Name.lower()] = Name
    
    # Parse and validate specified names...
    DescriptorNames = re.sub(" ", "", OptionsInfo["DescriptorNames"])
    if not DescriptorNames:
        MiscUtil.PrintError("No descriptor names specified for \"-d, --descriptorNames\" option")

    SMILESInfile = MiscUtil.CheckFileExt(Options["--infile"], "smi")
    Canonical3DNameMap = {}
    if SMILESInfile:
        for Name in DescriptorNamesMap["3D"]["Names"]:
            Canonical3DNameMap[Name.lower()] = Name
            
    SpecifiedDescriptorNames = []
    for Name in DescriptorNames.split(","):
        CanonicalName = Name.lower()
        if CanonicalName in CanonicalNameMap:
            SpecifiedDescriptorNames.append(CanonicalNameMap[CanonicalName])
        else:
            MiscUtil.PrintError("The descriptor name, %s, specified using \"-d, --descriptorNames\" option is not a valid name." % (Name))
        if SMILESInfile:
            if CanonicalName in Canonical3DNameMap:
                MiscUtil.PrintError("The 3D descriptor name, %s, specified using \"-d, --descriptorNames\" option is not a valid for SMILES input file." % (Name))
                
    if not len(SpecifiedDescriptorNames):
        MiscUtil.PrintError("No valid descriptor name specified for \"-d, --descriptorNames\" option")
    
    OptionsInfo["SpecifiedDescriptorNames"] = SpecifiedDescriptorNames

def RetrieveMolecularDescriptorsInfo():
    """Retrieve descriptors information."""

    MiscUtil.PrintInfo("\nRetrieving information for avalible molecular descriptors...")
    
    # Initialze data for 2D, FragmentCount and 3D descriptors...
    DescriptorNamesMap["Types"] = ["2D", "FragmentCount", "3D"]
    DescriptorNamesMap["ComputeFunction"] = {}

    Autocorr2DExclude = OptionsInfo["Autocorr2DExclude"]
    
    for Type in DescriptorNamesMap["Types"]:
        DescriptorNamesMap[Type] = {}
        DescriptorNamesMap[Type]["Names"] = []
    
    # Setup data for 2D and FragmentCount...
    DescriptorsInfo = Descriptors.descList
    for DescriptorInfo in DescriptorsInfo:
        Name = DescriptorInfo[0]
        ComputeFunction = DescriptorInfo[1]

        Type = "2D"
        if re.match("^fr_", Name, re.I):
            Type = "FragmentCount"
        elif re.match("^Autocorr2D$", Name, re.I):
            if Autocorr2DExclude:
                continue
        
        if Name in DescriptorNamesMap["ComputeFunction"]:
            MiscUtil.PrintWarning("Ignoring duplicate descriptor name: %s..." % Name)
        else:
            DescriptorNamesMap[Type]["Names"].append(Name)
            DescriptorNamesMap["ComputeFunction"][Name] = ComputeFunction

    # Add new 2D decriptor name to the list...
    Type = "2D"
    if not Autocorr2DExclude:
        try:
            Name = "Autocorr2D"
            ComputeFunction =  rdMolDescriptors.CalcAUTOCORR2D
            if not Name in DescriptorNamesMap["ComputeFunction"]:
                DescriptorNamesMap[Type]["Names"].append(Name)
                DescriptorNamesMap["ComputeFunction"][Name] = ComputeFunction
        except (AttributeError):
            MiscUtil.PrintInfo("2D descriptor, %s, is not available in your current version of RDKit." % Name)
        
    # Set data for 3D descriptors...
    Type = "3D"
    NameToComputeFunctionMap = {'PMI1' : Descriptors3D.PMI1, 'PMI2' : Descriptors3D.PMI2, 'PMI3' : Descriptors3D.PMI3, 'NPR1' : Descriptors3D.NPR1, 'NPR2' : Descriptors3D.NPR2, 'RadiusOfGyration' : Descriptors3D.RadiusOfGyration, 'InertialShapeFactor' :  Descriptors3D.InertialShapeFactor, 'Eccentricity' : Descriptors3D.Eccentricity, 'Asphericity' : Descriptors3D.Asphericity, 'SpherocityIndex' : Descriptors3D.SpherocityIndex}
    
    for Name in NameToComputeFunctionMap:
        ComputeFunction = NameToComputeFunctionMap[Name]
        if Name in DescriptorNamesMap["ComputeFunction"]:
            MiscUtil.PrintWarning("Ignoring duplicate descriptor name: %s..." % Name)
        else:
            DescriptorNamesMap[Type]["Names"].append(Name)
            DescriptorNamesMap["ComputeFunction"][Name] = ComputeFunction

    # Check and add new 3D descriptors not directly available through Descriptors3D module...
    Type = "3D"
    AvailableName3DMap = {}
    for Name in ['Autocorr3D', 'RDF', 'MORSE', 'WHIM', 'GETAWAY']:
        ComputeFunction = None
        try:
            if re.match("^Autocorr3D$", Name, re.I):
                ComputeFunction =  rdMolDescriptors.CalcAUTOCORR3D
            elif re.match("^RDF$", Name, re.I):
                ComputeFunction =  rdMolDescriptors.CalcRDF
            elif re.match("^MORSE$", Name, re.I):
                ComputeFunction =  rdMolDescriptors.CalcMORSE
            elif re.match("^WHIM$", Name, re.I):
                ComputeFunction =  rdMolDescriptors.CalcWHIM
            elif re.match("^GETAWAY$", Name, re.I):
                ComputeFunction =  rdMolDescriptors.CalcGETAWAY
            else:
                ComputeFunction = None
        except (AttributeError):
            MiscUtil.PrintWarning("3D descriptor, %s, is not available in your current version of RDKit" % Name)
        
        if ComputeFunction is not None:
            AvailableName3DMap[Name] = ComputeFunction

    for Name in AvailableName3DMap:
        ComputeFunction = AvailableName3DMap[Name]
        if not Name in DescriptorNamesMap["ComputeFunction"]:
            DescriptorNamesMap[Type]["Names"].append(Name)
            DescriptorNamesMap["ComputeFunction"][Name] = ComputeFunction

    Count = 0
    TypesCount = []
    for Type in DescriptorNamesMap["Types"]:
        TypeCount = len(DescriptorNamesMap[Type]["Names"])
        TypesCount.append(TypeCount)
        Count += TypeCount

    if not Count:
        MiscUtil.PrintError("Failed to retrieve any  molecular descriptors...")

    # Sort descriptor names...
    for Type in DescriptorNamesMap["Types"]:
        DescriptorNamesMap[Type]["Names"] = sorted(DescriptorNamesMap[Type]["Names"])
    
    MiscUtil.PrintInfo("\nTotal number of availble molecular descriptors: %d" % Count)
    for Index in range(0, len(DescriptorNamesMap["Types"])):
        Type = DescriptorNamesMap["Types"][Index]
        TypeCount = TypesCount[Index]
        MiscUtil.PrintInfo("Number of %s molecular descriptors: %d" % (Type, TypeCount))
    
    MiscUtil.PrintInfo("")
    
def ListMolecularDescriptorsInfo():
    """List descriptors information."""

    MiscUtil.PrintInfo("\nListing information for avalible molecular descriptors...")

    Delimiter = ", "
    for Type in DescriptorNamesMap["Types"]:
        Names = DescriptorNamesMap[Type]["Names"]
        MiscUtil.PrintInfo("\n%s descriptors: %s" % (Type, Delimiter.join(Names)))

    MiscUtil.PrintInfo("")

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()

    OptionsInfo["Autocorr2DExclude"] = True
    if not re.match("^Yes$", Options["--autocorr2DExclude"], re.I):
        OptionsInfo["Autocorr2DExclude"] = False
    
    OptionsInfo["FragmentCount"] = True
    if not re.match("^Yes$", Options["--fragmentCount"], re.I):
        OptionsInfo["FragmentCount"] = False
    
    OptionsInfo["DescriptorNames"] = Options["--descriptorNames"]
    OptionsInfo["Mode"] = Options["--mode"]
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"], Options["--infile"], Options["--outfile"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    
    TextOutFileMode = False
    TextOutFileDelim = ""
    if MiscUtil.CheckFileExt(Options["--outfile"], "csv"):
        TextOutFileMode = True
        TextOutFileDelim = ","
    elif MiscUtil.CheckFileExt(Options["--outfile"], "tsv txt"):
        TextOutFileMode = True
        TextOutFileDelim = "\t"
    OptionsInfo["TextOutFileMode"] = TextOutFileMode
    OptionsInfo["TextOutFileDelim"] = TextOutFileDelim
    
    OptionsInfo["Precision"] = int(Options["--precision"])
    
    OptionsInfo["SMILESOut"] = False
    if re.match("^Yes$", Options["--smilesOut"], re.I):
        OptionsInfo["SMILESOut"] = True

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
        
    # Handle listing of descriptor information...
    if  Options["--list"]:
        ProcessListMolecularDescriptorsOption()
        sys.exit(0)

def ProcessListMolecularDescriptorsOption():
    """Process list descriptors information."""

    RetrieveMolecularDescriptorsInfo()
    ListMolecularDescriptorsInfo()

def ValidateOptions():
    """Validate option values"""

    MiscUtil.ValidateOptionTextValue("-a, --autocorr2DExclude", Options["--autocorr2DExclude"], "yes no")
    MiscUtil.ValidateOptionTextValue("-f, --fragmentCount", Options["--fragmentCount"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "2D 3D All FragmentCountOnly Specify")
    
    if re.match("^Specify$", Options["--mode"], re.I):
        if re.match("^none$", Options["--descriptorNames"], re.I):
            MiscUtil.PrintError("The name(s) of molecular descriptors must be specified using \"-d, --descriptorNames\" option during \"Specify\" value of \"-m, --mode\" option.")
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi csv tsv txt")
    
    if re.match("^3D|All$", Options["--mode"], re.I):
        if MiscUtil.CheckFileExt(Options["--infile"], "smi"):
            MiscUtil.PrintError("The input SMILES file, %s, is not valid for  \"3D or All\" value of \"-m, --mode\" option.")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd csv tsv txt")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    
    MiscUtil.ValidateOptionIntegerValue("-p, --precision", Options["--precision"], {">": 0})
    MiscUtil.ValidateOptionTextValue("-s, --smilesOut", Options["--smilesOut"], "yes no")

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitCalculateMolecularDescriptors.py - Calculate 2D/3D molecular descriptors

Usage:
    RDKitCalculateMolecularDescriptors.py [--autocorr2DExclude <yes or no>] [--fragmentCount <yes or no>]
                                          [--descriptorNames <Name1,Name2,...>] [--infileParams <Name,Value,...>]
                                          [--mode <2D, 3D, All...>] [--outfileParams <Name,Value,...>]
                                          [--overwrite] [--precision <number>] [--smilesOut <yes or no>]
                                          [-w <dir>] -i <infile> -o <outfile> 
    RDKitCalculateMolecularDescriptors.py -l | --list
    RDKitCalculateMolecularDescriptors.py -h | --help | -e | --examples

Description:
    Calculate 2D/3D molecular descriptors for molecules and write them out to a SD or
    CSV/TSV text file.

    The complete list of currently available molecular descriptors may be obtained by
    using '-l, --list' option. The names of valid 2D, fragment count, and 3D molecular
    descriptors are shown below:

    2D descriptors: Autocorr2D, BalabanJ, BertzCT, Chi0, Chi1, Chi0n - Chi4n, Chi0v - Chi4v,
    EState_VSA1 - EState_VSA11, ExactMolWt, FpDensityMorgan1, FpDensityMorgan2, FpDensityMorgan3,
    FractionCSP3, HallKierAlpha, HeavyAtomCount, HeavyAtomMolWt, Ipc, Kappa1 - Kappa3,
    LabuteASA, MaxAbsEStateIndex, MaxAbsPartialCharge, MaxEStateIndex, MaxPartialCharge,
    MinAbsEStateIndex, MinAbsPartialCharge, MinEStateIndex, MinPartialCharge, MolLogP,
    MolMR, MolWt, NHOHCount, NOCount, NumAliphaticCarbocycles, NumAliphaticHeterocycles,
    NumAliphaticRings, NumAromaticCarbocycles, NumAromaticHeterocycles, NumAromaticRings,
    NumHAcceptors, NumHDonors, NumHeteroatoms, NumRadicalElectrons, NumRotatableBonds,
    NumSaturatedCarbocycles, NumSaturatedHeterocycles, NumSaturatedRings, NumValenceElectrons,
    PEOE_VSA1 - PEOE_VSA14,  RingCount, SMR_VSA1 - SMR_VSA10, SlogP_VSA1 - SlogP_VSA12,
    TPSA, VSA_EState1 - VSA_EState10, qed

    FragmentCount 2D descriptors: fr_Al_COO, fr_Al_OH, fr_Al_OH_noTert, fr_ArN, fr_Ar_COO,
    fr_Ar_N, fr_Ar_NH, fr_Ar_OH, fr_COO, fr_COO2, fr_C_O, fr_C_O_noCOO, fr_C_S, fr_HOCCN,
    fr_Imine, fr_NH0, fr_NH1, fr_NH2, fr_N_O, fr_Ndealkylation1, fr_Ndealkylation2, fr_Nhpyrrole,
    fr_SH, fr_aldehyde, fr_alkyl_carbamate, fr_alkyl_halide, fr_allylic_oxid, fr_amide, fr_amidine,
    fr_aniline, fr_aryl_methyl, fr_azide, fr_azo, fr_barbitur, fr_benzene, fr_benzodiazepine,
    fr_bicyclic, fr_diazo, fr_dihydropyridine, fr_epoxide, fr_ester, fr_ether, fr_furan, fr_guanido,
    fr_halogen, fr_hdrzine, fr_hdrzone, fr_imidazole, fr_imide, fr_isocyan, fr_isothiocyan, fr_ketone,
    fr_ketone_Topliss, fr_lactam, fr_lactone, fr_methoxy, fr_morpholine, fr_nitrile, fr_nitro,
    fr_nitro_arom, fr_nitro_arom_nonortho, fr_nitroso, fr_oxazole, fr_oxime, fr_para_hydroxylation,
    fr_phenol, fr_phenol_noOrthoHbond, fr_phos_acid, fr_phos_ester, fr_piperdine, fr_piperzine,
    fr_priamide, fr_prisulfonamd, fr_pyridine, fr_quatN, fr_sulfide, fr_sulfonamd, fr_sulfone,
    fr_term_acetylene, fr_tetrazole, fr_thiazole, fr_thiocyan, fr_thiophene, fr_unbrch_alkane, fr_urea

    3D descriptors: Asphericity, Autocorr3D, Eccentricity, GETAWAY, InertialShapeFactor, MORSE,
    NPR1, NPR2, PMI1, PMI2, PMI3, RDF, RadiusOfGyration, SpherocityIndex, WHIM

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd), SMILES (.smi,
    .txt, .csv, .tsv)

    The supported output file formats are: SD File (.sdf, .sd), CSV/TSV (.csv, .tsv, .txt)

Options:
    -a, --autocorr2DExclude <yes or no>  [default: yes]
        Exclude Autocorr2D descriptor from the calculation of 2D descriptors. 
    -f, --fragmentCount <yes or no>  [default: yes]
        Include 2D fragment count descriptors during the calculation. These descriptors are
        counted using SMARTS patterns specified in FragmentDescriptors.csv file distributed
        with RDKit. This option is only used during '2D' or 'All' value of '-m, --mode' option.
    -d, --descriptorNames <Name1,Name2,...>  [default: none]
        A comma delimited list of supported molecular descriptor names to calculate.
        This option is only used during 'Specify' value of '-m, --mode' option.
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
    -l, --list
        List molecular descriptors without performing any calculations.
    -m, --mode <2D, 3D, All, FragmentCountOnly, or Specify>  [default: 2D]
        Type of molecular descriptors to calculate. Possible values: 2D, 3D,
        All or Specify. The name of molecular descriptors must be specified using
        '-d, --descriptorNames' for 'Specify'. 2D descriptors also include 1D descriptors.
        The structure  of molecules must contain 3D coordinates for the  calculation
        of 3D descriptors.
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
        Floating point precision for writing the calculated descriptor values.
    -s, --smilesOut <yes or no>  [default: no]
        Write out SMILES string to CSV/TSV text output file.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To compute all available 2D descriptors except Autocorr2D descriptor and
    write out a CSV file, type:

        % RDKitCalculateMolecularDescriptors.py  -i Sample.smi -o SampleOut.csv

    To compute all available 2D descriptors including Autocorr2D descriptor and
    excluding fragment count descriptors, and write out a TSV file, type:

        % RDKitCalculateMolecularDescriptors.py  -m 2D -a no -f no
          -i Sample.smi -o SampleOut.tsv

    To compute all available 3D descriptors and write out a SD file, type:

        % RDKitCalculateMolecularDescriptors.py  -m 3D -i Sample3D.sdf
          -o Sample3DOut.sdf

    To compute only fragment count 2D descriptors and write out a SD
    file file, type:

        % RDKitCalculateMolecularDescriptors.py  -m FragmentCountOnly
          -i Sample.sdf -o SampleOut.sdf

    To compute all available 2D and 3D descriptors including fragment count and
    Autocorr2D and write out a CSV file, type:

        % RDKitCalculateMolecularDescriptors.py  -m All -a no -i Sample.sdf
          -o SampleOut.csv

    To compute a specific set of 2D and 3D descriptors and write out a
    write out a TSV file, type:

        % RDKitCalculateMolecularDescriptors.py  -m specify
          -d 'MolWt,MolLogP,NHOHCount, NOCount,RadiusOfGyration'
          -i Sample3D.sdf -o SampleOut.csv

    To compute all available 2D descriptors except Autocorr2D descriptor for 
    molecules in a CSV SMILES file, SMILES strings in column 1, name in
    column 2, and write out a SD file without calculation of 2D coordinates, type:

        % RDKitCalculateMolecularDescriptors.py --infileParams
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,no"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitCalculateRMSD.py, RDKitCompareMoleculeShapes.py, RDKitConvertFileFormat.py,
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
