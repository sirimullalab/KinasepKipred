#!/bin/env python
#
# File: RDKitEnumerateCompoundLibrary.py
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
    from rdkit.Chem import FunctionalGroups
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

RxnNamesMap = {}

def main():
    """Start execution of the script"""
    
    MiscUtil.PrintInfo("\n%s (RDK v%s; %s): Starting...\n" % (ScriptName, rdBase.rdkitVersion, time.asctime()))
    
    (WallClockTime, ProcessorTime) = MiscUtil.GetWallClockAndProcessorTime()
    
    # Retrieve command line arguments and options...
    RetrieveOptions()
    
    # Process and validate command line arguments and options...
    ProcessOptions()
    
    # Perform actions required by the script...
    PerformChemicalLibraryEnumeration()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformChemicalLibraryEnumeration():
    """Retrieve functional groups information and perform search."""

    ProcessReactionNamesInfo()
    PerformEnumeration()

def PerformEnumeration():
    """Enumerate virutal compound library."""

    ReactantFilesList = OptionsInfo["ReactantFilesList"]
    Outfile = OptionsInfo["Outfile"]

    RxnByNameMode = OptionsInfo["RxnByNameMode"]
    if RxnByNameMode:
        RxnSMIRKSPattern = OptionsInfo["RxnNameSMIRKS"]
    else:
        RxnSMIRKSPattern = OptionsInfo["SpecifiedSMIRKS"]

    # Set up a reaction and match number of reactants in rxn SMIRKS against number of
    # reactant files...
    Rxn = AllChem.ReactionFromSmarts(RxnSMIRKSPattern)
    RxnReactantsCount = Rxn.GetNumReactantTemplates()

    ReactantFilesList = OptionsInfo["ReactantFilesList"]
    ReactantFilesCount = len(ReactantFilesList)
    if  ReactantFilesCount != RxnReactantsCount:
        MiscUtil.PrintError("The number of specified reactant files, %d, must match number of rectants, %d, in reaction SMIRKS" % (ReactantFilesCount, RxnReactantsCount))
        
    # Retrieve reactant molecules...
    ReactantsMolsList = RetrieveReactantsMolecules()
    
    # Set up  a molecule writer...
    Writer = None
    Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)

    MiscUtil.PrintInfo("\nGenerating file %s..." % Outfile)

    # Set up reaction...
    ReturnReactants = False
    if OptionsInfo["UseReactantNames"]:
        ReturnReactants = True
    RxnProducts = AllChem.EnumerateLibraryFromReaction(Rxn, ReactantsMolsList, ReturnReactants)

    # Generate product molecules and write them out...
    
    Compute2DCoords = OptionsInfo["Compute2DCoords"]
    Sanitize = OptionsInfo["Sanitize"]
    
    ProdMolCount = 0
    ValidProdMolCount = 0
    
    if ReturnReactants:
        for Products, Reactants in list(RxnProducts):
            for ProdMol in Products:
                ProdMolCount += 1

                # Set product name...
                ReactantMolNames = [ReactantMol.GetProp("_Name") for ReactantMol in Reactants]
                Delimiter = "_"
                ProdMolName = Delimiter.join(ReactantMolNames) + "_Prod%d" % ProdMolCount
                ProdMol.SetProp("_Name", ProdMolName)

                Status = WriteProductMolecule(Writer, ProdMol, Sanitize, Compute2DCoords)
                if Status:
                    ValidProdMolCount += 1
    else:
        for Products in list(RxnProducts):
            for ProdMol in Products:
                ProdMolCount += 1

                # Set product name...
                ProdMolName = "Prod%d" % ProdMolCount
                ProdMol.SetProp("_Name", ProdMolName)
                
                Status = WriteProductMolecule(Writer, ProdMol, Sanitize, Compute2DCoords)
                if Status:
                    ValidProdMolCount += 1

    if Writer is not None:
        Writer.close()
    
    if ValidProdMolCount:
        MiscUtil.PrintInfo("\nTotal number of product molecules: %d" % ProdMolCount)
        MiscUtil.PrintInfo("Number of valid product molecules: %d" % ValidProdMolCount)
        MiscUtil.PrintInfo("Number of ignored product molecules: %d" % (ProdMolCount - ValidProdMolCount))
    else:
        MiscUtil.PrintInfo("\nThe compound library enumeration failed to generate any product molecules.\nCheck to make sure the reactants specified in input files match their corresponding specifications in reaction SMIRKS and try again.")

def WriteProductMolecule(Writer, ProdMol, Sanitize, Compute2DCoords):
    """Prepare and write out product  molecule."""

    try:
        if Sanitize:
            Chem.SanitizeMol(ProdMol)
    except (RuntimeError, ValueError):
        MiscUtil.PrintWarning("Ignoring product molecule: Failed to sanitize...\n")
        return False

    try:
        if Compute2DCoords:
            AllChem.Compute2DCoords(ProdMol)
    except (RuntimeError, ValueError):
        MiscUtil.PrintWarning("Ignoring product molecule: Failed to compute 2D coordinates...\n")
        return False

    Writer.write(ProdMol)

    return True

def RetrieveReactantsMolecules():
    """Retrieve reactant molecules from each reactant file and return a list containing lists of molecules
    for each reactant file."""

    MiscUtil.PrintInfo("\nProcessing reactant file(s)...")
    
    ReactantsMolsList = []
    ReactantFilesList = OptionsInfo["ReactantFilesList"]
    UseReactantNames = OptionsInfo["UseReactantNames"]
    ReactantCount = 0
    
    for FileIndex in range(0, len(ReactantFilesList)):
        ReactantCount += 1
        ReactantFile = ReactantFilesList[FileIndex]
        
        MiscUtil.PrintInfo("\nProcessing reactant file: %s..." % ReactantFile)

        Mols  = RDKitUtil.ReadMolecules(ReactantFile, **OptionsInfo["InfileParams"])
        
        ValidMols = []
        MolCount = 0
        ValidMolCount = 0
        
        for Mol in Mols:
            MolCount += 1
            if Mol is None:
                continue
            
            if RDKitUtil.IsMolEmpty(Mol):
                MolName = RDKitUtil.GetMolName(Mol, MolCount)
                MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
                continue
            
            ValidMolCount += 1

            # Check and set mol name...
            if UseReactantNames:
                MolName = RDKitUtil.GetMolName(Mol)
                if not len(MolName):
                    MolName = "React%dMol%d" % (ReactantCount, MolCount)
                    Mol.SetProp("_Name", MolName)
                
            ValidMols.append(Mol)

        ReactantsMolsList.append(ValidMols)
        
        MiscUtil.PrintInfo("Total number of molecules: %d" % MolCount)
        MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
        MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))
    
    return ReactantsMolsList
    
def ProcessReactionNamesInfo():
    """Process reaction names information."""
    
    if not OptionsInfo["RxnByNameMode"]:
        return

    RetrieveReactionNamesInfo()
    ProcessSpecifiedReactionName()

def ProcessSpecifiedReactionName():
    """Process and validate specified reaction name."""

    OptionsInfo["RxnNameSMIRKS"] = None
    
    # Set up a map of valid group rxn names for checking specified rxn names...
    CanonicalRxnNameMap = {}
    for Name in RxnNamesMap['Names']:
        CanonicalRxnNameMap[Name.lower()] = Name
    
    CanonicalRxnName = OptionsInfo["RxnName"].lower()
    if CanonicalRxnName in CanonicalRxnNameMap:
        Name = CanonicalRxnNameMap[CanonicalRxnName]
        OptionsInfo["RxnNameSMIRKS"] = RxnNamesMap['SMIRKSPattern'][Name]
    else:
        MiscUtil.PrintError("The rxn name name, %s, specified using \"-r, --rxnName\" option is not a valid name." % (OptionsInfo["RxnName"]))
    
def ProcessListReactionNamesOption():
    """Process list reaction names information."""

    # Validate and process dataFile option for listing reaction names information...
    OptionsInfo["RxnNamesFile"] = None
    if not re.match("^auto$", Options["--rxnNamesFile"], re.I):
        MiscUtil.ValidateOptionFilePath("--rxnNamesFile", Options["--rxnNamesFile"])
        OptionsInfo["RxnNamesFile"] = Options["--rxnNamesFile"]
    
    RetrieveReactionNamesInfo()
    ListReactionNamesInfo()

def RetrieveReactionNamesInfo():
    """Retrieve reaction names information."""

    RxnNamesFilePath = OptionsInfo["RxnNamesFile"]
    if RxnNamesFilePath is None:
        MayaChemToolsDataDir = MiscUtil.GetMayaChemToolsLibDataPath()
        RxnNamesFilePath = os.path.join(MayaChemToolsDataDir, "ReactionNamesAndSMIRKS.csv")
        
    MiscUtil.PrintInfo("\nRetrieving reaction names and SMIRKS patterns from file %s" % (RxnNamesFilePath))
    
    if not os.path.exists(RxnNamesFilePath):
        MiscUtil.PrintError("The reaction names file, %s, doesn't exist.\n" % (RxnNamesFilePath))

    Delimiter = ','
    QuoteChar = '"'
    IgnoreHeaderLine = True
    RxnLinesWords = MiscUtil.GetTextLinesWords(RxnNamesFilePath, Delimiter, QuoteChar, IgnoreHeaderLine)
    
    RxnNamesMap['Names'] = []
    RxnNamesMap['SMIRKSPattern'] = {}
    
    for LineWords in RxnLinesWords:
        Name = LineWords[0]
        SMIRKSPattern = LineWords[1]

        if Name in RxnNamesMap['SMIRKSPattern']:
            MiscUtil.PrintWarning("Ignoring duplicate reaction name: %s..." % Name)
        else:
            RxnNamesMap['Names'].append(Name)
            RxnNamesMap['SMIRKSPattern'][Name] = SMIRKSPattern
        
    if not len(RxnNamesMap['Names']):
        MiscUtil.PrintError("Failed to retrieve any reaction names and SMIRKS patterns...")
        
    MiscUtil.PrintInfo("Total number of reactions present in reaction names and SMIRKS file: %d" % (len(RxnNamesMap['Names'])))

def ListReactionNamesInfo():
    """List reaction names information"""

    MiscUtil.PrintInfo("\nListing available freaction names and SMIRKS patterns...")
    MiscUtil.PrintInfo("\nReactionName\tSMIRKSPattern")
    
    for Name in sorted(RxnNamesMap['Names']):
        SMIRKSPattern = RxnNamesMap['SMIRKSPattern'][Name]
        MiscUtil.PrintInfo("%s\t%s" % (Name, SMIRKSPattern))

    MiscUtil.PrintInfo("")

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    Compute2DCoords = True
    if not re.match("^yes$", Options["--compute2DCoords"], re.I):
        Compute2DCoords = False
    OptionsInfo["Compute2DCoords"]  = Compute2DCoords

    OptionsInfo["Mode"] = Options["--mode"]
    RxnByNameMode = True
    if not re.match("^RxnByName$", Options["--mode"], re.I):
        RxnByNameMode = False
    OptionsInfo["RxnByNameMode"] = RxnByNameMode

    OptionsInfo["ProdMolNamesMode"] = Options["--prodMolNames"]
    UseReactantNames = False
    if re.match("^UseReactants$", Options["--prodMolNames"], re.I):
        UseReactantNames = True
    OptionsInfo["UseReactantNames"] = UseReactantNames
    
    OptionsInfo["RxnName"] = Options["--rxnName"]
    OptionsInfo["RxnNameSMIRKS"] = None
    if OptionsInfo["RxnByNameMode"]:
        if not Options["--rxnName"]:
            MiscUtil.PrintError("No rxn name specified using \"-r, --rxnName\" option during \"RxnByName\" value of \"-m, --mode\" option")

    OptionsInfo["RxnNamesFile"] = None
    if not re.match("^auto$", Options["--rxnNamesFile"], re.I):
        OptionsInfo["RxnNamesFile"] = Options["--rxnNamesFile"]

    ReactantFiles = re.sub(" ", "", Options["--infiles"])
    ReactantFilesList = []
    ReactantFilesList = ReactantFiles.split(",")
    OptionsInfo["ReactantFiles"] = ReactantFiles
    OptionsInfo["ReactantFilesList"] = ReactantFilesList

    OptionsInfo["SpecifiedSMIRKS"] = Options["--smirksRxn"]
    if not OptionsInfo["RxnByNameMode"]:
        if not Options["--smirksRxn"]:
            MiscUtil.PrintError("No rxn SMIRKS pattern specified using \"-r, --rxnName\" option during \"RxnByName\" value of \"-m, --mode\" option")
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["Overwrite"] = Options["--overwrite"]

    # Use first reactant file as input file as all input files have the same format...
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], ReactantFilesList[0])

    # No need to pass any input or output file name due to absence of any auto parameter...
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"])
    
    Sanitize = True
    if not re.match("^yes$", Options["--sanitize"], re.I):
        Sanitize = False
    OptionsInfo["Sanitize"]  = Sanitize

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
    
    # Handle listing of functional group information...
    if  Options and Options["--list"]:
        ProcessListReactionNamesOption()
        sys.exit(0)

def ValidateOptions():
    """Validate option values"""
    
    MiscUtil.ValidateOptionTextValue("--compute2DCoords", Options["--compute2DCoords"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "RxnByName RxnBySMIRKS")
    MiscUtil.ValidateOptionTextValue("-p, --prodMolNames", Options["--prodMolNames"], "UseReactants Sequential")
    
    if not re.match("^auto$", Options["--rxnNamesFile"], re.I):
        MiscUtil.ValidateOptionFilePath("--rxnNamesFile", Options["--rxnNamesFile"])

    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd smi")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    
    ReactantFiles = re.sub(" ", "", Options["--infiles"])
    if not ReactantFiles:
        MiscUtil.PrintError("No reactant files specified for \"-i, --infiles\" option")

    # Validate file extensions...
    for ReactantFile in ReactantFiles.split(","):
        MiscUtil.ValidateOptionFilePath("-i, --infiles", ReactantFile)
        MiscUtil.ValidateOptionFileExt("-i, --infiles", ReactantFile, "sdf sd smi csv tsv txt")
        MiscUtil.ValidateOptionsDistinctFileNames("-i, --infiles", ReactantFile, "-o, --outfile", Options["--outfile"])
        MiscUtil.ValidateOptionsDistinctFileNames("-i, --infiles", ReactantFile, "-o, --outfile", Options["--outfile"])
        
    # Match file formats...
    FirstFile = True
    FirstFileFormat = ""
    for ReactantFile in ReactantFiles.split(","):
        FileFormat = ""
        if MiscUtil.CheckFileExt(ReactantFile, "sdf sd"):
            FileFormat = "SD"
        elif MiscUtil.CheckFileExt(ReactantFile, "smi csv tsv txt"):
            FileFormat = "SMILES"
        else:
            MiscUtil.PrintError("The file name specified , %s, for option \"-i, --infiles\" is not valid. Supported file formats: sdf sd smi csv tsv txt\n" % ReactantFile)
            
        if FirstFile:
            FirstFile = False
            FirstFileFormat = FileFormat
            continue
        
        if not re.match("^%s$" % FirstFileFormat, FileFormat, re.IGNORECASE):
            MiscUtil.PrintError("All reactant file names -  %s - specified using option \"-i, --infiles\" must have the same file format.\n" % ReactantFiles)
            

    MiscUtil.ValidateOptionTextValue("--sanitize", Options["--sanitize"], "yes no")
    
# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitEnumerateCompoundLibrary.py - Enumerate a virtual compound library

Usage:
    RDKitEnumerateCompoundLibrary.py  [--compute2DCoords <yes or no>] [--infileParams <Name,Value,...>]
                                      [--mode <RxnByName or RxnBySMIRKS>] [--outfileParams <Name,Value,...>] [--overwrite]
                                      [--prodMolNames <UseReactants or Sequential>] [--rxnName <text>]
                                      [--rxnNamesFile <FileName or auto>] [--smirksRxn <text>] [--sanitize <yes or no>]
                                      [-w <dir>] -i  <ReactantFile1,...> -o <outfile>
    RDKitEnumerateCompoundLibrary.py [--rxnNamesFile <FileName or auto>] -l | --list
    RDKitEnumerateCompoundLibrary.py -h | --help | -e | --examples

Description:
    Perform a combinatorial enumeration of a virtual library of molecules for a reaction specified
    using a reaction name or SMIRKS pattern and reactant input files.

    The SMIRKS patterns for supported reactions names [ Ref 134 ] are retrieved from file,
    ReactionNamesAndSMIRKS.csv, available in MayaChemTools data directory. The current
    list of supported reaction names is shown below:

    '1,2,4_triazole_acetohydrazide', '1,2,4_triazole_carboxylic_acid_ester', 3_nitrile_pyridine,
    Benzimidazole_derivatives_aldehyde, Benzimidazole_derivatives_carboxylic_acid_ester,
    Benzofuran, Benzothiazole, Benzothiophene, Benzoxazole_aromatic_aldehyde,
    Benzoxazole_carboxylic_acid, Buchwald_Hartwig, Decarboxylative_coupling, Fischer_indole,
    Friedlaender_chinoline, Grignard_alcohol, Grignard_carbonyl, Heck_non_terminal_vinyl,
    Heck_terminal_vinyl, Heteroaromatic_nuc_sub, Huisgen_Cu_catalyzed_1,4_subst,
    Huisgen_disubst_alkyne, Huisgen_Ru_catalyzed_1,5_subst, Imidazole, Indole, Mitsunobu_imide,
    Mitsunobu_phenole, Mitsunobu_sulfonamide, Mitsunobu_tetrazole_1, Mitsunobu_tetrazole_2,
    Mitsunobu_tetrazole_3, Mitsunobu_tetrazole_4, N_arylation_heterocycles, Negishi,
    Niementowski_quinazoline, Nucl_sub_aromatic_ortho_nitro, Nucl_sub_aromatic_para_nitro,
    Oxadiazole, Paal_Knorr_pyrrole, Phthalazinone, Pictet_Spengler, Piperidine_indole,
    Pyrazole, Reductive_amination, Schotten_Baumann_amide, Sonogashira, Spiro_chromanone,
    Stille, Sulfon_amide, Suzuki, Tetrazole_connect_regioisomer_1, Tetrazole_connect_regioisomer_2,
    Tetrazole_terminal, Thiazole, Thiourea, Triaryl_imidazole, Urea, Williamson_ether, Wittig 

    The supported input file formats are: SD (.sdf, .sd), SMILES (.smi, .csv, .tsv, .txt)

    The supported output file formats are:  SD (.sdf, .sd), SMILES (.smi)

Options:
    -c, --compute2DCoords <yes or no>  [default: yes]
        Compute 2D coordinates of product molecules before writing them out.
    -i, --infiles <ReactantFile1, ReactantFile2...>
        Comma delimited list of reactant file names for enumerating a compound library
        using reaction SMIRKS. The number of reactant files must match number of
        reaction components in reaction SMIRKS. All reactant input files must have
        the same format.
    --infileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for reading
        molecules from files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD, MOL: removeHydrogens,yes,sanitize,yes,strictParsing,yes
            SMILES: smilesColumn,1,smilesNameColumn,2,smilesDelimiter,space,
                smilesTitleLine,auto,sanitize,yes
            
        Possible values for smilesDelimiter: space, comma or tab. These parameters apply
        to all reactant input files, which must have the same file format.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -l, --list
        List available reaction names along with corresponding SMIRKS patterns without
        performing any enumeration.
    -m, --mode <RxnByName or RxnBySMIRKS>  [default: RxnByName]
        Indicate whether a reaction is specified by a reaction name or a SMIRKS pattern.
        Possible values: RxnByName or RxnBySMIRKS.
    -o, --outfile <outfile>
        Output file name.
    --outfileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for writing
        molecules to files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD: kekulize,no
            SMILES: kekulize,no,smilesDelimiter,space, smilesIsomeric,yes,
                smilesTitleLine,yes
            
    -p, --prodMolNames <UseReactants or Sequential>  [default: UseReactants]
        Generate names of product molecules using reactant names or assign names in
        a sequential order. Possible values: UseReactants or Sequential. Format of
        molecule names: UseReactants - <ReactName1>_<ReactName2>..._Prod<Num>;
        Sequential - Prod<Num>
    --overwrite
        Overwrite existing files.
    -r, --rxnName <text>
        Name of a reaction to use for enumerating a compound library. This option
        is only used during 'RxnByName' value of '-m, --mode' option.
    --rxnNamesFile <FileName or auto>  [default: auto]
        Specify a file name containing data for names of reactions and SMIRKS patterns or
        use default file, ReactionNamesAndSMIRKS.csv, available in MayaChemTools data
        directory.
        
        Reactions SMIRKS file format: RxnName,RxnSMIRKS.
        
        The format of data in local reaction names file must match format of the reaction
        SMIRKS file available in MayaChemTools data directory.
    -s, --smirksRxn <text>
        SMIRKS pattern of a reaction to use for enumerating a compound library. This
        option is only used during 'RxnBySMIRKS' value of '-m, --mode' option.
    --sanitize <yes or no>  [default: yes]
        Sanitize product molecules before writing them out.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To list all available reaction names along with their SMIRKS pattern, type:

         % RDKitEnumerateCompoundLibrary.py -l

    To perform a combinatorial enumeration of a virtual compound library corresponding
    to named amide reaction, Schotten_Baumann_amide and write out a SMILES file
    type:

        % RDKitEnumerateCompoundLibrary.py -r Schotten_Baumann_amide
          -i 'SampleAcids.smi,SampleAmines.smi' -o SampleOutCmpdLibrary.smi

    To perform a combinatorial enumeration of a virtual compound library corresponding
    to an amide reaction specified using a SMIRKS pattern and write out a SD file containing
    sanitized molecules, computed 2D coordinates, and generation of molecule names from
    reactant names, type:

        % RDKitEnumerateCompoundLibrary.py -m RxnBySMIRKS
          -s '[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]'
          -i 'SampleAcids.smi,SampleAmines.smi' -o SampleOutCmpdLibrary.sdf

    To perform a combinatorial enumeration of a virtual compound library corresponding
    to an amide reaction specified using a SMIRKS pattern  and write out a SD file containing
    unsanitized molecules, without generating 2D coordinates, and a sequential generation
    of molecule names, type:

        % RDKitEnumerateCompoundLibrary.py -m RxnBySMIRKS -c no -s no
          -p Sequential -s '[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]'
          -i 'SampleAcids.smi,SampleAmines.smi' -o SampleOutCmpdLibrary.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitFilterPAINS.py, RDKitSearchFunctionalGroups.py,
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
