#!/usr/bin/perl -w
#
# File: PathLengthFingerprints.pl
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2018 Manish Sud. All rights reserved.
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

use strict;
use FindBin; use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use File::Basename;
use Text::ParseWords;
use Benchmark;
use FileUtil;
use TextUtil;
use SDFileUtil;
use MoleculeFileIO;
use FileIO::FingerprintsSDFileIO;
use FileIO::FingerprintsTextFileIO;
use FileIO::FingerprintsFPFileIO;
use AtomTypes::AtomicInvariantsAtomTypes;
use AtomTypes::FunctionalClassAtomTypes;
use Fingerprints::PathLengthFingerprints;

my($ScriptName, %Options, $StartTime, $EndTime, $TotalTime);

# Autoflush STDOUT
$| = 1;

# Starting message...
$ScriptName = basename($0);
print "\n$ScriptName: Starting...\n\n";
$StartTime = new Benchmark;

# Get the options and setup script...
SetupScriptUsage();
if ($Options{help} || @ARGV < 1) {
  die GetUsageFromPod("$FindBin::Bin/$ScriptName");
}

my(@SDFilesList);
@SDFilesList = ExpandFileNames(\@ARGV, "sdf sd");

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Setup information about input files...
print "Checking input SD file(s)...\n";
my(%SDFilesInfo);
RetrieveSDFilesInfo();

# Process input files..
my($FileIndex);
if (@SDFilesList > 1) {
  print "\nProcessing SD files...\n";
}
for $FileIndex (0 .. $#SDFilesList) {
  if ($SDFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $SDFilesList[$FileIndex]...\n";
    GeneratePathLengthFingerprints($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Generate fingerprints for a SD file...
#
sub GeneratePathLengthFingerprints {
  my($FileIndex) = @_;
  my($CmpdCount, $IgnoredCmpdCount, $SDFile, $MoleculeFileIO, $Molecule, $PathLengthFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO);

  $SDFile = $SDFilesList[$FileIndex];

  # Setup output files...
  #
  ($NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO) = SetupAndOpenOutputFiles($FileIndex);

  $MoleculeFileIO = new MoleculeFileIO('Name' => $SDFile);
  $MoleculeFileIO->Open();

  $CmpdCount = 0;
  $IgnoredCmpdCount = 0;

  COMPOUND: while ($Molecule = $MoleculeFileIO->ReadMolecule()) {
    $CmpdCount++;

    # Filter compound data before calculating fingerprints...
    if ($OptionsInfo{Filter}) {
      if (CheckAndFilterCompound($CmpdCount, $Molecule)) {
	$IgnoredCmpdCount++;
	next COMPOUND;
      }
    }

    $PathLengthFingerprints = GenerateMoleculeFingerprints($Molecule);
    if (!$PathLengthFingerprints) {
      $IgnoredCmpdCount++;
      ProcessIgnoredCompound('FingerprintsGenerationFailed', $CmpdCount, $Molecule);
      next COMPOUND;
    }

    WriteDataToOutputFiles($FileIndex, $CmpdCount, $Molecule, $PathLengthFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO);
  }
  $MoleculeFileIO->Close();

  if ($NewFPSDFileIO) {
    $NewFPSDFileIO->Close();
  }
  if ($NewFPTextFileIO) {
    $NewFPTextFileIO->Close();
  }
  if ($NewFPFileIO) {
    $NewFPFileIO->Close();
  }

  WriteFingerprintsGenerationSummaryStatistics($CmpdCount, $IgnoredCmpdCount);
}

# Process compound being ignored due to problems in fingerprints geneation...
#
sub ProcessIgnoredCompound {
  my($Mode, $CmpdCount, $Molecule) = @_;
  my($CmpdID, $DataFieldLabelAndValuesRef);

  $DataFieldLabelAndValuesRef = $Molecule->GetDataFieldLabelAndValues();
  $CmpdID = SetupCmpdIDForOutputFiles($CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);

  MODE: {
    if ($Mode =~ /^ContainsNonElementalData$/i) {
      warn "\nWarning: Ignoring compound record number $CmpdCount with ID $CmpdID: Compound contains atom data corresponding to non-elemental atom symbol(s)...\n\n";
      next MODE;
    }

    if ($Mode =~ /^ContainsNoElementalData$/i) {
      warn "\nWarning: Ignoring compound record number $CmpdCount with ID $CmpdID: Compound contains no atom data...\n\n";
      next MODE;
    }

    if ($Mode =~ /^FingerprintsGenerationFailed$/i) {
      warn "\nWarning: Ignoring compound record number $CmpdCount with ID $CmpdID: Fingerprints generation didn't succeed...\n\n";
      next MODE;
    }
    warn "\nWarning: Ignoring compound record number $CmpdCount with ID $CmpdID: Fingerprints generation didn't succeed...\n\n";
  }
}

# Check and filter compounds....
#
sub CheckAndFilterCompound {
  my($CmpdCount, $Molecule) = @_;
  my($ElementCount, $NonElementCount);

  ($ElementCount, $NonElementCount) = $Molecule->GetNumOfElementsAndNonElements();

  if ($NonElementCount) {
    ProcessIgnoredCompound('ContainsNonElementalData', $CmpdCount, $Molecule);
    return 1;
  }

  if (!$ElementCount) {
    ProcessIgnoredCompound('ContainsNoElementalData', $CmpdCount, $Molecule);
    return 1;
  }

  return 0;
}

# Write out compounds fingerprints generation summary statistics...
#
sub WriteFingerprintsGenerationSummaryStatistics {
  my($CmpdCount, $IgnoredCmpdCount) = @_;
  my($ProcessedCmpdCount);

  $ProcessedCmpdCount = $CmpdCount - $IgnoredCmpdCount;

  print "\nNumber of compounds: $CmpdCount\n";
  print "Number of compounds processed successfully during fingerprints generation: $ProcessedCmpdCount\n";
  print "Number of compounds ignored during fingerprints generation: $IgnoredCmpdCount\n";
}

# Open output files...
#
sub SetupAndOpenOutputFiles {
  my($FileIndex) = @_;
  my($NewFPSDFile, $NewFPFile, $NewFPTextFile, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO, %FingerprintsFileIOParams);

  ($NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO) = (undef) x 3;

  # Setup common parameters for fingerprints file IO objects...
  #
  %FingerprintsFileIOParams = ();
  if ($OptionsInfo{Mode} =~ /^PathLengthBits$/i) {
    %FingerprintsFileIOParams = ('Mode' => 'Write', 'Overwrite' => $OptionsInfo{OverwriteFiles}, 'FingerprintsStringMode' => 'FingerprintsBitVectorString', 'BitStringFormat' => $OptionsInfo{BitStringFormat}, 'BitsOrder' => $OptionsInfo{BitsOrder});
  }
  elsif ($OptionsInfo{Mode} =~ /^PathLengthCount$/i) {
    %FingerprintsFileIOParams = ('Mode' => 'Write', 'Overwrite' => $OptionsInfo{OverwriteFiles}, 'FingerprintsStringMode' => 'FingerprintsVectorString', 'VectorStringFormat' => $OptionsInfo{VectorStringFormat});
  }

  if ($OptionsInfo{SDOutput}) {
    $NewFPSDFile = $SDFilesInfo{SDOutFileNames}[$FileIndex];
    print "Generating SD file $NewFPSDFile...\n";
    $NewFPSDFileIO = new FileIO::FingerprintsSDFileIO('Name' => $NewFPSDFile, %FingerprintsFileIOParams, 'FingerprintsFieldLabel' => $OptionsInfo{FingerprintsLabel});
    $NewFPSDFileIO->Open();
  }

  if ($OptionsInfo{FPOutput}) {
    $NewFPFile = $SDFilesInfo{FPOutFileNames}[$FileIndex];
    print "Generating FP file $NewFPFile...\n";
    $NewFPFileIO = new FileIO::FingerprintsFPFileIO('Name' => $NewFPFile, %FingerprintsFileIOParams);
    $NewFPFileIO->Open();
  }

  if ($OptionsInfo{TextOutput}) {
    my($ColLabelsRef);

    $NewFPTextFile = $SDFilesInfo{TextOutFileNames}[$FileIndex];
    $ColLabelsRef = SetupFPTextFileCoulmnLabels($FileIndex);

    print "Generating text file $NewFPTextFile...\n";
    $NewFPTextFileIO = new FileIO::FingerprintsTextFileIO('Name' => $NewFPTextFile, %FingerprintsFileIOParams, 'DataColLabels' => $ColLabelsRef, 'OutDelim' => $OptionsInfo{OutDelim}, 'OutQuote' => $OptionsInfo{OutQuote});
    $NewFPTextFileIO->Open();
  }

  return ($NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO);
}

# Write fingerpritns and other data to appropriate output files...
#
sub WriteDataToOutputFiles {
  my($FileIndex, $CmpdCount, $Molecule, $PathLengthFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO) = @_;
  my($DataFieldLabelAndValuesRef);

  $DataFieldLabelAndValuesRef = undef;
  if ($NewFPTextFileIO || $NewFPFileIO) {
    $DataFieldLabelAndValuesRef = $Molecule->GetDataFieldLabelAndValues();
  }

  if ($NewFPSDFileIO) {
    my($CmpdString);

    $CmpdString = $Molecule->GetInputMoleculeString();
    $NewFPSDFileIO->WriteFingerprints($PathLengthFingerprints, $CmpdString);
  }

  if ($NewFPTextFileIO) {
    my($ColValuesRef);

    $ColValuesRef = SetupFPTextFileCoulmnValues($FileIndex, $CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPTextFileIO->WriteFingerprints($PathLengthFingerprints, $ColValuesRef);
  }

  if ($NewFPFileIO) {
    my($CompoundID);

    $CompoundID = SetupCmpdIDForOutputFiles($CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPFileIO->WriteFingerprints($PathLengthFingerprints, $CompoundID);
  }
}

# Generate approriate column labels for FPText output file...
#
sub SetupFPTextFileCoulmnLabels {
  my($FileIndex) = @_;
  my($Line, @ColLabels);

  @ColLabels = ();
  if ($OptionsInfo{DataFieldsMode} =~ /^All$/i) {
    push @ColLabels, @{$SDFilesInfo{AllDataFieldsRef}[$FileIndex]};
  }
  elsif ($OptionsInfo{DataFieldsMode} =~ /^Common$/i) {
    push @ColLabels, @{$SDFilesInfo{CommonDataFieldsRef}[$FileIndex]};
  }
  elsif ($OptionsInfo{DataFieldsMode} =~ /^Specify$/i) {
    push @ColLabels, @{$OptionsInfo{SpecifiedDataFields}};
  }
  elsif ($OptionsInfo{DataFieldsMode} =~ /^CompoundID$/i) {
    push @ColLabels, $OptionsInfo{CompoundIDLabel};
  }
  # Add fingerprints label...
  push @ColLabels, $OptionsInfo{FingerprintsLabel};

  return \@ColLabels;
}

# Generate column values FPText output file..
#
sub SetupFPTextFileCoulmnValues {
  my($FileIndex, $CmpdCount, $Molecule, $DataFieldLabelAndValuesRef) = @_;
  my(@ColValues);

  @ColValues = ();
  if ($OptionsInfo{DataFieldsMode} =~ /^CompoundID$/i) {
    push @ColValues, SetupCmpdIDForOutputFiles($CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
  }
  elsif ($OptionsInfo{DataFieldsMode} =~ /^All$/i) {
    @ColValues = map { exists $DataFieldLabelAndValuesRef->{$_} ? $DataFieldLabelAndValuesRef->{$_} : ''} @{$SDFilesInfo{AllDataFieldsRef}[$FileIndex]};
  }
  elsif ($OptionsInfo{DataFieldsMode} =~ /^Common$/i) {
    @ColValues = map { exists $DataFieldLabelAndValuesRef->{$_} ? $DataFieldLabelAndValuesRef->{$_} : ''} @{$SDFilesInfo{CommonDataFieldsRef}[$FileIndex]};
  }
  elsif ($OptionsInfo{DataFieldsMode} =~ /^Specify$/i) {
    @ColValues = map { exists $DataFieldLabelAndValuesRef->{$_} ? $DataFieldLabelAndValuesRef->{$_} : ''} @{$OptionsInfo{SpecifiedDataFields}};
  }

  return \@ColValues;
}

# Generate compound ID for FP and FPText output files..
#
sub SetupCmpdIDForOutputFiles {
  my($CmpdCount, $Molecule, $DataFieldLabelAndValuesRef) = @_;
  my($CmpdID);

  $CmpdID = '';
  if ($OptionsInfo{CompoundIDMode} =~ /^MolNameOrLabelPrefix$/i) {
    my($MolName);
    $MolName = $Molecule->GetName();
    $CmpdID = $MolName ? $MolName : "$OptionsInfo{CompoundID}${CmpdCount}";
  }
  elsif ($OptionsInfo{CompoundIDMode} =~ /^LabelPrefix$/i) {
    $CmpdID = "$OptionsInfo{CompoundID}${CmpdCount}";
  }
  elsif ($OptionsInfo{CompoundIDMode} =~ /^DataField$/i) {
    my($SpecifiedDataField);
    $SpecifiedDataField = $OptionsInfo{CompoundID};
    $CmpdID = exists $DataFieldLabelAndValuesRef->{$SpecifiedDataField} ? $DataFieldLabelAndValuesRef->{$SpecifiedDataField} : '';
  }
  elsif ($OptionsInfo{CompoundIDMode} =~ /^MolName$/i) {
    $CmpdID = $Molecule->GetName();
  }
  return $CmpdID;
}

# Generate fingerprints for molecule...
#
sub GenerateMoleculeFingerprints {
  my($Molecule) = @_;
  my($PathLengthFingerprints);

  if ($OptionsInfo{KeepLargestComponent}) {
    $Molecule->KeepLargestComponent();
  }
  if ($OptionsInfo{IgnoreHydrogens}) {
    $Molecule->DeleteHydrogens();
  }

  if ($OptionsInfo{DetectAromaticity}) {
    if (!$Molecule->DetectRings()) {
      return undef;
    }
    $Molecule->SetAromaticityModel($OptionsInfo{AromaticityModel});
    $Molecule->DetectAromaticity();
  }

  $PathLengthFingerprints = undef;
  if ($OptionsInfo{Mode} =~ /^PathLengthBits$/i) {
    $PathLengthFingerprints = GeneratePathLengthBitsFingerprints($Molecule);
  }
  elsif ($OptionsInfo{Mode} =~ /^PathLengthCount$/i) {
    $PathLengthFingerprints = GeneratePathLengthCountFingerprints($Molecule);
  }
  else {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: PathLengthBits or PathLengthCount\n";
  }

  return $PathLengthFingerprints;
}

# Generate pathlength bits finerprints for molecule...
#
sub GeneratePathLengthBitsFingerprints {
  my($Molecule) = @_;
  my($PathLengthFingerprints);

  $PathLengthFingerprints = new Fingerprints::PathLengthFingerprints('Molecule' => $Molecule, 'Type' => 'PathLengthBits', 'AtomIdentifierType' => $OptionsInfo{AtomIdentifierType}, 'NumOfBitsToSetPerPath' => $OptionsInfo{NumOfBitsToSetPerPath}, 'Size' => $OptionsInfo{Size}, 'MinLength' => $OptionsInfo{MinPathLength}, 'MaxLength' => $OptionsInfo{MaxPathLength}, 'AllowRings' => $OptionsInfo{AllowRings}, 'AllowSharedBonds' => $OptionsInfo{AllowSharedBonds}, 'UseBondSymbols' => $OptionsInfo{UseBondSymbols}, 'UseUniquePaths' => $OptionsInfo{UseUniquePaths}, 'UsePerlCoreRandom' => $OptionsInfo{UsePerlCoreRandom});

  # Set atom identifier type...
  SetAtomIdentifierTypeValuesToUse($PathLengthFingerprints);

  # Generate fingerprints...
  $PathLengthFingerprints->GenerateFingerprints();

  # Make sure fingerprints generation is successful...
  if (!$PathLengthFingerprints->IsFingerprintsGenerationSuccessful()) {
    return undef;
  }

  if ($OptionsInfo{Fold}) {
    my($CheckSizeValue) = 0;
    $PathLengthFingerprints->FoldFingerprintsBySize($OptionsInfo{FoldedSize}, $CheckSizeValue);
  }

  return $PathLengthFingerprints;
}

# Generate pathlength count finerprints for molecule...
#
sub GeneratePathLengthCountFingerprints {
  my($Molecule) = @_;
  my($PathLengthFingerprints);

  $PathLengthFingerprints = new Fingerprints::PathLengthFingerprints('Molecule' => $Molecule, 'Type' => 'PathLengthCount', 'AtomIdentifierType' => $OptionsInfo{AtomIdentifierType}, 'MinLength' => $OptionsInfo{MinPathLength}, 'MaxLength' => $OptionsInfo{MaxPathLength}, 'AllowRings' => $OptionsInfo{AllowRings}, 'AllowSharedBonds' => $OptionsInfo{AllowSharedBonds}, 'UseBondSymbols' => $OptionsInfo{UseBondSymbols}, 'UseUniquePaths' => $OptionsInfo{UseUniquePaths});

  # Set atom identifier type...
  SetAtomIdentifierTypeValuesToUse($PathLengthFingerprints);

  # Generate fingerprints...
  $PathLengthFingerprints->GenerateFingerprints();

  # Make sure fingerprints generation is successful...
  if (!$PathLengthFingerprints->IsFingerprintsGenerationSuccessful()) {
    return undef;
  }
  return $PathLengthFingerprints;
}

# Set atom identifier type to use for generating path strings...
#
sub SetAtomIdentifierTypeValuesToUse {
  my($PathLengthFingerprints) = @_;

  if ($OptionsInfo{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
    $PathLengthFingerprints->SetAtomicInvariantsToUse(\@{$OptionsInfo{AtomicInvariantsToUse}});
  }
  elsif ($OptionsInfo{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
    $PathLengthFingerprints->SetFunctionalClassesToUse(\@{$OptionsInfo{FunctionalClassesToUse}});
  }
  elsif ($OptionsInfo{AtomIdentifierType} =~ /^(DREIDINGAtomTypes|EStateAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
    # Nothing to do for now...
  }
  else {
    die "Error: The value specified, $Options{atomidentifiertype}, for option \"-a, --AtomIdentifierType\" is not valid. Supported atom identifier types in current release of MayaChemTools: AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes\n";
  }
}

# Retrieve information about SD files...
#
sub RetrieveSDFilesInfo {
  my($SDFile, $Index, $FileDir, $FileExt, $FileName, $OutFileRoot, $TextOutFileExt, $SDOutFileExt, $FPOutFileExt, $NewSDFileName, $NewFPFileName, $NewTextFileName, $CheckDataField, $CollectDataFields, $AllDataFieldsRef, $CommonDataFieldsRef);

  %SDFilesInfo = ();
  @{$SDFilesInfo{FileOkay}} = ();
  @{$SDFilesInfo{OutFileRoot}} = ();
  @{$SDFilesInfo{SDOutFileNames}} = ();
  @{$SDFilesInfo{FPOutFileNames}} = ();
  @{$SDFilesInfo{TextOutFileNames}} = ();
  @{$SDFilesInfo{AllDataFieldsRef}} = ();
  @{$SDFilesInfo{CommonDataFieldsRef}} = ();

  $CheckDataField = ($OptionsInfo{TextOutput} && ($OptionsInfo{DataFieldsMode} =~ /^CompoundID$/i) && ($OptionsInfo{CompoundIDMode} =~ /^DataField$/i)) ? 1 : 0;
  $CollectDataFields = ($OptionsInfo{TextOutput} && ($OptionsInfo{DataFieldsMode} =~ /^(All|Common)$/i)) ? 1 : 0;

  FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFile = $SDFilesList[$Index];

    $SDFilesInfo{FileOkay}[$Index] = 0;
    $SDFilesInfo{OutFileRoot}[$Index] = '';
    $SDFilesInfo{SDOutFileNames}[$Index] = '';
    $SDFilesInfo{FPOutFileNames}[$Index] = '';
    $SDFilesInfo{TextOutFileNames}[$Index] = '';

    $SDFile = $SDFilesList[$Index];
    if (!(-e $SDFile)) {
      warn "Warning: Ignoring file $SDFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($SDFile, "sd sdf")) {
      warn "Warning: Ignoring file $SDFile: It's not a SD file\n";
      next FILELIST;
    }

    if ($CheckDataField) {
      # Make sure data field exists in SD file..
      my($CmpdString, $SpecifiedDataField, @CmpdLines, %DataFieldValues);

      @CmpdLines = ();
      open SDFILE, "$SDFile" or die "Error: Couldn't open $SDFile: $! \n";
      $CmpdString = ReadCmpdString(\*SDFILE);
      close SDFILE;
      @CmpdLines = split "\n", $CmpdString;
      %DataFieldValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
      $SpecifiedDataField = $OptionsInfo{CompoundID};
      if (!exists $DataFieldValues{$SpecifiedDataField}) {
	warn "Warning: Ignoring file $SDFile: Data field value, $SpecifiedDataField, using  \"--CompoundID\" option in \"DataField\" \"--CompoundIDMode\" doesn't exist\n";
	next FILELIST;
      }
    }

    $AllDataFieldsRef = '';
    $CommonDataFieldsRef = '';
    if ($CollectDataFields) {
      my($CmpdCount);
      open SDFILE, "$SDFile" or die "Error: Couldn't open $SDFile: $! \n";
      ($CmpdCount, $AllDataFieldsRef, $CommonDataFieldsRef) = GetAllAndCommonCmpdDataHeaderLabels(\*SDFILE);
      close SDFILE;
    }

    # Setup output file names...
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($SDFile);

    $TextOutFileExt = "csv";
    if ($Options{outdelim} =~ /^tab$/i) {
      $TextOutFileExt = "tsv";
    }
    $SDOutFileExt = $FileExt;
    $FPOutFileExt = "fpf";

    if ($OptionsInfo{OutFileRoot} && (@SDFilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($OptionsInfo{OutFileRoot});
      if ($RootFileName && $RootFileExt) {
	$FileName = $RootFileName;
      }
      else {
	$FileName = $OptionsInfo{OutFileRoot};
      }
      $OutFileRoot = $FileName;
    }
    else {
      $OutFileRoot = "${FileName}PathLengthFP";
    }

    $NewSDFileName = "${OutFileRoot}.${SDOutFileExt}";
    $NewFPFileName = "${OutFileRoot}.${FPOutFileExt}";
    $NewTextFileName = "${OutFileRoot}.${TextOutFileExt}";

    if ($OptionsInfo{SDOutput}) {
      if ($SDFile =~ /$NewSDFileName/i) {
	warn "Warning: Ignoring input file $SDFile: Same output, $NewSDFileName, and input file names.\n";
	print "Specify a different name using \"-r --root\" option or use default name.\n";
	next FILELIST;
      }
    }

    if (!$OptionsInfo{OverwriteFiles}) {
      # Check SD, FP and text outout files...
      if ($OptionsInfo{SDOutput}) {
	if (-e $NewSDFileName) {
	  warn "Warning: Ignoring file $SDFile: The file $NewSDFileName already exists\n";
	  next FILELIST;
	}
      }
      if ($OptionsInfo{FPOutput}) {
	if (-e $NewFPFileName) {
	  warn "Warning: Ignoring file $SDFile: The file $NewFPFileName already exists\n";
	  next FILELIST;
	}
      }
      if ($OptionsInfo{TextOutput}) {
	if (-e $NewTextFileName) {
	  warn "Warning: Ignoring file $SDFile: The file $NewTextFileName already exists\n";
	  next FILELIST;
	}
      }
    }

    $SDFilesInfo{FileOkay}[$Index] = 1;

    $SDFilesInfo{OutFileRoot}[$Index] = $OutFileRoot;
    $SDFilesInfo{SDOutFileNames}[$Index] = $NewSDFileName;
    $SDFilesInfo{FPOutFileNames}[$Index] = $NewFPFileName;
    $SDFilesInfo{TextOutFileNames}[$Index] = $NewTextFileName;

    $SDFilesInfo{AllDataFieldsRef}[$Index] = $AllDataFieldsRef;
    $SDFilesInfo{CommonDataFieldsRef}[$Index] = $CommonDataFieldsRef;
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};
  $OptionsInfo{AromaticityModel} = $Options{aromaticitymodel};
  $OptionsInfo{PathMode} = $Options{pathmode};

  ProcessAtomIdentifierTypeOptions();

  $OptionsInfo{BitsOrder} = $Options{bitsorder};
  $OptionsInfo{BitStringFormat} = $Options{bitstringformat};

  $OptionsInfo{CompoundIDMode} = $Options{compoundidmode};
  $OptionsInfo{CompoundIDLabel} = $Options{compoundidlabel};
  $OptionsInfo{DataFieldsMode} = $Options{datafieldsmode};

  my(@SpecifiedDataFields);
  @SpecifiedDataFields = ();

  @{$OptionsInfo{SpecifiedDataFields}} = ();
  $OptionsInfo{CompoundID} = '';

  if ($Options{datafieldsmode} =~ /^CompoundID$/i) {
    if ($Options{compoundidmode} =~ /^DataField$/i) {
      if (!$Options{compoundid}) {
	die "Error: You must specify a value for \"--CompoundID\" option in \"DataField\" \"--CompoundIDMode\". \n";
      }
      $OptionsInfo{CompoundID} = $Options{compoundid};
    }
    elsif ($Options{compoundidmode} =~ /^(LabelPrefix|MolNameOrLabelPrefix)$/i) {
      $OptionsInfo{CompoundID} = $Options{compoundid} ? $Options{compoundid} : 'Cmpd';
    }
  }
  elsif ($Options{datafieldsmode} =~ /^Specify$/i) {
    if (!$Options{datafields}) {
      die "Error: You must specify a value for \"--DataFields\" option in \"Specify\" \"-d, --DataFieldsMode\". \n";
    }
    @SpecifiedDataFields = split /\,/, $Options{datafields};
    push @{$OptionsInfo{SpecifiedDataFields}}, @SpecifiedDataFields;
  }

  if ($Options{atomidentifiertype} !~ /^AtomicInvariantsAtomTypes$/i) {
    if ($Options{detectaromaticity} =~ /^No$/i) {
      die "Error: The value specified, $Options{detectaromaticity}, for option \"--DetectAromaticity\" is not valid. No value is only allowed during AtomicInvariantsAtomTypes value for \"-a, --AtomIdentifierType\" \n";
    }
  }
  $OptionsInfo{DetectAromaticity} = ($Options{detectaromaticity} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{Filter} = ($Options{filter} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{FingerprintsLabel} = $Options{fingerprintslabel} ? $Options{fingerprintslabel} : 'PathLengthFingerprints';

  my($Size, $MinSize, $MaxSize);
  $MinSize = 32;
  $MaxSize = 2**32;
  $Size = $Options{size};
  if (!(IsPositiveInteger($Size) && $Size >= $MinSize && $Size <= $MaxSize && IsNumberPowerOfNumber($Size, 2))) {
    die "Error: Invalid size value, $Size, for \"-s, --size\" option. Allowed values: power of 2, >= minimum size of $MinSize, and <= maximum size of $MaxSize.\n";
  }
  $OptionsInfo{Size} = $Size;

  $OptionsInfo{Fold} = ($Options{fold} =~ /^Yes$/i) ? 1 : 0;
  my($FoldedSize);
  $FoldedSize = $Options{foldedsize};
  if ($Options{fold} =~ /^Yes$/i) {
    if (!(IsPositiveInteger($FoldedSize) && $FoldedSize < $Size && IsNumberPowerOfNumber($FoldedSize, 2))) {
      die "Error: Invalid folded size value, $FoldedSize, for \"--FoldedSize\" option. Allowed values: power of 2, >= minimum size of $MinSize, and < size value of $Size.\n";
    }
  }
  $OptionsInfo{FoldedSize} = $FoldedSize;

  $OptionsInfo{IgnoreHydrogens} = ($Options{ignorehydrogens} =~ /^Yes$/i) ? 1 : 0;
  $OptionsInfo{KeepLargestComponent} = ($Options{keeplargestcomponent} =~ /^Yes$/i) ? 1 : 0;

  my($MinPathLength, $MaxPathLength);
  $MinPathLength = $Options{minpathlength};
  $MaxPathLength = $Options{maxpathlength};
  if (!IsPositiveInteger($MinPathLength)) {
    die "Error: Invalid path length value, $MinPathLength, for \"--MinPathLength\" option. Allowed values: > 0\n";
  }
  if (!IsPositiveInteger($MaxPathLength)) {
    die "Error: Invalid path length value, $MaxPathLength, for \"--MinPathLength\" option. Allowed values: > 0\n";
  }
  if ($MinPathLength >= $MaxPathLength) {
    die "Error: Invalid minimum and maximum path length values, $MinPathLength and $MaxPathLength, for \"--MinPathLength\"  and \"--MaxPathLength\"options. Allowed values: minimum path length value must be smaller than maximum path length value.\n";
  }
  $OptionsInfo{MinPathLength} = $MinPathLength;
  $OptionsInfo{MaxPathLength} = $MaxPathLength;

  my($NumOfBitsToSetPerPath);
  $NumOfBitsToSetPerPath = $Options{numofbitstosetperpath};
  if (!IsPositiveInteger($MaxPathLength)) {
    die "Error: Invalid  value, $NumOfBitsToSetPerPath, for \"-n, --NumOfBitsToSetPerPath\" option. Allowed values: > 0\n";
  }
  if ($NumOfBitsToSetPerPath >= $Size) {
    die "Error: Invalid  value, $NumOfBitsToSetPerPath, for \"-n, --NumOfBitsToSetPerPath\" option. Allowed values: It must be less than the size, $Size, of the fingerprint bit-string.\n";
  }
  $OptionsInfo{NumOfBitsToSetPerPath} = $NumOfBitsToSetPerPath;

  $OptionsInfo{Output} = $Options{output};
  $OptionsInfo{SDOutput} = ($Options{output} =~ /^(SD|All)$/i) ? 1 : 0;
  $OptionsInfo{FPOutput} = ($Options{output} =~ /^(FP|All)$/i) ? 1 : 0;
  $OptionsInfo{TextOutput} = ($Options{output} =~ /^(Text|All)$/i) ? 1 : 0;

  $OptionsInfo{OutDelim} = $Options{outdelim};
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? 1 : 0;
  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : 0;

  $OptionsInfo{UseBondSymbols} = ($Options{usebondsymbols} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{UsePerlCoreRandom} = ($Options{useperlcorerandom} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{UseUniquePaths} = ($Options{useuniquepaths} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{VectorStringFormat} = $Options{vectorstringformat};

  # Setup parameters used during generation of fingerprints by PathLengthFingerprints class...
  my($AllowRings, $AllowSharedBonds);
  $AllowRings = 1;
  $AllowSharedBonds = 1;
  MODE: {
    if ($Options{pathmode} =~ /^AtomPathsWithoutRings$/i) { $AllowSharedBonds = 0; $AllowRings = 0; last MODE;}
    if ($Options{pathmode} =~ /^AtomPathsWithRings$/i) { $AllowSharedBonds = 0; $AllowRings = 1; last MODE;}
    if ($Options{pathmode} =~ /^AllAtomPathsWithoutRings$/i) { $AllowSharedBonds = 1; $AllowRings = 0; last MODE;}
    if ($Options{pathmode} =~ /^AllAtomPathsWithRings$/i) { $AllowSharedBonds = 1; $AllowRings = 1; last MODE;}
    die "Error: ProcessOptions: mode value, $Options{pathmode}, is not supported.\n";
  }
  $OptionsInfo{AllowRings} = $AllowRings;
  $OptionsInfo{AllowSharedBonds} = $AllowSharedBonds;
}

# Process atom identifier type and related options...
#
sub ProcessAtomIdentifierTypeOptions {

  $OptionsInfo{AtomIdentifierType} = $Options{atomidentifiertype};

  if ($Options{atomidentifiertype} =~ /^AtomicInvariantsAtomTypes$/i) {
    ProcessAtomicInvariantsToUseOption();
  }
  elsif ($Options{atomidentifiertype} =~ /^FunctionalClassAtomTypes$/i) {
    ProcessFunctionalClassesToUse();
  }
  elsif ($OptionsInfo{AtomIdentifierType} =~ /^(DREIDINGAtomTypes|EStateAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
    # Nothing to do for now...
  }
  else {
    die "Error: The value specified, $Options{atomidentifiertype}, for option \"-a, --AtomIdentifierType\" is not valid. Supported atom identifier types in current release of MayaChemTools: AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes\n";
  }
}

# Process specified atomic invariants to use...
#
sub ProcessAtomicInvariantsToUseOption {
  my($AtomicInvariant, $AtomSymbolSpecified, @AtomicInvariantsWords);

  @{$OptionsInfo{AtomicInvariantsToUse}} = ();
  if (IsEmpty($Options{atomicinvariantstouse})) {
    die "Error: Atomic invariants value specified using \"--AtomicInvariantsToUse\" option is empty\n";
  }
  $AtomSymbolSpecified = 0;
  @AtomicInvariantsWords = split /\,/, $Options{atomicinvariantstouse};
  for $AtomicInvariant (@AtomicInvariantsWords) {
    if (!AtomTypes::AtomicInvariantsAtomTypes::IsAtomicInvariantAvailable($AtomicInvariant)) {
      die "Error: Atomic invariant specified, $AtomicInvariant, using \"--AtomicInvariantsToUse\" option is not valid...\n ";
    }
    if ($AtomicInvariant =~ /^(AS|AtomSymbol)$/i) {
      $AtomSymbolSpecified = 1;
    }
    push @{$OptionsInfo{AtomicInvariantsToUse}}, $AtomicInvariant;
  }
  if (!$AtomSymbolSpecified) {
    die "Error: Atomic invariant, AS or AtomSymbol, must be specified as using \"--AtomicInvariantsToUse\" option...\n ";
  }
}

# Process specified functional classes invariants to use...
#
sub ProcessFunctionalClassesToUse {
  my($FunctionalClass, @FunctionalClassesToUseWords);

  @{$OptionsInfo{FunctionalClassesToUse}} = ();
  if (IsEmpty($Options{functionalclassestouse})) {
    die "Error: Functional classes value specified using \"--FunctionalClassesToUse\" option is empty\n";
  }
  @FunctionalClassesToUseWords = split /\,/, $Options{functionalclassestouse};
  for $FunctionalClass (@FunctionalClassesToUseWords) {
    if (!AtomTypes::FunctionalClassAtomTypes::IsFunctionalClassAvailable($FunctionalClass)) {
      die "Error: Functional class specified, $FunctionalClass, using \"--FunctionalClassesToUse\" option is not valid...\n ";
    }
    push @{$OptionsInfo{FunctionalClassesToUse}}, $FunctionalClass;
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{aromaticitymodel} = 'MayaChemToolsAromaticityModel';

  $Options{atomidentifiertype} = 'AtomicInvariantsAtomTypes';
  $Options{atomicinvariantstouse} = 'AS';

  $Options{functionalclassestouse} = 'HBD,HBA,PI,NI,Ar,Hal';

  $Options{bitsorder} = 'Ascending';
  $Options{bitstringformat} = 'HexadecimalString';

  $Options{compoundidmode} = 'LabelPrefix';
  $Options{compoundidlabel} = 'CompoundID';
  $Options{datafieldsmode} = 'CompoundID';
  $Options{detectaromaticity} = 'Yes';

  $Options{filter} = 'Yes';

  $Options{fold} = 'No';
  $Options{foldedsize} = 256;

  $Options{ignorehydrogens} = 'Yes';
  $Options{keeplargestcomponent} = 'Yes';

  $Options{mode} = 'PathLengthBits';
  $Options{pathmode} = 'AllAtomPathsWithRings';

  $Options{minpathlength} = 1;
  $Options{maxpathlength} = 8;

  $Options{numofbitstosetperpath} = 1;

  $Options{output} = 'text';
  $Options{outdelim} = 'comma';
  $Options{quote} = 'yes';

  $Options{size} = 1024;

  $Options{usebondsymbols} = 'yes';
  $Options{useperlcorerandom} = 'yes';
  $Options{useuniquepaths} = 'yes';

  $Options{vectorstringformat} = 'IDsAndValuesString';

  if (!GetOptions(\%Options, "aromaticitymodel=s", "atomidentifiertype|a=s", "atomicinvariantstouse=s", "functionalclassestouse=s", "bitsorder=s", "bitstringformat|b=s", "compoundid=s", "compoundidlabel=s", "compoundidmode=s", "datafields=s", "datafieldsmode|d=s", "detectaromaticity=s", "filter|f=s", "fingerprintslabel=s", "fold=s", "foldedsize=i", "help|h", "ignorehydrogens|i=s", "keeplargestcomponent|k=s", "mode|m=s", "minpathlength=i", "maxpathlength=i", "numofbitstosetperpath|n=i", "outdelim=s", "output=s", "overwrite|o", "pathmode|p=s", "quote|q=s", "root|r=s", "size|s=i", "usebondsymbols|u=s", "useperlcorerandom=s", "useuniquepaths=s", "vectorstringformat|v=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if (!Molecule::IsSupportedAromaticityModel($Options{aromaticitymodel})) {
    my(@SupportedModels) = Molecule::GetSupportedAromaticityModels();
    die "Error: The value specified, $Options{aromaticitymodel}, for option \"--AromaticityModel\" is not valid. Supported aromaticity models in current release of MayaChemTools: @SupportedModels\n";
  }
  if ($Options{atomidentifiertype} !~ /^(AtomicInvariantsAtomTypes|DREIDINGAtomTypes|EStateAtomTypes|FunctionalClassAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
    die "Error: The value specified, $Options{atomidentifiertype}, for option \"-a, --AtomIdentifierType\" is not valid. Supported atom identifier types in current release of MayaChemTools: AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes\n";
  }
  if ($Options{bitsorder} !~ /^(Ascending|Descending)$/i) {
    die "Error: The value specified, $Options{bitsorder}, for option \"--BitsOrder\" is not valid. Allowed values: Ascending or Descending\n";
  }
  if ($Options{bitstringformat} !~ /^(BinaryString|HexadecimalString)$/i) {
    die "Error: The value specified, $Options{bitstringformat}, for option \"-b, --bitstringformat\" is not valid. Allowed values: BinaryString or HexadecimalString\n";
  }
  if ($Options{compoundidmode} !~ /^(DataField|MolName|LabelPrefix|MolNameOrLabelPrefix)$/i) {
    die "Error: The value specified, $Options{compoundidmode}, for option \"--CompoundIDMode\" is not valid. Allowed values: DataField, MolName, LabelPrefix or MolNameOrLabelPrefix\n";
  }
  if ($Options{datafieldsmode} !~ /^(All|Common|Specify|CompoundID)$/i) {
    die "Error: The value specified, $Options{datafieldsmode}, for option \"-d, --DataFieldsMode\" is not valid. Allowed values: All, Common, Specify or CompoundID\n";
  }
  if ($Options{detectaromaticity} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{detectaromaticity}, for option \"--DetectAromaticity\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{filter} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{filter}, for option \"-f, --Filter\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{fold} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{fold}, for option \"--fold\" is not valid. Allowed values: Yes or No\n";
  }
  if (!IsPositiveInteger($Options{foldedsize})) {
    die "Error: The value specified, $Options{foldedsize}, for option \"--FoldedSize\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{ignorehydrogens} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{ignorehydrogens}, for option \"-i, --IgnoreHydrogens\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{keeplargestcomponent} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{keeplargestcomponent}, for option \"-k, --KeepLargestComponent\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{mode} !~ /^(PathLengthBits|PathLengthCount)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: PathLengthBits or PathLengthCount\n";
  }
  if (!IsPositiveInteger($Options{minpathlength})) {
    die "Error: The value specified, $Options{minpathlength}, for option \"--MinPathLength\" is not valid. Allowed values: > 0 \n";
  }
  if (!IsPositiveInteger($Options{numofbitstosetperpath})) {
    die "Error: The value specified, $Options{NumOfBitsToSetPerPath}, for option \"--NumOfBitsToSetPerPath\" is not valid. Allowed values: > 0 \n";
  }
  if (!IsPositiveInteger($Options{maxpathlength})) {
    die "Error: The value specified, $Options{maxpathlength}, for option \"--MaxPathLength\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{output} !~ /^(SD|FP|text|all)$/i) {
    die "Error: The value specified, $Options{output}, for option \"--output\" is not valid. Allowed values: SD, FP, text, or all\n";
  }
  if ($Options{outdelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{outdelim}, for option \"--outdelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if ($Options{pathmode} !~ /^(AtomPathsWithoutRings|AtomPathsWithRings|AllAtomPathsWithoutRings|AllAtomPathsWithRings)$/i) {
    die "Error: The value specified, $Options{pathmode}, for option \"-m, --PathMode\" is not valid. Allowed values: AtomPathsWithoutRings, AtomPathsWithRings, AllAtomPathsWithoutRings or AllAtomPathsWithRings\n";
  }
  if ($Options{quote} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{outdelim} =~ /semicolon/i && $Options{quote} =~ /^No$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not allowed with, semicolon value of \"--outdelim\" option: Fingerprints string use semicolon as delimiter for various data fields and must be quoted.\n";
  }

  if (!IsPositiveInteger($Options{size})) {
    die "Error: The value specified, $Options{size}, for option \"-s, --size\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{usebondsymbols} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{usebondsymbols}, for option \"-u, --UseBondSymbols\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{useperlcorerandom} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{useperlcorerandom}, for option \"--UsePerlCoreRandom\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{useuniquepaths} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{useuniquepaths}, for option \"--UseUniquePaths\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{vectorstringformat} !~ /^(IDsAndValuesString|IDsAndValuesPairsString|ValuesAndIDsString|ValuesAndIDsPairsString)$/i) {
    die "Error: The value specified, $Options{vectorstringformat}, for option \"-v, --VectorStringFormat\" is not valid. Allowed values: IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString or ValuesAndIDsPairsString\n";
  }
}

__END__

=head1 NAME

PathLengthFingerprints.pl - Generate atom path length based fingerprints for SD files

=head1 SYNOPSIS

PathLengthFingerprints.pl SDFile(s)...

PathLengthFingerprints.pl [B<--AromaticityModel> I<AromaticityModelType>]
[B<-a, --AtomIdentifierType> I<AtomicInvariantsAtomTypes>]
[B<--AtomicInvariantsToUse> I<"AtomicInvariant1,AtomicInvariant2...">]
[B<--FunctionalClassesToUse> I<"FunctionalClass1,FunctionalClass2...">]
[B<--BitsOrder> I<Ascending | Descending>] [B<-b, --BitStringFormat> I<BinaryString | HexadecimalString>]
[B<--CompoundID> I<DataFieldName or LabelPrefixString>] [B<--CompoundIDLabel> I<text>]
[B<--CompoundIDMode> I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>]
[B<--DataFields> I<"FieldLabel1,FieldLabel2,... ">] [B<-d, --DataFieldsMode> I<All | Common | Specify | CompoundID>]
[B<--DetectAromaticity> I<Yes | No>]  [B<-f, --Filter> I<Yes | No>] [B<--FingerprintsLabel> I<text>]
[B<--fold> I<Yes | No>] [B<--FoldedSize> I<number>] [B<-h, --help>]
[B<-i, --IgnoreHydrogens> I<Yes | No>] [B<-k, --KeepLargestComponent> I<Yes | No>]
[B<-m, --mode> I<PathLengthBits | PathLengthCount>]
[B<--MinPathLength> I<number>] [B<--MaxPathLength> I<number>] [B<-n, --NumOfBitsToSetPerPath> I<number>]
[B<--OutDelim> I<comma | tab | semicolon>]
[B<--output> I<SD | FP | text | all>] [B<-q, --quote> I<Yes | No>] [B<-r, --root> I<RootName>]
[B<-p, --PathMode> I<AtomPathsWithoutRings | AtomPathsWithRings | AllAtomPathsWithoutRings | AllAtomPathsWithRings>]
[B<-s, --size> I<number>] [B<-u, --UseBondSymbols> I<Yes | No>] [B<--UsePerlCoreRandom> I<Yes | No>]
[B<--UseUniquePaths> I<Yes | No>] [B<-q, --quote> I<Yes | No>] [B<-r, --root> I<RootName>]
[B<-v, --VectorStringFormat> I<IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>]
[B<-w, --WorkingDir> dirname] SDFile(s)...

=head1 DESCRIPTION

Generate atom path length fingerprints for I<SDFile(s)> and create appropriate SD, FP or
CSV/TSV text file(s) containing fingerprints bit-vector or vector strings corresponding to
molecular fingerprints.

Multiple SDFile names are separated by spaces. The valid file extensions are I<.sdf>
and I<.sd>. All other file names are ignored. All the SD files in a current directory
can be specified either by I<*.sdf> or the current directory name.

The current release of MayaChemTools supports generation of path length fingerprints
corresponding to following B<-a, --AtomIdentifierTypes>:

    AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
    FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
    SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes

Based on the values specified for B<-p, --PathMode>, B<--MinPathLength> and B<--MaxPathLength>,
all appropriate atom paths are generated for each atom in the molecule and collected in a list and
the list is filtered to remove any structurally duplicate paths as indicated by the value of
B<--UseUniquePaths> option.

For each atom path in the filtered atom paths list, an atom path string is created using value of
B<-a, --AtomIdentifierType> and specified values to use for a particular atom identifier type.
Value of B<-u, --UseBondSymbols> controls whether bond order symbols are used during generation
of atom path string. For each atom path, only lexicographically smaller atom path strings are kept.

For I<PathLengthBits> value of B<-m, --mode> option, each atom path is hashed to a 32 bit unsigned
integer key using B<TextUtil::HashCode> function. Using the hash key as a seed for a random number
generator, a random integer value between 0 and B<--Size> is used to set corresponding bits
in the fingerprint bit-vector string. Value of B<--NumOfBitsToSetPerPath> option controls the number
of time a random number is generated to set corresponding bits.

For I< PathLengthCount> value of B<-m, --mode> option, the number of times an atom path appears
is tracked and a fingerprints count-string corresponding to count of atom paths is generated.

Example of I<SD> file containing path length fingerprints string data:

    ... ...
    ... ...
    $$$$
    ... ...
    ... ...
    ... ...
    41 44  0  0  0  0  0  0  0  0999 V2000
     -3.3652    1.4499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    ... ...
    2  3  1  0  0  0  0
    ... ...
    M  END
    >  <CmpdID>
    Cmpd1

    >  <PathLengthFingerprints>
    FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes:MinLengt
    h1:MaxLength8;1024;HexadecimalString;Ascending;9c8460989ec8a49913991a66
    03130b0a19e8051c89184414953800cc2151082844a201042800130860308e8204d4028
    00831048940e44281c00060449a5000ac80c894114e006321264401600846c050164462
    08190410805000304a10205b0100e04c0038ba0fad0209c0ca8b1200012268b61c0026a
    aa0660a11014a011d46

    $$$$
    ... ...
    ... ...

Example of I<FP> file containing path length fingerprints string data:

    #
    # Package = MayaChemTools 7.4
    # ReleaseDate = Oct 21, 2010
    #
    # TimeStamp =  Mon Mar 7 15:14:01 2011
    #
    # FingerprintsStringType = FingerprintsBitVector
    #
    # Description = PathLengthBits:AtomicInvariantsAtomTypes:MinLength1:...
    # Size = 1024
    # BitStringFormat = HexadecimalString
    # BitsOrder = Ascending
    #
    Cmpd1 9c8460989ec8a49913991a6603130b0a19e8051c89184414953800cc21510...
    Cmpd2 000000249400840040100042011001001980410c000000001010088001120...
    ... ...
    ... ..

Example of CSV I<Text> file containing pathlength fingerprints string data:

    "CompoundID","PathLengthFingerprints"
    "Cmpd1","FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes
    :MinLength1:MaxLength8;1024;HexadecimalString;Ascending;9c8460989ec8a4
    9913991a6603130b0a19e8051c89184414953800cc2151082844a20104280013086030
    8e8204d402800831048940e44281c00060449a5000ac80c894114e006321264401..."
    ... ...
    ... ...

The current release of MayaChemTools generates the following types of path length
fingerprints bit-vector and vector strings:

    FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes:MinLeng
    th1:MaxLength8;1024;BinaryString;Ascending;001000010011010101011000110
    0100010101011000101001011100110001000010001001101000001001001001001000
    0010110100000111001001000001001010100100100000000011000000101001011100
    0010000001000101010100000100111100110111011011011000000010110111001101
    0101100011000000010001000011000010100011101100001000001000100000000...

    FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes:MinLeng
    th1:MaxLength8;1024;HexadecimalString;Ascending;48caa1315d82d91122b029
    42861c9409a4208182d12015509767bd0867653604481a8b1288000056090583603078
    9cedae54e26596889ab121309800900490515224208421502120a0dd9200509723ae89
    00024181b86c0122821d4e4880c38620dab280824b455404009f082003d52c212b4e6d
    6ea05280140069c780290c43

    FingerprintsVector;PathLengthCount:AtomicInvariantsAtomTypes:MinLength
    1:MaxLength8;432;NumericalValues;IDsAndValuesPairsString;C.X1.BO1.H3 2
    C.X2.BO2.H2 4 C.X2.BO3.H1 14 C.X3.BO3.H1 3 C.X3.BO4 10 F.X1.BO1 1 N.X
    2.BO2.H1 1 N.X3.BO3 1 O.X1.BO1.H1 3 O.X1.BO2 2 C.X1.BO1.H3C.X3.BO3.H1
    2 C.X2.BO2.H2C.X2.BO2.H2 1 C.X2.BO2.H2C.X3.BO3.H1 4 C.X2.BO2.H2C.X3.BO
    4 1 C.X2.BO2.H2N.X3.BO3 1 C.X2.BO3.H1:C.X2.BO3.H1 10 C.X2.BO3.H1:C....

    FingerprintsVector;PathLengthCount:DREIDINGAtomTypes:MinLength1:MaxLen
    gth8;410;NumericalValues;IDsAndValuesPairsString;C_2 2 C_3 9 C_R 22 F_
    1 N_3 1 N_R 1 O_2 2 O_3 3 C_2=O_2 2 C_2C_3 1 C_2C_R 1 C_2N_3 1 C_2O_3
    1 C_3C_3 7 C_3C_R 1 C_3N_R 1 C_3O_3 2 C_R:C_R 21 C_R:N_R 2 C_RC_R 2 C
    _RF_ 1 C_RN_3 1 C_2C_3C_3 1 C_2C_R:C_R 2 C_2N_3C_R 1 C_3C_2=O_2 1 C_3C
    _2O_3 1 C_3C_3C_3 5 C_3C_3C_R 2 C_3C_3N_R 1 C_3C_3O_3 4 C_3C_R:C_R ...

    FingerprintsVector;PathLengthCount:EStateAtomTypes:MinLength1:MaxLengt
    h8;454;NumericalValues;IDsAndValuesPairsString;aaCH 14 aasC 8 aasN 1 d
    O 2 dssC 2 sCH3 2 sF 1 sOH 3 ssCH2 4 ssNH 1 sssCH 3 aaCH:aaCH 10 aaCH:
    aasC 8 aasC:aasC 3 aasC:aasN 2 aasCaasC 2 aasCdssC 1 aasCsF 1 aasCssNH
    1 aasCsssCH 1 aasNssCH2 1 dO=dssC 2 dssCsOH 1 dssCssCH2 1 dssCssNH 1
    sCH3sssCH 2 sOHsssCH 2 ssCH2ssCH2 1 ssCH2sssCH 4 aaCH:aaCH:aaCH 6 a...

    FingerprintsVector;PathLengthCount:FunctionalClassAtomTypes:MinLength1
    :MaxLength8;404;NumericalValues;IDsAndValuesPairsString;Ar 22 Ar.HBA 1
    HBA 2 HBA.HBD 3 HBD 1 Hal 1 NI 1 None 10 Ar.HBA:Ar 2 Ar.HBANone 1 Ar:
    Ar 21 ArAr 2 ArHBD 1 ArHal 1 ArNone 2 HBA.HBDNI 1 HBA.HBDNone 2 HBA=NI
    1 HBA=None 1 HBDNone 1 NINone 1 NoneNone 7 Ar.HBA:Ar:Ar 2 Ar.HBA:ArAr
    1 Ar.HBA:ArNone 1 Ar.HBANoneNone 1 Ar:Ar.HBA:Ar 1 Ar:Ar.HBANone 2 ...

    FingerprintsVector;PathLengthCount:MMFF94AtomTypes:MinLength1:MaxLengt
    h8;463;NumericalValues;IDsAndValuesPairsString;C5A 2 C5B 2 C=ON 1 CB 1
    8 COO 1 CR 9 F 1 N5 1 NC=O 1 O=CN 1 O=CO 1 OC=O 1 OR 2 C5A:C5B 2 C5A:N
    5 2 C5ACB 1 C5ACR 1 C5B:C5B 1 C5BC=ON 1 C5BCB 1 C=ON=O=CN 1 C=ONNC=O 1
    CB:CB 18 CBF 1 CBNC=O 1 COO=O=CO 1 COOCR 1 COOOC=O 1 CRCR 7 CRN5 1 CR
    OR 2 C5A:C5B:C5B 2 C5A:C5BC=ON 1 C5A:C5BCB 1 C5A:N5:C5A 1 C5A:N5CR ...

    FingerprintsVector;PathLengthCount:SLogPAtomTypes:MinLength1:MaxLength
    8;518;NumericalValues;IDsAndValuesPairsString;C1 5 C10 1 C11 1 C14 1 C
    18 14 C20 4 C21 2 C22 1 C5 2 CS 2 F 1 N11 1 N4 1 O10 1 O2 3 O9 1 C10C1
    1 C10N11 1 C11C1 2 C11C21 1 C14:C18 2 C14F 1 C18:C18 10 C18:C20 4 C18
    :C22 2 C1C5 1 C1CS 4 C20:C20 1 C20:C21 1 C20:N11 1 C20C20 2 C21:C21 1
    C21:N11 1 C21C5 1 C22N4 1 C5=O10 1 C5=O9 1 C5N4 1 C5O2 1 CSO2 2 C10...

    FingerprintsVector;PathLengthCount:SYBYLAtomTypes:MinLength1:MaxLength
    8;412;NumericalValues;IDsAndValuesPairsString;C.2 2 C.3 9 C.ar 22 F 1
    N.am 1 N.ar 1 O.2 1 O.3 2 O.co2 2 C.2=O.2 1 C.2=O.co2 1 C.2C.3 1 C.2C.
    ar 1 C.2N.am 1 C.2O.co2 1 C.3C.3 7 C.3C.ar 1 C.3N.ar 1 C.3O.3 2 C.ar:C
    .ar 21 C.ar:N.ar 2 C.arC.ar 2 C.arF 1 C.arN.am 1 C.2C.3C.3 1 C.2C.ar:C
    .ar 2 C.2N.amC.ar 1 C.3C.2=O.co2 1 C.3C.2O.co2 1 C.3C.3C.3 5 C.3C.3...

    FingerprintsVector;PathLengthCount:TPSAAtomTypes:MinLength1:MaxLength8
    ;331;NumericalValues;IDsAndValuesPairsString;N21 1 N7 1 None 34 O3 2 O
    4 3 N21:None 2 N21None 1 N7None 2 None:None 21 None=O3 2 NoneNone 13 N
    oneO4 3 N21:None:None 2 N21:NoneNone 2 N21NoneNone 1 N7None:None 2 N7N
    one=O3 1 N7NoneNone 1 None:N21:None 1 None:N21None 2 None:None:None 20
    None:NoneNone 12 NoneN7None 1 NoneNone=O3 2 NoneNoneNone 8 NoneNon...

    FingerprintsVector;PathLengthCount:UFFAtomTypes:MinLength1:MaxLength8;
    410;NumericalValues;IDsAndValuesPairsString;C_2 2 C_3 9 C_R 22 F_ 1 N_
    3 1 N_R 1 O_2 2 O_3 3 C_2=O_2 2 C_2C_3 1 C_2C_R 1 C_2N_3 1 C_2O_3 1 C_
    3C_3 7 C_3C_R 1 C_3N_R 1 C_3O_3 2 C_R:C_R 21 C_R:N_R 2 C_RC_R 2 C_RF_
    1 C_RN_3 1 C_2C_3C_3 1 C_2C_R:C_R 2 C_2N_3C_R 1 C_3C_2=O_2 1 C_3C_2O_3
    1 C_3C_3C_3 5 C_3C_3C_R 2 C_3C_3N_R 1 C_3C_3O_3 4 C_3C_R:C_R 1 C_3...

=head1 OPTIONS

=over 4

=item B<--AromaticityModel> I<MDLAromaticityModel | TriposAromaticityModel | MMFFAromaticityModel | ChemAxonBasicAromaticityModel | ChemAxonGeneralAromaticityModel | DaylightAromaticityModel | MayaChemToolsAromaticityModel>

Specify aromaticity model to use during detection of aromaticity. Possible values in the current
release are: I<MDLAromaticityModel, TriposAromaticityModel, MMFFAromaticityModel,
ChemAxonBasicAromaticityModel, ChemAxonGeneralAromaticityModel, DaylightAromaticityModel
or MayaChemToolsAromaticityModel>. Default value: I<MayaChemToolsAromaticityModel>.

The supported aromaticity model names along with model specific control parameters
are defined in B<AromaticityModelsData.csv>, which is distributed with the current release
and is available under B<lib/data> directory. B<Molecule.pm> module retrieves data from
this file during class instantiation and makes it available to method B<DetectAromaticity>
for detecting aromaticity corresponding to a specific model.

This option is ignored during I<No> value of B<--DetectAromaticity> option.

=item B<-a, --AtomIdentifierType> I<AtomicInvariantsAtomTypes | DREIDINGAtomTypes | EStateAtomTypes | FunctionalClassAtomTypes | MMFF94AtomTypes | SLogPAtomTypes | SYBYLAtomTypes | TPSAAtomTypes | UFFAtomTypes>

Specify atom identifier type to use for assignment of atom types to hydrogen and/or
non-hydrogen atoms during calculation of atom types fingerprints. Possible values in the
current release are: I<AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes,
TPSAAtomTypes, UFFAtomTypes>. Default value: I<AtomicInvariantsAtomTypes>.

=item B<-a, --AtomIdentifierType> I<AtomicInvariantsAtomTypes | DREIDINGAtomTypes | EStateAtomTypes | FunctionalClassAtomTypes | MMFF94AtomTypes | SLogPAtomTypes | SYBYLAtomTypes | TPSAAtomTypes | UFFAtomTypes>

Specify atom identifier type to use during generation of atom path strings
corresponding to path length fingerprints. Possible values in the current release are:
I<AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes,
TPSAAtomTypes, UFFAtomTypes>. Default value: I<AtomicInvariantsAtomTypes>.

=item B<--AtomicInvariantsToUse> I<"AtomicInvariant1,AtomicInvariant2...">

This value is used during I<AtomicInvariantsAtomTypes> value of B<a, --AtomIdentifierType>
option. It's a list of comma separated valid atomic invariant atom types.

Possible values for atomic invariants are: I<AS, X, BO,  LBO, SB, DB, TB,
H, Ar, RA, FC, MN, SM>. Default value: I<AS>.

The atomic invariants abbreviations correspond to:

    AS = Atom symbol corresponding to element symbol

    X<n>   = Number of non-hydrogen atom neighbors or heavy atoms
    BO<n> = Sum of bond orders to non-hydrogen atom neighbors or heavy atoms
    LBO<n> = Largest bond order of non-hydrogen atom neighbors or heavy atoms
    SB<n> = Number of single bonds to non-hydrogen atom neighbors or heavy atoms
    DB<n> = Number of double bonds to non-hydrogen atom neighbors or heavy atoms
    TB<n> = Number of triple bonds to non-hydrogen atom neighbors or heavy atoms
    H<n>   = Number of implicit and explicit hydrogens for atom
    Ar     = Aromatic annotation indicating whether atom is aromatic
    RA     = Ring atom annotation indicating whether atom is a ring
    FC<+n/-n> = Formal charge assigned to atom
    MN<n> = Mass number indicating isotope other than most abundant isotope
    SM<n> = Spin multiplicity of atom. Possible values: 1 (singlet), 2 (doublet) or
            3 (triplet)

Atom type generated by AtomTypes::AtomicInvariantsAtomTypes class corresponds to:

    AS.X<n>.BO<n>.LBO<n>.<SB><n>.<DB><n>.<TB><n>.H<n>.Ar.RA.FC<+n/-n>.MN<n>.SM<n>

Except for AS which is a required atomic invariant in atom types, all other atomic invariants are
optional. Atom type specification doesn't include atomic invariants with zero or undefined values.

In addition to usage of abbreviations for specifying atomic invariants, the following descriptive words
are also allowed:

    X : NumOfNonHydrogenAtomNeighbors or NumOfHeavyAtomNeighbors
    BO : SumOfBondOrdersToNonHydrogenAtoms or SumOfBondOrdersToHeavyAtoms
    LBO : LargestBondOrderToNonHydrogenAtoms or LargestBondOrderToHeavyAtoms
    SB :  NumOfSingleBondsToNonHydrogenAtoms or NumOfSingleBondsToHeavyAtoms
    DB : NumOfDoubleBondsToNonHydrogenAtoms or NumOfDoubleBondsToHeavyAtoms
    TB : NumOfTripleBondsToNonHydrogenAtoms or NumOfTripleBondsToHeavyAtoms
    H :  NumOfImplicitAndExplicitHydrogens
    Ar : Aromatic
    RA : RingAtom
    FC : FormalCharge
    MN : MassNumber
    SM : SpinMultiplicity

Examples:

B<Benzene>: Using value of I<AS> for B<--AtomicInvariantsToUse>, I<Yes> for B<UseBondSymbols>,
and I< AllAtomPathsWithRings> for B<-p, --PathMode>, atom path strings generated are:

    C C:C C:C:C C:C:C:C C:C:C:C:C C:C:C:C:C:C C:C:C:C:C:C:C

And using I<AS,X,BO> for B<--AtomicInvariantsToUse> generates following atom path
strings:

    C.X2.BO3 C.X2.BO3:C.X2.BO3 C.X2.BO3:C.X2.BO3:C.X2.BO3
    C.X2.BO3:C.X2.BO3:C.X2.BO3:C.X2.BO3
    C.X2.BO3:C.X2.BO3:C.X2.BO3:C.X2.BO3:C.X2.BO3
    C.X2.BO3:C.X2.BO3:C.X2.BO3:C.X2.BO3:C.X2.BO3:C.X2.BO3
    C.X2.BO3:C.X2.BO3:C.X2.BO3:C.X2.BO3:C.X2.BO3:C.X2.BO3:C.X2.BO3

B<Urea>: Using value of I<AS> for B<--AtomicInvariantsToUse>, I<Yes> for B<UseBondSymbols>,
and I< AllAtomPathsWithRings> for B<-p, --PathMode>, atom path strings are:

    C N O C=O CN NC=O NCN

And using I<AS,X,BO> for B<--AtomicInvariantsToUse> generates following atom path
strings:

    C.X3.BO4 N.X1.BO1 O.X1.BO2 C.X3.BO4=O.X1.BO2
    C.X3.BO4N.X1.BO1 N.X1.BO1C.X3.BO4=O.X1.BO2
    N.X1.BO1C.X3.BO4N.X1.BO1

=item B<--FunctionalClassesToUse> I<"FunctionalClass1,FunctionalClass2...">

This value is used during I<FunctionalClassAtomTypes> value of B<a, --AtomIdentifierType>
option. It's a list of comma separated valid functional classes.

Possible values for atom functional classes are: I<Ar, CA, H, HBA, HBD, Hal, NI, PI, RA>.
Default value [ Ref 24 ]: I<HBD,HBA,PI,NI,Ar,Hal>.

The functional class abbreviations correspond to:

    HBD: HydrogenBondDonor
    HBA: HydrogenBondAcceptor
    PI :  PositivelyIonizable
    NI : NegativelyIonizable
    Ar : Aromatic
    Hal : Halogen
    H : Hydrophobic
    RA : RingAtom
    CA : ChainAtom

 Functional class atom type specification for an atom corresponds to:

    Ar.CA.H.HBA.HBD.Hal.NI.PI.RA

I<AtomTypes::FunctionalClassAtomTypes> module is used to assign functional class atom
types. It uses following definitions [ Ref 60-61, Ref 65-66 ]:

    HydrogenBondDonor: NH, NH2, OH
    HydrogenBondAcceptor: N[!H], O
    PositivelyIonizable: +, NH2
    NegativelyIonizable: -, C(=O)OH, S(=O)OH, P(=O)OH

=item B<--BitsOrder> I<Ascending | Descending>

Bits order to use during generation of fingerprints bit-vector string for I<PathLengthBits> value of
B<-m, --mode> option. Possible values: I<Ascending, Descending>. Default: I<Ascending>.

I<Ascending> bit order which corresponds to first bit in each byte as the lowest bit as
opposed to the highest bit.

Internally, bits are stored in I<Ascending> order using Perl vec function. Regardless
of machine order, big-endian or little-endian, vec function always considers first
string byte as the lowest byte and first bit within each byte as the lowest bit.

=item B<-b, --BitStringFormat> I<BinaryString | HexadecimalString>

Format of fingerprints bit-vector string data in output SD, FP or CSV/TSV text file(s) specified by
B<--output> used during I<PathLengthBits> value of B<-m, --mode> option. Possible
values: I<BinaryString, HexadecimalString>. Default value: I<HexadecimalString>.

I<BinaryString> corresponds to an ASCII string containing 1s and 0s. I<HexadecimalString>
contains bit values in ASCII hexadecimal format.

Examples:

    FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes:MinLeng
    th1:MaxLength8;1024;BinaryString;Ascending;001000010011010101011000110
    0100010101011000101001011100110001000010001001101000001001001001001000
    0010110100000111001001000001001010100100100000000011000000101001011100
    0010000001000101010100000100111100110111011011011000000010110111001101
    0101100011000000010001000011000010100011101100001000001000100000000...

    FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes:MinLeng
    th1:MaxLength8;1024;HexadecimalString;Ascending;48caa1315d82d91122b029
    42861c9409a4208182d12015509767bd0867653604481a8b1288000056090583603078
    9cedae54e26596889ab121309800900490515224208421502120a0dd9200509723ae89
    00024181b86c0122821d4e4880c38620dab280824b455404009f082003d52c212b4e6d
    6ea05280140069c780290c43

=item B<--CompoundID> I<DataFieldName or LabelPrefixString>

This value is B<--CompoundIDMode> specific and indicates how compound ID is generated.

For I<DataField> value of B<--CompoundIDMode> option, it corresponds to datafield label name
whose value is used as compound ID; otherwise, it's a prefix string used for generating compound
IDs like LabelPrefixString<Number>. Default value, I<Cmpd>, generates compound IDs which
look like Cmpd<Number>.

Examples for I<DataField> value of B<--CompoundIDMode>:

    MolID
    ExtReg

Examples for I<LabelPrefix> or I<MolNameOrLabelPrefix> value of B<--CompoundIDMode>:

    Compound

The value specified above generates compound IDs which correspond to Compound<Number>
instead of default value of Cmpd<Number>.

=item B<--CompoundIDLabel> I<text>

Specify compound ID column label for FP or CSV/TSV text file(s) used during I<CompoundID> value
of B<--DataFieldsMode> option. Default: I<CompoundID>.

=item B<--CompoundIDMode> I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>

Specify how to generate compound IDs and write to FP or CSV/TSV text file(s) along with generated
fingerprints for I<FP | text | all> values of B<--output> option: use a I<SDFile(s)> datafield value;
use molname line from I<SDFile(s)>; generate a sequential ID with specific prefix; use combination
of both MolName and LabelPrefix with usage of LabelPrefix values for empty molname lines.

Possible values: I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>.
Default: I<LabelPrefix>.

For I<MolNameAndLabelPrefix> value of B<--CompoundIDMode>, molname line in I<SDFile(s)> takes
precedence over sequential compound IDs generated using I<LabelPrefix> and only empty molname
values are replaced with sequential compound IDs.

This is only used for I<CompoundID> value of B<--DataFieldsMode> option.

=item B<--DataFields> I<"FieldLabel1,FieldLabel2,... ">

Comma delimited list of I<SDFiles(s)> data fields to extract and write to CSV/TSV text file(s) along
with generated fingerprints for I<text | all> values of B<--output> option.

This is only used for I<Specify> value of B<--DataFieldsMode> option.

Examples:

    Extreg
    MolID,CompoundName

=item B<-d, --DataFieldsMode> I<All | Common | Specify | CompoundID>

Specify how data fields in I<SDFile(s)> are transferred to output CSV/TSV text file(s) along
with generated fingerprints for I<text | all> values of B<--output> option: transfer all SD
data field; transfer SD data files common to all compounds; extract specified data fields;
generate a compound ID using molname line, a compound prefix, or a combination of both.
Possible values: I<All | Common | specify | CompoundID>. Default value: I<CompoundID>.

=item B<--DetectAromaticity> I<Yes | No>

Detect aromaticity before generating fingerprints. Possible values: I<Yes or No>.
Default value: I<Yes>.

I<No> B<--DetectAromaticity> forces usage of atom and bond aromaticity values
from I<SDFile(s)> and skips the step which detects and assigns aromaticity.

I<No> B<--DetectAromaticity> value is only allowed uring I<AtomicInvariantsAtomTypes>
value of B<-a, --AtomIdentifierType> options; for all possible values B<-a, --AtomIdentifierType>
values, it must be I<Yes>.

=item B<-f, --Filter> I<Yes | No>

Specify whether to check and filter compound data in SDFile(s). Possible values: I<Yes or No>.
Default value: I<Yes>.

By default, compound data is checked before calculating fingerprints and compounds containing
atom data corresponding to non-element symbols or no atom data are ignored.

=item B<--FingerprintsLabel> I<text>

SD data label or text file column label to use for fingerprints string in output SD or
CSV/TSV text file(s) specified by B<--output>. Default value: I<PathLenghFingerprints>.

=item B<--fold> I<Yes | No>

Fold fingerprints to increase bit density during I<PathLengthBits> value of
B<-m, --mode> option. Possible values: I<Yes or No>. Default value: I<No>.

=item B<--FoldedSize> I<number>

Size of folded fingerprint during I<PathLengthBits> value of B<-m, --mode> option. Default
value: I<256>. Valid values correspond to any positive integer which is less than
B<-s, --size> and meets the criteria for its value.

Examples:

   128
   512

=item B<-h, --help>

Print this help message

=item B<-i, --IgnoreHydrogens> I<Yes | No>

Ignore hydrogens during fingerprints generation. Possible values: I<Yes or No>.
Default value: I<Yes>.

For I<yes> value of B<-i, --IgnoreHydrogens>, any explicit hydrogens are also used for
generation of atoms path lengths and fingerprints; implicit hydrogens are still ignored.

=item B<-k, --KeepLargestComponent> I<Yes | No>

Generate fingerprints for only the largest component in molecule. Possible values:
I<Yes or No>. Default value: I<Yes>.

For molecules containing multiple connected components, fingerprints can be generated
in two different ways: use all connected components or just the largest connected
component. By default, all atoms except for the largest connected component are
deleted before generation of fingerprints.

=item B<-m, --mode> I<PathLengthBits | PathLengthCount>

Specify type of path length fingerprints to generate for molecules in I<SDFile(s)>. Possible
values: I<PathLengthBits, PathLengthCount>. Default value: I<PathLengthBits>.

For I<PathLengthBits> value of B<-m, --mode> option, a fingerprint bit-vector string containing
zeros and ones is generated and for I<PathLengthCount> value, a fingerprint vector string
corresponding to number of atom paths is generated.

=item B<--MinPathLength> I<number>

Minimum atom path length to include in fingerprints. Default value: I<1>. Valid values:
positive integers and less than B<--MaxPathLength>. Path length of 1 correspond to
a path containing only one atom.

=item B<--MaxPathLength> I<number>

Maximum atom path length to include in fingerprints. Default value: I<8>. Valid values:
positive integers and greater than B<--MinPathLength>.

=item B<-n, --NumOfBitsToSetPerPath> I<number>

Number of bits to set per path during generation of fingerprints bit-vector string for I<PathLengthBits>
value of B<-m, --mode> option. Default value: I<1>. Valid values: positive integers.

=item B<--OutDelim> I<comma | tab | semicolon>

Delimiter for output CSV/TSV text file(s). Possible values: I<comma, tab, or semicolon>
Default value: I<comma>.

=item B<--output> I<SD | FP | text | all>

Type of output files to generate. Possible values: I<SD, FP, text, or all>. Default value: I<text>.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-p, --PathMode> I<AtomPathsWithoutRings | AtomPathsWithRings | AllAtomPathsWithoutRings | AllAtomPathsWithRings>

Specify type of atom paths to use for generating  pathlength fingerprints for molecules in
I<SDFile(s)>. Possible values:I<AtomPathsWithoutRings, AtomPathsWithRings,
AllAtomPathsWithoutRings, AllAtomPathsWithRings>. Default value: I<AllAtomPathsWithRings>.

For molecules with no rings, first two and last two options are equivalent and generate
same set of atom paths starting from each atom with length between B<--MinPathLength>
and B<--MaxPathLength>. However, all these four options can result in the same set of
final atom paths for molecules containing fused, bridged or spiro rings.

For molecules containing rings, atom paths starting from each atom can be traversed in
four different ways:

I<AtomPathsWithoutRings> - Atom paths containing no rings and without sharing of bonds
in traversed paths.

I<AtomPathsWithRings> - Atom paths containing rings and without any sharing of bonds in
traversed paths.

I<AllAtomPathsWithoutRings> - All possible atom paths containing no rings and without any
sharing of bonds in traversed paths.

I<AllAtomPathsWithRings> - All possible atom paths containing rings and with sharing of
bonds in traversed paths.

Atom path traversal is terminated at the ring atom.

Based on values specified for for B<-p, --PathMode>, B<--MinPathLength> and
B<--MaxPathLength>, all appropriate atom paths are generated for each atom in the molecule
and collected in a list.

For each atom path in the filtered atom paths list, an atom path string is created using value of
B<-a, --AtomIdentifierType> and specified values to use for a particular atom identifier type.
Value of B<-u, --UseBondSymbols> controls whether bond order symbols are used during generation
of atom path string. Atom symbol corresponds to element symbol and characters used to represent
 bond order are: I<1 - None; 2 - '='; 3 - '#'; 1.5 or aromatic - ':'; others: bond order value>. By default,
bond symbols are included in atom path strings. Exclusion of bond symbols in atom path strings
results in fingerprints which correspond purely to atom paths without considering bonds.

B<UseUniquePaths> controls the removal of structurally duplicate atom path strings are removed
from the list.

For I<PathLengthBits> value of B<-m, --mode> option, each atom path is hashed to a 32 bit unsigned
integer key using B<TextUtil::HashCode> function. Using the hash key as a seed for a random number
generator, a random integer value between 0 and B<--Size> is used to set corresponding bits
in the fingerprint bit-vector string. Value of B<--NumOfBitsToSetPerPaths> option controls the number
of time a random number is generated to set corresponding bits.

For I< PathLengthCount> value of B<-m, --mode> option, the number of times an atom path appears
is tracked and a fingerprints count-string corresponding to count of atom paths is generated.

For molecule containing rings, combination of B<-p, --PathMode> and B<--UseBondSymbols> allows
generation of up to 8 different types of atom path length strings:

    AllowSharedBonds AllowRings UseBondSymbols

    0                0          1   - AtomPathsNoCyclesWithBondSymbols
    0                1          1   - AtomPathsWithCyclesWithBondSymbols

    1                0          1   - AllAtomPathsNoCyclesWithBondSymbols
    1                1          1   - AllAtomPathsWithCyclesWithBondSymbols
                                      [ DEFAULT ]

    0                0          0   - AtomPathsNoCyclesNoBondSymbols
    0                1          0   - AtomPathsWithCyclesNoBondSymbols

    1                0          0   - AllAtomPathsNoCyclesNoBondSymbols
    1                1          0   - AllAtomPathsWithCyclesNoWithBondSymbols

Default atom path length fingerprints generation for molecules containing rings with
I<AllAtomPathsWithRings> value for B<-p, --PathMode>, I<Yes> value for B<--UseBondSymbols>,
I<2> value for B<--MinPathLength> and I<8> value for B<--MaxPathLength> is the most time
consuming. Combinations of other options can substantially speed up fingerprint generation
for molecules containing complex ring systems.

Additionally, value for option B<-a, --AtomIdentifierType> in conjunction with corresponding specified
values for atom types changes the nature of atom path length strings and the fingerprints.

=item B<-q, --quote> I<Yes | No>

Put quote around column values in output CSV/TSV text file(s). Possible values:
I<Yes or No>. Default value: I<Yes>.

=item B<-r, --root> I<RootName>

New file name is generated using the root: <Root>.<Ext>. Default for new file
names: <SDFileName><PathLengthFP>.<Ext>. The file type determines <Ext> value.
The sdf, fpf, csv, and tsv <Ext> values are used for SD, FP, comma/semicolon, and tab
delimited text files, respectively.This option is ignored for multiple input files.

=item B<-s, --size> I<number>

Size of fingerprints. Default value: I<1024>. Valid values correspond to any positive
integer which satisfies the following criteria: power of 2, >= 32 and <= 2 ** 32.

Examples:

   256
   512
   2048

=item B<-u, --UseBondSymbols> I<Yes | No>

Specify whether to use bond symbols for atom paths during generation of atom path strings.
Possible values: I<Yes or No>. Default value: I<Yes>.

I<No> value option for B<-u, --UseBondSymbols> allows the generation of fingerprints corresponding
purely to atoms disregarding all bonds.

=item B<--UsePerlCoreRandom> I<Yes | No>

Specify whether to use Perl CORE::rand or MayaChemTools MathUtil::random function
during random number generation for setting bits in fingerprints bit-vector strings. Possible
values: I<Yes or No>. Default value: I<Yes>.

I<No> value option for B<--UsePerlCoreRandom> allows the generation of fingerprints
bit-vector strings which are same across different platforms.

The random number generator implemented in MayaChemTools is a variant of
linear congruential generator (LCG) as described by Miller et al. [ Ref 120 ].
It is also referred to as Lehmer random number generator or Park-Miller
random number generator.

Unlike Perl's core random number generator function rand, the random number
generator implemented in MayaChemTools, MathUtil::random,  generates consistent
random values across different platforms for a specific random seed and leads
to generation of portable fingerprints bit-vector strings.

=item B<--UseUniquePaths> I<Yes | No>

Specify whether to use structurally unique atom paths during generation of atom path strings.
Possible values: I<Yes or No>. Default value: I<Yes>.

I<No> value option for B<--UseUniquePaths> allows usage of all atom paths generated by
B<-p, --PathMode> option value for generation of atom path strings leading to duplicate
path count during I<PathLengthCount> value of B<-m, --mode> option. It doesn't change fingerprint
string generated during I<PathLengthBits> value of B<-m, --mode>.

For example, during I<AllAtomPathsWithRings> value of B<-p, --PathMode> option, benzene has
12 linear paths of length 2 and 12 cyclic paths length of 7, but only 6 linear paths of length 2 and
1 cyclic path of length 7 are structurally unique.

=item B<-v, --VectorStringFormat> I<IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>

Format of fingerprints vector string data in output SD, FP or CSV/TSV text file(s) specified by
B<--output> used during I<PathLengthCount> value of B<-m, --mode> option. Possible
values: I<IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString |
ValuesAndIDsPairsString>. Defaultvalue: I<IDsAndValuesString>.

Examples:

    FingerprintsVector;PathLengthCount:AtomicInvariantsAtomTypes:MinLength
    1:MaxLength8;432;NumericalValues;IDsAndValuesPairsString;C.X1.BO1.H3 2
     C.X2.BO2.H2 4 C.X2.BO3.H1 14 C.X3.BO3.H1 3 C.X3.BO4 10 F.X1.BO1 1 N.X
    2.BO2.H1 1 N.X3.BO3 1 O.X1.BO1.H1 3 O.X1.BO2 2 C.X1.BO1.H3C.X3.BO3.H1 
    2 C.X2.BO2.H2C.X2.BO2.H2 1 C.X2.BO2.H2C.X3.BO3.H1 4 C.X2.BO2.H2C.X3.BO
    4 1 C.X2.BO2.H2N.X3.BO3 1 C.X2.BO3.H1:C.X2.BO3.H1 10 C.X2.BO3.H1:C....

    FingerprintsVector;PathLengthCount:EStateAtomTypes:MinLength1:MaxLengt
    h8;454;NumericalValues;IDsAndValuesPairsString;aaCH 14 aasC 8 aasN 1 d
    O 2 dssC 2 sCH3 2 sF 1 sOH 3 ssCH2 4 ssNH 1 sssCH 3 aaCH:aaCH 10 aaCH:
    aasC 8 aasC:aasC 3 aasC:aasN 2 aasCaasC 2 aasCdssC 1 aasCsF 1 aasCssNH
    1 aasCsssCH 1 aasNssCH2 1 dO=dssC 2 dssCsOH 1 dssCssCH2 1 dssCssNH 1
    sCH3sssCH 2 sOHsssCH 2 ssCH2ssCH2 1 ssCH2sssCH 4 aaCH:aaCH:aaCH 6 a...

=item B<-w, --WorkingDir> I<DirName>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To generate path length fingerprints corresponding to all unique paths from length 1
through 8 in hexadecimal bit-vector string format of size 1024 and create a
SamplePLFPHex.csv file containing sequential compound IDs along with fingerprints
bit-vector strings data, type:

    % PathLengthFingerprints.pl -o -r SamplePLFPHex Sample.sdf

To generate path length fingerprints corresponding to all unique paths from length 1
through 8 in hexadecimal bit-vector string format of size 1024 and create SamplePLFPHex.sdf,
SamplePLFPHex.fpf, and SamplePLFPHex.csv files containing sequential compound IDs
in CSV file along with fingerprints bit-vector strings data, type:

    % PathLengthFingerprints.pl --output all -o -r SamplePLFPHex Sample.sdf

To generate path length fingerprints corresponding to all unique paths from length 1
through 8 in binary bit-vector string format of size 1024 and create a
SamplePLFPBin.csv file containing sequential compound IDs along with fingerprints
bit-vector strings data, type:

    % PathLengthFingerprints.pl --BitStringFormat BinaryString --size 2048
      -o -r SamplePLFPBin Sample.sdf

To generate path length fingerprints corresponding to count of all unique paths from
length 1 through 8 in IDsAndValuesString format and create a SamplePLFPCount.csv file
containing sequential compound IDs along with fingerprints vector strings data, type:

    % PathLengthFingerprints.pl -m PathLengthCount -o -r SamplePLFPCount
      Sample.sdf

To generate path length fingerprints corresponding to count of all unique paths from
length 1 through 8 in IDsAndValuesString format using  E-state atom types and
create a SamplePLFPCount.csv file containing sequential compound IDs along with fingerprints
vector strings data, type:

    % PathLengthFingerprints.pl -m  PathLengthCount --AtomIdentifierType
      EStateAtomTypes -o -r SamplePLFPCount Sample.sdf

To generate path length fingerprints corresponding to count of all unique paths from
length 1 through 8 in IDsAndValuesString format using  SLogP atom types and
create a SamplePLFPCount.csv file containing sequential compound IDs along with fingerprints
vector strings data, type:

    % PathLengthFingerprints.pl -m  PathLengthCount --AtomIdentifierType
      SLogPAtomTypes -o -r SamplePLFPCount Sample.sdf

To generate path length fingerprints corresponding to count of all unique paths from
length 1 through 8 in IDsAndValuesString format and create a SamplePLFPCount.csv file
containing sequential compound IDs along with fingerprints vector strings data, type:

    % PathLengthFingerprints.pl -m  PathLengthCount --VectorStringFormat
      ValuesAndIDsPairsString -o -r SamplePLFPCount Sample.sdf

To generate path length fingerprints corresponding to count of all unique paths from
length 1 through 8 in IDsAndValuesString format using AS,X,BO as atomic invariants and
create a SamplePLFPCount.csv file containing sequential compound IDs along with fingerprints
vector strings data, type:

    % PathLengthFingerprints.pl -m  PathLengthCount --AtomIdentifierType
      AtomicInvariantsAtomTypes --AtomicInvariantsToUse "AS,X,BO" -o
      -r SamplePLFPCount Sample.sdf

To generate path length fingerprints corresponding to count of all  paths from
length 1 through 8 in IDsAndValuesString format and create a SamplePLFPCount.csv file
containing compound IDs from MolName line along with fingerprints vector strings data, type:

    % PathLengthFingerprints.pl -m  PathLengthCount --UseUniquePaths No
      -o --CompoundIDMode MolName -r SamplePLFPCount --UseUniquePaths No
      Sample.sdf

To generate path length fingerprints corresponding to all unique paths from length 1
through 8 in hexadecimal bit-vector string format of size 512 after folding and create
SamplePLFPHex.sdf, SamplePLFPHex.fpf, and SamplePLFPHex.sdf  files containing sequential
compound IDs along with fingerprints bit-vector strings data, type:

    % PathLengthFingerprints.pl --output all --Fold Yes --FoldedSize 512
       -o -r SamplePLFPHex Sample.sdf

To generate path length fingerprints corresponding to all unique paths from length 1
through 8 containing no rings and  without sharing of bonds in hexadecimal bit-vector
string format of size 1024 and create a SamplePLFPHex.csv file containing sequential
compound IDs along with fingerprints bit-vector strings data and all data fields, type:

    % PathLengthFingerprints.pl -p AtomPathsWithoutRings --DataFieldsMode All
      -o -r SamplePLFPHex Sample.sdf

To generate path length fingerprints corresponding to all unique paths from length 1
through 8 containing rings and  without sharing of bonds in hexadecimal bit-vector
string format of size 1024 and create a SamplePLFPHex.tsv file containing compound IDs
derived from combination of molecule name line and an explicit compound prefix
along with fingerprints bit-vector strings data and all data fields, type:

    % PathLengthFingerprints.pl -p AtomPathsWithRings --DataFieldsMode
      CompoundID --CompoundIDMode MolnameOrLabelPrefix  --CompoundID Cmpd
      --CompoundIDLabel MolID --FingerprintsLabel PathLengthFP --OutDelim Tab
      -r SamplePLFPHex -o Sample.sdf

To generate path length fingerprints corresponding to count of all unique paths from
length 1 through 8 in IDsAndValuesString format and create a SamplePLFPCount.csv file
containing sequential compound IDs along with fingerprints vector strings data using
aromaticity specified in SD file, type:

    % PathLengthFingerprints.pl -m PathLengthCount --DetectAromaticity No
      -o -r SamplePLFPCount Sample.sdf

To generate path length fingerprints corresponding to all unique paths from length 2
through 6 in hexadecimal bit-vector string format of size 1024 and create a
SamplePLFPHex.csv file containing sequential compound IDs along with fingerprints
bit-vector strings data, type:

    % PathLengthFingerprints.pl --MinPathLength 2 --MaxPathLength 6
      -o -r SamplePLFPHex Sample.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoFingerprintsFiles.pl, SimilarityMatricesFingerprints.pl, AtomNeighborhoodsFingerprints.pl,
ExtendedConnectivityFingerprints.pl, MACCSKeysFingerprints.pl,
TopologicalAtomPairsFingerprints.pl, TopologicalAtomTorsionsFingerprints.pl,
TopologicalPharmacophoreAtomPairsFingerprints.pl, TopologicalPharmacophoreAtomTripletsFingerprints.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
