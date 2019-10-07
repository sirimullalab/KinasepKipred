#!/usr/bin/perl -w
#
# File: AtomTypesFingerprints.pl
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
use Fingerprints::AtomTypesFingerprints;

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
    GenerateAtomTypesFingerprints($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Generate fingerprints for a SD file...
#
sub GenerateAtomTypesFingerprints {
  my($FileIndex) = @_;
  my($CmpdCount, $IgnoredCmpdCount, $SDFile, $MoleculeFileIO, $Molecule, $AtomTypesFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO);

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

    $AtomTypesFingerprints = GenerateMoleculeFingerprints($Molecule);
    if (!$AtomTypesFingerprints) {
      $IgnoredCmpdCount++;
      ProcessIgnoredCompound('FingerprintsGenerationFailed', $CmpdCount, $Molecule);
      next COMPOUND;
    }

    WriteDataToOutputFiles($FileIndex, $CmpdCount, $Molecule, $AtomTypesFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO);
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
  if ($OptionsInfo{Mode} =~ /^AtomTypesBits$/i) {
    %FingerprintsFileIOParams = ('Mode' => 'Write', 'Overwrite' => $OptionsInfo{OverwriteFiles}, 'FingerprintsStringMode' => 'FingerprintsBitVectorString', 'BitStringFormat' => $OptionsInfo{BitStringFormat}, 'BitsOrder' => $OptionsInfo{BitsOrder});
  }
  elsif ($OptionsInfo{Mode} =~ /^AtomTypesCount$/i) {
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
  my($FileIndex, $CmpdCount, $Molecule, $AtomTypesFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO) = @_;
  my($DataFieldLabelAndValuesRef);

  $DataFieldLabelAndValuesRef = undef;
  if ($NewFPTextFileIO || $NewFPFileIO) {
    $DataFieldLabelAndValuesRef = $Molecule->GetDataFieldLabelAndValues();
  }

  if ($NewFPSDFileIO) {
    my($CmpdString);

    $CmpdString = $Molecule->GetInputMoleculeString();
    $NewFPSDFileIO->WriteFingerprints($AtomTypesFingerprints, $CmpdString);
  }

  if ($NewFPTextFileIO) {
    my($ColValuesRef);

    $ColValuesRef = SetupFPTextFileCoulmnValues($FileIndex, $CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPTextFileIO->WriteFingerprints($AtomTypesFingerprints, $ColValuesRef);
  }

  if ($NewFPFileIO) {
    my($CompoundID);

    $CompoundID = SetupCmpdIDForOutputFiles($CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPFileIO->WriteFingerprints($AtomTypesFingerprints, $CompoundID);
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
  my($AtomTypesFingerprints);

  if ($OptionsInfo{KeepLargestComponent}) {
    $Molecule->KeepLargestComponent();
  }
  if (!$Molecule->DetectRings()) {
    return undef;
  }
  $Molecule->SetAromaticityModel($OptionsInfo{AromaticityModel});
  $Molecule->DetectAromaticity();

  $AtomTypesFingerprints = undef;
  if ($OptionsInfo{Mode} =~ /^AtomTypesCount$/i) {
    $AtomTypesFingerprints = new Fingerprints::AtomTypesFingerprints('Molecule' => $Molecule, 'Type' => 'AtomTypesCount', 'AtomIdentifierType' => $OptionsInfo{AtomIdentifierType}, 'AtomTypesSetToUse' => $OptionsInfo{AtomTypesSetToUse}, 'IgnoreHydrogens' => $OptionsInfo{IgnoreHydrogens});

  }
  elsif ($OptionsInfo{Mode} =~ /^AtomTypesBits$/i) {
    $AtomTypesFingerprints = new Fingerprints::AtomTypesFingerprints('Molecule' => $Molecule, 'Type' => 'AtomTypesBits', 'AtomIdentifierType' => $OptionsInfo{AtomIdentifierType}, 'AtomTypesSetToUse' => 'FixedSize', 'IgnoreHydrogens' => $OptionsInfo{IgnoreHydrogens});
  }
  else {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: AtomTypesCount or AtomTypesBits\n";
  }

  SetAtomIdentifierTypeValuesToUse($AtomTypesFingerprints);

  # Generate atom types fingerprints...
  $AtomTypesFingerprints->GenerateFingerprints();

  # Make sure atom types fingerprints generation is successful...
  if (!$AtomTypesFingerprints->IsFingerprintsGenerationSuccessful()) {
    return undef;
  }

  return $AtomTypesFingerprints;
}

# Set atom identifier type to use for generating fingerprints...
#
sub SetAtomIdentifierTypeValuesToUse {
  my($AtomTypesFingerprints) = @_;

  if ($OptionsInfo{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
    $AtomTypesFingerprints->SetAtomicInvariantsToUse(\@{$OptionsInfo{AtomicInvariantsToUse}});
  }
  elsif ($OptionsInfo{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
    $AtomTypesFingerprints->SetFunctionalClassesToUse(\@{$OptionsInfo{FunctionalClassesToUse}});
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
      $OutFileRoot = $FileName . 'AtomTypesFP';
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
      # Check SD and text outout files...
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

  ProcessAtomIdentifierTypeOptions();

  my($AtomTypesSetToUse);
  $AtomTypesSetToUse = '';
  if ($Options{mode} =~ /^AtomTypesBits$/i) {
    if ($Options{atomtypessettouse} && $Options{atomtypessettouse} !~ /^FixedSize$/) {
      die "Error: The value specified, $Options{atomtypessettouse}, for option \"-e, --AtomTypesSetToUse\" is not valid. Allowed values for AtomTypesBits of \"-m, --mode\" option: FixedSize\n";
    }
    $AtomTypesSetToUse = 'FixedSize';
  }
  else {
    if ($Options{atomidentifiertype} =~ /^(AtomicInvariantsAtomTypes|FunctionalClassAtomTypes)$/i && $Options{atomtypessettouse} =~ /^FixedSize$/) {
      die "Error: The value specified, $Options{atomtypessettouse}, for option \"-e, --AtomTypesSetToUse\" is not valid during \"AtomicInvariantsAtomTypes or FunctionalClassAtomTypes\" value of \"-a, --AtomIdentifierType\". Allowed values: ArbitrarySize\n";
    }
    if ($Options{atomidentifiertype} =~ /^TPSAAtomTypes$/i && $Options{atomtypessettouse} =~ /^ArbitrarySize$/) {
      die "Error: The value specified, $Options{atomtypessettouse}, for option \"-e, --AtomTypesSetToUse\" is not valid during \"TPSAAtomTypes\" value of \"-a, --AtomIdentifierType\". Allowed values: FixedSize\n";
    }
    $AtomTypesSetToUse = $Options{atomtypessettouse} ? $Options{atomtypessettouse} : 'ArbitrarySize';
  }
  $OptionsInfo{AtomTypesSetToUse} = $AtomTypesSetToUse;

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

  $OptionsInfo{IgnoreHydrogens} = ($Options{ignorehydrogens} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{FingerprintsLabel} = $Options{fingerprintslabel} ? $Options{fingerprintslabel} : 'AtomTypesFingerprints';

  $OptionsInfo{Filter} = ($Options{filter} =~ /^Yes$/i) ? 1 : 0;

  if ($Options{fingerprintslabelmode} =~ /^FingerprintsLabelWithIDs$/) {
    if ($Options{mode} =~ /^(AtomTypesCount)$/i && $Options{atomtypessettouse} =~ /^FixedSize$/i) {
      # Append atom types to the fingerprints label...
      my($FixedSizeAtomTypesSetRef);
      $FixedSizeAtomTypesSetRef = GetFixedSizeAtomTypesSet();

      $OptionsInfo{FingerprintsLabel} .= "; AtomTypes: " . TextUtil::JoinWords($FixedSizeAtomTypesSetRef, " ", 0);
    }
  }
  $OptionsInfo{FingerprintsLabelMode} = $Options{fingerprintslabelmode};

  $OptionsInfo{KeepLargestComponent} = ($Options{keeplargestcomponent} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{Output} = $Options{output};
  $OptionsInfo{SDOutput} = ($Options{output} =~ /^(SD|All)$/i) ? 1 : 0;
  $OptionsInfo{FPOutput} = ($Options{output} =~ /^(FP|All)$/i) ? 1 : 0;
  $OptionsInfo{TextOutput} = ($Options{output} =~ /^(Text|All)$/i) ? 1 : 0;

  $OptionsInfo{OutDelim} = $Options{outdelim};
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? 1 : 0;
  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : 0;

  # Setup default vector string format...
  my($VectorStringFormat);
  $VectorStringFormat = '';
  if ($Options{vectorstringformat}) {
    $VectorStringFormat = $Options{vectorstringformat};
  }
  else {
    $VectorStringFormat = ($Options{atomtypessettouse} =~ /^FixedSize$/) ? "ValuesString" : "IDsAndValuesString";
  }
  $OptionsInfo{VectorStringFormat} = $VectorStringFormat;
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

# Get fixed size atom types set...
#
sub GetFixedSizeAtomTypesSet {
  my($AtomTypesRef);

  $AtomTypesRef = undef;

  IDENTIFIERTYPE: {
    if ($OptionsInfo{AtomIdentifierType} =~ /^DREIDINGAtomTypes$/i) {
      $AtomTypesRef = $OptionsInfo{IgnoreHydrogens} ? DREIDINGAtomTypes::GetAllPossibleDREIDINGNonHydrogenAtomTypes() : DREIDINGAtomTypes::GetAllPossibleDREIDINGAtomTypes();
      last IDENTIFIERTYPE;
    }

    if ($OptionsInfo{AtomIdentifierType} =~ /^EStateAtomTypes$/i) {
      $AtomTypesRef = $OptionsInfo{IgnoreHydrogens} ? EStateAtomTypes::GetAllPossibleEStateNonHydrogenAtomTypes() : EStateAtomTypes::GetAllPossibleEStateAtomTypes();
      last IDENTIFIERTYPE;
    }

    if ($OptionsInfo{AtomIdentifierType} =~ /^MMFF94AtomTypes$/i) {
      $AtomTypesRef = $OptionsInfo{IgnoreHydrogens} ? MMFF94AtomTypes::GetAllPossibleMMFF94NonHydrogenAtomTypes() : MMFF94AtomTypes::GetAllPossibleMMFF94AtomTypes();
      last IDENTIFIERTYPE;
    }

    if ($OptionsInfo{AtomIdentifierType} =~ /^SLogPAtomTypes$/i) {
      $AtomTypesRef = $OptionsInfo{IgnoreHydrogens} ? SLogPAtomTypes::GetAllPossibleSLogPNonHydrogenAtomTypes() : SLogPAtomTypes::GetAllPossibleSLogPAtomTypes();
      last IDENTIFIERTYPE;
    }

    if ($OptionsInfo{AtomIdentifierType} =~ /^SYBYLAtomTypes$/i) {
      $AtomTypesRef = $OptionsInfo{IgnoreHydrogens} ? SYBYLAtomTypes::GetAllPossibleSYBYLNonHydrogenAtomTypes() : SYBYLAtomTypes::GetAllPossibleSYBYLAtomTypes();
      last IDENTIFIERTYPE;
    }

    if ($OptionsInfo{AtomIdentifierType} =~ /^TPSAAtomTypes$/i) {
      $AtomTypesRef = TPSAAtomTypes::GetAllPossibleTPSAAtomTypes();
      last IDENTIFIERTYPE;
    }

    if ($OptionsInfo{AtomIdentifierType} =~ /^UFFAtomTypes$/i) {
      $AtomTypesRef = $OptionsInfo{IgnoreHydrogens} ? UFFAtomTypes::GetAllPossibleUFFNonHydrogenAtomTypes() : UFFAtomTypes::GetAllPossibleUFFAtomTypes();
      last IDENTIFIERTYPE;
    }
    die "Error: GetFixedSizeAtomTypesSet: Atom types set for atom indentifier type, $OptionsInfo{AtomIdentifierType}, is not available...";
  }

  return $AtomTypesRef;
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{aromaticitymodel} = 'MayaChemToolsAromaticityModel';

  $Options{atomidentifiertype} = 'AtomicInvariantsAtomTypes';
  $Options{atomicinvariantstouse} = 'AS,X,BO,H,FC';
  $Options{functionalclassestouse} = 'HBD,HBA,PI,NI,Ar,Hal';

  $Options{atomtypessettouse} = 'ArbitrarySize';

  $Options{bitsorder} = 'Ascending';
  $Options{bitstringformat} = 'BinaryString';

  $Options{compoundidmode} = 'LabelPrefix';
  $Options{compoundidlabel} = 'CompoundID';
  $Options{datafieldsmode} = 'CompoundID';

  $Options{filter} = 'Yes';

  $Options{fingerprintslabelmode} = 'FingerprintsLabelOnly';
  $Options{keeplargestcomponent} = 'Yes';

  $Options{mode} = 'AtomTypesCount';

  $Options{ignorehydrogens} = 'Yes';

  $Options{quote} = 'yes';

  $Options{output} = 'text';
  $Options{outdelim} = 'comma';
  $Options{quote} = 'yes';

  $Options{vectorstringformat} = '';

  if (!GetOptions(\%Options, "aromaticitymodel=s", "atomidentifiertype|a=s", "atomicinvariantstouse=s", "functionalclassestouse=s", "atomtypessettouse|e=s", "bitsorder=s", "bitstringformat|b=s", "compoundid=s", "compoundidlabel=s", "compoundidmode=s", "datafields=s", "datafieldsmode|d=s", "filter|f=s", "fingerprintslabelmode=s", "fingerprintslabel=s",  "help|h", "ignorehydrogens|i=s", "keeplargestcomponent|k=s", "mode|m=s", "outdelim=s", "output=s", "overwrite|o", "quote|q=s", "root|r=s", "vectorstringformat|v=s", "workingdir|w=s")) {
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
  if ($Options{atomtypessettouse} && $Options{atomtypessettouse} !~ /^(ArbitrarySize|FixedSize)$/) {
    die "Error: The value specified, $Options{atomtypessettouse}, for option \"--AtomTypesSetToUse\" is not valid. Allowed values: ArbitrarySize or FixedSize\n";
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
  if ($Options{filter} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{filter}, for option \"-f, --Filter\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{fingerprintslabelmode} !~ /^(FingerprintsLabelOnly|FingerprintsLabelWithIDs)$/i) {
    die "Error: The value specified, $Options{fingerprintslabelmode}, for option \"--FingerprintsLabelMode\" is not valid. Allowed values: FingerprintsLabelOnly or FingerprintsLabelWithIDs\n";
  }
  if ($Options{ignorehydrogens} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{ignorehydrogens}, for option \"-i, --IgnoreHydrogens\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{keeplargestcomponent} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{keeplargestcomponent}, for option \"-k, --KeepLargestComponent\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{mode} !~ /^(AtomTypesCount|AtomTypesBits)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: AtomTypesCount, or AtomTypesBits\n";
  }
  if ($Options{output} !~ /^(SD|FP|text|all)$/i) {
    die "Error: The value specified, $Options{output}, for option \"--output\" is not valid. Allowed values: SD, FP, text, or all\n";
  }
  if ($Options{outdelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{outdelim}, for option \"--outdelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if ($Options{quote} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{outdelim} =~ /semicolon/i && $Options{quote} =~ /^No$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not allowed with, semicolon value of \"--outdelim\" option: Fingerprints string use semicolon as delimiter for various data fields and must be quoted.\n";
  }
  if ($Options{vectorstringformat} && $Options{vectorstringformat} !~ /^(ValuesString|IDsAndValuesString|IDsAndValuesPairsString|ValuesAndIDsString|ValuesAndIDsPairsString)$/i) {
    die "Error: The value specified, $Options{vectorstringformat}, for option \"-v, --VectorStringFormat\" is not valid. Allowed values: ValuesString, IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString or ValuesAndIDsPairsString\n";
  }
}

__END__

=head1 NAME

AtomTypesFingerprints.pl - Generate atom types fingerprints for SD files

=head1 SYNOPSIS

AtomTypesFingerprints.pl SDFile(s)...

AtomTypesFingerprints.pl [B<--AromaticityModel> I<AromaticityModelType>]
[B<-a, --AtomIdentifierType> I<AtomicInvariantsAtomTypes |
DREIDINGAtomTypes | EStateAtomTypes | MMFF94AtomTypes | SLogPAtomTypes | SYBYLAtomTypes | TPSAAtomTypes | UFFAtomTypes>]
[B<--AtomicInvariantsToUse> I<"AtomicInvariant, AtomicInvariant...">]
[B<--FunctionalClassesToUse> I<"FunctionalClass1,FunctionalClass2...">]
[B<--AtomTypesSetToUse> I<ArbitrarySize | FixedSize>]
[B<--BitsOrder> I<Ascending | Descending>] [B<-b, --BitStringFormat> I<BinaryString | HexadecimalString>]
[B<--CompoundID> I<DataFieldName or LabelPrefixString>] [B<--CompoundIDLabel> I<text>]
[B<--CompoundIDMode> I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>]
[B<--DataFields> I<"FieldLabel1,FieldLabel2,...">] [B<-d, --DataFieldsMode> I<All | Common | Specify | CompoundID>]
[B<-f, --Filter> I<Yes | No>] [B<--FingerprintsLabelMode> I<FingerprintsLabelOnly | FingerprintsLabelWithIDs>] [B<--FingerprintsLabel> I<text>]
[B<-h, --help>] [B<-k, --KeepLargestComponent> I<Yes | No>]
[B<-m, --mode> I<AtomTypesCount | AtomTypesBits>] [B<-i, --IgnoreHydrogens> I<Yes | No>]
[B<--OutDelim> I<comma | tab | semicolon>] [B<--output> I<SD |FP | text | all>] [B<-o, --overwrite>]
[B<-q, --quote> I<Yes | No>] [B<-r, --root> I<RootName>] [B<-s, --size> I<number>] [B<--ValuesPrecision> I<number>]
[B<-v, --VectorStringFormat> I<IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>]
[B<-w, --WorkingDir> I<DirName>]

=head1 DESCRIPTION

Generate atom types fingerprints for I<SDFile(s)> and create appropriate SD, FP or
CSV/TSV text file(s) containing fingerprints bit-vector or vector strings corresponding to
molecular fingerprints.

Multiple SDFile names are separated by spaces. The valid file extensions are I<.sdf>
and I<.sd>. All other file names are ignored. All the SD files in a current directory
can be specified either by I<*.sdf> or the current directory name.

The current release of MayaChemTools supports generation of atom types fingerpritns
corresponding to following B<-a, --AtomIdentifierTypes>:

    AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
    FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
    SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes

Based on the values specified for B<-a, --AtomIdentifierType> along with other specified
parameters such as B<--AtomicInvariantsToUse> and B<--FunctionalClassesToUse>, initial
atom types are assigned to all non-hydrogen atoms or all atoms in a molecule

Using the assigned atom types and specified B<-m, --Mode>, one of the following types of
fingerprints are generated:

    AtomTypesCount - A vector containing count of atom types
    AtomTypesBits - A bit vector indicating presence/absence of atom types

For I<AtomTypesCount> fingerprints, two types of atom types set size are allowed as
value of B<--AtomTypesSetToUse> option:

    ArbitrarySize - Corresponds to only atom types detected in molecule
    FixedSize - Corresponds to fixed number of atom types previously defined

For I<AtomTypesBits> fingerprints, only I<FixedSize> atom type set is allowed.

I<ArbitrarySize> corresponds to atom types detected in a molecule where as I<FixedSize> implies
a fix number of all possible atom types previously defined for a specific B<-a, --AtomIdentifierType>.

Fix number of all possible atom types for supported I<AtomIdentifierTypes> in current release
of MayaChemTools are:

    AtomIdentifier       Total    TotalWithoutHydrogens

    DREIDINGAtomTypes    37       34
    EStateAtomTypes      109      87
    MMFF94AtomTypes      212      171
    SLogPAtomTypes       72       67
    SYBYLAtomTypes       45       44
    TPSAAtomTypes        47       47
    UFFAtomTypes         126      124

The current release of MayaChemTools generates the following atom types fingerprints
bit-vector and vector strings:

    FingerprintsVector;AtomTypesCount:AtomicInvariantsAtomTypes:ArbitraryS
    ize;10;NumericalValues;IDsAndValuesString;C.X1.BO1.H3 C.X2.BO2.H2 C.X2
    .BO3.H1 C.X3.BO3.H1 C.X3.BO4 F.X1.BO1 N.X2.BO2.H1 N.X3.BO3 O.X1.BO1.H1
    O.X1.BO2;2 4 14 3 10 1 1 1 3 2

    FingerprintsVector;AtomTypesCount:DREIDINGAtomTypes:ArbitrarySize;8;Nu
    mericalValues;IDsAndValuesString;C_2 C_3 C_R F_ N_3 N_R O_2 O_3;2 9 22
    1 1 1 2 3

    FingerprintsVector;AtomTypesCount:DREIDINGAtomTypes:FixedSize;34;Order
    edNumericalValues;IDsAndValuesString;B_3 B_2 C_3 C_R C_2 C_1 N_3 N_R N
    _2 N_1 O_3 O_R O_2 O_1 F_ Al3 Si3 P_3 S_3 Cl Ga3 Ge3 As3 Se3 Br In3 Sn
    3 Sb3 Te3 I_ Na Ca Fe Zn;0 0 9 22 2 0 1 1 0 0 3 0 2 0 1 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0

    FingerprintsBitVector;AtomTypesBits:DREIDINGAtomTypes:FixedSize;34;Bin
    aryString;Ascending;0011101100101010000000000000000000000000

    FingerprintsVector;AtomTypesCount:EStateAtomTypes:ArbitrarySize;11;Num
    ericalValues;IDsAndValuesString;aaCH aasC aasN dO dssC sCH3 sF sOH ssC
    H2 ssNH sssCH;14 8 1 2 2 2 1 3 4 1 3

    FingerprintsVector;AtomTypesCount:EStateAtomTypes:FixedSize;87;Ordered
    NumericalValues;IDsAndValuesString;sLi ssBe ssssBem sBH2 ssBH sssB sss
    sBm sCH3 dCH2 ssCH2 tCH dsCH aaCH sssCH ddC tsC dssC aasC aaaC ssssC s
    NH3p sNH2 ssNH2p dNH ssNH aaNH tN sssNHp dsN aaN sssN ddsN aasN ss...;
    0 0 0 0 0 0 0 2 0 4 0 0 14 3 0 0 2 8 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 3 2 0 0
    0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0...

    FingerprintsBitVector;AtomTypesBits:EStateAtomTypes:FixedSize;87;Binar
    yString;Ascending;0000000101001100110000001000000010110000100000000000
    000000000000000000000000000000000000

    FingerprintsVector;AtomTypesCount:FunctionalClassAtomTypes:ArbitrarySi
    ze;8;NumericalValues;IDsAndValuesString;Ar Ar.HBA HBA HBA.HBD HBD Hal 
    NI None;22 1 2 3 1 1 1 10

    FingerprintsVector;AtomTypesCount:MMFF94AtomTypes:ArbitrarySize;13;Num
    ericalValues;IDsAndValuesString;C5A C5B C=ON CB COO CR F N5 NC=O O=CN
    O=CO OC=O OR;2 2 1 18 1 9 1 1 1 1 1 1 2

    FingerprintsVector;AtomTypesCount:MMFF94AtomTypes:FixedSize;171;Ordere
    dNumericalValues;IDsAndValuesString;CR C=C CSP2 C=O C=N CGD C=OR C=ON
    CONN COO COON COOO C=OS C=S C=SN CSO2 CS=O CSS C=P CSP =C= OR OC=O OC=
    C OC=N OC=S ONO2 ON=O OSO3 OSO2 OSO OS=O -OS OPO3 OPO2 OPO -OP -O-...;
    9 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 2 1 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 1 0 1 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 0 0 0 0 0 0 0 0 ...

    FingerprintsBitVector;AtomTypesBits:MMFF94AtomTypes:FixedSize;171;Bina
    ryString;Ascending;100000010100000000000110000000000000000101000000100
    0100000000000000000000000000000000000000000100000000000000000000000000
    0000000011000000000000000001000000000000000000000000000

    FingerprintsVector;AtomTypesCount:SLogPAtomTypes:ArbitrarySize;16;Nume
    ricalValues;IDsAndValuesString;C1 C10 C11 C14 C18 C20 C21 C22 C5 CS F
    N11 N4 O10 O2 O9;5 1 1 1 14 4 2 1 2 2 1 1 1 1 3 1

    FingerprintsVector;AtomTypesCount:SLogPAtomTypes:FixedSize;67;OrderedN
    umericalValues;IDsAndValuesString;C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C
    12 C13 C14 C15 C16 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C27 CS N1 N
    2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12 N13 N14 NS O1 O2 O3 O4 O5 O6 O7 O8
    O9 O10 O11 O12 OS F Cl Br I Hal P S1 S2 S3 Me1 Me2;5 0 0 0 2 0 0 0 0 1
    1 0 0 1 0 0 0 14 0 4 2 1 0 0 0 0 0 2 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0...

    FingerprintsBitVector;AtomTypesBits:SLogPAtomTypes:FixedSize;67;Binary
    String;Ascending;10001000011001000101110000010001000000100000100000011
    0001000000000000000

    FingerprintsVector;AtomTypesCount:SYBYLAtomTypes:ArbitrarySize;9;Numer
    icalValues;IDsAndValuesString;C.2 C.3 C.ar F N.am N.ar O.2 O.3 O.co2;2
    9 22 1 1 1 1 2 2

    FingerprintsVector;AtomTypesCount:SYBYLAtomTypes:FixedSize;44;OrderedN
    umericalValues;IDsAndValuesString;C.3 C.2 C.1 C.ar C.cat N.3 N.2 N.1 N
    .ar N.am N.pl3 N.4 O.3 O.2 O.co2 S.3 S.2 S.o S.o2 P.3 F Cl Br I ANY HA
    L HET Li Na Mg Al Si K Ca Cr.th Cr.oh Mn Fe Co.oh Cu Zn Se Mo Sn;9 2 0
    22 0 0 0 0 1 1 0 0 2 1 2 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0

    FingerprintsBitVector;AtomTypesBits:SYBYLAtomTypes:FixedSize;44;Binary
    String;Ascending;110100001100111000001000000000000000000000000000

    FingerprintsVector;AtomTypesCount:TPSAAtomTypes:FixedSize;47;OrderedNu
    mericalValues;IDsAndValuesString;N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N1
    2 N13 N14 N15 N16 N17 N18 N19 N20 N21 N22 N23 N24 N25 N26 N O1 O2 O3 O
    4 O5 O6 O S1 S2 S3 S4 S5 S6 S7 S P1 P2 P3 P4 P;0 0 0 0 0 0 1 0 0 0 0 0
    0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

    FingerprintsBitVector;AtomTypesBits:TPSAAtomTypes:FixedSize;47;BinaryS
    tring;Ascending;000000100000000000001000000001100000000000000000

    FingerprintsVector;AtomTypesCount:UFFAtomTypes:ArbitrarySize;8;Numeric
    alValues;IDsAndValuesString;C_2 C_3 C_R F_ N_3 N_R O_2 O_3;2 9 22 1 1
    1 2 3

    FingerprintsVector;AtomTypesCount:UFFAtomTypes;124;OrderedNumerical
    Values;IDsAndValuesString;He4+4 Li Be3+2 B_3 B_2 C_3 C_R C_2 C_1 N_3 N_
    R N_2 N_1 O_3 O_3_z O_R O_2 O_1 F_ Ne4+4 Na Mg3+2 Al3 Si3 P_3+3 P_3+5 P
    _3+q S_3+2 S_3+4 S_3+6 S_R S_2 Cl Ar4+4 K_ Ca6+2 Sc3+3 Ti3+4 Ti6+4 V_3+
    ;0 0 0 0 0 12 0 3 0 3 0 1 0 2 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...

    FingerprintsVector;AtomTypesCount:UFFAtomTypes:FixedSize;124;OrderedNu
    mericalValues;IDsAndValuesString;He4+4 Li Be3+2 B_3 B_2 C_3 C_R C_2 C_
    1 N_3 N_R N_2 N_1 O_3 O_3_z O_R O_2 O_1 F_ Ne4+4 Na Mg3+2 Al3 Si3 P_3+
    3 P_3+5 P_3+q S_3+2 S_3+4 S_3+6 S_R S_2 Cl Ar4+4 K_ Ca6+2 Sc3+3 Ti...;
    0 0 0 0 0 9 22 2 0 1 1 0 0 3 0 0 2 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0...

    FingerprintsBitVector;AtomTypesBits:UFFAtomTypes:FixedSize;124;BinaryS
    tring;Ascending;000001110110010010100000000000000000000000000000000000
    0000000000000000000000000000000000000000000000000000000000000000000000

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

=item B<-a, --AtomIdentifierType> I<AtomicInvariantsAtomTypes | DREIDINGAtomTypes | EStateAtomTypes | FunctionalClassAtomTypes | MMFF94AtomTypes | SLogPAtomTypes | SYBYLAtomTypes | TPSAAtomTypes | UFFAtomTypes>

Specify atom identifier type to use for assignment of atom types to hydrogen and/or
non-hydrogen atoms during calculation of atom types fingerprints. Possible values in the
current release are: I<AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes,
TPSAAtomTypes, UFFAtomTypes>. Default value: I<AtomicInvariantsAtomTypes>.

=item B<--AtomicInvariantsToUse> I<"AtomicInvariant,AtomicInvariant...">

This value is used during I<AtomicInvariantsAtomTypes> value of B<a, --AtomIdentifierType>
option. It's a list of comma separated valid atomic invariant atom types.

Possible values for atomic invariants are: I<AS, X, BO,  LBO, SB, DB, TB,
H, Ar, RA, FC, MN, SM>. Default value: I<AS,X,BO,H,FC>.

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

I<AtomTypes::AtomicInvariantsAtomTypes> module is used to assign atomic invariant
atom types.

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


=item B<--AtomTypesSetToUse> I<ArbitrarySize | FixedSize>

Atom types set size to use during generation of atom types fingerprints.

Possible values for I<AtomTypesCount> values of B<-m, --mode> option: I<ArbitrarySize |
FixedSize>; Default value: I<ArbitrarySize>.

Possible values for I<AtomTypesBits> value of B<-m, --mode> option: I<FixedSize>;
Default value: I<FixedSize>.

I<FixedSize> value is not supported for I<AtomicInvariantsAtomTypes> value of
B<-a, --AtomIdentifierType> option.

I<ArbitrarySize> corresponds to only atom types detected in molecule; I<FixedSize> corresponds
to fixed number of previously defined atom types for specified B<-a, --AtomIdentifierType>.

=item B<--BitsOrder> I<Ascending | Descending>

Bits order to use during generation of fingerprints bit-vector string for I<AtomTypesBits> value of
=item B<--BitsOrder> I<Ascending | Descending>

Bits order to use during generation of fingerprints bit-vector string for I<AtomTypesBits> value of
B<-m, --mode> option. Possible values: I<Ascending, Descending>. Default: I<Ascending>.

I<Ascending> bit order which corresponds to first bit in each byte as the lowest bit as
opposed to the highest bit.

Internally, bits are stored in I<Ascending> order using Perl vec function. Regardless
of machine order, big-endian or little-endian, vec function always considers first
string byte as the lowest byte and first bit within each byte as the lowest bit.

=item B<-b, --BitStringFormat> I<BinaryString | HexadecimalString>

Format of fingerprints bit-vector string data in output SD, FP or CSV/TSV text file(s) specified by
B<--output> used during I<AtomTypesBits> value of B<-m, --mode> option. Possible
values: I<BinaryString, HexadecimalString>. Default value: I<BinaryString>.

I<BinaryString> corresponds to an ASCII string containing 1s and 0s. I<HexadecimalString>
contains bit values in ASCII hexadecimal format.

Examples:

    FingerprintsBitVector;AtomTypesBits:DREIDINGAtomTypes;34;BinaryString;
    Ascending;0010101010101000000000000000000000000000

    FingerprintsBitVector;AtomTypesBits:MMFF94AtomTypes;171;BinaryString;
    Ascending;1000010101000000000001100000000000000001010000101000000000000
    00000000000000000000000000000000000001000000000000000000000000000000000
    0000000000000000000000000000000000000000000

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

=item B<--DataFields> I<"FieldLabel1,FieldLabel2,...">

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

=item B<-f, --Filter> I<Yes | No>

Specify whether to check and filter compound data in SDFile(s). Possible values: I<Yes or No>.
Default value: I<Yes>.

By default, compound data is checked before calculating fingerprints and compounds containing
atom data corresponding to non-element symbols or no atom data are ignored.

=item B<--FingerprintsLabelMode> I<FingerprintsLabelOnly | FingerprintsLabelWithIDs>

Specify how fingerprints label is generated in conjunction with B<--FingerprintsLabel> option value:
use fingerprints label generated only by B<--FingerprintsLabel> option value or append atom type
value IDs to B<--FingerprintsLabel> option value.

Possible values: I<FingerprintsLabelOnly | FingerprintsLabelWithIDs>. Default value:
I<FingerprintsLabelOnly>.

This option is only used for I<FixedSize> value of B<-e, --AtomTypesSetToUse> option during
generation of I<AtomTypesCount> fingerprints and ignored for I<AtomTypesBits>.

Atom type IDs appended to B<--FingerprintsLabel> value during I<FingerprintsLabelWithIDs>
values of B<--FingerprintsLabelMode> correspond to fixed number of previously defined
atom types.

=item B<--FingerprintsLabel> I<text>

SD data label or text file column label to use for fingerprints string in output SD or
CSV/TSV text file(s) specified by B<--output>. Default value: I<AtomTypesFingerprints>.

=item B<-h, --help>

Print this help message.

=item B<-i, --IgnoreHydrogens> I<Yes | No>

Ignore hydrogens during fingerprints generation. Possible values: I<Yes or No>.
Default value: I<Yes>.

For I<yes> value of B<-i, --IgnoreHydrogens>, any explicit hydrogens are also used for
generation of atom type fingerprints; implicit hydrogens are still ignored.

=item B<-k, --KeepLargestComponent> I<Yes | No>

Generate fingerprints for only the largest component in molecule. Possible values:
I<Yes or No>. Default value: I<Yes>.

For molecules containing multiple connected components, fingerprints can be generated
in two different ways: use all connected components or just the largest connected
component. By default, all atoms except for the largest connected component are
deleted before generation of fingerprints.

=item B<-m, --mode> I<AtomTypesCount | AtomTypesBits>

Specify type of atom types fingerprints to generate for molecules in I<SDFile(s)>.
Possible values: I<AtomTypesCount or AtomTypesBits>. Default value: I<AtomTypesCount>.

For I<AtomTypesCount> values of B<-m, --mode> option, a fingerprint vector string is generated.
The vector string corresponding to I<AtomTypesCount> contains count of atom types.

For I<AtomTypesBits> value of B<-m, --mode> option, a fingerprint bit-vector string containing
zeros and ones indicating presence or absence of atom types is generated.

For I<AtomTypesCount> atom types fingerprints, two types of atom types set size can be specified
using B<-a, --AtomTypesSetToUse> option: I<ArbitrarySize or FixedSize>. I<ArbitrarySize> corrresponds
to only atom types detected in molecule; I<FixedSize> corresponds to fixed number of atom types
previously defined.

For I<AtomTypesBits> atom types fingeprints, only I<FixedSize> is allowed.

Combination of B<-m, --Mode> and B<--AtomTypesSetToUse> along with B<-a, --AtomtomIdentifierType>
allows generation of following different atom types fingerprints:

    Mode                  AtomIdentifierType           AtomTypesSetToUse

    AtomTypesCount        AtomicInvariantsAtomTypes    ArbitrarySize [ Default ]

    AtomTypesCount        DREIDINGAtomTypes            ArbitrarySize
    AtomTypesCount        DREIDINGAtomTypes            FixedSize
    AtomTypesBits         DREIDINGAtomTypes            FixedSize

    AtomTypesCount        EStateAtomTypes              ArbitrarySize
    AtomTypesCount        EStateAtomTypes              FixedSize
    AtomTypesBits         EStateAtomTypes              FixedSize

    AtomTypesCount        FunctionalClassAtomTypes    ArbitrarySize

    AtomTypesCount        MMFF94AtomTypes              ArbitrarySize
    AtomTypesCount        MMFF94AtomTypes              FixedSize
    AtomTypesBits         MMFF94AtomTypes              FixedSize

    AtomTypesCount        SLogPAtomTypes               ArbitrarySize
    AtomTypesCount        SLogPAtomTypes               FixedSize
    AtomTypesBits         SLogPAtomTypes               FixedSize

    AtomTypesCount        SYBYLAtomTypes               ArbitrarySize
    AtomTypesCount        SYBYLAtomTypes               FixedSize
    AtomTypesBits         SYBYLAtomTypes               FixedSize

    AtomTypesCount        TPSAAtomTypes                 FixedSize
    AtomTypesBits         TPSAAtomTypes                 FixedSize

    AtomTypesCount        UFFAtomTypes                 ArbitrarySize
    AtomTypesCount        UFFAtomTypes                 FixedSize
    AtomTypesBits         UFFAtomTypes                 FixedSize

The default is to generate I<AtomicInvariantAtomTypes> fingeprints corresponding to I<ArbitrarySize> as
value of B<--AtomTypesSetToUse> option.

=item B<--OutDelim> I<comma | tab | semicolon>

Delimiter for output CSV/TSV text file(s). Possible values: I<comma, tab, or semicolon>
Default value: I<comma>.

=item B<--output> I<SD | FP | text | all>

Type of output files to generate. Possible values: I<SD, FP, text, or all>. Default value: I<text>.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-q, --quote> I<Yes | No>

Put quote around column values in output CSV/TSV text file(s). Possible values:
I<Yes or No>. Default value: I<Yes>.

=item B<-r, --root> I<RootName>

New file name is generated using the root: <Root>.<Ext>. Default for new file
names: <SDFileName><AtomTypesFP>.<Ext>. The file type determines <Ext> value.
The sdf, fpf, csv, and tsv <Ext> values are used for SD, FP, comma/semicolon, and tab
delimited text files, respectively.This option is ignored for multiple input files.

=item B<-v, --VectorStringFormat> I<ValuesString | IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>

Format of fingerprints vector string data in output SD, FP or CSV/TSV text file(s) specified by
B<--output> used during <AtomTypesCount> value of B<-m, --mode> option. Possible values:
I<ValuesString, IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString |
ValuesAndIDsPairsString>.

Default value during I<ArbitrarySize> value of B<-e, --AtomTypesSetToUse>
option: I<IDsAndValuesString>. Default value during I<FixedSize> value of
B<-e, --AtomTypesSetToUse> option: I<ValuesString>.

Example of I<SD> file containing atom types fingerprints string data:

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

    >  <AtomTypesFingerprints>
    FingerprintsVector;AtomTypesCount:AtomicInvariantsAtomTypes:ArbitrarySi
    ze;10;NumericalValues;IDsAndValuesString;C.X1.BO1.H3 C.X2.BO2.H2 C.X2.B
    O3.H1 C.X3.BO3.H1 C.X3.BO4 F.X1.BO1 N.X2.BO2.H1 N.X3.BO3 O.X1.BO1.H1 O.
    X1.BO2;2 4 14 3 10 1 1 1 3 2

    $$$$
    ... ...
    ... ...

Example of I<FP> file containing atom types fingerprints string data:

    #
    # Package = MayaChemTools 7.4
    # Release Date = Oct 21, 2010
    #
    # TimeStamp = Fri Mar 11 14:28:07 2011
    #
    # FingerprintsStringType = FingerprintsVector
    #
    # Description = AtomTypesCount:AtomicInvariantsAtomTypes:ArbitrarySize
    # VectorStringFormat = IDsAndValuesString
    # VectorValuesType = NumericalValues
    #
    Cmpd1 10;C.X1.BO1.H3 C.X2.BO2.H2 C.X2.BO3.H1 C.X3.BO3.H1 C.X3.BO4 F...
    Cmpd2 9;C.X1.BO1.H3 C.X2.BO2.H2 C.X3.BO3.H1 C.X3.BO4 N.X1.BO1.H2 N....
    ... ...
    ... ..

Example of CSV I<Text> file atom types containing fingerprints string data:

    "CompoundID","AtomTypesFingerprints"
    "Cmpd1","FingerprintsVector;AtomTypesCount:AtomicInvariantsAtomTypes:Ar
    bitrarySize;10;NumericalValues;IDsAndValuesString;C.X1.BO1.H3 C.X2.BO2.
    H2 C.X2.BO3.H1 C.X3.BO3.H1 C.X3.BO4 F.X1.BO1 N.X2.BO2.H1 N.X3.BO3 O.X1.
    BO1.H1 O.X1.BO2;2 4 14 3 10 1 1 1 3 2"
    O.X1.BO2;3 3 6 3 1 1 2 2 2"
    ... ...
    ... ...

Examples:

    FingerprintsVector;AtomTypesCount:EStateAtomTypes:ArbitrarySize;11;Num
    ericalValues;IDsAndValuesString;aaCH aasC aasN dO dssC sCH3 sF sOH ssC
    H2 ssNH sssCH;14 8 1 2 2 2 1 3 4 1 3

    FingerprintsVector;AtomTypesCount:SYBYLAtomTypes:ArbitrarySize;9;Numer
    icalValues;IDsAndValuesString;C.2 C.3 C.ar F N.am N.ar O.2 O.3 O.co2;2
    9 22 1 1 1 1 2 2

    FingerprintsVector;AtomTypesCount:SYBYLAtomTypes:FixedSize;44;OrderedN
    umericalValues;IDsAndValuesString;C.3 C.2 C.1 C.ar C.cat N.3 N.2 N.1 N
    .ar N.am N.pl3 N.4 O.3 O.2 O.co2 S.3 S.2 S.o S.o2 P.3 F Cl Br I ANY HA
    L HET Li Na Mg Al Si K Ca Cr.th Cr.oh Mn Fe Co.oh Cu Zn Se Mo Sn;9 2 0
    22 0 0 0 0 1 1 0 0 2 1 2 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0

=item B<-w, --WorkingDir> I<DirName>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To generate atomic invariants atom types count fingerprints of arbitrary size in vector
string format and create a SampleATFP.csv file containing sequential compound IDs along
with fingerprints vector strings data, type:

    % AtomTypesFingerprints.pl -r SampleATFP -o Sample.sdf

To generate functional class atom types count fingerprints of arbitrary size in vector
string format and create a SampleATFP.csv file containing sequential compound IDs along
with fingerprints vector strings data, type:

    % AtomTypesFingerprints.pl -m AtomTypesCount -a FunctionalClassAtomTypes
      -r SampleATFP -o Sample.sdf

To generate E-state atom types count fingerprints of arbitrary size in vector string
format and create a SampleATFP.csv file containing sequential compound IDs along
with fingerprints vector strings data, type:

    % AtomTypesFingerprints.pl -m AtomTypesCount -a EStateAtomTypes
      --AtomTypesSetToUse ArbitrarySize -r SampleATFP -o Sample.sdf

To generate E-state atom types count fingerprints of fixed size in vector string
with IDsAndValues format and create a SampleATFP.csv file containing sequential
compound IDs along with fingerprints vector strings data, type:

    % AtomTypesFingerprints.pl -m AtomTypesCount -a EStateAtomTypes
      --AtomTypesSetToUse FixedSize -v IDsAndValuesString
      -r SampleATFP -o Sample.sdf

To generate E-state atom types bits fingerprints of fixed size in bit-vector string
format and create a SampleATFP.csv file containing sequential compound IDs along
with fingerprints vector strings data, type:

    % AtomTypesFingerprints.pl -m AtomTypesBits -a EStateAtomTypes
      --AtomTypesSetToUse FixedSize -r SampleATFP -o Sample.sdf

To generate MMFF94 atom types count fingerprints of arbitrary size in vector string
format and create a SampleATFP.csv file containing sequential compound IDs along
with fingerprints vector strings data, type:

    % AtomTypesFingerprints.pl -m AtomTypesCount -a MMFF94AtomTypes
      --AtomTypesSetToUse ArbitrarySize -r SampleATFP -o Sample.sdf

To generate MMFF94 atom types count fingerprints of fixed size in vector string
format and create a SampleATFP.csv file containing sequential compound IDs along
with fingerprints vector strings data, type:

    % AtomTypesFingerprints.pl -m AtomTypesCount -a MMFF94AtomTypes
      --AtomTypesSetToUse FixedSize -r SampleATFP -o Sample.sdf

To generate MMFF94 atom types count fingerprints of fixed size in vector string
with IDsAndValues format and create a SampleATFP.csv file containing sequential
compound IDs along with fingerprints vector strings data, type:

    % AtomTypesFingerprints.pl -m AtomTypesCount -a MMFF94AtomTypes
      --AtomTypesSetToUse FixedSize -v IDsAndValuesString
      -r SampleATFP -o Sample.sdf

To generate MMFF94 atom types bits fingerprints of fixed size in bit-vector string
format and create a SampleATFP.csv file containing sequential compound IDs along
with fingerprints vector strings data, type:

    % AtomTypesFingerprints.pl -m AtomTypesBits -a MMFF94AtomTypes
      --AtomTypesSetToUse FixedSize -r SampleATFP -o Sample.sdf

To generate MMFF94 atom types count fingerprints of arbitrary size in vector string
format and create a SampleATFP.csv file containing compound ID from molecule
name line along with fingerprints vector strings data, type

    % AtomTypesFingerprints.pl -m AtomTypesCount -a MMFF94AtomTypes
      --DataFieldsMode CompoundID --CompoundIDMode MolName
      -r SampleATFP -o Sample.sdf

To generate MMFF94 atom types count fingerprints of arbitrary size in vector string
format and create a SampleATFP.csv file containing compound IDs using specified
data field along with fingerprints vector strings data, type:

    % AtomTypesFingerprints.pl -m AtomTypesCount -a MMFF94AtomTypes
      --DataFieldsMode CompoundID --CompoundIDMode DataField --CompoundID
      Mol_ID -r SampleATFP -o Sample.sdf

To generate MMFF94 atom types count fingerprints of arbitrary size in vector string
format and create a SampleATFP.csv file containing compound ID using combination
of molecule name line and an explicit compound prefix along with fingerprints vector
strings data, type:

    % AtomTypesFingerprints.pl -m AtomTypesCount -a MMFF94AtomTypes
      --DataFieldsMode CompoundID --CompoundIDMode MolnameOrLabelPrefix
      --CompoundID Cmpd --CompoundIDLabel MolID -r SampleATFP -o Sample.sdf

To generate MMFF94 atom types count fingerprints of arbitrary size in vector string
format and create a SampleATFP.csv file containing specific data fields columns along
with fingerprints vector strings data, type:

    % AtomTypesFingerprints.pl -m AtomTypesCount -a MMFF94AtomTypes
      --DataFieldsMode Specify --DataFields Mol_ID -r SampleATFP
      -o Sample.sdf

To generate MMFF94 atom types count fingerprints of arbitrary size in vector string
format and create a SampleATFP.csv file containing common data fields columns along
with fingerprints vector strings data, type:

    % AtomTypesFingerprints.pl -m AtomTypesCount -a MMFF94AtomTypes
     --DataFieldsMode Common -r SampleATFP -o Sample.sdf

To generate MMFF94 atom types count fingerprints of arbitrary size in vector string
format and create SampleATFP.sdf,  SampleATFP.fpf and  SampleATFP.csv files containing
all data fields columns in CSV file along with fingerprints vector strings data, type:

    % AtomTypesFingerprints.pl -m AtomTypesCount -a MMFF94AtomTypes
     --DataFieldsMode All --output all -r SampleATFP -o Sample.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoFingerprintsFiles.pl, SimilarityMatricesFingerprints.pl, AtomNeighborhoodsFingerprints.pl,
ExtendedConnectivityFingerprints.pl, MACCSKeysFingeprints.pl, PathLengthFingerprints.pl,
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
