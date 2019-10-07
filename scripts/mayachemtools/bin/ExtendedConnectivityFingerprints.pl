#!/usr/bin/perl -w
#
# File: ExtendedConnectivityFingerprints.pl
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
use Fingerprints::ExtendedConnectivityFingerprints;

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
    GenerateExtendedConnectivityFingerprints($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Generate fingerprints for a SD file...
#
sub GenerateExtendedConnectivityFingerprints {
  my($FileIndex) = @_;
  my($CmpdCount, $IgnoredCmpdCount, $SDFile, $MoleculeFileIO, $Molecule, $ExtendedConnectivityFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO);

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

    $ExtendedConnectivityFingerprints = GenerateMoleculeFingerprints($Molecule);
    if (!$ExtendedConnectivityFingerprints) {
      $IgnoredCmpdCount++;
      ProcessIgnoredCompound('FingerprintsGenerationFailed', $CmpdCount, $Molecule);
      next COMPOUND;
    }

    WriteDataToOutputFiles($FileIndex, $CmpdCount, $Molecule, $ExtendedConnectivityFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO);
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
  if ($OptionsInfo{Mode} =~ /^(ExtendedConnectivity|ExtendedConnectivityCount)$/i) {
    %FingerprintsFileIOParams = ('Mode' => 'Write', 'Overwrite' => $OptionsInfo{OverwriteFiles}, 'FingerprintsStringMode' => 'FingerprintsVectorString', 'VectorStringFormat' => $OptionsInfo{VectorStringFormat});
  }
  elsif ($OptionsInfo{Mode} =~ /^ExtendedConnectivityBits$/i) {
    %FingerprintsFileIOParams = ('Mode' => 'Write', 'Overwrite' => $OptionsInfo{OverwriteFiles}, 'FingerprintsStringMode' => 'FingerprintsBitVectorString', 'BitStringFormat' => $OptionsInfo{BitStringFormat}, 'BitsOrder' => $OptionsInfo{BitsOrder});
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
  my($FileIndex, $CmpdCount, $Molecule, $ExtendedConnectivityFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO) = @_;
  my($DataFieldLabelAndValuesRef);

  $DataFieldLabelAndValuesRef = undef;
  if ($NewFPTextFileIO || $NewFPFileIO) {
    $DataFieldLabelAndValuesRef = $Molecule->GetDataFieldLabelAndValues();
  }

  if ($NewFPSDFileIO) {
    my($CmpdString);

    $CmpdString = $Molecule->GetInputMoleculeString();
    $NewFPSDFileIO->WriteFingerprints($ExtendedConnectivityFingerprints, $CmpdString);
  }

  if ($NewFPTextFileIO) {
    my($ColValuesRef);

    $ColValuesRef = SetupFPTextFileCoulmnValues($FileIndex, $CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPTextFileIO->WriteFingerprints($ExtendedConnectivityFingerprints, $ColValuesRef);
  }

  if ($NewFPFileIO) {
    my($CompoundID);

    $CompoundID = SetupCmpdIDForOutputFiles($CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPFileIO->WriteFingerprints($ExtendedConnectivityFingerprints, $CompoundID);
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
  my($ExtendedConnectivityFingerprints);

  if ($OptionsInfo{KeepLargestComponent}) {
    $Molecule->KeepLargestComponent();
  }
  if (!$Molecule->DetectRings()) {
    return undef;
  }
  $Molecule->SetAromaticityModel($OptionsInfo{AromaticityModel});
  $Molecule->DetectAromaticity();

  $ExtendedConnectivityFingerprints = undef;
  if ($OptionsInfo{Mode} =~ /^(ExtendedConnectivity|ExtendedConnectivityCount)$/i ) {
    $ExtendedConnectivityFingerprints = new Fingerprints::ExtendedConnectivityFingerprints('Type' => $OptionsInfo{Mode}, 'Molecule' => $Molecule, 'NeighborhoodRadius' => $OptionsInfo{NeighborhoodRadius}, 'AtomIdentifierType' => $OptionsInfo{AtomIdentifierType});
  }
  elsif ($OptionsInfo{Mode} =~ /^ExtendedConnectivityBits$/i) {
    $ExtendedConnectivityFingerprints = new Fingerprints::ExtendedConnectivityFingerprints('Type' => $OptionsInfo{Mode}, 'Molecule' => $Molecule, 'NeighborhoodRadius' => $OptionsInfo{NeighborhoodRadius}, 'AtomIdentifierType' => $OptionsInfo{AtomIdentifierType}, 'Size' => $OptionsInfo{Size}, 'UsePerlCoreRandom' => $OptionsInfo{UsePerlCoreRandom});
  }
  else {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: ExtendedConnectivity,  ExtendedConnectivityCount or ExtendedConnectivityBits\n";
  }
  SetAtomIdentifierTypeValuesToUse($ExtendedConnectivityFingerprints);

  # Generate fingerprints...
  $ExtendedConnectivityFingerprints->GenerateFingerprints();

  # Make sure fingerprints generation is successful...
  if (!$ExtendedConnectivityFingerprints->IsFingerprintsGenerationSuccessful()) {
    return undef;
  }

  return $ExtendedConnectivityFingerprints;
}

# Set atom identifier type to use for generating fingerprints...
#
sub SetAtomIdentifierTypeValuesToUse {
  my($ExtendedConnectivityFingerprints) = @_;

  if ($OptionsInfo{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
    $ExtendedConnectivityFingerprints->SetAtomicInvariantsToUse(\@{$OptionsInfo{AtomicInvariantsToUse}});
  }
  elsif ($OptionsInfo{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
    $ExtendedConnectivityFingerprints->SetFunctionalClassesToUse(\@{$OptionsInfo{FunctionalClassesToUse}});
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
      $OutFileRoot = "${FileName}ExtendedConnectivityFP";
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

  ProcessAtomIdentifierTypeOptions();

  $OptionsInfo{AromaticityModel} = $Options{aromaticitymodel};

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

  $OptionsInfo{FingerprintsLabel} = $Options{fingerprintslabel} ? $Options{fingerprintslabel} : 'ExtendedConnectivityFingerprints';

  $OptionsInfo{Filter} = ($Options{filter} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{KeepLargestComponent} = ($Options{keeplargestcomponent} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{NeighborhoodRadius} = $Options{neighborhoodradius};

  $OptionsInfo{UsePerlCoreRandom} = ($Options{useperlcorerandom} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{Output} = $Options{output};
  $OptionsInfo{SDOutput} = ($Options{output} =~ /^(SD|All)$/i) ? 1 : 0;
  $OptionsInfo{FPOutput} = ($Options{output} =~ /^(FP|All)$/i) ? 1 : 0;
  $OptionsInfo{TextOutput} = ($Options{output} =~ /^(Text|All)$/i) ? 1 : 0;

  $OptionsInfo{OutDelim} = $Options{outdelim};
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? 1 : 0;
  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : 0;

  my($Size, $MinSize, $MaxSize);
  $MinSize = 32;
  $MaxSize = 2**32;
  $Size = $Options{size};
  if (!(IsPositiveInteger($Size) && $Size >= $MinSize && $Size <= $MaxSize && IsNumberPowerOfNumber($Size, 2))) {
    die "Error: Invalid size value, $Size, for \"-s, --size\" option. Allowed values: power of 2, >= minimum size of $MinSize, and <= maximum size of $MaxSize.\n";
  }
  $OptionsInfo{Size} = $Size;

  # Setup default vector string format...
  #
  my($VectorStringFormat);
  $VectorStringFormat = '';
  if ($Options{vectorstringformat}) {
    $VectorStringFormat = $Options{vectorstringformat};
  }
  else {
    $VectorStringFormat = ($Options{mode} =~ /^ExtendedConnectivity$/) ? "ValuesString" : "IDsAndValuesString";
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
  elsif ($Options{atomidentifiertype} =~ /^(DREIDINGAtomTypes|EStateAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
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
  $Options{atomicinvariantstouse} = 'AS,X,BO,H,FC,MN';
  $Options{functionalclassestouse} = 'HBD,HBA,PI,NI,Ar,Hal';

  $Options{bitsorder} = 'Ascending';
  $Options{bitstringformat} = 'HexadecimalString';

  $Options{compoundidmode} = 'LabelPrefix';
  $Options{compoundidlabel} = 'CompoundID';
  $Options{datafieldsmode} = 'CompoundID';

  $Options{filter} = 'Yes';

  $Options{keeplargestcomponent} = 'Yes';

  $Options{mode} = 'ExtendedConnectivity';

  $Options{neighborhoodradius} = 2;

  $Options{useperlcorerandom} = 'yes';

  $Options{output} = 'text';
  $Options{outdelim} = 'comma';
  $Options{quote} = 'yes';

  $Options{size} = 1024;

  $Options{vectorstringformat} = '';

  if (!GetOptions(\%Options, "aromaticitymodel=s", "atomidentifiertype|a=s", "atomicinvariantstouse=s", "functionalclassestouse=s", "bitsorder=s", "bitstringformat|b=s", "compoundid=s", "compoundidlabel=s", "compoundidmode=s", "datafields=s", "datafieldsmode|d=s", "filter|f=s", "fingerprintslabel=s",  "help|h", "keeplargestcomponent|k=s",  "mode|m=s", "neighborhoodradius|n=s", "outdelim=s", "output=s", "overwrite|o", "quote|q=s", "root|r=s", "size|s=i", "useperlcorerandom=s", "vectorstringformat|v=s", "workingdir|w=s")) {
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
  if ($Options{atomidentifiertype} !~ /^(AtomicInvariantsAtomTypes|FunctionalClassAtomTypes|DREIDINGAtomTypes|EStateAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
    die "Error: The value specified, $Options{atomidentifiertype}, for option \"-a, --AtomIdentifierType\" is not valid. Supported atom identifier types in current release of MayaChemTools: AtomicInvariantsAtomTypes, FunctionalClassAtomTypes, DREIDINGAtomTypes, EStateAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes\n";
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
  if ($Options{keeplargestcomponent} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{keeplargestcomponent}, for option \"-k, --KeepLargestComponent\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{mode} !~ /^(ExtendedConnectivity|ExtendedConnectivityCount|ExtendedConnectivityBits)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: ExtendedConnectivity, ExtendedConnecticityCount, or ExtendedConnectivityBits\n";
  }
  if (!(IsInteger($Options{neighborhoodradius}) && ($Options{neighborhoodradius} >= 0))) {
    die "Error: The value specified, $Options{neighborhoodradius}, for option \"-n, --NeighborhoodRadius\" is not valid. Allowed values: >= 0 \n";
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
  if (!IsPositiveInteger($Options{size})) {
    die "Error: The value specified, $Options{size}, for option \"-s, --size\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{outdelim} =~ /semicolon/i && $Options{quote} =~ /^No$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not allowed with, semicolon value of \"--outdelim\" option: Fingerprints string use semicolon as delimiter for various data fields and must be quoted.\n";
  }
  if ($Options{useperlcorerandom} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{useperlcorerandom}, for option \"--UsePerlCoreRandom\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{vectorstringformat} && $Options{vectorstringformat} !~ /^(ValuesString|IDsAndValuesString|IDsAndValuesPairsString|ValuesAndIDsString|ValuesAndIDsPairsString)$/i) {
    die "Error: The value specified, $Options{vectorstringformat}, for option \"-v, --VectorStringFormat\" is not valid. Allowed values: ValuesString, IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString or ValuesAndIDsPairsString\n";
  }
}

__END__

=head1 NAME

ExtendedConnectivityFingerprints.pl - Generate extended connectivity fingerprints for SD files

=head1 SYNOPSIS

ExtendedConnectivityFingerprints.pl SDFile(s)...

ExtendedConnectivityFingerprints.pl [B<--AromaticityModel> I<AromaticityModelType>]
[B<-a, --AtomIdentifierType> I<AtomicInvariantsAtomTypes>]
[B<--AtomicInvariantsToUse> I<"AtomicInvariant,AtomicInvariant...">]
[B<--FunctionalClassesToUse> I<"FunctionalClass1,FunctionalClass2...">]
[B<--BitsOrder> I<Ascending | Descending>] [B<-b, --BitStringFormat> I<BinaryString | HexadecimalString>]
[B<--CompoundID> I<DataFieldName or LabelPrefixString>] [B<--CompoundIDLabel> I<text>]
[B<--CompoundIDMode>] [B<--DataFields> I<"FieldLabel1,FieldLabel2,...">]
[B<-d, --DataFieldsMode> I<All | Common | Specify | CompoundID>]  [B<-f, --Filter> I<Yes | No>]
[B<--FingerprintsLabel> I<text>] [B<-h, --help>] [B<-k, --KeepLargestComponent> I<Yes | No>]
[B<-m, --mode> I<ExtendedConnectivity | ExtendedConnecticityCount | ExtendedConnecticityBits>]
[B<-n, --NeighborhoodRadius> I<number>] [B<--OutDelim> I<comma | tab | semicolon>] [B<--output> I<SD | FP | text | all>]
[B<-o, --overwrite>] [B<-q, --quote> I<Yes | No>] [B<-r, --root> I<RootName>] [B<-s, --size> I<number>]
[B<--UsePerlCoreRandom> I<Yes | No>]
[B<-v, --VectorStringFormat> I<IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>]
[B<-w, --WorkingDir> dirname] SDFile(s)...

=head1 DESCRIPTION

Generate extended connectivity fingerprints [ Ref 48, Ref 52 ] for I<SDFile(s)> and create appropriate
SD, FP or CSV/TSV text file(s) containing fingerprints vector strings corresponding to molecular fingerprints.

Multiple SDFile names are separated by spaces. The valid file extensions are I<.sdf>
and I<.sd>. All other file names are ignored. All the SD files in a current directory
can be specified either by I<*.sdf> or the current directory name.

The current release of MayaChemTools supports generation of extended connectivity fingerprints
corresponding to following B<-a, --AtomIdentifierTypes>:

    AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
    FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
    SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes

Based on values specified for B<-a, --AtomIdentifierType>, B<--AtomicInvariantsToUse>
and B<--FunctionalClassesToUse>, initial atom types are assigned to all non-hydrogen atoms in
a molecule and these atom types strings are converted into initial atom identifier integers using
B<TextUtil::HashCode> function. The duplicate atom identifiers are removed.

For B<-n, --NeighborhoodRadius> value of I<0>, the initial set of unique atom identifiers comprises
the molecule fingerprints. Otherwise, atom neighborhoods are generated for each non-hydrogen
atom up to specified B<-n, --NeighborhoodRadius> value. For each non-hydrogen central atom
at a specific radius, its neighbors at next radius level along with their bond orders and previously
calculated atom identifiers are collected which in turn are used to generate a new integer
atom identifier; the bond orders and atom identifier pairs list is first sorted by bond order
followed by atom identifiers to make these values graph invariant.

After integer atom identifiers have been generated for all non-hydrogen atoms at all specified
neighborhood radii, the duplicate integer atom identifiers corresponding to same hash code
value generated using B<TextUtil::HashCode> are tracked by keeping the atom identifiers at
lower radius. Additionally, all structurally duplicate integer atom identifiers at each specified
radius are also tracked by identifying equivalent atoms and bonds corresponding to substructures
used for generating atom identifier and keeping integer atom identifier with lowest value.

For I<ExtendedConnnectivity> value of fingerprints B<-m, --mode>, the duplicate identifiers are
removed from the list and the unique atom identifiers constitute the extended connectivity
fingerprints of a molecule.

For I<ExtendedConnnectivityCount> value of fingerprints B<-m, --mode>, the occurrence of each
unique atom identifiers appears is counted and the unique atom identifiers along with their
count constitute the extended connectivity fingerprints of a molecule.

For I<ExtendedConnectivityBits> value of fingerprints B<-m, --mode>, the unique atom identifiers
are used as a random number seed to generate a random integer value between 0 and B<--Size> which
in turn is used to set corresponding bits in the fingerprint bit-vector string.

Example of I<SD> file containing extended connectivity fingerprints string data:

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

    >  <ExtendedConnectivityFingerprints>
    FingerprintsVector;ExtendedConnectivity:AtomicInvariantsAtomTypes:Radiu
    s2;60;AlphaNumericalValues;ValuesString;73555770 333564680 352413391 66
    6191900 1001270906 1371674323 1481469939 1977749791 2006158649 21414087
    99 49532520 64643108 79385615 96062769 273726379 564565671 855141035 90
    6706094 988546669 1018231313 1032696425 1197507444 1331250018 133853...

    $$$$
    ... ...
    ... ...

Example of I<FP> file containing extended connectivity fingerprints string data:

    #
    # Package = MayaChemTools 7.4
    # Release Date = Oct 21, 2010
    #
    # TimeStamp = Fri Mar 11 14:43:57 2011
    #
    # FingerprintsStringType = FingerprintsVector
    #
    # Description = ExtendedConnectivity:AtomicInvariantsAtomTypes:Radius2
    # VectorStringFormat = ValuesString
    # VectorValuesType = AlphaNumericalValues
    #
    Cmpd1 60;73555770 333564680 352413391 666191900 1001270906 137167432...
    Cmpd2 41;73555770 333564680 666191900 1142173602 1363635752 14814699...
    ... ...
    ... ..

Example of CSV I<Text> file containing extended connectivity fingerprints string data:

    "CompoundID","ExtendedConnectivityFingerprints"
    "Cmpd1","FingerprintsVector;ExtendedConnectivity:AtomicInvariantsAtomTy
    pes:Radius2;60;AlphaNumericalValues;ValuesString;73555770 333564680 352
    413391 666191900 1001270906 1371674323 1481469939 1977749791 2006158649
    2141408799 49532520 64643108 79385615 96062769 273726379 564565671 8551
    41035 906706094 988546669 1018231313 1032696425 1197507444 13312500..."
    ... ...
    ... ...

The current release of MayaChemTools generates the following types of extended connectivity
fingerprints vector strings:

    FingerprintsVector;ExtendedConnectivity:AtomicInvariantsAtomTypes:Radi
    us2;60;AlphaNumericalValues;ValuesString;73555770 333564680 352413391
    666191900 1001270906 1371674323 1481469939 1977749791 2006158649 21414
    08799 49532520 64643108 79385615 96062769 273726379 564565671 85514103
    5 906706094 988546669 1018231313 1032696425 1197507444 1331250018 1338
    532734 1455473691 1607485225 1609687129 1631614296 1670251330 17303...

    FingerprintsVector;ExtendedConnectivityCount:AtomicInvariantsAtomTypes
    :Radius2;60;NumericalValues;IDsAndValuesString;73555770 333564680 3524
    13391 666191900 1001270906 1371674323 1481469939 1977749791 2006158649
    2141408799 49532520 64643108 79385615 96062769 273726379 564565671...;
    3 2 1 1 14 1 2 10 4 3 1 1 1 1 2 1 2 1 1 1 2 3 1 1 2 1 3 3 8 2 2 2 6 2
    1 2 1 1 2 1 1 1 2 1 1 2 1 2 1 1 1 1 1 1 1 1 1 2 1 1

    FingerprintsBitVector;ExtendedConnectivityBits:AtomicInvariantsAtomTyp
    es:Radius2;1024;BinaryString;Ascending;0000000000000000000000000000100
    0000000001010000000110000011000000000000100000000000000000000000100001
    1000000110000000000000000000000000010011000000000000000000000000010000
    0000000000000000000000000010000000000000000001000000000000000000000000
    0000000000010000100001000000000000101000000000000000100000000000000...

    FingerprintsBitVector;ExtendedConnectivityBits:AtomicInvariantsAtomTyp
    es:Radius2;1024;HexadecimalString;Ascending;000000010050c0600800000803
    0300000091000004000000020000100000000124008200020000000040020000000000
    2080000000820040010020000000008040000000000080001000000000400000000000
    4040000090000061010000000800200000000000001400000000020080000000000020
    00008020200000408000

    FingerprintsVector;ExtendedConnectivity:FunctionalClassAtomTypes:Radiu
    s2;57;AlphaNumericalValues;ValuesString;24769214 508787397 850393286 8
    62102353 981185303 1231636850 1649386610 1941540674 263599683 32920567
    1 571109041 639579325 683993318 723853089 810600886 885767127 90326012
    7 958841485 981022393 1126908698 1152248391 1317567065 1421489994 1455
    632544 1557272891 1826413669 1983319256 2015750777 2029559552 20404...

    FingerprintsVector;ExtendedConnectivityCount:FunctionalClassAtomTypes:
    Radius2;57;NumericalValues;IDsAndValuesString;24769214 508787397 85039
    3286 862102353 981185303 1231636850 1649386610 1941540674 263599683 32
    9205671 571109041 639579325 683993318 723853089 810600886 885767127...;
    1 1 1 10 2 22 3 1 3 3 1 1 1 3 2 2 1 2 2 2 3 1 1 1 1 1 14 1 1 1 1 1 1 2
    1 2 1 1 2 2 1 1 2 1 1 1 2 1 1 2 1 1 1 1 1 1 1

    FingerprintsBitVector;ExtendedConnectivityBits:FunctionalClassAtomType
    s:Radius2;1024;BinaryString;Ascending;00000000000000000000100000000000
    0000000001000100000000001000000000000000000000000000000000101000000010
    0000001000000000010000000000000000000000000000000000000000000000000100
    0000000000001000000000000001000000000001001000000000000000000000000000
    0000000000000000100000000000001000000000000000000000000000000000000...

    FingerprintsVector;ExtendedConnectivity:DREIDINGAtomTypes:Radius2;56;A
    lphaNumericalValues;ValuesString;280305427 357928343 721790579 1151822
    898 1207111054 1380963747 1568213839 1603445250 4559268 55012922 18094
    0813 335715751 534801009 684609658 829361048 972945982 999881534 10076
    55741 1213692591 1222032501 1224517934 1235687794 1244268533 152812070
    0 1629595024 1856308891 1978806036 2001865095 2096549435 172675415 ...

    FingerprintsVector;ExtendedConnectivity:EStateAtomTypes:Radius2;62;Alp
    haNumericalValues;ValuesString;25189973 528584866 662581668 671034184
    926543080 1347067490 1738510057 1759600920 2034425745 2097234755 21450
    44754 96779665 180364292 341712110 345278822 386540408 387387308 50430
    1706 617094135 771528807 957666640 997798220 1158349170 1291258082 134
    1138533 1395329837 1420277211 1479584608 1486476397 1487556246 1566...

    FingerprintsVector;ExtendedConnectivity:MMFF94AtomTypes:Radius2;64;Alp
    haNumericalValues;ValuesString;224051550 746527773 998750766 103704190
    2 1239701709 1248384926 1259447756 1521678386 1631549126 1909437580 20
    37095052 2104274756 2117729376 8770364 31445800 81450228 314289324 344
    041929 581773587 638555787 692022098 811840536 929651561 936421792 988
    636432 1048624296 1054288509 1369487579 1454058929 1519352190 17271...

    FingerprintsVector;ExtendedConnectivity:SLogPAtomTypes:Radius2;71;Alph
    aNumericalValues;ValuesString;78989290 116507218 489454042 888737940 1
    162561799 1241797255 1251494264 1263717127 1471206899 1538061784 17654
    07295 1795036542 1809833874 2020454493 2055310842 2117729376 11868981
    56731842 149505242 184525155 196984339 288181334 481409282 556716568 6
    41915747 679881756 721736571 794256218 908276640 992898760 10987549...

    FingerprintsVector;ExtendedConnectivity:SYBYLAtomTypes:Radius2;58;Alph
    aNumericalValues;ValuesString;199957044 313356892 455463968 465982819
    1225318176 1678585943 1883366064 1963811677 2117729376 113784599 19153
    8837 196629033 263865277 416380653 477036669 681527491 730724924 90906
    5537 1021959189 1133014972 1174311016 1359441203 1573452838 1661585138
    1668649038 1684198062 1812312554 1859266290 1891651106 2072549404 ...

    FingerprintsVector;ExtendedConnectivity:TPSAAtomTypes:Radius2;47;Alpha
    NumericalValues;ValuesString;20818206 259344053 862102353 1331904542 1
    700688206 265614156 363161397 681332588 810600886 885767127 950172500
    951454814 1059668746 1247054493 1382302230 1399502637 1805025917 19189
    39561 2114677228 2126402271 8130483 17645742 32278373 149975755 160327
    654 256360355 279492740 291251259 317592700 333763396 972105960 101...

    FingerprintsVector;ExtendedConnectivity:UFFAtomTypes:Radius2;56;AlphaN
    umericalValues;ValuesString;280305427 357928343 721790579 1151822898 1
    207111054 1380963747 1568213839 1603445250 4559268 55012922 180940813
    335715751 534801009 684609658 829361048 972945982 999881534 1007655741
    1213692591 1222032501 1224517934 1235687794 1244268533 1528120700 162
    9595024 1856308891 1978806036 2001865095 2096549435 172675415 18344...

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

=item B<-a, --AtomIdentifierType> I<AtomicInvariantsAtomTypes | FunctionalClassAtomTypes | DREIDINGAtomTypes | EStateAtomTypes | MMFF94AtomTypes | SLogPAtomTypes | SYBYLAtomTypes | TPSAAtomTypes | UFFAtomTypes>

Specify atom identifier type to use for assignment of initial atom identifier to non-hydrogen
atoms during calculation of extended connectivity fingerprints [ Ref 48, Ref 52]. Possible values
in the current release are: I<AtomicInvariantsAtomTypes, FunctionalClassAtomTypes,
DREIDINGAtomTypes, EStateAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes,
TPSAAtomTypes, UFFAtomTypes>. Default value: I<AtomicInvariantsAtomTypes>.

=item B<--AtomicInvariantsToUse> I<"AtomicInvariant,AtomicInvariant...">

This value is used during I<AtomicInvariantsAtomTypes> value of B<a, --AtomIdentifierType>
option. It's a list of comma separated valid atomic invariant atom types.

Possible values for atomic invarians are: I<AS, X, BO,  LBO, SB, DB, TB,
H, Ar, RA, FC, MN, SM>. Default value [ Ref 24 ]: I<AS,X,BO,H,FC,MN>.

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

=item B<--BitsOrder> I<Ascending | Descending>

Bits order to use during generation of fingerprints bit-vector string for I<ExtendedConnectivityBits>
value of B<-m, --mode> option. Possible values: I<Ascending, Descending>. Default: I<Ascending>.

I<Ascending> bit order which corresponds to first bit in each byte as the lowest bit as
opposed to the highest bit.

Internally, bits are stored in I<Ascending> order using Perl vec function. Regardless
of machine order, big-endian or little-endian, vec function always considers first
string byte as the lowest byte and first bit within each byte as the lowest bit.

=item B<-b, --BitStringFormat> I<BinaryString | HexadecimalString>

Format of fingerprints bit-vector string data in output SD, FP or CSV/TSV text file(s) specified by
B<--output> used during I<ExtendedConnectivityBits> value of B<-m, --mode> option. Possible
values: I<BinaryString, HexadecimalString>. Default value: I<BinaryString>.

I<BinaryString> corresponds to an ASCII string containing 1s and 0s. I<HexadecimalString>
contains bit values in ASCII hexadecimal format.

Examples:

    FingerprintsBitVector;ExtendedConnectivityBits:AtomicInvariantsAtomTyp
    es:Radius2;1024;BinaryString;Ascending;0000000000000000000000000000100
    0000000001010000000110000011000000000000100000000000000000000000100001
    1000000110000000000000000000000000010011000000000000000000000000010000
    0000000000000000000000000010000000000000000001000000000000000000000000
    0000000000010000100001000000000000101000000000000000100000000000000...

    FingerprintsBitVector;ExtendedConnectivityBits:FunctionalClassAtomType
    s:Radius2;1024;BinaryString;Ascending;00000000000000000000100000000000
    0000000001000100000000001000000000000000000000000000000000101000000010
    0000001000000000010000000000000000000000000000000000000000000000000100
    0000000000001000000000000001000000000001001000000000000000000000000000
    0000000000000000100000000000001000000000000000000000000000000000000...

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

=item B<--FingerprintsLabel> I<text>

SD data label or text file column label to use for fingerprints string in output SD or
CSV/TSV text file(s) specified by B<--output>. Default value: I<ExtendedConnectivityFingerprints>.

=item B<-h, --help>

Print this help message.

=item B<-k, --KeepLargestComponent> I<Yes | No>

Generate fingerprints for only the largest component in molecule. Possible values:
I<Yes or No>. Default value: I<Yes>.

For molecules containing multiple connected components, fingerprints can be generated
in two different ways: use all connected components or just the largest connected
component. By default, all atoms except for the largest connected component are
deleted before generation of fingerprints.

=item B<-m, --mode> I<ExtendedConnectivity | ExtendedConnectivityCount | ExtendedConnectivityBits>

Specify type of extended connectivity fingerprints to generate for molecules in I<SDFile(s)>.
Possible values: I<ExtendedConnectivity, ExtendedConnecticityCount or
ExtendedConnectivityBits>. Default value: I<ExtendedConnectivity>.

For I<ExtendedConnnectivity> value of fingerprints B<-m, --mode>, a fingerprint vector
containing unique atom identifiers constitute the extended connectivity fingerprints
of a molecule.

For I<ExtendedConnnectivityCount> value of fingerprints B<-m, --mode>, a fingerprint vector
containing unique atom identifiers along with their count constitute the extended connectivity
fingerprints of a molecule.

For I<ExtendedConnnectivityBits> value of fingerprints B<-m, --mode>, a fingerprint bit vector
indicating presence/absence of structurally unique atom identifiers constitute the extended
connectivity fingerprints of a molecule.

=item B<-n, --NeighborhoodRadius> I<number>

Atomic neighborhood radius for generating extended connectivity neighborhoods. Default
value: I<2>.  Valid values: >= 0. Neighborhood radius of zero correspond to just the list
of non-hydrogen atoms.

Default value of I<2> for atomic neighborhood radius generates extended connectivity
fingerprints corresponding to path length or diameter value of I<4> [ Ref 52b ].

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

New file name is generated using the root: <Root>.<Ext>. Default for new file names:
<SDFileName><ExtendedConnectivityFP>.<Ext>. The file type determines <Ext>
value. The sdf, fpf, csv, and tsv <Ext> values are used for SD, FP, comma/semicolon, and tab
delimited text files, respectively.This option is ignored for multiple input files.

=item B<-s, --size> I<number>

Size of bit-vector to use during generation of fingerprints bit-vector string for
I<ExtendedConnectivityBits> value of B<-m, --mode>. Default value: I<1024>.
Valid values correspond to any positive integer which satisfies the following criteria:
power of 2, >= 32 and <= 2 ** 32.

Examples:

   512
   1024
   2048

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

=item B<-v, --VectorStringFormat> I<ValuesString | IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>

Format of fingerprints vector string data in output SD, FP or CSV/TSV text file(s) specified by
B<--output> used during <ExtendedConnectivityCount> value of B<-m, --mode> option. Possible
values: I<ValuesString, IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString |
ValuesAndIDsPairsString>.

Default value during <ExtendedConnectivityCount> value of B<-m, --mode> option:
I<IDsAndValuesString>.

Default value during <ExtendedConnectivity> value of B<-m, --mode> option: I<ValuesString>.

Examples:

    FingerprintsVector;ExtendedConnectivity:AtomicInvariantsAtomTypes:Radi
    us2;60;AlphaNumericalValues;ValuesString;73555770 333564680 352413391
    666191900 1001270906 1371674323 1481469939 1977749791 2006158649 21414
    08799 49532520 64643108 79385615 96062769 273726379 564565671 85514103
    5 906706094 988546669 1018231313 1032696425 1197507444 1331250018 1338
    532734 1455473691 1607485225 1609687129 1631614296 1670251330 17303...

    FingerprintsVector;ExtendedConnectivityCount:AtomicInvariantsAtomTypes
    :Radius2;60;NumericalValues;IDsAndValuesString;73555770 333564680 3524
    13391 666191900 1001270906 1371674323 1481469939 1977749791 2006158649
    2141408799 49532520 64643108 79385615 96062769 273726379 564565671...;
    3 2 1 1 14 1 2 10 4 3 1 1 1 1 2 1 2 1 1 1 2 3 1 1 2 1 3 3 8 2 2 2 6 2
    1 2 1 1 2 1 1 1 2 1 1 2 1 2 1 1 1 1 1 1 1 1 1 2 1 1

=item B<-w, --WorkingDir> I<DirName>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using atomic invariants atom types in vector string format and create a SampleECAIFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -r SampleECAIFP -o Sample.sdf

To generate extended connectivity count fingerprints corresponding to neighborhood radius up to
2 using atomic invariants atom types in vector string format and create a SampleECAIFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -m ExtendedConnectivityCount
      -r SampleECAIFP -o Sample.sdf

To generate extended connectivity bits fingerprints as hexadecimal bit-string corresponding to
neighborhood radius up to 2 using atomic invariants atom types in vector string format and
create a SampleECAIFP.csv file containing sequential compound IDs along with fingerprints
vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -m ExtendedConnectivityBits
      -r SampleECAIFP -o Sample.sdf

To generate extended connectivity bits fingerprints as binary bit-string corresponding to
neighborhood radius up to 2 using atomic invariants atom types in vector string format and
create a SampleECAIFP.csv file containing sequential compound IDs along with fingerprints
vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -m ExtendedConnectivityBits
      --BitStringFormat BinaryString -r SampleECAIFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using atomic invariants atom types in vector string format and create SampleECAIFP.sdf, SampleECAIFP.fpf
and SampleECAIFP.csv files containing sequential compound IDs in CSV file along with fingerprints
vector strings data, type:

    % ExtendedConnectivityFingerprints.pl --output all -r SampleECAIFP
      -o Sample.sdf

To generate extended connectivity count fingerprints corresponding to neighborhood radius up to
2 using atomic invariants atom types in vector string format and create SampleECAIFP.sdf, SampleECAIFP.fpf
and SampleECAIFP.csv files containing sequential compound IDs in CSV file along with fingerprints
vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -m ExtendedConnectivityCount
      --output all -r SampleECAIFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using functional class atom types in vector string format and create a SampleECFCFP.csv file
containing sequential compound IDs along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -a FunctionalClassAtomTypes
      -r SampleECFCFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using DREIDING atom types in vector string format and create a SampleECFP.csv file
containing sequential compound IDs along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -a DREIDINGAtomTypes
      -r SampleECFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using E-state atom types in vector string format and create a SampleECFP.csv file
containing sequential compound IDs along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -a EStateAtomTypes
      -r SampleECFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using MMFF94 atom types in vector string format and create a SampleECFP.csv file
containing sequential compound IDs along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -a MMFF94AtomTypes
      -r SampleECFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using SLogP atom types in vector string format and create a SampleECFP.csv file
containing sequential compound IDs along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -a SLogPAtomTypes
      -r SampleECFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using SYBYL atom types in vector string format and create a SampleECFP.csv file
containing sequential compound IDs along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -a SYBYLAtomTypes
      -r SampleECFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using TPSA atom types in vector string format and create a SampleECFP.csv file
containing sequential compound IDs along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -a TPSAAtomTypes
      -r SampleECFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using UFF atom types in vector string format and create a SampleECFP.csv file
containing sequential compound IDs along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -a UFFAtomTypes
      -r SampleECFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
3 using atomic invariants atom types in vector string format and create a SampleECAIFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -a AtomicInvariantsAtomTypes -n 3
      -r SampleECAIFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
3 using functional class atom types in vector string format and create a SampleECFCFP.csv file
containing sequential compound IDs along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -a FunctionalClassAtomTypes -n 3
      -r SampleECFCFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using only AS,X atomic invariants atom types in vector string format and create a
SampleECAIFP.csv file containing sequential compound IDs along with fingerprints vector
strings data, type:

    % ExtendedConnectivityFingerprints.pl -a AtomicInvariantsAtomTypes
      --AtomicInvariantsToUse "AS,X" -r SampleECAIFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using only HBD,HBA functional class atom types in vector string format and create a
SampleECFCFP.csv file containing sequential compound IDs along with fingerprints vector
strings data, type:

    % ExtendedConnectivityFingerprints.pl -a FunctionalClassAtomTypes
      --FunctionalClassesToUse "HBD,HBA" -r SampleECFCFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using atomic invariants atom types in vector string format and create a SampleECAIFP.csv
file containing compound ID from molecule name line along with fingerprints vector strings
data, type:

    % ExtendedConnectivityFingerprints.pl -a AtomicInvariantsAtomTypes
      --DataFieldsMode CompoundID -CompoundIDMode MolName
      -r SampleECAIFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using functional class atom types in vector string format and create a SampleECFCFP.csv
file containing compound IDs using specified data field along with fingerprints vector strings
data, type:

    % ExtendedConnectivityFingerprints.pl -a FunctionalClassAtomTypes
      --DataFieldsMode CompoundID -CompoundIDMode DataField --CompoundID Mol_ID
      -r SampleECFCFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using atomic invariants atom types in vector string format and create a SampleECAIFP.tsv
file containing compound ID using combination of molecule name line and an explicit compound
prefix along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -a AtomicInvariantsAtomTypes
      --DataFieldsMode CompoundID -CompoundIDMode MolnameOrLabelPrefix
      --CompoundID Cmpd --CompoundIDLabel MolID -r SampleECAIFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using functional class atom types in vector string format and create a SampleECFCFP.csv
file containing specific data fields columns along with fingerprints vector strings
data, type:

    % ExtendedConnectivityFingerprints.pl -a FunctionalClassAtomTypes
      --DataFieldsMode  Specify --DataFields Mol_ID -r SampleECFCFP
      -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using atomic invariants atom types in vector string format and create a SampleECAIFP.tsv
file containing common data fields columns along with fingerprints vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -a AtomicInvariantsAtomTypes
      --DataFieldsMode Common -r SampleECAIFP -o Sample.sdf

To generate extended connectivity fingerprints corresponding to neighborhood radius up to
2 using functional class atom types in vector string format and create SampleECFCFP.sdf, SampleECFCFP.fpf
and SampleECFCFP.csv files containing all data fields columns in CSV file along with fingerprints
vector strings data, type:

    % ExtendedConnectivityFingerprints.pl -a FunctionalClassAtomTypes
      --DataFieldsMode All  --output all -r SampleECFCFP
      -o Sample.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoFingerprintsFiles.pl, SimilarityMatricesFingerprints.pl, AtomNeighborhoodsFingerprints.pl,
MACCSKeysFingerprints.pl, PathLengthFingerprints.pl,
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
