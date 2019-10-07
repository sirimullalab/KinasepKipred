#!/usr/bin/perl -w
#
# File: TopologicalAtomTripletsFingerprints.pl
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
use Fingerprints::TopologicalAtomTripletsFingerprints;

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
    GenerateTopologicalAtomTripletsFingerprints($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Generate fingerprints for a SD file...
#
sub GenerateTopologicalAtomTripletsFingerprints {
  my($FileIndex) = @_;
  my($CmpdCount, $IgnoredCmpdCount, $SDFile, $MoleculeFileIO, $Molecule, $TopologicalAtomTripletsFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO);

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

    $TopologicalAtomTripletsFingerprints = GenerateMoleculeFingerprints($Molecule);
    if (!$TopologicalAtomTripletsFingerprints) {
      $IgnoredCmpdCount++;
      ProcessIgnoredCompound('FingerprintsGenerationFailed', $CmpdCount, $Molecule);
      next COMPOUND;
    }

    WriteDataToOutputFiles($FileIndex, $CmpdCount, $Molecule, $TopologicalAtomTripletsFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO);
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
  %FingerprintsFileIOParams = ('Mode' => 'Write', 'Overwrite' => $OptionsInfo{OverwriteFiles}, 'FingerprintsStringMode' => 'FingerprintsVectorString', 'VectorStringFormat' => $OptionsInfo{VectorStringFormat});

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
  my($FileIndex, $CmpdCount, $Molecule, $TopologicalAtomTripletsFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO) = @_;
  my($DataFieldLabelAndValuesRef);

  $DataFieldLabelAndValuesRef = undef;
  if ($NewFPTextFileIO || $NewFPFileIO) {
    $DataFieldLabelAndValuesRef = $Molecule->GetDataFieldLabelAndValues();
  }

  if ($NewFPSDFileIO) {
    my($CmpdString);

    $CmpdString = $Molecule->GetInputMoleculeString();
    $NewFPSDFileIO->WriteFingerprints($TopologicalAtomTripletsFingerprints, $CmpdString);
  }

  if ($NewFPTextFileIO) {
    my($ColValuesRef);

    $ColValuesRef = SetupFPTextFileCoulmnValues($FileIndex, $CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPTextFileIO->WriteFingerprints($TopologicalAtomTripletsFingerprints, $ColValuesRef);
  }

  if ($NewFPFileIO) {
    my($CompoundID);

    $CompoundID = SetupCmpdIDForOutputFiles($CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPFileIO->WriteFingerprints($TopologicalAtomTripletsFingerprints, $CompoundID);
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
  my($TopologicalAtomTripletsFingerprints);

  if ($OptionsInfo{KeepLargestComponent}) {
    $Molecule->KeepLargestComponent();
  }
  if (!$Molecule->DetectRings()) {
    return undef;
  }
  $Molecule->SetAromaticityModel($OptionsInfo{AromaticityModel});
  $Molecule->DetectAromaticity();

  $TopologicalAtomTripletsFingerprints = new Fingerprints::TopologicalAtomTripletsFingerprints('Molecule' => $Molecule, 'MinDistance' => $OptionsInfo{MinDistance},  'MaxDistance' => $OptionsInfo{MaxDistance}, 'UseTriangleInequality' => $OptionsInfo{UseTriangleInequality}, 'AtomIdentifierType' => $OptionsInfo{AtomIdentifierType});
  SetAtomIdentifierTypeValuesToUse($TopologicalAtomTripletsFingerprints);

  # Generate fingerprints...
  $TopologicalAtomTripletsFingerprints->GenerateFingerprints();

  # Make sure fingerprints generation is successful...
  if (!$TopologicalAtomTripletsFingerprints->IsFingerprintsGenerationSuccessful()) {
    return undef;
  }

  return $TopologicalAtomTripletsFingerprints;
}

# Set atom identifier type to use for generating fingerprints...
#
sub SetAtomIdentifierTypeValuesToUse {
  my($TopologicalAtomTripletsFingerprints) = @_;

  if ($OptionsInfo{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
    $TopologicalAtomTripletsFingerprints->SetAtomicInvariantsToUse(\@{$OptionsInfo{AtomicInvariantsToUse}});
  }
  elsif ($OptionsInfo{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
    $TopologicalAtomTripletsFingerprints->SetFunctionalClassesToUse(\@{$OptionsInfo{FunctionalClassesToUse}});
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
      $OutFileRoot = "${FileName}TopologicalAtomTripletsFP";
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

  $OptionsInfo{Filter} = ($Options{filter} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{FingerprintsLabel} = $Options{fingerprintslabel} ? $Options{fingerprintslabel} : 'TopologicalAtomTripletsFingerprints';

  $OptionsInfo{KeepLargestComponent} = ($Options{keeplargestcomponent} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{MinDistance} = $Options{mindistance};
  $OptionsInfo{MaxDistance} = $Options{maxdistance};

  $OptionsInfo{Output} = $Options{output};
  $OptionsInfo{SDOutput} = ($Options{output} =~ /^(SD|All)$/i) ? 1 : 0;
  $OptionsInfo{FPOutput} = ($Options{output} =~ /^(FP|All)$/i) ? 1 : 0;
  $OptionsInfo{TextOutput} = ($Options{output} =~ /^(Text|All)$/i) ? 1 : 0;

  $OptionsInfo{OutDelim} = $Options{outdelim};
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? 1 : 0;
  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : 0;

  $OptionsInfo{UseTriangleInequality} = ($Options{usetriangleinequality} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{VectorStringFormat} = $Options{vectorstringformat};
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
  $Options{atomicinvariantstouse} = 'AS,X,BO,H,FC';

  $Options{functionalclassestouse} = 'HBD,HBA,PI,NI,Ar,Hal';

  $Options{compoundidmode} = 'LabelPrefix';
  $Options{compoundidlabel} = 'CompoundID';
  $Options{datafieldsmode} = 'CompoundID';

  $Options{filter} = 'Yes';

  $Options{keeplargestcomponent} = 'Yes';

  $Options{mindistance} = 1;
  $Options{maxdistance} = 10;

  $Options{output} = 'text';
  $Options{outdelim} = 'comma';
  $Options{quote} = 'yes';

  $Options{usetriangleinequality} = 'No';

  $Options{vectorstringformat} = 'IDsAndValuesString';

  if (!GetOptions(\%Options, "aromaticitymodel=s", "atomidentifiertype|a=s", "atomicinvariantstouse=s", "functionalclassestouse=s", "compoundid=s", "compoundidlabel=s", "compoundidmode=s", "datafields=s", "datafieldsmode|d=s", "filter|f=s", "fingerprintslabel=s",  "help|h", "keeplargestcomponent|k=s",  "mindistance=s", "maxdistance=s", "outdelim=s", "output=s", "overwrite|o", "quote|q=s", "root|r=s", , "usetriangleinequality|u=s", "vectorstringformat|v=s", "workingdir|w=s")) {
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
  if (!IsPositiveInteger($Options{mindistance})) {
    die "Error: The value specified, $Options{mindistance}, for option \"--MinDistance\" is not valid. Allowed values: > 0 \n";
  }
  if (!IsPositiveInteger($Options{maxdistance})) {
    die "Error: The value specified, $Options{maxdistance}, for option \"--MaxDistance\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{mindistance} > $Options{maxdistance}) {
    die "Error: The value specified, specified, $Options{mindistance}, for option \"--MinDistance\" must be less than the value specified, $Options{maxdistance}, for option \"--MaxDistance\" \n";
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
  if ($Options{usetriangleinequality} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{usetriangleinequality}, for option \"-u, --UseTriangleInequality\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{vectorstringformat} !~ /^(IDsAndValuesString|IDsAndValuesPairsString|ValuesAndIDsString|ValuesAndIDsPairsString)$/i) {
    die "Error: The value specified, $Options{vectorstringformat}, for option \"-v, --VectorStringFormat\" is not valid. Allowed values: IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString or ValuesAndIDsPairsString\n";
  }
}

__END__

=head1 NAME

TopologicalAtomTripletsFingerprints.pl - Generate topological atom triplets fingerprints for SD files

=head1 SYNOPSIS

TopologicalAtomTripletsFingerprints.pl SDFile(s)...

TopologicalAtomTripletsFingerprints.pl [B<--AromaticityModel> I<AromaticityModelType>]
[B<-a, --AtomIdentifierType> I<AtomicInvariantsAtomTypes>]
[B<--AtomicInvariantsToUse> I<"AtomicInvariant,AtomicInvariant...">]
[B<--FunctionalClassesToUse> I<"FunctionalClass1,FunctionalClass2...">]
[B<--CompoundID> I<DataFieldName or LabelPrefixString>] [B<--CompoundIDLabel> I<text>]
[B<--CompoundIDMode>] [B<--DataFields> I<"FieldLabel1,FieldLabel2,...">]
[B<-d, --DataFieldsMode> I<All | Common | Specify | CompoundID>] [B<-f, --Filter> I<Yes | No>]
[B<--FingerprintsLabel> I<text>] [B<-h, --help>] [B<-k, --KeepLargestComponent> I<Yes | No>]
[B<--MinDistance> I<number>] [B<--MaxDistance> I<number>]
[B<--OutDelim> I<comma | tab | semicolon>] [B<--output> I<SD | FP | text | all>] [B<-o, --overwrite>]
[B<-q, --quote> I<Yes | No>] [B<-r, --root> I<RootName>]  [B<-u, --UseTriangleInequality> I<Yes | No>]
[B<-v, --VectorStringFormat> I<ValuesString, IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>]
[B<-w, --WorkingDir> dirname] SDFile(s)...

=head1 DESCRIPTION

Generate topological atom triplets fingerprints  for I<SDFile(s)> and create
appropriate SD, FP or CSV/TSV text file(s) containing fingerprints vector strings corresponding to
molecular fingerprints.

Multiple SDFile names are separated by spaces. The valid file extensions are I<.sdf>
and I<.sd>. All other file names are ignored. All the SD files in a current directory
can be specified either by I<*.sdf> or the current directory name.

The current release of MayaChemTools supports generation of topological atom triplets
fingerprints corresponding to following B<-a, --AtomIdentifierTypes>:

    AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
    FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
    SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes

Based on the values specified for B<-a, --AtomIdentifierType> and B<--AtomicInvariantsToUse>,
initial atom types are assigned to all non-hydrogen atoms in a molecule. Using the distance
matrix for the molecule and initial atom types assigned to non-hydrogen atoms, all unique atom
pairs within B<--MinDistance> and B<--MaxDistance> are identified and counted. An atom triplet
identifier is generated for each unique atom triplet; the format of the atom triplet identifier is:

    <ATx>-Dyz-<ATy>-Dxz-<ATz>-Dxy

    ATx, ATy, ATz: Atom types assigned to atom x, atom y, and atom z
    Dxy: Distance between atom x and atom y
    Dxz: Distance between atom x and atom z
    Dyz: Distance between atom y and atom z

    where <AT1>-D23 <= <AT2>-D13 <= <AT3>-D12

The atom triplet identifiers for all unique atom triplets corresponding to non-hydrogen atoms constitute
topological atom triplets fingerprints of the molecule.

Example of I<SD> file containing topological atom triplets fingerprints string data:

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

    >  <TopologicalAtomTripletsFingerprints>
    FingerprintsVector;TopologicalAtomTriplets:AtomicInvariantsAtomTypes:Mi
    nDistance1:MaxDistance10;3096;NumericalValues;IDsAndValuesString;C.X1.B
    O1.H3-D1-C.X1.BO1.H3-D1-C.X3.BO3.H1-D2 C.X1.BO1.H3-D1-C.X2.BO2.H2-D10-C
    .X3.BO4-D9 C.X1.BO1.H3-D1-C.X2.BO2.H2-D3-N.X3.BO3-D4 C.X1.BO1.H3-D1...;
    1 2 2 2 2 2 2 2 8 8 4 8 4 4 2 2 2 2 4 2 2 2 4 2 2 2 2 1 2 2 4 4 4 2 2 2
    4 4 4 8 4 4 2 4 4 4 2 4 4 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 8 8 ...

    $$$$
    ... ...
    ... ...

Example of I<FP> file containing topological atom triplets fingerprints string data:

    #
    # Package = MayaChemTools 7.4
    # Release Date = Oct 21, 2010
    #
    # TimeStamp = Fri Mar 11 15:24:01 2011
    #
    # FingerprintsStringType = FingerprintsVector
    #
    # Description = TopologicalAtomTriplets:AtomicInvariantsAtomTypes:Mi...
    # VectorStringFormat = IDsAndValuesString
    # VectorValuesType = NumericalValues
    #
    Cmpd1 3096;C.X1.BO1.H3-D1-C.X1.BO1.H3-D1-C.X3.BO3.H1-D2...;1 2 2 2 2...
    Cmpd2 1093;C.X1.BO1.H3-D1-C.X1.BO1.H3-D3-C.X2.BO2.H2-D4...;2 2 2 2 2... 
    ... ...
    ... ..

Example of CSV I<Text> file containing topological atom triplets fingerprints string data:

    "CompoundID","TopologicalAtomTripletsFingerprints"
    "Cmpd1","FingerprintsVector;TopologicalAtomTriplets:AtomicInvariantsAto
    mTypes:MinDistance1:MaxDistance10;3096;NumericalValues;IDsAndValuesStri
    ng;C.X1.BO1.H3-D1-C.X1.BO1.H3-D1-C.X3.BO3.H1-D2 C.X1.BO1.H3-D1-C.X2.BO2
    .H2-D10-C.X3.BO4-D9 C.X1.BO1.H3-D1-C.X2.BO2.H2-D3-N.X3.BO3-D4 C.X1....;
    1 2 2 2 2 2 2 2 8 8 4 8 4 4 2 2 2 2 4 2 2 2 4 2 2 2 2 1 2 2 4 4 4 2 2 2
    4 4 4 8 4 4 2 4 4 4 2 4 4 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 8 8 ...
    ... ...
    ... ...

The current release of MayaChemTools generates the following types of topological atom triplets
fingerprints vector strings:

    FingerprintsVector;TopologicalAtomTriplets:AtomicInvariantsAtomTypes:M
    inDistance1:MaxDistance10;3096;NumericalValues;IDsAndValuesString;C.X1
    .BO1.H3-D1-C.X1.BO1.H3-D1-C.X3.BO3.H1-D2 C.X1.BO1.H3-D1-C.X2.BO2.H2-D1
    0-C.X3.BO4-D9 C.X1.BO1.H3-D1-C.X2.BO2.H2-D3-N.X3.BO3-D4 C.X1.BO1.H3-D1
    -C.X2.BO2.H2-D4-C.X2.BO2.H2-D5 C.X1.BO1.H3-D1-C.X2.BO2.H2-D6-C.X3....;
    1 2 2 2 2 2 2 2 8 8 4 8 4 4 2 2 2 2 4 2 2 2 4 2 2 2 2 1 2 2 4 4 4 2 2
    2 4 4 4 8 4 4 2 4 4 4 2 4 4 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 8...

    FingerprintsVector;TopologicalAtomTriplets:AtomicInvariantsAtomTypes:M
    inDistance1:MaxDistance10;3096;NumericalValues;IDsAndValuesPairsString
    ;C.X1.BO1.H3-D1-C.X1.BO1.H3-D1-C.X3.BO3.H1-D2 1 C.X1.BO1.H3-D1-C.X2.BO
    2.H2-D10-C.X3.BO4-D9 2 C.X1.BO1.H3-D1-C.X2.BO2.H2-D3-N.X3.BO3-D4 2 C.X
    1.BO1.H3-D1-C.X2.BO2.H2-D4-C.X2.BO2.H2-D5 2 C.X1.BO1.H3-D1-C.X2.BO2.H2
    -D6-C.X3.BO3.H1-D5 2 C.X1.BO1.H3-D1-C.X2.BO2.H2-D6-C.X3.BO3.H1-D7 2...

    FingerprintsVector;TopologicalAtomTriplets:DREIDINGAtomTypes:MinDistan
    ce1:MaxDistance10;2377;NumericalValues;IDsAndValuesString;C_2-D1-C_2-D
    9-C_3-D10 C_2-D1-C_2-D9-C_R-D10 C_2-D1-C_3-D1-C_3-D2 C_2-D1-C_3-D10-C_
    3-D9 C_2-D1-C_3-D2-C_3-D3 C_2-D1-C_3-D2-C_R-D3 C_2-D1-C_3-D3-C_3-D4 C_
    2-D1-C_3-D3-N_R-D4 C_2-D1-C_3-D3-O_3-D2 C_2-D1-C_3-D4-C_3-D5 C_2-D...;
    1 1 1 2 1 1 3 1 1 2 2 1 1 1 1 1 1 1 1 2 1 3 4 5 1 1 6 4 2 2 3 1 1 1 2
    2 1 2 1 1 2 2 2 1 2 1 2 1 1 3 3 2 6 4 2 1 1 1 2 2 1 1 1 1 1 1 1 1 1...

    FingerprintsVector;TopologicalAtomTriplets:EStateAtomTypes:MinDistance
    1:MaxDistance10;3298;NumericalValues;IDsAndValuesString;aaCH-D1-aaCH-D
    1-aaCH-D2 aaCH-D1-aaCH-D1-aasC-D2 aaCH-D1-aaCH-D10-aaCH-D9 aaCH-D1-aaC
    H-D10-aasC-D9 aaCH-D1-aaCH-D2-aaCH-D3 aaCH-D1-aaCH-D2-aasC-D1 aaCH-D1-
    aaCH-D2-aasC-D3 aaCH-D1-aaCH-D3-aasC-D2 aaCH-D1-aaCH-D4-aasC-D5 aa...;
    6 4 24 4 16 8 8 4 8 8 8 12 10 14 4 16 24 4 12 2 2 4 1 10 2 2 15 2 2 2
    2 2 2 14 4 2 2 2 2 1 2 10 2 2 4 1 2 4 8 3 3 3 4 6 4 2 2 3 3 1 1 1 2 1
    2 2 4 2 3 2 1 2 4 5 3 2 2 1 2 4 3 2 8 12 6 2 2 4 4 7 1 4 2 4 2 2 2 ...

    FingerprintsVector;TopologicalAtomTriplets:FunctionalClassAtomTypes:Mi
    nDistance1:MaxDistance10;2182;NumericalValues;IDsAndValuesString;Ar-D1
    -Ar-D1-Ar-D2 Ar-D1-Ar-D1-Ar.HBA-D2 Ar-D1-Ar-D10-Ar-D9 Ar-D1-Ar-D10-Hal
    -D9 Ar-D1-Ar-D2-Ar-D2 Ar-D1-Ar-D2-Ar-D3 Ar-D1-Ar-D2-Ar.HBA-D1 Ar-D1-Ar
    -D2-Ar.HBA-D2 Ar-D1-Ar-D2-Ar.HBA-D3 Ar-D1-Ar-D2-HBD-D1 Ar-D1-Ar-D2...;
    27 1 32 2 2 63 3 2 1 2 1 2 3 1 1 40 3 1 2 2 2 2 4 2 2 47 4 2 2 1 2 1 5
    2 2 51 4 3 1 3 1 9 1 1 50 3 3 4 1 9 50 2 2 3 3 5 45 1 1 1 2 1 2 2 3 3
    4 4 3 2 1 1 3 4 5 5 3 1 2 3 2 3 5 7 2 7 3 7 1 1 2 2 2 2 3 1 4 3 1 2...

    FingerprintsVector;TopologicalAtomTriplets:MMFF94AtomTypes:MinDistance
    1:MaxDistance10;2966;NumericalValues;IDsAndValuesString;C5A-D1-C5A-D1-
    N5-D2 C5A-D1-C5A-D2-C5B-D2 C5A-D1-C5A-D3-CB-D2 C5A-D1-C5A-D3-CR-D2 C5A
    -D1-C5B-D1-C5B-D2 C5A-D1-C5B-D2-C=ON-D1 C5A-D1-C5B-D2-CB-D1 C5A-D1-C5B
    -D3-C=ON-D2 C5A-D1-C5B-D3-CB-D2 C5A-D1-C=ON-D3-NC=O-D2 C5A-D1-C=ON-D3-
    O=CN-D2 C5A-D1-C=ON-D4-NC=O-D3 C5A-D1-C=ON-D4-O=CN-D3 C5A-D1-CB-D1-...

    FingerprintsVector;TopologicalAtomTriplets:SLogPAtomTypes:MinDistance1
    :MaxDistance10;3710;NumericalValues;IDsAndValuesString;C1-D1-C1-D1-C11
    -D2 C1-D1-C1-D1-CS-D2 C1-D1-C1-D10-C5-D9 C1-D1-C1-D3-C10-D2 C1-D1-C1-D
    3-C5-D2 C1-D1-C1-D3-CS-D2 C1-D1-C1-D3-CS-D4 C1-D1-C1-D4-C10-D5 C1-D1-C
    1-D4-C11-D5 C1-D1-C1-D5-C10-D4 C1-D1-C1-D5-C5-D4 C1-D1-C1-D6-C11-D7 C1
    -D1-C1-D6-CS-D5 C1-D1-C1-D6-CS-D7 C1-D1-C1-D8-C11-D9 C1-D1-C1-D8-CS...

    FingerprintsVector;TopologicalAtomTriplets:SYBYLAtomTypes:MinDistance1
    :MaxDistance10;2332;NumericalValues;IDsAndValuesString;C.2-D1-C.2-D9-C
    .3-D10 C.2-D1-C.2-D9-C.ar-D10 C.2-D1-C.3-D1-C.3-D2 C.2-D1-C.3-D10-C.3-
    D9 C.2-D1-C.3-D2-C.3-D3 C.2-D1-C.3-D2-C.ar-D3 C.2-D1-C.3-D3-C.3-D4 C.2
    -D1-C.3-D3-N.ar-D4 C.2-D1-C.3-D3-O.3-D2 C.2-D1-C.3-D4-C.3-D5 C.2-D1-C.
    3-D5-C.3-D6 C.2-D1-C.3-D5-O.3-D4 C.2-D1-C.3-D6-C.3-D7 C.2-D1-C.3-D7...

    FingerprintsVector;TopologicalAtomTriplets:TPSAAtomTypes:MinDistance1:
    MaxDistance10;1007;NumericalValues;IDsAndValuesString;N21-D1-N7-D3-Non
    e-D4 N21-D1-N7-D5-None-D4 N21-D1-None-D1-None-D2 N21-D1-None-D2-None-D
    2 N21-D1-None-D2-None-D3 N21-D1-None-D3-None-D4 N21-D1-None-D4-None-D5
     N21-D1-None-D4-O3-D3 N21-D1-None-D4-O4-D3 N21-D1-None-D5-None-D6 N21-
    D1-None-D6-None-D7 N21-D1-None-D6-O4-D5 N21-D1-None-D7-None-D8 N21-...

    FingerprintsVector;TopologicalAtomTriplets:UFFAtomTypes:MinDistance1:M
    axDistance10;2377;NumericalValues;IDsAndValuesString;C_2-D1-C_2-D9-C_3
    -D10 C_2-D1-C_2-D9-C_R-D10 C_2-D1-C_3-D1-C_3-D2 C_2-D1-C_3-D10-C_3-D9 
    C_2-D1-C_3-D2-C_3-D3 C_2-D1-C_3-D2-C_R-D3 C_2-D1-C_3-D3-C_3-D4 C_2-D1-
    C_3-D3-N_R-D4 C_2-D1-C_3-D3-O_3-D2 C_2-D1-C_3-D4-C_3-D5 C_2-D1-C_3-D5-
    C_3-D6 C_2-D1-C_3-D5-O_3-D4 C_2-D1-C_3-D6-C_3-D7 C_2-D1-C_3-D7-C_3-...

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

Specify atom identifier type to use for assignment of initial atom identifier to non-hydrogen
atoms during calculation of topological atom triplets fingerprints. Possible values in the current
release are: I<AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
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

Specify compound ID column label for CSV/TSV text file(s) used during I<CompoundID> value
of B<--DataFieldsMode> option. Default value: I<CompoundID>.

=item B<--CompoundIDMode> I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>

Specify how to generate compound IDs and write to FP or CSV/TSV text file(s) along with generated
fingerprints for I<FP | text | all> values of B<--output> option: use a I<SDFile(s)> datafield value;
use molname line from I<SDFile(s)>; generate a sequential ID with specific prefix; use combination
of both MolName and LabelPrefix with usage of LabelPrefix values for empty molname lines.

Possible values: I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>.
Default value: I<LabelPrefix>.

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
CSV/TSV text file(s) specified by B<--output>. Default value: I<TopologicalAtomTripletsFingerprints>.

=item B<-h, --help>

Print this help message.

=item B<-k, --KeepLargestComponent> I<Yes | No>

Generate fingerprints for only the largest component in molecule. Possible values:
I<Yes or No>. Default value: I<Yes>.

For molecules containing multiple connected components, fingerprints can be generated
in two different ways: use all connected components or just the largest connected
component. By default, all atoms except for the largest connected component are
deleted before generation of fingerprints.

=item B<--MinDistance> I<number>

Minimum bond distance between atom triplets for generating topological atom triplets. Default value:
I<1>. Valid values: positive integers and less than B<--MaxDistance>.

=item B<--MaxDistance> I<number>

Maximum bond distance between atom triplets for generating topological atom triplets. Default value:
I<10>. Valid values: positive integers and greater than B<--MinDistance>.

=item B<--OutDelim> I<comma | tab | semicolon>

Delimiter for output CSV/TSV text file(s). Possible values: I<comma, tab, or semicolon>
Default value: I<comma>

=item B<--output> I<SD | FP | text | all>

Type of output files to generate. Possible values: I<SD, FP, text, or all>. Default value: I<text>.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-q, --quote> I<Yes | No>

Put quote around column values in output CSV/TSV text file(s). Possible values:
I<Yes or No>. Default value: I<Yes>.

=item B<-r, --root> I<RootName>

New file name is generated using the root: <Root>.<Ext>. Default for new file names:
<SDFileName><TopologicalAtomTripletsFP>.<Ext>. The file type determines <Ext> value.
The sdf, fpf, csv, and tsv <Ext> values are used for SD, FP, comma/semicolon, and tab
delimited text files, respectively.This option is ignored for multiple input files.

=item B<-u, --UseTriangleInequality> I<Yes | No>

Specify whether to imply triangle distance inequality test to distances between atom pairs in
atom triplets during generation of atom triplets generation. Possible values: I<Yes or No>.
Default value: I<No>.

Triangle distance inequality test implies that distance or binned distance between any two atom
pairs in an atom triplet must be less than the sum of distances or binned distances between other
two atoms pairs and greater than the difference of their distances.

    For atom triplet ATx-Dyz-ATy-Dxz-ATz-Dxy to satisfy triangle inequality:

    Dyz > |Dxz - Dxy| and Dyz < Dxz + Dxy
    Dxz > |Dyz - Dxy| and Dyz < Dyz + Dxy
    Dxy > |Dyz - Dxz| and Dxy < Dyz + Dxz

=item B<-v, --VectorStringFormat> I<IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>

Format of fingerprints vector string data in output SD, FP or CSV/TSV text file(s) specified by
B<--output> option. Possible values: I<IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString |
ValuesAndIDsPairsString>. Default value: I<IDsAndValuesString>.

Examples:


    FingerprintsVector;TopologicalAtomTriplets:AtomicInvariantsAtomTypes:M
    inDistance1:MaxDistance10;3096;NumericalValues;IDsAndValuesString;C.X1
    .BO1.H3-D1-C.X1.BO1.H3-D1-C.X3.BO3.H1-D2 C.X1.BO1.H3-D1-C.X2.BO2.H2-D1
    0-C.X3.BO4-D9 C.X1.BO1.H3-D1-C.X2.BO2.H2-D3-N.X3.BO3-D4 C.X1.BO1.H3-D1
    -C.X2.BO2.H2-D4-C.X2.BO2.H2-D5 C.X1.BO1.H3-D1-C.X2.BO2.H2-D6-C.X3....;
    1 2 2 2 2 2 2 2 8 8 4 8 4 4 2 2 2 2 4 2 2 2 4 2 2 2 2 1 2 2 4 4 4 2 2
    2 4 4 4 8 4 4 2 4 4 4 2 4 4 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 8...

    FingerprintsVector;TopologicalAtomTriplets:AtomicInvariantsAtomTypes:M
    inDistance1:MaxDistance10;3096;NumericalValues;IDsAndValuesPairsString
    ;C.X1.BO1.H3-D1-C.X1.BO1.H3-D1-C.X3.BO3.H1-D2 1 C.X1.BO1.H3-D1-C.X2.BO
    2.H2-D10-C.X3.BO4-D9 2 C.X1.BO1.H3-D1-C.X2.BO2.H2-D3-N.X3.BO3-D4 2 C.X
    1.BO1.H3-D1-C.X2.BO2.H2-D4-C.X2.BO2.H2-D5 2 C.X1.BO1.H3-D1-C.X2.BO2.H2
    -D6-C.X3.BO3.H1-D5 2 C.X1.BO1.H3-D1-C.X2.BO2.H2-D6-C.X3.BO3.H1-D7 2...

=item B<-w, --WorkingDir> I<DirName>

Location of working directory. Default value: current directory.

=back

=head1 EXAMPLES

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using atomic invariants atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using atomic invariants atom types in IDsAndValuesString format and create SampleTATFP.sdf,
SampleTATFP.fpf and SampleTATFP.csv files containing sequential compound IDs in CSV file along
with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl --output all -r SampleTATFP
      -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using atomic invariants atom types in IDsAndValuesPairsString format and create a SampleTATFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl --VectorStringFormat
      IDsAndValuesPairsString  -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using DREIDING atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl -a DREIDINGAtomTypes
      -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using E-state atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl -a EStateAtomTypes
      -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using functional class atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl -a FunctionalClassAtomTypes
      -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using DREIDING atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl -a DREIDINGAtomTypes
      -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using MM94 atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl -a MMFF94AtomTypes
      -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using SLogP atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl -a SLogPAtomTypes
      -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using SYBYL atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl -a SYBYLAtomTypes
      -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using TPSA atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl -a TPSAAtomTypes
      -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using UFF atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl -a UFFAtomTypes
      -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
6 using atomic invariants atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing sequential compound IDs along with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl -a AtomicInvariantsAtomTypes
     --MinDistance 1 --MaxDistance 6 -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using only AS,X atomic invariants atom types in IDsAndValuesString format and create a
SampleTATFP.csv file containing sequential compound IDs along with fingerprints vector strings
data, type:

    % TopologicalAtomTripletsFingerprints.pl -a AtomicInvariantsAtomTypes
      --AtomicInvariantsToUse "AS,X" --MinDistance 1 --MaxDistance 6
      -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using atomic invariants atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing compound ID from molecule name line along with fingerprints vector strings
data, type:

    % TopologicalAtomTripletsFingerprints.pl -a AtomicInvariantsAtomTypes
      --DataFieldsMode CompoundID -CompoundIDMode MolName
      -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using atomic invariants atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing compound IDs using specified data field along with fingerprints vector strings
data, type:

    % TopologicalAtomTripletsFingerprints.pl -a AtomicInvariantsAtomTypes
      --DataFieldsMode CompoundID -CompoundIDMode DataField --CompoundID
      Mol_ID -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using atomic invariants atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing compound ID using combination of molecule name line and an explicit compound
prefix along with fingerprints vector strings data, type:

    % TopologicalAtomTripletsFingerprints.pl -a AtomicInvariantsAtomTypes
      --DataFieldsMode CompoundID -CompoundIDMode MolnameOrLabelPrefix
      --CompoundID Cmpd --CompoundIDLabel MolID -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using atomic invariants atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing specific data fields columns along with fingerprints vector strings
data, type:

    % TopologicalAtomTripletsFingerprints.pl -a AtomicInvariantsAtomTypes
      --DataFieldsMode Specify --DataFields Mol_ID -r SampleTATFP
      -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using atomic invariants atom types in IDsAndValuesString format and create a SampleTATFP.csv
file containing common data fields columns along with fingerprints vector strings
data, type:

    % TopologicalAtomTripletsFingerprints.pl -a AtomicInvariantsAtomTypes
      --DataFieldsMode Common -r SampleTATFP -o Sample.sdf

To generate topological atom triplets fingerprints corresponding to bond distances from 1 through
10 using atomic invariants atom types in IDsAndValuesString format and create SampleTATFP.sdf,
SampleTATFP.fpf and SampleTATFP.csv files containing all data fields columns in CSV file along with
fingerprints data, type:

    % TopologicalAtomTripletsFingerprints.pl -a AtomicInvariantsAtomTypes
      --DataFieldsMode All  --output all -r SampleTATFP
      -o Sample.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoFingerprintsFiles.pl, SimilarityMatricesFingerprints.pl, AtomNeighborhoodsFingerprints.pl,
ExtendedConnectivityFingerprints.pl, MACCSKeysFingerprints.pl,
PathLengthFingerprints.pl, TopologicalAtomTorsionsFingerprints.pl,
TopologicalPharmacophoreAtomPairsFingerprints.pl, TopologicalPharmacophoreAtomTripletsFingerprints.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
