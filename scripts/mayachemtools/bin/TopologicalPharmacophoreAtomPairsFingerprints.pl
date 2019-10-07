#!/usr/bin/perl -w
#
# File: TopologicalPharmacophoreAtomPairsFingerprints.pl
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
use AtomTypes::FunctionalClassAtomTypes;
use Fingerprints::TopologicalPharmacophoreAtomPairsFingerprints;

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
    GenerateTopologicalPharmacophoreAtomPairsFingerprints($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Generate fingerprints for a SD file...
#
sub GenerateTopologicalPharmacophoreAtomPairsFingerprints {
  my($FileIndex) = @_;
  my($CmpdCount, $IgnoredCmpdCount, $SDFile, $MoleculeFileIO, $Molecule, $TopologicalPharmacophoreAtomPairsFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO, $SetupOutputFiles);

  $SDFile = $SDFilesList[$FileIndex];

  ($NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO) = (undef) x 3;
  $SetupOutputFiles = 1;

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

    $TopologicalPharmacophoreAtomPairsFingerprints = GenerateMoleculeFingerprints($Molecule);
    if (!$TopologicalPharmacophoreAtomPairsFingerprints) {
      $IgnoredCmpdCount++;
      ProcessIgnoredCompound('FingerprintsGenerationFailed', $CmpdCount, $Molecule);
      next COMPOUND;
    }

    if ($SetupOutputFiles) {
      $SetupOutputFiles = 0;
      SetupFingerprintsLabelValueIDs($TopologicalPharmacophoreAtomPairsFingerprints);
      ($NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO) = SetupAndOpenOutputFiles($FileIndex);
    }

    WriteDataToOutputFiles($FileIndex, $CmpdCount, $Molecule, $TopologicalPharmacophoreAtomPairsFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO);
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

# Append atom pair value IDs to fingerprint label...
#
sub SetupFingerprintsLabelValueIDs {
  my($TopologicalPharmacophoreAtomPairsFingerprints) = @_;

  if ($OptionsInfo{AtomPairsSetSizeToUse} =~ /^ArbitrarySize$/i ||
      $OptionsInfo{FingerprintsLabelMode} !~ /^FingerprintsLabelWithIDs$/i) {
    return;
  }

  $OptionsInfo{FingerprintsLabel} .= "; Value IDs: " . $TopologicalPharmacophoreAtomPairsFingerprints->GetFingerprintsVector->GetValueIDsString();
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
  my($FileIndex, $CmpdCount, $Molecule, $TopologicalPharmacophoreAtomPairsFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO) = @_;
  my($DataFieldLabelAndValuesRef);

  $DataFieldLabelAndValuesRef = undef;
  if ($NewFPTextFileIO || $NewFPFileIO) {
    $DataFieldLabelAndValuesRef = $Molecule->GetDataFieldLabelAndValues();
  }

  if ($NewFPSDFileIO) {
    my($CmpdString);

    $CmpdString = $Molecule->GetInputMoleculeString();
    $NewFPSDFileIO->WriteFingerprints($TopologicalPharmacophoreAtomPairsFingerprints, $CmpdString);
  }

  if ($NewFPTextFileIO) {
    my($ColValuesRef);

    $ColValuesRef = SetupFPTextFileCoulmnValues($FileIndex, $CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPTextFileIO->WriteFingerprints($TopologicalPharmacophoreAtomPairsFingerprints, $ColValuesRef);
  }

  if ($NewFPFileIO) {
    my($CompoundID);

    $CompoundID = SetupCmpdIDForOutputFiles($CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPFileIO->WriteFingerprints($TopologicalPharmacophoreAtomPairsFingerprints, $CompoundID);
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
  my($TopologicalPharmacophoreAtomPairsFingerprints);

  if ($OptionsInfo{KeepLargestComponent}) {
    $Molecule->KeepLargestComponent();
  }
  if (!$Molecule->DetectRings()) {
    return undef;
  }
  $Molecule->SetAromaticityModel($OptionsInfo{AromaticityModel});
  $Molecule->DetectAromaticity();

  if ($OptionsInfo{FuzzifyAtomPairsCount}) {
    $TopologicalPharmacophoreAtomPairsFingerprints = new Fingerprints::TopologicalPharmacophoreAtomPairsFingerprints('Molecule' => $Molecule, 'AtomPairsSetSizeToUse' => $OptionsInfo{AtomPairsSetSizeToUse}, 'MinDistance' => $OptionsInfo{MinDistance},  'MaxDistance' => $OptionsInfo{MaxDistance}, 'AtomTypesToUse' => \@{$OptionsInfo{AtomTypesToUse}}, , 'NormalizationMethodology' => $OptionsInfo{NormalizationMethodology}, , 'ValuesPrecision' => $OptionsInfo{ValuesPrecision}, 'FuzzifyAtomPairsCount' => $OptionsInfo{FuzzifyAtomPairsCount}, 'FuzzificationMode' =>  $OptionsInfo{FuzzificationMode}, 'FuzzificationMethodology' => $OptionsInfo{FuzzificationMethodology}, 'FuzzFactor' => $OptionsInfo{FuzzFactor});
  }
  else {
    $TopologicalPharmacophoreAtomPairsFingerprints = new Fingerprints::TopologicalPharmacophoreAtomPairsFingerprints('Molecule' => $Molecule, 'AtomPairsSetSizeToUse' => $OptionsInfo{AtomPairsSetSizeToUse}, 'MinDistance' => $OptionsInfo{MinDistance},  'MaxDistance' => $OptionsInfo{MaxDistance}, 'AtomTypesToUse' => \@{$OptionsInfo{AtomTypesToUse}}, 'NormalizationMethodology' => $OptionsInfo{NormalizationMethodology}, 'ValuesPrecision' => $OptionsInfo{ValuesPrecision});
  }

  # Set atom types weights...
  if ($OptionsInfo{UseAtomTypesWeight}) {
    $TopologicalPharmacophoreAtomPairsFingerprints->SetAtomTypesWeight(%{$OptionsInfo{AtomTypesWeight}});
  }

  # Generate fingerprints...
  $TopologicalPharmacophoreAtomPairsFingerprints->GenerateFingerprints();

  # Make sure fingerprints generation is successful...
  if (!$TopologicalPharmacophoreAtomPairsFingerprints->IsFingerprintsGenerationSuccessful()) {
    return undef;
  }

  return $TopologicalPharmacophoreAtomPairsFingerprints;
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
      $OutFileRoot = "${FileName}TopologicalPharmacophoreAtomPairsFP";
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

  ProcessAtomTypesToUseOption();
  ProcessAtomTypesWeightOption();

  $OptionsInfo{AromaticityModel} = $Options{aromaticitymodel};

  $OptionsInfo{AtomPairsSetSizeToUse} = $Options{atompairssetsizetouse};

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

  $OptionsInfo{FingerprintsLabelMode} = $Options{fingerprintslabelmode};
  $OptionsInfo{FingerprintsLabel} = $Options{fingerprintslabel} ? $Options{fingerprintslabel} : 'TopologicalPharmacophoreAtomPairsFingerprints';

  $OptionsInfo{FuzzifyAtomPairsCount} = ($Options{fuzzifyatompairscount} =~ /^Yes$/i) ? 1 : 0;
  $OptionsInfo{FuzzificationMode} = $Options{fuzzificationmode};
  $OptionsInfo{FuzzificationMethodology} = $Options{fuzzificationmethodology};
  $OptionsInfo{FuzzFactor} = $Options{fuzzfactor};

  $OptionsInfo{KeepLargestComponent} = ($Options{keeplargestcomponent} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{MinDistance} = $Options{mindistance};
  $OptionsInfo{MaxDistance} = $Options{maxdistance};

  $OptionsInfo{NormalizationMethodology} = $Options{normalizationmethodology};

  $OptionsInfo{Output} = $Options{output};
  $OptionsInfo{SDOutput} = ($Options{output} =~ /^(SD|All)$/i) ? 1 : 0;
  $OptionsInfo{FPOutput} = ($Options{output} =~ /^(FP|All)$/i) ? 1 : 0;
  $OptionsInfo{TextOutput} = ($Options{output} =~ /^(Text|All)$/i) ? 1 : 0;

  $OptionsInfo{OutDelim} = $Options{outdelim};
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? 1 : 0;
  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : 0;

  $OptionsInfo{ValuesPrecision} = $Options{valuesprecision};

  # Setup default vector string format...
  my($VectorStringFormat);
  $VectorStringFormat = '';

  if ($Options{vectorstringformat}) {
    $VectorStringFormat = $Options{vectorstringformat};

    if ($Options{atompairssetsizetouse} =~ /^ArbitrarySize$/i && $VectorStringFormat =~ /^ValuesString$/i) {
      die "Error: The value specified, $Options{vectorstringformat}, for option \"-v, --VectorStringFormat\" is not valid for $Options{atompairssetsizetouse} value of \"--AtomPairsSetSizeToUse\" option. Allowed values: IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString or ValuesAndIDsPairsString\n";
    }
  }
  else {
    $VectorStringFormat = ($Options{atompairssetsizetouse} =~ /^FixedSize$/) ? "ValuesString" : "IDsAndValuesString";
  }
  $OptionsInfo{VectorStringFormat} = $VectorStringFormat;
}

# Process atom type to use option...
#
sub ProcessAtomTypesToUseOption {
  my($AtomType, $SpecifiedAtomTypesToUse, @AtomTypesWords);

  @{$OptionsInfo{AtomTypesToUse}} = ();
  if (IsEmpty($Options{atomtypestouse})) {
    die "Error: Atom types value specified using \"-a, --AtomTypesToUse\" option is empty\n";
  }

  $SpecifiedAtomTypesToUse = $Options{atomtypestouse};
  $SpecifiedAtomTypesToUse =~ s/ //g;
  @AtomTypesWords = split /\,/, $SpecifiedAtomTypesToUse;

  for $AtomType (@AtomTypesWords) {
    if (!AtomTypes::FunctionalClassAtomTypes::IsFunctionalClassAvailable($AtomType)) {
      die "Error: Atomic type specified, $AtomType, using \"-a, --AtomTypesToUse\" option is not valid...\n ";
    }
    push @{$OptionsInfo{AtomTypesToUse}}, $AtomType;
  }
}

# Process atom types weight option...
#
sub ProcessAtomTypesWeightOption {
  my($Index, $AtomType, $AtomTypeWeight, $SpecifiedAtomTypesWeight, @AtomTypesWeightsPairs);

  %{$OptionsInfo{AtomTypesWeight}} = ();

  if (IsEmpty($Options{atomtypesweight})) {
    die "Error: Atom types weight value specified using \"--AtomTypesWeight\" option is empty\n";
  }
  $OptionsInfo{UseAtomTypesWeight} = ($Options{atomtypesweight} =~ /^None$/i) ? 0 : 1;
  if (!$OptionsInfo{UseAtomTypesWeight}) {
    return;
  }

  # Process specified atom type/weight pairs...
  $SpecifiedAtomTypesWeight = $Options{atomtypesweight};
  $SpecifiedAtomTypesWeight =~ s/ //g;
  @AtomTypesWeightsPairs = split /\,/, $SpecifiedAtomTypesWeight;

  if (@AtomTypesWeightsPairs % 2) {
    die "Error: Invalid number of values specified using \"--AtomTypesWeight\" option: It must contain even number of values.\n";
  }

  for ($Index = 0; $Index < @AtomTypesWeightsPairs; $Index += 2) {
    $AtomType = $AtomTypesWeightsPairs[$Index]; $AtomTypeWeight = $AtomTypesWeightsPairs[$Index + 1];
    if (!AtomTypes::FunctionalClassAtomTypes::IsFunctionalClassAvailable($AtomType)) {
      die "Error: Atom type specified, $AtomType, using \"--AtomTypesWeight\" option is not valid\n ";
    }
    if (!(IsFloat($AtomTypeWeight) && $AtomTypeWeight >= 0)) {
      die "Error: Atom type weight specified, $AtomTypeWeight, using option \"--AtomTypesWeight\" is not valid. Allowed values: real numbers >= 0 \n";
    }
    $OptionsInfo{AtomTypesWeight}{$AtomType} = $AtomTypeWeight;
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{aromaticitymodel} = 'MayaChemToolsAromaticityModel';

  $Options{atompairssetsizetouse} = 'ArbitrarySize';

  $Options{atomtypestouse} = 'HBD,HBA,PI,NI,H';
  $Options{atomtypesweight} = 'None';

  $Options{compoundidmode} = 'LabelPrefix';
  $Options{compoundidlabel} = 'CompoundID';
  $Options{datafieldsmode} = 'CompoundID';

  $Options{filter} = 'Yes';

  $Options{fingerprintslabelmode} = 'FingerprintsLabelOnly';

  $Options{fuzzifyatompairscount} = 'No';
  $Options{fuzzificationmode} = 'AfterNormalization';
  $Options{fuzzificationmethodology} = 'FuzzyBinning';
  $Options{fuzzfactor} = 0.15;

  $Options{keeplargestcomponent} = 'Yes';

  $Options{mindistance} = 1;
  $Options{maxdistance} = 10;

  $Options{normalizationmethodology} = 'None';

  $Options{output} = 'text';
  $Options{outdelim} = 'comma';
  $Options{quote} = 'yes';

  $Options{valuesprecision} = 2;

  $Options{vectorstringformat} = '';

  if (!GetOptions(\%Options, "aromaticitymodel=s", "atompairssetsizetouse=s", "atomtypestouse|a=s", "atomtypesweight=s", "compoundid=s", "compoundidlabel=s", "compoundidmode=s", "datafields=s", "datafieldsmode|d=s", "filter|f=s", "fingerprintslabelmode=s", "fingerprintslabel=s", "fuzzifyatompairscount=s", "fuzzificationmode=s", "fuzzificationmethodology=s", "fuzzfactor=s", "help|h", "keeplargestcomponent|k=s",  "mindistance=s", "maxdistance=s", "normalizationmethodology|n=s", "outdelim=s", "output=s", "overwrite|o", "quote|q=s", "root|r=s", "valuesprecision=s", "vectorstringformat|v=s", "workingdir|w=s")) {
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
  if ($Options{atompairssetsizetouse} !~ /^(ArbitrarySize|FixedSize)$/i) {
    die "Error: The value specified, $Options{atompairssetsizetouse}, for option \"--AtomPairsSetSizeToUse\" is not valid. Allowed values: ArbitrarySize or FixedSize\n";
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
  if ($Options{fuzzifyatompairscount} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{fuzzifyatompairscount}, for option \"--FuzzifyAtomPairsCount\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{fuzzificationmode} !~ /^(BeforeNormalization|AfterNormalization)$/i) {
    die "Error: The value specified, $Options{fuzzificationmode}, for option \"--FuzzificationMode\" is not valid. Allowed values: BeforeNormalization or AfterNormalization\n";
  }
  if ($Options{fuzzificationmethodology} !~ /^(FuzzyBinning|FuzzyBinSmoothing)$/i) {
    die "Error: The value specified, $Options{fuzzificationmethodology}, for option \"--FuzzificationMethodology\" is not valid. Allowed values: FuzzyBinning or FuzzyBinSmoothing\n";
  }
  if (!IsFloat($Options{fuzzfactor})) {
    die "Error: The value specified, $Options{fuzzfactor}, for option \"--FuzzFactor\" is not valid. Allowed values: real numbers >= 0 \n";
  }
  if ($Options{fuzzificationmethodology} !~ /^FuzzyBinning$/i) {
    if (!($Options{fuzzfactor} >=0 && $Options{fuzzfactor} <= 1.0)) {
      die "Error: The value specified, $Options{fuzzfactor}, for option \"--FuzzFactor\" during FuzzyBinning \"--FuzzificationMethodology\" is not valid. Allowed values: >= 0 and <= 1 \n";
    }
  }
  elsif ($Options{fuzzificationmethodology} !~ /^FuzzyBinSmoothing$/i) {
    if (!($Options{fuzzfactor} >=0 && $Options{fuzzfactor} <= 0.5)) {
      die "Error: The value specified, $Options{fuzzfactor}, for option \"--FuzzFactor\" during FuzzyBinSmoothing \"--FuzzificationMethodology\" is not valid. Allowed values: >= 0 and <= 0.5 \n";
    }
  }
  if ($Options{keeplargestcomponent} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{keeplargestcomponent}, for option \"-k, --KeepLargestComponent\" is not valid. Allowed values: Yes or No\n";
  }
  if (!IsInteger($Options{mindistance})) {
    die "Error: The value specified, $Options{mindistance}, for option \"--MinDistance\" is not valid. Allowed values: >= 0 \n";
  }
  if (!IsPositiveInteger($Options{maxdistance})) {
    die "Error: The value specified, $Options{maxdistance}, for option \"--MaxDistance\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{mindistance} > $Options{maxdistance}) {
    die "Error: The value specified, specified, $Options{mindistance}, for option \"--MinDistance\" must be less than the value specified, $Options{maxdistance}, for option \"--MaxDistance\" \n";
  }
  if ($Options{normalizationmethodology} !~ /^(None|ByHeavyAtomsCount|ByAtomTypesCount)$/i) {
    die "Error: The value specified, $Options{normalizationmethodology}, for option \"--NormalizationMethodology\" is not valid. Allowed values: None, ByHeavyAtomsCount, or ByAtomTypesCount\n";
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
  if (!IsPositiveInteger($Options{valuesprecision})) {
    die "Error: The value specified, $Options{valuesprecision}, for option \"--ValuesPrecision\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{vectorstringformat} && $Options{vectorstringformat} !~ /^(ValuesString|IDsAndValuesString|IDsAndValuesPairsString|ValuesAndIDsString|ValuesAndIDsPairsString)$/i) {
    die "Error: The value specified, $Options{vectorstringformat}, for option \"-v, --VectorStringFormat\" is not valid. Allowed values: ValuesString, IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString or ValuesAndIDsPairsString\n";
  }
}

__END__

=head1 NAME

TopologicalPharmacophoreAtomPairsFingerprints.pl - Generate topological pharmacophore atom pairs fingerprints for SD files

=head1 SYNOPSIS

TopologicalPharmacophoreAtomPairsFingerprints.pl SDFile(s)...

TopologicalPharmacophoreAtomPairsFingerprints.pl  [B<--AromaticityModel> I<AromaticityModelType>]
[B<--AtomPairsSetSizeToUse> I<ArbitrarySize | FixedSize>]
[B<-a, --AtomTypesToUse> I<"AtomType1, AtomType2...">]
[B<--AtomTypesWeight> I<"AtomType1, Weight1, AtomType2, Weight2...">]
[B<--CompoundID> I<DataFieldName or LabelPrefixString>] [B<--CompoundIDLabel> I<text>]
[B<--CompoundIDMode>] [B<--DataFields> I<"FieldLabel1, FieldLabel2,...">]
[B<-d, --DataFieldsMode> I<All | Common | Specify | CompoundID>] [B<-f, --Filter> I<Yes | No>]
[B<--FingerprintsLabelMode> I<FingerprintsLabelOnly | FingerprintsLabelWithIDs>] [B<--FingerprintsLabel> I<text>]
[B<--FuzzifyAtomPairsCount> I<Yes | No>] [B<--FuzzificationMode> I<FuzzyBinning | FuzzyBinSmoothing>]
[B<--FuzzificationMethodology> I<FuzzyBinning | FuzzyBinSmoothing>] [B<--FuzzFactor> I<number>]
[B<-h, --help>] [B<-k, --KeepLargestComponent> I<Yes | No>] [B<--MinDistance> I<number>]
[B<--MaxDistance> I<number>] [B<-n, --NormalizationMethodology> I<None | ByHeavyAtomsCount | ByAtomTypesCount>]
[B<--OutDelim> I<comma | tab | semicolon>] [B<--output> I<SD | FP | text | all>] [B<-o, --overwrite>]
[B<-q, --quote> I<Yes | No>] [B<-r, --root> I<RootName>] [B<--ValuesPrecision> I<number>]
[B<-v, --VectorStringFormat> I<ValuesString, IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>]
[B<-w, --WorkingDir> dirname] SDFile(s)...

=head1 DESCRIPTION

Generate topological pharmacophore atom pairs fingerprints [ Ref 60-62, Ref 65, Ref 68 ] for
I<SDFile(s)> and create appropriate SD, FP or CSV/TSV text file(s) containing fingerprints vector
strings corresponding to molecular fingerprints.

Multiple SDFile names are separated by spaces. The valid file extensions are I<.sdf>
and I<.sd>. All other file names are ignored. All the SD files in a current directory
can be specified either by I<*.sdf> or the current directory name.

Based on the values specified for B<--AtomTypesToUse>, pharmacophore atom types are
assigned to all non-hydrogen atoms in a molecule and a distance matrix is generated.
A pharmacophore atom pairs basis set is initialized for all unique possible pairs within
B<--MinDistance> and B<--MaxDistance> range.

    Let:

    P = Valid pharmacophore atom type

    Px = Pharmacophore atom type x
    Py = Pharmacophore atom type y

    Dmin = Minimum distance corresponding to number of bonds between
           two atoms
    Dmax = Maximum distance corresponding to number of bonds between
           two atoms
    D = Distance corresponding to number of bonds between two atoms

    Px-Dn-Py = Pharmacophore atom pair ID for atom types Px and Py at
               distance Dn

    P = Number of pharmacophore atom types to consider
    PPDn = Number of possible unique pharmacophore atom pairs at a distance Dn

    PPT = Total number of possible pharmacophore atom pairs at all distances
          between Dmin and Dmax

    Then:

    PPD =  (P * (P - 1))/2 + P

    PPT = ((Dmax - Dmin) + 1) * ((P * (P - 1))/2 + P)
        = ((Dmax - Dmin) + 1) * PPD

    So for default values of Dmin = 1, Dmax = 10 and P = 5,

    PPD =  (5 * (5 - 1))/2 + 5 = 15
    PPT = ((10 - 1) + 1) * 15 = 150

    The pharmacophore atom pairs bais set includes 150 values.

    The atom pair IDs correspond to:

    Px-Dn-Py = Pharmacophore atom pair ID for atom types Px and Py at
               distance Dn

    For example: H-D1-H, H-D2-HBA, PI-D5-PI and so on

Using distance matrix and pharmacohore atom types, occurrence of unique pharmacohore atom
pairs is counted. The contribution of each atom type to atom pair interaction is optionally
weighted by specified B<--AtomTypesWeight> before assigning its count to appropriate distance
bin. Based on B<--NormalizationMethodology> option, pharmacophore atom pairs count is optionally
normalized. Additionally, pharmacohore atom pairs count is optionally fuzzified before or after
the normalization controlled by values of B<--FuzzifyAtomPairsCount>, B<--FuzzificationMode>,
B<--FuzzificationMethodology> and B<--FuzzFactor> options.

The final pharmacophore atom pairs count along with atom pair identifiers involving all non-hydrogen
atoms, with optional normalization and fuzzification, constitute pharmacophore topological atom pairs
fingerprints of the molecule.

For I<ArbitrarySize> value of B<--AtomPairsSetSizeToUse> option, the fingerprint vector correspond to
only those topological pharmacophore atom pairs which are present and have non-zero count. However,
for I<FixedSize> value of B<--AtomPairsSetSizeToUse> option, the fingerprint vector contains all possible
valid topological pharmacophore atom pairs with both zero and non-zero count values.

Example of I<SD> file containing topological pharmacophore atom pairs fingerprints string data:

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

    >  <TopologicalPharmacophoreAtomPairsFingerprints>
    FingerprintsVector;TopologicalPharmacophoreAtomPairs:ArbitrarySize:Min
    Distance1:MaxDistance10;54;NumericalValues;IDsAndValuesString;H-D1-H H
    -D1-NI HBA-D1-NI HBD-D1-NI H-D2-H H-D2-HBA H-D2-HBD HBA-D2-HBA HBA-D2-
    HBD H-D3-H H-D3-HBA H-D3-HBD H-D3-NI HBA-D3-NI HBD-D3-NI H-D4-H H-D...;
    18 1 2 1 22 12 8 1 2 18 6 3 1 1 1 22 13 6 5 7 2 28 9 5 1 1 1 36 16 10 3
    4 1 37 10 8 1 35 10 9 3 3 1 28 7 7 4 18 16 12 5 1 2 1

    $$$$
    ... ...
    ... ...

Example of I<FP> file containing topological pharmacophore atom pairs fingerprints string data:

    #
    # Package = MayaChemTools 7.4
    # Release Date = Oct 21, 2010
    #
    # TimeStamp = Fri Mar 11 15:32:48 2011
    #
    # FingerprintsStringType = FingerprintsVector
    #
    # Description = TopologicalPharmacophoreAtomPairs:ArbitrarySize:MinDistance1:MaxDistance10
    # VectorStringFormat = IDsAndValuesString
    # VectorValuesType = NumericalValues
    #
    Cmpd1 54;H-D1-H H-D1-NI HBA-D1-NI HBD-D1-NI H-D2-H H-D2-HBA...;18 1 2...
    Cmpd2 61;H-D1-H H-D1-NI HBA-D1-NI HBD-D1-NI H-D2-H H-D2-HBA...;5 1 2 ...
    ... ...
    ... ..

Example of CSV I<Text> file containing topological pharmacophore atom pairs fingerprints string data:

    "CompoundID","TopologicalPharmacophoreAtomPairsFingerprints"
    "Cmpd1","FingerprintsVector;TopologicalPharmacophoreAtomPairs:Arbitrary
    Size:MinDistance1:MaxDistance10;54;NumericalValues;IDsAndValuesString;H
    -D1-H H-D1-NI HBA-D1-NI HBD-D1-NI H-D2-H H-D2-HBA H-D2-HBD HBA-D2-HBA H
    BA-D2-HBD H-D3-H H-D3-HBA H-D3-HBD H-D3-NI HBA-D3-NI HBD-D3-NI H-D4...;
    18 1 2 1 22 12 8 1 2 18 6 3 1 1 1 22 13 6 5 7 2 28 9 5 1 1 1 36 16 10 3
    4 1 37 10 8 1 35 10 9 3 3 1 28 7 7 4 18 16 12 5 1 2 1"
    ... ...
    ... ...

The current release of MayaChemTools generates the following types of topological pharmacophore
atom pairs fingerprints vector strings:

    FingerprintsVector;TopologicalPharmacophoreAtomPairs:ArbitrarySize:Min
    Distance1:MaxDistance10;54;NumericalValues;IDsAndValuesString;H-D1-H H
    -D1-NI HBA-D1-NI HBD-D1-NI H-D2-H H-D2-HBA H-D2-HBD HBA-D2-HBA HBA-D2-
    HBD H-D3-H H-D3-HBA H-D3-HBD H-D3-NI HBA-D3-NI HBD-D3-NI H-D4-H H-D4-H
    BA H-D4-HBD HBA-D4-HBA HBA-D4-HBD HBD-D4-HBD H-D5-H H-D5-HBA H-D5-...;
    18 1 2 1 22 12 8 1 2 18 6 3 1 1 1 22 13 6 5 7 2 28 9 5 1 1 1 36 16 10
    3 4 1 37 10 8 1 35 10 9 3 3 1 28 7 7 4 18 16 12 5 1 2 1

    FingerprintsVector;TopologicalPharmacophoreAtomPairs:FixedSize:MinDist
    ance1:MaxDistance10;150;OrderedNumericalValues;ValuesString;18 0 0 1 0
    0 0 2 0 0 1 0 0 0 0 22 12 8 0 0 1 2 0 0 0 0 0 0 0 0 18 6 3 1 0 0 0 1
    0 0 1 0 0 0 0 22 13 6 0 0 5 7 0 0 2 0 0 0 0 0 28 9 5 1 0 0 0 1 0 0 1 0
    0 0 0 36 16 10 0 0 3 4 0 0 1 0 0 0 0 0 37 10 8 0 0 0 0 1 0 0 0 0 0 0
    0 35 10 9 0 0 3 3 0 0 1 0 0 0 0 0 28 7 7 4 0 0 0 0 0 0 0 0 0 0 0 18...

    FingerprintsVector;TopologicalPharmacophoreAtomPairs:FixedSize:MinDist
    ance1:MaxDistance10;150;OrderedNumericalValues;IDsAndValuesString;H-D1
    -H H-D1-HBA H-D1-HBD H-D1-NI H-D1-PI HBA-D1-HBA HBA-D1-HBD HBA-D1-NI H
    BA-D1-PI HBD-D1-HBD HBD-D1-NI HBD-D1-PI NI-D1-NI NI-D1-PI PI-D1-PI H-D
    2-H H-D2-HBA H-D2-HBD H-D2-NI H-D2-PI HBA-D2-HBA HBA-D2-HBD HBA-D2...;
    18 0 0 1 0 0 0 2 0 0 1 0 0 0 0 22 12 8 0 0 1 2 0 0 0 0 0 0 0 0 18 6 3
    1 0 0 0 1 0 0 1 0 0 0 0 22 13 6 0 0 5 7 0 0 2 0 0 0 0 0 28 9 5 1 0 0 0
    1 0 0 1 0 0 0 0 36 16 10 0 0 3 4 0 0 1 0 0 0 0


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

=item B<--AtomPairsSetSizeToUse> I<ArbitrarySize | FixedSize>

Atom pairs set size to use during generation of topological pharmacophore atom pairs
fingerprints.

Possible values: I<ArbitrarySize | FixedSize>; Default value: I<ArbitrarySize>.

For I<ArbitrarySize> value of B<--AtomPairsSetSizeToUse> option, the fingerprint vector
correspond to only those topological pharmacophore atom pairs which are present and
have non-zero count. However, for I<FixedSize> value of B<--AtomPairsSetSizeToUse>
option, the fingerprint vector contains all possible valid topological pharmacophore atom
pairs with both zero and non-zero count values.

=item B<-a, --AtomTypesToUse> I<"AtomType1,AtomType2,...">

Pharmacophore atom types to use during generation of topological phramacophore
atom pairs. It's a list of comma separated valid pharmacophore atom types.

Possible values for pharmacophore atom types are: I<Ar, CA, H, HBA, HBD, Hal, NI, PI, RA>.
Default value [ Ref 60-62 ] : I<HBD,HBA,PI,NI,H>.

The pharmacophore atom types abbreviations correspond to:

    HBD: HydrogenBondDonor
    HBA: HydrogenBondAcceptor
    PI :  PositivelyIonizable
    NI : NegativelyIonizable
    Ar : Aromatic
    Hal : Halogen
    H : Hydrophobic
    RA : RingAtom
    CA : ChainAtom

I<AtomTypes::FunctionalClassAtomTypes> module is used to assign pharmacophore atom
types. It uses following definitions [ Ref 60-61, Ref 65-66 ]:

    HydrogenBondDonor: NH, NH2, OH
    HydrogenBondAcceptor: N[!H], O
    PositivelyIonizable: +, NH2
    NegativelyIonizable: -, C(=O)OH, S(=O)OH, P(=O)OH

=item B<--AtomTypesWeight> I<"AtomType1,Weight1,AtomType2,Weight2...">

Weights of specified pharmacophore atom types to use during calculation of their contribution
to atom pair count. Default value: I<None>. Valid values: real numbers greater than 0. In general
it's comma delimited list of valid atom type and its weight.

The weight values allow to increase the importance of specific pharmacophore atom type
in the generated fingerprints. A weight value of 0 for an atom type eliminates its contribution to
atom pair count where as weight value of 2 doubles its contribution.

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

=item B<--FingerprintsLabelMode> I<FingerprintsLabelOnly | FingerprintsLabelWithIDs>

Specify how fingerprints label is generated in conjunction with B<--FingerprintsLabel> option value:
use fingerprints label generated only by B<--FingerprintsLabel> option value or append topological
atom pair count value IDs to B<--FingerprintsLabel> option value.

Possible values: I<FingerprintsLabelOnly | FingerprintsLabelWithIDs>. Default value:
I<FingerprintsLabelOnly>.

Topological atom pairs IDs appended to B<--FingerprintsLabel> value during I<FingerprintsLabelWithIDs>
values of B<--FingerprintsLabelMode>  correspond to atom pair count values in fingerprint vector string.

I<FingerprintsLabelWithIDs> value of B<--FingerprintsLabelMode> is ignored during I<ArbitrarySize> value
of B<--AtomPairsSetSizeToUse> option and topological atom pairs IDs not appended to the label.

=item B<--FingerprintsLabel> I<text>

SD data label or text file column label to use for fingerprints string in output SD or
CSV/TSV text file(s) specified by B<--output>. Default value: I<TopologicalPharmacophoreAtomPairsFingerprints>.

=item B<--FuzzifyAtomPairsCount> I<Yes | No>

To fuzzify or not to fuzzify atom pairs count. Possible values: I<Yes or No>. Default value:
I<No>.

=item B<--FuzzificationMode> I<BeforeNormalization | AfterNormalization>

When to fuzzify atom pairs count. Possible values: I<BeforeNormalization | AfterNormalizationYes>.
Default value: I<AfterNormalization>.

=item B<--FuzzificationMethodology> I<FuzzyBinning | FuzzyBinSmoothing>

How to fuzzify atom pairs count. Possible values: I<FuzzyBinning | FuzzyBinSmoothing>.
Default value: I<FuzzyBinning>.

In conjunction with values for options B<--FuzzifyAtomPairsCount>, B<--FuzzificationMode> and
B<--FuzzFactor>, B<--FuzzificationMethodology> option is used to fuzzify pharmacophore atom
pairs count.

Let:

    Px = Pharmacophore atom type x
    Py = Pharmacophore atom type y
    PPxy = Pharmacophore atom pair between atom type Px and Py

    PPxyDn = Pharmacophore atom pairs count between atom type Px and Py
             at distance Dn
    PPxyDn-1 = Pharmacophore atom pairs count between atom type Px and Py
               at distance Dn - 1
    PPxyDn+1 = Pharmacophore atom pairs count between atom type Px and Py
               at distance Dn + 1

    FF = FuzzFactor for FuzzyBinning and FuzzyBinSmoothing

Then:

For I<FuzzyBinning>:

    PPxyDn = PPxyDn (Unchanged)

    PPxyDn-1 = PPxyDn-1 + PPxyDn * FF
    PPxyDn+1 = PPxyDn+1 + PPxyDn * FF

For I<FuzzyBinSmoothing>:

    PPxyDn = PPxyDn - PPxyDn * 2FF for Dmin < Dn < Dmax
    PPxyDn = PPxyDn - PPxyDn * FF for Dn = Dmin or Dmax

    PPxyDn-1 = PPxyDn-1 + PPxyDn * FF
    PPxyDn+1 = PPxyDn+1 + PPxyDn * FF

In both fuzzification schemes, a value of 0 for FF implies no fuzzification of occurrence counts.
A value of 1 during I<FuzzyBinning> corresponds to maximum fuzzification of occurrence counts;
however, a value of 1 during I<FuzzyBinSmoothing> ends up completely distributing the value over
the previous and next distance bins.

So for default value of B<--FuzzFactor> (FF) 0.15, the occurrence count of pharmacohore atom pairs
at distance Dn during FuzzyBinning is left unchanged and the counts at distances Dn -1 and Dn + 1
are incremented by PPxyDn * 0.15.

And during I<FuzzyBinSmoothing> the occurrence counts at Distance Dn is scaled back using multiplicative
factor of (1 - 2*0.15) and the occurrence counts at distances Dn -1 and Dn + 1 are incremented by
PPxyDn * 0.15. In otherwords, occurrence bin count is smoothed out by distributing it over the
previous and next distance value.

=item B<--FuzzFactor> I<number>

Specify by how much to fuzzify atom pairs count. Default value: I<0.15>. Valid values: For
I<FuzzyBinning> value of B<--FuzzificationMethodology> option: I<between 0 and 1.0>; For
I<FuzzyBinSmoothing> value of B<--FuzzificationMethodology> option: I<between 0 and 0.5>.

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

Minimum bond distance between atom pairs for generating topological pharmacophore atom
pairs. Default value: I<1>. Valid values: positive integers including 0 and less than B<--MaxDistance>.

=item B<--MaxDistance> I<number>

Maximum bond distance between atom pairs for generating topological pharmacophore atom
 pairs. Default value: I<10>. Valid values: positive integers and greater than B<--MinDistance>.

=item B<-n, --NormalizationMethodology> I<None | ByHeavyAtomsCount | ByAtomTypesCount>

Normalization methodology to use for scaling the occurrence count of pharmacophore atom
pairs within specified distance range. Possible values: I<None, ByHeavyAtomsCount or
ByAtomTypesCount>. Default value: I<None>.

=item B<--OutDelim> I<comma | tab | semicolon>

Delimiter for output CSV/TSV text file(s). Possible values: I<comma, tab, or semicolon>
Default value: I<comma>.

=item B<--output> I<SD | FP | text | all>

Type of output files to generate. Possible values: I<SD, FP, text, or all>. Default value: I<text>.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-q, --quote> I<Yes | No>

Put quote around column values in output CSV/TSV text file(s). Possible values:
I<Yes or No>. Default value: I<Yes>

=item B<-r, --root> I<RootName>

New file name is generated using the root: <Root>.<Ext>. Default for new file names:
<SDFileName><TopologicalPharmacophoreAtomPairsFP>.<Ext>. The file type determines <Ext> value.
The sdf, fpf, csv, and tsv <Ext> values are used for SD, FP, comma/semicolon, and tab
delimited text files, respectively.This option is ignored for multiple input files.

=item B<--ValuesPrecision> I<number>

Precision of atom pairs count real values which might be generated after normalization
or fuzzification. Default value: up to I<2> decimal places. Valid values: positive integers.

=item B<-v, --VectorStringFormat> I<ValuesString, IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>

Format of fingerprints vector string data in output SD, FP or CSV/TSV text file(s) specified by
B<--output> option. Possible values: I<ValuesString, IDsAndValuesString | IDsAndValuesPairsString
| ValuesAndIDsString | ValuesAndIDsPairsString>.

Default value during I<FixedSize> value of B<--AtomPairsSetSizeToUse> option: I<ValuesString>. Default
value during I<ArbitrarySize> value of B<--AtomPairsSetSizeToUse> option: I<IDsAndValuesString>.

I<ValuesString> option value is not allowed for I<ArbitrarySize> value of B<--AtomPairsSetSizeToUse>
option.

Examples:

    FingerprintsVector;TopologicalPharmacophoreAtomPairs:ArbitrarySize:Min
    Distance1:MaxDistance10;54;NumericalValues;IDsAndValuesString;H-D1-H H
    -D1-NI HBA-D1-NI HBD-D1-NI H-D2-H H-D2-HBA H-D2-HBD HBA-D2-HBA HBA-D2-
    HBD H-D3-H H-D3-HBA H-D3-HBD H-D3-NI HBA-D3-NI HBD-D3-NI H-D4-H H-D4-H
    BA H-D4-HBD HBA-D4-HBA HBA-D4-HBD HBD-D4-HBD H-D5-H H-D5-HBA H-D5-...;
    18 1 2 1 22 12 8 1 2 18 6 3 1 1 1 22 13 6 5 7 2 28 9 5 1 1 1 36 16 10
    3 4 1 37 10 8 1 35 10 9 3 3 1 28 7 7 4 18 16 12 5 1 2 1

    FingerprintsVector;TopologicalPharmacophoreAtomPairs:FixedSize:MinDist
    ance1:MaxDistance10;150;OrderedNumericalValues;ValuesString;18 0 0 1 0
    0 0 2 0 0 1 0 0 0 0 22 12 8 0 0 1 2 0 0 0 0 0 0 0 0 18 6 3 1 0 0 0 1
    0 0 1 0 0 0 0 22 13 6 0 0 5 7 0 0 2 0 0 0 0 0 28 9 5 1 0 0 0 1 0 0 1 0
    0 0 0 36 16 10 0 0 3 4 0 0 1 0 0 0 0 0 37 10 8 0 0 0 0 1 0 0 0 0 0 0
    0 35 10 9 0 0 3 3 0 0 1 0 0 0 0 0 28 7 7 4 0 0 0 0 0 0 0 0 0 0 0 18...

    FingerprintsVector;TopologicalPharmacophoreAtomPairs:FixedSize:MinDist
    ance1:MaxDistance10;150;OrderedNumericalValues;IDsAndValuesString;H-D1
    -H H-D1-HBA H-D1-HBD H-D1-NI H-D1-PI HBA-D1-HBA HBA-D1-HBD HBA-D1-NI H
    BA-D1-PI HBD-D1-HBD HBD-D1-NI HBD-D1-PI NI-D1-NI NI-D1-PI PI-D1-PI H-D
    2-H H-D2-HBA H-D2-HBD H-D2-NI H-D2-PI HBA-D2-HBA HBA-D2-HBD HBA-D2...;
    18 0 0 1 0 0 0 2 0 0 1 0 0 0 0 22 12 8 0 0 1 2 0 0 0 0 0 0 0 0 18 6 3
    1 0 0 0 1 0 0 1 0 0 0 0 22 13 6 0 0 5 7 0 0 2 0 0 0 0 0 28 9 5 1 0 0 0
    1 0 0 1 0 0 0 0 36 16 10 0 0 3 4 0 0 1 0 0 0 0

=item B<-w, --WorkingDir> I<DirName>

Location of working directory. Default value: current directory.

=back

=head1 EXAMPLES

To generate topological pharmacophore atom pairs fingerprints of arbitrary size corresponding to distances
from 1 through 10 using default atom types with no weighting, normalization, and fuzzification
of atom pairs count and create a SampleTPAPFP.csv file containing sequential compound IDs along
with fingerprints vector strings data in ValuesString format, type:

    % TopologicalPharmacophoreAtomPairsFingerprints.pl -r SampleTPAPFP
      -o Sample.sdf

To generate topological pharmacophore atom pairs fingerprints of fixed size corresponding to distances
from 1 through 10 using default atom types with no weighting, normalization, and fuzzification
of atom pairs count and create a SampleTPAPFP.csv file containing sequential compound IDs along
with fingerprints vector strings data in ValuesString format, type:

    % TopologicalPharmacophoreAtomPairsFingerprints.pl
       --AtomPairsSetSizeToUse FixedSize -r SampleTPAPFP-o Sample.sdf

To generate topological pharmacophore atom pairs fingerprints of arbitrary size corresponding to distances
from 1 through 10 using default atom types with no weighting, normalization, and fuzzification
of atom pairs count and create SampleTPAPFP.sdf, SampleTPAPFP.fpf and SampleTPAPFP.csv files containing
sequential compound IDs in CSV file along with fingerprints vector strings data in ValuesString
format, type:

    % TopologicalPharmacophoreAtomPairsFingerprints.pl --output all
      -r SampleTPAPFP -o Sample.sdf

To generate topological pharmacophore atom pairs fingerprints of arbitrary size corresponding to distances
from 1 through 10 using default atom types with no weighting, normalization, and fuzzification
of atom pairs count and create a SampleTPAPFP.csv file containing sequential compound IDs along
with fingerprints vector strings data in IDsAndValuesPairsString format, type:

    % TopologicalPharmacophoreAtomPairsFingerprints.pl --VectorStringFormat
      IDsAndValuesPairsString -r SampleTPAPFP -o Sample.sdf

To generate topological pharmacophore atom pairs fingerprints of arbitrary size corresponding to distances
from 1 through 6 using default atom types with no weighting, normalization, and fuzzification
of atom pairs count and create a SampleTPAPFP.csv file containing sequential compound IDs along
with fingerprints vector strings data in ValuesString format, type:

    % TopologicalPharmacophoreAtomPairsFingerprints.pl --MinDistance 1
      -MaxDistance 6 -r SampleTPAPFP -o Sample.sdf

To generate topological pharmacophore atom pairs fingerprints of arbitrary size corresponding to distances
from 1 through 10 using "HBD,HBA,PI,NI" atom types with double the weighting for "HBD,HBA" and
normalization by HeavyAtomCount but no fuzzification of atom pairs count and create a
SampleTPAPFP.csv file containing sequential compound IDs along with fingerprints vector strings
data in ValuesString format, type:

    % TopologicalPharmacophoreAtomPairsFingerprints.pl --MinDistance 1
      -MaxDistance 10  --AtomTypesToUse "HBD,HBA,PI, NI"  --AtomTypesWeight
      "HBD,2,HBA,2,PI,1,NI,1" --NormalizationMethodology ByHeavyAtomsCount
      --FuzzifyAtomPairsCount No -r SampleTPAPFP -o Sample.sdf

To generate topological pharmacophore atom pairs fingerprints of arbitrary size corresponding to
distances from 1 through 10 using "HBD,HBA,PI,NI,H" atom types with no weighting of atom types and
normalization but with fuzzification of atom pairs count using FuzzyBinning methodology
with FuzzFactor value 0.15 and create a SampleTPAPFP.csv file containing sequential compound
IDs along with fingerprints vector strings data in ValuesString format, type:

    % TopologicalPharmacophoreAtomPairsFingerprints.pl --MinDistance 1
      --MaxDistance 10  --AtomTypesToUse "HBD,HBA,PI, NI,H"  --AtomTypesWeight
      "HBD,1,HBA,1,PI,1,NI,1,H,1" --NormalizationMethodology None
      --FuzzifyAtomPairsCount Yes --FuzzificationMethodology FuzzyBinning
      --FuzzFactor  0.5 -r SampleTPAPFP -o Sample.sdf

To generate topological pharmacophore atom pairs fingerprints of arbitrary size corresponding to distances
distances from 1 through 10 using default atom types with no weighting,
normalization, and fuzzification of atom pairs count and create a SampleTPAPFP.csv
file containing compound ID from molecule name line along with fingerprints vector strings
data, type:

    % TopologicalPharmacophoreAtomPairsFingerprints.pl --DataFieldsMode
      CompoundID -CompoundIDMode MolName -r SampleTPAPFP -o Sample.sdf

To generate topological pharmacophore atom pairs fingerprints of arbitrary size corresponding
to distances from 1 through 10 using default atom types with no weighting,
normalization, and fuzzification of atom pairs count and create a SampleTPAPFP.csv
file containing compound IDs using specified data field along with fingerprints vector strings
data, type:

    % TopologicalPharmacophoreAtomPairsFingerprints.pl --DataFieldsMode
      CompoundID -CompoundIDMode DataField --CompoundID Mol_ID
      -r SampleTPAPFP -o Sample.sdf

To generate topological pharmacophore atom pairs fingerprints of arbitrary size corresponding
to distances from 1 through 10 using default atom types with no weighting,
normalization, and fuzzification of atom pairs count and create a SampleTPAPFP.csv
file containing compound ID using combination of molecule name line and an explicit compound
prefix along with fingerprints vector strings data, type:

    % TopologicalPharmacophoreAtomPairsFingerprints.pl --DataFieldsMode
      CompoundID -CompoundIDMode MolnameOrLabelPrefix
      --CompoundID Cmpd --CompoundIDLabel MolID -r SampleTPAPFP -o Sample.sdf

To generate topological pharmacophore atom pairs fingerprints of arbitrary size corresponding
to distances from 1 through 10 using default atom types with no weighting,
normalization, and fuzzification of atom pairs count and create a SampleTPAPFP.csv
file containing specific data fields columns along with fingerprints vector strings
data, type:

    % TopologicalPharmacophoreAtomPairsFingerprints.pl --DataFieldsMode
      Specify --DataFields Mol_ID -r SampleTPAPFP -o Sample.sdf

To generate topological pharmacophore atom pairs fingerprints of arbitrary size corresponding
to distances from 1 through 10 using default atom types with no weighting,
normalization, and fuzzification of atom pairs count and create a SampleTPAPFP.csv
file containing common data fields columns along with fingerprints vector strings
data, type:

    % TopologicalPharmacophoreAtomPairsFingerprints.pl --DataFieldsMode
      Common -r SampleTPAPFP -o Sample.sdf

To generate topological pharmacophore atom pairs fingerprints of arbitrary size corresponding
to distances from 1 through 10 using default atom types with no weighting,
normalization, and fuzzification of atom pairs count and create SampleTPAPFP.sdf, SampleTPAPFP.fpf,
and SampleTPAPFP.csv files containing all data fields columns in CSV file along with fingerprints
data, type:

    % TopologicalPharmacophoreAtomPairsFingerprints.pl --DataFieldsMode
      All  --output all -r SampleTPAPFP -o Sample.sdf


=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoFingerprintsFiles.pl, SimilarityMatricesFingerprints.pl, AtomNeighborhoodsFingerprints.pl,
ExtendedConnectivityFingerprints.pl, MACCSKeysFingerprints.pl, PathLengthFingerprints.pl,
TopologicalAtomPairsFingerprints.pl, TopologicalAtomTorsionsFingerprints.pl,
TopologicalPharmacophoreAtomTripletsFingerprints.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
