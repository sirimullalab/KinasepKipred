#!/usr/bin/perl -w
#
# File: TopologicalPharmacophoreAtomTripletsFingerprints.pl
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
use Fingerprints::TopologicalPharmacophoreAtomTripletsFingerprints;

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
    GenerateTopologicalPharmacophoreAtomTripletsFingerprints($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Generate fingerprints for a SD file...
#
sub GenerateTopologicalPharmacophoreAtomTripletsFingerprints {
  my($FileIndex) = @_;
  my($CmpdCount, $IgnoredCmpdCount, $SDFile, $MoleculeFileIO, $Molecule, $TopologicalPharmacophoreAtomTripletsFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO, $SetupOutputFiles);

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

    $TopologicalPharmacophoreAtomTripletsFingerprints = GenerateMoleculeFingerprints($Molecule);
    if (!$TopologicalPharmacophoreAtomTripletsFingerprints) {
      $IgnoredCmpdCount++;
      ProcessIgnoredCompound('FingerprintsGenerationFailed', $CmpdCount, $Molecule);
      next COMPOUND;
    }

    if ($SetupOutputFiles) {
      $SetupOutputFiles = 0;
      SetupFingerprintsLabelValueIDs($TopologicalPharmacophoreAtomTripletsFingerprints);
      ($NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO) = SetupAndOpenOutputFiles($FileIndex);
    }

    WriteDataToOutputFiles($FileIndex, $CmpdCount, $Molecule, $TopologicalPharmacophoreAtomTripletsFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO);
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
  my($TopologicalPharmacophoreAtomTripletsFingerprints) = @_;

  if ($OptionsInfo{AtomTripletsSetSizeToUse} =~ /^ArbitrarySize$/i ||
      $OptionsInfo{FingerprintsLabelMode} !~ /^FingerprintsLabelWithIDs$/i) {
    return;
  }
  $OptionsInfo{FingerprintsLabel} .= "; Value IDs: " . $TopologicalPharmacophoreAtomTripletsFingerprints->GetFingerprintsVector->GetValueIDsString();
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
  my($FileIndex, $CmpdCount, $Molecule, $TopologicalPharmacophoreAtomTripletsFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO) = @_;
  my($DataFieldLabelAndValuesRef);

  $DataFieldLabelAndValuesRef = undef;
  if ($NewFPTextFileIO || $NewFPFileIO) {
    $DataFieldLabelAndValuesRef = $Molecule->GetDataFieldLabelAndValues();
  }

  if ($NewFPSDFileIO) {
    my($CmpdString);

    $CmpdString = $Molecule->GetInputMoleculeString();
    $NewFPSDFileIO->WriteFingerprints($TopologicalPharmacophoreAtomTripletsFingerprints, $CmpdString);
  }

  if ($NewFPTextFileIO) {
    my($ColValuesRef);

    $ColValuesRef = SetupFPTextFileCoulmnValues($FileIndex, $CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPTextFileIO->WriteFingerprints($TopologicalPharmacophoreAtomTripletsFingerprints, $ColValuesRef);
  }

  if ($NewFPFileIO) {
    my($CompoundID);

    $CompoundID = SetupCmpdIDForOutputFiles($CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPFileIO->WriteFingerprints($TopologicalPharmacophoreAtomTripletsFingerprints, $CompoundID);
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
  my($TopologicalPharmacophoreAtomTripletsFingerprints);

  if ($OptionsInfo{KeepLargestComponent}) {
    $Molecule->KeepLargestComponent();
  }
  if (!$Molecule->DetectRings()) {
    return undef;
  }
  $Molecule->SetAromaticityModel($OptionsInfo{AromaticityModel});
  $Molecule->DetectAromaticity();

  $TopologicalPharmacophoreAtomTripletsFingerprints = new Fingerprints::TopologicalPharmacophoreAtomTripletsFingerprints('Molecule' => $Molecule, 'AtomTripletsSetSizeToUse' => $OptionsInfo{AtomTripletsSetSizeToUse}, 'MinDistance' => $OptionsInfo{MinDistance},  'MaxDistance' => $OptionsInfo{MaxDistance}, 'DistanceBinSize' => $OptionsInfo{DistanceBinSize}, 'UseTriangleInequality' => $OptionsInfo{UseTriangleInequality}, 'AtomTypesToUse' => \@{$OptionsInfo{AtomTypesToUse}});

  # Generate fingerprints...
  $TopologicalPharmacophoreAtomTripletsFingerprints->GenerateFingerprints();

  # Make sure fingerprints generation is successful...
  if (!$TopologicalPharmacophoreAtomTripletsFingerprints->IsFingerprintsGenerationSuccessful()) {
    return undef;
  }

  return $TopologicalPharmacophoreAtomTripletsFingerprints;
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
      $OutFileRoot = "${FileName}TopologicalPharmacophoreAtomTripletsFP";
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

  $OptionsInfo{AromaticityModel} = $Options{aromaticitymodel};

  $OptionsInfo{AtomTripletsSetSizeToUse} = $Options{atomtripletssetsizetouse};

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
  $OptionsInfo{FingerprintsLabel} = $Options{fingerprintslabel} ? $Options{fingerprintslabel} : 'TopologicalPharmacophoreAtomTripletsFingerprints';

  $OptionsInfo{KeepLargestComponent} = ($Options{keeplargestcomponent} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{DistanceBinSize} = $Options{distancebinsize};

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

  # Setup default vector string format...
  my($VectorStringFormat);
  $VectorStringFormat = '';

  if ($Options{vectorstringformat}) {
    $VectorStringFormat = $Options{vectorstringformat};

    if ($Options{atomtripletssetsizetouse} =~ /^ArbitrarySize$/i && $VectorStringFormat =~ /^ValuesString$/i) {
      die "Error: The value specified, $Options{vectorstringformat}, for option \"-v, --VectorStringFormat\" is not valid for $Options{atomtripletssetsizetouse} value of \"--AtomTripletsSetSizeToUse\" option. Allowed values: IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString or ValuesAndIDsPairsString\n";
    }
  }
  else {
    $VectorStringFormat = ($Options{atomtripletssetsizetouse} =~ /^FixedSize$/) ? "ValuesString" : "IDsAndValuesString";
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
      die "Error: Atom type specified, $AtomType, using \"-a, --AtomTypesToUse\" option is not valid...\n ";
    }
    push @{$OptionsInfo{AtomTypesToUse}}, $AtomType;
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{aromaticitymodel} = 'MayaChemToolsAromaticityModel';

  $Options{atomtripletssetsizetouse} = 'ArbitrarySize';

  $Options{atomtypestouse} = 'HBD,HBA,PI,NI,H,Ar';

  $Options{compoundidmode} = 'LabelPrefix';
  $Options{compoundidlabel} = 'CompoundID';
  $Options{datafieldsmode} = 'CompoundID';

  $Options{filter} = 'Yes';

  $Options{fingerprintslabelmode} = 'FingerprintsLabelOnly';

  $Options{keeplargestcomponent} = 'Yes';

  $Options{mindistance} = 1;
  $Options{maxdistance} = 10;

  $Options{distancebinsize} = 2;

  $Options{usetriangleinequality} = 'Yes';

  $Options{output} = 'text';
  $Options{outdelim} = 'comma';
  $Options{quote} = 'yes';

  $Options{vectorstringformat} = '';

  if (!GetOptions(\%Options, "aromaticitymodel=s", "atomtripletssetsizetouse=s", "atomtypestouse|a=s", "compoundid=s", "compoundidlabel=s", "compoundidmode=s", "datafields=s", "datafieldsmode|d=s", "distancebinsize=s", "filter|f=s", "fingerprintslabelmode=s", "fingerprintslabel=s", "help|h", "keeplargestcomponent|k=s",  "mindistance=s", "maxdistance=s", "outdelim=s", "output=s", "overwrite|o", "quote|q=s", "root|r=s", "usetriangleinequality|u=s", "vectorstringformat|v=s", "workingdir|w=s")) {
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
  if ($Options{atomtripletssetsizetouse} !~ /^(ArbitrarySize|FixedSize)$/i) {
    die "Error: The value specified, $Options{atomtripletssetsizetouse}, for option \"--AtomTripletsSetSizeToUse\" is not valid. Allowed values: ArbitrarySize or FixedSize\n";
  }
  if ($Options{compoundidmode} !~ /^(DataField|MolName|LabelPrefix|MolNameOrLabelPrefix)$/i) {
    die "Error: The value specified, $Options{compoundidmode}, for option \"--CompoundIDMode\" is not valid. Allowed values: DataField, MolName, LabelPrefix or MolNameOrLabelPrefix\n";
  }
  if ($Options{datafieldsmode} !~ /^(All|Common|Specify|CompoundID)$/i) {
    die "Error: The value specified, $Options{datafieldsmode}, for option \"-d, --DataFieldsMode\" is not valid. Allowed values: All, Common, Specify or CompoundID\n";
  }
  if (!IsPositiveInteger($Options{distancebinsize})) {
    die "Error: The value specified, $Options{distancebinsize}, for option \"--DistanceBinSize\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{filter} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{filter}, for option \"-f, --Filter\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{fingerprintslabelmode} !~ /^(FingerprintsLabelOnly|FingerprintsLabelWithIDs)$/i) {
    die "Error: The value specified, $Options{fingerprintslabelmode}, for option \"--FingerprintsLabelMode\" is not valid. Allowed values: FingerprintsLabelOnly or FingerprintsLabelWithIDs\n";
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
  if ($Options{vectorstringformat} && $Options{vectorstringformat} !~ /^(ValuesString|IDsAndValuesString|IDsAndValuesPairsString|ValuesAndIDsString|ValuesAndIDsPairsString)$/i) {
    die "Error: The value specified, $Options{vectorstringformat}, for option \"-v, --VectorStringFormat\" is not valid. Allowed values: ValuesString, IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString or ValuesAndIDsPairsString\n";
  }
}

__END__

=head1 NAME

TopologicalPharmacophoreAtomTripletsFingerprints.pl - Generate topological pharmacophore atom triplets fingerprints for SD files

=head1 SYNOPSIS

TopologicalPharmacophoreAtomTripletsFingerprints.pl SDFile(s)...

TopologicalPharmacophoreAtomTripletsFingerprints.pl [B<--AromaticityModel> I<AromaticityModelType>]
[B<--AtomTripletsSetSizeToUse> I<ArbitrarySize | FixedSize>]
[B<-a, --AtomTypesToUse> I<"AtomType1, AtomType2...">]
[B<--AtomTypesWeight> I<"AtomType1, Weight1, AtomType2, Weight2...">]
[B<--CompoundID> I<DataFieldName or LabelPrefixString>] [B<--CompoundIDLabel> I<text>]
[B<--CompoundIDMode>] [B<--DataFields> I<"FieldLabel1, FieldLabel2,...">]
[B<-d, --DataFieldsMode> I<All | Common | Specify | CompoundID>] [B<--DistanceBinSize> I<number>] [B<-f, --Filter> I<Yes | No>]
[B<--FingerprintsLabelMode> I<FingerprintsLabelOnly | FingerprintsLabelWithIDs>] [B<--FingerprintsLabel> I<text>]
[B<-h, --help>] [B<-k, --KeepLargestComponent> I<Yes | No>] [B<--MinDistance> I<number>] [B<--MaxDistance> I<number>]
[B<--OutDelim> I<comma | tab | semicolon>] [B<--output> I<SD | FP | text | all>] [B<-o, --overwrite>]
[B<-q, --quote> I<Yes | No>] [B<-r, --root> I<RootName>] [B<-u, --UseTriangleInequality> I<Yes | No>]
[B<-v, --VectorStringFormat> I<ValuesString, IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>]
[B<-w, --WorkingDir> dirname] SDFile(s)...

=head1 DESCRIPTION

Generate topological pharmacophore atom triplets fingerprints [ Ref 66, Ref 68-71 ] for
I<SDFile(s)> and create appropriate SD, FP or CSV/TSV text file(s) containing fingerprints vector
strings corresponding to molecular fingerprints.

Multiple SDFile names are separated by spaces. The valid file extensions are I<.sdf>
and I<.sd>. All other file names are ignored. All the SD files in a current directory
can be specified either by I<*.sdf> or the current directory name.

Based on the values specified for B<--AtomTypesToUse>, pharmacophore atom types are
assigned to all non-hydrogen atoms in a molecule and a distance matrix is generated.
Using B<--MinDistance>, B<--MaxDistance>, and B<--DistanceBinSize> values, a
binned distance matrix is generated with lower bound on the distance bin as the distance
in distance matrix; the lower bound on the distance bin is also used as the distance between
atom pairs for generation of atom triplet identifiers.

A pharmacophore atom triplets basis set is generated for all unique atom triplets constituting
atom pairs binned distances between B<--MinDistance> and B<--MaxDistance>. The value
of B<--UseTriangleInequality> determines whether the triangle inequality test is applied during
generation of atom triplets basis set. The lower distance bound, along with specified pharmacophore
types, is used during generation of atom triplet IDs.

    Let:

    P = Valid pharmacophore atom type

    Px = Pharmacophore atom x
    Py = Pharmacophore atom y
    Pz = Pharmacophore atom z

    Dmin = Minimum distance corresponding to number of bonds between two atoms
    Dmax = Maximum distance corresponding to number of bonds between two atoms
    D = Distance corresponding to number of bonds between two atom

    Bsize  = Distance bin size
    Nbins = Number of distance bins

    Dxy = Distance or lower bound of binned distance between Px and Py
    Dxz = Distance or lower bound of binned distance between Px and Pz
    Dyz = Distance or lower bound of binned distance between Py and Pz

    Then:

    PxDyz-PyDxz-PzDxy = Pharmacophore atom triplet IDs for atom types Px,
                        Py, and Pz

    For example: H1-H1-H1, H2-HBA-H2 and so on

    For default values of Dmin = 1 , Dmax = 10 and Bsize = 2:

    the number of distance bins, Nbins = 5, are:

    [1, 2] [3, 4] [5, 6] [7, 8] [9 10]

    and atom triplet basis set size is 2692.

    Atom triplet basis set size for various values of Dmin, Dmax and Bsize in
    conjunction with usage of triangle inequality is:

    Dmin    Dmax   Bsize   UseTriangleInequality   TripletBasisSetSize
    1       10     2       No                      4960
    1       10     2       Yes                     2692 [ Default ]
    2       12     2       No                      8436
    2       12     2       Yes                     4494


Using binned distance matrix and pharmacohore atom types, occurrence of unique pharmacohore
atom triplets is counted.

The final pharmacophore atom triples count along with atom pair identifiers involving all non-hydrogen
atoms constitute pharmacophore topological atom triplets fingerprints of the molecule.

For I<ArbitrarySize> value of B<--AtomTripletsSetSizeToUse> option, the fingerprint vector correspond to
only those topological pharmacophore atom triplets which are present and have non-zero count. However,
for I<FixedSize> value of B<--AtomTripletsSetSizeToUse> option, the fingerprint vector contains all possible
valid topological pharmacophore atom triplets with both zero and non-zero count values.

Example of I<SD> file containing topological pharmacophore atom triplets fingerprints string data:

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

    >  <TopologicalPharmacophoreAtomTripletsFingerprints>
    FingerprintsVector;TopologicalPharmacophoreAtomTriplets:ArbitrarySize:
    MinDistance1:MaxDistance10;696;NumericalValues;IDsAndValuesString;Ar1-
    Ar1-Ar1 Ar1-Ar1-H1 Ar1-Ar1-HBA1 Ar1-Ar1-HBD1 Ar1-H1-H1 Ar1-H1-HBA1 Ar1
    -H1-HBD1 Ar1-HBA1-HBD1 H1-H1-H1 H1-H1-HBA1 H1-H1-HBD1 H1-HBA1-HBA1 H1-
    HBA1-HBD1 H1-HBA1-NI1 H1-HBD1-NI1 HBA1-HBA1-NI1 HBA1-HBD1-NI1 Ar1-...;
    46 106 8 3 83 11 4 1 21 5 3 1 2 2 1 1 1 100 101 18 11 145 132 26 14 23
    28 3 3 5 4 61 45 10 4 16 20 7 5 1 3 4 5 3 1 1 1 1 5 4 2 1 2 2 2 1 1 1
    119 123 24 15 185 202 41 25 22 17 3 5 85 95 18 11 23 17 3 1 1 6 4 ...

    $$$$
    ... ...
    ... ...

Example of I<FP> file containing topological pharmacophore atom triplets fingerprints string data:

    #
    # Package = MayaChemTools 7.4
    # Release Date = Oct 21, 2010
    #
    # TimeStamp = Fri Mar 11 15:38:58 2011
    #
    # FingerprintsStringType = FingerprintsVector
    #
    # Description = TopologicalPharmacophoreAtomTriplets:ArbitrarySize:M...
    # VectorStringFormat = IDsAndValuesString
    # VectorValuesType = NumericalValues
    #
    Cmpd1 696;Ar1-Ar1-Ar1 Ar1-Ar1-H1 Ar1-Ar1-HBA1 Ar1-Ar1-HBD1...;;46 106...
    Cmpd2 251;H1-H1-H1 H1-H1-HBA1 H1-H1-HBD1 H1-H1-NI1...;4 1 3 1 1 2 2...
    ... ...
    ... ..

Example of CSV I<Text> file containing topological pharmacophore atom triplets fingerprints string data:

    "CompoundID","TopologicalPharmacophoreAtomTripletsFingerprints"
    "Cmpd1","FingerprintsVector;TopologicalPharmacophoreAtomTriplets:Arbitr
    arySize:MinDistance1:MaxDistance10;696;NumericalValues;IDsAndValuesStri
    ng;Ar1-Ar1-Ar1 Ar1-Ar1-H1 Ar1-Ar1-HBA1 Ar1-Ar1-HBD1 Ar1-H1-H1 Ar1-H1-HB
    A1 Ar1-H1-HBD1 Ar1-HBA1-HBD1 H1-H1-H1 H1-H1-HBA1 H1-H1-HBD1 H1-HBA1-HBA
    1 H1-HBA1-HBD1 H1-HBA1-NI1 H1-HBD1-NI1 HBA1-HBA1-NI1 HBA1-HBD1-NI1 A...;
    46 106 8 3 83 11 4 1 21 5 3 1 2 2 1 1 1 100 101 18 11 145 132 26 14 23
    28 3 3 5 4 61 45 10 4 16 20 7 5 1 3 4 5 3 1 1 1 1 5 4 2 1 2 2 2 1 1 1
    119 123 24 15 185 202 41 25 22 17 3 5 85 95 18 11 23 17 3 1 1 6 4 ...
    ... ...
    ... ...

The current release of MayaChemTools generates the following types of topological pharmacophore
atom triplets fingerprints vector strings:

    FingerprintsVector;TopologicalPharmacophoreAtomTriplets:ArbitrarySize:
    MinDistance1:MaxDistance10;696;NumericalValues;IDsAndValuesString;Ar1-
    Ar1-Ar1 Ar1-Ar1-H1 Ar1-Ar1-HBA1 Ar1-Ar1-HBD1 Ar1-H1-H1 Ar1-H1-HBA1 Ar1
    -H1-HBD1 Ar1-HBA1-HBD1 H1-H1-H1 H1-H1-HBA1 H1-H1-HBD1 H1-HBA1-HBA1 H1-
    HBA1-HBD1 H1-HBA1-NI1 H1-HBD1-NI1 HBA1-HBA1-NI1 HBA1-HBD1-NI1 Ar1-...;
    46 106 8 3 83 11 4 1 21 5 3 1 2 2 1 1 1 100 101 18 11 145 132 26 14 23
    28 3 3 5 4 61 45 10 4 16 20 7 5 1 3 4 5 3 1 1 1 1 5 4 2 1 2 2 2 1 1 1
    119 123 24 15 185 202 41 25 22 17 3 5 85 95 18 11 23 17 3 1 1 6 4 ...

    FingerprintsVector;TopologicalPharmacophoreAtomTriplets:FixedSize:MinD
    istance1:MaxDistance10;2692;OrderedNumericalValues;ValuesString;46 106
    8 3 0 0 83 11 4 0 0 0 1 0 0 0 0 0 0 0 0 21 5 3 0 0 1 2 2 0 0 1 0 0 0
    0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 101 18 11 0 0 145 132 26
    14 0 0 23 28 3 3 0 0 5 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 61 45 10 4 0
    0 16 20 7 5 1 0 3 4 5 3 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 5 ...

    FingerprintsVector;TopologicalPharmacophoreAtomTriplets:FixedSize:MinD
    istance1:MaxDistance10;2692;OrderedNumericalValues;IDsAndValuesString;
    Ar1-Ar1-Ar1 Ar1-Ar1-H1 Ar1-Ar1-HBA1 Ar1-Ar1-HBD1 Ar1-Ar1-NI1 Ar1-Ar1-P
    I1 Ar1-H1-H1 Ar1-H1-HBA1 Ar1-H1-HBD1 Ar1-H1-NI1 Ar1-H1-PI1 Ar1-HBA1-HB
    A1 Ar1-HBA1-HBD1 Ar1-HBA1-NI1 Ar1-HBA1-PI1 Ar1-HBD1-HBD1 Ar1-HBD1-...;
    46 106 8 3 0 0 83 11 4 0 0 0 1 0 0 0 0 0 0 0 0 21 5 3 0 0 1 2 2 0 0 1
    0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 101 18 11 0 0 145
    132 26 14 0 0 23 28 3 3 0 0 5 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 61 ...

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

=item B<--AtomTripletsSetSizeToUse> I<ArbitrarySize | FixedSize>

Atom triplets set size to use during generation of topological pharmacophore atom triplets
fingerprints.

Possible values: I<ArbitrarySize | FixedSize>; Default value: I<ArbitrarySize>.

For I<ArbitrarySize> value of B<--AtomTripletsSetSizeToUse> option, the fingerprint vector
correspond to only those topological pharmacophore atom triplets which are present and
have non-zero count. However, for I<FixedSize> value of B<--AtomTripletsSetSizeToUse>
option, the fingerprint vector contains all possible valid topological pharmacophore atom
triplets with both zero and non-zero count values.

=item B<-a, --AtomTypesToUse> I<"AtomType1,AtomType2,...">

Pharmacophore atom types to use during generation of topological phramacophore
atom triplets. It's a list of comma separated valid pharmacophore atom types.

Possible values for pharmacophore atom types are: I<Ar, CA, H, HBA, HBD, Hal, NI, PI, RA>.
Default value [ Ref 71 ] : I<HBD,HBA,PI,NI,H,Ar>.

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

=item B<--DistanceBinSize> I<number>

Distance bin size used to bin distances between atom pairs in atom triplets. Default value: I<2>.
Valid values: positive integers.

For default B<--MinDistance> and B<--MaxDistance> values of 1 and 10 with  B<--DistanceBinSize>
of 2 [ Ref 70 ], the following 5 distance bins are generated:

    [1, 2] [3, 4] [5, 6] [7, 8] [9 10]

The lower distance bound on the distance bin is uses to bin the distance between atom pairs in
atom triplets. So in the previous example, atom pairs with distances 1 and 2 fall in first distance
bin, atom pairs with distances 3 and 4  fall in second distance bin and so on.

In order to distribute distance bins of equal size, the last bin is allowed to go past B<--MaxDistance>
by up to distance bin size. For example, B<--MinDistance> and B<--MaxDistance> values of 2 and 10
with B<--DistanceBinSize> of 2 generates the following 6 distance bins:

    [2, 3] [4, 5] [6, 7] [8, 9] [10 11]

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
of B<--AtomTripletsSetSizeToUse> option and topological atom triplets IDs not appended to the label.

=item B<--FingerprintsLabel> I<text>

SD data label or text file column label to use for fingerprints string in output SD or
CSV/TSV text file(s) specified by B<--output>. Default value: I<TopologicalPharmacophoreAtomTripletsFingerprints>.

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

Minimum bond distance between atom pairs corresponding to atom triplets for generating
topological pharmacophore atom triplets. Default value: I<1>. Valid values: positive integers and
less than B<--MaxDistance>.

=item B<--MaxDistance> I<number>

Maximum bond distance between atom pairs corresponding to atom triplets for generating
topological pharmacophore atom triplets. Default value: I<10>. Valid values: positive integers and
greater than B<--MinDistance>.

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
<SDFileName><TopologicalPharmacophoreAtomTripletsFP>.<Ext>. The file type determines <Ext> value.
The sdf, fpf, csv, and tsv <Ext> values are used for SD, FP, comma/semicolon, and tab
delimited text files, respectively.This option is ignored for multiple input files.

=item B<-u, --UseTriangleInequality> I<Yes | No>

Specify whether to imply triangle distance inequality test to distances between atom pairs in
atom triplets during generation of atom triplets basis set generation. Possible values:
I<Yes or No>. Default value: I<Yes>.

Triangle distance inequality test implies that distance or binned distance between any two atom
pairs in an atom triplet must be less than the sum of distances or binned distances between other
two atoms pairs and greater than the difference of their distances.

    For atom triplet PxDyz-PyDxz-PzDxy to satisfy triangle inequality:

    Dyz > |Dxz - Dxy| and Dyz < Dxz + Dxy
    Dxz > |Dyz - Dxy| and Dyz < Dyz + Dxy
    Dxy > |Dyz - Dxz| and Dxy < Dyz + Dxz

=item B<-v, --VectorStringFormat> I<ValuesString, IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>

Format of fingerprints vector string data in output SD, FP or CSV/TSV text file(s) specified by
B<--output> option. Possible values: I<ValuesString, IDsAndValuesString | IDsAndValuesPairsString |
ValuesAndIDsString | ValuesAndIDsPairsString>. Defaultvalue: I<ValuesString>.

Default value during I<FixedSize> value of B<--AtomTripletsSetSizeToUse> option: I<ValuesString>. Default
value during I<ArbitrarySize> value of B<--AtomTripletsSetSizeToUse> option: I<IDsAndValuesString>.

I<ValuesString> option value is not allowed for I<ArbitrarySize> value of B<--AtomTripletsSetSizeToUse>
option.

Examples:

    FingerprintsVector;TopologicalPharmacophoreAtomTriplets:ArbitrarySize:
    MinDistance1:MaxDistance10;696;NumericalValues;IDsAndValuesString;Ar1-
    Ar1-Ar1 Ar1-Ar1-H1 Ar1-Ar1-HBA1 Ar1-Ar1-HBD1 Ar1-H1-H1 Ar1-H1-HBA1 Ar1
    -H1-HBD1 Ar1-HBA1-HBD1 H1-H1-H1 H1-H1-HBA1 H1-H1-HBD1 H1-HBA1-HBA1 H1-
    HBA1-HBD1 H1-HBA1-NI1 H1-HBD1-NI1 HBA1-HBA1-NI1 HBA1-HBD1-NI1 Ar1-...;
    46 106 8 3 83 11 4 1 21 5 3 1 2 2 1 1 1 100 101 18 11 145 132 26 14 23
    28 3 3 5 4 61 45 10 4 16 20 7 5 1 3 4 5 3 1 1 1 1 5 4 2 1 2 2 2 1 1 1
    119 123 24 15 185 202 41 25 22 17 3 5 85 95 18 11 23 17 3 1 1 6 4 ...

    FingerprintsVector;TopologicalPharmacophoreAtomTriplets:FixedSize:MinD
    istance1:MaxDistance10;2692;OrderedNumericalValues;ValuesString;46 106
    8 3 0 0 83 11 4 0 0 0 1 0 0 0 0 0 0 0 0 21 5 3 0 0 1 2 2 0 0 1 0 0 0
    0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 101 18 11 0 0 145 132 26
    14 0 0 23 28 3 3 0 0 5 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 61 45 10 4 0
    0 16 20 7 5 1 0 3 4 5 3 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 5 ...

    FingerprintsVector;TopologicalPharmacophoreAtomTriplets:FixedSize:MinD
    istance1:MaxDistance10;2692;OrderedNumericalValues;ValuesAndIDsPairsSt
    ring;46 Ar1-Ar1-Ar1 106 Ar1-Ar1-H1 8 Ar1-Ar1-HBA1 3 Ar1-Ar1-HBD1 0 Ar1
    -Ar1-NI1 0 Ar1-Ar1-PI1 83 Ar1-H1-H1 11 Ar1-H1-HBA1 4 Ar1-H1-HBD1 0 Ar1
    -H1-NI1 0 Ar1-H1-PI1 0 Ar1-HBA1-HBA1 1 Ar1-HBA1-HBD1 0 Ar1-HBA1-NI1 0
    Ar1-HBA1-PI1 0 Ar1-HBD1-HBD1 0 Ar1-HBD1-NI1 0 Ar1-HBD1-PI1 0 Ar1-NI...

=item B<-w, --WorkingDir> I<DirName>

Location of working directory. Default value: current directory.

=back

=head1 EXAMPLES

To generate topological pharmacophore atom triplets fingerprints  of arbitrary size corresponding to 5
distance bins spanning distances from 1 through 10 using default atoms with distances satisfying triangle
inequality and create a SampleTPATFP.csv file containing sequential compound IDs along
with fingerprints vector strings data in ValuesString format, type:

    % TopologicalPharmacophoreAtomTripletsFingerprints.pl -r SampleTPATFP
      -o Sample.sdf

To generate topological pharmacophore atom triplets fingerprints  of fixed size corresponding to 5
distance bins spanning distances from 1 through 10 using default atoms with distances satisfying triangle
inequality and create a SampleTPATFP.csv file containing sequential compound IDs along
with fingerprints vector strings data in ValuesString format, type:

    % TopologicalPharmacophoreAtomTripletsFingerprints.pl
      --AtomTripletsSetSizeToUse FixedSize -r SampleTPATFP -o Sample.sdf

To generate topological pharmacophore atom triplets fingerprints  of arbitrary size corresponding to 5
distance bins spanning distances from 1 through 10 using default atoms with distances satisfying triangle
inequality and create SampleTPATFP.sdf, SampleTPATFP.fpf and SampleTPATFP.csv files with CSV file containing
sequential compound IDs along with fingerprints vector strings data in ValuesString format, type:

    % TopologicalPharmacophoreAtomTripletsFingerprints.pl --output all
      -r SampleTPATFP -o Sample.sdf

To generate topological pharmacophore atom triplets fingerprints  of arbitrary size corresponding to 5
distance bins spanning distances from 1 through 10 using default atoms with distances satisfying triangle
inequality and create a SampleTPATFP.csv file containing sequential compound IDs along
with fingerprints vector strings data in ValuesString format and atom triplets IDs in the
fingerprint data column label starting with Fingerprints, type:

    % TopologicalPharmacophoreAtomTripletsFingerprints.pl
      --FingerprintsLabelMode FingerprintsLabelWithIDs --FingerprintsLabel
      Fingerprints -r SampleTPATFP -o Sample.sdf

To generate topological pharmacophore atom triplets fingerprints  of arbitrary size corresponding to 5
distance bins spanning distances from 1 through 10 using default atoms with distances not satisfying triangle
inequality and create a SampleTPATFP.csv file containing sequential compound IDs along
with fingerprints vector strings data in ValuesString format, type:

    % TopologicalPharmacophoreAtomTripletsFingerprints.pl
      --UseTriangleInequality No -r SampleTPATFP -o Sample.sdf

To generate topological pharmacophore atom triplets fingerprints  of arbitrary size corresponding to 6
distance bins spanning distances from 1 through 12 using default atoms with distances satisfying triangle
inequality and create a SampleTPATFP.csv file containing sequential compound IDs along
with fingerprints vector strings data in ValuesString format, type:

    % TopologicalPharmacophoreAtomTripletsFingerprints.pl
      --UseTriangleInequality Yes --MinDistance 1 --MaxDistance 12
      --DistanceBinSIze 2 -r SampleTPATFP -o Sample.sdf

To generate topological pharmacophore atom triplets fingerprints  of arbitrary size corresponding to 6
distance bins spanning distances from 1 through 12 using "HBD,HBA,PI, NI, H, Ar" atoms with distances
satisfying triangle inequality and create a SampleTPATFP.csv file containing sequential compound
IDs along with fingerprints vector strings data in ValuesString format, type:

    % TopologicalPharmacophoreAtomTripletsFingerprints.pl
      --AtomTypesToUse "HBD,HBA,PI,NI,H,Ar" --UseTriangleInequality Yes
      --MinDistance 1 --MaxDistance 12 --DistanceBinSIze 2
      --VectorStringFormat ValuesString -r SampleTPATFP -o Sample.sdf

To generate topological pharmacophore atom triplets fingerprints  of arbitrary size corresponding to 5
distance bins spanning distances from 1 through 10 using default atoms with distances satisfying triangle
inequality and create a SampleTPATFP.csv file containing sequential compound IDs from
molecule name line along with fingerprints vector strings data in ValuesString format, type:

    % TopologicalPharmacophoreAtomTripletsFingerprints.pl --DataFieldsMode
      CompoundID -CompoundIDMode MolName  -r SampleTPATFP -o Sample.sdf

To generate topological pharmacophore atom triplets fingerprints  of arbitrary size corresponding to 5
distance bins spanning distances from 1 through 10 using default atoms with distances satisfying triangle
inequality and create a SampleTPATFP.csv file containing sequential compound IDs using
specified data field along with fingerprints vector strings data in ValuesString format, type:

    % TopologicalPharmacophoreAtomTripletsFingerprints.pl --DataFieldsMode
      CompoundID -CompoundIDMode DataField --CompoundID Mol_ID
      -r SampleTPATFP -o Sample.sdf

To generate topological pharmacophore atom triplets fingerprints  of arbitrary size corresponding to 5
distance bins spanning distances from 1 through 10 using default atoms with distances satisfying triangle
inequality and create a SampleTPATFP.csv file containing sequential compound IDs using
combination of molecule name line and an explicit compound prefix along with fingerprints vector
strings data, type:

    % TopologicalPharmacophoreAtomTripletsFingerprints.pl --DataFieldsMode
      CompoundID -CompoundIDMode MolnameOrLabelPrefix
      --CompoundID Cmpd --CompoundIDLabel MolID -r SampleSampleTPATFP
      -o Sample.sdf

To generate topological pharmacophore atom triplets fingerprints  of arbitrary size corresponding to 5
distance bins spanning distances from 1 through 10 using default atoms with distances satisfying triangle
inequality and create a SampleTPATFP.csv file containing specific data fields columns along
with fingerprints vector strings data, type:

    % TopologicalPharmacophoreAtomTripletsFingerprints.pl --DataFieldsMode
      Specify --DataFields Mol_ID -r SampleTPATFP -o Sample.sdf

To generate topological pharmacophore atom triplets fingerprints  of arbitrary size corresponding to 5
distance bins spanning distances from 1 through 10 using default atoms with distances satisfying triangle
inequality and create a SampleTPATFP.csv file containing common data fields columns along
with fingerprints vector strings data, type:

    % TopologicalPharmacophoreAtomTripletsFingerprints.pl --DataFieldsMode
      Common -r SampleTPATFP -o Sample.sdf

To generate topological pharmacophore atom triplets fingerprints  of arbitrary size corresponding to 5
distance bins spanning distances from 1 through 10 using default atoms with distances satisfying triangle
inequality and create SampleTPATFP.sdf, SampleTPATFP.fpf and SampleTPATFP.csv files containing all
data fields columns in CSV file along with fingerprints data, type:

    % TopologicalPharmacophoreAtomTripletsFingerprints.pl --DataFieldsMode
      All  --output all -r SampleTPATFP -o Sample.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoFingerprintsFiles.pl, SimilarityMatricesFingerprints.pl, AtomNeighborhoodsFingerprints.pl,
ExtendedConnectivityFingerprints.pl, MACCSKeysFingerprints.pl, PathLengthFingerprints.pl,
TopologicalAtomPairsFingerprints.pl, TopologicalAtomTorsionsFingerprints.pl,
TopologicalPharmacophoreAtomPairsFingerprints.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
