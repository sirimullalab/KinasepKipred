#!/usr/bin/perl -w
#
# File: MACCSKeysFingerprints.pl
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
use Fingerprints::MACCSKeys;

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
    GenerateMACCSKeysFingerprints($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Generate fingerprints for a SD file...
#
sub GenerateMACCSKeysFingerprints {
  my($FileIndex) = @_;
  my($CmpdCount, $IgnoredCmpdCount, $SDFile, $MoleculeFileIO, $Molecule, $MACCSKeysFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO);

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

    $MACCSKeysFingerprints = GenerateMoleculeFingerprints($Molecule);
    if (!$MACCSKeysFingerprints) {
      $IgnoredCmpdCount++;
      ProcessIgnoredCompound('FingerprintsGenerationFailed', $CmpdCount, $Molecule);
      next COMPOUND;
    }

    WriteDataToOutputFiles($FileIndex, $CmpdCount, $Molecule, $MACCSKeysFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO);
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
  if ($OptionsInfo{Mode} =~ /^MACCSKeyBits$/i) {
    %FingerprintsFileIOParams = ('Mode' => 'Write', 'Overwrite' => $OptionsInfo{OverwriteFiles}, 'FingerprintsStringMode' => 'FingerprintsBitVectorString', 'BitStringFormat' => $OptionsInfo{BitStringFormat}, 'BitsOrder' => $OptionsInfo{BitsOrder});
  }
  elsif ($OptionsInfo{Mode} =~ /^MACCSKeyCount$/i) {
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
  my($FileIndex, $CmpdCount, $Molecule, $MACCSKeysFingerprints, $NewFPSDFileIO, $NewFPTextFileIO, $NewFPFileIO) = @_;
  my($DataFieldLabelAndValuesRef);

  $DataFieldLabelAndValuesRef = undef;
  if ($NewFPTextFileIO || $NewFPFileIO) {
    $DataFieldLabelAndValuesRef = $Molecule->GetDataFieldLabelAndValues();
  }

  if ($NewFPSDFileIO) {
    my($CmpdString);

    $CmpdString = $Molecule->GetInputMoleculeString();
    $NewFPSDFileIO->WriteFingerprints($MACCSKeysFingerprints, $CmpdString);
  }

  if ($NewFPTextFileIO) {
    my($ColValuesRef);

    $ColValuesRef = SetupFPTextFileCoulmnValues($FileIndex, $CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPTextFileIO->WriteFingerprints($MACCSKeysFingerprints, $ColValuesRef);
  }

  if ($NewFPFileIO) {
    my($CompoundID);

    $CompoundID = SetupCmpdIDForOutputFiles($CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    $NewFPFileIO->WriteFingerprints($MACCSKeysFingerprints, $CompoundID);
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
  my($MACCSKeysFingerprints);

  if ($OptionsInfo{KeepLargestComponent}) {
    $Molecule->KeepLargestComponent();
  }
  if (!$Molecule->DetectRings()) {
    return undef;
  }
  $Molecule->SetAromaticityModel($OptionsInfo{AromaticityModel});
  $Molecule->DetectAromaticity();

  $MACCSKeysFingerprints = undef;
  if ($OptionsInfo{Mode} =~ /^MACCSKeyBits$/i) {
    $MACCSKeysFingerprints = new Fingerprints::MACCSKeys('Molecule' => $Molecule, 'Type' => 'MACCSKeyBits', 'Size' => $OptionsInfo{Size});
  }
  elsif ($OptionsInfo{Mode} =~ /^MACCSKeyCount$/i) {
    $MACCSKeysFingerprints = new Fingerprints::MACCSKeys('Molecule' => $Molecule, 'Type' => 'MACCSKeyCount', 'Size' => $OptionsInfo{Size});
  }
  else {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: MACCSKeyBits or MACCSKeyCount\n";
  }
  $MACCSKeysFingerprints->GenerateMACCSKeys();

  return $MACCSKeysFingerprints;
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
      $OutFileRoot = "${FileName}MACCSKeysFP";
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

  $OptionsInfo{BitsOrder} = $Options{bitsorder};
  $OptionsInfo{BitStringFormat} = $Options{bitstringformat};

  $OptionsInfo{CompoundIDMode} = $Options{compoundidmode};
  $OptionsInfo{CompoundIDLabel} = $Options{compoundidlabel};
  $OptionsInfo{DataFieldsMode} = $Options{datafieldsmode};

  $OptionsInfo{Filter} = ($Options{filter} =~ /^Yes$/i) ? 1 : 0;

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

  $OptionsInfo{FingerprintsLabel} = $Options{fingerprintslabel} ? $Options{fingerprintslabel} : 'MACCSKeysFingerprints';

  $OptionsInfo{KeepLargestComponent} = ($Options{keeplargestcomponent} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{Output} = $Options{output};
  $OptionsInfo{SDOutput} = ($Options{output} =~ /^(SD|All)$/i) ? 1 : 0;
  $OptionsInfo{FPOutput} = ($Options{output} =~ /^(FP|All)$/i) ? 1 : 0;
  $OptionsInfo{TextOutput} = ($Options{output} =~ /^(Text|All)$/i) ? 1 : 0;

  $OptionsInfo{OutDelim} = $Options{outdelim};
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? 1 : 0;
  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : 0;

  $OptionsInfo{Size} = $Options{size};

  $OptionsInfo{VectorStringFormat} = $Options{vectorstringformat};
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{aromaticitymodel} = 'MayaChemToolsAromaticityModel';

  $Options{bitsorder} = 'Ascending';
  $Options{bitstringformat} = 'BinaryString';

  $Options{compoundidmode} = 'LabelPrefix';
  $Options{compoundidlabel} = 'CompoundID';
  $Options{datafieldsmode} = 'CompoundID';

  $Options{filter} = 'Yes';

  $Options{keeplargestcomponent} = 'Yes';

  $Options{mode} = 'MACCSKeyBits';

  $Options{output} = 'text';
  $Options{outdelim} = 'comma';
  $Options{quote} = 'yes';

  $Options{size} = 166;

  $Options{vectorstringformat} = 'ValuesString';

  if (!GetOptions(\%Options, "aromaticitymodel=s", "bitsorder=s", "bitstringformat|b=s", "compoundid=s", "compoundidlabel=s", "compoundidmode=s", "datafields=s", "datafieldsmode|d=s", "filter|f=s", "fingerprintslabel=s",  "help|h", "keeplargestcomponent|k=s", "mode|m=s", "outdelim=s", "output=s", "overwrite|o", "quote|q=s", "root|r=s", "size|s=i", "vectorstringformat|v=s", "workingdir|w=s")) {
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
  if ($Options{mode} !~ /^(MACCSKeyBits|MACCSKeyCount)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: MACCSKeyBits or MACCSKeyCount\n";
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
  if (!(IsPositiveInteger($Options{size}) && ($Options{size} == 166 || $Options{size} == 322))) {
    die "Error: The value specified, $Options{size}, for option \"-s, --size\" is not valid. Allowed values: 166 or 322 \n";
  }
  if ($Options{vectorstringformat} !~ /^(ValuesString|IDsAndValuesString|IDsAndValuesPairsString|ValuesAndIDsString|ValuesAndIDsPairsString)$/i) {
    die "Error: The value specified, $Options{vectorstringformat}, for option \"-v, --VectorStringFormat\" is not valid. Allowed values: ValuesString, IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString or ValuesAndIDsPairsString\n";
  }
}

__END__

=head1 NAME

MACCSKeysFingerprints.pl - Generate MACCS key fingerprints for SD files

=head1 SYNOPSIS

MACCSKeysFingerprints.pl SDFile(s)...

MACCSKeysFingerprints.pl [B<--AromaticityModel> I<AromaticityModelType>]
[B<--BitsOrder> I<Ascending | Descending>]
[B<-b, --BitStringFormat> I<BinaryString | HexadecimalString>]
[B<--CompoundID> I<DataFieldName or LabelPrefixString>] [B<--CompoundIDLabel> I<text>]
[B<--CompoundIDMode> I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>]
[B<--DataFields> I<"FieldLabel1,FieldLabel2,...">] [B<-d, --DataFieldsMode> I<All | Common | Specify | CompoundID>]
[B<-f, --Filter> I<Yes | No>] [B<--FingerprintsLabel> I<text>] [B<-h, --help>] [B<-k, --KeepLargestComponent> I<Yes | No>]
[B<-m, --mode> I<MACCSKeyBits | MACCSKeyCount>] [B<--OutDelim> I<comma | tab | semicolon>]
[B<--output> I<SD | FP | text | all>] [B<-o, --overwrite>]
[B<-q, --quote> I<Yes | No>] [B<-r, --root> I<RootName>] [B<-s, --size> I<number>]
[B<-v, --VectorStringFormat> I<IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>]
[B<-w, --WorkingDir> I<DirName>]

=head1 DESCRIPTION

Generate MACCS (Molecular ACCess System) keys fingerprints [ Ref 45-47 ] for I<SDFile(s)>
and create appropriate SD, FP or CSV/TSV text file(s) containing fingerprints bit-vector or
vector strings corresponding to molecular fingerprints.

Multiple SDFile names are separated by spaces. The valid file extensions are I<.sdf>
and I<.sd>. All other file names are ignored. All the SD files in a current directory
can be specified either by I<*.sdf> or the current directory name.

For each MACCS keys definition, atoms are processed to determine their membership to the key
and the appropriate molecular fingerprints strings are generated. An atom can belong to multiple
MACCS keys.

For I<MACCSKeyBits> value of B<-m, --mode> option, a fingerprint bit-vector string containing
zeros and ones is generated and for I<MACCSKeyCount> value, a fingerprint vector string
corresponding to number of MACCS keys [ Ref 45-47 ] is generated.

I<MACCSKeyBits | MACCSKeyCount> values for B<-m, --mode> option along with two possible
I<166 | 322>  values of B<-s, --size> supports generation of four different types of MACCS
keys fingerprint: I<MACCS166KeyBits, MACCS166KeyCount, MACCS322KeyBits, MACCS322KeyCount>.

Example of I<SD> file containing MAACS keys fingerprints string data:

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

    >  <MACCSKeysFingerprints>
    FingerprintsBitVector;MACCSKeyBits;166;BinaryString;Ascending;000000000
    00000000000000000000000000000000100100001001000000001001000000001110001
    00101010111100011011000100110110000011011110100110111111111111011111111
    11111111110111000

    $$$$
    ... ...
    ... ...

Example of I<FP> file containing MAACS keys fingerprints string data:

    #
    # Package = MayaChemTools 7.4
    # Release Date = Oct 21, 2010
    #
    # TimeStamp = Fri Mar 11 14:57:24 2011
    #
    # FingerprintsStringType = FingerprintsBitVector
    #
    # Description = MACCSKeyBits
    # Size = 166
    # BitStringFormat = BinaryString
    # BitsOrder = Ascending
    #
    Cmpd1 00000000000000000000000000000000000000000100100001001000000001...
    Cmpd2 00000000000000000000000010000000001000000010000000001000000000...
    ... ...
    ... ..

Example of CSV I<Text> file containing MAACS keys fingerprints string data:

    "CompoundID","MACCSKeysFingerprints"
    "Cmpd1","FingerprintsBitVector;MACCSKeyBits;166;BinaryString;Ascending;
    00000000000000000000000000000000000000000100100001001000000001001000000
    00111000100101010111100011011000100110110000011011110100110111111111111
    01111111111111111110111000"
    ... ...
    ... ...

The current release of MayaChemTools generates the following types of MACCS keys
fingerprints bit-vector and vector strings:

    FingerprintsBitVector;MACCSKeyBits;166;BinaryString;Ascending;00000000
    0000000000000000000000000000000001001000010010000000010010000000011100
    0100101010111100011011000100110110000011011110100110111111111111011111
    11111111111110111000

    FingerprintsBitVector;MACCSKeyBits;166;HexadecimalString;Ascending;000
    000000021210210e845f8d8c60b79dffbffffd1

    FingerprintsBitVector;MACCSKeyBits;322;BinaryString;Ascending;11101011
    1110011111100101111111000111101100110000000000000011100010000000000000
    0000000000000000000000000000000000000000000000101000000000000000000000
    0000000000000000000000000000000000000000000000000000000000000000000000
    0000000000000000000000000000000000000011000000000000000000000000000000
    0000000000000000000000000000000000000000

    FingerprintsBitVector;MACCSKeyBits;322;HexadecimalString;Ascending;7d7
    e7af3edc000c1100000000000000500000000000000000000000000000000300000000
    000000000

    FingerprintsVector;MACCSKeyCount;166;OrderedNumericalValues;ValuesStri
    ng;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 1 0 0 3 0 0 0 0 4 0 0 2 0 0 0 0 0 0 0 0 2 0 0 2 0 0 0 0
    0 0 0 0 1 1 8 0 0 0 1 0 0 1 0 1 0 1 0 3 1 3 1 0 0 0 1 2 0 11 1 0 0 0
    5 0 0 1 2 0 1 1 0 0 0 0 0 1 1 0 1 1 1 1 0 4 0 0 1 1 0 4 6 1 1 1 2 1 1
    3 5 2 2 0 5 3 5 1 1 2 5 1 2 1 2 4 8 3 5 5 2 2 0 3 5 4 1

    FingerprintsVector;MACCSKeyCount;322;OrderedNumericalValues;ValuesStri
    ng;14 8 2 0 2 0 4 4 2 1 4 0 0 2 5 10 5 2 1 0 0 2 0 5 13 3 28 5 5 3 0 0
    0 4 2 1 1 0 1 1 0 0 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 22 5 3 0 0 0 1 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 2 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...

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

=item B<--BitsOrder> I<Ascending | Descending>

Bits order to use during generation of fingerprints bit-vector string for I<MACCSKeyBits> value of
B<-m, --mode> option. Possible values: I<Ascending, Descending>. Default: I<Ascending>.

I<Ascending> bit order which corresponds to first bit in each byte as the lowest bit as
opposed to the highest bit.

Internally, bits are stored in I<Ascending> order using Perl vec function. Regardless
of machine order, big-endian or little-endian, vec function always considers first
string byte as the lowest byte and first bit within each byte as the lowest bit.

=item B<-b, --BitStringFormat> I<BinaryString | HexadecimalString>

Format of fingerprints bit-vector string data in output SD, FP or CSV/TSV text file(s) specified by
B<--output> used during I<MACCSKeyBits> value of B<-m, --mode> option. Possible
values: I<BinaryString, HexadecimalString>. Default value: I<BinaryString>.

I<BinaryString> corresponds to an ASCII string containing 1s and 0s. I<HexadecimalString>
contains bit values in ASCII hexadecimal format.

Examples:

    FingerprintsBitVector;MACCSKeyBits;166;BinaryString;Ascending;00000000
    0000000000000000000000000000000001001000010010000000010010000000011100
    0100101010111100011011000100110110000011011110100110111111111111011111
    11111111111110111000

    FingerprintsBitVector;MACCSKeyBits;166;HexadecimalString;Ascending;000
    000000021210210e845f8d8c60b79dffbffffd1

    FingerprintsBitVector;MACCSKeyBits;322;BinaryString;Ascending;11101011
    1110011111100101111111000111101100110000000000000011100010000000000000
    0000000000000000000000000000000000000000000000101000000000000000000000
    0000000000000000000000000000000000000000000000000000000000000000000000
    0000000000000000000000000000000000000011000000000000000000000000000000
    0000000000000000000000000000000000000000

    FingerprintsBitVector;MACCSKeyBits;322;HexadecimalString;Ascending;7d7
    e7af3edc000c1100000000000000500000000000000000000000000000000300000000
    000000000

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
CSV/TSV text file(s) specified by B<--output>. Default value: I<MACCSKeyFingerprints>.

=item B<-h, --help>

Print this help message.

=item B<-k, --KeepLargestComponent> I<Yes | No>

Generate fingerprints for only the largest component in molecule. Possible values:
I<Yes or No>. Default value: I<Yes>.

For molecules containing multiple connected components, fingerprints can be generated
in two different ways: use all connected components or just the largest connected
component. By default, all atoms except for the largest connected component are
deleted before generation of fingerprints.

=item B<-m, --mode> I<MACCSKeyBits | MACCSKeyCount>

Specify type of MACCS keys [ Ref 45-47 ] fingerprints to generate for molecules in I<SDFile(s)>.
Possible values: I<MACCSKeyBits, MACCSKeyCount>. Default value: I<MACCSKeyBits>.

For I<MACCSKeyBits> value of B<-m, --mode> option, a fingerprint bit-vector string containing
zeros and ones is generated and for I<MACCSKeyCount> value, a fingerprint vector string
corresponding to number of MACCS keys is generated.

I<MACCSKeyBits | MACCSKeyCount> values for B<-m, --mode> option along with two possible
I<166 | 322>  values of B<-s, --size> supports generation of four different types of MACCS
keys fingerprint: I<MACCS166KeyBits, MACCS166KeyCount, MACCS322KeyBits, MACCS322KeyCount>.

Definition of MACCS keys uses the following atom and bond symbols to define atom and
bond environments:

    Atom symbols for 166 keys [ Ref 47 ]:

    A : Any valid periodic table element symbol
    Q  : Hetro atoms; any non-C or non-H atom
    X  : Halogens; F, Cl, Br, I
    Z  : Others; other than H, C, N, O, Si, P, S, F, Cl, Br, I

    Atom symbols for 322 keys [ Ref 46 ]:

    A : Any valid periodic table element symbol
    Q  : Hetro atoms; any non-C or non-H atom
    X  : Others; other than H, C, N, O, Si, P, S, F, Cl, Br, I
    Z is neither defined nor used

    Bond types:

    -  : Single
    =  : Double
    T  : Triple
    #  : Triple
    ~  : Single or double query bond
    %  : An aromatic query bond

    None : Any bond type; no explicit bond specified

    $  : Ring bond; $ before a bond type specifies ring bond
    !  : Chain or non-ring bond; ! before a bond type specifies chain bond

    @  : A ring linkage and the number following it specifies the
         atoms position in the line, thus @1 means linked back to the first
         atom in the list.

    Aromatic: Kekule or Arom5

    Kekule: Bonds in 6-membered rings with alternate single/double bonds
            or perimeter bonds
    Arom5:  Bonds in 5-membered rings with two double bonds and a hetro
            atom at the apex of the ring.

MACCS 166 keys [ Ref 45-47 ] are defined as follows:

    Key Description

    1	ISOTOPE
    2	103 < ATOMIC NO. < 256
    3	GROUP IVA,VA,VIA PERIODS 4-6 (Ge...)
    4	ACTINIDE
    5	GROUP IIIB,IVB (Sc...)
    6	LANTHANIDE
    7	GROUP VB,VIB,VIIB (V...)
    8	QAAA@1
    9	GROUP VIII (Fe...)
    10	GROUP IIA (ALKALINE EARTH)
    11	4M RING
    12	GROUP IB,IIB (Cu...)
    13	ON(C)C
    14	S-S
    15	OC(O)O
    16	QAA@1
    17	CTC
    18	GROUP IIIA (B...)
    19	7M RING
    20	SI
    21	C=C(Q)Q
    22	3M RING
    23	NC(O)O
    24	N-O
    25	NC(N)N
    26	C$=C($A)$A
    27	I
    28	QCH2Q
    29	P
    30	CQ(C)(C)A
    31	QX
    32	CSN
    33	NS
    34	CH2=A
    35	GROUP IA (ALKALI METAL)
    36	S HETEROCYCLE
    37	NC(O)N
    38	NC(C)N
    39	OS(O)O
    40	S-O
    41	CTN
    42	F
    43	QHAQH
    44	OTHER
    45	C=CN
    46	BR
    47	SAN
    48	OQ(O)O
    49	CHARGE
    50	C=C(C)C
    51	CSO
    52	NN
    53	QHAAAQH
    54	QHAAQH
    55	OSO
    56	ON(O)C
    57	O HETEROCYCLE
    58	QSQ
    59	Snot%A%A
    60	S=O
    61	AS(A)A
    62	A$A!A$A
    63	N=O
    64	A$A!S
    65	C%N
    66	CC(C)(C)A
    67	QS
    68	QHQH (&...)
    69	QQH
    70	QNQ
    71	NO
    72	OAAO
    73	S=A
    74	CH3ACH3
    75	A!N$A
    76	C=C(A)A
    77	NAN
    78	C=N
    79	NAAN
    80	NAAAN
    81	SA(A)A
    82	ACH2QH
    83 	QAAAA@1
    84	NH2
    85	CN(C)C
    86	CH2QCH2
    87	X!A$A
    88	S
    89	OAAAO
    90	QHAACH2A
    91	QHAAACH2A
    92	OC(N)C
    93	QCH3
    94	QN
    95	NAAO
    96	5M RING
    97	NAAAO
    98	QAAAAA@1
    99	C=C
    100	ACH2N
    101	8M RING
    102	QO
    103	CL
    104	QHACH2A
    105	A$A($A)$A
    106	QA(Q)Q
    107	XA(A)A
    108	CH3AAACH2A
    109	ACH2O
    110	NCO
    111	NACH2A
    112	AA(A)(A)A
    113	Onot%A%A
    114	CH3CH2A
    115	CH3ACH2A
    116	CH3AACH2A
    117	NAO
    118	ACH2CH2A > 1
    119	N=A
    120	HETEROCYCLIC ATOM > 1 (&...)
    121	N HETEROCYCLE
    122	AN(A)A
    123	OCO
    124	QQ
    125	AROMATIC RING > 1
    126	A!O!A
    127	A$A!O > 1 (&...)
    128	ACH2AAACH2A
    129	ACH2AACH2A
    130	QQ > 1 (&...)
    131	QH > 1
    132	OACH2A
    133	A$A!N
    134	X (HALOGEN)
    135	Nnot%A%A
    136	O=A > 1
    137	HETEROCYCLE
    138	QCH2A > 1 (&...)
    139	OH
    140	O > 3 (&...)
    141	CH3 > 2 (&...)
    142	N > 1
    143	A$A!O
    144	Anot%A%Anot%A
    145	6M RING > 1
    146	O > 2
    147	ACH2CH2A
    148	AQ(A)A
    149	CH3 > 1
    150	A!A$A!A
    151	NH
    152	OC(C)C
    153	QCH2A
    154	C=O
    155	A!CH2!A
    156	NA(A)A
    157	C-O
    158	C-N
    159	O > 1
    160	CH3
    161	N
    162	AROMATIC
    163	6M RING
    164	O
    165	RING
    166 	FRAGMENTS

MACCS 322 keys set as defined in tables 1, 2 and 3 [ Ref 46 ] include:

    . 26 atom properties of type P, as listed in Table 1
    . 32 one-atom environments, as listed in Table 3
    . 264 atom-bond-atom combinations listed in Table 4

Total number of keys in three tables is : 322

Atom symbol, X, used for 322 keys [ Ref 46 ] doesn't refer to Halogens as it does for 166 keys. In
order to keep the definition of 322 keys consistent with the published definitions, the symbol X is
used to imply "others" atoms, but it's internally mapped to symbol X as defined for 166 keys
during the generation of key values.

Atom properties-based keys (26):

    Key   Description
    1     A(AAA) or AA(A)A - atom with at least three neighbors
    2     Q - heteroatom
    3     Anot%not-A - atom involved in one or more multiple bonds, not aromatic
    4     A(AAAA) or AA(A)(A)A - atom with at least four neighbors
    5     A(QQ) or QA(Q) - atom with at least two heteroatom neighbors
    6     A(QQQ) or QA(Q)Q - atom with at least three heteroatom neighbors
    7     QH - heteroatom with at least one hydrogen attached
    8     CH2(AA) or ACH2A - carbon with at least two single bonds and at least
          two hydrogens attached
    9     CH3(A) or ACH3 - carbon with at least one single bond and at least three
          hydrogens attached
    10    Halogen
    11    A(-A-A-A) or A-A(-A)-A - atom has at least three single bonds
    12    AAAAAA@1 > 2 - atom is in at least two different six-membered rings
    13    A($A$A$A) or A$A($A)$A - atom has more than two ring bonds
    14    A$A!A$A - atom is at a ring/chain boundary. When a comparison is done
          with another atom the path passes through the chain bond.
    15    Anot%A%Anot%A - atom is at an aromatic/nonaromatic boundary. When a
          comparison is done with another atom the path
          passes through the aromatic bond.
    16    A!A!A  - atom with more than one chain bond
    17    A!A$A!A - atom is at a ring/chain boundary. When a comparison is done
          with another atom the path passes through the ring bond.
    18    A%Anot%A%A - atom is at an aromatic/nonaromatic boundary. When a
          comparison is done with another atom the
          path passes through the nonaromatic bond.
    19    HETEROCYCLE - atom is a heteroatom in a ring.
    20    rare properties: atom with five or more neighbors, atom in
          four or more rings, or atom types other than
          H, C, N, O, S, F, Cl, Br, or I
    21    rare properties: atom has a charge, is an isotope, has two or
          more multiple bonds, or has a triple bond.
    22    N - nitrogen
    23    S - sulfur
    24    O - oxygen
    25    A(AA)A(A)A(AA) - atom has two neighbors, each with three or
          more neighbors (including the central atom).
    26    CHACH2 - atom has two hydrocarbon (CH2) neighbors

Atomic environments properties-based keys (32):

    Key   Description
    27    C(CC)
    28    C(CCC)
    29    C(CN)
    30    C(CCN)
    31    C(NN)
    32    C(NNC)
    33    C(NNN)
    34    C(CO)
    35    C(CCO)
    36    C(NO)
    37    C(NCO)
    38    C(NNO)
    39    C(OO)
    40    C(COO)
    41    C(NOO)
    42    C(OOO)
    43    Q(CC)
    44    Q(CCC)
    45    Q(CN)
    46    Q(CCN)
    47    Q(NN)
    48    Q(CNN)
    49    Q(NNN)
    50    Q(CO)
    51    Q(CCO)
    52    Q(NO)
    53    Q(CNO)
    54    Q(NNO)
    55    Q(OO)
    56    Q(COO)
    57    Q(NOO)
    58    Q(OOO)

Note: The first symbol is the central atom, with atoms bonded to the central atom listed in
parentheses. Q is any non-C, non-H atom. If only two atoms are in parentheses, there is
no implication concerning the other atoms bonded to the central atom.

Atom-Bond-Atom properties-based keys: (264)

    Key   Description
    59    C-C
    60    C-N
    61    C-O
    62    C-S
    63    C-Cl
    64    C-P
    65    C-F
    66    C-Br
    67    C-Si
    68    C-I
    69    C-X
    70    N-N
    71    N-O
    72    N-S
    73    N-Cl
    74    N-P
    75    N-F
    76    N-Br
    77    N-Si
    78    N-I
    79    N-X
    80    O-O
    81    O-S
    82    O-Cl
    83    O-P
    84    O-F
    85    O-Br
    86    O-Si
    87    O-I
    88    O-X
    89    S-S
    90    S-Cl
    91    S-P
    92    S-F
    93    S-Br
    94    S-Si
    95    S-I
    96    S-X
    97    Cl-Cl
    98    Cl-P
    99    Cl-F
    100   Cl-Br
    101   Cl-Si
    102   Cl-I
    103   Cl-X
    104   P-P
    105   P-F
    106   P-Br
    107   P-Si
    108   P-I
    109   P-X
    110   F-F
    111   F-Br
    112   F-Si
    113   F-I
    114   F-X
    115   Br-Br
    116   Br-Si
    117   Br-I
    118   Br-X
    119   Si-Si
    120   Si-I
    121   Si-X
    122   I-I
    123   I-X
    124   X-X
    125   C=C
    126   C=N
    127   C=O
    128   C=S
    129   C=Cl
    130   C=P
    131   C=F
    132   C=Br
    133   C=Si
    134   C=I
    135   C=X
    136   N=N
    137   N=O
    138   N=S
    139   N=Cl
    140   N=P
    141   N=F
    142   N=Br
    143   N=Si
    144   N=I
    145   N=X
    146   O=O
    147   O=S
    148   O=Cl
    149   O=P
    150   O=F
    151   O=Br
    152   O=Si
    153   O=I
    154   O=X
    155   S=S
    156   S=Cl
    157   S=P
    158   S=F
    159   S=Br
    160   S=Si
    161   S=I
    162   S=X
    163   Cl=Cl
    164   Cl=P
    165   Cl=F
    166   Cl=Br
    167   Cl=Si
    168   Cl=I
    169   Cl=X
    170   P=P
    171   P=F
    172   P=Br
    173   P=Si
    174   P=I
    175   P=X
    176   F=F
    177   F=Br
    178   F=Si
    179   F=I
    180   F=X
    181   Br=Br
    182   Br=Si
    183   Br=I
    184   Br=X
    185   Si=Si
    186   Si=I
    187   Si=X
    188   I=I
    189   I=X
    190   X=X
    191   C#C
    192   C#N
    193   C#O
    194   C#S
    195   C#Cl
    196   C#P
    197   C#F
    198   C#Br
    199   C#Si
    200   C#I
    201   C#X
    202   N#N
    203   N#O
    204   N#S
    205   N#Cl
    206   N#P
    207   N#F
    208   N#Br
    209   N#Si
    210   N#I
    211   N#X
    212   O#O
    213   O#S
    214   O#Cl
    215   O#P
    216   O#F
    217   O#Br
    218   O#Si
    219   O#I
    220   O#X
    221   S#S
    222   S#Cl
    223   S#P
    224   S#F
    225   S#Br
    226   S#Si
    227   S#I
    228   S#X
    229   Cl#Cl
    230   Cl#P
    231   Cl#F
    232   Cl#Br
    233   Cl#Si
    234   Cl#I
    235   Cl#X
    236   P#P
    237   P#F
    238   P#Br
    239   P#Si
    240   P#I
    241   P#X
    242   F#F
    243   F#Br
    244   F#Si
    245   F#I
    246   F#X
    247   Br#Br
    248   Br#Si
    249   Br#I
    250   Br#X
    251   Si#Si
    252   Si#I
    253   Si#X
    254   I#I
    255   I#X
    256   X#X
    257   C$C
    258   C$N
    259   C$O
    260   C$S
    261   C$Cl
    262   C$P
    263   C$F
    264   C$Br
    265   C$Si
    266   C$I
    267   C$X
    268   N$N
    269   N$O
    270   N$S
    271   N$Cl
    272   N$P
    273   N$F
    274   N$Br
    275   N$Si
    276   N$I
    277   N$X
    278   O$O
    279   O$S
    280   O$Cl
    281   O$P
    282   O$F
    283   O$Br
    284   O$Si
    285   O$I
    286   O$X
    287   S$S
    288   S$Cl
    289   S$P
    290   S$F
    291   S$Br
    292   S$Si
    293   S$I
    294   S$X
    295   Cl$Cl
    296   Cl$P
    297   Cl$F
    298   Cl$Br
    299   Cl$Si
    300   Cl$I
    301   Cl$X
    302   P$P
    303   P$F
    304   P$Br
    305   P$Si
    306   P$I
    307   P$X
    308   F$F
    309   F$Br
    310   F$Si
    311   F$I
    312   F$X
    313   Br$Br
    314   Br$Si
    315   Br$I
    316   Br$X
    317   Si$Si
    318   Si$I
    319   Si$X
    320   I$I
    321   I$X
    322   X$X

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
names: <SDFileName><MACCSKeysFP>.<Ext>. The file type determines <Ext> value.
The sdf, fpf, csv, and tsv <Ext> values are used for SD, FP, comma/semicolon, and tab
delimited text files, respectively.This option is ignored for multiple input files.

=item B<-s, --size> I<number>

Size of MACCS keys [ Ref 45-47 ] set to use during fingerprints generation. Possible values: I<166 or 322>.
Default value: I<166>.

=item B<-v, --VectorStringFormat> I<ValuesString | IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString | ValuesAndIDsPairsString>

Format of fingerprints vector string data in output SD, FP or CSV/TSV text file(s) specified by
B<--output> used during I<MACCSKeyCount> value of B<-m, --mode> option. Possible
values: I<ValuesString, IDsAndValuesString | IDsAndValuesPairsString | ValuesAndIDsString |
ValuesAndIDsPairsString>. Defaultvalue: I<ValuesString>.

Examples:

    FingerprintsVector;MACCSKeyCount;166;OrderedNumericalValues;ValuesStri
    ng;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 1 0 0 3 0 0 0 0 4 0 0 2 0 0 0 0 0 0 0 0 2 0 0 2 0 0 0 0
    0 0 0 0 1 1 8 0 0 0 1 0 0 1 0 1 0 1 0 3 1 3 1 0 0 0 1 2 0 11 1 0 0 0
    5 0 0 1 2 0 1 1 0 0 0 0 0 1 1 0 1 1 1 1 0 4 0 0 1 1 0 4 6 1 1 1 2 1 1
    3 5 2 2 0 5 3 5 1 1 2 5 1 2 1 2 4 8 3 5 5 2 2 0 3 5 4 1

    FingerprintsVector;MACCSKeyCount;322;OrderedNumericalValues;ValuesStri
    ng;14 8 2 0 2 0 4 4 2 1 4 0 0 2 5 10 5 2 1 0 0 2 0 5 13 3 28 5 5 3 0 0
    0 4 2 1 1 0 1 1 0 0 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 22 5 3 0 0 0 1 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 2 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...

=item B<-w, --WorkingDir> I<DirName>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To generate MACCS keys fingerprints of size 166 in binary bit-vector string format
and create a SampleMACCS166FPBin.csv file containing sequential compound IDs along with
fingerprints bit-vector strings data, type:

    % MACCSKeysFingerprints.pl -r SampleMACCS166FPBin -o Sample.sdf

To generate MACCS keys fingerprints of size 166 in binary bit-vector string format
and create SampleMACCS166FPBin.sdf, SampleMACCS166FPBin.csv and SampleMACCS166FPBin.csv
files containing sequential compound IDs in CSV file along with fingerprints bit-vector strings data, type:

    % MACCSKeysFingerprints.pl --output all -r SampleMACCS166FPBin
      -o Sample.sdf

To generate MACCS keys fingerprints of size 322 in binary bit-vector string format
and create a SampleMACCS322FPBin.csv file containing sequential compound IDs along with
fingerprints bit-vector strings data, type:

    % MACCSKeysFingerprints.pl -size 322 -r SampleMACCS322FPBin -o Sample.sdf

To generate MACCS keys fingerprints of size 166 corresponding to count of keys in
ValuesString format and create a SampleMACCS166FPCount.csv file containing sequential
compound IDs along with fingerprints vector strings data, type:

    % MACCSKeysFingerprints.pl -m MACCSKeyCount -r SampleMACCS166FPCount
      -o Sample.sdf

To generate MACCS keys fingerprints of size 322 corresponding to count of keys in
ValuesString format and create a SampleMACCS322FPCount.csv file containing sequential
compound IDs along with fingerprints vector strings data, type:

    % MACCSKeysFingerprints.pl -m MACCSKeyCount -size 322
      -r SampleMACCS322FPCount -o Sample.sdf

To generate MACCS keys fingerprints of size 166 in hexadecimal bit-vector string format with
ascending bits order and create a SampleMACCS166FPHex.csv file containing compound IDs
from MolName along with fingerprints bit-vector strings data, type:

    % MACCSKeysFingerprints.pl -m MACCSKeyBits --size 166 --BitStringFormat
      HexadecimalString --BitsOrder Ascending --DataFieldsMode CompoundID
      --CompoundIDMode MolName -r SampleMACCS166FPBin -o Sample.sdf

To generate MACCS keys fingerprints of size 166 corresponding to count of keys in
IDsAndValuesString format and create a SampleMACCS166FPCount.csv file containing
compound IDs from MolName line along with fingerprints vector strings data, type:

    % MACCSKeysFingerprints.pl -m MACCSKeyCount --size 166
      --VectorStringFormat IDsAndValuesString  --DataFieldsMode CompoundID
      --CompoundIDMode MolName -r SampleMACCS166FPCount -o Sample.sdf

To generate MACCS keys fingerprints of size 166 corresponding to count of keys in
IDsAndValuesString format and create a SampleMACCS166FPCount.csv file containing
compound IDs using specified data field along with fingerprints vector strings data, type:

    % MACCSKeysFingerprints.pl -m MACCSKeyCount --size 166
      --VectorStringFormat IDsAndValuesString  --DataFieldsMode CompoundID
      --CompoundIDMode DataField --CompoundID Mol_ID -r
      SampleMACCS166FPCount -o Sample.sdf

To generate MACCS keys fingerprints of size 322 corresponding to count of keys in
ValuesString format and create a SampleMACCS322FPCount.tsv file containing compound
IDs derived from combination of molecule name line and an explicit compound prefix
along with fingerprints vector strings data in a column labels MACCSKeyCountFP, type:

    % MACCSKeysFingerprints.pl -m MACCSKeyCount -size 322 --DataFieldsMode
      CompoundID --CompoundIDMode MolnameOrLabelPrefix --CompoundID Cmpd
      --CompoundIDLabel MolID --FingerprintsLabel MACCSKeyCountFP --OutDelim
      Tab -r SampleMACCS322FPCount -o Sample.sdf

To generate MACCS keys fingerprints of size 166 corresponding to count of keys in
ValuesString format and create a SampleMACCS166FPCount.csv file containing
specific data fields columns along with fingerprints vector strings data, type:

    % MACCSKeysFingerprints.pl -m MACCSKeyCount --size 166
      --VectorStringFormat ValuesString --DataFieldsMode Specify --DataFields
      Mol_ID  -r SampleMACCS166FPCount -o Sample.sdf

To generate MACCS keys fingerprints of size 322 corresponding to count of keys in
ValuesString format and create a SampleMACCS322FPCount.csv file containing
common data fields columns along with fingerprints vector strings data, type:

    % MACCSKeysFingerprints.pl -m MACCSKeyCount --size 322
      --VectorStringFormat ValuesString --DataFieldsMode Common -r
      SampleMACCS322FPCount -o Sample.sdf

To generate MACCS keys fingerprints of size 166 corresponding to count of keys in
ValuesString format and create SampleMACCS166FPCount.sdf, SampleMACCS166FPCount.fpf and
SampleMACCS166FPCount.csv files containing all data fields columns in CSV file
along with fingerprints vector strings data, type:

    % MACCSKeysFingerprints.pl -m MACCSKeyCount --size 166 --output all
      --VectorStringFormat ValuesString --DataFieldsMode All -r
      SampleMACCS166FPCount -o Sample.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoFingerprintsFiles.pl, SimilarityMatricesFingerprints.pl, AtomNeighborhoodsFingerprints.pl,
ExtendedConnectivityFingerprints.pl, PathLengthFingerprints.pl,
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
