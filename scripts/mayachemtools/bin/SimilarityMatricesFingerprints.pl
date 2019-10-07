#!/usr/bin/perl -w
#
# File: SimilarityMatricesFingerprints.pl
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
use File::Copy;
use Text::ParseWords;
use Benchmark;
use FileUtil;
use TextUtil;
use Fingerprints::FingerprintsFileUtil;
use Fingerprints::FingerprintsBitVector;
use Fingerprints::FingerprintsVector;

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

my(@FingerprintsFilesList);
@FingerprintsFilesList = ExpandFileNames(\@ARGV, "sdf sd fpf fp csv tsv");

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Setup information about input files...
print "Checking input fingerprints file(s)...\n";
my(%FingerprintsFilesInfo);
RetrieveFingerprintsFilesInfo();

# Process input files..
my($FileIndex);
if (@FingerprintsFilesList > 1) {
  print "\nProcessing fingerprints files...\n";
}
for $FileIndex (0 .. $#FingerprintsFilesList) {
  if ($FingerprintsFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $FingerprintsFilesList[$FileIndex]...\n";
    GenerateSimilarityMatrices($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Generate similarity matrices using fingerprints data in text file...
#
sub GenerateSimilarityMatrices {
  my($FileIndex) = @_;

  ProcessFingerprintsData($FileIndex);

  if ($FingerprintsFilesInfo{FingerprintsBitVectorStringMode}[$FileIndex]) {
    GenerateSimilarityMatricesForFingerprintsBitVectors($FileIndex);
  }
  elsif ($FingerprintsFilesInfo{FingerprintsVectorStringMode}[$FileIndex]) {
    GenerateSimilarityMatricesForFingerprintsVectors($FileIndex);
  }

  CleanupFingerprintsData($FileIndex);
}

# Generate bit vector similarity matrices...
#
sub GenerateSimilarityMatricesForFingerprintsBitVectors {
  my($FileIndex) = @_;
  my($SpecifiedComparisonMeasure, $ComparisonMeasure, $NewTextFile, $SimilarityMatrixRef, $MethodName, @MethodParameters);

  for $SpecifiedComparisonMeasure (@{$OptionsInfo{SpecifiedBitVectorComparisonsRef}}) {
    $ComparisonMeasure = $OptionsInfo{SpecifiedBitVectorComparisonsNameRef}->{lc($SpecifiedComparisonMeasure)};
    $NewTextFile = $FingerprintsFilesInfo{OutFileRoot}[$FileIndex] . "${ComparisonMeasure}." . $FingerprintsFilesInfo{OutFileExt}[$FileIndex];

    $MethodName = $OptionsInfo{SpecifiedBitVectorComparisonsMethodRef}->{lc($ComparisonMeasure)};

    @MethodParameters = ();
    @MethodParameters = @{$OptionsInfo{SpecifiedBitVectorComparisonsParameterRef}->{lc($ComparisonMeasure)}};

    GenerateSimilarityMatrix($FileIndex, $NewTextFile, $MethodName, \@MethodParameters);
  }
}

# Generate vector similarity and/or distance matrices...
#
sub GenerateSimilarityMatricesForFingerprintsVectors {
  my($FileIndex) = @_;
  my($SpecifiedComparisonMeasure, $ComparisonMode, $ComparisonMeasure, $NewTextFile, $MethodName, @MethodParameters);

  for $SpecifiedComparisonMeasure (@{$OptionsInfo{SpecifiedVectorComparisonsRef}}) {
    $ComparisonMeasure = $OptionsInfo{SpecifiedVectorComparisonsNameRef}->{lc($SpecifiedComparisonMeasure)};

    for $ComparisonMode (@{$OptionsInfo{SpecifiedVectorComparisonModesRef}}) {
      $NewTextFile = $FingerprintsFilesInfo{OutFileRoot}[$FileIndex] . "${ComparisonMeasure}${ComparisonMode}." . $FingerprintsFilesInfo{OutFileExt}[$FileIndex];

      $MethodName = $OptionsInfo{SpecifiedVectorComparisonsMethodRef}->{lc($ComparisonMeasure)};

      @MethodParameters = ();
      push @MethodParameters, $ComparisonMode;
      push @MethodParameters, @{$OptionsInfo{SpecifiedVectorComparisonsParameterRef}->{lc($ComparisonMeasure)}};

      GenerateSimilarityMatrix($FileIndex, $NewTextFile, $MethodName, \@MethodParameters);
    }
  }
}

# Calculate similarity matrix and write it out...
#
sub GenerateSimilarityMatrix {
  my($FileIndex, $NewTextFile, $MethodName, $MethodParametersRef) = @_;

  print "\nGenerating $NewTextFile...\n";

  # Open new file and write out column labels...
  open NEWTEXTFILE, ">$NewTextFile" or die "Error: Can't open $NewTextFile: $! \n";
  WriteColumnLabels($FileIndex, \*NEWTEXTFILE);

  # Calculate and write out similarity matrix values...
  if ($OptionsInfo{InputDataMode} =~ /^LoadInMemory$/i) {
    GenerateSimilarityMatrixUsingMemoryData($FileIndex, \*NEWTEXTFILE, $MethodName, $MethodParametersRef);
  }
  elsif ($OptionsInfo{InputDataMode} =~ /^ScanFile$/i) {
    GenerateSimilarityMatrixUsingFileData($FileIndex, \*NEWTEXTFILE, $MethodName, $MethodParametersRef);
  }
  else {
    warn "Warning: Input data mode, $OptionsInfo{InputDataMode}, is not supported.\n";
  }

  # Close new text file...
  close NEWTEXTFILE;

}

# Calculate and write out similarity values using fingerprints data already loaded in
# memory...
#
sub GenerateSimilarityMatrixUsingMemoryData {
  my($FileIndex, $NewTextFileRef, $MethodName, $MethodParametersRef) = @_;
  my($RowIndex, $ColIndex, $CmpdID1, $CmpdID2, $FingerprintsObject1, $FingerprintsObject2, $Value, $Line, @LineWords);

  for $RowIndex (0 .. $#{$FingerprintsFilesInfo{FingerprintsObjectsRef}}) {
    $FingerprintsObject1 = $FingerprintsFilesInfo{FingerprintsObjectsRef}->[$RowIndex];
    $CmpdID1 = $FingerprintsFilesInfo{CompundIDsRef}->[$RowIndex];

    if ($OptionsInfo{WriteRowsAndColumns}) {
      print $NewTextFileRef "$OptionsInfo{OutQuoteValue}${CmpdID1}$OptionsInfo{OutQuoteValue}";
    }

    COLINDEX: for $ColIndex (0 .. $#{$FingerprintsFilesInfo{FingerprintsObjectsRef}}) {
      if (SkipMatrixData($RowIndex, $ColIndex)) {
	next COLINDEX;
      }

      $FingerprintsObject2 = $FingerprintsFilesInfo{FingerprintsObjectsRef}->[$ColIndex];

      $Value = $FingerprintsObject1->$MethodName($FingerprintsObject2, @{$MethodParametersRef});
      $Value = (defined($Value) && length($Value)) ? (sprintf("%.$OptionsInfo{Precision}f", $Value) + 0) : '';

      if ($OptionsInfo{WriteRowsAndColumns}) {
	print $NewTextFileRef "$OptionsInfo{OutDelim}$OptionsInfo{OutQuoteValue}${Value}$OptionsInfo{OutQuoteValue}";
      }
      elsif ($OptionsInfo{WriteIDPairsAndValue}) {
	$CmpdID2 = $FingerprintsFilesInfo{CompundIDsRef}->[$ColIndex];

	@LineWords = ();
	push @LineWords,  ($CmpdID1, $CmpdID2, $Value);
	$Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
	print $NewTextFileRef "$Line\n";
      }
    }
    if ($OptionsInfo{WriteRowsAndColumns}) {
      print $NewTextFileRef "\n";
    }
  }
}

# Calculate and write out similarity values by retrieving and prcessing data
# from fingerprint file...
#
sub GenerateSimilarityMatrixUsingFileData {
  my($FileIndex, $NewTextFileRef, $MethodName, $MethodParametersRef) = @_;
  my($RowIndex, $ColIndex, $FingerprintsFileIO, $TmpFingerprintsFileIO, $FingerprintsObject1, $FingerprintsObject2, $CmpdID1, $CmpdID2, $FingerprintsCount, $IgnoredFingerprintsCount, $Value, $Line, @LineWords);

  print "\nReading and processing fingerprints data...\n";

  $FingerprintsFileIO = Fingerprints::FingerprintsFileUtil::NewFingerprintsFileIO(%{$FingerprintsFilesInfo{FingerprintsFileIOParameters}[$FileIndex]});
  $FingerprintsFileIO->Open();

  $RowIndex = 0; $ColIndex = 0;
  $FingerprintsCount = 0; $IgnoredFingerprintsCount = 0;

  FINGERPRINTSFILEIO: while ($FingerprintsFileIO->Read()) {
    $FingerprintsCount++;

    if (!$FingerprintsFileIO->IsFingerprintsDataValid()) {
      $IgnoredFingerprintsCount++;
      next FINGERPRINTSFILEIO;
    }
    $RowIndex++;
    $FingerprintsObject1 = $FingerprintsFileIO->GetFingerprints();
    $CmpdID1 = $FingerprintsFileIO->GetCompoundID();

    if ($OptionsInfo{WriteRowsAndColumns}) {
      print $NewTextFileRef "$OptionsInfo{OutQuoteValue}${CmpdID1}$OptionsInfo{OutQuoteValue}";
    }

    # Force detail level of 1 to avoid duplicate printing of diagnostic messages for invalid
    # fingerprints data...
    $TmpFingerprintsFileIO = Fingerprints::FingerprintsFileUtil::NewFingerprintsFileIO(%{$FingerprintsFilesInfo{TmpFingerprintsFileIOParameters}[$FileIndex]}, "DetailLevel" => 1);
    $TmpFingerprintsFileIO->Open();

    $ColIndex = 0;
    TMPFINGERPRINTSFILEIO: while ($TmpFingerprintsFileIO->Read()) {
      if (!$TmpFingerprintsFileIO->IsFingerprintsDataValid()) {
	next TMPFINGERPRINTSFILEIO;
      }
      $ColIndex++;

      if (SkipMatrixData($RowIndex, $ColIndex)) {
	next TMPFINGERPRINTSFILEIO;
      }

      $FingerprintsObject2 = $TmpFingerprintsFileIO->GetFingerprints();

      $Value = $FingerprintsObject1->$MethodName($FingerprintsObject2, @{$MethodParametersRef});
      $Value = (defined($Value) && length($Value)) ? (sprintf("%.$OptionsInfo{Precision}f", $Value) + 0) : '';

      if ($OptionsInfo{WriteRowsAndColumns}) {
	print $NewTextFileRef "$OptionsInfo{OutDelim}$OptionsInfo{OutQuoteValue}${Value}$OptionsInfo{OutQuoteValue}";
      }
      elsif ($OptionsInfo{WriteIDPairsAndValue}) {
	$CmpdID2 = $TmpFingerprintsFileIO->GetCompoundID();

	@LineWords = ();
	push @LineWords,  ($CmpdID1, $CmpdID2, $Value);
	$Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
	print $NewTextFileRef "$Line\n";
      }
    }
    $TmpFingerprintsFileIO->Close();

    if ($OptionsInfo{WriteRowsAndColumns}) {
      print $NewTextFileRef "\n";
    }
  }

  $FingerprintsFileIO->Close();

  print "Number of fingerprints data entries in database fingerprints file: $FingerprintsCount\n";
  print "Number of fingerprints date entries processed successfully: ", ($FingerprintsCount - $IgnoredFingerprintsCount)  , "\n";
  print "Number of fingerprints data entries ignored due to missing/invalid data: $IgnoredFingerprintsCount\n\n";
}

# Check whether matrix data need to be skipped...
#
sub SkipMatrixData {
  my($RowIndex, $ColIndex) = @_;

  if ($OptionsInfo{WriteFullMatrix}) {
    return 0;
  }
  elsif ($OptionsInfo{WriteUpperTriangularMatrix}) {
    return ($RowIndex > $ColIndex) ? 1 : 0;
  }
  elsif ($OptionsInfo{WriteLowerTriangularMatrix}) {
    return ($RowIndex < $ColIndex) ? 1 : 0;
  }

  return 0;
}

# Write out column labels...
#
sub WriteColumnLabels {
  my($FileIndex, $NewTextFileRef) = @_;
  my($Line, @LineWords);

  if ($OptionsInfo{OutMatrixFormat} =~ /^IDPairsAndValue$/i) {
    @LineWords = ();
    push @LineWords, ('CmpdID1', 'CmpdID2', 'Coefficient Value');
    $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print $NewTextFileRef "$Line\n";
  }
  elsif ($OptionsInfo{OutMatrixFormat} =~ /^RowsAndColumns$/i) {
    if ($OptionsInfo{InputDataMode} =~ /^LoadInMemory$/i) {
      @LineWords = ();
      push @LineWords, '';
      push @LineWords, @{$FingerprintsFilesInfo{CompundIDsRef}};
      $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
      print $NewTextFileRef "$Line\n";
    }
    elsif ($OptionsInfo{InputDataMode} =~ /^ScanFile$/i) {
      my( $FingerprintsFileIO, $CmpdID);

      # Scan file to retrieve compound IDs...
      #
      print "\nProcessing fingerprints file to generate compound IDs...\n";

      # Force detail level of 1 to avoid diagnostics messages for invalid fingeprints data during
      # retrieval of compound IDs as these get printed out during calculation of matrix...
      #
      $FingerprintsFileIO = Fingerprints::FingerprintsFileUtil::NewFingerprintsFileIO(%{$FingerprintsFilesInfo{FingerprintsFileIOParameters}[$FileIndex]}, "DetailLevel" => 1);
      $FingerprintsFileIO->Open();

      print $NewTextFileRef "$OptionsInfo{OutQuoteValue}$OptionsInfo{OutQuoteValue}";

      FINGERPRINTSFILEIO: while ($FingerprintsFileIO->Read()) {
	if (!$FingerprintsFileIO->IsFingerprintsDataValid()) {
	  next FINGERPRINTSFILEIO;
	}
	$CmpdID = $FingerprintsFileIO->GetCompoundID();
	print $NewTextFileRef "$OptionsInfo{OutDelim}$OptionsInfo{OutQuoteValue}${CmpdID}$OptionsInfo{OutQuoteValue}";
      }
      $FingerprintsFileIO->Close();

      print $NewTextFileRef "\n";

      print "Processing fingerprints file to generate matrix...\n";
    }
  }
  else {
    warn "Warning: Output matrix format, $OptionsInfo{OutMatrixFormat}, is not supported.\n";
  }
}

# Process fingerprints data...
#
sub ProcessFingerprintsData {
  my($FileIndex) = @_;
  my($FingerprintsFileIO);

  $FingerprintsFilesInfo{CompundIDsRef}  = undef;
  $FingerprintsFilesInfo{FingerprintsObjectsRef} = undef;

  if ($OptionsInfo{InputDataMode} =~ /^LoadInMemory$/i) {
    my($FingerprintsFileIO);

    $FingerprintsFileIO = Fingerprints::FingerprintsFileUtil::NewFingerprintsFileIO(%{$FingerprintsFilesInfo{FingerprintsFileIOParameters}[$FileIndex]});
    ($FingerprintsFilesInfo{CompundIDsRef}, $FingerprintsFilesInfo{FingerprintsObjectsRef}) = Fingerprints::FingerprintsFileUtil::ReadAndProcessFingerpritsData($FingerprintsFileIO);
  }
  elsif ($OptionsInfo{InputDataMode} =~ /^ScanFile$/i) {
    my($FingerprintsFile, $TmpFingerprintsFile);

    $FingerprintsFile = $FingerprintsFilesList[$FileIndex];
    $TmpFingerprintsFile = $FingerprintsFilesInfo{TmpFingerprintsFile}[$FileIndex];

    # Copy fingerprints file to a tmp file for calculating similarity matrix...
    print "\nCopying fingerprints file, $FingerprintsFile, to temporary fingperints file, $TmpFingerprintsFile...\n";
    copy $FingerprintsFile, $TmpFingerprintsFile or die "Error: Couldn't copy $FingerprintsFile to $TmpFingerprintsFile: $! \n";
  }
}

# Clean up fingerprints data...
#
sub CleanupFingerprintsData {
  my($FileIndex) = @_;

  if ($OptionsInfo{InputDataMode} =~ /^LoadInMemory$/i) {
    $FingerprintsFilesInfo{CompundIDsRef}  = undef;
    $FingerprintsFilesInfo{FingerprintsObjectsRef} = undef;
  }
  elsif ($OptionsInfo{InputDataMode} =~ /^ScanFile$/i) {
    my($TmpFingerprintsFile);

    # Delete temporary fingerprints file...
    $TmpFingerprintsFile = $FingerprintsFilesInfo{TmpFingerprintsFile}[$FileIndex];

    print "\nDeleting temporary fingerprints file $TmpFingerprintsFile...\n";
    unlink $TmpFingerprintsFile or die "Error: Couldn't unlink $TmpFingerprintsFile: $! \n";
  }
}

# Retrieve information about fingerprints files...
#
sub RetrieveFingerprintsFilesInfo {
  my($FingerprintsFile, $TmpFingerprintsFile, $FingerprintsFileIO, $FingerprintsBitVectorStringMode, $FingerprintsVectorStringMode, $FileType, $Index, $FileDir, $FileExt, $FileName, $InDelim, $OutFileRoot, $OutFileExt, %FingerprintsFileIOParameters);

  %FingerprintsFilesInfo = ();
  @{$FingerprintsFilesInfo{FileOkay}} = ();
  @{$FingerprintsFilesInfo{FileType}} = ();
  @{$FingerprintsFilesInfo{InDelim}} = ();
  @{$FingerprintsFilesInfo{OutFileRoot}} = ();
  @{$FingerprintsFilesInfo{OutFileExt}} = ();

  @{$FingerprintsFilesInfo{TmpFingerprintsFile}} = ();

  @{$FingerprintsFilesInfo{FingerprintsFileIOParameters}} = ();
  @{$FingerprintsFilesInfo{TmpFingerprintsFileIOParameters}} = ();

  @{$FingerprintsFilesInfo{FingerprintsBitVectorStringMode}} = ();
  @{$FingerprintsFilesInfo{FingerprintsVectorStringMode}} = ();

  FILELIST: for $Index (0 .. $#FingerprintsFilesList) {
    $FingerprintsFilesInfo{FileOkay}[$Index] = 0;
    $FingerprintsFilesInfo{FileType}[$Index] = '';
    $FingerprintsFilesInfo{InDelim}[$Index] = "";
    $FingerprintsFilesInfo{OutFileRoot}[$Index] = '';
    $FingerprintsFilesInfo{OutFileExt}[$Index] = '';

    %{$FingerprintsFilesInfo{FingerprintsFileIOParameters}[$Index]} = ();

    $FingerprintsFilesInfo{TmpFingerprintsFile}[$Index] = "";
    %{$FingerprintsFilesInfo{TmpFingerprintsFileIOParameters}[$Index]} = ();

    $FingerprintsFilesInfo{FingerprintsBitVectorStringMode}[$Index] = 0;
    $FingerprintsFilesInfo{FingerprintsVectorStringMode}[$Index] = 0;

    $FingerprintsFile = $FingerprintsFilesList[$Index];
    if (!(-e $FingerprintsFile)) {
      warn "Warning: Ignoring file $FingerprintsFile: It doesn't exist\n";
      next FILELIST;
    }

    $FileType = Fingerprints::FingerprintsFileUtil::GetFingerprintsFileType($FingerprintsFile);
    if (IsEmpty($FileType)) {
      warn "Warning: Ignoring file $FingerprintsFile: It's not a fingerprints file\n";
      next FILELIST;
    }

    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($FingerprintsFile);

    # Setup temporary fingerprints file name for scan file mode...
    $TmpFingerprintsFile = "${FileName}Tmp.${FileExt}";

    $InDelim = ($FileExt =~ /^tsv$/i) ? 'Tab' : $OptionsInfo{InDelim};

    # Setup output file names...
    $OutFileExt = "csv";
    if ($Options{outdelim} =~ /^tab$/i) {
      $OutFileExt = "tsv";
    }

    $OutFileRoot = $FileName;
    if ($OptionsInfo{OutFileRoot} && (@FingerprintsFilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($OptionsInfo{OutFileRoot});
      if ($RootFileName && $RootFileExt) {
	$FileName = $RootFileName;
      }
      else {
	$FileName = $OptionsInfo{OutFileRoot};
      }
      $OutFileRoot = $FileName;
    }

    if (!$Options{overwrite}) {
      # Similarity matrices output file names for bit-vector strings...
      my($SpecifiedComparisonMeasure, $ComparisonMeasure);
      for $SpecifiedComparisonMeasure (@{$OptionsInfo{SpecifiedBitVectorComparisonsRef}}) {
	$ComparisonMeasure = $OptionsInfo{SpecifiedBitVectorComparisonsNameRef}->{lc($SpecifiedComparisonMeasure)};
	if (-e "${OutFileRoot}${ComparisonMeasure}.${OutFileExt}") {
	  warn "Warning: Ignoring file $FingerprintsFile: The file ${OutFileRoot}${ComparisonMeasure}.${OutFileExt} already exists.\n";
	  next FILELIST;
	}
      }
      # Similarity matrices output file names for vector strings...
      my($ComparisonMode);
      for $SpecifiedComparisonMeasure (@{$OptionsInfo{SpecifiedVectorComparisonsRef}}) {
	$ComparisonMeasure = $OptionsInfo{SpecifiedVectorComparisonsNameRef}->{lc($SpecifiedComparisonMeasure)};
	for $ComparisonMode (@{$OptionsInfo{SpecifiedVectorComparisonModesRef}}) {
	  if (-e "${OutFileRoot}${ComparisonMeasure}${ComparisonMode}.${OutFileExt}") {
	    warn "Warning: Ignoring file $FingerprintsFile: The file ${OutFileRoot}${ComparisonMeasure}${ComparisonMode}.${OutFileExt} already exists.\n";
	    next FILELIST;
	  }
	}
      }
    }

    # Setup FingerprintsFileIO parameters...
    %FingerprintsFileIOParameters = ();
    FILEIOPARAMETERS: {
      if ($FileType =~ /^SD$/i) {
	%FingerprintsFileIOParameters = ('Name' => $FingerprintsFile, 'Mode' => 'Read', 'FingerprintsStringMode' => $OptionsInfo{Mode}, 'ValidateData' => $OptionsInfo{ValidateData}, 'DetailLevel' =>  $OptionsInfo{Detail}, 'FingerprintsFieldLabel' => $OptionsInfo{FingerprintsField}, 'CompoundIDMode' => $OptionsInfo{CompoundIDMode}, 'CompoundIDFieldLabel' => $OptionsInfo{CompoundIDField}, 'CompoundIDPrefix' => $OptionsInfo{CompoundIDPrefix});
	last FILEIOPARAMETERS;
      }
      if ($FileType =~ /^FP$/i) {
	%FingerprintsFileIOParameters = ('Name' => $FingerprintsFile, 'Mode' => 'Read', 'FingerprintsStringMode' => $OptionsInfo{Mode}, 'ValidateData' => $OptionsInfo{ValidateData}, 'DetailLevel' =>  $OptionsInfo{Detail});
	last FILEIOPARAMETERS;
      }
      if ($FileType =~ /^Text$/i) {
	%FingerprintsFileIOParameters = ('Name' => $FingerprintsFile, 'Mode' => 'Read', 'FingerprintsStringMode' => $OptionsInfo{Mode}, 'ValidateData' => $OptionsInfo{ValidateData}, 'DetailLevel' =>  $OptionsInfo{Detail}, 'FingerprintsCol' => $OptionsInfo{FingerprintsCol}, 'ColMode' => $OptionsInfo{ColMode}, 'CompoundIDCol' => $OptionsInfo{CompoundIDCol}, 'CompoundIDPrefix' => $OptionsInfo{CompoundIDPrefix}, 'InDelim' => $OptionsInfo{InDelim});
	last FILEIOPARAMETERS;
      }
      warn "Warning: File type for fingerprints file, $FingerprintsFile, is not valid. Supported file types: SD, FP or Text\n";
      next FILELIST;
    }

    # Retrieve fingerints file string mode information...
    $FingerprintsFileIO = Fingerprints::FingerprintsFileUtil::NewFingerprintsFileIO(%FingerprintsFileIOParameters);

    if (!$FingerprintsFileIO) {
      warn "Warning: Ignoring fingerprints file $FingerprintsFile: It contains invalid fingerprints data\n";
      next FILELIST;
    }
    if (!$FingerprintsFileIO->IsFingerprintsFileDataValid()) {
      warn "Warning: Ignoring fingerprints file $FingerprintsFile: It contains invalid fingerprints data\n";
      next FILELIST;
    }
    $FingerprintsBitVectorStringMode = $FingerprintsFileIO->GetFingerprintsBitVectorStringMode();
    $FingerprintsVectorStringMode = $FingerprintsFileIO->GetFingerprintsVectorStringMode();


    $FingerprintsFilesInfo{FileOkay}[$Index] = 1;
    $FingerprintsFilesInfo{FileType}[$Index] = $FileType;

    $FingerprintsFilesInfo{InDelim}[$Index] = $InDelim;

    $FingerprintsFilesInfo{OutFileRoot}[$Index] = $OutFileRoot;
    $FingerprintsFilesInfo{OutFileExt}[$Index] = $OutFileExt;

    %{$FingerprintsFilesInfo{FingerprintsFileIOParameters}[$Index]} = %FingerprintsFileIOParameters;

    $FingerprintsFilesInfo{TmpFingerprintsFile}[$Index] = $TmpFingerprintsFile;

    $FingerprintsFileIOParameters{Name} = $TmpFingerprintsFile;
    %{$FingerprintsFilesInfo{TmpFingerprintsFileIOParameters}[$Index]} = %FingerprintsFileIOParameters;

    $FingerprintsFilesInfo{FingerprintsBitVectorStringMode}[$Index] = $FingerprintsBitVectorStringMode;
    $FingerprintsFilesInfo{FingerprintsVectorStringMode}[$Index] = $FingerprintsVectorStringMode;
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{InputDataMode} = $Options{inputdatamode};

  ProcessBitVectorComparisonOptions();
  ProcessVectorComparisonOptions();

  $OptionsInfo{CompoundIDPrefix} = $Options{compoundidprefix} ? $Options{compoundidprefix} : 'Cmpd';

  # Compound ID and fingerprints column options for text files...
  $OptionsInfo{ColMode} = $Options{colmode};

  if (IsNotEmpty($Options{compoundidcol})) {
    if ($Options{colmode} =~ /^ColNum$/i) {
      if (!IsPositiveInteger($Options{compoundidcol})) {
	die "Error: Column value, $Options{compoundidcol}, specified using \"--CompoundIDCol\" is not valid: Allowed integer values: > 0\n";
      }
    }
    $OptionsInfo{CompoundIDCol} = $Options{compoundidcol};
  }
  else {
    $OptionsInfo{CompoundIDCol} = 'AutoDetect';
  }

  if (IsNotEmpty($Options{fingerprintscol})) {
    if ($Options{colmode} =~ /^ColNum$/i) {
      if (!IsPositiveInteger($Options{fingerprintscol})) {
	die "Error: Column value, $Options{fingerprintscol}, specified using \"--FingerprintsCol\" is not valid: Allowed integer values: > 0\n";
      }
    }
    $OptionsInfo{FingerprintsCol} = $Options{fingerprintscol};
  }
  else {
    $OptionsInfo{FingerprintsCol} = 'AutoDetect';
  }

  if (IsNotEmpty($Options{compoundidcol}) && IsNotEmpty($Options{fingerprintscol})) {
    if (IsPositiveInteger($Options{compoundidcol}) && IsPositiveInteger($Options{fingerprintscol})) {
      if (($Options{compoundidcol} == $Options{fingerprintscol})) {
	die "Error: Values specified using \"--CompoundIDCol\" and \"--FingerprintsCol\", $Options{compoundidcol}, must be different.\n";
      }
    }
    else {
      if (($Options{compoundidcol} eq $Options{fingerprintscol})) {
	die "Error: Values specified using \"--CompoundIDCol\" and \"--FingerprintsCol\", $Options{compoundidcol}, must be different.\n";
      }
    }
  }

  # Compound ID and fingerprints field options for SD files...
  $OptionsInfo{CompoundIDMode} = $Options{compoundidmode};
  $OptionsInfo{CompoundIDField} = '';

  if ($Options{compoundidmode} =~ /^DataField$/i) {
    if (!$Options{compoundidfield}) {
      die "Error: You must specify a value for \"--CompoundIDField\" option in \"DataField\" \"--CompoundIDMode\". \n";
    }
    $OptionsInfo{CompoundIDField} = $Options{compoundidfield};
  }


  if (IsNotEmpty($Options{fingerprintsfield})) {
    $OptionsInfo{FingerprintsField} = $Options{fingerprintsfield};
  }
  else {
    $OptionsInfo{FingerprintsField} = 'AutoDetect';
  }

  if ($Options{compoundidfield} && IsNotEmpty($Options{fingerprintsfield})) {
    if (($Options{compoundidfield} eq $Options{fingerprintsfield})) {
      die "Error: Values specified using \"--CompoundIDField\" and \"--Fingerprintsfield\", $Options{compoundidfield}, must be different.\n";
    }
  }

  $OptionsInfo{Detail} = $Options{detail};

  $OptionsInfo{InDelim} = $Options{indelim};
  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /tab/i ) ? "\t" : (($Options{outdelim} =~ /semicolon/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^Yes$/i) ? 1 : 0;
  $OptionsInfo{OutQuoteValue} = ($Options{quote} =~ /^Yes$/i) ? '"' : '';

  $OptionsInfo{OutMatrixFormat} = $Options{outmatrixformat};

  $OptionsInfo{WriteRowsAndColumns} = 0; $OptionsInfo{WriteIDPairsAndValue} = 0;
  OUTMATRIXFORMAT: {
    if ($OptionsInfo{OutMatrixFormat} =~ /^RowsAndColumns$/i) {
      $OptionsInfo{WriteRowsAndColumns} = 1; last OUTMATRIXFORMAT;
    }
    if ($OptionsInfo{OutMatrixFormat} =~ /^IDPairsAndValue$/i) {
      $OptionsInfo{WriteIDPairsAndValue} = 1; last OUTMATRIXFORMAT;
    }
    die "Error: The value specified, $Options{outmatrixformat}, for option \"--OutMatrixFormat\" is not valid. Allowed values: RowsAndColumns or IDPairsAndValue\n";
  }

  $OptionsInfo{OutMatrixType} = $Options{outmatrixtype};

  $OptionsInfo{WriteFullMatrix} = 0;
  $OptionsInfo{WriteUpperTriangularMatrix} = 0; $OptionsInfo{WriteLowerTriangularMatrix} = 0;
  OUTMATRIXTYPE: {
    if ($OptionsInfo{OutMatrixType} =~ /^FullMatrix$/i) {
      $OptionsInfo{WriteFullMatrix} = 1; last OUTMATRIXTYPE;
    }
    if ($OptionsInfo{OutMatrixType} =~ /^UpperTriangularMatrix$/i) {
      $OptionsInfo{WriteUpperTriangularMatrix} = 1; last OUTMATRIXTYPE;
    }
    if ($OptionsInfo{OutMatrixType} =~ /^LowerTriangularMatrix$/i) {
      $OptionsInfo{WriteLowerTriangularMatrix} = 1; last OUTMATRIXTYPE;
    }
    die "Error: The value specified, $Options{outmatrixtype}, for option \"--OutMatrixType\" is not valid. Allowed values: FullMatrix, UpperTriangularMatrix or LowerTriangularMatrix\n";
  }

  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? 1 : 0;
  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : 0;

  $OptionsInfo{Fast} = $Options{fast} ? 1 : 0;
  $OptionsInfo{ValidateData} = $Options{fast} ? 0 : 1;

  $OptionsInfo{Precision} = $Options{precision};

}

# Process options related to comparion of bit vector strings...
#
sub ProcessBitVectorComparisonOptions {
  # Setup supported bit vector similarity coefficients for bit vector strings...
  my($ComparisonMeasure, $SupportedComparisonMeasure, @SupportedComparisonMeasures, %SupportedComparisonMeasuresNameMap, %SupportedComparisonMeasuresMethodMap);

  @SupportedComparisonMeasures = ();
  %SupportedComparisonMeasuresNameMap = ();
  %SupportedComparisonMeasuresMethodMap = ();

  for $SupportedComparisonMeasure (Fingerprints::FingerprintsBitVector::GetSupportedSimilarityCoefficients()) {
    # Similarity coefficient function/method names contain "Coefficient" in their names.
    # So take 'em out and setup a map to original function/method name...
    $ComparisonMeasure = $SupportedComparisonMeasure;
    $ComparisonMeasure =~ s/Coefficient$//;

    push @SupportedComparisonMeasures, $ComparisonMeasure;
    $SupportedComparisonMeasuresNameMap{lc($ComparisonMeasure)} = $ComparisonMeasure;
    $SupportedComparisonMeasuresMethodMap{lc($ComparisonMeasure)} = $SupportedComparisonMeasure;
  }

  # Setup a list of similarity coefficients to use for calculating similarity matrices for bit vector strings...
  my($SpecifiedMeasure, @SpecifiedComparisonMeasures, %SpecifiedComparisonMeasuresNameMap, %SpecifiedComparisonMeasuresMethodMap, %SpecifiedComparisonMeasuresParameterMap);

  @SpecifiedComparisonMeasures = ();
  %SpecifiedComparisonMeasuresNameMap = ();
  %SpecifiedComparisonMeasuresMethodMap = ();
  %SpecifiedComparisonMeasuresParameterMap = ();

  if ($Options{bitvectorcomparisonmode} =~ /^All$/i) {
    push @SpecifiedComparisonMeasures, @SupportedComparisonMeasures;
  }
  else {
    # Comma delimited list of similarity coefficients...
    my($BitVectorComparisonMode, @SpecifiedMeasures, @UnsupportedSpecifiedMeasures);

    $BitVectorComparisonMode = $Options{bitvectorcomparisonmode};
    $BitVectorComparisonMode =~ s/ //g;
    @SpecifiedMeasures = split ",", $BitVectorComparisonMode;
    @UnsupportedSpecifiedMeasures = ();

    for $SpecifiedMeasure (@SpecifiedMeasures) {
      if (exists($SupportedComparisonMeasuresMethodMap{lc($SpecifiedMeasure)})) {
	push @SpecifiedComparisonMeasures, $SpecifiedMeasure;
      }
      else {
	push @UnsupportedSpecifiedMeasures, $SpecifiedMeasure;
      }
    }
    if (@UnsupportedSpecifiedMeasures) {
      if (@UnsupportedSpecifiedMeasures > 1) {
	warn "Error: The values specified - ", JoinWords(\@UnsupportedSpecifiedMeasures, ", ", 0)," - for option \"-b --BitVectorComparisonMode\" are not valid.\n";
      }
      else {
	warn "Error: The value specified, @UnsupportedSpecifiedMeasures, for option \"-b --BitVectorComparisonMode\" is not valid.\n";
      }
      die "Allowed values:", JoinWords(\@SupportedComparisonMeasures, ", ", 0), "\n";
    }
  }
  for $SpecifiedMeasure (@SpecifiedComparisonMeasures) {
    $SpecifiedComparisonMeasuresMethodMap{lc($SpecifiedMeasure)} = $SupportedComparisonMeasuresMethodMap{lc($SpecifiedMeasure)};
    $SpecifiedComparisonMeasuresNameMap{lc($SpecifiedMeasure)} = $SupportedComparisonMeasuresNameMap{lc($SpecifiedMeasure)};
  }

  $OptionsInfo{BitVectorComparisonMode} = $Options{bitvectorcomparisonmode};
  $OptionsInfo{SpecifiedBitVectorComparisonsRef} = \@SpecifiedComparisonMeasures;
  $OptionsInfo{SpecifiedBitVectorComparisonsNameRef} = \%SpecifiedComparisonMeasuresNameMap;
  $OptionsInfo{SpecifiedBitVectorComparisonsMethodRef} = \%SpecifiedComparisonMeasuresMethodMap;

  # Make sure valid alpha parameter is specified for Tversky calculation...
  my($SpecifiedMeasure1, $SpecifiedMeasure2);
  $OptionsInfo{Alpha} = '';
  $SpecifiedMeasure1 = 'TverskySimilarity';
  $SpecifiedMeasure2 = 'WeightedTverskySimilarity';
  if ($SpecifiedComparisonMeasuresNameMap{lc($SpecifiedMeasure1)} || $SpecifiedComparisonMeasuresNameMap{lc($SpecifiedMeasure2)}) {
    if (IsEmpty($Options{alpha})) {
      die "Error: You must specify a value for \"-a, --alpha\" option in \"$SpecifiedMeasure1, $SpecifiedMeasure2, or All\" \"-m --mode\". \n";
    }
    my($Alpha);
    $Alpha = $Options{alpha};
    if (!(IsFloat($Alpha) && $Alpha >=0 && $Alpha <= 1)) {
      die "Error: The value specified, $Options{alpha}, for option \"-a, --alpha\" is not valid. Allowed values: >= 0 and <= 1\n";
    }
    $OptionsInfo{Alpha} = $Alpha;
  }

  # Make sure valid beta parameter is specified for WeightedTanimoto and WeightedTversky
  # calculations...
  $OptionsInfo{Beta} = '';
  $SpecifiedMeasure1 = 'WeightedTverskySimilarity';
  $SpecifiedMeasure2 = 'WeightedTanimotoSimilarity';
  if ($SpecifiedComparisonMeasuresNameMap{lc($SpecifiedMeasure1)} || $SpecifiedComparisonMeasuresNameMap{lc($SpecifiedMeasure2)}) {
    if (IsEmpty($Options{beta})) {
      die "Error: You must specify a value for \"-b, --beta\" option in \"$SpecifiedMeasure1, $SpecifiedMeasure2, or All\" \"-m --mode\". \n";
    }
    my($Beta);
    $Beta = $Options{beta};
    if (!(IsFloat($Beta) && $Beta >=0 && $Beta <= 1)) {
      die "Error: The value specified, $Options{beta}, for option \"-b, --beta\" is not valid. Allowed values: >= 0 and <= 1\n";
    }
    $OptionsInfo{Beta} = $Beta;
  }

  # Setup any parameters required for specified comparison menthod...
  for $SpecifiedMeasure (@SpecifiedComparisonMeasures) {
    @{$SpecifiedComparisonMeasuresParameterMap{lc($SpecifiedMeasure)}} = ();
    if ($SpecifiedMeasure =~ /^TverskySimilarity$/i) {
      push @{$SpecifiedComparisonMeasuresParameterMap{lc($SpecifiedMeasure)}}, $OptionsInfo{Alpha};
    }
    elsif ($SpecifiedMeasure =~ /^WeightedTverskySimilarity$/i) {
      push @{$SpecifiedComparisonMeasuresParameterMap{lc($SpecifiedMeasure)}}, $OptionsInfo{Alpha};
      push @{$SpecifiedComparisonMeasuresParameterMap{lc($SpecifiedMeasure)}}, $OptionsInfo{Beta};
    }
    elsif ($SpecifiedMeasure =~ /^WeightedTanimotoSimilarity$/i) {
      push @{$SpecifiedComparisonMeasuresParameterMap{lc($SpecifiedMeasure)}}, $OptionsInfo{Beta};
    }
  }
  $OptionsInfo{SpecifiedBitVectorComparisonsParameterRef} = \%SpecifiedComparisonMeasuresParameterMap;
}

# Process options related to comparion of vector strings...
#
sub ProcessVectorComparisonOptions {
  # Setup specified similarity coefficients for vector strings..
  my($ComparisonMeasure, $SupportedComparisonMeasure, @SupportedComparisonMeasures, %SupportedComparisonMeasuresNameMap, %SupportedComparisonMeasuresMethodMap);

  @SupportedComparisonMeasures = ();
  %SupportedComparisonMeasuresNameMap = ();
  %SupportedComparisonMeasuresMethodMap = ();
  for $SupportedComparisonMeasure (Fingerprints::FingerprintsVector::GetSupportedDistanceAndSimilarityCoefficients()) {
    # Similarity and distance coefficient function/method names contain "Coefficient" in their names.
    # So take 'em out and setup a map to original function/method name...
    $ComparisonMeasure = $SupportedComparisonMeasure;
    if ($ComparisonMeasure =~ /Coefficient$/i) {
      $ComparisonMeasure =~ s/Coefficient$//i;
    }
    push @SupportedComparisonMeasures, $ComparisonMeasure;
    $SupportedComparisonMeasuresNameMap{lc($ComparisonMeasure)} = $ComparisonMeasure;
    $SupportedComparisonMeasuresMethodMap{lc($ComparisonMeasure)} = $SupportedComparisonMeasure;
  }

  # Setup a list of similarity coefficients to use for calculating similarity matrices for bit vector strings...
  my($SpecifiedMeasure, @SpecifiedComparisonMeasures, %SpecifiedComparisonMeasuresNameMap, %SpecifiedComparisonMeasuresMethodMap, %SpecifiedComparisonMeasuresParameterMap);

  @SpecifiedComparisonMeasures = ();
  %SpecifiedComparisonMeasuresNameMap = ();
  %SpecifiedComparisonMeasuresMethodMap = ();

  if ($Options{vectorcomparisonmode} =~ /^All$/i) {
    push @SpecifiedComparisonMeasures, @SupportedComparisonMeasures;
  }
  else {
    # Comma delimited list of similarity coefficients...
    my($VectorComparisonMode, @SpecifiedMeasures, @UnsupportedSpecifiedMeasures);

    $VectorComparisonMode = $Options{vectorcomparisonmode};
    $VectorComparisonMode =~ s/ //g;
    @SpecifiedMeasures = split ",", $VectorComparisonMode;
    @UnsupportedSpecifiedMeasures = ();

    for $SpecifiedMeasure (@SpecifiedMeasures) {
      if (exists($SupportedComparisonMeasuresMethodMap{lc($SpecifiedMeasure)})) {
	push @SpecifiedComparisonMeasures, $SpecifiedMeasure;
      }
      else {
	push @UnsupportedSpecifiedMeasures, $SpecifiedMeasure;
      }
    }
    if (@UnsupportedSpecifiedMeasures) {
      if (@UnsupportedSpecifiedMeasures > 1) {
	warn "Error: The values specified - ", JoinWords(\@UnsupportedSpecifiedMeasures, ", ", 0)," - for option \"-v --VectorComparisonMode\" are not valid.\n";
      }
      else {
	warn "Error: The value specified, @UnsupportedSpecifiedMeasures, for option \"-v --VectorComparisonMode\" is not valid.\n";
      }
      die "Allowed values:", JoinWords(\@SupportedComparisonMeasures, ", ", 0), "\n";
    }
  }
  for $SpecifiedMeasure (@SpecifiedComparisonMeasures) {
    $SpecifiedComparisonMeasuresMethodMap{lc($SpecifiedMeasure)} = $SupportedComparisonMeasuresMethodMap{lc($SpecifiedMeasure)};
    $SpecifiedComparisonMeasuresNameMap{lc($SpecifiedMeasure)} = $SupportedComparisonMeasuresNameMap{lc($SpecifiedMeasure)};
  }

  $OptionsInfo{VectorComparisonMode} = $Options{vectorcomparisonmode};
  $OptionsInfo{SpecifiedVectorComparisonsRef} = \@SpecifiedComparisonMeasures;
  $OptionsInfo{SpecifiedVectorComparisonsNameRef} = \%SpecifiedComparisonMeasuresNameMap;
  $OptionsInfo{SpecifiedVectorComparisonsMethodRef} = \%SpecifiedComparisonMeasuresMethodMap;

  # Setup specified vector comparison calculation modes...
  my(@SpecifiedVectorComparisonModes);
  @SpecifiedVectorComparisonModes = ();
  if ($Options{vectorcomparisonformulism} =~ /^All$/i) {
    push @SpecifiedVectorComparisonModes, ("AlgebraicForm", "BinaryForm", "SetTheoreticForm");
  }
  else {
    my($SpecifiedFormulism, @SpecifiedFormulismWords);

    @SpecifiedFormulismWords = split /\,/, $Options{vectorcomparisonformulism};
    for $SpecifiedFormulism (@SpecifiedFormulismWords) {
      if ($SpecifiedFormulism !~ /^(AlgebraicForm|BinaryForm|SetTheoreticForm)$/i) {
	die "Error: The value specified, $SpecifiedFormulism, for option \"--VectorComparisonFormulism\" is not valid. Allowed values: AlgebraicForm, BinaryForm or SetTheoreticForm\n";
      }
      push @SpecifiedVectorComparisonModes, $SpecifiedFormulism;
    }
  }
  $OptionsInfo{VectorComparisonFormulism} = $Options{vectorcomparisonformulism};
  $OptionsInfo{SpecifiedVectorComparisonModesRef} = \@SpecifiedVectorComparisonModes;

  # Setup any parameters required for specified comparison menthod...
  for $SpecifiedMeasure (@SpecifiedComparisonMeasures) {
    @{$SpecifiedComparisonMeasuresParameterMap{lc($SpecifiedMeasure)}} = ();
    push @{$SpecifiedComparisonMeasuresParameterMap{lc($SpecifiedMeasure)}}, ($Options{fast} ? 1 : 0);
  }
  $OptionsInfo{SpecifiedVectorComparisonsParameterRef} = \%SpecifiedComparisonMeasuresParameterMap;
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{alpha} = 0.5;
  $Options{beta} = 1;

  $Options{bitvectorcomparisonmode} = "TanimotoSimilarity";

  $Options{colmode} = 'colnum';

  $Options{compoundidprefix} = 'Cmpd';
  $Options{compoundidmode} = 'LabelPrefix';

  $Options{detail} = 1;

  $Options{indelim} = 'comma';
  $Options{outdelim} = 'comma';

  $Options{inputdatamode} = 'LoadInMemory';

  $Options{mode} = 'AutoDetect';

  $Options{outmatrixformat} = 'RowsAndColumns';

  $Options{outmatrixtype} = 'FullMatrix';

  $Options{quote} = 'yes';
  $Options{precision} = 2;

  $Options{vectorcomparisonmode} = "TanimotoSimilarity";
  $Options{vectorcomparisonformulism} = "AlgebraicForm";

  if (!GetOptions(\%Options, "alpha=f", "beta=f", "bitvectorcomparisonmode|b=s", "colmode|c=s", "compoundidcol=s", "compoundidprefix=s", "compoundidfield=s", "compoundidmode=s", "detail|d=i", "fast|f", "fingerprintscol=s", "fingerprintsfield=s", "help|h", "indelim=s", "inputdatamode=s", "mode|m=s", "outdelim=s", "overwrite|o", "outmatrixformat=s", "outmatrixtype=s", "precision|p=s", "quote|q=s", "root|r=s", "vectorcomparisonmode|v=s", "vectorcomparisonformulism=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{colmode} !~ /^(ColNum|ColLabel)$/i) {
    die "Error: The value specified, $Options{colmode}, for option \"-c, --ColMode\" is not valid. Allowed values: ColNum, or ColLabel\n";
  }
  if ($Options{compoundidmode} !~ /^(DataField|MolName|LabelPrefix|MolNameOrLabelPrefix)$/i) {
    die "Error: The value specified, $Options{compoundidmode}, for option \"--CompoundIDMode\" is not valid. Allowed values: DataField, MolName, LabelPrefix or MolNameOrLabelPrefix\n";
  }
  if (!IsPositiveInteger($Options{detail})) {
    die "Error: The value specified, $Options{detail}, for option \"-d, --detail\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{inputdatamode} !~ /^(LoadInMemory|ScanFile)$/i) {
    die "Error: The value specified, $Options{inputdatamode}, for option \"--InputDataMode\" is not valid. Allowed values: LoadInMemory or ScanFile\n";
  }
  if ($Options{mode} !~ /^(AutoDetect|FingerprintsBitVectorString|FingerprintsVectorString)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: AutoDetect, FingerprintsBitVectorString or FingerprintsVectorString \n";
  }
  if ($Options{indelim} !~ /^(comma|semicolon)$/i) {
    die "Error: The value specified, $Options{indelim}, for option \"--InDelim\" is not valid. Allowed values: comma, or semicolon\n";
  }
  if ($Options{outdelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{outdelim}, for option \"--OutDelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if ($Options{outmatrixformat} !~ /^(RowsAndColumns|IDPairsAndValue)$/i) {
    die "Error: The value specified, $Options{outmatrixformat}, for option \"--OutMatrixFormat\" is not valid. Allowed values: RowsAndColumns or IDPairsAndValue\n";
  }
  if ($Options{outmatrixtype} !~ /^(FullMatrix|UpperTriangularMatrix|LowerTriangularMatrix)$/i) {
    die "Error: The value specified, $Options{outmatrixtype}, for option \"--OutMatrixType\" is not valid. Allowed values: FullMatrix, UpperTriangularMatrix or LowerTriangularMatrix\n";
  }
  if ($Options{quote} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not valid. Allowed values: Yes or No\n";
  }
  if (!IsPositiveInteger($Options{precision})) {
    die "Error: The value specified, $Options{precision}, for option \"--precision\" is not valid. Allowed values: > 0 \n";
  }
}

__END__

=head1 NAME

SimilarityMatricesFingerprints.pl - Calculate similarity matrices using fingerprints strings data in SD, FP and CSV/TSV text file(s)

=head1 SYNOPSIS

SimilarityMatricesFingerprints.pl SDFile(s) FPFile(s) TextFile(s)...

SimilarityMatricesFingerprints.pl [B<--alpha> I<number>] [B<--beta> I<number>]
[B<-b, --BitVectorComparisonMode> I<All | "TanimotoSimilarity,[ TverskySimilarity, ... ]">]
[B<-c, --ColMode> I<ColNum | ColLabel>] [B<--CompoundIDCol> I<col number | col name>]
[B<--CompoundIDPrefix> I<text>] [B<--CompoundIDField> I<DataFieldName>]
[B<--CompoundIDMode> I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>]
[B<-d, --detail> I<InfoLevel>] [B<-f, --fast>] [B<--FingerprintsCol> I<col number | col name>]
[B<--FingerprintsField> I<FieldLabel>] [B<-h, --help>]  [B<--InDelim> I<comma | semicolon>]
[B<--InputDataMode> I<LoadInMemory | ScanFile>]
[B<-m, --mode> I<AutoDetect | FingerprintsBitVectorString | FingerprintsVectorString>]
[B<--OutDelim> I<comma | tab | semicolon>] [B<--OutMatrixFormat> I<RowsAndColumns | IDPairsAndValue>]
[B<--OutMatrixType> I<FullMatrix | UpperTriangularMatrix | LowerTriangularMatrix>]
[B<-o, --overwrite>] [B<-p, --precision> I<number>]
[B<-q, --quote> I<Yes | No>] [B<-r, --root> I<RootName>]
[B<-v, --VectorComparisonMode> I<All | "TanimotoSimilairy, [ ManhattanDistance, ...]">]
[B<--VectorComparisonFormulism> I<All | "AlgebraicForm, [BinaryForm, SetTheoreticForm]">]
[B<-w, --WorkingDir> dirname] SDFile(s) FPFile(s) TextFile(s)...

=head1 DESCRIPTION

Calculate similarity matrices using fingerprint bit-vector or vector strings data in I<SD, FP
and CSV/TSV> text file(s) and generate CSV/TSV text file(s) containing values for specified
similarity and distance coefficients.

The scripts SimilarityMatrixSDFiles.pl and SimilarityMatrixTextFiles.pl have been removed from the
current release of MayaChemTools and their functionality merged with this script.

The valid I<SDFile> extensions are I<.sdf> and I<.sd>. All SD files in a current directory
can be specified either by I<*.sdf> or the current directory name.

The valid I<FPFile> extensions are I<.fpf> and I<.fp>. All FP files in a current directory
can be specified either by I<*.fpf> or the current directory name.

The valid I<TextFile> extensions are I<.csv> and I<.tsv> for comma/semicolon and tab
delimited text files respectively. All other file names are ignored. All text files in a
current directory can be specified by I<*.csv>, I<*.tsv>, or the current directory
name. The B<--indelim> option determines the format of I<TextFile(s)>. Any file
which doesn't correspond to the format indicated by B<--indelim> option is ignored.

Example of I<FP> file containing fingerprints bit-vector string data:

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

Example of I<FP> file containing fingerprints vector string data:

    #
    # Package = MayaChemTools 7.4
    # ReleaseDate = Oct 21, 2010
    #
    # TimeStamp =  Mon Mar 7 15:14:01 2011
    #
    # FingerprintsStringType = FingerprintsVector
    #
    # Description = PathLengthBits:AtomicInvariantsAtomTypes:MinLength1:...
    # VectorStringFormat = IDsAndValuesString
    # VectorValuesType = NumericalValues
    #
    Cmpd1 338;C F N O C:C C:N C=O CC CF CN CO C:C:C C:C:N C:CC C:CF C:CN C:
    N:C C:NC CC:N CC=O CCC CCN CCO CNC NC=O O=CO C:C:C:C C:C:C:N C:C:CC...;
    33 1 2 5 21 2 2 12 1 3 3 20 2 10 2 2 1 2 2 2 8 2 5 1 1 1 19 2 8 2 2 2 2
    6 2 2 2 2 2 2 2 2 3 2 2 1 4 1 5 1 1 18 6 2 2 1 2 10 2 1 2 1 2 2 2 2 ...
    Cmpd2 103;C N O C=N C=O CC CN CO CC=O CCC CCN CCO CNC N=CN NC=O NCN O=C
    O C CC=O CCCC CCCN CCCO CCNC CNC=N CNC=O CNCN CCCC=O CCCCC CCCCN CC...;
    15 4 4 1 2 13 5 2 2 15 5 3 2 2 1 1 1 2 17 7 6 5 1 1 1 2 15 8 5 7 2 2 2 2
    1 2 1 1 3 15 7 6 8 3 4 4 3 2 2 1 2 3 14 2 4 7 4 4 4 4 1 1 1 2 1 1 1 ...
    ... ...
    ... ...

Example of I<SD> file containing fingerprints bit-vector string data:

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

Example of CSV I<Text> file containing fingerprints bit-vector string data:

    "CompoundID","PathLengthFingerprints"
    "Cmpd1","FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes
    :MinLength1:MaxLength8;1024;HexadecimalString;Ascending;9c8460989ec8a4
    9913991a6603130b0a19e8051c89184414953800cc2151082844a20104280013086030
    8e8204d402800831048940e44281c00060449a5000ac80c894114e006321264401..."
    ... ...
    ... ...

The current release of MayaChemTools supports the following types of fingerprint
bit-vector and vector strings:

    FingerprintsVector;AtomNeighborhoods:AtomicInvariantsAtomTypes:MinRadi
    us0:MaxRadius2;41;AlphaNumericalValues;ValuesString;NR0-C.X1.BO1.H3-AT
    C1:NR1-C.X3.BO3.H1-ATC1:NR2-C.X1.BO1.H3-ATC1:NR2-C.X3.BO4-ATC1 NR0-C.X
    1.BO1.H3-ATC1:NR1-C.X3.BO3.H1-ATC1:NR2-C.X1.BO1.H3-ATC1:NR2-C.X3.BO4-A
    TC1 NR0-C.X2.BO2.H2-ATC1:NR1-C.X2.BO2.H2-ATC1:NR1-C.X3.BO3.H1-ATC1:NR2
    -C.X2.BO2.H2-ATC1:NR2-N.X3.BO3-ATC1:NR2-O.X1.BO1.H1-ATC1 NR0-C.X2.B...

    FingerprintsVector;AtomTypesCount:AtomicInvariantsAtomTypes:ArbitraryS
    ize;10;NumericalValues;IDsAndValuesString;C.X1.BO1.H3 C.X2.BO2.H2 C.X2
    .BO3.H1 C.X3.BO3.H1 C.X3.BO4 F.X1.BO1 N.X2.BO2.H1 N.X3.BO3 O.X1.BO1.H1
    O.X1.BO2;2 4 14 3 10 1 1 1 3 2

    FingerprintsVector;AtomTypesCount:SLogPAtomTypes:ArbitrarySize;16;Nume
    ricalValues;IDsAndValuesString;C1 C10 C11 C14 C18 C20 C21 C22 C5 CS F
    N11 N4 O10 O2 O9;5 1 1 1 14 4 2 1 2 2 1 1 1 1 3 1

    FingerprintsVector;AtomTypesCount:SLogPAtomTypes:FixedSize;67;OrderedN
    umericalValues;IDsAndValuesString;C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C
    12 C13 C14 C15 C16 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C27 CS N1 N
    2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12 N13 N14 NS O1 O2 O3 O4 O5 O6 O7 O8
    O9 O10 O11 O12 OS F Cl Br I Hal P S1 S2 S3 Me1 Me2;5 0 0 0 2 0 0 0 0 1
    1 0 0 1 0 0 0 14 0 4 2 1 0 0 0 0 0 2 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0...

    FingerprintsVector;EStateIndicies:ArbitrarySize;11;NumericalValues;IDs
    AndValuesString;SaaCH SaasC SaasN SdO SdssC SsCH3 SsF SsOH SssCH2 SssN
    H SsssCH;24.778 4.387 1.993 25.023 -1.435 3.975 14.006 29.759 -0.073 3
    .024 -2.270

    FingerprintsVector;EStateIndicies:FixedSize;87;OrderedNumericalValues;
    ValuesString;0 0 0 0 0 0 0 3.975 0 -0.073 0 0 24.778 -2.270 0 0 -1.435
    4.387 0 0 0 0 0 0 3.024 0 0 0 0 0 0 0 1.993 0 29.759 25.023 0 0 0 0 1
    4.006 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0

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

    FingerprintsVector;ExtendedConnectivity:FunctionalClassAtomTypes:Radiu
    s2;57;AlphaNumericalValues;ValuesString;24769214 508787397 850393286 8
    62102353 981185303 1231636850 1649386610 1941540674 263599683 32920567
    1 571109041 639579325 683993318 723853089 810600886 885767127 90326012
    7 958841485 981022393 1126908698 1152248391 1317567065 1421489994 1455
    632544 1557272891 1826413669 1983319256 2015750777 2029559552 20404...

    FingerprintsVector;ExtendedConnectivity:EStateAtomTypes:Radius2;62;Alp
    haNumericalValues;ValuesString;25189973 528584866 662581668 671034184
    926543080 1347067490 1738510057 1759600920 2034425745 2097234755 21450
    44754 96779665 180364292 341712110 345278822 386540408 387387308 50430
    1706 617094135 771528807 957666640 997798220 1158349170 1291258082 134
    1138533 1395329837 1420277211 1479584608 1486476397 1487556246 1566...

    FingerprintsBitVector;MACCSKeyBits;166;BinaryString;Ascending;00000000
    0000000000000000000000000000000001001000010010000000010010000000011100
    0100101010111100011011000100110110000011011110100110111111111111011111
    11111111111110111000

    FingerprintsBitVector;MACCSKeyBits;322;BinaryString;Ascending;11101011
    1110011111100101111111000111101100110000000000000011100010000000000000
    0000000000000000000000000000000000000000000000101000000000000000000000
    0000000000000000000000000000000000000000000000000000000000000000000000
    0000000000000000000000000000000000000011000000000000000000000000000000
    0000000000000000000000000000000000000000

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

    FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes:MinLeng
    th1:MaxLength8;1024;BinaryString;Ascending;001000010011010101011000110
    0100010101011000101001011100110001000010001001101000001001001001001000
    0010110100000111001001000001001010100100100000000011000000101001011100
    0010000001000101010100000100111100110111011011011000000010110111001101
    0101100011000000010001000011000010100011101100001000001000100000000...

    FingerprintsVector;PathLengthCount:AtomicInvariantsAtomTypes:MinLength
    1:MaxLength8;432;NumericalValues;IDsAndValuesPairsString;C.X1.BO1.H3 2
    C.X2.BO2.H2 4 C.X2.BO3.H1 14 C.X3.BO3.H1 3 C.X3.BO4 10 F.X1.BO1 1 N.X
    2.BO2.H1 1 N.X3.BO3 1 O.X1.BO1.H1 3 O.X1.BO2 2 C.X1.BO1.H3C.X3.BO3.H1
    2 C.X2.BO2.H2C.X2.BO2.H2 1 C.X2.BO2.H2C.X3.BO3.H1 4 C.X2.BO2.H2C.X3.BO
    4 1 C.X2.BO2.H2N.X3.BO3 1 C.X2.BO3.H1:C.X2.BO3.H1 10 C.X2.BO3.H1:C....

    FingerprintsVector;PathLengthCount:MMFF94AtomTypes:MinLength1:MaxLengt
    h8;463;NumericalValues;IDsAndValuesPairsString;C5A 2 C5B 2 C=ON 1 CB 1
    8 COO 1 CR 9 F 1 N5 1 NC=O 1 O=CN 1 O=CO 1 OC=O 1 OR 2 C5A:C5B 2 C5A:N
    5 2 C5ACB 1 C5ACR 1 C5B:C5B 1 C5BC=ON 1 C5BCB 1 C=ON=O=CN 1 C=ONNC=O 1
    CB:CB 18 CBF 1 CBNC=O 1 COO=O=CO 1 COOCR 1 COOOC=O 1 CRCR 7 CRN5 1 CR
    OR 2 C5A:C5B:C5B 2 C5A:C5BC=ON 1 C5A:C5BCB 1 C5A:N5:C5A 1 C5A:N5CR ...

    FingerprintsVector;TopologicalAtomPairs:AtomicInvariantsAtomTypes:MinD
    istance1:MaxDistance10;223;NumericalValues;IDsAndValuesString;C.X1.BO1
    .H3-D1-C.X3.BO3.H1 C.X2.BO2.H2-D1-C.X2.BO2.H2 C.X2.BO2.H2-D1-C.X3.BO3.
    H1 C.X2.BO2.H2-D1-C.X3.BO4 C.X2.BO2.H2-D1-N.X3.BO3 C.X2.BO3.H1-D1-...;
    2 1 4 1 1 10 8 1 2 6 1 2 2 1 2 1 2 2 1 2 1 5 1 10 12 2 2 1 2 1 9 1 3 1
    1 1 2 2 1 3 6 1 6 14 2 2 2 3 1 3 1 8 2 2 1 3 2 6 1 2 2 5 1 3 1 23 1...

    FingerprintsVector;TopologicalAtomPairs:FunctionalClassAtomTypes:MinDi
    stance1:MaxDistance10;144;NumericalValues;IDsAndValuesString;Ar-D1-Ar
    Ar-D1-Ar.HBA Ar-D1-HBD Ar-D1-Hal Ar-D1-None Ar.HBA-D1-None HBA-D1-NI H
    BA-D1-None HBA.HBD-D1-NI HBA.HBD-D1-None HBD-D1-None NI-D1-None No...;
    23 2 1 1 2 1 1 1 1 2 1 1 7 28 3 1 3 2 8 2 1 1 1 5 1 5 24 3 3 4 2 13 4
    1 1 4 1 5 22 4 4 3 1 19 1 1 1 1 1 2 2 3 1 1 8 25 4 5 2 3 1 26 1 4 1 ...

    FingerprintsVector;TopologicalAtomTorsions:AtomicInvariantsAtomTypes;3
    3;NumericalValues;IDsAndValuesString;C.X1.BO1.H3-C.X3.BO3.H1-C.X3.BO4-
    C.X3.BO4 C.X1.BO1.H3-C.X3.BO3.H1-C.X3.BO4-N.X3.BO3 C.X2.BO2.H2-C.X2.BO
    2.H2-C.X3.BO3.H1-C.X2.BO2.H2 C.X2.BO2.H2-C.X2.BO2.H2-C.X3.BO3.H1-O...;
    2 2 1 1 2 2 1 1 3 4 4 8 4 2 2 6 2 2 1 2 1 1 2 1 1 2 6 2 4 2 1 3 1

    FingerprintsVector;TopologicalAtomTorsions:EStateAtomTypes;36;Numerica
    lValues;IDsAndValuesString;aaCH-aaCH-aaCH-aaCH aaCH-aaCH-aaCH-aasC aaC
    H-aaCH-aasC-aaCH aaCH-aaCH-aasC-aasC aaCH-aaCH-aasC-sF aaCH-aaCH-aasC-
    ssNH aaCH-aasC-aasC-aasC aaCH-aasC-aasC-aasN aaCH-aasC-ssNH-dssC a...;
    4 4 8 4 2 2 6 2 2 2 4 3 2 1 3 3 2 2 2 1 2 1 1 1 2 1 1 1 1 1 1 1 2 1 1 2

    FingerprintsVector;TopologicalAtomTriplets:AtomicInvariantsAtomTypes:M
    inDistance1:MaxDistance10;3096;NumericalValues;IDsAndValuesString;C.X1
    .BO1.H3-D1-C.X1.BO1.H3-D1-C.X3.BO3.H1-D2 C.X1.BO1.H3-D1-C.X2.BO2.H2-D1
    0-C.X3.BO4-D9 C.X1.BO1.H3-D1-C.X2.BO2.H2-D3-N.X3.BO3-D4 C.X1.BO1.H3-D1
    -C.X2.BO2.H2-D4-C.X2.BO2.H2-D5 C.X1.BO1.H3-D1-C.X2.BO2.H2-D6-C.X3....;
    1 2 2 2 2 2 2 2 8 8 4 8 4 4 2 2 2 2 4 2 2 2 4 2 2 2 2 1 2 2 4 4 4 2 2
    2 4 4 4 8 4 4 2 4 4 4 2 4 4 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 8...

    FingerprintsVector;TopologicalAtomTriplets:SYBYLAtomTypes:MinDistance1
    :MaxDistance10;2332;NumericalValues;IDsAndValuesString;C.2-D1-C.2-D9-C
    .3-D10 C.2-D1-C.2-D9-C.ar-D10 C.2-D1-C.3-D1-C.3-D2 C.2-D1-C.3-D10-C.3-
    D9 C.2-D1-C.3-D2-C.3-D3 C.2-D1-C.3-D2-C.ar-D3 C.2-D1-C.3-D3-C.3-D4 C.2
    -D1-C.3-D3-N.ar-D4 C.2-D1-C.3-D3-O.3-D2 C.2-D1-C.3-D4-C.3-D5 C.2-D1-C.
    3-D5-C.3-D6 C.2-D1-C.3-D5-O.3-D4 C.2-D1-C.3-D6-C.3-D7 C.2-D1-C.3-D7...

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

=head1 OPTIONS

=over 4

=item B<--alpha> I<number>

Value of alpha parameter for calculating I<Tversky> similarity coefficient specified for
B<-b, --BitVectorComparisonMode> option. It corresponds to weights assigned for bits set
to "1" in a pair of fingerprint bit-vectors during the calculation of similarity coefficient. Possible
values: I<0 to 1>. Default value: <0.5>.

=item B<--beta> I<number>

Value of beta parameter for calculating I<WeightedTanimoto> and  I<WeightedTversky>
similarity coefficients specified for B<-b, --BitVectorComparisonMode> option. It is used to
weight the contributions of bits set to "0" during the calculation of similarity coefficients. Possible
values: I<0 to 1>. Default value of <1> makes I<WeightedTanimoto> and  I<WeightedTversky>
equivalent to I<Tanimoto> and  I<Tversky>.

=item B<-b, --BitVectorComparisonMode> I<All | "TanimotoSimilarity,[TverskySimilarity,...]">

Specify what similarity coefficients to use for calculating similarity matrices for fingerprints bit-vector
strings data values in I<TextFile(s)>: calculate similarity matrices for all supported similarity
coefficients or specify a comma delimited list of similarity coefficients. Possible values:
I<All | "TanimotoSimilarity,[TverskySimilarity,...]>. Default: I<TanimotoSimilarity>

I<All> uses complete list of supported similarity coefficients: I<BaroniUrbaniSimilarity, BuserSimilarity,
CosineSimilarity, DiceSimilarity, DennisSimilarity, ForbesSimilarity, FossumSimilarity, HamannSimilarity, JacardSimilarity,
Kulczynski1Similarity, Kulczynski2Similarity, MatchingSimilarity, McConnaugheySimilarity, OchiaiSimilarity,
PearsonSimilarity, RogersTanimotoSimilarity, RussellRaoSimilarity, SimpsonSimilarity, SkoalSneath1Similarity,
SkoalSneath2Similarity, SkoalSneath3Similarity, TanimotoSimilarity, TverskySimilarity, YuleSimilarity,
WeightedTanimotoSimilarity, WeightedTverskySimilarity>. These similarity coefficients are described below.

For two fingerprint bit-vectors A and B of same size, let:

    Na = Number of bits set to "1" in A
    Nb = Number of bits set to "1" in B
    Nc = Number of bits set to "1" in both A and B
    Nd = Number of bits set to "0" in both A and B

    Nt = Number of bits set to "1" or "0" in A or B (Size of A or B)
    Nt = Na + Nb - Nc + Nd

    Na - Nc = Number of bits set to "1" in A but not in B
    Nb - Nc = Number of bits set to "1" in B but not in A

Then, various similarity coefficients [ Ref. 40 - 42 ] for a pair of bit-vectors A and B are
defined as follows:

I<BaroniUrbaniSimilarity>: ( SQRT( Nc * Nd ) + Nc ) / (  SQRT ( Nc * Nd ) + Nc + ( Na - Nc )  + ( Nb - Nc ) ) ( same as Buser )

I<BuserSimilarity>: ( SQRT ( Nc * Nd ) + Nc ) / (  SQRT ( Nc * Nd ) + Nc + ( Na - Nc )  + ( Nb - Nc ) ) ( same as BaroniUrbani )

I<CosineSimilarity>: Nc / SQRT ( Na * Nb ) (same as Ochiai)

I<DiceSimilarity>: (2 * Nc) / ( Na + Nb )

I<DennisSimilarity>: ( Nc * Nd - ( ( Na - Nc ) * ( Nb - Nc ) ) ) / SQRT ( Nt * Na * Nb)

I<ForbesSimilarity>: ( Nt * Nc ) / ( Na * Nb )

I<FossumSimilarity>: ( Nt * ( ( Nc - 1/2 ) ** 2 ) / ( Na * Nb )

I<HamannSimilarity>: ( ( Nc + Nd ) - ( Na - Nc ) - ( Nb - Nc ) ) / Nt

I<JaccardSimilarity>: Nc /  ( ( Na - Nc) + ( Nb - Nc ) + Nc ) = Nc / ( Na + Nb - Nc ) (same as Tanimoto)

I<Kulczynski1Similarity>: Nc / ( ( Na - Nc ) + ( Nb - Nc) ) = Nc / ( Na + Nb - 2Nc )

I<Kulczynski2Similarity>: ( ( Nc / 2 ) * ( 2 * Nc + ( Na - Nc ) + ( Nb - Nc) ) ) / ( ( Nc + ( Na - Nc ) ) * ( Nc + ( Nb - Nc ) ) ) = 0.5 * ( Nc / Na + Nc / Nb )

I<MatchingSimilarity>: ( Nc + Nd ) / Nt

I<McConnaugheySimilarity>: ( Nc ** 2 - ( Na - Nc ) * ( Nb - Nc) ) / (  Na * Nb )

I<OchiaiSimilarity>: Nc / SQRT ( Na * Nb ) (same as Cosine)

I<PearsonSimilarity>: ( ( Nc * Nd ) - ( ( Na - Nc ) * ( Nb - Nc ) ) / SQRT ( Na * Nb * (  Na - Nc + Nd ) * ( Nb - Nc + Nd ) )

I<RogersTanimotoSimilarity>: ( Nc + Nd ) / ( ( Na - Nc)  + ( Nb  - Nc) + Nt) = ( Nc + Nd ) / ( Na  + Nb  - 2Nc + Nt)

I<RussellRaoSimilarity>: Nc / Nt

I<SimpsonSimilarity>: Nc / MIN ( Na, Nb)

I<SkoalSneath1Similarity>: Nc / ( Nc + 2 * ( Na - Nc)  + 2 * ( Nb - Nc) ) = Nc / ( 2 * Na + 2 * Nb - 3 * Nc )

I<SkoalSneath2Similarity>: ( 2 * Nc + 2 * Nd ) / ( Nc + Nd + Nt )

I<SkoalSneath3Similarity>: ( Nc + Nd ) / ( ( Na - Nc ) + ( Nb - Nc ) ) = ( Nc + Nd ) / ( Na + Nb - 2 * Nc  )

I<TanimotoSimilarity>: Nc /  ( ( Na - Nc) + ( Nb - Nc ) + Nc ) = Nc / ( Na + Nb - Nc ) (same as Jaccard)

I<TverskySimilarity>: Nc / ( alpha * ( Na - Nc ) + ( 1 - alpha) * ( Nb - Nc) + Nc ) = Nc / ( alpha * ( Na - Nb )  + Nb)

I<YuleSimilarity>: ( ( Nc * Nd ) - ( ( Na - Nc ) * ( Nb - Nc ) ) ) / ( ( Nc * Nd ) + ( ( Na - Nc ) * ( Nb - Nc ) )  )

Values of Tanimoto/Jaccard and Tversky coefficients are dependent on only those bit which
are set to "1" in both A and B. In order to take into account all bit positions, modified versions
of Tanimoto [ Ref. 42 ] and Tversky [  Ref. 43 ] have been developed.

Let:

    Na' = Number of bits set to "0" in A
    Nb' = Number of bits set to "0" in B
    Nc' = Number of bits set to "0" in both A and B

Tanimoto': Nc' /  ( ( Na' - Nc') + ( Nb' - Nc' ) + Nc' ) = Nc' / ( Na' + Nb' - Nc' )

Tversky': Nc' / ( alpha * ( Na' - Nc' ) + ( 1 - alpha) * ( Nb' - Nc' ) + Nc' ) = Nc' / ( alpha * ( Na' - Nb' )  + Nb')

Then:

I<WeightedTanimotoSimilarity> = beta * Tanimoto + (1 - beta) * Tanimoto'

I<WeightedTverskySimilarity> = beta * Tversky + (1 - beta) * Tversky'

=item B<-c, --ColMode> I<ColNum | ColLabel>

Specify how columns are identified in I<TextFile(s)>: using column number or column
label. Possible values: I<ColNum or ColLabel>. Default value: I<ColNum>.

=item B<--CompoundIDCol> I<col number | col name>

This value is B<-c, --ColMode> mode specific. It specifies input I<TextFile(s)> column to use for
generating compound ID for similarity matrices  in output I<TextFile(s)>. Possible values: I<col number
or col label>. Default value: I<first column containing the word compoundID in its column label or sequentially
generated IDs>.

=item B<--CompoundIDPrefix> I<text>

Specify compound ID prefix to use during sequential generation of compound IDs for input I<SDFile(s)>
and I<TextFile(s)>. Default value: I<Cmpd>. The default value generates compound IDs which look
like Cmpd<Number>.

For input I<SDFile(s)>, this value is only used during I<LabelPrefix | MolNameOrLabelPrefix> values
of B<--CompoundIDMode> option; otherwise, it's ignored.

Examples for I<LabelPrefix> or I<MolNameOrLabelPrefix> value of B<--CompoundIDMode>:

    Compound

The values specified above generates compound IDs which correspond to Compound<Number>
instead of default value of Cmpd<Number>.

=item B<--CompoundIDField> I<DataFieldName>

Specify input I<SDFile(s)> datafield label for generating compound IDs. This value is only used
during I<DataField> value of B<--CompoundIDMode> option.

Examples for I<DataField> value of B<--CompoundIDMode>:

    MolID
    ExtReg

=item B<--CompoundIDMode> I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>

Specify how to generate compound IDs from input I<SDFile(s)> for similarity matrix CSV/TSV text
file(s): use a I<SDFile(s)> datafield value; use molname line from I<SDFile(s)>; generate a sequential ID
with specific prefix; use combination of both MolName and LabelPrefix with usage of LabelPrefix values
for empty molname lines.

Possible values: I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>.
Default: I<LabelPrefix>.

For I<MolNameAndLabelPrefix> value of B<--CompoundIDMode>, molname line in I<SDFile(s)> takes
precedence over sequential compound IDs generated using I<LabelPrefix> and only empty molname
values are replaced with sequential compound IDs.

=item B<-d, --detail> I<InfoLevel>

Level of information to print about lines being ignored. Default: I<1>. Possible values:
I<1, 2 or 3>.

=item B<-f, --fast>

In this mode, fingerprints columns specified using B<--FingerprintsCol> for I<TextFile(s)> and
B<--FingerprintsField> for I<SDFile(s)> are assumed to contain valid fingerprints data and no
checking is performed before calculating similarity matrices. By default, fingerprints data is
validated before computing pairwise similarity and distance coefficients.

=item B<--FingerprintsCol> I<col number | col name>

This value is B<-c, --colmode> specific. It specifies fingerprints column to use during
calculation similarity matrices for I<TextFile(s)>. Possible values: I<col number or col label>.
Default value: I<first column containing the word Fingerprints in its column label>.

=item B<--FingerprintsField> I<FieldLabel>

Fingerprints field label to use during calculation similarity matrices for I<SDFile(s)>.
Default value: I<first data field label containing the word Fingerprints in its label>

=item B<-h, --help>

Print this help message.

=item B<--InDelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile(s)>. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a
delimiter.

=item B<--InputDataMode> I<LoadInMemory | ScanFile>

Specify how fingerprints bit-vector or vector strings data from I<SD, FP and CSV/TSV>
fingerprint file(s) is processed: Retrieve, process and load all available fingerprints
data in memory; Retrieve and process data for fingerprints one at a time. Possible values
: I<LoadInMemory | ScanFile>. Default: I<LoadInMemory>.

During I<LoadInMemory> value of B<--InputDataMode>, fingerprints bit-vector or vector
strings data from input file is retrieved, processed, and loaded into memory all at once
as fingerprints objects for generation for similarity matrices.

During I<ScanFile> value of  B<--InputDataMode>, multiple passes over the input fingerprints
file are performed to retrieve and process fingerprints bit-vector or vector strings data one at
a time to generate fingerprints objects used during generation of similarity matrices. A temporary
copy of the input fingerprints file is made at the start and deleted after generating the matrices.

I<ScanFile> value of  B<--InputDataMode> allows processing of arbitrary large fingerprints files
without any additional memory requirement.

=item B<-m, --mode> I<AutoDetect | FingerprintsBitVectorString | FingerprintsVectorString>

Format of fingerprint strings data in  I<TextFile(s)>: automatically detect format of fingerprints
string created by MayaChemTools fingerprints generation scripts or explicitly specify its format.
Possible values: I<AutoDetect | FingerprintsBitVectorString | FingerprintsVectorString>. Default
value: I<AutoDetect>.

=item B<--OutDelim> I<comma | tab | semicolon>

Delimiter for output CSV/TSV text file(s). Possible values: I<comma, tab, or semicolon>
Default value: I<comma>.

=item B<--OutMatrixFormat> I<RowsAndColumns | IDPairsAndValue>

Specify how similarity or distance values calculated for fingerprints vector and bit-vector strings
are written to the output CSV/TSV text file(s): Generate text files containing rows and columns
with their labels corresponding to compound IDs and each matrix element value corresponding to
similarity or distance between corresponding compounds; Generate text files containing rows containing
compoundIDs for two compounds followed by similarity or distance value between these compounds.

Possible values: I<RowsAndColumns, or IDPairsAndValue>. Default value: I<RowsAndColumns>.

The value of B<--OutMatrixFormat> in conjunction with B<--OutMatrixType> determines type
of data written to output files and allows generation of up to 6 different output data formats:

    OutMatrixFormat OutMatrixType

    RowsAndColumns  FullMatrix   [ DEFAULT ]
    RowsAndColumns  UpperTriangularMatrix
    RowsAndColumns  LowerTriangularMatrix

    IDPairsAndValue FullMatrix
    IDPairsAndValue UpperTriangularMatrix
    IDPairsAndValue LowerTriangularMatrix

Example of data in output file for I<RowsAndColumns> B<--OutMatrixFormat> value for
I<FullMatrix> valueof B<--OutMatrixType>:

    "","Cmpd1","Cmpd2","Cmpd3","Cmpd4","Cmpd5","Cmpd6",... ...
    "Cmpd1","1","0.04","0.25","0.13","0.11","0.2",... ...
    "Cmpd2","0.04","1","0.06","0.05","0.19","0.07",... ...
    "Cmpd3","0.25","0.06","1","0.12","0.22","0.25",... ...
    "Cmpd4","0.13","0.05","0.12","1","0.11","0.13",... ...
    "Cmpd5","0.11","0.19","0.22","0.11","1","0.17",... ...
    "Cmpd6","0.2","0.07","0.25","0.13","0.17","1",... ...
    ... ... ..
    ... ... ..
    ... ... ..

Example of data in output file for I<RowsAndColumns> B<--OutMatrixFormat> value for
I<UpperTriangularMatrix> value of B<--OutMatrixType>:

    "","Cmpd1","Cmpd2","Cmpd3","Cmpd4","Cmpd5","Cmpd6",... ...
    "Cmpd1","1","0.04","0.25","0.13","0.11","0.2",... ...
    "Cmpd2","1","0.06","0.05","0.19","0.07",... ...
    "Cmpd3","1","0.12","0.22","0.25",... ...
    "Cmpd4","1","0.11","0.13",... ...
    "Cmpd5","1","0.17",... ...
    "Cmpd6","1",... ...
    ... ... ..
    ... ... ..
    ... ... ..

Example of data in output file for I<RowsAndColumns> B<--OutMatrixFormat> value for
I<LowerTriangularMatrix> value of B<--OutMatrixType>:

    "","Cmpd1","Cmpd2","Cmpd3","Cmpd4","Cmpd5","Cmpd6",... ...
    "Cmpd1","1"
    "Cmpd2","0.04","1"
    "Cmpd3","0.25","0.06","1"
    "Cmpd4","0.13","0.05","0.12","1"
    "Cmpd5","0.11","0.19","0.22","0.11","1"
    "Cmpd6","0.2","0.07","0.25","0.13","0.17","1"
    ... ... ..
    ... ... ..
    ... ... ..


Example of data in output file for I<IDPairsAndValue> B<--OutMatrixFormat> value for
<FullMatrix> value of  B<OutMatrixType>:

    "CmpdID1","CmpdID2","Coefficient Value"
    "Cmpd1","Cmpd1","1"
    "Cmpd1","Cmpd2","0.04"
    "Cmpd1","Cmpd3","0.25"
    "Cmpd1","Cmpd4","0.13"
    ... ... ...
    ... ... ...
    ... ... ...
    "Cmpd2","Cmpd1","0.04"
    "Cmpd2","Cmpd2","1"
    "Cmpd2","Cmpd3","0.06"
    "Cmpd2","Cmpd4","0.05"
    ... ... ...
    ... ... ...
    ... ... ...
    "Cmpd3","Cmpd1","0.25"
    "Cmpd3","Cmpd2","0.06"
    "Cmpd3","Cmpd3","1"
    "Cmpd3","Cmpd4","0.12"
    ... ... ...
    ... ... ...
    ... ... ...

Example of data in output file for I<IDPairsAndValue> B<--OutMatrixFormat> value for
<UpperTriangularMatrix> value of B<--OutMatrixType>:

    "CmpdID1","CmpdID2","Coefficient Value"
    "Cmpd1","Cmpd1","1"
    "Cmpd1","Cmpd2","0.04"
    "Cmpd1","Cmpd3","0.25"
    "Cmpd1","Cmpd4","0.13"
    ... ... ...
    ... ... ...
    ... ... ...
    "Cmpd2","Cmpd2","1"
    "Cmpd2","Cmpd3","0.06"
    "Cmpd2","Cmpd4","0.05"
    ... ... ...
    ... ... ...
    ... ... ...
    "Cmpd3","Cmpd3","1"
    "Cmpd3","Cmpd4","0.12"
    ... ... ...
    ... ... ...
    ... ... ...

Example of data in output file for I<IDPairsAndValue> B<--OutMatrixFormat> value for
<LowerTriangularMatrix> value of B<--OutMatrixType>:

    "CmpdID1","CmpdID2","Coefficient Value"
    "Cmpd1","Cmpd1","1"
    "Cmpd2","Cmpd1","0.04"
    "Cmpd2","Cmpd2","1"
    "Cmpd3","Cmpd1","0.25"
    "Cmpd3","Cmpd2","0.06"
    "Cmpd3","Cmpd3","1"
    "Cmpd4","Cmpd1","0.13"
    "Cmpd4","Cmpd2","0.05"
    "Cmpd4","Cmpd3","0.12"
    "Cmpd4","Cmpd4","1"
    ... ... ...
    ... ... ...
    ... ... ...

=item B<--OutMatrixType> I<FullMatrix | UpperTriangularMatrix | LowerTriangularMatrix>

Type of similarity or distance matrix to calculate for fingerprints vector and bit-vector strings:
Calculate full matrix; Calculate lower triangular matrix including diagonal; Calculate upper triangular
matrix including diagonal.

Possible values: I<FullMatrix, UpperTriangularMatrix, or LowerTriangularMatrix>. Default value:
I<FullMatrix>.

The value of B<--OutMatrixType> in conjunction with B<--OutMatrixFormat> determines type
of data written to output files.

=item B<-o, --overwrite>

Overwrite existing files

=item B<-p, --precision> I<number>

Precision of calculated values in the output file. Default: up to I<2> decimal places.
Valid values: positive integers.

=item B<-q, --quote> I<Yes | No>

Put quote around column values in output CSV/TSV text file(s). Possible values:
I<Yes or No>. Default value: I<Yes>.

=item B<-r, --root> I<RootName>

New file name is generated using the root: <Root><BitVectorComparisonMode>.<Ext> or
<Root><VectorComparisonMode><VectorComparisonFormulism>.<Ext>.
The csv, and tsv <Ext> values are used for comma/semicolon, and tab delimited text files
respectively. This option is ignored for multiple input files.

=item B<-v, --VectorComparisonMode> I<All | "TanimotoSimilarity,[ManhattanDistance,...]">

Specify what similarity or distance coefficients to use for calculating similarity matrices for
fingerprint vector strings data values in I<TextFile(s)>: calculate similarity matrices for all
supported similarity and distance coefficients or specify a comma delimited list of similarity
and distance coefficients. Possible values: I<All | "TanimotoSimilairy,[ManhattanDistance,..]">.
Default: I<TanimotoSimilarity>.

The value of B<-v, --VectorComparisonMode>, in conjunction with B<--VectorComparisonFormulism>,
decides which type of similarity and distance coefficient formulism gets used.

I<All> uses complete list of supported similarity and distance coefficients: I<CosineSimilarity,
CzekanowskiSimilarity, DiceSimilarity, OchiaiSimilarity, JaccardSimilarity, SorensonSimilarity, TanimotoSimilarity,
CityBlockDistance, EuclideanDistance, HammingDistance, ManhattanDistance, SoergelDistance>. These
similarity and distance coefficients are described below.

B<FingerprintsVector.pm> module, used to calculate similarity and distance coefficients,
provides support to perform comparison between vectors containing three different types of
values:

Type I: OrderedNumericalValues

    . Size of two vectors are same
    . Vectors contain real values in a specific order. For example: MACCS keys
      count, Topological pharmnacophore atom pairs and so on.

Type II: UnorderedNumericalValues

    . Size of two vectors might not be same
    . Vectors contain unordered real value identified by value IDs. For example:
      Toplogical atom pairs, Topological atom torsions and so on

Type III: AlphaNumericalValues

    . Size of two vectors might not be same
    . Vectors contain unordered alphanumerical values. For example: Extended
      connectivity fingerprints, atom neighborhood fingerprints.

Before performing similarity or distance calculations between vectors containing UnorderedNumericalValues
or AlphaNumericalValues, the vectors are transformed into vectors containing unique OrderedNumericalValues
using value IDs for UnorderedNumericalValues and values itself for AlphaNumericalValues.

Three forms of similarity and distance calculation between two vectors, specified using B<--VectorComparisonFormulism>
option, are supported: I<AlgebraicForm, BinaryForm or SetTheoreticForm>.

For I<BinaryForm>, the ordered list of processed final vector values containing the value or
count of each unique value type is simply converted into a binary vector containing 1s and 0s
corresponding to presence or absence of values before calculating similarity or distance between
two vectors.

For two fingerprint vectors A and B of same size containing OrderedNumericalValues, let:

    N = Number values in A or B

    Xa = Values of vector A
    Xb = Values of vector B

    Xai = Value of ith element in A
    Xbi = Value of ith element in B

   SUM = Sum of i over N values

For SetTheoreticForm of calculation between two vectors, let:

    SetIntersectionXaXb = SUM ( MIN ( Xai, Xbi ) )
    SetDifferenceXaXb = SUM ( Xai ) + SUM ( Xbi ) - SUM ( MIN ( Xai, Xbi ) )

For BinaryForm of calculation between two vectors, let:

    Na = Number of bits set to "1" in A = SUM ( Xai )
    Nb = Number of bits set to "1" in B = SUM ( Xbi )
    Nc = Number of bits set to "1" in both A and B = SUM ( Xai * Xbi )
    Nd = Number of bits set to "0" in both A and B
       = SUM ( 1 - Xai - Xbi + Xai * Xbi)

    N = Number of bits set to "1" or "0" in A or B = Size of A or B = Na + Nb - Nc + Nd

Additionally, for BinaryForm various values also correspond to:

    Na = | Xa |
    Nb = | Xb |
    Nc = | SetIntersectionXaXb |
    Nd = N - | SetDifferenceXaXb |

    | SetDifferenceXaXb | = N - Nd = Na + Nb - Nc + Nd - Nd = Na + Nb - Nc
                          =  | Xa | + | Xb | - | SetIntersectionXaXb |

Various similarity and distance coefficients [ Ref 40, Ref 62, Ref 64 ] for a pair of vectors A and B
in I<AlgebraicForm, BinaryForm and SetTheoreticForm> are defined as follows:

B<CityBlockDistance>: ( same as HammingDistance and ManhattanDistance)

I<AlgebraicForm>: SUM ( ABS ( Xai - Xbi ) )

I<BinaryForm>: ( Na - Nc ) + ( Nb - Nc ) = Na + Nb - 2 * Nc

I<SetTheoreticForm>: | SetDifferenceXaXb | - | SetIntersectionXaXb | = SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) )

B<CosineSimilarity>:  ( same as OchiaiSimilarityCoefficient)

I<AlgebraicForm>: SUM ( Xai * Xbi ) / SQRT ( SUM ( Xai ** 2) * SUM ( Xbi ** 2) )

I<BinaryForm>: Nc / SQRT ( Na * Nb)

I<SetTheoreticForm>: | SetIntersectionXaXb | / SQRT ( |Xa| * |Xb| ) = SUM ( MIN ( Xai, Xbi ) ) / SQRT ( SUM ( Xai ) * SUM ( Xbi ) )

B<CzekanowskiSimilarity>: ( same as DiceSimilarity and SorensonSimilarity)

I<AlgebraicForm>: ( 2 * ( SUM ( Xai * Xbi ) )  ) / ( SUM ( Xai ** 2) + SUM ( Xbi **2 ) )

I<BinaryForm>: 2 * Nc / ( Na + Nb )

I<SetTheoreticForm>: 2 * | SetIntersectionXaXb | / ( |Xa| + |Xb| ) = 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) / ( SUM ( Xai ) + SUM ( Xbi ) )

B<DiceSimilarity>: ( same as CzekanowskiSimilarity and SorensonSimilarity)

I<AlgebraicForm>: ( 2 * ( SUM ( Xai * Xbi ) )  ) / ( SUM ( Xai ** 2) + SUM ( Xbi **2 ) )

I<BinaryForm>: 2 * Nc / ( Na + Nb )

I<SetTheoreticForm>: 2 * | SetIntersectionXaXb | / ( |Xa| + |Xb| ) = 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) / ( SUM ( Xai ) + SUM ( Xbi ) )

B<EuclideanDistance>:

I<AlgebraicForm>: SQRT ( SUM ( ( ( Xai - Xbi ) ** 2 ) ) )

I<BinaryForm>: SQRT ( ( Na - Nc ) + ( Nb - Nc ) ) = SQRT ( Na + Nb - 2 * Nc )

I<SetTheoreticForm>: SQRT ( | SetDifferenceXaXb | - | SetIntersectionXaXb | ) = SQRT (  SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) )

B<HammingDistance>:  ( same as CityBlockDistance and ManhattanDistance)

I<AlgebraicForm>: SUM ( ABS ( Xai - Xbi ) )

I<BinaryForm>: ( Na - Nc ) + ( Nb - Nc ) = Na + Nb - 2 * Nc

I<SetTheoreticForm>: | SetDifferenceXaXb | - | SetIntersectionXaXb | = SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) )

B<JaccardSimilarity>: ( same as TanimotoSimilarity)

I<AlgebraicForm>:  SUM ( Xai * Xbi ) / ( SUM ( Xai ** 2 ) + SUM ( Xbi ** 2 ) - SUM ( Xai * Xbi ) )

I<BinaryForm>:  Nc / ( ( Na - Nc ) + ( Nb - Nc ) + Nc ) = Nc / ( Na + Nb - Nc )

I<SetTheoreticForm>: | SetIntersectionXaXb | / | SetDifferenceXaXb | = SUM ( MIN ( Xai, Xbi ) ) / (  SUM ( Xai ) + SUM ( Xbi ) - SUM ( MIN ( Xai, Xbi ) ) )

B<ManhattanDistance>:  ( same as CityBlockDistance and HammingDistance)

I<AlgebraicForm>: SUM ( ABS ( Xai - Xbi ) )

I<BinaryForm>: ( Na - Nc ) + ( Nb - Nc ) = Na + Nb - 2 * Nc

I<SetTheoreticForm>: | SetDifferenceXaXb | - | SetIntersectionXaXb | = SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) )

B<OchiaiSimilarity>:  ( same as CosineSimilarity)

I<AlgebraicForm>: SUM ( Xai * Xbi ) / SQRT ( SUM ( Xai ** 2) * SUM ( Xbi ** 2) )

I<BinaryForm>: Nc / SQRT ( Na * Nb)

I<SetTheoreticForm>: | SetIntersectionXaXb | / SQRT ( |Xa| * |Xb| ) = SUM ( MIN ( Xai, Xbi ) ) / SQRT ( SUM ( Xai ) * SUM ( Xbi ) )

B<SorensonSimilarity>: ( same as CzekanowskiSimilarity and DiceSimilarity)

I<AlgebraicForm>: ( 2 * ( SUM ( Xai * Xbi ) )  ) / ( SUM ( Xai ** 2) + SUM ( Xbi **2 ) )

I<BinaryForm>: 2 * Nc / ( Na + Nb )

I<SetTheoreticForm>: 2 * | SetIntersectionXaXb | / ( |Xa| + |Xb| ) = 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) / ( SUM ( Xai ) + SUM ( Xbi ) )

B<SoergelDistance>:

I<AlgebraicForm>:  SUM ( ABS ( Xai - Xbi ) ) / SUM ( MAX ( Xai, Xbi ) )

I<BinaryForm>: 1 - Nc / ( Na + Nb - Nc ) = ( Na + Nb - 2 * Nc ) / ( Na + Nb - Nc )

I<SetTheoreticForm>: ( | SetDifferenceXaXb | - | SetIntersectionXaXb | ) / | SetDifferenceXaXb | = ( SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) ) / ( SUM ( Xai ) + SUM ( Xbi ) - SUM ( MIN ( Xai, Xbi ) ) )

B<TanimotoSimilarity>:  ( same as JaccardSimilarity)

I<AlgebraicForm>:  SUM ( Xai * Xbi ) / ( SUM ( Xai ** 2 ) + SUM ( Xbi ** 2 ) - SUM ( Xai * Xbi ) )

I<BinaryForm>:  Nc / ( ( Na - Nc ) + ( Nb - Nc ) + Nc ) = Nc / ( Na + Nb - Nc )

I<SetTheoreticForm>: | SetIntersectionXaXb | / | SetDifferenceXaXb | = SUM ( MIN ( Xai, Xbi ) ) / (  SUM ( Xai ) + SUM ( Xbi ) - SUM ( MIN ( Xai, Xbi ) ) )

=item B<--VectorComparisonFormulism> I<All | "AlgebraicForm,[BinaryForm,SetTheoreticForm]">

Specify fingerprints vector comparison formulism to use for calculation similarity and distance
coefficients during B<-v, --VectorComparisonMode>: use all supported comparison formulisms
or specify a comma delimited. Possible values: I<All | "AlgebraicForm,[BinaryForm,SetTheoreticForm]">.
Default value: I<AlgebraicForm>.

I<All> uses all three forms of supported vector comparison formulism for values of B<-v, --VectorComparisonMode>
option.

For fingerprint vector strings containing B<AlphaNumericalValues> data values - B<ExtendedConnectivityFingerprints>,
B<AtomNeighborhoodsFingerprints> and so on - all three formulism result in same value during similarity and distance
calculations.

=item B<-w, --WorkingDir> I<DirName>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To generate a similarity matrix corresponding to Tanimoto similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints in text file present in a column
name containing Fingerprint substring by loading all fingerprints data into memory and create a
SampleFPHexTanimotoSimilarity.csv file containing compound IDs retrieved from column name
containing CompoundID substring, type:

    % SimilarityMatricesFingerprints.pl -o SampleFPHex.csv

To generate a similarity matrix corresponding to Tanimoto similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints in SD File present in a data field
with Fingerprint substring in its label by loading all fingerprints data into memory and create a
SampleFPHexTanimotoSimilarity.csv file containing sequentially generated compound IDs with
Cmpd prefix, type:

    % SimilarityMatricesFingerprints.pl -o SampleFPHex.sdf

To generate a similarity matrix corresponding to Tanimoto similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints in FP file by loading all fingerprints
data into memory and create a SampleFPHexTanimotoSimilarity.csv file along with compound IDs
retrieved from FP file, type:

    % SimilarityMatricesFingerprints.pl -o SampleFPHex.fpf

To generate a lower triangular similarity matrix corresponding to Tanimoto similarity coefficient for
fingerprints bit-vector strings data corresponding to supported fingerprints in text file present in a
column name containing Fingerprint substring by loading all fingerprints data into memory  and create
a SampleFPHexTanimotoSimilarity.csv file containing compound IDs retrieved from column name
containing CompoundID substring, type:

    % SimilarityMatricesFingerprints.pl -o --InputDataMode LoadInMemory
      --OutMatrixFormat RowsAndColumns --OutMatrixType LowerTriangularMatrix
      SampleFPHex.csv

To generate a upper triangular similarity matrix corresponding to Tanimoto similarity coefficient for
fingerprints bit-vector strings data corresponding to supported fingerprints in text file present in a
column name containing Fingerprint substring by loading all fingerprints data into memory  and create
a SampleFPHexTanimotoSimilarity.csv file in IDPairsAndValue format containing compound IDs retrieved
from column name containing CompoundID substring, type:

    % SimilarityMatricesFingerprints.pl -o --InputDataMode LoadInMemory
      --OutMatrixFormat IDPairsAndValue --OutMatrixType UpperTriangularMatrix
      SampleFPHex.csv

To generate a full similarity matrix corresponding to Tanimoto similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints in text file present in a column
name containing Fingerprint substring by scanning file without loading all fingerprints data into memory
and create a SampleFPHexTanimotoSimilarity.csv file containing compound IDs retrieved from
column name containing CompoundID substring, type:

    % SimilarityMatricesFingerprints.pl -o --InputDataMode ScanFile
      --OutMatrixFormat RowsAndColumns --OutMatrixType FullMatrix
      SampleFPHex.csv

To generate a lower triangular similarity matrix corresponding to Tanimoto similarity coefficient for
fingerprints bit-vector strings data corresponding to supported fingerprints in text file present in a
column name containing Fingerprint substring by scanning file without loading all fingerprints data into
memory and create a SampleFPHexTanimotoSimilarity.csv file in IDPairsAndValue format containing
compound IDs retrieved from column name containing CompoundID substring, type:

    % SimilarityMatricesFingerprints.pl -o --InputDataMode ScanFile
      --OutMatrixFormat IDPairsAndValue --OutMatrixType LowerTriangularMatrix
      SampleFPHex.csv

To generate a similarity matrix corresponding to Tanimoto similarity coefficient using algebraic formulism
for fingerprints vector strings data corresponding to supported fingerprints in text file present in a column name
containing Fingerprint substring and create a SampleFPCountTanimotoSimilarityAlgebraicForm.csv file
containing compound IDs retrieved from column name containing CompoundID substring, type:

    % SimilarityMatricesFingerprints.pl -o SampleFPCount.csv

To generate a similarity matrix corresponding to Tanimoto similarity coefficient using algebraic formulism
for fingerprints vector strings data corresponding to supported fingerprints in SD file present in a data field with
Fingerprint substring in its label and create a SampleFPCountTanimotoSimilarityAlgebraicForm.csv file
containing sequentially generated compound IDs with Cmpd prefix, type:

    % SimilarityMatricesFingerprints.pl -o SampleFPCount.sdf

To generate a similarity matrix corresponding to Tanimoto similarity coefficient using algebraic formulism
vector strings data corresponding to supported fingerprints in FP file and create a
SampleFPCountTanimotoSimilarityAlgebraicForm.csv file along with compound IDs retrieved from FP file, type:

    % SimilarityMatricesFingerprints.pl -o SampleFPCount.fpf

To generate a similarity matrix corresponding to Tanimoto similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints in text file present in a column name
containing Fingerprint substring and create a SampleFPHexTanimotoSimilarity.csv file in
IDPairsAndValue format containing compound IDs retrieved from column name containing
CompoundID substring, type:

    % SimilarityMatricesFingerprints.pl --OutMatrixFormat IDPairsAndValue -o
      SampleFPHex.csv

To generate a similarity matrix corresponding to Tanimoto similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints in SD file present in a data field with
Fingerprint substring in its label and create a SampleFPHexTanimotoSimilarity.csv file in
IDPairsAndValue format containing sequentially generated compound IDs with Cmpd prefix,
type:

    % SimilarityMatricesFingerprints.pl --OutMatrixFormat IDPairsAndValue -o
      SampleFPHex.sdf

To generate a similarity matrix corresponding to Tanimoto similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints in FP file and create a
SampleFPHexTanimotoSimilarity.csv file in IDPairsAndValue format along with compound IDs retrieved
from FP file, type:

    % SimilarityMatricesFingerprints.pl --OutMatrixFormat IDPairsAndValue -o
      SampleFPHex.fpf

To generate a similarity matrix corresponding to Tanimoto similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints in SD file present in a data field with
Fingerprint substring in its label and create a SampleFPHexTanimotoSimilarity.csv file
containing compound IDs from mol name line, type:

    % SimilarityMatricesFingerprints.pl --CompoundIDMode MolName -o
      SampleFPHex.sdf

To generate a similarity matrix corresponding to Tanimoto similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints present in a data field with
Fingerprint substring in its label and create a SampleFPHexTanimotoSimilarity.csv file
containing compound IDs from data field name Mol_ID, type:

    % SimilarityMatricesFingerprints.pl --CompoundIDMode DataField
      --CompoundIDField Mol_ID -o SampleFPBin.sdf

To generate similarity matrices corresponding to Buser, Dice and Tanimoto similarity coefficient
for fingerprints bit-vector strings data corresponding to supported fingerprints present in a column
name containing Fingerprint substring and create SampleFPBin[CoefficientName]Similarity.csv files
containing compound IDs retrieved from column name containing CompoundID substring, type:

    % SimilarityMatricesFingerprints.pl -b "BuserSimilarity,DiceSimilarity,
      TanimotoSimilarity" -o SampleFPBin.csv

To generate similarity matrices corresponding to Buser, Dice and Tanimoto similarity coefficient
for fingerprints bit-vector strings data corresponding to supported fingerprints present in a data field with
Fingerprint substring in its label and create SampleFPBin[CoefficientName]Similarity.csv files
containing sequentially generated compound IDs with Cmpd prefix, type:

    % SimilarityMatricesFingerprints.pl -b "BuserSimilarity,DiceSimilarity,
      TanimotoSimilarity" -o SampleFPBin.sdf

To generate similarity matrices corresponding to CityBlock distance and Tanimoto similarity coefficients using
algebraic formulism for fingerprints vector strings data corresponding to supported fingerprints present in
a column name containing Fingerprint substring and create  SampleFPCount[CoefficientName]AlgebraicForm.csv
files containing compound IDs retrieved from column name containing CompoundID substring, type:

    % SimilarityMatricesFingerprints.pl -v "CityBlockDistance,
      TanimotoSimilarity" -o SampleFPCount.csv

To generate similarity matrices corresponding to CityBlock distance and Tanimoto similarity coefficients using
algebraic formulism for fingerprints vector strings data corresponding to supported fingerprints present in
a data field with Fingerprint substring in its label and create SampleFPCount[CoefficientName]AlgebraicForm.csv
files containing sequentially generated compound IDs with Cmpd prefix, type:

    % SimilarityMatricesFingerprints.pl -v "CityBlockDistance,
      TanimotoSimilarity" -o SampleFPCount.sdf

To generate similarity matrices corresponding to CityBlock distance Tanimoto similarity coefficients using
binary formulism for fingerprints vector strings data corresponding to supported fingerprints present in
a column name containing Fingerprint substring and create SampleFPCount[CoefficientName]Binary.csv
files containing compound IDs retrieved from column name containing CompoundID substring, type:

    % SimilarityMatricesFingerprints.pl -v "CityBlockDistance,
      TanimotoSimilarity" --VectorComparisonFormulism BinaryForm -o
      SampleFPCount.csv

To generate similarity matrices corresponding to CityBlock distance Tanimoto similarity coefficients using
binary formulism for fingerprints vector strings data corresponding to supported fingerprints present in
a data field with Fingerprint substring in its label  and create SampleFPCount[CoefficientName]Binary.csv
files containing sequentially generated compound IDs with Cmpd prefix, type:

    % SimilarityMatricesFingerprints.pl -v "CityBlockDistance,
      TanimotoSimilarity" --VectorComparisonFormulism BinaryForm -o
      SampleFPCount.sdf

To generate similarity matrices corresponding to CityBlock distance Tanimoto similarity coefficients using
all supported comparison formulisms for fingerprints vector strings data corresponding to supported
fingerprints present in a column name containing Fingerprint substring and create
SampleFPCount[CoefficientName][FormulismName].csv files containing compound IDs retrieved from column
name containing CompoundID substring, type:

    % SimilarityMatricesFingerprints.pl -v "CityBlockDistance,
      TanimotoSimilarity" --VectorComparisonFormulism All -o SampleFPCount.csv

To generate similarity matrices corresponding to CityBlock distance Tanimoto similarity coefficients using
all supported comparison formulisms for fingerprints vector strings data corresponding to supported
fingerprints present in a data field with Fingerprint substring in its label and create
SampleFPCount[CoefficientName][FormulismName].csv files containing  sequentially generated
compound IDs with Cmpd prefix, type:

    % SimilarityMatricesFingerprints.pl -v "CityBlockDistance,TanimotoSimilarity"
      --VectorComparisonFormulism All -o SampleFPCount.sdf

To generate similarity matrices corresponding to all available similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints present in a column name
containing Fingerprint substring and create SampleFPHex[CoefficientName].csv files
containing compound IDs retrieved from column name containing CompoundID substring, type:

    % SimilarityMatricesFingerprints.pl -m AutoDetect --BitVectorComparisonMode
      All --alpha 0.5 -beta 0.5 -o SampleFPHex.csv

To generate similarity matrices corresponding to all available similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints present in a data field with Fingerprint
substring in its label and create SampleFPHex[CoefficientName].csv files containing  sequentially
generated compound IDs with Cmpd prefix, type

    % SimilarityMatricesFingerprints.pl -m AutoDetect --BitVectorComparisonMode
      All --alpha 0.5 -beta 0.5 -o SampleFPHex.sdf

To generate similarity matrices corresponding to all available similarity and distance coefficients using
all comparison formulism for fingerprints vector strings data corresponding to supported fingerprints
present in a column name containing Fingerprint substring and create
SampleFPCount[CoefficientName][FormulismName].csv files containing compound IDs
retrieved from column name containing CompoundID substring, type:

    % SimilarityMatricesFingerprints.pl -m AutoDetect --VectorComparisonMode
      All --VectorComparisonFormulism All -o SampleFPCount.csv

To generate similarity matrices corresponding to all available similarity and distance coefficients using
all comparison formulism for fingerprints vector strings data corresponding to supported fingerprints
present in a data field with Fingerprint substring in its label and create
SampleFPCount[CoefficientName][FormulismName].csv files containing  sequentially generated
compound IDs with Cmpd prefix, type:

    % SimilarityMatricesFingerprints.pl -m AutoDetect --VectorComparisonMode
      All --VectorComparisonFormulism All -o SampleFPCount.sdf

To generate a similarity matrix corresponding to Tanimoto similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints present in a column number 2
and create a SampleFPHexTanimotoSimilarity.csv file containing compound IDs retrieved column
number 1, type:

    % SimilarityMatricesFingerprints.pl --ColMode ColNum --CompoundIDCol 1
      --FingerprintsCol 2 -o SampleFPHex.csv

To generate a similarity matrix corresponding to Tanimoto similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints present in a data field name
Fingerprints and create a SampleFPHexTanimotoSimilarity.csv file containing compound IDs
present in data field name Mol_ID, type:

    % SimilarityMatricesFingerprints.pl --FingerprintsField Fingerprints
      --CompoundIDMode DataField --CompoundIDField Mol_ID -o SampleFPHex.sdf

To generate a similarity matrix corresponding to Tversky similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints present in a column named Fingerprints
and create a SampleFPHexTverskySimilarity.tsv file containing compound IDs retrieved column named
CompoundID, type:

    % SimilarityMatricesFingerprints.pl --BitVectorComparisonMode
      TverskySimilarity --alpha 0.5 --ColMode ColLabel --CompoundIDCol
      CompoundID --FingerprintsCol Fingerprints --OutDelim Tab --quote No
      -o SampleFPHex.csv

To generate a similarity matrix corresponding to Tanimoto similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints present in a data field
with Fingerprint substring in its label and create a SampleFPHexTanimotoSimilarity.csv file
containing compound IDs from molname line or sequentially generated compound IDs
with Mol prefix, type:

    % SimilarityMatricesFingerprints.pl --CompoundIDMode MolnameOrLabelPrefix
      --CompoundIDPrefix Mol -o SampleFPHex.sdf

To generate a similarity matrix corresponding to Tanimoto similarity coefficient for fingerprints
bit-vector strings data corresponding to supported fingerprints present in a data field with
Fingerprint substring in its label and create a SampleFPHexTanimotoSimilarity.tsv file
containing sequentially generated compound IDs with Cmpd prefix, type:

    % SimilarityMatricesFingerprints.pl -OutDelim Tab --quote No -o SampleFPHex.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoFingerprintsFiles.pl, SimilaritySearchingFingerprints.pl, AtomNeighborhoodsFingerprints.pl,
ExtendedConnectivityFingerprints.pl, MACCSKeysFingerprints.pl, PathLengthFingerprints.pl,
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
