#!/usr/bin/perl -w
#
# File: SimilaritySearchingFingerprints.pl
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
use StatisticsUtil;
use PseudoHeap;
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
if ($Options{help} || @ARGV != 2) {
  die GetUsageFromPod("$FindBin::Bin/$ScriptName");
}

# Process reference and database file names...
my(@FingerprintsFilesList);
ProcessFingerprintsFileNames();

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Setup information about fingerprints inut and SD/text output files...
my(%FingerprintsFilesInfo, %OutputFilesInfo, %SimilaritySearchInfo);
print "Checking and retrieving information from reference and database fingerprints files...\n";
RetrieveFingerprintsFilesInfo();

# Perform similarity search...
print "Performing similarity search...\n";
my(%SimilaritySearchResults, %DatabaseFingerprintsFileData);
PerformSimilaritySearch();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Perform similarity search using fingerprints data in reference and database text files...
#
sub PerformSimilaritySearch {

  print "\nProcessing fingerprints data for reference molecules...\n";
  ReadReferenceFingerprintsData();

  InitializeSimilaritySearchResults();
  GenerateSimilaritySearchResults();
  WriteSimilaritySearchResultFiles();
}

# Find similar molecules from database molecules for individual or multiple reference molecules...
#
sub GenerateSimilaritySearchResults {
  my($DatabaseFingerprintsFileIO, $FingerprintsCount, $IgnoredFingerprintsCount, $DatabaseFingerprintsObject, $DatabaseCmpdID, $ReferenceFingerprintsObject, $ReferenceIndex, $ReferenceCmpdID, $ComparisonValue, $FusedComparisonValue, @ComparisonValues);

  print "Processing fingerprints data for database molecules...\n";

  ($FingerprintsCount, $IgnoredFingerprintsCount) = (0) x 3;

  $DatabaseFingerprintsFileIO = Fingerprints::FingerprintsFileUtil::NewFingerprintsFileIO(%{$FingerprintsFilesInfo{Database}{FingerprintsFileIOParameters}});
  $DatabaseFingerprintsFileIO->Open();

  @ComparisonValues = ();

  DATABASEFP: while ($DatabaseFingerprintsFileIO->Read()) {
    $FingerprintsCount++;

    if (!$DatabaseFingerprintsFileIO->IsFingerprintsDataValid()) {
      $IgnoredFingerprintsCount++;
      next DATABASEFP;
    }
    $DatabaseFingerprintsObject = $DatabaseFingerprintsFileIO->GetFingerprints();
    $DatabaseCmpdID = $DatabaseFingerprintsFileIO->GetCompoundID();

    if ($SimilaritySearchInfo{MultipleReferencesMode}) {
      @ComparisonValues = ();
    }

    REFERENCEFP: for $ReferenceIndex (0 .. $#{$SimilaritySearchInfo{ReferenceCmpdIDsRef}}) {
      $ReferenceCmpdID = $SimilaritySearchInfo{ReferenceCmpdIDsRef}->[$ReferenceIndex];
      $ReferenceFingerprintsObject = $SimilaritySearchInfo{ReferenceFingerprintsObjectsRef}->[$ReferenceIndex];

      $ComparisonValue = CompareReferenceAndDatabaseFingerprintsPair($ReferenceFingerprintsObject, $DatabaseFingerprintsObject);
      if (!defined $ComparisonValue) {
	next REFERENCEFP;
      }

      if ($SimilaritySearchInfo{IndividualReferenceMode}) {
	CollectSimilaritySearchResults($DatabaseFingerprintsFileIO, $DatabaseCmpdID, $ComparisonValue, $ReferenceCmpdID);
      }
      elsif ($SimilaritySearchInfo{MultipleReferencesMode}) {
	push @ComparisonValues, $ComparisonValue;
      }
    }

    if ($SimilaritySearchInfo{MultipleReferencesMode}) {
      $FusedComparisonValue = CalculateGroupFusionComparisonValue(\@ComparisonValues);
      if (!defined $FusedComparisonValue) {
	next DATABASEFP;
      }
      CollectSimilaritySearchResults($DatabaseFingerprintsFileIO, $DatabaseCmpdID, $FusedComparisonValue);
    }
  }
  $DatabaseFingerprintsFileIO->Close();

  print "Number of fingerprints data entries in database fingerprints file: $FingerprintsCount\n";
  print "Number of fingerprints date entries processed successfully: ", ($FingerprintsCount - $IgnoredFingerprintsCount)  , "\n";
  print "Number of fingerprints data entries ignored due to missing/invalid data: $IgnoredFingerprintsCount\n\n";
}

# Compare a pair of reference and database fingerprints objects corresponding to bit-vector or
# vectors using specified comparison method and comparison cutoff...
#
sub CompareReferenceAndDatabaseFingerprintsPair {
  my($ReferenceFingerprintsObject, $DatabaseFingerprintsObject) = @_;
  my($ComparisonMethod, $ComparisonValue);

  $ComparisonMethod = $SimilaritySearchInfo{ComparisonMethod};
  $ComparisonValue = $ReferenceFingerprintsObject->$ComparisonMethod($DatabaseFingerprintsObject, @{$SimilaritySearchInfo{ComparisonMethodParameters}});

  if (!defined $ComparisonValue) {
    warn "Warning: Ignoring fingerprints data for reference compound ID ",  $ReferenceFingerprintsObject->GetID(), ": Its comparison with database compound ID, ", $DatabaseFingerprintsObject->GetID(), ", failed.\n";
    return undef;
  }

  $ComparisonValue = sprintf("%.$OptionsInfo{Precision}f", $ComparisonValue);

  # Apply any comparison cutoff...
  if ($SimilaritySearchInfo{ApplyComparisonCutoff}) {
    return $SimilaritySearchInfo{KeepTop} ? ($ComparisonValue >= $SimilaritySearchInfo{ComparisonCutoff} ? $ComparisonValue : undef) : ($ComparisonValue <= $SimilaritySearchInfo{ComparisonCutoff} ? $ComparisonValue : undef);
  }
  else {
    return $ComparisonValue;
  }
}

# Calculate group fusion comparison value...
#
sub CalculateGroupFusionComparisonValue {
  my($ComparisonValuesRef) = @_;
  my($FusedComparisonValue, @ComparisonValues);

  if (!@{$ComparisonValuesRef}) {
    return undef;
  }

  if ($SimilaritySearchInfo{SortComparisonValues}) {
    @ComparisonValues = sort { $SimilaritySearchInfo{KeepTop} ? ($b <=> $a) : ($a <=> $b) } @{$ComparisonValuesRef};
    if ($SimilaritySearchInfo{UsekNN} && ($OptionsInfo{kNN} < scalar @{$ComparisonValuesRef})) {
      # Keep only top kNN values for group fusion...
      splice @ComparisonValues, $OptionsInfo{kNN};
    }
    $ComparisonValuesRef = \@ComparisonValues;
  }

  $FusedComparisonValue = &{$SimilaritySearchInfo{GroupFusionMethodRef}}($ComparisonValuesRef);
  if ($SimilaritySearchInfo{ApplyPrecisionDuringFusion}) {
    $FusedComparisonValue = sprintf("%.$OptionsInfo{Precision}f", $FusedComparisonValue);
  }

  return $FusedComparisonValue;
}

# Collect similarity results for individual reference and multiple references search...
#
sub CollectSimilaritySearchResults {
  my($DatabaseFingerprintsFileIO, $DatabaseCmpdID, $ComparisonValue, $ReferenceCmpdID) = @_;

  if (defined $ReferenceCmpdID) {
    $SimilaritySearchResults{$ReferenceCmpdID}->AddKeyValuePair($ComparisonValue, $DatabaseCmpdID);
  }
  else {
    $SimilaritySearchResults{ResultsPseudoHeap}->AddKeyValuePair($ComparisonValue, $DatabaseCmpdID);
  }

  if ($FingerprintsFilesInfo{Database}{CollectInputFileData}) {
    CollectDatabaseFileData($DatabaseCmpdID, $DatabaseFingerprintsFileIO);
  }
}

# Initialize similarity results for individual or multiple reference molecules...
#
sub InitializeSimilaritySearchResults {
  my($ReferenceCmpdID);

  %SimilaritySearchResults = ();

  if ($SimilaritySearchInfo{IndividualReferenceMode}) {
    for $ReferenceCmpdID (@{$SimilaritySearchInfo{ReferenceCmpdIDsRef}}) {
      $SimilaritySearchResults{$ReferenceCmpdID} = new PseudoHeap('Type' => ($SimilaritySearchInfo{KeepTop} ? 'KeepTopN' : 'KeepBottomN'), 'KeyType' => 'Numeric', 'MaxSize' => $OptionsInfo{MaxSimilarMolecules});
    }
  }
  elsif ($SimilaritySearchInfo{MultipleReferencesMode}) {
    $SimilaritySearchResults{ResultsPseudoHeap} = new PseudoHeap('Type' => ($SimilaritySearchInfo{KeepTop} ? 'KeepTopN' : 'KeepBottomN'), 'KeyType' => 'Numeric', 'MaxSize' => $OptionsInfo{MaxSimilarMolecules});
  }

  %DatabaseFingerprintsFileData = ();
}

# Write out results SD and/or CSV/TSV text files for individual or multiple reference molecules...
#
sub WriteSimilaritySearchResultFiles {
  my($NewSDFileRef, $NewTextFileRef, $ReferenceCmpdID, $DatabaseCmpdID, $ComparisonValue);

  ($NewSDFileRef, $NewTextFileRef) = SetupAndOpenOutputFiles();

  if ($SimilaritySearchInfo{IndividualReferenceMode}) {
    for $ReferenceCmpdID (@{$SimilaritySearchInfo{ReferenceCmpdIDsRef}}) {
      for $ComparisonValue ($SimilaritySearchResults{$ReferenceCmpdID}->GetSortedKeys()) {
	for $DatabaseCmpdID ($SimilaritySearchResults{$ReferenceCmpdID}->GetKeyValues($ComparisonValue)) {
	  WriteDataToOutputFiles($NewSDFileRef, $NewTextFileRef, $ComparisonValue, $DatabaseCmpdID, $ReferenceCmpdID);
	}
      }
    }
  }
  elsif ($SimilaritySearchInfo{MultipleReferencesMode}) {
    for $ComparisonValue ($SimilaritySearchResults{ResultsPseudoHeap}->GetSortedKeys()) {
      for $DatabaseCmpdID ($SimilaritySearchResults{ResultsPseudoHeap}->GetKeyValues($ComparisonValue)) {
	WriteDataToOutputFiles($NewSDFileRef, $NewTextFileRef, $ComparisonValue, $DatabaseCmpdID);
      }
    }
  }

  if ($NewSDFileRef) {
    close $NewSDFileRef;
  }
  if ($NewTextFileRef) {
    close $NewTextFileRef;
  }
}

# Write individual reference or multiple references similarity results along with any other data to output files...
#
sub WriteDataToOutputFiles {
  my($NewSDFileRef, $NewTextFileRef, $ComparisonValue, $DatabaseCmpdID, $ReferenceCmpdID) = @_;

  if ($NewSDFileRef) {
    WriteMolStringDataToSDOutputFile($DatabaseCmpdID, $NewSDFileRef);
    if (defined $ReferenceCmpdID) {
      print $NewSDFileRef  ">  <ReferenceCmpdID>\n$ReferenceCmpdID\n\n";
    }
    print $NewSDFileRef  ">  <DatabaseCmpdID>\n$DatabaseCmpdID\n\n>  <ComparisonValue>\n$ComparisonValue\n\n";
    WriteDatabaseDataToSDOutputFile($DatabaseCmpdID, $NewSDFileRef);
    print $NewSDFileRef "\$\$\$\$\n";
  }

  if ($NewTextFileRef) {
    my(@LineWords);

    @LineWords = ();
    if (defined $ReferenceCmpdID) {
      push @LineWords, $ReferenceCmpdID;
    }
    push @LineWords, ($DatabaseCmpdID, $ComparisonValue);

    if ($FingerprintsFilesInfo{Database}{OutputDataFields} || $FingerprintsFilesInfo{Database}{OutputDataCols}) {
      push @LineWords, RetrieveDatabaseDataForTextOutputFile($DatabaseCmpdID);
    }
    print $NewTextFileRef JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote}), "\n";
  }
}

# Open output files...
#
sub SetupAndOpenOutputFiles {
  my($NewSDFileRef, $NewTextFileRef, $NewSDFile, $NewTextFile);

  ($NewSDFileRef, $NewTextFileRef) = (undef) x 2;

  if ($OptionsInfo{SDOutput}) {
    $NewSDFile = $OutputFilesInfo{SDOutFileName};
    print "Generating SD file $NewSDFile...\n";
    open NEWSDFILE, ">$NewSDFile" or die "Error: Couldn't open $NewSDFile: $! \n";
    $NewSDFileRef = \*NEWSDFILE;
  }

  if ($OptionsInfo{TextOutput}) {
    $NewTextFile = $OutputFilesInfo{TextOutFileName};
    print "Generating text file $NewTextFile...\n";
    open NEWTEXTFILE, ">$NewTextFile" or die "Error: Couldn't open $NewTextFile: $! \n";
    $NewTextFileRef = \*NEWTEXTFILE;

    WriteTextFileCoulmnLabels(\*NEWTEXTFILE);
  }

  return ($NewSDFileRef, $NewTextFileRef);
}

# Write out approriate column labels to text file...
#
sub WriteTextFileCoulmnLabels {
  my($NewTextFileRef) = @_;
  my($Line, @LineWords);

  @LineWords = ();

  if ($SimilaritySearchInfo{IndividualReferenceMode}) {
    push @LineWords, qw(ReferenceCompoundID DatabaseCompoundID ComparisonValue);
  }
  elsif ($SimilaritySearchInfo{MultipleReferencesMode}) {
    push @LineWords, qw(DatabaseCompoundID ComparisonValue);
  }

  # Add columns for other database fingerprints file data to be written to output file...
  if ($FingerprintsFilesInfo{Database}{OutputDataFields}) {
    push @LineWords, @{$FingerprintsFilesInfo{Database}{DataFieldsToOutput}};
  }
  elsif ($FingerprintsFilesInfo{Database}{OutputDataCols}) {
    push @LineWords, @{$FingerprintsFilesInfo{Database}{DataColLabelsToOutput}};
  }

  $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print $NewTextFileRef "$Line\n";
}

# Write molecule string data to SD output file...
#
sub WriteMolStringDataToSDOutputFile {
  my($DatabaseCmpdID, $NewSDFileRef) = @_;

  if ($FingerprintsFilesInfo{Database}{CollectCmpdStringData}) {
    my($MolString);

    ($MolString) = split /M  END/, $DatabaseFingerprintsFileData{$DatabaseCmpdID};
    print $NewSDFileRef "$MolString\nM  END\n";
  }
  else {
    # Just write out an empty molecule data string...
    print $NewSDFileRef SDFileUtil::GenerateEmptyCtabBlockLines(), "\n";
  }
}

# Write database data from SD or Text database file to SD output file...
#
sub WriteDatabaseDataToSDOutputFile {
  my($DatabaseCmpdID, $NewSDFileRef) = @_;

  if ($FingerprintsFilesInfo{Database}{OutputDataFields}) {
    my($DataFieldLabel, $DataFieldValue, @CmpdLines, %DataFieldLabelAndValues);

    @CmpdLines = split /\n/, $DatabaseFingerprintsFileData{$DatabaseCmpdID};
    %DataFieldLabelAndValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

    for $DataFieldLabel ($FingerprintsFilesInfo{Database}{OutputCurrentDataFields} ? GetCmpdDataHeaderLabels(\@CmpdLines) : @{$FingerprintsFilesInfo{Database}{DataFieldsToOutput}}) {
      $DataFieldValue = exists $DataFieldLabelAndValues{$DataFieldLabel} ? $DataFieldLabelAndValues{$DataFieldLabel} : '';
      print $NewSDFileRef  ">  <$DataFieldLabel>\n$DataFieldValue\n\n";
    }
  }
  elsif ($FingerprintsFilesInfo{Database}{OutputDataCols}) {
    my($DataColNum, $DataFieldLabel, $DataFieldValue);

    for $DataColNum (@{$FingerprintsFilesInfo{Database}{DataColNumsToOutput}}) {
      $DataFieldLabel = $FingerprintsFilesInfo{Database}{DataColNumToLabelMap}{$DataColNum};
      $DataFieldValue =  $DatabaseFingerprintsFileData{$DatabaseCmpdID}->[$DataColNum];
      print $NewSDFileRef  ">  <$DataFieldLabel>\n$DataFieldValue\n\n";
    }
  }
}

# Retriebe database data from SD or Text database file for text output file...
#
sub RetrieveDatabaseDataForTextOutputFile {
  my($DatabaseCmpdID) = @_;

  if ($FingerprintsFilesInfo{Database}{OutputDataFields}) {
    my(@CmpdLines, %DataFieldLabelAndValues);

    @CmpdLines = split /\n/, $DatabaseFingerprintsFileData{$DatabaseCmpdID};
    %DataFieldLabelAndValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

    return map { exists $DataFieldLabelAndValues{$_} ? $DataFieldLabelAndValues{$_} : ''} @{$FingerprintsFilesInfo{Database}{DataFieldsToOutput}};
  }
  elsif ($FingerprintsFilesInfo{Database}{OutputDataCols}) {
    if (exists $DatabaseFingerprintsFileData{$DatabaseCmpdID}) {
      return map { $DatabaseFingerprintsFileData{$DatabaseCmpdID}->[$_] } (0 .. $#{$FingerprintsFilesInfo{Database}{DataColNumsToOutput}});
    }
    else {
      return ('') x $#{$FingerprintsFilesInfo{Database}{DataColNumsToOutput}};
    }
  }
}

# Collect database file SD compound string or CSV/TSV data line for generating results
# files..
#
sub CollectDatabaseFileData {
  my($DatabaseCmpdID, $DatabaseFingerprintsFileIO) = @_;

  if (exists $DatabaseFingerprintsFileData{$DatabaseCmpdID}) {
    return;
  }

  if ($FingerprintsFilesInfo{Database}{CollectCmpdStringData}) {
    $DatabaseFingerprintsFileData{$DatabaseCmpdID} = $DatabaseFingerprintsFileIO->GetCompoundString();
  }

  if ($FingerprintsFilesInfo{Database}{CollectDataLine}) {
    my(@DataLineWords);
    @DataLineWords = $DatabaseFingerprintsFileIO->GetDataLineWords();
    $DatabaseFingerprintsFileData{$DatabaseCmpdID} = \@DataLineWords;
  }

}

# Read fingerprints data from reference fingerprints file...
#
sub ReadReferenceFingerprintsData {
  my($FingerprintsFileIO);

  $FingerprintsFileIO = Fingerprints::FingerprintsFileUtil::NewFingerprintsFileIO(%{$FingerprintsFilesInfo{Reference}{FingerprintsFileIOParameters}});
  ($SimilaritySearchInfo{ReferenceCmpdIDsRef}, $SimilaritySearchInfo{ReferenceFingerprintsObjectsRef}) = Fingerprints::FingerprintsFileUtil::ReadAndProcessFingerpritsData($FingerprintsFileIO);

}

# Retrieve information about fingerprints files...
#
sub RetrieveFingerprintsFilesInfo {

  %FingerprintsFilesInfo = ();
  %OutputFilesInfo = ();
  %SimilaritySearchInfo = ();

  %{$FingerprintsFilesInfo{Reference}} = ();
  %{$FingerprintsFilesInfo{Database}} = ();

  # Set up reference and database file names...
  $FingerprintsFilesInfo{Reference}{FileName} = $FingerprintsFilesList[0];
  $FingerprintsFilesInfo{Database}{FileName} = $FingerprintsFilesList[1];

  # Retrieve information about reference and database fingerprints file...
  RetrieveReferenceFingerprintsFileInfo();
  RetrieveDatabaseFingerprintsFileInfo();

  # Setup fingerprints comparison method and associated method parameters...
  SetupReferenceAndDatabaseFingerprintsComparisonInfo();

  # Retrieve information for output files...
  RetrieveOutputFilesInfo();
}

# Setup refrerence and database fingerprints comparison method and associated method parameters...
#
sub SetupReferenceAndDatabaseFingerprintsComparisonInfo {

  # Make sure reference and database fingerprints string match...
  if (($FingerprintsFilesInfo{Reference}{FirstFingerprintsStringType} !~ /^$FingerprintsFilesInfo{Database}{FirstFingerprintsStringType}$/i) ||
     ($FingerprintsFilesInfo{Reference}{FingerprintsBitVectorStringMode} != $FingerprintsFilesInfo{Database}{FingerprintsBitVectorStringMode}) ||
     ($FingerprintsFilesInfo{Reference}{FingerprintsVectorStringMode} != $FingerprintsFilesInfo{Database}{FingerprintsVectorStringMode}) ) {
    die "Error: First reference fingerprints string type, $FingerprintsFilesInfo{Reference}{FirstFingerprintsStringType}, must match first database fingerprints type, $FingerprintsFilesInfo{Database}{FirstFingerprintsStringType}.\n";
  }

  if ($FingerprintsFilesInfo{Reference}{FirstFingerprintsStringDescription} !~ /^$FingerprintsFilesInfo{Database}{FirstFingerprintsStringDescription}$/i) {
    warn "Warning: First reference fingerprints string description, $FingerprintsFilesInfo{Reference}{FirstFingerprintsStringDescription}, doesn't match first database fingerprints string description, $FingerprintsFilesInfo{Database}{FirstFingerprintsStringDescription}.\n";
  }

  # Setup individual reference and multiple references search mode...
  $SimilaritySearchInfo{IndividualReferenceMode} = undef;
  $SimilaritySearchInfo{MultipleReferencesMode} = undef;

  if ($OptionsInfo{Mode} =~ /^IndividualReference$/i) {
    $SimilaritySearchInfo{IndividualReferenceMode} = 1;
  }
  elsif ($OptionsInfo{Mode} =~ /^MultipleReferences$/i) {
    $SimilaritySearchInfo{MultipleReferencesMode} = 1;
  }
  else {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: IndividualReference, MultipleReferences\n";
  }

  # Set up reference and database fingerprints similarity search method and paramaters...
  my($ComparisonMeasure, $ComparisonMethod, $ApplyComparisonCutoff, $ComparisonCutoff, $KeepTop, @ComparisonMethodParameters);

  $SimilaritySearchInfo{ComparisonMethod} = '';
  @{$SimilaritySearchInfo{ComparisonMethodParameters}} = ();

  $SimilaritySearchInfo{ComparisonCutoff} = '';
  $SimilaritySearchInfo{KeepTop} = '';

  $ComparisonMeasure = ''; $ComparisonMethod = '';
  @ComparisonMethodParameters = ();

  FINGERPRINTSTYPE: {
    if ($FingerprintsFilesInfo{Reference}{FingerprintsBitVectorStringMode}) {
      $ComparisonMeasure = $OptionsInfo{SpecifiedBitVectorComparisonMeasure};
      $ComparisonMethod = $OptionsInfo{SpecifiedBitVectorComparisonMeasureMethod};

      if ($ComparisonMeasure =~ /^TverskySimilarity$/i) {
	push @ComparisonMethodParameters, $OptionsInfo{Alpha};
      }
      elsif ($ComparisonMeasure =~ /^WeightedTverskySimilarity$/i) {
	push @ComparisonMethodParameters, $OptionsInfo{Alpha};
	push @ComparisonMethodParameters, $OptionsInfo{Beta};
      }
      elsif ($ComparisonMeasure =~ /^WeightedTanimotoSimilarity$/i) {
	push @ComparisonMethodParameters, $OptionsInfo{Beta};
      }

      last FINGERPRINTSTYPE;
    }
    if ($FingerprintsFilesInfo{Reference}{FingerprintsVectorStringMode}) {
      my($SkipValuesCheck);

      $ComparisonMeasure = $OptionsInfo{SpecifiedVectorComparisonMeasure};
      $ComparisonMethod = $OptionsInfo{SpecifiedVectorComparisonMeasuresMethod};

      push @ComparisonMethodParameters, $OptionsInfo{SpecifiedVectorComparisonMode};

      $SkipValuesCheck = $OptionsInfo{Fast} ? 1 : 0;
      push @ComparisonMethodParameters, $SkipValuesCheck;

      last FINGERPRINTSTYPE;
    }
    die "Error: Uknown fingerprints string type. Supported values: FingerprintsBitVectorString or FingerprintsVectorString.\n";
  }

  $ApplyComparisonCutoff = $SimilaritySearchInfo{IndividualReferenceMode} ? 1 : (($SimilaritySearchInfo{MultipleReferencesMode} && $OptionsInfo{GroupFusionApplyCutoff}) ? 1 : 0);

  $ComparisonCutoff = ''; $KeepTop = '';
  if ($ComparisonMethod =~ /Distance/i) {
    $ComparisonCutoff = $OptionsInfo{DistanceCutoff};
    $KeepTop = ($OptionsInfo{SearchMode} =~ /^SimilaritySearch$/i) ? 0 : 1;
  }
  else {
    $ComparisonCutoff = $OptionsInfo{SimilarityCutoff};
    $KeepTop = ($OptionsInfo{SearchMode} =~ /^SimilaritySearch$/i) ? 1 : 0;
  }

  $SimilaritySearchInfo{ComparisonMethod} = $ComparisonMethod;
  @{$SimilaritySearchInfo{ComparisonMethodParameters}} = @ComparisonMethodParameters;

  $SimilaritySearchInfo{ComparisonCutoff} = $ComparisonCutoff;
  $SimilaritySearchInfo{KeepTop} = $KeepTop;
  $SimilaritySearchInfo{ApplyComparisonCutoff} = $ApplyComparisonCutoff;

  # Setup references to group fusion methods...
  $SimilaritySearchInfo{GroupFusionMethodRef} = undef;
  $SimilaritySearchInfo{ApplyPrecisionDuringFusion} = undef;

  FUSIONRULE: {
    if ($OptionsInfo{GroupFusionRule} =~ /^Max$/i) {
      # It's always the first value in the appropriated sorted list using value of KeepTop...
      $SimilaritySearchInfo{GroupFusionMethodRef} = sub { my($ComparisonValuesRef) = @_; return $ComparisonValuesRef->[0]; };
      last FUSIONRULE;
    }
    if ($OptionsInfo{GroupFusionRule} =~ /^Min$/i) {
      # It's always the last value in the appropriated sorted list using value of KeepTop...
      $SimilaritySearchInfo{GroupFusionMethodRef} = sub { my($ComparisonValuesRef) = @_; return $ComparisonValuesRef->[$#{$ComparisonValuesRef}]; };
      last FUSIONRULE;
    }
    if ($OptionsInfo{GroupFusionRule} =~ /^Mean$/i) {
      $SimilaritySearchInfo{GroupFusionMethodRef} = \&StatisticsUtil::Mean;
      $SimilaritySearchInfo{ApplyPrecisionDuringFusion} = 1;
      last FUSIONRULE;
    }
    if ($OptionsInfo{GroupFusionRule} =~ /^Median$/i) {
      $SimilaritySearchInfo{GroupFusionMethodRef} = \&StatisticsUtil::Median;
      $SimilaritySearchInfo{ApplyPrecisionDuringFusion} = 1;
      last FUSIONRULE;
    }
    if ($OptionsInfo{GroupFusionRule} =~ /^Sum$/i) {
      $SimilaritySearchInfo{GroupFusionMethodRef} = \&StatisticsUtil::Sum;
      $SimilaritySearchInfo{ApplyPrecisionDuringFusion} = 1;
      last FUSIONRULE;
    }
    if ($OptionsInfo{GroupFusionRule} =~ /^Euclidean$/i) {
      $SimilaritySearchInfo{GroupFusionMethodRef} = \&StatisticsUtil::Euclidean;
      $SimilaritySearchInfo{ApplyPrecisionDuringFusion} = 1;
      last FUSIONRULE;
    }
    die "Error: The value specified, $Options{groupfusionrule}, for option \"-g, --GroupFusionRule\" is not valid. Allowed values: Max, Min, Mean, Median, Sum, Euclidean\n";
  }

  $SimilaritySearchInfo{UsekNN} = ($OptionsInfo{kNN} !~ /^All$/i) ? 1 : 0;
  $SimilaritySearchInfo{SortComparisonValues} = (($OptionsInfo{GroupFusionRule} =~ /^(Max|Min)$/i) || $SimilaritySearchInfo{UsekNN}) ? 1 : 0;
}

# Retrieve information about reference fingerprints file...
#
sub RetrieveReferenceFingerprintsFileInfo {
  my($FingerprintsFile, $FileType, $InDelim, $FingerprintsFileIO, $FingerprintsStringMode, $FingerprintsBitVectorStringMode, $FingerprintsVectorStringMode, $FirstFingerprintsStringType, $FirstFingerprintsStringDescription);

  $FingerprintsFile = $FingerprintsFilesInfo{Reference}{FileName};
  ($FileType, $InDelim) =  RetrieveFingerprintsFileInfo($FingerprintsFile);

  $FingerprintsFilesInfo{Reference}{FileType} = $FileType;
  $FingerprintsFilesInfo{Reference}{InDelim} = $InDelim;

  # Setup reference FingerprintsFileIO parameters...
  %{$FingerprintsFilesInfo{Reference}{FingerprintsFileIOParameters}} = RetrieveFingerprintsFileIOParameters('Reference', $FileType, $FingerprintsFile);

  # Make sure reference fingerprints data file contains valid and retrieve fingerprints string mode information...
  ($FingerprintsStringMode, $FingerprintsBitVectorStringMode, $FingerprintsVectorStringMode, $FirstFingerprintsStringType, $FirstFingerprintsStringDescription) = RetrieveFingerprintsFileFingerprintsStringInfo('Reference', $FingerprintsFile);
  $FingerprintsFilesInfo{Reference}{FingerprintsStringMode} = $FingerprintsStringMode;
  $FingerprintsFilesInfo{Reference}{FingerprintsBitVectorStringMode} = $FingerprintsBitVectorStringMode;
  $FingerprintsFilesInfo{Reference}{FingerprintsVectorStringMode} = $FingerprintsVectorStringMode;
  $FingerprintsFilesInfo{Reference}{FirstFingerprintsStringType} = $FirstFingerprintsStringType;
  $FingerprintsFilesInfo{Reference}{FirstFingerprintsStringDescription} = $FirstFingerprintsStringDescription;

}

# Retrieve information about database fingerprints file...
#
sub RetrieveDatabaseFingerprintsFileInfo {
  my($FingerprintsFile, $FileType, $InDelim, $FingerprintsFileIO, $FingerprintsStringMode, $FingerprintsBitVectorStringMode, $FingerprintsVectorStringMode, $FirstFingerprintsStringType, $FirstFingerprintsStringDescription);

  $FingerprintsFile = $FingerprintsFilesInfo{Database}{FileName};
  ($FileType, $InDelim) =  RetrieveFingerprintsFileInfo($FingerprintsFile);

  $FingerprintsFilesInfo{Database}{FileType} = $FileType;
  $FingerprintsFilesInfo{Database}{InDelim} = $InDelim;

  # Setup reference FingerprintsFileIO parameters...
  %{$FingerprintsFilesInfo{Database}{FingerprintsFileIOParameters}} = RetrieveFingerprintsFileIOParameters('Database', $FileType, $FingerprintsFile);

  # Make sure database fingerprints data file contains valid and retrieve fingerprints string mode information...
  ($FingerprintsStringMode, $FingerprintsBitVectorStringMode, $FingerprintsVectorStringMode, $FirstFingerprintsStringType, $FirstFingerprintsStringDescription) = RetrieveFingerprintsFileFingerprintsStringInfo('Database', $FingerprintsFile);
  $FingerprintsFilesInfo{Database}{FingerprintsStringMode} = $FingerprintsStringMode;
  $FingerprintsFilesInfo{Database}{FingerprintsBitVectorStringMode} = $FingerprintsBitVectorStringMode;
  $FingerprintsFilesInfo{Database}{FingerprintsVectorStringMode} = $FingerprintsVectorStringMode;
  $FingerprintsFilesInfo{Database}{FirstFingerprintsStringType} = $FirstFingerprintsStringType;
  $FingerprintsFilesInfo{Database}{FirstFingerprintsStringDescription} = $FirstFingerprintsStringDescription;

  # Retrieve database fingerprints data field information for output file...
  #
  RetrieveDatabaseFingerprintsDataFieldsInfo($FingerprintsFile, $FileType, $InDelim);

  # Retrieve database fingerprints text file data columns information for output file...
  #
  RetrieveDatabaseFingerprintsDataColsInfo($FingerprintsFile, $FileType, $InDelim);

  # Any need to collect database compound string or data line for generation of results files...
  $FingerprintsFilesInfo{Database}{CollectCmpdStringData} = ($FileType =~ /^SD$/i) ? 1 : 0;
  $FingerprintsFilesInfo{Database}{CollectDataLine} = ($FileType =~ /^Text$/i && $OptionsInfo{DatabaseDataColsMode} =~ /^(All|Specify)$/i) ? 1 : 0;
  $FingerprintsFilesInfo{Database}{CollectInputFileData} = ($FingerprintsFilesInfo{Database}{CollectCmpdStringData} || $FingerprintsFilesInfo{Database}{CollectDataLine}) ? 1 : 0;

  # Set maximum number of similar compounds to find for individual reference of set of multiple
  # reference compounds...
  #
  SetMaximumSimilarMoleculesToRetrieve($FingerprintsFile, $FileType, $InDelim);
}

# Retrieve database fingerprints data field information...
#
sub RetrieveDatabaseFingerprintsDataFieldsInfo {
  my($FingerprintsFile, $FileType, $InDelim) = @_;
  my($CollectDataFields, $CmpdCount, $AllDataFieldsRef, $CommonDataFieldsRef, @DataFieldsToOutput);

  $FingerprintsFilesInfo{Database}{OutputDataFields} = 0;
  @{$FingerprintsFilesInfo{Database}{DataFieldsToOutput}} = ();

  $FingerprintsFilesInfo{Database}{OutputCurrentDataFields} = 0;

  @{$FingerprintsFilesInfo{Database}{AllDataFields}} = ();
  @{$FingerprintsFilesInfo{Database}{CommonDataFields}} = ();
  @{$FingerprintsFilesInfo{Database}{SpecifiedDatabaseDataFields}} = ();

  if ($FileType !~ /^SD$/i) {
    return;
  }

  # No need to go over SD file and collect data fields for SD file during All DatabaseDataFieldsMode as
  # they would be retrieved from database SD file compound string during generation of output files...
  #
  $CollectDataFields = (($OptionsInfo{TextOutput} && $OptionsInfo{DatabaseDataFieldsMode} =~ /^(All|Common)$/i) || ($OptionsInfo{SDOutput} && $OptionsInfo{DatabaseDataFieldsMode} =~ /^Common$/i)) ? 1 : 0;

  ($CmpdCount, $AllDataFieldsRef, $CommonDataFieldsRef) = (undef) x 2;

  if ($CollectDataFields) {
    open SDFILE, "$FingerprintsFile" or die "Error: Couldn't open $FingerprintsFile: $! \n";
    ($CmpdCount, $AllDataFieldsRef, $CommonDataFieldsRef) = GetAllAndCommonCmpdDataHeaderLabels(\*SDFILE);
    close SDFILE;
  }

  @DataFieldsToOutput = ();
  if ($OptionsInfo{DatabaseDataFieldsMode} =~ /^All$/i) {
    if (defined $AllDataFieldsRef) {
      push @DataFieldsToOutput, @{$AllDataFieldsRef};
      push @{$FingerprintsFilesInfo{Database}{AllDataFields}}, @{$AllDataFieldsRef};
    }
    else {
      # Retrieve and output data fields and values dynamically...
      $FingerprintsFilesInfo{Database}{OutputCurrentDataFields} = 1;
    }
  }
  elsif ($OptionsInfo{DatabaseDataFieldsMode} =~ /^Common$/i) {
    if (defined $CommonDataFieldsRef) {
      push @DataFieldsToOutput, @{$CommonDataFieldsRef};
      push @{$FingerprintsFilesInfo{Database}{CommonDataFields}}, @{$CommonDataFieldsRef};
    }
  }
  elsif ($OptionsInfo{DatabaseDataFieldsMode} =~ /^Specify$/i) {
    push @DataFieldsToOutput, @{$OptionsInfo{SpecifiedDatabaseDataFields}};
    push @{$FingerprintsFilesInfo{Database}{SpecifiedDatabaseDataFields}}, @{$OptionsInfo{SpecifiedDatabaseDataFields}};
  }

  if ($OptionsInfo{DatabaseDataFieldsMode} !~ /^CompoundID$/i) {
    $FingerprintsFilesInfo{Database}{OutputDataFields} = 1;
  }

  push @{$FingerprintsFilesInfo{Database}{DataFieldsToOutput}}, @DataFieldsToOutput;

}

# Retrieve database fingerprints data columns information...
#
sub RetrieveDatabaseFingerprintsDataColsInfo {
  my($FingerprintsFile, $FileType, $InDelim) = @_;
  my($Line, $ColNum, $ColLabel, $NumOfCols, @DataColLabels, @DataColLabelsToOutput, @DataColNumsToOutput, %DataColLabelToNumMap, %DataColNumToLabelMap);

  $FingerprintsFilesInfo{Database}{OutputDataCols} = 0;

  @{$FingerprintsFilesInfo{Database}{DataColLabels}} = ();
  %{$FingerprintsFilesInfo{Database}{DataColLabelToNumMap}} = ();
  %{$FingerprintsFilesInfo{Database}{DataColNumToLabelMap}} = ();

  @{$FingerprintsFilesInfo{Database}{DataColNumsToOutput}} = ();
  @{$FingerprintsFilesInfo{Database}{DataColLabelsToOutput}} = ();

  if ($FileType !~ /^Text$/i) {
    return;
  }

  @DataColLabels = ();
  @DataColLabelsToOutput = ();
  @DataColNumsToOutput = ();

  %DataColLabelToNumMap = ();
  %DataColNumToLabelMap = ();

  # Get column label line...
  open TEXTFILE, "$FingerprintsFile" or die "Error: Couldn't open $FingerprintsFile: $! \n";
  $Line = TextUtil::GetTextLine(\*TEXTFILE);
  close TEXTFILE;

  $InDelim = ($InDelim =~ /^Tab$/i) ? "\t" : ($InDelim =~ /semicolon/i ? "\;" : "\,");

  @DataColLabels = TextUtil::SplitWords($Line, $InDelim);
  $NumOfCols = scalar @DataColLabels;

  for $ColNum (0 .. $#DataColLabels) {
    $ColLabel = $DataColLabels[$ColNum];
    $DataColLabelToNumMap{$ColLabel} = $ColNum;
    $DataColNumToLabelMap{$ColNum} = $ColLabel;
  }

  if ($OptionsInfo{DatabaseDataColsMode} =~ /^Specify$/i) {
    if ($OptionsInfo{DatabaseColMode} =~ /^ColNum$/i) {
      for $ColNum (@{$OptionsInfo{SpecifiedDatabaseDataCols}}) {
	if ($ColNum > $NumOfCols) {
	  die "Error: Column number, $ColNum, specified using \"--DatabaseDataCols\" is not valid: It must be <= $NumOfCols\n";
	}
	push @DataColNumsToOutput, ($ColNum - 1);
      }
    }
    elsif ($OptionsInfo{DatabaseColMode} =~ /^ColLabel$/i) {
      for $ColLabel (@{$OptionsInfo{SpecifiedDatabaseDataCols}}) {
	if (!exists $DataColLabelToNumMap{$ColLabel}) {
	  die "Error: Column label, $ColLabel, specified using \"--DatabaseDataCols\" is not valid: It doesn't exist\n";
	}
	push @DataColNumsToOutput, $DataColLabelToNumMap{$ColLabel};
      }
    }
  }
  elsif ($OptionsInfo{DatabaseDataColsMode} =~ /^All$/i) {
    @DataColNumsToOutput = map { $_ } (0 .. $#DataColLabels);
  }

  # Setup data column labels to output...
  if (scalar @DataColNumsToOutput) {
    @DataColLabelsToOutput = map { $DataColNumToLabelMap{$_} } (0 .. $#DataColNumsToOutput);
  }

  $FingerprintsFilesInfo{Database}{OutputDataCols} = scalar @DataColNumsToOutput ? 1 : 0;

  @{$FingerprintsFilesInfo{Database}{DataColLabels}} = @DataColLabels;
  %{$FingerprintsFilesInfo{Database}{DataColLabelToNumMap}} = %DataColLabelToNumMap;
  %{$FingerprintsFilesInfo{Database}{DataColNumToLabelMap}} = %DataColNumToLabelMap;

  @{$FingerprintsFilesInfo{Database}{DataColNumsToOutput}} = @DataColNumsToOutput;
  @{$FingerprintsFilesInfo{Database}{DataColLabelsToOutput}} = @DataColLabelsToOutput;
}

# Set maximum number of similar compounds to find for individual reference of set of multiple
# reference compounds...
#
sub SetMaximumSimilarMoleculesToRetrieve {
  my($FingerprintsFile, $FileType, $InDelim) = @_;
  my($MaxSimilarMolecules, $NumOfDatabaseMolecules, $PercentSimilarMolecules, $Line);

  if ($OptionsInfo{SimilarCountMode} !~ /^PercentSimilar$/i) {
    return;
  }

  $PercentSimilarMolecules = $OptionsInfo{PercentSimilarMolecules};

  # Count database entries to figure out MaxSimilarMolecules using PercentSimilarMolecules
  # value...
  $NumOfDatabaseMolecules = 0;
  if ($FileType =~ /^SD$/i && exists($FingerprintsFilesInfo{Database}{NumOfDatabaseMolecules})) {
    # It might already be counted for SD file...
    $NumOfDatabaseMolecules = $FingerprintsFilesInfo{Database}{NumOfDatabaseMolecules};
  }
  else {
    print "Calculating maximum number of similar molecules to retrieve for \"PercentSimilar\" value of \"--SimilarCountMode\" option by counting number of molecules in database fingerprints file...\n";
    open FINGERPRINTSFILE, "$FingerprintsFile" or die "Error: Couldn't open $FingerprintsFile: $! \n";
    FILETYPE: {
      if ($FileType =~ /^SD$/i) {
	while ($Line = TextUtil::GetTextLine(\*FINGERPRINTSFILE)) {
	  if ($Line =~ /^\$\$\$\$/) {
	    $NumOfDatabaseMolecules++;
	  }
	}
	last FILETYPE;
      }
      if ($FileType =~ /^Text$/i) {
	# Ignore column label line...
	$Line = TextUtil::GetTextLine(\*FINGERPRINTSFILE);
	while ($Line = TextUtil::GetTextLine(\*FINGERPRINTSFILE)) {
	  $NumOfDatabaseMolecules++;
	}
	last FILETYPE;
      }
      if ($FileType =~ /^FP$/i) {
	while ($Line = TextUtil::GetTextLine(\*FINGERPRINTSFILE)) {
	  if ($Line !~ /^#/) {
	    $NumOfDatabaseMolecules++;
	  }
	}
	last FILETYPE;
      }
      $NumOfDatabaseMolecules = 0;
    }
    close FINGERPRINTSFILE;
    $FingerprintsFilesInfo{Database}{NumOfDatabaseMolecules} = $NumOfDatabaseMolecules;
  }

  $MaxSimilarMolecules = int (($NumOfDatabaseMolecules * $PercentSimilarMolecules)/100);
  if ($MaxSimilarMolecules < 1) {
    $MaxSimilarMolecules = 1;
  }

  $OptionsInfo{MaxSimilarMolecules} = $MaxSimilarMolecules;
}

# Retrieve information about fingerprints file...
#
sub RetrieveFingerprintsFileInfo {
  my($FingerprintsFile) = @_;
  my($FileType, $InDelim, $FileDir, $FileExt, $FileName);

  if (!(-e $FingerprintsFile)) {
    die "Error: Input fingerprints file, $FingerprintsFile, doesn't exist.\n";
  }

  $FileType = Fingerprints::FingerprintsFileUtil::GetFingerprintsFileType($FingerprintsFile);
  if (IsEmpty($FileType)) {
    die "Error: Input file, $FingerprintsFile, is not a fingerprints file.\n";
  }

  $InDelim = '';
  if ($FileType =~ /^Text$/i) {
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($FingerprintsFile);
    $InDelim = ($FileExt =~ /^tsv$/i) ? 'Tab' : $OptionsInfo{InDelim};
  }

  return ($FileType, $InDelim);
}

# Retrieve fingerprints file IO parameters...
#
sub RetrieveFingerprintsFileIOParameters {
  my($FingerprintsFileMode, $FileType, $FingerprintsFile) = @_;
  my(%FingerprintsFileIOParams);

  if ($FingerprintsFileMode !~ /^(Reference|Database)$/) {
    die "Error: Unknown fingerprints file mode: $FingerprintsFileMode. Supported values: Reference or Database\n";
  }

  %FingerprintsFileIOParams = ();

  FILETYPE: {
    if ($FileType =~ /^SD$/i) {
      %FingerprintsFileIOParams = ('Name' => $FingerprintsFile, 'Mode' => 'Read', 'FingerprintsStringMode' => $OptionsInfo{FingerprintsMode}, 'ValidateData' => $OptionsInfo{ValidateData}, 'DetailLevel' =>  $OptionsInfo{Detail}, 'FingerprintsFieldLabel' => $OptionsInfo{"${FingerprintsFileMode}FingerprintsField"}, 'CompoundIDMode' => $OptionsInfo{"${FingerprintsFileMode}CompoundIDMode"}, 'CompoundIDFieldLabel' => $OptionsInfo{"${FingerprintsFileMode}CompoundIDField"}, 'CompoundIDPrefix' => $OptionsInfo{"${FingerprintsFileMode}CompoundIDPrefix"});
      last FILETYPE;
    }
    if ($FileType =~ /^FP$/i) {
      %FingerprintsFileIOParams = ('Name' => $FingerprintsFile, 'Mode' => 'Read', 'FingerprintsStringMode' => $OptionsInfo{FingerprintsMode}, 'ValidateData' => $OptionsInfo{ValidateData}, 'DetailLevel' =>  $OptionsInfo{Detail});
      last FILETYPE;
    }
    if ($FileType =~ /^Text$/i) {
      %FingerprintsFileIOParams = ('Name' => $FingerprintsFile, 'Mode' => 'Read', 'FingerprintsStringMode' => $OptionsInfo{FingerprintsMode}, 'ValidateData' => $OptionsInfo{ValidateData}, 'DetailLevel' =>  $OptionsInfo{Detail}, 'FingerprintsCol' => $OptionsInfo{"${FingerprintsFileMode}FingerprintsCol"}, 'ColMode' => $OptionsInfo{"${FingerprintsFileMode}ColMode"}, 'CompoundIDCol' => $OptionsInfo{"${FingerprintsFileMode}CompoundIDCol"}, 'CompoundIDPrefix' => $OptionsInfo{"${FingerprintsFileMode}CompoundIDPrefix"}, 'InDelim' => $FingerprintsFilesInfo{$FingerprintsFileMode}{InDelim});
      last FILETYPE;
    }
    die "Error: Fingerprints file type, $FileType, is not valid. Supported file types: SD, FP or Text\n";
  }

  return %FingerprintsFileIOParams;
}

# Make sure fingerprints data file contains valid dta and retrieve fingerprints string mode information...
#
sub RetrieveFingerprintsFileFingerprintsStringInfo {
  my($FingerprintsFileMode, $FingerprintsFile) = @_;
  my($FingerprintsFileIO, $FingerprintsStringMode, $FingerprintsBitVectorStringMode, $FingerprintsVectorStringMode, $FirstFingerprintsStringType, $FirstFingerprintsStringDescription);

  $FingerprintsFileIO = Fingerprints::FingerprintsFileUtil::NewFingerprintsFileIO(%{$FingerprintsFilesInfo{$FingerprintsFileMode}{FingerprintsFileIOParameters}});
  if (!$FingerprintsFileIO) {
    die "Error: Reference fingerprints file, $FingerprintsFile, contains invalid fingerprints data.\n";
  }
  if (!$FingerprintsFileIO->IsFingerprintsFileDataValid()) {
    die "Error: Reference fingerprints file, $FingerprintsFile, contains invalid fingerprints data.\n";
  }

  $FingerprintsStringMode = $FingerprintsFileIO->GetFingerprintsStringMode();
  $FingerprintsBitVectorStringMode = $FingerprintsFileIO->GetFingerprintsBitVectorStringMode();
  $FingerprintsVectorStringMode = $FingerprintsFileIO->GetFingerprintsVectorStringMode();

  $FirstFingerprintsStringType = $FingerprintsFileIO->GetFirstFingerprintsStringType();
  $FirstFingerprintsStringDescription = $FingerprintsFileIO->GetFirstFingerprintsStringDescription();

  $FingerprintsFileIO->Close();

  return ($FingerprintsStringMode, $FingerprintsBitVectorStringMode, $FingerprintsVectorStringMode, $FirstFingerprintsStringType, $FirstFingerprintsStringDescription);
}

# Retrieve output files names using reference fingerprints file name...
#
sub RetrieveOutputFilesInfo {
  my($FingerprintsFile, $FileDir, $FileExt, $FileName, $OutFileRoot, $SDOutFileName, $TextOutFileName, $SDOutFileExt, $TextOutFileExt, $ReferenceFileName, $DatabaseFileName);

  $OutputFilesInfo{OutFileRoot} = '';
  $OutputFilesInfo{SDOutFileName} = '';
  $OutputFilesInfo{TextOutFileName} = '';

  $FingerprintsFile = $FingerprintsFilesInfo{Reference}{FileName};

  $FileDir = ""; $FileName = ""; $FileExt = "";
  ($FileDir, $FileName, $FileExt) = ParseFileName($FingerprintsFile);

  $SDOutFileExt = "sdf";
  $TextOutFileExt = ($Options{outdelim} =~ /^tab$/i) ? "tsv" : "csv";

  if ($OptionsInfo{OutFileRoot}) {
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
    $OutFileRoot = "${FileName}SimilaritySearching";
  }

  $SDOutFileName = "${OutFileRoot}.${SDOutFileExt}";
  $TextOutFileName = "${OutFileRoot}.${TextOutFileExt}";

  $ReferenceFileName = $FingerprintsFilesInfo{Reference}{FileName};
  $DatabaseFileName = $FingerprintsFilesInfo{Database}{FileName};

  if ($OptionsInfo{SDOutput}) {
    if ($SDOutFileName =~ /^$ReferenceFileName$/i) {
      die "Error: Same output, $SDOutFileName, and reference input file names.\nSpecify a different name using \"-r --root\" option or use default name.\n";
    }
    if ($SDOutFileName =~ /^$DatabaseFileName$/i) {
      die "Error: Same output, $SDOutFileName, and database input file names.\nSpecify a different name using \"-r --root\" option or use default name.\n";
    }
  }

  if ($OptionsInfo{TextOutput}) {
    if ($TextOutFileName =~ /^$ReferenceFileName$/i) {
      die "Error: Same output, $TextOutFileName, and reference input file names.\nSpecify a different name using \"-r --root\" option or use default name.\n";
    }
    if ($TextOutFileName =~ /^$DatabaseFileName$/i) {
      die "Error: Same output, $TextOutFileName, and database input file names.\nSpecify a different name using \"-r --root\" option or use default name.\n";
    }
  }

  if (!$OptionsInfo{OverwriteFiles}) {
    if ($OptionsInfo{SDOutput}) {
      if (-e $SDOutFileName) {
	die "Error: The output file $SDOutFileName already exists.\n";
      }
    }
    if ($OptionsInfo{TextOutput}) {
      if (-e $TextOutFileName) {
	die "Error: The output file $TextOutFileName already exists.\n";
      }
    }
  }

  $OutputFilesInfo{OutFileRoot} = $OutFileRoot;
  $OutputFilesInfo{SDOutFileName} = $SDOutFileName;
  $OutputFilesInfo{TextOutFileName} = $TextOutFileName;

}

# Process input fingerprints file names...
#
sub ProcessFingerprintsFileNames {
  @FingerprintsFilesList = ();

  if (@ARGV != 2) {
    die GetUsageFromPod("$FindBin::Bin/$ScriptName");
  }

  # Reference fingerprints file name...
  push @FingerprintsFilesList, $ARGV[0];

  # Database fingerprints file name...
  push @FingerprintsFilesList, $ARGV[1];

}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};
  $OptionsInfo{FingerprintsMode} = $Options{fingerprintsmode};

  $OptionsInfo{SearchMode} = $Options{searchmode};

  ProcessBitVectorComparisonOptions();
  ProcessVectorComparisonOptions();

  $OptionsInfo{GroupFusionRule} = $Options{groupfusionrule};
  $OptionsInfo{GroupFusionApplyCutoff} = ($Options{groupfusionapplycutoff} =~ /^Yes$/i) ? 1 : 0;;

  $OptionsInfo{SimilarCountMode} = $Options{similarcountmode};
  $OptionsInfo{NumOfSimilarMolecules} = $Options{numofsimilarmolecules};
  $OptionsInfo{PercentSimilarMolecules} = $Options{percentsimilarmolecules};

  # Set MaxSimilarMolecules to NumOfSimilarMolecules. For PercentSimilar value of SimilarCountMode,
  # it'll be overwritten using number of entries in database fingerprints file and value of PercentSimilarMolecules...
  #
  $OptionsInfo{MaxSimilarMolecules} = $OptionsInfo{NumOfSimilarMolecules};

  $OptionsInfo{SimilarityCutoff} = $Options{similaritycutoff};
  $OptionsInfo{DistanceCutoff} = $Options{distancecutoff};

  $OptionsInfo{kNN} = $Options{knn};
  if ($Options{knn} !~ /^All$/i) {
    if (!IsPositiveInteger($Options{knn})) {
      die "Error: The value specified, $Options{knn}, for option \"-k, --KNN\" is not valid. Allowed values: > 0 \n";
    }
  }

  ProcessReferenceFingerprintsDataOptions();
  ProcessDatabaseFingerprintsDataOptions();

  $OptionsInfo{Detail} = $Options{detail};

  $OptionsInfo{InDelim} = $Options{indelim};
  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /tab/i ) ? "\t" : (($Options{outdelim} =~ /semicolon/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{Output} = $Options{output};
  $OptionsInfo{SDOutput} = ($Options{output} =~ /^(SD|Both)$/i) ? 1 : 0;
  $OptionsInfo{TextOutput} = ($Options{output} =~ /^(Text|Both)$/i) ? 1 : 0;

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

  # Setup similarity coefficient to use for calculating similarity matrices for bit vector strings...
  my($SpecifiedMeasure, $SpecifiedComparisonMeasureName, $SpecifiedComparisonMeasureMethod);

  $SpecifiedComparisonMeasureName = '';
  $SpecifiedComparisonMeasureMethod = '';

  $SpecifiedMeasure = $Options{bitvectorcomparisonmode};

  if (! exists $SupportedComparisonMeasuresMethodMap{lc($SpecifiedMeasure)} )  {
      die "Error: The value specified, $SpecifiedMeasure, for option \"-b --BitVectorComparisonMode\" is not valid.\nAllowed values:", JoinWords(\@SupportedComparisonMeasures, ", ", 0), "\n";
  }

  $SpecifiedComparisonMeasureMethod = $SupportedComparisonMeasuresMethodMap{lc($SpecifiedMeasure)};
  $SpecifiedComparisonMeasureName = $SupportedComparisonMeasuresNameMap{lc($SpecifiedMeasure)};

  $OptionsInfo{BitVectorComparisonMode} = $Options{bitvectorcomparisonmode};

  $OptionsInfo{SpecifiedBitVectorComparisonMeasure} = $SpecifiedMeasure;
  $OptionsInfo{SpecifiedBitVectorComparisonMeasureName} = $SpecifiedComparisonMeasureName;
  $OptionsInfo{SpecifiedBitVectorComparisonMeasureMethod} = $SpecifiedComparisonMeasureMethod;

  # Make sure valid alpha parameter is specified for Tversky calculation...
  $OptionsInfo{Alpha} = '';
  if ($SpecifiedMeasure =~ /^(TverskySimilarity|WeightedTverskySimilarity)$/i) {
    if (IsEmpty($Options{alpha})) {
      die "Error: You must specify a value for \"-a, --alpha\" option in \"TverskySimilarity or WeightedTverskySimilarity\" \"-m --mode\". \n";
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
  if ($SpecifiedMeasure =~ /^(WeightedTverskySimilarity|WeightedTanimotoSimilarity)$/i) {
    if (IsEmpty($Options{beta})) {
      die "Error: You must specify a value for \"-b, --beta\" option in \"WeightedTverskySimilarity or WeightedTanimotoSimilarity\" \"-m --mode\". \n";
    }
    my($Beta);
    $Beta = $Options{beta};
    if (!(IsFloat($Beta) && $Beta >=0 && $Beta <= 1)) {
      die "Error: The value specified, $Options{beta}, for option \"-b, --beta\" is not valid. Allowed values: >= 0 and <= 1\n";
    }
    $OptionsInfo{Beta} = $Beta;
  }
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
  my($SpecifiedMeasure, $SpecifiedComparisonMeasureName, $SpecifiedComparisonMeasureMethod);

  $SpecifiedComparisonMeasureName = '';
  $SpecifiedComparisonMeasureMethod = '';

  $SpecifiedMeasure = $Options{vectorcomparisonmode};
  $SpecifiedMeasure =~ s/ //g;

  if (! exists($SupportedComparisonMeasuresMethodMap{lc($SpecifiedMeasure)})) {
    die "Error: The value specified, $SpecifiedMeasure, for option \"-v --VectorComparisonMode\" is not valid.\nAllowed values:", JoinWords(\@SupportedComparisonMeasures, ", ", 0), "\n";
  }

  $SpecifiedComparisonMeasureMethod = $SupportedComparisonMeasuresMethodMap{lc($SpecifiedMeasure)};
  $SpecifiedComparisonMeasureName = $SupportedComparisonMeasuresNameMap{lc($SpecifiedMeasure)};

  $OptionsInfo{VectorComparisonMode} = $Options{vectorcomparisonmode};

  $OptionsInfo{SpecifiedVectorComparisonMeasure} = $SpecifiedMeasure;
  $OptionsInfo{SpecifiedVectorComparisonMeasuresName} = $SpecifiedComparisonMeasureName;
  $OptionsInfo{SpecifiedVectorComparisonMeasuresMethod} = $SpecifiedComparisonMeasureMethod;

  # Setup specified vector comparison calculation modes...
  my($SpecifiedFormulism);

  $SpecifiedFormulism = $Options{vectorcomparisonformulism};
  $SpecifiedFormulism =~ s/ //g;
  if ($SpecifiedFormulism !~ /^(AlgebraicForm|BinaryForm|SetTheoreticForm)$/i) {
    die "Error: The value specified, $SpecifiedFormulism, for option \"--VectorComparisonFormulism\" is not valid. Allowed values: AlgebraicForm, BinaryForm or SetTheoreticForm\n";
  }

  $OptionsInfo{VectorComparisonFormulism} = $Options{vectorcomparisonformulism};
  $OptionsInfo{SpecifiedVectorComparisonMode} = $SpecifiedFormulism;

}

# Process options related to data retrieval from reference fingerprints SD and CSV/TSV
# text files...
#
sub ProcessReferenceFingerprintsDataOptions {

  $OptionsInfo{ReferenceCompoundIDPrefix} = $Options{referencecompoundidprefix} ? $Options{referencecompoundidprefix} : 'Cmpd';

  # Compound ID and fingerprints column options for text files...

  $OptionsInfo{ReferenceColMode} = $Options{referencecolmode};

  if (IsNotEmpty($Options{referencecompoundidcol})) {
    if ($Options{referencecolmode} =~ /^ColNum$/i) {
      if (!IsPositiveInteger($Options{referencecompoundidcol})) {
	die "Error: Column value, $Options{referencecompoundidcol}, specified using \"--ReferenceCompoundIDCol\" is not valid: Allowed integer values: > 0\n";
      }
    }
    $OptionsInfo{ReferenceCompoundIDCol} = $Options{referencecompoundidcol};
  }
  else {
    $OptionsInfo{ReferenceCompoundIDCol} = 'AutoDetect';
  }

  if (IsNotEmpty($Options{referencefingerprintscol})) {
    if ($Options{referencecolmode} =~ /^ColNum$/i) {
      if (!IsPositiveInteger($Options{referencefingerprintscol})) {
	die "Error: Column value, $Options{referencefingerprintscol}, specified using \"--ReferenceFingerprintsCol\" is not valid: Allowed integer values: > 0\n";
      }
    }
    $OptionsInfo{ReferenceFingerprintsCol} = $Options{referencefingerprintscol};
  }
  else {
    $OptionsInfo{ReferenceFingerprintsCol} = 'AutoDetect';
  }

  if (IsNotEmpty($Options{referencecompoundidcol}) && IsNotEmpty($Options{referencefingerprintscol})) {
    if (IsPositiveInteger($Options{referencecompoundidcol}) && IsPositiveInteger($Options{referencefingerprintscol})) {
      if (($Options{referencecompoundidcol} == $Options{referencefingerprintscol})) {
	die "Error: Values specified using \"--ReferenceCompoundIDCol\" and \"--ReferenceFingerprintsCol\", $Options{referencecompoundidcol}, must be different.\n";
      }
    }
    else {
      if (($Options{referencecompoundidcol} eq $Options{referencefingerprintscol})) {
	die "Error: Values specified using \"--ReferenceCompoundIDCol\" and \"--ReferenceFingerprintsCol\", $Options{referencecompoundidcol}, must be different.\n";
      }
    }
  }

  # Compound ID and fingerprints field options for SD files...

  $OptionsInfo{ReferenceCompoundIDMode} = $Options{referencecompoundidmode};
  $OptionsInfo{ReferenceCompoundIDField} = '';

  if ($Options{referencecompoundidmode} =~ /^DataField$/i && !$Options{referencecompoundidfield}) {
    die "Error: You must specify a value for \"--ReferenceCompoundIDField\" option in \"DataField\" \"--ReferenceCompoundIDMode\". \n";
  }
  if ($Options{referencecompoundidfield}) {
    $OptionsInfo{ReferenceCompoundIDField} = $Options{referencecompoundidfield};
  }

  if (IsNotEmpty($Options{referencefingerprintsfield})) {
    $OptionsInfo{ReferenceFingerprintsField} = $Options{referencefingerprintsfield};
  }
  else {
    $OptionsInfo{ReferenceFingerprintsField} = 'AutoDetect';
  }

  if ($Options{referencecompoundidfield} && IsNotEmpty($Options{referencefingerprintsfield})) {
    if (($Options{referencecompoundidfield} eq $Options{referencefingerprintsfield})) {
      die "Error: Values specified using \"--ReferenceCompoundIDField\" and \"--ReferenceFingerprintsfield\", $Options{referencecompoundidfield}, must be different.\n";
    }
  }

}

# Process options related to data retrieval from database fingerprints SD and CSV/TSV
# text files...
#
sub ProcessDatabaseFingerprintsDataOptions {

  $OptionsInfo{DatabaseCompoundIDPrefix} = $Options{databasecompoundidprefix} ? $Options{databasecompoundidprefix} : 'Cmpd';

  # Compound ID and fingerprints column options for text files...

  $OptionsInfo{DatabaseColMode} = $Options{databasecolmode};

  if (IsNotEmpty($Options{databasecompoundidcol})) {
    if ($Options{databasecolmode} =~ /^ColNum$/i) {
      if (!IsPositiveInteger($Options{databasecompoundidcol})) {
	die "Error: Column value, $Options{databasecompoundidcol}, specified using \"--DatabaseCompoundIDCol\" is not valid: Allowed integer values: > 0\n";
      }
    }
    $OptionsInfo{DatabaseCompoundIDCol} = $Options{databasecompoundidcol};
  }
  else {
    $OptionsInfo{DatabaseCompoundIDCol} = 'AutoDetect';
  }

  if (IsNotEmpty($Options{databasefingerprintscol})) {
    if ($Options{databasecolmode} =~ /^ColNum$/i) {
      if (!IsPositiveInteger($Options{databasefingerprintscol})) {
	die "Error: Column value, $Options{databasefingerprintscol}, specified using \"--DatabaseFingerprintsCol\" is not valid: Allowed integer values: > 0\n";
      }
    }
    $OptionsInfo{DatabaseFingerprintsCol} = $Options{databasefingerprintscol};
  }
  else {
    $OptionsInfo{DatabaseFingerprintsCol} = 'AutoDetect';
  }

  if (IsNotEmpty($Options{databasecompoundidcol}) && IsNotEmpty($Options{databasefingerprintscol})) {
    if (IsPositiveInteger($Options{databasecompoundidcol}) && IsPositiveInteger($Options{databasefingerprintscol})) {
      if (($Options{databasecompoundidcol} == $Options{databasefingerprintscol})) {
	die "Error: Values specified using \"--DatabaseCompoundIDCol\" and \"--DatabaseFingerprintsCol\", $Options{databasecompoundidcol}, must be different.\n";
      }
    }
    else {
      if (($Options{databasecompoundidcol} eq $Options{databasefingerprintscol})) {
	die "Error: Values specified using \"--DatabaseCompoundIDCol\" and \"--DatabaseFingerprintsCol\", $Options{databasecompoundidcol}, must be different.\n";
      }
    }
  }

  # Database data column options for text files...

  $OptionsInfo{DatabaseDataColsMode} = $Options{databasedatacolsmode};
  $OptionsInfo{DatabaseDataCols} = '';
  @{$OptionsInfo{SpecifiedDatabaseDataCols}} = ();

  if ($Options{databasedatacolsmode} =~ /^Specify$/i) {
    my($DatabaseDataCols, $DatabaseColNum, @SpecifiedDataCols);

    if (!$Options{databasedatacols}) {
      die "Error: You must specify a value for \"--DatabaseDataCols\" option in \"Specify\" \"--DatabaseDataColsMode\". \n";
    }
    $DatabaseDataCols = $Options{databasedatacols};

    if ($Options{databasecolmode} =~ /^ColNum$/i) {
      $DatabaseDataCols =~ s/ //g;
      @SpecifiedDataCols = split /\,/, $DatabaseDataCols;
      for $DatabaseColNum (@SpecifiedDataCols) {
	if (!IsPositiveInteger($DatabaseColNum)) {
	  die "Error: Column value, $DatabaseColNum, specified using \"--DatabaseDataCols\" is not valid: Allowed integer values: > 0\n";
	}
      }
    }
    else {
      @SpecifiedDataCols = split /\,/, $DatabaseDataCols;
    }
    $OptionsInfo{DatabaseDataCols} = $DatabaseDataCols;
    push @{$OptionsInfo{SpecifiedDatabaseDataCols}}, @SpecifiedDataCols;
  }
  elsif ($Options{databasedatacolsmode} =~ /^All$/i) {
    $OptionsInfo{DatabaseDataCols} = 'All';
  }

  if ($OptionsInfo{DatabaseDataColsMode} =~ /^Specify$/i && !$OptionsInfo{DatabaseDataCols}) {
    die "Error: You must specify a value for \"--DatabaseDataCols\" option in \"Specify\" \"--DatabaseDataColsMode\". \n";
  }

  # Compound ID and fingerprints field options for SD files...

  $OptionsInfo{DatabaseCompoundIDMode} = $Options{databasecompoundidmode};
  $OptionsInfo{DatabaseCompoundIDField} = $Options{databasecompoundidfield} ? $Options{databasecompoundidfield} : '';

  if ($Options{databasecompoundidmode} =~ /^DataField$/i) {
    if (!$Options{databasecompoundidfield}) {
      die "Error: You must specify a value for \"--DatabaseCompoundIDField\" option in \"DataField\" \"--DatabaseCompoundIDMode\". \n";
    }
    $OptionsInfo{DatabaseCompoundIDField} = $Options{databasecompoundidfield};
  }


  if (IsNotEmpty($Options{databasefingerprintsfield})) {
    $OptionsInfo{DatabaseFingerprintsField} = $Options{databasefingerprintsfield};
  }
  else {
    $OptionsInfo{DatabaseFingerprintsField} = 'AutoDetect';
  }

  if ($Options{databasecompoundidfield} && IsNotEmpty($Options{databasefingerprintsfield})) {
    if (($Options{databasecompoundidfield} eq $Options{databasefingerprintsfield})) {
      die "Error: Values specified using \"--DatabaseCompoundIDField\" and \"--DatabaseFingerprintsfield\", $Options{databasecompoundidfield}, must be different.\n";
    }
  }

  # Database data field options for SD files...

  $OptionsInfo{DatabaseDataFieldsMode} = $Options{databasedatafieldsmode};
  $OptionsInfo{DatabaseDataFields} = '';
  @{$OptionsInfo{SpecifiedDatabaseDataFields}} = ();

  if ($Options{databasedatafieldsmode} =~ /^Specify$/i && !$Options{databasedatafields}) {
    die "Error: You must specify a value for \"--DatabaseDataFields\" option in \"Specify\" \"--DatabaseDataFieldsMode\". \n";
  }
  if ($Options{databasedatafields}) {
    my(@SpecifiedDataFields);
    $OptionsInfo{DatabaseDataFields} = $Options{databasedatafields};

    @SpecifiedDataFields = split /\,/, $Options{databasedatafields};
    push @{$OptionsInfo{SpecifiedDatabaseDataFields}}, @SpecifiedDataFields;
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{alpha} = 0.5;
  $Options{beta} = 1;

  $Options{bitvectorcomparisonmode} = "TanimotoSimilarity";

  $Options{databasecolmode} = 'colnum';

  $Options{databasecompoundidprefix} = 'Cmpd';
  $Options{databasecompoundidmode} = 'LabelPrefix';

  $Options{databasedatacolsmode} = 'CompoundID';
  $Options{databasedatafieldsmode} = 'CompoundID';

  $Options{distancecutoff} = 10;

  $Options{referencecolmode} = 'colnum';

  $Options{referencecompoundidprefix} = 'Cmpd';
  $Options{referencecompoundidmode} = 'LabelPrefix';

  $Options{detail} = 1;

  $Options{fingerprintsmode} = 'AutoDetect';
  $Options{groupfusionrule} = 'Max';
  $Options{groupfusionapplycutoff} = 'Yes';

  $Options{knn} = 'All';

  $Options{mode} = 'MultipleReferences';

  $Options{numofsimilarmolecules} = 10;
  $Options{percentsimilarmolecules} = 1;

  $Options{indelim} = 'comma';
  $Options{outdelim} = 'comma';
  $Options{quote} = 'yes';

  $Options{output} = 'text';

  $Options{precision} = 2;

  $Options{searchmode} = 'SimilaritySearch';

  $Options{similarcountmode} = 'NumOfSimilar';

  $Options{similaritycutoff} = 0.75;

  $Options{vectorcomparisonmode} = 'TanimotoSimilarity';
  $Options{vectorcomparisonformulism} = 'AlgebraicForm';

  if (!GetOptions(\%Options, "alpha=f", "beta=f", "bitvectorcomparisonmode|b=s", "databasecolmode=s", "databasecompoundidcol=s", "databasecompoundidprefix=s", "databasecompoundidfield=s", "databasecompoundidmode=s", "databasedatacols=s", "databasedatacolsmode=s", "databasedatafields=s", "databasedatafieldsmode=s", "databasefingerprintscol=s", "databasefingerprintsfield=s", "distancecutoff=f", "detail|d=i", "fast|f", "fingerprintsmode=s", "groupfusionrule|g=s", , "groupfusionapplycutoff=s", "help|h", "indelim=s", "knn|k=s", "mode|m=s", "numofsimilarmolecules|n=i", "outdelim=s", "output=s", "overwrite|o", "percentsimilarmolecules|p=f", "precision=s", "quote|q=s", "referencecolmode=s", "referencecompoundidcol=s", "referencecompoundidprefix=s", "referencecompoundidfield=s", "referencecompoundidmode=s", "referencefingerprintscol=s", "referencefingerprintsfield=s", "root|r=s", "searchmode|s=s", "similarcountmode=s", "similaritycutoff=f", "vectorcomparisonmode|v=s", "vectorcomparisonformulism=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{databasecolmode} !~ /^(ColNum|ColLabel)$/i) {
    die "Error: The value specified, $Options{databasecolmode}, for option \"--DatabaseColMode\" is not valid. Allowed values: ColNum, or ColLabel\n";
  }
  if ($Options{databasecompoundidmode} !~ /^(DataField|MolName|LabelPrefix|MolNameOrLabelPrefix)$/i) {
    die "Error: The value specified, $Options{databasecompoundidmode}, for option \"--DatabaseCompoundIDMode\" is not valid. Allowed values: DataField, MolName, LabelPrefix or MolNameOrLabelPrefix\n";
  }
  if ($Options{databasedatacolsmode} !~ /^(All|Specify|CompoundID)$/i) {
    die "Error: The value specified, $Options{databasedatacolsmode}, for option \"--DatabaseDataColsMode\" is not valid. Allowed values: All, Specify, or CompoundID\n";
  }
  if ($Options{databasedatafieldsmode} !~ /^(All|Common|Specify|CompoundID)$/i) {
    die "Error: The value specified, $Options{databasedatafieldsmode}, for option \"--DatabaseDataFieldsMode\" is not valid. Allowed values: All, Common, Specify, or CompoundID\n";
  }
  if (!IsPositiveInteger($Options{detail})) {
    die "Error: The value specified, $Options{detail}, for option \"-d, --detail\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{fingerprintsmode} !~ /^(AutoDetect|FingerprintsBitVectorString|FingerprintsVectorString)$/i) {
    die "Error: The value specified, $Options{fingerprintsmode}, for option \"--FingerprintsMode\" is not valid. Allowed values: AutoDetect, FingerprintsBitVectorString or FingerprintsVectorString \n";
  }
  if ($Options{groupfusionrule} !~ /^(Max|Min|Mean|Median|Sum|Euclidean)$/i) {
    die "Error: The value specified, $Options{groupfusionrule}, for option \"-g, --GroupFusionRule\" is not valid. Allowed values: Max, Min, Mean, Median, Sum, Euclidean\n";
  }
  if ($Options{groupfusionapplycutoff} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"--GroupFusionApplyCutoff\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{indelim} !~ /^(comma|semicolon)$/i) {
    die "Error: The value specified, $Options{indelim}, for option \"--InDelim\" is not valid. Allowed values: comma, or semicolon\n";
  }
  if ($Options{mode} !~ /^(IndividualReference|MultipleReferences)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: IndividualReference, MultipleReferences\n";
  }
  if (!IsPositiveInteger($Options{numofsimilarmolecules})) {
    die "Error: The value specified, $Options{numofsimilarmolecules}, for option \"-n, --NumOfSimilarMolecules\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{outdelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{outdelim}, for option \"--OutDelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if ($Options{output} !~ /^(SD|text|both)$/i) {
    die "Error: The value specified, $Options{output}, for option \"--output\" is not valid. Allowed values: SD, text, or both\n";
  }
  if (!(IsFloat($Options{percentsimilarmolecules}) && $Options{percentsimilarmolecules} > 0 && $Options{percentsimilarmolecules} <= 100)) {
    die "Error: The value specified, $Options{percentsimilarmolecules}, for option \"-p, --PercentSimilarMolecules\" is not valid. Allowed values: > 0 and <= 100 \n";
  }
  if ($Options{quote} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not valid. Allowed values: Yes or No\n";
  }
  if (!IsPositiveInteger($Options{precision})) {
    die "Error: The value specified, $Options{precision}, for option \"--precision\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{referencecolmode} !~ /^(ColNum|ColLabel)$/i) {
    die "Error: The value specified, $Options{referencecolmode}, for option \"--ReferenceColMode\" is not valid. Allowed values: ColNum, or ColLabel\n";
  }
  if ($Options{referencecompoundidmode} !~ /^(DataField|MolName|LabelPrefix|MolNameOrLabelPrefix)$/i) {
    die "Error: The value specified, $Options{referencecompoundidmode}, for option \"--ReferenceCompoundIDMode\" is not valid. Allowed values: DataField, MolName, LabelPrefix or MolNameOrLabelPrefix\n";
  }
  if ($Options{searchmode} !~ /^(SimilaritySearch|DissimilaritySearch)$/i) {
    die "Error: The value specified, $Options{searchmode}, for option \"-s, --SearchMode\" is not valid. Allowed values: SimilaritySearch, DissimilaritySearch \n";
  }
  if ($Options{similarcountmode} !~ /^(NumOfSimilar|PercentSimilar)$/i) {
    die "Error: The value specified, $Options{similarcountmode}, for option \"--SimilarCountMode\" is not valid. Allowed values: NumOfSimilar, PercentSimilar \n";
  }
}

__END__

=head1 NAME

SimilaritySearchingFingerprints.pl - Perform similarity search using fingerprints strings data in SD, FP and CSV/TSV text file(s)

=head1 SYNOPSIS

SimilaritySearchingFingerprints.pl ReferenceFPFile DatabaseFPFile

SimilaritySearchingFingerprints.pl [B<--alpha> I<number>] [B<--beta> I<number>]
[B<-b, --BitVectorComparisonMode> I<TanimotoSimilarity | TverskySimilarity | ...>]
[B<--DatabaseColMode> I<ColNum | ColLabel>] [B<--DatabaseCompoundIDCol> I<col number | col name>]
[B<--DatabaseCompoundIDPrefix> I<text>] [B<--DatabaseCompoundIDField> I<DataFieldName>]
[B<--DatabaseCompoundIDMode> I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>]
[B<--DatabaseDataCols> I<"DataColNum1, DataColNum2,... " | DataColLabel1, DataCoLabel2,... ">]
[B<--DatabaseDataColsMode> I<All | Specify | CompoundID>] [B<--DatabaseDataFields> I<"FieldLabel1, FieldLabel2,... ">]
[B<--DatabaseDataFieldsMode> I<All | Common | Specify | CompoundID>]
[B<--DatabaseFingerprintsCol> I<col number | col name>] [B<--DatabaseFingerprintsField> I<FieldLabel>]
[]B<--DistanceCutoff> I<number>] [B<-d, --detail> I<InfoLevel>] [B<-f, --fast>]
[B<--FingerprintsMode> I<AutoDetect | FingerprintsBitVectorString | FingerprintsVectorString>]
[B<-g, --GroupFusionRule> I<Max, Mean, Median, Min, Sum, Euclidean>] [B<--GroupFusionApplyCutoff> I<Yes | No>]
[B<-h, --help>]  [B<--InDelim> I<comma | semicolon>] [B<-k, --KNN> I<all | number>]
[B<-m, --mode> I<IndividualReference | MultipleReferences>]
[B<-n, --NumOfSimilarMolecules> I<number>] [B<--OutDelim> I<comma | tab | semicolon>]
[B<--output> I<SD | text | both>] [B<-o, --overwrite>]
[B<-p, --PercentSimilarMolecules> I<number>] [B<--precision> I<number>] [B<-q, --quote> I<Yes | No>]
[B<--ReferenceColMode> I<ColNum | ColLabel>] [B<--ReferenceCompoundIDCol> I<col number | col name>]
[B<--ReferenceCompoundIDPrefix> I<text>] [B<--ReferenceCompoundIDField> I<DataFieldName>]
[B<--ReferenceCompoundIDMode> I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>]
[B<--ReferenceFingerprintsCol> I<col number | col name>] [B<--ReferenceFingerprintsField> I<FieldLabel>]
[B<-r, --root> I<RootName>] [B<-s, --SearchMode> I<SimilaritySearch | DissimilaritySearch>]
[B<--SimilarCountMode> I<NumOfSimilar | PercentSimilar>] [B<--SimilarityCutoff> I<number>]
[B<-v, --VectorComparisonMode> I<TanimotoSimilairy | ... | ManhattanDistance | ...>]
[B<--VectorComparisonFormulism> I<AlgebraicForm | BinaryForm | SetTheoreticForm>]
[B<-w, --WorkingDir> dirname] ReferenceFingerprintsFile DatabaseFingerprintsFile

=head1 DESCRIPTION

Perform molecular similarity search [ Ref 94-113 ] using fingerprint bit-vector or vector strings
data in I<SD, FP, or CSV/TSV text> files corresponding to I<ReferenceFingerprintsFile> and
I<DatabaseFingerprintsFile>, and generate SD and CSV/TSV text file(s) containing database
molecules which are similar to reference molecule(s). The reference molecules are also referred
to as query or seed molecules and database molecules as target molecules in the literature.

The current release of MayaChemTools supports two types of similarity search modes:
I<IndividualReference or MultipleReferences>. For default value of I<MultipleReferences> for B<-m, --mode>
option, reference molecules are considered as a set and B<-g, --GroupFusionRule> is used to calculate
similarity of a database molecule against reference molecules set. The group fusion rule is also
referred to as data fusion or consensus scoring in the literature. However, for I<IndividualReference>
value of B<-m, --mode> option, reference molecules are treated as individual molecules and each reference
molecule is compared against a database molecule by itself to identify similar molecules.

The molecular dissimilarity search can also be performed using I<DissimilaritySearch> value for
B<-s, --SearchMode> option. During dissimilarity search or usage of distance comparison coefficient
in similarity similarity search, the meaning of fingerprints comparison value is automatically reversed
as shown below:

    SeachMode      ComparisonCoefficient  ResultsSort   ComparisonValues

    Similarity     SimilarityCoefficient  Descending    Higher value imples
                                                        high similarity
    Similarity     DistanceCoefficient    Ascending     Lower value implies
                                                        high similarity

    Dissimilarity  SimilarityCoefficient  Ascending     Lower value implies
                                                        high dissimilarity
    Dissimilarity  DistanceCoefficient    Descending    Higher value implies
                                                        high dissimilarity

During I<IndividualReference> value of  B<-m, --Mode> option for similarity search, fingerprints bit-vector
or vector string of each reference molecule is compared with database molecules using specified
similarity or distance coefficients to identify most similar molecules for each reference molecule.
Based on value of B<--SimilarCountMode>, up to B<--n, --NumOfSimilarMolecules> or B<-p,
--PercentSimilarMolecules> at specified B<--SimilarityCutoff> or B<--DistanceCutoff> are
identified for each reference molecule.

During I<MultipleReferences> value B<-m, --mode> option for similarity search, all reference molecules
are considered as a set and B<-g, --GroupFusionRule> is used to calculate similarity of a database
molecule against reference molecules set either using all reference molecules or number of k-nearest
neighbors (k-NN) to a database molecule specified using B<-k, --kNN>. The fingerprints bit-vector
or vector string of each reference molecule in a set is compared with a database molecule using
a similarity or distance coefficient specified via B<-b, --BitVectorComparisonMode> or B<-v,
--VectorComparisonMode>. The reference molecules whose comparison values with a database
molecule fall outside specified B<--SimilarityCutoff> or B<--DistanceCutoff> are ignored during I<Yes>
value of B<--GroupFusionApplyCutoff>. The specified B<-g, --GroupFusionRule> is applied to
B<-k, --kNN> reference molecules to calculate final similarity value between a database molecule
and reference molecules set.

The input fingerprints I<SD, FP, or Text (CSV/TSV)> files for I<ReferenceFingerprintsFile> and
I<DatabaseTextFile> must contain valid fingerprint bit-vector or vector strings data corresponding to
same type of fingerprints.

The valid fingerprints I<SDFile> extensions are I<.sdf> and I<.sd>. The valid fingerprints I<FPFile>
extensions are I<.fpf> and I<.fp>. The valid fingerprints I<TextFile (CSV/TSV)> extensions are
I<.csv> and I<.tsv> for comma/semicolon and tab delimited text files respectively. The B<--indelim>
option determines the format of I<TextFile>. Any file which doesn't correspond to the format indicated
by B<--indelim> option is ignored.

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

Example of CSV I<TextFile> containing fingerprints bit-vector string data:

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

=item B<-b, --BitVectorComparisonMode> I<TanimotoSimilarity | TverskySimilarity | ...>

Specify what similarity coefficient to use for calculating similarity between fingerprints bit-vector
string data values in I<ReferenceFingerprintsFile> and I<DatabaseFingerprintsFile> during similarity
search. Possible values: I<TanimotoSimilarity | TverskySimilarity | ...>. Default: I<TanimotoSimilarity>

The current release supports the following similarity coefficients: I<BaroniUrbaniSimilarity, BuserSimilarity,
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

=item B<--DatabaseColMode> I<ColNum | ColLabel>

Specify how columns are identified in database fingerprints I<TextFile>: using column
number or column label. Possible values: I<ColNum or ColLabel>. Default value: I<ColNum>.

=item B<--DatabaseCompoundIDCol> I<col number | col name>

This value is B<--DatabaseColMode> mode specific. It specifies column to use for retrieving compound
ID from database fingerprints I<TextFile> during similarity and dissimilarity search for output SD and
CSV/TSV text files. Possible values: I<col number or col label>. Default value: I<first column containing
the word compoundID in its column label or sequentially generated IDs>.

This is only used for I<CompoundID> value of B<--DatabaseDataColsMode> option.

=item B<--DatabaseCompoundIDPrefix> I<text>

Specify compound ID prefix to use during sequential generation of compound IDs for database fingerprints
I<SDFile> and I<TextFile>. Default value: I<Cmpd>. The default value generates compound IDs which look
like Cmpd<Number>.

For database fingerprints I<SDFile>, this value is only used during I<LabelPrefix | MolNameOrLabelPrefix>
values of B<--DatabaseCompoundIDMode> option; otherwise, it's ignored.

Examples for I<LabelPrefix> or I<MolNameOrLabelPrefix> value of B<--DatabaseCompoundIDMode>:

    Compound

The values specified above generates compound IDs which correspond to Compound<Number>
instead of default value of Cmpd<Number>.

=item B<--DatabaseCompoundIDField> I<DataFieldName>

Specify database fingerprints I<SDFile> datafield label for generating compound IDs. This value is
only used during I<DataField> value of B<--DatabaseCompoundIDMode> option.

Examples for I<DataField> value of B<--DatabaseCompoundIDMode>:

    MolID
    ExtReg

=item B<--DatabaseCompoundIDMode> I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>

Specify how to generate compound IDs from database fingerprints I<SDFile> during similarity and
dissimilarity search for output SD and CSV/TSV text files: use a I<SDFile> datafield value; use
molname line from I<SDFile>; generate a sequential ID with specific prefix; use combination of both
MolName and LabelPrefix with usage of LabelPrefix values for empty molname lines.

Possible values: I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>.
Default: I<LabelPrefix>.

For I<MolNameAndLabelPrefix> value of B<--DatabaseCompoundIDMode>, molname line in I<SDFile> takes
precedence over sequential compound IDs generated using I<LabelPrefix> and only empty molname
values are replaced with sequential compound IDs.

This is only used for I<CompoundID> value of B<--DatabaseDataFieldsMode> option.

=item B<--DatabaseDataCols> I<"DataColNum1,DataColNum2,... " | DataColLabel1,DataCoLabel2,... ">

This value is B<--DatabaseColMode> mode specific. It is a comma delimited list of database fingerprints
I<TextFile> data column numbers or labels to extract and write to SD and CSV/TSV text files along with
other information for I<SD | text | both> values of B<--output> option.

This is only used for I<Specify> value of B<--DatabaseDataColsMode> option.

Examples:

    1,2,3
    CompoundName,MolWt

=item B<--DatabaseDataColsMode> I<All | Specify | CompoundID>

Specify how data columns from database fingerprints I<TextFile> are transferred to output SD and
CSV/TSV text files along with other information for I<SD | text | both> values of B<--output> option:
transfer all data columns; extract specified data columns; generate a compound ID database compound
prefix. Possible values: I<All | Specify | CompoundID>. Default value: I<CompoundID>.

=item B<--DatabaseDataFields> I<"FieldLabel1,FieldLabel2,... ">

Comma delimited list of database fingerprints I<SDFile> data fields to extract and write to SD
and CSV/TSV text files along with other information for I<SD | text | both> values of
B<--output> option.

This is only used for I<Specify> value of B<--DatabaseDataFieldsMode> option.

Examples:

    Extreg
    MolID,CompoundName

=item B<--DatabaseDataFieldsMode> I<All | Common | Specify | CompoundID>

Specify how data fields from database fingerprints I<SDFile> are transferred to output SD and
CSV/TSV text files along with other information for I<SD | text | both> values of B<--output>
option: transfer all SD data field; transfer SD data files common to all compounds; extract
specified data fields; generate a compound ID using molname line, a compound prefix, or a
combination of both. Possible values: I<All | Common | specify | CompoundID>. Default value:
I<CompoundID>.

=item B<--DatabaseFingerprintsCol> I<col number | col name>

This value is B<--DatabaseColMode> specific. It specifies fingerprints column to use during similarity
and dissimilarity search for database fingerprints I<TextFile>. Possible values: I<col number or col label>.
Default value: I<first column containing the word Fingerprints in its column label>.

=item B<--DatabaseFingerprintsField> I<FieldLabel>

Fingerprints field label to use during similarity and dissimilarity search for database fingerprints I<SDFile>.
Default value: I<first data field label containing the word Fingerprints in its label>

=item B<--DistanceCutoff> I<number>

Distance cutoff value to use during comparison of distance value between a pair of database
and reference molecule calculated by distance comparison methods for fingerprints vector
string data values. Possible values: I<Any valid number>. Default value: I<10>.

The comparison value between a pair of database and reference molecule must meet the cutoff
criterion as shown below:

    SeachMode      CutoffCriterion  ComparisonValues

    Similarity     <=               Lower value implies high similarity
    Dissimilarity  >=               Higher value implies high dissimilarity

This option is only used during distance coefficients values of B<-v, --VectorComparisonMode>
option.

This option is ignored during I<No> value of B<--GroupFusionApplyCutoff> for I<MultipleReferences>
B<-m, --mode>.

=item B<-d, --detail> I<InfoLevel>

Level of information to print about lines being ignored. Default: I<1>. Possible values:
I<1, 2 or 3>.

=item B<-f, --fast>

In this mode, fingerprints columns specified using B<--FingerprintsCol> for reference and database
fingerprints I<TextFile(s)>, and B<--FingerprintsField> for reference and database fingerprints I<SDFile(s)>
are assumed to contain valid fingerprints data and no checking is performed before performing similarity
and dissimilarity search. By default, fingerprints data is validated before computing pairwise similarity and
distance coefficients.

=item B<--FingerprintsMode> I<AutoDetect | FingerprintsBitVectorString | FingerprintsVectorString>

Format of fingerprint strings data in reference and database fingerprints I<SD, FP, or Text (CSV/TSV)>
files: automatically detect format of fingerprints string created by MayaChemTools fingerprints
generation scripts or explicitly specify its format. Possible values: I<AutoDetect | FingerprintsBitVectorString |
FingerprintsVectorString>. Default value: I<AutoDetect>.

=item B<-g, --GroupFusionRule> I<Max, Min, Mean, Median, Sum, Euclidean>

Specify what group fusion [ Ref 94-97, Ref 100, Ref 105 ] rule to use for calculating similarity of
a database molecule against a set of reference molecules during I<MultipleReferences> value of
similarity search B<-m, --mode>. Possible values: I<Max, Min, Mean, Median, Sum, Euclidean>. Default
value: I<Max>. I<Mean> value corresponds to average or arithmetic mean. The group fusion rule is
also referred to as data fusion of consensus scoring in the literature.

For a reference molecules set and a database molecule, let:

    N = Number of reference molecules in a set

    i = ith reference reference molecule in a set
    n = Nth reference reference molecule in a set

    d = dth database molecule

    Crd = Fingerprints comparison value between rth reference and dth database
          molecule - similarity/dissimilarity comparison using similarity or
          distance coefficient

Then, various group fusion rules to calculate fused similarity between a database molecule and
reference molecules set are defined as follows:

B<Max>: MAX ( C1d, C2d, ..., Cid, ..., Cnd )

B<Min>: MIN ( C1d, C2d, ..., Cid, ..., Cnd )

B<Mean>: SUM ( C1d, C2d, ..., Cid, ..., Cnd ) / N

B<Median>: MEDIAN (  C1d, C2d, ..., Cid, ..., Cnd )

B<Sum>: SUM (  C1d, C2d, ..., Cid, ..., Cnd )

B<Euclidean>: SQRT( SUM( C1d ** 2, C2d ** 2, ..., Cid ** 2, ..., Cnd *** 2) )

The fingerprints bit-vector or vector string of each reference molecule in a set is compared
with a database molecule using a similarity or distance coefficient specified via B<-b,
--BitVectorComparisonMode> or B<-v, --VectorComparisonMode>. The reference molecules
whose comparison values with a database molecule fall outside specified B<--SimilarityCutoff>
or B<--DistanceCutoff> are ignored during I<Yes> value of B<--GroupFusionApplyCutoff>. The
specified B<-g, --GroupFusionRule> is applied to B<-k, --kNN> reference molecules to calculate
final fused similarity value between a database molecule and reference molecules set.

During dissimilarity search or usage of distance comparison coefficient in similarity search,
the meaning of fingerprints comaprison value is automatically reversed as shown below:

    SeachMode      ComparisonCoefficient  ComparisonValues

    Similarity     SimilarityCoefficient  Higher value imples high similarity
    Similarity     DistanceCoefficient    Lower value implies high similarity

    Dissimilarity  SimilarityCoefficient  Lower value implies high
                                          dissimilarity
    Dissimilarity  DistanceCoefficient    Higher value implies high
                                          dissimilarity

Consequently, I<Max> implies highest and lowest comparison value for usage of similarity and
distance coefficient respectively during similarity search. And it corresponds to lowest and highest
comparison value for usage of similarity and distance coefficient respectively during dissimilarity
search. During I<Min> fusion rule, the highest and lowest comparison values are appropriately
reversed.

=item B<--GroupFusionApplyCutoff> I<Yes | No>

Specify whether to apply B<--SimilarityCutoff> or B<--DistanceCutoff> values during application
of B<-g, --GroupFusionRule> to reference molecules set. Possible values: I<Yes or No>. Default
value: I<Yes>.

During I<Yes> value of B<--GroupFusionApplyCutoff>, the reference molecules whose comparison
values with a database molecule fall outside specified B<--SimilarityCutoff> or B<--DistanceCutoff>
are not used to calculate final fused similarity value between a database molecule and reference
molecules set.

=item B<-h, --help>

Print this help message.

=item B<--InDelim> I<comma | semicolon>

Input delimiter for reference and database fingerprints CSV I<TextFile(s)>. Possible values:
I<comma or semicolon>. Default value: I<comma>. For TSV files, this option is ignored
and I<tab> is used as a delimiter.

=item B<-k, --kNN> I<all | number>

Number of k-nearest neighbors (k-NN) reference molecules to use during B<-g, --GroupFusionRule>
for calculating similarity of a database molecule against a set of reference molecules. Possible values:
I<all | positive integers>. Default: I<all>.

After ranking similarity values between a database molecule and reference molecules during
I<MultipleReferences> value of similarity search B<-m, --mode> option, a top B<-k, --KNN> reference
molecule are selected and used during B<-g, --GroupFusionRule>.

This option is B<-s, --SearchMode> dependent: It corresponds to dissimilar molecules during
I<DissimilaritySearch> value of B<-s, --SearchMode> option.

=item B<-m, --mode> I<IndividualReference | MultipleReferences>

Specify how to treat reference molecules in I<ReferenceFingerprintsFile> during similarity search:
Treat each reference molecule individually during similarity search or perform similarity
search by treating multiple reference molecules as a set. Possible values: I<IndividualReference
| MultipleReferences>. Default value: I<MultipleReferences>.

During I<IndividualReference> value of  B<-m, --Mode> for similarity search, fingerprints bit-vector
or vector string of each reference molecule is compared with database molecules using specified
similarity or distance coefficients to identify most similar molecules for each reference molecule.
Based on value of B<--SimilarCountMode>, upto B<--n, NumOfSimilarMolecules> or B<-p,
--PercentSimilarMolecules> at specified <--SimilarityCutoff> or B<--DistanceCutoff> are
identified for each reference molecule.

During I<MultipleReferences> value B<-m, --mode> for similarity search, all reference molecules
are considered as a set and B<-g, --GroupFusionRule> is used to calculate similarity of a database
molecule against reference molecules set either using all reference molecules or number of k-nearest
neighbors (k-NN) to a database molecule specified using B<-k, --kNN>. The fingerprints bit-vector
or vector string of each reference molecule in a set is compared with a database molecule using
a similarity or distance coefficient specified via B<-b, --BitVectorComparisonMode> or B<-v,
--VectorComparisonMode>. The reference molecules whose comparison values with a database
molecule fall outside specified B<--SimilarityCutoff> or B<--DistanceCutoff> are ignored. The
specified B<-g, --GroupFusionRule> is applied to rest of B<-k, --kNN> reference molecules to calculate
final similarity value between a database molecule and reference molecules set.

The meaning of similarity and distance is automatically reversed during I<DissimilaritySearch> value
of B<-s, --SearchMode> along with appropriate handling of B<--SimilarityCutoff> or
B<--DistanceCutoff> values.

=item B<-n, --NumOfSimilarMolecules> I<number>

Maximum number of most similar database molecules to find for each reference molecule or set of
reference molecules based on I<IndividualReference> or I<MultipleReferences> value of similarity
search B<-m, --mode> option. Default: I<10>. Valid values: positive integers.

This option is ignored during I<PercentSimilar> value of B<--SimilarCountMode> option.

This option is B<-s, --SearchMode> dependent: It corresponds to dissimilar molecules during
I<DissimilaritySearch> value of B<-s, --SearchMode> option.

=item B<--OutDelim> I<comma | tab | semicolon>

Delimiter for output CSV/TSV text file. Possible values: I<comma, tab, or semicolon>
Default value: I<comma>.

=item B<--output> I<SD | text | both>

Type of output files to generate. Possible values: I<SD, text, or both>. Default value: I<text>.

=item B<-o, --overwrite>

Overwrite existing files

=item B<-p, --PercentSimilarMolecules> I<number>

Maximum percent of mosy similar database molecules to find for each reference molecule or set of
reference molecules based on I<IndividualReference> or I<MultipleReferences> value of similarity
search B<-m, --mode> option. Default: I<1> percent of database molecules. Valid values: non-zero values
in between I<0 to 100>.

This option is ignored during I<NumOfSimilar> value of B<--SimilarCountMode> option.

During I<PercentSimilar> value of B<--SimilarCountMode> option, the number of molecules
in I<DatabaseFingerprintsFile> is counted and number of similar molecules correspond to
B<--PercentSimilarMolecules> of the total number of database molecules.

This option is B<-s, --SearchMode> dependent: It corresponds to dissimilar molecules during
I<DissimilaritySearch> value of B<-s, --SearchMode> option.

=item B<--precision> I<number>

Precision of calculated similarity values for comparison and generating output files. Default: up to I<2>
decimal places. Valid values: positive integers.

=item B<-q, --quote> I<Yes | No>

Put quote around column values in output CSV/TSV text file. Possible values:
I<Yes or No>. Default value: I<Yes>.

=item B<--ReferenceColMode> I<ColNum | ColLabel>

Specify how columns are identified in reference fingerprints I<TextFile>: using column
number or column label. Possible values: I<ColNum or ColLabel>. Default value: I<ColNum>.

=item B<--ReferenceCompoundIDCol> I<col number | col name>

This value is B<--ReferenceColMode> mode specific. It specifies column to use for retrieving compound
ID from reference fingerprints I<TextFile> during similarity and dissimilarity search for output SD and CSV/TSV
text files. Possible values: I<col number or col label>. Default value: I<first column containing the word compoundID
in its column label or sequentially generated IDs>.

=item B<--ReferenceCompoundIDPrefix> I<text>

Specify compound ID prefix to use during sequential generation of compound IDs for reference fingerprints
I<SDFile> and I<TextFile>. Default value: I<Cmpd>. The default value generates compound IDs which looks
like Cmpd<Number>.

For reference fingerprints I<SDFile>, this value is only used during I<LabelPrefix | MolNameOrLabelPrefix>
values of B<--ReferenceCompoundIDMode> option; otherwise, it's ignored.

Examples for I<LabelPrefix> or I<MolNameOrLabelPrefix> value of B<--DatabaseCompoundIDMode>:

    Compound

The values specified above generates compound IDs which correspond to Compound<Number>
instead of default value of Cmpd<Number>.

=item B<--ReferenceCompoundIDField> I<DataFieldName>

Specify reference fingerprints I<SDFile> datafield label for generating compound IDs.
This value is only used during I<DataField> value of B<--ReferenceCompoundIDMode> option.

Examples for I<DataField> value of B<--ReferenceCompoundIDMode>:

    MolID
    ExtReg

=item B<--ReferenceCompoundIDMode> I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>

Specify how to generate compound IDs from reference fingerprints I<SDFile> during similarity and
dissimilarity search for output SD and CSV/TSV text files: use a I<SDFile> datafield value; use
molname line from I<SDFile>; generate a sequential ID with specific prefix; use combination of both
MolName and LabelPrefix with usage of LabelPrefix values for empty molname lines.

Possible values: I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>.
Default: I<LabelPrefix>.

For I<MolNameAndLabelPrefix> value of B<--ReferenceCompoundIDMode>, molname line in I<SDFiles>
takes precedence over sequential compound IDs generated using I<LabelPrefix> and only empty molname
values are replaced with sequential compound IDs.

=item B<--ReferenceFingerprintsCol> I<col number | col name>

This value is B<--ReferenceColMode> specific. It specifies fingerprints column to use during similarity
and dissimilarity search for reference fingerprints I<TextFile>. Possible values: I<col number or col label>.
Default value: I<first column containing the word Fingerprints in its column label>.

=item B<--ReferenceFingerprintsField> I<FieldLabel>

Fingerprints field label to use during similarity and dissimilarity search for reference fingerprints I<SDFile>.
Default value: I<first data field label containing the word Fingerprints in its label>

=item B<-r, --root> I<RootName>

New file name is generated using the root: <Root>.<Ext>. Default for new file name:
<ReferenceFileName>SimilaritySearching.<Ext>. The output file type determines <Ext>
value. The sdf, csv, and tsv <Ext> values are used for SD, comma/semicolon, and tab delimited
text files respectively.

=item B<-s, --SearchMode> I<SimilaritySearch | DissimilaritySearch>

Specify how to find molecules from database molecules for individual reference molecules or
set of reference molecules: Find similar molecules or dissimilar molecules from database molecules.
Possible values: I<SimilaritySearch | DissimilaritySearch>. Default value: I<SimilaritySearch>.

During I<DissimilaritySearch> value of B<-s, --SearchMode> option, the meaning of the following
options is switched and they correspond to dissimilar molecules instead of similar molecules:
B<--SimilarCountMode>, B<-n, --NumOfSimilarMolecules>, B<--PercentSimilarMolecules>,
B<-k, --kNN>.

=item B<--SimilarCountMode> I<NumOfSimilar | PercentSimilar>

Specify method used to count similar molecules found from database molecules for individual
reference molecules or set of reference molecules: Find number of similar molecules or percent
of similar molecules from database molecules. Possible values: I<NumOfSimilar | PercentSimilar>.
Default value: I<NumOfSimilar>.

The values for number of similar molecules and percent similar molecules are specified
using options B<-n, NumOfSimilarMolecule> and B<--PercentSimilarMolecules>.

This option is B<-s, --SearchMode> dependent: It corresponds to dissimilar molecules during
I<DissimilaritySearch> value of B<-s, --SearchMode> option.

=item B<--SimilarityCutoff> I<number>

Similarity cutoff value to use during comparison of similarity value between a pair of database
and reference molecules calculated by similarity comparison methods for fingerprints bit-vector
vector strings data values. Possible values: I<Any valid number>. Default value: I<0.75>.

The comparison value between a pair of database and reference molecule must meet the cutoff
criterion as shown below:

    SeachMode      CutoffCriterion  ComparisonValues

    Similarity     >=               Higher value implies high similarity
    Dissimilarity  <=               Lower value implies high dissimilarity

This option is ignored during I<No> value of B<--GroupFusionApplyCutoff> for I<MultipleReferences>
B<-m, --mode>.

This option is B<-s, --SearchMode> dependent: It corresponds to dissimilar molecules during
I<DissimilaritySearch> value of B<-s, --SearchMode> option.

=item B<-v, --VectorComparisonMode> I<SupportedSimilarityName | SupportedDistanceName>

Specify what similarity or distance coefficient to use for calculating similarity between fingerprint
vector strings data values in I<ReferenceFingerprintsFile> and I<DatabaseFingerprintsFile> during
similarity search. Possible values:  I<TanimotoSimilairy | ... | ManhattanDistance | ...>. Default
value: I<TanimotoSimilarity>.

The value of B<-v, --VectorComparisonMode>, in conjunction with B<--VectorComparisonFormulism>,
decides which type of similarity and distance coefficient formulism gets used.

The current releases supports the following similarity and distance coefficients: I<CosineSimilarity,
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

=item B<--VectorComparisonFormulism> I<AlgebraicForm | BinaryForm | SetTheoreticForm>

Specify fingerprints vector comparison formulism to use for calculation similarity and distance
coefficients during B<-v, --VectorComparisonMode>. Possible values: I<AlgebraicForm | BinaryForm |
SetTheoreticForm>. Default value: I<AlgebraicForm>.

For fingerprint vector strings containing B<AlphaNumericalValues> data values - B<ExtendedConnectivityFingerprints>,
B<AtomNeighborhoodsFingerprints> and so on - all three formulism result in same value during similarity and distance
calculations.

=item B<-w, --WorkingDir> I<DirName>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To perform similarity search using Tanimoto coefficient by treating all reference molecules as a set
to find 10 most similar database molecules with application of Max group fusion rule and similarity
cutoff to supported fingerprints strings data in SD fingerprints files present in a data fields with
Fingerprint substring in their labels, and create a ReferenceFPHexSimilaritySearching.csv file containing
sequentially generated database compound IDs with Cmpd prefix, type:

    % SimilaritySearchingFingerprints.pl -o ReferenceSampleFPHex.sdf
      DatabaseSampleFPHex.sdf

To perform similarity search using Tanimoto coefficient by treating all reference molecules as a set
to find 10 most similar database molecules with application of Max group fusion rule and similarity
cutoff to supported fingerprints strings data in FP fingerprints files, and create a
SimilaritySearchResults.csv file containing database compound IDs retireved from FP file, type:

    % SimilaritySearchingFingerprints.pl -r SimilaritySearchResults -o
      ReferenceSampleFPBin.fpf DatabaseSampleFPBin.fpf

To perform similarity search using Tanimoto coefficient by treating all reference molecules as a set
to find 10 most similar database database molecules with application of Max group fusion rule and
similarity cutoff to supported fingerprints strings data in text fingerprints files present in a column
names containing Fingerprint substring in their names, and create a ReferenceFPHexSimilaritySearching.csv
file containing database compound IDs retireved column name containing CompoundID substring or
sequentially generated compound IDs, type:

    % SimilaritySearchingFingerprints.pl -o ReferenceSampleFPCount.csv
      DatabaseSampleFPCount.csv

To perform similarity search  using Tanimoto coefficient by treating reference molecules as individual molecules
to find 10 most similar database molecules for each reference molecule with application of similarity cutoff to
supported fingerprints strings data in SD fingerprints files present in a data fields with Fingerprint substring
in their labels, and create a ReferenceFPHexSimilaritySearching.csv file containing sequentially generated
reference and database compound IDs with Cmpd prefix, type:

    % SimilaritySearchingFingerprints.pl -mode IndividualReference -o
      ReferenceSampleFPHex.sdf DatabaseSampleFPHex.sdf

To perform similarity search  using Tanimoto coefficient by treating reference molecules as individual molecules
to find 10 most similar database molecules for each reference molecule with application of similarity cutoff to
supported fingerprints strings data in FP fingerprints files, and create a ReferenceFPHexSimilaritySearching.csv
file containing references and database compound IDs retireved from FP file, type:

    % SimilaritySearchingFingerprints.pl -mode IndividualReference -o
      ReferenceSampleFPHex.fpf DatabaseSampleFPHex.fpf

To perform similarity search  using Tanimoto coefficient by treating reference molecules as individual molecules
to find 10 most similar database molecules for each reference molecule with application of similarity cutoff to
supported fingerprints strings data in text fingerprints files present in a column names containing Fingerprint
substring in their names, and create a ReferenceFPHexSimilaritySearching.csv file containing reference and
database compound IDs retrieved column name containing CompoundID substring or sequentially generated
compound IDs, type:

    % SimilaritySearchingFingerprints.pl -mode IndividualReference -o
      ReferenceSampleFPHex.csv DatabaseSampleFPHex.csv

To perform dissimilarity search using Tanimoto coefficient by treating all reference molecules as a set
to find 10 most dissimilar database molecules with application of Max group fusion rule and similarity
cutoff to supported fingerprints strings data in SD fingerprints files present in a data fields with
Fingerprint substring in their labels, and create a ReferenceFPHexSimilaritySearching.csv file containing
sequentially generated database compound IDs with Cmpd prefix, type:

    % SimilaritySearchingFingerprints.pl --mode MultipleReferences --SearchMode
      DissimilaritySearch -o ReferenceSampleFPHex.sdf DatabaseSampleFPHex.sdf

To perform similarity search using CityBlock distance by treating reference molecules as individual molecules
to find 10 most similar database molecules for each reference molecule with application of distance cutoff
to supported vector fingerprints strings data in SD fingerprints files present in a data fields with Fingerprint
substring in their labels, and create a ReferenceFPHexSimilaritySearching.csv file containing sequentially generated
reference and database compound IDs with Cmpd prefix, type:

    % SimilaritySearchingFingerprints.pl -mode IndividualReference
      --VectorComparisonMode CityBlockDistance --VectorComparisonFormulism
      AlgebraicForm --DistanceCutoff 10 -o
      ReferenceSampleFPCount.sdf DatabaseSampleFPCount.sdf

To perform similarity search using Tanimoto coefficient by treating all reference molecules as a set
to find 100 most similar database molecules with application of Mean group fusion rule to to top 10
reference molecules with in similarity cutoff of 0.75 to supported fingerprints strings data in FP fingerprints
files, and create a ReferenceFPHexSimilaritySearching.csv file containing database compound IDs retrieved
from FP file, type:

    % SimilaritySearchingFingerprints.pl --mode MultipleReferences --SearchMode
      SimilaritySearch --BitVectorComparisonMode TanimotoSimilarity
      --GroupFusionRule Mean --GroupFusionApplyCutoff Yes --kNN 10
      --SimilarityCutoff 0.75 --SimilarCountMode NumOfSimilar
      --NumOfSimilarMolecules 100 -o
      ReferenceSampleFPHex.fpf DatabaseSampleFPHex.fpf

To perform similarity search  using Tanimoto coefficient by treating reference molecules as individual molecules
to find 2 percent of most similar database molecules for each reference molecule with application of similarity
cutoff of 0.85 to supported fingerprints strings data in text fingerprints files present in specific columns and
create a ReferenceFPHexSimilaritySearching.csv file containing reference and database compoundIDs retrieved
from specific columns, type:

    % SimilaritySearchingFingerprints.pl --mode IndividualReference --SearchMode
      SimilaritySearch --BitVectorComparisonMode TanimotoSimilarity
      --ReferenceColMode ColLabel --ReferenceFingerprintsCol Fingerprints
      --ReferenceCompoundIDCol CompoundID --DatabaseColMode Collabel
      --DatabaseCompoundIDCol CompoundID --DatabaseFingerprintsCol
      Fingerprints --SimilarityCutoff 0.85 --SimilarCountMode PercentSimilar
      --PercentSimilarMolecules 2 -o
      ReferenceSampleFPHex.csv DatabaseSampleFPHex.csv

To perform similarity search  using Tanimoto coefficient by treating reference molecules as individual molecules
to find top 50 most similar database molecules for each reference molecule with application of similarity
cutoff of 0.85 to supported fingerprints strings data in SD fingerprints files present in specific data fields and
create both ReferenceFPHexSimilaritySearching.csv and ReferenceFPHexSimilaritySearching.sdf files containing
reference and database compoundIDs retrieved from specific data fields, type:

    % SimilaritySearchingFingerprints.pl --mode IndividualReference --SearchMode
      SimilaritySearch --BitVectorComparisonMode TanimotoSimilarity
      --ReferenceFingerprintsField Fingerprints
      --DatabaseFingerprintsField Fingerprints
      --ReferenceCompoundIDMode DataField --ReferenceCompoundIDField CmpdID
      --DatabaseCompoundIDMode DataField --DatabaseCompoundIDField CmpdID
      --SimilarityCutoff 0.85 --SimilarCountMode NumOfSimilar
      --NumOfSimilarMolecules 50 --output both -o
      ReferenceSampleFPHex.sdf DatabaseSampleFPHex.sdf

To perform similarity search  using Tanimoto coefficient by treating reference molecules as individual molecules
to find 1 percent of  most similar database molecules for each reference molecule with application of similarity
cutoff to supported fingerprints strings data in SD fingerprints files present in specific data field labels, and create
both ReferenceFPHexSimilaritySearching.csv ReferenceFPHexSimilaritySearching.sdf files containing reference and
database compound IDs retrieved from specific data field labels along with other specific data for database
molecules, type:

    % SimilaritySearchingFingerprints.pl --mode IndividualReference --SearchMode
      SimilaritySearch --BitVectorComparisonMode TanimotoSimilarity
      --ReferenceFingerprintsField Fingerprints
      --DatabaseFingerprintsField Fingerprints
      --ReferenceCompoundIDMode DataField --ReferenceCompoundIDField CmpdID
      --DatabaseCompoundIDMode DataField --DatabaseCompoundIDField CmpdID
      --DatabaseDataFieldsMode Specify --DatabaseDataFields "TPSA,SLogP"
      --SimilarityCutoff 0.75 --SimilarCountMode PercentSimilar
      --PercentSimilarMolecules 1 --output both --OutDelim comma --quote Yes
      --precision 3 -o ReferenceSampleFPHex.sdf DatabaseSampleFPHex.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoFingerprintsFiles.pl, SimilarityMatricesFingerprints.pl, AtomNeighborhoodsFingerprints.pl,
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
