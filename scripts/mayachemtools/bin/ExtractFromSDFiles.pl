#!/usr/bin/perl -w
#
# File: ExtractFromSDFiles.pl
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
use SDFileUtil;
use FileUtil;
use TextUtil;

my($ScriptName, %Options, $StartTime, $EndTime, $TotalTime);

# Autoflush STDOUT
$| = 1;

# Starting message...
$ScriptName = basename($0);
print "\n$ScriptName:Starting...\n\n";
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

# Collect information about SD files...
print "Checking input SD file(s)...\n";
my(%SDFilesInfo);
RetrieveSDFilesInfo();

# Generate output files...
my($FileIndex);
if (@SDFilesList > 1) {
  print "\nProcessing SD files...\n";
}
for $FileIndex (0 .. $#SDFilesList) {
  if ($SDFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $SDFilesList[$FileIndex]...\n";
    ExtractFromSDFile($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Extract data from a SD file...
sub ExtractFromSDFile {
  my($FileIndex) = @_;

  OpenInputAndOutputFiles($FileIndex);

  MODE: {
    if ($OptionsInfo{Mode} =~ /^AllDataFields$/i) {
      ExtractAllDataFields($FileIndex);
      last MODE;
    }
    if ($OptionsInfo{Mode} =~ /^CommonDataFields$/i) {
      ExtractCommonDataFields($FileIndex);
      last MODE;
    }
    if ($OptionsInfo{Mode} =~ /^DataFields$/i) {
      ExtractDataFields($FileIndex);
      last MODE;
    }
    if ($OptionsInfo{Mode} =~ /^(DataFieldByList|DatafieldUniqueByList)$/i) {
      ExtractDataFieldByList($FileIndex);
      last MODE;
    }
    if ($OptionsInfo{Mode} =~ /^DataFieldNotByList$/i) {
      ExtractDataFieldNotByList($FileIndex);
      last MODE;
    }
    if ($OptionsInfo{Mode} =~ /^DataFieldsByValue$/i) {
      ExtractDataFieldsByValue($FileIndex);
      last MODE;
    }
    if ($OptionsInfo{Mode} =~ /^DataFieldsByRegex$/i) {
      ExtractDataFieldsByRegex($FileIndex);
      last MODE;
    }
    if ($OptionsInfo{Mode} =~ /^RandomCmpds$/i) {
      ExtractRandomCompounds($FileIndex);
      last MODE;
    }
    if ($OptionsInfo{Mode} =~ /^MolNames$/i) {
      ExtractMolNames($FileIndex);
      last MODE;
    }
    if ($OptionsInfo{Mode} =~ /^RecordNum$/i) {
      ExtractRecordNum($FileIndex);
      last MODE;
    }
    if ($OptionsInfo{Mode} =~ /^RecordNums$/i) {
      ExtractRecordNums($FileIndex);
      last MODE;
    }
    if ($OptionsInfo{Mode} =~ /^RecordRange$/i) {
      ExtractRecordRange($FileIndex);
      last MODE;
    }
    if ($OptionsInfo{Mode} =~ /^2DCmpdRecords$/i) {
      Extract2DCmpdRecords($FileIndex);
      last MODE;
    }
    if ($OptionsInfo{Mode} =~ /^3DCmpdRecords$/i) {
      Extract3DCmpdRecords($FileIndex);
      last MODE;
    }
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: alldatafields, commondatafields, datafields, datafieldsbyvalue, datafieldbylist, datafielduniquebylist, datafieldnotbylist, molnames, randomcmpds, recordnum, recordnums, recordrange, 2dcmpdrecords, 3dcmpdrecords\n";
  }

  CloseInputAndOutputFiles();
}

# Extract all data fields...
sub ExtractAllDataFields {
  my($FileIndex) = @_;
  my(@CmpdLines);

  @{$SDFilesInfo{DataLabels}} = @{$SDFilesInfo{AllDataFieldLabels}[$FileIndex]};
  WriteTextFileColLabels();

  while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    @CmpdLines = split "\n", $SDFilesInfo{CmpdString};
    %{$SDFilesInfo{DataFieldValues}} = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

    SetupDataValues();
    WriteTextFileCmpdData();
    WriteSDFileCmpdData();
  }
}

# Extract common data fields...
sub ExtractCommonDataFields {
  my($FileIndex) = @_;
  my(@CmpdLines);

  @{$SDFilesInfo{DataLabels}} = @{$SDFilesInfo{CommonDataFieldLabels}[$FileIndex]};
  WriteTextFileColLabels();

  while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    @CmpdLines = split "\n", $SDFilesInfo{CmpdString};
    %{$SDFilesInfo{DataFieldValues}} = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

    SetupDataValues();
    WriteTextFileCmpdData();
    WriteSDFileCmpdData();
  }
}

# Extract specified data fields...
sub ExtractDataFields {
  my($FileIndex) = @_;
  my(@CmpdLines);

  @{$SDFilesInfo{DataLabels}} = @{$OptionsInfo{SpecifiedDataFieldLabels}};
  WriteTextFileColLabels();

  while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    @CmpdLines = split "\n", $SDFilesInfo{CmpdString};
    %{$SDFilesInfo{DataFieldValues}} = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

    SetupDataValues();
    WriteTextFileCmpdData();
    WriteSDFileCmpdData();
  }
}

# Extract data fields using a list...
sub ExtractDataFieldByList {
  my($FileIndex) = @_;
  my($CmpdNum, $Value, $SpecifiedDataFieldValuesFoundCount, $CurrentValue, $SpecifiedDataFieldLabel, @CmpdLines);

  @{$SDFilesInfo{DataLabels}} = @{$SDFilesInfo{AllDataFieldLabels}[$FileIndex]};
  WriteTextFileColLabels();

  for $Value (keys %{$OptionsInfo{SpecifiedDataFieldValues}}) {
    $OptionsInfo{SpecifiedDataFieldValues}{$Value} = "NotFound";
  }
  $SpecifiedDataFieldValuesFoundCount = 0;
  $SpecifiedDataFieldLabel = $OptionsInfo{SpecifiedDataFieldLabel};

  CMPDSTRING: while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    $CmpdNum++;

    @CmpdLines = split "\n", $SDFilesInfo{CmpdString};
    %{$SDFilesInfo{DataFieldValues}} = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

    if (!exists $SDFilesInfo{DataFieldValues}{$SpecifiedDataFieldLabel}) {
      next CMPDSTRING;
    }

    SetupDataValues();

    $SpecifiedDataFieldLabel = $OptionsInfo{SpecifiedDataFieldLabel};
    $CurrentValue = $SDFilesInfo{DataFieldValues}{$SpecifiedDataFieldLabel};

    if (exists $OptionsInfo{SpecifiedDataFieldValues}{$CurrentValue}) {
      if ($SpecifiedDataFieldValuesFoundCount < $OptionsInfo{SpecifiedDataFieldValuesCount}) {
	if ($OptionsInfo{SpecifiedDataFieldValues}{$CurrentValue} eq "NotFound") {
	  $SpecifiedDataFieldValuesFoundCount++;
	  $OptionsInfo{SpecifiedDataFieldValues}{$CurrentValue} = "Found";
	  if ($OptionsInfo{Mode} =~ /^DataFieldUniqueByList$/i) {
	    WriteSDFileCmpdString();
	    WriteTextFileCmpdData();
	  }
	}
	if ($OptionsInfo{Mode} =~ /^DataFieldByList$/i) {
	  WriteSDFileCmpdString();
	  WriteTextFileCmpdData();
	}
      }
      if ($SpecifiedDataFieldValuesFoundCount >= $OptionsInfo{SpecifiedDataFieldValuesCount}) {
	last CMPDSTRING;
      }
    }
  }
}

# Extract data field whose values are not on the specified list...
sub ExtractDataFieldNotByList {
  my($FileIndex) = @_;
  my($CurrentValue, $SpecifiedDataFieldLabel, @CmpdLines);

  @{$SDFilesInfo{DataLabels}} = @{$SDFilesInfo{AllDataFieldLabels}[$FileIndex]};
  WriteTextFileColLabels();

  $SpecifiedDataFieldLabel = $OptionsInfo{SpecifiedDataFieldLabel};

  CMPDSTRING: while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    @CmpdLines = split "\n", $SDFilesInfo{CmpdString};
    %{$SDFilesInfo{DataFieldValues}} = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

    if (!exists $SDFilesInfo{DataFieldValues}{$SpecifiedDataFieldLabel}) {
      next CMPDSTRING;
    }

    SetupDataValues();

    $CurrentValue = $SDFilesInfo{DataFieldValues}{$SpecifiedDataFieldLabel};

    # Make sure the current value is not empty and is not only specified list of values...
    if (IsEmpty($CurrentValue) || exists $OptionsInfo{SpecifiedDataFieldValues}{$CurrentValue}) {
      next CMPDSTRING;
    }

    WriteSDFileCmpdString();
    WriteTextFileCmpdData();
  }
}

# Extract data fields by value...
sub ExtractDataFieldsByValue {
  my($FileIndex) = @_;
  my($Label, $CurrentValue, $SpecifiedCriterion, $SpecifiedValue, $ViolationCount, $Nothing, @CmpdLines);

  @{$SDFilesInfo{DataLabels}} = @{$SDFilesInfo{AllDataFieldLabels}[$FileIndex]};
  WriteTextFileColLabels();

  CMPDSTRING: while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    @CmpdLines = split "\n", $SDFilesInfo{CmpdString};
    %{$SDFilesInfo{DataFieldValues}} = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

    SetupDataValues();
    $ViolationCount = 0;

    for $Label (@{$OptionsInfo{SpecifiedDataFieldLabels}}) {
      if (exists $SDFilesInfo{DataFieldValues}{$Label}) {
	$CurrentValue = $SDFilesInfo{DataFieldValues}{$Label};
	$SpecifiedCriterion = $OptionsInfo{SpecifiedDataFieldCriteriaMap}{$Label};
	$SpecifiedValue = $OptionsInfo{SpecifiedDataFieldValuesMap}{$Label};

	if ($OptionsInfo{NumericalComparison}) {
	  CRITERION: {
	      if ($SpecifiedCriterion =~ /^eq$/i) { if ($CurrentValue != $SpecifiedValue) { $ViolationCount++; last CRITERION; } }
	      if ($SpecifiedCriterion =~ /^le$/i) { if ($CurrentValue > $SpecifiedValue) { $ViolationCount++; last CRITERION; } }
	      if ($SpecifiedCriterion =~ /^ge$/i) { if ($CurrentValue < $SpecifiedValue) { $ViolationCount++; last CRITERION; } }
	      $Nothing = 1;
	    }
	}
	else {
	  CRITERION: {
	      if ($SpecifiedCriterion =~ /^eq$/i) { if ($CurrentValue ne $SpecifiedValue) { $ViolationCount++; last CRITERION; } }
	      if ($SpecifiedCriterion =~ /^le$/i) { if ($CurrentValue gt $SpecifiedValue) { $ViolationCount++; last CRITERION; } }
	      if ($SpecifiedCriterion =~ /^ge$/i) { if ($CurrentValue lt $SpecifiedValue) { $ViolationCount++; last CRITERION; } }
	      $Nothing = 1;
	    }
	}
      }
    }
    if ($ViolationCount <= $OptionsInfo{Violations}) {
      WriteSDFileCmpdString();
      WriteTextFileCmpdData();
    }
  }
}

# Extract data fields by value using regular expression match...
sub ExtractDataFieldsByRegex {
  my($FileIndex) = @_;
  my($Label, $CurrentValue, $SpecifiedRegexCriterion, $SpecifiedRegex, $ViolationCount, $Nothing, @CmpdLines);

  @{$SDFilesInfo{DataLabels}} = @{$SDFilesInfo{AllDataFieldLabels}[$FileIndex]};
  WriteTextFileColLabels();

  CMPDSTRING: while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    @CmpdLines = split "\n", $SDFilesInfo{CmpdString};
    %{$SDFilesInfo{DataFieldValues}} = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

    SetupDataValues();
    $ViolationCount = 0;

    for $Label (@{$OptionsInfo{SpecifiedDataFieldLabels}}) {
      if (exists $SDFilesInfo{DataFieldValues}{$Label}) {
	$CurrentValue = $SDFilesInfo{DataFieldValues}{$Label};
           $SpecifiedRegexCriterion = $OptionsInfo{SpecifiedDataFieldRegexCriteriaMap}{$Label};
           $SpecifiedRegex = $OptionsInfo{SpecifiedDataFieldRegexMap}{$Label};

	if ($OptionsInfo{RegexIgnoreCase}) {
	  CRITERION: {
                 if ($SpecifiedRegexCriterion =~ /^eq$/i) { if ($CurrentValue !~ /$SpecifiedRegex/i) { $ViolationCount++; last CRITERION; } }
                 if ($SpecifiedRegexCriterion =~ /^ne$/i) { if ($CurrentValue =~ /$SpecifiedRegex/i) {  $ViolationCount++; last CRITERION; } }
	      $Nothing = 1;
	    }
	}
	else {
	  CRITERION: {
                 if ($SpecifiedRegexCriterion =~ /^eq$/i) { if ($CurrentValue !~ /$SpecifiedRegex/) { $ViolationCount++; last CRITERION; } }
                 if ($SpecifiedRegexCriterion =~ /^ne$/i) { if ($CurrentValue =~ /$SpecifiedRegex/) {  $ViolationCount++; last CRITERION; } }
	      $Nothing = 1;
	    }
	}
      }
    }
    if ($ViolationCount <= $OptionsInfo{Violations}) {
      WriteSDFileCmpdString();
      WriteTextFileCmpdData();
    }
  }
}

# Extract random compounds...
sub ExtractRandomCompounds {
  my($FileIndex) = @_;
  my($CmpdNum, $CmpdCount, $RandomCycleCount, $RandomIndex, @CmpdLines, %RandomCmpdIndexMap);

  @{$SDFilesInfo{DataLabels}} = @{$SDFilesInfo{AllDataFieldLabels}[$FileIndex]};
  WriteTextFileColLabels();

  $CmpdCount = $SDFilesInfo{CmpdCount}[$FileIndex];
  srand($OptionsInfo{Seed});
  $RandomCycleCount = 0;

  %RandomCmpdIndexMap = ();
  while ($RandomCycleCount <= $CmpdCount && $RandomCycleCount <= $OptionsInfo{NumOfCmpds}) {
    $RandomCycleCount++;
    $RandomIndex = int (rand $CmpdCount) + 1;
    $RandomCmpdIndexMap{$RandomIndex} = $RandomIndex;
  }

  $CmpdNum = 0;
  CMPDSTRING: while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    $CmpdNum++;
    if (!exists $RandomCmpdIndexMap{$CmpdNum}) {
      next CMPDSTRING;
    }

    @CmpdLines = split "\n", $SDFilesInfo{CmpdString};

    WriteSDFileCmpdString();

    if ($OptionsInfo{OutputTextFile}) {
      %{$SDFilesInfo{DataFieldValues}} = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
      SetupDataValues();
      WriteTextFileCmpdData();
    }
  }
}

# Extract mol names...
sub ExtractMolNames {
  my($FileIndex) = @_;
  my($MolName, $NewTextFileRef, @CmpdLines);

  push @{$SDFilesInfo{DataLabels}}, "MolName";
  WriteTextFileColLabels();

  $NewTextFileRef = $SDFilesInfo{NewTextFileRef};
  while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    @CmpdLines = split "\n", $SDFilesInfo{CmpdString};
    $MolName = QuoteAWord(ParseCmpdMolNameLine($CmpdLines[0]), $OptionsInfo{OutQuote});
    print $NewTextFileRef "$MolName\n";
  }
}

# Extract a specific compound record...
sub ExtractRecordNum {
  my($FileIndex) = @_;
  my($CmpdNum, @CmpdLines);

  @{$SDFilesInfo{DataLabels}} = @{$SDFilesInfo{AllDataFieldLabels}[$FileIndex]};
  WriteTextFileColLabels();

  $CmpdNum = 0;

  CMPDSTRING: while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    $CmpdNum++;
    if ($CmpdNum != $OptionsInfo{RecordNum}) {
      next CMPDSTRING;
    }

    @CmpdLines = split "\n", $SDFilesInfo{CmpdString};
    WriteSDFileCmpdString();

    if ($OptionsInfo{OutputTextFile}) {
      %{$SDFilesInfo{DataFieldValues}} = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
      SetupDataValues();
      WriteTextFileCmpdData();
    }
    last CMPDSTRING;
  }
}

# Extract a specific compound records...
sub ExtractRecordNums {
  my($FileIndex) = @_;
  my($CmpdNum, $CmpdCount, @CmpdLines);

  @{$SDFilesInfo{DataLabels}} = @{$SDFilesInfo{AllDataFieldLabels}[$FileIndex]};
  WriteTextFileColLabels();

  $CmpdNum = 0;
  $CmpdCount = 0;

  CMPDSTRING: while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    $CmpdNum++;

    if (exists $OptionsInfo{RecordNums}{$CmpdNum}) {
      $CmpdCount++;
      @CmpdLines = split "\n", $SDFilesInfo{CmpdString};

      WriteSDFileCmpdString();

      if ($OptionsInfo{OutputTextFile}) {
	%{$SDFilesInfo{DataFieldValues}} = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
	SetupDataValues();
	WriteTextFileCmpdData();
      }
    }
    elsif ($CmpdNum > $OptionsInfo{RecordNumsMax} || $CmpdCount >= $OptionsInfo{RecordNumsCount}) {
      last CMPDSTRING;
    }
  }
}


# Extract compounds in a specific record range...
sub ExtractRecordRange {
  my($FileIndex) = @_;
  my($CmpdNum, @CmpdLines);

  @{$SDFilesInfo{DataLabels}} = @{$SDFilesInfo{AllDataFieldLabels}[$FileIndex]};
  WriteTextFileColLabels();

  $CmpdNum = 0;
  CMPDSTRING: while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    $CmpdNum++;

    if ($CmpdNum >= $OptionsInfo{StartRecordNum} && $CmpdNum <= $OptionsInfo{EndRecordNum}) {
      @CmpdLines = split "\n", $SDFilesInfo{CmpdString};

      WriteSDFileCmpdString();

      if ($OptionsInfo{OutputTextFile}) {
	%{$SDFilesInfo{DataFieldValues}} = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
	SetupDataValues();
	WriteTextFileCmpdData();
      }
    }
    elsif ($CmpdNum > $OptionsInfo{EndRecordNum}) {
      last CMPDSTRING;
    }
  }
}

# Extract 2D compound records...
sub Extract2DCmpdRecords {
  my($FileIndex) = @_;
  my(@CmpdLines);

  @{$SDFilesInfo{DataLabels}} = @{$SDFilesInfo{AllDataFieldLabels}[$FileIndex]};
  WriteTextFileColLabels();


  CMPDSTRING: while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    @CmpdLines = split "\n", $SDFilesInfo{CmpdString};
    if (!IsCmpd2D(\@CmpdLines)) {
      next CMPDSTRING;
    }

    WriteSDFileCmpdString();

    if ($OptionsInfo{OutputTextFile}) {
      %{$SDFilesInfo{DataFieldValues}} = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
      SetupDataValues();
      WriteTextFileCmpdData();
    }
  }
}

# Extract 3D compound records...
sub Extract3DCmpdRecords {
  my($FileIndex) = @_;
  my(@CmpdLines);

  @{$SDFilesInfo{DataLabels}} = @{$SDFilesInfo{AllDataFieldLabels}[$FileIndex]};
  WriteTextFileColLabels();


  CMPDSTRING: while ($SDFilesInfo{CmpdString} = ReadCmpdString($SDFilesInfo{InputSDFileRef})) {
    @CmpdLines = split "\n", $SDFilesInfo{CmpdString};
    if (!IsCmpd3D(\@CmpdLines)) {
      next CMPDSTRING;
    }

    WriteSDFileCmpdString();

    if ($OptionsInfo{OutputTextFile}) {
      %{$SDFilesInfo{DataFieldValues}} = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
      SetupDataValues();
      WriteTextFileCmpdData();
    }
  }
}


# Open input and output files...
sub OpenInputAndOutputFiles {
  my($FileIndex) = @_;

  $SDFilesInfo{NewTextFileRef} = undef;
  $SDFilesInfo{NewSDFileRef} = undef;

  if ($OptionsInfo{OutputTextFile} && $OptionsInfo{OutputSDFile}) {
    print "Generating files $SDFilesInfo{NewSDFileName}[$FileIndex] and $SDFilesInfo{NewTextFileName}[$FileIndex]...\n";
  }
  elsif ($OptionsInfo{OutputSDFile}) {
    print "Generating file $SDFilesInfo{NewSDFileName}[$FileIndex]...\n";
  }
  else {
    print "Generating file $SDFilesInfo{NewTextFileName}[$FileIndex]...\n";
  }

  if ($OptionsInfo{OutputSDFile}) {
    open NEWSDFILE, ">$SDFilesInfo{NewSDFileName}[$FileIndex]" or die "Error: Couldn't open $SDFilesInfo{NewSDFileName}[$FileIndex]: $! \n";
    $SDFilesInfo{NewSDFileRef} = \*NEWSDFILE;
  }
  if ($OptionsInfo{OutputTextFile}) {
    open NEWTEXTFILE, ">$SDFilesInfo{NewTextFileName}[$FileIndex]" or die "Error: Couldn't open $SDFilesInfo{NewTextFileName}[$FileIndex]: $! \n";
    $SDFilesInfo{NewTextFileRef} = \*NEWTEXTFILE;
  }

  open SDFILE, "$SDFilesList[$FileIndex]" or die "Error: Couldn't open $SDFilesList[$FileIndex]: $! \n";
  $SDFilesInfo{InputSDFileRef} = \*SDFILE;

}

# Close open input and output files...
sub CloseInputAndOutputFiles {
  if ($SDFilesInfo{NewSDFileRef}) {
    close $SDFilesInfo{NewSDFileRef};
  }
  if ($SDFilesInfo{NewTextFileRef}) {
    close $SDFilesInfo{NewTextFileRef};
  }

  if ($SDFilesInfo{InputSDFileRef}) {
    close $SDFilesInfo{InputSDFileRef};
  }

  $SDFilesInfo{NewTextFileRef} = undef;
  $SDFilesInfo{NewSDFileRef} = undef;
  $SDFilesInfo{InputSDFileRef} = undef;
}

# Write out column labels for text file...
sub WriteTextFileColLabels {
  my($ColLabelsLine, $NewTextFileRef);

  if (!$OptionsInfo{OutputTextFile}) {
    return;
  }
  $NewTextFileRef = $SDFilesInfo{NewTextFileRef};

  if ($OptionsInfo{OutputStrDataString}) {
    # Append structure data string label...
    my(@DataLabels);

    @DataLabels = ();
    push @DataLabels, @{$SDFilesInfo{DataLabels}};
    push @DataLabels, "StructureDataString";

    $ColLabelsLine = JoinWords(\@DataLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  }
  else {
    $ColLabelsLine = JoinWords(\@{$SDFilesInfo{DataLabels}}, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  }
  print $NewTextFileRef "$ColLabelsLine\n";
}

# Setup values for data fields...
sub SetupDataValues {
  @{$SDFilesInfo{DataValues}} = map { exists $SDFilesInfo{DataFieldValues}{$_} ? $SDFilesInfo{DataFieldValues}{$_} : "" } @{$SDFilesInfo{DataLabels}};
}

# Write out structure data and specific data fields to SD file...
sub WriteSDFileCmpdData {
  my($MolString, $Count, $NewSDFileRef);

  if (!$OptionsInfo{OutputSDFile}) {
    return;
  }

  $NewSDFileRef = $SDFilesInfo{NewSDFileRef};

  ($MolString) = split "M  END", $SDFilesInfo{CmpdString};
  $MolString .= "M  END";
  print $NewSDFileRef "$MolString\n";

  for $Count (0 .. $#{$SDFilesInfo{DataLabels}}) {
    print $NewSDFileRef ">  <$SDFilesInfo{DataLabels}[$Count]>\n$SDFilesInfo{DataValues}[$Count]\n\n";
  }
  print $NewSDFileRef "\$\$\$\$\n";
}

# Write out compound string...
sub WriteSDFileCmpdString {
  my($NewSDFileRef);

  if (!$OptionsInfo{OutputSDFile}) {
    return;
  }

  $NewSDFileRef = $SDFilesInfo{NewSDFileRef};
  print $NewSDFileRef "$SDFilesInfo{CmpdString}\n";
}

# Write out data for text file...
sub WriteTextFileCmpdData {
  my($DataValuesLine, $NewTextFileRef);

  if (!$OptionsInfo{OutputTextFile}) {
    return;
  }

  $NewTextFileRef = $SDFilesInfo{NewTextFileRef};
  $DataValuesLine = JoinWords(\@{$SDFilesInfo{DataValues}}, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});

  # Handle multiple lines data values for data fields by joining 'em using semicolons...
  if ($DataValuesLine =~ /\n/) {
    $DataValuesLine =~ s/\n/;/g;
  }

  if ($OptionsInfo{OutputStrDataString}) {
    # Append structure data string...
    my($StrDataString, $OutQuoteValue, $OutDelim, $StrDataStringDelimiter);

    if ($OptionsInfo{StrDataStringWithFields}) {
      $StrDataString = $SDFilesInfo{CmpdString};
    }
    else {
      ($StrDataString) = split "M  END", $SDFilesInfo{CmpdString};
      $StrDataString .= "M  END";
    }
    $StrDataStringDelimiter = $OptionsInfo{StrDataStringDelimiter};
    $StrDataString =~ s/\n/$StrDataStringDelimiter/g;

    $OutDelim = $OptionsInfo{OutDelim};
    $OutQuoteValue = $OptionsInfo{OutQuote} ? "\"" : "";

    print $NewTextFileRef "$DataValuesLine${OutDelim}${OutQuoteValue}${StrDataString}${OutQuoteValue}\n";
  }
  else {
    print $NewTextFileRef "$DataValuesLine\n";
  }
}

# Retrieve information about input SD files...
sub RetrieveSDFilesInfo {
  my($SDFile, $Index, $FileDir, $FileExt, $FileName, $NewFileName, $NewSDFileName, $NewTextFileName, $CmpdCount);

  %SDFilesInfo = ();

  @{$SDFilesInfo{FileOkay}} = ();
  @{$SDFilesInfo{CmpdCount}} = ();
  @{$SDFilesInfo{NewTextFileName}} = ();
  @{$SDFilesInfo{NewSDFileName}} = ();

  @{$SDFilesInfo{AllDataFieldLabels}} = ();
  @{$SDFilesInfo{CommonDataFieldLabels}} = ();

  FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFile = $SDFilesList[$Index];

    $SDFilesInfo{FileOkay}[$Index] = 0;

    $SDFilesInfo{CmpdCount}[$Index] = 0;
    $SDFilesInfo{NewTextFileName}[$Index] = "";
    $SDFilesInfo{NewSDFileName}[$Index] = "";

    @{$SDFilesInfo{AllDataFieldLabels}[$Index]} = ();
    @{$SDFilesInfo{CommonDataFieldLabels}[$Index]} = ();

    if (!(-e $SDFile)) {
      warn "Warning: Ignoring file $SDFile: It doesn't exist\n";
      next FILELIST;
    }

    if (!CheckFileType($SDFile, "sd sdf")) {
      warn "Warning: Ignoring file $SDFile: It's not a SD file\n";
      next FILELIST;
    }

    # Generate appropriate name for the new output file.
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($SDFile);
    $NewFileName = $FileName;
    $NewFileName = $FileName  . $OptionsInfo{FileNameMode};
    if ($OptionsInfo{OutFileRoot} && (@SDFilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($OptionsInfo{OutFileRoot});
      if ($RootFileName && $RootFileExt) {
	$NewFileName = $RootFileName;
      }
      else {
	$NewFileName = $OptionsInfo{OutFileRoot};
      }
    }
    $NewSDFileName = $NewFileName . ".$OptionsInfo{SDFileExt}";
    $NewTextFileName = $NewFileName . ".$OptionsInfo{TextFileExt}";

    if ($OptionsInfo{OutputSDFile}) {
      if (lc($NewSDFileName) eq lc($SDFile)) {
	warn "Warning: Ignoring input file $SDFile: Same output, $NewSDFileName, and input file names.\n";
	print "Specify a different name using \"-r --root\" option or use default name.\n";
	next FILELIST;
      }
    }

    if (!$OptionsInfo{Overwrite}) {
      if ($OptionsInfo{OutputSDFile}) {
	if (-e $NewSDFileName) {
	  warn "Warning: Ignoring file $SDFile: New file, $NewSDFileName, already exists\n";
	  next FILELIST;
	}
      }
      if ($OptionsInfo{OutputTextFile}) {
	if (-e $NewTextFileName) {
	  warn "Warning: Ignoring file $SDFile: New file, $NewTextFileName, already exists\n";
	  next FILELIST;
	}
      }
    }

    if (!open SDFILE, "$SDFile") {
      warn "Warning: Ignoring file $SDFile: Couldn't open it: $! \n";
      next FILELIST;
    }

    my($CountCmpds, $CollectDataFields);
    my($CmpdString, @CmpdLines, @DataFieldLabels, %DataFieldLabelsMap,@CommonDataFieldLabels);

    $CountCmpds = ($OptionsInfo{Mode} =~ /^randomcmpds$/i) ? 1 : 0;

    $CollectDataFields = (($OptionsInfo{Mode} =~ /^(alldatafields|commondatafields|randomcmpds)$/i && $OptionsInfo{OutputTextFile}) || ($OptionsInfo{Mode} =~ /^(datafieldsbyvalue|datafieldsbyregex)$/i  && $OptionsInfo{OutputTextFile}) || ($OptionsInfo{Mode} =~ /^datafieldbylist$/i  && $OptionsInfo{OutputTextFile}) || ($OptionsInfo{Mode} =~ /^datafielduniquebylist$/i  && $OptionsInfo{OutputTextFile}) || ($OptionsInfo{Mode} =~ /^datafieldnotbylist$/i  && $OptionsInfo{OutputTextFile}) || ($OptionsInfo{Mode} =~ /^recordnum$/i && $OptionsInfo{OutputTextFile}) || ($OptionsInfo{Mode} =~ /^recordnums$/i && $OptionsInfo{OutputTextFile}) || ($OptionsInfo{Mode} =~ /^recordrange$/i && $OptionsInfo{OutputTextFile})) ? 1 : 0;

    $CmpdCount = 0;
    if ($CountCmpds || $CollectDataFields) {
      @DataFieldLabels = ();
      @CommonDataFieldLabels = ();
      %DataFieldLabelsMap = ();
      CMPDSTRING: while ($CmpdString = ReadCmpdString(\*SDFILE)) {
	$CmpdCount++;
	if ($OptionsInfo{Mode} =~ /^recordnum$/i) {
	  if ($CmpdCount == $OptionsInfo{RecordNum}) {
	    @CmpdLines = split "\n", $CmpdString;
	    @DataFieldLabels = GetCmpdDataHeaderLabels(\@CmpdLines);
	    last CMPDSTRING;
	  }
	}
	if ($CollectDataFields) {
	  my($Label);
	  @CmpdLines = split "\n", $CmpdString;
	  # Process compound data header labels and figure out which ones are present for
	  # all the compounds...
	  if (@DataFieldLabels) {
	    my (@CmpdDataFieldLabels) = GetCmpdDataHeaderLabels(\@CmpdLines);
	    my(%CmpdDataFieldLabelsMap) = ();
	    # Setup a map for the current labels...
	    for $Label (@CmpdDataFieldLabels) {
	      $CmpdDataFieldLabelsMap{$Label} = "PresentInSome";
	    }
	    # Check the presence old labels for this compound; otherwise, mark 'em new...
	    for $Label (@DataFieldLabels) {
	      if (!$CmpdDataFieldLabelsMap{$Label}) {
		$DataFieldLabelsMap{$Label} = "PresentInSome";
	      }
	    }
	    # Check the presence this compound in the old labels; otherwise, add 'em...
	    for $Label (@CmpdDataFieldLabels ) {
	      if (!$DataFieldLabelsMap{$Label}) {
		# It's a new label...
		push @DataFieldLabels, $Label;
		$DataFieldLabelsMap{$Label} = "PresentInSome";
	      }
	    }
	  }
	  else {
	    # Get the initial label set and set up a map...
	    @DataFieldLabels = GetCmpdDataHeaderLabels(\@CmpdLines);
	    for $Label (@DataFieldLabels) {
	      $DataFieldLabelsMap{$Label} = "PresentInAll";
	    }
	  }
	  # Identify the common data field labels...
	  if ($Options{mode} =~ /^commondatafields$/i) {
	    @CommonDataFieldLabels = ();
	    for $Label (@DataFieldLabels) {
	      if ($DataFieldLabelsMap{$Label} eq "PresentInAll") {
		push @CommonDataFieldLabels, $Label;
	      }
	    }
	  }
	}
      }
    }

    $SDFilesInfo{FileOkay}[$Index] = 1;

    $SDFilesInfo{NewTextFileName}[$Index] = $NewTextFileName;
    $SDFilesInfo{NewSDFileName}[$Index] = $NewSDFileName;

    $SDFilesInfo{CmpdCount}[$Index] = $CmpdCount;

    push @{$SDFilesInfo{AllDataFieldLabels}[$Index]}, @DataFieldLabels;
    push @{$SDFilesInfo{CommonDataFieldLabels}[$Index]}, @CommonDataFieldLabels;

    close SDFILE;
  }
}

# Process options...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{InDelim} = "\,";
  if ($Options{indelim} =~ /^semicolon$/i) {
    $OptionsInfo{InDelim} = "\;";
  }
  elsif ($Options{indelim} =~ /^tab$/i) {
    $OptionsInfo{InDelim} = "\t";
  }

  $OptionsInfo{OutDelim} = "\,";
  if ($Options{outdelim} =~ /^semicolon$/i) {
    $OptionsInfo{OutDelim} = "\;";
  }
  elsif ($Options{outdelim} =~ /^tab$/i) {
    $OptionsInfo{OutDelim} = "\t";
  }

  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^yes$/i) ? 1 : 0;

  $OptionsInfo{RegexIgnoreCase} = ($Options{regexignorecase} =~ /^yes$/i) ? 1 : 0;

  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : undef;
  $OptionsInfo{Overwrite} = $Options{overwrite} ? $Options{overwrite} : undef;

  $OptionsInfo{NumOfCmpds} = $Options{numofcmpds};

  $OptionsInfo{ValueComparisonMode} = $Options{valuecomparisonmode};
  $OptionsInfo{NumericalComparison} = ($Options{valuecomparisonmode} =~ /^Numeric$/i) ? 1 : 0;

  $OptionsInfo{Violations} = $Options{violations};
  $OptionsInfo{Seed} = $Options{seed};


  if ($Options{mode} =~ /^(datafields|datafieldsbyregex|datafieldsbyvalue|datafieldbylist|datafielduniquebylist|datafieldnotbylist)$/i) {
    if ($Options{datafields} || $Options{datafieldsfile}) {
      if ($Options{datafields} && $Options{datafieldsfile}) {
	die "Error: For \"-m --mode\" option values of datafields, datafieldsbyvalue, datafieldsbyregex, datafieldbylist, datafielduniquebylist, or datafieldnotbylist specify only one of the \"-d --datafields\" or \"--datafieldsfile\" option.\n";
      }
    }
    else {
      die "Error: For \"-m --mode\" option values of datafields, datafieldsbyvalue, datafieldsbyregex, datafieldbylist, datafielduniquebylist, or datafieldnotbylist specify one of the \"-d --datafields\" or \"--datafieldsfile\" option.\n";
    }
  }
  $OptionsInfo{DataFields} = $Options{datafields} ? $Options{datafields} : undef;
  $OptionsInfo{DataFieldsFile} = $Options{datafieldsfile} ? $Options{datafieldsfile} : undef;

  $OptionsInfo{RecordNum} = 0; $OptionsInfo{StartRecordNum} = 0; $OptionsInfo{EndRecordNum} = 0;

  %{$OptionsInfo{RecordNums}} = ();
  $OptionsInfo{RecordNumsMin} = 0; $OptionsInfo{RecordNumsMax} = 0; $OptionsInfo{RecordNumsCount} = 0;

  $OptionsInfo{Record} = $Options{record} ? $Options{record} : undef;

  if ($Options{mode} =~ /^(recordnum|recordnums|recordrange)$/i) {
    if ($Options{record}) {
      my($Record, @RecordSplit);

      $Record = $Options{record};
      $Record =~ s/ //g;

      @RecordSplit = split ",", $Record;

      if ($Options{mode} =~ /^recordnum$/i ) {
	if (@RecordSplit == 1) {
	  $OptionsInfo{RecordNum} = $RecordSplit[0];
	  if ($OptionsInfo{RecordNum} <= 0) {
	    die "Error: The value specified, $OptionsInfo{RecordNum},  for option \"--records\" is not valid. Allowed values: > 0 \n";
	  }
	}
	else {
	  die "Error: Invalid number of values, ", scalar(@RecordSplit), ", specified using \"--record\" option: only 1 value is allowed.\n";
	}
      }
      elsif ($Options{mode} =~ /^recordnums$/i ) {
	my($RecordNum, $RecordCount, @SortedRecordSplit);

	@SortedRecordSplit = sort { $a <=> $b } @RecordSplit;

	$RecordCount = 0;
	RECORDNUM: for $RecordNum (@SortedRecordSplit) {
	  if (exists $OptionsInfo{RecordNums}{$RecordNum}) {
	    next RECORDNUM;
	  }
	  $RecordCount++;
	  $OptionsInfo{RecordNums}{$RecordNum} = $RecordNum;
	}
	$OptionsInfo{RecordNumsCount} = $RecordCount;
	$OptionsInfo{RecordNumsMin} = $SortedRecordSplit[0];
	$OptionsInfo{RecordNumsMax} = $SortedRecordSplit[$#SortedRecordSplit];
      }
      else {
	if (@RecordSplit == 2) {
	  $OptionsInfo{StartRecordNum} = $RecordSplit[0];
	  $OptionsInfo{EndRecordNum} = $RecordSplit[1];
	  if ($OptionsInfo{StartRecordNum} <= 0 || $OptionsInfo{EndRecordNum} <= 0) {
	    die "Error: The value pair specified, $Options{record},  for option \"--records\" is not valid. Allowed values: > 0 \n";
	  }
	}
	else {
	  die "Error: Invalid number of values, ", scalar(@RecordSplit), ", specified using \"--record\" option: only 2 values is allowed.\n";
	}
	if ($OptionsInfo{StartRecordNum} > $OptionsInfo{EndRecordNum}) {
	  die "Error: Start record number, $OptionsInfo{StartRecordNum}, must be smaller than end record number, $OptionsInfo{EndRecordNum}.\nSpecify different values using \"--record\" option.\n";
	}
      }
    }
    else {
      die "Error: For \"-m --mode\" option values recordnum, recordnums or recordrange, specify \"--record\" option value.\n";
    }
  }

  @{$OptionsInfo{SpecifiedDataFieldLabels}} = ();

  my(@Words, $Line, $Value);
  if ($Options{mode} =~ /^datafields$/i) {
    @{$OptionsInfo{SpecifiedDataFieldLabels}} = ();
    if ($Options{datafields}) {
      @{$OptionsInfo{SpecifiedDataFieldLabels}} = split $OptionsInfo{InDelim}, $Options{datafields};
    }
    elsif ($Options{datafieldsfile}) {
      open DATAFIELDSFILE, "$Options{datafieldsfile}" or die "Error: Couldn't open $Options{datafieldsfile}: $! \n";
      while ($Line = GetTextLine(\*DATAFIELDSFILE)) {
	@Words = quotewords($OptionsInfo{InDelim}, 0, $Line);
	if (@Words) {
	  push @{$OptionsInfo{SpecifiedDataFieldLabels}}, @Words;
	}
      }
      close DATAFIELDSFILE;
    }
  }
  elsif ($Options{mode} =~ /^datafieldsbyvalue$/i) {
    my(@DataFieldsByValueTriplets);
    @DataFieldsByValueTriplets = ();
    if ($Options{datafields}) {
      @DataFieldsByValueTriplets = split $OptionsInfo{InDelim}, $Options{datafields};
    }
    elsif ($Options{datafieldsfile}) {
      open DATAFIELDSFILE, "$Options{datafieldsfile}" or die "Error: Couldn't open $Options{datafieldsfile}: $! \n";
      while ($Line = GetTextLine(\*DATAFIELDSFILE)) {
	@Words = quotewords($OptionsInfo{InDelim}, 0, $Line);
	if (@Words) {
	  push @DataFieldsByValueTriplets, @Words;
	}
      }
      close DATAFIELDSFILE;
    }
    if ((@DataFieldsByValueTriplets % 3)) {
      if ($Options{datafields}) {
	die "Error: Triplets not found in values specified by \"-d --datafields\" option\n";
      }
      elsif ($Options{datafieldsfile}) {
	die "Error: Triplets not found in values specified by \"--datafieldsfile\" option\n";
      }
    }
    my($Index, $Label, $Value, $Criterion);

    @{$OptionsInfo{SpecifiedDataFieldLabels}} = ();
    %{$OptionsInfo{SpecifiedDataFieldValuesMap}} = ();
    %{$OptionsInfo{SpecifiedDataFieldCriteriaMap}} = ();

    for ($Index = 0; $Index < @DataFieldsByValueTriplets; $Index = $Index + 3) {
      $Label = $DataFieldsByValueTriplets[$Index];
      $Value = $DataFieldsByValueTriplets[$Index + 1];
      $Criterion = $DataFieldsByValueTriplets[$Index + 2];

      if ($Criterion =~ /^(eq|le|ge)$/i) {
	push @{$OptionsInfo{SpecifiedDataFieldLabels}}, $Label;
	$OptionsInfo{SpecifiedDataFieldValuesMap}{$Label} = $Value;
	$OptionsInfo{SpecifiedDataFieldCriteriaMap}{$Label} = $Criterion;
      }
      else {
	warn "Warning: Ignoring triplet value, $Label $Value $Criterion , specified using \"-d --datafields\" or \"--datafieldsfile\" option: Invalid criterion value: $Criterion\n";
      }
    }
  }
  elsif ($Options{mode} =~ /^datafieldsbyregex$/i) {
    my(@DataFieldsByRegexTriplets);

    @DataFieldsByRegexTriplets = ();
    if ($Options{datafields}) {
      @DataFieldsByRegexTriplets = quotewords($OptionsInfo{InDelim}, 0, $Options{datafields});
    }
    elsif ($Options{datafieldsfile}) {
      open DATAFIELDSFILE, "$Options{datafieldsfile}" or die "Error: Couldn't open $Options{datafieldsfile}: $! \n";
      while ($Line = GetTextLine(\*DATAFIELDSFILE)) {
          @Words = quotewords($OptionsInfo{InDelim}, 0, $Line);
          if (@Words) {
            push @DataFieldsByRegexTriplets, @Words;
          }
      }
      close DATAFIELDSFILE;
    }
    if ((@DataFieldsByRegexTriplets % 3)) {
      if ($Options{datafields}) {
          die "Error: Triplet not found in values specified by \"-d --datafields\" option\n";
      }
      elsif ($Options{datafieldsfile}) {
          die "Error: Triplet not found in values specified by \"--datafieldsfile\" option\n";
      }
    }

    my($Index, $Label, $Value, $Criterion);

    @{$OptionsInfo{SpecifiedDataFieldLabels}} = ();
    %{$OptionsInfo{SpecifiedDataFieldRegexMap}} = ();
    %{$OptionsInfo{SpecifiedDataFieldRegexCriteriaMap}} = ();

    for ($Index = 0; $Index < @DataFieldsByRegexTriplets; $Index = $Index + 3) {
      $Label = $DataFieldsByRegexTriplets[$Index];
      $Value = $DataFieldsByRegexTriplets[$Index + 1];
      $Criterion = $DataFieldsByRegexTriplets[$Index + 2];

      if ($Criterion =~ /^(eq|ne)$/i) {
          push @{$OptionsInfo{SpecifiedDataFieldLabels}}, $Label;
          $OptionsInfo{SpecifiedDataFieldRegexMap}{$Label} = $Value;
          $OptionsInfo{SpecifiedDataFieldRegexCriteriaMap}{$Label} = $Criterion;
      }
      else {
          warn "Warning: Ignoring triplet value, $Label $Value $Criterion , specified using \"-d --datafields\" or \"--datafieldsfile\" option: Invalid criterion value: $Criterion; Supported values: eq or ne\n";
      }
    }
  }
  elsif ($Options{mode} =~ /^(datafieldbylist|datafielduniquebylist|datafieldnotbylist)$/i) {
    my($Index, @DataFieldAndValuesList);
    if ($Options{datafields}) {
      @DataFieldAndValuesList = split $OptionsInfo{InDelim}, $Options{datafields};
    }
    elsif ($Options{datafieldsfile}) {
      open DATAFIELDSFILE, "$Options{datafieldsfile}" or die "Error: Couldn't open $Options{datafieldsfile}: $! \n";
      while ($Line = GetTextLine(\*DATAFIELDSFILE)) {
	@Words = quotewords($OptionsInfo{InDelim}, 0, $Line);
	if (@Words) {
	  push @DataFieldAndValuesList, @Words;
	}
      }
      close DATAFIELDSFILE;
    }
    if (@DataFieldAndValuesList < 2) {
      if ($Options{datafields}) {
	die "Error: Invalid number of values specified by \"-d --datafields\" option\n";
      }
      elsif ($Options{datafieldsfile}) {
	die "Error: Invalid number values specified by \"--datafieldsfile\" option\n";
      }
    }

    $OptionsInfo{SpecifiedDataFieldLabel} = $DataFieldAndValuesList[0];
    $OptionsInfo{SpecifiedDataFieldValuesCount} = @DataFieldAndValuesList - 1;
    %{$OptionsInfo{SpecifiedDataFieldValues}} = ();

    for ($Index = 1; $Index < @DataFieldAndValuesList; $Index++) {
      $Value = $DataFieldAndValuesList[$Index];
      $OptionsInfo{SpecifiedDataFieldValues}{$Value} = "NotFound";
    }
  }

  $OptionsInfo{SDFileExt} = "sdf";
  $OptionsInfo{TextFileExt} = "csv";

  if ($Options{outdelim} =~ /^tab$/i) {
    $OptionsInfo{TextFileExt} = "tsv";
  }

  if ($Options{mode} =~ /^(alldatafields|molnames)$/i) {
    $OptionsInfo{OutputSDFile} = 0;
    $OptionsInfo{OutputTextFile} = 1;
  }
  else {
    $OptionsInfo{OutputSDFile} = ($Options{output} =~ /^(SD|both)$/i) ? 1 : 0;
    $OptionsInfo{OutputTextFile} = ($Options{output} =~ /^(text|both)$/i) ? 1 : 0;
  }

  $OptionsInfo{StrDataString} = $Options{strdatastring};
  $OptionsInfo{OutputStrDataString} = ($Options{strdatastring} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{StrDataStringDelimiter} = $Options{strdatastringdelimiter};

  if (IsEmpty($Options{strdatastringdelimiter})) {
    die "Error: No value specified for \"--StrDataStringDelimiter\" option.\n";
  }
  $OptionsInfo{StrDataStringMode} = $Options{strdatastringmode};
  $OptionsInfo{StrDataStringWithFields} = $Options{strdatastringmode} =~ /^StrAndDataFields$/i ? 1 : 0;

  MODE: {
    if ($Options{mode} =~ /^alldatafields$/i) { $OptionsInfo{FileNameMode} = "AllDataDields"; last MODE; }
    if ($Options{mode} =~ /^commondatafields$/i) { $OptionsInfo{FileNameMode} = "CommonDataDields"; last MODE; }
    if ($Options{mode} =~ /^datafields$/i) { $OptionsInfo{FileNameMode} = "SpecifiedDataFields"; last MODE; }
    if ($Options{mode} =~ /^datafieldsbyvalue$/i) { $OptionsInfo{FileNameMode} = "SpecifiedDataFieldsByValue"; last MODE; }
    if ($Options{mode} =~ /^datafieldsbyregex$/i) { $OptionsInfo{FileNameMode} = "SpecifiedDataFieldsByRegex"; last MODE; }
    if ($Options{mode} =~ /^datafieldbylist$/i) { $OptionsInfo{FileNameMode} = "SpecifiedDataField"; last MODE; }
    if ($Options{mode} =~ /^datafielduniquebylist$/i) { $OptionsInfo{FileNameMode} = "SpecifiedUniqueDataField"; last MODE; }
    if ($Options{mode} =~ /^datafieldnotbylist$/i) { $OptionsInfo{FileNameMode} = "SpecifiedDataFieldNotByList"; last MODE; }
    if ($Options{mode} =~ /^molnames$/i) { $OptionsInfo{FileNameMode} = "MolName"; last MODE; }
    if ($Options{mode} =~ /^randomcmpds$/i) { $OptionsInfo{FileNameMode} = "RandomCmpds"; last MODE; }
    if ($Options{mode} =~ /^recordnum$/i) { $OptionsInfo{FileNameMode} = "RecordNum$OptionsInfo{RecordNum}"; last MODE; }
    if ($Options{mode} =~ /^recordnums$/i) { $OptionsInfo{FileNameMode} = "RecordNums"; last MODE; }
    if ($Options{mode} =~ /^recordrange$/i) { $OptionsInfo{FileNameMode} = "RecordNum$OptionsInfo{StartRecordNum}" . "To" . "$OptionsInfo{EndRecordNum}"; last MODE; }
    if ($Options{mode} =~ /^2dcmpdrecords$/i) { $OptionsInfo{FileNameMode} = "2DCmpdRecords"; last MODE; }
    if ($Options{mode} =~ /^3dcmpdrecords$/i) { $OptionsInfo{FileNameMode} = "3DCmpdRecords"; last MODE; }
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: alldatafields, commondatafields, datafields, datafieldsbyvalue, datafieldbylist, datafielduniquebylist, , datafieldnotbylist, molnames, randomcmpds, recordnum, recordnums, recordrange, 2dcmpdrecords, 3dcmpdrecords\n";
  }

}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{numofcmpds} = 1;
  $Options{mode} = "alldatafields";
  $Options{indelim} = "comma";
  $Options{outdelim} = "comma";
  $Options{output} = "SD";
  $Options{quote} = "yes";
  $Options{regexignorecase} = "yes";
  $Options{valuecomparisonmode} = "numeric";
  $Options{violations} = 0;
  $Options{seed} = 123456789;

  $Options{strdatastring} = "no";
  $Options{strdatastringdelimiter} = "|";
  $Options{strdatastringmode} = "StrOnly";

  if (!GetOptions(\%Options, "help|h", "datafields|d=s", "datafieldsfile=s", "indelim=s", "mode|m=s", "numofcmpds|n=i", "outdelim=s", "output=s", "overwrite|o", "quote|q=s", "regexignorecase=s", "record=s", "root|r=s", "seed|s=i", "strdatastring=s", "strdatastringdelimiter=s", "strdatastringmode=s", "valuecomparisonmode=s", "violations|v=i", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{numofcmpds} < 1) {
    die "Error: The value specified, $Options{numofcmpds},  for option \"-n --numofcmpds\" is not valid. Allowed values: >= 1 \n";
  }
  if ($Options{valuecomparisonmode} !~ /^(Numeric|Alphanumeric)$/i) {
    die "Error: The value specified, $Options{valuecomparisonmode}, for option \"--ValueComparisonMode\" is not valid. Allowed values: Numeric or Alphanumeric\n";
  }
  if ($Options{violations} < 0) {
    die "Error: The value specified, $Options{violations},  for option \"-v --violations\" is not valid. Allowed values: >= 0 \n";
  }
  if ($Options{mode} !~ /^(alldatafields|commondatafields|datafields|datafieldsbyvalue|datafieldsbyregex|datafieldbylist|datafielduniquebylist|datafieldnotbylist|molnames|randomcmpds|recordnum|recordnums|recordrange|2dcmpdrecords|3dcmpdrecords)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: alldatafields, commondatafields, datafields, datafieldsbyvalue, datafieldbylist, datafielduniquebylist, datafieldnotbylist, molnames, randomcmpds, recordnum, recordnums, recordrange, 2dcmpdrecords, 3dcmpdrecords\n";
  }
  if ($Options{output} !~ /^(SD|text|both)$/i) {
    die "Error: The value specified, $Options{output}, for option \"--output\" is not valid. Allowed values: SD, text, or both\n";
  }
  if ($Options{indelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{indelim}, for option \"--indelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if ($Options{outdelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{outdelim}, for option \"--outdelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if ($Options{quote} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{regexignorecase} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{regexignorecase}, for option \"--regexignorecase\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{strdatastring} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{strdatastring}, for option \"--StrDataString\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{strdatastringmode} !~ /^(StrOnly|StrAndDataFields)$/i) {
    die "Error: The value specified, $Options{strdatastringmode}, for option \"--StrDataStringMode\" is not valid. Allowed values: StrOnly or StrAndDataFields\n";
  }
}

__END__

=head1 NAME

ExtractFromSDFiles.pl - Extract specific data from SDFile(s)

=head1 SYNOPSIS

ExtractFromSDFiles.pl SDFile(s)...

ExtractFromSDFiles.pl [B<-h, --help>]
[B<-d, --datafields> "fieldlabel,..." | "fieldlabel,value,criteria..." | "fieldlabel,value,value..."]
[B<--datafieldsfile> filename] [B<--indelim> comma | tab | semicolon] [B<-m, --mode> alldatafields |
commondatafields | | datafieldnotbylist | datafields | datafieldsbyvalue | datafieldsbyregex | datafieldbylist |
datafielduniquebylist | molnames | randomcmpds | recordnum | recordnums | recordrange | 2dcmpdrecords |
3dcmpdrecords ] [B<-n, --numofcmpds> number] [B<--outdelim> comma | tab | semicolon]
[B<--output> SD | text | both] [B<-o, --overwrite>] [B<-q, --quote> yes | no]
[B<--record> recnum | startrecnum,endrecnum] B<--RegexIgnoreCase> I<yes or no>
[B<-r, --root> rootname] [B<-s, --seed> number] [B<--StrDataString> yes | no]
[B<--StrDataStringDelimiter> text] [B<--StrDataStringMode> StrOnly | StrAndDataFields]
[B<--ValueComparisonMode> I<Numeric | Alphanumeric>]
[B<-v, --violations-> number] [B<-w, --workingdir> dirname] SDFile(s)...

=head1 DESCRIPTION

Extract specific data from I<SDFile(s)> and generate appropriate SD or CSV/TSV text
file(s). The structure data from SDFile(s) is not transferred to CSV/TSV text file(s).
Multiple SDFile names are separated by spaces. The valid file extensions are I<.sdf>
and I<.sd>. All other file names are ignored. All the SD files in a current directory
can be specified either by I<*.sdf> or the current directory name.

=head1 OPTIONS

=over 4

=item B<-h, --help>

Print this help message.

=item B<-d, --datafields> I<"fieldlabel,..." | "fieldlabel,value,criteria..." | "fieldlabel,value,value,...">

This value is mode specific. In general, it's a list of comma separated data field labels
and associated mode specific values.

For I<datafields> mode, input value format is: I<fieldlabel,...>. Examples:

    Extreg
    Extreg,CompoundName,ID

For I<datafieldsbyvalue> mode, input value format contains these triplets:
I<fieldlabel,value, criteria...>. Possible values for criteria: I<le, ge or eq>.
The values of B<--ValueComparisonMode> indicates whether values are
compared numerical or string comarison operators. Default is to consider
data field values as numerical values and use numerical comparison operators.
Examples:

    MolWt,450,le
    MolWt,450,le,LogP,5,le,SumNumNO,10,le,SumNHOH,5,le

For I<datafieldsbyregex> mode, input value format contains these triplets:
I<fieldlabel,regex, criteria...>. I<regex> corresponds to any valid regular expression
and is used to match the values for specified I<fieldlabel>. Possible values for criteria:
 I<eq or ne>. During I<eq> and I<ne> values, data field label value is matched with
regular expression using =~ and !~ respectively. B<--RegexIgnoreCase> option
value is used to determine whether to ignore letter upper/lower case during
regular expression match. Examples:

    Name,ol,eq
    Name,'^pat',ne

For I<datafieldbylist> and I<datafielduniquebylist> mode, input value format is:
I<fieldlabel,value1,value2...>. This is equivalent to I<datafieldsbyvalue> mode with
this input value format:I<fieldlabel,value1,eq,fieldlabel,value2,eq,...>. For
I<datafielduniquebylist> mode, only unique compounds identified by first occurrence
of I<value> associated with I<fieldlabel> in I<SDFile(s)> are kept; any subsequent compounds
are simply ignored.

For I<datafieldnotbylist> mode, input value format is: I<fieldlabel,value1,value2...>. In this
mode, the script behaves exactly opposite of I<datafieldbylist> mode, and only those compounds
are extracted whose data field values don't match any specified data field value.

=item B<--datafieldsfile> I<filename>

Filename which contains various mode specific values. This option provides a way
to specify mode specific values in a file instead of entering them on the command
line using B<-d --datafields>.

For I<datafields> mode, input file lines contain comma delimited field labels:
I<fieldlabel,...>. Example:

    Line 1:MolId
    Line 2:"Extreg",CompoundName,ID

For I<datafieldsbyvalue> mode, input file lines contains these comma separated triplets:
I<fieldlabel,value, criteria>. Possible values for criteria: I<le, ge or eq>. Examples:

    Line 1:MolWt,450,le

    Line 1:"MolWt",450,le,"LogP",5,le,"SumNumNO",10,le,"SumNHOH",5,le

    Line 1:MolWt,450,le
    Line 2:"LogP",5,le
    Line 3:"SumNumNO",10,le
    Line 4: SumNHOH,5,le

For I<datafieldbylist> and I<datafielduniquebylist> mode, input file line format is:

    Line 1:fieldlabel;
    Subsequent lines:value1,value2...

For I<datafieldbylist>, I<datafielduniquebylist>, and I<datafieldnotbylist> mode, input file
line format is:

    Line 1:fieldlabel;
    Subsequent lines:value1,value2...

For I<datafielduniquebylist> mode, only unique compounds identified by first occurrence
of I<value> associated with I<fieldlabel> in I<SDFile(s)> are kept; any subsequent compounds
are simply ignored. Example:

    Line 1: MolID
    Subsequent Lines:
    907508
    832291,4642
    "1254","907303"

=item B<--indelim> I<comma | tab | semicolon>

Delimiter used to specify text values for B<-d --datafields> and B<--datafieldsfile> options.
Possible values: I<comma, tab, or semicolon>. Default value: I<comma>.

=item B<-m, --mode> I<alldatafields | commondatafields | datafields | datafieldsbyvalue | datafieldsbyregex | datafieldbylist | datafielduniquebylist |  datafieldnotbylist | molnames | randomcmpds | recordnum | recordnums | recordrange | 2dcmpdrecords | 3dcmpdrecords>

Specify what to extract from I<SDFile(s)>. Possible values: I<alldatafields, commondatafields,
datafields, datafieldsbyvalue, datafieldsbyregex, datafieldbylist, datafielduniquebylist, datafieldnotbylist,
molnames, randomcmpds, recordnum, recordnums, recordrange, 2dcmpdrecords, 3dcmpdrecords>.
Default value: I<alldatafields>.

For I<alldatafields> and I<molnames> mode, only a CSV/TSV text file is generated; for all
other modes, however, a SD file is generated by default - you can change the behavior to genereate
text file using I<--output> option.

For I<3DCmpdRecords> mode, only those compounds with at least one non-zero value for Z atomic coordinates
are retrieved; however, during retrieval of compounds in I<2DCmpdRecords> mode, all Z atomic coordinates must
be zero.

=item B<-n, --numofcmpds> I<number>

Number of compouds to extract during I<randomcmpds> mode.

=item B<--outdelim> I<comma | tab | semicolon>

Delimiter for output CSV/TSV text file(s). Possible values: I<comma, tab, or semicolon>
Default value: I<comma>

=item B<--output> I<SD | text | both>

Type of output files to generate. Possible values: I<SD, text, or both>. Default value: I<SD>. For
I<alldatafields> and I<molnames> mode, this option is ingored and only a CSV/TSV text file is generated.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-q, --quote> I<yes | no>

Put quote around column values in output CSV/TSV text file(s). Possible values:
I<yes or no>. Default value: I<yes>.

=item B<--record> I<recnum | recnums | startrecnum,endrecnum>

Record number, record numbers or range of records to extract during I<recordnum>, I<recordnums>
and I<recordrange> mode. Input value format is: <num>, <num1,num2,...> and <startnum, endnum>
for I<recordnum>, I<recordnums> and I<recordrange> modes recpectively. Default value: none.

=item B<--RegexIgnoreCase> I<yes or no>

Specify whether to ingnore case during I<datafieldsbyregex> value of B<-m, --mode> option.
Possible values: I<yes or no>. Default value: I<yes>.

=item B<-r, --root> I<rootname>

New file name is generated using the root: <Root>.<Ext>. Default for new file
names: <SDFileName><mode>.<Ext>. The file type determines <Ext> value.
The sdf, csv, and tsv <Ext> values are used for SD, comma/semicolon, and tab
delimited text files respectively.This option is ignored for multiple input files.

=item B<-s, --seed> I<number>

Random number seed used for I<randomcmpds> mode. Default:123456789.

=item B<--StrDataString> I<yes | no>

Specify whether to write out structure data string to CSV/TSV text file(s). Possible values:
I<yes or no>. Default value: I<no>.

The value of B<StrDataStringDelimiter> option is used as a delimiter to join structure
data lines into a structure data string.

This option is ignored during generation of SD file(s).

=item B<--StrDataStringDelimiter> I<text>

Delimiter for joining multiple stucture data lines into a string before writing to CSV/TSV text
file(s). Possible values: I<any alphanumeric text>. Default value: I<|>.

This option is ignored during generation of SD file(s).

=item B<--StrDataStringMode> I<StrOnly | StrAndDataFields>

Specify whether to include SD data fields and values along with the structure data into structure
data string before writing it out to CSV/TSV text file(s). Possible values: I<StrOnly or StrAndDataFields>.
Default value: I<StrOnly>.

The value of B<StrDataStringDelimiter> option is used as a delimiter to join structure
data lines into a structure data string.

This option is ignored during generation of SD file(s).

=item B<--ValueComparisonMode> I<Numeric | Alphanumeric>

Specify how to compare data field values during I<datafieldsbyvalue> mode: Compare
values using either numeric or string ((eq, le, ge) comparison operators. Possible values:
I<Numeric or Alphanumeric>. Defaule value: I<Numeric>.

=item B<-v, --violations> I<number>

Number of criterion violations allowed for values specified during I<datafieldsbyvalue>
and I<datafieldsbyregex> mode. Default value: I<0>.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To retrieve all data fields from SD files and generate CSV text files, type:

    % ExtractFromSDFiles.pl -o Sample.sdf
    % ExtractFromSDFiles.pl -o *.sdf

To retrieve all data fields from SD file and generate CSV text files containing
a column with structure data as a string with | as line delimiter, type:

    % ExtractFromSDFiles.pl --StrDataString Yes -o Sample.sdf

To retrieve MOL_ID data fileld from SD file and generate CSV text files containing
a column with structure data along with all data fields as a string with | as line
delimiter, type:

    % ExtractFromSDFiles.pl -m datafields -d "Mol_ID" --StrDataString Yes
      --StrDataStringMode StrAndDataFields --StrDataStringDelimiter "|"
      --output text -o Sample.sdf

To retrieve common data fields which exists for all the compounds in
a SD file and generate a TSV text file NewSample.tsv, type:

    % ExtractFromSDFiles.pl -m commondatafields --outdelim tab -r NewSample
      --output Text -o Sample.sdf

To retrieve MolId, ExtReg, and CompoundName data field from a SD file and generate a
CSV text file NewSample.csv, type:

    % ExtractFromSDFiles.pl -m datafields -d "Mol_ID,MolWeight,
      CompoundName" -r NewSample --output Text -o Sample.sdf

To retrieve compounds from a SD which meet a specific set of criteria - MolWt <= 450,
LogP <= 5 and SumNO < 10 - from a SD file and generate a new SD file NewSample.sdf,
type:

    % ExtractFromSDFiles.pl -m datafieldsbyvalue -d "MolWt,450,le,LogP
      ,5,le,SumNO,10" -r NewSample -o Sample.sdf

To retrive compounds from a SD file with a specific set of values for MolID and
generate a new SD file NewSample.sdf, type:

    % ExtractFromSDFiles.pl -m datafieldbylist -d "Mol_ID,159,4509,4619"
      -r NewSample -o Sample.sdf

To retrive compounds from a SD file with values for MolID not on a list of specified
values and generate a new SD file NewSample.sdf, type:

    % ExtractFromSDFiles.pl -m datafieldnotbylist -d "Mol_ID,159,4509,4619"
      -r NewSample -o Sample.sdf

To retrive 10 random compounds from a SD file and generate a new SD file RandomSample.sdf, type:

    % ExtractFromSDFiles.pl -m randomcmpds -n 10 -r RandomSample
      -o Sample.sdf

To retrive compound record number 10 from a SD file and generate a new SD file NewSample.sdf, type:

    % ExtractFromSDFiles.pl -m recordnum --record 10 -r NewSample
      -o Sample.sdf

To retrive compound record numbers 10, 20 and 30  from a SD file and generate a new SD file
NewSample.sdf, type:

    % ExtractFromSDFiles.pl -m recordnums --record 10,20,30 -r NewSample
      -o Sample.sdf

To retrive compound records between 10 to 20 from  SD file and generate a new SD
file NewSample.sdf, type:

    % ExtractFromSDFiles.pl -m recordrange --record 10,20 -r NewSample
      -o Sample.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

FilterSDFiles.pl, InfoSDFiles.pl, SplitSDFiles.pl, MergeTextFilesWithSD.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
