#!/usr/bin/perl -w
#
# File: MergeTextFiles.pl
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
use FileHandle;
use FileUtil;
use TextUtil;

my($ScriptName, %Options, $StartTime, $EndTime, $TotalTime);

# Autoflush STDOUT
$| = 1;

# Starting message...
$ScriptName = basename $0;
print "\n$ScriptName:Starting...\n\n";
$StartTime = new Benchmark;

# Get the options and setup script...
SetupScriptUsage();
if ($Options{help} || @ARGV < 1) {
  die GetUsageFromPod("$FindBin::Bin/$ScriptName");
}

my(@TextFilesList);
@TextFilesList = ExpandFileNames(\@ARGV, "csv tsv");

if (@TextFilesList == 1) {
  die "Error: Specify more than one text file.\n";
}

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Setup information about input files...
my(%TextFilesInfo);
print "Checking input text files...\n";
RetrieveTextFilesInfo();
RetrieveColumnsAndKeysInfo();

# Merge files...
print "\nGenerating new text file $OptionsInfo{NewTextFile}...\n";
MergeTextFiles();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Merge all valid Text files...
sub MergeTextFiles {
  my($Index);

  open NEWTEXTFILE, ">$OptionsInfo{NewTextFile}" or die "Error: Couldn't open $OptionsInfo{NewTextFile}: $! \n";

  WriteNewTextFileColumnLabels(\*NEWTEXTFILE);

  #Open up all the files and skip coumn label line...
  @{$TextFilesInfo{FileHandle}} = ();
  for $Index (0 .. $#TextFilesList) {
    $TextFilesInfo{FileHandle}[$Index] = new FileHandle;

    open $TextFilesInfo{FileHandle}[$Index], "$TextFilesList[$Index]" or die "Error: Couldn't open $TextFilesList[$Index]: $! \n";
    GetTextLine($TextFilesInfo{FileHandle}[$Index]);
  }

  # Merge files...
  if ($OptionsInfo{Keys}) {
    MergeColumnValuesUsingKeys(\*NEWTEXTFILE);
  }
  else {
    MergeColumnValues(\*NEWTEXTFILE);
  }

  # Close all opened files...
  close NEWTEXTFILE;
  for $Index (0 .. $#TextFilesList) {
    close $TextFilesInfo{FileHandle}[$Index];
  }

}

# Merge all the column values...
sub MergeColumnValues {
  my($NewTextFileRef) = @_;
  my($Index, $Line, $InDelim, $Value, $ColNum, @LineWords, @File1LineWords, @ColValues);

  while ($Line = GetTextLine($TextFilesInfo{FileHandle}[0])) {
    $InDelim = $TextFilesInfo{InDelim}[0];
    @ColValues = ();

    #Collect column values from first file before the merge point...
    @File1LineWords = quotewords($InDelim, 0, $Line);
    for $ColNum (@{$TextFilesInfo{File1Part1ColNums}}) {
      $Value = ($ColNum < @File1LineWords) ? $File1LineWords[$ColNum] : "";
      push @ColValues, $Value;
    }

    #Collect column values from other text files...
    for $Index (1 .. $#TextFilesList) {
      $InDelim = $TextFilesInfo{InDelim}[$Index];
      if ($Line = GetTextLine($TextFilesInfo{FileHandle}[$Index])) {
	@LineWords = quotewords($InDelim, 0, $Line);
	for $ColNum (@{$TextFilesInfo{ColToMerge}[$Index]}) {
	  $Value = ($ColNum < @LineWords) ? $LineWords[$ColNum] : "";
	  push @ColValues, $Value;
	}
      }
    }

    #Collect column labels from first file after the merge point...
    for $ColNum (@{$TextFilesInfo{File1Part2ColNums}}) {
      $Value = ($ColNum < @File1LineWords) ? $File1LineWords[$ColNum] : "";
      push @ColValues, $Value;
    }

    # Write it out...
    $Line = JoinWords(\@ColValues, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print $NewTextFileRef "$Line\n";
  }

}

# Merge column values using keys...
sub MergeColumnValuesUsingKeys {
  my($NewTextFileRef) = @_;
  my($Index, $InDelim, $Line, $Value, $ColNum, $KeyColNum, $KeyColValue, @LineWords, @ColValues, @File1LineWords, @TextFilesKeysToLinesMap);

  @TextFilesKeysToLinesMap = ();

  # Retrieve text lines from all the files except for the first file...
  for $Index (1 .. $#TextFilesList) {
    %{$TextFilesKeysToLinesMap[$Index]} = ();

    $InDelim = $TextFilesInfo{InDelim}[$Index];
    $KeyColNum = $TextFilesInfo{KeysToUse}[$Index];

    while ($Line = GetTextLine($TextFilesInfo{FileHandle}[$Index])) {
      @LineWords = quotewords($InDelim, 0, $Line);
      if ($KeyColNum < @LineWords) {
	$KeyColValue = $LineWords[$KeyColNum];
	if (length($KeyColValue)) {
	  if (exists($TextFilesKeysToLinesMap[$Index]{$KeyColValue})) {
	    warn "Warning: Ignoring line, $Line, in text file $TextFilesList[$Index]: Column key value, $KeyColValue, already exists\n";
	  }
	  else {
	    @{$TextFilesKeysToLinesMap[$Index]{$KeyColValue}} = ();
	    push @{$TextFilesKeysToLinesMap[$Index]{$KeyColValue}}, @LineWords;
	  }
	}
      }
    }
  }

  while ($Line = GetTextLine($TextFilesInfo{FileHandle}[0])) {
    $InDelim = $TextFilesInfo{InDelim}[0];

    @ColValues = ();
    @File1LineWords = quotewords($InDelim, 0, $Line);

    $KeyColNum = $TextFilesInfo{KeysToUse}[0];
    $KeyColValue = $File1LineWords[$KeyColNum];

    #Collect column values from first file before the merge point...
    for $ColNum (@{$TextFilesInfo{File1Part1ColNums}}) {
      $Value = ($ColNum < @File1LineWords) ? $File1LineWords[$ColNum] : "";
      push @ColValues, $Value;
    }

    #Collect column values from other text files...
    for $Index (1 .. $#TextFilesList) {
      @LineWords = ();
      if (exists($TextFilesKeysToLinesMap[$Index]{$KeyColValue})) {
	push @LineWords, @{$TextFilesKeysToLinesMap[$Index]{$KeyColValue}};
      }
      for $ColNum (@{$TextFilesInfo{ColToMerge}[$Index]}) {
	$Value = ($ColNum < @LineWords) ? $LineWords[$ColNum] : "";
	push @ColValues, $Value;
      }
    }

    #Collect column labels from first file after the merge point...
    for $ColNum (@{$TextFilesInfo{File1Part2ColNums}}) {
      $Value = ($ColNum < @File1LineWords) ? $File1LineWords[$ColNum] : "";
      push @ColValues, $Value;
    }

    # Write it out...
    $Line = JoinWords(\@ColValues, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print $NewTextFileRef "$Line\n";
  }

}

# Write out column labels for new merged text file...
sub WriteNewTextFileColumnLabels {
  my($NewTextFileRef) = @_;
  my($Index, $Line, $ColNum, @ColLabels);

  #Write out column labels for the merged text file...
  @ColLabels = ();

  #Collect column labels from first file before the merge point...
  for $ColNum (@{$TextFilesInfo{File1Part1ColNums}}) {
    push @ColLabels, $TextFilesInfo{ColToMergeNumToLabelMap}[0]{$ColNum};
  }

  #Collect column labels from other text files...
  for $Index (1 .. $#TextFilesList) {
    for $ColNum (@{$TextFilesInfo{ColToMerge}[$Index]}) {
      push @ColLabels, $TextFilesInfo{ColToMergeNumToLabelMap}[$Index]{$ColNum};
    }
  }

  #Collect column labels from first file after the merge point...
  for $ColNum (@{$TextFilesInfo{File1Part2ColNums}}) {
    push @ColLabels, $TextFilesInfo{ColToMergeNumToLabelMap}[0]{$ColNum};
  }

  #Write it out...
  $Line = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print NEWTEXTFILE "$Line\n";
}

# Retrieve text file columns and keys information for specified options...
sub RetrieveColumnsAndKeysInfo {
  ProcessColumnsInfo();

  if ($OptionsInfo{Keys}) {
    ProcessKeysInfo();
  }

  ProcessStartColInfo();
}

# Process specified columns...
sub ProcessColumnsInfo {
  my($Index, $SpecifiedColNum, $Values, $ColIndex, $ColNum, $ColLabel, @Words);

  @{$TextFilesInfo{ColSpecified}} = ();
  @{$TextFilesInfo{ColToMerge}} = ();
  @{$TextFilesInfo{ColToMergeNumToLabelMap}} = ();

  for $Index (0 .. $#TextFilesList) {

    @{$TextFilesInfo{ColSpecified}[$Index]} = ();

    $Values = "all";
    if ($OptionsInfo{Columns}) {
      $Values = $OptionsInfo{ColValues}[$Index];
    }

    if ($Values =~ /all/i) {
      if ($OptionsInfo{Mode} =~ /^colnum$/i) {
	for $ColNum (1 .. $TextFilesInfo{ColCount}[$Index]) {
	  push @{$TextFilesInfo{ColSpecified}[$Index]}, $ColNum;
	}
      }
      else {
	push @{$TextFilesInfo{ColSpecified}[$Index]}, @{$TextFilesInfo{ColLabels}[$Index]};
      }
    }
    else {
      @Words = split ",", $Values;
      push @{$TextFilesInfo{ColSpecified}[$Index]}, @Words;
    }

    @{$TextFilesInfo{ColToMerge}[$Index]} = ();
    %{$TextFilesInfo{ColToMergeNumToLabelMap}[$Index]} = ();

    if ($OptionsInfo{Mode} =~ /^collabel$/i) {
      for $ColIndex (0 .. $#{$TextFilesInfo{ColSpecified}[$Index]}) {
	$ColLabel = $TextFilesInfo{ColSpecified}[$Index][$ColIndex];
	if (exists($TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel})) {
	  $ColNum = $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel};
	  push @{$TextFilesInfo{ColToMerge}[$Index]}, $ColNum;
	  $TextFilesInfo{ColToMergeNumToLabelMap}[$Index]{$ColNum} = $ColLabel;
	}
	else {
	  warn "Warning: Ignoring value, $ColLabel, specified using \"-c --column\" option: column name doesn't exist in  $TextFilesList[$Index]  \n";
	}
      }
    }
    else {
      for $ColIndex (0 .. $#{$TextFilesInfo{ColSpecified}[$Index]}) {
	$SpecifiedColNum = $TextFilesInfo{ColSpecified}[$Index][$ColIndex];
	if ($SpecifiedColNum > 0 && $SpecifiedColNum <= $TextFilesInfo{ColCount}[$Index]) {
	  $ColNum = $SpecifiedColNum - 1;
	  push @{$TextFilesInfo{ColToMerge}[$Index]}, $ColNum;
	  $TextFilesInfo{ColToMergeNumToLabelMap}[$Index]{$ColNum} = $TextFilesInfo{ColLabels}[$Index][$ColNum];
	}
	else {
	  warn "Warning: Ignoring value, $SpecifiedColNum, specified using \"-c --column\" option: column number doesn't exist in  $TextFilesList[$Index]  \n";
	}
      }
    }
    my (@ColToMergeSorted) = sort { $a <=> $b } @{$TextFilesInfo{ColToMerge}[$Index]};
    @{$TextFilesInfo{ColToMerge}[$Index]} = ();
    push @{$TextFilesInfo{ColToMerge}[$Index]}, @ColToMergeSorted;
  }
}

# Process specified key column values...
sub ProcessKeysInfo {
  my($Index, $Key, $ColLabel, $ColNum);

  @{$TextFilesInfo{KeysSpecified}} = ();
  @{$TextFilesInfo{KeysToUse}} = ();

  for $Index (0 .. $#TextFilesList) {
    $Key = $OptionsInfo{KeyValues}[$Index];

    $TextFilesInfo{KeysSpecified}[$Index] = $Key;
    $TextFilesInfo{KeysToUse}[$Index] = -1;

    if ($OptionsInfo{Mode} =~ /^collabel$/i) {
      $ColLabel = $Key;
      if (exists($TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel})) {
	$TextFilesInfo{KeysToUse}[$Index] =  $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel};
      }
      else {
	warn "Warning: Ignoring value, $ColLabel, specified using \"-k --keys\" option: column name doesn't exist in  $TextFilesList[$Index]  \n";
      }
    }
    else {
      $ColNum = $Key;
      if ($ColNum > 0 && $ColNum <= $TextFilesInfo{ColCount}[$Index]) {
	$TextFilesInfo{KeysToUse}[$Index] = $ColNum - 1;
      }
      else {
	warn "Warning: Ignoring value, $ColNum, specified using \"-k --keys\" option: column number doesn't exist in  $TextFilesList[$Index]  \n";
      }
    }
  }

  # Modify columns to merge list to make sure the columns identified by key are taken off the list
  # except for the first text file...
  my(@ColToMergeFiltered);

  for $Index (1 .. $#TextFilesList) {
    @ColToMergeFiltered = ();
    for $ColNum (@{$TextFilesInfo{ColToMerge}[$Index]}) {
      if ($TextFilesInfo{KeysToUse}[$Index] != $ColNum) {
	push @ColToMergeFiltered, $ColNum;
      }
    }
    @{$TextFilesInfo{ColToMerge}[$Index]} = ();
    push @{$TextFilesInfo{ColToMerge}[$Index]}, @ColToMergeFiltered;
  }
}

# Process specified start column value...
sub ProcessStartColInfo {
  my($Index, $ColIndex, $ColNum, $StartColNum, $Part1StartColNum, $Part1EndColNum, $Part2StartColNum, $Part2EndColNum, $BeforeStartColNum, $AfterStartColNum, $FirstColNum, $LastColNum, $FirstIndex, $LastIndex);

  @{$TextFilesInfo{File1Part1ColNums}} = ();
  @{$TextFilesInfo{File1Part2ColNums}} = ();

  $StartColNum = "last";
  if ($OptionsInfo{StartCol}) {
    if (length($OptionsInfo{StartCol})) {
      $StartColNum = $OptionsInfo{StartCol}
    }
  }

  if ($StartColNum !~ /^last$/i) {
    if ($OptionsInfo{Mode} =~ /^collabel$/i) {
      if (exists($TextFilesInfo{ColLabelToNumMap}[0]{$StartColNum})) {
	$StartColNum = $TextFilesInfo{ColLabelToNumMap}[0]{$StartColNum};
      }
      else {
	die "Error: Invalid value $StartColNum specified using \"-s --startcol\" option: column name doesn't exist in  $TextFilesList[0]  \n";
      }
    }
    else {
      if ($StartColNum > 0 && $StartColNum <= $TextFilesInfo{ColCount}[0]) {
	$StartColNum -= 1;
      }
      else {
	die "Error: Invalid value $StartColNum specified using \"-s --startcol\" option: column number doesn't exist in  $TextFilesList[0]  \n";
      }
    }
  }
  else {
    $StartColNum = $TextFilesInfo{ColCount}[0] - 1;
  }

  # Make sure StartColNum is present on the list of columns to merge for the first text file...
  if (!exists($TextFilesInfo{ColToMergeNumToLabelMap}[0]{$StartColNum})) {
    die "Error: Invalid value $StartColNum specified using \"-s --startcol\" option: doesn't exist in the specified lists of columns to merge for  $TextFilesList[0]  \n";
  }

  # Find out the column number before and after StartColNum in first text file...
  $BeforeStartColNum = $StartColNum;
  $AfterStartColNum = $StartColNum;

  $FirstIndex = 0; $LastIndex = $#{$TextFilesInfo{ColToMerge}[0]};

  $FirstColNum = $TextFilesInfo{ColToMerge}[0][$FirstIndex];
  $LastColNum = $TextFilesInfo{ColToMerge}[0][$LastIndex];

  for $Index (0 .. $LastIndex) {
    if ($TextFilesInfo{ColToMerge}[0][$Index] == $StartColNum) {
      $BeforeStartColNum = (($Index -1) >= $FirstIndex) ? $TextFilesInfo{ColToMerge}[0][$Index - 1] : ($FirstColNum - 1);
      $AfterStartColNum = (($Index + 1) <= $LastIndex) ? $TextFilesInfo{ColToMerge}[0][$Index + 1] : ($LastColNum + 1);
    }
  }

  if ($OptionsInfo{StartColMode} =~ /^after$/i) {
    $Part1StartColNum = $FirstColNum; $Part1EndColNum = $StartColNum;
    $Part2StartColNum = $AfterStartColNum; $Part2EndColNum = $LastColNum;
  }
  else {
    $Part1StartColNum = $FirstColNum; $Part1EndColNum = $BeforeStartColNum;
    $Part2StartColNum = $StartColNum; $Part2EndColNum = $LastColNum;
  }

  @{$TextFilesInfo{File1Part1ColNums}} = ();
  @{$TextFilesInfo{File1Part2ColNums}} = ();

  for $ColIndex (0 .. $#{$TextFilesInfo{ColToMerge}[0]}) {
    $ColNum = $TextFilesInfo{ColToMerge}[0][$ColIndex];
    if ($ColNum >= $Part1StartColNum && $ColNum <= $Part1EndColNum) {
      push @{$TextFilesInfo{File1Part1ColNums}}, $ColNum;
    }
  }

  for $ColIndex (0 .. $#{$TextFilesInfo{ColToMerge}[0]}) {
    $ColNum = $TextFilesInfo{ColToMerge}[0][$ColIndex];
    if ($ColNum >= $Part2StartColNum && $ColNum <= $Part2EndColNum) {
      push @{$TextFilesInfo{File1Part2ColNums}}, $ColNum;
    }
  }

}

# Retrieve information about input text files...
sub RetrieveTextFilesInfo {
  my($Index, $TextFile, $FileDir, $FileName, $FileExt, $InDelim, $Line, $ColNum, $ColLabel, $FileNotOkayCount, @ColLabels,);

  %TextFilesInfo = ();

  @{$TextFilesInfo{FileOkay}} = ();
  @{$TextFilesInfo{ColCount}} = ();
  @{$TextFilesInfo{ColLabels}} = ();
  @{$TextFilesInfo{ColLabelToNumMap}} = ();
  @{$TextFilesInfo{InDelim}} = ();

  $FileNotOkayCount = 0;

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];

    $TextFilesInfo{FileOkay}[$Index] = 0;
    $TextFilesInfo{ColCount}[$Index] = 0;
    $TextFilesInfo{InDelim}[$Index] = "";

    @{$TextFilesInfo{ColLabels}[$Index]} = ();
    %{$TextFilesInfo{ColLabelToNumMap}[$Index]} = ();

    if (!(-e $TextFile)) {
      warn "Warning: Ignoring file $TextFile: It doesn't exist\n";
      $FileNotOkayCount++;
      next FILELIST;
    }
    if (!CheckFileType($TextFile, "csv tsv")) {
      warn "Warning: Ignoring file $TextFile: It's not a csv or tsv file\n";
      $FileNotOkayCount++;
      next FILELIST;
    }
    ($FileDir, $FileName, $FileExt) = ParseFileName($TextFile);
    if ($FileExt =~ /^tsv$/i) {
      $InDelim = "\t";
    }
    else {
      $InDelim = "\,";
      if ($OptionsInfo{InDelim} !~ /^(comma|semicolon)$/i) {
	warn "Warning: Ignoring file $TextFile: The value specified, $OptionsInfo{InDelim}, for option \"--indelim\" is not valid for csv files\n";
	$FileNotOkayCount++;
	next FILELIST;
      }
      if ($OptionsInfo{InDelim} =~ /^semicolon$/i) {
	$InDelim = "\;";
      }
    }

    if (!open TEXTFILE, "$TextFile") {
      warn "Warning: Ignoring file $TextFile: Couldn't open it: $! \n";
      $FileNotOkayCount++;
      next FILELIST;
    }

    $Line = GetTextLine(\*TEXTFILE);
    @ColLabels = quotewords($InDelim, 0, $Line);
    close TEXTFILE;

    $TextFilesInfo{FileOkay}[$Index] = 1;
    $TextFilesInfo{InDelim}[$Index] = $InDelim;

    $TextFilesInfo{ColCount}[$Index] = @ColLabels;
    push @{$TextFilesInfo{ColLabels}[$Index]}, @ColLabels;
    for $ColNum (0 .. $#ColLabels) {
      $ColLabel = $ColLabels[$ColNum];
      $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel} = $ColNum;
    }
  }
  # Make sure all specified files are valid for merging to work properly...
  if ($FileNotOkayCount) {
    die "Error: Problems with input text file(s)...\n";
  }
}

# Process option values...
sub ProcessOptions {
  my($Index, $FileDir, $FileName, $FileExt, $NewTextFile, @ColValues, @KeyValues);

  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{Columns} = $Options{columns};
  @{$OptionsInfo{ColValues}} = ();

  if ($Options{columns}) {
    @ColValues = split ";", $Options{columns};
    if (@ColValues != @TextFilesList) {
      die "Error: Invalid number of values specified by \"-c --columns\" option: it must be equal to number of input text files.\n";
    }
    for $Index (0 .. $#ColValues) {
      if (!length($ColValues[$Index])) {
	die "Error: Invalid value specified by \"-c --columns\" option: empty values are not allowed.\n";
      }
    }
    @{$OptionsInfo{ColValues}} = @ColValues;
  }

  $OptionsInfo{Keys} = $Options{keys};
  @{$OptionsInfo{KeyValues}} = ();

  if ($Options{keys}) {
    @KeyValues = split ";", $Options{keys};
    if (@KeyValues != @TextFilesList) {
      die "Error: Invalid number of values specified by \"-k --keys\" option: it must be equal to number of input text files.\n";
    }
    for $Index (0 .. $#KeyValues) {
      if (!length($KeyValues[$Index])) {
	die "Error: Invalid value specified by \"-k --keys\" option: empty values are not allowed.\n";
      }
    }
    @{$OptionsInfo{KeyValues}} = @KeyValues;
  }

  $OptionsInfo{InDelim} = $Options{indelim};

  $OptionsInfo{StartCol} = $Options{startcol} ? $Options{startcol} : undef;
  $OptionsInfo{StartColMode} = $Options{startcolmode};

  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : undef;
  $OptionsInfo{Overwrite} = $Options{overwrite} ? $Options{overwrite} : undef;

  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /^tab$/i ) ? "\t" : (($Options{outdelim} =~ /^semicolon$/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^yes$/i) ? 1 : 0;

  if ($Options{root}) {
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($Options{root});
    if ($FileName && $FileExt) {
      $NewTextFile = $FileName;
    } else {
      $NewTextFile =  $Options{root};
    }
  } else {
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($TextFilesList[0]);
    $NewTextFile = $FileName . "1To" . @TextFilesList . "Merged";
  }
  if ($Options{outdelim} =~ /^tab$/i) {
    $NewTextFile .= ".tsv";
  } else {
    $NewTextFile .= ".csv";
  }
  if (!$Options{overwrite}) {
    if (-e $NewTextFile) {
      die "Error: The file $NewTextFile already exists.\n";
    }
  }
  if ($Options{root}) {
    for $Index (0 .. $#TextFilesList) {
      if (lc($NewTextFile) eq lc($TextFilesList[$Index])) {
	die "Error: Output filename, $NewTextFile, is similar to a input file name.\nSpecify a different name using \"-r --root\" option or use default name.\n";
      }
    }
  }

  $OptionsInfo{NewTextFile} = $NewTextFile;
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{mode} = "colnum";
  $Options{indelim} = "comma";
  $Options{outdelim} = "comma";
  $Options{quote} = "yes";
  $Options{startcolmode} = "after";

  if (!GetOptions(\%Options, "help|h", "indelim=s", "columns|c=s", "keys|k=s", "mode|m=s", "outdelim=s", "overwrite|o", "quote|q=s", "root|r=s", "startcol|s=s", "startcolmode=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{mode} !~ /^(colnum|collabel)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: colnum, or collabel\n";
  }
  if ($Options{indelim} !~ /^(comma|semicolon)$/i) {
    die "Error: The value specified, $Options{indelim}, for option \"--indelim\" is not valid. Allowed values: comma or semicolon\n";
  }
  if ($Options{outdelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{outdelim}, for option \"--outdelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if ($Options{quote} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{startcolmode} !~ /^(before|after)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"--startcolmode\" is not valid. Allowed values: before or after\n";
  }
}

__END__

=head1 NAME

MergeTextFiles.pl - Merge multiple CSV or TSV text files into a single text file

=head1 SYNOPSIS

MergeTextFiles.pl TextFiles...

MergeTextFiles.pl [B<-h, --help>] [B<--indelim> comma | semicolon] [B<-c, --columns> colnum,...;... | collabel,...;...]
[B<-k, --keys> colnum,...;... | collabel,...;...] [B<-m, --mode> colnum | collabel]
[B<-o, --overwrite>] [B<--outdelim> comma | tab | semicolon] [B<-q, --quote> yes | no]
[B<-r, --root> rootname] [B<-s, --startcol> colnum | collabel] [B<--startcolmode> before | after]
[B<-w, --workingdir> dirname] TextFiles...

=head1 DESCRIPTION

Merge multiple CSV or TSV I<TextFiles> into first I<TextFile> to generate a single
text file. Unless B<-k --keys> option is used, data rows from other I<TextFiles>
are added to first I<TextFile> in a sequential order, and the number of rows in first
I<TextFile> is used to determine how many rows of data are added from other
I<TextFiles>.

Multiple I<TextFiles> names are separated by space. The valid file extensions are I<.csv> and
I<.tsv> for comma/semicolon and tab delimited text files respectively. All other file names
are ignored. All the text files in a current directory can be specified by I<*.csv>,
I<*.tsv>, or the current directory name. The B<--indelim> option determines the
format of I<TextFiles>. Any file which doesn't correspond to the format indicated
by B<--indelim> option is ignored.

=head1 OPTIONS

=over 4

=item B<-h, --help>

Print this help message.

=item B<--indelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile(s)>. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a
delimiter.

=item B<-c, --columns> I<colnum,...;... | collabel,...;...>

This value is mode specific. It is a list of columns to merge into first
text file specified by column numbers or labels for each text file
delimited by ";". All specified text files are merged into first text file.

Default value: I<all;all;...>. By default, all columns from specified text files are
merged into first text file.

For I<colnum> mode, input value format is: I<colnum,...;colnum,...;...>. Example:

    "1,2;1,3,4;7,8,9"

For I<collabel> mode, input value format is: I<collabel,...;collabel,...;...>. Example:

    "MW,SumNO;SumNHOH,ClogP,PSA;MolName,Mol_Id,Extreg"

=item B<-k, --keys> I<colnum,...;... | collabel,...;...>

This value is mode specific. It specifies column keys to use for merging
all specified text files into first text file. The column keys are specified by
column numbers or labels for each text file delimited by ";".

By default, data rows from text files are merged into first file in the order they appear.

For I<colnum> mode, input value format is:I<colkeynum, colkeynum;...>. Example:

    "1;3;7"

For I<collabel> mode, input value format is:I<colkeylabel, colkeylabel;...>. Example:

    "Mol_Id;Mol_Id;Cmpd_Id"

=item B<-m, --mode> I<colnum | collabel>

Specify how to merge text files: using column numbers or column labels.
Possible values: I<colnum or collabel>. Default value: I<colnum>.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<--outdelim> I<comma | tab | semicolon>

Output text file delimiter. Possible values: I<comma, tab, or semicolon>
Default value: I<comma>.

=item B<-q, --quote> I<yes | no>

Put quotes around column values in output text file. Possible values: I<yes or
no>. Default value: I<yes>.

=item B<-r, --root> I<rootname>

New text file name is generated using the root: <Root>.<Ext>. Default file
name: <FirstTextFileName>1To<Count>Merged.<Ext>. The csv, and tsv
<Ext> values are used for comma/semicolon, and tab delimited text files
respectively.

=item B<-s, --startcol> I<colnum | collabel>

This value is mode specific. It specifies the column in first text file which is
used for start merging other text files.For I<colnum> mode, specify column
number and for I<collabel> mode, specify column label.

Default value: I<last>. Start merge after the last column.

=item B<--startcolmode> I<before | after>

Start the merge before or after the B<-s, --startcol> value. Possible values: I<before or after>
Default value: I<after>.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To merge Sample2.csv and Sample3.csv into Sample1.csv and generate
NewSample.csv, type:

    % MergeTextFiles.pl -r NewSample -o Sample1.csv Sample2.csv
      Sample3.csv

To merge all Sample*.tsv and generate NewSample.tsv file, type:

    % MergeTextFiles.pl -r NewSample --indelim comma --outdelim tab -o
      Sample*.csv

To merge column numbers "1,2" and "3,4,5" from Sample2.csv and Sample3.csv
into Sample1.csv starting before column number 3 in Sample1.csv and to generate
NewSample.csv without quoting column data, type:

    % MergeTextFiles.pl -s 3 --startcolmode before -r NewSample -q no
      -m colnum -c "all;1,2;3,4,5" -o Sample1.csv Sample2.csv
      Sample3.csv

To merge column "Mol_ID,Formula,MolWeight" and "Mol_ID,NAME,ChemBankID"
from Sample2.csv and Sample3.csv into Sample1.csv using "Mol_ID" as a column keys
starting after the last column and to generate NewSample.tsv, type:

    % MergeTextFiles.pl -r NewSample --outdelim tab -k "Mol_ID;Mol_ID;
      Mol_ID" -m collabel -c "all;Mol_ID,Formula,MolWeight;Mol_ID,NAME,
      ChemBankID" -o Sample1.csv Sample2.csv Sample3.csv

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

JoinTextFiles.pl, MergeTextFilesWithSD.pl, ModifyTextFilesFormat.pl, SplitTextFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
