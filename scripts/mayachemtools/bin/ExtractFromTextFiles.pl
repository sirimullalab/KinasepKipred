#!/usr/bin/perl -w
#
# File: ExtractFromTextFiles.pl
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
use FileHandle;
use Benchmark;
use FileUtil;
use TextUtil;

my($ScriptName, %Options, $StartTime, $EndTime, $TotalTime);

# Autoflush STDOUT
$| = 1;

$StartTime = new Benchmark;

# Starting message...
$ScriptName = basename $0;
print "\n$ScriptName:Starting...\n\n";

# Get the options and setup script...
SetupScriptUsage();
if ($Options{help} || @ARGV < 1) {
  die GetUsageFromPod("$FindBin::Bin/$ScriptName");
}

my(@TextFilesList);
@TextFilesList = ExpandFileNames(\@ARGV, "csv tsv");

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Collect column information for all the text files...
print "Checking input text file(s)...\n";
my(%TextFilesInfo);
RetrieveTextFilesInfo();
RetrieveColumnsAndRowsInfo();

# Generate output files...
my($FileIndex);
if (@TextFilesList > 1) {
  print "\nProcessing text files...\n";
}
for $FileIndex (0 .. $#TextFilesList) {
  if ($TextFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $TextFilesList[$FileIndex]...\n";
    ExtractFromTextFile($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Extract appropriate data from text file...
sub ExtractFromTextFile {
  my($Index) = @_;

  if ($OptionsInfo{Mode} =~ /^categories$/i) {
    ExtractCategoryData($Index);
  }
  elsif ($OptionsInfo{Mode} =~ /^rows$/i){
    ExtractRowsData($Index);
  }
  else {
    ExtractColumnData($Index);
  }
}

# Geneate category files...
sub ExtractCategoryData {
  my($Index) = @_;
  my($TextFile, $CategoryCol, $NewTextFile, $InDelim, @ColLabels);

  $TextFile = $TextFilesList[$Index];

  $NewTextFile = $TextFilesInfo{OutFile}[$Index];
  $CategoryCol = $TextFilesInfo{CategoryColNum}[$Index];
  $InDelim = $TextFilesInfo{InDelim}[$Index];
  @ColLabels = @{$TextFilesInfo{ColLabels}[$Index]};

  my($Line, @LineWords, $CategoryName, $CategoryCount, %CategoriesNameToCountMap, %CategoriesNameToLinesMap);
  # Collect category data...
  open TEXTFILE, "$TextFile" or die "Couldn't open $TextFile: $! \n";
  # Skip label line...
  $_ = <TEXTFILE>;

  %CategoriesNameToCountMap = ();
  %CategoriesNameToLinesMap = ();

  while ($Line = GetTextLine(\*TEXTFILE)) {
    @LineWords = quotewords($InDelim, 0, $Line);
    $CategoryName = ($CategoryCol <= @LineWords) ? $LineWords[$CategoryCol] : "";
    if (exists($CategoriesNameToCountMap{$CategoryName})) {
      $CategoriesNameToCountMap{$CategoryName} += 1;
      push @{$CategoriesNameToLinesMap{$CategoryName}}, $Line;
    }
    else {
      $CategoriesNameToCountMap{$CategoryName} = 1;
      @{$CategoriesNameToLinesMap{$CategoryName}} = ();
      push @{$CategoriesNameToLinesMap{$CategoryName}}, $Line;
    }
  }
  close TEXTFILE;

  # Setup file names for individual category files...
  my(%CategoriesNameToFileHandleMap, %CategoriesNameToFileNameMap, $CategoryFile, $CategoryFileHandle);

  %CategoriesNameToFileHandleMap = ();
  %CategoriesNameToFileNameMap = ();

  for $CategoryName (keys %CategoriesNameToCountMap) {
    $CategoryFile = $TextFilesInfo{CategoryOutFileRoot}[$Index] . "$CategoryName" . ".$TextFilesInfo{OutFileExt}[$Index]";;
    $CategoryFile =~ s/ //g;
    $CategoryFileHandle = new FileHandle;
    open $CategoryFileHandle, ">$CategoryFile" or die "Couldn't open $CategoryFile: $! \n";
    $CategoriesNameToFileNameMap{$CategoryName} = $CategoryFile;
    $CategoriesNameToFileHandleMap{$CategoryName} = $CategoryFileHandle;
  }

  # Write out summary file...
  print "Generating file $NewTextFile...\n";
  open NEWTEXTFILE, ">$NewTextFile" or die "Couldn't open $NewTextFile: $! \n";

  # Write out column labels...
  @LineWords = ("Category","Count");
  $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print NEWTEXTFILE "$Line\n";

  # Write out the category names and count...
  for $CategoryName (sort { lc($a) cmp lc($b) } keys %CategoriesNameToCountMap) {
    $CategoryCount = $CategoriesNameToCountMap{$CategoryName};
    @LineWords = ("$CategoryName","$CategoryCount");
    $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print NEWTEXTFILE "$Line\n";
  }
  close NEWTEXTFILE;

  # Write out a file for each category...
  my($ColLabelLine, $LineIndex);

  $ColLabelLine = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print "\nGenerating text files for each category...\n";

  for $CategoryName (sort { lc($a) cmp lc($b) } keys %CategoriesNameToCountMap) {
    print "Generating file $CategoriesNameToFileNameMap{$CategoryName}...\n";
    $CategoryFileHandle = $CategoriesNameToFileHandleMap{$CategoryName};
    print $CategoryFileHandle "$ColLabelLine\n";
    for $LineIndex (0 .. $#{$CategoriesNameToLinesMap{$CategoryName}}) {
      $Line = ${$CategoriesNameToLinesMap{$CategoryName}}[$LineIndex];
      @LineWords = quotewords($InDelim, 0, $Line);
      $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
      print $CategoryFileHandle "$Line\n";
    }
    close $CategoryFileHandle;
  }
}

# Extract data for specific columns...
sub ExtractColumnData {
  my($Index) = @_;
  my($TextFile, @ColNumsToExtract, $NewTextFile, $InDelim);

  $TextFile = $TextFilesList[$Index];
  $NewTextFile =$TextFilesInfo{OutFile}[$Index];
  $InDelim = $TextFilesInfo{InDelim}[$Index];
  @ColNumsToExtract = @{$TextFilesInfo{ColNumsToExtract}[$Index]};

  print "Generating file $NewTextFile...\n";
  open TEXTFILE, "$TextFile" or die "Couldn't open $TextFile: $! \n";
  open NEWTEXTFILE, ">$NewTextFile" or die "Couldn't open $NewTextFile: $! \n";

  $_ = <TEXTFILE>;
  # Write out column labels...
  my($Line, @LineWords, @ColLabels, $ColLabelLine, @ColValues, $ColValuesLine, $ColNum, $ColValue);
  @ColLabels = (); $ColLabelLine = "";
  for $ColNum (@ColNumsToExtract) {
    push @ColLabels, $TextFilesInfo{ColLabels}[$Index][$ColNum];
  }
  $ColLabelLine = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print NEWTEXTFILE "$ColLabelLine\n";

  while ($Line = GetTextLine(\*TEXTFILE)) {
    @LineWords = quotewords($InDelim, 0, $Line);
    @ColValues = (); $ColValuesLine = "";
    for $ColNum (@ColNumsToExtract) {
      $ColValue = "";
      if ($ColNum < @LineWords) {
	$ColValue = (defined $LineWords[$ColNum]) ? $LineWords[$ColNum] : "";
      }
      push @ColValues, $ColValue;
    }
    $ColValuesLine = JoinWords(\@ColValues, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print NEWTEXTFILE "$ColValuesLine\n";
  }
  close NEWTEXTFILE;
  close TEXTFILE;
}

# Extract data for specific rows...
sub ExtractRowsData {
  my($Index) = @_;
  my($TextFile, $NewTextFile, $InDelim, $SpecifiedRowsMode);

  $TextFile = $TextFilesList[$Index];
  $NewTextFile =$TextFilesInfo{OutFile}[$Index];
  $InDelim = $TextFilesInfo{InDelim}[$Index];

  $SpecifiedRowsMode = $OptionsInfo{SpecifiedRowsMode};

  print "Generating file $NewTextFile...\n";
  open TEXTFILE, "$TextFile" or die "Couldn't open $TextFile: $! \n";
  open NEWTEXTFILE, ">$NewTextFile" or die "Couldn't open $NewTextFile: $! \n";

  my($Line, $RowCount, @LineWords, @ColLabels);

  # Write out column labels...
  $Line = <TEXTFILE>;
  push @ColLabels, @{$TextFilesInfo{ColLabels}[$Index]};
  $Line = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print NEWTEXTFILE "$Line\n";

  if ($SpecifiedRowsMode =~ /^rowsbycolvalue$/i) {
    ExtractRowsByColValue($Index, \*TEXTFILE, \*NEWTEXTFILE);
  }
  elsif ($SpecifiedRowsMode =~ /^rowsbycolvaluelist$/i) {
    ExtractRowsByColValueList($Index, \*TEXTFILE, \*NEWTEXTFILE);
  }
  elsif ($SpecifiedRowsMode =~ /^rowsbycolvaluerange$/i) {
    ExtractRowsByColValueRange($Index, \*TEXTFILE, \*NEWTEXTFILE);
  }
  elsif ($SpecifiedRowsMode =~ /^(rowbymincolvalue|rowbymaxcolvalue)$/i) {
    ExtractRowByMinOrMaxColValue($Index, \*TEXTFILE, \*NEWTEXTFILE);
  }
  elsif ($SpecifiedRowsMode =~ /^rownums$/i) {
    ExtractRowsByRowNums($Index, \*TEXTFILE, \*NEWTEXTFILE);
  }
  elsif ($SpecifiedRowsMode =~ /^rownumrange$/i) {
    ExtractRowsByRowNumRange($Index, \*TEXTFILE, \*NEWTEXTFILE);
  }

  close NEWTEXTFILE;
  close TEXTFILE;
}

# Extract rows by column value...
sub ExtractRowsByColValue {
  my($Index, $TextFileRef, $NewTextFileRef) = @_;
  my($Line, $ColNum, $ColValue, $Criterion, $Value, $ValueIndex, $InDelim, @LineWords);

  $InDelim = $TextFilesInfo{InDelim}[$Index];

  LINE: while ($Line = GetTextLine($TextFileRef)) {
    @LineWords = quotewords($InDelim, 0, $Line);
    for ($ValueIndex = 0; $ValueIndex < @{$TextFilesInfo{RowValues}[$Index]}; $ValueIndex = $ValueIndex + 3) {
      $ColNum = $TextFilesInfo{RowValues}[$Index][$ValueIndex];
      $ColValue = $TextFilesInfo{RowValues}[$Index][$ValueIndex + 1];
      $Criterion = $TextFilesInfo{RowValues}[$Index][$ValueIndex + 2];
      if ($ColNum > $#LineWords) {
	next LINE;
      }
      $Value = $LineWords[$ColNum];
      if ($Criterion =~ /^le$/i) {
	if ($Value > $ColValue) {
	  next LINE;
	}
      }
      elsif ($Criterion =~ /^ge$/i) {
	if ($Value < $ColValue) {
	  next LINE;
	}
      }
      elsif ($Criterion =~ /^eq$/i) {
	if ($Value ne $ColValue) {
	  next LINE;
	}
      }
    }
    # Write it out...
    $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print $NewTextFileRef "$Line\n";
  }
}
# Extract rows by column value list...
sub ExtractRowsByColValueList {
  my($Index, $TextFileRef, $NewTextFileRef) = @_;
  my($Line, $ColNum, $ColValue, $ValueIndex, $Value, $InDelim, %ColValueMap, @LineWords);

  $InDelim = $TextFilesInfo{InDelim}[$Index];
  $ColNum = $TextFilesInfo{RowValues}[$Index][0];

  # Setup a col value map...
  %ColValueMap = ();
  for $ValueIndex (1 .. $#{$TextFilesInfo{RowValues}[$Index]}) {
    $Value = $TextFilesInfo{RowValues}[$Index][$ValueIndex];
    $ColValueMap{$Value} = $Value;
  }

  LINE: while ($Line = GetTextLine($TextFileRef)) {
    @LineWords = quotewords($InDelim, 0, $Line);
    if ($ColNum > $#LineWords) {
      next LINE;
    }
    $ColValue = $LineWords[$ColNum];
    if (exists $ColValueMap{$ColValue}) {
      $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
      print $NewTextFileRef "$Line\n";
    }
  }
}

# Extract row by minimum column value...
sub ExtractRowByMinOrMaxColValue {
  my($Index, $TextFileRef, $NewTextFileRef) = @_;
  my($Line, $ColNum, $ColValue, $FirstValue, $ValueLine, $InDelim, @LineWords);

  $InDelim = $TextFilesInfo{InDelim}[$Index];
  $ColNum = $TextFilesInfo{RowValues}[$Index][0];

  $ValueLine = ''; $ColValue = ''; $FirstValue = 1;
  LINE: while ($Line = GetTextLine($TextFileRef)) {
    @LineWords = quotewords($InDelim, 0, $Line);
    if ($ColNum > $#LineWords) {
      next LINE;
    }
    if ($FirstValue) {
      $FirstValue = 0;
      $ColValue = $LineWords[$ColNum];
      $ValueLine = $Line;
      next LINE;
    }
    if ($OptionsInfo{SpecifiedRowsMode} =~ /^rowbymaxcolvalue$/i) {
      if ($LineWords[$ColNum] > $ColValue) {
	$ColValue = $LineWords[$ColNum];
	$ValueLine = $Line;
      }
    }
    else {
      if ($LineWords[$ColNum] < $ColValue) {
	$ColValue = $LineWords[$ColNum];
	$ValueLine = $Line;
      }
    }
  }
  if ($ValueLine) {
    @LineWords = quotewords($InDelim, 0, $ValueLine);
    $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print $NewTextFileRef "$Line\n";
  }
}

# Extract rows by column value range...
sub ExtractRowsByColValueRange {
  my($Index, $TextFileRef, $NewTextFileRef) = @_;
  my($Line, $ColNum, $ColValue, $MinValue, $MaxValue, $InDelim, @LineWords);

  $InDelim = $TextFilesInfo{InDelim}[$Index];
  $ColNum = $TextFilesInfo{RowValues}[$Index][0];
  $MinValue = $TextFilesInfo{RowValues}[$Index][1];
  $MaxValue = $TextFilesInfo{RowValues}[$Index][2];

  LINE: while ($Line = GetTextLine($TextFileRef)) {
    @LineWords = quotewords($InDelim, 0, $Line);
    if ($ColNum > $#LineWords) {
      next LINE;
    }
    $ColValue = $LineWords[$ColNum];
    if ($ColValue >= $MinValue && $ColValue <= $MaxValue) {
      $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
      print $NewTextFileRef "$Line\n";
    }
  }
}

# Extract rows by row number range...
sub ExtractRowsByRowNumRange {
  my($Index, $TextFileRef, $NewTextFileRef) = @_;

  my($Line, $MinRowNum, $MaxRowNum, $RowCount, $InDelim, @LineWords);
  $InDelim = $TextFilesInfo{InDelim}[$Index];
  $MinRowNum = $TextFilesInfo{RowValues}[$Index][0];
  $MaxRowNum = $TextFilesInfo{RowValues}[$Index][1];

  $RowCount = 1;
  LINE: while ($Line = GetTextLine($TextFileRef)) {
    $RowCount++;
    if ($RowCount >= $MinRowNum && $RowCount <= $MaxRowNum) {
      @LineWords = quotewords($InDelim, 0, $Line);
      $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
      print $NewTextFileRef "$Line\n";
    }
    elsif ($RowCount > $MaxRowNum) {
      last LINE;
    }
  }
}

# Extract rows by row numbers...
sub ExtractRowsByRowNums {
  my($Index, $TextFileRef, $NewTextFileRef) = @_;
  my($Line, $RowNum, $MaxRowNum, $RowCount, $InDelim, %RowNumMap, @LineWords);

  $InDelim = $TextFilesInfo{InDelim}[$Index];

  # Setup a row nums map...
  %RowNumMap = ();
  $MaxRowNum = $TextFilesInfo{RowValues}[$Index][0];
  for $RowNum (@{$TextFilesInfo{RowValues}[$Index]}) {
    if ($RowNum > $MaxRowNum) {
      $MaxRowNum = $RowNum;
    }
    $RowNumMap{$RowNum} = $RowNum;
  }

  $RowCount = 1;
  LINE: while ($Line = GetTextLine($TextFileRef)) {
    $RowCount++;
    if (exists $RowNumMap{$RowCount}) {
      @LineWords = quotewords($InDelim, 0, $Line);
      $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
      print $NewTextFileRef "$Line\n";
    }
    elsif ($RowCount > $MaxRowNum) {
      last LINE;
    }
  }
}

# Retrieve text file columns and rows information for specified options...
sub RetrieveColumnsAndRowsInfo {
  ProcessColumnsInfo();
  ProcessRowsInfo();
}

# Make sure the specified columns exists in text files...
sub ProcessColumnsInfo {
  my($Index, $SpecifiedCategoryCol, $TextFile, @ColNumsToExtract);

  @{$TextFilesInfo{CategoryColNum}} = ();
  @{$TextFilesInfo{ColNumsToExtract}} = ();

  $SpecifiedCategoryCol = $OptionsInfo{SpecifiedCategoryCol};

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];

    $TextFilesInfo{CategoryColNum}[$Index] = 0;
    @{$TextFilesInfo{ColNumsToExtract}[$Index]} = ();

    if ($TextFilesInfo{FileOkay}[$Index]) {
      if ($OptionsInfo{Mode} =~ /^categories$/i) {
	my($CategoryColNum, $CategoryColValid);

	$CategoryColNum = 0;
	$CategoryColValid = 1;
	if ($SpecifiedCategoryCol) {
	  if ($OptionsInfo{ColMode} =~ /^colnum$/i) {
	    if ($SpecifiedCategoryCol <= $TextFilesInfo{ColCount}[$Index]) {
	      $CategoryColNum = $SpecifiedCategoryCol - 1;
	    }
	    else {
	      $CategoryColValid = 0;
	    }
	  }
	  else {
	    if (exists($TextFilesInfo{ColLabelToNumMap}[$Index]{$SpecifiedCategoryCol})) {
	      $CategoryColNum =  $TextFilesInfo{ColLabelToNumMap}[$Index]{$SpecifiedCategoryCol};
	    }
	    else {
	      $CategoryColValid = 0;
	    }
	  }
	}
	if ($CategoryColValid) {
	  $TextFilesInfo{CategoryColNum}[$Index] = $CategoryColNum;
	}
	else {
	  warn "Warning: Ignoring file $TextFile: Category column specified, $SpecifiedCategoryCol, using \"--categorycol\" option doesn't exist\n";
	  $TextFilesInfo{FileOkay}[$Index] = 0;
	}
      }
      elsif ($OptionsInfo{Mode} =~ /^columns$/i) {
	my($SpecifiedColNum, $ColNum);

	$ColNum = 0;
	@ColNumsToExtract = ();

	if (@{$OptionsInfo{SpecifiedColumns}}) {
	  if ($OptionsInfo{ColMode} =~ /^colnum$/i) {
	    for $SpecifiedColNum (@{$OptionsInfo{SpecifiedColumns}}) {
	      if ($SpecifiedColNum >=1 && $SpecifiedColNum <= $TextFilesInfo{ColCount}[$Index]) {
		$ColNum = $SpecifiedColNum - 1;
		push @ColNumsToExtract, $ColNum;
	      }
	    }
	  }
	  else {
	    my($ColLabel);
	    for $ColLabel (@{$OptionsInfo{SpecifiedColumns}}) {
	      if (exists($TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel})) {
		push @ColNumsToExtract, $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel};
	      }
	    }
	  }
	}
	else {
	  push @ColNumsToExtract, $ColNum;
	}
	if (@ColNumsToExtract) {
	  push @{$TextFilesInfo{ColNumsToExtract}[$Index]}, @ColNumsToExtract;
	}
	else {
	  warn "Warning: Ignoring file $TextFile: None of the columns specified, @{$OptionsInfo{SpecifiedColumns}}, using \"--columns\" option exist\n";
	  $TextFilesInfo{FileOkay}[$Index] = 0;
	}
      }
    }
  }
}

# Process specified rows info...
sub ProcessRowsInfo {
  my($Index, $TextFile, $ColID, $ColIDOkay, $Value, $Criterion, $ColNum, @RowValues);

  @{$TextFilesInfo{RowValues}} = ();

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];
    @{$TextFilesInfo{RowValues}[$Index]} = ();

    if ($OptionsInfo{Mode} !~ /^rows$/i) {
      next FILELIST;
    }
    if (!$TextFilesInfo{FileOkay}[$Index]) {
      next FILELIST;
    }

    @RowValues = ();

    if ($OptionsInfo{RowsMode} =~ /^rowsbycolvalue$/i) {
      my($ValueIndex);
      for ($ValueIndex = 0; $ValueIndex < @{$OptionsInfo{SpecifiedRowValues}}; $ValueIndex = $ValueIndex + 3) {
	$ColID = $OptionsInfo{SpecifiedRowValues}[$ValueIndex];
	$Value = $OptionsInfo{SpecifiedRowValues}[$ValueIndex + 1];
	$Criterion = $OptionsInfo{SpecifiedRowValues}[$ValueIndex + 2];

	$ColIDOkay = 0;
	if ($OptionsInfo{ColMode} =~ /^collabel$/i) {
	  if (exists $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColID}) {
	    $ColIDOkay = 1;
	    $ColNum = $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColID};
	  }
	}
	else {
	  if ($ColID >=1 && $ColID <= $TextFilesInfo{ColCount}[$Index]) {
	    $ColNum = $ColID - 1;
	    $ColIDOkay = 1;
	  }
	}
	if ($ColIDOkay) {
	  push @RowValues, ($ColNum, $Value, $Criterion);
	}
      }
    }
    elsif ($OptionsInfo{RowsMode} =~ /^(rowsbycolvaluelist|rowsbycolvaluerange|rowbymincolvalue|rowbymaxcolvalue)$/i) {
      # Process coulumn id...
      $ColID = $OptionsInfo{SpecifiedRowValues}[0];
      $ColIDOkay = 0;

      if ($OptionsInfo{ColMode} =~ /^collabel$/i) {
	if (exists $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColID}) {
	  $ColIDOkay = 1;
	  $ColNum = $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColID};
	}
      }
      else {
	if ($ColID >=1 && $ColID <= $TextFilesInfo{ColCount}[$Index]) {
	  $ColIDOkay = 1;
	  $ColNum = $ColID - 1;
	}
      }
      if ($ColIDOkay) {
	push @RowValues, $ColNum;
	# Get rest of the specified values...
	if (@{$OptionsInfo{SpecifiedRowValues}} > 1) {
	  for $Index (1 .. $#{$OptionsInfo{SpecifiedRowValues}}) {
	    push @RowValues, $OptionsInfo{SpecifiedRowValues}[$Index];
	  }
	}
      }
    }
    elsif ($OptionsInfo{RowsMode} =~ /^(rownums|rownumrange)$/i) {
      push @RowValues, @{$OptionsInfo{SpecifiedRowValues}};
    }

    if (@RowValues) {
      push @{$TextFilesInfo{RowValues}[$Index]}, @RowValues;
    }
    else {
      warn "Warning: Ignoring file $TextFile: Column specified, $ColID, using \"--rows\" option doesn't exist\n";
      $TextFilesInfo{FileOkay}[$Index] = 0;
    }
  }
}

# Retrieve information about input text files...
sub RetrieveTextFilesInfo {
  my($Index, $TextFile, $FileDir, $FileName, $FileExt, $InDelim, $Line, @ColLabels, $OutFileRoot, $CategoryOutFileRoot, $OutFile, $ColNum, $ColLabel);

  %TextFilesInfo = ();

  @{$TextFilesInfo{FileOkay}} = ();
  @{$TextFilesInfo{ColCount}} = ();
  @{$TextFilesInfo{ColLabels}} = ();
  @{$TextFilesInfo{ColLabelToNumMap}} = ();
  @{$TextFilesInfo{InDelim}} = ();
  @{$TextFilesInfo{OutFile}} = ();
  @{$TextFilesInfo{OutFileExt}} = ();
  @{$TextFilesInfo{CategoryOutFileRoot}} = ();

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];

    $TextFilesInfo{FileOkay}[$Index] = 0;
    $TextFilesInfo{ColCount}[$Index] = 0;
    $TextFilesInfo{InDelim}[$Index] = "";
    $TextFilesInfo{OutFile}[$Index] = "";
    $TextFilesInfo{OutFileExt}[$Index] = "";
    $TextFilesInfo{CategoryOutFileRoot}[$Index] = "";

    @{$TextFilesInfo{ColLabels}[$Index]} = ();
    %{$TextFilesInfo{ColLabelToNumMap}[$Index]} = ();

    if (!(-e $TextFile)) {
      warn "Warning: Ignoring file $TextFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($TextFile, "csv tsv")) {
      warn "Warning: Ignoring file $TextFile: It's not a csv or tsv file\n";
      next FILELIST;
    }

    ($FileDir, $FileName, $FileExt) = ParseFileName($TextFile);
    if ($FileExt =~ /^tsv$/i) {
      $InDelim = "\t";
    }
    else {
      $InDelim = "\,";
      if (!($OptionsInfo{InDelim} =~ /^(comma|semicolon)$/i)) {
	warn "Warning: Ignoring file $TextFile: The value specified, $OptionsInfo{InDelim}, for option \"--indelim\" is not valid for csv files\n";
	next FILELIST;
      }
      if ($OptionsInfo{InDelim} =~ /^semicolon$/i) {
	$InDelim = "\;";
      }
    }

    if (!open TEXTFILE, "$TextFile") {
      warn "Warning: Ignoring file $TextFile: Couldn't open it: $! \n";
      next FILELIST;
    }

    $Line = GetTextLine(\*TEXTFILE);
    @ColLabels = quotewords($InDelim, 0, $Line);
    close TEXTFILE;

    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($TextFile);
    $FileExt = "csv";
    if ($Options{outdelim} =~ /^tab$/i) {
      $FileExt = "tsv";
    }

    if ($OptionsInfo{OutFileRoot} && (@TextFilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($OptionsInfo{OutFileRoot});
      if ($RootFileName && $RootFileExt) {
	$FileName = $RootFileName;
      }
      else {
	$FileName = $OptionsInfo{OutFileRoot};
      }
      $OutFileRoot .= $FileName;
    }
    else {
      $OutFileRoot = $FileName;
      $OutFileRoot .= ($OptionsInfo{Mode} =~ /^categories$/i) ? "CategoriesSummary" : (($OptionsInfo{Mode} =~ /^rows$/i) ? "ExtractedRows" : "ExtractedColumns");
    }
    $CategoryOutFileRoot = "$FileName" . "Category";

    $OutFile = $OutFileRoot . ".$FileExt";
    if (lc($OutFile) eq lc($TextFile)) {
      warn "Warning: Ignoring file $TextFile:Output file name, $OutFile, is same as input text file name, $TextFile\n";
      next FILELIST;
    }

    if (!$OptionsInfo{Overwrite}) {
      if (-e $OutFile) {
	warn "Warning: Ignoring file $TextFile: The file $OutFile already exists\n";
	next FILELIST;
      }
    }

    $TextFilesInfo{FileOkay}[$Index] = 1;
    $TextFilesInfo{InDelim}[$Index] = $InDelim;
    $TextFilesInfo{CategoryOutFileRoot}[$Index] = $CategoryOutFileRoot;
    $TextFilesInfo{OutFile}[$Index] = "$OutFile";
    $TextFilesInfo{OutFileExt}[$Index] = "$FileExt";

    $TextFilesInfo{ColCount}[$Index] = @ColLabels;
    push @{$TextFilesInfo{ColLabels}[$Index]}, @ColLabels;

    for $ColNum (0 .. $#ColLabels) {
      $ColLabel = $ColLabels[$ColNum];
      $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel} = $ColNum;
    }
  }
}

# Process option values...
sub ProcessOptions {
  my(@SpecifiedColumns, @SpecifiedRowValues);

  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{ColMode} = $Options{colmode};

  $OptionsInfo{CategoryCol} = defined $Options{categorycol} ? $Options{categorycol} : undef;
  $OptionsInfo{SpecifiedCategoryCol} = "";

  if (defined $Options{categorycol}) {
    my(@SpecifiedValues) = split ",", $Options{categorycol};
    if (@SpecifiedValues != 1) {
      die "Error: Invalid number of values, ",scalar(@SpecifiedValues), " using \"--categorycol\" option: Only one value is allowed.\n";
    }
    $OptionsInfo{SpecifiedCategoryCol} = $SpecifiedValues[0];
    if ($Options{colmode} =~ /^colnum$/i) {
      if (!IsPositiveInteger($OptionsInfo{SpecifiedCategoryCol})) {
	die "Error: Category column value, $OptionsInfo{SpecifiedCategoryCol}, specified using \"--categorycol\" is not valid. Allowed integer values: > 0.\n";
      }
    }
  }

  $OptionsInfo{Columns} = defined $Options{columns} ? $Options{columns} : undef;
  @{$OptionsInfo{SpecifiedColumns}} = ();
  @SpecifiedColumns = ();

  if (defined $Options{columns}) {
    my(@SpecifiedValues) = split ",", $Options{columns};
    if ($Options{colmode} =~ /^colnum$/i) {
      my($ColValue);
      for $ColValue (@SpecifiedValues) {
	if (!IsPositiveInteger($ColValue)) {
	  die "Error: Column value, $ColValue, specified using \"--columns\" is not valid: Allowed integer values: > 0.\n";
	}
      }
    }
    push @SpecifiedColumns, @SpecifiedValues;
  }
  @{$OptionsInfo{SpecifiedColumns}} = @SpecifiedColumns;

  $OptionsInfo{InDelim} = $Options{indelim};

  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /^tab$/i ) ? "\t" : (($Options{outdelim} =~ /^semicolon$/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^yes$/i) ? 1 : 0;
  $OptionsInfo{Overwrite} = defined $Options{overwrite} ? $Options{overwrite} : undef;

  $OptionsInfo{OutFileRoot} = defined $Options{root} ? $Options{root} : undef;

  # Process any specified rows values...
  @SpecifiedRowValues = ();
  @{$OptionsInfo{SpecifiedRowValues}} = ();

  $OptionsInfo{RowsMode} = $Options{rowsmode};
  $OptionsInfo{Rows} = defined $Options{rows} ? $Options{rows} : undef;

  $OptionsInfo{SpecifiedRowsMode} = $Options{rowsmode};

  if (defined $Options{rows}) {
    (@SpecifiedRowValues) = split ",", $Options{rows};
  }
  else {
    if ($Options{rowsmode} !~ /^rownums$/i) {
      die "Error: Specify value for \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode}.\n";
    }
    push @SpecifiedRowValues, "1";
  }
  @{$OptionsInfo{SpecifiedRowValues}} = @SpecifiedRowValues;

  my($SpecifiedColID, $SpecifiedRowID);
  # Make sure specified values are okay...
  if ($Options{rowsmode} =~ /^rowsbycolvalue$/i) {
    if (@SpecifiedRowValues % 3) {
      die "Error: Invalid number of values, ", scalar(@SpecifiedRowValues) , ", specified by \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode}.\nIt must contain triplets.\n";
    }
    # Triplet format: colid,value,criteria. Criterion: le,ge,eq
    my($Index, $ColID, $Criterion, $Value);
    for ($Index = 0; $Index < @SpecifiedRowValues; $Index = $Index + 3) {
      $ColID = $SpecifiedRowValues[$Index];
      $Value = $SpecifiedRowValues[$Index + 1];
      $Criterion = $SpecifiedRowValues[$Index + 2];
      if ($Options{colmode} =~ /^colnum$/i) {
	if (!IsPositiveInteger($ColID)) {
	  die "Error: Invalid column id, $ColID, specified in triplet, \"$ColID,$Criterion,$Value\", using \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode} is not valid. Allowed integer values: > 0.\n";
	}
      }
      if ($Criterion !~ /^(eq|le|ge)$/i) {
	die "Error: Invalid criterion value, $Criterion, specified in triplet, \"$ColID,$Criterion,$Value\", using \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode} is not valid. Allowed values: le, ge, or eq.\n";
      }
    }
  }
  elsif ($Options{rowsmode} =~ /^rowsbycolvaluelist$/i) {
    ($SpecifiedColID) = $SpecifiedRowValues[0];
    if ($Options{colmode} =~ /^colnum$/i) {
      if (!IsPositiveInteger($SpecifiedColID)) {
	die "Error: Rows value, $SpecifiedColID, specified using \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode} is not valid. Allowed integer values: > 0.\n";
      }
    }
    if (@SpecifiedRowValues == 1) {
      die "Error: Invalid number of values, ", scalar(@SpecifiedRowValues) , ", specified by \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode}.\nIt must contain more than one value\n";
    }
  }
  elsif ($Options{rowsmode} =~ /^rowsbycolvaluerange$/i) {
    if (@SpecifiedRowValues != 3) {
      die "Error: Invalid number of values, ", scalar(@SpecifiedRowValues) , ", specified by \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode}.\nIt must contain three values\n";
    }
    ($SpecifiedColID) = $SpecifiedRowValues[0];
    if ($Options{colmode} =~ /^colnum$/i) {
      if (!IsPositiveInteger($SpecifiedColID)) {
	die "Error: Rows value, $SpecifiedColID, specified using \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode} is not valid. Allowed integer values: > 0.\n";
      }
    }
    if ($SpecifiedRowValues[1] >= $SpecifiedRowValues[2]) {
      die "Error: Invalid value triplet - ", JoinWords(\@SpecifiedRowValues, ',', 0) , " - specified by \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode}.\nAllowed values: second value < third value\n";
    }
  }
  elsif ($Options{rowsmode} =~ /^(rowbymincolvalue|rowbymaxcolvalue)$/i) {
    if (@SpecifiedRowValues != 1) {
      die "Error: Invalid number of values, ", scalar(@SpecifiedRowValues) , ", specified by \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode}.\nOnly one value is allowed.\n";
    }
    ($SpecifiedColID) = $SpecifiedRowValues[0];
    if ($Options{colmode} =~ /^colnum$/i) {
      if (!IsPositiveInteger($SpecifiedColID)) {
	die "Error: Rows value, $SpecifiedColID, specified using \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode} is not valid. Allowed integer values: > 0.\n";
      }
    }
  }
  elsif ($Options{rowsmode} =~ /^rownums$/i) {
    for $SpecifiedRowID (@SpecifiedRowValues) {
      if (!IsPositiveInteger($SpecifiedRowID)) {
	die "Error: Rows value, $SpecifiedRowID, specified using \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode} is not valid. Allowed integer values: > 0.\n";
      }
    }
  }
  elsif ($Options{rowsmode} =~ /^rownumrange$/i) {
    if (@SpecifiedRowValues != 2) {
      die "Error: Invalid number of values, ", scalar(@SpecifiedRowValues) , ", specified by \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode}.\nIt must contain only two values.\n";
    }
    for $SpecifiedRowID (@SpecifiedRowValues) {
      if (!IsPositiveInteger($SpecifiedRowID)) {
	die "Error: Rows value, $SpecifiedRowID, specified using \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode} is not valid. Allowed integer values: > 0.\n";
      }
    }
    if ($SpecifiedRowValues[0] >= $SpecifiedRowValues[1]) {
      die "Error: Invalid value pair -  ", JoinWords(\@SpecifiedRowValues, ',', 0) , " - specified by \"--rows\" option with \"--rowsmode\" value of $Options{rowsmode}.\nAllowed values: First value < second value\n";
    }
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Setup default and retrieve all the options...
  %Options = ();
  $Options{colmode} = "colnum";
  $Options{indelim} = "comma";
  $Options{mode} = "columns";
  $Options{outdelim} = "comma";
  $Options{quote} = "yes";
  $Options{rowsmode} = "rownums";

  if (!GetOptions(\%Options, "categorycol=s", "columns=s", "colmode|c=s", "help|h", "indelim=s", "mode|m=s", "outdelim=s", "overwrite|o", "quote|q=s", "root|r=s", "rows=s", "rowsmode=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} || die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{mode} !~ /^(columns|rows|categories)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: columns, rows or categories \n";
  }
  if ($Options{colmode} !~ /^(colnum|collabel)$/i) {
    die "Error: The value specified, $Options{colmode}, for option \"--colmode\" is not valid. Allowed values: colnum or collabel \n";
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
  if ($Options{rowsmode} !~ /^(rowsbycolvalue|rowsbycolvaluelist|rowsbycolvaluerange|rowbymincolvalue|rowbymaxcolvalue|rownums|rownumrange)$/i) {
    die "Error: The value specified, $Options{rowsmode}, for option \"--rowsmode\" is not valid. Allowed values: rowsbycolvalue, rowsbycolvaluelist, rowsbycolvaluerange, rowbymincolvalue, rowbymaxcolvalue, rownum, rownumrange\n";
  }
}
__END__


=head1 NAME

ExtractFromTextFiles.pl - Extract specific data from TextFile(s)

=head1 SYNOPSIS

ExtractFromTextFiles.pl TextFile(s)...

ExtractFromTextFiles.pl [B<-c, --colmode> colnum | collabel] [B<--categorycol > number | string]
[B<--columns> "colnum,[colnum]..." | "collabel,[collabel]..."] [B<-h, --help>]
[B<--indelim> I<comma | semicolon>] [B<-m, --mode > I<columns | rows | categories>]
[B<-o, --overwrite>] [B<--outdelim> I<comma | tab | semicolon>] [B<-q, --quote> I<yes | no>]
[B<--rows> "colid,value,criteria..." | "colid,value..." | "colid,mincolvalue,maxcolvalue" | "rownum,rownum,..." | colid | "minrownum,maxrownum"]
[ B<--rowsmode> rowsbycolvalue | rowsbycolvaluelist | rowsbycolvaluerange | rowbymincolvalue | rowbymaxcolvalue | rownums | rownumrange]
[B<-r, --root> I<rootname>] [B<-w, --workingdir> I<dirname>] TextFile(s)...

=head1 DESCRIPTION

Extract column(s)/row(s) data from I<TextFile(s)> identified by column numbers or labels. Or categorize
data using a specified column category. During categorization, a summary text file is
generated containing category name and count; an additional text file, containing data for
for each category, is also generated. The file names are separated by space. The
valid file extensions are I<.csv> and I<.tsv> for comma/semicolon and tab delimited
text files respectively. All other file names are ignored. All the text files in a
current directory can be specified by I<*.csv>, I<*.tsv>, or the current directory
name. The B<--indelim> option determines the format of I<TextFile(s)>. Any file
which doesn't correspond to the format indicated by B<--indelim> option is ignored.

=head1 OPTIONS

=over 4

=item B<-c, --colmode> I<colnum | collabel>

Specify how columns are identified in I<TextFile(s)>: using column number or column
label. Possible values: I<colnum or collabel>. Default value: I<colnum>.

=item B<--categorycol > I<number | string>

Column used to categorize data. Default value: First column.

For I<colnum> value of B<-c, --colmode> option, input value is a column number.
Example: I<1>.

For I<collabel> value of B<-c, --colmode> option, input value is a column label.
Example: I<Mol_ID>.

=item B<--columns> I<"colnum,[colnum]..." | "collabel,[collabel]...">

List of comma delimited columns to extract. Default value: First column.

For I<colnum> value of B<-c, --colmode> option, input values format is:
I<colnum,colnum,...>. Example: I<1,3,5>

For I<collabel> value of B<-c, --colmode> option, input values format is:
I<collabel,collabel,..>. Example: I<Mol_ID,MolWeight>

=item B<-h, --help>

Print this help message.

=item B<--indelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile(s)>. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a
delimiter.

=item B<-m, --mode > I<columns | rows | categories>

Specify what to extract from I<TextFile(s)>. Possible values: I<columns, rows,
or categories>. Default value: I<columns>.

For I<columns> mode, data for appropriate columns specified by B<--columns> option
is extracted from I<TextFile(s)> and placed into new text files.

For I<rows> mode, appropriate rows specified in conjuction with B<--rowsmode> and
B<rows> options are extracted from I<TextFile(s)> and placed into new text files.

For I<categories> mode, coulmn specified by B<--categorycol> is
used to categorize data, and a summary text file is generated
containing category name and count;  an additional text file, containing data for
for each category, is also generated.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<--outdelim> I<comma | tab | semicolon>.

Output text file delimiter. Possible values: I<comma, tab, or semicolon>.
Default value: I<comma>

=item B<-q, --quote> I<yes | no>

Put quotes around column values in output text file. Possible values: I<yes or
no>. Default value: I<yes>.

=item B<-r, --root> I<rootname>

New file name is generated using the root: <Root>.<Ext>. Default for new file
names: <TextFile>CategoriesSummary.<Ext>, <TextFile>ExtractedColumns.<Ext>, and
<TextFile>ExtractedRows.<Ext> for I<categories>, I<columns>, and I<rows> mode
respectively. And <TextFile>Category<CategoryName>.<Ext>
for each category retrieved from each text file. The output file type determines <Ext>
value: csv and tsv for CSV, and TSV files respectively.

This option is ignored for multiple input files.

=item B<--rows> I<"colid,value,criteria..." | "colid,value..." | "colid,mincolvalue,maxcolvalue" | "rownum,rownum,..." | colid | "minrownum,maxrownum">

This value is B<--rowsmode> specific. In general, it's a list of comma separated column ids and
associated mode specific value. Based on Column ids specification, column label or number, is
controlled by B<-c, --colmode> option.

First line containing column labels is always written out. And value comparisons assume
numerical column data.

For I<rowsbycolvalue> mode, input value format contains these triplets:
I<colid,value, criteria...>. Possible values for criteria: I<le, ge or eq>.
Examples:

    MolWt,450,le
    MolWt,450,le,LogP,5,le,SumNumNO,10,le,SumNHOH,5,le

For I<rowsbycolvaluelist> mode, input value format is: I<colid,value...>. Examples:

    Mol_ID,20
    Mol_ID,20,1002,1115

For I<rowsbycolvaluerange> mode, input value format is: I<colid,mincolvalue,maxcolvalue>. Examples:

    MolWt,100,450

For I<rowbymincolvalue, rowbymaxcolvalue> modes, input value format is: I<colid>.

For I<rownum> mode, input value format is: I<rownum>. Default value: I<2>.

For I<rownumrange> mode, input value format is: I<minrownum, maxrownum>. Examples:

    10,40

=item B<--rowsmode> I<rowsbycolvalue | rowsbycolvaluelist | rowsbycolvaluerange | rowbymincolvalue | rowbymaxcolvalue | rownums | rownumrange>

Specify how to extract rows from I<TextFile(s)>. Possible values: I<rowsbycolvalue, rowsbycolvaluelist, rowsbycolvaluerange,
rowbymincolvalue, rowbymaxcolvalue, rownum, rownumrange>. Default value: I<rownum>.

Use B<--rows> option to list rows criterion used for extraction of rows from
I<TextFile(s)>.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To extract first column from a text file and generate a new CSV text file NewSample1.csv,
type:

    % ExtractFromTextFiles.pl -r NewSample1 -o Sample1.csv

To extract columns Mol_ID, MolWeight, and NAME from Sample1.csv and generate a new
textfile NewSample1.tsv with no quotes, type:

    % ExtractFromTextFiles.pl -m columns -c collabel --columns "Mol_ID,
      MolWeight,NAME" --outdelim tab --quote no -r NewSample1
      -o Sample1.csv

To extract rows containing values for MolWeight column of less than 450 from
Sample1.csv and generate a new textfile NewSample1.csv, type:

    % ExtractFromTextFiles.pl -m rows --rowsmode rowsbycolvalue
      -c collabel --rows MolWeight,450,le -r NewSample1
      -o Sample1.csv

To extract rows containing values for MolWeight column between 400 and 500 from
Sample1.csv and generate a new textfile NewSample1.csv, type:

    % ExtractFromTextFiles.pl -m rows --rowsmode rowsbycolvaluerange
      -c collabel --rows MolWeight,450,500 -r NewSample1
      -o Sample1.csv

To extract a row containing minimum value for column MolWeight from Sample1.csv and generate
a new textfile NewSample1.csv, type:

    % ExtractFromTextFiles.pl -m rows --rowsmode rowbymincolvalue
      -c collabel --rows MolWeight -r NewSample1
      -o Sample1.csv

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
