#!/usr/bin/perl -w
#
# File: AnalyzeTextFilesData.pl
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
use StatisticsUtil;

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

my(@TextFilesList);
@TextFilesList = ExpandFileNames(\@ARGV, "csv tsv");

print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Collect column information for all the text files...
print "Checking input text file(s)...\n";
my(%TextFilesInfo);
RetrieveTextFilesInfo();
ProcessColumnsInfo();

# Generate output files...
my($FileIndex);
if (@TextFilesList > 1) {
  print "\nProcessing text files...\n";
}
for $FileIndex (0 .. $#TextFilesList) {
  if ($TextFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $TextFilesList[$FileIndex]...\n";
    AnalyzeTextFile($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Analyze data...
sub AnalyzeTextFile {
  my($Index) = @_;
  my($TextFile, $Line, $InDelim, $ColNum, $Value, @LineWords, @ColNumsToAnalyze, %ColValuesToAnalyzeMap);

  $TextFile = $TextFilesList[$Index];
  $InDelim = $TextFilesInfo{InDelim}[$Index];
  @ColNumsToAnalyze = @{$TextFilesInfo{UniqueColNumsToAnalyze}[$Index]};
  %ColValuesToAnalyzeMap = ();
  for $ColNum (@ColNumsToAnalyze) {
    @{$ColValuesToAnalyzeMap{$ColNum}} = ();
  }

  my($LineCount, $InvalidLineCount, @InvalidColLabels);

  open TEXTFILE, "$TextFile" or die "Error: Can't open $TextFile: $! \n";
  # Skip over column labels line in text file and collect appropriate column data
  # for analysis...
  $Line = GetTextLine(\*TEXTFILE);
  $LineCount = 1;
  $InvalidLineCount = 0;
  while ($Line = GetTextLine(\*TEXTFILE)) {
    $LineCount++;
    @LineWords = quotewords($InDelim, 0, $Line);
    @InvalidColLabels = ();
    COLNUM: for $ColNum (@ColNumsToAnalyze) {
      $Value = $LineWords[$ColNum];
      if ($OptionsInfo{CheckData}) {
	if (!IsNumerical($Value)) {
	  push @InvalidColLabels, $TextFilesInfo{ColLabels}[$Index][$ColNum];
	  next COLNUM;
	}
      }
      push @{$ColValuesToAnalyzeMap{$ColNum}}, $Value;
    }
    if (@InvalidColLabels) {
      $InvalidLineCount++;
      if ($OptionsInfo{DetailLevel} >=4 ) {
	print "Line number $LineCount contains ", scalar(@InvalidColLabels)," non-numerical or empty value(s) for column(s) - ", JoinWords(\@InvalidColLabels, ", ", 0)," - to be analyzed: $Line \n";
      }
      elsif ($OptionsInfo{DetailLevel} >= 3) {
	print "Line number $LineCount contains ", scalar(@InvalidColLabels)," non-numerical or empty value(s) for column(s) - ", JoinWords(\@InvalidColLabels, ", ", 0)," - to be analyzed...\n";
      }
      elsif ($OptionsInfo{DetailLevel} >= 2) {
	print "Line number $LineCount contains ", scalar(@InvalidColLabels)," non-numerical or empty value(s) for columns to be analyzed...\n";
      }
    }
  }
  if ($InvalidLineCount && ($OptionsInfo{DetailLevel} >= 1)) {
    print "Non-numerical or empty data present in $InvalidLineCount line(s)...\n";
  }
  close TEXTFILE;

  # Perform the analysis...
  my(@SpecifiedFunctionNames, $SpecifiedFunction);
  @SpecifiedFunctionNames = ();

  for $SpecifiedFunction (@{$OptionsInfo{SpecifiedStatisticalFunctions}}) {
    if ($SpecifiedFunction !~ /^(Covariance|Correlation|Frequency|Rsquare|StandardScores|StandardScoresN)$/i) {
      push @SpecifiedFunctionNames, $OptionsInfo{SpecifiedStatisticalFunctionsMap}{lc($SpecifiedFunction)};
    }
  }
  if (@SpecifiedFunctionNames) {
    PerformAnalysis($Index, \@SpecifiedFunctionNames, \%ColValuesToAnalyzeMap)
  }
  if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{covariance}) || exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{correlation}) || exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{rsquare})) {
    if ($OptionsInfo{AllColumnPairs}) {
      PerformMatrixAnalysis($Index, \%ColValuesToAnalyzeMap);
    }
    else {
      # Perform pairwise analysis for specified columns and write out calculated values - correlation
      # rsquare, or covariance - in the same file.
      PerformColumnPairAnalysis($Index, \%ColValuesToAnalyzeMap);
    }
  }
  if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{standardscores}) || exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{standardscoresn}) ) {
    PerformStandardScoresAnalysis($Index, \%ColValuesToAnalyzeMap);
  }
  if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{frequency})) {
    PerformFrequencyAnalysis($Index, \%ColValuesToAnalyzeMap);
  }
}

# Calculate values for various statistical functions...
sub PerformAnalysis {
  my($Index, $SpecifiedFunctionNamesRef, $ColValuesToAnalyzeMapRef) = @_;
  my($NewTextFile, $Line, $SpecifiedFunction, $Label, @ColLabels, @ColNumsToAnalyze);

  $NewTextFile = $TextFilesInfo{OutFileRoot}[$Index] . $OptionsInfo{FileNameMode} . "." . $TextFilesInfo{OutFileExt}[$Index];

  print "Generating new text file $NewTextFile...\n";
  open NEWTEXTFILE, ">$NewTextFile" or die "Error: Can't open $NewTextFile: $! \n";

  # Write out column labels...
  @ColLabels = ();
  push @ColLabels, "ColumnID";
  for $SpecifiedFunction (@{$SpecifiedFunctionNamesRef}) {
    $Label = $SpecifiedFunction;
    if ($SpecifiedFunction =~ /^(KLargest|KSmallest)$/i) {
      my($KthValue);
      $KthValue = ($SpecifiedFunction =~ /^KLargest$/i) ? $OptionsInfo{KLargest} : $OptionsInfo{KSmallest};
      $Label = AddNumberSuffix($KthValue) . "$SpecifiedFunction";
      $Label =~ s/K//g;
    }
    elsif ($SpecifiedFunction =~ /^TrimMean$/i) {
      $Label = "${SpecifiedFunction}($OptionsInfo{TrimFraction})";
    }
    push @ColLabels, $Label;
  }
  $Line = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print NEWTEXTFILE "$Line\n";

  # Go over each column to be analyzed...
  @ColNumsToAnalyze = @{$TextFilesInfo{ColNumsToAnalyze}[$Index]};

  # Turn off "strict"; otherwise, invoking statistical functions using function name string
  # is problematic.
  no strict;

  my($ColValuesRef, $ColNum, $Value, @RowValues, %CalculatedValues);
  %CalculatedValues = ();
  for $ColNum (@ColNumsToAnalyze) {
    @RowValues = ();
    # Setup column id...
    push @RowValues, $TextFilesInfo{ColLabels}[$Index][$ColNum];
    $ColValuesRef =  \@{$ColValuesToAnalyzeMapRef->{$ColNum}};
    FUNCTIONNAME: for $SpecifiedFunction (@{$SpecifiedFunctionNamesRef}) {
      $Value = "";
      if (!@{$ColValuesToAnalyzeMapRef->{$ColNum}}) {
	# Invalid column values...
	push @RowValues, $Value;
	next FUNCTIONNAME;
      }
      if ($SpecifiedFunction =~ /^Count$/i) {
	$Value = @{$ColValuesToAnalyzeMapRef->{$ColNum}};
      }
      elsif ($SpecifiedFunction =~ /^KLargest$/i) {
	$Value = &$SpecifiedFunction($ColValuesRef, $OptionsInfo{KLargest});
      }
      elsif ($SpecifiedFunction =~ /^KSmallest$/i) {
	$Value = &$SpecifiedFunction($ColValuesRef, $OptionsInfo{KSmallest});
      }
      elsif ($SpecifiedFunction =~ /^StandardDeviation$/i) {
	if (exists($CalculatedValues{$ColNum}{StandardDeviation})) {
	  $Value = $CalculatedValues{$ColNum}{StandardDeviation};
	}
	else {
	  $Value = &$SpecifiedFunction($ColValuesRef);
	  $CalculatedValues{$ColNum}{StandardDeviation} = $Value;
	}
      }
      elsif ($SpecifiedFunction =~ /^StandardError$/i) {
	if (!exists($CalculatedValues{$ColNum}{StandardDeviation})) {
	  $Value = StandardDeviation($ColValuesRef);
	  $CalculatedValues{$ColNum}{StandardDeviation} = $Value;
	}
	if (defined $CalculatedValues{$ColNum}{StandardDeviation}) {
	  $Value = &$SpecifiedFunction($CalculatedValues{$ColNum}{StandardDeviation}, @{$ColValuesToAnalyzeMapRef->{$ColNum}});
	}
      }
      elsif ($SpecifiedFunction =~ /^TrimMean$/i) {
	$Value = &$SpecifiedFunction($ColValuesRef, $OptionsInfo{TrimFraction});
      }
      else {
	$Value = &$SpecifiedFunction($ColValuesRef);
      }
      # Format the output value. And add zero to get rid of tariling zeros...
      $Value = (defined($Value) && length($Value)) ? (sprintf("%.$OptionsInfo{Precision}f", $Value) + 0) : "";
      push @RowValues, $Value;
    }
    $Line = JoinWords(\@RowValues, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print NEWTEXTFILE "$Line\n";
  }
  close NEWTEXTFILE;
}

# Calculate covariance, correlation, rsquare for specified column pairs....
sub PerformColumnPairAnalysis {
  my($Index, $ColValuesToAnalyzeMapRef) = @_;
  my($NewTextFile, @ColLabels, $Line, $CalculateCorrelation, $CalculateRSquare, $CalculateCovariance);
  $CalculateCorrelation = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{correlation}) ? 1 : 0;
  $CalculateRSquare = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{rsquare}) ? 1 : 0;
  $CalculateCovariance = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{covariance}) ? 1 : 0;

  $NewTextFile = $TextFilesInfo{OutFileRoot}[$Index] . "ColumnPairsAnalysis." .  $TextFilesInfo{OutFileExt}[$Index];
  print "Generating new text file $NewTextFile...\n";
  open NEWTEXTFILE, ">$NewTextFile" or die "Error: Can't open $NewTextFile: $! \n";

  # Write out the column labels...
  @ColLabels = ();
  push @ColLabels, ("ColumnID1", "ColumnID2");
  if ($CalculateCorrelation || $CalculateRSquare) {
    push @ColLabels, "Correlation";
    if ($CalculateRSquare) {
      push @ColLabels, "RSquare";
    }
  }
  if ($CalculateCovariance) {
    push @ColLabels, "Covariance";
  }
  $Line = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print NEWTEXTFILE "$Line\n";

  # Go over each column pair...
  my($CorrelationValue, $RSquareValue, $CovarianceValue,  $ColIndex, $ColNum1, $ColNum2, $ColValuesRef1, $ColValuesRef2, @ColPairs1ToAnalyze, @ColPairs2ToAnalyze, @RowValues, $Value);

  @ColPairs1ToAnalyze = @{$TextFilesInfo{ColPairs1ToAnalyze}[$Index]};
  @ColPairs2ToAnalyze = @{$TextFilesInfo{ColPairs2ToAnalyze}[$Index]};
  for $ColIndex (0 .. $#ColPairs1ToAnalyze) {
    @RowValues = ();
    $ColNum1 = $ColPairs1ToAnalyze[$ColIndex];
    $ColNum2 = $ColPairs2ToAnalyze[$ColIndex];
    $ColValuesRef1 =  \@{$ColValuesToAnalyzeMapRef->{$ColNum1}};
    $ColValuesRef2 =  \@{$ColValuesToAnalyzeMapRef->{$ColNum2}};

    # Setup column ids...
    push @RowValues, $TextFilesInfo{ColLabels}[$Index][$ColNum1];
    push @RowValues, $TextFilesInfo{ColLabels}[$Index][$ColNum2];

    if (@$ColValuesRef1 != @$ColValuesRef2) {
      # Print a warning...
      warn "Warning: Skipping analysis for column pair $TextFilesInfo{ColLabels}[$Index][$ColNum1], $TextFilesInfo{ColLabels}[$Index][$ColNum2]: Number of valid data values must be same.\n";
      if ($CalculateCorrelation || $CalculateRSquare) {
	push @RowValues, "";
	if ($CalculateRSquare) {
	  push @RowValues, "";
	}
      }
      if ($CalculateCovariance) {
	push @RowValues, "";
      }
    }
    else {
      # Calculate appropriate value...
      if ($CalculateCorrelation || $CalculateRSquare) {
	$CorrelationValue = Correlation($ColValuesRef1, $ColValuesRef2);
	$Value = (defined($CorrelationValue) && length($CorrelationValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $CorrelationValue) + 0) : "";
	push @RowValues, $Value;
	if ($CalculateRSquare) {
	  $RSquareValue = (defined($CorrelationValue) && length($CorrelationValue)) ? ($CorrelationValue ** 2) : "";
	  $Value = (length($RSquareValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $RSquareValue) + 0) : "";
	  push @RowValues, $Value;
	}
      }
      if ($CalculateCovariance) {
	$CovarianceValue = Covariance($ColValuesRef1, $ColValuesRef2);
	$Value = (defined($CovarianceValue) && length($CovarianceValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $CovarianceValue) + 0) : "";
	push @RowValues, $Value;
      }
    }
    $Line = JoinWords(\@RowValues, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print NEWTEXTFILE "$Line\n";
  }
  close NEWTEXTFILE;
}

# Generate histogram numbers...
sub PerformFrequencyAnalysis {
  my($Index, $ColValuesToAnalyzeMapRef) = @_;
  my($NewTextFile, $ColLabel, @ColLabels, @RowValues, $Line, $ColNum, @ColNumsToAnalyze, $ColValuesRef, $BinValue, $FrequencyValue, $Value, %FrequencyMap);

  @ColNumsToAnalyze = @{$TextFilesInfo{ColNumsToAnalyze}[$Index]};
  for $ColNum (@ColNumsToAnalyze) {
    $NewTextFile = $TextFilesInfo{OutFileRoot}[$Index] . $TextFilesInfo{ColLabels}[$Index][$ColNum] . "FrequencyAnalysis." .  $TextFilesInfo{OutFileExt}[$Index];
    print "Generating new text file $NewTextFile...\n";
    open NEWTEXTFILE, ">$NewTextFile" or die "Error: Can't open $NewTextFile: $! \n";

    # Write out the column labels...
    @ColLabels = ();
    push @ColLabels , ("Bins", "Frequency");
    $Line = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print NEWTEXTFILE "$Line\n";

    #Calculate and write out frequency values...
    %FrequencyMap = ();
    $ColValuesRef =  \@{$ColValuesToAnalyzeMapRef->{$ColNum}};
    if (@$ColValuesRef) {
      if (@{$OptionsInfo{BinRange}}) {
	%FrequencyMap = Frequency($ColValuesRef, \@{$OptionsInfo{BinRange}});
      }
      else {
	%FrequencyMap = Frequency($ColValuesRef, $OptionsInfo{NumOfBins});
      }
    }
    for $BinValue (sort { $a <=> $b }  keys %FrequencyMap) {
      $FrequencyValue = $FrequencyMap{$BinValue};

      @RowValues = ();
      $Value = (length($BinValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $BinValue) + 0) : "";
      push @RowValues, $Value;
      $Value = (length($FrequencyValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $FrequencyValue) + 0) : "";
      push @RowValues, $Value;

      $Line = JoinWords(\@RowValues, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
      print NEWTEXTFILE "$Line\n";
    }
    close NEWTEXTFILE;
  }
}

# Calculate covariance, correlation/rsquare matrices....
sub PerformMatrixAnalysis {
  my($Index, $ColValuesToAnalyzeMapRef) = @_;
  my($CorrelationTextFile, $CovarianceTextFile, $RSquareTextFile, $CalculateCorrelation, $CalculateRSquare, $CalculateCovariance);

  $CalculateCorrelation = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{correlation}) ? 1 : 0;
  $CalculateRSquare = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{rsquare}) ? 1 : 0;
  $CalculateCovariance = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{covariance}) ? 1 : 0;

  $CorrelationTextFile = $TextFilesInfo{OutFileRoot}[$Index] . "CorrelationMatrix." .  $TextFilesInfo{OutFileExt}[$Index];
  $RSquareTextFile = $TextFilesInfo{OutFileRoot}[$Index] . "RSquareMatrix." .  $TextFilesInfo{OutFileExt}[$Index];
  $CovarianceTextFile = $TextFilesInfo{OutFileRoot}[$Index] . "CovarianceMatrix." .  $TextFilesInfo{OutFileExt}[$Index];

  my($TextFilesList, $Delimiter);
  $TextFilesList =  "";
  if ($CalculateCorrelation || $CalculateRSquare) {
    $TextFilesList = $CorrelationTextFile;
    if ($CalculateRSquare) {
      $TextFilesList .= ", $CorrelationTextFile";
    }
  }
  $Delimiter = length($TextFilesList) ? "," : "";
  if ($CalculateCovariance) {
    $TextFilesList .= "${Delimiter} ${CorrelationTextFile}";
  }
  if ($TextFilesList =~ /\,/) {
    print "Generating new text files $TextFilesList...\n"
  }
  else {
    print "Generating new text file $TextFilesList...\n"
  }
  if ($CalculateCorrelation || $CalculateRSquare) {
    open CORRELATIONTEXTFILE, ">$CorrelationTextFile" or die "Error: Can't open $CorrelationTextFile: $! \n";
    if ($CalculateRSquare) {
      open RSQUARETEXTFILE, ">$RSquareTextFile" or die "Error: Can't open $RSquareTextFile: $! \n";
    }
  }
  if ($CalculateCovariance) {
    open COVARIANCETEXTFILE, ">$CovarianceTextFile" or die "Error: Can't open $CovarianceTextFile: $! \n";
  }

  my($Line, $Value, $CorrelationValue, $RSquareValue, $CovarianceValue, $ColNum, $ColNum1, $ColNum2, $ColValuesRef1, $ColValuesRef2, @ColLabels, @CovarianceRowValues, @CorrelationRowValues, @RSquareRowValues);

  # Write out the column labels...
  @ColLabels = ();
  push @ColLabels, "";
  for $ColNum (0 .. ($TextFilesInfo{ColCount}[$Index] - 1)) {
    push @ColLabels, $TextFilesInfo{ColLabels}[$Index][$ColNum];
  }
  $Line = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  if ($CalculateCorrelation || $CalculateRSquare) {
    print CORRELATIONTEXTFILE "$Line\n";
    if ($CalculateRSquare) {
      print RSQUARETEXTFILE "$Line\n";
    }
  }
  if ($CalculateCovariance) {
    print COVARIANCETEXTFILE "$Line\n";
  }

  # Due to symmetric nature of these matrices, only one half needs to be
  # calculated. So, just calculate the lower half and copy it to upper half...
  my(%CorrelationMatrixMap, %RSquareMatrixMap, %CovarianceMatrixMap);

  %CorrelationMatrixMap = (); %RSquareMatrixMap = (); %CovarianceMatrixMap = ();
  for $ColNum1 (0 .. ($TextFilesInfo{ColCount}[$Index] - 1)) {
    for $ColNum2 (0 .. $ColNum1) {
      $ColValuesRef1 =  \@{$ColValuesToAnalyzeMapRef->{$ColNum1}};
      $ColValuesRef2 =  \@{$ColValuesToAnalyzeMapRef->{$ColNum2}};
      if ($CalculateCorrelation || $CalculateRSquare) {
	$CorrelationValue = Correlation($ColValuesRef1, $ColValuesRef2);
	$CorrelationValue = (defined($CorrelationValue) && length($CorrelationValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $CorrelationValue) + 0) : "";
	$CorrelationMatrixMap{$ColNum1}{$ColNum2} = $CorrelationValue;
	if ($ColNum1 != $ColNum2) {
	  $CorrelationMatrixMap{$ColNum2}{$ColNum1} = $CorrelationValue;
	}
	if ($CalculateRSquare) {
	  $RSquareValue = (defined($CorrelationValue) && length($CorrelationValue)) ? ($CorrelationValue ** 2) : "";
	  $RSquareValue = (length($RSquareValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $RSquareValue) + 0) : "";
	  $RSquareMatrixMap{$ColNum1}{$ColNum2} = $RSquareValue;
	  if ($ColNum1 != $ColNum2) {
	    $RSquareMatrixMap{$ColNum2}{$ColNum1} = $RSquareValue;
	  }
	}
      }
      if ($CalculateCovariance) {
	$CovarianceValue = Covariance($ColValuesRef1, $ColValuesRef2);
	$CovarianceValue = (defined($CovarianceValue) && length($CovarianceValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $CovarianceValue) + 0) : "";
	$CovarianceMatrixMap{$ColNum1}{$ColNum2} = $CovarianceValue;
	if ($ColNum1 != $ColNum2) {
	  $CovarianceMatrixMap{$ColNum2}{$ColNum1} = $CovarianceValue;
	}
      }
    }
  }

  # Write out the matrices...
  for $ColNum1 (0 .. ($TextFilesInfo{ColCount}[$Index] - 1)) {
    @CorrelationRowValues = ();
    @RSquareRowValues = ();
    @CovarianceRowValues = ();
    if ($CalculateCorrelation || $CalculateRSquare) {
      push @CorrelationRowValues, $TextFilesInfo{ColLabels}[$Index][$ColNum1];
      if ($CalculateRSquare) {
	push @RSquareRowValues, $TextFilesInfo{ColLabels}[$Index][$ColNum1];
      }
    }
    if ($CalculateCovariance) {
      push @CovarianceRowValues, $TextFilesInfo{ColLabels}[$Index][$ColNum1];
    }
    for $ColNum2 (0 .. ($TextFilesInfo{ColCount}[$Index] - 1)) {
      if ($CalculateCorrelation || $CalculateRSquare) {
	push @CorrelationRowValues, $CorrelationMatrixMap{$ColNum1}{$ColNum2};
	if ($CalculateRSquare) {
	  push @RSquareRowValues, $RSquareMatrixMap{$ColNum1}{$ColNum2};
	}
      }
      if ($CalculateCovariance) {
	push @CovarianceRowValues, $CovarianceMatrixMap{$ColNum1}{$ColNum2};
      }
    }
    if ($CalculateCorrelation || $CalculateRSquare) {
      $Line = JoinWords(\@CorrelationRowValues, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
      print CORRELATIONTEXTFILE "$Line\n";
      if ($CalculateRSquare) {
	$Line = JoinWords(\@RSquareRowValues, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
	print RSQUARETEXTFILE "$Line\n";
      }
    }
    if ($CalculateCovariance) {
      $Line = JoinWords(\@CovarianceRowValues, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
      print COVARIANCETEXTFILE "$Line\n";
    }
  }
  if ($CalculateCorrelation || $CalculateRSquare) {
    close CORRELATIONTEXTFILE;
    if ($CalculateRSquare) {
      close RSQUARETEXTFILE;
    }
  }
  if ($CalculateCovariance) {
    close COVARIANCETEXTFILE;
  }
}

# Calculate standard scores...
sub PerformStandardScoresAnalysis {
  my($Index, $ColValuesToAnalyzeMapRef) = @_;
  my($StandardScores, $StandardScoresN, $NewTextFile, @ColLabels, $Label, $NewLine);

  $StandardScores = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{standardscores}) ? 1 : 0;
  $StandardScoresN = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{standardscoresn}) ? 1 : 0;

  $NewTextFile = $TextFilesInfo{OutFileRoot}[$Index] . "StandardScores." .  $TextFilesInfo{OutFileExt}[$Index];
  print "Generating new text file $NewTextFile...\n";
  open NEWTEXTFILE, ">$NewTextFile" or die "Error: Can't open $NewTextFile: $! \n";

  my($ColValuesRef, $ColNum, @ColNumsToAnalyze);
  # Write out column labels...
  @ColLabels = ();
  @ColNumsToAnalyze = @{$TextFilesInfo{ColNumsToAnalyze}[$Index]};
  for $ColNum (@ColNumsToAnalyze) {
    $Label = $TextFilesInfo{ColLabels}[$Index][$ColNum];
    if ($StandardScores) {
      push @ColLabels, "${Label}\(StandardScores)";
    }
    if ($StandardScoresN) {
      push @ColLabels, "${Label}\(StandardScoresN)";
    }
  }
  $NewLine = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print NEWTEXTFILE "$NewLine\n";

  # Go over each column to be analyzed and calculate standard deviation
  # and mean values...
  my(%StandardDeviationMap, %StandardDeviationNMap, %MeanMap);
  %StandardDeviationMap = ();
  %StandardDeviationNMap = ();
  %MeanMap = ();
  for $ColNum (@ColNumsToAnalyze) {
    $ColValuesRef =  \@{$ColValuesToAnalyzeMapRef->{$ColNum}};
    if (!exists($MeanMap{$ColNum})) {
      $MeanMap{$ColNum} = Mean($ColValuesRef);
    }
    if ($StandardScores) {
      if (!exists($StandardDeviationMap{$ColNum})) {
	$StandardDeviationMap{$ColNum} = StandardDeviation($ColValuesRef);
      }
    }
    if ($StandardScoresN) {
      if (!exists($StandardDeviationNMap{$ColNum})) {
	$StandardDeviationNMap{$ColNum} = StandardDeviationN($ColValuesRef);
      }
    }
  }
  #
  # Go over each row and calculate standard scores for each column
  # using (x[i] - mean) / (n - 1) for StandardScores and (x[i] - mean) / n
  # for StandardScoresN; write out the calculated values as well...

  my($TextFile, $InDelim, $Line, $Value, $ValueOkay, $ScoreValue, @RowValues, @LineWords);
  $TextFile = $TextFilesList[$Index];
  $InDelim = $TextFilesInfo{InDelim}[$Index];

  open TEXTFILE, "$TextFile" or die "Error: Can't open $TextFile: $! \n";
  $Line = GetTextLine(\*TEXTFILE);
  while ($Line = GetTextLine(\*TEXTFILE)) {
    @LineWords = quotewords($InDelim, 0, $Line);
    @RowValues = ();
    COLNUM: for $ColNum (@ColNumsToAnalyze) {
      $Value = $LineWords[$ColNum];
      $ValueOkay = ($OptionsInfo{CheckData} && !IsNumerical($Value)) ? 0 : 1;
      if ($StandardScores) {
	$ScoreValue = $ValueOkay ? (($Value - $MeanMap{$ColNum})/$StandardDeviationMap{$ColNum}) : "";
	$ScoreValue = (defined($ScoreValue) && length($ScoreValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $ScoreValue) + 0) : "";
	push @RowValues, $ScoreValue;
      }
      if ($StandardScoresN) {
	$ScoreValue = $ValueOkay ? (($Value - $MeanMap{$ColNum})/$StandardDeviationNMap{$ColNum}) : "";
	$ScoreValue = (defined($ScoreValue) && length($ScoreValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $ScoreValue) + 0) : "";
	push @RowValues, $ScoreValue;
      }
    }
    $NewLine = JoinWords(\@RowValues, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print NEWTEXTFILE "$NewLine\n";
  }
  close TEXTFILE;
  close NEWTEXTFILE;
}

# Make sure the specified columns exists in text files...
sub ProcessColumnsInfo {
  my($Index, $TextFile, $ColNum, $NewColNum, $ColIndex, @ColNumsToAnalyze, %UniqueColNumsToAnalyzeMap);

  @{$TextFilesInfo{ColNumsToAnalyze}} = ();
  @{$TextFilesInfo{ColPairs1ToAnalyze}} = ();
  @{$TextFilesInfo{ColPairs2ToAnalyze}} = ();
  @{$TextFilesInfo{UniqueColNumsToAnalyze}} = ();

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];

    @{$TextFilesInfo{ColNumsToAnalyze}[$Index]} = ();
    @{$TextFilesInfo{ColPairs1ToAnalyze}[$Index]} = ();
    @{$TextFilesInfo{ColPairs2ToAnalyze}[$Index]} = ();
    @{$TextFilesInfo{UniqueColNumsToAnalyze}[$Index]} = ();

    %UniqueColNumsToAnalyzeMap = ();

    if ($TextFilesInfo{FileOkay}[$Index]) {
      @ColNumsToAnalyze = ();
      if (@{$OptionsInfo{SpecifiedColumns}}) {
	if ($OptionsInfo{ColMode} =~ /^colnum$/i) {
	  for $ColNum (@{$OptionsInfo{SpecifiedColumns}}) {
	    if ($ColNum >=1 && $ColNum <= $TextFilesInfo{ColCount}[$Index]) {
	      $NewColNum = $ColNum -1;
	      push @ColNumsToAnalyze, $NewColNum;
	    }
	  }
	}
	else {
	  my($ColLabel);
	  for $ColLabel (@{$OptionsInfo{SpecifiedColumns}}) {
	    if (exists($TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel})) {
	      push @ColNumsToAnalyze, $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel};
	    }
	  }
	}
      }
      elsif (defined  $OptionsInfo{Columns} && $OptionsInfo{Columns} =~ /^All$/i) {
	for $ColNum (0 .. ($TextFilesInfo{ColCount}[$Index] - 1)) {
	  push @ColNumsToAnalyze, $ColNum;
	}
      }
      else {
	push @ColNumsToAnalyze, 0;
      }
      if (@ColNumsToAnalyze) {
	push @{$TextFilesInfo{ColNumsToAnalyze}[$Index]}, @ColNumsToAnalyze;
	# Set up unique columns map as well...
	for $ColNum (@ColNumsToAnalyze) {
	  if (!exists $UniqueColNumsToAnalyzeMap{$ColNum}) {
	    $UniqueColNumsToAnalyzeMap{$ColNum} = $ColNum;
	  }
	}
      }
      else {
	warn "Warning: Ignoring file $TextFile: None of the columns specified, @{$OptionsInfo{SpecifiedColumns}}, using \"--columns\" option exist.\n";
	$TextFilesInfo{FileOkay}[$Index] = 0;
	next FILELIST;
      }
      if (!$OptionsInfo{Overwrite} && exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{frequency})) {
	# Make sure specific frequency files don't exist...
	my($FrequencyFile);
	for $ColNum (@ColNumsToAnalyze) {
	  $FrequencyFile = $TextFilesInfo{OutFileRoot}[$Index] . $TextFilesInfo{ColLabels}[$Index][$ColNum] . "FrequencyAnalysis." .  $TextFilesInfo{OutFileExt}[$Index];
	  if (-e $FrequencyFile) {
	    warn "Warning: Ignoring file $TextFile: The file $FrequencyFile already exists.\n";
	    $TextFilesInfo{FileOkay}[$Index] = 0;
	    next FILELIST;
	  }
	}
      }
      # Setup specified column pairs...
      if (exists $OptionsInfo{SpecifiedStatisticalFunctionsMap}{correlation} || exists $OptionsInfo{SpecifiedStatisticalFunctionsMap}{covariance} || exists $OptionsInfo{SpecifiedStatisticalFunctionsMap}{rsquare}) {
	my(@ColPairsToAnalyze, $ColNum1, $ColNum2);
	if (@{$OptionsInfo{SpecifiedColumnPairs}}) {
	  # Make sure both columns exist...
	  if ($OptionsInfo{ColMode} =~ /^colnum$/i) {
	    for ($ColIndex = 0; (($ColIndex + 1) < @{$OptionsInfo{SpecifiedColumnPairs}}); $ColIndex += 2 ) {
	      $ColNum1 = $OptionsInfo{SpecifiedColumnPairs}[$ColIndex];
	      $ColNum2 = $OptionsInfo{SpecifiedColumnPairs}[$ColIndex + 1];
	      if ($ColNum1 >=1 && $ColNum1 <= $TextFilesInfo{ColCount}[$Index] && $ColNum2 >=1 && $ColNum2 <= $TextFilesInfo{ColCount}[$Index]) {
		$ColNum1 -= 1;
		$ColNum2 -= 1;
		push @ColPairsToAnalyze, ($ColNum1, $ColNum2);
	      }
	    }
	  }
	  else {
	    my($ColLabel1, $ColLabel2);
	    for ($ColIndex = 0; (($ColIndex + 1) < @{$OptionsInfo{SpecifiedColumnPairs}}); $ColIndex += 2 ) {
	      $ColLabel1 = $OptionsInfo{SpecifiedColumnPairs}[$ColIndex];
	      $ColLabel2 = $OptionsInfo{SpecifiedColumnPairs}[$ColIndex + 1];
	      if (exists($TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel1}) && exists($TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel2})) {
		$ColNum1 = $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel1};
		$ColNum2 = $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel2};
		push @ColPairsToAnalyze, ($ColNum1, $ColNum2);
	      }
	    }
	  }
	}
	elsif ($OptionsInfo{AllColumnPairs}) {
	  for $ColNum1 (0 .. ($TextFilesInfo{ColCount}[$Index] - 1)) {
	    for $ColNum2 (0 .. ($TextFilesInfo{ColCount}[$Index] - 1)) {
	      push @ColPairsToAnalyze, ($ColNum1, $ColNum2);
	    }
	  }
	}
	else {
	  if ($TextFilesInfo{ColCount}[$Index] >= 2) {
	    push @ColPairsToAnalyze, (0,1);
	  }
	}
	if (@ColPairsToAnalyze) {
	  if (@ColPairsToAnalyze % 2) {
	    warn "Warning: Ignoring file $TextFile: Invalid number  values specified using \"--columnpairs\" option: It must contain even number of valid values.\n";
	    $TextFilesInfo{FileOkay}[$Index] = 0;
	    next FILELIST;
	  }
	  else {
	    for ($ColIndex = 0; $ColIndex < @ColPairsToAnalyze; $ColIndex += 2) {
	      push @{$TextFilesInfo{ColPairs1ToAnalyze}[$Index]}, $ColPairsToAnalyze[$ColIndex];
	      push @{$TextFilesInfo{ColPairs2ToAnalyze}[$Index]}, $ColPairsToAnalyze[$ColIndex + 1];
	    }
	    # Set up unique columns map as well...
	    for $ColNum (@ColPairsToAnalyze) {
	      if (!exists $UniqueColNumsToAnalyzeMap{$ColNum}) {
		$UniqueColNumsToAnalyzeMap{$ColNum} = $ColNum;
	      }
	    }
	  }
	}
      }
      # Setup uniques columns array...
      push @{$TextFilesInfo{UniqueColNumsToAnalyze}[$Index]}, (sort keys %UniqueColNumsToAnalyzeMap);
    }
  }
}

# Retrieve information about input text files...
sub RetrieveTextFilesInfo {
  my($Index, $TextFile, $FileDir, $FileName, $FileExt, $InDelim, $Line, @ColLabels, $OutFileRoot,  $OutFile, $OutFileExt, $ColNum, $ColLabel);

  %TextFilesInfo = ();

  @{$TextFilesInfo{FileOkay}} = ();
  @{$TextFilesInfo{ColCount}} = ();
  @{$TextFilesInfo{ColLabels}} = ();
  @{$TextFilesInfo{ColLabelToNumMap}} = ();
  @{$TextFilesInfo{InDelim}} = ();
  @{$TextFilesInfo{OutFileRoot}} = ();
  @{$TextFilesInfo{OutFileExt}} = ();

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];

    $TextFilesInfo{FileOkay}[$Index] = 0;
    $TextFilesInfo{ColCount}[$Index] = 0;
    $TextFilesInfo{InDelim}[$Index] = "";
    $TextFilesInfo{OutFileRoot}[$Index] = "";
    $TextFilesInfo{OutFileExt}[$Index] = "";

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
      if ($Options{indelim} !~ /^(comma|semicolon)$/i) {
	warn "Warning: Ignoring file $TextFile: The value specified, $Options{indelim}, for option \"--indelim\" is not valid for csv files\n";
	next FILELIST;
      }
      if ($Options{indelim} =~ /^semicolon$/i) {
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
    $OutFileExt = $FileExt;
    if ($Options{root} && (@TextFilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($Options{root});
      if ($RootFileName && $RootFileExt) {
	$FileName = $RootFileName;
      }
      else {
	$FileName = $Options{root};
      }
      $OutFileRoot = $FileName;
    }
    else {
      $OutFileRoot = $FileName;
    }
    $OutFile = $OutFileRoot . $OptionsInfo{FileNameMode} . ".$OutFileExt";

    if (lc($OutFile) eq lc($TextFile)) {
      warn "Warning: Ignoring file $TextFile:Output file name, $OutFile, is same as input text file name, $TextFile\n";
      next FILELIST;
    }
    if (!$Options{overwrite}) {
      if (-e $OutFile) {
	warn "Warning: Ignoring file $TextFile: The file $OutFile already exists\n";
	next FILELIST;
      }
      if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{covariance}) || exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{correlation}) || exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{rsquare})) {
	if ($OptionsInfo{AllColumnPairs}) {
	  if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{covariance}) && (-e "${OutFileRoot}CovarianceMatrix.${FileExt}")) {
	    warn "Warning: Ignoring file $TextFile: The file ${OutFileRoot}Covariance.${FileExt} already exists.\n";
	    next FILELIST;
	  }
	  if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{correlation}) && (-e "${OutFileRoot}CorrelationMatrix.${FileExt}")) {
	    warn "Warning: Ignoring file $TextFile: The file ${OutFileRoot}CorrelationMatrix.${FileExt} already exists.\n";
	    next FILELIST;
	  }
	  if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{rsquare}) && (-e "${OutFileRoot}RSquareMatrix.${FileExt}")) {
	    warn "Warning: Ignoring file $TextFile: The file ${OutFileRoot}RSquareMatrix.${FileExt} already exists.\n";
	    next FILELIST;
	  }
	}
	else {
	  if (-e "${OutFileRoot}ColumnPairsAnalysis.${FileExt}") {
	    warn "Warning: Ignoring file $TextFile: The file ${OutFileRoot}ColumnPairsAnalysis.${FileExt} already exists.\n";
	    next FILELIST;
	  }
	}
      }
      if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{standardscores}) && (-e "${OutFileRoot}StandardScores.${FileExt}")) {
	warn "Warning: Ignoring file $TextFile: The file ${OutFileRoot}StandardScores.${FileExt} already exists.\n";
	next FILELIST;
      }
    }

    $TextFilesInfo{FileOkay}[$Index] = 1;
    $TextFilesInfo{InDelim}[$Index] = $InDelim;
    $TextFilesInfo{OutFileRoot}[$Index] = "$OutFileRoot";
    $TextFilesInfo{OutFileExt}[$Index] = "$OutFileExt";

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
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{DetailLevel} = $Options{detail};

  # Setup supported statistical functions...
  my($SupportedFunction, @SupportedStatisticaFunctions, %SupportedStatisticaFunctionsMap);
  %SupportedStatisticaFunctionsMap = ();
  @SupportedStatisticaFunctions = qw(Average AverageDeviation Correlation Count Covariance GeometricMean Frequency HarmonicMean KLargest KSmallest Kurtosis Maximum Minimum Mean Median Mode RSquare Skewness Sum SumOfSquares StandardDeviation StandardDeviationN StandardError StandardScores StandardScoresN TrimMean Variance VarianceN);

  for $SupportedFunction (@SupportedStatisticaFunctions) {
    $SupportedStatisticaFunctionsMap{lc($SupportedFunction)} = $SupportedFunction;
  }

  # Setup a list of functions to use for analysis...
  my($SpecifiedFunction);
  %{$OptionsInfo{SpecifiedStatisticalFunctionsMap}} = ();
  @{$OptionsInfo{SpecifiedStatisticalFunctions}} = ();
  # Check mode values...
  if ($Options{mode} =~ /^DescriptiveStatisticsBasic$/i ) {
    $OptionsInfo{FileNameMode} = "DescriptiveStatisticsBasic";
    @{$OptionsInfo{SpecifiedStatisticalFunctions}} = qw(Count Maximum Minimum Mean Median StandardDeviation StandardError Variance Sum);
  }
  elsif ($Options{mode} =~ /^DescriptiveStatisticsAll$/i ) {
    $OptionsInfo{FileNameMode} = "DescriptiveStatisticsAll";
    @{$OptionsInfo{SpecifiedStatisticalFunctions}} = qw(Count Maximum Minimum Mean GeometricMean HarmonicMean TrimMean Median Mode StandardDeviation Kurtosis Skewness StandardError Variance  RSquare Frequency  KLargest KSmallest Sum);
  }
  elsif ($Options{mode} =~ /^All$/i ) {
    $OptionsInfo{FileNameMode} = "AllStatistics";
    @{$OptionsInfo{SpecifiedStatisticalFunctions}} = @SupportedStatisticaFunctions;
  }
  else {
    $OptionsInfo{FileNameMode} = "SpecifiedStatistics";
    # Comma delimited list of functions...
    my($Mode, @SpecifiedFunctions, @UnsupportedSpecifiedFunctions);
    $Mode = $Options{mode};
    $Mode =~ s/ //g;
    @SpecifiedFunctions = split ",", $Mode;
    @UnsupportedSpecifiedFunctions = ();
    for $SpecifiedFunction (@SpecifiedFunctions) {
      if (exists($SupportedStatisticaFunctionsMap{lc($SpecifiedFunction)})) {
	push @{$OptionsInfo{SpecifiedStatisticalFunctions}}, $SpecifiedFunction;
      }
      else {
	push @UnsupportedSpecifiedFunctions, $SpecifiedFunction;
      }
    }
    if (@UnsupportedSpecifiedFunctions) {
      if (@UnsupportedSpecifiedFunctions > 1) {
	warn "Error: The values specified - ", JoinWords(\@UnsupportedSpecifiedFunctions, ", ", 0)," - for option \"-m --mode\" are not valid.\n";
      }
      else {
	warn "Error: The value specified, @UnsupportedSpecifiedFunctions , for option \"-m --mode\" is not valid.\n";
      }
      die "Allowed values:", JoinWords(\@SupportedStatisticaFunctions, ", ", 0), "\n";
    }
  }
  FUNCTION: for $SpecifiedFunction (@{$OptionsInfo{SpecifiedStatisticalFunctions}}) {
    if (exists $OptionsInfo{SpecifiedStatisticalFunctionsMap}{lc($SpecifiedFunction)} ) {
      next FUNCTION;
    }
    $OptionsInfo{SpecifiedStatisticalFunctionsMap}{lc($SpecifiedFunction)} = $SupportedStatisticaFunctionsMap{lc($SpecifiedFunction)};
  }

  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /tab/i ) ? "\t" : (($Options{outdelim} =~ /semicolon/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /yes/i ) ? 1 : 0;

  $OptionsInfo{Overwrite} = defined $Options{overwrite} ? $Options{overwrite} : undef;
  $OptionsInfo{Root} = defined $Options{root} ? $Options{root} : undef;

  $OptionsInfo{CheckData} = $Options{fast} ? 0 : 1;
  $OptionsInfo{Precision} = $Options{precision};

  $OptionsInfo{KLargest} = $Options{klargest};
  $OptionsInfo{KSmallest} = $Options{ksmallest};

  $OptionsInfo{TrimFraction} = $Options{trimfraction};

  # Setup frequency bin values...
  $OptionsInfo{NumOfBins} = 10;
  @{$OptionsInfo{BinRange}} = ();
  if ($Options{frequencybins} =~ /\,/) {
    my($BinValue, @SpecifiedBinRange);
    @SpecifiedBinRange = split /\,/,  $Options{frequencybins};
    if (@SpecifiedBinRange < 2) {
      die "Error: The value specified, $Options{frequencybins}, for option \"--frequencybins\" is not valid: Must contain at least two values. \n";
    }
    for $BinValue (@SpecifiedBinRange) {
      if (!IsNumerical($BinValue)) {
	die "Error: The value specified, $Options{frequencybins}, for option \"--frequencybins\" is not valid: Contains non numeric values. \n";
      }
    }
    my($Index1, $Index2);
    for $Index1 (0 .. $#SpecifiedBinRange) {
      for $Index2 (($Index1 + 1) .. $#SpecifiedBinRange) {
	if ($SpecifiedBinRange[$Index1] >= $SpecifiedBinRange[$Index2]) {
	  die "Error: The value specified, $Options{frequencybins}, for option \"--frequencybins\" is not valid: Must contain values in ascending order. \n";
	}
      }
    }
    push @{$OptionsInfo{BinRange}}, @SpecifiedBinRange;
  }
  else {
    $OptionsInfo{NumOfBins} = $Options{frequencybins};
    if (!IsPositiveInteger($OptionsInfo{NumOfBins})) {
      die "Error: The value specified, $Options{frequencybins}, for option \"--frequencybins\" is not valid. Allowed values: positive integer or \"number,number,[number]...\". \n";
    }
  }

  # Setup specified columns...
  $OptionsInfo{ColMode} = $Options{colmode};
  $OptionsInfo{Columns} = defined $Options{columns} ? $Options{columns} : undef;

  @{$OptionsInfo{SpecifiedColumns}} = ();
  if (defined $Options{columns} && $Options{columns} !~ /^All$/i) {
    my(@SpecifiedValues) = split ",", $Options{columns};
    if ($Options{colmode} =~ /^colnum$/i) {
      my($ColValue);
      for $ColValue (@SpecifiedValues) {
	if (!IsPositiveInteger($ColValue)) {
	  die "Error: Column value, $ColValue, specified using \"--columns\" is not valid: Allowed integer values: > 0.\n";
	}
      }
    }
    push @{$OptionsInfo{SpecifiedColumns}}, @SpecifiedValues;
  }
  @{$OptionsInfo{SpecifiedColumnPairs}} = ();
  $OptionsInfo{AllColumnPairs} = (defined($Options{columnpairs}) && $Options{columnpairs} =~ /^AllPairs$/i) ? 1 : 0;
  if (defined($Options{columnpairs}) && !$OptionsInfo{AllColumnPairs}) {
    my(@SpecifiedValues) = split ",", $Options{columnpairs};
    if (@SpecifiedValues % 2) {
      die "Error: Invalid number of values specified using \"--columnpairs\" option: It must contain even number of values.\n";
    }
    if ($Options{colmode} =~ /^colnum$/i) {
      my($ColValue);
      for $ColValue (@SpecifiedValues) {
	if (!IsPositiveInteger($ColValue)) {
	  die "Error: Column value, $ColValue, specified using \"--columnpairs\" is not valid: Allowed integer values: > 0.\n";
	}
      }
    }
    push @{$OptionsInfo{SpecifiedColumnPairs}}, @SpecifiedValues;
  }

}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{colmode} = "colnum";
  $Options{detail} = 1;
  $Options{indelim} = "comma";
  $Options{frequencybins} = 10;
  $Options{klargest} = 2;
  $Options{ksmallest} = 2;
  $Options{mode} = "DescriptiveStatisticsBasic";
  $Options{outdelim} = "comma";
  $Options{precision} = 2;
  $Options{quote} = "yes";
  $Options{trimfraction} = 0.1;

  if (!GetOptions(\%Options, "colmode|c=s", "columns=s", "columnpairs=s", "detail|d=i", "frequencybins=s", "fast|f", "help|h", "indelim=s", "klargest=i", "ksmallest=i", "mode|m=s", "outdelim=s", "overwrite|o", "precision|p=i", "quote|q=s", "root|r=s", "trimfraction=f", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{colmode} !~ /^(colnum|collabel)$/i) {
    die "Error: The value specified, $Options{colmode}, for option \"-c --colmode\" is not valid. Allowed values: colnum or collabel\n";
  }
  if (!IsPositiveInteger($Options{detail})) {
    die "Error: The value specified, $Options{detail}, for option \"-d --detail\" is not valid. Allowed values: > 0\n";
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
  if (!IsPositiveInteger($Options{precision})) {
    die "Error: The value specified, $Options{precision}, for option \"-p --precision\" is not valid. Allowed values: > 0 \n";
  }
  if (!IsPositiveInteger($Options{klargest})) {
    die "Error: The value specified, $Options{klargest}, for option \"--klargest\" is not valid. Allowed values: > 0 \n";
  }
  if (!IsPositiveInteger($Options{ksmallest})) {
    die "Error: The value specified, $Options{ksmallest}, for option \"--ksmallest\" is not valid. Allowed values: > 0 \n";
  }
  if (IsFloat($Options{trimfraction})) {
    if ($Options{trimfraction} <= 0 || $Options{trimfraction} >= 1.0) {
      die "Error: The value specified, $Options{trimfraction}, for option \"--trimfraction\" is not valid. Allowed values: > 0 and < 1.0\n";
    }
  }
  else {
    die "Error: The value specified, $Options{trimfraction}, for option \"--trimfraction\" is not valid. Allowed values: > 0 and < 1.0\n";
  }
}

__END__

=head1 NAME

AnalyzeTextFilesData.pl - Analyze numerical coulmn data in TextFile(s)

=head1 SYNOPSIS

AnalyzeTextFilesData.pl TextFile(s)...

AnalyzeTextFilesData.pl [B<-c, --colmode> colnum | collabel] [B<--columns> "colnum,[colnum,...]" | "collabel,[collabel,...]" | All]
[B<--columnpairs> "colnum,colnum,[colnum,colnum]..." | "collabel,collabel,[collabel,collabel]..." | AllPairs]
[B<-d, --detail> infolevel] [B<-f, --fast>] [B<--frequencybins> number | "number,number,[number,...]"] [B<-h, --help>]
[B<--indelim> comma | semicolon] [B<--klargest> number] [B<--ksmallest> number]
[B<-m, --mode> DescriptiveStatisticsBasic | DescriptiveStatisticsAll | All | "function1, [function2,...]"]
[B<-o, --overwrite>] [B<--outdelim> comma | tab | semicolon] [B<-p, --precision> number]
[B<-q, --quote> yes | no] [B<-r, --root> rootname] [B<--trimfraction> number] [B<-w, --workingdir> dirname] TextFiles(s)...

=head1 DESCRIPTION

Anaylze numerical column data in I<TextFile(s)> using a combination of various statistical
functions; Non-numerical values are simply ignored. For I<Correlation, RSquare, and Covariance>
analysis, the count of valid values in specifed column pair must be same; otherwise, column
pair is ignored. The file names are separated by space. The valid file extensions are I<.csv>
and I<.tsv> for comma/semicolon and tab delimited text files respectively. All other
file names are ignored. All the text files in a current directory can be specified by
I<*.csv>, I<*.tsv>, or the current directory name. The B<--indelim> option determines
the format of I<TextFile(s)>. Any file which doesn't correspond to the format indicated
by B<--indelim> option is ignored.

=head1 OPTIONS

=over 4

=item B<-c, --colmode> I<colnum | collabel>

Specify how columns are identified in TextFile(s): using column number or column
label. Possible values: I<colnum or collabel>. Default value: I<colnum>.

=item B<--columns> I<"colnum,[colnum,...]" | "collabel,[collabel]..." | All>

This value is mode specific. It's a list of comma delimited columns to use
for data analysis. Default value: I<First column>.

This value is ignored during I<Correlation/Pearson Correlation> and I<Covariance>
data analysis; B<-coulmnparis> option is used instead.

For I<colnum> value of B<-c, --colmode> option, input values format is:
I<colnum,colnum,...>. Example:

   1,3,5

For I<collabel> value of B<-c, --colmode> option, input values format is:
I<collabel,collabel,..>. Example:

    ALogP,MolWeight,EC50

=item B<--columnpairs> I<"colnum,colnum,[colnum,colnum,...]" | "collabel,collabel,[collabel,collabel,...]" | AllPairs>

This value is mode specific and is only used for I<Correlation, PearsonCorrelation, or
Covariance> value of B<-m, --mode> option. It is a comma delimited list of column pairs
to use for data analysis during I<Correlation> and I<Covariance> calculations. Default value:
I<First column, Second column>.

For I<colnum> value of B<-c, --colmode> option, input values format is:
I<colnum,colnum,[colnum,colnum]...>. Example:

    1,3,5,6,1,6

For I<collabel> value of B<-c, --colmode> option, input values format is:
I<collabel,collabel,[collabel,collabel]..>. Example:

    MolWeight,EC50,NumN+O,PSA

For I<AllPairs> value of B<--columnparis> option, all column pairs are used for I<Correlation>
and I<Covariance> calculations.

=item B<-d, --detail> I<infolevel>

Level of information to print about column values being ignored. Default: I<1>. Possible values:
1, 2, 3, or 4.

=item B<-f, --fast>

In this mode, all the columns specified for analysis are assumed to contain numerical
data and no checking is performed before analysis. By default, only numerical data is
used for analysis; other types of column data is ignored.

=item B<--frequencybins> I<number | "number,number,[number,...]">

Specify number of bins or bin range to use for frequency analysis. Default value: I<10>

Number of bins value along with the smallest and largest value for a column is used to
group the column values into different groups.

The bin range list is used to group values for a column into different groups; It must contain
values in ascending order. Examples:

    10,20,30
    0.1,0.2,0.3,0.4,0.5

The frequency value calculated for a specific bin corresponds to all the column values
which are greater than the previous bin value and less than or equal to the current bin value.

=item B<-h, --help>

Print this help message.

=item B<--indelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile(s)>. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a
delimiter.

=item B<--klargest> I<number>

Kth largest value to find by I<KLargest> function. Default value: I<2> Valid values: positive
integers.

=item B<--ksmallest> I<number>

Kth smallest value to find by I<KSmallest> function. Default value: I<2>. Valid values: positive
integers.

=item B<-m, --mode> I<DescriptiveStatisticsBasic | DescriptiveStatisticsAll | All | "function1, [function2,...]">

Specify how to analyze data in TextFile(s): calculate basic or all descriptive statistics; or
use a comma delimited list of supported statistical functions. Possible values:
I<DescriptiveStatisticsBasic | DescriptiveStatisticsAll | "function1,[function2]...">. Default
value: I<DescriptiveStatisticsBasic>

I<DescriptiveStatisticsBasic> includes these functions: I<Count, Maximum, Minimum, Mean,
Median, Sum, StandardDeviation, StandardError, Variance>.

I<DescriptiveStatisticsAll>, in addition to  I<DescriptiveStatisticsBasic> functions, includes:
I<GeometricMean, Frequency, HarmonicMean, KLargest, KSmallest, Kurtosis, Mode, RSquare,
Skewness, TrimMean>.

I<All> uses complete list of supported functions: I<Average, AverageDeviation, Correlation,
Count, Covariance, GeometricMean, Frequency, HarmonicMean, KLargest, KSmallest, Kurtosis,
Maximum, Minimum, Mean, Median, Mode, RSquare, Skewness, Sum,
SumOfSquares, StandardDeviation, StandardDeviationN, StandardError, StandardScores,
StandardScoresN, TrimMean, Variance, VarianceN>. The function names ending with N
calculate corresponding values assuming an entire population instead of a population sample.

Here are the formulas for these functions:

Average: See Mean

AverageDeviation: SUM( ABS(x[i] - Xmean) ) / n

Correlation: See Pearson Correlation

Covariance: SUM( (x[i] - Xmean)(y[i] - Ymean) ) / n

GeometricMean: NthROOT( PRODUCT(x[i]) )

HarmonicMean: 1 / ( SUM(1/x[i]) / n )

Mean: SUM( x[i] ) / n

Median: Xsorted[(n - 1)/2 + 1] for even values of n; (Xsorted[n/2] + Xsorted[n/2 + 1])/2
for odd values of n.

Kurtosis: [ {n(n + 1)/(n - 1)(n - 2)(n - 3)}  SUM{ ((x[i] - Xmean)/STDDEV)^4 } ] -
{3((n - 1)^2)}/{(n - 2)(n-3)}

PearsonCorrelation: SUM( (x[i] - Xmean)(y[i] - Ymean) ) / SQRT( SUM( (x[i] - Xmean)^2 )
(SUM( (y[i] - Ymean)^2 ))   )

RSquare: PearsonCorrelation^2

Skewness: {n/(n - 1)(n - 2)} SUM{ ((x[i] - Xmean)/STDDEV)^3 }

StandardDeviation: SQRT ( SUM( (x[i] - Mean)^2 ) / (n - 1) )

StandardDeviationN: SQRT ( SUM( (x[i] - Mean)^2 ) / n )

StandardError: StandardDeviation / SQRT( n )

StandardScore: (x[i] - Mean) / (n - 1)

StandardScoreN: (x[i] - Mean) / n

Variance: SUM( (x[i] - Xmean)^2  / (n - 1) )

VarianceN: SUM( (x[i] - Xmean)^2  / n )

=item B<-o, --overwrite>

Overwrite existing files.

=item B<--outdelim> I<comma | tab | semicolon>

Output text file delimiter. Possible values: I<comma, tab, or semicolon>
Default value: I<comma>.

=item B<-p, --precision> I<number>

Precision of calculated values in the output file. Default: up to I<2> decimal places.
Valid values: positive integers.

=item B<-q, --quote> I<yes | no>

Put quotes around column values in output text file. Possible values: I<yes or
no>. Default value: I<yes>.

=item B<-r, --root> I<rootname>

New text file name is generated using the root: <Root>.<Ext>. Default new file
name: <InitialTextFileName><Mode>.<Ext>. Based on the specified analysis,
<Mode> corresponds to one of these values: DescriptiveStatisticsBasic,
DescriptiveStatisticsAll, AllStatistics, SpecifiedStatistics, Covariance, Correlation,
Frequency, or StandardScores. The csv, and tsv <Ext> values are used for
comma/semicolon, and tab delimited text files respectively. This option is ignored for
multiple input files.

=item B<--trimfraction> I<number>

Fraction of data to exclude from the top and bottom of the data set during
I<TrimMean> calculation. Default value: I<0.1>. Valid values: > 0 and < 1.

=item B<-w --workingdir> I<text>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To calculate basic statistics for data in first column and generate a
NewSample1DescriptiveStatisticsBasic.csv file, type:

    % AnalyzeTextFilesData.pl -o -r NewSample1 Sample1.csv

To calculate basic statistics for data in third column and generate a
NewSample1DescriptiveStatisticsBasic.csv file, type:

    % AnalyzeTextFilesData.pl --columns 3 -o -r NewSample1 Sample1.csv

To calculate basic statistics for data in MolWeight column and generate a
NewSample1DescriptiveStatisticsBasic.csv file, type:

    % AnalyzeTextFilesData.pl -colmode collabel --columns MolWeight -o
    -r NewSample1 Sample1.csv

To calculate all available statistics for data in third column and all column pairs,
and generate NewSample1DescriptiveStatisticsAll.csv, NewSample1CorrelationMatrix.csv,
NewSample1CorrelationMatrix.csv, and NewSample1MolWeightFrequencyAnalysis.csv files,
type:

    % AnalyzeTextFilesData.pl -m DescriptiveStatisticsAll --columns 3 -o
    --columnpairs AllPairs -r NewSample1 Sample1.csv

To compute frequency distribution of data in third column into five bins and
generate NewSample1MolWeightFrequencyAnalysis.csv, type:

    % AnalyzeTextFilesData.pl -m Frequency --frequencybins 5 --columns 3
    -o -r NewSample1 Sample1.csv

To compute frequency distribution of data in third column into specified bin range
values, and generate NewSample1MolWeightFrequencyAnalysis.csv, type:

    % AnalyzeTextFilesData.pl -m Frequency --frequencybins "100,200,400"
    --columns 3 -o -r NewSample1 Sample1.csv

To calculate all available statistics for data in all columns and column pairs, type:

    % AnalyzeTextFilesData.pl -m All --columns  All --columnpairs
    AllPairs -o -r NewSample1 Sample1.csv

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

JoinTextFiles.pl, MergeTextFilesWithSD.pl, ModifyTextFilesFormat.pl, SplitTextFiles.pl, TextFilesToHTML.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
