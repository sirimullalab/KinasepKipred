#!/usr/bin/perl -w
#
# File: AnalyzeSDFilesData.pl
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
use SDFileUtil;
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

my(@SDFilesList);
@SDFilesList = ExpandFileNames(\@ARGV, "sd sdf");

print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Collect information about SD files...
print "Checking input SD file(s)...\n";
my(%SDFilesInfo);
RetrieveSDFilesInfo();
ProcessSDFilesDataLabelsInfo();

# Generate output files...
my($FileIndex);
if (@SDFilesList > 1) {
  print "\nProcessing SD files...\n";
}
for $FileIndex (0 .. $#SDFilesList) {
  if ($SDFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $SDFilesList[$FileIndex]...\n";
    AnalyzeSDFile($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Analyze data...
sub AnalyzeSDFile {
  my($Index) = @_;
  my($SDFile, $DataLabel, $DataValue, @DataLabelsToAnalyze, %DataFieldValuesToAnalyzeMap);

  $SDFile = $SDFilesList[$Index];
  @DataLabelsToAnalyze = @{$SDFilesInfo{UniqueDataLabelsToAnalyze}[$Index]};
  %DataFieldValuesToAnalyzeMap = ();
  for $DataLabel (@DataLabelsToAnalyze) {
    @{$DataFieldValuesToAnalyzeMap{$DataLabel}} = ();
  }

  # Collect appropriate data field label values for analysis...
  my($CmpdString, @CmpdLines, %DataFieldValues, $CmpdCount, $InvalidCmpdCount, @InvalidCmpdDataLabels);
  open SDFILE, "$SDFile" or die "Error: Can't open $SDFile: $! \n";
  $CmpdCount = 0;
  $InvalidCmpdCount = 0;
  while ($CmpdString = ReadCmpdString(\*SDFILE)) {
    $CmpdCount++;
    @CmpdLines = split "\n", $CmpdString;
    %DataFieldValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
    @InvalidCmpdDataLabels = ();
    DATALABEL: for $DataLabel (@DataLabelsToAnalyze) {
      if (exists $DataFieldValues{$DataLabel}) {
	$DataValue = $DataFieldValues{$DataLabel};
	if ($OptionsInfo{CheckData}) {
	  if (!IsNumerical($DataValue)) {
	    push @InvalidCmpdDataLabels, $DataLabel;
	    next DATALABEL;
	  }
	}
	push @{$DataFieldValuesToAnalyzeMap{$DataLabel}}, $DataValue;
      }
    }
    if (@InvalidCmpdDataLabels) {
      $InvalidCmpdCount++;
      if ($OptionsInfo{DetailLevel} >=4 ) {
	print "Compound record $CmpdCount contains ", scalar(@InvalidCmpdDataLabels)," non-numerical or empty value(s) for data field(s) - ", JoinWords(\@InvalidCmpdDataLabels, ", ", 0)," - to be analyzed:\n$CmpdString \n";
      }
      elsif ($OptionsInfo{DetailLevel} >= 3) {
	print "Compound record $CmpdCount contains ", scalar(@InvalidCmpdDataLabels)," non-numerical or empty value(s) for data field(s) - ", JoinWords(\@InvalidCmpdDataLabels, ", ", 0)," - to be analyzed...\n";
      }
      elsif ($OptionsInfo{DetailLevel} >= 2) {
	print "Compound record $CmpdCount contains ", scalar(@InvalidCmpdDataLabels)," non-numerical or empty value(s) for data field to be analyzed...\n";
      }
    }
  }
  if ($InvalidCmpdCount && ($OptionsInfo{DetailLevel} >= 1)) {
    print "Non-numerical or empty data present in $InvalidCmpdCount compound record(s)...\n";
  }
  close SDFILE;

  # Perform the analysis...
  my(@SpecifiedFunctionNames, $SpecifiedFunction);
  @SpecifiedFunctionNames = ();

  for $SpecifiedFunction (@{$OptionsInfo{SpecifiedStatisticalFunctions}}) {
    if ($SpecifiedFunction !~ /^(Covariance|Correlation|Frequency|Rsquare|StandardScores|StandardScoresN)$/i) {
      push @SpecifiedFunctionNames, $OptionsInfo{SpecifiedStatisticalFunctionsMap}{lc($SpecifiedFunction)};
    }
  }
  if (@SpecifiedFunctionNames) {
    PerformAnalysis($Index, \@SpecifiedFunctionNames, \%DataFieldValuesToAnalyzeMap)
  }
  if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{covariance}) || exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{correlation}) || exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{rsquare})) {
    if ($OptionsInfo{AllDataLabelPairs} || $OptionsInfo{CommonDataLabelPairs}) {
      PerformMatrixAnalysis($Index, \%DataFieldValuesToAnalyzeMap);
    }
    else {
      # Perform pairwise analysis for specified columns and write out calculated values - correlation
      # rsquare, or covariance - in the same file.
      PerformDataLabelPairAnalysis($Index, \%DataFieldValuesToAnalyzeMap);
    }
  }
  if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{standardscores}) || exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{standardscoresn}) ) {
    PerformStandardScoresAnalysis($Index, \%DataFieldValuesToAnalyzeMap);
  }
  if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{frequency})) {
    PerformFrequencyAnalysis($Index, \%DataFieldValuesToAnalyzeMap);
  }

}

# Calculate values for various statistical functions...
sub PerformAnalysis {
  my($Index, $SpecifiedFunctionNamesRef, $DataValuesToAnalyzeMapRef) = @_;
  my($NewTextFile, $Line, $SpecifiedFunction, $Label, @ColLabels, @DataLabelsToAnalyze);

  $NewTextFile = $SDFilesInfo{NewTextFileRoot}[$Index] . $OptionsInfo{FileNameMode} . "." . $SDFilesInfo{NewTextFileExt}[$Index];

  print "Generating new text file $NewTextFile...\n";
  open NEWTEXTFILE, ">$NewTextFile" or die "Error: Can't open $NewTextFile: $! \n";

  # Write out column labels...
  @ColLabels = ();
  push @ColLabels, "DataLabel";
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
  @DataLabelsToAnalyze = @{$SDFilesInfo{DataLabelsToAnalyze}[$Index]};

  # Turn off "strict"; otherwise, invoking statistical functions using function name string
  # is problematic.
  no strict;

  my($DataValuesRef, $DataLabel, $Value, @RowValues, %CalculatedValues);
  %CalculatedValues = ();
  for $DataLabel (@DataLabelsToAnalyze) {
    @RowValues = ();
    # Setup column id...
    push @RowValues, $DataLabel;
    $DataValuesRef =  \@{$DataValuesToAnalyzeMapRef->{$DataLabel}};
    FUNCTIONNAME: for $SpecifiedFunction (@{$SpecifiedFunctionNamesRef}) {
      $Value = "";
      if (!@{$DataValuesToAnalyzeMapRef->{$DataLabel}}) {
	# Invalid column values...
	push @RowValues, $Value;
	next FUNCTIONNAME;
      }
      if ($SpecifiedFunction =~ /^Count$/i) {
	$Value = @{$DataValuesToAnalyzeMapRef->{$DataLabel}};
      }
      elsif ($SpecifiedFunction =~ /^KLargest$/i) {
	$Value = &$SpecifiedFunction($DataValuesRef, $OptionsInfo{KLargest});
      }
      elsif ($SpecifiedFunction =~ /^KSmallest$/i) {
	$Value = &$SpecifiedFunction($DataValuesRef, $OptionsInfo{KSmallest});
      }
      elsif ($SpecifiedFunction =~ /^StandardDeviation$/i) {
	if (exists($CalculatedValues{$DataLabel}{StandardDeviation})) {
	  $Value = $CalculatedValues{$DataLabel}{StandardDeviation};
	}
	else {
	  $Value = &$SpecifiedFunction($DataValuesRef);
	  $CalculatedValues{$DataLabel}{StandardDeviation} = $Value;
	}
      }
      elsif ($SpecifiedFunction =~ /^StandardError$/i) {
	if (!exists($CalculatedValues{$DataLabel}{StandardDeviation})) {
	  $Value = StandardDeviation($DataValuesRef);
	  $CalculatedValues{$DataLabel}{StandardDeviation} = $Value;
	}
	if (defined $CalculatedValues{$DataLabel}{StandardDeviation}) {
	  $Value = &$SpecifiedFunction($CalculatedValues{$DataLabel}{StandardDeviation}, @{$DataValuesToAnalyzeMapRef->{$DataLabel}});
	}
      }
      elsif ($SpecifiedFunction =~ /^TrimMean$/i) {
	$Value = &$SpecifiedFunction($DataValuesRef, $OptionsInfo{TrimFraction});
      }
      else {
	$Value = &$SpecifiedFunction($DataValuesRef);
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

# Calculate covariance, correlation, rsquare for specified data field label pairs....
sub PerformDataLabelPairAnalysis {
  my($Index, $DataValuesToAnalyzeMapRef) = @_;
  my($NewTextFile, @ColLabels, $Line, $CalculateCorrelation, $CalculateRSquare, $CalculateCovariance);

  $CalculateCorrelation = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{correlation}) ? 1 : 0;
  $CalculateRSquare = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{rsquare}) ? 1 : 0;
  $CalculateCovariance = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{covariance}) ? 1 : 0;

  $NewTextFile = $SDFilesInfo{NewTextFileRoot}[$Index] . "DataFieldPairsAnalysis." .  $SDFilesInfo{NewTextFileExt}[$Index];
  print "Generating new text file $NewTextFile...\n";
  open NEWTEXTFILE, ">$NewTextFile" or die "Error: Can't open $NewTextFile: $! \n";

  # Write out the column labels...
  @ColLabels = ();
  push @ColLabels, ("DataLabel1", "DataLabel2");
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

  # Go over each data field pair...
  my($CorrelationValue, $RSquareValue, $CovarianceValue,  $LabelIndex, $DataLabel1, $DataLabel2, $DataValues1, $DataValues2, @DataLabelPairs1ToAnalyze, @DataLabelPairs2ToAnalyze, @RowValues, $Value);

  @DataLabelPairs1ToAnalyze = @{$SDFilesInfo{DataLabelPairs1ToAnalyze}[$Index]};
  @DataLabelPairs2ToAnalyze = @{$SDFilesInfo{DataLabelPairs2ToAnalyze}[$Index]};
  for $LabelIndex (0 .. $#DataLabelPairs1ToAnalyze) {
    @RowValues = ();
    $DataLabel1 = $DataLabelPairs1ToAnalyze[$LabelIndex];
    $DataLabel2 = $DataLabelPairs2ToAnalyze[$LabelIndex];
    $DataValues1 =  \@{$DataValuesToAnalyzeMapRef->{$DataLabel1}};
    $DataValues2 =  \@{$DataValuesToAnalyzeMapRef->{$DataLabel2}};

    # Setup column ids...
    push @RowValues, $DataLabel1;
    push @RowValues, $DataLabel2;

    if (@$DataValues1 != @$DataValues2) {
      # Print a warning...
      warn "Warning: Skipping analysis for data field pair $DataLabel1, $DataLabel2: Number of valid data values must be same.\n";
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
	$CorrelationValue = Correlation($DataValues1, $DataValues2);
	$Value = (defined($CorrelationValue) && length($CorrelationValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $CorrelationValue) + 0) : "";
	push @RowValues, $Value;
	if ($CalculateRSquare) {
	  $RSquareValue = (defined($CorrelationValue) && length($CorrelationValue)) ? ($CorrelationValue ** 2) : "";
	  $Value = (length($RSquareValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $RSquareValue) + 0) : "";
	  push @RowValues, $Value;
	}
      }
      if ($CalculateCovariance) {
	$CovarianceValue = Covariance($DataValues1, $DataValues2);
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
  my($Index, $DataValuesToAnalyzeMapRef) = @_;
  my($NewTextFile, $ColLabel, @ColLabels, @RowValues, $Line, $DataLabel, @DataLabelsToAnalyze, $DataValuesRef, $BinValue, $FrequencyValue, $Value, %FrequencyMap);

  @DataLabelsToAnalyze = @{$SDFilesInfo{DataLabelsToAnalyze}[$Index]};
  for $DataLabel (@DataLabelsToAnalyze) {
    $NewTextFile = $SDFilesInfo{NewTextFileRoot}[$Index] . $DataLabel . "FrequencyAnalysis." .  $SDFilesInfo{NewTextFileExt}[$Index];
    print "Generating new text file $NewTextFile...\n";
    open NEWTEXTFILE, ">$NewTextFile" or die "Error: Can't open $NewTextFile: $! \n";

    # Write out the column labels...
    @ColLabels = ();
    push @ColLabels , ("Bins", "Frequency");
    $Line = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print NEWTEXTFILE "$Line\n";

    #Calculate and write out frequency values...
    %FrequencyMap = ();
    $DataValuesRef =  \@{$DataValuesToAnalyzeMapRef->{$DataLabel}};
    if (@$DataValuesRef) {
      if (@{$OptionsInfo{BinRange}}) {
	%FrequencyMap = Frequency($DataValuesRef, \@{$OptionsInfo{BinRange}});
      }
      else {
	%FrequencyMap = Frequency($DataValuesRef, $OptionsInfo{NumOfBins});
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
  my($Index, $DataValuesToAnalyzeMapRef) = @_;
  my($CorrelationTextFile, $CovarianceTextFile, $RSquareTextFile, $CalculateCorrelation, $CalculateRSquare, $CalculateCovariance);

  $CalculateCorrelation = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{correlation}) ? 1 : 0;
  $CalculateRSquare = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{rsquare}) ? 1 : 0;
  $CalculateCovariance = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{covariance}) ? 1 : 0;

  $CorrelationTextFile = $SDFilesInfo{NewTextFileRoot}[$Index] . "CorrelationMatrix." .  $SDFilesInfo{NewTextFileExt}[$Index];
  $RSquareTextFile = $SDFilesInfo{NewTextFileRoot}[$Index] . "RSquareMatrix." .  $SDFilesInfo{NewTextFileExt}[$Index];
  $CovarianceTextFile = $SDFilesInfo{NewTextFileRoot}[$Index] . "CovarianceMatrix." .  $SDFilesInfo{NewTextFileExt}[$Index];

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

  my($Line, $Value, $CorrelationValue, $RSquareValue, $CovarianceValue, $DataLabel, $DataLabel1, $DataLabel2, $DataValuesRef1, $DataValuesRef2, @ColLabels, @CovarianceRowValues, @CorrelationRowValues, @RSquareRowValues);

  # Write out the column labels...
  @ColLabels = ();
  push @ColLabels, @{$SDFilesInfo{AllDataLabels}[$Index]};
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
  my(%CorrelationMatrixMap, %RSquareMatrixMap, %CovarianceMatrixMap, $LabelIndex1, $LabelIndex2, @DataLabelsToAnalyze);

  %CorrelationMatrixMap = (); %RSquareMatrixMap = (); %CovarianceMatrixMap = ();
  @DataLabelsToAnalyze = ();
  @DataLabelsToAnalyze = $OptionsInfo{AllDataLabelPairs} ? @{$SDFilesInfo{AllDataLabels}[$Index]} : @{$SDFilesInfo{CommonDataLabels}[$Index]};

  for $LabelIndex1 (0 .. (@DataLabelsToAnalyze - 1)) {
    $DataLabel1 = $DataLabelsToAnalyze[$LabelIndex1];
    for $LabelIndex2 (0 .. $LabelIndex1) {
      $DataLabel2 = $DataLabelsToAnalyze[$LabelIndex2];
      $DataValuesRef1 =  \@{$DataValuesToAnalyzeMapRef->{$DataLabel1}};
      $DataValuesRef2 =  \@{$DataValuesToAnalyzeMapRef->{$DataLabel2}};
      if ($CalculateCorrelation || $CalculateRSquare) {
	$CorrelationValue = Correlation($DataValuesRef1, $DataValuesRef2);
	$CorrelationValue = (defined($CorrelationValue) && length($CorrelationValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $CorrelationValue) + 0) : "";
	$CorrelationMatrixMap{$DataLabel1}{$DataLabel2} = $CorrelationValue;
	if ($DataLabel1 ne $DataLabel2) {
	  $CorrelationMatrixMap{$DataLabel2}{$DataLabel1} = $CorrelationValue;
	}
	if ($CalculateRSquare) {
	  $RSquareValue = (defined($CorrelationValue) && length($CorrelationValue)) ? ($CorrelationValue ** 2) : "";
	  $RSquareValue = (length($RSquareValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $RSquareValue) + 0) : "";
	  $RSquareMatrixMap{$DataLabel1}{$DataLabel2} = $RSquareValue;
	  if ($DataLabel1 ne $DataLabel2) {
	    $RSquareMatrixMap{$DataLabel2}{$DataLabel1} = $RSquareValue;
	  }
	}
      }
      if ($CalculateCovariance) {
	$CovarianceValue = Covariance($DataValuesRef1, $DataValuesRef2);
	$CovarianceValue = (defined($CovarianceValue) && length($CovarianceValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $CovarianceValue) + 0) : "";
	$CovarianceMatrixMap{$DataLabel1}{$DataLabel2} = $CovarianceValue;
	if ($DataLabel1 ne $DataLabel2) {
	  $CovarianceMatrixMap{$DataLabel2}{$DataLabel1} = $CovarianceValue;
	}
      }
    }
  }

  # Write out the matrices...
  for $LabelIndex1 (0 .. (@DataLabelsToAnalyze - 1)) {
    $DataLabel1 = $DataLabelsToAnalyze[$LabelIndex1];
    @CorrelationRowValues = ();
    @RSquareRowValues = ();
    @CovarianceRowValues = ();
    if ($CalculateCorrelation || $CalculateRSquare) {
      push @CorrelationRowValues, $DataLabel1;
      if ($CalculateRSquare) {
	push @RSquareRowValues, $DataLabel1;
      }
    }
    if ($CalculateCovariance) {
      push @CovarianceRowValues, $DataLabel;
    }
    for $LabelIndex2 (0 .. (@DataLabelsToAnalyze - 1)) {
      $DataLabel2 = $DataLabelsToAnalyze[$LabelIndex2];
      if ($CalculateCorrelation || $CalculateRSquare) {
	push @CorrelationRowValues, $CorrelationMatrixMap{$DataLabel1}{$DataLabel2};
	if ($CalculateRSquare) {
	  push @RSquareRowValues, $RSquareMatrixMap{$DataLabel1}{$DataLabel2};
	}
      }
      if ($CalculateCovariance) {
	push @CovarianceRowValues, $CovarianceMatrixMap{$DataLabel1}{$DataLabel2};
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
  my($Index, $DataValuesToAnalyzeMapRef) = @_;
  my($StandardScores, $StandardScoresN, $NewTextFile, @ColLabels, $Label, $NewLine);

  $StandardScores = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{standardscores}) ? 1 : 0;
  $StandardScoresN = exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{standardscoresn}) ? 1 : 0;

  $NewTextFile = $SDFilesInfo{NewTextFileRoot}[$Index] . "StandardScores." .  $SDFilesInfo{NewTextFileExt}[$Index];
  print "Generating new text file $NewTextFile...\n";
  open NEWTEXTFILE, ">$NewTextFile" or die "Error: Can't open $NewTextFile: $! \n";

  my($DataLabel, @DataLabelsToAnalyze);
  # Write out column labels...
  @ColLabels = ();
  @DataLabelsToAnalyze = @{$SDFilesInfo{DataLabelsToAnalyze}[$Index]};
  for $DataLabel (@DataLabelsToAnalyze) {
    if ($StandardScores) {
      push @ColLabels, "${DataLabel}\(StandardScores)";
    }
    if ($StandardScoresN) {
      push @ColLabels, "${DataLabel}\(StandardScoresN)";
    }
  }
  $NewLine = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print NEWTEXTFILE "$NewLine\n";

  # Go over each column to be analyzed and calculate standard deviation
  # and mean values...
  my($DataValuesRef, %StandardDeviationMap, %StandardDeviationNMap, %MeanMap);
  %StandardDeviationMap = ();
  %StandardDeviationNMap = ();
  %MeanMap = ();
  for $DataLabel (@DataLabelsToAnalyze) {
    $DataValuesRef =  \@{$DataValuesToAnalyzeMapRef->{$DataLabel}};
    if (!exists($MeanMap{$DataLabel})) {
      $MeanMap{$DataLabel} = Mean($DataValuesRef);
    }
    if ($StandardScores) {
      if (!exists($StandardDeviationMap{$DataLabel})) {
	$StandardDeviationMap{$DataLabel} = StandardDeviation($DataValuesRef);
      }
    }
    if ($StandardScoresN) {
      if (!exists($StandardDeviationNMap{$DataLabel})) {
	$StandardDeviationNMap{$DataLabel} = StandardDeviationN($DataValuesRef);
      }
    }
  }
  #
  # Go over each data field and calculate standard scores for each column
  # using (x[i] - mean) / (n - 1) for StandardScores and (x[i] - mean) / n
  # for StandardScoresN; write out the calculated values as well...

  my($SDFile, $Value, $ValueOkay, $ScoreValue, @RowValues, $CmpdString, @CmpdLines, %DataFieldValues);
  $SDFile = $SDFilesList[$Index];

  open SDFILE, "$SDFile" or die "Error: Can't open $SDFile: $! \n";
  while ($CmpdString = ReadCmpdString(\*SDFILE)) {
    @CmpdLines = split "\n", $CmpdString;
    %DataFieldValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
    @RowValues = ();
    for $DataLabel (@DataLabelsToAnalyze) {
      $Value = "";
      if (exists $DataFieldValues{$DataLabel}) {
	$Value = $DataFieldValues{$DataLabel};
      }
      $ValueOkay = ($OptionsInfo{CheckData} && !IsNumerical($Value)) ? 0 : 1;
      if ($StandardScores) {
	$ScoreValue = $ValueOkay ? (($Value - $MeanMap{$DataLabel})/$StandardDeviationMap{$DataLabel}) : "";
	$ScoreValue = (defined($ScoreValue) && length($ScoreValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $ScoreValue) + 0) : "";
	push @RowValues, $ScoreValue;
      }
      if ($StandardScoresN) {
	$ScoreValue = $ValueOkay ? (($Value - $MeanMap{$DataLabel})/$StandardDeviationNMap{$DataLabel}) : "";
	$ScoreValue = (defined($ScoreValue) && length($ScoreValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $ScoreValue) + 0) : "";
	push @RowValues, $ScoreValue;
      }
    }
    $NewLine = JoinWords(\@RowValues, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print NEWTEXTFILE "$NewLine\n";
  }
  close SDFILE;
  close NEWTEXTFILE;

}

# Make sure the specified data field labels exists in SD files...
sub ProcessSDFilesDataLabelsInfo {
  my($Index, $DataFieldIndex, $SDFile, $DataLabel, @DataLabelsToAnalyze, %UniqueDataLabelsToAnalyzeMap);

  @{$SDFilesInfo{DataLabelsToAnalyze}} = ();
  @{$SDFilesInfo{DataLabelPairs1ToAnalyze}} = ();
  @{$SDFilesInfo{DataLabelPairs2ToAnalyze}} = ();
  @{$SDFilesInfo{UniqueDataLabelsToAnalyze}} = ();

  FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFile = $SDFilesList[$Index];

    @{$SDFilesInfo{DataLabelsToAnalyze}[$Index]} = ();
    @{$SDFilesInfo{DataLabelPairs1ToAnalyze}[$Index]} = ();
    @{$SDFilesInfo{DataLabelPairs2ToAnalyze}[$Index]} = ();
    @{$SDFilesInfo{UniqueDataLabelsToAnalyze}[$Index]} = ();

    %UniqueDataLabelsToAnalyzeMap = ();

    if ($SDFilesInfo{FileOkay}[$Index]) {
      @DataLabelsToAnalyze = ();
      if (@{$OptionsInfo{SpecifiedDataLabels}}) {
	for $DataLabel (@{$OptionsInfo{SpecifiedDataLabels}}) {
	  if (exists($SDFilesInfo{AllDataLabelsMap}[$Index]{$DataLabel})) {
	    push @DataLabelsToAnalyze, $DataLabel;
	  }
	}
      }
      elsif (defined($OptionsInfo{DataFields}) && $OptionsInfo{DataFields} =~ /^All$/i) {
	push @DataLabelsToAnalyze, @{$SDFilesInfo{AllDataLabels}[$Index]};
      }
      else {
	push @DataLabelsToAnalyze, @{$SDFilesInfo{CommonDataLabels}[$Index]};
      }
      if (@DataLabelsToAnalyze) {
	push @{$SDFilesInfo{DataLabelsToAnalyze}[$Index]}, @DataLabelsToAnalyze;
	# Set up unique data field label map as well...
	for $DataLabel (@DataLabelsToAnalyze) {
	  if (!exists $UniqueDataLabelsToAnalyzeMap{$DataLabel}) {
	    $UniqueDataLabelsToAnalyzeMap{$DataLabel} = $DataLabel;
	  }
	}
      }
      else {
	warn "Warning: Ignoring file $SDFile: None of the data field labels specified, @{$OptionsInfo{SpecifiedDataLabels}}, using \"--datafields\" option exist.\n";
	$SDFilesInfo{FileOkay}[$Index] = 0;
	next FILELIST;
      }
      if (!$OptionsInfo{Overwrite} && exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{frequency})) {
	# Make sure specific frequency files don't exist...
	my($FrequencyFile);
	for $DataLabel (@DataLabelsToAnalyze) {
	  $FrequencyFile = $SDFilesInfo{NewTextFileRoot}[$Index] . $SDFilesInfo{AllDataLabelsMap}[$Index]{$DataLabel} . "FrequencyAnalysis." .  $SDFilesInfo{NewTextFileExt}[$Index];
	  if (-e $FrequencyFile) {
	    warn "Warning: Ignoring file $SDFile: The file $FrequencyFile already exists.\n";
	    $SDFilesInfo{FileOkay}[$Index] = 0;
	    next FILELIST;
	  }
	}
      }
      # Setup specified data field label pairs...
      if (exists $OptionsInfo{SpecifiedStatisticalFunctionsMap}{correlation} || exists $OptionsInfo{SpecifiedStatisticalFunctionsMap}{covariance} || exists $OptionsInfo{SpecifiedStatisticalFunctionsMap}{rsquare}) {
	my(@DataLabelPairsToAnalyze, $DataLabel1, $DataLabel2);
	if (@{$OptionsInfo{SpecifiedDataLabelPairs}}) {
	  # Make sure both data field labels exist...
	  my($DataFieldIndex);
	  for ($DataFieldIndex = 0; (($DataFieldIndex + 1) < @{$OptionsInfo{SpecifiedDataLabelPairs}}); $DataFieldIndex += 2 ) {
	    $DataLabel1 = $OptionsInfo{SpecifiedDataLabelPairs}[$DataFieldIndex];
	    $DataLabel2 = $OptionsInfo{SpecifiedDataLabelPairs}[$DataFieldIndex + 1];
	    if (exists($SDFilesInfo{AllDataLabelsMap}[$Index]{$DataLabel1}) && exists($SDFilesInfo{AllDataLabelsMap}[$Index]{$DataLabel2})) {
	      push @DataLabelPairsToAnalyze, ($DataLabel1, $DataLabel2);
	    }
	  }
	}
	elsif ($OptionsInfo{AllDataLabelPairs}) {
	  for $DataLabel1 (@{$SDFilesInfo{AllDataLabels}[$Index]}) {
	    for $DataLabel2 (@{$SDFilesInfo{AllDataLabels}[$Index]}) {
	      push @DataLabelPairsToAnalyze, ($DataLabel1, $DataLabel2);
	    }
	  }
	}
	else {
	  for $DataLabel1 (@{$SDFilesInfo{CommonDataLabels}[$Index]}) {
	    for $DataLabel2 (@{$SDFilesInfo{CommonDataLabels}[$Index]}) {
	      push @DataLabelPairsToAnalyze, ($DataLabel1, $DataLabel2);
	    }
	  }
	}
	if (@DataLabelPairsToAnalyze) {
	  if (@DataLabelPairsToAnalyze % 2) {
	    warn "Warning: Ignoring file $SDFile: Invalid number  values specified using \"--datafieldpairs\" option: It must contain even number of valid values.\n";
	    $SDFilesInfo{FileOkay}[$Index] = 0;
	    next FILELIST;
	  }
	  else {
	    for ($DataFieldIndex = 0; $DataFieldIndex < @DataLabelPairsToAnalyze; $DataFieldIndex += 2) {
	      push @{$SDFilesInfo{DataLabelPairs1ToAnalyze}[$Index]}, $DataLabelPairsToAnalyze[$DataFieldIndex];
	      push @{$SDFilesInfo{DataLabelPairs2ToAnalyze}[$Index]}, $DataLabelPairsToAnalyze[$DataFieldIndex + 1];
	    }
	    # Set up unique data field labe map as well...
	    for $DataLabel (@DataLabelPairsToAnalyze) {
	      if (!exists $UniqueDataLabelsToAnalyzeMap{$DataLabel}) {
		$UniqueDataLabelsToAnalyzeMap{$DataLabel} = $DataLabel;
	      }
	    }
	  }
	}
      }
      # Setup unique data field label array...
      push @{$SDFilesInfo{UniqueDataLabelsToAnalyze}[$Index]}, (sort keys %UniqueDataLabelsToAnalyzeMap);
    }
  }
}

# Retrieve information about input SD files...
sub RetrieveSDFilesInfo {
  my($SDFile, $Index, $FileDir, $FileExt, $FileName, $OutFile, $OutFileRoot, $OutFileExt, $CmpdCount);

  %SDFilesInfo = ();

  @{$SDFilesInfo{FileOkay}} = ();
  @{$SDFilesInfo{CmpdCount}} = ();
  @{$SDFilesInfo{NewTextFileRoot}} = ();
  @{$SDFilesInfo{NewTextFileExt}} = ();

  @{$SDFilesInfo{AllDataFieldLabels}} = ();
  @{$SDFilesInfo{AllDataFieldLabelsMap}} = ();
  @{$SDFilesInfo{CommonDataLabels}} = ();

  FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFile = $SDFilesList[$Index];

    $SDFilesInfo{FileOkay}[$Index] = 0;

    $SDFilesInfo{CmpdCount}[$Index] = 0;
    $SDFilesInfo{NewTextFileRoot}[$Index] = "";
    $SDFilesInfo{NewTextFileExt}[$Index] = "";

    @{$SDFilesInfo{AllDataLabels}[$Index]} = ();
    %{$SDFilesInfo{AllDataLabelsMap}[$Index]} = ();
    @{$SDFilesInfo{CommonDataLabels}[$Index]} = ();

    if (!(-e $SDFile)) {
      warn "Warning: Ignoring file $SDFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($SDFile, "sd sdf")) {
      warn "Warning: Ignoring file $SDFile: It's not a SD file\n";
      next FILELIST;
    }

    # Generate appropriate name for the new text files...
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($SDFile);
    $OutFileExt = "csv";
    if ($Options{outdelim} =~ /^tab$/i) {
      $OutFileExt = "tsv";
    }
    if ($Options{root} && (@SDFilesList == 1)) {
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

    if (!$OptionsInfo{Overwrite}) {
      if (-e $OutFile) {
	warn "Warning: Ignoring file $SDFile: The file $OutFile already exists\n";
	next FILELIST;
      }
      if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{covariance}) || exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{correlation}) || exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{rsquare})) {
	if ($OptionsInfo{AllDataLabelPairs}) {
	  if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{covariance}) && (-e "${OutFileRoot}CovarianceMatrix.${FileExt}")) {
	    warn "Warning: Ignoring file $SDFile: The file ${OutFileRoot}Covariance.${FileExt} already exists.\n";
	    next FILELIST;
	  }
	  if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{correlation}) && (-e "${OutFileRoot}CorrelationMatrix.${FileExt}")) {
	    warn "Warning: Ignoring file $SDFile: The file ${OutFileRoot}CorrelationMatrix.${FileExt} already exists.\n";
	    next FILELIST;
	  }
	  if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{rsquare}) && (-e "${OutFileRoot}RSquareMatrix.${FileExt}")) {
	    warn "Warning: Ignoring file $SDFile: The file ${OutFileRoot}RSquareMatrix.${FileExt} already exists.\n";
	    next FILELIST;
	  }
	}
	else {
	  if (-e "${OutFileRoot}ColumnPairsAnalysis.${FileExt}") {
	    warn "Warning: Ignoring file $SDFile: The file ${OutFileRoot}ColumnPairsAnalysis.${FileExt} already exists.\n";
	    next FILELIST;
	  }
	}
      }
      if (exists($OptionsInfo{SpecifiedStatisticalFunctionsMap}{standardscores}) && (-e "${OutFileRoot}StandardScores.${FileExt}")) {
	warn "Warning: Ignoring file $SDFile: The file ${OutFileRoot}StandardScores.${FileExt} already exists.\n";
	next FILELIST;
      }
    }

    if (!open SDFILE, "$SDFile") {
      warn "Warning: Ignoring file $SDFile: Couldn't open it: $! \n";
      next FILELIST;
    }

    my($CmpdCount, $Label, $DataFieldLabelsRef, $CommonDataFieldLabelsRef, @DataFieldLabels, @CommonDataFieldLabels);
    $CmpdCount = 0;
    @DataFieldLabels = ();
    @CommonDataFieldLabels = ();
    ($CmpdCount, $DataFieldLabelsRef, $CommonDataFieldLabelsRef) = GetAllAndCommonCmpdDataHeaderLabels(\*SDFILE);
    push @DataFieldLabels, @{$DataFieldLabelsRef};
    push @CommonDataFieldLabels, @{$CommonDataFieldLabelsRef};
    close SDFILE;

    $SDFilesInfo{FileOkay}[$Index] = 1;
    $SDFilesInfo{NewTextFileRoot}[$Index] = "$OutFileRoot";
    $SDFilesInfo{NewTextFileExt}[$Index] = "$OutFileExt";

    $SDFilesInfo{CmpdCount}[$Index] = $CmpdCount;
    push @{$SDFilesInfo{AllDataLabels}[$Index]}, @DataFieldLabels;
    push @{$SDFilesInfo{CommonDataLabels}[$Index]}, @CommonDataFieldLabels;
    for $Label (@DataFieldLabels) {
      $SDFilesInfo{AllDataLabelsMap}[$Index]{$Label} = $Label;
    }
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{DataFields} = defined $Options{datafields} ? $Options{datafields} : undef;

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

  # Setup delimiter and quotes...
  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /tab/i ) ? "\t" : (($Options{outdelim} =~ /semicolon/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /yes/i ) ? 1 : 0;

  $OptionsInfo{Overwrite} = defined $Options{overwrite} ? $Options{overwrite} : undef;
  $OptionsInfo{Root} = defined $Options{root} ? $Options{root} : undef;

  # Setup miscellaneous options...
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

  # Setup specified data field labels...
  @{$OptionsInfo{SpecifiedDataLabels}} = ();
  if (defined $Options{datafields} && $Options{datafields} !~ /^(All|Common)$/i ) {
    my(@SpecifiedValues) = split ",", $Options{datafields};
    push @{$OptionsInfo{SpecifiedDataLabels}}, @SpecifiedValues;
  }
  @{$OptionsInfo{SpecifiedDataLabelPairs}} = ();
  $OptionsInfo{AllDataLabelPairs} = (defined($Options{datafieldpairs}) && $Options{datafieldpairs} =~ /^AllPairs$/i) ? 1 : 0;
  $OptionsInfo{CommonDataLabelPairs} = (defined($Options{datafieldpairs}) && $Options{datafieldpairs} =~ /^CommonPairs$/i) ? 1 : 0;
  if (defined($Options{datafieldpairs}) && !$OptionsInfo{AllDataLabelPairs} && !$OptionsInfo{CommonDataLabelPairs}) {
    my(@SpecifiedValues) = split ",", $Options{datafieldpairs};
    if (@SpecifiedValues % 2) {
      die "Error: Invalid number of values specified using \"--datafieldpairs\" option: It must contain even number of values.\n";
    }
    push @{$OptionsInfo{SpecifiedDataLabelPairs}}, @SpecifiedValues;
  }

}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{detail} = 0;
  $Options{datafields} = "Common";
  $Options{datafieldpairs} = "CommonPairs";
  $Options{frequencybins} = 10;
  $Options{klargest} = 2;
  $Options{ksmallest} = 2;
  $Options{mode} = "DescriptiveStatisticsBasic";
  $Options{outdelim} = "comma";
  $Options{precision} = 2;
  $Options{quote} = "yes";
  $Options{trimfraction} = 0.1;

  if (!GetOptions(\%Options, "datafields=s", "datafieldpairs=s", "detail|d=i", "frequencybins=s", "fast|f", "help|h", "klargest=i", "ksmallest=i", "mode|m=s", "outdelim=s", "overwrite|o", "precision|p=i", "quote|q=s", "root|r=s", "trimfraction=f", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if (!IsInteger($Options{detail})) {
    die "Error: The value specified, $Options{detail}, for option \"-d --detail\" is not valid. Allowed values: >= 0\n";
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

AnalyzeSDFilesData.pl - Analyze numerical data field values in SDFile(s)

=head1 SYNOPSIS

AnalyzeSDFilesData.pl SDFile(s)...

AnalyzeSDFilesData.pl [B<--datafields> "fieldlabel,[fieldlabel,...]" | All]
[B<--datafieldpairs> "fieldlabel,fieldlabel,[fieldlabel,fieldlabel,...]" | AllPairs] [B<-d, --detail> infolevel]
[B<-f, --fast>] [B<--frequencybins> number | "number,number,[number,...]"]
[B<-h, --help>] [B<--klargest> number] [B<--ksmallest> number]
[B<-m, --mode> DescriptiveStatisticsBasic | DescriptiveStatisticsAll | All | "function1, [function2,...]"]
[B<--trimfraction> number] [B<-w, --workingdir> dirname] SDFiles(s)...

=head1 DESCRIPTION

Analyze numerical data field values in I<SDFile(s)> using a combination of various statistical
functions; Non-numerical values are simply ignored. For I<Correlation, RSquare, and
Covariance> analysis, the count of valid values in specified data field pairs must be same;
otherwise, column data field pair is ignored. The file names are separated by space.The valid file
extensions are I<.sdf> and I<.sd>. All other file names are ignored. All the SD files in a
current directory can be specified either by I<*.sdf> or the current directory name.

=head1 OPTIONS

=over 4

=item B<--datafields> I<"fieldlabel,[fieldlabel,...]" | Common | All>

Data fields to use for analysis. Possible values: list of comma separated data field
labels, data fields common to all records, or all data fields. Default value: I<Common>.
Examples:

    ALogP,MolWeight,EC50
    "MolWeight,PSA"

=item B<--datafieldpairs> I<"fieldlabel,fieldlabel,[fieldlabel,fieldlabel,...]" | CommonPairs | AllPairs>

This value is mode specific and is only used for I<Correlation, PearsonCorrelation, or
Covariance> value of B<-m, --mode> option. It specifies data field label pairs to use
for data analysis during I<Correlation> and I<Covariance> calculations. Possible values:
comma delimited list of data field label pairs, data field label pairs common to all records,
or all data field pairs. Default value:I<CommonPairs>. Example:

    MolWeight,EC50,NumN+O,PSA

For I<AllPairs> value of B<--datafieldpairs> option, all data field label pairs are used for
I<Correlation> and I<Covariance> calculations.

=item B<-d, --detail> I<infolevel>

Level of information to print about column values being ignored. Default: I<0>. Possible values:
0, 1, 2, 3, or 4.

=item B<-f, --fast>

In this mode, all the data field values specified for analysis are assumed to contain numerical
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

=item B<--klargest> I<number>

Kth largest value to find by I<KLargest> function. Default value: I<2>. Valid values: positive
integers.

=item B<--ksmallest> I<number>

Kth smallest value to find by I<KSmallest> function. Default values: I<2>. Valid values: positive
integers.

=item B<-m, --mode> I<DescriptiveStatisticsBasic | DescriptiveStatisticsAll | All | "function1, [function2,...]">

Specify how to analyze data in SDFile(s): calculate basic or all descriptive statistics; or
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
name: <InitialSDFileName><Mode>.<Ext>. Based on the specified analysis,
<Mode> corresponds to one of these values: DescriptiveStatisticsBasic,
DescriptiveStatisticsAll, AllStatistics, SpecifiedStatistics, Covariance, Correlation,
Frequency, or StandardScores. The csv, and tsv <Ext> values are used for
comma/semicolon, and tab delimited text files respectively. This option is ignored for
multiple input files.

=item B<--trimfraction> I<number>

Fraction of data to exclude from the top and bottom of the data set during
I<TrimMean> calculation. Default value: I<0.1> Valid values: > 0 and < 1.

=item B<-w --workingdir> I<text>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To calculate basic statistics for data in all common data fields and generate a
NewSample1DescriptiveStatisticsBasic.csv file, type:

    % AnalyzeSDFilesData.pl -o -r NewSample1 Sample1.sdf

To calculate basic statistics for MolWeight data field and generate a
NewSample1DescriptiveStatisticsBasic.csv file, type:

    % AnalyzeSDFilesData.pl --datafields MolWeight -o -r NewSample1
    Sample1.sdf

To calculate all available statistics for MolWeight data field and all data field pairs,
and generate NewSample1DescriptiveStatisticsAll.csv, NewSample1CorrelationMatrix.csv,
NewSample1CorrelationMatrix.csv, and NewSample1MolWeightFrequencyAnalysis.csv
files, type:

    % AnalyzeSDFilesData.pl -m DescriptiveStatisticsAll --datafields
    MolWeight -o --datafieldpairs AllPairs -r NewSample1 Sample1.sdf

To compute frequency distribution of MolWeight data field into five bins and
generate NewSample1MolWeightFrequencyAnalysis.csv, type:

    % AnalyzeSDFilesData.pl -m Frequency --frequencybins 5 --datafields
    MolWeight -o -r NewSample1 Sample1.sdf

To compute frequency distribution of data in MolWeight data field into specified bin range
values, and generate NewSample1MolWeightFrequencyAnalysis.csv, type:

    % AnalyzeSDFilesData.pl -m Frequency --frequencybins "100,200,400"
    --datafields MolWeight -o -r NewSample1 Sample1.sdf

To calculate all available statistics for data in all data fields and pairs, type:

    % AnalyzeSDFilesData.pl -m All --datafields  All --datafieldpairs
    AllPairs -o -r NewSample1 Sample1.sdf

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
