#!/usr/bin/perl -w
#
# File: ElementalAnalysisTextFiles.pl
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
use MolecularFormula;

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

# Process options...
my(%OptionsInfo);
print "Processing options...\n";
ProcessOptions();

print "Checking input text file(s)...\n";
my(%TextFilesInfo);
RetrieveTextFilesInfo();
RetrieveColumnsAndLabelsInfo();

# Generate output files...
my($FileIndex);
if (@TextFilesList > 1) {
  print "\nProcessing text files...\n";
}
for $FileIndex (0 .. $#TextFilesList) {
  if ($TextFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $TextFilesList[$FileIndex]...\n";
    PerformElementalAnalysis($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Perform elemental analysis...
sub PerformElementalAnalysis {
  my($Index) = @_;
  my($TextFile, $NewTextFile, $FormulaCol, $Line, $NewLine, $FormulaColValue, $InDelim, $ColNum, $Value, $Status, $ErrorMsg, @ColLabels, @LineWords, @ColNumsBeforeNew, @ColNumsAfterNew);

  $TextFile = $TextFilesList[$Index];
  $InDelim = $TextFilesInfo{InDelim}[$Index];
  $NewTextFile = $TextFilesInfo{OutFile}[$Index];
  $FormulaCol = $TextFilesInfo{FormulaColNum}[$Index];

  @ColNumsBeforeNew = @{$TextFilesInfo{ColNumsBeforeNew}[$Index]};
  @ColNumsAfterNew = @{$TextFilesInfo{ColNumsAfterNew}[$Index]};

  print "Generating new Text file $NewTextFile...\n";
  open NEWTEXTFILE, ">$NewTextFile" or die "Error: Couldn't open $NewTextFile: $! \n";
  open TEXTFILE, "$TextFile" or die "Error: Can't open $TextFile: $! \n";

  # Skip over column labels from old file...
  $Line = GetTextLine(\*TEXTFILE);

  # Add column lablels in new file...
  @ColLabels = ();
  for $ColNum (@ColNumsBeforeNew) {
    push @ColLabels, $TextFilesInfo{ColLabels}[$Index][$ColNum];
  }
  for $Value (@{$OptionsInfo{SpecifiedCalculations}}) {
    push @ColLabels, $TextFilesInfo{ValueLabelsMap}[$Index]{$Value};
  }
  for $ColNum (@ColNumsAfterNew) {
    push @ColLabels, $TextFilesInfo{ColLabels}[$Index][$ColNum];
  }
  $NewLine = '';
  $NewLine = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print NEWTEXTFILE "$NewLine\n";

  # Go over all rows...
  my($LineCount, $ElementsRef, $ElementCompositionRef, $CalculationType, $CalculatedValue, @CalculatedValues);

  $LineCount = 1;
  TEXTLINE: while ($Line = GetTextLine(\*TEXTFILE)) {
    @LineWords = quotewords($InDelim, 0, $Line);
    $LineCount++;

    @CalculatedValues = ();
    for $Value (@{$OptionsInfo{SpecifiedCalculations}}) {
      push @CalculatedValues, '';
    }
    if ($FormulaCol > @LineWords) {
      $ErrorMsg = "Ignoring line $LineCount: Formula column $ColLabels[$FormulaCol] not found";
      PrintErrorMsg($Line, $ErrorMsg);
      ComposeAndWriteNewLine(\*NEWTEXTFILE, \@LineWords, \@ColNumsBeforeNew, \@ColNumsAfterNew, \@CalculatedValues);
      next TEXTLINE;
    }

    # Make sure it's a valid molecular formula...
    $FormulaColValue = $LineWords[$FormulaCol];
    if ($OptionsInfo{CheckFormula}) {
      ($Status, $ErrorMsg) = MolecularFormula::IsMolecularFormula($FormulaColValue);
      if (!$Status) {
	$ErrorMsg = "Ignoring line $LineCount: Formula column $ColLabels[$FormulaCol] value is not valid: $ErrorMsg";
	PrintErrorMsg($Line, $ErrorMsg);
	ComposeAndWriteNewLine(\*NEWTEXTFILE, \@LineWords, \@ColNumsBeforeNew, \@ColNumsAfterNew, \@CalculatedValues);
	next TEXTLINE;
      }
    }

    # Calculate appropriate values and write 'em out...
    @CalculatedValues = ();
    for $CalculationType (@{$OptionsInfo{SpecifiedCalculations}}) {
      if ($CalculationType =~ /^ElementalAnalysis$/i) {
	($ElementsRef, $ElementCompositionRef) = MolecularFormula::CalculateElementalComposition($FormulaColValue);
	$CalculatedValue = (defined($ElementsRef) && defined($ElementCompositionRef)) ? MolecularFormula::FormatCompositionInfomation($ElementsRef, $ElementCompositionRef, $OptionsInfo{Precision}) : '';
      }
      elsif ($CalculationType =~ /^MolecularWeight$/i) {
	$CalculatedValue = MolecularFormula::CalculateMolecularWeight($FormulaColValue);
	$CalculatedValue = (defined($CalculatedValue) && length($CalculatedValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $CalculatedValue)) : "";
      }
      elsif ($CalculationType =~ /^ExactMass$/i) {
	$CalculatedValue = MolecularFormula::CalculateExactMass($FormulaColValue);
	$CalculatedValue = (defined($CalculatedValue) && length($CalculatedValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $CalculatedValue)) : "";
      }
      else {
	$CalculatedValue = '';
      }
      push @CalculatedValues, $CalculatedValue;
    }

    ComposeAndWriteNewLine(\*NEWTEXTFILE, \@LineWords, \@ColNumsBeforeNew, \@ColNumsAfterNew, \@CalculatedValues);
  }
  close NEWTEXTFILE;
  close TEXTFILE;

}

# Write out new line using old and new calculated data...
sub ComposeAndWriteNewLine {
  my($NewTextFileRef, $OldLineWordsRef, $ColNumsBeforeNewRef, $ColNumsAfterNewRef, $CalculatedValuesRef) = @_;
  my($NewLine, $ColNum, $Value, @NewLineWords);

  @NewLineWords = ();
  for $ColNum (@{$ColNumsBeforeNewRef}) {
    push @NewLineWords, $OldLineWordsRef->[$ColNum];
  }
  for $Value (@{$CalculatedValuesRef}) {
    push @NewLineWords, $Value;
  }
  for $ColNum (@{$ColNumsAfterNewRef}) {
    push @NewLineWords, $OldLineWordsRef->[$ColNum];
  }
  $NewLine = JoinWords(\@NewLineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print $NewTextFileRef "$NewLine\n";
}

# Print out error message...
sub PrintErrorMsg {
  my($Line, $ErrorMsg) = @_;

  if ($OptionsInfo{DetailLevel} >= 2 ) {
    print "$ErrorMsg: $Line\n";
  }
  elsif ($OptionsInfo{DetailLevel} >= 1) {
    print "$ErrorMsg\n";
  }
}

# Process formula columns and other information...
sub RetrieveColumnsAndLabelsInfo {
  RetrieveFormulaColumnsInfo();
  RetrieveStartColumnsAndValueLabelsInfo();
}

# Make sure specified formula column are okay...
sub RetrieveFormulaColumnsInfo {
  my($Index, $TextFile);

  @{$TextFilesInfo{FormulaColNum}} = ();

 FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];

    $TextFilesInfo{FormulaColNum}[$Index] = 0;

    if ($TextFilesInfo{FileOkay}[$Index]) {
      my($FormulaColNum, $FormulaColValid);

      $FormulaColNum = 0;
      $FormulaColValid = 0;
      if ($OptionsInfo{SpecifiedFormulaCol}) {
	if ($OptionsInfo{ColMode} =~ /^colnum$/i) {
	  if ($OptionsInfo{SpecifiedFormulaCol} <= $TextFilesInfo{ColCount}[$Index]) {
	    $FormulaColNum = $OptionsInfo{SpecifiedFormulaCol} - 1;
	    $FormulaColValid = 1;
	  }
	}
	else {
	  if (exists($TextFilesInfo{ColLabelToNumMap}[$Index]{$OptionsInfo{SpecifiedFormulaCol}})) {
	    $FormulaColNum =  $TextFilesInfo{ColLabelToNumMap}[$Index]{$OptionsInfo{SpecifiedFormulaCol}};
	    $FormulaColValid = 1;
	  }
	}
      }
      else {
	# Grab the first column with the word Formula in its label...
	my($ColLabel);
	LABEL: for $ColLabel (@{$TextFilesInfo{ColLabels}[$Index]}) {
	  if ($ColLabel =~ /Formula/i) {
	    $FormulaColNum =  $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel};
	    $FormulaColValid = 1;
	    last LABEL;
	  }
	}
      }
      if ($FormulaColValid) {
	$TextFilesInfo{FormulaColNum}[$Index] = $FormulaColNum;
      }
      else {
	if ($OptionsInfo{SpecifiedFormulaCol}) {
	  warn "Warning: Ignoring file $TextFile: Formula column specified, $OptionsInfo{SpecifiedFormulaCol}, using \"f --formulacol\" option doesn't exist\n";
	}
	else {
	  warn "Warning: Ignoring file $TextFile: Column label containing the word Formula doesn't exist\n";
	}
	$TextFilesInfo{FileOkay}[$Index] = 0;
      }
    }
  }
}

# Setup starting column number for adding calculated values and
# column lables to use for these values...
sub RetrieveStartColumnsAndValueLabelsInfo {
  my($Index, $TextFile, $SpecifiedStartColNum, $StartColNum, $Label, $Value, $NewLabel, $Count, $BeforeStartColNum, $AfterStartColNum, $FirstColNum, $LastColNum, $ColNum, $Part1StartColNum, $Part1EndColNum, $Part2StartColNum, $Part2EndColNum, @Part1ColNums, @Part2ColNums);

  # Start column number for inserting new values...
  $SpecifiedStartColNum = "last";
  if (defined($OptionsInfo{StartCol})) {
    if (length($OptionsInfo{StartCol})) {
      $SpecifiedStartColNum = $OptionsInfo{StartCol}
    }
  }

  # Column labels for for new calculated values...
  my(%NewValueLabels) = (ElementalAnalysis => 'ElementalAnalysis', MolecularWeight => 'MolecularWeight', ExactMass => 'ExactMass');
  if (@{$OptionsInfo{SpecifiedValueLabels}}) {
    for ($Index = 0; $Index < @{$OptionsInfo{SpecifiedValueLabels}}; $Index +=2) {
      $Value = $OptionsInfo{SpecifiedValueLabels}[$Index];
      $Label = $OptionsInfo{SpecifiedValueLabels}[$Index + 1];
      if (exists $NewValueLabels{$Value}) {
	$NewValueLabels{$Value} = $Label;
      }
    }
  }

  @{$TextFilesInfo{ColNumsBeforeNew}} = ();
  @{$TextFilesInfo{ColNumsAfterNew}} = ();
  @{$TextFilesInfo{ValueLabelsMap}} = ();

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];

    @{$TextFilesInfo{ColNumsBeforeNew}[$Index]} = ();
    @{$TextFilesInfo{ColNumsAfterNew}[$Index]} = ();
    %{$TextFilesInfo{ValueLabelsMap}[$Index]} = ();

    if (!$TextFilesInfo{FileOkay}[$Index]) {
      next FILELIST;
    }

    if ($SpecifiedStartColNum !~ /^last$/i) {
      if ($OptionsInfo{ColMode} =~ /^collabel$/i) {
	if (exists($TextFilesInfo{ColLabelToNumMap}[$Index]{$SpecifiedStartColNum})) {
	  $StartColNum = $TextFilesInfo{ColLabelToNumMap}[$Index]{$SpecifiedStartColNum};
	}
	else {
	  die "Error: Invalid value $SpecifiedStartColNum specified using \"-s --startcol\" option: column name doesn't exist in  $TextFile  \n";
	}
      }
      else {
	if ($SpecifiedStartColNum > 0 && $SpecifiedStartColNum <= $TextFilesInfo{ColCount}[$Index]) {
	  $StartColNum -= 1;
	}
	else {
	  die "Error: Invalid value $SpecifiedStartColNum specified using \"-s --startcol\" option: column number doesn't exist in  $TextFile  \n";
	}
      }
    }
    else {
      $StartColNum = $TextFilesInfo{ColCount}[$Index] - 1;
    }
    # Set up columns lists for before and after the addition of calculated column values
    # for each text file...
    my($BeforeStartColNum, $AfterStartColNum, $FirstColNum, $LastColNum, $ColNum, $Part1StartColNum, $Part1EndColNum, $Part2StartColNum, $Part2EndColNum, @Part1ColNums, @Part2ColNums);

    $FirstColNum = 0; $LastColNum = $TextFilesInfo{ColCount}[$Index] - 1;

    $BeforeStartColNum = $StartColNum - 1;
    $AfterStartColNum = $StartColNum + 1;

    if ($OptionsInfo{StartColMode} =~ /^after$/i) {
      $Part1StartColNum = $FirstColNum; $Part1EndColNum = $StartColNum;
      $Part2StartColNum = $AfterStartColNum; $Part2EndColNum = $LastColNum;
    }
    else {
      $Part1StartColNum = $FirstColNum; $Part1EndColNum = $BeforeStartColNum;
      $Part2StartColNum = $StartColNum; $Part2EndColNum = $LastColNum;
    }
    @Part1ColNums = (); @Part2ColNums = ();
    for $ColNum (0 .. $TextFilesInfo{ColCount}[$Index]) {
      if ($ColNum >= $Part1StartColNum && $ColNum <= $Part1EndColNum) {
	push @Part1ColNums, $ColNum;
      }
    }
    for $ColNum (0 .. $TextFilesInfo{ColCount}[$Index]) {
      if ($ColNum >= $Part2StartColNum && $ColNum <= $Part2EndColNum) {
	push @Part2ColNums, $ColNum;
      }
    }
    push @{$TextFilesInfo{ColNumsBeforeNew}[$Index]}, @Part1ColNums;
    push @{$TextFilesInfo{ColNumsAfterNew}[$Index]}, @Part2ColNums;

    # Setup column labels for calculated values...
    for $Value (keys %NewValueLabels) {
      $Label = $NewValueLabels{$Value};

      # Make sure it doesn't already exists...
      $Count = 1;
      $NewLabel = $Label;
      while (exists $TextFilesInfo{ColLabelToNumMap}[$Index]{$NewLabel}) {
	$Count++;
	$NewLabel = $Label . $Count;
      }
      $TextFilesInfo{ValueLabelsMap}[$Index]{$Value} = $NewLabel;
    }
  }
}

# Retrieve information about input text files...
sub RetrieveTextFilesInfo {
  my($Index, $TextFile, $FileDir, $FileName, $FileExt, $InDelim, $Line, @ColLabels, $OutFileRoot,  $OutFile, $ColNum, $ColLabel);

  %TextFilesInfo = ();

  @{$TextFilesInfo{FileOkay}} = ();
  @{$TextFilesInfo{ColCount}} = ();
  @{$TextFilesInfo{ColLabels}} = ();
  @{$TextFilesInfo{ColLabelToNumMap}} = ();
  @{$TextFilesInfo{InDelim}} = ();
  @{$TextFilesInfo{OutFile}} = ();

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];

    $TextFilesInfo{FileOkay}[$Index] = 0;
    $TextFilesInfo{ColCount}[$Index] = 0;
    $TextFilesInfo{InDelim}[$Index] = "";
    $TextFilesInfo{OutFile}[$Index] = "";

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
      $OutFileRoot = $FileName . "ElementalAnalysis";
    }

    $OutFile = $OutFileRoot . ".$FileExt";
    if (lc($OutFile) eq lc($TextFile)) {
      warn "Warning: Ignoring file $TextFile:Output file name, $OutFile, is same as input text file name, $TextFile\n";
      next FILELIST;
    }
    if (!$Options{overwrite}) {
      if (-e $OutFile) {
	warn "Warning: Ignoring file $TextFile: The file $OutFile already exists\n";
	next FILELIST;
      }
    }

    $TextFilesInfo{FileOkay}[$Index] = 1;
    $TextFilesInfo{InDelim}[$Index] = $InDelim;
    $TextFilesInfo{OutFile}[$Index] = "$OutFile";

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

  $OptionsInfo{ColMode} = $Options{colmode};
  $OptionsInfo{StartColMode} = $Options{startcolmode};

  $OptionsInfo{Fast} = defined $Options{fast} ? $Options{fast} : undef;

  $OptionsInfo{DetailLevel} = $Options{detail};
  $OptionsInfo{CheckFormula} = $Options{fast} ? 0 : 1;
  $OptionsInfo{Precision} = $Options{precision};

  $OptionsInfo{InDelim} = $Options{indelim};

  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /^tab$/i ) ? "\t" : (($Options{outdelim} =~ /^semicolon$/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^yes$/i) ? 1 : 0;

  $OptionsInfo{Overwrite} = defined $Options{overwrite} ? $Options{overwrite} : undef;
  $OptionsInfo{OutFileRoot} = defined $Options{root} ? $Options{root} : undef;

  $OptionsInfo{StartCol} = defined $Options{startcol} ? $Options{startcol} : undef;

  $OptionsInfo{FormulaCol} = defined $Options{formulacol} ? $Options{formulacol} : undef;
  $OptionsInfo{SpecifiedFormulaCol} = "";

  if (defined $Options{formulacol}) {
    $OptionsInfo{SpecifiedFormulaCol} = $Options{formulacol};
    if ($Options{colmode} =~ /^colnum$/i) {
      if (!IsPositiveInteger($OptionsInfo{SpecifiedFormulaCol})) {
	die "Error: Invalid value $Options{formulacol} specified using \"-f -formulacol\" option: Allowed values: > 0\n";
      }
    }
  }

  # Setup what to calculate...
  @{$OptionsInfo{SpecifiedCalculations}} = ();
  if ($Options{mode} =~ /^All$/i) {
    @{$OptionsInfo{SpecifiedCalculations}} = qw(ElementalAnalysis MolecularWeight ExactMass);
  }
  else {
    my($Mode, $ModeValue, @SpecifiedModeValues);
    $Mode = $Options{mode};
    $Mode =~ s/ //g;
    @SpecifiedModeValues = split /\,/, $Mode;
    for $ModeValue (@SpecifiedModeValues) {
      if ($ModeValue !~ /^(ElementalAnalysis|MolecularWeight|ExactMass)$/i) {
	if ($ModeValue =~ /^All$/i) {
	  die "Error: All value for option \"-m --mode\" is not allowed with other valid values.\n";
	}
	else {
	  die "Error: The value specified, $ModeValue, for option \"-m --mode\" is not valid. Allowed values: ElementalAnalysis, MolecularWeight, or ExactMass\n";
	}
      }
      push @{$OptionsInfo{SpecifiedCalculations}}, $ModeValue;
    }
  }

  $OptionsInfo{ValueColLabels} = defined $Options{valuecollabels} ? $Options{valuecollabels} : undef;
  @{$OptionsInfo{SpecifiedValueLabels}} = ();

  if ($Options{valuecollabels}) {
    my($Value, $Label, @ValueLabels);
    @ValueLabels = split /\,/, $Options{valuecollabels};
    if (@ValueLabels % 2) {
      die "Error: The value specified, $Options{valuecollabels}, for option \"-v --valuecollabels\" is not valid: It must contain even number of comma delimited values\n";
    }
    my($Index);
    for ($Index = 0; $Index < @ValueLabels; $Index +=2) {
      $Value = $ValueLabels[$Index];
      $Value =~ s/ //g;
      $Label = $ValueLabels[$Index + 1];
      if ($Value !~ /^(ElementalAnalysis|MolecularWeight|ExactMass)$/i) {
	die "Error: The value specified, $Value, using option \"-v --valuecollabels\" is not valid. Allowed values: ElementalAnalysis, MolecularWeight, or ExactMass\n";
      }
      push @{$OptionsInfo{SpecifiedValueLabels}}, ($Value, $Label);
    }
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{colmode} = "colnum";
  $Options{detail} = 1;
  $Options{mode} = "All";
  $Options{indelim} = "comma";
  $Options{outdelim} = "comma";
  $Options{precision} = 2;
  $Options{quote} = "yes";
  $Options{startcolmode} = "after";

  if (!GetOptions(\%Options, "colmode|c=s", "detail|d=i", "fast", "formulacol|f=s", "help|h", "indelim=s", "mode|m=s", "outdelim=s", "overwrite|o", "precision|p=i", "quote|q=s", "root|r=s", "startcol|s=s", "startcolmode=s", "valuecollabels|v=s", "workingdir|w=s")) {
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
  if ($Options{indelim} !~ /^(comma|semicolon)$/i) {
    die "Error: The value specified, $Options{indelim}, for option \"--indelim\" is not valid. Allowed values: comma or semicolon\n";
  }
  if ($Options{outdelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{outdelim}, for option \"--outdelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if (!IsPositiveInteger($Options{precision})) {
    die "Error: The value specified, $Options{precision}, for option \"-p --precision\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{quote} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{startcolmode} !~ /^(before|after)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"--startcolmode\" is not valid. Allowed values: before or after\n";
  }
  if (!IsPositiveInteger($Options{detail})) {
    die "Error: The value specified, $Options{detail}, for option \"-d --detail\" is not valid. Allowed values: > 0\n";
  }
}

__END__

=head1 NAME

ElementalAnalysisTextFiles.pl - Perform elemental analysis using formula column in TextFile(s)

=head1 SYNOPSIS

ElementalAnalysisTextFiles.pl TextFile(s)...

ElementalAnalysisTextFiles.pl [B<-c, --colmode> colnum | collabel] [B<-d, --detail> infolevel] [B<-f, --fast>]
[B<-f, --formulacol> colnum | collabel] [B<-h, --help>] [B<--indelim> comma | semicolon]
[B<-m, --mode> All | "ElementalAnysis, [MolecularWeight, ExactMass]"] [B<-o, --overwrite>]
[B<--outdelim> comma | tab | semicolon] [B<-p, --precision> number] [B<-q, --quote> yes | no]
[B<-r, --root> rootname] [B<-s, --startcol> colnum | collabel] [B<--startcolmode> before | after]
B<-v --valuecollabels> [Name, Label, [Name, Label,...]] [B<-w, --workingdir> dirname] TextFile(s)...

=head1 DESCRIPTION

Perform elemental analysis using molecular formula column specified by a column number or label in
I<TextFile(s)>.

In addition to straightforward molecular formulas - H2O, HCl, C3H7O2N -
other supported variations are: Ca3(PO4)2, [PCl4]+, [Fe(CN)6]4-, C37H42N2O6+2, Na2CO3.10H2O,
8H2S.46H2O, and so on. Charges are simply ignored. Isotope symbols in formulas specification, including
D and T, are not supported.

The valid file extensions are I<.csv> and I<.tsv> for comma/semicolon and tab delimited
text files respectively. All other file names are ignored. All the text files in a
current directory can be specified by I<*.csv>, I<*.tsv>, or the current directory
name. The B<--indelim> option determines the format of I<TextFile(s)>. Any file
which doesn't correspond to the format indicated by B<--indelim> option is ignored.

=head1 OPTIONS

=over 4

=item B<-c, --colmode> I<colnum | collabel>

Specify how columns are identified in I<TextFile(s)>: using column number or column
label. Possible values: I<colnum or collabel>. Default value: I<colnum>.

=item B<-d, --detail> I<infolevel>

Level of information to print about lines being ignored. Default: I<1>. Possible values:
I<1, 2 or 3>.

=item B<-h, --help>

Print this help message.

=item B<--fast>

In this mode, the formula column specified using B<-f, --formulacol> option is assumed
to contain valid molecular formula data and initial formula validation check is skipped.

=item B<-f, --formulacol> I<col number | col name>

This value is mode specific. It specifies molecular formula column to use for performing
elemental analysis on I<TextFile(s)>. Possible values: I<col number or col label>.
Default value: I<first column containing the word formula in its column label>.

=item B<-m, --mode> I<All | "ElementalAnalysis,[MolecularWeight,ExactMass]">

Specify what values to calculate using molecular formula in I<TextFile(s)>: calculate all supported
values or specify a comma delimited list of values. Possible values: I<All | "ElementalAnalysis,
[MolecularWeight, ExactMass]">. Default: I<All>

=item B<--indelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile(s)>. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a
delimiter.

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
name: <InitialTextFileName>ElementalAnalysis.<Ext>. The csv, and tsv
<Ext> values are used for comma/semicolon, and tab delimited text files
respectively. This option is ignored for multiple input files.

=item B<-s, --startcol> I<colnum | collabel>

This value is mode specific. It specifies the column in text files which is
used for start adding calculated column values. For I<colnum> mode, specify
column number and for I<collabel> mode, specify column label.

Default value: I<last>. Start merge after the last column.

=item B<--startcolmode> I<before | after>

Start adding calculated column values after the B<-s, --startcol> value. Possible values: I<before or after>.
Default value: I<after>.

=item B<-v --valuecollabels> I<Name,Label,[Name,Label,...]>

Specify column labels to use for calculated values. In general, it's a comma delimited
list of value name and column label pairs. Supported value names: I<ElementalAnalysis,
MolecularWeight,  and ExactMass>. Default labels: I<ElementalAnalysis, MolecularWeight,
and ExactMass>.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To perform elemental analysis, calculate molecular weight and exact mass using formulas
in a column with the word Formula in its column label and generate a new CSV text
file NewSample1.csv, type:

    % ElementalAnalysisTextFiles.pl -o -r NewSample1 Sample1.csv

To perform elemental analysis using formulas in column number two, use column label
Analysis for calculated data, and generate a new CSV text file NewSample1.csv, type:

    % ElementalAnalysisTextFiles.pl --m ElementalAnalysis --formulacol 2
      --valuecollabels "ElementalAnalysis,Analysis" -o -r NewSample1
      Sample1.csv

To calculate molecular weight using formula in column label Formula with four decimal
precision and generate a new CSV text file NewSample1.csv, type

    % ElementalAnalysisTextFiles.pl --m MolecularWeight --colmode collabel
      --formulacol Formula --precision 4 -o -r NewSample1 Sample1.csv

To calculate exact mass  using formula in column label Formula with four decimal
precision, adding column for exact mass right after Formula column, and generate a
new CSV text file NewSample1.csv, type

    % ElementalAnalysisTextFiles.pl --m ExactMass --colmode collabel
      --formulacol Formula --precision 4 --startcolmode after
      --startcol Formula -o -r NewSample1 Sample1.csv


=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AnalyzeTextFilesData.pl, InfoTextFiles.pl, ExtractFromTextFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
