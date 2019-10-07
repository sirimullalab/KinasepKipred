#!/usr/bin/perl -w
#
# File: MergeTextFilesWithSD.pl
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
use SDFileUtil;
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

my($SDFile, @TextFilesList);
@TextFilesList = ExpandFileNames(\@ARGV, "csv tsv");

if (@TextFilesList < 2) {
  die "Error: Specify one or more text files.\n";
}
$SDFile = shift @TextFilesList;

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Setup information about input files...
print "Checking input SD and text files...\n";
my(%TextFilesInfo);
ProcessSDFileInfo();
RetrieveTextFilesInfo();
RetrieveColumnsAndKeysInfo();

# Merge files...
print "\nGenerating new SD file $OptionsInfo{NewSDFile}...\n";
MergeTextFilesWithSD();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Merge all valid Text files with SD file...
sub MergeTextFilesWithSD {
  my($Index, $Line);

  open NEWSDFILE, ">$OptionsInfo{NewSDFile}" or die "Error: Couldn't open $OptionsInfo{NewSDFile}: $! \n";

  open SDFILE, "$SDFile" or die "Error: Couldn't open $SDFile: $! \n";

  @{$TextFilesInfo{FileHandle}} = ();
  for $Index (0 .. $#TextFilesList) {
    $TextFilesInfo{FileHandle}[$Index] = new FileHandle;

    open $TextFilesInfo{FileHandle}[$Index], "$TextFilesList[$Index]" or die "Error: Couldn't open $TextFilesList[$Index]: $! \n";
    GetTextLine($TextFilesInfo{FileHandle}[$Index]);
  }

  if ($OptionsInfo{Keys}) {
    MergeTextColumnValuesUsingKeys(\*NEWSDFILE, \*SDFILE);
  }
  else {
    MergeTextColumnValues(\*NEWSDFILE, \*SDFILE);
  }

  # Close all opened files...
  close NEWSDFILE;
  close SDFILE;
  for $Index (0 .. $#TextFilesList) {
    close $TextFilesInfo{FileHandle}[$Index];
  }
}

# Merge the specified text columns into SD file...
sub MergeTextColumnValues {
  my($NewSDFileRef, $SDFileRef) = @_;
  my($Index, $Value, $CmpdString, $Line, $InDelim, $ColNum, $ColIndex, @ColLabels, @ColValues, @LineWords);

  while ($CmpdString = ReadCmpdString($SDFileRef)) {
    $CmpdString =~ s/\$\$\$\$$//g;
    print $NewSDFileRef "$CmpdString";

    # Merge coulmn values from other text files...
    @ColLabels = (); @ColValues = ();
    for $Index (0 .. $#TextFilesList) {
      push @ColLabels, @{$TextFilesInfo{ColToMergeLabels}[$Index]};
      $InDelim = $TextFilesInfo{InDelim}[$Index];

      if ($Line = GetTextLine($TextFilesInfo{FileHandle}[$Index])) {
	@LineWords = quotewords($InDelim, 0, $Line);

	for $ColNum (@{$TextFilesInfo{ColToMerge}[$Index]}) {
	  $Value = ($ColNum < @LineWords) ? $LineWords[$ColNum] : "";
	  push @ColValues, $Value;
	}
      }
    }

    for $ColIndex (0 .. $#ColLabels) {
      print $NewSDFileRef "> <$ColLabels[$ColIndex]>\n$ColValues[$ColIndex]\n\n";
    }
    print $NewSDFileRef "\$\$\$\$\n";
  }
}

# Merge the specified text columns into SD file using keys...
sub MergeTextColumnValuesUsingKeys {
  my($NewSDFileRef, $SDFileRef) = @_;
  my($Index, $CmpdString, $Value, $InDelim, $KeyColNum, $KeyColValue, $Line, $ColIndex, $ColNum, @ColLabels, @ColValues, @LineWords, @CmpdLines, @TextFilesKeysToLinesMap, %DataFieldValues);

  # Retrieve text lines from all the text files...
  @TextFilesKeysToLinesMap = ();

  for $Index (0 .. $#TextFilesList) {
    $InDelim = $TextFilesInfo{InDelim}[$Index];
    %{$TextFilesKeysToLinesMap[$Index]} = ();
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

  while ($CmpdString = ReadCmpdString($SDFileRef)) {
    @CmpdLines = split "\n", $CmpdString;
    %DataFieldValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

    if (exists($DataFieldValues{$OptionsInfo{SDKey}})) {
      @ColLabels = (); @ColValues = ();
      $CmpdString =~ s/\$\$\$\$$//g;
      print $NewSDFileRef "$CmpdString";

      $KeyColValue = $DataFieldValues{$OptionsInfo{SDKey}};

      # Merge coulmn values from other text files...
      for $Index (0 .. $#TextFilesList) {
	push @ColLabels, @{$TextFilesInfo{ColToMergeLabels}[$Index]};
	@LineWords = ();

	if (exists($TextFilesKeysToLinesMap[$Index]{$KeyColValue})) {
	  push @LineWords, @{$TextFilesKeysToLinesMap[$Index]{$KeyColValue}};
	}

	for $ColNum (@{$TextFilesInfo{ColToMerge}[$Index]}) {
	  $Value = ($ColNum < @LineWords) ? $LineWords[$ColNum] : "";
	  push @ColValues, $Value;
	}
      }

      for $ColIndex (0 .. $#ColLabels) {
	$Value = (($ColIndex < @ColValues) && IsNotEmpty($ColValues[$ColIndex]) ) ? $ColValues[$ColIndex] : "";
	print $NewSDFileRef "> <$ColLabels[$ColIndex]>\n$Value\n\n";
      }
      print $NewSDFileRef "\$\$\$\$\n";
    }
  }
}

# Retrieve text file columns and keys information for specified options...
sub RetrieveColumnsAndKeysInfo {
  ProcessColumnsInfo();

  if ($OptionsInfo{Keys}) {
    ProcessKeysInfo();
  }
}

# Process specified columns...
sub ProcessColumnsInfo {
  my($Index, $Values, $ColIndex, $ColNum, $ColLabel, @Words);

  @{$TextFilesInfo{ColSpecified}} = ();
  @{$TextFilesInfo{ColToMerge}} = ();
  @{$TextFilesInfo{ColToMergeLabels}} = ();
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
      } else {
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
	$ColNum = $TextFilesInfo{ColSpecified}[$Index][$ColIndex];

	# Make sure it's a numeric value...
	if (!IsPositiveInteger($ColNum)) {
	  warn "Warning: Ignoring value, $ColNum, specified using \"-c --column\" option:  Allowed integer values: > 0\n";
	}
	else {
	  if ($ColNum > 0 && $ColNum <= $TextFilesInfo{ColCount}[$Index]) {
	    $ColNum -= 1;
	    push @{$TextFilesInfo{ColToMerge}[$Index]}, $ColNum;
	    $TextFilesInfo{ColToMergeNumToLabelMap}[$Index]{$ColNum} = $TextFilesInfo{ColLabels}[$Index][$ColNum];
	  }
	  else {
	    warn "Warning: Ignoring value, $ColNum, specified using \"-c --column\" option: column number doesn't exist in  $TextFilesList[$Index]  \n";
	  }
	}
      }
    }

    my (@TextFilesColToMergeSorted) = sort { $a <=> $b } @{$TextFilesInfo{ColToMerge}[$Index]};

    @{$TextFilesInfo{ColToMerge}[$Index]} = ();
    push @{$TextFilesInfo{ColToMerge}[$Index]}, @TextFilesColToMergeSorted;

    # Set up the labels...
    @{$TextFilesInfo{ColToMergeLabels}[$Index]} = ();
    for $ColNum (@TextFilesColToMergeSorted) {
      push @{$TextFilesInfo{ColToMergeLabels}[$Index]}, $TextFilesInfo{ColToMergeNumToLabelMap}[$Index]{$ColNum};
    }
  }
}

# Process specified keys....
sub ProcessKeysInfo {
  my($Index, $ColNum, $ColLabel, $Key);

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
      if (!IsPositiveInteger($ColNum)) {
	warn "Warning: Ignoring value, $ColNum, specified using \"-k --keys\" option: Allowed integer values: > 0  \n";
      }
      else {
	if ($ColNum > 0 && $ColNum <= $TextFilesInfo{ColCount}[$Index]) {
	  $TextFilesInfo{KeysToUse}[$Index] = $ColNum - 1;
	}
	else {
	  warn "Warning: Ignoring value, $ColNum, specified using \"-k --keys\" option: column number doesn't exist in  $TextFilesList[$Index]  \n";
	}
      }
    }
  }

  # Modify columns to merge list to make sure the columns identified by key are taken off the list
  my(@TextFilesColToMergeFiltered, @TextFilesColToMergeLabelsFiltered);

  for $Index (0 .. $#TextFilesList) {
    @TextFilesColToMergeFiltered = ();
    @TextFilesColToMergeLabelsFiltered = ();

    for $ColNum (@{$TextFilesInfo{ColToMerge}[$Index]}) {
      if ($TextFilesInfo{KeysToUse}[$Index] != $ColNum) {
	push @TextFilesColToMergeFiltered, $ColNum;
	push @TextFilesColToMergeLabelsFiltered, $TextFilesInfo{ColToMergeNumToLabelMap}[$Index]{$ColNum};
      }
    }

    @{$TextFilesInfo{ColToMerge}[$Index]} = ();
    push @{$TextFilesInfo{ColToMerge}[$Index]}, @TextFilesColToMergeFiltered;

    @{$TextFilesInfo{ColToMergeLabels}[$Index]} = ();
    push @{$TextFilesInfo{ColToMergeLabels}[$Index]}, @TextFilesColToMergeLabelsFiltered;
  }
}

# Check SD file...
sub ProcessSDFileInfo {
  if (!CheckFileType($SDFile, "sd sdf")) {
    die "Error: Invalid first file $SDFile: It's not a SD file\n";
  }
  if (!(-e $SDFile)) {
    die "Error: SDFile $SDFile doesn't exist\n";
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
  my($Index, $FileDir, $FileName, $FileExt, $NewSDFile, @ColValues, @KeyValues);

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

  $OptionsInfo{OutFileRoot} = defined $Options{root} ? $Options{root} : undef;
  $OptionsInfo{Overwrite} = defined $Options{overwrite} ? $Options{overwrite} : undef;

  $OptionsInfo{SDKey} = defined $Options{sdkey} ? $Options{sdkey} : undef;

  # Setup new SD file...
  if ($Options{root}) {
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($Options{root});
    if ($FileName && $FileExt) {
      $NewSDFile = $FileName;
    }
    else {
      $NewSDFile =  $Options{root};
    }
  }
  else {
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($SDFile);

    $NewSDFile = $FileName;
    ($FileDir, $FileName, $FileExt) = ParseFileName($TextFilesList[0]);

    $NewSDFile = $NewSDFile . "MergedWith" . $FileName . "1To" . @TextFilesList;
  }

  $NewSDFile = $NewSDFile . ".sdf";
  if (!$Options{overwrite}) {
    if (-e $NewSDFile) {
      die "Error: The file $NewSDFile already exists.\n";
    }
  }
  if ($Options{root}) {
    if (lc($NewSDFile) eq lc($SDFile)) {
      die "Error: Output filename, $NewSDFile, is similar to a input file name.\nSpecify a different name using \"-r --root\" option or use default name.\n";
    }
  }
  $OptionsInfo{NewSDFile} = $NewSDFile;
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{mode} = "colnum";
  $Options{indelim} = "comma";

  if (!GetOptions(\%Options, "help|h", "indelim=s", "columns|c=s", "keys|k=s", "mode|m=s", "overwrite|o", "root|r=s", "sdkey|s=s", "workingdir|w=s")) {
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
  if ($Options{sdkey} && !$Options{keys}) {
    die "Error: The option \"-s --sdkey\" can't be specified without the \"-k --keys\" option.\n";
  }
  elsif (!$Options{sdkey} && $Options{keys}) {
    die "Error: The option \"-k --keys\" can't be specified without the \"-s --sdkey\" option.\n";
  }
}

__END__

=head1 NAME

MergeTextFilesWithSD.pl - Merge CSV or TSV TextFile(s) into SDFile

=head1 SYNOPSIS

MergeTextFilesWithSD.pl  SDFile TextFile(s)...

MergeTextFilesWithSD.pl  [B<-h, --help>] [B<--indelim> comma | semicolon]
[B<-c, --columns> colnum,...;... | collabel,...;...] [B<-k, --keys> colkeynum;... | colkeylabel;...]
[B<-m, --mode> colnum | collabel] [B<-o, --overwrite>] [B<-r, --root> rootname]
[B<-s, --sdkey> sdfieldname] [B<-w, --workingdir> dirname] SDFile TextFile(s)...

=head1 DESCRIPTION

Merge multiple CSV or TSV I<TextFile(s)> into I<SDFile>. Unless B<-k --keys>
option is used, data rows from all I<TextFile(s)> are added to I<SDFile> in a
sequential order, and the number of compounds in I<SDFile> is used to determine
how many rows of data are added from I<TextFile(s)>.

Multiple I<TextFile(s)> names are separated by spaces. The valid file extensions are I<.csv> and
I<.tsv> for comma/semicolon and tab delimited text files respectively. All other file names
are ignored. All the text files in a current directory can be specified by I<*.csv>,
I<*.tsv>, or the current directory name. The B<--indelim> option determines the
format of I<TextFile(s)>. Any file which doesn't correspond to the format indicated
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

This value is mode specific. It is a list of columns to merge into I<SDFile>
specified by column numbers or labels for each text file delimited by ";".
All I<TextFile(s)> are merged into I<SDFile>.

Default value: I<all;all;...>. By default, all columns from TextFile(s) are
merged into I<SDFile>.

For I<colnum> mode, input value format is: I<colnum,...;colnum,...;...>. Example:

    "1,2;1,3,4;7,8,9"

For I<collabel> mode, input value format is: I<collabel,...;collabel,...;...>. Example:

    "MW,SumNO;SumNHOH,ClogP,PSA;MolName,Mol_Id,Extreg"

=item B<-k, --keys> I<colkeynum;... | colkeylabel;...>

This value is mode specific. It specifies column keys to use for merging
I<TextFile(s)> into I<SDFile>. The column keys, delimited by ";", are specified by column
numbers or labels for I<TextFile(s)>.

By default, data rows from I<TextFile(s)> are merged into I<SDFile> in the order they appear.

For I<colnum> mode, input value format is:I<colkeynum, colkeynum;...>. Example:

    "1;3;7"

For I<collabel> mode, input value format is:I<colkeylabel, colkeylabel;...>. Example:

    "Mol_Id;Mol_Id;Cmpd_Id"

=item B<-m, --mode> I<colnum | collabel>

Specify how to merge I<TextFile(s)> into I<SDFile>: using column numbers or column labels.
Possible values: I<colnum or collabel>. Default value: I<colnum>.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

New SD file name is generated using the root: <Root>.sdf. Default file name:
<InitialSDFileName>MergedWith<FirstTextFileName>1To<Count>.sdf.

=item B<-s, --sdkey> I<sdfieldname>

I<SDFile> data field name used as a key to merge data from TextFile(s). By default,
data rows from I<TextFile(s)> are merged into I<SDFile> in the order they appear.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To merge Sample1.csv and Sample2.csv into Sample.sdf and generate
NewSample.sdf, type:

    % MergeTextFileswithSD.pl -r NewSample -o Sample.sdf
      Sample1.csv Sample2.csv

To merge all Sample*.tsv into Sample.sdf and generate NewSample.sdf file, type:

    % MergeTextFilesWithSD.pl -r NewSample -o Sample.sdf
      Sample*.tsv

To merge column numbers "1,2" and "3,4,5" from Sample2.csv and Sample3.csv
into Sample.sdf and to generate NewSample.sdf, type:

    % MergeTextFilesWithSD.pl -r NewSample -m colnum -c "1,2;3,4,5"
      -o Sample.sdf Sample1.csv Sample2.csv

To merge column "Mol_ID,Formula,MolWeight" and "Mol_ID,ChemBankID,NAME"
from Sample1.csv and Sample2.csv into Sample.sdf using "Mol_ID" as SD and column keys
to generate NewSample.sdf, type:

    % MergeTextFilesWithSD.pl -r NewSample -s Mol_ID -k "Mol_ID;Mol_ID"
      -m collabel -c "Mol_ID,Formula,MolWeight;Mol_ID,ChemBankID,NAME"
      -o Sample1.sdf Sample1.csv Sample2.csv

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

ExtractFromSDFiles.pl, FilterSDFiles.pl, InfoSDFiles.pl, JoinSDFiles.pl, JoinTextFiles.pl,
MergeTextFiles.pl, ModifyTextFilesFormat.pl, SplitSDFiles.pl, SplitTextFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
