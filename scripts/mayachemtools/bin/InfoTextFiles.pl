#!/usr/bin/perl -w
#
# File: InfoTextFiles.pl
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
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

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
    ListTextFileInfo($FileIndex);
  }
}
ListTotalSizeOfFiles();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# List appropriate information...
sub ListTextFileInfo {
  my($Index) = @_;
  my($TextFile,  $Line, $InDelim, $LineCount, $EmptyLinesCount, $EmptyColDataLinesCount, $GreaterThanMaxColLinesCount, $Label, $Value, $ColNum, $EmptyColValueFound, $PrintTextLine, $NonNumericalDataFound, @ColLabels, @LineWords, %EmptyColValuesCountMap, %NonEmptyColValuesCountMap, %SpecifiedNonNumericalColValuesCountMap, %NonNumericalColValuesCountMap, %NumericalColValuesCountMap,);

  $TextFile = $TextFilesList[$Index];
  $InDelim = $TextFilesInfo{InDelim}[$Index];
  @ColLabels = @{$TextFilesInfo{ColLabels}[$Index]};

  open TEXTFILE, "$TextFile" or die "Error: Can't open $TextFile: $! \n";

  $LineCount = 0;
  $EmptyLinesCount = 0;
  $EmptyColDataLinesCount = 0;
  $GreaterThanMaxColLinesCount = 0;

  %EmptyColValuesCountMap = ();
  %NonEmptyColValuesCountMap = ();
  %SpecifiedNonNumericalColValuesCountMap = ();
  %NonNumericalColValuesCountMap = ();
  %NumericalColValuesCountMap = ();

  if ($OptionsInfo{ParseLines}) {
    # Skip over column labels from old file...
    if (<TEXTFILE>) {
      $LineCount++;
      LINE: while ($Line = <TEXTFILE>) {
	$LineCount++;
	$PrintTextLine = 0;
	$Line =~ s/(\r\n)|(\r)|\n//g;
	@LineWords = quotewords($InDelim, 0, $Line);
	if ($OptionsInfo{CountEmpty}) {
	  # Count lines with no data...
	  if (!@LineWords) {
	    $EmptyLinesCount++;
	    if ($OptionsInfo{DetailLevel} >= 2) {
	      print "Line number $LineCount is empty...\n";
	    }
	    next LINE;
	  }
	  # Count lines with empty data for some columns...
	  $EmptyColValueFound = 0;
	  VALUE: for $Value (@LineWords) {
	      if (!IsNotEmpty($Value)) {
		$EmptyColValueFound = 1;
		next VALUE;
	      }
	  }
	  if ($EmptyColValueFound) {
	    $EmptyColDataLinesCount++;
	    if ($OptionsInfo{DetailLevel} >= 2) {
	      print "Line number $LineCount contains empty column value(s)...\n";
	    }
	    $PrintTextLine = ($OptionsInfo{DetailLevel} >= 3) ? 1 : 0;
	  }
	  # Count lines with columns greater than the column label line...
	  if (@LineWords > @ColLabels) {
	    $GreaterThanMaxColLinesCount++;
	    if ($OptionsInfo{DetailLevel} >= 2) {
	      print "Line number $LineCount contains more than ", scalar(@ColLabels), " columns...\n";
	    }
	    $PrintTextLine = ($OptionsInfo{DetailLevel} >= 3) ? 1 : 0;
	  }
	  # Count empty values for each coulmn...
	  for $ColNum (0 .. $#LineWords) {
	    if ($ColNum < @ColLabels) {
	      $Label = $ColLabels[$ColNum];
	      if (IsNotEmpty($LineWords[$ColNum])) {
		if (exists($NonEmptyColValuesCountMap{$Label})) {
		  $NonEmptyColValuesCountMap{$Label} += 1;
		}
		else {
		  $NonEmptyColValuesCountMap{$Label} = 1;
		}
	      }
	      else {
		$PrintTextLine = ($OptionsInfo{DetailLevel} >= 3) ? 1 : 0;
		if (exists($EmptyColValuesCountMap{$Label})) {
		  $EmptyColValuesCountMap{$Label} += 1;
		}
		else {
		  $EmptyColValuesCountMap{$Label} = 1;
		}
	      }
	    }
	  }
	}
	if ($OptionsInfo{CheckData}) {
	  for $ColNum (0 .. $#LineWords) {
	    if ($ColNum < @ColLabels) {
	      if (IsNumerical($LineWords[$ColNum])) {
		$Label = $ColLabels[$ColNum];
		if (exists($NumericalColValuesCountMap{$Label})) {
		  $NumericalColValuesCountMap{$Label} += 1;
		}
		else {
		  $NumericalColValuesCountMap{$Label} = 1;
		}
	      }
	      else {
		$Label = $ColLabels[$ColNum];
		if (IsNotEmpty($LineWords[$ColNum])) {
		  if (exists($NonNumericalColValuesCountMap{$Label})) {
		    $NonNumericalColValuesCountMap{$Label} += 1;
		  }
		  else {
		    $NonNumericalColValuesCountMap{$Label} = 1;
		  }
		}
	      }
	    }
	  }
	}
	if ($OptionsInfo{CheckNumericalData}) {
	  $NonNumericalDataFound = 0;
	  for $ColNum (@{$TextFilesInfo{NumericalDataColNums}[$Index]}) {
	    if ($ColNum < @LineWords) {
	      if (!IsNumerical($LineWords[$ColNum])) {
		$NonNumericalDataFound = 1;
		$Label = $ColLabels[$ColNum];
		if (exists($SpecifiedNonNumericalColValuesCountMap{$Label})) {
		  $SpecifiedNonNumericalColValuesCountMap{$Label} += 1;
		}
		else {
		  $SpecifiedNonNumericalColValuesCountMap{$Label} = 1;
		}
	      }
	    }
	  }
	  if ($NonNumericalDataFound) {
	    $PrintTextLine = ($OptionsInfo{DetailLevel} >= 3) ? 1 : 0;
	    if ($OptionsInfo{DetailLevel} >=2 ) {
	      print "Line number $LineCount contains non-numerical data for some specified column(s)...\n";
	    }
	  }
	}
	if ($PrintTextLine) {
	  print "Line $LineCount: $Line\n\n";
	}
      }
    }
  }
  else {
    while (<TEXTFILE>) {
      $LineCount++;
    }
  }
  close TEXTFILE;

  print "\nNumber of lines: $LineCount\n";
  print "Number of columns: $TextFilesInfo{ColCount}[$Index]\n";
  print "Column labels: ", JoinWords(\@ColLabels, ", ", 1), "\n";

  if ($OptionsInfo{CountEmpty}) {
    print "\nNumber of lines with no data: $EmptyLinesCount\n";
    print "Number of lines with some missing column data: $EmptyColDataLinesCount\n";
    print "Number of lines containing greater than ", scalar(@ColLabels), " columns: $GreaterThanMaxColLinesCount\n";
    PrintDataInformation("Number of non-empty values for each column(s)", \@ColLabels, \%NonEmptyColValuesCountMap);
    PrintDataInformation("Number of empty values for each column(s)", \@ColLabels, \%EmptyColValuesCountMap);
  }

  if ($OptionsInfo{CheckData}) {
    print "\n";
    PrintDataInformation("Number of non-numerical data values for each column(s)", \@ColLabels, \%NonNumericalColValuesCountMap);
    PrintDataInformation("Number of numerical data values for each column(s)", \@ColLabels, \%NumericalColValuesCountMap);
    print "\n";
  }

  if ($OptionsInfo{CheckNumericalData} && @{$TextFilesInfo{NumericalDataColLabels}[$Index]}) {
    PrintDataInformation("Number of non-numerical data values for each column(s)", \@{$TextFilesInfo{NumericalDataColLabels}[$Index]}, \%SpecifiedNonNumericalColValuesCountMap);
  }

  # File size and modification information...
  print "\nFile size: ", FormatFileSize($TextFilesInfo{FileSize}[$Index]), " \n";
  print "Last modified: ", $TextFilesInfo{FileLastModified}[$Index], " \n";
}

# Total size of all the fiels...
sub ListTotalSizeOfFiles {
  my($FileOkayCount, $TotalSize, $Index);

  $FileOkayCount = 0;
  $TotalSize = 0;

  for $Index (0 .. $#TextFilesList) {
    if ($TextFilesInfo{FileOkay}[$Index]) {
      $FileOkayCount++;
      $TotalSize += $TextFilesInfo{FileSize}[$Index];
    }
  }
  if ($FileOkayCount > 1) {
    print "\nTotal size of $FileOkayCount files: ", FormatFileSize($TotalSize), "\n";
  }
}

# List data information...
sub PrintDataInformation {
  my($InfoLabel, $DataLabelRef, $DataLabelToValueMapRef) = @_;
  my($Line, $Label);

  $Line = "";
  for $Label (@{$DataLabelRef}) {
    $Line .= " \"$Label\" - " . (exists($DataLabelToValueMapRef->{$Label}) ? $DataLabelToValueMapRef->{$Label} : 0) . ",";
  }
  $Line =~ s/\,$//g;
  print "$InfoLabel: $Line\n";
}

# Retrieve information about input text files...
sub RetrieveTextFilesInfo {
  my($Index, $TextFile, $FileDir, $FileName, $FileExt, $InDelim, $Line, @ColLabels,  $ColNum, $ColLabel, $ModifiedTimeString, $ModifiedDateString);

  %TextFilesInfo = ();
  @{$TextFilesInfo{FileOkay}} = ();
  @{$TextFilesInfo{ColCount}} = ();
  @{$TextFilesInfo{ColLabels}} = ();
  @{$TextFilesInfo{ColLabelToNumMap}} = ();
  @{$TextFilesInfo{InDelim}} = ();
  @{$TextFilesInfo{FileSize}} = ();
  @{$TextFilesInfo{FileLastModified}} = ();

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];

    $TextFilesInfo{FileOkay}[$Index] = 0;
    $TextFilesInfo{ColCount}[$Index] = 0;
    $TextFilesInfo{InDelim}[$Index] = "";
    $TextFilesInfo{FileSize}[$Index] = 0;
    $TextFilesInfo{FileLastModified}[$Index] = '';
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
      if ($OptionsInfo{InDelim} !~ /^(comma|semicolon)$/i) {
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

    $TextFilesInfo{FileOkay}[$Index] = 1;
    $TextFilesInfo{InDelim}[$Index] = $InDelim;

    $TextFilesInfo{ColCount}[$Index] = @ColLabels;
    push @{$TextFilesInfo{ColLabels}[$Index]}, @ColLabels;
    for $ColNum (0 .. $#ColLabels) {
      $ColLabel = $ColLabels[$ColNum];
      $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel} = $ColNum;
    }
    $TextFilesInfo{FileSize}[$Index] = FileSize($TextFile);
    ($ModifiedTimeString, $ModifiedDateString) = FormattedFileModificationTimeAndDate($TextFile);
    $TextFilesInfo{FileLastModified}[$Index] = "$ModifiedTimeString; $ModifiedDateString";
  }

}

# Make sure specified numerical data columns are okay...
sub ProcessColumnsInfo {
  my($Index, $TextFile);

  @{$TextFilesInfo{NumericalDataColNums}} = ();
  @{$TextFilesInfo{NumericalDataColLabels}} = ();

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];
    @{$TextFilesInfo{NumericalDataColNums}[$Index]} = ();
    @{$TextFilesInfo{NumericalDataColLabels}[$Index]} = ();

    if ($TextFilesInfo{FileOkay}[$Index]) {
      my($SpecifiedColNum, $ColNum, $ColLabel, @SpecifiedColNums, @SpecifiedColLabels);
      @SpecifiedColNums = ();
      if ($OptionsInfo{Mode} =~ /^colnum$/i) {
	for $SpecifiedColNum (@{$OptionsInfo{SpecifiedNumericalDataCols}}) {
	  if ($SpecifiedColNum <= $TextFilesInfo{ColCount}[$Index]) {
	    $ColNum = $SpecifiedColNum - 1;
	    push @SpecifiedColNums, $ColNum;
	    push @SpecifiedColLabels, $TextFilesInfo{ColLabels}[$Index][$ColNum];
	  }
	}
      }
      else {
	for $ColLabel (@{$OptionsInfo{SpecifiedNumericalDataCols}}) {
	  if (exists($TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel})) {
	    $ColNum = $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel};
	    push @SpecifiedColNums, $ColNum;
	    push @SpecifiedColLabels, $ColLabel;
	  }
	}
      }
      if (@SpecifiedColNums) {
	push @{$TextFilesInfo{NumericalDataColNums}[$Index]}, @SpecifiedColNums;
	push @{$TextFilesInfo{NumericalDataColLabels}[$Index]}, @SpecifiedColLabels;
      }
    }
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{All} = $Options{all} ? $Options{all} : 0;
  $OptionsInfo{Count} = $Options{count} ? $Options{count} : 0;

  $OptionsInfo{DetailLevel} = $Options{detail} ? $Options{detail} : 1;

  $OptionsInfo{Empty} = $Options{empty} ? $Options{empty} : 0;

  $OptionsInfo{InDelim} = $Options{indelim};
  $OptionsInfo{NumericalDataCols} = $Options{numericaldatacols} ? $Options{numericaldatacols} : 0;

  $OptionsInfo{ParseLines} = ($Options{all} || $Options{empty} || $Options{numericaldatacols}) ? 1 : 0;
  $OptionsInfo{CountEmpty} = ($Options{all} || $Options{empty}) ? 1 : 0;
  $OptionsInfo{CheckData} = ($Options{all} || $Options{datacheck}) ? 1 : 0;
  $OptionsInfo{CheckNumericalData} = ($Options{all} || $Options{numericaldatacols}) ? 1 : 0;

  @{$OptionsInfo{SpecifiedNumericalDataCols}} = ();
  if ($Options{numericaldatacols}) {
    @{$OptionsInfo{SpecifiedNumericalDataCols}} = split ",", $Options{numericaldatacols};
    if ($Options{mode} =~ /^colnum$/i) {
      my($ColNum);
      for $ColNum (@{$OptionsInfo{SpecifiedNumericalDataCols}}) {
	if (!IsPositiveInteger($ColNum)) {
	  die "Error: Invalid value $ColNum specified using \"--numericaldatacols\" option: Allowed values: > 0\n";
	}
      }
    }
  }

}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{detail} = 1;
  $Options{mode} = "colnum";
  $Options{indelim} = "comma";
  if (!GetOptions(\%Options, "all|a", "count|c", "datacheck", "detail|d=i", "empty|e", "help|h", "indelim=s", "mode|m=s", "numericaldatacols|n=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{mode} !~ /^(colnum|collabel)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: colnum or collabel\n";
  }
  if ($Options{indelim} !~ /^(comma|semicolon)$/i) {
    die "Error: The value specified, $Options{indelim}, for option \"--indelim\" is not valid. Allowed values: comma or semicolon\n";
  }
  if (!IsPositiveInteger($Options{detail})) {
    die "Error: The value specified, $Options{detail}, for option \"-d --detail\" is not valid. Allowed values: > 0\n";
  }
}

__END__

=head1 NAME

InfoTextFiles.pl - List information about TextFile(s)

=head1 SYNOPSIS

InfoTextFiles.pl TextFile(s)...

InfoTextFiles.pl [B<-a, --all>] [B<-c, --count>] [B<--datacheck>] [B<-d, --detail> infolevel] [B<-e, --empty>]
[B<-h, --help>] [B<--indelim> comma | semicolon] [B<-m, --mode> colnum | collabel]
[B<-n, --numericaldatacols> colnum,[colnum,...] | collabel,[collabel,...]]
[B<-w, --workingdir> dirname] TextFile(s)...

=head1 DESCRIPTION

List information about I<TextFile(s)> contents: number of lines and columns, empty
column values, and so on. The file names are separated by spaces.
The valid file extensions are I<.csv> and I<.tsv> for comma/semicolon and tab delimited
text files respectively. All other file names are ignored. All the text files in a
current directory can be specified by I<*.csv>, I<*.tsv>, or the current directory
name. The B<--indelim> option determines the format of I<TextFile(s)>. Any file
which doesn't correspond to the format indicated by B<--indelim> option is ignored.

=head1 OPTIONS

=over 4

=item B<-a, --all>

List all the available information.

=item B<-c, --count>

List number of rows and columns. This is B<default behavior>.

=item B<--datacheck>

List number of numerical and non-numerical values for each column.

=item B<-d, --detail> I<infolevel>

Level of information to print about lines being ignored. Default: I<1>. Possible values:
I<1, 2 or 3>.

=item B<-e, --empty>

List number of empty row and column values.

=item B<-h, --help>

Print this help message.

=item B<--indelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile(s)>. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a
delimiter.

=item B<-m, --mode> I<colnum | collabel>

Specify how to identify numerical data columns: using column number or column label.
Possible values: I<colnum or collabel>. Default value: I<colnum>.

=item B<-n, --numericaldatacols> I<colnum,[colnum,...] | collabel,[collabel,...]>

This value is mode specific. It is a list of column number or labels to check for
presence of numerical data only; otherwise, the value is flagged. Default value: I<all;all;...>.

For I<colnum> mode, input value format is: I<colnum,...;colnum,...;...>. Example:

    1,3,5
    "2,4,6"

For I<collabel> mode, input value format is: I<collabel,...;collabel,...;...>. Example:

    "MW,SumNO,SumNHOH"


=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To count number of lines and columns in Text file(s), type:

    % InfoTextFiles.pl Sample1.csv
    % InfoTextFiles.pl Sample1.csv Sample1.tsv
    % InfoTextFiles.pl *.csv *.tsv

To count number of lines, columns and empty values in Sample1.csv file and print
detailed information, type:

    % InfoTextFiles.pl -d 3 -e Sample1.csv

To track all available information and non-numerical values for Mol_ID and MolWeight
columns in Sample1.csv file and print detailed information, type:

    % InfoTextFiles.pl -d 3 -a -m collabel -n Mol_ID,MolWeight Sample1.csv

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
