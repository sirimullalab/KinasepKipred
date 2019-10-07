#!/usr/bin/perl -w
#
# File: SortTextFiles.pl
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
    SortTextFile($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Sort it out...
sub SortTextFile {
  my($Index) = @_;
  my($TextFile, $NewTextFile, $KeyCol, $Line, $KeyColValue, $InDelim, @ColLabels, @LineWords);

  $TextFile = $TextFilesList[$Index];
  $InDelim = $TextFilesInfo{InDelim}[$Index];
  $NewTextFile = $TextFilesInfo{OutFile}[$Index];
  $KeyCol = $TextFilesInfo{KeyColNum}[$Index];
  @ColLabels = @{$TextFilesInfo{ColLabels}[$Index]};

  print "Generating new Text file $NewTextFile...\n";
  open NEWTEXTFILE, ">$NewTextFile" or die "Error: Couldn't open $NewTextFile: $! \n";
  open TEXTFILE, "$TextFile" or die "Error: Can't open $TextFile: $! \n";

  # Skip over column labels from old file...
  $Line = GetTextLine(\*TEXTFILE);

  # Add column lablels in new file...
  $Line = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print NEWTEXTFILE "$Line\n";

  # Go over all rows and store the lines using key value as hash...
  my(%KeyToLinesMap, @InvalidDataLines, $LineCount);

  %KeyToLinesMap = ();
  @InvalidDataLines = ();
  $LineCount = 1;
  TEXTLINE: while ($Line = GetTextLine(\*TEXTFILE)) {
    @LineWords = quotewords($InDelim, 0, $Line);
    $LineCount++;
    if ($KeyCol < @LineWords) {
      $KeyColValue = $LineWords[$KeyCol];
      if (!IsNotEmpty($KeyColValue)) {
	$Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
	push @InvalidDataLines, $Line;
	if ($OptionsInfo{DetailLevel} >= 3 ) {
	  print "Ignoring line $LineCount: Contains empty value for key column $ColLabels[$KeyCol]: $Line\n";
	}
	elsif ($OptionsInfo{DetailLevel} >= 2) {
	  print "Ignoring line $LineCount: Contains empty value for key column $ColLabels[$KeyCol]...\n";
	}
	next TEXTLINE;
      }
      if ($OptionsInfo{KeyData} =~ /^numeric$/i) {
	if (!IsFloat($KeyColValue)) {
	  $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
	  push @InvalidDataLines, $Line;
	  if ($OptionsInfo{DetailLevel} >= 3 ) {
	    print "Line number $LineCount: Contains non-numerical value for key column $ColLabels[$KeyCol]: $Line\n";
	  }
	  elsif ($OptionsInfo{DetailLevel} >= 2) {
	    print "Line number $LineCount: Contains non-numerical value for key column $ColLabels[$KeyCol]...\n";
	  }
	  next TEXTLINE;
	}
      }
      if (exists($KeyToLinesMap{$KeyColValue})) {
	# Append to existing line...
	$Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
	$KeyToLinesMap{$KeyColValue} .= "\n" . $Line;
      }
      else {
	$Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
	$KeyToLinesMap{$KeyColValue} = $Line;
      }
    }
  }
  if ($OptionsInfo{Sort} =~ /^ascending$/i) {
    if ($OptionsInfo{KeyData} =~ /^alphanumeric$/i) {
      for $KeyColValue (sort { lc($a) cmp lc($b) } keys %KeyToLinesMap ) {
	print NEWTEXTFILE "$KeyToLinesMap{$KeyColValue}\n";
      }
    }
    else {
      for $KeyColValue (sort { $a <=> $b } keys %KeyToLinesMap ) {
	print NEWTEXTFILE "$KeyToLinesMap{$KeyColValue}\n";
      }
    }
  }
  else {
    if ($OptionsInfo{KeyData} =~ /^alphanumeric$/i) {
      for $KeyColValue (sort { lc($b) cmp lc($a) } keys %KeyToLinesMap ) {
	print NEWTEXTFILE "$KeyToLinesMap{$KeyColValue}\n";
      }
    }
    else {
      for $KeyColValue (sort { $b <=> $a } keys %KeyToLinesMap ) {
	print NEWTEXTFILE "$KeyToLinesMap{$KeyColValue}\n";
      }
    }
  }
  # Write out the lines with invalid data...
  if (@InvalidDataLines) {
    print "Placing ", scalar(@InvalidDataLines)," line(s) with invalid column key data at the end...\n";
    for $Line (@InvalidDataLines) {
      print NEWTEXTFILE "$Line\n";
    }
  }
  close NEWTEXTFILE;
  close TEXTFILE;

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
      $OutFileRoot = $FileName . "SortedByColumn";
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

# Make sure specified key column are okay...
sub ProcessColumnsInfo {
  my($Index, $TextFile, $SpecifiedKeyCol);

  @{$TextFilesInfo{KeyColNum}} = ();

  $SpecifiedKeyCol = $OptionsInfo{SpecifiedKeyCol};

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];
    $TextFilesInfo{KeyColNum}[$Index] = 0;

    if ($TextFilesInfo{FileOkay}[$Index]) {
      my($KeyColNum, $KeyColValid);

      $KeyColNum = 0;
      $KeyColValid = 1;
      if ($SpecifiedKeyCol) {
	if ($OptionsInfo{Mode} =~ /^colnum$/i) {
	  if ($SpecifiedKeyCol <= $TextFilesInfo{ColCount}[$Index]) {
	    $KeyColNum = $SpecifiedKeyCol - 1;
	  }
	  else {
	    $KeyColValid = 0;
	  }
	}
	else {
	  if (exists($TextFilesInfo{ColLabelToNumMap}[$Index]{$SpecifiedKeyCol})) {
	    $KeyColNum =  $TextFilesInfo{ColLabelToNumMap}[$Index]{$SpecifiedKeyCol};
	  }
	  else {
	    $KeyColValid = 0;
	  }
	}
      }
      if ($KeyColValid) {
	$TextFilesInfo{KeyColNum}[$Index] = $KeyColNum;
      }
      else {
	warn "Warning: Ignoring file $TextFile: Column key specified, $SpecifiedKeyCol, using \"k --key\" option doesn't exist\n";
	$TextFilesInfo{FileOkay}[$Index] = 0;
      }
    }
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{DetailLevel} = $Options{detail};

  $OptionsInfo{Sort} = $Options{sort};
  $OptionsInfo{KeyData} = $Options{keydata};

  $OptionsInfo{InDelim} = $Options{indelim};

  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /^tab$/i ) ? "\t" : (($Options{outdelim} =~ /^semicolon$/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^yes$/i) ? 1 : 0;

  $OptionsInfo{Overwrite} = defined $Options{overwrite} ? $Options{overwrite} : undef;
  $OptionsInfo{Root} = defined $Options{root} ? $Options{root} : undef;

  $OptionsInfo{Key} = defined $Options{key} ? $Options{key} : undef;
  $OptionsInfo{SpecifiedKeyCol} = "";

  if (defined $Options{key}) {
    $OptionsInfo{SpecifiedKeyCol} = $Options{key};
    if ($Options{mode} =~ /^colnum$/i) {
      if (!IsPositiveInteger($OptionsInfo{SpecifiedKeyCol})) {
	die "Error: Invalid value $Options{key} specified using \"-k --key\" option: Allowed values: > 0\n";
      }
    }
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {
  %Options = ();

  $Options{detail} = 1;
  $Options{mode} = "colnum";
  $Options{sort} = "ascending";
  $Options{keydata} = "numeric";
  $Options{indelim} = "comma";
  $Options{outdelim} = "comma";
  $Options{quote} = "yes";
  if (!GetOptions(\%Options, "detail|d=i", "help|h", "indelim=s", "key|k=s", "keydata=s", "mode|m=s", "outdelim=s", "overwrite|o", "quote|q=s", "root|r=s", "sort|s=s", "workingdir|w=s")) {
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
  if ($Options{keydata} !~ /^(numeric|alphanumeric)$/i) {
    die "Error: The value specified, $Options{keydata}, for option \"--keydata\" is not valid. Allowed values: numeric or alphanumeric\n";
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
  if ($Options{sort} !~ /^(ascending|descending)$/i) {
    die "Error: The value specified, $Options{sort}, for option \"-s --sort\" is not valid. Allowed values: ascending or descending\n";
  }
  if (!IsPositiveInteger($Options{detail})) {
    die "Error: The value specified, $Options{detail}, for option \"-d --detail\" is not valid. Allowed values: > 0\n";
  }
}

__END__

=head1 NAME

SortTextFiles.pl - Sort TextFile(s) using values for a column

=head1 SYNOPSIS

SortTextFiles.pl TextFile(s)...

SortTextFiles.pl [B<-d, --detail> infolevel] [B<-h, --help>] [B<--indelim> comma | semicolon] [B<-k, --key> colnum | collabel]
[B<--keydata> numeric | alphanumeric] [B<-m, --mode> colnum | collabel] [B<-o, --overwrite>]
[B<--outdelim> comma | tab | semicolon] [B<-q, --quote> yes | no] [B<-r, --root> rootname]
[B<-s, --sort> ascending | descending] [B<-w, --workingdir> dirname] TextFile(s)...

=head1 DESCRIPTION

Sort I<TextFile(s)> using values for a key column specified by a column number or label.
Only one column key can be specified for sorting. In an event of conflict during sorting
process, two similar values for a column key are simply transferred to output files in
order of their presence in input files. Additionally, rows with empty or inappropriate
values for column key are simply placed at the end. The file names are separated by space.
The valid file extensions are I<.csv> and I<.tsv> for comma/semicolon and tab delimited
text files respectively. All other file names are ignored. All the text files in a
current directory can be specified by I<*.csv>, I<*.tsv>, or the current directory
name. The B<--indelim> option determines the format of I<TextFile(s)>. Any file
which doesn't correspond to the format indicated by B<--indelim> option is ignored.

=head1 OPTIONS

=over 4

=item B<-d, --detail> I<infolevel>

Level of information to print about lines being ignored. Default: I<1>. Possible values:
I<1, 2 or 3>.

=item B<-h, --help>

Print this help message.

=item B<--indelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile(s)>. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a
delimiter.

=item B<-k, --key> I<col number | col name>

This value is mode specific. It specifies which column to use for sorting I<TextFile(s)>.
Possible values: I<col number or col label>. Default value: I<first column>.

=item B<--keydata> I<numeric | alphanumeric>

Data type for column key. Possible values: I<numeric or alphanumeric>. Default value:
I<numeric>. For I<alphanumeric> data values, comparison is case insensitive.

=item B<-m, --mode> I<colnum | collabel>

Specify how to sort text files: using column number or column label.
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

New text file name is generated using the root: <Root>.<Ext>. Default new file
name: <InitialTextFileName>SortedByColumn.<Ext>. The csv, and tsv
<Ext> values are used for comma/semicolon, and tab delimited text files
respectively. This option is ignored for multiple input files.

=item B<-s, --sort> I<ascending | descending>

Sorting order for column values. Possible values: I<ascending or descending>.
Default value: I<ascending>.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To perform numerical sort in ascending order using first column values and generate
a new CSV text file NewSample1.csv, type:

    % SortTextFiles.pl -o -r NewSample1 Sample1.csv

To perform numerical sort in descending order using MolWeight column  and generate
a new CSV text file NewSample1.csv, type:

    % SortTextFiles.pl -m collabel -k MolWeight --keydata numeric
      -s descending -r NewSample1 -o Sample1.csv

To perform numerical sort in ascending order using column number 1 and generate
a new TSV text file NewSample1.csv, type:

    % SortTextFiles.pl -m colnum -k 1 --keydata numeric -s ascending
      -r NewSample1 --outdelim tab -o Sample1.csv

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
