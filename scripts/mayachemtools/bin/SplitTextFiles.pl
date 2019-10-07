#!/usr/bin/perl -w
#
# File: SplitTextFiles.pl
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

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

print "Checking input text file(s)...\n";
my(%TextFilesInfo);
RetrieveTextFilesInfo();

# Generate output files...
my($FileIndex);
if (@TextFilesList > 1) {
  print "\nProcessing text files...\n";
}
for $FileIndex (0 .. $#TextFilesList) {
  if ($TextFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $TextFilesList[$FileIndex]...\n";
    SplitTextFile($FileIndex);
  }
}

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Split a Text file...
#
sub SplitTextFile {
  my($FileIndex) = @_;
  my($TextFile, $LineCount, $MaxLinesPerFile, $MaxNumOfFiles);

  $TextFile = $TextFilesList[$FileIndex];

  if (!open TEXTFILE, "$TextFile") {
    warn "Warning: Ignoring file $TextFile: Couldn't open it: $! \n";
    return;
  }

  $MaxNumOfFiles = $OptionsInfo{NumOfFiles};

  # Count number of lines to figure out maximum number of lines per file...
  $LineCount = 0;
  while (<TEXTFILE>) {
      $LineCount++;
  }
  close TEXTFILE;

  if ($LineCount < $MaxNumOfFiles) {
    warn "Warning: Ignoring file $TextFile: Total number of lines, $LineCount, is smaller than\nnumber of new files, $MaxNumOfFiles\n";
    return;
  }

  $MaxLinesPerFile = int $LineCount / $MaxNumOfFiles;

  GenerateTextFiles($FileIndex, $MaxNumOfFiles, $MaxLinesPerFile);
}

# Generate new Text files...
#
sub GenerateTextFiles {
  my($FileIndex, $NumOfFiles, $NumOfLinesPerFile) = @_;
  my($TextFile, $LineCount, $NewFileIndex, $NewFileName, $MaxLinesCount, $InDelim, $OutDelim, $OutQuote, $ColLabelsLine, $Line, @ColLabels, @Words, @NewTextFilesList);

  # Setup new file names list...
  @NewTextFilesList = ();
  for $NewFileIndex (1 .. $NumOfFiles) {
    $NewFileName = $TextFilesInfo{OutFileRoot}[$FileIndex] . "Part${NewFileIndex}." . $TextFilesInfo{OutFileExt}[$FileIndex];
    if (!$OptionsInfo{OverwriteFiles}) {
      if (-e $NewFileName) {
	warn "Warning: Ignoring file $TextFile: New Text file, $NewFileName, already exists\n";
	return;
      }
    }
    push @NewTextFilesList, $NewFileName;
  }

  $TextFile = $TextFilesList[$FileIndex];

  if (!open TEXTFILE, "$TextFile") {
    warn "Warning: Ignoring file $TextFile: Couldn't open it: $! \n";
    return;
  }

  $InDelim = $TextFilesInfo{InDelim}[$FileIndex];

  $OutDelim = $OptionsInfo{OutDelim};
  $OutQuote = $OptionsInfo{OutQuote};

  $MaxLinesCount = $NumOfLinesPerFile;
  $LineCount = 0;
  $NewFileIndex = 1;

  open NEWTEXTFILE, ">$NewTextFilesList[$NewFileIndex - 1]" or die "Error: Can't open $NewTextFilesList[$NewFileIndex -1]: $! \n";
  print "Generating $NewTextFilesList[$NewFileIndex - 1] file...\n";

  if ($OptionsInfo{Label}) {
    if ($OptionsInfo{Fast}) {
      $ColLabelsLine = GetTextLine(\*TEXTFILE);
    }
    else {
      $Line = GetTextLine(\*TEXTFILE);
      @ColLabels = quotewords($InDelim, 0, $Line);
      $ColLabelsLine = JoinWords(\@ColLabels, $OutDelim, $OutQuote);
    }
    print NEWTEXTFILE "$ColLabelsLine\n";
  }

  while ($Line = GetTextLine(\*TEXTFILE)) {
    $LineCount++;

    if (!$Options{fast}) {
      @Words = quotewords($InDelim, 0, $Line);
      $Line = JoinWords(\@Words, $OutDelim, $OutQuote);
    }
    print NEWTEXTFILE "$Line\n";

    if ($NewFileIndex <= $NumOfFiles) {
      if ($LineCount >= $MaxLinesCount) {
	if ($NewFileIndex < $NumOfFiles) {
	  close NEWTEXTFILE;
	}
	$NewFileIndex++;
	$MaxLinesCount = $NumOfLinesPerFile * $NewFileIndex;

	if ($NewFileIndex <= $NumOfFiles) {
	  open NEWTEXTFILE, ">$NewTextFilesList[$NewFileIndex - 1]" or die "Error: Can't open $NewTextFilesList[$NewFileIndex -1]: $! \n";
	  print "Generating $NewTextFilesList[$NewFileIndex - 1] file...\n";

	  if ($OptionsInfo{Label}) {
	    print NEWTEXTFILE "$ColLabelsLine\n";
	  }
	}
      }
    }
  }
  close NEWTEXTFILE;
  close TEXTFILE;
}

# Retrieve information about Text files...
sub RetrieveTextFilesInfo {
  my($Index, $TextFile, $InDelim, $FileDir, $FileName, $FileExt, $OutFileRoot, $OutFileExt);

  %TextFilesInfo = ();
  @{$TextFilesInfo{FileOkay}} = ();
  @{$TextFilesInfo{InDelim}} = ();
  @{$TextFilesInfo{OutFileRoot}} = ();
  @{$TextFilesInfo{OutFileExt}} = ();

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFilesInfo{FileOkay}[$Index] = 0;
    $TextFilesInfo{InDelim}[$Index] = "";
    $TextFilesInfo{OutFileRoot}[$Index] = "";
    $TextFilesInfo{OutFileExt}[$Index] = "";

    $TextFile = $TextFilesList[$Index];
    if (!(-e $TextFile)) {
      warn "Warning: Ignoring file $TextFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($TextFile, "csv tsv")) {
      warn "Warning: Ignoring file $TextFile: It's not a Text file\n";
      next FILELIST;
    }
    if (! open TEXTFILE, "$TextFile") {
      warn "Warning: Ignoring file $TextFile: Couldn't open it: $! \n";
      next FILELIST;
    }
    close TEXTFILE;

    ($FileDir, $FileName, $FileExt) = ParseFileName($TextFile);

    # Setup input delimiter...
    $InDelim = '';
    if (!$OptionsInfo{Fast}) {
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
    }

    # Setup output file root...
    $OutFileExt = $OptionsInfo{Fast} ? $FileExt : (($Options{outdelim} =~ /^tab$/i ) ? "tsv" : "csv");

    if ($OptionsInfo{OutFileRoot} && (@TextFilesList == 1)) {
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
      $OutFileRoot = $FileName;
    }

    $TextFilesInfo{FileOkay}[$Index] = 1;
    $TextFilesInfo{InDelim}[$Index] = $InDelim;
    $TextFilesInfo{OutFileRoot}[$Index] = $OutFileRoot;
    $TextFilesInfo{OutFileExt}[$Index] = $OutFileExt;
  }
}

# Process option values...
sub ProcessOptions {

  %OptionsInfo = ();

  $OptionsInfo{Fast} = defined $Options{fast} ? $Options{fast} : undef;

  $OptionsInfo{InDelim} = $Options{indelim};
  $OptionsInfo{Label} = ($Options{label} =~ /^yes$/i) ? 1 : 0;

  $OptionsInfo{NumOfFiles} = $Options{numfiles};

  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : undef;
  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? $Options{overwrite} : undef;

  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /^tab$/i ) ? "\t" : (($Options{outdelim} =~ /^semicolon$/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^yes$/i) ? 1 : 0;
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{label} = "yes";
  $Options{numfiles} = 2;
  $Options{indelim} = "comma";
  $Options{outdelim} = "comma";
  $Options{quote} = "yes";
  if (!GetOptions(\%Options, "fast|f", "help|h", "indelim=s", "label|l=s", "numfiles|n=i", "outdelim=s", "overwrite|o", "quote|q=s", "root|r=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir},  for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{numfiles} < 2) {
    die "Error: The value specified, $Options{numfiles},  for option \"-n --numfiles\" is not valid. Allowed values: >= 2 \n";
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
  if ($Options{label} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{label}, for option \"-l --label\" is not valid. Allowed values: yes or no\n";
  }
}

__END__

=head1 NAME

SplitTextFiles.pl - Split CSV or TSV TextFile(s) into multiple text files

=head1 SYNOPSIS

SplitTextFiles.pl TextFile(s)...

SplitTextFiles.pl [B<-f, --fast>] [B<-h, --help>] [B<--indelim> comma | semicolon]
[B<-l, --label> yes | no] [B<-n, --numfiles> number] [B<-o, --overwrite>]
[B<--outdelim> comma | tab | semicolon] [B<-q, --quote> yes | no]
[B<-r, --root> rootname] [B<-w, --workingdir> dirname] TextFile(s)...

=head1 DESCRIPTION

Split CSV or TSV I<TextFile(s)> into multiple text files. Each new text file contains
a subset of similar number of lines from the initial file. The file names are separated
by space. The valid file extensions are I<.csv> and I<.tsv> for comma/semicolon and
tab delimited text files respectively. All other file names are ignored. All the text
files in a current directory can be specified by I<*.csv>, I<*.tsv>, or the current
directory name. The B<--indelim> option determines the format of I<TextFile(s)>.
Any file which doesn't correspond to the format indicated by B<--indelim> option
is ignored.

=head1 OPTIONS

=over 4

=item B<-f, --fast>

In this mode, B<--indelim, --outdelim>, and B<-q --quote> options are ignored. The
format of input and output file(s) are assumed to be similar. And the text lines
from input I<TextFile(s)> are just transferred to output file(s) without any processing.

=item B<-h, --help>

Print this help message.

=item B<--indelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile(s)>. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a
delimiter.

=item B<-l, --label> I<yes | no>

First line contains column labels. Possible values: I<yes or no>. Default value: I<yes>.

=item B<-n, --numfiles> I<number>

Number of new files to generate for each TextFile(s). Default: I<2>.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<--outdelim> I<comma | tab | semicolon>

Output text file delimiter. Possible values: I<comma, tab, or semicolon>.
Default value: I<comma>

=item B<-q, --quote> I<yes | no>

Put quotes around column values in output text file. Possible values: I<yes or
no>. Default value: I<yes>.

=item B<-r, --root> I<rootname>

New text file names are generated using the root: <Root>Part<Count>.<Ext>.
Default new file names: <InitialTextFileName>Part<Count>.<Ext>. The csv, and tsv
<Ext> values are used for comma/semicolon, and tab delimited text files
respectively.This option is ignored for multiple input files.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To split each CSV text files into 4 different text files type:

    % SplitTextFiles.pl -n 5 -o Sample1.csv Sample2.csv
    % SplitTextFiles.pl -n 5 -o *.csv

To split Sample1.tsv into 10 different CSV text files, type:

    % SplitTextFiles.pl -n 10 --outdelim comma -o Sample1.tsv

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

JoinTextFiles.pl, MergeTextFiles.pl, ModifyTextFilesFormat.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
