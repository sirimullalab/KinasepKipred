#!/usr/bin/perl -w
#
# File: JoinTextFiles.pl
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

if (@TextFilesList == 1) {
  die "Error: Specify more than one Text file.\n";
}

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Setup information about input files...
print "Checking input text files...\n";
my(%TextFilesInfo);
RetrieveTextFilesInfo();

# Join files...
print "\nGenerating new text file $OptionsInfo{NewTextFile}...\n";
JoinTextFiles();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Join all valid Text files...
sub JoinTextFiles {
  my($FileIndex, $TextFile, $NewTextFile, $Line, $FirstColLabelsLine, $OutDelim, $OutQuote, $InDelim, @Words, @ColLabels);

  $NewTextFile = $OptionsInfo{NewTextFile};

  $FirstColLabelsLine = '';

  $OutDelim = $OptionsInfo{OutDelim}; $OutQuote = $OptionsInfo{OutQuote};

  open NEWTEXTFILE, ">$NewTextFile" or die "Error: Couldn't open $NewTextFile: $! \n";
  FILELIST: for $FileIndex (0 .. $#TextFilesList) {
    if (!$TextFilesInfo{FileOkay}[$FileIndex]) {
      next FILELIST;
    }

    $TextFile = $TextFilesList[$FileIndex];
    $InDelim = $TextFilesInfo{InDelim}[$FileIndex];

    print "\nProcessing file $TextFile...\n";

    open TEXTFILE, "$TextFile" or die "Error: Couldn't open $TextFile: $! \n";

    if ($OptionsInfo{Label}) {
      if ($OptionsInfo{Fast}) {
	$Line = GetTextLine(\*TEXTFILE);
	if (!$FirstColLabelsLine) {
	  $FirstColLabelsLine = $Line;
	  print NEWTEXTFILE "$FirstColLabelsLine\n";
	}
      }
      else {
	$Line = GetTextLine(\*TEXTFILE);
	if (!$FirstColLabelsLine) {
	  @ColLabels = quotewords($InDelim, 0, $Line);
	  $FirstColLabelsLine = JoinWords(\@ColLabels, $OutDelim, $OutQuote);
	  print NEWTEXTFILE "$FirstColLabelsLine\n";
	}
      }
    }

    while ($Line = GetTextLine(\*TEXTFILE)) {
      if (!$OptionsInfo{Fast}) {
	@Words = quotewords($InDelim, 0, $Line);
	$Line = JoinWords(\@Words, $OutDelim, $OutQuote);
      }
      print NEWTEXTFILE "$Line\n";
    }
    close TEXTFILE;
  }

  close NEWTEXTFILE;
}

# Retrieve information about Text files...
sub RetrieveTextFilesInfo {
  my($Index, $TextFile, $InDelim, $Line, $FirstColLabelsLine, $ColLabelsLine, $FileDir, $FileName, $FileExt, @FirstColLabels, @ColLabels);

  %TextFilesInfo = ();
  @{$TextFilesInfo{FileOkay}} = ();
  @{$TextFilesInfo{InDelim}} = ();

  $FirstColLabelsLine = ''; $ColLabelsLine = '';
  @FirstColLabels = (); @ColLabels = ();

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFilesInfo{FileOkay}[$Index] = 0;
    $TextFilesInfo{InDelim}[$Index] = "";

    $TextFile = $TextFilesList[$Index];
    if (!(-e $TextFile)) {
      warn "Warning: Ignoring file $TextFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($TextFile, "csv tsv")) {
      warn "Warning: Ignoring file $TextFile: It's not a Text file\n";
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

    if (! open TEXTFILE, "$TextFile") {
      warn "Warning: Ignoring file $TextFile: Couldn't open it: $! \n";
      next FILELIST;
    }
    $Line = GetTextLine(\*TEXTFILE);
    close TEXTFILE;

    if ($OptionsInfo{Label}) {
      if (!$OptionsInfo{Fast}) {
	@ColLabels = quotewords($InDelim, 0, $Line);
	if ($FirstColLabelsLine) {
	  if (@ColLabels != @FirstColLabels) {
	    warn "Warning: Ignoring file $TextFile: The number of columns in this file, ", scalar(@ColLabels), ", is different from the number of columns, ", scalar(@FirstColLabels), ", in the first valid text file. \n";
	    next FILELIST;
	  }
	  $ColLabelsLine = JoinWords(\@ColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
	  if ($ColLabelsLine ne $FirstColLabelsLine) {
	    warn "Warning: Ignoring file $TextFile: The column names in this file are different from those in first valid text file.\nColumnlabels in first valid text file: $FirstColLabelsLine \nColumnlabels in current text file: $ColLabelsLine\n";
	    next FILELIST;
	  }
	}
	else {
	  @FirstColLabels = @ColLabels;
	  $FirstColLabelsLine = JoinWords(\@FirstColLabels, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
	}
      }
    }

    $TextFilesInfo{FileOkay}[$Index] = 1;
    $TextFilesInfo{InDelim}[$Index] = $InDelim;
  }
}

# Process option values...
sub ProcessOptions {
  my($FileDir, $FileName, $FileExt, $NewTextFile, $Index);

  %OptionsInfo = ();

  $OptionsInfo{Fast} = $Options{fast} ? $Options{fast} : 0;

  $OptionsInfo{InDelim} = $Options{indelim};
  $OptionsInfo{Label} = ($Options{label} =~ /^yes$/i) ? 1 : 0;

  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : undef;
  $OptionsInfo{Overwrite} = $Options{overwrite} ? $Options{overwrite} : undef;

  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /^tab$/i ) ? "\t" : (($Options{outdelim} =~ /^semicolon$/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^yes$/i) ? 1 : 0;

  if ($Options{root}) {
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($Options{root});
    if ($FileName && $FileExt) {
      $NewTextFile = $FileName;
    }
    else {
      $NewTextFile =  $Options{root};
    }
  }
  else {
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($TextFilesList[0]);
    $NewTextFile = $FileName . "1To" . @TextFilesList . "Joined";
  }

  if ($Options{outdelim} =~ /^tab$/i) {
    $NewTextFile .= ".tsv";
  }
  else {
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
  $Options{label} = "yes";
  $Options{indelim} = "comma";
  $Options{outdelim} = "comma";
  $Options{quote} = "yes";
  if (!GetOptions(\%Options, "fast|f", "help|h", "indelim=s", "label|l=s", "outdelim=s", "overwrite|o", "quote|q=s", "root|r=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir},  for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
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

JoinTextFiles.pl - Join multiple CSV or TSV text files into a single text file

=head1 SYNOPSIS

JoinTextFiles.pl TextFiles...

JoinTextFiles.pl  [B<-f, --fast>] [B<-h, --help>] [B<--indelim> comma | semicolon]
[B<-l, --label> yes | no] [B<-o, --overwrite>] [B<--outdelim> comma | tab | semicolon]
[B<-q, --quote> yes | no] [B<-r, --root> rootname] [B<-w, --workingdir> dirname] TextFiles...

=head1 DESCRIPTION

Multiple CSV or TSV I<TextFiles> are joined to generate a single text file. The
file names are separated by spaces. The valid file extensions are I<.csv> and I<.tsv>
for comma/semicolon and tab delimited text files respectively. All other file names
are ignored. All the text files in a current directory can be specified by I<*.csv>,
I<*.tsv>, or the current directory name. The B<--indelim> option determines the
format of I<TextFiles>. Any file which doesn't correspond to the format indicated
by B<--indelim> option is ignored.

=head1 OPTIONS

=over 4

=item B<-f, --fast>

In this mode, B<--indelim> and B<-q --quote> options are ignored. The format of
input and output file(s) are assumed to be similar. And the text lines from I<TextFiles>
are simply transferred to output file without any processing.

=item B<-h, --help>

Print this help message.

=item B<--indelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile(s)>. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a
delimiter.

=item B<-l, --label> I<yes | no>

First line contains column labels. Possible values: I<yes or no>. Default value: I<yes>.

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
name: <FirstTextFileName>1To<Count>Joined.<Ext>. The csv, and tsv
<Ext> values are used for comma/semicolon, and tab delimited text files
respectively.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To join CSV text files, type:

    % JoinTextFiles.pl -o Sample1.csv Sample2.csv
    % JoinTextFiles.pl -o *.csv

To join Sample*.tsv TSV text files into a NewSample.tsv file, type:

    % JoinTextFiles.pl -o -r NewSample Sample*.tsv

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

MergeTextFiles.pl, ModifyTextFilesFormat.pl, SplitTextFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
