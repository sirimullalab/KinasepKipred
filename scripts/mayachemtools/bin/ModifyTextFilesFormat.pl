#!/usr/bin/perl -w
#
# File: ModifyTextFilesFormat.pl
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
    ModifyTextFileFormat($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Modify text file format...
sub ModifyTextFileFormat {
  my($Index) = @_;
  my($TextFile, $NewTextFile, $InDelim, $Line, @Words);

  $TextFile = $TextFilesList[$Index];

  $NewTextFile = $TextFilesInfo{OutFile}[$Index];
  $InDelim = $TextFilesInfo{InDelim}[$Index];

  print "Generating new $NewTextFile file...\n";

  open NEWTEXTFILE, ">$NewTextFile" or die "Error: Can't open $NewTextFile !$ \n";
  open TEXTFILE, "$TextFile" or die "Error: Can't open $TextFile: $! \n";

  while ($Line = GetTextLine(\*TEXTFILE)) {
    @Words = quotewords($InDelim, 0, $Line);
    $Line = JoinWords(\@Words, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print NEWTEXTFILE "$Line\n";
  }

  close NEWTEXTFILE;
  close TEXTFILE;
}

# Retrieve information about input text files...
sub RetrieveTextFilesInfo {
  my($Index, $TextFile, $NewTextFile, $FileDir, $FileName, $FileExt, $InDelim, $OutFileExt);

  %TextFilesInfo = ();
  @{$TextFilesInfo{FileOkay}} = ();
  @{$TextFilesInfo{InDelim}} = ();
  @{$TextFilesInfo{OutFile}} = ();

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];

    $TextFilesInfo{FileOkay}[$Index] = 0;
    $TextFilesInfo{InDelim}[$Index] = "";
    $TextFilesInfo{OutFile}[$Index] = "";

    if (!(-e $TextFile)) {
      warn "Warning: Ignoring file $TextFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($TextFile, "csv tsv")) {
      warn "Warning: Ignoring file $TextFile: It's not a csv or tsv file\n";
      next FILELIST;
    }
    if (!open TEXTFILE, "$TextFile") {
      warn "Warning: Ignoring file $TextFile: Couldn't open it: $! \n";
      next FILELIST;
    }
    close TEXTFILE;

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

    if (lc($InDelim) eq lc($OptionsInfo{OutDelim})) {
      warn "Warning: Ignoring file $TextFile: The value specified, $Options{outdelim}, for option \"--outdelim\" is same as input delimiter\n";
      next FILELIST;
    }

    $OutFileExt = ($Options{outdelim} =~ /^tab$/i ) ? "tsv" : "csv";

    $NewTextFile = $FileName;
    if ($OptionsInfo{OutFileRoot} && (@TextFilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($OptionsInfo{OutFileRoot});
      if ($RootFileName && $RootFileExt) {
	$NewTextFile = $RootFileName;
      }
      else {
	$NewTextFile = $OptionsInfo{OutFileRoot};
      }
      $NewTextFile .= ".$OutFileExt";
    }
    else {
      $NewTextFile .= "FormatModified" . ".$OutFileExt";
    }

    if (!$OptionsInfo{Overwrite}) {
      if (-e $NewTextFile) {
	warn "Warning: Ignoring file $TextFile: New Text file, $NewTextFile, already exists\n";
	next FILELIST;
      }
    }

    $TextFilesInfo{FileOkay}[$Index] = 1;
    $TextFilesInfo{InDelim}[$Index] = $InDelim;
    $TextFilesInfo{OutFile}[$Index] = $NewTextFile;
  }

}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{InDelim} = $Options{indelim};

  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : undef;
  $OptionsInfo{Overwrite} = $Options{overwrite} ? $Options{overwrite} : undef;

  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /^tab$/i ) ? "\t" : (($Options{outdelim} =~ /^semicolon$/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^yes$/i) ? 1 : 0;

}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{indelim} = "comma";
  $Options{outdelim} = "tab";
  $Options{quote} = "yes";

  if (!GetOptions(\%Options, "help|h", "indelim=s", "outdelim=s", "overwrite|o", "quote|q=s", "root|r=s", "workingdir|w=s")) {
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
}

__END__

=head1 NAME

ModifyTextFilesFormat.pl - Change CSV Textfile(s) into TSV Textfile(s) and vice versa

=head1 SYNOPSIS

ModifyTextFilesFormat.pl TextFile(s)...

ModifyTextFilesFormat.pl [B<-h, --help>] [B<--indelim> comma | semicolon]
[B<--outdelim> comma | tab | semicolon] [B<-q, --quote> yes | no] [B<-r, --root> rootname]
[B<-w, --workingdir> dirname] TextFile(s)...

=head1 DESCRIPTION

Interchange CSV and TSV I<TextFile(s)> format. Mutiple file names are separated by spaces.
The valid file extensions are I<.csv> and I<.tsv> for comma/semicolon and tab delimited
text files respectively. All other file names are ignored. All the text files in a current
directory can be specified by I<*.csv>, I<*.tsv>, or the current directory name. The
B<--indelim> option determines the format of I<TextFile(s)>. Any file which doesn't
correspond to the format indicated by B<--indelim> option is ignored.

=head1 OPTIONS

=over 4

=item B<-h, --help>

Print this help message.

=item B<--indelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile(s)>. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a
delimiter.

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
name: <InitialTextFileName>FormatModified.<Ext>. The csv, and tsv
<Ext> values are used for comma/semicolon, and tab delimited text files
respectively. This option is ignored for multiple input files.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To convert Sample*.csv into TSV files, type:

    % ModifyTextFilesFormat.pl --outdelim tab -q no -o Sample*.csv

To convert Sample1.tsv into NewSample1.csv without any quotes around column
data values, type:

    % ModifyTextFilesFormat.pl --outdelim comma - q no
      -r NewSample1 -o Sample1.tsv

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

ModifyNewLineChar.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
