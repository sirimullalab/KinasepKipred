#!/usr/bin/perl -w
#
# File: TextFilesToSDFiles.pl
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
use SDFileUtil;

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
    ConvertTextFile($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Convert text file to SD file...
sub ConvertTextFile {
  my($Index) = @_;
  my($TextFile, $SDFile, $Line, $InDelim, $Label, $Value, $ColIndex, $ColCount, @ColLabels, @LineWords);

  $TextFile = $TextFilesList[$Index];
  $InDelim = $TextFilesInfo{InDelim}[$Index];
  $SDFile = $TextFilesInfo{OutSDFile}[$Index];
  @ColLabels = @{$TextFilesInfo{ColLabels}[$Index]};
  $ColCount = @ColLabels;

  print "Generating SD file $SDFile...\n";
  open SDFILE, ">$SDFile" or die "Error: Couldn't open $SDFile: $! \n";
  open TEXTFILE, "$TextFile" or die "Error: Can't open $TextFile: $! \n";
  if ($OptionsInfo{ColLabelsPresent}) {
    # Skip over column labels from old file...
    $Line = GetTextLine(\*TEXTFILE);
  }
  my($Date) = GenerateMiscLineDateStamp();
  while ($Line = GetTextLine(\*TEXTFILE)) {
    @LineWords = quotewords($InDelim, 0, $Line);

    # Write out empty CTAB block...
    print SDFILE GenerateEmptyCtabBlockLines($Date), "\n";

    # Write out data fields and values...
    for $ColIndex (0 .. $#LineWords) {
      if ($ColIndex < $ColCount) {
	$Label = $ColLabels[$ColIndex];
	$Value = $LineWords[$ColIndex];
	print SDFILE "> <$Label>\n$Value\n\n";
      }
    }
    print SDFILE "\$\$\$\$\n";
  }
  close SDFILE;
  close TEXTFILE;
}

# Retrieve information about input text files...
sub RetrieveTextFilesInfo {
  my($Index, $TextFile, $FileDir, $FileName, $FileExt, $InDelim, $Line, @LineWords, @ColLabels, $OutFileRoot,  $OutFile, $ColNum, $ColLabel);

  %TextFilesInfo = ();

  @{$TextFilesInfo{FileOkay}} = ();
  @{$TextFilesInfo{ColCount}} = ();
  @{$TextFilesInfo{ColLabels}} = ();
  @{$TextFilesInfo{InDelim}} = ();
  @{$TextFilesInfo{OutSDFile}} = ();


  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];

    $TextFilesInfo{FileOkay}[$Index] = 0;
    $TextFilesInfo{ColCount}[$Index] = 0;
    $TextFilesInfo{InDelim}[$Index] = "";
    $TextFilesInfo{OutSDFile}[$Index] = "";

    @{$TextFilesInfo{ColLabels}[$Index]} = ();

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
    @LineWords = quotewords($InDelim, 0, $Line);
    @ColLabels = ();
    if ($OptionsInfo{ColLabelsPresent}) {
      push @ColLabels, @LineWords;
    }
    else {
      for $ColNum (1 .. @LineWords) {
	$ColLabel = "Column${ColNum}Data";
	push @ColLabels, $ColLabel;
      }
    }
    close TEXTFILE;

    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($TextFile);
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
      $OutFileRoot = "${FileName}WithNoStrData";
    }

    $OutFile = "${OutFileRoot}.sdf";
    if (!$Options{overwrite}) {
      if (-e $OutFile) {
	warn "Warning: Ignoring file $TextFile: The file $OutFile already exists\n";
	next FILELIST;
      }
    }
    $TextFilesInfo{FileOkay}[$Index] = 1;
    $TextFilesInfo{InDelim}[$Index] = $InDelim;
    $TextFilesInfo{OutSDFile}[$Index] = "$OutFile";

    $TextFilesInfo{ColCount}[$Index] = @ColLabels;
    push @{$TextFilesInfo{ColLabels}[$Index]}, @ColLabels;
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Label} = $Options{label};
  $OptionsInfo{ColLabelsPresent} = ($Options{label} =~ /^yes$/i) ? 1 : 0;

  $OptionsInfo{InDelim} = $Options{indelim};
  $OptionsInfo{Overwrite} = defined $Options{overwrite} ? $Options{overwrite} : undef;

  $OptionsInfo{OutFileRoot} = defined $Options{root} ? $Options{root} : undef;

}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{label} = "yes";
  $Options{indelim} = "comma";
  if (!GetOptions(\%Options, "help|h", "indelim=s", "label|l=s", "overwrite|o", "root|r=s", "workingdir|w=s")) {
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
  if ($Options{label} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{label}, for option \"-l --label\" is not valid. Allowed values: yes or no\n";
  }
}

__END__

=head1 NAME

TextFilesToSDFiles.pl - Generate SD files from CSV or TSV TextFile(s)

=head1 SYNOPSIS

TextFilesToSDFiles.pl TextFile(s)...

TextFilesToSDFiles.pl [B<-h, --help>] [B<--indelim> comma | semicolon] [B<-l, --label> yes | no]
[B<-o, --overwrite>] [B<-r, --root> rootname] [B<-w, --workingdir> dirname] TextFile(s)...

=head1 DESCRIPTION

Generate SD files from CSV or TSV I<TextFile(s)>. The new SD files contain no structure data
as indicated by empty structure data block; Data fields and values in SD files are generated
using I<TextFile(s)> column labels and corresponding data values.

Multiple I<TextFile(s)> names are separated by space. The valid file extensions are
I<.csv> and I<.tsv> for comma/semicolon and tab delimited text files respectively.
All other file names are ignored. All the text files in a current directory can be specified
by I<*.csv>, I<*.tsv>, or the current directory name. The B<--indelim> option
determines the format of I<TextFile(s)>. Any file which doesn't correspond to
the format indicated by B<--indelim> option is ignored.

=head1 OPTIONS

=over 4

=item B<-h, --help>

Print this help message.

=item B<--indelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile(s)>. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a
delimiter.

=item B<-l, --label> I<yes | no>

First line contains column labels. Possible values: I<yes or no>. Default value: I<yes>.
Column labels are used to create SD data field labels; otherwise, data field labels look
like Column<colnumber>Data.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

New SD file names are generated using the root: <Root>.sdf. Default new file names:
<TextFileName>WithNoStrData.sdf. This option is ignored for multiple input files.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To generate NewSample1.sdf file from Sample1.csv file, type:

    % TextFilesToSDFiles.pl -o -r NewSample1 Sample1.csv

To generate NewSample1.sdf file from Sample1.tsv file which doesn't
contain column labels line, type:

    % TextFilesToSDFiles.pl --label no -o -r NewSample1 Sample1.tsv

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

JoinTextFiles.pl, MergeTextFiles.pl, ModifySDFilesDataFields.pl, ModifyTextFilesFormat.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
