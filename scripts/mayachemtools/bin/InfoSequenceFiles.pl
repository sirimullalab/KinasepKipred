#!/usr/bin/perl -w
#
# File: InfoSequenceFiles.pl
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
use SequenceFileUtil;
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

my(@SequenceFilesList);
@SequenceFilesList = ExpandFileNames(\@ARGV, "aln msf fasta fta pir");

print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

print "Checking input sequence file(s)...\n";
my(%SequenceFilesInfo);
RetrieveSequenceFilesInfo();

my($FileIndex);
if (@SequenceFilesList > 1) {
  print "\nProcessing sequence files...\n";
}
for $FileIndex (0 .. $#SequenceFilesList) {
  if ($SequenceFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $SequenceFilesList[$FileIndex]...\n";
    ListSequenceFileInfo($FileIndex);
  }
}
ListTotalSizeOfFiles();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# List appropriate information...
sub ListSequenceFileInfo {
  my($Index) = @_;
  my($SequenceFile, $SequenceDataRef);

  $SequenceFile = $SequenceFilesList[$Index];

  $SequenceDataRef = ReadSequenceFile($SequenceFile);

  my($SequencesCount) = $SequenceDataRef->{Count};
  print "\nNumber of sequences: $SequencesCount\n";

  if ($OptionsInfo{ListShortestSequence} && ($SequencesCount > 1)) {
    my($ShortestSeqID, $ShortestSeq, $ShortestSeqLen, $Description) = GetShortestSequence($SequenceDataRef, $OptionsInfo{IgnoreGaps});
    print "\nShortest sequence information:\nID: $ShortestSeqID; Length:$ShortestSeqLen\n";
    if ($OptionsInfo{DetailLevel} >= 2) {
      print "Description: $Description\n";
    }
    if ($OptionsInfo{DetailLevel} >= 3) {
      print "Sequence: $ShortestSeq\n";
    }
  }
  if ($OptionsInfo{ListLongestSequence} && ($SequencesCount > 1)) {
    my($LongestSeqID, $LongestSeq, $LongestSeqLen, $Description) = GetLongestSequence($SequenceDataRef, $OptionsInfo{IgnoreGaps});
    print "\nLongest sequence information:\nID: $LongestSeqID; Length: $LongestSeqLen\n";
    if ($OptionsInfo{DetailLevel} >= 2) {
      print "Description: $Description\n";
    }
    if ($OptionsInfo{DetailLevel} >= 3) {
      print "Sequence: $LongestSeq\n";
    }
  }
  if ($OptionsInfo{FrequencyAnalysis} && ($SequencesCount > 1)) {
    PerformLengthFrequencyAnalysis($SequenceDataRef);
  }
  if ($OptionsInfo{ListSequenceLengths}) {
    ListSequenceLengths($SequenceDataRef);
  }

  # File size and modification information...
  print "\nFile size: ", FormatFileSize($SequenceFilesInfo{FileSize}[$Index]), " \n";
  print "Last modified: ", $SequenceFilesInfo{FileLastModified}[$Index], " \n";
}

# List information about sequence lengths...
sub ListSequenceLengths {
  my($SequenceDataRef) = @_;
  my($ID, $SeqLen, $Sequence, $Description);

  print "\nSequence lengths information:\n";
  for $ID (@{$SequenceDataRef->{IDs}}) {
    $Sequence = $SequenceDataRef->{Sequence}{$ID};
    $Description = $SequenceDataRef->{Description}{$ID};
    $SeqLen = GetSequenceLength($Sequence, $OptionsInfo{IgnoreGaps});
    if ($OptionsInfo{IgnoreGaps}) {
      $Sequence = RemoveSequenceGaps($Sequence);
    }
    print "ID: $ID; Length:$SeqLen\n";
    if ($OptionsInfo{DetailLevel} >= 2) {
      print "Description: $Description\n";
    }
    if ($OptionsInfo{DetailLevel} >= 3) {
      print "Sequence: $Sequence\n";
    }
    if ($OptionsInfo{DetailLevel} >= 2) {
      print "\n";
    }
  }
}

# Total size of all the fiels...
sub ListTotalSizeOfFiles {
  my($FileOkayCount, $TotalSize, $Index);

  $FileOkayCount = 0;
  $TotalSize = 0;

  for $Index (0 .. $#SequenceFilesList) {
    if ($SequenceFilesInfo{FileOkay}[$Index]) {
      $FileOkayCount++;
      $TotalSize += $SequenceFilesInfo{FileSize}[$Index];
    }
  }
  if ($FileOkayCount > 1) {
    print "\nTotal size of $FileOkayCount files: ", FormatFileSize($TotalSize), "\n";
  }
}


# Perform frequency analysis of sequence lengths
sub PerformLengthFrequencyAnalysis {
  my($SequenceDataRef, $SequenceLengthsRef) = @_;
  my ($ID, $SeqLen, $Sequence, $SequenceLenBin, $LenBin, $SequenceLenCount, @SequenceLengths, %SequenceLenFrequency);

  @SequenceLengths = ();
  %SequenceLenFrequency = ();
  for $ID (@{$SequenceDataRef->{IDs}}) {
    $Sequence = $SequenceDataRef->{Sequence}{$ID};
    $SeqLen = GetSequenceLength($Sequence, $OptionsInfo{IgnoreGaps});
    push @SequenceLengths, $SeqLen;
  }
  if (@{$OptionsInfo{BinRange}}) {
    %SequenceLenFrequency = Frequency(\@SequenceLengths, \@{$OptionsInfo{BinRange}});
  }
  else {
    %SequenceLenFrequency = Frequency(\@SequenceLengths, $OptionsInfo{NumOfBins});
  }
  print "\nDistribution of sequence lengths (LengthBin => Count):\n";
  for $SequenceLenBin (sort { $a <=> $b} keys %SequenceLenFrequency) {
    $SequenceLenCount = $SequenceLenFrequency{$SequenceLenBin};
    $LenBin = sprintf("%.1f", $SequenceLenBin) + 0;
    print "$LenBin => $SequenceLenCount; ";
  }
  print "\n";
}

# Retrieve information about sequence files...
sub RetrieveSequenceFilesInfo {
  my($Index, $SequenceFile, $FileSupported, $FileFormat, $ModifiedTimeString, $ModifiedDateString);

  %SequenceFilesInfo = ();
  @{$SequenceFilesInfo{FileOkay}} = ();
  @{$SequenceFilesInfo{FileFormat}} = ();
  @{$SequenceFilesInfo{FileSize}} = ();
  @{$SequenceFilesInfo{FileLastModified}} = ();

  FILELIST: for $Index (0 .. $#SequenceFilesList) {
    $SequenceFile = $SequenceFilesList[$Index];

    if (! open SEQUENCEFILE, "$SequenceFile") {
      warn "Warning: Ignoring file $SequenceFile: Couldn't open it: $! \n";
      next FILELIST;
    }
    close SEQUENCEFILE;

    $SequenceFilesInfo{FileOkay}[$Index] = 0;
    $SequenceFilesInfo{FileFormat}[$Index] = 'NotSupported';
    $SequenceFilesInfo{FileSize}[$Index] = 0;
    $SequenceFilesInfo{FileLastModified}[$Index] = '';

    ($FileSupported, $FileFormat) = IsSupportedSequenceFile($SequenceFile);
    if (!$FileSupported) {
      warn "Warning: Ignoring file $SequenceFile: Sequence file format is not supported.\n";
      next FILELIST;
    }

    $SequenceFilesInfo{FileOkay}[$Index] = 1;
    $SequenceFilesInfo{FileFormat}[$Index] = $FileFormat;
    $SequenceFilesInfo{FileSize}[$Index] = FileSize($SequenceFile);

    ($ModifiedTimeString, $ModifiedDateString) = FormattedFileModificationTimeAndDate($SequenceFile);
    $SequenceFilesInfo{FileLastModified}[$Index] = "$ModifiedTimeString; $ModifiedDateString";
  }
}

# Process option values...
sub ProcessOptions {

  $OptionsInfo{All} = defined $Options{all} ? $Options{all} : undef;

  $OptionsInfo{Count} = defined $Options{count} ? $Options{count} : undef;
  $OptionsInfo{DetailLevel} = $Options{detail};
  $OptionsInfo{Frequency} = defined $Options{frequency} ? $Options{frequency} : undef;
  $OptionsInfo{FrequencyBins} = defined $Options{frequencybins} ? $Options{frequencybins} : undef;
  $OptionsInfo{IgnoreGaps} = defined $Options{ignoregaps} ? $Options{ignoregaps} : undef;
  $OptionsInfo{Longest} = defined $Options{longest} ? $Options{longest} : undef;
  $OptionsInfo{Shortest} = defined $Options{shortest} ? $Options{shortest} : undef;
  $OptionsInfo{SequenceLengths} = defined $Options{sequencelengths} ? $Options{sequencelengths} : undef;

  $OptionsInfo{FrequencyAnalysis} = ($Options{all} || $Options{frequency}) ? 1 : 0;
  $OptionsInfo{ListLongestSequence} = ($Options{all} || $Options{longest}) ? 1 : 0;
  $OptionsInfo{ListShortestSequence} = ($Options{all} || $Options{shortest}) ? 1 : 0;
  $OptionsInfo{ListSequenceLengths} = ($Options{all} || $Options{sequencelengths}) ? 1 : 0;
  $OptionsInfo{IgnoreGaps} = ($Options{ignoregaps} =~ /Yes/i) ? 1 : 0;

  # Setup frequency bin values...
  $OptionsInfo{NumOfBins} = 4;
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
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{detail} = 1;
  $Options{ignoregaps} = 'no';
  $Options{frequencybins} = 10;

  if (!GetOptions(\%Options, "all|a", "count|c", "detail|d=i", "frequency|f", "frequencybins=s", "help|h", "ignoregaps|i=s", "longest|l", "shortest|s", "sequencelengths", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if (!IsPositiveInteger($Options{detail})) {
    die "Error: The value specified, $Options{detail}, for option \"-d --detail\" is not valid. Allowed values: > 0\n";
  }
  if ($Options{ignoregaps} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{ignoregaps}, for option \"-i --IgnoreGaps\" is not valid. Allowed values: yes or no\n";
  }
}

__END__

=head1 NAME

InfoSequenceFiles.pl - List information about sequence and alignment files

=head1 SYNOPSIS

InfoSequenceFiles.pl SequenceFile(s) AlignmentFile(s)...

InfoSequenceFiles.pl [B<-a, --all>] [B<-c, --count>] [B<-d, --detail> infolevel]
[B<-f, --frequency>] [B<--FrequencyBins> number | "number, number, [number,...]"]
[B<-h, --help>] [B<-i, --IgnoreGaps> yes | no] [B<-l, --longest>] [B<-s, --shortest>]
[B<--SequenceLengths>] [B<-w, --workingdir> dirname] SequenceFile(s)...

=head1 DESCRIPTION

List information about contents of I<SequenceFile(s) and AlignmentFile(s)>: number of sequences,
shortest and longest sequences, distribution of sequence lengths and so on. The file names are
separated by spaces. All the sequence files in a current directory can be specified by I<*.aln>,
I<*.msf>, I<*.fasta>, I<*.fta>, I<*.pir> or any other supported formats; additionally, I<DirName>
corresponds to all the sequence files in the current directory with any of the supported file
extension: I<.aln, .msf, .fasta, .fta, and .pir>.

Supported sequence formats are: I<ALN/CLustalW>, I<GCG/MSF>, I<PILEUP/MSF>, I<Pearson/FASTA>,
and I<NBRF/PIR>. Instead of using file extensions, file formats are detected by parsing the contents
of I<SequenceFile(s) and AlignmentFile(s)>.

=head1 OPTIONS

=over 4

=item B<-a, --all>

List all the available information.

=item B<-c, --count>

List number of of sequences. This is B<default behavior>.

=item B<-d, --detail> I<InfoLevel>

Level of information to print about sequences during various options. Default: I<1>.
Possible values: I<1, 2 or 3>.

=item B<-f, --frequency>

List distribution of sequence lengths using the specified number of bins or bin range specified
using B<FrequencyBins> option.

This option is ignored for input files containing only single sequence.

=item B<--FrequencyBins> I<number | "number,number,[number,...]">

This value is used with B<-f, --frequency> option to list distribution of sequence lengths using
the specified number of bins or bin range. Default value: I<10>.

The bin range list is used to group sequence lengths  into different groups; It must contain
values in ascending order. Examples:

    100,200,300,400,500,600
    200,400,600,800,1000

The frequency value calculated for a specific bin corresponds to all the sequence lengths
which are greater than the previous bin value and less than or equal to the current bin value.

=item B<-h, --help>

Print this help message.

=item B<-i, --IgnoreGaps> I<yes | no>

Ignore gaps during calculation of sequence lengths. Possible values: I<yes or
no>. Default value: I<no>.

=item B<-l, --longest>

List information about longest sequence: ID, sequence and sequence length. This option
is ignored for input files containing only single sequence.

=item B<-s, --shortest>

List information about shortest sequence: ID, sequence and sequence length. This option
is ignored for input files containing only single sequence.

=item B<--SequenceLengths>

List information about sequence lengths.

=item B<-w, --WorkingDir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To count number of sequences in sequence files, type:

    % InfoSequenceFiles.pl Sample1.fasta
    % InfoSequenceFiles.pl Sample1.msf Sample1.aln Sample1.pir
    % InfoSequenceFiles.pl *.fasta *.fta *.msf *.pir *.aln

To list all available information with maximum level of available detail for a sequence
alignment file Sample1.msf, type:

    % InfoSequenceFiles.pl -a -d 3 Sample1.msf

To list sequence length information after ignoring sequence gaps in Sample1.aln file, type:

    % InfoSequenceFiles.pl --SequenceLengths --IgnoreGaps Yes
      Sample1.aln

To list shortest and longest sequence length information after ignoring sequence
gaps in Sample1.aln file, type:

    % InfoSequenceFiles.pl --longest --shortest --IgnoreGaps Yes
      Sample1.aln

To list distribution of sequence lengths after ignoring sequence gaps in Sample1.aln file and
report the frequency distribution into 10 bins, type:

    % InfoSequenceFiles.pl --frequency --FrequencyBins 10
      --IgnoreGaps Yes Sample1.aln

To list distribution of sequence lengths after ignoring sequence gaps in Sample1.aln file and
report the frequency distribution into specified bin range, type:

    % InfoSequenceFiles.pl --frequency --FrequencyBins
      "150,200,250,300,350" --IgnoreGaps Yes Sample1.aln

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AnalyzeSequenceFilesData.pl, ExtractFromSequenceFiles.pl, InfoAminoAcids.pl, InfoNucleicAcids.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
