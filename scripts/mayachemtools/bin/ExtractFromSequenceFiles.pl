#!/usr/bin/perl -w
#
# File: ExtractFromSequenceFiles.pl
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

my($ScriptName, %Options, $StartTime, $EndTime, $TotalTime);

# Autoflush STDOUT
$| = 1;

# Starting message...
$ScriptName = basename($0);
print "\n$ScriptName: Starting...\n\n";
$StartTime = new Benchmark;

# Setup script usage message...
SetupScriptUsage();
if ($Options{help} || @ARGV < 1) {
  die GetUsageFromPod("$FindBin::Bin/$ScriptName");
}

# Expand wild card file names...
my(@SequenceFilesList);
@SequenceFilesList = ExpandFileNames(\@ARGV, "aln msf fasta fta pir");

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Set up information about input files...
print "Checking input sequence file(s)...\n";
my(%SequenceFilesInfo);
RetrieveSequenceFilesInfo();

# Process input files..
my($FileIndex);
if (@SequenceFilesList > 1) {
  print "\nProcessing sequence files...\n";
}
for $FileIndex (0 .. $#SequenceFilesList) {
  if ($SequenceFilesInfo{FilesOkay}[$FileIndex]) {
    print "\nProcessing file $SequenceFilesList[$FileIndex]...\n";
    ExtractFromSequenceFiles($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Extract from sequence files...
sub ExtractFromSequenceFiles {
  my($FileIndex) = @_;
  my($OutSequenceFile, $SequenceFile, $SequenceDataRef, $SpecifiedSequenceDataRef);

  # Read sequence file...
  $SequenceFile = $SequenceFilesList[$FileIndex];
  open SEQUENCEFILE, "$SequenceFile" or die "Error: Can't open $SequenceFile: $! \n";
  $SequenceDataRef = ReadSequenceFile($SequenceFile);
  close SEQUENCEFILE;

  $OutSequenceFile = $SequenceFilesInfo{OutFile}[$FileIndex];
  print "Generating sequence file $OutSequenceFile...\n";

  # Retrieve sequence data for specified sequences...
  $SpecifiedSequenceDataRef = GetSpecifiedSequenceData($SequenceDataRef);

  # Handle gaps...
  if ($OptionsInfo{IgnoreGaps}) {
    if (@{$SpecifiedSequenceDataRef->{IDs}} > 1) {
      if (AreSequenceLengthsIdentical($SpecifiedSequenceDataRef)) {
	$SpecifiedSequenceDataRef = RemoveSequenceAlignmentGapColumns($SpecifiedSequenceDataRef);
      }
    }
    else {
      # Remove the gaps from the sequence...
      my($ID, $Sequence);
      $ID = $SpecifiedSequenceDataRef->{IDs}[0];
      $Sequence = $SpecifiedSequenceDataRef->{Sequence}{$ID};
      $SpecifiedSequenceDataRef->{Sequence}{$ID} = RemoveSequenceGaps($Sequence);
    }
  }

  # Write out the file...
  WritePearsonFastaSequenceFile($OutSequenceFile, $SpecifiedSequenceDataRef, $OptionsInfo{MaxSequenceLength});
}

# Get specified sequence data...
sub GetSpecifiedSequenceData {
  my($SequenceDataRef) = @_;

  if ($OptionsInfo{Mode} =~ /^SequenceID$/i) {
    return GetDataBySequenceIDs($SequenceDataRef);
  }
  elsif ($Options{mode} =~ /^SequenceNum$/i) {
    return GetDataBySequenceNums($SequenceDataRef);
  }
  elsif ($Options{mode} =~ /^SequenceNumRange$/i) {
    return GetDataBySequenceNumRange($SequenceDataRef);
  }
  else {
    return undef;
  }
}

# Get specified sequence data...
sub GetDataBySequenceIDs {
  my($SequenceDataRef) = @_;
  my($ID, $SequenceCount, $IDMatched, $SpecifiedID, %SpecifiedSequenceDataMap);

  # Go over sequences and collect sequences for writing out a new sequence file...
  %SpecifiedSequenceDataMap = ();
  @{$SpecifiedSequenceDataMap{IDs}} = ();
  %{$SpecifiedSequenceDataMap{Description}} = ();
  %{$SpecifiedSequenceDataMap{Sequence}} = ();

  $SequenceCount = 0;
  ID: for $ID (@{$SequenceDataRef->{IDs}}) {
    if ($OptionsInfo{MatchExactSequenceIDs}) {
      if (!exists $OptionsInfo{SpecifiedSequenceIDsMap}{lc($ID)}) {
	next ID;
      }
      if ($SequenceCount >= scalar @{$OptionsInfo{SpecifiedSequenceIDs}}) {
	last ID;
      }
      $SequenceCount++;
    }
    else {
      # Does this ID contains specified ID as substring...
      $IDMatched = 0;
      SPECIFIEDID: for $SpecifiedID (@{$OptionsInfo{SpecifiedSequenceIDs}}) {
	if ($ID =~ /$SpecifiedID/i) {
	  $IDMatched = 1;
	  last SPECIFIEDID;
	}
      }
      if (!$IDMatched) {
	next ID;
      }
      $SequenceCount++;
    }
    # Collect sequence data...
    push @{$SpecifiedSequenceDataMap{IDs}}, $ID;
    $SpecifiedSequenceDataMap{Description}{$ID} = $SequenceDataRef->{Description}{$ID};
    $SpecifiedSequenceDataMap{Sequence}{$ID} = $SequenceDataRef->{Sequence}{$ID};
  }

  return \%SpecifiedSequenceDataMap;
}

# Get specified sequence data...
sub GetDataBySequenceNums {
  my($SequenceDataRef) = @_;
  my($ID, $SequenceNum, $SequenceCount, %SpecifiedSequenceDataMap);

  # Go over sequences and collect sequences for writing out a new sequence file...
  %SpecifiedSequenceDataMap = ();
  @{$SpecifiedSequenceDataMap{IDs}} = ();
  %{$SpecifiedSequenceDataMap{Description}} = ();
  %{$SpecifiedSequenceDataMap{Sequence}} = ();

  $SequenceNum = 0;
  $SequenceCount = 0;
  ID: for $ID (@{$SequenceDataRef->{IDs}}) {
    $SequenceNum++;
    if (!exists $OptionsInfo{SpecifiedSequenceIDsMap}{$SequenceNum}) {
      next ID;
    }
    if ($SequenceCount >= scalar @{$OptionsInfo{SpecifiedSequenceIDs}}) {
      last ID;
    }
    $SequenceCount++;

    # Collect sequence data...
    push @{$SpecifiedSequenceDataMap{IDs}}, $ID;
    $SpecifiedSequenceDataMap{Description}{$ID} = $SequenceDataRef->{Description}{$ID};
    $SpecifiedSequenceDataMap{Sequence}{$ID} = $SequenceDataRef->{Sequence}{$ID};
  }

  return \%SpecifiedSequenceDataMap;
}

# Get specified sequence data...
sub GetDataBySequenceNumRange {
  my($SequenceDataRef) = @_;
  my($ID, $SequenceNum, $SequenceCount, %SpecifiedSequenceDataMap);

  # Go over sequences and collect sequences for writing out a new sequence file...
  %SpecifiedSequenceDataMap = ();
  @{$SpecifiedSequenceDataMap{IDs}} = ();
  %{$SpecifiedSequenceDataMap{Description}} = ();
  %{$SpecifiedSequenceDataMap{Sequence}} = ();

  $SequenceNum = 0;
  $SequenceCount = 0;
  ID: for $ID (@{$SequenceDataRef->{IDs}}) {
    $SequenceNum++;

    if (!($SequenceNum >= $OptionsInfo{SpecifiedSequenceIDs}[0] && $SequenceNum <= $OptionsInfo{SpecifiedSequenceIDs}[1])) {
      next ID;
    }
    if ($SequenceNum > $OptionsInfo{SpecifiedSequenceIDs}[1]) {
      last ID;
    }
    $SequenceCount++;
    # Collect sequence data...
    push @{$SpecifiedSequenceDataMap{IDs}}, $ID;
    $SpecifiedSequenceDataMap{Description}{$ID} = $SequenceDataRef->{Description}{$ID};
    $SpecifiedSequenceDataMap{Sequence}{$ID} = $SequenceDataRef->{Sequence}{$ID};
  }

  return \%SpecifiedSequenceDataMap;
}


# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  # Miscellaneous options...
  $OptionsInfo{IgnoreGaps} = ($Options{ignoregaps} =~ /Yes/i) ? 1 : 0;

  $OptionsInfo{Mode} = $Options{mode};
  $OptionsInfo{MatchExactSequenceIDs} = $Options{sequenceidmatch} =~ /Exact/i ? 1 :0;

  # Check specified sequences value...
  $OptionsInfo{SpecifiedSequences} = $Options{sequences};
  @{$OptionsInfo{SpecifiedSequenceIDs}} = ();
  %{$OptionsInfo{SpecifiedSequenceIDsMap}} = ();

  my(@SpecifiedSequenceIDs) = ();
  if ($Options{mode} =~ /^SequenceID$/i) {
    if (!$Options{sequences}) {
      die "Error: No value specified for option \"-s, --Sequences\" during \"SequenceID\" of \"-m, --mode\" option\n";
    }
    @SpecifiedSequenceIDs = split /\,/, $Options{sequences};
  }
  elsif ($Options{mode} =~ /^SequenceNum$/i) {
    if ($Options{sequences}) {
      @SpecifiedSequenceIDs = split /\,/, $Options{sequences};
      my($SequenceNum);
      for $SequenceNum (@SpecifiedSequenceIDs) {
	if (!IsPositiveInteger($SequenceNum)) {
	  die "Error: The value specified, $SequenceNum, in \"$Options{sequences}\" for option \"-s, --Sequences\" is not valid: Valid values: > 0\n";
	}
      }
    }
    else {
      push @SpecifiedSequenceIDs, "1";
    }
  }
  elsif ($Options{mode} =~ /^SequenceNumRange$/i) {
    if (!$Options{sequences}) {
      die "Error: No value specified for option \"-s, --Sequences\" during \"SequenceNumRange\" of \"-m, --mode\" option\n";
    }
    @SpecifiedSequenceIDs = split /\,/, $Options{sequences};
    if (@SpecifiedSequenceIDs != 2) {
      die "Error: The number of values", scalar @SpecifiedSequenceIDs, " specified, $Options{sequences}, for option \"-s, --Sequences\" are not valid. Number of values must be 2 to indicate starting and ending sequence number.\n";
    }
    my($SequenceNum);
    for $SequenceNum (@SpecifiedSequenceIDs) {
      if (!IsPositiveInteger($SequenceNum)) {
	die "Error: The value specified, $SequenceNum, in \"$Options{sequences}\" for option \"-s, --Sequences\" is not valid: Valid values: > 0\n";
      }
    }
    if ($SpecifiedSequenceIDs[0] > $SpecifiedSequenceIDs[1]) {
      die "Error: The value specified \"$Options{sequences}\" for option \"-s, --Sequences\" are not valid: Starting sequence number $SpecifiedSequenceIDs[0] must be smaller than ending sequence number $SpecifiedSequenceIDs[1]\n";
    }
  }
  push @{$OptionsInfo{SpecifiedSequenceIDs}}, @SpecifiedSequenceIDs;
  my($SequenceID);
  for $SequenceID (@SpecifiedSequenceIDs) {
    if ($Options{mode} =~ /^SequenceID$/i) {
      $OptionsInfo{SpecifiedSequenceIDsMap}{lc($SequenceID)} = $SequenceID;
    }
    else {
      $OptionsInfo{SpecifiedSequenceIDsMap}{$SequenceID} = $SequenceID;
    }
  }

  $OptionsInfo{MaxSequenceLength} = $Options{sequencelength};
  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? 1 : 0;
  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : 0;
}

# Retrieve information about sequence files...
sub RetrieveSequenceFilesInfo {
  my($Index, $SequenceFile, $FileSupported, $FileFormat, $SequenceCount, $FileDir, $FileName, $FileExt, $OutFileRoot, $OutFileExt, $OutFileMode, $SequenceDataRef);

  %SequenceFilesInfo = ();
  @{$SequenceFilesInfo{FilesOkay}} = ();
  @{$SequenceFilesInfo{OutFileRoot}} = ();
  @{$SequenceFilesInfo{OutFileExt}} = ();
  @{$SequenceFilesInfo{OutFile}} = ();
  @{$SequenceFilesInfo{Format}} = ();
  @{$SequenceFilesInfo{SequenceCount}} = ();

  FILELIST: for $Index (0 .. $#SequenceFilesList) {
    $SequenceFile = $SequenceFilesList[$Index];
    $SequenceFilesInfo{FilesOkay}[$Index] = 0;
    $SequenceFilesInfo{OutFileRoot}[$Index] = '';
    $SequenceFilesInfo{OutFileExt}[$Index] = '';
    $SequenceFilesInfo{OutFile}[$Index] = '';
    $SequenceFilesInfo{Format}[$Index] = 'NotSupported';
    $SequenceFilesInfo{SequenceCount}[$Index] = 0;

    if (! open SEQUENCEFILE, "$SequenceFile") {
      warn "Warning: Ignoring file $SequenceFile: Couldn't open it: $! \n";
      next FILELIST;
    }
    close SEQUENCEFILE;

    ($FileSupported, $FileFormat) = IsSupportedSequenceFile($SequenceFile);
    if (!$FileSupported) {
      warn "Warning: Ignoring file $SequenceFile: Sequence file format is not supported.\n";
      next FILELIST;
    }
    $SequenceDataRef = ReadSequenceFile($SequenceFile);

    $SequenceCount = $SequenceDataRef->{Count};
    if (!$SequenceCount) {
      warn "Warning: Ignoring file $SequenceFile: Sequence data is missing.\n";
      next FILELIST;
    }

    # Setup output file names...
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($SequenceFile);
    $OutFileExt = 'fasta';
    if ($OptionsInfo{OutFileRoot} && (@SequenceFilesList == 1)) {
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
    MODE: {
	if ($OptionsInfo{Mode} =~ /^SequenceID$/i) { $OutFileMode = 'SequenceID'; last MODE;}
	if ($OptionsInfo{Mode} =~ /^SequenceNum$/i) { $OutFileMode = 'SequenceNum'; last MODE;}
	if ($OptionsInfo{Mode} =~ /^SequenceNumRange$/i) { $OutFileMode = 'SequenceNumRange'; last MODE;}
	$OutFileMode = '';
    }
    if (!$OptionsInfo{OverwriteFiles}) {
      if (-e "${OutFileRoot}${OutFileMode}.${OutFileExt}") {
	warn "Warning: Ignoring file $SequenceFile: The file ${OutFileRoot}${OutFileMode}.${OutFileExt} already exists\n";
	next FILELIST;
      }
    }

    $SequenceFilesInfo{FilesOkay}[$Index] = 1;
    $SequenceFilesInfo{OutFileRoot}[$Index] = $OutFileRoot;
    $SequenceFilesInfo{OutFileExt}[$Index] = $OutFileExt;
    $SequenceFilesInfo{OutFile}[$Index] = "${OutFileRoot}${OutFileMode}.${OutFileExt}";

    $SequenceFilesInfo{Format}[$Index] = $FileFormat;
    $SequenceFilesInfo{SequenceCount}[$Index] = $SequenceCount;
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{ignoregaps} = 'Yes';
  $Options{mode} = 'SequenceNum';
  $Options{sequenceidmatch} = 'Relaxed';
  $Options{sequencelength} = 80;

  if (!GetOptions(\%Options, "help|h", "ignoregaps|i=s", "mode|m=s", "overwrite|o", "root|r=s", "sequences|s=s", "sequenceidmatch=s", "sequencelength=i", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{ignoregaps} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{ignoregaps}, for option \"-q --quote\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{mode} !~ /^(SequenceID|SequenceNum|SequenceNumRange)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: SequenceID, SequenceNum, or SequenceNumRange\n";
  }
  if ($Options{sequenceidmatch} !~ /^(Exact|Relaxed)$/i) {
    die "Error: The value specified, $Options{sequenceidmatch}, for option \"--SequenceIDMatch\" is not valid. Allowed values: Exact or Relaxed\n";
  }
  if (!IsPositiveInteger($Options{sequencelength})) {
    die "Error: The value specified, $Options{sequencelength}, for option \"--SequenceLength\" is not valid. Allowed values: >0\n";
  }
}

__END__

=head1 NAME

ExtractFromSequenceFiles.pl - Extract data from sequence and alignment files

=head1 SYNOPSIS

ExtractFromSequenceFiles.pl SequenceFile(s) AlignmentFile(s)...

ExtractFromSequenceFiles.pl [B<-h, --help>] [B<-i, --IgnoreGaps> yes | no]
[B<-m, --mode> SequenceID | SequenceNum | SequenceNumRange] [B<-o, --overwrite>]
[B<-r, --root> rootname] [B<-s, --Sequences> "SequenceID, [SequenceID,...]" | "SequenceNum, [SequenceNum,...]" |
"StartingSeqNum, EndingSeqNum"] [B<--SequenceIDMatch> Exact | Relaxed]
[B<-w, --WorkingDir> dirname] SequenceFile(s) AlignmentFile(s)...

=head1 DESCRIPTION

Extract specific data from I<SequenceFile(s) and AlignmentFile(s)> and generate
FASTA files. You can extract sequences using sequence IDs or sequence numbers.

The file names are separated by spaces. All the sequence files in a current directory can
be specified by I<*.aln>, I<*.msf>, I<*.fasta>, I<*.fta>, I<*.pir> or any other supported
formats; additionally, I<DirName> corresponds to all the sequence files in the current directory
with any of the supported file extension: I<.aln, .msf, .fasta, .fta, and .pir>.

Supported sequence formats are: I<ALN/CLustalW>, I<GCG/MSF>, I<PILEUP/MSF>, I<Pearson/FASTA>,
and I<NBRF/PIR>. Instead of using file extensions, file formats are detected by parsing the contents
of I<SequenceFile(s) and AlignmentFile(s)>.

=head1 OPTIONS

=over 4

=item B<-h, --help>

Print this help message.

=item B<-i, --IgnoreGaps> I<yes | no>

Ignore gaps or gap columns during during generation of new sequence or alignment file(s).
Possible values: I<yes or no>. Default value: I<yes>.

In order to remove gap columns, length of all the sequence must be same; otherwise,
this option is ignored.

=item B<-m, --mode> I<SequenceID | SequenceNum | SequenceNumRange>

Specify how to extract data from sequence files: extract sequences using sequence
IDs or sequence numbers. Possible values: I<SequenceID | SequenceNum
| SequenceNumRange>. Default: I<SequenceNum> with value of 1.

The sequence numbers correspond to position of sequences starting from 1 for first sequence
in I<SequenceFile(s) and AlignmentFile(s)>.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

New sequence file name is generated using the root: <Root><Mode>.<Ext>. Default new file:
<SequenceFileName><Mode>.<Ext>. This option is ignored for multiple input files.

=item B<-s, --Sequences> I<"SequenceID,[SequenceID,...]" | "SequenceNum,[SequenceNum,...]" | "StartingSeqNum,EndingSeqNum">

This value is B<-m, --mode> specific. In general, it's a comma delimites list of sequence IDs or sequence
numbers.

For I<SequenceID> value of B<-m, --mode> option, input value format is: I<SequenceID,...>. Examples:

    ACHE_BOVIN
    ACHE_BOVIN,ACHE_HUMAN

For I<SequenceNum> value of B<-m, --mode> option, input value format is: I<SequenceNum,...>. Examples:

    2
    1,5

For I<SequenceNum> value of B<-m, --mode> option, input value format is: I<StaringSeqNum,EndingSeqNum>. Examples:

    2,4

=item B<--SequenceIDMatch> I<Exact | Relaxed>

Sequence IDs matching criterion during I<SequenceID> value of B<-m, --mode> option: match
specified sequence ID exactly or as sub string against sequence IDs in the files. Possible
values: I<Exact | Relaxed>. Default: I<Relaxed>. Sequence ID match is case insenstitive
during both options.

=item B<--SequenceLength> I<number>

Maximum sequence length per line in sequence file(s). Default: I<80>.

=item B<-w --WorkingDir> I<text>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To extract first sequence from Sample1.fasta sequence file and generate Sample1SequenceNum.fasta
sequence file, type:

    % ExtractFromSequenceFiles.pl -o Sample1.fasta

To extract first sequence from Sample1.aln alignment file and generate Sample1SequenceNum.fasta
sequence file without any column gaps, type:

    % ExtractFromSequenceFiles.pl -o Sample1.aln

To extract first sequence from Sample1.aln alignment file and generate Sample1SequenceNum.fasta
sequence file with column gaps, type:

    % ExtractFromSequenceFiles.pl --IgnroreGaps No -o Sample1.aln

To extract sequence number 1 and 4 from Sample1.fasta sequence file and generate
Sample1SequenceNum.fasta sequence file, type:

    % ExtractFromSequenceFiles.pl -o -m SequenceNum --Sequences 1,4
      -o Sample1.fasta

To extract sequences from sequence  number 1 to 4 from Sample1.fasta sequence file and generate
Sample1SequenceNumRange.fasta sequence file, type:

    % ExtractFromSequenceFiles.pl -o -m SequenceNumRange --Sequences
      1,4 -o Sample1.fasta

To extract sequence ID "Q9P993/104-387" from sequence  from Sample1.fasta sequence file and generate
Sample1SequenceID.fasta sequence file, type:

    % ExtractFromSequenceFiles.pl -o -m SequenceID --Sequences
      "Q9P993/104-387" --SequenceIDMatch Exact -o Sample1.fasta

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AnalyzeSequenceFilesData.pl, InfoSequenceFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
