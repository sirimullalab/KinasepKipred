#!/usr/bin/perl -w
#
# File: DownloadPDBFiles.pl
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
use File::Fetch;
use File::Copy;
use Text::ParseWords;
use Benchmark;
use FileUtil;
use TextUtil;
use PDBFileUtil;

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

# Process options...
print "Processing options...\n";
my(%OptionsInfo, %PDBIDsFileInfo);
ProcessOptions();

# Collect PDB IDs and download corresponding files...
my(@PDBIDs);
RetrievePDBIDs();
DownloadPDBFiles();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Download appropriate PDB fies...
sub DownloadPDBFiles {
  my($PDBFilesLocationURL, $PDBIDCount, $PDBIDOkayCount, $PDBIDFailedCount, $PDBIDsIgnoredCount, $PDBID, $PDBRemoteFileURL, $PDBLocalFile, $DownloadEDSMap, $EDSMapFilesLocaltionURL, $EDSMapRemoteFileURL, $EDSMapLocalFile, $Status, $EDSMapOkayCount, $EDSMapFailedCount, $EDSMapPDBID, %PDBIDsMap);

  print "\nDownloading PDB files...\n";

  $PDBFilesLocationURL = $OptionsInfo{DataLocationURL};

  $DownloadEDSMap = $OptionsInfo{EDSMap};
  $EDSMapFilesLocaltionURL = $OptionsInfo{EDSMapLocationURL};

  ($PDBIDCount, $PDBIDOkayCount, $PDBIDFailedCount, $PDBIDsIgnoredCount, $EDSMapOkayCount, $EDSMapFailedCount) = (0) x 6;

  # Turn off warnings from File::Fetch
  $File::Fetch::WARN = 0;

  %PDBIDsMap = ();

  PDBID: for $PDBID (@PDBIDs) {
    $PDBIDCount++;

    print "\nProcessing PDB ID $PDBID...\n";

    if (exists $PDBIDsMap{$PDBID}) {
      $PDBIDsIgnoredCount++;
      warn "Warning: Ignoring duplicate PDB ID $PDBID\n";
      next PDBID;
    }
    else {
      $PDBIDsMap{$PDBID} = $PDBID;
    }

    if ($PDBID =~ /\./) {
      $PDBIDsIgnoredCount++;
      warn "Warning: Ignoring invalid PDB ID $PDBID\n";
      next PDBID;
    }

    # Download PDB file...
    $PDBRemoteFileURL = "${PDBFilesLocationURL}${PDBID}.pdb";
    $PDBLocalFile = "${PDBID}.pdb";

    print "Downloading PDB file: $PDBRemoteFileURL\n";
    $Status = DownloadFile($PDBRemoteFileURL, $PDBLocalFile);

    if (!$Status) {
      $PDBIDFailedCount++;
      next PDBID;
    }
    $PDBIDOkayCount++;

    # Download EDS map file...
    if (!$DownloadEDSMap) {
      next PDBID;
    }

    $EDSMapPDBID = lc $PDBID;
    $EDSMapRemoteFileURL = "${EDSMapFilesLocaltionURL}${EDSMapPDBID}.ccp4";
    $EDSMapLocalFile = "${EDSMapPDBID}.ccp4";

    print "Downloading PDB EDS map file: $EDSMapRemoteFileURL\n";
    $Status = DownloadFile($EDSMapRemoteFileURL, $EDSMapLocalFile );

    if (!$Status) {
      $EDSMapFailedCount++;
      next PDBID;
    }
    print "Moving file from ${EDSMapPDBID}.ccp4 to ${PDBID}.ccp4\n";
    move "${EDSMapPDBID}.ccp4", "${EDSMapPDBID}Tmp.ccp4" or "Warning: Couldn't move file ${EDSMapPDBID}.ccp4 to ${EDSMapPDBID}Tmp.ccp4\n";
    move "${EDSMapPDBID}Tmp.ccp4", "${PDBID}.ccp4" or "Warning: Couldn't move file ${EDSMapPDBID}Tmp.ccp4 to ${PDBID}.ccp4\n";
    $EDSMapOkayCount++;

  }

  print "\nTotal number of PDB IDs:  $PDBIDCount\n";
  print "Number of duplicate PDB IDs ignored:  $PDBIDsIgnoredCount\n";

  print "\nNumber of successful downloads:  $PDBIDOkayCount\n";
  print "Number of failed downloads:  $PDBIDFailedCount\n";

  if ($DownloadEDSMap) {
    print "\nNumber of successful EDS map downloads:  $EDSMapOkayCount\n";
    print "Number of failed EDS map downloads:  $EDSMapFailedCount\n";
  }
}

# Download specified file...
sub DownloadFile {
  my($RemoteFileURL, $LocalFileName) = @_;
  my($Status, $FileFetch, $FetchedFilePath);

  $Status = 1;

  # Setup a fetch object...
  $FileFetch = File::Fetch->new(uri => $RemoteFileURL);

  # Fetch file to the CWD...
  $FetchedFilePath = $FileFetch->fetch();

  if (IsEmpty($FetchedFilePath)) {
    warn "Warning: Download failed for file $RemoteFileURL: " . $FileFetch->error() . "\n";
    if (-e $LocalFileName) {
      warn "Warning: Deleting empty file $LocalFileName\n";
      unlink $LocalFileName or warn "Warning: Couldn't delete file $LocalFileName\n";
    }
    $Status = 0;
  }
  return $Status;
}

# Collect PDB IDs...
#
sub RetrievePDBIDs {
  @PDBIDs = ();

  if ($OptionsInfo{Mode} =~ /^IDsOnCmdLine$/i) {
    RetriveCommandLinePDBIDs();
  }
  elsif ($OptionsInfo{Mode} =~ /^IDsInFile$/i) {
    RetriveTextFilePDBIDs();
  }
}

# Collect PDB IDs specified on the command line...
#
sub RetriveCommandLinePDBIDs {
  my($SpecifiedPDBID, @ProcessedPDBIDs );

  print "\nProcessing PDB ID(s) from command line...\n";

  for $SpecifiedPDBID (@{$OptionsInfo{CmdLinePDBIDs}}) {
    @ProcessedPDBIDs = ProcessSpecifiedPDBIDs($SpecifiedPDBID);
    if (@ProcessedPDBIDs) {
      push @PDBIDs, @ProcessedPDBIDs;
    }
  }
}

# Collect PDB IDs specified in the text file...
#
sub RetriveTextFilePDBIDs {
  my($TextFile, $InDelim, $IDsColIndex, $LineCount, $ProcessedLineCount, $IgnoredLineCount, $PDBID, $Line, @ProcessedPDBIDs, @LineWords);

  $TextFile = $PDBIDsFileInfo{Name};

  $IDsColIndex = $PDBIDsFileInfo{IDsColIndex};
  $InDelim = $PDBIDsFileInfo{InDelim} ;

  ($LineCount, $ProcessedLineCount, $IgnoredLineCount) = (0) x 3;

  print "\nProcessing PDB ID(s) from PDB IDs  file $TextFile...\n";

  open TEXTFILE, "$TextFile" or die "Couldn't open $TextFile: $! \n";
  # Skip label line...
  $_ = <TEXTFILE>;

  LINE: while ($Line = GetTextLine(\*TEXTFILE)) {
    $LineCount++;
    @LineWords = quotewords($InDelim, 0, $Line);

    if ($IDsColIndex >= scalar @LineWords) {
      $IgnoredLineCount++;
      warn "Warning: Ignoring line number $LineCount: PDB IDs column number, ". ($IDsColIndex + 1) . ", doesn't exist in the line containing, " . (scalar @LineWords) .  ", columns.\nLine: $Line\n";
      next LINE;
    }
    $PDBID = $LineWords[$IDsColIndex];
    if (IsEmpty($PDBID )) {
      $IgnoredLineCount++;
      warn "Warning: Ignoring line number $LineCount: PDB ID value is empty.\nLine: $Line\n";
      next LINE;
    }
    $ProcessedLineCount++;

    @ProcessedPDBIDs = ProcessSpecifiedPDBIDs($PDBID);
    if (@ProcessedPDBIDs) {
      push @PDBIDs, @ProcessedPDBIDs;
    }
  }
  close TEXTFILE;

  print "\nTotal number of lines in PDB IDs text file: $LineCount\n";
  print "Total number of lines processed: $ProcessedLineCount\n";
  print "Total number of lines ignored: $IgnoredLineCount\n";
}

# Process specified PDB IDs...
#
# Notes:
#   . Commas and spaces in the specification of PBD IDs are allowed.
#   . All PDB IDs are turned into uppercase letters.
#
sub ProcessSpecifiedPDBIDs {
  my($SpecifiedPDBID) = @_;
  my($PDBID, @PDBIDWords, @PDBIDs);

  $SpecifiedPDBID = RemoveLeadingAndTrailingWhiteSpaces($SpecifiedPDBID);
   if ($SpecifiedPDBID =~ / /) {
    @PDBIDWords = split " ",  $SpecifiedPDBID;
  }
  elsif ($SpecifiedPDBID =~ /,/) {
    @PDBIDWords = split ",",  $SpecifiedPDBID;
  }
  else {
    push @PDBIDWords, $SpecifiedPDBID;
  }

  @PDBIDs = ();
  for $PDBID (@PDBIDWords) {
    $PDBID =~ s/( |,)//g;
    push @PDBIDs, uc $PDBID;
  }
  return @PDBIDs;
}

# Process option values...
sub ProcessOptions {

  %OptionsInfo = ();
  %PDBIDsFileInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{DataLocationURL} = $Options{datalocationurl};
  if (IsEmpty($OptionsInfo{DataLocationURL} )) {
    die "Error: PDB data location URL specified using \"-d, --dataLocationURL\" is empty. Allowed value: Non empty string\n";
  }

  $OptionsInfo{EDSMap} = $Options{edsmap} =~ /^Yes$/i ? 1 : 0;
  $OptionsInfo{EDSMapLocationURL} = $Options{edsmaplocationurl};

  $OptionsInfo{ColMode} = $Options{colmode};
  $OptionsInfo{InDelim} = $Options{indelim};

  $OptionsInfo{PDBIDsCol } = defined $Options{pdbidscol} ? $Options{pdbidscol} : '';

  @{$OptionsInfo{CmdLinePDBIDs}} = ();
  $OptionsInfo{PDBIDsFile} = "";

  if ($OptionsInfo{Mode} =~ /^IDsOnCmdLine$/i) {
    push @{$OptionsInfo{CmdLinePDBIDs}}, @ARGV;
  }
  elsif ($OptionsInfo{Mode} =~ /^IDsInFile$/i) {
    if (@ARGV != 1) {
      die "Error: Invalid number of PDB IDs text files, ". (scalar @ARGV) . ",specified on the command line for \"IDsInFile\" value of  $Options{mode}, for option \"-m --mode\". Allowed value: Only one text file\n";
    }
    $OptionsInfo{PDBIDsFile} = $ARGV[0];

    RetrievePDBIDsTextFileInfo();
  }
  else {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: IDsOnCmdLine or IDsInFile\n";
  }
}

# Retrieve information for PDB IDs text file...
#
sub RetrievePDBIDsTextFileInfo {
  my($TextFile, $FileDir, $FileName, $FileExt, $InDelim, $Line, $ColNum, $ColLabel, $PDBIDsColIndex, $ColMode, $PDBIDsCol , $ColCount, $ColIndex, @ColLabels);

  $TextFile = $OptionsInfo{PDBIDsFile};

  %PDBIDsFileInfo = ();
  $PDBIDsFileInfo{Name} = $TextFile;

  $PDBIDsFileInfo{ColCount} = 0;
  @{$PDBIDsFileInfo{ColLabels}} = ();
  %{$PDBIDsFileInfo{ColLabelToNumMap}} = ();
  $PDBIDsFileInfo{InDelim} = "";

  $PDBIDsFileInfo{IDsColIndex} = "";

  if (!-e $TextFile) {
    die "Error: PDBIDs text file, $TextFile, doesn't exist\n";
  }

  if (!CheckFileType($TextFile, "csv tsv")) {
    die "Error: Ignoring file $TextFile: It's not a csv or tsv file\n";
  }

  ($FileDir, $FileName, $FileExt) = ParseFileName($TextFile);
  if ($FileExt =~ /^tsv$/i) {
    $InDelim = "\t";
  }
  else {
    $InDelim = "\,";
    if (!($OptionsInfo{InDelim} =~ /^(comma|semicolon)$/i)) {
      die "Error: Ignoring file $TextFile: The value specified, $OptionsInfo{InDelim}, for option \"--indelim\" is not valid for textfile\n";
    }
    if ($OptionsInfo{InDelim} =~ /^semicolon$/i) {
      $InDelim = "\;";
    }
  }

  if (!open TEXTFILE, "$TextFile") {
    die "Error: Ignoring file $TextFile: Couldn't open it: $! \n";
  }

  $Line = GetTextLine(\*TEXTFILE);
  @ColLabels = quotewords($InDelim, 0, $Line);
  close TEXTFILE;
  $ColCount = scalar @ColLabels;

  push @{$PDBIDsFileInfo{ColLabels}}, @ColLabels;
  $PDBIDsFileInfo{ColCount} = $ColCount ;
  $PDBIDsFileInfo{InDelim} = $InDelim;

  # Setup collabel to colnum map...
  %{$PDBIDsFileInfo{ColLabelToNumMap}} = ();
  for $ColNum (0 .. $#ColLabels) {
    $ColLabel = $ColLabels[$ColNum];
    $PDBIDsFileInfo{ColLabelToNumMap}{lc $ColLabel} = $ColNum;
  }

  # Identify column containing PDB IDs...
  $PDBIDsColIndex = "";

  $ColMode = $OptionsInfo{ColMode};
  $PDBIDsCol = $OptionsInfo{PDBIDsCol };

  if (IsNotEmpty($PDBIDsCol )) {
    if ($ColMode =~ /^collabel$/i) {
      $ColLabel = lc $PDBIDsCol;
      if (!exists $PDBIDsFileInfo{ColLabelToNumMap}{$ColLabel} ) {
	die "Error: Ignoring file $TextFile: The column name, $PDBIDsCol, specified for option \"-p, --PDBIDsCol \" is not valid for text file\n";
      }
      $PDBIDsColIndex = $PDBIDsFileInfo{ColLabelToNumMap}{$ColLabel};
    }
    else {
      $ColNum = $PDBIDsCol;
      if ($ColNum <= 0 || $ColNum > $ColCount) {
	die "Error: Ignoring file $TextFile: The column number, $PDBIDsCol, specified for option \"-p, --PDBIDsCol \" is not valid for text file. It must be > 0 and <= $ColCount\n";
      }
      $PDBIDsColIndex = $ColNum - 1;
    }
  }
  else {
    # Look for column name containing PDB_ID or PDBID text string...
    $PDBIDsCol = "";
    $ColIndex = 0;
    COLLABEL: for $ColLabel (@ColLabels) {
      if ($ColLabel =~ /(PDB_ID|PDBID)/i) {
	$PDBIDsCol = $ColLabel;
	$PDBIDsColIndex = $ColIndex;
	last COLLABEL;
      }
      $ColIndex++;
    }
    if (IsEmpty($PDBIDsCol)) {
      die "Error: Ignoring file $TextFile: Couldn't find PDB IDs default column containing text string PDB_ID or PDBID in its name\n";
    }
  }
  $PDBIDsFileInfo{IDsColIndex} = $PDBIDsColIndex;
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{colmode} = "colnum";
  $Options{datalocationurl} = "http://www.rcsb.org/pdb/files/";
  $Options{edsmap} = "no";
  $Options{edsmaplocationurl} = "http://www.ebi.ac.uk/pdbe/coordinates/files/";
  $Options{indelim} = "comma";
  $Options{mode} = "IDsOnCmdLine";

  if (!GetOptions(\%Options, "colmode|c=s", "datalocationurl|d=s", "edsmap|e=s", "edsmaplocationurl=s", "help|h",  "indelim=s", "mode|m=s", "pdbidscol|s=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{colmode} !~ /^(colnum|collabel)$/i) {
    die "Error: The value specified, $Options{colmode}, for option \"-c, --colmode\" is not valid. Allowed values: colnum or collabel\n";
  }
  if ($Options{edsmap} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{edsmap}, for option \"-e --edsmap\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{indelim} !~ /^(comma|semicolon)$/i) {
    die "Error: The value specified, $Options{indelim}, for option \"--indelim\" is not valid. Allowed values: comma or semicolon\n";
  }
  if ($Options{mode} !~ /^(IDsOnCmdLine|IDsInFile)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: IDsOnCmdLine or IDsInFile\n";
  }
}

__END__

=head1 NAME

DownloadPDBFiles.pl - Download PDB files for PDB ID(s)

=head1 SYNOPSIS

DownloadPDBFiles.pl PDBID(s) or PDBIDsTextFile

DownloadPDBFiles.pl [B<-c, --colmode> I<colnum | collabel>]
[B<-d, --dataLocationURL> I<PDB URL>] [B<-h, --help>] [B<-d, --dataLocationURL> I<PDB URL>]
[B<-e, --EDSMap > I<yes | no>] [B<--EDSMapLocationURL> I<EDS Map URL>]
[B<-m, --mode> <IDsOnCmdLine | IDsInFile>] [B<--PDBIDsCol > I<number | string>]
[B<-w, --WorkingDir> dirname] PDBID(s) or PDBIDsTextFile

=head1 DESCRIPTION

Download PDB files corresponding to PDB IDs specified in a column in a CSV/TSV text file or
on the command line as space delimited parameters.

=head1 OPTIONS

=over 4

=item B<-c, --colmode> I<colnum | collabel>

Specify how columns are identified in a I<TextFile> containing PDB IDs: using column number
or column label. Possible values: I<colnum or collabel>. Default value: I<colnum>.

=item B<-d, --dataLocationURL> I<PDB URL>

Specify location of PDB URL where data files are available for download. Default value:
I<http://www.rcsb.org/pdb/files/>.

=item B<-e, --EDSMap > I<yes | no>

Download Electron Density Map (EDS) in CCP4 format. Possible values: I<Yes or No>.
Default value: I<no>.

=item B<--EDSMapLocationURL> I<EDS Map URL>

Specify location of EDS Map URL where data files are available for download. Default value:
I<http://www.ebi.ac.uk/pdbe/coordinates/files/>.

=item B<-h, --help>

Print this help message.

=item B<--indelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile> containing PDB IDs. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a delimiter.

=item B<-m, --mode> <IDsOnCmdLine | IDsInFile>

Indicate how PDB IDs are specified: PDB IDs are either present as space delimited command line
parameters or in a specific column in a CSV/TSV text file. Possible values: I<IDsOnCmdLine or  IDsInFile>.
Default: I<IDsOnCmdLine>.

=item B<-p, --PDBIDsCol > I<number | string>

Column used to identify PDB ID(s) in a text file. Default value: First column containing text
string B<PDB_ID> or <PDBID>.

For I<colnum> value of B<-c, --colmode> option, input value is a column number.
Example: I<1>.

For I<collabel> value of B<-c, --colmode> option, input value is a column label.
Example: I<PDB_ID>.

This option is ignored during I<IDsOnCmdLine> value of B<m, --mode> option.

=item B<-w, --WorkingDir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To retrieve a PDB files for PDB ID 2HYY and generate a local 2HYY.pdb file, type:

    % DownloadPDBFiles.pl 2HYY

To retrieve PDB files for multiple PDB IDs 2HYY and 1KV2 and generate corresponding
local PDB files, type:

    % DownloadPDBFiles.pl 2HYY 1KV2

To download PDB files for PDB IDs present in column name PDB_ID or PDBID in
SamplePDBIDs.csv file and generate correponding PDB files, type

    % DownloadPDBFiles.pl -m IDsInFile SamplePDBIDs.csv

To download PDB files for PDB IDs present in a specific column name in
SamplePDBIDs.csv file and generate correponding PDB files, type

    % DownloadPDBFiles.pl -m IDsInFile -c collabel -p PDB_ID SamplePDBIDs.csv

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

ExtractFromPDBFiles.pl,  InfoPDBFiles.pl, ModifyPDBFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
