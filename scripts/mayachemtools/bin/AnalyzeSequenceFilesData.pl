#!/usr/bin/perl -w
#
# File: AnalyzeSequenceFilesData.pl
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
use AminoAcids;
use NucleicAcids;

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

print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Set up information about input files...
print "Checking input sequence file(s)...\n";
my(%SequenceFilesInfo);
RetrieveSequenceFilesInfo();
SetupSequenceRegionsData();

# Process input files..
my($FileIndex);
if (@SequenceFilesList > 1) {
  print "\nProcessing sequence files...\n";
}
for $FileIndex (0 .. $#SequenceFilesList) {
  if ($SequenceFilesInfo{FilesOkay}[$FileIndex]) {
    print "\nProcessing file $SequenceFilesList[$FileIndex]...\n";
    AnalyzeSequenceFileData($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Analyze sequence file...
sub AnalyzeSequenceFileData {
  my($FileIndex) = @_;
  my($SequenceFile, $SequenceDataRef);

  $SequenceFile = $SequenceFilesList[$FileIndex];

  open SEQUENCEFILE, "$SequenceFile" or die "Error: Can't open $SequenceFile: $! \n";
  $SequenceDataRef = ReadSequenceFile($SequenceFile);
  close SEQUENCEFILE;

  if ($OptionsInfo{CalculatePercentIdentityMatrix}) {
    CalculatePercentIdentityMatrix($FileIndex, $SequenceDataRef);
  }
  if ($OptionsInfo{PerformResidueFrequencyAnalysis}) {
    PerformResidueFrequencyAnalysis($FileIndex, $SequenceDataRef);
  }
}

# Calculate percent identity matrix...
sub CalculatePercentIdentityMatrix {
  my($FileIndex, $SequenceDataRef) = @_;
  my($PercentIdentity, $PercentIdentityMatrixFile, $PercentIdentityMatrixRef, $RowID, $ColID, $Line, @LineWords);

  $PercentIdentityMatrixFile = $SequenceFilesInfo{OutFileRoot}[$FileIndex] . 'PercentIdentityMatrix.' . $SequenceFilesInfo{OutFileExt}[$FileIndex];
  $PercentIdentityMatrixRef = CalculatePercentSequenceIdentityMatrix($SequenceDataRef, $OptionsInfo{IgnoreGaps}, $OptionsInfo{Precision});

  print "Generating percent identity matrix file $PercentIdentityMatrixFile...\n";
  open OUTFILE, ">$PercentIdentityMatrixFile" or die "Can't open $PercentIdentityMatrixFile: $!\n";

  # Write out column labels...
  @LineWords = ();
  push @LineWords, '';
  for $ColID (@{$PercentIdentityMatrixRef->{IDs}}) {
    push @LineWords, $ColID;
  }
  $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print OUTFILE "$Line\n";

  # Write out rows...
  for $RowID (@{$PercentIdentityMatrixRef->{IDs}}) {
    @LineWords = ();
    push @LineWords, $RowID;
    for $ColID (@{$PercentIdentityMatrixRef->{IDs}}) {
      $PercentIdentity = $PercentIdentityMatrixRef->{PercentIdentity}{$RowID}{$ColID};
      push @LineWords, $PercentIdentity;
    }
    $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print OUTFILE "$Line\n";
  }
  close OUTFILE;
}

# Perform frequency analysis...
sub PerformResidueFrequencyAnalysis {
  my($FileIndex, $SequenceDataRef) = @_;

  CountResiduesInRegions($FileIndex, $SequenceDataRef);
  CalculatePercentResidueFrequencyInRegions($FileIndex, $SequenceDataRef);
  GeneratePercentResidueFrequencyOutFilesForRegions($FileIndex, $SequenceDataRef);
}

# Count residues...
sub CountResiduesInRegions {
  my($FileIndex, $SequenceDataRef) = @_;

  # Setup rerfernce sequence data...
  my($RefereceSequenceID, $RefereceSequence);
  $RefereceSequenceID = $SequenceFilesInfo{RefereceSequenceID}[$FileIndex];
  $RefereceSequence = $SequenceFilesInfo{RefereceSequence}[$FileIndex];

  # Count residues...
  my($RegionID, $StartResNum, $EndResNum, $ResNum, $ResIndex, $ID, $Sequence, $Residue);
  for $RegionID (@{$SequenceFilesInfo{RegionsData}[$FileIndex]{RegionIDs}}) {
    $StartResNum = $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{StartResNum};
    $EndResNum = $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{EndResNum};
    RESNUM: for $ResNum ($StartResNum .. $EndResNum) {
      $ResIndex = $ResNum - 1;
      if ($OptionsInfo{IgnoreGaps} && $SequenceFilesInfo{RefereceSequenceResNums}[$FileIndex]{IsGap}{$ResNum}) {
	next RESNUM;
      }
      # Go over residues in column $ResNum in all the sequences...
      ID: for $ID (@{$SequenceDataRef->{IDs}}) {
	$Sequence = $SequenceDataRef->{Sequence}{$ID};
	$Residue = substr($Sequence, $ResIndex, 1);
	if (IsGapResidue($Residue)) {
	  $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{Count}{$ResNum}{Gap} += 1;
	}
	else {
	  if (exists $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{Count}{$ResNum}{$Residue}) {
	    $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{Count}{$ResNum}{$Residue} += 1;
	  }
	  else {
	    # Internal error...
	    print "Warning: Ignoring residue $Residue in sequence $ID during ResidueFrequencyAnalysis calculation: Unknown residue...\n";
	  }
	}
      }
    }
  }
}

# Calculate percent frequency for various residues in the sequence regions...
sub CalculatePercentResidueFrequencyInRegions {
  my($FileIndex, $SequenceDataRef) = @_;
  my($RegionID, $StartResNum, $EndResNum, $ResNum, $Residue, $Count, $PercentCount, $MaxResiduesCount, $Precision);

  $MaxResiduesCount = $SequenceDataRef->{Count};
  $Precision = $OptionsInfo{Precision};
  for $RegionID (@{$SequenceFilesInfo{RegionsData}[$FileIndex]{RegionIDs}}) {
    $StartResNum = $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{StartResNum};
    $EndResNum = $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{EndResNum};
    RESNUM: for $ResNum ($StartResNum .. $EndResNum) {
      if ($OptionsInfo{IgnoreGaps} && $SequenceFilesInfo{RefereceSequenceResNums}[$FileIndex]{IsGap}{$ResNum}) {
	next RESNUM;
      }
      for $Residue (keys %{$SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{Count}{$ResNum}}) {
	$Count = $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{Count}{$ResNum}{$Residue};
	$PercentCount = ($Count / $MaxResiduesCount) * 100;
	$PercentCount = sprintf("%.${Precision}f", $PercentCount) + 0;
	$SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{PercentCount}{$ResNum}{$Residue} = $PercentCount;
      }
    }
  }
}

# Generate output files...
sub GeneratePercentResidueFrequencyOutFilesForRegions {
  my($FileIndex, $SequenceDataRef) = @_;

  # Setup rerfernce sequence data...
  my($RefereceSequenceID, $RefereceSequence);
  $RefereceSequenceID = $SequenceFilesInfo{RefereceSequenceID}[$FileIndex];
  $RefereceSequence = $SequenceFilesInfo{RefereceSequence}[$FileIndex];

  my($RegionID, $StartResNum, $EndResNum, $ResNum, $Count, $PercentCount, $Residue, $RegionNum, $RegionOutFile, $PercentRegionOutFile, $OutFileRoot, $OutFileExt, $Line, @LineWords, @PercentLineWords);

  $OutFileRoot = $SequenceFilesInfo{OutFileRoot}[$FileIndex];
  $OutFileExt = $SequenceFilesInfo{OutFileExt}[$FileIndex];
  $RegionNum = 0;
  for $RegionID (@{$SequenceFilesInfo{RegionsData}[$FileIndex]{RegionIDs}}) {
    $RegionNum++;
    $StartResNum = $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{StartResNum};
    $EndResNum = $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{EndResNum};

    $RegionOutFile = "${OutFileRoot}ResidueFrequencyAnalysisRegion${RegionNum}.${OutFileExt}";
    $PercentRegionOutFile = "${OutFileRoot}PercentResidueFrequencyAnalysisRegion${RegionNum}.${OutFileExt}";

    print "Generating $RegionOutFile and $PercentRegionOutFile...\n";
    open REGIONOUTFILE, ">$RegionOutFile" or die "Error: Can't open $RegionOutFile: $! \n";
    open PERCENTREGIONOUTFILE, ">$PercentRegionOutFile" or die "Error: Can't open $PercentRegionOutFile: $! \n";

    # Write out reference residue positions as column values....
    @LineWords = ();
    push @LineWords, '';
    RESNUM: for $ResNum ($StartResNum .. $EndResNum) {
      if ($OptionsInfo{IgnoreGaps} && $SequenceFilesInfo{RefereceSequenceResNums}[$FileIndex]{IsGap}{$ResNum}) {
	next RESNUM;
      }
      push @LineWords, $ResNum;
    }
    $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print REGIONOUTFILE "$Line\n";
    print PERCENTREGIONOUTFILE "$Line\n";


    # Write out row data for each residue; Gap residue is written last...
    RESIDUE: for $Residue (sort @{$SequenceFilesInfo{ResidueCodes}[$FileIndex]}) {
      if ($Residue =~ /^Gap$/i) {
	next RESIDUE;
      }
      @LineWords = ();
      push @LineWords, $Residue;
      @PercentLineWords = ();
      push @PercentLineWords, $Residue;

      RESNUM: for $ResNum ($StartResNum .. $EndResNum) {
	if ($OptionsInfo{IgnoreGaps} && $SequenceFilesInfo{RefereceSequenceResNums}[$FileIndex]{IsGap}{$ResNum}) {
	  next RESNUM;
	}
	$Count = $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{Count}{$ResNum}{$Residue};
	push @LineWords, $Count;
	$PercentCount = $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{PercentCount}{$ResNum}{$Residue};
	push @PercentLineWords, $PercentCount;
      }
      $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
      print REGIONOUTFILE "$Line\n";

      $Line = JoinWords(\@PercentLineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
      print PERCENTREGIONOUTFILE "$Line\n";
    }

    # Write out data for gap...
    $Residue = 'Gap';
    @LineWords = ();
    push @LineWords, $Residue;
    @PercentLineWords = ();
    push @PercentLineWords, $Residue;

    RESNUM: for $ResNum ($StartResNum .. $EndResNum) {
      if ($OptionsInfo{IgnoreGaps} && $SequenceFilesInfo{RefereceSequenceResNums}[$FileIndex]{IsGap}{$ResNum}) {
	next RESNUM;
      }
      $Count = $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{Count}{$ResNum}{$Residue};
      push @LineWords, $Count;

      $PercentCount = $SequenceFilesInfo{RegionsData}[$FileIndex]{$RegionID}{PercentCount}{$ResNum}{$Residue};
      push @PercentLineWords, $PercentCount;
    }
    $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print REGIONOUTFILE "$Line\n";

    $Line = JoinWords(\@PercentLineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print PERCENTREGIONOUTFILE "$Line\n";

    close REGIONOUTFILE;
    close PERCENTREGIONOUTFILE;
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  # Setup analysis mode...
  $OptionsInfo{CalculatePercentIdentityMatrix} = ($Options{mode} =~ /^(PercentIdentityMatrix|All)$/i) ? 1 : 0;
  $OptionsInfo{PerformResidueFrequencyAnalysis} = ($Options{mode} =~ /^(ResidueFrequencyAnalysis|All)$/i) ? 1 : 0;

  # Setup delimiter and quotes...
  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /tab/i ) ? "\t" : (($Options{outdelim} =~ /semicolon/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /yes/i ) ? 1 : 0;

  # Setup reference sequence and regions for residue frequence analysis...
  $OptionsInfo{SpecifiedRefereceSequence} = $Options{referencesequence};
  $OptionsInfo{SpecifiedRegion} = $Options{region};
  @{$OptionsInfo{SpecifiedRegions}} = ();

  my(@SpecifiedRegions);
  @SpecifiedRegions = ();
  if ($Options{region} =~ /\,/) {
    @SpecifiedRegions = split /\,/, $OptionsInfo{SpecifiedRegion};
    if (@SpecifiedRegions % 2) {
      die "Error: The value specified, $Options{region}, for option \"--region\" is not valid. Allowed values: \"StartResNum,EndResNum,[StartResNum,EndResNum...]\" or UseCompleteSequence\n";
    }
    # Make sure EndResNum > StartResNum...
    my($StartResNum, $EndResNum, $Index, $RegionNum);
    $RegionNum = 0;
    for ($Index = 0; $Index <= $#SpecifiedRegions; $Index += 2) {
      $StartResNum = $SpecifiedRegions[$Index];
      $EndResNum = $SpecifiedRegions[$Index + 1];
      $RegionNum++;
      if (!IsPositiveInteger($StartResNum)) {
	die "Error: The value specified, $Options{region}, for option \"--region\" is not valid: The start residue number, $StartResNum, must be a positive integer for region $RegionNum.\n";
      }
      if (!IsPositiveInteger($EndResNum)) {
	die "Error: The value specified, $Options{region}, for option \"--region\" is not valid: The start residue number, $EndResNum, must be a positive integer for region $RegionNum.\n";
      }
      if ($StartResNum >= $EndResNum) {
	die "Error: The value specified, $Options{region}, for option \"--region\" is not valid: The start residue number, $StartResNum, must be smaller than end residue number, $EndResNum, for region $RegionNum.\n";
      }
    }
  }
  else {
    if ($Options{region} !~ /^UseCompleteSequence$/i) {
      die "Error: The value specified, $Options{region}, for option \"--region\" is not valid. Allowed values: \"StartResNum,EndResNum,[StartResNum,EndResNum...]\" or UseCompleteSequence\n";
    }
  }
  push @{$OptionsInfo{SpecifiedRegions}}, @SpecifiedRegions;

  # Miscellaneous options...
  $OptionsInfo{Precision} = $Options{precision};
  $OptionsInfo{IgnoreGaps} = ($Options{ignoregaps} =~ /Yes/i) ? 1 : 0;
  $OptionsInfo{RegionResiduesMode} = $Options{regionresiduesmode};

  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? 1 : 0;
  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : 0;
}

# Retrieve information about sequence files...
sub RetrieveSequenceFilesInfo {
  my($Index, $SequenceFile, $FileSupported, $FileFormat, $SequenceCount, $RefereceSequence, $RefereceSequenceID, $RefereceSequenceLen, $RefereceSequenceWithNoGaps, $RefereceSequenceWithNoGapsLen, $RefereceSequenceRegionCount, $FileDir, $FileName, $FileExt, $OutFileRoot, $OutFileExt, $SequenceDataRef, $SpecifiedRefereceSequence, @SpecifiedRegions, @RefereceSequenceRegions);

  %SequenceFilesInfo = ();
  @{$SequenceFilesInfo{FilesOkay}} = ();
  @{$SequenceFilesInfo{OutFileRoot}} = ();
  @{$SequenceFilesInfo{OutFileExt}} = ();
  @{$SequenceFilesInfo{Format}} = ();
  @{$SequenceFilesInfo{SequenceCount}} = ();
  @{$SequenceFilesInfo{RefereceSequenceID}} = ();
  @{$SequenceFilesInfo{RefereceSequence}} = ();
  @{$SequenceFilesInfo{RefereceSequenceLen}} = ();
  @{$SequenceFilesInfo{RefereceSequenceWithNoGaps}} = ();
  @{$SequenceFilesInfo{RefereceSequenceWithNoGapsLen}} = ();
  @{$SequenceFilesInfo{RefereceSequenceRegions}} = ();
  @{$SequenceFilesInfo{RefereceSequenceRegionCount}} = ();
  @{$SequenceFilesInfo{ResidueCodes}} = ();

  FILELIST: for $Index (0 .. $#SequenceFilesList) {
    $SequenceFile = $SequenceFilesList[$Index];
    $SequenceFilesInfo{FilesOkay}[$Index] = 0;
    $SequenceFilesInfo{OutFileRoot}[$Index] = '';
    $SequenceFilesInfo{OutFileExt}[$Index] = '';
    $SequenceFilesInfo{Format}[$Index] = 'NotSupported';
    $SequenceFilesInfo{SequenceCount}[$Index] = 0;
    $SequenceFilesInfo{RefereceSequenceID}[$Index] = '';
    $SequenceFilesInfo{RefereceSequence}[$Index] = '';
    $SequenceFilesInfo{RefereceSequenceLen}[$Index] = '';
    $SequenceFilesInfo{RefereceSequenceWithNoGaps}[$Index] = '';
    $SequenceFilesInfo{RefereceSequenceWithNoGapsLen}[$Index] = '';
    @{$SequenceFilesInfo{RefereceSequenceRegions}[$Index]} = ();
    $SequenceFilesInfo{RefereceSequenceRegionCount}[$Index] = 0;
    @{$SequenceFilesInfo{ResidueCodes}[$Index]} = ();

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

    # Make sure all sequence lengths are identical...
    if (!AreSequenceLengthsIdentical($SequenceDataRef)) {
      warn "Warning: Ignoring file $SequenceFile: Sequence legths are not identical.\n";
      next FILELIST;
    }
    $SpecifiedRefereceSequence = $OptionsInfo{SpecifiedRefereceSequence};
    # Make sure reference sequence ID is valid...
    if ($SpecifiedRefereceSequence =~ /^UseFirstSequenceID$/i) {
      $RefereceSequenceID = $SequenceDataRef->{IDs}[0];
    }
    else {
      if (!exists($SequenceDataRef->{Sequence}{$SpecifiedRefereceSequence})) {
	warn "Warning: Ignoring file $SequenceFile: Rreference sequence ID, $SpecifiedRefereceSequence, specified using option \"--referencesequence\" is missing.\n";
	next FILELIST;
      }
      $RefereceSequenceID = $SpecifiedRefereceSequence;
    }

    # Make sure sequence regions corresponding to reference sequence are valid...
    @RefereceSequenceRegions = ();
    $RefereceSequenceRegionCount = 0;
    $RefereceSequence = $SequenceDataRef->{Sequence}{$RefereceSequenceID};
    $RefereceSequenceLen = length $RefereceSequence;

    $RefereceSequenceWithNoGaps = RemoveSequenceGaps($RefereceSequence);
    $RefereceSequenceWithNoGapsLen = length $RefereceSequenceWithNoGaps;

    @SpecifiedRegions = ();
    push @SpecifiedRegions, @{$OptionsInfo{SpecifiedRegions}};
    if (@SpecifiedRegions) {
      # Make sure specified regions are valid...
      my($StartResNum, $EndResNum, $RegionIndex, $RegionNum);
      $RegionNum = 0;
      for ($RegionIndex = 0; $RegionIndex <= $#SpecifiedRegions; $RegionIndex += 2) {
	$StartResNum = $SpecifiedRegions[$RegionIndex];
	$EndResNum = $SpecifiedRegions[$RegionIndex + 1];
	$RegionNum++;
	if ($OptionsInfo{IgnoreGaps}) {
	  if ($StartResNum > $RefereceSequenceWithNoGapsLen) {
	    warn "Warning: Ignoring file $SequenceFile: The value specified, $Options{region}, for option \"--region\" is not valid: The start residue number, $StartResNum, must be smaller the sequence length, $RefereceSequenceWithNoGapsLen, of reference sequence ID,  $RefereceSequenceID, in region $RegionNum. The reference sequence residue numbers correspond to the sequence with no gaps. Specify \"No\" value for \"-i, --ignoregaps\" option to use residue numbers corresponding to reference sequence including gaps.\n";
	    next FILELIST;
	  }
	  if ($EndResNum > $RefereceSequenceWithNoGapsLen) {
	    warn "Warning: Ignoring file $SequenceFile: The value specified, $Options{region}, for option \"--region\" is not valid: The end residue number, $EndResNum, must be smaller the sequence length, $RefereceSequenceWithNoGapsLen, of reference sequence ID,  $RefereceSequenceID, in region $RegionNum. The reference sequence residue numbers correspond to the sequence with no gaps. Specify \"No\" value for \"-i, --ignoregaps\" option to use residue numbers corresponding to reference sequence including gaps.\n";
	    next FILELIST;
	  }
	}
	else {
	  if ($StartResNum > $RefereceSequenceLen) {
	    warn "Warning: Ignoring file $SequenceFile: The value specified, $Options{region}, for option \"--region\" is not valid: The start residue number, $StartResNum, must be smaller the sequence length, $RefereceSequenceLen, of reference sequence ID,  $RefereceSequenceID, in region $RegionNum.\n";
	    next FILELIST;
	  }
	  if ($EndResNum > $RefereceSequenceLen) {
	    warn "Warning: Ignoring file $SequenceFile: The value specified, $Options{region}, for option \"--region\" is not valid: The end residue number, $EndResNum, must be smaller the sequence length, $RefereceSequenceLen, of reference sequence ID,  $RefereceSequenceID, in region $RegionNum.\n";
	    next FILELIST;
	  }
	}
	push @RefereceSequenceRegions, ($StartResNum, $EndResNum);
      }
      $RefereceSequenceRegionCount = $RegionNum;
    }
    else {
      # Use complete sequence corresponding to reference sequence ID...
      if ($OptionsInfo{IgnoreGaps}) {
	push @RefereceSequenceRegions, (1, $RefereceSequenceWithNoGapsLen);
      }
      else {
	push @RefereceSequenceRegions, (1, $RefereceSequenceLen);
      }
      $RefereceSequenceRegionCount = 1;
    }
    # Setup output file names...
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($SequenceFile);
    $FileExt = "csv";
    if ($Options{outdelim} =~ /^tab$/i) {
      $FileExt = "tsv";
    }
    $OutFileExt = $FileExt;
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
    if (!$OptionsInfo{OverwriteFiles}) {
      if ($OptionsInfo{CalculatePercentIdentityMatrix}) {
	if (-e "${OutFileRoot}PercentIdentityMatrix.${OutFileExt}") {
	  warn "Warning: Ignoring file $SequenceFile: The file ${OutFileRoot}PercentIdentityMatrix.${OutFileExt} already exists\n";
	  next FILELIST;
	}
      }
      if ($OptionsInfo{PerformResidueFrequencyAnalysis}) {
	my($RegionNum);
	for $RegionNum (1 .. $RefereceSequenceRegionCount) {
	  if (-e "${OutFileRoot}ResidueFrequencyAnalysisRegion${RegionNum}.${OutFileExt}") {
	    warn "Warning: Ignoring file $SequenceFile: The file ${OutFileRoot}ResidueFrequencyAnalysisRegion${RegionNum}.${OutFileExt} already exists\n";
	    next FILELIST;
	  }
	  if (-e "${OutFileRoot}PercentResidueFrequencyAnalysisRegion${RegionNum}.${OutFileExt}") {
	    warn "Warning: Ignoring file $SequenceFile: The file ${OutFileRoot}PercentResidueFrequencyAnalysisRegion${RegionNum}.${OutFileExt} already exists\n";
	    next FILELIST;
	  }
	}
      }
    }

    $SequenceFilesInfo{FilesOkay}[$Index] = 1;
    $SequenceFilesInfo{OutFileRoot}[$Index] = $OutFileRoot;
    $SequenceFilesInfo{OutFileExt}[$Index] = $OutFileExt;

    $SequenceFilesInfo{Format}[$Index] = $FileFormat;
    $SequenceFilesInfo{SequenceCount}[$Index] = $SequenceCount;
    $SequenceFilesInfo{RefereceSequenceID}[$Index] = $RefereceSequenceID;
    $SequenceFilesInfo{RefereceSequence}[$Index] = $RefereceSequence;
    $SequenceFilesInfo{RefereceSequenceLen}[$Index] = $RefereceSequenceLen;
    $SequenceFilesInfo{RefereceSequenceWithNoGaps}[$Index] = $RefereceSequenceWithNoGaps;
    $SequenceFilesInfo{RefereceSequenceWithNoGapsLen}[$Index] = $RefereceSequenceWithNoGapsLen;
    push @{$SequenceFilesInfo{RefereceSequenceRegions}[$Index]}, @RefereceSequenceRegions;
    $SequenceFilesInfo{RefereceSequenceRegionCount}[$Index] = $RefereceSequenceRegionCount;

    # Setup residue codes...
    SetupSequenceFileResidueCodes($SequenceDataRef, $Index);
  }
}

sub SetupSequenceFileResidueCodes {
  my($SequenceDataRef, $FileIndex) = @_;
  my($Residue, @ResidueCodesList, %ResidueCodesMap);

  # Initialize
  @{$SequenceFilesInfo{ResidueCodes}[$FileIndex]} = ();

  %ResidueCodesMap = ();
  @ResidueCodesList = ();

  # Setup standard amino acids and nucleic acids one letter codes...
  if ($OptionsInfo{RegionResiduesMode} =~ /^AminoAcids$/i) {
    @ResidueCodesList = AminoAcids::GetAminoAcids('OneLetterCode');
  }
  elsif ($OptionsInfo{RegionResiduesMode} =~ /^NucleicAcids$/i) {
    push @ResidueCodesList, ('A', 'G', 'T', 'U', 'C');
  }
  push @ResidueCodesList, 'Gap';
  for $Residue (@ResidueCodesList) {
    $ResidueCodesMap{$Residue} = $Residue;
  }

  # Go over all the residues in all the sequences and add missing ones to the list...
  my($ID, $Sequence, $ResIndex);
  for $ID (@{$SequenceDataRef->{IDs}}) {
    $Sequence = $SequenceDataRef->{Sequence}{$ID};
    RES: for $ResIndex (0 .. (length($Sequence) - 1)) {
      $Residue = substr($Sequence, $ResIndex, 1);
      if (IsGapResidue($Residue)) {
	next RES;
      }
      if (exists $ResidueCodesMap{$Residue}) {
	next RES;
      }
      push @ResidueCodesList, $Residue;
      $ResidueCodesMap{$Residue} = $Residue;
    }
  }
  push @{$SequenceFilesInfo{ResidueCodes}[$FileIndex]}, @ResidueCodesList;
}

# Setup regions data for performing residue frequency analysis...
sub SetupSequenceRegionsData {
  my($Index, $RefereceSequence, $RefereceSequenceLen, $RegionID, $StartResNum, $EndResNum, $RegionIndex, $RegionNum, $NoGapResNum, $ResNum, $ResIndex, $Residue, $ResidueCode, @RefereceSequenceRegions);


  @{$SequenceFilesInfo{RefereceSequenceResNums}} = ();
  @{$SequenceFilesInfo{RegionsData}} = ();

  FILELIST: for $Index (0 .. $#SequenceFilesList) {
    %{$SequenceFilesInfo{RefereceSequenceResNums}[$Index]{IsGap}} = ();
    %{$SequenceFilesInfo{RefereceSequenceResNums}[$Index]{NoGapToGap}} = ();
    %{$SequenceFilesInfo{RefereceSequenceResNums}[$Index]{GapToNoGap}} = ();
    %{$SequenceFilesInfo{RegionsData}[$Index]} = ();
    @{$SequenceFilesInfo{RegionsData}[$Index]{RegionIDs}} = ();

    if (!$SequenceFilesInfo{FilesOkay}[$Index]) {
      next FILELIST;
    }
    if (!$OptionsInfo{PerformResidueFrequencyAnalysis}) {
      next FILELIST;
    }

    $RefereceSequence = $SequenceFilesInfo{RefereceSequence}[$Index];
    $RefereceSequenceLen = $SequenceFilesInfo{RefereceSequenceLen}[$Index];

    # Setup residue number mapping and gap status for refernece sequence...
    $NoGapResNum = 0;
    $ResNum = 0;
    for $ResIndex (0 .. ($RefereceSequenceLen - 1)) {
      $ResNum++;
      $Residue = substr($RefereceSequence, $ResIndex, 1);
      $SequenceFilesInfo{RefereceSequenceResNums}[$Index]{IsGap}{$ResNum} = 1;
      if (!IsGapResidue($Residue)) {
	$NoGapResNum++;
	$SequenceFilesInfo{RefereceSequenceResNums}[$Index]{IsGap}{$ResNum} = 0;
	$SequenceFilesInfo{RefereceSequenceResNums}[$Index]{NoGapToGap}{$NoGapResNum} = $ResNum;
	$SequenceFilesInfo{RefereceSequenceResNums}[$Index]{GapToNoGap}{$ResNum} = $NoGapResNum;
      }
    }
    # Map residue numbers for specified regions to the reference residue in input sequence/alignment files
    $RegionNum = 0;
    @RefereceSequenceRegions = ();
    push @RefereceSequenceRegions, @{$SequenceFilesInfo{RefereceSequenceRegions}[$Index]};
    for ($RegionIndex = 0; $RegionIndex <= $#RefereceSequenceRegions; $RegionIndex += 2) {
      $StartResNum = $RefereceSequenceRegions[$RegionIndex];
      $EndResNum = $RefereceSequenceRegions[$RegionIndex + 1];
      $RegionNum++;
      $RegionID = "Region${RegionNum}";
      if ($OptionsInfo{IgnoreGaps}) {
	# Map residue numbers to the reference sequence residue numbers to account for any ignored gaps...
	$StartResNum = $SequenceFilesInfo{RefereceSequenceResNums}[$Index]{NoGapToGap}{$StartResNum};
	$EndResNum = $SequenceFilesInfo{RefereceSequenceResNums}[$Index]{NoGapToGap}{$EndResNum};
      }
      push @{$SequenceFilesInfo{RegionsData}[$Index]{RegionIDs}}, $RegionID;
      $SequenceFilesInfo{RegionsData}[$Index]{$RegionID}{StartResNum} = $StartResNum;
      $SequenceFilesInfo{RegionsData}[$Index]{$RegionID}{EndResNum} = $EndResNum;

      # Initialize data for residue codes...
      for $ResNum ($StartResNum .. $EndResNum) {
	for $ResidueCode (@{$SequenceFilesInfo{ResidueCodes}[$Index]}) {
	  $SequenceFilesInfo{RegionsData}[$Index]{$RegionID}{Count}{$ResNum}{$ResidueCode} = 0;
	}
      }
    }
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{ignoregaps} = 'yes';
  $Options{regionresiduesmode} = 'None';
  $Options{mode} = 'PercentIdentityMatrix';
  $Options{outdelim} = 'comma';
  $Options{precision} = 2;
  $Options{quote} = 'yes';
  $Options{referencesequence} = 'UseFirstSequenceID';
  $Options{region} = 'UseCompleteSequence';

  if (!GetOptions(\%Options, "help|h", "ignoregaps|i=s", "mode|m=s", "outdelim=s", "overwrite|o", "precision|p=i", "quote|q=s", "referencesequence=s", "region=s", "regionresiduesmode=s", "root|r=s", "workingdir|w=s")) {
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
  if ($Options{regionresiduesmode} !~ /^(AminoAcids|NucleicAcids|None)$/i) {
    die "Error: The value specified, $Options{regionresiduesmode}, for option \"--regionresiduesmode\" is not valid. Allowed values: AminoAcids, NucleicAcids or None\n";
  }
  if ($Options{mode} !~ /^(PercentIdentityMatrix|ResidueFrequencyAnalysis|All)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: PercentIdentityMatrix, ResidueFrequencyAnalysis  or All\n";
  }
  if ($Options{outdelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{outdelim}, for option \"--outdelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if ($Options{quote} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not valid. Allowed values: yes or no\n";
  }
  if (!IsPositiveInteger($Options{precision})) {
    die "Error: The value specified, $Options{precision}, for option \"-p --precision\" is not valid. Allowed values: > 0 \n";
  }
}

__END__

=head1 NAME

AnalyzeSequenceFilesData.pl - Analyze sequence and alignment files

=head1 SYNOPSIS

AnalyzeSequenceFilesData.pl SequenceFile(s) AlignmentFile(s)...

AnalyzeSequenceFilesData.pl [B<-h, --help>] [B<-i, --IgnoreGaps> yes | no]
[B<-m, --mode> PercentIdentityMatrix | ResidueFrequencyAnalysis | All]
[B<--outdelim> comma | tab | semicolon] [B<-o, --overwrite>] [B<-p, --precision> number] [B<-q, --quote> yes | no]
[B<--ReferenceSequence> SequenceID | UseFirstSequenceID]
[B<--region> "StartResNum, EndResNum, [StartResNum, EndResNum...]" | UseCompleteSequence]
[B<--RegionResiduesMode> AminoAcids | NucleicAcids | None]
[B<-w, --WorkingDir> dirname] SequenceFile(s) AlignmentFile(s)...

=head1 DESCRIPTION

Analyze I<SequenceFile(s) and AlignmentFile(s)> data: calculate pairwise percent identity matrix or
calculate percent occurrence of various residues in specified sequence regions. All the sequences
in the input file must have the same sequence lengths; otherwise, the sequence file is ignored.

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

Ignore gaps during calculation of sequence lengths and specification of regions during residue
frequency analysis. Possible values: I<yes or no>. Default value: I<yes>.

=item B<-m, --mode> I<PercentIdentityMatrix | ResidueFrequencyAnalysis | All>

Specify how to analyze data in sequence files: calculate percent identity matrix or calculate
frequency of occurrence of residues in specific regions. During I<ResidueFrequencyAnalysis> value
of B<-m, --mode> option, output files are generated for both the residue count and percent residue
count. Possible values: I<PercentIdentityMatrix, ResidueFrequencyAnalysis, or All>. Default value:
I<PercentIdentityMatrix>.

=item B<--outdelim> I<comma | tab | semicolon>

Output text file delimiter. Possible values: I<comma, tab, or semicolon>.
Default value: I<comma>.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-p, --precision> I<number>

Precision of calculated values in the output file. Default: up to I<2> decimal places.
Valid values: positive integers.

=item B<-q, --quote> I<yes | no>

Put quotes around column values in output text file. Possible values: I<yes or
no>. Default value: I<yes>.

=item B<--ReferenceSequence> I<SequenceID | UseFirstSequenceID>

Specify reference sequence ID to identify regions for performing I<ResidueFrequencyAnalysis> specified
using B<-m, --mode> option. Default: I<UseFirstSequenceID>.

=item B<--region> I<StartResNum,EndResNum,[StartResNum,EndResNum...] | UseCompleteSequence>

Specify how to perform frequency of occurrence analysis for residues: use specific regions
indicated by starting and ending residue numbers in reference sequence or use the whole reference
sequence as one region. Default: I<UseCompleteSequence>.

Based on the value of B<-i, --IgnoreGaps> option, specified residue numbers I<StartResNum,EndResNum>
correspond to the positions in the reference sequence without gaps or with gaps.

For residue numbers corresponding to the reference sequence including gaps, percent occurrence
of various residues corresponding to gap position in reference sequence is also calculated.

=item B<--RegionResiduesMode> I<AminoAcids | NucleicAcids | None>

Specify how to process residues in the regions specified using B<--region> option during
I<ResidueFrequencyAnalysis> calculation: categorize residues as amino acids, nucleic acids, or simply
ignore residue category during the calculation. Possible values: I<AminoAcids, NucleicAcids or None>.
Default value: I<None>.

For I<AminoAcids> or I<NucleicAcids> values of B<--RegionResiduesMode> option, all the standard amino
acids or nucleic acids are listed in the output file for each region; Any gaps and other non standard residues
are added to the list as encountered.

For I<None> value of B<--RegionResiduesMode> option, no assumption is made about type of residues.
Residue and gaps are added to the list as encountered.

=item B<-r, --root> I<rootname>

New sequence file name is generated using the root: <Root><Mode>.<Ext> and
<Root><Mode><RegionNum>.<Ext>. Default new file
name: <SequenceFileName><Mode>.<Ext> for I<PercentIdentityMatrix> value B<m, --mode> option
and <SequenceFileName><Mode><RegionNum>.<Ext>  for I<ResidueFrequencyAnalysis>.
The csv, and tsv <Ext> values are used for comma/semicolon, and tab delimited text
files respectively. This option is ignored for multiple input files.

=item B<-w --WorkingDir> I<text>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To calculate percent identity matrix for all sequences in Sample1.msf file and generate
Sample1PercentIdentityMatrix.csv, type:

    % AnalyzeSequenceFilesData.pl Sample1.msf

To perform residue frequency analysis for all sequences in Sample1.aln file corresponding to
non-gap positions in the first sequence and generate Sample1ResidueFrequencyAnalysisRegion1.csv
and Sample1PercentResidueFrequencyAnalysisRegion1.csv files, type:

    % AnalyzeSequenceFilesData.pl -m ResidueFrequencyAnalysis -o
      Sample1.aln

To perform residue frequency analysis for all sequences in Sample1.aln file corresponding to
all positions in the first sequence and generate TestResidueFrequencyAnalysisRegion1.csv
and TestPercentResidueFrequencyAnalysisRegion1.csv files, type:

    % AnalyzeSequenceFilesData.pl -m ResidueFrequencyAnalysis --IgnoreGaps
      No -o -r Test Sample1.aln

To perform residue frequency analysis for all sequences in Sample1.aln file corresponding to
non-gap residue positions 5 to 10, and 30 to 40 in sequence ACHE_BOVIN and generate
Sample1ResidueFrequencyAnalysisRegion1.csv, Sample1ResidueFrequencyAnalysisRegion2.csv,
SamplePercentResidueFrequencyAnalysisRegion1.csv, and
SamplePercentResidueFrequencyAnalysisRegion2.csv files, type:

    % AnalyzeSequenceFilesData.pl -m ResidueFrequencyAnalysis
      --ReferenceSequence ACHE_BOVIN --region "5,15,30,40" -o Sample1.msf


=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

ExtractFromSequenceFiles.pl, InfoSequenceFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
