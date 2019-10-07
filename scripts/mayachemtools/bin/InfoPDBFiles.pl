#!/usr/bin/perl -w
#
# File: InfoPDBFiles.pl
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

my(@PDBFilesList);
@PDBFilesList = ExpandFileNames(\@ARGV, "pdb");

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Setup information about input files...
my(%PDBFilesInfo);
print "Checking input PDB file(s)...\n";
RetrievePDBFilesInfo();

# Process input files..
my($FileIndex);
if (@PDBFilesList > 1) {
  print "\nProcessing PDB files...\n";
}
for $FileIndex (0 .. $#PDBFilesList) {
  if ($PDBFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $PDBFilesList[$FileIndex]...\n";
    ListPDBFileInfo($FileIndex);
  }
}
ListTotalSizeOfFiles();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# List appropriate information...
sub ListPDBFileInfo {
  my($Index) = @_;
  my($PDBFile, $PDBRecordLinesRef);

  $PDBFile = $PDBFilesList[$Index];
  $PDBRecordLinesRef = ReadPDBFile($PDBFile);

  # Header informaton...
  if ($OptionsInfo{ListHeaderInfo}) {
    ListHeaderInfo($PDBRecordLinesRef);
  }

  # Experiment informaton...
  if ($OptionsInfo{ListExperimentalTechniqueInfo}) {
    ListExperimentalTechniqueInfo($PDBRecordLinesRef);
  }

  # Total number of records...
  my($TotalRecordsCount) = scalar @{$PDBRecordLinesRef};
  print "\nTotal number of records: $TotalRecordsCount\n";

  # List record type count information...
  ListRecordTypesInfo($PDBRecordLinesRef);

  if ($OptionsInfo{CountChains} || $OptionsInfo{CountResiduesInChains} || $OptionsInfo{ResiduesFrequencyInChains}) {
    ListChainsAndResiduesInfo($PDBRecordLinesRef);
  }
  if ($OptionsInfo{CountResiduesAll} || $OptionsInfo{ResiduesFrequencyAll}) {
    ListAllResiduesInfo($PDBRecordLinesRef);
  }
  if ($OptionsInfo{ResidueNumbersInfo}) {
    ListResidueNumbersInfo($PDBRecordLinesRef);
  }
  if ($OptionsInfo{CalculateBoundingBox}) {
    ListBoundingBox($PDBRecordLinesRef);
  }

  # File size and modification information...
  print "\nFile size: ", FormatFileSize($PDBFilesInfo{FileSize}[$Index]), " \n";
  print "Last modified: ", $PDBFilesInfo{FileLastModified}[$Index], " \n";
}

sub ListHeaderInfo {
  my($PDBRecordLinesRef) = @_;
  my($HeaderRecordLine, $Classification, $DepositionDate, $IDCode);

  ($Classification, $DepositionDate, $IDCode) = (undef) x 3;
  $HeaderRecordLine = $PDBRecordLinesRef->[0];
  if (IsHeaderRecordType($HeaderRecordLine)) {
    ($Classification, $DepositionDate, $IDCode) = ParseHeaderRecordLine($HeaderRecordLine);
  }

  $Classification = IsEmpty($Classification) ? 'Not available' : $Classification;
  $DepositionDate = IsEmpty($DepositionDate) ? 'Not available' : $DepositionDate;
  $IDCode = IsEmpty($IDCode) ? 'Not available' : $IDCode;

  print "\nClassification: $Classification\nID: $IDCode\nDeposition date: $DepositionDate\n";
}

# List experimental technique information info...
sub ListExperimentalTechniqueInfo {
  my($PDBRecordLinesRef) = @_;
  my($ExperimentalTechnique, $Resolution, $ResolutionUnits);

  $ExperimentalTechnique = GetExperimentalTechnique($PDBRecordLinesRef);
  print "\nExperimental technique: " . ($ExperimentalTechnique ? $ExperimentalTechnique : "Not available") . "\n";

  ($Resolution, $ResolutionUnits) = GetExperimentalTechniqueResolution($PDBRecordLinesRef);
  print "Resolution: " . ($Resolution ? "$Resolution $ResolutionUnits" : "Not available") . "\n";

}

# List record type info...
sub ListRecordTypesInfo {
  my($PDBRecordLinesRef) = @_;
  my($RecordType, $RecordCount, $RecordTypesCountRef, @RecordTypeCountInfo);

  $RecordTypesCountRef = GetRecordTypesCount($PDBRecordLinesRef);

  @RecordTypeCountInfo = ();
  if ($OptionsInfo{CountRecordType} =~ /^All$/i) {
    for $RecordType (@{$RecordTypesCountRef->{RecordTypes}}) {
      $RecordCount = $RecordTypesCountRef->{Count}{$RecordType};
      push @RecordTypeCountInfo, "$RecordType - $RecordCount";
    }
  }
  else {
    for $RecordType (@{$OptionsInfo{SpecifiedRecordTypes}}) {
      $RecordCount = (exists $RecordTypesCountRef->{Count}{$RecordType}) ? ($RecordTypesCountRef->{Count}{$RecordType}) : 0;
      push @RecordTypeCountInfo, "$RecordType - $RecordCount";
    }
  }
  print "Number of individual records: ", JoinWords(\@RecordTypeCountInfo, '; ', 0), "\n";

  if ($OptionsInfo{CheckMasterRecord}) {
    CheckMasterRecord($RecordTypesCountRef, $PDBRecordLinesRef);
  }
}

# List information about residues and chains...
sub ListChainsAndResiduesInfo {
  my($PDBRecordLinesRef) = @_;
  my($ResidueName, $ResidueCount, $ChainCount, $ChainID, $CollectChainResiduesBeyondTER, $ChainsAndResiduesInfoRef);

  $CollectChainResiduesBeyondTER = 1;
  $ChainsAndResiduesInfoRef = GetChainsAndResidues($PDBRecordLinesRef, 'AtomAndHetatm', $CollectChainResiduesBeyondTER);
  $ChainCount = @{$ChainsAndResiduesInfoRef->{ChainIDs}};
  if ($OptionsInfo{CountChains}) {
    print "\nNumber of chains: $ChainCount \n";
    my($ChainID, @ChainIDsList);
    @ChainIDsList = ();
    for $ChainID (@{$ChainsAndResiduesInfoRef->{ChainIDs}}) {
      push @ChainIDsList, CleanupChainID($ChainID);
    }
    print "Chain IDs: ", JoinWords(\@ChainIDsList, ', ', 0),"\n";
  }

  if ($OptionsInfo{CountResiduesInChains}) {
    my($TotalResiduesCount, $ResidueCountInfo, @ResiduesCountInfo);
    @ResiduesCountInfo = ();
    $TotalResiduesCount = 0;
    for $ChainID (@{$ChainsAndResiduesInfoRef->{ChainIDs}}) {
      $ResidueCount = @{$ChainsAndResiduesInfoRef->{Residues}{$ChainID}};
      $TotalResiduesCount += $ResidueCount;
      $ResidueCountInfo =  "Chain ${ChainID} - ${ResidueCount}";
      push @ResiduesCountInfo, $ResidueCountInfo;
    }
    print "\nNumber of residues in chain(s): ";
    if ($ChainCount > 1) {
      print "Total - $TotalResiduesCount; ", JoinWords(\@ResiduesCountInfo, '; ', 0),"\n";
    }
    else {
      print "$TotalResiduesCount\n";
    }

    # List of residues in each chain...
    if ($OptionsInfo{DetailLevel} >= 3) {
      print "List of residues in chain(s): \n";
      for $ChainID (@{$ChainsAndResiduesInfoRef->{ChainIDs}}) {
	if ($ChainCount > 1) {
	  print "Chain ", CleanupChainID($ChainID), ": ";
	}
	print JoinWords(\@{$ChainsAndResiduesInfoRef->{Residues}{$ChainID}}, ', ', 0),"\n";
      }
    }
  }
  if ($OptionsInfo{ResiduesFrequencyInChains}) {
    # Setup a hash using residue count as key for sorting the values...
    my(%ResiduesCountToNameMap);
    %ResiduesCountToNameMap = ();
    @{$ResiduesCountToNameMap{ChainIDs}} = ();
    %{$ResiduesCountToNameMap{ResidueNames}} = ();

    for $ChainID (@{$ChainsAndResiduesInfoRef->{ChainIDs}}) {
      push @{$ResiduesCountToNameMap{ChainIDs}}, $ChainID;
      %{$ResiduesCountToNameMap{ResidueNames}{$ChainID}} = ();

      for $ResidueName (sort keys %{$ChainsAndResiduesInfoRef->{ResidueCount}{$ChainID}}) {
	$ResidueCount = $ChainsAndResiduesInfoRef->{ResidueCount}{$ChainID}{$ResidueName};
	# Setup count value for each chain...
	if (exists $ResiduesCountToNameMap{ResidueNames}{$ChainID}{$ResidueCount}) {
	  $ResiduesCountToNameMap{ResidueNames}{$ChainID}{$ResidueCount} .= "~${ResidueName}";
	}
	else {
	  $ResiduesCountToNameMap{ResidueNames}{$ChainID}{$ResidueCount} = $ResidueName;
	}
      }
    }
    # Collect data for all the residues in all the chains...
    my(%AllResiduesNameToCountMap, %AllResiduesCountToNameMap);
    %AllResiduesNameToCountMap = ();
    %AllResiduesCountToNameMap = ();
    if ($ChainCount > 1) {
      for $ChainID (@{$ChainsAndResiduesInfoRef->{ChainIDs}}) {
	for $ResidueName (keys %{$ChainsAndResiduesInfoRef->{ResidueCount}{$ChainID}}) {
	  $ResidueCount = $ChainsAndResiduesInfoRef->{ResidueCount}{$ChainID}{$ResidueName};
	  if (exists $AllResiduesNameToCountMap{$ResidueName}) {
	    $AllResiduesNameToCountMap{$ResidueName} += $ResidueCount;
	  }
	  else {
	    $AllResiduesNameToCountMap{$ResidueName} = $ResidueCount;
	  }
	}
      }
      for $ResidueName (keys %AllResiduesNameToCountMap) {
	$ResidueCount = $AllResiduesNameToCountMap{$ResidueName};
	if (exists $AllResiduesCountToNameMap{$ResidueCount}) {
	  $AllResiduesCountToNameMap{$ResidueCount} .= "~${ResidueName}";
	}
	else {
	  $AllResiduesCountToNameMap{$ResidueCount} = $ResidueName;
	}
      }
    }

    # Setup distribution data for individual chains and the grand total as well...
    my($ChainResidueCount, $PercentResidueCount, $TotalResidueCount, $ResidueNames, @ResidueNamesList, %ResiduesFrequencyInfoMap);
    @{$ResiduesFrequencyInfoMap{ChainIDs}} = ();
    %{$ResiduesFrequencyInfoMap{Frequency}} = ();
    %{$ResiduesFrequencyInfoMap{PercentFrequency}} = ();

    @{$ResiduesFrequencyInfoMap{AllChainsFrequency}} = ();
    @{$ResiduesFrequencyInfoMap{AllChainsPercentFrequency}} = ();

    $TotalResidueCount = 0;

    for $ChainID (@{$ResiduesCountToNameMap{ChainIDs}}) {
      push @{$ResiduesFrequencyInfoMap{ChainIDs}}, $ChainID;
      @{$ResiduesFrequencyInfoMap{Frequency}{$ChainID}} = ();
      @{$ResiduesFrequencyInfoMap{PercentFrequency}{$ChainID}} = ();

      $ChainResidueCount = @{$ChainsAndResiduesInfoRef->{Residues}{$ChainID}};
      $TotalResidueCount += $ChainResidueCount;

      for $ResidueCount (sort {$b <=> $a} keys %{$ResiduesCountToNameMap{ResidueNames}{$ChainID}}) {
	$ResidueNames = $ResiduesCountToNameMap{ResidueNames}{$ChainID}{$ResidueCount};
	@ResidueNamesList = split /~/, $ResidueNames;
	for $ResidueName (@ResidueNamesList) {
	  push @{$ResiduesFrequencyInfoMap{Frequency}{$ChainID}}, "${ResidueName} - ${ResidueCount}";
	  $PercentResidueCount = sprintf("%.1f", (($ResidueCount/$ChainResidueCount)*100)) + 0;
	  push @{$ResiduesFrequencyInfoMap{PercentFrequency}{$ChainID}}, "${ResidueName} - ${PercentResidueCount}%";
	}
      }
    }
    if ($ChainCount > 1) {
      for $ResidueCount (sort {$b <=> $a} keys %AllResiduesCountToNameMap) {
	$ResidueNames = $AllResiduesCountToNameMap{$ResidueCount};
	@ResidueNamesList = split /~/, $ResidueNames;
	for $ResidueName (@ResidueNamesList) {
	  push @{$ResiduesFrequencyInfoMap{AllChainsFrequency}}, "${ResidueName} - ${ResidueCount}";
	  $PercentResidueCount = sprintf("%.1f", (($ResidueCount/$TotalResidueCount)*100)) + 0;
	  push @{$ResiduesFrequencyInfoMap{AllChainsPercentFrequency}}, "${ResidueName} - ${PercentResidueCount}%";
	}
      }
    }

    # List distribution of residues
    print "\nDistribution of residues in chain(s): \n";
    for $ChainID (@{$ResiduesFrequencyInfoMap{ChainIDs}}) {
      if ($ChainCount > 1) {
	print "Chain ", CleanupChainID($ChainID), ": ";
      }
      print JoinWords(\@{$ResiduesFrequencyInfoMap{Frequency}{$ChainID}}, '; ', 0), "\n";
    }
    if ($OptionsInfo{DetailLevel} >= 2) {
      print "\nPercent distribution of residues in chain(s): \n";
      for $ChainID (@{$ResiduesFrequencyInfoMap{ChainIDs}}) {
	if ($ChainCount > 1) {
	  print "Chain ", CleanupChainID($ChainID), ": ";
	}
	print JoinWords(\@{$ResiduesFrequencyInfoMap{PercentFrequency}{$ChainID}}, '; ', 0), "\n";
      }
    }
    if ($ChainCount > 1) {
      print "\nDistribution of residues across all chains: \n";
      print JoinWords(\@{$ResiduesFrequencyInfoMap{AllChainsFrequency}}, '; ', 0), "\n";

      if ($OptionsInfo{DetailLevel} >= 2) {
	print "\nPercent distribution of residues across all chains: \n";
	print JoinWords(\@{$ResiduesFrequencyInfoMap{AllChainsPercentFrequency}}, '; ', 0), "\n";
      }
    }
  }
}

# List information about all the residues...
sub ListAllResiduesInfo {
  my($PDBRecordLinesRef) = @_;
  my($TotalResidueCount, $AtomResiduesCount, $HetatmResiduesCount, $ResiduesInfoRef);

  $ResiduesInfoRef = GetAllResidues($PDBRecordLinesRef);
  $TotalResidueCount = @{$ResiduesInfoRef->{ResidueNames}};
  $AtomResiduesCount = @{$ResiduesInfoRef->{AtomResidueNames}};
  $HetatmResiduesCount = @{$ResiduesInfoRef->{HetatmResidueNames}};

  if ($OptionsInfo{CountResiduesAll}) {
    print "\nTotal number of residues: Total - $TotalResidueCount; ATOM residues - $AtomResiduesCount; HETATM residues - $HetatmResiduesCount\n";

    if ($OptionsInfo{DetailLevel} >= 3) {
      print "List of residues: \n";
      if ($AtomResiduesCount) {
	print "ATOM residues: ", JoinWords(\@{$ResiduesInfoRef->{AtomResidueNames}}, ', ', 0), "\n";
      }
      if ($HetatmResiduesCount) {
	print "HETATM residues: ", JoinWords(\@{$ResiduesInfoRef->{HetatmResidueNames}}, ', ', 0), "\n";
      }
    }
  }

  if ($OptionsInfo{ResiduesFrequencyAll}) {
    my($ResidueName, $ResidueCount);

    # Setup a hash using residue count as key for sorting the values...
    my(%ResiduesCountToNameMap, %AtomResiduesCountToNameMap, %HetatmResiduesCountToNameMap);
    %ResiduesCountToNameMap = ();
    %{$ResiduesCountToNameMap{ResidueNames}} = ();

    %AtomResiduesCountToNameMap = ();
    %{$AtomResiduesCountToNameMap{ResidueNames}} = ();

    %HetatmResiduesCountToNameMap = ();
    %{$HetatmResiduesCountToNameMap{ResidueNames}} = ();

    for $ResidueName (keys %{$ResiduesInfoRef->{ResidueCount}}) {
      $ResidueCount = $ResiduesInfoRef->{ResidueCount}{$ResidueName};
      if (exists $ResiduesCountToNameMap{ResidueNames}{$ResidueCount}) {
	$ResiduesCountToNameMap{ResidueNames}{$ResidueCount} .= "~${ResidueName}";
      }
      else {
	$ResiduesCountToNameMap{ResidueNames}{$ResidueCount} = $ResidueName;
      }
    }

    if ($OptionsInfo{DetailLevel} >= 1) {
      for $ResidueName (keys %{$ResiduesInfoRef->{AtomResidueCount}}) {
	$ResidueCount = $ResiduesInfoRef->{AtomResidueCount}{$ResidueName};
	if (exists $AtomResiduesCountToNameMap{ResidueNames}{$ResidueCount}) {
	  $AtomResiduesCountToNameMap{ResidueNames}{$ResidueCount} .= "~${ResidueName}";
	}
	else {
	  $AtomResiduesCountToNameMap{ResidueNames}{$ResidueCount} = $ResidueName;
	}
      }
      for $ResidueName (keys %{$ResiduesInfoRef->{HetatmResidueCount}}) {
	$ResidueCount = $ResiduesInfoRef->{HetatmResidueCount}{$ResidueName};
	if (exists $HetatmResiduesCountToNameMap{ResidueNames}{$ResidueCount}) {
	  $HetatmResiduesCountToNameMap{ResidueNames}{$ResidueCount} .= "~${ResidueName}";
	}
	else {
	  $HetatmResiduesCountToNameMap{ResidueNames}{$ResidueCount} = $ResidueName;
	}
      }
    }

    # Setup distribution of residues info...
    my($ResidueNames, $PercentResidueCount, @ResidueNamesList, %ResiduesCountInfoMap, %AtomResiduesCountInfoMap, %HetatmResiduesCountInfoMap);

    @{$ResiduesCountInfoMap{Frequency}} = ();
    @{$ResiduesCountInfoMap{PercentFrequency}} = ();
    for $ResidueCount (sort {$b <=> $a} keys %{$ResiduesCountToNameMap{ResidueNames}}) {
      $PercentResidueCount = sprintf("%.1f", (($ResidueCount/$TotalResidueCount)*100)) + 0;
      $ResidueNames = $ResiduesCountToNameMap{ResidueNames}{$ResidueCount};
      @ResidueNamesList = split /~/, $ResidueNames;
      for $ResidueName (@ResidueNamesList) {
	push @{$ResiduesCountInfoMap{Frequency}}, "${ResidueName} - ${ResidueCount}";
	push @{$ResiduesCountInfoMap{PercentFrequency}}, "${ResidueName} - ${PercentResidueCount}";
      }
    }
    if ($OptionsInfo{DetailLevel} >= 1) {
      @{$AtomResiduesCountInfoMap{Frequency}} = ();
      @{$AtomResiduesCountInfoMap{PercentFrequency}} = ();
      for $ResidueCount (sort {$b <=> $a} keys %{$AtomResiduesCountToNameMap{ResidueNames}}) {
	$PercentResidueCount = sprintf("%.1f", (($ResidueCount/$TotalResidueCount)*100)) + 0;
	$ResidueNames = $AtomResiduesCountToNameMap{ResidueNames}{$ResidueCount};
	@ResidueNamesList = split /~/, $ResidueNames;
	for $ResidueName (@ResidueNamesList) {
	  push @{$AtomResiduesCountInfoMap{Frequency}}, "${ResidueName} - ${ResidueCount}";
	  push @{$AtomResiduesCountInfoMap{PercentFrequency}}, "${ResidueName} - ${PercentResidueCount}";
	}
      }
      @{$HetatmResiduesCountInfoMap{Frequency}} = ();
      @{$HetatmResiduesCountInfoMap{PercentFrequency}} = ();
      for $ResidueCount (sort {$b <=> $a} keys %{$HetatmResiduesCountToNameMap{ResidueNames}}) {
	$PercentResidueCount = sprintf("%.1f", (($ResidueCount/$TotalResidueCount)*100)) + 0;
	$ResidueNames = $HetatmResiduesCountToNameMap{ResidueNames}{$ResidueCount};
	@ResidueNamesList = split /~/, $ResidueNames;
	for $ResidueName (@ResidueNamesList) {
	  push @{$HetatmResiduesCountInfoMap{Frequency}}, "${ResidueName} - ${ResidueCount}";
	  push @{$HetatmResiduesCountInfoMap{PercentFrequency}}, "${ResidueName} - ${PercentResidueCount}";
	}
      }
    }

    # List distribution of residues
    print "\nDistribution of residues: ", JoinWords(\@{$ResiduesCountInfoMap{Frequency}},'; ', 0), "\n";
    if ($OptionsInfo{DetailLevel} >= 2) {
      print "\nPercent distribution of residues: ", JoinWords(\@{$ResiduesCountInfoMap{PercentFrequency}},'; ', 0), "\n";
    }

    if ($OptionsInfo{DetailLevel} >= 1) {
      print "\nDistribution of ATOM residues: ", JoinWords(\@{$AtomResiduesCountInfoMap{Frequency}},'; ', 0), "\n";
      if ($OptionsInfo{DetailLevel} >= 2) {
	print "\nPercent distribution of ATOM residues: ", JoinWords(\@{$AtomResiduesCountInfoMap{PercentFrequency}},'; ', 0), "\n";
      }

      print "\nDistribution of HETATM residues: ", JoinWords(\@{$HetatmResiduesCountInfoMap{Frequency}},'; ', 0), "\n";
      if ($OptionsInfo{DetailLevel} >= 2) {
	print "\nPercent distribution of HETATM residues: ", JoinWords(\@{$HetatmResiduesCountInfoMap{PercentFrequency}},'; ', 0), "\n";
      }
    }
  }
}

# List information about residue numbers for each chain...
sub ListResidueNumbersInfo {
  my($PDBRecordLinesRef) = @_;
  my($Index, $ResidueCount, $StartResidueNum, $EndResidueNum, $ChainID, $CollectChainResiduesBeyondTER, $ChainsAndResiduesInfoRef, $ResidueNum, $PreviousResidueNum, $ResidueName, $PreviousResidueName, $GapResiduePairsCount, $GapLength, $DescendingOrderResiduePairsCount, @DescendingOrderResiduePairs, @GapResiduePairs);

  $CollectChainResiduesBeyondTER = 0;
  $ChainsAndResiduesInfoRef = GetChainsAndResidues($PDBRecordLinesRef, 'AtomAndHetatm', $CollectChainResiduesBeyondTER);

  print "\nATOM/HETATM residue numbers information for chains:\n";

  for $ChainID (@{$ChainsAndResiduesInfoRef->{ChainIDs}}) {
    print "\nChain ID - ",  CleanupChainID($ChainID), "";

    $ResidueCount = @{$ChainsAndResiduesInfoRef->{ResidueNumbers}{$ChainID}};

    # Start and end residue numbers...
    $StartResidueNum = $ChainsAndResiduesInfoRef->{ResidueNumbers}{$ChainID}[0];
    $EndResidueNum = $ChainsAndResiduesInfoRef->{ResidueNumbers}{$ChainID}[$ResidueCount - 1];
    print "; Number of residues: $ResidueCount; Start residue number - $StartResidueNum; End residue number - $EndResidueNum\n";

    # Identify any gaps in residue numbers or non-ascending order residue numbers...
    $GapResiduePairsCount = 0;
    $DescendingOrderResiduePairsCount = 0;

    @DescendingOrderResiduePairs = ();
    @GapResiduePairs = ();

    RESIDUE: for $Index (1 .. ($ResidueCount - 1)) {
      $ResidueNum = $ChainsAndResiduesInfoRef->{ResidueNumbers}{$ChainID}[$Index];
      $PreviousResidueNum = $ChainsAndResiduesInfoRef->{ResidueNumbers}{$ChainID}[$Index - 1];

      $ResidueName = $ChainsAndResiduesInfoRef->{Residues}{$ChainID}[$Index];
      $PreviousResidueName = $ChainsAndResiduesInfoRef->{Residues}{$ChainID}[$Index - 1];

      if ($ResidueNum == ($PreviousResidueNum + 1)) {
	# All is good...
	next RESIDUE;
      }

      # Are residue in descending order?
      if ($ResidueNum < $PreviousResidueNum) {
	$DescendingOrderResiduePairsCount++;
	push @DescendingOrderResiduePairs, "<${PreviousResidueName}${PreviousResidueNum} - ${ResidueName}${ResidueNum}>";
      }

      # Track gaps in residue pairs...
      $GapResiduePairsCount++;
      $GapLength = abs($ResidueNum - $PreviousResidueNum) - 1;

      push @GapResiduePairs, "<${PreviousResidueName}${PreviousResidueNum} - ${ResidueName}${ResidueNum}; $GapLength>";
    }

    # Print gaps information...
    print "Gaps in residue numbers: ", $GapResiduePairsCount ? "Yes" : "None";
    if ($GapResiduePairsCount) {
      print "; Number of gap residue number pairs: $GapResiduePairsCount; Gap residue pairs: <StartRes-EndRes; GapLength> - ", JoinWords(\@GapResiduePairs, "; ", 0);
    }
    print "\n";

    # Print descending residue order information...
    print "Residue numbers in descending order: ", $DescendingOrderResiduePairsCount ? "Yes" : "None";
    if ($DescendingOrderResiduePairsCount) {
      print "; Number of descending residue number pairs: $DescendingOrderResiduePairsCount; Descending residue number pairs: <StartRes-EndRes> ", JoinWords(\@DescendingOrderResiduePairs, "; ", 0);
    }
    print "\n";
  }
}

# List min/max XYZ coordinates for ATOM/HETATM records...
sub ListBoundingBox {
  my($PDBRecordLinesRef) = @_;
  my($XMin, $YMin, $ZMin, $XMax, $YMax, $ZMax, $XSize, $YSize, $ZSize);

  ($XMin, $YMin, $ZMin, $XMax, $YMax, $ZMax) = GetMinMaxCoords($PDBRecordLinesRef);
  $XSize = abs($XMax - $XMin);
  $YSize = abs($YMax - $YMin);
  $ZSize = abs($ZMax - $ZMin);

  $XMin = sprintf("%.3f", $XMin) + 0; $XMax = sprintf("%.3f", $XMax) + 0;
  $YMin = sprintf("%.3f", $YMin) + 0; $YMax = sprintf("%.3f", $YMax) + 0;
  $ZMin = sprintf("%.3f", $ZMin) + 0; $ZMax = sprintf("%.3f", $ZMax) + 0;

  $XSize = sprintf("%.3f", $XSize) + 0;
  $YSize = sprintf("%.3f", $YSize) + 0;
  $ZSize = sprintf("%.3f", $ZSize) + 0;

  print "\nBounding box coordinates: <XMin, XMax> - <$XMin, $XMax>; <YMin, YMax> - <$YMin, $YMax>; <ZMin, ZMax> - <$ZMin, $ZMax>;\n";
  print "Bounding box size in angstroms: XSize - $XSize; YSize - $YSize; ZSize - $ZSize\n";

}

# Check master record counts against actual record counts...
sub CheckMasterRecord {
  my($RecordTypesCountRef, $PDBRecordLinesRef) = @_;

  # Get master record information...
  my($NumOfRemarkRecords, $NumOfHetRecords, $NumOfHelixRecords, $NumOfSheetRecords, $NumOfTurnRecords, $NumOfSiteRecords, $NumOfTransformationsRecords, $NumOfAtomAndHetatmRecords, $NumOfTerRecords, $NumOfConectRecords, $NumOfSeqresRecords) = (undef) x 11;
  my($RecordLine, $MasterRecordFound);
  $MasterRecordFound = 0;

  LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
      if (IsMasterRecordType($RecordLine)) {
	($NumOfRemarkRecords, $NumOfHetRecords, $NumOfHelixRecords, $NumOfSheetRecords, $NumOfTurnRecords, $NumOfSiteRecords, $NumOfTransformationsRecords, $NumOfAtomAndHetatmRecords, $NumOfTerRecords, $NumOfConectRecords, $NumOfSeqresRecords) = ParseMasterRecordLine($RecordLine);
	$MasterRecordFound = 1;
	last LINE;
      }
  }
  if (!$MasterRecordFound) {
    print "\nWarning: MASTER record is missing.\n";
    return;
  }
  my(@MasterRecordValidationInfo);
  @MasterRecordValidationInfo = ();
  $NumOfRemarkRecords += 0;
  if (exists($RecordTypesCountRef->{Count}{REMARK}) && $NumOfRemarkRecords != $RecordTypesCountRef->{Count}{REMARK}) {
    push @MasterRecordValidationInfo, "Number of REMARK records, $NumOfRemarkRecords, specified in MASTER record doen't match its explict count, $RecordTypesCountRef->{Count}{REMARK}.";
  }
  $NumOfHetRecords += 0;
  if (exists($RecordTypesCountRef->{Count}{HET}) && $NumOfHetRecords != $RecordTypesCountRef->{Count}{HET}) {
    push @MasterRecordValidationInfo, "Number of HET records, $NumOfHetRecords, specified in MASTER record doen't match its explict count, $RecordTypesCountRef->{Count}{HET}.";
  }
  $NumOfHelixRecords += 0;
  if (exists($RecordTypesCountRef->{Count}{HELIX}) && $NumOfHelixRecords != $RecordTypesCountRef->{Count}{HELIX}) {
    push @MasterRecordValidationInfo, "Number of HELIX records, $NumOfHelixRecords, specified in MASTER record doen't match its explict count, $RecordTypesCountRef->{Count}{HELIX}.";
  }
  $NumOfSheetRecords += 0;
  if (exists($RecordTypesCountRef->{Count}{SHEET}) && $NumOfSheetRecords != $RecordTypesCountRef->{Count}{SHEET}) {
    push @MasterRecordValidationInfo, "Number of SHEET records, $NumOfSheetRecords, specified in MASTER record doen't match its explict count, $RecordTypesCountRef->{Count}{SHEET}.";
  }
  $NumOfTurnRecords += 0;
  if (exists($RecordTypesCountRef->{Count}{TURN}) && $NumOfTurnRecords != $RecordTypesCountRef->{Count}{TURN}) {
    push @MasterRecordValidationInfo, "Number of TURN records, $NumOfTurnRecords, specified in MASTER record doen't match its explict count, $RecordTypesCountRef->{Count}{REMARK}.";
  }
  $NumOfSiteRecords += 0;
  if (exists($RecordTypesCountRef->{Count}{SITE}) && $NumOfSiteRecords != $RecordTypesCountRef->{Count}{SITE}) {
    push @MasterRecordValidationInfo, "Number of SITE records, $NumOfSiteRecords, specified in MASTER record doen't match its explict count, $RecordTypesCountRef->{Count}{SITE}.";
  }

  $NumOfTransformationsRecords += 0;
  my($RecordsCount, $ID, $RecordID, $RecordLabel);
  $RecordsCount = 0;
  for $RecordLabel ('ORIGX', 'SCALE', 'MTRIX') {
    for $ID (1 .. 3) {
      $RecordID = "${RecordLabel}${ID}";
      if (exists $RecordTypesCountRef->{Count}{$RecordID}) {
	$RecordsCount += $RecordTypesCountRef->{Count}{$RecordID};
      }
    }
  }
  if ($NumOfTransformationsRecords != $RecordsCount) {
    push @MasterRecordValidationInfo, "Number of transformation records (ORIGXn+SCALEn+MTRIXn), $NumOfTransformationsRecords, specified in MASTER record doen't match its explict count, $RecordsCount.";
  }

  $RecordsCount = 0;
  for $RecordLabel ('ATOM', 'HETATM') {
      if (exists $RecordTypesCountRef->{Count}{$RecordLabel}) {
	$RecordsCount += $RecordTypesCountRef->{Count}{$RecordLabel};
      }
  }
  $NumOfAtomAndHetatmRecords += 0;
  if ($NumOfAtomAndHetatmRecords != $RecordsCount) {
    push @MasterRecordValidationInfo, "Number of ATOM + HETATM records, $NumOfAtomAndHetatmRecords, specified in MASTER record doen't match its explict count, $RecordsCount.";
  }
  $NumOfTerRecords += 0;
  if (exists($RecordTypesCountRef->{Count}{TER}) && $NumOfTerRecords != $RecordTypesCountRef->{Count}{TER}) {
    push @MasterRecordValidationInfo, "Number of TER records, $NumOfTerRecords, specified in MASTER record doen't match its explict count, $RecordTypesCountRef->{Count}{TER}.";
  }
  $NumOfConectRecords += 0;
  if (exists($RecordTypesCountRef->{Count}{CONECT}) && $NumOfConectRecords != $RecordTypesCountRef->{Count}{CONECT}) {
    push @MasterRecordValidationInfo, "Number of CONECT records, $NumOfConectRecords, specified in MASTER record doen't match its explict count, $RecordTypesCountRef->{Count}{CONECT}.";
  }
  $NumOfSeqresRecords += 0;
  if (exists($RecordTypesCountRef->{Count}{SEQRES}) && $NumOfSeqresRecords != $RecordTypesCountRef->{Count}{SEQRES}) {
    push @MasterRecordValidationInfo, "Number of SITE records, $NumOfSeqresRecords, specified in MASTER record doen't match its explict count, $RecordTypesCountRef->{Count}{SEQRES}.";
  }

  if (@MasterRecordValidationInfo) {
    print "\nMASTER record validation: Count mismatches found:\n";
    print JoinWords(\@MasterRecordValidationInfo, "\n", 0), "\n";
  }
  else {
    print "\nMASTER record validation: Count values match with the explicit count of the corresponding records.\n";
  }
}

# Total size of all the files...
sub ListTotalSizeOfFiles {
  my($FileOkayCount, $TotalSize, $Index);

  $FileOkayCount = 0;
  $TotalSize = 0;

  for $Index (0 .. $#PDBFilesList) {
    if ($PDBFilesInfo{FileOkay}[$Index]) {
      $FileOkayCount++;
      $TotalSize += $PDBFilesInfo{FileSize}[$Index];
    }
  }
  if ($FileOkayCount > 1) {
    print "\nTotal size of $FileOkayCount files: ", FormatFileSize($TotalSize), "\n";
  }

}

# Empty chain IDs are replaced with "None[1-9]". But for displaying purposes, take out any
# numbers from label...
sub CleanupChainID {
  my($ChainID) = @_;

  if ($ChainID =~ /^None/i) {
    return 'None';
  }
  return $ChainID;
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  # Setup record types to count...
  if ($Options{count}) {
    $OptionsInfo{CountRecordType} = $Options{count};
  }
  else {
    $OptionsInfo{CountRecordType} = $Options{all} ? 'All' : 'ATOM,HETATM';
  }
  @{$OptionsInfo{SpecifiedRecordTypes}} =();
  if ($OptionsInfo{CountRecordType} !~ /^All$/i) {
    my(@RecordTypes);
    @RecordTypes = split /\,/, $OptionsInfo{CountRecordType};
    push @{$OptionsInfo{SpecifiedRecordTypes}}, @RecordTypes;
  }
  $OptionsInfo{CountChains} = ($Options{chains} || $Options{all}) ? 1 : 0;
  $OptionsInfo{CheckMasterRecord} = ($Options{mastercheck} || $Options{all}) ? 1 : 0;

  # Residue count is the default. So $Options{residues} is simply ignored.
  my($CountResidues) = 1;
  $OptionsInfo{CountResiduesInChains} = (($CountResidues || $Options{all}) && $Options{residuesmode} =~ /^(InChains|Both)$/i) ? 1 : 0;
  $OptionsInfo{CountResiduesAll} = (($CountResidues || $Options{all}) && $Options{residuesmode} =~ /^(All|Both)$/i) ? 1 : 0;

  $OptionsInfo{ResiduesFrequencyInChains} = (($Options{frequency} || $Options{all}) && $Options{residuesmode} =~ /^(InChains|Both)$/i) ? 1 : 0;
  $OptionsInfo{ResiduesFrequencyAll} = (($Options{frequency} || $Options{all}) && $Options{residuesmode} =~ /^(All|Both)$/i) ? 1 : 0;

  $OptionsInfo{ResidueNumbersInfo} = ($Options{residuenumbers} || $Options{all})  ? 1 : 0;

  $OptionsInfo{CalculateBoundingBox} = ($Options{boundingbox} || $Options{all}) ? 1 : 0;

  $OptionsInfo{ListHeaderInfo} = ($Options{header} || $Options{all}) ? 1 : 0;
  $OptionsInfo{DetailLevel} = $Options{detail};

  $OptionsInfo{ListExperimentalTechniqueInfo} = ($Options{experiment} || $Options{all}) ? 1 : 0;

}

# Retrieve information about PDB files...
sub RetrievePDBFilesInfo {
  my($Index, $PDBFile, $ModifiedTimeString, $ModifiedDateString);

  %PDBFilesInfo = ();
  @{$PDBFilesInfo{FileOkay}} = ();
  @{$PDBFilesInfo{FileSize}} = ();
  @{$PDBFilesInfo{FileLastModified}} = ();

  FILELIST: for $Index (0 .. $#PDBFilesList) {
    $PDBFilesInfo{FileOkay}[$Index] = 0;
    $PDBFilesInfo{FileSize}[$Index] = 0;
    $PDBFilesInfo{FileLastModified}[$Index] = '';

    $PDBFile = $PDBFilesList[$Index];
    if (!(-e $PDBFile)) {
      warn "Warning: Ignoring file $PDBFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($PDBFile, "pdb")) {
      warn "Warning: Ignoring file $PDBFile: It's not a PDB file\n";
      next FILELIST;
    }
    if (! open PDBFILE, "$PDBFile") {
      warn "Warning: Ignoring file $PDBFile: Couldn't open it: $! \n";
      next FILELIST;
    }
    close PDBFILE;

    $PDBFilesInfo{FileOkay}[$Index] = 1;
    $PDBFilesInfo{FileSize}[$Index] = FileSize($PDBFile);
    ($ModifiedTimeString, $ModifiedDateString) = FormattedFileModificationTimeAndDate($PDBFile);
    $PDBFilesInfo{FileLastModified}[$Index] = "$ModifiedTimeString; $ModifiedDateString";
  }
}


# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{count} = '';
  $Options{detail} = 1;
  $Options{residuesmode} = 'Both';

  if (!GetOptions(\%Options, "all|a", "boundingbox|b", "count|c=s", "chains", "detail|d=i", "experiment|e", "frequency|f", "mastercheck|m", "header", "help|h", "residues", "residuesmode=s", "residuenumbers", "workingdir|w=s")) {
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
  if ($Options{residuesmode} !~ /^(InChains|All|Both)$/i) {
    die "Error: The value specified, $Options{residuesmode}, for option \"--ResiduesMode\" is not valid. Allowed values: InChains, All, or Both\n";
  }
}

__END__

=head1 NAME

InfoPDBFiles.pl - List information about PDB files

=head1 SYNOPSIS

InfoPDBFiles.pl PDBFile(s) PDB(s)...

InfoPDBFiles.pl [B<-a, --all>] [B<-b, --BoundingBox>]
[B<-c, --count> "RecordType, [RecordType,...]" | All] [B<--chains>]
[B<-d, --detail> infolevel] [B<-e, --experiment>] [B<-f, --frequency>]
[B<-h, --help>] [B<--header>] [B<m, --MasterCheck>] [B<--residues>]
[B<--ResiduesMode> InChains | All | Both] [B<--ResidueNumbers>]
[B<-w, --WorkingDir> dirname] PDBFile(s)...

=head1 DESCRIPTION

List information about contents of I<PDBFile(s)>: number of each record type, number of chains,
count and percent distribution of residues in each chain, bounding box and so on.
Multiple PDBFile names are separated by spaces. The valid file extension is I<.pdb>.
All other file name extensions are ignored during the wild card expansion. All the PDB files
in a current directory can be specified either by I<*.pdb> or the current directory name.

In PDB files containing data for multiple models, all ATOM/HETAM records for chains after the first model
are ignored.

=head1 OPTIONS

=over 4

=item B<-a, --all>

List all the available information.

=item B<-b, --BoundingBox>

List min/max XYZ coordiates of ATOM/HETATM records.

=item B<-c, --count> I<RecordType,[RecordType,...]|All>

Types of PDB records to count in I<PDBFile(s)>. You can specify a list of any valid PDB
record type or count all record types found in the files. Possible values: Comma delimited list
of valid I<RecordTypes> or I<All>. Default: I<ATOM,HETATM>. And this is also B<default behavior>.

The list of valid PDB record types includes: I<HEADER, OBSLTE, TITLE, CAVEAT, COMPND, SOURCE, KEYWDS,
EXPDTA, AUTHOR, REVDAT, SPRSDE, JRN, REMARK, DBRE, SEQADV, SEQRES, MODRES, HET, HETNAM, HETSYN,
FORMUL, HELIX, SHEET, TURN, SSBOND, LINK, HYDBND, SLTBRG, CISPEP, SITE, CRYST1, ORIGX1, ORIGX2, ORIGX3,
SCALE1, SCALE2, SCALE3, MTRIX1 MTRIX2 MTRIX3, TVECT, MODEL, ATOM, SIGATM, ANISOU, SIGUIJ, TER,
HETATM, ENDMDL, CONECT, MASTER, END>.

=item B<--chains>

Count number of chains.

=item B<-d, --detail> I<infolevel>

Level of information to print about PDB during various options. Default: I<1>.
Possible values: I<1, 2 or 3>.

=item B<-e, --experiment>

List experimental technique information along with any applicable resolution.

=item B<-f, --frequency>

List distribution of residues: report count and percent of residues in individual chains and
across all the chains, or for all the residues in the file. The value of option B<--residuesmode>
determines how residues are counted and what is listed. The list is sorted by frequency in
descending order. By default, only residue count values are reported. To list percent distribution
of residues, specify B<-d, --detail> value of I<2> or higher.

=item B<-h, --help>

Print this help message.

=item B<--header>

List header information.

=item B<m, --MasterCheck>

Check master record by explicitly counting the number of REMARK, HET, HELIX, SHEET, TURN, SITE,
ORIGX, SCALE, MTRIX, ATOM, HETATM, TER, CONECT and SEQRES records and comparing their
values against contents of master record.

=item B<--residues>

Count residues in I<PDBFile(s)>. This is also B<default behavior>.

By default, only residue count values are reported. To list percent distribution of residues,
specify B<-d, --detail> value of I<2> or higher.

=item B<--ResiduesMode> <InChains | All | Both>

Specify how to count residues in I<PDBFile(s)>: Count residue in each chain and across all the chains,
list count iof all the residues in the file, or list both. Possible values: I<InChains, All, or Both>.
Default: I<Both>.

=item B<--ResidueNumbers>

List information about ATOM residue numbers in each chain before TER record: start and end residue
number; gaps in residue numbers corresponding to non-sequential residue numbers; residue
numbers not in ascending order.

=item B<-w, --WorkingDir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To list total number of records and number of chain(s) residues in PDB files, type:

    % InfoPDBFiles.pl Sample1.pdb
    % InfoPDBFiles.pl Sample2.pdb

To list all available information for PDB file Sample2.pdb, type:

    % InfoPDBFiles.pl -a Sample2.pdb

To list all available information for PDB file Sample2.pdb with all available details, type:

    % InfoPDBFiles.pl -a -d Sample2.pdb

To count ATOM and HETATM records in Sample2.pdb file, type:

    % InfoPDBFiles.pl -c "ATOM,HETATM" Sample2.pdb

To list distribution of residues in chains across the whole PDB file Sample2.pdb along with
percent distribution, type

    % InfoPDBFiles.pl --frequency -d 2 Sample2.pdb

To list distribution of residues only across chains in PDB file Sample2.pdb along with
percent distribution, type

    % InfoPDBFiles.pl --frequency -d 2 --ResiduesMode InChains Sample2.pdb

To list min/max coordinates of the bounding box which encompasses the structure in Sample1.pdb
file, type:

    % InfoPDBFiles.pl -b Sample1.pdb

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

ExtractFromPDBFiles.pl, InfoAminoAcids.pl, InfoNucleicAcids.pl, InfoSequenceFiles.pl, ModifyPDBFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
