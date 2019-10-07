#!/usr/bin/perl -w
#
# File: ExtractFromPDBFiles.pl
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
use AminoAcids;
use SequenceFileUtil;

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
print "Checking input PDB file(s)...\n";
my(%PDBFilesInfo);
RetrievePDBFilesInfo();

# Process input files..
my($FileIndex);
if (@PDBFilesList > 1) {
  print "\nProcessing PDB files...\n";
}
for $FileIndex (0 .. $#PDBFilesList) {
  if ($PDBFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $PDBFilesList[$FileIndex]...\n";
    ExtractFromPDBFiles($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Extract appropriate information...
sub ExtractFromPDBFiles {
  my($FileIndex) = @_;
  my($PDBFile, $PDBRecordLinesRef);

  # Get PDB data...
  $PDBFile = $PDBFilesList[$FileIndex];
  $PDBRecordLinesRef = ReadPDBFile($PDBFile);

  if ($OptionsInfo{Mode} =~ /Chains/i) {
    ExtractChains($FileIndex, $PDBRecordLinesRef);
  }
  elsif ($OptionsInfo{Mode} =~ /Sequences/i) {
    ExtractSequences($FileIndex, $PDBRecordLinesRef);
  }
  elsif ($OptionsInfo{Mode} =~ /^(Atoms|CAlphas|AtomNums|AtomsRange|AtomNames)$/i) {
    ExtractByAtoms($FileIndex, $PDBRecordLinesRef);
  }
  elsif ($OptionsInfo{Mode} =~ /^(ResidueNums|ResiduesRange|ResidueNames)$/i) {
    ExtractByResidues($FileIndex, $PDBRecordLinesRef);
  }
  elsif ($OptionsInfo{Mode} =~ /Distance/i) {
    ExtractByDistance($FileIndex, $PDBRecordLinesRef);
  }
  elsif ($OptionsInfo{Mode} =~ /NonWater/i) {
    ExtractNonWaterRecords($FileIndex, $PDBRecordLinesRef);
  }
  elsif ($OptionsInfo{Mode} =~ /NonHydrogens/i) {
    ExtractNonHydrogenRecords($FileIndex, $PDBRecordLinesRef);
  }
}

# Extract chains and generate new PDB files...
#
sub ExtractChains {
  my($FileIndex, $PDBRecordLinesRef) = @_;
  my($ChainIndex, $ChainID, $ChainLabel, $PDBFileName, $RecordLine, $ChainsAndResiduesInfoRef, $AtomNumber, $AtomName, $ResidueName, $AtomChainID, $ResidueNumber, $AlternateLocation, $InsertionCode, $ConectRecordLinesRef, $CollectChainResiduesBeyondTER, %ChainAtomNumbersMap);

  $CollectChainResiduesBeyondTER = ($OptionsInfo{ChainsRecordMode} =~ /^AcrossTER$/i) ? 1 : 0;

  # Get chains and residues data...
  $ChainsAndResiduesInfoRef = GetChainsAndResidues($PDBRecordLinesRef, 'AtomAndHetatm', $CollectChainResiduesBeyondTER, 1);

  if ($OptionsInfo{CombineChains}) {
    $PDBFileName = $PDBFilesInfo{OutFileNames}[$FileIndex][0];
    print "Generating PDBFileName file $PDBFileName...\n";

    open OUTFILE, ">$PDBFileName" or die "Error: Can't open $PDBFileName: $! \n";

    # Write out header and other older records...
    WriteHeaderAndOlderRecords(\*OUTFILE, $PDBRecordLinesRef);
  }

  for $ChainIndex (0 .. $#{$PDBFilesInfo{SpecifiedChains}[$FileIndex]}) {
    $ChainID = $PDBFilesInfo{SpecifiedChains}[$FileIndex][$ChainIndex];
    $ChainLabel = $PDBFilesInfo{ChainLabels}[$FileIndex][$ChainIndex];

    if (!$OptionsInfo{CombineChains}) {
      $PDBFileName = $PDBFilesInfo{OutFileNames}[$FileIndex][$ChainIndex];
      print "Generating PDBFileName file $PDBFileName...\n";

      open OUTFILE, ">$PDBFileName" or die "Error: Can't open $PDBFileName: $! \n";

      # Write out header and other older recors...
      WriteHeaderAndOlderRecords(\*OUTFILE, $PDBRecordLinesRef);
    }

    # Write out ATOM/HETATM line for chain and collect all ATOM/HETATM serial numbers
    # for writing out appropriate CONECT records...
    %ChainAtomNumbersMap = ();
    for $RecordLine (@{$ChainsAndResiduesInfoRef->{Lines}{$ChainID}}) {
      if (CheckRecordType($RecordLine)) {
	print OUTFILE "$RecordLine\n";
	($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $AtomChainID, $ResidueNumber, $InsertionCode) = ParseAtomRecordLine($RecordLine);
	$AtomNumber = int $AtomNumber;
	$ChainAtomNumbersMap{$AtomNumber} = $AtomName;
      }
    }
    # Write out TER record using information from last chain record...
    $AtomNumber += 1;
    print OUTFILE GenerateTerRecordLine($AtomNumber, $ResidueName, $AtomChainID, $ResidueNumber, $InsertionCode), "\n";

    # Write out CONECT records...
    $ConectRecordLinesRef = GetConectRecordLines($PDBRecordLinesRef, \%ChainAtomNumbersMap);

    for $RecordLine (@{$ConectRecordLinesRef}) {
      print OUTFILE "$RecordLine\n";
    }

    if (!$OptionsInfo{CombineChains}) {
      # Write out END record...
      print OUTFILE GenerateEndRecordLine(), "\n";

      close OUTFILE;
    }
  }

  if ($OptionsInfo{CombineChains}) {
    # Write out END record...
    print OUTFILE GenerateEndRecordLine(), "\n";

    close OUTFILE;
  }

}

# Extract sequences for individual chains or combine all the chains...
sub ExtractSequences {
  my($FileIndex, $PDBRecordLinesRef) = @_;
  my($ChainIndex, $ChainID, $ChainLabel, $SequenceFileName, $Residue, $ResidueCode, $StandardResidue, $ChainSequence, $WrappedChainSequence, $ChainSequenceID, $ChainsAndResiduesInfoRef, $ChainResiduesRef, %ChainSequencesDataMap);

  if ($OptionsInfo{SequenceRecordSource} =~ /^SeqRes$/i) {
    $ChainsAndResiduesInfoRef = GetChainsAndResidues($PDBRecordLinesRef, 'SeqRes');
  }
  else {
    $ChainsAndResiduesInfoRef = GetChainsAndResidues($PDBRecordLinesRef);
  }

  # Generate sequence data for all the chains...
  %ChainSequencesDataMap = ();
  @{$ChainSequencesDataMap{IDs}} = ();
  %{$ChainSequencesDataMap{Sequence}} = ();
  %{$ChainSequencesDataMap{Description}} = ();

  for $ChainIndex (0 .. $#{$PDBFilesInfo{SpecifiedChains}[$FileIndex]}) {
    $ChainID = $PDBFilesInfo{SpecifiedChains}[$FileIndex][$ChainIndex];
    $ChainLabel = $PDBFilesInfo{ChainLabels}[$FileIndex][$ChainIndex];

    # Setup sequence ID...
    $ChainSequenceID = $PDBFilesInfo{ChainSequenceIDs}[$FileIndex][$ChainIndex];
    push @{$ChainSequencesDataMap{IDs}}, $ChainSequenceID;
    $ChainSequencesDataMap{Description}{$ChainID} = $ChainSequenceID;

    # Collect sequence data for the chain...
    if ($OptionsInfo{SequenceRecordSource} =~ /^SeqRes/i) {
      $ChainResiduesRef = \@{$ChainsAndResiduesInfoRef->{Residues}{$ChainID}};
    }
    else {
      $ChainResiduesRef = \@{$ChainsAndResiduesInfoRef->{Residues}{$ChainID}};
    }
    # Setup sequence data...
    $ChainSequence = '';
    RESIDUE: for $Residue (@{$ChainResiduesRef}) {
      ($ResidueCode, $StandardResidue) = GetResidueCode($Residue);
      if (!$StandardResidue) {
	if ($OptionsInfo{KeepNonStandardSequences}) {
	  $ResidueCode = $OptionsInfo{NonStandardSequenceCode};
	  warn "Warning: Keeping nonstandard residue $Residue in $ChainLabel...\n";
	}
	else {
	  warn "Warning: Ignoring nonstandard residue $Residue in $ChainLabel...\n";
	  next RESIDUE;
	}
      }
      $ChainSequence .= $ResidueCode;
    }
    $ChainSequencesDataMap{Sequence}{$ChainID} = $ChainSequence;

  }

  # Write out the sequence files...
  my($SequenceID, $SequenceDescription, $Sequence, %SequencesDataMap );
  if ($OptionsInfo{CombineChainSequences}) {
    # Combine all the chain sequences...
    $Sequence = '';
    for $ChainIndex (0 .. $#{$PDBFilesInfo{SpecifiedChains}[$FileIndex]}) {
      $ChainID = $PDBFilesInfo{SpecifiedChains}[$FileIndex][$ChainIndex];

      $Sequence .= $ChainSequencesDataMap{Sequence}{$ChainID};
    }
    $SequenceID = $PDBFilesInfo{ChainSequenceIDsPrefix}[$FileIndex][0] . "_CombinedChains|PDB";;
    $SequenceDescription = $SequenceID;
    $SequenceFileName = $PDBFilesInfo{OutFileNames}[$FileIndex][0];

    print "Generating sequence file $SequenceFileName...\n";
    %SequencesDataMap = ();
    @{$SequencesDataMap{IDs}} = ();
    %{$SequencesDataMap{Sequence}} = ();
    %{$SequencesDataMap{Description}} = ();

    push @{$SequencesDataMap{IDs}}, $SequenceID;
    $SequencesDataMap{Description}{$SequenceID} = $SequenceDescription;
    $SequencesDataMap{Sequence}{$SequenceID} = $Sequence;

    WritePearsonFastaSequenceFile($SequenceFileName, \%SequencesDataMap, $OptionsInfo{MaxSequenceLength});
  }
  else {
    # For each specifed chain, write out the sequences...
    for $ChainIndex (0 .. $#{$PDBFilesInfo{SpecifiedChains}[$FileIndex]}) {
      $ChainID = $PDBFilesInfo{SpecifiedChains}[$FileIndex][$ChainIndex];

      $SequenceFileName = $PDBFilesInfo{OutFileNames}[$FileIndex][$ChainIndex];

      $SequenceID = $ChainSequencesDataMap{IDs}[$ChainIndex];
      $SequenceDescription = $ChainSequencesDataMap{Description}{$ChainID};
      $Sequence = $ChainSequencesDataMap{Sequence}{$ChainID};

      print "Generating sequence file $SequenceFileName...\n";
      %SequencesDataMap = ();
      @{$SequencesDataMap{IDs}} = ();
      %{$SequencesDataMap{Sequence}} = ();
      %{$SequencesDataMap{Description}} = ();

      push @{$SequencesDataMap{IDs}}, $SequenceID;
      $SequencesDataMap{Description}{$SequenceID} = $SequenceDescription;
      $SequencesDataMap{Sequence}{$SequenceID} = $Sequence;

      WritePearsonFastaSequenceFile($SequenceFileName, \%SequencesDataMap, $OptionsInfo{MaxSequenceLength});
    }
  }
}

# Extract atoms...
sub ExtractByAtoms {
  my($FileIndex, $PDBRecordLinesRef) = @_;
  my($PDBFileName,  $RecordLine, $ChainRecordCount, $AtomNumber, $AtomName, $IgnoreRecord, $ConectRecordLinesRef, %AtomNumbersMap);

  $PDBFileName = $PDBFilesInfo{OutFileNames}[$FileIndex][0];
  print "Generating PDBFileName file $PDBFileName...\n";
  open OUTFILE, ">$PDBFileName" or die "Error: Can't open $PDBFileName: $! \n";

  # Write out header and other older recors...
  WriteHeaderAndOlderRecords(\*OUTFILE, $PDBRecordLinesRef);

  # Write out all ATOM records along with TER and model records to indicate
  # chains and multiple models..
  %AtomNumbersMap = ();
  $ChainRecordCount = 0;
  for $RecordLine (@{$PDBRecordLinesRef}) {
    if (CheckRecordType($RecordLine)) {
      ($AtomNumber, $AtomName) = ParseAtomRecordLine($RecordLine);

      # Check atoms...
      $IgnoreRecord = 1;
      if ($OptionsInfo{Mode} =~ /^Atoms$/i) {
	$IgnoreRecord = 0;
      }
      elsif ($OptionsInfo{Mode} =~ /^(CAlphas|AtomNames)$/i) {
	if (exists $OptionsInfo{SpecifiedAtomNamesMap}{lc $AtomName}) {
	  $IgnoreRecord = 0;
	}
      }
      elsif ($OptionsInfo{Mode} =~ /^AtomNums$/i) {
	if (exists $OptionsInfo{SpecifiedAtomNumsMap}{$AtomNumber}) {
	  $IgnoreRecord = 0;
	}
      }
      elsif ($OptionsInfo{Mode} =~ /^AtomsRange$/i) {
	if ($AtomNumber >= $OptionsInfo{SpecifiedStartAtomNum} && $AtomNumber <= $OptionsInfo{SpecifiedEndAtomNum}) {
	  $IgnoreRecord = 0;
	}
      }

      if (!$IgnoreRecord) {
	$ChainRecordCount++;
	print OUTFILE "$RecordLine\n";

	$AtomNumber = int $AtomNumber;
	$AtomNumbersMap{$AtomNumber} = $AtomName;
      }
    }
    elsif (IsTerRecordType($RecordLine)) {
      if ($ChainRecordCount) {
	print OUTFILE GenerateTerRecordLine(), "\n";
      }
      $ChainRecordCount = 0;
    }
    elsif (IsModelRecordType($RecordLine) || IsEndmdlRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
  }

  # Write out appropriate CONECT records...
  $ConectRecordLinesRef = GetConectRecordLines($PDBRecordLinesRef, \%AtomNumbersMap);
  for $RecordLine (@{$ConectRecordLinesRef}) {
    print OUTFILE "$RecordLine\n";
  }

  # Write out END record...
  print OUTFILE GenerateEndRecordLine(), "\n";

  close OUTFILE;
}

# Extract residues...
sub ExtractByResidues {
  my($FileIndex, $PDBRecordLinesRef) = @_;
  my($PDBFileName,  $RecordLine, $ChainRecordCount, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $ConectRecordLinesRef, $IgnoreRecord, %AtomNumbersMap);

  $PDBFileName = $PDBFilesInfo{OutFileNames}[$FileIndex][0];
  print "Generating PDBFileName file $PDBFileName...\n";
  open OUTFILE, ">$PDBFileName" or die "Error: Can't open $PDBFileName: $! \n";

  # Write out header and other older recors...
  WriteHeaderAndOlderRecords(\*OUTFILE, $PDBRecordLinesRef);

  # Write out all ATOM records for specified residues with TER and model records to indicate
  # chains and multiple models...
  %AtomNumbersMap = ();
  $ChainRecordCount = 0;
  for $RecordLine (@{$PDBRecordLinesRef}) {
    if (CheckRecordType($RecordLine)) {
      ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber) = ParseAtomRecordLine($RecordLine);

      # Check residues...
      $IgnoreRecord = 1;
      if ($OptionsInfo{Mode} =~ /^ResidueNums$/i) {
	if (exists $OptionsInfo{SpecifiedResidueNumsMap}{$ResidueNumber}) {
	  $IgnoreRecord = 0;
	}
      }
      elsif ($OptionsInfo{Mode} =~ /^ResiduesRange$/i) {
	if ($ResidueNumber >= $OptionsInfo{SpecifiedStartResidueNum} && $ResidueNumber <= $OptionsInfo{SpecifiedEndResidueNum}) {
	  $IgnoreRecord = 0;
	}
      }
      elsif ($OptionsInfo{Mode} =~ /^ResidueNames$/i) {
	if (exists $OptionsInfo{SpecifiedResidueNamesMap}{lc $ResidueName}) {
	  $IgnoreRecord = 0;
	}
      }
      if (!$IgnoreRecord) {
	$ChainRecordCount++;
	print OUTFILE "$RecordLine\n";
	$AtomNumber = int $AtomNumber;
	$AtomNumbersMap{$AtomNumber} = $AtomName;
      }
    }
    elsif (IsTerRecordType($RecordLine)) {
      if ($ChainRecordCount) {
	print OUTFILE GenerateTerRecordLine(), "\n";
      }
      $ChainRecordCount = 0;
    }
    elsif (IsModelRecordType($RecordLine) || IsEndmdlRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
  }

  # Write out appropriate CONECT records...
  $ConectRecordLinesRef = GetConectRecordLines($PDBRecordLinesRef, \%AtomNumbersMap);
  for $RecordLine (@{$ConectRecordLinesRef}) {
    print OUTFILE "$RecordLine\n";
  }
  # Write out END record...
  print OUTFILE GenerateEndRecordLine(), "\n";

  close OUTFILE;
}

# Extract non water records...
sub ExtractNonWaterRecords {
  my($FileIndex, $PDBRecordLinesRef) = @_;
  my($PDBFileName,  $RecordLine, $ChainRecordCount, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ConectRecordLinesRef, %AtomNumbersMap);

  $PDBFileName = $PDBFilesInfo{OutFileNames}[$FileIndex][0];
  print "Generating PDBFileName file $PDBFileName...\n";
  open OUTFILE, ">$PDBFileName" or die "Error: Can't open $PDBFileName: $! \n";

  # Write out header and other older recors...
  WriteHeaderAndOlderRecords(\*OUTFILE, $PDBRecordLinesRef);

  # Write out all ATOM/HETATM non water records along with TER and model records to indicate
  # chains and multiple models..
  %AtomNumbersMap = ();
  $ChainRecordCount = 0;
  for $RecordLine (@{$PDBRecordLinesRef}) {
    if (CheckRecordType($RecordLine)) {
      ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName) = ParseAtomRecordLine($RecordLine);
      if (! exists $OptionsInfo{SpecifiedWaterResiduesMap}{$ResidueName} ) {
	$ChainRecordCount++;
	print OUTFILE "$RecordLine\n";
	$AtomNumber = int $AtomNumber;
	$AtomNumbersMap{$AtomNumber} = $AtomName;
      }
    }
    elsif (IsTerRecordType($RecordLine)) {
      if ($ChainRecordCount) {
	print OUTFILE GenerateTerRecordLine(), "\n";
      }
      $ChainRecordCount = 0;
    }
    elsif (IsModelRecordType($RecordLine) || IsEndmdlRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
  }

  # Write out appropriate CONECT records...
  $ConectRecordLinesRef = GetConectRecordLines($PDBRecordLinesRef, \%AtomNumbersMap);
  for $RecordLine (@{$ConectRecordLinesRef}) {
    print OUTFILE "$RecordLine\n";
  }
  # Write out END record...
  print OUTFILE GenerateEndRecordLine(), "\n";

  close OUTFILE;
}

# Extract non hydrogen records...
sub ExtractNonHydrogenRecords {
  my($FileIndex, $PDBRecordLinesRef) = @_;
  my($PDBFileName,  $RecordLine, $ChainRecordCount, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge, $ConectRecordLinesRef, %AtomNumbersMap);

  $PDBFileName = $PDBFilesInfo{OutFileNames}[$FileIndex][0];
  print "Generating PDBFileName file $PDBFileName...\n";
  open OUTFILE, ">$PDBFileName" or die "Error: Can't open $PDBFileName: $! \n";

  # Write out header and other older recors...
  WriteHeaderAndOlderRecords(\*OUTFILE, $PDBRecordLinesRef);

  # Write out all ATOM/HETATM non hydrogen records along with TER and model records to indicate
  # chains and multiple models..
  %AtomNumbersMap = ();
  $ChainRecordCount = 0;
  for $RecordLine (@{$PDBRecordLinesRef}) {
    if (CheckRecordType($RecordLine)) {
      ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = ParseAtomRecordLine($RecordLine);
      if ($ElementSymbol !~ /^H$/i) {
	$ChainRecordCount++;
	print OUTFILE "$RecordLine\n";
	$AtomNumber = int $AtomNumber;
	$AtomNumbersMap{$AtomNumber} = $AtomName;
      }
    }
    elsif (IsTerRecordType($RecordLine)) {
      if ($ChainRecordCount) {
	print OUTFILE GenerateTerRecordLine(), "\n";
      }
      $ChainRecordCount = 0;
    }
    elsif (IsModelRecordType($RecordLine) || IsEndmdlRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
  }

  # Write out appropriate CONECT records...
  $ConectRecordLinesRef = GetConectRecordLines($PDBRecordLinesRef, \%AtomNumbersMap);
  for $RecordLine (@{$ConectRecordLinesRef}) {
    print OUTFILE "$RecordLine\n";
  }
  # Write out END record...
  print OUTFILE GenerateEndRecordLine(), "\n";

  close OUTFILE;
}

# Extract ATOM/HETATM records by distance...
sub ExtractByDistance {
  my($FileIndex, $PDBRecordLinesRef) = @_;
  my($PDBFileName,  $RecordLine, $RecordLineNum, $ChainRecordCount, $ConectRecordLinesRef, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $IgnoreRecord, $ResidueID, @OriginCoords, @Coords, %AtomNumbersMap, %ResiduesDataMap);

  $PDBFileName = $PDBFilesInfo{OutFileNames}[$FileIndex][0];
  print "Generating PDBFileName file $PDBFileName...\n";
  open OUTFILE, ">$PDBFileName" or die "Error: Can't open $PDBFileName: $! \n";

  # Write out header and other older recors...
  WriteHeaderAndOlderRecords(\*OUTFILE, $PDBRecordLinesRef);

  # Setup coordinates of origin to calculate distance...
  @OriginCoords = ();
  push @OriginCoords, @{$PDBFilesInfo{DistanceOrigin}[$FileIndex]};

  # Write out all ATOM records for which meet specified criteria along with TER and model records to indicate
  # chains and multiple models...
  %AtomNumbersMap = ();

  %ResiduesDataMap = ();
  %{$ResiduesDataMap{ID}} = ();
  %{$ResiduesDataMap{Status}} = ();

  $ChainRecordCount = 0;
  $RecordLineNum = 0;

  for $RecordLine (@{$PDBRecordLinesRef}) {
    $RecordLineNum++;
    if (CheckRecordType($RecordLine)) {
      ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z) = ParseAtomRecordLine($RecordLine);
      @Coords = (); push @Coords, ($X, $Y, $Z);

      $IgnoreRecord = 1;
      if ($OptionsInfo{DistanceSelectionMode} =~ /^ByResidue$/i) {
	$ResidueID = "${ResidueName}_${ResidueNumber}_${ChainID}";
	if (exists $ResiduesDataMap{ID}{$ResidueID}) {
	  # Residue data has been processed; check its selection status...
	  if ($ResiduesDataMap{Status}{$ResidueID}) {
	    $IgnoreRecord = 0;
	  }
	}
	else {
	  # Residue hasn't been processed...
	  $ResiduesDataMap{ID}{$ResidueID} = $ResidueID;
	  $ResiduesDataMap{Status}{$ResidueID} = 0;
	  if (CheckResidueDistance($ResidueID, $RecordLineNum, $PDBRecordLinesRef, \@OriginCoords)) {
	    $IgnoreRecord = 0;
	    $ResiduesDataMap{Status}{$ResidueID} = 1;
	  }
	}
      }
      elsif ($OptionsInfo{DistanceSelectionMode} =~ /^ByAtom$/i) {
	if (CheckDistance(\@Coords, \@OriginCoords)) {
	  $IgnoreRecord = 0;
	}
      }

      if (!$IgnoreRecord) {
	$ChainRecordCount++;
	print OUTFILE "$RecordLine\n";
	$AtomNumber = int $AtomNumber;
	$AtomNumbersMap{$AtomNumber} = $AtomName;
      }
    }
    elsif (IsTerRecordType($RecordLine)) {
      if ($ChainRecordCount) {
	print OUTFILE GenerateTerRecordLine(), "\n";
      }
      $ChainRecordCount = 0;
    }
    elsif (IsModelRecordType($RecordLine) || IsEndmdlRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
  }

  # Write out appropriate CONECT records...
  $ConectRecordLinesRef = GetConectRecordLines($PDBRecordLinesRef, \%AtomNumbersMap);
  for $RecordLine (@{$ConectRecordLinesRef}) {
    print OUTFILE "$RecordLine\n";
  }

  # Write out END record...
  print OUTFILE GenerateEndRecordLine(), "\n";

  close OUTFILE;
}

# Does record type correspond to the specified record type?
sub CheckRecordType {
  my($RecordLine) = @_;
  my($Status);

  $Status = 0;
  if ($OptionsInfo{RecordMode} =~ /^Atom$/i) {
    $Status = IsAtomRecordType($RecordLine) ? 1 : 0;
  }
  elsif ($OptionsInfo{RecordMode} =~ /^Hetatm$/i) {
    $Status = IsHetatmRecordType($RecordLine) ? 1 : 0;
  }
  elsif ($OptionsInfo{RecordMode} =~ /^AtomAndHetatm$/i) {
    $Status = (IsAtomRecordType($RecordLine) || IsHetatmRecordType($RecordLine)) ? 1 : 0;
  }

  return $Status;
}

# Does record meets distance citerion specified by the user?
sub CheckResidueDistance {
  my($SpecifiedResidueID, $StartingLineNum, $PDBRecordLinesRef, $OriginCoordsRef) = @_;
  my($Status, $RecordLine, $RecordLineIndex, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $ResidueID, @Coords);

  $Status = 0;

  RECORDLINE: for $RecordLineIndex (($StartingLineNum - 1) .. $#{$PDBRecordLinesRef}) {
    $RecordLine = $PDBRecordLinesRef->[$RecordLineIndex];
    if (!CheckRecordType($RecordLine)) {
      next RECORDLINE;
    }
    ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z) = ParseAtomRecordLine($RecordLine);
    $ResidueID = "${ResidueName}_${ResidueNumber}_${ChainID}";

    if ($ResidueID !~ /^$SpecifiedResidueID$/i) {
      # It's a new residue line...
      last RECORDLINE;
    }

    # Check distance...
    @Coords = (); push @Coords, ($X, $Y, $Z);
    if (CheckDistance(\@Coords, $OriginCoordsRef)) {
      # Distance criterion is met for at least one record in the residue...
      $Status = 1;
      last RECORDLINE;
    }
  }
  return $Status;
}

# Does record meets distance citerion specified by the user?
sub CheckDistance {
  my($CoordsRef, $OriginCoordsRef) = @_;
  my($Status, $Index, $Distance, $DistanceSquare);

  $Status = 0;

  if ($OptionsInfo{ExtractionDistanceMode} =~ /^Residue$/i) {
    # Go over coordinates of all the atoms in the residue...
    my($ResidueCoordsCount) = scalar @{$OriginCoordsRef};
    INDEX: for ($Index = 0; $Index < $ResidueCoordsCount; $Index += 3) {
      $DistanceSquare = ($CoordsRef->[0] - $OriginCoordsRef->[$Index])**2 + ($CoordsRef->[1] - $OriginCoordsRef->[$Index + 1])**2 + ($CoordsRef->[2] - $OriginCoordsRef->[$Index + 2])**2;
      $Distance = sqrt $DistanceSquare;
      if ($Distance <= $OptionsInfo{MaxExtractionDistance}) {
	$Status = 1;
	last INDEX;
      }
    }
  }
  else {
    $DistanceSquare = 0;
    for $Index (0 .. 2) {
      $DistanceSquare += ($CoordsRef->[$Index] - $OriginCoordsRef->[$Index])**2;
    }
    $Distance = sqrt $DistanceSquare;
    $Status = ($Distance <= $OptionsInfo{MaxExtractionDistance}) ? 1 : 0;
  }

  return $Status;
}

# Write out modifed header and other older records...
sub WriteHeaderAndOlderRecords {
  my($OutFileRef, $PDBRecordLinesRef) = @_;

  if ($OptionsInfo{ModifyHeaderRecord}) {
    # Write out modified HEADER record...
    my($Classification, $DepositionDate, $IDCode) = GetHeaderRecordInformation($PDBRecordLinesRef);
    $Classification = 'Data extracted using MayaChemTools';
    print $OutFileRef GenerateHeaderRecordLine($IDCode, $Classification), "\n";
  }
  else {
    print $OutFileRef $PDBRecordLinesRef->[0], "\n";
  }

  # Write out any old records...
  if ($OptionsInfo{KeepOldRecords}) {
    my($RecordLineIndex, $RecordLine);
    # Skip HEADER record and write out older records all the way upto first MODEL/ATOM/HETATM records from input file...
    RECORDLINE: for $RecordLineIndex (1 .. $#{$PDBRecordLinesRef}) {
      $RecordLine = $PDBRecordLinesRef->[$RecordLineIndex];
      if (IsModelRecordType($RecordLine) || IsAtomRecordType($RecordLine) || IsHetatmRecordType($RecordLine)) {
	last RECORDLINE;
      }
      print $OutFileRef "$RecordLine\n";
    }
  }
}

# Get header record information assuming it's the first record...
sub GetHeaderRecordInformation {
  my($PDBRecordLinesRef) = @_;
  my($Classification, $DepositionDate, $IDCode, $HeaderRecordLine);

  ($Classification, $DepositionDate, $IDCode) = ('') x 3;
  $HeaderRecordLine = $PDBRecordLinesRef->[0];
  if (IsHeaderRecordType($HeaderRecordLine)) {
    ($Classification, $DepositionDate, $IDCode) = ParseHeaderRecordLine($HeaderRecordLine);
  }
  return ($Classification, $DepositionDate, $IDCode);
}

# Get one letter residue code...
sub GetResidueCode {
  my($ResidueName) = @_;
  my($ResidueCode, $StandardResidue);

  $ResidueCode = $OptionsInfo{NonStandardSequenceCode};
  $StandardResidue = 0;

  if (length($ResidueName) == 3) {
    # Assume it's an amino acid...
    if (AminoAcids::IsAminoAcid($ResidueName)) {
      # Standard amino acid...
      $ResidueCode = AminoAcids::GetAminoAcidOneLetterCode($ResidueName);
      $StandardResidue = 1;
    }
  }
  elsif (length($ResidueName) == 1) {
    # Assume it's a nucleic acid...
    if ($ResidueName =~ /^(A|G|T|U|C)$/i) {
      $ResidueCode = $ResidueName;
      $StandardResidue = 1;
    }
  }

  return ($ResidueCode, $StandardResidue);
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();
  $OptionsInfo{Mode} = $Options{mode};

  my(@SpecifiedChains) = ();
  if ($Options{chains} =~ /^(First|All)$/i) {
    $OptionsInfo{ChainsToExtract} = $Options{chains};
  }
  else {
    @SpecifiedChains = split /\,/, $Options{chains};
    $OptionsInfo{ChainsToExtract} = 'Specified';
  }
  @{$OptionsInfo{SpecifiedChains}} = ();
  push @{$OptionsInfo{SpecifiedChains}}, @SpecifiedChains;

  $OptionsInfo{ChainsRecordMode} = $Options{chainsrecordmode};

  $OptionsInfo{CombineChains} = ($Options{combinechains} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{CombineChainSequences} = ($Options{combinechains} =~ /^Yes$/i) ? 1 : 0;

  ProcessResiduesOptions();
  ProcessAtomsOptions();
  ProcessDistanceOptions();

  $OptionsInfo{WaterResidueNames} = $Options{waterresiduenames};
  @{$OptionsInfo{SpecifiedWaterResiduesList}} = ();
  %{$OptionsInfo{SpecifiedWaterResiduesMap}} = ();

  my(@SpecifiedWaterResiduesList);
  @SpecifiedWaterResiduesList = ();

  if ($OptionsInfo{Mode} =~ /^NonWater$/i) {
    my($WaterResidueName);
    if ($OptionsInfo{WaterResidueNames} =~ /Automatic/i) {
      push @SpecifiedWaterResiduesList, ('HOH', 'WAT', 'H2O');
    }
    else {
      @SpecifiedWaterResiduesList = split /\,/, $Options{waterresiduenames};
    }
    for $WaterResidueName (@SpecifiedWaterResiduesList) {
      $OptionsInfo{SpecifiedWaterResiduesMap}{$WaterResidueName} = $WaterResidueName;
    }
  }
  push @{$OptionsInfo{SpecifiedWaterResiduesList}}, @SpecifiedWaterResiduesList;

  $OptionsInfo{RecordMode} = $Options{recordmode} ? $Options{recordmode} : ($Options{mode} =~ /^(Atoms|CAlphas|AtomNums|AtomsRange|AtomNames)$/i ? "Atom" : "AtomAndHetatm");

  $OptionsInfo{KeepOldRecords} = ($Options{keepoldrecords} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{ModifyHeaderRecord} = ($Options{modifyheader} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{KeepNonStandardSequences} = ($Options{nonstandardkeep} =~ /^Yes$/i) ? 1 : 0;
  $OptionsInfo{NonStandardSequenceCode} = $Options{nonstandardcode};
  $OptionsInfo{MaxSequenceLength} = $Options{sequencelength};
  $OptionsInfo{SequenceRecordSource} = $Options{sequencerecords};
  $OptionsInfo{SequenceIDPrefixSource} = $Options{sequenceidprefix};

  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? 1 : 0;
  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : 0;
}

# Process specified residue options...
sub ProcessResiduesOptions {
  my($ResidueNum, $StartResidueNum, $EndResNum, $ResidueName, @SpecifiedResidueNumsList, @SpecifiedResidueNamesList);

  @SpecifiedResidueNumsList = ();
  ($StartResidueNum, $EndResNum) = (0, 0);

  @SpecifiedResidueNamesList = ();

  if ($OptionsInfo{Mode} =~ /^(ResidueNums|ResiduesRange|ResidueNames)$/i) {
    if (!$Options{residues}) {
      die "Error: You must specify a value for \"--Residues\" option in \"ResidueNums, ResiduesRange, or ResidueNames\" \"-m, --mode\". \n";
    }
    $OptionsInfo{Residues} = $Options{residues};
    $OptionsInfo{Residues} =~ s/ //g;

    if ($OptionsInfo{Mode} =~ /^ResidueNames$/i) {
      @SpecifiedResidueNamesList = split /\,/, $OptionsInfo{Residues};
    }
    else {
      @SpecifiedResidueNumsList = split /\,/, $OptionsInfo{Residues};
      for $ResidueNum (@SpecifiedResidueNumsList) {
	if (!IsPositiveInteger($ResidueNum)) {
	  die "Error: Invalid residue number value, $ResidueNum, for \"--Residues\" option during \"ResidueNumes\" or \"ResiduesRange\"value of \"-m --mode\" option: Residue number must be a positive integer.\n";
	}
      }
      if ($OptionsInfo{Mode} =~ /^ResiduesRange$/i) {
	if (@SpecifiedResidueNumsList != 2) {
	  die "Error: Invalid number of residue number values, ", scalar(@SpecifiedResidueNumsList), ", for \"--Residues\" option during \"ResiduesRange\" value of \"-m --mode\" option: The number of values must be 2 corresponding to start and end residue numbers.\n";
	}
	if ($SpecifiedResidueNumsList[0] > $SpecifiedResidueNumsList[1]) {
	  die "Error: Invalid residue number values, @SpecifiedResidueNumsList, for \"--Residues\" option during \"ResiduesRange\" value of \"-m --mode\" option: The start residue number must be less than end residue number.\n";
	}
	($StartResidueNum, $EndResNum) = @SpecifiedResidueNumsList;
      }
    }
  }

  @{$OptionsInfo{SpecifiedResidueNumsList}} = ();
  push @{$OptionsInfo{SpecifiedResidueNumsList}}, @SpecifiedResidueNumsList;

  $OptionsInfo{SpecifiedStartResidueNum} = $StartResidueNum;
  $OptionsInfo{SpecifiedEndResidueNum} = $EndResNum;

  @{$OptionsInfo{SpecifiedResidueNamesList}} = ();
  push @{$OptionsInfo{SpecifiedResidueNamesList}}, @SpecifiedResidueNamesList;

  # Set up a specified residue numbers map...
  %{$OptionsInfo{SpecifiedResidueNumsMap}} = ();
  for $ResidueNum (@{$OptionsInfo{SpecifiedResidueNumsList}}) {
    $OptionsInfo{SpecifiedResidueNumsMap}{$ResidueNum} = $ResidueNum;
  }

  # Set up a specified residue names map...
  %{$OptionsInfo{SpecifiedResidueNamesMap}} = ();
  for $ResidueName (@{$OptionsInfo{SpecifiedResidueNamesList}}) {
    $OptionsInfo{SpecifiedResidueNamesMap}{lc $ResidueName} = lc $ResidueName;
  }

}

# Process specified atom options...
sub ProcessAtomsOptions {
  my($AtomNum, $StartAtomNum, $EndAtomNum, $AtomName, @SpecifiedAtomNumsList, @SpecifiedAtomNamesList);

  @SpecifiedAtomNumsList = ();
  ($StartAtomNum, $EndAtomNum) = (0, 0);

  @SpecifiedAtomNamesList = ();

  if ($OptionsInfo{Mode} =~ /^(AtomNums|AtomsRange|AtomNames)$/i) {
    if (!$Options{atoms}) {
      die "Error: You must specify a value for \"--Atoms\" option in \"AtomNums, AtomsRange, or AtomNames\" \"-m, --mode\". \n";
    }
    $OptionsInfo{Atoms} = $Options{atoms};
    $OptionsInfo{Atoms} =~ s/ //g;

    if ($OptionsInfo{Mode} =~ /^AtomNames$/i) {
      @SpecifiedAtomNamesList = split /\,/, $OptionsInfo{Atoms};
    }
    else {
      @SpecifiedAtomNumsList = split /\,/, $OptionsInfo{Atoms};
      for $AtomNum (@SpecifiedAtomNumsList) {
	if (!IsPositiveInteger($AtomNum)) {
	  die "Error: Invalid atom number value, $AtomNum, for \"--Atoms\" option during \"AtomNums\" or \"AtomsRange\"value of \"-m --mode\" option: Atom number must be a positive integer.\n";
	}
      }
      if ($OptionsInfo{Mode} =~ /^AtomsRange$/i) {
	if (@SpecifiedAtomNumsList != 2) {
	  die "Error: Invalid number of atom number values, ", scalar(@SpecifiedAtomNumsList), ", for \"--Atoms\" option during \"AtomsRange\" value of \"-m --mode\" option: The number of values must be 2 corresponding to start and end atom numbers.\n";
	}
	if ($SpecifiedAtomNumsList[0] > $SpecifiedAtomNumsList[1]) {
	  die "Error: Invalid atom number values, @SpecifiedAtomNumsList, for \"--Atoms\" option during \"AtomsRange\" value of \"-m --mode\" option: The start atom number must be less than end atom number.\n";
	}
	($StartAtomNum, $EndAtomNum) = @SpecifiedAtomNumsList;
      }
    }
  }
  elsif ($OptionsInfo{Mode} =~ /^CAlphas$/i) {
    @SpecifiedAtomNamesList = ("CA");
  }

  @{$OptionsInfo{SpecifiedAtomNumsList}} = ();
  push @{$OptionsInfo{SpecifiedAtomNumsList}}, @SpecifiedAtomNumsList;

  $OptionsInfo{SpecifiedStartAtomNum} = $StartAtomNum;
  $OptionsInfo{SpecifiedEndAtomNum} = $EndAtomNum;

  @{$OptionsInfo{SpecifiedAtomNamesList}} = ();
  push @{$OptionsInfo{SpecifiedAtomNamesList}}, @SpecifiedAtomNamesList;

  # Set up a specified residue numbers map...
  %{$OptionsInfo{SpecifiedAtomNumsMap}} = ();
  for $AtomNum (@{$OptionsInfo{SpecifiedAtomNumsList}}) {
    $OptionsInfo{SpecifiedAtomNumsMap}{$AtomNum} = $AtomNum;
  }

  # Set up a specified residue names map...
  %{$OptionsInfo{SpecifiedAtomNamesMap}} = ();
  for $AtomName (@{$OptionsInfo{SpecifiedAtomNamesList}}) {
    $OptionsInfo{SpecifiedAtomNamesMap}{lc $AtomName} = lc $AtomName;
  }

}

# Process specified distance options...
sub ProcessDistanceOptions {
  my(@SpecifiedDistanceOrigin) = ();

  $OptionsInfo{MaxExtractionDistance} = $Options{distance};
  $OptionsInfo{ExtractionDistanceMode} = $Options{distancemode};
  $OptionsInfo{ExtractionDistanceOrigin} = $Options{distanceorigin} ? $Options{distanceorigin} : '';
  $OptionsInfo{DistanceSelectionMode} = $Options{distanceselectionmode};

  if ($OptionsInfo{Mode} =~ /^Distance$/i) {
    if (!$Options{distanceorigin}) {
      die "Error: You must specify a value for \"--distanceorigin\" option in \"Distance\" \"-m, --mode\". \n";
    }
    @SpecifiedDistanceOrigin = split /\,/, $Options{distanceorigin};
    if ($OptionsInfo{ExtractionDistanceMode} =~ /^Atom$/i) {
      if (@SpecifiedDistanceOrigin != 2) {
	die "Error: Invalid number of values, ", scalar(@SpecifiedDistanceOrigin), " for option \"distanceorigin\" option during \"Atom\" value of \"--distancemode\" : The number of values must be 2.\n";
      }
      if (!IsPositiveInteger($SpecifiedDistanceOrigin[0])) {
	die "Error: Invalid atom number value, ", $SpecifiedDistanceOrigin[0], ", for option \"distanceorigin\" option during \"Atom\" value of \"--distancemode\". Allowed values: > 0\n";
      }
    }
    elsif ($OptionsInfo{ExtractionDistanceMode} =~ /^Hetatm$/i) {
      if (@SpecifiedDistanceOrigin != 2) {
	die "Error: Invalid number of values, ", scalar(@SpecifiedDistanceOrigin), " for option \"distanceorigin\" option during \"Hetatm\" value of \"--distancemode\" : The number of values must be 2.\n";
      }
      if (!IsPositiveInteger($SpecifiedDistanceOrigin[0])) {
	die "Error: Invalid hetatm number value, ", $SpecifiedDistanceOrigin[0], ", for option \"distanceorigin\" option during \"Hetatm\" value of \"--distancemode\". Allowed values: > 0\n";
      }
    }
    elsif ($OptionsInfo{ExtractionDistanceMode} =~ /^Residue$/i) {
      if (!(@SpecifiedDistanceOrigin == 2 || @SpecifiedDistanceOrigin == 3)) {
	die "Error: Invalid number of values, ", scalar(@SpecifiedDistanceOrigin), " for option \"distanceorigin\" option during \"Residue\" value of \"--distancemode\" : The number of values must be either 2 or 3.\n";
      }
      if (!IsPositiveInteger($SpecifiedDistanceOrigin[0])) {
	die "Error: Invalid residue number value, ", $SpecifiedDistanceOrigin[0], ", for option \"distanceorigin\" option during \"Residue\" value of \"--distancemode\". Allowed values: > 0\n";
      }
    }
    elsif ($OptionsInfo{ExtractionDistanceMode} =~ /^XYZ$/i) {
      if (@SpecifiedDistanceOrigin != 3) {
	die "Error: Invalid number of values, ", scalar(@SpecifiedDistanceOrigin), " for option \"distanceorigin\" option during \"XYZ\" value of \"--distancemode\" : The number of values must be 3.\n";
      }
      my($Value);
      for $Value (@SpecifiedDistanceOrigin) {
	if (!IsNumerical($Value)) {
	  die "Error: Invalid coordinate value, ", $SpecifiedDistanceOrigin[0], ", for option \"distanceorigin\" option during \"XYZ\" value of \"--distancemode\". Allowed values: numerical\n";
	}
      }
    }
  }
  @{$OptionsInfo{SpecifiedExtractionDistanceOrigin}} = ();
  push @{$OptionsInfo{SpecifiedExtractionDistanceOrigin}}, @SpecifiedDistanceOrigin;

}

# Retrieve information about PDB files...
sub RetrievePDBFilesInfo {
  my($Index, $PDBFile, $PDBRecordLinesRef, $ChainID, $ChainLabel, $ChainsAndResiduesInfoRef, $Mode, $FileDir, $FileName, $FileExt, $OutFileName, $OutFileRoot, @SpecifiedChains, @DistanceOrigin, @OutFileNames, @ChainLabels, @ChainSequenceIDs, @ChainSequenceIDsPrefix);

  %PDBFilesInfo = ();
  @{$PDBFilesInfo{FileOkay}} = ();
  @{$PDBFilesInfo{OutFileRoot}} = ();
  @{$PDBFilesInfo{OutFileNames}} = ();
  @{$PDBFilesInfo{ChainLabels}} = ();
  @{$PDBFilesInfo{ChainSequenceIDs}} = ();
  @{$PDBFilesInfo{ChainSequenceIDsPrefix}} = ();
  @{$PDBFilesInfo{SpecifiedChains}} = ();
  @{$PDBFilesInfo{DistanceOrigin}} = ();

  FILELIST: for $Index (0 .. $#PDBFilesList) {
    $PDBFilesInfo{FileOkay}[$Index] = 0;

    $PDBFilesInfo{OutFileRoot}[$Index] = '';
    @{$PDBFilesInfo{OutFileNames}[$Index]} = ();
    @{$PDBFilesInfo{OutFileNames}[$Index]} = ();
    @{$PDBFilesInfo{ChainLabels}[$Index]} = ();
    @{$PDBFilesInfo{ChainSequenceIDs}[$Index]} = ();
    @{$PDBFilesInfo{ChainSequenceIDsPrefix}[$Index]} = ();
    @{$PDBFilesInfo{SpecifiedChains}[$Index]} = ();
    @{$PDBFilesInfo{DistanceOrigin}[$Index]} = ();

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

    # Get PDB data...
    $PDBRecordLinesRef = ReadPDBFile($PDBFile);
    if ($OptionsInfo{Mode} =~ /^Sequences$/i && $OptionsInfo{SequenceRecordSource} =~ /^SeqRes$/i) {
      $ChainsAndResiduesInfoRef = GetChainsAndResidues($PDBRecordLinesRef, 'SeqRes');
    }
    else {
      $ChainsAndResiduesInfoRef = GetChainsAndResidues($PDBRecordLinesRef);
    }
    if (!scalar @{$ChainsAndResiduesInfoRef->{ChainIDs}}) {
      warn "Warning: Ignoring file $PDBFile: No chains found \n";
      next FILELIST;
    }

    # Make sure specified chains exist in PDB file...
    @SpecifiedChains = ();
    if ($OptionsInfo{ChainsToExtract} =~ /^Specified$/i) {
      for $ChainID (@{$OptionsInfo{SpecifiedChains}}) {
	if (exists $ChainsAndResiduesInfoRef->{Residues}{$ChainID}) {
	  push @SpecifiedChains, $ChainID;
	}
	else {
	  warn "Warning: Ignoring file $PDBFile: Specified chain, $ChainID, in \"-c, --chains\" option doesn't exist.\n";
	  next FILELIST;
	}
      }
    }
    elsif ($OptionsInfo{ChainsToExtract} =~ /^First$/i) {
      push @SpecifiedChains, $ChainsAndResiduesInfoRef->{ChainIDs}[0];
    }
    elsif ($OptionsInfo{ChainsToExtract} =~ /^All$/i) {
      push @SpecifiedChains, @{$ChainsAndResiduesInfoRef->{ChainIDs}};
    }
    # Setup chain labels to use for sequence IDs and generating output files...
    @ChainLabels = ();
    for $ChainID (@SpecifiedChains) {
      $ChainLabel = $ChainID; $ChainLabel =~ s/^None//ig;
      $ChainLabel = "Chain${ChainLabel}";
      push @ChainLabels, $ChainLabel;
    }

    # Make sure specified distance origin is valid...
    @DistanceOrigin = ();
    if ($OptionsInfo{Mode} =~ /^Distance$/i) {
      if ($OptionsInfo{ExtractionDistanceMode} =~ /^(Atom|Hetatm)$/i) {
	my($RecordType, $SpecifiedAtomName, $SpecifiedAtomNumber, $RecordFound, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $RecordLine);
	$RecordType = $OptionsInfo{ExtractionDistanceMode};
	($SpecifiedAtomNumber, $SpecifiedAtomName) = @{$OptionsInfo{SpecifiedExtractionDistanceOrigin}};
	$RecordFound = 0;
	LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
	  if (!(IsAtomRecordType($RecordLine) || IsHetatmRecordType($RecordLine))) {
	    next LINE;
	  }
	  ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z) = ParseAtomRecordLine($RecordLine);
	  $AtomName = RemoveLeadingAndTrailingWhiteSpaces($AtomName);
	  if (($RecordType =~ /^Atom$/i && IsAtomRecordType($RecordLine)) || ($RecordType =~ /^Hetatm$/i && IsHetatmRecordType($RecordLine))) {
	    if ($AtomNumber == $SpecifiedAtomNumber && $AtomName eq $SpecifiedAtomName) {
	      $RecordFound = 1;
	      last LINE;
	    }
	  }
	}
	if (!$RecordFound) {
	  warn "Warning: Ignoring file $PDBFile: ", uc($RecordType), " record corresponding to \"--distanceorigin\" option value, $OptionsInfo{ExtractionDistanceOrigin}, doesn't exist.\n";
	  next FILELIST;
	}
	push @DistanceOrigin, ($X, $Y, $Z);
      }
      elsif ($OptionsInfo{ExtractionDistanceMode} =~ /^Residue$/i) {
	my($SpecifiedResidueNumber, $SpecifiedResidueName, $SpecifiedChainID, $RecordFound, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $RecordLine);
	$SpecifiedChainID = '';
	if (@{$OptionsInfo{SpecifiedExtractionDistanceOrigin}} == 3) {
	  ($SpecifiedResidueNumber, $SpecifiedResidueName, $SpecifiedChainID) = @{$OptionsInfo{SpecifiedExtractionDistanceOrigin}};
	}
	else {
	  ($SpecifiedResidueNumber, $SpecifiedResidueName) = @{$OptionsInfo{SpecifiedExtractionDistanceOrigin}};
	}
	$RecordFound = 0;
	LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
	  if (!(IsAtomRecordType($RecordLine) || IsHetatmRecordType($RecordLine))) {
	    next LINE;
	  }
	  ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z) = ParseAtomRecordLine($RecordLine);
	  $ResidueName = RemoveLeadingAndTrailingWhiteSpaces($ResidueName);
	  $ChainID = RemoveLeadingAndTrailingWhiteSpaces($ChainID);
	  if ($SpecifiedChainID && ($SpecifiedChainID ne $ChainID)) {
	    next LINE;
	  }
	  if ($ResidueNumber == $SpecifiedResidueNumber && $ResidueName eq $SpecifiedResidueName) {
	    # Store coordinates for all the atoms...
	    $RecordFound = 1;
	    push @DistanceOrigin, ($X, $Y, $Z);
	    next LINE;
	  }
	}
	if (!$RecordFound) {
	  warn "Warning: Ignoring file $PDBFile: ATOM/HETATM record corresponding to \"--distanceorigin\" option value, $OptionsInfo{ExtractionDistanceOrigin}, doesn't exist.\n";
	  next FILELIST;
	}
      }
      elsif ($OptionsInfo{ExtractionDistanceMode} =~ /^XYZ$/i) {
	push @DistanceOrigin, @{$OptionsInfo{SpecifiedExtractionDistanceOrigin}};
      }
    }
    # Setup output file names...
    @OutFileNames = ();
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($PDBFile);
    if ($OptionsInfo{OutFileRoot} && (@PDBFilesList == 1)) {
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
    $Mode = $OptionsInfo{Mode};
    if ($Mode =~ /^(Atoms|CAlphas|AtomNums|AtomsRange|AtomNames|ResidueNums|ResiduesRange|ResidueNames|Distance|NonWater|NonHydrogens)$/i) {
      $OutFileName = '';
      if ($Mode =~ /^CAlphas$/i) {
	$OutFileName = "${OutFileRoot}CAlphas.pdb";
      }
      elsif ($Mode =~ /^Atoms$/i) {
	$OutFileName = "${OutFileRoot}Atoms.pdb";
      }
      elsif ($Mode =~ /^AtomNums$/i) {
	$OutFileName = "${OutFileRoot}AtomNums.pdb";
      }
      elsif ($Mode =~ /^AtomsRange$/i) {
	$OutFileName = "${OutFileRoot}AtomsRange.pdb";
      }
      elsif ($Mode =~ /^AtomNames$/i) {
	$OutFileName = "${OutFileRoot}AtomNames.pdb";
      }
      elsif ($Mode =~ /^ResidueNums$/i) {
	$OutFileName = "${OutFileRoot}ResidueNums.pdb";
      }
      elsif ($Mode =~ /^ResiduesRange$/i) {
	$OutFileName = "${OutFileRoot}ResiduesRange.pdb";
      }
      elsif ($Mode =~ /^ResidueNames$/i) {
	$OutFileName = "${OutFileRoot}ResidueNames.pdb";
      }
      elsif ($Mode =~ /^NonWater$/i) {
	$OutFileName = "${OutFileRoot}NonWater.pdb";
      }
      elsif ($Mode =~ /^NonHydrogens$/i) {
	$OutFileName = "${OutFileRoot}NonHydrogens.pdb";
      }
      elsif ($Mode =~ /^Distance$/i) {
	my($DistanceMode) = '';
	if ($OptionsInfo{ExtractionDistanceMode} =~ /^Atom$/i) {
	  $DistanceMode = 'Atom';
	}
	elsif ($OptionsInfo{ExtractionDistanceMode} =~ /^Hetatm$/i) {
	  $DistanceMode = 'Hetatm';
	}
	elsif ($OptionsInfo{ExtractionDistanceMode} =~ /^Residue$/i) {
	  $DistanceMode = 'Residue';
	}
	elsif ($OptionsInfo{ExtractionDistanceMode} =~ /^XYZ$/i) {
	  $DistanceMode = 'XYZ';
	}
	$OutFileName = "${OutFileRoot}DistanceBy${DistanceMode}.pdb";
      }
      push @OutFileNames, $OutFileName;
      if (!$OptionsInfo{OverwriteFiles} && (-e $OutFileName)) {
	warn "Warning: Ignoring file $PDBFile: The file $OutFileName already exists\n";
	next FILELIST;
      }
    }
    elsif ($Mode =~ /^(Chains|Sequences)$/i) {
      if ($OptionsInfo{CombineChainSequences}) {
	$OutFileName = ($Mode =~ /^Chains$/i) ? "${OutFileRoot}ExtractedChains.pdb" : "${OutFileRoot}SequencesChainsCombined.fasta";
	push @OutFileNames, $OutFileName;
	if (!$OptionsInfo{OverwriteFiles} && (-e $OutFileName)) {
	  warn "Warning: Ignoring file $PDBFile: The file $OutFileName already exists\n";
	  next FILELIST;
	}
      }
      else {
	for $ChainLabel (@ChainLabels) {
	  $OutFileName = ($Mode =~ /^Chains$/i) ? "${OutFileRoot}${ChainLabel}.pdb" : "${OutFileRoot}Sequences${ChainLabel}.fasta";
	  push @OutFileNames, $OutFileName;
	  if (!$OptionsInfo{OverwriteFiles} && (-e $OutFileName)) {
	    warn "Warning: Ignoring file $PDBFile: The file $OutFileName already exists\n";
	    next FILELIST;
	  }
	}
      }
    }
    @ChainSequenceIDs = ();
    @ChainSequenceIDsPrefix = ();
    if ($Mode =~ /^Sequences$/i) {
      my($HeaderRecordLine, $Classification, $DepositionDate, $IDCode, $IDPrefix);
      ($Classification, $DepositionDate, $IDCode) = GetHeaderRecordInformation($PDBRecordLinesRef);

      if ($OptionsInfo{SequenceIDPrefixSource} =~ /^FileName$/i) {
	$IDPrefix = $FileName;
      }
      elsif ($OptionsInfo{SequenceIDPrefixSource} =~ /^HeaderRecord$/i) {
	$IDPrefix = IsNotEmpty($IDCode) ? $IDCode : '';
      }
      else {
	$IDPrefix = IsNotEmpty($IDCode) ? $IDCode : $FileName;
      }

      for $ChainLabel (@ChainLabels) {
	push @ChainSequenceIDsPrefix, $IDPrefix;
	push @ChainSequenceIDs, "${IDPrefix}_${ChainLabel}|PDB";
      }
    }

    $PDBFilesInfo{FileOkay}[$Index] = 1;
    $PDBFilesInfo{OutFileRoot}[$Index] = $OutFileRoot;

    push @{$PDBFilesInfo{OutFileNames}[$Index]}, @OutFileNames;
    push @{$PDBFilesInfo{ChainLabels}[$Index]}, @ChainLabels;
    push @{$PDBFilesInfo{ChainSequenceIDsPrefix}[$Index]}, @ChainSequenceIDsPrefix;
    push @{$PDBFilesInfo{ChainSequenceIDs}[$Index]}, @ChainSequenceIDs;
    push @{$PDBFilesInfo{SpecifiedChains}[$Index]}, @SpecifiedChains;
    push @{$PDBFilesInfo{DistanceOrigin}[$Index]}, @DistanceOrigin;
  }
}


# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{chains} = 'First';
  $Options{chainsrecordmode} = 'NotAcrossTER';
  $Options{combinechains} = 'no';
  $Options{distance} = 10.0;
  $Options{distancemode} = 'XYZ';
  $Options{distanceselectionmode} = 'ByAtom';
  $Options{keepoldrecords} = 'no';
  $Options{mode} = 'NonWater';
  $Options{modifyheader} = 'yes';
  $Options{nonstandardkeep} = 'yes';
  $Options{nonstandardcode} = 'X';
  $Options{sequencelength} = 80;
  $Options{sequenceidprefix} = 'Automatic';
  $Options{sequencerecords} = 'Atom';
  $Options{waterresiduenames} = 'Automatic';

  if (!GetOptions(\%Options, "atoms|a=s", "chains|c=s", "chainsrecordmode=s", "combinechains=s", "distance|d=f", "distancemode=s", "distanceorigin=s", "distanceselectionmode=s", "help|h", "keepoldrecords|k=s", "mode|m=s", "modifyheader=s", "nonstandardkeep=s", "nonstandardcode=s", "overwrite|o", "root|r=s", "recordmode=s", "residues=s", "sequencelength=i", "sequenceidprefix=s", "sequencerecords=s", "waterresiduenames=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{combinechains} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{combinechains}, for option \"--CombineChains\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{chainsrecordmode} !~ /^(AcrossTER|NotAcrossTER)$/i) {
    die "Error: The value specified, $Options{chainsrecordmode}, for option \"--ChainsRecordMode\" is not valid. Allowed values: AcrossTER or NotAcrossTER\n";
  }
  if ($Options{distancemode} !~ /^(Atom|Hetatm|Residue|XYZ)$/i) {
    die "Error: The value specified, $Options{distancemode}, for option \"--DistanceMode\" is not valid. Allowed values: Atom, Hetatm, Residue, or XYZ\n";
  }
  if ($Options{distanceselectionmode} !~ /^(ByAtom|ByResidue)$/i) {
    die "Error: The value specified, $Options{distanceselectionmode}, for option \"--DistanceSelectionMode\" is not valid. Allowed values: ByAtom or ByResidue\n";
  }
  if ($Options{keepoldrecords} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{keepoldrecords}, for option \"--KeepOldRecords\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{mode} !~ /^(Chains|Sequences|Atoms|CAlphas|AtomNums|AtomsRange|AtomNames|ResidueNums|ResidueNames|ResiduesRange|Distance|NonWater|NonHydrogens)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"m, --mode\" is not valid. Allowed values: Chains, Sequences, Atoms, CAlphas, AtomNums, AtomsRange, AtomNames, ResidueNums, ResiduesRange, ResidueNames, Distance, NonWater, NonHydrogens\n";
  }
  if ($Options{modifyheader} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{modifyheader}, for option \"--ModifyHeader\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{nonstandardkeep} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{nonstandardkeep}, for option \"--NonStandardKeep\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{nonstandardcode} !~ /^(\?|\-|X)$/i) {
    die "Error: The value specified, $Options{nonstandardcode}, for option \"--NonStandardCode\" is not valid. Allowed values: ?, -, or X\n";
  }
  if ($Options{recordmode} && $Options{recordmode} !~ /^(Atom|Hetatm|AtomAndHetatm)$/i) {
    die "Error: The value specified, $Options{recordmode}, for option \"--RecordMode\" is not valid. Allowed values: Atom, Hetatm, AtomAndHetatm\n";
  }
  if (!IsPositiveInteger($Options{sequencelength})) {
    die "Error: The value specified, $Options{sequencelength}, for option \"--SequenceLength\" is not valid. Allowed values: >0\n";
  }
  if ($Options{sequencerecords} !~ /^(Atom|SeqRes)$/i) {
    die "Error: The value specified, $Options{sequencerecords}, for option \"--SequenceRecords\" is not valid. Allowed values: Atom or SeqRes\n";
  }
  if ($Options{sequenceidprefix} !~ /^(FileName|HeaderRecord|Automatic)$/i) {
    die "Error: The value specified, $Options{sequenceidprefix}, for option \"--SequenceIDPrefix\" is not valid. Allowed values: FileName, HeaderRecord, or AutomaticAtom\n";
  }
}

__END__

=head1 NAME

ExtractFromPDBFiles.pl - Extract specific data from PDBFile(s)

=head1 SYNOPSIS

ExtractFromPDBFiles.pl PDBFile(s)...

ExtractFromPDBFiles.pl [B<-a, --Atoms> "AtomNum, [AtomNum...]" | "StartAtomNum, EndAtomNum" |
"AtomName, [AtomName...]"] [B<-c, --chains> First | All | "ChainID, [ChainID,...]"] [B<--ChainsRecordMode> I<AcrossTER | NotAcrossTER>]
[<--CombineChains> yes | no] [B<-d, --distance> number] [B<--DistanceMode> Atom | Hetatm | Residue | XYZ]
[B<--DistanceOrigin> "AtomNumber, AtomName" | "HetatmNumber, HetAtmName" | "ResidueNumber, ResidueName, [ChainID]" | "X,Y,Z">]
[<--DistanceSelectionMode> ByAtom | ByResidue] [B<-h, --help>] [B<-k, --KeepOldRecords> yes | no]
[B<-m, --mode > Chains | Sequences | Atoms | CAlphas | AtomNums | AtomsRange | AtomNames |
ResidueNums | ResiduesRange | ResidueNames | Distance | NonWater | NonHydrogens]
[B<--ModifyHeader> yes | no] [B<--NonStandardKeep> yes | no] [B<--NonStandardCode> character]
[B<-o, --overwrite>] [B<-r, --root> rootname] B<--RecordMode> I<Atom | Hetatm | AtomAndHetatm>]
[B<--Residues> "ResidueNum,[ResidueNum...]" | StartResidueNum,EndResiduNum ]
[B<--SequenceLength> number] [B<--SequenceRecords> Atom | SeqRes]
[B<--SequenceIDPrefix> FileName | HeaderRecord | Automatic]
[B<--WaterResidueNames> Automatic | "ResidueName, [ResidueName,...]"]
[B<-w, --WorkingDir> dirname] PDBFile(s)...

=head1 DESCRIPTION

Extract specific data from I<PDBFile(s)> and generate appropriate PDB or sequence file(s).
Multiple PDBFile names are separated by spaces. The valid file extension is I<.pdb>.
All other file name extensions are ignored during the wild card expansion. All the PDB files
in a current directory can be specified either by I<*.pdb> or the current directory name.

During I<Chains> and I<Sequences> values of B<-m, --mode> option, all ATOM/HETAM records
for chains after the first model in PDB fils containing data for multiple models are ignored.

=head1 OPTIONS

=over 4

=item B<-a, --Atoms> I<"AtomNum,[AtomNum...]" | "StartAtomNum,EndAtomNum" | "AtomName,[AtomName...]">

Specify which atom records to extract from I<PDBFiles(s)> during I<AtomNums>,
I<AtomsRange>, and I<AtomNames> value of B<-m, --mode> option: extract records
corresponding to atom numbers specified in a comma delimited list of atom numbers/names,
or with in the range of start and end atom numbers. Possible values: I<"AtomNum[,AtomNum,..]">,
I<StartAtomNum,EndAtomNum>, or I<"AtomName[,AtomName,..]">. Default: I<None>. Examples:

    10
    15,20
    N,CA,C,O

=item B<-c, --chains> I<First | All | ChainID,[ChainID,...]>

Specify which chains to extract from I<PDBFile(s)> during I<Chains | Sequences> value of
B<-m, --mode> option: first chain, all chains, or a specific list of comma delimited chain IDs.
Possible values: I<First | All | ChainID,[ChainID,...]>. Default: I<First>. Examples:

    A
    A,B
    All

=item B<--ChainsRecordMode> I<AcrossTER | NotAcrossTER>

Specify whether to extract ATOM and HETATM record lines across TER records from
I<PDBFile(s)> during I<Chains> value of B<-m, --mode> option. Possible values:
I<AcrossTER | NotAcrossTER>. Defaul value: I<NotAcrossTER>.

This option allows retrieval ATOM and HETATM record lines for a specific chain
which spread across TER record in I<PDBFile(s)>.

=item B<--CombineChains> I<yes | no>

Specify whether to combine extracted chains data into a single file during I<Chains> or
I<Sequences> value of B<-m, --mode> option. Possible values: I<yes | no>. Default: I<no>.

During I<Chains> value of <-m, --mode> option with I<Yes> value of <--CombineChains>,
extracted data for specified chains is written into a single file instead of individual file for each
chain.

During I<Sequences> value of <-m, --mode> option with I<Yes> value of <--CombineChains>,
residues sequences for specified chains are extracted and concatenated into a single sequence
file instead of  individual file for each chain.

=item B<-d, --distance> I<number>

Specify distance used to extract ATOM/HETATM recods during I<Distance> value of
B<-m, --mode> option. Default: I<10.0> angstroms.

B<--RecordMode> option controls type of record lines to extract from I<PDBFile(s)>:
ATOM, HETATM or both.

=item B<--DistanceMode> I<Atom | Hetatm | Residue | XYZ>

Specify how to extract ATOM/HETATM records from I<PDBFile(s)> during I<Distance> value of
B<-m, --mode> option: extract all the records within a certain distance specifed by B<-d, --distance>
from an atom or hetro atom record, a residue, or any artbitrary point. Possible values: I<Atom |
Hetatm | Residue | XYZ>. Default: I<XYZ>.

During I<Residue> value of B<--distancemode>, distance of ATOM/HETATM records is calculated from
all the atoms in the residue and the records are selected as long as any atom of the residue lies with
in the distace specified using B<-d, --distance> option.

B<--RecordMode> option controls type of record lines to extract from I<PDBFile(s)>:
ATOM, HETATM or both.

=item B<--DistanceSelectionMode> I<ByAtom | ByResidue>

Specify how how to extract ATOM/HETATM records from I<PDBFile(s)> during I<Distance> value of
B<-m, --mode> option for all values of B<--DistanceMode> option: extract only those ATOM/HETATM
records that meet specified distance criterion; extract all records corresponding to a residue as
long as one of the ATOM/HETATM record in the residue satisfies specified distance criterion. Possible
values: I<ByAtom, ByResidue>. Default value: I<ByAtom>.

=item B<--DistanceOrigin> I<"AtomNumber,AtomName" | "HetatmNumber,HetAtmName" | "ResidueNumber,ResidueName[,ChainID]" | "X,Y,Z">

This value is B<--distancemode> specific. In general, it identifies a point used to select
other ATOM/HETATMS with in a specific distance from this point.

For I<Atom> value of B<--distancemode>, this option corresponds to an atom specification.
Format: I<AtomNumber,AtomName>. Example:

    455,CA

For I<Hetatm> value of B<--distancemode>, this option corresponds to a hetatm specification.
Format: I<HetatmNumber,HetAtmName>. Example:

    5295,C1

For I<Residue> value of B<--distancemode>, this option corresponds to a residue specification.
Format: I<ResidueNumber, ResidueName[,ChainID]>. Example:

    78,MSE
    977,RET,A
    978,RET,B

For I<XYZ> value of B<--distancemode>, this option corresponds to a coordinate of an
arbitrary point. Format: I<X,Y,X>. Example:

    10.044,19.261,-4.292

B<--RecordMode> option controls type of record lines to extract from I<PDBFile(s)>:
ATOM, HETATM or both.

=item B<-h, --help>

Print this help message.

=item B<-k, --KeepOldRecords> I<yes | no>

Specify whether to transfer old non ATOM and HETATM records from input PDBFile(s) to new
PDBFile(s) during I<Chains | Atoms | HetAtms | CAlphas | Distance| NonWater | NonHydrogens>
value of B<-m --mode> option. By default, except for the HEADER record, all
other unnecessary non ATOM/HETATM records are dropped during the
generation of new PDB files. Possible values: I<yes | no>. Default: I<no>.

=item B<-m, --mode > I<Chains | Sequences | Atoms | CAlphas | AtomNums | AtomsRange | AtomNames | ResidueNums | ResiduesRange | ResidueNames | Distance | NonWater | NonHydrogens>

Specify what to extract from I<PDBFile(s)>: I<Chains> - retrieve records for
specified chains; I<Sequences> - generate sequence files for specific chains;
I<Atoms> - extract atom records; I<CAlphas> - extract atom records for alpha
carbon atoms; I<AtomNums> - extract atom records for specified atom numbers;
I<AtomsRange> - extract atom records between specified atom number range;
I<AtomNames> - extract atom records for specified atom names; I<ResidueNums>
- extract records for specified residue numbers; I<ResiduesRange> - extract records
for residues between specified residue number range; I<ResidueNames> - extract
records for specified residue names; I<Distance> - extract records with in a
certain distance from a specific position; I<NonWater> - extract records corresponding
to residues other than water; I<NonHydrogens> - extract non-hydrogen records.

Possible values: I<Chains, Sequences Atoms, CAlphas, AtomNums, AtomsRange,
AtomNames, ResidueNums, ResiduesRange, ResidueNames, Distance, NonWater,
NonHydrogens>. Default value: I<NonWater>

During the generation of new PDB files, unnecessay CONECT records are dropped.

For I<Chains> mode, data for appropriate chains specified by B<--c --chains> option
is extracted from I<PDBFile(s)> and placed into new PDB file(s).

For I<Sequences> mode, residues names using various sequence related options are
extracted for chains specified by B<--c --chains> option from I<PDBFile(s)> and
FASTA sequence file(s) are generated.

For I<Distance> mode, all ATOM/HETATM records with in a distance specified
by B<-d --distance> option from a specific atom, residue or a point indicated by
B<--distancemode> are extracted and placed into new PDB file(s).

For I<NonWater> mode, non water ATOM/HETATM record lines, identified using value of
B<--WaterResidueNames>, are extracted and written to new PDB file(s).

For I<NonHydrogens> mode, ATOM/HETATOM record lines containing element symbol
other than I<H> are extracted and written to new PDB file(s).

For all other options, appropriate ATOM/HETATM records are extracted to generate new
PDB file(s).

B<--RecordMode> option controls type of record lines to extract and process from
I<PDBFile(s)>: ATOM, HETATM or both.

=item B<--ModifyHeader> I<yes | no>

Specify whether to modify HEADER record during the generation of new PDB files
for B<-m, --mode> values of I<Chains | Atoms | CAlphas | Distance>. Possible values:
I<yes | no>.  Default: I<yes>. By default, Classification data is replaced by I<Data extracted
using MayaChemTools> before writing out HEADER record.

=item B<--NonStandardKeep> I<yes | no>

Specify whether to include and convert non-standard three letter residue codes into
a code specified using B<--nonstandardcode> option and include them into sequence file(s)
generated during I<Sequences> value of B<-m, --mode> option. Possible values: I<yes | no>.
Default: I<yes>.

A warning is also printed about the presence of non-standard residues. Any residue other
than standard 20 amino acids and 5 nucleic acid is considered non-standard; additionally,
HETATM residues in chains also tagged as non-standard.

=item B<--NonStandardCode> I<character>

A single character code to use for non-standard residues. Default: I<X>. Possible values:
I<?, -, or X>.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

New PDB and sequence file name is generated using the root: <Root><Mode>.<Ext>.
Default new file name: <PDBFileName>Chain<ChainID>.pdb for I<Chains> B<mode>;
<PDBFileName>SequenceChain<ChainID>.fasta for I<Sequences> B<mode>;
<PDBFileName>DistanceBy<DistanceMode>.pdb for I<Distance> B<-m, --mode>
<PDBFileName><Mode>.pdb for I<Atoms | CAlphas | NonWater | NonHydrogens> B<-m, --mode>
values. This option is ignored for multiple input files.

=item B<--RecordMode> I<Atom | Hetatm | AtomAndHetatm>

Specify type of record lines to extract and process from I<PDBFile(s)> during various
values of B<-m, --mode> option: extract only ATOM record lines; extract only HETATM
record lines; extract both ATOM and HETATM lines. Possible values: I<Atom | Hetatm
| AtomAndHetatm | XYZ>. Default during I<Atoms, CAlphas, AtomNums, AtomsRange,
AtomNames> values of B<-m, --mode> option: I<Atom>; otherwise: I<AtomAndHetatm>.

This option is ignored during I<Sequences> values of B<-m, --mode> option.

=item B<--Residues> I<"ResidueNum,[ResidueNum...]" | "StartResidueNum,EndResiduNum" | "ResidueName,[ResidueName...]">

Specify which resiude records to extract from I<PDBFiles(s)> during I<ResidueNums>,
I<ResiduesRange>,and I<ResidueNames> value of B<-m, --mode> option: extract records
corresponding to residue numbers specified in a comma delimited list of residue numbers/names,
or with in the range of start and end residue numbers. Possible values: I<"ResidueNum[,ResidueNum,..]">,
I<StartResidueNum,EndResiduNum>, or I<<"ResidueName[,ResidueName,..]">. Default: I<None>. Examples:

    20
    5,10
    TYR,SER,THR

B<--RecordMode> option controls type of record lines to extract from I<PDBFile(s)>:
ATOM, HETATM or both.

=item B<--SequenceLength> I<number>

Maximum sequence length per line in sequence file(s). Default: I<80>.

=item B<--SequenceRecords> I<Atom | SeqRes>

Specify which records to use for extracting residue names from I<PDBFiles(s)> during
I<Sequences> value of B<-m, --mode> option: use ATOM records to compile a list
of residues in a chain or parse SEQRES record to get a list of residues. Possible values:
I<Atom | SeqRes>. Default: I<Atom>.

=item B<--SequenceIDPrefix> I<FileName | HeaderRecord | Automatic>

Specify how to generate a prefix for sequence IDs during I<Sequences> value
of B<-m, --mode> option: use input file name prefix; retrieve PDB ID from HEADER record;
or automatically decide the method for generating the prefix. The chain IDs are also
appended to the prefix. Possible values: I<FileName | HeaderRecord | Automatic>.
Default: I<Automatic>

=item B<--WaterResidueNames> I<Automatic | "ResidueName,[ResidueName,...]">

Identification of water residues during I<NonWater> value of B<-m, --mode> option. Possible values:
I<Automatic | "ResidueName,[ResidueName,...]">. Default: I<Automatic> - corresponds
to "HOH,WAT,H20". You can also specify a different comma delimited list of residue names
to use for water.

=item B<-w, --WorkingDir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To extract non-water records from Sample2.pdb file and generate Sample2NonWater.pdb
file, type:

    % ExtractFromPDBFiles.pl Sample2.pdb

To extract non-water records corresponding to only ATOM records from Sample2.pdb file
and generate Sample2NonWater.pdb file, type:

    % ExtractFromPDBFiles.pl --RecordMode Atom Sample2.pdb

To extract non-water records from Sample2.pdb file using HOH or WAT residue name for water along
with all old non-coordinate records and generate Sample2NewNonWater.pdb file, type:

    % ExtractFromPDBFiles.pl -m NonWater --WaterResidueNames "HOH,WAT"
      -KeepOldRecords Yes -r Sample2New -o Sample2.pdb

To extract non-hydrogens records from Sample2.pdb file and generate Sample2NonHydrogen.pdb
file, type:

    % ExtractFromPDBFiles.pl -m NonHydrogens Sample2.pdb

To extract data for first chain in Sample2.pdb and generate Sample2ChainA.pdb, type
file, type:

    % ExtractFromPDBFiles.pl -m chains -o Sample2.pdb

To extract data for both chains in Sample2.pdb and generate Sample2ChainA.pdb and
Sample2ChainB.pdb, type:

    % ExtractFromPDBFiles.pl -m chains -c All -o Sample2.pdb

To extract data for alpha carbons in Sample2.pdb and generate Sample2CAlphas.pdb, type:

    % ExtractFromPDBFiles.pl -m CAlphas -o Sample2.pdb

To extract records for specific residue numbers in all chains from Sample2.pdb file and generate
Sample2ResidueNums.pdb file, type:

    % ExtractFromPDBFiles.pl -m ResidueNums --Residues "3,6"
      Sample2.pdb

To extract records for a specific range of residue number in all chains from Sample2.pdb
file and generate Sample2ResiduesRange.pdb file, type:

    % ExtractFromPDBFiles.pl -m ResiduesRange --Residues "10,30"
      Sample2.pdb

To extract data for all ATOM and HETATM records with in 10 angstrom of an atom specifed by
atom serial number and name "1,N" in Sample2.pdb file and generate Sample2DistanceByAtom.pdb,
type:

    % ExtractFromPDBFiles.pl -m Distance --DistanceMode Atom
      --DistanceOrigin "1,N" -k No --distance 10 -o Sample2.pdb

To extract data for all ATOM and HETATM records for complete residues with any atom or hetatm
less than 10 angstrom of an atom specifed by atom serial number and name "1,N" in Sample2.pdb
file and generate Sample2DistanceByAtom.pdb, type:

    % ExtractFromPDBFiles.pl -m Distance --DistanceMode Atom
      --DistanceOrigin "1,N" --DistanceSelectionMode ByResidue
      -k No --distance 10 -o Sample2.pdb

To extract data for all ATOM and HETATM records with in 25 angstrom of an arbitrary point "0,0,0"
in Sample2.pdb file and generate Sample2DistanceByXYZ.pdb, type:

    % ExtractFromPDBFiles.pl -m Distance --DistanceMode XYZ
      --DistanceOrigin "0,0,0" -k No --distance 25 -o Sample2.pdb

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoPDBFiles.pl, ModifyPDBFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
