#!/usr/bin/perl -w
#
# File: ModifyPDBFiles.pl
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
    ModifyPDBFiles($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Modify appropriate information...
sub ModifyPDBFiles {
  my($FileIndex) = @_;
  my($PDBFile, $PDBRecordLinesRef);

  # Get PDB data...
  $PDBFile = $PDBFilesList[$FileIndex];
  $PDBRecordLinesRef = ReadPDBFile($PDBFile);

  if ($OptionsInfo{Mode} =~ /^RenumberAtoms$/i) {
    RenumberAtoms($FileIndex, $PDBRecordLinesRef);
  }
  elsif ($OptionsInfo{Mode} =~ /^RenumberResidues$/i) {
    RenumberResidues($FileIndex, $PDBRecordLinesRef);
  }
  elsif ($OptionsInfo{Mode} =~ /^RenumberWaters$/i) {
    RenumberWaters($FileIndex, $PDBRecordLinesRef);
  }
  elsif ($OptionsInfo{Mode} =~ /^RenameChainIDs$/i) {
    RenameChainsIDs($FileIndex, $PDBRecordLinesRef);
  }
}

# Renumber atom and hetro atom numbers...
sub RenumberAtoms {
  my($FileIndex, $PDBRecordLinesRef) = @_;
  my($PDBFileName,  $RecordLine, $ConectRecordLinesRef, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge, $NewAtomNumber, $RecordType, %OldToNewAtomNumbersMap);

  $PDBFileName = $PDBFilesInfo{OutFileNames}[$FileIndex][0];
  print "Generating PDBFileName file $PDBFileName...\n";
  open OUTFILE, ">$PDBFileName" or die "Error: Can't open $PDBFileName: $! \n";

  # Write out header and other older recors...
  WriteHeaderAndOlderRecords(\*OUTFILE, $PDBRecordLinesRef);

  # Write out all ATOM records along with TER and model records to indicate
  # chains and multiple models..
  %OldToNewAtomNumbersMap = ();
  $NewAtomNumber = $OptionsInfo{StartingAtomNumber};
  for $RecordLine (@{$PDBRecordLinesRef}) {
    if (IsAtomRecordType($RecordLine) || IsHetatmRecordType($RecordLine)) {
      $RecordType = GetPDBRecordType($RecordLine);

      ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = ParseAtomOrHetatmRecordLine($RecordLine);

      print OUTFILE GenerateAtomOrHetatmRecordLine($RecordType, $NewAtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge), "\n";

      $OldToNewAtomNumbersMap{$AtomNumber} = $NewAtomNumber;
      $NewAtomNumber++;
    }
    elsif (IsTerRecordType($RecordLine)) {
      $NewAtomNumber++;
      print OUTFILE GenerateTerRecordLine($NewAtomNumber, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode), "\n";
    }
    elsif (IsModelRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
    elsif (IsEndmdlRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
      # Restart numbering...
      $NewAtomNumber = $OptionsInfo{StartingAtomNumber};
    }
  }

  # Write out modified CONECT records...
  my($ModifiedConectAtomNum, $ConectAtomNum, @ConectAtomNums, @ModifiedConectAtomNums);
  LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
    if (!IsConectRecordType($RecordLine)) {
      next LINE;
    }
    @ConectAtomNums = ();
    @ModifiedConectAtomNums = ();
    push @ConectAtomNums, ParseConectRecordLine($RecordLine);
    ATOMNUMBER: for $ConectAtomNum (@ConectAtomNums) {
      $ModifiedConectAtomNum = $ConectAtomNum;
      if (defined($ConectAtomNum)) {
	$AtomNumber = $ConectAtomNum;
	if ($AtomNumber) {
	  if (exists $OldToNewAtomNumbersMap{$AtomNumber}) {
	    $ModifiedConectAtomNum = $OldToNewAtomNumbersMap{$AtomNumber};
	  }
	}
      }
      push @ModifiedConectAtomNums, $ModifiedConectAtomNum;
    }
    # Write out the record...
    print OUTFILE GenerateConectRecordLine(@ModifiedConectAtomNums), "\n";
  }

  # Write out END record...
  print OUTFILE GenerateEndRecordLine(), "\n";

  close OUTFILE;
}

# Renumber residues...
sub RenumberResidues {
  my($FileIndex, $PDBRecordLinesRef) = @_;
  my($PDBFileName,  $RecordLine, $ConectRecordLinesRef, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge, $NewResidueNumber, $NewHetatmResidueNumber, $TERCount, $TotalTERCount, $PreviousResidueNumber, $PreviousHetatmResidueNumber, $RecordType);

  $PDBFileName = $PDBFilesInfo{OutFileNames}[$FileIndex][0];
  print "Generating PDBFileName file $PDBFileName...\n";
  open OUTFILE, ">$PDBFileName" or die "Error: Can't open $PDBFileName: $! \n";

  # Write out header and other older recors...
  WriteHeaderAndOlderRecords(\*OUTFILE, $PDBRecordLinesRef);

  # Do a quick count of all TER records...
  $TotalTERCount = 0;
  for $RecordLine (@{$PDBRecordLinesRef}) {
    if (IsTerRecordType($RecordLine)) {
      $TotalTERCount++;
    }
  }

  # Write out all ATOM records along with TER and model records to indicate
  # chains and multiple models..
  $NewResidueNumber = $OptionsInfo{StartingResidueNumber};
  $NewHetatmResidueNumber = $OptionsInfo{StartingHetatmResidueNumber};

  $TERCount = 0;
  $PreviousResidueNumber = 0;
  $PreviousHetatmResidueNumber = 0;

  for $RecordLine (@{$PDBRecordLinesRef}) {
    if (IsAtomRecordType($RecordLine) || (IsHetatmRecordType($RecordLine) && ($TERCount < $TotalTERCount || $OptionsInfo{HetatmResidueNumberMode} =~ /^Automatic$/i))) {
      $RecordType = GetPDBRecordType($RecordLine);
      ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = ParseAtomOrHetatmRecordLine($RecordLine);

      if ($PreviousResidueNumber && $PreviousResidueNumber != $ResidueNumber) {
	$PreviousResidueNumber = $ResidueNumber;
	$NewResidueNumber++;
      }
      else {
	# First residue in a chain...
	$PreviousResidueNumber = $ResidueNumber;
      }
      print OUTFILE GenerateAtomOrHetatmRecordLine($RecordType, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $NewResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge), "\n";

    }
    elsif (IsHetatmRecordType($RecordLine)) {
      ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = ParseHetatmRecordLine($RecordLine);

      # User HETATM residue numbers...
      if ($PreviousHetatmResidueNumber && $PreviousHetatmResidueNumber != $ResidueNumber) {
	$PreviousHetatmResidueNumber = $ResidueNumber;
	$NewHetatmResidueNumber++;
      }
      else {
	# First HETATM residue outside a chain...
	$PreviousHetatmResidueNumber = $ResidueNumber;
      }

      print OUTFILE GenerateHetatmRecordLine($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $NewHetatmResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge), "\n";
    }
    elsif (IsTerRecordType($RecordLine)) {
      $TERCount++;
      $AtomNumber++;
      print OUTFILE GenerateTerRecordLine($AtomNumber, $ResidueName, $ChainID, $NewResidueNumber, $InsertionCode), "\n";
      # For per chain numbering, start over again...
      if ($OptionsInfo{ResidueNumberMode} =~ /^PerChain$/i) {
	if ($TERCount < $TotalTERCount ) {
	  $NewResidueNumber = $OptionsInfo{StartingResidueNumber};
	}
	$PreviousResidueNumber = 0;
      }
    }
    elsif (IsModelRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
    elsif (IsEndmdlRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
  }

  # Write out CONECT records...
  for $RecordLine (@{$PDBRecordLinesRef}) {
    if (IsConectRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
  }

  # Write out END record...
  print OUTFILE GenerateEndRecordLine(), "\n";

  close OUTFILE;
}

# Renumber water residues...
sub RenumberWaters {
  my($FileIndex, $PDBRecordLinesRef) = @_;
  my($PDBFileName,  $RecordLine, $ConectRecordLinesRef, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge, $NewResidueNumber, $RecordType);

  $PDBFileName = $PDBFilesInfo{OutFileNames}[$FileIndex][0];
  print "Generating PDBFileName file $PDBFileName...\n";
  open OUTFILE, ">$PDBFileName" or die "Error: Can't open $PDBFileName: $! \n";

  # Write out header and other older recors...
  WriteHeaderAndOlderRecords(\*OUTFILE, $PDBRecordLinesRef);

  # Write out all ATOM records along with TER and model records to indicate
  # chains and multiple models..
  $NewResidueNumber = $OptionsInfo{StartingWaterResidueNumber};
  for $RecordLine (@{$PDBRecordLinesRef}) {
    if (IsAtomRecordType($RecordLine) || IsHetatmRecordType($RecordLine)) {
      $RecordType = GetPDBRecordType($RecordLine);

      ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = ParseAtomOrHetatmRecordLine($RecordLine);

      if (exists $OptionsInfo{SpecifiedWaterResiduesMap}{$ResidueName}) {
	$ResidueNumber = $NewResidueNumber;
	print OUTFILE GenerateAtomOrHetatmRecordLine($RecordType, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge), "\n";
	$NewResidueNumber++;
      }
      else {
	print OUTFILE "$RecordLine\n";
      }
    }
    elsif (IsTerRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
    elsif (IsModelRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
    elsif (IsEndmdlRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
  }

  # Write out CONECT records...
  for $RecordLine (@{$PDBRecordLinesRef}) {
    if (IsConectRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
  }

  # Write out END record...
  print OUTFILE GenerateEndRecordLine(), "\n";

  close OUTFILE;
}

# Rename chain IDs...
sub RenameChainsIDs {
  my($FileIndex, $PDBRecordLinesRef) = @_;
  my($PDBFileName,  $RecordLine, $ConectRecordLinesRef, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge, $RecordType, $PreviousChainID, $FirstChainID, $NewChainID, $NewChainIDCounter, %OldToNewChainIDsMap);

  $PDBFileName = $PDBFilesInfo{OutFileNames}[$FileIndex][0];
  print "Generating PDBFileName file $PDBFileName...\n";
  open OUTFILE, ">$PDBFileName" or die "Error: Can't open $PDBFileName: $! \n";

  # Write out header and other older recors...
  WriteHeaderAndOlderRecords(\*OUTFILE, $PDBRecordLinesRef);

  # Write out all ATOM records along with TER and model records to indicate
  # chains and multiple models..
  %OldToNewChainIDsMap = ();
  $NewChainIDCounter = $OptionsInfo{StartingChainID};
  $FirstChainID = 1;
  $PreviousChainID = '';
  LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
    if (IsAtomRecordType($RecordLine) || IsHetatmRecordType($RecordLine)) {
      $RecordType = GetPDBRecordType($RecordLine);

      ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = ParseAtomOrHetatmRecordLine($RecordLine);

      if (exists $OptionsInfo{SpecifiedWaterResiduesMap}{$ResidueName}) {
	# Chain IDs are not assigned to water residues...
	print OUTFILE "$RecordLine\n";
	next LINE;
      }

      if ($FirstChainID) {
	$FirstChainID = 0;
	$PreviousChainID = $ChainID;
	if ($ChainID || (!$ChainID && $OptionsInfo{RenameEmptyChainIDs})) {
	  $NewChainID = $NewChainIDCounter;
	  $OldToNewChainIDsMap{$ChainID} = $NewChainID;
	}
	else {
	  $NewChainID = '';
	}
      }
      elsif ($PreviousChainID ne $ChainID) {
	if ($ChainID || (!$ChainID && $OptionsInfo{RenameEmptyChainIDs})) {
	  $PreviousChainID = $ChainID;
	  if (exists $OldToNewChainIDsMap{$ChainID}) {
	    $NewChainID = $OldToNewChainIDsMap{$ChainID};
	  }
	  else {
	    $NewChainIDCounter++;
	    $NewChainID = $NewChainIDCounter;
	    $OldToNewChainIDsMap{$ChainID} = $NewChainID;
	  }
	}
	else {
	  $NewChainID = '';
	}
      }

      print OUTFILE GenerateAtomOrHetatmRecordLine($RecordType, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $NewChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge), "\n";
    }
    elsif (IsTerRecordType($RecordLine)) {
      $AtomNumber++;
      print OUTFILE GenerateTerRecordLine($AtomNumber, $ResidueName, $NewChainID, $ResidueNumber, $InsertionCode), "\n";
    }
    elsif (IsModelRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
    elsif (IsEndmdlRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
  }

  # Write out CONECT records...
  for $RecordLine (@{$PDBRecordLinesRef}) {
    if (IsConectRecordType($RecordLine)) {
      print OUTFILE "$RecordLine\n";
    }
  }

  # Write out END record...
  print OUTFILE GenerateEndRecordLine(), "\n";

  close OUTFILE;
}


# Write out modifed header and other older records...
sub WriteHeaderAndOlderRecords {
  my($OutFileRef, $PDBRecordLinesRef) = @_;

  if ($OptionsInfo{ModifyHeaderRecord}) {
    # Write out modified HEADER record...
    my($Classification, $DepositionDate, $IDCode) = GetHeaderRecordInformation($PDBRecordLinesRef);
    $Classification = 'Data modified using MayaChemTools';
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


# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();
  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{StartingAtomNumber} = $Options{atomnumberstart};
  $OptionsInfo{StartingChainID} = $Options{chainidstart};
  $OptionsInfo{RenameEmptyChainIDs} = ($Options{chainidrenameempty} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{KeepOldRecords} = ($Options{keepoldrecords} =~ /^Yes$/i) ? 1 : 0;
  $OptionsInfo{ModifyHeaderRecord} = ($Options{modifyheader} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{ResidueNumberMode} = $Options{residuenumbermode};
  $OptionsInfo{StartingResidueNumber} = $Options{residuenumberstart};

  $OptionsInfo{HetatmResidueNumberMode} = $Options{residuenumberhetatmmode};
  $OptionsInfo{StartingHetatmResidueNumber} = $Options{residuenumberstarthetatm};

  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? 1 : 0;
  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : 0;

  $OptionsInfo{WaterResidueNames} = $Options{waterresiduenames};
  $OptionsInfo{StartingWaterResidueNumber} = $Options{waterresiduestart};
  @{$OptionsInfo{SpecifiedWaterResiduesList}} = ();
  %{$OptionsInfo{SpecifiedWaterResiduesMap}} = ();

  my(@SpecifiedWaterResiduesList);
  @SpecifiedWaterResiduesList = ();
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
  push @{$OptionsInfo{SpecifiedWaterResiduesList}}, @SpecifiedWaterResiduesList;
}

# Retrieve information about PDB files...
sub RetrievePDBFilesInfo {
  my($Index, $PDBFile, $PDBRecordLinesRef, $ChainsAndResiduesInfoRef, $FileDir, $FileName, $FileExt, $OutFileName, $OutFileRoot,  $Mode, $OutFileMode, @OutFileNames);

  %PDBFilesInfo = ();
  @{$PDBFilesInfo{FileOkay}} = ();
  @{$PDBFilesInfo{OutFileRoot}} = ();
  @{$PDBFilesInfo{OutFileNames}} = ();

  FILELIST: for $Index (0 .. $#PDBFilesList) {
    $PDBFilesInfo{FileOkay}[$Index] = 0;

    $PDBFilesInfo{OutFileRoot}[$Index] = '';
    @{$PDBFilesInfo{OutFileNames}[$Index]} = ();
    @{$PDBFilesInfo{OutFileNames}[$Index]} = ();

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
    $ChainsAndResiduesInfoRef = GetChainsAndResidues($PDBRecordLinesRef);
    if (!scalar @{$ChainsAndResiduesInfoRef->{ChainIDs}}) {
      warn "Warning: Ignoring file $PDBFile: No chains found \n";
      next FILELIST;
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
    MODE: {
	if ($Mode =~ /^RenumberAtoms$/i) { $OutFileMode = 'RenumberAtoms'; last MODE;}
	if ($Mode =~ /^RenumberResidues$/i) { $OutFileMode = 'RenumberResidues'; last MODE;}
	if ($Mode =~ /^RenumberWaters$/i) { $OutFileMode = 'RenumberWaters'; last MODE;}
	if ($Mode =~ /^RenameChainIDs$/i) { $OutFileMode = 'RenameChainIDs'; last MODE;}
	$OutFileMode = '';
    }
    $OutFileName = "${OutFileRoot}${OutFileMode}.pdb";
    push @OutFileNames, $OutFileName;

    $PDBFilesInfo{FileOkay}[$Index] = 1;
    $PDBFilesInfo{OutFileRoot}[$Index] = $OutFileRoot;

    push @{$PDBFilesInfo{OutFileNames}[$Index]}, @OutFileNames;
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{atomnumberstart} = 1;
  $Options{chainidstart} = 'A';
  $Options{chainidrenameempty} = 'No';
  $Options{keepoldrecords} = 'no';
  $Options{mode} = 'RenumberResidues';
  $Options{modifyheader} = 'yes';
  $Options{residuenumbermode} = 'PerChain';
  $Options{residuenumberstart} = 1;
  $Options{residuenumberhetatmmode} = 'Automatic';
  $Options{residuenumberstarthetatm} = 6000;
  $Options{waterresiduenames} = 'Automatic';
  $Options{waterresiduestart} = 8000;

  if (!GetOptions(\%Options, "help|h", "atomnumberstart|a=i", "chainidstart|c=s", "chainidrenameempty=s", "keepoldrecords|k=s", "mode|m=s", "modifyheader=s", "overwrite|o", "residuenumbermode=s", "residuenumberstart=i", "residuenumberhetatmmode=s", "residuenumberstarthetatm=i", "root|r=s", "sequencelength=i", "waterresiduenames=s", "waterresiduestart=i", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if (!IsPositiveInteger($Options{atomnumberstart})) {
    die "Error: The value specified, $Options{atomnumberstart}, for option \"-a, --AtomNumberStart\" is not valid. Allowed values: >0\n";
  }
  if ((length($Options{chainidstart}) > 1) || ($Options{chainidstart} !~ /[A-Z]/i)) {
    die "Error: The value specified, $Options{chainidstart}, for option \"-c, --ChainIDStart\" is not valid. Allowed values: a single character from A to Z\n";
  }
  if ($Options{chainidrenameempty} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{chainidrenameempty}, for option \"--chainidrenameempty\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{keepoldrecords} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{keepoldrecords}, for option \"--KeepOldRecords\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{mode} !~ /^(RenumberAtoms|RenumberResidues|RenumberWaters|RenameChainIDs)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: RenumberAtoms, RenumberResidues, RenumberWaters or RenameChainIDs\n";
  }
  if ($Options{modifyheader} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{modifyheader}, for option \"--ModifyHeader\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{residuenumbermode} !~ /^(Sequential|PerChain)$/i) {
    die "Error: The value specified, $Options{residuenumbermode}, for option \"--ResidueNumberMode\" is not valid. Allowed values: Sequential or PerChain\n";
  }
  if (!IsPositiveInteger($Options{residuenumberstart})) {
    die "Error: The value specified, $Options{residuenumberstart}, for option \"--ResidueNumberStart\" is not valid. Allowed values: >0\n";
  }
  if ($Options{residuenumberhetatmmode} !~ /^(automatic|specify)$/i) {
    die "Error: The value specified, $Options{residuenumberhetatmmode}, for option \"--residuenumbermode\" is not valid. Allowed values: automatic or specify\n";
  }
  if (!IsPositiveInteger($Options{residuenumberstarthetatm})) {
    die "Error: The value specified, $Options{residuenumberstarthetatm}, for option \"--residuenumberstartHetatm\" is not valid. Allowed values: >0\n";
  }
  if (!IsPositiveInteger $Options{waterresiduestart}) {
    die "Error: The value specified, $Options{waterresiduestart}, for option \"--waterresiduestart\" is not valid. Allowed values: >0\n";
  }
}

__END__

=head1 NAME

ModifyPDBFiles.pl - Modify data in PDBFile(s)

=head1 SYNOPSIS

ModifyPDBFiles.pl PDBFile(s)...

ModifyPDBFiles.pl [B<-a, --AtomNumberStart> number] [B<-c, --ChainIDStart> character]
[B<--ChainIDRenameEmpty> yes | no] [B<-h, --help>] [B<-k, --KeepOldRecords> yes | no]
[B<-m, --mode > RenumberAtoms | RenumberResidues | RenumberWaters | RenameChainIDs]
[B<--ModifyHeader> yes | no] [B<-o, --overwrite>] [B<--ResidueNumberMode> Sequential | PerChain]
[B<--ResidueNumberStart> number] [B<--ResidueNumberHetatmMode> automatic | specify]
[B<--ResidueNumberStarHetatm> number] [B<-r, --root> rootname]
[B<--WaterResidueNames> Automatic | "ResidueName, [ResidueName,...]"] [B<--WaterResidueStart> number]
[B<-w, --WorkingDir> dirname] PDBFile(s)...

=head1 DESCRIPTION

Modify data in I<PDBFile(s)>: renumber atoms, residues, and water residues or assign new
chain IDs. Multiple PDBFile names are separated by spaces. The valid file extension is I<.pdb>.
All other file name extensions are ignored during the wild card expansion. All the PDB files
in a current directory can be specified either by I<*.pdb> or the current directory name.

=head1 OPTIONS

=over 4

=item B<-a, --AtomNumberStart> I<number>

Starting atom number to use during I<RenumberAtoms> value of B<-m, --mode> option. Default: I<1>.
Valid values: positive integers.

=item B<-c, --ChainIDStart> I<character>

A single character to use for starting IDs for chains during I<RenameChainIDs> value of B<-m, --mode> option.
Default: I<A>. Valid values: I<A to Z>.

=item B<--ChainIDRenameEmpty> I<Yes | No>

Specify whether to rename empty chain IDs during I<RenameChainIDs> B<-m, --mode> value. By
default, ATOM and HETATM records with no chain IDs are left unchanged. Possible values:
I<yes | no>. Default: I<No>.

=item B<-h, --help>

Print this help message.

=item B<-k, --KeepOldRecords> I<yes | no>

Specify whether to transfer old non ATOM and HETATM records from input PDBFile(s) to new
PDBFile(s). By default, except for the HEADER record, all records other than ATOM/HETATM
are dropped during the generation of new PDB files. Possible values: I<yes | no>.
Default: I<no>.

=item B<-m, --mode > I<RenumberAtoms | RenumberResidues | RenumberWaters | RenameChainIDs>

Specify how to modify I<PDBFile(s)>. Possible values: I<RenumberAtoms | RenumberResidues
| RenumberWaters | RenameChainIDs>. Default: I<RenumberResidues>.

For I<RenumberAtoms> mode, residue number in ATOM and HETATM records are reassigned
sequentially starting using value of B<-a, --AtomNumberStart> option.

For I<RenumberResidues> mode, serial number in ATOM and HETATM records are reassigned
either sequentially or statring from specified values for ATOM and HETATM records in each
chain.

For I<RenumberWaters> mode, residue number for waters are reassigned starting from a specific
value.

For I<RenameChainIDs> mode, all the chain IDs are reassigned starting from a specific chain ID.

During the generation of new PDB files, unnecessary CONECT records are dropped.

=item B<--ModifyHeader> I<yes | no>

Specify whether to modify HEADER record during the generation of new PDB files
Possible values: I<yes | no>.  Default: I<yes>. By defailt, Classification data is replaced
by I<Data modified using MayaChemTools> before writing out HEADER record.

=item B<-o, --overwrite>

Overwrite existing files

=item B<--ResidueNumberMode> I<Sequential | PerChain>

Specify how to renumber residues: renumber residues sequentially across all the chains
or start from the begining for each chain. Possible values: I<Sequential | PerChain>. Default:
I<PerChain>.

=item B<--ResidueNumberStart> I<number>

Starting residue number to use for ATOM records in chains. Default: I<1>. Valid values
positive integers.

For I<Sequential> value of B<--ResidueNumberMode> option, residue numbers are
assigned sequentially across all the chains starting from the specified value.

For I<PerChain> value of B<--ResidueNumberMode> option, residue numbers are
starting again from the specified value for each chain.

HETATM residues with in the chains are numbered using this value as well

=item B<--ResidueNumberHetatmMode> I<automatic | specify>

Specify how to start residue number for HETATM records: use the next sequential
residue number after the last residue number from ATOM records or start from a
specific residue number. Possible values: I<automatic | specify>. Default:
I<automatic>

For I<automatic> , residue number after highest residue number of ATOM
records is used as the starting residue number for HETATM records.

For I<specify>,  value of option B<--ResidueNumberStarHetatm> is used as the
starting residue number for HETATM records.

This option along with B<--ResidueNumberStartHetatm> only applies to HETATM records
outside the chains.

=item B<--ResidueNumberStartHetatm> I<number>

Starting residue number to use for HETATM records. Default: I<6000>. Valid values
positive integers.

=item B<-r, --root> I<rootname>

New PDB and sequence file name is generated using the root: <Root><Mode>.<Ext>.
Default new file name: <PDBFileName><Mode>.pdb. This option is ignored for multiple
input files.

=item B<--WaterResidueNames> I<Automatic | "ResidueName,[ResidueName,...]">

Identification of water residues during I<RenumberWaters> value of B<-m, --mode> option. Possible
values: I<Automatic | "ResidueName,[ResidueName,...]">. Default: I<Automatic> which corresponds
to "HOH,WAT,H20". You can also specify a different comma delimited list of residue names
to use for water.

=item B<--WaterResidueStart> I<number>

Starting water residue number to use during I<RenumberWaters> B<-m, --mode> value.
Default: I<8000>. Valid values: positive integers.

=item B<-w, --WorkingDir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To renumber ATOM and HETATM residues starting from 1 for each chain with continuation to
HETATM residues outside TER records in Sample2.pdb and generate
Sample2RenumberResidues.pdb file, type:

    % ModifyPDBFiles.pl Sample1.pdb

To renumber ATOM and HETATM residues sequentially across all chains starting from 1 with
continuation to HETATM residues outside TER records in Sample2.pdb and generate
Sample2RenumberResidues.pdb file, type:

    % ModifyPDBFiles.pl --ResidueNumberMode Sequential -o Sample1.pdb

To renumber ATOM and HETATM residues sequentially across all chains starting from 1 and
HETATM residues outside TER records starting from 6000 in Sample2.pdb and generate
Sample2RenumberResidues.pdb file, type:

    % ModifyPDBFiles.pl --ResidueNumberMode Sequential
      --ResidueNumberHetatmMode Specify  -o Sample1.pdb


To renumber ATOM and HETATM residues sequentially across all chains starting from 100 for
ATOM/HETATM  residues with in TER records and starting from 999 for HETATM residues
outside TER records in Sample2.pdb and generate Sample2RenumberResidues.pdb file, type:

    % ModifyPDBFiles.pl --ResidueNumberMode Sequential
      --ResidueNumberHetatmMode Specify --ResidueNumberStart 100
      --ResidueNumberStartHetatm 999 -o Sample2.pdb

To renumber ATOM and HETATM residues from 100 for each chain and starting from 999 for
HETATM  residues outside TER records in Sample2.pdb and generate Sample2RenumberResidues.pdb
file, type:

    % ModifyPDBFiles.pl --ResidueNumberMode PerChain
      --ResidueNumberHetatmMode Specify --ResidueNumberStart 100
      --ResidueNumberStartHetatm 999 -o Sample2.pdb

To renumber ATOM serial numbers sequentially starting from 100 in Sample1.pdb file and generate
Sample1RenumberAtoms.pdb file, type:

    % ModifyPDBFiles.pl -m RenumberAtoms --AtomNumberStart 100
      -o Sample1.pdb

To renumber water residues identified by "HOH,WAT" starting from residue number 1000
in Sample2.pdb file and generate Sample2RenumberWaters.pdb file, type:

    % ModifyPDBFiles.pl -m RenumberWaters --WaterResidueNames "HOH,WAT"
      -o --WaterResidueStart 950 Sample2.pdb

To rename all chain IDs starting from A in Sample1.pdb file and generate
Sample1RenameChainIDs.pdb file, type:

    % ModifyPDBFiles.pl -m RenameChainIDs -o Sample1.pdb

To rename all chain IDs starting from B without assigning any chain IDs to ATOM/HETATOM
with no chain IDs in Sample2.pdb file and generate Sample2RenameChainIDs.pdb file, type:

    % ModifyPDBFiles.pl l -m RenameChainIDs -c B --ChainIDRenameEmpty No
      -o Sample2.pdb


=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

ExtractFromPDBFiles.pl, InfoPDBFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
