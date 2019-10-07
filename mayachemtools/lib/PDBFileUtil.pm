package PDBFileUtil;
#
# File: PDBFileUtil.pm
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
use Exporter;
use Text::ParseWords;
use TextUtil;
use FileUtil;
use TimeUtil ();

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw(GetPDBRecordType GetRecordTypesCount GetAllResidues GetConectRecordLines GetChainsAndResidues GetExperimentalTechnique GetExperimentalTechniqueResolution GetMinMaxCoords IsPDBFile IsAtomRecordType IsConectRecordType IsHeaderRecordType IsHetatmRecordType IsSeqresRecordType IsModelRecordType IsEndmdlRecordType IsTerRecordType IsMasterRecordType ReadPDBFile ParseHeaderRecordLine GenerateHeaderRecordLine GenerateHeaderRecordTimeStamp ParseAtomRecordLine GenerateAtomRecordLine ParseAtomOrHetatmRecordLine GenerateAtomOrHetatmRecordLine GenerateHetatmRecordLine ParseHetatmRecordLine ParseConectRecordLine GenerateConectRecordLine ParseExpdtaRecordLine ParseRemark2ResolutionRecordLine ParseSeqresRecordLine ParseTerRecordLine GenerateTerRecordLine ParseMasterRecordLine GenerateEndRecordLine);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Get PDB record type...
sub GetPDBRecordType {
  my($Line) = @_;

  return _GetRecordType($Line);
}

# Is it a PDB file?
sub IsPDBFile {
  my($PDBFile) = @_;
  my($Line, $Status);

  $Status = 0;
  open PDBFILE, "$PDBFile" or die "Can't open $PDBFile: $!\n";
  $Line = GetTextLine(\*PDBFILE);
  $Status = ($Line =~ /^HEADER/i) ? 1 : 0;
  close PDBFILE;

  return $Status;
}

# Is it a atom record type?
sub IsAtomRecordType {
  my($Line) = @_;

  return _IsRecordType($Line, 'ATOM');
}

# Is it a connect record type?
sub IsConectRecordType {
  my($Line) = @_;

  return _IsRecordType($Line, 'CONECT');
}

# Is it a header atom record type?
sub IsHeaderRecordType {
  my($Line) = @_;

  return _IsRecordType($Line, 'HEADER');
}

# Is it a hetro atom record type?
sub IsHetatmRecordType {
  my($Line) = @_;

  return _IsRecordType($Line, 'HETATM');
}

# Is it a seqres record type?
sub IsSeqresRecordType {
  my($Line) = @_;

  return _IsRecordType($Line, 'SEQRES');
}

# Is it a MODEL record type?
sub IsModelRecordType {
  my($Line) = @_;

  return _IsRecordType($Line, 'MODEL');
}

# Is it a ENDMDL record type?
sub IsEndmdlRecordType {
  my($Line) = @_;

  return _IsRecordType($Line, 'ENDMDL');
}

# Is it a TER record type?
sub IsTerRecordType {
  my($Line) = @_;

  return _IsRecordType($Line, 'TER');
}

# Is it a MASTER record type?
sub IsMasterRecordType {
  my($Line) = @_;

  return _IsRecordType($Line, 'MASTER');
}

# Count the number of each record type and a reference to data type with these key/value pairs:
# {RecordTypes} - An array of unique record types in order of their presence in the file
# {Count}{$RecordType} - Count of each record type
# {Lines}{$RecordType} - Optional lines data for a specific record type.
#
sub GetRecordTypesCount {
  my($PDBRecordLinesRef, $SpecifiedRecordType, $GetRecordLinesFlag, $RecordType, $RecordLine, %RecordTypeDataMap);

  %RecordTypeDataMap = ();
  @{$RecordTypeDataMap{RecordTypes}} = ();
  %{$RecordTypeDataMap{Count}} = ();
  %{$RecordTypeDataMap{Lines}} = ();

  $SpecifiedRecordType = '';
  $GetRecordLinesFlag = 0;
  if (@_ == 3) {
    ($PDBRecordLinesRef, $SpecifiedRecordType, $GetRecordLinesFlag) = @_;
    $SpecifiedRecordType = uc $SpecifiedRecordType;
  }
  elsif (@_ == 2) {
    ($PDBRecordLinesRef, $SpecifiedRecordType) = @_;
    $SpecifiedRecordType = uc $SpecifiedRecordType;
  }
  else {
    ($PDBRecordLinesRef) = @_;
  }
  LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
    $RecordType = _GetRecordType($RecordLine);
    if ($SpecifiedRecordType && ($SpecifiedRecordType ne $RecordType)) {
      next LINE;
    }
    if (exists $RecordTypeDataMap{Count}{$RecordType}) {
      # Update count...
      $RecordTypeDataMap{Count}{$RecordType} += 1;

      if ($GetRecordLinesFlag) {
	push @{$RecordTypeDataMap{Lines}{$RecordType}}, $RecordLine;
      }
    }
    else {
      # New record type...
      push @{$RecordTypeDataMap{RecordTypes}}, $RecordType;
      $RecordTypeDataMap{Count}{$RecordType} = 1;

      if ($GetRecordLinesFlag) {
	@{$RecordTypeDataMap{Lines}{$RecordType}} = ();
	push @{$RecordTypeDataMap{Lines}{$RecordType}}, $RecordLine;
      }
    }
  }
  return (\%RecordTypeDataMap);
}

# Collect CONECT record lines for specific atom number, modified specified data to exclude any atom
# number not present in the list of specified atom numbers and return a reference to list of
# CONECT record lines.
#
sub GetConectRecordLines {
  my($PDBRecordLinesRef, $AtomNumbersMapRef) = @_;
  my($AtomNumber, $ConectAtomNumber, $RecordLine, @ConectRecordAtomNums, @ConectRecordLines);

  @ConectRecordLines = ();
  LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
    if (!IsConectRecordType($RecordLine)) {
      next LINE;
    }
    @ConectRecordAtomNums = ();
    push @ConectRecordAtomNums, ParseConectRecordLine($RecordLine);
    ATOMNUMBER: for $ConectAtomNumber (@ConectRecordAtomNums) {
      if (defined $ConectAtomNumber) {
	$AtomNumber = $ConectAtomNumber;
	if ($AtomNumber) {
	  if (! exists $AtomNumbersMapRef->{$AtomNumber}) {
	    next LINE;
	  }
	}
      }
    }
    push @ConectRecordLines, $RecordLine;
  }
  return \@ConectRecordLines;
}

# Get chains and residue information using ATOM/HETATM or SEQRES records. And return a reference to a
#  hash with these keys:
#
# @{$ChainsDataMap{ChainIDs}} - List of chain IDs with 'None' for no chain identification
# @{$ChainsDataMap{Residues}{$ChainID}} - List of residues in order of their appearance in a chain
# @{$ChainsDataMap{ResidueNumbers}{$ChainID}} - List of residue numbers in order of their appearance in a chain
# %{$ChainsDataMap{ResidueCount}{$ChainID}{$ResidueName}} - Count of specific residues in a chain
#
# Notes:
#  . Chains and residue data can be extacted using either ATOM/HETATM records or SEQRES records.
#  . In addition to a different chain ID in ATOM/HETATM a TER record also indicates end of an existing chain
#    and start of a new one: ChainID in ATOM/HETATM records might still be emtpy.
#   . ATOM/HETATM records after the first ENDMDL records are simply ingnored.
#
sub GetChainsAndResidues {
  my($PDBRecordLinesRef, $RecordsSource, $GetChainResiduesBeyondTERFlag, $GetRecordLinesFlag);

  $RecordsSource = 'AtomAndHetatm';
  $GetChainResiduesBeyondTERFlag = 0;
  $GetRecordLinesFlag = 0;

  if (@_ == 4) {
    ($PDBRecordLinesRef, $RecordsSource, $GetChainResiduesBeyondTERFlag, $GetRecordLinesFlag) = @_;
  }
  elsif (@_ == 3) {
    ($PDBRecordLinesRef, $RecordsSource, $GetChainResiduesBeyondTERFlag) = @_;
  }
  elsif (@_ == 2) {
    ($PDBRecordLinesRef, $RecordsSource) = @_;
  }
  else {
    ($PDBRecordLinesRef) = @_;
  }

  if ($RecordsSource =~ /^AtomAndHetatm$/i) {
    return _GetChainsAndResiduesFromAtomHetatmRecords($PDBRecordLinesRef, $GetChainResiduesBeyondTERFlag, $GetRecordLinesFlag);
  }
  elsif ($RecordsSource =~ /^Seqres$/i) {
    return _GetChainsAndResiduesFromSeqresRecords($PDBRecordLinesRef);
  }
  else {
    my(%ChainsDataMap);
    %ChainsDataMap = ();
    @{$ChainsDataMap{ChainIDs}} = ();
    %{$ChainsDataMap{Residues}} = ();
    %{$ChainsDataMap{ResidueNumbers}} = ();
    %{$ChainsDataMap{ResidueCount}} = ();

    return \%ChainsDataMap;
  }
}


# Get residue information using ATOM/HETATM records and return a reference to a hash with
# these keys:
#
# @{$ResiduesDataMap{ResidueNames}} - List of all the residues
# %{$ResiduesDataMap{ResidueCount}{$ResidueName}} - Count of residues
# @{$ResiduesDataMap{AtomResidueNames}} - List of all the residues
# %{$ResiduesDataMap{AtomResidueCount}{$ResidueName}} - Count of residues in ATOM records
# @{$ResiduesDataMap{HetatomResidueNames}} - List of all the residues
# %{$ResiduesDataMap{HetatmResidueCount}{$ResidueName}} - Count of residues HETATM records
#
# Notes:
#  . ATOM/HETATM records after the first ENDMDL records are simply ingnored.
#
sub GetAllResidues {
  my($PDBRecordLinesRef) = @_;

  my($PreviousChainID, $PreviousResidueNumber, $RecordLine, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge, %ResiduesDataMap);

  %ResiduesDataMap = ();
  @{$ResiduesDataMap{ResidueNames}} = ();
  %{$ResiduesDataMap{ResidueCount}} = ();
  @{$ResiduesDataMap{AtomResidueNames}} = ();
  %{$ResiduesDataMap{AtomResidueCount}} = ();
  @{$ResiduesDataMap{HetatmResidueNames}} = ();
  %{$ResiduesDataMap{HetatmResidueCount}} = ();

  $PreviousChainID = '';
  $PreviousResidueNumber = 0;

  LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
    if (IsEndmdlRecordType($RecordLine)) {
      last LINE;
    }
    if (!(IsAtomRecordType($RecordLine) || IsHetatmRecordType($RecordLine))) {
      next LINE;
    }
    ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = ParseAtomRecordLine($RecordLine);

    if ($PreviousChainID eq $ChainID) {
      if ($ResidueNumber == $PreviousResidueNumber) {
	next LINE;
      }
      $PreviousResidueNumber = $ResidueNumber;
    }
    else {
      # New chain...
      $PreviousChainID = $ChainID;
      $PreviousResidueNumber = $ResidueNumber;
    }

    # Store the residue and update its count...
    push @{$ResiduesDataMap{ResidueNames}}, $ResidueName;
    if (exists $ResiduesDataMap{ResidueCount}{$ResidueName}) {
      $ResiduesDataMap{ResidueCount}{$ResidueName} += 1;
    }
    else {
      $ResiduesDataMap{ResidueCount}{$ResidueName} = 1;
    }
    # Update ATOM residue data...
    if (IsAtomRecordType($RecordLine)) {
      push @{$ResiduesDataMap{AtomResidueNames}}, $ResidueName;
      if (exists $ResiduesDataMap{AtomResidueCount}{$ResidueName}) {
	$ResiduesDataMap{AtomResidueCount}{$ResidueName} += 1;
      }
      else {
	$ResiduesDataMap{AtomResidueCount}{$ResidueName} = 1;
      }
    }
    # Update HETATM residue data...
    if (IsHetatmRecordType($RecordLine)) {
      push @{$ResiduesDataMap{HetatmResidueNames}}, $ResidueName;
      if (exists $ResiduesDataMap{HetatmResidueCount}{$ResidueName}) {
	$ResiduesDataMap{HetatmResidueCount}{$ResidueName} += 1;
      }
      else {
	$ResiduesDataMap{HetatmResidueCount}{$ResidueName} = 1;
      }
    }
  }

  return \%ResiduesDataMap;
}

# Return min/max XYZ coordinates for ATOM/HETATM records...
sub GetMinMaxCoords {
  my($PDBRecordLinesRef) = @_;

  my($XMin, $YMin, $ZMin, $XMax, $YMax, $ZMax, $RecordLine, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge);

  ($XMin, $YMin, $ZMin) = (99999) x 3;
  ($XMax, $YMax, $ZMax) = (-99999) x 3;

  LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
    if (!(IsAtomRecordType($RecordLine) || IsHetatmRecordType($RecordLine))) {
      next LINE;
    }
    ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = ParseAtomRecordLine($RecordLine);

    $XMin = ($X < $XMin) ? $X : $XMin;
    $YMin = ($Y < $YMin) ? $Y : $YMin;
    $ZMin = ($Z < $ZMin) ? $Z : $ZMin;

    $XMax = ($X > $XMax) ? $X : $XMax;
    $YMax = ($Y > $YMax) ? $Y : $YMax;
    $ZMax = ($Z > $ZMax) ? $Z : $ZMax;
  }

  if ($XMin == 99999) { $XMin = undef; }
  if ($YMin == 99999) { $YMin = undef; }
  if ($ZMin == 99999) { $ZMin = undef; }
  if ($XMax == -99999) { $XMax = undef; }
  if ($YMax == -99999) { $YMax = undef; }
  if ($ZMax == -99999) { $ZMax = undef; }

  return ($XMin, $YMin, $ZMin, $XMax, $YMax, $ZMax);
}

# Read PDB file and return reference to record lines..
sub ReadPDBFile {
  my($PDBFile) = @_;

  my($Line, @PDBRecordLines);

  @PDBRecordLines = ();
  open PDBFILE, "$PDBFile" or die "Can't open $PDBFile: $!\n";
  while ($Line = GetTextLine(\*PDBFILE)) {
    push @PDBRecordLines, $Line;
  }

  close PDBFILE;

  return (\@PDBRecordLines);
}

#
# Get experimental technique information...
#
sub GetExperimentalTechnique {
  my($PDBRecordLinesRef) = @_;
  my($RecordLine, $ContinuationNum, $ExperimentalTechnique);

  $ExperimentalTechnique = undef;

  LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
      if (_IsRecordType($RecordLine, 'EXPDTA')) {
	($ContinuationNum, $ExperimentalTechnique) = ParseExpdtaRecordLine($RecordLine);
	last LINE;
      }
  }

  return $ExperimentalTechnique;
}

#
# Get experimental technique resolution information...
#
sub GetExperimentalTechniqueResolution {
  my($PDBRecordLinesRef) = @_;
  my($RecordLine, $Resolution, $ResolutionUnits);

  ($Resolution, $ResolutionUnits) = ((undef) x 2);

  LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
    if ($RecordLine =~ /^REMARK   2 RESOLUTION./i) {
      ($Resolution, $ResolutionUnits) = ParseRemark2ResolutionRecordLine($RecordLine);
      last LINE;
    }
  }

  return ($Resolution, $ResolutionUnits);
}

#
# Parse HEADER record line...
sub ParseHeaderRecordLine {
  my($Line) = @_;
  my($Classification, $DepositionDate, $IDCode) = (undef, undef, undef);

  if ($Line !~ /^HEADER/i) {
    return ($Classification, $DepositionDate, $IDCode);
  }
  my($Length);

  ($Classification, $DepositionDate, $IDCode) = ('') x 3;

  $Length = length $Line;

  if ($Length <= 62) {
    ($Classification, $DepositionDate) = unpack("x10A40A9", $Line);
  }
  else {
    ($Classification, $DepositionDate, $IDCode) = unpack("x10A40A9x3A4", $Line);
  }

  $Classification = RemoveLeadingAndTrailingWhiteSpaces($Classification);
  $DepositionDate =~ s/ //g;
  $IDCode =~ s/ //g;

  return ($Classification, $DepositionDate, $IDCode);
}

#
# Generate HEADER record line...
sub GenerateHeaderRecordLine {
  my($Classification, $Date, $IDCode, $Line);

  $Classification = "Created using MayaChemTools";
  $Date = GenerateHeaderRecordTimeStamp();
  if (@_ == 3) {
    ($IDCode, $Classification, $Date) = @_;
  }
  elsif (@_ == 2) {
    ($IDCode, $Classification) = @_;
  }
  elsif (@_ == 1) {
    ($IDCode) = @_;
  }

  $Line = sprintf "HEADER    %-40.40s%9.9s   %4.4s", $Classification, $Date, $IDCode;
  return $Line;
}

# Generate PDB header time stamp...
sub GenerateHeaderRecordTimeStamp {
  return TimeUtil::PDBFileTimeStamp();
}

#
# Parse ATOM record line.
#
sub ParseAtomRecordLine {
  my($Line) = @_;

  return _ParseAtomOrHetatmRecordLine($Line);
}

# Generate ATOM record line...
sub GenerateAtomRecordLine {
  my($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = @_;

  return _GenerateAtomOrHetatmRecordLine('ATOM', $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge);
}

#
# Parse ATOM/HETATm record line.
#
sub ParseAtomOrHetatmRecordLine {
  my($Line) = @_;

  return _ParseAtomOrHetatmRecordLine($Line);
}

# Generate ATOM/HETATM record line...
sub GenerateAtomOrHetatmRecordLine {
  my($RecordType, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = @_;

  return _GenerateAtomOrHetatmRecordLine($RecordType, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge);
}
#
# Parse HETATM record line...
#
sub ParseHetatmRecordLine {
  my($Line) = @_;

  return _ParseAtomOrHetatmRecordLine($Line);
}

# Generate HETATM record line...
sub GenerateHetatmRecordLine {
  my($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = @_;

  return _GenerateAtomOrHetatmRecordLine('HETATM', $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge);
}

# Parse EXPDTA record line...
#
# EXPDTA format:
#
#1 - 6 Record name "EXPDTA"
# 9 - 10 Continuation continuation Allows concatenation of multiple records.
# 11 - 70 SList technique The experimental technique(s) with optional comment describing the sample or experiment.
#
# The EXPDTA record identifies the experimental technique used. This may refer to the type of radiation and
# sample, or include the spectroscopic or modeling technique. Permitted values include:
#
# ELECTRON DIFFRACTION
# FIBER DIFFRACTION
# FLUORESCENCE TRANSFER
# NEUTRON DIFFRACTION
# NMR
# THEORETICAL MODEL
# X-RAY DIFFRACTION
#
sub ParseExpdtaRecordLine {
  my($Line) = @_;

  if ($Line !~ /^EXPDTA/i) {
    return ((undef) x 2);
  }

  my($ContinuationNum, $ExperimentalTechnique) = unpack("x8A2A60" , $Line);

  $ContinuationNum =~ s/ //g;
  $ExperimentalTechnique = RemoveLeadingAndTrailingWhiteSpaces($ExperimentalTechnique);

  return ($ContinuationNum, $ExperimentalTechnique);
}

# Parse REMARK 2 record line...
#
# REMARK 2 format:
#
# The second REMARK 2 record has one of two formats. The first is used for diffraction studies, the second
# for other types of experiments in which resolution is not relevant, e.g., NMR and theoretical modeling.
#
#For diffraction experiments:
#
# 1 - 6 Record name "REMARK"
# 10 LString(1) "2"
# 12 - 22 LString(11) "RESOLUTION."
# 23 - 27 Real(5.2) resolution Resolution.
# 29 - 38 LString(10) "ANGSTROMS."
#
# REMARK 2 when not a diffraction experiment:
#
# 1 - 6 Record name "REMARK"
# 10 LString(1) "2"
# 12 - 38 LString(28) "RESOLUTION. NOT APPLICABLE."
# 41 - 70 String comment Comment.
#
sub ParseRemark2ResolutionRecordLine {
  my($Line) = @_;

  if ($Line !~ /^REMARK   2 RESOLUTION./i) {
    return ((undef) x 2);
  }

  my($Resolution, $ResolutionUnits);

  if ($Line =~ /NOT APPLICABLE/i) {
    ($Resolution, $ResolutionUnits) = ("NOT APPLICABLE", "");
  }
  else {
    ($Resolution, $ResolutionUnits) = unpack("x22A5x1A10" , $Line);
  }

  $Resolution = RemoveLeadingAndTrailingWhiteSpaces($Resolution);

  $ResolutionUnits = RemoveLeadingAndTrailingWhiteSpaces($ResolutionUnits);
  $ResolutionUnits =~ s/\.$//;

  return ($Resolution, $ResolutionUnits);
}

#
# Parse SEQRES record line...
#
# SEQRES format:
#
# 1 - 6 Record name "SEQRES"
# 9 - 10 Serial number of the SEQRES record for the current chain. Starts at 1 and increments by one each line. Reset to 1 for each chain.
# 12 - Chain identifier
# 14 - 17 Integer numRes Number of residues in the chain
# 20 - 22 24 -26 ... ... 68 - 70 Residue name resName Residue name.
#
sub ParseSeqresRecordLine {
  my($Line) = @_;

  if ($Line !~ /^SEQRES/i) {
    return ((undef) x 5);
  }
  my($RecordSerialNumber, $ChainID, $NumOfResidues, $ResidueNames) = unpack("x8A2x1A1x1A4x2A51" , $Line);
  $RecordSerialNumber =~ s/ //g;
  $ChainID =~ s/ //g;
  $NumOfResidues =~ s/ //g;
  $ResidueNames = RemoveLeadingAndTrailingWhiteSpaces($ResidueNames);

  return ($RecordSerialNumber, $ChainID, $NumOfResidues, $ResidueNames);
}

#
# Parse CONECT record line...
#
# CONECT format:
#
# 1 - 6 Record name "CONECT"
# 7 - 11 Atom number
# 12 - 16, 17 - 21, 22 - 26, 27 - 31 Atom number of bonded atom
#
# 32 - 36, 37 - 41 Atom number of hydrogen bonded atom
# 42 - 46 Atom number of salt bridged atom
# 47 - 51, 52 -56 Atom number of hydrogen bonded atom
# 57 - 61 Atom number of salt bridged atom
#
sub ParseConectRecordLine {
  my($Line) = @_;

  if ($Line !~ /^CONECT/i) {
    return ((undef) x 11);
  }
  my($AtomNum, $BondedAtomNum1, $BondedAtomNum2, $BondedAtomNum3, $BondedAtomNum4, $HBondedAtomNum1, $HBondedAtomNum2, $SaltBridgedAtomNum1, $HBondedAtomNum3, $HBondedAtomNum4, $SaltBridgedAtomNum2) = map {s/ //g; $_} unpack("x6A5A5A5A5A5A5A5A5A5A5A5", $Line);

  return ($AtomNum, $BondedAtomNum1, $BondedAtomNum2, $BondedAtomNum3, $BondedAtomNum4, $HBondedAtomNum1, $HBondedAtomNum2, $SaltBridgedAtomNum1, $HBondedAtomNum3, $HBondedAtomNum4, $SaltBridgedAtomNum2);
}

# Generate CONECT record line...
sub GenerateConectRecordLine {
  my($AtomNum, $BondedAtomNum1, $BondedAtomNum2, $BondedAtomNum3, $BondedAtomNum4, $HBondedAtomNum1, $HBondedAtomNum2, $SaltBridgedAtomNum1, $HBondedAtomNum3, $HBondedAtomNum4, $SaltBridgedAtomNum2) = @_;
  my($Line);

  $Line = sprintf "CONECT%5.5s%5.5s%5.5s%5.5s%5.5s%5.5s%5.5s%5.5s%5.5s%5.5s%5.5s", $AtomNum, $BondedAtomNum1, $BondedAtomNum2, $BondedAtomNum3, $BondedAtomNum4, $HBondedAtomNum1, $HBondedAtomNum2, $SaltBridgedAtomNum1, $HBondedAtomNum3, $HBondedAtomNum4, $SaltBridgedAtomNum2;

  return $Line;
}

#
# Parse TER record line...
#
# TER format:
#
#1 - 6 Record name "TER "
# 7 - 11 Serial number
# 18 - 20 Residue name
# 22 Chain identifier
# 23 - 26 Residue sequence number
# 27 Insertion code
#
sub ParseTerRecordLine {
  my($Line) = @_;

  if ($Line !~ /^TER/i) {
    return ((undef) x 5);
  }
  my($SerialNumber, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $Length);

  ($SerialNumber, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode) = ('') x 5;

  $Length = length $Line;

  if ($Length <= 17) {
    ($SerialNumber, $ResidueName) = map {s/ //g; $_} unpack("x6A5", $Line);
  }
  elsif ($Length <= 21) {
    ($SerialNumber, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode) = map {s/ //g; $_} unpack("x6A5x6A3", $Line);
  }
  else {
    ($SerialNumber, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode) = map {s/ //g; $_} unpack("x6A5x6A3xA1A4A1", $Line);
  }

  return ($SerialNumber, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode);
}

# Generate TER record line...
sub GenerateTerRecordLine {
  my($SerialNumber, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $Line) = ('') x 6;

  if (@_ == 5) {
    ($SerialNumber, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode) = @_;
  }
  elsif (@_ == 4) {
    ($SerialNumber, $ResidueName, $ChainID, $ResidueNumber) = @_;
  }
  elsif (@_ == 3) {
    ($SerialNumber, $ResidueName) = @_;
  }
  elsif (@_ == 2) {
    ($SerialNumber, $ResidueName) = @_;
  }
  elsif (@_ == 1) {
    ($SerialNumber) = @_;
  }
  $Line = sprintf "TER   %5.5s      %-3.3s %1.1s%4.4s%1.1s", $SerialNumber, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode;

  return $Line;
}

#
# Parse MASTER record line...
#
# MASTER record format:
#
#1 - 6 Record name "MASTER"
# 11 - 15 Number of REMARK records
# 16 - 20 "0"
# 21 - 25 Number of HET records
# 26 - 30 Number of HELIX records
# 31 - 35 Number of SHEET records
# 36 - 40 Number of TURN records
# 41 - 45 Number of SITE records
# 46 - 50 Number of coordinate transformation records (ORIGXn+SCALEn+MTRIXn)
# 51 - 55 Number of atomic coordinate records (ATOM+HETATM)
# 56 - 60 Number of TER records
# 61 - 65 Number of CONECT records
# 66 - 70 Number of SEQRES records
#
sub ParseMasterRecordLine {
  my($Line) = @_;

  if ($Line !~ /^MASTER/i) {
    return ((undef) x 11);
  }
  my($NumOfRemarkRecords, $NumOfHetRecords, $NumOfHelixRecords, $NumOfSheetRecords, $NumOfTurnRecords, $NumOfSiteRecords, $NumOfTransformationsRecords, $NumOfAtomAndHetatmRecords, $NumOfTerRecords, $NumOfConectRecords, $NumOfSeqresRecords) = map {s/ //g; $_} unpack("x6x4A5x5A5A5A5A5A5A5A5A5A5A5", $Line);

  return ($NumOfRemarkRecords, $NumOfHetRecords, $NumOfHelixRecords, $NumOfSheetRecords, $NumOfTurnRecords, $NumOfSiteRecords, $NumOfTransformationsRecords, $NumOfAtomAndHetatmRecords, $NumOfTerRecords, $NumOfConectRecords, $NumOfSeqresRecords);
}

# End record...
sub GenerateEndRecordLine {
  my($Line);
  $Line = 'END   ';
  return $Line;
}

# ATOM/HETATM record format:
#
# 1 - 6 Record name
# 7 - 11  Atom serial number - right justified
# 13 - 16 Atom name
# 17 Alternate location indicator.
# 18 - 20 Residue name - right justified
# 22 Chain identifier.
# 23 - 26 Residue sequence number - right justified
# 27 Code for insertion of residues.
# 31 - 38 Real(8.3), Orthogonal coordinates for X in Angstroms.
# 39 - 46 Real(8.3), Orthogonal coordinates for Y in Angstroms.
# 47 - 54 Real(8.3), Orthogonal coordinates for Z in Angstroms.
# 55 - 60 Real(6.2), Occupancy
# 61 - 66 Real(6.2), Temperature factor
# 73 - 76 LString(4), Segment identifier, left-justified.
# 77 - 78 LString(2), Element symbol, right-justified.
#79 - 80 LString(2), Charge on the atom.
#
# Notes:
#  . Atom names starting with C, N, O and S are left justified starting with column 14
#    and others are left justified starting with column 13.
#
#  . Six characters (columns) are reserved for atom names, assigned as follows:
#
#   13 - 14 Chemical symbol - right justified, except for hydrogen atoms
#
#   And for amino acids:
#
#   15 Remoteness indicator (alphabetic) (A, B, G, D, E, Z and so on)
#   16 Branch designator (numeric)
#
sub _ParseAtomOrHetatmRecordLine {
  my($Line) = @_;

  if ($Line !~ /^(ATOM|HETATM)/i) {
    return ((undef) x 15);
  }
  my($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge, $Length);

  ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = ('') x 15;

  $Length = length $Line;

  if ($Length <= 72) {
    ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor) = map {s/ //g; $_} unpack("x6A5xA4A1A3xA1A4A1x3A8A8A8A6A6", $Line);
  }
  else {
    ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = map {s/ //g; $_} unpack("x6A5xA4A1A3xA1A4A1x3A8A8A8A6A6x6A4A2A2", $Line);
  }
  return($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge);
}

# Generate ATOM/HETATM record line...
sub _GenerateAtomOrHetatmRecordLine {
  my($RecordType, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = @_;
  my($Line, $AtomNameFormat);

  if (length($AtomName) >= 4) {
    # Left justified starting at column 13 for all atom names of length 4...
    $AtomNameFormat = "%-4.4s";
  }
  elsif (IsEmpty($ElementSymbol)) {
    # No element symbol specified; just guess from atom name to cover most likely cases...
    $AtomNameFormat = ($AtomName =~ /^(C|N|O|S)/i) ? " %-3.3s" : "%-4.4s";
  }
  else {
    # Element symbol specified...
    if ($ElementSymbol =~ /^H$/i) {
      # Hydrogen atom name with <=3 characters is left justified starting at column 14;
      # Otherwise, left justified starting at column 13.
      $AtomNameFormat = (length($AtomName) <= 3) ? " %-3.3s" : "%-4.4s";
    }
    else {
      # Non-hydrogen atom name...
      $AtomNameFormat = (length($ElementSymbol) == 1) ? " %-3.3s" : "%-4.4s";
    }
  }

  $Line = sprintf "%-6.6s%5.5s ${AtomNameFormat}%1.1s%3.3s %1.1s%4.4s%1.1s   %8.8s%8.8s%8.8s%6.6s%6.6s      %-4.4s%2.2s%2.2s", $RecordType, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge;

  return $Line;
}

# Check record type...
sub _IsRecordType {
  my($Line, $SpecifiedType) = @_;
  my($Type, $Status);

  ($Type) = map {s/ //g; $_} unpack("A6", $Line);

  $Status = ($SpecifiedType eq $Type) ? 1 : 0;

  return $Status;
}

# Get record type...
sub _GetRecordType {
  my($Line) = @_;
  my($Type);

  ($Type) = map {s/ //g; $_} unpack("A6", $Line);

  return $Type;
}

# Get chains and residues data using ATOM/HETATM records...
#
sub _GetChainsAndResiduesFromAtomHetatmRecords {
  my($PDBRecordLinesRef, $GetChainResiduesBeyondTERFlag, $GetRecordLinesFlag) = @_;

  my($LineCount, $TotalChainCount, $PreviousResidueNumber, $ChainCount, $DefaultChainID, $DefaultChainLabel, $RecordLine, $AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge, %ChainsDataMap);

  # Do a quick chain count using TER record...
  $TotalChainCount = 0;
  LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
    if (IsEndmdlRecordType($RecordLine)) {
      last LINE;
    }
    if (IsTerRecordType($RecordLine)) {
      $TotalChainCount++;
    }
  }

  %ChainsDataMap = ();
  @{$ChainsDataMap{ChainIDs}} = ();
  %{$ChainsDataMap{Residues}} = ();
  %{$ChainsDataMap{ResidueNumbers}} = ();
  %{$ChainsDataMap{Lines}} = ();
  %{$ChainsDataMap{ResidueCount}} = ();

  $LineCount = 0;
  $ChainCount = 0;
  $DefaultChainLabel = 'None';
  $DefaultChainID = $DefaultChainLabel . ($ChainCount + 1);
  $PreviousResidueNumber = 0;

  LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
      $LineCount++;
      if (IsTerRecordType($RecordLine)) {
	$DefaultChainID = $DefaultChainLabel . ($ChainCount + 1);
	$ChainCount++;
	if ($ChainCount == $TotalChainCount) {
	  last LINE;
	}
	else {
	  next LINE;
	}
      }
      elsif (!(IsAtomRecordType($RecordLine) || IsHetatmRecordType($RecordLine))) {
	next LINE;
      }
      ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = ParseAtomRecordLine($RecordLine);

      if (IsEmpty($ChainID)) {
	$ChainID = $DefaultChainID;
      }
      if (exists $ChainsDataMap{Residues}{$ChainID}) {
	# Data for existing chain...
	if ($GetRecordLinesFlag) {
	  push @{$ChainsDataMap{Lines}{$ChainID}}, $RecordLine;
	}

	if ($ResidueNumber != $PreviousResidueNumber) {
	  # Next residue with in the chain...
	  push @{$ChainsDataMap{Residues}{$ChainID}}, $ResidueName;
	  push @{$ChainsDataMap{ResidueNumbers}{$ChainID}}, $ResidueNumber;

	  if (exists $ChainsDataMap{ResidueCount}{$ChainID}{$ResidueName}) {
	    $ChainsDataMap{ResidueCount}{$ChainID}{$ResidueName} += 1;
	  }
	  else {
	    $ChainsDataMap{ResidueCount}{$ChainID}{$ResidueName} = 1;
	  }
	  $PreviousResidueNumber = $ResidueNumber;
	}
      }
      else {
	# Data for new chain...
	push @{$ChainsDataMap{ChainIDs}}, $ChainID;

	@{$ChainsDataMap{Residues}{$ChainID}} = ();
	push @{$ChainsDataMap{Residues}{$ChainID}}, $ResidueName;

	@{$ChainsDataMap{ResidueNumbers}{$ChainID}} = ();
	push @{$ChainsDataMap{ResidueNumbers}{$ChainID}}, $ResidueNumber;

	@{$ChainsDataMap{Lines}{$ChainID}} = ();
	if ($GetRecordLinesFlag) {
	  push @{$ChainsDataMap{Lines}{$ChainID}}, $RecordLine;
	}

	%{$ChainsDataMap{ResidueCount}{$ChainID}} = ();
	$ChainsDataMap{ResidueCount}{$ChainID}{$ResidueName} = 1;
	$PreviousResidueNumber = $ResidueNumber;
      }
  }
  if (!$GetChainResiduesBeyondTERFlag) {
    return \%ChainsDataMap;
  }
  # Look for any HETATM residues specified outside TER records which could belong to an existing chain...
  my($LineIndex, $PreviousChainID);
  $PreviousChainID = '';
  $PreviousResidueNumber = 0;
  LINE: for $LineIndex (($LineCount - 1) .. $#{$PDBRecordLinesRef}) {
    $RecordLine = $PDBRecordLinesRef->[$LineIndex];
    if (IsEndmdlRecordType($RecordLine)) {
      last LINE;
    }
    if (!(IsAtomRecordType($RecordLine) || IsHetatmRecordType($RecordLine))) {
      next LINE;
    }
    ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) = ParseAtomRecordLine($RecordLine);
    if (IsEmpty($ChainID)) {
      # Ignore the chains with no ids...
      next LINE;
    }
    if (! exists($ChainsDataMap{Residues}{$ChainID})) {
      # Don't collect any new chains after TER record...
      next LINE;
    }
    if ($GetRecordLinesFlag) {
      push @{$ChainsDataMap{Lines}{$ChainID}}, $RecordLine;
    }
    if ($ResidueNumber != $PreviousResidueNumber || $ChainID ne $PreviousChainID) {

      push @{$ChainsDataMap{Residues}{$ChainID}}, $ResidueName;
      push @{$ChainsDataMap{ResidueNumbers}{$ChainID}}, $ResidueNumber;

      if (exists $ChainsDataMap{ResidueCount}{$ChainID}{$ResidueName}) {
	$ChainsDataMap{ResidueCount}{$ChainID}{$ResidueName} += 1;
      }
      else {
	$ChainsDataMap{ResidueCount}{$ChainID}{$ResidueName} = 1;
      }
      $PreviousChainID = $ChainID;
      $PreviousResidueNumber = $ResidueNumber;
    }
  }
  return \%ChainsDataMap;
}

# Get chains and residues data using SEQRES records...
#
sub _GetChainsAndResiduesFromSeqresRecords {
  my($PDBRecordLinesRef) = @_;

  my($ChainCount, $DefaultChainLabel, $DefaultChainID, $RecordLine, $RecordSerialNumber, $ChainID, $NumOfResidues, $ResidueName, $ResidueNamesString, @ResidueNamesList, %ChainsDataMap);

  %ChainsDataMap = ();
  @{$ChainsDataMap{ChainIDs}} = ();
  %{$ChainsDataMap{Residues}} = ();
  %{$ChainsDataMap{ResidueNumbers}} = ();
  %{$ChainsDataMap{ResidueCount}} = ();

  $ChainCount = 0;
  $DefaultChainLabel = 'None';
  $DefaultChainID = $DefaultChainLabel . ($ChainCount + 1);

  LINE: for $RecordLine (@{$PDBRecordLinesRef}) {
      if (!IsSeqresRecordType($RecordLine)) {
	next LINE;
      }
      ($RecordSerialNumber, $ChainID, $NumOfResidues, $ResidueNamesString) = ParseSeqresRecordLine($RecordLine);
      if ($RecordSerialNumber == 1) {
	# Indicates start of a new chain...
	$DefaultChainID = $DefaultChainLabel . ($ChainCount + 1);
	$ChainCount++;
      }
      if (IsEmpty($ChainID)) {
	$ChainID = $DefaultChainID;
      }
      # Process the residues...
      @ResidueNamesList = split /[ ]+/, $ResidueNamesString;

      if (exists $ChainsDataMap{Residues}{$ChainID}) {
	# Data for existing chain...
	push @{$ChainsDataMap{Residues}{$ChainID}}, @ResidueNamesList;
      }
      else {
	# Data for new chain...
	push @{$ChainsDataMap{ChainIDs}}, $ChainID;
	@{$ChainsDataMap{Residues}{$ChainID}} = ();
	push @{$ChainsDataMap{Residues}{$ChainID}}, @ResidueNamesList;
      }

      # Setup residue count...
      for $ResidueName (@ResidueNamesList) {
	if (exists $ChainsDataMap{ResidueCount}{$ChainID}{$ResidueName}) {
	  $ChainsDataMap{ResidueCount}{$ChainID}{$ResidueName} += 1;
	}
	else {
	  $ChainsDataMap{ResidueCount}{$ChainID}{$ResidueName} = 1;
	}
      }
  }
  return \%ChainsDataMap;
}

1;

__END__

=head1 NAME

PDBFileUtil

=head1 SYNOPSIS

use PDBFileUtil ;

use PDBFileUtil qw(:all);

=head1 DESCRIPTION

B<PDBFileUtil> module provides the following functions:

GenerateAtomOrHetatmRecordLine, GenerateAtomRecordLine, GenerateConectRecordLine,
GenerateEndRecordLine, GenerateHeaderRecordLine, GenerateHeaderRecordTimeStamp,
GenerateHetatmRecordLine, GenerateTerRecordLine, GetAllResidues,
GetChainsAndResidues, GetConectRecordLines, GetExperimentalTechnique,
GetExperimentalTechniqueResolution, GetMinMaxCoords, GetPDBRecordType,
GetRecordTypesCount, IsAtomRecordType, IsConectRecordType, IsEndmdlRecordType,
IsHeaderRecordType, IsHetatmRecordType, IsMasterRecordType, IsModelRecordType,
IsPDBFile, IsSeqresRecordType, IsTerRecordType, ParseAtomOrHetatmRecordLine,
ParseAtomRecordLine, ParseConectRecordLine, ParseExpdtaRecordLine,
ParseHeaderRecordLine, ParseHetatmRecordLine, ParseMasterRecordLine,
ParseRemark2ResolutionRecordLine, ParseSeqresRecordLine, ParseTerRecordLine,
ReadPDBFile

=head1 METHODS

=over 4

=item B<GenerateAtomOrHetatmRecordLine>

    $RecordLine = GenerateAtomOrHetatmRecordLine($RecordType,
      $AtomNumber, $AtomName, $AlternateLocation, $ResidueName,
      $ChainID, $ResidueNumber, $InsertionCode, $X, $Y, $Z,
      $Occupancy, $TemperatureFactor, $SegmentID,
      $ElementSymbol, $AtomCharge);

Returns ATOM or HETATM record line.

=item B<GenerateAtomRecordLine>

    $RecordLine = GenerateAtomRecordLine($AtomNumber,
      $AtomName, $AlternateLocation, $ResidueName, $ChainID,
      $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy,
      $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge);

Returns ATOM record line.

=item B<GenerateConectRecordLine>

    $RecordLine = GenerateConectRecordLine($AtomNum, $BondedAtomNum1,
      $BondedAtomNum2, $BondedAtomNum3, $BondedAtomNum4,
      $HBondedAtomNum1, $HBondedAtomNum2, $SaltBridgedAtomNum1,
      $HBondedAtomNum3, $HBondedAtomNum4, $SaltBridgedAtomNum2);

Returns CONECT record line.

=item B<GenerateHeaderRecordLine>

    $RecordLine = GenerateHeaderRecordLine($IDCode, [$Classification,
      $Date]);

Returns HEADER record line.

=item B<GenerateHeaderRecordTimeStamp>

    $Date = GenerateHeaderRecordTimeStamp();

Returns PDB header time stamp.

=item B<GenerateHetatmRecordLine>

    $RecordLine = GenerateHetatmRecordLine($AtomNumber, $AtomName,
    $AlternateLocation, $ResidueName, $ChainID, $ResidueNumber,
    $InsertionCode, $X, $Y, $Z, $Occupancy, $TemperatureFactor,
    $SegmentID, $ElementSymbol, $AtomCharge);

Returns HETATM record line.

=item B<GenerateEndRecordLine>

    $RecordLine = GenerateEndRecordLine();

Returns END record line.

=item B<GenerateTerRecordLine>

    $RecordLine = GenerateTerRecordLine($SerialNumber, [$ResidueName,
      $ChainID, $ResidueNumber, $InsertionCode]);

Returns TER record line.

=item B<GetAllResidues>

    $ResiduesDataRef = GetAllResidues($PDBRecordLinesRef);

Gets residue information using ATOM/HETATM records and returns a reference to a hash with
following key/value pairs:

    $ResiduesDataRef->{ResidueNames} - Array of all the residues
    $ResiduesDataRef->{ResidueCount}{$ResidueName} - Count of residues
    $ResiduesDataRef->{AtomResidueNames}} - Array of all ATOM residues
    $ResiduesDataRef->{AtomResidueCount}{$ResidueName} - Count of
       residues in ATOM records
    $ResiduesDataRef->{HetatomResidueNames} - List of all HETATM
       residues
    $ResiduesDataRef->{HetatmResidueCount}{$ResidueName} - Count of
      residues HETATM records

ATOM/HETATM records after the first ENDMDL records are simply ingnored.

=item B<GetChainsAndResidues>

    $ChainsDataRef = GetChainsAndResidues($PDBRecordLinesRef,
      [$RecordsSource, $GetChainResiduesBeyondTERFlag,
      $GetRecordLinesFlag]);

Gets chains and residue information using ATOM/HETATM or SEQRES records and returns a reference to a
hash with these keys:

    $ChainsDataRef->{ChainIDs} - List of chain IDs with 'None' for
      no IDs
    $ChainsDataRef->{Residues}{$ChainID} - List of residues in order
      of their appearance in a chain
    $ChainsDataRef->{ResidueCount}{$ChainID}{$ResidueName} - Count of
      residues in a chain

Chains and residue data can be extacted using either ATOM/HETATM records or SEQRES records.
ATOM/HETATM records after the first ENDMDL records are simply ingnored.

=item B<GetConectRecordLines>

    $ConectRecordLinesRef = GetConectRecordLines($PDBRecordLinesRef,
      $AtomNumbersMapRef);

Collects CONECT record lines for specific atom number, modified specified data to exclude any atom
number not present in the list of specified atom numbers and returns a reference to list of
CONECT record lines.

=item B<GetExperimentalTechnique>

    $ExperimentalTechnique = GetExperimentalTechnique($PDBRecordLinesRef);

Returns I<ExperimentalTechnique> value retrieved from EXPDATA record line.

=item B<GetExperimentalTechniqueResolution>

    ($Resolution, $ResolutionUnits) = GetExperimentalTechniqueResolution(
                                      $PDBRecordLinesRef);

Returns I<Resolution> and I<ResolutionUnits> values from REMARK 2 RESOLUTION
record line.

=item B<GetMinMaxCoords>

    ($XMin, $YMin, $ZMin, $XMax, $YMax, $ZMax) =
      GetMinMaxCoords($PDBRecordLinesRef);

Returns minimum and maximum XYZ coordinates for ATOM/HETATM records.

=item B<GetPDBRecordType>

    $RecordType = GetPDBRecordType($RecordLine);

Returns type of I<RecordLine>.

=item B<GetRecordTypesCount>

    $RecordTypeDataRef = GetRecordTypesCount($PDBRecordLinesRef,
      [$SpecifiedRecordType, $GetRecordLinesFlag]);

Counts the number of each record type or a $SpecifiedRecordType and returns a reference to data
type with following key/value pairs:

    $RecordTypeDataRef->{RecordTypes} - An array of unique record types
       in order of their presence in the file
    $RecordTypeDataRef->{Count}{$RecordType} - Count of each record type
    $RecordTypeDataRef->{Lines}{$RecordType} - Optional lines data for a
      specific record type.

=item B<IsAtomRecordType>

    $Status = IsAtomRecordType($RecordLine);

Returns 1 or 0 based on whether it's a ATOM record line.

=item B<IsConectRecordType>

    $Status = IsAtomConectType($RecordLine);

Returns 1 or 0 based on whether it's a CONECT record line.

=item B<IsEndmdlRecordType>

    $Status = IsEndmdlRecordType($RecordLine);

Returns 1 or 0 based on whether it's a ENDMDL a record line.

=item B<IsHeaderRecordType>

    $Status = IsHeaderRecordType($RecordLine);

Returns 1 or 0 based on whether it's a HEADER a record line.

=item B<IsHetatmRecordType>

    $Status = IsHetatmRecordType($RecordLine);

Returns 1 or 0 based on whether it's a HETATM a record line.

=item B<IsMasterRecordType>

    $Status = IsMasterRecordType($RecordLine);

Returns 1 or 0 based on whether it's a MASTER a record line.

=item B<IsModelRecordType>

    $Status = IsModelRecordType($RecordLine);

Returns 1 or 0 based on whether it's a MODEL record line.

=item B<IsPDBFile>

    $Status = IsPDBFile($PDBFile);

Returns 1 or 0 based on whether it's a PDB file.

=item B<IsSeqresRecordType>

    $Status = IsSeqresRecordType($RecordLine);

Returns 1 or 0 based on whether it's SEQRES a record line.

=item B<IsTerRecordType>

    $Status = IsTerRecordType($RecordLine);

Returns 1 or 0 based on whether it's a TER record line.

=item B<ParseAtomOrHetatmRecordLine>

    ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID,
      $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy,
      $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) =
      ParseAtomOrHetatmRecordLine($RecordLine);

Parses ATOM or HETATM record line.

=item B<ParseAtomRecordLine>

    ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID,
      $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy,
      $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) =
      ParseAtomRecordLine($RecordLine);

Parses ATOM record line.

=item B<ParseConectRecordLine>

    ($AtomNum, $BondedAtomNum1, $BondedAtomNum2, $BondedAtomNum3,
       $BondedAtomNum4, $HBondedAtomNum1, $HBondedAtomNum2,
       $SaltBridgedAtomNum1, $HBondedAtomNum3, $HBondedAtomNum4,
       $SaltBridgedAtomNum2) = ParseConectRecordLine($RecordLine);

Parses CONECT record line.

=item B<ParseExpdtaRecordLine>

    ($ContinuationNum, $ExperimentalTechnique) = ParseExpdtaRecordLine($Line);

Parses EXPDTA record line.

=item B<ParseHeaderRecordLine>

    ($Classification, $DepositionDate, $IDCode) = ParseHeaderRecordLine($RecordLine);

Parses HEADER record line

=item B<ParseHetatmRecordLine>

    ($AtomNumber, $AtomName, $AlternateLocation, $ResidueName, $ChainID,
      $ResidueNumber, $InsertionCode, $X, $Y, $Z, $Occupancy,
      $TemperatureFactor, $SegmentID, $ElementSymbol, $AtomCharge) =
      ParseHetatmRecordLine($RecordLine);

Parses HETATM record line.

=item B<ParseMasterRecordLine>

    ($NumOfRemarkRecords, $NumOfHetRecords, $NumOfHelixRecords,
      $NumOfSheetRecords, $NumOfTurnRecords, $NumOfSiteRecords,
      $NumOfTransformationsRecords, $NumOfAtomAndHetatmRecords,
      $NumOfTerRecords, $NumOfConectRecords, $NumOfSeqresRecords) =
      ParseMasterRecordLine($RecordLine);

Parses MASTER ecord line.

=item B<ParseRemark2ResolutionRecordLine>

    ($Resolution, $ResolutionUnits) = ParseRemark2ResolutionRecordLine(
                                      $RecordLine);

Parses REMARK 2 RESOLUTION record line.

=item B<ParseSeqresRecordLine>

    ($RecordSerialNumber, $ChainID, $NumOfResidues, $ResidueNames) =
      ParseSeqresRecordLine($RecordLine);

Parses SEQRES record line.

=item B<ParseTerRecordLine>

    ($SerialNumber, $ResidueName, $ChainID, $ResidueNumber, $InsertionCode) =
      ParseTerRecordLine($RecordLine);

Parses TER record line.

=item B<ReadPDBFile>

    $PDBRecordLinesRef = ReadPDBFile($PDBFile);

Reads PDB file and returns reference to record lines.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

FileUtil.pm, SequenceFileUtil.pm, TextUtil.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
