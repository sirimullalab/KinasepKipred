package SDFileUtil;
#
# File: SDFileUtil.pm
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
use Carp;
use PeriodicTable qw(IsElement);
use TimeUtil ();

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw(GenerateCmpdAtomLine GenerateCmpdBondLine GenerateCmpdChargePropertyLines GenerateCmpdCommentsLine GenerateCmpdCountsLine GenerateCmpdAtomAliasPropertyLines GenerateCmpdIsotopePropertyLines GenerateCmpdDataHeaderLabelsAndValuesLines GenerateCmpdMiscInfoLine GenerateCmpdRadicalPropertyLines GenerateCmpdMolNameLine GenerateEmptyCtabBlockLines GenerateMiscLineDateStamp GetAllAndCommonCmpdDataHeaderLabels GetCmpdDataHeaderLabels GetCmpdDataHeaderLabelsAndValues GetCmpdFragments GetCtabLinesCount GetUnknownAtoms GetInvalidAtomNumbers MDLChargeToInternalCharge InternalChargeToMDLCharge MDLBondTypeToInternalBondOrder InternalBondOrderToMDLBondType MDLBondStereoToInternalBondStereochemistry InternalBondStereochemistryToMDLBondStereo InternalSpinMultiplicityToMDLRadical MDLRadicalToInternalSpinMultiplicity IsCmpd3D IsCmpd2D ParseCmpdAtomLine ParseCmpdBondLine ParseCmpdCommentsLine ParseCmpdCountsLine ParseCmpdMiscInfoLine ParseCmpdMolNameLine ParseCmpdAtomAliasPropertyLine ParseCmpdChargePropertyLine ParseCmpdIsotopePropertyLine ParseCmpdRadicalPropertyLine ReadCmpdString RemoveCmpdDataHeaderLabelAndValue WashCmpd);
@EXPORT_OK = qw();
%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Format data for compounds count line...
sub GenerateCmpdCountsLine {
  my($AtomCount, $BondCount, $ChiralFlag, $PropertyCount, $Version, $Line);

  if (@_ == 5) {
    ($AtomCount, $BondCount, $ChiralFlag, $PropertyCount, $Version) = @_;
  }
  elsif (@_ == 3) {
    ($AtomCount, $BondCount, $ChiralFlag) = @_;
    $PropertyCount = 999;
    $Version = "V2000";
  }
  else {
    ($AtomCount, $BondCount) = @_;
    $ChiralFlag = 0;
    $PropertyCount = 999;
    $Version = "V2000";
  }
  if ($AtomCount > 999) {
    croak "Error: SDFileUtil::GenerateCmpdCountsLine: The atom count, $AtomCount, exceeds maximum of 999 allowed for CTAB version 2000. The Extended Connection Table (V3000) format in MDL MOL and SD files is not supported by the current release of MayaChemTools...";
  }
  $Line = sprintf "%3i%3i%3i%3i%3i%3i%3i%3i%3i%3i%3i%6s", $AtomCount, $BondCount, 0, 0, $ChiralFlag, 0, 0, 0, 0, 0, $PropertyCount, $Version;

  return ($Line);
}

# Generate comments line...
sub GenerateCmpdCommentsLine {
  my($Comments) = @_;
  my($Line);

  $Line = (length($Comments) > 80) ? substr($Comments, 0, 80) : $Comments;

  return $Line;
}

# Generate molname line...
sub GenerateCmpdMolNameLine {
  my($MolName) = @_;
  my($Line);

  $Line = (length($MolName) > 80) ? substr($MolName, 0, 80) : $MolName;

  return $Line;
}

# Generate data for compounds misc info line...
sub GenerateCmpdMiscInfoLine {
  my($ProgramName, $UserInitial, $Code) = @_;
  my($Date, $Line);

  if (!(defined($ProgramName) && $ProgramName)) {
    $ProgramName = "MayaChem";
  }
  if (!(defined($UserInitial) && $UserInitial)) {
    $UserInitial = "  ";
  }
  if (!(defined($Code) && $Code)) {
    $Code = "2D";
  }

  if (length($ProgramName) > 8) {
    $ProgramName = substr($ProgramName, 0, 8);
  }
  if (length($UserInitial) > 2) {
    $UserInitial = substr($UserInitial, 0, 2);
  }
  if (length($Code) > 2) {
    $Code = substr($Code, 0, 2);
  }
  $Date = GenerateMiscLineDateStamp();

  $Line = "${UserInitial}${ProgramName}${Date}${Code}";

  return $Line;
}

# Generate data for compounds misc info line...
sub GenerateEmptyCtabBlockLines {
  my($Date, $Lines);

  if (@_ == 1) {
    ($Date) = @_;
  }
  else {
    $Date = GenerateMiscLineDateStamp();
  }
  # First line: Blank molname line...
  # Second line: Misc info...
  # Third line: Blank comments line...
  # Fourth line: Counts line reflecting empty structure data block...
  $Lines = "\n";
  $Lines .= "  MayaChem${Date}2D\n";
  $Lines .= "\n";
  $Lines .= GenerateCmpdCountsLine(0, 0, 0) . "\n";
  $Lines .= "M  END";

  return $Lines;
}

# Generate SD file data stamp...
sub GenerateMiscLineDateStamp {
  return TimeUtil::SDFileTimeStamp();
}

# Generate data for compound atom line...
#
sub GenerateCmpdAtomLine {
  my($AtomSymbol, $AtomX, $AtomY, $AtomZ, $MassDifference, $Charge, $StereoParity) = @_;
  my($Line);

  if (!defined $MassDifference) {
    $MassDifference = 0;
  }
  if (!defined $Charge) {
    $Charge = 0;
  }
  if (!defined $StereoParity) {
    $StereoParity = 0;
  }
  $Line = sprintf "%10.4f%10.4f%10.4f %-3s%2i%3i%3i  0  0  0  0  0  0  0  0  0", $AtomX, $AtomY, $AtomZ, $AtomSymbol, $MassDifference, $Charge, $StereoParity;

  return $Line
}

# Generate data for compound bond line...
#
sub GenerateCmpdBondLine {
  my($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo) = @_;
  my($Line);

  if (!defined $BondStereo) {
    $BondStereo = 0;
  }
  $Line = sprintf "%3i%3i%3i%3i  0  0  0", $FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo;

  return $Line
}

# Generate charge property lines for CTAB block...
#
sub GenerateCmpdChargePropertyLines {
  my($ChargeValuePairsRef) = @_;

  return _GenerateCmpdGenericPropertyLines('Charge', $ChargeValuePairsRef);
}

# Generate isotope property lines for CTAB block...
#
sub GenerateCmpdIsotopePropertyLines {
  my($IsotopeValuePairsRef) = @_;

  return _GenerateCmpdGenericPropertyLines('Isotope', $IsotopeValuePairsRef);
}

# Generate radical property line property lines for CTAB block...
#
sub GenerateCmpdRadicalPropertyLines {
  my($RadicalValuePairsRef) = @_;

  return _GenerateCmpdGenericPropertyLines('Radical', $RadicalValuePairsRef);
}

# Generate atom alias property line property lines for CTAB block...
#
# Atom alias property line format:
#
# A  aaa
# x...
#
#    aaa: Atom number
#    x: Atom alias in next line
#
sub GenerateCmpdAtomAliasPropertyLines {
  my($PropertyValuePairsRef) = @_;
  my($Index, $AtomNum, $AtomAlias, $Line, @PropertyLines);

  @PropertyLines = ();

  for ($Index = 0; $Index < $#{$PropertyValuePairsRef}; $Index += 2) {
    $AtomNum = $PropertyValuePairsRef->[$Index];
    $AtomAlias = $PropertyValuePairsRef->[$Index + 1];

    $Line = "A  " . sprintf "%3i", $AtomNum;

    push @PropertyLines, $Line;
    push @PropertyLines, $AtomAlias;
  }

  return @PropertyLines;
}

# Generate data header labels and values lines...
#
sub GenerateCmpdDataHeaderLabelsAndValuesLines {
  my($DataHeaderLabelsRef, $DataHeaderLabelsAndValuesRef, $SortDataLabels) = @_;
  my($DataLabel, $DataValue, @DataLabels, @DataLines);

  if (!defined $SortDataLabels) {
    $SortDataLabels = 0;
  }

  @DataLines = ();
  @DataLabels = ();
  if ($SortDataLabels) {
    push @DataLabels, sort @{$DataHeaderLabelsRef};
  }
  else {
    push @DataLabels,  @{$DataHeaderLabelsRef};
  }
  for $DataLabel (@DataLabels) {
    $DataValue = '';
    if (exists $DataHeaderLabelsAndValuesRef->{$DataLabel}) {
      $DataValue = $DataHeaderLabelsAndValuesRef->{$DataLabel};
    }
    push @DataLines, (">  <${DataLabel}>", "$DataValue", "");
  }
  return @DataLines;
}

# Parse data field header in SD file and return lists of all and common data field
# labels.
sub GetAllAndCommonCmpdDataHeaderLabels {
  my($SDFileRef) = @_;
  my($CmpdCount, $CmpdString, $Label, @CmpdLines, @DataFieldLabels, @CommonDataFieldLabels, %DataFieldLabelsMap);

  $CmpdCount = 0;
  @DataFieldLabels = ();
  @CommonDataFieldLabels = ();
  %DataFieldLabelsMap = ();

  while ($CmpdString = ReadCmpdString($SDFileRef)) {
    $CmpdCount++;
    @CmpdLines = split "\n", $CmpdString;
    # Process compound data header labels and figure out which ones are present for
    # all the compounds...
    if (@DataFieldLabels) {
      my (@CmpdDataFieldLabels) = GetCmpdDataHeaderLabels(\@CmpdLines);
      my(%CmpdDataFieldLabelsMap) = ();
      # Setup a map for the current labels...
      for $Label (@CmpdDataFieldLabels) {
	$CmpdDataFieldLabelsMap{$Label} = "PresentInSome";
      }
      # Check the presence old labels for this compound; otherwise, mark 'em new...
      for $Label (@DataFieldLabels) {
	if (!$CmpdDataFieldLabelsMap{$Label}) {
	  $DataFieldLabelsMap{$Label} = "PresentInSome";
	}
      }
      # Check the presence this compound in the old labels; otherwise, add 'em...
      for $Label (@CmpdDataFieldLabels ) {
	if (!$DataFieldLabelsMap{$Label}) {
	  # It's a new label...
	  push @DataFieldLabels, $Label;
	  $DataFieldLabelsMap{$Label} = "PresentInSome";
	}
      }
    }
    else {
      # Get the initial label set and set up a map...
      @DataFieldLabels = GetCmpdDataHeaderLabels(\@CmpdLines);
      for $Label (@DataFieldLabels) {
	$DataFieldLabelsMap{$Label} = "PresentInAll";
      }
    }
  }
  # Identify the common data field labels...
  @CommonDataFieldLabels = ();
  for $Label (@DataFieldLabels) {
    if ($DataFieldLabelsMap{$Label} eq "PresentInAll") {
      push @CommonDataFieldLabels, $Label;
    }
  }
  return ($CmpdCount, \@DataFieldLabels, \@CommonDataFieldLabels);
}

# Parse all the data header labels and return 'em as an list...
#
# Format:
#
#> Data header line
#Data line(s)
#Blank line
#
# [Data Header] (one line) precedes each item of data, starts with a greater than (>) sign, and
# contains at least one of the following:
#  The field name enclosed in angle brackets. For example: <melting.point>
#  The field number, DTn , where n represents the number assigned to the field in a MACCS-II database
#
#Optional information for the data header includes:
#  The compound’s external and internal registry numbers. External registry numbers must be enclosed in parentheses.
#  Any combination of information
#
#The following are examples of valid data headers:
#> <MELTING.POINT>
#> 55 (MD-08974) <BOILING.POINT> DT12
#> DT12 55
#> (MD-0894) <BOILING.POINT> FROM ARCHIVES
#
#Notes: Sometimes last blank line is missing and can be just followed by $$$$
#
sub GetCmpdDataHeaderLabels {
  my($CmpdLines) = @_;
  my($CmpdLine, $Label, @Labels);

  @Labels = ();
  CMPDLINE: for $CmpdLine (@$CmpdLines) {
    if ($CmpdLine !~ /^>/) {
      next CMPDLINE;
    }
    # Does the line contains field name enclosed in angular brackets?
    ($Label) = $CmpdLine =~ /<.*?>/g;
    if (!defined($Label)) {
      next CMPDLINE;
    }
    $Label =~ s/(<|>)//g;
    push @Labels, $Label;
  }
  return (@Labels);
}

# Parse all the data header labels and values
sub GetCmpdDataHeaderLabelsAndValues {
  my($CmpdLines) = @_;
  my($CmpdLine, $CurrentLabel, $Label, $Value, $ValueCount, $ProcessingLabelData, @Values, %DataFields);

  %DataFields = ();
  $ProcessingLabelData = 0;
  $ValueCount = 0;
  CMPDLINE: for $CmpdLine (@$CmpdLines) {
    if ($CmpdLine =~ /^\$\$\$\$/) {
      last CMPDLINE;
    }
    if ($CmpdLine =~ /^>/) {
      # Does the line contains field name enclosed in angular brackets?
      ($Label) = $CmpdLine =~ /<.*?>/g;
      if (defined $Label) {
	$CurrentLabel = $Label;
	$CurrentLabel =~ s/(<|>)//g;
	$ProcessingLabelData = 0;
	$ValueCount = 0;

	if ($CurrentLabel) {
	  $ProcessingLabelData = 1;
	  $DataFields{$CurrentLabel} = '';
	  next CMPDLINE;
	}
      }
      else {
	if (!$ProcessingLabelData) {
	  # Data line containing no <label> as allowed by SDF format. Just ignore it...
	  next CMPDLINE;
	}
      }
    }
    if (!$ProcessingLabelData) {
      next CMPDLINE;
    }
    if (!(defined($CmpdLine) && length($CmpdLine))) {
      # Blank line terminates value for a label...
      $CurrentLabel = '';
      $ValueCount = 0;
      $ProcessingLabelData = 0;
      next CMPDLINE;
    }
    $ValueCount++;
    $Value = $CmpdLine;

    if ($ValueCount > 1) {
      $DataFields{$CurrentLabel} .= "\n" . $Value;
    }
    else {
      $DataFields{$CurrentLabel} = $Value;
    }
  }
  return (%DataFields);
}

# Return an updated compoud string after removing  data header label along with its
# value from the specified compound string...
#
sub RemoveCmpdDataHeaderLabelAndValue {
  my($CmpdString, $DataHeaderLabel) = @_;
  my($Line, $PorcessingDataHeaderLabel, @CmpdLines);

  @CmpdLines = ();
  $PorcessingDataHeaderLabel = 0;

  CMPDLINE: for $Line (split "\n", $CmpdString) {
    if ($Line =~ /^>/ && $Line =~ /<$DataHeaderLabel>/i) {
      $PorcessingDataHeaderLabel = 1;
      next CMPDLINE;
    }

    if ($PorcessingDataHeaderLabel) {
      # Blank line indicates end of fingerprints data value...
      if ($Line =~ /^\$\$\$\$/) {
	push @CmpdLines, $Line;
	$PorcessingDataHeaderLabel = 0;
      }
      elsif (!length($Line)) {
	$PorcessingDataHeaderLabel = 0;
      }
      next CMPDLINE;
    }

    # Track compound lines without fingerprints data...
    push @CmpdLines, $Line;
  }

  return join "\n", @CmpdLines;
}

#
# Using bond blocks, figure out the number of disconnected fragments  and
# return their values along with the atom numbers in a string delimited by new
# line character.
#
sub GetCmpdFragments {
  my($CmpdLines) = @_;
  my($AtomCount, $BondCount, $FirstAtomNum, $SecondAtomNum, @AtomConnections, $BondType, $FragmentString, $FragmentCount, $LineIndex, $Index, $AtomNum, $NbrAtomNum, @ProcessedAtoms, $ProcessedAtomCount, $ProcessAtomNum, @ProcessingAtoms, @ConnectedAtoms, %Fragments, $FragmentNum, $AFragmentString);

  # Setup the connection table for each atom...
  @AtomConnections = ();
  ($AtomCount, $BondCount) = ParseCmpdCountsLine(@$CmpdLines[3]);
  for $AtomNum (1 .. $AtomCount) {
    %{$AtomConnections[$AtomNum]} = ();
  }
  for ($LineIndex = 4 + $AtomCount; $LineIndex < (4 + $AtomCount + $BondCount); $LineIndex++) {
    ($FirstAtomNum, $SecondAtomNum, $BondType) = ParseCmpdBondLine(@$CmpdLines[$LineIndex]);
    if (!$AtomConnections[$FirstAtomNum]{$SecondAtomNum}) {
      $AtomConnections[$FirstAtomNum]{$SecondAtomNum} = $BondType;
    }
    if (!$AtomConnections[$SecondAtomNum]{$FirstAtomNum}) {
      $AtomConnections[$SecondAtomNum]{$FirstAtomNum} = $BondType;
    }
  }

  #Get set to count fragments...
  $ProcessedAtomCount = 0;
  $FragmentNum = 0;
  %Fragments = ();
  @ProcessedAtoms = ();
  for $AtomNum (1 .. $AtomCount) {
    $ProcessedAtoms[$AtomNum] = 0;
  }
  while ($ProcessedAtomCount < $AtomCount) {
    @ProcessingAtoms = ();
    @ConnectedAtoms = ();
    ATOMNUM: for $AtomNum (1 .. $AtomCount) {
      if (!$ProcessedAtoms[$AtomNum]) {
	$ProcessedAtomCount++;
	$ProcessedAtoms[$AtomNum] = 1;
	push @ProcessingAtoms, $AtomNum;
	$FragmentNum++;
	@{$Fragments{$FragmentNum} } = ();
	push @{$Fragments{$FragmentNum} }, $AtomNum;
	last ATOMNUM;
      }
    }

    # Go over the neighbors and follow the connection trail while collecting the
    # atoms numbers present in the connected fragment...
    while (@ProcessingAtoms) {
      for ($Index = 0; $Index < @ProcessingAtoms; $Index++) {
	$ProcessAtomNum = $ProcessingAtoms[$Index];
	for $NbrAtomNum (keys %{$AtomConnections[$ProcessAtomNum]})  {
	  if (!$ProcessedAtoms[$NbrAtomNum]) {
	    $ProcessedAtomCount++;
	    $ProcessedAtoms[$NbrAtomNum] = 1;
	    push @ConnectedAtoms, $NbrAtomNum;
	    push @{ $Fragments{$FragmentNum} }, $NbrAtomNum;
	  }
	}
      }
      @ProcessingAtoms = ();
      @ProcessingAtoms = @ConnectedAtoms;
      @ConnectedAtoms = ();
    }
  }
  $FragmentCount = $FragmentNum;
  $FragmentString = "";

  # Sort out the fragments by size...
  for $FragmentNum (sort { @{$Fragments{$b}} <=> @{$Fragments{$a}}  } keys %Fragments ) {
    # Sort the atoms in a fragment by their numbers...
    $AFragmentString = join " ", sort { $a <=> $b } @{ $Fragments{$FragmentNum} };
    if ($FragmentString) {
      $FragmentString .=  "\n" . $AFragmentString;
    }
    else {
      $FragmentString = $AFragmentString;
    }
  }
  return ($FragmentCount, $FragmentString);
}

# Count number of lines present in between 4th and line containg "M END"
sub GetCtabLinesCount {
  my($CmpdLines) = @_;
  my($LineIndex, $CtabLinesCount);

  $CtabLinesCount = 0;
 LINE: for ($LineIndex = 4; $LineIndex < @$CmpdLines; $LineIndex++) {
    #
    # Any line after atom and bond data starting with anything other than space or
    # a digit indicates end of Ctab atom/bond data block...
    #
    if (@$CmpdLines[$LineIndex] !~ /^[0-9 ]/) {
      $CtabLinesCount = $LineIndex - 4;
      last LINE;
    }
  }
  return $CtabLinesCount;
}

# Using atom blocks, count the number of atoms which contain special element
# symbols not present in the periodic table.
sub GetUnknownAtoms {
  my($CmpdLines) = @_;
  my($UnknownAtomCount, $UnknownAtoms, $UnknownAtomLines, $LineIndex, $AtomCount, $AtomSymbol);

  $UnknownAtomCount = 0;
  $UnknownAtoms = "";
  $UnknownAtomLines = "";
  ($AtomCount) = ParseCmpdCountsLine(@$CmpdLines[3]);
  for ($LineIndex = 4; $LineIndex < (4 + $AtomCount); $LineIndex++) {
    ($AtomSymbol) = ParseCmpdAtomLine(@$CmpdLines[$LineIndex]);
    if (!IsElement($AtomSymbol)) {
      $UnknownAtomCount++;
      $UnknownAtoms .= " $AtomSymbol";
      if ($UnknownAtomLines) {
	$UnknownAtomLines .= "\n" . @$CmpdLines[$LineIndex];
      }
      else {
	$UnknownAtomLines = @$CmpdLines[$LineIndex];
      }
    }
  }
  return ($UnknownAtomCount, $UnknownAtoms, $UnknownAtomLines);
}

# Check z coordinates of all atoms to see whether any of them is non-zero
# which makes the compound geometry three dimensional...
#
sub IsCmpd3D {
  my($CmpdLines) = @_;
  my($LineIndex, $AtomCount, $AtomSymbol, $AtomX, $AtomY, $AtomZ);

  ($AtomCount) = ParseCmpdCountsLine(@$CmpdLines[3]);
  for ($LineIndex = 4; $LineIndex < (4 + $AtomCount); $LineIndex++) {
    ($AtomSymbol, $AtomX, $AtomY, $AtomZ) = ParseCmpdAtomLine(@$CmpdLines[$LineIndex]);
    if ($AtomZ != 0) {
      return 1;
    }
  }
  return 0;
}

# Check whether it's a 2D compound...
#
sub IsCmpd2D {
  my($CmpdLines) = @_;

  return IsCmpd3D($CmpdLines) ? 0 : 1;
}

# Using bond blocks, count the number of bond lines which contain atom numbers
# greater than atom count specified in compound count line...
#
sub GetInvalidAtomNumbers {
  my($CmpdLines) = @_;
  my($LineIndex, $AtomCount, $BondCount, $FirstAtomNum, $SecondAtomNum, $InvalidAtomNumbersCount, $InvalidAtomNumbers, $InvalidAtomNumberLines, $Line, $InvalidAtomPropertyLine, $ValuePairIndex, $AtomNum, $Value, @ValuePairs);

  ($AtomCount, $BondCount) = ParseCmpdCountsLine(@$CmpdLines[3]);

  $InvalidAtomNumbersCount = 0;
  $InvalidAtomNumbers = "";
  $InvalidAtomNumberLines = "";

  # Go over bond block lines...
  LINE: for ($LineIndex = 4 + $AtomCount; $LineIndex < (4 + $AtomCount + $BondCount); $LineIndex++) {
    ($FirstAtomNum, $SecondAtomNum) = ParseCmpdBondLine(@$CmpdLines[$LineIndex]);
    if ($FirstAtomNum <= $AtomCount && $SecondAtomNum <= $AtomCount) {
      next LINE;
    }
    if ($FirstAtomNum > $AtomCount) {
      $InvalidAtomNumbersCount++;
      $InvalidAtomNumbers .= " $FirstAtomNum";
    }
    if ($SecondAtomNum > $AtomCount) {
      $InvalidAtomNumbersCount++;
      $InvalidAtomNumbers .= " $SecondAtomNum";
    }
    if ($InvalidAtomNumberLines) {
      $InvalidAtomNumberLines .= "\n" . @$CmpdLines[$LineIndex];
    }
    else {
      $InvalidAtomNumberLines = @$CmpdLines[$LineIndex];
    }
  }
  # Go over property lines before M  END...
  #
  LINE: for ($LineIndex = (4 + $AtomCount + $BondCount); $LineIndex < @$CmpdLines; $LineIndex++) {
    $Line = @$CmpdLines[$LineIndex];
    @ValuePairs = ();
    if ($Line =~ /^M  END/i) {
      last LINE;
    }
    @ValuePairs = ();
    if ($Line =~ /^M  CHG/i) {
      @ValuePairs = ParseCmpdChargePropertyLine($Line);
    }
    elsif ($Line =~ /^M  RAD/i) {
      @ValuePairs = ParseCmpdRadicalPropertyLine($Line);
    }
    elsif ($Line =~ /^M  ISO/i) {
      @ValuePairs = ParseCmpdIsotopePropertyLine($Line);
    }
    elsif ($Line =~ /^A  /i) {
      my($NextLine);
      $LineIndex++;
      $NextLine = @$CmpdLines[$LineIndex];
      @ValuePairs = ParseCmpdAtomAliasPropertyLine($Line, $NextLine);
    }
    else {
      next LINE;
    }

    $InvalidAtomPropertyLine = 0;
    for ($ValuePairIndex = 0; $ValuePairIndex < $#ValuePairs; $ValuePairIndex += 2) {
      $AtomNum = $ValuePairs[$ValuePairIndex]; $Value = $ValuePairs[$ValuePairIndex + 1];
      if ($AtomNum > $AtomCount) {
	$InvalidAtomPropertyLine = 1;
	$InvalidAtomNumbersCount++;
	$InvalidAtomNumbers .= " $AtomNum";
      }
    }
    if ($InvalidAtomPropertyLine) {
      if ($InvalidAtomNumberLines) {
	$InvalidAtomNumberLines .= "\n" . $Line;
      }
      else {
	$InvalidAtomNumberLines = $Line;
      }
    }
  }

  return ($InvalidAtomNumbersCount, $InvalidAtomNumbers, $InvalidAtomNumberLines);
}

# Ctab lines: Atom block
#
# Format: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
#         A10       A10       A10       xA3 A2A3 A3 A3 A3 A3 A3 A3 A3 A3 A3 A3
# x,y,z: Atom coordinates
# aaa: Atom symbol. Entry in periodic table or L for atom list, A, Q, * for unspecified
#      atom, and LP for lone pair, or R# for Rgroup label
# dd: Mass difference. -3, -2, -1, 0, 1, 2, 3, 4 (0 for value beyond these limits)
# ccc: Charge. 0 = uncharged or value other than these, 1 = +3, 2 = +2, 3 = +1,
#      4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
# sss: Atom stereo parity. 0 = not stereo, 1 = odd, 2 = even, 3 = either or unmarked stereo center
# hhh: Hydrogen count + 1. 1 = H0, 2 = H1, 3 = H2, 4 = H3, 5 = H4
# bbb: Stereo care box. 0 = ignore stereo configuration of this double bond atom, 1 = stereo
#      configuration of double bond atom must match
# vvv: Valence. 0 = no marking (default)(1 to 14) = (1 to 14) 15 = zero valence
# HHH: H0 designator. 0 = not specified, 1 = no H atoms allowed (redundant due to hhh)
# rrr: Not used
# iii: Not used
# mmm: Atom-atom mapping number. 1 - number of atoms
# nnn: Inversion/retention flag. 0 = property not applied, 1 = configuration is inverted,
#      2 = configuration is retained.
# eee: Exact change flag. 0 = property not applied, 1 = change on atom must be
#      exactly as shown
#
# Notes:
#  . StereoParity: 1 - ClockwiseStereo, 2 - AntiClockwiseStereo; 3 - Either; 0 - none. These
#    values determine chirailty around the chiral center; a non zero value indicates atom
#    has been marked as chiral center.
#
sub ParseCmpdAtomLine {
  my($Line) = @_;
  my ($LineIndex, $AtomX, $AtomY, $AtomZ, $AtomSymbol, $MassDifference, $Charge, $StereoParity);

  ($AtomX, $AtomY, $AtomZ, $AtomSymbol, $MassDifference, $Charge, $StereoParity) = ('') x 7;
  if (length($Line) > 31) {
    ($AtomX, $AtomY, $AtomZ, $AtomSymbol, $MassDifference, $Charge, $StereoParity) = unpack("A10A10A10xA3A2A3A3", $Line);
  }
  else {
    ($AtomX, $AtomY, $AtomZ, $AtomSymbol) = unpack("A10A10A10", $Line);
  }
  return ($AtomSymbol, $AtomX, $AtomY, $AtomZ, $MassDifference, $Charge, $StereoParity);
}

# Map MDL charge value used in SD and MOL files to internal charge used by MayaChemTools.
#
sub MDLChargeToInternalCharge {
  my($MDLCharge) = @_;
  my($InternalCharge);

  CHARGE: {
    if ($MDLCharge == 0) { $InternalCharge = 0; last CHARGE;}
    if ($MDLCharge == 1) { $InternalCharge = 3; last CHARGE;}
    if ($MDLCharge == 2) { $InternalCharge = 2; last CHARGE;}
    if ($MDLCharge == 3) { $InternalCharge = 1; last CHARGE;}
    if ($MDLCharge == 5) { $InternalCharge = -1; last CHARGE;}
    if ($MDLCharge == 6) { $InternalCharge = -2; last CHARGE;}
    if ($MDLCharge == 7) { $InternalCharge = -3; last CHARGE;}
    # All other MDL charge values, including 4 corresponding to "doublet radical",
    # are assigned internal value of 0.
    $InternalCharge = 0;
    if ($MDLCharge != 4) {
      carp "Warning: MDLChargeToInternalCharge: MDL charge value, $MDLCharge, is not supported: An internal charge value, 0, has been assigned...";
    }
  }
  return $InternalCharge;
}

# Map internal charge used by MayaChemTools to MDL charge value used in SD and MOL files.
#
sub InternalChargeToMDLCharge {
  my($InternalCharge) = @_;
  my($MDLCharge);

  CHARGE: {
    if ($InternalCharge == 3) { $MDLCharge = 1; last CHARGE;}
    if ($InternalCharge == 2) { $MDLCharge = 2; last CHARGE;}
    if ($InternalCharge == 1) { $MDLCharge = 3; last CHARGE;}
    if ($InternalCharge == -1) { $MDLCharge = 5; last CHARGE;}
    if ($InternalCharge == -2) { $MDLCharge = 6; last CHARGE;}
    if ($InternalCharge == -3) { $MDLCharge = 7; last CHARGE;}
    # All other MDL charge values, including 4 corresponding to "doublet radical",
    # are assigned internal value of 0.
    $MDLCharge = 0;
  }
  return $MDLCharge;
}

# Ctab lines: Bond block
#
# Format: 111222tttsssxxxrrrccc
#
# 111: First atom number.
# 222: Second atom number.
# ttt: Bond type. 1 = Single, 2 = Double, 3 = Triple, 4 = Aromatic, 5 = Single or Double,
#      6 = Single or Aromatic, 7 = Double or Aromatic, 8 = Any
# sss: Bond stereo. Single bonds: 0 = not stereo, 1 = Up, 4 = Either, 6 = Down,
#      Double bonds: 0 = Use x-, y-, z-coords from atom block to determine cis or trans,
#      3 = Cis or trans (either) double bond
# xxx: Not used
# rrr: Bond topology. 0 = Either, 1 = Ring, 2 = Chain
# ccc: Reacting center status. 0 = unmarked, 1 = a center, -1 = not a center,
#      Additional: 2 = no change,4 = bond made/broken, 8 = bond order changes 12 = 4+8
#      (both made/broken and changes); 5 = (4 + 1), 9 = (8 + 1), and 13 = (12 + 1) are also possible
#
sub ParseCmpdBondLine {
  my($Line) = @_;
  my($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo);

  ($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo) = map {s/ //g; $_} unpack("A3A3A3A3", $Line);
  return ($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo);
}

# Map MDL bond type value used in SD and MOL files to internal bond order  and bond types
# values used by MayaChemTools...
#
sub MDLBondTypeToInternalBondOrder {
  my($MDLBondType) = @_;
  my($InternalBondOrder, $InternalBondType);

  $InternalBondType = '';

  BONDTYPE: {
    if ($MDLBondType == 1) { $InternalBondOrder = 1; $InternalBondType = 'Single'; last BONDTYPE;}
    if ($MDLBondType == 2) { $InternalBondOrder = 2; $InternalBondType = 'Double'; last BONDTYPE;}
    if ($MDLBondType == 3) { $InternalBondOrder = 3; $InternalBondType = 'Triple'; last BONDTYPE;}
    if ($MDLBondType == 4) { $InternalBondOrder = 1.5; $InternalBondType = 'Aromatic'; last BONDTYPE;} # Aromatic
    if ($MDLBondType == 5) { $InternalBondOrder = 1; $InternalBondType = 'SingleOrDouble'; last BONDTYPE;} # Aromatic
    if ($MDLBondType == 6) { $InternalBondOrder = 1; $InternalBondType = 'SingleOrAromatic'; last BONDTYPE;} # Aromatic
    if ($MDLBondType == 7) { $InternalBondOrder = 2; $InternalBondType = 'DoubleOrAromatic'; last BONDTYPE;} # Aromatic
    if ($MDLBondType == 8) { $InternalBondOrder = 1; $InternalBondType = 'Any'; last BONDTYPE;} # Aromatic
    #
    # Although MDL aromatic bond values are used for query only and explicit Kekule bond order
    # values must be assigned, internal value of 1.5 is allowed to indicate aromatic bond orders.
    #
    # All other MDL bond type values -  5 = Single or Double, 6 = Single or Aromatic, 7 = Double or Aromatic,
    # 8 = Any - are also assigned appropriate internal value of 1: These are meant to be used for
    # structure queries by MDL products.
    #
    $InternalBondOrder = 1;
    $InternalBondType = 'Single';

    carp "Warning: MDLBondTypeToInternalBondOrder: MDL bond type value, $MDLBondType, is not supported: An internal bond order value, 0, has been assigned...";
  }
  return ($InternalBondOrder, $InternalBondType);
}

# Map internal bond order  and bond type values used by MayaChemTools to MDL bond type value used
# in SD and MOL files...
#
sub InternalBondOrderToMDLBondType {
  my($InternalBondOrder, $InternalBondType) = @_;
  my($MDLBondType);

  BONDTYPE: {
    if ($InternalBondOrder == 1) {
      if ($InternalBondType =~ /^SingleOrDouble$/i) {
	$MDLBondType = 5;
      }
      elsif ($InternalBondType =~ /^SingleOrAromatic$/i) {
	$MDLBondType = 6;
      }
      elsif ($InternalBondType =~ /^Any$/i) {
	$MDLBondType = 8;
      }
      else {
	$MDLBondType = 1;
      }
      $MDLBondType = 1;
      last BONDTYPE;
    }
    if ($InternalBondOrder == 2) {
      if ($InternalBondType =~ /^DoubleOrAromatic$/i) {
	$MDLBondType = 7;
      }
      else {
	$MDLBondType = 2;
      }
      last BONDTYPE;
    }
    if ($InternalBondOrder == 3) { $MDLBondType = 3; last BONDTYPE;}
    if ($InternalBondOrder == 1.5) { $MDLBondType = 4; last BONDTYPE;}
    if ($InternalBondType =~ /^Any$/i) { $MDLBondType = 8; last BONDTYPE;}

    $MDLBondType = 1;

    carp "Warning: InternalBondOrderToMDLBondType: Internal bond order and type values, $InternalBondOrder and $InternalBondType, don't match any valid MDL bond type: MDL bond type value, 1, has been assigned...";
  }
  return $MDLBondType;
}

# Third line: Comments - A blank line is also allowed.
sub ParseCmpdCommentsLine {
  my($Line) = @_;
  my($Comments);

  $Comments = unpack("A80", $Line);

  return ($Comments);
}

# Map MDL bond stereo value used in SD and MOL files to internal bond stereochemistry values used by MayaChemTools...
#
sub MDLBondStereoToInternalBondStereochemistry {
  my($MDLBondStereo) = @_;
  my($InternalBondStereo);

  $InternalBondStereo = '';

  BONDSTEREO: {
    if ($MDLBondStereo == 1) { $InternalBondStereo = 'Up'; last BONDSTEREO;}
    if ($MDLBondStereo == 4) { $InternalBondStereo = 'UpOrDown'; last BONDSTEREO;}
    if ($MDLBondStereo == 6) { $InternalBondStereo = 'Down'; last BONDSTEREO;}
    if ($MDLBondStereo == 3) { $InternalBondStereo = 'CisOrTrans'; last BONDSTEREO;}
    if ($MDLBondStereo == 0) { $InternalBondStereo = 'None'; last BONDSTEREO;}

    $InternalBondStereo = '';
    carp "Warning: MDLBondStereoToInternalBondType: MDL bond stereo value, $MDLBondStereo, is not supported: It has been ignored and bond order would be used to determine bond type...";
  }
  return $InternalBondStereo;
}

# Map internal bond stereochemistry values used by MayaChemTools to MDL bond stereo value used in SD and MOL files...
#
sub InternalBondStereochemistryToMDLBondStereo {
  my($InternalBondStereo) = @_;
  my($MDLBondStereo);

  $MDLBondStereo = 0;

  BONDSTEREO: {
    if ($InternalBondStereo =~ /^Up$/i) { $MDLBondStereo = 1; last BONDSTEREO;}
    if ($InternalBondStereo =~ /^UpOrDown$/i) { $MDLBondStereo = 4; last BONDSTEREO;}
    if ($InternalBondStereo =~ /^Down$/) { $MDLBondStereo = 6; last BONDSTEREO;}
    if ($InternalBondStereo =~ /^CisOrTrans$/) { $MDLBondStereo = 3; last BONDSTEREO;}

    $MDLBondStereo = 0;
  }
  return $MDLBondStereo;
}

# Fourth line: Counts
#
# Format: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
#
# aaa: number of atoms; bbb: number of bonds; lll: number of atom lists; fff: (obsolete)
# ccc: chiral flag: 0=not chiral, 1=chiral; sss: number of stext entries; xxx,rrr,ppp,iii:
# (obsolete); mmm: number of lines of additional properties, including the M END line, No
# longer supported, default is set to 999; vvvvvv: version

sub ParseCmpdCountsLine {
  my($Line) = @_;
  my($AtomCount, $BondCount, $ChiralFlag, $PropertyCount, $Version);

  if (length($Line) >= 39) {
    ($AtomCount, $BondCount, $ChiralFlag, $PropertyCount, $Version) = unpack("A3A3x3x3A3x3x3x3x3x3A3A6", $Line);
  }
  elsif (length($Line) >= 15) {
    ($PropertyCount, $Version) = ("999", "v2000");
    ($AtomCount, $BondCount, $ChiralFlag) = unpack("A3A3x3x3A3", $Line);
  }
  else {
    ($ChiralFlag, $PropertyCount, $Version) = ("0", "999", "v2000");
    ($AtomCount, $BondCount) = unpack("A3A3", $Line);
  }

  if ($Version =~ /V3000/i) {
    # Current version of MayaChemTools modules and classes for processing MDL MOL and SD don't support
    # V3000. So instead of relying on callers, just exit with an error to disable any processing of V3000
    # format.
    croak "Error: SDFileUtil::ParseCmpdCountsLine: The Extended Connection Table (V3000) format in MDL MOL and SD files is not supported by the current release of MayaChemTools...";
  }

  return ($AtomCount, $BondCount, $ChiralFlag, $PropertyCount, $Version);
}

# Second line: Misc info
#
# Format: IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR
#         A2A8      A10       A2I2A10       A12         A6
# User's first and last initials (I), program name (P), date/time (M/D/Y,H:m),
# dimensional codes - 2D or 3D (d),scaling factors (S, s), energy (E) if modeling program input,
# internal registry number (R) if input through MDL form. A blank line is also allowed.
sub ParseCmpdMiscInfoLine {
  my($Line) = @_;
  my($UserInitial, $ProgramName, $Date, $Code, $ScalingFactor1, $ScalingFactor2, $Energy, $RegistryNum);

  ($UserInitial, $ProgramName, $Date, $Code, $ScalingFactor1, $ScalingFactor2, $Energy, $RegistryNum) = unpack("A2A8A10A2A2A10A12A6", $Line);
  return ($UserInitial, $ProgramName, $Date, $Code, $ScalingFactor1, $ScalingFactor2, $Energy, $RegistryNum);
}

# First line: Molecule name. This line is unformatted, but like all other lines in a
# molfile may not extend beyond column 80. A blank line is also allowed.
sub ParseCmpdMolNameLine {
  my($Line) = @_;
  my($MolName);

  $MolName = unpack("A80", $Line);

  return ($MolName);
}

# Parse atom alias property line in CTAB generic properties block.
#
# Atom alias property line format:
#
# A  aaa
# x...
#
#    aaa: Atom number
#    x: Atom alias in next line
#
sub ParseCmpdAtomAliasPropertyLine {
  my($Line, $NextLine) = @_;
  my($Label, $AtomNumber, $AtomAlias);

  ($Label, $AtomNumber) = split(' ', $Line);
  $AtomAlias = $NextLine;

  if (!$AtomAlias) {
    carp "Warning: _ParseCmpdAtomAliasPropertyLine: No atom alias value specified on the line following atom alias property line...";
  }

  return ($AtomNumber, $AtomAlias);
}

# Parse charge property line in CTAB generic properties block.
#
# Charge property line format:
#
# M  CHGnn8 aaa vvv ...
#
#    nn8: Number of value pairs. Maximum of 8 pairs allowed.
#    aaa: Atom number
#    vvv: -15 to +15. Default of 0 = uncharged atom. When present, this property supersedes
#    all charge and radical values in the atom block, forcing a 0 charge on all atoms not
#    listed in an M  CHG or M  RAD line.
#
sub ParseCmpdChargePropertyLine {
  my($Line) = @_;

  return _ParseCmpdGenericPropertyLine('Charge', $Line);
}


# Parse isotope property line in CTAB generic properties block.
#
# Isoptope property line format:
#
# M  ISOnn8 aaa vvv ...
#
#    nn8: Number of value paris. Maximum of 8 pairs allowed.
#    aaa: Atom number
#    vvv: Absolute mass of the atom isotope as a positive integer. When present, this property
#    supersedes all isotope values in the atom block. Default (no entry) means natural
#    abundance. The difference between this absolute mass value and the natural
#    abundance value specified in the PTABLE.DAT file must be within the range of -18
#    to +12
#
# Notes:
#  . Values correspond to mass numbers...
#
sub ParseCmpdIsotopePropertyLine {
  my($Line) = @_;

  return _ParseCmpdGenericPropertyLine('Isotope', $Line);
}

# Parse radical property line in CTAB generic properties block.
#
# Radical property line format:
#
# M  RADnn8 aaa vvv ...
#
#    nn8: Number of value paris. Maximum of 8 pairs allowed.
#    aaa: Atom number
#    vvv: Default of 0 = no radical, 1 = singlet, 2 = doublet, 3 = triplet . When
#    present, this property supersedes all charge and radical values in the atom block,
#    forcing a 0 (zero) charge and radical on all atoms not listed in an M  CHG or
#    M  RAD line.
#
sub ParseCmpdRadicalPropertyLine {
  my($Line) = @_;

  return _ParseCmpdGenericPropertyLine('Radical', $Line);
}

# Map MDL radical stereo value used in SD and MOL files to internal spin multiplicity values used by MayaChemTools...
#
sub MDLRadicalToInternalSpinMultiplicity {
  my($MDLRadical) = @_;
  my($InternalSpinMultiplicity);

  $InternalSpinMultiplicity = '';

  SPINMULTIPLICITY: {
    if ($MDLRadical == 0) { $InternalSpinMultiplicity = 0; last SPINMULTIPLICITY;}
    if ($MDLRadical == 1) { $InternalSpinMultiplicity = 1; last SPINMULTIPLICITY;}
    if ($MDLRadical == 2) { $InternalSpinMultiplicity = 2; last SPINMULTIPLICITY;}
    if ($MDLRadical == 3) { $InternalSpinMultiplicity = 3; last SPINMULTIPLICITY;}
    $InternalSpinMultiplicity = '';
    carp "Warning: MDLRadicalToInternalSpinMultiplicity: MDL radical value, $MDLRadical, specifed on line M  RAD is not supported...";
  }
  return $InternalSpinMultiplicity;
}

# Map internal spin multiplicity values used by MayaChemTools to MDL radical stereo value used in SD and MOL files...
#
sub InternalSpinMultiplicityToMDLRadical {
  my($InternalSpinMultiplicity) = @_;
  my($MDLRadical);

  $MDLRadical = 0;

  SPINMULTIPLICITY: {
    if ($InternalSpinMultiplicity == 1) { $MDLRadical = 1; last SPINMULTIPLICITY;}
    if ($InternalSpinMultiplicity == 2) { $MDLRadical = 2; last SPINMULTIPLICITY;}
    if ($InternalSpinMultiplicity == 3) { $MDLRadical = 3; last SPINMULTIPLICITY;}
    $MDLRadical = 0;
  }
  return $MDLRadical;
}

# Process generic CTAB property line...
sub _ParseCmpdGenericPropertyLine {
  my($PropertyName, $Line) = @_;

  my($Label, $PropertyLabel, $ValuesCount, $ValuePairsCount, @ValuePairs);

  @ValuePairs = ();
  ($Label, $PropertyLabel, $ValuesCount, @ValuePairs) = split(' ', $Line);
  $ValuePairsCount = (scalar @ValuePairs)/2;
  if ($ValuesCount != $ValuePairsCount) {
    carp "Warning: _ParseCmpdGenericPropertyLine: Number of atom number and $PropertyName value paris specified on $Label $PropertyLabel property line, $ValuePairsCount, does not match expected value of $ValuesCount...";
  }

  return (@ValuePairs);
}

# Generic CTAB property lines for charge, istope and radical properties...
#
sub _GenerateCmpdGenericPropertyLines {
  my($PropertyName, $PropertyValuePairsRef) = @_;
  my($Index, $PropertyLabel, $Line, $PropertyCount, $AtomNum, $PropertyValue, @PropertyLines);

  @PropertyLines = ();
  NAME: {
    if ($PropertyName =~ /^Charge$/i) { $PropertyLabel = "M  CHG"; last NAME; }
    if ($PropertyName =~ /^Isotope$/i) { $PropertyLabel = "M  ISO"; last NAME; }
    if ($PropertyName =~ /^Radical$/i) { $PropertyLabel = "M  RAD"; last NAME; }
    carp "Warning: _GenerateCmpdGenericPropertyLines: Unknown property name, $PropertyName, specified...";
    return @PropertyLines;
  }

  # A maximum of 8 property pair values allowed per line...
  $PropertyCount = 0;
  $Line = '';
  for ($Index = 0; $Index < $#{$PropertyValuePairsRef}; $Index += 2) {
    if ($PropertyCount > 8) {
      # Setup property line...
      $Line = "${PropertyLabel}  8${Line}";
      push @PropertyLines, $Line;

      $PropertyCount = 0;
      $Line = '';
    }
    $PropertyCount++;
    $AtomNum = $PropertyValuePairsRef->[$Index];
    $PropertyValue = $PropertyValuePairsRef->[$Index + 1];
    $Line .= sprintf " %3i %3i", $AtomNum, $PropertyValue;
  }
  if ($Line) {
    $Line = "${PropertyLabel}  ${PropertyCount}${Line}";
    push @PropertyLines, $Line;
  }
  return @PropertyLines;
}

#
# Read compound data into a string and return its value
sub ReadCmpdString {
  my($SDFileRef) = @_;
  my($CmpdString);

  $CmpdString = "";
  LINE: while (defined($_ = <$SDFileRef>)) {
    # Change Windows and Mac new line char to UNIX...
    s/(\r\n)|(\r)/\n/g;

    if (/^\$\$\$\$/) {
      # Take out any new line char at the end by explicitly removing it instead of using
      # chomp, which might not always work correctly on files generated on a system
      # with a value of input line separator different from the current system...
      s/\n$//g;

      # Doesn't hurt to chomp...
      chomp;

      $CmpdString .=  $_;
      last LINE;
    }
    else {
      $CmpdString .=  $_;
    }
  }
  return $CmpdString;
}

# Find out the number of fragements in the compounds. And for the compound with
# more than one fragment, remove all the others besides the largest one.
sub WashCmpd {
  my($CmpdLines) = @_;
  my($WashedCmpdString, $FragmentCount, $Fragments);

  $WashedCmpdString = "";
  ($FragmentCount, $Fragments) = GetCmpdFragments($CmpdLines);
  if ($FragmentCount > 1) {
    # Go over the compound data for the largest fragment including property
    # data...
    my (@AllFragments, @LargestFragment, %LargestFragmentAtoms, @WashedCmpdLines, $Index, $LineIndex, $AtomCount, $BondCount, $NewAtomCount, $NewBondCount, $FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo, $FirstNewAtomNum, $SecondNewAtomNum, $AtomNum, $ChiralFlag, $BondLine, $MENDLineIndex, $Line, $Value, @ValuePairs, @NewValuePairs, $ValuePairIndex, $NewAtomNum, @NewPropertyLines);

    @AllFragments = (); @LargestFragment = ();
    %LargestFragmentAtoms = ();
    @AllFragments = split "\n", $Fragments;
    @LargestFragment = split " ", $AllFragments[0];
    for $Index (0 .. $#LargestFragment) {
      # Map old atom numbers to new atom numbers as the fragment atom numbers are sorted
      # from lowest to highest old atom numbers...
      $LargestFragmentAtoms{$LargestFragment[$Index]} = $Index + 1;
    }
    @WashedCmpdLines = ();
    push @WashedCmpdLines, @$CmpdLines[0], @$CmpdLines[1], @$CmpdLines[2], @$CmpdLines[3];
    ($AtomCount, $BondCount, $ChiralFlag) = ParseCmpdCountsLine(@$CmpdLines[3]);
    $NewAtomCount = @LargestFragment;
    $NewBondCount = 0;
    $AtomNum = 0;
    # Retrieve the largest fragment atom lines...
    for ($LineIndex = 4; $LineIndex < (4 + $AtomCount); $LineIndex++) {
      $AtomNum++;
      if ($LargestFragmentAtoms{$AtomNum}) {
	push @WashedCmpdLines, @$CmpdLines[$LineIndex];
      }
    }
    # Retrieve the largest fragment bond lines...
    for ($LineIndex = 4 + $AtomCount; $LineIndex < (4 + $AtomCount + $BondCount); $LineIndex++) {
      ($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo) = ParseCmpdBondLine(@$CmpdLines[$LineIndex]);
      if ($LargestFragmentAtoms{$FirstAtomNum} && $LargestFragmentAtoms{$SecondAtomNum}) {
	$NewBondCount++;
	# Set up bond line with new atom number mapping...
	$FirstNewAtomNum =  $LargestFragmentAtoms{$FirstAtomNum};
	$SecondNewAtomNum =  $LargestFragmentAtoms{$SecondAtomNum};
	$BondLine = GenerateCmpdBondLine($FirstNewAtomNum, $SecondNewAtomNum, $BondType, $BondStereo);
	push @WashedCmpdLines, $BondLine;
      }
    }
    # Get property lines for CHG, ISO and RAD label and map the old atom numbers to new
    # atom numners; Others, property lines before M  END line are skipped as atom numbers for
    # other properties might not valid anymore...
    #
    $MENDLineIndex = $LineIndex;
    LINE: for ($LineIndex = (4 + $AtomCount + $BondCount); $LineIndex < @$CmpdLines; $LineIndex++) {
      $Line = @$CmpdLines[$LineIndex];
      if ($Line =~ /^M  END/i) {
	push @WashedCmpdLines, "M  END";
	$MENDLineIndex = $LineIndex;
	last LINE;
      }

      @ValuePairs = ();
      if ($Line =~ /^M  CHG/i) {
	@ValuePairs = ParseCmpdChargePropertyLine($Line);
      }
      elsif ($Line =~ /^M  RAD/i) {
	@ValuePairs = ParseCmpdRadicalPropertyLine($Line);
      }
      elsif ($Line =~ /^M  ISO/i) {
	@ValuePairs = ParseCmpdIsotopePropertyLine($Line);
      }
      elsif ($Line =~ /^A  /i) {
	my($NextLine);
	$LineIndex++;
	$NextLine = @$CmpdLines[$LineIndex];
	@ValuePairs = ParseCmpdAtomAliasPropertyLine($Line, $NextLine);
      }
      else {
	next LINE;
      }

      if (!@ValuePairs) {
	next LINE;
      }

      # Collect values for valid atom numbers with mapping to new atom numbers...
      @NewValuePairs = ();
      VALUEINDEX: for ($ValuePairIndex = 0; $ValuePairIndex < $#ValuePairs; $ValuePairIndex += 2) {
	$AtomNum = $ValuePairs[$ValuePairIndex]; $Value = $ValuePairs[$ValuePairIndex + 1];
	if (!exists $LargestFragmentAtoms{$AtomNum}) {
	  next VALUEINDEX;
	}
	$NewAtomNum = $LargestFragmentAtoms{$AtomNum};
	push @NewValuePairs, ($NewAtomNum, $Value)
      }
      if (!@NewValuePairs) {
	next LINE;
      }
      @NewPropertyLines = ();
      if ($Line =~ /^M  CHG/i) {
	@NewPropertyLines = GenerateCmpdChargePropertyLines(\@NewValuePairs);
      }
      elsif ($Line =~ /^M  RAD/i) {
	@NewPropertyLines = GenerateCmpdRadicalPropertyLines(\@NewValuePairs);
      }
      elsif ($Line =~ /^M  ISO/i) {
	@NewPropertyLines = GenerateCmpdIsotopePropertyLines(\@NewValuePairs);
      }
      elsif ($Line =~ /^A  /i) {
	@NewPropertyLines = GenerateCmpdAtomAliasPropertyLines(\@NewValuePairs);
      }
      push @WashedCmpdLines, @NewPropertyLines;
    }

    # Retrieve rest of the data label and value property data...
    for ($LineIndex = (1 + $MENDLineIndex); $LineIndex < @$CmpdLines; $LineIndex++) {
      push @WashedCmpdLines, @$CmpdLines[$LineIndex];
    }
    # Update atom and bond count line...
    $WashedCmpdLines[3] = GenerateCmpdCountsLine($NewAtomCount, $NewBondCount, $ChiralFlag);

    $WashedCmpdString = join "\n", @WashedCmpdLines;
  }
  return ($FragmentCount, $Fragments, $WashedCmpdString);
}

1;

__END__

=head1 NAME

SDFileUtil

=head1 SYNOPSIS

use SDFileUtil ;

use SDFileUtil qw(:all);

=head1 DESCRIPTION

B<SDFileUtil> module provides the following functions:

GenerateCmpdAtomAliasPropertyLines, GenerateCmpdAtomLine, GenerateCmpdBondLine,
GenerateCmpdChargePropertyLines, GenerateCmpdCommentsLine, GenerateCmpdCountsLine,
GenerateCmpdDataHeaderLabelsAndValuesLines, GenerateCmpdIsotopePropertyLines,
GenerateCmpdMiscInfoLine, GenerateCmpdMolNameLine,
GenerateCmpdRadicalPropertyLines, GenerateEmptyCtabBlockLines,
GenerateMiscLineDateStamp, GetAllAndCommonCmpdDataHeaderLabels,
GetCmpdDataHeaderLabels, GetCmpdDataHeaderLabelsAndValues, GetCmpdFragments,
GetCtabLinesCount, GetInvalidAtomNumbers, GetUnknownAtoms,
InternalBondOrderToMDLBondType, InternalBondStereochemistryToMDLBondStereo,
InternalChargeToMDLCharge, InternalSpinMultiplicityToMDLRadical, IsCmpd2D,
IsCmpd3D, MDLBondStereoToInternalBondStereochemistry,
MDLBondTypeToInternalBondOrder, MDLChargeToInternalCharge,
MDLRadicalToInternalSpinMultiplicity, ParseCmpdAtomAliasPropertyLine,
ParseCmpdAtomLine, ParseCmpdBondLine, ParseCmpdChargePropertyLine,
ParseCmpdCommentsLine, ParseCmpdCountsLine, ParseCmpdIsotopePropertyLine,
ParseCmpdMiscInfoLine, ParseCmpdMolNameLine, ParseCmpdRadicalPropertyLine,
ReadCmpdString, RemoveCmpdDataHeaderLabelAndValue, WashCmpd

=head1 METHODS

=over 4

=item B<GenerateCmpdAtomAliasPropertyLines>

    @Lines = GenerateCmpdAtomAliasPropertyLines($AliasValuePairsRef);

Returns a formatted atom alias property lines corresponding to successive pairs
of atom number and alias values specified by a refernce to an array. Two lines
are generate for each atom number and alias value pairs: First line - A  <AtomNum>;
Second line:<AtomAlias>.

=item B<GenerateCmpdAtomLine>

    $Line = GenerateCmpdAtomLine($AtomSymbol, $AtomX, $AtomY,
            $AtomZ, [$MassDifference, $Charge, $StereoParity]);

Returns a formatted atom data line containing all the input values.

=item B<GenerateCmpdBondLine>

    $Line = GenerateCmpdBondLine($FirstAtomNum, $SecondAtomNum,
            $BondType, [$BondStereo]);

Returns a formatted bond data line containing all the input values.

=item B<GenerateCmpdChargePropertyLines>

    @Lines = GenerateCmpdChargePropertyLines($ChargeValuePairsRef);

Returns a formatted M  CHG property lines corresponding to successive pairs of
atom number and charge values specified by a refernce to an array.

=item B<GenerateCmpdCommentsLine>

    $Line = GenerateCmpdCommentsLine($Comments);

Returns a formatted comments data line.

=item B<GenerateCmpdCountsLine>

    $Line = GenerateCmpdCountsLine($AtomCount, $BondCount,
            $ChiralFlag, [$PropertyCount, $Version]);

Returns a formatted line containing all the input values. The default values of 999
and  V2000 are used for I<PropertyCount> and I<Version>.

=item B<GenerateCmpdDataHeaderLabelsAndValuesLines>

    @Lines = GenerateCmpdDataHeaderLabelsAndValuesLines(
             $DataHeaderLabelsRef, $DataHeaderLabelsAndValuesRef,
             [$SortDataLabels]);

Returns formatted data lines containing header label and values lines corresponding to
all data header labels in array reference I<DataHeaderLabelsRef> with values in hash
reference I<DataHeaderLabelsAndValuesRef>. By default, data header labels are
not sorted and correspond to the label order in array reference I<DataHeaderLabelsRef>.

=item B<GenerateCmpdIsotopePropertyLines>

    @Lines = GenerateCmpdIsotopePropertyLines($IsotopeValuePairsRef);

Returns a formatted M ISO property lines corresponding to successive pairs of
atom number and isotope values specified by a refernce to an array.

=item B<GenerateCmpdMiscInfoLine>

    $Line = GenerateCmpdMiscInfoLine([$ProgramName, $UserInitial,
            $Code]);

Returns a formatted line containing specified user initial, program name, date and code.
Default values are: I<ProgramName - MayaChem; UserInitial - NULL; Code - 2D>.

=item B<GenerateCmpdMolNameLine>

    $Line = GenerateCmpdMolNameLine($MolName);

Returns a formatted molecule name data line.

=item B<GenerateCmpdRadicalPropertyLines>

    @Lines = GenerateCmpdRadicalPropertyLines($RadicalValuePairsRef);

Returns a formatted M  CHG property lines corresponding to successive pairs of
atom number and multiplicity values specified by a refernce to an array.

=item B<GenerateEmptyCtabBlockLines>

    $Lines = GenerateCmpdMiscInfoLine([$Date]);

Returns formatted lines representing empty CTAB block.

=item B<GenerateMiscLineDateStamp>

    $Line = GenerateMiscLineDateStamp();

Returns date stamp for misc line.

=item B<GetAllAndCommonCmpdDataHeaderLabels>

    ($CmpdCount, $DataFieldLabelsArrayRef,
       $CommonDataFieldLabelsArrayRef) =
          GetAllAndCommonCmpdDataHeaderLabels(\*SDFILE);

Returns number of comopunds, a reference to an array containing all unique data header
label and a reference to an array containing common data field labels for all compounds
in SD file.

=item B<GetCmpdDataHeaderLabels>

    (@Labels) = GetCmpdDataHeaderLabels(\@CmpdLines);

Returns an array containg data header labels for a compound

=item B<GetCmpdDataHeaderLabelsAndValues>

    (%DataValues) = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

Returns a hash conating data header labes and values for a compound.

=item B<GetCmpdFragments>

    ($FragmentCount, $FragmentString) = GetCmpdFragments(\@CmpLines);

Figures out the number of disconnected fragments and return their values along
with the atom numbers in a string delimited by new line character. Fragment data
in B<FragmentString> is sorted on based on its size.

=item B<GetCtabLinesCount>

    $CtabLinesCount = GetCtabLinesCount(\@CmpdLines);

Returns number of lines present between the 4th line and the line containg "M END".

=item B<GetInvalidAtomNumbers>

    ($InvalidAtomNumbersCount, $InvalidAtomNumbers, $InvalidAtomNumberLines) =
       GetInvalidAtomNumbers(\@CmpdLines);

Returns a list of values containing information about invalid atom numbers present
in block or atom property lines.

=item B<GetUnknownAtoms>

    ($UnknownAtomCount, $UnknownAtoms, $UnknownAtomLines) =
       GetUnknownAtoms(\@CmpdLines);

Returns a list of values containing information about atoms which contain special element
symbols not present in the periodic table.

=item B<InternalBondOrderToMDLBondType>

    $MDLBondType = InternalBondOrderToMDLBondType($InternalBondOrder);

Returns value of I<MDLBondType> corresponding to I<InternalBondOrder>.

    InternalBondOrder  MDLBondType

     1                  1
     2                  2
     3                  3
     1.5                4

=item B<InternalBondStereochemistryToMDLBondStereo>

    $MDLBondStereo = InternalBondStereochemistryToMDLBondStereo(
                     $InternalBondStereo);

Returns value of I<MDLBondStereo> corresponding to I<InternalBondStereo> using following
mapping:

    InternalBondStereo  MDLBondStereo

     Up          1
     UpOrDown    4
     Down        6
     CisOrTrans  3
     Other       0

=item B<InternalChargeToMDLCharge>

    $MDLCharge = InternalChargeToMDLCharge($InternalCharge);

Returns value of I<MDLCharge> corresponding to I<InternalCharge> using following
mapping:

    InternalCharge  MDLCharge

     3               1
     2               2
     1               3
    -1              5
    -2              6
    -3              7

=item B<InternalSpinMultiplicityToMDLRadical>

    $MDLRadical = InternalSpinMultiplicityToMDLRadical(
                  $InternalSpinMultiplicity);

Returns value of I<MDLRadical> corresponding to I<InternalSpinMultiplicity>. These
value are equivalent.

=item B<MDLBondStereoToInternalBondType>

    $InternalBondType = MDLBondStereoToInternalBondType($MDLBondStereo);

Returns value of I<InternalBondType> corresponding to I<MDLBondStereo> using
mapping shown for B<InternalBondTypeToMDLBondStereo> function.

=item B<IsCmpd2D>

    $Status = IsCmpd2D();

Returns 1 or 0 based on whether z-coordinate of any atom is non-zero.

=item B<IsCmpd3D>

    $Status = IsCmpd3D();

Returns 1 or 0 based on whether z-coordinate of any atom is non-zero.

=item B<MDLBondStereoToInternalBondStereochemistry>

    $InternalBondStereo = MDLBondStereoToInternalBondStereochemistry(
                          $MDLBondStereo);

Returns value of I<InternalBondStereo> corresponding to I<MDLBondStereo> using
mapping shown for B<InternalBondStereochemistryToMDLBondStereo> function.

=item B<MDLBondTypeToInternalBondOrder>

    $InternalBondOrder = MDLBondTypeToInternalBondOrder($MDLBondType);

Returns value of I<InternalBondOrder> corresponding to I<MDLBondType> using
mapping shown for B<InternalBondOrderToMDLBondType> function.

=item B<MDLChargeToInternalCharge>

    $InternalCharge = MDLChargeToInternalCharge($MDLCharge);

Returns value of I<$InternalCharge> corresponding to I<MDLCharge> using
mapping shown for B<InternalChargeToMDLCharge> function.

=item B<MDLRadicalToInternalSpinMultiplicity>

    $InternalSpinMultiplicity = MDLRadicalToInternalSpinMultiplicity(
                                $MDLRadical);

Returns value of I<InternalSpinMultiplicity> corresponding to I<MDLRadical>. These
value are equivalent.

=item B<ParseCmpdAtomAliasPropertyLine>

    @AtomNumAndValuePairs = ParseCmpdAtomAliasPropertyLine(
                            $CurrentLine, $NexLine);

Parses atom alias propery lines in CTAB generic properties block and returns an array
with successive pairs of values corresponding to atom number and its alias.

=item B<ParseCmpdAtomLine>

    ($AtomSymbol, $AtomX, $AtomY, $AtomZ, $MassDifference, $Charge,
       $StereoParity) = ParseCmpdAtomLine($AtomDataLine);

Parses compound data line containing atom information and returns a list
of values.

=item B<ParseCmpdBondLine>

    ($FirstAtomNum, $SecondAtomNum, $BondType) =
       ParseCmpdBondLine($BondDataLine);

Parses compound data line containing bond information and returns a list of
values.

=item B<ParseCmpdCommentsLine>

    $Comments = ParseCmpdCommentsLine($CommentsDataLine);

Returns the comment string.

=item B<ParseCmpdChargePropertyLine>

    @AtomNumAndValuePairs = ParseCmpdChargePropertyLine(
                            $ChargeDataLine);

Parses charge propery line in CTAB generic properties block and returns an array
with successive pairs of values corresponding to atom number and its charge.

=item B<ParseCmpdCountsLine>

    ($AtomCount, $BondCount, $ChiralFlag, $PropertyCount, $Version) =
       ParseCmpdCountsLine(\@CountDataLines);

Returns a list of values containing count information.

=item B<ParseCmpdMiscInfoLine>

    ($UserInitial, $ProgramName, $Date, $Code, $ScalingFactor1, $ScalingFactor2,
       $Energy, $RegistryNum) =  ParseCmpdMiscInfoLine($Line);

Returns a list of values containing miscellaneous information.

=item B<ParseCmpdIsotopePropertyLine>

    @AtomNumAndValuePairs = ParseCmpdIsotopePropertyLine(
                            $IsotopeDataLine);

Parses isotopic propery line in CTAB generic properties block and returns an array
with successive pairs of values corresponding to atom number and absolute mass of
atom isotope.

=item B<ParseCmpdMolNameLine>

    $MolName = ParseCmpdMolNameLine($Line);

Returns a string containing molecule name.

=item B<ParseCmpdRadicalPropertyLine>

    @AtomNumAndValuePairs = ParseCmpdRadicalPropertyLine(
                            $RadicalDataLine);

Parses radical propery line in CTAB generic properties block and returns an array
with successive pairs of values corresponding to atom number and radical number
value.

=item B<RemoveCmpdDataHeaderLabelAndValue>

    $NewCmpdString = RemoveCmpdDataHeaderLabelAndValue($CmpdString,
                                                       $DataHeaderLabel);

Returns a B<NewCmpdString> after removing  I<DataHeaderLabel> along with its
value from I<CmpdString>.

=item B<ReadCmpdString>

    $CmpdString = ReadCmpdString(\*SDFILEHANDLE);

Returns a string containing all the data lines for the next available compound
in an already open file indicated by SDFILEHANDLE. A NULL string is returned
on EOF.

=item B<WashCmpd>

    ($FragmentCount, $Fragments, $WashedCmpdString) = 
       WashCmpd(\@CmpdLines);

Figures out the number of disconnected fragments and return their values along
with the atom numbers in a string delimited by new line character. Fragment data
in B<FragmentString> is sorted on based on its size.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

TextUtil.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
