package AminoAcids;
#
# File: AminoAcids.pm
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
use Carp;
use Text::ParseWords;
use TextUtil;
use FileUtil;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GetAminoAcids GetAminoAcidPropertiesData GetAminoAcidPropertiesNames IsAminoAcid IsAminoAcidProperty);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

#
# Load amino acids data...
#
my(%AminoAcidDataMap, %AminoAcidThreeLetterCodeMap, %AminoAcidOneLetterCodeMap, %AminoAcidNameMap, @AminoAcidPropertyNames, %AminoAcidPropertyNamesMap, );
_LoadAminoAcidsData();

#
# Get a list of all known amino acids as one of these values:
# one letter code, three letter code, or amino acid name...
#
sub GetAminoAcids {
  my($NameType, $ThreeLetterCode, $Name, @AminoAcidNames, %AminoAcidNamesMap);

  $NameType = 'ThreeLetterCode';
  if (@_ >= 1) {
    ($NameType) = @_;
  }

  # Collect names...
  %AminoAcidNamesMap = ();
  for $ThreeLetterCode (keys %AminoAcidDataMap) {
    NAME : {
      if ($NameType =~ /^OneLetterCode$/i) {$Name = $AminoAcidDataMap{$ThreeLetterCode}{OneLetterCode}; last NAME; }
      if ($NameType =~ /^AminoAcid$/i) {$Name = $AminoAcidDataMap{$ThreeLetterCode}{AminoAcid}; last NAME; }
      $Name = $ThreeLetterCode;
    }
    $AminoAcidNamesMap{$Name} = $Name;
  }

  # Sort 'em out
  @AminoAcidNames = ();
  for $Name (sort keys %AminoAcidNamesMap) {
    push @AminoAcidNames, $Name;
  }

  return (wantarray ? @AminoAcidNames : \@AminoAcidNames);
}


#
# Get all available properties data for an amino acid using any of these symbols:
# three letter code; one letter code; name.
#
# A reference to a hash array is returned with keys and values representing property
# name and its values respectively.
#
sub GetAminoAcidPropertiesData {
  my($AminoAcidID) = @_;
  my($ThreeLetterCode);

  if ($ThreeLetterCode = _ValidateAminoAcidID($AminoAcidID)) {
    return \%{$AminoAcidDataMap{$ThreeLetterCode}};
  }
  else {
    return undef;
  }
}

#
# Get names of all available amino acid properties. A reference to  an array containing
# names of all available properties is returned.
#
sub GetAminoAcidPropertiesNames {
  my($Mode);
  my($PropertyName, @PropertyNames);

  $Mode = 'ByGroup';
  if (@_ == 1) {
    ($Mode) = @_;
  }

  @PropertyNames = ();
  if ($Mode =~ /^Alphabetical$/i) {
    my($PropertyName);
    # ThreeLetterCode, OneLetterCode, and AminoAcid are always listed first...
    push @PropertyNames, qw(ThreeLetterCode OneLetterCode AminoAcid);
    for $PropertyName (sort keys %AminoAcidPropertyNamesMap) {
      if ($PropertyName !~ /^(ThreeLetterCode|OneLetterCode|AminoAcid)$/) {
	push @PropertyNames, $PropertyName;
      }
    }
  }
  else {
    push @PropertyNames, @AminoAcidPropertyNames;
  }
  return (wantarray ? @PropertyNames : \@PropertyNames);
}

#
# Is it a known amino acid? Input is either an one/three letter code or a name.
#
sub IsAminoAcid {
  my($AminoAcidID) = @_;
  my($Status);

  $Status = (_ValidateAminoAcidID($AminoAcidID)) ? 1 : 0;

  return $Status;
}


#
# Is it an available amino acid property?
#
sub IsAminoAcidProperty {
  my($PropertyName) = @_;
  my($Status);

  $Status = (exists($AminoAcidPropertyNamesMap{$PropertyName})) ? 1 : 0;

  return $Status;
}

#
# Implents GetAminoAcid<PropertyName> for a valid proprty name.
#
sub AUTOLOAD {
  my($AminoAcidID) = @_;
  my($FunctionName, $PropertyName, $PropertyValue, $ThreeLetterCode);

  $PropertyValue = undef;

  use vars qw($AUTOLOAD);
  $FunctionName = $AUTOLOAD;
  $FunctionName =~ s/.*:://;

  # Only Get<PropertyName> functions are supported...
  if ($FunctionName !~ /^Get/) {
    croak "Error: Function, AminoAcid::$FunctionName, is not supported by AUTOLOAD in AminoAcid module: Only Get<PropertyName> functions are implemented...";
  }

  $PropertyName = $FunctionName;
  $PropertyName =~  s/^GetAminoAcid//;
  if (!exists $AminoAcidPropertyNamesMap{$PropertyName}) {
    croak "Error: Function, AminoAcid::$FunctionName, is not supported by AUTOLOAD in AminoAcid module: Unknown amino acid property name, $PropertyName, specified...";
  }

  if (!($ThreeLetterCode = _ValidateAminoAcidID($AminoAcidID))) {
    return undef;
  }
  $PropertyValue = $AminoAcidDataMap{$ThreeLetterCode}{$PropertyName};
  return $PropertyValue;
}


#
# Load AminoAcidsData.csv files from <MayaChemTools>/lib directory...
#
sub _LoadAminoAcidsData {
  my($AminoAcidsDataFile, $MayaChemToolsLibDir);

  $MayaChemToolsLibDir = GetMayaChemToolsLibDirName();

  $AminoAcidsDataFile =  "$MayaChemToolsLibDir" . "/data/AminoAcidsData.csv";

  if (! -e "$AminoAcidsDataFile") {
    croak "Error: MayaChemTools package file, $AminoAcidsDataFile, is missing: Possible installation problems...";
  }

  _LoadData($AminoAcidsDataFile);
}

#
# Load AminoAcidsData.csv file from <MayaChemTools>/lib directory...
#
sub _LoadData {
  my($AminoAcidsDataFile) = @_;

  %AminoAcidDataMap = ();
  @AminoAcidPropertyNames = ();
  %AminoAcidPropertyNamesMap = ();
  %AminoAcidThreeLetterCodeMap = ();
  %AminoAcidOneLetterCodeMap = ();
  %AminoAcidNameMap = ();

  # Load property data for all amino acids...
  #
  # File Format:
  #"ThreeLetterCode","OneLetterCode","AminoAcid","AcidicBasic","PolarNonpolar","Charged","Aromatic","HydrophobicHydophilic","IsoelectricPoint","pKCOOH","pKNH3+","MolecularWeight","MolecularWeightMinusH2O(18.01524)","ExactMass","ExactMassMinusH2O(18.01056)","vanderWaalsVolume","%AccessibleResidues","%BuriedResidues","AlphaHelixChouAndFasman","AlphaHelixDeleageAndRoux","AlphaHelixLevitt","AminoAcidsComposition","AminoAcidsCompositionInSwissProt","AntiparallelBetaStrand","AverageAreaBuried","AverageFlexibility","BetaSheetChouAndFasman","BetaSheetDeleageAndRoux","BetaSheetLevitt","BetaTurnChouAndFasman","BetaTurnDeleageAndRoux","BetaTurnLevitt","Bulkiness","CoilDeleageAndRoux","HPLCHFBARetention","HPLCRetentionAtpH2.1","HPLCRetentionAtpH7.4","HPLCTFARetention","HydrophobicityAbrahamAndLeo","HydrophobicityBlack","HydrophobicityBullAndBreese","HydrophobicityChothia","HydrophobicityEisenbergAndOthers","HydrophobicityFauchereAndOthers","HydrophobicityGuy","HydrophobicityHPLCAtpH3.4Cowan","HydrophobicityHPLCAtpH7.5Cowan","HydrophobicityHPLCParkerAndOthers","HydrophobicityHPLCWilsonAndOthers","HydrophobicityHoppAndWoods","HydrophobicityJanin","HydrophobicityKyteAndDoolittle","HydrophobicityManavalanAndOthers","HydrophobicityMiyazawaAndOthers","HydrophobicityOMHSweetAndOthers","HydrophobicityRaoAndArgos","HydrophobicityRfMobility","HydrophobicityRoseAndOthers","HydrophobicityRoseman","HydrophobicityWellingAndOthers","HydrophobicityWolfendenAndOthers","MolecularWeight","NumberOfCodons","ParallelBetaStrand","PolarityGrantham","PolarityZimmerman","RatioHeteroEndToSide","RecognitionFactors","Refractivity","RelativeMutability","TotalBetaStrand","LinearStructure","LinearStructureAtpH7.4"
  #
  #
  my($ThreeLetterCode, $OneLetterCode, $AminoAcidName, $Line, $NumOfCols, $InDelim, $Index, $Name, $Value, $Units, @LineWords, @ColLabels);

  $InDelim = "\,";
  open AMINOACIDSDATAFILE, "$AminoAcidsDataFile" or croak "Couldn't open $AminoAcidsDataFile: $! ...";

  # Skip lines up to column labels...
  LINE: while ($Line = GetTextLine(\*AMINOACIDSDATAFILE)) {
    if ($Line !~ /^#/) {
      last LINE;
    }
  }
  @ColLabels= quotewords($InDelim, 0, $Line);
  $NumOfCols = @ColLabels;

  # Extract property names from column labels...
  @AminoAcidPropertyNames = ();
  for $Index (0 .. $#ColLabels) {
    $Name = $ColLabels[$Index];
    push @AminoAcidPropertyNames, $Name;

    # Store property names...
    $AminoAcidPropertyNamesMap{$Name} = $Name;
  }

  # Process amino acid data...
  LINE: while ($Line = GetTextLine(\*AMINOACIDSDATAFILE)) {
    if ($Line =~ /^#/) {
      next LINE;
    }
    @LineWords = ();
    @LineWords = quotewords($InDelim, 0, $Line);
    if (@LineWords != $NumOfCols) {
      croak "Error: The number of data fields, @LineWords, in $AminoAcidsDataFile must be $NumOfCols.\nLine: $Line...";
    }
    $ThreeLetterCode = $LineWords[0]; $OneLetterCode = $LineWords[1]; $AminoAcidName = $LineWords[3];
    if (exists $AminoAcidDataMap{$ThreeLetterCode}) {
      carp "Warning: Ignoring data for amino acid $ThreeLetterCode: It has already been loaded.\nLine: $Line....";
      next LINE;
    }

    # Store all the values...
    %{$AminoAcidDataMap{$ThreeLetterCode}} = ();
    for $Index (0 .. $#LineWords) {
      $Name = $AminoAcidPropertyNames[$Index];
      $Value = $LineWords[$Index];
      $AminoAcidDataMap{$ThreeLetterCode}{$Name} = $Value;
    }
  }
  close AMINOACIDSDATAFILE;

  # Setup one letter and amino acid name maps...
  _SetupAminoAcidIDMap();
}


#
# Setup lowercase three/one letter code and name maps pointing
# to three letter code as show in data file.
#
sub _SetupAminoAcidIDMap {
  my($ThreeLetterCode, $OneLetterCode, $AminoAcidName);

  %AminoAcidThreeLetterCodeMap = ();
  %AminoAcidOneLetterCodeMap = ();
  %AminoAcidNameMap = ();

  for $ThreeLetterCode (keys %AminoAcidDataMap) {
    $OneLetterCode = $AminoAcidDataMap{$ThreeLetterCode}{OneLetterCode};
    $AminoAcidName = $AminoAcidDataMap{$ThreeLetterCode}{AminoAcid};

    $AminoAcidThreeLetterCodeMap{lc($ThreeLetterCode)} = $ThreeLetterCode;
    $AminoAcidOneLetterCodeMap{lc($OneLetterCode)} = $ThreeLetterCode;
    $AminoAcidNameMap{lc($AminoAcidName)} = $ThreeLetterCode;
  }
}

# Validate amino acid ID...
sub _ValidateAminoAcidID {
  my($AminoAcidID) = @_;
  my($ThreeLetterCode);


  if (length($AminoAcidID) == 3) {
    if (! exists $AminoAcidThreeLetterCodeMap{lc($AminoAcidID)}) {
      return undef;
    }
    $ThreeLetterCode = $AminoAcidThreeLetterCodeMap{lc($AminoAcidID)};
  }
  elsif (length($AminoAcidID) == 1) {
    if (! exists $AminoAcidOneLetterCodeMap{lc($AminoAcidID)}) {
      return undef;
    }
    $ThreeLetterCode = $AminoAcidOneLetterCodeMap{lc($AminoAcidID)};
  }
  else {
    if (! exists $AminoAcidNameMap{lc($AminoAcidID)}) {
      return undef;
    }
    $ThreeLetterCode = $AminoAcidNameMap{lc($AminoAcidID)};
  }
  return $ThreeLetterCode;
}


1;

__END__

=head1 NAME

AminoAcids

=head1 SYNOPSIS

use AminoAcids;

use AminoAcids qw(:all);

=head1 DESCRIPTION

B<AminoAcids> module provides the following functions:

GetAminoAcidPropertiesData, GetAminoAcidPropertiesNames, GetAminoAcid<PropertyName>,
GetAminoAcids, IsAminoAcid, IsAminoAcidProperty

=head1 FUNCTIONS

=over 4

=item B<GetAminoAcidPropertiesData>

    $DataHashRef = GetAminoAcidPropertiesData($AminoAcidID);

Returns a reference to hash containing property names and values for a specified
amino acid.

=item B<GetAminoAcidPropertiesNames>

    @Names = GetAminoAcidPropertiesNames([$Mode]);
    $NamesRef = GetAminoAcidPropertiesNames([$Mode]);

Returns an array or a reference to an array containing names of amino acids
properties. Order of amino acids properties is controlled by optional parameter
I<Mode>. Possible values for I<Mode>: I<Alphabetical or  ByGroup>; Default: I<ByGroup>

=item B<GetAminoAcidPropertyName>

    $Value = GetAminoAcid<PropertyName>($AminoAcidID);

Returns amino acid property value for a specified amino acid. These functions are
not defined in this modules; these are implemented on the fly using Perl's AUTOLOAD
funcion. Here is the list of known amino acids I<property names>: DNACodons, RNACodons,
AcidicBasic, PolarNonpolar, Charged, Aromatic, HydrophobicHydophilic, IsoelectricPoint,
pKCOOH, pKNH3+, ChemicalFormula, MolecularWeight, ExactMass, ChemicalFormulaMinusH2O,
MolecularWeightMinusH2O(18.01524), ExactMassMinusH2O(18.01056), vanderWaalsVolume,
%AccessibleResidues, %BuriedResidues, AlphaHelixChouAndFasman,
AlphaHelixDeleageAndRoux, AlphaHelixLevitt, AminoAcidsComposition,
AminoAcidsCompositionInSwissProt, AntiparallelBetaStrand, AverageAreaBuried, AverageFlexibility,
BetaSheetChouAndFasman, BetaSheetDeleageAndRoux, BetaSheetLevitt,
BetaTurnChouAndFasman, BetaTurnDeleageAndRoux, BetaTurnLevitt, Bulkiness,
CoilDeleageAndRoux, HPLCHFBARetention, HPLCRetentionAtpH2.1, HPLCRetentionAtpH7.4,
HPLCTFARetention, HydrophobicityAbrahamAndLeo, HydrophobicityBlack,
HydrophobicityBullAndBreese, HydrophobicityChothia, HydrophobicityEisenbergAndOthers,
HydrophobicityFauchereAndOthers, HydrophobicityGuy, HydrophobicityHPLCAtpH3.4Cowan,
HydrophobicityHPLCAtpH7.5Cowan, HydrophobicityHPLCParkerAndOthers,
HydrophobicityHPLCWilsonAndOthers, HydrophobicityHoppAndWoods, HydrophobicityJanin,
HydrophobicityKyteAndDoolittle, HydrophobicityManavalanAndOthers,
HydrophobicityMiyazawaAndOthers, HydrophobicityOMHSweetAndOthers,
HydrophobicityRaoAndArgos, HydrophobicityRfMobility, HydrophobicityRoseAndOthers,
HydrophobicityRoseman, HydrophobicityWellingAndOthers, HydrophobicityWolfendenAndOthers,
ParallelBetaStrand, PolarityGrantham, PolarityZimmerman, RatioHeteroEndToSide,
RecognitionFactors, Refractivity, RelativeMutability, TotalBetaStrand, LinearStructure,
LinearStructureAtpH7.4

=item B<GetAminoAcids>

    $NamesRef = GetAminoAcids([$NameType]);
    (@Names) = GetAminoAcids([$NameType]);

Returns an array or a reference to an array containing names of amino acids
as one letter code, three letter code, or amino acid name controlled by optional
parameter $NameType. By default, amino acids names are returned as three
letter code. Possible values for I<NameType>: I<ThreeLetterCode, OneLetterCode, or
AminoAcid>.

=item B<IsAminoAcid>

    $Status = IsAminoAcid($AminoAcidID);

Returns a flag indicating whether or not its a known amino acid ID.

=item B<IsAminoAcidProperty>

    $Status = IsAminoAcid($PropertyName);

Returns a flag indicating whether or not its a known amino acid property name.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

NucleicAcids.pm, PeriodicTable.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
