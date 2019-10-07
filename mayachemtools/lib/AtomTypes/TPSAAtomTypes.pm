package AtomTypes::TPSAAtomTypes;
#
# File: TPSAAtomTypes.pm
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
use Exporter;
use Scalar::Util ();
use Text::ParseWords;
use AtomTypes::AtomTypes;
use Molecule;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(AtomTypes::AtomTypes Exporter);
@EXPORT = qw(GetTPSAAtomTypesData GetAllPossibleTPSAAtomTypes);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, %TPSAAtomTypesDataMap);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyTPSAAtomTypes';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeTPSAAtomTypes();

  $This->_InitializeTPSAAtomTypesProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Initialize the data hash. It'll be loaded on demand later...
  %TPSAAtomTypesDataMap = ();
}

# Initialize object data...
#
sub _InitializeTPSAAtomTypes {
  my($This) = @_;

  # Type of AtomTypes...
  $This->{Type} = 'TPSA';

  # Besides polar atoms - N, O, P, S - no TPSA atom types are assigned to any other
  # atoms.
  #
  # By default, TPSA atom types are not assigned to Phosphorus and Sulfur atoms.
  #
  $This->{IgnorePhosphorus} = 0;
  $This->{IgnoreSulfur} = 0;

  return $This;
}

# Initialize object properties...
#
sub _InitializeTPSAAtomTypesProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  # Make sure molecule object was specified...
  if (!exists $NamesAndValues{Molecule}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying molecule...";
  }

  return $This;
}

# Get TPSA atom types and associated data loaded from TPSA data file as
# a reference to hash with the following hash data format:
#
# @{$TPSAAtomTypesDataMap{AtomTypes}} - Array of all possible atom types for all atoms
# @{$TPSAAtomTypesDataMap->{ColLabels}} - Array of column labels
# %{$TPSAAtomTypesDataMap->{DataCol<Num>}} - Hash keys pair: <DataCol<Num>, AtomType>
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetTPSAAtomTypesData {

  # Make sure data is loaded...
  _CheckAndLoadTPSAAtomTypesData();

  return \%TPSAAtomTypesDataMap;
}

# Get all possible TPSA atom types atoms as an array reference...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAllPossibleTPSAAtomTypes {
  return _GetAllPossibleTPSAAtomTypes();
}

# Are all atoms types successfully assigned?
#
# Notes:
#   . Dynamic checking of atom types assignment for atoms eliminates the need
#     to check and synchronize valid atom types during SetAtomType.
#   . Base class method is overrided to check atom assignment to nitrogen and
#     oxygen atom with optional check for phosphorus and sulfur atoms.
#
sub IsAtomTypesAssignmentSuccessful {
  my($This) = @_;
  my($Atom, $AtomType);

  ATOM: for $Atom ($This->{Molecule}->GetAtoms()) {
    if (!($Atom->IsNitrogen() || $Atom->IsOxygen() ||
        ($Atom->IsPhosphorus() && !$This->{IgnorePhosphorus}) ||
        ($Atom->IsSulfur() && !$This->{IgnoreSulfur}))) {
      next ATOM;
    }
    $AtomType = $This->GetAtomType($Atom);
    if ($AtomType =~ /^None$/i) {
      return 0;
    }
  }

  return 1;
}

# Get all possible TPSA atom types as an array reference...
#
sub _GetAllPossibleTPSAAtomTypes {
  my($TPSAAtomTypesDataRef);

  $TPSAAtomTypesDataRef = GetTPSAAtomTypesData();

  return \@{$TPSAAtomTypesDataRef->{AtomTypes}};
}

# Assign Topological Polar Surface Area (TPSA) atom types [ Ref 90-91 ] to Nitrogen and Oxygen
# atoms with optional assignment to Phosphorus and Sulfur atoms.
#
# Notes:
#     o Number of atom type symbols for:
#         o N: 27
#         o O: 7
#         o P: 5
#         o S: 8
#
sub AssignAtomTypes {
  my($This) = @_;
  my($Atom, $AtomType);

  ATOM: for $Atom ($This->GetMolecule()->GetAtoms()) {
    $AtomType = $This->_GetAtomType($Atom);
    $This->SetAtomType($Atom, $AtomType);
  }

  return $This;
}

# Get TPSA atom type for atom...
#
sub _GetAtomType {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOM: {
    if ($Atom->IsNitrogen()) {
      $AtomType = $This->_GetAtomTypeForNitrogen($Atom);
      last ATOM;
    }

    if ($Atom->IsOxygen()) {
      $AtomType = $This->_GetAtomTypeForOxygen($Atom);
      last ATOM;
    }

    if ($Atom->IsPhosphorus() && !$This->{IgnorePhosphorus}) {
      $AtomType = $This->_GetAtomTypeForPhosphorus($Atom);
      last ATOM;
    }

    if ($Atom->IsSulfur() && !$This->{IgnoreSulfur}) {
      $AtomType = $This->_GetAtomTypeForSulfur($Atom);
      last ATOM;
    }
    $AtomType = 'None';
  }

  return $AtomType;
}


# Get TPSA atom type for Nitrogen atom...
#
# 27 AtomTypeSymbols for element N:
#
# AtomTypeSymbol - SMARTS - Comments
# N1 - '[N](-*)(-*)-*'
# N2 - '[N](-*)=*'
# N3 - '[N]#*'
# N4 - '[N](-*)(=*)=*' - As in nitro group
# N5 - '[N](=*)#*' - Middle nitrogen in azide group
# N6 - '[N]1(-*)-*-*-1' - Atom in a 3 membered ring
# N7 - '[NH](-*)-*'
# N8 - '[NH]1-*-*-1' - Atom in a 3 membered ring
# N9 - '[NH]=*'
# N10 - '[NH2]-*'
# N11 - '[N+](-*)(-*)(-*)-*'
# N12 - '[N+](-*)(-*)=*'
# N13 - '[N+](-*)#*' - Nitrogen in isocyano group
# N14 - '[NH+](-*)(-*)-*'
# N15 - '[NH+](-*)=*'
# N16 - '[NH2+](-*)-*'
# N17 - '[NH2+]=*'
# N18 - '[NH3+]-*'
# N19 - '[n](:*):*'
# N20 - '[n](:*)(:*):*'
# N21 - '[n](-*)(:*):*'
# N22 - '[n](=*)(:*):*' - As in pyridine N-oxide
# N23 - '[nH](:*):*'
# N24 - '[n+](:*)(:*):*'
# N25 - '[n+](-*)(:*):*'
# N26 - '[nH+](:*):*'
# N - '[#7]' - Any other Nitrogen; Contribution: 30.5 - X*8.2 + H*1.5 or 0.0 for negative value
#
sub _GetAtomTypeForNitrogen {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'None';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {

    # Aromatic Nitrogens...
    if ($Atom->IsAromatic()) {
      $AtomType = $This->_GetAtomTypeForAromaticNitrogen($Atom);
      last ATOMTYPE;
    }

    # Only single bonds...
    if ($NumOfPiBonds == 0) {
      $AtomType = $This->_GetAtomTypeForNitrogenWithOnlySigmaBonds($Atom);
      last ATOMTYPE;
    }

    # One double bond...
    if ($NumOfPiBonds == 1) {
      $AtomType = $This->_GetAtomTypeForNitrogenWithOnePiBond($Atom);
      last ATOMTYPE;
    }

    # One triple bond or two double bonds...
    if ($NumOfPiBonds == 2) {
      $AtomType = $This->_GetAtomTypeForNitrogenWithTwoPiBonds($Atom);
      last ATOMTYPE;
    }

    # One triple bond and a double bond...
    if ($NumOfPiBonds == 3) {
      $AtomType = $This->_GetAtomTypeForNitrogenWithThreePiBonds($Atom);
      last ATOMTYPE;
    }

    $AtomType = 'N';
  }
  return $AtomType;
}

# Get TPSA atom type for Oxygen atom...
#
# AtomTypeSymbol - SMARTS - Comments
# O1 - '[O](-*)-*'
# O2 - '[O]1-*-*-1' - Atom in a 3 membered ring
# O3 - '[O]=*'
# O4 - '[OH]-*'
# O5 - '[O-]-*'
# O6 - '[o](:*):*'
# O - '[#8]' - Any other Oxygen; Contribution: 28.5 - X*8.6 + H*1.5 or 0.0 for negative value
#
sub _GetAtomTypeForOxygen {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # O6 - '[o](:*):*'
    if ($This->_IsO6Oxygen($Atom)) {
      $AtomType = 'O6';
      last ATOMTYPE;
    }

    # O3 - '[O]=*'
    if ($This->_IsO3Oxygen($Atom)) {
      $AtomType = 'O3';
      last ATOMTYPE;
    }

    # O4 - '[OH]-*'
    if ($This->_IsO4Oxygen($Atom)) {
      $AtomType = 'O4';
      last ATOMTYPE;
    }

    # O5 - '[O-]-*'
    if ($This->_IsO5Oxygen($Atom)) {
      $AtomType = 'O5';
      last ATOMTYPE;
    }

    # O2 - '[O]1-*-*-1' - Atom in a 3 membered ring
    if ($This->_IsO2Oxygen($Atom)) {
      $AtomType = 'O2';
      last ATOMTYPE;
    }

    # O1 - '[O](-*)-*'
    if ($This->_IsO1Oxygen($Atom)) {
      $AtomType = 'O1';
      last ATOMTYPE;
    }

    # Any other Oxygen...
    $AtomType = 'O';
  }

  return $AtomType;
}

# Get TPSA atom type for Phosphorus atom...
#
# 4 AtomTypeSymbols for element P:
#
# AtomTypeSymbol - SMARTS - Comments
# P1 - '[P](-*)(-*)-*'
# P2 - '[P](-*)=*'
# P3 - '[P](-*)(-*)(-*)=*'
# P4 - '[PH](-*)(-*)=*'
# P - '[#15]' - Any other Sulfur
#
sub _GetAtomTypeForPhosphorus {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # P1 - '[P](-*)(-*)-*'
    if ($This->_IsP1Phosphorus($Atom)) {
      $AtomType = 'P1';
      last ATOMTYPE;
    }

    # P2 - '[P](-*)=*'
    if ($This->_IsP2Phosphorus($Atom)) {
      $AtomType = 'P2';
      last ATOMTYPE;
    }

    # P3 - '[P](-*)(-*)(-*)=*'
    if ($This->_IsP3Phosphorus($Atom)) {
      $AtomType = 'P3';
      last ATOMTYPE;
    }

    # P4 - '[PH](-*)(-*)=*'
    if ($This->_IsP4Phosphorus($Atom)) {
      $AtomType = 'P4';
      last ATOMTYPE;
    }

    # Any other Phosphorus...
    $AtomType = 'P';
  }

  return $AtomType;
}

# Get TPSA atom type for Sulfur atom...
#
# 7 AtomTypeSymbols for element S:
#
# AtomTypeSymbol - SMARTS - Comments
# S1 - '[S](-*)-*'
# S2 - '[S]=*'
# S3 - '[S](-*)(-*)=*'
# S4 - '[S](-*)(-*)(=*)=*'
# S5 - '[SH]-*'
# S6 - '[s](:*):*'
# S - '[#16]' - Any other Phosphorus
#
sub _GetAtomTypeForSulfur {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # S6 - '[s](:*):*'
    if ($This->_IsS6Sulfur($Atom)) {
      $AtomType = 'S6';
      last ATOMTYPE;
    }

    # S4 - '[S](-*)(-*)(=*)=*'
    if ($This->_IsS4Sulfur($Atom)) {
      $AtomType = 'S4';
      last ATOMTYPE;
    }

    # S1 - '[S](-*)-*'
    if ($This->_IsS1Sulfur($Atom)) {
      $AtomType = 'S1';
      last ATOMTYPE;
    }

    # S2 - '[S]=*'
    if ($This->_IsS2Sulfur($Atom)) {
      $AtomType = 'S2';
      last ATOMTYPE;
    }

    # S3 - '[S](-*)(-*)=*'
    if ($This->_IsS3Sulfur($Atom)) {
      $AtomType = 'S3';
      last ATOMTYPE;
    }

    # S5 - '[SH]-*'
    if ($This->_IsS5Sulfur($Atom)) {
      $AtomType = 'S5';
      last ATOMTYPE;
    }

    # Any other Sulfur...
    $AtomType = 'S';
  }

  return $AtomType;
}

# Get TPSA atom type for aromatic Nitrogen...
#
sub _GetAtomTypeForAromaticNitrogen {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # N19 - '[n](:*):*'
    if ($This->_IsN19Nitrogen($Atom)) {
      $AtomType = 'N19';
      last ATOMTYPE;
    }

    # N20 - '[n](:*)(:*):*'
    if ($This->_IsN20Nitrogen($Atom)) {
      $AtomType = 'N20';
      last ATOMTYPE;
    }

    # N21 - '[n](-*)(:*):*'
    if ($This->_IsN21Nitrogen($Atom)) {
      $AtomType = 'N21';
      last ATOMTYPE;
    }

    # N22 - '[n](=*)(:*):*' - As in pyridine N-oxide
    if ($This->_IsN22Nitrogen($Atom)) {
      $AtomType = 'N22';
      last ATOMTYPE;
    }

    # N23 - '[nH](:*):*'
    if ($This->_IsN23Nitrogen($Atom)) {
      $AtomType = 'N23';
      last ATOMTYPE;
    }

    # N24 - '[n+](:*)(:*):*'
    if ($This->_IsN24Nitrogen($Atom)) {
      $AtomType = 'N24';
      last ATOMTYPE;
    }

    # N25 - '[n+](-*)(:*):*'
    if ($This->_IsN25Nitrogen($Atom)) {
      $AtomType = 'N25';
      last ATOMTYPE;
    }

    # N26 - '[nH+](:*):*'
    if ($This->_IsN26Nitrogen($Atom)) {
      $AtomType = 'N26';
      last ATOMTYPE;
    }

    $AtomType = 'N';
  }

  return $AtomType;
}

# Get TPSA atom type for Nitrogen with only sigma bonds...
#
sub _GetAtomTypeForNitrogenWithOnlySigmaBonds {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

   # N6 - '[N]1(-*)-*-*-1' - Atom in a 3 membered ring
    if ($This->_IsN6Nitrogen($Atom)) {
      $AtomType = 'N6';
      last ATOMTYPE;
    }

    # N1 - '[N](-*)(-*)-*'
    if ($This->_IsN1Nitrogen($Atom)) {
      $AtomType = 'N1';
      last ATOMTYPE;
    }

    # N8 - '[NH]1-*-*-1' - Atom in a 3 membered ring
    if ($This->_IsN8Nitrogen($Atom)) {
      $AtomType = 'N8';
      last ATOMTYPE;
    }

    # N7 - '[NH](-*)-*'
    if ($This->_IsN7Nitrogen($Atom)) {
      $AtomType = 'N7';
      last ATOMTYPE;
    }

    # N10 - '[NH2]-*'
    if ($This->_IsN10Nitrogen($Atom)) {
      $AtomType = 'N10';
      last ATOMTYPE;
    }

    # N11 - '[N+](-*)(-*)(-*)-*'
    if ($This->_IsN11Nitrogen($Atom)) {
      $AtomType = 'N11';
      last ATOMTYPE;
    }

    # N14 - '[NH+](-*)(-*)-*'
    if ($This->_IsN14Nitrogen($Atom)) {
      $AtomType = 'N14';
      last ATOMTYPE;
    }

    # N16 - '[NH2+](-*)-*'
    if ($This->_IsN16Nitrogen($Atom)) {
      $AtomType = 'N16';
      last ATOMTYPE;
    }

    # N18 - '[NH3+]-*'
    if ($This->_IsN18Nitrogen($Atom)) {
      $AtomType = 'N18';
      last ATOMTYPE;
    }

    $AtomType = 'N';
  }

  return $AtomType;
}

# Get TPSA atom type for Nitrogen with one pi bonds...
#
sub _GetAtomTypeForNitrogenWithOnePiBond {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # N2 - '[N](-*)=*'
    if ($This->_IsN2Nitrogen($Atom)) {
      $AtomType = 'N2';
      last ATOMTYPE;
    }

    # N9 - '[NH]=*'
    if ($This->_IsN9Nitrogen($Atom)) {
      $AtomType = 'N9';
      last ATOMTYPE;
    }

    # N12 - '[N+](-*)(-*)=*'
    if ($This->_IsN12Nitrogen($Atom)) {
      $AtomType = 'N12';
      last ATOMTYPE;
    }

    # N15 - '[NH+](-*)=*'
    if ($This->_IsN15Nitrogen($Atom)) {
      $AtomType = 'N15';
      last ATOMTYPE;
    }

    # N17 - '[NH2+]=*'
    if ($This->_IsN17Nitrogen($Atom)) {
      $AtomType = 'N17';
      last ATOMTYPE;
    }

    $AtomType = 'N';
  }

  return $AtomType;
}

# Get TPSA atom type for Nitrogen with two pi bonds...
#
sub _GetAtomTypeForNitrogenWithTwoPiBonds {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # N3 - '[N]#*'
    if ($This->_IsN3Nitrogen($Atom)) {
      $AtomType = 'N3';
      last ATOMTYPE;
    }

    # N4 - '[N](-*)(=*)=*' - As in nitro group
    if ($This->_IsN4Nitrogen($Atom)) {
      $AtomType = 'N4';
      last ATOMTYPE;
    }

    # N13 - '[N+](-*)#*'- Nitrogen in isocyano group
    if ($This->_IsN13Nitrogen($Atom)) {
      $AtomType = 'N13';
      last ATOMTYPE;
    }

    $AtomType = 'N';
  }

  return $AtomType;
}

# Get TPSA atom type for Nitrogen with three pi bonds...
#
sub _GetAtomTypeForNitrogenWithThreePiBonds {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # N5 - '[N](=*)#*' - Middle nitrogen in azide group
    if ($This->_IsN5Nitrogen($Atom)) {
      $AtomType = 'N5';
      last ATOMTYPE;
    }

    $AtomType = 'N';
  }

  return $AtomType;
}

# N1 - '[N](-*)(-*)-*'
#
sub _IsN1Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.!RA3.X3.SB3.H0.FC0') ? 1 : 0;
}

# N2 - '[N](-*)=*'
#
sub _IsN2Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X2.SB1.DB1.H0.FC0') ? 1 : 0;
}

# N3 - '[N]#*'
#
sub _IsN3Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X1.TB1.H0.FC0') ? 1 : 0;
}

# N4 - '[N](-*)(=*)=*' - As in nitro group
#
sub _IsN4Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X3.SB1.DB2.H0.FC0') ? 1 : 0;
}

# N5 - '[N](=*)#*' - Middle nitrogen in azide group
#
sub _IsN5Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X2.DB1.TB1.H0.FC0') ? 1 : 0;
}

# N6 - '[N]1(-*)-*-*-1' - Atom in a 3 membered ring
#
sub _IsN6Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.RA3.X3.SB3.H0.FC0') ? 1 : 0;
}

# N7 - '[NH](-*)-*'
#
sub _IsN7Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.!RA3.X2.SB2.H1.FC0') ? 1 : 0;
}

# N8 - '[NH]1-*-*-1' - Atom in a 3 membered ring
#
sub _IsN8Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.RA3.X2.SB2.H1.FC0') ? 1 : 0;
}

# N9 - '[NH]=*'
#
sub _IsN9Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X1.DB1.H1.FC0') ? 1 : 0;
}

# N10 - '[NH2]-*'
#
sub _IsN10Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X1.SB1.H2.FC0') ? 1 : 0;
}

# N11 - '[N+](-*)(-*)(-*)-*'
#
sub _IsN11Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X4.SB4.H0.FC+1') ? 1 : 0;
}

# N12 - '[N+](-*)(-*)=*'
#
sub _IsN12Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X3.SB2.DB1.H0.FC+1') ? 1 : 0;
}

# N13 - '[N+](-*)#*'- Nitrogen in isocyano group
#
sub _IsN13Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X2.SB1.TB1.H0.FC+1') ? 1 : 0;
}

# N14 - '[NH+](-*)(-*)-*'
#
sub _IsN14Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X3.SB3.H1.FC+1') ? 1 : 0;
}

# N15 - '[NH+](-*)=*'
#
sub _IsN15Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X2.SB1.DB1.H1.FC+1') ? 1 : 0;
}

# N16 - '[NH2+](-*)-*'
#
sub _IsN16Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X2.SB2.H2.FC+1') ? 1 : 0;
}

# N17 - '[NH2+]=*'
#
sub _IsN17Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X1.DB1.H2.FC+1') ? 1 : 0;
}

# N18 - '[NH3+]-*'
#
sub _IsN18Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X1.SB1.H3.FC+1') ? 1 : 0;
}

# N19 - '[n](:*):*'
#
sub _IsN19Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.Ar.X2.AB2.H0.FC0') ? 1 : 0;
}

# N20 - '[n](:*)(:*):*'
#
sub _IsN20Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.Ar.X3.AB3.H0.FC0') ? 1 : 0;
}

# N21 - '[n](-*)(:*):*'
#
sub _IsN21Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.Ar.X3.AB2.H0.FC0', ['*', '*', '*'], [':', ':', '-']) ? 1 : 0;
}

# N22 - '[n](=*)(:*):*' - As in pyridine N-oxide
#
sub _IsN22Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.Ar.X3.AB2.H0.FC0', ['*', '*', '*'], [':', ':', '=']) ? 1 : 0;
}

# N23 - '[nH](:*):*'
#
sub _IsN23Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.Ar.X2.AB2.H1.FC0') ? 1 : 0;
}

# N24 - '[n+](:*)(:*):*'
#
sub _IsN24Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.Ar.X3.AB3.H0.FC+1') ? 1 : 0;
}

# N25 - '[n+](-*)(:*):*'
#
sub _IsN25Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.Ar.X3.AB2.H0.FC+1', ['*', '*', '*'], [':', ':', '-']) ? 1 : 0;
}

# N26 - '[nH+](:*):*'
#
sub _IsN26Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.Ar.X2.AB2.H1.FC+1') ? 1 : 0;
}

# O1 - '[O](-*)-*'
#
sub _IsO1Oxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.!RA.X2.SB2.H0.FC0') ? 1 : 0;
}

# O2 - '[O]1-*-*-1' - Atom in a 3 membered ring
#
sub _IsO2Oxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.RA3.X2.SB2.H0.FC0') ? 1 : 0;
}

# O3 - '[O]=*'
#
sub _IsO3Oxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.X1.DB1.H0.FC0') ? 1 : 0;
}

# O4 - '[OH]-*'
#
sub _IsO4Oxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.X1.SB1.H1.FC0') ? 1 : 0;
}

# O5 - '[O-]-*'
#
sub _IsO5Oxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.X1.SB1.H0.FC-1') ? 1 : 0;
}

# O6 - '[o](:*):*'
#
sub _IsO6Oxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.Ar.X2.AB2.H0.FC0') ? 1 : 0;
}

# P1 - '[P](-*)(-*)-*'
#
sub _IsP1Phosphorus {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('P.X3.SB3.H0.FC0') ? 1 : 0;
}

# P2 - '[P](-*)=*'
#
sub _IsP2Phosphorus {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('P.X2.SB1.DB1.H0.FC0') ? 1 : 0;
}

# P3 - '[P](-*)(-*)(-*)=*'
#
sub _IsP3Phosphorus {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('P.X4.SB3.DB1.H0.FC0') ? 1 : 0;
}

# P4 - '[PH](-*)(-*)=*'
#
sub _IsP4Phosphorus {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('P.X3.SB2.DB1.H1.FC0') ? 1 : 0;
}

# S1 - '[S](-*)-*'
#
sub _IsS1Sulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.X2.SB2.H0.FC0') ? 1 : 0;
}

# S2 - '[S]=*'
#
sub _IsS2Sulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.X1.DB1.H0.FC0') ? 1 : 0;
}

# S3 - '[S](-*)(-*)=*'
#
sub _IsS3Sulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.X3.SB2.DB1.H0.FC0') ? 1 : 0;
}

# S4 - '[S](-*)(-*)(=*)=*'
#
sub _IsS4Sulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.X4.SB2.DB2.H0.FC0') ? 1 : 0;
}

# S5 - '[SH]-*'
#
sub _IsS5Sulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.X1.SB1.H1.FC0') ? 1 : 0;
}

# S6 - '[s](:*):*'
#
sub _IsS6Sulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.Ar.X2.AB2.H0.FC0') ? 1 : 0;
}

# Return a string containg data for TPSAAtomTypes object...
#
sub StringifyTPSAAtomTypes {
  my($This) = @_;
  my($AtomTypesString);

  # Type of AtomTypes...
  $AtomTypesString = "AtomTypes: $This->{Type}; IgnorePhosphorus: " . ($This->{IgnorePhosphorus} ? "Yes" : "No") . "; IgnoreSulfur: " .  ($This->{IgnoreSulfur} ? "Yes" : "No");

  # Setup atom types information...
  my($AtomID, $AtomType, @AtomTypesInfo, %AssignedAtomTypes);

  @AtomTypesInfo = ();
  %AssignedAtomTypes = $This->GetAtomTypes();

  for $AtomID (sort { $a <=> $b } keys %AssignedAtomTypes) {
    $AtomType = $AssignedAtomTypes{$AtomID} ? $AssignedAtomTypes{$AtomID} : 'None';
    push @AtomTypesInfo, "$AtomID:$AtomType";
  }
  $AtomTypesString .= "; AtomIDs:AtomTypes: <" . TextUtil::JoinWords(\@AtomTypesInfo, ", ", 0) . ">";

  return $AtomTypesString;
}

# Is it a TPSAAtomTypes object?
sub _IsTPSAAtomTypes {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Check and load TPSA atom types data...
#
sub _CheckAndLoadTPSAAtomTypesData {

  # Is it already loaded?
  if (exists $TPSAAtomTypesDataMap{AtomTypes}) {
    return;
  }

  _LoadTPSAAtomTypesData();
}

# Load TPSA atom types data from the file assuming first column to be atom type symbol..
#
# Format:
#
# "AtomType","SMARTS","TPSAContribution","Comments"
# "N1","[N](-*)(-*)-*","3.24",""
# "N2","[N](-*)=*","12.36",""
#
sub _LoadTPSAAtomTypesData {
  my($AtomTypesDataFile, $MayaChemToolsLibDir);

  $MayaChemToolsLibDir = FileUtil::GetMayaChemToolsLibDirName();

  $AtomTypesDataFile =  "$MayaChemToolsLibDir" . "/data/TPSAAtomTypes.csv";
  if (! -e "$AtomTypesDataFile") {
    croak "Error: MayaChemTools package file, $AtomTypesDataFile, is missing: Possible installation problems...";
  }

  %TPSAAtomTypesDataMap = ();
  AtomTypes::AtomTypes::LoadAtomTypesData($AtomTypesDataFile, \%TPSAAtomTypesDataMap);
}

1;

__END__

=head1 NAME

TPSAAtomTypes

=head1 SYNOPSIS

use AtomTypes::TPSAAtomTypes;

use AtomTypes::TPSAAtomTypes qw(:all);

=head1 DESCRIPTION

B<TPSAAtomTypes> class provides the following methods:

new, AssignAtomTypes, GetAllPossibleTPSAAtomTypes, GetTPSAAtomTypesData,
StringifyTPSAAtomTypes

The following functions are available:

GetAllPossibleTPSAAtomTypes, GetTPSAAtomTypesData


B<TPSAAtomTypes> is derived from B<AtomTypes> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<TPSAAtomTypes>, B<AtomTypes> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The data file TPSAAomTypes.csv distributed with MayaChemTools release contains
all possible topological surface area (TPSA) [ Ref 90-91 ] atom types.

The atom type symbols assigned by MayaChemTools are not used in the original publication
and the numbers in atom symbols for element types simply correspond to their order of
appearance in Table 1 [ Ref 90 ].

Examples of TPSA atom types:

    N1, N2, N3, O1, O2, O3, S1, S2, P1, P2 and so on

=head2 METHODS

=over 4

=item B<new>

    $NewTPSAAtomTypes = new AtomTypes::TPSAAtomTypes(%NamesAndValues);

Using specified I<TPSAAtomTypes> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<TPSAAtomTypes>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'TPSA'
    IgnorePhosphorus = 0
    IgnoreSulfur = 0

Examples:

    $TPSAAtomTypes = new AtomTypes::TPSAAtomTypes(
                              'Molecule' => $Molecule,
                              'IgnorePhosphorus' => 0,
                              'IgnoreSulfur' => 0);

=item B<AssignAtomTypes>

    $TPSAAtomTypes->AssignAtomTypes();

Assigns TPSA atom types to all the atoms in a molecule and returns
I<TPSAAtomTypes>.

=item B<GetAllPossibleTPSAAtomTypes>

    $AllAtomTypesDataRef = $TPSAAtomTypes->
                           GetAllPossibleTPSAAtomTypes();
    $AllAtomTypesDataRef = AtomTypes::TPSAAtomTypes::
                           GetAllPossibleTPSAAtomTypes();

Returns all possible TPSA atom types corresponding to hydrogen and non-hydrogen
atoms as an array reference.

=item B<GetTPSAAtomTypesData>

    $AtomTypesDataMapRef = $TPSAAtomTypes->GetTPSAAtomTypesData();
    $AtomTypesDataMapRef = AtomTypes::TPSAAtomTypes::GetTPSAAtomTypesData();

Returns TPSA atom types and associated data loaded from TPSA data file as
a reference to hash with the following hash data format:

    @{$TPSAAtomTypesDataMap{AtomTypes}} - Array of all possible atom
                                          types for all atoms
    @{$TPSAAtomTypesDataMap->{ColLabels}} - Array of column labels
    %{$TPSAAtomTypesDataMap->{DataCol<Num>}} - Hash keys pair:
                                               DataCol<Num>, AtomType

=item B<IsAtomTypesAssignmentSuccessful>

    $Status = $AtomTypes->IsAtomTypesAssignmentSuccessful();

Returns 1 or 0 based on whether atom types assignment was successfully performed.
This method overrides the same method available in the base class AtomTypes.pm used
to derived this class. It checks for successful atom type assignments for nitrogen and
oxygen atoms with an optional check for phosphorous and sulfur atoms.

=item B<StringifyTPSAAtomTypes>

    $String = $TPSAAtomTypes->StringifyTPSAAtomTypes();

Returns a string containing information about I<TPSAAtomTypes> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AtomTypes.pm, AtomicInvariantsAtomTypes.pm, DREIDINGAtomTypes.pm,
EStateAtomTypes.pm, FunctionalClassAtomTypes.pm, MMFF94AtomTypes.pm,
SLogPAtomTypes.pm, SYBYLAtomTypes.pm, UFFAtomTypes.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
