package AtomTypes::SLogPAtomTypes;
#
# File: SLogPAtomTypes.pm
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
use FileUtil ();
use AtomTypes::AtomTypes;
use Molecule;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(AtomTypes::AtomTypes Exporter);
@EXPORT = qw(GetSLogPAtomTypesData GetAllPossibleSLogPAtomTypes GetAllPossibleSLogPNonHydrogenAtomTypes);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, %SLogPAtomTypesDataMap);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifySLogPAtomTypes';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeSLogPAtomTypes();

  $This->_InitializeSLogPAtomTypesProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Initialize the data hash. It'll be loaded on demand later...
  %SLogPAtomTypesDataMap = ();
}

# Initialize object data...
#
sub _InitializeSLogPAtomTypes {
  my($This) = @_;

  # Type of AtomTypes...
  $This->{Type} = 'SLogP';

  # By default, SLogP atom types are also assigned to hydrogens...
  $This->{IgnoreHydrogens} = 0;

  return $This;
}

# Initialize object properties...
#
sub _InitializeSLogPAtomTypesProperties {
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

# Get SLogP atom types and associated data loaded from SLogP data file as
# a reference to hash with the following hash data format:
#
# @{$SLogPAtomTypesDataMap{AtomTypes}} - Array of all possible atom types for all atoms
# @{$SLogPAtomTypesDataMap{NonHydrogenAtomTypes}} - Array of all possible atom types for non-hydrogen atoms
# @{$SLogPAtomTypesDataMap->{ColLabels}} - Array of column labels
# %{$SLogPAtomTypesDataMap->{DataCol<Num>}} - Hash keys pair: <DataCol<Num>, AtomType>
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetSLogPAtomTypesData {

  # Make sure data is loaded...
  _CheckAndLoadSLogPAtomTypesData();

  return \%SLogPAtomTypesDataMap;
}

# Get all possible SLogP atom types corresponding to hydrogen and non-hydrogen
# atoms as an array reference...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAllPossibleSLogPAtomTypes {
  return _GetAllPossibleSLogPAtomTypes();
}

# Get all possible SLogP atom types corresponding to non-hydrogen atoms
# as an array reference...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAllPossibleSLogPNonHydrogenAtomTypes {
  my($NonHydrogensOnly);

  $NonHydrogensOnly = 1;
  return _GetAllPossibleSLogPAtomTypes($NonHydrogensOnly);
}

# Get all possible SLogP atom types as an array reference...
#
sub _GetAllPossibleSLogPAtomTypes {
  my($NonHydrogensOnly) = @_;
  my($SLogPAtomTypesDataRef);

  $NonHydrogensOnly = defined $NonHydrogensOnly ? $NonHydrogensOnly : 0;

  $SLogPAtomTypesDataRef = GetSLogPAtomTypesData();

  return $NonHydrogensOnly ? \@{$SLogPAtomTypesDataRef->{NonHydrogenAtomTypes}}: \@{$SLogPAtomTypesDataRef->{AtomTypes}};
}

# Assign SLogP [ Ref 89 ] atom types to all atoms...
#
# Notes:
#     o 72 SLogP atom type symbols are listed
#     o Number of atom types symbols for:
#         o C: 28
#         o N: 15
#         o O: 13
#         o P: 1
#         o S: 3
#         o H: 5
#         o F, Cl, Br, I: 1 each
#         o Ionic halogen: 1
#         o p-block elements: 1
#         o d-block elements: 1
#
sub AssignAtomTypes {
  my($This) = @_;
  my($Atom, $AtomType);

  ATOM: for $Atom ($This->GetMolecule()->GetAtoms()) {
    if ($This->{IgnoreHydrogens} && $Atom->IsHydrogen()) {
      next ATOM;
    }
    $AtomType = $This->_GetAtomType($Atom);
    $This->SetAtomType($Atom, $AtomType);
  }

  return $This;
}

# Get SLogP atom type for atom...
#
sub _GetAtomType {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOM: {
    if ($Atom->IsCarbon()) {
      $AtomType = $This->_GetAtomTypeForCarbon($Atom);
      last ATOM;
    }

    if ($Atom->IsNitrogen()) {
      $AtomType = $This->_GetAtomTypeForNitrogen($Atom);
      last ATOM;
    }

    if ($Atom->IsOxygen()) {
      $AtomType = $This->_GetAtomTypeForOxygen($Atom);
      last ATOM;
    }

    if ($Atom->IsPhosphorus()) {
      $AtomType = $This->_GetAtomTypeForPhosphorus($Atom);
      last ATOM;
    }

    if ($Atom->IsSulfur()) {
      $AtomType = $This->_GetAtomTypeForSulfur($Atom);
      last ATOM;
    }

    if ($Atom->IsHydrogen()) {
      $AtomType = $This->_GetAtomTypeForHydrogen($Atom);
      last ATOM;
    }

    $AtomType = $This->_GetAtomTypeForOtherAtoms($Atom);
  }

  return $AtomType;
}

# Get SLogP atom type for Carbon atom...
#
# 28 AtomTypeSymbols for element C:
#
# AtomTypeSymbol - Description - SMARTS
# C1 - primary, secondary aliphatic - '[CH4]','[CH3]C','[CH2](C)C'
# C2 - tertiary, quaternary aliphatic - '[CH](C)(C)C','[C](C)(C)(C)C'
# C3 - primary, secondary heteroatom - '[CH3][(N,O,P,S,F,Cl,Br,I)]','[CH2X4](N,O,P,S,F,Cl,Br,I)]'
# C4 - tertiary, quaternary heteroatom - '[CH1X4][(N,O,P,S,F,Cl,Br,I)]','[CH0X4][(N,O,P,S,F,Cl,Br,I)]'
# C5 - C = heteroatom - '[C]=[A#X]'
# C6 - C = C aliphatic - '[CH2]=C','[CH1](=C)A','[CH0](=C)(A)A','[C](=C)=C'
# C7 - acetylene, nitrile - '[CX2]#A'
# C8 - primary aromatic carbon - '[CH3]c'
# C9 - primary aromatic heteroatom - '[CH3][a#X]'
# C10 - secondary aromatic - '[CH2X4]a'
# C11 - tertiary aromatic - '[CHX4]a'
# C12 - quaternary aromatic - '[CH0X4]a'
# C13 - aromatic heteroatom - '[cH0]-[!(C,N,O,S,F,Cl,Br,I)]'
# C14 - aromatic halide - '[c][#9]'
# C15 - aromatic halide - '[c][#17]'
# C16 - aromatic halide - '[c][#35]'
# C17 - aromatic halide - '[c][#53]'
# C18 - aromatic - '[cH]'
# C19 - aromatic bridgehead - '[c](:a)(:a):a'
# C20 - quaternary aromatic - '[c](:a)(:a)-a'
# C21 - quaternary aromatic - '[c](:a)(:a)-C'
# C22 - quaternary aromatic - '[c](:a)(:a)-N'
# C23 - quaternary aromatic - '[c](:a)(:a)-O'
# C24 - quaternary aromatic - '[c](:a)(:a)-S'
# C25 - quaternary aromatic - '[c](:a)(:a)=C','[c](:a)(:a)=N','[c](:a)(:a)=O'
# C26 - C = C aromatic - '[C](=C)(a)A','[C](=C)(c)a','[CH](=C)a','[C]=c'
# C27 - aliphatic heteroatom - '[CX4][!(C,N,O,P,S,F,Cl,Br,I)]'
# CS - carbon supplemental not matching any basic C type - '[#6]'
#
sub _GetAtomTypeForCarbon {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'None';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {

    # Only single bonds...
    if ($NumOfPiBonds == 0) {
      $AtomType = $This->_GetAtomTypeForCarbonWithOnlySigmaBonds($Atom);
      last ATOMTYPE;
    }

    # One double bond...
    if ($NumOfPiBonds == 1) {
      $AtomType = $This->_GetAtomTypeForCarbonWithOnePiBond($Atom);
      last ATOMTYPE;
    }

    # One triple bond or two double bonds...
    if ($NumOfPiBonds == 2) {
      $AtomType = $This->_GetAtomTypeForCarbonWithTwoPiBonds($Atom);
      last ATOMTYPE;
    }

    $AtomType = 'CS';
  }
  return $AtomType;
}

# Get SLogP atom type for Nitrogen atom...
#
# 15 AtomTypeSymbols for element N:
#
# AtomTypeSymbol - Description - SMARTS
# N1 - primary amine - '[NH2+0]A'
# N2 - secondary amine - '[NH+0](A)A'
# N3 - primary aromatic amine - '[NH2+0]a'
# N4 - secondary aromatic amine - '[NH+0](A)a','[NH+0](a)a'
# N5 - imine - '[NH+0]=A','[NH+0]=a'
# N6 - substituted imine - '[N+0](=A)A','[N+0](=A)a','[N+0](=a)A','[N+0](=a)a'
# N7 - tertiary amine - '[N+0](A)(A)A'
# N8 - tertiary aromatic amine - '[N+0](a)(A)A','[N+0](a)(a)A','[N+0](a)(a)a'
# N9 - nitrile - '[N+0]#A'
# N10 - protonated amine - '[NH3+*]','[NH2+*]','[NH+*]'
# N11 - unprotonated aromatic - '[n+0]'
# N12 - protonated aromatic - '[n+*]'
# N13 - quaternary amine - '[NH0+*](A)(A)(A)A','[NH0+*](=A)(A)A','[NH0+*](=A)(A)a','[NH0+*](=[#6])=[#7]'
# N14 - other ionized nitrogen - '[N+*]#A','[N-*]','[N+*](=[N-*])=N'
# NS - nitrogen supplemental not matching any basic N type - '[#7]'
#
sub _GetAtomTypeForNitrogen {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'None';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {

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

    $AtomType = 'NS';
  }
  return $AtomType;
}

# Get SLogP atom type for Oxygen atom...
#
# 13 AtomTypeSymbols for element O:
#
# AtomTypeSymbol - Description - SMARTS
# O1 - aromatic - '[o]'
# O2 - alcohol - '[OH]','[OH2]'
# O3 - aliphatic ether - '[O](C)C','[O](C)[A#X]','[O]([A#X])[A#X]'
# O4 - aromatic ether - '[O](A)a','[O](a)a'
# O5 - oxide - '[O]=[#8]','[O]=[#7]','[OX1-*][#7]'
# O6 - oxide - '[OX1-*][#16]'
# O7 - oxide - '[OX1-*][!(N,S)]'
# O8 - aromatic carbonyl - '[O]=c'
# O9 - carbonyl aliphatic - '[O]=[CH]C','[O]=C(C)C','[O]=C(C)[A#X]','[O]=[CH]N','[O]=[CH]O','[O]=[CH2]','[O]=[CX2]=O'
# O10 - carbonyl aromatic - '[O]=[CH]c','[O]=C(C)c','[O]=C(c)c','[O]=C(c)[a#X]','[O]=C(c)[A#X]','[O]=C(C)[a#X]'
# O11 - carbonyl heteroatom - '[O]=C([A#X])[A#X]','[O]=C([A#X])[a#X]','[O]=C([a#X])[a#X]'
# O12 - acid - '[O-1]C(=O)'
# OS - oxygen supplemental not matching any basic O type - '[#8]'
#
sub _GetAtomTypeForOxygen {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'None';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {

    # Only single bonds...
    if ($NumOfPiBonds == 0) {
      $AtomType = $This->_GetAtomTypeForOxygenWithOnlySigmaBonds($Atom);
      last ATOMTYPE;
    }

    # One double bond...
    if ($NumOfPiBonds == 1) {
      $AtomType = $This->_GetAtomTypeForOxygenWithOnePiBond($Atom);
      last ATOMTYPE;
    }

    # OS - oxygen supplemental not matching any basic O type - '[#8]'
    $AtomType = 'OS';
  }
  return $AtomType;
}

# Get SLogP atom type for Phosphorus atom...
#
# 1 AtomTypeSymbols for element P:
#
# AtomTypeSymbol - Description - SMARTS
# P - phosphorous - '[#15]'
#
sub _GetAtomTypeForPhosphorus {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'P';

  return $AtomType;
}

# Get SLogP atom type for Sulfur atom...
#
# 3 AtomTypeSymbols for element S:
#
# AtomTypeSymbol - Description - SMARTS
# S1 - aliphatic - '[S-0]'
# S2 - ionic sulfur - '[S-*]','[S+*]'
# S3 - aromatic - '[s]'
#
sub _GetAtomTypeForSulfur {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # S1 - aliphatic - '[S-0]'
    if ($This->_IsS1Sulfur($Atom)) {
      $AtomType = 'S1';
      last ATOMTYPE;
    }

    # S2 - ionic sulfur - '[S-*]','[S+*]'
    if ($This->_IsS2Sulfur($Atom)) {
      $AtomType = 'S2';
      last ATOMTYPE;
    }

    # S3 - aromatic - '[s]'
    if ($This->_IsS3Sulfur($Atom)) {
      $AtomType = 'S3';
      last ATOMTYPE;
    }

    # S1 - aliphatic - '[S-0]'
    $AtomType = 'S1';
  }

  return $AtomType;
}

# Get SLogP atom type for Hydrogen atom...
#
# 5 AtomTypeSymbols for element H:
#
# AtomTypeSymbol - Description - SMARTS
# H1 - hydrocarbon - '[#1][#6]','[#1][#1]'
# H2 - alcohol - '[#1]O[CX4]','[#1]Oc','[#1]O[!(C,N,O,S)]','[#1][!C,N,O)]'
# H3 - amine - '[#1][#7]','[#1]O[#7]'
# H4 - acid - '[#1]OC=[#6]','[#1]OC=[#7]','[#1]OC=O','[#1]OC=S','[#1]OO','[#1]OS'
# HS - hydrogen supplemental not matching any basic H type - '[#1]'
#
sub _GetAtomTypeForHydrogen {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # H1 - hydrocarbon - '[#1][#6]','[#1][#1]'
    if ($This->_IsH1Hydrogen($Atom)) {
      $AtomType = 'H1';
      last ATOMTYPE;
    }

    # H2 - alcohol - '[#1]O[CX4]','[#1]Oc','[#1]O[!(C,N,O,S)]','[#1][!C,N,O)]'
    if ($This->_IsH2Hydrogen($Atom)) {
      $AtomType = 'H2';
      last ATOMTYPE;
    }

    # H3 - amine - '[#1][#7]','[#1]O[#7]'
    if ($This->_IsH3Hydrogen($Atom)) {
      $AtomType = 'H3';
      last ATOMTYPE;
    }

    # H4 - acid - '[#1]OC=[#6]','[#1]OC=[#7]','[#1]OC=O','[#1]OC=S','[#1]OO','[#1]OS'
    if ($This->_IsH4Hydrogen($Atom)) {
      $AtomType = 'H4';
      last ATOMTYPE;
    }

    $AtomType = 'HS';
  }
  return $AtomType;
}

# Get SLogP atom type for atoms other than Carbon, Nitrogen, Oxygen, Phosporus,
# Sulfur and Hydrogen...
#
# AtomTypeSymbol - Description - SMARTS
# F - fluorine - '[#9-0]'
# Cl - chlorine - '[#17-0]'
# Br - bromine - '[#35-0]'
# I - iodine - '[#53-0]'
# Hal - ionic halogens - '[#9-*]','[#17-*]','[#35-*]',[#53-*]','[#53+*]'
# Hal - all remaining s-block elements
# Me1 - all remaining p-block elements
# Me2 - all remaining d-block elements
#
sub _GetAtomTypeForOtherAtoms {
  my($This, $Atom) = @_;
  my($AtomType, $AtomSymbol);

  $AtomType = 'None';
  $AtomSymbol = $Atom->GetAtomSymbol();

  ATOMTYPE: {

    # Halogens...
    if ($AtomSymbol =~ /^(F|Cl|Br|I)$/i) {
      $AtomType = $Atom->GetFormalCharge() ? 'Hal' : $AtomSymbol;
      last ATOMTYPE;
    }

    # Me1 - all remaining p-block elements
    if ($This->_IsPBlockElement($Atom)) {
      $AtomType = 'Me1';
      last ATOMTYPE;
    }

    # Me2 - all remaining d-block elements
    if ($This->_IsDBlockElement($Atom)) {
      $AtomType = 'Me2';
      last ATOMTYPE;
    }

    # Hal - all remaining s-block elements
    if ($This->_IsSBlockElement($Atom)) {
      $AtomType = 'Hal';
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForOtherAtoms: SLogP atom type for $AtomSymbol cann't be assigned...";
  }
  return $AtomType;
}

# Get SLogP atom type for Carbon with only sigma bonds...
#
sub _GetAtomTypeForCarbonWithOnlySigmaBonds {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # C1 - primary, secondary aliphatic - '[CH4]','[CH3]C','[CH2](C)C'
    if ($This->_IsC1Carbon($Atom)) {
      $AtomType = 'C1';
      last ATOMTYPE;
    }

    # C2 - tertiary, quaternary aliphatic - '[CH](C)(C)C','[C](C)(C)(C)C'
    if ($This->_IsC2Carbon($Atom)) {
      $AtomType = 'C2';
      last ATOMTYPE;
    }

    # C3 - primary, secondary heteroatom - '[CH3][(N,O,P,S,F,Cl,Br,I)]','[CH2X4](N,O,P,S,F,Cl,Br,I)]'
    if ($This->_IsC3Carbon($Atom)) {
      $AtomType = 'C3';
      last ATOMTYPE;
    }

    # C4 - tertiary, quaternary heteroatom - '[CH1X4][(N,O,P,S,F,Cl,Br,I)]','[CH0X4][(N,O,P,S,F,Cl,Br,I)]'
    if ($This->_IsC4Carbon($Atom)) {
      $AtomType = 'C4';
      last ATOMTYPE;
    }

    # C8 - primary aromatic carbon - '[CH3]c'
    if ($This->_IsC8Carbon($Atom)) {
      $AtomType = 'C8';
      last ATOMTYPE;
    }

    # C9 - primary aromatic heteroatom - '[CH3][a#X]'
    if ($This->_IsC9Carbon($Atom)) {
      $AtomType = 'C9';
      last ATOMTYPE;
    }

    # C10 - secondary aromatic - '[CH2X4]a'
    if ($This->_IsC10Carbon($Atom)) {
      $AtomType = 'C10';
      last ATOMTYPE;
    }

    # C11 - tertiary aromatic - '[CHX4]a'
    if ($This->_IsC11Carbon($Atom)) {
      $AtomType = 'C11';
      last ATOMTYPE;
    }

    # C12 - quaternary aromatic - '[CH0X4]a'
    if ($This->_IsC12Carbon($Atom)) {
      $AtomType = 'C12';
      last ATOMTYPE;
    }

    # C27 - aliphatic heteroatom - '[CX4][!(C,N,O,P,S,F,Cl,Br,I)]'
    if ($This->_IsC27Carbon($Atom)) {
      $AtomType = 'C27';
      last ATOMTYPE;
    }

    # CS - carbon supplemental not matching any basic C type - '[#6]'
    $AtomType = 'CS';
  }

  return $AtomType;
}

# Get SLogP atom type for Carbon with one pi bond...
#
sub _GetAtomTypeForCarbonWithOnePiBond {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # C5 - C = heteroatom - '[C]=[A#X]'
    if ($This->_IsC5Carbon($Atom)) {
      $AtomType = 'C5';
      last ATOMTYPE;
    }

    # C6 - C = C aliphatic - '[CH2]=C','[CH1](=C)A','[CH0](=C)(A)A','[C](=C)=C'
    if ($This->_IsC6Carbon($Atom)) {
      $AtomType = 'C6';
      last ATOMTYPE;
    }

    # C13 - aromatic heteroatom - '[cH0]-[!(C,N,O,S,F,Cl,Br,I)]'
    if ($This->_IsC13Carbon($Atom)) {
      $AtomType = 'C13';
      last ATOMTYPE;
    }

    # C14 - aromatic halide - '[c][#9]'
    if ($This->_IsC14Carbon($Atom)) {
      $AtomType = 'C14';
      last ATOMTYPE;
    }

    # C15 - aromatic halide - '[c][#17]'
    if ($This->_IsC15Carbon($Atom)) {
      $AtomType = 'C15';
      last ATOMTYPE;
    }

    # C16 - aromatic halide - '[c][#35]'
    if ($This->_IsC16Carbon($Atom)) {
      $AtomType = 'C16';
      last ATOMTYPE;
    }

    # C17 - aromatic halide - '[c][#53]'
    if ($This->_IsC17Carbon($Atom)) {
      $AtomType = 'C17';
      last ATOMTYPE;
    }

    # C18 - aromatic - '[cH]'
    if ($This->_IsC18Carbon($Atom)) {
      $AtomType = 'C18';
      last ATOMTYPE;
    }

    # C19 - aromatic bridgehead - '[c](:a)(:a):a'
    if ($This->_IsC19Carbon($Atom)) {
      $AtomType = 'C19';
      last ATOMTYPE;
    }

    # C20 - quaternary aromatic - '[c](:a)(:a)-a'
    if ($This->_IsC20Carbon($Atom)) {
      $AtomType = 'C20';
      last ATOMTYPE;
    }

    # C21 - quaternary aromatic - '[c](:a)(:a)-C'
    if ($This->_IsC21Carbon($Atom)) {
      $AtomType = 'C21';
      last ATOMTYPE;
    }

    # C22 - quaternary aromatic - '[c](:a)(:a)-N'
    if ($This->_IsC22Carbon($Atom)) {
      $AtomType = 'C22';
      last ATOMTYPE;
    }

    # C23 - quaternary aromatic - '[c](:a)(:a)-O'
    if ($This->_IsC23Carbon($Atom)) {
      $AtomType = 'C23';
      last ATOMTYPE;
    }

    # C24 - quaternary aromatic - '[c](:a)(:a)-S'
    if ($This->_IsC24Carbon($Atom)) {
      $AtomType = 'C24';
      last ATOMTYPE;
    }

    # C26 - C = C aromatic - '[C](=C)(a)A','[C](=C)(c)a','[CH](=C)a','[C]=c'
    if ($This->_IsC26Carbon($Atom)) {
      $AtomType = 'C26';
      last ATOMTYPE;
    }

    # CS - carbon supplemental not matching any basic C type - '[#6]'
    $AtomType = 'CS';
  }

  return $AtomType;
}

# Get SLogP atom type for Carbon with two pi bonds...
#
sub _GetAtomTypeForCarbonWithTwoPiBonds {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # C6 - C = C aliphatic - '[CH2]=C','[CH1](=C)A','[CH0](=C)(A)A','[C](=C)=C'
    if ($This->_IsC6Carbon($Atom)) {
      $AtomType = 'C6';
      last ATOMTYPE;
    }

    # C7 - acetylene, nitrile - '[CX2]#A'
    if ($This->_IsC7Carbon($Atom)) {
      $AtomType = 'C7';
      last ATOMTYPE;
    }

    # CS - carbon supplemental not matching any basic C type - '[#6]'
    $AtomType = 'CS';
  }

  return $AtomType;
}

# C1 - primary, secondary aliphatic - '[CH4]','[CH3]C','[CH2](C)C'
#
sub _IsC1Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T4.!Ar', ['C.!Ar,H', 'C.!Ar,H', 'H', 'H'], ['-', '-', '-', '-']) ? 1 : 0;
}

# C2 - tertiary, quaternary aliphatic - '[CH](C)(C)C','[C](C)(C)(C)C'
#
sub _IsC2Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T4.!Ar', ['C.!Ar', 'C.!Ar', 'C.!Ar', 'C.!Ar,H'], ['-', '-', '-', '-']) ? 1 : 0;
}

# C3 - primary, secondary heteroatom - '[CH3][(N,O,P,S,F,Cl,Br,I)]','[CH2X4](N,O,P,S,F,Cl,Br,I)]'
#
sub _IsC3Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T4.!Ar', ['N.!Ar,O.!Ar,P.!Ar,S.!Ar,F,Cl,Br,I', 'N.!Ar,O.!Ar,P.!Ar,S.!Ar,F,Cl,Br,I,H', 'H', 'H'], ['-', '-', '-', '-']) ? 1 : 0;
}

# C4 - tertiary, quaternary heteroatom - '[CH1X4][(N,O,P,S,F,Cl,Br,I)]','[CH0X4][(N,O,P,S,F,Cl,Br,I)]'
#
sub _IsC4Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T4.!Ar', ['N.!Ar,O.!Ar,P.!Ar,.!ArS,F,Cl,Br,I', 'N.!Ar,O.!Ar,P.!Ar,S.!Ar,F,Cl,Br,I', 'N.!Ar,O.!Ar,P.!Ar,S.!Ar,F,Cl,Br,I', 'N.!Ar,O.!Ar,P.!Ar,S.!Ar,F,Cl,Br,I,H'], ['-', '-', '-', '-']) ? 1 : 0;
}

# C5 - C = heteroatom - '[C]=[A#X]'
#
sub _IsC5Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.!Ar', ['!C.!Ar'], ['=']) ? 1 : 0;
}

# C6 - C = C aliphatic - '[CH2]=C','[CH1](=C)A','[CH0](=C)(A)A','[C](=C)=C'
#
sub _IsC6Carbon {
  my($This, $Atom) = @_;

  if ($Atom->DoesAtomNeighborhoodMatch('C.DB1.!Ar', ['C.!Ar', 'H,!H.!Ar', 'H,!H.!Ar'], ['=', '-', '-']) ||
     $Atom->DoesAtomNeighborhoodMatch('C.H0.DB2.!Ar', ['C.!Ar', 'C.!Ar'], ['=', '='])) {
    return 1;
  }
  return 0;
}

# C7 - acetylene, nitrile - '[CX2]#A'
#
sub _IsC7Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T2.TB1', ['*', '*'], ['#', '-']) ? 1 : 0;
}

# C8 - primary aromatic carbon - '[CH3]c'
#
sub _IsC8Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T4.!Ar', ['C.Ar', 'H', 'H', 'H'], ['-', '-', '-', '-']) ? 1 : 0;
}

# C9 - primary aromatic heteroatom - '[CH3][a#X]'
#
sub _IsC9Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T4.!Ar', ['!C.Ar', 'H', 'H', 'H'], ['-', '-', '-', '-']) ? 1 : 0;
}

# C10 - secondary aromatic - '[CH2X4]a'
#
sub _IsC10Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T4.!Ar', ['!H.Ar', '!H', 'H', 'H'], ['-', '-', '-', '-']) ? 1 : 0;
}

# C11 - tertiary aromatic - '[CHX4]a'
#
sub _IsC11Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T4.!Ar', ['!H.Ar', '!H', '!H', 'H'], ['-', '-', '-', '-']) ? 1 : 0;
}

# C12 - quaternary aromatic - '[CH0X4]a'
#
sub _IsC12Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T4.!Ar', ['!H.Ar', '!H', '!H', '!H'], ['-', '-', '-', '-']) ? 1 : 0;
}

# C13 - aromatic heteroatom - '[cH0]-[!(C,N,O,S,F,Cl,Br,I)]' or matching [ Table 1 annotations, Ref 89 ]
# '[cH0]-[B,Si,P,As,Se,Sn,Hg]'
#
#
sub _IsC13Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.Ar.H0', ['B,Si,P,As,Se,Sn,Hg'], ['-']) ? 1 : 0;
}

# C14 - aromatic halide - '[c][#9]'
#
sub _IsC14Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.Ar', ['F'], ['-']) ? 1 : 0;
}

# C15 - aromatic halide - '[c][#17]'
#
sub _IsC15Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.Ar', ['Cl'], ['-']) ? 1 : 0;
}

# C16 - aromatic halide - '[c][#35]'
#
sub _IsC16Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.Ar', ['Br'], ['-']) ? 1 : 0;
}

# C17 - aromatic halide - '[c][#53]'
#
sub _IsC17Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.Ar', ['I'], ['-']) ? 1 : 0;
}

# C18 - aromatic - '[cH]'
#
sub _IsC18Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.Ar', ['H'], ['-']) ? 1 : 0;
}

# C19 - aromatic bridgehead - '[c](:a)(:a):a'
#
sub _IsC19Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.Ar', ['!H.Ar', '!H.Ar', '!H.Ar'], [':', ':', ':']) ? 1 : 0;
}

# C20 - quaternary aromatic - '[c](:a)(:a)-a'
#
sub _IsC20Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.Ar', ['!H.Ar', '!H.Ar', '!H.Ar'], [':', ':', '-.!:']) ? 1 : 0;
}

# C21 - quaternary aromatic - '[c](:a)(:a)-C'
#
sub _IsC21Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.Ar', ['!H.Ar', '!H.Ar', 'C.!Ar'], [':', ':', '-']) ? 1 : 0;
}

# C22 - quaternary aromatic - '[c](:a)(:a)-N'
#
sub _IsC22Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.Ar', ['!H.Ar', '!H.Ar', 'N.!Ar'], [':', ':', '-']) ? 1 : 0;
}

# C23 - quaternary aromatic - '[c](:a)(:a)-O'
#
sub _IsC23Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.Ar', ['!H.Ar', '!H.Ar', 'O.!Ar'], [':', ':', '-']) ? 1 : 0;
}

# C24 - quaternary aromatic - '[c](:a)(:a)-S'
#
sub _IsC24Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.Ar', ['!H.Ar', '!H.Ar', 'S.!Ar'], [':', ':', '-']) ? 1 : 0;
}

# C25 - quaternary aromatic - '[c](:a)(:a)=C','[c](:a)(:a)=N','[c](:a)(:a)=O'
#
sub _IsC25Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.Ar', ['!H.Ar', '!H.Ar', 'C.!Ar,N.!Ar,O.!Ar'], [':', ':', '=']) ? 1 : 0;
}

# C26 - C = C aromatic - '[C](=C)(a)A','[C](=C)(c)a','[CH](=C)a','[C]=c'
#
sub _IsC26Carbon {
  my($This, $Atom) = @_;

  if ($Atom->DoesAtomNeighborhoodMatch('C.!Ar.T3', ['C.!Ar', '!H.Ar', '!H.!Ar'], ['=', '-', '-']) ||
     $Atom->DoesAtomNeighborhoodMatch('C.!Ar.T3', ['C.!Ar', 'C.Ar', '!H.Ar'], ['=', '-', '-']) ||
     $Atom->DoesAtomNeighborhoodMatch('C.!Ar.T3.H1', ['C.!Ar', '!H.Ar', 'H'], ['=', '-', '-']) ||
     $Atom->DoesAtomNeighborhoodMatch('C.!Ar', ['C.Ar'], ['='])) {
    return 1;
  }
  return 0;
}

# C27 - aliphatic heteroatom - '[CX4][!(C,N,O,P,S,F,Cl,Br,I)]'
#
# Notes:
#   . X4 implies four neighbors including Hydrogen
#   . For C27 match, at least one of the neighbors must be a hetro atom other
#     than C,N,O,P,S,F,Cl,Br,I. In other words, it's primary, secondary, tertiary,
#     or queaternary aliphatic heteroatom not in (C,N,O,P,S,F,Cl,Br,I) group defined
#     using C3 and C4.
#
sub _IsC27Carbon {
  my($This, $Atom) = @_;
  my($AtomNeighbor, $AtomNeighborSymbol);

  if (!$Atom->DoesAtomNeighborhoodMatch('C.T4.!Ar')) {
    return 0;
  }

  ATOMNEIGHBOR: for $AtomNeighbor ($Atom->GetNonHydrogenAtomNeighbors()) {
    $AtomNeighborSymbol = $AtomNeighbor->GetAtomSymbol();
    # Is it a heteroatom?
    if ($AtomNeighborSymbol =~ /^(C|N|O|P|S|F|Cl|Br|I|H)$/) {
      next ATOMNEIGHBOR;
    }
    # Is it aromatic?
    if ($AtomNeighbor->IsAromatic()) {
      next ATOMNEIGHBOR;
    }
    return 1;
  }
  return 0;
}

# CS - carbon supplemental not matching any basic C type - '[#6]'
#
sub _IsCSCarbon {
  my($This, $Atom) = @_;

  return $Atom->IsCarbon() ? 1 : 0;
}

# Get SLogP atom type for Nitrogen with only sigma bonds...
#
sub _GetAtomTypeForNitrogenWithOnlySigmaBonds {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # N1 - primary amine - '[NH2+0]A'
    if ($This->_IsN1Nitrogen($Atom)) {
      $AtomType = 'N1';
      last ATOMTYPE;
    }

    # N2 - secondary amine - '[NH+0](A)A'
    if ($This->_IsN2Nitrogen($Atom)) {
      $AtomType = 'N2';
      last ATOMTYPE;
    }

    # N3 - primary aromatic amine - '[NH2+0]a'
    if ($This->_IsN3Nitrogen($Atom)) {
      $AtomType = 'N3';
      last ATOMTYPE;
    }

    # N4 - secondary aromatic amine - '[NH+0](A)a','[NH+0](a)a'
    if ($This->_IsN4Nitrogen($Atom)) {
      $AtomType = 'N4';
      last ATOMTYPE;
    }

    # N7 - tertiary amine - '[N+0](A)(A)A'
    if ($This->_IsN7Nitrogen($Atom)) {
      $AtomType = 'N7';
      last ATOMTYPE;
    }

    # N8 - tertiary aromatic amine - '[N+0](a)(A)A','[N+0](a)(a)A','[N+0](a)(a)a'
    if ($This->_IsN8Nitrogen($Atom)) {
      $AtomType = 'N8';
      last ATOMTYPE;
    }

    # N10 - protonated amine - '[NH3+*]','[NH2+*]','[NH+*]'
    if ($This->_IsN10Nitrogen($Atom)) {
      $AtomType = 'N10';
      last ATOMTYPE;
    }

    # N11 - unprotonated aromatic - '[n+0]'
    if ($This->_IsN11Nitrogen($Atom)) {
      $AtomType = 'N11';
      last ATOMTYPE;
    }

    # N12 - protonated aromatic - '[n+*]'
    if ($This->_IsN12Nitrogen($Atom)) {
      $AtomType = 'N12';
      last ATOMTYPE;
    }

    # N13 - quaternary amine - '[NH0+*](A)(A)(A)A','[NH0+*](=A)(A)A','[NH0+*](=A)(A)a','[NH0+*](=[#6])=[#7]'
    if ($This->_IsN13Nitrogen($Atom)) {
      $AtomType = 'N13';
      last ATOMTYPE;
    }

    # N14 - other ionized nitrogen - '[N+*]#A','[N-*]','[N+*](=[N-*])=N'
    if ($This->_IsN14Nitrogen($Atom)) {
      $AtomType = 'N14';
      last ATOMTYPE;
    }

    $AtomType = 'NS';
  }

  return $AtomType;
}

# Get SLogP atom type for Nitrogen with one pi bond...
#
sub _GetAtomTypeForNitrogenWithOnePiBond {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # N5 - imine - '[NH+0]=A','[NH+0]=a'
    if ($This->_IsN5Nitrogen($Atom)) {
      $AtomType = 'N5';
      last ATOMTYPE;
    }

    # N6 - substituted imine - '[N+0](=A)A','[N+0](=A)a','[N+0](=a)A','[N+0](=a)a'
    if ($This->_IsN6Nitrogen($Atom)) {
      $AtomType = 'N6';
      last ATOMTYPE;
    }

    # N11 - unprotonated aromatic - '[n+0]'
    if ($This->_IsN11Nitrogen($Atom)) {
      $AtomType = 'N11';
      last ATOMTYPE;
    }

    # N12 - protonated aromatic - '[n+*]'
    if ($This->_IsN12Nitrogen($Atom)) {
      $AtomType = 'N12';
      last ATOMTYPE;
    }

    # N13 - quaternary amine - '[NH0+*](A)(A)(A)A','[NH0+*](=A)(A)A','[NH0+*](=A)(A)a','[NH0+*](=[#6])=[#7]'
    if ($This->_IsN13Nitrogen($Atom)) {
      $AtomType = 'N13';
      last ATOMTYPE;
    }

    # N14 - other ionized nitrogen - '[N+*]#A','[N-*]','[N+*](=[N-*])=N'
    if ($This->_IsN14Nitrogen($Atom)) {
      $AtomType = 'N14';
      last ATOMTYPE;
    }

    $AtomType = 'NS';
  }

  return $AtomType;
}

# Get SLogP atom type for Nitrogen with two pi bonds...
#
sub _GetAtomTypeForNitrogenWithTwoPiBonds {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # N9 - nitrile - '[N+0]#A'
    if ($This->_IsN9Nitrogen($Atom)) {
      $AtomType = 'N9';
      last ATOMTYPE;
    }

    # N13 - quaternary amine - '[NH0+*](A)(A)(A)A','[NH0+*](=A)(A)A','[NH0+*](=A)(A)a','[NH0+*](=[#6])=[#7]'
    if ($This->_IsN13Nitrogen($Atom)) {
      $AtomType = 'N13';
      last ATOMTYPE;
    }

    # N14 - other ionized nitrogen - '[N+*]#A','[N-*]','[N+*](=[N-*])=N'
    if ($This->_IsN14Nitrogen($Atom)) {
      $AtomType = 'N14';
      last ATOMTYPE;
    }

    $AtomType = 'NS';
  }

  return $AtomType;
}

# N1 - primary amine - '[NH2+0]A'
#
sub _IsN1Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.!Ar.T3.FC0', ['!H.!Ar', 'H', 'H'], ['-', '-', '-']) ? 1 : 0;
}

# N2 - secondary amine - '[NH+0](A)A'
#
sub _IsN2Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.!Ar.T3.FC0', ['!H.!Ar', '!H.!Ar', 'H'], ['-', '-', '-']) ? 1 : 0;
}

# N3 - primary aromatic amine - '[NH2+0]a'
#
sub _IsN3Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.!Ar.T3.FC0', ['!H.Ar', 'H', 'H'], ['-', '-', '-']) ? 1 : 0;
}

# N4 - secondary aromatic amine - '[NH+0](A)a','[NH+0](a)a'
#
sub _IsN4Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.!Ar.T3.FC0', ['!H.Ar', '!H', 'H'], ['-', '-', '-']) ? 1 : 0;
}

# N5 - imine - '[NH+0]=A','[NH+0]=a'
#
sub _IsN5Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.!Ar.T2.FC0', ['!H', 'H'], ['=', '-']) ? 1 : 0;
}

# N6 - substituted imine - '[N+0](=A)A','[N+0](=A)a','[N+0](=a)A','[N+0](=a)a'
#
sub _IsN6Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.!Ar.T2.FC0', ['!H', '!H'], ['=', '-']) ? 1 : 0;
}

# N7 - tertiary amine - '[N+0](A)(A)A'
#
sub _IsN7Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.!Ar.T3.FC0', ['!H.!Ar', '!H.!Ar', '!H.!Ar'], ['-', '-', '-']) ? 1 : 0;
}

# N8 - tertiary aromatic amine - '[N+0](a)(A)A','[N+0](a)(a)A','[N+0](a)(a)a'
#
sub _IsN8Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.!Ar.T3.FC0', ['!H.Ar', '!H', '!H'], ['-', '-', '-']) ? 1 : 0;
}

# N9 - nitrile - '[N+0]#A'
#
sub _IsN9Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.!Ar.T1.FC0', ['!H.!Ar'], ['#']) ? 1 : 0;
}

# N10 - protonated amine - '[NH3+*]','[NH2+*]','[NH+*]'
#
sub _IsN10Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.!Ar.FC+*.H3,N.!Ar.FC+*.H2,N.!Ar.FC+*.H1') ? 1 : 0;
}

# N11 - unprotonated aromatic - '[n+0]'
#
sub _IsN11Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.Ar.FC0.H0') ? 1 : 0;
}

# N12 - protonated aromatic - '[n+*]'
#
sub _IsN12Nitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.Ar.FC+*.!H0') ? 1 : 0;
}

# N13 - quaternary amine - '[NH0+*](A)(A)(A)A','[NH0+*](=A)(A)A','[NH0+*](=A)(A)a','[NH0+*](=[#6])=[#7]'
#
sub _IsN13Nitrogen {
  my($This, $Atom) = @_;

  if ($Atom->DoesAtomNeighborhoodMatch('N.!Ar.FC+*.T4.H0', ['!H.!Ar', '!H.!Ar', '!H.!Ar', '!H.!Ar'], ['-', '-', '-', '-']) ||
     $Atom->DoesAtomNeighborhoodMatch('N.!Ar.FC+*.T3.H0', ['!H.!Ar', '!H.!Ar', '!H.!Ar'], ['=', '-', '-']) ||
     $Atom->DoesAtomNeighborhoodMatch('N.!Ar.FC+*.T3.H0', ['!H.!Ar', '!H.!Ar', '!H.Ar'], ['=', '-', '-']) ||
     $Atom->DoesAtomNeighborhoodMatch('N.!Ar.FC+*.T2.H0', ['C.!Ar', 'N.!Ar'], ['=', '=']) ) {
    return 1;
  }
  return 0;
}

# N14 - other ionized nitrogen - '[N+*]#A','[N-*]','[N+*](=[N-*])=N'
#
sub _IsN14Nitrogen {
  my($This, $Atom) = @_;

  if ($Atom->DoesAtomNeighborhoodMatch('N.!Ar.FC+*', ['!H.!Ar'], ['#']) ||
     $Atom->DoesAtomNeighborhoodMatch('N.!Ar.FC-*') ||
     $Atom->DoesAtomNeighborhoodMatch('N.!Ar.FC+*', ['N.!Ar.FC-*', 'N.!Ar'], ['=', '=']) ) {
    return 1;
  }
  return 0;
}

# NS - nitrogen supplemental not matching any basic N type - '[#7]'
#
sub _IsNSNitrogen {
  my($This, $Atom) = @_;

  return $Atom->IsNitrogen() ? 1 : 0;
}

# Get SLogP atom type for Oxygen with only sigma bonds...
#
sub _GetAtomTypeForOxygenWithOnlySigmaBonds {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # O1 - aromatic - '[o]'
    if ($This->_IsO1Oxygen($Atom)) {
      $AtomType = 'O1';
      last ATOMTYPE;
    }

    # O2 - alcohol - '[OH]','[OH2]'
    if ($This->_IsO2Oxygen($Atom)) {
      $AtomType = 'O2';
      last ATOMTYPE;
    }

    # O3 - aliphatic ether - '[O](C)C','[O](C)[A#X]','[O]([A#X])[A#X]'
    if ($This->_IsO3Oxygen($Atom)) {
      $AtomType = 'O3';
      last ATOMTYPE;
    }

    # O4 - aromatic ether - '[O](A)a','[O](a)a'
    if ($This->_IsO4Oxygen($Atom)) {
      $AtomType = 'O4';
      last ATOMTYPE;
    }

    # O5 - oxide - '[O]=[#8]','[O]=[#7]','[OX1-*][#7]'
    if ($This->_IsO5Oxygen($Atom)) {
      $AtomType = 'O5';
      last ATOMTYPE;
    }

    # O6 - oxide - '[OX1-*][#16]'
    if ($This->_IsO6Oxygen($Atom)) {
      $AtomType = 'O6';
      last ATOMTYPE;
    }

    # O7 - oxide - '[OX1-*][!(N,S)]'
    if ($This->_IsO7Oxygen($Atom)) {
      $AtomType = 'O7';
      last ATOMTYPE;
    }

    # O12 - acid - '[O-1]C(=O)'
    if ($This->_IsO12Oxygen($Atom)) {
      $AtomType = 'O12';
      last ATOMTYPE;
    }

    $AtomType = 'OS';
  }

  return $AtomType;
}

# Get SLogP atom type for Oxygen with only sigma bonds...
#
sub _GetAtomTypeForOxygenWithOnePiBond {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # O1 - aromatic - '[o]'
    if ($This->_IsO1Oxygen($Atom)) {
      $AtomType = 'O1';
      last ATOMTYPE;
    }

    # O5 - oxide - '[O]=[#8]','[O]=[#7]','[OX1-*][#7]'
    if ($This->_IsO5Oxygen($Atom)) {
      $AtomType = 'O5';
      last ATOMTYPE;
    }

    # O8 - aromatic carbonyl - '[O]=c'
    if ($This->_IsO8Oxygen($Atom)) {
      $AtomType = 'O8';
      last ATOMTYPE;
    }

    # O9 - carbonyl aliphatic - '[O]=[CH]C','[O]=C(C)C','[O]=C(C)[A#X]','[O]=[CH]N','[O]=[CH]O','[O]=[CH2]','[O]=[CX2]=O'
    if ($This->_IsO9Oxygen($Atom)) {
      $AtomType = 'O9';
      last ATOMTYPE;
    }

    # O10 - carbonyl aromatic - '[O]=[CH]c','[O]=C(C)c','[O]=C(c)c','[O]=C(c)[a#X]','[O]=C(c)[A#X]','[O]=C(C)[a#X]'
    if ($This->_IsO10Oxygen($Atom)) {
      $AtomType = 'O10';
      last ATOMTYPE;
    }

    # O11 - carbonyl heteroatom - '[O]=C([A#X])[A#X]','[O]=C([A#X])[a#X]','[O]=C([a#X])[a#X]'
    if ($This->_IsO11Oxygen($Atom)) {
      $AtomType = 'O11';
      last ATOMTYPE;
    }

    $AtomType = 'OS';
  }

  return $AtomType;
}

# O1 - aromatic - '[o]'
#
sub _IsO1Oxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.Ar') ? 1 : 0;
}

# O2 - alcohol - '[OH]','[OH2]'
#
sub _IsO2Oxygen {
  my($This, $Atom) = @_;

  if ($Atom->DoesAtomNeighborhoodMatch('O.TSB2.!Ar', ['C', 'H'], ['-', '-']) ||
     $Atom->DoesAtomNeighborhoodMatch('O.TSB2.!Ar', ['H', 'H'])) {
    return 1;
  }
  return 0;
}

# O3 - aliphatic ether - '[O](C)C','[O](C)[A#X]','[O]([A#X])[A#X]'
#
sub _IsO3Oxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.!Ar.X2', ['!H.!Ar', '!H.!Ar'], ['-', '-']) ? 1 : 0;
}

# O4 - aromatic ether - '[O](A)a','[O](a)a'
#
sub _IsO4Oxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.X2.!Ar', ['!H', '!H.Ar'], ['-', '-']) ? 1 : 0;
}

# O5 - oxide - '[O]=[#8]','[O]=[#7]','[OX1-*][#7]'
#
sub _IsO5Oxygen {
  my($This, $Atom) = @_;

  if ($Atom->DoesAtomNeighborhoodMatch('O.DB1.FC0', ['N,O'], ['=']) ||
     $Atom->DoesAtomNeighborhoodMatch('O.T1.FC-*', ['N'], ['-']))  {
    return 1;
  }
  return 0;
}

# O6 - oxide - '[OX1-*][#16]'
#
sub _IsO6Oxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.T1.FC-*', ['S'], ['-']) ? 1 : 0;
}

# O7 - oxide - '[OX1-*][!(N,S)]'  or matching [ Table 1 annotations, Ref 89 ]
#  '[OX1-*](P,As,Tc,I)
#
sub _IsO7Oxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.T1.FC-*', ['P,As,Tc,I'], ['-']) ? 1 : 0;
}

# O8 - aromatic carbonyl - '[O]=c'
#
sub _IsO8Oxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.DB1.!Ar', ['C.Ar'], ['=']) ? 1 : 0;
}

# O9 - carbonyl aliphatic - '[O]=[CH]C','[O]=C(C)C','[O]=C(C)[A#X]','[O]=[CH]N','[O]=[CH]O','[O]=[CH2]','[O]=[CX2]=O'
#
sub _IsO9Oxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  # Is it a doubly bonded non-aromatic Oxygen?
  if (!$Atom->DoesAtomNeighborhoodMatch('O.DB1.!Ar.FC0')) {
    return 0;
  }

  # Is it attached to appopriate Carbon?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('C.DB1.T3', ['O.!Ar', 'C.!Ar', 'H,C.!Ar'], ['=', '-', '-'])  ||
       $AtomNeighbor->DoesAtomNeighborhoodMatch('C.DB1.T3.H0', ['O.!Ar', '!C.!Ar', 'C.!Ar'], ['=', '-', '-'])  ||
       $AtomNeighbor->DoesAtomNeighborhoodMatch('C.DB1.T3', ['O.!Ar', 'N.!Ar,O.!Ar', 'H'], ['=', '-', '-']) ||
       $AtomNeighbor->DoesAtomNeighborhoodMatch('C.DB1.T3', ['O.!Ar', 'H', 'H'], ['=', '-', '-']) ||
       $AtomNeighbor->DoesAtomNeighborhoodMatch('C.DB2.T2', ['O.!Ar', 'O.!Ar'], ['=', '='])) {
      return 1;
    }
  }
  return 0;
}

# O10 - carbonyl aromatic - '[O]=[CH]c','[O]=C(C)c','[O]=C(c)c','[O]=C(c)[a#X]','[O]=C(c)[A#X]','[O]=C(C)[a#X]'
#
sub _IsO10Oxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  # Is it a doubly bonded non-aromatic Oxygen?
  if (!$Atom->DoesAtomNeighborhoodMatch('O.DB1.!Ar.FC0')) {
    return 0;
  }

  # Is it attached to appopriate Carbon?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('C.DB1.T3', ['O.!Ar', 'C.Ar', 'H'], ['=', '-', '-'])  ||
       $AtomNeighbor->DoesAtomNeighborhoodMatch('C.DB1.T3', ['O.!Ar', 'C.Ar', 'C.!Ar'], ['=', '-', '-'])  ||
       $AtomNeighbor->DoesAtomNeighborhoodMatch('C.DB1.T3', ['O.!Ar', 'C.Ar', 'C.Ar'], ['=', '-', '-'])  ||
       $AtomNeighbor->DoesAtomNeighborhoodMatch('C.DB1.T3.H0', ['O.!Ar', '!C', 'C.Ar'], ['=', '-', '-'])  ||
       $AtomNeighbor->DoesAtomNeighborhoodMatch('C.DB1.T3.H0', ['O.!Ar', '!C.Ar', 'C.!Ar'], ['=', '-', '-'])) {
      return 1;
    }
  }
  return 0;
}

# O11 - carbonyl heteroatom - '[O]=C([A#X])[A#X]','[O]=C([A#X])[a#X]','[O]=C([a#X])[a#X]'
#
sub _IsO11Oxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  # Is it a doubly bonded non-aromatic Oxygen?
  if (!$Atom->DoesAtomNeighborhoodMatch('O.DB1.!Ar.FC0')) {
    return 0;
  }

  # Is it attached to appopriate Carbon?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('C.DB1.T3.H0', ['O.!Ar', '!C', '!C'], ['=', '-', '-'])) {
      return 1;
    }
  }
  return 0;
}

# O12 - acid - '[O-1]C(=O)'
#
sub _IsO12Oxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  # Is it a acid Oxygen?
  if (!$Atom->DoesAtomNeighborhoodMatch('O.DB0.!Ar.FC-1')) {
    return 0;
  }

  # Is it attached to appopriate Carbon?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('C.DB1.T3', ['O.FC0', 'O.FC-1'], ['=', '-'])) {
      return 1;
    }
  }
  return 0;
}

# OS - oxygen supplemental not matching any basic O type - '[#8]'
#
sub _IsOSOxygen {
  my($This, $Atom) = @_;

  return $Atom->IsOxygen() ? 1 : 0;
}

# S1 - aliphatic - '[S-0]'
#
sub _IsS1Sulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.!Ar.FC0') ? 1 : 0;
}

# S2 - ionic sulfur - '[S-*]','[S+*]'
#
sub _IsS2Sulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.!Ar.FC-*,S.!Ar.FC+*') ? 1 : 0;
}

# S3 - aromatic - '[s]'
#
sub _IsS3Sulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.Ar') ? 1 : 0;
}

# Hal - all remaining s-block elements
#
sub _IsSBlockElement {
  my($This, $Atom) = @_;
  my($GroupNumber);

  $GroupNumber = $Atom->GetGroupNumber();

  return ($GroupNumber >= 1 && $GroupNumber <= 2) ? 1 : 0;
}

# Me1 - all remaining p-block elements
#
sub _IsPBlockElement {
  my($This, $Atom) = @_;
  my($GroupNumber);

  $GroupNumber = $Atom->GetGroupNumber();

  return ($GroupNumber >= 13 && $GroupNumber <= 18) ? 1 : 0;
}

# Me2 - all remaining d-block elements
#
sub _IsDBlockElement {
  my($This, $Atom) = @_;
  my($GroupNumber);

  $GroupNumber = $Atom->GetGroupNumber();

  return ($GroupNumber >= 3 && $GroupNumber <= 12) ? 1 : 0;
}

# H1 - hydrocarbon - '[#1][#6]','[#1][#1]'
#
sub _IsH1Hydrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('H', ['C,H']) ? 1 : 0;
}

# H2 - alcohol - '[#1]O[CX4]','[#1]Oc','[#1]O[!(C,N,O,S)]','[#1][!C,N,O)]' or matching [ Table 1 annotations, Ref 89 ]
# '[H]O[CX4]', '[H]Oc', '[H]O[H,B,Si,P,As,Sn]','[H][B,Si,P,S,Sn]'
#
sub _IsH2Hydrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  # Is Hydrogen directly attached to B,Si,P,S,Sn?
  if ($Atom->DoesAtomNeighborhoodMatch('H', ['B,Si,P,S,Sn'], ['-'])) {
    return 1;
  }

  # Is Hydrogen directly attached to Oxygen?
  if (!$Atom->DoesAtomNeighborhoodMatch('H', ['O'])) {
    return 0;
  }

  # Is it attached to appropriate Oxygen?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('O')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('O.T2', ['C.T4', 'H'], ['-', '-'])  ||
       $AtomNeighbor->DoesAtomNeighborhoodMatch('O.T2', ['C.Ar', 'H'], ['-', '-'])  ||
       $AtomNeighbor->DoesAtomNeighborhoodMatch('O.T2', ['B,Si,P,As,Sn', 'H'], ['-', '-']) ||
       $AtomNeighbor->DoesAtomNeighborhoodMatch('O.T2', ['H', 'H'], ['-', '-'])) {
      return 1;
    }
  }
  return 0;
}

# H3 - amine - '[#1][#7]','[#1]O[#7]'
#
sub _IsH3Hydrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  # Is Hydrogen directly attached to Nitrogen?
  if ($Atom->DoesAtomNeighborhoodMatch('H', ['N'])) {
    return 1;
  }

  # Is Hydrogen directly attached to Oxygen?
  if (!$Atom->DoesAtomNeighborhoodMatch('H', ['O'])) {
    return 0;
  }

  # Is it attached to appropriate Oxygen?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('O')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('O.T2', ['N', 'H'], ['-', '-'])) {
      return 1;
    }
  }
  return 0;
}

# H4 - acid - '[#1]OC=[#6]','[#1]OC=[#7]','[#1]OC=O','[#1]OC=S','[#1]OO','[#1]OS'
#
sub _IsH4Hydrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor, $AtomNbrOfNbr);

  # Get non-hydrogen atom neighbor...
  $AtomNeighbor = $Atom->GetNonHydrogenNeighborOfHydrogenAtom();
  if (!$AtomNeighbor) {
    return 0;
  }

  # Is it Oxygen neighbor?
  if (!$AtomNeighbor->IsOxygen()) {
    return 0;
  }

  # '[#1]OO','[#1]OS'
  if ($AtomNeighbor->DoesAtomNeighborhoodMatch('O.TSB2', ['O,S', 'H'], ['-', '-'])) {
    return 1;
  }

  # '[#1]OC=[#6]','[#1]OC=[#7]','[#1]OC=O','[#1]OC=S'
  for $AtomNbrOfNbr ($AtomNeighbor->GetNeighborsUsingAtomSpecification('C.DB1')) {
    if ($AtomNbrOfNbr->DoesAtomNeighborhoodMatch('C.T3.DB1', ['O.H1', 'O,C,N,S', '*'], ['-', '=', '-'])) {
      return 1;
    }
  }
  return 0;
}

# HS - hydrogen supplemental not matching any basic H type - '[#1]'
#
sub _IsHSHydrogen {
  my($This, $Atom) = @_;

  return $Atom->IsHydrogen() ? 1 : 0;
}

# Return a string containg data for SLogPAtomTypes object...
#
sub StringifySLogPAtomTypes {
  my($This) = @_;
  my($AtomTypesString);

  # Type of AtomTypes...
  $AtomTypesString = "AtomTypes: $This->{Type}; IgnoreHydrogens: " . ($This->{IgnoreHydrogens} ? "Yes" : "No");

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

# Is it a SLogPAtomTypes object?
sub _IsSLogPAtomTypes {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Check and load SLogP atom types data...
#
sub _CheckAndLoadSLogPAtomTypesData {

  # Is it already loaded?
  if (exists $SLogPAtomTypesDataMap{AtomTypes}) {
    return;
  }

  _LoadSLogPAtomTypesData();
}

# Load SLogP atom types data from the file assuming first column to be atom type symbol..
#
# Format:
#
# "AtomTypeSymbol","Description","SMARTS","LogPContribution","MRContribution"
# "C1","primary, secondary aliphatic","'[CH4]','[CH3]C','[CH2](C)C'","0.1441","2.503"
# "C2","tertiary, quaternary aliphatic","'[CH](C)(C)C','[C](C)(C)(C)C'","0.0000","2.433"
#
sub _LoadSLogPAtomTypesData {
  my($AtomTypesDataFile, $MayaChemToolsLibDir);

  $MayaChemToolsLibDir = FileUtil::GetMayaChemToolsLibDirName();

  $AtomTypesDataFile =  "$MayaChemToolsLibDir" . "/data/SLogPAtomTypes.csv";
  if (! -e "$AtomTypesDataFile") {
    croak "Error: MayaChemTools package file, $AtomTypesDataFile, is missing: Possible installation problems...";
  }

  %SLogPAtomTypesDataMap = ();
  AtomTypes::AtomTypes::LoadAtomTypesData($AtomTypesDataFile, \%SLogPAtomTypesDataMap);
}

1;

__END__

=head1 NAME

SLogPAtomTypes

=head1 SYNOPSIS

use AtomTypes::SLogPAtomTypes;

use AtomTypes::SLogPAtomTypes qw(:all);

=head1 DESCRIPTION

B<SLogPAtomTypes> class provides the following methods:

new, AssignAtomTypes, GetAllPossibleSLogPAtomTypes,
GetAllPossibleSLogPNonHydrogenAtomTypes, GetSLogPAtomTypesData,
StringifySLogPAtomTypes

The following functions are available:

GetAllPossibleSLogPAtomTypes,
GetAllPossibleSLogPNonHydrogenAtomTypes, GetSLogPAtomTypesData

B<SLogPAtomTypes> is derived from B<AtomTypes> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<SLogPAtomTypes>, B<AtomTypes> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The data file SLogPAomTypes.csv distributed with MayaChemTools release contains
all possible SLogP [ Ref 89 ] atom types.

Examples of SLogP atom types:

    C1, C2, C3, N1, N2, O1, O2 and so on

=head2 METHODS

=over 4

=item B<new>

    $NewSLogPAtomTypes = new AtomTypes::SLogPAtomTypes(%NamesAndValues);

Using specified I<SLogPAtomTypes> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<SLogPAtomTypes>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'SLogP'
    IgnoreHydrogens = 0

Examples:

    $SLogPAtomTypes = new AtomTypes::SLogPAtomTypes(
                              'Molecule' => $Molecule,
                              'IgnoreHydrogens' => 0);

=item B<AssignAtomTypes>

    $SLogPAtomTypes->AssignAtomTypes();

Assigns SLogP atom types to all the atoms in a molecule and returns
I<SLogPAtomTypes>.

=item B<GetAllPossibleSLogPAtomTypes>

    $AllAtomTypesDataRef = $SLogPAtomTypes->
                           GetAllPossibleSLogPAtomTypes();
    $AllAtomTypesDataRef = AtomTypes::SLogPAtomTypes::
                           GetAllPossibleSLogPAtomTypes();

Returns all possible SLogP atom types corresponding to hydrogen and non-hydrogen
atoms as an array reference.

=item B<GetAllPossibleSLogPNonHydrogenAtomTypes>

    $AtomTypesDataRef = $SLogPAtomTypes->
                        GetAllPossibleSLogPNonHydrogenAtomTypes();
    $AtomTypesDataRef = AtomTypes::SLogPAtomTypes::
                        GetAllPossibleSLogPNonHydrogenAtomTypes();

Returns all possible SLogP atom types corresponding to non-hydrogen atoms as
an array reference.

=item B<GetSLogPAtomTypesData>

    $AtomTypesDataMapRef = $SLogPAtomTypes->GetSLogPAtomTypesData();
    $AtomTypesDataMapRef = AtomTypes::SLogPAtomTypes::GetSLogPAtomTypesData();

Returns SLogP atom types and associated data loaded from SLogP data file as
a reference to hash with the following hash data format:

    @{$SLogPAtomTypesDataMap{AtomTypes}} - Array of all possible atom
                              types for all atoms
    @{$SLogPAtomTypesDataMap{NonHydrogenAtomTypes}} - Array of all
                              possible atom types for non-hydrogen atoms
    @{$SLogPAtomTypesDataMap->{ColLabels}} - Array of column labels
    %{$SLogPAtomTypesDataMap->{DataCol<Num>}} - Hash keys pair:
                                                DataCol<Num>, AtomType

=item B<StringifySLogPAtomTypes>

    $String = $SLogPAtomTypes->StringifySLogPAtomTypes();

Returns a string containing information about I<SLogPAtomTypes> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AtomTypes.pm, AtomicInvariantsAtomTypes.pm, DREIDINGAtomTypes.pm,
EStateAtomTypes.pm, FunctionalClassAtomTypes.pm, MMFF94AtomTypes.pm,
SYBYLAtomTypes.pm, TPSAAtomTypes.pm, UFFAtomTypes.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
