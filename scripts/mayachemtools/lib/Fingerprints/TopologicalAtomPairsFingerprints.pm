package Fingerprints::TopologicalAtomPairsFingerprints;
#
# File: TopologicalAtomPairsFingerprints.pm
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
use Fingerprints::Fingerprints;
use TextUtil ();
use Molecule;
use AtomTypes::AtomicInvariantsAtomTypes;
use AtomTypes::DREIDINGAtomTypes;
use AtomTypes::EStateAtomTypes;
use AtomTypes::FunctionalClassAtomTypes;
use AtomTypes::MMFF94AtomTypes;
use AtomTypes::SLogPAtomTypes;
use AtomTypes::SYBYLAtomTypes;
use AtomTypes::TPSAAtomTypes;
use AtomTypes::UFFAtomTypes;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Fingerprints::Fingerprints Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyTopologicalAtomPairsFingerprints';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeTopologicalAtomPairsFingerprints();

  $This->_InitializeTopologicalAtomPairsFingerprintsProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeTopologicalAtomPairsFingerprints {
  my($This) = @_;

  # Type of fingerprint...
  $This->{Type} = 'TopologicalAtomPairs';

  # Type of vector...
  $This->{VectorType} = 'FingerprintsVector';

  # Type of FingerprintsVector...
  $This->{FingerprintsVectorType} = 'NumericalValues';

  # Minimum and maximum bond distance between atom paris...
  $This->{MinDistance} = 1;
  $This->{MaxDistance} = 10;

  # Atom identifier type to use for atom IDs in atom pairs...
  #
  # Currently supported values are: AtomicInvariantsAtomTypes, DREIDINGAtomTypes,
  # EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
  # SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes
  #
  $This->{AtomIdentifierType} = '';

  # Atom types assigned to each heavy atom...
  #
  %{$This->{AssignedAtomTypes}} = ();

  # All atom pairs between minimum and maximum distance...
  #
  @{$This->{AtomPairsIDs}} = ();
  %{$This->{AtomPairsCount}} = ();
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize object properties....
sub _InitializeTopologicalAtomPairsFingerprintsProperties {
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
  if (!exists $NamesAndValues{AtomIdentifierType}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying AtomIdentifierType...";
  }

  $This->_InitializeFingerprintsVector();

  return $This;
}

# Set minimum distance for atom pairs...
#
sub SetMinDistance {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetMinDistance: MinDistance value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{MinDistance} = $Value;

  return $This;
}

# Set maximum distance for atom pairs...
#
sub SetMaxDistance {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetMaxDistance: MaxDistance value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{MaxDistance} = $Value;

  return $This;
}

# Set atom identifier type..
#
sub SetAtomIdentifierType {
  my($This, $IdentifierType) = @_;

  if ($IdentifierType !~ /^(AtomicInvariantsAtomTypes|DREIDINGAtomTypes|EStateAtomTypes|FunctionalClassAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
    croak "Error: ${ClassName}->SetAtomIdentifierType: Specified value, $IdentifierType, for AtomIdentifierType is not vaild. Supported types in current release of MayaChemTools: AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes, and UFFAtomTypes.";
  }

  if ($This->{AtomIdentifierType}) {
    croak "Error: ${ClassName}->SeAtomIdentifierType: Can't change intial atom identifier type:  It's already set...";
  }

  $This->{AtomIdentifierType} = $IdentifierType;

  # Initialize atom identifier type information...
  $This->_InitializeAtomIdentifierTypeInformation();

  return $This;
}

# Generate fingerprints description...
#
sub GetDescription {
  my($This) = @_;

  # Is description explicity set?
  if (exists $This->{Description}) {
    return $This->{Description};
  }

  # Generate fingerprints description...

  return "$This->{Type}:$This->{AtomIdentifierType}:MinDistance$This->{MinDistance}:MaxDistance$This->{MaxDistance}";
}

# Generate topological atom pairs [ Ref 57, Ref 59, Ref 72 ] fingerprints...
#
# Methodology:
#   . Generate a distance matrix.
#   . Assign atom types to all the atoms.
#   . Using distance matrix and atom types, count occurrence of
#     unique atom pairs within specified distance range - It corresponds to the
#     correlation-vector for the atom pairs.
#
# Notes:
#   . Hydrogen atoms are ignored during the fingerprint generation.
#
sub GenerateFingerprints {
  my($This) = @_;

  if ($This->{MinDistance} > $This->{MaxDistance}) {
    croak "Error: ${ClassName}->GenerateTopologicalAtomPairsFingerprints: No fingerpritns generated: MinDistance, $This->{MinDistance}, must be <= MaxDistance, $This->{MaxDistance}...";
  }

  # Cache appropriate molecule data...
  $This->_SetupMoleculeDataCache();

  # Generate distance matrix...
  if (!$This->_SetupDistanceMatrix()) {
    carp "Warning: ${ClassName}->GenerateFingerprints: $This->{AtomIdentifierType} fingerprints generation didn't succeed: Couldn't generate distance matrix...";
    return $This;
  }

  # Assign atom types to all heavy atoms...
  if (!$This->_AssignAtomTypes()) {
    carp "Warning: ${ClassName}->GenerateFingerprints: $This->{AtomIdentifierType} fingerprints generation didn't succeed: Couldn't assign valid $This->{AtomIdentifierType} to all atoms...";
    return $This;
  }

  # Intialize values of toplogical atom pairs...
  $This->_InitializeToplogicalAtomPairs();

  # Count atom pairs...
  $This->_GenerateAndCountAtomPairs();

  # Set final fingerprints...
  $This->_SetFinalFingerprints();

  # Clear cached molecule data...
  $This->_ClearMoleculeDataCache();

  return $This;
}

# Setup distance matrix...
#
sub _SetupDistanceMatrix {
  my($This) = @_;

  $This->{DistanceMatrix} = $This->GetMolecule()->GetDistanceMatrix();

  if (!$This->{DistanceMatrix}) {
    return undef;
  }

  return $This;
}

# Assign appropriate atom types to all heavy atoms...
#
sub _AssignAtomTypes {
  my($This) = @_;
  my($SpecifiedAtomTypes, $Atom, $AtomID, $IgnoreHydrogens);

  %{$This->{AssignedAtomTypes}} = ();
  $IgnoreHydrogens = 1;

  $SpecifiedAtomTypes = undef;

  IDENTIFIERTYPE: {
    if ($This->{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::AtomicInvariantsAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $IgnoreHydrogens, 'AtomicInvariantsToUse' => $This->{AtomicInvariantsToUse});
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^DREIDINGAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::DREIDINGAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $IgnoreHydrogens);
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^EStateAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::EStateAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $IgnoreHydrogens);
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::FunctionalClassAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $IgnoreHydrogens, 'FunctionalClassesToUse' => $This->{FunctionalClassesToUse});
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^MMFF94AtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::MMFF94AtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $IgnoreHydrogens);
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^SLogPAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::SLogPAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $IgnoreHydrogens);
      last IDENTIFIERTYPE;
    }
    if ($This->{AtomIdentifierType} =~ /^SYBYLAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::SYBYLAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $IgnoreHydrogens);
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^TPSAAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::TPSAAtomTypes('Molecule' => $This->{Molecule}, 'IgnorePhosphorus' => 0, 'IgnoreSulfur' => 0);
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^UFFAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::UFFAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $IgnoreHydrogens);
      last IDENTIFIERTYPE;
    }

    croak "Error: ${ClassName}->_AssignAtomTypes: Unknown atom indentifier type $This->{AtomIdentifierType}...";
  }

  # Assign atom types...
  $SpecifiedAtomTypes->AssignAtomTypes();

  # Make sure atom types assignment is successful...
  if (!$SpecifiedAtomTypes->IsAtomTypesAssignmentSuccessful()) {
    return undef;
  }

  # Collect assigned atom types...
  ATOM: for $Atom (@{$This->{Atoms}}) {
    if ($Atom->IsHydrogen()) {
      next ATOM;
    }
    $AtomID = $Atom->GetID();
    $This->{AssignedAtomTypes}{$AtomID} = $SpecifiedAtomTypes->GetAtomType($Atom);
  }

  return $This;
}

# Initialize topological atom pairs between specified distance range...
#
sub _InitializeToplogicalAtomPairs {
  my($This) = @_;
  my($Distance);

  @{$This->{AtomPairsIDs}} = ();
  %{$This->{AtomPairsCount}} = ();

  for $Distance ($This->{MinDistance} .. $This->{MaxDistance}) {
    %{$This->{AtomPairsCount}{$Distance}} = ();
  }

  return $This;
}

# Count atom pairs between mininum and maximum distance at each
# distance using distance matrix and atom types assiged to each heavy
# atom.
#
# Notes:
#   . The row and column indices of distance matrix correspond to atom indices.
#   . Distance value of BigNumber implies the atom is not connected to any other atom.
#   . Due to symmetric nature of distance matrix, only upper or lower triangular matrix
#     needs to be processed during identification and count of atom pairs.
#
sub _GenerateAndCountAtomPairs {
  my($This) = @_;

  my($NumOfRows, $NumOfCols, $RowIndex, $ColIndex, $DistanceMatrix, $Distance, $AtomID1, $AtomID2, $AtomType1, $AtomType2, $SkipIndexCheck, $CountIncrement);

  $DistanceMatrix = $This->{DistanceMatrix};
  ($NumOfRows, $NumOfCols) = $DistanceMatrix->GetSize();
  $SkipIndexCheck = 0;

  ROWINDEX: for $RowIndex (0 .. ($NumOfRows - 1) ) {
    $AtomID1 = $This->{AtomIndexToID}{$RowIndex};
    if ( !(exists($This->{AssignedAtomTypes}{$AtomID1})) ) {
      next ROWINDEX;
    }
    $AtomType1 = $This->{AssignedAtomTypes}{$AtomID1};

    COLINDEX: for $ColIndex (($RowIndex + 1) .. ($NumOfCols - 1) ) {
      $AtomID2 = $This->{AtomIndexToID}{$ColIndex};
      if ( !(exists($This->{AssignedAtomTypes}{$AtomID2})) ) {
	next COLINDEX;
      }
      $Distance = $DistanceMatrix->GetValue($RowIndex, $ColIndex, $SkipIndexCheck);
      if ($Distance < $This->{MinDistance} || $Distance > $This->{MaxDistance}) {
	next COLINDEX;
      }
      $AtomType2 = $This->{AssignedAtomTypes}{$AtomID2};

      if ($AtomType1 le $AtomType2) {
	$This->_SetAtomPairsCount($Distance, $AtomType1, $AtomType2);
      }
      else {
	$This->_SetAtomPairsCount($Distance, $AtomType2, $AtomType1);
      }
    }
  }
  return $This;
}

# Set atom paris count for a specific atom ID pair at a specific distance...
#
sub _SetAtomPairsCount {
  my($This, $Distance, $AtomType1, $AtomType2) = @_;

  if (! exists $This->{AtomPairsCount}{$Distance}{$AtomType1}) {
    %{$This->{AtomPairsCount}{$Distance}{$AtomType1}} = ();
    $This->{AtomPairsCount}{$Distance}{$AtomType1}{$AtomType2} = 1;
    return $This;
  }

  if (exists $This->{AtomPairsCount}{$Distance}{$AtomType1}{$AtomType2}) {
    $This->{AtomPairsCount}{$Distance}{$AtomType1}{$AtomType2} += 1;
  }
  else {
    $This->{AtomPairsCount}{$Distance}{$AtomType1}{$AtomType2} = 1;
  }

  return $This;
}

# Set final fingerpritns vector...
#
sub _SetFinalFingerprints {
  my($This) = @_;
  my($Distance, $AtomType1, $AtomType2, $Value, @Values);

  # Mark successful generation of fingerprints...
  $This->{FingerprintsGenerated} = 1;

  @Values = ();
  @{$This->{AtomPairsIDs}} = ();

  for $Distance ($This->{MinDistance} .. $This->{MaxDistance}) {
    for $AtomType1 (sort keys %{$This->{AtomPairsCount}{$Distance}} ) {
      for $AtomType2 (sort keys %{$This->{AtomPairsCount}{$Distance}{$AtomType1}} ) {
	push @{$This->{AtomPairsIDs}}, "${AtomType1}-D${Distance}-${AtomType2}";
	$Value = $This->{AtomPairsCount}{$Distance}{$AtomType1}{$AtomType2};
	push @Values, $Value;
      }
    }
  }

  # Add AtomPairsIDs and values to fingerprint vector...
  $This->{FingerprintsVector}->AddValueIDs(\@{$This->{AtomPairsIDs}});
  $This->{FingerprintsVector}->AddValues(\@Values);

  return $This;
}

# Get atom pair IDs corresponding to atom pairs count values in fingerprint
# vector as an array or reference to an array...
#
# AtomPairIDs list differes in molecules and is generated during finalization
# of fingerprints to make sure the fingerprint vector containing count values
# matches the atom pairs array.
#
sub GetAtomPairIDs {
  my($This) = @_;

  return wantarray ? @{$This->{AtomPairsIDs}} : \@{$This->{AtomPairsIDs}};
}

# Cache  appropriate molecule data...
#
sub _SetupMoleculeDataCache {
  my($This) = @_;

  # Get all atoms including hydrogens to correctly map atom indices to atom IDs for
  # usage of distance matrix. The hydrogen atoms are ignored during processing...
  #
  @{$This->{Atoms}} = $This->GetMolecule()->GetAtoms();

  # Get all atom IDs...
  my(@AtomIDs);
  @AtomIDs = ();
  @AtomIDs =  map { $_->GetID() } @{$This->{Atoms}};

  # Set AtomIndex to AtomID hash...
  %{$This->{AtomIndexToID}} = ();
  @{$This->{AtomIndexToID}}{ (0 .. $#AtomIDs) } = @AtomIDs;

  return $This;
}

# Clear cached molecule data...
#
sub _ClearMoleculeDataCache {
  my($This) = @_;

  @{$This->{Atoms}} = ();

  return $This;
}

# Set atomic invariants to use for atom identifiers...
#
sub SetAtomicInvariantsToUse {
  my($This, @Values) = @_;
  my($FirstValue, $TypeOfFirstValue, $AtomicInvariant, $SpecifiedAtomicInvariant, $AtomicInvariantValue, @SpecifiedAtomicInvariants, @AtomicInvariantsToUse);

  if (!@Values) {
    carp "Warning: ${ClassName}->SetAtomicInvariantsToUse: No values specified...";
    return;
  }

  $FirstValue = $Values[0];
  $TypeOfFirstValue = ref $FirstValue;

  @SpecifiedAtomicInvariants = ();
  @AtomicInvariantsToUse = ();

  if ($TypeOfFirstValue =~ /^ARRAY/) {
    push @SpecifiedAtomicInvariants, @{$FirstValue};
  }
  else {
    push @SpecifiedAtomicInvariants, @Values;
  }

  # Make sure specified AtomicInvariants are valid...
  for $SpecifiedAtomicInvariant (@SpecifiedAtomicInvariants) {
    if (!AtomTypes::AtomicInvariantsAtomTypes::IsAtomicInvariantAvailable($SpecifiedAtomicInvariant)) {
      croak "Error: ${ClassName}->SetAtomicInvariantsToUse: Specified atomic invariant, $SpecifiedAtomicInvariant, is not supported...\n ";
    }
    $AtomicInvariant = $SpecifiedAtomicInvariant;
    push @AtomicInvariantsToUse, $AtomicInvariant;
  }

  # Set atomic invariants to use...
  @{$This->{AtomicInvariantsToUse}} = ();
  push @{$This->{AtomicInvariantsToUse}}, @AtomicInvariantsToUse;

  return $This;
}

# Set functional classes to use for atom identifiers...
#
sub SetFunctionalClassesToUse {
  my($This, @Values) = @_;
  my($FirstValue, $TypeOfFirstValue, $FunctionalClass, $SpecifiedFunctionalClass, @SpecifiedFunctionalClasses, @FunctionalClassesToUse);

  if (!@Values) {
    carp "Warning: ${ClassName}->SetFunctionalClassesToUse: No values specified...";
    return;
  }

  if ($This->{AtomIdentifierType} !~ /^FunctionalClassAtomTypes$/i) {
    carp "Warning: ${ClassName}->SetFunctionalClassesToUse: FunctionalClassesToUse can't be set for InitialAtomIdentifierType of $This->{AtomIdentifierType}...";
    return;
  }

  $FirstValue = $Values[0];
  $TypeOfFirstValue = ref $FirstValue;

  @SpecifiedFunctionalClasses = ();
  @FunctionalClassesToUse = ();

  if ($TypeOfFirstValue =~ /^ARRAY/) {
    push @SpecifiedFunctionalClasses, @{$FirstValue};
  }
  else {
    push @SpecifiedFunctionalClasses, @Values;
  }

  # Make sure specified FunctionalClasses are valid...
  for $SpecifiedFunctionalClass (@SpecifiedFunctionalClasses) {
    if (!AtomTypes::FunctionalClassAtomTypes::IsFunctionalClassAvailable($SpecifiedFunctionalClass)) {
      croak "Error: ${ClassName}->SetFunctionalClassesToUse: Specified functional class, $SpecifiedFunctionalClass, is not supported...\n ";
    }
    push @FunctionalClassesToUse, $SpecifiedFunctionalClass;
  }

  # Set functional classes to use...
  @{$This->{FunctionalClassesToUse}} = ();
  push @{$This->{FunctionalClassesToUse}}, @FunctionalClassesToUse;

  return $This;
}

# Initialize atom indentifier type information...
#
# Current supported values:
#
# AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes, FunctionalClassAtomTypes,
# MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes
#
sub _InitializeAtomIdentifierTypeInformation {
  my($This) = @_;

  if ($This->{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
    $This->_InitializeAtomicInvariantsAtomTypesInformation();
  }
  elsif ($This->{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
    $This->_InitializeFunctionalClassAtomTypesInformation();
  }
  elsif ($This->{AtomIdentifierType} =~ /^(DREIDINGAtomTypes|EStateAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
    # Nothing to do for now...
  }
  else {
    croak "Error: ${ClassName}->_InitializeAtomIdentifierTypeInformation: Unknown atom indentifier type $This->{AtomIdentifierType}...";
  }

  return $This;
}

# Initialize atomic invariants atom types to use for generating atom identifiers...
#
# Let:
#   AS = Atom symbol corresponding to element symbol
#
#   X<n>   = Number of non-hydrogen atom neighbors or heavy atoms attached to atom
#   BO<n> = Sum of bond orders to non-hydrogen atom neighbors or heavy atoms attached to atom
#   LBO<n> = Largest bond order of non-hydrogen atom neighbors or heavy atoms attached to atom
#   SB<n> = Number of single bonds to non-hydrogen atom neighbors or heavy atoms attached to atom
#   DB<n> = Number of double bonds to non-hydrogen atom neighbors or heavy atoms attached to atom
#   TB<n> = Number of triple bonds to non-hydrogen atom neighbors or heavy atoms attached to atom
#   H<n>   = Number of implicit and explicit hydrogens for atom
#   Ar     = Aromatic annotation indicating whether atom is aromatic
#   RA     = Ring atom annotation indicating whether atom is a ring
#   FC<+n/-n> = Formal charge assigned to atom
#   MN<n> = Mass number indicating isotope other than most abundant isotope
#   SM<n> = Spin multiplicity of atom. Possible values: 1 (singlet), 2 (doublet) or 3 (triplet)
#
#   AtomTypeIDx = Atomic invariants atom type for atom x
#   AtomTypeIDy = Atomic invariants atom type for atom y
#   Dn   = Topological distance between atom x and y
#
# Then:
#
#   Atom pair AtomID generated by AtomTypes::AtomicInvariantsAtomTypes class corresponds to:
#
#     AS.X<n>.BO<n>.LBO<n>.<SB><n>.<DB><n>.<TB><n>.H<n>.Ar.RA.FC<+n/-n>.MN<n>.SM<n>
#
#  AtomPairID corresponds to:
#
#    AtomTypeIDx-D<n>-AtomTypeIDy
#
# Except for AS which is a required atomic invariant in atom pair AtomIDs, all other atomic invariants are
# optional. Default atomic invariants used for AtomID are: AS, X<n>, BO<n>, H<n>, FC<+n/-n>.
# AtomID specification doesn't include atomic invariants with zero or undefined values.
#
# Examples of atom pair AtomIDs:
#
#   O.X1.BO1.H1 - Hydroxyl oxygen in carboxylate with attached hydrogen and no explicit charge
#   O.X1.BO1.FC-1 - Hydroxyl ozygen in carboxylate with explicit negative charge
#   O.X1.BO2 - Carbonyl oxygen in carboxylate with double bond to carbon
#   O.X2.BO2 - Hydroxyl ozygen in carboxylate attached to carbonyl carbon and another heavy atom
#
#   C.X2.BO3.H1.Ar - Aromatic carbon
#
# Examples of AtomPairIDs:
#
#   C.X2.BO2.H3-D1-O.X1.BO1 - Carbon with two heavy atom neighbors attached to oxygen at bond distance 1(methanol)
#
#   C.X2.BO3.H1.Ar-D3-C.X2.BO3.H1.Ar  - Two aromatic carbons at bond distance 3 where each carbon has
#                                       two heavy atom neighbors and bond order of 3 (benzene)
#
sub _InitializeAtomicInvariantsAtomTypesInformation {
  my($This) = @_;

  # Default atomic invariants to use for generating atom neighborhood atom IDs: AS, X, BO, H, FC
  #
  @{$This->{AtomicInvariantsToUse}} = ();
  @{$This->{AtomicInvariantsToUse}} = ('AS', 'X', 'BO', 'H', 'FC');

  return $This;
}

# Initialize functional class atom types, generated by AtomTypes::FunctionalClassAtomTypes
# class, to use for generating atom identifiers...
#
# Let:
#   HBD: HydrogenBondDonor
#   HBA: HydrogenBondAcceptor
#   PI :  PositivelyIonizable
#   NI : NegativelyIonizable
#   Ar : Aromatic
#   Hal : Halogen
#   H : Hydrophobic
#   RA : RingAtom
#   CA : ChainAtom
#
# Then:
#
#   Functiononal class atom type specification for an atom corresponds to:
#
#     Ar.CA.H.HBA.HBD.Hal.NI.PI.RA
#
#   Default functional classes used are: HBD, HBA, PI, NI, Ar, Hal
#
#   FunctionalAtomTypes are assigned using the following definitions [ Ref 60-61, Ref 65-66 ]:
#
#     HydrogenBondDonor: NH, NH2, OH
#     HydrogenBondAcceptor: N[!H], O
#     PositivelyIonizable: +, NH2
#     NegativelyIonizable: -, C(=O)OH, S(=O)OH, P(=O)OH
#
sub _InitializeFunctionalClassAtomTypesInformation {
  my($This) = @_;

  # Default functional class atom typess to use for generating atom identifiers
  # are: HBD, HBA, PI, NI, Ar, Hal
  #
  @{$This->{FunctionalClassesToUse}} = ();
  @{$This->{FunctionalClassesToUse}} = ('HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal');

  return $This;
}

# Return a string containg data for TopologicalAtomPairsFingerprints object...
#
sub StringifyTopologicalAtomPairsFingerprints {
  my($This) = @_;
  my($FingerprintsString);

  # Type of fingerprint...
  $FingerprintsString = "Fingerprint type: $This->{Type}; AtomIdentifierType: $This->{AtomIdentifierType}";

  # Min and max distance...
  $FingerprintsString .= "; MinDistance:  $This->{MinDistance}; MaxDistance: $This->{MaxDistance}";

  if ($This->{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
    my($AtomicInvariant, @AtomicInvariants, @AtomicInvariantsOrder, %AvailableAtomicInvariants);

    @AtomicInvariantsOrder = AtomTypes::AtomicInvariantsAtomTypes::GetAtomicInvariantsOrder();
    %AvailableAtomicInvariants = AtomTypes::AtomicInvariantsAtomTypes::GetAvailableAtomicInvariants();

    for $AtomicInvariant (@AtomicInvariantsOrder) {
      push @AtomicInvariants, "$AtomicInvariant: $AvailableAtomicInvariants{$AtomicInvariant}";
    }

    $FingerprintsString .= "; AtomicInvariantsToUse: <" . TextUtil::JoinWords(\@{$This->{AtomicInvariantsToUse}}, ", ", 0) . ">";
    $FingerprintsString .= "; AtomicInvariantsOrder: <" . TextUtil::JoinWords(\@AtomicInvariantsOrder, ", ", 0) . ">";
    $FingerprintsString .= "; AvailableAtomicInvariants: <" . TextUtil::JoinWords(\@AtomicInvariants, ", ", 0) . ">";
  }
  elsif ($This->{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
    my($FunctionalClass, @FunctionalClasses, @FunctionalClassesOrder, %AvailableFunctionalClasses);

    @FunctionalClassesOrder = AtomTypes::FunctionalClassAtomTypes::GetFunctionalClassesOrder();
    %AvailableFunctionalClasses = AtomTypes::FunctionalClassAtomTypes::GetAvailableFunctionalClasses();

    for $FunctionalClass (@FunctionalClassesOrder) {
      push @FunctionalClasses, "$FunctionalClass: $AvailableFunctionalClasses{$FunctionalClass}";
    }

    $FingerprintsString .= "; FunctionalClassesToUse: <" . TextUtil::JoinWords(\@{$This->{FunctionalClassesToUse}}, ", ", 0) . ">";
    $FingerprintsString .= "; FunctionalClassesOrder: <" . TextUtil::JoinWords(\@FunctionalClassesOrder, ", ", 0) . ">";
    $FingerprintsString .= "; AvailableFunctionalClasses: <" . TextUtil::JoinWords(\@FunctionalClasses, ", ", 0) . ">";
  }

  # Total number of atom pairs...
  $FingerprintsString .= "; NumOfAtomPairs: " . $This->{FingerprintsVector}->GetNumOfValues();

  # FingerprintsVector...
  $FingerprintsString .= "; FingerprintsVector: < $This->{FingerprintsVector} >";

  return $FingerprintsString;
}

1;

__END__

=head1 NAME

TopologicalAtomPairsFingerprints

=head1 SYNOPSIS

use Fingerprints::TopologicalAtomPairsFingerprints;

use Fingerprints::TopologicalAtomPairsFingerprints qw(:all);

=head1 DESCRIPTION

B<TopologicalAtomPairsFingerprints>  [ Ref 57, Ref 59, Ref 72 ] class provides the following methods:

new, GenerateFingerprints, GetAtomPairIDs, GetDescription, SetAtomIdentifierType,
SetAtomicInvariantsToUse, SetFunctionalClassesToUse, SetMaxDistance,
SetMinDistance, StringifyTopologicalAtomPairsFingerprints

B<TopologicalAtomPairsFingerprints> is derived from B<Fingerprints> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<TopologicalAtomPairsFingerprints>, B<Fingerprints> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The current release of MayaChemTools supports generation of B<AtomTypesFingerpritns>
corresponding to following B<AtomtomIdentifierTypes>:

    AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
    FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
    SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes

Based on the values specified for B<AtomIdentifierType> along with other specified
parameters such as B<AtomicInvariantsToUse> and B<FunctionalClassesToUse>, initial
atom types are assigned to all non-hydrogen atoms in a molecule. Using the distance
matrix for the molecule and initial atom types assigned to non-hydrogen atoms, all unique atom
pairs within B<MinDistance> and B<MaxDistance> are identified and counted. An atom pair
identifier is generated for each unique atom pair; the format of atom pair identifier is:

    <AtomType1>-D<n>-<AtomType2>

    AtomType1, AtomType2: Atom types assigned to atom1 and atom2
    D: Distance between atom1 and atom2

    where AtomType1 <= AtomType2

The atom pair identifiers for all unique atom pairs corresponding to non-hydrogen atoms constitute
topological atom pairs fingerprints of the molecule.

The current release of MayaChemTools generates the following types of topological atom pairs
fingerprints vector strings:

    FingerprintsVector;TopologicalAtomPairs:AtomicInvariantsAtomTypes:MinD
    istance1:MaxDistance10;223;NumericalValues;IDsAndValuesString;C.X1.BO1
    .H3-D1-C.X3.BO3.H1 C.X2.BO2.H2-D1-C.X2.BO2.H2 C.X2.BO2.H2-D1-C.X3.BO3.
    H1 C.X2.BO2.H2-D1-C.X3.BO4 C.X2.BO2.H2-D1-N.X3.BO3 C.X2.BO3.H1-D1-...;
    2 1 4 1 1 10 8 1 2 6 1 2 2 1 2 1 2 2 1 2 1 5 1 10 12 2 2 1 2 1 9 1 3 1
    1 1 2 2 1 3 6 1 6 14 2 2 2 3 1 3 1 8 2 2 1 3 2 6 1 2 2 5 1 3 1 23 1...

    FingerprintsVector;TopologicalAtomPairs:AtomicInvariantsAtomTypes:MinD
    istance1:MaxDistance10;223;NumericalValues;IDsAndValuesPairsString;C.X
    1.BO1.H3-D1-C.X3.BO3.H1 2 C.X2.BO2.H2-D1-C.X2.BO2.H2 1 C.X2.BO2.H2-D1-
    C.X3.BO3.H1 4 C.X2.BO2.H2-D1-C.X3.BO4 1 C.X2.BO2.H2-D1-N.X3.BO3 1 C.X2
    .BO3.H1-D1-C.X2.BO3.H1 10 C.X2.BO3.H1-D1-C.X3.BO4 8 C.X3.BO3.H1-D1-C.X
    3.BO4 1 C.X3.BO3.H1-D1-O.X1.BO1.H1 2 C.X3.BO4-D1-C.X3.BO4 6 C.X3.BO...

    FingerprintsVector;TopologicalAtomPairs:DREIDINGAtomTypes:MinDistance1
    :MaxDistance10;157;NumericalValues;IDsAndValuesString;C_2-D1-C_3 C_2-D
    1-C_R C_2-D1-N_3 C_2-D1-O_2 C_2-D1-O_3 C_3-D1-C_3 C_3-D1-C_R C_3-D1-N_
    R C_3-D1-O_3 C_R-D1-C_R C_R-D1-F_ C_R-D1-N_3 C_R-D1-N_R C_2-D2-C_3 C_2
    1 1 1 2 1 7 1 1 2 23 1 1 2 1 3 5 5 2 1 5 28 2 3 3 1 1 1 2 4 1 1 4 9 3
    1 4 24 2 4 3 3 4 5 5 14 1 1 2 3 22 1 3 4 4 1 1 1 1 2 2 5 1 4 21 3 1...

    FingerprintsVector;TopologicalAtomPairs:EStateAtomTypes:MinDistance1:M
    axDistance10;251;NumericalValues;IDsAndValuesString;aaCH-D1-aaCH aaCH-
    D1-aasC aasC-D1-aasC aasC-D1-aasN aasC-D1-dssC aasC-D1-sF aasC-D1-ssNH
    aasC-D1-sssCH aasN-D1-ssCH2 dO-D1-dssC dssC-D1-sOH dssC-D1-ssCH2 d...;
    10 8 5 2 1 1 1 1 1 2 1 1 1 2 2 1 4 10 12 2 2 6 3 1 3 2 2 1 1 1 1 1 1 1
    1 1 5 2 1 1 6 12 2 2 2 2 6 1 3 2 2 5 2 2 1 2 1 1 1 1 1 1 3 1 3 19 2...

    FingerprintsVector;TopologicalAtomPairs:FunctionalClassAtomTypes:MinDi
    stance1:MaxDistance10;144;NumericalValues;IDsAndValuesString;Ar-D1-Ar
    Ar-D1-Ar.HBA Ar-D1-HBD Ar-D1-Hal Ar-D1-None Ar.HBA-D1-None HBA-D1-NI H
    BA-D1-None HBA.HBD-D1-NI HBA.HBD-D1-None HBD-D1-None NI-D1-None No...;
    23 2 1 1 2 1 1 1 1 2 1 1 7 28 3 1 3 2 8 2 1 1 1 5 1 5 24 3 3 4 2 13 4
    1 1 4 1 5 22 4 4 3 1 19 1 1 1 1 1 2 2 3 1 1 8 25 4 5 2 3 1 26 1 4 1 ...

    FingerprintsVector;TopologicalAtomPairs:MMFF94AtomTypes:MinDistance1:M
    axDistance10;227;NumericalValues;IDsAndValuesPairsString;C5A-D1-C5B 2 
    C5A-D1-CB 1 C5A-D1-CR 1 C5A-D1-N5 2 C5B-D1-C5B 1 C5B-D1-C=ON 1 C5B-D1-
    CB 1 C=ON-D1-NC=O 1 C=ON-D1-O=CN 1 CB-D1-CB 18 CB-D1-F 1 CB-D1-NC=O 1
    COO-D1-CR 1 COO-D1-O=CO 1 COO-D1-OC=O 1 CR-D1-CR 7 CR-D1-N5 1 CR-D1-OR
    2 C5A-D2-C5A 1 C5A-D2-C5B 2 C5A-D2-C=ON 1 C5A-D2-CB 3 C5A-D2-CR 4 ...

    FingerprintsVector;TopologicalAtomPairs:SLogPAtomTypes:MinDistance1:Ma
    xDistance10;329;NumericalValues;IDsAndValuesPairsString;C1-D1-C10 1 C1
    -D1-C11 2 C1-D1-C5 1 C1-D1-CS 4 C10-D1-N11 1 C11-D1-C21 1 C14-D1-C18 2
    C14-D1-F 1 C18-D1-C18 10 C18-D1-C20 4 C18-D1-C22 2 C20-D1-C20 3 C20-D
    1-C21 1 C20-D1-N11 1 C21-D1-C21 1 C21-D1-C5 1 C21-D1-N11 1 C22-D1-N4 1
    C5-D1-N4 1 C5-D1-O10 1 C5-D1-O2 1 C5-D1-O9 1 CS-D1-O2 2 C1-D2-C1 3...

    FingerprintsVector;TopologicalAtomPairs:SYBYLAtomTypes:MinDistance1:Ma
    xDistance10;159;NumericalValues;IDsAndValuesPairsString;C.2-D1-C.3 1 C
    .2-D1-C.ar 1 C.2-D1-N.am 1 C.2-D1-O.2 1 C.2-D1-O.co2 2 C.3-D1-C.3 7 C.
    3-D1-C.ar 1 C.3-D1-N.ar 1 C.3-D1-O.3 2 C.ar-D1-C.ar 23 C.ar-D1-F 1 C.a
    r-D1-N.am 1 C.ar-D1-N.ar 2 C.2-D2-C.3 1 C.2-D2-C.ar 3 C.3-D2-C.3 5 C.3
    -D2-C.ar 5 C.3-D2-N.ar 2 C.3-D2-O.3 4 C.3-D2-O.co2 2 C.ar-D2-C.ar 2...

    FingerprintsVector;TopologicalAtomPairs:TPSAAtomTypes:MinDistance1:Max
    Distance10;64;NumericalValues;IDsAndValuesPairsString;N21-D1-None 3 N7
    -D1-None 2 None-D1-None 34 None-D1-O3 2 None-D1-O4 3 N21-D2-None 5 N7-
    D2-None 3 N7-D2-O3 1 None-D2-None 44 None-D2-O3 2 None-D2-O4 5 O3-D2-O
    4 1 N21-D3-None 7 N7-D3-None 4 None-D3-None 45 None-D3-O3 4 None-D3-O4
    5 N21-D4-N7 1 N21-D4-None 5 N21-D4-O3 1 N21-D4-O4 1 N7-D4-None 4 N...

    FingerprintsVector;TopologicalAtomPairs:UFFAtomTypes:MinDistance1:MaxD
    istance10;157;NumericalValues;IDsAndValuesPairsString;C_2-D1-C_3 1 C_2
    -D1-C_R 1 C_2-D1-N_3 1 C_2-D1-O_2 2 C_2-D1-O_3 1 C_3-D1-C_3 7 C_3-D1-C
    _R 1 C_3-D1-N_R 1 C_3-D1-O_3 2 C_R-D1-C_R 23 C_R-D1-F_ 1 C_R-D1-N_3 1 
    C_R-D1-N_R 2 C_2-D2-C_3 1 C_2-D2-C_R 3 C_3-D2-C_3 5 C_3-D2-C_R 5 C_3-D
    2-N_R 2 C_3-D2-O_2 1 C_3-D2-O_3 5 C_R-D2-C_R 28 C_R-D2-F_ 2 C_R-D2-...

=head2 METHODS

=over 4

=item B<new>

    $NewTopologicalAtomPairsFingerprints = new TopologicalAtomPairsFingerprints(
                                                   %NamesAndValues);

Using specified I<TopologicalAtomPairsFingerprints> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<TopologicalAtomPairsFingerprints>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'TopologicalAtomPairs'
    MinDistance = 1
    MaxDistance = 10
    AtomIdentifierType = ''
    AtomicInvariantsToUse = ['AS', 'X', 'BO', 'H', 'FC']
    FunctionalClassesToUse = ['HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal']

Examples:

    $TopologicalAtomPairsFingerprints = new TopologicalAtomPairsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes');

    $TopologicalAtomPairsFingerprints = new TopologicalAtomPairsFingerprints(
                              'Molecule' => $Molecule,
                              'MinDistance' => 1,
                              'MaxDistance' => 10,
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes',
                              'AtomicInvariantsToUse' =>
                                              ['AS', 'X', 'BO', 'H', 'FC'] );

    $TopologicalAtomPairsFingerprints = new TopologicalAtomPairsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'EStateAtomTypes');

    $TopologicalAtomPairsFingerprints = new TopologicalAtomPairsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'SLogPAtomTypes');

    $TopologicalAtomPairsFingerprints = new TopologicalAtomPairsFingerprints(
                              'Molecule' => $Molecule,
                              'MinDistance' => 1,
                              'MaxDistance' => 10,
                              'AtomIdentifierType' =>
                                              'FunctionalClassAtomTypes',
                              'FunctionalClassesToUse' =>
                                              ['HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal']);


    $TopologicalAtomPairsFingerprints->GenerateFingerprints();
    print "$TopologicalAtomPairsFingerprints\n";

=item B<GetDescription>

    $Description = $TopologicalAtomPairsFingerprints->GetDescription();

Returns a string containing description of topological atom pairs fingerprints fingerprints.

=item B<GenerateFingerprints>

    $TopologicalAtomPairsFingerprints->GenerateFingerprints();

Generates topological atom pairs fingerprints and returns I<TopologicalAtomPairsFingerprints>.

=item B<GetAtomPairIDs>

    $AtomPairIDsRef = $TopologicalAtomPairsFingerprints->GetAtomPairIDs();
    @AtomPairIDs = $TopologicalAtomPairsFingerprints->GetAtomPairIDs();

Returns atom pair IDs corresponding to atom pairs count values in topological atom pairs
fingerprints vector as an array or reference to an array.

=item B<SetAtomIdentifierType>

    $TopologicalAtomPairsFingerprints->SetAtomIdentifierType($IdentifierType);

Sets atom I<IdentifierType> to use during atom pairs fingerprints generation and
returns I<TopologicalAtomPairsFingerprints>.

Possible values: I<AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes,
TPSAAtomTypes, UFFAtomTypes>.

=item B<SetAtomicInvariantsToUse>

    $TopologicalAtomPairsFingerprints->SetAtomicInvariantsToUse($ValuesRef);
    $TopologicalAtomPairsFingerprints->SetAtomicInvariantsToUse(@Values);

Sets atomic invariants to use during I<AtomicInvariantsAtomTypes> value of I<AtomIdentifierType>
for topological atom pairs fingerprints generation and returns I<TopologicalAtomPairsFingerprints>.

Possible values for atomic invariants are: I<AS, X, BO,  LBO, SB, DB, TB,
H, Ar, RA, FC, MN, SM>. Default value: I<AS,X,BO,H,FC>.

The atomic invariants abbreviations correspond to:

    AS = Atom symbol corresponding to element symbol

    X<n>   = Number of non-hydrogen atom neighbors or heavy atoms
    BO<n> = Sum of bond orders to non-hydrogen atom neighbors or heavy atoms
    LBO<n> = Largest bond order of non-hydrogen atom neighbors or heavy atoms
    SB<n> = Number of single bonds to non-hydrogen atom neighbors or heavy atoms
    DB<n> = Number of double bonds to non-hydrogen atom neighbors or heavy atoms
    TB<n> = Number of triple bonds to non-hydrogen atom neighbors or heavy atoms
    H<n>   = Number of implicit and explicit hydrogens for atom
    Ar     = Aromatic annotation indicating whether atom is aromatic
    RA     = Ring atom annotation indicating whether atom is a ring
    FC<+n/-n> = Formal charge assigned to atom
    MN<n> = Mass number indicating isotope other than most abundant isotope
    SM<n> = Spin multiplicity of atom. Possible values: 1 (singlet), 2 (doublet) or
            3 (triplet)

Atom type generated by AtomTypes::AtomicInvariantsAtomTypes class corresponds to:

    AS.X<n>.BO<n>.LBO<n>.<SB><n>.<DB><n>.<TB><n>.H<n>.Ar.RA.FC<+n/-n>.MN<n>.SM<n>

Except for AS which is a required atomic invariant in atom types, all other atomic invariants are
optional. Atom type specification doesn't include atomic invariants with zero or undefined values.

In addition to usage of abbreviations for specifying atomic invariants, the following descriptive words
are also allowed:

    X : NumOfNonHydrogenAtomNeighbors or NumOfHeavyAtomNeighbors
    BO : SumOfBondOrdersToNonHydrogenAtoms or SumOfBondOrdersToHeavyAtoms
    LBO : LargestBondOrderToNonHydrogenAtoms or LargestBondOrderToHeavyAtoms
    SB :  NumOfSingleBondsToNonHydrogenAtoms or NumOfSingleBondsToHeavyAtoms
    DB : NumOfDoubleBondsToNonHydrogenAtoms or NumOfDoubleBondsToHeavyAtoms
    TB : NumOfTripleBondsToNonHydrogenAtoms or NumOfTripleBondsToHeavyAtoms
    H :  NumOfImplicitAndExplicitHydrogens
    Ar : Aromatic
    RA : RingAtom
    FC : FormalCharge
    MN : MassNumber
    SM : SpinMultiplicity

I<AtomTypes::AtomicInvariantsAtomTypes> module is used to assign atomic invariant
atom types.

=item B<SetFunctionalClassesToUse>

    $TopologicalAtomPairsFingerprints->SetFunctionalClassesToUse($ValuesRef);
    $TopologicalAtomPairsFingerprints->SetFunctionalClassesToUse(@Values);

Sets functional classes invariants to use during I<FunctionalClassAtomTypes> value of I<AtomIdentifierType>
for topological atom pairs fingerprints generation and returns I<TopologicalAtomPairsFingerprints>.

Possible values for atom functional classes are: I<Ar, CA, H, HBA, HBD, Hal, NI, PI, RA>.
Default value [ Ref 24 ]: I<HBD,HBA,PI,NI,Ar,Hal>.

The functional class abbreviations correspond to:

    HBD: HydrogenBondDonor
    HBA: HydrogenBondAcceptor
    PI :  PositivelyIonizable
    NI : NegativelyIonizable
    Ar : Aromatic
    Hal : Halogen
    H : Hydrophobic
    RA : RingAtom
    CA : ChainAtom

 Functional class atom type specification for an atom corresponds to:

    Ar.CA.H.HBA.HBD.Hal.NI.PI.RA or None

I<AtomTypes::FunctionalClassAtomTypes> module is used to assign functional class atom
types. It uses following definitions [ Ref 60-61, Ref 65-66 ]:

    HydrogenBondDonor: NH, NH2, OH
    HydrogenBondAcceptor: N[!H], O
    PositivelyIonizable: +, NH2
    NegativelyIonizable: -, C(=O)OH, S(=O)OH, P(=O)OH

=item B<SetMaxDistance>

    $TopologicalAtomPairsFingerprints->SetMaxDistance($Distance);

Sets maximum distance to use during topological atom pairs fingerprints generation and
returns I<TopologicalAtomPairsFingerprints>.

=item B<SetMinDistance>

    $TopologicalAtomPairsFingerprints->SetMinDistance($Distance);

Sets minimum distance to use during topological atom pairs fingerprints generation and
returns I<TopologicalAtomPairsFingerprints>.

=item B<StringifyTopologicalAtomPairsFingerprints>

    $String = $TopologicalAtomPairsFingerprints->
                  StringifyTopologicalAtomPairsFingerprints();

Returns a string containing information about I<TopologicalAtomPairsFingerprints> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Fingerprints.pm, FingerprintsStringUtil.pm, AtomNeighborhoodsFingerprints.pm,
AtomTypesFingerprints.pm, EStateIndiciesFingerprints.pm, ExtendedConnectivityFingerprints.pm,
MACCSKeys.pm, PathLengthFingerprints.pm, TopologicalAtomTripletsFingerprints.pm,
TopologicalAtomTorsionsFingerprints.pm, TopologicalPharmacophoreAtomPairsFingerprints.pm,
TopologicalPharmacophoreAtomTripletsFingerprints.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
