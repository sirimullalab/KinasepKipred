package Fingerprints::TopologicalAtomTripletsFingerprints;
#
# File: TopologicalAtomTripletsFingerprints.pm
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
use overload '""' => 'StringifyTopologicalAtomTripletsFingerprints';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeTopologicalAtomTripletsFingerprints();

  $This->_InitializeTopologicalAtomTripletsFingerprintsProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeTopologicalAtomTripletsFingerprints {
  my($This) = @_;

  # Type of fingerprint...
  $This->{Type} = 'TopologicalAtomTriplets';

  # Type of vector...
  $This->{VectorType} = 'FingerprintsVector';

  # Type of FingerprintsVector...
  $This->{FingerprintsVectorType} = 'NumericalValues';

  # Minimum and maximum bond distance between atom paris...
  $This->{MinDistance} = 1;
  $This->{MaxDistance} = 10;

  # Determines whether to apply triangle inequality to distance triplets...
  #
  $This->{UseTriangleInequality} = 0;

  # Atom identifier type to use for atom IDs in atom triplets...
  #
  # Currently supported values are: AtomicInvariantsAtomTypes, DREIDINGAtomTypes,
  # EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
  # SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes
  #
  $This->{AtomIdentifierType} = '';

  # Atom types assigned to each heavy atom...
  #
  %{$This->{AssignedAtomTypes}} = ();

  # All atom triplets between minimum and maximum distance...
  #
  @{$This->{AtomTripletsIDs}} = ();
  %{$This->{AtomTripletsCount}} = ();
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize object properties....
sub _InitializeTopologicalAtomTripletsFingerprintsProperties {
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

# Set minimum distance for atom triplets...
#
sub SetMinDistance {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetMinDistance: MinDistance value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{MinDistance} = $Value;

  return $This;
}

# Set maximum distance for atom triplets...
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

# Generate topological atom triplets fingerprints...
#
# Let:
#
#   AT = Any of the supported atom types
#
#   ATx = Atom type for  atom x
#   ATy = Atom type for  atom y
#   ATz = Atom type for  atom z
#
#   Dxy = Distance between Px and Py
#   Dxz = Distance between Px and Pz
#   Dyz = Distance between Py and Pz
#
# Then:
#
#   ATx-Dyz-ATy-Dxz-ATz-Dxy = Atom triplet ID for atom types ATx, ATy and Atz
#
# Methodology:
#   . Generate a distance matrix.
#   . Assign atom types to all the atoms.
#   . Using distance matrix and atom types, count occurrence of unique atom triplets
#     within specified distance range along with optional trinagle inequality
#
# Notes:
#   . Hydrogen atoms are ignored during the fingerprint generation.
#   . For a molecule containing N atoms with all different atom type, the total number of
#     possible unique atom triplets without applying triangle inquality check corresponds to:
#
#     Factorial( N ) / ( Factorial( N - 3 ) * Factorial (3) )
#
#     However, due to similar atom types assigned to atoms in a molecule for a specific atom
#     typing methodology and specified distance range used during fingerprints generation, the
#     actual number of unique triplets is usually smaller than the theoretical limit.
#
sub GenerateFingerprints {
  my($This) = @_;

  if ($This->{MinDistance} > $This->{MaxDistance}) {
    croak "Error: ${ClassName}->GenerateTopologicalAtomTripletsFingerprints: No fingerpritns generated: MinDistance, $This->{MinDistance}, must be <= MaxDistance, $This->{MaxDistance}...";
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

  # Intialize values of toplogical atom triplets...
  $This->_InitializeToplogicalAtomTriplets();

  # Count atom triplets...
  $This->_GenerateAndCountAtomTriplets();

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

# Initialize topological atom triplets between specified distance range...
#
sub _InitializeToplogicalAtomTriplets {
  my($This) = @_;
  my($Distance);

  @{$This->{AtomTripletsIDs}} = ();
  %{$This->{AtomTripletsCount}} = ();

  return $This;
}

# Count atom triplets between mininum and maximum distance at each
# distance using distance matrix and atom types assiged to each heavy
# atom.
#
sub _GenerateAndCountAtomTriplets {
  my($This) = @_;
  my($NumOfAtoms, $AtomIndex1, $AtomIndex2, $AtomIndex3, $AtomID1, $AtomID2, $AtomID3, $AtomType1, $AtomType2, $AtomType3, $Distance12, $Distance13, $Distance23, $SkipIndexCheck, $DistanceMatrix, $AtomTripletID);

  $NumOfAtoms = @{$This->{Atoms}};
  $DistanceMatrix = $This->{DistanceMatrix};
  $SkipIndexCheck = 0;

  ATOMINDEX1: for $AtomIndex1 (0 .. ($NumOfAtoms - 1)) {
    $AtomID1 = $This->{AtomIndexToID}{$AtomIndex1};
    if (!exists($This->{AssignedAtomTypes}{$AtomID1})) {
      next ATOMINDEX1;
    }
    $AtomType1 = $This->{AssignedAtomTypes}{$AtomID1};

    ATOMINDEX2: for $AtomIndex2 (($AtomIndex1 + 1) .. ($NumOfAtoms - 1)) {
      $AtomID2 = $This->{AtomIndexToID}{$AtomIndex2};
      if (!exists($This->{AssignedAtomTypes}{$AtomID2})) {
	next ATOMINDEX2;
      }
      $AtomType2 = $This->{AssignedAtomTypes}{$AtomID2};

      $Distance12 = $DistanceMatrix->GetValue($AtomIndex1, $AtomIndex2, $SkipIndexCheck);
      if ($Distance12 < $This->{MinDistance} || $Distance12 > $This->{MaxDistance}) {
	next ATOMINDEX2;
      }

      ATOMINDEX3: for $AtomIndex3 (($AtomIndex2 + 1) .. ($NumOfAtoms - 1)) {
	$AtomID3 = $This->{AtomIndexToID}{$AtomIndex3};
	if (!exists($This->{AssignedAtomTypes}{$AtomID3})) {
	  next ATOMINDEX3;
	}
	$AtomType3 = $This->{AssignedAtomTypes}{$AtomID3};

	$Distance13 = $DistanceMatrix->GetValue($AtomIndex1, $AtomIndex3, $SkipIndexCheck);
	$Distance23 = $DistanceMatrix->GetValue($AtomIndex2, $AtomIndex3, $SkipIndexCheck);

	if ($Distance13 < $This->{MinDistance} || $Distance13 > $This->{MaxDistance}) {
	  next ATOMINDEX3;
	}
	if ($Distance23 < $This->{MinDistance} || $Distance23 > $This->{MaxDistance}) {
	  next ATOMINDEX3;
	}
	if ($This->{UseTriangleInequality} && !$This->_DoDistancesSatisfyTriangleInequality($Distance12, $Distance13, $Distance23)) {
	  next ATOMINDEX3;
	}

	$AtomTripletID = $This->_GetAtomTripletID($AtomType1, $Distance23, $AtomType2, $Distance13, $AtomType3, $Distance12);
	if (!exists $This->{AtomTripletsCount}{$AtomTripletID}) {
	  $This->{AtomTripletsCount}{$AtomTripletID} = 0;
	}
	$This->{AtomTripletsCount}{$AtomTripletID} += 1;
      }
    }
  }
  return $This;
}

# Check triangle inequality...
#
sub _DoDistancesSatisfyTriangleInequality {
  my($This, $Distance1, $Distance2, $Distance3) = @_;

  if ( !($Distance1 > abs($Distance2 - $Distance3) && $Distance1 < ($Distance2 + $Distance3)) ) {
    return 0;
  }
  if ( !($Distance2 > abs($Distance1 - $Distance3) && $Distance2 < ($Distance1 + $Distance3)) ) {
    return 0;
  }
  if ( !($Distance3 > abs($Distance1 - $Distance2) && $Distance3 < ($Distance1 + $Distance2)) ) {
    return 0;
  }
  return 1;
}

# Get atom triplet ID corresponding to atom types and distances corresponding to atom triplet...
#
sub _GetAtomTripletID {
  my($This, $ATx, $Dyz, $ATy, $Dxz, $ATz, $Dxy) = @_;
  my($AtomTripletID, @AtomIDs);

  @AtomIDs = ();

  @AtomIDs = sort("${ATx}-D${Dyz}", "${ATy}-D${Dxz}", "${ATz}-D${Dxy}");
  $AtomTripletID = join "-", @AtomIDs;

  return $AtomTripletID;
}

# Set final fingerpritns vector...
#
sub _SetFinalFingerprints {
  my($This) = @_;
  my($AtomTripletID, $Value, @Values);

  # Mark successful generation of fingerprints...
  $This->{FingerprintsGenerated} = 1;

  @Values = ();
  @{$This->{AtomTripletsIDs}} = ();

  for $AtomTripletID (sort keys %{$This->{AtomTripletsCount}}) {
    push @{$This->{AtomTripletsIDs}}, $AtomTripletID;
    $Value = $This->{AtomTripletsCount}{$AtomTripletID};
    push @Values, $Value;
  }

  # Add AtomTripletsIDs and values to fingerprint vector...
  $This->{FingerprintsVector}->AddValueIDs(\@{$This->{AtomTripletsIDs}});
  $This->{FingerprintsVector}->AddValues(\@Values);

  return $This;
}

# Get atom triplet IDs corresponding to atom triplets count values in fingerprint
# vector as an array or reference to an array...
#
# AtomTripletIDs list differes in molecules and is generated during finalization
# of fingerprints to make sure the fingerprint vector containing count values
# matches the atom triplets array.
#
sub GetAtomTripletIDs {
  my($This) = @_;

  return wantarray ? @{$This->{AtomTripletsIDs}} : \@{$This->{AtomTripletsIDs}};
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

# Initialize atomic invariants atom types to use for generating atom IDs in atom triplets...
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
#   ATx = Atomic invariants atom type for atom x
#   ATy = Atomic invariants atom type for atom y
#   ATz = Atomic invariants atom type for atom z
#
#   Dxy = Distance between Px and Py
#   Dxz = Distance between Px and Pz
#   Dyz = Distance between Py and Pz
#
# Then:
#
#   Atom triplet AtomID generated by AtomTypes::AtomicInvariantsAtomTypes class corresponds to:
#
#     AS.X<n>.BO<n>.LBO<n>.<SB><n>.<DB><n>.<TB><n>.H<n>.Ar.RA.FC<+n/-n>.MN<n>.SM<n>
#
#  Toplogical atom triplet ID between atom IDs ATx, ATy and ATz corresponds to:
#
#    ATx-Dyz-ATy-Dxz-ATz-Dxy
#
# Except for AS which is a required atomic invariant in atom triplet AtomIDs, all other atomic invariants are
# optional. Default atomic invariants used for AtomID are: AS, X<n>, BO<n>, H<n>, FC<+n/-n>.
# AtomID specification doesn't include atomic invariants with zero or undefined values.
#
# Examples of atom triplet AtomIDs:
#
#   O.X1.BO1.H1 - Hydroxyl oxygen in carboxylate with attached hydrogen and no explicit charge
#   O.X1.BO1.FC-1 - Hydroxyl ozygen in carboxylate with explicit negative charge
#   O.X1.BO2 - Carbonyl oxygen in carboxylate with double bond to carbon
#   O.X2.BO2 - Hydroxyl ozygen in carboxylate attached to carbonyl carbon and another heavy atom
#
#   C.X2.BO3.H1.Ar - Aromatic carbon
#
sub _InitializeAtomicInvariantsAtomTypesInformation {
  my($This) = @_;

  # Default atomic invariants to use for generating atom triplet atom IDs: AS, X, BO, H, FC
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

# Clear cached molecule data...
#
sub _ClearMoleculeDataCache {
  my($This) = @_;

  @{$This->{Atoms}} = ();

  return $This;
}

# Return a string containg data for TopologicalAtomTripletsFingerprints object...
#
sub StringifyTopologicalAtomTripletsFingerprints {
  my($This) = @_;
  my($FingerprintsString);

  # Type of fingerprint...
  $FingerprintsString = "Fingerprint type: $This->{Type}; AtomIdentifierType: $This->{AtomIdentifierType}";

  # Min and max distance...
  $FingerprintsString .= "; MinDistance:  $This->{MinDistance}; MaxDistance: $This->{MaxDistance}; UseTriangleInequality: " . ($This->{UseTriangleInequality} ? "Yes" : "No");

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

  # Total number of atom triplets...
  $FingerprintsString .= "; NumOfAtomTriplets: " . $This->{FingerprintsVector}->GetNumOfValues();

  # FingerprintsVector...
  $FingerprintsString .= "; FingerprintsVector: < $This->{FingerprintsVector} >";

  return $FingerprintsString;
}

1;

__END__

=head1 NAME

TopologicalAtomTripletsFingerprints

=head1 SYNOPSIS

use Fingerprints::TopologicalAtomTripletsFingerprints;

use Fingerprints::TopologicalAtomTripletsFingerprints qw(:all);

=head1 DESCRIPTION

B<TopologicalAtomTripletsFingerprints>  [ Ref 57, Ref 59, Ref 72 ] class provides the following methods:

new, GenerateFingerprints, GetAtomTripletIDs, GetDescription,
SetAtomIdentifierType, SetAtomicInvariantsToUse, SetFunctionalClassesToUse,
SetMaxDistance, SetMinDistance, StringifyTopologicalAtomTripletsFingerprints

B<TopologicalAtomTripletsFingerprints> is derived from B<Fingerprints> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<TopologicalAtomTripletsFingerprints>, B<Fingerprints> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The current release of MayaChemTools supports generation of B<TopologicalAtomTripletsFingerprints>
corresponding to following B<AtomtomIdentifierTypes>:

    AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
    FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
    SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes

Based on the values specified for B<AtomIdentifierType> along with other specified
parameters such as B<AtomicInvariantsToUse> and B<FunctionalClassesToUse>, initial
atom types are assigned to all non-hydrogen atoms in a molecule. Using the distance
matrix for the molecule and initial atom types assigned to non-hydrogen atoms, all unique atom
triplets within B<MinDistance> and B<MaxDistance> are identified and counted. An atom triplet
identifier is generated for each unique atom triplet; the format of atom triplet identifier is:

    <ATx>-Dyz-<ATy>-Dxz-<ATz>-Dxy

    ATx, ATy, ATz: Atom types assigned to atom x, atom y, and atom z
    Dxy: Distance between atom x and atom y
    Dxz: Distance between atom x and atom z
    Dyz: Distance between atom y and atom z

    where <AT1>-D23 <= <AT2>-D13 <= <AT3>-D12

The atom triplet identifiers for all unique atom triplets corresponding to non-hydrogen atoms constitute
topological atom triplets fingerprints of the molecule.

The current release of MayaChemTools generates the following types of topological atom triplets
fingerprints vector strings:

    FingerprintsVector;TopologicalAtomTriplets:AtomicInvariantsAtomTypes:M
    inDistance1:MaxDistance10;3096;NumericalValues;IDsAndValuesString;C.X1
    .BO1.H3-D1-C.X1.BO1.H3-D1-C.X3.BO3.H1-D2 C.X1.BO1.H3-D1-C.X2.BO2.H2-D1
    0-C.X3.BO4-D9 C.X1.BO1.H3-D1-C.X2.BO2.H2-D3-N.X3.BO3-D4 C.X1.BO1.H3-D1
    -C.X2.BO2.H2-D4-C.X2.BO2.H2-D5 C.X1.BO1.H3-D1-C.X2.BO2.H2-D6-C.X3....;
    1 2 2 2 2 2 2 2 8 8 4 8 4 4 2 2 2 2 4 2 2 2 4 2 2 2 2 1 2 2 4 4 4 2 2
    2 4 4 4 8 4 4 2 4 4 4 2 4 4 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 8...

    FingerprintsVector;TopologicalAtomTriplets:AtomicInvariantsAtomTypes:M
    inDistance1:MaxDistance10;3096;NumericalValues;IDsAndValuesPairsString
    ;C.X1.BO1.H3-D1-C.X1.BO1.H3-D1-C.X3.BO3.H1-D2 1 C.X1.BO1.H3-D1-C.X2.BO
    2.H2-D10-C.X3.BO4-D9 2 C.X1.BO1.H3-D1-C.X2.BO2.H2-D3-N.X3.BO3-D4 2 C.X
    1.BO1.H3-D1-C.X2.BO2.H2-D4-C.X2.BO2.H2-D5 2 C.X1.BO1.H3-D1-C.X2.BO2.H2
    -D6-C.X3.BO3.H1-D5 2 C.X1.BO1.H3-D1-C.X2.BO2.H2-D6-C.X3.BO3.H1-D7 2...

    FingerprintsVector;TopologicalAtomTriplets:DREIDINGAtomTypes:MinDistan
    ce1:MaxDistance10;2377;NumericalValues;IDsAndValuesString;C_2-D1-C_2-D
    9-C_3-D10 C_2-D1-C_2-D9-C_R-D10 C_2-D1-C_3-D1-C_3-D2 C_2-D1-C_3-D10-C_
    3-D9 C_2-D1-C_3-D2-C_3-D3 C_2-D1-C_3-D2-C_R-D3 C_2-D1-C_3-D3-C_3-D4 C_
    2-D1-C_3-D3-N_R-D4 C_2-D1-C_3-D3-O_3-D2 C_2-D1-C_3-D4-C_3-D5 C_2-D...;
    1 1 1 2 1 1 3 1 1 2 2 1 1 1 1 1 1 1 1 2 1 3 4 5 1 1 6 4 2 2 3 1 1 1 2
    2 1 2 1 1 2 2 2 1 2 1 2 1 1 3 3 2 6 4 2 1 1 1 2 2 1 1 1 1 1 1 1 1 1...

    FingerprintsVector;TopologicalAtomTriplets:EStateAtomTypes:MinDistance
    1:MaxDistance10;3298;NumericalValues;IDsAndValuesString;aaCH-D1-aaCH-D
    1-aaCH-D2 aaCH-D1-aaCH-D1-aasC-D2 aaCH-D1-aaCH-D10-aaCH-D9 aaCH-D1-aaC
    H-D10-aasC-D9 aaCH-D1-aaCH-D2-aaCH-D3 aaCH-D1-aaCH-D2-aasC-D1 aaCH-D1-
    aaCH-D2-aasC-D3 aaCH-D1-aaCH-D3-aasC-D2 aaCH-D1-aaCH-D4-aasC-D5 aa...;
    6 4 24 4 16 8 8 4 8 8 8 12 10 14 4 16 24 4 12 2 2 4 1 10 2 2 15 2 2 2
    2 2 2 14 4 2 2 2 2 1 2 10 2 2 4 1 2 4 8 3 3 3 4 6 4 2 2 3 3 1 1 1 2 1
    2 2 4 2 3 2 1 2 4 5 3 2 2 1 2 4 3 2 8 12 6 2 2 4 4 7 1 4 2 4 2 2 2 ...

    FingerprintsVector;TopologicalAtomTriplets:FunctionalClassAtomTypes:Mi
    nDistance1:MaxDistance10;2182;NumericalValues;IDsAndValuesString;Ar-D1
    -Ar-D1-Ar-D2 Ar-D1-Ar-D1-Ar.HBA-D2 Ar-D1-Ar-D10-Ar-D9 Ar-D1-Ar-D10-Hal
    -D9 Ar-D1-Ar-D2-Ar-D2 Ar-D1-Ar-D2-Ar-D3 Ar-D1-Ar-D2-Ar.HBA-D1 Ar-D1-Ar
    -D2-Ar.HBA-D2 Ar-D1-Ar-D2-Ar.HBA-D3 Ar-D1-Ar-D2-HBD-D1 Ar-D1-Ar-D2...;
    27 1 32 2 2 63 3 2 1 2 1 2 3 1 1 40 3 1 2 2 2 2 4 2 2 47 4 2 2 1 2 1 5
    2 2 51 4 3 1 3 1 9 1 1 50 3 3 4 1 9 50 2 2 3 3 5 45 1 1 1 2 1 2 2 3 3
    4 4 3 2 1 1 3 4 5 5 3 1 2 3 2 3 5 7 2 7 3 7 1 1 2 2 2 2 3 1 4 3 1 2...

    FingerprintsVector;TopologicalAtomTriplets:MMFF94AtomTypes:MinDistance
    1:MaxDistance10;2966;NumericalValues;IDsAndValuesString;C5A-D1-C5A-D1-
    N5-D2 C5A-D1-C5A-D2-C5B-D2 C5A-D1-C5A-D3-CB-D2 C5A-D1-C5A-D3-CR-D2 C5A
    -D1-C5B-D1-C5B-D2 C5A-D1-C5B-D2-C=ON-D1 C5A-D1-C5B-D2-CB-D1 C5A-D1-C5B
    -D3-C=ON-D2 C5A-D1-C5B-D3-CB-D2 C5A-D1-C=ON-D3-NC=O-D2 C5A-D1-C=ON-D3-
    O=CN-D2 C5A-D1-C=ON-D4-NC=O-D3 C5A-D1-C=ON-D4-O=CN-D3 C5A-D1-CB-D1-...

    FingerprintsVector;TopologicalAtomTriplets:SLogPAtomTypes:MinDistance1
    :MaxDistance10;3710;NumericalValues;IDsAndValuesString;C1-D1-C1-D1-C11
    -D2 C1-D1-C1-D1-CS-D2 C1-D1-C1-D10-C5-D9 C1-D1-C1-D3-C10-D2 C1-D1-C1-D
    3-C5-D2 C1-D1-C1-D3-CS-D2 C1-D1-C1-D3-CS-D4 C1-D1-C1-D4-C10-D5 C1-D1-C
    1-D4-C11-D5 C1-D1-C1-D5-C10-D4 C1-D1-C1-D5-C5-D4 C1-D1-C1-D6-C11-D7 C1
    -D1-C1-D6-CS-D5 C1-D1-C1-D6-CS-D7 C1-D1-C1-D8-C11-D9 C1-D1-C1-D8-CS...

    FingerprintsVector;TopologicalAtomTriplets:SYBYLAtomTypes:MinDistance1
    :MaxDistance10;2332;NumericalValues;IDsAndValuesString;C.2-D1-C.2-D9-C
    .3-D10 C.2-D1-C.2-D9-C.ar-D10 C.2-D1-C.3-D1-C.3-D2 C.2-D1-C.3-D10-C.3-
    D9 C.2-D1-C.3-D2-C.3-D3 C.2-D1-C.3-D2-C.ar-D3 C.2-D1-C.3-D3-C.3-D4 C.2
    -D1-C.3-D3-N.ar-D4 C.2-D1-C.3-D3-O.3-D2 C.2-D1-C.3-D4-C.3-D5 C.2-D1-C.
    3-D5-C.3-D6 C.2-D1-C.3-D5-O.3-D4 C.2-D1-C.3-D6-C.3-D7 C.2-D1-C.3-D7...

    FingerprintsVector;TopologicalAtomTriplets:TPSAAtomTypes:MinDistance1:
    MaxDistance10;1007;NumericalValues;IDsAndValuesString;N21-D1-N7-D3-Non
    e-D4 N21-D1-N7-D5-None-D4 N21-D1-None-D1-None-D2 N21-D1-None-D2-None-D
    2 N21-D1-None-D2-None-D3 N21-D1-None-D3-None-D4 N21-D1-None-D4-None-D5
     N21-D1-None-D4-O3-D3 N21-D1-None-D4-O4-D3 N21-D1-None-D5-None-D6 N21-
    D1-None-D6-None-D7 N21-D1-None-D6-O4-D5 N21-D1-None-D7-None-D8 N21-...

    FingerprintsVector;TopologicalAtomTriplets:UFFAtomTypes:MinDistance1:M
    axDistance10;2377;NumericalValues;IDsAndValuesString;C_2-D1-C_2-D9-C_3
    -D10 C_2-D1-C_2-D9-C_R-D10 C_2-D1-C_3-D1-C_3-D2 C_2-D1-C_3-D10-C_3-D9 
    C_2-D1-C_3-D2-C_3-D3 C_2-D1-C_3-D2-C_R-D3 C_2-D1-C_3-D3-C_3-D4 C_2-D1-
    C_3-D3-N_R-D4 C_2-D1-C_3-D3-O_3-D2 C_2-D1-C_3-D4-C_3-D5 C_2-D1-C_3-D5-
    C_3-D6 C_2-D1-C_3-D5-O_3-D4 C_2-D1-C_3-D6-C_3-D7 C_2-D1-C_3-D7-C_3-...

=head2 METHODS

=over 4

=item B<new>

    $NewTopologicalAtomTripletsFingerprints = new TopologicalAtomTripletsFingerprints(
                                                   %NamesAndValues);

Using specified I<TopologicalAtomTripletsFingerprints> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<TopologicalAtomTripletsFingerprints>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'TopologicalAtomTriplets'
    MinDistance = 1
    MaxDistance = 10
    UseTriangleInequality = 1
    AtomIdentifierType = ''
    AtomicInvariantsToUse = ['AS', 'X', 'BO', 'H', 'FC']
    FunctionalClassesToUse = ['HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal']

Examples:

    $TopologicalAtomTripletsFingerprints = new TopologicalAtomTripletsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes');

    $TopologicalAtomTripletsFingerprints = new TopologicalAtomTripletsFingerprints(
                              'Molecule' => $Molecule,
                              'MinDistance' => 1,
                              'MaxDistance' => 10,
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes',
                              'AtomicInvariantsToUse' =>
                                              ['AS', 'X', 'BO', 'H', 'FC'] );

    $TopologicalAtomTripletsFingerprints = new TopologicalAtomTripletsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'DREIDINGAtomTypes');

    $TopologicalAtomTripletsFingerprints = new TopologicalAtomTripletsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'MMFF94AtomTypes');

    $TopologicalAtomTripletsFingerprints = new TopologicalAtomTripletsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'TPSAAtomTypes');

    $TopologicalAtomTripletsFingerprints = new TopologicalAtomTripletsFingerprints(
                              'Molecule' => $Molecule,
                              'MinDistance' => 1,
                              'MaxDistance' => 10,
                              'AtomIdentifierType' =>
                                              'FunctionalClassAtomTypes',
                              'FunctionalClassesToUse' =>
                                              ['HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal']);

    $TopologicalAtomTripletsFingerprints->GenerateFingerprints();
    print "$TopologicalAtomTripletsFingerprints\n";

=item B<GetDescription>

    $Return = $TopologicalAtomTripletsFingerprints->GetDescription();

Returns a string containing description of topological atom triplets fingerprints.

=item B<GenerateFingerprints>

    $TopologicalAtomTripletsFingerprints->GenerateFingerprints();

Generates topological atom triplets fingerprints and returns I<TopologicalAtomTripletsFingerprints>.

=item B<GetAtomTripletIDs>

    $AtomTripletIDsRef = $TopologicalAtomTripletsFingerprints->GetAtomTripletIDs();
    @AtomTripletIDs = $TopologicalAtomTripletsFingerprints->GetAtomTripletIDs();

Returns atom triplet IDs corresponding to atom triplets count values in topological atom triplets
fingerprints vector as an array or reference to an array.

=item B<SetAtomIdentifierType>

    $TopologicalAtomTripletsFingerprints->SetAtomIdentifierType($IdentifierType);

Sets atom I<IdentifierType> to use during atom triplets fingerprints generation and
returns I<TopologicalAtomTripletsFingerprints>.

Possible values: I<AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes,
TPSAAtomTypes, UFFAtomTypes>.

=item B<SetAtomicInvariantsToUse>

    $TopologicalAtomTripletsFingerprints->SetAtomicInvariantsToUse($ValuesRef);
    $TopologicalAtomTripletsFingerprints->SetAtomicInvariantsToUse(@Values);

Sets atomic invariants to use during I<AtomicInvariantsAtomTypes> value of I<AtomIdentifierType>
for topological atom triplets fingerprints generation and returns I<TopologicalAtomTripletsFingerprints>.

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

    $TopologicalTripletsFingerprints->SetFunctionalClassesToUse($ValuesRef);
    $TopologicalTripletsFingerprints->SetFunctionalClassesToUse(@Values);

Sets functional classes invariants to use during I<FunctionalClassAtomTypes> value of I<AtomIdentifierType>
for topological atom triplets fingerprints generation and returns I<TopologicalAtomTripletsFingerprints>.

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

    $TopologicalAtomTripletsFingerprints->SetMaxDistance($Distance);

Sets maximum distance to use during topological atom triplets fingerprints generation and
returns I<TopologicalAtomTripletsFingerprints>.

=item B<SetMinDistance>

    $TopologicalAtomTripletsFingerprints->SetMinDistance($Distance);

Sets minimum distance to use during topological atom triplets fingerprints generation and
returns I<TopologicalAtomTripletsFingerprints>.

=item B<StringifyTopologicalAtomTripletsFingerprints>

    $String = $TopologicalAtomTripletsFingerprints->
                  StringifyTopologicalAtomTripletsFingerprints();

Returns a string containing information about I<TopologicalAtomTripletsFingerprints> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Fingerprints.pm, FingerprintsStringUtil.pm, AtomNeighborhoodsFingerprints.pm,
AtomTypesFingerprints.pm, EStateIndiciesFingerprints.pm, ExtendedConnectivityFingerprints.pm,
MACCSKeys.pm, PathLengthFingerprints.pm, TopologicalAtomPairsFingerprints.pm,
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
