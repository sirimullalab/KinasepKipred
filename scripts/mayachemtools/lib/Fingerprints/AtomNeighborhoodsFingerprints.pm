package Fingerprints::AtomNeighborhoodsFingerprints;
#
# File: AtomNeighborhoodsFingerprints.pm
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
use overload '""' => 'StringifyAtomNeighborhoodsFingerprints';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeAtomNeighborhoodsFingerprints();

  $This->_InitializeAtomNeighborhoodsFingerprintsProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeAtomNeighborhoodsFingerprints {
  my($This) = @_;

  # Type of fingerprint...
  $This->{Type} = 'AtomNeighborhoods';

  # Type of vector...
  $This->{VectorType} = 'FingerprintsVector';

  # Type of FingerprintsVector...
  $This->{FingerprintsVectorType} = 'AlphaNumericalValues';

  # Minimum and maximum atomic neighborhoods radii...
  $This->{MinNeighborhoodRadius} = 0;
  $This->{MaxNeighborhoodRadius} = 2;

  # Atom identifier type to use for atom IDs in atom neighborhood atoms...
  #
  # Currently supported values are: AtomicInvariantsAtomTypes, DREIDINGAtomTypes,
  # EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
  # SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes
  #
  $This->{AtomIdentifierType} = '';

  # Atom types assigned to each heavy atom...
  %{$This->{AssignedAtomTypes}} = ();

  # Atom neighorhoods with in specified atom radii..
  %{$This->{AtomNeighborhoods}} = ();

  # Atom neighborhoods atom types count at different neighborhoods...
  %{$This->{NeighborhoodAtomTypesCount}} = ();

  # Atom neighborhood identifiers using specified atom identifier types methodology...
  @{$This->{AtomNeighborhoodsIdentifiers}} = ();
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize object properties....
sub _InitializeAtomNeighborhoodsFingerprintsProperties {
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
  if (exists $NamesAndValues{Size}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated with a user specified size: It's an arbitrary length vector...";
  }
  if (!exists $NamesAndValues{AtomIdentifierType}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying AtomIdentifierType...";
  }

  $This->_InitializeFingerprintsVector();

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
    croak "Error: ${ClassName}->SetAtomIdentifierType: Can't change intial atom identifier type:  It's already set...";
  }

  $This->{AtomIdentifierType} = $IdentifierType;

  # Initialize atom identifier type information...
  $This->_InitializeAtomIdentifierTypeInformation();

  return $This;
}

# Set minimum atom neighborhood radius...
#
sub SetMinNeighborhoodRadius {
  my($This, $Value) = @_;

  if (!TextUtil::IsInteger($Value)) {
    croak "Error: ${ClassName}->SetMinNeighborhoodRadius: MinNeighborhoodRadius value, $Value, is not valid:  It must be an  integer...";
  }

  if ($Value < 0 ) {
    croak "Error: ${ClassName}->SetMinNeighborhoodRadius: MinNeighborhoodRadius value, $Value, is not valid:  It must be >= 0...";
  }
  $This->{MinNeighborhoodRadius} = $Value;

  return $This;
}

# Set maximum atom neighborhood radius...
#
sub SetMaxNeighborhoodRadius {
  my($This, $Value) = @_;

  if (!TextUtil::IsInteger($Value)) {
    croak "Error: ${ClassName}->SetMaxNeighborhoodRadius: MaxNeighborhoodRadius value, $Value, is not valid:  It must be an  integer...";
  }

  if ($Value < 0 ) {
    croak "Error: ${ClassName}->SetMaxNeighborhoodRadius: MaxNeighborhoodRadius value, $Value, is not valid:  It must be >= 0...";
  }
  $This->{MaxNeighborhoodRadius} = $Value;

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

  return "$This->{Type}:$This->{AtomIdentifierType}:MinRadius$This->{MinNeighborhoodRadius}:MaxRadius$This->{MaxNeighborhoodRadius}";
}

# Generate atom neighborhood [ Ref 53-56, Ref 73 ] fingerprints...
#
# Methodology:
#   . Assign atom types to all non-hydrogen atoms in the molecule
#   . Get atom neighborhoods up to MaxNeighborhoodRadis
#   . Count unqiue atom types at each neighborhood radii for all heavy atoms
#   . Generate neighborhood identifiers for all neighborhoods around central
#     heavy atom
#      . Atom neighborhood identifier for a specific radii is generated using neighborhood
#        radius, assigned atom type and its count as follows:
#
#            NR<n>-<AtomType>-ATC<n>
#
#      . Atom neighborhood identifier for a central atom at all specified radii is generated
#        by concatenating neighborhood identifiers at each radii by colon:
#
#            NR<n>-<AtomType>-ATC<n>:NR<n>-<AtomType>-ATC<n>:
#
#   . Set final fingerprints as list of neighborhood atom indentifiers
#
sub GenerateFingerprints {
  my($This) = @_;

  if ($This->{MinNeighborhoodRadius} > $This->{MaxNeighborhoodRadius}) {
    croak "Error: ${ClassName}->GenerateFingerprints: No fingerpritns generated: MinLength, $This->{MinNeighborhoodRadius}, must be less than MaxLength, $This->{MaxNeighborhoodRadius}...";
  }

  # Cache appropriate molecule data...
  $This->_SetupMoleculeDataCache();

  # Assign atom types to all heavy atoms...
  if (!$This->_AssignAtomTypes()) {
    carp "Warning: ${ClassName}->GenerateFingerprints: $This->{AtomIdentifierType} fingerprints generation didn't succeed: Couldn't assign valid $This->{AtomIdentifierType} to all atoms...";
    return $This;
  }

  # Intialize atom neighborhoods information...
  $This->_InitializeAtomNeighborhoods();

  # Identify atom neighborhoods with in specified radii...
  $This->_GetAtomNeighborhoods();

  # Count atom neighborhoods atom types...
  $This->_CountAtomNeighborhoodsAtomTypes();

  # Genenerate atom neighborhood identifiers...
  $This->_GenerateAtomNeighborhoodIdentifiers();

  # Set final fingerprints...
  $This->_SetFinalFingerprints();

  # Clear cached molecule data...
  $This->_ClearMoleculeDataCache();

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
sub _InitializeAtomNeighborhoods {
  my($This) = @_;
  my($Radius);

  # Initialize atom neighborhood count information between specified radii...
  %{$This->{NeighborhoodAtomTypesCount}} = ();

  for $Radius ($This->{MinNeighborhoodRadius} .. $This->{MaxNeighborhoodRadius}) {
    %{$This->{NeighborhoodAtomTypesCount}{$Radius}} = ();
  }

  # Initialize atom neighborhoods atoms information at all specified radii...
  #
  %{$This->{AtomNeighborhoods}} = ();

  for $Radius (0 .. $This->{MaxNeighborhoodRadius}) {
    %{$This->{AtomNeighborhoods}{$Radius}} = ();
  }

  return $This;
}

# Collect atom neighborhoods upto maximum neighborhood radius...
#
# Notes:
#  . Fingerprints are only generated for neighborhoods between specified minimum
#    and maximum neighborhood radii.
#
sub _GetAtomNeighborhoods {
  my($This) = @_;
  my($Atom, $AtomID, $MaxRadius, $Radius, $Molecule);

  $MaxRadius = $This->{MaxNeighborhoodRadius};
  $Molecule = $This->GetMolecule();

  # Collect atom neighborhoods...

  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomID = $Atom->GetID();
    $Radius = 0;

    if ($MaxRadius == 0) {
      # Atom is its own neighborhood at 0 radius...
      my(@AtomNeighborhoodsAtoms);

      @AtomNeighborhoodsAtoms = ($Atom);
      $This->{AtomNeighborhoods}{$Radius}{$AtomID} = \@AtomNeighborhoodsAtoms;

      next ATOM;
    }

    # Collect available atom neighborhoods at different neighborhood radii levels...
    my($AtomNeighborhoodAtomsRef);

    for $AtomNeighborhoodAtomsRef ($Molecule->GetAtomNeighborhoodsWithRadiusUpto($Atom, $MaxRadius)) {
      $This->{AtomNeighborhoods}{$Radius}{$AtomID} = $AtomNeighborhoodAtomsRef;
      $Radius++;
    }
  }
  return $This;
}

# Count atom neighborhoods atom types for each non-hydrogen central atoms with
# neighborhoods in specified radii range...
#
sub _CountAtomNeighborhoodsAtomTypes {
  my($This) = @_;
  my($AtomID, $NeighborhoodAtomID, $Radius, $NeighborhoodAtom, $NeighborhoodAtomType, $AtomNeighborhoodAtomsRef);

  RADIUS: for $Radius (sort { $a <=> $b } keys %{$This->{AtomNeighborhoods}} ) {
    if ($Radius < $This->{MinNeighborhoodRadius} || $Radius > $This->{MaxNeighborhoodRadius}) {
      next RADIUS;
    }
    # Go over the neighborhoods of each atom at the current radius...
    for $AtomID (keys %{$This->{AtomNeighborhoods}{$Radius}}) {
      $AtomNeighborhoodAtomsRef = $This->{AtomNeighborhoods}{$Radius}{$AtomID};
      NEIGHBORHOODATOM: for $NeighborhoodAtom (@{$AtomNeighborhoodAtomsRef}) {
	if ($NeighborhoodAtom->IsHydrogen()) {
	  next NEIGHBORHOODATOM;
	}
	$NeighborhoodAtomID = $NeighborhoodAtom->GetID();
	$NeighborhoodAtomType = $This->{AssignedAtomTypes}{$NeighborhoodAtomID};

	# Count neighbothood atom types for each atom at different radii...
	if (!exists $This->{NeighborhoodAtomTypesCount}{$Radius}{$AtomID}) {
	  %{$This->{NeighborhoodAtomTypesCount}{$Radius}{$AtomID}} = ();
	}
	if (exists $This->{NeighborhoodAtomTypesCount}{$Radius}{$AtomID}{$NeighborhoodAtomType}) {
	  $This->{NeighborhoodAtomTypesCount}{$Radius}{$AtomID}{$NeighborhoodAtomType} += 1;
	}
	else {
	  $This->{NeighborhoodAtomTypesCount}{$Radius}{$AtomID}{$NeighborhoodAtomType} = 1;
	}
      }
    }
  }
  return $This;
}

# Generate atom neighborhood identifiers for each non-hydrogen atom using atom
# neighborhood atom types and their count information...
#
# Let:
#   NR<n> = Neighborhood radius
#   AtomType = Assigned atom type
#   ATC<n> = AtomType count
#
# Then:
#
#   AtomNeighborhoodAtomIdentifier for a neighborhood atom generated for
#   AtomTypes::AtomicInvariantsAtomTypes class corresponds to:
#
#     NR<n>-<AtomType>-ATC<n>
#
#   AtomNeighborhoodsIdentifier for all specified atom neighbothoods of an atom generated for
#   AtomTypes::AtomicInvariantsAtomTypes class corresponds to:
#
#     NR<n>-<AtomType>-ATC<n>;NR<n>-<AtomType>-ATC<n>;...
#
sub _GenerateAtomNeighborhoodIdentifiers {
  my($This) = @_;
  my($Atom, $AtomID, $Radius, $AtomType, $AtomTypeCount, $AtomNeighborhoodIdentifier, @AtomNeighborhoodIdentifiers);

  @{$This->{AtomNeighborhoodsIdentifiers}} = ();

  for $Atom (@{$This->{Atoms}}) {
    $AtomID = $Atom->GetID();
    @AtomNeighborhoodIdentifiers = ();
    RADIUS: for $Radius ($This->{MinNeighborhoodRadius} .. $This->{MaxNeighborhoodRadius}) {
      if (!exists $This->{NeighborhoodAtomTypesCount}{$Radius}{$AtomID}) {
	next RADIUS;
      }
      for $AtomType (sort keys %{$This->{NeighborhoodAtomTypesCount}{$Radius}{$AtomID}}) {
	$AtomTypeCount = $This->{NeighborhoodAtomTypesCount}{$Radius}{$AtomID}{$AtomType};
	push @AtomNeighborhoodIdentifiers, "NR${Radius}-${AtomType}-ATC${AtomTypeCount}";
      }
    }
    $AtomNeighborhoodIdentifier = join(":", @AtomNeighborhoodIdentifiers);
    push @{$This->{AtomNeighborhoodsIdentifiers}}, $AtomNeighborhoodIdentifier;
  }

  return $This;
}

# Set final fingerprits vector...
#
sub _SetFinalFingerprints {
  my($This) = @_;

  # Mark successful generation of fingerprints...
  $This->{FingerprintsGenerated} = 1;

  # Sort AtomNeighborhoodsIdentifiers..
  #
  @{$This->{AtomNeighborhoodsIdentifiers}} = sort @{$This->{AtomNeighborhoodsIdentifiers}};

  # Add sorted atom neighborhood identifiers to FingerprintsVector which is already defined
  # during initialization containing AlphaNumericalValues...
  #
  $This->{FingerprintsVector}->AddValues(\@{$This->{AtomNeighborhoodsIdentifiers}});

  return $This;
}

# Cache  appropriate molecule data...
#
sub _SetupMoleculeDataCache {
  my($This) = @_;

  # Get all non-hydrogen atoms...
  my($NegateAtomCheckMethod);
  $NegateAtomCheckMethod = 1;
  @{$This->{Atoms}} = $This->GetMolecule()->GetAtoms("IsHydrogen", $NegateAtomCheckMethod);

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
# Then:
#
#   Atom type generated by AtomTypes::AtomicInvariantsAtomTypes class corresponds to:
#
#     AS.X<n>.BO<n>.LBO<n>.<SB><n>.<DB><n>.<TB><n>.H<n>.Ar.RA.FC<+n/-n>.MN<n>.SM<n>
#
# Except for AS which is a required atomic invariant in atom types, all other atomic invariants are
# optional. Default atomic invariants used for AtomID are: AS, X<n>, BO<n>, H<n>, FC<+n/-n>.
# AtomID specification doesn't include atomic invariants with zero or undefined values.
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

# Return a string containg data for AtomNeighborhoodsFingerprints object...
#
sub StringifyAtomNeighborhoodsFingerprints {
  my($This) = @_;
  my($FingerprintsString);

  # Type of fingerprint...
  $FingerprintsString = "Fingerprint type: $This->{Type}; AtomIdentifierType: $This->{AtomIdentifierType}; MinNeighborhoodRadius: $This->{MinNeighborhoodRadius}; MaxNeighborhoodRadius: $This->{MaxNeighborhoodRadius}";

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

  # Total number of atom neighborhood atom IDs...
  $FingerprintsString .= "; NumOfAtomNeighborhoodAtomIdentifiers: " . $This->{FingerprintsVector}->GetNumOfValues();

  # FingerprintsVector...
  $FingerprintsString .= "; FingerprintsVector: < $This->{FingerprintsVector} >";

  return $FingerprintsString;
}

1;

__END__

=head1 NAME

AtomNeighborhoodsFingerprints

=head1 SYNOPSIS

use Fingerprints::AtomNeighborhoodsFingerprints;

use Fingerprints::AtomNeighborhoodsFingerprints qw(:all);

=head1 DESCRIPTION

B<AtomNeighborhoodsFingerprints> [ Ref 53-56, Ref 73 ] class provides the following methods:

new, GenerateFingerprints, GetDescription, SetAtomIdentifierType,
SetAtomicInvariantsToUse, SetFunctionalClassesToUse, SetMaxNeighborhoodRadius,
SetMinNeighborhoodRadius, StringifyAtomNeighborhoodsFingerprints

B<AtomNeighborhoodsFingerprints> is derived from B<Fingerprints> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<AtomNeighborhoodsFingerprints>, B<Fingerprints> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The current release of MayaChemTools supports generation of B<AtomNeighborhoodsFingerprints>
corresponding to following B<AtomIdentifierTypes>:

    AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
    FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
    SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes

Based on the values specified for B<AtomIdentifierType> along with other specified
sucb as B<AtomicInvariantsToUse> and B<FunctionalClassesToUse>, initial atom types are
assigned to all non-hydrogen atoms in a molecule. Using atom neighborhoods
around each non-hydrogen central atom corresponding to radii between specified values
B<MinNeighborhoodRadius> and B<MaxNeighborhoodRadius>, unique atom types at each radii
level are counted and an atom neighborhood identifier is generated.

The format of an atom neighborhood identifier around a central non-hydrogen atom at a
specific radius is:

    NR<n>-<AtomType>-ATC<n>

    NR: Neighborhood radius
    AtomType: Assigned atom type
    ATC: Atom type count

The atom neighborhood identifier for non-hydrogen central atom corresponding to all specified radii
is generated by concatenating neighborhood identifiers at each radii by colon as a delimiter:

    NR<n>-<AtomType>-ATC<n>:NR<n>-<AtomType>-ATC<n>:...

The atom neighborhood identifiers for all non-hydrogen central atoms at all specified radii are
concatenated using space as a delimiter and constitute atom neighborhood fingerprint of the molecule.

The current release of MayaChemTools generates the following types of atom neighborhoods
fingerprints vector strings:

    FingerprintsVector;AtomNeighborhoods:AtomicInvariantsAtomTypes:MinRadi
    us0:MaxRadius2;41;AlphaNumericalValues;ValuesString;NR0-C.X1.BO1.H3-AT
    C1:NR1-C.X3.BO3.H1-ATC1:NR2-C.X1.BO1.H3-ATC1:NR2-C.X3.BO4-ATC1 NR0-C.X
    1.BO1.H3-ATC1:NR1-C.X3.BO3.H1-ATC1:NR2-C.X1.BO1.H3-ATC1:NR2-C.X3.BO4-A
    TC1 NR0-C.X2.BO2.H2-ATC1:NR1-C.X2.BO2.H2-ATC1:NR1-C.X3.BO3.H1-ATC1:NR2
    -C.X2.BO2.H2-ATC1:NR2-N.X3.BO3-ATC1:NR2-O.X1.BO1.H1-ATC1 NR0-C.X2.B...

    FingerprintsVector;AtomNeighborhoods:DREIDINGAtomTypes:MinRadius0:MaxR
    adius2;41;AlphaNumericalValues;ValuesString;NR0-C_2-ATC1:NR1-C_3-ATC1:
    NR1-O_2-ATC1:NR1-O_3-ATC1:NR2-C_3-ATC1 NR0-C_2-ATC1:NR1-C_R-ATC1:NR1-N
    _3-ATC1:NR1-O_2-ATC1:NR2-C_R-ATC3 NR0-C_3-ATC1:NR1-C_2-ATC1:NR1-C_3-AT
    C1:NR2-C_3-ATC1:NR2-O_2-ATC1:NR2-O_3-ATC2 NR0-C_3-ATC1:NR1-C_3-ATC1:NR
    1-N_R-ATC1:NR2-C_3-ATC1:NR2-C_R-ATC2 NR0-C_3-ATC1:NR1-C_3-ATC1:NR2-...

    FingerprintsVector;AtomNeighborhoods:EStateAtomTypes:MinRadius0:MaxRad
    ius2;41;AlphaNumericalValues;ValuesString;NR0-aaCH-ATC1:NR1-aaCH-ATC1:
    NR1-aasC-ATC1:NR2-aaCH-ATC1:NR2-aasC-ATC1:NR2-sF-ATC1 NR0-aaCH-ATC1:NR
    1-aaCH-ATC1:NR1-aasC-ATC1:NR2-aaCH-ATC1:NR2-aasC-ATC1:NR2-sF-ATC1 NR0-
    aaCH-ATC1:NR1-aaCH-ATC1:NR1-aasC-ATC1:NR2-aaCH-ATC1:NR2-aasC-ATC2 NR0-
    aaCH-ATC1:NR1-aaCH-ATC1:NR1-aasC-ATC1:NR2-aaCH-ATC1:NR2-aasC-ATC2 N...

    FingerprintsVector;AtomNeighborhoods:FunctionalClassAtomTypes:MinRadiu
    s0:MaxRadius2;41;AlphaNumericalValues;ValuesString;NR0-Ar-ATC1:NR1-Ar-
    ATC1:NR1-Ar.HBA-ATC1:NR1-None-ATC1:NR2-Ar-ATC2:NR2-None-ATC4 NR0-Ar-AT
    C1:NR1-Ar-ATC2:NR1-Ar.HBA-ATC1:NR2-Ar-ATC5:NR2-None-ATC1 NR0-Ar-ATC1:N
    R1-Ar-ATC2:NR1-HBD-ATC1:NR2-Ar-ATC2:NR2-None-ATC1 NR0-Ar-ATC1:NR1-Ar-A
    TC2:NR1-Hal-ATC1:NR2-Ar-ATC2 NR0-Ar-ATC1:NR1-Ar-ATC2:NR1-None-ATC1:...

    FingerprintsVector;AtomNeighborhoods:MMFF94AtomTypes:MinRadius0:MaxRad
    ius2;41;AlphaNumericalValues;ValuesString;NR0-C5A-ATC1:NR1-C5B-ATC1:NR
    1-CB-ATC1:NR1-N5-ATC1:NR2-C5A-ATC1:NR2-C5B-ATC1:NR2-CB-ATC3:NR2-CR-ATC
    1 NR0-C5A-ATC1:NR1-C5B-ATC1:NR1-CR-ATC1:NR1-N5-ATC1:NR2-C5A-ATC1:NR2-C
    5B-ATC1:NR2-C=ON-ATC1:NR2-CR-ATC3 NR0-C5B-ATC1:NR1-C5A-ATC1:NR1-C5B-AT
    C1:NR1-C=ON-ATC1:NR2-C5A-ATC1:NR2-CB-ATC1:NR2-CR-ATC1:NR2-N5-ATC1:N...

    FingerprintsVector;AtomNeighborhoods:SLogPAtomTypes:MinRadius0:MaxRadi
    us2;41;AlphaNumericalValues;ValuesString;NR0-C1-ATC1:NR1-C10-ATC1:NR1-
    CS-ATC1:NR2-C1-ATC1:NR2-N11-ATC1:NR2-O2-ATC1 NR0-C1-ATC1:NR1-C11-ATC1:
    NR2-C1-ATC1:NR2-C21-ATC1 NR0-C1-ATC1:NR1-C11-ATC1:NR2-C1-ATC1:NR2-C21-
    ATC1 NR0-C1-ATC1:NR1-C5-ATC1:NR1-CS-ATC1:NR2-C1-ATC1:NR2-O2-ATC2:NR2-O
    9-ATC1 NR0-C1-ATC1:NR1-CS-ATC2:NR2-C1-ATC2:NR2-O2-ATC2 NR0-C10-ATC1...

    FingerprintsVector;AtomNeighborhoods:SYBYLAtomTypes:MinRadius0:MaxRadi
    us2;41;AlphaNumericalValues;ValuesString;NR0-C.2-ATC1:NR1-C.3-ATC1:NR1
    -O.co2-ATC2:NR2-C.3-ATC1 NR0-C.2-ATC1:NR1-C.ar-ATC1:NR1-N.am-ATC1:NR1-
    O.2-ATC1:NR2-C.ar-ATC3 NR0-C.3-ATC1:NR1-C.2-ATC1:NR1-C.3-ATC1:NR2-C.3-
    ATC1:NR2-O.3-ATC1:NR2-O.co2-ATC2 NR0-C.3-ATC1:NR1-C.3-ATC1:NR1-N.ar-AT
    C1:NR2-C.3-ATC1:NR2-C.ar-ATC2 NR0-C.3-ATC1:NR1-C.3-ATC1:NR2-C.3-ATC...

    FingerprintsVector;AtomNeighborhoods:TPSAAtomTypes:MinRadius0:MaxRadiu
    s2;41;AlphaNumericalValues;ValuesString;NR0-N21-ATC1:NR1-None-ATC3:NR2
    -None-ATC5 NR0-N7-ATC1:NR1-None-ATC2:NR2-None-ATC3:NR2-O3-ATC1 NR0-Non
    e-ATC1:NR1-N21-ATC1:NR1-None-ATC1:NR2-None-ATC3 NR0-None-ATC1:NR1-N21-
    ATC1:NR1-None-ATC2:NR2-None-ATC6 NR0-None-ATC1:NR1-N21-ATC1:NR1-None-A
    TC2:NR2-None-ATC6 NR0-None-ATC1:NR1-N7-ATC1:NR1-None-ATC1:NR1-O3-AT...

    FingerprintsVector;AtomNeighborhoods:UFFAtomTypes:MinRadius0:MaxRadius
    2;41;AlphaNumericalValues;ValuesString;NR0-C_2-ATC1:NR1-C_3-ATC1:NR1-O
    _2-ATC1:NR1-O_3-ATC1:NR2-C_3-ATC1 NR0-C_2-ATC1:NR1-C_R-ATC1:NR1-N_3-AT
    C1:NR1-O_2-ATC1:NR2-C_R-ATC3 NR0-C_3-ATC1:NR1-C_2-ATC1:NR1-C_3-ATC1:NR
    2-C_3-ATC1:NR2-O_2-ATC1:NR2-O_3-ATC2 NR0-C_3-ATC1:NR1-C_3-ATC1:NR1-N_R
    -ATC1:NR2-C_3-ATC1:NR2-C_R-ATC2 NR0-C_3-ATC1:NR1-C_3-ATC1:NR2-C_3-A...

=head2 METHODS

=over 4

=item B<new>

    $NewAtomNeighborhoodsFingerprints = new AtomNeighborhoodsFingerprints(
                                                   %NamesAndValues);

Using specified I<AtomNeighborhoodsFingerprints> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<AtomNeighborhoodsFingerprints>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'AtomNeighborhoods'
    MinNeighborhoodRadius = 0
    MaxNeighborhoodRadius = 2
    AtomIdentifierType = ''
    AtomicInvariantsToUse = ['AS', 'X', 'BO', 'H', 'FC', 'MN']
    FunctionalClassesToUse = ['HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal']

Examples:

    $AtomNeighborhoodsFingerprints = new AtomNeighborhoodsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              "AtomicInvariantsAtomTypes");

    $AtomNeighborhoodsFingerprints = new AtomNeighborhoodsFingerprints(
                              'Molecule' => $Molecule,
                              'MinNeighborhoodRadius' => 0,
                              'MaxNeighborhoodRadius' => 2,
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes',
                              'AtomicInvariantsToUse' =>
                                              ['AS', 'X', 'BO', 'H', 'FC'] );

    $AtomNeighborhoodsFingerprints = new AtomNeighborhoodsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'SYBYLAtomTypes');

    $AtomNeighborhoodsFingerprints = new AtomNeighborhoodsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'MMFF94AtomTypes');

    $AtomNeighborhoodsFingerprints = new AtomNeighborhoodsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes');

    $AtomNeighborhoodsFingerprints = new AtomNeighborhoodsFingerprints(
                              'Molecule' => $Molecule,
                              'MinNeighborhoodRadius' => 0,
                              'MaxNeighborhoodRadius' => 2,
                              'AtomIdentifierType' =>
                                              'FunctionalClassAtomTypes',
                              'FunctionalClassesToUse' =>
                                          ['HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal'] );

    $AtomNeighborhoodsFingerprints->GenerateFingerprints();
    print "$AtomNeighborhoodsFingerprints\n";

=item B<GenerateFingerprints>

    $AtomNeighborhoodsFingerprints->GenerateFingerprints();

Generates atom neighborhood fingerprints and returns I<AtomNeighborhoodsFingerprints>.

=item B<GetDescription>

    $Description = $AtomNeighborhoodsFingerprints->GetDescription();

Returns a string containing description of atom neighborhood fingerprints.

=item B<SetAtomIdentifierType>

    $AtomNeighborhoodsFingerprints->SetAtomIdentifierType($IdentifierType);

Sets atom I<IdentifierType> to use during atom neighborhood fingerprints generation and
returns I<AtomNeighborhoodsFingerprints>.

Possible values: I<AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes,
TPSAAtomTypes, UFFAtomTypes>.

=item B<SetAtomicInvariantsToUse>

    $AtomNeighborhoodsFingerprints->SetAtomicInvariantsToUse($ValuesRef);
    $AtomNeighborhoodsFingerprints->SetAtomicInvariantsToUse(@Values);

Sets atomic invariants to use during I<AtomicInvariantsAtomTypes> value of I<AtomIdentifierType>
for atom neighborhood fingerprints generation and returns I<AtomNeighborhoodsFingerprints>.

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

    $AtomNeighborhoodsFingerprints->SetFunctionalClassesToUse($ValuesRef);
    $AtomNeighborhoodsFingerprints->SetFunctionalClassesToUse(@Values);

Sets functional classes invariants to use during I<FunctionalClassAtomTypes> value of I<AtomIdentifierType>
for atom neighborhoods fingerprints generation and returns I<AtomNeighborhoodsFingerprints>.

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

=item B<SetMaxNeighborhoodRadius>

    $AtomNeighborhoodsFingerprints->SetMaxNeighborhoodRadius($Radius);

Sets maximum neighborhood radius to use during atom neighborhood fingerprints generation and
returns I<AtomNeighborhoodsFingerprints>.

=item B<SetMinNeighborhoodRadius>

    $AtomNeighborhoodsFingerprints->SetMinNeighborhoodRadius($Radius);

Sets minimum neighborhood radius to use during atom neighborhood fingerprints generation and
returns I<AtomNeighborhoodsFingerprints>.

=item B<StringifyAtomNeighborhoodsFingerprints>

    $String = $Fingerprints->StringifyAtomNeighborhoodsFingerprints();

Returns a string containing information about I<AtomNeighborhoodsFingerprints> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Fingerprints.pm, FingerprintsStringUtil.pm, AtomTypesFingerprints.pm, EStateIndiciesFingerprints.pm,
ExtendedConnectivityFingerprints.pm, MACCSKeys.pm, PathLengthFingerprints.pm,
TopologicalAtomPairsFingerprints.pm, TopologicalAtomTripletsFingerprints.pm,
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
