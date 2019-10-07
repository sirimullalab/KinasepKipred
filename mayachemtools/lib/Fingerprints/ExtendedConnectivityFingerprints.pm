package Fingerprints::ExtendedConnectivityFingerprints;
#
# File: ExtendedConnectivityFingerprints.pm
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
use TextUtil ();
use MathUtil ();
use Fingerprints::Fingerprints;
use Molecule;
use AtomTypes::AtomicInvariantsAtomTypes;
use AtomTypes::FunctionalClassAtomTypes;
use AtomTypes::DREIDINGAtomTypes;
use AtomTypes::EStateAtomTypes;
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
use overload '""' => 'StringifyExtendedConnectivityFingerprints';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeExtendedConnectivityFingerprints();

  $This->_InitializeExtendedConnectivityFingerprintsProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeExtendedConnectivityFingerprints {
  my($This) = @_;

  # Type of fingerprint to generate:
  #
  # ExtendedConnectivity - Set of integer identifiers corresponding to structurally unique features
  # ExtendedConnectivityCount - Set of integer identifiers corresponding to structurally unique features and their count
  # ExtendedConnectivityBits - A bit vector indicating presence/absence of structurally unique features
  #
  $This->{Type} = 'ExtendedConnectivity';

  # Atomic neighborhoods radius for extended connectivity...
  $This->{NeighborhoodRadius} = 2;

  # Size of bit bector to use during generation of ExtendedConnectivityBits fingerprints...
  $This->{Size} = 1024;

  # Min and max size of bit bector to use during generation of ExtendedConnectivityBits fingerprints...
  $This->{MinSize} = 32;
  $This->{MaxSize} = 2**32;

  # Type of atom attributes to use for initial identifier assignment to non-hydrogen atoms
  # during the calculation of extended connectivity fingerprints [ Ref 48, Ref 52 ]...
  #
  # Currently supported values are: AtomicInvariantsAtomTypes, FunctionalClassAtomTypes,
  # DREIDINGAtomTypes, EStateAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes,
  # TPSAAtomTypes, UFFAtomTypes
  #
  $This->{AtomIdentifierType} = '';

  # Random number generator to use during generation of fingerprints bit-vector
  # string: Perl CORE::rand or MayaChemTools MathUtil::random function.
  #
  # The random number generator implemented in MayaChemTools is a variant of
  # linear congruential generator (LCG) as described by Miller et al. [ Ref 120 ].
  # It is also referred to as Lehmer random number generator or Park-Miller
  # random number generator.
  #
  # Unlike Perl's core random number generator function rand, the random number
  # generator implemented in MayaChemTools, MathUtil::random,  generates consistent
  # random values across different platformsfor a specific random seed and leads
  # to generation of portable fingerprints bit-vector strings.
  #
  $This->{UsePerlCoreRandom} = 1;

  # Atom neighorhoods up to specified neighborhood radius...
  %{$This->{AtomNeighborhoods}} = ();

  # Atom identifiers at different neighborhoods up to specified neighborhood radius...
  %{$This->{AtomIdentifiers}} = ();

  # Structurally unique atom identifiers at different neighborhoods up to specified neighborhood radius...
  %{$This->{UniqueAtomIdentifiers}} = ();
  %{$This->{UniqueAtomIdentifiersCount}} = ();

  # Unique atom identifiers at different neighborhoods up to specified neighborhood radius...
  %{$This->{StructurallyUniqueAtomIdentifiers}} = ();
  %{$This->{StructurallyUniqueAtomIdentifiersCount}} = ();

  # Structure feature  information at different neighborhoods up to specified neighborhood
  # radius used during removal of atom indentifiers which are structually equivalent...
  %{$This->{StructureFeatures}} = ();
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize object properties....
sub _InitializeExtendedConnectivityFingerprintsProperties {
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

  # Make sure AtomIdentifierType was specified...
  if (!exists $NamesAndValues{AtomIdentifierType}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying AtomIdentifierType...";
  }

  # Make sure it's power of 2...
  if (exists $NamesAndValues{Size}) {
    if (!TextUtil::IsNumberPowerOfNumber($NamesAndValues{Size}, 2)) {
      croak "Error: ${ClassName}->New: Specified size value, $NamesAndValues{Size}, must be power of 2...";
    }
  }

  if ($This->{Type} =~ /^ExtendedConnectivity$/i) {
    $This->_InitializeExtendedConnectivityFingerprintsVector();
  }
  elsif ($This->{Type} =~ /^ExtendedConnectivityCount$/i) {
    $This->_InitializeExtendedConnectivityCountFingerprintsVector();
  }
  elsif ($This->{Type} =~ /^ExtendedConnectivityBits$/i) {
    $This->_InitializeExtendedConnectivityBitsFingerprintsBitVector();
  }
  else {
    croak "Error: ${ClassName}->_InitializeExtendedConnectivityFingerprintsProperties: Unknown ExtendedConnectivity fingerprints type: $This->{Type}; Supported fingerprints types: ExtendedConnectivity, ExtendedConnectivityCount or ExtendedConnectivityBits...";
  }

  return $This;
}

# Initialize extended connectivity fingerprints vector...
#
sub _InitializeExtendedConnectivityFingerprintsVector {
  my($This) = @_;

  # Type of vector...
  $This->{VectorType} = 'FingerprintsVector';

  # Type of FingerprintsVector...
  $This->{FingerprintsVectorType} = 'AlphaNumericalValues';

  $This->_InitializeFingerprintsVector();

  return $This;
}

# Initialize extended connectivity count fingerprints vector...
#
sub _InitializeExtendedConnectivityCountFingerprintsVector {
  my($This) = @_;

  # Type of vector...
  $This->{VectorType} = 'FingerprintsVector';

  # Type of FingerprintsVector...
  $This->{FingerprintsVectorType} = 'NumericalValues';

  $This->_InitializeFingerprintsVector();

  return $This;
}

# Initialize extended connectivity bit fingerprints vector...
#
sub _InitializeExtendedConnectivityBitsFingerprintsBitVector {
  my($This) = @_;

  # Type of vector...
  $This->{VectorType} = 'FingerprintsBitVector';

  $This->_InitializeFingerprintsBitVector();

  return $This;
}

# Set type...
#
sub SetType {
  my($This, $Type) = @_;

  if ($Type =~ /^ExtendedConnectivity$/i) {
    $This->{Type} = 'ExtendedConnectivity';;
  }
  elsif ($Type =~ /^ExtendedConnectivityCount$/i) {
    $This->{Type} = 'ExtendedConnectivityCount';;
  }
  elsif ($Type =~ /^ExtendedConnectivityBits$/i) {
    $This->{Type} = 'ExtendedConnectivityBits';;
  }
  else {
    croak "Error: ${ClassName}->SetType: Unknown ExtendedConnectivity fingerprints type: $This->{Type}; Supported fingerprints types: ExtendedConnectivity, ExtendedConnectivityCount or ExtendedConnectivityBits...";
  }
  return $This;
}

# Disable vector type change...
#
sub SetVectorType {
  my($This, $Type) = @_;

  croak "Error: ${ClassName}->SetVectorType: Can't change vector type...";

  return $This;
}

# Disable vector type change...
#
sub SetFingerprintsVectorType {
  my($This, $Type) = @_;

  croak "Error: ${ClassName}->SetFingerprintsVectorType: Can't change fingerprints vector type...";

  return $This;
}

# Set intial atom identifier type..
#
sub SetAtomIdentifierType {
  my($This, $IdentifierType) = @_;

  if ($IdentifierType !~ /^(AtomicInvariantsAtomTypes|FunctionalClassAtomTypes|DREIDINGAtomTypes|EStateAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
    croak "Error: ${ClassName}->SetAtomIdentifierType: Specified value, $IdentifierType, for AtomIdentifierType is not vaild. Supported types in current release of MayaChemTools: AtomicInvariantsAtomTypes, FunctionalClassAtomTypes, DREIDINGAtomTypes, EStateAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes and UFFAtomTypes.";
  }

  if ($This->{AtomIdentifierType}) {
    croak "Error: ${ClassName}->SetAtomIdentifierType: Can't change intial atom identifier type:  It's already set...";
  }

  $This->{AtomIdentifierType} = $IdentifierType;

  # Initialize identifier type information...
  $This->_InitializeAtomIdentifierTypeInformation();

  return $This;
}

# Set atom neighborhood radius...
#
sub SetNeighborhoodRadius {
  my($This, $Value) = @_;

  if (!TextUtil::IsInteger($Value)) {
    croak "Error: ${ClassName}->SetNeighborhoodRadius: NeighborhoodRadius value, $Value, is not valid:  It must be an  integer...";
  }

  if ($Value < 0 ) {
    croak "Error: ${ClassName}->SetNeighborhoodRadius: NeighborhoodRadius value, $Value, is not valid:  It must be >= 0...";
  }
  $This->{NeighborhoodRadius} = $Value;

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

  return "$This->{Type}:$This->{AtomIdentifierType}:Radius$This->{NeighborhoodRadius}";
}

# Generate fingerprints...
#
# Methodology:
#   . Assign initial atom identfiers to all non-hydrogen atoms in the molecule
#
#   . Remove duplicates from the initial identifiers and add them to list corresponding
#     to molecule fingerprint
#
#   . For NeighborhoodRadius value of 0, just return the molecule fingerprint list
#
#   . For each NeighborhoodRadius level
#      . For each non-hydrogen CentralAtom at this NeighborhoodRadius level
#         . For each non-hydrogen SuccessorNeighborAtom
#           . Collect (BondOrder AtomIdentifier) pair of values corresponding to
#             (CentralAtom SuccessorNeighborAtom)  and add it to a list
#
#         . Sort list containing (BondOrder AtomIdentifier) pairs first by BondOrder followed
#            by AtomIdendifiers to make these values graph invariant
#         . Generate a hash code for the values in the list
#         . Assign hash code as new atom identifier at the current NeighborhoodRadius level
#         . Save all atoms and bonds corresponding to the substructure involved in
#           generating the hash code to be used for identifying structural duplicate hash code
#
#         . Add the new identifier to the molecule fingerprint list making sure it's not a duplicate
#           identifier
#
#   Hash code atom identifier deduplication:
#     . Track/remove the identifier generated at higher neighborhood radius level
#
#  Structural atom identifier deduplication:
#    . For equivalent atoms and bonds corresponding to substructure at a NeighborhoodRadius level,
#      track/remove the atom identifier with largest value
#
#
sub GenerateFingerprints {
  my($This) = @_;

  # Cache appropriate molecule data...
  $This->_SetupMoleculeDataCache();

  # Assign intial atom identifers...
  if (!$This->_AssignInitialAtomIdentifiers()) {
    carp "Warning: ${ClassName}->GenerateFingerprints: $This->{AtomIdentifierType} fingerprints generation didn't succeed: Couldn't assign valid $This->{AtomIdentifierType} to all atoms...";
    return $This;
  }

  # Identify atom neighborhoods up to specified radius...
  $This->_GetAtomNeighborhoods();

  # Assign atom identifiers to central atoms considering atom neighborhoods at each
  # radius level...
  $This->_AssignAtomIdentifiersToAtomNeighborhoods();

  # Remove duplicates identifiers...
  $This->_RemoveDuplicateAtomIdentifiers();

  # Set final fingerprints...
  $This->_SetFinalFingerprints();

  # Clear cached molecule data...
  $This->_ClearMoleculeDataCache();

  return $This;
}

# Assign appropriate initial atom identifiers...
#
#   Generation of initial identifier for a specific atom involves:
#     . Values of the specified atom attributes are appended in a specific order to
#       generate an initial atom identifier string
#     . A 32 bit unsigned integer hash key, using TextUtil::HashCode function,  is
#       generated for the atom indentifier and assigned to the atom as initial
#       atom identifier.
#
sub _AssignInitialAtomIdentifiers {
  my($This) = @_;
  my($Atom, $AtomID, $Radius, $SpecifiedAtomTypes, $IgnoreHydrogens, $AtomType, $InitialAtomTypeString, $InitialAtomIdentifier);

  # Initialize atom identifiers...
  $This->_InitializeAtomIdentifiers();

  # Set up atom types...
  $IgnoreHydrogens = 1;
  $SpecifiedAtomTypes = undef;

  IDENTIFIERTYPE: {
    if ($This->{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::AtomicInvariantsAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $IgnoreHydrogens, 'AtomicInvariantsToUse' => $This->{AtomicInvariantsToUse});
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::FunctionalClassAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $IgnoreHydrogens, 'FunctionalClassesToUse' => $This->{FunctionalClassesToUse});
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

    croak "Error: ${ClassName}->_AssignInitialAtomIdentifiers: Couldn't assign intial atom identifiers: InitialAtomIdentifierType $This->{AtomIdentifierType} is not supported...";
  }

  # Assign atom types...
  $SpecifiedAtomTypes->AssignAtomTypes();

  # Make sure atom types assignment is successful...
  if (!$SpecifiedAtomTypes->IsAtomTypesAssignmentSuccessful()) {
    return undef;
  }

  # Assign atom identifiers at radius 0...
  $Radius = 0;
  for $Atom (@{$This->{Atoms}}) {
    $AtomID = $Atom->GetID();

    $AtomType = $SpecifiedAtomTypes->GetAtomType($Atom);
    $InitialAtomTypeString = $AtomType ? $AtomType : 'None';

    $InitialAtomIdentifier = TextUtil::HashCode($InitialAtomTypeString);
    $This->{AtomIdentifiers}{$Radius}{$AtomID} = $InitialAtomIdentifier;
  }

  return $This;
}

# Initialize atom identifiers...
#
sub _InitializeAtomIdentifiers {
  my($This) = @_;
  my($Radius, $CurrentRadius);

  $Radius = $This->{NeighborhoodRadius};

  %{$This->{AtomIdentifiers}} = ();
  for $CurrentRadius (0 .. $Radius) {
    # Atom idenfiers key and value correspond to AtomID and AtomIdentifier
    %{$This->{AtomIdentifiers}{$CurrentRadius}} = ();

    # Unique and strcuturally unique idenfiers key and value correspond to AtomIdentifier and AtomID
    %{$This->{UniqueAtomIdentifiers}{$CurrentRadius}} = ();
    %{$This->{UniqueAtomIdentifiersCount}{$CurrentRadius}} = ();

    %{$This->{StructurallyUniqueAtomIdentifiers}{$CurrentRadius}} = ();
    %{$This->{StructurallyUniqueAtomIdentifiersCount}{$CurrentRadius}} = ();
  }

}

# Collect atom neighborhoods upto specified neighborhood radius...
#
sub _GetAtomNeighborhoods {
  my($This) = @_;
  my($Atom, $AtomID, $Radius, $CurrentRadius, $Molecule);

  %{$This->{AtomNeighborhoods}} = ();

  $Radius = $This->{NeighborhoodRadius};
  if ($Radius < 1) {
    # At radius level 0, it's just the atoms...
    return;
  }

  # Initialize neighborhood at different radii...
  for $CurrentRadius (0 .. $Radius) {
    %{$This->{AtomNeighborhoods}{$CurrentRadius}} = ();
  }

  $Molecule = $This->GetMolecule();

  # Collect available atom neighborhoods at different at different neighborhood level for each atom...
  my($AtomsNeighborhoodWithSuccessorAtomsRef);

  for $Atom (@{$This->{Atoms}}) {
    $AtomID = $Atom->GetID();
    $CurrentRadius = 0;
    for $AtomsNeighborhoodWithSuccessorAtomsRef ($Molecule->GetAtomNeighborhoodsWithSuccessorAtomsAndRadiusUpto($Atom, $Radius)) {
      $This->{AtomNeighborhoods}{$CurrentRadius}{$AtomID} = $AtomsNeighborhoodWithSuccessorAtomsRef;
      $CurrentRadius++;
    }
  }
  return $This;
}

# Assign atom identifiers to central atom at each neighborhood radius level...
#
sub _AssignAtomIdentifiersToAtomNeighborhoods {
  my($This) = @_;
  my($Radius, $NextRadius, $Atom, $AtomID, $NeighborhoodAtom, $SuccessorAtom, $SuccessorAtomID, $NeighborhoodAtomSuccessorAtomsRef, $NeighborhoodAtomsWithSuccessorAtomsRef, $Bond, $BondOrder, $SuccessorAtomCount);

  if ($This->{NeighborhoodRadius} < 1) {
    return;
  }

  # Go over the atom neighborhoods at each radius upto specified radius and assign atom
  # indentifiers using their connected successor atoms and their identifiers.
  #
  # For a neighborhood atom at a specified radius, the successor connected atoms correpond
  # to next radius level and the last set of neighorhood atoms don't have any successor connected
  # atoms. Additionally, radius level 0 just correspond to initial atom identifiers.
  #
  # So in order to process atom neighborhood upto specified radius level, the last atom neighborhood
  # doesn't need to be processed: it gets processed at previous radius level as successor connected
  # atoms.
  #
  RADIUS: for $Radius (0 .. ($This->{NeighborhoodRadius} - 1)) {
    ATOM: for $Atom (@{$This->{Atoms}}) {
      $AtomID = $Atom->GetID();

      # Are there any available atom neighborhoods at this radius?
      if (!exists $This->{AtomNeighborhoods}{$Radius}{$AtomID}) {
	next ATOM;
      }
      $NextRadius = $Radius + 1;

      # Go over neighborhood atoms and their successor connected atoms at this radius and collect
      # (BondOrder AtomIdentifier) values for bonded atom pairs. Additionally, keep track of atom and bonds
      # for the neighorhoods to remove identifieres generated from structurally duplicate features.
      #
      my(%BondOrdersAndAtomIdentifiers);

      %BondOrdersAndAtomIdentifiers = ();
      $SuccessorAtomCount = 0;

      NEIGHBORHOODS: for $NeighborhoodAtomsWithSuccessorAtomsRef (@{$This->{AtomNeighborhoods}{$Radius}{$AtomID}}) {
	($NeighborhoodAtom, $NeighborhoodAtomSuccessorAtomsRef) = @{$NeighborhoodAtomsWithSuccessorAtomsRef};

	# Any connected successors for the NeighborhoodAtom?
	if (!@{$NeighborhoodAtomSuccessorAtomsRef}) {
	  next NEIGHBORHOODS;
	}
	SUCCESSORATOM: for $SuccessorAtom (@{$NeighborhoodAtomSuccessorAtomsRef}) {
	  if ($SuccessorAtom->IsHydrogen()) {
	    # Skip successor hydrogen atom...
	    next SUCCESSORATOM;
	  }
	  $SuccessorAtomID = $SuccessorAtom->GetID();
	  $SuccessorAtomCount++;

	  $Bond = $NeighborhoodAtom->GetBondToAtom($SuccessorAtom);
	  $BondOrder = $Bond->IsAromatic() ? "1.5" : $Bond->GetBondOrder();

	  if (!exists $BondOrdersAndAtomIdentifiers{$BondOrder}) {
	    @{$BondOrdersAndAtomIdentifiers{$BondOrder}} = ();
	  }
	  push @{$BondOrdersAndAtomIdentifiers{$BondOrder}}, $This->{AtomIdentifiers}{$Radius}{$SuccessorAtomID};
	}
      }
      if (!$SuccessorAtomCount) {
	next ATOM;
      }
      # Assign a new atom identifier at the NextRadius level...
      $This->_AssignAtomIdentifierToAtomNeighborhood($AtomID, $Radius, \%BondOrdersAndAtomIdentifiers);
    }
 }
  return $This;
}

# Generate and assign atom indentifier for AtomID using atom neighborhood at next radius level...
#
sub _AssignAtomIdentifierToAtomNeighborhood {
  my($This, $AtomID, $Radius, $BondOrdersAndAtomIdentifiersRef) = @_;
  my($NextRadius, $AtomIdentifier,  $SuccessorAtomIdentifier, $BondOrder, $AtomIdentifierString, @AtomIndentifiersInfo);

  $NextRadius = $Radius + 1;

  @AtomIndentifiersInfo = ();

  $AtomIdentifier = $This->{AtomIdentifiers}{$Radius}{$AtomID};
  push @AtomIndentifiersInfo, ($NextRadius, $AtomIdentifier);

  # Sort out successor atom bond order and identifier pairs by bond order followed by atom identifiers
  # in order to make the final atom identifier graph invariant...
  #
  for $BondOrder (sort { $a <=> $b } keys %{$BondOrdersAndAtomIdentifiersRef}) {
    for $SuccessorAtomIdentifier (sort { $a <=> $b } @{$BondOrdersAndAtomIdentifiersRef->{$BondOrder}}) {
      push @AtomIndentifiersInfo, ($BondOrder, $SuccessorAtomIdentifier);
    }
  }
  $AtomIdentifierString = join("", @AtomIndentifiersInfo);
  $AtomIdentifier = TextUtil::HashCode($AtomIdentifierString);

  # Assign atom identifier to the atom at next radius level...
  $This->{AtomIdentifiers}{$NextRadius}{$AtomID} = $AtomIdentifier;

  return $This;
}

# Remove duplicates atom identifiers...
#
sub _RemoveDuplicateAtomIdentifiers {
  my($This) = @_;

  $This->_RemoveDuplicateIdentifiersByValue();
  $This->_RemoveStructurallyDuplicateIdenfiers();

  return $This;
}

# Remove duplicate identifiers at each radius level by just using their value...
#
sub _RemoveDuplicateIdentifiersByValue {
  my($This) = @_;
  my($Radius, $Atom, $AtomID, $AtomIdentifier);

  for $Radius (0 .. $This->{NeighborhoodRadius}) {
    ATOM: for $Atom (@{$This->{Atoms}}) {
      $AtomID = $Atom->GetID();
      if (!exists $This->{AtomIdentifiers}{$Radius}{$AtomID}) {
	next ATOM;
      }
      $AtomIdentifier = $This->{AtomIdentifiers}{$Radius}{$AtomID};
      if (exists $This->{UniqueAtomIdentifiers}{$Radius}{$AtomIdentifier}) {
	# It's a duplicate atom idenfier at this radius level...
	$This->{UniqueAtomIdentifiersCount}{$Radius}{$AtomIdentifier} += 1;
	next ATOM;
      }
      $This->{UniqueAtomIdentifiers}{$Radius}{$AtomIdentifier} = $AtomID;
      $This->{UniqueAtomIdentifiersCount}{$Radius}{$AtomIdentifier} = 1;
    }
  }
  return $This;
}

# Remove structurally duplicate identifiers at each radius level...
#
# Methodology:
#   . For unquie atom identifiers at each radius level, assign complete structure features
#     in terms all the bonds involved to generate that identifier
#   . Use the complete structure features to remover atom identifiers which are
#     structurally equivalent which can also be at earlier radii levels
#
#
sub _RemoveStructurallyDuplicateIdenfiers {
  my($This) = @_;
  my($Radius, $AtomID, $AtomIdentifier, $SimilarAtomIdentifierRadius, $SimilarAtomIdentifier);

  # Setup structure features...
  $This->_SetupStructureFeaturesForAtomIDsInvolvedInUniqueIdentifiers();

  # Identify structurally unqiue identifiers...
  for $Radius (0 .. $This->{NeighborhoodRadius}) {
    ATOMIDENTIFIER: for $AtomIdentifier (sort { $a <=> $b } keys %{$This->{UniqueAtomIdentifiers}{$Radius}}) {
      $AtomID = $This->{UniqueAtomIdentifiers}{$Radius}{$AtomIdentifier};

      ($SimilarAtomIdentifierRadius, $SimilarAtomIdentifier) = $This->_FindStructurallySimilarAtomIdentifier($Radius, $AtomID, $AtomIdentifier);
      if ($SimilarAtomIdentifier) {
	# Current atom identifier is similar to an earlier structurally unique atom identifier...
	$This->{StructurallyUniqueAtomIdentifiersCount}{$SimilarAtomIdentifierRadius}{$SimilarAtomIdentifier} += $This->{UniqueAtomIdentifiersCount}{$Radius}{$AtomIdentifier};
	next ATOMIDENTIFIER;
      }
      $This->{StructurallyUniqueAtomIdentifiers}{$Radius}{$AtomIdentifier} = $AtomID;

      # Set structurally unique atom identifier count to the unique atom identifiers count...
      $This->{StructurallyUniqueAtomIdentifiersCount}{$Radius}{$AtomIdentifier} = $This->{UniqueAtomIdentifiersCount}{$Radius}{$AtomIdentifier};
    }
  }
  return $This;
}

# Set final fingerpritns vector...
#
sub _SetFinalFingerprints {
  my($This) = @_;

  # Mark successful generation of fingerprints...
  $This->{FingerprintsGenerated} = 1;

  if ($This->{Type} =~ /^ExtendedConnectivity$/i) {
    $This->_SetFinalExtendedConnectivityFingerprints();
  }
  elsif ($This->{Type} =~ /^ExtendedConnectivityCount$/i) {
    $This->_SetFinalExtendedConnectivityCountFingerprints();
  }
  elsif ($This->{Type} =~ /^ExtendedConnectivityBits$/i) {
    $This->_SetFinalExtendedConnectivityBitsFingerprints();
  }

  return $This;
}

# Set final extended connectivity fingerpritns vector...
#
sub _SetFinalExtendedConnectivityFingerprints {
  my($This) = @_;
  my($Radius, $AtomIdentifier, @AtomIdentifiers);

  @AtomIdentifiers = ();

  for $Radius (0 .. $This->{NeighborhoodRadius}) {
    for $AtomIdentifier (sort { $a <=> $b } keys %{$This->{StructurallyUniqueAtomIdentifiers}{$Radius}}) {
      push @AtomIdentifiers, $AtomIdentifier;
    }
  }
  # Add atom identifiers to fingerprint vector...
  $This->{FingerprintsVector}->AddValues(\@AtomIdentifiers);

  return $This;
}

# Set final extended connectivity count fingerpritns vector...
#
sub _SetFinalExtendedConnectivityCountFingerprints {
  my($This) = @_;
  my($Radius, $AtomIdentifier, $AtomIdentifierCount, @AtomIdentifiers, @AtomIdentifiersCount);

  @AtomIdentifiers = (); @AtomIdentifiersCount = ();

  for $Radius (0 .. $This->{NeighborhoodRadius}) {
    for $AtomIdentifier (sort { $a <=> $b } keys %{$This->{StructurallyUniqueAtomIdentifiers}{$Radius}}) {
      $AtomIdentifierCount = $This->{StructurallyUniqueAtomIdentifiersCount}{$Radius}{$AtomIdentifier};
      push @AtomIdentifiers, $AtomIdentifier;
      push @AtomIdentifiersCount, $AtomIdentifierCount;
    }
  }
  # Add atom identifiers to fingerprint vector as value IDs...
  $This->{FingerprintsVector}->AddValueIDs(\@AtomIdentifiers);

  # Add atom identifiers to count to fingerprint vector as values...
  $This->{FingerprintsVector}->AddValues(\@AtomIdentifiersCount);

  return $This;
}

# Set final extended connectivity bits fingerpritns vector...
#
sub _SetFinalExtendedConnectivityBitsFingerprints {
  my($This) = @_;
  my($Radius, $AtomIdentifier, $FingerprintsBitVector, $Size, $SkipBitPosCheck, $AtomIdentifierBitPos, $SetBitNum);

  $FingerprintsBitVector = $This->{FingerprintsBitVector};

  $Size = $This->{Size};

  $SkipBitPosCheck = 1;

  for $Radius (0 .. $This->{NeighborhoodRadius}) {
    for $AtomIdentifier (keys %{$This->{StructurallyUniqueAtomIdentifiers}{$Radius}}) {
      # Set random number seed...
      if ($This->{UsePerlCoreRandom}) {
	CORE::srand($AtomIdentifier);
      }
      else {
	MathUtil::srandom($AtomIdentifier);
      }

      # Set bit position...
      $AtomIdentifierBitPos = $This->{UsePerlCoreRandom} ? int(CORE::rand($Size)) : int(MathUtil::random($Size));
      $FingerprintsBitVector->SetBit($AtomIdentifierBitPos, $SkipBitPosCheck);
    }
  }
  return $This;
}


# Identify structurally unique identifiers by comparing structure features involved in
# generating identifiear by comparing it agains all the previous structurally unique
# identifiers...
#
sub _FindStructurallySimilarAtomIdentifier {
  my($This, $SpecifiedRadius, $SpecifiedAtomID, $SpecifiedAtomIdentifier) = @_;
  my($Radius, $AtomID, $AtomIdentifier, $FeatureAtomCount, $FeatureAtomIDsRef,  $SpecifiedFeatureAtomID, $SpecifiedFeatureAtomCount, $SpecifiedFeatureAtomIDsRef);

  if ($SpecifiedRadius == 0) {
    # After duplicate removal by value, all identifier at radius level 0 would be structurally unique...
    return (undef, undef);
  }

  $SpecifiedFeatureAtomCount = $This->{StructureFeatures}{AtomCount}{$SpecifiedRadius}{$SpecifiedAtomID};
  $SpecifiedFeatureAtomIDsRef = $This->{StructureFeatures}{AtomIDs}{$SpecifiedRadius}{$SpecifiedAtomID};

  # No need to compare features at radius 0...
  for $Radius (1 .. $SpecifiedRadius) {
    ATOMIDENTIFIER: for $AtomIdentifier (keys %{$This->{StructurallyUniqueAtomIdentifiers}{$Radius}}) {
      $AtomID = $This->{StructurallyUniqueAtomIdentifiers}{$Radius}{$AtomIdentifier};

      $FeatureAtomCount = $This->{StructureFeatures}{AtomCount}{$Radius}{$AtomID};
      $FeatureAtomIDsRef = $This->{StructureFeatures}{AtomIDs}{$Radius}{$AtomID};

      if ($SpecifiedFeatureAtomCount != $FeatureAtomCount) {
	# Couldn't be structurally equivalent...
	next ATOMIDENTIFIER;
      }
      for $SpecifiedFeatureAtomID (keys % {$SpecifiedFeatureAtomIDsRef}) {
	if (! exists $FeatureAtomIDsRef->{$SpecifiedFeatureAtomID}) {
	  # For structural equivalency, all atom in specified feature must also be present in a previously
	  # identified structurally unique structure feature...
	  next ATOMIDENTIFIER;
	}
      }
      # Found structurally equivalent feature...
      return ($Radius, $AtomIdentifier);
    }
  }
  return (undef, undef);
}

# Setup structure features for atom IDs involved in unique atom identifiers at all
# radii level...
#
sub _SetupStructureFeaturesForAtomIDsInvolvedInUniqueIdentifiers {
  my($This) = @_;
  my($Radius, $PreviousRadius, $Atom, $AtomID, $AtomIdentifier, $NeighborhoodAtomID, $NeighborhoodAtomsWithSuccessorAtomsRef, $NeighborhoodAtom, $NeighborhoodAtomSuccessorAtomsRef, %AtomIDs);

  $This->_InitializeStructureFeatures();

  # Collect atom IDs involved in unique atom identifiers...
  %AtomIDs = ();
  for $Radius (0 .. $This->{NeighborhoodRadius}) {
    for $AtomIdentifier (keys %{$This->{UniqueAtomIdentifiers}{$Radius}}) {
      $AtomID = $This->{UniqueAtomIdentifiers}{$Radius}{$AtomIdentifier};
      $AtomIDs{$AtomID} = $AtomID;
    }
  }

  # Setup structure features...
  for $Radius (0 .. $This->{NeighborhoodRadius}) {
    for $AtomID (keys %AtomIDs) {
      my($StructureFeatureAtomCount, %StructureFeatureAtomIDs);

      $StructureFeatureAtomCount = 0;
      %StructureFeatureAtomIDs = ();

      # Get partial structure features for the atom at previous radius level...
      $PreviousRadius = $Radius - 1;
      if ($PreviousRadius >= 0) {
	$StructureFeatureAtomCount += $This->{StructureFeatures}{AtomCount}{$PreviousRadius}{$AtomID};
	%StructureFeatureAtomIDs = %{$This->{StructureFeatures}{AtomIDs}{$PreviousRadius}{$AtomID}};
      }

      # Get all neighborhood atom at this radius level...
      if (exists($This->{AtomNeighborhoods}{$Radius}) && exists($This->{AtomNeighborhoods}{$Radius}{$AtomID})) {
	NEIGHBORHOODS: for $NeighborhoodAtomsWithSuccessorAtomsRef (@{$This->{AtomNeighborhoods}{$Radius}{$AtomID}}) {
	  ($NeighborhoodAtom, $NeighborhoodAtomSuccessorAtomsRef) = @{$NeighborhoodAtomsWithSuccessorAtomsRef};
	  if ($NeighborhoodAtom->IsHydrogen()) {
	    next NEIGHBORHOODS;
	  }
	  $NeighborhoodAtomID = $NeighborhoodAtom->GetID();
	  $StructureFeatureAtomCount++;
	  $StructureFeatureAtomIDs{$NeighborhoodAtomID} = $NeighborhoodAtomID;
	}
      }

      # Assign structure features to atom at this radius level...
      $This->{StructureFeatures}{AtomCount}{$Radius}{$AtomID} = $StructureFeatureAtomCount;
      $This->{StructureFeatures}{AtomIDs}{$Radius}{$AtomID} = \%StructureFeatureAtomIDs;
    }
  }
  return $This;
}

# Intialize structure features at each radius level...
#
sub _InitializeStructureFeatures {
  my($This) = @_;
  my($Radius, $CurrentRadius, $Atom, $AtomID);

  # Initialize all structure features...

  %{$This->{StructureFeatures}} = ();
  %{$This->{StructureFeatures}{AtomCount}} = ();
  %{$This->{StructureFeatures}{AtomIDs}} = ();

  $Radius = $This->{NeighborhoodRadius};
  for $CurrentRadius (0 .. $Radius) {
    # Structure features for at specific radii accessed using atom IDs...
    %{$This->{StructureFeatures}{AtomCount}{$CurrentRadius}} = ();
    %{$This->{StructureFeatures}{AtomIDs}{$CurrentRadius}} = ();
  }
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

# Initialize atom indentifier type information...
#
# Current supported values:
#
# AtomicInvariantsAtomTypes, FunctionalClassAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
# MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes
#
sub _InitializeAtomIdentifierTypeInformation {
  my($This) = @_;

  IDENTIFIERTYPE: {
    if ($This->{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
      $This->_InitializeAtomicInvariantsAtomTypesInformation();
      last IDENTIFIERTYPE;
    }
    if ($This->{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
      $This->_InitializeFunctionalClassAtomTypesInformation();
      last IDENTIFIERTYPE;
    }
    if ($This->{AtomIdentifierType} =~ /^(DREIDINGAtomTypes|EStateAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
      # Nothing to do for now...
      last IDENTIFIERTYPE;
    }
    carp "Warning: ${ClassName}->_InitializeAtomIdentifierTypeInformation: Unknown atom indentifier type $This->{AtomIdentifierType}...";
  }
  return $This;
}

# Initialize atomic invariants atom types, generated by AtomTypes::AtomicInvariantsAtomTypes
# class, to use for generating initial atom identifiers...
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
# optional.
#
# Default atomic invariants used for generating inital atom identifiers are [ Ref 24 ]:
#
#   AS, X<n>, BO<n>, H<n>, FC<+n/-n>, MN<n>
#
# In addition to usage of abbreviations for specifying atomic invariants, the following descriptive words
# are also allowed:
#
# X : NumOfNonHydrogenAtomNeighbors or NumOfHeavyAtomNeighbors
# BO : SumOfBondOrdersToNonHydrogenAtoms or SumOfBondOrdersToHeavyAtoms
# LBO : LargestBondOrderToNonHydrogenAtoms or LargestBondOrderToHeavyAtoms
# SB :  NumOfSingleBondsToNonHydrogenAtoms or NumOfSingleBondsToHeavyAtoms
# DB : NumOfDoubleBondsToNonHydrogenAtoms or NumOfDoubleBondsToHeavyAtoms
# TB : NumOfTripleBondsToNonHydrogenAtoms or NumOfTripleBondsToHeavyAtoms
# H :  NumOfImplicitAndExplicitHydrogens
# Ar : Aromatic
# RA : RingAtom
# FC : FormalCharge
# MN : MassNumber
# SM : SpinMultiplicity
#
sub _InitializeAtomicInvariantsAtomTypesInformation {
  my($This) = @_;

  # Default atomic invariants to use for generating initial atom identifiers are: AS, X, BO, LBO, H, FC
  #
  @{$This->{AtomicInvariantsToUse}} = ();
  @{$This->{AtomicInvariantsToUse}} = ('AS', 'X', 'BO', 'H', 'FC', 'MN');

  return $This;
}

# Initialize functional class atom types, generated by AtomTypes::FunctionalClassAtomTypes
# class, to use for generating initial atom identifiers...
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

  # Default functional class atom typess to use for generating initial atom identifiers
  # are: HBD, HBA, PI, NI, Ar, Hal
  #
  @{$This->{FunctionalClassesToUse}} = ();
  @{$This->{FunctionalClassesToUse}} = ('HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal');

  return $This;
}

# Set atomic invariants to use for generation of intial atom indentifiers...
#
sub SetAtomicInvariantsToUse {
  my($This, @Values) = @_;
  my($FirstValue, $TypeOfFirstValue, $AtomicInvariant, $SpecifiedAtomicInvariant, @SpecifiedAtomicInvariants, @AtomicInvariantsToUse);

  if (!@Values) {
    carp "Warning: ${ClassName}->SetAtomicInvariantsToUse: No values specified...";
    return;
  }

  if ($This->{AtomIdentifierType} !~ /^AtomicInvariantsAtomTypes$/i) {
    carp "Warning: ${ClassName}->SetFunctionalAtomTypesToUse: AtomicInvariantsToUse can't be set for InitialAtomIdentifierType of $This->{AtomIdentifierType}...";
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

# Set functional classes to use for generation of intial atom indentifiers...
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

# Return a string containg data for ExtendedConnectivityFingerprints object...
sub StringifyExtendedConnectivityFingerprints {
  my($This) = @_;
  my($ExtendedConnectivityFingerprintsString);

  $ExtendedConnectivityFingerprintsString = "InitialAtomIdentifierType: $This->{AtomIdentifierType}; NeighborhoodRadius: $This->{NeighborhoodRadius}";

  if ($This->{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
    my($AtomicInvariant, @AtomicInvariants, @AtomicInvariantsOrder, %AvailableAtomicInvariants);

    @AtomicInvariantsOrder = AtomTypes::AtomicInvariantsAtomTypes::GetAtomicInvariantsOrder();
    %AvailableAtomicInvariants = AtomTypes::AtomicInvariantsAtomTypes::GetAvailableAtomicInvariants();

    for $AtomicInvariant (@AtomicInvariantsOrder) {
      push @AtomicInvariants, "$AtomicInvariant: $AvailableAtomicInvariants{$AtomicInvariant}";
    }

    $ExtendedConnectivityFingerprintsString .= "; AtomicInvariantsToUse: <" . TextUtil::JoinWords(\@{$This->{AtomicInvariantsToUse}}, ", ", 0) . ">";
    $ExtendedConnectivityFingerprintsString .= "; AtomicInvariantsOrder: <" . TextUtil::JoinWords(\@AtomicInvariantsOrder, ", ", 0) . ">";
    $ExtendedConnectivityFingerprintsString .= "; AvailableAtomicInvariants: <" . TextUtil::JoinWords(\@AtomicInvariants, ", ", 0) . ">";
  }
  elsif ($This->{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
    my($FunctionalClass, @FunctionalClasses, @FunctionalClassesOrder, %AvailableFunctionalClasses);

    @FunctionalClassesOrder = AtomTypes::FunctionalClassAtomTypes::GetFunctionalClassesOrder();
    %AvailableFunctionalClasses = AtomTypes::FunctionalClassAtomTypes::GetAvailableFunctionalClasses();

    for $FunctionalClass (@FunctionalClassesOrder) {
      push @FunctionalClasses, "$FunctionalClass: $AvailableFunctionalClasses{$FunctionalClass}";
    }

    $ExtendedConnectivityFingerprintsString .= "; FunctionalClassesToUse: <" . TextUtil::JoinWords(\@{$This->{FunctionalClassesToUse}}, ", ", 0) . ">";
    $ExtendedConnectivityFingerprintsString .= "; FunctionalClassesOrder: <" . TextUtil::JoinWords(\@FunctionalClassesOrder, ", ", 0) . ">";
    $ExtendedConnectivityFingerprintsString .= "; AvailableFunctionalClasses: <" . TextUtil::JoinWords(\@FunctionalClasses, ", ", 0) . ">";
  }

  if ($This->{Type} =~ /^ExtendedConnectivityBits$/i) {
    # Size...
    $ExtendedConnectivityFingerprintsString .= "; Size: $This->{Size}; MinSize: $This->{MinSize}; MaxSize: $This->{MaxSize}";

    # Fingerprint bit density and num of bits set...
    my($NumOfSetBits, $BitDensity);
    $NumOfSetBits = $This->{FingerprintsBitVector}->GetNumOfSetBits();
    $BitDensity = $This->{FingerprintsBitVector}->GetFingerprintsBitDensity();
    $ExtendedConnectivityFingerprintsString .= "; NumOfOnBits: $NumOfSetBits; BitDensity: $BitDensity";

    $ExtendedConnectivityFingerprintsString .= "; FingerprintsBitVector: < $This->{FingerprintsBitVector} >";
  }
  else {
    # Number of identifiers...
    $ExtendedConnectivityFingerprintsString .= "; NumOfIdentifiers: " . $This->{FingerprintsVector}->GetNumOfValues();

    # FingerprintsVector...
    $ExtendedConnectivityFingerprintsString .= "; FingerprintsVector: < $This->{FingerprintsVector} >";
  }

  return $ExtendedConnectivityFingerprintsString;
}

1;

__END__

=head1 NAME

ExtendedConnectivityFingerprints

=head1 SYNOPSIS

use Fingerprints::ExtendedConnectivityFingerprints;

use Fingerprints::ExtendedConnectivityFingerprints qw(:all);

=head1 DESCRIPTION

ExtendedConnectivityFingerprints  [ Ref 48, Ref 52 ] class provides the following methods:

new, GenerateFingerprints, GetDescription, SetAtomIdentifierType,
SetAtomicInvariantsToUse, SetFunctionalClassesToUse, SetNeighborhoodRadius,
StringifyExtendedConnectivityFingerprints

B<ExtendedConnectivityFingerprints> is derived from B<Fingerprints> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<ExtendedConnectivityFingerprints>, B<Fingerprints> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The current release of MayaChemTools supports generation of B<ExtendedConnectivityFingerprints>
corresponding to following B<AtomtomIdentifierTypes>:

    AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
    FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
    SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes

Based on the values specified for B<AtomIdentifierType>, B<AtomicInvariantsToUse>
and B<FunctionalClassesToUse>, initial atom types are assigned to all non-hydrogen atoms in
a molecule and these atom types strings are converted into initial atom identifier integers using
B<TextUtil::HashCode> function. The duplicate atom identifiers are removed.

For B<NeighborhoodRadius> value of I<0>, the initial set of unique atom identifiers comprises
the molecule fingerprints. Otherwise, atom neighborhoods are generated for each non-hydrogen
atom up-to specified B<NeighborhoodRadius> value. For each non-hydrogen central atom
at a specific radius, its neighbors at next radius level along with their bond orders and previously
calculated atom identifiers are collected which in turn are used to generate a new integer
atom identifier; the bond orders and atom identifier pairs list is first sorted by bond order
followed by atom identifiers to make these values graph invariant.

After integer atom identifiers have been generated for all non-hydrogen atoms at all specified
neighborhood radii, the duplicate integer atom identifiers corresponding to same hash code
value generated using B<TextUtil::HashCode> are tracked by keeping the atom identifiers at
lower radius. Additionally, all structurally duplicate integer atom identifiers at each specified
radius are also tracked by identifying equivalent atom and bonds corresponding to substructures
used for generating atom identifier and keeping integer atom identifier with lowest value.

For I<ExtendedConnnectivity> value of fingerprints B<Type>, the duplicate identifiers are
removed from the list and the unique atom identifiers constitute the extended connectivity
fingerprints of a molecule.

For I<ExtendedConnnectivityCount> value of fingerprints B<Type>, the occurrence of each
unique atom identifiers appears is counted and the unique atom identifiers along with their
count constitute the extended connectivity fingerprints of a molecule.

For I<ExtendedConnectivityBits> value of fingerprints B<-m, --mode>, the unique atom identifiers
are used as a random number seed to generate a random integer value between 0 and B<--Size> which
in turn is used to set corresponding bits in the fingerprint bit-vector string.

The current release of MayaChemTools generates the following types of extended connectivity
fingerprints vector strings:

    FingerprintsVector;ExtendedConnectivity:AtomicInvariantsAtomTypes:Radi
    us2;60;AlphaNumericalValues;ValuesString;73555770 333564680 352413391
    666191900 1001270906 1371674323 1481469939 1977749791 2006158649 21414
    08799 49532520 64643108 79385615 96062769 273726379 564565671 85514103
    5 906706094 988546669 1018231313 1032696425 1197507444 1331250018 1338
    532734 1455473691 1607485225 1609687129 1631614296 1670251330 17303...

    FingerprintsVector;ExtendedConnectivityCount:AtomicInvariantsAtomTypes
    :Radius2;60;NumericalValues;IDsAndValuesString;73555770 333564680 3524
    13391 666191900 1001270906 1371674323 1481469939 1977749791 2006158649
    2141408799 49532520 64643108 79385615 96062769 273726379 564565671...;
    3 2 1 1 14 1 2 10 4 3 1 1 1 1 2 1 2 1 1 1 2 3 1 1 2 1 3 3 8 2 2 2 6 2
    1 2 1 1 2 1 1 1 2 1 1 2 1 2 1 1 1 1 1 1 1 1 1 2 1 1

    FingerprintsBitVector;ExtendedConnectivityBits:AtomicInvariantsAtomTyp
    es:Radius2;1024;BinaryString;Ascending;0000000000000000000000000000100
    0000000001010000000110000011000000000000100000000000000000000000100001
    1000000110000000000000000000000000010011000000000000000000000000010000
    0000000000000000000000000010000000000000000001000000000000000000000000
    0000000000010000100001000000000000101000000000000000100000000000000...

    FingerprintsBitVector;ExtendedConnectivityBits:AtomicInvariantsAtomTyp
    es:Radius2;1024;HexadecimalString;Ascending;000000010050c0600800000803
    0300000091000004000000020000100000000124008200020000000040020000000000
    2080000000820040010020000000008040000000000080001000000000400000000000
    4040000090000061010000000800200000000000001400000000020080000000000020
    00008020200000408000

    FingerprintsVector;ExtendedConnectivity:FunctionalClassAtomTypes:Radiu
    s2;57;AlphaNumericalValues;ValuesString;24769214 508787397 850393286 8
    62102353 981185303 1231636850 1649386610 1941540674 263599683 32920567
    1 571109041 639579325 683993318 723853089 810600886 885767127 90326012
    7 958841485 981022393 1126908698 1152248391 1317567065 1421489994 1455
    632544 1557272891 1826413669 1983319256 2015750777 2029559552 20404...

    FingerprintsVector;ExtendedConnectivityCount:FunctionalClassAtomTypes:
    Radius2;57;NumericalValues;IDsAndValuesString;24769214 508787397 85039
    3286 862102353 981185303 1231636850 1649386610 1941540674 263599683 32
    9205671 571109041 639579325 683993318 723853089 810600886 885767127...;
    1 1 1 10 2 22 3 1 3 3 1 1 1 3 2 2 1 2 2 2 3 1 1 1 1 1 14 1 1 1 1 1 1 2
    1 2 1 1 2 2 1 1 2 1 1 1 2 1 1 2 1 1 1 1 1 1 1

    FingerprintsBitVector;ExtendedConnectivityBits:FunctionalClassAtomType
    s:Radius2;1024;BinaryString;Ascending;00000000000000000000100000000000
    0000000001000100000000001000000000000000000000000000000000101000000010
    0000001000000000010000000000000000000000000000000000000000000000000100
    0000000000001000000000000001000000000001001000000000000000000000000000
    0000000000000000100000000000001000000000000000000000000000000000000...

    FingerprintsVector;ExtendedConnectivity:DREIDINGAtomTypes:Radius2;56;A
    lphaNumericalValues;ValuesString;280305427 357928343 721790579 1151822
    898 1207111054 1380963747 1568213839 1603445250 4559268 55012922 18094
    0813 335715751 534801009 684609658 829361048 972945982 999881534 10076
    55741 1213692591 1222032501 1224517934 1235687794 1244268533 152812070
    0 1629595024 1856308891 1978806036 2001865095 2096549435 172675415 ...

    FingerprintsVector;ExtendedConnectivity:EStateAtomTypes:Radius2;62;Alp
    haNumericalValues;ValuesString;25189973 528584866 662581668 671034184
    926543080 1347067490 1738510057 1759600920 2034425745 2097234755 21450
    44754 96779665 180364292 341712110 345278822 386540408 387387308 50430
    1706 617094135 771528807 957666640 997798220 1158349170 1291258082 134
    1138533 1395329837 1420277211 1479584608 1486476397 1487556246 1566...

    FingerprintsVector;ExtendedConnectivity:MMFF94AtomTypes:Radius2;64;Alp
    haNumericalValues;ValuesString;224051550 746527773 998750766 103704190
    2 1239701709 1248384926 1259447756 1521678386 1631549126 1909437580 20
    37095052 2104274756 2117729376 8770364 31445800 81450228 314289324 344
    041929 581773587 638555787 692022098 811840536 929651561 936421792 988
    636432 1048624296 1054288509 1369487579 1454058929 1519352190 17271...

    FingerprintsVector;ExtendedConnectivity:SLogPAtomTypes:Radius2;71;Alph
    aNumericalValues;ValuesString;78989290 116507218 489454042 888737940 1
    162561799 1241797255 1251494264 1263717127 1471206899 1538061784 17654
    07295 1795036542 1809833874 2020454493 2055310842 2117729376 11868981
    56731842 149505242 184525155 196984339 288181334 481409282 556716568 6
    41915747 679881756 721736571 794256218 908276640 992898760 10987549...

    FingerprintsVector;ExtendedConnectivity:SYBYLAtomTypes:Radius2;58;Alph
    aNumericalValues;ValuesString;199957044 313356892 455463968 465982819
    1225318176 1678585943 1883366064 1963811677 2117729376 113784599 19153
    8837 196629033 263865277 416380653 477036669 681527491 730724924 90906
    5537 1021959189 1133014972 1174311016 1359441203 1573452838 1661585138
    1668649038 1684198062 1812312554 1859266290 1891651106 2072549404 ...

    FingerprintsVector;ExtendedConnectivity:TPSAAtomTypes:Radius2;47;Alpha
    NumericalValues;ValuesString;20818206 259344053 862102353 1331904542 1
    700688206 265614156 363161397 681332588 810600886 885767127 950172500
    951454814 1059668746 1247054493 1382302230 1399502637 1805025917 19189
    39561 2114677228 2126402271 8130483 17645742 32278373 149975755 160327
    654 256360355 279492740 291251259 317592700 333763396 972105960 101...

    FingerprintsVector;ExtendedConnectivity:UFFAtomTypes:Radius2;56;AlphaN
    umericalValues;ValuesString;280305427 357928343 721790579 1151822898 1
    207111054 1380963747 1568213839 1603445250 4559268 55012922 180940813
    335715751 534801009 684609658 829361048 972945982 999881534 1007655741
    1213692591 1222032501 1224517934 1235687794 1244268533 1528120700 162
    9595024 1856308891 1978806036 2001865095 2096549435 172675415 18344...

=head2 METHODS

=over 4

=item B<new>

    $NewExtendedConnectivityFingerprints = new ExtendedConnectivityFingerprints(
                                                   %NamesAndValues);

Using specified I<ExtendedConnectivityFingerprints> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<ExtendedConnectivityFingerprints>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'ExtendedConnectivity'
    NeighborhoodRadius = 2
    AtomIdentifierType = ''
    AtomicInvariantsToUse = ['AS', 'X', 'BO', 'H', 'FC', 'MN']
    FunctionalClassesToUse = ['HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal']

Examples:

    $ExtendedConnectivityFingerprints = new ExtendedConnectivityFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes');

    $ExtendedConnectivityFingerprints = new ExtendedConnectivityFingerprints(
                              'Type' => 'ExtendedConnectivityCount',
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes');

    $ExtendedConnectivityFingerprints = new ExtendedConnectivityFingerprints(
                              'Type' => 'ExtendedConnectivityBits',
                              'Molecule' => $Molecule,
                              'Size' => 1024,
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes');

    $ExtendedConnectivityFingerprints = new ExtendedConnectivityFingerprints(
                              'Type' => 'ExtendedConnectivity',
                              'Molecule' => $Molecule,
                              'NeighborhoodRadius' => 2,
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes',
                              'AtomicInvariantsToUse' =>
                                              ['AS', 'X', 'BO', 'H', 'FC', 'MN'] );

    $ExtendedConnectivityFingerprints = new ExtendedConnectivityFingerprints(
                              'Type' => 'ExtendedConnectivity',
                              'Molecule' => $Molecule,
                              'NeighborhoodRadius' => 2,
                              'AtomIdentifierType' =>
                                          'FunctionalClassAtomTypes',
                              'FunctionalClassesToUse' =>
                                          ['HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal'] );

    $ExtendedConnectivityFingerprints = new ExtendedConnectivityFingerprints(
                              'Type' => 'ExtendedConnectivity',
                              'Molecule' => $Molecule,;
                              'AtomIdentifierType' =>
                                              'MMFF94AtomTypes');

    $ExtendedConnectivityFingerprints = new ExtendedConnectivityFingerprints(
                              'Type' => 'ExtendedConnectivityCount',
                              'Molecule' => $Molecule,;
                              'AtomIdentifierType' =>
                                              'MMFF94AtomTypes');

    $ExtendedConnectivityFingerprints = new ExtendedConnectivityFingerprints(
                              'Type' => 'ExtendedConnectivityCount',
                              'Molecule' => $Molecule,;
                              'AtomIdentifierType' =>
                                              'SLogPAtomTypes');

    $ExtendedConnectivityFingerprints = new ExtendedConnectivityFingerprints(
                              'Type' => 'ExtendedConnectivity',
                              'Molecule' => $Molecule,;
                              'AtomIdentifierType' =>
                                              'SLogPAtomTypes');

    $ExtendedConnectivityFingerprints = new ExtendedConnectivityFingerprints(
                              'Type' => 'ExtendedConnectivity',
                              'Molecule' => $Molecule,;
                              'AtomIdentifierType' =>
                                              'SYBYLAtomTypes');

    $ExtendedConnectivityFingerprints->GenerateFingerprints();
    print "$ExtendedConnectivityFingerprints\n";

=item B<GenerateFingerprints>

    $ExtendedConnectivityFingerprints->GenerateFingerprints();

Generates extended connectivity fingerprints and returns I<ExtendedConnectivityFingerprints>.

=item B<GetDescription>

    $Description = $ExtendedConnectivityFingerprints->GetDescription();

Returns a string containing description of extended connectivity fingerprints
fingerprints.

=item B<SetAtomIdentifierType>

    $ExtendedConnectivityFingerprints->SetAtomIdentifierType($IdentifierType);

Sets atom I<IdentifierType> to use during extended connectivity fingerprints generation and
returns I<ExtendedConnectivityFingerprints>.

Possible values: I<AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes,
TPSAAtomTypes, UFFAtomTypes>.

=item B<SetAtomicInvariantsToUse>

    $ExtendedConnectivityFingerprints->SetAtomicInvariantsToUse($ValuesRef);
    $ExtendedConnectivityFingerprints->SetAtomicInvariantsToUse(@Values);

Sets atomic invariants to use during I<AtomicInvariantsAtomTypes> value of I<AtomIdentifierType>
for extended connectivity fingerprints generation and returns I<ExtendedConnectivityFingerprints>.

Possible values for atomic invariants are: I<AS, X, BO,  LBO, SB, DB, TB,
H, Ar, RA, FC, MN, SM>. Default value [ Ref 24 ]: I<AS,X,BO,H,FC,MN>.

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

    $ExtendedConnectivityFingerprints->SetFunctionalClassesToUse($ValuesRef);
    $ExtendedConnectivityFingerprints->SetFunctionalClassesToUse(@Values);

Sets functional classes invariants to use during I<FunctionalClassAtomTypes> value of I<AtomIdentifierType>
for extended connectivity fingerprints generation and returns I<ExtendedConnectivityFingerprints>.

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

=item B<SetNeighborhoodRadius>

    $ExtendedConnectivityFingerprints->SetNeighborhoodRadius($Radius);

Sets neighborhood radius to use during extended connectivity fingerprints generation and
returns I<ExtendedConnectivityFingerprints>.

=item B<StringifyExtendedConnectivityFingerprints>

    $String = $Fingerprints->StringifyExtendedConnectivityFingerprints();

Returns a string containing information about I<ExtendedConnectivityFingerprints> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Fingerprints.pm, FingerprintsStringUtil.pm, AtomNeighborhoodsFingerprints.pm,
AtomTypesFingerprints.pm, EStateIndiciesFingerprints.pm, MACCSKeys.pm,
PathLengthFingerprints.pm, TopologicalAtomPairsFingerprints.pm,
TopologicalAtomTripletsFingerprints.pm, TopologicalAtomTorsionsFingerprints.pm,
TopologicalPharmacophoreAtomPairsFingerprints.pm,
TopologicalPharmacophoreAtomTripletsFingerprints.pm


=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
