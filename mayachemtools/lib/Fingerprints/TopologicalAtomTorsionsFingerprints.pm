package Fingerprints::TopologicalAtomTorsionsFingerprints;
#
# File: TopologicalAtomTorsionsFingerprints.pm
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
use overload '""' => 'StringifyTopologicalAtomTorsionsFingerprints';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeTopologicalAtomTorsionsFingerprints();

  $This->_InitializeTopologicalAtomTorsionsFingerprintsProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeTopologicalAtomTorsionsFingerprints {
  my($This) = @_;

  # Type of fingerprint...
  $This->{Type} = 'TopologicalAtomTorsions';

  # Type of vector...
  $This->{VectorType} = 'FingerprintsVector';

  # Type of FingerprintsVector...
  $This->{FingerprintsVectorType} = 'NumericalValues';

  # Atom identifier type to use for atom IDs in atom torsions...
  #
  # Currently supported values are: AtomicInvariantsAtomTypes, DREIDINGAtomTypes,
  # EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
  # SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes
  #
  $This->{AtomIdentifierType} = '';

  # Atom types assigned to each heavy atom...
  #
  %{$This->{AssignedAtomTypes}} = ();

  # Final unique atom torsions...
  #
  @{$This->{AtomTorsionsIDs}} = ();
  %{$This->{AtomTorsionsCount}} = ();
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize object properties....
sub _InitializeTopologicalAtomTorsionsFingerprintsProperties {
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

  return "$This->{Type}:$This->{AtomIdentifierType}";
}

# Generate topological atom torsions [ Ref 58, Ref 72 ] fingerprints...
#
# Methodology:
#   . Assign atom types to all the atoms.
#   . Generate and count atom torsions.
#
# Notes:
#   . Hydrogen atoms are ignored during the fingerprint generation.
#
sub GenerateFingerprints {
  my($This) = @_;

  # Cache appropriate molecule data...
  $This->_SetupMoleculeDataCache();

  # Assign atom types to all heavy atoms...
  if (!$This->_AssignAtomTypes()) {
    carp "Warning: ${ClassName}->GenerateFingerprints: $This->{AtomIdentifierType} fingerprints generation didn't succeed: Couldn't assign valid $This->{AtomIdentifierType} to all atoms...";
    return $This;
  }

  # Count atom torsions...
  $This->_GenerateAndCountAtomTorsions();

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

# Count atom torsions involving non-hydrogen atoms by going over the structurally
# unique atom torsions...
#
sub _GenerateAndCountAtomTorsions {
  my($This) = @_;
  my($Atom1, $Atom2, $Atom3, $Atom4, $AtomID1, $AtomID2, $AtomID3, $AtomID4, $AtomTorsionID, @Atom1Neighbors, @Atom2Neighbors, @Atom3Neighbors);

  # Setup a hash to track structurally unique atom torsions by atom IDs...
  %{$This->{StructurallyUniqueAtomTorsions}} = ();

  ATOM1: for $Atom1 (@{$This->{Atoms}}) {
    if ($Atom1->IsHydrogen()) {
      next ATOM1;
    }
    $AtomID1 = $Atom1->GetID();
    # Go over Atom1 neighbors other than Atom1...
    @Atom1Neighbors = $Atom1->GetNeighbors($Atom1);
    ATOM2: for $Atom2 (@Atom1Neighbors) {
      if ($Atom2->IsHydrogen()) {
	next ATOM2;
      }
      $AtomID2 = $Atom2->GetID();
      # Go over Atom2 neighbors other than Atom1 and Atom2...
      @Atom2Neighbors = $Atom2->GetNeighbors($Atom1, $Atom2);
      ATOM3: for $Atom3 (@Atom2Neighbors) {
	if ($Atom3->IsHydrogen()) {
	  next ATOM3;
	}
	$AtomID3 = $Atom3->GetID();
	@Atom3Neighbors = $Atom3->GetNeighbors($Atom1, $Atom2, $Atom3);
	# Go over Atom3 neighbors other than Atom1, Atom2 and Atom3...
	ATOM4: for $Atom4 (@Atom3Neighbors) {
	  if ($Atom4->IsHydrogen()) {
	    next ATOM4;
	  }
	  $AtomID4 = $Atom4->GetID();

	  # Is it a structurally unique torsion?
	  if (!$This->_IsStructurallyUniqueTorsion($AtomID1, $AtomID2, $AtomID3, $AtomID4)) {
	    next ATOM4;
	  }

	  # Track structurally unique torsions...
	  $AtomTorsionID = $This->_GetAtomTorsionID($AtomID1, $AtomID2, $AtomID3, $AtomID4);
	  if (exists $This->{AtomTorsionsCount}{$AtomTorsionID}) {
	    $This->{AtomTorsionsCount}{$AtomTorsionID} += 1;
	  }
	  else {
	    $This->{AtomTorsionsCount}{$AtomTorsionID} = 1;
	  }
	}
      }
    }
  }

  return $This;
}

# Is it a structurally unique torsions?
#
# Notes:
#  . For a torsion to be structurally unique which hasn't already been encountered,
#    all the four atoms involved in the torsion must be new atoms. And this can be
#    simply implemented by tracking the torsions using atom IDs and maintaining a
#    hash of already encountered torsions using lexicographically smaller torsion ID
#    consisting of four atom IDs.
#
sub _IsStructurallyUniqueTorsion {
  my($This, @AtomIDs) = @_;
  my($TorsionID, $ReverseTorsionID);

  $TorsionID = join "-", @AtomIDs;
  $ReverseTorsionID = join "-", reverse @AtomIDs;

  # Use lexicographically smaller string...
  if ($ReverseTorsionID lt $TorsionID) {
    $TorsionID = $ReverseTorsionID;
  }

  if (exists $This->{StructurallyUniqueAtomTorsions}{$TorsionID}) {
    return 0;
  }

  # Keep track...
  $This->{StructurallyUniqueAtomTorsions}{$TorsionID} = 1;

  return 1;
}

# Get atom torsion ID corresponding to atom types involved in torsion...
#
# Notes:
#  . TorsionID corresponds to assigned atom types of all the four torsion atoms
#    concatenated by hyphen.
#  . TorsionIDs are generated for both forward and backward sequence of atoms
#    in the torsion and keeping the lexicographically smaller TorsionID to keep TorsionID
#    independent of atom ordering.
#
sub _GetAtomTorsionID {
  my($This, @AtomIDs) = @_;
  my($AtomTorsionID, $ReverseAtomTorsionID, @AtomTypes);

  @AtomTypes = ();
  @AtomTypes = map { $This->{AssignedAtomTypes}{$_} } @AtomIDs;

  $AtomTorsionID = join "-", @AtomTypes;
  $ReverseAtomTorsionID = join "-", reverse @AtomTypes;

  # Use lexicographically smaller string as ID...
  return ($ReverseAtomTorsionID lt $AtomTorsionID) ? $ReverseAtomTorsionID : $AtomTorsionID;
}

# Set final fingerpritns vector...
#
sub _SetFinalFingerprints {
  my($This) = @_;
  my($AtomTorsionID, $Value, @Values);

  # Mark successful generation of fingerprints...
  $This->{FingerprintsGenerated} = 1;

  @Values = ();
  @{$This->{AtomTorsionsIDs}} = ();

  for $AtomTorsionID (sort keys %{$This->{AtomTorsionsCount}}) {
    $Value = $This->{AtomTorsionsCount}{$AtomTorsionID};
    push @{$This->{AtomTorsionsIDs}}, $AtomTorsionID;
    push @Values, $Value;
  }

  # Add AtomPairsIDs and values to fingerprint vector...
  $This->{FingerprintsVector}->AddValueIDs(\@{$This->{AtomTorsionsIDs}});
  $This->{FingerprintsVector}->AddValues(\@Values);

  return $This;
}

# Get atom torsions IDs corresponding to atom torsions count values in fingerprint
# vector as an array or reference to an array...
#
# AtomTorsionsIDs list differes in molecules and is generated during finalization
# of fingerprints to make sure the fingerprint vector containing count values
# matches the atom torsions array.
#
sub GetAtomTorsionsIDs {
  my($This) = @_;

  return wantarray ? @{$This->{AtomTorsionsIDs}} : \@{$This->{AtomTorsionsIDs}};
}

# Cache  appropriate molecule data...
#
sub _SetupMoleculeDataCache {
  my($This) = @_;

  # Get all atoms including hydrogens. The hydrogen atoms are ignored during processing...
  @{$This->{Atoms}} = $This->GetMolecule()->GetAtoms();

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

# Initialize atomic invariants to use for generating atom IDs in atom torsions...
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
#   AtomTypeIDz = Atomic invariants atom type for atom z
#   AtomTypeIDw = Atomic invariants atom type for atom w
#
# Then:
#
#   Atom torsion AtomID generated by AtomTypes::AtomicInvariantsAtomTypes class corresponds to:
#
#     AS.X<n>.BO<n>.LBO<n>.<SB><n>.<DB><n>.<TB><n>.H<n>.Ar.RA.FC<+n/-n>.MN<n>.SM<n>
#
#  AtomTorsion ID corresponds to:
#
#    AtomTypeIDx-AtomTypeIDy-AtomTypeIDz-AtomTypeIDw
#
# Except for AS which is a required atomic invariant in atom torsions AtomIDs, all other atomic invariants are
# optional. Default atomic invariants used for AtomID are: AS, X<n>, BO<n>, H<n>, FC<+n/-n>.
# AtomID specification doesn't include atomic invariants with zero or undefined values.
#
# Examples of atom torsion AtomIDs in Aspirin using default atomic invariants:
#
#  C.X1.BO1.H3-C.X3.BO4-O.X2.BO2-C.X3.BO4
#  C.X2.BO3.H1-C.X2.BO3.H1-C.X2.BO3.H1-C.X2.BO3.H1
#  C.X3.BO4-C.X3.BO4-O.X2.BO2-C.X3.BO4
#  C.X3.BO4-O.X2.BO2-C.X3.BO4-O.X1.BO2
#
# Examples of atom torsion AtomIDs in Aspirin using AS, X and BO atomic invariants:
#
#  C.X1.BO1-C.X3.BO4-O.X2.BO2-C.X3.BO4
#  C.X2.BO3-C.X2.BO3-C.X2.BO3-C.X2.BO3
#  C.X3.BO4-C.X3.BO4-O.X2.BO2-C.X3.BO4
#  C.X3.BO4-O.X2.BO2-C.X3.BO4-O.X1.BO2
#
sub _InitializeAtomicInvariantsAtomTypesInformation {
  my($This) = @_;

  # Default atomic invariants to use for generating atom torsions atom IDs: AS, X, BO, H, FC
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

# Return a string containg data for TopologicalAtomTorsionsFingerprints object...
#
sub StringifyTopologicalAtomTorsionsFingerprints {
  my($This) = @_;
  my($FingerprintsString);

  # Type of fingerprint...
  $FingerprintsString = "Fingerprint type: $This->{Type}; AtomIdentifierType: $This->{AtomIdentifierType}";

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

  # Total number of atom torsions...
  $FingerprintsString .= "; NumOfAtomTorsions: " . $This->{FingerprintsVector}->GetNumOfValues();

  # FingerprintsVector...
  $FingerprintsString .= "; FingerprintsVector: < $This->{FingerprintsVector} >";

  return $FingerprintsString;
}

1;

__END__

=head1 NAME

TopologicalAtomTorsionsFingerprints

=head1 SYNOPSIS

use Fingerprints::TopologicalAtomTorsionsFingerprints;

use Fingerprints::TopologicalAtomTorsionsFingerprints qw(:all);

=head1 DESCRIPTION

B<TopologicalAtomTorsionsFingerprints> class provides the following methods:

new, GenerateFingerprints, GetAtomTorsionsIDs, GetDescription,
SetAtomIdentifierType, SetAtomicInvariantsToUse, SetFunctionalClassesToUse,
StringifyTopologicalAtomTorsionsFingerprints

B<TopologicalAtomTorsionsFingerprints> is derived from B<Fingerprints> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<TopologicalAtomTorsionsFingerprints>, B<Fingerprints> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The current release of MayaChemTools supports generation of B<TopologicalAtomTorsionsFingerprints>
corresponding to following B<AtomtomIdentifierTypes>:

    AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
    FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
    SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes

Based on the values specified for B<AtomIdentifierType> along with other specified
parameters such as B<AtomicInvariantsToUse> and B<FunctionalClassesToUse>, initial
atom types are assigned to all non-hydrogen  in a molecule. All unique atom torsions
are identified and an atom torsion identifier is generated; the format of atom torsion identifier is:

    <AtomType1>-<AtomType2>-<AtomType3>-<AtomType4>

    AtomType1, AtomType2, AtomType3, AtomTyp4: Assigned atom types

    where AtomType1 <= AtomType2 <= AtomType3 <= AtomType4

The atom torsion identifiers for all unique atom torsions corresponding to non-hydrogen atoms constitute
topological atom torsions fingerprints of the molecule.

The current release of MayaChemTools generates the following types of topological atom torsions
fingerprints vector strings:

    FingerprintsVector;TopologicalAtomTorsions:AtomicInvariantsAtomTypes;3
    3;NumericalValues;IDsAndValuesString;C.X1.BO1.H3-C.X3.BO3.H1-C.X3.BO4-
    C.X3.BO4 C.X1.BO1.H3-C.X3.BO3.H1-C.X3.BO4-N.X3.BO3 C.X2.BO2.H2-C.X2.BO
    2.H2-C.X3.BO3.H1-C.X2.BO2.H2 C.X2.BO2.H2-C.X2.BO2.H2-C.X3.BO3.H1-O...;
    2 2 1 1 2 2 1 1 3 4 4 8 4 2 2 6 2 2 1 2 1 1 2 1 1 2 6 2 4 2 1 3 1

    FingerprintsVector;TopologicalAtomTorsions:AtomicInvariantsAtomTypes;3
    3;NumericalValues;IDsAndValuesPairsString;C.X1.BO1.H3-C.X3.BO3.H1-C.X3
    .BO4-C.X3.BO4 2 C.X1.BO1.H3-C.X3.BO3.H1-C.X3.BO4-N.X3.BO3 2 C.X2.BO2.H
    2-C.X2.BO2.H2-C.X3.BO3.H1-C.X2.BO2.H2 1 C.X2.BO2.H2-C.X2.BO2.H2-C.X3.B
    O3.H1-O.X1.BO1.H1 1 C.X2.BO2.H2-C.X2.BO2.H2-N.X3.BO3-C.X3.BO4 2 C.X2.B
    O2.H2-C.X3.BO3.H1-C.X2.BO2.H2-C.X3.BO3.H1 2 C.X2.BO2.H2-C.X3.BO3.H1...

    FingerprintsVector;TopologicalAtomTorsions:DREIDINGAtomTypes;27;Numeri
    calValues;IDsAndValuesString;C_2-C_3-C_3-C_3 C_2-C_3-C_3-O_3 C_2-C_R-C
    _R-C_3 C_2-C_R-C_R-C_R C_2-C_R-C_R-N_R C_2-N_3-C_R-C_R C_3-C_3-C_2-O_2
    C_3-C_3-C_2-O_3 C_3-C_3-C_3-C_3 C_3-C_3-C_3-N_R C_3-C_3-C_3-O_3 C_...;
    1 1 1 2 1 2 1 1 3 1 3 2 2 2 1 1 1 3 1 2 2 32 2 2 5 3 1

    FingerprintsVector;TopologicalAtomTorsions:EStateAtomTypes;36;Numerica
    lValues;IDsAndValuesString;aaCH-aaCH-aaCH-aaCH aaCH-aaCH-aaCH-aasC aaC
    H-aaCH-aasC-aaCH aaCH-aaCH-aasC-aasC aaCH-aaCH-aasC-sF aaCH-aaCH-aasC-
    ssNH aaCH-aasC-aasC-aasC aaCH-aasC-aasC-aasN aaCH-aasC-ssNH-dssC a...;
    4 4 8 4 2 2 6 2 2 2 4 3 2 1 3 3 2 2 2 1 2 1 1 1 2 1 1 1 1 1 1 1 2 1 1 2

    FingerprintsVector;TopologicalAtomTorsions:FunctionalClassAtomTypes;26
    ;NumericalValues;IDsAndValuesString;Ar-Ar-Ar-Ar Ar-Ar-Ar-Ar.HBA Ar-Ar-
    Ar-HBD Ar-Ar-Ar-Hal Ar-Ar-Ar-None Ar-Ar-Ar.HBA-Ar Ar-Ar-Ar.HBA-None Ar
    -Ar-HBD-None Ar-Ar-None-HBA Ar-Ar-None-HBD Ar-Ar-None-None Ar-Ar.H...;
    32 5 2 2 3 3 3 2 2 2 2 1 2 1 1 1 2 1 1 1 1 3 1 1 1 3

    FingerprintsVector;TopologicalAtomTorsions:MMFF94AtomTypes;43;Numerica
    lValues;IDsAndValuesString;C5A-C5B-C5B-C5A C5A-C5B-C5B-C=ON C5A-C5B-C5
    B-CB C5A-C5B-C=ON-NC=O C5A-C5B-C=ON-O=CN C5A-C5B-CB-CB C5A-CB-CB-CB C5
    A-N5-C5A-C5B C5A-N5-C5A-CB C5A-N5-C5A-CR C5A-N5-CR-CR C5B-C5A-CB-C...;
    1 1 1 1 1 2 2 2 1 1 2 2 2 2 1 1 2 1 1 2 1 2 1 1 1 2 1 1 1 2 18 2 2 1 1
    1 1 2 1 1 3 1 3

    FingerprintsVector;TopologicalAtomTorsions:SLogPAtomTypes;49;Numerical
    Values;IDsAndValuesPairsString;C1-C10-N11-C20 1 C1-C10-N11-C21 1 C1-C1
    1-C21-C21 2 C1-C11-C21-N11 2 C1-CS-C1-C10 1 C1-CS-C1-C5 1 C1-CS-C1-CS
    2 C10-C1-CS-O2 1 C10-N11-C20-C20 2 C10-N11-C21-C11 1 C10-N11-C21-C21 1
    C11-C21-C21-C20 1 C11-C21-C21-C5 1 C11-C21-N11-C20 1 C14-C18-C18-C20
    2 C18-C14-C18-C18 2 C18-C18-C14-F 2 C18-C18-C18-C18 4 C18-C18-C18-C...

    FingerprintsVector;TopologicalAtomTorsions:SYBYLAtomTypes;26;Numerical
    Values;IDsAndValuesPairsString;C.2-C.3-C.3-C.3 1 C.2-C.3-C.3-O.3 1 C.2
    -C.ar-C.ar-C.3 1 C.2-C.ar-C.ar-C.ar 2 C.2-C.ar-C.ar-N.ar 1 C.2-N.am-C.
    ar-C.ar 2 C.3-C.3-C.2-O.co2 2 C.3-C.3-C.3-C.3 3 C.3-C.3-C.3-N.ar 1 C.3
    -C.3-C.3-O.3 3 C.3-C.3-C.ar-C.ar 2 C.3-C.3-C.ar-N.ar 2 C.3-C.3-N.ar-C.
    ar 2 C.3-C.ar-C.ar-C.ar 1 C.3-C.ar-N.ar-C.3 1 C.3-C.ar-N.ar-C.ar 1 ...

    FingerprintsVector;TopologicalAtomTorsions:TPSAAtomTypes;8;NumericalVa
    lues;IDsAndValuesPairsString;N21-None-None-None 9 N7-None-None-None 4
    None-N21-None-None 10 None-N7-None-None 3 None-N7-None-O3 1 None-None-
    None-None 44 None-None-None-O3 3 None-None-None-O4 5

    FingerprintsVector;TopologicalAtomTorsions:UFFAtomTypes;27;NumericalVa
    lues;IDsAndValuesPairsString;C_2-C_3-C_3-C_3 1 C_2-C_3-C_3-O_3 1 C_2-C
    _R-C_R-C_3 1 C_2-C_R-C_R-C_R 2 C_2-C_R-C_R-N_R 1 C_2-N_3-C_R-C_R 2 C_3
    -C_3-C_2-O_2 1 C_3-C_3-C_2-O_3 1 C_3-C_3-C_3-C_3 3 C_3-C_3-C_3-N_R 1 C
    _3-C_3-C_3-O_3 3 C_3-C_3-C_R-C_R 2 C_3-C_3-C_R-N_R 2 C_3-C_3-N_R-C_R 2
     C_3-C_R-C_R-C_R 1 C_3-C_R-N_R-C_3 1 C_3-C_R-N_R-C_R 1 C_3-N_R-C_R-...

=head2 METHODS

=over 4

=item B<new>

    $NewTopologicalAtomTorsionsFingerprints = new TopologicalAtomTorsionsFingerprints(
                                                   %NamesAndValues);

Using specified I<TopologicalAtomTorsionsFingerprints> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<TopologicalAtomTorsionsFingerprints>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'TopologicalAtomTorsions'
    AtomIdentifierType = ''
    AtomicInvariantsToUse = ['AS', 'X', 'BO', 'H', 'FC']
    FunctionalClassesToUse = ['HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal']

Examples:

    $TopologicalAtomTorsionsFingerprints = new TopologicalAtomTorsionsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes');

    $TopologicalAtomTorsionsFingerprints = new TopologicalAtomTorsionsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes',
                              'AtomicInvariantsToUse' =>
                                              ['AS', 'X', 'BO', 'H', 'FC'] );

    $TopologicalAtomTorsionsFingerprints = new TopologicalAtomTorsionsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'DREIDINGAtomTypes');

    $TopologicalAtomTorsionsFingerprints = new TopologicalAtomTorsionsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'SYBYLAtomTypes');

    $TopologicalAtomTorsionsFingerprints = new TopologicalAtomTorsionsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'SLogPAtomTypes');

    $TopologicalAtomTorsionsFingerprints = new TopologicalAtomTorsionsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomIdentifierType' =>
                                              'FunctionalClassAtomTypes',
                              'FunctionalClassesToUse' =>
                                              ['HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal'] );


    $TopologicalAtomTorsionsFingerprints->GenerateFingerprints();
    print "$TopologicalAtomTorsionsFingerprints\n";

=item B<GetDescription>

    $Description = $TopologicalAtomTorsionsFingerprints->GetDescription();

Returns a string containing description of topological atom torsions fingerprints.

=item B<GenerateFingerprints>

    $TopologicalAtomTorsionsFingerprints->GenerateFingerprints();

Generates topological atom torsions fingerprints and returns I<TopologicalAtomTorsionsFingerprints>.

=item B<GetAtomTorsionsIDs>

    $AtomPairIDsRef = $TopologicalAtomTorsionsFingerprints->GetAtomTorsionsIDs();
    @AtomPairIDs = $TopologicalAtomTorsionsFingerprints->GetAtomTorsionsIDs();

Returns atom torsion IDs corresponding to atom torsion count values in topological atom torsions
fingerprints vector as an array or reference to an array.

=item B<SetAtomIdentifierType>

    $TopologicalAtomTorsionsFingerprints->SetAtomIdentifierType($IdentifierType);

Sets atom I<IdentifierType> to use during atom torsions fingerprints generation and
returns I<TopologicalAtomTorsionsFingerprints>.

Possible values: I<AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes,
TPSAAtomTypes, UFFAtomTypes>.

=item B<SetAtomicInvariantsToUse>

    $TopologicalAtomTorsionsFingerprints->SetAtomicInvariantsToUse($ValuesRef);
    $TopologicalAtomTorsionsFingerprints->SetAtomicInvariantsToUse(@Values);

Sets atomic invariants to use during I<AtomicInvariantsAtomTypes> value of I<AtomIdentifierType>
for topological atom torsions fingerprints generation and returns I<TopologicalAtomTorsionsFingerprints>.

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

    $TopologicalTorsionsFingerprints->SetFunctionalClassesToUse($ValuesRef);
    $TopologicalTorsionsFingerprints->SetFunctionalClassesToUse(@Values);

Sets functional classes invariants to use during I<FunctionalClassAtomTypes> value of I<AtomIdentifierType>
for topological atom torsions fingerprints generation and returns I<TopologicalAtomTorsionsFingerprints>.

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

=item B<StringifyTopologicalAtomTorsionsFingerprints>

    $String = $TopologicalAtomTorsionsFingerprints->
                  StringifyTopologicalAtomTorsionsFingerprints();

Returns a string containing information about I<TopologicalAtomTorsionsFingerprints> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Fingerprints.pm, FingerprintsStringUtil.pm, AtomNeighborhoodsFingerprints.pm,
AtomTypesFingerprints.pm, EStateIndiciesFingerprints.pm, ExtendedConnectivityFingerprints.pm,
MACCSKeys.pm, PathLengthFingerprints.pm, TopologicalAtomPairsFingerprints.pm,
TopologicalAtomTripletsFingerprints.pm, TopologicalPharmacophoreAtomPairsFingerprints.pm,
TopologicalPharmacophoreAtomTripletsFingerprints.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
