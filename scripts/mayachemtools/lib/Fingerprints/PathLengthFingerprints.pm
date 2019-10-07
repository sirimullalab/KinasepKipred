package Fingerprints::PathLengthFingerprints;
#
# File: PathLengthFingerprints.pm
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
use overload '""' => 'StringifyPathLengthFingerprints';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializePathLengthFingerprints();

  $This->_InitializePathLengthFingerprintsProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializePathLengthFingerprints {
  my($This) = @_;

  # Type of fingerprint to generate...
  #
  # PathLengthBits - A bit vector indicating presence/absence of atom paths
  # PathLengthCount - A vector containing count of atom paths
  #
  $This->{Type} = '';

  # Type of vector: FingerprintsBitVector or FingerprintsVector
  $This->{VectorType} = '';

  # Set default mininum, maximum, and default size. Although any arbitrary size can
  # be specified, bit vector used to store bits work on a vector size which is
  # power of 2 and additonal bits are automatically added and cleared.
  #
  $This->{Size} = 1024;

  $This->{MinSize} = 32;
  $This->{MaxSize} = 2**32;

  # Minimum and maximum path lengths to use for fingerprints generation...
  $This->{MinLength} = 1;
  $This->{MaxLength} = 8;

  # Numner of bits to set for each atom path for FingerprintsBitVector...
  $This->{NumOfBitsToSetPerPath} = 1;

  # Atom identifier type to use for path atoms during fingerprints generation...
  #
  # Currently supported values are: AtomicInvariantsAtomTypes, DREIDINGAtomTypes,
  # EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
  # SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes
  #
  $This->{AtomIdentifierType} = '';

  # Atom types assigned to atoms...
  %{$This->{AssignedAtomTypes}} = ();

  # For molecules containing rings, atom paths starting from each atom can be traversed in four
  # different ways:
  #
  # . Atom paths without any rings and sharing of bonds in traversed paths.
  # . Atom paths containing rings and without any sharing of bonds in traversed paths
  # . All possible atom paths without any rings and sharing of bonds in traversed paths
  # . All possible atom paths containing rings and with sharing of bonds in traversed paths.
  #
  # Atom path traversal is terminated at the last ring atom. For molecules containing no rings,
  # first two and last two types described above are equivalent.
  #
  # AllowSharedBonds and AllowRings variables allow generation of differen types of paths
  # to be used for fingerprints generation.
  #
  # In addition to atom symbols, bond symbols are also used to generate a string
  # for atom paths. These atom paths strings are hased to a 32 bit integer key which
  # in turn is used as a seed for a random number generation in range of 1 to fingerprint
  # size for setting corresponding bit in bit vector.
  #
  # UseBondSymbols variable allow generation of atom path strings and consequently fingerprints.
  #
  # Combination of AllowSharedBonds, AllowRings, and UseBondSymbols allow generation of
  # 8 different types of path length fingerprints:
  #
  # AllowSharedBonds    AllowRings    UseBondSymbols    PathLengthFingerprintsType
  #
  # No                  No            Yes                AtomPathsNoCyclesWithBondSymbols
  # No                  Yes           Yes                AtomPathsWithCyclesWithBondSymbols
  #
  # Yes                 No            Yes                AllAtomPathsNoCyclesWithBondSymbols
  # Yes                 Yes           Yes                AllAtomPathsWithCyclesWithBondSymbols [ DEFAULT ]
  #
  # No                  No            No                 AtomPathsNoCyclesNoBondSymbols
  # No                  Yes           No                 AtomPathsWithCyclesNoBondSymbols
  #
  # Yes                 No            No                 AllAtomPathsNoCyclesNoBondSymbols
  # Yes                 Yes           No                 AllAtomPathsWithCyclesNoWithBondSymbols
  #
  #

  # By default, atom paths starting from atoms are allowed to share bonds already traversed...
  $This->{AllowSharedBonds} = 1;

  # By default rings are included in paths...
  $This->{AllowRings} = 1;

  # By default bond symbols are included in atom path strings...
  $This->{UseBondSymbols} = 1;

  # By default only structurally unique atom paths are used for generation
  # atom path strings...
  $This->{UseUniquePaths} = 1;

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

  # Bond symbols to use during generation of atom path strings...
  %{$This->{BondOrderToSymbol}} = ();
  %{$This->{BondOrderToSymbol}} = ('1' => '', '1.5' => ':', '2' => '=', '3' => '#');

  # BondSymbols map to use for bonded atom IDs to use during atom path strings...
  %{$This->{BondSymbols}} = ();

  # Path atom IDs to remove duplicate paths...
  %{$This->{UniqueLinearAtomPathsIDs}} = ();
  %{$This->{UniqueCyclicAtomPathsIDs}} = ();

  # Reference to all the atom paths upto specified path length...
  $This->{AtomPathsRef} = '';

  # Atom paths strings created using specified atom types and bond symbols...
  %{$This->{AtomPathsStrings}} = ();
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize object properties....
sub _InitializePathLengthFingerprintsProperties {
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

  if (!exists $NamesAndValues{Type}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying Type...";
  }

  if (!exists $NamesAndValues{AtomIdentifierType}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying AtomIdentifierType...";
  }

  # Make sure it's power of 2...
  if (exists $NamesAndValues{Size}) {
    if (!TextUtil::IsNumberPowerOfNumber($NamesAndValues{Size}, 2)) {
      croak "Error: ${ClassName}->New: Specified size value, $NamesAndValues{Size}, must be power of 2...";
    }
  }

  if ($This->{Type} =~ /^PathLengthBits$/i) {
    $This->_InitializePathLengthBits();
  }
  elsif ($This->{Type} =~ /^PathLengthCount$/i) {
    $This->_InitializePathLengthCount();
  }
  else {
    croak "Error: ${ClassName}->_InitializePathLengthFingerprintsProperties: Unknown PathLength type: $This->{Type}; Supported PathLength type : PathLengthBits or PathLengthCount......";
  }

  return $This;
}

# Initialize PathLength bits...
#
sub _InitializePathLengthBits {
  my($This) = @_;

  # Vector type...
  $This->{VectorType} = 'FingerprintsBitVector';

  $This->_InitializeFingerprintsBitVector();

  return $This;
}

# Initialize PathLength key count...
#
sub _InitializePathLengthCount {
  my($This) = @_;

  # Vector type and type of values...
  $This->{VectorType} = 'FingerprintsVector';
  $This->{FingerprintsVectorType} = 'NumericalValues';

  $This->_InitializeFingerprintsVector();

  return $This;
}

# Set type...
#
sub SetType {
  my($This, $Type) = @_;

  if ($This->{Type}) {
    croak "Error: ${ClassName}->SetType: Can't change type:  It's already set...";
  }

  if ($Type =~ /^PathLengthBits$/i) {
    $This->{Type} = 'PathLengthBits';;
  }
  elsif ($Type =~ /^PathLengthCount$/i) {
    $This->{Type} = 'PathLengthCount';;
  }
  else {
    croak "Error: ${ClassName}->SetType: Unknown PathLength keys: $Type; Supported PathLength types: PathLengthBits or PathLengthCount...";
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

# Set atom identifier type to use for path length atom identifiers...
#
sub SetAtomIdentifierType {
  my($This, $IdentifierType) = @_;

  if ($IdentifierType !~ /^(AtomicInvariantsAtomTypes|DREIDINGAtomTypes|EStateAtomTypes|FunctionalClassAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
    croak "Error: ${ClassName}->SetAtomIdentifierType: Specified value, $IdentifierType, for AtomIdentifierType is not vaild. Supported types in current release of MayaChemTools: AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes, and UFFAtomTypes.";
  }

  if ($This->{AtomIdentifierType}) {
    croak "Error: ${ClassName}->SetAtomIdentifierType: Can't change atom identifier type:  It's already set...";
  }

  $This->{AtomIdentifierType} = $IdentifierType;

  # Initialize atom identifier type information...
  $This->_InitializeAtomIdentifierTypeInformation();

  return $This;
}

# Set minimum path length...
#
sub SetMinLength {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetMinLength: MinLength value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{MinLength} = $Value;

  return $This;
}

# Set maximum path length...
#
sub SetMaxLength {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetMaxLength: MaxLength value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{MaxLength} = $Value;

  return $This;
}

# Set number of bits to set for each path...
#
sub SetNumOfBitsToSetPerPath {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetNumOfBitsToSetPerPath: NumOfBitsToSetPerPath value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{NumOfBitsToSetPerPath} = $Value;

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

  return "$This->{Type}:$This->{AtomIdentifierType}:MinLength$This->{MinLength}:MaxLength$This->{MaxLength}";
}

# Generate path length fingerprints...
#
sub GenerateFingerprints {
  my($This) = @_;

  if ($This->{MinLength} > $This->{MaxLength}) {
    croak "Error: ${ClassName}->GenerateFingerprints: No fingerpritns generated: MinLength, $This->{MinLength}, must be <= MaxLength, $This->{MaxLength}...";
  }

  # Cache appropriate molecule data...
  $This->_SetupMoleculeDataCache();

  # Assign atom types to all atoms...
  if (!$This->_AssignAtomTypes()) {
    carp "Warning: ${ClassName}->GenerateFingerprints: $This->{AtomIdentifierType} fingerprints generation didn't succeed: Couldn't assign valid $This->{AtomIdentifierType} to all atoms...";
    return $This;
  }

  # Setup bond symbol map...
  if ($This->{UseBondSymbols}) {
    $This->_InitializeBondSymbols();
  }

  # Generate appropriate atom paths...
  $This->_GenerateAtomPathsUpToMaxLength();

  # Initialize atom path strings...
  $This->_InitializeAtomPathsStrings();

  # Generate appropriate atom path strings for unique atom paths...
  $This->_GenerateAtomPathsStrings();

  # Set final fingerprints...
  $This->_SetFinalFingerprints();

  # Clear cached molecule data...
  $This->_ClearMoleculeDataCache();

  return $This;
}

# Assign appropriate atom types to all atoms...
#
sub _AssignAtomTypes {
  my($This) = @_;
  my($SpecifiedAtomTypes, $Atom, $AtomID, $IgnoreHydrogens);

  %{$This->{AssignedAtomTypes}} = ();
  $IgnoreHydrogens = 0;

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
    $AtomID = $Atom->GetID();
    $This->{AssignedAtomTypes}{$AtomID} = $SpecifiedAtomTypes->GetAtomType($Atom);
  }

  return $This;
}

# Setup bond symbol map for atoms to speed up generation of path length identifiers
# during fingerprints generation...
#
sub _InitializeBondSymbols {
  my($This) = @_;
  my($Atom1, $Atom2, $AtomID1, $AtomID2, $Bond, $BondSymbol, $BondOrder);

  %{$This->{BondSymbols}} = ();

  if (!$This->{UseBondSymbols}) {
    return $This;
  }

  for $Bond ($This->{Molecule}->GetBonds()) {
    $BondOrder = $Bond->GetBondOrder();
    $BondSymbol = $Bond->IsAromatic() ? ':' : (exists($This->{BondOrderToSymbol}{$BondOrder}) ? $This->{BondOrderToSymbol}{$BondOrder} : $BondOrder);
    ($Atom1, $Atom2) = $Bond->GetAtoms();
    $AtomID1 = $Atom1->GetID(); $AtomID2 = $Atom2->GetID();
    if ($AtomID1 > $AtomID2) {
      ($AtomID1, $AtomID2) =  ($AtomID2, $AtomID1);
    }

    if (!exists $This->{BondSymbols}{$AtomID1}) {
      %{$This->{BondSymbols}{$AtomID1}} = ();
    }
    $This->{BondSymbols}{$AtomID1}{$AtomID2} = $BondSymbol;
  }
  return $This;
}

# Get appropriate atom paths with length up to MaxLength...
#
sub _GenerateAtomPathsUpToMaxLength {
  my($This) = @_;
  my($PathLength, $AllowRings, $Molecule, $AtomPathsRef);

  $PathLength = $This->{MaxLength};
  $AllowRings = $This->{AllowRings};
  $Molecule = $This->{Molecule};

  if ($This->{AllowSharedBonds}) {
    $AtomPathsRef =  $Molecule->GetAllAtomPathsWithLengthUpto($PathLength, $AllowRings);
  }
  else {
    $AtomPathsRef = $Molecule->GetAtomPathsWithLengthUpto($PathLength, $AllowRings);
  }
  $This->{AtomPathsRef} = $AtomPathsRef;

  return $This;
}

# Initialize atom paths strings at various pathlength levels...
#
sub _InitializeAtomPathsStrings {
  my($This) = @_;
  my($PathLength);

  %{$This->{AtomPathsStrings}} = ();

  for $PathLength ($This->{MinLength} .. $This->{MaxLength}) {
    %{$This->{AtomPathsStrings}{$PathLength}} = ();
  }

  return $This;
}

# Generate appropriate atom path strings for unique atom paths...
#
sub _GenerateAtomPathsStrings {
  my($This, $PathAtomsRef) = @_;
  my($PathLength, $MinPathLength, $UseUniquePaths);

  $MinPathLength = $This->{MinLength};
  $UseUniquePaths = $This->{UseUniquePaths};

  PATHATOMS: for $PathAtomsRef (@{$This->{AtomPathsRef}}) {
    $PathLength = scalar @{$PathAtomsRef};
    if ($PathLength < $MinPathLength) {
      next PATHATOMS;
    }
    if ($UseUniquePaths) {
      $This->_GenerateAtomPathStringUsingUniquePath($PathAtomsRef);
    }
    else {
      $This->_GenerateAtomPathString($PathAtomsRef);
    }
  }
  return $This;
}

# Generate atom path string using unique path...
#
sub _GenerateAtomPathStringUsingUniquePath {
  my($This, $PathAtomsRef) = @_;

  if ($This->{AllowRings} && $This->_DoesAtomPathContainsCycle($PathAtomsRef)) {
    $This->_GenerateAtomPathStringUsingUniquePathContainingCycle($PathAtomsRef);
  }
  else {
    $This->_GenerateAtomPathStringUsingUniqueLinearPath($PathAtomsRef);
  }
  return $This;
}

# Generate atom path string for specified path containing no cycle...
#
sub _GenerateAtomPathStringUsingUniqueLinearPath {
  my($This, $PathAtomsRef) = @_;

  # Is it a unique linear atom path?
  #
  if (!$This->_IsUniqueLinearAtomPath($PathAtomsRef)) {
    return $This;
  }
  $This->_GenerateAtomPathString($PathAtomsRef);

  return $This;
}

# Is it a structurally unique linear path?
#
# For a path to be structurally unique, all of its atom IDs must be diffferent from any
# earlier path atom IDs. In order to generate atom path atom ID invariant of the atom
# order in the molecule, atom IDs are sorted numerically before generating the path ID.
#
# Notes:
#   . Atom path ID doesn't reflect the order of atoms in the atom path.
#
sub _IsUniqueLinearAtomPath {
  my($This, $PathAtomsRef) = @_;
  my($AtomPathID, $PathLength, @PathAtomIDs);

  @PathAtomIDs = ();
  @PathAtomIDs = map { $_->GetID(); } @{$PathAtomsRef};

  $AtomPathID = join '-', sort { $a <=> $b } @PathAtomIDs;
  if (exists $This->{UniqueLinearAtomPathsIDs}{$AtomPathID}) {
    return 0;
  }

  # It's a unique atom path...
  $This->{UniqueLinearAtomPathsIDs}{$AtomPathID} = 1;

  return 1;
}

# Generate atom path string for specified path containing a cycle...
#
sub _GenerateAtomPathStringUsingUniquePathContainingCycle {
  my($This, $PathAtomsRef) = @_;

  # Is it a unique atom path containing a cycle?
  #
  if (!$This->_IsUniqueAtomPathContainingCycle($PathAtomsRef)) {
    return $This;
  }

  my($CycleClosingPathAtomIndex);
  ($CycleClosingPathAtomIndex) = $This->_GetAtomPathCycleClosingAtomIndex($PathAtomsRef);

  if ($CycleClosingPathAtomIndex == 0) {
    $This->_GenerateUniqueAtomPathStringForPathCycle($PathAtomsRef);
  }
  else {
    $This->_GenerateUniqueAtomPathStringForPathContainingCycle($PathAtomsRef, $CycleClosingPathAtomIndex);
  }
  return $This;
}

# Generate a unique atom path string for a cyclic path by generating atom path
# strings for all possible paths in the cycle and keeping the lexicographically smallest
# one.
#
# Although all the paths enumerated during atom path string generation are also
# present in the intial paths list, but structural uniqueness check would detect
# 'em earlier and this method ends being invoked only once for the first cyclic path.
#
# For atom paths containg same atom types and bond symbols, atom path strings
# would be same for the paths.
#
sub _GenerateUniqueAtomPathStringForPathCycle {
  my($This, $PathAtomsRef) = @_;

  if ($This->_AreAllPathAtomsSymbolsSame($PathAtomsRef) && $This->_AreAllPathBondSymbolsSame($PathAtomsRef)) {
    return $This->_GenerateAtomPathString($PathAtomsRef);
  }

  # Generate all possible atom path strings and select the lexicographically smallest one...
  my($Index, $PathLength, $FinalAtomPathString, $FirstAtomPathString, $LastIndex, $FirstPartIndex, $FirstPartStartIndex, $FirstPartEndIndex, $SecondPartIndex, $SecondPartStartIndex, $SecondPartEndIndex, $AtomPathSymbolsRef, $AtomPathString, $ReverseAtomPathString, @FirstPartPathAtoms, @SecondPartPathAtoms, @PathAtoms);

  $PathLength = scalar @{$PathAtomsRef};
  $LastIndex = $PathLength - 1;

  $FinalAtomPathString = '';
  $FirstAtomPathString = 1;

  @FirstPartPathAtoms = (); @SecondPartPathAtoms = (); @PathAtoms = ();

  for $Index (0 .. ($LastIndex - 1)) {
    @FirstPartPathAtoms = (); @SecondPartPathAtoms = (); @PathAtoms = ();

    $FirstPartStartIndex = 0; $FirstPartEndIndex = $Index - 1;
    $SecondPartStartIndex = $Index; $SecondPartEndIndex = $LastIndex - 1;

    # Get first part atoms...
    for $FirstPartIndex ($FirstPartStartIndex .. $FirstPartEndIndex) {
      push @FirstPartPathAtoms, $PathAtomsRef->[$FirstPartIndex];
    }

    # Get second part atoms...
    for $SecondPartIndex ($SecondPartStartIndex .. $SecondPartEndIndex) {
      push @SecondPartPathAtoms, $PathAtomsRef->[$SecondPartIndex];
    }

    # Get final list of path atoms...
    if (@SecondPartPathAtoms) {
      push @PathAtoms, @SecondPartPathAtoms;
    }
    if (@FirstPartPathAtoms) {
      push @PathAtoms, @FirstPartPathAtoms;
    }

    # Complete the cycle by adding first atom as the last atom...
    push @PathAtoms, $PathAtomsRef->[$SecondPartStartIndex];

    # Generate atom path string...
    $AtomPathSymbolsRef = $This->_GenerateAtomPathSymbols(\@PathAtoms);

    $AtomPathString = join '', @{$AtomPathSymbolsRef};
    $ReverseAtomPathString = join '', reverse @{$AtomPathSymbolsRef};

    if ($ReverseAtomPathString le $AtomPathString) {
      $AtomPathString = $ReverseAtomPathString;
    }

    # Update final atom path string...

    if ($FirstAtomPathString) {
      $FirstAtomPathString = 0;
      $FinalAtomPathString = $AtomPathString;
    }
    else {
      if ($AtomPathString le $FinalAtomPathString) {
	$FinalAtomPathString = $AtomPathString;
      }
    }
  }

  # Set final atom path string...
  #
  if (exists $This->{AtomPathsStrings}{$PathLength}{$FinalAtomPathString}) {
    $This->{AtomPathsStrings}{$PathLength}{$FinalAtomPathString} += 1;
  }
  else {
    $This->{AtomPathsStrings}{$PathLength}{$FinalAtomPathString} = 1;
  }

  return $This;
}

#
# Generate a unique atom path string for paths containing a cycle closed by
# the specified atom index and the last atom index.
#
# The following methodology is used to generate atom path string which is
# independemt of initial atom ordering:
#   . Generate atom paths string from first atom to the atom before the first cycle
#     closing atom.
#   . Generate atom path string from atoms from first cycle closing atom index to
#     the last path atom in both forward and reverse order. And select the lexicographically
#     smallest atom path string.
#   . Combine atom path string generated in first step with second step to generate
#     final atom path string.
#
sub _GenerateUniqueAtomPathStringForPathContainingCycle {
  my($This, $PathAtomsRef, $CycleClosingAtomIndex) = @_;
  my($Index, $PathLength, $LastIndex, $LinearPartStartIndex, $LinearPartEndIndex, $CyclicPartStartIndex, $CyclicPartEndIndex, $CyclicPartAtomPathSymbolsRef, $CyclicPartAtomPathString, $ReverseCyclicPartAtomPathString, $AtomPathString, $AtomPathSymbolsRef, @CyclicPartPathAtoms, @PathAtoms);

  $PathLength = scalar @{$PathAtomsRef};
  $LastIndex = $PathLength - 1;

  @PathAtoms = ();

  # Get path atoms corresponding to linear  part of the path...
  $LinearPartStartIndex = 0; $LinearPartEndIndex = $CycleClosingAtomIndex - 1;

  for $Index ($LinearPartStartIndex .. $LinearPartEndIndex) {
    push @PathAtoms, $PathAtomsRef->[$Index];
  }

  # Get atoms correcponding to cyclic part of the path...
  @CyclicPartPathAtoms = ();
  $CyclicPartStartIndex = $CycleClosingAtomIndex; $CyclicPartEndIndex = $LastIndex;

  for $Index ($CyclicPartStartIndex .. $CyclicPartEndIndex) {
    push @CyclicPartPathAtoms, $PathAtomsRef->[$Index];
  }

  # Setup a lexicographically smaller atom path string for cyclic part...

  $CyclicPartAtomPathSymbolsRef = $This->_GenerateAtomPathSymbols(\@CyclicPartPathAtoms);
  $CyclicPartAtomPathString = join '', @{$CyclicPartAtomPathSymbolsRef};
  $ReverseCyclicPartAtomPathString = join '', reverse @{$CyclicPartAtomPathSymbolsRef};

  # Setup atom path corresponding to linear part and lexigraphicall smaller cyclic part...

  if ($ReverseCyclicPartAtomPathString le $CyclicPartAtomPathString) {
    push @PathAtoms, reverse @CyclicPartPathAtoms;
  }
  else {
    push @PathAtoms, @CyclicPartPathAtoms;
  }

  # Setup final atom path string...

  $AtomPathSymbolsRef = $This->_GenerateAtomPathSymbols(\@PathAtoms);
  $AtomPathString = join '', @{$AtomPathSymbolsRef};

  if (exists $This->{AtomPathsStrings}{$PathLength}{$AtomPathString}) {
    $This->{AtomPathsStrings}{$PathLength}{$AtomPathString} += 1;
  }
  else {
    $This->{AtomPathsStrings}{$PathLength}{$AtomPathString} = 1;
  }

  return $This;
}

# Does atom path contain a cycle?
#
# For an atom path to contain cycle, it must satisfy the following conditions:
#   . Pathlength >= 3
#   . Last atom ID is equal to first atom ID or some other atom ID besides itself
#
sub _DoesAtomPathContainsCycle {
  my($This, $PathAtomsRef) = @_;
  my($PathLength);

  $PathLength = scalar @{$PathAtomsRef};
  if ($PathLength <= 2) {
    return 0;
  }

  my($AtomIndex, $LastAtomIndex, $Atom, $AtomID, $LastAtom, $LastAtomID);

  $LastAtomIndex = $PathLength - 1;
  $LastAtom = $PathAtomsRef->[$LastAtomIndex];
  $LastAtomID = $LastAtom->GetID();

  # Look for atomID similar to last atom ID...
  for $AtomIndex (0 .. ($LastAtomIndex - 1)) {
    $Atom =  $PathAtomsRef->[$AtomIndex];
    $AtomID = $Atom->GetID();

    if ($AtomID == $LastAtomID) {
      # It's a cycle...
      return 1;
    }
  }
  return 0;
}

# Get atom path cycle closing atom index...
#
sub _GetAtomPathCycleClosingAtomIndex {
  my($This, $PathAtomsRef) = @_;
  my($AtomIndex, $LastAtomIndex, $Atom, $AtomID, $LastAtom, $LastAtomID, $PathLength);

  $PathLength = scalar @{$PathAtomsRef};

  $LastAtomIndex = $PathLength - 1;
  $LastAtom = $PathAtomsRef->[$LastAtomIndex]; $LastAtomID = $LastAtom->GetID();

  # Look for atomID similar to last atom ID...
  for $AtomIndex (0 .. ($LastAtomIndex - 1)) {
    $Atom =  $PathAtomsRef->[$AtomIndex]; $AtomID = $Atom->GetID();

    if ($AtomID == $LastAtomID) {
      # It's a cycle closing atom...
      return $AtomIndex;
    }
  }
  return undef;
}

# Is it a structurally unique path containing a cycle?
#
# For atom paths containing cycles, last atom ID is either equal to first atom ID or
# some other atom ID besides itself.
#
# In order to determine its structurally unqiue independent of initial atom ordering,
# the following methodolgy is used:
#
#   . For paths with same first and atom IDs:
#      . Remove the last atom ID from atom path
#      . Sort atom IDs in the path
#      . Add first atom ID from the sorted list to the end of list to complete the cycle
#      . Generate a atom path ID
#      . Use final path ID to track uniqueness of path containing cycle.
#
#   . For paths with last atom ID equal to some other atom ID besidies itself:
#      . Sort atom IDs in atom path
#      . Generate atom path ID and use it to track unqiueness of atom paths.
#
sub _IsUniqueAtomPathContainingCycle {
  my($This, $PathAtomsRef) = @_;
  my($PathLength, $AtomPathID, $FirstAtom, $LastAtom, $FirstAtomID, $LastAtomID, @PathAtomIDs, @SortedPathAtomIDs);

  @PathAtomIDs = ();
  @PathAtomIDs = map { $_->GetID(); } @{$PathAtomsRef};

  $PathLength = scalar @{$PathAtomsRef};

  $FirstAtom = $PathAtomsRef->[0]; $FirstAtomID = $FirstAtom->GetID();
  $LastAtom = $PathAtomsRef->[$PathLength - 1]; $LastAtomID = $LastAtom->GetID();

  if ($FirstAtomID == $LastAtomID) {
    pop @PathAtomIDs;

    @SortedPathAtomIDs = ();
    @SortedPathAtomIDs = sort { $a <=> $b } @PathAtomIDs;

    push @SortedPathAtomIDs, $SortedPathAtomIDs[0];

    $AtomPathID = join '-', @SortedPathAtomIDs;
  }
  else {
    $AtomPathID = join '-', sort { $a <=> $b } @PathAtomIDs;
  }

  if (exists $This->{UniqueCyclicAtomPathsIDs}{$AtomPathID}) {
    return 0;
  }

  # It's a unique atom path containing a cycle...
  $This->{UniqueCyclicAtomPathsIDs}{$AtomPathID} = 1;

  return 1;
}

# Generate atom path string for specified atom path...
#
sub _GenerateAtomPathString {
  my($This, $PathAtomsRef) = @_;
  my($PathLength, $AtomPathString, $ReverseAtomPathString, $AtomPathSymbolsRef);

  $PathLength = scalar @{$PathAtomsRef};

  # Generate path atom and bond symbols...
  #
  $AtomPathSymbolsRef = $This->_GenerateAtomPathSymbols($PathAtomsRef);

  # Check presence of path using path ID created by atom path symbols...
  $AtomPathString = join '', @{$AtomPathSymbolsRef};
  if (exists $This->{AtomPathsStrings}{$PathLength}{$AtomPathString}) {
    $This->{AtomPathsStrings}{$PathLength}{$AtomPathString} += 1;
    return $This;
  }

  # Check presence of reverse path using path ID created by atom path symbols...
  #
  $ReverseAtomPathString = join '', reverse @{$AtomPathSymbolsRef};
  if (exists $This->{AtomPathsStrings}{$PathLength}{$ReverseAtomPathString}) {
    $This->{AtomPathsStrings}{$PathLength}{$ReverseAtomPathString} += 1;
    return $This;
  }

  # Use lexicographically smaller atom path string as PathID...
  #
  if ($AtomPathString le $ReverseAtomPathString) {
    $This->{AtomPathsStrings}{$PathLength}{$AtomPathString} = 1;
  }
  else {
    $This->{AtomPathsStrings}{$PathLength}{$ReverseAtomPathString} = 1;
  }
  return $This;
}

#  Are atom types for all path atoms same?
#
sub _AreAllPathAtomsSymbolsSame {
  my($This, $PathAtomsRef) = @_;
  my($Index, $Atom, $AtomID, $AtomType, $FirstAtomType);

  $Atom = $PathAtomsRef->[0]; $AtomID = $Atom->GetID();
  $FirstAtomType = $This->{AssignedAtomTypes}{$AtomID};

  for $Index (1 .. $#{$PathAtomsRef}) {
    $Atom = $PathAtomsRef->[$Index]; $AtomID = $Atom->GetID();
    $AtomType = $This->{AssignedAtomTypes}{$AtomID};

    if ($AtomType ne $FirstAtomType) {
      return 0;
    }
  }
  return 1;
}

#  Are bond symbols for all path bonds same?
#
sub _AreAllPathBondSymbolsSame {
  my($This, $PathAtomsRef) = @_;
  my($Index, $Atom, $BondedAtom, $AtomID, $BondedAtomID, $BondAtomID1, $BondAtomID2, $FirstBondSymbol, $BondSymbol);

  # During no usage of bond symbols, just ignore them and assume they are same...
  if (!$This->{UseBondSymbols}) {
    return 1;
  }

  $Atom = $PathAtomsRef->[0]; $BondedAtom = $PathAtomsRef->[1];
  $AtomID = $Atom->GetID(); $BondedAtomID = $BondedAtom->GetID();

  ($BondAtomID1, $BondAtomID2) = ($AtomID < $BondedAtomID) ? ($AtomID, $BondedAtomID) : ($BondedAtomID, $AtomID);
  $FirstBondSymbol = $This->{BondSymbols}{$BondAtomID1}{$BondAtomID2};

  for $Index (1 .. ($#{$PathAtomsRef} - 1)) {
    $Atom = $PathAtomsRef->[$Index]; $BondedAtom = $PathAtomsRef->[$Index + 1];
    $AtomID = $Atom->GetID(); $BondedAtomID = $BondedAtom->GetID();

    ($BondAtomID1, $BondAtomID2) = ($AtomID < $BondedAtomID) ? ($AtomID, $BondedAtomID) : ($BondedAtomID, $AtomID);
    $BondSymbol = $This->{BondSymbols}{$BondAtomID1}{$BondAtomID2};

    if ($BondSymbol ne $FirstBondSymbol) {
      return 0;
    }
  }
  return 1;
}

# Generate atom path symbols...
#
sub _GenerateAtomPathSymbols {
  my($This, $PathAtomsRef) = @_;
  my($Atom, $AtomID, @AtomPathSymbols);

  @AtomPathSymbols = ();

  if (@{$PathAtomsRef} == 1) {
    $Atom = $PathAtomsRef->[0]; $AtomID = $Atom->GetID();
    push @AtomPathSymbols, $This->{AssignedAtomTypes}{$AtomID};
    return \@AtomPathSymbols;
  }

  # Ignore bond information...
  if (!$This->{UseBondSymbols}) {
    for $Atom (@{$PathAtomsRef}) {
      $AtomID = $Atom->GetID();
      push @AtomPathSymbols, $This->{AssignedAtomTypes}{$AtomID};
    }
    return \@AtomPathSymbols;
  }

  # Use atoms and bonds to generate atom path string...
  my($Index, $BondedAtom, $BondedAtomID, $BondAtomID1, $BondAtomID2);

  # Process atom type of first atom in path...
  $Atom = $PathAtomsRef->[0]; $AtomID = $Atom->GetID();
  push @AtomPathSymbols, $This->{AssignedAtomTypes}{$AtomID};

  for $Index (0 .. ($#{$PathAtomsRef} - 1)) {
    $Atom = $PathAtomsRef->[$Index]; $BondedAtom = $PathAtomsRef->[$Index + 1];
    $AtomID = $Atom->GetID(); $BondedAtomID = $BondedAtom->GetID();

    ($BondAtomID1, $BondAtomID2) = ($AtomID < $BondedAtomID) ? ($AtomID, $BondedAtomID) : ($BondedAtomID, $AtomID);
    push @AtomPathSymbols, $This->{BondSymbols}{$BondAtomID1}{$BondAtomID2};

    # Process atom type of next atom in path...
    push @AtomPathSymbols, $This->{AssignedAtomTypes}{$BondedAtomID};
  }
  return \@AtomPathSymbols;
}

# Set final fingerprits...
#
sub _SetFinalFingerprints {
  my($This) = @_;

  # Mark successful generation of fingerprints...
  $This->{FingerprintsGenerated} = 1;

  if ($This->{Type} =~ /^PathLengthBits$/i) {
    $This->_SetFinalFingerprintsBitVector();
  }
  elsif ($This->{Type} =~ /^PathLengthCount$/i) {
    $This->_SetFinalFingerprintsVector();
  }

  return $This;
}

# Set final fingerprits bit vector...
#
sub _SetFinalFingerprintsBitVector {
  my($This) = @_;
  my($PathLength, $Size, $AtomPathString, $AtomPathHashCode, $AtomPathBitPos, $FingerprintsBitVector, $SkipBitPosCheck, $NumOfBitsToSetPerPath, $SetBitNum);

  $FingerprintsBitVector = $This->{FingerprintsBitVector};

  $Size = $This->{Size};

  $SkipBitPosCheck = 1;
  $NumOfBitsToSetPerPath = $This->{NumOfBitsToSetPerPath};

  for $PathLength (keys %{$This->{AtomPathsStrings}}) {
    for $AtomPathString (keys %{$This->{AtomPathsStrings}{$PathLength}}) {
      $AtomPathHashCode = TextUtil::HashCode($AtomPathString);

      # Set random number seed...
      if ($This->{UsePerlCoreRandom}) {
	CORE::srand($AtomPathHashCode);
      }
      else {
	MathUtil::srandom($AtomPathHashCode);
      }

      for $SetBitNum (1 .. $NumOfBitsToSetPerPath) {
	$AtomPathBitPos = $This->{UsePerlCoreRandom} ? int(CORE::rand($Size)) : int(MathUtil::random($Size));
	$FingerprintsBitVector->SetBit($AtomPathBitPos, $SkipBitPosCheck);
      }
    }
  }
  return $This;
}

# Set final fingerprits vector...
#
sub _SetFinalFingerprintsVector {
  my($This) = @_;
  my($PathLength, $AtomPathString, $FingerprintsVector, $AtomPathCount, @Values, @ValueIDs);

  @Values = ();
  @ValueIDs = ();

  for $PathLength (sort { $a <=> $b } keys %{$This->{AtomPathsStrings}}) {
    for $AtomPathString (sort keys %{$This->{AtomPathsStrings}{$PathLength}}) {
      $AtomPathCount = $This->{AtomPathsStrings}{$PathLength}{$AtomPathString};

      push @Values, $AtomPathCount;
      push @ValueIDs, $AtomPathString;
    }
  }

  # Add PathLengthIDs and values to fingerprint vector...
  $This->{FingerprintsVector}->AddValueIDs(\@ValueIDs);
  $This->{FingerprintsVector}->AddValues(\@Values);

  return $This;
}

# Cache  appropriate molecule data...
#
sub _SetupMoleculeDataCache {
  my($This) = @_;

  # Get all atoms...
  @{$This->{Atoms}} = $This->GetMolecule()->GetAtoms();

  return $This;
}

# Clear cached molecule data...
#
sub _ClearMoleculeDataCache {
  my($This) = @_;

  # Clear atoms...
  @{$This->{Atoms}} = ();

  # Clear path atoms..
  $This->{AtomPathsRef} = '';

  return $This;
}

# Set atomic invariants to use atom identifiers...
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

# Return a string containg data for PathLengthFingerprints object...
#
sub StringifyPathLengthFingerprints {
  my($This) = @_;
  my($PathLengthsFingerprintsString);

  # Type of fingerprint...
  $PathLengthsFingerprintsString = "Fingerprint type: $This->{Type}; AtomIdentifierType: $This->{AtomIdentifierType}";

  # Path length...
  $PathLengthsFingerprintsString .= "; MinPathLength: $This->{MinLength}; MaxPathLength: $This->{MaxLength}";

  # Fingerprint generation control...
  my($AllowSharedBonds, $AllowRings, $UseBondSymbols, $UseUniquePaths);

  $AllowSharedBonds = $This->{AllowSharedBonds} ? "Yes" : "No";
  $AllowRings = $This->{AllowRings} ? "Yes" : "No";
  $UseBondSymbols = $This->{UseBondSymbols} ? "Yes" : "No";
  $UseUniquePaths = $This->{UseBondSymbols} ? "Yes" : "No";

  $PathLengthsFingerprintsString .= "; UseUniquePaths: $UseUniquePaths; AllowSharedBonds: $AllowSharedBonds; AllowRings: $AllowRings; UseBondSymbols: $UseBondSymbols";

  if ($This->{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
    my($AtomicInvariant, @AtomicInvariants, @AtomicInvariantsOrder, %AvailableAtomicInvariants);

    @AtomicInvariantsOrder = AtomTypes::AtomicInvariantsAtomTypes::GetAtomicInvariantsOrder();
    %AvailableAtomicInvariants = AtomTypes::AtomicInvariantsAtomTypes::GetAvailableAtomicInvariants();

    for $AtomicInvariant (@AtomicInvariantsOrder) {
      push @AtomicInvariants, "$AtomicInvariant: $AvailableAtomicInvariants{$AtomicInvariant}";
    }

    $PathLengthsFingerprintsString .= "; AtomicInvariantsToUse: <" . TextUtil::JoinWords(\@{$This->{AtomicInvariantsToUse}}, ", ", 0) . ">";
    $PathLengthsFingerprintsString .= "; AtomicInvariantsOrder: <" . TextUtil::JoinWords(\@AtomicInvariantsOrder, ", ", 0) . ">";
    $PathLengthsFingerprintsString .= "; AvailableAtomicInvariants: <" . TextUtil::JoinWords(\@AtomicInvariants, ", ", 0) . ">";
  }
  elsif ($This->{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
    my($FunctionalClass, @FunctionalClasses, @FunctionalClassesOrder, %AvailableFunctionalClasses);

    @FunctionalClassesOrder = AtomTypes::FunctionalClassAtomTypes::GetFunctionalClassesOrder();
    %AvailableFunctionalClasses = AtomTypes::FunctionalClassAtomTypes::GetAvailableFunctionalClasses();

    for $FunctionalClass (@FunctionalClassesOrder) {
      push @FunctionalClasses, "$FunctionalClass: $AvailableFunctionalClasses{$FunctionalClass}";
    }

    $PathLengthsFingerprintsString .= "; FunctionalClassesToUse: <" . TextUtil::JoinWords(\@{$This->{FunctionalClassesToUse}}, ", ", 0) . ">";
    $PathLengthsFingerprintsString .= "; FunctionalClassesOrder: <" . TextUtil::JoinWords(\@FunctionalClassesOrder, ", ", 0) . ">";
    $PathLengthsFingerprintsString .= "; AvailableFunctionalClasses: <" . TextUtil::JoinWords(\@FunctionalClasses, ", ", 0) . ">";
  }

  if ($This->{Type} =~ /^PathLengthBits$/i) {
    # Size...
    $PathLengthsFingerprintsString .= "; Size: $This->{Size}; MinSize: $This->{MinSize}; MaxSize: $This->{MaxSize}";

    # NumOfBitsToSetPerPath...
    $PathLengthsFingerprintsString .= "; NumOfBitsToSetPerPath: $This->{NumOfBitsToSetPerPath}";

    # Fingerprint bit density and num of bits set...
    my($NumOfSetBits, $BitDensity);
    $NumOfSetBits = $This->{FingerprintsBitVector}->GetNumOfSetBits();
    $BitDensity = $This->{FingerprintsBitVector}->GetFingerprintsBitDensity();
    $PathLengthsFingerprintsString .= "; NumOfOnBits: $NumOfSetBits; BitDensity: $BitDensity";

    $PathLengthsFingerprintsString .= "; FingerprintsBitVector: < $This->{FingerprintsBitVector} >";
  }
  elsif ($This->{Type} =~ /^PathLengthCount$/i) {
    $PathLengthsFingerprintsString .= "; FingerprintsVector: < $This->{FingerprintsVector} >";
  }

  return $PathLengthsFingerprintsString;
}

1;

__END__

=head1 NAME

PathLengthFingerprints

=head1 SYNOPSIS

use Fingerprints::PathLengthFingerprints;

use Fingerprints::PathLengthFingerprints qw(:all);

=head1 DESCRIPTION

B<PathLengthFingerprints> class provides the following methods:

new, GenerateFingerprints, , GetDescription, SetAtomIdentifierType,
SetAtomicInvariantsToUse, SetFunctionalClassesToUse, SetMaxLength,
SetMinLength, SetNumOfBitsToSetPerPath, SetType,
StringifyPathLengthFingerprints

B<PathLengthFingerprints> is derived from B<Fingerprints> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<PathLengthFingerprints>, B<Fingerprints> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The current release of MayaChemTools supports generation of B<AtomTypesFingerpritns>
corresponding to following B<AtomtomIdentifierTypes>:

    AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
    FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
    SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes

Based on the values specified for B<Type>, B<AtomtomIdentifierTypes>, B<MinPathLength> and
B<MaxPathLength>, all appropriate atom paths are generated for each atom in the molecule
and collected in a list and the list is filtered to remove any structurally duplicate paths as
indicated by the value of B<UseUniquePaths>.

For molecules containing rings, atom paths starting from each atom can be traversed in four
different ways:

    o Atom paths without any rings and sharing of bonds in traversed paths.
    o Atom paths containing rings and without any sharing of bonds in
      traversed paths
    o All possible atom paths without any rings and sharing of bonds in
      traversed paths
    o All possible atom paths containing rings and with sharing of bonds in
      traversed paths.

Atom path traversal is terminated at the last ring atom. For molecules containing no rings,
first two and last two types described above are equivalent.

B<AllowSharedBonds> and B<AllowRings> allow generation of different types of paths
to be used for fingerprints generation.

The combination of B<AllowSharedBonds>, B<AllowRings>, and B<UseBondSymbols> allows generation of
8 different types of path length fingerprints:

    AllowSharedBonds AllowRings UseBondSymbols

    0                0          1   - AtomPathsNoCyclesWithBondSymbols
    0                1          1   - AtomPathsWithCyclesWithBondSymbols

    1                0          1   - AllAtomPathsNoCyclesWithBondSymbols
    1                1          1   - AllAtomPathsWithCyclesWithBondSymbols
                                      [ DEFAULT ]

    0                0          0   - AtomPathsNoCyclesNoBondSymbols
    0                1          0   - AtomPathsWithCyclesNoBondSymbols

    1                0          0   - AllAtomPathsNoCyclesNoBondSymbols
    1                1          0   - AllAtomPathsWithCyclesNoWithBondSymbols

Additionally, possible values for option B<--AtomIdentifierType> in conjunction with corresponding
specified values for B<AtomicInvariantsToUse> and B<FunctionalClassesToUse > changes the nature
of atom path length strings and the fingerprints.

For each atom path in the filtered atom paths list, an atom path string is created using value of
B<AtomIdentifierType> and specified values to use for a particular atom identifier type.
Value of B<UseBondSymbols> controls whether bond order symbols are used during generation
of atom path string. Atom symbol corresponds to element symbol and characters used to represent
 bond order are: I<1 - None; 2 - '='; 3 - '#'; 1.5 or aromatic - ':'; others: bond order value>. By default,
bond symbols are included in atom path strings. Exclusion of bond symbols in atom path strings
results in fingerprints which correspond purely to atom paths without considering bonds.

B<UseUniquePaths> controls the removal of structurally duplicate atom path strings are removed
from the list.

For I<PathLengthBits> value of B<Type>, each atom path is hashed to a 32 bit unsigned
integer key using B<TextUtil::HashCode> function. Using the hash key as a seed for a random number
generator, a random integer value between 0 and B<Size> is used to set corresponding bits
in the fingerprint bit-vector string. Value of B<NumOfBitsToSetPerPaths> option controls the number
of time a random number is generated to set corresponding bits.

For I< PathLengthCount> value of B<Type>n, the number of times an atom path appears
is tracked and a fingerprints count-string corresponding to count of atom paths is generated.

The current release of MayaChemTools generates the following types of path length
fingerprints bit-vector and vector strings:

    FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes:MinLeng
    th1:MaxLength8;1024;BinaryString;Ascending;001000010011010101011000110
    0100010101011000101001011100110001000010001001101000001001001001001000
    0010110100000111001001000001001010100100100000000011000000101001011100
    0010000001000101010100000100111100110111011011011000000010110111001101
    0101100011000000010001000011000010100011101100001000001000100000000...

    FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes:MinLeng
    th1:MaxLength8;1024;HexadecimalString;Ascending;48caa1315d82d91122b029
    42861c9409a4208182d12015509767bd0867653604481a8b1288000056090583603078
    9cedae54e26596889ab121309800900490515224208421502120a0dd9200509723ae89
    00024181b86c0122821d4e4880c38620dab280824b455404009f082003d52c212b4e6d
    6ea05280140069c780290c43

    FingerprintsVector;PathLengthCount:AtomicInvariantsAtomTypes:MinLength
    1:MaxLength8;432;NumericalValues;IDsAndValuesPairsString;C.X1.BO1.H3 2
    C.X2.BO2.H2 4 C.X2.BO3.H1 14 C.X3.BO3.H1 3 C.X3.BO4 10 F.X1.BO1 1 N.X
    2.BO2.H1 1 N.X3.BO3 1 O.X1.BO1.H1 3 O.X1.BO2 2 C.X1.BO1.H3C.X3.BO3.H1
    2 C.X2.BO2.H2C.X2.BO2.H2 1 C.X2.BO2.H2C.X3.BO3.H1 4 C.X2.BO2.H2C.X3.BO
    4 1 C.X2.BO2.H2N.X3.BO3 1 C.X2.BO3.H1:C.X2.BO3.H1 10 C.X2.BO3.H1:C....

    FingerprintsVector;PathLengthCount:DREIDINGAtomTypes:MinLength1:MaxLen
    gth8;410;NumericalValues;IDsAndValuesPairsString;C_2 2 C_3 9 C_R 22 F_
    1 N_3 1 N_R 1 O_2 2 O_3 3 C_2=O_2 2 C_2C_3 1 C_2C_R 1 C_2N_3 1 C_2O_3
    1 C_3C_3 7 C_3C_R 1 C_3N_R 1 C_3O_3 2 C_R:C_R 21 C_R:N_R 2 C_RC_R 2 C
    _RF_ 1 C_RN_3 1 C_2C_3C_3 1 C_2C_R:C_R 2 C_2N_3C_R 1 C_3C_2=O_2 1 C_3C
    _2O_3 1 C_3C_3C_3 5 C_3C_3C_R 2 C_3C_3N_R 1 C_3C_3O_3 4 C_3C_R:C_R ...

    FingerprintsVector;PathLengthCount:EStateAtomTypes:MinLength1:MaxLengt
    h8;454;NumericalValues;IDsAndValuesPairsString;aaCH 14 aasC 8 aasN 1 d
    O 2 dssC 2 sCH3 2 sF 1 sOH 3 ssCH2 4 ssNH 1 sssCH 3 aaCH:aaCH 10 aaCH:
    aasC 8 aasC:aasC 3 aasC:aasN 2 aasCaasC 2 aasCdssC 1 aasCsF 1 aasCssNH
    1 aasCsssCH 1 aasNssCH2 1 dO=dssC 2 dssCsOH 1 dssCssCH2 1 dssCssNH 1
    sCH3sssCH 2 sOHsssCH 2 ssCH2ssCH2 1 ssCH2sssCH 4 aaCH:aaCH:aaCH 6 a...

    FingerprintsVector;PathLengthCount:FunctionalClassAtomTypes:MinLength1
    :MaxLength8;404;NumericalValues;IDsAndValuesPairsString;Ar 22 Ar.HBA 1
    HBA 2 HBA.HBD 3 HBD 1 Hal 1 NI 1 None 10 Ar.HBA:Ar 2 Ar.HBANone 1 Ar:
    Ar 21 ArAr 2 ArHBD 1 ArHal 1 ArNone 2 HBA.HBDNI 1 HBA.HBDNone 2 HBA=NI
    1 HBA=None 1 HBDNone 1 NINone 1 NoneNone 7 Ar.HBA:Ar:Ar 2 Ar.HBA:ArAr
    1 Ar.HBA:ArNone 1 Ar.HBANoneNone 1 Ar:Ar.HBA:Ar 1 Ar:Ar.HBANone 2 ...

    FingerprintsVector;PathLengthCount:MMFF94AtomTypes:MinLength1:MaxLengt
    h8;463;NumericalValues;IDsAndValuesPairsString;C5A 2 C5B 2 C=ON 1 CB 1
    8 COO 1 CR 9 F 1 N5 1 NC=O 1 O=CN 1 O=CO 1 OC=O 1 OR 2 C5A:C5B 2 C5A:N
    5 2 C5ACB 1 C5ACR 1 C5B:C5B 1 C5BC=ON 1 C5BCB 1 C=ON=O=CN 1 C=ONNC=O 1
    CB:CB 18 CBF 1 CBNC=O 1 COO=O=CO 1 COOCR 1 COOOC=O 1 CRCR 7 CRN5 1 CR
    OR 2 C5A:C5B:C5B 2 C5A:C5BC=ON 1 C5A:C5BCB 1 C5A:N5:C5A 1 C5A:N5CR ...

    FingerprintsVector;PathLengthCount:SLogPAtomTypes:MinLength1:MaxLength
    8;518;NumericalValues;IDsAndValuesPairsString;C1 5 C10 1 C11 1 C14 1 C
    18 14 C20 4 C21 2 C22 1 C5 2 CS 2 F 1 N11 1 N4 1 O10 1 O2 3 O9 1 C10C1
    1 C10N11 1 C11C1 2 C11C21 1 C14:C18 2 C14F 1 C18:C18 10 C18:C20 4 C18
    :C22 2 C1C5 1 C1CS 4 C20:C20 1 C20:C21 1 C20:N11 1 C20C20 2 C21:C21 1
    C21:N11 1 C21C5 1 C22N4 1 C5=O10 1 C5=O9 1 C5N4 1 C5O2 1 CSO2 2 C10...

    FingerprintsVector;PathLengthCount:SYBYLAtomTypes:MinLength1:MaxLength
    8;412;NumericalValues;IDsAndValuesPairsString;C.2 2 C.3 9 C.ar 22 F 1
    N.am 1 N.ar 1 O.2 1 O.3 2 O.co2 2 C.2=O.2 1 C.2=O.co2 1 C.2C.3 1 C.2C.
    ar 1 C.2N.am 1 C.2O.co2 1 C.3C.3 7 C.3C.ar 1 C.3N.ar 1 C.3O.3 2 C.ar:C
    .ar 21 C.ar:N.ar 2 C.arC.ar 2 C.arF 1 C.arN.am 1 C.2C.3C.3 1 C.2C.ar:C
    .ar 2 C.2N.amC.ar 1 C.3C.2=O.co2 1 C.3C.2O.co2 1 C.3C.3C.3 5 C.3C.3...

    FingerprintsVector;PathLengthCount:TPSAAtomTypes:MinLength1:MaxLength8
    ;331;NumericalValues;IDsAndValuesPairsString;N21 1 N7 1 None 34 O3 2 O
    4 3 N21:None 2 N21None 1 N7None 2 None:None 21 None=O3 2 NoneNone 13 N
    oneO4 3 N21:None:None 2 N21:NoneNone 2 N21NoneNone 1 N7None:None 2 N7N
    one=O3 1 N7NoneNone 1 None:N21:None 1 None:N21None 2 None:None:None 20
    None:NoneNone 12 NoneN7None 1 NoneNone=O3 2 NoneNoneNone 8 NoneNon...

    FingerprintsVector;PathLengthCount:UFFAtomTypes:MinLength1:MaxLength8;
    410;NumericalValues;IDsAndValuesPairsString;C_2 2 C_3 9 C_R 22 F_ 1 N_
    3 1 N_R 1 O_2 2 O_3 3 C_2=O_2 2 C_2C_3 1 C_2C_R 1 C_2N_3 1 C_2O_3 1 C_
    3C_3 7 C_3C_R 1 C_3N_R 1 C_3O_3 2 C_R:C_R 21 C_R:N_R 2 C_RC_R 2 C_RF_
    1 C_RN_3 1 C_2C_3C_3 1 C_2C_R:C_R 2 C_2N_3C_R 1 C_3C_2=O_2 1 C_3C_2O_3
    1 C_3C_3C_3 5 C_3C_3C_R 2 C_3C_3N_R 1 C_3C_3O_3 4 C_3C_R:C_R 1 C_3...

=head2 METHODS

=over 4

=item B<new>

    $NewPathLengthFingerprints = new PathLengthFingerprints(
                                                   %NamesAndValues);

Using specified I<PathLengthFingerprints> property names and values hash, B<new> method creates a new object
and returns a reference to newly created B<PathLengthFingerprints> object. By default, the following properties are
initialized:

    Molecule = '';
    Type = ''
    Size = 1024
    MinSize = 32
    MaxSize = 2**32
    NumOfBitsToSetPerPath = 1
    MinLength = 1
    MaxLength = 8
    AllowSharedBonds = 1
    AllowRings = 1
    UseBondSymbols = 1
    UseUniquePaths = ''
    AtomIdentifierType = ''
    SetAtomicInvariantsToUse = ['AS']
    FunctionalClassesToUse = ['HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal']

Examples:

    $PathLengthFingerprints = new PathLengthFingerprints(
                              'Molecule' => $Molecule,
                               'Type' => 'PathLengthBits',
                               'AtomIdentifierType' =
                                              'AtomicInvariantsAtomTypes');

    $PathLengthFingerprints = new PathLengthFingerprints(
                               'Molecule' => $Molecule,
                               'Type' => 'PathLengthBits',
                               'Size' => 1024,
                               'MinLength' => 1,
                               'MaxLength' => 8,
                               'AllowRings' => 1,
                               'AllowSharedBonds' => 1,
                               'UseBondSymbols' => 1,
                               'UseUniquePaths' => 1,
                               'AtomIdentifierType' =
                                              'AtomicInvariantsAtomTypes',
                               'AtomicInvariantsToUse' => ['AS']);

    $PathLengthFingerprints = new PathLengthFingerprints(
                               'Molecule' => $Molecule,
                               'Type' => 'PathLengthCount',
                               'MinLength' => 1,
                               'MaxLength' => 8,
                               'AllowRings' => 1,
                               'AllowSharedBonds' => 1,
                               'UseBondSymbols' => 1,
                               'UseUniquePaths' => 1,
                               'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes',
                               'AtomicInvariantsToUse' => ['AS']);

    $PathLengthFingerprints = new PathLengthFingerprints(
                              'Molecule' => $Molecule,
                               'Type' => 'PathLengthBits',
                               'AtomIdentifierType' =
                                              'SLogPAtomTypes');

    $PathLengthFingerprints = new PathLengthFingerprints(
                              'Molecule' => $Molecule,
                               'Type' => 'PathLengthCount',
                               'AtomIdentifierType' =
                                              'SYBYLAtomTypes');

    $PathLengthFingerprints = new PathLengthFingerprints(
                               'Molecule' => $Molecule,
                               'Type' => 'PathLengthBits',
                               'AtomIdentifierType' =
                                              'FunctionalClassAtomTypes',
                               'FunctionalClassesToUse' => ['HBD', 'HBA', 'Ar']);

    $PathLengthFingerprints->GenerateFingerprints();
    print "$PathLengthFingerprints\n";

=item B<GetDescription>

    $Description = $PathLengthFingerprints->GetDescription();

Returns a string containing description of path length fingerprints.

=item B<GenerateFingerprints>

    $PathLengthFingerprints->GenerateFingerprints();

Generates path length fingerprints and returns I<PathLengthFingerprints>.

=item B<SetMaxLength>

    $PathLengthFingerprints->SetMaxLength($Length);

Sets maximum value of atom path length to be used during atom path length fingerprints
generation and returns I<PathLengthFingerprints>

=item B<SetAtomIdentifierType>

    $PathLengthFingerprints->SetAtomIdentifierType();

Sets atom I<IdentifierType> to use during path length fingerprints generation and
returns I<PathLengthFingerprints>.

Possible values: I<AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes,
TPSAAtomTypes, UFFAtomTypes>.

=item B<SetAtomicInvariantsToUse>

    $PathLengthFingerprints->SetAtomicInvariantsToUse($ValuesRef);
    $PathLengthFingerprints->SetAtomicInvariantsToUse(@Values);

Sets atomic invariants to use during I<AtomicInvariantsAtomTypes> value of I<AtomIdentifierType>
for path length fingerprints generation and returns I<PathLengthFingerprints>.

Possible values for atomic invariants are: I<AS, X, BO,  LBO, SB, DB, TB,
H, Ar, RA, FC, MN, SM>. Default value: I<AS>.

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

    $PathLengthFingerprints->SetFunctionalClassesToUse($ValuesRef);
    $PathLengthFingerprints->SetFunctionalClassesToUse(@Values);

Sets functional classes invariants to use during I<FunctionalClassAtomTypes> value of I<AtomIdentifierType>
for path length fingerprints generation and returns I<PathLengthFingerprints>.

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

=item B<SetMinLength>

    $PathLengthFingerprints->SetMinLength($Length);

Sets minimum value of atom path length to be used during atom path length fingerprints
generation and returns I<PathLengthFingerprints>.

=item B<SetMaxLength>

    $PathLengthFingerprints->SetMaxLength($Length);

Sets maximum value of atom path length to be used during atom path length fingerprints
generation and returns I<PathLengthFingerprints>.

=item B<SetNumOfBitsToSetPerPath>

    $PathLengthFingerprints->SetNumOfBitsToSetPerPath($NumOfBits);

Sets number of bits to set for each path during I<PathLengthBits> B<Type > during path length fingerprints
generation and returns I<PathLengthFingerprints>.

=item B<SetType>

    $PathLengthFingerprints->SetType($Type);

Sets type of path length fingerprints and returns I<PathLengthFingerprints>. Possible values:
I<PathLengthBits or PathLengthCount>.

=item B<StringifyPathLengthFingerprints>

    $String = $PathLengthFingerprints->StringifyPathLengthFingerprints();

Returns a string containing information about I<PathLengthFingerprints> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Fingerprints.pm, FingerprintsStringUtil.pm, AtomNeighborhoodsFingerprints.pm,
AtomTypesFingerprints.pm, EStateIndiciesFingerprints.pm, ExtendedConnectivityFingerprints.pm,
MACCSKeys.pm, TopologicalAtomPairsFingerprints.pm, TopologicalAtomTripletsFingerprints.pm,
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
