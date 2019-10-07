package Fingerprints::AtomTypesFingerprints;
#
# File: AtomTypesFingerprints.pm
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
use Text::ParseWords;
use TextUtil ();
use FileUtil ();
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
use overload '""' => 'StringifyAtomTypesFingerprints';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeAtomTypesFingerprints();

  $This->_InitializeAtomTypesFingerprintsProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeAtomTypesFingerprints {
  my($This) = @_;

  # Type of atom type fingerprint to generate:
  #
  # AtomTypesCount - A vector containing count of atom types
  # AtomTypesBits - A bit vector indicating presence/absence of atom types
  #
  $This->{Type} = '';

  # AtomTypes to use for generating fingerprints...
  #
  # Currently supported values are: AtomicInvariantsAtomTypes, DREIDINGAtomTypes,
  # EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
  # SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes
  #
  $This->{AtomIdentifierType} = '';

  # AtomTypesSetToUse for AtomTypesCount:
  #
  # ArbitrarySize - Corrresponds to only AtomTypes atom types detected in molecule
  # FixedSize - Corresponds to fixed number of atom types previously defined for
  #             specific atom types.
  #
  # The default AtomTypesSetToUse value for AtomTypesCount fingerprints type: ArbitrarySize.
  #
  # Possible values: ArbitrarySize or FixedSize. However, for AtomTypesBits fingerprints type, only FixedSize
  # value is allowed.
  #
  $This->{AtomTypesSetToUse} = '';

  # By default, hydrogens are ignored during fingerprint generation...
  $This->{IgnoreHydrogens} = 1;

  # Assigned AtomTypes atom types...
  %{$This->{AtomTypes}} = ();

  # AtomTypes atom types count for generating atom types count and bits fingerprints...
  %{$This->{AtomTypesCount}} = ();
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize object properties....
sub _InitializeAtomTypesFingerprintsProperties {
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

  # Make sure type and identifier type were specified...
  if (!exists $NamesAndValues{Type}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying type...";
  }
  if (!exists $NamesAndValues{AtomIdentifierType}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying AtomIdentifierType...";
  }

  if ($This->{Type} =~ /^AtomTypesCount$/i) {
    $This->_InitializeAtomTypesCount();
  }
  elsif ($This->{Type} =~ /^AtomTypesBits$/i) {
    $This->_InitializeAtomTypesBits();
  }
  else {
    croak "Error: ${ClassName}->_InitializeAtomTypesFingerprintsProperties: Unknown AtomTypes fingerprints type: $This->{Type}; Supported fingerprints types: AtomTypesCount or AtomTypesBits...";
  }

  return $This;
}

# Initialize atom type counts...
#
sub _InitializeAtomTypesCount {
  my($This) = @_;

  # Set default AtomTypesSetToUse...
  if (!$This->{AtomTypesSetToUse}) {
    $This->{AtomTypesSetToUse} = ($This->{AtomIdentifierType} =~ /^TPSAAtomTypes$/i) ? 'FixedSize' : 'ArbitrarySize';
  }

  # Make sure AtomTypesSetToUse value is okay...
  $This->_ValidateAtomTypesSetToUse($This->{AtomTypesSetToUse});

  # Vector type and type of values...
  $This->{VectorType} = 'FingerprintsVector';

  if ($This->{AtomTypesSetToUse} =~ /^FixedSize$/i) {
    $This->{FingerprintsVectorType} = 'OrderedNumericalValues';
  }
  else {
    $This->{FingerprintsVectorType} = 'NumericalValues';
  }

  $This->_InitializeFingerprintsVector();

  return $This;
}

# Initialize atom types bits...
#
sub _InitializeAtomTypesBits {
  my($This) = @_;

  # Set default AtomTypesSetToUse...
  $This->{AtomTypesSetToUse} = 'FixedSize';

  # Make sure AtomTypesSetToUse value is okay...
  $This->_ValidateAtomTypesSetToUse($This->{AtomTypesSetToUse});

  # Vector type...
  $This->{VectorType} = 'FingerprintsBitVector';

  # Vector size...
  $This->{Size} = $This->_GetFixedSizeAtomTypesSetSize();

  $This->_InitializeFingerprintsBitVector();

  return $This;
}

# Set type...
#
sub SetType {
  my($This, $Type) = @_;

  if ($This->{Type}) {
    croak "Error: ${ClassName}->SetType: Can't change type:  It's already set...";
  }

  if ($Type =~ /^AtomTypesCount$/i) {
    $This->{Type} = 'AtomTypesCount';;
  }
  elsif ($Type =~ /^AtomTypesBits$/i) {
    $This->{Type} = 'AtomTypesBits';;
  }
  else {
    croak "Error: ${ClassName}->SetType: Unknown AtomTypes fingerprints type: $Type; Supported fingerprints types: AtomTypesCount or AtomTypesBit...";
  }
  return $This;
}

# Disable set size method...
#
sub SetSize {
  my($This, $Type) = @_;

  croak "Error: ${ClassName}->SetSize: Can't change size:  It's not allowed...";
}

# Set atom types set to use...
#
sub SetAtomTypesSetToUse {
  my($This, $Value) = @_;

  if ($This->{AtomTypesSetToUse}) {
    croak "Error: ${ClassName}->SetAtomTypesSetToUse: Can't change size:  It's already set...";
  }

  $This->_ValidateAtomTypesSetToUse($Value);

  $This->{AtomTypesSetToUse} = $Value;

  return $This;
}

# Validate AtomTypesSetToUse value...
#
sub _ValidateAtomTypesSetToUse {
  my($This, $Value) = @_;

  if ($Value !~ /^(ArbitrarySize|FixedSize)/i) {
    croak "Error: ${ClassName}->_ValidateAtomTypesSetToUse: Unknown AtomTypesSetToUse value: $Value; Supported values: ArbitrarySize or FixedSize";
  }

  if ($Value =~ /^ArbitrarySize$/i && $This->{Type} =~ /^AtomTypesBits$/i) {
    croak "Error: ${ClassName}->_ValidateAtomTypesSetToUse: Specified AtomTypesSetToUse value, $Value, is not allowed for AtomTypesBits fingerprints...";
  }

  if ($Value =~ /^FixedSize$/i && $This->{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
    croak "Error: ${ClassName}->_ValidateAtomTypesSetToUse: Specified AtomTypesSetToUse value, $Value, is not allowed for AtomicInvariantsAtomTypes fingerprints...";
  }

  if ($Value =~ /^FixedSize$/i && $This->{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
    croak "Error: ${ClassName}->_ValidateAtomTypesSetToUse: Specified AtomTypesSetToUse value, $Value, is not allowed for FunctionalClassAtomTypes fingerprints...";
  }

  if ($Value =~ /^ArbitrarySize$/i && $This->{AtomIdentifierType} =~ /^TPSAAtomTypes$/i) {
    croak "Error: ${ClassName}->_ValidateAtomTypesSetToUse: Specified AtomTypesSetToUse value, $Value, is not allowed for TPSAAtomTypes fingerprints...";
  }

  return $This;
}

# Set atom identifier type...
#
sub SetAtomIdentifierType {
  my($This, $IdentifierType) = @_;

  if ($IdentifierType !~ /^(AtomicInvariantsAtomTypes|DREIDINGAtomTypes|EStateAtomTypes|FunctionalClassAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
    croak "Error: ${ClassName}->SetAtomIdentifierType: Specified value, $IdentifierType, for AtomIdentifierType is not vaild. Supported types in current release of MayaChemTools: AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes, and UFFAtomTypes.";
  }

  if ($This->{AtomIdentifierType}) {
    croak "Error: ${ClassName}->SetAtomIdentifierType: Can't change intial atom identifier type:  It's already set...";
  }

  if ($This->{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i && $This->{AtomTypesSetToUse} =~ /^FixedSize$/i) {
    croak "Error: ${ClassName}->SetAtomIdentifierType: Specified AtomTypesSetToUse value, $IdentifierType, is not allowed for AtomicInvariantsAtomTypes fingerprints...";
  }

  if ($This->{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i && $This->{AtomTypesSetToUse} =~ /^FixedSize$/i) {
    croak "Error: ${ClassName}->SetAtomIdentifierType: Specified AtomTypesSetToUse value, $IdentifierType, is not allowed for FunctionalClassAtomTypes fingerprints...";
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

  return "$This->{Type}:$This->{AtomIdentifierType}:$This->{AtomTypesSetToUse}";
}

# Generate atom types fingerprints...
#
# The current release of MayaChemTools supports generation of two types of AtomTypes
# fingerprints corresponding to non-hydrogen and/or hydrogen atoms:
#
# AtomTypesCount - A vector containing count of  atom types
# AtomTypesBits - A bit vector indicating presence/absence of atom types
#
# For AtomTypesCount fingerprints, two types of atom types set size is allowed:
#
# ArbitrarySize - Corrresponds to only atom types detected in molecule
# FixedSize - Corresponds to fixed number of atom types previously defined
#
# For AtomTypesBits fingeprints, only FixedSize atom type set is allowed.
#
# The fixed size atom type set size used during generation of fingerprints corresponding
# to FixedSize value of AtomTypesSetToUse contains all possible atom types in datafiles
# distributed with MayaChemTools release for each supported type.
#
# Combination of Type and AtomTypesSetToUse allow generation of 21 different types of
# AtomTypes fingerprints:
#
# Type                  AtomIdentifierType           AtomTypesSetToUse
#
# AtomTypesCount        AtomicInvariantsAtomTypes    ArbitrarySize
#
# AtomTypesCount        DREIDINGAtomTypes            ArbitrarySize
# AtomTypesCount        DREIDINGAtomTypes            FixedSize
# AtomTypesBits         DREIDINGAtomTypes            FixedSize
#
# AtomTypesCount        EStateAtomTypes              ArbitrarySize
# AtomTypesCount        EStateAtomTypes              FixedSize
# AtomTypesBits         EStateAtomTypes              FixedSize
#
# AtomTypesCount        FunctionalClassAtomTypes    ArbitrarySize
#
# AtomTypesCount        MMFF94AtomTypes              ArbitrarySize
# AtomTypesCount        MMFF94AtomTypes              FixedSize
# AtomTypesBits         MMFF94AtomTypes              FixedSize
#
# AtomTypesCount        SLogPAtomTypes               ArbitrarySize
# AtomTypesCount        SLogPAtomTypes               FixedSize
# AtomTypesBits         SLogPAtomTypes               FixedSize
#
# AtomTypesCount        SYBYLAtomTypes               ArbitrarySize
# AtomTypesCount        SYBYLAtomTypes               FixedSize
# AtomTypesBits         SYBYLAtomTypes               FixedSize
#
# AtomTypesCount        TPSAAtomTypes                 FixedSize
# AtomTypesBits         TPSAAtomTypes                 FixedSize
#
# AtomTypesCount        UFFAtomTypes                 ArbitrarySize
# AtomTypesCount        UFFAtomTypes                 FixedSize
# AtomTypesBits         UFFAtomTypes                 FixedSize
#
sub GenerateFingerprints {
  my($This) = @_;

  # Cache appropriate molecule data...
  $This->_SetupMoleculeDataCache();

  # Check and assign appropriate atom types...
  if (!$This->_AssignAtomTypes()) {
    carp "Warning: ${ClassName}->GenerateFingerprints: $This->{AtomIdentifierType} fingerprints generation didn't succeed: Couldn't assign valid $This->{AtomIdentifierType} to all atoms...";
    return $This;
  }

  # Count atom types...
  $This->_CountAtomTypes();

  # Set final fingerprints...
  $This->_SetFinalFingerprints();

  # Clear cached molecule data...
  $This->_ClearMoleculeDataCache();

  return $This;
}

# Assign appropriate atom types...
#
sub _AssignAtomTypes {
  my($This) = @_;
  my($SpecifiedAtomTypes, $Atom, $AtomID);

  %{$This->{AtomTypes}} = ();
  $SpecifiedAtomTypes = undef;

  IDENTIFIERTYPE: {
    if ($This->{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::AtomicInvariantsAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $This->{IgnoreHydrogens}, 'AtomicInvariantsToUse' => $This->{AtomicInvariantsToUse});
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^DREIDINGAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::DREIDINGAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $This->{IgnoreHydrogens});
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^EStateAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::EStateAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $This->{IgnoreHydrogens});
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::FunctionalClassAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $This->{IgnoreHydrogens}, 'FunctionalClassesToUse' => $This->{FunctionalClassesToUse});
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^MMFF94AtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::MMFF94AtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $This->{IgnoreHydrogens});
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^SLogPAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::SLogPAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $This->{IgnoreHydrogens});
      last IDENTIFIERTYPE;
    }
    if ($This->{AtomIdentifierType} =~ /^SYBYLAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::SYBYLAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $This->{IgnoreHydrogens});
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^TPSAAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::TPSAAtomTypes('Molecule' => $This->{Molecule}, 'IgnorePhosphorus' => 0, 'IgnoreSulfur' => 0);
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^UFFAtomTypes$/i) {
      $SpecifiedAtomTypes = new AtomTypes::UFFAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => $This->{IgnoreHydrogens});
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
  for $Atom (@{$This->{Atoms}}) {
    $AtomID = $Atom->GetID();
    $This->{AtomTypes}{$AtomID} = $SpecifiedAtomTypes->GetAtomType($Atom);
  }

  return $This;
}

# Count atom types...
#
sub _CountAtomTypes {
  my($This) = @_;
  my($Atom, $AtomID, $AtomType);

  %{$This->{AtomTypesCount}} = ();

  for $Atom (@{$This->{Atoms}}) {
    $AtomID = $Atom->GetID();
    $AtomType = $This->{AtomTypes}{$AtomID};

    if (!exists $This->{AtomTypesCount}{$AtomType}) {
      $This->{AtomTypesCount}{$AtomType} = 0;
    }

    $This->{AtomTypesCount}{$AtomType} += 1;
  }
  return $This;
}

# Set final fingerprints...
#
sub _SetFinalFingerprints {
  my($This) = @_;

  # Mark successful generation of fingerprints...
  $This->{FingerprintsGenerated} = 1;

  if ($This->{Type} =~ /^AtomTypesCount$/i) {
    $This->_SetFinalAtomTypesCountFingerprints();
  }
  elsif ($This->{Type} =~ /^AtomTypesBits$/i) {
    $This->_SetFinalAtomTypesBitsFingerprints();
  }
  return $This;
}

# Set final final fingerpritns for atom types count...
#
sub _SetFinalAtomTypesCountFingerprints {
  my($This) = @_;
  my($AtomType, @Values, @IDs);

  @Values = ();
  @IDs = ();

  if ($This->{AtomTypesSetToUse} =~ /^FixedSize$/i) {
    for $AtomType (@{$This->_GetFixedSizeAtomTypesSet()}) {
      push @IDs, $AtomType;
      push @Values, exists($This->{AtomTypesCount}{$AtomType}) ? $This->{AtomTypesCount}{$AtomType} : 0;
    }
  }
  else {
    for $AtomType (sort keys %{$This->{AtomTypesCount}}) {
      push @IDs, $AtomType;
      push @Values, $This->{AtomTypesCount}{$AtomType};
    }
  }

  # Add IDs and values to fingerprint vector...
  if (@IDs) {
    $This->{FingerprintsVector}->AddValueIDs(\@IDs);
  }
  $This->{FingerprintsVector}->AddValues(\@Values);

  return $This;
}

# Set final final fingerpritns for atom types count bits...
#
sub _SetFinalAtomTypesBitsFingerprints {
  my($This) = @_;
  my($AtomType, $SkipPosCheck, $AtomTypeNum, $AtomTypeBitIndex);

  $SkipPosCheck = 1;
  $AtomTypeNum = 0;

  ATOMTYPE: for $AtomType (@{$This->_GetFixedSizeAtomTypesSet()}) {
    $AtomTypeNum++;
    if (!(exists($This->{AtomTypesCount}{$AtomType}) && $This->{AtomTypesCount}{$AtomType})) {
      next ATOMTYPE;
    }
    $AtomTypeBitIndex = $AtomTypeNum - 1;
    $This->{FingerprintsBitVector}->SetBit($AtomTypeBitIndex, $SkipPosCheck);
  }

  return $This;
}

# Cache  appropriate molecule data...
#
sub _SetupMoleculeDataCache {
  my($This) = @_;

  if ($This->{IgnoreHydrogens}) {
    # Get all non-hydrogen atoms...
    my($NegateAtomCheckMethod);
    $NegateAtomCheckMethod = 1;

    @{$This->{Atoms}} = $This->GetMolecule()->GetAtoms("IsHydrogen", $NegateAtomCheckMethod);
  }
  else {
    @{$This->{Atoms}} = $This->GetMolecule()->GetAtoms();
  }

  return $This;
}

# Clear cached molecule data...
#
sub _ClearMoleculeDataCache {
  my($This) = @_;

  @{$This->{Atoms}} = ();

  return $This;
}

# Get fixed size atom types set size...
#
sub _GetFixedSizeAtomTypesSetSize {
  my($This) = @_;
  my($Size);

  $Size = 0;

  IDENTIFIERTYPE: {
    if ($This->{AtomIdentifierType} =~ /^DREIDINGAtomTypes$/i) {
      $Size = $This->{IgnoreHydrogens} ? scalar @{AtomTypes::DREIDINGAtomTypes::GetAllPossibleDREIDINGNonHydrogenAtomTypes()} : scalar @{AtomTypes::DREIDINGAtomTypes::GetAllPossibleDREIDINGAtomTypes()};
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^EStateAtomTypes$/i) {
      $Size = $This->{IgnoreHydrogens} ? scalar @{AtomTypes::EStateAtomTypes::GetAllPossibleEStateNonHydrogenAtomTypes()} : scalar @{AtomTypes::EStateAtomTypes::GetAllPossibleEStateAtomTypes()};
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^MMFF94AtomTypes$/i) {
      $Size = $This->{IgnoreHydrogens} ? scalar @{AtomTypes::MMFF94AtomTypes::GetAllPossibleMMFF94NonHydrogenAtomTypes()} : scalar @{AtomTypes::MMFF94AtomTypes::GetAllPossibleMMFF94AtomTypes()};
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^SLogPAtomTypes$/i) {
      $Size = $This->{IgnoreHydrogens} ? scalar @{AtomTypes::SLogPAtomTypes::GetAllPossibleSLogPNonHydrogenAtomTypes()} : scalar @{AtomTypes::SLogPAtomTypes::GetAllPossibleSLogPAtomTypes()};
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^SYBYLAtomTypes$/i) {
      $Size = $This->{IgnoreHydrogens} ? scalar @{AtomTypes::SYBYLAtomTypes::GetAllPossibleSYBYLNonHydrogenAtomTypes()} : scalar @{AtomTypes::SYBYLAtomTypes::GetAllPossibleSYBYLAtomTypes()};
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^TPSAAtomTypes$/i) {
      $Size =  scalar @{AtomTypes::TPSAAtomTypes::GetAllPossibleTPSAAtomTypes()};
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^UFFAtomTypes$/i) {
      $Size = $This->{IgnoreHydrogens} ? scalar @{AtomTypes::UFFAtomTypes::GetAllPossibleUFFNonHydrogenAtomTypes()} : scalar @{AtomTypes::UFFAtomTypes::GetAllPossibleUFFAtomTypes()};
      last IDENTIFIERTYPE;
    }

    croak "Error: ${ClassName}->_GetFixedSizeAtomTypesSetSize: Atom types set size for atom indentifier type, $This->{AtomIdentifierType}, is not available...";
  }

  return $Size;
}

# Get fixed size atom types set...
#
sub _GetFixedSizeAtomTypesSet {
  my($This) = @_;
  my($AtomTypesRef);

  $AtomTypesRef = undef;

  IDENTIFIERTYPE: {
    if ($This->{AtomIdentifierType} =~ /^DREIDINGAtomTypes$/i) {
      $AtomTypesRef = $This->{IgnoreHydrogens} ? AtomTypes::DREIDINGAtomTypes::GetAllPossibleDREIDINGNonHydrogenAtomTypes() : AtomTypes::DREIDINGAtomTypes::GetAllPossibleDREIDINGAtomTypes();
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^EStateAtomTypes$/i) {
      $AtomTypesRef = $This->{IgnoreHydrogens} ? AtomTypes::EStateAtomTypes::GetAllPossibleEStateNonHydrogenAtomTypes() : AtomTypes::EStateAtomTypes::GetAllPossibleEStateAtomTypes();
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^MMFF94AtomTypes$/i) {
      $AtomTypesRef = $This->{IgnoreHydrogens} ? AtomTypes::MMFF94AtomTypes::GetAllPossibleMMFF94NonHydrogenAtomTypes() : AtomTypes::MMFF94AtomTypes::GetAllPossibleMMFF94AtomTypes();
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^SLogPAtomTypes$/i) {
      $AtomTypesRef = $This->{IgnoreHydrogens} ? AtomTypes::SLogPAtomTypes::GetAllPossibleSLogPNonHydrogenAtomTypes() : AtomTypes::SLogPAtomTypes::GetAllPossibleSLogPAtomTypes();
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^SYBYLAtomTypes$/i) {
      $AtomTypesRef = $This->{IgnoreHydrogens} ? AtomTypes::SYBYLAtomTypes::GetAllPossibleSYBYLNonHydrogenAtomTypes() : AtomTypes::SYBYLAtomTypes::GetAllPossibleSYBYLAtomTypes();
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^TPSAAtomTypes$/i) {
      $AtomTypesRef = AtomTypes::TPSAAtomTypes::GetAllPossibleTPSAAtomTypes();
      last IDENTIFIERTYPE;
    }

    if ($This->{AtomIdentifierType} =~ /^UFFAtomTypes$/i) {
      $AtomTypesRef = $This->{IgnoreHydrogens} ? AtomTypes::UFFAtomTypes::GetAllPossibleUFFNonHydrogenAtomTypes() : AtomTypes::UFFAtomTypes::GetAllPossibleUFFAtomTypes();
      last IDENTIFIERTYPE;
    }

    croak "Error: ${ClassName}->_GetFixedSizeAtomTypesSet: Atom types set for atom indentifier type, $This->{AtomIdentifierType}, is not available...";
  }

  return $AtomTypesRef;
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

# Initialize atomic invariants atom types to use for generating atom IDs in atom pairs...
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
#   AtomID generated by AtomTypes::AtomicInvariantsAtomTypes class corresponds to:
#
#     AS.X<n>.BO<n>.LBO<n>.<SB><n>.<DB><n>.<TB><n>.H<n>.Ar.RA.FC<+n/-n>.MN<n>.SM<n>
#
# Except for AS which is a required atomic invariant atom types AtomIDs, all other atomic invariants are
# optional. Default atomic invariants used for AtomID are: AS, X<n>, BO<n>, H<n>, FC<+n/-n>.
# AtomID specification doesn't include atomic invariants with zero or undefined values.
#
# Examples of  AtomIDs:
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

  # Default atomic invariants to use for generating atom pair atom IDs: AS, X, BO, H, FC
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

# Set atomic invariants to use for atom IDs...
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

# Return a string containg data for AtomTypesFingerprints object...
sub StringifyAtomTypesFingerprints {
  my($This) = @_;
  my($FingerprintsString, $IgnoreHydrogens);

  $FingerprintsString = "Type: $This->{Type}; AtomIdentifierType: $This->{AtomIdentifierType}; AtomTypesSetToUse: $This->{AtomTypesSetToUse}";

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


  $IgnoreHydrogens = $This->{IgnoreHydrogens} ? "Yes" : "No";
  $FingerprintsString .= "; IgnoreHydrogens: $IgnoreHydrogens";

  if ($This->{Type} =~ /^AtomTypesCount$/i) {
    $FingerprintsString .= "; FingerprintsVector: < $This->{FingerprintsVector} >";
  }
  elsif ($This->{Type} =~ /^AtomTypesBits$/i) {
    $FingerprintsString .= "; FingerprintsBitVector: < $This->{FingerprintsBitVector} >";
  }

  return $FingerprintsString;
}

1;

__END__

=head1 NAME

AtomTypesFingerprints

=head1 SYNOPSIS

use Fingerprints::AtomTypesFingerprints;

use Fingerprints::AtomTypesFingerprints qw(:all);

=head1 DESCRIPTION

B<AtomTypesFingerprints> class provides the following methods:

new, GenerateFingerprints, GetDescription, SetAtomIdentifierType,
SetAtomTypesSetToUse, SetAtomicInvariantsToUse, SetFunctionalClassesToUse,
SetType, StringifyAtomTypesFingerprints

B<AtomTypesFingerprints> is derived from B<Fingerprints> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<AtomNeighborhoodsFingerprints>, B<Fingerprints> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The current release of MayaChemTools supports generation of B<AtomTypesFingerpritns>
corresponding to following B<AtomIdentifierTypes>:

    AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
    FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
    SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes

Based on the values specified for B<AtomIdentifierType> along with other specified
parameters such as B<AtomicInvariantsToUse> and B<FunctionalClassesToUse>, initial
atom types are assigned to all non-hydrogen atoms or all atoms in a molecule.

Using the assigned atom types and specified B<Type>, one of the following types of
fingerprints are generated:

    AtomTypesCount - A vector containing count of atom types
    AtomTypesBits - A bit vector indicating presence/absence of atom types

For I<AtomTypesCount> fingerprints, two types of atom types set size is allowed:

    ArbitrarySize - Corresponds to only atom types detected in molecule
    FixedSize - Corresponds to fixed number of atom types previously defined

For I<AtomTypesBits> fingerprints, only I<FixedSize> atom type set is allowed.

I<ArbitrarySize> corresponds to atom types detected in a molecule where as I<FixedSize> implies
a fix number of all possible atom types previously defined for a specific I<AtomIdentifierType>.

Fix number of all possible atom types for supported I<AtomIdentifierTypes> in current release
of MayaChemTools are:

    AtomIdentifier       Total    TotalWithoutHydrogens

    DREIDINGAtomTypes    37       34
    EStateAtomTypes      109      87
    MMFF94AtomTypes      212      171
    SLogPAtomTypes       72       67
    SYBYLAtomTypes       45       44
    TPSAAtomTypes        47       47
    UFFAtomTypes         126      124

Combination of B<Type> and B<AtomTypesSetToUse> along with B<AtomtomIdentifierType>
allows generation of following different atom types fingerprints:

    Type                  AtomIdentifierType           AtomTypesSetToUse

    AtomTypesCount        AtomicInvariantsAtomTypes    ArbitrarySize

    AtomTypesCount        DREIDINGAtomTypes            ArbitrarySize
    AtomTypesCount        DREIDINGAtomTypes            FixedSize
    AtomTypesBits         DREIDINGAtomTypes            FixedSize

    AtomTypesCount        EStateAtomTypes              ArbitrarySize
    AtomTypesCount        EStateAtomTypes              FixedSize
    AtomTypesBits         EStateAtomTypes              FixedSize

    AtomTypesCount        FunctionalClassAtomTypes     ArbitrarySize

    AtomTypesCount        MMFF94AtomTypes              ArbitrarySize
    AtomTypesCount        MMFF94AtomTypes              FixedSize
    AtomTypesBits         MMFF94AtomTypes              FixedSize

    AtomTypesCount        SLogPAtomTypes               ArbitrarySize
    AtomTypesCount        SLogPAtomTypes               FixedSize
    AtomTypesBits         SLogPAtomTypes               FixedSize

    AtomTypesCount        SYBYLAtomTypes               ArbitrarySize
    AtomTypesCount        SYBYLAtomTypes               FixedSize
    AtomTypesBits         SYBYLAtomTypes               FixedSize

    AtomTypesCount        TPSAAtomTypes                FixedSize
    AtomTypesBits         TPSAAtomTypes                FixedSize

    AtomTypesCount        UFFAtomTypes                 ArbitrarySize
    AtomTypesCount        UFFAtomTypes                 FixedSize
    AtomTypesBits         UFFAtomTypes                 FixedSize

The current release of MayaChemTools generates the following types of atom types
fingerprints bit-vector and vector strings:

    FingerprintsVector;AtomTypesCount:AtomicInvariantsAtomTypes:ArbitraryS
    ize;10;NumericalValues;IDsAndValuesString;C.X1.BO1.H3 C.X2.BO2.H2 C.X2
    .BO3.H1 C.X3.BO3.H1 C.X3.BO4 F.X1.BO1 N.X2.BO2.H1 N.X3.BO3 O.X1.BO1.H1
    O.X1.BO2;2 4 14 3 10 1 1 1 3 2

    FingerprintsVector;AtomTypesCount:DREIDINGAtomTypes:ArbitrarySize;8;Nu
    mericalValues;IDsAndValuesString;C_2 C_3 C_R F_ N_3 N_R O_2 O_3;2 9 22
    1 1 1 2 3

    FingerprintsVector;AtomTypesCount:DREIDINGAtomTypes:FixedSize;34;Order
    edNumericalValues;IDsAndValuesString;B_3 B_2 C_3 C_R C_2 C_1 N_3 N_R N
    _2 N_1 O_3 O_R O_2 O_1 F_ Al3 Si3 P_3 S_3 Cl Ga3 Ge3 As3 Se3 Br In3 Sn
    3 Sb3 Te3 I_ Na Ca Fe Zn;0 0 9 22 2 0 1 1 0 0 3 0 2 0 1 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0

    FingerprintsBitVector;AtomTypesBits:DREIDINGAtomTypes:FixedSize;34;Bin
    aryString;Ascending;0011101100101010000000000000000000000000

    FingerprintsVector;AtomTypesCount:EStateAtomTypes:ArbitrarySize;11;Num
    ericalValues;IDsAndValuesString;aaCH aasC aasN dO dssC sCH3 sF sOH ssC
    H2 ssNH sssCH;14 8 1 2 2 2 1 3 4 1 3

    FingerprintsVector;AtomTypesCount:EStateAtomTypes:FixedSize;87;Ordered
    NumericalValues;IDsAndValuesString;sLi ssBe ssssBem sBH2 ssBH sssB sss
    sBm sCH3 dCH2 ssCH2 tCH dsCH aaCH sssCH ddC tsC dssC aasC aaaC ssssC s
    NH3p sNH2 ssNH2p dNH ssNH aaNH tN sssNHp dsN aaN sssN ddsN aasN ss...;
    0 0 0 0 0 0 0 2 0 4 0 0 14 3 0 0 2 8 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 3 2 0 0
    0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0...

    FingerprintsBitVector;AtomTypesBits:EStateAtomTypes:FixedSize;87;Binar
    yString;Ascending;0000000101001100110000001000000010110000100000000000
    000000000000000000000000000000000000

    FingerprintsVector;AtomTypesCount:FunctionalClassAtomTypes:ArbitrarySi
    ze;8;NumericalValues;IDsAndValuesString;Ar Ar.HBA HBA HBA.HBD HBD Hal 
    NI None;22 1 2 3 1 1 1 10

    FingerprintsVector;AtomTypesCount:MMFF94AtomTypes:ArbitrarySize;13;Num
    ericalValues;IDsAndValuesString;C5A C5B C=ON CB COO CR F N5 NC=O O=CN
    O=CO OC=O OR;2 2 1 18 1 9 1 1 1 1 1 1 2

    FingerprintsVector;AtomTypesCount:MMFF94AtomTypes:FixedSize;171;Ordere
    dNumericalValues;IDsAndValuesString;CR C=C CSP2 C=O C=N CGD C=OR C=ON
    CONN COO COON COOO C=OS C=S C=SN CSO2 CS=O CSS C=P CSP =C= OR OC=O OC=
    C OC=N OC=S ONO2 ON=O OSO3 OSO2 OSO OS=O -OS OPO3 OPO2 OPO -OP -O-...;
    9 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 2 1 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 1 0 1 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 0 0 0 0 0 0 0 0 ...

    FingerprintsBitVector;AtomTypesBits:MMFF94AtomTypes:FixedSize;171;Bina
    ryString;Ascending;100000010100000000000110000000000000000101000000100
    0100000000000000000000000000000000000000000100000000000000000000000000
    0000000011000000000000000001000000000000000000000000000

    FingerprintsVector;AtomTypesCount:SLogPAtomTypes:ArbitrarySize;16;Nume
    ricalValues;IDsAndValuesString;C1 C10 C11 C14 C18 C20 C21 C22 C5 CS F
    N11 N4 O10 O2 O9;5 1 1 1 14 4 2 1 2 2 1 1 1 1 3 1

    FingerprintsVector;AtomTypesCount:SLogPAtomTypes:FixedSize;67;OrderedN
    umericalValues;IDsAndValuesString;C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C
    12 C13 C14 C15 C16 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C27 CS N1 N
    2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12 N13 N14 NS O1 O2 O3 O4 O5 O6 O7 O8
    O9 O10 O11 O12 OS F Cl Br I Hal P S1 S2 S3 Me1 Me2;5 0 0 0 2 0 0 0 0 1
    1 0 0 1 0 0 0 14 0 4 2 1 0 0 0 0 0 2 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0...

    FingerprintsBitVector;AtomTypesBits:SLogPAtomTypes:FixedSize;67;Binary
    String;Ascending;10001000011001000101110000010001000000100000100000011
    0001000000000000000

    FingerprintsVector;AtomTypesCount:SYBYLAtomTypes:ArbitrarySize;9;Numer
    icalValues;IDsAndValuesString;C.2 C.3 C.ar F N.am N.ar O.2 O.3 O.co2;2
    9 22 1 1 1 1 2 2

    FingerprintsVector;AtomTypesCount:SYBYLAtomTypes:FixedSize;44;OrderedN
    umericalValues;IDsAndValuesString;C.3 C.2 C.1 C.ar C.cat N.3 N.2 N.1 N
    .ar N.am N.pl3 N.4 O.3 O.2 O.co2 S.3 S.2 S.o S.o2 P.3 F Cl Br I ANY HA
    L HET Li Na Mg Al Si K Ca Cr.th Cr.oh Mn Fe Co.oh Cu Zn Se Mo Sn;9 2 0
    22 0 0 0 0 1 1 0 0 2 1 2 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0

    FingerprintsBitVector;AtomTypesBits:SYBYLAtomTypes:FixedSize;44;Binary
    String;Ascending;110100001100111000001000000000000000000000000000

    FingerprintsVector;AtomTypesCount:TPSAAtomTypes:FixedSize;47;OrderedNu
    mericalValues;IDsAndValuesString;N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N1
    2 N13 N14 N15 N16 N17 N18 N19 N20 N21 N22 N23 N24 N25 N26 N O1 O2 O3 O
    4 O5 O6 O S1 S2 S3 S4 S5 S6 S7 S P1 P2 P3 P4 P;0 0 0 0 0 0 1 0 0 0 0 0
    0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

    FingerprintsBitVector;AtomTypesBits:TPSAAtomTypes:FixedSize;47;BinaryS
    tring;Ascending;000000100000000000001000000001100000000000000000

    FingerprintsVector;AtomTypesCount:UFFAtomTypes:ArbitrarySize;8;Numeric
    alValues;IDsAndValuesString;C_2 C_3 C_R F_ N_3 N_R O_2 O_3;2 9 22 1 1
    1 2 3

    FingerprintsVector;AtomTypesCount:UFFAtomTypes;124;OrderedNumerical
    Values;IDsAndValuesString;He4+4 Li Be3+2 B_3 B_2 C_3 C_R C_2 C_1 N_3 N_
    R N_2 N_1 O_3 O_3_z O_R O_2 O_1 F_ Ne4+4 Na Mg3+2 Al3 Si3 P_3+3 P_3+5 P
    _3+q S_3+2 S_3+4 S_3+6 S_R S_2 Cl Ar4+4 K_ Ca6+2 Sc3+3 Ti3+4 Ti6+4 V_3+
    ;0 0 0 0 0 12 0 3 0 3 0 1 0 2 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...

    FingerprintsVector;AtomTypesCount:UFFAtomTypes:FixedSize;124;OrderedNu
    mericalValues;IDsAndValuesString;He4+4 Li Be3+2 B_3 B_2 C_3 C_R C_2 C_
    1 N_3 N_R N_2 N_1 O_3 O_3_z O_R O_2 O_1 F_ Ne4+4 Na Mg3+2 Al3 Si3 P_3+
    3 P_3+5 P_3+q S_3+2 S_3+4 S_3+6 S_R S_2 Cl Ar4+4 K_ Ca6+2 Sc3+3 Ti...;
    0 0 0 0 0 9 22 2 0 1 1 0 0 3 0 0 2 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0...

    FingerprintsBitVector;AtomTypesBits:UFFAtomTypes:FixedSize;124;BinaryS
    tring;Ascending;000001110110010010100000000000000000000000000000000000
    0000000000000000000000000000000000000000000000000000000000000000000000

=head2 METHODS

=over 4

=item B<new>

    $NewAtomTypesFingerprints = new AtomTypesFingerprints(%NamesAndValues);

Using specified I<AtomTypesFingerprints> property names and values hash, B<new> method creates a new object
and returns a reference to newly created B<PathLengthFingerprints> object. By default, the
following properties are initialized:

    Molecule = '';
    Type = ''
    AtomIdentifierType = ''
    AtomTypesSetToUse = ''
    IgnoreHydrogens = 1
    AtomicInvariantsToUse = ['AS', 'X', 'BO', 'H', 'FC', 'MN']
    FunctionalClassesToUse = ['HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal']

Examples:

    $AtomTypesFingerprints = new AtomTypesFingerprints(
                              'Molecule' => $Molecule,
                              'Type' => 'AtomTypesCount',
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes');

    $AtomTypesFingerprints = new AtomTypesFingerprints(
                              'Molecule' => $Molecule,
                              'Type' => 'AtomTypesCount',
                              'AtomIdentifierType' =>
                                              'AtomicInvariantsAtomTypes',
                              'AtomicInvariantsToUse' =>
                                              ['AS', 'X', 'BO', 'H', 'FC'] );

    $AtomTypesFingerprints = new AtomTypesFingerprints(
                              'Molecule' => $Molecule,
                              'Type' => 'AtomTypesCount',
                              'AtomIdentifierType' =>
                                              'DREIDINGAtomTypes');

    $AtomTypesFingerprints = new AtomTypesFingerprints(
                              'Molecule' => $Molecule,
                              'Type' => 'AtomTypesCount',
                              'AtomIdentifierType' =>
                                              'EStateAtomTypes',
                              'AtomTypesSetToUse' =>
                                              'ArbitrarySize');

    $AtomTypesFingerprints = new AtomTypesFingerprints(
                              'Molecule' => $Molecule,
                              'Type' => 'AtomTypesCount',
                              'AtomIdentifierType' =>
                                              'SLogPAtomTypes',
                              'AtomTypesSetToUse' =>
                                              'FixedSize');

    $AtomTypesFingerprints = new AtomTypesFingerprints(
                              'Molecule' => $Molecule,
                              'Type' => 'AtomTypesBits',
                              'AtomIdentifierType' =>
                                              'MMFF94AtomTypes',
                              'AtomTypesSetToUse' =>
                                              'FixedSize');

    $AtomTypesFingerprints->GenerateFingerprints();
    print "$AtomTypesFingerprints\n";

=item B<GenerateFingerprints>

    $AtomTypesFingerprints->GenerateFingerprints();

Generates atom types fingerprints and returns I<AtomTypesFingerprints>.

=item B<GetDescription>

    $Description = $AtomTypesFingerprints->GetDescription();

Returns a string containing description of atom types fingerprints.

=item B<SetAtomIdentifierType>

    $AtomTypesFingerprints->SetAtomIdentifierType($IdentifierType);

Sets atom I<IdentifierType> to use during atom types fingerprints generation and
returns I<AtomTypesFingerprints>.

Possible values: I<AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes,
TPSAAtomTypes, UFFAtomTypes>.

=item B<SetAtomTypesSetToUse>

    $AtomTypesFingerprints->SetAtomTypesSetToUse($Value);

Sets I<Value> of I<AtomTypesSetToUse> and returns I<AtomTypesFingerprints>. Possible
values: I<ArbitrarySize or FixedSize>. Default for I<AtomTypesCount> value of
B<AtomTypesSetToUse>: I<ArbitrarySize>.

=item B<SetAtomicInvariantsToUse>

    $AtomTypesFingerprints->SetAtomicInvariantsToUse($ValuesRef);
    $AtomTypesFingerprints->SetAtomicInvariantsToUse(@Values);

Sets atomic invariants to use during I<AtomicInvariantsAtomTypes> value of I<AtomIdentifierType>
for atom neighborhood fingerprints generation and returns I<AtomTypesFingerprints>.

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

    $AtomTypesFingerprints->SetFunctionalClassesToUse($ValuesRef);
    $AtomTypesFingerprints->SetFunctionalClassesToUse(@Values);

Sets functional classes invariants to use during I<FunctionalClassAtomTypes> value of I<AtomIdentifierType>
for atom types fingerprints generation and returns I<AtomTypesFingerprints>.

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

=item B<SetType>

    $AtomTypesFingerprints->SetType($Type);

Sets type of AtomTypes fingerpritns and returns I<AtomTypesFingerprints>. Possible values: I<AtomTypesFingerprintsBits or
AtomTypesFingerprintsCount>.


=item B<StringifyAtomTypesFingerprints>

    $String = $AtomTypesFingerprints->StringifyAtomTypesFingerprints();

Returns a string containing information about I<AtomTypesFingerprints> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Fingerprints.pm, FingerprintsStringUtil.pm, AtomNeighborhoodsFingerprints.pm,
EStateIndiciesFingerprints.pm, ExtendedConnectivityFingerprints.pm, MACCSKeys.pm,
PathLengthFingerprints.pm, TopologicalAtomPairsFingerprints.pm, TopologicalAtomTripletsFingerprints.pm,
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
