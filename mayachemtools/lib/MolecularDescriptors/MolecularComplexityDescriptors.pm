package MolecularDescriptors::MolecularComplexityDescriptors;
#
# File: MolecularComplexityDescriptors.pm
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
use TextUtil ();
use MathUtil ();
use Atom;
use Molecule;
use MolecularDescriptors::MolecularDescriptors;
use AtomTypes::AtomicInvariantsAtomTypes;
use AtomTypes::FunctionalClassAtomTypes;
use Fingerprints::AtomTypesFingerprints;
use Fingerprints::ExtendedConnectivityFingerprints;
use Fingerprints::MACCSKeys;
use Fingerprints::PathLengthFingerprints;
use Fingerprints::TopologicalAtomPairsFingerprints;
use Fingerprints::TopologicalAtomTripletsFingerprints;
use Fingerprints::TopologicalAtomTorsionsFingerprints;
use Fingerprints::TopologicalPharmacophoreAtomPairsFingerprints;
use Fingerprints::TopologicalPharmacophoreAtomTripletsFingerprints;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(MolecularDescriptors::MolecularDescriptors Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GetDescriptorNames GetMolecularComplexityTypeAbbreviation);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, @DescriptorNames);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyMolecularComplexityDescriptors';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeMolecularComplexityDescriptors();

  $This->_InitializeMolecularComplexityDescriptorsProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Descriptor names...
  @DescriptorNames = ('MolecularComplexity');

}

# Get descriptor names as an array.
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetDescriptorNames {
  return @DescriptorNames;
}

# Initialize object data...
#
sub _InitializeMolecularComplexityDescriptors {
  my($This) = @_;

  # Type of MolecularDescriptor...
  $This->{Type} = 'MolecularComplexity';

  #
  # The current release of MayaChemTools supports calculation of molecular complexity
  # corresponding to number of bits-set or unique keys [ Ref 117-119 ] in molecular
  # fingerprints. The following types of fingerprints based molecular complexity measures
  # are supported:
  #
  # AtomTypesFingerprints
  # ExtendedConnectivityFingerprints
  # MACCSKeys
  # PathLengthFingerprints
  # TopologicalAtomPairsFingerprints
  # TopologicalAtomTripletsFingerprints
  # TopologicalAtomTorsionsFingerprints
  # TopologicalPharmacophoreAtomPairsFingerprints
  # TopologicalPharmacophoreAtomTripletsFingerprints
  #
  # Default: MACCSKeys
  #
  $This->{MolecularComplexityType} = '';

  # Atom types to use for generating fingerprints...
  #
  # Currently supported values are: AtomicInvariantsAtomTypes, DREIDINGAtomTypes,
  # EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes,
  # SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes
  #
  # Notes:
  #   . AtomicInvariantsAtomTypes for all supported MolecularComplexityType except for
  #     TopologicalPharmacophoreAtomPairsFingerprints and TopologicalPharmacophoreAtomTripletsFingerprints
  #   . This value is not used for MACCSKeys
  #   . FunctionalClassAtomTypes is the only valid value during topological pharmacophore fingerprints.
  #
  #   . Default values for AtomicInvariantsToUse and FunctionalClassesToUse are set appropriately
  #     for different types of fingerprints as shown below.
  #
  #     MolecularComplexityType              AtomicInvariantsToUse
  #
  #     AtomTypesFingerprints                AS, X, BO, H, FC
  #     TopologicalAtomPairsFingerprints     AS, X, BO, H, FC
  #     TopologicalAtomTripletsFingerprints  AS, X, BO, H, FC
  #     TopologicalAtomTorsionsFingerprints  AS, X, BO, H, FC
  #
  #     ExtendedConnectivityFingerprints     AS, X, BO, H, FC, MN
  #     PathLengthFingerprints               AS
  #
  #     Default for FunctionalClassesToUse for all fingerprints is set to:
  #
  #     HBD, HBA, PI, NI, Ar, Hal
  #
  #     except for the following two MolecularComplexityType fingerprints:
  #
  #     TopologicalPharmacophoreAtomPairsFingerprints     HBD, HBA, PI, NI, H
  #     TopologicalPharmacophoreAtomTripletsFingerprints  HBD, HBA, PI, NI, H, Ar
  #
  $This->{AtomIdentifierType} = '';

  # Size of MACCS key set: 166 or 322...
  #
  $This->{MACCSKeysSize} = 166;

  # Atomic neighborhoods radius for extended connectivity fingerprints...
  $This->{NeighborhoodRadius} = 2;

  # Minimum and maximum path lengths to use for path length fingerprints...
  $This->{MinPathLength} = 1;
  $This->{MaxPathLength} = 8;

  # By default bond symbols are included in atom path strings used to generate path length
  # fingerprints... ...
  $This->{UseBondSymbols} = 1;

  # Minimum and maximum bond distance between atom pairs during topological
  # atom pairs/triplets fingerprints...
  $This->{MinDistance} = 1;
  $This->{MaxDistance} = 10;

  # Determines whether to apply triangle inequality to distance triplets...
  #
  # Default for TopologicalAtomTripletsFingerprints: 0
  # Default for TopologicalPharmacophoreAtomTripletsFingerprints: 1
  #
  $This->{UseTriangleInequality} = '';

  # Distance bin size used for binning distances during generation of
  # topological pharmacophore atom triplets fingerprints...
  #
  $This->{DistanceBinSize} = 2;

  # Normalization methodology to use for scaling the number of bits-set or unique keys
  # for:
  #
  # ExtendedConnectivityFingerprints
  # TopologicalPharmacophoreAtomPairsFingerprints
  # TopologicalPharmacophoreAtomTripletsFingerprints
  #
  # This option is gnored for all other types of fingerprints.
  #
  # Possible values during extended connectivity fingerprints: None or ByHeavyAtomsCount. Default:
  # None.
  #
  # Possible values during topological pharmacophore atom pairs and tripletes fingerprints: None,
  # or ByPossibleKeysCount. Default: None. ByPossibleKeysCount corresponds to total number of
  # possible topological pharmacophore atom pairs or triplets in a molecule.
  #
  #
  $This->{NormalizationMethodology} = 'None';

  # Intialize descriptor names and values...
  $This->_InitializeDescriptorNamesAndValues(@DescriptorNames);

  return $This;
}

# Initialize object properties...
#
sub _InitializeMolecularComplexityDescriptorsProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  # Make sure MolecularComplexityType is set...
  if (!exists $NamesAndValues{MolecularComplexityType}) {
    $This->{MolecularComplexityType} = 'MACCSKeys';
  }

  # Make sure AtomIdentifierType is set...
  if ($This->{MolecularComplexityType} !~ /^MACCSKeys$/i) {
    if (!exists $NamesAndValues{AtomIdentifierType}) {
      $This->_InitializeAtomIdentifierType();
    }
  }

  # Make sure UseTriangleInequality is set...
  if ($This->{MolecularComplexityType} =~ /^(TopologicalAtomTripletsFingerprints|TopologicalPharmacophoreAtomTripletsFingerprints)$/i) {
    if (!exists $NamesAndValues{UseTriangleInequality}) {
      $This->{UseTriangleInequality} =  ($This->{MolecularComplexityType} =~ /^TopologicalPharmacophoreAtomTripletsFingerprints$/i) ? 1 : 0;
    }
  }

  return $This;
}

# Initialize atom identifer type...
#
sub _InitializeAtomIdentifierType {
  my($This) = @_;
  my($AtomIdentifierType);

  if ($This->{MolecularComplexityType} =~ /^MACCSKeys$/i) {
    return $This;
  }

  $AtomIdentifierType = ($This->{MolecularComplexityType} =~ /^(TopologicalPharmacophoreAtomPairsFingerprints|TopologicalPharmacophoreAtomTripletsFingerprints)$/i) ? 'FunctionalClassAtomTypes' : 'AtomicInvariantsAtomTypes';

  $This->SetAtomIdentifierType($AtomIdentifierType);

  return $This;
}

# Get abbreviation for specified molecular complexity type or using descriptors object...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetMolecularComplexityTypeAbbreviation {
  my($FirstParameter) = @_;
  my($This, $ComplexityType, %ComplexityTypeToAbbrev);

  if (_IsMolecularComplexityDescriptors($FirstParameter)) {
    $This = $FirstParameter;
    $ComplexityType = $This->{MolecularComplexityType};
  }
  else {
    $ComplexityType = $FirstParameter;
  }

  %ComplexityTypeToAbbrev = (lc 'AtomTypesFingerprints' => 'ATFP', lc 'ExtendedConnectivityFingerprints' => 'ECFP',
			     lc 'MACCSKeys' => 'MACCSKeys', lc 'PathLengthFingerprints' => 'PLFP',
			     lc 'TopologicalAtomPairsFingerprints' => 'TAPFP', lc 'TopologicalAtomTripletsFingerprints' => 'TATFP',
			     lc 'TopologicalAtomTorsionsFingerprints' => 'TATFP',
			     lc 'TopologicalPharmacophoreAtomPairsFingerprints' => 'TPAPFP',
			     lc 'TopologicalPharmacophoreAtomTripletsFingerprints' => 'TPATFP');

  return exists $ComplexityTypeToAbbrev{lc $ComplexityType} ? $ComplexityTypeToAbbrev{lc $ComplexityType} : '';
}

# Set MACCS key set size...
#
sub SetMACCSKeysSize {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetMACCSKeysSize: Size value, $Value, is not valid:  It must be a positive integer...";
  }
  if ($Value !~ /^(166|322)/i) {
    croak "Error: ${ClassName}->SetMACCSKeysSize: The current release of MayaChemTools doesn't support MDL MACCS $Value keys...";
  }
  $This->{MACCSKeysSize} = $Value;

  return $This;
}

# Set minimum path length...
#
sub SetMinPathLength {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetMinPathLength: MinPathLength value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{MinPathLength} = $Value;

  return $This;
}

# Set maximum path length...
#
sub SetMaxPathLength {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetMaxPathLength: MaxPathLength value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{MaxPathLength} = $Value;

  return $This;
}

# Set minimum  bond distance between atom pairs during topological and topological
# pharmacophore atom pairs/triplets fingerprints...
#
sub SetMinDistance {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetMinDistance: MinDistance value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{MinDistance} = $Value;

  return $This;
}

# Set maximum  bond distance between atom pairs during topological and topological
# pharmacophore atom pairs/triplets fingerprints...
#
sub SetMaxDistance {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetMaxDistance: MaxDistance value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{MaxDistance} = $Value;

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

# Set molecular complexity type...
#
sub SetMolecularComplexityType {
  my($This, $Value) = @_;

  if ($Value !~ /^(AtomTypesFingerprints|ExtendedConnectivityFingerprints|MACCSKeys|PathLengthFingerprints|TopologicalAtomPairsFingerprints|TopologicalAtomTripletsFingerprints|TopologicalAtomTorsionsFingerprints|TopologicalPharmacophoreAtomPairsFingerprints|TopologicalPharmacophoreAtomTripletsFingerprints)$/i) {
    croak "Error: ${ClassName}->SetMolecularComplexityType: MolecularComplexityType value, $Value, is not valid. Supported values: AtomTypesFingerprints, ExtendedConnectivityFingerprints, MACCSKeys, PathLengthFingerprints, TopologicalAtomPairsFingerprints, TopologicalAtomTripletsFingerprints, TopologicalAtomTorsionsFingerprints, TopologicalPharmacophoreAtomPairsFingerprints, or TopologicalPharmacophoreAtomTripletsFingerprints...";
  }

  $This->{MolecularComplexityType} = $Value;

  return $This;
}

# Set distance bin size for binning pharmacophore atom pair distances in atom triplets...
#
sub SetDistanceBinSize {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetDistanceBinSize: DistanceBinSize value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{DistanceBinSize} = $Value;

  return $This;
}

# Set normalization methodology to use for scaling the number of bits-set or unique keys
# in fingerprints...
#
sub SetNormalizationMethodology {
  my($This, $Value) = @_;

  if ($Value !~ /^(ByHeavyAtomsCount|ByPossibleKeysCount|None)$/i) {
    croak "Error: ${ClassName}->SetNormalizationMethodology: NormalizationMethodology value, $Value, is not valid. Supported values: None, ByHeavyAtomsCount or ByPossibleKeysCount...";
  }

  if ($This->{MolecularComplexityType}) {
    if ($This->{MolecularComplexityType} !~ /^(ExtendedConnectivityFingerprints|TopologicalPharmacophoreAtomPairsFingerprints|TopologicalPharmacophoreAtomTripletsFingerprints)$/i) {
      croak "Error: ${ClassName}->SetNormalizationMethodology: Normalization is not supported for MolecularComplexityType: $This->{MolecularComplexityType}. Valid MolecularComplexityType values: ExtendedConnectivityFingerprints, TopologicalPharmacophoreAtomPairsFingerprints, or TopologicalPharmacophoreAtomTripletsFingerprints...\n";
    }

    if ($This->{MolecularComplexityType} =~ /^ExtendedConnectivityFingerprints$/i && $Value !~ /^(ByHeavyAtomsCount|None)$/i) {
      croak "Error: ${ClassName}->SetNormalizationMethodology: NormalizationMethodology value, $Value, is not valid for MolecularComplexityType: $This->{MolecularComplexityType}. Supported values: None or ByHeavyAtomsCount...";
    }

    if ($This->{MolecularComplexityType} =~ /^(TopologicalPharmacophoreAtomPairsFingerprints|TopologicalPharmacophoreAtomTripletsFingerprints)$/i && $Value !~ /^(ByPossibleKeysCount|None)$/i) {
      croak "Error: ${ClassName}->SetNormalizationMethodology: NormalizationMethodology value, $Value, is not valid for MolecularComplexityType: $This->{MolecularComplexityType}. Supported values: None or ByPossibleKeysCount...";
    }
  }

  $This->{NormalizationMethodology} = $Value;

  return $This;
}

# Set intial atom identifier type..
#
sub SetAtomIdentifierType {
  my($This, $IdentifierType) = @_;

  if ($IdentifierType !~ /^(AtomicInvariantsAtomTypes|FunctionalClassAtomTypes|DREIDINGAtomTypes|EStateAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
    croak "Error: ${ClassName}->SetAtomIdentifierType: Specified value, $IdentifierType, for AtomIdentifierType is not vaild. Supported types in current release of MayaChemTools: AtomicInvariantsAtomTypes, FunctionalClassAtomTypes, DREIDINGAtomTypes, EStateAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes and UFFAtomTypes.";
  }

  # FunctionalClassAtomTypes is the only valid atom identifier type for pharmacophore fingerprints...
  if ($This->{MolecularComplexityType} =~ /^(TopologicalPharmacophoreAtomPairsFingerprints|TopologicalPharmacophoreAtomTripletsFingerprints)$/i) {
    if ($IdentifierType !~ /^FunctionalClassAtomTypes$/i) {
      croak "Error: ${ClassName}->SetAtomIdentifierType: Specified value, $IdentifierType, for AtomIdentifierType is not vaild. Supported type for $This->{MolecularComplexityType} complexity type: FunctionalClassAtomTypes.";
    }
  }

  if ($This->{AtomIdentifierType}) {
    croak "Error: ${ClassName}->SetAtomIdentifierType: Can't change intial atom identifier type:  It's already set...";
  }

  $This->{AtomIdentifierType} = $IdentifierType;

  # Initialize identifier type information...
  $This->_InitializeAtomIdentifierTypeInformation();

  return $This;
}

# Calculate molecular complexity [ Ref 117-119 ] of a molecule using its fingerprints.
#
# The current release of MayaChemTools supports calculation of molecular complexity
# corresponding to the number of bits-set or unique keys in molecular fingerprints. The
# following types of fingerprints based molecular complexity measures are supported:
#
# AtomTypesFingerprints
# ExtendedConnectivityFingerprints
# MACCSKeys
# PathLengthFingerprints
# TopologicalAtomPairsFingerprints
# TopologicalAtomTripletsFingerprints
# TopologicalAtomTorsionsFingerprints
# TopologicalPharmacophoreAtomPairsFingerprints
# TopologicalPharmacophoreAtomTripletsFingerprints
#
# After the molecular complexity value has been calculated, it can also be normalized by
# by scaling the number of bits-set or unique keys for following types of fingerprints:
#
# ExtendedConnectivityFingerprints
# TopologicalPharmacophoreAtomPairsFingerprints
# TopologicalPharmacophoreAtomTripletsFingerprints
#
# Two types of normalization methodologies are supported: by heavy atoms count for
# extended connectivity fingerprints; by possible keys count for topological pharmacophore
# atom pairs and triplets fingerprints.
#
sub GenerateDescriptors {
  my($This) = @_;

  # Initialize descriptor values...
  $This->_InitializeDescriptorValues();

  # Check availability of molecule...
  if (!$This->{Molecule}) {
    carp "Warning: ${ClassName}->GenerateDescriptors: $This->{Type} molecular descriptors generation didn't succeed: Molecule data is not available: Molecule object hasn't been set...";
    return undef;
  }

  # Calculate descriptor values...
  if (!$This->_CalculateDescriptorValues()) {
    carp "Warning: ${ClassName}->GenerateDescriptors: $This->{Type} molecular descriptors generation didn't succeed: Couldn't calculate MolecularComplexity values corresponding to assigned MolecularComplexity atom types...";
    return undef;
  }

  # Set final descriptor values...
  $This->_SetFinalDescriptorValues();

  return $This;
}

# Calculate molecular complexity value...
#
sub _CalculateDescriptorValues {
  my($This) = @_;
  my($FingerprintsObject, $MethodName);

  # Setup fingerprints object and generate fingerprints...
  $MethodName = "_Setup" . $This->{MolecularComplexityType};
  $FingerprintsObject = $This->$MethodName();

  $FingerprintsObject->GenerateFingerprints();

  # Make sure atom types fingerprints generation is successful...
  if (!$FingerprintsObject->IsFingerprintsGenerationSuccessful()) {
    return undef;
  }

  if (!$This->_CalculateMolecularComplexity($FingerprintsObject)) {
    return undef;
  }

  # Normalize molecular complexity...
  if ($This->{NormalizationMethodology} !~ /^None$/i) {
    if (!$This->_NormalizeMolecularComplexity($FingerprintsObject)) {
      return undef;
    }
  }

  return $This;
}

# Setup atom types fingerprints...
#
sub _SetupAtomTypesFingerprints {
  my($This) = @_;
  my($FingerprintsObject);

  $FingerprintsObject = new Fingerprints::AtomTypesFingerprints('Molecule' => $This->{Molecule}, 'Type' => 'AtomTypesCount', 'AtomIdentifierType' => $This->{AtomIdentifierType},  'IgnoreHydrogens' => 1);
  $This->_SetAtomIdentifierTypeValuesToUse($FingerprintsObject);

  return $FingerprintsObject;
}

# Setup extended connectivity fingerprints...
#
sub _SetupExtendedConnectivityFingerprints {
  my($This) = @_;
  my($FingerprintsObject);

  $FingerprintsObject = new Fingerprints::ExtendedConnectivityFingerprints('Molecule' => $This->{Molecule}, 'Type' => 'ExtendedConnectivity', 'NeighborhoodRadius' => $This->{NeighborhoodRadius}, 'AtomIdentifierType' => $This->{AtomIdentifierType});
  $This->_SetAtomIdentifierTypeValuesToUse($FingerprintsObject);

  return $FingerprintsObject;
}

# Setup MACCS keys...
#
sub _SetupMACCSKeys {
  my($This) = @_;
  my($FingerprintsObject);

  $FingerprintsObject = new Fingerprints::MACCSKeys('Molecule' => $This->{Molecule}, 'Type' => 'MACCSKeyBits', 'Size' => $This->{MACCSKeysSize});

  return $FingerprintsObject;
}

# Set up path length fingerprints...
#
sub _SetupPathLengthFingerprints {
  my($This) = @_;
  my($FingerprintsObject);

  $FingerprintsObject = new Fingerprints::PathLengthFingerprints('Molecule' => $This->{Molecule}, 'Type' => 'PathLengthCount', 'AtomIdentifierType' => $This->{AtomIdentifierType}, 'MinLength' => $This->{MinPathLength}, 'MaxLength' => $This->{MaxPathLength}, 'AllowRings' => 1, 'AllowSharedBonds' => 1, 'UseBondSymbols' => $This->{UseBondSymbols}, 'UseUniquePaths' => 1);
  $This->_SetAtomIdentifierTypeValuesToUse($FingerprintsObject);

  return $FingerprintsObject;
}

# Setup topological atom pairs fingerprints...
#
sub _SetupTopologicalAtomPairsFingerprints {
  my($This) = @_;
  my($FingerprintsObject);

  $FingerprintsObject = new Fingerprints::TopologicalAtomPairsFingerprints('Molecule' => $This->{Molecule}, 'MinDistance' => $This->{MinDistance}, 'MaxDistance' => $This->{MaxDistance}, 'AtomIdentifierType' => $This->{AtomIdentifierType});
  $This->_SetAtomIdentifierTypeValuesToUse($FingerprintsObject);

  return $FingerprintsObject;
}

# Setup topological atom triplets fingerprints...
#
sub _SetupTopologicalAtomTripletsFingerprints {
  my($This) = @_;
  my($FingerprintsObject);

  $FingerprintsObject = new Fingerprints::TopologicalAtomTripletsFingerprints('Molecule' => $This->{Molecule}, 'MinDistance' => $This->{MinDistance}, 'MaxDistance' => $This->{MaxDistance}, 'UseTriangleInequality' => $This->{UseTriangleInequality}, 'AtomIdentifierType' => $This->{AtomIdentifierType});
  $This->_SetAtomIdentifierTypeValuesToUse($FingerprintsObject);

  return $FingerprintsObject;
}

# Setup topological atom torsions fingerprints...
#
sub _SetupTopologicalAtomTorsionsFingerprints {
  my($This) = @_;
  my($FingerprintsObject);

  $FingerprintsObject = new Fingerprints::TopologicalAtomTorsionsFingerprints('Molecule' => $This->{Molecule},  'AtomIdentifierType' => $This->{AtomIdentifierType});

  $This->_SetAtomIdentifierTypeValuesToUse($FingerprintsObject);

  return $FingerprintsObject;
}

# Setup TopologicalPharmacophoreAtomPairsFingerprints...
#
sub _SetupTopologicalPharmacophoreAtomPairsFingerprints {
  my($This) = @_;
  my($FingerprintsObject, $AtomPairsSetSizeToUse);

  # Use fixed size to get total number of possible keys for normalization...
  $AtomPairsSetSizeToUse = ($This->{NormalizationMethodology} =~ /^ByPossibleKeysCount$/i) ? 'FixedSize' : 'ArbitrarySize';

  $FingerprintsObject = new Fingerprints::TopologicalPharmacophoreAtomPairsFingerprints('Molecule' => $This->{Molecule}, 'AtomPairsSetSizeToUse' => $AtomPairsSetSizeToUse, 'MinDistance' => $This->{MinDistance}, 'MaxDistance' => $This->{MaxDistance}, 'AtomTypesToUse' => \@{$This->{FunctionalClassesToUse}}, 'NormalizationMethodology' => 'None', 'ValuesPrecision' => 2);

  return $FingerprintsObject;
}

# Setup TopologicalPharmacophoreAtomTripletsFingerprints...
#
sub _SetupTopologicalPharmacophoreAtomTripletsFingerprints {
  my($This) = @_;
  my($FingerprintsObject, $AtomTripletsSetSizeToUse);

  # Use fixed size to get total number of possible keys for normalization...
  $AtomTripletsSetSizeToUse = ($This->{NormalizationMethodology} =~ /^ByPossibleKeysCount$/i) ? 'FixedSize' : 'ArbitrarySize';

  $FingerprintsObject = new Fingerprints::TopologicalPharmacophoreAtomTripletsFingerprints('Molecule' => $This->{Molecule}, 'AtomTripletsSetSizeToUse' => $AtomTripletsSetSizeToUse, 'MinDistance' => $This->{MinDistance}, 'MaxDistance' => $This->{MaxDistance}, 'DistanceBinSize' => $This->{DistanceBinSize}, 'UseTriangleInequality' => $This->{UseTriangleInequality}, 'AtomTypesToUse' => \@{$This->{FunctionalClassesToUse}});

  return $FingerprintsObject;
}

# Normalize molecular complexity value...
#
sub _NormalizeMolecularComplexity {
  my($This, $FingerprintsObject) = @_;

  if ($This->{MolecularComplexityType} =~ /^ExtendedConnectivityFingerprints$/i && $This->{NormalizationMethodology} =~ /^ByHeavyAtomsCount$/i) {
    return $This->_NormalizeMolecularComplexityByHeavyAtomsCount($FingerprintsObject);
  }
  elsif ($This->{MolecularComplexityType} =~ /^(TopologicalPharmacophoreAtomPairsFingerprints|TopologicalPharmacophoreAtomTripletsFingerprints)$/i && $This->{NormalizationMethodology} =~ /^ByPossibleKeysCount$/i) {
    return $This->_NormalizeMolecularComplexityByPossibleKeysCount($FingerprintsObject);
  }
  else {
    warn "Warning: ${ClassName}->_NormalizeMolecularComplexity: NormalizationMethodology value, $This->{NormalizationMethodology}, is not valid. Supported values: ByHeavyAtomsCount or ByPossibleKeysCount...";
  }
  return undef;
}

# Normalize molecular complexity value by heavy atom count...
#
sub _NormalizeMolecularComplexityByHeavyAtomsCount {
  my($This, $FingerprintsObject) = @_;
  my($NumOfHeavyAtoms, $NormalizedComplexity);

  $NumOfHeavyAtoms = $This->{Molecule}->GetNumOfHeavyAtoms();
  if (!$NumOfHeavyAtoms) {
    return $This;
  }

  $NormalizedComplexity = $This->{MolecularComplexity} / $NumOfHeavyAtoms;
  $This->{MolecularComplexity} = MathUtil::round($NormalizedComplexity, 2) + 0;

  return $This;
}

# Normalize molecular complexity value by possible keys count...
#
sub _NormalizeMolecularComplexityByPossibleKeysCount {
  my($This, $FingerprintsObject) = @_;
  my($NumOfPossibleKeys, $NormalizedComplexity);

  $NumOfPossibleKeys = $FingerprintsObject->GetFingerprintsVector()->GetNumOfValues();
  if (!$NumOfPossibleKeys) {
    return $This;
  }

  $NormalizedComplexity = $This->{MolecularComplexity} / $NumOfPossibleKeys;
  $This->{MolecularComplexity} = MathUtil::round($NormalizedComplexity, 2) + 0;

  return $This;
}

# Calculate molecular complexity value using fingerprints objects...
#
sub _CalculateMolecularComplexity {
  my($This, $FingerprintsObject) = @_;

  if ($FingerprintsObject->GetVectorType() =~ /^FingerprintsBitVector$/i) {
    return $This->_CalculateMolecularComplexityUsingFingerprintsBitVector($FingerprintsObject->GetFingerprintsBitVector());
  }
  elsif ($FingerprintsObject->GetVectorType() =~ /^FingerprintsVector$/i) {
    return $This->_CalculateMolecularComplexityUsingFingerprintsVector($FingerprintsObject->GetFingerprintsVector());
  }
  else {
    warn "Warning: ${ClassName}->_CalculateMolecularComplexity: Fingerprints vector type  is not valid. Supported values: FingerprintsBitVector or FingerprintsVector...";
  }

  return undef;
}

# Calculate molecular complexity value using fingerprints vector...
#
sub _CalculateMolecularComplexityUsingFingerprintsVector {
  my($This, $FingerprintsVector) = @_;

  $This->{MolecularComplexity} = ($FingerprintsVector->GetType() =~ /^(OrderedNumericalValues|NumericalValues)$/i) ? $FingerprintsVector->GetNumOfNonZeroValues() : $FingerprintsVector->GetNumOfValues();

  return $This;
}

# Calculate molecular complexity value using fingerprints vector...
#
sub _CalculateMolecularComplexityUsingFingerprintsBitVector {
  my($This, $FingerprintsBitVector) = @_;

  $This->{MolecularComplexity} = $FingerprintsBitVector->GetNumOfSetBits();

  return $This;
}

# Setup final descriptor values...
#
sub _SetFinalDescriptorValues {
  my($This) = @_;

  $This->{DescriptorsGenerated} = 1;

  $This->SetDescriptorValues($This->{MolecularComplexity});

  return $This;
}

# Set atom identifier type to use for generating fingerprints...
#
sub _SetAtomIdentifierTypeValuesToUse {
  my($This, $FingerprintsObject) = @_;

  if ($This->{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
    $FingerprintsObject->SetAtomicInvariantsToUse(\@{$This->{AtomicInvariantsToUse}});
  }
  elsif ($This->{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
    $FingerprintsObject->SetFunctionalClassesToUse(\@{$This->{FunctionalClassesToUse}});
  }
  elsif ($This->{AtomIdentifierType} =~ /^(DREIDINGAtomTypes|EStateAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
    # Nothing to do for now...
  }
  else {
    croak "Error: The value specified, $This->{AtomIdentifierType}, for option \"-a, --AtomIdentifierType\" is not valid. Supported atom identifier types in current release of MayaChemTools: AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes\n";
  }
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

  @{$This->{AtomicInvariantsToUse}} = ();

  if ($This->{MolecularComplexityType} =~ /^(AtomTypesFingerprints|TopologicalAtomPairsFingerprints|TopologicalAtomTripletsFingerprints|TopologicalAtomTorsionsFingerprints)$/i) {
    @{$This->{AtomicInvariantsToUse}} = ('AS', 'X', 'BO', 'H', 'FC');
  }
  elsif ($This->{MolecularComplexityType} =~ /^ExtendedConnectivityFingerprints$/i) {
    @{$This->{AtomicInvariantsToUse}} = ('AS', 'X', 'BO', 'H', 'FC', 'MN');
  }
  elsif ($This->{MolecularComplexityType} =~ /^PathLengthFingerprints$/i) {
    @{$This->{AtomicInvariantsToUse}} = ('AS');
  }

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

  @{$This->{FunctionalClassesToUse}} = ();

  if ($This->{MolecularComplexityType} =~ /^(AtomTypesFingerprints|ExtendedConnectivityFingerprints|PathLengthFingerprints|TopologicalAtomPairsFingerprints|TopologicalAtomTripletsFingerprints|TopologicalAtomTorsionsFingerprints)$/i) {
    @{$This->{FunctionalClassesToUse}} = ('HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal');
  }
  elsif ($This->{MolecularComplexityType} =~ /^TopologicalPharmacophoreAtomPairsFingerprints$/i) {
    @{$This->{FunctionalClassesToUse}} = ('HBD', 'HBA', 'PI', 'NI', 'H');
  }
  elsif ($This->{MolecularComplexityType} =~ /^TopologicalPharmacophoreAtomTripletsFingerprints$/i) {
    @{$This->{FunctionalClassesToUse}} = ('HBD', 'HBA', 'PI', 'NI', 'H', 'Ar');
  }

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
    carp "Warning: ${ClassName}->SetAtomicInvariantsToUse: AtomicInvariantsToUse can't be set for InitialAtomIdentifierType of $This->{AtomIdentifierType}...";
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

# Return a string containg data for MolecularComplexityDescriptors object...
#
sub StringifyMolecularComplexityDescriptors {
  my($This) = @_;
  my($ComplexityDescriptorsString, $Nothing);

  $ComplexityDescriptorsString = "MolecularDescriptorType: $This->{Type}; MolecularComplexityType: $This->{MolecularComplexityType}; " . $This->_StringifyDescriptorNamesAndValues();

  # Setup fingerprints specific information...
  if ($This->{MolecularComplexityType} =~ /^MACCSKeys$/i) {
    $ComplexityDescriptorsString .= "; MACCSKeysSize = $This->{MACCSKeysSize}";
  }
  elsif ($This->{MolecularComplexityType} =~ /^ExtendedConnectivityFingerprints$/i) {
    $ComplexityDescriptorsString .= "; NeighborhoodRadius = $This->{NeighborhoodRadius}; NormalizationMethodology = $This->{NormalizationMethodology}";
  }
  elsif ($This->{MolecularComplexityType} =~ /^PathLengthFingerprints$/i) {
    $ComplexityDescriptorsString .= "; MinPathLength = $This->{MinPathLength}; MaxPathLength = $This->{MaxPathLength}; UseBondSymbols: " . ($This->{UseBondSymbols} ? "Yes" : "No");
  }
  elsif ($This->{MolecularComplexityType} =~ /^TopologicalAtomPairsFingerprints$/i) {
    $ComplexityDescriptorsString .= "; MinDistance = $This->{MinDistance}; MaxDistance = $This->{MaxDistance}";
  }
  elsif ($This->{MolecularComplexityType} =~ /^TopologicalAtomTripletsFingerprints$/i) {
    $ComplexityDescriptorsString .= "; MinDistance = $This->{MinDistance}; MaxDistance = $This->{MaxDistance}; UseTriangleInequality: " . ($This->{UseTriangleInequality} ? "Yes" : "No");
  }
  elsif ($This->{MolecularComplexityType} =~ /^TopologicalAtomTorsionsFingerprints$/i) {
    $ComplexityDescriptorsString .= "; MinDistance = $This->{MinDistance}; MaxDistance = $This->{MaxDistance}";
  }
  elsif ($This->{MolecularComplexityType} =~ /^TopologicalPharmacophoreAtomPairsFingerprints$/i) {
    $ComplexityDescriptorsString .= "; MinDistance = $This->{MinDistance}; MaxDistance = $This->{MaxDistance}; NormalizationMethodology = $This->{NormalizationMethodology}";
  }
  elsif ($This->{MolecularComplexityType} =~ /^TopologicalPharmacophoreAtomTripletsFingerprints$/i) {
    $ComplexityDescriptorsString .= "; MinDistance = $This->{MinDistance}; MaxDistance = $This->{MaxDistance}; NormalizationMethodology = $This->{NormalizationMethodology};  DistanceBinSize: $This->{DistanceBinSize}; UseTriangleInequality: " . ($This->{UseTriangleInequality} ? "Yes" : "No");
  }

  # Setup atom identifier information...
  if ($This->{MolecularComplexityType} =~ /^(AtomTypesFingerprints|ExtendedConnectivityFingerprints|PathLengthFingerprints|TopologicalAtomPairsFingerprints|TopologicalAtomTripletsFingerprints|TopologicalAtomTorsionsFingerprints|TopologicalPharmacophoreAtomPairsFingerprints|TopologicalPharmacophoreAtomTripletsFingerprints)$/i) {
    $ComplexityDescriptorsString .= "; AtomIdentifierType = $This->{AtomIdentifierType}";

    if ($This->{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
      my($AtomicInvariant, @AtomicInvariants, @AtomicInvariantsOrder, %AvailableAtomicInvariants);

      @AtomicInvariantsOrder = AtomTypes::AtomicInvariantsAtomTypes::GetAtomicInvariantsOrder();
      %AvailableAtomicInvariants = AtomTypes::AtomicInvariantsAtomTypes::GetAvailableAtomicInvariants();

      for $AtomicInvariant (@AtomicInvariantsOrder) {
	push @AtomicInvariants, "$AtomicInvariant: $AvailableAtomicInvariants{$AtomicInvariant}";
      }

      $ComplexityDescriptorsString .= "; AtomicInvariantsToUse: <" . TextUtil::JoinWords(\@{$This->{AtomicInvariantsToUse}}, ", ", 0) . ">";
      $ComplexityDescriptorsString .= "; AtomicInvariantsOrder: <" . TextUtil::JoinWords(\@AtomicInvariantsOrder, ", ", 0) . ">";
      $ComplexityDescriptorsString .= "; AvailableAtomicInvariants: <" . TextUtil::JoinWords(\@AtomicInvariants, ", ", 0) . ">";
    }
    elsif ($This->{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
      my($FunctionalClass, @FunctionalClasses, @FunctionalClassesOrder, %AvailableFunctionalClasses);

      @FunctionalClassesOrder = AtomTypes::FunctionalClassAtomTypes::GetFunctionalClassesOrder();
      %AvailableFunctionalClasses = AtomTypes::FunctionalClassAtomTypes::GetAvailableFunctionalClasses();

      for $FunctionalClass (@FunctionalClassesOrder) {
	push @FunctionalClasses, "$FunctionalClass: $AvailableFunctionalClasses{$FunctionalClass}";
      }

      $ComplexityDescriptorsString .= "; FunctionalClassesToUse: <" . TextUtil::JoinWords(\@{$This->{FunctionalClassesToUse}}, ", ", 0) . ">";
      $ComplexityDescriptorsString .= "; FunctionalClassesOrder: <" . TextUtil::JoinWords(\@FunctionalClassesOrder, ", ", 0) . ">";
      $ComplexityDescriptorsString .= "; AvailableFunctionalClasses: <" . TextUtil::JoinWords(\@FunctionalClasses, ", ", 0) . ">";
    }
  }
  return $ComplexityDescriptorsString;
}

# Is it a MolecularComplexityDescriptors object?
sub _IsMolecularComplexityDescriptors {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

MolecularComplexityDescriptors

=head1 SYNOPSIS

use MolecularDescriptors::MolecularComplexityDescriptors;

use MolecularDescriptors::MolecularComplexityDescriptors qw(:all);

=head1 DESCRIPTION

B<MolecularComplexityDescriptors> class provides the following methods:

new, GenerateDescriptors, GetDescriptorNames,
GetMolecularComplexityTypeAbbreviation, MACCSKeysSize, SetAtomIdentifierType,
SetAtomicInvariantsToUse, SetDistanceBinSize, SetFunctionalClassesToUse,
SetMaxDistance, SetMaxPathLength, SetMinDistance, SetMinPathLength,
SetMolecularComplexityType, SetNeighborhoodRadius, SetNormalizationMethodology,
StringifyMolecularComplexityDescriptors

B<MolecularComplexityDescriptors> is derived from B<MolecularDescriptors> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<MolecularComplexityDescriptors>, B<MolecularDescriptors> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The current release of MayaChemTools supports calculation of molecular complexity using
I<MolecularComplexityType> parameter corresponding to number of bits-set or unique
keys [ Ref 117-119 ] in molecular  fingerprints. The valid values for I<MolecularComplexityType>
are:

    AtomTypesFingerprints
    ExtendedConnectivityFingerprints
    MACCSKeys
    PathLengthFingerprints
    TopologicalAtomPairsFingerprints
    TopologicalAtomTripletsFingerprints
    TopologicalAtomTorsionsFingerprints
    TopologicalPharmacophoreAtomPairsFingerprints
    TopologicalPharmacophoreAtomTripletsFingerprints

Default value for I<MolecularComplexityType>: I<MACCSKeys>.

I<AtomIdentifierType> parameter name corresponds to atom types used during generation of
fingerprints. The valid values for I<AtomIdentifierType> are: I<AtomicInvariantsAtomTypes,
DREIDINGAtomTypes, EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes,
SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes>. I<AtomicInvariantsAtomTypes>
is not supported for following values of I<MolecularComplexityType>: I<MACCSKeys,
TopologicalPharmacophoreAtomPairsFingerprints, TopologicalPharmacophoreAtomTripletsFingerprints>.
I<FunctionalClassAtomTypes> is the only valid value of I<AtomIdentifierType> for topological
pharmacophore fingerprints.

Default value for I<AtomIdentifierType>: I<AtomicInvariantsAtomTypes> for all fingerprints;
I<FunctionalClassAtomTypes> for topological pharmacophore fingerprints.

I<AtomicInvariantsToUse> parameter name and values are used during I<AtomicInvariantsAtomTypes>
value of parameter I<AtomIdentifierType>. It's a list of space separated valid atomic invariant atom types.

Possible values for atomic invariants are: I<AS, X, BO,  LBO, SB, DB, TB, H, Ar, RA, FC, MN, SM>.
Default value for I<AtomicInvariantsToUse> parameter are set differently for different fingerprints
using I<MolecularComplexityType> parameter as shown below:

    MolecularComplexityType              AtomicInvariantsToUse

    AtomTypesFingerprints                AS X BO H FC
    TopologicalAtomPairsFingerprints     AS X BO H FC
    TopologicalAtomTripletsFingerprints  AS X BO H FC
    TopologicalAtomTorsionsFingerprints  AS X BO H FC

    ExtendedConnectivityFingerprints     AS X  BO H FC MN
    PathLengthFingerprints               AS

I<FunctionalClassesToUse> parameter name and values are used during I<FunctionalClassAtomTypes>
value of parameter I<AtomIdentifierType>. It's a list of space separated valid atomic invariant atom types.

Possible values for atom functional classes are: I<Ar, CA, H, HBA, HBD, Hal, NI, PI, RA>.

Default value for I<FunctionalClassesToUse> parameter is set to:

    HBD HBA PI NI Ar Hal

for all fingerprints except for the following two I<MolecularComplexityType> fingerints:

    MolecularComplexityType                           FunctionalClassesToUse

    TopologicalPharmacophoreAtomPairsFingerprints     HBD HBA P, NI H
    TopologicalPharmacophoreAtomTripletsFingerprints  HBD HBA PI NI H Ar

I<MACCSKeysSize> parameter name is only used during I<MACCSKeys> value of
I<MolecularComplexityType> and corresponds to size of MACCS key set. Possible
values: I<166 or 322>. Default value: I<166>.

I<NeighborhoodRadius> parameter name is only used during I<ExtendedConnectivityFingerprints>
value of I<MolecularComplexityType> and corresponds to atomic neighborhoods radius for
generating extended connectivity fingerprints. Possible values: positive integer. Default value:
I<2>.

I<MinPathLength> and I<MaxPathLength> parameters are only used during I<PathLengthFingerprints>
value of I<MolecularComplexityType> and correspond to minimum and maximum path lengths to use
for generating path length fingerprints. Possible values: positive integers. Default value: I<MinPathLength - 1>;
I<MaxPathLength - 8>.

I<UseBondSymbols> parameter is only used during I<PathLengthFingerprints> value of
I<MolecularComplexityType> and indicates whether bond symbols are included in atom path
strings used to generate path length fingerprints. Possible value: I<Yes or No>. Default value:
I<Yes>.

I<MinDistance> and I<MaxDistance> parameters are only used during I<TopologicalAtomPairsFingerprints>
and I<TopologicalAtomTripletsFingerprints> values of I<MolecularComplexityType> and correspond to
minimum and maximum bond distance between atom pairs during topological pharmacophore fingerprints.
Possible values: positive integers. Default value: I<MinDistance - 1>; I<MaxDistance - 10>.

I<UseTriangleInequality> parameter is used during these values for I<MolecularComplexityType>:
I<TopologicalAtomTripletsFingerprints> and I<TopologicalPharmacophoreAtomTripletsFingerprints>.
Possible values: I<Yes or No>. It determines wheter to apply triangle inequality to distance triplets.
Default value: I<TopologicalAtomTripletsFingerprints - No>;
I<TopologicalPharmacophoreAtomTripletsFingerprints - Yes>.

I<DistanceBinSize> parameter is used during I<TopologicalPharmacophoreAtomTripletsFingerprints>
value of I<MolecularComplexityType> and corresponds to distance bin size used for binning
distances during generation of topological pharmacophore atom triplets fingerprints. Possible
value: positive integer. Default value: I<2>.

I<NormalizationMethodology> is only used for these values for I<MolecularComplexityType>:
I<ExtendedConnectivityFingerprints>, I<TopologicalPharmacophoreAtomPairsFingerprints>
and I<TopologicalPharmacophoreAtomTripletsFingerprints>. It corresponds to normalization
methodology to use for scaling the number of bits-set or unique keys during generation of
fingerprints. Possible values during I<ExtendedConnectivityFingerprints>: I<None or
ByHeavyAtomsCount>; Default value: I<None>. Possible values during topological
pharmacophore atom pairs and triplets fingerprints: I<None or ByPossibleKeysCount>;
Default value: I<None>. I<ByPossibleKeysCount> corresponds to total number of
possible topological pharmacophore atom pairs or triplets in a molecule.

=head2 METHODS

=over 4

=item B<new>

    $NewMolecularComplexityDescriptors = new MolecularDescriptors::
                                         MolecularComplexityDescriptors(
                                             %NamesAndValues);

Using specified I<MolecularComplexityDescriptors> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<MolecularComplexityDescriptors>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'MolecularComplexity'
    MolecularComplexityType = 'MACCSKeys'
    AtomIdentifierType = ''
    MACCSKeysSize = 166
    NeighborhoodRadius = 2
    MinPathLength = 1
    MaxPathLength = 8
    UseBondSymbols = 1
    MinDistance = 1
    MaxDistance = 10
    UseTriangleInequality = ''
    DistanceBinSize = 2
    NormalizationMethodology = 'None'
    @DescriptorNames = ('MolecularComplexity')
    @DescriptorValues = ('None')

Examples:

    $MolecularComplexityDescriptors = new MolecularDescriptors::
                                      MolecularComplexityDescriptors(
                                      'Molecule' => $Molecule);

    $MolecularComplexityDescriptors = new MolecularDescriptors::
                                      MolecularComplexityDescriptors();

    $MolecularComplexityDescriptors->SetMolecule($Molecule);
    $MolecularComplexityDescriptors->GenerateDescriptors();
    print "MolecularComplexityDescriptors: $MolecularComplexityDescriptors\n";


=item B<GenerateDescriptors>

    $MolecularComplexityDescriptors->GenerateDescriptors();

Calculates MolecularComplexity value for a molecule and returns I<MolecularComplexityDescriptors>.

=item B<GetDescriptorNames>

    @DescriptorNames = $MolecularComplexityDescriptors->GetDescriptorNames();
    @DescriptorNames = MolecularDescriptors::MolecularComplexityDescriptors::
                       GetDescriptorNames();

Returns all available descriptor names as an array.

=item B<GetMolecularComplexityTypeAbbreviation>

    $Abbrev = $MolecularComplexityDescriptors->
                  GetMolecularComplexityTypeAbbreviation();
    $Abbrev = MolecularDescriptors::MolecularComplexityDescriptors::
                  GetMolecularComplexityTypeAbbreviation($ComplexityType);

Returns abbreviation for a specified molecular complexity type or corresponding to
I<MolecularComplexityDescriptors> object.

=item B<SetMACCSKeysSize>

    $MolecularComplexityDescriptors->MACCSKeysSize($Size);

Sets MACCS keys size and returns I<MolecularComplexityDescriptors>.

=item B<SetAtomIdentifierType>

    $MolecularComplexityDescriptors->SetAtomIdentifierType($IdentifierType);

Sets atom I<IdentifierType> to use during fingerprints generation corresponding to
I<MolecularComplexityType> and returns I<MolecularComplexityDescriptors>.

Possible values: I<AtomicInvariantsAtomTypes, DREIDINGAtomTypes, EStateAtomTypes,
FunctionalClassAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes,
TPSAAtomTypes, UFFAtomTypes>.

=item B<SetAtomicInvariantsToUse>

    $MolecularComplexityDescriptors->SetAtomicInvariantsToUse($ValuesRef);
    $MolecularComplexityDescriptors->SetAtomicInvariantsToUse(@Values);

Sets atomic invariants to use during I<AtomicInvariantsAtomTypes> value of I<AtomIdentifierType>
for fingerprints generation and returns I<MolecularComplexityDescriptors>.

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

=item B<SetDistanceBinSize>

    $MolecularComplexityDescriptors->SetDistanceBinSize($BinSize);

Sets distance bin size used to bin distances between atom pairs in atom triplets for
topological pharmacophore atom triplets fingerprints generation and returns
I<MolecularComplexityDescriptors>.

=item B<SetFunctionalClassesToUse>

    $MolecularComplexityDescriptors->SetFunctionalClassesToUse($ValuesRef);
    $MolecularComplexityDescriptors->SetFunctionalClassesToUse(@Values);

Sets functional classes invariants to use during I<FunctionalClassAtomTypes> value of I<AtomIdentifierType>
for fingerprints generation and returns I<MolecularComplexityDescriptors>.

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

    $MolecularComplexityDescriptors->SetMaxDistance($MaxDistance);

Sets maximum distance to use during topological atom pairs and triplets fingerprints
generation and returns I<MolecularComplexityDescriptors>.

=item B<SetMaxPathLength>

    $MolecularComplexityDescriptors->SetMaxPathLength($Length);

Sets maximum path length to use during path length fingerprints generation and returns
I<MolecularComplexityDescriptors>.

=item B<SetMinDistance>

    $MolecularComplexityDescriptors->SetMinDistance($MinDistance);

Sets minimum distance to use during topological atom pairs and triplets fingerprints
generation and returns I<MolecularComplexityDescriptors>.

=item B<SetMinPathLength>

    $MolecularComplexityDescriptors->SetMinPathLength($MinPathLength);

Sets minimum path length to use during path length fingerprints generation and returns
I<MolecularComplexityDescriptors>.

=item B<SetMolecularComplexityType>

    $MolecularComplexityDescriptors->SetMolecularComplexityType($ComplexityType);

Sets molecular complexity type to use for calculating its value and returns
I<MolecularComplexityDescriptors>.

=item B<SetNeighborhoodRadius>

    $MolecularComplexityDescriptors->SetNeighborhoodRadius($Radius);

Sets neighborhood radius to use during extended connectivity fingerprints generation and
returns I<MolecularComplexityDescriptors>.

=item B<SetNormalizationMethodology>

    $MolecularComplexityDescriptors->SetNormalizationMethodology($Methodology);

Sets normalization methodology to use during calculation of molecular complexity
corresponding to extended connectivity, topological pharmacophore atom pairs and
tripletes fingerprints returns I<MolecularComplexityDescriptors>.

=item B<StringifyMolecularComplexityDescriptors>

    $String = $MolecularComplexityDescriptors->
                  StringifyMolecularComplexityDescriptors();

Returns a string containing information about I<MolecularComplexityDescriptors> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

MolecularDescriptors.pm, MolecularDescriptorsGenerator.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
