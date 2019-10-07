package Fingerprints::TopologicalPharmacophoreAtomPairsFingerprints;
#
# File: TopologicalPharmacophoreAtomPairsFingerprints.pm
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
use MathUtil ();
use Molecule;
use AtomTypes::FunctionalClassAtomTypes;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Fingerprints::Fingerprints Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyTopologicalPharmacophoreAtomPairsFingerprints';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeTopologicalPharmacophoreAtomPairsFingerprints();

  $This->_InitializeTopologicalPharmacophoreAtomPairsFingerprintsProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeTopologicalPharmacophoreAtomPairsFingerprints {
  my($This) = @_;

  # Type of fingerprint...
  $This->{Type} = 'TopologicalPharmacophoreAtomPairs';

  # Type of vector...
  $This->{VectorType} = 'FingerprintsVector';

  # AtomPairsSetSizeToUse...
  #
  # ArbitrarySize - Corrresponds to atom pairs with non-zero count
  # FixedSize - Corresponds to all atom pairs with zero and non-zero count
  #
  # Possible values: ArbitrarySize or FixedSize. Default: ArbitrarySize
  #
  $This->{AtomPairsSetSizeToUse} = '';

  # Type of FingerprintsVector...
  #
  # OrderedNumericalValues - For ArbitrarySize value of AtomPairsSetSizeToUse
  # NumericalValues - For FixedSize value of AtomPairsSetSizeToUse
  #
  # Possible values: OrderedNumericalValues or NumericalValues. Default: NumericalValues
  #
  $This->{FingerprintsVectorType} = '';

  # Vector values precision for real values which might be generated after
  # normalization and fuzzification...
  $This->{ValuesPrecision} = 2;

  # Minimum and maximum bond distance between pharmacophore atom paris...
  $This->{MinDistance} = 1;
  $This->{MaxDistance} = 10;

  # Initialize atom types and weight information...
  $This->_InitializePharmacophoreAtomTypesAndWeightInformation();

  # Normalization methodology to use for scaling the occurance count of pharmacophore atom
  # pairs at various distances.
  #
  # Possible values: None, ByHeavyAtomsCount, ByAtomTypesCount. Default: None
  #
  $This->{NormalizationMethodology} = 'None';

  # Initialize fuzzification parameters...
  #
  $This->_InitializeFuzzificationInformation();

  # Pharmacophore types assigned to each heavy atom...
  #
  %{$This->{AssignedAtomTypes}} = ();

  # Assigned Atom types count of each type in the molecule...
  #
  %{$This->{AssignedAtomTypesCount}} = ();

  # All pharmacophore atom pairs between minimum and maximum distance...
  #
  @{$This->{AtomPairsIDs}} = ();
  %{$This->{AtomPairsCount}} = ();
}

# Inialize pharmacophore atom types and weight information...
#
sub _InitializePharmacophoreAtomTypesAndWeightInformation {
  my($This) = @_;

  # Default pharmacophore atom types to use for atom pairs fingerprint generation
  # are: HBD, HBA, PI, NI, H
  #
  @{$This->{AtomTypesToUse}} = ();
  @{$This->{AtomTypesToUse}} = sort ('HBD', 'HBA', 'PI', 'NI', 'H');

  # Weight of the various pharmacophore atom types to use for their contribution to atom
  # pair interaction. It allows to increase the importance of specific pharmacophore atom
  # types in the generted fingerprints.
  #
  # A value of 0 eliminates the contribution by a particular pharmacophore atom
  # type and 2 doubles its contribution.
  #
  my($AtomType, %AvailableAtomTypes);

  %AvailableAtomTypes = AtomTypes::FunctionalClassAtomTypes::GetAvailableFunctionalClasses();

  %{$This->{AtomTypesWeight}} = ();
  for $AtomType (keys %AvailableAtomTypes) {
    $This->{AtomTypesWeight}{$AtomType} = 1;
  }
  return $This;
}

# Initialize fuzzification information...
#
sub _InitializeFuzzificationInformation {
  my($This) = @_;

  # To fuzz or not to fuzz atom pairs count. Default: No fuzzication
  #
  $This->{FuzzifyAtomPairsCount} = 0;

  # When to fuzz atom pair count...
  #
  # Possible values: BeforeNormalization or AfterNormalization. Default: AfterNormalization
  #
  $This->{FuzzificationMode} = 'AfterNormalization';

  # How to fuzz atom pair count...
  #
  # Possible values: FuzzyBinning or FuzzyBinSmoothing. Default: FuzzyBinning
  #
  $This->{FuzzificationMethodology} = 'FuzzyBinning';

  # By how much to fuzz atom pairs count...
  #
  $This->{FuzzFactor} = 0.15;

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize object properties....
sub _InitializeTopologicalPharmacophoreAtomPairsFingerprintsProperties {
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

  $This->_InitializeTopologicalPharmacophoreAtomPairsFingerprintsVector();

  return $This;
}

# Initialize fingerprints vector...
#
sub _InitializeTopologicalPharmacophoreAtomPairsFingerprintsVector {
  my($This) = @_;

  if (!$This->{AtomPairsSetSizeToUse}) {
    $This->{AtomPairsSetSizeToUse} =  'ArbitrarySize';
  }

  # Vector type and type of values...
  $This->{VectorType} = 'FingerprintsVector';

  if ($This->{AtomPairsSetSizeToUse} =~ /^FixedSize$/i) {
    $This->{FingerprintsVectorType} = 'OrderedNumericalValues';
  }
  else {
    $This->{FingerprintsVectorType} = 'NumericalValues';
  }

  $This->_InitializeFingerprintsVector();
}

# Set atom parits set size to use...
#
sub SetAtomPairsSetSizeToUse {
  my($This, $Value) = @_;

  if ($This->{AtomPairsSetSizeToUse}) {
    croak "Error: ${ClassName}->SetAtomPairsSetSizeToUse: Can't change size:  It's already set...";
  }

  if ($Value !~ /^(ArbitrarySize|FixedSize)$/i) {
    croak "Error: ${ClassName}->SetAtomPairsSetSizeToUse: Unknown AtomPairsSetSizeToUse value: $Value; Supported values: ArbitrarySize or FixedSize";
  }

  $This->{AtomPairsSetSizeToUse} = $Value;

  return $This;
}

# Disable change of AvailableAtomTypes...
#
sub SetAvailableAtomTypes {
  my($This) = @_;

  carp "Warning: ${ClassName}->SetAvailableAtomTypes: AvailableAtomTypes value can't be set...";

  return $This;
}

# Set atom types to use for atom pairs...
#
sub SetAtomTypesToUse {
  my($This, @Values) = @_;
  my($FirstValue, $TypeOfFirstValue, $AtomType, $SpecifiedAtomType, @SpecifiedAtomTypes, @AtomTypesToUse);

  if (!@Values) {
    carp "Warning: ${ClassName}->SetAtomTypesToUse: No values specified...";
    return;
  }

  $FirstValue = $Values[0];
  $TypeOfFirstValue = ref $FirstValue;

  @SpecifiedAtomTypes = ();
  @AtomTypesToUse = ();

  if ($TypeOfFirstValue =~ /^ARRAY/) {
    push @SpecifiedAtomTypes, @{$FirstValue};
  }
  else {
    push @SpecifiedAtomTypes, @Values;
  }

  # Make sure specified AtomTypes are valid...
  for $SpecifiedAtomType (@SpecifiedAtomTypes) {
    if (!AtomTypes::FunctionalClassAtomTypes::IsFunctionalClassAvailable($SpecifiedAtomType)) {
      croak "Error: ${ClassName}->SetAtomTypesToUse: Specified atom type, $SpecifiedAtomType, is not supported...\n ";
    }
    $AtomType = $SpecifiedAtomType;
    push @AtomTypesToUse, $AtomType;
  }

  # Set atom types to use...
  @{$This->{AtomTypesToUse}} = ();
  push @{$This->{AtomTypesToUse}}, sort @AtomTypesToUse;

  return $This;
}

# Set vector values precision for real values which might be generated after
# normalization and fuzzification...
#
sub SetValuesPrecision {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetValuesPrecision: ValuesPrecision value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{ValuesPrecision} = $Value;

  return $This;
}

# Set minimum distance for pharmacophore atom pairs...
#
sub SetMinDistance {
  my($This, $Value) = @_;

  if (!TextUtil::IsInteger($Value)) {
    croak "Error: ${ClassName}->SetMinDistance: MinDistance value, $Value, is not valid:  It must be an integer...";
  }
  $This->{MinDistance} = $Value;

  return $This;
}

# Set maximum distance for pharmacophore atom pairs...
#
sub SetMaxDistance {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetMaxDistance: MaxDistance value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{MaxDistance} = $Value;

  return $This;
}

# Set normalization methodology to use for scaling the occurance count of pharmacophore atom
# pairs over distance range beween minimum and maximum distance.
#
sub SetNormalizationMethodology {
  my($This, $Value) = @_;

  if ($Value !~ /^(ByHeavyAtomsCount|ByAtomTypesCount|None)$/i) {
    croak "Error: ${ClassName}->SetNormalizationMethodology: NormalizationMethodology value, $Value, is not valid. Supported values: None, ByHeavyAtomsCount or ByAtomTypesCount...";
  }

  $This->{NormalizationMethodology} = $Value;

  return $This;
}

# Set weight of the various pharmacophore atom types to use for their contribution to atom
# pair interaction using atom types label and value hash.
#
# It allows to increase the importance of specific pharmacophore atom
# types in the generted fingerprints.
#
# A value of 0 eliminates the contribution by a particular pharmacophore atom
# type and 2 doubles its contribution.
#
sub SetAtomTypesWeight {
  my($This, %AtomTypesWeight) = @_;
  my($AtomType, $Weight);

  while (($AtomType, $Weight) = each %AtomTypesWeight) {
    if (!exists $This->{AtomTypesWeight}{$AtomType}) {
      croak "Error: ${ClassName}->SetAtomTypesWeight: AtomTypeWeight for $AtomType couldn't be set: Unknown atom type...";
    }
    if (!(TextUtil::IsFloat($Weight) && ($Weight >= 0))) {
      croak "Error: ${ClassName}->SetAtomTypesWeight: Specified weight value, $Weight, for AtomType, $AtomType, muts be >= 0...";
    }
    $This->{AtomTypesWeight}{$AtomType}  = $Weight;
  }
}

# Set fuzzification methodology to use for fuzzifying atom pairs count...
#
sub SetFuzzificationMethodology {
  my($This, $Value) = @_;

  if ($Value !~ /^(FuzzyBinning|FuzzyBinSmoothing)$/i) {
    croak "Error: ${ClassName}->SetFuzzificationMethodology: FuzzificationMethodology value, $Value, is not valid. Supported values: FuzzyBinning or FuzzyBinSmoothing...";
  }

  $This->{FuzzificationMethodology} = $Value;

  return $This;
}

# Set fuzzification mode for fuzzifying atom pairs count...
#
sub SetFuzzificationMode {
  my($This, $Value) = @_;

  if ($Value !~ /^(BeforeNormalization|AfterNormalization)$/i) {
    croak "Error: ${ClassName}->SetFuzzificationMode: FuzzificationMode value, $Value, is not valid. Supported values: BeforeNormalization or AfterNormalization...";
  }

  $This->{FuzzificationMode} = $Value;

  return $This;
}

# Set fuzz factor values used for fuzzifying atom pairs count...
#
sub SetFuzzFactor {
  my($This, $Value) = @_;

  if ($This->{FuzzificationMethodology} =~ /^FuzzyBinning$/i) {
    if (!(TextUtil::IsFloat($Value) && $Value >=0 && $Value <= 1.0)) {
      croak "Error: ${ClassName}->SetFuzzFactor: Specified fuzz factor value, $Value, must be >= 0 and <= 1...";
    }
  }
  elsif ($This->{FuzzificationMethodology} =~ /^FuzzyBinSmoothing$/i) {
    if (!(TextUtil::IsFloat($Value) && $Value >=0 && $Value <= 0.5)) {
      croak "Error: ${ClassName}->SetFuzzFactor: Specified fuzz factor value, $Value, must be >= 0 and <= 0.5...";
    }
  }
  else {
    croak "Error: ${ClassName}->SetFuzzFactor: Fuzz factor value can't be changed: Uknown FuzzificationMethodology: $This->{FuzzificationMethodology}...";
  }

  $This->{FuzzFactor} = $Value;

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

  return "$This->{Type}:$This->{AtomPairsSetSizeToUse}:MinDistance$This->{MinDistance}:MaxDistance$This->{MaxDistance}";
}

# Generate topological pharmacophore atom pairs [ Ref 60-62, Ref 65, Ref 68 ] fingerprints...
#
# Methodology:
#   . Generate a distance matrix.
#   . Assign pharmacophore atom types to all the atoms.
#   . Initialize pharmacophore atom pairs basis set for all unique pairs between
#     minimum and maximum distance.
#   . Using distance matrix and pharmacophore atom types, count occurance of
#     unique atom pairs between specified distance range - It corresponds to the
#     correlation-vector for the atom pairs.
#       . Weigh contribution of each atom type to atom pair interaction by its specified
#         weight during occurance count.
#       . Assign count to appropriate distance bin for a specific atom pair
#
#   . Normalize occurance count of pharmacophore atom pairs by heavy atom count
#     or sum of AtomTypeCounts of each pharmacophore atom type in the atom pair
#     at a specific distance.
#
#   . Fuzzify occurance count of pharmacophore atom pairs using FuzzyBinning or
#     FuzzySmothing methodology.
#
# Notes:
#   . Hydrogen atoms are ignored during the fingerprint generation.
#
sub GenerateFingerprints {
  my($This) = @_;

  if ($This->{MinDistance} > $This->{MaxDistance}) {
    croak "Error: ${ClassName}->GenerateTopologicalPharmacophoreAtomPairsFingerprints: No fingerpritns generated: MinDistance, $This->{MinDistance}, must be <= MaxDistance, $This->{MaxDistance}...";
  }

  # Cache appropriate molecule data...
  $This->_SetupMoleculeDataCache();

  # Generate distance matrix...
  if (!$This->_SetupDistanceMatrix()) {
    carp "Warning: ${ClassName}->GenerateFingerprints: Fingerprints generation didn't succeed: Couldn't generate distance matrix...";
    return $This;
  }

  # Assign pharmacohore atom types to all heavy atoms...
  $This->_AssignPharmacophoreAtomTypes();

  # Initialize values of all possible pharmacohore atom pairs...
  $This->_InitializePharmacophoreAtomPairs();

  # Count atom pairs...
  $This->_CountPharmacohoreAtomPairs();

  # Fuzzify atom pairs count...
  if ($This->{FuzzificationMode} =~ /^BeforeNormalization$/i) {
    $This->_FuzzifyPharmacohoreAtomPairsCount();
  }

  # Normalize atom pairs count...
  $This->_NormalizePharmacohoreAtomPairsCount();

  # Fuzzify atom pairs count...
  if ($This->{FuzzificationMode} =~ /^AfterNormalization$/i) {
    $This->_FuzzifyPharmacohoreAtomPairsCount();
  }

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

# Assign pharmacohore atom types to all heavy atoms and count each atom
# types assigned...
#
sub _AssignPharmacophoreAtomTypes {
  my($This) = @_;
  my($Atom, $AtomID, $AtomType, $AssignedAtomType, $FunctionalClassAtomTypes);

  # Assign topological pharmacophore atom types...
  $FunctionalClassAtomTypes = new AtomTypes::FunctionalClassAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => 1, 'FunctionalClassesToUse' => $This->{AtomTypesToUse});
  $FunctionalClassAtomTypes->AssignAtomTypes();

  %{$This->{AssignedAtomTypes}} = ();

  # Initialize assigned atom types count...
  %{$This->{AssignedAtomTypesCount}} = ();
  for $AtomType (@{$This->{AtomTypesToUse}}) {
    $This->{AssignedAtomTypesCount}{$AtomType} = 0;
  }

  $This->{HeavyAtomCount} = 0;

  ATOM: for $Atom (@{$This->{Atoms}}) {
    if ($Atom->IsHydrogen()) {
      next ATOM;
    }
    $This->{HeavyAtomCount} += 1;

    $AtomID = $Atom->GetID();

    # Collect all possible pharmacophore atom types which could be assigned to atom...
    my(@AtomTypes);

    @AtomTypes = ();
    $AssignedAtomType = $FunctionalClassAtomTypes->GetAtomType($Atom);
    if ($AssignedAtomType && $AssignedAtomType !~ /^None$/i) {
      push @AtomTypes, split /\./, $AssignedAtomType;
      for $AtomType (@AtomTypes) {
	$This->{AssignedAtomTypesCount}{$AtomType} += 1;
      }
    }

    # Assign phramacophore types to atom...
    $AtomID = $Atom->GetID();
    $This->{AssignedAtomTypes}{$AtomID} = \@AtomTypes;
  }
  return $This;
}

# Initialize values of all possible pharmacohore atom pairs...
#
# Let:
#   Dmin = Minimum distance correspoding to number of bonds between two atoms
#   Dmax = Maximum distance correspoding to number of bonds between two atoms
#   D = Distance correspoding to number of bonds between two atoms
#
#   P = Number of pharmacophore atom types to consider
#   PPDn = Number of possible unique pharmacophore atom pairs at a distance Dn
#
#   PPT = Total number of possible pharmacophore atom pairs at all distances between Dmin and Dmax
#
# Then:
#
#   PPD =  (P * (P - 1))/2 + P
#
#   PPT = ((Dmax - Dmin) + 1) * ((P * (P - 1))/2 + P)
#       = ((Dmax - Dmin) + 1) * PPD
#
#
# So for default values of Dmin = 1, Dmax = 10 and P = 5,
#
#   PPD =  (5 * (5 - 1))/2 + 5 = 15
#   PPT = ((10 - 1) + 1) * 15 = 150
#
# the pharmacophore atom pairs bais set includes 150 values.
#
sub _InitializePharmacophoreAtomPairs {
  my($This) = @_;
  my($Distance, $Index1, $Index2, $AtomType1, $AtomType2);

  %{$This->{AtomPairsCount}} = ();

  for $Distance ($This->{MinDistance} .. $This->{MaxDistance}) {
    %{$This->{AtomPairsCount}{$Distance}} = ();

    for $Index1 (0 .. $#{$This->{AtomTypesToUse}}) {
      $AtomType1 = $This->{AtomTypesToUse}[$Index1];
      %{$This->{AtomPairsCount}{$Distance}{$AtomType1}} = ();

      for $Index2 ($Index1 .. $#{$This->{AtomTypesToUse}}) {
	$AtomType2 = $This->{AtomTypesToUse}[$Index2];
	$This->{AtomPairsCount}{$Distance}{$AtomType1}{$AtomType2} = 0;
      }
    }
  }
  return $This;
}

# Count pharmacophore atom pairs between mininum and maximum distance at each
# distance using distance matrix and pharmacophore atom types assiged to each heavy
# atom.
#
# Let:
#   Px = Pharmacophore atom type x
#   Py = Pharmacophore atom type y
#   Dn = Distance between Px and Py in specified distance range
#
# Then:
#   Px-Dn-Py = Pharmacophore atom pair ID for atom types Px and Py at distance Dn
#
# For example: H-D1-H, H-D2-HBA, PI-D5-PI and so on
#
# Notes:
#   . The row and column indices of distance matrix correspond to atom indices.
#   . Distance value of BigNumber implies the atom is not connected to any other atom.
#   . Due to symmetric nature of distance matrix, only upper or lower triangular matrix
#     needs to be processed during identification and count of pharmacophore atom pairs.
#
sub _CountPharmacohoreAtomPairs {
  my($This) = @_;
  my($NumOfRows, $NumOfCols, $RowIndex, $ColIndex, $DistanceMatrix, $Distance, $AtomID1, $AtomID2, $AtomType1, $AtomType2, $SkipIndexCheck, $CountIncrement);

  $DistanceMatrix = $This->{DistanceMatrix};
  ($NumOfRows, $NumOfCols) = $DistanceMatrix->GetSize();
  $SkipIndexCheck = 0;

  ROWINDEX: for $RowIndex (0 .. ($NumOfRows - 1) ) {
    $AtomID1 = $This->{AtomIndexToID}{$RowIndex};
    if ( !((exists($This->{AssignedAtomTypes}{$AtomID1}) && @{$This->{AssignedAtomTypes}{$AtomID1}})) ) {
      next ROWINDEX;
    }

    COLINDEX: for $ColIndex ($RowIndex .. ($NumOfCols - 1) ) {
      $AtomID2 = $This->{AtomIndexToID}{$ColIndex};
      if ( !((exists($This->{AssignedAtomTypes}{$AtomID2}) && @{$This->{AssignedAtomTypes}{$AtomID2}})) ) {
	next COLINDEX;
      }

      $Distance = $DistanceMatrix->GetValue($RowIndex, $ColIndex, $SkipIndexCheck);
      if ($Distance < $This->{MinDistance} || $Distance > $This->{MaxDistance}) {
	next COLINDEX;
      }

      ATOMTYPE1: for $AtomType1 (@{$This->{AssignedAtomTypes}{$AtomID1}}) {
	if ($This->{AtomTypesWeight}{$AtomType1} == 0) {
	  next ATOMTYPE1;
	}
	ATOMTYPE2: for $AtomType2 (@{$This->{AssignedAtomTypes}{$AtomID2}}) {
	  if ($This->{AtomTypesWeight}{$AtomType2} == 0) {
	    next ATOMTYPE2;
	  }
	  $CountIncrement = $This->{AtomTypesWeight}{$AtomType1} * $This->{AtomTypesWeight}{$AtomType2};
	  if ($AtomType1 le $AtomType2) {
	    $This->{AtomPairsCount}{$Distance}{$AtomType1}{$AtomType2} += $CountIncrement;
	  }
	  else {
	    $This->{AtomPairsCount}{$Distance}{$AtomType2}{$AtomType1} += $CountIncrement;
	  }
	}
      }
    }
  }
  return $This;
}

# Normalize the occurance count of pharmacophore atom pairs over the specified distance
# range...
#
sub _NormalizePharmacohoreAtomPairsCount {
  my($This) = @_;

  METHODOLOGY: {
    if ($This->{NormalizationMethodology} =~ /^None$/i) {
      last METHODOLOGY;
    }
    if ($This->{NormalizationMethodology} =~ /^ByHeavyAtomsCount$/i) {
      $This->_NormalizeAtomPairsCountByHeavyAtomsCount();
      last METHODOLOGY;
    }
    if ($This->{NormalizationMethodology} =~ /^ByAtomTypesCount$/i) {
      $This->_NormalizeAtomPairsCountByAtomTypesCount();
      last METHODOLOGY;
    }
    croak "Error: ${ClassName}->_NormalizePharmacohoreAtomPairsCount: Unknown NormalizationMethodology: $This->{NormalizationMethodology}...";
  }
  return $This;
}


# Normalize the occurance count of pharmacophore atom pairs at various distances by
# heavy atom count...
#
sub _NormalizeAtomPairsCountByHeavyAtomsCount {
  my($This) = @_;
  my($Distance, $AtomType1, $AtomType2);

  if ($This->{HeavyAtomCount} == 0) {
    return $This;
  }

  for $Distance (keys %{$This->{AtomPairsCount}} ) {
    for $AtomType1 (keys %{$This->{AtomPairsCount}{$Distance}} ) {
      ATOMTYPE2: for $AtomType2 (keys %{$This->{AtomPairsCount}{$Distance}{$AtomType1}} ) {
	if ($This->{AtomPairsCount}{$Distance}{$AtomType1}{$AtomType2} == 0) {
	  next ATOMTYPE2;
	}
	$This->{AtomPairsCount}{$Distance}{$AtomType1}{$AtomType2} /= $This->{HeavyAtomCount};
      }
    }
  }
  return $This;
}

# Normalize the occurance count of pharmacophore atom pairs at various distances by
# dividing it using sum of the count of each pharmacophore atom type present in the
# molecule for the corresponding atom pair.
#
sub _NormalizeAtomPairsCountByAtomTypesCount {
  my($This) = @_;
  my($Distance, $AtomType1, $AtomType2, $AtomType1Count, $AtomType2Count, $NormalizationFactor);

  for $Distance (keys %{$This->{AtomPairsCount}} ) {
    for $AtomType1 (keys %{$This->{AtomPairsCount}{$Distance}} ) {
      ATOMTYPE2: for $AtomType2 (keys %{$This->{AtomPairsCount}{$Distance}{$AtomType1}} ) {
	if ($This->{AtomPairsCount}{$Distance}{$AtomType1}{$AtomType2} == 0) {
	  next ATOMTYPE2;
	}
	$NormalizationFactor = $This->{AssignedAtomTypesCount}{$AtomType1} + $This->{AssignedAtomTypesCount}{$AtomType2};
	$This->{AtomPairsCount}{$Distance}{$AtomType1}{$AtomType2} /= $NormalizationFactor;
      }
    }
  }
  return $This;
}

# Fuzzify pharmacophore atom pairs count...
#
# Let:
#   Px = Pharmacophore atom type x
#   Py = Pharmacophore atom type y
#
#   PPxy = Pharmacophore atom pair between atom type Px and Py
#
#   PPxyDn = Pharmacophore atom pairs count between atom type Px and Py at distance Dn
#   PPxyDn-1 = Pharmacophore atom pairs count between atom type Px and Py at distance Dn - 1
#   PPxyDn+1 = Pharmacophore atom pairs count between atom type Px and Py at distance Dn + 1
#
#   FF = FuzzFactor for FuzzyBinning and FuzzyBinSmoothing
#
# Then:
#
# For FuzzyBinning:
#
#   PPxyDn = PPxyDn (Unchanged)
#
#   PPxyDn-1 = PPxyDn-1 + PPxyDn * FF
#   PPxyDn+1 = PPxyDn+1 + PPxyDn * FF
#
# For FuzzyBinSmoothing:
#
#   PPxyDn = PPxyDn - PPxyDn * 2FF for Dmin < Dn < Dmax
#   PPxyDn = PPxyDn - PPxyDn * FF for Dn = Dmin or Dmax
#
#   PPxyDn-1 = PPxyDn-1 + PPxyDn * FF
#   PPxyDn+1 = PPxyDn+1 + PPxyDn * FF
#
# In both fuzzification schemes, a value of 0 for FF implies no fuzzification of occurance counts.
# A value of 1 during FuzzyBinning corresponds to maximum fuzzification of occurance counts;
# however, a value of 1 during FuzzyBinSmoothing ends up completely distributing the value over
# the previous and next distance bins.
#
# So for default value of FuzzFactor (FF) 0.15, the occurance count of pharmacohore atom pairs
# at distance Dn during FuzzyBinning is left unchanged and the counts at distances Dn -1 and Dn + 1
# are incremened by PPxyDn * 0.15.
#
# And during FuzzyBinSmoothing the occurance counts at Distance Dn is scaled back using multiplicate
# factor of (1 - 2*0.15) and the occurance counts at distances Dn -1 and Dn + 1 are incremened by
# PPxyDn * 0.15. In otherwords, occurance bin count is smoothed out by distributing it over the
# previous and next distance value.
#
sub _FuzzifyPharmacohoreAtomPairsCount {
  my($This) = @_;
  my($Index1, $Index2, $AtomType1, $AtomType2, $CurrentDistance, $CurrentCount, $NextDistance, $NextCount, $PreviousDistance, $ModifyCurrentCount, $ChangeInCountValue);

  if (!($This->{FuzzifyAtomPairsCount} && $This->{FuzzFactor} > 0)) {
    return $This;
  }

  $ModifyCurrentCount = ($This->{FuzzificationMethodology} =~ /^FuzzyBinSmoothing$/i) ? 1 : 0;

  for $Index1 (0 .. $#{$This->{AtomTypesToUse}}) {
    $AtomType1 = $This->{AtomTypesToUse}[$Index1];
    for $Index2 ($Index1 .. $#{$This->{AtomTypesToUse}}) {
      $AtomType2 = $This->{AtomTypesToUse}[$Index2];

      $CurrentCount = 0; $NextCount = 0;

      $NextDistance = $This->{MinDistance};
      $NextCount = $This->{AtomPairsCount}{$NextDistance}{$AtomType1}{$AtomType2};

      DISTANCE: for $CurrentDistance ($This->{MinDistance} .. $This->{MaxDistance}) {
	$NextDistance = $CurrentDistance + 1;
	$PreviousDistance = $CurrentDistance - 1;

	$CurrentCount = $NextCount;
	$NextCount = ($CurrentDistance < $This->{MaxDistance}) ? $This->{AtomPairsCount}{$NextDistance}{$AtomType1}{$AtomType2} : 0;

	if ($CurrentCount == 0) {
	  # No contribution to fuzzy binning from this distance...
	  next DISTANCE;
	}

	$ChangeInCountValue = $CurrentCount * $This->{FuzzFactor};

	if ($CurrentDistance > $This->{MinDistance}) {
	  # Increment count at previous distance...
	  $This->{AtomPairsCount}{$PreviousDistance}{$AtomType1}{$AtomType2} += $ChangeInCountValue;
	}

	if ($ModifyCurrentCount) {
	  # Decrement count at current distance for FuzzyBinSmoothing...
	  if ($CurrentDistance > $This->{MinDistance} && $CurrentDistance < $This->{MaxDistance}) {
	    $This->{AtomPairsCount}{$CurrentDistance}{$AtomType1}{$AtomType2} -= 2 * $ChangeInCountValue;
	  }
	  else {
	    $This->{AtomPairsCount}{$CurrentDistance}{$AtomType1}{$AtomType2} -= $ChangeInCountValue;
	  }
	}

	if ($CurrentDistance < $This->{MaxDistance}) {
	  # Increment count at next distance...
	  $This->{AtomPairsCount}{$NextDistance}{$AtomType1}{$AtomType2} += $ChangeInCountValue;
	}
      }
    }
  }
  return $This;
}

# Set final fingerpritns vector...
#
sub _SetFinalFingerprints {
  my($This) = @_;
  my($Distance, $Index1, $Index2, $AtomType1, $AtomType2, $Value, $RoundOffValues, $ValuesPrecision, $UseArbitrarySetSize, @Values);

  # Mark successful generation of fingerprints...
  $This->{FingerprintsGenerated} = 1;

  @Values = ();
  @{$This->{AtomPairsIDs}} = ();

  # Do values need to be rounded off?
  $RoundOffValues = (($This->{NormalizationMethodology} !~ /^None$/i) || ($This->{FuzzifyAtomPairsCount})) ? 1 : 0;
  $ValuesPrecision = $This->{ValuesPrecision};

  # Is it an ArbitraySize atom pairs set size?
  $UseArbitrarySetSize = $This->{AtomPairsSetSizeToUse} =~ /^ArbitrarySize$/i ? 1 : 0;

  # Collect all atom paris count values...
  for $Distance ($This->{MinDistance} .. $This->{MaxDistance}) {
    for $Index1 (0 .. $#{$This->{AtomTypesToUse}}) {
      $AtomType1 = $This->{AtomTypesToUse}[$Index1];
      INDEX2: for $Index2 ($Index1 .. $#{$This->{AtomTypesToUse}}) {
	$AtomType2 = $This->{AtomTypesToUse}[$Index2];

	# Atom pair count...
	$Value = $This->{AtomPairsCount}{$Distance}{$AtomType1}{$AtomType2};
	if ($RoundOffValues) {
	  $Value = MathUtil::round($Value, $This->{ValuesPrecision}) + 0;
	}

	# Ignore or not to ignore...
	if ($UseArbitrarySetSize && $Value == 0) {
	  next INDEX2;
	}

	push @{$This->{AtomPairsIDs}}, "${AtomType1}-D${Distance}-${AtomType2}";
	push @Values, $Value;
      }
    }
  }

  # Add AtomPairsIDs and count values to fingerprint vector...
  $This->{FingerprintsVector}->AddValueIDs(\@{$This->{AtomPairsIDs}});
  $This->{FingerprintsVector}->AddValues(\@Values);

  return $This;
}

# Get pharmacophore atom pair IDs corresponding to atom pairs count values in
# fingerprint vector as an array or reference to an array...
#
# AtomPairIDs list  is generated during finalization  of fingerprints  and the fingerprint
# vector containing count values matches the atom pairs array.
#
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


# Return a string containg data for TopologicalPharmacophoreAtomPairsFingerprints object...
sub StringifyTopologicalPharmacophoreAtomPairsFingerprints {
  my($This) = @_;
  my($FingerprintsString);

  # Type of fingerprint...
  $FingerprintsString = "Fingerprint type: $This->{Type}; AtomPairsSetSizeToUse: $This->{AtomPairsSetSizeToUse}";

  # Min and max distance...
  $FingerprintsString .= "; MinDistance:  $This->{MinDistance}; MaxDistance: $This->{MaxDistance}";

  # Pharmacophore type labels and description...
  my($AtomType, @AtomTypes, @AtomTypesOrder, %AvailableAtomTypes);

  @AtomTypesOrder = AtomTypes::FunctionalClassAtomTypes::GetFunctionalClassesOrder();
  %AvailableAtomTypes = AtomTypes::FunctionalClassAtomTypes::GetAvailableFunctionalClasses();

  @AtomTypes = ();
  for $AtomType (@AtomTypesOrder) {
    push @AtomTypes, "$AtomType: $AvailableAtomTypes{$AtomType}";
  }

  $FingerprintsString .= "; AtomTypesToUse: <" . TextUtil::JoinWords(\@{$This->{AtomTypesToUse}}, ", ", 0) . ">";
  $FingerprintsString .= "; AtomTypesOrder: <" . TextUtil::JoinWords(\@AtomTypesOrder, ", ", 0) . ">";
  $FingerprintsString .= "; AvailableAtomTypes: <" . TextUtil::JoinWords(\@AtomTypes, ", ", 0) . ">";

  # Normalization method...
  $FingerprintsString .= "; NormalizationMethodology: $This->{NormalizationMethodology}";

  # Weights...
  my($FirstLabel, $Label, $Weight);

  $FingerprintsString .= "; AtomTypesWeight <Labels: Weight>: <";
  $FirstLabel = 1;
  for $Label (sort @{$This->{AtomTypesToUse}}) {
    $Weight = $This->{AtomTypesWeight}{$Label};
    if ($FirstLabel) {
      $FirstLabel = 0;
      $FingerprintsString .= " ${Label}: ${Weight}";
    }
    else {
      $FingerprintsString .= "; ${Label}: ${Weight}";
    }
  }
  $FingerprintsString .= ">";

  # Fuzzification of count...
  my($FuzzifyFlag);
  $FuzzifyFlag = $This->{FuzzifyAtomPairsCount} ? "Yes" : "No";
  $FingerprintsString .= "; FuzzifyAtomPairsCount: $FuzzifyFlag; FuzzificationMode: $This->{FuzzificationMode}; FuzzificationMethodology: $This->{FuzzificationMethodology}; FuzzFactor: $This->{FuzzFactor}";

  # Total number of pharmacophore atom pairs...
  $FingerprintsString .= "; NumOfAtomPairs: " . $This->{FingerprintsVector}->GetNumOfValues();

  # FingerprintsVector...
  $FingerprintsString .= "; FingerprintsVector: < $This->{FingerprintsVector} >";

  return $FingerprintsString;
}

1;

__END__

=head1 NAME

TopologicalPharmacophoreAtomPairsFingerprints

=head1 SYNOPSIS

use Fingerprints::TopologicalPharmacophoreAtomPairsFingerprints;

use Fingerprints::TopologicalPharmacophoreAtomPairsFingerprints qw(:all);

=head1 DESCRIPTION

B<TopologicalPharmacophoreAtomPairsFingerprints> [ Ref 60-62, Ref 65, Ref 68 ] class provides
the following methods:

new, GenerateFingerprints, GetDescription, GetAtomPairIDs, SetAtomTypesToUse,
SetAtomTypesWeight, SetFuzzFactor, SetFuzzificationMethodology,
SetFuzzificationMode, SetMaxDistance, SetMinDistance,
SetNormalizationMethodology, SetValuesPrecision,
StringifyTopologicalPharmacophoreAtomPairsFingerprints

B<TopologicalPharmacophoreAtomPairsFingerprints> is derived from B<Fingerprints> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<TopologicalPharmacophoreAtomPairsFingerprints>, B<Fingerprints> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

Based on the values specified for B<AtomTypesToUse>, pharmacophore atom types are
assigned to all non-hydrogen atoms in a molecule and a distance matrix is generated.
A pharmacophore atom pairs basis set is initialized for all unique possible pairs within
B<MinDistance> and B<MaxDistance> range.

    Let:

    P = Valid pharmacophore atom type

    Px = Pharmacophore atom type x
    Py = Pharmacophore atom type y

    Dmin = Minimum distance corresponding to number of bonds between two atoms
    Dmax = Maximum distance corresponding to number of bonds between two atoms
    D = Distance corresponding to number of bonds between two atoms

    Px-Dn-Py = Pharmacophore atom pair ID for atom types Px and Py at distance Dn

    P = Number of pharmacophore atom types to consider
    PPDn = Number of possible unique pharmacophore atom pairs at a distance Dn

    PPT = Total number of possible pharmacophore atom pairs at all distances between Dmin and Dmax

    Then:

    PPD =  (P * (P - 1))/2 + P

    PPT = ((Dmax - Dmin) + 1) * ((P * (P - 1))/2 + P)
        = ((Dmax - Dmin) + 1) * PPD

    So for default values of Dmin = 1, Dmax = 10 and P = 5,

    PPD =  (5 * (5 - 1))/2 + 5 = 15
    PPT = ((10 - 1) + 1) * 15 = 150

    The pharmacophore atom pairs bais set includes 150 values.

    The atom pair IDs correspond to:

    Px-Dn-Py = Pharmacophore atom pair ID for atom types Px and Py at distance Dn

    For example: H-D1-H, H-D2-HBA, PI-D5-PI and so on

Using distance matrix and pharmacohore atom types, occurrence of unique pharmacohore atom
pairs is counted. The contribution of each atom type to atom pair interaction is optionally
weighted by specified B<AtomTypesWeight> before assigning its count to appropriate distance
bin. Based on B<NormalizationMethodology> option, pharmacophore atom pairs count is optionally
normalized. Additionally, pharmacohore atom pairs count is optionally fuzzified before or after
the normalization controlled by values of B<FuzzifyAtomPairsCount>, B<FuzzificationMode>,
B<FuzzificationMethodology> and B<FuzzFactor>.

The final pharmacophore atom pairs count along with atom pair identifiers involving all non-hydrogen
atoms, with optional normalization and fuzzification, constitute pharmacophore topological atom pairs
fingerprints of the molecule.

For I<ArbitrarySize> value of B<AtomPairsSetSizeToUse>, the fingerprint vector correspond to
only those topological pharmacophore atom pairs which are present and have non-zero count. However,
for I<FixedSize> value of B<AtomPairsSetSizeToUse>, the fingerprint vector contains all possible
valid topological pharmacophore atom pairs with both zero and non-zero count values.

The current release of MayaChemTools generates the following types of topological pharmacophore
atom pairs fingerprints vector strings:

    FingerprintsVector;TopologicalPharmacophoreAtomPairs:ArbitrarySize:Min
    Distance1:MaxDistance10;54;NumericalValues;IDsAndValuesString;H-D1-H H
    -D1-NI HBA-D1-NI HBD-D1-NI H-D2-H H-D2-HBA H-D2-HBD HBA-D2-HBA HBA-D2-
    HBD H-D3-H H-D3-HBA H-D3-HBD H-D3-NI HBA-D3-NI HBD-D3-NI H-D4-H H-D4-H
    BA H-D4-HBD HBA-D4-HBA HBA-D4-HBD HBD-D4-HBD H-D5-H H-D5-HBA H-D5-...;
    18 1 2 1 22 12 8 1 2 18 6 3 1 1 1 22 13 6 5 7 2 28 9 5 1 1 1 36 16 10
    3 4 1 37 10 8 1 35 10 9 3 3 1 28 7 7 4 18 16 12 5 1 2 1

    FingerprintsVector;TopologicalPharmacophoreAtomPairs:FixedSize:MinDist
    ance1:MaxDistance10;150;OrderedNumericalValues;ValuesString;18 0 0 1 0
    0 0 2 0 0 1 0 0 0 0 22 12 8 0 0 1 2 0 0 0 0 0 0 0 0 18 6 3 1 0 0 0 1
    0 0 1 0 0 0 0 22 13 6 0 0 5 7 0 0 2 0 0 0 0 0 28 9 5 1 0 0 0 1 0 0 1 0
    0 0 0 36 16 10 0 0 3 4 0 0 1 0 0 0 0 0 37 10 8 0 0 0 0 1 0 0 0 0 0 0
    0 35 10 9 0 0 3 3 0 0 1 0 0 0 0 0 28 7 7 4 0 0 0 0 0 0 0 0 0 0 0 18...

    FingerprintsVector;TopologicalPharmacophoreAtomPairs:FixedSize:MinDist
    ance1:MaxDistance10;150;OrderedNumericalValues;IDsAndValuesString;H-D1
    -H H-D1-HBA H-D1-HBD H-D1-NI H-D1-PI HBA-D1-HBA HBA-D1-HBD HBA-D1-NI H
    BA-D1-PI HBD-D1-HBD HBD-D1-NI HBD-D1-PI NI-D1-NI NI-D1-PI PI-D1-PI H-D
    2-H H-D2-HBA H-D2-HBD H-D2-NI H-D2-PI HBA-D2-HBA HBA-D2-HBD HBA-D2...;
    18 0 0 1 0 0 0 2 0 0 1 0 0 0 0 22 12 8 0 0 1 2 0 0 0 0 0 0 0 0 18 6 3
    1 0 0 0 1 0 0 1 0 0 0 0 22 13 6 0 0 5 7 0 0 2 0 0 0 0 0 28 9 5 1 0 0 0
    1 0 0 1 0 0 0 0 36 16 10 0 0 3 4 0 0 1 0 0 0 0

=head2 METHODS

=over 4

=item B<new>

    $TPAPFP = new TopologicalPharmacophoreAtomPairsFingerprints(
                                                   %NamesAndValues);

Using specified I<TopologicalPharmacophoreAtomPairsFingerprints> property names and
values hash, B<new> method creates a new object and returns a reference to newly created
B<TopologicalPharmacophoreAtomPairsFingerprints> object. By default, the following properties
are initialized:

    Molecule = ''
    Type = 'TopologicalPharmacophoreAtomPairs'
    MinDistance = 1
    MaxDistance = 10
    NormalizationMethodology = 'None'
    AtomTypesToUse = ['HBD', 'HBA', 'PI', 'NI', 'H']

    FuzzifyAtomPairsCount = 0
    FuzzificationMode = 'AfterNormalization'
    FuzzificationMethodology =  'FuzzyBinning'
    FuzzFactor = 0.15

    ValuesPrecision = 2

Examples:

    $TPAPFP = new TopologicalPharmacophoreAtomPairsFingerprints(
                              'Molecule' => $Molecule);

    $TPAPFP = new TopologicalPharmacophoreAtomPairsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomPairsSetSizeToUse' => 'ArbitrarySize',
                              'MinDistance' => 1,
                              'MaxDistance' => 10,
                              'NormalizationMethodology' => 'None',
                              'AtomTypesToUse' => ['HBD', 'HBA', 'PI', 'NI', 'H'],
                              'FuzzifyAtomPairsCount' => 0);

    $TPAPFP = new TopologicalPharmacophoreAtomPairsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomPairsSetSizeToUse' => 'FizedSize',
                              'MinDistance' => 1,
                              'MaxDistance' => 10,
                              'NormalizationMethodology' => 'None',
                              'AtomTypesToUse' => ['HBD', 'HBA', 'PI', 'NI', 'H'],
                              'FuzzifyAtomPairsCount' => 1,
                              'FuzzificationMethodology' => 'FuzzyBinning',
                              'FuzzFactor' => 0.15,
                              'ValuesPrecision' => 2);

    $TPAPFP->GenerateFingerprints();
    print "$TPAPFP\n";

=item B<GetDescription>

    $Description = $TopologicalPharmacophoreAtomPairsFP->GetDescription();

Returns a string containing description of topological pharmacophore atom pairs fingerprints.

=item B<GenerateFingerprints>

    $TopologicalPharmacophoreAtomPairsFP->GenerateFingerprints();

Generates topological pharmacophore atom pairs fingerprints and returns
I<TopologicalPharmacophoreAtomPairsFP>.

=item B<GetAtomPairIDs>

    $AtomPairIDsRef = $TopologicalPharmacophoreAtomPairsFP->GetAtomPairIDs();
    @AtomPairIDs = $TopologicalPharmacophoreAtomPairsFP->GetAtomPairIDs();

Returns atom pair IDs corresponding to atom pairs count values in topological pharmacophore
atom pairs fingerprints vector as an array or reference to an array.

=item B<SetAtomPairsSetSizeToUse>

    $TopologicalPharmacophoreAtomPairsFP->SetAtomPairsSetSizeToUse($Values);

Sets pharmacophore atom pairs set size to use for topological pharmacophore fingerprints
generation and returns I<TopologicalPharmacophoreAtomPairsFingerprints>.

Possible values for pharmacophore atom pairs set size are: I<ArbitrarySize, FizedSize>.
Default value: I<ArbitrarySize>.

For I<ArbitrarySize> value of B<AtomPairsSetSizeToUse>, the fingerprint vector correspond to
only those topological pharmacophore atom pairs which are present and have non-zero count. However,
for I<FixedSize> value of B<AtomPairsSetSizeToUse>, the fingerprint vector contains all possible
valid topological pharmacophore atom pairs with both zero and non-zero count values.

=item B<SetAtomTypesToUse>

    $TopologicalPharmacophoreAtomPairsFP->SetAtomTypesToUse($ValuesRef);
    $TopologicalPharmacophoreAtomPairsFP->SetAtomTypesToUse(@Values);

Sets pharmacophore atom types to use for topological pharmacophore fingerprints
generation and returns I<TopologicalPharmacophoreAtomPairsFingerprints>.

Possible values for pharmacophore atom types are: I<Ar, CA, H, HBA, HBD, Hal, NI, PI, RA>.
Default value [ Ref 60-62 ] : I<HBD,HBA,PI,NI,H>.

The pharmacophore atom types abbreviations correspond to:

    HBD: HydrogenBondDonor
    HBA: HydrogenBondAcceptor
    PI :  PositivelyIonizable
    NI : NegativelyIonizable
    Ar : Aromatic
    Hal : Halogen
    H : Hydrophobic
    RA : RingAtom
    CA : ChainAtom

I<AtomTypes::FunctionalClassAtomTypes> module is used to assign pharmacophore atom
types. It uses following definitions [ Ref 60-61, Ref 65-66 ]:

    HydrogenBondDonor: NH, NH2, OH
    HydrogenBondAcceptor: N[!H], O
    PositivelyIonizable: +, NH2
    NegativelyIonizable: -, C(=O)OH, S(=O)OH, P(=O)OH

=item B<SetAtomTypesWeight>

    $TopologicalPharmacophoreAtomPairsFP->SetAtomTypesWeight(
        %AtomTypesToWeight);

Sets weights of specified pharmacophore atom types to use during calculation of their contribution
to atom pair count and returns I<TopologicalPharmacophoreAtomPairsFP>.  Default values: I<1 for
each atom type>.

The weight values allow to increase the importance of specific pharmacophore atom type
in the generated fingerprints. A weight value of 0 for an atom type eliminates its contribution to
atom pair count where as weight value of 2 doubles its contribution.

=item B<SetFuzzFactor>

    $TopologicalPharmacophoreAtomPairsFP->SetFuzzFactor($Value);

Sets fuzz factor value to use during fuzzification of atom pairs count and returns
I<TopologicalPharmacophoreAtomPairsFP>. Default value: I<0.15>.

Valid values: For I<FuzzyBinning> value of B<FuzzificationMethodology>: I<between 0 and 1.0>; For
I<FuzzyBinSmoothing> value of B<FuzzificationMethodology>: I<between 0 and 0.5>.

=item B<SetFuzzificationMethodology>

    $TopologicalPharmacophoreAtomPairsFP->SetFuzzificationMethodology($Value);

Sets fuzzification methodology to use for fuzzification of atom pairs count and returns
I<TopologicalPharmacophoreAtomPairsFP>. Default value: I<FuzzyBinning>.  Possible values:
I<FuzzyBinning | FuzzyBinSmoothing>.

In conjunction with values for options B<FuzzifyAtomPairsCount>, B<FuzzificationMode> and
B<FuzzFactor>, B<FuzzificationMethodology> option is used to fuzzify pharmacophore atom
pairs count.

Let:

    Px = Pharmacophore atom type x
    Py = Pharmacophore atom type y
    PPxy = Pharmacophore atom pair between atom type Px and Py

    PPxyDn = Pharmacophore atom pairs count between atom type Px and Py
             at distance Dn
    PPxyDn-1 = Pharmacophore atom pairs count between atom type Px and Py
               at distance Dn - 1
    PPxyDn+1 = Pharmacophore atom pairs count between atom type Px and Py
               at distance Dn + 1

    FF = FuzzFactor for FuzzyBinning and FuzzyBinSmoothing

Then:

For I<FuzzyBinning>:

    PPxyDn = PPxyDn (Unchanged)

    PPxyDn-1 = PPxyDn-1 + PPxyDn * FF
    PPxyDn+1 = PPxyDn+1 + PPxyDn * FF

For I<FuzzyBinSmoothing>:

    PPxyDn = PPxyDn - PPxyDn * 2FF for Dmin < Dn < Dmax
    PPxyDn = PPxyDn - PPxyDn * FF for Dn = Dmin or Dmax

    PPxyDn-1 = PPxyDn-1 + PPxyDn * FF
    PPxyDn+1 = PPxyDn+1 + PPxyDn * FF

In both fuzzification schemes, a value of 0 for FF implies no fuzzification of occurrence counts.
A value of 1 during I<FuzzyBinning> corresponds to maximum fuzzification of occurrence counts;
however, a value of 1 during I<FuzzyBinSmoothing> ends up completely distributing the value over
the previous and next distance bins.

So for default value of B<FuzzFactor> (FF) 0.15, the occurrence count of pharmacohore atom pairs
at distance Dn during FuzzyBinning is left unchanged and the counts at distances Dn -1 and Dn + 1
are incremented by PPxyDn * 0.15.

And during I<FuzzyBinSmoothing> the occurrence counts at Distance Dn is scaled back using multiplicative
factor of (1 - 2*0.15) and the occurrence counts at distances Dn -1 and Dn + 1 are incremented by
PPxyDn * 0.15. In other words, occurrence bin count is smoothed out by distributing it over the
previous and next distance value.

=item B<SetFuzzificationMode>

    $TopologicalPharmacophoreAtomPairsFP->SetFuzzificationMode($Value);

Sets fuzzification mode to use for fuzzification of atom pairs count and returns
I<TopologicalPharmacophoreAtomPairsFP>. Default value: I<AfterNormalization>.  Possible values:
I<BeforeNormalization | AfterNormalization>.

=item B<SetMaxDistance>

    $TopologicalPharmacophoreAtomPairsFP->SetMaxDistance($Value);

Sets maximum bond distance between atom pairs for generating topological pharmacophore atom
pairs fingerprints and returns I<TopologicalPharmacophoreAtomPairsFP>.

=item B<SetMinDistance>

    $TopologicalPharmacophoreAtomPairsFP->SetMinDistance($Value);

Sets minimum bond distance between atom pairs for generating topological pharmacophore atom
pairs fingerprints and returns I<TopologicalPharmacophoreAtomPairsFP>.

=item B<SetNormalizationMethodology>

    $TopologicalPharmacophoreAtomPairsFP->SetNormalizationMethodology($Value);

Sets normalization methodology to use for scaling the occurrence count of pharmacophore atom
pairs within specified distance range and returns I<TopologicalPharmacophoreAtomPairsFP>.
Default value: I<None>. Possible values: I<None, ByHeavyAtomsCount or ByAtomTypesCount>.

=item B<SetValuesPrecision>

    $TopologicalPharmacophoreAtomPairsFP->SetValuesPrecision($Value);

Sets precision of atom pairs count real values which might be generated after normalization
or fuzzification  and returns I<TopologicalPharmacophoreAtomPairsFP>. Default: up to I<2> decimal
places.

=item B<StringifyTopologicalPharmacophoreAtomPairsFingerprints>

    $String = $TopologicalPharmacophoreAtomPairsFP->
                  StringifyTopologicalPharmacophoreAtomPairsFingerprints();

Returns a string containing information about I<TopologicalPharmacophoreAtomPairsFingerprints> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Fingerprints.pm, FingerprintsStringUtil.pm, AtomNeighborhoodsFingerprints.pm,
AtomTypesFingerprints.pm, EStateIndiciesFingerprints.pm, ExtendedConnectivityFingerprints.pm,
MACCSKeys.pm, PathLengthFingerprints.pm, TopologicalAtomPairsFingerprints.pm,
TopologicalAtomTripletsFingerprints.pm, TopologicalAtomTorsionsFingerprints.pm,
TopologicalPharmacophoreAtomTripletsFingerprints.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
