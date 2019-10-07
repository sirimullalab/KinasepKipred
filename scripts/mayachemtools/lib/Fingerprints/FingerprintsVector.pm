package Fingerprints::FingerprintsVector;
#
# File: FingerprintsVector.pm
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
use MathUtil ();
use TextUtil ();
use StatisticsUtil ();
use BitVector;
use Vector;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);

# Distance coefficients
my(@DistanceCoefficients) = qw(CityBlockDistanceCoefficient EuclideanDistanceCoefficient HammingDistanceCoefficient ManhattanDistanceCoefficient SoergelDistanceCoefficient);

# Similarity coefficients...
my(@SimilarityCoefficients) = qw(CosineSimilarityCoefficient CzekanowskiSimilarityCoefficient DiceSimilarityCoefficient OchiaiSimilarityCoefficient JaccardSimilarityCoefficient SorensonSimilarityCoefficient TanimotoSimilarityCoefficient);

# New from string...
my(@NewFromString) = qw(NewFromValuesString NewFromValuesAndIDsString NewFromIDsAndValuesString NewFromValuesAndIDsPairsString NewFromIDsAndValuesPairsString);

@EXPORT = qw(IsFingerprintsVector);
@EXPORT_OK = qw(GetSupportedDistanceCoefficients GetSupportedSimilarityCoefficients GetSupportedDistanceAndSimilarityCoefficients @DistanceCoefficients @SimilarityCoefficients);

%EXPORT_TAGS = (
		new => [@NewFromString],
		distancecoefficients => [@DistanceCoefficients],
		similaritycoefficients => [@SimilarityCoefficients],
		all  => [@EXPORT, @EXPORT_OK]
	       );

# Setup class variables...
my($ClassName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyFingerprintsVector';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;

  $This->_InitializeFingerprintsVector();

  $This->_InitializeFingerprintsVectorProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeFingerprintsVector {
  my($This) = @_;

  # Type of fingerprint vector...
  $This->{Type} = '';

  # Fingerprint vector values...
  @{$This->{Values}} = ();

  # Fingerprint vector value IDs...
  @{$This->{ValueIDs}} = ();

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize object properties....
sub _InitializeFingerprintsVectorProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  if (!exists $NamesAndValues{Type}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying type...";
  }
  return $This;
}

# Create a new fingerprints vector using space delimited values string. This functionality can be
# either invoked as a class function or an object method.
#
sub NewFromValuesString ($$;$) {
  my($FirstParameter, $SecondParameter, $ThirdParamater) = @_;
  my($This, $Type, $ValuesString);

  if (@_ == 3) {
    ($This, $Type, $ValuesString) = ($FirstParameter, $SecondParameter, $ThirdParamater);
  }
  else {
    ($This, $Type, $ValuesString) = (undef, $FirstParameter, $SecondParameter);
  }
  my($FingerprintsVector, @Values);

  @Values = ();
  if (defined($ValuesString) && length($ValuesString) && $ValuesString !~ /^None$/i) {
    @Values = split(' ', $ValuesString);
  }

  $FingerprintsVector = new Fingerprints::FingerprintsVector('Type' => $Type, 'Values' => \@Values);

  return $FingerprintsVector;
}

# Create a new fingerprints vector using values and IDs string containing semicolon
# delimited value string and value IDs strings. The values within value and value IDs
# string are delimited by spaces.
#
# This functionality can be either invoked as a class function or an object method.
#
sub NewFromValuesAndIDsString ($$;$) {
  my($FirstParameter, $SecondParameter, $ThirdParamater) = @_;
  my($This, $Type, $ValuesAndIDsString);

  if (@_ == 3) {
    ($This, $Type, $ValuesAndIDsString) = ($FirstParameter, $SecondParameter, $ThirdParamater);
  }
  else {
    ($This, $Type, $ValuesAndIDsString) = (undef, $FirstParameter, $SecondParameter);
  }
  my($FingerprintsVector, $ValuesString, $ValueIDsString, @Values, @ValueIDs);

  ($ValuesString, $ValueIDsString) = split(';', $ValuesAndIDsString);

  @Values = ();
  if (defined($ValuesString) && length($ValuesString) && $ValuesString !~ /^None$/i) {
    @Values = split(' ', $ValuesString);
  }
  @ValueIDs = ();
  if (defined($ValueIDsString) && length($ValueIDsString) && $ValueIDsString !~ /^None$/i) {
    @ValueIDs = split(' ', $ValueIDsString);
  }

  if (@Values != @ValueIDs ) {
    carp "Warning: ${ClassName}->NewFromValuesAndIDsString: Object can't be instantiated: Number specified values, " . scalar @Values . ", must be equal to number of specified value IDs, " . scalar @ValueIDs .  "...";
    return undef;
  }

  $FingerprintsVector = new Fingerprints::FingerprintsVector('Type' => $Type, 'Values' => \@Values, 'ValueIDs' => \@ValueIDs);

  return $FingerprintsVector;
}

# Create a new fingerprints vector using IDs and values string containing semicolon
# delimited value IDs string and values strings. The values within value and value IDs
# string are delimited by spaces.
#
# This functionality can be either invoked as a class function or an object method.
#
sub NewFromIDsAndValuesString ($$;$) {
  my($FirstParameter, $SecondParameter, $ThirdParamater) = @_;
  my($This, $Type, $IDsAndValuesString);

  if (@_ == 3) {
    ($This, $Type, $IDsAndValuesString) = ($FirstParameter, $SecondParameter, $ThirdParamater);
  }
  else {
    ($This, $Type, $IDsAndValuesString) = (undef, $FirstParameter, $SecondParameter);
  }
  my($FingerprintsVector, $ValuesString, $ValueIDsString, @Values, @ValueIDs);

  ($ValueIDsString, $ValuesString) = split(';', $IDsAndValuesString);

  @Values = ();
  if (defined($ValuesString) && length($ValuesString) && $ValuesString !~ /^None$/i) {
    @Values = split(' ', $ValuesString);
  }
  @ValueIDs = ();
  if (defined($ValueIDsString) && length($ValueIDsString) && $ValueIDsString !~ /^None$/i) {
    @ValueIDs = split(' ', $ValueIDsString);
  }

  if (@Values != @ValueIDs ) {
    carp "Warning: ${ClassName}->NewFromIDsAndValuesString: Object can't be instantiated: Number specified values, " . scalar @Values . ", must be equal to number of specified value IDs, " . scalar @ValueIDs .  "...";
    return undef;
  }

  $FingerprintsVector = new Fingerprints::FingerprintsVector('Type' => $Type, 'Values' => \@Values, 'ValueIDs' => \@ValueIDs);

  return $FingerprintsVector;
}

# Create a new fingerprints vector using values and IDs pairs string containing space
# value and value IDs pairs.
#
# This functionality can be either invoked as a class function or an object method.
#
sub NewFromValuesAndIDsPairsString ($$;$) {
  my($FirstParameter, $SecondParameter, $ThirdParamater) = @_;
  my($This, $Type, $ValuesAndIDsPairsString);

  if (@_ == 3) {
    ($This, $Type, $ValuesAndIDsPairsString) = ($FirstParameter, $SecondParameter, $ThirdParamater);
  }
  else {
    ($This, $Type, $ValuesAndIDsPairsString) = (undef, $FirstParameter, $SecondParameter);
  }
  my($FingerprintsVector, $Index, @Values, @ValueIDs, @ValuesAndIDsPairs);

  @ValuesAndIDsPairs = split(' ', $ValuesAndIDsPairsString);
  if (@ValuesAndIDsPairs % 2) {
    carp "Warning: ${ClassName}->NewFromValuesAndIDsPairsString: No fingerprint vector created: Invalid values and IDs pairs data: Input list must contain even number of values and IDs pairs...";
    return undef;
  }

  @Values = (); @ValueIDs = ();
  if (!(@ValuesAndIDsPairs == 2 && $ValuesAndIDsPairs[0] =~ /^None$/i && $ValuesAndIDsPairs[1] =~ /^None$/i)) {
    for ($Index = 0; $Index < $#ValuesAndIDsPairs; $Index += 2) {
      push @Values, $ValuesAndIDsPairs[$Index];
      push @ValueIDs, $ValuesAndIDsPairs[$Index + 1];
    }
  }
  $FingerprintsVector = new Fingerprints::FingerprintsVector('Type' => $Type, 'Values' => \@Values, 'ValueIDs' => \@ValueIDs);

  return $FingerprintsVector;
}

# Create a new fingerprints vector using IDs and values pairs string containing space
# value IDs and valus pairs.
#
# This functionality can be either invoked as a class function or an object method.
#
sub NewFromIDsAndValuesPairsString ($$;$) {
  my($FirstParameter, $SecondParameter, $ThirdParamater) = @_;
  my($This, $Type, $IDsAndValuesPairsString);

  if (@_ == 3) {
    ($This, $Type, $IDsAndValuesPairsString) = ($FirstParameter, $SecondParameter, $ThirdParamater);
  }
  else {
    ($This, $Type, $IDsAndValuesPairsString) = (undef, $FirstParameter, $SecondParameter);
  }
  my($FingerprintsVector, $Index, @Values, @ValueIDs, @IDsAndValuesPairs);

  @IDsAndValuesPairs = split(' ', $IDsAndValuesPairsString);
  if (@IDsAndValuesPairs % 2) {
    croak "Error: ${ClassName}->NewFromIDsAndValuesPairsString: No fingerprint vector created: Invalid values and IDs pairs data: Input list must contain even number of values and IDs pairs...";
    return undef;
  }

  @Values = (); @ValueIDs = ();
  if (!(@IDsAndValuesPairs == 2 && $IDsAndValuesPairs[0] =~ /^None$/i && $IDsAndValuesPairs[1] =~ /^None$/i)) {
    for ($Index = 0; $Index < $#IDsAndValuesPairs; $Index += 2) {
      push @ValueIDs, $IDsAndValuesPairs[$Index];
      push @Values, $IDsAndValuesPairs[$Index + 1];
    }
  }
  $FingerprintsVector = new Fingerprints::FingerprintsVector('Type' => $Type, 'Values' => \@Values, 'ValueIDs' => \@ValueIDs);

  return $FingerprintsVector;
}

# Set type of fingerprint vector. Supported types are: OrderedNumericalValues, NumericalValues, and
# AlphaNumericalValues
#
#  .  For OrderedNumericalValues type, both vectors must be of the same size and contain similar
#     types of numerical values in the same order.
#
#  .  For NumericalValues type, vector value IDs for both vectors must be specified; however, their
#     size and order of IDs and numerical values may be different. For each vector, value IDs must
#     correspond to vector values.
#
#  .  For AlphaNumericalValues type, vectors may contain both numerical and alphanumerical values
#     and their sizes may be different.
#
sub SetType {
  my($This, $Type) = @_;

  if ($Type !~ /^(OrderedNumericalValues|NumericalValues|AlphaNumericalValues)$/i) {
    croak "Error: ${ClassName}->SetType: Specified value, $Type, for Type is not vaild. Supported types in current release of MayaChemTools: OrderedNumericalValues, NumericalValues or AlphaNumericalValues";
  }

  if ($This->{Type}) {
    croak "Error: ${ClassName}->SetType: Can't change intial fingerprints vector type:  It's already set...";
  }
  $This->{Type} = $Type;

  return $This;
}

# Get fingerpints vector type...
#
sub GetType {
  my($This) = @_;

  return $This->{Type};
}

# Set ID...
sub SetID {
  my($This, $Value) = @_;

  $This->{ID} = $Value;

  return $This;
}

# Get ID...
sub GetID {
  my($This) = @_;

  return exists $This->{ID} ? $This->{ID} : 'None';
}

# Set description...
sub SetDescription {
  my($This, $Value) = @_;

  $This->{Description} = $Value;

  return $This;
}

# Get description...
sub GetDescription {
  my($This) = @_;

  return exists $This->{Description} ? $This->{Description} : 'No description available';
}

# Set vector type...
sub SetVectorType {
  my($This, $Value) = @_;

  $This->{VectorType} = $Value;

  return $This;
}

# Get vector type...
sub GetVectorType {
  my($This) = @_;

  return exists $This->{VectorType} ? $This->{VectorType} : 'FingerprintsVector';
}

# Set values of a fingerprint vector using a vector, reference to an array or an array...
#
sub SetValues {
  my($This, @Values) = @_;

  $This->_SetOrAddValuesOrValueIDs("SetValues", @Values);

  return $This;
}

# Set value IDs of a fingerprint vector using a vector, reference to an array or an array...
#
sub SetValueIDs {
  my($This, @Values) = @_;

  $This->_SetOrAddValuesOrValueIDs("SetValueIDs", @Values);

  return $This;
}

# Add values to a fingerprint vector using a vector, reference to an array or an array...
#
sub AddValues {
  my($This, @Values) = @_;

  $This->_SetOrAddValuesOrValueIDs("AddValues", @Values);

  return $This;
}

# Add value IDs to a fingerprint vector using a vector, reference to an array or an array...
#
sub AddValueIDs {
  my($This, @Values) = @_;

  $This->_SetOrAddValuesOrValueIDs("AddValueIDs", @Values);

  return $This;
}

# Set or add values or value IDs using:
#
#    o List of values or ValueIDs
#    o Reference to an list of values or ValuesIDs
#    o A vector containing values or ValueIDs
#
sub _SetOrAddValuesOrValueIDs {
  my($This, $Mode, @Values) = @_;

  if (!@Values) {
    return;
  }

  # Collect specified values or valueIDs...
  my($FirstValue, $TypeOfFirstValue, $ValuesRef);

  $FirstValue = $Values[0];
  $TypeOfFirstValue = ref $FirstValue;
  if ($TypeOfFirstValue =~ /^(SCALAR|HASH|CODE|REF|GLOB)/) {
    croak "Error: ${ClassName}-> _SetOrAddValuesOrValueIDs: Trying to add values to vector object with a reference to unsupported value format...";
  }

  if (Vector::IsVector($FirstValue)) {
    # It's a vector...
    $ValuesRef = $FirstValue->GetValues();
  }
  elsif ($TypeOfFirstValue =~ /^ARRAY/) {
    # It's an array refernce...
    $ValuesRef = $FirstValue;
  }
  else {
    # It's a list of values...
    $ValuesRef = \@Values;
  }

  # Set or add values or value IDs...
  MODE: {
    if ($Mode =~ /^SetValues$/i) { @{$This->{Values}} = (); push @{$This->{Values}}, @{$ValuesRef}; last MODE; }
    if ($Mode =~ /^SetValueIDs$/i) { @{$This->{ValueIDs}} = (); push @{$This->{ValueIDs}}, @{$ValuesRef}; last MODE; }
    if ($Mode =~ /^AddValues$/i) { push @{$This->{Values}}, @{$ValuesRef}; last MODE; }
    if ($Mode =~ /^AddValueIDs$/i) { push @{$This->{ValueIDs}}, @{$ValuesRef}; last MODE; }
    croak "Error: ${ClassName}-> _SetOrAddValuesOrValueIDs: Unknown mode $Mode...";
  }
  return $This;
}

# Set a specific value in fingerprint vector with indicies starting from 0..
#
sub SetValue {
  my($This, $Index, $Value, $SkipCheck) = @_;

  # Just set it...
  if ($SkipCheck) {
    return $This->_SetValue($Index, $Value);
  }

  # Check and set...
  if ($Index < 0) {
    croak "Error: ${ClassName}->SetValue: Index value must be a positive number...";
  }
  if ($Index >= $This->GetNumOfValues()) {
    croak "Error: ${ClassName}->SetValue: Index vaue must be less than number of values...";
  }

  return $This->_SetValue($Index, $Value);
}

# Set a fingerprint vector value...
#
sub _SetValue {
  my($This, $Index, $Value) = @_;

  $This->{Values}[$Index] = $Value;

  return $This;
}

# Get a specific value from fingerprint vector with indicies starting from 0...
#
sub GetValue {
  my($This, $Index) = @_;

  if ($Index < 0) {
    croak "Error: ${ClassName}->GetValue: Index value must be a positive number...";
  }
  if ($Index >= $This->GetNumOfValues()) {
    croak "Error: ${ClassName}->GetValue: Index value must be less than number of values...";
  }
  return $This->_GetValue($Index);
}

# Get a fingerprint vector value...
sub _GetValue {
  my($This, $Index) = @_;

  return $This->{Values}[$Index];
}

# Return vector values as an array or reference to an array...
#
sub GetValues {
  my($This) = @_;

  return wantarray ? @{$This->{Values}} : \@{$This->{Values}};
}

# Set a specific value ID in fingerprint vector with indicies starting from 0..
#
sub SetValueID {
  my($This, $Index, $Value, $SkipCheck) = @_;

  # Just set it...
  if ($SkipCheck) {
    return $This->_SetValueID($Index, $Value);
  }

  # Check and set...
  if ($Index < 0) {
    croak "Error: ${ClassName}->SetValueID: Index value must be a positive number...";
  }
  if ($Index >= $This->GetNumOfValueIDs()) {
    croak "Error: ${ClassName}->SetValueID: Index vaue must be less than number of value IDs...";
  }

  return $This->_SetValueID($Index, $Value);
}

# Set a fingerprint vector value ID...
#
sub _SetValueID {
  my($This, $Index, $Value) = @_;

  $This->{ValueIDs}[$Index] = $Value;

  return $This;
}

# Get a specific value ID from fingerprint vector with indicies starting from 0...
#
sub GetValueID {
  my($This, $Index) = @_;

  if ($Index < 0) {
    croak "Error: ${ClassName}->GetValueID: Index value must be a positive number...";
  }
  if ($Index >= $This->GetNumOfValueIDs()) {
    croak "Error: ${ClassName}->GetValueID: Index value must be less than number of value IDs...";
  }
  return $This->_GetValueID($Index);
}

# Get a fingerprint vector value ID...
#
sub _GetValueID {
  my($This, $Index) = @_;

  return $This->{ValueIDs}[$Index];
}

# Return vector value IDs as an array or reference to an array...
#
sub GetValueIDs {
  my($This) = @_;

  return wantarray ? @{$This->{ValueIDs}} : \@{$This->{ValueIDs}};
}

# Get fingerprints vector string containing values and/or IDs string in a specifed format...
#
sub GetFingerprintsVectorString {
  my($This, $Format) = @_;

  FORMAT : {
    if ($Format =~ /^(IDsAndValuesString|IDsAndValues)$/i) { return $This->GetIDsAndValuesString(); last FORMAT; }
    if ($Format =~ /^(IDsAndValuesPairsString|IDsAndValuesPairs)$/i) { return $This->GetIDsAndValuesPairsString(); last FORMAT; }
    if ($Format =~ /^(ValuesAndIDsString|ValuesAndIDs)$/i) { return $This->GetValuesAndIDsString(); last FORMAT; }
    if ($Format =~ /^(ValuesAndIDsPairsString|ValuesAndIDsPairs)$/i) { return $This->GetValuesAndIDsPairsString(); last FORMAT;}
    if ($Format =~ /^(ValueIDsString|ValueIDs)$/i) { return $This->GetValueIDsString(); last FORMAT; }
    if ($Format =~ /^(ValuesString|Values)$/i) { return $This->GetValuesString(); last FORMAT; }
    croak "Error: ${ClassName}->GetFingerprintsVectorString: Specified vector string format, $Format, is not supported. Value values: IDsAndValuesString, IDsAndValues, IDsAndValuesPairsString, IDsAndValuesPairs, ValuesAndIDsString, ValuesAndIDs, ValuesAndIDsPairsString, ValuesAndIDsPairs, ValueIDsString, ValueIDs, ValuesString, Values...";
  }
  return '';
}
# Get vector value IDs and values string as space delimited ASCII string separated
# by semicolon...
#
sub GetIDsAndValuesString {
  my($This) = @_;

  if (@{$This->{ValueIDs}} && @{$This->{Values}}) {
    # Both IDs and values are available...
    return join(' ', @{$This->{ValueIDs}}) . ";" . join(' ', @{$This->{Values}});
  }
  elsif (@{$This->{Values}}) {
    # Only values are available...
    return "None;" . join(' ', @{$This->{Values}});
  }
  else {
    # Values are not available...
    return "None;None";
  }
}

# Get vector value IDs and value pairs string as space delimited ASCII string...
#
sub GetIDsAndValuesPairsString {
  my($This) = @_;
  my($Index, $ValueIDsPresent, @IDsAndValuesPairs);

  if (!@{$This->{Values}}) {
    # Values are unavailable...
    return "None None";
  }

  $ValueIDsPresent = @{$This->{ValueIDs}} ? 1 : 0;

  @IDsAndValuesPairs = ();
  for $Index (0 .. $#{$This->{Values}}) {
    if ($ValueIDsPresent) {
      push @IDsAndValuesPairs, ($This->{ValueIDs}->[$Index], $This->{Values}->[$Index]);
    }
    else {
      push @IDsAndValuesPairs, ('None', $This->{Values}->[$Index]);
    }
  }
  return join(' ', @IDsAndValuesPairs);
}

# Get vector value and value IDs string as space delimited ASCII string separated
# by semicolon...
#
sub GetValuesAndIDsString {
  my($This) = @_;

  if (@{$This->{ValueIDs}} && @{$This->{Values}}) {
    # Both IDs and values are available...
    return join(' ', @{$This->{Values}}) . ";" . join(' ', @{$This->{ValueIDs}});
  }
  elsif (@{$This->{Values}}) {
    # Only values are available...
    return join(' ', @{$This->{Values}}) . ";None";
  }
  else {
    # Values are not available...
    return "None;None";
  }
}

# Get vector value and value ID pairs string as space delimited ASCII string...
#
sub GetValuesAndIDsPairsString {
  my($This) = @_;
  my($Index, $ValueIDsPresent, @ValuesAndIDsPairs);

  if (!@{$This->{Values}}) {
    # Values are unavailable...
    return "None None";
  }

  $ValueIDsPresent = @{$This->{ValueIDs}} ? 1 : 0;

  @ValuesAndIDsPairs = ();
  for $Index (0 .. $#{$This->{Values}}) {
    if ($ValueIDsPresent) {
      push @ValuesAndIDsPairs, ($This->{Values}->[$Index], $This->{ValueIDs}->[$Index]);
    }
    else {
      push @ValuesAndIDsPairs, ($This->{Values}->[$Index], 'None');
    }
  }
  return join(' ', @ValuesAndIDsPairs);
}

# Get vector value IDs string as space delimited ASCII string...
#
sub GetValueIDsString {
  my($This) = @_;

  return @{$This->{ValueIDs}} ? join(' ', @{$This->{ValueIDs}}) : 'None';
}

# Get vector value string as space delimited ASCII string...
#
sub GetValuesString {
  my($This) = @_;

  return @{$This->{Values}} ? join(' ', @{$This->{Values}}) : 'None';
}

# Get number of values...
sub GetNumOfValues {
  my($This) = @_;

  return scalar @{$This->{Values}};
}

# Get number of non-zero values...
sub GetNumOfNonZeroValues {
  my($This) = @_;
  my($Count, $Index, $Size);

  $Count = 0;
  $Size = $This->GetNumOfValues();

  for $Index (0 .. ($Size -1)) {
    if ($This->{Values}[$Index] != 0) {
      $Count++;
    }
  }
  return $Count;
}

# Get number of value IDs...
sub GetNumOfValueIDs {
  my($This) = @_;

  return scalar @{$This->{ValueIDs}};
}

# FinegerprintsVectors class provides methods to calculate similarity between vectors
# containing three different types of values:
#
# Type I: OrderedNumericalValues
#
#   . Size of two vectors are same
#   . Vectors contain real values in a specific order. For example: MACCS keys count, Topological
#     pharnacophore atom pairs and so on.
#   . Option to calculate similarity value using continious values or binary values
#
# Type II: UnorderedNumericalValues
#
#   . Size of two vectors might not be same
#   . Vectors contain unordered real value identified by value IDs. For example: Toplogical atom pairs,
#     Topological atom torsions and so on
#   . Option to calculate similarity value using continous values or binary values
#
# Type III: AlphaNumericalValues
#
#   . Size of two vectors might not be same
#   . Vectors contain unordered alphanumerical values. For example: Extended connectivity fingerprints,
#     atom neighbothood fingerpritns.
#   . The vector values are treated as keys or bit indices and similarity value is calculated accordingly.
#
# Before performing similarity or distance calculations between vectors containing UnorderedNumericalValues
# or AlphaNumericalValues, the vectors are tranformed into vectors containing unique OrderedNumericalValues
# using value IDs for UnorderedNumericalValues and values itself for AlphaNumericalValues.
#
# Three forms similarity or distance calculation between two vectors: AlgebraicForm, BinaryForm or
# SetTheoreticForm.
#
# The value of an extra paramter, CalculationMode, passed to each similarity or distance function
# controls the calculation. Supported values for CalculationMode: AlgebraicForm, BinaryForm and
# SetTheoreticForm. Default: AlgebraicForm.
#
# For BinaryForm CalculationMode, the ordered list of processed final vector values containing the value or
# count of each unique value type is simply converted into a binary vector containing 1s and 0s
# corresponding to presence or absence of values before calculating similarity or distance between
# two vectors.
#
# For two fingerprint vectors A and B of same size containing OrderedNumericalValues, let:
#
#  N = Number values in A or B
#
#  Xa = Values of vector A
#  Xb = Values of vector B
#
#  Xai = Value of ith element in A
#  Xbi = Value of ith element in B
#
#  SUM = Sum of i over N values
#
# For SetTheoreticForm of calculation between two vectors, let:
#
#  SetIntersectionXaXb = SUM ( MIN ( Xai, Xbi ) )
#  SetDifferenceXaXb = SUM ( Xai ) + SUM ( Xbi ) - SUM ( MIN ( Xai, Xbi ) )
#
# For BinaryForm of calculation between two vectors, let:
#
#  Na = Number of bits set to "1" in A = SUM ( Xai )
#  Nb = Number of bits set to "1" in B = SUM ( Xbi )
#  Nc = Number of bits set to "1" in both A and B = SUM ( Xai * Xbi )
#  Nd = Number of bits set to "0" in both A and B = SUM ( 1 - Xai - Xbi + Xai * Xbi)
#
#  N = Number of bits set to "1" or "0" in A or B = Size of A or B = Na + Nb - Nc + Nd
#
# Additionally, for BinaryForm various values also correspond to:
#
#  Na = | Xa |
#  Nb = | Xb |
#  Nc = | SetIntersectionXaXb |
#  Nd = N - | SetDifferenceXaXb |
#
#  | SetDifferenceXaXb | = N - Nd = Na + Nb - Nc + Nd - Nd = Na + Nb - Nc
#                        =  | Xa | + | Xb | - | SetIntersectionXaXb |
#
# Various distance coefficients and similarity coefficients [ Ref 40, Ref 62, Ref 64 ] for a pair vectors A and B
# in AlgebraicForm and BinaryForm are defined as follows:
#
# . CityBlockDistanceCoefficient: ( same as HammingDistanceCoefficient and ManhattanDistanceCoefficient)
#
#     . AlgebraicForm: SUM ( ABS ( Xai - Xbi ) )
#
#     . BinaryForm: ( Na - Nc ) + ( Nb - Nc ) = Na + Nb - 2 * Nc
#
#     . SetTheoreticForm: | SetDifferenceXaXb | - | SetIntersectionXaXb |
#                        = SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) )
#
# . CosineSimilarityCoefficient:  ( same as OchiaiSimilarityCoefficient)
#
#     . AlgebraicForm: SUM ( Xai * Xbi ) / SQRT ( SUM ( Xai ** 2) * SUM ( Xbi ** 2) )
#
#     . BinaryForm: Nc / SQRT ( Na * Nb)
#
#     . SetTheoreticForm: | SetIntersectionXaXb | / SQRT ( |Xa| * |Xb| )
#                        = SUM ( MIN ( Xai, Xbi ) ) / SQRT ( SUM ( Xai ) * SUM ( Xbi ) )
#
# . CzekanowskiSimilarityCoefficient: ( same as DiceSimilarityCoefficient and SorensonSimilarityCoefficient)
#
#     . AlgebraicForm: ( 2 * ( SUM ( Xai * Xbi ) )  ) / ( SUM ( Xai ** 2) + SUM ( Xbi **2 ) )
#
#     . BinaryForm: 2 * Nc / ( Na + Nb )
#
#     . SetTheoreticForm: 2 * | SetIntersectionXaXb | / ( |Xa| + |Xb| )
#                        = 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) / ( SUM ( Xai ) + SUM ( Xbi ) )
#
# . DiceSimilarityCoefficient: ( same as CzekanowskiSimilarityCoefficient and SorensonSimilarityCoefficient)
#
#     . AlgebraicForm: ( 2 * ( SUM ( Xai * Xbi ) )  ) / ( SUM ( Xai ** 2) + SUM ( Xbi **2 ) )
#
#     . BinaryForm: 2 * Nc / ( Na + Nb )
#
#     . SetTheoreticForm: 2 * | SetIntersectionXaXb | / ( |Xa| + |Xb| )
#                        = 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) / ( SUM ( Xai ) + SUM ( Xbi ) )
#
# . EuclideanDistanceCoefficient:
#
#     . AlgebraicForm: SQRT ( SUM ( ( ( Xai - Xbi ) ** 2 ) ) )
#
#     . BinaryForm: SQRT ( ( Na - Nc ) + ( Nb - Nc ) ) = SQRT ( Na + Nb - 2 * Nc )
#
#     . SetTheoreticForm: SQRT ( | SetDifferenceXaXb | - | SetIntersectionXaXb | )
#                        = SQRT (  SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) )
#
# . HammingDistanceCoefficient:  ( same as CityBlockDistanceCoefficient and ManhattanDistanceCoefficient)
#
#     . AlgebraicForm: SUM ( ABS ( Xai - Xbi ) )
#
#     . BinaryForm: ( Na - Nc ) + ( Nb - Nc ) = Na + Nb - 2 * Nc
#
#     . SetTheoreticForm: | SetDifferenceXaXb | - | SetIntersectionXaXb |
#                        = SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) )
#
# . JaccardSimilarityCoefficient: ( same as TanimotoSimilarityCoefficient)
#
#     . AlgebraicForm:  SUM ( Xai * Xbi ) / ( SUM ( Xai ** 2 ) + SUM ( Xbi ** 2 ) - SUM ( Xai * Xbi ) )
#
#     . BinaryForm:  Nc / ( ( Na - Nc ) + ( Nb - Nc ) + Nc ) = Nc / ( Na + Nb - Nc )
#
#     . SetTheoreticForm: | SetIntersectionXaXb | / | SetDifferenceXaXb |
#                        = SUM ( MIN ( Xai, Xbi ) ) / (  SUM ( Xai ) + SUM ( Xbi ) - SUM ( MIN ( Xai, Xbi ) ) )
#
# . ManhattanDistanceCoefficient:  ( same as CityBlockDistanceCoefficient and HammingDistanceCoefficient)
#
#     . AlgebraicForm: SUM ( ABS ( Xai - Xbi ) )
#
#     . BinaryForm: ( Na - Nc ) + ( Nb - Nc ) = Na + Nb - 2 * Nc
#
#     . SetTheoreticForm: | SetDifferenceXaXb | - | SetIntersectionXaXb |
#                        = SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) )
#
# . OchiaiSimilarityCoefficient:  ( same as CosineSimilarityCoefficient)
#
#     . AlgebraicForm: SUM ( Xai * Xbi ) / SQRT ( SUM ( Xai ** 2) * SUM ( Xbi ** 2) )
#
#     . BinaryForm: Nc / SQRT ( Na * Nb)
#
#     . SetTheoreticForm: | SetIntersectionXaXb | / SQRT ( |Xa| * |Xb| )
#                        = SUM ( MIN ( Xai, Xbi ) ) / SQRT ( SUM ( Xai ) * SUM ( Xbi ) )
#
# . SorensonSimilarityCoefficient: ( same as CzekanowskiSimilarityCoefficient and DiceSimilarityCoefficient)
#
#     . AlgebraicForm: ( 2 * ( SUM ( Xai * Xbi ) )  ) / ( SUM ( Xai ** 2) + SUM ( Xbi **2 ) )
#
#     . BinaryForm: 2 * Nc / ( Na + Nb )
#
#     . SetTheoreticForm: 2 * | SetIntersectionXaXb | / ( |Xa| + |Xb| )
#                        = 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) / ( SUM ( Xai ) + SUM ( Xbi ) )
#
# . SoergelDistanceCoefficient:
#
#     . AlgebraicForm:  SUM ( ABS ( Xai - Xbi ) ) / SUM ( MAX ( Xai, Xbi ) )
#
#     . BinaryForm: 1 - Nc / ( Na + Nb - Nc ) = ( Na + Nb - 2 * Nc ) / ( Na + Nb - Nc )
#
#     . SetTheoreticForm: ( | SetDifferenceXaXb | - | SetIntersectionXaXb | ) / | SetDifferenceXaXb |
#                        = ( SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) ) / ( SUM ( Xai ) + SUM ( Xbi ) - SUM ( MIN ( Xai, Xbi ) ) )
#
# . TanimotoSimilarityCoefficient:  ( same as JaccardSimilarityCoefficient)
#
#     . AlgebraicForm:  SUM ( Xai * Xbi ) / ( SUM ( Xai ** 2 ) + SUM ( Xbi ** 2 ) - SUM ( Xai * Xbi ) )
#
#     . BinaryForm:  Nc / ( ( Na - Nc ) + ( Nb - Nc ) + Nc ) = Nc / ( Na + Nb - Nc )
#
#     . SetTheoreticForm: | SetIntersectionXaXb | / | SetDifferenceXaXb |
#                        = SUM ( MIN ( Xai, Xbi ) ) / (  SUM ( Xai ) + SUM ( Xbi ) - SUM ( MIN ( Xai, Xbi ) ) )
#
#

# Calculate Hamming distance coefficient between two fingerprint vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub HammingDistanceCoefficient ($$;$$) {
  my($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck) = @_;

  return CityBlockDistanceCoefficient($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck);
}

# Calculate Hamming distance coefficient between two fingerprint vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub ManhattanDistanceCoefficient ($$;$$) {
  my($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck) = @_;

  return CityBlockDistanceCoefficient($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck);
}

# Calculate CityBlock distance coefficient between two fingerprint vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub CityBlockDistanceCoefficient ($$;$$) {
  my($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck) = @_;

  $CalculationMode = defined $CalculationMode ? $CalculationMode : 'AlgebraicForm';
  $SkipValuesCheck = defined $SkipValuesCheck ? $SkipValuesCheck : 0;

  # Validate and process fingerprints vectors for similarity calculations...
  #
  _ValidateAndProcessFingerprintsVectorsForSimilarityCalculation("CityBlockDistanceCoefficient: Calculation failed", $FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck);

  # Perform the calculation...
  if ($CalculationMode =~ /^AlgebraicForm$/i) {
    return _CityBlockDistanceCoefficientUsingAlgebraicForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($CalculationMode =~ /^BinaryForm$/i) {
    return _CityBlockDistanceCoefficientUsingBinaryForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($CalculationMode =~ /^SetTheoreticForm$/i) {
    return _CityBlockDistanceCoefficientUsingSetTheoreticForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  else {
    return undef;
  }
}

# Calculate CityBlock distance coefficient using algebraic form...
#
sub _CityBlockDistanceCoefficientUsingAlgebraicForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($SumAbsSubtractionXaiXbi);

  $SumAbsSubtractionXaiXbi = _GetSumOfAbsoluteValueOfSubtractionOfFingerprintsOrderedValues($FingerprintsVectorA, $FingerprintsVectorB);

  return $SumAbsSubtractionXaiXbi;
}

# Calculate CityBlock distance coefficient using binary form...
#
sub _CityBlockDistanceCoefficientUsingBinaryForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($Na, $Nb, $Nc);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsVectorA, $FingerprintsVectorB);

  return  ($Na + $Nb - 2 * $Nc);
}

# Calculate  CityBlock distance coefficient using set theoretic form...
#
sub _CityBlockDistanceCoefficientUsingSetTheoreticForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($SumMinXaiXbi, $SumXai, $SumXbi);

  $SumXai = _GetSumOfFingerprintsOrderedValues($FingerprintsVectorA);
  $SumXbi = _GetSumOfFingerprintsOrderedValues($FingerprintsVectorB);
  $SumMinXaiXbi = _GetSumOfMinimumOfFingerprintsOrderdedValues($FingerprintsVectorA, $FingerprintsVectorB);

  return  ($SumXai + $SumXbi - 2 * $SumMinXaiXbi);
}

# Calculate Ochiai similarity cofficient between two fingerprint vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub OchiaiSimilarityCoefficient ($$;$$) {
  my($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck) = @_;

  return CosineSimilarityCoefficient($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck);
}

# Calculate Cosine similarity cofficient between two fingerprint vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub CosineSimilarityCoefficient ($$;$$) {
  my($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck) = @_;

  $CalculationMode = defined $CalculationMode ? $CalculationMode : 'AlgebraicForm';
  $SkipValuesCheck = defined $SkipValuesCheck ? $SkipValuesCheck : 0;

  # Validate and process fingerprints vectors for similarity calculations...
  #
  _ValidateAndProcessFingerprintsVectorsForSimilarityCalculation("CosineSimilarityCoefficient: Calculation failed", $FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck);

  # Perform the calculation...
  if ($CalculationMode =~ /^AlgebraicForm$/i) {
    return _CosineSimilarityCoefficientUsingAlgebraicForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($CalculationMode =~ /^BinaryForm$/i) {
    return _CosineSimilarityCoefficientUsingBinaryForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($CalculationMode =~ /^SetTheoreticForm$/i) {
    return _CosineSimilarityCoefficientUsingSetTheoreticForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  else {
    return undef;
  }
}

# Calculate Cosine similarity coefficient using algebraic form...
#
sub _CosineSimilarityCoefficientUsingAlgebraicForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($SumProductXaiXbi, $SumXai2, $SumXbi2, $Numerator, $Denominator);

  $SumXai2 = _GetSumOfSquaresOfFingerprintsOrderedValues($FingerprintsVectorA);
  $SumXbi2 = _GetSumOfSquaresOfFingerprintsOrderedValues($FingerprintsVectorB);
  $SumProductXaiXbi = _GetSumOfProductOfFingerprintsOrderedValues($FingerprintsVectorA, $FingerprintsVectorB);

  $Numerator = $SumProductXaiXbi;
  $Denominator = sqrt($SumXai2 * $SumXbi2);

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# CalculateCosine similarity coefficient using binary form...
#
sub _CosineSimilarityCoefficientUsingBinaryForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($Na, $Nb, $Nc, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsVectorA, $FingerprintsVectorB);

  $Numerator = $Nc;
  $Denominator = sqrt($Na * $Nb);

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Cosine similarity coefficient using set theoretic form...
#
sub _CosineSimilarityCoefficientUsingSetTheoreticForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($SumMinXaiXbi, $SumXai, $SumXbi, $Numerator, $Denominator);

  $SumXai = _GetSumOfFingerprintsOrderedValues($FingerprintsVectorA);
  $SumXbi = _GetSumOfFingerprintsOrderedValues($FingerprintsVectorB);
  $SumMinXaiXbi = _GetSumOfMinimumOfFingerprintsOrderdedValues($FingerprintsVectorA, $FingerprintsVectorB);

  $Numerator = $SumMinXaiXbi;
  $Denominator = sqrt($SumXai * $SumXbi);

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Czekanowski similarity cofficient between two fingerprint vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub CzekanowskiSimilarityCoefficient ($$;$$) {
  my($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck) = @_;

  return DiceSimilarityCoefficient($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck);
}

# Calculate Sorenson similarity cofficient between two fingerprint vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub SorensonSimilarityCoefficient ($$;$$) {
  my($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck) = @_;

  return DiceSimilarityCoefficient($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck);
}

# Calculate Dice similarity cofficient between two fingerprint vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub DiceSimilarityCoefficient ($$;$$) {
  my($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck) = @_;

  $CalculationMode = defined $CalculationMode ? $CalculationMode : 'AlgebraicForm';
  $SkipValuesCheck = defined $SkipValuesCheck ? $SkipValuesCheck : 0;

  # Validate and process fingerprints vectors for similarity calculations...
  #
  _ValidateAndProcessFingerprintsVectorsForSimilarityCalculation("DiceSimilarityCoefficient: Calculation failed", $FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck);

  # Perform the calculation...
  if ($CalculationMode =~ /^AlgebraicForm$/i) {
    return _DiceSimilarityCoefficientUsingAlgebraicForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($CalculationMode =~ /^BinaryForm$/i) {
    return _DiceSimilarityCoefficientUsingBinaryForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($CalculationMode =~ /^SetTheoreticForm$/i) {
    return _DiceSimilarityCoefficientUsingSetTheoreticForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  else {
    return undef;
  }
}

# Calculate Dice similarity coefficient using algebraic form...
#
sub _DiceSimilarityCoefficientUsingAlgebraicForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($SumProductXaiXbi, $SumXai2, $SumXbi2, $Numerator, $Denominator);

  $SumXai2 = _GetSumOfSquaresOfFingerprintsOrderedValues($FingerprintsVectorA);
  $SumXbi2 = _GetSumOfSquaresOfFingerprintsOrderedValues($FingerprintsVectorB);
  $SumProductXaiXbi = _GetSumOfProductOfFingerprintsOrderedValues($FingerprintsVectorA, $FingerprintsVectorB);

  $Numerator = 2 * $SumProductXaiXbi;
  $Denominator = $SumXai2 + $SumXbi2;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Dice similarity coefficient using binary form...
#
sub _DiceSimilarityCoefficientUsingBinaryForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($Na, $Nb, $Nc, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsVectorA, $FingerprintsVectorB);

  $Numerator = 2 * $Nc;
  $Denominator = $Na + $Nb;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Dice similarity coefficient using set theoretic form...
#
sub _DiceSimilarityCoefficientUsingSetTheoreticForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($SumMinXaiXbi, $SumXai, $SumXbi, $Numerator, $Denominator);

  $SumXai = _GetSumOfFingerprintsOrderedValues($FingerprintsVectorA);
  $SumXbi = _GetSumOfFingerprintsOrderedValues($FingerprintsVectorB);
  $SumMinXaiXbi = _GetSumOfMinimumOfFingerprintsOrderdedValues($FingerprintsVectorA, $FingerprintsVectorB);

  $Numerator = 2 * $SumMinXaiXbi;
  $Denominator = $SumXai + $SumXbi;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}


# Calculate Euclidean distance coefficient between two fingerprint vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub EuclideanDistanceCoefficient ($$;$$) {
  my($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck) = @_;

  $CalculationMode = defined $CalculationMode ? $CalculationMode : 'AlgebraicForm';
  $SkipValuesCheck = defined $SkipValuesCheck ? $SkipValuesCheck : 0;

  # Validate and process fingerprints vectors for similarity calculations...
  #
  _ValidateAndProcessFingerprintsVectorsForSimilarityCalculation("EuclideanDistanceCoefficient: Calculation failed", $FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck);

  # Perform the calculation...
  if ($CalculationMode =~ /^AlgebraicForm$/i) {
    return _EuclideanDistanceCoefficientUsingAlgebraicForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($CalculationMode =~ /^BinaryForm$/i) {
    return _EuclideanDistanceCoefficientUsingBinaryForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($CalculationMode =~ /^SetTheoreticForm$/i) {
    return _EuclideanDistanceCoefficientUsingSetTheoreticForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  else {
    return undef;
  }
}

# Calculate Euclidean distance coefficient using algebraic form...
#
sub _EuclideanDistanceCoefficientUsingAlgebraicForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($SumSquaresSubtractionXaiXbi);

  $SumSquaresSubtractionXaiXbi = _GetSumOfSquaresOfSubtractionOfFingerprintsOrderedValues($FingerprintsVectorA, $FingerprintsVectorB);

  return sqrt($SumSquaresSubtractionXaiXbi);
}

# Calculate Euclidean distance coefficient using binary form...
#
sub _EuclideanDistanceCoefficientUsingBinaryForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($Na, $Nb, $Nc);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsVectorA, $FingerprintsVectorB);

  return  (sqrt($Na + $Nb - 2 * $Nc));
}

# Calculate Euclidean distance coefficient using set theoretic form...
#
sub _EuclideanDistanceCoefficientUsingSetTheoreticForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($SumMinXaiXbi, $SumXai, $SumXbi);

  $SumXai = _GetSumOfFingerprintsOrderedValues($FingerprintsVectorA);
  $SumXbi = _GetSumOfFingerprintsOrderedValues($FingerprintsVectorB);
  $SumMinXaiXbi = _GetSumOfMinimumOfFingerprintsOrderdedValues($FingerprintsVectorA, $FingerprintsVectorB);

  return  (sqrt($SumXai + $SumXbi - 2 * $SumMinXaiXbi));
}

# Calculate Jaccard similarity cofficient between two fingerprint vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub JaccardSimilarityCoefficient ($$;$$) {
  my($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck) = @_;

  return TanimotoSimilarityCoefficient($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck);
}

# Calculate Tanimoto similarity cofficient between two fingerprint vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub TanimotoSimilarityCoefficient ($$;$$) {
  my($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck) = @_;

  $CalculationMode = defined $CalculationMode ? $CalculationMode : 'AlgebraicForm';
  $SkipValuesCheck = defined $SkipValuesCheck ? $SkipValuesCheck : 0;

  # Validate and process fingerprints vectors for similarity calculations...
  #
  _ValidateAndProcessFingerprintsVectorsForSimilarityCalculation("TanimotoSimilarityCoefficient: Calculation failed", $FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck);

  # Perform the calculation...
  if ($CalculationMode =~ /^AlgebraicForm$/i) {
    return _TanimotoSimilarityCoefficientUsingAlgebraicForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($CalculationMode =~ /^BinaryForm$/i) {
    return _TanimotoSimilarityCoefficientUsingBinaryForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($CalculationMode =~ /^SetTheoreticForm$/i) {
    return _TanimotoSimilarityCoefficientUsingSetTheoreticForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  else {
    return undef;
  }
}

# Calculate Tanimoto similarity coefficient using algebraic form...
#
sub _TanimotoSimilarityCoefficientUsingAlgebraicForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($SumProductXaiXbi, $SumXai2, $SumXbi2, $Numerator, $Denominator);

  $SumXai2 = _GetSumOfSquaresOfFingerprintsOrderedValues($FingerprintsVectorA);
  $SumXbi2 = _GetSumOfSquaresOfFingerprintsOrderedValues($FingerprintsVectorB);
  $SumProductXaiXbi = _GetSumOfProductOfFingerprintsOrderedValues($FingerprintsVectorA, $FingerprintsVectorB);

  $Numerator = $SumProductXaiXbi;
  $Denominator = $SumXai2 + $SumXbi2 - $SumProductXaiXbi;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Tanimoto similarity coefficient using binary form...
#
sub _TanimotoSimilarityCoefficientUsingBinaryForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($Na, $Nb, $Nc, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsVectorA, $FingerprintsVectorB);

  $Numerator = $Nc;
  $Denominator = $Na + $Nb - $Nc;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Tanimoto similarity coefficient using set theoretic form...
#
sub _TanimotoSimilarityCoefficientUsingSetTheoreticForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($SumMinXaiXbi, $SumXai, $SumXbi, $Numerator, $Denominator);

  $SumXai = _GetSumOfFingerprintsOrderedValues($FingerprintsVectorA);
  $SumXbi = _GetSumOfFingerprintsOrderedValues($FingerprintsVectorB);
  $SumMinXaiXbi = _GetSumOfMinimumOfFingerprintsOrderdedValues($FingerprintsVectorA, $FingerprintsVectorB);

  $Numerator = $SumMinXaiXbi;
  $Denominator = $SumXai + $SumXbi - $SumMinXaiXbi;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}


# Calculate Soergel distance coefficient between two fingerprint vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub SoergelDistanceCoefficient ($$;$$) {
  my($FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck) = @_;

  $CalculationMode = defined $CalculationMode ? $CalculationMode : 'AlgebraicForm';
  $SkipValuesCheck = defined $SkipValuesCheck ? $SkipValuesCheck : 0;

  # Validate and process fingerprints vectors for similarity calculations...
  #
  _ValidateAndProcessFingerprintsVectorsForSimilarityCalculation("SoergelDistanceCoefficient: Calculation failed", $FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck);

  # Perform the calculation...
  if ($CalculationMode =~ /^AlgebraicForm$/i) {
    return _SoergelDistanceCoefficientUsingAlgebraicForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($CalculationMode =~ /^BinaryForm$/i) {
    return _SoergelDistanceCoefficientUsingBinaryForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($CalculationMode =~ /^SetTheoreticForm$/i) {
    return _SoergelDistanceCoefficientUsingSetTheoreticForm($FingerprintsVectorA, $FingerprintsVectorB);
  }
  else {
    return undef;
  }
}

# Calculate Soergel distance coefficientusing algebraic form...
#
sub _SoergelDistanceCoefficientUsingAlgebraicForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($SumAbsSubtractionXaiXbi, $SumMaxXaiXbi, $Numerator, $Denominator);

  $SumAbsSubtractionXaiXbi = _GetSumOfAbsoluteValueOfSubtractionOfFingerprintsOrderedValues($FingerprintsVectorA, $FingerprintsVectorB);
  $SumMaxXaiXbi = _GetSumOfMaximumOfFingerprintsOrderdedValues($FingerprintsVectorA, $FingerprintsVectorB);

  $Numerator = $SumAbsSubtractionXaiXbi;
  $Denominator = $SumMaxXaiXbi;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Soergel distance coefficient using binary form...
#
sub _SoergelDistanceCoefficientUsingBinaryForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($Na, $Nb, $Nc, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsVectorA, $FingerprintsVectorB);

  $Numerator = $Na + $Nb - 2 * $Nc;
  $Denominator = $Na + $Nb - $Nc;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate SoergelDistanceCoefficient using set theoretic form...
#
sub _SoergelDistanceCoefficientUsingSetTheoreticForm {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($SumMinXaiXbi, $SumXai, $SumXbi, $Numerator, $Denominator);

  $SumXai = _GetSumOfFingerprintsOrderedValues($FingerprintsVectorA);
  $SumXbi = _GetSumOfFingerprintsOrderedValues($FingerprintsVectorB);
  $SumMinXaiXbi = _GetSumOfMinimumOfFingerprintsOrderdedValues($FingerprintsVectorA, $FingerprintsVectorB);

  $Numerator = $SumXai + $SumXbi - 2 * $SumMinXaiXbi;
  $Denominator = $SumXai + $SumXbi - $SumMinXaiXbi;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Validate and process fingerprints vectors for similarity calculations...
#
sub _ValidateAndProcessFingerprintsVectorsForSimilarityCalculation {
  my($ErrorMsg, $FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode, $SkipValuesCheck) = @_;

  $CalculationMode = defined $CalculationMode ? $CalculationMode : 'AlgebraicForm';
  $SkipValuesCheck = defined $SkipValuesCheck ? $SkipValuesCheck : 0;

  if (!$SkipValuesCheck) {
    _ValidateFingerprintsVectorsForSimilarityCalculation($ErrorMsg, $FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode);
  }
  _ProcessFingerprintsVectorsForSimilarityCalculation($ErrorMsg, $FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode);
}

# Make sure fingerprint vectors are good for performing similarity/distance calculation...
#
sub _ValidateFingerprintsVectorsForSimilarityCalculation {
  my($ErrorMsg, $FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode) = @_;

  # Make sure both are fingerprint vectors..
  if (!(IsFingerprintsVector($FingerprintsVectorA) && IsFingerprintsVector($FingerprintsVectorB))) {
    croak "Error: ${ClassName}->${ErrorMsg}: Both objects must be fingerprint vectors...";
  }

  # Check types...
  if ($FingerprintsVectorA->{Type} ne $FingerprintsVectorB->{Type}) {
    croak "Error: ${ClassName}->${ErrorMsg}: Type of first fingerprint vector, $FingerprintsVectorA->{Type}, must be same as type of second fingerprint vector, $FingerprintsVectorB->{Type}...";
  }

  # Check calculation mode...
  if ($CalculationMode !~ /^(AlgebraicForm|BinaryForm|SetTheoreticForm)$/i) {
    croak "Error: ${ClassName}->${ErrorMsg}: Specified similarity calculation mode, $CalculationMode, is not valid. Supported values: AlgebraicForm, BinaryForm, and SetTheoreticForm...";
  }

  # Check values and value IDs...
  my($Na, $Nb, $NIDa, $NIDb);
  $Na = $FingerprintsVectorA->GetNumOfValues(); $Nb = $FingerprintsVectorB->GetNumOfValues();
  $NIDa = $FingerprintsVectorA->GetNumOfValueIDs(); $NIDb = $FingerprintsVectorB->GetNumOfValueIDs();

  if ($Na == 0) {
    croak "Error: ${ClassName}->${ErrorMsg}: Number of values in first fingerprint vector, $Na, must be > 0 for fingerprint vector type $FingerprintsVectorA->{Type} ...";
  }
  if ($Nb == 0) {
    croak "Error: ${ClassName}->${ErrorMsg}: Number of values in second fingerprint vector, $Nb, must be > 0 for fingerprint vector type $FingerprintsVectorB->{Type} ...";
  }

  if ($FingerprintsVectorA->{Type} =~ /^OrderedNumericalValues$/i) {
    if ($Na != $Nb) {
      croak "Error: ${ClassName}->${ErrorMsg}: Number of values in first fingerprint vector, $Na, must be equal to number of values, $Nb, in second fingerprint vector for fingerprint vector types $FingerprintsVectorA->{Type} ...";
    }
  }
  elsif ($FingerprintsVectorA->{Type} =~ /^NumericalValues$/i) {
    if ($NIDa == 0) {
      croak "Error: ${ClassName}->${ErrorMsg}: Number of value IDs in first fingerprint vector, $NIDa, must be > 0 for fingerprint vector type $FingerprintsVectorA->{Type} ...";
    }
    if ($NIDb == 0) {
      croak "Error: ${ClassName}->${ErrorMsg}: Number of value IDs in first fingerprint vector, $NIDb, must be > 0 for fingerprint vector type $FingerprintsVectorB->{Type} ...";
    }

    if ($NIDa != $Na) {
      croak "Error: ${ClassName}->${ErrorMsg}: Number of value IDs in first fingerprint vector, $NIDa, must be equal to its number of values, $Na, for fingerprint vector type $FingerprintsVectorA->{Type} ...";
    }
    if ($NIDb != $Nb) {
      croak "Error: ${ClassName}->${ErrorMsg}: Number of value IDs in second fingerprint vector, $NIDb, must be equal to its number of values, $Nb, for fingerprint vector type $FingerprintsVectorA->{Type} ...";
    }
  }
  elsif ($FingerprintsVectorA->{Type} =~ /^AlphaNumericalValues$/i) {
    if ($NIDa || $NIDb) {
      croak "Error: ${ClassName}->${ErrorMsg}: ValueIDs cann't be specified for fingerprint vector types $FingerprintsVectorA->{Type} ...";
    }
  }
  else {
    croak "Error: ${ClassName}->${ErrorMsg}: Fingerprint vector types $FingerprintsVectorA->{Type} is not valid...";
  }
}

# Process fingerprints vectors for similarity calculation by generating vectors
# containing ordered list of values...
#
sub _ProcessFingerprintsVectorsForSimilarityCalculation {
  my($ErrorMsg, $FingerprintsVectorA, $FingerprintsVectorB, $CalculationMode) = @_;

  $FingerprintsVectorA->{OrderedValuesRef} = undef; $FingerprintsVectorB->{OrderedValuesRef} = undef;
  $FingerprintsVectorA->{BitVector} = undef; $FingerprintsVectorB->{BitVector} = undef;

  if ($FingerprintsVectorA->{Type} =~ /^OrderedNumericalValues$/i) {
    _ProcessOrderedNumericalValuesFingerprintsVectorsForSimilarityCalculation($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($FingerprintsVectorA->{Type} =~ /^NumericalValues$/i) {
    _ProcessNumericalValuesFingerprintsVectorsForSimilarityCalculation($FingerprintsVectorA, $FingerprintsVectorB);
  }
  elsif ($FingerprintsVectorA->{Type} =~ /^AlphaNumericalValues$/i) {
    _ProcessAlphaNumericalValuesFingerprintsVectorsForSimilarityCalculation($FingerprintsVectorA, $FingerprintsVectorB);
  }
  else {
    croak "Error: ${ClassName}->${ErrorMsg}: Fingerprint vector types $FingerprintsVectorA->{Type} is not valid...";
  }
  if ($CalculationMode =~ /^BinaryForm$/i) {
    _TransformFinalOrderedValuesIntoBitVectorsForSimilarityCalculation($FingerprintsVectorA, $FingerprintsVectorB);
  }
}

# Process fingerprints vectors with ordered numerical values for similarity calculations...
#
sub _ProcessOrderedNumericalValuesFingerprintsVectorsForSimilarityCalculation {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;

  $FingerprintsVectorA->{OrderedValuesRef} = \@{$FingerprintsVectorA->{Values}};
  $FingerprintsVectorB->{OrderedValuesRef} = \@{$FingerprintsVectorB->{Values}};
}

# Process fingerprints vectors with numerical values for similarity calculations...
#
sub _ProcessNumericalValuesFingerprintsVectorsForSimilarityCalculation {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;

  # Set up unique IDs and values map for each fingerprint vector...
  my($Index, $Value, $ValueID, %UniqueFingerprintsVectorAValueIDValues, %UniqueFingerprintsVectorBValueIDValues, %UniqueFingerprintsVectorsValueIDs);

  %UniqueFingerprintsVectorAValueIDValues = ();
  %UniqueFingerprintsVectorBValueIDValues = ();
  %UniqueFingerprintsVectorsValueIDs = ();

  # Go over first vector...
  for $Index (0 .. $#{$FingerprintsVectorA->{ValueIDs}}) {
    $ValueID = $FingerprintsVectorA->{ValueIDs}[$Index];
    $Value = $FingerprintsVectorA->{Values}[$Index];
    if (exists $UniqueFingerprintsVectorAValueIDValues{$ValueID}) {
      $UniqueFingerprintsVectorAValueIDValues{$ValueID} += $Value;
    }
    else {
      $UniqueFingerprintsVectorAValueIDValues{$ValueID} = $Value;
    }
    if (!exists $UniqueFingerprintsVectorsValueIDs{$ValueID}) {
      $UniqueFingerprintsVectorsValueIDs{$ValueID} = 1;
    }
  }

  # Go over second vector...
  for $Index (0 .. $#{$FingerprintsVectorB->{ValueIDs}}) {
    $ValueID = $FingerprintsVectorB->{ValueIDs}[$Index];
    $Value = $FingerprintsVectorB->{Values}[$Index];
    if (exists $UniqueFingerprintsVectorBValueIDValues{$ValueID}) {
      $UniqueFingerprintsVectorBValueIDValues{$ValueID} += $Value;
    }
    else {
      $UniqueFingerprintsVectorBValueIDValues{$ValueID} = $Value;
    }
    if (!exists $UniqueFingerprintsVectorsValueIDs{$ValueID}) {
      $UniqueFingerprintsVectorsValueIDs{$ValueID} = 1;
    }
  }

  # Setup ordered values...
  my(@UniqueOrderedValueIDs, @OrderedValuesA, @OrderedValuesB);

  @UniqueOrderedValueIDs = ();
  @UniqueOrderedValueIDs = sort keys %UniqueFingerprintsVectorsValueIDs;

  @OrderedValuesA = ();
  @OrderedValuesA = map { exists $UniqueFingerprintsVectorAValueIDValues{$_} ? $UniqueFingerprintsVectorAValueIDValues{$_} : 0 } @UniqueOrderedValueIDs;

  @OrderedValuesB = ();
  @OrderedValuesB = map { exists $UniqueFingerprintsVectorBValueIDValues{$_} ? $UniqueFingerprintsVectorBValueIDValues{$_} : 0 } @UniqueOrderedValueIDs;

  $FingerprintsVectorA->{OrderedValuesRef} = \@OrderedValuesA;
  $FingerprintsVectorB->{OrderedValuesRef} = \@OrderedValuesB;
}

# Process fingerprints vectors with allpha numerical values for similarity calculations...
#
sub _ProcessAlphaNumericalValuesFingerprintsVectorsForSimilarityCalculation {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;

  # Set up unique IDs and values map for each vector...
  my($Index, $Value, $ValueID, %UniqueFingerprintsVectorAValuesCount, %UniqueFingerprintsVectorBValuesCount, %UniqueFingerprintsVectorsValues);

  %UniqueFingerprintsVectorAValuesCount = ();
  %UniqueFingerprintsVectorBValuesCount = ();
  %UniqueFingerprintsVectorsValues = ();

  # Go over first vector...
  for $Value (@{$FingerprintsVectorA->{Values}}) {
    if (exists $UniqueFingerprintsVectorAValuesCount{$Value}) {
      $UniqueFingerprintsVectorAValuesCount{$Value} += 1;
    }
    else {
      $UniqueFingerprintsVectorAValuesCount{$Value} = 1;
    }
    if (!exists $UniqueFingerprintsVectorsValues{$Value}) {
      $UniqueFingerprintsVectorsValues{$Value} = 1;
    }
  }

  # Go over second vector...
  for $Value (@{$FingerprintsVectorB->{Values}}) {
    if (exists $UniqueFingerprintsVectorBValuesCount{$Value}) {
      $UniqueFingerprintsVectorBValuesCount{$Value} += 1;
    }
    else {
      $UniqueFingerprintsVectorBValuesCount{$Value} = 1;
    }
    if (!exists $UniqueFingerprintsVectorsValues{$Value}) {
      $UniqueFingerprintsVectorsValues{$Value} = 1;
    }
  }

  # Setup ordered values...
  my(@UniqueOrderedValueIDs, @OrderedValuesA, @OrderedValuesB);

  @UniqueOrderedValueIDs = ();
  @UniqueOrderedValueIDs = sort keys %UniqueFingerprintsVectorsValues;

  @OrderedValuesA = ();
  @OrderedValuesA = map { exists $UniqueFingerprintsVectorAValuesCount{$_} ? $UniqueFingerprintsVectorAValuesCount{$_} : 0 } @UniqueOrderedValueIDs;

  @OrderedValuesB = ();
  @OrderedValuesB = map { exists $UniqueFingerprintsVectorBValuesCount{$_} ? $UniqueFingerprintsVectorBValuesCount{$_} : 0 } @UniqueOrderedValueIDs;

  $FingerprintsVectorA->{OrderedValuesRef} = \@OrderedValuesA;
  $FingerprintsVectorB->{OrderedValuesRef} = \@OrderedValuesB;

}

# Transform final ordered values array into a BitVector for similarity calculation...
#
sub _TransformFinalOrderedValuesIntoBitVectorsForSimilarityCalculation {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($Index, $Size, $BitVectorA, $BitVectorB, $SkipCheck);

  # Create bit vectors...
  $Size = scalar @{$FingerprintsVectorA->{OrderedValuesRef}};

  $FingerprintsVectorA->{BitVector} = new BitVector($Size);
  $FingerprintsVectorB->{BitVector} = new BitVector($Size);

  # Set bits...
  $SkipCheck = 1;
  for $Index (0 .. ($Size - 1)) {
    if ($FingerprintsVectorA->{OrderedValuesRef}[$Index]) {
      $FingerprintsVectorA->{BitVector}->SetBit($Index, $SkipCheck);
    }
    if ($FingerprintsVectorB->{OrderedValuesRef}[$Index]) {
      $FingerprintsVectorB->{BitVector}->SetBit($Index, $SkipCheck);
    }
  }
}

# Return sum of ordered vector values...
#
sub _GetSumOfFingerprintsOrderedValues {
  my($FingerprintVector) = @_;

  return StatisticsUtil::Sum($FingerprintVector->{OrderedValuesRef});
}

# Return sum of squared ordered vector values...
#
sub _GetSumOfSquaresOfFingerprintsOrderedValues {
  my($FingerprintVector) = @_;

  return StatisticsUtil::SumOfSquares($FingerprintVector->{OrderedValuesRef});
}

# Return sum of product of correponding ordered vector values...
#
sub _GetSumOfProductOfFingerprintsOrderedValues {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($Index, $SumProductXaiXbi);

  $SumProductXaiXbi = 0;
  for $Index (0 .. $#{$FingerprintsVectorA->{OrderedValuesRef}}) {
    $SumProductXaiXbi += $FingerprintsVectorA->{OrderedValuesRef}[$Index] * $FingerprintsVectorB->{OrderedValuesRef}[$Index];
  }
  return $SumProductXaiXbi;
}

# Return sum of absolute value of subtraction of correponding ordered vector values...
#
sub _GetSumOfAbsoluteValueOfSubtractionOfFingerprintsOrderedValues {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($Index, $SumAbsSubtractionXaiXbi);

  $SumAbsSubtractionXaiXbi = 0;
  for $Index (0 .. $#{$FingerprintsVectorA->{OrderedValuesRef}}) {
    $SumAbsSubtractionXaiXbi += abs($FingerprintsVectorA->{OrderedValuesRef}[$Index] - $FingerprintsVectorB->{OrderedValuesRef}[$Index]);
  }
  return $SumAbsSubtractionXaiXbi;
}

# Return sum of squares of subtraction of correponding ordered vector values...
#
sub _GetSumOfSquaresOfSubtractionOfFingerprintsOrderedValues {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($Index, $SumSquaresSubtractionXaiXbi);

  $SumSquaresSubtractionXaiXbi = 0;
  for $Index (0 .. $#{$FingerprintsVectorA->{OrderedValuesRef}}) {
    $SumSquaresSubtractionXaiXbi += ($FingerprintsVectorA->{OrderedValuesRef}[$Index] - $FingerprintsVectorB->{OrderedValuesRef}[$Index])**2;
  }
  return $SumSquaresSubtractionXaiXbi;
}

# Return sum of minimum of correponding ordered vector values...
#
sub _GetSumOfMinimumOfFingerprintsOrderdedValues {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($Index, $SumMinXaiXbi);

  $SumMinXaiXbi = 0;
  for $Index (0 .. $#{$FingerprintsVectorA->{OrderedValuesRef}}) {
    $SumMinXaiXbi += MathUtil::min($FingerprintsVectorA->{OrderedValuesRef}[$Index], $FingerprintsVectorB->{OrderedValuesRef}[$Index]);
  }
  return $SumMinXaiXbi;
}

# Return sum of maximum of correponding ordered vector values...
#
sub _GetSumOfMaximumOfFingerprintsOrderdedValues {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($Index, $SumMaxXaiXbi);

  $SumMaxXaiXbi = 0;
  for $Index (0 .. $#{$FingerprintsVectorA->{OrderedValuesRef}}) {
    $SumMaxXaiXbi += MathUtil::max($FingerprintsVectorA->{OrderedValuesRef}[$Index], $FingerprintsVectorB->{OrderedValuesRef}[$Index]);
  }
  return $SumMaxXaiXbi;
}

# Get number of Na, Nb and Nc bits in vector A and B for BinaryForm calculation...
#
sub _GetNumOfIndividualAndCommonSetBits ($$) {
  my($FingerprintsVectorA, $FingerprintsVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd, $FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $FingerprintsBitVectorA = $FingerprintsVectorA->{BitVector};
  $FingerprintsBitVectorB = $FingerprintsVectorB->{BitVector};

  # Number of bits set to "1" in A
  $Na = $FingerprintsBitVectorA->GetNumOfSetBits();

  # Number of bits set to "1" in B
  $Nb = $FingerprintsBitVectorB->GetNumOfSetBits();

  # Number of bits set to "1" in both A and B
  my($NcBitVector);
  $NcBitVector = $FingerprintsBitVectorA & $FingerprintsBitVectorB;
  $Nc = $NcBitVector->GetNumOfSetBits();

  return ($Na, $Nb, $Nc);
}

# Return a list of supported distance coefficients...
#
sub GetSupportedDistanceCoefficients () {

  return @DistanceCoefficients;
}

# Return a list of supported similarity coefficients...
#
sub GetSupportedSimilarityCoefficients () {

  return @SimilarityCoefficients;
}

# Return a list of supported distance and similarity coefficients...
#
sub GetSupportedDistanceAndSimilarityCoefficients () {
  my(@DistanceAndSimilarityCoefficients);

  @DistanceAndSimilarityCoefficients = ();
  push @DistanceAndSimilarityCoefficients, @DistanceCoefficients;
  push @DistanceAndSimilarityCoefficients, @SimilarityCoefficients;

  return sort @DistanceAndSimilarityCoefficients;
}

# Is it a fingerprints vector object?
sub IsFingerprintsVector ($) {
  my($Object) = @_;

  return _IsFingerprintsVector($Object);
}

# Is it a fingerprints vector object?
sub _IsFingerprintsVector {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Return a string containing vector values...
sub StringifyFingerprintsVector {
  my($This) = @_;
  my($FingerprintsVectorString);

  # Set type, values and value IDs...
  my($NumOfValues, $ValuesString, $NumOfValueIDs, $ValueIDsString, $MaxValuesToStringify);

  $NumOfValues = $This->GetNumOfValues();
  $MaxValuesToStringify = 500;

  if ($NumOfValues < $MaxValuesToStringify) {
    # Append all values...
    $ValuesString = $NumOfValues ? join ' ', @{$This->{Values}} : 'None';
  }
  else {
    # Truncate values...
    my($Index, @Values);
    for $Index (0 .. ($MaxValuesToStringify - 1)) {
      push @Values, $This->{Values}[$Index];
    }
    $ValuesString = join(' ', @Values) . " ...";
  }

  $NumOfValueIDs = $This->GetNumOfValueIDs();
  if ($NumOfValueIDs < $MaxValuesToStringify) {
    # Append all valueIDs...
    $ValueIDsString = $NumOfValueIDs ? join ' ', @{$This->{ValueIDs}} : 'None';
  }
  else {
    # Truncate value IDs...
    my($Index, @ValueIDs);
    @ValueIDs = ();
    for $Index (0 .. ($MaxValuesToStringify - 1)) {
      push @ValueIDs, $This->{ValueIDs}[$Index];
    }
    $ValueIDsString = join(' ', @ValueIDs) . " ...";
  }

  $FingerprintsVectorString = "Type: $This->{Type}; NumOfValues: $NumOfValues";
  if ($This->{Type} =~ /^(OrderedNumericalValues|NumericalValues)$/i) {
    my($NumOfNonZeroValues);
    $NumOfNonZeroValues = $This->GetNumOfNonZeroValues();
    $FingerprintsVectorString .= "; NumOfNonZeroValues: $NumOfNonZeroValues";
  }

  # Append all the values and value IDs...
  if ($NumOfValues < $MaxValuesToStringify) {
    $FingerprintsVectorString .= "; Values: <$ValuesString>; NumOfValueIDs: $NumOfValueIDs; ValueIDs: <$ValueIDsString>";
  }
  else {
    $FingerprintsVectorString .= "; Values (Truncated after $MaxValuesToStringify): <$ValuesString>; NumOfValueIDs: $NumOfValueIDs; ValueIDs (Truncated after $MaxValuesToStringify): <$ValueIDsString>";
  }

  return $FingerprintsVectorString;
}

1;

__END__

=head1 NAME

FingerprintsVector

=head1 SYNOPSIS

use Fingerprints::FingerprintsVector;

use Fingerprints::FingerprintsVector qw(:all);

=head1 DESCRIPTION

B<FingerprintsVector> class provides the following methods:

new, AddValueIDs, AddValues, CityBlockDistanceCoefficient,
CosineSimilarityCoefficient, CzekanowskiSimilarityCoefficient,
DiceSimilarityCoefficient, EuclideanDistanceCoefficient, GetDescription,
GetFingerprintsVectorString, GetID, GetIDsAndValuesPairsString,
GetIDsAndValuesString, GetNumOfNonZeroValues, GetNumOfValueIDs, GetNumOfValues,
GetSupportedDistanceAndSimilarityCoefficients, GetSupportedDistanceCoefficients,
GetSupportedSimilarityCoefficients, GetType, GetValue, GetValueID, GetValueIDs,
GetValueIDsString, GetValues, GetValuesAndIDsPairsString, GetValuesAndIDsString,
GetValuesString, GetVectorType, HammingDistanceCoefficient, IsFingerprintsVector,
JaccardSimilarityCoefficient, ManhattanDistanceCoefficient,
NewFromIDsAndValuesPairsString, NewFromIDsAndValuesString,
NewFromValuesAndIDsPairsString, NewFromValuesAndIDsString, NewFromValuesString,
OchiaiSimilarityCoefficient, SetDescription, SetID, SetType, SetValue, SetValueID,
SetValueIDs, SetValues, SetVectorType, SoergelDistanceCoefficient,
SorensonSimilarityCoefficient, StringifyFingerprintsVector,
TanimotoSimilarityCoefficient

The methods available to create fingerprints vector from strings and to calculate similarity
and distance coefficients between two vectors can also be invoked as class functions.

B<FingerprintsVector> class provides support to perform comparison between vectors
containing three different types of values:

Type I: OrderedNumericalValues

    o Size of two vectors are same
    o Vectors contain real values in a specific order. For example: MACCS keys
      count, Topological pharmacophore atom pairs and so on.

Type II: UnorderedNumericalValues

    o Size of two vectors might not be same
    o Vectors contain unordered real value identified by value IDs. For example:
      Topological atom pairs, Topological atom torsions and so on

Type III: AlphaNumericalValues

    o Size of two vectors might not be same
    o Vectors contain unordered alphanumerical values. For example: Extended
      connectivity fingerprints, atom neighborhood fingerprints.

Before performing similarity or distance calculations between vectors containing UnorderedNumericalValues
or AlphaNumericalValues, the vectors are transformed into vectors containing unique OrderedNumericalValues
using value IDs for UnorderedNumericalValues and values itself for AlphaNumericalValues.

Three forms of similarity and distance calculation between two vectors, specified using B<CalculationMode>
option, are supported: I<AlgebraicForm, BinaryForm or SetTheoreticForm>.

For I<BinaryForm>, the ordered list of processed final vector values containing the value or
count of each unique value type is simply converted into a binary vector containing 1s and 0s
corresponding to presence or absence of values before calculating similarity or distance between
two vectors.

For two fingerprint vectors A and B of same size containing OrderedNumericalValues, let:

    N = Number values in A or B

    Xa = Values of vector A
    Xb = Values of vector B

    Xai = Value of ith element in A
    Xbi = Value of ith element in B

   SUM = Sum of i over N values

For SetTheoreticForm of calculation between two vectors, let:

    SetIntersectionXaXb = SUM ( MIN ( Xai, Xbi ) )
    SetDifferenceXaXb = SUM ( Xai ) + SUM ( Xbi ) - SUM ( MIN ( Xai, Xbi ) )

For BinaryForm of calculation between two vectors, let:

    Na = Number of bits set to "1" in A = SUM ( Xai )
    Nb = Number of bits set to "1" in B = SUM ( Xbi )
    Nc = Number of bits set to "1" in both A and B = SUM ( Xai * Xbi )
    Nd = Number of bits set to "0" in both A and B
       = SUM ( 1 - Xai - Xbi + Xai * Xbi)

    N = Number of bits set to "1" or "0" in A or B = Size of A or B = Na + Nb - Nc + Nd

Additionally, for BinaryForm various values also correspond to:

    Na = | Xa |
    Nb = | Xb |
    Nc = | SetIntersectionXaXb |
    Nd = N - | SetDifferenceXaXb |

    | SetDifferenceXaXb | = N - Nd = Na + Nb - Nc + Nd - Nd = Na + Nb - Nc
                          =  | Xa | + | Xb | - | SetIntersectionXaXb |

Various similarity and distance coefficients [ Ref 40, Ref 62, Ref 64 ] for a pair of vectors A and B
in I<AlgebraicForm, BinaryForm and SetTheoreticForm> are defined as follows:

B<CityBlockDistance>: ( same as HammingDistance and ManhattanDistance)

I<AlgebraicForm>: SUM ( ABS ( Xai - Xbi ) )

I<BinaryForm>: ( Na - Nc ) + ( Nb - Nc ) = Na + Nb - 2 * Nc

I<SetTheoreticForm>: | SetDifferenceXaXb | - | SetIntersectionXaXb | = SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) )

B<CosineSimilarity>:  ( same as OchiaiSimilarityCoefficient)

I<AlgebraicForm>: SUM ( Xai * Xbi ) / SQRT ( SUM ( Xai ** 2) * SUM ( Xbi ** 2) )

I<BinaryForm>: Nc / SQRT ( Na * Nb)

I<SetTheoreticForm>: | SetIntersectionXaXb | / SQRT ( |Xa| * |Xb| ) = SUM ( MIN ( Xai, Xbi ) ) / SQRT ( SUM ( Xai ) * SUM ( Xbi ) )

B<CzekanowskiSimilarity>: ( same as DiceSimilarity and SorensonSimilarity)

I<AlgebraicForm>: ( 2 * ( SUM ( Xai * Xbi ) )  ) / ( SUM ( Xai ** 2) + SUM ( Xbi **2 ) )

I<BinaryForm>: 2 * Nc / ( Na + Nb )

I<SetTheoreticForm>: 2 * | SetIntersectionXaXb | / ( |Xa| + |Xb| ) = 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) / ( SUM ( Xai ) + SUM ( Xbi ) )

B<DiceSimilarity>: ( same as CzekanowskiSimilarity and SorensonSimilarity)

I<AlgebraicForm>: ( 2 * ( SUM ( Xai * Xbi ) )  ) / ( SUM ( Xai ** 2) + SUM ( Xbi **2 ) )

I<BinaryForm>: 2 * Nc / ( Na + Nb )

I<SetTheoreticForm>: 2 * | SetIntersectionXaXb | / ( |Xa| + |Xb| ) = 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) / ( SUM ( Xai ) + SUM ( Xbi ) )

B<EuclideanDistance>:

I<AlgebraicForm>: SQRT ( SUM ( ( ( Xai - Xbi ) ** 2 ) ) )

I<BinaryForm>: SQRT ( ( Na - Nc ) + ( Nb - Nc ) ) = SQRT ( Na + Nb - 2 * Nc )

I<SetTheoreticForm>: SQRT ( | SetDifferenceXaXb | - | SetIntersectionXaXb | ) = SQRT (  SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) )

B<HammingDistance>:  ( same as CityBlockDistance and ManhattanDistance)

I<AlgebraicForm>: SUM ( ABS ( Xai - Xbi ) )

I<BinaryForm>: ( Na - Nc ) + ( Nb - Nc ) = Na + Nb - 2 * Nc

I<SetTheoreticForm>: | SetDifferenceXaXb | - | SetIntersectionXaXb | = SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) )

B<JaccardSimilarity>: ( same as TanimotoSimilarity)

I<AlgebraicForm>:  SUM ( Xai * Xbi ) / ( SUM ( Xai ** 2 ) + SUM ( Xbi ** 2 ) - SUM ( Xai * Xbi ) )

I<BinaryForm>:  Nc / ( ( Na - Nc ) + ( Nb - Nc ) + Nc ) = Nc / ( Na + Nb - Nc )

I<SetTheoreticForm>: | SetIntersectionXaXb | / | SetDifferenceXaXb | = SUM ( MIN ( Xai, Xbi ) ) / (  SUM ( Xai ) + SUM ( Xbi ) - SUM ( MIN ( Xai, Xbi ) ) )

B<ManhattanDistance>:  ( same as CityBlockDistance and HammingDistance)

I<AlgebraicForm>: SUM ( ABS ( Xai - Xbi ) )

I<BinaryForm>: ( Na - Nc ) + ( Nb - Nc ) = Na + Nb - 2 * Nc

I<SetTheoreticForm>: | SetDifferenceXaXb | - | SetIntersectionXaXb | = SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) )

B<OchiaiSimilarity>:  ( same as CosineSimilarity)

I<AlgebraicForm>: SUM ( Xai * Xbi ) / SQRT ( SUM ( Xai ** 2) * SUM ( Xbi ** 2) )

I<BinaryForm>: Nc / SQRT ( Na * Nb)

I<SetTheoreticForm>: | SetIntersectionXaXb | / SQRT ( |Xa| * |Xb| ) = SUM ( MIN ( Xai, Xbi ) ) / SQRT ( SUM ( Xai ) * SUM ( Xbi ) )

B<SorensonSimilarity>: ( same as CzekanowskiSimilarity and DiceSimilarity)

I<AlgebraicForm>: ( 2 * ( SUM ( Xai * Xbi ) )  ) / ( SUM ( Xai ** 2) + SUM ( Xbi **2 ) )

I<BinaryForm>: 2 * Nc / ( Na + Nb )

I<SetTheoreticForm>: 2 * | SetIntersectionXaXb | / ( |Xa| + |Xb| ) = 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) / ( SUM ( Xai ) + SUM ( Xbi ) )

B<SoergelDistance>:

I<AlgebraicForm>:  SUM ( ABS ( Xai - Xbi ) ) / SUM ( MAX ( Xai, Xbi ) )

I<BinaryForm>: 1 - Nc / ( Na + Nb - Nc ) = ( Na + Nb - 2 * Nc ) / ( Na + Nb - Nc )

I<SetTheoreticForm>: ( | SetDifferenceXaXb | - | SetIntersectionXaXb | ) / | SetDifferenceXaXb | = ( SUM ( Xai ) + SUM ( Xbi ) - 2 * ( SUM ( MIN ( Xai, Xbi ) ) ) ) / ( SUM ( Xai ) + SUM ( Xbi ) - SUM ( MIN ( Xai, Xbi ) ) )

B<TanimotoSimilarity>:  ( same as JaccardSimilarity)

I<AlgebraicForm>:  SUM ( Xai * Xbi ) / ( SUM ( Xai ** 2 ) + SUM ( Xbi ** 2 ) - SUM ( Xai * Xbi ) )

I<BinaryForm>:  Nc / ( ( Na - Nc ) + ( Nb - Nc ) + Nc ) = Nc / ( Na + Nb - Nc )

I<SetTheoreticForm>: | SetIntersectionXaXb | / | SetDifferenceXaXb | = SUM ( MIN ( Xai, Xbi ) ) / (  SUM ( Xai ) + SUM ( Xbi ) - SUM ( MIN ( Xai, Xbi ) ) )

=head2 METHODS

=over 4

=item B<new>

    $FPVector = new Fingerprints::FingerprintsVector(%NamesAndValues);

Using specified I<FingerprintsVector> property names and values hash, B<new> method creates
a new object and returns a reference to newly created B<FingerprintsVectorsVector>
object. By default, the following properties are initialized:

    Type = ''
    @{Values} = ()
    @{ValuesIDs} = ()

Examples:

    $FPVector = new Fingerprints::FingerprintsVector('Type' => 'OrderedNumericalValues',
                                       'Values' => [1, 2, 3, 4]);
    $FPVector = new Fingerprints::FingerprintsVector('Type' => 'NumericalValues',
                                       'Values' => [10, 22, 33, 44],
                                       'ValueIDs' => ['ID1', 'ID2', 'ID3', 'ID4']);
    $FPVector = new Fingerprints::FingerprintsVector('Type' => 'AlphaNumericalValues',
                                       'Values' => ['a1', 2, 'a3', 4]);

=item B<AddValueIDs>

    $FingerprintsVector->AddValueIDs($ValueIDsRef);
    $FingerprintsVector->AddValueIDs(@ValueIDs);

Adds specified I<ValueIDs> to I<FingerprintsVector> and returns I<FingerprintsVector>.

=item B<AddValues>

    $FingerprintsVector->AddValues($ValuesRef);
    $FingerprintsVector->AddValues(@Values);
    $FingerprintsVector->AddValues($Vector);

Adds specified I<Values> to I<FingerprintsVector> and returns I<FingerprintsVector>.

=item B<CityBlockDistanceCoefficient>

    $Value = $FingerprintsVector->CityBlockDistanceCoefficient(
              $OtherFingerprintVector, [$CalculationMode, $SkipValuesCheck]);
    $Value = Fingerprints::FingerprintsVector::CityBlockDistanceCoefficient(
              $FingerprintsVectorA, $FingerprintVectorB,
              [$CalculationMode, $SkipValuesCheck]);

Returns value of I<CityBlock> distance coefficient between two I<FingerprintsVectors> using
optionally specified I<CalculationMode> and optional checking of vector values.

Possible I<CalculationMode> values: I<AlgebraicForm, BinaryForm or SetTheoreticForm>. Default
I<CalculationMode> value: I<AlgebraicForm>. Default I<SkipValuesCheck> value: I<0>.

=item B<CosineSimilarityCoefficient>

    $Value = $FingerprintsVector->CosineSimilarityCoefficient(
              $OtherFingerprintVector, [$CalculationMode, $SkipValuesCheck]);
    $Value = Fingerprints::FingerprintsVector::CosineSimilarityCoefficient(
              $FingerprintsVectorA, $FingerprintVectorB,
              [$CalculationMode, $SkipValuesCheck]);

Returns value of I<Cosine> similarity coefficient between two I<FingerprintsVectors> using
optionally specified I<CalculationMode> and optional checking of vector values.

Possible I<CalculationMode> values: I<AlgebraicForm, BinaryForm or SetTheoreticForm>. Default
I<CalculationMode> value: I<AlgebraicForm>. Default I<SkipValuesCheck> value: I<0>.

=item B<CzekanowskiSimilarityCoefficient>

    $Value = $FingerprintsVector->CzekanowskiSimilarityCoefficient(
              $OtherFingerprintVector, [$CalculationMode, $SkipValuesCheck]);
    $Value = Fingerprints::FingerprintsVector::CzekanowskiSimilarityCoefficient(
              $FingerprintsVectorA, $FingerprintVectorB,
              [$CalculationMode, $SkipValuesCheck]);

Returns value of I<Czekanowski> similarity coefficient between two I<FingerprintsVectors> using
optionally specified I<CalculationMode> and optional checking of vector values.

Possible I<CalculationMode> values: I<AlgebraicForm, BinaryForm or SetTheoreticForm>. Default
I<CalculationMode> value: I<AlgebraicForm>. Default I<SkipValuesCheck> value: I<0>.

=item B<DiceSimilarityCoefficient>

    $Value = $FingerprintsVector->DiceSimilarityCoefficient(
              $OtherFingerprintVector, [$CalculationMode, $SkipValuesCheck]);
    $Value = Fingerprints::FingerprintsVector::DiceSimilarityCoefficient(
              $FingerprintsVectorA, $FingerprintVectorB,
              [$CalculationMode, $SkipValuesCheck]);

Returns value of I<Dice> similarity coefficient between two I<FingerprintsVectors> using
optionally specified I<CalculationMode> and optional checking of vector values.

Possible I<CalculationMode> values: I<AlgebraicForm, BinaryForm or SetTheoreticForm>. Default
I<CalculationMode> value: I<AlgebraicForm>. Default I<SkipValuesCheck> value: I<0>.

=item B<EuclideanDistanceCoefficient>

    $Value = $FingerprintsVector->EuclideanDistanceCoefficient(
              $OtherFingerprintVector, [$CalculationMode, $SkipValuesCheck]);
    $Value = Fingerprints::FingerprintsVector::EuclideanDistanceCoefficient(
              $FingerprintsVectorA, $FingerprintVectorB,
              [$CalculationMode, $SkipValuesCheck]);

Returns value of I<Euclidean> distance coefficient between two I<FingerprintsVectors> using
optionally specified I<CalculationMode> and optional checking of vector values.

Possible I<CalculationMode> values: I<AlgebraicForm, BinaryForm or SetTheoreticForm>. Default
I<CalculationMode> value: I<AlgebraicForm>. Default I<SkipValuesCheck> value: I<0>.

=item B<GetDescription>

    $Description = $FingerprintsVector->GetDescription();

Returns a string containing description of fingerprints vector.

=item B<GetFingerprintsVectorString>

    $FPString = $FingerprintsVector->GetFingerprintsVectorString($Format);

Returns a B<FingerprintsString> containing vector values and/or IDs in I<FingerprintsVector>
corresponding to specified I<Format>.

Possible I<Format> values: I<IDsAndValuesString, IDsAndValues, IDsAndValuesPairsString,
IDsAndValuesPairs, ValuesAndIDsString, ValuesAndIDs, ValuesAndIDsPairsString, ValuesAndIDsPairs,
ValueIDsString, ValueIDs, ValuesString, or Values>.

=item B<GetID>

    $ID = $FingerprintsVector->GetID();

Returns I<ID> of I<FingerprintsVector>.

=item B<GetVectorType>

    $VectorType = $FingerprintsVector->GetVectorType();

Returns I<VectorType> of I<FingerprintsVector>.

=item B<GetIDsAndValuesPairsString>

    $IDsValuesPairsString = $FingerprintsVector->GetIDsAndValuesPairsString();

Returns I<FingerprintsVector> value IDs and values as space delimited ID/value pair
string.

=item B<GetIDsAndValuesString>

    $IDsValuesString = $FingerprintsVector->GetIDsAndValuesString();

Returns I<FingerprintsVector> value IDs and values as string containing space delimited IDs followed by
values with semicolon as IDs and values delimiter.

=item B<GetNumOfNonZeroValues>

    $NumOfNonZeroValues = $FingerprintsVector->GetNumOfNonZeroValues();

Returns number of non-zero values in I<FingerprintsVector>.

=item B<GetNumOfValueIDs>

    $NumOfValueIDs = $FingerprintsVector->GetNumOfValueIDs();

Returns number of value IDs I<FingerprintsVector>.

=item B<GetNumOfValues>

    $NumOfValues = $FingerprintsVector->GetNumOfValues();

Returns number of values I<FingerprintsVector>.

=item B<GetSupportedDistanceAndSimilarityCoefficients>

    @SupportedDistanceAndSimilarityCoefficientsReturn =
        Fingerprints::FingerprintsVector::GetSupportedDistanceAndSimilarityCoefficients();

Returns an array containing names of supported distance and similarity coefficients.

=item B<GetSupportedDistanceCoefficients>

    @SupportedDistanceCoefficientsReturn =
        Fingerprints::FingerprintsVector::GetSupportedDistanceCoefficients();

Returns an array containing names of supported disyance coefficients.

=item B<GetSupportedSimilarityCoefficients>

    @SupportedSimilarityCoefficientsReturn =
        Fingerprints::FingerprintsVector::GetSupportedSimilarityCoefficients();

Returns an array containing names of supported similarity coefficients.

=item B<GetType>

    $VectorType = $FingerprintsVector->GetType();

Returns I<FingerprintsVector> vector type.

=item B<GetValue>

    $Value = $FingerprintsVector->GetValue($Index);

Returns fingerprints vector B<Value> specified using I<Index> starting at 0.

=item B<GetValueID>

    $ValueID = $FingerprintsVector->GetValueID();

Returns fingerprints vector B<ValueID> specified using I<Index> starting at 0.

=item B<GetValueIDs>

    $ValueIDs = $FingerprintsVector->GetValueIDs();
    @ValueIDs = $FingerprintsVector->GetValueIDs();

Returns fingerprints vector B<ValueIDs> as an array or reference to an array.

=item B<GetValueIDsString>

    $ValueIDsString = $FingerprintsVector->GetValueIDsString();

Returns fingerprints vector B<ValueIDsString> with value IDs delimited by space.

=item B<GetValues>

    $ValuesRef = $FingerprintsVector->GetValues();
    @Values = $FingerprintsVector->GetValues();

Returns fingerprints vector B<Values> as an array or reference to an array.

=item B<GetValuesAndIDsPairsString>

    $ValuesIDsPairsString = $FingerprintsVector->GetValuesAndIDsPairsString();

Returns I<FingerprintsVector> value and value IDs as space delimited ID/value pair
string.

=item B<GetValuesAndIDsString>

    $ValuesIDsString = $FingerprintsVector->GetValuesAndIDsString();

Returns I<FingerprintsVector> values and value IDs as string containing space delimited IDs followed by
values with semicolon as IDs and values delimiter.

=item B<GetValuesString>

    $Return = $FingerprintsVector->GetValuesString();

Returns I<FingerprintsVector> values as space delimited string.

=item B<HammingDistanceCoefficient>

    $Value = $FingerprintsVector->HammingDistanceCoefficient(
              $OtherFingerprintVector, [$CalculationMode, $SkipValuesCheck]);
    $Value = Fingerprints::FingerprintsVector::HammingDistanceCoefficient(
              $FingerprintsVectorA, $FingerprintVectorB,
              [$CalculationMode, $SkipValuesCheck]);

Returns value of I<Hamming> distance coefficient between two I<FingerprintsVectors> using
optionally specified I<CalculationMode> and optional checking of vector values.

Possible I<CalculationMode> values: I<AlgebraicForm, BinaryForm or SetTheoreticForm>. Default
I<CalculationMode> value: I<AlgebraicForm>. Default I<SkipValuesCheck> value: I<0>.

=item B<IsFingerprintsVector>

    $Status = Fingerprints::FingerprintsVector::IsFingerprintsVector($Object);

Returns 1 or 0 based on whether I<Object> is a I<FingerprintsVector>.

=item B<JaccardSimilarityCoefficient>

    $Value = $FingerprintsVector->JaccardSimilarityCoefficient(
              $OtherFingerprintVector, [$CalculationMode, $SkipValuesCheck]);
    $Value = Fingerprints::FingerprintsVector::JaccardSimilarityCoefficient(
              $FingerprintsVectorA, $FingerprintVectorB,
              [$CalculationMode, $SkipValuesCheck]);

Returns value of I<Jaccard> similarity coefficient between two I<FingerprintsVectors> using
optionally specified I<CalculationMode> and optional checking of vector values.

Possible I<CalculationMode> values: I<AlgebraicForm, BinaryForm or SetTheoreticForm>. Default
I<CalculationMode> value: I<AlgebraicForm>. Default I<SkipValuesCheck> value: I<0>.

=item B<ManhattanDistanceCoefficient>

    $Value = $FingerprintsVector->ManhattanDistanceCoefficient(
              $OtherFingerprintVector, [$CalculationMode, $SkipValuesCheck]);
    $Value = Fingerprints::FingerprintsVector::ManhattanDistanceCoefficient(
              $FingerprintsVectorA, $FingerprintVectorB,
              [$CalculationMode, $SkipValuesCheck]);

Returns value of I<Manhattan> distance coefficient between two I<FingerprintsVectors> using
optionally specified I<CalculationMode> and optional checking of vector values.

Possible I<CalculationMode> values: I<AlgebraicForm, BinaryForm or SetTheoreticForm>. Default
I<CalculationMode> value: I<AlgebraicForm>. Default I<SkipValuesCheck> value: I<0>.

=item B<NewFromIDsAndValuesPairsString>

    $FingerprintsVector = $FingerprintsVector->NewFromIDsAndValuesPairsString(
                          $ValuesType, $IDsAndValuesPairsString);
    $FingerprintsVector = Fingerprints::FingerprintsVector::NewFromIDsAndValuesPairsString(
                          $ValuesType, $IDsAndValuesPairsString);

Creates a new I<FingerprintsVector> of I<ValuesType> using I<IDsAndValuesPairsString> containing
space delimited value IDs and values pairs and returns new B<FingerprintsVector> object.
Possible I<ValuesType> values: I<OrderedNumericalValues, NumericalValues, or AlphaNumericalValues>.

=item B<NewFromIDsAndValuesString>

    $FingerprintsVector = $FingerprintsVector->NewFromIDsAndValuesString(
                          $ValuesType, $IDsAndValuesString);
    $FingerprintsVector = Fingerprints::FingerprintsVector::NewFromIDsAndValuesString(
                          $ValuesType, $IDsAndValuesString);

Creates a new I<FingerprintsVector> of I<ValuesType> using I<IDsAndValuesString> containing
semicolon delimited value IDs string followed by values strings and returns new B<FingerprintsVector>
object. The values within value and value IDs tring are delimited by spaces. Possible I<ValuesType>
values: I<OrderedNumericalValues, NumericalValues, or AlphaNumericalValues>.

=item B<NewFromValuesAndIDsPairsString>

    $FingerprintsVector = $FingerprintsVector->NewFromValuesAndIDsPairsString(
                          $ValuesType, $ValuesAndIDsPairsString);
    $FingerprintsVector = Fingerprints::FingerprintsVector::NewFromValuesAndIDsPairsString(
                          $ValuesType, $ValuesAndIDsPairsString);

Creates a new I<FingerprintsVector> of I<ValuesType> using I<ValuesAndIDsPairsString> containing
space delimited value and value IDs pairs and returns new B<FingerprintsVector> object.
Possible I<ValuesType> values: I<OrderedNumericalValues, NumericalValues, or AlphaNumericalValues>.

=item B<NewFromValuesAndIDsString>

    $FingerprintsVector = $FingerprintsVector->NewFromValuesAndIDsString(
                          $ValuesType, $IDsAndValuesString);
    $FingerprintsVector = Fingerprints::FingerprintsVector::NewFromValuesAndIDsString(
                          $ValuesType, $IDsAndValuesString);

Creates a new I<FingerprintsVector> of I<ValuesType> using I<ValuesAndIDsString> containing
semicolon delimited values string followed by value IDs strings and returns new B<FingerprintsVector>
object. The values within values and value IDs tring are delimited by spaces. Possible I<ValuesType>
values: I<OrderedNumericalValues, NumericalValues, or AlphaNumericalValues>.

=item B<NewFromValuesString>

    $FingerprintsVector = $FingerprintsVector->NewFromValuesString(
                          $ValuesType, $ValuesString);
    $FingerprintsVector = Fingerprints::FingerprintsVector::NewFromValuesString(
                          $ValuesType, $ValuesString);

Creates a new I<FingerprintsVector> of I<ValuesType> using I<ValuesString> containing space
delimited values string and returns new B<FingerprintsVector> object. The values within values
and value IDs tring are delimited by spaces. Possible I<ValuesType> values: I<OrderedNumericalValues,
NumericalValues, or AlphaNumericalValues>.

=item B<OchiaiSimilarityCoefficient>

    $Value = $FingerprintsVector->OchiaiSimilarityCoefficient(
              $OtherFingerprintVector, [$CalculationMode, $SkipValuesCheck]);
    $Value = Fingerprints::FingerprintsVector::OchiaiSimilarityCoefficient(
              $FingerprintsVectorA, $FingerprintVectorB,
              [$CalculationMode, $SkipValuesCheck]);

Returns value of I<Ochiai> similarity coefficient between two I<FingerprintsVectors> using
optionally specified I<CalculationMode> and optional checking of vector values.

Possible I<CalculationMode> values: I<AlgebraicForm, BinaryForm or SetTheoreticForm>. Default
I<CalculationMode> value: I<AlgebraicForm>. Default I<SkipValuesCheck> value: I<0>.

=item B<SetDescription>

    $FingerprintsVector->SetDescription($Description);

Sets I<Description> of fingerprints vector and returns I<FingerprintsVector>.

=item B<SetID>

    $FingerprintsVector->SetID($ID);

Sets I<ID> of fingerprints vector and returns I<FingerprintsVector>.

=item B<SetVectorType>

    $FingerprintsVector->SetVectorType($VectorType);

Sets I<VectorType> of fingerprints vector and returns I<FingerprintsVector>.

=item B<SetType>

    $FingerprintsVector->SetType($Type);

Sets I<FingerprintsVector> values I<Type> and returns I<FingerprintsVector>. Possible I<Type>
values: I<OrderedNumericalValues, NumericalValues, or AlphaNumericalValues>.

During calculation of similarity and distance coefficients between two I<FingerprintsVectors>, the
following conditions apply to vector type, size, value and value IDs:

    o For OrderedNumericalValues type, both vectors must be of the same size
      and contain similar types of numerical values in the same order.

    o For NumericalValues type, vector value IDs for both vectors must be
      specified; however, their size and order of IDs and numerical values may
      be different. For each vector, value IDs must correspond to vector values.

    o For AlphaNumericalValues type, vectors may contain both numerical and
      alphanumerical values and their sizes may be different.

=item B<SetValue>

    $FingerprintsVector->SetValue($Index, $Value, [$SkipIndexCheck]);

Sets a I<FingerprintsVector> value specified by I<Index> starting at 0 to I<Value> along with
optional index range check and returns I<FingerprintsVector>.

=item B<SetValueID>

    $FingerprintsVector->SetValueID($Index, $ValueID, [$SkipIndexCheck]);

Sets a I<FingerprintsVector> value ID specified by I<Index> starting at 0 to I<ValueID> along with
optional index range check and returns I<FingerprintsVector>.

=item B<SetValueIDs>

    $FingerprintsVector->SetValueIDs($ValueIDsRef);
    $FingerprintsVector->SetValueIDs(@ValueIDs);

Sets I<FingerprintsVector> value IDs to specified I<ValueIDs> and returns I<FingerprintsVector>.

=item B<SetValues>

    $FingerprintsVector->SetValues($ValuesRef);
    $FingerprintsVector->SetValues(@Values);

Sets I<FingerprintsVector> value to specified I<Values> and returns I<FingerprintsVector>.

=item B<SoergelDistanceCoefficient>

    $Value = $FingerprintsVector->SoergelDistanceCoefficient(
              $OtherFingerprintVector, [$CalculationMode, $SkipValuesCheck]);
    $Value = Fingerprints::FingerprintsVector::SoergelDistanceCoefficient(
              $FingerprintsVectorA, $FingerprintVectorB,
              [$CalculationMode, $SkipValuesCheck]);

Returns value of I<Soergel> distance coefficient between two I<FingerprintsVectors> using
optionally specified I<CalculationMode> and optional checking of vector values.

Possible I<CalculationMode> values: I<AlgebraicForm, BinaryForm or SetTheoreticForm>. Default
I<CalculationMode> value: I<AlgebraicForm>. Default I<SkipValuesCheck> value: I<0>.

=item B<SorensonSimilarityCoefficient>

    $Value = $FingerprintsVector->SorensonSimilarityCoefficient(
              $OtherFingerprintVector, [$CalculationMode, $SkipValuesCheck]);
    $Value = Fingerprints::FingerprintsVector::SorensonSimilarityCoefficient(
              $FingerprintsVectorA, $FingerprintVectorB,
              [$CalculationMode, $SkipValuesCheck]);

Returns value of I<Sorenson> similarity coefficient between two I<FingerprintsVectors> using
optionally specified I<CalculationMode> and optional checking of vector values.

Possible I<CalculationMode> values: I<AlgebraicForm, BinaryForm or SetTheoreticForm>. Default
I<CalculationMode> value: I<AlgebraicForm>. Default I<SkipValuesCheck> value: I<0>.

=item B<TanimotoSimilarityCoefficient>

    $Value = $FingerprintsVector->TanimotoSimilarityCoefficient(
              $OtherFingerprintVector, [$CalculationMode, $SkipValuesCheck]);
    $Value = Fingerprints::FingerprintsVector::TanimotoSimilarityCoefficient(
              $FingerprintsVectorA, $FingerprintVectorB,
              [$CalculationMode, $SkipValuesCheck]);

Returns value of I<Tanimoto> similarity coefficient between two I<FingerprintsVectors> using
optionally specified I<CalculationMode> and optional checking of vector values.

Possible I<CalculationMode> values: I<AlgebraicForm, BinaryForm or SetTheoreticForm>. Default
I<CalculationMode> value: I<AlgebraicForm>. Default I<SkipValuesCheck> value: I<0>.

=item B<StringifyFingerprintsVector>

    $String = $FingerprintsVector->StringifyFingerprintsVector();

Returns a string containing information about I<FingerprintsVector> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

BitVector.pm, FingerprintsStringUtil.pm, FingerprintsBitVector.pm, Vector.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
