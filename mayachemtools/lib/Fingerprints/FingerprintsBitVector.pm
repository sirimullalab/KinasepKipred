package Fingerprints::FingerprintsBitVector;
#
# File: FingerprintsBitVector.pm
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
use BitVector;
use MathUtil;
use TextUtil ();

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(BitVector Exporter);

# Similiarity coefficients...
my(@SimilarityCoefficients) = qw(BaroniUrbaniSimilarityCoefficient BuserSimilarityCoefficient CosineSimilarityCoefficient DiceSimilarityCoefficient DennisSimilarityCoefficient ForbesSimilarityCoefficient FossumSimilarityCoefficient HamannSimilarityCoefficient JacardSimilarityCoefficient Kulczynski1SimilarityCoefficient Kulczynski2SimilarityCoefficient MatchingSimilarityCoefficient McConnaugheySimilarityCoefficient OchiaiSimilarityCoefficient PearsonSimilarityCoefficient RogersTanimotoSimilarityCoefficient RussellRaoSimilarityCoefficient SimpsonSimilarityCoefficient SkoalSneath1SimilarityCoefficient SkoalSneath2SimilarityCoefficient SkoalSneath3SimilarityCoefficient TanimotoSimilarityCoefficient TverskySimilarityCoefficient YuleSimilarityCoefficient WeightedTanimotoSimilarityCoefficient WeightedTverskySimilarityCoefficient);

# New from string...
my(@NewFromString) = qw(NewFromBinaryString NewFromHexadecimalString NewFromRawBinaryString);

@EXPORT = qw(IsFingerprintsBitVector);
@EXPORT_OK = qw(GetSupportedSimilarityCoefficients @NewFromString @SimilarityCoefficients);

%EXPORT_TAGS = (
		new => [@NewFromString],
		coefficients => [@SimilarityCoefficients],
		all  => [@EXPORT, @EXPORT_OK]
	       );

# Setup class variables...
my($ClassName);
_InitializeClass();

use overload '""' => 'StringifyFingerprintsBitVector';

# Class constructor...
sub new {
  my($Class, $Size) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new($Size);
  bless $This, ref($Class) || $Class;
  $This->_InitializeFingerprintsBitVector($Size);

  return $This;
}

# Initialize object data...
#
# Note:
#  . The class, BitVector, used to derive this class provides all the functionality to
#    manipulate bits.
#  . Irrespective of specified size, Perl functions used to handle bit data in
#    BitVector class automatically sets the size to the next nearest power of 2.
#    SpecifiedSize is used by this class to process any aribitray size during similarity
#    coefficient calculations.
#
sub _InitializeFingerprintsBitVector {
  my($This, $Size) = @_;

  if (!defined $Size) {
    croak "Error: ${ClassName}->new: FingerprintsBitVector object instantiated without specifying its size ...";
  }
  if ($Size <=0) {
    croak "Error: ${ClassName}->new: Fingerprints bit vector size, $Size, must be a positive integer...";
  }

  # Specified size of fingerprints...
  $This->{SpecifiedSize} = $Size;

}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Set specified size...
#
# Notes:
#   Irrespective of specified size, Perl functions used to handle bit data in
#   BitVector class automatically sets the size to the next nearest power of 2.
#   SpecifiedSize is used by this class to process any aribitray size during similarity
#   coefficient calculations.
#
sub SetSpecifiedSize {
  my($This, $SpecifiedSize) = @_;

  if (!($SpecifiedSize > 0 && $SpecifiedSize <= $This->{Size})) {
    croak "Error: ${ClassName}->SetSpecifiedSize: Specified size, $SpecifiedSize, is not valid:  It must be > 0 && <= ", $This->GetSize()," ...";
  }
  $This->{SpecifiedSize} = $SpecifiedSize;
}

# Get specified size...
sub GetSpecifiedSize {
  my($This) = @_;

  return $This->{SpecifiedSize};
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

  return exists $This->{VectorType} ? $This->{VectorType} : 'FingerprintsBitVector';
}

# Create a new fingerprints bit vector using binary string. This functionality can be
# either invoked as a class function or an object method.
#
sub NewFromBinaryString ($;$) {
  my($FirstParameter, $SecondParameter, $ThirdParameter) = @_;

  if (_IsFingerprintsBitVector($FirstParameter)) {
    return _NewFingerptinsBitVectorFromString('Binary', $SecondParameter, $ThirdParameter);
  }
  else {
    return _NewFingerptinsBitVectorFromString( 'Binary', $FirstParameter, $SecondParameter);
  }
}

# Create a new fingerprints bit vector using hexadecimal string. This functionality can be
# either invoked as a class function or an object method.
#
sub NewFromHexadecimalString ($;$) {
  my($FirstParameter, $SecondParameter, $ThirdParameter) = @_;

  if (_IsFingerprintsBitVector($FirstParameter)) {
    return _NewFingerptinsBitVectorFromString('Hexadecimal', $SecondParameter, $ThirdParameter);
  }
  else {
    return _NewFingerptinsBitVectorFromString( 'Hexadecimal', $FirstParameter, $SecondParameter);
  }
}

# Create a new fingerprints bit vector using octal string. This functionality can be
# either invoked as a class function or an object method.
#
#
sub NewFromOctalString ($) {
  croak "Error: ${ClassName}->NewFromOctalString: Creation of fingerprits bit vector from an octal string is not supported ...";
}

# Create a new fingerprints bit vector using decimal string. This functionality can be
# either invoked as a class function or an object method.
#
sub NewFromDecimalString ($;$) {
  croak "Error: ${ClassName}->NewFromDecimalString: Creation of fingerprits bit vector from a decimal string is not supported ...";
}

# Create a new fingerprints bit vector using raw binary string. This functionality can be
# either invoked as a class function or an object method.
#
sub NewFromRawBinaryString ($;$) {
  my($FirstParameter, $SecondParameter, $ThirdParameter) = @_;

  if (_IsFingerprintsBitVector($FirstParameter)) {
    return _NewFingerptinsBitVectorFromString('RawBinary', $SecondParameter, $ThirdParameter);
  }
  else {
    return _NewFingerptinsBitVectorFromString( 'RawBinary', $FirstParameter, $SecondParameter);
  }
}

# Create a new fingerprints bit vector from a string...
#
#
sub _NewFingerptinsBitVectorFromString ($$;$) {
  my($Format, $String, $BitsOrder) = @_;
  my($FingerprintsBitVector, $Size);

  $Size = BitVector::_CalculateStringSizeInBits($Format, $String);

  $FingerprintsBitVector = new Fingerprints::FingerprintsBitVector($Size);
  $FingerprintsBitVector->_SetBitsAsString($Format, $String, $BitsOrder);

  return $FingerprintsBitVector;
}

# Get fingerprint bits as a hexadecimal string...
#
sub GetBitsAsHexadecimalString {
  my($This, $BitsOrder) = @_;

  return $This->_GetFingerprintBitsAsString('Hexadecimal', $BitsOrder);
}

# Get fingerprint bits as an octal string...
#
sub GetBitsAsOctalString {
  my($This, $BitsOrder) = @_;

  croak "Error: ${ClassName}->GetBitsAsOctalString: Retrieval of fingerprits bits as an octal string is not supported ...";
}

# Get fingerprint bits as an decimal string...
#
sub GetBitsAsDecimalString {
  my($This, $BitsOrder) = @_;

  croak "Error: ${ClassName}->GetBitsAsOctalString: Retrieval of fingerprits bits as a decimal string is not supported ...";
}

# Get fingerprint bits as a binary string conatning 1s and 0s...
#
sub GetBitsAsBinaryString {
  my($This, $BitsOrder) = @_;

  return $This->_GetFingerprintBitsAsString('Binary', $BitsOrder);
}

# Get fingerprint bits as a binary string conatning 1s and 0s...
#
sub GetBitsAsRawBinaryString {
  my($This) = @_;

  return $This->_GetFingerprintBitsAsString('RawBinary');
}

# Return fingerprint bits as a string...
#
sub _GetFingerprintBitsAsString {
  my($This, $Format, $BitsOrder) = @_;

  $BitsOrder = (defined($BitsOrder) && $BitsOrder) ? $BitsOrder : 'Ascending';

  return $This->_GetBitsAsString($Format, $BitsOrder);
}

# Is it a fingerprints bit vector object?
sub IsFingerprintsBitVector ($) {
  my($Object) = @_;

  return _IsFingerprintsBitVector($Object);
}

# Is it a fingerprints bit vector object?
sub _IsFingerprintsBitVector {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Return a list of supported similarity coefficients...
sub GetSupportedSimilarityCoefficients () {

  return @SimilarityCoefficients;
}

# Get bit density for fingerprints bit vector corresponding to on bits...
#
sub GetFingerprintsBitDensity {
  my($This) = @_;
  my($BitDensity);

  $BitDensity = $This->GetDensityOfSetBits();

  return round($BitDensity, 2);
}

# Fold fingerprints bit vector by recursively reducing its size by half untill size is less than or equal to
# specified size...
#
sub FoldFingerprintsBitVectorBySize {
  my($This, $Size) = @_;

  if (!($Size > 0 && $Size <= $This->GetSize())) {
    croak "Error: ${ClassName}->FoldFingerprintsBitVectorBySize: Specified size, $Size, is not valid:  It must be > 0 && <= ", $This->GetSize()," ...";
  }

  if ($This->GetSize() <= $Size) {
    return $This;
  }
  return $This->_FoldFingerprintsBitVector('BySize', $Size);
}

# Fold fingerprints bit vector by recursively reducing its size by half untill bit density of set bits is greater than
#  or equal to specified density...
#
sub FoldFingerprintsBitVectorByDensity {
  my($This, $Density) = @_;

  if (!($Density > 0 && $Density <= 1)) {
    croak "Error: ${ClassName}->FoldFingerprintsBitVectorByDensity: Specified bit density, $Density, is not valid:  It must be > 0 && <= 1 ...";
  }

  if ($This->GetDensityOfSetBits() >= $Density) {
    return $This;
  }
  return $This->_FoldFingerprintsBitVector('ByDensity', $Density);
}

# Fold fingerprints bit vector using size or density and return folded fingerprint bit vector...
#
sub _FoldFingerprintsBitVector {
  my($This, $Mode, $Value) = @_;

  # Fold upto size of 8 bits...
  if ($This->GetSize() <= 8) {
    return $This;
  }

  # Check size or density....
  if ($Mode =~ /^BySize$/i) {
    if ($This->GetSize() <= $Value) {
      return $This;
    }
  }
  elsif ($Mode =~ /^ByDensity$/i) {
    if ($This->GetDensityOfSetBits() >= $Value) {
      return $This;
    }
  }
  else {
    return $This;
  }

  # Recursively reduce its size by half...
  my($FirstHalfBinaryString, $SecondHalfBinaryString, $FirstHalfFingerprintsBitVector, $SecondHalfFingerprintsBitVector, $FoldedFingerprintsBitVector, $BinaryString, $StringLength);

  $BinaryString = $This->GetBitsAsBinaryString();
  $StringLength = length $BinaryString;

  $FirstHalfBinaryString = substr($BinaryString, 0, $StringLength/2);
  $SecondHalfBinaryString = substr($BinaryString, $StringLength/2);

  $FirstHalfFingerprintsBitVector = NewFromBinaryString($FirstHalfBinaryString);
  $SecondHalfFingerprintsBitVector = NewFromBinaryString($SecondHalfBinaryString);

  $FoldedFingerprintsBitVector = $FirstHalfFingerprintsBitVector | $SecondHalfFingerprintsBitVector;

  return $FoldedFingerprintsBitVector->_FoldFingerprintsBitVector($Mode, $Value);
}

# Is first bit vector subset of second bit vector?
#
# For a bit vector to be a subset of another bit vector, both vectors must be of
# the same size and the bit positions set in first vector must also be set in the
# secons bit vector.
#
# This functionality can be either invoked as a class function or an object method.
#
sub IsSubSet ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;

  if ($FingerprintsBitVectorA->GetSize() != $FingerprintsBitVectorB->GetSize()) {
    return 0;
  }
  my($AndFingerprintsBitVector);

  $AndFingerprintsBitVector = $FingerprintsBitVectorA & $FingerprintsBitVectorB;

  return ($FingerprintsBitVectorA->GetNumOfSetBits() == $AndFingerprintsBitVector->GetNumOfSetBits()) ? 1 : 0;
}

# Return a string containing vector values...
sub StringifyFingerprintsBitVector {
  my($This) = @_;
  my($FingerprintsBitVectorString);

  # BitVector size information...
  #
  if ($This->{SpecifiedSize} != $This->GetSize()) {
    $FingerprintsBitVectorString = "SpecifiedSize: " . $This->{SpecifiedSize} . "; BitVectorSize: " . $This->GetSize();
  }
  else {
    $FingerprintsBitVectorString = "BitVectorSize: " . $This->GetSize();
  }
  my($NumOfSetBits, $BitDensity);
  $NumOfSetBits = $This->GetNumOfSetBits();
  $BitDensity = $This->GetFingerprintsBitDensity();

  $FingerprintsBitVectorString .= "; NumOfOnBits: $NumOfSetBits; BitDensity: $BitDensity";

  # BitVector values...
  $FingerprintsBitVectorString .= "; BitVector: " . $This->StringifyBitVector();

  return $FingerprintsBitVectorString;
}

# For two fingerprints bit vectors A and B of same size, let:
#
#  Na = Number of bits set to "1" in A
#  Nb = Number of bits set to "1" in B
#  Nc = Number of bits set to "1" in both A and B
#  Nd = Number of bits set to "0" in both A and B
#
#  Nt = Number of bits set to "1" or "0" in A or B = Size of A or B = Na + Nb - Nc + Nd
#
#  Na - Nc = Number of bits set to "1" in A but not in B
#  Nb - Nc = Number of bits set to "1" in B but not in A
#
# Various similarity coefficients [ Ref 40 - 42 ] for a pair of bit vectors A and B are
# defined as follows:
#
# . BaroniUrbani: ( SQRT( Nc * Nd ) + Nc ) / (  SQRT ( Nc * Nd ) + Nc + ( Na - Nc )  + ( Nb - Nc ) ) ( same as Buser )
#
# . Buser: ( SQRT ( Nc * Nd ) + Nc ) / (  SQRT ( Nc * Nd ) + Nc + ( Na - Nc )  + ( Nb - Nc ) ) ( same as BaroniUrbani )
#
# . Cosine: Nc / SQRT ( Na * Nb ) (same as Ochiai)
#
# . Dice: (2 * Nc) / ( Na + Nb )
#
# . Dennis: ( Nc * Nd - ( ( Na - Nc ) * ( Nb - Nc ) ) ) / SQRT ( Nt * Na * Nb)
#
# . Forbes: ( Nt * Nc ) / ( Na * Nb )
#
# . Fossum: ( Nt * ( ( Nc - 1/2 ) ** 2 ) / ( Na * Nb )
#
# . Hamann: ( ( Nc + Nd ) - ( Na - Nc ) - ( Nb - Nc ) ) / Nt
#
# . Jaccard: Nc /  ( ( Na - Nc) + ( Nb - Nc ) + Nc ) = Nc / ( Na + Nb - Nc ) (same as Tanimoto)
#
# . Kulczynski1: Nc / ( ( Na - Nc ) + ( Nb - Nc) ) = Nc / ( Na + Nb - 2Nc )
#
# . Kulczynski2: ( ( Nc / 2 ) * ( 2 * Nc + ( Na - Nc ) + ( Nb - Nc) ) ) / ( ( Nc + ( Na - Nc ) ) * ( Nc + ( Nb - Nc ) ) ) = 0.5 * ( Nc / Na + Nc / Nb )
#
# . Matching: ( Nc + Nd ) / Nt
#
# . McConnaughey: ( Nc ** 2 - ( Na - Nc ) * ( Nb - Nc) ) / (  Na * Nb )
#
# . Ochiai: Nc / SQRT ( Na * Nb ) (same as Cosine)
#
# . Pearson: ( ( Nc * Nd ) - ( ( Na - Nc ) * ( Nb - Nc ) ) / SQRT ( Na * Nb * (  Na - Nc + Nd ) * ( Nb - Nc + Nd ) )
#
# . RogersTanimoto: ( Nc + Nd ) / ( ( Na - Nc)  + ( Nb  - Nc) + Nt) = ( Nc + Nd ) / ( Na  + Nb  - 2Nc + Nt)
#
# . RussellRao: Nc / Nt
#
# . Simpson: Nc / MIN ( Na, Nb)
#
# . SkoalSneath1: Nc / ( Nc + 2 * ( Na - Nc)  + 2 * ( Nb - Nc) ) = Nc / ( 2 * Na + 2 * Nb - 3 * Nc )
#
# . SkoalSneath2: ( 2 * Nc + 2 * Nd ) / ( Nc + Nd + Nt )
#
# . SkoalSneath3: ( Nc + Nd ) / ( ( Na - Nc ) + ( Nb - Nc ) ) = ( Nc + Nd ) / ( Na + Nb - 2 * Nc  )
#
# . Tanimoto: Nc /  ( ( Na - Nc) + ( Nb - Nc ) + Nc ) = Nc / ( Na + Nb - Nc ) (same as Jaccard)
#
# . Tversky: Nc / ( alpha * ( Na - Nc ) + ( 1 - alpha) * ( Nb - Nc) + Nc ) = Nc / ( alpha * ( Na - Nb )  + Nb)
#
# . Yule: ( ( Nc * Nd ) - ( ( Na - Nc ) * ( Nb - Nc ) ) ) / ( ( Nc * Nd ) + ( ( Na - Nc ) * ( Nb - Nc ) )  )
#
#
# Values of Tanimoto/Jaccard and Tversky coefficients are dependent on only those bit which
# are set to "1" in both A and B. In order to take into account all bit positions, modified versions
# of Tanimoto [ Ref. 42 ] and Tversky [  Ref. 43 ] have been developed.
#
# Let:
#
#  Na' = Number of bits set to "0" in A
#  Nb' = Number of bits set to "0" in B
#  Nc' = Number of bits set to "0" in both A and B
#
# . Tanimoto': Nc' /  ( ( Na' - Nc') + ( Nb' - Nc' ) + Nc' ) = Nc' / ( Na' + Nb' - Nc' )
#
# . Tversky': Nc' / ( alpha * ( Na' - Nc' ) + ( 1 - alpha) * ( Nb' - Nc' ) + Nc' ) = Nc' / ( alpha * ( Na' - Nb' )  + Nb')
#
# Then:
#
# . WeightedTanimoto = beta * Tanimoto + (1 - beta) * Tanimoto'
#
# . WeightedTversky = beta * Tversky + (1 - beta) * Tversky'
#
#

# Calculate BaroniUrbani similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub BaroniUrbaniSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;

  return BuserSimilarityCoefficient($FingerprintsBitVectorA, $FingerprintsBitVectorB);
}

# Calculate Buser similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub BuserSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd, $Nt, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nd = _GetNumOfCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nt = $Na + $Nb - $Nc + $Nd;

  $Numerator = sqrt($Nc*$Nd) + $Nc;
  $Denominator = sqrt($Nc*$Nd) + ($Na - $Nc)  + ($Nb - $Nc ) + $Nc;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Cosine similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub CosineSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $Numerator = $Nc;
  $Denominator = sqrt($Na*$Nb);

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Dice similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub DiceSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $Numerator = 2*$Nc;
  $Denominator = $Na + $Nb;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Dennis similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub DennisSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd, $Nt, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nd = _GetNumOfCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nt = $Na + $Nb - $Nc + $Nd;

  $Numerator = $Nc*$Nd - (($Na - $Nc)*($Nb - $Nc));
  $Denominator = sqrt($Nt*$Na*$Nb);

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Forbes similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub ForbesSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd, $Nt, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nd = _GetNumOfCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nt = $Na + $Nb - $Nc + $Nd;

  $Numerator = $Nt*$Nc;
  $Denominator = $Na*$Nb;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Fossum similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub FossumSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd, $Nt, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nd = _GetNumOfCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nt = $Na + $Nb - $Nc + $Nd;

  $Numerator =  $Nt*(($Nc - 0.5)** 2);
  $Denominator =  $Na*$Nb ;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Hamann similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub HamannSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd, $Nt, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nd = _GetNumOfCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nt = $Na + $Nb - $Nc + $Nd;

  $Numerator =  ($Nc + $Nd ) - ($Na - $Nc) - ($Nb - $Nc) ;
  $Denominator = $Nt;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Jacard similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub JacardSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;

  return TanimotoSimilarityCoefficient($FingerprintsBitVectorA, $FingerprintsBitVectorB);
}

# Calculate Kulczynski1 similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub Kulczynski1SimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $Numerator = $Nc;
  $Denominator =  $Na + $Nb - 2*$Nc;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Kulczynski2 similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub Kulczynski2SimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $Numerator = 0.5*($Na*$Nc + $Nb*$Nc);
  $Denominator = $Na*$Nb;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Matching similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub MatchingSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd, $Nt, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nd = _GetNumOfCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nt = $Na + $Nb - $Nc + $Nd;

  $Numerator =  $Nc + $Nd;
  $Denominator = $Nt;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate McConnaughey similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub McConnaugheySimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $Numerator =  $Nc**2 - (($Na - $Nc)*($Nb - $Nc));
  $Denominator = $Na*$Nb ;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Ochiai similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub OchiaiSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;

  return CosineSimilarityCoefficient($FingerprintsBitVectorA, $FingerprintsBitVectorB);
}

# Calculate Pearson similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub PearsonSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd, $Nt, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nd = _GetNumOfCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nt = $Na + $Nb - $Nc + $Nd;

  $Numerator = ($Nc*$Nd ) - (($Na - $Nc)*($Nb - $Nc));
  $Denominator =  sqrt($Na*$Nb*($Na - $Nc + $Nd )*($Nb - $Nc + $Nd));

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate RogersTanimoto similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub RogersTanimotoSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd, $Nt, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nd = _GetNumOfCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nt = $Na + $Nb - $Nc + $Nd;

  $Numerator = $Nc + $Nd;
  $Denominator =  ($Na - $Nc)  + ($Nb  - $Nc) + $Nt;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate RussellRao similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub RussellRaoSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd, $Nt, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nd = _GetNumOfCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nt = $Na + $Nb - $Nc + $Nd;

  $Numerator = $Nc;
  $Denominator = $Nt;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Simpson similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub SimpsonSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $Numerator = $Nc;
  $Denominator =  min($Na, $Nb);

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate SkoalSneath1 similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub SkoalSneath1SimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $Numerator = $Nc;
  $Denominator = $Nc + 2*($Na - $Nc)  + 2*($Nb - $Nc);

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate SkoalSneath2 similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub SkoalSneath2SimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd, $Nt, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nd = _GetNumOfCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nt = $Na + $Nb - $Nc + $Nd;

  $Numerator = 2*$Nc + 2*$Nd  ;
  $Denominator = $Nc + $Nd + $Nt ;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate SkoalSneath3 similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub SkoalSneath3SimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd, $Nt, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nd = _GetNumOfCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nt = $Na + $Nb - $Nc + $Nd;

  $Numerator =  $Nc + $Nd;
  $Denominator = ($Na - $Nc) + ($Nb - $Nc ) ;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Tanimoto similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub TanimotoSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $Numerator = $Nc;
  $Denominator = $Na + $Nb - $Nc;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Tversky similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub TverskySimilarityCoefficient ($$$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB, $Alpha) = @_;
  my($Na, $Nb, $Nc, $Numerator, $Denominator);

  if (!(defined($Alpha) && ($Alpha >= 0 && $Alpha <= 1))) {
    croak "Error: ${ClassName}->TverskySimilarityCoefficient: Alpha parameters must be defined and its value must be >=0 and <=1 ...";
  }

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $Numerator = $Nc;
  $Denominator =  $Alpha*($Na - $Nb )  + $Nb;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate Yule similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub YuleSimilarityCoefficient ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd, $Nt, $Numerator, $Denominator);

  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nd = _GetNumOfCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);
  $Nt = $Na + $Nb - $Nc + $Nd;

  $Numerator = ($Nc*$Nd) - (($Na - $Nc)*($Nb - $Nc)) ;
  $Denominator = ($Nc*$Nd) + (($Na - $Nc)*($Nb - $Nc))  ;

  return  $Denominator ? ($Numerator/$Denominator) : 0;
}

# Calculate WeightedTanimoto similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub WeightedTanimotoSimilarityCoefficient ($$$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB, $Beta) = @_;
  my($Na, $Nb, $Nc, $TanimotoForSetBits, $TanimotoForClearBits, $Numerator, $Denominator, $WeightedTanimoto);

  if (!(defined($Beta) && ($Beta >= 0 && $Beta <= 1))) {
    croak "Error: ${ClassName}->WeightedTanimotoSimilarityCoefficient: Beta parameters must be defined and its value must be >=0 and <=1 ...";
  }

  # Get Tanimoto for set bits...
  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $Numerator = $Nc;
  $Denominator = $Na + $Nb - $Nc;
  $TanimotoForSetBits = $Denominator ? ($Numerator/$Denominator) : 0;

  # Get Tanimoto for clear bits...
  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $Numerator = $Nc;
  $Denominator = $Na + $Nb - $Nc;
  $TanimotoForClearBits = $Denominator ? ($Numerator/$Denominator) : 0;

  $WeightedTanimoto = $Beta*$TanimotoForSetBits + (1 - $Beta)*$TanimotoForClearBits;

  return  $WeightedTanimoto;
}

# Calculate WeightedTversky similarity coefficient for two same size bit vectors.
#
# This functionality can be either invoked as a class function or an object method.
#
sub WeightedTverskySimilarityCoefficient ($$$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB, $Alpha, $Beta) = @_;
  my($Na, $Nb, $Nc, $TverskyForSetBits, $TverskyForClearBits, $Numerator, $Denominator, $WeightedTversky);

  if (!(defined($Alpha) && ($Alpha >= 0 && $Alpha <= 1))) {
    croak "Error: ${ClassName}->WeightedTverskySimilarityCoefficient: Alpha parameters must be defined and its value must be >=0 and <=1 ...";
  }
  if (!(defined($Beta) && ($Beta >= 0 && $Beta <= 1))) {
    croak "Error: ${ClassName}->WeightedTverskySimilarityCoefficient: Beta parameters must be defined and its value must be >=0 and <=1 ...";
  }

  # Get Tversky for set bits...
  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonSetBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $Numerator = $Nc;
  $Denominator =  $Alpha*($Na - $Nb )  + $Nb;
  $TverskyForSetBits =  $Denominator ? ($Numerator/$Denominator) : 0;

  # Get Tversky for clear bits...
  ($Na, $Nb, $Nc) = _GetNumOfIndividualAndCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  $Numerator = $Nc;
  $Denominator =  $Alpha*($Na - $Nb )  + $Nb;
  $TverskyForClearBits =  $Denominator ? ($Numerator/$Denominator) : 0;

  $WeightedTversky = $Beta*$TverskyForSetBits + (1 - $Beta)*$TverskyForClearBits;

  return  $WeightedTversky;
}

# Get number of Na, Nb and Nc bits in bit vector A and B to be used for similarity coefficient calculations...
#
sub _GetNumOfIndividualAndCommonSetBits ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd);

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

# Get number of Nd bits in bit vector A and B to be used for similarity coefficient calculations...
#
sub _GetNumOfCommonClearBits ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Nd, $NdBitVector);

  #  Number of bits set to "0" in both A and B
  $NdBitVector = ~$FingerprintsBitVectorA & ~$FingerprintsBitVectorB;
  $Nd = $NdBitVector->GetNumOfSetBits();

  # Correct for number of clear bits used for padding...
  if (_IsNumOfClearBitsCorrectionRequired($FingerprintsBitVectorA)) {
    $Nd = $Nd - _GetNumOfClearBitsCorrection($FingerprintsBitVectorA);
  }
  elsif (_IsNumOfClearBitsCorrectionRequired($FingerprintsBitVectorB)) {
    $Nd = $Nd - _GetNumOfClearBitsCorrection($FingerprintsBitVectorB);
  }

  return $Nd;
}

# Get number of Na, Nb and Nc bits in bit vector A and B to be used for similarity coefficient calculations...
#
sub _GetNumOfIndividualAndCommonClearBits ($$) {
  my($FingerprintsBitVectorA, $FingerprintsBitVectorB) = @_;
  my($Na, $Nb, $Nc, $Nd);

  # Number of bits set to "0" in A
  $Na = $FingerprintsBitVectorA->GetNumOfClearBits();

  # Correct for number of clear bits used for padding...
  if (_IsNumOfClearBitsCorrectionRequired($FingerprintsBitVectorA)) {
    $Na = $Na - _GetNumOfClearBitsCorrection($FingerprintsBitVectorA);
  }

  # Number of bits set to "0" in B
  $Nb = $FingerprintsBitVectorB->GetNumOfClearBits();

  # Correct for number of clear bits used for padding...
  if (_IsNumOfClearBitsCorrectionRequired($FingerprintsBitVectorB)) {
    $Nb = $Nb - _GetNumOfClearBitsCorrection($FingerprintsBitVectorB);
  }

  # Number of bits set to "0" in both A and B
  $Nc = _GetNumOfCommonClearBits($FingerprintsBitVectorA, $FingerprintsBitVectorB);

  return ($Na, $Nb, $Nc);
}

# Irrespective of specified size, Perl functions used to handle bit data data in
# BitVector class automatically sets the size to the next nearest power of 2
# and clear the extra bits.
#
# SpecifiedSize is used by this class to process any aribitray size during similarity
# coefficient calculations.
#
# Assuming the FingerprintsBitBector class only manipulates bits upto specified
# size, a correction for the extra bits added by BitVector class needs to be applied
# to number of clear bits.
#
sub _GetNumOfClearBitsCorrection {
  my($FingerprintsBitVector) = @_;

  return ($FingerprintsBitVector->{Size} - $FingerprintsBitVector->{SpecifiedSize});
}

# Is number of clear bits correction required?
#
sub _IsNumOfClearBitsCorrectionRequired {
  my($FingerprintsBitVector) = @_;

  return ($FingerprintsBitVector->{Size} > $FingerprintsBitVector->{SpecifiedSize}) ? 1 : 0;
}


1;

__END__

=head1 NAME

FingerprintsBitVector

=head1 SYNOPSIS

use Fingerprints::FingerprintsBitVector;

use Fingerprints::FingerprintsBitVector qw(:coefficients);

use Fingerprints::FingerprintsBitVector qw(:all);

=head1 DESCRIPTION

B<FingerprintsBitVector> class provides the following methods:

new, BaroniUrbaniSimilarityCoefficient, BuserSimilarityCoefficient,
CosineSimilarityCoefficient, DennisSimilarityCoefficient,
DiceSimilarityCoefficient, FoldFingerprintsBitVectorByDensity,
FoldFingerprintsBitVectorBySize, ForbesSimilarityCoefficient,
FossumSimilarityCoefficient, GetBitsAsBinaryString, GetBitsAsDecimalString,
GetBitsAsHexadecimalString, GetBitsAsOctalString, GetBitsAsRawBinaryString,
GetDescription, GetFingerprintsBitDensity, GetID, GetSpecifiedSize,
GetSupportedSimilarityCoefficients, GetVectorType, HamannSimilarityCoefficient,
IsFingerprintsBitVector, IsSubSet, JacardSimilarityCoefficient,
Kulczynski1SimilarityCoefficient, Kulczynski2SimilarityCoefficient,
MatchingSimilarityCoefficient, McConnaugheySimilarityCoefficient,
NewFromBinaryString, NewFromDecimalString, NewFromHexadecimalString,
NewFromOctalString, NewFromRawBinaryString, OchiaiSimilarityCoefficient,
PearsonSimilarityCoefficient, RogersTanimotoSimilarityCoefficient,
RussellRaoSimilarityCoefficient, SetDescription, SetID, SetSpecifiedSize,
SetVectorType, SimpsonSimilarityCoefficient, SkoalSneath1SimilarityCoefficient,
SkoalSneath2SimilarityCoefficient, SkoalSneath3SimilarityCoefficient,
StringifyFingerprintsBitVector, TanimotoSimilarityCoefficient,
TverskySimilarityCoefficient, WeightedTanimotoSimilarityCoefficient,
WeightedTverskySimilarityCoefficient, YuleSimilarityCoefficient

The methods available to create fingerprints bit vector from strings and to calculate similarity
coefficient between two bit vectors can also be invoked as class functions.

B<FingerprintsBitVector> class is derived from B<BitVector> class which provides the functionality
to manipulate bits.

For two fingerprints bit vectors A and B of same size, let:

    Na = Number of bits set to "1" in A
    Nb = Number of bits set to "1" in B
    Nc = Number of bits set to "1" in both A and B
    Nd = Number of bits set to "0" in both A and B

    Nt = Number of bits set to "1" or "0" in A or B (Size of A or B)
    Nt = Na + Nb - Nc + Nd

    Na - Nc = Number of bits set to "1" in A but not in B
    Nb - Nc = Number of bits set to "1" in B but not in A

Then, various similarity coefficients [ Ref. 40 - 42 ] for a pair of bit vectors A and B are
defined as follows:

BaroniUrbani: ( SQRT( Nc * Nd ) + Nc ) / (  SQRT ( Nc * Nd ) + Nc + ( Na - Nc )  + ( Nb - Nc ) ) ( same as Buser )

Buser: ( SQRT ( Nc * Nd ) + Nc ) / (  SQRT ( Nc * Nd ) + Nc + ( Na - Nc )  + ( Nb - Nc ) ) ( same as BaroniUrbani )

Cosine: Nc / SQRT ( Na * Nb ) (same as Ochiai)

Dice: (2 * Nc) / ( Na + Nb )

Dennis: ( Nc * Nd - ( ( Na - Nc ) * ( Nb - Nc ) ) ) / SQRT ( Nt * Na * Nb)

Forbes: ( Nt * Nc ) / ( Na * Nb )

Fossum: ( Nt * ( ( Nc - 1/2 ) ** 2 ) / ( Na * Nb )

Hamann: ( ( Nc + Nd ) - ( Na - Nc ) - ( Nb - Nc ) ) / Nt

Jaccard: Nc /  ( ( Na - Nc) + ( Nb - Nc ) + Nc ) = Nc / ( Na + Nb - Nc ) (same as Tanimoto)

Kulczynski1: Nc / ( ( Na - Nc ) + ( Nb - Nc) ) = Nc / ( Na + Nb - 2Nc )

Kulczynski2: ( ( Nc / 2 ) * ( 2 * Nc + ( Na - Nc ) + ( Nb - Nc) ) ) / ( ( Nc + ( Na - Nc ) ) * ( Nc + ( Nb - Nc ) ) )
= 0.5 * ( Nc / Na + Nc / Nb )

Matching: ( Nc + Nd ) / Nt

McConnaughey: ( Nc ** 2 - ( Na - Nc ) * ( Nb - Nc) ) / (  Na * Nb )

Ochiai: Nc / SQRT ( Na * Nb ) (same as Cosine)

Pearson: ( ( Nc * Nd ) - ( ( Na - Nc ) * ( Nb - Nc ) ) / SQRT ( Na * Nb * (  Na - Nc + Nd ) * ( Nb - Nc + Nd ) )

RogersTanimoto: ( Nc + Nd ) / ( ( Na - Nc)  + ( Nb  - Nc) + Nt) = ( Nc + Nd ) / ( Na  + Nb  - 2Nc + Nt)

RussellRao: Nc / Nt

Simpson: Nc / MIN ( Na, Nb)

SkoalSneath1: Nc / ( Nc + 2 * ( Na - Nc)  + 2 * ( Nb - Nc) ) = Nc / ( 2 * Na + 2 * Nb - 3 * Nc )

SkoalSneath2: ( 2 * Nc + 2 * Nd ) / ( Nc + Nd + Nt )

SkoalSneath3: ( Nc + Nd ) / ( ( Na - Nc ) + ( Nb - Nc ) ) = ( Nc + Nd ) / ( Na + Nb - 2 * Nc  )

Tanimoto: Nc /  ( ( Na - Nc) + ( Nb - Nc ) + Nc ) = Nc / ( Na + Nb - Nc ) (same as Jaccard)

Tversky: Nc / ( alpha * ( Na - Nc ) + ( 1 - alpha) * ( Nb - Nc) + Nc ) = Nc / ( alpha * ( Na - Nb )  + Nb)

Yule: ( ( Nc * Nd ) - ( ( Na - Nc ) * ( Nb - Nc ) ) ) / ( ( Nc * Nd ) + ( ( Na - Nc ) * ( Nb - Nc ) )  )

The values of Tanimoto/Jaccard and Tversky coefficients are dependent on only those bit which
are set to "1" in both A and B. In order to take into account all bit positions, modified versions
of Tanimoto [ Ref. 42 ] and Tversky [  Ref. 43 ] have been developed.

Let:

    Na' = Number of bits set to "0" in A
    Nb' = Number of bits set to "0" in B
    Nc' = Number of bits set to "0" in both A and B

Tanimoto': Nc' /  ( ( Na' - Nc') + ( Nb' - Nc' ) + Nc' ) = Nc' / ( Na' + Nb' - Nc' )

Tversky': Nc' / ( alpha * ( Na' - Nc' ) + ( 1 - alpha) * ( Nb' - Nc' ) + Nc' ) = Nc' / ( alpha * ( Na' - Nb' )  + Nb')

Then:

WeightedTanimoto = beta * Tanimoto + (1 - beta) * Tanimoto'

WeightedTversky = beta * Tversky + (1 - beta) * Tversky'

=head2 METHODS

=over 4

=item B<new>

    $NewFPBitVector = new Fingerprints::FingerprintsBitVector($Size);

Creates a new I<FingerprintsBitVector> object of size I<Size> and returns  newly created
B<FingerprintsBitVector>. Bit numbers range from 0 to 1 less than I<Size>.

=item B<BaroniUrbaniSimilarityCoefficient>

    $Value = $FingerprintsBitVector->BaroniUrbaniSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::
              BaroniUrbaniSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<BaroniUrbani> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<BuserSimilarityCoefficient>

    $Value = $FingerprintsBitVector->BuserSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::BuserSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Buser> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<CosineSimilarityCoefficient>

    $Value = $FingerprintsBitVector->CosineSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::CosineSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Cosine> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<DennisSimilarityCoefficient>

    $Value = $FingerprintsBitVector->DennisSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::DennisSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Dennis> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<DiceSimilarityCoefficient>

    $Value = $FingerprintsBitVector->DiceSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::DiceSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Dice> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<FoldFingerprintsBitVectorByDensity>

    $FingerprintsBitVector->FoldFingerprintsBitVectorByDensity($Density);

Folds I<FingerprintsBitVector> by recursively reducing its size by half until bit density of set bits is
greater than or equal to specified I<Density> and returns folded I<FingerprintsBitVector>.

=item B<FoldFingerprintsBitVectorBySize>

    $FingerprintsBitVector->FoldFingerprintsBitVectorBySize($Size);

Folds I<FingerprintsBitVector> by recursively reducing its size by half until size is less than or equal to
specified I<Size> and returns folded I<FingerprintsBitVector>.

=item B<ForbesSimilarityCoefficient>

    $Value = $FingerprintsBitVector->ForbesSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::ForbesSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Forbes> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<FossumSimilarityCoefficient>

    $Value = $FingerprintsBitVector->FossumSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::FossumSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Fossum> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<GetBitsAsBinaryString>

    $BinaryASCIIString = $FingerprintsBitVector->GetBitsAsBinaryString();

Returns fingerprints as a binary ASCII string containing 0s and 1s.

=item B<GetBitsAsHexadecimalString>

    $HexadecimalString = $FingerprintsBitVector->GetBitsAsHexadecimalString();

Returns fingerprints as a hexadecimal string.

=item B<GetBitsAsRawBinaryString>

    $RawBinaryString = $FingerprintsBitVector->GetBitsAsRawBinaryString();

Returns fingerprints as a raw binary string containing packed bit values for each byte.

=item B<GetDescription>

    $Description = $FingerprintsBitVector->GetDescription();

Returns a string containing description of fingerprints bit vector.

=item B<GetFingerprintsBitDensity>

    $BitDensity = $FingerprintsBitVector->GetFingerprintsBitDensity();

Returns I<BitDensity> of I<FingerprintsBitVector> corresponding to bits set to 1s.

=item B<GetID>

    $ID = $FingerprintsBitVector->GetID();

Returns I<ID> of I<FingerprintsBitVector>.

=item B<GetVectorType>

    $VectorType = $FingerprintsBitVector->GetVectorType();

Returns I<VectorType> of I<FingerprintsBitVector>.

=item B<GetSpecifiedSize>

    $Size = $FingerprintsBitVector->GetSpecifiedSize();

Returns value of specified size for bit vector.

=item B<GetSupportedSimilarityCoefficients>

    @SimilarityCoefficient =
       Fingerprints::FingerprintsBitVector::GetSupportedSimilarityCoefficients();

Returns an array containing names of supported similarity coefficients.

=item B<HamannSimilarityCoefficient>

    $Value = $FingerprintsBitVector->HamannSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::HamannSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Hamann> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<IsFingerprintsBitVector>

    $Status = Fingerprints::FingerprintsBitVector::
              IsFingerprintsBitVector($Object);

Returns 1 or 0 based on whether I<Object> is a B<FingerprintsBitVector> object.

=item B<IsSubSet>

    $Status = $FingerprintsBitVector->IsSubSet($OtherFPBitVector);
    $Status = Fingerprints::FingerprintsBitVector::IsSubSet(
              $FPBitVectorA, $FPBitVectorB);

Returns 1 or 0 based on whether first firngerprints bit vector is a subset of second
fingerprints bit vector.

For a bit vector to be a subset of another bit vector, both vectors must be of
the same size and the bit positions set in first vector must also be set in the
second bit vector.

=item B<JacardSimilarityCoefficient>

    $Value = $FingerprintsBitVector->JacardSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::JacardSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Jacard> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<Kulczynski1SimilarityCoefficient>

    $Value = $FingerprintsBitVector->Kulczynski1SimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::
              Kulczynski1SimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Kulczynski1> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<Kulczynski2SimilarityCoefficient>

    $Value = $FingerprintsBitVector->Kulczynski2SimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::
              Kulczynski2SimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Kulczynski2> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<MatchingSimilarityCoefficient>

    $Value = $FingerprintsBitVector->MatchingSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::
              MatchingSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Matching> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<McConnaugheySimilarityCoefficient>

    $Value = $FingerprintsBitVector->McConnaugheySimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::
              McConnaugheySimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<McConnaughey> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<NewFromBinaryString>

    $NewFPBitVector = $FingerprintsBitVector->NewFromBinaryString(
                      $BinaryString);
    $NewFPBitVector = Fingerprints::FingerprintsBitVector::NewFromBinaryString(
                      $BinaryString);

Creates a new I<FingerprintsBitVector> using I<BinaryString> and returns new
B<FingerprintsBitVector> object.

=item B<NewFromHexadecimalString>

    $NewFPBitVector = $FingerprintsBitVector->NewFromHexadecimalString(
                      $HexdecimalString);
    $NewFPBitVector = Fingerprints::FingerprintsBitVector::
                    NewFromHexadecimalString(
                      $HexdecimalString);

Creates a new I<FingerprintsBitVector> using I<HexdecimalString> and returns new
B<FingerprintsBitVector> object.

=item B<NewFromRawBinaryString>

    $NewFPBitVector = $FingerprintsBitVector->NewFromRawBinaryString(
                      $RawBinaryString);
    $NewFPBitVector = Fingerprints::FingerprintsBitVector::
                      NewFromRawBinaryString(
                      $RawBinaryString);

Creates a new I<FingerprintsBitVector> using I<RawBinaryString> and returns new
B<FingerprintsBitVector> object.

=item B<OchiaiSimilarityCoefficient>

    $Value = $FingerprintsBitVector->OchiaiSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::OchiaiSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Ochiai> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<PearsonSimilarityCoefficient>

    $Value = $FingerprintsBitVector->PearsonSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::PearsonSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Pearson> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<RogersTanimotoSimilarityCoefficient>

    $Value = $FingerprintsBitVector->RogersTanimotoSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::
              RogersTanimotoSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<RogersTanimoto> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<RussellRaoSimilarityCoefficient>

    $Value = $FingerprintsBitVector->RussellRaoSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::
              RussellRaoSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<RussellRao> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<SetSpecifiedSize>

    $FingerprintsBitVector->SetSpecifiedSize($Size);

Sets specified size for fingerprints bit vector.

Irrespective of specified size, Perl functions used to handle bit data in B<BitVector> class
automatically sets the size to the next nearest power of 2. I<SpecifiedSize> is used by
B<FingerprintsBitVector> class to process any aribitray size during similarity coefficient calculations.

=item B<SetDescription>

    $FingerprintsBitVector->SetDescription($Description);

Sets I<Description> of fingerprints bit vector and returns I<FingerprintsBitVector>.

=item B<SetID>

    $FingerprintsBitVector->SetID($ID);

Sets I<ID> of fingerprints bit vector and returns I<FingerprintsBitVector>.

=item B<SetVectorType>

    $FingerprintsBitVector->SetVectorType($VectorType);

Sets I<VectorType> of fingerprints bit vector and returns I<FingerprintsBitVector>.

=item B<SimpsonSimilarityCoefficient>

    $Value = $FingerprintsBitVector->SimpsonSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::SimpsonSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Simpson> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<SkoalSneath1SimilarityCoefficient>

    $Value = $FingerprintsBitVector->SkoalSneath1SimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::
              SkoalSneath1SimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<SkoalSneath1> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<SkoalSneath2SimilarityCoefficient>

    $Value = $FingerprintsBitVector->SkoalSneath2SimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::
              SkoalSneath2SimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<SkoalSneath2> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<SkoalSneath3SimilarityCoefficient>

    $Value = $FingerprintsBitVector->SkoalSneath3SimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::
              SkoalSneath3SimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<SkoalSneath3> similarity coefficient between two same size I<FingerprintsBitVectors>

=item B<StringifyFingerprintsBitVector>

    $String = $FingerprintsBitVector->StringifyFingerprintsBitVector();

Returns a string containing information about I<FingerprintsBitVector> object.

=item B<TanimotoSimilarityCoefficient>

    $Value = $FingerprintsBitVector->TanimotoSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::
              TanimotoSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Tanimoto> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<TverskySimilarityCoefficient>

    $Value = $FingerprintsBitVector->TverskySimilarityCoefficient(
              $OtherFingerprintBitVector, $Alpha);
    $Value = Fingerprints::FingerprintsBitVector::
              TverskySimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB, $Alpha);

Returns value of I<Tversky> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<WeightedTanimotoSimilarityCoefficient>

    $Value =
       $FingerprintsBitVector->WeightedTanimotoSimilarityCoefficient(
         $OtherFingerprintBitVector, $Beta);
    $Value =
       Fingerprints::FingerprintsBitVector::
         WeightedTanimotoSimilarityCoefficient(
         $FingerprintsBitVectorA, $FingerprintBitVectorB, $Beta);

Returns value of I<WeightedTanimoto> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<WeightedTverskySimilarityCoefficient>

    $Value =
       $FingerprintsBitVector->WeightedTverskySimilarityCoefficient(
          $OtherFingerprintBitVector, $Alpha, $Beta);
    $Value =
      Fingerprints::FingerprintsBitVector::
        WeightedTverskySimilarityCoefficient(
        $FingerprintsBitVectorA, $FingerprintBitVectorB, $Alpha, $Beta);

Returns value of I<WeightedTversky> similarity coefficient between two same size I<FingerprintsBitVectors>.

=item B<YuleSimilarityCoefficient>

    $Value = $FingerprintsBitVector->YuleSimilarityCoefficient(
              $OtherFingerprintBitVector);
    $Value = Fingerprints::FingerprintsBitVector::YuleSimilarityCoefficient(
              $FingerprintsBitVectorA, $FingerprintBitVectorB);

Returns value of I<Yule> similarity coefficient between two same size I<FingerprintsBitVectors>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

BitVector.pm, FingerprintsStringUtil.pm, FingerprintsVector.pm, Vector.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
