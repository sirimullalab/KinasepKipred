package BitVector;
#
# File: BitVector.pm
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
use ConversionsUtil ();
use MathUtil;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw(IsBitVector);
@EXPORT_OK = qw(NewFromBinaryString NewFromDecimalString NewFromHexadecimalString NewFromOctalString NewFromRawBinaryString);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, $ValueFormat, $ValueBitOrder);
_InitializeClass();

#
# Overload bitwise and some logical operators...
#
# 'fallback' is set to 'false' to raise exception for all other operators.
#
use overload '""' => 'StringifyBitVector',
  '&' => '_BitVectorAndOperator',
  '|' => '_BitVectorOrOperator',
  '^' => '_BitVectorExclusiveOrOperator',

  '~' => '_BitVectorNegationOperator',

  '==' => '_BitVectorEqualOperator',
  '!=' => '_BitVectorNotEqualOperator',

  'fallback' => undef;

# Class constructor...
#
sub new {
  my($Class, $Size) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializeBitVector($Size);

  return $This;
}

# Initialize object data...
#
# Note:
#  . Perl pack function used to initialize vector automatically sets its size to
#    nearest power of 2 value.
#
sub _InitializeBitVector {
  my($This, $Size) = @_;

  if (!defined $Size) {
    croak "Error: ${ClassName}->new: BitVector object instantiated without specifying its size ...";
  }
  if ($Size <=0) {
    croak "Error: ${ClassName}->new: Bit vector size, $Size, must be a positive integer...";
  }

  # Initialize vector with zeros...
  $This->{BitValues} = pack("b*", "0" x $Size);

  # Size to automatically set to nearest power of 2 by Perl pack function. So use the length
  # of packed vector to set size...
  $This->{Size} = length($This->GetBitsAsBinaryString());

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Print format for bit vectore values...
  $ValueFormat = "Binary";

  # Bit ordering for printing bit vector value strings. Default is to print lowest bit of each
  # byte on the left.
  #
  # Internally, bits are stored in ascending order using Perl vec function. Regardless
  # of machine order, big-endian or little-endian, vec function always considers first
  # string byte as the lowest byte and first bit within each byte as the lowest bit.
  #
  # Possible values: Ascending or Descending
  #
  $ValueBitOrder = 'Ascending';
}

# Create a new bit vector using binary string. This functionality can be
# either invoked as a class function or an object method.
#
# The size of bit vector is automatically set to reflect the string.
#
sub NewFromBinaryString ($;$) {
  my($FirstParameter, $SecondParameter, $ThirdParameter) = @_;

  if (_IsBitVector($FirstParameter)) {
    return _NewBitVectorFromString('Binary', $SecondParameter, $ThirdParameter);
  }
  else {
    return _NewBitVectorFromString( 'Binary', $FirstParameter, $SecondParameter);
  }
}

# Create a new bit vector using hexadecimal string. This functionality can be
# either invoked as a class function or an object method.
#
# The size of bit vector is automatically set to reflect the string.
#
sub NewFromHexadecimalString ($;$) {
  my($FirstParameter, $SecondParameter, $ThirdParameter) = @_;

  if (_IsBitVector($FirstParameter)) {
    return _NewBitVectorFromString('Hexadecimal', $SecondParameter, $ThirdParameter);
  }
  else {
    return _NewBitVectorFromString( 'Hexadecimal', $FirstParameter, $SecondParameter);
  }
}

# Create a new bit vector using octal string. This functionality can be
# either invoked as a class function or an object method.
#
# The size of bit vector is automatically set to reflect the string.
#
sub NewFromOctalString ($;$) {
  my($FirstParameter, $SecondParameter, $ThirdParameter) = @_;

  if (_IsBitVector($FirstParameter)) {
    return _NewBitVectorFromString('Octal', $SecondParameter, $ThirdParameter);
  }
  else {
    return _NewBitVectorFromString( 'Octal', $FirstParameter, $SecondParameter);
  }
}

# Create a new bit vector using decimal string. This functionality can be
# either invoked as a class function or an object method.
#
# The size of bit vector is automatically set to reflect the string.
#
sub NewFromDecimalString ($;$) {
  my($FirstParameter, $SecondParameter, $ThirdParameter) = @_;

  if (_IsBitVector($FirstParameter)) {
    return _NewBitVectorFromString('Decimal', $SecondParameter, $ThirdParameter);
  }
  else {
    return _NewBitVectorFromString( 'Decimal', $FirstParameter, $SecondParameter);
  }
}

# Create a new bit vector using raw binary string. This functionality can be
# either invoked as a class function or an object method.
#
# The size of bit vector is automatically set to reflect the string.
#
sub NewFromRawBinaryString ($;$) {
  my($FirstParameter, $SecondParameter, $ThirdParameter) = @_;

  if (_IsBitVector($FirstParameter)) {
    return _NewBitVectorFromString('RawBinary', $SecondParameter, $ThirdParameter);
  }
  else {
    return _NewBitVectorFromString( 'RawBinary', $FirstParameter, $SecondParameter);
  }
}

# Create a new bit vector from a string...
#
sub _NewBitVectorFromString ($$;$) {
  my($Format, $String, $BitOrder) = @_;
  my($Size, $BitVector);

  $Size = _CalculateStringSizeInBits($Format, $String);

  $BitVector = new BitVector($Size);
  $BitVector->_SetBitsAsString($Format, $String, $BitOrder);

  return $BitVector;
}

# Copy bit vector...
sub Copy {
  my($This) = @_;
  my($BitVector);

  # Make a new bit vector...
  $BitVector = (ref $This)->new($This->{Size});

  # Copy bit values...
  $BitVector->{BitValues} = $This->{BitValues};

  # Copy value format for stringification...
  if (exists $This->{ValueFormat}) {
    $BitVector->{ValueFormat} = $This->{ValueFormat};
  }
  # Copy value bit order for stringification...
  if (exists $This->{ValueBitOrder}) {
    $BitVector->{ValueBitOrder} = $This->{ValueBitOrder};
  }
  return $BitVector;
}

# Reverse bit values in bit vector...
sub Reverse {
  my($This) = @_;
  my($BitNum, $ReverseBitNum, $BitValue, $ReverseBitValue);

  $BitNum = 0; $ReverseBitNum = $This->{Size} - 1;

  while ($BitNum < $ReverseBitNum) {
    $BitValue = $This->_GetBitValue($BitNum);
    $ReverseBitValue = $This->_GetBitValue($ReverseBitNum);

    $This->_SetBitValue($BitNum, $ReverseBitValue);
    $This->_SetBitValue($ReverseBitNum, $BitValue);

    $BitNum++; $ReverseBitNum--;
  }
  return $This;
}

# Is it a bit vector object?
sub IsBitVector ($) {
  my($Object) = @_;

  return _IsBitVector($Object);
}

# Get size...
sub GetSize {
  my($This) = @_;

  return $This->{Size};
}

# Set a bit...
#
sub SetBit {
  my($This, $BitNum, $SkipCheck) = @_;

  # Just set it...
  if ($SkipCheck) {
    return $This->_SetBitValue($BitNum, 1);
  }

  # Check and set...
  $This->_ValidateBitNumber("SetBit", $BitNum);

  return $This->_SetBitValue($BitNum, 1);
}

# Set arbitrary bits specified as a list of bit numbers...
#
sub SetBits {
  my($This, @BitNums) = @_;
  my($BitNum);

  for $BitNum (@BitNums) {
    $This->SetBit($BitNum);
  }
  return $This;
}

# Set bits in a specified range...
#
sub SetBitsRange {
  my($This, $MinBitNum, $MaxBitNum) = @_;
  my($BitNum);

  $This->_ValidateBitNumber("SetBitsRange", $MinBitNum);
  $This->_ValidateBitNumber("SetBitsRange", $MaxBitNum);

  for $BitNum ($MinBitNum .. $MaxBitNum) {
    $This->_SetBitValue($BitNum, 1);
  }
  return $This;
}

# Set all bits...
#
sub SetAllBits {
  my($This) = @_;

  $This->{BitValues} = pack("b*", "1" x $This->{Size});
}

# Clear a bit...
#
sub ClearBit {
  my($This, $BitNum) = @_;

  $This->_ValidateBitNumber("ClearBit", $BitNum);

  return $This->_SetBitValue($BitNum, 0);
}

# Clear arbitrary bits specified as a list of bit numbers...
#
sub ClearBits {
  my($This, @BitNums) = @_;
  my($BitNum);

  for $BitNum (@BitNums) {
    $This->ClearBit($BitNum);
  }
  return $This;
}

# Clear bits in a specified range...
#
sub ClearBitsRange {
  my($This, $MinBitNum, $MaxBitNum) = @_;
  my($BitNum);

  $This->_ValidateBitNumber("ClearBitsRange", $MinBitNum);
  $This->_ValidateBitNumber("ClearBitsRange", $MaxBitNum);

  for $BitNum ($MinBitNum .. $MaxBitNum) {
    $This->_SetBitValue($BitNum, 0);
  }
  return $This;
}

# Clear all bits...
#
sub ClearAllBits {
  my($This) = @_;

  $This->{BitValues} = pack("b*", "0" x $This->{Size});

  return $This;
}

# Set or clear bit...
#
sub SetBitValue {
  my($This, $BitNum, $BitValue) = @_;

 BITVALUE: {
    if ($BitValue == 1) { return $This->SetBit($BitNum); last BITVALUE; }
    if ($BitValue == 0) { return $This->ClearBit($BitNum); last BITVALUE; }
    croak "Error: ${ClassName}->SetBit: Specified bit value, $BitValue, must be 0 or 1...";
  }
  return $This;
}

# Flip bit value...
#
sub FlipBit {
  my($This, $BitNum) = @_;

  $This->_ValidateBitNumber("FlipBit", $BitNum);
  return $This->_FlipBit($BitNum);
}

# Flip arbitrary bits specified as a list of bit numbers...
#
sub FlipBits {
  my($This, @BitNums) = @_;
  my($BitNum);

  for $BitNum (@BitNums) {
    $This->FlipBit();
  }
  return $This;
}

# Flip bit value in a specified bit range...
#
sub FlipBitsRange {
  my($This, $MinBitNum, $MaxBitNum) = @_;
  my($BitNum);

  $This->_ValidateBitNumber("FlipBitsRange", $MinBitNum);
  $This->_ValidateBitNumber("FlipBitsRange", $MaxBitNum);

  for $BitNum ($MinBitNum .. $MaxBitNum) {
    $This->_FlipBit();
  }
  return $This;
}

# Flip all bit valus...
#
sub FlipAllBits {
  my($This) = @_;

  return $This->FlipBits(0, ($This->{Size} - 1));
}

# Flip bit value...
sub _FlipBit {
  my($This, $BitNum) = @_;

  if ($This->_GetBitValue($BitNum)) {
    return $This->_SetBitValue($BitNum, 0);
  }
  else {
    return $This->_SetBitValue($BitNum, 1);
  }
}

# Get bit value...
#
sub GetBit {
  my($This, $BitNum) = @_;

  $This->_ValidateBitNumber("GetBit", $BitNum);

  return $This->_GetBitValue($BitNum);
}

# Is a specific bit set?
#
sub IsBitSet {
  my($This, $BitNum) = @_;

  if (!(defined($BitNum) && ($BitNum >= 0) && ($BitNum < $This->{Size}))) {
    return undef;
  }

  return $This->_GetBitValue($BitNum) ? 1 : 0;
}

# Is a specific bit clear?
#
sub IsBitClear {
  my($This, $BitNum) = @_;

  if (!(defined($BitNum) && ($BitNum >= 0) && ($BitNum < $This->{Size}))) {
    return undef;
  }

  return $This->_GetBitValue($BitNum) ? 0 : 1;
}

# Get number of set bits...
#
sub GetNumOfSetBits {
  my($This) = @_;

  return unpack("%b*", $This->{BitValues});
}

# Get number of clear bits...
#
sub GetNumOfClearBits {
  my($This) = @_;

  return ($This->{Size} - $This->GetNumOfSetBits());
}

# Get density of set bits...
#
sub GetDensityOfSetBits {
  my($This) = @_;

  return $This->{Size} ? ($This->GetNumOfSetBits()/$This->{Size}) : 0;
}

# Get density of clear bits...
#
sub GetDensityOfClearBits {
  my($This) = @_;

  return $This->GetNumOfClearBits()/$This->{Size};
}

# Convert internal bit values stored using Perl vec function with first string byte
# as the lowest byte and first bit within each byte as the lowest bit into a binary
# string with ascending or descending bit order within each byte. The internal
# bit order corresponds to ascending bit order within each byte.
#
sub GetBitsAsBinaryString {
  my($This, $BitOrder) = @_;

  return $This->_GetBitsAsString('Binary', $BitOrder);
}

# Convert internal bit values stored using Perl vec function with first string byte
# as the lowest byte and first bit within each byte as the lowest bit into a hexadecimal
# string with ascending or descending bit order within each byte. The internal
# bit order corresponds to ascending bit order within each byte.
#
#
sub GetBitsAsHexadecimalString {
  my($This, $BitOrder) = @_;

  return $This->_GetBitsAsString('Hexadecimal', $BitOrder);
}

# Convert bit values into a octal string value...
#
sub GetBitsAsOctalString {
  my($This, $BitOrder) = @_;

  return $This->_GetBitsAsString('Octal', $BitOrder);
}

# Convert bit values into a decimal string value...
#
sub GetBitsAsDecimalString {
  my($This, $BitOrder) = @_;

  return $This->_GetBitsAsString('Decimal', $BitOrder);
}

# Return packed bit values which also contains nonprintable characters...
#
sub GetBitsAsRawBinaryString {
  my($This) = @_;

  return $This->_GetBitsAsString('RawBinary');
}

# Convert internal bit values stored using Perl vec function with first string byte
# as the lowest byte and first bit within each byte as the lowest bit into a
# string with ascending or descending bit order within each byte. The internal
# bit order corresponds to ascending bit order within each byte.
#
#
sub _GetBitsAsString {
  my($This, $Format, $BitOrder) = @_;
  my($BinaryTemplate, $HexadecimalTemplate);

  ($BinaryTemplate, $HexadecimalTemplate) = $This->_SetupBitsPackUnpackTemplate($BitOrder);

  FORMAT : {
    if ($Format =~ /^(Hexadecimal|Hex|HexadecimalString)$/i) { return unpack($HexadecimalTemplate, $This->{BitValues}); last FORMAT; }
    if ($Format =~ /^(Octal|Oct|OctalString)$/i) { return ConversionsUtil::HexadecimalToOctal(unpack($HexadecimalTemplate, $This->{BitValues})); last FORMAT; }
    if ($Format =~ /^(Decimal|Dec|DecimalString)$/i) { return ConversionsUtil::HexadecimalToDecimal(unpack($HexadecimalTemplate, $This->{BitValues})); last FORMAT; }
    if ($Format =~ /^(Binary|Bin|BinaryString)$/i) { return unpack($BinaryTemplate, $This->{BitValues}); last FORMAT; }
    if ($Format =~ /^(RawBinary|RawBin|RawBinaryString)$/i) { return $This->{BitValues}; last FORMAT; }
    croak "Error: ${ClassName}->_GetBitsAsString: Specified bit vector string format, $Format, is not supported. Value values: Binary, Bin, BinaryString, Hexdecimal, Hex, HexadecimalString, Decimal, Dec, DecimalString, Octal, Oct, OctalString, RawBinary, RawBin, RawBinaryString...";
  }
}

# Setup templates to unpack bits...
#
sub _SetupBitsPackUnpackTemplate {
  my($This, $BitOrder) = @_;
  my($BinaryTemplate, $HexadecimalTemplate);

  $BitOrder = (defined($BitOrder) && $BitOrder) ? $BitOrder : 'Ascending';

  if ($BitOrder =~ /^Ascending$/i) {
    $BinaryTemplate = "b*";
    $HexadecimalTemplate = "h*";
  }
  elsif ($BitOrder =~ /^Descending$/i) {
    $BinaryTemplate = "B*";
    $HexadecimalTemplate = "H*";
  }
  else {
    croak "Warning: ${ClassName}::_SetupBitsPackUnpackTemplate: Specified bit order value, $BitOrder, is not supported. Supported values: Ascending, Descending...";
  }
  return ($BinaryTemplate, $HexadecimalTemplate);
}

# Set bit values using hexadecimal string. The initial size of bit vector is not changed.
#
sub SetBitsAsHexadecimalString {
  my($This, $Hexadecimal, $BitOrder) = @_;

  if ($Hexadecimal =~ /^0x/i) {
    $Hexadecimal =~ s/^0x//i;
  }
  return $This->_SetBitsAsString('Hexadecimal', $Hexadecimal, $BitOrder);
}

# Set bit values using octal string. The initial size of bit vector is not changed.
#
sub SetBitsAsOctalString {
  my($This, $Octal, $BitOrder) = @_;

  if ($Octal =~ /^0/i) {
    $Octal =~ s/^0//i;
  }
  return $This->_SetBitsAsString('Octal', $Octal, $BitOrder);
}

# Set bit values using a decimal number. The initial size of bit vector is not changed.
#
sub SetBitsAsDecimalString {
  my($This, $Decimal, $BitOrder) = @_;

  if (!TextUtil::IsPositiveInteger($Decimal)) {
    croak "Error: ${ClassName}->SetBitsAsDecimalString: Specified decimal value, $Decimal, must be a positive integer...";
  }
  if ($Decimal =~ /[+]/) {
    $Decimal =~ s/[+]//;
  }
  return $This->_SetBitsAsString('Decimal', $Decimal, $BitOrder);
}

# Set bit values using hexadecimal string. The initial size of bit vector is not changed.
#
sub SetBitsAsBinaryString {
  my($This, $Binary, $BitOrder) = @_;

  if ($Binary =~ /^0b/i) {
    $Binary =~ s/^0b//i;
  }
  return $This->_SetBitsAsString('Binary', $Binary, $BitOrder);
}

# Set bit values using packed binary string. The size of bit vector is changed to reflect
# the input raw string...
#
sub SetBitsAsRawBinaryString {
  my($This, $RawBinary) = @_;

  return $This->_SetBitsAsString('RawBinary', $RawBinary);
}

# Set bits using string in a specified format. This size of bit vector is not changed except for
# RawBinary string type...
#
sub _SetBitsAsString {
  my($This, $Format, $String, $BitOrder) = @_;
  my($Size, $BinaryTemplate, $HexadecimalTemplate);

  ($BinaryTemplate, $HexadecimalTemplate) = $This->_SetupBitsPackUnpackTemplate($BitOrder);

  $Size = $This->{Size};
  FORMAT : {
    if ($Format =~ /^(Hexadecimal|Hex|HexadecimalString)$/i) { $This->{BitValues} = pack($HexadecimalTemplate, $String); last FORMAT; }
    if ($Format =~ /^(Octal|Oct|OctalString)$/i) { vec($This->{BitValues}, 0, $Size) = ConversionsUtil::OctalToDecimal($String); last FORMAT; }
    if ($Format =~ /^(Decimal|Dec|DecimalString)$/i) { vec($This->{BitValues}, 0, $Size) = $String; last FORMAT; }
    if ($Format =~ /^(Binary|Bin|BinaryString)$/i) { $This->{BitValues} = pack($BinaryTemplate, $String); last FORMAT; }
    if ($Format =~ /^(RawBinary|RawBin|RawBinaryString)$/i) { $This->{BitValues} = $String; last FORMAT; }
    croak "Error: ${ClassName}->_SetBitsAsString: Specified bit vector string format, $Format, is not supported. Value values: Binary, Bin, BinaryString, Hexdecimal, Hex, HexadecimalString, Decimal, Dec, DecimalString, Octal, Oct, OctalString, RawBinary, RawBin, RawBinaryString...";
  }

  # Set size using packed string...
  $Size = length($This->GetBitsAsBinaryString());
  if ($Size <=0) {
    croak "Error: ${ClassName}->_SetBitsAsString: Bit vector size, $Size, must be a positive integer...";
  }
  $This->{Size} = $Size;

  return $This;
}

# Calculate string size in bits...
#
sub _CalculateStringSizeInBits ($$;$) {
  my($FirstParameter, $SecondParameter, $ThisParameter) = @_;
  my($This, $Format, $String, $Size);

  if ((@_ == 3) && (_IsBitVector($FirstParameter))) {
    ($This, $Format, $String) = ($FirstParameter, $SecondParameter, $ThisParameter);
  }
  else {
    ($This, $Format, $String) = (undef, $FirstParameter, $SecondParameter);
  }

  FORMAT : {
    if ($Format =~ /^(Hexadecimal|Hex|HexadecimalString)$/i) { $Size = length($String) * 4; last FORMAT; }
    if ($Format =~ /^(Octal|Oct|OctalString)$/i) { $Size = length($String) * 3; last FORMAT; }
    if ($Format =~ /^(Decimal|Dec|DecimalString)$/i) { $Size = length(ConversionsUtil::DecimalToHexadecimal($String)) * 4; last FORMAT; }
    if ($Format =~ /^(Binary|Bin|BinaryString)$/i) { $Size = length($String); last FORMAT; }
    if ($Format =~ /^(RawBinary|RawBin|RawBinaryString)$/i) { $Size = length(unpack("B*", $String)); last FORMAT; }
    croak "Error: ${ClassName}::_CalculateStringSizeInBits: Specified bit vector string format, $Format, is not supported. Value values: Binary, Bin, BinaryString, Hexdecimal, Hex, HexadecimalString, Decimal, Dec, DecimalString, Octal, Oct, OctalString, RawBinary, RawBin, RawBinaryString...";
  }
  return $Size;
}

# Set bit value using Perl vec function with bit numbers going from left to right.
# First bit number corresponds to 0.
#
sub _SetBitValue {
  my($This, $BitNum, $BitValue) = @_;
  my($Offset, $Width);

  $Offset = $BitNum;
  $Width = 1;

  vec($This->{BitValues}, $Offset, $Width) = $BitValue;

  return $This;
}

# Get bit value Perl vec function with bit numbers going from left to right.
# First bit number corresponds to 0.
#
sub _GetBitValue {
  my($This, $BitNum) = @_;
  my($Offset, $Width, $BitValue);

  $Offset = $BitNum;
  $Width = 1;

  $BitValue = vec($This->{BitValues}, $Offset, $Width);

  return $BitValue;
}

# Check to make sure it's a valid bit number...
#
sub _ValidateBitNumber {
  my($This, $CallerName, $BitNum) = @_;

  if (!defined $BitNum) {
    croak "Error: ${ClassName}->${CallerName}: Bit number is not defined...";
  }
  if ($BitNum < 0) {
    croak "Error: ${ClassName}->${CallerName}: Bit number value, $BitNum, must be >= 0 ...";
  }
  if ($BitNum >= $This->{Size}) {
    croak "Error: ${ClassName}->${CallerName}: Bit number number value, $BitNum, must be less than the size of bit vector, ", $This->{Size}, "...";
  }

  return $This;
}

# Set bit values print format for an individual object or the whole class...
#
sub SetBitValuePrintFormat ($;$) {
  my($FirstParameter, $SecondParameter) = @_;

  if ((@_ == 2) && (_IsBitVector($FirstParameter))) {
    # Set bit values print format for the specific object...
    my($This, $ValuePrintFormat) = ($FirstParameter, $SecondParameter);

    if (!_ValidateBitValuePrintFormat($ValuePrintFormat)) {
      return;
    }

    $This->{ValueFormat} = $ValuePrintFormat;
  }
  else {
    # Set value print format for the class...
    my($ValuePrintFormat) = ($FirstParameter);

    if (!_ValidateBitValuePrintFormat($ValuePrintFormat)) {
      return;
    }

    $ValueFormat = $ValuePrintFormat;
  }
}

# Set bit values bit order for an individual object or the whole class...
#
sub SetBitValueBitOrder ($;$) {
  my($FirstParameter, $SecondParameter) = @_;

  if ((@_ == 2) && (_IsBitVector($FirstParameter))) {
    # Set bit value bit order for the specific object...
    my($This, $BitOrder) = ($FirstParameter, $SecondParameter);

    if (!_ValidateBitValueBitOrder($BitOrder)) {
      return;
    }

    $This->{ValueBitOrder} = $BitOrder;
  }
  else {
    # Set bit value bit order for the class...
    my($BitOrder) = ($FirstParameter);

    if (!_ValidateBitValueBitOrder($BitOrder)) {
      return;
    }

    $ValueBitOrder = $BitOrder;
  }
}

# Validate print format for bit values...
sub _ValidateBitValueBitOrder {
  my($BitOrder) = @_;

  if ($BitOrder !~ /^(Ascending|Descending)$/i) {
    carp "Warning: ${ClassName}::_ValidateBitValueBitOrder: Specified bit order value, $BitOrder, is not supported. Supported values: Ascending, Descending...";
    return 0;
  }
  return 1;
}

# Validate print format for bit values...
sub _ValidateBitValuePrintFormat {
  my($ValuePrintFormat) = @_;

  if ($ValuePrintFormat !~ /^(Binary|Bin||BinaryString|Hexadecimal|Hex||HexadecimalString|Decimal|Dec||DecimalString|Octal|Oct||OctalString|RawBinary|RawBin|RawBinaryString)$/i) {
    carp "Warning: ${ClassName}::_ValidateBitValuePrintFormat: Specified bit vector print format value, $ValuePrintFormat, is not supported. Supported values: Binary, Bin, BinaryString, Hexdecimal, Hex, HexadecimalString, Decimal, Dec, DecimalString, Octal, Oct, OctalString, RawBinary, RawBin, RawBinaryString...";
    return 0;
  }
  return 1;
}

# Bitwise AND operation for BitVectors...
#
sub _BitVectorAndOperator {
  my($This, $Other, $OrderFlipped, $OtherIsBitVector, $ErrorMsg, $CheckBitVectorSizes);

  $ErrorMsg = "_BitVectorAndOperator: Bitwise AND oparation failed";
  $CheckBitVectorSizes = 1;
  ($This, $Other, $OrderFlipped, $OtherIsBitVector) = _ProcessOverloadedOperatorParameters($ErrorMsg, @_, $CheckBitVectorSizes);

  if (!$OtherIsBitVector) {
    if ($OrderFlipped) {
      croak "Error: ${ClassName}->${ErrorMsg}: First object must be a bit vector...";
    }
  }
  my($BitVector);
  $BitVector = (ref $This)->new($This->{Size});
  $BitVector->{BitValues} = $This->{BitValues} & $Other->{BitValues};

  return $BitVector;
}

# Bitwise OR operation for BitVectors...
#
sub _BitVectorOrOperator {
  my($This, $Other, $OrderFlipped, $OtherIsBitVector, $ErrorMsg, $CheckBitVectorSizes);

  $ErrorMsg = "_BitVectorAndOperator: Bitwise OR oparation failed";
  $CheckBitVectorSizes = 1;
  ($This, $Other, $OrderFlipped, $OtherIsBitVector) = _ProcessOverloadedOperatorParameters($ErrorMsg, @_, $CheckBitVectorSizes);

  if (!$OtherIsBitVector) {
    if ($OrderFlipped) {
      croak "Error: ${ClassName}->${ErrorMsg}: First object must be a bit vector...";
    }
  }
  my($BitVector);
  $BitVector = (ref $This)->new($This->{Size});
  $BitVector->{BitValues} = $This->{BitValues} | $Other->{BitValues};

  return $BitVector;
}

# Bitwise XOR operation for BitVectors...
#
sub _BitVectorExclusiveOrOperator {
  my($This, $Other, $OrderFlipped, $OtherIsBitVector, $ErrorMsg, $CheckBitVectorSizes);

  $ErrorMsg = "_BitVectorAndOperator: Bitwise XOR oparation failed";
  $CheckBitVectorSizes = 1;
  ($This, $Other, $OrderFlipped, $OtherIsBitVector) = _ProcessOverloadedOperatorParameters($ErrorMsg, @_, $CheckBitVectorSizes);

  if (!$OtherIsBitVector) {
    if ($OrderFlipped) {
      croak "Error: ${ClassName}->${ErrorMsg}: First object must be a bit vector...";
    }
  }
  my($BitVector);
  $BitVector = (ref $This)->new($This->{Size});
  $BitVector->{BitValues} = $This->{BitValues} ^ $Other->{BitValues};

  return $BitVector;
}

# Bitwise negation operation for BitVectors...
#
sub _BitVectorNegationOperator {
  my($This, $Other, $OrderFlipped, $OtherIsBitVector, $ErrorMsg, $CheckBitVectorSizes);

  $ErrorMsg = "_BitVectorAndOperator: Bitwise negation oparation failed";
  $CheckBitVectorSizes = 1;
  ($This, $Other, $OrderFlipped, $OtherIsBitVector) = _ProcessOverloadedOperatorParameters($ErrorMsg, @_, $CheckBitVectorSizes);

  my($BitVector);
  $BitVector = (ref $This)->new($This->{Size});
  $BitVector->{BitValues} = ~ $This->{BitValues};

  return $BitVector;
}

# Bit vector equla operator. Two bit vectors are considered equal assuming their size
# is same and bits are on at the same positions...
#
sub _BitVectorEqualOperator {
  my($This, $Other, $OrderFlipped, $OtherIsBitVector, $ErrorMsg, $CheckBitVectorSizes);

  $ErrorMsg = "_BitVectorEqualOperator: BitVector == oparation failed";
  $CheckBitVectorSizes = 0;
  ($This, $Other, $OrderFlipped, $OtherIsBitVector) = _ProcessOverloadedOperatorParameters($ErrorMsg, @_, $CheckBitVectorSizes);

  if (!$OtherIsBitVector) {
    if ($OrderFlipped) {
      croak "Error: ${ClassName}->${ErrorMsg}: First object must be a bit vector...";
    }
  }
  if ($This->GetSize() != $Other->GetSize()) {
    return 0;
  }
  if ($This->GetNumOfSetBits() != $Other->GetNumOfSetBits()) {
    return 0;
  }
  # Check number of On bits only in This vector. It must be zero for vectors to be equal...
  my($BitVector);
  $BitVector = $This & ~$Other;

  return $BitVector->GetNumOfSetBits() ? 0 : 1;
}

# Bit vector not equal operator. Two bit vectors are considered not equal when their size
# is different or bits are on at the same positions...
#
sub _BitVectorNotEqualOperator {
  my($This, $Other, $OrderFlipped, $OtherIsBitVector, $ErrorMsg, $CheckBitVectorSizes);

  $ErrorMsg = "_BitVectorEqualOperator: BitVector != oparation failed";
  $CheckBitVectorSizes = 0;
  ($This, $Other, $OrderFlipped, $OtherIsBitVector) = _ProcessOverloadedOperatorParameters($ErrorMsg, @_, $CheckBitVectorSizes);

  if (!$OtherIsBitVector) {
    if ($OrderFlipped) {
      croak "Error: ${ClassName}->${ErrorMsg}: First object must be a bit vector...";
    }
  }
  if ($This->GetSize() != $Other->GetSize()) {
    return 1;
  }
  if ($This->GetNumOfSetBits() != $Other->GetNumOfSetBits()) {
    return 1;
  }
  # Check number of On bits only in This vector. It must be zero for vectors to be equal...
  my($BitVector);
  $BitVector = $This & ~$Other;

  return $BitVector->GetNumOfSetBits() ? 1 : 0;
}

# Process parameters passed to overloaded operators...
#
# For uninary operators, $SecondParameter is not defined.
sub _ProcessOverloadedOperatorParameters {
  my($ErrorMsg, $FirstParameter, $SecondParameter, $ParametersOrderStatus, $CheckBitVectorSizesStatus) = @_;
  my($This, $Other, $OrderFlipped, $OtherIsBitVector, $CheckBitVectorSizes);

  ($This, $Other) =  ($FirstParameter, $SecondParameter);
  $OrderFlipped = (defined($ParametersOrderStatus) && $ParametersOrderStatus) ? 1 : 0;
  $CheckBitVectorSizes = (defined $CheckBitVectorSizesStatus) ? $CheckBitVectorSizesStatus : 1;

  _ValidateBitVector($ErrorMsg, $This);

  $OtherIsBitVector = 0;
  if (defined($Other) && (ref $Other)) {
    # Make sure $Other is a vector...
    _ValidateBitVector($ErrorMsg, $Other);
    if ($CheckBitVectorSizes) {
      _ValidateBitVectorSizesAreEqual($ErrorMsg, $This, $Other);
    }
    $OtherIsBitVector = 1;
  }
  return ($This, $Other, $OrderFlipped, $OtherIsBitVector);
}

# Is it a bit vector object?
sub _IsBitVector {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Make sure it's a bit vector reference...
sub _ValidateBitVector {
  my($ErrorMsg, $Vector) = @_;

  if (!_IsBitVector($Vector)) {
    croak "Error: ${ClassName}->${ErrorMsg}: Object must be a bit vector...";
  }
}

# Make sure size of the two bit vectors are equal...
sub _ValidateBitVectorSizesAreEqual {
  my($ErrorMsg, $BitVector1, $BitVector2) = @_;

  if ($BitVector1->GetSize() != $BitVector2->GetSize()) {
    croak "Error: ${ClassName}->${ErrorMsg}: Size of the bit vectors must be same...";
  }
}

# Return a string containing vector values...
sub StringifyBitVector {
  my($This) = @_;
  my($BitVectorString, $PrintFormat, $BitOrder, $BitsValue);

  $PrintFormat = (exists $This->{ValueFormat}) ? $This->{ValueFormat} : $ValueFormat;
  $BitOrder = (exists $This->{ValueBitOrder}) ? $This->{ValueBitOrder} : $ValueBitOrder;
  $BitVectorString = '';

  FORMAT: {
      if ($PrintFormat =~ /^(Hexadecimal|Hex|HexadecimalString)$/i) { $BitsValue = $This->_GetBitsAsString('Hexadecimal', $BitOrder);  last FORMAT; }
      if ($PrintFormat =~ /^(Octal|Oct|OctalString)$/i) { $BitsValue = $This->_GetBitsAsString('Octal', $BitOrder);  last FORMAT; }
      if ($PrintFormat =~ /^(Decimal|Dec|DecimalString)$/i) { $BitsValue = $This->_GetBitsAsString('Decimal', $BitOrder);  last FORMAT; }
      if ($PrintFormat =~ /^(RawBinary|RawBin|RawBinaryString)$/i) { $BitsValue = $This->_GetBitsAsString('RawBinary');  last FORMAT; }
      # Default is bninary format...
      $BitsValue = $This->_GetBitsAsString('Binary', $BitOrder);
  }
  $BitVectorString = "<Size: ". $This->GetSize() . ";BitOrder: $BitOrder; Value: " . $BitsValue . ">";

  return $BitVectorString;
}

1;

__END__

=head1 NAME

BitVector

=head1 SYNOPSIS

use BitVector;

use BitVector ();

use BitVector qw(:all);

=head1 DESCRIPTION

B<BitVector> class provides the following methods:

new, ClearAllBits, ClearBit, ClearBits, ClearBitsRange, Copy, FlipAllBits,
FlipBit, FlipBits, FlipBitsRange, GetBit, GetBitsAsBinaryString,
GetBitsAsDecimalString, GetBitsAsHexadecimalString, GetBitsAsOctalString,
GetBitsAsRawBinaryString, GetDensityOfClearBits, GetDensityOfSetBits,
GetNumOfClearBits, GetNumOfSetBits, GetSize, IsBitClear, IsBitSet, IsBitVector,
NewFromBinaryString, NewFromDecimalString, NewFromHexadecimalString,
NewFromOctalString, NewFromRawBinaryString, Reverse, SetAllBits, SetBit,
SetBitValue, SetBitValueBitOrder, SetBitValuePrintFormat, SetBits,
SetBitsAsBinaryString, SetBitsAsDecimalString, SetBitsAsHexadecimalString,
SetBitsAsOctalString, SetBitsAsRawBinaryString, SetBitsRange, StringifyBitVector

The following methods can also be used as functions:

IsBitVector, NewFromBinaryString, NewFromDecimalString, NewFromHexadecimalString,
NewFromOctalString, NewFromRawBinaryString

The following operators are overloaded:

    "" & | ^ ~ == !=

Internally, bits are stored in ascending order using Perl vec function. Regardless
of machine order, big-endian or little-endian, vec function always considers first
string byte as the lowest byte and first bit within each byte as the lowest bit.

Things to keep in mind:

    o Bit numbers range from 0 to (Size - 1).
    o Bit data retieval methods provide options to data in ascending or
      descending bit order. Default is ascending bit order.
    o Stringyfy method provides an option to print data in ascending or
      descending bit order. Default is ascending bit order.

=head2 METHODS

=over 4

=item B<new>

    $NewBitVector = new BitVector($Size);

Create a new I<BitVector> object of size I<Size> and return  newly created
B<BitVector>. Bit numbers range from 0 to 1 less than I<Size>.

=item B<ClearAllBits>

    $BitVector->ClearAllBits();

Set all bit values to 0 in I<BitVector> object and return I<BitVector>.

=item B<ClearBit>

    $BitVector->ClearBit($BitNum);

Set specified bit number I<BitNum> to 0 in I<BitVector> object and return I<BitVector>.

=item B<ClearBits>

    $BitVector->ClearBits(@BitNums);

Set specified bit numbers I<BitNums> to 0 in I<BitVector> object and return I<BitVector>.

=item B<ClearBitsRange>

    $BitVector->ClearBitsRange($MinBitNum, $MaxBitNum);

Set specified bit numbers between I<MinBitNum> and I<MaxBitNum> to 0 in I<BitVector>
object and return I<BitVector>.

=item B<Copy>

    $NewBitVector = $BitVector->Copy();

Copy I<BitVector> and its associated data to a new B<BitVector> and return a new
B<BitVector>.

=item B<FlipAllBits>

    $BitVector->FlipAllBits();

Flip values of all bits in I<BitVector> and its associated data to a new B<BitVector> and return
I<BitVector>.

=item B<FlipBit>

    $BitVector->FlipBit($BitNum);

Flip value of specified I<BitNum> of in I<BitVector> and return I<BitVector>.

=item B<FlipBits>

    $BitVector->FlipBits(@BitNums);

Flip values of specified bit numbers I<BitNums> in I<BitVector> object and return I<BitVector>.

=item B<FlipBitsRange>

    $BitVector->FlipBitsRange($MinBitNum, $MaxBitNum);

Flip values of specified bit numbers between I<MinBitNum> and I<MaxBitNum> in I<BitVector>
object and return I<BitVector>.

=item B<GetBit>

    $BitValue = $BitVector->GetBit($BitNum);

Returns value of bit number I<BitNum> in I<BitVector> object.

=item B<GetBitsAsBinaryString>

    $BitString = $BitVector->GetBitsAsBinaryString([$BitOrder]);

Returns values of bits in I<BitVector> as an ascii bit string containing 0s and 1s.

Default I<BitOrder> is I<Ascending> bit order which corresponds to first bit in each
byte as the loweset bit as opposed to the higest bit.

=item B<GetBitsAsDecimalString>

    $BitString = $BitVector->GetBitsAsDecimalString([$BitOrder]);

Returns values of bits in I<BitVector> as a decimal bit string containing values from 0 to
9.

Default I<BitOrder> is I<Ascending> bit order which corresponds to first bit in each
byte as the loweset bit as opposed to the higest bit.

=item B<GetBitsAsHexadecimalString>

    $BitString = $BitVector->GetBitsAsHexadecimalString([$BitOrder]);

Returns values of bits in I<BitVector> as a hexadecimal bit string containing values from 0 to 9
and a to f.

Default I<BitOrder> is I<Ascending> bit order which corresponds to first bit in each
byte as the loweset bit as opposed to the higest bit.

=item B<GetBitsAsOctalString>

    $BitString = $BitVector->GetBitsAsOctalString([$BitOrder]);

Returns values of bits in I<BitVector> as an octal bit string containing values form 0 to
7.

Default I<BitOrder> is I<Ascending> bit order which corresponds to first bit in each
byte as the loweset bit as opposed to the higest bit.

=item B<GetBitsAsRawBinaryString>

    $BitString = $BitVector->GetBitsAsRawBinaryString();

Returns values of bits in I<BitVector> as an string corresponding to packed bit values
used by Perl vec function without perfoming any unpacking.

=item B<GetDensityOfClearBits>

    $ClearBitsDensity = $BitVector->GetDensityOfClearBits();

Returns density of clear bits in I<BitVector> which corresponds to number of bits set to 0
I<BitVector> divided by its size.

=item B<GetDensityOfSetBits>

    $SetBitsDensity = $BitVector->GetDensityOfSetBits();

Returns density of set bits in I<BitVector> which corresponds to number of bits set to 1 in
I<BitVector> divided by its size.

=item B<GetNumOfClearBits>

    $NumOfClearBits = $BitVector->GetNumOfClearBits();

Returns number of bits set to 0 in I<BitVector>.

=item B<GetNumOfSetBits>

    $NumOfSetBits = $BitVector->GetNumOfSetBits();

Returns number of bits set to 1 in I<BitVector>.

=item B<GetSize>

    $Size = $BitVector->GetSize();

Returns size of I<BitVector>.

=item B<IsBitClear>

    $Status = $BitVector->IsBitClear();

Returns 1 or 0 based on whether I<BitNum> is set to 0 in I<BitVector>.

=item B<IsBitSet>

    $Status = $BitVector->IsBitSet($BitNum);

Returns 1 or 0 based on whether I<BitNum> is set to 1 in I<BitVector>.

=item B<IsBitVector>

    $Status = BitVector::IsBitVector($Object);

Returns 1 or 0 based on whether I<Object> is a B<BitVector> object.

=item B<NewFromBinaryString>

    $NewBitVector = BitVector::NewFromBinaryString($BinaryString,
                    [$BitOrder]);
    $NewBitVector = $BitVector->NewFromBinaryString($BinaryString,
                    [$BitOrder]);

Creates a new I<BitVector> using I<BinaryString> and returns new B<BitVector> object.

Default I<BitOrder> is I<Ascending> bit order which corresponds to first bit in each
byte as the loweset bit as opposed to the higest bit.

=item B<NewFromDecimalString>

    $NewBitVector = BitVector::NewFromDecimalString($DecimalString,
                    [$BitOrder]);
    $NewBitVector = $BitVector->NewFromDecimalString($DecimalString,
                    [$BitOrder]);

Creates a new I<BitVector> using I<DecimalString> and returns new B<BitVector> object.

Default I<BitOrder> is I<Ascending> bit order which corresponds to first bit in each
byte as the loweset bit as opposed to the higest bit.

=item B<NewFromHexadecimalString>

    $NewBitVector = BitVector::NewFromHexadecimalString(
                    $HexadecimalString, [$BitOrder]);
    $NewBitVector = $BitVector->NewFromHexadecimalString(
                    $HexadecimalString, [$BitOrder]);

Creates a new I<BitVector> using I<HexadecimalString> and returns new B<BitVector> object.

Default I<BitOrder> is I<Ascending> bit order which corresponds to first bit in each
byte as the loweset bit as opposed to the higest bit.

=item B<NewFromOctalString>

    $NewBitVector = BitVector::NewFromOctalString($OctalString, [$BitOrder]);
    $NewBitVector = $BitVector->NewFromOctalString($OctalString, [$BitOrder]);

Creates a new I<BitVector> using I<OctalString> and returns new B<BitVector> object.

Default I<BitOrder> is I<Ascending> bit order which corresponds to first bit in each
byte as the loweset bit as opposed to the higest bit.

=item B<NewFromRawBinaryString>

    $NewBitVector = BitVector::NewFromRawBinaryString(
                    $RawBinaryString);
    $NewBitVector = $BitVector->NewFromRawBinaryString(
                    $RawBinaryString);

Creates a new I<BitVector> using I<RawBinaryString> and returns new B<BitVector> object.

=item B<Reverse>

    $BitVector->Reverse();

Reverses values of bits in I<BitVector> and returns I<BitVector>. First bit number ends up with
value of last bit number.

=item B<SetAllBits>

    $BitVector->SetAllBits();

Sets values of all bits in I<BitVector> to 1 and returns I<BitVector>.

=item B<SetBit>

    $BitVector->SetBit($BitNum);

Sets value of I<BitNum> to 1 in I<BitVector> and returns I<BitVector>.

=item B<SetBitValue>

    $BitVector->SetBitValue($BitNum, $BitValue);

Sets value of I<BitNum> to I<BitValue> in I<BitVector> and returns I<BitVector>.

=item B<SetBitValueBitOrder>

    BitVector::SetBitValueBitOrder($BitOrder);
    $BitVector->SetBitValueBitOrder($BitOrder);

Set bit order for printing B<BitVector> values during stringification of B<BitVector> object.
Possible bit order values: I<Ascending or Descending>.

Bit order can be set for either an individual B<BitVector> object or the class. Default is
to print bits in each byte in I<Asscending> bit order.

Internally, bits are stored in I<Ascending> bit order using Perl vec function. Regardless
of machine order, big-endian or little-endian, vec function always considers first
string byte as the lowest byte and first bit within each byte as the lowest bit.

=item B<SetBitValuePrintFormat>

    BitVector::SetBitValuePrintFormat($PrintValueFormat);
    $BitVector->SetBitValuePrintFormat($PrintValueFormat);

Set bit values print format for printing B<BitVector> values during stringification of B<BitVector>
object. Possible print format values: I<Binary, Bin, Hexadecimal, Hex, Decimal, Dec, Octal,
Oct, RawBinary, RawBin>. Default: I<Binary>.

Bit values print format can be set for either an individual B<BitVector> object or the class.

=item B<SetBits>

    $BitVector->SetBits(@BitNums);

Set specified bit numbers I<BitNums> to 1 in I<BitVector> object and return I<BitVector>.

=item B<SetBitsAsBinaryString>

    $BitVector->SetBitsAsBinaryString($BinaryString);

Set bit values in I<BitVector> using specified I<BinaryString> and return I<BitVector>. The
size of I<BitVector> is not changed.

=item B<SetBitsAsDecimalString>

    $BitVector->SetBitsAsDecimalString($DecimalString, [$BitOrder]);

Set bit values in I<BitVector> using specified I<DecimalString> and return I<BitVector>. The
size of I<BitVector> is not changed.

=item B<SetBitsAsHexadecimalString>

    $BitVector->SetBitsAsHexadecimalString($HexadecimalString, [$BitOrder]);

Set bit values in I<BitVector> using specified I<HexadecimalString> and return I<BitVector>. The
size of I<BitVector> is not changed.

=item B<SetBitsAsOctalString>

    $BitVector->SetBitsAsOctalString($OctalString, [$BitOrder]);

Set bit values in I<BitVector> using specified I<OctalString> and return I<BitVector>. The
size of I<BitVector> is not changed.

=item B<SetBitsAsRawBinaryString>

    $BitVector->SetBitsAsRawBinaryString($RawBinaryString);

Set bit values in I<BitVector> using specified I<RawBinaryString> and return I<BitVector>. The
size of I<BitVector> is not changed.

=item B<SetBitsRange>

    $BitVector->SetBitsRange($MinBitNum, $MaxBitNum);

Set specified bit numbers between I<MinBitNum> and I<MaxBitNum> to 1 in I<BitVector>
object and return I<BitVector>.

=item B<StringifyBitVector>

    $String = $BitVector->StringifyBitVector();

Returns a string containing information about I<BitVector> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Vector.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
