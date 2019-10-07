package ConversionsUtil;
#
# File: ConversionsUtil.pm
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
use Exporter;
use Constants;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);

# Groups of conversion functions...
my(@MathConversions) = qw(DegreesToRadians RadiansToDegrees);
my(@NumericBaseConversions) = qw(BinaryToDecimal DecimalToBinary HexadecimalToDecimal DecimalToHexadecimal OctalToDecimal DecimalToOctal BinaryToHexadecimal HexadecimalToBinary HexadecimalToOctal OctalToHexadecimal StringToBinary StringToHexadecimal);

# Export all conversion functions...
@EXPORT = (@MathConversions, @NumericBaseConversions);
@EXPORT_OK = qw();

%EXPORT_TAGS = (
		math => [@MathConversions],
		all  => [@EXPORT, @EXPORT_OK]
	       );


# Degrees to radians...
sub DegreesToRadians ($;$) {
  my($Degrees, $IgnoreWrap) = @_;
  my($Radians, $WrapValue);

  $WrapValue = (defined($IgnoreWrap) && $IgnoreWrap) ? 0 : 1;
  if ($Degrees > 360 && $WrapValue) {
    $Degrees = $Degrees % 360;
  }
  $Radians = ($Degrees * TwoPi) / 360;

  return $Radians;
}

# Radians to degrees...
sub RadiansToDegrees ($;$) {
  my($Radians, $IgnoreWrap) = @_;
  my($Degrees, $WrapValue);

  $WrapValue = (defined($IgnoreWrap) && $IgnoreWrap) ? 0 : 1;
  $Degrees = ($Radians * 360) / TwoPi;
  if ($Degrees > 360 && $WrapValue) {
    $Degrees = $Degrees % 360;
  }

  return $Degrees;
}

# Convert a binary string to a decimal number...
sub BinaryToDecimal ($) {
  my($Value) = @_;

  if ($Value !~ /^0b/) {
    $Value = "0b${Value}";
  }
  return _ConvertToDecimal($Value);
}

# Convert a decimal number into a binary string...
sub DecimalToBinary ($) {
  my($Value) = @_;

  return sprintf("%b", $Value);
}

# Convert a hexadecimal string to a decimal number...
sub HexadecimalToDecimal ($) {
  my($Value) = @_;

  if ($Value !~ /^0x/) {
    $Value = "0x${Value}";
  }
  return _ConvertToDecimal($Value);
}

# Convert a decimal number into a hexadecimal string...
sub DecimalToHexadecimal ($) {
  my($Value) = @_;

  return sprintf("%x", $Value);
}

# Convert an octal string to a decimal number...
sub OctalToDecimal ($) {
  my($Value) = @_;

  if ($Value !~ /^0/) {
    $Value = "0${Value}";
  }
  return _ConvertToDecimal($Value);
}

# Convert a decimal number into an octal string...
sub DecimalToOctal ($) {
  my($Value) = @_;

  return sprintf("%o", $Value);
}

# Convert string into a binary string. Going from left to right, two ways of arranging bits
# inside each byte are available: Most Significat Bits (MSB) first or Least Significat Bits (LSB)
# first. Default is MSB corresponding to  descending bits order (PerlSpeak) inside each
# each packed byte (Most singificat bits first).
#
sub StringToBinary ($;$) {
  my($Value, $UseReverseBitOrder) = @_;
  my($BinTemplate);

  $BinTemplate = (defined($UseReverseBitOrder) && $UseReverseBitOrder) ? 'b*' : 'B*';
  return unpack($BinTemplate, $Value);
}

# Convert string into a hexadecimal string. Two ways of arranging nybbles (pair of 4 bits in each
# byte) are available: high nybbles first or low nybbles first. Default is MSB corresponding to high
# nybbles (PerlSpeak) first. Low and high nybbles correspond to pair of a low and high four bits in a byte.
#
sub StringToHexadecimal ($;$) {
  my($Value, $UseReverseBitOrder) = @_;
  my($HexTemplate);

  $HexTemplate = (defined($UseReverseBitOrder) && $UseReverseBitOrder) ? 'h*' : 'H*';
  return unpack($HexTemplate, $Value);
}

# Convert a binary string into a hexadecimal string...
#
sub BinaryToHexadecimal ($;$) {
  my($Value, $UseReverseBitOrder) = @_;
  my($BinTemplate, $HexTemplate);

  $BinTemplate = (defined($UseReverseBitOrder) && $UseReverseBitOrder) ? 'b*' : 'B*';
  $HexTemplate = (defined($UseReverseBitOrder) && $UseReverseBitOrder) ? 'h*' : 'H*';

  return unpack($HexTemplate, pack($BinTemplate, $Value));
}

# Convert a hexadecimal string into a binary string...
#
sub HexadecimalToBinary ($;$) {
  my($Value, $UseReverseBitOrder) = @_;
  my($BinTemplate, $HexTemplate);

  $BinTemplate = (defined($UseReverseBitOrder) && $UseReverseBitOrder) ? 'b*' : 'B*';
  $HexTemplate = (defined($UseReverseBitOrder) && $UseReverseBitOrder) ? 'h*' : 'H*';

  return unpack($BinTemplate, pack($HexTemplate, $Value));
}

# Convert a hexadecimal string into a octal string...
#
sub HexadecimalToOctal {
  my($Hexadecimal) = @_;

  return DecimalToOctal(HexadecimalToDecimal($Hexadecimal));
}

# Convert a octal string into a hexadecimal string...
#
sub OctalToHexadecimal {
  my($Octal) = @_;

  return DecimalToHexadecimal(OctalToDecimal($Octal));
}

# Use Perl oct function to convert binary, octal or hexadecimal strings into decimal numbers.
sub _ConvertToDecimal ($) {
  my($Value) = @_;

  return ($Value =~ /^0/) ? oct($Value) : $Value;
}

1;

__END__

=head1 NAME

ConversionsUtil

=head1 SYNOPSIS

use ConversionsUtil;

use ConversionsUtil qw(:math);

use ConversionsUtil qw(:all);

=head1 DESCRIPTION

B<ConversionsUtil> module provides the following functions:

BinaryToDecimal, BinaryToHexadecimal, DecimalToBinary, DecimalToHexadecimal,
DecimalToOctal, DegreesToRadians, HexadecimalToBinary, HexadecimalToDecimal,
HexadecimalToOctal, OctalToDecimal, OctalToHexadecimal, RadiansToDegrees,
StringToBinary, StringToHexadecimal

=head2 FUNCTIONS

=over 4

=item B<BinaryToDecimal>

    $Decimal = BinaryToDecimal($Binary);

Converts a I<Binary> string to B<Decimal> string.

=item B<BinaryToHexadecimal>

    $Hexadecimal = BinaryToHexadecimal($Binary);

Converts a I<Binary> string to B<Hexadecimal> string.

=item B<DecimalToBinary>

    $Binary = DecimalToBinary($Decimal);

Converts a I<Decimal> string to B<Binary> string.

=item B<DecimalToHexadecimal>

    $Hexadecimal = DecimalToHexadecimal($Decimal);

Converts a I<Decimal> string to B<Hexadecimal> string.

=item B<DecimalToOctal>

    $Octal = DecimalToOctal($Decimal);

Converts a I<Decimal> string to B<Octal> string.

=item B<DegreesToRadians>

    $Radians = DegreesToRadians($Degrees, [$DoNotWrapValue]);

Converts degrees to radians in the range from 0 to 2PI or to corresponding radians without
wrapping the converted value to 0 to 2PI. Default is to wrap the converted value.

=item B<HexadecimalToBinary>

    $Binary = HexadecimalToBinary($Hexadecimal);

Converts a I<Hexadecimal> string to B<Binary> string.

=item B<HexadecimalToDecimal>

    $Decimal = HexadecimalToDecimal($Hexadecimal);

Converts a I<Hexadecimal> string to B<Decimal> string.

=item B<HexadecimalToOctal>

    $Octal = HexadecimalToOctal($Hexadecimal);

Converts a I<Hexadecimal> string to B<Octal> string.

=item B<OctalToDecimal>

    $Decimal = OctalToDecimal($Octal);

Converts a I<Octal> string to B<Decimal> string.

=item B<OctalToHexadecimal>

    $Hexadecimal = OctalToHexadecimal($Octal);

Converts a I<Octal> string to B<Hexadecimal> string.

=item B<RadiansToDegrees>

    $Degrees = RadiansToDegrees($Radians, [$DoNotWrapValue]);

Converts radians to degrees in the range from 0 to 360 or to corresponding degrees without
wrapping the converted value to 0 to 360. Default is to wrap the converted value.

=item B<StringToBinary>

    $BinaryString = StringToBinary($String, [$UseReverseBitOrder]);

Converts specified I<String> into a B<Binarystring>. Going from left to right, two ways of arranging
bits inside each byte are available: Most Significat Bits (MSB) first or Least Significat Bits (LSB) first.
Default is MSB corresponding to  descending bits order (PerlSpeak) inside each each packed byte
(Most singificat bits first).

=item B<StringToHexadecimal>

    $HexadecimalString = StringToHexadecimal($String,
                         [$UseReverseBitOrder]);

Convert string into a hexadecimal string. Two ways of arranging nybbles (pair of 4 bits in each
byte) are available: high nybbles first or low nybbles first. Default is MSB corresponding to high
nybbles (PerlSpeak) first. Low and high nybbles correspond to pair of a low and high four bits in a byte.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Constants.pm, MathUtil.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
