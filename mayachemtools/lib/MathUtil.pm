package MathUtil;
#
# File: MathUtil.pm
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
use Math::Trig ();
use POSIX ();

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw(acos asin atan tan ceil floor log10 min max srandom random round GeneratePrimeNumbersUpToLimit GeneratePrimeNumbersUpToCount);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]
	       );


# Return next largest integer...
sub ceil ($) {
  my($Value) = @_;

  return POSIX::ceil($Value);
}

# Return previous smallest integer...
sub floor ($) {
  my($Value) = @_;

  return POSIX::floor($Value);
}

# Calculate log value using base 10...
sub log10 ($) {
  my($Value) = @_;

  return CORE::log($Value)/CORE::log(10);
}

# Return the smaller of two numbers...
sub min ($$) {
  my($Value1, $Value2) = @_;

  return ($Value1 <= $Value2) ? $Value1 : $Value2;
}

# Return the larger of two numbers...
sub max ($$) {
  my($Value1, $Value2) = @_;

  return ($Value1 >= $Value2) ? $Value1 : $Value2;
}

# The random number generator implemented in MayaChemTools is a variant of linear
# congruential generator (LCG) as described by Miller et al. [ Ref 120 ]. It is
# also referred to as Lehmer random number generator or Park-Miller random number
# generator.
#
# Unlike Perl's core random number generator function rand, the random number
# generator implemented in MayaChemTools generates consistent random values
# across different platforms - Windows, CygWin, Linux, Unix - for a specific random
# seed.
#

# $RandomModulus = 2**31 - 1;
# $RandomMultiplier = 16807;
# $RandomQuotient = $RandomModulus / $RandomMultiplier;
# $RandomRemainder = $RandomModulus % $RandomMultiplier
#
# $MaxRandomSeed = 2*31 -2
#
my($MaxRandomSeed, $RandomSeed, $RandomModulus, $RandomMultiplier, $RandomQuotient, $RandomRemainder);

$MaxRandomSeed = 2147483646;
$RandomSeed = 123456789;

$RandomModulus = 2147483647;
$RandomMultiplier = 16807;
$RandomQuotient = 127773;
$RandomRemainder = 2836;

# Set random number seed...
#
# The intial value of random number seed is recommeded to be an integer between 1
# and 2**31 - 2 [Ref 120] which translates to be 1 and 2147483646
#
sub srandom ($) {
  my($Seed) = @_;

  if ($Seed <= 0 ) {
    die "Error: srandom: Specified seed value must be greater than 0...";
  }

  $RandomSeed = ($Seed > $MaxRandomSeed) ? ($Seed % $MaxRandomSeed) : $Seed;

  return $RandomSeed;
}

# Retrun a random number between 0 and less than 1 or specified size...
#
sub random (;$) {
  my($Size) = @_;
  my($Value, $LowValue, $HighValue);

  $Size = defined $Size ? $Size : 1.0;

  $HighValue = $RandomSeed / $RandomQuotient;
  $LowValue = $RandomSeed % $RandomQuotient;

  $Value = $RandomMultiplier * $LowValue - $RandomRemainder * $HighValue;

  $RandomSeed = ($Value > 0) ? $Value : ($Value + $RandomModulus);

  return ($RandomSeed / $RandomModulus) * $Size;
}

# Round a integer/real number to:
# . A nearest integer
# . Specified number of decimal places
#
sub round ($;$) {
  my($Value, $DecimalPlaces) = @_;
  my($RoundedValue);

  if (defined($DecimalPlaces) && $DecimalPlaces > 0) {
    $RoundedValue = sprintf "%.${DecimalPlaces}f", $Value;
  }
  else {
    if ($Value < 0) {
      $RoundedValue = int($Value - 0.5);
    }
    else {
      $RoundedValue = int($Value + 0.5);
    }
  }
  return $RoundedValue;
}

# Return tangent of an angle expressed in radians.
sub tan {
  my($Value) = @_;

  return (CORE::sin($Value)/CORE::cos($Value));
}

# Return inverse sine of an angle expressed in radians.
#
# For a right angle triangle defined by sides X and Y in a unit circle, Pythagorean theorem implies
# X**2 + Y**2 = 1 and sin value corresponds to Y. So asin is equivalent to atan2(Y, sqrt(1-Y**2)).
# However, taking sqrt of negative numbers is problematic; Math::Trig::asin handles it using complex
# numbers.
#
sub asin ($) {
  my($Value) = @_;

  return Math::Trig::asin($Value);
}

# Return inverse cosine of an angle expressed in radians.
#
# For a right angle triangle defined by sides X and Y in a unit circle, Pythagorean theorem implies
# X**2 + Y**2 = 1 and cos value corresponds to X. So asin is equivalent to atan2(sqrt(1-X**2), X)
# However, taking sqrt of negative numbers is problematic; Math::Trig::acos handles it using complex
# numbers.
#
sub acos ($) {
  my($Value) = @_;

  return Math::Trig::acos($Value);
}

# Generate prime numbers up to a specified limit and return a reference to an
# array containing the prime numbers.
#
# By default, the first 1000 prime numbers are generated. The 1000th prime
# number is 7919 and that's why default limit is set to 7920.
#
sub GeneratePrimeNumbersUpToLimit (;$) {
  my($Limit) = @_;

  $Limit = defined $Limit ? $Limit : 7920;

  return _GeneratePrimeNumbers('ByLimit', $Limit)
}

# Generate prime numbers up to specified count of prime numbers and return a
# reference to an array containing the prime numbers.
#
# By default, the first 1000 prime numbers are generated. The 1000th prime
# number is 7919.
#
sub GeneratePrimeNumbersUpToCount (;$) {
  my($Count) = @_;

  $Count = defined $Count ? $Count : 1000;

  return _GeneratePrimeNumbers('ByCount', $Count)
}

# Generate prime numbers up to specified limit or count and return a reference
# to an array containing the prime numbers.
#
# The algorithm to generate prime numbers is a modification of  Sieve of Erastothenes
# prime number generator.
#
sub _GeneratePrimeNumbers {
  my($Mode, $Value) = @_;
  my($ByLimit, $PrimeNumber, $Number, $SqrtOfNumber, $NumberIsPrime, @PrimeNumbers);

  $ByLimit = ($Mode =~ /^ByLimit$/i) ? 1 : 0;

  @PrimeNumbers = (2, 3);
  $Number = 3;

  # while ($Number <= $Limit) {
  while ($ByLimit ? ($Number < $Value) : (@PrimeNumbers < $Value)) {
    $Number += 2;
    $SqrtOfNumber = sqrt $Number;

    $NumberIsPrime = 1;
    PRIMENUMBER: for $PrimeNumber (@PrimeNumbers) {
      if ($PrimeNumber > $SqrtOfNumber) {
	last PRIMENUMBER;
      }
      if (!($Number % $PrimeNumber)) {
	$NumberIsPrime = 0;
	last PRIMENUMBER;
      }
    }
    if ($NumberIsPrime) {
      push @PrimeNumbers, $Number;
    }
  }
  return \@PrimeNumbers;
}

1;

__END__

=head1 NAME

MathUtil

=head1 SYNOPSIS

use MathUtil;

use MathUtil qw(:all);

=head1 DESCRIPTION

B<MathUtil> module provides a variety of common math functions not available in core
Perl package or some other useful math utilities. In order to be consistent with other
Perl functions, name of all the functions in this package are in lowercase which differs
from MayaChemTools naming convention for function names.

B<MathUtil> module provides the following functions:

GeneratePrimeNumbersUpToCount, GeneratePrimeNumbersUpToLimit, acos, asin, ceil,
floor, log10, max, min, random, round, srandom, tan

=head2 FUNCTIONS

=over 4

=item B<GeneratePrimeNumbersUpToCount>

    $PrimesRef = GeneratePrimeNumbersUpToCount();
    $PrimesRef = GeneratePrimeNumbersUpToCount($Count);

Generate prime numbers up to specified I<Count> of prime numbers and return a
reference to an array containing the prime numbers.

By default, the first 1000 prime numbers are generated. The 1000th prime
number is 7919.

The algorithm to generate prime numbers is a modification of  Sieve of Erastothenes
prime number generator.

=item B<GeneratePrimeNumbersUpToLimit>

    $PrimesRef = GeneratePrimeNumbersUpToLimit();
    $PrimesRef = GeneratePrimeNumbersUpToLimit($Limit);

Generate prime numbers up to a specified I<Limit> and return a reference to an
array containing the prime numbers.

By default, the first 1000 prime numbers are generated. The 1000th prime
number is 7919 and that's why default limit is set to 7920.

The algorithm to generate prime numbers is a modification of  Sieve of Erastothenes
prime number generator.

=item B<acos>

    $Value = acos($AngleInRadians);

Returns the nverse cosine of an angle expressed in I<Radians> using Math::Trig::acos
function.

=item B<asin>

    $Value = asin($AngleInRadians);

Returns the inverse sine of an angle expressed in I<Radians> using Math::Trig::asin
function.

=item B<ceil>

    $IntegerValue = ceil($Value);

Returns the next largest integer for I<Value> using POSIX::ceil function.

=item B<floor>

    $IntegerValue = floor($Value);

Returns the previous smallest integer for I<Value> using POSIX::floor function.

=item B<log10>

    $Log10Value = log10($Value);

Returns the log of I<Value> using base 10.

=item B<max>

    $Number = max($Number1, $Number2);

Returns a B<Number> corresponding to the maximum of I<Number1> and I<Number2>.

=item B<min>

    $Number = min($Number1, $Number2);

Returns a B<Number> corresponding to the minimum of I<Number1> and I<Number2>.

=item B<round>

    $RoundedValue = round($Number);
    $RoundedValue = round($Number, $DecimalPlaces);

Returns a value corresponding to a nearst ingeter for I<Number> or formatted to I<DecimalPlaces>.

=item B<random>

    $RandomNumber = random();
    $RandomNumber = random($Size);

Returns a random number between 0 and less than 1 or specified size.

The random number generator implemented in MayaChemTools is a variant of linear
congruential generator (LCG) as described by Miller et al. [ Ref 120 ]. It is
also referred to as Lehmer random number generator or Park-Miller random number
generator.

Unlike Perl's core random number generator function rand, the random number
generator implemented in MayaChemTools generates consistent random values
across different platforms - Windows, CygWin, Linux, Unix - for a specific random
seed.

=item B<srandom>

    $Seed = srandom($Seed);

Sets random number seed to be used by <random> function and returns seed value.

The random number seed is recommeded to be an integer between 1 and 2**31 - 2
[Ref 120] which translates to be 1 and 2147483646.

The default seed is set to 123456789.

=item B<tan>

    $Value = tan($AngleInRadians);

Returns the tangent of an angle expressed in I<Radians>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Constants.pm, ConversionsUtil.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
