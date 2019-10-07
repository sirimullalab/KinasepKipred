package StatisticsUtil;
#
# File: StatisticsUtil.pm
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

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw(Average AverageDeviation Covariance Correlation Euclidean Factorial FactorialDivison GeometricMean Frequency HarmonicMean KLargest KSmallest Kurtosis Maximum Minimum Mean Median Mode PearsonCorrelation Permutations Product Range RSquare Skewness Sum SumOfSquares StandardDeviation StandardDeviationN  StandardError Standardize StandardScores StandardScoresN TrimMean Variance VarianceN);
@EXPORT_OK = qw();
%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Compute the mean of an array of numbers
sub Average {
  my($XArrayRef) = @_;
  return Mean($XArrayRef);
}

# Compute the average of the absolute deviation of an array of numbers: SUM( ABS(x[i] - Xmean) ) / n
sub AverageDeviation {
  my($XArrayRef) = @_;

  if (!@$XArrayRef) {
    return undef;
  }
  my($AverageDeviation, $Mean, $Value, $SumOfDeviations);

  $AverageDeviation = 0;
  $Mean = Mean($XArrayRef);
  foreach $Value (@$XArrayRef) {
    $SumOfDeviations += abs($Value - $Mean);
  }
  $AverageDeviation = $SumOfDeviations / @$XArrayRef;

  return $AverageDeviation;
}

# Compute correlation coefficient between two arrays of numbers
sub Correlation {
  my($XArrayRef, $YArrayRef) = @_;
  return PearsonCorrelation($XArrayRef, $YArrayRef);
}

# Compute the covariance between two arrays of numbers: SUM( (x[i] - Xmean) (y[i] - Ymean) ) / n
sub Covariance {
  my($XArrayRef, $YArrayRef) = @_;

  if (!(@$XArrayRef && @$YArrayRef && (@$XArrayRef == @$YArrayRef))) {
    return undef;
  }
  my($Covariance, $XMean, $YMean, $Index, $ProductOfDeviations);

  $Covariance = 0;
  $XMean = Mean($XArrayRef);
  $YMean = Mean($YArrayRef);
  $ProductOfDeviations = 0;
  for $Index (0 .. $#{@$XArrayRef}) {
    $ProductOfDeviations += (($XArrayRef->[$Index] - $XMean) * ($YArrayRef->[$Index] - $YMean));
  }
  $Covariance = $ProductOfDeviations / @$XArrayRef;
  return $Covariance;
}

# Compute the euclidean distance of an array of numbers: SQRT( SUM( x[i] ** 2) )
sub Euclidean {
  my($XArrayRef) = @_;

  if (!@$XArrayRef) {
    return undef;
  }
  my($SumOfSquares);

  $SumOfSquares = SumOfSquares($XArrayRef);

  return sqrt $SumOfSquares;
}

# Compute factorial of a number...
sub Factorial {
  my($Num) = @_;

  return _Factorial($Num, 1);
}

# Perform factorial division of two numbers...
sub FactorialDivison {
  my($Numerator, $Denominator) = @_;

  # Only works for integer numbers...
  if ($Numerator <= 0 || ($Numerator != int($Numerator)) ||
     $Denominator <= 0 || ($Denominator != int($Denominator)) ) {
    return undef;
  }
  my($LargerNum, $SmallerNum, $Result);
  $LargerNum = ($Numerator > $Denominator) ? $Numerator : $Denominator;
  $SmallerNum = ($Numerator < $Denominator) ? $Numerator : $Denominator;

  $Result = _Factorial($LargerNum, $SmallerNum);
  if ($Numerator < $Denominator) {
    $Result = 1/$Result;
  }
  return $Result;
}

# Calculate factorial of a number upto a specific limit...
sub _Factorial {
  my($Num, $Limit) = @_;

  # Only works for integer numbers...
  if ($Num <= 0 || ($Num != int($Num)) || $Limit < 1) {
    return undef;
  }

  my($Result) = 1;

  while ($Num > $Limit) {
    $Result *= $Num;
    $Num--;
  }
  return $Result;
}

# Generate all possible permuations or a specific permutations of items in an array
# and return a reference to an array containing array references to generated permuations...
#
# This alogrithm is based on the example provided by Mark Jason-Dominus, and is available
# at CPAN as mjd_permute standalone script.
#
sub Permutations {
  my(@DataToPermute) = @_;
  my($PermutationNum, $NumOfPermutations, @Permutations);

  if (!@DataToPermute) {
    return undef;
  }

  @Permutations = ();
  $NumOfPermutations = Factorial(scalar @DataToPermute);

  for ($PermutationNum = 0; $PermutationNum < $NumOfPermutations; $PermutationNum++) {
    my @Permutation = @DataToPermute[_PermutationNumToPermutation($PermutationNum, $#DataToPermute)];
    push @Permutations, \@Permutation;
  }

  return \@Permutations;
}

# Generte Nth permutation for a collection of specific size...
#
sub _PermutationNumToPermutation {
  my($Num, $Size) = @_;

  return _PatternToPermutation(_PermutationNumToPattern($Num, $Size));
}

# Generate Nth pattern for a collection of specific size...
#
sub _PermutationNumToPattern {
  my($Num, $Size) = @_;
  my($Index, @Pattern);

  $Index = 1;

  while ($Index <= $Size + 1) {
    push @Pattern, $Num % $Index;
    $Num = int($Num/$Index);
    $Index++;
  }

  return @Pattern;
}

# Generate permutation of integers from pattern...
#
sub _PatternToPermutation {
  my(@Pattern) = @_;
  my(@Source, @Permutation);

  @Source = (0 .. $#Pattern);

  while (@Pattern) {
    push @Permutation, splice(@Source, (pop @Pattern), 1);
  }

  return @Permutation;
}

# Compute the frequency of occurance of values in an array of numbers. Three different
# invocation methods are supported:
#
# Frequency(\@ArrayRef) : Using the smallest and largest values, group the numbers into
# 10 bins.
#
# Frequency(\@ArrayRef, $NumOfBins) : Using the smallest and largest values, group the
# numbers into specified bins.
#
# Frequency(\@ArrayRef, \@BinRange): Use bin range to goup the values into different bins.
#
# A hash array is returned with keys and values representing range and frequency values respectively.
# The frequency value for a specific key corresponds to all the values which are greater than
# the previous key and less than or equal to the current key. A key value representing maximum value is
# added for generating frequency distribution for specific number of bins, and whenever the maximum
# array value is greater than the maximum specified in bin range, it is also added to bin range.
#
sub Frequency {
  my($XArrayRef) = @_;

  if (!@$XArrayRef) {
    return undef;
  }

  my($BinRange, $NumOfBins, $BinRangeSpecified);

  $BinRangeSpecified = 0;
  $NumOfBins = 10;
  if (@_ == 2) {
    if (ref($_[1]) eq 'ARRAY') {
      $BinRange = $_[1];
      if (!(@$BinRange && (@$BinRange > 1))) {
	return undef;
      }
      # Make sure the bin range contains values in increasing order...
      my($Index1, $Index2);
      for $Index1 (0 .. $#{@$BinRange}) {
	for $Index2 (($Index1 + 1) .. $#{@$BinRange}) {
	  if ($BinRange->[$Index1] >= $BinRange->[$Index2]) {
	    return undef;
	  }
	}
      }
      $BinRangeSpecified = 1;
    }
    else {
      $NumOfBins = $_[1];
      if ($NumOfBins <= 1) {
	return undef;
      }
    }
  }

  # Setup range keys...
  my(@RangeKeys);
  @RangeKeys = ();

  my($MinValue, $MaxValue) = Range($XArrayRef);
  if ($BinRangeSpecified) {
    push @RangeKeys, @$BinRange;
    if ($MaxValue > $RangeKeys[$#RangeKeys]) {
      push @RangeKeys, $MaxValue;
    }
  }
  else {
    my($MinValue, $MaxValue) = Range($XArrayRef);
    my($Interval) = ($MaxValue - $MinValue)/$NumOfBins;
    my($KeyValue) = $MinValue + $Interval;
    while ($KeyValue  < $MaxValue) {
      push @RangeKeys, $KeyValue;
      $KeyValue += $Interval;
    }
    push @RangeKeys, $MaxValue;
  }

  #Setup frequency hash array...
  my(%FrequencyMap);
  %FrequencyMap = ();

  %FrequencyMap = map { $_ => 0 } @RangeKeys;

  # Count values...
  my($Key, $Value);

  VALUE: for $Value (@$XArrayRef) {
      for $Key (@RangeKeys) {
	if ($Value <= $Key) {
	  $FrequencyMap{$Key} += 1;
	  next VALUE;
	}
      }
  }
  return (%FrequencyMap);
}

# Compute the geometric mean of an array of numbers: NthROOT( PRODUCT(x[i]) )
sub GeometricMean {
  my($XArrayRef) = @_;

  if (!@$XArrayRef) {
    return undef;
  }
  my($Mean, $Product, $Value);
  $Product = 1;
  foreach $Value (@$XArrayRef) {
    if ($Value <= 0 ) {
      return undef;
    }
    $Product *= $Value;
  }
  $Mean = $Product ** (1 / @$XArrayRef);
  return $Mean;
}

# Compute the harmonic mean of an array of numbers: 1 / ( SUM(1/x[i]) / n )
sub HarmonicMean {
  my($XArrayRef) = @_;

  if (!@$XArrayRef) {
    return undef;
  }
  my($Mean, $Sum, $Value);
  $Sum = 0;
  foreach $Value (@$XArrayRef) {
    if ($Value <= 0 ) {
      return undef;
    }
    $Sum += 1/$Value;
  }
  $Mean = 1/($Sum/@$XArrayRef);
  return $Mean;
}

# Return the k-largest value from an array of numbers
sub KLargest {
  my($XArrayRef, $K) = @_;

  if (!(@$XArrayRef && ($K > 0) && ($K <= @$XArrayRef))) {
    return undef;
  }
  my($KLargest, @SortedXArray);
  @SortedXArray = sort { $b <=> $a } @$XArrayRef;
  $KLargest = $SortedXArray[$K - 1];
  return $KLargest;
}

# Return the k-smallest value from an array of numbers
sub KSmallest {
  my($XArrayRef, $K) = @_;

  if (!(@$XArrayRef && ($K > 0) && ($K <= @$XArrayRef))) {
    return undef;
  }
  my($KSmallest, @SortedXArray);
  @SortedXArray = sort { $a <=> $b } @$XArrayRef;
  $KSmallest = $SortedXArray[$K - 1];
  return $KSmallest;
}

# Compute the kurtosis of an array of numbers:
# [ {n(n + 1)/(n - 1)(n - 2)(n - 3)}  SUM{ ((x[i] - Xmean)/STDDEV)^4 } ] - {3((n - 1)^2)}/{(n - 2)(n-3)}
#
sub Kurtosis {
  my($XArrayRef) = @_;

  if (!@$XArrayRef || ((@$XArrayRef - 3) <= 0)) {
    return undef;
  }
  my($Kurtosis, $Mean, $StandardDeviation, $Value);
  $Mean = Mean($XArrayRef);
  if (!defined $Mean) {
    return undef;
  }
  $StandardDeviation = StandardDeviation($XArrayRef);
  if (!(defined $StandardDeviation &&  $StandardDeviation != 0)) {
    return undef;
  }

  my($SumOfScores, $SampleSize);
  $SumOfScores = 0;
  for $Value (@$XArrayRef) {
    $SumOfScores += (($Value - $Mean)/$StandardDeviation) ** 4;
  }
  $SampleSize = @$XArrayRef;
  $Kurtosis = ((($SampleSize * ($SampleSize + 1))/(($SampleSize - 1) * ($SampleSize - 2) * ($SampleSize - 3))) * $SumOfScores) - ((3 * (($SampleSize - 1) ** 2))/(($SampleSize - 2) * ($SampleSize - 3)));
  return $Kurtosis;
}

# Return the smallest value from an array of numbers
sub Minimum {
  my($XArrayRef) = @_;
  return KSmallest($XArrayRef, 1);
}

# Return the largest value from an array of numbers
sub Maximum {
  my($XArrayRef) = @_;
  return KLargest($XArrayRef, 1);
}

# Compute the mean of an array of numbers: SUM( x[i] ) / n
sub Mean {
  my($XArrayRef) = @_;

  if (!@$XArrayRef) {
    return undef;
  }
  my($Mean, $Sum, $Value);
  $Sum = 0;
  foreach $Value (@$XArrayRef) {
    $Sum += $Value;
  }
  $Mean = $Sum / @$XArrayRef;
  return $Mean;
}

# Compute the median value of an array of numbers. For an even number array, it's
# the average of two middle values.
#
# For even values of n: Xsorted[(n - 1)/2 + 1]
# For odd values of n: (Xsorted[n/2] + Xsorted[n/2 + 1])/2
#
sub Median {
  my($XArrayRef) = @_;

  if (!@$XArrayRef) {
    return undef;
  }
  my($Median, @SortedXArray);
  $Median = 0;
  @SortedXArray = sort { $a <=> $b } @$XArrayRef;
  if (@$XArrayRef % 2) {
    my($MidIndex);
    $MidIndex = int(@SortedXArray - 1)/2;
    $Median = $SortedXArray[$MidIndex];
  }
  else {
    # Even number array...
    my($MidPosition);
    $MidPosition = int(@SortedXArray / 2);
    $Median = ($SortedXArray[$MidPosition - 1] + $SortedXArray[$MidPosition]) / 2;
  }
  return $Median;
}

# Return the most frequently occuring value in an array of numbers
sub Mode {
  my($XArrayRef) = @_;

  if (!@$XArrayRef) {
    return undef;
  }
  my($Value, %ValueToCountMap, @CountList, @SortedCountList);
  %ValueToCountMap = ();
  @CountList = ();
  @SortedCountList = ();
  for $Value (@$XArrayRef) {
    if (exists $ValueToCountMap{$Value}) {
      $ValueToCountMap{$Value} += 1;
    }
    else {
      $ValueToCountMap{$Value} = 1;
    }
  }
  for $Value (keys %ValueToCountMap) {
    push @CountList, $ValueToCountMap{$Value};
  }
  @SortedCountList = sort { $b <=> $a } @CountList;

  # Make sure the frequency of mode value is greater than one and check for
  # multiple modes as well...
  #
  my($ModeCount, $ModeValue);
  $ModeCount = $SortedCountList[0];
  if ($ModeCount <= 1) {
    return undef;
  }
  # Get the first mode value...
  VALUE: for $Value (keys %ValueToCountMap) {
    if ($ValueToCountMap{$Value} == $ModeCount) {
      $ModeValue = $Value;
      # Set it to zero to skip it next time...
      $ValueToCountMap{$Value} = 0;
      last VALUE;
    }
  }

  if (wantarray) {
    # Retrieve all the modes...
    my(@Modes, $Count);
    @Modes = ();
    push @Modes, $ModeValue;
    for $Count (@SortedCountList) {
      if ($Count == $ModeCount) {
      VALUE: for $Value (keys %ValueToCountMap) {
	  if ($ValueToCountMap{$Value} == $ModeCount) {
	    push @Modes, $Value;
	    # Set it to zero to skip it next time...
	    $ValueToCountMap{$Value} = 0;
	    last VALUE;
	  }
	}
      }
    }
    return sort {$b <=> $a} @Modes;
  }
  else {
    return $ModeValue;
  }
}


# Compute the Pearson correlation coefficient between two arrays of numbers:
#
#  SUM( (x[i] - Xmean)(y[i] - Ymean) ) / SQRT( SUM( (x[i] - Xmean)^2 )(SUM( (y[i] - Ymean)^2 ))   )
#
# It returns values in the range from -1.0 to 1.0
sub PearsonCorrelation {
  my($XArrayRef, $YArrayRef) = @_;

  if (!(@$XArrayRef && @$YArrayRef && (@$XArrayRef == @$YArrayRef))) {
    return undef;
  }
  my($Correlation, $XMean, $YMean, $Index, $XValueDeviation, $YValueDeviation, $SquareOfXDeviations, $SquareOfYDeviations, $ProductOfDeviations);

  $Correlation = 0;
  $XMean = Mean($XArrayRef);
  $YMean = Mean($YArrayRef);
  $ProductOfDeviations = 0; $SquareOfXDeviations = 0; $SquareOfYDeviations = 0;
  for $Index (0 .. $#{@$XArrayRef}) {
    $XValueDeviation = $XArrayRef->[$Index] - $XMean;
    $YValueDeviation = $YArrayRef->[$Index] - $YMean;
    $ProductOfDeviations += ($XValueDeviation * $YValueDeviation);
    $SquareOfXDeviations += $XValueDeviation ** 2;
    $SquareOfYDeviations += $YValueDeviation ** 2;
  }
  $Correlation = $ProductOfDeviations / sqrt($SquareOfXDeviations *  $SquareOfYDeviations);
  return $Correlation;
}

# Return the smallest and largest values from an array of numbers
sub Range {
  my($XArrayRef) = @_;

  if (!@$XArrayRef) {
    return (undef, undef);
  }
  my($Smallest, $Largest, @SortedXArray);
  @SortedXArray = sort { $a <=> $b } @$XArrayRef;
  $Smallest = $SortedXArray[0];
  $Largest = $SortedXArray[$#SortedXArray];
  return ($Smallest, $Largest);
}

# Compute square of the Pearson correlation coefficient between two arrays of numbers.
#
sub RSquare {
  my($XArrayRef, $YArrayRef) = @_;
  my($RSquare, $Correlation);

  $RSquare = undef;
  $Correlation = PearsonCorrelation($XArrayRef, $YArrayRef);
  if (defined $Correlation) {
    $RSquare = $Correlation ** 2;
  }
  return $RSquare;
}

# Compute the skewness of an array of numbers:
#  {n/(n - 1)(n - 2)} SUM{ ((x[i] - Xmean)/STDDEV)^3 }
#
sub Skewness {
  my($XArrayRef) = @_;

  if (!@$XArrayRef || ((@$XArrayRef - 2) <= 0)) {
    return undef;
  }
  my($Skewness, $Mean, $StandardDeviation, $Value);
  $Mean = Mean($XArrayRef);
  if (!defined $Mean) {
    return undef;
  }
  $StandardDeviation = StandardDeviation($XArrayRef);
  if (!(defined $StandardDeviation &&  $StandardDeviation != 0)) {
    return undef;
  }

  my($SumOfScores, $SampleSize);
  $SumOfScores = 0;
  for $Value (@$XArrayRef) {
    $SumOfScores += (($Value - $Mean)/$StandardDeviation) ** 3;
  }
  $SampleSize = @$XArrayRef;
  $Skewness = ($SampleSize/(($SampleSize - 1) * ($SampleSize - 2) )) * $SumOfScores;
  return $Skewness;
}

# Compute the standard deviation of an array of numbers
sub StandardDeviation {
  my($XArrayRef) = @_;
  return _CalculateStandardDeviation($XArrayRef, 2);
}

# Compute the standard deviation of an array of numbers representing entire population
sub StandardDeviationN {
  my($XArrayRef) = @_;
  return _CalculateStandardDeviation($XArrayRef, 1);
}

# Compute the standard deviation of an array of numbers.
# Mode 1: SQRT ( SUM( (x[i] - mean)^2 ) / n )
# Mode 2: SQRT ( SUM( (x[i] - mean)^2 ) / (n - 1) )
#
sub _CalculateStandardDeviation {
  my($XArrayRef, $Mode) = @_;

  if (!@$XArrayRef) {
    return undef;
  }
  my($StandardDeviation, $Value, $SquareOfDeviations, $Mean, $N);

  $StandardDeviation = 0;
  $Mean = Mean($XArrayRef);
  $SquareOfDeviations = 0;
  foreach $Value (@$XArrayRef) {
    $SquareOfDeviations += ($Value - $Mean) ** 2;
  }
  $N = ($Mode == 1) ? @$XArrayRef : (@$XArrayRef - 1);
  $StandardDeviation =  sqrt($SquareOfDeviations / $N);

  return $StandardDeviation;
}

# Compute the standard error using standard deviation and sample size
sub StandardError {
  my($StandardDeviation, $Count) = @_;

  if ($Count <= 0) {
    return undef;
  }
  my($StandardError);
  $StandardError = $StandardDeviation / sqrt($Count);

  return $StandardError;
}

# Standardize the value using mean and standard deviation
sub Standardize {
  my($Value, $Mean, $StandardDeviation) = @_;

  if ($StandardDeviation <= 0) {
    return undef;
  }
  my($StandardizedValue);
  $StandardizedValue = ($Value - $Mean)/$StandardDeviation;

  return $StandardizedValue;
}

# Compute the standard deviation above the mean for an array of numbers.
sub StandardScores {
  my($XArrayRef) = @_;
  return _CalculateStandardScores($XArrayRef, 2);
}

# Compute the standard deviation above the mean for an array of numbers representing entire population
sub StandardScoresN {
  my($XArrayRef) = @_;
  return _CalculateStandardScores($XArrayRef, 1);
}

# Compute the standard deviation above the mean for an array of numbers.
# Mode 1:  (x[i] - mean) / n
# Mode 2:  (x[i] - mean) / (n - 1)
#
sub _CalculateStandardScores {
  my($XArrayRef, $Mode) = @_;

  if (!@$XArrayRef) {
    return undef;
  }
  my(@StandardScores, $Mean, $StandardDeviation, $Value);

  $Mean = Mean($XArrayRef);
  $StandardDeviation = _CalculateStandardDeviation($XArrayRef, $Mode);
  if (!(defined($StandardDeviation) && $StandardDeviation > 0)) {
    return undef;
  }
  @StandardScores = ();
  for $Value (@$XArrayRef) {
    push @StandardScores, ($Value - $Mean)/$StandardDeviation;
  }

  return @StandardScores;
}

# Compute the product of an array of numbers
sub Product {
  my($XArrayRef) = @_;

  if (!@$XArrayRef) {
    return undef;
  }
  my($Product, $Value);
  $Product = 1;
  foreach $Value (@$XArrayRef) {
    $Product *= $Value;
  }
  return $Product;
}

# Compute the sum of an array of numbers
sub Sum {
  my($XArrayRef) = @_;

  if (!@$XArrayRef) {
    return undef;
  }
  my($Sum, $Value);
  $Sum = 0;
  foreach $Value (@$XArrayRef) {
    $Sum += $Value;
  }
  return $Sum;
}

# Compute the sum of squares of an array of numbers
sub SumOfSquares {
  my($XArrayRef) = @_;

  if (!@$XArrayRef) {
    return undef;
  }
  my($SumOfSquares, $Value);
  $SumOfSquares = 0;
  foreach $Value (@$XArrayRef) {
    $SumOfSquares += $Value ** 2;
  }
  return $SumOfSquares;
}

# Compute the mean of an array of numbers by excluding a fraction of
# numbers from the top and bottom of the data set.
sub TrimMean {
  my($XArrayRef, $FractionToExclude) = @_;

  if (!(@$XArrayRef && $FractionToExclude > 0 && $FractionToExclude <= 1)) {
    return undef;
  }
  my($NumberToExclude);
  $NumberToExclude = int(@$XArrayRef * $FractionToExclude);
  $NumberToExclude = ($NumberToExclude % 2) ? ($NumberToExclude - 1) : $NumberToExclude;
  if ($NumberToExclude == @$XArrayRef) {
    return undef;
  }
  my($Mean, $Sum, $Index, $FirstIndex, $LastIndex);
  $FirstIndex = $NumberToExclude/2;
  $LastIndex = @$XArrayRef - ($NumberToExclude/2) - 1;
  $Sum = 0;
  my(@SortedXArray);
  @SortedXArray = sort { $a <=> $b } @$XArrayRef;
  for $Index ($FirstIndex .. $LastIndex) {
    $Sum += $SortedXArray[$Index];
  }
  $Mean = $Sum/(@SortedXArray - $NumberToExclude);
  return $Mean;
}

# Compute the variance of an array of numbers
sub Variance {
  my($XArrayRef) = @_;
  return _CalculateVariance($XArrayRef, 2);
}

# Compute the variance of an array of numbers representing entire population
sub VarianceN {
  my($XArrayRef) = @_;
  return _CalculateVariance($XArrayRef, 1);
}

# Compute the variance of an array of numbers:
# Mode 1: SUM( (x[i] - Xmean)^2  / n )
# Mode 2: SUM( (x[i] - Xmean)^2  / (n - 1) )
#
sub _CalculateVariance {
  my($XArrayRef, $Mode) = @_;

  if (!@$XArrayRef) {
    return undef;
  }
  my($Variance, $Value, $SquareOfDeviations, $Mean, $N);

  $Variance = 0;
  $Mean = Mean($XArrayRef);
  $SquareOfDeviations = 0;
  foreach $Value (@$XArrayRef) {
    $SquareOfDeviations += ($Value - $Mean) ** 2;
  }
  $N = ($Mode == 1) ? @$XArrayRef : (@$XArrayRef - 1);
  $Variance =  $SquareOfDeviations / $N;

  return $Variance;
}

1;

__END__

=head1 NAME

StatisticsUtil

=head1 SYNOPSIS

use StatisticsUtil;

use Statistics qw(:all);

=head1 DESCRIPTION

B<StatisticsUtil> module provides the following functions:

Average, AverageDeviation, Correlation, Covariance, Euclidean, Factorial,
FactorialDivison, Frequency, GeometricMean, HarmonicMean, KLargest, KSmallest,
Kurtosis, Maximum, Mean, Median, Minimum, Mode, PearsonCorrelation, Permutations,
Product, RSquare, Range, Skewness, StandardDeviation, StandardDeviationN,
StandardError, StandardScores, StandardScoresN, Standardize, Sum, SumOfSquares,
TrimMean, Variance, VarianceN

=head2 METHODS

=over 4

=item B<Average>

    $Value = Average(\@DataArray);

Computes the mean of an array of numbers: SUM( x[i] ) / n

=item B<AverageDeviation>

    $Value = AverageDeviation(\@DataArray);

Computes the average of the absolute deviation of an array of numbers: SUM( ABS(x[i] - Xmean) ) / n

=item B<Correlation>

    $Value = Correlation(\@XDataArray, \@YDataArray);

Computes the Pearson correlation coefficient between two arrays of numbers:
SUM( (x[i] - Xmean)(y[i] - Ymean) ) / SQRT( SUM( (x[i] - Xmean)^2 )(SUM( (y[i] - Ymean)^2 ))   )

=item B<Euclidean>

    $Return = Euclidean(\@DataArray);

Computes the euclidean distance of an array of numbers: SQRT( SUM( x[i] ** 2) )

=item B<Covariance>

    $Value = Covariance(\@XDataArray, \@YDataArray);

Computes the covariance between two arrays of numbers: SUM( (x[i] - Xmean) (y[i] - Ymean) ) / n

=item B<Factorial>

    $Value = Factorial($Num);

Computes the factorial of a positive integer.

=item B<FactorialDivison>

    $Value = FactorialDivision($Numerator, $Denominator);

Compute the factorial divison of two positive integers.

=item B<Frequency>

    %FrequencyValues = Frequency(\@DataArray, [$NumOfBins]);
    %FrequencyValues = Frequency(\@DataArray, [\@BinRange]);

A hash array is returned with keys and values representing range and frequency values, respectively.
The frequency value for a specific key corresponds to all the values which are greater than
the previous key and less than or equal to the current key. A key value representing maximum value is
added for generating frequency distribution for specific number of bins, and whenever the maximum
array value is greater than the maximum specified in bin range, it is also added to bin range.

=item B<GeometricMean>

    $Value = GeometricMean(\@DataArray);

Computes the geometric mean of an array of numbers: NthROOT( PRODUCT(x[i]) )

=item B<HarmonicMean>

    $Value = HarmonicMean(\@DataArray);

Computes the harmonic mean of an array of numbers: 1 / ( SUM(1/x[i]) / n )

=item B<KLargest>

    $Value = KLargest(\@DataArray, $KthNumber);

Returns the k-largest value from an array of numbers.

=item B<KSmallest>

    $Value = KSmallest(\@DataArray, $KthNumber);

Returns the k-smallest value from an array of numbers.

=item B<Kurtosis>

    $Value = Kurtosis(\@DataArray);

Computes the kurtosis of an array of numbers:
[ {n(n + 1)/(n - 1)(n - 2)(n - 3)}  SUM{ ((x[i] - Xmean)/STDDEV)^4 } ] - {3((n - 1)^2)}/{(n - 2)(n-3)}

=item B<Maximum>

    $Value = Maximum(\@DataArray);

Returns the largest value from an array of numbers.

=item B<Minimum>

    $Value = Minimum(\@DataArray);

Returns the smallest value from an array of numbers.

=item B<Mean>

    $Value = Mean(\@DataArray);

Computes the mean of an array of numbers: SUM( x[i] ) / n

=item B<Median>

    $Value = Median(\@DataArray);

Computes the median value of an array of numbers. For an even number array, it's
the average of two middle values.

For even values of n: Xsorted[(n - 1)/2 + 1]
For odd values of n: (Xsorted[n/2] + Xsorted[n/2 + 1])/2

=item B<Mode>

    $Value = Mode(\@DataArray);

Returns the most frequently occuring value in an array of numbers.

=item B<PearsonCorrelation>

    $Value = Correlation(\@XDataArray, \@YDataArray);

Computes the Pearson correlation coefficient between two arrays of numbers:
SUM( (x[i] - Xmean)(y[i] - Ymean) ) / SQRT( SUM( (x[i] - Xmean)^2 )(SUM( (y[i] - Ymean)^2 ))   )

=item B<Permutations>

    $PermutationsRef = Permutations(@DataToPermute);

Generate all possible permuations or a specific permutations of items in an array
and return a reference to an array containing array references to generated permuations.

This alogrithm is based on the example provided by Mark Jason-Dominus, and is available
at CPAN as mjd_permute standalone script.

=item B<Product>

    $Value = Product(\@DataArray);

Compute the product of an array of numbers.

=item B<Range>

    ($Smallest, $Largest) = Range(\@DataArray);

Return the smallest and largest values from an array of numbers.

=item B<RSquare>

    $Value = RSquare(\@XDataArray, \@YDataArray);

Computes square of the Pearson correlation coefficient between two arrays of numbers.

=item B<Skewness>

    $Value = Skewness(\@DataArray);

Computes the skewness of an array of numbers:
{n/(n - 1)(n - 2)} SUM{ ((x[i] - Xmean)/STDDEV)^3 }

=item B<StandardDeviation>

    $Value = StandardDeviation(\@DataArray);

Computes the standard deviation of an array of numbers.
SQRT ( SUM( (x[i] - mean)^2 ) / (n - 1) )

=item B<StandardDeviationN>

    $Value = StandardDeviationN(\@DataArray);

Computes the standard deviation of an array of numbers representing entire population:
SQRT ( SUM( (x[i] - mean)^2 ) / n )

=item B<StandardError>

    $Value = StandardError($StandardDeviation, $Count);

Computes the standard error using standard deviation and sample size.

=item B<Standardize>

    $Value = Standardize($Value, $Mean, $StandardDeviation);

Standardizes the value using mean and standard deviation.

=item B<StandardScores>

    @Values = StandardScores(\@DataArray);

Computes the standard deviation above the mean for an array of numbers:
(x[i] - mean) / (n - 1)

=item B<StandardScoresN>

    @Values = StandardScoresN(\@DataArray);

Computes the standard deviation above the mean for an array of numbers representing entire population:
(x[i] - mean) / n

=item B<Sum>

    $Value = Sum(\@DataArray);

Compute the sum of an array of numbers.

=item B<SumOfSquares>

    $Value = SumOfSquares(\@DataArray);

Computes the sum of an array of numbers.

=item B<TrimMean>

    $Value = TrimMean(\@DataArray, $FractionToExclude));

Computes the mean of an array of numbers by excluding a fraction of
numbers from the top and bottom of the data set.

=item B<Variance>

    $Value = Variance(\@DataArray);

Computes the variance of an array of numbers: SUM( (x[i] - Xmean)^2  / (n - 1) )

=item B<VarianceN>

    $Value = Variance(\@DataArray);

Compute the variance of an array of numbers representing entire population:
SUM( (x[i] - Xmean)^2  / n )

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Constants.pm, ConversionsUtil.pm, MathUtil.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
