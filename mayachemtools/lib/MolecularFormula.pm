package MolecularFormula;
#
# File: MolecularFormula.pm
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
use Text::ParseWords;
use TextUtil;
use PeriodicTable;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(CalculateMolecularWeight CalculateExactMass CalculateElementalComposition FormatCompositionInfomation GetElementsAndCount IsMolecularFormula);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

#
# Calculate molecular weight assuming its a valid molecular formula...
#
sub CalculateMolecularWeight {
  my($MolecularFormula) = @_;
  my($Index, $MolecularWeight, $ElementSymbol, $ElementCount, $AtomicWeight, $FormulaElementsRef, $FormulaElementCountRef);

  ($FormulaElementsRef, $FormulaElementCountRef) = _ProcessMolecularFormula($MolecularFormula);
  if (!(defined($FormulaElementsRef) && defined($FormulaElementCountRef))) {
    return undef;
  }

  $MolecularWeight = 0;

  for $Index (0 .. $#{$FormulaElementsRef}) {
    $ElementSymbol = $FormulaElementsRef->[$Index];
    $ElementCount = $FormulaElementCountRef->[$Index];
    $AtomicWeight = PeriodicTable::GetElementAtomicWeight($ElementSymbol);
    $MolecularWeight += $AtomicWeight * $ElementCount;
  }
  return $MolecularWeight;
}

#
# Calculate exact mass assuming it's a valid formula...
#
sub CalculateExactMass {
  my($MolecularFormula) = @_;
  my($Index, $ElementSymbol, $ElementCount, $ExactMass, $RelativeAtomicMass, $FormulaElementsRef, $FormulaElementCountRef);

  ($FormulaElementsRef, $FormulaElementCountRef) = _ProcessMolecularFormula($MolecularFormula);
  if (!(defined($FormulaElementsRef) && defined($FormulaElementCountRef))) {
    return undef;
  }
  $ExactMass = 0;

  for $Index (0 .. $#{$FormulaElementsRef}) {
    $ElementSymbol = $FormulaElementsRef->[$Index];
    $ElementCount = $FormulaElementCountRef->[$Index];
    $RelativeAtomicMass = PeriodicTable::GetElementMostAbundantNaturalIsotopeMass($ElementSymbol);
    if (!defined($RelativeAtomicMass)) {
      next ELEMENT;
    }
    $ExactMass += $RelativeAtomicMass * $ElementCount;
  }
  return $ExactMass;
}


#
# Calculate elemental composition and return reference to arrays
# containing elements and their percent composition...
#
sub CalculateElementalComposition {
  my($MolecularFormula) = @_;
  my($Index, $MolecularWeight, $ElementSymbol, $ElementCount, $AtomicWeight, $Composition, $CompositionMultiplier, $FormulaElementsRef, $FormulaElementCountRef, @FormulaElements, @FormulaElementComposition);

  $MolecularWeight = CalculateMolecularWeight($MolecularFormula);
  if (! defined $MolecularWeight) {
    return (undef, undef);
  }
  ($FormulaElementsRef, $FormulaElementCountRef) = _ProcessMolecularFormula($MolecularFormula);

  @FormulaElements = ();
  @FormulaElementComposition = ();

  if (!$MolecularWeight) {
    return ( \@FormulaElements, \@FormulaElementComposition);
  }

  $CompositionMultiplier = 100 / $MolecularWeight;

  for $Index (0 .. $#{$FormulaElementsRef}) {
    $ElementSymbol = $FormulaElementsRef->[$Index];
    $ElementCount = $FormulaElementCountRef->[$Index];
    $AtomicWeight = PeriodicTable::GetElementAtomicWeight($ElementSymbol);
    $Composition = ($AtomicWeight * $ElementCount) * $CompositionMultiplier;

    push @FormulaElements, $ElementSymbol;
    push @FormulaElementComposition, $Composition;
  }

  return ( \@FormulaElements, \@FormulaElementComposition);
}

# Using refernece to element and its composition arrays, format composition information
# as: Element: Composition;...
#
sub FormatCompositionInfomation {
  my($Index, $ElementSymbol, $ElementComposition, $ElementsRef, $ElementCompositionRef, $Precision, $Composition);

  $Precision = 2;
  if (@_ == 3) {
    ($ElementsRef, $ElementCompositionRef, $Precision) = @_;
  }
  else {
    ($ElementsRef, $ElementCompositionRef) = @_;
  }

  $Composition = '';
  for $Index (0 .. $#{$ElementsRef}) {
    $ElementSymbol = $ElementsRef->[$Index];
    $ElementComposition = $ElementCompositionRef->[$Index];
    $ElementComposition = sprintf("%.${Precision}f", $ElementComposition);

    $Composition .= ($Composition) ? '; ' : '';
    $Composition .=  "${ElementSymbol}: ${ElementComposition}%";
  }

  return $Composition;
}

#
# Get elements and their count...
#
sub GetElementsAndCount {
  my($MolecularFormula) = @_;
  my($FormulaElementsRef, $FormulaElementCountRef, $ErrorMsg);

  ($FormulaElementsRef, $FormulaElementCountRef, $ErrorMsg) = _ProcessMolecularFormula($MolecularFormula);

  return ($FormulaElementsRef, $FormulaElementCountRef);
}

#
# Is it a valid molecular formula?
#
sub IsMolecularFormula {
  my($MolecularFormula, $PrintErrorMsg, $Status, $FormulaElementsRef, $FormulaElementCountRef, $ErrorMsg);

  ($MolecularFormula) = @_;

  ($FormulaElementsRef, $FormulaElementCountRef, $ErrorMsg) = _ProcessMolecularFormula($MolecularFormula);
  $Status = (defined($FormulaElementsRef) && defined($FormulaElementCountRef)) ? 1 : 0;

  return (wantarray ? ($Status, $ErrorMsg) : $Status);
}

#
# Process molecular formula. For a valid formula, return references to arrays conatining elements
# and element count; otherwsie, return undef.
#
sub _ProcessMolecularFormula {
  my($MolecularFormula) = @_;
  my($ErrorMsg) = '';

  $MolecularFormula = _CleanUpFormula($MolecularFormula);

  # Make sure it only contains numbers and letters...
  if ($MolecularFormula =~ /[^a-zA-Z0-9\(\)\[\]]/) {
    $ErrorMsg = 'Molecular formula contains characters other than a-zA-Z0-9';
    return (undef, undef, $ErrorMsg);
  }

  # Parse the formula...
  my($ElementSpec, $FormulaElementSpec, $Spec, $ElementSymbol, $ElementCount,  @FormulaElements, @ElementCount, %FormulaElementsToCountMap, @SubFormulaElements, %SubFormulaElementsToCountMap);

  @FormulaElements = (); @ElementCount = ();
  %FormulaElementsToCountMap = ();

# Setup element symbol and count regular expression...
# IUPAC: http://www.iupac.org/reports/provisional/abstract04/RB-prs310804/Chap4-3.04.pdf
#

  $FormulaElementSpec = qr/
                   \G(                         # $1
                         (?:
                           ([A-Z][a-z]?)   # Two or one letter element symbol; $2
                           ([0-9]*)          # Optionally followed by element count; $3
                         )
                         | \( | \[
                         | \)[0-9]* | \][0-9]*
                         | .
                      )
                   /x;

  my($ProcessingParenthesis);
  $ProcessingParenthesis = 0;
  # Go over the formula...
  FORMULA: while ($MolecularFormula =~ /$FormulaElementSpec/gx) {
    ($Spec, $ElementSymbol, $ElementCount) = ($1, $2, $3);

    # Handle parenthesis in formula to indicate repeating units...
    if ($Spec =~ /^(\(|\[)/) {
      if ($ProcessingParenthesis) {
	$ErrorMsg = "Molecular formula contains multiple level of () or []";
	return (undef, undef, $ErrorMsg);
      }
      $ProcessingParenthesis = 1;
      @SubFormulaElements = ();
      %SubFormulaElementsToCountMap = ();
      next FORMULA;
    }
    elsif ($Spec =~ /^(\)|\])/) {
      $ProcessingParenthesis = 0;

      # Retrieve repeat count and move data to @FormulaElements and %FormulaElementsToCountMap;
      my($RepeatCount, $Symbol, $Count);
      $RepeatCount = $Spec;
      $RepeatCount =~  s/(\)|\])//g;
      if (!$RepeatCount) {
	$RepeatCount = 1;
      }
      # Copy data...
      for $Symbol (@SubFormulaElements) {
	$Count = $SubFormulaElementsToCountMap{$Symbol} * $RepeatCount;
	_SetupFormulaElementData(\@FormulaElements, \%FormulaElementsToCountMap, $Symbol, $Count);
      }

      # Get ready again...
      @SubFormulaElements = ();
      %SubFormulaElementsToCountMap = ();

      next FORMULA;
    }

    # Retrieve element symbol and count...
    $ElementSymbol = ($Spec && !$ElementSymbol) ? $Spec : ($ElementSymbol ? $ElementSymbol : '');
    $ElementCount = $ElementCount ? $ElementCount : 1;
    if (!PeriodicTable::IsElement($ElementSymbol)) {
      $ErrorMsg = "Molecular formula contains unknown elemental symbol $ElementSymbol";
      return (undef, undef, $ErrorMsg);
    }

    if ($ProcessingParenthesis) {
      _SetupFormulaElementData(\@SubFormulaElements, \%SubFormulaElementsToCountMap, $ElementSymbol, $ElementCount);
    }
    else {
      _SetupFormulaElementData(\@FormulaElements, \%FormulaElementsToCountMap, $ElementSymbol, $ElementCount);
    }
  }

  # Setup element count array...
  for $ElementSymbol (@FormulaElements) {
    $ElementCount = $FormulaElementsToCountMap{$ElementSymbol};
    push @ElementCount, $ElementCount;
  }

  # Make sure it all adds up to 100%; otherwise, adjust the last value..

  return (\@FormulaElements, \@ElementCount, $ErrorMsg);
}

# Clean it up...
sub _CleanUpFormula {
  my($MolecularFormula) = @_;
  #Take out any spaces...
  $MolecularFormula =~ s/ //g;

  # Eliminate any charge specifications: +, - or [1-9]+[+-]
  # e.g NO+ [Al(H2O)6]3+ [H2NO3]+
  if ($MolecularFormula =~ /[\+\-]/) {
    if ($MolecularFormula =~ /\][0-9]+[\+\-]/) {
      # Bracket followed optionally by number and then, +/- ...
      # [Al(H2O)6]3+ ...
      $MolecularFormula =~ s/\][0-9]+[\+\-]/\]/g;
    }
    elsif ($MolecularFormula =~ /[\+\-][0-9]*/) {
      # +/- followed optionally by a number...
      # C37H42N2O6+2, Cu+
      $MolecularFormula =~ s/[\+\-][0-9]*//g;
    }
  }

  # Eliminate any brackets - ] or ) - not followed by numbers:
  # e.g. Li[H2PO4]
  if ($MolecularFormula !~ /\][0-9]+/) {
    $MolecularFormula =~ s/[\[\]]//g;
  }
  if ($MolecularFormula !~ /\)[0-9]+/) {
    $MolecularFormula =~ s/[\(\)]//g;
  }
  # Change adducts to parenthesis format...
  # Na2CO3.10H2O -> Na2CO3(H2O)10
  # 3CdSO4.8H2O -> (CdSO4)3(H2O)8
  if ($MolecularFormula =~ /\./) {
    my($SubFormula, $Count, $Spec);
    my(@MolecularFormulaSplits) = split /\./, $MolecularFormula;
    $MolecularFormula = '';
    for $SubFormula (@MolecularFormulaSplits) {
      ($Count, $Spec) = $SubFormula =~ /^([0-9]*)(.*?)$/;
      if ($Count) {
	$MolecularFormula .= "(${Spec})${Count}";
      }
      else {
	$MolecularFormula .= $Spec;
      }
    }
  }

  return $MolecularFormula;
}

# Store the element and count...
sub _SetupFormulaElementData {
  my($ElementsRef, $ElementsToCountMapRef, $Element, $Count) = @_;

  if (exists $ElementsToCountMapRef->{$Element}) {
    $ElementsToCountMapRef->{$Element} += $Count;
  }
  else {
    push @{$ElementsRef}, $Element;
    $ElementsToCountMapRef->{$Element} = $Count;
  }
}

1;

__END__

=head1 NAME

MolecularFormula

=head1 SYNOPSIS

use MolecularFormula;

use MolecularFormula qw(:all);

=head1 DESCRIPTION

B<MolecularFormula> module provides the following functions:

CalculateElementalComposition, CalculateExactMass, CalculateMolecularWeight,
FormatCompositionInfomation, GetElementsAndCount, IsMolecularFormula

=head1 FUNCTIONS

=over 4

=item B<CalculateMolecularWeight>

    $MolecularWeight = CalculateMolecularWeight($MolecularFormula);

Calculates and returns the molecular weight for a specified I<MolecularFormula>.

=item B<CalculateElementalComposition>

    ($ElementsRef, $ElementCompositionRef) =
       CalculateElementalComposition($MolecularFormula);

Calculates the percent composition in a specified I<MolecularFormula> and returns references
to arrays containing elements and their percent composition.

=item B<CalculateExactMass>

    $ExactMass = CalculateMolecularWeight($MolecularFormula);

Calculates and returns the exact mass for a specified I<MolecularFormula>.

=item B<FormatCompositionInfomation>

    $FormattedString = FormatCompositionInfomation($ElementsRef,
                       $ElementCompositionRef, [$Precision]);

Returns a formatted elemental composition string using references to elements and elemental
composition arrays. Precision is an optional parameter; its default value is I<2>.

=item B<GetElementsAndCount>

    ($ElementsRef, $ElementCountRef) = GetElementsAndCount(
                                       $MolecularFormula);

Retrieves elements and their count composition in a specified I<MolecularFormula> and
returns references to arrays containing elements and their count.

=item B<IsMolecularFormula>

    $Status = IsMolecularFormula($MolecularFormula);
    ($Status, $ErrorMsg) = IsMolecularFormula($MolecularFormula);

Returns 1 or 0 a based on whether it's a valid I<MolecularFormula>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Molecule.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
