package Fingerprints::Fingerprints;
#
# File: Fingerprints.pm
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
use ObjectProperty;
use MathUtil ();
use TextUtil ();
use Fingerprints::FingerprintsBitVector;
use Fingerprints::FingerprintsVector;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(ObjectProperty Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializeFingerprints();

  $This->_InitializeFingerprintsProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeFingerprints {
  my($This) = @_;

  # Molecule object...
  $This->{Molecule} = '';

  # Type of fingerprints...
  $This->{Type} = '';

  # Type of fingerprints vector: FingerprintsBitVector or FingerprintsVector...
  $This->{VectorType} = '';

  # Marks successful generation of fingerprints...
  $This->{FingerprintsGenerated} = 0;

  # Initialize values for FingerprintsBitVector...
  _InitializeClassValuesForFingerprintsBitVector();

  # Initialize values for FingerprintsVector...
  _InitializeClassValuesForFingerprintsVector();
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize class values specific to FingerprintsBitVector...
#
sub _InitializeClassValuesForFingerprintsBitVector {
  my($This) = @_;

  # Size of FingerprintsBitVector...
  $This->{Size} = '';

  # Min/Max sizes used for folding FingerprintsBitVector...
  $This->{MinSize} = '';
  $This->{MaxSize} = '';

  # FingerprintsBitVector...
  $This->{FingerprintsBitVector} = '';
}

# Initialize class values specific to FingerprintsVector...
#
sub _InitializeClassValuesForFingerprintsVector {
  my($This) = @_;

  # Types of FingerprintsVector values: OrderedNumericalValues, NumericalValues or AlphaNumericalValues.
  $This->{FingerprintsVectorType} = '';

  # Fingerprints vector...
  $This->{FingerprintsVector} = '';
}

# Initialize object properties....
sub _InitializeFingerprintsProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  return $This;
}

# Set molecule object and make sure it's not already set...
#
sub SetMolecule {
  my($This, $Molecule) = @_;

  if ($This->{Molecule}) {
    croak "Error: ${ClassName}->SetMolecule: Can't change molecule object:  It's already set...";
  }
  $This->{Molecule} = $Molecule;

  # Weaken the reference to disable increment of reference count...
  Scalar::Util::weaken($This->{Molecule});

  return $This;
}

# Set type and make sure it's not already set...
#
sub SetType {
  my($This, $Type) = @_;

  if ($This->{Type}) {
    croak "Error: ${ClassName}->SetType: Can't change fingerprint type:  It's already set...";
  }
  $This->{Type} = $Type;

  return $This;
}

# Is fingerprints generation successful?
#
# Notes:
#   . The specific fingerprints generation class sets the value of FingerprinsGenerated
#     to 1 after the successful generation of fingerprints; otherwise, it's set to 0.
#
sub IsFingerprintsGenerationSuccessful {
  my($This) = @_;

  return $This->{FingerprintsGenerated} ? 1 : 0;
}

# Set vector type and make sure it's not already set...
#
sub SetVectorType {
  my($This, $VectorType) = @_;

  if ($This->{VectorType}) {
    croak "Error: ${ClassName}->SetVectorType: Can't change fingerprint vector type:  It's already set...";
  }
  if ($VectorType !~ /^(FingerprintsBitVector|FingerprintsVector)$/i) {
    croak "Error: ${ClassName}->SetVectorType: Specified value, $VectorType, for Type is not vaild. Supported types in current release of MayaChemTools: FingerprintsBitVector and FingerprintsVector...";
  }
  $This->{VectorType} = $VectorType;

  return $This;
}

#
# Methods/functions for FingerprintsBitVector...
#

# Initialize fingerprint bit vector...
sub _InitializeFingerprintsBitVector {
  my($This) = @_;

  if ($This->{VectorType} !~ /^FingerprintsBitVector$/i) {
    croak "Error: ${ClassName}->_InitializeFingerprintsBitVector: Can't initialize fingerprints bit vector: VectorType must be FingerprintsBitVector...";
  }

  if ($This->{Size}) {
    $This->{FingerprintsBitVector} = new Fingerprints::FingerprintsBitVector($This->{Size});
  }
  return $This;
}

# Set size...
#
sub SetSize {
  my($This, $Size) = @_;

  if ($This->{MinSize} && $Size < $This->{MinSize}) {
    croak "Error: ${ClassName}->SetSize: Fingerprint size value, $Size, is not valid :  It must be >= $This->{MinSize}...";
  }
  if ($This->{MaxSize} && $Size > $This->{MaxSize}) {
    croak "Error: ${ClassName}->SetSize: Fingerprint size value, $Size, is not valid:  It must be <= $This->{MaxSize}...";
  }

  $This->{Size} = $Size;

  return $This;
}

# Set FingerprintsBitVector object and make sure it's not already set...
#
sub SetFingerprintsBitVector {
  my($This, $FingerprintsBitVector) = @_;

  if ($This->{FingerprintsBitVector}) {
    croak "Error: ${ClassName}->SetFingerprintsBitVector: Can't change FingerprintsBitVector object:  It's already set...";
  }
  $This->{FingerprintsBitVector} = $FingerprintsBitVector;

  return $This;
}

# Fold fingerprints by recursively reducing its size by half untill bit density is greater than or equal to
# specified bit density...
#
sub FoldFingerprintsByBitDensity {
  my($This, $BitDensity) = @_;

  if (!($BitDensity >= 0 && $BitDensity <= 1)) {
    croak "Error: ${ClassName}->FoldFingerprintsByBitDensity: Specified bit density, $BitDensity, is not valid:  It must be > 0 && <= 1 ...";
  }

  return $This->_FoldFingerprintsBitVector('ByDensity', $BitDensity);
}

# Fold fingerprints by recursively reducing its size by half untill size is less than or equal to
# specified size...
#
sub FoldFingerprintsBySize {
  my($This, $Size, $CheckSizeValue) = @_;

  if (!defined $CheckSizeValue) {
    $CheckSizeValue = 1;
  }

  if ($CheckSizeValue) {
    if (!TextUtil::IsPositiveInteger($Size)) {
      croak "Error: ${ClassName}->FoldFingerprintsBySize: Specified size, $Size, is not valid:  It must be a positive integer";
    }
    if (!($Size >= $This->{MinSize} && $Size < $This->{Size})) {
      croak "Error: ${ClassName}->FoldFingerprintsBySize: Specified size, $Size, is not valid:  It must be greater than or equal to minimum size of $This->{MinSize} and less than current size of $This->{Size}...";
    }
    if (!TextUtil::IsNumberPowerOfNumber($Size, 2)) {
      croak "Error: ${ClassName}->FoldFingerprintsBySize: Specified size value, $Size, must be power of 2...";
    }
  }

  return $This->_FoldFingerprintsBitVector('BySize', $Size);
}

# Fold fingerprints bit vector using specified size of bit density...
#
sub _FoldFingerprintsBitVector {
  my($This, $Mode, $Value) = @_;

  if (!$This->{FingerprintsBitVector}) {
    return $This;
  }
  my($FingerprintsBitVector, $FoldedFingerprintsBitVector);

  $FoldedFingerprintsBitVector = '';
  $FingerprintsBitVector = $This->{FingerprintsBitVector};
  MODE: {
    if ($Mode =~ /^BySize$/i) { $FoldedFingerprintsBitVector = $FingerprintsBitVector->FoldFingerprintsBitVectorBySize($Value); last MODE; }
    if ($Mode =~ /^ByDensity$/i) { $FoldedFingerprintsBitVector = $FingerprintsBitVector->FoldFingerprintsBitVectorByDensity($Value); last MODE; }
    $FoldedFingerprintsBitVector = '';
  }
  if ($FoldedFingerprintsBitVector) {
    $This->{FingerprintsBitVector} = $FoldedFingerprintsBitVector;
    $This->{Size} = $FoldedFingerprintsBitVector->GetSize();
  }
  return $This;
}

# Get fingerprints as a binary ascii string containing 0s and 1s...
#
sub GetFingerprintBitsAsBinaryString {
  my($This, $BitOrder) = @_;

  return $This->_GetFingerprintBitsAsString('Binary', $BitOrder);
}

# Get fingerprints as a hexadecimal string...
#
sub GetFingerprintBitsAsHexadecimalString {
  my($This, $BitOrder) = @_;

  return $This->_GetFingerprintBitsAsString('Hexadecimal', $BitOrder);
}

# Get fingerprints as a raw binary string containing packed bit values for each
# byte...
#
sub GetFingerprintBitsAsRawBinaryString {
  my($This, $BitOrder) = @_;

  return $This->_GetFingerprintBitsAsString('RawBinary', $BitOrder);
}

# Get fingerprint bits as a string...
#
sub _GetFingerprintBitsAsString {
  my($This, $Format) = @_;

  if (!$This->{FingerprintsBitVector}) {
    return undef;
  }
  FORMAT : {
    if ($Format =~ /^(Binary|Bin|BinaryString)$/i) { return $This->{FingerprintsBitVector}->GetBitsAsBinaryString(); last FORMAT; }
    if ($Format =~ /^(Hexadecimal|Hex|HexadecimalString)$/i) { return $This->{FingerprintsBitVector}->GetBitsAsHexadecimalString(); last FORMAT; }
    if ($Format =~ /^(RawBinary|RawBin|RawBinaryString)$/i) { return $This->{FingerprintsBitVector}->GetBitsAsRawBinaryString(); last FORMAT; }
    croak "Error: ${ClassName}->_GetFingerprintBitsAsString: Specified bit vector string format, $Format, is not supported. Value values: Binary, Bin, BinaryString, Hexdecimal, Hex, HexadecimalString, RawBinary, RawBin, RawBinaryString...";
  }
  return undef;
}

#
# Methods/functions for FingerprintsVector...
#

# Initialize fingerprint vector...
sub _InitializeFingerprintsVector {
  my($This) = @_;

  if ($This->{VectorType} !~ /^FingerprintsVector$/i) {
    croak "Error: ${ClassName}->_InitializeFingerprintsVector: Can't initialize fingerprints vector: VectorType must be FingerprintsVector...";
  }

  if ($This->{FingerprintsVectorType}) {
    $This->{FingerprintsVector} = new Fingerprints::FingerprintsVector('Type' => $This->{FingerprintsVectorType});
  }
  return $This;
}

# Set FingerprintsVector object and make sure it's not already set...
#
sub SetFingerprintsVector {
  my($This, $FingerprintsVector) = @_;

  if ($This->{FingerprintsVector}) {
    croak "Error: ${ClassName}->SetFingerprintsVector: Can't change FingerprintsVector object:  It's already set...";
  }
  $This->{FingerprintsVector} = $FingerprintsVector;

  return $This;
}

# Types of FingerprintsVector values: OrderedNumericalValues, NumericalValues or AlphaNumericalValues.
#
sub SetFingerprintsVectorType {
  my($This, $FingerprintsVectorType) = @_;

  if ($This->{FingerprintsVectorType}) {
    croak "Error: ${ClassName}->SetFingerprintsVector: Can't change FingerprintsVectorType:  It's already set...";
  }
  if ($FingerprintsVectorType !~ /^(OrderedNumericalValues|NumericalValues|AlphaNumericalValues)$/i) {
    croak "Error: ${ClassName}->SetFingerprintsVectorType: Specified value, $FingerprintsVectorType, for Type is not vaild. Supported types in current release of MayaChemTools: OrderedNumericalValues, NumericalValues or AlphaNumericalValues";
  }
  $This->{FingerprintsVectorType} = $FingerprintsVectorType;

  return $This;
}

# Get fingerprints vector values as an array or reference to an array...
#
sub GetFingerprintsVectorValues {
  my($This) = @_;

  if (!$This->{FingerprintsVector}) {
    return undef;
  }
  return $This->{FingerprintsVector}->GetValues();
}

# Get fingerprints vector value IDs as an array or reference to an array...
#
sub GetFingerprintsVectorValueIDs {
  my($This) = @_;

  if (!$This->{FingerprintsVector}) {
    return undef;
  }
  return $This->{FingerprintsVector}->GetValueIDs();
}

1;

__END__

=head1 NAME

Fingerprints - Fingerprints class

=head1 SYNOPSIS

use Fingerprints::Fingerprints;

use Fingerprints::Fingerprints qw(:all);

=head1 DESCRIPTION

B<Fingerprints> class provides the following methods:

new, FoldFingerprintsByBitDensity, FoldFingerprintsBySize,
GetFingerprintBitsAsBinaryString, GetFingerprintBitsAsHexadecimalString,
GetFingerprintBitsAsRawBinaryString, GetFingerprintsVectorValueIDs,
GetFingerprintsVectorValues, IsFingerprintsGenerationSuccessful,
SetFingerprintsBitVector, SetFingerprintsVector, SetFingerprintsVectorType,
SetMolecule, SetSize, SetType, SetVectorType

B<Fingerprints> class is used as a base class for various specific fingerprint classes such as
B<AtomNeighborhoodsFingerprints>, B<AtomTypesFingerprints>, B<EStateIndiciesFingerprints>,
B<PathLengthFingerprints>, B<ExtendedConnectivityFingerprints>, B<MACCSKeys> and so on.
It implements functionality common to fingerprint classes.

B<Fingerprints> class is  derived from B<ObjectProperty> base class which provides methods not
explicitly defined in B<Fingerprints> or B<ObjectProperty> classes using Perl's AUTOLOAD functionality.
These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

B<Fingerprints> class uses B<FingerprintsBitVector> class to provide bits manipulation functionality.

=head2 METHODS

=over 4

=item B<new>

    $NewFingerprints = new Fingerprints(%NamesAndValues);

Using specified I<Fingerprints> property names and values hash, B<new> method creates a new object
and returns a reference to newly created B<Fingerprints> object. By default, following properties are
initialized:

    Molecule = '';
    Type = '';
    VectorType = '';
    Size = '';
    MinSize = '';
    MaxSize = '';
    FingerprintsBitVector = '';
    FingerprintsVectorType = '';
    FingerprintsVector = '';

=item B<FoldFingerprintsByBitDensity>

    $Fingerprints->FoldFingerprintsByBitDensity($BitDensity);

Folds fingerprints by recursively reducing its size by half until bit density is greater than or equal to
specified I<BitDensity> and returns I<Fingerprints>.

=item B<FoldFingerprintsBySize>

    $Fingerprints->FoldFingerprintsBySize($Size, [$CheckSizeValue]);

Fold fingerprints by recursively reducing its size by half until size is less than or equal to specified
I<Size> and returns I<Fingerprints>. By default, value I<Size> is checked to make sure it's:

    >= MinSize and < Size and IsPowerOfTwo

=item B<GetFingerprintBitsAsBinaryString>

    $BinaryASCIIString =
       $Fingerprints->GetFingerprintBitsAsBinaryString();

Returns fingerprints as a binary ASCII string containing 0s and 1s.

=item B<GetFingerprintBitsAsHexadecimalString>

    $HexadecimalString =
       $Fingerprints->GetFingerprintBitsAsHexadecimalString();

Returns fingerprints as a hexadecimal string

=item B<GetFingerprintBitsAsRawBinaryString>

    $RawBinaryString =
       $Fingerprints->GetFingerprintBitsAsRawBinaryString();

Returns fingerprints as a raw binary string containing packed bit values for each byte.

=item B<GetFingerprintsVectorValueIDs>

    $ValueIDsRef = $Fingerprints->GetFingerprintsVectorValueIDs();
    @ValueIDs = $Fingerprints->GetFingerprintsVectorValueIDs();

Returns fingerprints vector value IDs as an array or reference to an array.

=item B<GetFingerprintsVectorValues>

    $ValuesRef = $Fingerprints->GetFingerprintsVectorValues();
    @Values = $Fingerprints->GetFingerprintsVectorValues();

Returns fingerprints vector values as an array or reference to an array.

=item B<IsFingerprintsGenerationSuccessful>

    $Return = $Fingerprints->IsFingerprintsGenerationSuccessful();

Returns 1 or 0 based on whether fingerprints were successfully generated.

=item B<SetFingerprintsBitVector>

    $Fingerprints->SetFingerprintsBitVector($FingerprintsBitVector);

Sets I<FingerprintsBitVector> object for I<Fingerprints> and returns I<Fingerprints>.

=item B<SetFingerprintsVector>

    $Fingerprints->SetFingerprintsVector();

Sets I<FingerprintsVector> object for I<Fingerprints> and returns I<Fingerprints>.

=item B<SetFingerprintsVectorType>

    $Fingerprints->SetFingerprintsVectorType($VectorType);

Sets I<FingerprintsVector> type for I<Fingerprints> and returns I<Fingerprints>. Possible
I<VectorType> values: I<OrderedNumericalValues, NumericalValues or AlphaNumericalValues>.

=item B<SetMolecule>

    $Fingerprints->SetMolecule($Molecule);

Sets I<Molecule> object for I<Fingerprints> and returns I<Fingerprints>.

=item B<SetSize>

    $Fingerprints->SetSize($Size);

Sets I<Size> of fingerprints and returns I<Fingerprints>.

=item B<SetType>

    $Fingerprints->SetType($Type);

Sets I<Type> of fingerprints and returns I<Fingerprints>.

=item B<SetVectorType>

    $Fingerprints->SetVectorType($Type);

Sets I<Type> of fingerprints vector and returns I<Fingerprints>. Possible I<Type> values:
I<FingerprintsBitVector or FingerprintsVector>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

FingerprintsStringUtil.pm, AtomNeighborhoodsFingerprints.pm, AtomTypesFingerprints.pm,
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
