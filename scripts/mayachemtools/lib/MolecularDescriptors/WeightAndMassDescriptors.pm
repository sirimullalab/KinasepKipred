package MolecularDescriptors::WeightAndMassDescriptors;
#
# File: WeightAndMassDescriptors.pm
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

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(MolecularDescriptors::MolecularDescriptors Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GetDescriptorNames);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, @DescriptorNames);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyWeightAndMassDescriptors';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeWeightAndMassDescriptors();

  $This->_InitializeWeightAndMassDescriptorsProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Descriptor names...
  @DescriptorNames = ('MolecularWeight', 'ExactMass');

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
sub _InitializeWeightAndMassDescriptors {
  my($This) = @_;

  # Type of MolecularDescriptor...
  $This->{Type} = 'WeightAndMass';

  # Precision for molecular weight and exact mass values...
  $This->{WeightPrecision} = 2;
  $This->{MassPrecision} = 4;

  # Intialize descriptor names and values...
  $This->_InitializeDescriptorNamesAndValues(@DescriptorNames);

  return $This;
}

# Initialize object properties...
#
sub _InitializeWeightAndMassDescriptorsProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  return $This;
}

# Set weight precision for moelcular weight...
#
sub SetWeightPrecision {
  my($This, $Value) = @_;

  if (!TextUtil::IsInteger($Value)) {
    croak "Error: ${ClassName}->SetWeightPrecision: WeightPrecision value, $Value, is not valid:  It must be a an integer...";
  }
  $This->{WeightPrecision} = $Value;

  return $This;
}

# Set mass precision for exact weight...
#
sub SetMassPrecision {
  my($This, $Value) = @_;

  if (!TextUtil::IsInteger($Value)) {
    croak "Error: ${ClassName}->SetMassPrecision: MassPrecision value, $Value, is not valid:  It must be a an integer...";
  }
  $This->{MassPrecision} = $Value;

  return $This;
}

# Generate molecular weight and exact mass values...
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
    carp "Warning: ${ClassName}->GenerateDescriptors: $This->{Type} molecular descriptors generation didn't succeed...";
    return undef;
  }

  # Set final descriptor values...
  $This->_SetFinalDescriptorValues();

  return $This;
}

# Calculate molecular weight and exact mass values..
#
sub _CalculateDescriptorValues {
  my($This) = @_;
  my($MolecularWeight, $ExactMass);

  $MolecularWeight = $This->{Molecule}->GetMolecularWeight();
  $ExactMass = $This->{Molecule}->GetExactMass();

  # Track values...
  $This->{MolecularWeight} = MathUtil::round($MolecularWeight, $This->{WeightPrecision});
  $This->{ExactMass} = MathUtil::round($ExactMass, $This->{MassPrecision});

  return $This;
}

# Setup final descriptor values...
#
sub _SetFinalDescriptorValues {
  my($This) = @_;

  $This->{DescriptorsGenerated} = 1;

  $This->SetDescriptorValues($This->{MolecularWeight}, $This->{ExactMass});

  return $This;
}

# Return a string containg data for WeightAndMassDescriptors object...
#
sub StringifyWeightAndMassDescriptors {
  my($This) = @_;
  my($TheString);

  $TheString = "MolecularDescriptorType: $This->{Type}; " . $This->_StringifyDescriptorNamesAndValues();

  return $TheString;
}

# Is it a WeightAndMassDescriptors object?
sub _IsWeightAndMassDescriptors {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

WeightAndMassDescriptors

=head1 SYNOPSIS

use MolecularDescriptors::WeightAndMassDescriptors;

use MolecularDescriptors::WeightAndMassDescriptors qw(:all);

=head1 DESCRIPTION

B<WeightAndMassDescriptors> class provides the following methods:

new, GenerateDescriptors, GetDescriptorNames, StringifyWeightAndMassDescriptors

B<WeightAndMassDescriptors> is derived from B<MolecularDescriptors> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<WeightAndMassDescriptors>, B<MolecularDescriptors> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

B<WeightAndMassDescriptors> calculates molecular weight and exact mass descriptors using
methods available in B<Molecule> class.

=head2 METHODS

=over 4

=item B<new>

    $NewWeightAndMassDescriptors = new MolecularDescriptors::
                                   WeightAndMassDescriptors(
                                   %NamesAndValues);

Using specified I<WeightAndMassDescriptors> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<WeightAndMassDescriptors>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'WeightAndMass'
    WeightPrecision = 2;
    MassPrecision = 4;

    @DescriptorNames = ('MolecularWeight', 'ExactMass')
    @DescriptorValues = ('None', 'None')

Examples:

    $WeightAndMassDescriptors = new MolecularDescriptors::
                                WeightAndMassDescriptors(
                                'Molecule' => $Molecule);

    $WeightAndMassDescriptors = new MolecularDescriptors::
                                WeightAndMassDescriptors();

    $WeightAndMassDescriptors->SetMolecule($Molecule);
    $WeightAndMassDescriptors->GenerateDescriptors();
    print "WeightAndMassDescriptors: $WeightAndMassDescriptors\n";

=item B<GenerateDescriptors>

    $WeightAndMassDescriptors->GenerateDescriptors();

Calculates molecular weight and exact mass of a molecule and returns
I<WeightAndMassDescriptors>.

=item B<GetDescriptorNames>

    @DescriptorNames = $WeightAndMassDescriptors->GetDescriptorNames();
    @DescriptorNames = MolecularDescriptors::WeightAndMassDescriptors::
                       GetDescriptorNames();

Returns all available descriptor names as an array.

=item B<SetMassPrecision>

    $WeightAndMassDescriptors->SetMassPrecision($Precision);

Sets precision for calculated exact mass value and returns
I<WeightAndMassDescriptors>.

=item B<SetWeightPrecision>

    $WeightAndMassDescriptors->SetWeightPrecision($Precision);

Sets precision for calculated molecular weight value and returns
I<WeightAndMassDescriptors>.

=item B<StringifyWeightAndMassDescriptors>

    $String = $WeightAndMassDescriptors->StringifyWeightAndMassDescriptors();

Returns a string containing information about I<WeightAndMassDescriptors> object.

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
