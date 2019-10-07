package MolecularDescriptors::Fsp3CarbonsDescriptors;
#
# File: Fsp3CarbonsDescriptors.pm
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
use overload '""' => 'StringifyFsp3CarbonsDescriptors';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeFsp3CarbonsDescriptors();

  $This->_InitializeFsp3CarbonsDescriptorsProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Descriptor names...
  @DescriptorNames = ('Fsp3Carbons', 'Sp3Carbons');

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
sub _InitializeFsp3CarbonsDescriptors {
  my($This) = @_;

  # Type of MolecularDescriptor...
  $This->{Type} = 'Fsp3Carbons';

  # Intialize descriptor names and values...
  $This->_InitializeDescriptorNamesAndValues(@DescriptorNames);

  return $This;
}

# Initialize object properties...
#
sub _InitializeFsp3CarbonsDescriptorsProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  return $This;
}

# Calculate fraction of SP3 carbons (Fsp3Carbons)  [ Ref 115-116, Ref 119 ] in a molecule...
#
# It is defined as follows:
#
# Fsp3 = Number of SP3 carbons/Total number of carbons
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
    carp "Warning: ${ClassName}->GenerateDescriptors: $This->{Type} molecular descriptors generation didn't succeed: Couldn't calculate Fsp3Carbons values corresponding to assigned Fsp3Carbons atom types...";
    return undef;
  }

  # Set final descriptor values...
  $This->_SetFinalDescriptorValues();

  return $This;
}

# Calculate Fsp3Carbons value...
#
sub _CalculateDescriptorValues {
  my($This) = @_;
  my($Atom, $AtomID, $TotalCarbons, $NumOfSp3Carbons, $Fsp3Carbons);

  $TotalCarbons = 0;
  $NumOfSp3Carbons = 0;

  ATOM: for $Atom ($This->{Molecule}->GetAtoms()) {
    if (!$Atom->IsCarbon()) {
      next ATOM;
    }
    $TotalCarbons += 1;

    if ($Atom->DoesAtomNeighborhoodMatch('C.T4.TSB4')) {
      $NumOfSp3Carbons += 1;
    }
  }

  $Fsp3Carbons = $NumOfSp3Carbons ? $NumOfSp3Carbons/$TotalCarbons : 0;

  # Track values...
  $This->{Fsp3Carbons} = MathUtil::round($Fsp3Carbons, 2);
  $This->{Sp3Carbons} = $NumOfSp3Carbons;

  return $This;
}

# Setup final descriptor values...
#
sub _SetFinalDescriptorValues {
  my($This) = @_;

  $This->{DescriptorsGenerated} = 1;

  $This->SetDescriptorValues($This->{Fsp3Carbons}, $This->{Sp3Carbons});

  return $This;
}

# Return a string containg data for Fsp3CarbonsDescriptors object...
#
sub StringifyFsp3CarbonsDescriptors {
  my($This) = @_;
  my($Fsp3CarbonsDescriptorsString);

  $Fsp3CarbonsDescriptorsString = "MolecularDescriptorType: $This->{Type}; " . $This->_StringifyDescriptorNamesAndValues();

  return $Fsp3CarbonsDescriptorsString;
}

# Is it a Fsp3CarbonsDescriptors object?
sub _IsFsp3CarbonsDescriptors {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

Fsp3CarbonsDescriptors

=head1 SYNOPSIS

use MolecularDescriptors::Fsp3CarbonsDescriptors;

use MolecularDescriptors::Fsp3CarbonsDescriptors qw(:all);

=head1 DESCRIPTION

B<Fsp3CarbonsDescriptors> class provides the following methods:

new, GenerateDescriptors, GetDescriptorNames, StringifyFsp3CarbonsDescriptors

B<Fsp3CarbonsDescriptors> is derived from B<MolecularDescriptors> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<Fsp3CarbonsDescriptors>, B<MolecularDescriptors> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

Fraction sp3 carbons (Fsp3Carbons) [ Ref 115-116, Ref 119 ] value is calculated by dividing the number of sp3
carbons (Sp3Carbons) with the total number of carbons in a molecule.

=head2 METHODS

=over 4

=item B<new>

    $NewFsp3CarbonsDescriptors = new MolecularDescriptors::
                                 Fsp3CarbonsDescriptors(%NamesAndValues);

Using specified I<Fsp3CarbonsDescriptors> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<Fsp3CarbonsDescriptors>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'Fsp3Carbons'
    @DescriptorNames = ('Fsp3Carbons', 'Sp3Carbons')
    @DescriptorValues = ('None', 'None')

Examples:

    $Fsp3CarbonsDescriptors = new MolecularDescriptors::Fsp3CarbonsDescriptors(
                              'Molecule' => $Molecule);

    $Fsp3CarbonsDescriptors = new MolecularDescriptors::Fsp3CarbonsDescriptors();

    $Fsp3CarbonsDescriptors->SetMolecule($Molecule);
    $Fsp3CarbonsDescriptors->GenerateDescriptors();
    print "Fsp3CarbonsDescriptors: $Fsp3CarbonsDescriptors\n";

=item B<GenerateDescriptors>

    $Fsp3CarbonsDescriptors->GenerateDescriptors();

Calculates Fsp3Carbons and Sp3Carbons values for a molecule and returns I<Fsp3CarbonsDescriptors>.

=item B<GetDescriptorNames>

    @DescriptorNames = $Fsp3CarbonsDescriptors->GetDescriptorNames();
    @DescriptorNames = MolecularDescriptors::Fsp3CarbonsDescriptors::
                       GetDescriptorNames();

Returns all available descriptor names as an array.

=item B<StringifyFsp3CarbonsDescriptors>

    $String = $Fsp3CarbonsDescriptors->StringifyFsp3CarbonsDescriptors();

Returns a string containing information about I<Fsp3CarbonsDescriptors> object.

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
