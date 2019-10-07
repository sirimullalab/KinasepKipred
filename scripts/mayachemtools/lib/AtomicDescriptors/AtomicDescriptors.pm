package AtomicDescriptors::AtomicDescriptors;
#
# File: AtomicDescriptors.pm
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
use TextUtil ();

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
  my($Class, %PropertyNamesAndValues) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializeAtomicDescriptors();

  $This->_InitializeAtomicDescriptorsProperties(%PropertyNamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeAtomicDescriptors {
  my($This) = @_;

  # Molecule object...
  $This->{Molecule} = '';

  # Type of atomic descriptors...
  $This->{Type} = '';

  # By default, atomic decriptor values are also calculated for hydrogens...
  $This->{IgnoreHydrogens} = 0;

  # Calculated atomic descriptor values hash. Instead of assigning the calculated values to Atom
  # objects, these values are stored in the current object in a hash with atom ID and atomic descriptor
  # values as key/value pairs.
  #
  # Unlike molecular descriptors, no descriptor names are assigned to individual atomic descriptor
  # values.
  #
  %{$This->{DescriptorValues}} = ();
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize object properties....
sub _InitializeAtomicDescriptorsProperties {
  my($This, %PropertiesNamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %PropertiesNamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  return $This;
}

# Initialize descriptor values for all atoms in a molecule...
#
sub _InitializeDescriptorValues {
  my($This) = @_;

  if (!$This->{Molecule}) {
    return $This;
  }

  # Assign 'None' to all atomic descriptor values...
  #
  my($Atom, $AtomID);

  ATOM: for $Atom ($This->{Molecule}->GetAtoms()) {
    $AtomID = $Atom->GetID();
    $This->{DescriptorValues}{$AtomID} = 'None';
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
    croak "Error: ${ClassName}->SetType: Can't change AtomicDescriptors type:  It's already set...";
  }
  $This->{Type} = $Type;

  return $This;
}

# Set specific atomic descriptor value...
#
sub SetDescriptorValue {
  my($This, $Atom, $AtomicDescriptor) = @_;
  my($AtomID);

  $AtomID = $Atom->GetID();
  $This->{DescriptorValues}{$AtomID} = $AtomicDescriptor;

  return $This;
}

# Get specific atomic descriptor value...
#
sub GetDescriptorValue {
  my($This, $Atom) = @_;
  my($AtomID);

  $AtomID = $Atom->GetID();

  return exists $This->{DescriptorValues}{$AtomID} ? $This->{DescriptorValues}{$AtomID} : 'None';
}

# Get calculated atomic descriptor values as a  hash with atom ID and atomic descriptor
# values as key/value pairs...
#
sub GetDescriptorValues {
  my($This) = @_;

  return %{$This->{DescriptorValues}};
}

# Are all atomic descriptor values successfully calculated?
#
# Notes:
#   . Dynamic checking of calculated descriptor values for atoms eliminates the need
#     to check and synchronize valid descriptor values during SetDescriptorValue.
#
sub IsDescriptorsGenerationSuccessful {
  my($This) = @_;
  my($Atom, $DescriptorValue, @Atoms);

  ATOM: for $Atom ($This->{Molecule}->GetAtoms()) {
    if ($This->{IgnoreHydrogens} && $Atom->IsHydrogen()) {
      next ATOM;
    }
    $DescriptorValue = $This->GetDescriptorValue($Atom);
    if ($DescriptorValue =~ /^None$/i) {
      return 0;
    }
  }

  return 1;
}

1;

__END__

=head1 NAME

AtomicDescriptors - AtomicDescriptors class

=head1 SYNOPSIS

use AtomicDescriptors::AtomicDescriptors;

use AtomicDescriptors::AtomicDescriptors qw(:all);

=head1 DESCRIPTION

B<AtomicDescriptors> base class used to derive all other atomic descriptors classes provides the following methods:

new, GetDescriptorValue, GetDescriptorValues,
IsDescriptorsGenerationSuccessful, SetDescriptorValue

B<AtomicDescriptors> class is  derived from B<ObjectProperty> base class which provides methods not
explicitly defined in B<Fingerprints> or B<ObjectProperty> classes using Perl's AUTOLOAD functionality.
These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

=head2 METHODS

=over 4

=item B<new>

    $NewAtomicDescriptors = new AtomicDescriptors::
                            AtomicDescriptors(%NamesAndValues);

Using specified I<AtomicDescriptors> property names and values hash, B<new> method creates a new object
and returns a reference to newly created B<AtomicDescriptors> object. By default, following properties are
initialized:

    Molecule = '';
    Type = '';
    IgnoreHydrogens = 0;

=item B<GetDescriptorValue>

    $Value = $AtomicDescriptors->GetDescriptorValue($Atom);

Returns calculated atomic descriptor I<Value> for specified I<Atom>.

=item B<GetDescriptorValues>

    %Values = $AtomicDescriptors->GetDescriptorValues();

Returns calculated atomic descriptor values for all atoms as a hash with atom ID
and atomic descriptor values as key/value pairs.

=item B<IsDescriptorsGenerationSuccessful>

    $Status = $AtomicDescriptors->
              IsDescriptorsGenerationSuccessful();

Returns 1 or 0 based on whether atomic desctiptors calculations was successful.
For a successful atomic descriptors calculation, all atoms must have a value of other
than a string I<None>.

=item B<SetDescriptorValue>

    $AtomicDescriptors->SetDescriptorValue($Atom, $Value);

Sets specified atomic descriptor I<Value> for I<Atom> and returns I<$AtomicDescriptors>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

demo

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
