package MolecularDescriptors::AtomsCountDescriptors;
#
# File: AtomsCountDescriptors.pm
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
use overload '""' => 'StringifyAtomsCountDescriptors';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeAtomsCountDescriptors();

  $This->_InitializeAtomsCountDescriptorsProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Descriptor names...
  @DescriptorNames = ('Atoms', 'HeavyAtoms');

}

# Get descriptor names as an array...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetDescriptorNames {
  return @DescriptorNames;
}

# Initialize object data...
#
sub _InitializeAtomsCountDescriptors {
  my($This) = @_;

  # Type of MolecularDescriptor...
  $This->{Type} = 'AtomsCount';

  # Intialize descriptor names and values...
  $This->_InitializeDescriptorNamesAndValues(@DescriptorNames);

  return $This;
}

# Initialize object properties...
#
sub _InitializeAtomsCountDescriptorsProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

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

  $This->{Atoms} = $This->{Molecule}->GetNumOfAtoms();
  $This->{HeavyAtoms} = $This->{Molecule}->GetNumOfHeavyAtoms();

  return $This;
}

# Setup final descriptor values...
#
sub _SetFinalDescriptorValues {
  my($This) = @_;

  $This->{DescriptorsGenerated} = 1;

  $This->SetDescriptorValues($This->{Atoms}, $This->{HeavyAtoms});

  return $This;
}

# Return a string containg data for AtomsCountDescriptors object...
#
sub StringifyAtomsCountDescriptors {
  my($This) = @_;
  my($TheString);

  $TheString = "MolecularDescriptorType: $This->{Type}; " . $This->_StringifyDescriptorNamesAndValues();

  return $TheString;
}

# Is it a AtomsCountDescriptors object?
sub _IsAtomsCountDescriptors {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

AtomsCountDescriptors

=head1 SYNOPSIS

use MolecularDescriptors::AtomsCountDescriptors;

use MolecularDescriptors::AtomsCountDescriptors qw(:all);

=head1 DESCRIPTION

B<AtomsCountDescriptors> class provides the following methods:

new, GenerateDescriptors, GetDescriptorNames, StringifyAtomsCountDescriptors

B<AtomsCountDescriptors> is derived from B<MolecularDescriptors> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<AtomsCountDescriptors>, B<MolecularDescriptors> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

B<AtomsCountDescriptors> class counts the number of atoms and heavy atoms in a molecule
corresponding to total number of atom and non-hydrogen atoms respectively.

=head2 METHODS

=over 4

=item B<new>

    $NewAtomsCountDescriptors = new MolecularDescriptors::
                                AtomsCountDescriptors(
                                %NamesAndValues);

Using specified I<AtomsCountDescriptors> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<AtomsCountDescriptors>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'AtomsCount'

    @DescriptorNames = ('Atoms', 'HeavyAtoms')
    @DescriptorValues = ('None', 'None')

Examples:

    $AtomsCountDescriptors = new AtomsCountDescriptors(
                              'Molecule' => $Molecule);

    $AtomsCountDescriptors = new AtomsCountDescriptors();

    $AtomsCountDescriptors->SetMolecule($Molecule);
    $AtomsCountDescriptors->GenerateDescriptors();
    print "AtomsCountDescriptors: $AtomsCountDescriptors\n";


=item B<GenerateDescriptors>

    $AtomsCountDescriptors->GenerateDescriptors();

Calculates number of atoms and heavy atoms in a molecule and returns I<AtomsCountDescriptors>.

=item B<GetDescriptorNames>

    @DescriptorNames = $AtomsCountDescriptors->GetDescriptorNames();
    @DescriptorNames = MolecularDescriptors::AtomsCountDescriptors::
                       GetDescriptorNames();

Returns all available descriptor names as an array.

=item B<StringifyAtomsCountDescriptors>

    $String = $AtomsCountDescriptors->StringifyAtomsCountDescriptors();

Returns a string containing information about I<AtomsCountDescriptors> object.

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
