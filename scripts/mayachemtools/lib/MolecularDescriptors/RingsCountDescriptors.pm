package MolecularDescriptors::RingsCountDescriptors;
#
# File: RingsCountDescriptors.pm
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
use overload '""' => 'StringifyRingsCountDescriptors';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeRingsCountDescriptors();

  $This->_InitializeRingsCountDescriptorsProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Descriptor names...
  @DescriptorNames = ('Rings', 'AromaticRings');

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
sub _InitializeRingsCountDescriptors {
  my($This) = @_;

  # Type of MolecularDescriptor...
  $This->{Type} = 'RingsCount';

  # Intialize descriptor names and values...
  $This->_InitializeDescriptorNamesAndValues(@DescriptorNames);

  return $This;
}

# Initialize object properties...
#
sub _InitializeRingsCountDescriptorsProperties {
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

  $This->{Rings} = $This->{Molecule}->GetNumOfRings();
  $This->{AromaticRings} = $This->{Molecule}->GetNumOfAromaticRings();

  return $This;
}

# Setup final descriptor values...
#
sub _SetFinalDescriptorValues {
  my($This) = @_;

  $This->{DescriptorsGenerated} = 1;

  $This->SetDescriptorValues($This->{Rings}, $This->{AromaticRings});

  return $This;
}

# Return a string containg data for RingsCountDescriptors object...
#
sub StringifyRingsCountDescriptors {
  my($This) = @_;
  my($TheString);

  $TheString = "MolecularDescriptorType: $This->{Type}; " . $This->_StringifyDescriptorNamesAndValues();

  return $TheString;
}

# Is it a RingsCountDescriptors object?
sub _IsRingsCountDescriptors {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

RingsCountDescriptors

=head1 SYNOPSIS

use MolecularDescriptors::RingsCountDescriptors;

use MolecularDescriptors::RingsCountDescriptors qw(:all);

=head1 DESCRIPTION

B<RingsCountDescriptors> class provides the following methods:

new, GenerateDescriptors, GetDescriptorNames, StringifyRingsCountDescriptors

B<RingsCountDescriptors> is derived from B<MolecularDescriptors> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<RingsCountDescriptors>, B<MolecularDescriptors> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

B<RingsCountDescriptors> class doesn't perform any ring or aromaticity detection before
counting their number in a molecule. Instead, it assumes ring and aromaticity detection have
been performed by caller using B<DetectRings> [Ref 31] and B<DetectAromaticity> methods
available in B<Molecule>.

B<DetectAromaticity> method available in B<Molecule> class assigns aromaticity to rings
using Huckel rule as explained below:

o Ring aromaticity is determined using Huckel's rule: a ring containing 4n + 2 pi electrons is
considered aromatic.

o Hetrocyclic rings containing N, O and S atoms fall into two classes: Basic aromatic and
Non-basic aromatic. In Basic aromatic hetrocyclic rings, heteroatom itself is involved in a
double bond. (e.g. Pyridine) However, in non-basic hetrocyclic rings, heteroatom might have
an attached hydrogen atom and the remaining lone pair contribute to electron delocalization
and contributes to 4n + 2 electrons. (e.g. Pyrrole, Furan)

o For molecules containing fused rings, each fused ring set is considered as one aromatic
system for counting pi electrons to satisfy Huckel's rule; In case of a failure, rings in
fused set are treated individually for aromaticity detection. Additionally, non-fused
rings are handled on their own during aromaticity detection.

=head2 METHODS

=over 4

=item B<new>

    $NewRingsCountDescriptors = new MolecularDescriptors::
                                RingsCountDescriptors(
                                %NamesAndValues);

Using specified I<RingsCountDescriptors> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<RingsCountDescriptors>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'RingsCount'

    @DescriptorNames = ('Rings', 'AromaticRings')
    @DescriptorValues = ('None', 'None')

Examples:

    $RingsCountDescriptors = new MolecularDescriptors::RingsCountDescriptors(
                              'Molecule' => $Molecule);

    $RingsCountDescriptors = new MolecularDescriptors::RingsCountDescriptors();

    $RingsCountDescriptors->SetMolecule($Molecule);
    $RingsCountDescriptors->GenerateDescriptors();
    print "RingsCountDescriptors: $RingsCountDescriptors\n";

=item B<GenerateDescriptors>

    $RingsCountDescriptors->GenerateDescriptors();

Calculate number of rings and aromatic rings in a molecule and returns
I<RingsCountDescriptors>.

=item B<GetDescriptorNames>

    @DescriptorNames = $RingsCountDescriptors->GetDescriptorNames();
    @DescriptorNames = MolecularDescriptors::RingsCountDescriptors::
                       GetDescriptorNames();

Returns all available descriptor names as an array.

=item B<StringifyRingsCountDescriptors>

    $String = $RingsCountDescriptors->
                              StringifyRingsCountDescriptors();

Returns a string containing information about I<RingsCountDescriptors> object.

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
