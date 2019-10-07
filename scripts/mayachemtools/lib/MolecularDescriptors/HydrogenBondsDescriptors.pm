package MolecularDescriptors::HydrogenBondsDescriptors;
#
# File: HydrogenBondsDescriptors.pm
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
use overload '""' => 'StringifyHydrogenBondsDescriptors';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeHydrogenBondsDescriptors();

  $This->_InitializeHydrogenBondsDescriptorsProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Descriptor names...
  @DescriptorNames = ('HydrogenBondDonors', 'HydrogenBondAcceptors');

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
sub _InitializeHydrogenBondsDescriptors {
  my($This) = @_;

  # Type of MolecularDescriptor...
  $This->{Type} = 'HydrogenBonds';

  # The currrent release of MayaChemTools supports identification of two types of
  # hydrogen bond donor and acceptor atoms with these names:
  #
  # HBondsType1 or HydrogenBondsType1
  # HBondsType2 or HydrogenBondsType2
  #
  # The names of these hydrogen bond types are rather arbirary. However, their
  # definitions have specific meaning and are as follows:
  #
  # HydrogenBondsType1 [ Ref 60-61, Ref 65-66 ]:
  #   . Donor: NH, NH2, OH - Any N and O with available H
  #   . Acceptor: N[!H], O - Any N without available H and any O
  #
  # HydrogenBondsType2 [ Ref 91 ]:
  #   . Donor: NH, NH2, OH - N and O with availabe H
  #   . Acceptor: N, O - Add N and O
  #
  # Note:
  #   . HydrogenBondsType2 definition corresponds to Rule of 5.
  #
  $This->{HydrogenBondsType} = 'HBondsType2';

  # Intialize descriptor names and values...
  $This->_InitializeDescriptorNamesAndValues(@DescriptorNames);

  return $This;
}

# Initialize object properties...
#
sub _InitializeHydrogenBondsDescriptorsProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  return $This;
}

# Set hydrogen bonds type...
#
sub SetHydrogenBondsType {
  my($This, $HydrogenBondsType) = @_;

  if ($HydrogenBondsType !~ /^(HBondsType1|HBondsType2|HydrogenBondsType1|HydrogenBondsType2)$/i) {
    croak "Error: ${ClassName}->SetHydrogenBondsType: Specified hydrogen bonds type, $HydrogenBondsType, is not supported. Valid values: HBondsType1, HBondsType2, HydrogenBondsType1, HydrogenBondsType2 ...\n ";
  }

  $This->{HydrogenBondsType} = $HydrogenBondsType;

  return $This;
}

# Calculate number of hydrogen bond donors and acceptors in a molecule...
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
    carp "Warning: ${ClassName}->CalculateDescriptorValues: $This->{Type} molecular descriptors generation didn't succeed: Couldn't calculate number of hydrogen bond donor and accepror values...";
    return undef;
  }

  # Set final descriptor values...
  $This->_SetFinalDescriptorValues();

}

# Calculate number of hydrogen bond donors and acceptors...
#
sub _CalculateDescriptorValues {
  my($This) = @_;
  my($HydrogenBondDonors, $HydrogenBondAcceptors, $Atom);

  $HydrogenBondDonors = 0;
  $HydrogenBondAcceptors = 0;

  for $Atom ($This->{Molecule}->GetAtoms()) {
    if ($Atom->IsHydrogenBondDonor($This->{HydrogenBondsType})) {
      $HydrogenBondDonors++;
    }
    if ($Atom->IsHydrogenBondAcceptor($This->{HydrogenBondsType})) {
      $HydrogenBondAcceptors++;
    }
  }

  # Track the calculated values...
  $This->{HydrogenBondDonors} = $HydrogenBondDonors;
  $This->{HydrogenBondAcceptors} = $HydrogenBondAcceptors;

  return $This;
}

# Setup final descriptor values...
#
sub _SetFinalDescriptorValues {
  my($This) = @_;

  $This->{DescriptorsGenerated} = 1;

  $This->SetDescriptorValues($This->{HydrogenBondDonors}, $This->{HydrogenBondAcceptors});

  return $This;
}

# Return a string containg data for HydrogenBondsDescriptors object...
#
sub StringifyHydrogenBondsDescriptors {
  my($This) = @_;
  my($HydrogenBondsDescriptorsString);

  $HydrogenBondsDescriptorsString = "MolecularDescriptorType: $This->{Type}; HydrogenBondsType: $This->{HydrogenBondsType}; " . $This->_StringifyDescriptorNamesAndValues();

  return $HydrogenBondsDescriptorsString;
}

# Is it a HydrogenBondsDescriptors object?
sub _IsHydrogenBondsDescriptors {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

HydrogenBondsDescriptors

=head1 SYNOPSIS

use MolecularDescriptors::HydrogenBondsDescriptors;

use MolecularDescriptors::HydrogenBondsDescriptors qw(:all);

=head1 DESCRIPTION

B<HydrogenBondsDescriptors> class provides the following methods:

new, GenerateDescriptors, GetDescriptorNames, SetHydrogenBondsType,
StringifyHydrogenBondsDescriptors

B<HydrogenBondsDescriptors> is derived from B<MolecularDescriptors> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<HydrogenBondsDescriptors>, B<MolecularDescriptors> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The current release of MayaChemTools supports identification of two types of hydrogen bond
donor and acceptor atoms with these names:

    HBondsType1 or HydrogenBondsType1
    HBondsType2 or HydrogenBondsType2

The names of these hydrogen bond types are rather arbitrary. However, their definitions have
specific meaning and are as follows:

    HydrogenBondsType1 [ Ref 60-61, Ref 65-66 ]:

        Donor: NH, NH2, OH - Any N and O with available H
        Acceptor: N[!H], O - Any N without available H and any O

    HydrogenBondsType2 [ Ref 91 ]:

        Donor: NH, NH2, OH - N and O with available H
        Acceptor: N, O - And N and O

By default, I<HydrogenBondsType2> is used to calculate number hydrogen bond donor
and acceptor atoms. This corresponds to B<RuleOf5> definition of hydrogen bond donors
and acceptors.

=head2 METHODS

=over 4

=item B<new>

    $HydrogenBondsDescriptors = new MolecularDescriptors::
                                HydrogenBondsDescriptors(%NamesAndValues);

Using specified I<HydrogenBondsDescriptors> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<HydrogenBondsDescriptors>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'HydrogenBonds'
    HydrogenBondsType = 'HBondsType2'
    @DescriptorNames = ('HydrogenBondDonors', 'HydrogenBondAcceptors')
    @DescriptorValues = ('None', 'None')

Examples:

    $HydrogenBondsDescriptors = new MolecularDescriptors::
                                HydrogenBondsDescriptors();

    $HydrogenBondsDescriptors = new MolecularDescriptors::
                                HydrogenBondsDescriptors(
                                'HydrogenBondsType' => 'HBondsType2');

    $HydrogenBondsDescriptors->SetMolecule($Molecule);
    $HydrogenBondsDescriptors->GenerateDescriptors();
    print "HydrogenBondsDescriptors: $HydrogenBondsDescriptors\n";

=item B<GenerateDescriptors>

    $HydrogenBondsDescriptors->GenerateDescriptors();

Calculates number of hydrogen bond donors and acceptors a molecule and returns
I<HydrogenBondsDescriptors>.

=item B<GetDescriptorNames>

    @DescriptorNames = $HydrogenBondsDescriptors->GetDescriptorNames();
    @DescriptorNames = MolecularDescriptors::HydrogenBondsDescriptors::
                       GetDescriptorNames();

Returns all available descriptor names as an array.

=item B<SetHydrogenBondsType>

    $HydrogenBondsDescriptors->SetHydrogenBondsType($HBondsType);

Sets value of hydrogen bonds type to use during calculation of descriptors and returns
I<HydrogenBondsDescriptors>. Possible values: I<HBondsType1, HydrogenBondsType1,
HBondsType2, HydrogenBondsType2>.

=item B<StringifyHydrogenBondsDescriptors>

    $String = $HydrogenBondsDescriptors->
                              StringifyHydrogenBondsDescriptors();

Returns a string containing information about I<HydrogenBondsDescriptors> object.

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
