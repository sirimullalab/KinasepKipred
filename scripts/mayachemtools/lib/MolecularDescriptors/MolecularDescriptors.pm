package MolecularDescriptors::MolecularDescriptors;
#
# File: MolecularDescriptors.pm
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
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializeMolecularDescriptors();

  $This->_InitializeMolecularDescriptorsProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeMolecularDescriptors {
  my($This) = @_;

  # Molecule object...
  $This->{Molecule} = '';

  # Type of molecular descriptors...
  $This->{Type} = '';

  # Names and calculated value of molecular descriptors...
  #
  # The specific descriptor class, derived from this base class, populate descriptor names and values
  # arrays...
  #
  @{$This->{DescriptorNames}} = ();
  @{$This->{DescriptorValues}} = ();

  # Marks successful generation of descriptors...
  $This->{DescriptorsGenerated} = 0;

}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}


# Initialize object properties....
sub _InitializeMolecularDescriptorsProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  return $This;
}

# Initialize descriptor names and values...
#
sub _InitializeDescriptorNamesAndValues {
  my($This, @Names) = @_;

  @{$This->{DescriptorNames}} = @Names;

  $This->_InitializeDescriptorValues();

  return $This;
}

# Initialize descriptor values...
#
sub _InitializeDescriptorValues {
  my($This) = @_;

  $This->{DescriptorsGenerated} = 0;

  @{$This->{DescriptorValues}} = ();

  return $This;
}

# Set molecule object...
#
sub SetMolecule {
  my($This, $Molecule) = @_;

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
    croak "Error: ${ClassName}->SetType: Can't change MolecularDescriptors type:  It's already set...";
  }
  $This->{Type} = $Type;

  return $This;
}

# Get molecular descriptor names as an array...
#
sub GetDescriptorNames {
  my($This) = @_;

  return @{$This->{DescriptorNames}};
}

# Set descriptor names...
#
sub SetDescriptorNames {
  my($This, @Names) = @_;

  @{$This->{DescriptorNames}} = @Names;

  return $This;
}

# Add descriptor names...
#
sub AddDescriptorNames {
  my($This, @Names) = @_;

  push @{$This->{DescriptorNames}}, @Names;

  return $This;
}

# Set descriptor values...
#
sub SetDescriptorValues {
  my($This, @Values) = @_;

  @{$This->{DescriptorValues}} = @Values;

  return $This;
}

# Add descriptor values...
#
sub AddDescriptorValues {
  my($This, @Values) = @_;

  push @{$This->{DescriptorValues}}, @Values;

  return $This;
}

# Is descriptors generation successful?
#
# Notes:
#   . The specific molecular descriptor class generation class sets the value of
#     DescriptorsCalculated  to 1 after the successful generation of descriptors;
#     otherwise, it's set to 0.
#
sub IsDescriptorsGenerationSuccessful {
  my($This) = @_;

  return $This->{DescriptorsGenerated} ? 1 : 0;
}

# Get all descriptor values as an array...
#
sub GetDescriptorValues {
  my($This) = @_;

  if ($This->{DescriptorsGenerated}) {
    return wantarray ? @{$This->{DescriptorValues}} : scalar @{$This->{DescriptorValues}};
  }
  else {
    my(@DescriptorValues);

    @DescriptorValues = ('None') x scalar @{$This->{DescriptorNames}};

    return wantarray ? @DescriptorValues : scalar @DescriptorValues;
  }
}

# Get descriptor value for a specified descriptor name...
#
sub GetDescriptorValueByName {
  my($This, $Name) = @_;
  my(%NamesAndValues);

  %NamesAndValues = $This->GetDescriptorNamesAndValues();

  return exists $NamesAndValues{$Name} ? $NamesAndValues{$Name} : 'None';

}

# Get calculated molecular descriptor names sand values as a to a hash with names
# and values as key/value pairs...
#
sub GetDescriptorNamesAndValues {
  my($This) = @_;
  my(%NamesAndValues);

  %NamesAndValues = ();
  @NamesAndValues{ @{$This->{DescriptorNames}} } = $This->GetDescriptorValues();

  return %NamesAndValues;
}

# Return a string containing descriptor names and values...
#
sub _StringifyDescriptorNamesAndValues {
  my($This) = @_;
  my($NamesAndValuesString, $Name, $Value, @NamesAndValuesInfo, %NamesAndValues);

  @NamesAndValuesInfo = ();
  %NamesAndValues = $This->GetDescriptorNamesAndValues();

  for $Name (@{$This->{DescriptorNames}}) {
    $Value = $NamesAndValues{$Name};
    $Value = (TextUtil::IsEmpty($Value) || $Value =~ /^None$/i) ? 'None' : $Value;
    push @NamesAndValuesInfo, "$Name - $Value";
  }
  if (@NamesAndValuesInfo) {
    $NamesAndValuesString = "Names - Values: <" . TextUtil::JoinWords(\@NamesAndValuesInfo, ", ", 0) . ">";
  }
  else {
    $NamesAndValuesString = "Names - Values: < None>";
  }
  return $NamesAndValuesString;
}

1;

__END__

=head1 NAME

MolecularDescriptors - MolecularDescriptors class

=head1 SYNOPSIS

use MolecularDescriptors::MolecularDescriptors;

use MolecularDescriptors::MolecularDescriptors qw(:all);

=head1 DESCRIPTION

B<MolecularDescriptors> base class used to derive all other molecular descriptors classes provides the following methods:

new, AddDescriptorNames, AddDescriptorValues, GetDescriptorNames,
GetDescriptorNamesAndValues, GetDescriptorValueByName, GetDescriptorValues,
IsDescriptorsGenerationSuccessful, SetDescriptorNames, SetDescriptorValues,
SetMolecule, SetType

B<MolecularDescriptors> class is  derived from B<ObjectProperty> base class which provides methods not
explicitly defined in B<Fingerprints> or B<ObjectProperty> classes using Perl's AUTOLOAD functionality.
These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

=head2 METHODS

=over 4

=item B<new>

    $NewMolecularDescriptors = new MolecularDescriptors::
                               MolecularDescriptors(%NamesAndValues);

Using specified I<MolecularDescriptors> property names and values hash, B<new> method creates a new object
and returns a reference to newly created B<MolecularDescriptors> object. By default, following properties are
initialized:

    Molecule = '';
    Type = '';

=item B<AddDescriptorNames>

    $MolecularDescriptors->AddDescriptorNames(@Name);

Adds specified descriptor I<Names> to the list of available descriptor names and returns
I<MolecularDescriptors>.

=item B<AddDescriptorValues>

    $MolecularDescriptors->AddDescriptorValues(@Values);

Adds specified descriptor I<Values> to the list of calculated descriptor values and returns
I<MolecularDescriptors>.

=item B<GetDescriptorNames>

    @Names = $MolecularDescriptors->GetDescriptorNames();

Returns an array  containing all available descriptor names.

=item B<GetDescriptorNamesAndValues>

    %NamesAndValuesReturn = $MolecularDescriptors->
                              GetDescriptorNamesAndValues();

Returns a hash containing all available descriptor names and calculated values.

=item B<GetDescriptorValueByName>

    $Value = $MolecularDescriptors->
                              GetDescriptorValueByName($Name);

Returns calculated value for a specified descriptor name. A string B<None> is returned for
unknown descriptor names or for those descriptors whose values haven't been calculated.

=item B<GetDescriptorValues>

    @Values = $MolecularDescriptors->GetDescriptorValues();

Returns an array containing calculated descriptor values for all available descriptors.
Unless B<CalculateDescriptorsValues> method has been successfully invoked on a
I<MolecularDescriptors> object, value of each descriptor corresponds to string B<None>.

=item B<IsDescriptorsGenerationSuccessful>

    $Status = $MolecularDescriptors->
                              IsDescriptorsGenerationSuccessful();

Returns 1 or 0 based on whether molecular descriptors generation was successful.
For a successful molecular descriptors calculation, all available descriptors must have
a values other than a string I<None> which are set by B<CalculateDescriptorsValues>
method after successful completion of descriptors calculation.

=item B<SetDescriptorNames>

    $MolecularDescriptors->SetDescriptorNames(@Names);

Sets names of available descriptors to specified names and returns I<MolecularDescriptors>.

=item B<SetDescriptorValues>

    $MolecularDescriptors->SetDescriptorValues(@Values);

Sets values of available descriptors to specified values and returns I<MolecularDescriptors>.

=item B<SetMolecule>

    $MolecularDescriptors->SetMolecule($Molecule);

Sets molecule to use during calculation of molecular descriptors and returns I<MolecularDescriptors>.

=item B<SetType>

    $MolecularDescriptors->SetType($Type);

Sets I<Type> for I<MolecularDescriptors> object and returns I<MolecularDescriptors>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

MolecularDescriptorsGenerator.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
