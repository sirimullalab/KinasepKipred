package ObjectProperty;
#
# File: ObjectProperty.pm
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

use vars qw($AUTOLOAD);

# Set property for an object...
sub SetProperty {
  my($This, $Name, $Value) = @_;

  if (!(defined($Name) && defined($Value))) {
    return undef;
  }
  return $This->_SetProperty($Name, $Value);
}

# Set properties for an object...
sub SetProperties {
  my($This, %NamesAndValues) = @_;
  my($Name, $Value);

  while (($Name, $Value) = each  %NamesAndValues) {
    $This->_SetProperty($Name, $Value);
  }

  return $This;
}

# Set object property...
sub _SetProperty {
  my($This, $Name, $Value) = @_;

  $This->{$Name} = $Value;
}

# Get property for an object...
sub GetProperty {
  my($This, $Name) = @_;

  if (!defined $Name) {
    return undef;
  }
  return $This->_GetProperty($Name);
}

# Get object property...
sub _GetProperty {
  my($This, $Name) = @_;

  if (exists $This->{$Name}) {
    return $This->{$Name};
  }
  else {
    return undef;
  }
}

# Does this property exist?
sub HasProperty {
  my($This, $Name) = @_;

  if (!defined $Name) {
    return 0;
  }
  return (exists $This->{$Name}) ? 1 : 0;
}

# Delete object property...
sub DeleteProperty {
  my($This, $Name) = @_;

  if (!defined $Name) {
    return undef;
  }
  return $This->_DeleteProperty($Name);
}

# Delete object property...
sub _DeleteProperty {
  my($This, $Name) = @_;

  if (exists $This->{$Name}) {
    delete $This->{$Name};
  }
  return $This;
}

# Implements Set<PropertyName> and Get<PropertyName> methods...
sub AUTOLOAD {
  my($This, $PropertyValue) = @_;
  my($PackageName, $MethodName, $PropertyName, $ThisType);

  # Do a greedy match to make sure package name and method names are
  # picked up correctly from invocation names containing multiple occurences
  # of ::. For example: FileIO::SDFileIO::GetFileHandle and so on.
  #
  ($PackageName, $MethodName) = $AUTOLOAD =~ /^(.*)::(.*)$/;

  if ($MethodName =~ /^(BEGIN|DESTROY)$/) {
    return;
  }

  $ThisType = ref($This) or croak "Error: Invocation of function ${PackageName}::${MethodName} invocation is not supported: It must be invoked using an object reference...";

  if (!($MethodName =~ /^Get/ || $MethodName =~ /^Set/ || $MethodName =~ /^Delete/)) {
    croak "Error: Can't locate object method \"$MethodName\" via package \"$ThisType\": This method is not automatically implemented by AUTOLOAD: Only Get<PropertyName>, Set<PropertyName> and Delete<PropertyName> functions are implemented via AUTOLOAD...";
  }
  if ($MethodName =~ /^Delete/) {
    ($PropertyName) = $MethodName =~ /^Delete(.*?)$/;
  }
  else {
    ($PropertyName) = $MethodName =~ /^[SG]et(.*?)$/;
  }
  if ($MethodName =~ /^Set/ && !defined($PropertyValue)) {
    carp "Warning:  ${PackageName}::${MethodName}: Didn't set value for property $PropertyName: Property value for must be specified...\n";
    return undef;
  }

  if ($MethodName =~ /^Get/) {
    return $This->_GetProperty($PropertyName);
  }
  elsif ($MethodName =~ /^Set/) {
    return $This->_SetProperty($PropertyName, $PropertyValue);
  }
  elsif ($MethodName =~ /^Delete/) {
    return $This->_DeleteProperty($PropertyName);
  }

}

1;

__END__

=head1 NAME

ObjectProperty

=head1 SYNOPSIS

use ObjectProperty;

=head1 DESCRIPTION

B<ObjectProperty> is an abstract base class which implements methods not explicitly defined
in classed derived from this class using Perl's AUTOLOAD functionality. These methods are generated
on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

This class uses its parent class hash to set, get, and delete  propery names and values.

ObjectProperty module provides the following methods to be used in context of its parent class:

DeleteProperty, GetProperty, HasProperty, SetProperties, SetProperty

=head2 METHODS

=over 4

=item B<DeleteProperty>

    DeleteProperty($Name);

Deletes specified property I<Name>

=item B<GetProperty>

    GetProperty($Name);

Returns value associated with specified property I<Name>.

=item B<HasProperty>

    HasProperty($Name);

Returns 1 or 0 based on whether specified property I<Name> associated with an object.

=item B<SetProperties>

    SetProperties(%NamesAndValues);

Using specified property name and value hash I<NamesAndValues>, associates each
property I<Name> and I<Values> to an object.

=item B<SetProperty>

    SetProperty($Name, $Value);

Associate property I<Name> and I<Value> to an object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
