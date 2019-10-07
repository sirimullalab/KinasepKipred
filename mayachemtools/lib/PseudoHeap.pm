package PseudoHeap;
#
# File: PseudoHeap.pm
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
use TextUtil ();

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK]
	       );

# Setup class variables...
my($ClassName);
_InitializeClass();

use overload '""' => 'StringifyPseudoHeap';

# PseudoHeap is designed to support tracking of a specific number of largest or smallest key/value
# pairs with numeric or alphanumeric keys along with corresponding scalar or reference values.
#
# Although PseudoHeap is similar to a heap, it lacks number of key properties of a traditional heap data
# structure: no concept of root, parent and child nodes; no ordering of keys in any particular order; no
# specific localtion greatest or smallest key.
#
# The keys are simply stored in a hash with each key poining to an array containing specified values.
# The min/max keys are updated during addition and deletion of key/value pairs; these can be retrieved
# by accessing corresponding hash.
#
# Addition and deletion of key/value is also straightforward using hashes. However, min/max keys
# need to be identified which is done using Perl sort on the keys.
#
#
# Class constructor...
#
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializePseudoHeap();

  $This->_InitializePseudoHeapProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializePseudoHeap {
  my($This) = @_;

  # Type of pseudo heap:
  #
  # KeepTopN - Keep track of a specified number largest of key/value pairs
  # KeepBottomN - Keep track of a specified number smallest of key/value pairs
  #
  $This->{Type} = undef;

  # Type of keys: Numeric or Alphanumeric
  #
  # The value of KeyType determines comparison function used to sort and
  # and compare keys for a specific heap type as shown below:
  #
  # Type             KeyType       Comp  Sorting
  #
  # KeepTopN      Numeric       <     Descending
  # KeepTopN      AlphaNumeric  lt    Descending
  # KeepBottomN  Numeric        >     Ascending
  # KeepBottomN  AlphaNumeric   gt    Ascending
  #
  $This->{KeyType} = undef;

  # Maximum number of largest or smallest key/value pairs to keep...
  #
  $This->{MaxSize} = 10;

  # Keys and values associated with each key as an array...
  %{$This->{Keys}} = ();

  # Min and max keys...
  $This->{MinKey} = undef;
  $This->{MaxKey} = undef;

  # Number of key/valur pairs currently present...
  $This->{CurrentSize} = 0;

  # Number of keys currently present where each key correspond to multiple values...
  $This->{KeysCount} = 0;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

}

# Initialize object properties....
#
sub _InitializePseudoHeapProperties {
  my($This, %NamesAndValues) = @_;
  my($Name, $Value, $MethodName);

  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  if (!exists $NamesAndValues{Type}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying Type...";
  }

  if (!exists $NamesAndValues{KeyType}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying KeyType...";
  }
}

# Set heap type...
#
sub SetType {
  my($This, $Type) = @_;

  if (defined $This->{Type}) {
    croak "Error: ${ClassName}->SetType: Can't change Type...";
  }

  if ($Type !~ /^(KeepTopN|KeepBottomN)$/i) {
    croak "Error: ${ClassName}->SetType: Unknown PseudoHeap type: $Type; Supported types: KeepTopN or KeepBottomN...";
  }
  $This->{Type} = $Type;

  return $This;
}

# Get heap type..
#
sub GetType {
  my($This) = @_;

  return defined $This->{Type} ? $This->{Type} : 'None';
}

# Set key type...
#
sub SetKeyType {
  my($This, $KeyType) = @_;

  if (defined $This->{KeyType}) {
    croak "Error: ${ClassName}->SetType: Can't change KeyType...";
  }

  if ($KeyType !~ /^(Numeric|Alphanumeric)$/i) {
    croak "Error: ${ClassName}->SetType: Unknown PseudoHeap key type: $KeyType; Supported key types: Numeric or Alphanumeric...";
  }
  $This->{KeyType} = $KeyType;

  return $This;
}

# Get key type..
#
sub GetKeyType {
  my($This) = @_;

  return defined $This->{KeyType} ? $This->{KeyType} : 'None';
}

# Add a key/value pair...
#
sub AddKeyValuePair {
  my($This, $Key, $Value) = @_;

  if (!(defined($Key) && defined($Value))) {
    carp "Warning: ${ClassName}->AddKeyValuePair: No key added: Both key and value must be defined...";
    return undef;
  }

  $This->_AddKeyValuePair($Key, $Value);

  return $This;
}

# Add multiple key/value pairs...
#
sub AddKeyValuePairs {
  my($This, @KeyValuePairs) = @_;

  if (!@KeyValuePairs) {
    carp "Warning: ${ClassName}->AddKeyValuePairs: No keys added: Key/Value pairs list is empty...";
    return undef;
  }
  if (@KeyValuePairs % 2) {
    carp "Warning: ${ClassName}->AddKeyValuePairs: No keys pairs added: Invalid key/value pairs data: Input list must contain even number of values...";
    return undef;
  }

  my($Key, $Value, $Index);
  for ($Index = 0; $Index < $#KeyValuePairs; $Index += 2) {
    $Key = $KeyValuePairs[$Index]; $Value = $KeyValuePairs[$Index + 1];
    $This->AddKeyValuePair($Key, $Value);
  }

  return $This;
}

# Delete specified keys along with all associated values for each key...
#
sub DeleteKeys {
  my($This, @Keys) = @_;

  if (!@Keys) {
    carp "Warning: ${ClassName}->DeleteKeys: No keys deleted: Keys list is empty...";
    return undef;
  }
  my($Key);
  for $Key (@Keys) {
    $This->DeleteKey($Key);
  }

  return $This;
}

# Delete a sepcified key along with all of its associated values...
#
sub DeleteKey {
  my($This, $Key) = @_;

  if (!defined $Key ) {
    carp "Warning: ${ClassName}->DeleteKey: No key deleted: Key must be specified...";
    return undef;
  }

  return $This->_DeleteKey($Key);
}

# Delete min key along with all of its associated values...
#
sub DeleteMinKey {
  my($This) = @_;

  return $This->DeleteKey($This->{MinKey});
}

# Delete max key along with all of its associated values...
#
sub DeleteMaxKey {
  my($This) = @_;

  return $This->DeleteKey($This->{MaxKey});
}

# Set max size...
#
sub SetMaxSize {
  my($This, $Size) = @_;

  if (!TextUtil::IsPositiveInteger($Size)) {
    croak "Error: ${ClassName}->SetMaxSize: Max size value, $Size, is not valid: It must be a positive  integer...";
  }

  if (defined($This->{MinKey}) || defined($This->{MaxKey})) {
    croak "Error: ${ClassName}->SetMaxSize: Can't change max size: Keys are already present...";
  }

  $This->{MaxSize} = $Size;

  return $This;
}

# Get max size...
#
sub GetMaxSize {
  my($This) = @_;

  return $This->{MaxMaxSize};
}

# Get current size...
#
sub GetCurrentSize {
  my($This) = @_;

  return $This->{CurrentSize};
}

# Get min key...
#
sub GetMinKey {
  my($This) = @_;

  return defined $This->{MinKey} ? $This->{MinKey} : 'None';
}

# Get max key...
#
sub GetMaxKey {
  my($This) = @_;

  return defined $This->{MaxKey} ? $This->{MaxKey} : 'None';
}

# Get keys...
#
sub GetKeys {
  my($This) = @_;

  return wantarray ? keys %{$This->{Keys}} : scalar keys %{$This->{Keys}};
}

# Get sorted keys...
#
sub GetSortedKeys {
  my($This) = @_;
  my(@SortedKeys);

  @SortedKeys = ();
  if ($This->{Type} =~ /^KeepTopN$/i) {
    @SortedKeys = ($This->{KeyType} =~ /^Numeric$/i) ? (sort { $b <=> $a } keys %{$This->{Keys}}) : (sort { $b cmp $a } keys %{$This->{Keys}});
  }
  elsif ($This->{Type} =~ /^KeepBottomN$/i) {
    @SortedKeys = ($This->{KeyType} =~ /^Numeric$/i) ? (sort { $a <=> $b } keys %{$This->{Keys}}) : (sort { $a cmp $b } keys %{$This->{Keys}});
  }

  return wantarray ? @SortedKeys : scalar @SortedKeys;
}

# Get values associated with a specified key...
sub GetKeyValues {
  my($This, $Key) = @_;
  my(@KeyValues);

  @KeyValues = ();
  if (defined($Key) && exists($This->{Keys}{$Key})) {
    @KeyValues = @{$This->{Keys}{$Key}};
  }
  return wantarray ? @KeyValues : scalar @KeyValues;
}

#  Add key/value pair...
#
sub _AddKeyValuePair{
  my($This, $Key, $Value) = @_;

  if ($This->{CurrentSize} < $This->{MaxSize}) {
    return $This->_AppendKeyValuePair($Key, $Value);
  }
  else {
    return $This->_InsertKeyValuePair($Key, $Value);
  }
}

# Append key/value pair...
#
sub _AppendKeyValuePair {
  my($This, $Key, $Value) = @_;

  if (!exists $This->{Keys}{$Key}) {
    @{$This->{Keys}{$Key}} = ();
    $This->{KeysCount} += 1;

    $This->_CompareAndSetMinKey($Key);
    $This->_CompareAndSetMaxKey($Key);
  }

  push @{$This->{Keys}{$Key}}, $Value;
  $This->{CurrentSize} += 1;

  return $This;
}

# Insert key/value pair...
#
sub _InsertKeyValuePair {
  my($This, $Key, $Value) = @_;

  # Is this key need to be inserted?
  if (!$This->_IsKeyNeedToBeInserted($Key)) {
    return $This;
  }

  # Insert key/value pair...
  if (!exists $This->{Keys}{$Key}) {
    @{$This->{Keys}{$Key}} = ();
    $This->{KeysCount} += 1;
  }
  push @{$This->{Keys}{$Key}}, $Value;
  $This->{CurrentSize} += 1;

  # Remove min or max key/value pair along with its update...
  my($KeyToDetele);

  $KeyToDetele = ($This->{Type} =~ /^KeepTopN$/i) ? $This->{MinKey} : $This->{MaxKey};
  $This->_DeleteKeyValuePair($KeyToDetele);

  return $This;
}

# Check whether it makes sense to insert specified key...
#
sub _IsKeyNeedToBeInserted {
  my($This, $Key) = @_;

  if ($This->{Type} =~ /^KeepTopN$/i) {
    if ($This->{KeyType} =~ /^Numeric$/i) {
      return ($Key < $This->{MinKey}) ? 0 : ((($This->{KeysCount} == 1) && ($This->{MinKey} == $Key)) ? 0 : 1);
    }
    else {
      return ($Key lt $This->{MinKey}) ? 0 : ((($This->{KeysCount} == 1) && ($This->{MinKey} eq $Key)) ? 0 : 1);
    }
  }
  elsif ($This->{Type} =~ /^KeepBottomN$/i) {
    if ($This->{KeyType} =~ /^Numeric$/i) {
      return ($Key > $This->{MaxKey}) ? 0 : ((($This->{KeysCount} == 1) && ($This->{MaxKey} == $Key)) ? 0 : 1);
    }
    else {
      return ($Key gt $This->{MaxKey}) ? 0 : ((($This->{KeysCount} == 1) && ($This->{MaxKey} eq $Key)) ? 0 : 1);
    }
  }

  return 1;
}

# Set min key...
#
sub _CompareAndSetMinKey {
  my($This, $Key) = @_;

  if (!defined $This->{MinKey}) {
    $This->{MinKey} = $Key;
    return $This;
  }

  if ($This->{KeyType} =~ /^Numeric$/i) {
    if ($Key < $This->{MinKey}) {
      $This->{MinKey} = $Key;
    }
  }
  else {
    if ($Key lt $This->{MinKey}) {
      $This->{MinKey} = $Key;
    }
  }

  return $This;
}

# Set max key...
#
sub _CompareAndSetMaxKey {
  my($This, $Key) = @_;

  if (!defined $This->{MaxKey}) {
    $This->{MaxKey} = $Key;
    return $This;
  }

  if ($This->{KeyType} =~ /^Numeric$/i) {
    if ($Key > $This->{MaxKey}) {
      $This->{MaxKey} = $Key;
    }
  }
  else {
    if ($Key gt $This->{MaxKey}) {
      $This->{MaxKey} = $Key;
    }
  }

  return $This;
}

# Delete a sepcified key along with all of its values added to the list...
#
sub _DeleteKey {
  my($This, $Key) = @_;
  my($NumOfValues);

  if (!exists $This->{Keys}{$Key}) {
    return undef;
  }

  # Delete all key values...
  $NumOfValues = scalar @{$This->{Keys}{$Key}};
  @{$This->{Keys}{$Key}} = ();
  $This->{CurrentSize} -= $NumOfValues;

  # Delete key...
  delete $This->{Keys}{$Key};
  $This->{KeysCount} -= 1;

  # Set min and max keys...
  $This->_FindAndSetMinAndMaxKeys();

  return $This;
}

# Delete a sepcified key along with its most recent value added to the list...
#
sub _DeleteKeyValuePair {
  my($This, $Key) = @_;

  if (!exists $This->{Keys}{$Key}) {
    return undef;
  }

  # Delete value...
  pop @{$This->{Keys}{$Key}};
  $This->{CurrentSize} -= 1;

  # Delete key...
  if (!@{$This->{Keys}{$Key}}) {
    delete $This->{Keys}{$Key};
    $This->{KeysCount} -= 1;
  }

  # Set min and max keys...
  $This->_FindAndSetMinAndMaxKeys();

  return $This;
}

# Set min and max key...
#
sub _FindAndSetMinAndMaxKeys {
  my($This) = @_;
  my(@SortedKeys);

  @SortedKeys = ($This->{KeyType} =~ /^Numeric$/i) ? (sort { $a <=> $b } keys %{$This->{Keys}}) : (sort { $a cmp $b } keys %{$This->{Keys}});

  if (@SortedKeys) {
    $This->{MinKey} = $SortedKeys[0];
    $This->{MaxKey} = $SortedKeys[$#SortedKeys];
  }
  else {
    $This->{MinKey} = undef;
    $This->{MaxKey} = undef;
  }

  return $This;
}

# Return a string containing vector values...
sub StringifyPseudoHeap {
  my($This) = @_;
  my($PseudoHeapString, $Key, $Value, $KeyValuesString, @KeysAndValues);

  $PseudoHeapString = "PseudoHeap: Type: " . $This->GetType() . "; KeyType: " . $This->GetKeyType() . "; MaxSize: $This->{MaxSize}; CurrentSize: $This->{CurrentSize}; MinKey: " . $This->GetMinKey() .  "; MaxKey: " . $This->GetMaxKey() . "; NumOfUniqueKeys: $This->{KeysCount}";

  @KeysAndValues = ();
  for $Key ($This->GetSortedKeys()) {
    for $Value ($This->GetKeyValues($Key)) {
      push @KeysAndValues, "$Key - $Value";
    }
  }
  if (@KeysAndValues) {
    $KeyValuesString = TextUtil::JoinWords(\@KeysAndValues, "; ", 0);
  }
  else {
    $KeyValuesString = "None";
  }

  $PseudoHeapString .= "; Sorted Key - Value pairs: [$KeyValuesString]";

  return $PseudoHeapString;
}

1;

__END__

=head1 NAME

PseudoHeap

=head1 SYNOPSIS

use PseudoHeap;

use PseudoHeap qw(:all);

=head1 DESCRIPTION

B<PseudoHeap> class provides the following methods:

new, AddKeyValuePair, AddKeyValuePairs, DeleteKey, DeleteKeys, DeleteMaxKey,
DeleteMinKey, GetCurrentSize, GetKeyType, GetKeyValues, GetKeys, GetMaxKey,
GetMaxSize, GetMinKey, GetSortedKeys, GetType, SetKeyType, SetMaxSize, SetType,
StringifyPseudoHeap

PseudoHeap is designed to support tracking of a specific number of largest or smallest key/value
pairs with numeric or alphanumeric keys along with corresponding scalar or reference values.

Although PseudoHeap is conceptually similar to a heap, it lacks number of key properties of a traditional
heap data structure: no concept of root, parent and child nodes; no ordering of keys in any particular
order; no specific location greatest or smallest key.

The keys are simply stored in a hash with each key pointing to an array containing specified values.
The min/max keys are updated during addition and deletion of key/value pairs; these can be retrieved
by accessing corresponding hash.

Addition and deletion of key/value is also straightforward using hashes. However, min/max keys
need to be identified which is done using Perl sort function on the keys.

=head2 FUNCTIONS

=over 4

=item B<new>

    $NewPseudoHeap = new PseudoHeap(%NamesAndValues);

Using specified parameters I<NamesAndValues> names and values hash, B<new> method creates
a new object and returns a reference to a newly created B<NewPseudoHeap> object. By default,
the following property names are initialized:

    Type = undef;
    KeyType = undef;
    MaxSize = 10;

Examples:

    $NewPseudoHeap = new PseudoHeap(
                               'Type' => 'KeepTopN',
                               'KeyType' => 'Numeric');

    $NewPseudoHeap = new PseudoHeap(
                               'Type' => 'KeepTopN',
                               'KeyType' => 'AlphaNumeric',
                               'MaxSize' => '20');

    $NewPseudoHeap = new PseudoHeap(
                               'Type' => 'KeepBottomN',
                               'KeyType' => 'AlphaNumeric',
                               'MaxSize' => '20');

=item B<AddKeyValuePair>

    $PseudoHeap->AddKeyValuePair($Key, $Value);

Add specified I<Key> and I<Value> pair to pseudo heap using a new or an existing
key and returns B<PseudoHeap>.

=item B<AddKeyValuePairs>

    $PseudoHeap->AddKeyValuePairs(@KeyValuePairs);

Adds multiple key and value pairs specified in array I<KeyValuePairs> to pseudo heap
using a new or existing keys and returns B<PseudoHeap>.

=item B<DeleteKey>

    $PseudoHeap->DeleteKey($Key);

Deletes a specified I<Key> from pseudo heap and returns B<PseudoHeap>.

=item B<DeleteKeys>

    $PseudoHeap->DeleteKeys(@Keys);

Deletes a specified I<Keys> from pseudo heap and returns B<PseudoHeap>.

=item B<DeleteMaxKey>

    $PseudoHeap->DeleteMaxKey();

Deletes a I<MaxKey> along with its associated values from pseudo heap and returns
B<PseudoHeap>.

=item B<DeleteMinKey>

    $PseudoHeap->DeleteMinKey();

Deletes a I<MinKey> along with its associated values from pseudo heap and returns
B<PseudoHeap>.

=item B<GetCurrentSize>

    $Size = $PseudoHeap->GetCurrentSize();

Returns current I<Size> of pseudo heap corresponding to number to keys in heap.

=item B<GetKeyType>

    $KeyType = $PseudoHeap->GetKeyType();

Returns I<KeyType> of pseudo heap. Possible B<KeyType> values: I<Numeric or Alphanumeric>.

=item B<GetKeyValues>

    @Values = $PseudoHeap->GetKeyValues($Key);
    $NumOfValues = $PseudoHeap->GetKeyValues($Key);

Returns an array containing B<Values> associated with a specified I<Key> in pseudo heap. In
scalar context, it returns number of values associated with a key.

=item B<GetKeys>

    @Keys = $PseudoHeap->GetKeys();
    $NumOfKeys = $PseudoHeap->GetKeys();

Returns an array containing all B<Keys> in pseudo heap. In scalar context, it returns total
number of keys.

=item B<GetMaxKey>

    $MaxKey = $PseudoHeap->GetMaxKey();

Returns I<MaxKey> present in pseudo heap.

=item B<GetMaxSize>

    $MaxSize = $PseudoHeap->GetMaxSize();

Returns I<MaxSize> of pseudo heap.

=item B<GetMinKey>

    $MinKey = $PseudoHeap->GetMinKey();

Returns I<MinKey> present in pseudo heap.

=item B<GetSortedKeys>

    @Keys = $PseudoHeap->GetSortedKeys();
    $NumOfKeys = $PseudoHeap->GetSortedKeys();

Returns an array containing all sorted B<Keys> in pseudo heap. In scalar context, it retruns
total number of keys.

Keys are sorted based on values of B<Type> and B<KeyType> for pseudo heap:

    Type          KeyType       SortOrder   SortOperator
    KeepTopN      Numeric       Descending  <=>
    KeepTopN      Alphanumeric  Descending  cmp
    KeepBottomN   Numeric       Ascending    <=>
    KeepBottomN   Alphanumeric  Ascending   cmp

=item B<GetType>

    $Type = $PseudoHeap->GetType();

Returns I<Type> of pseudo heap.

=item B<SetKeyType>

    $PseudoHeap->SetKeyType($KeyType);

Sets I<KeyType> of pseudo heap and returns B<PseudoHeap>.

=item B<SetMaxSize>

    $PseudoHeap->SetMaxSize($MaxSize);

Sets I<MaxSize> of pseudo heap and returns B<PseudoHeap>.

=item B<SetType>

    $PseudoHeap->SetType($Type);

Sets I<Type> of pseudo heap and returns B<PseudoHeap>.

=item B<StringifyPseudoHeap>

    $PseudoHeapString = $PseudoHeap->StringifyPseudoHeap();

Returns a string containing information about I<PseudoHeap> object

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
