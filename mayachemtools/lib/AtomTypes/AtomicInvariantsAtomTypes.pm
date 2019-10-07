package AtomTypes::AtomicInvariantsAtomTypes;
#
# File: AtomicInvariantsAtomTypes.pm
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
use AtomTypes::AtomTypes;
use Molecule;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(AtomTypes::AtomTypes Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(IsAtomicInvariantAvailable GetAvailableAtomicInvariants);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, @AtomicInvariantsOrder, %AvailableAtomicInvariants, %AvailableAtomicInvariantsByDescription);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyAtomicInvariantsAtomTypes';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeAtomicInvariantsAtomTypes();

  $This->_InitializeAtomicInvariantsAtomTypesProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Initialize class atomic invariants...
  _InitializeClassAtomicInvariants();
}

# Initialize class level atomic invariants information which doesn't change during
# instantiations of objects...
#
sub _InitializeClassAtomicInvariants {
  # Available atomic invariants for generating atom types...
  #
  %AvailableAtomicInvariants = ();
  %AvailableAtomicInvariants = ('AS' => 'AtomSymbol|ElementSymbol',
				'X' => 'NumOfNonHydrogenAtomNeighbors|NumOfHeavyAtomNeighbors',
				'BO' => 'SumOfBondOrdersToNonHydrogenAtoms|SumOfBondOrdersToHeavyAtoms',
				'LBO' => 'LargestBondOrderToNonHydrogenAtoms|LargestBondOrderToHeavyAtoms',
				'SB' => 'NumOfSingleBondsToNonHydrogenAtoms|NumOfSingleBondsToHeavyAtoms',
				'DB' => 'NumOfDoubleBondsToNonHydrogenAtoms|NumOfDoubleBondsToHeavyAtoms',
				'TB' => 'NumOfTripleBondsToNonHydrogenAtoms|NumOfTripleBondsToHeavyAtoms',
				'H' => 'NumOfImplicitAndExplicitHydrogens',
				'Ar' => 'Aromatic',
				'RA' => 'RingAtom',
				'FC' => 'FormalCharge',
				'MN' => 'MassNumber',
				'SM' => 'SpinMultiplicity');

  # Setup available atomic invariants description to abbreviation map...
  #
  my($Key, $Value, $Description, @Descriptions);
  %AvailableAtomicInvariantsByDescription = ();
  while (($Key, $Value) = each %AvailableAtomicInvariants) {
    @Descriptions = ($Value =~ /|/) ? (split /\|/, $Value) : ($Value);
    for $Description (@Descriptions) {
      $AvailableAtomicInvariantsByDescription{$Description} = $Key;
    }
  }

  # Atomic invariants order used for generating atom types...
  #
  @AtomicInvariantsOrder = ();
  @AtomicInvariantsOrder = ('AS', 'X', 'BO', 'LBO', 'SB', 'DB', 'TB', 'H', 'Ar', 'RA', 'FC', 'MN', 'SM');
}

# Initialize object data...
#
sub _InitializeAtomicInvariantsAtomTypes {
  my($This) = @_;

  # Type of AtomTypes...
  $This->{Type} = 'AtomicInvariants';

  # By default hydrogens are also assigned atom types...
  $This->{IgnoreHydrogens} = 0;

  # Initialize atom types information...
  $This->_InitializeAtomTypesInformation();

  return $This;
}

# Inialize atomic invariants information used for generating atom types...
#
sub _InitializeAtomTypesInformation {
  my($This) = @_;

  # Default atomic invariants to use for generating atom types: AS, X, BO, H, FC
  #
  %{$This->{AtomicInvariantsToUse}} = ();
  %{$This->{AtomicInvariantsToUse}} = ('AS' => 1, 'X' => 1, 'BO' => 1, 'LBO' => 0,
				       'SB' => 0, 'DB' => 0, 'TB' => 0,
				       'H' => 1, 'Ar' => 0, 'RA' => 0, 'FC' => 1, 'MN' => 0, 'SM' => 0);

  return $This;
}

# Initialize object properties...
#
sub _InitializeAtomicInvariantsAtomTypesProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  # Make sure molecule object was specified...
  if (!exists $NamesAndValues{Molecule}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying molecule...";
  }

  return $This;
}

# Disable change of AvailableAtomicInvariants...
#
sub SetAvailableAtomicInvariants {
  my($This) = @_;

  carp "Warning: ${ClassName}->SetAtomicInvariantsOrder: Available atomic invariants can't be changed...";

  return $This;
}

# Disable change of atomic invariants order used for generation of atom types...
#
sub SetAtomicInvariantsOrder {
  my($This) = @_;

  carp "Warning: ${ClassName}->SetAtomicInvariantsOrder: Atomic invariants order can't be changed...";

  return $This;
}

# Set atom invariants to use for atom types...
#
sub SetAtomicInvariantsToUse {
  my($This, @Values) = @_;
  my($FirstValue, $TypeOfFirstValue, $AtomicInvariant, $SpecifiedAtomicInvariant, $AtomicInvariantValue, @SpecifiedAtomicInvariants, %AtomicInvariantsToUse);

  if (!@Values) {
    carp "Warning: ${ClassName}->SetAtomicInvariantsToUse: No values specified...";
    return;
  }

  $FirstValue = $Values[0];
  $TypeOfFirstValue = ref $FirstValue;
  @SpecifiedAtomicInvariants = ();

  if ($TypeOfFirstValue =~ /^ARRAY/) {
    push @SpecifiedAtomicInvariants, @{$FirstValue};
  }
  else {
    push @SpecifiedAtomicInvariants, @Values;
  }

  # Make sure specified AtomicInvariants are valid...
  for $SpecifiedAtomicInvariant (@SpecifiedAtomicInvariants) {
    if (exists $AvailableAtomicInvariants{$SpecifiedAtomicInvariant}) {
      $AtomicInvariant = $SpecifiedAtomicInvariant;
    }
    elsif ($AvailableAtomicInvariantsByDescription{$SpecifiedAtomicInvariant}) {
      $AtomicInvariant = $AvailableAtomicInvariantsByDescription{$SpecifiedAtomicInvariant};
    }
    else {
      croak "Error: ${ClassName}->SetAtomicInvariantsToUse: Specified atomic invariant, $SpecifiedAtomicInvariant, is not supported...\n ";
    }
    $AtomicInvariantsToUse{$AtomicInvariant} = 1;
  }

  # Make sure AtomSymbol is always used...
  if (!(exists($AtomicInvariantsToUse{AS}) && $AtomicInvariantsToUse{AS} == 1)) {
    croak "Error: ${ClassName}->SetAtomicInvariantsToUse: AtomicInvariant AtomSymbol must be specified...\n ";
  }

  # Set atomic invariants...
  for $AtomicInvariant (keys %{$This->{AtomicInvariantsToUse}}) {
    $This->{AtomicInvariantsToUse}{$AtomicInvariant} = 0;
    if (exists $AtomicInvariantsToUse{$AtomicInvariant}) {
      $This->{AtomicInvariantsToUse}{$AtomicInvariant} = 1;
    }
  }

  return $This;
}

# Is it an available AtomicInvariant?
#
sub IsAtomicInvariantAvailable {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $AtomicInvariant, $Status);

  if ((@_ == 2) && (_IsAtomicInvariantsAtomTypes($FirstParameter))) {
    ($This, $AtomicInvariant) = ($FirstParameter, $SecondParameter);
  }
  else {
    $AtomicInvariant = $FirstParameter;
  }
  $Status = exists($AvailableAtomicInvariants{$AtomicInvariant}) || exists($AvailableAtomicInvariantsByDescription{$AtomicInvariant}) ? 1 : 0;

  return $Status;
}

# Get a hash containing available atomic invariants and their description
# as key/value pairs.
#
sub GetAvailableAtomicInvariants {
  return %AvailableAtomicInvariants;
}

# Get an array containing order of atomic invariants used to generate atom types...
#
sub GetAtomicInvariantsOrder {
  return @AtomicInvariantsOrder;
}

# Assign atom types to all atoms...
#
# Let:
#   AS = Atom symbol corresponding to element symbol
#
#   X<n>   = Number of non-hydrogen atom neighbors or heavy atoms attached to atom
#   BO<n> = Sum of bond orders to non-hydrogen atom neighbors or heavy atoms attached to atom
#   LBO<n> = Largest bond order of non-hydrogen atom neighbors or heavy atoms attached to atom
#   SB<n> = Number of single bonds to non-hydrogen atom neighbors or heavy atoms attached to atom
#   DB<n> = Number of double bonds to non-hydrogen atom neighbors or heavy atoms attached to atom
#   TB<n> = Number of triple bonds to non-hydrogen atom neighbors or heavy atoms attached to atom
#   H<n>   = Number of implicit and explicit hydrogens for atom
#   Ar     = Aromatic annotation indicating whether atom is aromatic
#   RA     = Ring atom annotation indicating whether atom is a ring
#   FC<+n/-n> = Formal charge assigned to atom
#   MN<n> = Mass number indicating isotope other than most abundant isotope
#   SM<n> = Spin multiplicity of atom. Possible values: 1 (singlet), 2 (doublet) or 3 (triplet)
#
# Then:
#
#   AtomType specification corresponds to:
#
#     AS.X<n>.BO<n>.LBO<n>.<SB><n>.<DB><n>.<TB><n>.H<n>.Ar.RA.FC<+n/-n>.MN<n>.SM<n>
#
# Except for AS which is a required atomic invariant in atom types, all other atomic invariants are
# optional. Default atomic invariants used for AtomID are: AS, X<n>, BO<n>, H<n>, FC<+n/-n>.
# AtomID specification doesn't include atomic invariants with zero or undefined values.
#
# Notes:
#   . AtomicInvariants with zero or undefined values are not shown.
#   . LBO with value of 1 is not shown. And absence of LBO in AtomTypes implies the largest
#     bond order value is one.
#   . SB, DB and TB with values of zero are not shown.
#   . The difference in BO and X values corresponds to numbed of pi electrons [ Ref 57 ].
#
# Examples of atomic invariant atom types:
#
#   O.X1.BO1.H1 - Hydroxyl oxygen in carboxylate with attached hydrogen and no explicit charge
#   O.X1.BO1.FC-1 - Hydroxyl ozygen in carboxylate with explicit negative charge
#   O.X1.BO2 - Carbonyl oxygen in carboxylate with double bond to carbon
#   O.X2.BO2 - Hydroxyl ozygen in carboxylate attached to carbonyl carbon and another heavy atom
#
#   C.X2.BO3.H1.Ar - Aromatic carbon
#
sub AssignAtomTypes {
  my($This) = @_;
  my($Atom, $AtomType, $AtomicInvariant, $AtomicInvariantValue, @AtomicInvariants);

  ATOM: for $Atom ($This->GetMolecule()->GetAtoms()) {
    if ($This->{IgnoreHydrogens} && $Atom->IsHydrogen()) {
      next ATOM;
    }
    @AtomicInvariants = ();

    # Go over atomic invariants...
    ATOMICINVARIANT: for $AtomicInvariant (@AtomicInvariantsOrder) {
      if (!$This->{AtomicInvariantsToUse}{$AtomicInvariant}) {
	next ATOMICINVARIANT;
      }
      $AtomicInvariantValue = $Atom->GetAtomicInvariantValue($AtomicInvariant);
      if (!(defined($AtomicInvariantValue) && $AtomicInvariantValue)) {
	next ATOMICINVARIANT;
      }
      if ($AtomicInvariant =~ /^AS$/i) {
	push @AtomicInvariants, $AtomicInvariantValue;
      }
      elsif ($AtomicInvariant =~ /^Ar$/i) {
	push @AtomicInvariants, "Ar";
      }
      elsif ($AtomicInvariant =~ /^RA$/i) {
	push @AtomicInvariants, "RA";
      }
      elsif ($AtomicInvariant =~ /^FC$/i) {
	push @AtomicInvariants, ($AtomicInvariantValue > 0) ? "FC+${AtomicInvariantValue}" : "FC${AtomicInvariantValue}";
      }
      elsif ($AtomicInvariant =~ /^LBO$/i) {
	if ($AtomicInvariantValue > 1) {
	  push @AtomicInvariants, "${AtomicInvariant}${AtomicInvariantValue}";
	}
      }
      else {
	push @AtomicInvariants, "${AtomicInvariant}${AtomicInvariantValue}";
      }
    }
    # Create and assign atom type to atom...
    $AtomType = TextUtil::JoinWords(\@AtomicInvariants, ".", 0);
    $This->SetAtomType($Atom, $AtomType);
  }
  return $This;
}

# Are all atoms types successfully assigned?
#
# Notes:
#   . Base class method is overridden to always return 1: An appropriate value, atomic invariant
#     atom types delimited by dot, is always assigned to atoms.
#
sub IsAtomTypesAssignmentSuccessful {
  my($This) = @_;

  return 1;
}

# Return a string containg data for AtomicInvariantsAtomTypes object...
#
sub StringifyAtomicInvariantsAtomTypes {
  my($This) = @_;
  my($AtomTypesString);

  # Type of AtomTypes...
  $AtomTypesString = "AtomTypes: $This->{Type}; IgnoreHydrogens: " . ($This->{IgnoreHydrogens} ? "Yes" : "No");

  # AvailableAtomicInvariants and AtomicInvariantsToUse...
  my($AtomicInvariant, @AtomicInvariants, @AtomicInvariantsToUse);

  @AtomicInvariantsToUse = ();
  @AtomicInvariants = ();
  for $AtomicInvariant (@AtomicInvariantsOrder) {
    push @AtomicInvariants, "$AtomicInvariant: $AvailableAtomicInvariants{$AtomicInvariant}";
    if ($This->{AtomicInvariantsToUse}{$AtomicInvariant}) {
      push @AtomicInvariantsToUse, $AtomicInvariant;
    }
  }
  $AtomTypesString .= "; AtomicInvariantsToUse: <" . TextUtil::JoinWords(\@AtomicInvariantsToUse, ", ", 0) . ">";
  $AtomTypesString .= "; AtomicInvariantsOrder: <" . TextUtil::JoinWords(\@AtomicInvariantsOrder, ", ", 0) . ">";
  $AtomTypesString .= "; AvailableAtomicInvariants: <" . TextUtil::JoinWords(\@AtomicInvariants, ", ", 0) . ">";

  # Setup atom types information...
  my($AtomID, $AtomType, @AtomTypesInfo, %AssignedAtomTypes);

  @AtomTypesInfo = ();
  %AssignedAtomTypes = $This->GetAtomTypes();

  for $AtomID (sort { $a <=> $b } keys %AssignedAtomTypes) {
    $AtomType = $AssignedAtomTypes{$AtomID} ? $AssignedAtomTypes{$AtomID} : 'None';
    push @AtomTypesInfo, "$AtomID:$AtomType";
  }
  $AtomTypesString .= "; AtomIDs:AtomTypes: <" . TextUtil::JoinWords(\@AtomTypesInfo, ", ", 0) . ">";

  return $AtomTypesString;
}

# Is it a AtomicInvariantsAtomTypes object?
sub _IsAtomicInvariantsAtomTypes {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

AtomicInvariantsAtomTypes

=head1 SYNOPSIS

use AtomTypes::AtomicInvariantsAtomTypes;

use AtomTypes::AtomicInvariantsAtomTypes qw(:all);

=head1 DESCRIPTION

B<AtomicInvariantsAtomTypes> class provides the following methods:

new, AssignAtomTypes, GetAtomicInvariantsOrder, GetAvailableAtomicInvariants,
IsAtomicInvariantAvailable, SetAtomicInvariantsToUse, StringifyAtomicInvariantsAtomTypes

The following functions are available:

GetAvailableAtomicInvariants, IsAtomicInvariantAvailable

B<AtomicInvariantsAtomTypes> is derived from B<AtomTypes> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<AtomicInvariantsAtomTypes>, B<AtomTypes> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

Possible values for atomic invariants are: I<AS, X, BO,  LBO, SB, DB, TB,
H, Ar, RA, FC, MN, SM>. Default atom invariants values: I<AS,X,BO,H,FC>.

The atomic invariants abbreviations correspond to:

    AS = Atom symbol corresponding to element symbol

    X<n>   = Number of non-hydrogen atom neighbors or heavy atoms
    BO<n> = Sum of bond orders to non-hydrogen atom neighbors or heavy atoms
    LBO<n> = Largest bond order of non-hydrogen atom neighbors or heavy atoms
    SB<n> = Number of single bonds to non-hydrogen atom neighbors or heavy atoms
    DB<n> = Number of double bonds to non-hydrogen atom neighbors or heavy atoms
    TB<n> = Number of triple bonds to non-hydrogen atom neighbors or heavy atoms
    H<n>   = Number of implicit and explicit hydrogens for atom
    Ar     = Aromatic annotation indicating whether atom is aromatic
    RA     = Ring atom annotation indicating whether atom is a ring
    FC<+n/-n> = Formal charge assigned to atom
    MN<n> = Mass number indicating isotope other than most abundant isotope
    SM<n> = Spin multiplicity of atom. Possible values: 1 (singlet), 2 (doublet) or
            3 (triplet)

Atom type generated by AtomTypes::AtomTypes::AtomicInvariantsAtomTypes class corresponds to:

    AS.X<n>.BO<n>.LBO<n>.<SB><n>.<DB><n>.<TB><n>.H<n>.Ar.RA.FC<+n/-n>.MN<n>.SM<n>

Except for AS which is a required atomic invariant in atom types, all other atomic invariants are
optional. Atom type specification doesn't include atomic invariants with zero or undefined values.

In addition to usage of abbreviations for specifying atomic invariants, the following descriptive words
are also allowed:

    X : NumOfNonHydrogenAtomNeighbors or NumOfHeavyAtomNeighbors
    BO : SumOfBondOrdersToNonHydrogenAtoms or SumOfBondOrdersToHeavyAtoms
    LBO : LargestBondOrderToNonHydrogenAtoms or LargestBondOrderToHeavyAtoms
    SB :  NumOfSingleBondsToNonHydrogenAtoms or NumOfSingleBondsToHeavyAtoms
    DB : NumOfDoubleBondsToNonHydrogenAtoms or NumOfDoubleBondsToHeavyAtoms
    TB : NumOfTripleBondsToNonHydrogenAtoms or NumOfTripleBondsToHeavyAtoms
    H :  NumOfImplicitAndExplicitHydrogens
    Ar : Aromatic
    RA : RingAtom
    FC : FormalCharge
    MN : MassNumber
    SM : SpinMultiplicity

 Notes:

    . AtomicInvariants with zero or undefined values are not shown.
    . LBO with value of 1 is not shown. And absence of LBO in AtomTypes
      implies the largest bond order value is one.
    . SB, DB and TB with values of zero are not shown.
    . The difference in BO and X values corresponds to numbed of pi electrons [ Ref 57 ].

Examples of atomic invariant atom types:

    . O.X1.BO1.H1 - Hydroxyl oxygen in carboxylate with attached hydrogen
      and no explicit charge
    . O.X1.BO1.FC-1 - Hydroxyl ozygen in carboxylate with explicit negative
      charge
    . O.X1.BO2 - Carbonyl oxygen in carboxylate with double bond to carbon
    . O.X2.BO2 - Hydroxyl ozygen in carboxylate attached to carbonyl carbon
      and another heavy atom
    . C.X2.BO3.H1.Ar - Aromatic carbon

=head2 METHODS

=over 4

=item B<new>

    $NewAtomicInvariantsAtomTypes = new AtomTypes::AtomicInvariantsAtomTypes(
                                                   %NamesAndValues);

Using specified I<AtomicInvariantsAtomTypes> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<AtomicInvariantsAtomTypes>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'AtomicInvariants'
    IgnoreHydrogens = 0
    AtomicInvariantsToUse = AS,X,BO,H,FC

Examples:

    $AtomicInvariantsAtomTypes = new AtomTypes::AtomicInvariantsAtomTypes(
                              'Molecule' => $Molecule,
                              'IgnoreHydrogens' => 0,
                              'AtomicInvariantsToUse' =>
                                         ['AS', 'X', 'BO', 'H', 'FC']);

=item B<AssignAtomTypes>

    $AtomicInvariantsAtomTypes->AssignAtomTypes();

Assigns atomic invariant atom types to all the atoms in a molecule and returns
I<AtomicInvariantsAtomTypes>.

=item B<GetAtomicInvariantsOrder>

    @AtomicInvariantsOrder = $AtomicInvariantsAtomTypes->
                             GetAtomicInvariantsOrder();

Returns an array obtaining order of atomic invariants used to generate atom types.

=item B<GetAvailableAtomicInvariants>

    %AvailableAtomicInvariants = $AtomicInvariantsAtomTypes->
                                 GetAvailableAtomicInvariants();

Returns available atomic invariants as a hash containing available atomic invariants
and their description as key/value pairs.

=item B<IsAtomTypesAssignmentSuccessful>

    $Status = $AtomTypes->IsAtomTypesAssignmentSuccessful();

Returns 1 or 0 based on whether atom types assignment was successfully performed.
This method overrides the same method available in the base class AtomTypes.pm used
to derived this class.

=item B<IsAtomicInvariantAvailable>

    $Status = $AtomicInvariantsAtomTypes->
              IsAtomicInvariantAvailable($AtomicInvariant);
    $Status = AtomTypes::AtomicInvariantsAtomTypes::
              IsAtomicInvariantAvailable($AtomicInvariant);

Returns 1 or 0 based on whether I<AtomicInvariant> is valid.

=item B<SetAtomicInvariantsToUse>

    $AtomicInvariantsAtomTypes->SetAtomicInvariantsToUse($ValuesRef);
    $AtomicInvariantsAtomTypes->SetAtomicInvariantsToUse(@Values);

Sets atomic invariants to use for generating and assigning atom types and returns
I<AtomicInvariantsAtomTypes>.

=item B<StringifyAtomicInvariantsAtomTypes>

    $String = $AtomicInvariantsAtomTypes->StringifyAtomicInvariantsAtomTypes();

Returns a string containing information about I<AtomicInvariantsAtomTypes> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AtomTypes.pm, DREIDINGAtomTypes.pm, EStateAtomTypes.pm,
FunctionalClassAtomTypes.pm, MMFF94AtomTypes.pm, SLogPAtomTypes.pm,
SYBYLAtomTypes.pm, TPSAAtomTypes.pm, UFFAtomTypes.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
