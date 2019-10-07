package AtomTypes::FunctionalClassAtomTypes;
#
# File: FunctionalClassAtomTypes.pm
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
@EXPORT_OK = qw(IsFunctionalClassAvailable GetAvailableFunctionalClasses GetFunctionalClassesOrder);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, @FunctionalClassesOrder, %AvailableFunctionalClasses, %AvailableFunctionalClassesByDescription);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyFunctionalClassAtomTypes';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeFunctionalClassAtomTypes();

  $This->_InitializeFunctionalClassAtomTypesProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Initialize class level functional classes information...
  _InitializeFunctionalClasses();
}

# Initialize class level functional class information which doesn't change during
# instantiations of objects...
#
sub _InitializeFunctionalClasses {
  # Available functional classes for generating atom types...
  #
  %AvailableFunctionalClasses = ();
  %AvailableFunctionalClasses = ('HBD' => 'HydrogenBondDonor',
				     'HBA' => 'HydrogenBondAcceptor',
				     'PI' => 'PositivelyIonizable',
				     'NI' => 'NegativelyIonizable',
				     'Ar' => 'Aromatic',
				     'Hal' => 'Halogen',
				     'H' => 'Hydrophobic',
				     'RA' => 'RingAtom',
				     'CA' => 'ChainAtom');

  # Setup available functional classe description to abbreviation map...
  #
  my($Key, $Value, $Description, @Descriptions);
  %AvailableFunctionalClassesByDescription = ();
  while (($Key, $Value) = each %AvailableFunctionalClasses) {
    @Descriptions = ($Value =~ /|/) ? (split /\|/, $Value) : ($Value);
    for $Description (@Descriptions) {
      $AvailableFunctionalClassesByDescription{$Description} = $Key;
    }
  }

  # Functional classes order used for generating atom types...
  #
  @FunctionalClassesOrder = ();
  @FunctionalClassesOrder = sort keys %AvailableFunctionalClasses;
}

# Initialize object data...
#
sub _InitializeFunctionalClassAtomTypes {
  my($This) = @_;

  # Type of AtomTypes...
  $This->{Type} = 'FunctionalClass';

  # By default hydrogens are also assigned atom types...
  $This->{IgnoreHydrogens} = 0;

  # Initialize atom types information...
  $This->_InitializeAtomTypesInformation();

  return $This;
}

# Inialize functional class information used for generating atom types...
#
sub _InitializeAtomTypesInformation {
  my($This) = @_;

  # Default functional classes to use for generating atom types: HBD, HBA, PI, NI, Ar, Hal
  #
  %{$This->{FunctionalClassesToUse}} = ();
  %{$This->{FunctionalClassesToUse}} = ('HBD' => 1, 'HBA' => 1,
					'PI' => 1, 'NI' => 1,
					'Ar' => 1,
					'Hal' => 1,
					'H' => 0,
					'RA' => 0, 'CA' => 0);
  return $This;
}

# Initialize object properties...
#
sub _InitializeFunctionalClassAtomTypesProperties {
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

# Disable change of AvailableFunctionalClasses...
#
sub SetAvailableFunctionalClasses {
  my($This) = @_;

  carp "Warning: ${ClassName}->SetFunctionalClassesOrder: Available functional classes can't be changed...";

  return $This;
}

# Disable change of functional classes order used for generation of atom types...
#
sub SetFunctionalClassesOrder {
  my($This) = @_;

  carp "Warning: ${ClassName}->SetFunctionalClassesOrder: functional classes order can't be changed...";

  return $This;
}

# Set functional classes to use for atom types...
#
sub SetFunctionalClassesToUse {
  my($This, @Values) = @_;
  my($FirstValue, $TypeOfFirstValue, $FunctionalClass, $SpecifiedFunctionalClass, $FunctionalClassValue, @SpecifiedFunctionalClasses, %FunctionalClassesToUse);

  if (!@Values) {
    carp "Warning: ${ClassName}->SetFunctionalClassesToUse: No values specified...";
    return;
  }

  $FirstValue = $Values[0];
  $TypeOfFirstValue = ref $FirstValue;
  @SpecifiedFunctionalClasses = ();

  if ($TypeOfFirstValue =~ /^ARRAY/) {
    push @SpecifiedFunctionalClasses, @{$FirstValue};
  }
  else {
    push @SpecifiedFunctionalClasses, @Values;
  }

  # Make sure specified FunctionalClasses are valid...
  for $SpecifiedFunctionalClass (@SpecifiedFunctionalClasses) {
    if (exists $AvailableFunctionalClasses{$SpecifiedFunctionalClass}) {
      $FunctionalClass = $SpecifiedFunctionalClass;
    }
    elsif ($AvailableFunctionalClassesByDescription{$SpecifiedFunctionalClass}) {
      $FunctionalClass = $AvailableFunctionalClassesByDescription{$SpecifiedFunctionalClass};
    }
    else {
      croak "Error: ${ClassName}->SetFunctionalClassesToUse: Specified functional class, $SpecifiedFunctionalClass, is not supported...\n ";
    }
    $FunctionalClassesToUse{$FunctionalClass} = 1;
  }

  # Set functional classes...
  for $FunctionalClass (keys %{$This->{FunctionalClassesToUse}}) {
    $This->{FunctionalClassesToUse}{$FunctionalClass} = 0;
    if (exists $FunctionalClassesToUse{$FunctionalClass}) {
      $This->{FunctionalClassesToUse}{$FunctionalClass} = 1;
    }
  }

  return $This;
}

# Is it an available FunctionalClass?
#
sub IsFunctionalClassAvailable {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $FunctionalClass, $Status);

  if ((@_ == 2) && (_IsFunctionalClassAtomTypes($FirstParameter))) {
    ($This, $FunctionalClass) = ($FirstParameter, $SecondParameter);
  }
  else {
    $FunctionalClass = $FirstParameter;
  }
  $Status = exists($AvailableFunctionalClasses{$FunctionalClass}) || exists($AvailableFunctionalClassesByDescription{$FunctionalClass}) ? 1 : 0;

  return $Status;
}

# Get a hash containing available functional classes and their description
# as key/value pairs.
#
sub GetAvailableFunctionalClasses {
  return %AvailableFunctionalClasses;
}

# Get an array containing order of functional classes used to generate atom types...
#
sub GetFunctionalClassesOrder {
  return @FunctionalClassesOrder;
}

# Assign functional class atom types to all atoms...
#
# Let:
#   HBD: HydrogenBondDonor
#   HBA: HydrogenBondAcceptor
#   PI :  PositivelyIonizable
#   NI : NegativelyIonizable
#   Ar : Aromatic
#   Hal : Halogen
#   H : Hydrophobic
#   RA : RingAtom
#   CA : ChainAtom
#
# Then:
#
#   Function class atom type specification for an atom corresponds to:
#
#     Ar.CA.H.HBA.HBD.Hal.NI.PI.RA
#
#   Default functional classes used are: HBD, HBA, PI, NI, Ar, Hal
#
#   FunctionalAtomTypes are assigned using the following definitions [ Ref 60-61, Ref 65-66 ]:
#
#     HydrogenBondDonor: NH, NH2, OH
#     HydrogenBondAcceptor: N[!H], O
#     PositivelyIonizable: +, NH2
#     NegativelyIonizable: -, C(=O)OH, S(=O)OH, P(=O)OH
#
# Notes:
#   . Final functional class atom type shows only those functional
#     classes to which an atom belongs; others are not shown.
#   . A null string is assigned as final atom type to those atom which
#     don't belong to any of the specified functional classes.
#
# Examples of functional class atom types:
#
# HBD.HBA - Hydrogen bond donor and acceptor
# HBD.RA - Hydrogen bond donor in a ring
#
#
sub AssignAtomTypes {
  my($This) = @_;
  my($Atom, $AtomType, $FunctionalClass, $FunctionalClassValue, @FunctionalClasses);

  ATOM: for $Atom ($This->GetMolecule()->GetAtoms()) {
    if ($This->{IgnoreHydrogens} && $Atom->IsHydrogen()) {
      next ATOM;
    }
    @FunctionalClasses = ();

    # Go over functional classes...
    ATOMICINVARIANT: for $FunctionalClass (@FunctionalClassesOrder) {
      if (!$This->{FunctionalClassesToUse}{$FunctionalClass}) {
	next ATOMICINVARIANT;
      }
      if ($Atom->IsFunctionalClassType($FunctionalClass)) {
	push @FunctionalClasses, $FunctionalClass;
      }
    }
    # Create and assign atom type to atom...
    $AtomType = (scalar @FunctionalClasses) ? TextUtil::JoinWords(\@FunctionalClasses, ".", 0) : "None";
    $This->SetAtomType($Atom, $AtomType);
  }
  return $This;
}

# Are all atoms types successfully assigned?
#
# Notes:
#   . Base class method is overridden to always return 1: An appropriate value, functional
#     class types delmited by dot or None, is always assigned to atoms.
#
sub IsAtomTypesAssignmentSuccessful {
  my($This) = @_;

  return 1;
}

# Return a string containg data for FunctionalClassAtomTypes object...
#
sub StringifyFunctionalClassAtomTypes {
  my($This) = @_;
  my($AtomTypesString);

  # Type of AtomTypes...
  $AtomTypesString = "AtomTypes: $This->{Type}; IgnoreHydrogens: " . ($This->{IgnoreHydrogens} ? "Yes" : "No");

  # AvailableFunctionalClasses and FunctionalClassesToUse...
  my($FunctionalClass, @FunctionalClasses, @FunctionalClassesToUse);

  @FunctionalClassesToUse = ();
  @FunctionalClasses = ();
  for $FunctionalClass (@FunctionalClassesOrder) {
    push @FunctionalClasses, "$FunctionalClass: $AvailableFunctionalClasses{$FunctionalClass}";
    if ($This->{FunctionalClassesToUse}{$FunctionalClass}) {
      push @FunctionalClassesToUse, $FunctionalClass;
    }
  }
  $AtomTypesString .= "; FunctionalClassesToUse: <" . TextUtil::JoinWords(\@FunctionalClassesToUse, ", ", 0) . ">";
  $AtomTypesString .= "; FunctionalClassesOrder: <" . TextUtil::JoinWords(\@FunctionalClassesOrder, ", ", 0) . ">";
  $AtomTypesString .= "; AvailableFunctionalClasses: <" . TextUtil::JoinWords(\@FunctionalClasses, ", ", 0) . ">";

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

# Is it a FunctionalClassAtomTypes object?
sub _IsFunctionalClassAtomTypes {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

FunctionalClassAtomTypes

=head1 SYNOPSIS

use AtomTypes::FunctionalClassAtomTypes;

use AtomTypes::FunctionalClassAtomTypes qw(:all);

=head1 DESCRIPTION

B<FunctionalClassAtomTypes> class provides the following methods:

new, AssignAtomTypes, GetAvailableFunctionalClasses, GetFunctionalClassesOrder,
IsFunctionalClassAvailable, SetFunctionalClassesToUse, StringifyFunctionalClassAtomTypes

B<FunctionalClassAtomTypes> is derived from B<AtomTypes> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<FunctionalClassAtomTypes>, B<AtomTypes> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

Possible values for functional clas atom types are: I<Ar, CA, H, HBA, HBD, Hal, NI, PI, RA>.
Default value: I<HBD, HBA, PI, NI, Ar, Hal>.

The functional calss atom types abbreviations correspond to:

    HBD: HydrogenBondDonor
    HBA: HydrogenBondAcceptor
    PI :  PositivelyIonizable
    NI : NegativelyIonizable
    Ar : Aromatic
    Hal : Halogen
    H : Hydrophobic
    RA : RingAtom
    CA : ChainAtom

FunctionalAtomTypes are assigned using the following definitions [ Ref 60-61, Ref 65-66 ]:

    HydrogenBondDonor: NH, NH2, OH
    HydrogenBondAcceptor: N[!H], O
    PositivelyIonizable: +, NH2
    NegativelyIonizable: -, C(=O)OH, S(=O)OH, P(=O)OH

Notes:

    o Final functional class atom type shows only those functional
      classes to which an atom belongs; others are not shown.
    o A null string is assigned as final atom type to those atom which
      don't belong to any of the specified functional classes.

 Examples of functional class atom types:

    HBD.HBA - Hydrogen bond donor and acceptor
    HBD.RA - Hydrogen bond donor in a ring

=head2 METHODS

=over 4

=item B<new>

    $NewFunctionalClassAtomTypes = new AtomTypes::FunctionalClassAtomTypes(
                                                   %NamesAndValues);

Using specified I<FunctionalClassAtomTypes> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<FunctionalClassAtomTypes>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'FunctionalClass'
    IgnoreHydrogens = 0
    FunctionalClassesToUse = HBD, HBA, PI, NI, Ar, Hal

Examples:

    $FunctionalClassAtomTypes = new AtomTypes::FunctionalClassAtomTypes(
                              'Molecule' => $Molecule,
                              'IgnoreHydrogens' => 0,
                              'FunctionalClassesToUse' =>
                                         ['HBD', 'HBA', 'PI', 'NI', 'Ar', 'Hal']);

=item B<AssignAtomTypes>

    $FunctionalClassAtomTypes->AssignAtomTypes();

Assigns functional class atom types to all the atoms in a molecule and returns
I<FunctionalClassAtomTypes>.

=item B<GetAvailableFunctionalClasses>

    %AvailableFunctionalClasses = $FunctionalClassAtomTypes->
                                 GetAvailableFunctionalClasses();

Returns available functional classes as a hash containing available functional classes
and their description as key/value pairs.

=item B<GetFunctionalClassesOrder>

    @FunctionalClassesOrder = $FunctionalClassAtomTypes->
                             GetFunctionalClassesOrder();

Returns an array obtaining order of functional classes used to generate atom types.

=item B<IsAtomTypesAssignmentSuccessful>

    $Status = $AtomTypes->IsAtomTypesAssignmentSuccessful();

Returns 1 or 0 based on whether atom types assignment was successfully performed.
This method overrides the same method available in the base class AtomTypes.pm used
to derived this class.

=item B<IsFunctionalClassAvailable>

    $Status = $FunctionalClassAtomTypes->
              IsFunctionalClassAvailable($FunctionalClass);
    $Status = AtomTypes::FunctionalClassAtomTypes::
              IsFunctionalClassAvailable($FunctionalClass);

Returns 1 or 0 based on whether I<FunctionalClass> is valid.

=item B<SetFunctionalClassesToUse>

    $FunctionalClassAtomTypes->SetFunctionalClassesToUse($ValuesRef);
    $FunctionalClassAtomTypes->SetFunctionalClassesToUse(@Values);

Set functional classes to use for generating and assigning atom types and returns
I<FunctionalClassAtomTypes>.

=item B<StringifyFunctionalClassAtomTypes>

    $String = $FunctionalClassAtomTypes->StringifyFunctionalClassAtomTypes();

Returns a string containing information about I<FunctionalClassAtomTypes> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AtomTypes.pm, AtomicInvariantsAtomTypes.pm, DREIDINGAtomTypes.pm,
EStateAtomTypes.pm, MMFF94AtomTypes.pm, SLogPAtomTypes.pm,
SYBYLAtomTypes.pm, TPSAAtomTypes.pm, UFFAtomTypes.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
