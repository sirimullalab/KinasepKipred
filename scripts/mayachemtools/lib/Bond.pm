package Bond;
#
# File: Bond.pm
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
use Storable ();
use Scalar::Util ();
use ObjectProperty;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(ObjectProperty Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, $ObjectID);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyBond';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializeBond();

  $This->_InitializeBondProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeBond {
  my($This) = @_;
  my($ObjectID) = _GetNewObjectID();

  # All other property names and values along with all Set/Get<PropertyName> methods
  # are implemented on-demand using ObjectProperty class.
  $This->{ID} = $ObjectID;

  # Bond from and to atoms indicate begining and ending bond atoms...
  $This->{BondFromAtom} = undef;
  $This->{BondToAtom} = undef;

  $This->{BondOrder} = '';
  $This->{BondType} = '';
  $This->{BondStereochemistry} = '';
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # ID to keep track of objects...
  $ObjectID = 0;
}

# Initialize bond properties...
sub _InitializeBondProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  if (!exists $NamesAndValues{'Atoms'}) {
    carp "Warning: ${ClassName}->new: Bond object instantiated without specifying atoms list...";
  }
  if (!exists $NamesAndValues{'BondOrder'}) {
    carp "Warning: ${ClassName}->new: Bond object instantiated without setting bond order...";
  }
  return $This;
}

# Setup an explicit SetID method to block setting of ID by AUTOLOAD function...
sub SetID {
  my($This, $Value) = @_;

  carp "Warning: ${ClassName}->SetID: Object ID can't be changed: it's used for internal tracking...";

  return $This;
}

# Setup an explicit SetMolecule method to block setting of ID by AUTOLOAD function...
sub SetMolecule {
  my($This, $Value) = @_;

  carp "Warning: ${ClassName}->SetMolecule: Molecule property can't be changed: it's used for internal tracking...";

  return $This;
}

# Assign bond to  molecule...
sub _SetMolecule {
  my($This, $Molecule) = @_;

  $This->{Molecule} = $Molecule;

  # Weaken the reference to disable increment of reference count; otherwise,
  # it it becomes a circular reference and destruction of Molecule object doesn't
  # get initiated which in turn disables destruction of bond object.
  #
  Scalar::Util::weaken($This->{Molecule});

  return $This;
}

# Set bond atoms...
sub SetAtoms {
  my($This, @Values) = @_;

  if (!@Values) {
    croak "Error: ${ClassName}->SetAtoms: No atoms specified...";
  }

  my($FirstValue, $TypeOfFirstValue, $Atom1, $Atom2, $AtomID1, $AtomID2);
  $FirstValue = $Values[0];
  $TypeOfFirstValue = ref $FirstValue;

  if ($TypeOfFirstValue =~ /^ARRAY/) {
    # Initialize using array refernce...
    if (@{$FirstValue} != 2) {
      croak "Warning: ${ClassName}->SetAtoms: Number of atoms specified in bond object is not equal to 2...";
    }
    ($Atom1, $Atom2) = @{$FirstValue};
  }
  else {
    # It's a list of values...
    if (@Values != 2) {
      croak "Warning: ${ClassName}->SetAtoms: Number of atoms specified in bond object is not equal to 2...";
    }
    ($Atom1, $Atom2) = @Values;
  }

  $AtomID1 = $Atom1->GetID(); $AtomID2 = $Atom2->GetID();
  if ($AtomID1 == $AtomID2) {
      croak "Warning: ${ClassName}->SetAtoms: Can't specify same atoms to create a bond...";
  }

  $This->{BondFromAtom} = $Atom1;
  $This->{BondToAtom} = $Atom2;

  return $This;
}

# Get bond atoms as array. In scalar context, return number of atoms involved in
# bond...
#
sub GetAtoms {
  my($This) = @_;

  return wantarray ? ($This->{BondFromAtom}, $This->{BondToAtom}) : 2;
}

# Get atom bonded to specified atom...
sub GetBondedAtom {
  my($This, $Atom) = @_;
  my($Atom1, $Atom2, $AtomID1, $AtomID2, $AtomID);

  ($Atom1, $Atom2) = $This->GetAtoms();
  $AtomID1 = $Atom1->GetID(); $AtomID2 = $Atom2->GetID(); $AtomID = $Atom->GetID();

  return ($AtomID1 == $AtomID) ? $Atom2 : (($AtomID2 == $AtomID) ? $Atom1 : undef) ;
}

# Get common atom between two bonds...
sub GetCommonAtom {
  my($This, $Other) = @_;
  my($Atom1, $Atom2, $AtomID1, $AtomID2, $OtherAtom1, $OtherAtom2, $OtherAtomID1, $OtherAtomID2);

  ($Atom1, $Atom2) = $This->GetAtoms();
  $AtomID1 = $Atom1->GetID(); $AtomID2 = $Atom2->GetID();

  ($OtherAtom1, $OtherAtom2) = $Other->GetAtoms();
  $OtherAtomID1 = $OtherAtom1->GetID(); $OtherAtomID2 = $OtherAtom2->GetID();

  return ($AtomID1 == $OtherAtomID1 || $AtomID1 == $OtherAtomID2) ? $Atom1 : (($AtomID2 == $OtherAtomID1 || $AtomID2 == $OtherAtomID2) ? $Atom2 : undef) ;
}

# Get bond from atom indicating begining atom of a bond...
sub GetBondFromAtom {
  my($This) = @_;

  return $This->{BondFromAtom};
}

# Get begin bond atom...
sub GetBondBeginAtom {
  my($This) = @_;

  return $This->GetBondFromAtom();
}

# Get bond to atom indicating ending atom of a bond...
sub GetBondToAtom {
  my($This) = @_;

  return $This->{BondToAtom};
}

# Get begin bond atom...
sub GetBondEndAtom {
  my($This) = @_;

  return $This->GetBondToAtom();
}

# Switch bond from and to atoms...
sub SwitchBondFromAndToAtoms {
  my($This) = @_;
  my($FromAtom, $ToAtom);

  ($FromAtom, $ToAtom) = $This->GetAtoms();

  $This->{BondFromAtom} = $ToAtom;
  $This->{BondToAtom} = $FromAtom;

  return $This;
}

# Set bond order...
#
# Possible BondOrder are: 0 = None, 1 = Single, 1.5 = Aromatic, 2 = Double, 3 = Triple, 4 = Quadruple,
# 5 = Quintuple, 6 = Sextuple, 7 = Septuple
#
# Possible BondType for different BondOrders are:
#
# 0 : None, Ionic, Unspecified
# 1 : Single, Dative, Coordinate, SingleOrDouble, SingleOrAromatic, Tautomeric
# 1.5 : Aromatic, Resonance, SingleOrAromatic, DoubleOrAromatic
# 2 : Double, Tautomeric, SingleOrDouble, DoubleOrAromatic
# 3 : Triple
# 4 : Quadruple
# 5 : Quintuple
# 6 : Sextuple
# 7 : Septuple
#
# Note: BondType Any is valid for all BondOrders.
#
# Possible BondStereochemistry values for different BondOrders are:
#
# 0 : None, Unspecified
# 1 : Wedge, Up, Hash, Down, Wavy, WedgeOrHash, UpOrDown, Upward, Downward, None, Unspecified
# 2 : Cis, Trans, Z, E, DoubleCross, CisOrTrans, None, Unspecified
#
# Notes:
#   . BondType property is automatically assigned using default BondType values for
#     specified BondOrder.
#   . BondType values can also be explicit set.
#   . To make bonds aromatic in a ring, explicitly set "Aromatic" property for bond/atoms and make sure
#     appropriate BondOrder values are assigned.
#   . Dative or coordinate bond types are treated as single bond types with explicit formal charge
#     of + and - on first and second bond atoms.
#
sub SetBondOrder {
  my($This, $BondOrder) = @_;

  if ($BondOrder !~ /^(0|1|1.5|2|3|4|5|6|7)$/) {
    croak "Error: ${ClassName}->SetBondOrder: BondOrder value $BondOrder is not valid. Supported values: 0, 1, 1.5, 2, 3, 4, 5, 6, 7";
  }

  # Set bond order and type...
  my($BondType);

  BONDORDER: {
      if ($BondOrder == 1) {$BondType = 'Single'; last BONDORDER; }
      if ($BondOrder == 1.5) {$BondType = 'Aromatic'; last BONDORDER; }
      if ($BondOrder == 2) {$BondType = 'Double'; last BONDORDER; }
      if ($BondOrder == 3) {$BondType = 'Triple'; last BONDORDER; }
      if ($BondOrder == 4) {$BondType = 'Quadruple'; last BONDORDER; }
      if ($BondOrder == 5) {$BondType = 'Quintuple'; last BONDORDER; }
      if ($BondOrder == 6) {$BondType = 'Sextuple'; last BONDORDER; }
      if ($BondOrder == 7) {$BondType = 'Septuple'; last BONDORDER; }
      $BondType = '';
      $BondOrder = 0;
  }
  $This->{BondType} = $BondType;
  $This->{BondOrder} = $BondOrder;

  return $This;

}

# Set bond type for a specific bond...
#
sub SetBondType {
  my($This, $BondType) = @_;

  if ($BondType !~ /^(None|Ionic|Unspecified|Single|Dative|Coordinate|SingleOrDouble|SingleOrAromatic|Aromatic|Tautomeric|Resonance|DoubleOrAromatic|Double|Triple|Quadruple|Any)$/i) {
    croak "Error: ${ClassName}->SetBondType: BondType value $BondType is not valid. Supported values: None, Ionic, Unspecified, Single, Dative, Coordinate, SingleOrDouble, SingleOrAromatic, Aromatic, Tautomeric, Resonance, DoubleOrAromatic, Double, Triple, Quadruple, Any...";
  }

  # Make sure its a valid BondType value for BondOrder...
  my($BondOrder, $ValidBondType);

  $ValidBondType = 0;
  $BondOrder = $This->{BondOrder};

  BONDORDER: {
      if ($BondOrder == 0 && $BondType =~ /^(None|Ionic|Unspecified|Any)$/i) {$ValidBondType = 1; last BONDORDER; }
      if ($BondOrder == 1 && $BondType =~ /^(Single|Dative|Coordinate|SingleOrDouble|SingleOrAromatic|Tautomeric|Any)$/i) {$ValidBondType = 1; last BONDORDER; }
      if ($BondOrder == 1.5 && $BondType =~ /^(Aromatic|Resonance|SingleOrAromatic|DoubleOrAromatic|Any)$/i) {$ValidBondType = 1; last BONDORDER; }
      if ($BondOrder == 2 && $BondType =~ /^(Double|SingleOrDouble|DoubleOrAromatic|Tautomeric|Any)$/i ) {$ValidBondType = 1; last BONDORDER; }
      if ($BondOrder == 3 && $BondType =~ /^(Triple|Any)$/i) {$ValidBondType = 1; last BONDORDER; }
      if ($BondOrder == 4 && $BondType =~ /^(Quadruple|Any)$/i) {$ValidBondType = 1; last BONDORDER; }
      if ($BondOrder == 5 && $BondType =~ /^(Quintuple|Any)$/i) {$ValidBondType = 1; last BONDORDER; }
      if ($BondOrder == 6 && $BondType =~ /^(Sextuple|Any)$/i) {$ValidBondType = 1; last BONDORDER; }
      if ($BondOrder == 7 && $BondType =~ /^(Septuple|Any)$/i) {$ValidBondType = 1; last BONDORDER; }
      $ValidBondType = 0;
  }

  if (!$ValidBondType) {
    carp "Warning: ${ClassName}->SetBondType: Assigning invalid BondType value, $BondType, for BondOrder $BondOrder...";
  }

  $This->{BondType} = $BondType;

  # Set explicit formal charge for atoms involved in Dative or coordinate bonds...
  if ($BondType =~ /^(Dative|Coordinate)$/i) {
    my($Atom1, $Atom2) = $This->GetAtoms();
    $Atom1->SetFormalCharge('1');
    $Atom2->SetFormalCharge('-1');
  }

  return $This;
}

# Set bond stereochemistry...
#
# Single bond stereochemistry:
#
# None, Unspecified: Not a stereo bond or unspecified
#
# Wedge, Up : Wedge end pointing up
# Hash, Down: Wedge end pointing down
# Wavy, WedgeOrHash, UpOrDown: Wedge end up or down
#
# Note: Wedge starts at begin atom of bond making wedge pointed end always at this
#       atom.
#
# Upward: Single bond around cis/trans double bonds pointing upward
# Downward: Single bond around cis/trans double bonds pointing upward
#
# Note: Upward/downward bonds start at atoms involved in cis/trans double bond.
#
# Double bond stereochemistry:
#
# None, Unspecified: Not a stereo bond or unspecified
#
# Z, cis: Similar groups on same side of double bond
# E, trans: Similar groups on different side of double bond
#
# CisOrTrans, DoubleCross: cis or trans
#
sub SetBondStereochemistry {
  my($This, $BondStereochemistry) = @_;

  if ($BondStereochemistry !~ /^(None|Unspecified|Wedge|Hash|Up|Down|Wavy|WedgeOrHash|UpOrDown|Upward|Downward|Cis|Trans|Z|E|DoubleCross|CisOrTrans)$/i) {
    croak "Error: ${ClassName}->SetBondStereochemistry: BondStereochemistry value $BondStereochemistry is not valid. Supported values: None, Unspecified, Wedge, Hash, Up, Down, Wavy, WedgeOrHash, UpOrDown, Upward, Downward, Cis, Trans, Z, E, DoubleCross, CisOrTrans...";
  }

  # Make sure its a valid BondStereochemistry value for BondOrder...
  my($BondOrder, $ValidBondType);

  $ValidBondType = 0;
  $BondOrder = $This->{BondOrder};

  BONDORDER: {
      if ($BondOrder == 0 && $BondStereochemistry =~ /^(None|Unspecified)$/i) {$ValidBondType = 1; last BONDORDER; }
      if ($BondOrder == 1 && $BondStereochemistry =~ /^(Wedge|Hash|Up|Down|Wavy|WedgeOrHash|UpOrDown|Upward|Downward|None|Unspecified)$/i) {$ValidBondType = 1; last BONDORDER; }
      if ($BondOrder == 2 && $BondStereochemistry =~ /^(Cis|Trans|Z|E|DoubleCross|CisOrTrans|None|Unspecified)$/i ) {$ValidBondType = 1; last BONDORDER; }
      $ValidBondType = 0;
  }

  if (!$ValidBondType) {
    carp "Warning: ${ClassName}->SetBondStereochemistry: Assigning invalid BondStereochemistry value, $BondStereochemistry, for BondOrder $BondOrder...";
  }

  $This->{BondStereochemistry} = $BondStereochemistry;

  return $This;
}

# Is it a single bond?
sub IsSingle {
  my($This) = @_;

  return ($This->{BondOrder} == 1) ? 1 : 0;
}

# Is it a double bond?
sub IsDouble {
  my($This) = @_;

  return ($This->{BondOrder} == 2) ? 1 : 0;
}

# Is it a triple bond?
sub IsTriple {
  my($This) = @_;

  return ($This->{BondOrder} == 3) ? 1 : 0;
}

# Is it a quadruple bond?
sub IsQuadruple {
  my($This) = @_;

  return ($This->{BondOrder} == 4) ? 1 : 0;
}

# Is aromatic property set for the bond?
sub IsAromatic {
  my($This) = @_;
  my($Aromatic);

  $Aromatic = $This->GetAromatic();

  return ((defined($Aromatic) && $Aromatic) || $This->{BondOrder} == 1.5) ? 1 : 0;
}

# Is it a quintuple bond?
sub IsQuintuple {
  my($This) = @_;

  return ($This->{BondOrder} == 5) ? 1 : 0;
}

# Is it a sextuple bond?
sub IsSextuple {
  my($This) = @_;

  return ($This->{BondOrder} == 6) ? 1 : 0;
}

# Is bond type specified?
sub IsBondTypeSpecified {
  my($This) = @_;

  return ($This->{BondType} && $This->{BondType} !~ /^(None|Unspecified)$/i) ? 1 : 0;
}

# Is it a dative or coordinate bond?
sub IsDative {
  my($This) = @_;

  return ($This->{BondType} =~ /^(Dative|Coordinate)$/i) ? 1 : 0;
}

# Is it a dative or coordinate bond?
sub IsCoordinate {
  my($This) = @_;

  return $This->IsDative();
}

# Is it a ionic bond?
sub IsIonic {
  my($This) = @_;

  return ($This->{BondType} =~ /^Ionic$/i) ? 1 : 0;
}

# Is it a tautomeric bond?
sub IsTautomeric {
  my($This) = @_;

  return ($This->{BondType} =~ /^Tautomeric$/i) ? 1 : 0;
}

# Is bond stereochemistry specified?
sub IsBondStereochemistrySpecified {
  my($This) = @_;

  return ($This->{BondStereochemistry} && $This->{BondStereochemistry} !~ /^(None|Unspecified)$/i) ? 1 : 0;
}

# Is it a Wedge or Up single bond?
sub IsWedge {
  my($This) = @_;

  return $This->IsUp();
}

# Is it a Wedge or Up single bond?
sub IsUp {
  my($This) = @_;

  return ($This->{BondStereochemistry} =~ /^(Wedge|Up)$/i) ? 1 : 0;
}

# Is it a Hash or Down single bond?
sub IsHash {
  my($This) = @_;

  return $This->IsDown();
}

# Is it a Hash or Down single bond?
sub IsDown {
  my($This) = @_;

  return ($This->{BondStereochemistry} =~ /^(Hash|Down)$/i) ? 1 : 0;
}

# Is it a Wavy, WedgeOrHash or UpOrDown single bond?
sub IsWedgeOrHash {
  my($This) = @_;

  return $This->IsUpOrDown();
}

# Is it a Wavy, WedgeOrHash or UpOrDown single bond?
sub IsUpOrDown {
  my($This) = @_;

  return ($This->{BondStereochemistry} =~ /^(Wavy|WedgeOrHash|UpOrDown)$/i) ? 1 : 0;
}

# Is it a upward single bond?
sub IsUpward {
  my($This) = @_;

  return ($This->{BondStereochemistry} =~ /^Upward$/i) ? 1 : 0;
}

# Is it a Downward single bond?
sub IsDownward {
  my($This) = @_;

  return ($This->{BondStereochemistry} =~ /^Downward$/i) ? 1 : 0;
}
# Is it a cis or Z double bond?
sub IsCis {
  my($This) = @_;

  return ($This->{BondStereochemistry} =~ /^(Cis|Z)$/i) ? 1 : 0;
}

# Is it a trans or E double bond?
sub IsTrans {
  my($This) = @_;

  return ($This->{BondStereochemistry} =~ /^(Trans|E)$/i) ? 1 : 0;
}

# Is it a cis or trans double bond?
sub IsCisOrTrans {
  my($This) = @_;

  return ($This->{BondStereochemistry} =~ /^(CisOrTrans|DoubleCross)$/i) ? 1 : 0;
}

# Delete bond...
sub DeleteBond {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    # Nothing to do...
    return $This;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');
  $Molecule->DeleteBond($This);

  return $This;
}

# Copy bond and all its associated data...
sub Copy {
  my($This) = @_;
  my($Bond);

  $Bond = Storable::dclone($This);

  return $Bond;
}

# Is bond in a ring?
#
sub IsInRing {
  my($This) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_IsBondInRing($This);
}

# Is bond not in a ring?
#
sub IsNotInRing {
  my($This) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_IsBondNotInRing($This);
}

# Is bond only in one ring?
#
sub IsOnlyInOneRing {
  my($This) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_IsBondInOnlyOneRing($This);
}

# Is bond in a ring of specific size?
#
sub IsInRingOfSize {
  my($This, $RingSize) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_IsBondInRingOfSize($This, $RingSize);
}

# Get size of smallest ring containing the bond...
#
sub GetSizeOfSmallestRing {
  my($This) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetSizeOfSmallestBondRing($This);
}

# Get size of largest ring containing the bond...
#
sub GetSizeOfLargestRing {
  my($This) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetSizeOfLargestBondRing($This);
}

# Get number of  rings containing the bond...
#
sub GetNumOfRings {
  my($This) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetNumOfBondRings($This);
}

# Get number of  rings with odd size containing the bond...
#
sub GetNumOfRingsWithOddSize {
  my($This) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetNumOfBondRingsWithOddSize($This);
}

# Get number of  rings with even size containing the bond...
#
sub GetNumOfRingsWithEvenSize {
  my($This) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetNumOfBondRingsWithEvenSize($This);
}

# Get number of  rings with specified size containing the bond...
#
sub GetNumOfRingsWithSize {
  my($This, $RingSize) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetNumOfBondRingsWithSize($This, $RingSize);
}

# Get number of  rings with size less than specified containing the bond...
#
sub GetNumOfRingsWithSizeLessThan {
  my($This, $RingSize) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetNumOfBondRingsWithSizeLessThan($This, $RingSize);
}

# Get number of  rings with size greater than specified size containing the bond...
#
sub GetNumOfRingsWithSizeGreaterThan {
  my($This, $RingSize) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetNumOfBondRingsWithSizeGreaterThan($This, $RingSize);
}

# Get all rings as an array of references to arrays containing ring atoms...
#
sub GetRings {
  my($This) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetBondRings($This);
}

# Get smallest ring as an array containing ring atoms...
#
sub GetSmallestRing {
  my($This) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetSmallestBondRing($This);
}

# Get largest ring as an array containing ring atoms...
#
sub GetLargestRing {
  my($This) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetLargestBondRing($This);
}

# Get odd size rings an array of references to arrays containing ring atoms...
#
sub GetRingsWithOddSize {
  my($This) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetBondRingsWithOddSize($This);
}

# Get even size rings an array of references to arrays containing ring atoms...
#
sub GetRingsWithEvenSize {
  my($This) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetBondRingsWithEvenSize($This);
}

# Get rings with specified size  an array of references to arrays containing ring atoms...
#
sub GetRingsWithSize {
  my($This, $RingSize) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetBondRingsWithSize($This, $RingSize);
}

# Get rings with size less than specfied size as an array of references to arrays containing ring atoms...
#
sub GetRingsWithSizeLessThan {
  my($This, $RingSize) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetBondRingsWithSizeLessThan($This, $RingSize);
}

# Get rings with size greater than specfied size as an array of references to arrays containing ring atoms...
#
sub GetRingsWithSizeGreaterThan {
  my($This, $RingSize) = @_;

  # Is this bond in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetBondRingsWithSizeGreaterThan($This, $RingSize);
}

# Get next object ID...
sub _GetNewObjectID {
  $ObjectID++;
  return $ObjectID;
}

# Return a string containing bond and other properties...
sub StringifyBond {
  my($This) = @_;
  my($BondString, $ID, $BondOrder, $BondType, $Stereochemistry, $AtomsString, $RingBond, $AromaticBond, $NumOfRings, $Atom1, $Atom2);

  $ID = $This->GetID();
  $BondOrder = $This->GetBondOrder();
  if (!defined $BondOrder) {
    $BondOrder = "undefined";
  }
  $BondType = $This->GetBondType();
  if (!defined $BondOrder) {
    $BondType = "undefined";
  }
  if (defined($BondOrder) && $BondOrder == 2) {
    $Stereochemistry = $This->GetStereochemistry();
    if (!defined $Stereochemistry) {
      $Stereochemistry = "undefined";
    }
  }
  $RingBond = $This->IsInRing();
  if (defined $RingBond) {
    $RingBond = $RingBond  ? 'Yes' : 'No';
    $NumOfRings = $This->GetNumOfRings();
  }
  else {
    $RingBond = 'undefined';
    $NumOfRings = 'undefined';
  }

  $AromaticBond = $This->GetAromatic();
  if (defined $AromaticBond) {
    $AromaticBond = $AromaticBond  ? 'Yes' : 'No';
  }
  else {
    $AromaticBond = 'undefined';
  }

  ($Atom1, $Atom2) = $This->GetAtoms();
  $AtomsString = "Atoms: undefined";
  if (defined($Atom1) && defined($Atom2)) {
    my($Atom, $AtomID, $AtomCount, $AtomName, $AtomSymbol, $AtomicNumber, @BondAtoms);
    @BondAtoms = ();
    push @BondAtoms, ($Atom1, $Atom2);
    $AtomCount = 0;
    $AtomsString = "";
    for $Atom (@BondAtoms) {
      $AtomCount++;
      $AtomID = $Atom->GetID();
      $AtomName = $Atom->GetName();
      $AtomSymbol = $Atom->GetAtomSymbol();
      $AtomicNumber = $Atom->GetAtomicNumber();
      if ($AtomCount == 1) {
	$AtomsString .= "FirstAtom:";
      }
      else {
	$AtomsString .= "; SecondAtom:";
      }
      $AtomsString .= " [ ID: $AtomID; Name: \"$AtomName\"; AtomSymbol: \"$AtomSymbol\"; AtomicNumber: $AtomicNumber ]";
    }
  }
  $BondString = "Bond: ID: $ID; BondOrder: $BondOrder; BondType: $BondType; RingBond: $RingBond; NumOfBondRings: $NumOfRings; AromaticBond: $AromaticBond;";
  if (defined($BondOrder) && $BondOrder == 2) {
  $BondString .= " Stereochemistry: $Stereochemistry;";
  }
  $BondString .= " $AtomsString";

  return $BondString;
}

1;

__END__

=head1 NAME

Bond

=head1 SYNOPSIS

use Bond;

use Bond qw(:all);

=head1 DESCRIPTION

B<Bond> class provides the following methods:

new, Copy, DeleteBond, GetAtoms, GetBondBeginAtom, GetBondEndAtom,
GetBondFromAtom, GetBondToAtom, GetBondedAtom, GetCommonAtom, GetLargestRing,
GetNumOfRings, GetNumOfRingsWithEvenSize, GetNumOfRingsWithOddSize,
GetNumOfRingsWithSize, GetNumOfRingsWithSizeGreaterThan,
GetNumOfRingsWithSizeLessThan, GetRings, GetRingsWithEvenSize,
GetRingsWithOddSize, GetRingsWithSize, GetRingsWithSizeGreaterThan,
GetRingsWithSizeLessThan, GetSizeOfLargestRing, GetSizeOfSmallestRing,
GetSmallestRing, IsAromatic, IsBondStereochemistrySpecified, IsBondTypeSpecified,
IsCis, IsCisOrTrans, IsCoordinate, IsDative, IsDouble, IsDown, IsDownward, IsHash,
IsInRing, IsInRingOfSize, IsIonic, IsNotInRing, IsOnlyInOneRing, IsQuadruple,
IsQuintuple, IsSextuple, IsSingle, IsTautomeric, IsTrans, IsTriple, IsUp,
IsUpOrDown, IsUpward, IsWedge, IsWedgeOrHash, SetAtoms, SetBondOrder,
SetBondStereochemistry, SetBondType, StringifyBond, SwitchBondFromAndToAtoms

B<Bond> class is derived from B<ObjectProperty> base class which provides methods not explicitly
defined in B<Atom> or B<ObjectProperty> class using Perl's AUTOLOAD functionality. These methods
are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

=head2 METHODS

=over 4

=item B<new>

    $NewBond = new Bond([%PropertyNameAndValues]);

Using specified I<Bond> property names and values hash, B<new> method creates a new object
and returns a reference to newly created B<Bond> object. By default, following properties are
initialized:

    ID = SequentialObjectID
    @Atoms = ();
    BondType = ""
    BondOrder = ""

Except for I<ID> property, all other default properties and other additional properties can
be set during invocation of this method.

Examples:

    $Bond = new Bond();
    $DoubleBond = new Bond('Atoms' => [$Atom2, $Atom1],
                           'BondOrder' => 2);

=item B<Copy>

    $BondCopy = $Bond->Copy();

Copy I<Bond> and its associated data using B<Storable::dclone> and return a new
B<Bond> object.

=item B<DeleteBond>

    $Bond->DeleteBond();

Delete I<Bond> between atoms in from a molecule.

=item B<GetAtoms>

    @BondedAtoms = $Bond->GetAtoms();

Returns an array containing I<Atoms> invoved in I<Bond>.

=item B<GetBondedAtom>

    $BondedAtom = $Bond->GetBondedAtom($Atom);

Returns B<BondedAtom> bonded to I<Atom> in  I<Bond>.

=item B<GetBondBeginAtom>

    $BeginAtom = $Bond->GetBondBeginAtom();

Returns B<BeginAtom> corresponding to bond starting atom in I<Bond>.

=item B<GetBondEndAtom>

    $EndAtom = $Bond->GetBondEndAtom();

Returns B<EndAtom> corresponding to bond ending atom in I<Bond>.

=item B<GetBondFromAtom>

    $FromAtom = $Bond->GetBondFromAtom();

Returns B<FromAtom> corresponding to bond starting atom in I<Bond>.

=item B<GetBondToAtom>

    $ToAotm = $Bond->GetBondToAtom();

Returns B<ToAtom> corresponding to bond ending atom in I<Bond>.

=item B<GetCommonAtom>

    $CommonAtom = $Bond->GetCommonAtom($OtherBond);

Returns B<Atom> common to both I<Bond> and I<$OtherBond>.

=item B<GetLargestRing>

    @RingAtoms = $Bond->GetLargestRing();

Returns an array of ring I<Atoms> corresponding to the largest ring containing I<Bond>.
in a molecule

=item B<GetNumOfRings>

    $NumOfRings = $Bond->GetNumOfRings();

Returns number of rings containing I<Bond> in a molecule.

=item B<GetNumOfRingsWithEvenSize>

    $NumOfRings = $Bond->GetNumOfRingsWithEvenSize();

Returns number of rings with even size containing I<Bond> in a molecule.

=item B<GetNumOfRingsWithOddSize>

    $NumOfRings = $Bond->GetNumOfRingsWithOddSize();

Returns number of rings with odd size containing I<Bond> in a molecule.

=item B<GetNumOfRingsWithSize>

    $NumOfRings = $Bond->GetNumOfRingsWithSize($RingSize);

Returns number of rings with specific I<RingSize> containing I<Bond> in a molecule.

=item B<GetNumOfRingsWithSizeGreaterThan>

    $NumOfRings = $Bond->GetNumOfRingsWithSizeGreaterThan($RingSize);

Returns number of rings with size greater than specific I<RingSize> containing
I<Bond> in a molecule.

=item B<GetNumOfRingsWithSizeLessThan>

    $NumOfRings = $Bond->GetNumOfRingsWithSizeLessThan($RingSize);

Returns number of rings with size less than specific I<RingSize> containing I<Bond>
in a molecule.

=item B<GetRings>

    @Rings = $Bond->GetRings();

Returns an array of references to arrays containing ring atoms corressponding
to all rings containing I<Bond> in a molecule.

=item B<GetRingsWithEvenSize>

    @Rings = $Bond->GetRingsWithEvenSize();

Returns an array of references to arrays containing ring atoms corressponding to all rings with even size
containing I<Bond> in a molecule.

=item B<GetRingsWithOddSize>

    @Rings = $Bond->GetRingsWithOddSize();

Returns an array of references to arrays containing ring atoms corressponding to all rings with odd size
containing I<Bond> in a molecule.

=item B<GetRingsWithSize>

    @Rings = $Bond->GetRingsWithSize($RingSize);

Returns an array of references to arrays containing ring atoms corressponding to all rings with specific
I<RingSize >containing I<Bond> in a molecule.

=item B<GetRingsWithSizeGreaterThan>

    @Rings = $Bond->GetRingsWithSizeGreaterThan($RingSize);

Returns an array of references to arrays containing ring atoms corressponding to all rings with size
greater than specific I<RingSize >containing I<Bond> in a molecule.

=item B<GetRingsWithSizeLessThan>

    @Rings = $Bond->GetRingsWithSizeLessThan($RingSize);

Returns an array of references to arrays containing ring atoms corressponding to all rings with size
less than specific I<RingSize >containing I<Bond> in a molecule.

=item B<GetSizeOfLargestRing>

    $Size = $Bond->GetSizeOfLargestRing();

Returns size of the largest ring containing I<Bond> in a molecule.

=item B<GetSizeOfSmallestRing>

    $Size = $Bond->GetSizeOfSmallestRing();

Returns size of the smallest ring containing I<Bond> in a molecule.

=item B<GetSmallestRing>

    @RingAtoms = $Bond->GetSmallestRing();

Returns an array of ring I<Atoms> corresponding to the largest ring containing I<Bond>
in a molecule.

=item B<IsAromatic>

    $Status = $Bond->IsAromatic();

Returns 1 or 0 based on whether it's an aromatic I<Bond>.

=item B<IsBondStereochemistrySpecified>

    $Status = $Bond->IsBondStereochemistrySpecified();

Returns 1 or 0 based on whether I<Bond>'s sterochemistry is specified.

=item B<IsBondTypeSpecified>

    $Status = $Bond->IsBondTypeSpecified();

Returns 1 or 0 based on whether I<Bond>'s type is specified.

=item B<IsCis>

    $Status = $Bond->IsCis();

Returns 1 or 0 based on whether it's a cis I<Bond>.

=item B<IsCisOrTrans>

    $Status = $Bond->IsCisOrTrans();

Returns 1 or 0 based on whether it's a cis or trans I<Bond>.

=item B<IsCoordinate>

    $Status = $Bond->IsCoordinate();

Returns 1 or 0 based on whether it's a coordinate or dative  I<Bond>.

=item B<IsDative>

    $Status = $Bond->IsDative();

Returns 1 or 0 based on whether it's a coordinate or dative  I<Bond>.

=item B<IsDouble>

    $Status =$Bond->IsDouble();

Returns 1 or 0 based on whether it's a double I<Bond>.

=item B<IsDown>

    $Status = $Bond->IsDown();

Returns 1 or 0 based on whether it's a hash or down single I<Bond>.

=item B<IsDownward>

    $Return = $Bond->IsDownward();

Returns 1 or 0 based on whether it's a downward I<Bond>.

=item B<IsHash>

    $Status = $Bond->IsHash();

Returns 1 or 0 based on whether it's a hash or down single I<Bond>.

=item B<IsInRing>

    $Status = $Bond->IsInRing();

Returns 1 or 0 based on whether I<Bond> is present in a ring.

=item B<IsInRingOfSize>

    $Status = $Bond->IsInRingOfSize($Size);

Returns 1 or 0 based on whether I<Bond> is present in a ring of specific I<Size>.

=item B<IsIonic>

    $Status = $Bond->IsIonic();

Returns 1 or 0 based on whether it's an ionic I<Bond>.

=item B<IsNotInRing>

    $Status = $Bond->IsNotInRing();

Returns 1 or 0 based on whether I<Bond> is not present in a ring.

=item B<IsOnlyInOneRing>

    $Status = $Bond->IsOnlyInOneRing();

Returns 1 or 0 based on whether I<Bond> is only present in one ring.

=item B<IsQuadruple>

    $Status = $Bond->IsQuadruple();

Returns 1 or 0 based on whether it's a quadruple I<Bond>.

=item B<IsQuintuple>

    $Status = $Bond->IsQuintuple();

Returns 1 or 0 based on whether it's a quintuple I<Bond>.

=item B<IsSextuple>

    $Status = $Bond->IsSextuple();

Returns 1 or 0 based on whether it's a sextuple I<Bond>.

=item B<IsSingle>

    $Status =$Bond->IsSingle();

Returns 1 or 0 based on whether it's a single I<Bond>.

=item B<IsTriple>

    $Status =$Bond->IsTriple();

Returns 1 or 0 based on whether it's a triple I<Bond>.

=item B<IsTautomeric>

    $Status = $Bond->IsTautomeric();

Returns 1 or 0 based on whether it's a I<Bond>.

=item B<IsTrans>

    $Status = $Bond->IsTrans();

Returns 1 or 0 based on whether it's a trans I<Bond>.

=item B<IsUp>

    $Status = $Bond->IsUp();

Returns 1 or 0 based on whether it's a up I<Bond>.

=item B<IsUpOrDown>

    $Status = $Bond->IsUpOrDown();

Returns 1 or 0 based on whether it's an up or down I<Bond>.

=item B<IsUpward>

    $Status = $Bond->IsUpward();

Returns 1 or 0 based on whether it's an upward I<Bond>.

=item B<IsWedge>

    $Status = $Bond->IsWedge();

Returns 1 or 0 based on whether it's a wedge I<Bond>.

=item B<IsWedgeOrHash>

    $Status = $Bond->IsWedgeOrHash();

Returns 1 or 0 based on whether it's a wedge or hash I<Bond>.

=item B<SetAtoms>

    $Bond->SetAtoms($AtomsRef);
    $Bond->SetAtoms(@Atoms);

Set atoms of I<Bond> to atoms in I<Atoms> array or in a reference to an array of atoms
and return I<Bond>.

=item B<SetBondOrder>

    $Bond->SetBondOrder($BondOrder);

Sets bond order of I<Bond> to specified I<BondOrder> and returns I<Bond>. Possible bond order
values: 1 = Single, 1.5 = Aromatic, 2 = Double, 3 = Triple, 4 = Quadruple, 5 = Quintuple,
6 = Sextuple, 7 = Septuple

Notes:

    . BondType property is automatically assigned using default BondType
      values for specified BondOrder.
    . BondType values can also be explicit set.
    . To make bonds aromatic in a ring, explicitly set "Aromatic"
      property for bond/atoms and make sure appropriate BondOrder
      values are assigned.
    . Dative or coordinate bond types are treated as single bond types with
      explicit formal charge of + and - on first and second bond atoms.

=item B<SetBondType>

    $Bond->SetBondType($BondType);

Sets bond type for I<Bond> to specified I<BondType> and returns I<Bond>. Possible bond type
values for different bond orders are:

    0: None, Ionic, Unspecified
    1 : Single, Dative, Coordinate, SingleOrDouble, SingleOrAromatic, Tautomeric
    2 : Double, SingleOrDouble, DoubleOrAromatic, Tautomeric
    3 : Triple
    4 : Quadruple
    5 : Quintuple
    6 : Sextuple
    7 : Septuple
    1.5 : Aromatic, Resonance, SingleOrAromatic, DoubleOrAromatic

Notes:

    o BondType Any is valid for all BondOrders.
    o BondOrder property is automatically assigned using default BondOrder
      values for specified BondType.

Possible bond stereochemistry values for different bond orders are:

    0 : None, Unspecified
    1 : Wedge, Up, Hash, Down, Wavy, WedgeOrHash, UpOrDown, Upward, Downward,
        None, Unspecified
    2 : Cis, Trans, Z, E, DoubleCross, CisOrTrans, None, Unspecified

=item B<SetBondStereochemistry>

    $Bond = $Bond->SetBondStereochemistry($BondStereochemistry);

Sets bond stereochemistry of I<Bond> to specified I<BondStereochemistry> and
returns I<Bond>. Possible I<BondStereoChemistry> values for different bond orders
are:

BondOrder: 1

    None, Unspecified: Not a stereo bond or unspecified

    Wedge, Up : Wedge end pointing up
    Hash, Down: Wedge end pointing down
    Wavy, WedgeOrHash, UpOrDown: Wedge end up or down

    Upward: Single bond around cis/trans double bonds pointing upward
    Downward: Single bond around cis/trans double bonds pointing upward

Notes:

    o Wedge starts at begin atom of a bond making wedge pointed end always
      at this atom.
    o Upward/downward bonds start at atoms involved in cis/trans double bonds.

BondOrder: 2

    None, Unspecified: Not a stereo bond or unspecified

    Z, cis: Similar groups on same side of double bond
    E, trans: Similar groups on different side of double bond

    CisOrTrans, DoubleCross: cis or trans

=item B<StringifyBond>

    $BondString = $Bond->StringifyBond();

Returns a string containing information about I<bond> object.

=item B<SwitchBondFromAndToAtoms>

    $Bond = $Bond->SwitchBondFromAndToAtoms();

Swaps bond from and to atoms in I<Bond> and returns I<Bond>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Atom.pm, Molecule.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
