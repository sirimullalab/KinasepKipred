package Atom;
#
# File: Atom.pm
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
use PeriodicTable;
use Vector;
use MathUtil;
use Text::ParseWords;
use TextUtil;
use FileUtil;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(ObjectProperty Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, $ObjectID, %MDLValenceModelDataMap, %DaylightValenceModelDataMap);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyAtom';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializeAtom();

  $This->_InitializeAtomProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
sub _InitializeAtom {
  my($This) = @_;
  my($ObjectID) = _GetNewObjectID();

  # All other property names and values along with all Set/Get<PropertyName> methods
  # are implemented on-demand using ObjectProperty class.
  $This->{ID} = $ObjectID;
  $This->{Name} = "Atom ${ObjectID}";
  $This->{AtomSymbol} = '';
  $This->{AtomicNumber} = 0;
  $This->{XYZ} = Vector::ZeroVector;
}

# Initialize atom properties...
sub _InitializeAtomProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }
  if (!exists $NamesAndValues{'AtomSymbol'}) {
    carp "Warning: ${ClassName}->new: Atom object instantiated without setting atom symbol...";
  }

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # ID to keep track of objects...
  $ObjectID = 0;

  # Load atom class data...
  _LoadAtomClassData();
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

# Assign atom to  molecule...
sub _SetMolecule {
  my($This, $Molecule) = @_;

  $This->{Molecule} = $Molecule;

  # Weaken the reference to disable increment of reference count; otherwise,
  # it it becomes a circular reference and destruction of Molecule object doesn't
  # get initiated which in turn disables destruction of atom object.
  #
  Scalar::Util::weaken($This->{Molecule});

  return $This;
}

# Setup atom symbol and atomic number for the element...
#
# Possible atom symbol values:
#    . An element symbol or some other type of atom: L - Atom list; LP - Lone pair; R# - R group;
#       A, Q, * - unknown atom; or something else?
#
# Default mass number corresponds to the most abundant natural isotope unless it's explicity
# set using "MassNumber" property.
#
sub SetAtomSymbol {
  my($This, $AtomSymbol) = @_;
  my($AtomicNumber);

  $This->{AtomSymbol} = $AtomSymbol;

  $AtomicNumber = PeriodicTable::GetElementAtomicNumber($AtomSymbol);
  $This->{AtomicNumber} = (defined $AtomicNumber) ? $AtomicNumber : 0;

  return $This;
}

# Setup atom symbol and atomic number for the element...
sub SetAtomicNumber {
  my($This, $AtomicNumber) = @_;
  my($AtomSymbol);

  $AtomSymbol = PeriodicTable::GetElementAtomSymbol($AtomicNumber);
  if (!defined $AtomSymbol) {
    carp "Warning: ${ClassName}->SetAtomicNumber: Didn't set atomic number: Invalid atomic number, $AtomicNumber, specified...";
    return;
  }
  $This->{AtomicNumber} = $AtomicNumber;
  $This->{AtomSymbol} = $AtomSymbol;

  return $This;
}

# Set atom as stereo center...
#
sub SetStereoCenter {
  my($This, $StereoCenter) = @_;

  $This->SetProperty('StereoCenter', $StereoCenter);

  return $This;
}

# Is it a stereo center?
#
sub IsStereoCenter {
  my($This) = @_;
  my($StereoCenter);

  $StereoCenter = $This->GetProperty('StereoCenter');

  return (defined($StereoCenter) && $StereoCenter) ? 1 : 0;
}

# Set atom stereochemistry.
#
# Supported values are: R, S.
#
# Notes:
#
# . After the ligands around a central stereocenter has been ranked using CIP priority scheme and
# the lowest ranked ligand lies behind the center atom, then R and S values correspond to:
#
# R: Clockwise arrangement of remaining ligands around the central atom going from highest to lowest ranked ligand
# S: CounterClockwise arrangement of remaining ligands around the central atom going from highest to lowest ranked ligand
#
# . Assignment of any other arbitray values besides R and S is also allowed; however, a warning is printed.
#
sub SetStereochemistry {
  my($This, $Stereochemistry) = @_;

  if ($Stereochemistry !~ /^(R|S)$/i) {
    carp "Warning: ${ClassName}->SetStereochemistry: Assigning non-supported Stereochemistry value of $Stereochemistry. Supported values: R, S...";
  }

  $This->SetProperty('StereoCenter', 1);
  $This->SetProperty('Stereochemistry', $Stereochemistry);

  return $This;
}

# Setup mass number for atom...
sub SetMassNumber {
  my($This, $MassNumber) = @_;
  my($AtomicNumber, $AtomSymbol);

  $AtomicNumber = $This->{AtomicNumber};
  $AtomSymbol = $This->{AtomSymbol};
  if (!$AtomicNumber) {
    carp "Warning: ${ClassName}->SetMassNumber: Didn't set mass number: Non standard atom with atomic number, $AtomicNumber, and atomic symbol, $AtomSymbol...";
    return;
  }
  if (!PeriodicTable::IsElementNaturalIsotopeMassNumber($AtomicNumber, $MassNumber)) {
    carp "Warning: ${ClassName}->SetMassNumber: Unknown mass number, $MassNumber, specified for atom with atomic number, $AtomicNumber, and atomic symbol, $AtomSymbol. Don't forget to Set ExactMass property explicitly; otherwise, GetExactMass method would return mass of most abundant isotope...";
  }
  $This->SetProperty('MassNumber', $MassNumber);

  return $This;
}

# Get mass number...
#
sub GetMassNumber {
  my($This) = @_;

  # Is mass number explicity set?
  if ($This->HasProperty('MassNumber')) {
    return $This->GetProperty('MassNumber');
  }

  # Is it an element symbol?
  my($AtomicNumber) = $This->{AtomicNumber};
  if (!$AtomicNumber) {
    return 0;
  }

  # Return most abundant mass number...
  return PeriodicTable::GetElementMostAbundantNaturalIsotopeMassNumber($AtomicNumber);
}

# Get atomic weight:
#   . Explicitly set by the caller
#   . Using atomic number
#
sub GetAtomicWeight {
  my($This) = @_;

  # Is atomic weight explicity set?
  if ($This->HasProperty('AtomicWeight')) {
    return $This->GetProperty('AtomicWeight');
  }

  # Is it an element symbol?
  my($AtomicNumber) = $This->{AtomicNumber};
  if (!$AtomicNumber) {
    return 0;
  }

  # Return its atomic weight...
  return PeriodicTable::GetElementAtomicWeight($AtomicNumber);
}

# Get exact mass weight:
#   . Explicitly set by the caller
#   . Using atomic number and mass number explicity set by the caller
#   . Using atomic number and most abundant isotope
#
sub GetExactMass {
  my($This) = @_;

  # Is exact mass explicity set?
  if ($This->HasProperty('ExactMass')) {
    return $This->GetProperty('ExactMass');
  }

  # Is it an element symbol?
  my($AtomicNumber) = $This->{AtomicNumber};
  if (!$AtomicNumber) {
    return 0;
  }

  # Is mass number explicitly set?
  if ($This->HasProperty('MassNumber')) {
    my($MassNumber) = $This->GetProperty('MassNumber');
    if (PeriodicTable::IsElementNaturalIsotopeMassNumber($AtomicNumber, $MassNumber)) {
      return PeriodicTable::GetElementNaturalIsotopeMass($AtomicNumber, $MassNumber);
    }
  }

  # Return most abundant isotope mass...
  return PeriodicTable::GetElementMostAbundantNaturalIsotopeMass($AtomicNumber);
}

# Get formal charge:
#   . Explicitly set by the caller
#   . Or return zero insetad of undef
#
sub GetFormalCharge {
  my($This) = @_;
  my($FormalCharge);

  $FormalCharge = 0;
  if ($This->HasProperty('FormalCharge')) {
    $FormalCharge = $This->GetProperty('FormalCharge');
  }

  return defined($FormalCharge) ? $FormalCharge : 0;
}

# Get spin multiplicity:
#   . Explicitly set by the caller
#   . From FreeRadicalElectrons value explicitly set by the caller
#   . Or return zero insetad of undef
#
sub GetSpinMultiplicity {
  my($This) = @_;
  my($SpinMultiplicity);

  $SpinMultiplicity = 0;
  if ($This->HasProperty('SpinMultiplicity')) {
    $SpinMultiplicity = $This->GetProperty('SpinMultiplicity');
    return defined($SpinMultiplicity) ? $SpinMultiplicity : 0;
  }

  if ($This->HasProperty('FreeRadicalElectrons')) {
    my($FreeRadicalElectrons);
    $FreeRadicalElectrons = $This->GetProperty('FreeRadicalElectrons');

    SPINMULTIPLICITY: {
      if ($FreeRadicalElectrons == 1) { $SpinMultiplicity = 2; last SPINMULTIPLICITY;}
      if ($FreeRadicalElectrons == 2) { $SpinMultiplicity = 1; last SPINMULTIPLICITY;}
      carp "Warning: ${ClassName}->GetSpinMultiplicity: It's not possible to determine spin multiplicity from the specified free radical electrons value, $FreeRadicalElectrons. It has been set to 0...";
      $SpinMultiplicity = 0;
    }
  }

  return $SpinMultiplicity;
}

# Get number of free radical electrons:
#   . Explicitly set by the caller
#   . From SpinMultiplicity value explicitly set by the caller
#   . Or return zero insetad of undef
#
# Notes:
#  . For atoms with explicit assignment of SpinMultiplicity property values corresponding to
#    Singlet (two unpaired electrons corresponding to one spin state), Doublet (free radical; an unpaired
#    electron corresponding to two spin states), and Triplet (two unparied electrons corresponding to
#    three spin states; divalent carbon atoms (carbenes)), FreeRadicalElectrons are calculated as follows:
#
#       SpinMultiplicity: Doublet(2); FreeRadicalElectrons: 1
#       SpinMultiplicity: Singlet(1)/Triplet(3); FreeRadicalElectrons: 2
#
sub GetFreeRadicalElectrons {
  my($This) = @_;
  my($FreeRadicalElectrons);

  $FreeRadicalElectrons = 0;

  if ($This->HasProperty('FreeRadicalElectrons')) {
    $FreeRadicalElectrons = $This->GetProperty('FreeRadicalElectrons');
    return defined($FreeRadicalElectrons) ? $FreeRadicalElectrons : 0;
  }

  if ($This->HasProperty('SpinMultiplicity')) {
    my($SpinMultiplicity);
    $SpinMultiplicity = $This->GetProperty('SpinMultiplicity');

    SPINMULTIPLICITY: {
      if ($SpinMultiplicity == 1) { $FreeRadicalElectrons = 2; last SPINMULTIPLICITY;}
      if ($SpinMultiplicity == 2) { $FreeRadicalElectrons = 1; last SPINMULTIPLICITY;}
      if ($SpinMultiplicity == 3) { $FreeRadicalElectrons = 2; last SPINMULTIPLICITY;}
      carp "Warning: ${ClassName}->GetFreeRadicalElectrons: It's not possible to determine free radical electrons from the specified spin multiplicity value, $FreeRadicalElectrons. It has been set to 0...";
      $FreeRadicalElectrons = 0;
    }
  }

  return $FreeRadicalElectrons;
}

# Set atom coordinates using:
# . An array reference with three values
# . An array containg three values
# . A 3D vector
#
sub SetXYZ {
  my($This, @Values) = @_;

  if (!@Values) {
    carp "Warning: ${ClassName}->SetXYZ: No values specified...";
    return;
  }

  $This->{XYZ}->SetXYZ(@Values);
  return $This;
}

# Set X value...
sub SetX {
  my($This, $Value) = @_;

  if (!defined $Value) {
    carp "Warning: ${ClassName}->SetX: Undefined X value...";
    return;
  }
  $This->{XYZ}->SetX($Value);
  return $This;
}

# Set Y value...
sub SetY {
  my($This, $Value) = @_;

  if (!defined $Value) {
    carp "Warning: ${ClassName}->SetY: Undefined Y value...";
    return;
  }
  $This->{XYZ}->SetY($Value);
  return $This;
}

# Set Z value...
sub SetZ {
  my($This, $Value) = @_;

  if (!defined $Value) {
    carp "Warning: ${ClassName}->SetZ: Undefined Z value...";
    return;
  }
  $This->{XYZ}->SetZ($Value);
  return $This;
}

# Return XYZ as:
# . Reference to an array
# . An array
#
sub GetXYZ {
  my($This) = @_;

  return $This->{XYZ}->GetXYZ();
}

# Return XYZ as a vector object...
#
sub GetXYZVector {
  my($This) = @_;

  return $This->{XYZ};
}

# Get X value...
sub GetX {
  my($This) = @_;

  return $This->{XYZ}->GetX();
}

# Get Y value...
sub GetY {
  my($This) = @_;

  return $This->{XYZ}->GetY();
}

# Get Z value...
sub GetZ {
  my($This) = @_;

  return $This->{XYZ}->GetZ();
}

# Delete atom...
sub DeleteAtom {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    # Nothing to do...
    return $This;
  }
  my($Molecule) = $This->GetProperty('Molecule');

  return $Molecule->_DeleteAtom($This);
}

# Get atom neighbor objects as array. In scalar conetxt, return number of neighbors...
sub GetNeighbors {
  my($This, @ExcludeNeighbors) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule) = $This->GetProperty('Molecule');

  if (@ExcludeNeighbors) {
    return $This->_GetAtomNeighbors(@ExcludeNeighbors);
  }
  else {
    return $This->_GetAtomNeighbors();
  }
}

# Get atom neighbor objects as array. In scalar conetxt, return number of neighbors...
sub _GetAtomNeighbors {
  my($This, @ExcludeNeighbors) = @_;
  my($Molecule) = $This->GetProperty('Molecule');

  if (!@ExcludeNeighbors) {
    return $Molecule->_GetAtomNeighbors($This);
  }

  # Setup a map for neigbhors to exclude...
  my($ExcludeNeighbor, $ExcludeNeighborID, %ExcludeNeighborsIDsMap);

  %ExcludeNeighborsIDsMap = ();
  for $ExcludeNeighbor (@ExcludeNeighbors) {
    $ExcludeNeighborID = $ExcludeNeighbor->GetID();
    $ExcludeNeighborsIDsMap{$ExcludeNeighborID} = $ExcludeNeighborID;
  }

  # Generate a filtered neighbors list...
  my($Neighbor, $NeighborID, @FilteredAtomNeighbors);
  @FilteredAtomNeighbors = ();
  NEIGHBOR: for $Neighbor ($Molecule->_GetAtomNeighbors($This)) {
      $NeighborID = $Neighbor->GetID();
      if (exists $ExcludeNeighborsIDsMap{$NeighborID}) {
	next NEIGHBOR;
      }
    push @FilteredAtomNeighbors, $Neighbor;
  }

  return wantarray ? @FilteredAtomNeighbors : scalar @FilteredAtomNeighbors;
}

# Get specific atom neighbor objects as array. In scalar conetxt, return number of neighbors.
#
# Notes:
#   . AtomSpecification correspond to any valid AtomicInvariant based atomic specifications
#     as implemented in DoesAtomNeighborhoodMatch method.
#   . Multiple atom specifications can be used in a string delimited by comma.
#
sub GetNeighborsUsingAtomSpecification {
  my($This, $AtomSpecification, @ExcludeNeighbors) = @_;
  my(@AtomNeighbors);

  @AtomNeighbors = ();
  @AtomNeighbors = $This->GetNeighbors(@ExcludeNeighbors);

  # Does atom has any neighbors and do they need to be filtered?
  if (!(@AtomNeighbors && defined($AtomSpecification) && $AtomSpecification)) {
    return wantarray ? @AtomNeighbors : scalar @AtomNeighbors;
  }

  # Filter neighbors using atom specification...
  my($AtomNeighbor, @FilteredAtomNeighbors);

  @FilteredAtomNeighbors = ();
  NEIGHBOR: for $AtomNeighbor (@AtomNeighbors) {
    if (!$AtomNeighbor->_DoesAtomSpecificationMatch($AtomSpecification)) {
      next NEIGHBOR;
    }
    push @FilteredAtomNeighbors, $AtomNeighbor;
  }

  return wantarray ? @FilteredAtomNeighbors : scalar @FilteredAtomNeighbors;
}


# Get non-hydrogen atom neighbor objects as array. In scalar context, return number of neighbors...
sub GetHeavyAtomNeighbors {
  my($This) = @_;

  return $This->GetNonHydrogenAtomNeighbors();
}

# Get non-hydrogen atom neighbor objects as array. In scalar context, return number of neighbors...
sub GetNonHydrogenAtomNeighbors {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($NonHydrogenAtomsOnly, $HydrogenAtomsOnly) = (1, 0);

  return $This->_GetFilteredAtomNeighbors($NonHydrogenAtomsOnly, $HydrogenAtomsOnly);
}

# Get hydrogen atom neighbor objects as array. In scalar context, return numbe of neighbors...
sub GetHydrogenAtomNeighbors {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($NonHydrogenAtomsOnly, $HydrogenAtomsOnly) = (0, 1);

  return $This->_GetFilteredAtomNeighbors($NonHydrogenAtomsOnly, $HydrogenAtomsOnly);
}

# Get non-hydrogen neighbor of hydrogen atom...
#
sub GetNonHydrogenNeighborOfHydrogenAtom {
  my($This) = @_;

  # Is it Hydrogen?
  if (!$This->IsHydrogen()) {
    return undef;
  }
  my(@Neighbors);

  @Neighbors = $This->GetNonHydrogenAtomNeighbors();

  return (@Neighbors == 1) ? $Neighbors[0] : undef;
}

# Get filtered atom atom neighbors
sub _GetFilteredAtomNeighbors {
  my($This, $NonHydrogenAtomsOnly, $HydrogenAtomsOnly) = @_;

  # Check flags...
  if (!defined $NonHydrogenAtomsOnly) {
    $NonHydrogenAtomsOnly = 0;
  }
  if (!defined $HydrogenAtomsOnly) {
    $HydrogenAtomsOnly = 0;
  }
  my($Neighbor, @FilteredAtomNeighbors);

  @FilteredAtomNeighbors = ();
  NEIGHBOR: for $Neighbor ($This->GetNeighbors()) {
    if ($NonHydrogenAtomsOnly && $Neighbor->IsHydrogen()) {
      next NEIGHBOR;
    }
    if ($HydrogenAtomsOnly && (!$Neighbor->IsHydrogen())) {
      next NEIGHBOR;
    }
    push @FilteredAtomNeighbors, $Neighbor;
  }

  return wantarray ? @FilteredAtomNeighbors : scalar @FilteredAtomNeighbors;
}

# Get number of neighbors...
#
sub GetNumOfNeighbors {
  my($This) = @_;
  my($NumOfNeighbors);

  $NumOfNeighbors = $This->GetNeighbors();

  return (defined $NumOfNeighbors) ? $NumOfNeighbors : undef;
}

# Get number of neighbors which are non-hydrogen atoms...
sub GetNumOfHeavyAtomNeighbors {
  my($This) = @_;

  return $This->GetNumOfNonHydrogenAtomNeighbors();
}

# Get number of neighbors which are non-hydrogen atoms...
sub GetNumOfNonHydrogenAtomNeighbors {
  my($This) = @_;
  my($NumOfNeighbors);

  $NumOfNeighbors = $This->GetNonHydrogenAtomNeighbors();

  return (defined $NumOfNeighbors) ? $NumOfNeighbors : undef;
}

# Get number of neighbors which are hydrogen atoms...
sub GetNumOfHydrogenAtomNeighbors {
  my($This) = @_;
  my($NumOfNeighbors);

  $NumOfNeighbors = $This->GetHydrogenAtomNeighbors();

  return (defined $NumOfNeighbors) ? $NumOfNeighbors : undef;
}

# Get bond objects as array. In scalar context, return number of bonds...
sub GetBonds {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule) = $This->GetProperty('Molecule');

  return $Molecule->_GetAtomBonds($This);
}

# Get bond to specified atom...
sub GetBondToAtom {
  my($This, $Other) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule) = $This->GetProperty('Molecule');

  return $Molecule->_GetBondToAtom($This, $Other);
}

# It it bonded to a specified atom?
sub IsBondedToAtom {
  my($This, $Other) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule) = $This->GetProperty('Molecule');

  return $Molecule->_IsBondedToAtom($This, $Other);
}

# Get bond objects to non-hydrogen atoms as array. In scalar context, return number of bonds...
sub GetBondsToHeavyAtoms {
  my($This) = @_;

  return $This->GetBondsToNonHydrogenAtoms();
}

# Get bond objects to non-hydrogen atoms as array. In scalar context, return number of bonds...
sub GetBondsToNonHydrogenAtoms {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($BondsToNonHydrogenAtomsOnly, $BondsToHydrogenAtomsOnly) = (1, 0);

  return $This->_GetFilteredBonds($BondsToNonHydrogenAtomsOnly, $BondsToHydrogenAtomsOnly);
}

# Get bond objects to hydrogen atoms as array. In scalar context, return number of bonds...
sub GetBondsToHydrogenAtoms {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($BondsToNonHydrogenAtomsOnly, $BondsToHydrogenAtomsOnly) = (0, 1);

  return $This->_GetFilteredBonds($BondsToNonHydrogenAtomsOnly, $BondsToHydrogenAtomsOnly);
}

# Get filtered bonds...
sub _GetFilteredBonds {
  my($This, $BondsToNonHydrogenAtomsOnly, $BondsToHydrogenAtomsOnly) = @_;

  # Check flags...
  if (!defined $BondsToNonHydrogenAtomsOnly) {
    $BondsToNonHydrogenAtomsOnly = 0;
  }
  if (!defined $BondsToHydrogenAtomsOnly) {
    $BondsToHydrogenAtomsOnly = 0;
  }

  my($Bond, $BondedAtom, @FilteredBonds);

  @FilteredBonds = ();
  BOND: for $Bond ($This->GetBonds()) {
    $BondedAtom = $Bond->GetBondedAtom($This);
    if ($BondsToNonHydrogenAtomsOnly && $BondedAtom->IsHydrogen()) {
      next BOND;
    }
    if ($BondsToHydrogenAtomsOnly && (!$BondedAtom->IsHydrogen())) {
      next BOND;
    }
    push @FilteredBonds, $Bond;
  }

  return wantarray ? @FilteredBonds : (scalar @FilteredBonds);
}

# Get number of bonds...
#
sub GetNumOfBonds {
  my($This) = @_;
  my($NumOfBonds);

  $NumOfBonds = $This->GetBonds();

  return (defined $NumOfBonds) ? ($NumOfBonds) : undef;
}

# Get number of bonds to non-hydrogen atoms...
sub GetNumOfBondsToHeavyAtoms {
  my($This) = @_;

  return $This->GetNumOfBondsToNonHydrogenAtoms();
}

# Get number of bonds to non-hydrogen atoms...
sub GetNumOfBondsToNonHydrogenAtoms {
  my($This) = @_;
  my($NumOfBonds);

  $NumOfBonds = $This->GetBondsToNonHydrogenAtoms();

  return (defined $NumOfBonds) ? ($NumOfBonds) : undef;
}

# Get number of single bonds to heavy atoms...
sub GetNumOfSingleBondsToHeavyAtoms {
  my($This) = @_;

  return $This->GetNumOfSingleBondsToNonHydrogenAtoms();
}

# Get number of single bonds to non-hydrogen atoms...
sub GetNumOfSingleBondsToNonHydrogenAtoms {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  return $This->_GetNumOfBondsWithSpecifiedBondOrderToNonHydrogenAtoms(1);
}

# Get number of double bonds to heavy atoms...
sub GetNumOfDoubleBondsToHeavyAtoms {
  my($This) = @_;

  return $This->GetNumOfDoubleBondsToNonHydrogenAtoms();
}

# Get number of double bonds to non-hydrogen atoms...
sub GetNumOfDoubleBondsToNonHydrogenAtoms {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  return $This->_GetNumOfBondsWithSpecifiedBondOrderToNonHydrogenAtoms(2);
}

# Get number of triple bonds to heavy atoms...
sub GetNumOfTripleBondsToHeavyAtoms {
  my($This) = @_;

  return $This->GetNumOfTripleBondsToNonHydrogenAtoms();
}

# Get number of triple bonds to non-hydrogen atoms...
sub GetNumOfTripleBondsToNonHydrogenAtoms {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  return $This->_GetNumOfBondsWithSpecifiedBondOrderToNonHydrogenAtoms(3);
}

# Get number of bonds of specified bond order to non-hydrogen atoms...
sub _GetNumOfBondsWithSpecifiedBondOrderToNonHydrogenAtoms {
  my($This, $SpecifiedBondOrder) = @_;
  my($NumOfBonds, $Bond, $BondOrder, @Bonds);

  $NumOfBonds = 0;
  @Bonds = $This->GetBondsToNonHydrogenAtoms();
  for $Bond (@Bonds) {
    $BondOrder = $Bond->GetBondOrder();
    if ($SpecifiedBondOrder == $BondOrder) {
      $NumOfBonds++;
    }
  }
  return $NumOfBonds;
}

# Get number of aromatic bonds to heavy atoms...
sub GetNumOfAromaticBondsToHeavyAtoms {
  my($This) = @_;

  return $This->GetNumOfAromaticBondsToNonHydrogenAtoms();
}

# Get number of aromatic bonds to non-hydrogen atoms...
sub GetNumOfAromaticBondsToNonHydrogenAtoms {
  my($This) = @_;
  my($NumOfBonds, $Bond, @Bonds);

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }

  $NumOfBonds = 0;
  @Bonds = $This->GetBondsToNonHydrogenAtoms();
  for $Bond (@Bonds) {
    if ($Bond->IsAromatic()) { $NumOfBonds++; }
  }
  return $NumOfBonds;
}

# Get number of different bond types to non-hydrogen atoms...
#
sub GetNumOfBondTypesToHeavyAtoms {
  my($This, $CountAromaticBonds) = @_;

  return $This->GetNumOfBondTypesToNonHydrogenAtoms($CountAromaticBonds);
}

# Get number of single, double, triple, and aromatic bonds from an atom to all other
# non-hydrogen atoms. Value of CountAtomaticBonds parameter controls whether
# number of aromatic bonds is returned; default is not to count aromatic bonds. During
# counting of aromatic bonds, the bond marked aromatic is not included in the count
# of other bond types.
#
sub GetNumOfBondTypesToNonHydrogenAtoms {
  my($This, $CountAromaticBonds) = @_;
  my($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds, $None, $Bond, @Bonds);

  $CountAromaticBonds = defined($CountAromaticBonds) ? $CountAromaticBonds : 0;

  ($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds) = ('0') x 3;
  $NumOfAromaticBonds = $CountAromaticBonds ? 0 : undef;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return ($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds);
  }

  @Bonds = $This->GetBondsToNonHydrogenAtoms();

  for $Bond (@Bonds) {
    BONDTYPE: {
      if ($CountAromaticBonds) {
	if ($Bond->IsAromatic()) { $NumOfAromaticBonds++; last BONDTYPE; }
      }
      if ($Bond->IsSingle()) { $NumOfSingleBonds++; last BONDTYPE; }
      if ($Bond->IsDouble()) { $NumOfDoubleBonds++; last BONDTYPE; }
      if ($Bond->IsTriple()) { $NumOfTripleBonds++; last BONDTYPE; }
      $None = 1;
    }
  }
  return ($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds);
}

# Get number of sigma and pi bonds to heavy atoms...
#
sub GetNumOfSigmaAndPiBondsToHeavyAtoms {
  my($This) = @_;

  return $This->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
}

# Get number of sigma and pi bonds from an atom to all other non-hydrogen atoms.
# Sigma and pi bonds are counted using the following methodology: a single bond
# correspond to one sigma bond; a double bond contributes one to sigma bond count
# and one to pi bond count; a triple bond contributes one to sigma bond count and
# two to pi bond count.
#
sub GetNumOfSigmaAndPiBondsToNonHydrogenAtoms {
  my($This) = @_;
  my($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfSigmaBonds, $NumOfPiBonds);

  ($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds) = $This->GetNumOfBondTypesToNonHydrogenAtoms();

  $NumOfSigmaBonds = $NumOfSingleBonds + $NumOfDoubleBonds + $NumOfTripleBonds;
  $NumOfPiBonds = $NumOfDoubleBonds + 2*$NumOfTripleBonds;

  return ($NumOfSigmaBonds, $NumOfPiBonds);
}

# Get information related to atoms for all heavy atoms attached to an atom..
#
sub GetHeavyAtomNeighborsAtomInformation {
  my($This) = @_;

  return $This->GetNonHydrogenAtomNeighborsAtomInformation();
}

# Get information related to atoms for all non-hydrogen atoms attached to an atom..
#
# The following values are returned:
#   . Number of non-hydrogen atom neighbors
#   . A reference to an array containing atom objects correpsonding to non-hydrogen
#     atom neighbors
#   . Number of different types of non-hydrogen atom neighbors
#   . A reference to a hash containing atom symbol as key with value corresponding
#     to its count for non-hydrogen atom neighbors
#
sub GetNonHydrogenAtomNeighborsAtomInformation {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return (undef, undef, undef, undef);
  }
  my($AtomSymbol, $AtomNeighbor, $NumOfAtomNeighbors, $NumOfAtomNeighborsType, @AtomNeighbors, %AtomNeighborsTypeMap);

  $NumOfAtomNeighbors = 0; @AtomNeighbors = ();
  $NumOfAtomNeighborsType = 0; %AtomNeighborsTypeMap = ();

  @AtomNeighbors = $This->GetNonHydrogenAtomNeighbors();
  $NumOfAtomNeighbors = scalar @AtomNeighbors;

  for $AtomNeighbor (@AtomNeighbors) {
    $AtomSymbol = $AtomNeighbor->{AtomSymbol};
    if (exists $AtomNeighborsTypeMap{$AtomSymbol}) {
      $AtomNeighborsTypeMap{$AtomSymbol} += 1;
    }
    else {
      $AtomNeighborsTypeMap{$AtomSymbol} = 1;
      $NumOfAtomNeighborsType++;
    }
  }

  return ($NumOfAtomNeighbors, \@AtomNeighbors, $NumOfAtomNeighborsType, \%AtomNeighborsTypeMap);
}

# Get information related to bonds for all heavy atoms attached to an atom..
#
sub GetHeavyAtomNeighborsBondformation {
  my($This) = @_;

  return $This->GetNonHydrogenAtomNeighborsBondInformation();
}

# Get information related to bonds for all non-hydrogen atoms attached to an atom..
#
# The following values are returned:
#   . Number of bonds to non-hydrogen atom neighbors
#   . A reference to an array containing bond objects correpsonding to non-hydrogen
#     atom neighbors
#   . A reference to a hash containing bond type as key with value corresponding
#     to its count for non-hydrogen atom neighbors. Bond types are: Single, Double or Triple
#   . A reference to a hash containing atom symbol as key pointing to bond type as second
#     key with values correponding to count of bond types for atom symbol for non-hydrogen
#     atom neighbors
#   . A reference to a hash containing atom symbol as key pointing to bond type as second
#     key with values correponding to atom objects array involved in corresponding bond type for
#     atom symbol for non-hydrogen atom neighbors
#
sub GetNonHydrogenAtomNeighborsBondInformation {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return (undef, undef, undef, undef, undef);
  }
  my($BondedAtom, $BondedAtomSymbol, $BondType, $None, $Bond, $NumOfBonds, @Bonds, %BondTypeCountMap, %AtomsBondTypesCountMap, %AtomsBondTypeAtomsMap);

  $NumOfBonds = 0; @Bonds = ();
  %BondTypeCountMap = ();
  %AtomsBondTypesCountMap = (); %AtomsBondTypeAtomsMap = ();

  $BondTypeCountMap{Single} = 0;
  $BondTypeCountMap{Double} = 0;
  $BondTypeCountMap{Triple} = 0;

  @Bonds = $This->GetBondsToNonHydrogenAtoms();
  $NumOfBonds = scalar @Bonds;

  BOND: for $Bond (@Bonds) {
    $BondType = $Bond->IsSingle() ? "Single" : ($Bond->IsDouble() ? "Double" : ($Bond->IsTriple() ? "Triple" : ""));
    if (!$BondType) {
      next BOND;
    }

    # Track bond types...
    if (exists $BondTypeCountMap{$BondType}) {
      $BondTypeCountMap{$BondType} += 1;
    }
    else {
      $BondTypeCountMap{$BondType} = 1;
    }

    $BondedAtom = $Bond->GetBondedAtom($This);
    $BondedAtomSymbol = $BondedAtom->{AtomSymbol};

    # Track bond types count for atom types involved in specific bond types...
    if (!exists $AtomsBondTypesCountMap{$BondedAtomSymbol}) {
      %{$AtomsBondTypesCountMap{$BondedAtomSymbol}} = ();
    }
    if (exists $AtomsBondTypesCountMap{$BondedAtomSymbol}{$BondType}) {
      $AtomsBondTypesCountMap{$BondedAtomSymbol}{$BondType} += 1;
    }
    else {
      $AtomsBondTypesCountMap{$BondedAtomSymbol}{$BondType} = 1;
    }

    # Track atoms involved in specific bond types for specific atom types...
    if (!exists $AtomsBondTypeAtomsMap{$BondedAtomSymbol}) {
      %{$AtomsBondTypeAtomsMap{$BondedAtomSymbol}} = ();
    }
    if (!exists $AtomsBondTypeAtomsMap{$BondedAtomSymbol}{$BondType}) {
      @{$AtomsBondTypeAtomsMap{$BondedAtomSymbol}{$BondType}} = ();
    }
    push @{$AtomsBondTypeAtomsMap{$BondedAtomSymbol}{$BondType}}, $BondedAtom;
  }

  return ($NumOfBonds, \@Bonds, \%BondTypeCountMap, \%AtomsBondTypesCountMap, \%AtomsBondTypeAtomsMap);
}

# Get number of bonds to hydrogen atoms...
sub GetNumOfBondsToHydrogenAtoms {
  my($This) = @_;
  my($NumOfBonds);

  $NumOfBonds = $This->GetBondsToHydrogenAtoms();

  return (defined $NumOfBonds) ? ($NumOfBonds) : undef;
}

# Get sum of bond orders to all bonded atoms...
#
sub GetSumOfBondOrders {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }

  return $This->_GetSumOfBondOrders();
}

# Get sum of bond orders to non-hydrogen atoms only...
#
sub GetSumOfBondOrdersToHeavyAtoms {
  my($This) = @_;

  return $This->GetSumOfBondOrdersToNonHydrogenAtoms();
}

# Get sum of bond orders to non-hydrogen atoms only...
#
sub GetSumOfBondOrdersToNonHydrogenAtoms {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($ToNonHydrogenAtomsOnly, $ToHydrogenAtomsOnly) = (1, 0);

  return $This->_GetSumOfBondOrders($ToNonHydrogenAtomsOnly, $ToHydrogenAtomsOnly);
}

# Get sum of bond orders to hydrogen atoms only...
#
sub GetSumOfBondOrdersToHydrogenAtoms {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($ToNonHydrogenAtomsOnly, $ToHydrogenAtomsOnly) = (0, 1);

  return $This->_GetSumOfBondOrders($ToNonHydrogenAtomsOnly, $ToHydrogenAtomsOnly);
}

# Get sum of bond orders to all bonded atoms,  non-hydrogen or hydrogen bonded atoms...
#
sub _GetSumOfBondOrders {
  my($This, $ToNonHydrogenAtomsOnly, $ToHydrogenAtomsOnly) = @_;

  # Check flags...
  if (!defined $ToNonHydrogenAtomsOnly) {
    $ToNonHydrogenAtomsOnly = 0;
  }
  if (!defined $ToHydrogenAtomsOnly) {
    $ToHydrogenAtomsOnly = 0;
  }
  my($Bond, $SumOfBondOrders, @Bonds);
  @Bonds = ();

  if ($ToNonHydrogenAtomsOnly) {
    @Bonds = $This->GetBondsToNonHydrogenAtoms();
  }
  elsif ($ToHydrogenAtomsOnly) {
    @Bonds = $This->GetBondsToHydrogenAtoms();
  }
  else {
    # All bonds...
    @Bonds = $This->GetBonds();
  }

  $SumOfBondOrders = 0;
  for $Bond (@Bonds) {
    $SumOfBondOrders += $Bond->GetBondOrder();
  }

  if ($SumOfBondOrders =~ /\./) {
    #
    # Change any fractional bond order to next largest integer...
    #
    # As long as aromatic bond orders in a ring are correctly using using 4n + 2 Huckel rule
    # (BondOrder: 1.5) or explicity set as Kekule  bonds (alternate single/double),
    # SumOfBondOrders should add up to an integer.
    #
    $SumOfBondOrders = ceil($SumOfBondOrders);
  }

  return $SumOfBondOrders;
}

# Get largest bond order to any bonded atoms...
#
sub GetLargestBondOrder {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }

  return $This->_GetLargestBondOrder();
}

# Get largest bond order to bonded non-hydrogen atoms...
#
sub GetLargestBondOrderToHeavyAtoms {
  my($This) = @_;

  return $This->GetLargestBondOrderToNonHydrogenAtoms();
}

# Get largest bond order to bonded non-hydrogen atoms...
#
sub GetLargestBondOrderToNonHydrogenAtoms {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }

  my($ToNonHydrogenAtomsOnly) = (1);

  return $This->_GetLargestBondOrder($ToNonHydrogenAtomsOnly);
}

# Get largest bond order to all bonded atoms, non-hydrogen or hydrogen bonded atoms...
#
sub _GetLargestBondOrder {
  my($This, $ToNonHydrogenAtomsOnly, $ToHydrogenAtomsOnly) = @_;

  # Check flags...
  if (!defined $ToNonHydrogenAtomsOnly) {
    $ToNonHydrogenAtomsOnly = 0;
  }
  if (!defined $ToHydrogenAtomsOnly) {
    $ToHydrogenAtomsOnly = 0;
  }
  my($Bond, $LargestBondOrder, $BondOrder, @Bonds);
  @Bonds = ();

  if ($ToNonHydrogenAtomsOnly) {
    @Bonds = $This->GetBondsToNonHydrogenAtoms();
  }
  elsif ($ToHydrogenAtomsOnly) {
    @Bonds = $This->GetBondsToHydrogenAtoms();
  }
  else {
    # All bonds...
    @Bonds = $This->GetBonds();
  }

  $LargestBondOrder = 0;
  for $Bond (@Bonds) {
    $BondOrder = $Bond->GetBondOrder();
    if ($BondOrder > $LargestBondOrder) {
      $LargestBondOrder = $BondOrder;
    }
  }

  return $LargestBondOrder;
}

# Get number of implicit hydrogen for atom...
#
sub GetImplicitHydrogens {
  my($This) = @_;

  return $This->GetNumOfImplicitHydrogens();
}

# Get number of implicit hydrogen for atom...
#
sub GetNumOfImplicitHydrogens {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }

  # Is ImplicitHydrogens property explicitly set?
  if ($This->HasProperty('ImplicitHydrogens')) {
    return $This->GetProperty('ImplicitHydrogens');
  }

  # Is it an element symbol?
  if (!$This->{AtomicNumber}) {
    return 0;
  }

  my($ImplicitHydrogens, $PotentialTotalValence, $SumOfBondOrders);

  $ImplicitHydrogens = 0;
  $SumOfBondOrders = $This->GetSumOfBondOrders();
  $PotentialTotalValence = $This->GetPotentialTotalCommonValence();

  if (defined($PotentialTotalValence) && defined($SumOfBondOrders)) {
    # Subtract sum of bond orders to non-hydrogen and hydrogen atom neighbors...
    $ImplicitHydrogens = $PotentialTotalValence - $SumOfBondOrders;
  }

  return $ImplicitHydrogens > 0 ? $ImplicitHydrogens : 0;
}

# Get number of bonds available to form additional bonds with heavy atoms, excluding
# any implicit bonds to hydrogens set using ImplicitHydrogens property.
#
# It's different from number of implicit or missing hydrogens, both of which are equivalent.
#
# For example, in a SMILES string, [nH] ring atom corresponds to an aromatic nitrogen.
# Although the hydrogen specified for n is treated internally as implicit hydrogen and shows
# up in missing hydrogen count, it's not available to participate in double bonds to additional
# heavy atoms.
#
sub GetNumOfBondsAvailableForHeavyAtoms {
  my($This) = @_;

  return $This->GetNumOfBondsAvailableForNonHydrogenAtoms();
}

# It's another name for GetNumOfBondsAvailableForHeavyAtoms
#
sub GetNumOfBondsAvailableForNonHydrogenAtoms {
  my($This) = @_;
  my($NumOfAvailableBonds, $PotentialTotalValence, $SumOfBondOrders);

  $NumOfAvailableBonds = 0;

  $SumOfBondOrders = $This->GetSumOfBondOrders();
  $PotentialTotalValence = $This->GetPotentialTotalCommonValence();

  if (defined($PotentialTotalValence) && defined($SumOfBondOrders)) {
    # Subtract sum of bond orders to non-hydrogen and hydrogen atom neighbors...
    $NumOfAvailableBonds = $PotentialTotalValence - $SumOfBondOrders;
  }

  if ($This->HasProperty('ImplicitHydrogens')) {
    $NumOfAvailableBonds -= $This->GetProperty('ImplicitHydrogens');
  }

  return $NumOfAvailableBonds > 0 ? $NumOfAvailableBonds : 0;
}

# Disable setting of explicit hydrogens property...
sub SetExplicitHydrogens {
  my($This, $Value) = @_;

  carp "Warning: ${ClassName}->SetExplicitHydrogens: Setting of explicit hydrogens is not supported...";

  return $This;
}

# Get number of explicit hydrogens for atom...
#
sub GetExplicitHydrogens {
  my($This) = @_;

  return $This->GetNumOfExplicitHydrogens();
}

# Get number of explicit hydrogens for atom...
#
sub GetNumOfExplicitHydrogens {
  my($This) = @_;
  my($HydrogenAtomNbrs);

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }

  $HydrogenAtomNbrs = $This->GetNumOfHydrogenAtomNeighbors();

  return defined $HydrogenAtomNbrs ? $HydrogenAtomNbrs : 0;
}

# Get num of missing hydrogens...
#
sub GetMissingHydrogens {
  my($This) = @_;

  return $This->GetNumOfMissingHydrogens();
}

# Get num of missing hydrogens...
#
sub GetNumOfMissingHydrogens {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }

  return $This->GetNumOfImplicitHydrogens();
}

# Get total number of hydrogens...
#
sub GetHydrogens {
  my($This) = @_;

  return $This->GetNumOfHydrogens();
}

# Get total number of hydrogens...
#
sub GetNumOfHydrogens {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }

  return $This->GetNumOfImplicitHydrogens() + $This->GetNumOfExplicitHydrogens();
}

# Valence corresponds to the number of electrons used by an atom in bonding:
#
#   Valence = ValenceElectrons - ValenceFreeElectrons = BondingElectrons
#
# Single, double, triple bonds with bond orders of 1, 2 and 3 correspond to contribution of
# 1, 2, and 3 bonding electrons. So Valence can be computed using:
#
#   Valence = SumOfBondOrders + NumOfMissingHydrogens + FormalCharge
#
# where positive and negative values of FormalCharge increase and decrease the number
# of bonding electrons respectively.
#
# The current release of MayaChemTools supports the following three valence models, which
# are used during calculation of implicit hydrogens: MDLValenceModel, DaylightValenceModel,
# InternalValenceModel or MayaChemToolsValenceModel.
#
# Notes:
#  . This doesn't always corresponds to explicit valence.
#  . Missing hydrogens are included in the valence.
#  . For neutral molecules, valence and sum of bond order are equal.
#  . For molecules containing only single bonds, SumOfBondOrders and NumOfBonds are equal.
#  . Free radical electrons lead to the decrease in valence. For atoms with explicit assignment
#    of SpinMultiplicity property values corresponding to Singlet (two unparied electrons
#    corresponding to one spin state), Doublet (free radical; an unpaired electron corresponding
#    to two spin states), and Triplet (two unparied electrons corresponding to three spin states;
#    divalent carbon atoms (carbenes)), FreeRadicalElectrons are calculated as follows:
#
#       SpinMultiplicity: Doublet(2); FreeRadicalElectrons: 1 (one valence electron not available for bonding)
#       SpinMultiplicity: Singlet(1)/Triplet(3); FreeRadicalElectrons: 2 (two valence electrons not available for bonding)
#
sub GetValence {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }

  # Is Valence property explicitly set?
  if ($This->HasProperty('Valence')) {
    return $This->GetProperty('Valence');
  }
  my($Valence);

  $Valence = $This->GetSumOfBondOrders() + $This->GetNumOfMissingHydrogens() + $This->GetFormalCharge();

  return $Valence > 0 ? $Valence : 0;
}

# Get free non-bodning valence electrons left on atom after taking into account
# sum of bond orders, missing hydrogens and formal charged on the atom. Free
# radical electrons are included in the valence free electrons count by default.
#
# Valence corresponds to number of electrons used by atom in bonding:
#
#   Valence = ValenceElectrons - ValenceFreeElectrons
#
# Additionally, valence can also be calculated by:
#
#   Valence = SumOfBondOrders + NumOfMissingHydrogens + FormalCharge
#
# Valence and SumOfBondOrders are equal for neutral molecules.
#
# From two formulas for Valence described above, non-bonding free electrons
# left can be computed by:
#
#  ValenceFreeElectrons = ValenceElectrons - Valence
#                       = ValenceElectrons - SumOfBondOrders -
#                         NumOfMissingHydrogens - FormalCharge
#
# . Notes:
#    . Missing hydrogens are excluded from the valence free electrons.
#    . Any free radical electrons are considered part of the valence free electrons
#      by default.
#
# Examples:
#
# o NH3: ValenceFreeElectrons = 5 - 3 = 5 - 3 - 0 - 0 = 2
# o NH2: ValenceFreeElectrons = 5 - 3 = 5 - 2 - 1 - 0 = 2
# o NH4+; ValenceFreeElectrons = 5 - 5 = 5 - 4 - 0 - 1 = 0
# o NH3+; ValenceFreeElectrons = 5 - 5 = 5 - 3 - 1 - 1 = 0
# o C(=O)O- : ValenceFreeElectrons on O- = 6 - 0 = 6 - 1 - 0 - (-1) = 6
# o C(=O)O- : ValenceFreeElectrons on =O = 6 - 2 = 6 - 2 - 0 - 0 = 4
#
#
sub GetValenceFreeElectrons {
  my($This, $ExcludeFreeRadicalElectrons) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }

  # Is ValenceFreeElectrons property explicitly set?
  if ($This->HasProperty('ValenceFreeElectrons')) {
    return $This->GetProperty('ValenceFreeElectrons');
  }

  if (!$This->{AtomicNumber}) {
    return 0;
  }

  my($ValenceFreeElectrons);

  $ValenceFreeElectrons = $This->GetValenceElectrons() - $This->GetValence();
  if ($ExcludeFreeRadicalElectrons) {
    $ValenceFreeElectrons -= $This->GetFreeRadicalElectrons();
  }

  return $ValenceFreeElectrons > 0 ? $ValenceFreeElectrons : 0;
}

# Get potential total common valence for calculating the number of implicit hydrogens
# using the specified common valence model or default internal model for a molecule...
#
sub GetPotentialTotalCommonValence {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($PotentialTotalValence, $ValenceModel);

  $PotentialTotalValence = 0;
  $ValenceModel = $This->GetProperty('Molecule')->GetValenceModel();

  VALENCEMODEL: {
    if ($ValenceModel =~ /^MDLValenceModel$/i) {
      $PotentialTotalValence = $This->_GetPotentialTotalCommonValenceUsingMDLValenceModel();
      last VALENCEMODEL;
    }
    if ($ValenceModel =~ /^DaylightValenceModel$/i) {
       $PotentialTotalValence = $This->_GetPotentialTotalCommonValenceUsingDaylightValenceModel();
      last VALENCEMODEL;
    }
    if ($ValenceModel !~ /^(InternalValenceModel|MayaChemToolsValenceModel)$/i) {
      carp "Warning: ${ClassName}->GetPotentialTotalCommonValence: The current release of MayaChemTools doesn't support the specified valence model $ValenceModel. Supported valence models: MDLValenceModel, DaylightValenceModel, InternalValenceModel. Using internal valence model...";
    }
    # Use internal valence model as the default valence model...
    $PotentialTotalValence = $This->_GetPotentialTotalCommonValenceUsingInternalValenceModel();
  }

  return $PotentialTotalValence;
}

# Get potential total common valence using data for MDL valence model available in file,
# lib/data/MDLValenceModelData.csv, distributed with the package...
#
sub _GetPotentialTotalCommonValenceUsingMDLValenceModel {
  my($This) = @_;

  return $This->_GetPotentialTotalCommonValenceUsingValenceModelData(\%MDLValenceModelDataMap);

}

# Get potential total common valence using data for Daylight valence model available in file,
# lib/data/DaylightValenceModelData.csv, distributed with the release...
#
sub _GetPotentialTotalCommonValenceUsingDaylightValenceModel {
  my($This) = @_;

  return $This->_GetPotentialTotalCommonValenceUsingValenceModelData(\%DaylightValenceModelDataMap);
}

# Get potential total common valence using data for a specific valence model...
#
sub _GetPotentialTotalCommonValenceUsingValenceModelData {
  my($This, $ValenceModelDataRef) = @_;
  my($AtomicNumber, $FormalCharge);

  $AtomicNumber = $This->{AtomicNumber};
  if (!$AtomicNumber) {
    return 0;
  }

  $FormalCharge = $This->GetFormalCharge();

  # Is any valence model data available for atomic number and formal charge?
  if (!exists $ValenceModelDataRef->{$AtomicNumber}) {
    return 0;
  }
  if (!exists $ValenceModelDataRef->{$AtomicNumber}{$FormalCharge}) {
    return 0;
  }

  my($PotentialTotalValence, $SumOfBondOrders, $CurrentEffectiveValence, $AvailableCommonValence);

  $SumOfBondOrders = $This->GetSumOfBondOrders();
  if (!defined $SumOfBondOrders) {
    $SumOfBondOrders = 0;
  }
  $CurrentEffectiveValence = $SumOfBondOrders + $This->GetFreeRadicalElectrons();

  $PotentialTotalValence = 0;
  VALENCE: for $AvailableCommonValence (@{$ValenceModelDataRef->{$AtomicNumber}{$FormalCharge}{CommonValences}}) {
      if ($CurrentEffectiveValence <= $AvailableCommonValence) {
	$PotentialTotalValence = $AvailableCommonValence;
	last VALENCE;
      }
  }

  return $PotentialTotalValence;
}

#
# For elements with one one common valence, potential total common valence used
# during the calculation for number of implicit hydrogens during InternalValenceMode
# corresponds to:
#
#   CommonValence + FormalCharge - FreeRadicalElectrons
#
# For elements with multiple common valences, each common valence is used to
# calculate total potential common valence as shown above, and the first total potential
# common valence gerater than the sum of bond orderes is selected as the final total
# common valence.
#
# Group numbers > 14 - Group numbers 15 (N), 16 (O), 17 (F), 18 (He)
#
# Formal charge sign is not adjusted. Positive and negative values result in the
# increase and decrease of valence.
#
# Group 14 containing C, Si, Ge, Sn, Pb...
#
# Formal charge sign is reversed for positive values. Both positive and negative
# values result in the decrease of valence.
#
# Group 13 containing B, Al, Ga, In, Tl...
#
# Formal charge sign is always reversed. Positive and negative values result in the
# decrease and increase of valence.
#
# Groups 1 (H) through 12 (Zn)...
#
# Formal charge sign is reversed for positive values. Both positive and negative
# values result in the decrease of valence.
#
# Lanthanides and actinides...
#
# Formal charge sign is reversed for positive values. Both positive and negative
# values result in the decrease of valence.
#
# Notes:
#  . CommonValence and HighestCommonValence available from PeriodicTable module
#    are equivalent to most common and highest sum of bond orders for an element. For
#    neutral atoms involved only in single bonds, it corresponds to highest number of
#    allowed bonds for the atom.
#  . FormalCharge sign is reversed for electropositive elements with positive formal charge
#    during common valence calculations. Electropositive elements, metals and transition elements,
#    have usually plus formal charge and it leads to decrease in common valence; the negative
#    formal charge should result in the decrease of common valence.
#  . For carbon, both plus/minus formal charge cause decrease in common valence
#  . For elements on the right of carbon in periodic table, electronegative elements, plus formal
#    charge causes common valence to increase and minus formal charge cause it to decrease.
#
sub _GetPotentialTotalCommonValenceUsingInternalValenceModel {
  my($This) = @_;
  my($AtomicNumber, $CommonValences);

  $AtomicNumber = $This->{AtomicNumber};
  if (!$AtomicNumber) {
    return 0;
  }

  $CommonValences = PeriodicTable::GetElementCommonValences($AtomicNumber);
  if (!$CommonValences) {
    return 0;
  }

  my($PotentialTotalValence, $AdjustedFormalCharge, $FreeRadicalElectrons, $SumOfBondOrders, $AvailableCommonValence, @AvailableCommonValences);

  $AdjustedFormalCharge = $This->_GetFormalChargeAdjustedForInternalValenceModel();
  $FreeRadicalElectrons = $This->GetFreeRadicalElectrons();

  $SumOfBondOrders = $This->GetSumOfBondOrders();
  if (!defined $SumOfBondOrders) {
    $SumOfBondOrders = 0;
  }

  @AvailableCommonValences = split /\,/, $CommonValences;

  if (@AvailableCommonValences == 1) {
    # Calculate potential total valence using the only available common valence...
    $PotentialTotalValence = $AvailableCommonValences[0] + $AdjustedFormalCharge - $FreeRadicalElectrons;
  }
  else {
    # Calculate potential total valence using common valence from a list of available valences
    # that makes it higher than sum of bond orders or using the highest common valence...
    VALENCE: for $AvailableCommonValence (@AvailableCommonValences) {
      $PotentialTotalValence = $AvailableCommonValence + $AdjustedFormalCharge - $FreeRadicalElectrons;

      if ($PotentialTotalValence < 0 || $PotentialTotalValence >= $SumOfBondOrders) {
	last VALENCE;
      }
    }
  }

  return $PotentialTotalValence > 0 ? $PotentialTotalValence : 0;
}

# Adjust sign of the formal charge for potential total common valence calculation
# used during internal valence model to figure out number of implicit hydrogens.
#
sub _GetFormalChargeAdjustedForInternalValenceModel {
  my($This) = @_;
  my($FormalCharge, $GroupNumber, $SwitchSign);

  $FormalCharge = $This->GetFormalCharge();
  if ($FormalCharge == 0) {
    return 0;
  }

  $GroupNumber = $This->GetGroupNumber();
  if (!defined $GroupNumber) {
    return $FormalCharge;
  }

  # Group numbers > 14 - Group numbers 15 (N), 16 (O), 17 (F), 18 (He)
  #
  # Formal charge sign is not adjusted. Positive and negative values result in the
  # increase and decrease of valence.
  #
  # Group 14 containing C, Si, Ge, Sn, Pb...
  #
  # Formal charge sign is reversed for positive values. Both positive and negative
  # values result in the decrease of valence.
  #
  # Group 13 containing B, Al, Ga, In, Tl...
  #
  # Formal charge sign is always reversed. Positive and negative values result in the
  # decrease and increase of valence.
  #
  # Groups 1 (H) through 12 (Zn)...
  #
  # Formal charge sign is reversed for positive values. Both positive and negative
  # values result in the decrease of valence.
  #
  # Lanthanides and actinides...
  #
  # Formal charge sign is reversed for positive values. Both positive and negative
  # values result in the decrease of valence.
  #

  $SwitchSign = 0;
  if (length $GroupNumber) {
    GROUPNUMBER: {
      if ($GroupNumber > 14) {
	# Groups on the right side of C group in the periodic table...
	$SwitchSign = 0;
	last GROUPNUMBER;
      }
      if ($GroupNumber == 14) {
	# Group containing C, Si, Ge, Sn, Pb...
	$SwitchSign = ($FormalCharge > 0) ? 1 : 0;
	last GROUPNUMBER;
      }
      if ($GroupNumber == 13) {
	# Group containing B, Al, Ga, In, Tl...
	$SwitchSign = 1;
	last GROUPNUMBER;
      }
      # Groups 1 (H) through 12 (Zn)...
      if ($GroupNumber >=1 && $GroupNumber <= 12) {
	# Groups 1 (H) through 12 (Zn)...
	$SwitchSign = ($FormalCharge > 0) ? 1 : 0;
	last GROUPNUMBER;
      }
    }
  }
  else {
    # Lanthanides and actinides...
    $SwitchSign = ($FormalCharge > 0) ? 1 : 0;
  }

  if ($SwitchSign) {
    $FormalCharge *= -1.0;
  }

  return $FormalCharge;
}

# Get lowest common valence...
sub GetLowestCommonValence {
  my($This) = @_;

  # Is LowestCommonValence property explicitly set?
  if ($This->HasProperty('LowestCommonValence')) {
    return $This->GetProperty('LowestCommonValence');
  }
  my($AtomicNumber, $LowestCommonValence);

  $AtomicNumber = $This->{AtomicNumber};
  if (!$AtomicNumber) {
    return 0;
  }
  # Any need to differentiate between internal and other valence models...

  # LowestCommonValence is not set for all elements...
  $LowestCommonValence = PeriodicTable::GetElementLowestCommonValence($AtomicNumber);
  if (!$LowestCommonValence) {
    $LowestCommonValence = undef;
  }

  return $LowestCommonValence;
}

# Get highest common valence...
sub GetHighestCommonValence {
  my($This) = @_;

  # Is HighestCommonValence property explicitly set?
  if ($This->HasProperty('HighestCommonValence')) {
    return $This->GetProperty('HighestCommonValence');
  }
  my($AtomicNumber, $HighestCommonValence);

  $AtomicNumber = $This->{AtomicNumber};
  if (!$AtomicNumber) {
    return 0;
  }

  # Any need to differentiate between internal and other valence models...

  # HighestCommonValence is not set for all elements...
  $HighestCommonValence = PeriodicTable::GetElementHighestCommonValence($AtomicNumber);
  if (!$HighestCommonValence) {
    $HighestCommonValence = undef;
  }

  return $HighestCommonValence;
}

# Get valence electrons...
sub GetValenceElectrons {
  my($This) = @_;

  # Is ValenceElectrons property explicitly set?
  if ($This->HasProperty('ValenceElectrons')) {
    return $This->GetProperty('ValenceElectrons');
  }
  my($AtomicNumber, $ValenceElectrons);

  $AtomicNumber = $This->{AtomicNumber};
  if (!$AtomicNumber) {
    return 0;
  }

  $ValenceElectrons = PeriodicTable::GetElementValenceElectrons($AtomicNumber);

  return $ValenceElectrons;
}

# Add hydrogens to specified atom in molecule and return number of hydrogens added:
#
#   o HydrogensToAdd = ImplicitHydrogenCount - ExplicitHydrogenCount
#
#   o XYZ are set to ZeroVector
#
sub AddHydrogens {
  my($This, $HydrogenPositionsWarning) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  if (!defined $HydrogenPositionsWarning) {
    $HydrogenPositionsWarning = 1;
  }
  if ($HydrogenPositionsWarning) {
    carp "Warning: ${ClassName}->AddHydrogens: The current release of MayaChemTools doesn't assign any hydrogen positions...";
  }

  # Is it an element symbol?
  if (!$This->{AtomicNumber}) {
    return 0;
  }

  my($Molecule, $HydrogensAdded, $HydrogensToAdd);

  $Molecule = $This->GetProperty('Molecule');
  $HydrogensAdded = 0;
  $HydrogensToAdd = $This->GetNumOfMissingHydrogens();
  if ($HydrogensToAdd <= 0) {
    return $HydrogensAdded;
  }

  my($Count, $Hydrogen);

  for $Count (1 .. $HydrogensToAdd) {
    $HydrogensAdded++;

    $Hydrogen = $Molecule->NewAtom('AtomSymbol' => 'H', 'XYZ' => [0, 0, 0]);
    $Molecule->NewBond('Atoms' => [$This, $Hydrogen], 'BondOrder' => 1);
  }

  return $HydrogensAdded;
}

# Delete hydrogens attached to atom in molecule and return total number of hydrogens deleted...
sub DeleteHydrogens {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }

  # Is it an element symbol?
  if (!$This->{AtomicNumber}) {
    return 0;
  }

  my($Molecule, $Neighbor, $HydrogensDeleted, @Neighbors);

  $Molecule = $This->GetProperty('Molecule');
  $HydrogensDeleted = 0;
  @Neighbors = $This->GetNeighbors();

  NEIGHBOR: for $Neighbor (@Neighbors) {
    if (!$Neighbor->IsHydrogen()) {
      next NEIGHBOR;
    }
    $Molecule->_DeleteAtom($Neighbor);
    $HydrogensDeleted++;
  }

  return $HydrogensDeleted;
}

# Copy atom and all its associated data...
sub Copy {
  my($This) = @_;
  my($Atom);

  $Atom = Storable::dclone($This);

  return $Atom;
}

# Get atomic invariant value...
#
sub GetAtomicInvariantValue {
  my($This, $AtomicInvariant) = @_;
  my($Value);

  $Value = "";

  ATOMICVARIANT: {
    if ($AtomicInvariant =~ /^(AS|AtomSymbol|ElementSymbol)$/i) {
      $Value = $This->GetAtomSymbol();
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(X|NumOfNonHydrogenAtomNeighbors|NumOfHeavyAtomNeighbors)$/i) {
      $Value = $This->GetNumOfNonHydrogenAtomNeighbors();
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(BO|SumOfBondOrdersToNonHydrogenAtoms|SumOfBondOrdersToHeavyAtoms)$/i) {
      $Value = $This->GetSumOfBondOrdersToNonHydrogenAtoms();
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(LBO|LargestBondOrderToNonHydrogenAtoms|LargestBondOrderToHeavyAtoms)$/i) {
      $Value = $This->GetLargestBondOrderToNonHydrogenAtoms();
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(H|NumOfImplicitAndExplicitHydrogens)$/i) {
      $Value = $This->GetNumOfHydrogens();
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(SB|NumOfSingleBondsToNonHydrogenAtoms|NumOfSingleBondsToHeavyAtoms)$/i) {
      $Value = $This->GetNumOfSingleBondsToNonHydrogenAtoms();
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(DB|NumOfDoubleBondsToNonHydrogenAtoms|NumOfDoubleBondsToHeavyAtoms)$/i) {
      $Value = $This->GetNumOfDoubleBondsToNonHydrogenAtoms();
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(TB|NumOfTripleBondsToNonHydrogenAtoms|NumOfTripleBondsToHeavyAtoms)$/i) {
      $Value = $This->GetNumOfTripleBondsToNonHydrogenAtoms();
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(AB|NumOfAromaticBondsToNonHydrogenAtoms|NumOfAromaticBondsToHeavyAtoms)$/i) {
      $Value = $This->GetNumOfAromaticBondsToNonHydrogenAtoms();
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(FC|FormalCharge)$/i) {
      $Value = $This->GetFormalCharge();
      $Value = defined $Value ? $Value : 0;
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(T|TotalNumOfAtomNeighbors)$/i) {
      $Value = $This->GetNumOfNonHydrogenAtomNeighbors() + $This->GetNumOfHydrogens();
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(TSB|TotalNumOfSingleBonds)$/i) {
      $Value = $This->GetNumOfSingleBondsToNonHydrogenAtoms() + $This->GetNumOfHydrogens();
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(Ar|Aromatic)$/i) {
      $Value = $This->IsAromatic() ? 1 : 0;
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(RA|RingAtom)$/i) {
      $Value = $This->IsInRing() ? 1 : 0;
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(Str|Stereochemistry)$/i) {
      $Value = $This->GetStereochemistry();
      $Value= (defined($Value) && ($Value =~ /^(R|S)$/i)) ? $Value : '';
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(AN|AtomicNumber)$/i) {
      $Value = $This->GetAtomicNumber();
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(AM|AtomicMass)$/i) {
      $Value = round($This->GetExactMass(), 4) + 0;
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(MN|MassNumber)$/i) {
      $Value = $This->GetMassNumber();
      last ATOMICVARIANT;
    }
    if ($AtomicInvariant =~ /^(SM|SpinMultiplicity)$/i) {
      $Value = $This->GetSpinMultiplicity();
      $Value = defined $Value ? $Value : '';
      last ATOMICVARIANT;
    }
    $Value = "";
    carp "Warning: ${ClassName}->GetAtomicInvariantValue: Unknown atomic invariant $AtomicInvariant...";
  }

  return $Value;
}

# Get period number of the atom..
#
sub GetPeriodNumber {
  my($This) = @_;

  # Is PeriodNumber property explicitly set?
  if ($This->HasProperty('PeriodNumber')) {
    return $This->GetProperty('PeriodNumber');
  }
  my($AtomicNumber, $PeriodNumber);

  $AtomicNumber = $This->{AtomicNumber};
  if (!$AtomicNumber) {
    return 0;
  }

  $PeriodNumber = PeriodicTable::GetElementPeriodNumber($AtomicNumber);

  return $PeriodNumber;
}

# Get group number of the atom..
#
sub GetGroupNumber {
  my($This) = @_;

  # Is GroupNumber property explicitly set?
  if ($This->HasProperty('GroupNumber')) {
    return $This->GetProperty('GroupNumber');
  }
  my($AtomicNumber, $GroupNumber);

  $AtomicNumber = $This->{AtomicNumber};
  if (!$AtomicNumber) {
    return 0;
  }

  $GroupNumber = PeriodicTable::GetElementGroupNumber($AtomicNumber);

  return $GroupNumber;
}

# Is it a specified topological pharmacophore atom type?
#
sub IsTopologicalPharmacophoreType {
  my($This, $Type) = @_;

  return $This->_IsFunctionalClassType($Type);
}

# Is it a specified functional class atom type?
#
sub IsFunctionalClassType {
  my($This, $Type) = @_;

  return $This->_IsFunctionalClassType($Type);
}

# Is it a specified functional/topological pharmacophore atom type?
#
sub _IsFunctionalClassType {
  my($This, $Type) = @_;
  my($Value);

  $Value = 0;

  TYPE: {
    if ($Type =~ /^(HBD|HydrogenBondDonor)$/i) {
      $Value = $This->IsHydrogenBondDonor();
      last TYPE;
    }
    if ($Type =~ /^(HBA|HydrogenBondAcceptor)$/i) {
      $Value = $This->IsHydrogenBondAcceptor();
      last TYPE;
    }
    if ($Type =~ /^(PI|PositivelyIonizable)$/i) {
      $Value = $This->IsPositivelyIonizable();
      last TYPE;
    }
    if ($Type =~ /^(NI|NegativelyIonizable)$/i) {
      $Value = $This->IsNegativelyIonizable();
      last TYPE;
    }
    if ($Type =~ /^(H|Hydrophobic)$/i) {
      $Value = $This->IsHydrophobic();
      last TYPE;
    }
    if ($Type =~ /^(Ar|Aromatic)$/i) {
      $Value = $This->IsAromatic();
      last TYPE;
    }
    if ($Type =~ /^(Hal|Halogen)$/i) {
      $Value = $This->IsHalogen();
      last TYPE;
    }
    if ($Type =~ /^(RA|RingAtom)$/i) {
      $Value = $This->IsInRing();
      last TYPE;
    }
    if ($Type =~ /^(CA|ChainAtom)$/i) {
      $Value = $This->IsNotInRing();
      last TYPE;
    }
    $Value = 0;
    carp "Warning: ${ClassName}->_IsType: Unknown functional/pharmacohore type $Type...";
  }
  return $Value;
}

# Is it a Hydrogen atom?
sub IsHydrogen {
  my($This) = @_;

  return ($This->{AtomicNumber} == 1) ? 1 : 0;
}

# Is it a Carbon atom?
sub IsCarbon {
  my($This) = @_;

  return ($This->{AtomicNumber} == 6) ? 1 : 0;
}

# Is it a Nitrogen atom?
sub IsNitrogen {
  my($This) = @_;

  return ($This->{AtomicNumber} == 7) ? 1 : 0;
}

# Is it a Oxygen atom?
sub IsOxygen {
  my($This) = @_;

  return ($This->{AtomicNumber} == 8) ? 1 : 0;
}

# Is it a Fluorine atom?
sub IsFluorine {
  my($This) = @_;

  return ($This->{AtomicNumber} == 9) ? 1 : 0;
}

# Is it a Silicon atom?
sub IsSilicon {
  my($This) = @_;

  return ($This->{AtomicNumber} == 14) ? 1 : 0;
}

# Is it a Phosphorus atom?
sub IsPhosphorus {
  my($This) = @_;

  return ($This->{AtomicNumber} == 15) ? 1 : 0;
}

# Is it a Sulphur atom?
sub IsSulphur {
  my($This) = @_;

  return $This->IsSulfur();
}

# Is it a Sulfur atom?
sub IsSulfur {
  my($This) = @_;

  return ($This->{AtomicNumber} == 16) ? 1 : 0;
}

# Is it a Chlorine atom?
sub IsChlorine {
  my($This) = @_;

  return ($This->{AtomicNumber} == 17) ? 1 : 0;
}

# Is it a Arsenic atom?
sub IsArsenic {
  my($This) = @_;

  return ($This->{AtomicNumber} == 33) ? 1 : 0;
}

# Is it a Selenium atom?
sub IsSelenium {
  my($This) = @_;

  return ($This->{AtomicNumber} == 34) ? 1 : 0;
}

# Is it a Bromine atom?
sub IsBromine {
  my($This) = @_;

  return ($This->{AtomicNumber} == 35) ? 1 : 0;
}

# Is it a Tellurium atom?
sub IsTellurium {
  my($This) = @_;

  return ($This->{AtomicNumber} == 52) ? 1 : 0;
}

# Is it a Iodine atom?
sub IsIodine {
  my($This) = @_;

  return ($This->{AtomicNumber} == 53) ? 1 : 0;
}

# Is it a hetro atom? (N, O, F, P, S, Cl, Br, I)
sub IsHeteroAtom {
  my($This) = @_;

  return ($This->{AtomicNumber} =~ /^(7|8|9|15|16|17|35|53)$/) ? 1 : 0;
}

# Is it a halogen atom? (F, Cl, Br, I)
sub IsHalogen {
  my($This) = @_;

  return ($This->{AtomicNumber} =~ /^(9|17|35|53)$/) ? 1 : 0;
}

# Is it classified as metallic?
sub IsMetallic {
  my($This) = @_;
  my($Classification);

  $Classification = PeriodicTable::GetElementClassification($This->{AtomicNumber});

  return ($Classification =~ /^Metallic$/i) ? 1 : 0;
}

# Is it a non carbon or hydrogen atom? (C, H)
sub IsNonCarbonOrHydrogen {
  my($This) = @_;

  return ($This->{AtomicNumber} =~ /^(1|6)$/) ? 0 : 1;
}

# Is it a polar atom? ( N, O,  P, S)
sub IsPolarAtom {
  my($This) = @_;

  return ($This->{AtomicNumber} =~ /^(7|8|15|16)$/) ? 1 : 0;
}

# Is it an isotope?
sub IsIsotope {
  my($This) = @_;

  my($AtomicNumber) = $This->{AtomicNumber};
  if (!$AtomicNumber) {
    return 0;
  }

  if (!$This->HasProperty('MassNumber')) {
    return 0;
  }
  my($MassNumber, $MostAbundantMassNumber);

  $MassNumber = $This->GetProperty('MassNumber');
  $MostAbundantMassNumber = PeriodicTable::GetElementMostAbundantNaturalIsotopeMassNumber($AtomicNumber);

  return ($MassNumber == $MostAbundantMassNumber) ? 0 : 1;
}

# Is it a terminal atom?
sub IsTerminal {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }

  return ($This->GetNumOfNonHydrogenAtomNeighbors() <= 1) ? 1 : 0

}

# Is aromatic property set for the atom?
sub IsAromatic {
  my($This) = @_;
  my($Aromatic);

  $Aromatic = $This->GetAromatic();

  return (defined($Aromatic) && $Aromatic) ? 1 : 0;
}

# Is this a hydrogen atom and attached to one of these atoms: N, O, P, S
sub IsPolarHydrogen {
  my($This) = @_;

  if (!$This->IsHydrogen()) {
    return 0;
  }

  my(@Bonds);
  @Bonds = $This->GetBonds();
  if (@Bonds > 1) {
    return 0;
  }

  my($Bond, $BondedAtom);
  ($Bond) = @Bonds;
  $BondedAtom = $Bond->GetBondedAtom($This);

  return $BondedAtom->IsPolarAtom() ? 1 : 0;
}

# Is it a hydrogen bond donor atom?
#
sub IsHBondDonor {
  my($This, $HydrogenBondsType) = @_;

  return $This->IsHydrogenBondDonor($HydrogenBondsType);
}

# The currrent release of MayaChemTools supports identification of two types of
# hydrogen bond donor and acceptor atoms with these names:
#
# HBondsType1 or HydrogenBondsType1
# HBondsType2 or HydrogenBondsType2
#
# The names of these hydrogen bond types are rather arbitrary. However, their
# definitions have specific meaning and are as follows:
#
# HydrogenBondsType1 [ Ref 60-61, Ref 65-66 ]:
#   . Donor: NH, NH2, NH3, OH - Any N and O with available H
#   . Acceptor: N[!H], O - Any N without available H and any O
#
# HydrogenBondsType2 [ Ref 91 ]:
#   . Donor: NH, NH2, NH3, OH - Any N and O with availabe H
#   . Acceptor: N, O - Any N and O
#
# Note:
#   . HydrogenBondsType2 definition corresponds to Rule of 5.
#

# Is it a hydrogen bond donor atom?
#
# The currrent release of MayaChemTools supports identification of two types of
sub IsHydrogenBondDonor {
  my($This, $HydrogenBondsType) = @_;
  my($Status);

  $HydrogenBondsType = defined $HydrogenBondsType ? $HydrogenBondsType : 'HBondsType1';
  $Status = 0;

  HYDROGENBONDSTYPE: {

      if ($HydrogenBondsType =~ /^(HBondsType1|HydrogenBondsType1)$/i) {
	$Status = $This->_IsHydrogenBondDonorOfType1();
	last HYDROGENBONDSTYPE;
      }

      if ($HydrogenBondsType =~ /^(HBondsType2|HydrogenBondsType2)$/i) {
	$Status = $This->_IsHydrogenBondDonorOfType2();
	last HYDROGENBONDSTYPE;
      }

      $Status = 0;
      carp "Warning: ${ClassName}->IsHydrogenBondDonor: The current release of MayaChemTools doesn't support specified value, $HydrogenBondsType, for HydrogenBondsType. Valid values: HBondsType1, HydrogenBondsType1, HBondsType2 HydrogenBondsType2 ...";
  }

  return $Status;
}

# Is it a MayaChemTools HBondType1 hydrogen bond donor atom?
#
sub _IsHydrogenBondDonorOfType1 {
  my($This) = @_;

  return $This->_IsHydrogenBondDonorOfType1OrType2();
}

# Is it a MayaChemTools HBondType2 hydrogen bond donor atom?
#
sub _IsHydrogenBondDonorOfType2 {
  my($This) = @_;

  return $This->_IsHydrogenBondDonorOfType1OrType2();
}

# Is it a hydrogen bond donor atom of MayaChemTools Type1 or Type2?
#
# HydrogenBondDonor definition [ Ref 60-61, Ref 65-66, Ref 91 ]: NH, NH2, OH
#
# In other words:
#   . NH, NH2 - Nitrogen atom with available hydrogen
#   . OH - Oxygen atom with avilable hydrogen
#
sub _IsHydrogenBondDonorOfType1OrType2 {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return 0;
  }

  # Is it N or O?
  if ($This->{AtomicNumber} !~ /^(7|8)$/) {
    return 0;
  }

  # Any explicitly attached hydrogens?
  if ($This->GetExplicitHydrogens()) {
    return 1;
  }

  # Any missing hydrogens?
  return $This->GetNumOfMissingHydrogens() ? 1 : 0;
}

# Is it a hydrogen bond acceptor atom?
#
sub IsHBondAcceptor {
  my($This, $HydrogenBondsType) = @_;

  return $This->IsHydrogenBondAcceptor($HydrogenBondsType);
}

# Is it a hydrogen bond acceptor atom?
#
sub IsHydrogenBondAcceptor {
  my($This, $HydrogenBondsType) = @_;
  my($Status);

  $HydrogenBondsType = defined $HydrogenBondsType ? $HydrogenBondsType : 'HBondsType1';
  $Status = 0;

  HYDROGENBONDSTYPE: {

      if ($HydrogenBondsType =~ /^(HBondsType1|HydrogenBondsType1)$/i) {
	$Status = $This->_IsHydrogenBondAcceptorOfType1();
	last HYDROGENBONDSTYPE;
      }

      if ($HydrogenBondsType =~ /^(HBondsType2|HydrogenBondsType2)$/i) {
	$Status = $This->_IsHydrogenBondAcceptorOfType2();
	last HYDROGENBONDSTYPE;
      }

      $Status = 0;
      carp "Warning: ${ClassName}->IsHydrogenBondAcceptor: The current release of MayaChemTools doesn't support specified value, $HydrogenBondsType, for HydrogenBondsType. Valid values: HBondsType1, HydrogenBondsType1, HBondsType2 HydrogenBondsType2 ...";
  }

  return $Status;
}

# Is it a MayaChemTools HBondType1 hydrogen bond acceptor atom?
#
# HydrogenBondAcceptor definition [ Ref 60-61, Ref 65-66 ]: N[!H], O
#
# In other words:
#   . N[!H] - Nitrogen atom with no hydrogen
#   . O - Oxygen atom
#
sub _IsHydrogenBondAcceptorOfType1 {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return 0;
  }

  # Is it N or O?
  if ($This->{AtomicNumber} !~ /^(7|8)$/) {
    return 0;
  }

  # Is it O?
  if ($This->{AtomicNumber} == 8 ) {
    return 1;
  }

  # Any explicitly attached hydrogens?
  if ($This->GetExplicitHydrogens()) {
    return 0;
  }

  # Any missing hydrogens?
  return $This->GetNumOfMissingHydrogens() ? 0 : 1;
}

# Is it a MayaChemTools HBondType2 hydrogen bond acceptor atom?
#
# HydrogenBondAcceptor definition [ Ref 91 ]: N, O
#
# In other words:
#   . Any Nitrogen or Oxygen atom
#
# Note:
#   . HydrogenBondsType2 definition corresponds to Rule of 5.
#
sub _IsHydrogenBondAcceptorOfType2 {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return 0;
  }

  return ($This->{AtomicNumber} =~ /^(7|8)$/) ? 1 : 0;
}

# Is it a positively ionizable atom?
#
# PositivelyIonizable defintion [ Ref 60-61, Ref 65-66 ]: +, NH2
#
# In other words:
#   . Any atom with positve formal charge
#   . NH2 - Nitogen atom in amino group
#
sub IsPositivelyIonizable {
  my($This) = @_;
  my($FormalCharge);

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return 0;
  }

  # Any explicit positive formal charge?
  $FormalCharge = $This->GetFormalCharge();
  if (defined($FormalCharge) && $FormalCharge > 0) {
    return 1;
  }

  # Is it  N?
  if ($This->{AtomicNumber} != 7 ) {
    return 0;
  }

  return ($This->GetNumOfHydrogens() == 2) ? 1 : 0;
}

# Is it a negatively ionizable atom?
#
# NegativelyIonizable definition [ Ref 60-61, Ref 65-66 ]: -, C(=O)OH, S(=O)OH, P(=O)OH
#
# In other words:
#   . Any atom with negative formal charge
#   . Carbon atom in C(=O)OH group
#   . Phosphorous in P(=O)OH group
#   . Sulfur atom in S(=O)OH group
#
sub IsNegativelyIonizable {
  my($This) = @_;
  my($FormalCharge);

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return 0;
  }

  # Any explicit negative formal charge?
  $FormalCharge = $This->GetFormalCharge();
  if (defined($FormalCharge) && $FormalCharge < 0) {
    return 1;
  }

  # Is it C, P or S?
  if ($This->{AtomicNumber} !~ /^(6|15|16)$/ ) {
    return 0;
  }

  # Collect oxygens connected to C, P or S with single or double bonds and not connected to
  # any other heavy atom...
  my($Neighbor, $NeighborOxygenBondOrder, $NumOfNeighborOxygensWithSingleBonds, $NumOfNeighborOxygensWithDoubleBonds);

  $NumOfNeighborOxygensWithSingleBonds = 0; $NumOfNeighborOxygensWithDoubleBonds = 0;

  NEIGHBOR: for $Neighbor ($This->GetNeighbors()) {
    # Is it an oxygen?
    if ($Neighbor->{AtomicNumber} != 8) {
      next NEIGHBOR;
    }
    # Is oxygent connected to only heavy atom?
    if ($Neighbor->GetNumOfHeavyAtomNeighbors() != 1) {
      next NEIGHBOR;
    }
    $NeighborOxygenBondOrder = $This->GetBondToAtom($Neighbor)->GetBondOrder();

    if ($NeighborOxygenBondOrder == 2) {
      $NumOfNeighborOxygensWithDoubleBonds++;
    }
    elsif ($NeighborOxygenBondOrder == 1) {
      $NumOfNeighborOxygensWithSingleBonds++;
    }
  }
  return ($NumOfNeighborOxygensWithDoubleBonds >= 1 && $NumOfNeighborOxygensWithSingleBonds >= 1) ? 1 : 0;
}

# Is it a liphophilic atom?
#
# Lipophilic definition [ Ref 60-61, Ref 65-66 ]: C(C)(C)(C)(C), Cl, Br, I, S(C)(C)
#
# In other words:
#   . C(C)(C)(C)(C) - Carbon atom connected to only other carbons
#   . Chlorine, Bromine or Iodine atom
#   . S(C)(C) - Sulfur connected to two carbons
#
sub IsLipophilic {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return 0;
  }

  # Is it Cl, Br, I?
  if ($This->{AtomicNumber} =~ /^(17|35|53)$/) {
    return 1;
  }

  # Is it C, S?
  if ($This->{AtomicNumber} !~ /^(6|16)$/) {
    return 0;
  }

  # Are all heavy atom neighbors Carbons?
  my($HeavyAtomNeighbor, @HeavyAtomNeighbors);
  @HeavyAtomNeighbors = ();
  @HeavyAtomNeighbors = $This->GetHeavyAtomNeighbors();

  for $HeavyAtomNeighbor (@HeavyAtomNeighbors) {
    if ($HeavyAtomNeighbor->{AtomicNumber} != 6) {
      return 0;
    }
  }

  # Does sulfur has two carbon neighbors?
  if ($This->{AtomicNumber} == 16) {
    if (@HeavyAtomNeighbors != 2) {
      return 0;
    }
  }
  return 1;
}

# Is it hydrophobic?
#
sub IsHydrophobic {
  my($This) = @_;

  return $This->IsLipophilic();
}

# Is it a Nitrogen atom in Guadinium group?
#
sub IsGuadiniumNitrogen {
  my($This) = @_;

  # Is it Nitrogen?
  if (!$This->IsNitrogen()) {
    return 0;
  }

  # Is it connected to a Guadinium Carbon?
  my($AtomNeighbor);

  for $AtomNeighbor ($This->GetNonHydrogenAtomNeighbors()) {
    if ($AtomNeighbor->IsGuadiniumCarbon()) {
      return 1;
    }
  }

  return 0;
}

# Is it a Carbon atom in Guadinium group?
#
# Guadinium group definition:
#
#   R2N-C(=NR)-(NR2) or R2N-C(=NR2+)-(NR2)
#
#   where:
#      . R = Hydrogens or group of atoms attached through Carbon
#      . Only one of the three Nitrogens has a double bond to Carbon and has optional
#        formal charge allowing it to be neutral or charged state
#
sub IsGuadiniumCarbon {
  my($This) = @_;

  # Is it Carbon?
  if (!$This->IsCarbon()) {
    return 0;
  }

  # Match atom neighborhood...
  my($CentralAtomSpec, @NbrAtomSpecsRef, @NbrBondSpecsRef, @NbrOfNbrAtomSpecsRef);

  $CentralAtomSpec = 'C.X3.BO4';
  @NbrAtomSpecsRef = ('N.FC0', 'N.FC0', 'N.FC0,N.FC+1');
  @NbrBondSpecsRef = ('-', '-', '=');
  @NbrOfNbrAtomSpecsRef = ('C,H', 'C,H', 'C,H');

  if ($This->DoesAtomNeighborhoodMatch($CentralAtomSpec, \@NbrAtomSpecsRef, \@NbrBondSpecsRef, \@NbrOfNbrAtomSpecsRef)) {
    return 1;
  }

  return 0;
}

# Is it a Nitrogen atom in Amide group?
#
sub IsAmideNitrogen {
  my($This) = @_;

  # Is it Nitrogen?
  if (!$This->IsNitrogen()) {
    return 0;
  }

  # Is it connected to a Amide Carbon?
  my($AtomNeighbor);

  for $AtomNeighbor ($This->GetNonHydrogenAtomNeighbors()) {
    if ($AtomNeighbor->IsAmideCarbon()) {
      return 1;
    }
  }

  return 0;
}

# Is it a Carbon atom in Amide group?
#
# Amide group definition: R-C(=O)-N(-R')-R''
#
#   where:
#      . R = Hydrogen or groups of atoms attached through Carbon
#      . R' = Hydrogens or groups of atoms attached through Carbon or hetro atoms
#      . R'' = Hydrogens or groups of atoms attached through Carbon or hetro atoms
#
sub IsAmideCarbon {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return 0;
  }

  # Is it Carbon?
  if (!$This->IsCarbon()) {
    return 0;
  }

  # Match atom neighborhood...
  my($CentralAtomSpec, @NbrAtomSpecsRef, @NbrBondSpecsRef, @NbrOfNbrAtomSpecsRef);

  $CentralAtomSpec = 'C.X3.BO4,C.X2.BO3';
  @NbrAtomSpecsRef = ('C,H', 'O', 'N');
  @NbrBondSpecsRef = ('-', '=', '-');
  @NbrOfNbrAtomSpecsRef = ('C,H', 'C', 'C,H,N,O,S,P,F,Cl,Br,I');

  if ($This->DoesAtomNeighborhoodMatch($CentralAtomSpec, \@NbrAtomSpecsRef, \@NbrBondSpecsRef, \@NbrOfNbrAtomSpecsRef)) {
    return 1;
  }

  return 0;
}

# Is it a Oxygen atom in Carboxylate group?
#
sub IsCarboxylateOxygen {
  my($This) = @_;

  return $This->_MatchCarboxylateAndOrCarboxylOxygen('Carboxylate');
}

# Is it a Carbon atom in Carboxylate group?
#
# Carboxyl group definition: R-C(=O)-O-
#
sub IsCarboxylateCarbon {
  my($This) = @_;

  return $This->_MatchCarboxylateAndOrCarboxylCarbon('Carboxylate');
}

# Is it a Oxygen atom in Carboxyl group?
#
sub IsCarboxylOxygen {
  my($This) = @_;

  return $This->_MatchCarboxylateAndOrCarboxylOxygen('Carboxyl');
}

# Is it a Carbon atom in Carboxyl group?
#
# Carboxyl group definition: R-C(=O)-OH
#
sub IsCarboxylCarbon {
  my($This) = @_;

  return $This->_MatchCarboxylateAndOrCarboxylCarbon('Carboxyl');
}

# Match Carboxylate and/or Carboxyl oxygen...
#
sub _MatchCarboxylateAndOrCarboxylOxygen {
  my($This, $Mode) = @_;

  # Is it Oxygen?
  if (!$This->IsOxygen()) {
    return 0;
  }

  # Is it connected to a Carboxylate Carbon?
  my($AtomNeighbor);

  for $AtomNeighbor ($This->GetNonHydrogenAtomNeighbors()) {
    if ($AtomNeighbor->_MatchCarboxylateAndOrCarboxylCarbon($Mode)) {
      return 1;
    }
  }

  return 0;
}

# Match Carboxylate and Carboxyl Carbon
#
# Carboxylate group definition: R-C(=O)-O-
# Carboxyl group definition: R-C(=O)-OH
#
#   where:
#      . R = Hydrogens or groups of atoms attached through Carbon
#
sub _MatchCarboxylateAndOrCarboxylCarbon {
  my($This, $Mode) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return 0;
  }

  # Is it Carbon?
  if (!$This->IsCarbon()) {
    return 0;
  }

  # Match atom neighborhood...
  my($CentralAtomSpec, @NbrAtomSpecsRef, @NbrBondSpecsRef, @NbrOfNbrAtomSpecsRef);

  $CentralAtomSpec = 'C.X3.BO4,C.X2.BO3';
  MODE: {
    if ($Mode =~ /^Carboxylate$/i) {
      @NbrAtomSpecsRef = ('C,H', 'O', 'O.X1.FC-1');
      last MODE;
    }
    if ($Mode =~ /^Carboxyl$/i) {
      @NbrAtomSpecsRef = ('C,H', 'O', 'O.X1.FC0');
      last MODE;
    }
    if ($Mode =~ /^CarboxylateOrCarboxyl$/i) {
      @NbrAtomSpecsRef = ('C,H', 'O', 'O.X1.FC-1,O.X1.FC0');
      last MODE;
    }
    carp "Warning: ${ClassName}->_MatchCarboxylateAndCarboxylCarbon.: Unknown mode $Mode...";
    return 0;
  }
  @NbrBondSpecsRef = ('-', '=', '-');
  @NbrOfNbrAtomSpecsRef = ('C,H', 'C', 'C');

  if ($This->DoesAtomNeighborhoodMatch($CentralAtomSpec, \@NbrAtomSpecsRef, \@NbrBondSpecsRef, \@NbrOfNbrAtomSpecsRef)) {
    return 1;
  }

  return 0;
}

# Is it a Oxygen atom in Phosphate group?
#
sub IsPhosphateOxygen {
  my($This) = @_;

  # Is it Oxygen?
  if (!$This->IsOxygen()) {
    return 0;
  }

  # Is it connected to a Phosphate Phosphorus?
  my($AtomNeighbor);

  for $AtomNeighbor ($This->GetNonHydrogenAtomNeighbors()) {
    if ($AtomNeighbor->IsPhosphatePhosphorus()) {
      return 1;
    }
  }

  return 0;
}

# Is it a Phosphorus atom in Phosphate group?
#
# Phosphate group definition: AO-(O=)P(-OA)-OA
#
#   where:
#      . A = Any Groups of atoms including hydrogens
#
sub IsPhosphatePhosphorus {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return 0;
  }

  # Is it Phosphorus?
  if (!$This->IsPhosphorus()) {
    return 0;
  }

  # Match atom neighborhood...
  my($CentralAtomSpec, @NbrAtomSpecsRef, @NbrBondSpecsRef, @NbrOfNbrAtomSpecsRef);

  $CentralAtomSpec = 'P.X4.BO5';
  @NbrAtomSpecsRef = ('O', 'O', 'O', 'O');
  @NbrBondSpecsRef = ('-', '=', '-', '-');
  @NbrOfNbrAtomSpecsRef = (undef, undef, undef, undef);

  if ($This->DoesAtomNeighborhoodMatch($CentralAtomSpec, \@NbrAtomSpecsRef, \@NbrBondSpecsRef, \@NbrOfNbrAtomSpecsRef)) {
    return 1;
  }

  return 0;
}


# Match central atom and its neighborhood using specified atom and bonds specifications...
#
# Let:
#   AS = Atom symbol corresponding to element symbol, atomic number (#n) or any
#        atom (A)
#
#   X<n>   = Number of non-hydrogen atom neighbors or heavy atoms attached to atom
#   T<n>   = Total number of atom neighbors including implcit and explicit hydrogens
#   BO<n> = Sum of bond orders to non-hydrogen atom neighbors or heavy atoms attached to atom
#   LBO<n> = Largest bond order of non-hydrogen atom neighbors or heavy atoms attached to atom
#   SB<n> = Number of single bonds to non-hydrogen atom neighbors or heavy atoms attached to atom
#   TSB<n> = Total number of single bonds to atom neighbors including implcit and explicit hydrogens
#   DB<n> = Number of double bonds to non-hydrogen atom neighbors or heavy atoms attached to atom
#   TB<n> = Number of triple bonds to non-hydrogen atom neighbors or heavy atoms attached to atom
#   H<n>   = Number of implicit and explicit hydrogens for atom
#   Ar     = Aromatic annotation indicating whether atom is aromatic
#   RA or RA<n>  = Ring atom annotation indicating whether atom is a ring
#   TR<n>  = Total number of rings containing atom
#   FC<+n/-n> = Formal charge assigned to atom
#   MN<n> = Mass number indicating isotope other than most abundant isotope
#   SM<n> = Spin multiplicity of atom. Possible values: 1 (singlet), 2 (doublet) or 3 (triplet)
#
# Then:
#
#   Atom specification corresponds to:
#
#     AS.X<n>.T<n>.BO<n>.LBO<n>.<SB><n>.TSB<n>.<DB><n>.<TB><n>.H<n>.Ar.RA<n>.TR<n>FC<+n/-n>.MN<n>.SM<n>
#
# Except for AS which is a required atomic invariant in atom specification, all other atomic invariants are
# optional. For an atom specification to match an atom, the values of all specified atomic invariants must
# match. Exclamation in from of atomic invariant can be used to negate its effect during the match.
#
# A comma delimited atom specification string is used to match any one of the specifed atom specification.
#
# Notes:
#   . During atom specification match to an atom, the first atomic invariant is always assumed to
#     atom symbol.
#   . Atom match specfication is based on AtomicInvariantAtomTypes implemented in
#     AotmTypes::AtomicInvariantAtomType.pm module
#
# Examples:
#     . ('N', 'N', 'N')
#     . ('N.FC0', 'N.FC0', 'N,N.FC+1.H1')
#     . ('N.H2', 'N.H2', 'N.H1')
#     . ('C,N', '!N', '!H')
#     . ('C,N', 'N.Ar', 'N.R5')
#
# Let:
#   -|1|s|Single = Single bond
#   =|2|d|Double = Double bond
#   #|3|t|Triple  = Triple bond
#   :|1.5|a|Ar|Aromatic = Aromatic bond
#
#   @|RB|Ring = Ring bond
#   ~|*|Any = Any bond
#
# Then:
#
#   Bond specification corresponds to:
#
#     -.:
#     =.@
#     Double.Aromatic
#
# For a bond specification to match bond between two atoms, the values of all specified bond symbols must
# match. Exclamation in from of bond symbol can be used to negate its effect during the match.
#
# A comma delimited bond specification string is used to match any one of the specifed atom specification.
#
sub DoesAtomNeighborhoodMatch {
  my($CentralAtom, $CentralAtomSpec, $NbrAtomSpecsRef, $NbrBondSpecsRef, $NbrOfNbrAtomSpecsRef) = @_;
  my($NumOfNbrAtomSpecs, $NumOfNbrBondSpecs, $NumOfNbrOfNbrAtomSpecs);

  # Is this atom in a molecule?
  if (!$CentralAtom->HasProperty('Molecule')) {
    return 0;
  }

  $NumOfNbrAtomSpecs = defined $NbrAtomSpecsRef ? scalar @{$NbrAtomSpecsRef} : 0;
  $NumOfNbrBondSpecs = defined $NbrBondSpecsRef ? scalar @{$NbrBondSpecsRef} : 0;
  $NumOfNbrOfNbrAtomSpecs = defined $NbrOfNbrAtomSpecsRef ? scalar @{$NbrOfNbrAtomSpecsRef} : 0;

  # Validate number of specifications...
  if ($NumOfNbrBondSpecs && ($NumOfNbrAtomSpecs != $NumOfNbrBondSpecs)) {
    carp "Warning: ${ClassName}->DoesAtomNeighborhoodMatch: Number of specified central atom, $NumOfNbrAtomSpecs, and bond, $NumOfNbrBondSpecs, specifications must be same; No neighborhood match performed ...";
    return 0;
  }

  if ($NumOfNbrOfNbrAtomSpecs && ($NumOfNbrOfNbrAtomSpecs != $NumOfNbrAtomSpecs)) {
    carp "Warning: ${ClassName}->DoesAtomNeighborhoodMatch: Number of specified central atom, $NumOfNbrAtomSpecs, and neighbor of neighbor atoms specifications, $NumOfNbrOfNbrAtomSpecs, must be same; No neighborhood match performed ...";
    return 0;
  }

  # Sort atom and bond specifications in terms of being most specific to least specific..
  ($NbrAtomSpecsRef, $NbrBondSpecsRef, $NbrOfNbrAtomSpecsRef) = $CentralAtom->_SortSpecificationsForAtomNeighborhoodMatch($NbrAtomSpecsRef, $NbrBondSpecsRef, $NbrOfNbrAtomSpecsRef);

  # Does central atom specification match?
  if (!$CentralAtom->_DoesAtomSpecificationMatch($CentralAtomSpec)) {
    return 0;
  }

  # No neighbors to match...
  if (!$NumOfNbrAtomSpecs) {
    return 1;
  }

  # Match neighbors...
  my($NbrSpecsMatched, $NbrSpecCount, $NbrSpecMatchCount, %NbrSpecAlreadyMatchedMap);

  $NbrSpecCount = $NumOfNbrAtomSpecs;
  $NbrSpecMatchCount = 0;

  %NbrSpecAlreadyMatchedMap = ();
  ($NbrSpecsMatched, $NbrSpecMatchCount) = $CentralAtom->_MatchAtomNeighborhoodUsingAtomBondSpecs($NbrSpecCount, $NbrSpecMatchCount, \%NbrSpecAlreadyMatchedMap, $NbrAtomSpecsRef, $NbrBondSpecsRef, $NbrOfNbrAtomSpecsRef);

  if ($NbrSpecsMatched) {
    # It's match...
    return 1;
  }

  # Match central atom's missing hydrogens with any unmatched atom
  # and bond specifications...
  #
  ($NbrSpecsMatched, $NbrSpecMatchCount) = $CentralAtom->_MatchAtomNeighborhoodUsingMissingHydrogens($NbrSpecCount, $NbrSpecMatchCount, \%NbrSpecAlreadyMatchedMap, $NbrAtomSpecsRef, $NbrBondSpecsRef, $NbrOfNbrAtomSpecsRef);

  if ($NbrSpecsMatched) {
    # It's match...
    return 1;
  }

  # No match...
  return 0;
}

# Match central atom neighborhood atom and bond specifications...
#
sub _MatchAtomNeighborhoodUsingAtomBondSpecs {
  my($CentralAtom, $NbrSpecCount, $NbrSpecMatchCount, $NbrSpecAlreadyMatchedRef, $NbrAtomSpecsRef, $NbrBondSpecsRef, $NbrOfNbrAtomSpecsRef) = @_;
  my($Index, $NbrAtom, $NbrAtomSpec, $NbrBondSpec, $NbrOfNbrAtom, $NbrOfNbrAtomSpec, $MatchNbrOfNbrAtomSpecs, $NbrSpecsMatched);

  $MatchNbrOfNbrAtomSpecs = (defined $NbrOfNbrAtomSpecsRef && scalar @{$NbrOfNbrAtomSpecsRef}) ? 1 : 0;

  $NbrSpecsMatched = 0;

  # Match central atom's  immediate neighbors atom and bond specifications...
  NBRATOM:  for $NbrAtom ($CentralAtom->GetNeighbors()) {
    NBRATOMSPEC: for $Index (0 .. ($NbrSpecCount - 1)) {
      if (exists $NbrSpecAlreadyMatchedRef->{$Index}) {
	next NBRATOMSPEC;
      }
      $NbrAtomSpec = $NbrAtomSpecsRef->[$Index];
      $NbrBondSpec = $NbrBondSpecsRef->[$Index];

      $NbrOfNbrAtomSpec = $MatchNbrOfNbrAtomSpecs ? $NbrOfNbrAtomSpecsRef->[$Index] : undef;

      # Match neighbor atom specification...
      if (!$NbrAtom->_DoesAtomSpecificationMatch($NbrAtomSpec)) {
	next NBRATOMSPEC;
      }

      # Match central atom to neighbor atom bond specification...
      if (!$CentralAtom->_DoesBondSpecificationMatch($NbrAtom, $NbrBondSpec)) {
	next NBRATOMSPEC;
      }

      # Match any neighbor of neighbor atom specifications...
      if (defined $NbrOfNbrAtomSpec) {
	# Go over the neighbors of central atom skipping the central atom...
	for $NbrOfNbrAtom ($NbrAtom->GetNeighbors($CentralAtom)) {
	  if (!$NbrOfNbrAtom->_DoesAtomSpecificationMatch($NbrOfNbrAtomSpec)) {
	    next NBRATOMSPEC;
	  }
	}
      }

      # It's a match for a neighbor atom specification...
      $NbrSpecAlreadyMatchedRef->{$Index} = $Index;
      $NbrSpecMatchCount++;

      if ($NbrSpecMatchCount == $NbrSpecCount) {
	# It's match...
	$NbrSpecsMatched = 1;
	last NBRATOM;
      }
      # Match next neighbor atom...
      next NBRATOM;
    }
  }
  return ($NbrSpecsMatched, $NbrSpecMatchCount);
}

# Match central atom's missing hydrogens with any unmatched atom and bond
# specifications...
#
sub _MatchAtomNeighborhoodUsingMissingHydrogens {
  my($CentralAtom, $NbrSpecCount, $NbrSpecMatchCount, $NbrSpecAlreadyMatchedRef, $NbrAtomSpecsRef, $NbrBondSpecsRef, $NbrOfNbrAtomSpecsRef) = @_;
  my($Index, $NbrAtom, $NbrAtomSpec, $NbrBondSpec, $NumOfMissingHydrogens, $MissingHydrogensIndex, $NbrSpecsMatched, $AtomSpecMatched, $AtomSpec, $AtomSymbol);

  $NbrSpecsMatched = 0;

  $NumOfMissingHydrogens = $CentralAtom->GetNumOfMissingHydrogens();
  if (($NbrSpecCount - $NbrSpecMatchCount) > $NumOfMissingHydrogens) {
    # It won't match...
    return ($NbrSpecsMatched, $NbrSpecMatchCount);
  }

  MISSINGHYDROGENNBR: for $MissingHydrogensIndex (0 .. ($NumOfMissingHydrogens - 1)) {
    NBRATOMSPEC: for $Index (0 .. ($NbrSpecCount - 1)) {
      if (exists $NbrSpecAlreadyMatchedRef->{$Index}) {
	next NBRATOMSPEC;
      }
      $NbrAtomSpec = $NbrAtomSpecsRef->[$Index];
      $NbrBondSpec = $NbrBondSpecsRef->[$Index];

      $NbrAtomSpec =~ s/ //g;

      # Match neighbor atom specification hydrogen atom symbol...
      $AtomSpecMatched = 0;
      ATOMSPEC: for $AtomSpec (split /\,/, $NbrAtomSpec) {
	($AtomSymbol) = split /\./, $AtomSpec;
	if ($AtomSymbol =~ /^(H|A|\*)$/i) {
	  $AtomSpecMatched = 1;
	  last ATOMSPEC;
	}
      }
      if (!$AtomSpecMatched) {
	next NBRATOMSPEC;
      }

      # Match neighbor atom bond specification to singal bond...
      if (defined $NbrBondSpec) {
	$NbrBondSpec =~ s/ //g;
	if ($NbrBondSpec !~ /^(-|1|s|Single|\~|\*|Any)/i) {
	  next NBRATOMSPEC;
	}
      }

      # It's a match for a neighbor atom specification...
      $NbrSpecAlreadyMatchedRef->{$Index} = $Index;
      $NbrSpecMatchCount++;

      if ($NbrSpecMatchCount == $NbrSpecCount) {
	# It's match...
	$NbrSpecsMatched = 1;
	last MISSINGHYDROGENNBR;
      }
      # Match next missing hydrogen neighbor...
      next MISSINGHYDROGENNBR;
    }
  }

  return ($NbrSpecsMatched, $NbrSpecMatchCount);
}

# Sort atom and bond specifications base on neighbor atom specifications going
# from most to least specific atom specifications.
#
# Atom specifications are sorted at the following two levels:
#
#   o By atom specification count with in each specification going from most specific
#      to least specific, where count is determined by the number of "," in each
#      specification. Wild card containing specifications are considered least specific
#      and end up at the end of the sorted list.
#   o By atomic invariant count with in each sorted list going from most specific to
#      least specific, where count is determined by the number of "." in each atom
#      specification.
#
#A single atom specification,
# without any commas in atom specification, is is considered most specific...
#
sub _SortSpecificationsForAtomNeighborhoodMatch {
  my($This, $NbrAtomSpecsRef, $NbrBondSpecsRef, $NbrOfNbrAtomSpecsRef) = @_;
  my($Index, $NeedToSort, $NumOfNbrAtomSpecs, $NbrAtomSpecCount,  $NbrAtomSpecAtomicInvarintCount, $NbrAtomSpecToMatch, $NbrAtomSpec, $FirstAtomicInvariant, $WildCardInNbrAtomSpec, @NbrAtomSpecs, @NbrAtomSpecAtomicInvariants, @SortedNbrAtomSpecs, @SortedNbrBondSpecs, @SortedNbrOfNbrAtomSpecs, %NbrAtomSpecDataMap);

  $NumOfNbrAtomSpecs = defined $NbrAtomSpecsRef ? scalar @{$NbrAtomSpecsRef} : 0;

  # Figure out whether sorting is necessary...
  $NeedToSort = 0;
  if ($NumOfNbrAtomSpecs > 1) {
    ATOMSPEC: for $NbrAtomSpecToMatch (@{$NbrAtomSpecsRef}) {
      if ($NbrAtomSpecToMatch =~ /(,|\.|A|\*)/i) {
	$NeedToSort = 1;
	last ATOMSPEC;
      }
    }
  }
  if (!$NeedToSort) {
    # Nothing to do...
    return ($NbrAtomSpecsRef, $NbrBondSpecsRef, $NbrOfNbrAtomSpecsRef);
  }

  %NbrAtomSpecDataMap = ();

  for $Index (0 .. ($NumOfNbrAtomSpecs - 1)) {
    $NbrAtomSpecToMatch = $NbrAtomSpecsRef->[$Index];
    $NbrAtomSpecToMatch =~ s/ //g;

    @NbrAtomSpecs = split /\,/, $NbrAtomSpecToMatch;
    $NbrAtomSpecCount = scalar @NbrAtomSpecs;

    # Does neighbor specification contains a wild card in atom symbol specification?
    #
    if ($NbrAtomSpecToMatch =~ /(A|\*)/i) {
      $WildCardInNbrAtomSpec = 0;
      NBRATOMSPEC: for $NbrAtomSpec (@NbrAtomSpecs) {
	($FirstAtomicInvariant) = split /\./, $NbrAtomSpec;
	if ($FirstAtomicInvariant =~ /^!/) {
	  $FirstAtomicInvariant =~ s/^!//;
	}
	$WildCardInNbrAtomSpec = $This->_IsWildCardAtomSymbolAtomicInvariant($FirstAtomicInvariant);
	if ($WildCardInNbrAtomSpec) {
	  last NBRATOMSPEC;
	}
      }
      if ($WildCardInNbrAtomSpec) {
	# Set NbrAtomSpecCount arbitrarily high to make the spec containing wild
	# card last on the sorted list while maintaining its original order in the list...
	$NbrAtomSpecCount = 999;
      }
    }

    if (!exists $NbrAtomSpecDataMap{$NbrAtomSpecCount}) {
      %{$NbrAtomSpecDataMap{$NbrAtomSpecCount}} = ();
    }

    # Use first NbrAtomSpec available in @NbrAtomSpecs to determine atomic invariant count
    # with in each NbrAtomSpecToMatch, as @NbrAtomSpecs derived from $NbrAtomSpecToMatch
    # simply corresponds to a list of possible matches...
    #
    ($NbrAtomSpec) = @NbrAtomSpecs;
    @NbrAtomSpecAtomicInvariants = split /\./, $NbrAtomSpec;
    $NbrAtomSpecAtomicInvarintCount = scalar @NbrAtomSpecAtomicInvariants;

    if (!exists $NbrAtomSpecDataMap{$NbrAtomSpecCount}{$NbrAtomSpecAtomicInvarintCount}) {
      @{$NbrAtomSpecDataMap{$NbrAtomSpecCount}{$NbrAtomSpecAtomicInvarintCount}} = ();
    }
    push @{$NbrAtomSpecDataMap{$NbrAtomSpecCount}{$NbrAtomSpecAtomicInvarintCount}}, $Index;

  }

  @SortedNbrAtomSpecs = (); @SortedNbrBondSpecs = ();
  @SortedNbrOfNbrAtomSpecs = ();

  for $NbrAtomSpecCount ( sort { $a <=> $b } keys %NbrAtomSpecDataMap) {
    for $NbrAtomSpecAtomicInvarintCount ( sort { $b <=> $a } keys %{$NbrAtomSpecDataMap{$NbrAtomSpecCount}}) {
      for $Index (@{$NbrAtomSpecDataMap{$NbrAtomSpecCount}{$NbrAtomSpecAtomicInvarintCount}}) {
	push @SortedNbrAtomSpecs, $NbrAtomSpecsRef->[$Index];
	if (defined $NbrBondSpecsRef) {
	  push @SortedNbrBondSpecs, $NbrBondSpecsRef->[$Index];
	}
	if (defined $NbrOfNbrAtomSpecsRef) {
	  push @SortedNbrOfNbrAtomSpecs, $NbrOfNbrAtomSpecsRef->[$Index];
	}
      }
    }
  }

  return (\@SortedNbrAtomSpecs, defined $NbrBondSpecsRef ? \@SortedNbrBondSpecs : undef, defined $NbrOfNbrAtomSpecsRef ? \@SortedNbrOfNbrAtomSpecs : undef);
}

# Check whether atom matches supported atom specification...
#
sub _DoesAtomSpecificationMatch {
  my($This, $AtomSpecificationToMatch) = @_;
  my($AtomSpecification, $AtomicInvariant, $AtomSpecificationMatched, $AtomicInvariantMatched, $FirstMatch);

  # Anything to match...
  if (!(defined($AtomSpecificationToMatch) && $AtomSpecificationToMatch)) {
    return 1;
  }

  # Take out any spaces...
  $AtomSpecificationToMatch =~ s/ //g;

  # Match specified atom specifications. For multiple atom specifications in a comma delimited string,
  # only one atom specification needs to match for a successful match. It's up to the caller to make
  # sure that the specificaton list is ordered from least to most specific atom specification...
  #
  for $AtomSpecification (split /\,/, $AtomSpecificationToMatch) {
    $AtomSpecificationMatched = 1;
    $FirstMatch = 1;

    # Match all atom symbol atomic invariants...
    ATOMICINVARIANT: for $AtomicInvariant (split /\./, $AtomSpecification) {
      if ($FirstMatch) {
	# Match atom symbol atomic invariant...
	$FirstMatch = 0;
	$AtomicInvariantMatched = $This->_MatchAtomSymbolAtomicInvariant($AtomicInvariant);
      }
      else {
	# Match non atom symbol atomic invariant...
	$AtomicInvariantMatched = $This->_MatchNonAtomSymbolAtomicInvariant($AtomicInvariant);
      }

      if (!$AtomicInvariantMatched) {
	# No need to match other atomic invariants...
	$AtomSpecificationMatched = 0;
	last ATOMICINVARIANT;
      }
    }

    if ($AtomSpecificationMatched) {
      # No need to match other atom specifications...
      return 1;
    }
  }

  # Nothing matched...
  return 0;
}

# Check whether atom matches atom symbol atomic invariant...
#
sub _MatchAtomSymbolAtomicInvariant {
  my($This, $AtomicInvariant) = @_;
  my($NegateMatch, $Status, $AtomicNumber);

  $Status = 0;
  $NegateMatch = 0;

  # Does match needs to be negated?
  if ($AtomicInvariant =~ /^!/) {
    $NegateMatch = 1;
    $AtomicInvariant =~ s/^!//;
  }

  ATOMICINVARIANT: {
    # Any atom match...
    if ($This->_IsWildCardAtomSymbolAtomicInvariant($AtomicInvariant)) {
      $Status = 1;
      last ATOMICINVARIANT;
    }

    # Atomic number match...
    if ($AtomicInvariant =~ /^#/) {
      $AtomicNumber = $AtomicInvariant; $AtomicNumber =~ s/^#//;
      $Status = ($This->{AtomicNumber} == $AtomicNumber) ? 1 : 0;
      last ATOMICINVARIANT;
    }

    # Atom symbol match...
    $Status = ($This->{AtomSymbol} =~ /^$AtomicInvariant$/i) ? 1 : 0;
  }

  if ($NegateMatch) {
    $Status = $Status ? 0 : 1;
  }

  return $Status;
}

# Is it a wild card atom symbol atomic invariant?
#
sub _IsWildCardAtomSymbolAtomicInvariant {
  my($This, $AtomicInvariant) = @_;

  return ($AtomicInvariant =~ /^(A|\*)$/i) ? 1 : 0;
}

# Check whether atom matches non atom symbol atomic invariants...
#
sub _MatchNonAtomSymbolAtomicInvariant {
  my($This, $AtomicInvariant) = @_;
  my($NegateMatch, $Status, $Name, $Value, $UnknownName);

  ($Status, $NegateMatch, $UnknownName) = ('0') x 3;

  # Does match needs to be negated?
  if ($AtomicInvariant =~ /^!/) {
    $NegateMatch = 1;
    $AtomicInvariant =~ s/^!//;
  }

  # Extract atomic invariant name and any value...
  if ($AtomicInvariant =~ /[0-9\*]+/) {
    ($Name, $Value) = $AtomicInvariant =~ /^([a-zA-Z]+)([0-9\-\+\*\>\<\=]+)$/;
  }
  else {
    ($Name, $Value) = ($AtomicInvariant, undef);
  }

  NAME: {
    # Match number of non-hydrogen atom neighbors
    if ($Name =~ /^X$/i) {
      $Status = (defined($Value) && $This->GetNumOfNonHydrogenAtomNeighbors() == $Value) ? 1 : 0;
      last NAME;
    }

    # Match total number of atom neighbors including missing hydrogens...
    if ($Name =~ /^T$/i) {
      $Status = (defined($Value) && ($This->GetNumOfNonHydrogenAtomNeighbors() + $This->GetNumOfHydrogens()) == $Value) ? 1 : 0;
      last NAME;
    }

    # Match formal charge...
    if ($Name =~ /^FC$/i) {
      my $FormalCharge = $This->GetFormalCharge();
      $Status = $This->_MatchNonAtomSymbolAtomicInvariantValue($FormalCharge, $Value);
      last NAME;
    }

    # Match aromatic annotation indicating whether atom is aromatic...
    if ($Name =~ /^Ar$/i) {
      $Status = $This->IsAromatic() ? 1 : 0;
      last NAME;
    }

    # Match number of implicit and explicit hydrogens...
    if ($Name =~ /^H$/i) {
      $Status = (defined($Value) && ($This->GetNumOfHydrogens() == $Value)) ? 1 : 0;
      last NAME;
    }

    # Match ring atom annotation indicating whether atom is in ring...
    if ($Name =~ /^RA$/i) {
      $Status = defined($Value) ? $This->IsInRingOfSize($Value) : ($This->IsInRing() ? 1 : 0);
      last NAME;
    }

    # Match number of rings for atom..
    if ($Name =~ /^TR$/i) {
      $Status = (defined($Value) && ($Value == $This->GetNumOfRings())) ? 1 : 0;
      last NAME;
    }

    # Match sum of bond orders to non-hydrogen atom neighbors...
    if ($Name =~ /^BO$/i) {
      $Status = (defined($Value) && $This->GetSumOfBondOrdersToNonHydrogenAtoms() == $Value) ? 1 : 0;
      last NAME;
    }

    # Match largest bond order of non-hydrogen atom neighbors...
    if ($Name =~ /^LBO$/i) {
      $Status = (defined($Value) && $This->GetLargestBondOrderToNonHydrogenAtoms() == $Value) ? 1 : 0;
      last NAME;
    }

    # Match number of single bonds to non-hydrogen atom neighbors...
    if ($Name =~ /^SB$/i) {
      $Status = (defined($Value) && $This->GetNumOfSingleBondsToNonHydrogenAtoms() == $Value) ? 1 : 0;
      last NAME;
    }

    # Match total number of single bonds to atom neighbors including missing and explicit hydrogens...
    if ($Name =~ /^TSB$/i) {
      $Status = (defined($Value) && ($This->GetNumOfSingleBondsToNonHydrogenAtoms() + $This->GetNumOfHydrogens()) == $Value) ? 1 : 0;
      last NAME;
    }

    # Match number of double bonds to non-hydrogen atom neighbors...
    if ($Name =~ /^DB$/i) {
      $Status = (defined($Value) && $This->GetNumOfDoubleBondsToNonHydrogenAtoms() == $Value) ? 1 : 0;
      last NAME;
    }

    # Match number of triple bonds to non-hydrogen atom neighbors...
    if ($Name =~ /^TB$/i) {
      $Status = (defined($Value) && $This->GetNumOfTripleBondsToNonHydrogenAtoms() == $Value) ? 1 : 0;
      last NAME;
    }

    # Match number of aromatic bonds to non-hydrogen atom neighbors...
    if ($Name =~ /^AB$/i) {
      $Status = (defined($Value) && $This->GetNumOfAromaticBondsToNonHydrogenAtoms() == $Value) ? 1 : 0;
      last NAME;
    }


    # Match mass number indicating isotope other than most abundant isotope...
    if ($Name =~ /^MN$/i) {
      $Status = (defined($Value) && $This->GetMassNumber() == $Value) ? 1 : 0;
      last NAME;
    }

    # Match spin multiplicity...
    if ($Name =~ /^SM$/i) {
      my $SpinMultiplicity = $This->GetSpinMultiplicity();
      if (!defined $SpinMultiplicity) { $SpinMultiplicity = 0; }
      $Status = (defined($Value) && defined($SpinMultiplicity) && $Value == $SpinMultiplicity) ? 1 : 0;
      last NAME;
    }

    $UnknownName = 1;
    carp "Warning: ${ClassName}->_MatchNonAtomSymbolAtomicInvariant: Unknown atomic invariant $AtomicInvariant...";
  }

  if (!$UnknownName) {
    if ($NegateMatch) {
      $Status = $Status ? 0 : 1;
    }
  }

  return $Status;
}

# Match atomic invariant value...
#
# Specified value format:
#   . +* : Any positive value
#   . -* : Any negative value
#   . >ValidNumber or >=ValidNumber
#   . <ValidNumber or <=ValidNumber
#   . Any valid number
#
sub _MatchNonAtomSymbolAtomicInvariantValue {
  my($This, $TargetValue, $SpecifiedValue) = @_;
  my($Status);

  $Status = 0;

  if (!(defined($TargetValue) && defined($SpecifiedValue))) {
    return $Status;
  }

  VALUE: {
    if ($SpecifiedValue =~ /^\+\*/) {
      $Status = ($TargetValue > 0) ? 1 : 0;
      last VALUE;
    }
    if ($SpecifiedValue =~ /^\-\*/) {
      $Status = ($TargetValue < 0) ? 1 : 0;
      last VALUE;
    }
    if ($SpecifiedValue =~ /^>/) {
      if ($SpecifiedValue =~ /^>=/) {
	$SpecifiedValue =~ s/^>=//;
	$Status = ($SpecifiedValue >= $TargetValue) ? 1 : 0;
      }
      else {
	$SpecifiedValue =~ s/^>//;
	$Status = ($SpecifiedValue > $TargetValue) ? 1 : 0;
      }
      last VALUE;
    }
    if ($SpecifiedValue =~ /^</) {
      if ($SpecifiedValue =~ /^<=/) {
	$SpecifiedValue =~ s/^<=//;
	$Status = ($SpecifiedValue <= $TargetValue) ? 1 : 0;
      }
      else {
	$SpecifiedValue =~ s/^<//;
	$Status = ($SpecifiedValue < $TargetValue) ? 1 : 0;
      }
      last VALUE;
    }
    # Default is do perform an equality match...
    $Status = ($SpecifiedValue == $TargetValue) ? 1 : 0;
  }

  return $Status;
}

# Check whether atoms match  bond specifications...
#
sub _DoesBondSpecificationMatch {
  my($This, $BondedAtom, $BondSpecificationToMatch) = @_;
  my($BondSpecification, $BondSymbolSpecification, $BondSpecificationMatched);

  # Anything to match...
  if (!(defined($BondSpecificationToMatch) && $BondSpecificationToMatch)) {
    return 1;
  }

  # Take out any spaces...
  $BondSpecificationToMatch =~ s/ //g;

  # Match specified bond specifications. For multiple bond specifications in a comma delimited string,
  # only one bond specification needs to match for a successful match...
  #
  for $BondSpecification (split /\,/, $BondSpecificationToMatch) {
    $BondSpecificationMatched = 1;

    # Match all specified bond symbol specifications...
    BONDSYMBOL: for $BondSymbolSpecification (split /\./, $BondSpecification) {
      if (!$This->_MatchBondSymbolSpecification($BondedAtom, $BondSymbolSpecification)) {
	# No need to match other bond symbol specifications...
	$BondSpecificationMatched = 0;
	last BONDSYMBOL;
      }
    }
    if ($BondSpecificationMatched) {
      # No need to try matching other bond specifications...
      return 1;
    }
  }

  # Nothing matched...
  return 0;
}

# Check whether atoms match  bond symbol specification...
#
sub _MatchBondSymbolSpecification {
  my($This, $BondedAtom, $BondSymbolSpecification) = @_;
  my($NegateMatch, $Status, $Bond, $BondSymbol, $UnknownBondSymbol);

  ($Status, $NegateMatch, $UnknownBondSymbol) = ('0') x 3;

  # Does match needs to be negated?
  if ($BondSymbolSpecification =~ /^!/) {
    $NegateMatch = 1;
    $BondSymbolSpecification =~ s/^!//;
  }
  $BondSymbol = $BondSymbolSpecification;
  $Bond = $This->GetBondToAtom($BondedAtom);

  BONDSYMBOL: {
    if ($BondSymbol =~ /^(-|1|s|Single)$/i) { $Status = $Bond->IsSingle() ? 1 : 0; last BONDSYMBOL; }
    if ($BondSymbol =~ /^(=|2|d|Double)$/i) { $Status = $Bond->IsDouble() ? 1 : 0; last BONDSYMBOL; }
    if ($BondSymbol =~ /^(#|3|t|Triple)$/i) { $Status = $Bond->IsTriple() ? 1 : 0; last BONDSYMBOL; }
    if ($BondSymbol =~ /^(:|a|Ar|Aromatic)$/i) { $Status = $Bond->IsAromatic() ? 1 : 0; last BONDSYMBOL; }

    if ($BondSymbol =~ /^(\@|RB|Ring)$/i) { $Status = $Bond->IsInRing() ? 1 : 0; last BONDSYMBOL; }

    if ($BondSymbol =~ /^(\~|\*|Any)$/i) { $Status = 1; last BONDSYMBOL; }

    $UnknownBondSymbol = 1;
    carp "Warning: ${ClassName}->_MatchBondSpecification: Unknown bond specification $BondSymbolSpecification...";
  }

  if (!$UnknownBondSymbol) {
    if ($NegateMatch) {
      $Status = $Status ? 0 : 1;
    }
  }

  return $Status;
}

# Is it a saturated atom?
#
sub IsSaturated {
  my($This) = @_;

  return !$This->IsUnsaturated();
}

# Is it an unsaturated atom containing at least one non-single bond?
#
sub IsUnsaturated {
  my($This) = @_;
  my($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds);

  ($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds) = $This->GetNumOfBondTypesToNonHydrogenAtoms();

  return ($NumOfDoubleBonds || $NumOfTripleBonds || $NumOfAromaticBonds) ? 1 : 0;
}

# Is atom in a ring?
#
sub IsInRing {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_IsAtomInRing($This);
}

# Is atom not in a ring?
#
sub IsNotInRing {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_IsAtomNotInRing($This);
}

# Is atom only in one ring?
#
sub IsOnlyInOneRing {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_IsAtomInOnlyOneRing($This);
}

# Is atom in a ring of specific size?
#
sub IsInRingOfSize {
  my($This, $RingSize) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_IsAtomInRingOfSize($This, $RingSize);
}

# Get size of smallest ring containing the atom...
#
sub GetSizeOfSmallestRing {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetSizeOfSmallestAtomRing($This);
}

# Get size of largest ring containing the atom...
#
sub GetSizeOfLargestRing {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetSizeOfLargestAtomRing($This);
}

# Get number of  rings containing the atom...
#
sub GetNumOfRings {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetNumOfAtomRings($This);
}

# Get number of  rings with odd size containing the atom...
#
sub GetNumOfRingsWithOddSize {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetNumOfAtomRingsWithOddSize($This);
}

# Get number of  rings with even size containing the atom...
#
sub GetNumOfRingsWithEvenSize {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetNumOfAtomRingsWithEvenSize($This);
}

# Get number of  rings with specified size containing the atom...
#
sub GetNumOfRingsWithSize {
  my($This, $RingSize) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetNumOfAtomRingsWithSize($This, $RingSize);

}

# Get number of  rings with size less than specified containing the atom...
#
sub GetNumOfRingsWithSizeLessThan {
  my($This, $RingSize) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetNumOfAtomRingsWithSizeLessThan($This, $RingSize);
}

# Get number of  rings with size greater than specified size containing the atom...
#
sub GetNumOfRingsWithSizeGreaterThan {
  my($This, $RingSize) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetNumOfAtomRingsWithSizeGreaterThan($This, $RingSize);
}

# Get all rings an array of references to arrays containing ring atoms...
#
sub GetRings {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetAtomRings($This);
}

# Get smallest ring as an array containing ring atoms...
#
sub GetSmallestRing {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetSmallestAtomRing($This);
}

# Get largest ring as an array containing ring atoms...
#
sub GetLargestRing {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetLargestAtomRing($This);
}

# Get odd size rings an array of references to arrays containing ring atoms...
#
sub GetRingsWithOddSize {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetAtomRingsWithOddSize($This);
}

# Get even size rings an array of references to arrays containing ring atoms...
#
sub GetRingsWithEvenSize {
  my($This) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetAtomRingsWithEvenSize($This);
}

# Get rings with specified size as an array of references to arrays containing ring atoms...
#
sub GetRingsWithSize {
  my($This, $RingSize) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetAtomRingsWithSize($This, $RingSize);
}

# Get rings with size less than specfied size as an array of references to arrays containing ring atoms...
#
sub GetRingsWithSizeLessThan {
  my($This, $RingSize) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetAtomRingsWithSizeLessThan($This, $RingSize);
}

# Get rings with size greater than specfied size as an array of references to arrays containing ring atoms...
#
sub GetRingsWithSizeGreaterThan {
  my($This, $RingSize) = @_;

  # Is this atom in a molecule?
  if (!$This->HasProperty('Molecule')) {
    return undef;
  }
  my($Molecule);
  $Molecule = $This->GetProperty('Molecule');

  return $Molecule->_GetAtomRingsWithSizeGreaterThan($This, $RingSize);
}

# Get next object ID...
sub _GetNewObjectID {
  $ObjectID++;
  return $ObjectID;
}

# Return a string containing vertices, edges and other properties...
sub StringifyAtom {
  my($This) = @_;
  my($AtomString, $ID, $Name, $AtomSymbol, $AtomicNumber, $XYZVector, $AtomicWeight, $ExactMass, $NumOfNeighbors, $NumOfBonds, $Valence, $MissingHydrogens, $TotalHydrogens, $ImplicitHydrogens, $ExplicitHydrogens, $FormalCharge, $Charge, $SpinMultiplicity, $FreeRadicalElectrons, $StereoCenter, $StereoCenterStatus, $StereoChemistry, $StereochemistryString, $RingAtom, $NumOfRings, $AromaticAtom);

  $ID = $This->GetID();
  $Name = $This->GetName();
  $AtomSymbol = $This->GetAtomSymbol();
  $AtomicNumber = $This->GetAtomicNumber();
  $XYZVector = $This->GetXYZVector();

  $AtomicWeight = $This->GetAtomicWeight();
  if (!defined $AtomicWeight) {
    $AtomicWeight = 'undefined';
  }
  $ExactMass = $This->GetExactMass();
  if (!defined $ExactMass) {
    $ExactMass = 'undefined';
  }
  $NumOfNeighbors = $This->GetNumOfNeighbors();
  if (!defined $NumOfNeighbors) {
    $NumOfNeighbors = 'undefined';
  }
  $NumOfBonds = $This->GetNumOfBonds();
  if (!defined $NumOfBonds) {
    $NumOfBonds = 'undefined';
  }
  $Valence = $This->GetValence();
  if (!defined $Valence) {
    $Valence = 'undefined';
  }

  $MissingHydrogens = $This->GetNumOfMissingHydrogens();
  if (!defined $MissingHydrogens) {
    $MissingHydrogens = 'undefined';
  }
  $TotalHydrogens = $This->GetNumOfHydrogens();
  if (!defined $TotalHydrogens) {
    $TotalHydrogens = 'undefined';
  }
  $ImplicitHydrogens = $This->GetNumOfImplicitHydrogens();
  if (!defined $ImplicitHydrogens) {
    $ImplicitHydrogens = 'undefined';
  }
  $ExplicitHydrogens = $This->GetNumOfExplicitHydrogens();
  if (!defined $ExplicitHydrogens) {
    $ExplicitHydrogens = 'undefined';
  }

  $FormalCharge = $This->GetFormalCharge();
  if (!defined $FormalCharge) {
    $FormalCharge = 'undefined';
  }
  $Charge = $This->GetCharge();
  if (!defined $Charge) {
    $Charge = 'undefined';
  }

  $SpinMultiplicity = $This->GetSpinMultiplicity();
  if (!defined $SpinMultiplicity) {
    $SpinMultiplicity = 'undefined';
  }
  $FreeRadicalElectrons = $This->GetFreeRadicalElectrons();
  if (!defined $FreeRadicalElectrons) {
    $FreeRadicalElectrons = 'undefined';
  }

  $RingAtom = $This->IsInRing();
  if (defined $RingAtom) {
    $RingAtom = $RingAtom  ? 'Yes' : 'No';
    $NumOfRings = $This->GetNumOfRings();
  }
  else {
    $RingAtom = 'undefined';
    $NumOfRings = 'undefined';
  }

  $AromaticAtom = $This->GetAromatic();
  if (defined $AromaticAtom) {
    $AromaticAtom = $AromaticAtom  ? 'Yes' : 'No';
  }
  else {
    $AromaticAtom = 'undefined';
  }

  $StereochemistryString = '';
  $StereoCenter = $This->GetStereoCenter();
  if (defined $StereoCenter) {
    $StereoCenterStatus = $This->IsStereoCenter() ? 'Yes' : 'No';
    $StereoChemistry = $This->GetStereochemistry();
    if (!defined $StereoChemistry) {
      $StereoChemistry = 'undefined';
    }
    $StereochemistryString = "StereoCenter: $StereoCenterStatus; Stereochemistry: $StereoChemistry";
  }

  $AtomString = "Atom: ID: $ID; Name: \"$Name\"; AtomSymbol: \"$AtomSymbol\"; AtomicNumber: $AtomicNumber; XYZ: $XYZVector; AtomicWeight: $AtomicWeight; ExactMass: $ExactMass; NumOfNeighbors: $NumOfNeighbors;  NumOfBonds: $NumOfBonds; Valence: $Valence; MissingHydrogens: $MissingHydrogens; TotalHydrogens: $TotalHydrogens; ImplicitHydrogens: $ImplicitHydrogens; ExplicitHydrogens: $ExplicitHydrogens; FormalCharge: $FormalCharge; Charge: $Charge; SpinMultiplicity: $SpinMultiplicity; FreeRadicalElectrons: $FreeRadicalElectrons; RingAtom: $RingAtom; NumOfAtomRings: $NumOfRings; AromaticAtom: $AromaticAtom";

  if ($StereochemistryString) {
    $AtomString .= "; $StereochemistryString";
  }

  return $AtomString;
}

# Load appropriate atom data files from <MayaChemTools>/lib directory used by various
# object methods in the current class...
#
sub _LoadAtomClassData {
  my($MayaChemToolsLibDir);

  $MayaChemToolsLibDir = GetMayaChemToolsLibDirName();

  _LoadAtomValenceModelData($MayaChemToolsLibDir);

}

#
# Load data for supported valence models...
#
sub _LoadAtomValenceModelData {
  my($MayaChemToolsLibDir) = @_;
  my($MDLValenceModelDataFile, $DaylightValenceModelDataFile);

  %MDLValenceModelDataMap = ();
  %DaylightValenceModelDataMap = ();

  $MDLValenceModelDataFile = $MayaChemToolsLibDir . "/data/MDLValenceModelData.csv";
  $DaylightValenceModelDataFile = $MayaChemToolsLibDir . "/data/DaylightValenceModelData.csv";

  if (! -e "$MDLValenceModelDataFile") {
    croak "Error: ${ClassName}::_LoadAtomValenceModelData: MayaChemTools package file, $MDLValenceModelDataFile, is missing: Possible installation problems...";
  }

  if (! -e "$DaylightValenceModelDataFile") {
    croak "Error: ${ClassName}::_LoadAtomValenceModelData: MayaChemTools package file, $DaylightValenceModelDataFile, is missing: Possible installation problems...";
  }

  _LoadValenceModelDataFile($MDLValenceModelDataFile, \%MDLValenceModelDataMap);
  _LoadValenceModelDataFile($DaylightValenceModelDataFile, \%DaylightValenceModelDataMap);

}

#
# Load valence model data file...
#
sub _LoadValenceModelDataFile {
  my($DataFile, $DataMapRef) = @_;

  # File format:
  #
  # "AtomicNumber","ElementSymbol","FormalCharge","CommomValences"
  #
  my($InDelim, $Line, $NumOfCols, @ColLabels, @LineWords);

  $InDelim = "\,";
  open DATAFILE, "$DataFile" or croak "Couldn't open $DataFile: $! ...";

  # Skip lines up to column labels...
  LINE: while ($Line = GetTextLine(\*DATAFILE)) {
    if ($Line !~ /^#/) {
      last LINE;
    }
  }
  @ColLabels= quotewords($InDelim, 0, $Line);
  $NumOfCols = @ColLabels;

  my($AtomicNumber, $ElementSymbol, $FormalCharge, $CommonValences);

  # Process element data...
  LINE: while ($Line = GetTextLine(\*DATAFILE)) {
    if ($Line =~ /^#/) {
      next LINE;
    }
    @LineWords = ();
    @LineWords = quotewords($InDelim, 0, $Line);
    if (@LineWords != $NumOfCols) {
      croak "Error: ${ClassName}::_LoadValenceModelDataFile: The number of data fields, @LineWords, in $DataFile must be $NumOfCols.\nLine: $Line...";
    }

    ($AtomicNumber, $ElementSymbol, $FormalCharge, $CommonValences) = @LineWords;

    if (exists $DataMapRef->{$AtomicNumber}) {
      # Additional data for an element...
      if (exists $DataMapRef->{$AtomicNumber}{$FormalCharge}) {
	# Duplicate data entries for an element...
	carp "Warning: ${ClassName}::_LoadValenceModelDataFile: Ignoring valence data for element with atomic number $AtomicNumber and formal charge $FormalCharge in data file $DataFile: It has already been loaded.\nLine: $Line...";
	next LINE;
      }
    }
    else {
      # Data for a new element...
      %{$DataMapRef->{$AtomicNumber}} = ();
    }

    %{$DataMapRef->{$AtomicNumber}{$FormalCharge}} = ();
    $DataMapRef->{$AtomicNumber}{$FormalCharge}{ElementSymbol} = $ElementSymbol;

    @{$DataMapRef->{$AtomicNumber}{$FormalCharge}{CommonValences}} = ();
    @{$DataMapRef->{$AtomicNumber}{$FormalCharge}{CommonValences}} = sort { $a <=> $b } split /\,/, $CommonValences;
  }
  close DATAFILE;
}

1;

__END__

=head1 NAME

Atom

=head1 SYNOPSIS

use Atom;

=head1 DESCRIPTION

B<Atom> class provides the following methods:

new, AddHydrogens, Copy, DeleteAtom, DeleteHydrogens, DoesAtomNeighborhoodMatch,
GetAtomicInvariantValue, GetAtomicWeight, GetBondToAtom, GetBonds,
GetBondsToHeavyAtoms, GetBondsToHydrogenAtoms, GetBondsToNonHydrogenAtoms,
GetExactMass, GetExplicitHydrogens, GetFormalCharge, GetFreeRadicalElectrons,
GetGroupNumber, GetHeavyAtomNeighbors, GetHeavyAtomNeighborsAtomInformation,
GetHeavyAtomNeighborsBondformation, GetHighestCommonValence,
GetHydrogenAtomNeighbors, GetHydrogens, GetImplicitHydrogens, GetLargestBondOrder,
GetLargestBondOrderToHeavyAtoms, GetLargestBondOrderToNonHydrogenAtoms,
GetLargestRing, GetLowestCommonValence, GetMassNumber, GetMissingHydrogens,
GetNeighbors, GetNeighborsUsingAtomSpecification, GetNonHydrogenAtomNeighbors,
GetNonHydrogenAtomNeighborsAtomInformation,
GetNonHydrogenAtomNeighborsBondInformation, GetNonHydrogenNeighborOfHydrogenAtom,
GetNumOfAromaticBondsToHeavyAtoms, GetNumOfAromaticBondsToNonHydrogenAtoms,
GetNumOfBondTypesToHeavyAtoms, GetNumOfBondTypesToNonHydrogenAtoms, GetNumOfBonds,
GetNumOfBondsToHeavyAtoms, GetNumOfBondsToHydrogenAtoms,
GetNumOfBondsToNonHydrogenAtoms, GetNumOfDoubleBondsToHeavyAtoms,
GetNumOfBondsAvailableForHeavyAtoms, GetNumOfBondsAvailableForNonHydrogenAtoms,
GetNumOfDoubleBondsToNonHydrogenAtoms, GetNumOfExplicitHydrogens,
GetNumOfHeavyAtomNeighbors, GetNumOfHydrogenAtomNeighbors, GetNumOfHydrogens,
GetNumOfImplicitHydrogens, GetNumOfMissingHydrogens, GetNumOfNeighbors,
GetNumOfNonHydrogenAtomNeighbors, GetNumOfRings, GetNumOfRingsWithEvenSize,
GetNumOfRingsWithOddSize, GetNumOfRingsWithSize, GetNumOfRingsWithSizeGreaterThan,
GetNumOfRingsWithSizeLessThan, GetNumOfSigmaAndPiBondsToHeavyAtoms,
GetNumOfSigmaAndPiBondsToNonHydrogenAtoms, GetNumOfSingleBondsToHeavyAtoms,
GetNumOfSingleBondsToNonHydrogenAtoms, GetNumOfTripleBondsToHeavyAtoms,
GetNumOfTripleBondsToNonHydrogenAtoms, GetPeriodNumber,
GetPotentialTotalCommonValence, GetRings, GetRingsWithEvenSize,
GetRingsWithOddSize, GetRingsWithSize, GetRingsWithSizeGreaterThan,
GetRingsWithSizeLessThan, GetSizeOfLargestRing, GetSizeOfSmallestRing,
GetSmallestRing, GetSpinMultiplicity, GetSumOfBondOrders,
GetSumOfBondOrdersToHeavyAtoms, GetSumOfBondOrdersToHydrogenAtoms,
GetSumOfBondOrdersToNonHydrogenAtoms, GetValence, GetValenceElectrons,
GetValenceFreeElectrons, GetX, GetXYZ, GetXYZVector, GetY, GetZ, IsAmideCarbon,
IsAmideNitrogen, IsAromatic, IsArsenic, IsBondedToAtom, IsBromine, IsCarbon, IsCarboxylCarbon,
IsCarboxylOxygen, IsCarboxylateCarbon, IsCarboxylateOxygen, IsChlorine,
IsFluorine, IsFunctionalClassType, IsGuadiniumCarbon, IsGuadiniumNitrogen,
IsHBondAcceptor, IsHBondDonor, IsHalogen, IsHeteroAtom, IsHydrogen,
IsHydrogenBondAcceptor, IsHydrogenBondDonor, IsHydrophobic, IsInRing,
IsInRingOfSize, IsIodine, IsIsotope, IsLipophilic, IsMetallic,
IsNegativelyIonizable, IsNitrogen, IsNonCarbonOrHydrogen, IsNotInRing,
IsOnlyInOneRing, IsOxygen, IsPhosphateOxygen, IsPhosphatePhosphorus, IsPhosphorus,
IsPolarAtom, IsPolarHydrogen, IsPositivelyIonizable, IsSaturated, IsSelenium,
IsSilicon, IsStereoCenter, IsSulfur, IsSulphur, IsTellurium, IsTerminal,
IsTopologicalPharmacophoreType, IsUnsaturated, SetAtomSymbol, SetAtomicNumber,
SetExplicitHydrogens, SetMassNumber, SetStereoCenter, SetStereochemistry,
SetX, SetXYZ, SetY, SetZ, StringifyAtom

B<Atom> class is derived from B<ObjectProperty> base class which provides methods not explicitly
defined in B<Atom> or B<ObjectProperty> class using Perl's AUTOLOAD functionality. These methods
are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

=head2 METHODS

=over 4

=item B<new>

    $NewAtom = new Atom([%PropertyNameAndValues]);

Using specified I<Atom> property names and values hash, B<new> method creates a new object
and returns a reference to newly created B<Atom> object. By default, following properties are
initialized:

    ID = SequentialObjectID
    Name = "Atom <SequentialObjectID>"
    AtomSymbol = ""
    AtomicNumber = 0
    XYZ = ZeroVector

Except for I<ID> property, all other default properties and other additional properties can
be set during invocation of this method.

Examples:

    $Atom = new Atom();
    $CarbonAtom = new Atom('AtomSymbol' => 'C', 'XYZ' => (0.0, 1.0,
                  0.0));
    $OxygenAtom = new Atom('AtomName' => 'Oxygen', AtomSymbol' => 'O',
                  'XYZ' => (1.0, 1.0, 1.0));

=item B<AddHydrogens>

    $NumOfHydrogensAdded = $Atom->AddHydrogens();

Adds hydrogens to an B<Atom> present in a B<Molecule> object and returns
the number of added hydrogens. The current release of MayaChemTools doesn't
assign hydrogen positions.

=item B<Copy>

    $AtomCopy = $Atom->Copy();

Copy I<Atom> and its associated data using B<Storable::dclone> and return a new
B<Atom> object.

=item B<DeleteAtom>

    $Atom->DeleteAtom();

Delete I<Atom> from a molecule.

=item B<DoesAtomNeighborhoodMatch>

    $Status = $Atom->DoesAtomNeighborhoodMatch($CentralAtomSpec);
    $Status = $Atom->DoesAtomNeighborhoodMatch($CentralAtomSpec,
                              $NbrAtomSpecsRef);
    $Status = $Atom->DoesAtomNeighborhoodMatch($CentralAtomSpec,
                              $NbrAtomSpecsRef, $AllowedNbrBondSpecsRef);
    $Status = $Atom->DoesAtomNeighborhoodMatch($CentralAtomSpec,
                              $NbrAtomSpecsRef, $NbrBondSpecsRef,
                              $AllowedNbrOfNbrAtomSpecsRef);

Returns 1 or 0 based on whether atom matches central atom and its neighborhood
using specified atom and bonds specifications. Neighborhood atom and bond specifications
are specified as array references containing neighbor atom and bond specifications.

Let:

    AS = Atom symbol corresponding to element symbol, atomic number (#n)
         or any atom (A)

    X<n>  = Number of non-hydrogen atom neighbors or heavy atoms
            attached to atom
    T<n>  = Total number of atom neighbors including implicit and explicit
            hydrogens
    BO<n> = Sum of bond orders to non-hydrogen atom neighbors or heavy
            atoms attached to atom
    LBO<n> = Largest bond order of non-hydrogen atom neighbors or heavy
             atoms attached to atom
    SB<n> = Number of single bonds to non-hydrogen atom neighbors or
            heavy atoms attached to atom
    TSB<n> = Total number of single bonds to atom neighbors including implicit
             and explicit hydrogens
    DB<n> = Number of double bonds to non-hydrogen atom neighbors or
            heavy atoms attached to atom
    TB<n> = Number of triple bonds to non-hydrogen atom neighbors or
            heavy atoms attached to atom
    AB<n> = Number of aromatic bonds to non-hydrogen atom neighbors or
            heavy atoms attached to atom
    H<n>   = Number of implicit and explicit hydrogens for atom
    Ar     = Aromatic annotation indicating whether atom is aromatic
    RA or RA<n>  = Ring atom annotation indicating whether atom
                   is a ring
    TR<n>  = Total number of rings containing atom
    FC<+n/-n> = Formal charge assigned to atom
    MN<n> = Mass number indicating isotope other than most abundant isotope
    SM<n> = Spin multiplicity of atom. Possible values: 1 (singlet),
            2 (doublet) or 3 (triplet)

Then, atom specification corresponds to:

    AS.X<n>.T<n>.BO<n>.LBO<n>.<SB><n>.TSB<n>.<DB><n>.<TB><n>.AB<n>.H<n>.Ar.
    RA<n>.TR<n>FC<+n/-n>.MN<n>.SM<n>

Except for AS which is a required atomic invariant in atom specification, all other atomic invariants are
optional. For an atom specification to match an atom, the values of all specified atomic invariants must
match. Exclamation in from of atomic invariant can be used to negate its effect during the match.

For I<FC> value matching, the following value operators are also supported:

    o +* : Any positive value
    o -* : Any negative value
    o > ValidNumber or >= ValidNumber
    o < ValidNumber or <= ValidNumber

A comma delimited atom specification string is used to match any one of the specified atom specification.

Notes:

    o During atom specification match to an atom, the first atomic invariant is always assumed to
      atom symbol.

Examples:

    o ('N', 'N', 'N')
    o ('N.FC0', 'N.FC0', 'N,N.FC+1.H1')
    o ('N.H2', 'N.H2', 'N.H1')
    o ('C,N', '!N', '!H')
    o ('C,N', 'N.Ar', 'N.R5')

Let:

    -|1|s|Single = Single bond
    =|2|d|Double = Double bond
    #|3|t|Triple  = Triple bond
    :|1.5|a|Ar|Aromatic = Aromatic bond

    @|RB|Ring = Ring bond
    ~|*|Any = Any bond

Then, bond specification corresponds to:

    -.:
    =.@
    Double.Aromatic

For a bond specification to match bond between two atoms, the values of all specified bond symbols must
match. Exclamation in from of bond symbol can be used to negate its effect during the match.

A comma delimited bond specification string is used to match any one of the specified atom specification.

Notes:

    o During atom neighborhood match for central atom neighborhood atom and bond specifications,
      implicit or missing hydrogens are automatically checked for any matches to unmatched
      specifications.

Examples:


    Aromatic carbon in a 5 membered ring:
                              $Atom->DoesAtomNeighborhoodMatch('C.Ar.RA5');

    AcetylenicCarbon: $Atom->DoesAtomNeighborhoodMatch('C.T2.TB1'); or
                       $Atom->DoesAtomNeighborhoodMatch('C.T2.TB1',
                              ['*', '*'], ['#', '-']);

    GuadiniumCarbon: $Atom->DoesAtomNeighborhoodMatch('C.X3.BO4',
                              ['N.FC0', 'N.FC0', 'N.FC0,N.FC+1'],
                              ['-', '-', '='],
                              ['C,H', 'C,H', 'C,H']);

    AmideCarbon: $Atom->DoesAtomNeighborhoodMatch('C.X3.BO4,C.X2.BO3',
                              ['C,H', 'O', 'N'],
                              ['-', '=', '-'],
                              ['C,H', 'C', 'C,H,N,O,S,P,F,Cl,Br,I']);

    CarboxylCarbon: $Atom->DoesAtomNeighborhoodMatch('C.X3.BO4,C.X2.BO3',
                              ['C,H', 'O', 'O.X1.FC0'],
                              ['-', '=', '-'],
                              ['C,H', 'C', 'C']);

    CarboxylateCarbon: $Atom->DoesAtomNeighborhoodMatch('C.X3.BO4,C.X2.BO3',
                              ['C,H', 'O', 'O.X1.FC-1'],
                              ['-', '=', '-'],
                              ['C,H', 'C', 'C']);


=item B<DeleteHydrogens>

    $NumOfHydrogensDeleted = $Atom->AddHydrogens();

Delete hydrogens from an B<Atom> present in a B<Molecule> object and returns
the number of deleted hydrogens.

=item B<GetAtomicInvariantValue>

    $Value = $Atom->GetAtomicInvariantValue($AtomicInvariant);

Returns atomic invariant value for a specified I<AtomicInvariant>. The current release
of MayaChemTools supports following abbreviations and descriptive names for
I<AtomicInvarints>:

    AS : Atom or element symbol
    X : NumOfNonHydrogenAtomNeighbors or NumOfHeavyAtomNeighbors
    T : TotalNumOfAtomNeighbors
    BO : SumOfBondOrdersToNonHydrogenAtoms or SumOfBondOrdersToHeavyAtoms
    LBO : LargestBondOrderToNonHydrogenAtoms or LargestBondOrderToHeavyAtoms
    SB :  NumOfSingleBondsToNonHydrogenAtoms or NumOfSingleBondsToHeavyAtoms
    TSB :  TotalNumOfSingleBonds
    DB : NumOfDoubleBondsToNonHydrogenAtoms or NumOfDoubleBondsToHeavyAtoms
    TB : NumOfTripleBondsToNonHydrogenAtoms or NumOfTripleBondsToHeavyAtoms
    AB : NumOfAromaticBondsToNonHydrogenAtoms or NumOfAromaticBondsToHeavyAtoms
    H :  NumOfImplicitAndExplicitHydrogens
    Ar : Aromatic
    Str : Stereochemistry
    RA : RingAtom
    FC : FormalCharge
    AN : AtomicNumber
    AM : AtomicMass
    MN : MassNumber
    SM : SpinMultiplicity

=item B<GetAtomicWeight>

    $Value = $Aom->GetAtomicWeight();

Returns atomic weight of an B<Atom> which corresponds to either explicity set I<AtomicWeight>
atom property or atomic weight of the corresponding element in the periodic table available by
B<PeriodicTable> module.

=item B<GetBondToAtom>

    $Bond = $Atom->GetBondToAtom($OtherAtom);

Returns a B<Bond> object corresponding to bond between I<Atom> and I<OtherAtom> in
a molecule.

=item B<GetBonds>

    @Bonds = $Aoto->GetBonds();

Returns an array of B<Bond> objects corresponding to all bonds from I<Atom> to other
bonded atoms in a molecule.

=item B<GetBondsToHeavyAtoms>

    @Bonds = $Atom->GetBondsToHeavyAtoms();

Returns an array of B<Bond> objects corresponding to bonds from I<Atom> to other bonded
non-hydrogen atoms in a molecule.

=item B<GetBondsToHydrogenAtoms>

    @Bonds = $Atom->GetBondsToHydrogenAtoms();

Returns an array of B<Bond> objects corresponding to bonds from I<Atom> to any other
hydrogen atom in a molecule.

=item B<GetBondsToNonHydrogenAtoms>

    @Bonds = $Atom->GetBondsToNonHydrogenAtoms();

Returns an array of B<Bond> objects corresponding to bonds from I<Atom> to other bonded
non-hydrogen atoms in a molecule.

=item B<GetExactMass>

    $ExactMass = $Atom->GetExactMass();

Returns exact mass of an I<Atom> which correspond to one of these three values: explicity set
I<ExactMass> property; mass of natural isotope for an explicty set value of I<MassNumber>; most
abundant natural isotope mass for I<Atom> with valid atomic number value available by
B<PerodicTable> module.

=item B<GetExplicitHydrogens>

    $NumOfExplicitHydrogens = $Atom->GetExplicitHydrogens();

Returns number of hydrogens explicity bonded to an I<Atom> in a molecule.

=item B<GetFormalCharge>

    $FormalCharge = $Atom->GetFormalCharge();

Returns formal charge of an I<Atom> in a molecule.

=item B<GetFreeRadicalElectrons>

    $FreeRadicalElectrons = $Atom->GetFreeRadicalElectrons();

Returns number of free radical electrons corresponding to to one of these
three values: I<FreeRadicalElectrons> property; I<SpinMultiplicity> property; value
of 0.

For atoms with explicit assignment of I<SpinMultiplicity> atom property values,

    Singlet  - two unparied electrons corresponding to one spin state
    Doublet - free radical; an unpaired electron corresponding to two
              spin states
    Triplet - two unparied electrons corresponding to three spin states
              (divalent carbon atoms: carbenes)

B<FreeRadicalElectrons> are calculated as follows:

    Doublet: 1 (one valence electron not available for bonding)
    Singlet: 2 (two valence electrons not available for bonding)
    Triplet: 2 (two valence electrons not available for bonding)

=item B<GetGroupNumber>

    $GroupNumber = $Atom->GetGroupNumber();

Returns group number of an I<Atom> in a molecule with a valid atomic number.

=item B<GetHeavyAtomNeighbors>

    $NumOfHeavyAtoms = $Atom->GetHeavyAtomNeighbors();
    @HeavyAtoms = $Atom->GetHeavyAtomNeighbors();

Return number of heavy atoms or an array of B<Atom> objects corresponding to heavy atoms
bonded to an I<Atom> in a molecule.

=item B<GetHeavyAtomNeighborsAtomInformation>

    ($NumOfAtomNeighbors, $AtomNeighborsRef,
     $NumOfAtomNeighborsType, $AtomNeighborsTypeMapRef) = $Atom->
                              GetHeavyAtomNeighborsAtomInformation();

Returns atoms information for all non-hydrogen atoms attached to an I<Atom>
in a molecule.

The following values are returned:

    o Number of non-hydrogen atom neighbors
    o A reference to an array containing atom objects corresponding to
      non-hydrogen atom neighbors
    o Number of different types of non-hydrogen atom neighbors
    o A reference to a hash containing atom symbol as key with value
      corresponding to its count for non-hydrogen atom neighbors

=item B<GetHeavyAtomNeighborsBondformation>

    ($NumOfBonds, $BondTypeCountMapRef,
    $AtomsBondTypesCountMapRef,
    $AtomsBondTypeAtomsMap) = $Atom->
                              GetHeavyAtomNeighborsBondformation();

Returns bonds information for all non-hydrogen atoms attached to an I<Atom>
in a molecule.

The following values are returned:

    o Number of bonds to non-hydrogen atom neighbors
    o A reference to an array containing bond objects corresponding to
      non-hydrogen atom neighbors
    o A reference to a hash containing bond type as key with value
      corresponding to its count for non-hydrogen atom neighbors. Bond
      types are: Single, Double or Triple
    o A reference to a hash containing atom symbol as key pointing to bond
      type as second key with values corresponding to count of bond types for atom
      symbol for non-hydrogen atom neighbors
    o A reference to a hash containing atom symbol as key pointing to bond
      type as second key with values corresponding to atom objects array involved
      in corresponding bond type for atom symbol for non-hydrogen atom neighbors

=item B<GetHighestCommonValence>

    $HighestCommonValence = $Atom->GetHighestCommonValence();

Returns highest common valence of an I<Atom> which corresponds to either explicity set
I<HighestCommonValence> atom property or highest common valence of the corresponding
element in the periodic table available by B<PerodicTable> module.

=item B<GetHydrogens>

    $NumOfHydrogens = $Atom->GetHydrogens();

Returns total number of hydrogens for an I<Atom> in a molecule including both hydrogen atom
neighbors and implicit hydrogens.

=item B<GetHydrogenAtomNeighbors>

    $NumOfHydrogenAtomNeighbors = $Atom->GetHydrogenAtomNeighbors();
    @HydrogenAtomNeighbors = $Atom->GetHydrogenAtomNeighbors();

Return number of hydrogen atoms or an array of I<Atom> objects corresponding to hydrogen
atoms bonded to an I<Atom> in a molecule.

=item B<GetImplicitHydrogens>

    $NumOfImplicitHydrogens = $Atom->GetImplicitHydrogens();

Returns number of implicit hydrogens for an I<Atom> in a molecule. This value either
corresponds to explicitly set I<ImplicitHydrogens> atom property or calculated as the
difference between the value of potential total valence and sum of bond orders to
both hydrogen and non-hydrogen atom neighbors.

=item B<GetPotentialTotalCommonValence>

    $PotentialTotalValence = $Atom->GetPotentialTotalCommonValence();

Returns potential total common valence of an I<Atom> in a molecule corresponding
to a specific valence model set for the molecule using its B<SetValenceModel> method
or default internal valence model. It is used during the calculation of missing or
implicit hydrogens.

The current release of MayaChemTools supports three valence models: I<MDLValenceModel,
DaylightValenceModel, InternalValenceModel or MayaChemToolsValenceModel>.

For I<MDLValenceModel> and I<DaylightValenceModel>, the following data files, distributed
with the package, are used to calculate potential total valence:

    lib/data/MDLValenceModelData.csv
    lib/data/DaylightValenceModelData.csv

The calculation of potential total common valence for these two models is performed as
follows: Calculate current effective total valence of the I<Atom> by adding up the bond
order of its neighbors and number of free radical electrons; Find available common valence
for the I<Atom>, corresponding to any specified formal charge, higher than the effective
total valence, and return it as I<PotentialTotalValence>.

The calculation of potential total common valence For I<InternalValenceModel> or
I<MayaChenToolsValenceModel> doesn't uses B<PeriodicTable> module to retrieve values
for common valence, which in turn reads in PeriodicTableElements.csv file distributed with
the package.

For elements with one one common valence, potential total common valence corresponds
to:

    CommonValence + FormalCharge - FreeRadicalElectrons

For elements with multiple common valences, each common valence is used to
calculate total potential common valence as shown above, and the first total potential
common valence greater than the sum of bond orders to all neighbors is selected as
the final total common valence.

FormalCharge sign is reversed for electropositive elements with positive formal charge
during common valence calculations. Electropositive elements, metals and transition elements,
have usually plus formal charge and it leads to decrease in common valence; the negative
formal charge should result in the decrease of common valence.

The sign of formal charge is adjusted as follows.

Group numbers > 14 - Group numbers 15 (N), 16 (O), 17 (F), 18 (He):

Formal charge sign is not adjusted. Positive and negative values result in the
increase and decrease of valence.

Group 14 containing C, Si, Ge, Sn, Pb...:

Formal charge sign is reversed for positive values. Both positive and negative
values result in the decrease of valence.

Group 13 containing B, Al, Ga, In, Tl...:

Formal charge sign is always reversed. Positive and negative values result in the
decrease and increase of valence.

Groups 1 (H) through 12 (Zn)...:

Formal charge sign is reversed for positive values. Both positive and negative
values result in the decrease of valence.

Lanthanides and actinides:

Formal charge sign is reversed for positive values. Both positive and negative
values result in the decrease of valence.

=item B<GetLargestBondOrder>

    $LargestBO =$Atom->GetLargestBondOrder();

Returns largest bond order for an I<Atom> among the bonds to other atoms in a molecule.

=item B<GetLargestBondOrderToHeavyAtoms>

    $LargestBO =$Atom->GetLargestBondOrderToHeavyAtoms();

Returns largest bond order for an I<Atom> among the bonds to other heavy atoms in a molecule.

=item B<GetLargestBondOrderToNonHydrogenAtoms>

    $LargestBO =$Atom->GetLargestBondOrderToNonHydrogenAtoms();

Returns largest bond order for an I<Atom> among the bonds to other non-hydrogen atoms
in a molecule.

=item B<GetLargestRing>

    @RingAtoms = $Atom->GetLargestRing();

Returns an array of ring I<Atom> objects corresponding to the largest ring containing I<Atom>
in a molecule.

=item B<GetLowestCommonValence>

    $LowestCommonValence = $Atom->GetLowestCommonValence();

Returns lowest common valence of an I<Atom> which corresponds to either explicity set
I<LowestCommonValence> atom property or highest common valence of the corresponding
element in the periodic table available by B<PerodicTable> module.

=item B<GetMassNumber>

    $MassNumber = $Aom->GetMassNumber();

Returns atomic weight of an B<Atom> which corresponds to either explicity set I<MassNumber>
atom property or mass number of the most abundant natural isotope of the corresponding element
in the periodic table available by B<PeriodicTable> module.

=item B<GetMissingHydrogens>

    $NumOfMissingHydrogens = $Atom->GetMissingHydrogens();

Returns number of missing hydrogens for an I<Atom> in a molecule. This value either
corresponds to explicitly set I<ImplicitHydrogens> atom property or calculated as the
difference between the value of potential total valence and sum of bond orders to
both hydrogen and non-hydrogen atom neighbors.

=item B<GetNeighbors>

    $NumOfNeighbors = $Atom->GetNeighbors();
    @Neighbors = $Atom->GetNeighbors();

Returns number of neighbor atoms or an array of I<Atom> objects corresponding to all
atoms bonded to an I<Atom> in a molecule.

=item B<GetNeighborsUsingAtomSpecification>

    @AtomNeighbors = $Atom->GetNeighborsUsingAtomSpecification($AtomSpec);
    $NumOfNeighbors = $Atom->GetNeighborsUsingAtomSpecification($AtomSpec);

    @AtomNeighbors = $Atom->GetNeighborsUsingAtomSpecification($AtomSpec,
                     @ExcludeNeighbors);

Returns number of neighbor atoms or an array of I<Atom> objects matching atom
specification corresponding to atom neighbors of an I<Atom> in a molecule. Optionally,
I<Atom> neighbors can be excluded from the neighbors list using I<ExcludeNeighbors>.

Notes:

    o AtomSpecification correspond to any valid AtomicInvariant based atomic specifications
      as supported by DoesAtomNeighborhoodMatch method
    o Multiple atom specifications can be used in a string delimited by comma

=item B<GetNonHydrogenAtomNeighbors>

    $NumOfNeighbors = $Atom->GetNonHydrogenAtomNeighbors();
    @Neighbors = $Atom->GetNonHydrogenAtomNeighbors();

Returns number of non-hydrogen atoms or an array of B<Atom> objects corresponding to non-hydrogen
atoms bonded to an I<Atom> in a molecule.

=item B<GetNonHydrogenAtomNeighborsAtomInformation>

    ($NumOfAtomNeighbors, $AtomNeighborsRef,
     $NumOfAtomNeighborsType, $AtomNeighborsTypeMapRef) = $Atom->
                              GetNonHydrogenAtomNeighborsAtomInformation();

Returns atoms information for all non-hydrogen atoms attached to an I<Atom>
in a molecule.

The following values are returned:

    o Number of non-hydrogen atom neighbors
    o A reference to an array containing atom objects corresponding to
      non-hydrogen atom neighbors
    o Number of different types of non-hydrogen atom neighbors
    o A reference to a hash containing atom symbol as key with value
      corresponding to its count for non-hydrogen atom neighbors


=item B<GetNonHydrogenAtomNeighborsBondInformation>

    ($NumOfBonds, $BondTypeCountMapRef,
    $AtomsBondTypesCountMapRef,
    $AtomsBondTypeAtomsMap) = $Atom->
                              GetNonHydrogenAtomNeighborsBondInformation();

Returns bonds information for all non-hydrogen atoms attached to an I<Atom>
in a molecule.

The following values are returned:

    o Number of bonds to non-hydrogen atom neighbors
    o A reference to an array containing bond objects corresponding to
      non-hydrogen atom neighbors
    o A reference to a hash containing bond type as key with value
      corresponding to its count for non-hydrogen atom neighbors. Bond
      types are: Single, Double or Triple
    o A reference to a hash containing atom symbol as key pointing to bond
      type as second key with values corresponding to count of bond types for atom
      symbol for non-hydrogen atom neighbors
    o A reference to a hash containing atom symbol as key pointing to bond
      type as second key with values corresponding to atom objects array involved
      in corresponding bond type for atom symbol for non-hydrogen atom neighbors

=item B<GetNonHydrogenNeighborOfHydrogenAtom>

    $Atom = $Atom->GetNonHydrogenNeighborOfHydrogenAtom();

Returns non-hydrogen or heavy atom neighbor of a hydrogen atom in a molecule..

=item B<GetNumOfAromaticBondsToHeavyAtoms>

    $NumOfBonds = $Atom->GetNumOfAromaticBondsToHeavyAtoms();

Returns number of aromatic bonds from an I<Atom> to other non-hydrogen or heavy atoms in
a molecule.

=item B<GetNumOfAromaticBondsToNonHydrogenAtoms>

    $NumOfBonds = $Atom->GetNumOfAromaticBondsToNonHydrogenAtoms();

Returns number of aromatic bonds from an I<Atom> to other non-hydrogen or heavy atoms in
a molecule.

=item B<GetNumOfBonds>

    $NumOfBonds = $Atom->GetNumOfBonds();

Returns number of bonds from an I<Atom> to other atoms in a molecule.

=item B<GetNumOfBondsAvailableForHeavyAtoms>

    $NumOfBonds = $Atom->GetNumOfBondsAvailableForHeavyAtoms();

Get number of bonds available to form additional bonds with heavy atoms, excluding
any implicit bonds to hydrogens set using I<ImplicitHydrogens> property.

It's different from number of implicit or missing hydrogens, both of which are equivalent.

For example, in a SMILES string, [nH] ring atom corresponds to an aromatic nitrogen.
Although the hydrogen specified for n is treated internally as implicit hydrogen and shows
up in missing hydrogen count, it's not available to participate in double bonds to additional
heavy atoms.

=item B<GetNumOfBondsAvailableForNonHydrogenAtoms>

    $NumOfBonds = $Atom->GetNumOfBondsAvailableForNonHydrogenAtoms();

Get number of bonds available to form additional bonds with heavy atoms, excluding
any implicit bonds to hydrogens set using ImplicitHydrogens property.

=item B<GetNumOfBondsToHeavyAtoms>

    $NumOfBondsToHeavyAtoms = $Atom->GetNumOfBondsToHeavyAtoms();

Returns number of bonds from an I<Atom> to other heavy atoms in a molecule.

=item B<GetNumOfBondsToHydrogenAtoms>

    $NumOfBonds = $Atom->GetNumOfBondsToHydrogenAtoms();

Returns number of bonds from an I<Atom> to other hydrogen atoms in a molecule.

=item B<GetNumOfBondsToNonHydrogenAtoms>

    $NumOfBonds = $Atom->GetNumOfBondsToNonHydrogenAtoms();

Returns number of bonds from an I<Atom> to other non-hydrogen atoms in a molecule.

=item B<GetNumOfBondTypesToHeavyAtoms>

    ($NumOfSingleBonds, $NumOfDoubleBonds,
     $NumOfTripleBonds, $NumOfAromaticBonds) = $Atom->
                   GetNumOfBondTypesToHeavyAtoms($CountAromaticBonds);

Get number of single, double, triple, and aromatic bonds from an I<Atom> to all other
non-hydrogen atoms in a molecule.

Value of I<CountAtomaticBonds> parameter controls whether number of aromatic
bonds is returned; default is not to count aromatic bonds. During  counting of
aromatic bonds, the bond marked aromatic is not included in the count
of other bond types.

=item B<GetNumOfBondTypesToNonHydrogenAtoms>

    ($NumOfSingleBonds, $NumOfDoubleBonds,
     $NumOfTripleBonds, $NumOfAromaticBonds) = $Atom->
             GetNumOfBondTypesToNonHydrogenAtoms($CountAromaticBonds);

Get number of single, double, triple, and aromatic bonds from an I<Atom> to all other
non-hydrogen atoms in a molecule.

Value of I<CountAtomaticBonds> parameter controls whether number of aromatic
bonds is returned; default is not to count aromatic bonds. During  counting of
aromatic bonds, the bond marked aromatic is not included in the count
of other bond types.

=item B<GetNumOfDoubleBondsToHeavyAtoms>

    $NumOfDoubleBonds = $Atom->GetNumOfDoubleBondsToHeavyAtoms();

Returns number of double bonds from an I<Atom> to other heavy atoms or non-hydrogen
atoms in a molecule.

=item B<GetNumOfDoubleBondsToNonHydrogenAtoms>

    $NumOfDoubleBonds =$Atom->GetNumOfDoubleBondsToNonHydrogenAtoms();

Returns number of double bonds from an I<Atom> to other heavy atoms or non-hydrogen
atoms in a molecule.

=item B<GetNumOfHeavyAtomNeighbors>

    $NumOfNeighbors = $Atom->GetNumOfHeavyAtomNeighbors();

Returns number heavy atom neighbors for an I<Atom> in a molecule.

=item B<GetNumOfHydrogenAtomNeighbors>

    $NumOfNeighbors = $Atom->GetNumOfHydrogenAtomNeighbors();

Returns number hydrogens atom neighbors for an I<Atom> in a molecule.

=item B<GetNumOfMissingHydrogens>

    $NumOfMissingHydrogens = $Atom->GetNumOfMissingHydrogens();

Returns number of implicit hydrogens for an I<Atom> in a molecule. This value either
corresponds to explicitly set I<ImplicitHydrogens> atom property or calculated as the
difference between the value of potential total valence and sum of bond orders to
both hydrogen and non-hydrogen atom neighbors.

=item B<GetNumOfExplicitHydrogens>

    $NumOfExplicitHydrogens = $Atom->GetNumOfExplicitHydrogens();

Returns number hydrogens atom neighbors for an I<Atom> in a molecule.

=item B<GetNumOfHydrogens>

    $NumOfHydrogens = $Atom->GetNumOfHydrogens();

Returns total number of hydrogens for an I<Atom> in a molecule including both hydrogen atom
neighbors and implicit hydrogens.

=item B<GetNumOfImplicitHydrogens>

    $NumOfImplicitHydrogens = $Atom->GetNumOfImplicitHydrogens();

Returns number of implicit hydrogens for an I<Atom> in a molecule. This value either
corresponds to explicitly set I<ImplicitHydrogens> atom property or calculated as the
difference between the value of potential total valence and sum of bond orders to
both hydrogen and non-hydrogen atom neighbors.

=item B<GetNumOfNeighbors>

    $NumOfNeighbors = $Atom->GetNumOfNeighbors();

Returns number atom neighbors for an I<Atom> in a molecule.

=item B<GetNumOfNonHydrogenAtomNeighbors>

    $NumNeighbors = $This->GetNumOfNonHydrogenAtomNeighbors();

Returns number non-hydrogens atom neighbors for an I<Atom> in a molecule.

=item B<GetNumOfRings>

    $NumOfRings = $Atom->GetNumOfRings();

Returns number of rings containing I<Atom> in a molecule.

=item B<GetNumOfRingsWithEvenSize>

    $NumOfRings = $Atom->GetNumOfRingsWithEvenSize();

Returns number of rings with even size containing I<Atom> in a molecule.

=item B<GetNumOfRingsWithOddSize>

    $NumOfRings = $Atom->GetNumOfRingsWithOddSize();

Returns number of rings with odd size containing I<Atom> in a molecule.

=item B<GetNumOfRingsWithSize>

    $NumOfRings = $Atom->GetNumOfRingsWithSize($RingSize);

Returns number of rings with specific I<RingSize> containing I<Atom> in a molecule.

=item B<GetNumOfRingsWithSizeGreaterThan>

    $NumOfRings = $Atom->GetNumOfRingsWithSizeGreaterThan($RingSize);

Returns number of rings with size greater than specific I<RingSize> containing I<Atom>
in a molecule.

=item B<GetNumOfRingsWithSizeLessThan>

    $NumOfRings = $Atom->GetNumOfRingsWithSizeLessThan($RingSize);

Returns number of rings with size less than specific I<RingSize> containing I<Atom> in a molecule.

=item B<GetNumOfSigmaAndPiBondsToHeavyAtoms>

    ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->
                              GetNumOfSigmaAndPiBondsToHeavyAtoms();

Get number of sigma and pi bonds from an I<Atom> to all other non-hydrogen
atoms in a molecule.

Sigma and pi bonds are counted using the following methodology: a single bond
correspond to one sigma bond; a double bond contributes one to sigma bond count
and one to pi bond count; a triple bond contributes one to sigma bond count and
two to pi bond count.


=item B<GetNumOfSigmaAndPiBondsToNonHydrogenAtoms>

    ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->
                              GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();

Get number of sigma and pi bonds from an I<Atom> to all other non-hydrogen
atoms in a molecule.

Sigma and pi bonds are counted using the following methodology: a single bond
correspond to one sigma bond; a double bond contributes one to sigma bond count
and one to pi bond count; a triple bond contributes one to sigma bond count and
two to pi bond count.

=item B<GetNumOfSingleBondsToNonHydrogenAtoms>

    $NumOfSingleBonds =$Atom->GetNumOfSingleBondsToNonHydrogenAtoms();

Returns number of single bonds from an I<Atom> to other heavy atoms or non-hydrogen
atoms in a molecule.

=item B<GetNumOfSingleBondsToHeavyAtoms>

    $NumOfSingleBonds = $Atom->GetNumOfSingleBondsToHeavyAtoms();

Returns number of single bonds from an I<Atom> to other heavy atoms or non-hydrogen
atoms in a molecule.

=item B<GetNumOfTripleBondsToNonHydrogenAtoms>

    $NumOfTripleBonds =$Atom->GetNumOfTripleBondsToNonHydrogenAtoms();

Returns number of triple bonds from an I<Atom> to other heavy atoms or non-hydrogen
atoms in a molecule.

=item B<GetNumOfTripleBondsToHeavyAtoms>

    $NumOfTripleBonds = $Atom->GetNumOfTripleBondsToHeavyAtoms();

Returns number of triple bonds from an I<Atom> to other heavy atoms or non-hydrogen
atoms in a molecule.

=item B<GetPeriodNumber>

    $PeriodNumber = $Atom->GetPeriodNumber();

Returns periodic table period number for an I<Atom> in a molecule with a valid atomic number .

=item B<GetRings>

    @Rings = $Aotm->GetRings();

Returns an array of references to arrays containing ring atoms corressponding to all rings containing
I<Atom> in a molecule.

=item B<GetRingsWithEvenSize>

    @Rings = $Aotm->GetRingsWithEvenSize();

Returns an array of references to arrays containing ring atoms corressponding to all rings with even size
containing I<Atom> in a molecule.

=item B<GetRingsWithOddSize>

    @Rings = $Aotm->GetRingsWithOddSize();

Returns an array of references to arrays containing ring atoms corressponding to all rings with odd size
containing I<Atom> in a molecule.

=item B<GetRingsWithSize>

    @Rings = $Aotm->GetRingsWithSize($RingSize);

Returns an array of references to arrays containing ring atoms corressponding to all rings with specific
I<RingSize >containing I<Atom> in a molecule.

=item B<GetRingsWithSizeGreaterThan>

    @Rings = $Aotm->GetRingsWithSizeGreaterThan($RingSize);

Returns an array of references to arrays containing ring atoms corressponding to all rings with size
greater than specific I<RingSize >containing I<Atom> in a molecule.

=item B<GetRingsWithSizeLessThan>

    @Rings = $Aotm->GetRingsWithSizeLessThan($RingSize);

Returns an array of references to arrays containing ring atoms corressponding to all rings with size
less than specific I<RingSize >containing I<Atom> in a molecule.

=item B<GetSizeOfLargestRing>

    $Size = $Atom->GetSizeOfLargestRing();

Returns size of the largest ring containing I<Atom> in a molecule.

=item B<GetSizeOfSmallestRing>

    $Size = $Atom->GetSizeOfSmallestRing();

Returns size of the smallest ring containing I<Atom> in a molecule.

=item B<GetSmallestRing>

    @RingAtoms = $Atom->GetSmallestRing();

Returns an array of ring I<Atom> objects corresponding to the largest ring containing I<Atom>
in a molecule.

=item B<GetSpinMultiplicity>

    $SpinMultiplicity = $Atom->GetSpinMultiplicity();

Returns spin multiplicity of an I<Atom> corresponding to one of these three
values: explicitly set B<SpinMultiplicity> property value; calculated from
B<FreeRadicalElectrons> property; value of 0.

The B<SpinMultiplicity> is calculate from I<FreeRadicalElectrons> property as
follows:

    FreeRadicalElectrons: 1; SpinMultiplicity: 2
    FreeRadicalElectrons: 2; SpinMultiplicity: 1
    FreeRadicalElectrons: other; SpinMultiplicity: 0

=item B<GetSumOfBondOrders>

    $SumBondOrders = $Atom->GetSumOfBondOrders();

Returns sum of bond orders corresponding to all atoms bonded to an I<Atom> in a molecule.

=item B<GetSumOfBondOrdersToHeavyAtoms>

    $SumBondOrders = $Atom->GetSumOfBondOrdersToHeavyAtoms();

Returns sum of bond orders corresponding to all heavy atoms bonded to an I<Atom> in a molecule.

=item B<GetSumOfBondOrdersToHydrogenAtoms>

    $SumBondOrders = $Atom->GetSumOfBondOrdersToHydrogenAtoms();

Returns sum of bond orders corresponding to all hydrogen atoms bonded to an I<Atom> in a molecule.

=item B<GetSumOfBondOrdersToNonHydrogenAtoms>

    $SumBondOrders = $Atom->GetSumOfBondOrdersToNonHydrogenAtoms();

Returns sum of bond orders corresponding to all non-hydrogen atoms bonded to an I<Atom>
in a molecule.

=item B<GetValence>

    $Valence = $Atom->GetValence();

Returns valence of an I<Atom> in a molecule. Valence corresponds to number of electrons used
by an atom in bonding:

    Valence = ValenceElectrons - ValenceFreeElectrons = BondingElectrons

Single, double and triple bonds with bond orders of 1, 2, and 3 correspond to
contribution of 1, 2, and 3 bonding electrons. So:

    Valence = SumOfBondOrders + NumOfMissingHydrogens + FormalCharge

where positive and negative values of FormalCharge increase and decrease the number of bonding
electrons, respectively.

The current release of MayaChemTools supports the following three valence models, which
are used during calculation of implicit hydrogens: MDLValenceModel, DaylightValenceModel,
InternalValenceModel or MayaChemToolsValenceModel.

Notes:

    . Missing hydrogens are included in the valence.
    . For neutral molecules, valence and sum of bond orders are equal.
    . For molecules containing only single bonds, SumOfBondOrders and
      NumOfBonds are equal.
    . Free radical electrons lead to the decrease in valence. For atoms with
      explicit assignment of SpinMultiplicity property values corresponding to
      Singlet (two unparied electrons corresponding to one spin state), Doublet
      (free radical; an unpaired electron corresponding to two spin states),
      and Triplet (two unparied electrons corresponding to three spin states;
      divalent carbon atoms (carbenes)), FreeRadicalElectrons are calculated as follows:

       SpinMultiplicity: Doublet(2); FreeRadicalElectrons: 1 (one valence
           electron not available for bonding)
       SpinMultiplicity: Singlet(1)/Triplet(3); FreeRadicalElectrons: 2 (two
           valence electrons not available for bonding)

=item B<GetValenceElectrons>

    $ValenceElectrons = $Atom->GetValenceElectrons();

Returns valence electrons for an B<Atom> which corresponds to either explicity set I<ValenceElectrons>
atom property or valence electrons for the corresponding element in the periodic table available by
B<PeriodicTable> module.

=item B<GetValenceFreeElectrons>

    $ValenceFreeElectrons = $Atom->GetValenceFreeElectrons();
    $ValenceFreeElectrons = $Atom->GetValenceFreeElectrons(
                            $ExcludeFreeRadicalElectrons);

Returns valence frees electrons for an B<Atom> in a molecule. It corresponds to:

    ValenceElectrons - Valence
    or
    ValenceElectrons - NumOfMissingHydrogens - SumOfBondOrders - FormalCharge

Free radical electrons are included in the valence free electrons count by default.

Examples:

    NH3: ValenceFreeElectrons = 5 - 3 = 5 - 3 - 0 - 0 = 2
    NH2: ValenceFreeElectrons = 5 - 3 = 5 - 2 - 1 - 0 = 2
    NH4+; ValenceFreeElectrons = 5 - 5 = 5 - 4 - 0 - 1 = 0
    NH3+; ValenceFreeElectrons = 5 - 5 = 5 - 3 - 1 - 1 = 0
    C(=O)O- : ValenceFreeElectrons on O- = 6 - 0 = 6 - 1 - 0 - (-1) = 6
    C(=O)O- : ValenceFreeElectrons on =O = 6 - 2 = 6 - 2 - 0 - 0 = 4

=item B<GetX>

    $X = $Atom->GetX();

Returns value of X-coordinate for an I<Atom>.

=item B<GetXYZ>

    @XYZ = $Atom->GetXYZ();
    $XYZRef = $Atom->GetXYZ();

Returns an array or a reference to an array containing values for I<Atom> coordinates.

=item B<GetXYZVector>

    $XYZVector = $Atom->GetXYZVector();

Returns a I<Vector> object containing values for I<Atom> coordinates

=item B<GetY>

    $Y = $Atom->GetY();

Returns value of Y-coordinate for an I<Atom>.

=item B<GetZ>

    $Z = $Atom->GetZ();

Returns value of Z-coordinate for an I<Atom>.

=item B<IsAmideCarbon>

    $Status = $Atom->IsAmideCarbon();

Returns 1 or 0 based on whether it's amide carbon I<Atom>.

An amide group is defineds as:

    R-C(=O)-N(-R')-R''

where:

    o R = Hydrogen or groups of atoms attached through carbon
    o R' = Hydrogens or groups of atoms attached through carbon or
           hetro atoms
    o R'' = Hydrogens or groups of atoms attached through carbon or
           hetro atoms

=item B<IsAmideNitrogen>

    $Status = $Atom->IsAmideNitrogen();

Returns 1 or 0 based on whether it's amide nitrogen I<Atom>.

=item B<IsAromatic>

    $Status = $Atom->IsAromatic();

Returns 1 or 0 based on whether it's an aromatic I<Atom>.

=item B<IsArsenic>

    $Status = $Atom->IsArsenic();

Returns 1 or 0 based on whether it's an arsenic I<Atom>.

=item B<IsBondedToAtom>

    $Status = $Atom->IsBondedToAtom($OtherAtom);

Returns 1 or 0 based on whether I<Atom> is bonded to I<OtherAtom>.

=item B<IsBromine>

    $Status = $Atom->IsBromine();

Returns 1 or 0 based on whether it's a bromine I<Atom>.

=item B<IsCarbon>

    $Status = $Atom->IsCarbon();

Returns 1 or 0 based on whether it's a carbon I<Atom>.

=item B<IsCarboxylCarbon>

    $Status = $Atom->IsCarboxylCarbon();

Returns 1 or 0 based on whether it's a carboxyl carbon atom in carboxyl group:
R-C(=O)-OH.

=item B<IsCarboxylOxygen>

    $Status = $Atom->IsCarboxylOxygen();

Returns 1 or 0 based on whether it's a carboxyl oxygen atom in carboxyl group:
R-C(=O)-OH.

=item B<IsCarboxylateCarbon>

    $Status = $Atom->IsCarboxylateCarbon();

Returns 1 or 0 based on whether it's a carboxylate carbon atom in carboxyl group:
R-C(=O)-O-.

=item B<IsCarboxylateOxygen>

    $Status = $Atom->IsCarboxylateOxygen();

Returns 1 or 0 based on whether it's a carboxylate oxygen atom in carboxyl group:
R-C(=O)-O-.

=item B<IsChlorine>

    $Status = $Atom->IsChlorine();

Returns 1 or 0 based on whether it's a chlorine I<Atom>.

=item B<IsFluorine>

    $Status = $Atom->IsFluorine();

Returns 1 or 0 based on whether it's a fluorine I<Atom>.

=item B<IsFunctionalClassType>

    $Status =$Atom->IsFunctionalClassType($Type);

Returns 1 or 0 based on whether it's a specified functional class I<Type>.

The current release of MayaChemTools supports following abbreviations and descriptive
names for I<FunctionalClassType>:

    HBD: HydrogenBondDonor
    HBA: HydrogenBondAcceptor
    PI :  PositivelyIonizable
    NI : NegativelyIonizable
    Ar : Aromatic
    Hal : Halogen
    H : Hydrophobic
    RA : RingAtom
    CA : ChainAtom

The following definitions are used to determine functional class types: [ Ref 60-61, Ref 65-66 ]:

    HydrogenBondDonor: NH, NH2, OH
    HydrogenBondAcceptor: N[!H], O
    PositivelyIonizable: +, NH2
    NegativelyIonizable: -, C(=O)OH, S(=O)OH, P(=O)OH

=item B<IsGuadiniumCarbon>

    $Status = $Atom->IsGuadiniumCarbon();

Returns 1 or 0 based on whether it's a guadinium carbon in guadinium group by
checking its neighbors for a nitrogen in guadinium group.

=item B<IsGuadiniumNitrogen>

    $Status = $Atom->IsGuadiniumNitrogen();

Returns 1 or 0 based on whether it's a guadinium nitrogen in guadinium group.

A guadinium group is defined as:

    R2N-C(=NR)-(NR2) or R2N-C(=NR2+)-(NR2)

where:

    o R = Hydrogens or group of atoms attached through carbon
    o Only one of the three nitrogens has a double bond to carbon
      and has optional formal charge allowing it to be neutral or charged state

=item B<IsHBondAcceptor>

    $Status =$Atom->IsHBondAcceptor();
    $Status =$Atom->IsHBondAcceptor($HydrogenBondsType);

Returns 1 or 0 based on whether it's a hydrogen bond acceptor I<Atom>.

=item B<IsHBondDonor>

    $Status =$Atom->IsHBondDonor();
    $Status =$Atom->IsHBondDonor($HydrogenBondsType);

Returns 1 or 0 based on whether it's a hydrogen bond donor I<Atom>.

=item B<IsHydrogenBondAcceptor>

    $Status =$Atom->IsHydrogenBondAcceptor();
    $Status =$Atom->IsHydrogenBondAcceptor($HydrogenBondsType);

Returns 1 or 0 based on whether it's a hydrogen bond acceptor I<Atom>.

=item B<IsHydrogenBondDonor>

    $Status =$Atom->IsHydrogenBondDonor();
    $Status =$Atom->IsHydrogenBondDonor($HydrogenBondsType);

Returns 1 or 0 based on whether it's a hydrogen bond donor I<Atom>.

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

By default, I<HydrogenBondsType1> is used to calculate number hydrogen bond donor
and acceptor atoms. I<HydrogenBondsType2> corresponds to B<RuleOf5> definition
of hydrogen bond donors and acceptors.

=item B<IsHalogen>

    $Status =$Atom->IsHalogen();

Returns 1 or 0 based on whether it's a halogen I<Atom>.

=item B<IsHeteroAtom>

    $Status = $Atom->IsHeteroAtom();

Returns 0 or 1 based on whether it's a hetro I<Atom>. Following atoms are considered hetro atoms:
B<N, O, F, P, S, Cl, Br, I>.

=item B<IsHydrogen>

    $Status = $Atom->IsHydrogen();

Returns 1 or 0 based on whether it's a hydrogen I<Atom>.

=item B<IsHydrophobic>

    $Status =$Atom->IsHydrophobic();

Returns 1 or 0 based on whether it's a hydrophobic I<Atom>.

=item B<IsInRing>

    $Status = $Atom->IsInRing();

Returns 1 or 0 based on whether I<Atom> is present in a ring.

=item B<IsInRingOfSize>

    $Status = $Atom->IsInRingOfSize($Size);

Returns 1 or 0 based on whether I<Atom> is present in a ring of specific I<Size>.

=item B<IsIodine>

    $Status = $Atom->IsIodine();

Returns 1 or 0 based on whether it's an iodine I<Atom>.

=item B<IsIsotope>

    $Status =$Atom->IsIsotope();

Returns 1 or 0 based on whether it's an isotope I<Atom>.

=item B<IsLipophilic>

    $Status =$Atom->IsLipophilic();

Returns 1 or 0 based on whether it's a lipophilic I<Atom>.

=item B<IsMetallic>

    $Status = $Atom->IsMetallic();

Returns 1 or 0 based on whether it's a metallic I<Atom>.

=item B<IsNegativelyIonizable>

    $Status =$Atom->IsNegativelyIonizable();

Returns 1 or 0 based on whether it's a negatively ionizable atom I<Atom>.

=item B<IsNitrogen>

    $Status = $Atom->IsNitrogen();

Returns 1 or 0 based on whether it's a nitrogen I<Atom>.

=item B<IsNonCarbonOrHydrogen>

    $Status =$Atom->IsNonCarbonOrHydrogen();

Returns 1 or 0 based on whether it's not a carbon or hydrogen I<Atom>.

=item B<IsNotInRing>

    $Status = $Atom->IsNotInRing();

Returns 1 or 0 based on whether I<Atom> is not present in a ring.

=item B<IsOnlyInOneRing>

    $Status = $Atom->IsOnlyInOneRing();

Returns 1 or 0 based on whether I<Atom> is only present in one ring.

=item B<IsOxygen>

    $Status = $Atom->IsOxygen();

Returns 0 or 1 based on whether it's an oxygen I<Atom>.

=item B<IsPhosphorus>

    $Status = $Atom->IsPhosphorus();

Returns 0 or 1 based on whether it's a phosphorus I<Atom>.

=item B<IsPhosphateOxygen>

    $Status = $Atom->IsPhosphateOxygen();

Returns 1 or 0 based on whether it's a phosphate oxygen in phosphate group.

A phosphate group is defined as:

    AO-(O=)P(-OA)-OA

Where:

   A - Any group of atoms including hydrogens

=item B<IsPhosphatePhosphorus>

    $Status = $Atom->IsPhosphatePhosphorus();

Returns 1 or 0 based on whether it's a phosphate phosphorus in phosphate group.

=item B<IsPolarAtom>

    $Status = $Atom->IsPolarAtom();

Returns 0 or 1 based on whether it's a polar I<Atom>. Following atoms are considered polar atoms:
B<N, O, P, S>.

=item B<IsPolarHydrogen>

    $Status = $Atom->IsPolarHydrogen();

Returns 0 or 1 based on whether it's a hydrogen I<Atom> bonded to a polar atom.

=item B<IsPositivelyIonizable>

    $Status =$Atom->IsPositivelyIonizable();

Returns 1 or 0 based on whether it's a positively ionizable I<Atom>.

=item B<IsSaturated>

    $Status = $Atom->IsSaturated();

Returns 1 or 0 based on whether it's a saturated I<Atom>. An atom attached
to other atoms with only single bonds is considered a saturated atom.

=item B<IsSelenium>

    $Status = $Atom->IsSelenium();

Returns 0 or 1 based on whether it's a selenium I<Atom>.

=item B<IsStereoCenter>

    $Status = $Atom->IsStereoCenter();

Returns 0 or 1 based on whether it's marked as a stero center I<Atom> by explicit setting
of I<StereoCenter> atom propert to value of I<1>.

=item B<IsSilicon>

    $Status = $Atom->IsSilicon();

Returns 0 or 1 based on whether it's a silicon I<Atom>.

=item B<IsSulfur>

    $Status = $Atom->IsSulfur();

Returns 0 or 1 based on whether it's a sulfur I<Atom>.

=item B<IsSulphur>

    $Status = $Atom->IsSulphur();

Returns 0 or 1 based on whether it's a sulfur I<Atom>.

=item B<IsTellurium>

    $Status = $Atom->IsTellurium();

Returns 0 or 1 based on whether it's a tellurium I<Atom>.

=item B<IsTerminal>

    $Status = $Atom->IsTerminal();

Returns 0 or 1 based on whether it's a terminal I<Atom> attached to no
more than one non-hydrogen atom.

=item B<IsUnsaturated>

    $Status = $Atom->IsUnsaturated();

Returns 1 or 0 based on whether it's as unsaturated I<Atom>. An atom attached
to other atoms with at least one non-single bond is considered an unsaturated atom.

=item B<IsTopologicalPharmacophoreType>

    $Status =$Atom->IsTopologicalPharmacophoreType();

Returns 1 or 0 based on whether it's any of the supportyed topological pharmacophore
I<Atom> type. See I<IsFunctionalClassType> for a list of supported types.

=item B<SetAtomSymbol>

    $Atom->SetAtomSymbol($AtomicSymbol);

Sets atom symbol for I<Atom> and returns I<Atom> object. The appropriate atomic number is also
set automatically.

=item B<SetAtomicNumber>

    $Atom->SetAtomicNumber($AtomicNumber);

Sets atomic number for I<Atom> and returns I<Atom> object. The appropriate atom symbol is also
set automatically.

=item B<SetMassNumber>

    $Atom->SetMassNumber($MassNumber);

Sets mass number for I<Atom> and returns I<Atom> object.

=item B<SetStereoCenter>

    $Atom->SetStereoCenter($StereoCenter);

Sets stereo center for I<Atom> and returns I<Atom> object.

=item B<SetStereochemistry>

    $Atom->SetStereochemistry($Stereochemistry);

Sets stereo chemistry for I<Atom> and returns I<Atom> object.

=item B<SetX>

    $Atom->SetX($Value);

Sets X-coordinate value for I<Atom> and returns I<Atom> object.

=item B<SetXYZ>

    $Atom->SetXYZ(@XYZValues);
    $Atom->SetXYZ($XYZValuesRef);
    $Atom->SetXYZ($XYZVector);

Sets I<Atom> coordinates using an array, reference to an array or a I<Vector> object and
returns I<Atom> object.

=item B<SetY>

    $Atom->SetY($Value);

Sets Y-coordinate value for I<Atom> and returns I<Atom> object.

=item B<SetZ>

    $Atom->SetZ($Value);

Sets Z-coordinate value for I<Atom> and returns I<Atom> object.

=item B<StringifyAtom>

    $AtomString = $Atom->StringifyAtom();

Returns a string containing information about I<Atom> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Bond.pm, Molecule.pm, MoleculeFileIO.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
