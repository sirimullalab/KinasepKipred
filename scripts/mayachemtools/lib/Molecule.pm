package Molecule;
#
# File: Molecule.pm
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
use MathUtil;
use PeriodicTable;
use Text::ParseWords;
use TextUtil;
use FileUtil;
use Graph;
use Atom;
use Bond;
use MolecularFormula;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Graph ObjectProperty Exporter);
@EXPORT = qw(IsMolecule);
@EXPORT_OK = qw(FormatElementalCompositionInformation GetSupportedAromaticityModels IsSupportedAromaticityModel);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, $ObjectID, %AromaticityModelsDataMap, %CanonicalAromaticityModelNamesMap);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyMolecule';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeMolecule();

  if (keys %NamesAndValues) { $This->_InitializeMoleculeProperties(%NamesAndValues); }

  return $This;
}

# Initialize object data...
sub _InitializeMolecule {
  my($This) = @_;
  my($ObjectID) = _GetNewObjectID();

  # All other property names and values along with all Set/Get<PropertyName> methods
  # are implemented on-demand using ObjectProperty class.
  $This->{ID} = $ObjectID;
  $This->{Name} = "Molecule ${ObjectID}";
}

# Initialize molecule properties...
sub _InitializeMoleculeProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # ID to keep track of objects...
  $ObjectID = 0;

  # Load molecule class data...
  _LoadMoleculeClassData();
}

# Setup an explicit SetID method to block setting of ID by AUTOLOAD function...
sub SetID {
  my($This, $Value) = @_;

  carp "Warning: ${ClassName}->SetID: Object ID can't be changed: it's used for internal tracking...";

  return $This;
}

# Add an atom...
sub AddAtom {
  my($This, $Atom) = @_;

  if (!defined $Atom) {
    carp "Warning: ${ClassName}->AddAtom: No atom added: Atom must be specified...";
    return undef;
  }
  if ($This->HasAtom($Atom)) {
    carp "Warning: ${ClassName}->AddAtom: No atom added: Atom already exists...";
    return undef;
  }
  return $This->_AddAtom($Atom);
}

# Add an atom...
sub _AddAtom {
  my($This, $Atom) = @_;

  # Assign atom to this molecule...
  $Atom->_SetMolecule($This);

  # Add it to the graph as a vertex...
  my($AtomID);
  $AtomID = $Atom->GetID();
  $This->AddVertex($AtomID);
  $This->SetVertexProperty('Atom', $Atom, $AtomID);

  return $This;
}

# Add atoms...
sub AddAtoms {
  my($This, @Atoms) = @_;

  if (!@Atoms) {
    carp "Warning: ${ClassName}->AddAtoms: No atoms added: Atoms list must be specified...";
    return undef;
  }
  my($Atom);
  for $Atom (@Atoms) {
    $This->AddAtom($Atom);
  }
  return $This;
}

# Create an atom and add it to molecule...
sub NewAtom {
  my($This, %NamesAndValues) = @_;
  my($Atom);

  $Atom = new Atom(%NamesAndValues);
  $This->AddAtom($Atom);

  return $Atom;
}

# Delete an atom...
sub DeleteAtom {
  my($This, $Atom) = @_;

  if (!defined $Atom) {
    carp "Warning: ${ClassName}->DeleteAtom: No atom deleted: Atom must be specified...";
    return undef;
  }
  # Does the atom exist in  molecule?
  if (!$This->HasAtom($Atom)) {
    carp "Warning: ${ClassName}->DeleteAtom: No atom deleted: Atom doesn't exist...";
    return undef;
  }
  return $This->_DeleteAtom($Atom);
}

# Delete atom...
sub _DeleteAtom {
  my($This, $Atom) = @_;

  my($AtomID);
  $AtomID = $Atom->GetID();
  $This->DeleteVertex($AtomID);

  return $This;
}

# Delete atoms...
sub DeleteAtoms {
  my($This, @Atoms) = @_;

  if (!@Atoms) {
    carp "Warning: ${ClassName}->DeleteAtoms: No atoms added: Atoms list must be specified...";
    return undef;
  }
  my($Atom);
  for $Atom (@Atoms) {
    $This->DeleteAtom($Atom);
  }

  return $This;
}

# Is this atom present?
sub HasAtom {
  my($This, $Atom) = @_;

  if (!defined $Atom) {
    return 0;
  }
  if (!$Atom->HasProperty('Molecule')) {
    # It's not in molecule...
    return 0;
  }
  my($AtomID);
  $AtomID = $Atom->GetID();
  if (!$This->HasVertex($AtomID)) {
    # It's not in molecule...
    return 0;
  }
  my($Molecule);
  $Molecule = $Atom->GetProperty('Molecule');

  return ($This->HasVertex($AtomID) && $This->GetID() == $Molecule->GetID()) ? 1 : 0;
}

# Return an array of atoms. In scalar context,  it returns number of atoms. Additionally,
# atoms array can be filtered by any user specifiable Atom class method...
#
sub GetAtoms {
  my($This, $AtomCheckMethodName, $NegateMethodResult) = @_;
  my(@Atoms, @AtomIDs);

  @Atoms = (); @AtomIDs = ();

  @AtomIDs = $This->GetVertices();
  if (!@AtomIDs) {
    return wantarray ? @Atoms : scalar @Atoms;
  }

  @Atoms = $This->GetVerticesProperty('Atom', @AtomIDs);

  if (!defined $AtomCheckMethodName) {
    return wantarray ? @Atoms : scalar @Atoms;
  }
  $NegateMethodResult = (defined($NegateMethodResult) &&  $NegateMethodResult) ? 1 : 0;

  # Filter out atoms...
  my($Atom, $KeepAtom, @FilteredAtoms);
  @FilteredAtoms = ();
  for $Atom (@Atoms) {
    $KeepAtom = $NegateMethodResult ? (!$Atom->$AtomCheckMethodName()) : $Atom->$AtomCheckMethodName();
    if ($KeepAtom) {
      push @FilteredAtoms, $Atom;
    }
  }
  return wantarray ? @FilteredAtoms : scalar @FilteredAtoms;
}

# Return an array of bonds. In scalar context, it returns number of bonds...
sub GetBonds {
  my($This) = @_;
  my(@Bonds, @EdgesAtomsIDs);

  @Bonds = (); @EdgesAtomsIDs = ();

  @EdgesAtomsIDs = $This->GetEdges();
  if (@EdgesAtomsIDs) {
    @Bonds = $This->GetEdgesProperty('Bond', @EdgesAtomsIDs);
  }
  return wantarray ? @Bonds : scalar @Bonds;
}

# Get number of atoms in molecule...
sub GetNumOfAtoms {
  my($This) = @_;
  my($NumOfAtoms);

  $NumOfAtoms = $This->GetAtoms();

  return $NumOfAtoms;
}

# Get number of bonds in molecule...
sub GetNumOfBonds {
  my($This) = @_;
  my($NumOfBonds);

  $NumOfBonds = $This->GetBonds();

  return $NumOfBonds;
}

# Get number of heavy atoms in molecule...
sub GetNumOfHeavyAtoms {
  my($This) = @_;

  return $This->GetNumOfNonHydrogenAtoms();
}

# Get number of non-hydrogen atoms in molecule...
sub GetNumOfNonHydrogenAtoms {
  my($This) = @_;
  my($NumOfNonHydrogenAtoms, $Atom, @Atoms);

  @Atoms = $This->GetAtoms();
  $NumOfNonHydrogenAtoms = 0;
  for $Atom (@Atoms) {
    if (!$Atom->IsHydrogen()) {
      $NumOfNonHydrogenAtoms++;
    }
  }
  return $NumOfNonHydrogenAtoms;
}

# Get number of hydrogen atoms in molecule...
sub GetNumOfHydrogenAtoms {
  my($This) = @_;
  my($NumOfHydrogenAtoms, $Atom, @Atoms);

  @Atoms = $This->GetAtoms();
  $NumOfHydrogenAtoms = 0;
  for $Atom (@Atoms) {
    if ($Atom->IsHydrogen()) {
      $NumOfHydrogenAtoms++;
    }
  }
  return $NumOfHydrogenAtoms;
}

# Get number of missing hydrogen atoms in molecule...
sub GetNumOfMissingHydrogenAtoms {
  my($This) = @_;
  my($NumOfMissingHydrogenAtoms, $Atom, @Atoms);

  @Atoms = $This->GetAtoms();
  $NumOfMissingHydrogenAtoms = 0;
  for $Atom (@Atoms) {
    if (!$Atom->IsHydrogen()) {
      $NumOfMissingHydrogenAtoms += $Atom->GetNumOfMissingHydrogens();
    }
  }
  return $NumOfMissingHydrogenAtoms;
}

# Add bond...
sub AddBond {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();
  if (!(defined($Atom1) && defined($Atom2))) {
    carp "Warning: ${ClassName}->AddBond: No bond added: Both atoms must be specified...";
    return undef;
  }
  if (!($This->HasAtom($Atom1) && $This->HasAtom($Atom2))) {
    carp "Warning: ${ClassName}->AddBond: No bond added: Both atoms must be present...";
    return undef;
  }
  if ($This->HasBond($Bond)) {
    carp "Warning: ${ClassName}->AddBond: No bond added: Bond already exists...";
    return undef;
  }
  return $This->_AddBond($Bond);
}

# Add bond...
sub _AddBond {
  my($This, $Bond) = @_;

  # Assign bond to this molecule...
  $Bond->_SetMolecule($This);

  # Add it to the graph as an edge...
  my($Atom1, $Atom2, $AtomID1, $AtomID2);
  ($Atom1, $Atom2) = $Bond->GetAtoms();
  $AtomID1 = $Atom1->GetID(); $AtomID2 = $Atom2->GetID();
  $This->AddEdge($AtomID1, $AtomID2);
  $This->SetEdgeProperty('Bond', $Bond, $AtomID1, $AtomID2);

  return $This;
}

# Add Bonds...
sub AddBonds {
  my($This, @Bonds) = @_;

  if (!@Bonds) {
    carp "Warning: ${ClassName}->AddBonds: No bonds added: Bonds list must be specified...";
    return undef;
  }
  my($Bond);
  for $Bond (@Bonds) {
    $This->AddBond($Bond);
  }
  return $This;
}

# Create a bond and add it to molecule...
sub NewBond {
  my($This, %NamesAndValues) = @_;
  my($Bond);

  $Bond = new Bond(%NamesAndValues);
  $This->AddBond($Bond);

  return $Bond;
}

# Delete a bond...
sub DeleteBond {
  my($This, $Bond) = @_;

  if (!defined $Bond) {
    carp "Warning: ${ClassName}->DeleteBond: No bond deleted: Bond must be specified...";
    return undef;
  }
  # Does the bond exist in molecule?
  if (!$This->HasBond($Bond)) {
    carp "Warning: ${ClassName}->DeleteBond: No bond deleted: Bond doesn't exist...";
    return undef;
  }
  return $This->_DeleteBond($Bond);
}

# Delete bond...
sub _DeleteBond {
  my($This, $Bond) = @_;

  my($Atom1, $Atom2, $AtomID1, $AtomID2);
  ($Atom1, $Atom2) = $Bond->GetAtoms();
  $AtomID1 = $Atom1->GetID(); $AtomID2 = $Atom2->GetID();
  $This->DeleteEdge($AtomID1, $AtomID2);

  return $This;
}

# Delete bonds...
sub DeleteBonds {
  my($This, @Bonds) = @_;

  if (!@Bonds) {
    carp "Warning: ${ClassName}->DeleteBonds: No bonds added: Bonds list must be specified...";
    return undef;
  }
  my($Bond);
  for $Bond (@Bonds) {
    $This->DeleteBond($Bond);
  }

  return $This;
}

# Has bond...
sub HasBond {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();
  if (!($This->HasAtom($Atom1) && $This->HasAtom($Atom2))) {
    return 0;
  }
  if (!$Bond->HasProperty('Molecule')) {
    # It's not in molecule...
    return 0;
  }
  my($AtomID1, $AtomID2, $Molecule);
  $AtomID1 = $Atom1->GetID(); $AtomID2 = $Atom2->GetID();
  $Molecule = $Bond->GetMolecule();

  return ($This->HasEdge($AtomID1, $AtomID2) && $This->GetID() == $Molecule->GetID()) ? 1 : 0;
}

# Get atom neighbors...
sub _GetAtomNeighbors {
  my($This, $Atom) = @_;

  my($AtomID, @Atoms, @AtomIDs);

  @Atoms = (); @AtomIDs = ();
  $AtomID = $Atom->GetID();
  @AtomIDs = $This->GetNeighbors($AtomID);
  if (@AtomIDs) {
    @Atoms = $This->GetVerticesProperty('Atom', @AtomIDs);
  }
  return wantarray ? @Atoms : scalar @Atoms;
}

# Get atom bonds...
sub _GetAtomBonds {
  my($This, $Atom) = @_;
  my($AtomID, @AtomIDs, @Bonds);

  @Bonds = (); @AtomIDs = ();
  $AtomID = $Atom->GetID();
  @AtomIDs = $This->GetEdges($AtomID);
  if (@AtomIDs) {
    @Bonds = $This->GetEdgesProperty('Bond', @AtomIDs);
  }
  return wantarray ? @Bonds : scalar @Bonds;
}

# Get bond to specified atom...
sub _GetBondToAtom {
  my($This, $Atom1, $Atom2) = @_;
  my($AtomID1, $AtomID2);

  $AtomID1 = $Atom1->GetID();
  $AtomID2 = $Atom2->GetID();

  return $This->GetEdgeProperty('Bond', $AtomID1, $AtomID2);
}

# Are two specified atoms bonded?
sub _IsBondedToAtom {
  my($This, $Atom1, $Atom2) = @_;
  my($AtomID1, $AtomID2);

  $AtomID1 = $Atom1->GetID();
  $AtomID2 = $Atom2->GetID();

  return $This->HasEdgeProperty('Bond', $AtomID1, $AtomID2);
}

# Add hydrogens to each atoms in molecule and return total number of hydrogens added...
sub AddHydrogens {
  my($This) = @_;

  return $This->_AddHydrogens();
}

# Add hydrogens to polar atoms (N, O, P, S) in molecule and return total number of hydrogens added...
sub AddPolarHydrogens {
  my($This) = @_;
  my($PolarHydrogensOnly) = 1;

  return $This->_AddHydrogens($PolarHydrogensOnly);
}

# Add all the hydrogens or hydrogens for polar atoms only...
#
# Note:
#   . The current release of MayaChemTools doesn't assign any hydrogen positions.
#
sub _AddHydrogens {
  my($This, $PolarHydrogensOnly) = @_;
  my($Atom, $NumOfHydrogensAdded, $HydrogenPositionsWarning, @Atoms);

  if (! defined $PolarHydrogensOnly) {
    $PolarHydrogensOnly = 0;
  }

  $NumOfHydrogensAdded = 0;
  @Atoms = $This->GetAtoms();
  $HydrogenPositionsWarning = 0;

  ATOM: for $Atom (@Atoms) {
    if ($PolarHydrogensOnly) {
      if (!$Atom->IsPolarAtom()) {
	next ATOM;
      }
    }
    $NumOfHydrogensAdded += $Atom->AddHydrogens($HydrogenPositionsWarning);
  }
  return $NumOfHydrogensAdded;
}

# Delete all hydrogens atoms in molecule and return total number of hydrogens removed...
sub DeleteHydrogens {
  my($This) = @_;

  return $This->_DeleteHydrogens();
}

# Delete hydrogens to polar atoms (N, O, P, S) in molecule and return total number of hydrogens removed...
sub DeletePolarHydrogens {
  my($This) = @_;
  my($PolarHydrogensOnly) = 1;

  return $This->_DeleteHydrogens($PolarHydrogensOnly);
}

# Delete all hydrogens atoms in molecule and return total number of hydrogens removed...
sub _DeleteHydrogens {
  my($This, $PolarHydrogensOnly) = @_;
  my($Atom, $NumOfHydrogensRemoved, @Atoms);

  if (! defined $PolarHydrogensOnly) {
    $PolarHydrogensOnly = 0;
  }

  $NumOfHydrogensRemoved = 0;
  @Atoms = $This->GetAtoms();

  ATOM: for $Atom (@Atoms) {
    if ($PolarHydrogensOnly) {
      if (!$Atom->IsPolarHydrogen()) {
	next ATOM;
      }
    }
    elsif (!$Atom->IsHydrogen()) {
      next ATOM;
    }
    $This->_DeleteAtom($Atom);
    $NumOfHydrogensRemoved++;
  }
  return $NumOfHydrogensRemoved;
}

# Get molecular weight by summing up atomic weights of all the atoms...
sub GetMolecularWeight {
  my($This, $IncludeMissingHydrogens) = @_;
  my($MolecularWeight, $AtomicWeight, @Atoms, $Atom);

  $IncludeMissingHydrogens = defined($IncludeMissingHydrogens) ? $IncludeMissingHydrogens : 1;

  $MolecularWeight = 0;
  @Atoms = $This->GetAtoms();
  for $Atom (@Atoms) {
    $AtomicWeight = $Atom->GetAtomicWeight();
    if (defined $AtomicWeight) {
      $MolecularWeight += $AtomicWeight;
    }
  }

  if (!$IncludeMissingHydrogens) {
    return $MolecularWeight;
  }

  # Account for missing hydrogen atoms...
  my($NumOfMissingHydrogenAtoms);

  $NumOfMissingHydrogenAtoms = $This->GetNumOfMissingHydrogenAtoms();
  if ($NumOfMissingHydrogenAtoms) {
    $MolecularWeight += $NumOfMissingHydrogenAtoms * PeriodicTable::GetElementAtomicWeight('H');
  }

  return $MolecularWeight;
}

# Get exact mass by summing up exact masses of all the atoms...
sub GetExactMass {
  my($This, $IncludeMissingHydrogens) = @_;
  my($ExactMass, $AtomicMass, @Atoms, $Atom);

  $IncludeMissingHydrogens = defined($IncludeMissingHydrogens) ? $IncludeMissingHydrogens : 1;

  $ExactMass = 0;
  @Atoms = $This->GetAtoms();
  for $Atom (@Atoms) {
    $AtomicMass = $Atom->GetExactMass();
    if (defined $AtomicMass) {
      $ExactMass += $AtomicMass;
    }
  }

  if (!$IncludeMissingHydrogens) {
    return $ExactMass;
  }

  # Account for missing hydrogen atoms...
  my($NumOfMissingHydrogenAtoms);

  $NumOfMissingHydrogenAtoms = $This->GetNumOfMissingHydrogenAtoms();
  if ($NumOfMissingHydrogenAtoms) {
    $ExactMass += $NumOfMissingHydrogenAtoms * PeriodicTable::GetElementMostAbundantNaturalIsotopeMass('H');
  }

  return $ExactMass;
}

# Get net formal charge on the molecule using one of the following two methods:
#   . Using explicitly set FormalCharge property
#   . Adding up formal charge on each atom in the molecule
#
# Caveats:
#   . FormalCharge is different from Charge property of the molecule which corresponds to
#     sum of partial atomic charges explicitly set for each atom using a specific methodology.
#
sub GetFormalCharge {
  my($This) = @_;

  # Is FormalCharge property explicitly set?
  if ($This->HasProperty('FormalCharge')) {
    return $This->GetProperty('FormalCharge');
  }
  my($FormalCharge, $AtomicFormalCharge, @Atoms, $Atom);

  $FormalCharge = 0;
  @Atoms = $This->GetAtoms();
  for $Atom (@Atoms) {
    $AtomicFormalCharge = $Atom->GetFormalCharge();
    if (defined $AtomicFormalCharge) {
      $FormalCharge += $AtomicFormalCharge;
    }
  }
  return $FormalCharge;
}

# Get net charge on the molecule using one of the following two methods:
#   . Using explicitly set Charge property
#   . Adding up charge on each atom in the molecule
#
# Caveats:
#   . FormalCharge is different from Charge property of the molecule which corresponds to
#     sum of partial atomic charges explicitly set for each atom using a specific methodology.
#
sub GetCharge {
  my($This) = @_;

  # Is Charge property explicitly set?
  if ($This->HasProperty('Charge')) {
    return $This->GetProperty('Charge');
  }
  my($Charge, $AtomicCharge, @Atoms, $Atom);

  $Charge = 0;
  @Atoms = $This->GetAtoms();
  for $Atom (@Atoms) {
    $AtomicCharge = $Atom->GetCharge();
    if (defined $AtomicCharge) {
      $Charge += $AtomicCharge;
    }
  }
  return $Charge;
}

# Get total SpinMultiplicity for the molecule using one of the following two methods:
#   . Using explicitly set SpinMultiplicity property
#   . Adding up SpinMultiplicity on each atom in the molecule
#
#
sub GetSpinMultiplicity {
  my($This) = @_;

  # Is SpinMultiplicity property explicitly set?
  if ($This->HasProperty('SpinMultiplicity')) {
    return $This->GetProperty('SpinMultiplicity');
  }
  my($AtomicSpinMultiplicity, $SpinMultiplicity, @Atoms, $Atom);

  $SpinMultiplicity = 0;
  @Atoms = $This->GetAtoms();
  for $Atom (@Atoms) {
    $AtomicSpinMultiplicity = $Atom->GetSpinMultiplicity();
    if (defined $AtomicSpinMultiplicity) {
      $SpinMultiplicity += $AtomicSpinMultiplicity;
    }
  }
  return $SpinMultiplicity;
}

# Get total FreeRadicalElectrons for the molecule using one of the following two methods:
#   . Using explicitly set FreeRadicalElectrons property
#   . Adding up FreeRadicalElectrons on each atom in the molecule
#
#
sub GetFreeRadicalElectrons {
  my($This) = @_;

  # Is FreeRadicalElectrons property explicitly set?
  if ($This->HasProperty('FreeRadicalElectrons')) {
    return $This->GetProperty('FreeRadicalElectrons');
  }
  my($AtomicFreeRadicalElectrons, $FreeRadicalElectrons, @Atoms, $Atom);

  $FreeRadicalElectrons = 0;
  @Atoms = $This->GetAtoms();
  for $Atom (@Atoms) {
    $AtomicFreeRadicalElectrons = $Atom->GetFreeRadicalElectrons();
    if (defined $AtomicFreeRadicalElectrons) {
      $FreeRadicalElectrons += $AtomicFreeRadicalElectrons;
    }
  }
  return $FreeRadicalElectrons;
}

# Set valence model...
#
sub SetValenceModel {
  my($This, $ValenceModel) = @_;

  if ($ValenceModel !~ /^(MDLValenceModel|DaylightValenceModel|InternalValenceModel|MayaChemToolsValenceModel)$/i) {
    carp "Warning: ${ClassName}->SetValenceModel: The current release of MayaChemTools doesn't support the specified valence model $ValenceModel. Supported valence models: MDLValenceModel, DaylightValenceModel, InternalValenceModel or MayaChemToolsValenceModel. Using internal valence model...";
    $ValenceModel = 'InternalValenceModel';
  }

  $This->SetProperty('ValenceModel', $ValenceModel);

  return $This;
}

# Get valence model...
#
sub GetValenceModel {
  my($This) = @_;

  # Is ValenceModel property explicitly set?
  if ($This->HasProperty('ValenceModel')) {
    return $This->GetProperty('ValenceModel');
  }

  # Used internal valence model as default model...
  return 'InternalValenceModel';
}

# Get molecular formula by collecting information about all atoms in the molecule and
# composing the formula using Hills ordering system:
#   . C shows up first and H follows assuming C is present
#   . All other standard elements are sorted alphanumerically
#   . All other non-stanard atom symbols are also sorted alphanumerically and follow standard elements
#
# Caveats:
#   . By default, missing hydrogens and nonelements are also included
#   . Elements for disconnected fragments are combined into the same formula
#
# Handle formula generation for disconnected structures. e.g: molecule generated by
# [Na+].[O-]c1ccccc1
#
sub GetMolecularFormula {
  my($This, $IncludeMissingHydrogens, $IncludeNonElements) = @_;
  my($MolecularFormula, $AtomSymbol, $ElementsCountRef, $NonElementsCountRef);

  $IncludeMissingHydrogens = defined($IncludeMissingHydrogens) ? $IncludeMissingHydrogens : 1;
  $IncludeNonElements = defined($IncludeNonElements) ? $IncludeNonElements : 1;

  # Get elements count and setup molecular formula...
  ($ElementsCountRef, $NonElementsCountRef) = $This->GetElementsAndNonElements($IncludeMissingHydrogens);
  $MolecularFormula = '';

  # Count C and H first...
  if (exists $ElementsCountRef->{C} ) {
    $MolecularFormula .= 'C' . ($ElementsCountRef->{C} > 1 ? $ElementsCountRef->{C} : '');
    delete $ElementsCountRef->{C};

    if (exists $ElementsCountRef->{H} ) {
      $MolecularFormula .= 'H' . ($ElementsCountRef->{H} > 1 ? $ElementsCountRef->{H} : '');
      delete $ElementsCountRef->{H};
    }
  }

  # Sort elements...
  for $AtomSymbol (sort {$a cmp $b} keys %{$ElementsCountRef}) {
    $MolecularFormula .= $AtomSymbol . ($ElementsCountRef->{$AtomSymbol} > 1 ? $ElementsCountRef->{$AtomSymbol} : '');
  }

  # Sort non-elements...
  if ($IncludeNonElements) {
    for $AtomSymbol (sort {$a cmp $b} keys %{$NonElementsCountRef}) {
      $MolecularFormula .= $AtomSymbol . ($NonElementsCountRef->{$AtomSymbol} > 1 ? $NonElementsCountRef->{$AtomSymbol} : '');
    }
  }

  # Check formal charge...
  my($FormalCharge);
  $FormalCharge = $This->GetFormalCharge();
  if ($FormalCharge) {
    # Setup formal charge string...
    my($FormalChargeString);
    if ($FormalCharge == 1 ) {
      $FormalChargeString =  "+";
    }
    elsif ($FormalCharge == -1 ) {
      $FormalChargeString =  "-";
    }
    else {
      $FormalChargeString = ($FormalCharge > 0) ? ("+" . abs($FormalCharge)) : ("-" . abs($FormalCharge));
    }
    $MolecularFormula = "${MolecularFormula}${FormalChargeString}";
  }

  return $MolecularFormula;
}

# Count elements and non-elements in molecule and return references to hashes
# containing count of elements and non-elements respectively. By default, missing
# hydrogens are not added to the element hash.
#
#
sub GetElementsAndNonElements {
  my($This, $IncludeMissingHydrogens) = @_;
  my($Atom, $AtomicNumber, $AtomSymbol, $NumOfMissingHydrogens, @Atoms, %ElementsCount, %NonElementsCount);

  $IncludeMissingHydrogens = defined($IncludeMissingHydrogens) ? $IncludeMissingHydrogens : 0;

  %ElementsCount = (); %NonElementsCount = ();
  $NumOfMissingHydrogens = 0;

  # Count elements and non elements...
  @Atoms = $This->GetAtoms();
  for $Atom (@Atoms) {
    $AtomicNumber = $Atom->GetAtomicNumber();
    $AtomSymbol = $Atom->GetAtomSymbol();
    if ($AtomicNumber) {
      if (exists $ElementsCount{$AtomSymbol}) {
	$ElementsCount{$AtomSymbol} += 1;
      }
      else {
	$ElementsCount{$AtomSymbol} = 1;
      }
      if ($IncludeMissingHydrogens) {
	$NumOfMissingHydrogens += $Atom->GetNumOfMissingHydrogens();
      }
    }
    else {
      if (exists $NonElementsCount{$AtomSymbol}) {
	$NonElementsCount{$AtomSymbol} += 1;
      }
      else {
	$NonElementsCount{$AtomSymbol} = 1;
      }
    }
  }
  if ($IncludeMissingHydrogens && $NumOfMissingHydrogens) {
    $AtomSymbol = 'H';
    if (exists $ElementsCount{$AtomSymbol}) {
      $ElementsCount{$AtomSymbol} += $NumOfMissingHydrogens;
    }
    else {
      $ElementsCount{$AtomSymbol} = $NumOfMissingHydrogens;
    }
  }

  return (\%ElementsCount, \%NonElementsCount);
}

# Get number of element and non-elements in molecule. By default, missing
# hydrogens are not added to element count.
#
sub GetNumOfElementsAndNonElements {
  my($This, $IncludeMissingHydrogens) = @_;
  my($ElementCount, $NonElementCount, $Atom);

  $IncludeMissingHydrogens = defined($IncludeMissingHydrogens) ? $IncludeMissingHydrogens : 0;

  ($ElementCount, $NonElementCount) = (0, 0);

  ATOM: for $Atom ($This->GetAtoms()) {
    if (!$Atom->GetAtomicNumber()) {
      $NonElementCount++;
      next ATOM;
    }
    # Process element...
    $ElementCount++;
    if ($IncludeMissingHydrogens) {
      if (!$Atom->IsHydrogen()) {
	$ElementCount += $Atom->GetNumOfMissingHydrogens();
      }
    }
  }

  return ($ElementCount, $NonElementCount);
}

# Calculate elemental composition and return reference to arrays
# containing elements and their percent composition.
#
# Caveats:
#   . By default, missing hydrogens are included
#   . Non elemnents are ignored
#   . Mass number are ignored
#
sub GetElementalComposition {
  my($This, $IncludeMissingHydrogens) = @_;
  my($MolecularFormula, $IncludeNonElements, $ElementsCountRef, $NonElementsCountRef, $ElementsRef, $ElementsCompositionRef);

  $IncludeMissingHydrogens = defined($IncludeMissingHydrogens) ? $IncludeMissingHydrogens : 1;

  $IncludeNonElements = 0;
  ($ElementsCountRef, $NonElementsCountRef) = $This->GetElementsAndNonElements($IncludeMissingHydrogens);

  $MolecularFormula = $This->GetMolecularFormula($IncludeMissingHydrogens, $IncludeNonElements);

  ($ElementsRef, $ElementsCompositionRef) = MolecularFormula::CalculateElementalComposition($MolecularFormula);

  return ($ElementsRef, $ElementsCompositionRef);
}

# Using refernece to element and its composition arrays, format composition information
# as: Element: Composition;...
#
sub FormatElementalCompositionInformation {
  my($FirstParameter, $SecondParameter, $ThirdParameter, $FourthParameter) = @_;
  my($This, $ElementsRef, $ElementCompositionRef, $Precision);

  if (_IsMolecule($FirstParameter)) {
    ($This, $ElementsRef, $ElementCompositionRef, $Precision) = ($FirstParameter, $SecondParameter, $ThirdParameter, $FourthParameter);
  }
  else {
    ($ElementsRef, $ElementCompositionRef, $Precision) = ($FirstParameter, $SecondParameter, $ThirdParameter);
  }
  my($FormattedInfo) = '';

  if (!(defined($ElementsRef) && @{$ElementsRef})) {
    carp "Warning: ${ClassName}->FormatElementalCompositionInformation: Elements list is not defined or empty...";
    return undef;
  }
  if (!(defined($ElementCompositionRef) && @{$ElementCompositionRef})) {
    carp "Warning: ${ClassName}->FormatElementalCompositionInformation: Elements composition list is not defined or empty...";
    return undef;
  }

  if (!defined $Precision) {
    $Precision = 2;
  }

  $FormattedInfo = MolecularFormula::FormatCompositionInfomation($ElementsRef, $ElementCompositionRef, $Precision);

  return $FormattedInfo;
}

# Copy molecule and its associated data using Storable::dclone and update:
#
#   o Atom references corresponding atoms and bonds objects
#   o Bond object references
#
# Object IDs for Atoms and Bonds don't get changed. So there is no need to clear
# up any exisiting ring data attached to molecule via graph as vertex IDs.
#
sub Copy {
  my($This) = @_;
  my($NewMolecule, $Atom, $NewAtom, $AtomID, @Atoms, @AtomIDs, %AtomsIDsToNewAtoms);

  $NewMolecule = Storable::dclone($This);

  # Update atom references stored as vertex property...

  @Atoms = (); @AtomIDs = ();
  %AtomsIDsToNewAtoms = ();

  @AtomIDs = $This->GetVertices();
  if (@AtomIDs) {
    @Atoms = $This->GetVerticesProperty('Atom', @AtomIDs);
  }

  for $Atom (@Atoms) {
    $AtomID = $Atom->GetID();

    # Setup a reference to copied atom object...
    $NewAtom = $Atom->Copy();
    $AtomsIDsToNewAtoms{$AtomID} = $NewAtom;

    # Update atom reference to new atom object...
    $NewMolecule->UpdateVertexProperty('Atom', $NewAtom, $AtomID);
  }

  # Update bond object and bond atom references stored as edge property...
  my($Index, $AtomID1, $AtomID2, $Bond, $NewBond, $NewAtom1, $NewAtom2, @EdgesAtomsIDs);
  @EdgesAtomsIDs = ();
  @EdgesAtomsIDs = $This->GetEdges();
  for ($Index = 0; $Index < $#EdgesAtomsIDs; $Index += 2) {
    $AtomID1 = $EdgesAtomsIDs[$Index]; $AtomID2 = $EdgesAtomsIDs[$Index + 1];

    # Get reference to bond object...
    $Bond = $This->GetEdgeProperty('Bond', $AtomID1, $AtomID2);

    # Make a new bond object and update its atom references...
    $NewBond = $Bond->Copy();
    $NewAtom1 = $AtomsIDsToNewAtoms{$AtomID1};
    $NewAtom2 = $AtomsIDsToNewAtoms{$AtomID2};
    $NewBond->SetAtoms($NewAtom1, $NewAtom2);

    # Update bond object reference in the new molecule...
    $NewMolecule->UpdateEdgeProperty('Bond', $NewBond, $AtomID1, $AtomID2);
  }

  return $NewMolecule;
}

# Get number of connected components...
#
sub GetNumOfConnectedComponents {
  my($This) = @_;
  my($NumOfComponents);

  $NumOfComponents = $This->GetConnectedComponentsVertices();

  return $NumOfComponents;
}

# Return a reference to an array containing molecules corresponding
# to connected components sorted in decreasing order of component size...
#
sub GetConnectedComponents {
  my($This) = @_;
  my($Index, @ComponentMolecules, @ConnectedComponents);

  @ConnectedComponents = ();
  @ConnectedComponents = $This->GetConnectedComponentsVertices();
  @ComponentMolecules = ();

  for $Index (0 .. $#ConnectedComponents) {
    push @ComponentMolecules, $This->_GetConnectedComponent(\@ConnectedComponents, $Index);
  }
  return @ComponentMolecules;
}

# Return a reference to largest connected component as a molecule object...
#
sub GetLargestConnectedComponent {
  my($This) = @_;
  my($LargestComponentIndex, @ConnectedComponents);

  $LargestComponentIndex = 0;
  @ConnectedComponents = ();
  @ConnectedComponents = $This->GetConnectedComponentsVertices();

  return $This->_GetConnectedComponent(\@ConnectedComponents, $LargestComponentIndex);
}

# Return connected component as a molecule...
#
sub _GetConnectedComponent {
  my($This, $ConnectedComponentsRef, $ComponentIndex) = @_;
  my($ComponentMolecule);

  # Copy existing molecule...
  $ComponentMolecule = $This->Copy();

  # Delete all atoms besides the ones in specified component...
  $ComponentMolecule->_DeleteConnectedComponents($ConnectedComponentsRef, $ComponentIndex);

  # Clear any deteced rings...
  if ($ComponentMolecule->HasRings()) {
    $ComponentMolecule->ClearRings();
  }
  return $ComponentMolecule;
}

# Delete atoms corresponding to all connected components except the one specified...
#
sub _DeleteConnectedComponents {
  my($This, $ConnectedComponentsRef, $KeepComponentIndex) = @_;
  my($Index, $AtomID);

  INDEX: for $Index (0 .. $#{$ConnectedComponentsRef}) {
    if ($Index == $KeepComponentIndex) {
      next INDEX;
    }
    for $AtomID (@{$ConnectedComponentsRef->[$Index]}) {
      $This->DeleteVertex($AtomID);
    }
  }
  return $This;
}

# Return an array containing references to atom arrays corresponding to atoms of
# connected components sorted in order of their decreasing size...
#
sub GetConnectedComponentsAtoms {
  my($This) = @_;
  my($Index, @ComponentsAtoms, @ConnectedComponents);

  @ConnectedComponents = ();
  @ConnectedComponents = $This->GetConnectedComponentsVertices();

  @ComponentsAtoms = ();
  for $Index (0 .. $#ConnectedComponents) {
    my(@ComponentAtoms);

    @ComponentAtoms = ();
    @ComponentAtoms = $This->_GetConnectedComponentAtoms(\@ConnectedComponents, $Index);
    push @ComponentsAtoms, \@ComponentAtoms;
  }
  return @ComponentsAtoms;
}

# Return an array containing atoms correspondig to largest connected component...
#
sub GetLargestConnectedComponentAtoms {
  my($This) = @_;
  my($LargestComponentIndex, @ConnectedComponents);

  $LargestComponentIndex = 0;
  @ConnectedComponents = ();
  @ConnectedComponents = $This->GetConnectedComponentsVertices();

  return $This->_GetConnectedComponentAtoms(\@ConnectedComponents, $LargestComponentIndex);
}

# Return an array containing atoms corresponding to specified connected component...
#
sub _GetConnectedComponentAtoms {
  my($This, $ConnectedComponentsRef, $ComponentIndex) = @_;
  my($AtomID, @AtomIDs, @ComponentAtoms);

  @ComponentAtoms = ();
  @AtomIDs = ();

  for $AtomID (@{$ConnectedComponentsRef->[$ComponentIndex]}) {
    push @AtomIDs, $AtomID;
  }
  @ComponentAtoms = $This->_GetAtomsFromAtomIDs(@AtomIDs);

  return @ComponentAtoms;
}

# Except for the largest connected component, delete atoms corresponding to all other
# connected components...
#
sub KeepLargestComponent {
  my($This) = @_;
  my($LargestComponentIndex, @ConnectedComponents);

  @ConnectedComponents = ();
  @ConnectedComponents = $This->GetConnectedComponentsVertices();
  if (@ConnectedComponents == 1) {
    return $This;
  }
  $LargestComponentIndex = 0;
  $This->_DeleteConnectedComponents(\@ConnectedComponents, $LargestComponentIndex);

  # Clear any deteced rings...
  if ($This->HasRings()) {
    $This->ClearRings();
  }

  return $This;
}

# Get an array of topologically sorted atoms starting from a specified atom or
# an arbitrary atom in the molecule...
#
sub GetTopologicallySortedAtoms {
  my($This, $StartAtom) = @_;
  my(@SortedAtoms);

  @SortedAtoms = ();
  if (defined($StartAtom) && !$This->HasAtom($StartAtom)) {
    carp "Warning: ${ClassName}->_GetTopologicallySortedAtoms: No atoms retrieved: Start atom doesn't exist...";
    return @SortedAtoms;
  }
  my($StartAtomID, @AtomIDs);

  @AtomIDs = ();
  $StartAtomID = defined($StartAtom) ? $StartAtom->GetID() : undef;

  @AtomIDs = $This->GetTopologicallySortedVertices($StartAtomID);
  @SortedAtoms = $This->_GetAtomsFromAtomIDs(@AtomIDs);

  return @SortedAtoms;
}

# Detect rings in molecule...
#
sub DetectRings {
  my($This) = @_;

  # Use graph method to detect all cycles and associate 'em to graph as graph
  # and vertex properties...
  return $This->DetectCycles();
}

# Clear rings in molecule...
#
sub ClearRings {
  my($This) = @_;

  # Use graph method to clear all cycles...
  $This->ClearCycles();

  return $This;
}

# Setup rings type paths to use during all ring related methods. Possible values:
# Independent or All. Default is to use Independent rings.
#
sub SetActiveRings {
  my($This, $RingsType) = @_;

  if (!defined $This->SetActiveCyclicPaths($RingsType)) {
    return undef;
  }
  return $This;
}

# Is it a supported aromaticity model?
#
sub IsSupportedAromaticityModel {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $AromaticityModel);

  if (_IsMolecule($FirstParameter)) {
    ($This, $AromaticityModel) = ($FirstParameter, $SecondParameter);
  }
  else {
    ($This, $AromaticityModel) = (undef, $FirstParameter);
  }

  return exists $CanonicalAromaticityModelNamesMap{lc($AromaticityModel)} ? 1 : 0;
}

# Get a list of supported aromaticity model names...
#
sub GetSupportedAromaticityModels {
  return (sort values %CanonicalAromaticityModelNamesMap);
}

# Set aromaticity model...
#
sub SetAromaticityModel {
  my($This, $AromaticityModel) = @_;

  if (!$This->IsSupportedAromaticityModel($AromaticityModel)) {
    my(@SupportedModels) = $This->GetSupportedAromaticityModels();

    carp "Warning: ${ClassName}->SetAromaticityModel: The current release of MayaChemTools doesn't support the specified aromaticity model $AromaticityModel. Supported aromaticity models defined in AromaticityModelsData.csv file are: @SupportedModels . Using MayaChemToolsAromaticityModel...";
    $AromaticityModel = 'MayaChemToolsAromaticityModel';
  }

  $This->SetProperty('AromaticityModel', $AromaticityModel);

  return $This;
}

# Get aromaticity model...
#
sub GetAromaticityModel {
  my($This) = @_;

  # Is ValenceModel property explicitly set?
  if ($This->HasProperty('AromaticityModel')) {
    return $This->GetProperty('AromaticityModel');
  }

  # Used internal aromaticity model as default model...
  return 'MayaChemToolsAromaticityModel';
}

# Identify aromatic rings and ring systems in a molecule and set aromaticity for
# corresponding atoms and bonds.
#
# What is aromaticity? [ Ref 124 ] It's in the code of the implementer, did you
# say? Agree. The implementation of aromaticity varies widely across different
# packages [ Ref 125 ]; additionally, the implementation details are not always
# completely available, and it's not possible to figure out the exact implementation
# of aromaticity across various packages. Using the publicly available information,
# however, one can try to reproduce the available results to the extent possible,
# along with parameterizing all the control parameters used to implement different
# aromaticity models, and that's exactly what the current release of MayaChemTools
# does.
#
# The implementation of aromaticity corresponding to various aromaticity models in
# MayaChemTools package is driven by an external CSV file AromaticityModelsData.csv,
# which is distributed with the package and is available in lib/data directory. The CSV
# files contains names of supported aromaticity models, along with various control
# parameters and their values. This file is loaded and processed during instantiation
# of Molecule class and data corresponding to specific aromaticity model are used
# to detect aromaticity for that model. Any new aromaticity model added to the
# aromaticity data file, using different combinations of values for existing control
# parameters would work without any changes to the code; the addition of any new
# control parameters, however, requires its implementation in the code used to
# calculate number of pi electrons available towards delocalization in a ring or ring
# systems.
#
# The current release of MayaChemTools package supports these aromaticity
# models: MDLAromaticityModel, TriposAromaticityModel, MMFFAromaticityModel,
# ChemAxonBasicAromaticityModel, ChemAxonGeneralAromaticityModel,
# DaylightAromaticityModel, MayaChemToolsAromaticityModel.
#
# The current list of control parameters available to detect aromaticity corresponding
# to different aromaticity models are: AllowHeteroRingAtoms, HeteroRingAtomsList,
# AllowExocyclicDoubleBonds, AllowHomoNuclearExocyclicDoubleBonds,
# AllowElectronegativeRingAtomExocyclicDoubleBonds, AllowRingAtomFormalCharge,
# AllowHeteroRingAtomFormalCharge, MinimumRingSize. The values for these control
# parameters are specified in AromaticityModelsData.csv file.
#
# Although definition of aromaticity differs across various aromaticity models, a ring
# or a ring system containing 4n + 2 pi electrons (Huckel's rule) corresponding to
# alternate single and double bonds, in general, is considered aromatic.
#
# The available valence free electrons on heterocyclic ring atoms, involved in two single
# ring bonds, are also allowed to participate in pi electron delocalizaiton for most of
# the supported aromaticity models.
#
# The presence of exocyclic terminal double bond on ring atoms involved in pi electron
# delocalization is only allowed for some of the aromaticity models. Additionally, the type
# atoms involved in exocyclic terminal double bonds may result in making a ring or ring
# system non-aromatic.
#
sub DetectAromaticity {
  my($This) = @_;

  # Delete aromaticity property for atoms and bonds...
  $This->DeleteAromaticity();

  # Any ring out there...
  if (!$This->HasRings()) {
    return $This;
  }

  if ($This->HasFusedRings()) {
    $This->_DetectAromaticityUsingFusedRingSets();
  }
  else {
    $This->_DetectAromaticityUsingIndividualRings();
  }
  return $This;
}

# Go over all rings and set aromaticity property for corresponding ring atoms
# and bonds involved in aromatic rings...
#
sub _DetectAromaticityUsingIndividualRings {
  my($This) = @_;

  return $This->_DetectRingsAromaticity($This->GetRings());
}

# For each fused ring set, detect aromaticity by treating all of its ring as one aromatic
# system for counting pi electrons to satisfy Huckel's rule; In case of a failure, rings in
# fused set are treated individually for aromaticity detection. Additionally, non-fused
# rings are handled on their own during aromaticity detection.
#
# Note:
#   . pi electrons in common bonds involved in fused ring sets are only counted once.
#
#
sub _DetectAromaticityUsingFusedRingSets {
  my($This) = @_;
  my($Index, $RingAtomsRef, $FusedRingSetRef, $FusedRingSetsRef, $NonFusedRingsRef, @FusedRingSetIsAromatic);

  ($FusedRingSetsRef, $NonFusedRingsRef) = $This->GetFusedAndNonFusedRings();

  @FusedRingSetIsAromatic = ();
  RINGSET: for $Index (0 .. $#{$FusedRingSetsRef}) {
    $FusedRingSetRef = $FusedRingSetsRef->[$Index];
    $FusedRingSetIsAromatic[$Index] = 0;

    my($NumOfPiElectronsInRingSet, $NumOfPiElectronsInRing, $Bond, $BondID, %FusedRingSetsBondsMap, %FusedRingSetsBondsVisitedMap, %FusedRingBondsMap);

    $NumOfPiElectronsInRingSet = 0;

    %FusedRingSetsBondsMap = ();
    %FusedRingSetsBondsVisitedMap = ();
    %FusedRingBondsMap = ();

    # Setup a bond ID map for all bonds in fused ring set and another one
    # for bonds involved in more than one ring...
    #
    for $RingAtomsRef (@{$FusedRingSetRef}) {
      for $Bond ($This->GetRingBonds(@{$RingAtomsRef})) {
	$BondID = $Bond->GetID();
	$FusedRingSetsBondsMap{$BondID} = $BondID;

	if ($Bond->GetNumOfRings() == 2) {
	  $FusedRingBondsMap{$BondID} = $BondID;
	}
      }
    }

    for $RingAtomsRef (@{$FusedRingSetRef}) {
      my(@RingBonds);

      @RingBonds = ();
      @RingBonds = $This->GetRingBonds(@{$RingAtomsRef});
      $NumOfPiElectronsInRing = $This->_GetNumOfPiElectronsAvailableForDelocalization($RingAtomsRef, \@RingBonds, \%FusedRingSetsBondsMap, \%FusedRingSetsBondsVisitedMap, \%FusedRingBondsMap);

      if (!$NumOfPiElectronsInRing) {
	next RINGSET;
      }
      $NumOfPiElectronsInRingSet += $NumOfPiElectronsInRing;
    }
    if ($This->_DoPiElectronSatifyHuckelsRule($NumOfPiElectronsInRingSet)) {
      $FusedRingSetIsAromatic[$Index] = 1;
    }
  }

  # Set atom and bond aromatic flags for ring sets whose pi electrons satisfy Huckel's rule; otherwise,
  # treat rings in a ring set as individual rings for detecting aromaticity...
  for $Index (0 .. $#{$FusedRingSetsRef}) {
    $FusedRingSetRef = $FusedRingSetsRef->[$Index];
    if ($FusedRingSetIsAromatic[$Index]) {
      $This->_SetRingsAromaticity(@{$FusedRingSetRef});
    }
    else {
      $This->_DetectRingsAromaticity(@{$FusedRingSetRef});
    }
  }

  $This->_DetectRingsAromaticity(@{$NonFusedRingsRef});

  return $This;
}

# Detect and set aromaticity for rings...
#
sub _DetectRingsAromaticity {
  my($This, @Rings) = @_;
  my($RingAtom, $RingBond, $RingAtomsRef);

  RING: for $RingAtomsRef (@Rings) {
    if (!$This->_CheckRingAromaticity(@{$RingAtomsRef})) {
      next RING;
    }
    $This->_SetRingAromaticity(@{$RingAtomsRef});
  }
  return $This;
}

# Set aromatic property for all all atoms and bonds involved in all specified rings..
#
sub _SetRingsAromaticity {
  my($This, @Rings) = @_;
  my($RingAtomsRef );

  for $RingAtomsRef (@Rings) {
    $This->_SetRingAromaticity(@{$RingAtomsRef});
  }
  return $This;
}

# Set aromatic property for all all atoms and bonds involved in ring..
#
sub _SetRingAromaticity {
  my($This, @RingAtoms) = @_;
  my($RingAtom, $RingBond);

  for $RingAtom (@RingAtoms) {
    $RingAtom->SetAromatic(1);
  }
  for $RingBond ($This->GetRingBonds(@RingAtoms)) {
    $RingBond->SetAromatic(1);
  }
  return $This;
}


# For a ring to be an aromatic ring, all of its atoms must have aromatic property
# set.
#
sub IsRingAromatic {
  my($This, @RingAtoms) = @_;
  my($RingAtom);

  for $RingAtom (@RingAtoms) {
    if (!$RingAtom->IsAromatic()) {
      return 0;
    }
  }
  return 1;
}

# Delete aromatic property for all atoms and bonds...
#
sub DeleteAromaticity {
  my($This) = @_;

  return $This->_DeleteAtomsAndBondsAromaticity();
}

# Check ring aromaticity...
#
sub _CheckRingAromaticity {
  my($This, @RingAtoms) = @_;
  my($NumOfPiElectrons, $BondID, @RingBonds);

  @RingBonds = ();
  @RingBonds = $This->GetRingBonds(@RingAtoms);

  $NumOfPiElectrons = $This->_GetNumOfPiElectronsAvailableForDelocalization(\@RingAtoms, \@RingBonds);

  return $This->_DoPiElectronSatifyHuckelsRule($NumOfPiElectrons);
}

# Get number of pi electrons available for delocalizaiton in a ring or ring system...
#
sub _GetNumOfPiElectronsAvailableForDelocalization {
  my($This, $RingAtomsRef, $RingBondsRef, $FusedRingSetsBondsMapRef, $FusedRingSetsBondsVisitedMapRef, $FusedRingBondsMapRef) = @_;
  my($AromaticityModelName, $AromaticityModelDataRef, $ExocyclicDoubleBondsDataMapRef, $NumOfConjugatedDoubleBonds, $NumOfExplicitAromaticBonds, $NumOfRingAtomElectronPairs, $NumOfRingBondsProcessed, $NumOfPiElectrons, $Index, $RingBond, $RingAtom, $RingAtomSymbol, $BondOrder, $RingAtomID, $RingBondID, $PreviousIndex, $PreviousRingBond, $ExcludeFreeRadicalElectrons, %ElectronPairContributionProcessedMap);

  $AromaticityModelName = $CanonicalAromaticityModelNamesMap{lc($This->GetAromaticityModel())};
  $AromaticityModelDataRef = \%{$AromaticityModelsDataMap{$AromaticityModelName}};

  # Perform an intial check for potential atomatic ring atoms..
  if (!$This->_CheckRingAtomsForPotentialAromaticity($RingAtomsRef, $RingBondsRef, $AromaticityModelDataRef)) {
    return 0;
  }

  $ExocyclicDoubleBondsDataMapRef = undef;
  $ExcludeFreeRadicalElectrons = 1;

  %ElectronPairContributionProcessedMap = ();

  ($NumOfPiElectrons, $NumOfRingBondsProcessed, $NumOfConjugatedDoubleBonds, $NumOfExplicitAromaticBonds, $NumOfRingAtomElectronPairs) = ('0') x 5;

  # Go over ring atoms and bonds to check their participation in aromaticity and count
  # pi electrons available for delocalization corresponding to various aromaticity models...
  #
  RINGBOND: for $Index (0 .. $#{$RingBondsRef}) {
    $RingBond = $RingBondsRef->[$Index];
    $RingAtom = $RingAtomsRef->[$Index];
    $BondOrder = $RingBond->GetBondOrder();

    # Is this ring bond part of a fused ring system which has been already processed?
    if (defined($FusedRingSetsBondsVisitedMapRef) && $RingBond->GetNumOfRings() == 2) {
      $RingBondID = $RingBond->GetID();
      if (exists $FusedRingSetsBondsVisitedMapRef->{$RingBondID}) {
	next RINGBOND;
      }
      $FusedRingSetsBondsVisitedMapRef->{$RingBondID} = $RingBondID;
    }
    $NumOfRingBondsProcessed++;

    # For first ring, previous ring bond corrresponds to last ring bond...
    $PreviousIndex = $Index ? ($Index -1) : $#{$RingBondsRef};
    $PreviousRingBond = $RingBondsRef->[$PreviousIndex];

    # Check for presence of alternate single/double bond configuration, and pesence of
    # hetero atoms with two single ring bonds along with any exocyclic double bonds...
    #
    BONDORDER: {
      # Is current ring double bond in an alternate single/double bond configuration?
      if ($BondOrder == 2) {
	if ($PreviousRingBond->GetBondOrder() != 1) {
	  return 0;
	}
	$NumOfConjugatedDoubleBonds += 1;
	last BONDORDER;
      }

      # Is current ring bond order correspond to an explicit aromatic bond?
      if ($BondOrder == 1.5) {
	if ($PreviousRingBond->GetBondOrder() != 1.5) {
	  return 0;
	}
	$NumOfExplicitAromaticBonds += 1;
	last BONDORDER;
      }

      # Check for potential hetero atoms involved in two single ring bonds along
      # with any terminal exocyclic bonds...
      if ($BondOrder == 1) {
	if ($PreviousRingBond->GetBondOrder() != 1) {
	  # Part of a conjugated system...
	  last BONDORDER;
	}

	# Identify any exocylic bonds on rings atoms...
	if (!defined $ExocyclicDoubleBondsDataMapRef) {
	  $ExocyclicDoubleBondsDataMapRef = $This->_IdentifyRingAtomsInvolvedInExocyclicDoubleBonds($RingAtomsRef, $RingBondsRef, $FusedRingSetsBondsMapRef);
	}

	# Is current ring atom part of an allowed exocyclic terminal bond?
	if (!$This->_CheckPotentialAromaticRingAtomForExocylicDoubleBonds($RingAtom, $AromaticityModelDataRef, $ExocyclicDoubleBondsDataMapRef)) {
	  return 0;
	}

	# Is it allowed to have any formal charge?
	if (!$This->_CheckPotentialAromaticRingAtomForFormalCharge($RingAtom, $AromaticityModelDataRef)) {
	  return 0;
	}

	# It it an allowed hetero ring atom or a carbon atom?
	if (!$This->_CheckPotentialAromaticRingAtomForAllowedHeteroAtoms($RingAtom, $AromaticityModelDataRef)) {
	  return 0;
	}

	$RingAtomID = $RingAtom->GetID();
	$ElectronPairContributionProcessedMap{$RingAtomID} = $RingAtomID;

	# Is it able to donate a pair for electrons towards pi electron delocalization?
	if ($RingAtom->GetValenceFreeElectrons($ExcludeFreeRadicalElectrons) >= 2) {
	  # Possibilites:
	  #   . Hetero atom with or without formal charge and an available electron pair
	  #   . Carbon atom with -ve formal charge and with an available electron pair
	  #
	  $NumOfRingAtomElectronPairs += 1;
	}
	else {
	  # Is ring atom involved in two single bonds without any electron pair allowed?
	  if (!$This->_AllowRingAtomInTwoSingleBondsWithoutElectronPair($RingAtom, $RingBond, $PreviousRingBond, $ExocyclicDoubleBondsDataMapRef, $FusedRingBondsMapRef)) {
	    return 0;
	  }
	}
	last BONDORDER;
      }

      # Any other type of ring atom/bond is not allowed to contribute towards pi electron count
      # and caused loss of aromaticity...
      return 0;
    }
  }

  # Check for any electron pair contributions towards pi electron delocalization due to
  # -ve formal charge on ring atoms which haven't been already processed and part of
  # conjugated single/double bond system...
  #
  $NumOfRingAtomElectronPairs += $This->_GetElectronPairsContributionFromConjugatedRingAtoms($RingAtomsRef, $RingBondsRef, $ExcludeFreeRadicalElectrons, $AromaticityModelDataRef, \%ElectronPairContributionProcessedMap);

  # Setup pi electron count available for delocalization...
  COUNT: {
    if ($NumOfExplicitAromaticBonds == $NumOfRingBondsProcessed) {
      # Each aromatic bond contribute one electron towards pi electron delocalization...
      $NumOfPiElectrons = $NumOfExplicitAromaticBonds;
      last COUNT;
    }

    # Each conjugated double bond contribute two electrons towards pi electron delocalization...
    $NumOfPiElectrons = 2*$NumOfConjugatedDoubleBonds + 2*$NumOfRingAtomElectronPairs;
  }

  return $NumOfPiElectrons;
}

# Check ring atoms for their potential participation in aromatic systems..
#
sub _CheckRingAtomsForPotentialAromaticity {
  my($This, $RingAtomsRef, $RingBondsRef, $AromaticityModelDataRef) = @_;
  my($Index, $RingBond, $RingAtom);

  # Check availability of ring atoms and bonds...
  if (!(defined($RingAtomsRef) && @{$RingBondsRef})) {
    return 0;
  }

  # Is there any minimum ring size limit?
  if ($AromaticityModelDataRef->{MinimumRingSize}) {
    if (@{$RingAtomsRef} < $AromaticityModelDataRef->{MinimumRingSize}) {
      return 0;
    }
  }

  # Make sure ring bond order is not greater than 2 and ring atom is not connected to more
  # than 3 other atoms to eliminate any non sp2 carbon atoms and still allow for hetero atoms
  # to contrbute towards electron delocalization...
  #
  for $Index (0 .. $#{$RingBondsRef}) {
    $RingBond = $RingBondsRef->[$Index];
    $RingAtom = $RingAtomsRef->[$Index];

    if (($RingBond->GetBondOrder() > 2) || ($RingAtom->GetNumOfBonds() + $RingAtom->GetNumOfMissingHydrogens()) > 3) {
      return 0;
    }
  }

  return 1;
}

# Identify any exocylic double bonds on ring atoms...
#
sub _IdentifyRingAtomsInvolvedInExocyclicDoubleBonds {
  my($This, $RingAtomsRef, $RingBondsRef, $FusedRingSetsBondsMapRef) = @_;
  my($Index, $RingAtom, $RingBond, $RingAtomID, $Bond, $BondID, $BondedAtom, $RingBondsMapRef, %RingBondsMap, %ExocyclicDoubleBondsDataMap);

  # Setup a ring bond map to process exocyclic bonds...
  $RingBondsMapRef = undef;
  %RingBondsMap = ();

  if (defined $FusedRingSetsBondsMapRef) {
    $RingBondsMapRef = $FusedRingSetsBondsMapRef;
  }
  else {
    for $BondID (map { $_->GetID() } @{$RingBondsRef}) {
      $RingBondsMap{$BondID} = $BondID;
    }
    $RingBondsMapRef = \%RingBondsMap;
  }

  # Intialize exocyclic terminal double bond data...
  %ExocyclicDoubleBondsDataMap = ();
  %{$ExocyclicDoubleBondsDataMap{RingAtomID}} = ();

  for $Index (0 .. $#{$RingBondsRef}) {
    $RingBond = $RingBondsRef->[$Index];
    $RingAtom = $RingAtomsRef->[$Index];

    $RingAtomID = $RingAtom->GetID();

    BOND: for $Bond ($RingAtom->GetBonds()) {
      if ($Bond->GetBondOrder != 2) {
	next BOND;
      }

      # Is it part of ring or ring system under consideration?
      if (exists $RingBondsMapRef->{$Bond->GetID()}) {
	next BOND;
      }

      # Is bonded atom in a ring or a non-terminal atom?
      $BondedAtom = $Bond->GetBondedAtom($RingAtom);
      if ($BondedAtom->IsInRing() || !$BondedAtom->IsTerminal() ) {
	next BOND;
      }

      # Track exocyclic terminal double bond information...
      if (!exists $ExocyclicDoubleBondsDataMap{RingAtomID}{$RingAtomID}) {
	@{$ExocyclicDoubleBondsDataMap{RingAtomID}{$RingAtomID}} = ();
      }
      push @{$ExocyclicDoubleBondsDataMap{RingAtomID}{$RingAtomID}}, $BondedAtom;
    }
  }

  return \%ExocyclicDoubleBondsDataMap;
}

# Check to see whether ring atoms are allowed to participate in exocyclic terminal double
# bonds...
#
sub _CheckPotentialAromaticRingAtomForExocylicDoubleBonds {
  my($This, $RingAtom, $AromaticityModelDataRef, $ExocyclicDoubleBondsDataMapRef) = @_;
  my($RingAtomID, $ExocyclicTerminalAtom, $RingAtomElectronegativity, $TerminalAtomElectronagativity);

  $RingAtomID = $RingAtom->GetID();

  # Is it part of an exocyclic terminal double bond?
  if (!exists $ExocyclicDoubleBondsDataMapRef->{RingAtomID}{$RingAtomID}) {
    return 1;
  }

  # Are exocyclic terminal double bonds allowed?
  if (!$AromaticityModelDataRef->{AllowExocyclicDoubleBonds}) {
    return 0;
  }

  # Are there multiple exocyclic double bonds?
  if (@{$ExocyclicDoubleBondsDataMapRef->{RingAtomID}{$RingAtomID}} > 1) {
    return 0;
  }
  ($ExocyclicTerminalAtom) = @{$ExocyclicDoubleBondsDataMapRef->{RingAtomID}{$RingAtomID}};

  # Are homo nuclear exocyclic terminal double bonds allowed?
  if (!$AromaticityModelDataRef->{AllowHomoNuclearExocyclicDoubleBonds}) {
    if ($RingAtom->GetAtomicNumber() == $ExocyclicTerminalAtom->GetAtomicNumber()) {
      return 0;
    }
  }

  # Are ring atoms with higher electronegativity allowed in exocyclic double bonds?
  if (!$AromaticityModelDataRef->{AllowElectronegativeRingAtomExocyclicDoubleBonds}) {
    $RingAtomElectronegativity = PeriodicTable::GetElementPaulingElectronegativity($RingAtom->GetAtomicNumber());
    $TerminalAtomElectronagativity = PeriodicTable::GetElementPaulingElectronegativity($ExocyclicTerminalAtom->GetAtomicNumber());

    if ($RingAtomElectronegativity && $TerminalAtomElectronagativity) {
      if ($RingAtomElectronegativity > $TerminalAtomElectronagativity) {
	return 0;
      }
    }
  }

  return 1;
}

#
# Check for any formal charge participation into electron delocalization...
#
sub _CheckPotentialAromaticRingAtomForFormalCharge {
  my($This, $RingAtom, $AromaticityModelDataRef) = @_;
  my($FormalCharge);

  # Does atom has any formal charge?
  $FormalCharge = $RingAtom->GetFormalCharge();
  if (!$FormalCharge) {
    return 1;
  }

  # Are ring atoms with formal charge allowed to participate in electron delocalization?
  if (!$AromaticityModelDataRef->{AllowRingAtomFormalCharge}) {
    return 0;
  }

  # Are hetero ring atoms with formal charge allowed to participate in electron delocalization?
  if (!$RingAtom->IsCarbon()) {
    if (!$AromaticityModelDataRef->{AllowHeteroRingAtomFormalCharge}) {
      return 0;
    }
  }

  return 1;
}

#
# Check ring atoms for allowed hetero atoms...
#
sub _CheckPotentialAromaticRingAtomForAllowedHeteroAtoms {
  my($This, $RingAtom, $AromaticityModelDataRef) = @_;
  my($RingAtomSymbol);

  # Is it a Carbon atom?
  if ($RingAtom->IsCarbon()) {
    return 1;
  }

  # Are heteroatoms allowed?
  if (!$AromaticityModelDataRef->{AllowHeteroRingAtoms}) {
    return 0;
  }

  # Is it an allowed hetero atom?
  $RingAtomSymbol = $RingAtom->GetAtomSymbol();
  if (!exists $AromaticityModelDataRef->{HeteroRingAtomsListMapRef}->{$RingAtomSymbol}) {
    return 0;
  }

  return 1;
}

# Check for any electron pair contributions toward pi electron delocalization due to
# -ve formal charge on ring atoms which haven't been already processed and part of
# conjugated single/double bond system...
#
sub _GetElectronPairsContributionFromConjugatedRingAtoms {
  my($This, $RingAtomsRef, $RingBondsRef, $ExcludeFreeRadicalElectrons, $AromaticityModelDataRef, $ElectronPairContributionProcessedMapRef) = @_;
  my($Index, $RingBond, $RingAtom, $NumOfRingAtomElectronPairs, $RingAtomID);

  # Is formal charge allowed on ring atoms?
  if (!$AromaticityModelDataRef->{AllowRingAtomFormalCharge}) {
    return 0;
  }

  $NumOfRingAtomElectronPairs = 0;

  # Process ring atoms...
  RINGBOND: for $Index (0 .. $#{$RingBondsRef}) {
    $RingBond = $RingBondsRef->[$Index];
    $RingAtom = $RingAtomsRef->[$Index];
    $RingAtomID = $RingAtom->GetID();

    # Is is already processed?
    if (exists $ElectronPairContributionProcessedMapRef->{$RingAtomID}) {
      next RINGBOND;
    }
    $ElectronPairContributionProcessedMapRef->{$RingAtomID} = $RingAtomID;

    # Is it allowed to have any formal charge?
    if (!$This->_CheckPotentialAromaticRingAtomForFormalCharge($RingAtom, $AromaticityModelDataRef)) {
      next RINGBOND;
    }

    # It it an allowed hetero ring atom or a carbon atom?
    if (!$This->_CheckPotentialAromaticRingAtomForAllowedHeteroAtoms($RingAtom, $AromaticityModelDataRef)) {
      next RINGBOND;
    }

    # It is an atom with -ve formal charge?
    if ($RingAtom->GetFormalCharge() >= 0) {
      next RINGBOND;
    }

    # Is it able to donate a pair for electrons towards pi electron delocalization?
    if ($RingAtom->GetValenceFreeElectrons($ExcludeFreeRadicalElectrons) < 2) {
      next RINGBOND;
    }
    $NumOfRingAtomElectronPairs += 1;
  }

  return $NumOfRingAtomElectronPairs;
}

# Check for ring atoms involved in two single ring bonds without any available electron
# pair which are allowed to participate in aromatic system, after all other checks
# corresponding to specified aromaticity models have already been performed...
#
sub _AllowRingAtomInTwoSingleBondsWithoutElectronPair {
  my($This, $RingAtom, $RingBond, $PreviousRingBond, $ExocyclicDoubleBondsDataMapRef, $FusedRingBondsMapRef) = @_;

  ALLOWRINGATOM: {
    if (exists $ExocyclicDoubleBondsDataMapRef->{RingAtomID}{$RingAtom->GetID()}) {
      # Ring atom in an exocylic terminal double bond without any available electron pair...
      last ALLOWRINGATOM;
    }

    if ($RingAtom->GetFormalCharge() > 0) {
      # Ring atom with positive formal charge without any available electron pair...
      last ALLOWRINGATOM;
    }

    if (defined $FusedRingBondsMapRef && (exists $FusedRingBondsMapRef->{$RingBond->GetID()} || exists $FusedRingBondsMapRef->{$PreviousRingBond->GetID()})) {
      # Ring atom involved in fused ring bond, which might end up being part of a conjugated
      # system in another fused ring...
      last ALLOWRINGATOM;
    }

    # Ring atom in any other environment is not allowed...
    return 0;
  }

  return 1;
}

# Do pi electrons satify huckel's rule: Number of pi electrons correspond to 4n + 2 where
# n is a positive integer...
#
sub _DoPiElectronSatifyHuckelsRule {
  my($This, $NumOfPiElectrons) = @_;

  $NumOfPiElectrons = $NumOfPiElectrons - 2;

  return ($NumOfPiElectrons > 0) ? (($NumOfPiElectrons % 4) ? 0 : 1) : 0;
}

# Delete aromatic property for all atoms and bonds...
#
sub _DeleteAtomsAndBondsAromaticity {
  my($This) = @_;
  my($Atom, $Bond);

  for $Atom ($This->GetAtoms()) {
    $Atom->DeleteAromatic();
  }
  for $Bond ($This->GetBonds()) {
    $Bond->DeleteAromatic();
  }
  return $This;
}

# Kekulize marked ring and non-ring aromatic atoms in a molecule...
#
sub KekulizeAromaticAtoms {
  my($This) = @_;

  if (!$This->_KekulizeAromaticAtomsInRings()) {
    return 0;
  }

  if (!$This->_KekulizeAromaticAtomsNotInRings()) {
    return 0;
  }

  return 1;
}

# Kekulize marked aromatic atoms in rings and fused ring sets...
#
sub _KekulizeAromaticAtomsInRings {
  my($This) = @_;

  if (!$This->HasRings()) {
    # Nothing to do...
    return 1;
  }

  if (!$This->HasAromaticAtomsInRings()) {
    # Nothing to do...
    return 1;
  }

  # Identify fully aromatic fused and individual rings along with any partially aromatic ring components
  # using marked aromatic atoms in a molecule and kekulize them as individual stes...
  #
  my($AromaticFusedRingSetsRef, $AromaticRingsRef, $PartiallyAromaticRingComponentsRef) = (undef) x 3;
  if ($This->HasFusedRings()) {
    ($AromaticFusedRingSetsRef, $AromaticRingsRef, $PartiallyAromaticRingComponentsRef) = $This->_GetFusedAndNonFusedRingsContainingAromaticAtoms();
  }
  else {
    ($AromaticRingsRef, $PartiallyAromaticRingComponentsRef) = $This->_GetIndividualRingsContainingAromaticAtoms();
  }

  return $This->_KekulizeCompleteAndPartialAromaticRings($AromaticFusedRingSetsRef, $AromaticRingsRef, $PartiallyAromaticRingComponentsRef);
}

# Identify fully aromatic fused and individual rings along with any partially aromatic ring components
# using marked aromatic atoms in a molecule...
#
sub _GetFusedAndNonFusedRingsContainingAromaticAtoms {
  my($This)  = @_;
  my($Index, $SetAtomsCount, $SetAromaticAtomsCount, $FusedRingSetRef, $FusedRingSetsRef, $NonFusedRingsRef, $IndividualRingsRef, $RingAtomsRef, $RingAtomsCount, $AromaticAtomsCount, $RingAtom, $NonFusedFullyAromaticRingsRef, $NonFusedPartiallyAromaticRingComponentsRef, $PartiallyAromaticRingComponentsRef, @FullyAromaticFusedRingSets, @PotentialFullyAromaticRings, @FullyAromaticRings, @PotentialPartiallyAromaticRings, @PartiallyAromaticRingComponents);

  @FullyAromaticFusedRingSets = ();

  @PotentialFullyAromaticRings = ();
  @FullyAromaticRings = ();

  @PotentialPartiallyAromaticRings = ();
  @PartiallyAromaticRingComponents = ();

  ($FusedRingSetsRef, $NonFusedRingsRef) = $This->GetFusedAndNonFusedRings();

  # Go over fused ring sets...
  RINGSET: for $Index (0 .. $#{$FusedRingSetsRef}) {
    $FusedRingSetRef = $FusedRingSetsRef->[$Index];

    $SetAtomsCount = 0;
    $SetAromaticAtomsCount = 0;

    for $RingAtomsRef (@{$FusedRingSetRef}) {
      $SetAtomsCount += scalar @{$RingAtomsRef};

      for $RingAtom (@{$RingAtomsRef}) {
	if ($RingAtom->IsAromatic()) {
	  $SetAromaticAtomsCount += 1;
	}
      }
    }

    if (!($SetAtomsCount && $SetAromaticAtomsCount)) {
      next RINGSET;
    }

    if ($SetAromaticAtomsCount == $SetAtomsCount) {
      push @FullyAromaticFusedRingSets, $FusedRingSetRef;
    }
    else {
      # Identify any individual rings in partial aromatic fused ring sets which might be
      # fully or partially aromatic...
      #
      RING: for $RingAtomsRef (@{$FusedRingSetRef}) {
	$RingAtomsCount = scalar @{$RingAtomsRef};
	$AromaticAtomsCount = 0;

	RINGATOM: for $RingAtom (@{$RingAtomsRef}) {
	  if (!$RingAtom->IsAromatic()) {
	    next RINGATOM;
	  }
	  $AromaticAtomsCount += 1;
	}

	if (!($RingAtomsCount && $AromaticAtomsCount)) {
	  next RING;
	}

	if ($RingAtomsCount == $AromaticAtomsCount) {
	  push @PotentialFullyAromaticRings, $RingAtomsRef;
	}
	else {
	  #  Track partially aromatic rings in an different list before removing them for
	  #  any overlap with other rings and then add to fully aromatic rings...
	  push @PotentialPartiallyAromaticRings, $RingAtomsRef;
	}
      }
    }
  }

  if (@PotentialFullyAromaticRings > 1) {
    # Get any fully aromatic fused ring subsets from potentially fully aromatic rings...
    my($FullyAromaticFusedRingSetsRefs, $FullyAromaticNonFusedRingsRef);
    ($FullyAromaticFusedRingSetsRefs, $FullyAromaticNonFusedRingsRef) = $This->_GetFullyAromaticFusedAndNonFusedRingsInFusedSubset(\@PotentialFullyAromaticRings);

    if (@{$FullyAromaticFusedRingSetsRefs}) {
      push @FullyAromaticFusedRingSets, @{$FullyAromaticFusedRingSetsRefs};
    }
    if (@{$FullyAromaticNonFusedRingsRef}) {
      push @FullyAromaticRings, @{$FullyAromaticNonFusedRingsRef};
    }
  }
  else {
    push @FullyAromaticRings, @PotentialFullyAromaticRings;
  }

  # Go over partial aromatic ring components...
  if (@PotentialPartiallyAromaticRings) {
    $PartiallyAromaticRingComponentsRef = $This->_GetPartiallyAromaticRingComponents(\@PotentialPartiallyAromaticRings, \@PotentialFullyAromaticRings);
    if (@{$PartiallyAromaticRingComponentsRef}) {
      push @PartiallyAromaticRingComponents, @{$PartiallyAromaticRingComponentsRef};
    }
  }

  # Go over non-fused rings...
  if (@{$NonFusedRingsRef}) {
    ($NonFusedFullyAromaticRingsRef, $NonFusedPartiallyAromaticRingComponentsRef) = $This->_GetRingsContainingAromaticAtoms(@{$NonFusedRingsRef});

    if (@{$NonFusedFullyAromaticRingsRef}) {
      push @FullyAromaticRings, @{$NonFusedFullyAromaticRingsRef};
    }
    if (@{$NonFusedPartiallyAromaticRingComponentsRef}) {
      push @PartiallyAromaticRingComponents, @{$NonFusedPartiallyAromaticRingComponentsRef};
    }
  }

  return (\@FullyAromaticFusedRingSets, \@FullyAromaticRings, \@PartiallyAromaticRingComponents);
}

# Identify fully aromatic fused sets and non-fused rings in potentially fully aromatic
# rings in fused ring sets...
#
# Fully aromatic rings in fused ring sets might contain fully aromatic fused subsets. These
# fused subets need to be tracked and treated as fused sets.
#
# Note:
#   . Fused ring sets share at least one common bond, which could be used to identify
#     any multiple fully aromatic fused rings sets; absence of a shared ring bond implies
#     there are no fused ring sets.
#
#
sub _GetFullyAromaticFusedAndNonFusedRingsInFusedSubset {
  my($This, $PotentialFullyAromaticFusedRingsRef) = @_;
  my($RingIndex, $RingIndex1, $RingIndex2, $RingAtom, $RingAtomID, $RingIsFuesd, $RingIndicesGraph, $FusedRingSetIndicesRef, @RingIndices, @FusedRingPairIndices, @FusedRingSetIndicesRefs, @FullyAromaticFusedRingSets, @FullyAromaticRings, %RingIndexToAtomIDMap, %FullyAromaticFusedRingIndexMap);

  @FullyAromaticFusedRingSets = ();
  @FullyAromaticRings = ();

  # Setup a ring index map for ring atoms...
  #
  %RingIndexToAtomIDMap = ();
  for $RingIndex (0 .. $#{$PotentialFullyAromaticFusedRingsRef}) {
    %{$RingIndexToAtomIDMap{$RingIndex}} = ();
    for $RingAtom (@{$PotentialFullyAromaticFusedRingsRef->[$RingIndex]}) {
      $RingAtomID = $RingAtom->GetID();
      $RingIndexToAtomIDMap{$RingIndex}{$RingAtomID} = $RingAtomID;
    }
  }

  # Identify fused ring pairs...
  #
  @RingIndices = ();
  @FusedRingPairIndices = ();

  for $RingIndex1 (0 .. $#{$PotentialFullyAromaticFusedRingsRef}) {
    push @RingIndices, $RingIndex1;
    for $RingIndex2 (($RingIndex1 + 1) .. $#{$PotentialFullyAromaticFusedRingsRef}) {
      $RingIsFuesd = 0;
      RINGATOM: for $RingAtom (@{$PotentialFullyAromaticFusedRingsRef->[$RingIndex2]}) {
	$RingAtomID = $RingAtom->GetID();
	if (exists $RingIndexToAtomIDMap{$RingIndex1}{$RingAtomID}) {
	  $RingIsFuesd = 1;
	  last RINGATOM;
	}
      }
      if ($RingIsFuesd) {
	push @FusedRingPairIndices, ($RingIndex1, $RingIndex2);
      }
    }
  }

  if (!@FusedRingPairIndices) {
    # No fused ring subset out there...
    push @FullyAromaticRings, @{$PotentialFullyAromaticFusedRingsRef};

    return (\@FullyAromaticFusedRingSets, \@FullyAromaticRings);
  }

  # Identify fused ring sets...
  #
  $RingIndicesGraph = new Graph(@RingIndices);
  $RingIndicesGraph->AddEdges(@FusedRingPairIndices);
  @FusedRingSetIndicesRefs = $RingIndicesGraph->GetConnectedComponentsVertices();

  # Collect fully aromatic fused ring sets...
  #
  %FullyAromaticFusedRingIndexMap = ();
  for $FusedRingSetIndicesRef (@FusedRingSetIndicesRefs) {
    my(@FullyAromaticFusedRingSet) = ();
    for $RingIndex (@{$FusedRingSetIndicesRef}) {
      $FullyAromaticFusedRingIndexMap{$RingIndex} = $RingIndex;
      push @FullyAromaticFusedRingSet, $PotentialFullyAromaticFusedRingsRef->[$RingIndex];
    }
    if (@FullyAromaticFusedRingSet) {
      # Sort rings by size with in the fused ring set...
      @FullyAromaticFusedRingSet = sort { scalar @$a <=> scalar @$b } @FullyAromaticFusedRingSet;
      push @FullyAromaticFusedRingSets, \@FullyAromaticFusedRingSet;
    }
  }

  # Collect fully aromatic non-fused rings...
  #
  RINGINDEX: for $RingIndex (0 .. $#{$PotentialFullyAromaticFusedRingsRef}) {
    if (exists $FullyAromaticFusedRingIndexMap{$RingIndex}) {
      next RINGINDEX;
    }
    push @FullyAromaticRings, $PotentialFullyAromaticFusedRingsRef->[$RingIndex];
  }

  return (\@FullyAromaticFusedRingSets, \@FullyAromaticRings);
}

# Identify individual non-fused rings containing aromatic atoms...
#
sub _GetIndividualRingsContainingAromaticAtoms {
  my($This)  = @_;

  return $This->_GetRingsContainingAromaticAtoms($This->GetRings());
}

# Identify individual non-fused rings containing aromatic atoms...
#
sub _GetRingsContainingAromaticAtoms {
  my($This, @Rings)  = @_;
  my($RingAtom, $RingAtomsRef, $RingAtomsCount, $AromaticAtomsCount, $PartiallyAromaticRingComponentsRef, @FullyAromaticRings, @PartiallyAromaticRings);

  @FullyAromaticRings = ();
  @PartiallyAromaticRings = ();

  RING: for $RingAtomsRef (@Rings) {
    $RingAtomsCount = scalar @{$RingAtomsRef};
    $AromaticAtomsCount = 0;

    for $RingAtom (@{$RingAtomsRef}) {
      if ($RingAtom->IsAromatic()) {
	$AromaticAtomsCount += 1;
      }
    }

    if (!($AromaticAtomsCount && $RingAtomsCount)) {
      next RING;
    }

    if ($AromaticAtomsCount == $RingAtomsCount) {
      push @FullyAromaticRings, $RingAtomsRef;
    }
    else {
      push @PartiallyAromaticRings, $RingAtomsRef;
    }
  }

  $PartiallyAromaticRingComponentsRef = $This->_GetPartiallyAromaticRingComponents(\@PartiallyAromaticRings);

  return (\@FullyAromaticRings, $PartiallyAromaticRingComponentsRef);
}

# Get connected aromatic components with in partially aromatic rings...
#
sub _GetPartiallyAromaticRingComponents {
  my($This, $PotentialPartiallyAromaticRingsRef, $FullyAromaticRingsRef) = @_;
  my($RingAtomsRef, $RingAtom, $RingAtomID, $Index, @PartiallyAromaticRingComponents, %FullyAromaticRingAtomsMap);

  @PartiallyAromaticRingComponents = ();

  # Setup a map for atoms involve in fully aromatic rings to remove remove partial rings
  # containing only those atoms which are already part of some other fully aromatic ring
  # in fused ring scenarios or some other partially aromatic ring...
  #
  %FullyAromaticRingAtomsMap = ();
  if (defined $FullyAromaticRingsRef) {
    for $RingAtomsRef (@{$FullyAromaticRingsRef}) {
      for $RingAtom (@{$RingAtomsRef}) {
	$RingAtomID = $RingAtom->GetID();
	$FullyAromaticRingAtomsMap{$RingAtomID} = $RingAtomID;
      }
    }
  }

  # . Identify any connected components with in each partially aromatic ring.
  # . Use ring atom indices to figure out connnected components in rings: All ring atoms
  #   in a connected component have sequential indices and a difference by more than
  #   1 indicates a new component in the list.
  #
  RING: for $RingAtomsRef (@{$PotentialPartiallyAromaticRingsRef}) {
    my(@AromaticRingAtoms, @AromaticRingAtomsIndices);

    @AromaticRingAtoms = ();
    @AromaticRingAtomsIndices = ();

    RINGATOM: for $Index (0 .. $#{$RingAtomsRef}) {
      $RingAtom = $RingAtomsRef->[$Index];
      $RingAtomID = $RingAtom->GetID();

      if (defined $FullyAromaticRingsRef && exists $FullyAromaticRingAtomsMap{$RingAtomID}) {
	next RINGATOM;
      }
      if (!$RingAtom->IsAromatic()) {
	next RINGATOM;
      }
      push @AromaticRingAtoms, $RingAtom;
      push @AromaticRingAtomsIndices, $Index;

    }
    if (!@AromaticRingAtoms) {
      next RING;
    }

    # Start off with a new connected component...
    #
    my($ComponentNum);
    $ComponentNum = scalar @PartiallyAromaticRingComponents;
    @{$PartiallyAromaticRingComponents[$ComponentNum]} = ();

    $Index = 0;
    push @{$PartiallyAromaticRingComponents[$ComponentNum]}, $AromaticRingAtoms[$Index];

    for $Index (1 .. $#AromaticRingAtoms) {
      if (($AromaticRingAtomsIndices[$Index] - $AromaticRingAtomsIndices[$Index -1]) > 1) {
	# New connected component...
	$ComponentNum += 1;
	@{$PartiallyAromaticRingComponents[$ComponentNum]} = ();
      }
      push @{$PartiallyAromaticRingComponents[$ComponentNum]}, $AromaticRingAtoms[$Index];
    }
  }

  return (\@PartiallyAromaticRingComponents);
}

# Kekulize fully aromatic fused and individual rings along with any partially aromatic ring
# components...
#
sub _KekulizeCompleteAndPartialAromaticRings {
  my($This, $AromaticFusedRingSetsRef, $AromaticRingsRef, $PartiallyAromaticRingComponentsRef) = @_;
  my($Status, $ConnectedPathsAtomsSetsRef, $ConnectedPathsBondsSetsRef, $ConnectdPathsSetsTypesRef, $PathSetIndex, $PathAtom, $AtomID, $BondID, $PathBondsRef, $DeleteAtomsAromaticity, $DeleteBondsAromaticity, %PathAtomsProcessingStatusMap, %PathBondsProcessingStatusMap);

  ($ConnectedPathsAtomsSetsRef, $ConnectedPathsBondsSetsRef, $ConnectdPathsSetsTypesRef) = $This->_SetupCompleteAndPartialAromaticRingsForKekulizaiton($AromaticFusedRingSetsRef, $AromaticRingsRef, $PartiallyAromaticRingComponentsRef);

  if (!@{$ConnectedPathsAtomsSetsRef}) {
    # Nothing to do...
    return 1;
  }

  # Delete any aromaticity property set for non-ring bonds connected any two ring
  # aromatic atoms...
  #
  $This->_ProcessNonRingAromaticBondsBetweenAromaticRingAtoms();

  %PathAtomsProcessingStatusMap = ();
  %PathBondsProcessingStatusMap = ();

  $Status = 1;

  PATHSET: for $PathSetIndex (0 .. $#{$ConnectedPathsAtomsSetsRef}) {
    my($AtomProcessingStatusMapRef, $BondProcessingStatusMapRef, $PathSetProcessingStatusRef);

    ($AtomProcessingStatusMapRef, $BondProcessingStatusMapRef, $PathSetProcessingStatusRef) = (undef) x 3;

    if ($ConnectdPathsSetsTypesRef->[$PathSetIndex] =~ /^FusedAromatic$/i) {
      # Fused set of connected paths...
      #
      my($FusedConnectedPathAtomsSetRef, $FusedConnectedPathBondsSetRef);

      $FusedConnectedPathAtomsSetRef = $ConnectedPathsAtomsSetsRef->[$PathSetIndex];
      $FusedConnectedPathBondsSetRef = $ConnectedPathsBondsSetsRef->[$PathSetIndex];

      # Prepare for kekulization...
      ($AtomProcessingStatusMapRef, $BondProcessingStatusMapRef, $PathSetProcessingStatusRef) = $This->_SetupConnectedPathSetsForKekulization($FusedConnectedPathAtomsSetRef, $FusedConnectedPathBondsSetRef);

      # Perform kekulization starting with the first path set...
      $PathSetProcessingStatusRef->[0] = 'Processed';
      if (!$This->_KekulizeConnectedPathSets($FusedConnectedPathAtomsSetRef->[0],  $FusedConnectedPathBondsSetRef->[0], $AtomProcessingStatusMapRef, $BondProcessingStatusMapRef, $FusedConnectedPathAtomsSetRef, $FusedConnectedPathBondsSetRef, $PathSetProcessingStatusRef)) {
	# Kekulization failed for the current fused paths set...
	$Status = 0;
      }
    }
    else {
      # An individual connected path...
      #
      my(@ConnectedPathAtomsSet, @ConnectedPathBondsSet);

      @ConnectedPathAtomsSet = ($ConnectedPathsAtomsSetsRef->[$PathSetIndex]);
      @ConnectedPathBondsSet = ($ConnectedPathsBondsSetsRef->[$PathSetIndex]);

      # Prepare for kekulization...
      ($AtomProcessingStatusMapRef, $BondProcessingStatusMapRef) = $This->_SetupConnectedPathSetsForKekulization(\@ConnectedPathAtomsSet, \@ConnectedPathBondsSet);

      # Perform kekulization...
      if (!$This->_KekulizeConnectedPathSets($ConnectedPathsAtomsSetsRef->[$PathSetIndex], $ConnectedPathsBondsSetsRef->[$PathSetIndex], $AtomProcessingStatusMapRef, $BondProcessingStatusMapRef)) {
	# Kekulization failed for the current path...
	$Status = 0;
      }
    }

    # Did kekulization succeed for the current path or path set?
    if (!$Status) {
      last PATHSET;
    }

    # Track atom and bond processing state for final assignment after kekulization
    # is successfully completed for all the paths and fused path sets...
    #
    for $AtomID (keys %{$AtomProcessingStatusMapRef}) {
      $PathAtomsProcessingStatusMap{$AtomID} = $AtomProcessingStatusMapRef->{$AtomID};
    }

    for $BondID (keys %{$BondProcessingStatusMapRef}) {
      $PathBondsProcessingStatusMap{$BondID} = $BondProcessingStatusMapRef->{$BondID};
    }
  }

  if (!$Status) {
    carp "Warning: ${ClassName}->_KekulizeCompleteAndPartialAromaticRings: Couldn't perform kekulization for marked ring aromatic atoms...";
    return 0;
  }

  # Use PathAtomsProcessingStatusMap and PathBondsProcessingStatusMap to set
  # single/double bonds in the molecule after successful kekulization along with modification of
  # any aromatic flags...

  for $PathSetIndex (0 .. $#{$ConnectedPathsAtomsSetsRef}) {
    $DeleteAtomsAromaticity = 0; $DeleteBondsAromaticity = 0;

    if ($ConnectdPathsSetsTypesRef->[$PathSetIndex] =~ /^FusedAromatic$/i) {
      for $PathBondsRef (@{$ConnectedPathsBondsSetsRef->[$PathSetIndex]}) {
	$This->_ProcessBondOrdersAssignedDuringSuccessfulKekulization($PathBondsRef, \%PathBondsProcessingStatusMap, $DeleteBondsAromaticity);
      }
    }
    else {
      if ($ConnectdPathsSetsTypesRef->[$PathSetIndex] =~ /^PartiallyAromatic$/i ) {
	$DeleteBondsAromaticity = 1; $DeleteAtomsAromaticity = 1;
      }

      if ($DeleteAtomsAromaticity) {
	for $PathAtom (@{$ConnectedPathsAtomsSetsRef->[$PathSetIndex]}) {
	  $PathAtom->DeleteAromatic();
	}
      }

      $This->_ProcessBondOrdersAssignedDuringSuccessfulKekulization($ConnectedPathsBondsSetsRef->[$PathSetIndex], \%PathBondsProcessingStatusMap, $DeleteBondsAromaticity);
    }
  }

  return 1;
}

# Look for any aromatic bonds outside the rings between two ring aromatic atoms
# and turn them into single non-aromatic bonds before kekulization; otherwise, kekulization
# fails.
#
# Note:
#   . Two atoms marked as aromatic atoms in two different rings, such as two rings
#     connected through a single bond, are still aromatic, but the bond is outside
#     the ring and shouldn't be marked as aromatic. It should be set to single bond without
#     any aromatic property for kekulization to succeed.
#
#     For example, the molecule  generated by SMILES parser for biphenyl SMILES string
#     "c1ccccc1c2ccccc2" sets up an aromatic bond between the two phenyl rings, as
#     it's connected to two aromatic atoms.
#
sub _ProcessNonRingAromaticBondsBetweenAromaticRingAtoms {
  my($This) = @_;
  my($Bond, $Atom1, $Atom2);

  BOND: for $Bond ($This->GetBonds()) {
    if (!($Bond->IsAromatic() && $Bond->IsNotInRing())) {
      next BOND;
    }

    ($Atom1, $Atom2) = $Bond->GetAtoms();
    if (!($Atom1->IsAromatic() && $Atom2->IsAromatic() && $Atom1->IsInRing() && $Atom2->IsInRing())) {
      next BOND;
    }

    $Bond->SetBondOrder(1);
    $Bond->DeleteAromatic();
  }

  return $This;
}

# Setup completelty aromatic fused and individual rings along with partially aromatic ring
# components as sets of connected paths...
#
sub _SetupCompleteAndPartialAromaticRingsForKekulizaiton {
  my($This, $AromaticFusedRingSetsRef, $AromaticRingsRef, $PartiallyAromaticRingComponentsRef) = @_;
  my(@ConnectedPathsSets, @ConnectedPathsBondsSets, @ConnectdPathsSetsTypes);

  @ConnectedPathsSets = ();
  @ConnectedPathsBondsSets = ();
  @ConnectdPathsSetsTypes = ();

  # Setup atoms and bonds for connected paths in fused aromatic ring sets...
  #
  if (defined $AromaticFusedRingSetsRef && @{$AromaticFusedRingSetsRef}) {
    my($RingSetIndex);

    push @ConnectdPathsSetsTypes, ('FusedAromatic') x scalar @{$AromaticFusedRingSetsRef};
    push @ConnectedPathsSets, @{$AromaticFusedRingSetsRef};

    for $RingSetIndex (0 .. $#{$AromaticFusedRingSetsRef}) {
      my(@AromaticFusedRingBondsSet);

      # Get ring bonds for each ring set...
      #
      @AromaticFusedRingBondsSet = $This->GetRingBondsFromRings(@{$AromaticFusedRingSetsRef->[$RingSetIndex]});
      push @ConnectedPathsBondsSets, \@AromaticFusedRingBondsSet;
    }
  }

  # Set up atoms and bonds for connected paths in aromatic rings...
  #
  if (defined $AromaticRingsRef && @{$AromaticRingsRef}) {
    my(@AromaticRingBondsSets);

    push @ConnectdPathsSetsTypes, ('Aromatic') x scalar @{$AromaticRingsRef};
    push @ConnectedPathsSets, @{$AromaticRingsRef};

    # Get ring bonds for each ring...
    @AromaticRingBondsSets = $This->GetRingBondsFromRings(@{$AromaticRingsRef});
    push @ConnectedPathsBondsSets, @AromaticRingBondsSets;
  }

  # Set up atoms and bonds for connected paths in partially aromatic rings...
  #
  if (defined $PartiallyAromaticRingComponentsRef && @{$PartiallyAromaticRingComponentsRef}) {
    my($ComponentIndex);

    push @ConnectedPathsSets, @{$PartiallyAromaticRingComponentsRef};
    push @ConnectdPathsSetsTypes, ('PartiallyAromatic') x scalar @{$PartiallyAromaticRingComponentsRef};

    for $ComponentIndex (0 .. $#{$PartiallyAromaticRingComponentsRef}) {
      my(@ComponentBonds);
      @ComponentBonds = $This->_GetPathBonds($This->_GetAtomsIDsFromAtoms(@{$PartiallyAromaticRingComponentsRef->[$ComponentIndex]}));
      push @ConnectedPathsBondsSets, \@ComponentBonds;
    }
  }

  return (\@ConnectedPathsSets, \@ConnectedPathsBondsSets, \@ConnectdPathsSetsTypes);
}

# Process non-ring connected atoms which are marked aromatic and set connected
# bonds as alternate single/double bonds...
#
# Notes:
#   . Atom and bond aromaticity is deleted during kekulization of non-ring atoms.
#
sub _KekulizeAromaticAtomsNotInRings {
  my($This) = @_;
  my($Status, $PathIndex, $PathAtom, $PathAtomID, $PathBondID, $ConnectedPathsAtomsRef, $ConnectedPathsBondsRef, $DeleteAtomsAromaticity, $DeleteBondsAromaticity, %PathAtomsProcessingStatusMap, %PathBondsProcessingStatusMap);

  if (!$This->HasAromaticAtomsNotInRings()) {
    # Nothing to do...
    return 1;
  }

  # Identify paths for connected components containing non-ring aromatic atoms...
  ($ConnectedPathsAtomsRef, $ConnectedPathsBondsRef) = $This->_GetConnectedComponentsPathsForNonRingAromaticAtoms();

  if (!@{$ConnectedPathsAtomsRef}) {
    carp "Warning: ${ClassName}->_KekulizeAromaticAtomsNotInRings: Couldn't perform kekulization for marked non-ring aromatic atoms...";
    return 0;
  }

  %PathAtomsProcessingStatusMap = ();
  %PathBondsProcessingStatusMap = ();

  $Status = 1;

  PATH: for $PathIndex (0 .. $#{$ConnectedPathsAtomsRef}) {
    my($AtomProcessingStatusMapRef, $BondProcessingStatusMapRef, @ConnectedPathAtomsSet, @ConnectedPathBondsSet);

    ($AtomProcessingStatusMapRef, $BondProcessingStatusMapRef) = (undef) x 2;

    @ConnectedPathAtomsSet = ($ConnectedPathsAtomsRef->[$PathIndex]);
    @ConnectedPathBondsSet = ($ConnectedPathsBondsRef->[$PathIndex]);

    # Prepare for kekulization...
    ($AtomProcessingStatusMapRef, $BondProcessingStatusMapRef) = $This->_SetupConnectedPathSetsForKekulization(\@ConnectedPathAtomsSet, \@ConnectedPathBondsSet);

    # Perform kekulization...
    if (!$This->_KekulizeConnectedPathSets($ConnectedPathsAtomsRef->[$PathIndex], $ConnectedPathsBondsRef->[$PathIndex], $AtomProcessingStatusMapRef, $BondProcessingStatusMapRef)) {
      # Kekulization failed for the current path...
      $Status = 0;
      last PATH;
    }

    # Track atom and bond processing state for final assignment after kekulization
    # is successfully completed for all the paths and fused path sets...
    #
    for $PathAtomID (keys %{$AtomProcessingStatusMapRef}) {
      $PathAtomsProcessingStatusMap{$PathAtomID} = $AtomProcessingStatusMapRef->{$PathAtomID};
    }

    for $PathBondID (keys %{$BondProcessingStatusMapRef}) {
      $PathBondsProcessingStatusMap{$PathBondID} = $BondProcessingStatusMapRef->{$PathBondID};
    }
  }

  if (!$Status) {
    carp "Warning: ${ClassName}->_KekulizeAromaticAtomsNotInRings: Couldn't perform kekulization for marked non-ring aromatic atoms...";
    return 0;
  }

  $DeleteAtomsAromaticity = 1; $DeleteBondsAromaticity = 1;
  for $PathIndex (0 .. $#{$ConnectedPathsAtomsRef}) {
    if ($DeleteAtomsAromaticity) {
      for $PathAtom (@{$ConnectedPathsAtomsRef->[$PathIndex]}) {
	$PathAtom->DeleteAromatic();
      }
    }
    $This->_ProcessBondOrdersAssignedDuringSuccessfulKekulization($ConnectedPathsBondsRef->[$PathIndex], \%PathBondsProcessingStatusMap, $DeleteBondsAromaticity);
  }

  return 1;
}

# Collect path atoms for connected components paths containing non-ring aromatic atoms...
#
sub _GetConnectedComponentsPathsForNonRingAromaticAtoms {
  my($This) = @_;
  my($ComponentRef, $AtomIDsRef, $AtomIDsMapRef, $ConnectedComponentsAtomIDsRef, $ConnectedComponentsAtomIDsMapRef, $ConnectedComponentsPathsAtomIDsRef, $ConnectedComponentsPathsAtomsRef, $ConnectedComponentsPathsBondsRef);

  # Retrieve information for marked aromatic atoms not in the rings...
  ($AtomIDsRef, $AtomIDsMapRef) = $This->_GetNonRingAromaticAtomIDs();

  # Identify connected components containing marked aromatic atoms not in the rings...
  ($ConnectedComponentsAtomIDsRef, $ConnectedComponentsAtomIDsMapRef) = $This->_GetConnectedComponentsForNonRingAromaticAtoms($AtomIDsRef);

  # Identify paths for connected components containing non-ring aromatic atoms...
  ($ConnectedComponentsPathsAtomsRef, $ConnectedComponentsPathsBondsRef) = $This->_GetConnectedComponentsPathsAtomsAndBondsForNonRingAromaticAtoms($AtomIDsMapRef, $ConnectedComponentsAtomIDsRef, $ConnectedComponentsAtomIDsMapRef);

  return ($ConnectedComponentsPathsAtomsRef, $ConnectedComponentsPathsBondsRef);
}

# Collect information for marked aromatic atoms not in the rings...
#
sub _GetNonRingAromaticAtomIDs {
  my($This) = @_;
  my($Atom, $AtomID, @AtomIDs, %AtomIDsMap);

  @AtomIDs = ();
  %AtomIDsMap = ();

  ATOM: for $Atom ($This->GetAtoms()) {
    if (!$Atom->IsAromatic()) {
      next ATOM;
    }
    if ($Atom->IsInRing()) {
      next ATOM;
    }
    $AtomID = $Atom->GetID();

    push @AtomIDs, $AtomID;
    $AtomIDsMap{$AtomID} = $Atom;
  }

  return (\@AtomIDs, \%AtomIDsMap);
}

# Retrieve connected non-ring atom components as a reference to an array of references
# containing atom IDs of connecnted components...
#
sub _GetConnectedComponentsForNonRingAromaticAtoms {
  my($This, $AtomIDsRef) = @_;
  my($Index, $AtomID, $AtomIDsGraph, @BondedAtomPairIDs, @ComponentsAtomIDsRefs, @ComponentsAtomIDsMapRefs);

  @ComponentsAtomIDsRefs = ();
  @ComponentsAtomIDsMapRefs = ();

  # Get bonded atom pair IDs...
  @BondedAtomPairIDs = $This->_GetBondedAtomPairAtomIDsFromAtomIDs(@{$AtomIDsRef});

  if (!@BondedAtomPairIDs) {
    return (\@ComponentsAtomIDsRefs, \@ComponentsAtomIDsMapRefs);
  }

  $AtomIDsGraph = new Graph(@{$AtomIDsRef});
  $AtomIDsGraph->AddEdges(@BondedAtomPairIDs);

  @ComponentsAtomIDsRefs = $AtomIDsGraph->GetConnectedComponentsVertices();

  # Setup atom IDs map for each component...
  for $Index (0 .. $#ComponentsAtomIDsRefs) {
    %{$ComponentsAtomIDsMapRefs[$Index]} = ();

    for $AtomID (@{$ComponentsAtomIDsRefs[$Index]}) {
      $ComponentsAtomIDsMapRefs[$Index]{$AtomID} = $AtomID;
    }
  }

  return (\@ComponentsAtomIDsRefs, \@ComponentsAtomIDsMapRefs);
}

# Get linear paths for connected components starting and ending at terminal aromatic atoms,
# which are connected to only one other aromatic atom in the connected component..
#
sub _GetConnectedComponentsPathsAtomsAndBondsForNonRingAromaticAtoms {
  my($This, $AtomIDsMapRef, $ComponentsAtomIDsRef, $ComponentsAtomIDsMapRef) = @_;
  my($Index, $AtomID, $Atom, $AtomNbr, $AtomNbrID, $NumOfNonRingAromaticNbrs, $AtomIndex1, $AtomIndex2, $AtomID1, $AtomID2, $Atom1, $Atom2, $AtomIDsGraph, $StartTerminalAtomID, $EndTerminalAtomID, @Paths, @PathAtomIDs, @PathsAtoms, @PathsBonds, @TerminalAtomIDs, @AtomIDs, @BondedAtomPairIDs);

  @PathsAtoms = ();
  @PathsBonds = ();

  @TerminalAtomIDs = ();

  $Index = 0;
  COMPONENT: for $Index (0 .. $#{$ComponentsAtomIDsRef}) {
    @{$TerminalAtomIDs[$Index]} = ();

    # Identify terminal atoms for connected components...
    #
    # Notes:
    #   . Terminal atoms are defined as atoms connected to only one marked
    #     aromatic atom.
    #   . Linear connected compoents contain only two terminal atoms.
    #
    ATOM: for $AtomID (@{$ComponentsAtomIDsRef->[$Index]}) {
      $Atom = $AtomIDsMapRef->{$AtomID};
      $NumOfNonRingAromaticNbrs = 0;

      ATOMNBRID: for $AtomNbr ($Atom->GetNeighbors()) {
	$AtomNbrID = $AtomNbr->GetID();

	# Is neighbor in the same connected components containing aromatic atoms?
	if (!exists $ComponentsAtomIDsMapRef->[$Index]{$AtomNbrID}) {
	  next ATOMNBRID;
	}
	$NumOfNonRingAromaticNbrs++;
      }

      # Is it a terminal atom?
      if ($NumOfNonRingAromaticNbrs != 1) {
	next ATOM;
      }
      push @{$TerminalAtomIDs[$Index]}, $AtomID;
    }

    if (@{$TerminalAtomIDs[$Index]} != 2) {
      next COMPONENT;
    }

    # Setup bonded atom pair IDs for connected component...
    #
    @AtomIDs = @{$ComponentsAtomIDsRef->[$Index]};
    @BondedAtomPairIDs = ();

    for $AtomIndex1 ( 0 .. $#AtomIDs) {
      $AtomID1 = $AtomIDs[$AtomIndex1];
      $Atom1 = $AtomIDsMapRef->{$AtomID1};

      for $AtomIndex2 ( ($AtomIndex1 + 1) .. $#AtomIDs) {
	$AtomID2 = $AtomIDs[$AtomIndex2];
	$Atom2 = $AtomIDsMapRef->{$AtomID2};

	if ($Atom1->IsBondedToAtom($Atom2)) {
	  push @BondedAtomPairIDs, ($AtomID1, $AtomID2);
	}
      }
    }

    if (!@BondedAtomPairIDs) {
      next COMPONENT;
    }

    # Get path for connected component...
    $AtomIDsGraph = new Graph(@AtomIDs);
    $AtomIDsGraph->AddEdges(@BondedAtomPairIDs);

    ($StartTerminalAtomID, $EndTerminalAtomID) = sort { $a <=> $b }  @{$TerminalAtomIDs[$Index]};
    @Paths = $AtomIDsGraph->GetPathsBetween($StartTerminalAtomID, $EndTerminalAtomID);

    if (@Paths != 1) {
      next COMPONENT;
    }

    @PathAtomIDs = $Paths[0]->GetVertices();

    my(@PathAtoms);
    @PathAtoms = $This->_GetAtomsFromAtomIDs(@PathAtomIDs);
    push @PathsAtoms, \@PathAtoms;

    my(@PathBonds);
    @PathBonds = $This->_GetPathBonds(@PathAtomIDs);
    push @PathsBonds, \@PathBonds;

  }

  return (\@PathsAtoms, \@PathsBonds);
}

# Setup initial processing status of atoms and bonds involved in connected paths
# before starting kekulization...
#
# Possible atom processing status: DoubleBondPossible, DoubleBondAssigned, DoubleBondNotPossible
# Initial status: DoubleBondPossible or DoubleBondNotPossible
#
# Possible bond processing status: DoubleBondAssigned, SingleBondAssigned, NotProcessed
#
# Possible paths processing status: Processed, NotProcessed
# Initial status: NotProcessed
#
sub _SetupConnectedPathSetsForKekulization {
  my($This, $PathAtomsSetsRef, $PathBondsSetsRef) = @_;
  my($PathIndex, $PathAtomsRef, $PathBondsRef, $Atom, $AtomID, $Bond, $BondID, %AtomProcessingStatusMap, %BondProcessingStatusMap, @PathsProcessingStatus, %InitialPathBondOrderMap);

  # Possible path set status values: Processed, NotProcessed
  # Initial value: NotProcessed
  #
  @PathsProcessingStatus = ('NotProcessed') x scalar @{$PathAtomsSetsRef};

  # Collect initial bond order of path bonds before setting bond orders to 1
  # and use it to set the bond order back to intial value after it has been processed for
  # availability of double bonds...
  #
  %InitialPathBondOrderMap = ();
  for $PathBondsRef (@{$PathBondsSetsRef}) {
    BOND: for $Bond (@{$PathBondsRef}) {
      $BondID = $Bond->GetID();
      if (exists $InitialPathBondOrderMap{$BondID}) {
	next BOND;
      }
      $InitialPathBondOrderMap{$BondID} = $Bond->GetBondOrder();
      $Bond->SetBondOrder(1);
    }
  }

  %AtomProcessingStatusMap = ();
  %BondProcessingStatusMap = ();

  for $PathIndex (0 .. $#{$PathAtomsSetsRef}) {

    $PathAtomsRef = $PathAtomsSetsRef->[$PathIndex];
    ATOM: for $Atom (@{$PathAtomsRef}) {
      $AtomID = $Atom->GetID();
      if (exists $AtomProcessingStatusMap{$AtomID}) {
	next ATOM;
      }
      $AtomProcessingStatusMap{$AtomID} = ($Atom->GetNumOfBondsAvailableForNonHydrogenAtoms() >= 1) ? 'DoubleBondPossible' : 'DoubleBondNotPossible';
    }

    $PathBondsRef = $PathBondsSetsRef->[$PathIndex];
    BOND: for $Bond (@{$PathBondsRef}) {
      $BondID = $Bond->GetID();
      if (exists $BondProcessingStatusMap{$BondID}) {
	next BOND;
      }
      $BondProcessingStatusMap{$BondID} = 'NotProcessed';
    }
  }

  # Set bond orders back to initial bond orders...
  for $PathIndex (0 .. $#{$PathAtomsSetsRef}) {
    $PathBondsRef = $PathBondsSetsRef->[$PathIndex];

    for $Bond (@{$PathBondsRef}) {
      $BondID = $Bond->GetID();
      if (exists $InitialPathBondOrderMap{$BondID}) {
	$Bond->SetBondOrder($InitialPathBondOrderMap{$BondID});
      }
    }
  }

  return (\%AtomProcessingStatusMap, \%BondProcessingStatusMap, \@PathsProcessingStatus);
}

# Kekulize connected path sets corresponding to fused rings, individual rings, or any other
# connected path...
#
# Note:
#   . PathAtomsRef and PathBondsRef contain paths and bonds corresponding to path
#     under consideration for kekulization
#   . PathAtomsSetsRef and PathBondsSetsRef contain any other available paths fused
#     to the path being kekulized
#   . _KekulizeConnectedPathSets is invoked recursively to kekulize all available paths
#
sub _KekulizeConnectedPathSets {
  my($This, $PathAtomsRef, $PathBondsRef, $AtomProcessingStatusMapRef, $BondProcessingStatusMapRef, $PathAtomsSetsRef, $PathBondsSetsRef, $PathsProcessingStatusRef) = @_;
  my($PathBond);

  # Get next available path bond...
  $PathBond = $This->_GetNextAvailablePathBondForKekulization($PathBondsRef, $BondProcessingStatusMapRef);

  if ($PathBond) {
    return $This->_ProcessNextAvailablePathBondForKekulization($PathBond, $PathAtomsRef, $PathBondsRef, $AtomProcessingStatusMapRef, $BondProcessingStatusMapRef, $PathAtomsSetsRef, $PathBondsSetsRef, $PathsProcessingStatusRef);
  }

  # Did kekulization succeed for the current path bonds?
  if (!$This->_DidKekulizationSucceedForPathBonds($PathAtomsRef, $PathBondsRef, $AtomProcessingStatusMapRef, $BondProcessingStatusMapRef)) {
    return 0;
  }

  # Is there any other path available for kekulization?
  ($PathAtomsRef, $PathBondsRef) = $This->_GetNextAvailablePathForKekulization($PathAtomsSetsRef, $PathBondsSetsRef, $PathsProcessingStatusRef, $AtomProcessingStatusMapRef, $BondProcessingStatusMapRef);

  if ($PathAtomsRef && $PathBondsRef) {
    # Recursively call itself to kekulize next path, which could either be a new path or part
    # of a fused paths corresponding to fused ring sets...
    #
    return $This->_KekulizeConnectedPathSets($PathAtomsRef, $PathBondsRef, $AtomProcessingStatusMapRef, $BondProcessingStatusMapRef, $PathAtomsSetsRef, $PathBondsSetsRef, $PathsProcessingStatusRef);
  }

  return 1;
}

# Get next available path bond in a list of path bonds...
#
sub _GetNextAvailablePathBondForKekulization {
  my($This, $PathBondsRef, $BondProcessingStatusMapRef) = @_;
  my($AvailablePathBond, $PathBond, $PathBondID);

  $AvailablePathBond = undef;

  BOND: for $PathBond (@{$PathBondsRef}) {
    $PathBondID = $PathBond->GetID();
    if (!exists $BondProcessingStatusMapRef->{$PathBondID}) {
      next BOND;
    }
    if ($BondProcessingStatusMapRef->{$PathBondID} =~ /^NotProcessed$/i) {
      $AvailablePathBond = $PathBond;
      last BOND;
    }
  }

  return ($AvailablePathBond);
}

# Process next available path bond for kekulizaiton...
#
sub _ProcessNextAvailablePathBondForKekulization {
  my($This, $PathBond, $PathAtomsRef, $PathBondsRef, $AtomProcessingStatusMapRef, $BondProcessingStatusMapRef, $PathAtomsSetsRef, $PathBondsSetsRef, $PathsProcessingStatusRef) = @_;
  my($PathBondID, $PathAtom1, $PathAtom2, $PathAtomID1, $PathAtomID2, %CurrentAtomProcessingStatusMap, %CurrentBondProcessingStatusMap);

  $PathBondID = $PathBond->GetID();

  ($PathAtom1, $PathAtom2) = $PathBond->GetAtoms();
  ($PathAtomID1, $PathAtomID2) = ($PathAtom1->GetID(), $PathAtom2->GetID());

  %CurrentAtomProcessingStatusMap = %{$AtomProcessingStatusMapRef};
  %CurrentBondProcessingStatusMap = %{$BondProcessingStatusMapRef};

  # Is it possible to assign a double bond to the current path bond?
  if ($AtomProcessingStatusMapRef->{$PathAtomID1} =~ /^DoubleBondPossible$/i && $AtomProcessingStatusMapRef->{$PathAtomID2} =~ /^DoubleBondPossible$/i ) {
    # Set current bond to double bond by appropriately marking atom and bond process status...
    $AtomProcessingStatusMapRef->{$PathAtomID1} = 'DoubleBondAssigned';
    $AtomProcessingStatusMapRef->{$PathAtomID2} = 'DoubleBondAssigned';

    $BondProcessingStatusMapRef->{$PathBondID} = 'DoubleBondAssigned';

    # Recursively call  _KekulizeConnectedPathSets to kekulize next available bond...
    if ($This->_KekulizeConnectedPathSets($PathAtomsRef, $PathBondsRef, $AtomProcessingStatusMapRef, $BondProcessingStatusMapRef, $PathAtomsSetsRef, $PathBondsSetsRef, $PathsProcessingStatusRef)) {
      return 1;
    }

    # Double bond at the current ring bond position didn't lead to successful kekulization...
    %{$AtomProcessingStatusMapRef} = %CurrentAtomProcessingStatusMap;
    %{$BondProcessingStatusMapRef} = %CurrentBondProcessingStatusMap;
  }

  # Try single bond at the current ring bond position and recursively call _KekulizeConnectedPathSets to kekulize
  # rest of the ring bonds...
  #
  $BondProcessingStatusMapRef->{$PathBondID} = 'SingleBondAssigned';

  if ($This->_KekulizeConnectedPathSets($PathAtomsRef, $PathBondsRef, $AtomProcessingStatusMapRef, $BondProcessingStatusMapRef, $PathAtomsSetsRef, $PathBondsSetsRef, $PathsProcessingStatusRef)) {
    return 1;
  }

  %{$AtomProcessingStatusMapRef} = %CurrentAtomProcessingStatusMap;
  %{$BondProcessingStatusMapRef} = %CurrentBondProcessingStatusMap;

  # Kekulization didn't work out for path bonds...

  return 0;

}

# Get next available path for kekulization from a set of fused ring paths...
#
sub _GetNextAvailablePathForKekulization {
  my($This, $PathAtomsSetsRef, $PathBondsSetsRef, $PathsProcessingStatusRef, $AtomProcessingStatusMapRef, $BondProcessingStatusMapRef) = @_;
  my($PathIndex, $AvailablePathIndex, $PathAtomsRef, $PathBondsRef, $PathBond, $PathBondID, $MaxNumOfPathBondsProcessed, $NumOfPathBondsProcessed);

  ($PathAtomsRef, $PathBondsRef, $AvailablePathIndex) = (undef) x 3;

  if (!(defined($PathAtomsSetsRef)  && defined($PathBondsSetsRef) && defined($PathsProcessingStatusRef))) {
    return ($PathAtomsRef, $PathBondsRef);
  }

  $MaxNumOfPathBondsProcessed = -999;
  $AvailablePathIndex = undef;

  PATHINDEX: for $PathIndex (0 .. $#{$PathsProcessingStatusRef}) {
    if ($PathsProcessingStatusRef->[$PathIndex] =~ /^Processed$/i) {
      next PATHINDEX;
    }

    # Count of already processed bonds in an unprocessed path bonds through
    # their participation in any fused bonds sets...
    #
    $NumOfPathBondsProcessed = 0;
    PATHBOND: for $PathBond (@{$PathBondsSetsRef->[$PathIndex]}) {
      $PathBondID = $PathBond->GetID();
      if ($BondProcessingStatusMapRef->{$PathBondID} =~ /^NotProcessed$/i) {
	next PATHBOND;
      }
      $NumOfPathBondsProcessed++;
    }

    if ($NumOfPathBondsProcessed > $MaxNumOfPathBondsProcessed) {
      $AvailablePathIndex = $PathIndex;
      $MaxNumOfPathBondsProcessed = $NumOfPathBondsProcessed;
    }

  }

  # Is any path available?
  if (!$AvailablePathIndex) {
    return ($PathAtomsRef, $PathBondsRef);
  }

  $PathsProcessingStatusRef->[$AvailablePathIndex] = 'Processed';

  $PathAtomsRef = $PathAtomsSetsRef->[$AvailablePathIndex];
  $PathBondsRef = $PathBondsSetsRef->[$AvailablePathIndex];

  return ($PathAtomsRef, $PathBondsRef);
}

# Check for kekulization in a specific set of path bonds. For successful kekulization, all
# all path atoms marked with DoubleBondPossible must be involved in a path double bond...
#
sub _DidKekulizationSucceedForPathBonds {
  my($This, $PathAtomsRef, $PathBondsRef, $AtomProcessingStatusMapRef, $BondProcessingStatusMapRef) = @_;
  my($PathAtom, $PathAtomID);

  for $PathAtom (@{$PathAtomsRef}) {
    $PathAtomID = $PathAtom->GetID();
    if (exists $AtomProcessingStatusMapRef->{$PathAtomID} && $AtomProcessingStatusMapRef->{$PathAtomID} =~ /^DoubleBondPossible$/i) {
      return 0;
    }
  }
  return 1;
}

# Assign bond orders to the bonds in a molecule which have been successfully
# kekulized along with optional clearing of aromaticty property...
#
sub _ProcessBondOrdersAssignedDuringSuccessfulKekulization {
  my($This, $BondsRef, $BondsProcessingStatusMapRef, $DeleteBondsAromaticity) = @_;
  my($Bond, $BondID, $BondOrder);

  $DeleteBondsAromaticity = defined $DeleteBondsAromaticity ? $DeleteBondsAromaticity : 0;

  BOND: for $Bond (@{$BondsRef}) {
    $BondID = $Bond->GetID();

    if (!exists $BondsProcessingStatusMapRef->{$BondID}) {
      carp "Warning: ${ClassName}->_ProcessBondOrdersAssignedDuringSuccessfulKekulization: Couldn't process bond with bond ID, $BondID: It's not available in the list of bonds processed for kekulization...";
      next BOND;
    }

    $BondOrder = ($BondsProcessingStatusMapRef->{$BondID} =~ /^DoubleBondAssigned$/i) ? 2 : 1;
    $Bond->SetBondOrder($BondOrder);

    if ($DeleteBondsAromaticity) {
      $Bond->DeleteAromatic();
    }
  }
  return $This;
}

# Does molecule contains aromatic rings?
#
sub HasAromaticRings {
  my($This) = @_;

  return $This->GetNumOfAromaticRings() ? 1 : 0;
}

# Does molecule contains any aromatic atom in a ring?
#
sub HasAromaticAtomsInRings {
  my($This) = @_;
  my($Atom);

  ATOM: for $Atom ($This->GetAtoms()) {
    if (!$Atom->IsAromatic()) {
      next ATOM;
    }
    if ($Atom->IsInRing()) {
      return 1;
    }
  }
  return 0;
}

# Does molecule contains any aromatic atom not in a ring?
#
sub HasAromaticAtomsNotInRings {
  my($This) = @_;
  my($Atom);

  ATOM: for $Atom ($This->GetAtoms()) {
    if (!$Atom->IsAromatic()) {
      next ATOM;
    }
    if ($Atom->IsNotInRing()) {
      return 1;
    }
  }
  return 0;
}

# Does molecule contains rings?
#
sub HasRings {
  my($This) = @_;

  return $This->IsCyclic();
}

# Does molecule contains only one ring?
#
sub HasOnlyOneRing {
  my($This) = @_;

  return $This->IsUnicyclic();
}

# Does molecule contains any rings?
#
sub HasNoRings {
  my($This) = @_;

  return $This->IsAcyclic();
}

# Get size of smallest ring...
#
sub GetSizeOfSmallestRing {
  my($This) = @_;

  return $This->GetSizeOfSmallestCycle();
}

# Get size of largest ring...
#
sub GetSizeOfLargestRing {
  my($This) = @_;

  return $This->GetSizeOfLargestCycle();
}

# Get number of rings...
#
sub GetNumOfRings {
  my($This) = @_;

  return $This->GetNumOfCycles();
}

# Get number of aromatic rings...
#
sub GetNumOfAromaticRings {
  my($This) = @_;
  my($NumOfRings);

  $NumOfRings = scalar $This->GetAromaticRings();

  return $NumOfRings;
}

# Get num of rings with odd size...
#
sub GetNumOfRingsWithOddSize {
  my($This) = @_;

  return $This->GetNumOfCyclesWithOddSize();
}

# Get num of rings with even size...
#
sub GetNumOfRingsWithEvenSize {
  my($This) = @_;

  return $This->GetNumOfCyclesWithEvenSize();
}

# Get num of rings with specified size...
#
sub GetNumOfRingsWithSize {
  my($This, $RingSize) = @_;

  return $This->GetNumOfCyclesWithSize($RingSize);
}

# Get num of rings with size less than a specified size...
#
sub GetNumOfRingsWithSizeLessThan {
  my($This, $RingSize) = @_;

  return $This->GetNumOfCyclesWithSizeLessThan($RingSize);
}

# Get num of rings with size greater than a specified size...
#
sub GetNumOfRingsWithSizeGreaterThan {
  my($This, $RingSize) = @_;

  return $This->GetNumOfCyclesWithSizeGreaterThan($RingSize);
}

# Get largest ring as an array containing ring atoms...
#
sub GetLargestRing {
  my($This) = @_;

  return $This->_GetRing($This->GetLargestCycle());
}

# Get smallest ring as an array containing ring atoms...
#
sub GetSmallestRing {
  my($This) = @_;

  return $This->_GetRing($This->GetSmallestCycle());
}

# Get rings as an array containing references to arrays with ring atoms...
#
sub GetRings {
  my($This) = @_;

  return $This->_GetRings($This->GetCycles());
}

# Get aromatic rings as an array containing references to arrays with ring atoms...
#
sub GetAromaticRings {
  my($This) = @_;

  return $This->_GetAromaticRings($This->GetCycles());
}

# Get odd size rings as an array containing references to arrays with ring atoms...
#
sub GetRingsWithOddSize {
  my($This) = @_;

  return $This->_GetRings($This->GetCyclesWithOddSize());
}

# Get even size rings as an array containing references to arrays with ring atoms...
#
sub GetRingsWithEvenSize {
  my($This) = @_;

  return $This->_GetRings($This->GetCyclesWithEvenSize());
}

# Get rings with a specific size as an array containing references to arrays with ring atoms...
#
sub GetRingsWithSize {
  my($This, $RingSize) = @_;

  return $This->_GetRings($This->GetCyclesWithSize($RingSize));
}

# Get rings with size less than a specific size as an array containing references to arrays with ring atoms...
#
sub GetRingsWithSizeLessThan {
  my($This, $RingSize) = @_;

  return $This->_GetRings($This->GetCyclesWithSizeLessThan($RingSize));
}

# Get rings with size greater than a specific size as an array containing references to arrays with ring atoms...
#
sub GetRingsWithSizeGreaterThan {
  my($This, $RingSize) = @_;

  return $This->_GetRings($This->GetCyclesWithSizeGreaterThan($RingSize));
}

# Generate an array of bond objects for an array of ring atoms and return an array
# of bond objects...
#
sub GetRingBonds {
  my($This, @RingAtoms) = @_;
  my(@Bonds);

  @Bonds = ();
  if (!@RingAtoms) {
    # Return an empty ring bonds list...
    return @Bonds;
  }

  my(@RingAtomIDs);

  @RingAtomIDs = ();
  @RingAtomIDs = $This->_GetAtomsIDsFromAtoms(@RingAtoms);
  if (!@RingAtomIDs) {
    carp "Warning: ${ClassName}->GetRingBonds: No ring bonds retrieved: Atom IDs couldn't be retrieved for specified atoms...";
    return @Bonds;
  }

  # Add start atom to the end to make it a cyclic path for ring: It's taken out during conversion
  # of cyclic path to a ring...
  push @RingAtomIDs, $RingAtomIDs[0];

  return $This->_GetPathBonds(@RingAtomIDs);
}

# Generate an array containing references to arrays of ring bond objects for rings specified
# in an array of references to ring atoms...
#
sub GetRingBondsFromRings {
  my($This, @RingAtomsSets) = @_;
  my($RingAtomsRef, @RingBondsSets);

  @RingBondsSets = ();
  for $RingAtomsRef  (@RingAtomsSets) {
    my(@RingBonds);
    @RingBonds = $This->GetRingBonds(@{$RingAtomsRef});

    push @RingBondsSets, \@RingBonds;
  }

  return @RingBondsSets;
}

# Does molecule has any fused rings?
#
sub HasFusedRings {
  my($This) = @_;

  return $This->HasFusedCycles();
}

# Get references to array of fused ring sets and non-fused rings. Fused ring sets array reference
# contains refernces to arrays of rings; Non-fused rings array reference contains references to
# arrays of ring atoms...
# rings.
#
sub GetFusedAndNonFusedRings {
  my($This) = @_;
  my($FusedCyclesSetsRef, $NonFusedCyclesRef, @FusedRingSets, @NonFusedRings);

  @FusedRingSets = (); @NonFusedRings = ();
  ($FusedCyclesSetsRef, $NonFusedCyclesRef) = $This->GetFusedAndNonFusedCycles();
  if (!(defined($FusedCyclesSetsRef) && defined($NonFusedCyclesRef))) {
    return (\@FusedRingSets, \@NonFusedRings);
  }
  my($FusedCyclesSetRef);

  for $FusedCyclesSetRef (@{$FusedCyclesSetsRef}) {
    my(@FusedRingSet);
    @FusedRingSet = ();
    @FusedRingSet = $This->_GetRings(@{$FusedCyclesSetRef});
    push @FusedRingSets, \@FusedRingSet;
  }

  @NonFusedRings = $This->_GetRings(@{$NonFusedCyclesRef});

  return (\@FusedRingSets, \@NonFusedRings);
}

# Get rings as an array containing references to arrays with ring atoms...
#
sub _GetRings {
  my($This, @CyclicPaths) = @_;
  my($CyclicPath, @Rings);

  @Rings = ();
  if (!@CyclicPaths) {
    return @Rings;
  }
  if (!@CyclicPaths) {
    # Return an empty ring list...
    return @Rings;
  }

  for $CyclicPath (@CyclicPaths) {
    my(@RingAtoms);
    @RingAtoms = ();
    push @RingAtoms, $This->_GetRing($CyclicPath);

    push @Rings, \@RingAtoms;
  }
  return @Rings;
}

# Get aromatic rings as an array containing references to arrays with ring atoms...
#
sub _GetAromaticRings {
  my($This, @CyclicPaths) = @_;
  my($RingAtomsRef, @Rings, @AromaticRings);

  @AromaticRings = ();
  @Rings = $This->_GetRings(@CyclicPaths);

  if (!@Rings) {
    return @AromaticRings;
  }
  RING: for $RingAtomsRef (@Rings) {
    if (!$This->IsRingAromatic(@{$RingAtomsRef})) {
      next RING;
    }
    my(@RingAtoms);
    @RingAtoms = ();
    push @RingAtoms, @{$RingAtomsRef};

    push @AromaticRings, \@RingAtoms;
  }
  return @AromaticRings;
}

# Map atom IDs in cyclic path to atoms and return a reference to an array containing ring atoms...
#
# Note:
#   . Start and end vertex is same for cyclic paths. So end atom is removed before
#     returning atoms array as ring atoms...
#
sub _GetRing {
  my($This, $CyclicPath) = @_;
  my(@RingAtoms);

  @RingAtoms = ();
  if (!defined $CyclicPath) {
    # Return an empty atoms list...
    return @RingAtoms;
  }

  @RingAtoms = $This->_GetPathAtoms($CyclicPath);
  if (@RingAtoms) {
    pop @RingAtoms;
  }
  return @RingAtoms;
}

# Map atom IDs to atoms and return a reference to an array containing these atoms...
#
sub _GetPathAtoms {
  my($This, $Path) = @_;
  my(@PathAtoms);

  @PathAtoms = ();
  if (!defined $Path) {
    carp "Warning: ${ClassName}->_GetPathAtoms: No path atoms retrieved: Path must be defined...";
    return @PathAtoms;
  }
  my(@AtomIDs);

  @AtomIDs = ();
  @AtomIDs = $Path->GetVertices();

  @PathAtoms = $This->_GetAtomsFromAtomIDs(@AtomIDs);

  return @PathAtoms;
}

# Get bonds for a path specified by atom IDs...
#
sub _GetPathBonds {
  my($This, @AtomIDs) = @_;
  my($Index, $AtomID1, $AtomID2, @Bonds, @EdgesAtomIDs);

  @Bonds = (); @EdgesAtomIDs = ();

  if (!@AtomIDs || @AtomIDs == 1) {
    return @Bonds;
  }

  # Setup edges...
  for $Index (0 .. ($#AtomIDs - 1) ) {
    $AtomID1 = $AtomIDs[$Index];
    $AtomID2 = $AtomIDs[$Index + 1];
    push @EdgesAtomIDs, ($AtomID1, $AtomID2);
  }
  @Bonds =  $This->GetEdgesProperty('Bond', @EdgesAtomIDs);

  return @Bonds;
}

# Map atom ID to an atom...
#
sub _GetAtomFromAtomID {
  my($This, $AtomID) = @_;

  return $This->GetVertexProperty('Atom', $AtomID);
}

# Map atom IDs to atoms and return an array containing these atoms...
#
sub _GetAtomsFromAtomIDs {
  my($This, @AtomIDs) = @_;

  return $This->GetVerticesProperty('Atom', @AtomIDs);
}

# Map atoms to atom IDs and return an array containing these atoms...
#
sub _GetAtomsIDsFromAtoms {
  my($This, @Atoms) = @_;

  return map { $_->GetID() } @Atoms;
}

# Get bonded atom pair atom IDs for specified list of atom IDs...
#
sub _GetBondedAtomPairAtomIDsFromAtomIDs {
  my($This, @AtomIDs) = @_;
  my($AtomIndex1, $AtomID1, $Atom1, $AtomIndex2, $AtomID2, $Atom2, @Atoms, @BondedAtomPairIDs);

  @BondedAtomPairIDs = ();
  @Atoms = $This->_GetAtomsFromAtomIDs(@AtomIDs);

  for $AtomIndex1 ( 0 .. $#Atoms) {
    $Atom1 = $Atoms[$AtomIndex1];
    $AtomID1 = $Atom1->GetID();

    ATOMINDEX2: for $AtomIndex2 ( ($AtomIndex1 + 1) .. $#Atoms) {
      $Atom2 = $Atoms[$AtomIndex2];
      if (!$Atom1->IsBondedToAtom($Atom2)) {
	next ATOMINDEX2;
       }
      $AtomID2 = $Atom2->GetID();

      push @BondedAtomPairIDs, ($AtomID1, $AtomID2);
    }
  }

  return @BondedAtomPairIDs;
}

# Get bonded atom pair atoms for specified list of atoms...
#
sub _GetBondedAtomPairAtomsFromAtoms {
  my($This, @Atoms) = @_;
  my($AtomIndex1, $Atom1, $AtomIndex2, $Atom2, @BondedAtomPairAtoms);

  @BondedAtomPairAtoms = ();

  for $AtomIndex1 ( 0 .. $#Atoms) {
    $Atom1 = $Atoms[$AtomIndex1];

    ATOMINDEX2: for $AtomIndex2 ( ($AtomIndex1 + 1) .. $#Atoms) {
      $Atom2 = $Atoms[$AtomIndex2];
      if ($Atom1->IsBondedToAtom($Atom2)) {
	next ATOMINDEX2;
       }

      push @BondedAtomPairAtoms, ($Atom1, $Atom2);
    }
  }

  return @BondedAtomPairAtoms;
}

# Is atom in a ring?
#
sub _IsAtomInRing {
  my($This, $Atom) = @_;

  return $This->IsCyclicVertex($Atom->GetID());
}

# Is atom not in a ring?
#
sub _IsAtomNotInRing {
  my($This, $Atom) = @_;

  return $This->IsAcyclicVertex($Atom->GetID());
}

# Is atom only in one ring?
#
sub _IsAtomInOnlyOneRing {
  my($This, $Atom) = @_;

  return $This->IsUnicyclicVertex($Atom->GetID());
}

# Is atom in a ring of specified size?
#
sub _IsAtomInRingOfSize {
  my($This, $Atom, $RingSize) = @_;

  return $This->GetNumOfVertexCyclesWithSize($Atom->GetID(), $RingSize) ? 1 : 0;
}

# Get size of smallest ring containing specified atom...
#
sub _GetSizeOfSmallestAtomRing {
  my($This, $Atom) = @_;

  return $This->GetSizeOfSmallestVertexCycle($Atom->GetID());
}

# Get size of largest ring containing specified atom...
#
sub _GetSizeOfLargestAtomRing {
  my($This, $Atom) = @_;

  return $This->GetSizeOfLargestVertexCycle($Atom->GetID());
}

# Get number of  rings containing specified atom...
#
sub _GetNumOfAtomRings {
  my($This, $Atom) = @_;

  return $This->GetNumOfVertexCycles($Atom->GetID());
}

# Get number of  rings with odd size containing specified atom...
#
sub _GetNumOfAtomRingsWithOddSize {
  my($This, $Atom) = @_;

  return $This->GetNumOfVertexCyclesWithOddSize($Atom->GetID());
}

# Get number of  rings with even size containing specified atom...
#
sub _GetNumOfAtomRingsWithEvenSize {
  my($This, $Atom) = @_;

  return $This->GetNumOfVertexCyclesWithEvenSize($Atom->GetID());
}

# Get number of  rings with specified size containing specified atom...
#
sub _GetNumOfAtomRingsWithSize {
  my($This, $Atom, $RingSize) = @_;

  return $This->GetNumOfVertexCyclesWithSize($Atom->GetID(), $RingSize);
}

# Get number of  rings with size less than specified containing specified atom...
#
sub _GetNumOfAtomRingsWithSizeLessThan {
  my($This, $Atom, $RingSize) = @_;

  return $This->GetNumOfVertexCyclesWithSizeLessThan($Atom->GetID(), $RingSize);
}

# Get number of  rings with size greater than specified containing specified atom...
#
sub _GetNumOfAtomRingsWithSizeGreaterThan {
  my($This, $Atom, $RingSize) = @_;

  return $This->GetNumOfVertexCyclesWithSizeGreaterThan($Atom->GetID(), $RingSize);
}

# Get smallest ring as an array containing ring atoms...
#
sub _GetSmallestAtomRing {
  my($This, $Atom) = @_;

  return $This->_GetRing($This->GetSmallestVertexCycle($Atom->GetID()));
}

# Get odd size rings an array of references to arrays containing ring atoms...
#
sub _GetLargestAtomRing {
  my($This, $Atom) = @_;

  return $This->_GetRing($This->GetLargestVertexCycle($Atom->GetID()));
}

# Get all rings an array of references to arrays containing ring atoms...
#
sub _GetAtomRings {
  my($This, $Atom) = @_;

  return $This->_GetRings($This->GetVertexCycles($Atom->GetID()));
}

# Get odd size rings an array of references to arrays containing ring atoms...
#
sub _GetAtomRingsWithOddSize {
  my($This, $Atom) = @_;

  return $This->_GetRings($This->GetVertexCyclesWithOddSize($Atom->GetID()));
}

# Get even size rings an array of references to arrays containing ring atoms...
#
sub _GetAtomRingsWithEvenSize {
  my($This, $Atom) = @_;

  return $This->_GetRings($This->GetVertexCyclesWithEvenSize($Atom->GetID()));
}

# Get rings with specified size  an array of references to arrays containing ring atoms...
#
sub _GetAtomRingsWithSize {
  my($This, $Atom, $RingSize) = @_;

  return $This->_GetRings($This->GetVertexCyclesWithSize($Atom->GetID(), $RingSize));
}

# Get rings with size less than specfied size as an array of references to arrays containing ring atoms...
#
sub _GetAtomRingsWithSizeLessThan {
  my($This, $Atom, $RingSize) = @_;

  return $This->_GetRings($This->GetVertexCyclesWithSizeLessThan($Atom->GetID(), $RingSize));
}

# Get rings with size less than specfied size as an array of references to arrays containing ring atoms...
#
sub _GetAtomRingsWithSizeGreaterThan {
  my($This, $Atom, $RingSize) = @_;

  return $This->_GetRings($This->GetVertexCyclesWithSizeGreaterThan($Atom->GetID(), $RingSize));
}

# Is bond in a ring?
#
sub _IsBondInRing {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->IsCyclicEdge($Atom1->GetID(), $Atom2->GetID());
}

# Is bond not in a ring?
#
sub _IsBondNotInRing {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->IsAcyclicEdge($Atom1->GetID(), $Atom2->GetID());
}

# Is bond only in one ring?
#
sub _IsBondInOnlyOneRing {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->IsUnicyclicEdge($Atom1->GetID(), $Atom2->GetID());
}

# Is bond in a ring of specified size?
#
sub _IsBondInRingOfSize {
  my($This, $Bond, $RingSize) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->GetNumOfEdgeCyclesWithSize($Atom1->GetID(), $Atom2->GetID(), $RingSize) ? 1 : 0;
}

# Get size of smallest ring containing specified bond...
#
sub _GetSizeOfSmallestBondRing {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->GetSizeOfSmallestEdgeCycle($Atom1->GetID(), $Atom2->GetID());
}

# Get size of largest ring containing specified bond...
#
sub _GetSizeOfLargestBondRing {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->GetSizeOfLargestEdgeCycle($Atom1->GetID(), $Atom2->GetID());
}

# Get number of  rings containing specified bond...
#
sub _GetNumOfBondRings {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->GetNumOfEdgeCycles($Atom1->GetID(), $Atom2->GetID());
}

# Get number of  rings with odd size containing specified bond...
#
sub _GetNumOfBondRingsWithOddSize {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->GetNumOfEdgeCyclesWithOddSize($Atom1->GetID(), $Atom2->GetID());
}

# Get number of  rings with even size containing specified bond...
#
sub _GetNumOfBondRingsWithEvenSize {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->GetNumOfEdgeCyclesWithEvenSize($Atom1->GetID(), $Atom2->GetID());
}

# Get number of  rings with specified size containing specified bond...
#
sub _GetNumOfBondRingsWithSize {
  my($This, $Bond, $RingSize) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->GetNumOfEdgeCyclesWithSize($Atom1->GetID(), $Atom2->GetID(), $RingSize);
}

# Get number of  rings with size less than specified containing specified bond...
#
sub _GetNumOfBondRingsWithSizeLessThan {
  my($This, $Bond, $RingSize) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->GetNumOfEdgeCyclesWithSizeLessThan($Atom1->GetID(), $Atom2->GetID(), $RingSize);
}

# Get number of  rings with size greater than specified containing specified bond...
#
sub _GetNumOfBondRingsWithSizeGreaterThan {
  my($This, $Bond, $RingSize) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->GetNumOfEdgeCyclesWithSizeGreaterThan($Atom1->GetID(), $Atom2->GetID(), $RingSize);
}

# Get smallest ring as an array containing ring atoms...
#
sub _GetSmallestBondRing {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->_GetRing($This->GetSmallestEdgeCycle($Atom1->GetID(), $Atom2->GetID()));
}

# Get odd size rings an array of references to arrays containing ring atoms...
#
sub _GetLargestBondRing {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->_GetRing($This->GetLargestEdgeCycle($Atom1->GetID(), $Atom2->GetID()));
}

# Get all rings an array of references to arrays containing ring atoms...
#
sub _GetBondRings {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->_GetRings($This->GetEdgeCycles($Atom1->GetID(), $Atom2->GetID()));
}

# Get odd size rings an array of references to arrays containing ring atoms...
#
sub _GetBondRingsWithOddSize {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->_GetRings($This->GetEdgeCyclesWithOddSize($Atom1->GetID(), $Atom2->GetID()));
}

# Get even size rings an array of references to arrays containing ring atoms...
#
sub _GetBondRingsWithEvenSize {
  my($This, $Bond) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->_GetRings($This->GetEdgeCyclesWithEvenSize($Atom1->GetID(), $Atom2->GetID()));
}

# Get rings with specified size  an array of references to arrays containing ring atoms...
#
sub _GetBondRingsWithSize {
  my($This, $Bond, $RingSize) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->_GetRings($This->GetEdgeCyclesWithSize($Atom1->GetID(), $Atom2->GetID(), $RingSize));
}

# Get rings with size less than specfied size as an array of references to arrays containing ring atoms...
#
sub _GetBondRingsWithSizeLessThan {
  my($This, $Bond, $RingSize) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->_GetRings($This->GetEdgeCyclesWithSizeLessThan($Atom1->GetID(), $Atom2->GetID(), $RingSize));
}

# Get rings with size less than specfied size as an array of references to arrays containing ring atoms...
#
sub _GetBondRingsWithSizeGreaterThan {
  my($This, $Bond, $RingSize) = @_;
  my($Atom1, $Atom2);

  ($Atom1, $Atom2) = $Bond->GetAtoms();

  return $This->_GetRings($This->GetEdgeCyclesWithSizeGreaterThan($Atom1->GetID(), $Atom2->GetID(), $RingSize));
}


# Get atom paths starting from a specified atom as a reference to an array containing references
# to arrays with path atoms.
#
# Path atoms atoms correspond to to all possible paths for specified atom in molecule with length
# upto a specified length and sharing of bonds in paths traversed. By default, rings are
# included in paths. A path containing a ring is terminated at an atom completing the ring.
#
# Note:
#    . For molecule without any rings, this method returns the same set of atom paths
#      as GetAtomPathsStartingAtWithLengthUpto method.
#
sub GetAllAtomPathsStartingAtWithLengthUpto {
  my($This, $StartAtom, $Length, $AllowCycles) = @_;

  return $This->_GetAtomPathsStartingAt('AllAtomPathsWithLengthUpto', $StartAtom, $Length, $AllowCycles);
}

# Get atom paths starting from a specified atom as a reference to an array containing references
# to arrays with path atoms.
#
# Path atoms atoms correspond to to all possible paths for specified atom in molecule with
# specified length and sharing of bonds in paths traversed. By default, rings are
# included in paths. A path containing a ring is terminated at an atom completing the ring.
#
# Note:
#    . For molecule without any rings, this method returns the same set of atom paths
#      as GetAtomPathsStartingAtWithLengthUpto method.
#
sub GetAllAtomPathsStartingAtWithLength {
  my($This, $StartAtom, $Length, $AllowCycles) = @_;

  return $This->_GetAtomPathsStartingAt('AllAtomPathsWithLength', $StartAtom, $Length, $AllowCycles);
}

# Get atom paths starting from a specified atom as a reference to an array containing references
# to arrays with path atoms.
#
# Path atoms atoms correspond to to all possible paths for specified atom in molecule with all
# possible lengths and sharing of bonds in paths traversed. By default, rings are
# included in paths. A path containing a ring is terminated at an atom completing the ring.
#
# Note:
#    . For molecule without any rings, this method returns the same set of atom paths
#      as GetAtomPathsStartingAt method.
#
sub GetAllAtomPathsStartingAt {
  my($This, $StartAtom, $AllowCycles) = @_;

  return $This->_GetAtomPathsStartingAt('AllAtomPathsWithAllLengths', $StartAtom, undef, $AllowCycles);
}

# Get atom paths starting from a specified atom as a reference to an array containing references
# to arrays with path atoms.
#
# Path atoms atoms correspond to to all possible paths for specified atom in molecule with length
# upto a specified length and no sharing of bonds in paths traversed. By default, rings are
# included in paths. A path containing a ring is terminated at an atom completing the ring.
#
sub GetAtomPathsStartingAtWithLengthUpto {
  my($This, $StartAtom, $Length, $AllowCycles) = @_;

  return $This->_GetAtomPathsStartingAt('AtomPathsWithLengthUpto', $StartAtom, $Length, $AllowCycles);
}

# Get atom paths starting from a specified atom as a reference to an array containing references
# to arrays with path atoms.
#
# Path atoms atoms correspond to to all possible paths for specified atom in molecule with
# specified length and no sharing of bonds in paths traversed. By default, rings are
# included in paths. A path containing a ring is terminated at an atom completing the ring.
#
sub GetAtomPathsStartingAtWithLength {
  my($This, $StartAtom, $Length, $AllowCycles) = @_;

  return $This->_GetAtomPathsStartingAt('AtomPathsWithLength', $StartAtom, $Length, $AllowCycles);
}

# Get atom paths starting from a specified atom as a reference to an array containing references
# to arrays with path atoms.
#
# Path atoms atoms correspond to to all possible paths for specified atom in molecule with all
# possible lengths and no sharing of bonds in paths traversed. By default, rings are
# included in paths. A path containing a ring is terminated at an atom completing the ring.
#
#
sub GetAtomPathsStartingAt {
  my($This, $StartAtom, $AllowCycles) = @_;

  return $This->_GetAtomPathsStartingAt('AtomPathsWithAllLengths', $StartAtom, undef, $AllowCycles);
}

# Get atom paths as an array containing references to arrays with path atoms...
#
sub _GetAtomPathsStartingAt {
  my($This, $Mode, $StartAtom, $Length, $AllowCycles) = @_;
  my(@AtomPaths);

  @AtomPaths = ();
  if (!defined $StartAtom) {
    carp "Warning: ${ClassName}->_GetAtomPathsStartingAt: No atom paths retrieved: Start atom is not defined...";
    return @AtomPaths;
  }
  if (!$This->HasAtom($StartAtom)) {
    carp "Warning: ${ClassName}->_GetAtomPathsStartingAt: No atom paths retrieved: Start atom doesn't exist...";
    return @AtomPaths;
  }
  my($StartAtomID, @Paths);

  $StartAtomID = $StartAtom->GetID();
  @Paths = ();

  # Collect appropriate atom paths...
  MODE: {
    if ($Mode =~ /^AtomPathsWithLengthUpto$/i) { @Paths = $This->GetPathsStartingAtWithLengthUpto($StartAtomID, $Length, $AllowCycles); last MODE; }
    if ($Mode =~ /^AtomPathsWithLength$/i) { @Paths = $This->GetPathsStartingAtWithLength($StartAtomID, $Length, $AllowCycles); last MODE; }
    if ($Mode =~ /^AtomPathsWithAllLengths$/i) { @Paths = $This->GetPathsStartingAt($StartAtomID, $AllowCycles); last MODE; }

    if ($Mode =~ /^AllAtomPathsWithLengthUpto$/i) { @Paths = $This->GetAllPathsStartingAtWithLengthUpto($StartAtomID, $Length, $AllowCycles); last MODE; }
    if ($Mode =~ /^AllAtomPathsWithLength$/i) { @Paths = $This->GetAllPathsStartingAtWithLength($StartAtomID, $Length, $AllowCycles); last MODE; }
    if ($Mode =~ /^AllAtomPathsWithAllLengths$/i) { @Paths = $This->GetAllPathsStartingAt($StartAtomID, $AllowCycles); last MODE; }

    print "Warn: ${ClassName}->_GetAtomPathsStartingAt: No atom paths retrieved: Mode, $Mode, is not supported...";
    return @AtomPaths;
  }
  return $This->_GetAtomPathsFromPaths(\@Paths);
}

# Get atom paths for all atoms as a reference to an array containing references to arrays with
# path atoms.
#
# Path atoms correspond to to all possible paths for each atom in molecule with length
# upto a specified length and sharing of bonds in paths traversed. By default, rings are
# included in paths. A path containing a ring is terminated at an atom completing the ring.
#
# Notes:
#    . For molecule without any rings, this method returns the same set of atom paths
#      as GetAtomPathsWithLengthUpto method.
#
sub GetAllAtomPathsWithLengthUpto {
  my($This, $Length, $AllowCycles) = @_;

  return $This->_GetAtomPaths('AllAtomPathsWithLengthUpto', $Length, $AllowCycles);
}

# Get atom paths for all atoms as a reference to an array containing references to arrays with
# path atoms.
#
# Path atoms correspond to to all possible paths for each atom in molecule with
# a specified length and sharing of bonds in paths traversed. By default, rings are
# included in paths. A path containing a ring is terminated at an atom completing the ring.
#
# Notes:
#    . For molecule without any rings, this method returns the same set of atom paths
#      as GetAtomPathsWithLengthUpto method.
#
sub GetAllAtomPathsWithLength {
  my($This, $Length, $AllowCycles) = @_;

  return $This->_GetAtomPaths('AllAtomPathsWithLength', $Length, $AllowCycles);
}

# Get atom paths for all atoms as a reference to an array containing references to arrays with
# path atoms.
#
# Path atoms correspond to to all possible paths for each atom in molecule with all
# possible lengths and sharing of bonds in paths traversed. By default, rings are
# included in paths. A path containing a ring is terminated at an atom completing the ring.
#
# Notes:
#    . For molecule without any rings, this method returns the same set of atom paths
#      as GetAtomPaths method.
#
sub GetAllAtomPaths {
  my($This, $AllowCycles) = @_;

  return $This->_GetAtomPaths('AllAtomPathsWithAllLengths', undef, $AllowCycles);
}

# Get atom paths for all atoms as a reference to an array containing references to arrays with
# path atoms.
#
# Path atoms correspond to to all possible paths for each atom in molecule with length
# upto a specified length and no sharing of bonds in paths traversed. By default, rings are
# included in paths. A path containing a ring is terminated at an atom completing the ring.
#
sub GetAtomPathsWithLengthUpto {
  my($This, $Length, $AllowCycles) = @_;

  return $This->_GetAtomPaths('AtomPathsWithLengthUpto', $Length, $AllowCycles);
}

# Get atom paths for all atoms as a reference to an array containing references to arrays with
# path atoms.
#
# Path atoms correspond to to all possible paths for each atom in molecule with
# a specified length and no sharing of bonds in paths traversed. By default, rings are
# included in paths. A path containing a ring is terminated at an atom completing the ring.
#
sub GetAtomPathsWithLength {
  my($This, $Length, $AllowCycles) = @_;

  return $This->_GetAtomPaths('AtomPathsWithLength', $Length, $AllowCycles);
}


# Get atom paths for all atoms as a reference to an array containing references to arrays with
# path atoms.
#
# Path atoms correspond to to all possible paths for each atom in molecule with all
# possible lengths and no sharing of bonds in paths traversed. By default, rings are
# included in paths. A path containing a ring is terminated at an atom completing the ring.
#
sub GetAtomPaths {
  my($This, $AllowCycles) = @_;

  return $This->_GetAtomPaths('AtomPathsWithAllLengths', undef, $AllowCycles);
}

# Get atom paths for all atoms as a reference to an array containing references to arrays with
# path atoms.
#
sub _GetAtomPaths {
  my($This, $Mode, $Length, $AllowCycles) = @_;
  my($PathsRef, @AtomPaths);

  @AtomPaths = ();
  # Collect appropriate atom paths...
  MODE: {
    if ($Mode =~ /^AtomPathsWithLengthUpto$/i) { $PathsRef = $This->GetPathsWithLengthUpto($Length, $AllowCycles); last MODE; }
    if ($Mode =~ /^AtomPathsWithLength$/i) { $PathsRef = $This->GetPathsWithLength($Length, $AllowCycles); last MODE; }
    if ($Mode =~ /^AtomPathsWithAllLengths$/i) { $PathsRef = $This->GetPaths($AllowCycles); last MODE; }

    if ($Mode =~ /^AllAtomPathsWithLengthUpto$/i) { $PathsRef = $This->GetAllPathsWithLengthUpto($Length, $AllowCycles); last MODE; }
    if ($Mode =~ /^AllAtomPathsWithLength$/i) { $PathsRef = $This->GetAllPathsWithLength($Length, $AllowCycles); last MODE; }
    if ($Mode =~ /^AllAtomPathsWithAllLengths$/i) { $PathsRef = $This->GetAllPaths($AllowCycles); last MODE; }

    print "Warn: ${ClassName}->_GetAtomPaths: No atom paths retrieved: Mode, $Mode, is not supported...";
    return \@AtomPaths;
  }
  return $This->_GetAtomPathsFromPaths($PathsRef);
}

# Get atom paths as an array reference containing references to arrays with path atoms...
#
sub _GetAtomPathsFromPaths {
  my($This, $PathsRef) = @_;
  my($Path, @AtomPaths);

  @AtomPaths = ();
  if (!defined $PathsRef) {
    return \@AtomPaths;
  }
  if (!@{$PathsRef}) {
    # Return an empty atom paths list...
    return \@AtomPaths;
  }
  for $Path (@{$PathsRef}) {
    my(@PathAtoms);
    @PathAtoms = ();
    @PathAtoms = $This->_GetAtomPathFromPath($Path);

    push @AtomPaths, \@PathAtoms;
  }
  return \@AtomPaths;
}

# Generate an array of bond objects for an array of path atoms and return an array
# of bond objects...
#
sub GetAtomPathBonds {
  my($This, @PathAtoms) = @_;
  my(@Bonds);

  if (!@PathAtoms) {
    # Return an empty ring bonds list...
    return @Bonds;
  }
  my(@PathAtomIDs);

  @PathAtomIDs = ();
  @PathAtomIDs = $This->_GetAtomsIDsFromAtoms(@PathAtoms);

  return $This->_GetPathBonds(@PathAtomIDs);
}

# Map atom IDs in path to atoms and return a reference to an array containing ring atoms...
#
sub _GetAtomPathFromPath {
  my($This, $Path) = @_;
  my(@PathAtoms);

  @PathAtoms = ();
  if (!defined $Path) {
    # Return an empty atoms list...
    return @PathAtoms;
  }

  return $This->_GetPathAtoms($Path);
}

# Get atom paths between two specified atoms as a reference to an array containing references
# to arrays with path atoms. For molecules with rings, atom paths array contains may contain
# two paths.
#
sub GetAtomPathsBetween {
  my($This, $StartAtom, $EndAtom) = @_;
  my(@AtomPaths);

  @AtomPaths = ();
  if (!(defined($StartAtom) && $This->HasAtom($StartAtom))) {
    carp "Warning: ${ClassName}->_GetAtomPathsBetween: No atom paths retrieved: Start atom is not defined  or it doesn't exist...";
    return @AtomPaths;
  }
  if (!(defined($EndAtom) && $This->HasAtom($EndAtom))) {
    carp "Warning: ${ClassName}->_GetAtomPathsBetween: No atom paths retrieved: End atom is not defined  or it doesn't exist...";
    return @AtomPaths;
  }
  return $This->_GetAtomPathsBetween($StartAtom, $EndAtom);
}

# Get atom paths between two specified atoms as a reference to an array containing references
# to arrays with path atoms.
#
sub _GetAtomPathsBetween {
  my($This, $StartAtom, $EndAtom) = @_;
  my($StartAtomID, $EndAtomID, @Paths);

  $StartAtomID = $StartAtom->GetID();
  $EndAtomID = $EndAtom->GetID();

  @Paths = ();
  @Paths = $This->GetPathsBetween($StartAtomID, $EndAtomID);

  return $This->_GetAtomPathsFromPaths(\@Paths);
}

# Get atom neighborhoods around a specified atom as an array containing references
# to arrays with neighborhood atoms at different radii upto specified radius...
#
sub GetAtomNeighborhoodsWithRadiusUpto {
  my($This, $StartAtom, $Radius) = @_;

  return $This->_GetAtomNeighborhoods('RadiusUpto', $StartAtom, $Radius);
}

# Get atom neighborhoods around a specified atom as an array containing references
# to arrays with neighborhood atoms at possible radii...
#
sub GetAtomNeighborhoods {
  my($This, $StartAtom) = @_;

  return $This->_GetAtomNeighborhoods('AllRadii', $StartAtom, undef);
}

# Get atom neighborhood around a specified atom, along with their successor connected atoms, collected
# with in a specified radius as a list containing references to lists with first value corresponding to neighborhood
# atom at a specific radius and second value as reference to a list containing its successor connected atoms.
#
# For a neighborhood atom at each radius level, the successor connected atoms correspond to the
# neighborhood atoms at the next radius level. Consequently, the neighborhood atoms at the last
# radius level don't contain any successor atoms which fall outside the range of specified radius.
#
sub GetAtomNeighborhoodsWithSuccessorAtomsAndRadiusUpto {
  my($This, $StartAtom, $Radius) = @_;

  return $This->_GetAtomNeighborhoods('WithSuccessorsAndRadiusUpto', $StartAtom, $Radius);
}

# Get atom neighborhood around a specified atom, along with their successor connected atoms, collected
# at all radii as a list containing references to lists with first value corresponding to neighborhood
# atom at a specific radius and second value as reference to a list containing its successor connected atoms.
#
# For a neighborhood atom at each radius level, the successor connected atoms correspond to the
# neighborhood atoms at the next radius level. Consequently, the neighborhood atoms at the last
# radius level don't contain any successor atoms which fall outside the range of specified radius.
#
#
sub GetAtomNeighborhoodsWithSuccessorAtoms {
  my($This, $StartAtom) = @_;

  return $This->_GetAtomNeighborhoods('WithSuccessorsAndAllRadii', $StartAtom, undef);
}

# Get atom neighborhoods...
#
sub _GetAtomNeighborhoods {
  my($This, $Mode, $StartAtom, $Radius) = @_;
  my(@AtomNeighborhoods);

  @AtomNeighborhoods = ();

  if (!(defined($StartAtom) && $This->HasAtom($StartAtom))) {
    carp "Warning: ${ClassName}->_GetAtomNeighborhoods: No atom neighborhoods retrieved: Start atom is not defined  or it doesn't exist...";
    return @AtomNeighborhoods;
  }
  if ($Mode =~ /^(RadiusUpto|WithSuccessorsAndRadiusUpto)$/i) {
    if (!(defined($Radius) && $Radius > 0)) {
      carp "Warning: ${ClassName}->_GetAtomNeighborhoods: No atom neighborhoods retrieved: Radius is not defined or it's <= 0 ...";
      return @AtomNeighborhoods;
    }
  }

  # Collect neighborhood atom IDs...
  my($StartAtomID, @NeighborhoodAtomIDs, @NeighborhoodAtomIDsWithSuccessors);

  @NeighborhoodAtomIDs = (); @NeighborhoodAtomIDsWithSuccessors = ();
  $StartAtomID = $StartAtom->GetID();

  MODE: {
    if ($Mode =~ /^RadiusUpto$/i) { @NeighborhoodAtomIDs = $This->GetNeighborhoodVerticesWithRadiusUpto($StartAtomID, $Radius); last MODE; }
    if ($Mode =~ /^AllRadii$/i) { @NeighborhoodAtomIDs = $This->GetNeighborhoodVertices($StartAtomID); last MODE; }

    if ($Mode =~ /^WithSuccessorsAndRadiusUpto$/i) { @NeighborhoodAtomIDsWithSuccessors = $This->GetNeighborhoodVerticesWithSuccessorsAndRadiusUpto($StartAtomID, $Radius); last MODE; }
    if ($Mode =~ /^WithSuccessorsAndAllRadii$/i) { @NeighborhoodAtomIDsWithSuccessors = $This->GetNeighborhoodVerticesWithSuccessors($StartAtomID); last MODE; }

    print "Warn: ${ClassName}->_GetAtomNeighborhood: No atom neighborhoods retrieved: Mode, $Mode, is not supported...";
    return @AtomNeighborhoods;
  }
  if ($Mode =~ /^(RadiusUpto|AllRadii)$/i) {
    return $This->_GetNeighborhoodAtomsFromAtomIDs(\@NeighborhoodAtomIDs);
  }
  elsif ($Mode =~ /^(WithSuccessorsAndRadiusUpto|WithSuccessorsAndAllRadii)$/i) {
    return $This->_GetNeighborhoodAtomsWithSuccessorsFromAtomIDs(\@NeighborhoodAtomIDsWithSuccessors);
  }

  return @AtomNeighborhoods;
}

# Map neighborhood atom IDs to atoms...
#
sub _GetNeighborhoodAtomsFromAtomIDs {
  my($This, $NeighborhoodsAtomIDsRef) = @_;
  my($NeighborhoodAtomIDsRef, @AtomNeighborhoods);

  @AtomNeighborhoods = ();
  for $NeighborhoodAtomIDsRef (@{$NeighborhoodsAtomIDsRef}) {
    my(@AtomNeighborhood);

    @AtomNeighborhood = ();
    @AtomNeighborhood = $This->_GetAtomsFromAtomIDs(@{$NeighborhoodAtomIDsRef});
    push @AtomNeighborhoods, \@AtomNeighborhood;
  }
  return @AtomNeighborhoods;
}

# Map neighborhood atom IDs with successors to atoms...
#
sub _GetNeighborhoodAtomsWithSuccessorsFromAtomIDs {
  my($This, $NeighborhoodsAtomIDsWithSuccessorsRef) = @_;
  my($Depth, $NeighborhoodAtomIDsWithSuccessorsRef, $NeighborhoodAtomIDWithSuccessorsRef, $NeighborhoodAtomID, $NeighborhoodAtomSuccessorsIDsRef, @AtomNeighborhoods);

  $Depth = 0;
  @AtomNeighborhoods = ();

  # Go over neighborhoods at each level...
  for $NeighborhoodAtomIDsWithSuccessorsRef (@{$NeighborhoodsAtomIDsWithSuccessorsRef}) {
    @{$AtomNeighborhoods[$Depth]} = ();

    # Go over the neighborhood atoms and their successors at a specific level..
    for $NeighborhoodAtomIDWithSuccessorsRef (@{$NeighborhoodAtomIDsWithSuccessorsRef}) {
      my($NeighborhoodAtom, @NeighborhoodAtomWithSuccessors, @NeighborhoodAtomSuccessorAtoms);

      @NeighborhoodAtomWithSuccessors = (); @NeighborhoodAtomSuccessorAtoms = ();
      ($NeighborhoodAtomID, $NeighborhoodAtomSuccessorsIDsRef) = @{$NeighborhoodAtomIDWithSuccessorsRef};

      # Map atom IDs to atoms...
      $NeighborhoodAtom = $This->_GetAtomFromAtomID($NeighborhoodAtomID);
      if (@{$NeighborhoodAtomSuccessorsIDsRef}) {
	@NeighborhoodAtomSuccessorAtoms = $This->_GetAtomsFromAtomIDs(@{$NeighborhoodAtomSuccessorsIDsRef});
      }

      # Store an atom and its successors at each level in an array...
      push @NeighborhoodAtomWithSuccessors, ($NeighborhoodAtom, \@NeighborhoodAtomSuccessorAtoms);

      push @{$AtomNeighborhoods[$Depth]} , \@NeighborhoodAtomWithSuccessors;
    }
      $Depth++;
  }
  return @AtomNeighborhoods;
}

# Get next object ID...
sub _GetNewObjectID {
  $ObjectID++;
  return $ObjectID;
}

# Is aromatic property set for the molecule?
sub IsAromatic {
  my($This) = @_;
  my($Aromatic);

  $Aromatic = $This->GetAromatic();

  return (defined($Aromatic) && $Aromatic) ? 1 : 0;
}

# Does molecule contains any atoms with non-zero Z coordiantes?
sub IsThreeDimensional {
  my($This) = @_;
  my($Atom, @Atoms);

  @Atoms = $This->GetAtoms();
  ATOM: for $Atom (@Atoms) {
      if ($Atom->GetZ() != 0) {
	return 1;
      }
  }
  return 0;
}

# Does molecule contains any atoms with non-zero X or Y coordinates
# and only zero Z-coordinates?
sub IsTwoDimensional {
  my($This) = @_;
  my($Atom, @Atoms);

  @Atoms = $This->GetAtoms();
  ATOM: for $Atom (@Atoms) {
      if ($Atom->GetZ() != 0) {
	return 0;
      }
      if ($Atom->GetX() != 0 || $Atom->GetY() != 0) {
	return 1;
      }
  }
  return 0;
}

# Get dimensionality of the molecule using one of the following two methods:
#   . Using explicitly set Dimensionality
#   . Going over atomic coordinates
#
# The valid dimensionality values are:
#   . 3D - Three dimensional: One of X, Y or Z coordinate is non-zero
#   . 2D - Two dimensional: One of X or Y coordinate is non-zero; All Z coordinates are zero
#   . 0D - Zero dimensional: All atomic coordinates are zero
#
sub GetDimensionality {
  my($This) = @_;

  # Is Dimensionality property explicitly set?
  if ($This->HasProperty('Dimensionality')) {
    return $This->GetProperty('Dimensionality');
  }
  my($Atom, @Atoms);

  @Atoms = $This->GetAtoms();
  ATOM: for $Atom (@Atoms) {
      if ($Atom->GetZ() != 0) {
	return '3D';
      }
      if ($Atom->GetX() != 0 || $Atom->GetY() != 0) {
	return '2D';
      }
  }
  return '0D';
}

# Is it a molecule object?
sub IsMolecule ($) {
  my($Object) = @_;

  return _IsMolecule($Object);
}

# Return a string containing vertices, edges and other properties...
sub StringifyMolecule {
  my($This) = @_;
  my($MoleculeString, $ID, $Name, $NumOfAtoms, $NumOfBonds, $MolecularFormula, $NumOfRings, $MolecularWeight, $ExactMass, $FormalCharge, $SpinMultiplicity, $FreeRadicalElectrons, $Charge, $ElementsRef, $ElementsCompositionRef, $ElementalComposition);

  $ID = $This->GetID();
  $Name = $This->GetName();
  $NumOfAtoms = $This->GetNumOfAtoms();
  $NumOfBonds = $This->GetNumOfBonds();

  $NumOfRings = $This->GetNumOfRings();
  if (!defined $NumOfRings) {
    $NumOfRings = 'undefined';
  }

  $MolecularFormula = $This->GetMolecularFormula();

  $MolecularWeight = $This->GetMolecularWeight();
  $MolecularWeight = round($MolecularWeight, 4) + 0;

  $ExactMass = $This->GetExactMass();
  $ExactMass = round($ExactMass, 4) + 0;

  $FormalCharge = $This->GetFormalCharge();
  $Charge = $This->GetCharge();

  $SpinMultiplicity = $This->GetSpinMultiplicity();
  $FreeRadicalElectrons = $This->GetFreeRadicalElectrons();

  ($ElementsRef, $ElementsCompositionRef) = $This->GetElementalComposition();
  $ElementalComposition = 'None';
  if (defined($ElementsRef) && @{$ElementsRef}) {
    $ElementalComposition = "[ " . FormatElementalCompositionInformation($ElementsRef, $ElementsCompositionRef) . " ]";
  }

  $MoleculeString = "Molecule: ID: $ID; Name: \"$Name\"; NumOfAtoms: $NumOfAtoms; NumOfBonds: $NumOfBonds; NumOfRings: $NumOfRings; MolecularFormula: $MolecularFormula; MolecularWeight: $MolecularWeight; ExactMass: $ExactMass; FormalCharge: $FormalCharge; Charge: $Charge; SpinMultiplicity: $SpinMultiplicity; FreeRadicalElectrons: $FreeRadicalElectrons; ElementalComposition: $ElementalComposition";

  return $MoleculeString;
}

# Load appropriate atom data files from <MayaChemTools>/lib directory used by various
# object methods in the current class...
#
sub _LoadMoleculeClassData {
  my($MayaChemToolsLibDir);

  $MayaChemToolsLibDir = GetMayaChemToolsLibDirName();

  # Load and process data for aromaticity models...
  _LoadAromaticityModelsData($MayaChemToolsLibDir);
  _ProcessAromaticityModelsData();
}

#
# Load data for supported aromaticity models...
#
sub _LoadAromaticityModelsData {
  my($MayaChemToolsLibDir) = @_;
  my($DataFile, $Index, $InDelim, $Line, $NumOfCols, $ParameterName, $ParameterValue, $ModelName, @ColLabels, @LineWords, %ParameterNames, %ColIndexToModelName, %SupportedParameterNames);

  %AromaticityModelsDataMap = ();
  %CanonicalAromaticityModelNamesMap = ();

  # File format:
  #
  # "ParameterName","MDLAromaticityModel","TriposAromaticityModel","MMFFAromaticityModel","ChemAxonBasicAromaticityModel","ChemAxonGeneralAromaticityModel","DaylightAromaticityModel","MayaChemToolsAromaticityModel"
  # "AllowHeteroRingAtoms","No","No","Yes","Yes","Yes","Yes","Yes"
  #
  $DataFile = $MayaChemToolsLibDir . "/data/AromaticityModelsData.csv";
  if (! -e "$DataFile") {
    croak "Error: ${ClassName}::_LoadAromaticityModelsData: MayaChemTools package file, $DataFile, is missing: Possible installation problems...";
  }

  # Setup a list of currently supported aromaticity parameters...
  #
  my(@KnownNames);
  @KnownNames = qw(AllowHeteroRingAtoms HeteroRingAtomsList AllowExocyclicDoubleBonds AllowHomoNuclearExocyclicDoubleBonds AllowElectronegativeRingAtomExocyclicDoubleBonds AllowRingAtomFormalCharge AllowHeteroRingAtomFormalCharge MinimumRingSize);

  %SupportedParameterNames = ();
  for $ParameterName (@KnownNames) {
    $SupportedParameterNames{$ParameterName} = $ParameterName;
  }

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

  %ColIndexToModelName = ();

  # Process names of aromaticity models...
  for $Index (1 .. $#ColLabels) {
    $ModelName = $ColLabels[$Index];
    $ModelName =~ s/ //g;

    if (exists $AromaticityModelsDataMap{$ModelName}) {
      croak "Error: ${ClassName}::_LoadAromaticityModelsData: The aromaticity model name, $ModelName, in $DataFile has already exists.\nLine: $Line...";
    }
    %{$AromaticityModelsDataMap{$ModelName}} = ();

    # Cannonicalize aromatic model name by converting into all lowercase...
    $CanonicalAromaticityModelNamesMap{lc($ModelName)} = $ModelName;

    $ColIndexToModelName{$Index} = $ModelName;
  }

  # Process paramater name and their values for specified aromaticity models...
  #
  %ParameterNames = ();
  LINE: while ($Line = GetTextLine(\*DATAFILE)) {
    if ($Line =~ /^#/) {
      next LINE;
    }
    @LineWords = ();
    @LineWords = quotewords($InDelim, 0, $Line);
    if (@LineWords != $NumOfCols) {
      croak "Error: ${ClassName}::_LoadAromaticityModelsData: The number of data fields, @LineWords, in $DataFile must be $NumOfCols.\nLine: $Line...";
    }

    # Process parameter name and values for aromaticity models...
    #
    $ParameterName = $LineWords[0];

    if (!exists $SupportedParameterNames{$ParameterName}) {
      carp "Warning: ${ClassName}::_LoadAromaticityModelsData: The current release of MayaChemTools doesn't support aromaticity model parameter name, $ParameterName, specified in $DataFile. It would be ignore during aromaticity detection.\nLine: $Line...";
    }

    if (exists $ParameterNames{$ParameterName}) {
      carp "Warning: ${ClassName}::_LoadAromaticityModelsData: Ignoring aromaticity model data for parameter name, $ParameterName, in $DataFile. It has already been loaded.\nLine: $Line...";
      next LINE;
    }
    $ParameterNames{$ParameterName} = $ParameterName;

    for $Index (1 .. $#LineWords) {
      $ModelName = $ColIndexToModelName{$Index};
      $ParameterValue = $LineWords[$Index];
      $AromaticityModelsDataMap{$ModelName}{$ParameterName} = $ParameterValue;
    }
  }
  close DATAFILE;
}

# Process already loaded aromaticity model data...
#
sub _ProcessAromaticityModelsData {
  my($ParameterName, $ParameterValue, $ModelName, $NewParameterValue);

  for $ModelName (keys %AromaticityModelsDataMap) {
    for $ParameterName (keys %{$AromaticityModelsDataMap{$ModelName}}) {
      $ParameterValue = $AromaticityModelsDataMap{$ModelName}{$ParameterName};
      $ParameterValue =~ s/ //g;

      VALUE: {
	if ($ParameterValue =~ /^Yes$/i) {
	  $NewParameterValue = 1;
	  last VALUE;
	}
	if ($ParameterValue =~ /^(NA|No)$/i) {
	  $NewParameterValue = 0;
	  last VALUE;
	}
	if ($ParameterValue =~ /^None$/i) {
	  $NewParameterValue = '';
	  last VALUE;
	}
	$NewParameterValue = $ParameterValue;
      }
      $AromaticityModelsDataMap{$ModelName}{$ParameterName} = $NewParameterValue;

      if ($ParameterName =~ /List/i) {
	# Setup a new parameter conatining a reference to a hash for the specified values...
	my($DataMapRefName, $DataValue, %DataMap);

	$DataMapRefName = "${ParameterName}MapRef";

	%DataMap = ();
	for $DataValue (split /\,/, $NewParameterValue) {
	  $DataMap{$DataValue} = $DataValue;
	}
	$AromaticityModelsDataMap{$ModelName}{$DataMapRefName} = \%DataMap;
      }
    }
  }
}

# Is it a molecule object?
sub _IsMolecule {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

Molecule - Molecule class

=head1 SYNOPSIS

use Molecule;

use Molecule qw(:all);

=head1 DESCRIPTION

B<Molecule> class provides the following methods:

new, AddAtom, AddAtoms, AddBond, AddBonds, AddHydrogens, AddPolarHydrogens,
ClearRings, Copy, DeleteAromaticity, DeleteAtom, DeleteAtoms, DeleteBond,
DeleteBonds, DeleteHydrogens, DeletePolarHydrogens, DetectAromaticity,
DetectRings, FormatElementalCompositionInformation, GetAllAtomPaths,
GetAllAtomPathsStartingAt, GetAllAtomPathsStartingAtWithLength,
GetAllAtomPathsStartingAtWithLengthUpto, GetAllAtomPathsWithLength,
GetAllAtomPathsWithLengthUpto, GetAromaticRings, GetAromaticityModel,
GetAtomNeighborhoods, GetAtomNeighborhoodsWithRadiusUpto,
GetAtomNeighborhoodsWithSuccessorAtoms,
GetAtomNeighborhoodsWithSuccessorAtomsAndRadiusUpto, GetAtomPathBonds,
GetAtomPaths, GetAtomPathsBetween, GetAtomPathsStartingAt,
GetAtomPathsStartingAtWithLength, GetAtomPathsStartingAtWithLengthUpto,
GetAtomPathsWithLength, GetAtomPathsWithLengthUpto, GetAtoms, GetBonds, GetCharge,
GetConnectedComponents, GetConnectedComponentsAtoms, GetDimensionality,
GetElementalComposition, GetElementsAndNonElements, GetExactMass, GetFormalCharge,
GetFreeRadicalElectrons, GetFusedAndNonFusedRings, GetLargestConnectedComponent,
GetLargestConnectedComponentAtoms, GetLargestRing, GetMolecularFormula,
GetMolecularWeight, GetNumOfAromaticRings, GetNumOfAtoms, GetNumOfBonds,
GetNumOfConnectedComponents, GetNumOfElementsAndNonElements, GetNumOfHeavyAtoms,
GetNumOfHydrogenAtoms, GetNumOfMissingHydrogenAtoms, GetNumOfNonHydrogenAtoms,
GetNumOfRings, GetNumOfRingsWithEvenSize, GetNumOfRingsWithOddSize,
GetNumOfRingsWithSize, GetNumOfRingsWithSizeGreaterThan,
GetNumOfRingsWithSizeLessThan, GetRingBonds, GetRingBondsFromRings, GetRings,
GetRingsWithEvenSize, GetRingsWithOddSize, GetRingsWithSize,
GetRingsWithSizeGreaterThan, GetRingsWithSizeLessThan, GetSizeOfLargestRing,
GetSizeOfSmallestRing, GetSmallestRing, GetSpinMultiplicity,
GetSupportedAromaticityModels, GetTopologicallySortedAtoms, GetValenceModel,
HasAromaticAtomsInRings, HasAromaticAtomsNotInRings, HasAromaticRings, HasAtom,
HasBond, HasFusedRings, HasNoRings, HasOnlyOneRing, HasRings, IsAromatic,
IsMolecule, IsRingAromatic, IsSupportedAromaticityModel, IsThreeDimensional,
IsTwoDimensional, KeepLargestComponent, KekulizeAromaticAtoms, NewAtom, NewBond,
SetActiveRings, SetAromaticityModel, SetID, SetValenceModel, StringifyMolecule

The following methods can also be used as functions:

FormatElementalCompositionInformation, IsMolecule

B<Molecule> class is derived from B<ObjectProperty> base class which provides methods not explicitly
defined in B<Molecule> or B<ObjectProperty> class using Perl's AUTOLOAD functionality. These methods
are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

=head2 METHODS

=over 4

=item B<new>

    $NewMolecule = new Molecule([%PropertyNameAndValues]);

Using specified I<Atom> property names and values hash, B<new> method creates a new object
and returns a reference to newly created B<Atom> object. By default, the following properties are
initialized:

    ID = SequentialObjectID
    Name = "Molecule <SequentialObjectID>"

Examples:

    $Molecule = new Molecule();

    $WaterMolecule = new Molecule('Name' => 'Water');

    $Oxygen = new Atom('AtomSymbol' => 'O', 'XYZ' => [0, 0, 0]);
    $Hydrogen1 = new Atom('AtomSymbol' => 'H',
                          'XYZ' => [0.7144, 0.4125, 0]);
    $Hydrogen2 = new Atom('AtomSymbol' => 'H',
                          'XYZ' => [1.1208, -0.2959, 0]);
    $WaterMolecule->AddAtoms($Oxygen, $Hydrogen1, $Hydrogen2);

    $Bond1 = new Bond('Atoms' => [$Oxygen, $Hydrogen1],
                      'BondOrder' => 1);
    $Bond2 = new Bond('Atoms' => [$Oxygen, $Hydrogen2],
                      'BondOrder' => 1);
    $WaterMolecule->AddBonds($Bond1, $Bond2);

=item B<AddAtom>

    $Molecule->AddAtom($Atom);

Adds an I<Atom> to a I<Molecule> and returns I<Molecule>.

=item B<AddAtoms>

    $Molecule->AddAtoms(@Atoms);

Adds I<Atoms> to a I<Molecule> and returns I<Molecule>.

=item B<AddBond>

    $Molecule->AddBond($Bond);

Adds a I<Bond> to a I<Molecule> and returns I<Molecule>.

=item B<AddBonds>

    $Molecule->AddBonds(@Bonds);

Adds I<Bonds> to a I<Molecule> and returns I<Molecule>.

=item B<AddHydrogens>

    $NumOfHydrogensAdded = $Molecule->AddHydrogens();

Adds hydrogens to each atom in a I<Molecule> and returns total number of hydrogens
added. The current release of MayaChemTools doesn't assign hydrogen positions.

=item B<AddPolarHydrogens>

    $NumOfHydrogensAdded = $Molecule->AddPolarHydrogens();

Adds hydrogens to each polar atom - N, O, P or S - in a I<Molecule> and returns total
number of polar hydrogens added. The current release of MayaChemTools doesn't
assign hydrogen positions.

=item B<ClearRings>

    $Molecule->ClearRings();

Deletes all rings associated with I<Molecule> and returns I<Molecule>.

=item B<Copy>

    $MoleculeCopy = $Molecule->Copy();

Copies I<Molecule> and its associated data using B<Storable::dclone> and returns a new
B<Molecule> object.

=item B<DeleteAromaticity>

    $Molecule->DeleteAromaticity();

Deletes aromatic property associated with all atoms and bonds in a I<Molecule> and returns
I<Molecule>.

=item B<DeleteAtom>

    $Molecule->DeleteAtom($Atom);

Deletes I<Atom> from a I<Molecule> and returns I<Molecule>.

=item B<DeleteAtoms>

    $Molecule->DeleteAtoms(@Atoms);

Deletes I<Atoms> from a I<Molecule> and returns I<Molecule>.

=item B<DeleteBond>

    $Molecule->DeleteBond($Bond);

Deletes I<Bond> from a I<Molecule> and returns I<Molecule>.

=item B<DeleteBonds>

    $Molecule->DeleteBonds(@Bonds);

Deletes I<Bonds> from a I<Molecule> and returns I<Molecule>.

=item B<DeleteHydrogens>

    $NumOfHydrogensDeleted = $Molecule->DeleteHydrogens();

Removes hydrogens from each atom in a I<Molecule> and returns total number of hydrogens
deleted.

=item B<DeletePolarHydrogens>

    $NumOfHydrogensDeleted = $Molecule->DeletePolarHydrogens();

Removes hydrogens to each polar atom - N, O, P or S - in a I<Molecule> and returns total
number of polar hydrogens deleted.

=item B<DetectAromaticity>

    $Molecule->DetectAromaticity();

Associates I<Aromatic> property to atoms and bonds involved in aromatic rings or ring
systems in a I<Molecule> and returns I<Molecule>.

This method assumes the ring detection has already been perfomed using B<DetectRings>.
And any existing I<Aromatic> property associated with atoms and bonds is deleted before
performing aromaticity detection.

What is aromaticity? [ Ref 124 ] It's in the code of the implementer, did you
say? Agree. The implementation of aromaticity varies widely across different
packages [ Ref 125 ]; additionally, the implementation details are not always
completely available, and it's not possible to figure out the exact implementation
of aromaticity across various packages. Using the publicly available information,
however, one can try to reproduce the available results to the extent possible,
along with parameterizing all the control parameters used to implement different
aromaticity models, and that's exactly what the current release of MayaChemTools
does.

The implementation of aromaticity corresponding to various aromaticity models in
MayaChemTools package is driven by an external CSV file AromaticityModelsData.csv,
which is distributed with the package and is available in lib/data directory. The CSV
files contains names of supported aromaticity models, along with various control
parameters and their values. This file is loaded and processed during instantiation
of Molecule class and data corresponding to specific aromaticity model are used
to detect aromaticity for that model. Any new aromaticity model added to the
aromaticity data file, using different combinations of values for existing control
parameters, would work without any changes to the code; the addition of any new
control parameters, however, requires its implementation in the code used to
calculate number of pi electrons available towards delocalization in a ring or ring
systems.

The current release of MayaChemTools package supports these aromaticity
models: MDLAromaticityModel, TriposAromaticityModel, MMFFAromaticityModel,
ChemAxonBasicAromaticityModel, ChemAxonGeneralAromaticityModel,
DaylightAromaticityModel, MayaChemToolsAromaticityModel.

The current list of control parameters available to detect aromaticity corresponding
to different aromaticity models are: AllowHeteroRingAtoms, HeteroRingAtomsList,
AllowExocyclicDoubleBonds, AllowHomoNuclearExocyclicDoubleBonds,
AllowElectronegativeRingAtomExocyclicDoubleBonds, AllowRingAtomFormalCharge,
AllowHeteroRingAtomFormalCharge, MinimumRingSize. The values for these control
parameters are specified in AromaticityModelsData.csv file.

Although definition of aromaticity differs across various aromaticity models, a ring
or a ring system containing 4n + 2 pi electrons (Huckel's rule) corresponding to
alternate single and double bonds, in general, is considered aromatic.

The available valence free electrons on heterocyclic ring atoms, involved in two single
ring bonds, are also allowed to participate in pi electron delocalizaiton for most of
the supported aromaticity models.

The presence of exocyclic terminal double bond on ring atoms involved in pi electron
delocalization is only allowed for some of the aromaticity models. Additionally, the type
atoms involved in exocyclic terminal double bonds may result in making a ring or ring
system non-aromatic.

For molecules containing fused rings, each fused ring set is considered as one aromatic
system for counting pi electrons to satisfy Huckel's rule; In case of a failure, rings in
fused set are treated individually for aromaticity detection. Additionally, non-fused
rings are handled on their own during aromaticity detection.

=item B<DetectRings>

    $Molecule->DetectRings();

Detects rings in a I<Molecule> and returns I<Molecule>. Ring detection is performed using
B<DetectCycles> method avaible in B<Graph> class which in turn uses methods available
B<Graph::CyclesDetection> class. B<Graph::CyclesDetection> class implements collapsing path graph
[Ref 31] methodology to detect all cycles in a graph.

=item B<FormatElementalCompositionInformation>

    $FormattedInfo = $Molecule->FormatElementalCompositionInformation(
                     $ElementsRef, $ElementCompositionRef,
                     [$Precision]);
    $FormattedInfo = Molecule::FormatElementalCompositionInformation(
                     $ElementsRef, $ElementCompositionRef,
                     [$Precision]);

Using I<ElementsRef> and I<ElementCompositionRef> arrays referneces containg informatio
about elements and their composition, formats elemental composition information and returns
a I<FormattedInfo> string. Defaule I<Precision> value: I<2>.

=item B<GetAromaticityModel>

    $AromaticityModel = $Molecule->GetAromaticityModel();

Returns name of B<AromaticityModel> set for I<Molecule> corresponding to B<AromaticityModel>
property or default model name of B<MayaChemToolsAromaticityModel>.

=item B<GetAllAtomPaths>

    $AtomPathsRef = $Molecule->GetAllAtomPaths([$AllowCycles]);

Returns all paths as a reference to an array containing reference to arrays with path
B<Atom> objects.

Path atoms correspond to to all possible paths for each atom in molecule with all
possible lengths and sharing of bonds in paths traversed. By default, rings are
included in paths. A path containing a ring is terminated at an atom completing the ring.

For molecule without any rings, this method returns the same set of atom paths as
B<GetAtomPaths> method.

=item B<GetAllAtomPathsStartingAt>

    $AtomPathsRef = $Molecule->GetAllAtomPathsStartingAt($StartAtom,
                    [$AllowCycles]);

Returns all atom paths starting from I<StartAtom> as a reference to an array containing
reference to arrays with path B<Atom> objects.

Path atoms atoms correspond to to all possible paths for specified atom in molecule with all
possible lengths and sharing of bonds in paths traversed. By default, rings are
included in paths. A path containing a ring is terminated at an atom completing the ring.

For molecule without any rings, this method returns the same set of atom paths as
B<GetAtomPathsStartingAt>  method.

=item B<GetAllAtomPathsStartingAtWithLength>

    $AtomPathsRef = $Molecule->GetAllAtomPathsStartingAtWithLength(
                    $StartAtom, $Length, [$AllowCycles]);

Returns all atom paths starting from I<StartAtom> with specified I<Length>as a reference
to an array containing reference to arrays with path B<Atom> objects.

Path atoms atoms correspond to to all possible paths for specified atom in molecule with all
possible lengths and sharing of bonds in paths traversed. By default, rings are
included in paths. A path containing a ring is terminated at an atom completing the ring.

For molecule without any rings, this method returns the same set of atom paths as
B<GetAtomPathsStartingAtWithLength>  method.

=item B<GetAllAtomPathsStartingAtWithLengthUpto>

    $AtomPathsRef = $Molecule->GetAllAtomPathsStartingAtWithLengthUpto(
                    $StartAtom, $Length, [$AllowCycles]);

Returns atom paths starting from I<StartAtom> with length up to I<Length> as a reference
to an array containing reference to arrays with path B<Atom> objects.

Path atoms atoms correspond to all possible paths for specified atom in molecule with length
up to a specified length and sharing of bonds in paths traversed. By default, rings are
included in paths. A path containing a ring is terminated at an atom completing the ring.

For molecule without any rings, this method returns the same set of atom paths as
I<GetAtomPathsStartingAtWithLengthUpto> method.

=item B<GetAllAtomPathsWithLength>

    $AtomPathsRef = $Molecule->GetAllAtomPathsWithLength($Length,
                    [$AllowCycles]);

Returns all atom paths with specified I<Length> as a reference to an array containing
reference to arrays with path B<Atom> objects.

Path atoms correspond to to all possible paths for each atom in molecule with length
up to a specified length and sharing of bonds in paths traversed. By default, rings are
included in paths. A path containing a ring is terminated at an atom completing the ring.

For molecule without any rings, this method returns the same set of atom paths as
as I<GetAtomPathsWithLength> method.

=item B<GetAllAtomPathsWithLengthUpto>

    $AtomPathsRef = $Molecule->GetAllAtomPathsWithLengthUpto($Length,
                    [$AllowCycles]);

Returns all atom paths with length up to I<Length> as a reference to an array containing
reference to arrays with path B<Atom> objects.

Path atoms correspond to to all possible paths for each atom in molecule with length
up to a specified length and sharing of bonds in paths traversed. By default, rings are
included in paths. A path containing a ring is terminated at an atom completing the ring.

For molecule without any rings, this method returns the same set of atom paths as
as I<GetAtomPathsWithLengthUpto> method.

=item B<GetAromaticRings>

    @AtomaticRings = $Molecule->GetAromaticRings();

Returns aromatic rings as an array containing references to arrays of ring I<Atom> objects
in a I<Molecule>.

=item B<GetAtomNeighborhoods>

    @Neighborhoods = $Molecule->GetAtomNeighborhoods($StartAtom);

Returns atom neighborhoods around a I<StartAtom> as an array containing references
to arrays with neighborhood I<Atom> objects at possible radii.

=item B<GetAtomNeighborhoodsWithRadiusUpto>

    @Neighborhoods = $Molecule->GetAtomNeighborhoodsWithRadiusUpto($StartAtom,
                     $Radius);

Returns atom neighborhoods around a I<StartAtom> as an array containing references
to arrays with neighborhood I<Atom> objects up to I<Radius>.

=item B<GetAtomNeighborhoodsWithSuccessorAtoms>

    @Neighborhoods = $Molecule->GetAtomNeighborhoodsWithSuccessorAtoms(
                     $StartAtom);

Returns atom neighborhood around a specified I<StartAtom>, along with their successor
connected atoms, collected at all radii as an array containing references to arrays with first
value corresponding to neighborhood atom at a specific radius and second value as reference
to an array containing its successor connected atoms.

For a neighborhood atom at each radius level, the successor connected atoms correspond to the
neighborhood atoms at the next radius level. Consequently, the neighborhood atoms at the last
radius level don't contain any successor atoms which fall outside the range of specified radius.

=item B<GetAtomNeighborhoodsWithSuccessorAtomsAndRadiusUpto>

    @Neighborhoods = $Molecule->GetAtomNeighborhoodsWithSuccessorAtomsAndRadiusUpto(
                     $StartAtom, $Radius);

Returns atom neighborhood around a specified I<StartAtom>, along with their successor
connected atoms, collected upto specified I<Radiud> as an array containing references to arrays
with first value corresponding to neighborhood atom at a specific radius and second value as
reference to an array containing its successor connected atoms.

For a neighborhood atom at each radius level, the successor connected atoms correspond to the
neighborhood atoms at the next radius level. Consequently, the neighborhood atoms at the last
radius level don't contain any successor atoms which fall outside the range of specified radius.

=item B<GetAtomPathBonds>

    $Return = $Molecule->GetAtomPathBonds(@PathAtoms);

Returns an array containing B<Bond> objects corresponding to successive pair of
atoms in I<PathAtoms>

=item B<GetAtomPaths>

    $AtomPathsRef = $Molecule->GetAtomPaths([$AllowCycles]);

Returns all paths as a reference to an array containing reference to arrays with path
B<Atom> objects.

Path atoms correspond to to all possible paths for each atom in molecule with all
possible lengths and no sharing of bonds in paths traversed. By default, rings are
included in paths. A path containing a ring is terminated at an atom completing the ring.

=item B<GetAtomPathsBetween>

    $AtomPathsRef = $Molecule->GetAtomPathsBetween($StartAtom, $EndAtom);

Returns all paths as between I<StartAtom> and I<EndAtom> as a reference to an array
containing reference to arrays with path B<Atom> objects.

For molecules with rings, atom paths array contains may contain two paths.

=item B<GetAtomPathsStartingAt>

    $AtomPathsRef = $Molecule->GetAtomPathsStartingAt($StartAtom, [$AllowCycles]);

Returns paths starting at I<StartAtom> as a reference to an array containing reference to
arrays with path B<Atom> objects.

Path atoms correspond to all possible paths for specified atom in molecule with all
possible lengths and no sharing of bonds in paths traversed. By default, rings are
included in paths. A path containing a ring is terminated at an atom completing the ring.

=item B<GetAtomPathsStartingAtWithLength>

    $AtomPathsRef = $Molecule->GetAtomPathsStartingAtWithLength($StartAtom,
	            $Length, [$AllowCycles]);

Returns paths starting at I<StartAtom> with length  I<Length> as a reference to an array
containing reference to arrays with path B<Atom> objects.

Path atoms correspond to all possible paths for specified atom in molecule with length
upto a specified length and no sharing of bonds in paths traversed. By default, rings are
included in paths. A path containing a ring is terminated at an atom completing the ring.

=item B<GetAtomPathsStartingAtWithLengthUpto>

    $AtomPathsRef = $Molecule->GetAtomPathsStartingAtWithLengthUpto($StartAtom,
	            $Length, [$AllowCycles]);

Returns paths starting at I<StartAtom> with length up to I<Length> as a reference to an array
containing reference to arrays with path B<Atom> objects.

Path atoms correspond to all possible paths for specified atom in molecule with length
upto a specified length and no sharing of bonds in paths traversed. By default, rings are
included in paths. A path containing a ring is terminated at an atom completing the ring.

=item B<GetAtomPathsWithLength>

    $AtomPathsRef = $Molecule->GetAtomPathsWithLength($Length, [$AllowCycles]);

Returns all paths with specified I<Length> as a reference to an array containing reference
to arrays with path B<Atom> objects.

Path atoms correspond to all possible paths for each atom in molecule with length
upto a specified length and no sharing of bonds in paths traversed. By default, rings are
included in paths. A path containing a ring is terminated at an atom completing the ring.

=item B<GetAtomPathsWithLengthUpto>

    $AtomPathsRef = $Molecule->GetAtomPathsWithLengthUpto($Length, [$AllowCycles]);

Returns all paths with length up to I<Length> as a reference to an array containing reference
to arrays with path B<Atom> objects.

Path atoms correspond to all possible paths for each atom in molecule with length
upto a specified length and no sharing of bonds in paths traversed. By default, rings are
included in paths. A path containing a ring is terminated at an atom completing the ring.

=item B<GetAtoms>

    @AllAtoms = $Molecule->GetAtoms();
    @PolarAtoms = $Molecule->GetAtoms('IsPolarAtom');

    $NegateMethodResult = 1;
    @NonHydrogenAtoms = $Molecule->GetAtoms('IsHydrogenAtom',
                        $NegateMethodResult);

    $AtomsCount = $Molecule->GetAtoms();

Returns an array of I<Atoms> in a I<Molecule>. In scalar context,  it returns number of atoms.
Additionally, B<Atoms> array can be filtered by any user specifiable valid B<Atom> class method
and the result of the B<Atom> class method used to filter the atoms can also be negated by
an optional negate results flag as third parameter.

=item B<GetBonds>

    @Bonds = $Molecule->GetBonds();
    $BondsCount = $Molecule->GetBonds();

Returns an array of I<Bonds> in a I<Molecule>. In scalar context,  it returns number of bonds.

=item B<GetCharge>

    $Charge = $Molecule->GetCharge();

Returns net charge on a I<Molecule> using one of the following two methods: explicitly
set B<Charge> property or sum of partial atomic charges on each atom.

=item B<GetConnectedComponents>

    @ConnectedComponents = $Molecule->GetConnectedComponents();

Returns a reference to an array containing I<Molecule> objects corresponding
to connected components sorted in decreasing order of component size in a I<Molecule>.

=item B<GetConnectedComponentsAtoms>

    @ConnectedComponentsAtoms =
      $Molecule->GetConnectedComponentsAtoms();

Returns an array containing references to arrays with I<Atom> objects corresponding to
atoms of connected components sorted in order of component decreasing size in a
I<Molecule>.

=item B<GetDimensionality>

    $Dimensionality = $Molecule->GetDimensionality();

Returns I<Dimensionality> of a I<Molecule> corresponding to explicitly set
I<Dimensionality> property value or by processing atomic.

The I<Dimensionality> value from atomic coordinates is calculated as follows:

    3D - Three dimensional: One of X, Y or Z coordinate is non-zero
    2D - Two dimensional: One of X or Y coordinate is non-zero; All Z
         coordinates are zero
    0D - Zero dimensional: All atomic coordinates are zero

=item B<GetElementalComposition>

    ($ElementsRef, $CompositionRef) =
      $Molecule->GetElementalComposition([$IncludeMissingHydrogens]);

Calculates elemental composition and returns references to arrays containing elements
and their percent composition in a I<Molecule>. By default, missing hydrogens are included
during the calculation.

=item B<GetElementsAndNonElements>

    ($ElementsRef, $NonElementsRef) =
     $Molecule->GetElementsAndNonElements([$IncludeMissingHydrogens]);

Counts elements and non-elements in a I<Molecule> and returns references to hashes
containing element and non-element as hash keys with values corresponding to their
count. By default, missing hydrogens are not added to the element hash.

=item B<GetExactMass>

    $ExactMass = $Molecule->GetExactMass();

Returns exact mass of a I<Molecule> corresponding to sum of exact masses of all
the atoms.

=item B<GetFormalCharge>

    $FormalCharge = $Molecule->GetFormalCharge();

Returns net formal charge on a I<Molecule> using one of the following two methods: explicitly
set B<FormalCharge> property or sum of formal charges on each atom.

B<FormalCharge> is different from B<Charge> property of the molecule which corresponds to
sum of partial atomic charges explicitly set for each atom using a specific methodology.

=item B<GetFreeRadicalElectrons>

    $FreeRadicalElectrons = $Molecule->GetFreeRadicalElectrons();

Returns total number of free radical electrons available in a I<Molecule> using one of the
following two methods: explicitly set B<FreeRadicalElectrons> property or sum of available
free radical electrons on each atom.

=item B<GetFusedAndNonFusedRings>

    ($FusedRingSetRef, $NonFusedRingsRef) =
       $Molecule->GetFusedAndNonFusedRings();

Returns references to array of fused ring sets and non-fused rings in a I<Molecule>. Fused ring sets
array reference contains refernces to arrays of rings corresponding to ring I<Atom> objects;
Non-fused rings array reference contains references to arrays of ring I<Atom> objects.

=item B<GetLargestConnectedComponent>

    $ComponentMolecule = $Molecule->GetLargestConnectedComponent();

Returns a reference to B<Molecule> object corresponding to a largest connected component
in a I<Molecule>.

=item B<GetLargestConnectedComponentAtoms>

    @ComponentAtoms = $Molecule->GetLargestConnectedComponentAtoms();

Returns a reference to an array of B<Atom> objects corresponding to a largest connected
component in a I<Molecule>.

=item B<GetLargestRing>

    @RingAtoms = $Molecule->GetLargestRing();

Returns an array of I<Atoms> objects corresponding to a largest ring in a I<Molecule>.

=item B<GetMolecularFormula>

    $FormulaString = $Molecule->GetMolecularFormula(
                     [$IncludeMissingHydrogens,
                     $IncludeNonElements]);

Returns molecular formula of a I<Molecule> by collecting information about all atoms in
the molecule and composing the formula using Hills ordering system:

    o C shows up first and H follows assuming C is present.
    o All other standard elements are sorted alphanumerically.
    o All other non-stanard atom symbols are also sorted
      alphanumerically and follow standard elements.

Notes:

    o By default, missing hydrogens and nonelements are also included.
    o Elements for disconnected fragments are combined into the same
      formula.
    o Formal charge is also used during compoisiton of molecular formula.

=item B<GetMolecularWeight>

    $MolWeight = $Molecule->GetMolecularWeight();

Returns molecular weight of a I<Molecule> corresponding to sum of atomic weights of all
the atoms.

=item B<GetNumOfAromaticRings>

    $NumOfAromaticRings = $Molecule->GetNumOfAromaticRings();

Returns number of aromatic rings in a I<Molecule>.

=item B<GetNumOfAtoms>

    $NumOfAtoms = $Molecule->GetNumOfAtoms();

Returns number of atoms in a I<Molecule>.

=item B<GetNumOfBonds>

    $NumOfBonds = $Molecule->GetNumOfBonds();

Returns number of bonds in a I<Molecule>.

=item B<GetNumOfConnectedComponents>

    $NumOfComponents = $Molecule->GetNumOfConnectedComponents();

Returns number of connected components in a I<Molecule>.

=item B<GetNumOfElementsAndNonElements>

    ($NumOfElements, $NumOfNonElements) = $Molecule->
                              GetNumOfElementsAndNonElements();
    ($NumOfElements, $NumOfNonElements) = $Molecule->
                   GetNumOfElementsAndNonElements($IncludeMissingHydrogens);

Returns number of elements and non-elements in a I<Molecule>. By default, missing
hydrogens are not added to element count.

=item B<GetNumOfHeavyAtoms>

    $NumOfHeavyAtoms = $Molecule->GetNumOfHeavyAtoms();

Returns number of heavy atoms, non-hydrogen atoms, in a I<Molecule>.

=item B<GetNumOfHydrogenAtoms>

    $NumOfHydrogenAtoms = $Molecule->GetNumOfHydrogenAtoms();

Returns number of hydrogen atoms in a I<Molecule>.

=item B<GetNumOfMissingHydrogenAtoms>

    $NumOfMissingHydrogenAtoms = $Molecule->GetNumOfMissingHydrogenAtoms();

Returns number of hydrogen atoms in a I<Molecule>.

=item B<GetNumOfNonHydrogenAtoms>

    $NumOfNonHydrogenAtoms = $Molecule->GetNumOfNonHydrogenAtoms();

Returns number of non-hydrogen atoms in a I<Molecule>.

=item B<GetNumOfRings>

    $RingCount = $Molecule->GetNumOfRings();

Returns number of rings in a I<Molecule>.

=item B<GetNumOfRingsWithEvenSize>

    $RingCount = $Molecule->GetNumOfRingsWithEvenSize();

Returns number of rings with even size in a I<Molecule>.

=item B<GetNumOfRingsWithOddSize>

    $RingCount = $Molecule->GetNumOfRingsWithOddSize();

Returns number of rings with odd size in a I<Molecule>.

=item B<GetNumOfRingsWithSize>

    $RingCount = $Molecule->GetNumOfRingsWithSize($Size);

Returns number of rings with I<Size> in a I<Molecule>.

=item B<GetNumOfRingsWithSizeGreaterThan>

    $RingCount = $Molecule->GetNumOfRingsWithSizeGreaterThan($Size);

Returns number of rings with size greater than I<Size> in a I<Molecule>.

=item B<GetNumOfRingsWithSizeLessThan>

    $RingCount = $Molecule->GetNumOfRingsWithSizeLessThan($Size);

Returns number of rings with size less than I<Size> in a I<Molecule>.

=item B<GetRingBonds>

    @RingBonds = $Molecule->GetRingBonds(@RingAtoms);

Returns an array of ring B<Bond> objects correponding to an array of ring I<Atoms> in a
I<Molecule>.

=item B<GetRingBondsFromRings>

    @RingBondsSets = $Molecule->GetRingBondsFromRings(@RingAtomsSets);

Returns an array containing references to arrays of ring B<Bond> objects for rings specified
in an array of references to ring I<Atom> objects.

=item B<GetRings>

    @Rings = $Molecule->GetRings();

Returns rings as an array containing references to arrays of ring I<Atom> objects in a I<Molecule>.

=item B<GetRingsWithEvenSize>

    @Rings = $Molecule->GetRingsWithEvenSize();

Returns even size rings as an array containing references to arrays of ring I<Atom> objects in
a I<Molecule>.

=item B<GetRingsWithOddSize>

    @Rings = $Molecule->GetRingsWithOddSize();

Returns odd size rings as an array containing references to arrays of ring I<Atom> objects in
a I<Molecule>.

=item B<GetRingsWithSize>

    @Rings = $Molecule->GetRingsWithSize($Size);

Returns rings with I<Size> as an array containing references to arrays of ring I<Atom> objects in
a I<Molecule>.

=item B<GetRingsWithSizeGreaterThan>

    @Rings = $Molecule->GetRingsWithSizeGreaterThan($Size);

Returns rings with size greater than I<Size> as an array containing references to arrays of
ring I<Atom> objects in a I<Molecule>.

=item B<GetRingsWithSizeLessThan>

    @Rings = $Molecule->GetRingsWithSizeLessThan($Size);

Returns rings with size less than I<Size> as an array containing references to arrays of
ring I<Atom> objects in a I<Molecule>.

=item B<GetSizeOfLargestRing>

    $Size = $Molecule->GetSizeOfLargestRing();

Returns size of the largest ring in a I<Molecule>.

=item B<GetSizeOfSmallestRing>

    $Size = $Molecule->GetSizeOfSmallestRing();

Returns size of the smalles ring in a I<Molecule>.

=item B<GetSmallestRing>

    @RingAtoms = $Molecule->GetSmallestRing();

Returns an array containing I<Atom> objects corresponding to the smallest ring in
a I<Molecule>.

=item B<GetSpinMultiplicity>

    $SpinMultiplicity = $Molecule->GetSpinMultiplicity();

Returns net spin multiplicity of a I<Molecule> using one of the following two methods: explicitly
set B<SpinMultiplicity> property or sum of spin multiplicity on each atom.

=item B<GetSupportedAromaticityModels>

    @SupportedModels = $Molecule->GetSupportedAromaticityModels();

Returns an array containing a list of supported aromaticity models.

=item B<GetValenceModel>

    $ValenceModel = $Molecule->GetValenceModel();

Returns valence model for I<Molecule> using one of the following two methods: explicitly
set B<ValenceModel> property or defaul value of I<InternalValenceModel>.

=item B<GetTopologicallySortedAtoms>

    @SortedAtoms = $Molecule->GetTopologicallySortedAtoms([$StartAtom]);

Returns an array of topologically sorted I<Atom> objects starting from I<StartAtom> or
an arbitrary atom in a I<Molecule>.

=item B<HasAromaticRings>

    $Status = $Molecule->HasAromaticRings();

Returns 1 or 0 based on whether any aromatic ring is present in a I<Molecule>.

=item B<HasAromaticAtomsInRings>

    $Status = $Molecule->HasAromaticAtomsInRings();

Returns 1 or 0 based on whether any aromatic ring atom is present in a I<Molecule>.

=item B<HasAromaticAtomsNotInRings>

    $Status = $Molecule->HasAromaticAtomsNotInRings();

Returns 1 or 0 based on whether any non-ring atom is marked aromatic in a I<Molecule>.

=item B<HasAtom>

    $Status = $Molecule->HasAtom($Atom);

Returns 1 or 0 based on whether I<Atom> is present in a I<Molecule>.

=item B<HasBond>

    $Status = $Molecule->HasBond($Bond);

Returns 1 or 0 based on whether I<Bond> is present in a I<Molecule>.

=item B<HasFusedRings>

    $Status = $Molecule->HasFusedRings();

Returns 1 or 0 based on whether any fused rings set is present in a I<Molecule>.

=item B<HasNoRings>

    $Status = $Molecule->HasNoRings();

Returns 0 or 1 based on whether any ring is present in a I<Molecule>.

=item B<HasOnlyOneRing>

    $Status = $Molecule->HasOnlyOneRing();

Returns 1 or 0 based on whether only one ring is present in a I<Molecule>.

=item B<HasRings>

    $Status = $Molecule->HasRings();

Returns 1 or 0 based on whether rings are present in a I<Molecule>.

=item B<IsAromatic>

    $Status = $Molecule->IsAromatic();

Returns 1 or 0 based on whether I<Molecule> is aromatic.

=item B<IsMolecule>

    $Status = Molecule::IsMolecule();

Returns 1 or 0 based on whether I<Object> is a B<Molecule> object.

=item B<IsRingAromatic>

    $Status = $Molecule->IsRingAromatic(@RingAtoms);

Returns 1 or 0 based on whether all I<RingAtoms> are aromatic.

=item B<IsSupportedAromaticityModel>

    $Status = $Molecule->IsSupportedAromaticityModel($AromaticityModel);
    $Status = Molecule::IsSupportedAromaticityModel($AromaticityModel);

Returns 1 or 0 based on whether specified I<AromaticityModel> is supported.

=item B<IsTwoDimensional>

    $Status = $Molecule->IsTwoDimensional();

Returns 1 or 0 based on whether any atom in I<Molecule> has a non-zero value
for X or Y coordinate and all atoms have zero value for Z coordinates.

=item B<IsThreeDimensional>

    $Status = $Molecule->IsThreeDimensional();

Returns 1 or 0 based on whether any atom in I<Molecule> has a non-zero value
for Z coordinate.

=item B<KeepLargestComponent>

    $Molecule->KeepLargestComponent();

Deletes atoms corresponding to all other connected components Except for the largest
connected component in a I<Molecule> and returns I<Molecule>.

=item B<KekulizeAromaticAtoms>

    $Status = $Molecule->KekulizeAromaticAtoms();

Kekulize marked ring and non-ring aromatic atoms in a molecule and return 1 or 1 based
on whether the kekulization succeeded.

=item B<NewAtom>

    $NewAtom = $Molecule->NewAtom(%AtomPropertyNamesAndValues);

Creates a new atom using I<AtomPropertyNamesAndValues>, add its to I<Molecule>, and returns
new B<Atom> object.

=item B<NewBond>

    $NewBond = $Molecule->NewBond(%BondPropertyNamesAndValues);

Creates a new bond using I<AtomPropertyNamesAndValues>, add its to I<Molecule>, and returns
new B<Bond> object.

=item B<SetActiveRings>

    $Molecule->SetActiveRings($RingsType);

Sets up type of detected ring sets to use during all ring related methods and returns I<Molecule>.
Possible I<RingType> values: I<Independent or All>. By default, I<Independent> ring set is used
during all ring methods.

=item B<SetAromaticityModel>

    $Molecule = $Molecule->SetAromaticityModel($AromaticityModel);

Sets up I<AromaticityModel> property value for I<Molecule> and retrurns I<Molecule>.

=item B<SetValenceModel>

    $Molecule = $Molecule->SetValenceModel(ValenceModel);

Sets up I<ValenceModel> property value for I<Molecule> and retrurns I<Molecule>.

=item B<StringifyMolecule>

    $MoleculeString = $Molecule->StringifyMolecule();

Returns a string containing information about I<Molecule> object

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Atom.pm, Bond.pm, MoleculeFileIO.pm, MolecularFormula.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
