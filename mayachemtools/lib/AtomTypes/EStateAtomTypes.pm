package AtomTypes::EStateAtomTypes;
#
# File: EStateAtomTypes.pm
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
@EXPORT = qw(GetEStateAtomTypesData GetAllPossibleEStateAtomTypes GetAllPossibleEStateNonHydrogenAtomTypes);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, %EStateAtomTypesDataMap);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyEStateAtomTypes';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeEStateAtomTypes();

  $This->_InitializeEStateAtomTypesProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Initialize the data hash. It'll be loaded on demand later...
  %EStateAtomTypesDataMap = ();
}

# Initialize object data...
#
sub _InitializeEStateAtomTypes {
  my($This) = @_;

  # Type of AtomTypes...
  $This->{Type} = 'EState';

  # By default, EState atom types are also assigned to hydrogens...
  $This->{IgnoreHydrogens} = 0;

  return $This;
}

# Initialize object properties...
#
sub _InitializeEStateAtomTypesProperties {
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

# Get EState atom types and associated data loaded from EState data file as
# a reference to hash with the following hash data format:
#
# @{$EStateAtomTypesDataMap{AtomTypes}} - Array of all possible atom types for all atoms
# @{$EStateAtomTypesDataMap{NonHydrogenAtomTypes}} - Array of all possible atom types for non-hydrogen atoms
# @{$EStateAtomTypesDataMap->{ColLabels}} - Array of column labels
# %{$EStateAtomTypesDataMap->{DataCol<Num>}} - Hash keys pair: <DataCol<Num>, AtomType>
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetEStateAtomTypesData {

  # Make sure data is loaded...
  _CheckAndLoadEStateAtomTypesData();

  return \%EStateAtomTypesDataMap;
}

# Get all possible E-state atom types corresponding to hydrogen and non-hydrogen
# atoms as an array reference...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAllPossibleEStateAtomTypes {
  return _GetAllPossibleEStateAtomTypes();
}

# Get all possible E-state atom types corresponding to non-hydrogen atoms
# as an array reference...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAllPossibleEStateNonHydrogenAtomTypes {
  my($NonHydrogensOnly);

  $NonHydrogensOnly = 1;
  return _GetAllPossibleEStateAtomTypes($NonHydrogensOnly);
}

# Get all possible E-state atom types as an array reference...
#
sub _GetAllPossibleEStateAtomTypes {
  my($NonHydrogensOnly) = @_;
  my($EStateAtomTypesDataRef);

  $NonHydrogensOnly = defined $NonHydrogensOnly ? $NonHydrogensOnly : 0;

  $EStateAtomTypesDataRef = GetEStateAtomTypesData();

  return $NonHydrogensOnly ? \@{$EStateAtomTypesDataRef->{NonHydrogenAtomTypes}}: \@{$EStateAtomTypesDataRef->{AtomTypes}};
}

# Assign electrotopological state (E-state) [ Ref 75-78 ] atom types to all atoms...
#
# E-state atom types for various different atom groups [Appendix Table 1 in Ref 76, Appendix III
# in Ref 77 ] are defined using central atom environments indicating its topological and valence state
# along with bonded hydrogens.
#
# The current release of MayaChemTools implements an extended E-state atom assignment
# methodology which is able to assign atom types to any valid non-hydrogen atom in any
# atom group instead of a fixed set of E-state atom types [ Ref 77].
#
# Let:
#   As = Atom symbol corresponding to element symbol
#
#   H<n>   = Number of implicit and explicit hydrogens for atom
#
#   s = Single bond to non-hydrogen atom neighbors or heavy atoms attached to atom
#   s<x> = Symbol s repeated x times to indicate multiple single bonds
#
#   d = Double bond to non-hydrogen atom neighbors or heavy atoms attached to atom
#   d<x> = Symbol d repeated x times to indicate multiple double bonds
#
#   t = Triple bond to non-hydrogen atom neighbors or heavy atoms attached to atom
#   t<x> = Symbol t repeated x times to indicate multiple triple bonds
#
#   a = Aromatic to bond non-hydrogen atom neighbors or heavy atoms attached to atom
#   t<x> = Symbol a repeated x times to indicate multiple aromatic bonds
#
#   p = Plus or positive formal charge
#   m = Minus or negative formal charge
#
# Then:
#
#   AtomType specification corresponds to:
#
#     t<x>d<x>a<x>s<x>AsH<n>p or t<x>d<x>a<x>s<x>AsH<n>m
#
# Notes:
#   . p and n with values of 0 are not shown.
#   . s, d, t, and a bond symbol with values of zero are not shown.
#   . s and d bonds which are also aromatic don't contribute to the count of single and
#     double bonds; instead, aromatic bond count reflect these bonds.
#   . The E-state atom type assignment scheme implemented in the current release of
#     MayaChemToools package supports assignment of atom types to all the periodic tab'e
#     element.
#
# Hydrogen E-state [ Ref 76-77 ] atom type definitions:
#
# HGroup    AtomType
#
#   -OH        HsOH
#   -SH        HsSH
#
#   -NH2       HsNH2
#   >NH        HssNH
#   =NH        HdNH
#   :NH:       HaaNH
#   -NH3+      HsNH3p
#   >NH2+     HssNH2p
#   >NH-+      HsssNHp
#
#   #CH        HtCH
#   =CH2       HdCH2 - H attached to a terminal vinyl group
#   =CH-       HdsCH - H attached a non-terminal vinyl group
#   :CH:       HaaCH
#
#   >CHF       HCHF
#   -CH2F      HCH2F
#   >CHCl      HCHCl
#   -CH2Cl     HCH2Cl
#
#   CHn (saturated)      HCsats   - H attached to sp3 carbon attached to saturated carbon(s)
#   CHn (unsatd.)        HCsatu   - H attached to sp3 carbon attached to unsaturated carbon(s)
#
#   CHn (aromatic)       Havin    - H attached to a non-terminal vinyl group, =CH-, attached to an aromatic carbon
#
#   CHn        Hother    - H attached to any other type of C, N, O or S
#   AHn        Hmisc     - H not attached to C, N, O or  S
#
# Notes:
#   . - : Single bond; = : Double bond; # : Triple bond
#   . Hother atom type capture Hydrogen atom groups not explicitly defined.
#   . HGroup doesn't explicitly corresponds to functional groups
#     . -OH group could be a hydroxyl group or part of carboxylic acid group and so on
#     . -NH2 group could be primary amine or part of an amide group and so on
#
sub AssignAtomTypes {
  my($This) = @_;
  my($Atom);

  ATOM: for $Atom ($This->GetMolecule()->GetAtoms()) {
    if ($Atom->IsHydrogen()) {
      if (!$This->{IgnoreHydrogens}) {
	$This->_AssignAtomTypeToHydrogenAtom($Atom);
      }
      next ATOM;
    }

    # Handle non-hydrogen atoms..
    $This->_AssignAtomTypeToNonHydrogenAtom($Atom);
  }
  return $This;
}

# Assign E-State atom type to non-hydrogen atom...
#
sub _AssignAtomTypeToNonHydrogenAtom {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = $This->_GetAtomTypeForNonHydrogenAtom($Atom);
  $This->SetAtomType($Atom, $AtomType);

  return $This;
}

# Get E-State atom type for non-hydrogen atom...
#
sub _GetAtomTypeForNonHydrogenAtom {
  my($This, $Atom) = @_;
  my($AtomType, $AtomSymbol, $NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds, $CountAromaticBonds, $NumOfHydrogens, $FormalCharge, @EStateAtomInvariants);

  @EStateAtomInvariants = ();

  $AtomSymbol = $Atom->GetAtomicInvariantValue('AS');
  $NumOfHydrogens = $Atom->GetAtomicInvariantValue('H');
  $FormalCharge = $Atom->GetAtomicInvariantValue('FC');

  $CountAromaticBonds = 1;
  ($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds) = $Atom->GetNumOfBondTypesToNonHydrogenAtoms($CountAromaticBonds);

  # Set up E-state atom invariants symbols...
  if ($NumOfTripleBonds > 0) {
    push @EStateAtomInvariants, "t" x $NumOfTripleBonds;
  }
  if ($NumOfDoubleBonds > 0) {
    push @EStateAtomInvariants, "d" x $NumOfDoubleBonds;
  }
  if ($NumOfAromaticBonds > 0) {
    push @EStateAtomInvariants, "a" x $NumOfAromaticBonds;
  }
  if ($NumOfSingleBonds > 0) {
    push @EStateAtomInvariants, "s" x $NumOfSingleBonds;
  }
  push @EStateAtomInvariants, $AtomSymbol;
  if ($NumOfHydrogens > 0) {
    push @EStateAtomInvariants, ($NumOfHydrogens > 1) ? "H${NumOfHydrogens}" : "H";
  }
  if ($FormalCharge > 0) {
    push @EStateAtomInvariants, "p";
  }
  elsif ($FormalCharge < 0) {
    push @EStateAtomInvariants, "m";
  }

  $AtomType = TextUtil::JoinWords(\@EStateAtomInvariants, "", 0);

  return $AtomType;
}

# Assign E-State atom type to hydrogen atom...
#
sub _AssignAtomTypeToHydrogenAtom {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = $This->_GetAtomTypeForHydrogenAtom($Atom);
  $This->SetAtomType($Atom, $AtomType);

  return $This;
}

# Get E-State atom type for hydrogen atom...
#
sub _GetAtomTypeForHydrogenAtom {
  my($This, $Atom) = @_;
  my($AtomType, $AtomNeighbor);

  $AtomType = "Hmisc";

  # Get non-hydrogen atom neighbor...
  $AtomNeighbor = $Atom->GetNonHydrogenNeighborOfHydrogenAtom();
  if (!$AtomNeighbor) {
    return $AtomType;
  }

  ATOMNEIGHBOR: {
    if ($AtomNeighbor->IsCarbon()) {
      $AtomType = $This->_GetAtomTypeForHydrogenAttachedToCarbon($AtomNeighbor);
      last ATOMNEIGHBOR;
    }
    if ($AtomNeighbor->IsNitrogen()) {
      $AtomType = $This->_GetAtomTypeForHydrogenAttachedToNitrogen($AtomNeighbor);
      last ATOMNEIGHBOR;
    }
    if ($AtomNeighbor->IsOxygen()) {
      $AtomType = $This->_GetAtomTypeForHydrogenAttachedToOxygen($AtomNeighbor);
      last ATOMNEIGHBOR;
    }
    if ($AtomNeighbor->IsSulfur()) {
      $AtomType = $This->_GetAtomTypeForHydrogenAttachedToSulfur($AtomNeighbor);
      last ATOMNEIGHBOR;
    }
    $AtomType = "Hmisc";
  }
  return $AtomType;
}

# Get E-state atom type for Hydrogen attached to Carbon...
#
# HGroup    AtomType
#
#   #CH        HtCH
#   =CH2       HdCH2 - H attached to a terminal vinyl group
#   =CH-       HdsCH - H attached a non-terminal vinyl group
#   :CH:       HaaCH
#
#   >CHF       HCHF
#   -CH2F      HCH2F
#   >CHCl      HCHCl
#   -CH2Cl     HCH2Cl
#
#   CHn (saturated)      HCsats   - H attached to sp3 carbon attached to saturated carbon(s)
#   CHn (unsatd.)        HCsatu   - H attached to sp3 carbon attached to unsaturated carbon(s)
#
#   CHn (aromatic)       Havin    - H attached to a non-terminal vinyl group, =CH-, attached to an aromatic carbon
#
#
sub _GetAtomTypeForHydrogenAttachedToCarbon {
  my($This, $CarbonAtom) = @_;
  my($AtomType, $AtomNeighbor, $NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds, $CountAromaticBonds, $NumOfHydrogens, $NumOfFluorines, $NumOfChlorines);

  $AtomType = "Hother";

  $NumOfHydrogens = $CarbonAtom->GetAtomicInvariantValue('H');

  $CountAromaticBonds = 1;
  ($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds) = $CarbonAtom->GetNumOfBondTypesToNonHydrogenAtoms($CountAromaticBonds);

  ($NumOfFluorines,  $NumOfChlorines) = $This->_GetNumOfFluorineAndChlorineNeighbors($CarbonAtom);

  ATOMTYPE: {
    if ($NumOfHydrogens == 1 && $NumOfSingleBonds == 0 && $NumOfDoubleBonds == 0 && $NumOfTripleBonds == 1 && $NumOfAromaticBonds == 0) {
      $AtomType = "HtCH";
      last ATOMTYPE;
    }
    if ($NumOfHydrogens == 2 && $NumOfSingleBonds == 0 && $NumOfDoubleBonds == 1 && $NumOfTripleBonds == 0 && $NumOfAromaticBonds == 0) {
      $AtomType = "HdCH2";
      last ATOMTYPE;
    }
    if ($NumOfHydrogens == 1 && $NumOfSingleBonds == 1 && $NumOfDoubleBonds == 1 && $NumOfTripleBonds == 0 && $NumOfAromaticBonds == 0) {
      if ($This->_IsAttachedToAromaticCarbon($CarbonAtom)) {
	$AtomType = "Havin";
      }
      else {
	$AtomType = "HdsCH";
      }
      last ATOMTYPE;
    }
    if ($NumOfHydrogens == 1 && $NumOfSingleBonds == 0 && $NumOfDoubleBonds == 0 && $NumOfTripleBonds == 0 && $NumOfAromaticBonds == 2) {
      $AtomType = "HaaCH";
      last ATOMTYPE;
    }

    if ($NumOfFluorines == 1 && $NumOfHydrogens == 1 && $NumOfSingleBonds == 3 && $NumOfDoubleBonds == 0 && $NumOfTripleBonds == 0 && $NumOfAromaticBonds == 0) {
      $AtomType = "HCHF";
      last ATOMTYPE;
    }
    if ($NumOfFluorines == 1 && $NumOfHydrogens == 2 && $NumOfSingleBonds == 2 && $NumOfDoubleBonds == 0 && $NumOfTripleBonds == 0 && $NumOfAromaticBonds == 0) {
      $AtomType = "HCH2F";
      last ATOMTYPE;
    }
    if ($NumOfChlorines == 1 && $NumOfHydrogens == 1 && $NumOfSingleBonds == 3 && $NumOfDoubleBonds == 0 && $NumOfTripleBonds == 0 && $NumOfAromaticBonds == 0) {
      $AtomType = "HCHCl";
      last ATOMTYPE;
    }
    if ($NumOfChlorines == 1 && $NumOfHydrogens == 2 && $NumOfSingleBonds == 2 && $NumOfDoubleBonds == 0 && $NumOfTripleBonds == 0 && $NumOfAromaticBonds == 0) {
      $AtomType = "HCH2Cl";
      last ATOMTYPE;
    }
    if ($NumOfDoubleBonds == 0 && $NumOfTripleBonds == 0 && $NumOfAromaticBonds == 0) {
      # H  attached to sp3 carbon...
      $AtomType = $This->_GetAtomTypeForHydrogenAttachedToSaturatedOrUnsaturatedCarbons($CarbonAtom);
      last ATOMTYPE;
    }
    $AtomType = "Hother";
  }

  return $AtomType;
}

# Get number of fluorines and chlorines atached to an atom...
#
sub _GetNumOfFluorineAndChlorineNeighbors {
  my($This, $Atom) = @_;
  my($NumOfFluorines, $NumOfChlorines, $AtomNeighbor);

  $NumOfFluorines = 0; $NumOfChlorines = 0;
  ATOMNEIGHBOR: for $AtomNeighbor ($Atom->GetNonHydrogenAtomNeighbors()) {
    if ($AtomNeighbor->IsFluorine()) {
      $NumOfFluorines++;
      next ATOMNEIGHBOR;
    }
    if ($AtomNeighbor->IsChlorine()) {
      $NumOfChlorines++;
      next ATOMNEIGHBOR;
    }
  }
  return ($NumOfFluorines, $NumOfChlorines);
}

# Get atom type of hydrogen atom attached to a sp3 carbon attached to saturated or
# unsaturated carbons...
#
# HGroup    AtomType
#
#   CHn (saturated)      HCsats   - H attached to sp3 carbon attached to saturated carbon(s)
#   CHn (unsatd.)        HCsatu   - H attached to sp3 carbon attached to unsaturated carbon(s)
#
sub _GetAtomTypeForHydrogenAttachedToSaturatedOrUnsaturatedCarbons {
  my($This, $CarbonAtom) = @_;
  my($AtomType, $AtomNeighbor, @AtomNeighbors);

  $AtomType = "Hother";
  @AtomNeighbors = $CarbonAtom->GetNonHydrogenAtomNeighbors();

  # Make sure all neighbors are carbon atoms...
  for $AtomNeighbor (@AtomNeighbors) {
    if (!$AtomNeighbor->IsCarbon()) {
      return $AtomType;
    }
  }

  $AtomType = "HCsats";
  ATOMNEIGHBOR: for $AtomNeighbor ($CarbonAtom->GetNonHydrogenAtomNeighbors()) {
    if ($AtomNeighbor->IsUnsaturated()) {
      $AtomType = "HCsatu";
      last ATOMNEIGHBOR;
    }
  }

  return $AtomType;
}

# Is vinyl carbon attached to an aromatic carbon?
#
sub _IsAttachedToAromaticCarbon {
  my($This, $CarbonAtom) = @_;
  my($Status, $AtomNeighbor, @AtomNeighbors);

  $Status = 0;

  @AtomNeighbors = $CarbonAtom->GetNonHydrogenAtomNeighbors();
  ATOMNEIGHBOR: for $AtomNeighbor (@AtomNeighbors) {
    if ($AtomNeighbor->IsCarbon() && $AtomNeighbor->IsAromatic()) {
      $Status = 1;
      last ATOMNEIGHBOR;
    }
  }
  return $Status;
}


# Get E-state atom type for Hydrogen attached to Nitrogen...
#
# HGroup    AtomType
#
#   -NH2       HsNH2
#   >NH        HssNH
#   =NH        HdNH
#   :NH:       HaaNH
#   -NH3+      HsNH3p
#   >NH2+     HssNH2p
#   >NH-+      HsssNHp
#
sub _GetAtomTypeForHydrogenAttachedToNitrogen {
  my($This, $NitrogenAtom) = @_;
  my($AtomType, $NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds, $CountAromaticBonds, $NumOfHydrogens, $FormalCharge);

  $AtomType = "Hother";

  $NumOfHydrogens = $NitrogenAtom->GetAtomicInvariantValue('H');
  $FormalCharge = $NitrogenAtom->GetAtomicInvariantValue('FC');

  $CountAromaticBonds = 1;
  ($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds) = $NitrogenAtom->GetNumOfBondTypesToNonHydrogenAtoms($CountAromaticBonds);

  if (!($NumOfTripleBonds == 0)) {
    return $AtomType;
  }

  ATOMTYPE: {
    if ($FormalCharge == 0) {
      if ($NumOfHydrogens == 2 && $NumOfSingleBonds == 1 &&  $NumOfDoubleBonds == 0 && $NumOfAromaticBonds == 0) {
	$AtomType = "HsNH2";
	last ATOMTYPE;
      }
      if ($NumOfHydrogens == 1 && $NumOfSingleBonds == 2 && $NumOfDoubleBonds == 0 && $NumOfAromaticBonds == 0) {
	$AtomType = "HssNH";
	last ATOMTYPE;
      }
      if ($NumOfHydrogens == 1 && $NumOfSingleBonds == 0 && $NumOfDoubleBonds == 1 && $NumOfAromaticBonds == 0) {
	$AtomType = "HdNH";
	last ATOMTYPE;
      }
      if ($NumOfHydrogens == 1 && $NumOfSingleBonds == 0 && $NumOfDoubleBonds == 0 && $NumOfAromaticBonds == 2) {
	$AtomType = "HaaNH";
	last ATOMTYPE;
      }
      $AtomType = "Hother";
      last ATOMTYPE;
    }
    if ($FormalCharge == 1) {
      if ($NumOfHydrogens == 3 && $NumOfSingleBonds == 1 && $NumOfDoubleBonds == 0 && $NumOfAromaticBonds == 0) {
	$AtomType = "HsNH3p";
	last ATOMTYPE;
      }
      if ( $NumOfHydrogens == 2 && $NumOfSingleBonds == 2 && $NumOfDoubleBonds == 0 && $NumOfAromaticBonds == 0) {
	$AtomType = "HssNH2p";
	last ATOMTYPE;
      }
      if ($NumOfHydrogens == 1 && $NumOfSingleBonds == 3 && $NumOfDoubleBonds == 0 && $NumOfAromaticBonds == 0) {
	$AtomType = "HsssNHp";
	last ATOMTYPE;
      }
      $AtomType = "Hother";
      last ATOMTYPE;
    }
    $AtomType = "Hother";
  }

  return $AtomType;
}

# Get E-state atom type for Hydrogen attached to Oxygen...
#
# HGroup    AtomType
#
#   -OH        HsOH
#
sub _GetAtomTypeForHydrogenAttachedToOxygen {
  my($This, $OxygenAtom) = @_;
  my($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds, $CountAromaticBonds, $NumOfHydrogens);

  $NumOfHydrogens = $OxygenAtom->GetAtomicInvariantValue('H');

  $CountAromaticBonds = 1;
  ($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds) = $OxygenAtom->GetNumOfBondTypesToNonHydrogenAtoms($CountAromaticBonds);

  return ($NumOfSingleBonds == 1 && $NumOfHydrogens == 1 && $NumOfDoubleBonds == 0 && $NumOfTripleBonds == 0 && $NumOfAromaticBonds == 0) ? "HsOH" : "Hother";
}

# Get E-state atom type for Hydrogen attached to Sulfur...
#
# HGroup    AtomType
#
#   -SH        HsSH
#
sub _GetAtomTypeForHydrogenAttachedToSulfur {
  my($This, $SulfurAtom) = @_;
  my($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds, $CountAromaticBonds, $NumOfHydrogens);

  $NumOfHydrogens = $SulfurAtom->GetAtomicInvariantValue('H');

  $CountAromaticBonds = 1;
  ($NumOfSingleBonds, $NumOfDoubleBonds, $NumOfTripleBonds, $NumOfAromaticBonds) = $SulfurAtom->GetNumOfBondTypesToNonHydrogenAtoms($CountAromaticBonds);

  return ($NumOfSingleBonds == 1 && $NumOfHydrogens == 1 && $NumOfDoubleBonds == 0 && $NumOfTripleBonds == 0 && $NumOfAromaticBonds == 0) ? "HsSH" : "Hother";
}

# Return a string containg data for EStateAtomTypes object...
#
sub StringifyEStateAtomTypes {
  my($This) = @_;
  my($AtomTypesString);

  # Type of AtomTypes...
  $AtomTypesString = "AtomTypes: $This->{Type}; IgnoreHydrogens: " . ($This->{IgnoreHydrogens} ? "Yes" : "No");

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

# Is it a EStateAtomTypes object?
sub _IsEStateAtomTypes {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Check and load EState atom types data...
#
sub _CheckAndLoadEStateAtomTypesData {

  # Is it already loaded?
  if (exists $EStateAtomTypesDataMap{AtomTypes}) {
    return;
  }

  _LoadEStateAtomTypesData();
}

# Load EState atom types data from the file assuming first column to be atom type symbol..
#
# Format:
#
# "AtomType","AtomGroup"
# "sLi","-Li"
# "ssBe","-Be-"
# "ssssBem",">Be<-2"
#
sub _LoadEStateAtomTypesData {
  my($AtomTypesDataFile, $MayaChemToolsLibDir);

  $MayaChemToolsLibDir = FileUtil::GetMayaChemToolsLibDirName();

  $AtomTypesDataFile =  "$MayaChemToolsLibDir" . "/data/EStateAtomTypes.csv";
  if (! -e "$AtomTypesDataFile") {
    croak "Error: MayaChemTools package file, $AtomTypesDataFile, is missing: Possible installation problems...";
  }

  %EStateAtomTypesDataMap = ();
  AtomTypes::AtomTypes::LoadAtomTypesData($AtomTypesDataFile, \%EStateAtomTypesDataMap);
}

1;

__END__

=head1 NAME

EStateAtomTypes

=head1 SYNOPSIS

use AtomTypes::EStateAtomTypes;

use AtomTypes::EStateAtomTypes qw(:all);

=head1 DESCRIPTION

B<EStateAtomTypes> class provides the following methods:

new, AssignAtomTypes, GetAllPossibleEStateAtomTypes,
GetAllPossibleEStateNonHydrogenAtomTypes, GetEStateAtomTypesData,
StringifyEStateAtomTypes

The following functions are available:

GetAllPossibleEStateAtomTypes,
GetAllPossibleEStateNonHydrogenAtomTypes, GetEStateAtomTypesData

B<EStateAtomTypes> is derived from B<AtomTypes> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<EStateAtomTypes>, B<AtomTypes> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The data file EStateAtomTypes.csv distributed with MayaChemTools release contains
all possible electrotopological state (E-state) [ Ref 75-78 ] atom types.

E-state atom types for various different atom groups [Appendix Table 1 in Ref 76, Appendix III
in Ref 77 ] are defined using central atom environments indicating its topological and valence state
along with bonded hydrogens.

The current release of MayaChemTools implements an extended E-state atom assignment
methodology which is able to assign atom types to any valid non-hydrogen atom in any
atom group instead of a fixed set of E-state atom types [ Ref 77].

Let:

    As = Atom symbol corresponding to element symbol

    H<n>   = Number of implicit and explicit hydrogens for atom

    s = Single bond to non-hydrogen atoms attached to atom
    s<x> = Symbol s repeated x times to indicate multiple single bonds

    d = Double bond to non-hydrogen atoms attached to atom
    d<x> = Symbol d repeated x times to indicate multiple double bonds

    t = Triple bond to non-hydrogen atoms attached to atom
    t<x> = Symbol t repeated x times to indicate multiple triple bonds

    a = Aromatic to bond non-hydrogen atoms attached to atom
    t<x> = Symbol a repeated x times to indicate multiple aromatic bonds

    p = Plus or positive formal charge
    m = Minus or negative formal charge

Then, E-state atom type specification for non-hydrogen or heavy atoms corresponds to:

    t<x>d<x>a<x>s<x>AsH<n>p or t<x>d<x>a<x>s<x>AsH<n>m

 Notes:

    o p and n with values of 0 are not shown.
    o s, d, t, and a bond symbol with values of zero are not shown.
    o s and d bonds which are also aromatic don't contribute to the count
      of single and double bonds; instead, aromatic bond count reflect these
      bonds.

Hydrogen E-state [ Ref 76-77 ] atom type definitions are:

HGroup         AtomType

    -OH        HsOH
    -SH        HsSH

    -NH2       HsNH2
    >NH        HssNH
    =NH        HdNH
    :NH:       HaaNH
    -NH3+      HsNH3p
    >NH2+     HssNH2p
    >NH-+      HsssNHp

    #CH        HtCH
    =CH2       HdCH2 - H attached to a terminal vinyl group
    =CH-       HdsCH - H attached a non-terminal vinyl group
    :CH:       HaaCH

    >CHF       HCHF
    -CH2F      HCH2F
    >CHCl      HCHCl
    -CH2Cl     HCH2Cl

    CHn (saturated)      HCsats - H attached to sp3 carbon attached
                                  to saturated carbon(s)
    CHn (unsatd.)        HCsatu - H attached to sp3 carbon attached
                                  to unsaturated carbon(s)

    CHn (aromatic)       Havin -  H attached to a non-terminal vinyl
                                  group, =CH-, attached to an aromatic carbon

    CHn        Hother    - H attached to any other type of C, N, O or S
    AHn        Hmisc     - H not attached to C, N, O or  S

 Notes:

    o - : Single bond; = : Double bond; # : Triple bond
    o Hother atom type capture Hydrogen atom groups not explicitly defined.
    o HGroup doesn't explicitly corresponds to functional groups
    o -OH group could be a hydroxyl group or part of carboxylic acid group and so on
    o -NH2 group could be primary amine or part of an amide group and so on

Examples of E-state atom types for non-hydrogen or heavy atoms:

    sCH3, dCH2, dsCH, ddC, aasC, sNH2 and so on

=head2 METHODS

=over 4

=item B<new>

    $NewEStateAtomTypes = new AtomTypes::EStateAtomTypes(%NamesAndValues);

Using specified I<EStateAtomTypes> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<EStateAtomTypes>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'EState'
    IgnoreHydrogens = 0

Examples:

    $EStateAtomTypes = new AtomTypes::EStateAtomTypes(
                              'Molecule' => $Molecule,
                              'IgnoreHydrogens' => 0);

=item B<AssignAtomTypes>

    $EStateAtomTypes->AssignAtomTypes();

Assigns E-state atom types to all the atoms in a molecule and returns
I<EStateAtomTypes>.

=item B<GetAllPossibleEStateAtomTypes>

    $AllAtomTypesDataRef = $EStateAtomTypes->
                           GetAllPossibleEStateAtomTypes();
    $AllAtomTypesDataRef = AtomTypes::EStateAtomTypes::
                           GetAllPossibleEStateAtomTypes();

Returns all possible EState atom types corresponding to hydrogen and non-hydrogen
atoms as an array reference.

=item B<GetAllPossibleEStateNonHydrogenAtomTypes>

    $AtomTypesDataRef = $EStateAtomTypes->
                        GetAllPossibleEStateNonHydrogenAtomTypes();
    $AtomTypesDataRef = AtomTypes::EStateAtomTypes::
                        GetAllPossibleEStateNonHydrogenAtomTypes();

Returns all possible EState atom types corresponding to non-hydrogen atoms as
an array reference.

=item B<GetEStateAtomTypesData>

    $AtomTypesDataMapRef = $EStateAtomTypes->GetEStateAtomTypesData();
    $AtomTypesDataMapRef = AtomTypes::EStateAtomTypes::
                           GetEStateAtomTypesData();

Returns EState atom types and associated data loaded from EState data file as
a reference to hash with the following hash data format:

    @{$EStateAtomTypesDataMap{AtomTypes}} - Array of all possible atom
                              types for all atoms
    @{$EStateAtomTypesDataMap{NonHydrogenAtomTypes}} - Array of all
                              possible atom types for non-hydrogen atoms
    @{$EStateAtomTypesDataMap->{ColLabels}} - Array of column labels
    %{$EStateAtomTypesDataMap->{DataCol<Num>}} - Hash keys pair:
                                                 DataCol<Num>, AtomType

=item B<StringifyEStateAtomTypes>

    $String = $EStateAtomTypes->StringifyEStateAtomTypes();

Returns a string containing information about I<EStateAtomTypes> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AtomTypes.pm, AtomicInvariantsAtomTypes.pm, DREIDINGAtomTypes.pm,
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
