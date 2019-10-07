package AtomTypes::UFFAtomTypes;
#
# File: UFFAtomTypes.pm
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
@EXPORT = qw(GetUFFAtomTypesData GetAllPossibleUFFAtomTypes GetAllPossibleUFFNonHydrogenAtomTypes);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, %UFFAtomTypesDataMap);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyUFFAtomTypes';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeUFFAtomTypes();

  $This->_InitializeUFFAtomTypesProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Initialize the data hash. It'll be loaded on demand later...
  %UFFAtomTypesDataMap = ();
}

# Initialize object data...
#
sub _InitializeUFFAtomTypes {
  my($This) = @_;

  # Type of AtomTypes...
  $This->{Type} = 'UFF';

  # By default, UFF atom types are also assigned to hydrogens...
  $This->{IgnoreHydrogens} = 0;

  return $This;
}

# Initialize object properties...
#
sub _InitializeUFFAtomTypesProperties {
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

# Get UFF atom types and associated data loaded from UFF data file as
# a reference to hash with the following hash data format:
#
# @{$UFFAtomTypesDataMap{AtomTypes}} - Array of all possible atom types for all atoms
# @{$UFFAtomTypesDataMap{NonHydrogenAtomTypes}} - Array of all possible atom types for non-hydrogen atoms
# @{$UFFAtomTypesDataMap->{ColLabels}} - Array of column labels
# %{$UFFAtomTypesDataMap->{DataCol<Num>}} - Hash keys pair: <DataCol<Num>, AtomType>
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetUFFAtomTypesData {

  # Make sure data is loaded...
  _CheckAndLoadUFFAtomTypesData();

  return \%UFFAtomTypesDataMap;
}

# Get all possible UFF atom types corresponding to hydrogen and non-hydrogen
# atoms as an array reference...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAllPossibleUFFAtomTypes {
  return _GetAllPossibleUFFAtomTypes();
}

# Get all possible UFF atom types corresponding to non-hydrogen atoms
# as an array reference...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAllPossibleUFFNonHydrogenAtomTypes {
  my($NonHydrogensOnly);

  $NonHydrogensOnly = 1;
  return _GetAllPossibleUFFAtomTypes($NonHydrogensOnly);
}

# Get all possible UFF atom types as an array reference...
#
sub _GetAllPossibleUFFAtomTypes {
  my($NonHydrogensOnly) = @_;
  my($UFFAtomTypesDataRef);

  $NonHydrogensOnly = defined $NonHydrogensOnly ? $NonHydrogensOnly : 0;

  $UFFAtomTypesDataRef = GetUFFAtomTypesData();

  return $NonHydrogensOnly ? \@{$UFFAtomTypesDataRef->{NonHydrogenAtomTypes}}: \@{$UFFAtomTypesDataRef->{AtomTypes}};
}
# Assign UFF [ Ref 81-82 ] atom types to all atoms...
#
# Notes:
#   . Some listed atom types - O_3_z,
#     are not assigned to any atom
#     o 126 UFF atom types are listed for elements with atomic number upto 103
#     o AtomTypes::AtomTypes::UFFAtomTypes.pm module is used to assign UFF atom types
#     o Units:
#         o ValenceBondRadius and NonBondRadius: Angstroms
#         o ValenceAngle: Degrees
#         o NonBondEnergy and SP3TorsionalBarrier: kcal/mol
#     o Five-character mnemonic label for UFF atom types
#         o First two characters correspond to chemical symbol with an underscore as second
#           character for elements with one character symbol
#         o Third character describes hybridization or geometry: 1 - linear; 2 - trigonal; R - resonant;
#           3 = tetrahedral; 4 - square planar; 5 - trigonal bipyramidal; 6 - octahedral
#         o Fourth and fifth characters are used as indicators of alternate parameters: formal oxidation
#           state, bridging hydrogens and so on.
#
#
sub AssignAtomTypes {
  my($This) = @_;
  my($Atom, $AtomType);

  ATOM: for $Atom ($This->GetMolecule()->GetAtoms()) {
    if ($This->{IgnoreHydrogens} && $Atom->IsHydrogen()) {
      next ATOM;
    }
    $AtomType = $This->_GetAtomType($Atom);
    $This->SetAtomType($Atom, $AtomType);
  }
  return $This;
}

# Get UFF atom type for atom...
#
sub _GetAtomType {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = '';

  ATOM: {
    if ($Atom->IsCarbon()) {
      $AtomType = $This->_GetAtomTypeForCarbon($Atom);
      last ATOM;
    }
    if ($Atom->IsNitrogen()) {
      $AtomType = $This->_GetAtomTypeForNitrogen($Atom);
      last ATOM;
    }
    if ($Atom->IsOxygen()) {
      $AtomType = $This->_GetAtomTypeForOxygen($Atom);
      last ATOM;
    }
    if ($Atom->IsPhosphorus()) {
      $AtomType = $This->_GetAtomTypeForPhosphorus($Atom);
      last ATOM;
    }
    if ($Atom->IsSulfur()) {
      $AtomType = $This->_GetAtomTypeForSulfur($Atom);
      last ATOM;
    }
    if ($Atom->IsHydrogen()) {
      $AtomType = $This->_GetAtomTypeForHydrogen($Atom);
      last ATOM;
    }
    $AtomType = $This->_GetAtomTypeForOtherAtoms($Atom);
  }

  return $AtomType;
}

# Get UFF atom type for Carbon atom...
#
sub _GetAtomTypeForCarbon {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'None';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {
    if ($Atom->IsAromatic()) {
      $AtomType = 'C_R';
      last ATOMTYPE;
    }

    # Only single bonds...
    if ($NumOfPiBonds == 0) {
      $AtomType = 'C_3';
      last ATOMTYPE;
    }

    # One double bond...
    if ($NumOfPiBonds == 1) {
      $AtomType = 'C_2';
      last ATOMTYPE;
    }

    # One triple bond or two double bonds...
    if ($NumOfPiBonds == 2) {
      $AtomType = 'C_1';
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForCarbon: UFF atom types for Carbon cann't be assigned...";
  }

  return $AtomType;
}

# Get UFF atom type for Nitrogen atom...
#
sub _GetAtomTypeForNitrogen {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'None';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {
    if ($Atom->IsAromatic()) {
      $AtomType = 'N_R';
      last ATOMTYPE;
    }

    # Only single bonds...
    if ($NumOfPiBonds == 0) {
      $AtomType = 'N_3';
      last ATOMTYPE;
    }

    # One double bond...
    if ($NumOfPiBonds == 1) {
      $AtomType = 'N_2';
      last ATOMTYPE;
    }

    # One triple bond or two double bonds...
    if ($NumOfPiBonds == 2) {
      $AtomType = 'N_1';
      last ATOMTYPE;
    }
    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForNitrogen: UFF atom types for Nitrogen cann't be assigned...";
  }

  return $AtomType;
}

# Get UFF atom type for Oxygen atom...
#
sub _GetAtomTypeForOxygen {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'None';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {
    if ($Atom->IsAromatic()) {
      $AtomType = 'O_R';
      last ATOMTYPE;
    }

    # Only single bonds...
    if ($NumOfPiBonds == 0) {
      $AtomType = 'O_3';
      last ATOMTYPE;
    }

    # One double bond...
    if ($NumOfPiBonds == 1) {
      $AtomType = 'O_2';
      last ATOMTYPE;
    }

    # One triple bond or two double bonds...
    if ($NumOfPiBonds == 2) {
      $AtomType = 'O_1';
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForOxygen: UFF atom types for Oxygen cann't be assigned...";
  }

  return $AtomType;
}

# Get UFF atom type for Phosphorus atom...
#
sub _GetAtomTypeForPhosphorus {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'None';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {
    # Is it a four-coordinated Phosphorus for describing organometallic coordinated phosphines?
    if ($This->_IsFourCoordinatedOrganometallicPhosphorus($Atom)) {
      $AtomType = 'P_3+q';
      last ATOMTYPE;
    }

    # -P(-)-
    if ($NumOfSigmaBonds == 3 && $NumOfPiBonds == 0) {
      $AtomType = 'P_3+3';
      last ATOMTYPE;
    }

    # =P(-)(-)-
    if ($NumOfSigmaBonds == 4 && $NumOfPiBonds == 1) {
      $AtomType = 'P_3+5';
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForPhosphorus: UFF atom types for Phosphorus cann't be assigned...";
  }

  return $AtomType;
}

# Get UFF atom type for Sulfur atom...
#
sub _GetAtomTypeForSulfur {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'None';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {
    if ($Atom->IsAromatic()) {
      $AtomType = 'S_R';
      last ATOMTYPE;
    }

    # -S-
    if ($NumOfSigmaBonds == 2 && $NumOfPiBonds == 0) {
      $AtomType = 'S_3+2';
      last ATOMTYPE;
    }

    # -S(=)-
    if ($NumOfSigmaBonds == 3 && $NumOfPiBonds == 1) {
      $AtomType = 'S_3+4';
      last ATOMTYPE;
    }

    # -S(=)(=)-
    if ($NumOfSigmaBonds == 4 && $NumOfPiBonds == 2) {
      $AtomType = 'S_3+6';
      last ATOMTYPE;
    }

    # S=
    if ($NumOfSigmaBonds == 1 && $NumOfPiBonds == 1) {
      $AtomType = 'S_2';
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForSulfur: UFF atom types for Sulfur cann't be assigned...";
  }
  return $AtomType;
}

# Get UFF atom type for Hydrogen atom...
#
sub _GetAtomTypeForHydrogen {
  my($This, $Atom) = @_;
  my($AtomType);

  if ($Atom->GetNumOfHeavyAtomNeighbors() > 1) {
    # Bridging hydrogen as in B2H6
    $AtomType = 'H___b';
  }
  else {
    $AtomType = 'H_';
  }

  return $AtomType;
}

# Get UFF atom type for atoms other than Carbon, Nitrogen, Oxygen, Phosporus,
# Sulfur and Hydrogen...
#
sub _GetAtomTypeForOtherAtoms {
  my($This, $Atom) = @_;
  my($AtomType, $AtomicNumber, $AtomSymbol, $GroupNumber, $MethodName);

  $AtomType = 'None';

  $AtomicNumber = $Atom->GetAtomicNumber();
  $AtomSymbol = $Atom->GetAtomSymbol();
  $GroupNumber = $Atom->GetGroupNumber();

  ATOMTYPE: {
    # Get atom types for atoms in a valid periodic table group number...
    if (defined($GroupNumber) && $GroupNumber) {
      $MethodName = "_GetAtomTypeForOtherAtomsInGroupNumber${GroupNumber}";
      $AtomType = $This->$MethodName($Atom);
      last ATOMTYPE;
    }

    # Get atom types for Lanthanidies...
    if ($AtomicNumber >= 57 && $AtomicNumber <= 71) {
      $AtomType = $This->_GetAtomTypeForOtherAtomsInLanthanoidGroup($Atom);
      last ATOMTYPE;
    }

    # Get atom types for Actinides...
    if ($AtomicNumber >= 89 && $AtomicNumber <= 103) {
      $AtomType = $This->_GetAtomTypeForOtherAtomsInActinoidGroup($Atom);
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForOtherAtoms: UFF atom types for atom, $AtomSymbol, with atomic number, $AtomicNumber, cann't be assigned...";
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 1...
#
# Group number 1 contains: H, Li, Na, K, Rb, Cs, Fr
#
# And corresponding UFF atom types are: Li, Na, K_, Rb, Cs, Fr
#
# Notes:
#   . This method doesn't assign UFF atom type for H.
#
sub _GetAtomTypeForOtherAtomsInGroupNumber1 {
  my($This, $Atom) = @_;
  my($AtomType, $AtomSymbol);

  $AtomSymbol = $Atom->GetAtomSymbol();
  $AtomType = (length($AtomSymbol) == 1) ? "${AtomSymbol}_" : $AtomSymbol;

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 2...
#
# Group number 2 contains: Be, Mg, Ca, Sr, Ba, Ra
#
# And corresponding UFF atom types are: Be3+2, Mg3+2, Ca6+2, Sr6+2, Ba6+2, Ra6+2
#
# Notes:
#   . Although the number of valence electrons is two, the tetrahedral and octahedral
#     geometry is attributed to coordination bonds from other atoms.
#
sub _GetAtomTypeForOtherAtomsInGroupNumber2 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber2';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^(Be|Mg)$/) {
      $AtomType = "${AtomSymbol}3+2";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'tetrahedral', 4, 2);
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^(Ca|Sr|Ba|Ra)$/) {
      $AtomType = "${AtomSymbol}6+2";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'octahedral', 6, 2);
      last ATOMSYMBOL;
    }
    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 3...
#
# Group number 3 contains: Sc, Y, Lu, Lr
#
# And corresponding UFF atom types are: Sc3+3, Y_3+3, Lu6+3, Lr6+3
#
sub _GetAtomTypeForOtherAtomsInGroupNumber3 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber3';

  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  $AtomSymbol = $Atom->GetAtomSymbol();

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^(Sc|Y)$/) {
      $AtomType = (length($AtomSymbol) == 1) ? "${AtomSymbol}_3+3" : "${AtomSymbol}3+3";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'tetrahedral', 4, 3);
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^(Lu|Lr)$/) {
      $AtomType = "${AtomSymbol}6+3";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'octahedral', 6, 3);
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 4...
#
# Group number 4 contains: Ti, Zr, Hf, Rf
#
# And corresponding UFF atom types are: Ti3+4, Ti3+6, Zr3+4, Hf3+4
#
# Notes:
#   . No UFF atom type for Rf
#
sub _GetAtomTypeForOtherAtomsInGroupNumber4 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber4';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^Ti$/) {
      TI: {
	if ($NumOfNeighbors == 4 && $FormalOxidationState == 4) {
	  $AtomType = "Ti3+4";
	  $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'tetrahedral', 4, 4);
	  last TI;
	}

	if ($NumOfNeighbors == 4 && $FormalOxidationState == 6) {
	  $AtomType = "Ti3+6";
	  $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'tetrahedral', 4, 6);
	  last TI;
	}

	# Assign default value...
	$AtomType = "Ti3+6";
	$This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'tetrahedral', 4, 6);
      }
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^(Zr|Hf)$/) {
      $AtomType = "${AtomSymbol}3+4";
      $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'tetrahedral', 4, 4);
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 5...
#
# Group number 5 contains: V, Nb, Ta, Db
#
# And corresponding UFF atom types are: V_3+5, Nb3+5, Ta3+5
#
# Notes:
#   . No UFF atom type for Db
#
sub _GetAtomTypeForOtherAtomsInGroupNumber5 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber5';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^(V|Nb|Ta)$/) {
      $AtomType = (length($AtomSymbol) == 1) ? "${AtomSymbol}_3+5" : "${AtomSymbol}3+5";
      $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'tetrahedral', 4, 5);
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 6...
#
# Group number 6 contains: Cr, Mo, W, Sg
#
# And corresponding UFF atom types are: Cr6+3, Mo6+6, Mo3+6, W_6+6, W_3+4, W_3+6
#
# Notes:
#   . No UFF atom type for Sg
#
sub _GetAtomTypeForOtherAtomsInGroupNumber6 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber6';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^Cr$/) {
      $AtomType = "Cr6+3";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'octahedral', 6, 3);
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^Mo$/) {
      MO: {
	if ($NumOfNeighbors == 6 && $FormalOxidationState == 6) {
	  $AtomType = "Mo6+6";
	  $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'octahedral', 6, 6);
	  last MO;
	}

	if ($NumOfNeighbors == 4 && $FormalOxidationState == 6) {
	  $AtomType = "Mo3+6";
	  $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'tetrahedral', 4, 6);
	  last MO;
	}

	# Assign default value...
	$AtomType = "Mo6+6";
	$This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'octahedral', 6, 6);
      }
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^W$/) {
      W: {
	if ($NumOfNeighbors == 4 && $FormalOxidationState == 4) {
	  $AtomType = "W_3+4";
	  $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'tetrahedral', 4, 4);
	  last W;
	}

	if ($NumOfNeighbors == 4 && $FormalOxidationState == 6) {
	  $AtomType = "W_3+6";
	  $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'tetrahedral', 4, 6);
	  last W;
	}

	if ($NumOfNeighbors == 6 && $FormalOxidationState == 6) {
	  $AtomType = "W_6+6";
	  $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'octahedral', 6, 6);
	  last W;
	}

	# Assign default value...
	$AtomType = "W_6+6";
	$This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'octahedral', 6, 6);
      }
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 7...
#
# Group number 7 contains: Mn, Tc, Re, Bh
#
# And corresponding UFF atom types are: Mn6+2, Tc6+5, Re6+5, Re3+7
#
# Notes:
#   . No UFF atom type for Bh
#
sub _GetAtomTypeForOtherAtomsInGroupNumber7 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber7';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^Mn$/) {
      $AtomType = "Mn6+2";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'octahedral', 6, 2);
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^Tc$/) {
      $AtomType = "Tc6+5";
      $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'octahedral', 6, 6);
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^Re$/) {
      RE: {
	if ($NumOfNeighbors == 6) {
	  $AtomType = "Re6+5";
	  $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'octahedral', 6, 6);
	  last RE;
	}

	if ($NumOfNeighbors == 4 && $FormalOxidationState == 7) {
	  $AtomType = "Re3+7";
	  $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'tetrahedral', 4, 7);
	  last RE;
	}

	# Assign default value...
	$AtomType = "Re6+5";
	$This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'octahedral', 6, 6);
      }
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 8...
#
# Group number 8 contains: Fe, Ru, Os, Hs
#
# And corresponding UFF atom types are: Fe6+2, Ru6+2, Ru6+3, Os6+6
#
# Notes:
#   . No UFF atom type for Hs
#
sub _GetAtomTypeForOtherAtomsInGroupNumber8 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber8';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^Fe$/) {
      $AtomType = "Fe6+2";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'octahedral', 6, 2);
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^Ru$/) {
      RU: {
	if ($NumOfNeighbors == 6 && $FormalCharge == 2) {
	  $AtomType = "Ru6+2";
	  $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'octahedral', 6, 2);
	  last RU;
	}

	if ($NumOfNeighbors == 6 && $FormalCharge == 3) {
	  $AtomType = "Ru6+3";
	  $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'octahedral', 6, 3);
	  last RU;
	}

	# Assign default value...
	$AtomType = "Ru6+3";
	$This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'octahedral', 6, 2);
      }
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^Os$/) {
      $AtomType = "Os6+6";
      $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'octahedral', 6, 6);
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 9...
#
# Group number 9 contains: Co, Rh, Ir, Mt
#
# And corresponding UFF atom types are: Co6+3, Rh6+3, Ir6+3
#
# Notes:
#   . No UFF atom type for Mt
#
sub _GetAtomTypeForOtherAtomsInGroupNumber9 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber9';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^(Co|Rh|Ir)$/) {
      $AtomType = "${AtomSymbol}6+3";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'octahedral', 6, 3);
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 10...
#
# Group number 10 contains: Ni, Pd, Pt
#
# And corresponding UFF atom types are: Ni4+2, Pd4+2, Pt4+2
#
sub _GetAtomTypeForOtherAtomsInGroupNumber10 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber10';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^(Ni|Pd|Pt)$/) {
      $AtomType = "${AtomSymbol}4+2";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'planar', 4, 2);
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 11...
#
# Group number 11 contains: Cu, Ag, Au
#
# And corresponding UFF atom types are: Cu3+1, Ag1+1, Au4+3
#
sub _GetAtomTypeForOtherAtomsInGroupNumber11 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber11';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^Cu$/) {
      $AtomType = "Cu3+1";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'tetrahedral', 4, 1);
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^Ag$/) {
      $AtomType = "Ag1+1";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'linear', 1, 1);
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^Au$/) {
      $AtomType = "Au4+3";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'planar', 4, 3);
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 102..
#
# Group number 12 contains: Zn, Cd, Hg
#
# And corresponding UFF atom types are: Zn3+2, Cd3+2, Hg1+2
#
sub _GetAtomTypeForOtherAtomsInGroupNumber12 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber12';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^(Zn|Cd)$/) {
      $AtomType = "${AtomSymbol}3+2";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'tetrahedral', 4, 2);
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^Hg$/) {
      $AtomType = "Hg1+2";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'linear', 1, 2);
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 13...
#
# Group number 13 contains: B, Al, Ga, In, Tl
#
# And corresponding UFF atom types are: B_3, B_2, Al3, Ga3+3, In3+3, Tl3+3
#
sub _GetAtomTypeForOtherAtomsInGroupNumber13 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber13';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^B$/) {
      B: {
	if ($NumOfNeighbors == 4) {
	  $AtomType = "B_3";
	  $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'tetrahedral', 4, 0);
	  last B;
	}

	if ($NumOfNeighbors == 3) {
	  $AtomType = "B_2";
	  $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'trigonal', 3, 0);
	  last B;
	}

	# Assign default value...
	$AtomType = "B_2";
	$This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'trigonal', 3, 0);
      }

      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^Al$/) {
      $AtomType = "Al3";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'tetrahedral', 4, 0);
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^(Ga|In|Tl)$/) {
      $AtomType = "${AtomSymbol}3+3";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'tetrahedral', 4, 3);
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 14...
#
# Group number 14 contains: C, Si, Ge, Sn, Pb
#
# And corresponding UFF atom types are: Si3, Ge3, Sn3, Pb3
#
# Notes:
#   . This method doesn't assign UFF atom type for C.
#
sub _GetAtomTypeForOtherAtomsInGroupNumber14 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber14';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^(Si|Ge|Sn|Pb)$/) {
      $AtomType = "${AtomSymbol}3";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'tetrahedral', 4, 0);
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 15...
#
# Group number 15 contains: N, P, As, Sb, Bi
#
# And corresponding UFF atom types are: As3+3, Sb3+3, Bi3+3
#
# Notes:
#   . This method doesn't assign UFF atom type for N and P.
#
sub _GetAtomTypeForOtherAtomsInGroupNumber15 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber15';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^(As|Sb|Bi)$/) {
      $AtomType = "${AtomSymbol}3+3";
      $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'tetrahedral', 3, 3);
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 16...
#
# Group number 16 contains: O, S, Se, Te, Po
#
# And corresponding UFF atom types are: Se3+2, Te3+2, Po3+2
#
# Notes:
#   . This method doesn't assign UFF atom type for O and S.
#
sub _GetAtomTypeForOtherAtomsInGroupNumber16 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInGroupNumber16';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^(Se|Te|Po)$/) {
      $AtomType = "${AtomSymbol}3+2";
      $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'tetrahedral', 2, 2);
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 17...
#
# Group number 17 contains: F, Cl, Br, I, At
#
# And corresponding UFF atom types are: F_, Cl, Br, I_, At
#
sub _GetAtomTypeForOtherAtomsInGroupNumber17 {
  my($This, $Atom) = @_;
  my($AtomType, $AtomSymbol);

  $AtomSymbol = $Atom->GetAtomSymbol();
  $AtomType = (length($AtomSymbol) == 1) ? "${AtomSymbol}_" : $AtomSymbol;

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group number 18...
#
# Group number 18 contains: He, Ne, Ar, Kr, Xe, Rn
#
# And corresponding UFF atom types are: He4+4, Ne4+4, Ar4+4, Kr4+4, Xe4+4, Rn4+4
#
sub _GetAtomTypeForOtherAtomsInGroupNumber18 {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType);

  $AtomType = 'None';

  $AtomSymbol = $Atom->GetAtomSymbol();
  $AtomType = "${AtomSymbol}4+4";

  return $AtomType;
}

# Get UFF atom type for atoms in periodic table group name Lanthanoid...
#
# Group name Lanthanoid contains: La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu
#
# And corresponding UFF atom types are: La3+3, Ce6+3, Pr6+3, Nd6+3, Pm6+3, Sm6+3,
# Eu6+3, Gd6+3, Tb6+3, Dy6+3, Ho6+3, Er6+3, Tm6+3, Yb6+3, Lu6+3
#
sub _GetAtomTypeForOtherAtomsInLanthanoidGroup {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInLanthanoidGroup';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^La$/) {
      $AtomType = "La3+3";
      $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'tetrahedral', 4, 3);
      last ATOMSYMBOL;
    }

    $AtomType = "${AtomSymbol}6+3";
    $This->_CheckGeometryAndFormalChargeMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalCharge, 'octahedral', 6, 3);
  }

  return $AtomType;

}

# Get UFF atom type for atoms in periodic table group name Actinoid...
#
# Group name Actinoid contains: Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr
#
# And corresponding UFF atom types are: Ac6+3, Th6+4, Pa6+4, U_6+4, Np6+4, Pu6+4,
# Am6+4, Cm6+3, Bk6+3, Cf6+3, Es6+3, Fm6+3, Md6+3, No6+3, Lr6+3
#
sub _GetAtomTypeForOtherAtomsInActinoidGroup {
  my($This, $Atom) = @_;
  my($AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, $FormalCharge, $MethodName);

  $AtomType = 'None';
  $MethodName = '_GetAtomTypeForOtherAtomsInActinoidGroup';

  $AtomSymbol = $Atom->GetAtomSymbol();
  ($NumOfNeighbors, $FormalOxidationState, $FormalCharge) = $This->_GetAtomEnvironmentInfoForUFFAtomTypes($Atom);

  ATOMSYMBOL: {
    if ($AtomSymbol =~ /^(Ac|Cm|Bk|Cf|Es|Fm|Md|No|Lr)$/) {
      $AtomType = "${AtomSymbol}6+3";
      $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'octahedral', 6, 3);
      last ATOMSYMBOL;
    }

    if ($AtomSymbol =~ /^(Th|Pa|U|Np|Pu|Am)$/) {
      $AtomType = length($AtomSymbol) == 1 ? "${AtomSymbol}_6+4" : "${AtomSymbol}6+4";
      $This->_CheckGeometryAndOxidationStateMismatch($MethodName, $AtomSymbol, $AtomType, $NumOfNeighbors, $FormalOxidationState, 'octahedral', 6, 4);
      last ATOMSYMBOL;
    }

    $AtomType = 'None';
  }

  return $AtomType;
}


# Is it a four-coordinated Phosphorus for describing organometallic coordinated phosphines?
#
sub _IsFourCoordinatedOrganometallicPhosphorus {
  my($This, $Atom) = @_;
  my($AtomNeighbor, $NumOfNeighbors, $MetalNeighborFound, @AtomNeighbors);

  @AtomNeighbors = $Atom->GetHeavyAtomNeighbors();

  # Is it attached to a metallic atom?
  $MetalNeighborFound = 0;
  NEIGHBOR: for $AtomNeighbor (@AtomNeighbors) {
    if ($AtomNeighbor->IsMetallic()) {
      $MetalNeighborFound = 1;
      last NEIGHBOR;
    }
  }

  if (!$MetalNeighborFound) {
    return 0;
  }

  # Is it four coordinated Phosphorus?
  $NumOfNeighbors = scalar @AtomNeighbors;
  if ($NumOfNeighbors <= 4) {
    # As long as total number of heavy atom neighbors, including attached
    # metal atom is <= 4, missing hydrogens would make it a tetra coordinated
    # Phosphorous...
    return 1;
  }

  return 0;
}

# Get UFF atom environment information, number of neighbors and formal oxidatiion state,
#  for assigning UFF atom types...
#
sub _GetAtomEnvironmentInfoForUFFAtomTypes {
  my($This, $Atom) = @_;
  my($NumOfNeighbors, $NumOfHydrogens, $FormalOxidationState, $FormalCharge);

  $NumOfHydrogens = $Atom->GetNumOfMissingHydrogens() + $Atom->GetExplicitHydrogens();

  # Total number of neighbor atoms...
  $NumOfNeighbors = $Atom->GetNumOfNonHydrogenAtomNeighbors() + $NumOfHydrogens;

  # UFF formal oxidation state appears to just the sum of bond orders to all attched non-hyrdogen
  # atoms and all hydrogen atoms...
  #
  $FormalOxidationState = $Atom->GetSumOfBondOrdersToNonHydrogenAtoms() + $NumOfHydrogens;

  # Any explicit formal charge...
  $FormalCharge = $Atom->GetFormalCharge();

  return ($NumOfNeighbors, $FormalOxidationState, $FormalCharge);
}

# Check and warn about any geometry and/or formal charge mismatch...
#
sub _CheckGeometryAndFormalChargeMismatch {
  my($This, $CallingMethod, $AtomSymbol, $AssignedAtomType, $NumOfNeighbors, $FormalCharge, $ExpectedGeometryType, $ExpectedNumOfNeighbors, $ExpectedFormalCharge) = @_;
  my($MsgHeader);

  $MsgHeader = "Warning: ${ClassName}->${CallingMethod}:_CheckGeometryAndFormalChargeMismatch";

  if ($NumOfNeighbors !=$ExpectedNumOfNeighbors  && $FormalCharge != $ExpectedFormalCharge) {
    carp "\n${MsgHeader}: UFF atom type for $AtomSymbol corresponding to $ExpectedGeometryType geometry and formal charge $ExpectedFormalCharge cann't be assigned; Number of neighbors, $NumOfNeighbors, is different from $ExpectedNumOfNeighbors and formal charge, $FormalCharge, is different from $ExpectedFormalCharge. Default UFF atom type, $AssignedAtomType, has been assigned...";
  }
  elsif ($NumOfNeighbors !=$ExpectedNumOfNeighbors) {
    carp "\n${MsgHeader}: UFF atom type for $AtomSymbol corresponding to $ExpectedGeometryType geometry and formal charge $ExpectedFormalCharge cann't be assigned; Number of neighbors, $NumOfNeighbors, is different from $ExpectedNumOfNeighbors. Default UFF atom type, $AssignedAtomType, has been assigned...";
  }
  elsif ($FormalCharge != $ExpectedFormalCharge) {
    carp "\n${MsgHeader}: UFF atom type for $AtomSymbol corresponding to $ExpectedGeometryType geometry and formal charge $ExpectedFormalCharge cann't be assigned; Formal charge, $FormalCharge, is different from $ExpectedFormalCharge. Default UFF atom type, $AssignedAtomType, has been assigned...";
  }
}

# Check and warn about any geometry and/or formal oxidation state mismatch...
#
sub _CheckGeometryAndOxidationStateMismatch {
  my($This, $CallingMethod, $AtomSymbol, $AssignedAtomType, $NumOfNeighbors, $FormalOxidationState, $ExpectedGeometryType, $ExpectedNumOfNeighbors, $ExpectedFormalOxidationState) = @_;
  my($MsgHeader);

  $MsgHeader = "Warning: ${ClassName}->${CallingMethod}:_CheckGeometryAndOxidationStateMismatch";

  if ($NumOfNeighbors !=$ExpectedNumOfNeighbors  && $FormalOxidationState != $ExpectedFormalOxidationState) {
    carp "\n${MsgHeader}: UFF atom type for $AtomSymbol corresponding to $ExpectedGeometryType geometry and formal oxidation state $ExpectedFormalOxidationState cann't be assigned; Number of neighbors, $NumOfNeighbors, is different from $ExpectedNumOfNeighbors and formal oxidation state, $FormalOxidationState, is different from $ExpectedFormalOxidationState. Default UFF atom type, $AssignedAtomType, has been assigned...";
  }
  elsif ($NumOfNeighbors !=$ExpectedNumOfNeighbors) {
    carp "\n${MsgHeader}: UFF atom type for $AtomSymbol corresponding to $ExpectedGeometryType geometry and formal oxidation state $ExpectedFormalOxidationState cann't be assigned; Number of neighbors, $NumOfNeighbors, is different from $ExpectedNumOfNeighbors. Default UFF atom type, $AssignedAtomType, has been assigned...";
  }
  elsif ($FormalOxidationState != $ExpectedFormalOxidationState) {
    carp "\n${MsgHeader}: UFF atom type for $AtomSymbol corresponding to $ExpectedGeometryType geometry and formal oxidation state $ExpectedFormalOxidationState cann't be assigned; Formal oxidation state, $FormalOxidationState, is different from $ExpectedFormalOxidationState. Default UFF atom type, $AssignedAtomType, has been assigned...";
  }
}


# Return a string containg data for UFFAtomTypes object...
#
sub StringifyUFFAtomTypes {
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

# Is it a UFFAtomTypes object?
sub _IsUFFAtomTypes {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Check and load UFF atom types data...
#
sub _CheckAndLoadUFFAtomTypesData {

  # Is it already loaded?
  if (exists $UFFAtomTypesDataMap{AtomTypes}) {
    return;
  }

  _LoadUFFAtomTypesData();
}

# Load UFF atom types data from the file assuming first column to be atom type symbol..
#
# Format:
#
# "AtomType","Mass","ValenceBondRadius","ValenceAngle","NonBondRadius","NonBondEnergy","NonBondScale","EffectiveCharge","SP3TorsionalBarrier"
# "H_","1.0080","0.354","180.000","2.886","0.044","12.000","0.733","0.000"
# "C_3","12.0110","0.757","109.471","3.851","0.105","12.730","1.967","2.119"
# "C_R","12.0110","0.729","120.000","3.851","0.105","12.730","1.967","0.000"
# "C_2","12.0110","0.732","120.000","3.851","0.105","12.730","1.967","0.000"
# "C_1","12.0110","0.711","180.000","3.851","0.105","12.730","1.967","0.000"
#
sub _LoadUFFAtomTypesData {
  my($AtomTypesDataFile, $MayaChemToolsLibDir);

  $MayaChemToolsLibDir = FileUtil::GetMayaChemToolsLibDirName();

  $AtomTypesDataFile =  "$MayaChemToolsLibDir" . "/data/UFFAtomTypes.csv";
  if (! -e "$AtomTypesDataFile") {
    croak "Error: MayaChemTools package file, $AtomTypesDataFile, is missing: Possible installation problems...";
  }

  %UFFAtomTypesDataMap = ();
  AtomTypes::AtomTypes::LoadAtomTypesData($AtomTypesDataFile, \%UFFAtomTypesDataMap);
}

1;

__END__

=head1 NAME

UFFAtomTypes

=head1 SYNOPSIS

use AtomTypes::UFFAtomTypes;

use AtomTypes::UFFAtomTypes qw(:all);

=head1 DESCRIPTION

B<UFFAtomTypes> class provides the following methods:

new, AssignAtomTypes, GetAllPossibleUFFAtomTypes,
GetAllPossibleUFFNonHydrogenAtomTypes, GetUFFAtomTypesData, StringifyUFFAtomTypes

The following functions are available:

GetAllPossibleUFFAtomTypes,
GetAllPossibleUFFNonHydrogenAtomTypes, GetUFFAtomTypesData

B<UFFAtomTypes> is derived from B<AtomTypes> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<UFFAtomTypes>, B<AtomTypes> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The data file UFFAomTypes.csv distributed with MayaChemTools release contains
all possible UFF [ Ref 81-82 ] atom types.

Format of a Five-character mnemonic label used for UFF atom types:

    o First two characters correspond to chemical symbol with an underscore
      as second character for elements with one character symbol
    o Third character describes hybridization or geometry: 1 - linear;
      2 - trigonal; R - resonant; 3 = tetrahedral; 4 - square planar;
      5 - trigonal bipyramidal; 6 - octahedral
    o Fourth and fifth characters are used as indicators of alternate
      parameters: formal oxidation state, bridging hydrogens and so on.

Examples of UFF atom types:

    C_3, C_2, C_R, N_3, N_R, O_3, O_2, and so on

=head2 METHODS

=over 4

=item B<new>

    $NewUFFAtomTypes = new AtomTypes::UFFAtomTypes(%NamesAndValues);

Using specified I<UFFAtomTypes> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<UFFAtomTypes>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'UFF'
    IgnoreHydrogens = 0

Examples:

    $UFFAtomTypes = new AtomTypes::UFFAtomTypes(
                              'Molecule' => $Molecule,
                              'IgnoreHydrogens' => 0);

=item B<AssignAtomTypes>

    $UFFAtomTypes->AssignAtomTypes();

Assigns UFF atom types to all the atoms in a molecule and returns
I<UFFAtomTypes>.

=item B<GetAllPossibleUFFAtomTypes>

    $AllAtomTypesDataRef = $UFFAtomTypes->
                           GetAllPossibleUFFAtomTypes();
    $AllAtomTypesDataRef = AtomTypes::UFFAtomTypes::
                           GetAllPossibleUFFAtomTypes();

Returns all possible UFF atom types corresponding to hydrogen and non-hydrogen
atoms as an array reference.

=item B<GetAllPossibleUFFNonHydrogenAtomTypes>

    $AtomTypesDataRef = $UFFAtomTypes->
                        GetAllPossibleUFFNonHydrogenAtomTypes();
    $AtomTypesDataRef = AtomTypes::UFFAtomTypes::
                        GetAllPossibleUFFNonHydrogenAtomTypes();

Returns all possible UFF atom types corresponding to non-hydrogen atoms as
an array reference.

=item B<GetUFFAtomTypesData>

    $AtomTypesDataMapRef = $UFFAtomTypes->GetUFFAtomTypesData();
    $AtomTypesDataMapRef = AtomTypes::UFFAtomTypes::GetUFFAtomTypesData();

Returns UFF atom types and associated data loaded from UFF data file as
a reference to hash with the following hash data format:

    @{$UFFAtomTypesDataMap{AtomTypes}} - Array of all possible atom
                              types for all atoms
    @{$UFFAtomTypesDataMap{NonHydrogenAtomTypes}} - Array of all
                              possible atom types for non-hydrogen atoms
    @{$UFFAtomTypesDataMap->{ColLabels}} - Array of column labels
    %{$UFFAtomTypesDataMap->{DataCol<Num>}} - Hash keys pair:
                                              DataCol<Num>, AtomType


=item B<StringifyUFFAtomTypes>

    $String = $UFFAtomTypes->StringifyUFFAtomTypes();

Returns a string containing information about I<UFFAtomTypes> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AtomTypes.pm, AtomicInvariantsAtomTypes.pm, DREIDINGAtomTypes.pm,
EStateAtomTypes.pm, FunctionalClassAtomTypes.pm, MMFF94AtomTypes.pm,
SLogPAtomTypes.pm, SYBYLAtomTypes.pm, TPSAAtomTypes.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
