package AtomTypes::SYBYLAtomTypes;
#
# File: SYBYLAtomTypes.pm
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
@EXPORT = qw(GetSYBYLAtomTypesData GetAllPossibleSYBYLAtomTypes GetAllPossibleSYBYLNonHydrogenAtomTypes);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, %SYBYLAtomTypesDataMap);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifySYBYLAtomTypes';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeSYBYLAtomTypes();

  $This->_InitializeSYBYLAtomTypesProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Initialize the data hash. It'll be loaded on demand later...
  %SYBYLAtomTypesDataMap = ();
}


# Initialize object data...
#
sub _InitializeSYBYLAtomTypes {
  my($This) = @_;

  # Type of AtomTypes...
  $This->{Type} = 'SYBYL';

  # By default, SYBYL atom types are also assigned to hydrogens...
  $This->{IgnoreHydrogens} = 0;

  return $This;
}

# Initialize object properties...
#
sub _InitializeSYBYLAtomTypesProperties {
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

# Get SYBYL atom types and associated data loaded from SYBYL data file as
# a reference to hash with the following hash data format:
#
# @{$SYBYLAtomTypesDataMap{AtomTypes}} - Array of all possible atom types for all atoms
# @{$SYBYLAtomTypesDataMap{NonHydrogenAtomTypes}} - Array of all possible atom types for non-hydrogen atoms
# @{$SYBYLAtomTypesDataMap->{ColLabels}} - Array of column labels
# %{$SYBYLAtomTypesDataMap->{DataCol<Num>}} - Hash keys pair: <DataCol<Num>, AtomType>
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetSYBYLAtomTypesData {

  # Make sure data is loaded...
  _CheckAndLoadSYBYLAtomTypesData();

  return \%SYBYLAtomTypesDataMap;
}

# Get all possible SYBYL atom types corresponding to hydrogen and non-hydrogen
# atoms as an array reference...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAllPossibleSYBYLAtomTypes {
  return _GetAllPossibleSYBYLAtomTypes();
}

# Get all possible SYBYL atom types corresponding to non-hydrogen atoms
# as an array reference...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAllPossibleSYBYLNonHydrogenAtomTypes {
  my($NonHydrogensOnly);

  $NonHydrogensOnly = 1;
  return _GetAllPossibleSYBYLAtomTypes($NonHydrogensOnly);
}

# Get all possible SYBYL atom types as an array reference...
#
sub _GetAllPossibleSYBYLAtomTypes {
  my($NonHydrogensOnly) = @_;
  my($SYBYLAtomTypesDataRef);

  $NonHydrogensOnly = defined $NonHydrogensOnly ? $NonHydrogensOnly : 0;

  $SYBYLAtomTypesDataRef = GetSYBYLAtomTypesData();

  return $NonHydrogensOnly ? \@{$SYBYLAtomTypesDataRef->{NonHydrogenAtomTypes}}: \@{$SYBYLAtomTypesDataRef->{AtomTypes}};
}
# Assign Tripos SYBYL [ Ref 79-80 ] atom types to all atoms...
#
# Notes:
#   . 8 SYBYL listed atom types - O.spc, O.t3p, H.spc, H.t3p, Du, Du.C, HEV, LP -
#     are not assigned to any atom
#   . N.pl3 atom type is assigned to Nitrogens in a guadinium group attached
#     to Carbon C.cat
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

# Get SYBYL atom type for atom...
#
sub _GetAtomType {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'Any';

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

# Get SYBYL atom type for Carbon atom...
#
sub _GetAtomTypeForCarbon {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'Any';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {
    if ($Atom->IsAromatic()) {
      $AtomType = 'C.ar';
      last ATOMTYPE;
    }

    if ($Atom->IsGuadiniumCarbon()) {
      $AtomType = 'C.cat';
      last ATOMTYPE;
    }

    # Only single bonds...
    if ($NumOfPiBonds == 0) {
      $AtomType = 'C.3';
      last ATOMTYPE;
    }

    # One double bond...
    if ($NumOfPiBonds == 1) {
      $AtomType = 'C.2';
      last ATOMTYPE;
    }

    # One triple bond or two double bonds...
    if ($NumOfPiBonds == 2) {
      $AtomType = 'C.1';
      last ATOMTYPE;
    }

    $AtomType = 'Any';
  }

  return $AtomType;
}

# Get SYBYL atom type for Nitrogen atom...
#
sub _GetAtomTypeForNitrogen {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'Het';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {
    if ($Atom->IsAromatic()) {
      $AtomType = 'N.ar';
      last ATOMTYPE;
    }

    if ($Atom->IsGuadiniumNitrogen()) {
      $AtomType = 'N.pl3';
      last ATOMTYPE;
    }

    if ($Atom->IsAmideNitrogen()) {
      $AtomType = 'N.am';
      last ATOMTYPE;
    }

    # Only single bonds...
    if ($NumOfSigmaBonds == 4 && $NumOfPiBonds == 0) {
      $AtomType = 'N.4';
      last ATOMTYPE;
    }

    # Only single bonds...
    if ($NumOfSigmaBonds == 3 && $NumOfPiBonds == 0) {
      $AtomType = $This->_IsN3NitrogenPlanar($Atom) ? 'N.pl3' : 'N.3';
      last ATOMTYPE;
    }

    # One double bond...
    if ($NumOfPiBonds == 1) {
      $AtomType = 'N.2';
      last ATOMTYPE;
    }

    # One triple bond or two double bonds...
    if ($NumOfPiBonds == 2) {
      $AtomType = 'N.1';
      last ATOMTYPE;
    }

    $AtomType = 'Het';
  }

  return $AtomType;
}

# Get SYBYL atom type for Oxygen atom...
#
sub _GetAtomTypeForOxygen {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'Het';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {
    if ($Atom->IsCarboxylateOxygen() || $Atom->IsCarboxylOxygen() || $Atom->IsPhosphateOxygen()) {
      $AtomType = 'O.co2';
      last ATOMTYPE;
    }

    # Only single bonds...
    if ($NumOfPiBonds == 0) {
      $AtomType = 'O.3';
      last ATOMTYPE;
    }

    # One double bond...
    if ($NumOfPiBonds == 1) {
      $AtomType = 'O.2';
      last ATOMTYPE;
    }

    $AtomType = 'Het';
  }

  return $AtomType;
}

# Get SYBYL atom type for Phosphorus atom...
#
sub _GetAtomTypeForPhosphorus {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'Het';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {
    # -P(-)-, =P(-)(-)-
    if (($NumOfSigmaBonds == 3 && $NumOfPiBonds == 0) || ($NumOfSigmaBonds == 4 && $NumOfPiBonds == 1)) {
      $AtomType = 'P.3';
      last ATOMTYPE;
    }

    $AtomType = 'Het';
  }

  return $AtomType;
}

# Get SYBYL atom type for Sulfur atom...
#
sub _GetAtomTypeForSulfur {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'Het';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {
    if ($This->_IsSulfoneSulfur($Atom)) {
      $AtomType = 'S.O2';
      last ATOMTYPE;
    }

    if ($This->_IsSulfoxideSulfur($Atom)) {
      $AtomType = 'S.O';
      last ATOMTYPE;
    }

    # -S-
    if ($NumOfSigmaBonds == 2 && $NumOfPiBonds == 0) {
      $AtomType = 'S.3';
      last ATOMTYPE;
    }

    # S=
    if ($NumOfSigmaBonds == 1 && $NumOfPiBonds == 1) {
      $AtomType = 'S.2';
      last ATOMTYPE;
    }

    $AtomType = 'Het';
  }

  return $AtomType;
}

# Get SYBYL atom type for Hydrogen atom...
#
sub _GetAtomTypeForHydrogen {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'H';

  return $AtomType;
}

# Get SYBYL atom type for atoms other than Carbon, Nitrogen, Oxygen, Phosporus
# and Sulfur...
#
sub _GetAtomTypeForOtherAtoms {
  my($This, $Atom) = @_;
  my($AtomType, $AtomicNumber, $AtomSymbol);

  $AtomType = 'Any';

  $AtomicNumber = $Atom->GetAtomicNumber();
  $AtomSymbol = $Atom->GetAtomSymbol();

  ATOMICNUMBER: {
    if ($AtomicNumber =~ /^(9|17|35|53)$/i) {
      # F, Cl, Br, I
      $AtomType = $AtomSymbol;
      last ATOMICNUMBER;
    }

    if ($AtomicNumber =~ /^(3|11|12|13|14)$/i) {
      # Li, Na, Mg, Al, Si
      $AtomType = $AtomSymbol;
      last ATOMICNUMBER;
    }

    if ($AtomicNumber =~ /^(19|20|25|26|29|30|34)$/i) {
      # K, Ca, Mn, Fe, Cu, Zn, Se
      $AtomType = $AtomSymbol;
      last ATOMICNUMBER;
    }

    if ($AtomicNumber =~ /^24$/i) {
      $AtomType = $This->_GetAtomTypeForChromium($Atom);
      last ATOMICNUMBER;
    }

    if ($AtomicNumber =~ /^27$/i) {
      $AtomType = $This->_GetAtomTypeForCobalt($Atom);
      last ATOMICNUMBER;
    }

    if ($AtomicNumber =~ /^(42|50)$/i) {
      # Mo, Sn
      $AtomType = $AtomSymbol;
      last ATOMICNUMBER;
    }

    $AtomType = 'Any';
  }

  return $AtomType;
}

# Get SYBYL atom type for Chromium atom...
#
sub _GetAtomTypeForChromium {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfNeighbors);

  $AtomType = 'Any';
  $NumOfNeighbors = $Atom->GetNumOfNeighbors();

  NUMOFNEIGHBORS: {
    if ($NumOfNeighbors == 4) {
      $AtomType = 'Cr.th';
      last NUMOFNEIGHBORS;
    }

    if ($NumOfNeighbors == 6) {
      $AtomType = 'Cr.oh';
      last NUMOFNEIGHBORS;
    }

    $AtomType = 'Cr.oh';
    carp "Warning: ${ClassName}->_GetAtomTypeForChromium: SYBYL atom types for Cromimum, Co.th or Cr.oh, corresponding to tetrahedral or octahedral geometry cann't be assigned; Number of neighbors, $NumOfNeighbors, is different from 4 or 6. Default SYBYL atom type, Cr.oh, has been assigned...";
  }

  return $AtomType;
}

# Get SYBYL atom type for Cobalt atom...
#
sub _GetAtomTypeForCobalt {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfNeighbors);

  $AtomType = 'Any';

  $NumOfNeighbors = $Atom->GetNumOfNeighbors();

  if ($NumOfNeighbors == 6) {
    $AtomType = 'Co.oh';
  }
  else {
    $AtomType = 'Co.oh';
    carp "Warning: ${ClassName}->_GetAtomTypeForCobalt: SYBYL atom type for Cobalt, Co.oh, corresponding to octahedral geometry cann't be assigned; Number of neighbors, $NumOfNeighbors, is different from 6. Default SYBYL atom type, Co.oh, has been assigned...";
  }

  return $AtomType;
}

# Is it N.3 Nitrogen a planar Nitrogen?
#
# A N.3 Nitrogen is a planar Nitrogen when it meets any of the following constraints:
#
#   . Nitrogen atom is in a ring
#   . Attached to only one heavy atom which is an aromatic atom or in a ring
#   . Attached to two or more heavy atom which are aromatic atoms or in a ring
#
sub _IsN3NitrogenPlanar {
  my($This, $Atom) = @_;

  # Is it a ring Nitrogen atom?
  if ($Atom->IsInRing()) {
    return 1;
  }

  # Count number of ring and aromatic heavy atoms attached to Nitrogen...
  my($AtomNeighbor, $NumOfAromaticAtomNeighbors, $NumOfRingAotmNeighbors, $NumOfNeighbors, @AtomNeighbors);

  @AtomNeighbors = $Atom->GetHeavyAtomNeighbors();
  $NumOfNeighbors = scalar @AtomNeighbors;

  $NumOfAromaticAtomNeighbors = 0;  $NumOfRingAotmNeighbors = 0;

  for $AtomNeighbor (@AtomNeighbors) {
    if ($AtomNeighbor->IsAromatic()) {
      $NumOfAromaticAtomNeighbors++;
    }
    if ($AtomNeighbor->IsInRing()) {
      $NumOfRingAotmNeighbors++;
    }
  }

  # Is attached to only one heavy atom which is in a ring or aromatic?
  if ($NumOfNeighbors == 1) {
    if ($NumOfAromaticAtomNeighbors || $NumOfRingAotmNeighbors) {
      return 1;
    }
  }

  # Is attached to more than heavy atoms which are in a ring or aromatic?
  if ($NumOfAromaticAtomNeighbors >= 2 || $NumOfRingAotmNeighbors >= 2) {
    return 1;
  }

  return 0;
}

# Is it a Sulfur atom in Sulfoxide group?
#
# SYBYL Sulfoxide group definition: A-S(=O)-A
#
#   where:
#      . A = Any atom
#
sub _IsSulfoxideSulfur {
  my($This, $Atom) = @_;

  # Is it Sulfur?
  if (!$Atom->IsSulfur()) {
    return 0;
  }
  # Match atom neighborhood...
  my($CentralAtomSpec, @NbrAtomSpecsRef, @NbrBondSpecsRef);

  $CentralAtomSpec = 'S.X3.BO4';
  @NbrAtomSpecsRef = ('O', '*', '*');
  @NbrBondSpecsRef = ('=', '-', '-');

  if ($Atom->DoesAtomNeighborhoodMatch($CentralAtomSpec, \@NbrAtomSpecsRef, \@NbrBondSpecsRef)) {
    return 1;
  }

  return 0;
}

# Is it a Sulfur atom in Sulfone group?
#
# Sulfoxide group definition: A-(O=)S(=O)-A
#
#   where:
#      . A = Any atom
#
sub _IsSulfoneSulfur {
  my($This, $Atom) = @_;

  # Is it Sulfur?
  if (!$Atom->IsSulfur()) {
    return 0;
  }

  # Match atom neighborhood...
  my($CentralAtomSpec, @NbrAtomSpecsRef, @NbrBondSpecsRef);

  $CentralAtomSpec = 'S.X4.BO6';
  @NbrAtomSpecsRef = ('O', 'O', '*', '*');
  @NbrBondSpecsRef = ('=', '=', '-', '-');

  if ($Atom->DoesAtomNeighborhoodMatch($CentralAtomSpec, \@NbrAtomSpecsRef, \@NbrBondSpecsRef)) {
    return 1;
  }

  return 0;
}

# Return a string containg data for SYBYLAtomTypes object...
#
sub StringifySYBYLAtomTypes {
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

# Is it a SYBYLAtomTypes object?
sub _IsSYBYLAtomTypes {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Check and load SYBYL atom types data...
#
sub _CheckAndLoadSYBYLAtomTypesData {

  # Is it already loaded?
  if (exists $SYBYLAtomTypesDataMap{AtomTypes}) {
    return;
  }

  _LoadSYBYLAtomTypesData();
}

# Load SYBYL atom types data from the file assuming first column to be atom type symbol..
#
# Format:
#
# "AtomType","Description"
# "C.3","sp3 carbon"
# "C.2","sp2 carbon"
#
sub _LoadSYBYLAtomTypesData {
  my($AtomTypesDataFile, $MayaChemToolsLibDir);

  $MayaChemToolsLibDir = FileUtil::GetMayaChemToolsLibDirName();

  $AtomTypesDataFile =  "$MayaChemToolsLibDir" . "/data/SYBYLAtomTypes.csv";
  if (! -e "$AtomTypesDataFile") {
    croak "Error: MayaChemTools package file, $AtomTypesDataFile, is missing: Possible installation problems...";
  }

  %SYBYLAtomTypesDataMap = ();
  AtomTypes::AtomTypes::LoadAtomTypesData($AtomTypesDataFile, \%SYBYLAtomTypesDataMap);
}


1;

__END__

=head1 NAME

SYBYLAtomTypes

=head1 SYNOPSIS

use AtomTypes::SYBYLAtomTypes;

use AtomTypes::SYBYLAtomTypes qw(:all);

=head1 DESCRIPTION

B<SYBYLAtomTypes> class provides the following methods:

new, AssignAtomTypes, GetAllPossibleSYBYLAtomTypes,
GetAllPossibleSYBYLNonHydrogenAtomTypes, GetSYBYLAtomTypesData,
StringifySYBYLAtomTypes

The following functions are available:

GetAllPossibleSYBYLAtomTypes,
GetAllPossibleSYBYLNonHydrogenAtomTypes, GetSYBYLAtomTypesData

B<SYBYLAtomTypes> is derived from B<AtomTypes> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<SYBYLAtomTypes>, B<AtomTypes> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The data file SYBYLAomTypes.csv distributed with MayaChemTools release contains
all possible Triops SYBYL [ Ref 79-80 ] atom types.

Examples of SYBYL atom types:

    C.3,C.2, C.ar, N.3, N.2, N.ar and so on

=head2 METHODS

=over 4

=item B<new>

    $NewSYBYLAtomTypes = new AtomTypes::SYBYLAtomTypes(%NamesAndValues);

Using specified I<SYBYLAtomTypes> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<SYBYLAtomTypes>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'SYBYL'
    IgnoreHydrogens = 0

Examples:

    $SYBYLAtomTypes = new AtomTypes::SYBYLAtomTypes(
                              'Molecule' => $Molecule,
                              'IgnoreHydrogens' => 0);

=item B<AssignAtomTypes>

    $SYBYLAtomTypes->AssignAtomTypes();

Assigns SYBYL atom types to all the atoms in a molecule and returns
I<SYBYLAtomTypes>.

=item B<GetAllPossibleSYBYLAtomTypes>

    $AllAtomTypesDataRef = $SYBYLAtomTypes->
                           GetAllPossibleSYBYLAtomTypes();
    $AllAtomTypesDataRef = AtomTypes::SYBYLAtomTypes::
                           GetAllPossibleSYBYLAtomTypes();

Returns all possible SYBYL atom types corresponding to hydrogen and non-hydrogen
atoms as an array reference.

=item B<GetAllPossibleSYBYLNonHydrogenAtomTypes>

    $AtomTypesDataRef = $SYBYLAtomTypes->
                        GetAllPossibleSYBYLNonHydrogenAtomTypes();
    $AtomTypesDataRef = AtomTypes::SYBYLAtomTypes::
                        GetAllPossibleSYBYLNonHydrogenAtomTypes();

Returns all possible SYBYL atom types corresponding to non-hydrogen atoms as
an array reference.

=item B<GetSYBYLAtomTypesData>

    $AtomTypesDataMapRef = $SYBYLAtomTypes->GetSYBYLAtomTypesData();
    $AtomTypesDataMapRef = AtomTypes::SYBYLAtomTypes::GetSYBYLAtomTypesData();

Returns SYBYL atom types and associated data loaded from SYBYL data file as
a reference to hash with the following hash data format:

    @{$SYBYLAtomTypesDataMap{AtomTypes}} - Array of all possible atom
                              types for all atoms
    @{$SYBYLAtomTypesDataMap{NonHydrogenAtomTypes}} - Array of all
                              possible atom types for non-hydrogen atoms
    @{$SYBYLAtomTypesDataMap->{ColLabels}} - Array of column labels
    %{$SYBYLAtomTypesDataMap->{DataCol<Num>}} - Hash keys pair:
                                                DataCol<Num>, AtomType


=item B<StringifySYBYLAtomTypes>

    $String = $SYBYLAtomTypes->StringifySYBYLAtomTypes();

Returns a string containing information about I<SYBYLAtomTypes> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AtomTypes.pm, AtomicInvariantsAtomTypes.pm, DREIDINGAtomTypes.pm,
EStateAtomTypes.pm, FunctionalClassAtomTypes.pm, MMFF94AtomTypes.pm,
SLogPAtomTypes.pm, TPSAAtomTypes.pm, UFFAtomTypes.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
