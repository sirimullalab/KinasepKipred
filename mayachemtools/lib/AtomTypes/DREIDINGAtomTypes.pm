package AtomTypes::DREIDINGAtomTypes;
#
# File: DREIDINGAtomTypes.pm
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
@EXPORT = qw(GetDREIDINGAtomTypesData GetAllPossibleDREIDINGAtomTypes GetAllPossibleDREIDINGNonHydrogenAtomTypes);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, %DREIDINGAtomTypesDataMap);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyDREIDINGAtomTypes';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeDREIDINGAtomTypes();

  $This->_InitializeDREIDINGAtomTypesProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Initialize the data hash. It'll be loaded on demand later...
  %DREIDINGAtomTypesDataMap = ();
}


# Initialize object data...
#
sub _InitializeDREIDINGAtomTypes {
  my($This) = @_;

  # Type of AtomTypes...
  $This->{Type} = 'DREIDING';

  # By default, DREIDING atom types are also assigned to hydrogens...
  $This->{IgnoreHydrogens} = 0;

  return $This;
}

# Initialize object properties...
#
sub _InitializeDREIDINGAtomTypesProperties {
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

# Get DREIDING atom types and associated data loaded from DREIDING data file as
# a reference to hash with the following hash data format:
#
# @{$DREIDINGAtomTypesDataMap{AtomTypes}} - Array of all possible atom types for all atoms
# @{$DREIDINGAtomTypesDataMap{NonHydrogenAtomTypes}} - Array of all possible atom types for non-hydrogen atoms
# @{$DREIDINGAtomTypesDataMap->{ColLabels}} - Array of column labels
# %{$DREIDINGAtomTypesDataMap->{DataCol<Num>}} - Hash keys pair: <DataCol<Num>, AtomType>
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetDREIDINGAtomTypesData {

  # Make sure data is loaded...
  _CheckAndLoadDREIDINGAtomTypesData();

  return \%DREIDINGAtomTypesDataMap;
}

# Get all possible DREIDING atom types corresponding to hydrogen and non-hydrogen
# atoms as an array reference...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAllPossibleDREIDINGAtomTypes {
  return _GetAllPossibleDREIDINGAtomTypes();
}

# Get all possible DREIDING atom types corresponding to non-hydrogen atoms
# as an array reference...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAllPossibleDREIDINGNonHydrogenAtomTypes {
  my($NonHydrogensOnly);

  $NonHydrogensOnly = 1;
  return _GetAllPossibleDREIDINGAtomTypes($NonHydrogensOnly);
}

# Get all possible DREIDING atom types as an array reference...
#
sub _GetAllPossibleDREIDINGAtomTypes {
  my($NonHydrogensOnly) = @_;
  my($DREIDINGAtomTypesDataRef);

  $NonHydrogensOnly = defined $NonHydrogensOnly ? $NonHydrogensOnly : 0;

  $DREIDINGAtomTypesDataRef = GetDREIDINGAtomTypesData();

  return $NonHydrogensOnly ? \@{$DREIDINGAtomTypesDataRef->{NonHydrogenAtomTypes}}: \@{$DREIDINGAtomTypesDataRef->{AtomTypes}};
}

# Assign DREIDING [ Ref 88 ] atom types to all atoms...
#
# Notes:
#     o 37 DREIDING atom types are listed
#     o AtomTypes::DREIDINGAtomTypes.pm module is used to assign DREIDING atom types
#     o Units:
#         o ValenceBondRadius and NonBondRadius: Angstroms
#         o ValenceAngle: Degrees
#     o Five-character mnemonic label for DREIDING atom types
#         o First two characters correspond to chemical symbol with an underscore as second
#           character for elements with one character symbol
#         o Third character describes hybridization: 1 - linear (sp); 2 - trigonal (sp2);
#           3 = tetrahedral (sp3); R - sp2 involved in resonance situation
#         o Fourth character used to indicate number of implicit hydrogens
#         o Fourth and fifth chracters are used as indicators of alternate parameters: formal oxidation
#           state, bridging hydrogens and so on. The _HB type denotes a hydrogen atom capable
#           of forming hdyrogen bonds attached to (N, O, F). The H_b is the bridging hydrogen
#           of diborane.
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

# Get DREIDING atom type for atom...
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

# Get DREIDING atom type for Carbon atom...
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
    carp "Warning: ${ClassName}->_GetAtomTypeForCarbon: DREIDING atom types for Carbon cann't be assigned...";
  }

  return $AtomType;
}

# Get DREIDING atom type for Nitrogen atom...
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
    carp "Warning: ${ClassName}->_GetAtomTypeForNitrogen: DREIDING atom types for Nitrogen cann't be assigned...";
  }

  return $AtomType;
}

# Get DREIDING atom type for Oxygen atom...
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
    carp "Warning: ${ClassName}->_GetAtomTypeForOxygen: DREIDING atom types for Oxygen cann't be assigned...";
  }

  return $AtomType;
}

# Get DREIDING atom type for Phosphorus atom...
#
sub _GetAtomTypeForPhosphorus {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'P_3';

  return $AtomType;
}

# Get DREIDING atom type for Sulfur atom...
#
sub _GetAtomTypeForSulfur {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'S_3';

  return $AtomType;
}

# Get DREIDING atom type for Hydrogen atom...
#
sub _GetAtomTypeForHydrogen {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfNeighbors, $NeighborAtom, @NonHydrogenAtomNeighbors);

  @NonHydrogenAtomNeighbors = $Atom->GetNonHydrogenAtomNeighbors();

  $NumOfNeighbors = scalar @NonHydrogenAtomNeighbors;
  $NeighborAtom = $NonHydrogenAtomNeighbors[0];

  ATOMTYPE: {
    if ($NumOfNeighbors > 1) {
      # Bridging hydrogen as in B2H6
      $AtomType = 'H___b';
      last ATOMTYPE;
    }

    if ($NeighborAtom->GetAtomicNumber() =~ /^(7|8|9)$/) {
      # Involved in hydrogen bonding due to its attachment to N, O, or F
      $AtomType = 'H__HB';
      last ATOMTYPE;
    }
    $AtomType = 'H_';
  }

  return $AtomType;
}

# Get DREIDING atom type for atoms other than Carbon, Nitrogen, Oxygen, Phosporus,
# Sulfur and Hydrogen...
#
sub _GetAtomTypeForOtherAtoms {
  my($This, $Atom) = @_;
  my($AtomType, $AtomicNumber, $AtomSymbol);

  $AtomType = 'None';

  $AtomicNumber = $Atom->GetAtomicNumber();
  $AtomSymbol = $Atom->GetAtomSymbol();

  ATOMICNUMBER: {
    if ($AtomicNumber =~ /^(9|17|35|53)$/i) {
      # F, Cl, Br, I
      $AtomType = length($AtomSymbol) == 1 ? "${AtomSymbol}_" : $AtomSymbol;
      last ATOMICNUMBER;
    }

    if ($AtomicNumber =~ /^5$/i) {
      # B: B_2 and B_3
      $AtomType = (($Atom->GetNumOfNonHydrogenAtomNeighbors() + $Atom->GetAtomicInvariantValue('H')) == 4) ? "B_3" : "B_2";
      last ATOMICNUMBER;
    }

    if ($AtomicNumber =~ /^(13|14|31|32|33|34|49|50|51|52)$/i) {
      # Al, Si, Ga, Ge, As, Se, In, Sn, Sb, Te
      $AtomType = "${AtomSymbol}3";
      last ATOMICNUMBER;
    }

    if ($AtomicNumber =~ /^(11|20|26|30)$/i) {
      # Na, Ca, Fe, Zn
      $AtomType = $AtomSymbol;
      last ATOMICNUMBER;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForOtherAtoms: DREIDING atom types for atom, $AtomSymbol, with atomic number, $AtomicNumber, cann't be assigned...";
  }

  return $AtomType;
}

# Return a string containg data for DREIDINGAtomTypes object...
#
sub StringifyDREIDINGAtomTypes {
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

# Is it a DREIDINGAtomTypes object?
sub _IsDREIDINGAtomTypes {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Check and load DREIDING atom types data...
#
sub _CheckAndLoadDREIDINGAtomTypesData {

  # Is it already loaded?
  if (exists $DREIDINGAtomTypesDataMap{AtomTypes}) {
    return;
  }

  _LoadDREIDINGAtomTypesData();
}

# Load DREIDING atom types data from the file assuming first column to be atom type symbol..
#
# Format:
#
# "AtomType","ValenceBondRadius","ValenceAngle"
# "H_","0.330","180.0"
# "C_3","0.770","109.471"
# "C_R","0.700","120.0"
# "C_2","0.670","120.0"
# "C_1","0.602","180.0"
# "N_3","0.702","106.7"
#
sub _LoadDREIDINGAtomTypesData {
  my($AtomTypesDataFile, $MayaChemToolsLibDir);

  $MayaChemToolsLibDir = FileUtil::GetMayaChemToolsLibDirName();

  $AtomTypesDataFile =  "$MayaChemToolsLibDir" . "/data/DREIDINGAtomTypes.csv";
  if (! -e "$AtomTypesDataFile") {
    croak "Error: MayaChemTools package file, $AtomTypesDataFile, is missing: Possible installation problems...";
  }

  %DREIDINGAtomTypesDataMap = ();
  AtomTypes::AtomTypes::LoadAtomTypesData($AtomTypesDataFile, \%DREIDINGAtomTypesDataMap);
}

1;

__END__

=head1 NAME

DREIDINGAtomTypes

=head1 SYNOPSIS

use AtomTypes::DREIDINGAtomTypes;

use AtomTypes::DREIDINGAtomTypes qw(:all);

=head1 DESCRIPTION

B<DREIDINGAtomTypes> class provides the following methods:

new, AssignAtomTypes, GetAllPossibleDREIDINGAtomTypes,
GetAllPossibleDREIDINGNonHydrogenAtomTypes, GetDREIDINGAtomTypesData,
StringifyDREIDINGAtomTypes

The following functions are available:

GetAllPossibleDREIDINGAtomTypes,
GetAllPossibleDREIDINGNonHydrogenAtomTypes, GetDREIDINGAtomTypesData

B<DREIDINGAtomTypes> is derived from B<AtomTypes> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<DREIDINGAtomTypes>, B<AtomTypes> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The data file DREIDINGAtomTypes.csv distributed with MayaChemTools release contains
all possible DREIDING [ Ref 88 ] atom types.

Format of a Five-character mnemonic label used for DREIDING atom types:

    o First two characters correspond to chemical symbol with an
      underscore as second character for elements with one character symbol
    o Third character describes hybridization: 1 - linear (sp);
      2 - trigonal (sp2); 3 = tetrahedral (sp3); R - sp2 involved in
      resonance situation
    o Fourth character used to indicate number of implicit hydrogens
    o Fourth and fifth characters are used as indicators of alternate
      parameters: formal oxidation state, bridging hydrogens and so on.
      The _HB type denotes a hydrogen atom capable of forming hydrogen
      bonds attached to (N, O, F). The H_b is the bridging hydrogen
      of diborane.

Examples of DREIDING atom types:

    H_, C_3, C_R, C_2, N_3, N_R, O_3, O_R and so on

=head2 METHODS

=over 4

=item B<new>

    $NewDREIDINGAtomTypes = new AtomTypes::DREIDINGAtomTypes(%NamesAndValues);

Using specified I<DREIDINGAtomTypes> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<DREIDINGAtomTypes>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'DREIDING'
    IgnoreHydrogens = 0

Examples:

    $DREIDINGAtomTypes = new AtomTypes::DREIDINGAtomTypes(
                              'Molecule' => $Molecule,
                              'IgnoreHydrogens' => 0);

=item B<AssignAtomTypes>

    $DREIDINGAtomTypes->AssignAtomTypes();

Assigns DREIDING atom types to all the atoms in a molecule and returns
I<DREIDINGAtomTypes>.

=item B<GetAllPossibleDREIDINGAtomTypes>

    $AllAtomTypesDataRef = $DREIDINGAtomTypes->
                           GetAllPossibleDREIDINGAtomTypes();
    $AllAtomTypesDataRef = AtomTypes::DREIDINGAtomTypes::
                           GetAllPossibleDREIDINGAtomTypes();

Returns all possible DREIDING atom types corresponding to hydrogen and non-hydrogen
atoms as an array reference.

=item B<GetAllPossibleDREIDINGNonHydrogenAtomTypes>

    $AtomTypesDataRef = $DREIDINGAtomTypes->
                        GetAllPossibleDREIDINGNonHydrogenAtomTypes();
    $AtomTypesDataRef = AtomTypes::DREIDINGAtomTypes::
                        GetAllPossibleDREIDINGNonHydrogenAtomTypes();

Returns all possible DREIDING atom types corresponding to non-hydrogen atoms as
an array reference.

=item B<GetDREIDINGAtomTypesData>

    $AtomTypesDataMapRef = $DREIDINGAtomTypes->GetDREIDINGAtomTypesData();
    $AtomTypesDataMapRef = AtomTypes::DREIDINGAtomTypes::
                           GetDREIDINGAtomTypesData();

Returns DREIDING atom types and associated data loaded from DREIDING data file as
a reference to hash with the following hash data format:

    @{$DREIDINGAtomTypesDataMap{AtomTypes}} - Array of all possible atom
                              types for all atoms
    @{$DREIDINGAtomTypesDataMap{NonHydrogenAtomTypes}} - Array of all
                              possible atom types for non-hydrogen atoms
    @{$DREIDINGAtomTypesDataMap->{ColLabels}} - Array of column labels
    %{$DREIDINGAtomTypesDataMap->{DataCol<Num>}} - Hash keys pair:
                                                   DataCol<Num>, AtomType

=item B<StringifyDREIDINGAtomTypes>

    $String = $DREIDINGAtomTypes->StringifyDREIDINGAtomTypes();

Returns a string containing information about I<DREIDINGAtomTypes> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AtomTypes.pm, AtomicInvariantsAtomTypes.pm, EStateAtomTypes.pm,
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
