package AtomTypes::MMFF94AtomTypes;
#
# File: MMFF94AtomTypes.pm
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
use Text::ParseWords;
use AtomTypes::AtomTypes;
use Molecule;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(AtomTypes::AtomTypes Exporter);
@EXPORT = qw(GetMMFF94AtomTypesData GetAllPossibleMMFF94AtomTypes GetAllPossibleMMFF94NonHydrogenAtomTypes);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, %MMFF94AtomTypesDataMap);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyMMFF94AtomTypes';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeMMFF94AtomTypes();

  $This->_InitializeMMFF94AtomTypesProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Initialize the data hash. It'll be loaded on demand later...
  %MMFF94AtomTypesDataMap = ();
}


# Initialize object data...
#
sub _InitializeMMFF94AtomTypes {
  my($This) = @_;

  # Type of AtomTypes...
  $This->{Type} = 'MMFF94';

  # By default, MMFF94 atom types are also assigned to hydrogens...
  $This->{IgnoreHydrogens} = 0;

  return $This;
}

# Initialize object properties...
#
sub _InitializeMMFF94AtomTypesProperties {
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

# Get MMFF94 atom types and associated data loaded from MMFF94 data file as
# a reference to hash with the following hash data format:
#
# @{$MMFF94AtomTypesDataMap{AtomTypes}} - Array of all possible atom types for all atoms
# @{$MMFF94AtomTypesDataMap{NonHydrogenAtomTypes}} - Array of all possible atom types for non-hydrogen atoms
# @{$MMFF94AtomTypesDataMap->{ColLabels}} - Array of column labels
# %{$MMFF94AtomTypesDataMap->{DataCol<Num>}} - Hash keys pair: <DataCol<Num>, AtomType>
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetMMFF94AtomTypesData {

  # Make sure data is loaded...
  _CheckAndLoadMMFF94AtomTypesData();

  return \%MMFF94AtomTypesDataMap;
}

# Get all possible MMFF94 atom types corresponding to hydrogen and non-hydrogen
# atoms as an array reference...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAllPossibleMMFF94AtomTypes {
  return _GetAllPossibleMMFF94AtomTypes();
}

# Get all possible MMFF94 atom types corresponding to non-hydrogen atoms
# as an array reference...
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAllPossibleMMFF94NonHydrogenAtomTypes {
  my($NonHydrogensOnly);

  $NonHydrogensOnly = 1;
  return _GetAllPossibleMMFF94AtomTypes($NonHydrogensOnly);
}

# Get all possible MMFF94 atom types as an array reference...
#
sub _GetAllPossibleMMFF94AtomTypes {
  my($NonHydrogensOnly) = @_;
  my($MMFF94AtomTypesDataRef);

  $NonHydrogensOnly = defined $NonHydrogensOnly ? $NonHydrogensOnly : 0;

  $MMFF94AtomTypesDataRef = GetMMFF94AtomTypesData();

  return $NonHydrogensOnly ? \@{$MMFF94AtomTypesDataRef->{NonHydrogenAtomTypes}}: \@{$MMFF94AtomTypesDataRef->{AtomTypes}};
}

# Assign MMFF94 [ Ref 83-87 ] atom types to all atoms...
#
# Notes:
#     o 212 MMFF94 atom type symbols are listed
#     o 95 MMFF94 atom type numbers are listed
#     o Atom type numbers from 83 to 86 are not used
#     o Number of atom type symbols for:
#         o C: 34
#         o N: 47
#         o O: 45
#         o P: 7
#         o S: 18
#         o F, Br: 2
#         o Cl: 3
#         o I: 1
#         o H: 41
#         o Fe,Cu, Zn: 2
#         o Li, Na, L, K, Mg, Si, : 1
#
sub AssignAtomTypes {
  my($This) = @_;

  $This->_AssignAtomTypesToNonHydrogenAtoms();

  if (!$This->{IgnoreHydrogens}) {
    $This->_AssignAtomTypesToHydrogenAtoms();
  }

  return $This;
}

# Assign atom types to all non-Hydrogen atoms..
#
sub _AssignAtomTypesToNonHydrogenAtoms {
  my($This) = @_;
  my($Atom, $AtomType);

  ATOM: for $Atom ($This->GetMolecule()->GetAtoms()) {
    if ($Atom->IsHydrogen()) {
      next ATOM;
    }
    $AtomType = $This->_GetAtomType($Atom);
    $This->SetAtomType($Atom, $AtomType);
  }
  return $This;
}

# Assign atom types to Hydrogen atoms..
#
sub _AssignAtomTypesToHydrogenAtoms {
  my($This) = @_;
  my($Atom, $AtomType);

  if ($This->{IgnoreHydrogens}) {
    return $This;
  }

  ATOM: for $Atom ($This->GetMolecule()->GetAtoms()) {
    if (!$Atom->IsHydrogen()) {
      next ATOM;
    }
    $AtomType = $This->_GetAtomTypeForHydrogen($Atom);
    $This->SetAtomType($Atom, $AtomType);
  }
  return $This;
}

# Get MMFF94 atom type for atom...
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

# Get MMFF94 atom type for Carbon atom...
#
# 34 AtomTypeSymbols for element C:
#
# AtomTypeSymbol   AtomTypeNum   AtomTypeDefinition
#   CR       1     ALKYL CARBON, SP3
#   C=C      2     VINYLIC CARBON, SP2
#   CSP2     2     GENERIC SP2 CARBON
#   C=O      3     GENERAL CARBONYL CARBON
#   C=N      3     SP2 CARBON IN C=N
#   CGD      3     GUANIDINE CARBON, DOUBLY BONDED TO N
#   C=OR     3     KETONE OR ALDEHYDE CARBONYL CARBON
#   C=ON     3     AMIDE CARBONYL CARBON
#   CONN     3     UREA CARBONYL CARBON
#   COO      3     CARBOXYLIC ACID OR ESTER CARBONYL CARBON
#   COON     3     CARBAMATE CARBONYL CARBON
#   COOO     3     CARBONIC ACID OR ESTER CARBONYL CARBON
#   C=OS     3     THIOESTER CARBONYL CARBON, DOUBLE BONDED TO O
#   C=S      3     THIOESTER CARBON, DOUBLY BONDED TO S
#   C=SN     3     THIOAMIDE, CARBON, DOUBLY BONDED TO S
#   CSO2     3     CARBON IN >C=SO2
#   CS=O     3     CARBON IN >C=S=O (SULFINYL GROUP)
#   CSS      3     THIOCARBOXYLIC ACID OR ESTER CARBONYL CARBON
#   C=P      3     CARBON DOUBLE BONDED TO PHOSPHOROUS
#   CSP      4     ACETYLENIC CARBON
#   =C=      4     ALLENIC CARBON
#   CR4R     20    CARBON IN 4-MEMBERED RINGS
#   CR3R     22    CARBON IN A 3-MEMBERED RING
#   CE4R     30    OLEFINIC CARBON IN 4-MEMBERED RINGS
#   CB       37    CARBON AS IN BENZENE, PYRROLE
#   CO2M     41    CARBOXYLATE ANION CARBON
#   CS2M     41    CARBON IN THIOCARBOXYLATE ANION
#   CGD+     57    GUANIDINIUM CARBON
#   CNN+     57    C IN +N=C-N RESONANCE STRUCTURES
#   C%       60    ISONITRILE CARBON
#   C5A      63    ALPHA CARBON IN 5-MEMBERED HETEROAROMATIC RING
#   C5B      64    BETA CARBON IN 5-MEMBERED HETEROAROMATIC RING
#   C5       78    GENERAL CARBON IN 5-MEMBERED HETEROAROMATIC RING
#   CIM+     80    C IN N-C-N IN IMIDAZOLIUM ION
#
# Notes:
#     . During atom type assignments, matches are performed starting from specific to generic.
#
sub _GetAtomTypeForCarbon {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'None';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {

    # Only single bonds...
    if ($NumOfPiBonds == 0) {
      $AtomType = $This->_GetAtomTypeForCarbonWithOnlySigmaBonds($Atom);
      last ATOMTYPE;
    }

    # One double bond...
    if ($NumOfPiBonds == 1) {
      $AtomType = $This->_GetAtomTypeForCarbonWithOnePiBond($Atom);
      last ATOMTYPE;
    }

    # One triple bond or two double bonds...
    if ($NumOfPiBonds == 2) {
      $AtomType = $This->_GetAtomTypeForCarbonWithTwoPiBonds($Atom);
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForCarbon: MMFF94 atom type for Carbon cann't be assigned...";
  }
  return $AtomType;
}

# Get MMFF94 atom type for Nitrogen atom...
#
# 47 AtomTypeSymbols for element N:
#
# AtomTypeSymbol   AtomTypeNum   AtomTypeDefinition
#   NR       8     NITROGEN IN ALIPHATIC AMINES
#   N=C      9     NITROGEN IN IMINES
#   N=N      9     NITROGEN IN AZO COMPOUNDS
#   NC=O     10    NITROGEN IN AMIDES
#   NC=S     10    NITROGEN IN N-C=S, THIOAMIDE
#   NN=C     10    NITROGEN IN N-N=C
#   NN=N     10    NITROGEN IN N-N=N
#   NR+      34    QUATERNARY NITROGEN, SP3, POSITIVELY CHARGED
#   NPYD     38    NITROGEN, AS IN PYRIDINE
#   NPYL     39    NITROGEN, AS IN PYRROLE
#   NC=C     40    NITROGEN ON N-C=C
#   NC=N     40    NITROGEN IN N-C=N
#   NC=P     40    NITROGEN IN N-C=P
#   NC%C     40    NITROGEN ATTACHED TO C-C TRIPLE BOND
#   NSP      42    NITROGEN, TRIPLE BONDED
#   NSO2     43    NITROGEN IN SULFONAMIDES
#   NSO3     43    NITROGEN IN SULFONAMIDES, THREE Os ON S
#   NPO2     43    NITROGEN IN PHOSPHONAMIDES
#   NPO3     43    NITROGEN IN PHOSPHONAMIDES, THREE Os ON P
#   NC%N     43    NITROGEN ATTACHED TO CYANO GROUP
#   NO2      45    NITRO GROUP NITROGEN
#   NO3      45    NITRATE GROUP NITROGEN
#   N=O      46    NITROSO NITROGEN
#   NAZT     47    TERMINAL NITROGEN IN AZIDO OR DIAZO GROUP
#   NSO      48    DIVALENT NITROGEN REPLACING MONOVALENT O IN SO2 GROUP
#   =N=      53    NITROGEN IN C=N=N OR -N=N=N
#   N+=C     54    POSITIVELY CHARGED IMINIUM NITROGEN
#   N+=N     54    POSITIVELY CHARGED NITROGEN DOUBLE-BONDED TO N
#   NCN+     55    N IN +N=C-N RESONANCE STRUCTURES - FORMAL CHARGE=1/2
#   NGD+     56    GUANIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1/3
#   NPD+     58    PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1
#   NR%      61    ISONITRILE NITROGEN [FC = 0] OR DIAZO NITROGEN [FC = 1]
#   NM       62    DEPROTONATED SULFONAMIDE N-; FORMAL CHARGE=-1
#   N5A      65    ALPHA AROM HETEROCYCLIC 5-RING  NITROGEN
#   N5B      66    BETA AROM HETEROCYCLIC 5-RING  NITROGEN
#   N2OX     67    SP2-HYDRIDIZED N-OXIDE NITROGEN
#   N3OX     68    SP3-HYDRIDIZED N-OXIDE NITROGEN
#   NPOX     69    PYRIDINE N-OXIDE NITROGEN
#   N5M      76    NEGATIVELY CHARGED N IN, E.G, TRI- OR TETRAZOLE ANION
#   N5       79    GENERAL NITROGEN IN 5-MEMBERED HETEROCYCLIC RING
#   NIM+     81    IMIDAZOLIUM-TYPE NITROGEN - FORMAL CHARGE=1/2
#   N5A+     81    POSITIVE N5A NITROGEN - FORMAL CHARGE=1
#   N5B+     81    POSITIVE N5B NITROGEN - FORMAL CHARGE=1
#   N5+      81    POSITIVE N5 NITROGEN - FORMAL CHARGE=1
#   N5AX     82    N-OXIDE NITROGEN IN 5-RING ALPHA POSITION
#   N5BX     82    N-OXIDE NITROGEN IN 5-RING BETA POSITION
#   N5OX     82    N-OXIDE NITROGEN IN GENERAL 5-RING POSITION
#
# Notes:
#   . The current release of MayaChemTools assigns "None" to Oxygens in the following environment
#     as no generic or specific MMFF94 atom types exists to handle them:
#
#      . Terminal Nitrogens attched to Sulfur in >S=N
#
sub _GetAtomTypeForNitrogen {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'None';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  ATOMTYPE: {

    # Nitrogens in five membered rings...
    if ($Atom->IsInRingOfSize(5)) {
      $AtomType = $This->_GetAtomTypeForFiveMemberedRingNitrogen($Atom);
      last ATOMTYPE;
    }

    # -N(-)-, -N+(-)(-)-, -(N-1)(-)
    if ($NumOfPiBonds == 0) {
      $AtomType = $This->_GetAtomTypeForNitrogenWithOnlySigmaBonds($Atom);
      last ATOMTYPE;
    }

    # -N=, and -N+(=)-
    if ($NumOfPiBonds == 1) {
      $AtomType = $This->_GetAtomTypeForNitrogenWithOnePiBond($Atom);
      last ATOMTYPE;
    }

    # #N, #N+-, and =N+=
    if ($NumOfPiBonds == 2) {
      $AtomType = $This->_GetAtomTypeForNitrogenWithTwoPiBonds($Atom);
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForNitrogen: MMFF94 atom type for Nitrogen cann't be assigned...";
  }
  return $AtomType;
}

# Get MMFF94 atom type for Oxygen atom...
#
# 45 AtomTypeSymbols for element O:
#
# AtomTypeSymbol   AtomTypeNum   AtomTypeDefinition
#   OR       6     ALCOHOL OR ETHER OXYGEN
#   OC=O     6     ESTER OR CARBOXYLIC ACID -O-
#   OC=C     6     ENOLIC OR PHENOLIC OXYGEN
#   OC=N     6     DIVALENT OXYGEN
#   OC=S     6     THIOESTER OR THIOACID -O-
#   ONO2     6     DIVALENT NITRATE ETHER OXYGEN
#   ON=O     6     DIVALENT NITRITE ETHER OXYGEN
#   OSO3     6     DIVALENT OXYGEN ATTACHED TO SULFUR
#   OSO2     6     DIVALENT OXYGEN ATTACHED TO SULFUR
#   OSO      6     DIVALENT OXYGEN ATTACHED TO SULFUR
#   OS=O     6     DIVALENT OXYGEN ATTACHED TO SULFOXIDE SULFUR
#   -OS      6     GENERAL DIVALENT OXYGEN ATTACHED TO S
#   OPO3     6     DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
#   OPO2     6     DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
#   OPO      6     DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
#   -OP      6     DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
#   -O-      6     GENERAL DIVALENT O
#   O=C      7     GENERAL C=O
#   O=CN     7     CARBONYL OXYGEN, AMIDES
#   O=CR     7     CARBONYL OXYGEN, ALDEHYDES AND KETONES
#   O=CO     7     CARBONYL OXYGEN, CARBOXYLIC ACIDS AND ESTERS
#   O=N      7     NITROSO OXYGEN
#   O=S      7     O=S IN SULFOXIDES
#   O=S=     7     O=S ON SULFUR DOUBLY BONDED TO, E.G., CARBON
#   O2CM     32    OXYGEN IN CARBOXYLATE ANION
#   OXN      32    N-OXIDE OXYGEN
#   O2N      32    NITRO OXYGEN
#   O2NO     32    NITRO-GROUP OXYGEN IN NITRATE
#   O3N      32    NITRATE ANION OXYGEN
#   O-S      32    SINGLE TERMINAL OXYGEN ON TETRACOORD SULFUR
#   O2S      32    TERMINAL O-S IN SULFONES AND SULFONAMIDES
#   O3S      32    TERMINAL O IN SULFONATES
#   O4S      32    TERMINAL O IN SO4(-3)
#   OSMS     32    TERM O IN THIOSULFINATE ANION - FORMAL CHARGE=-0.5
#   OP       32    TERMINAL O IN PHOSPHOXIDES
#   O2P      32    TERMINAL O IN PHOSPHINATES
#   O3P      32    TERMINAL OXYGEN IN PHOSPHONATES
#   O4P      32    TERMINAL OXYGEN IN PHOSPHATES AND PHOSPHODIESTERS
#   O4CL     32    OXYGEN IN CLO4(-) ANION - FORMAL CHARGE=-0.25
#   OM       35    ALKOXIDE OXYGEN, NEGATIVELY CHARGED
#   OM2      35    OXIDE OXYGEN ON SP2 CARBON, NEGATIVELY CHARGED
#   O+       49    POSITIVELY CHARGED OXONIUM (TRICOORDINATE) OXYGEN
#   O=+      51    POSITIVELY CHARGED OXENIUM (DICOORDINATE) OXYGEN
#   OFUR     59    AROMATIC OXYGEN AS IN FURAN
#   OH2      70    OXYGEN ON WATER
#
# Notes:
#   . The current release of MayaChemTools assigns "None" to Oxygens in the following environment
#     as no generic or specific MMFF94 atom types exists to handle them:
#
#      . Terminal anion Oxygen corresponding to divalent Oxygen attached to Sulfoxide Sulfur in
#        OS=O and divalent Oxygen attached to Sulfur in -OS
#      . Terminal Oxygens attched to Sulfur in =SO2
#      . Terminal anion Oxygen attched to Sulfur in SO2M which is same as OS=O
#
sub _GetAtomTypeForOxygen {
  my($This, $Atom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds, $OxygenAttachedToSulfur, $OxygenAttachedToPhosphorus);

  $AtomType = 'None';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $Atom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $Atom->GetAtomicInvariantValue('H');

  $OxygenAttachedToSulfur = $Atom->GetNeighborsUsingAtomSpecification('S');
  $OxygenAttachedToPhosphorus = $Atom->GetNeighborsUsingAtomSpecification('P');

  ATOMTYPE: {

   # Divalent or terminal Oxygen attached to Sulfur...
    if ($OxygenAttachedToSulfur) {
      $AtomType = $This->_GetAtomTypeForOxygenAttachedToSulfur($Atom);
      last ATOMTYPE;
    }

   # Divalent or terminal Oxygen attached to Phosphorous...
    if ($OxygenAttachedToPhosphorus) {
      $AtomType = $This->_GetAtomTypeForOxygenAttachedToPhosphorus($Atom);
      last ATOMTYPE;
    }

    # Only single bonds...
    if ($NumOfPiBonds == 0) {
      $AtomType = $This->_GetAtomTypeForOxygenWithOnlySigmaBonds($Atom);
      last ATOMTYPE;
    }

    # One double bond...
    if ($NumOfPiBonds == 1) {
      $AtomType = $This->_GetAtomTypeForOxygenWithOnePiBond($Atom);
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForOxygen: MMFF94 atom type for Oxygen cann't be assigned...";
  }
  return $AtomType;
}

# Get MMFF94 atom type for Phosphorus atom...
#
# 7 AtomTypeSymbols for element P:
#
# AtomTypeSymbol   AtomTypeNum   AtomTypeDefinition
#   PO4      25    PHOSPHOROUS IN PHOSPHATES AND PHOSPHODIESTERS
#   PO3      25    TETRACOORDINATE P WITH THREE ATTACHED OXYGENS
#   PO2      25    TETRACOORDINATE P WITH TWO ATTACHED OXYGENS
#   PO       25    TETRACOORDINATE P WITH ONE ATTACHED OXYGEN
#   PTET     25    GENERAL TETRACOORDINATE PHOSPHORUS
#   P        26    TRICOORDINATE P, AS IN PHOSPHINES
#   -P=C     75    PHOSPHOROUS DOUBLY BONDED TO CARBON
#
sub _GetAtomTypeForPhosphorus {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # PO4 : PHOSPHOROUS IN PHOSPHATES AND PHOSPHODIESTERS
    if ($This->_IsPhosphateOrPhosphodiesterPhosphorus($Atom)) {
      $AtomType = 'PO4';
      last ATOMTYPE;
    }

    # PO3 : TETRACOORDINATE P WITH THREE ATTACHED OXYGENS
    if ($This->_IsPhosphonatePhosphorus($Atom)) {
      $AtomType = 'PO3';
      last ATOMTYPE;
    }

    # PO2 : TETRACOORDINATE P WITH TWO ATTACHED OXYGENS
    if ($This->_IsPhosphinatePhosphorus($Atom)) {
      $AtomType = 'PO2';
      last ATOMTYPE;
    }

    # PO : TETRACOORDINATE P WITH ONE ATTACHED OXYGEN
    if ($This->_IsPhosphoxidePhosphorus($Atom)) {
      $AtomType = 'PO';
      last ATOMTYPE;
    }

    # -P=C : PHOSPHOROUS DOUBLY BONDED TO CARBON
    if ($This->_IsDoublyBondedToCarbonPhosphorous($Atom)) {
      $AtomType = '-P=C';
      last ATOMTYPE;
    }

    # PTET : GENERAL TETRACOORDINATE PHOSPHORUS
    if ($This->_IsTetraCoordinatedPhosphorus($Atom)) {
      $AtomType = 'PTET';
      last ATOMTYPE;
    }

    # P : TRICOORDINATE P, AS IN PHOSPHINES
    if ($This->_IsTriCoordinatedPhosphorus($Atom)) {
      $AtomType = 'P';
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForPhosphorus: MMFF94 atom type for Phosphorous cann't be assigned...";
  }

  return $AtomType;
}

# Get MMFF94 atom type for Sulfur atom...
#
# 18 AtomTypeSymbols for element S:
#
# AtomTypeSymbol   AtomTypeNum   AtomTypeDefinition
#   S        15    SULFUR IN THIOETHERS AND MERCAPTANS
#   S=C      16    TERMINAL SULFUR DOUBLY BONDED TO CARBON
#   S=O      17    SULFUR IN SULFOXIDES
#   >S=N     17    SULFUR, TRICOORD, DOUBLY BONDED TO N
#   SO2      18    SULFUR IN SULFONES
#   SO2N     18    SULFUR IN SULFONAMIDES
#   SO3      18    SULFONATE SULFUR
#   SO4      18    SULFATE SULFUR
#   =SO2     18    SULFONE SULPHER DOUBLY BONDED TO CARBON
#   SNO      18    SULFUR IN NITROGEN ANALOG OF A SULFONE
#   STHI     44    SULFUR AS IN THIOPHENE
#   S-P      72    TERMINAL SULFUR BONDED TO PHOSPHORUS
#   S2CM     72    TERMINAL SULFUR IN THIOCARBOXYLATE ANION
#   SM       72    TERMINAL SULFUR - FORMAL CHARGE=-1
#   SSMO     72    TERMINAL SULFUR IN THIOSULFINATE GROUP
#   SO2M     73    SULFUR IN NEGATIVELY CHARGED SULFINATE GROUP
#   SSOM     73    TRICOORD SULFUR IN THIOSULFINATE GROUP
#   =S=O     74    SULFINYL SULFUR, EG. IN C=S=O
#
sub _GetAtomTypeForSulfur {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # SO4 : SULFATE SULFUR
    if ($This->_IsSulfateSulfur($Atom)) {
      $AtomType = 'SO4';
      last ATOMTYPE;
    }

    # SO3 : SULFONATE SULFUR
    if ($This->_IsSulfonateSulfur($Atom)) {
      $AtomType = 'SO3';
      last ATOMTYPE;
    }

    # SO2N : SULFUR IN SULFONAMIDES
    if ($This->_IsSulfonamideSulfur($Atom)) {
      $AtomType = 'SO2N';
      last ATOMTYPE;
    }

    # SO2 : SULFUR IN SULFONES
    if ($This->_IsSulfoneSulfur($Atom)) {
      $AtomType = 'SO2';
      last ATOMTYPE;
    }

    #  =SO2: SULFONE SULPHER DOUBLY BONDED TO CARBON
    if ($This->_IsDoublyBondedToCarbonSulfoneSulfur($Atom)) {
      $AtomType = '=SO2';
      last ATOMTYPE;
    }

    # SO2M: SULFUR IN NEGATIVELY CHARGED SULFINATE GROUP
    if ($This->_IsNegativelyChargedSulfinateSulfur($Atom)) {
      $AtomType = 'SO2M';
      last ATOMTYPE;
    }

    # SNO : SULFUR IN NITROGEN ANALOG OF A SULFONE
    if ($This->_IsNitrogenAnalogOfSulfoneSulfur($Atom)) {
      $AtomType = 'SNO';
      last ATOMTYPE;
    }

    # S=O : SULFUR IN SULFOXIDES
    if ($This->_IsSulfoxideSulfur($Atom)) {
      $AtomType = 'S=O';
      last ATOMTYPE;
    }

    # >S=N : SULFUR, TRICOORD, DOUBLY BONDED TO N
    if ($This->_IsSNTricoordinatedSulfur($Atom)) {
      $AtomType = '>S=N';
      last ATOMTYPE;
    }

    # STHI : SULFUR AS IN THIOPHENE
    if ($This->_IsSTHISulfur($Atom)) {
      $AtomType = 'STHI';
      last ATOMTYPE;
    }

    # S2CM : TERMINAL SULFUR IN THIOCARBOXYLATE ANION
    if ($This->_IsThioCarboxylateAnionTerminalSulfur($Atom)) {
      $AtomType = 'S2CM';
      last ATOMTYPE;
    }

    # SSMO : TERMINAL SULFUR IN THIOSULFINATE GROUP
    if ($This->_IsThioSulfinateTerminalSulfur($Atom)) {
      $AtomType = 'SSMO';
      last ATOMTYPE;
    }

    # SSOM : TRICOORD SULFUR IN THIOSULFINATE GROUP
    if ($This->_IsTriCoordinatedThioSulfinateSulfur($Atom)) {
      $AtomType = 'SSOM';
      last ATOMTYPE;
    }

    #   =S=O:  SULFINYL SULFUR, EG. IN C=S=O
    if ($This->_IsSulfinylSulfur($Atom)) {
      $AtomType = '=S=O';
      last ATOMTYPE;
    }

    # S-P : TERMINAL SULFUR BONDED TO PHOSPHORUS
    if ($This->_IsSPTerminalSulfur($Atom)) {
      $AtomType = 'S-P';
      last ATOMTYPE;
    }

    # S=C : TERMINAL SULFUR DOUBLY BONDED TO CARBON
    if ($This->_IsSCTerminalSulfur($Atom)) {
      $AtomType = 'S=C';
      last ATOMTYPE;
    }

    # SM : TERMINAL SULFUR - FORMAL CHARGE=-1
    if ($This->_IsNegativelyChargedTerminalSulfur($Atom)) {
      $AtomType = 'SM';
      last ATOMTYPE;
    }

    # S : SULFUR IN THIOETHERS AND MERCAPTANS
    if ($This->_IsThioEthersOrMercaptansSulfur($Atom)) {
      $AtomType = 'S';
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForSulfur: MMFF94 atom type for Sulfur cann't be assigned...";
  }
  return $AtomType;
}

# Get MMFF94 atom type for Hydrogen atom...
#
# 41 AtomTypeSymbols for element H:
#
# AtomTypeSymbol   AtomTypeNum   AtomTypeDefinition
#   HC       5     H  ATTACHED TO C
#   HSI      5     H ATTACHED TO SI
#   HOR      21    HYDROGEN IN ALCOHOLS
#   HO       21    GENERAL H ON OXYGEN
#   HOM      21    HYDROGEN IN HYDROXIDE ANION
#   HNR      23    H-N(SP3)
#   H3N      23    H-N(SP3), AMMONIA
#   HPYL     23    H-N IN PYRROLE
#   HNOX     23    H-N IN IN A N-OXIDE
#   HNM      23    H ON DICOORD, NEGATIVELY CHARGED NITROGEN
#   HN       23    GENERAL H ON NITROGEN
#   HOCO     24    H-O IN CARBOXYLIC ACIDS
#   HOP      24    HYDROGEN ON OXYGEN ATTACHED TO PHOSPHOROUS
#   HN=N     27    AZO HYDROGEN
#   HN=C     27    IMINE HYDROGEN
#   HNCO     28    AMIDE HYDROGEN
#   HNCS     28    THIOAMIDE HYDROGEN
#   HNCC     28    H-N IN ENAMINES
#   HNCN     28    H-N IN H-N-C=N
#   HNNC     28    H-N IN H-N-N=C
#   HNNN     28    H-N IN H-N-N=N
#   HNSO     28    H-N IN SULFONAMIDE
#   HNPO     28    H-N IN PHOSPHONAMIDE
#   HNC%     28    HYDROGEN ON N ATTACHED TO TRIPLY BONDED CARBON
#   HSP2     28    GENERAL H ON SP2 NITROGEN
#   HOCC     29    H-O IN ENOLS AND PHENOLS
#   HOCN     29    H-O IN HO-C=N
#   HOH      31    HYDROGEN IN H2O
#   HOS      33    H ON OXYGEN ATTACHED TO SULFUR
#   HNR+     36    H ON QUATERNARY NITROGEN
#   HIM+     36    H ON IMIDAZOLIUM-TYPE NITROGEN
#   HPD+     36    H ON PROTONATED PYRIDINE NITROGEN
#   HNN+     36    H ON AMIDINIUM-TYPE NITROGEN
#   HNC+     36    H ON PROTONATED IMINE NITROGEN
#   HGD+     36    H ON GUANIDINIUM-TYPE NITROGEN
#   HN5+     36    H ON N5+, N5A+ OR N5B+
#   HO+      50    HYDROGEN ON O+ OXYGEN
#   HO=+     52    HYDROGEN ON OXENIUM OXYGEN
#   HS       71    H ATTACHED TO DIVALENT, DICOORDINATE S
#   HS=N     71    H ATTACHED TO TETRAVALENT, TRICOODR S DBL BONDED TO N
#   HP       71    H ATTACHED TO TRI- OR TETRACOORDINATE PHOSPHORUS
#
sub _GetAtomTypeForHydrogen {
  my($This, $Atom) = @_;
  my($AtomType, $AtomNeighbor);

  $AtomType = 'None';

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
    if ($AtomNeighbor->IsPhosphorus()) {
      $AtomType = $This->_GetAtomTypeForHydrogenAttachedToPhosphorus($AtomNeighbor);
      last ATOMNEIGHBOR;
    }
    if ($AtomNeighbor->IsSulfur()) {
      $AtomType = $This->_GetAtomTypeForHydrogenAttachedToSulfur($AtomNeighbor);
      last ATOMNEIGHBOR;
    }
    if ($AtomNeighbor->IsSilicon()) {
      $AtomType = $This->_GetAtomTypeForHydrogenAttachedToSilicon($AtomNeighbor);
      last ATOMNEIGHBOR;
    }
    $AtomType = "None";
    carp "Warning: ${ClassName}->_GetAtomTypeForHydrogen: MMFF94 atom type for Hydrogen cann't be assigned...";
  }
  return $AtomType;
}

# Get MMFF94 atom type for atoms other than Carbon, Nitrogen, Oxygen, Phosporus,
# Sulfur and Hydrogen...
#
# AtomTypeSymbol   AtomTypeNum   AtomTypeDefinition
#   LI+      92    LITHIUM CATION
#   F        11    FLUORINE
#   F-       89    FLUORIDE ANION
#   NA+      93    SODIUM CATION
#   MG+2     99    DIPOSITIVE MAGNESIUM CATION
#   SI       19    SILICON
#   CL       12    CHLORINE
#   CLO4     77    CHLORINE IN PERCHLORATE ANION, CLO4(-)
#   CL-      90    CHLORIDE ANION
#   K+       94    POTASSIUM CATION
#   CA+2     96    DIPOSITIVE CALCIUM
#   FE+2     87    IRON +2 CATION
#   FE+3     88    IROM +3 CATION
#   CU+1     97    MONOPOSITIVE COPPER
#   CU+2     98    DIPOSITIVE COPPER
#   ZINC     95    DIPOSITIVE ZINC
#   ZN+2     95    DIPOSITIVE ZINC
#   BR       13    BROMINE
#   BR-      91    BROMIDE ANION
#   I        14    IODINE
#
#
sub _GetAtomTypeForOtherAtoms {
  my($This, $Atom) = @_;
  my($AtomType, $AtomSymbol, $FormalCharge, $CallingMethod, @AllowedFormalCharges);

  $AtomType = 'None';

  $AtomSymbol = $Atom->GetAtomSymbol();
  $FormalCharge = $Atom->GetFormalCharge();

  $CallingMethod = "_GetAtomTypeForOtherAtoms";
  @AllowedFormalCharges = ();

  ATOMSYMBOL: {
    # FLUORINE
    if ($AtomSymbol =~ /^F$/i) {
      #  F : FLUORINE;  F- :  FLUORIDE ANION
      $AtomType = ($FormalCharge == -1) ? 'F-' : 'F';
      @AllowedFormalCharges = ('0', '-1');
      last ATOMSYMBOL;
    }
    # CHLORINE
    if ($AtomSymbol =~ /^Cl$/i) {
      # CL : CHLORINE; CLO4 : CHLORINE IN PERCHLORATE ANION, CLO4(-);
      # CL- : CHLORIDE ANION
      $AtomType = $This->_IsPerChlorateAnionChlorine($Atom) ? 'CLO4' : (($FormalCharge == -1) ? 'CL-' : 'CL');
      @AllowedFormalCharges = ('0', '-1');
      last ATOMSYMBOL;
    }
    # BROMINE
    if ($AtomSymbol =~ /^Br$/i) {
      # BR : BROMINE; BR- : BROMIDE ANION
      $AtomType = ($FormalCharge == -1) ? 'BR-' : 'BR';
      @AllowedFormalCharges = ('0', '-1');
      last ATOMSYMBOL;
    }
    # IODINE
    if ($AtomSymbol =~ /^I$/i) {
      $AtomType = 'I';
      @AllowedFormalCharges = ('0');
      last ATOMSYMBOL;
    }
    # LI+ : LITHIUM CATION
    if ($AtomSymbol =~ /^Li$/i) {
      $AtomType = 'LI+';
      @AllowedFormalCharges = ('+1');
      last ATOMSYMBOL;
    }
    # NA+ : SODIUM CATION
    if ($AtomSymbol =~ /^Na$/i) {
      $AtomType = 'NA+';
      @AllowedFormalCharges = ('+1');
      last ATOMSYMBOL;
    }
    # MG+2 : DIPOSITIVE MAGNESIUM CATION
    if ($AtomSymbol =~ /^Mg$/i) {
      $AtomType = 'MG+2';
      @AllowedFormalCharges = ('+2');
      last ATOMSYMBOL;
    }
    # SI : SILICON
    if ($AtomSymbol =~ /^Si$/i) {
      $AtomType = 'SI';
      @AllowedFormalCharges = ('0');
      last ATOMSYMBOL;
    }
    # K+ : POTASSIUM CATION
    if ($AtomSymbol =~ /^K$/i) {
      $AtomType = 'K+';
      @AllowedFormalCharges = ('+1');
      last ATOMSYMBOL;
    }
    # CA+2 : DIPOSITIVE CALCIUM
    if ($AtomSymbol =~ /^Ca$/i) {
      $AtomType = 'CA+2';
      @AllowedFormalCharges = ('+2');
      last ATOMSYMBOL;
    }
    # IRON
    if ($AtomSymbol =~ /^Fe$/i) {
      # FE+2 : IRON +2 CATION; FE+3 : IROM +3 CATION
      $AtomType = ($FormalCharge == 3) ? 'FE+3' : 'FE+2';
      @AllowedFormalCharges = ('+2', '+3');
      last ATOMSYMBOL;
    }
    # COPPER
    if ($AtomSymbol =~ /^Cu$/i) {
      # CU+1 : MONOPOSITIVE COPPER; CU+2 : DIPOSITIVE COPPER
      $AtomType = ($FormalCharge == 1) ? 'CU+1' : 'CU+2';
      @AllowedFormalCharges = ('+1', '+2');
      last ATOMSYMBOL;
    }
    # ZINC
    if ($AtomSymbol =~ /^Zn$/i) {
      # ZN+2 : DIPOSITIVE ZINC
      $AtomType = 'Zn+2';
      @AllowedFormalCharges = ('+2');
      last ATOMSYMBOL;
    }
    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForOtherAtoms: MMFF94 atom type for $AtomSymbol cann't be assigned...";
  }
  if (@AllowedFormalCharges) {
    $This->_CheckFormalChargeMismatch($CallingMethod, $AtomSymbol, $AtomType, $FormalCharge, \@AllowedFormalCharges);
  }
  return $AtomType;
}

# Check any formal charge mismatches...
#
sub _CheckFormalChargeMismatch {
  my($This, $CallingMethod, $AtomSymbol, $AssignedAtomType, $FormalCharge, $AllowedFormalChargesRef) = @_;
  my($AllowedFormalCharge, $FormalChargeMismatch);

  $FormalChargeMismatch = 1;

  FORMALCHARGE: for $AllowedFormalCharge (@{$AllowedFormalChargesRef}) {
    if ($AllowedFormalCharge == $FormalCharge) {
      $FormalChargeMismatch = 0;
      last FORMALCHARGE;
    }
  }
  if ($FormalChargeMismatch) {
    my($AllowedFormalCharges);
    $AllowedFormalCharges = TextUtil::JoinWords($AllowedFormalChargesRef, ",", 0);

    carp "\nWarning: ${ClassName}->${CallingMethod}:_CheckFormalChargeMismatch: MMFF94 atom for $AtomSymbol with formal charge, $FormalCharge, cann't be assigned: Formal charge, $FormalCharge, is different from allowed formal charge(s): $AllowedFormalCharges. Default UFF atom type, $AssignedAtomType, has been assigned...";
  }

  return $This;
}

# Get MMFF94 atom type for Carbon with only sigma bonds...
#
sub _GetAtomTypeForCarbonWithOnlySigmaBonds {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {
    # CR4R : CARBON IN 4-MEMBERED RINGS
    if ($This->_IsFourMemberedRingCarbon($Atom)) {
      $AtomType = 'CR4R';
      last ATOMTYPE;
    }

    # CR3R : CARBON IN A 3-MEMBERED RING
    if ($This->_IsThreeMemberedRingCarbon($Atom)) {
      $AtomType = 'CR3R';
      last ATOMTYPE;
    }

    #  CR : ALKYL CARBON, SP3
    if ($This->_IsAlkylCarbon($Atom)) {
      $AtomType = 'CR';
      last ATOMTYPE;
    }

    # As far as I can tell, MMFF94 doesn't have a generic atom type for SP3 carbon. So the current
    # release of MayaChemTools package used CR as defaul SP3 carbon.
    #
    $AtomType = 'CR';
  }

  return $AtomType;
}

# Get MMFF94 atom type for Carbon with one pi bond...
#
sub _GetAtomTypeForCarbonWithOnePiBond {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {
    # CR3R :  CARBON IN 3-MEMBERED RINGS
    if ($This->_IsThreeMemberedRingOlefinicCarbon($Atom)) {
      $AtomType = 'CR3R';
      last ATOMTYPE;
    }

    # CE4R :  OLEFINIC CARBON IN 4-MEMBERED RINGS
    if ($This->_IsFourMemberedRingOlefinicCarbon($Atom)) {
      $AtomType = 'CE4R';
      last ATOMTYPE;
    }

    # CIM+ : C IN N-C-N IN IMIDAZOLIUM ION
    if ($This->_IsImidazoliumCarbon($Atom)) {
      $AtomType = 'CIM+';
      last ATOMTYPE;
    }

    # C5A : ALPHA CARBON IN 5-MEMBERED HETEROAROMATIC RING
    if ($This->_IsFiveMemberedHeteroAromaticRingAlphaCarbon($Atom)) {
      $AtomType = 'C5A';
      last ATOMTYPE;
    }

    # C5B : BETA CARBON IN 5-MEMBERED HETEROAROMATIC RING
    if ($This->_IsFiveMemberedHeteroAromaticRingBetaCarbon($Atom)) {
      $AtomType = 'C5B';
      last ATOMTYPE;
    }

    # C5 : GENERAL CARBON IN 5-MEMBERED HETEROAROMATIC RING
    if ($This->_IsFiveMemberedHeteroAromaticRingCarbon($Atom)) {
      $AtomType = 'C5';
      last ATOMTYPE;
    }

    # CB : CARBON AS IN BENZENE, PYRROLE
    if ($This->_IsRingAromaticCarbon($Atom)) {
      $AtomType = 'CB';
      last ATOMTYPE;
    }

    # C=C : VINYLIC CARBON, SP2
    if ($This->_IsVinylicCarbon($Atom)) {
      $AtomType = 'C=C';
      last ATOMTYPE;
    }

    # C=OR : KETONE OR ALDEHYDE CARBONYL CARBON
    if ($This->_IsKetoneOrAldehydeCarbonylCarbon($Atom)) {
      $AtomType = 'C=OR';
      last ATOMTYPE;
    }

    # C=ON : AMIDE CARBONYL CARBON
    if ($This->_IsAmideCarbonylCarbon($Atom)) {
      $AtomType = 'C=ON';
      last ATOMTYPE;
    }

    # CONN : UREA CARBONYL CARBON
    if ($This->_IsUreaCarbonylCarbon($Atom)) {
      $AtomType = 'CONN';
      last ATOMTYPE;
    }

    # COO : CARBOXYLIC ACID OR ESTER CARBONYL CARBON
    if ($This->_IsCarboxylicAcidOrEsterCarbonylCarbon($Atom)) {
      $AtomType = 'COO';
      last ATOMTYPE;
    }

    # CO2M : CARBOXYLATE ANION CARBON
    if ($This->_IsCarboxylateAnionCarbon($Atom)) {
      $AtomType = 'CO2M';
      last ATOMTYPE;
    }

    # COON: CARBAMATE CARBONYL CARBON
    if ($This->_IsCarbamateCarbonylCarbon($Atom)) {
      $AtomType = 'COON';
      last ATOMTYPE;
    }

    # COOO: CARBONIC ACID OR ESTER CARBONYL CARBON
    if ($This->_IsCarbonicAcidOrEsterCarbonylCarbon($Atom)) {
      $AtomType = 'COOO';
      last ATOMTYPE;
    }

    # C=OS : THIOESTER CARBONYL CARBON, DOUBLE BONDED TO O
    if ($This->_IsThioEsterCarbonylCarbon($Atom)) {
      $AtomType = 'C=OS';
      last ATOMTYPE;
    }

    # C=O : GENERAL CARBONYL CARBON
    if ($This->_IsGeneralCarbonylCarbon($Atom)) {
      $AtomType = 'C=O';
      last ATOMTYPE;
    }

    # C=SN : THIOAMIDE, CARBON, DOUBLY BONDED TO S
    if ($This->_IsThioAmideCarbon($Atom)) {
      $AtomType = 'C=SN';
      last ATOMTYPE;
    }

    # CSO2 : CARBON IN >C=SO2
    if ($This->_IsSulfonylCarbon($Atom)) {
      $AtomType = 'CSO2';
      last ATOMTYPE;
    }

    # CS=O :  CARBON IN >C=S=O (SULFINYL GROUP)
    if ($This->_IsSulfinylCarbon($Atom)) {
      $AtomType = 'CS=O';
      last ATOMTYPE;
    }

    # CSS : THIOCARBOXYLIC ACID OR ESTER CARBONYL CARBON
    if ($This->_IsThioCarboxylicAcidOrEsterCarbonylCarbon($Atom)) {
      $AtomType = 'CSS';
      last ATOMTYPE;
    }

    # CS2M : CARBON IN THIOCARBOXYLATE ANION
    if ($This->_IsThioCarboxylateAnionCarbon($Atom)) {
      $AtomType = 'CS2M';
      last ATOMTYPE;
    }

    # C=S : THIOESTER CARBON, DOUBLY BONDED TO S
    if ($This->_IsThioEsterCarbon($Atom)) {
      $AtomType = 'C=S';
      last ATOMTYPE;
    }

    # C=P : CARBON DOUBLE BONDED TO PHOSPHOROUS
    if ($This->_IsDoublyBondedToPhosphorousCarbon($Atom)) {
      $AtomType = 'C=P';
      last ATOMTYPE;
    }

    # CGD : GUANIDINE CARBON, DOUBLY BONDED TO N
    if ($This->_IsGuandineCarbon($Atom)) {
      $AtomType = 'CGD';
      last ATOMTYPE;
    }

    # CGD+ : GUANIDINIUM CARBON
    if ($This->_IsGuandiniumCarbon($Atom)) {
      $AtomType = 'CGD+';
      last ATOMTYPE;
    }

    #CNN+ : C IN +N=C-N RESONANCE STRUCTURES
    if ($This->_IsNCNResonaceStructuresCarbon($Atom)) {
      $AtomType = 'CNN+';
      last ATOMTYPE;
    }

    # C=N : SP2 CARBON IN C=N
    if ($This->_IsDoublyBondedToNitrogenSP2Carbon($Atom)) {
      $AtomType = 'C=N';
      last ATOMTYPE;
    }

    # CSP2 : GENERIC SP2 CARBON
    if ($This->_IsGenericSP2Carbon($Atom)) {
      $AtomType = 'CSP2';
      last ATOMTYPE;
    }

    $AtomType = 'CSP2';
    carp "Warning: ${ClassName}->_GetAtomTypeForCarbonWithOnePiBond: MMFF94 atom type for Carbon cann't be assigned; Default MMFF94 atom type, $AtomType, has been assigned...\n";
  }

  return $AtomType;
}

# Get MMFF94 atom type for Carbon with two pi bonds...
#
sub _GetAtomTypeForCarbonWithTwoPiBonds {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {
    # =C= : ALLENIC CARBON
    if ($This->_IsAllenicCarbon($Atom)) {
      $AtomType = '=C=';
      last ATOMTYPE;
    }

    # C% : ISONITRILE CARBON ( R-N+#C- )
    if ($This->_IsIsoNitrileCarbon($Atom)) {
      $AtomType = 'C%';
      last ATOMTYPE;
    }

    # CSP : ACETYLENIC CARBON
    if ($This->_IsAcetylenicCarbon($Atom)) {
      $AtomType = 'CSP';
      last ATOMTYPE;
    }

    $AtomType = 'CSP';
    carp "Warning: ${ClassName}->_GetAtomTypeForCarbonWithTwoPiBonds: MMFF94 atom type for Carbon cann't be assigned; Default MMFF94 atom type, $AtomType, has been assigned...\n";
  }

  return $AtomType;
}


# CR4R : CARBON IN 4-MEMBERED RINGS
#
sub _IsFourMemberedRingCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.RA4.T4.TSB4') ? 1 : 0;
}

# CR4R : CARBON IN 3-MEMBERED RINGS
#
sub _IsThreeMemberedRingCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.RA3.T4.TSB4') ? 1 : 0;
}

#  CR : ALKYL CARBON, SP3
#
sub _IsAlkylCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T4.TSB4', ['C,H', 'C,H', 'C,H', 'C,H'], ['-', '-', '-', '-']) ? 1 : 0;
}

# C5A : ALPHA CARBON IN 5-MEMBERED HETEROAROMATIC RING
#
sub _IsFiveMemberedHeteroAromaticRingAlphaCarbon {
  my($This, $Atom) = @_;

  # Is it an aromatic atom in a five membered ring?
  if (!$Atom->DoesAtomNeighborhoodMatch('C.Ar.RA5')) {
    return 0;
  }

  # Is it part of a five membered ring containing hetero atom at alpha position?
  my($RingAtomsRef);

  for $RingAtomsRef ($Atom->GetRingsWithSize(5)) {
    if ($This->_IsAtomPositionAlphaInHeteroAromaticRing($Atom, $RingAtomsRef)) {
      return 1;
    }
  }
  return 0;
}

# C5B : BETA CARBON IN 5-MEMBERED HETEROAROMATIC RING
#
sub _IsFiveMemberedHeteroAromaticRingBetaCarbon {
  my($This, $Atom) = @_;

  # Is it an aromatic atom in a five membered ring?
  if (!$Atom->DoesAtomNeighborhoodMatch('C.Ar.RA5')) {
    return 0;
  }

  # Is it part of five membered rings containing hetero atom at beta position?
  my($RingAtomsRef);

  for $RingAtomsRef ($Atom->GetRingsWithSize(5)) {
    if ($This->_IsAtomPositionBetaInHeteroAromaticRing($Atom, $RingAtomsRef)) {
      return 1;
    }
  }
  return 0;
}


# C5 : GENERAL CARBON IN 5-MEMBERED HETEROAROMATIC RING
#
sub _IsFiveMemberedHeteroAromaticRingCarbon {
  my($This, $Atom) = @_;

  # Is it an aromatic atom in a five membered ring?
  if (!$Atom->DoesAtomNeighborhoodMatch('C.Ar.RA5')) {
    return 0;
  }

  # Is it part of five membered rings containing at least one hetero atom?
  my($RingAtomsRef, $RingIsAromatic, $NumOfHeteroAtoms);
  for $RingAtomsRef ($Atom->GetRingsWithSize(5)) {
    ($RingIsAromatic, $NumOfHeteroAtoms) = $This->_GetHeteroAtomsInformationInRing($RingAtomsRef);
    if ($RingIsAromatic && $NumOfHeteroAtoms >= 1) {
      return 1;
    }
  }
  return 0;
}

# CR3R :  OLEFINIC CARBON IN 3-MEMBERED RINGS
#
# Notes:
#    . MMFF94 atom type for olefinic Carbon in 3-membered rings is same as SP3 Carbon.
#
sub _IsThreeMemberedRingOlefinicCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.RA3.T3.DB1', ['*', '*', '*'], ['=', '-', '-']) ? 1 : 0;
}

# CE4R :  OLEFINIC CARBON IN 4-MEMBERED RINGS
#
sub _IsFourMemberedRingOlefinicCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.RA4.T3.DB1', ['*', '*', '*'], ['=', '-', '-']) ? 1 : 0;
}

# CB : CARBON AS IN BENZENE, PYRROLE
#
# Notes:
#    . MayaChemTools assigns all aromatic carbons, other than five membered
#      hetero aromatic ring Carbons, as CB.
#
sub _IsRingAromaticCarbon {
  my($This, $Atom) = @_;

  return ($Atom->DoesAtomNeighborhoodMatch('C.Ar.RA')) ? 1 : 0;
}

# C=C : VINYLIC CARBON, SP2
#
sub _IsVinylicCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['C', '*', '*'], ['=', '-', '-'], ['C,H', 'C,H', 'C,H']) ? 1 : 0;
}

# C=OR : KETONE OR ALDEHYDE CARBONYL CARBON
#
sub _IsKetoneOrAldehydeCarbonylCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['O', 'C,H', 'C,H'], ['=', '-', '-']) ? 1 : 0;
}

# C=ON : AMIDE CARBONYL CARBON
#
sub _IsAmideCarbonylCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['O', 'N', 'C,H'], ['=', '-', '-']) ? 1 : 0;
}

# CONN : UREA CARBONYL CARBON : R-(R'-)N-C(=O)-N(-R")-R'''
#
sub _IsUreaCarbonylCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.X3.DB1', ['O', 'N', 'N'], ['=', '-', '-']) ? 1 : 0;
}

# COO : CARBOXYLIC ACID OR ESTER CARBONYL CARBON
#
sub _IsCarboxylicAcidOrEsterCarbonylCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['O', 'O.X1.FC0,O.X2', 'C,H'], ['=', '-', '-'], ['C', 'C,H', 'C,H']) ? 1 : 0;
}

# CO2M : CARBOXYLATE ANION CARBON
#
sub _IsCarboxylateAnionCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['O', 'O.X1.FC-1', 'C,H'], ['=', '-', '-'], ['C', 'C,H', 'C,H']) ? 1 : 0;
}

# COON: CARBAMATE CARBONYL CARBON : R-O-C(=O)-N(-R')-R"
#
sub _IsCarbamateCarbonylCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.X3.DB1', ['O', 'O', 'N'], ['=', '-', '-']) ? 1 : 0;
}

# COOO: CARBONIC ACID OR ESTER CARBONYL CARBON:  R-O-C(=O)-O-R'
#
sub _IsCarbonicAcidOrEsterCarbonylCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.X3.DB1', ['O', 'O', 'O'], ['=', '-', '-']) ? 1 : 0;
}

# C=OS : THIOESTER CARBONYL CARBON, DOUBLE BONDED TO O
#
sub _IsThioEsterCarbonylCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['O', 'S', 'C,H'], ['=', '-', '-']) ? 1 : 0;
}

# C=O : GENERAL CARBONYL CARBON
#
sub _IsGeneralCarbonylCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['O', '*', '*'], ['=', '-', '-']) ? 1 : 0;
}

# C=SN : THIOAMIDE, CARBON, DOUBLY BONDED TO S
#
sub _IsThioAmideCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['S', 'N', 'C,H'], ['=', '-', '-']) ? 1 : 0;
}

# CSO2 : CARBON IN >C=SO2
#
sub _IsSulfonylCarbon {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  # Is it connected to a sulfone Sulfur?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S')) {
    if ($This->_IsDoublyBondedToCarbonSulfoneSulfur($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}


# CS=O :  CARBON IN >C=S=O (SULFINYL GROUP)
#
sub _IsSulfinylCarbon {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  # Is it connected to a sulfinyl Sulfur?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S')) {
    if ($This->_IsSulfinylSulfur($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# CSS : THIOCARBOXYLIC ACID OR ESTER CARBONYL CARBON
#
sub _IsThioCarboxylicAcidOrEsterCarbonylCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['S', 'S.X1.FC0,S.X2', 'C,H'], ['=', '-', '-'], ['C', 'C,H', 'C,H']) ? 1 : 0;
}

# CS2M : CARBON IN THIOCARBOXYLATE ANION
#
sub _IsThioCarboxylateAnionCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['S', 'S.X1.FC-1', 'C,H'], ['=', '-', '-'], ['C', 'C,H', 'C,H']) ? 1 : 0;
}

# C=S : THIOESTER CARBON, DOUBLY BONDED TO S
#
sub _IsThioEsterCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['S.X1', '*', '*'], ['=', '-', '-']) ? 1 : 0;
}

# C=P : CARBON DOUBLE BONDED TO PHOSPHOROUS
#
sub _IsDoublyBondedToPhosphorousCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['P.X1', '*', '*'], ['=', '-', '-']) ? 1 : 0;
}

# CIM+ : C IN N-C-N IN IMIDAZOLIUM ION
#
# Notes:
#    . Imidazole is five membered ring contain  C*=C-N-C=N* (* - Ring bond)
#    . Imidazolium ion is where N in C=N has a formal charge of +1 and has an extra substituent
#
sub _IsImidazoliumCarbon {
  my($This, $Atom) = @_;

  # Is it in a 5 membered aromatic ring surrounded by two aromatic Nitrogens?
  #
  if (!($Atom->DoesAtomNeighborhoodMatch('C.Ar.RA5.TR1', ['N.Ar.RA5.TR1.FC0', 'N.Ar.RA5.TR1.FC+1'], ['-', '=']) ||
       $Atom->DoesAtomNeighborhoodMatch('C.Ar.RA5.TR1', ['N.Ar.RA5.TR1.!FC0', 'N.Ar.RA5.TR1.!FC0'], ['-', '=']))) {
    return 0;
  }

  # Check to make sure ring contains only two Nitrogen hetero atom...
  my($RingAtomsRef, $RingIsAromatic, $NumOfHeteroAtoms, $HeteroAtomSymbolsRef);

  ($RingAtomsRef) = $Atom->GetRingsWithSize(5);
  ($RingIsAromatic, $NumOfHeteroAtoms, $HeteroAtomSymbolsRef) = $This->_GetHeteroAtomsInformationInRing($RingAtomsRef);

  if ($NumOfHeteroAtoms != 2) {
    return 0;
  }

  return (exists($HeteroAtomSymbolsRef->{N}) && ($HeteroAtomSymbolsRef->{N} == 2)) ? 1 : 0;
}

# CNN+ : C IN +N=C-N RESONANCE STRUCTURES
#
sub _IsNCNResonaceStructuresCarbon {
  my($This, $Atom) = @_;

  # Match +1 formal charge on Nitrogen...
  if ($Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['N.FC+1', 'N.FC0'], ['=', '-'])) {
    return 1;
  }

  # Match non-zero (+1/2) formal charge on Nitrogen with single, double or resonnace bonds...
  if ($Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['N.!FC0', 'N.!FC0'], ['-,=,:', '-,=,:'])) {
    return 1;
  }

  return 0;
}

# CGD : GUANIDINE CARBON, DOUBLY BONDED TO N
#
sub _IsGuandineCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.X3.DB1', ['N.FC0', 'N.FC0', 'N.FC0'], ['=', '-', '-']) ? 1 : 0;
}

# CGD+ : GUANIDINIUM CARBON
#
sub _IsGuandiniumCarbon {
  my($This, $Atom) = @_;

  # Match +1 formal charge on a Nitrogen with explicitly single and double bonds...
  if ($Atom->DoesAtomNeighborhoodMatch('C.X3.DB1', ['N.FC+1', 'N.FC0', 'N.FC0'], ['=', '-', '-'])) {
    return 1;
  }

  # Match +1/3 formal charge on each Nitrogen with single, double or resonance bonds...
  if ($Atom->DoesAtomNeighborhoodMatch('C.X3.DB1', ['N.!FC0', 'N.!FC0', 'N.!FC0'], ['-,=,:', '-,=,:', '-,=,:'])) {
    return 1;
  }

  return 0;
}

# C=N : SP2 CARBON IN C=N
#
sub _IsDoublyBondedToNitrogenSP2Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['N', '*', '*'], ['=', '-', '-']) ? 1 : 0;
}

# CSP2 : GENERIC SP2 CARBON
#
sub _IsGenericSP2Carbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['*', '*', '*'], ['=', '-', '-']) ? 1 : 0;
}

# =C= : ALLENIC CARBON
#
sub _IsAllenicCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.X2.DB2', ['*', '*'], ['=', '=']) ? 1 : 0;
}

# C% : ISONITRILE CARBON ( R-N+#C- )
#
sub _IsIsoNitrileCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.X1.TB1.FC-1', ['N.T2.TB1.FC+1'], ['#'], ['C,H']) ? 1 : 0;
}

# CSP : ACETYLENIC CARBON
#
sub _IsAcetylenicCarbon {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('C.T2.TB1', ['*', '*'], ['#', '-']) ? 1 : 0;
}

# Get MMFF94 atom type Nitrogen in five membered ring...
#
# Note:
#    . The method handles both SP3 and SP2 five membered ring Nitrogens.
#
sub _GetAtomTypeForFiveMemberedRingNitrogen {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # NIM+ : IMIDAZOLIUM-TYPE NITROGEN - FORMAL CHARGE=1/2
    if ($This->_IsImidazoliumNitrogen($Atom)) {
      $AtomType = 'NIM+';
      last ATOMTYPE;
    }

    # N5A+ : POSITIVE N5A NITROGEN - FORMAL CHARGE=1
    if ($This->_IsPositivelyChargedFiveMemberedHeteroAromaticRingAlphaNitrogen($Atom)) {
      $AtomType = 'N5A+';
      last ATOMTYPE;
    }

    # N5A : ALPHA AROM HETEROCYCLIC 5-RING  NITROGEN
    if ($This->_IsFiveMemberedHeteroAromaticRingAlphaNitrogen($Atom)) {
      $AtomType = 'N5A';
      last ATOMTYPE;
    }

    # N5B+ : POSITIVE N5B NITROGEN - FORMAL CHARGE=1
    if ($This->_IsPositivelyChargedFiveMemberedHeteroAromaticRingBetaNitrogen($Atom)) {
      $AtomType = 'N5B+';
      last ATOMTYPE;
    }

    # N5B : BETA AROM HETEROCYCLIC 5-RING  NITROGEN
    if ($This->_IsFiveMemberedHeteroAromaticRingBetaNitrogen($Atom)) {
      $AtomType = 'N5B';
      last ATOMTYPE;
    }

    # N5AX : N-OXIDE NITROGEN IN 5-RING ALPHA POSITION
    if ($This->_IsNOxideFiveMemberedHeteroCyclicRingAlphaNitrogen($Atom)) {
      $AtomType = 'N5AX';
      last ATOMTYPE;
    }

    # N5BX : N-OXIDE NITROGEN IN 5-RING BETA POSITION
    if ($This->_IsNOxideFiveMemberedHeteroCyclicRingBetaNitrogen($Atom)) {
      $AtomType = 'N5BX';
      last ATOMTYPE;
    }

    # N5OX : N-OXIDE NITROGEN IN GENERAL 5-RING POSITION
    if ($This->_IsNOxideFiveMemberedHeteroCyclicRingNitrogen($Atom)) {
      $AtomType = 'N5OX';
      last ATOMTYPE;
    }

    # N5M : NEGATIVELY CHARGED N IN, E.G, TRI- OR TETRAZOLE ANION
    if ($This->_IsNegativelyChargedFiveMemberedHeteroCyclicRingNitrogen($Atom)) {
      $AtomType = 'N5M';
      last ATOMTYPE;
    }

    # N5+ : POSITIVE N5 NITROGEN - FORMAL CHARGE=1
    if ($This->_IsPositivelyChargedFiveMemberedHeteroCyclicRingNitrogen($Atom)) {
      $AtomType = 'N5+';
      last ATOMTYPE;
    }

    # N5 : GENERAL NITROGEN IN 5-MEMBERED HETEROCYCLIC RING
    if ($This->_IsFiveMemberedHeteroCyclicRingNitrogen($Atom)) {
      $AtomType = 'N5';
      last ATOMTYPE;
    }

    $AtomType = 'N5';
  }

  return $AtomType;
}

# Get MMFF94 atom type for Nitrogen with only sigma bonds...
#
sub _GetAtomTypeForNitrogenWithOnlySigmaBonds {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # NPYL : NITROGEN, AS IN PYRROLE
    if ($This->_IsPyrroleNitrogen($Atom)) {
      $AtomType = 'NPYL';
      last ATOMTYPE;
    }

    # NC=O : NITROGEN IN AMIDES
    if ($This->_IsAmideNitrogen($Atom)) {
      $AtomType = 'NC=O';
      last ATOMTYPE;
    }

    # NC=S : NITROGEN IN N-C=S, THIOAMIDE
    if ($This->_IsThioAmideNitrogen($Atom)) {
      $AtomType = 'NC=S';
      last ATOMTYPE;
    }

    # NN=C : NITROGEN IN N-N=C
    if ($This->_IsNNCNitrogen($Atom)) {
      $AtomType = 'NN=C';
      last ATOMTYPE;
    }

    # NGD+ : GUANIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1/3
    if ($This->_IsGuanidiniumNitrogen($Atom)) {
      $AtomType = 'NGD+';
      last ATOMTYPE;
    }

    # NCN+ :  N IN +N=C-N RESONANCE STRUCTURES - FORMAL CHARGE=1/2
    if ($This->_IsNCNResonaceStructuresNitrogen($Atom)) {
      $AtomType = 'NCN+';
      last ATOMTYPE;
    }

    # NN=N : NITROGEN IN N-N=N
    if ($This->_IsNNNNitrogen($Atom)) {
      $AtomType = 'NN=N';
      last ATOMTYPE;
    }

    #  NC=C : NITROGEN ON N-C=C
    if ($This->_IsNCCNitrogen($Atom)) {
      $AtomType = 'NC=C';
      last ATOMTYPE;
    }

    # NC=N : NITROGEN IN N-C=N
    if ($This->_IsNCNNitrogen($Atom)) {
      $AtomType = 'NC=N';
      last ATOMTYPE;
    }

    # NC=P : NITROGEN IN N-C=P
    if ($This->_IsNCPNitrogen($Atom)) {
      $AtomType = 'NC=P';
      last ATOMTYPE;
    }

    # NM : DEPROTONATED SULFONAMIDE N-; FORMAL CHARGE=-1
    if ($This->_IsDeprotonatedSulfonamideNitrogen($Atom)) {
      $AtomType = 'NM';
      last ATOMTYPE;
    }

    # NSO2 : NITROGEN IN SULFONAMIDES
    if ($This->_IsNSO2SulfonamideNitrogen($Atom)) {
      $AtomType = 'NSO2';
      last ATOMTYPE;
    }

    # NSO3 : NITROGEN IN SULFONAMIDES, THREE Os ON S
    if ($This->_IsNSO3SulfonamideNitrogen($Atom)) {
      $AtomType = 'NSO3';
      last ATOMTYPE;
    }

    # NPO2 : NITROGEN IN PHOSPHONAMIDES
    if ($This->_IsNPO2PhosphonamideNitrogen($Atom)) {
      $AtomType = 'NPO2';
      last ATOMTYPE;
    }

    # NPO3 : NITROGEN IN PHOSPHONAMIDES, THREE Os ON P
    if ($This->_IsNPO3PhosphonamideNitrogen($Atom)) {
      $AtomType = 'NPO3';
      last ATOMTYPE;
    }

    #  NC%N : NITROGEN ATTACHED TO CYANO GROUP
    if ($This->_IsAttachedToCyanoNitrogen($Atom)) {
      $AtomType = 'NC%N';
      last ATOMTYPE;
    }

    # N3OX : SP3-HYDRIDIZED N-OXIDE NITROGEN
    if ($This->_IsSP3NOxideNitrogen($Atom)) {
      $AtomType = 'N3OX';
      last ATOMTYPE;
    }

    #  NC%C : NITROGEN ATTACHED TO C-C TRIPLE BOND
    if ($This->_IsAttchedToCCTripleBondNitrogen($Atom)) {
      $AtomType = 'NC%C';
      last ATOMTYPE;
    }

    # NR : NITROGEN IN ALIPHATIC AMINES
    if ($This->_IsAliphaticAmineNitrogen($Atom)) {
      $AtomType = 'NR';
      last ATOMTYPE;
    }

    # NR+ : QUATERNARY NITROGEN, SP3, POSITIVELY CHARGED
    if ($This->_IsAliphaticAmineQuaternaryNitrogen($Atom)) {
      $AtomType = 'NR+';
      last ATOMTYPE;
    }

    $AtomType = 'NR';
  }

  return $AtomType;
}

# Get MMFF94 atom type for Nitrogen with one pi bond...
#
sub _GetAtomTypeForNitrogenWithOnePiBond {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # NPOX : PYRIDINE N-OXIDE NITROGEN
    if ($This->_IsPyridineNOxideNitrogen($Atom)) {
      $AtomType = 'NPOX';
      last ATOMTYPE;
    }

    # NPD+ : PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1
    if ($This->_IsPyridiniumNitrogen($Atom)) {
      $AtomType = 'NPD+';
      last ATOMTYPE;
    }

    # NPYD : NITROGEN, AS IN PYRIDINE
    if ($This->_IsPyridineNitrogen($Atom)) {
      $AtomType = 'NPYD';
      last ATOMTYPE;
    }

    # N2OX : SP2-HYDRIDIZED N-OXIDE NITROGEN
    if ($This->_IsSP2NOxideNitrogen($Atom)) {
      $AtomType = 'N2OX';
      last ATOMTYPE;
    }

    # NAZT : TERMINAL NITROGEN IN AZIDO OR DIAZO GROUP
    if ($This->_IsNAZTNitrogen($Atom)) {
      $AtomType = 'NAZT';
      last ATOMTYPE;
    }

    # NR% :  ISONITRILE NITROGEN [FC = 0] OR DIAZO NITROGEN [FC = 1]
    if ($This->_IsIsoNitrileOrDiAzoNitrogen($Atom)) {
      $AtomType = 'NR%';
      last ATOMTYPE;
    }

    # NGD+ : GUANIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1/3
    if ($This->_IsGuanidiniumNitrogen($Atom)) {
      $AtomType = 'NGD+';
      last ATOMTYPE;
    }

    # NCN+ :  N IN +N=C-N RESONANCE STRUCTURES - FORMAL CHARGE=1/2
    if ($This->_IsNCNResonaceStructuresNitrogen($Atom)) {
      $AtomType = 'NCN+';
      last ATOMTYPE;
    }

    # N=C : NITROGEN IN IMINES
    if ($This->_IsImineNitrogen($Atom)) {
      $AtomType = 'N=C';
      last ATOMTYPE;
    }

    # N+=C : POSITIVELY CHARGED IMINIUM NITROGEN
    if ($This->_IsPositivelyChargedIminiumNitrogen($Atom)) {
      $AtomType = 'N+=C';
      last ATOMTYPE;
    }

    # N=N :  NITROGEN IN AZO COMPOUNDS
    if ($This->_IsAzoNitrogen($Atom)) {
      $AtomType = 'N=N';
      last ATOMTYPE;
    }

    # N+=N : POSITIVELY CHARGED NITROGEN DOUBLE-BONDED TO N
    if ($This->_IsPositivelyChargedAzoNitrogen($Atom)) {
      $AtomType = 'N+=N';
      last ATOMTYPE;
    }

    # NO2 : NITRO GROUP NITROGEN
    if ($This->_IsNitroNitrogen($Atom)) {
      $AtomType = 'NO2';
      last ATOMTYPE;
    }

    # NO3 : NITRATE GROUP NITROGEN
    if ($This->_IsNitrateNitrogen($Atom)) {
      $AtomType = 'NO3';
      last ATOMTYPE;
    }

    # N=O : NITROSO NITROGEN
    if ($This->_IsNitrosoNitrogen($Atom)) {
      $AtomType = 'N=O';
      last ATOMTYPE;
    }

    # NSO : DIVALENT NITROGEN REPLACING MONOVALENT O IN SO2 GROUP
    if ($This->_IsNSONitrogen($Atom)) {
      $AtomType = 'NSO';
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForNitrogenWithOnePiBond: MMFF94 atom type for Nitrogen cann't be assigned; Default MMFF94 atom type, $AtomType, has been assigned...\n";
  }

  return $AtomType;
}

# Get MMFF94 atom type for Nitrogen with two pi bonds...
#
sub _GetAtomTypeForNitrogenWithTwoPiBonds {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # NAZT : TERMINAL NITROGEN IN AZIDO OR DIAZO GROUP
    if ($This->_IsNAZTNitrogen($Atom)) {
      $AtomType = 'NAZT';
      last ATOMTYPE;
    }

    # NR% :  ISONITRILE NITROGEN [FC = 0] OR DIAZO NITROGEN [FC = 1]
    if ($This->_IsIsoNitrileOrDiAzoNitrogen($Atom)) {
      $AtomType = 'NR%';
      last ATOMTYPE;
    }

    #  =N= : NITROGEN IN C=N=N OR -N=N=N
    if ($This->_IsTwoDoublyBondedNitrogen($Atom)) {
      $AtomType = '=N=';
      last ATOMTYPE;
    }

    # NSP : NITROGEN, TRIPLE BONDED
    if ($This->_IsTriplyBondedSPNitrogen($Atom)) {
      $AtomType = 'NSP';
      last ATOMTYPE;
    }

    $AtomType = 'NSP';
  }

  return $AtomType;
}

# NR : NITROGEN IN ALIPHATIC AMINES
#
sub _IsAliphaticAmineNitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.T3', ['C,H', 'C,H', 'C,H'], ['-', '-','-']) ? 1 : 0;
}

# NR+ : QUATERNARY NITROGEN, SP3, POSITIVELY CHARGED
#
sub _IsAliphaticAmineQuaternaryNitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.T4.FC+1', ['C,H', 'C,H', 'C,H', 'C,H'], ['-', '-','-','-']) ? 1 : 0;
}

# NC=O : NITROGEN IN AMIDES
#
sub _IsAmideNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T3.TSB3')) {
    return 0;
  }

  # Is it attached to amide carbonyl Carbon?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.T3.DB1')) {
    if ($This->_IsAmideCarbonylCarbon($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# NC=S : NITROGEN IN N-C=S, THIOAMIDE
#
sub _IsThioAmideNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T3.TSB3')) {
    return 0;
  }

  # Is it attached to thioamide Carbon?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.T3.DB1')) {
    if ($This->_IsThioAmideCarbon($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# NN=C : NITROGEN IN N-N=C
#
sub _IsNNCNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T3.TSB3')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('N.DB1')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('N.T2.DB1', ['C.T3', 'N.T3'], ['=', '-'])) {
      return 1;
    }
  }
  return 0;
}

# NN=N : NITROGEN IN N-N=N
#
sub _IsNNNNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T3.TSB3')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('N.DB1')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('N.T2.DB1', ['N.T2', 'N.T3'], ['=', '-'])) {
      return 1;
    }
  }
  return 0;
}

# NPYD : NITROGEN, AS IN PYRIDINE
#
sub _IsPyridineNitrogen {
  my($This, $Atom) = @_;

  # Is it an aromatic Nitrogen in only one six membered ring?
  if (!$Atom->DoesAtomNeighborhoodMatch('N.Ar.RA6.TR1')) {
    return 0;
  }

  # Is it part of six membered ring containing only one hetero atom?
  my($RingAtomsRef, $RingIsAromatic, $NumOfHeteroAtoms);
  for $RingAtomsRef ($Atom->GetRingsWithSize(6)) {
    ($RingIsAromatic, $NumOfHeteroAtoms) = $This->_GetHeteroAtomsInformationInRing($RingAtomsRef);
    if ($RingIsAromatic && $NumOfHeteroAtoms == 1) {
      return 1;
    }
  }
  return 0;
}

# NPYL : NITROGEN, AS IN PYRROLE
#
sub _IsPyrroleNitrogen {
  my($This, $Atom) = @_;

  # Is it an aromatic Nitrogen in only one five membered ring?
  if (!$Atom->DoesAtomNeighborhoodMatch('N.Ar.RA5.TR1')) {
    return 0;
  }

  # Is it part of five membered ring containing only one hetero atom?
  my($RingAtomsRef, $RingIsAromatic, $NumOfHeteroAtoms);
  for $RingAtomsRef ($Atom->GetRingsWithSize(5)) {
    ($RingIsAromatic, $NumOfHeteroAtoms) = $This->_GetHeteroAtomsInformationInRing($RingAtomsRef);
    if ($RingIsAromatic && $NumOfHeteroAtoms == 1) {
      return 1;
    }
  }
  return 0;
}

#  NC=C : NITROGEN ON N-C=C
#
sub _IsNCCNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T3.TSB3')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.DB1')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('C.T3.DB1', ['C.T3', 'N.T3', '*'], ['=', '-','-'])) {
      return 1;
    }
  }
  return 0;
}

# NC=N : NITROGEN IN N-C=N
#
sub _IsNCNNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T3.TSB3')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.DB1')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('C.T3.DB1', ['N.T2', 'N.T3', '*'], ['=', '-','-'])) {
      return 1;
    }
  }
  return 0;
}

# NC=P : NITROGEN IN N-C=P
#
sub _IsNCPNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T3.TSB3')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.DB1')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('C.T3.DB1', ['P.X1', 'N.T3', '*'], ['=', '-', '-'])) {
      return 1;
    }
  }
  return 0;
}

#  NC%C : NITROGEN ATTACHED TO C-C TRIPLE BOND
#
sub _IsAttchedToCCTripleBondNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T3.TSB3')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.TB1')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('C.X2.TB1', ['C.T2', 'N.T3'], ['#', '-'])) {
      return 1;
    }
  }
  return 0;
}

# NSO2 : NITROGEN IN SULFONAMIDES
#
sub _IsNSO2SulfonamideNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T3.TSB3')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.T4.DB2')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('S.T4.DB2', ['N', 'O', 'O', '!O'], ['-', '=', '=', '-'])) {
      return 1;
    }
  }
  return 0;
}

# NSO3 : NITROGEN IN SULFONAMIDES, THREE Os ON S
#
sub _IsNSO3SulfonamideNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T3.TSB3')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.T4.DB2')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('S.X4.DB2', ['N', 'O', 'O', 'O'], ['-', '=', '=', '-'])) {
      return 1;
    }
  }
  return 0;
}

# NPO2 : NITROGEN IN PHOSPHONAMIDES
#
sub _IsNPO2PhosphonamideNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T3.TSB3')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('P.DB1')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('P.T4.DB1', ['N', 'O', 'O', '!O'], ['-', '=', '-', '-'])) {
      return 1;
    }
  }
  return 0;
}

# NPO3 : NITROGEN IN PHOSPHONAMIDES, THREE Os ON P
#
sub _IsNPO3PhosphonamideNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T3.TSB3')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('P.DB1')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('P.X4.DB1', ['N', 'O', 'O', 'O'], ['-', '=', '-', '-'])) {
      return 1;
    }
  }
  return 0;
}

#  NC%N : NITROGEN ATTACHED TO CYANO GROUP
#
sub _IsAttachedToCyanoNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T3.TSB3')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.TB1')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('C.T2.TB1', ['N', 'N'], ['-', '#'])) {
      return 1;
    }
  }
  return 0;
}

# NM : DEPROTONATED SULFONAMIDE N-; FORMAL CHARGE=-1
#
sub _IsDeprotonatedSulfonamideNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T2.FC-1.TSB2')) {
    return 0;
  }

  # Is it attached to sulfonamide Sulfur?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.T4.DB2')) {
    if ($This->_IsSulfonamideSulfur($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# N=C : NITROGEN IN IMINES
#
sub _IsImineNitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.T2.DB1', ['C', '*'], ['=', '-']) ? 1 : 0;
}

# N+=C : POSITIVELY CHARGED IMINIUM NITROGEN
#
sub _IsPositivelyChargedIminiumNitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.DB1.FC+1', ['C', '*', '*'], ['=', '-', '-']) ? 1 : 0;
}

# N=N :  NITROGEN IN AZO COMPOUNDS
#
sub _IsAzoNitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.T2.DB1', ['N', '*'], ['=', '-']) ? 1 : 0;
}

# N+=N : POSITIVELY CHARGED NITROGEN DOUBLE-BONDED TO N
#
sub _IsPositivelyChargedAzoNitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.DB1.FC+1', ['N', '*', '*'], ['=', '-', '-']) ? 1 : 0;
}

# NO2 : NITRO GROUP NITROGEN
#
sub _IsNitroNitrogen {
  my($This, $Atom) = @_;

  # R-N+(=O)-(O-)
  return $Atom->DoesAtomNeighborhoodMatch('N.T3.DB1.FC+1', ['O', 'O.FC-1', '!O'], ['=', '-', '-']) ? 1 : 0;
}

# NO3 : NITRATE GROUP NITROGEN
#
sub _IsNitrateNitrogen {
  my($This, $Atom) = @_;

  # (O-)-N+(=O)-(O-) or R-O-N+(=O)-(O-)
  #
  return $Atom->DoesAtomNeighborhoodMatch('N.T3.DB1.FC+1', ['O', 'O.FC-1', 'O.T2.FC0'], ['=', '-', '-']) ? 1 : 0;
}

# N=O : NITROSO NITROGEN
#
sub _IsNitrosoNitrogen {
  my($This, $Atom) = @_;

  # R-N=O
  return $Atom->DoesAtomNeighborhoodMatch('N.T2.DB1', ['O'], ['=']) ? 1 : 0;
}

# NAZT : TERMINAL NITROGEN IN AZIDO OR DIAZO GROUP
#
sub _IsNAZTNitrogen {
  my($This, $Atom) = @_;

  return ($This->_IsAzidoTerminalNitrogen($Atom) || $This->_IsDiazoTerminalNitrogen($Atom)) ? 1 : 0;
}

# TERMINAL NITROGEN IN AZIDO GROUP
#
sub _IsAzidoTerminalNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.X1.DB1')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('N.DB2')) {
    if ($This->_IsAzidoGroupMiddleNitrogen($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# TERMINAL NITROGEN IN  DIAZO GROUP
#
sub _IsDiazoTerminalNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.X1.TB1,N.X1.DB1')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('N.DB2,N.TB1')) {
    if ($This->_IsDiazoGroupMiddleNitrogen($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# Azido group: R-N=N+=N-
#
sub _IsAzidoGroupMiddleNitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.X2.DB2.FC+1', ['N.X1.FC-1', 'N.T2.FC0'], ['=', '=']) ? 1 : 0;
}

# Diazo group: R'-(C-)(-R")-(N+)#N or R'-C(-R")=(N+)=(N-)
#
sub _IsDiazoGroupMiddleNitrogen {
  my($This, $Atom) = @_;

  # Diazo group: R'-C(-R")=(N+)=(N-)
  if ($Atom->DoesAtomNeighborhoodMatch('N.X2.DB2.FC+1', ['C.T3.FC0', 'N.X1.FC-1'], ['=', '='])) {
    return 1;
  }

  # Diazo group: R'-(C-)(-R")-(N+)#N
  if ($Atom->DoesAtomNeighborhoodMatch('N.X2.TB1.FC+1', ['C.T3.FC-1', 'N.X1.FC0'], ['-', '#'])) {
    return 1;
  }

  return 0;
}

# NSO : DIVALENT NITROGEN REPLACING MONOVALENT O IN SO2 GROUP
#
sub _IsNSONitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.T2.DB1')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.DB2')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('S.DB2', ['N', 'O'], ['=', '='])) {
      return 1;
    }
  }
  return 0;
}

# NCN+ :  N IN +N=C-N RESONANCE STRUCTURES - FORMAL CHARGE=1/2
#
# Notes:
#    . This method is called invoked for both SP and SP2 Nitrogens to check and assign
#      NCN+ to both Nitrogens.
#
sub _IsNCNResonaceStructuresNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.DB1')) {
    if ($This->_IsNCNResonaceStructuresCarbon($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# NGD+ : GUANIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1/3
#
# Notes:
#    . This method is called invoked for both SP and SP2 Nitrogens to check and assign
#      NGD+ to all three Guanidinium Nitrogens.
#
sub _IsGuanidiniumNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N')) {
    return 0;
  }

  # Is it attached to Guandinium Carbon?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.DB1')) {
    if ($This->_IsGuandiniumCarbon($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# NPD+ : PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1
#
sub _IsPyridiniumNitrogen {
  my($This, $Atom) = @_;

  return ($This->_IsPyridineNitrogen($Atom) && $Atom->DoesAtomNeighborhoodMatch('N.FC+1')) ? 1 : 0;
}

# NR% :  ISONITRILE NITROGEN [FC = 0] OR DIAZO NITROGEN [FC = 1]
#
# Notes:
#    . This method is called invoked for both SP and SP2 Nitrogens to check and assign
#      NR%.
#
sub _IsIsoNitrileOrDiAzoNitrogen {
  my($This, $Atom) = @_;

  return ($This->_IsIsoNitrileNitrogen($Atom) || $This->_IsDiAzoNitrogen($Atom)) ? 1 : 0;
}

# Isonitrile nitrogen: R-N+#C-
#
# Notes:
#    . MayaChemTools assumes isonitrile Nitrogen to have a formal charge of +1; otherwise
#      it won't be an isonitrile group. It's not clear why it would be zero as indicated in
#      MMFF documentation.
#
sub _IsIsoNitrileNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('N.FC+1.TB1')) {
    return 0;
  }

  # Is it attached to isonitrile Carbon?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.TB1')) {
    if ($This->_IsIsoNitrileCarbon($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# Diazo group: R1-(C-)(-R2)-(N+)#N or R1-C(-R2)=(N+)=(N-)
#
sub _IsDiAzoNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  return $This->_IsDiazoGroupMiddleNitrogen($Atom) ? 1 : 0;
}

# N5A : ALPHA AROM HETEROCYCLIC 5-RING  NITROGEN
#
sub _IsFiveMemberedHeteroAromaticRingAlphaNitrogen {
  my($This, $Atom) = @_;

  # Is it an aromatic atom in a five membered ring?
  if (!$Atom->DoesAtomNeighborhoodMatch('N.Ar.RA5.FC0')) {
    return 0;
  }

  # Is it part of five membered rings containing hetero atom at alpha position?
  my($RingAtomsRef);

  for $RingAtomsRef ($Atom->GetRingsWithSize(5)) {
    if ($This->_IsAtomPositionAlphaInHeteroAromaticRing($Atom, $RingAtomsRef)) {
      return 1;
    }
  }
  return 0;
}

# N5B : BETA AROM HETEROCYCLIC 5-RING  NITROGEN
#
sub _IsFiveMemberedHeteroAromaticRingBetaNitrogen {
  my($This, $Atom) = @_;

  # Is it an aromatic atom in a five membered ring?
  if (!$Atom->DoesAtomNeighborhoodMatch('N.Ar.RA5.FC0')) {
    return 0;
  }

  # Is it part of five membered rings containing hetero atom at beta position?
  my($RingAtomsRef);

  for $RingAtomsRef ($Atom->GetRingsWithSize(5)) {
    if ($This->_IsAtomPositionBetaInHeteroAromaticRing($Atom, $RingAtomsRef)) {
      return 1;
    }
  }
  return 0;
}

# N5A+ : POSITIVE N5A NITROGEN - FORMAL CHARGE=1
#
sub _IsPositivelyChargedFiveMemberedHeteroAromaticRingAlphaNitrogen {
  my($This, $Atom) = @_;

  # Is it an aromatic atom in a five membered ring?
  if (!$Atom->DoesAtomNeighborhoodMatch('N.Ar.RA5.FC+1')) {
    return 0;
  }

  # Is it part of five membered rings containing hetero atom at alpha position?
  my($RingAtomsRef);

  for $RingAtomsRef ($Atom->GetRingsWithSize(5)) {
    if ($This->_IsAtomPositionAlphaInHeteroAromaticRing($Atom, $RingAtomsRef)) {
      return 1;
    }
  }
  return 0;
}

# N5B+ : POSITIVE N5B NITROGEN - FORMAL CHARGE=1
#
sub _IsPositivelyChargedFiveMemberedHeteroAromaticRingBetaNitrogen {
  my($This, $Atom) = @_;

  # Is it an aromatic atom in a five membered ring?
  if (!$Atom->DoesAtomNeighborhoodMatch('N.Ar.RA5.FC+1')) {
    return 0;
  }

  # Is it part of five membered rings containing hetero atom at beta position?
  my($RingAtomsRef);

  for $RingAtomsRef ($Atom->GetRingsWithSize(5)) {
    if ($This->_IsAtomPositionBetaInHeteroAromaticRing($Atom, $RingAtomsRef)) {
      return 1;
    }
  }
  return 0;
}

# N2OX : SP2-HYDRIDIZED N-OXIDE NITROGEN
#
sub _IsSP2NOxideNitrogen {
  my($This, $Atom) = @_;

  # R=(R'-)(N+)-(O-)
  #
  return $Atom->DoesAtomNeighborhoodMatch('N.FC+1.DB1.T3', ['O.X1.FC-1', 'C', '*'], ['-', '=', '-']) ? 1 : 0;
}

# N3OX : SP3-HYDRIDIZED N-OXIDE NITROGEN
#
sub _IsSP3NOxideNitrogen {
  my($This, $Atom) = @_;

  # R-(R'-)(R"-)(N+)-(O-)
  #
  return $Atom->DoesAtomNeighborhoodMatch('N.FC+1.T4', ['O.X1.FC-1', '*', '*', '*'], ['-', '-', '-', '-']) ? 1 : 0;
}

# NPOX : PYRIDINE N-OXIDE NITROGEN
#
sub _IsPyridineNOxideNitrogen {
  my($This, $Atom) = @_;

  return ($This->_IsPyridineNitrogen($Atom) && $This->_IsSP2NOxideNitrogen($Atom)) ? 1 : 0;
}

# N5M : NEGATIVELY CHARGED N IN, E.G, TRI- OR TETRAZOLE ANION
#
sub _IsNegativelyChargedFiveMemberedHeteroCyclicRingNitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.RA5.FC-1') ? 1 : 0;
}

# N5 : GENERAL NITROGEN IN 5-MEMBERED HETEROCYCLIC RING
#
sub _IsFiveMemberedHeteroCyclicRingNitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.RA5') ? 1 : 0;
}

# N5+ : POSITIVE N5 NITROGEN - FORMAL CHARGE=1
#
sub _IsPositivelyChargedFiveMemberedHeteroCyclicRingNitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.RA5.FC+1') ? 1 : 0;
}

# N5AX : N-OXIDE NITROGEN IN 5-RING ALPHA POSITION
#
sub _IsNOxideFiveMemberedHeteroCyclicRingAlphaNitrogen {
  my($This, $Atom) = @_;

  # Is it a N-oxide Nitrogen atom in a five membered ring?
  if (!$Atom->DoesAtomNeighborhoodMatch('N.RA5.FC+1', ['O.X1.FC-1.!RA'], ['-'])) {
    return 0;
  }

  # Is it part of five membered rings containing hetero atom at alpha position?
  my($RingAtomsRef);

  for $RingAtomsRef ($Atom->GetRingsWithSize(5)) {
    if ($This->_IsAtomPositionAlphaInHeteroRing($Atom, $RingAtomsRef)) {
      return 1;
    }
  }
  return 0;
}

# N5BX : N-OXIDE NITROGEN IN 5-RING BETA POSITION
#
sub _IsNOxideFiveMemberedHeteroCyclicRingBetaNitrogen {
  my($This, $Atom) = @_;

  # Is it a N-oxide Nitrogen atom in a five membered ring?
  if (!$Atom->DoesAtomNeighborhoodMatch('N.RA5.FC+1', ['O.X1.FC-1.!RA'], ['-'])) {
    return 0;
  }

  # Is it part of five membered rings containing hetero atom at alpha position?
  my($RingAtomsRef);

  for $RingAtomsRef ($Atom->GetRingsWithSize(5)) {
    if ($This->_IsAtomPositionBetaInHeteroRing($Atom, $RingAtomsRef)) {
      return 1;
    }
  }
  return 0;
}

# N5OX : N-OXIDE NITROGEN IN GENERAL 5-RING POSITION
#
sub _IsNOxideFiveMemberedHeteroCyclicRingNitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.FC+1.RA5', ['O.X1.FC-1.!RA'], ['-']) ? 1 : 0;
}

# NIM+ : IMIDAZOLIUM-TYPE NITROGEN - FORMAL CHARGE=1/2
#
# Notes:
#    . MayaChemTools assigns NIM+ to both the Nitrogens around Imidazolium Carbon.
#
sub _IsImidazoliumNitrogen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  # Is it an aromatic Nitrogen in only one five membered ring?
  if (!$Atom->DoesAtomNeighborhoodMatch('N.Ar.RA5.TR1')) {
    return 0;
  }

  # Is it attached to Imidazolium Carbon?
  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.RA5')) {
    if ($This->_IsImidazoliumCarbon($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# =N= : NITROGEN IN C=N=N OR -N=N=N
#
# Notes:
#    . MayaChemTools assign =N= to all doubly bonded Nitrogen atoms
#      irrespective of whether it has explcit formal charge of +1.
#
sub _IsTwoDoublyBondedNitrogen {
  my($This, $Atom) = @_;

  if ($Atom->DoesAtomNeighborhoodMatch('N.DB2', ['N', 'N'], ['=', '='])) {
    return 1;
  }
  if ($Atom->DoesAtomNeighborhoodMatch('N.DB2', ['C', 'N'], ['=', '='])) {
    return 1;
  }
  return 0;
}

# NSP : NITROGEN, TRIPLE BONDED
#
sub _IsTriplyBondedSPNitrogen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('N.TB1') ? 1 : 0;
}


# Get MMFF94 atom type for divalent or terminal Oxygen attached to Sulfur...
#
sub _GetAtomTypeForOxygenAttachedToSulfur {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # OSMS : TERM O IN THIOSULFINATE ANION - FORMAL CHARGE=-0.5
    if ($This->_IsThioSulfinateTerminalOxygen($Atom)) {
      $AtomType = 'OSMS';
      last ATOMTYPE;
    }

    # OSO3 : DIVALENT OXYGEN ATTACHED TO SULFUR
    if ($This->_IsOSO3DivalentOxygen($Atom)) {
      $AtomType = 'OSO3';
      last ATOMTYPE;
    }

    # OSO2 : DIVALENT OXYGEN ATTACHED TO SULFUR
    if ($This->_IsOSO2DivalentOxygen($Atom)) {
      $AtomType = 'OSO2';
      last ATOMTYPE;
    }

    # OSO : DIVALENT OXYGEN ATTACHED TO SULFUR
    if ($This->_IsOSODivalentOxygen($Atom)) {
      $AtomType = 'OSO';
      last ATOMTYPE;
    }

    # OS=O : DIVALENT OXYGEN ATTACHED TO SULFOXIDE SULFUR
    if ($This->_IsSulfoxideDivalentOxygen($Atom)) {
      $AtomType = 'OS=O';
      last ATOMTYPE;
    }

    # -OS : GENERAL DIVALENT OXYGEN ATTACHED TO S
    if ($This->_IsOSDivalentOxygen($Atom)) {
      $AtomType = ' -OS';
      last ATOMTYPE;
    }

    # O4S : TERMINAL O IN SO4(-3)
    if ($This->_IsSulfateTerminalOxygen($Atom)) {
      $AtomType = 'O4S';
      last ATOMTYPE;
    }

    # O3S : TERMINAL O IN SULFONATES
    if ($This->_IsSulfonateTerminalOxygen($Atom)) {
      $AtomType = 'O3S';
      last ATOMTYPE;
    }

    # O2S : TERMINAL O-S IN SULFONES AND SULFONAMIDES
    if ($This->_IsSulfoneOrSulfonamideTerminalOxygen($Atom)) {
      $AtomType = 'O2S';
      last ATOMTYPE;
    }

    # O=S : O=S IN SULFOXIDES
    if ($This->_IsSulfoxideOxygen($Atom)) {
      $AtomType = 'O=S';
      last ATOMTYPE;
    }

    # O=S= :  O=S ON SULFUR DOUBLY BONDED TO, E.G., CARBON
    if ($This->_IsDoublyBondedOSOxygen($Atom)) {
      $AtomType = 'O=S=';
      last ATOMTYPE;
    }

    # O-S : SINGLE TERMINAL OXYGEN ON TETRACOORD SULFUR
    if ($This->_IsOSTerminalOxygen($Atom)) {
      $AtomType = 'O-S';
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForOxygenAttachedToSulfur: MMFF94 atom type for Oxygen cann't be assigned...";
  }

  return $AtomType;
}

# Get MMFF94 atom type for divalent or terminal Oxygen attached to Phosphorus...
#
sub _GetAtomTypeForOxygenAttachedToPhosphorus {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # OPO3 : DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
    if ($This->_IsOPO3DivalentOxygen($Atom)) {
      $AtomType = 'OPO3';
      last ATOMTYPE;
    }

    # OPO2 : DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
    if ($This->_IsOPO2DivalentOxygen($Atom)) {
      $AtomType = 'OPO2';
      last ATOMTYPE;
    }

    # OPO : DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
    if ($This->_IsOPODivalentOxygen($Atom)) {
      $AtomType = 'OPO';
      last ATOMTYPE;
    }

    # -OP : DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
    if ($This->_IsOPDivalentOxygen($Atom)) {
      $AtomType = '-OP';
      last ATOMTYPE;
    }

    # O4P : TERMINAL OXYGEN IN PHOSPHATES AND PHOSPHODIESTERS
    if ($This->_IsPhosphateOrPhosphodiesterTerminalOxygen($Atom)) {
      $AtomType = 'O4P';
      last ATOMTYPE;
    }

    # O3P : TERMINAL OXYGEN IN PHOSPHONATES
    if ($This->_IsPhosphonateTerminalOxygen($Atom)) {
      $AtomType = 'O3P';
      last ATOMTYPE;
    }

    # O2P : TERMINAL O IN PHOSPHINATES
    if ($This->_IsPhosphinateTerminalOxygen($Atom)) {
      $AtomType = 'O2P';
      last ATOMTYPE;
    }

    # OP : TERMINAL O IN PHOSPHOXIDES
    if ($This->_IsPhosphoxideTerminalOxygen($Atom)) {
      $AtomType = 'OP';
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForOxygenAttachedToPhosphorus: MMFF94 atom type for Oxygen cann't be assigned...";
  }

  return $AtomType;
}

# Get MMFF94 atom type for Oxygen with only sigma bonds...
#
sub _GetAtomTypeForOxygenWithOnlySigmaBonds {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # OFUR : AROMATIC OXYGEN AS IN FURAN
    if ($This->_IsAromaticOxygen($Atom)) {
      $AtomType = 'OFUR';
      last ATOMTYPE;
    }

    # OC=O : ESTER OR CARBOXYLIC ACID -O-
    if ($This->_IsEsterOrCarboxylicAcidOxygen($Atom)) {
      $AtomType = 'OC=O';
      last ATOMTYPE;
    }

    # O2CM : OXYGEN IN CARBOXYLATE ANION
    if ($This->_IsCarboxylateAnionOxygen($Atom)) {
      $AtomType = 'O2CM';
      last ATOMTYPE;
    }

    # OC=C : ENOLIC OR PHENOLIC OXYGEN
    if ($This->_IsEnolicOrPhenolicOxygen($Atom)) {
      $AtomType = 'OC=C ';
      last ATOMTYPE;
    }

    # OC=N : DIVALENT OXYGEN
    if ($This->_IsOCNDivalentOxygen($Atom)) {
      $AtomType = 'OC=N';
      last ATOMTYPE;
    }

    # OC=S : THIOESTER OR THIOACID -O-
    if ($This->_IsThioEsterOrThioAcidOxygen($Atom)) {
      $AtomType = 'OC=S';
      last ATOMTYPE;
    }

    # ON=O : DIVALENT NITRITE ETHER OXYGEN
    if ($This->_IsDivalentNitriteEtherOxygen($Atom)) {
      $AtomType = 'ON=O';
      last ATOMTYPE;
    }

    # ONO2 : DIVALENT NITRATE ETHER OXYGEN
    if ($This->_IsDivalentNitrateEtherOxygen($Atom)) {
      $AtomType = 'ONO2';
      last ATOMTYPE;
    }

    # O3N : NITRATE ANION OXYGEN
    if ($This->_IsNitrateAnionOxygen($Atom)) {
      $AtomType = 'O3N';
      last ATOMTYPE;
    }

    # O2N : NITRO OXYGEN
    if ($This->_IsNitroOxygen($Atom)) {
      $AtomType = 'O2N';
      last ATOMTYPE;
    }

    # OXN : N-OXIDE OXYGEN
    if ($This->_IsNOxideOxygen($Atom)) {
      $AtomType = 'OXN';
      last ATOMTYPE;
    }

    # OM2 : OXIDE OXYGEN ON SP2 CARBON, NEGATIVELY CHARGED
    if ($This->_IsNegativelyChargedSP2OxideOxygen($Atom)) {
      $AtomType = 'OM2';
      last ATOMTYPE;
    }

    # OM : ALKOXIDE OXYGEN, NEGATIVELY CHARGED
    if ($This->_IsNegativelyChargedAlkoxideOxygen($Atom)) {
      $AtomType = 'OM';
      last ATOMTYPE;
    }

    # O4CL : OXYGEN IN CLO4(-) ANION - FORMAL CHARGE=-0.25
    if ($This->_IsPerChlorateAnionOxygen($Atom)) {
      $AtomType = 'O4CL';
      last ATOMTYPE;
    }

    # O+ : POSITIVELY CHARGED OXONIUM (TRICOORDINATE) OXYGEN
    if ($This->_IsPositivelyChargedOxoniumOxygen($Atom)) {
      $AtomType = 'O+';
      last ATOMTYPE;
    }

    # OH2 : OXYGEN ON WATER
    if ($This->_IsWaterOxygen($Atom)) {
      $AtomType = 'OH2';
      last ATOMTYPE;
    }

    # OR : ALCOHOL OR ETHER OXYGEN
    if ($This->_IsAlcoholOrEtherOxygen($Atom)) {
      $AtomType = 'OR';
      last ATOMTYPE;
    }

    # GENERAL DIVALENT O
    if ($This->_IsDivalentOxygen($Atom)) {
      $AtomType = '-O-';
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForOxygenWithOnlySigmaBonds: MMFF94 atom type for Oxygen cann't be assigned...";
  }

  return $AtomType;
}

# Get MMFF94 atom type for Oxygen with only sigma bonds...
#
sub _GetAtomTypeForOxygenWithOnePiBond {
  my($This, $Atom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # OFUR : AROMATIC OXYGEN AS IN FURAN
    if ($This->_IsAromaticOxygen($Atom)) {
      $AtomType = 'OFUR';
      last ATOMTYPE;
    }

    # O=+ : POSITIVELY CHARGED OXENIUM (DICOORDINATE) OXYGEN
    if ($This->_IsPositivelyChargedOxeniumOxygen($Atom)) {
      $AtomType = 'O=+';
      last ATOMTYPE;
    }

    # O=CN : CARBONYL OXYGEN, AMIDES
    if ($This->_IsAmideCarbonylOxygen($Atom)) {
      $AtomType = 'O=CN';
      last ATOMTYPE;
    }

    # O=CO : CARBONYL OXYGEN, CARBOXYLIC ACIDS AND ESTERS
    if ($This->_IsCarbobylCarboxylicAcidsOrEstersOxygen($Atom)) {
      $AtomType = 'O=CO';
      last ATOMTYPE;
    }

    # O=CR : CARBONYL OXYGEN, ALDEHYDES AND KETONES
    if ($This->_IsCarbonylAldehydeOrKetoneOxygen($Atom)) {
      $AtomType = 'O=CR';
      last ATOMTYPE;
    }

    # O=C : GENERAL C=O
    if ($This->_IsCarbonylOxygen($Atom)) {
      $AtomType = 'O=C';
      last ATOMTYPE;
    }

    # O2NO : NITRO-GROUP OXYGEN IN NITRATE
    if ($This->_IsNitroGroupNitrateOxygen($Atom)) {
      $AtomType = 'O2NO';
      last ATOMTYPE;
    }

    # O2N : NITRO OXYGEN
    if ($This->_IsNitroOxygen($Atom)) {
      $AtomType = 'O2N';
      last ATOMTYPE;
    }

    # O=N : NITROSO OXYGEN
    if ($This->_IsNitrosoOxygen($Atom)) {
      $AtomType = 'O=N';
      last ATOMTYPE;
    }

    # O4CL : OXYGEN IN CLO4(-) ANION - FORMAL CHARGE=-0.25
    if ($This->_IsPerChlorateAnionOxygen($Atom)) {
      $AtomType = 'O4CL';
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForOxygenWithOnePiBond: MMFF94 atom type for Oxygen cann't be assigned...";
  }

  return $AtomType;
}

# OR : ALCOHOL OR ETHER OXYGEN
#
sub _IsAlcoholOrEtherOxygen {
  my($This, $Atom) = @_;

  return ($This->_IsAlcoholOxygen($Atom) || $This->_IsEtherOxygen($Atom)) ? 1 : 0;
}

# Alcohol Oxygen..
#
sub _IsAlcoholOxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.TSB2', ['C', 'H'], ['-', '-']) ? 1 : 0;
}

# Ether Oxygen..
#
sub _IsEtherOxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.TSB2', ['C', 'C'], ['-', '-']) ? 1 : 0;
}


# OC=O : ESTER OR CARBOXYLIC ACID -O-
#
# Notes:
#    . Carboxylate anion Oxygen is matched using O2CM atom type.
#
sub _IsEsterOrCarboxylicAcidOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.TSB2')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.DB1')) {
    if ($This->_IsCarboxylicAcidOrEsterCarbonylCarbon($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# OC=C : ENOLIC OR PHENOLIC OXYGEN
#
sub _IsEnolicOrPhenolicOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.TSB2')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.DB1')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('C.T3.DB1', ['C', 'O', '*'], ['=', '-', '-'])) {
      return 1;
    }
  }
  return 0;
}

# OC=N : DIVALENT OXYGEN
#
sub _IsOCNDivalentOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.TSB2')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.DB1')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('C.T3.DB1', ['N', 'O', '*'], ['=', '-', '-'])) {
      return 1;
    }
  }
  return 0;
}

# OC=S : THIOESTER OR THIOACID -O-
#
sub _IsThioEsterOrThioAcidOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.TSB2,O.X1.FC-1')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.DB1')) {
    # Thio ester, acid or anion Oxygen...
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('C.T3.DB1', ['S', 'O.X1.FC0,O.X2,O.X1.FC-1', 'C,H'], ['=', '-', '-'], ['C', 'C,H', 'C,H'])) {
      return 1;
    }
  }
  return 0;
}

# ONO2 : DIVALENT NITRATE ETHER OXYGEN
#
# Notes:
#    . Nitrate anion Oxygen is matched using O3N atom type.
#    . Divalent Oxygen in Nitrate with one Hydrogen atom is not treated as ether ONO2.
#
sub _IsDivalentNitrateEtherOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X2')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('N.DB1')) {
    if ($This->_IsNitrateNitrogen($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# ON=O : DIVALENT NITRITE ETHER OXYGEN
#
# Notes:
#    . Divalent Oxygen in Nitrite with one Hydrogen atom is not treated as ether ON=O.
#
sub _IsDivalentNitriteEtherOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X2')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('N.DB1')) {
    # R-O-N=O
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('N.T2.DB1.FC0', ['O', 'O.FC0'], ['=', '-'])) {
      return 1;
    }
  }
  return 0;
}

# OSO3 : DIVALENT OXYGEN ATTACHED TO SULFUR
#
# Notes:
#    . It corresponds to divalent Oxygen in Sulfates.
#    . Oxygen attached to one heavy atom and one Hydrogen atom is treated as divalent
#      Oxygen and matched to OSO3.
#    . Anion Oxygen attached to one heavy atom with no Hydrogen atom is treated as terminal
#      Oxygen and matched to O4S.
#
sub _IsOSO3DivalentOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.TSB2', ['S.T4'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.T4')) {
    # R'-O-S(=O)(=O)-O-R"
    if ($This->_IsSulfateSulfur($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# OSO2 : DIVALENT OXYGEN ATTACHED TO SULFUR
#
# Notes:
#    . It corresponds to divalent Oxygen in Sulfonates, Sulfones and so on.
#    . Oxygen attached to one heavy atom and one Hydrogen atom is treated as divalent
#      Oxygen and matched to OSO2.
#    . Anion Oxygen attached to one heavy atom with no Hydrogen atom is treated as terminal
#      Oxygen and matched to O3S.
#
sub _IsOSO2DivalentOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.TSB2', ['S.T4'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.T4')) {
    # R'-O-S(=O)(=O)-R"
    if ($This->_IsSulfonateSulfur($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# OSO : DIVALENT OXYGEN ATTACHED TO SULFUR
#
# Notes:
#    . It corresponds to divalent Oxygen attached to Sulfur with single bond.
#    . Oxygen attached to one heavy atom and one Hydrogen atom is treated as divalent
#      Oxygen and matched to OSO.
#
sub _IsOSODivalentOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.TSB2', ['S.T2'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.T2.DB0')) {
    # R'-O-S-O-R"
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('S', ['O', 'O'], ['-', '-'])) {
      return 1;
    }
  }
  return 0;
}

# OS=O : DIVALENT OXYGEN ATTACHED TO SULFOXIDE SULFUR
#
# Notes:
#    . It corresponds to divalent Oxygen in Sulfoxides.
#    . Oxygen attached to one heavy atom and one Hydrogen atom is treated as divalent
#      Oxygen and matched to OS=O.
#
sub _IsSulfoxideDivalentOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.TSB2', ['S.T3'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.T3.DB1')) {
    # R'-O-S(=O)-R"
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('S.T3.DB1', ['O', 'O', '*'], ['=', '-', '-'])) {
      return 1;
    }
  }
  return 0;
}

# -OS : GENERAL DIVALENT OXYGEN ATTACHED TO S
#
# Notes:
#    . It covers any divalent Oxygen attached to Sulfur.
#
sub _IsOSDivalentOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.TSB2', ['S'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('S', ['O'], ['-'])) {
      return 1;
    }
  }
  return 0;
}

# OPO3 : DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
#
# Notes:
#    . It corresponds to divalent Oxygen in Phopsphates or Phosphodiesters.
#
sub _IsOPO3DivalentOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.TSB2')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('P.T4')) {
    if ($This->_IsPhosphateOrPhosphodiesterPhosphorus($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# OPO2 : DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
#
# Notes:
#    . It corresponds to divalent Oxygen in Phopsphonates or their esters.
#
sub _IsOPO2DivalentOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.TSB2')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('P.T4')) {
    if ($This->_IsPhosphonatePhosphorus($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# OPO : DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
#
# Notes:
#    . It corresponds to divalent Oxygen in Phopsphinates.
#
sub _IsOPODivalentOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.TSB2')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('P.T4')) {
    if ($This->_IsPhosphinatePhosphorus($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# -OP : DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
#
sub _IsOPDivalentOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.TSB2')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('P')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('P', ['O'], ['-'])) {
      return 1;
    }
  }
  return 0;
}

# -O- : GENERAL DIVALENT O
#
sub _IsDivalentOxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.TSB2') ? 1 : 0;
}

# O=C : GENERAL C=O
#
sub _IsCarbonylOxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.T1.DB1', ['C'], ['=']) ? 1 : 0;
}

# O=CN : CARBONYL OXYGEN, AMIDES
#
sub _IsAmideCarbonylOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.T1.DB1')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.DB1')) {
    if ($This->_IsAmideCarbonylCarbon($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# O=CR : CARBONYL OXYGEN, ALDEHYDES AND KETONES
#
sub _IsCarbonylAldehydeOrKetoneOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.T1.DB1')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.DB1')) {
    if ($This->_IsKetoneOrAldehydeCarbonylCarbon($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# O=CO : CARBONYL OXYGEN, CARBOXYLIC ACIDS AND ESTERS
#
sub _IsCarbobylCarboxylicAcidsOrEstersOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.T1.DB1')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.DB1')) {
    if ($This->_IsCarboxylicAcidOrEsterCarbonylCarbon($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# O=N : NITROSO OXYGEN
#
sub _IsNitrosoOxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.T1.DB1', ['N'], ['=']) ? 1 : 0;
}

# O=S : O=S IN SULFOXIDES
#
sub _IsSulfoxideOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.T1.DB1')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.DB1')) {
    if ($This->_IsSulfoxideSulfur($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# O=S= :  O=S ON SULFUR DOUBLY BONDED TO, E.G., CARBON
#
sub _IsDoublyBondedOSOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.T1.DB1')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.DB2')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('S', ['O', '!O'], ['=', '='])) {
      return 1;
    }
  }
  return 0;
}

# O2CM : OXYGEN IN CARBOXYLATE ANION
#
sub _IsCarboxylateAnionOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X1.FC-1')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.T3.DB1')) {
    if ($This->_IsCarboxylateAnionCarbon($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# OXN : N-OXIDE OXYGEN
#
sub _IsNOxideOxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.X1.FC-1', ['N.FC+1'], ['-']) ? 1 : 0;
}

# O2N : NITRO OXYGEN
#
sub _IsNitroOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X1.DB1,O.X1.FC-1')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('N.T3.DB1')) {
    if ($This->_IsNitroNitrogen($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# O2NO : NITRO-GROUP OXYGEN IN NITRATE
#
sub _IsNitroGroupNitrateOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X1.DB1')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('N.T3.DB1')) {
    if ($This->_IsNitrateNitrogen($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# O3N : NITRATE ANION OXYGEN
#
sub _IsNitrateAnionOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X1.FC-1')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('N.T3.DB1')) {
    if ($This->_IsNitrateNitrogen($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# O-S : SINGLE TERMINAL OXYGEN ON TETRACOORD SULFUR
#
sub _IsOSTerminalOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X1', ['S.T4'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.T4')) {
    if ($AtomNeighbor->DoesAtomNeighborhoodMatch('S', ['O.X1', '!O', '!O', '!O'])) {
      return 1;
    }
  }
  return 0;
}

# O2S : TERMINAL O-S IN SULFONES AND SULFONAMIDES
#
sub _IsSulfoneOrSulfonamideTerminalOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X1', ['S.T4'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.T4')) {
    if ($This->_IsSulfonamideSulfur($AtomNeighbor) || $This->_IsSulfoneSulfur($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# O3S : TERMINAL O IN SULFONATES
#
# Notes:
#    . It corresponds to monovalent Oxygen in Sulfonates.
#    . Anion Oxygen attached to one heavy atom with no Hydrogen atom is treated as terminal
#      Oxygen and matched to O3S.
#
sub _IsSulfonateTerminalOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X1', ['S.T4'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.T4')) {
    # R'-O-S(=O)(=O)-R", [O-]-S(=O)(=O)-R",
    if ($This->_IsSulfonateSulfur($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# O4S : TERMINAL O IN SO4(-3)
#
# Notes:
#    . It corresponds to monovalent Oxygen in Sulfates.
#    . As far I can tell, SO4 should have a formal charge of -2.
#    . Anion Oxygen attached to one heavy atom with no Hydrogen atom is treated as terminal
#      Oxygen and matched to O4S.
#
sub _IsSulfateTerminalOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X1', ['S.T4'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.T4')) {
    # R'-O-S(=O)(=O)-O-R", R'-O-S(=O)(=O)-[O-]
    if ($This->_IsSulfateSulfur($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# OSMS : TERM O IN THIOSULFINATE ANION - FORMAL CHARGE=-0.5
#
sub _IsThioSulfinateTerminalOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X1', ['S.FC+1'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.FC+1')) {
    if ($This->_IsThioSulfinateSulfur($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# OP : TERMINAL O IN PHOSPHOXIDES
#
sub _IsPhosphoxideTerminalOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X1', ['P.T4'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('P.T4')) {
    if ($This->_IsPhosphoxidePhosphorus($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# O2P : TERMINAL O IN PHOSPHINATES
#
sub _IsPhosphinateTerminalOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X1', ['P.T4'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('P.T4')) {
    if ($This->_IsPhosphinatePhosphorus($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# O3P : TERMINAL OXYGEN IN PHOSPHONATES
#
sub _IsPhosphonateTerminalOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X1', ['P.T4'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('P.T4')) {
    if ($This->_IsPhosphonatePhosphorus($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# O4P : TERMINAL OXYGEN IN PHOSPHATES AND PHOSPHODIESTERS
#
sub _IsPhosphateOrPhosphodiesterTerminalOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('O.X1', ['P.T4'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('P.T4')) {
    if ($This->_IsPhosphateOrPhosphodiesterPhosphorus($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# O4CL : OXYGEN IN CLO4(-) ANION - FORMAL CHARGE=-0.25
#
sub _IsPerChlorateAnionOxygen {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  # All Oxygens in PerChlorate anion and ester are matched to O4Cl...
  if (!$Atom->DoesAtomNeighborhoodMatch('O', ['Cl'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('Cl')) {
    if ($This->_IsPerChlorateAnionChlorine($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# OM : ALKOXIDE OXYGEN, NEGATIVELY CHARGED
#
sub _IsNegativelyChargedAlkoxideOxygen {
  my($This, $Atom) = @_;

  # R-(O-)

  return $Atom->DoesAtomNeighborhoodMatch('O.FC-1.T1', ['C,H'], ['-']) ? 1 : 0;
}

# OM2 : OXIDE OXYGEN ON SP2 CARBON, NEGATIVELY CHARGED
#
sub _IsNegativelyChargedSP2OxideOxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.FC-1.T1', ['C.DB1.T3'], ['-']) ? 1 : 0;
}

# O+ : POSITIVELY CHARGED OXONIUM (TRICOORDINATE) OXYGEN
#
sub _IsPositivelyChargedOxoniumOxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.FC+1.T3.DB0') ? 1 : 0;
}

# O=+ : POSITIVELY CHARGED OXENIUM (DICOORDINATE) OXYGEN
#
sub _IsPositivelyChargedOxeniumOxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.FC+1.T2.DB1') ? 1 : 0;
}

# OFUR : AROMATIC OXYGEN AS IN FURAN
#
sub _IsAromaticOxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.Ar.RA5') ? 1 : 0;
}

# OH2 : OXYGEN ON WATER
#
sub _IsWaterOxygen {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('O.T2.TSB2', ['H', 'H'], ['-', '-']) ? 1 : 0;
}

# PO4 : PHOSPHOROUS IN PHOSPHATES AND PHOSPHODIESTERS
#
sub _IsPhosphateOrPhosphodiesterPhosphorus {
  my($This, $Atom) = @_;

  # Phosphate: R-O-P(=O)(-O)-O; Phosphate diester: R'-O-P(=O)(-O)-O-R"

  return $Atom->DoesAtomNeighborhoodMatch('P.T4.DB1', ['O', 'O', 'O', 'O'], ['=', '-', '-', '-']) ? 1 : 0;
}

# PO3 : TETRACOORDINATE P WITH THREE ATTACHED OXYGENS
#
sub _IsPhosphonatePhosphorus {
  my($This, $Atom) = @_;

  # Phosphonate: R-P(=O)(-O)-O; Phosphonate ester: R'-P(=O)(-O)-O-R"

  return $Atom->DoesAtomNeighborhoodMatch('P.T4.DB1', ['O', 'O', 'O', '!O'], ['=', '-', '-', '-']) ? 1 : 0;
}

# PO2 : TETRACOORDINATE P WITH TWO ATTACHED OXYGENS
#
sub _IsPhosphinatePhosphorus {
  my($This, $Atom) = @_;

  # Phosphinate: R-P(=O)(-O)-R"; Phosphinate ester: R'-P(=O)(-O)-R"

  return $Atom->DoesAtomNeighborhoodMatch('P.T4.DB1', ['O', 'O', '!O', '!O'], ['=', '-', '-', '-']) ? 1 : 0;
}

# PO : TETRACOORDINATE P WITH ONE ATTACHED OXYGEN
#
sub _IsPhosphoxidePhosphorus {
  my($This, $Atom) = @_;

  # Phosphoxide: R-P(=O)(-R")-R'''

  return $Atom->DoesAtomNeighborhoodMatch('P.T4.DB1', ['O', '!O', '!O', '!O'], ['=', '-', '-', '-']) ? 1 : 0;
}

# PTET : GENERAL TETRACOORDINATE PHOSPHORUS
#
sub _IsTetraCoordinatedPhosphorus {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('P.T4') ? 1 : 0;
}

# P : TRICOORDINATE P, AS IN PHOSPHINES
#
sub _IsTriCoordinatedPhosphorus {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('P.T3') ? 1 : 0;
}

# -P=C : PHOSPHOROUS DOUBLY BONDED TO CARBON
#
sub _IsDoublyBondedToCarbonPhosphorous {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('P.DB1', ['C'], ['=']) ? 1 : 0;
}

# S : SULFUR IN THIOETHERS AND MERCAPTANS
#
sub _IsThioEthersOrMercaptansSulfur {
  my($This, $Atom) = @_;

  return ($This->_IsThioEtherSulfur($Atom) || $This->_IsMercaptansSulfur($Atom)) ? 1 : 0;

  return 0;
}

# Thioethers: R'-S-R"
#
sub _IsThioEtherSulfur {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  return $Atom->DoesAtomNeighborhoodMatch('S.X2.T2', ['C', 'C'], ['-', '-']) ? 1 : 0;
}

# Mercaptans or Thiols: R-S-H
#
sub _IsMercaptansSulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.T2', ['C', 'C,H'], ['-', '-']) ? 1 : 0;
}

# Is it a divalent dicoordinated Sulfur...
#
sub _IsDivalentDiCoordinatedSulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.T2.FC0', ['*', '*'], ['-', '-']) ? 1 : 0;
}

# S=C : TERMINAL SULFUR DOUBLY BONDED TO CARBON
#
sub _IsSCTerminalSulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.X1.DB1', ['C'], ['=']) ? 1 : 0;
}

# >S=N : SULFUR, TRICOORD, DOUBLY BONDED TO N
#
sub _IsSNTricoordinatedSulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.T3.DB1', ['N', '*', '*'], ['=', '-', '-']) ? 1 : 0;
}

# S=O : SULFUR IN SULFOXIDES
#
sub _IsSulfoxideSulfur {
  my($This, $Atom) = @_;

  # Sulfone: R'-S(=O)-R", R'-S(=O)-O-R", and so on...

  return $Atom->DoesAtomNeighborhoodMatch('S.T3.DB1', ['O', '*', '*'], ['=', '-', '-']) ? 1 : 0;
}

# SO2 : SULFUR IN SULFONES
#
sub _IsSulfoneSulfur {
  my($This, $Atom) = @_;

  # Sulfone: R'-S(=O)(=O)-R"

  return $Atom->DoesAtomNeighborhoodMatch('S.T4.DB2', ['O', 'O', '!O', '!O'], ['=', '=', '-', '-']) ? 1 : 0;
}

#  =SO2: SULFONE SULPHER DOUBLY BONDED TO CARBON
#
sub _IsDoublyBondedToCarbonSulfoneSulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.T3.DB3', ['C', 'O', 'O'], ['=', '=', '=']) ? 1 : 0;
}

# SO2N : SULFUR IN SULFONAMIDES
#
sub _IsSulfonamideSulfur {
  my($This, $Atom) = @_;

  # Sulfonamide: R-S(=O)(=O)-N(-R)(-R")

  return $Atom->DoesAtomNeighborhoodMatch('S.T4.DB2', ['O', 'O', 'N', '!O'], ['=', '=', '-', '-']) ? 1 : 0;
}

# SO3 : SULFONATE SULFUR
#
sub _IsSulfonateSulfur {
  my($This, $Atom) = @_;

  # Sulfonate ion: R'-S(=O)(=O)-(O-); Sulfonate ester: R'-S(=O)(=O)-O-R"

  return $Atom->DoesAtomNeighborhoodMatch('S.T4.DB2', ['O', 'O', 'O', '!O'], ['=', '=', '-', '-']) ? 1 : 0;
}

# SO4 : SULFATE SULFUR
#
sub _IsSulfateSulfur {
  my($This, $Atom) = @_;

  # Sulfate ion: (O-)-S(=O)(=O)-(O-); Sulfate esters: R-O-S(=O)(=O)-O-R"

  return $Atom->DoesAtomNeighborhoodMatch('S.T4.DB2', ['O', 'O', 'O', 'O'], ['=', '=', '-', '-']) ? 1 : 0;
}

# SNO : SULFUR IN NITROGEN ANALOG OF A SULFONE
#
sub _IsNitrogenAnalogOfSulfoneSulfur {
  my($This, $Atom) = @_;

  # Sulfone: R'-S(=N)(=O)-R"

  return $Atom->DoesAtomNeighborhoodMatch('S.T4.DB2', ['N', 'O', '!O', '!O'], ['=', '=', '-', '-']) ? 1 : 0;
}

# STHI : SULFUR AS IN THIOPHENE
#
sub _IsSTHISulfur {
  my($This, $Atom) = @_;

  # Is it an aromatic atom in a five membered ring?
  if (!$Atom->DoesAtomNeighborhoodMatch('S.Ar.RA5.FC0')) {
    return 0;
  }

  # Is it part of five membered ring containing only one hetero atom?
  my($RingAtomsRef, $RingIsAromatic, $NumOfHeteroAtoms);

  for $RingAtomsRef ($Atom->GetRingsWithSize(5)) {
    ($RingIsAromatic, $NumOfHeteroAtoms) = $This->_GetHeteroAtomsInformationInRing($RingAtomsRef);
    if ($RingIsAromatic && $NumOfHeteroAtoms == 1) {
      return 1;
    }
  }
  return 0;
}

# S-P : TERMINAL SULFUR BONDED TO PHOSPHORUS
#
sub _IsSPTerminalSulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.X1.T2', ['P', '*'], ['-', '-']) ? 1 : 0;
}

# S2CM : TERMINAL SULFUR IN THIOCARBOXYLATE ANION
#
sub _IsThioCarboxylateAnionTerminalSulfur {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('S.X1.FC-1, S.X1.DB1.FC0')) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('C.DB1')) {
    if ($This->_IsThioCarboxylateAnionCarbon($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# SM : TERMINAL SULFUR - FORMAL CHARGE=-1
#
sub _IsNegativelyChargedTerminalSulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.X1.FC-1', ['*'], ['-']) ? 1 : 0;
}

# SSMO : TERMINAL SULFUR IN THIOSULFINATE GROUP
#
sub _IsThioSulfinateTerminalSulfur {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('S.X1', ['S.FC+1'])) {
    return 0;
  }

  for $AtomNeighbor ($Atom->GetNeighborsUsingAtomSpecification('S.FC+1')) {
    if ($This->_IsThioSulfinateSulfur($AtomNeighbor)) {
      return 1;
    }
  }
  return 0;
}

# SO2M : SULFUR IN NEGATIVELY CHARGED SULFINATE GROUP
#
sub _IsNegativelyChargedSulfinateSulfur {
  my($This, $Atom) = @_;

  # Sulfinate ion: R'-S(=O)-(O-)

  return $Atom->DoesAtomNeighborhoodMatch('S.T3.DB1', ['O', 'O.X1.FC-1', '!O'], ['=', '-', '-']) ? 1 : 0;
}

# Sulfinate: R'-S(=O)-OH; Sulfinate esters: R'-S(=O)-O-R
#
sub _IsSulfinateSulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.T3.DB1', ['O', 'O.X1.FC0,O.X2', '!O'], ['=', '-', '-']) ? 1 : 0;
}

# SSOM : TRICOORD SULFUR IN THIOSULFINATE GROUP
#
sub _IsTriCoordinatedThioSulfinateSulfur {
  my($This, $Atom) = @_;
  my($AtomNeighbor);

  if (!$Atom->DoesAtomNeighborhoodMatch('S.X3', ['S'])) {
    return 0;
  }

  if ($This->_IsThioSulfinateSulfur($Atom)) {
      return 1;
  }

  return 0;
}

# Is it a Thiosulfinate group?
#
sub _IsThioSulfinateSulfur {
  my($This, $Atom) = @_;

  # R'-[S+1](-[O-1])-S-R"

  return $Atom->DoesAtomNeighborhoodMatch('S.FC+1', ['O.FC-1', 'S', '!O'], ['-', '-', '-']) ? 1 : 0;
}

#   =S=O:  SULFINYL SULFUR, EG. IN C=S=O
#
sub _IsSulfinylSulfur {
  my($This, $Atom) = @_;

  return $Atom->DoesAtomNeighborhoodMatch('S.X2.DB2', ['C', 'O'], ['=', '=']) ? 1 : 0;
}

# CLO4 : CHLORINE IN PERCHLORATE ANION, CLO4(-)
#
sub _IsPerChlorateAnionChlorine {
  my($This, $Atom) = @_;

  # R-O-Cl(=O)(=O)=O or (O-)-Cl(=O)(=O)=O
  if ($Atom->DoesAtomNeighborhoodMatch('Cl.X4.DB3', ['O.X1.FC-1,O.T2.FC0', 'O', 'O', 'O'], ['-', '=', '=', '='])) {
    return 1;
  }

  # Match distributed formal charge of -0.25 on each Oxygen?
  if ($Atom->DoesAtomNeighborhoodMatch('Cl.X4.DB3', ['O.FC-0.25', 'O.FC-0.25', 'O.FC-0.25', 'O.FC-0.25'], ['*', '*', '*', '*'])) {
    return 1;
  }

  return 0;
}

# Get MMFF94 atom type for Hydrogen attached to Carbon...
#
sub _GetAtomTypeForHydrogenAttachedToCarbon {
  my($This, $CarbonAtom) = @_;
  my($AtomType);

  # HC :  H  ATTACHED TO C
  $AtomType = 'HC';

  return $AtomType;
}

# Get MMFF94 atom type for Hydrogen attached to Nitrogen...
#
# 25 AtomTypeSymbols for element H:
#
# AtomTypeSymbol   AtomTypeNum   AtomTypeDefinition
#   HNR      23    H-N(SP3)
#   H3N      23    H-N(SP3), AMMONIA
#   HPYL     23    H-N IN PYRROLE
#   HNOX     23    H-N IN IN A N-OXIDE
#   HNM      23    H ON DICOORD, NEGATIVELY CHARGED NITROGEN
#   HN       23    GENERAL H ON NITROGEN
#   HN=N     27    AZO HYDROGEN
#   HN=C     27    IMINE HYDROGEN
#   HNCO     28    AMIDE HYDROGEN
#   HNCS     28    THIOAMIDE HYDROGEN
#   HNCC     28    H-N IN ENAMINES
#   HNCN     28    H-N IN H-N-C=N
#   HNNC     28    H-N IN H-N-N=C
#   HNNN     28    H-N IN H-N-N=N
#   HNSO     28    H-N IN SULFONAMIDE
#   HNPO     28    H-N IN PHOSPHONAMIDE
#   HNC%     28    HYDROGEN ON N ATTACHED TO TRIPLY BONDED CARBON
#   HSP2     28    GENERAL H ON SP2 NITROGEN
#   HNR+     36    H ON QUATERNARY NITROGEN
#   HIM+     36    H ON IMIDAZOLIUM-TYPE NITROGEN
#   HPD+     36    H ON PROTONATED PYRIDINE NITROGEN
#   HNN+     36    H ON AMIDINIUM-TYPE NITROGEN
#   HNC+     36    H ON PROTONATED IMINE NITROGEN
#   HGD+     36    H ON GUANIDINIUM-TYPE NITROGEN
#   HN5+     36    H ON N5+, N5A+ OR N5B+
#
sub _GetAtomTypeForHydrogenAttachedToNitrogen {
  my($This, $NitrogenAtom) = @_;
  my($AtomType, $NumOfSigmaBonds, $NumOfPiBonds);

  $AtomType = 'None';

  ($NumOfSigmaBonds, $NumOfPiBonds) = ('0') x 2;

  ($NumOfSigmaBonds, $NumOfPiBonds) = $NitrogenAtom->GetNumOfSigmaAndPiBondsToNonHydrogenAtoms();
  $NumOfSigmaBonds += $NitrogenAtom->GetAtomicInvariantValue('H');

  ATOMTYPE: {

    if ($NumOfPiBonds == 0) {
      $AtomType = $This->_GetAtomTypeForHydrogenAttachedToNitrogenWithOnlySigmaBonds($NitrogenAtom);
      last ATOMTYPE;
    }

    if ($NumOfPiBonds == 1) {
      $AtomType = $This->_GetAtomTypeForHydrogenAttachedToNitrogenWithOnePiBond($NitrogenAtom);
      last ATOMTYPE;
    }

    if ($NumOfPiBonds == 2) {
      $AtomType = $This->_GetAtomTypeForHydrogenAttachedToNitrogenWithTwoPiBonds($NitrogenAtom);
      last ATOMTYPE;
    }

    # HN : GENERAL H ON NITROGEN
    $AtomType = 'HN';
  }
  return $AtomType;
}

# Get atom type for Hydrogen attached to a Nitrogen with only sigma bonds...
#
sub _GetAtomTypeForHydrogenAttachedToNitrogenWithOnlySigmaBonds {
  my($This, $NitrogenAtom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # HNC% : HYDROGEN ON N ATTACHED TO TRIPLY BONDED CARBON
    if ($This->_IsHydrogenAttachedToTriplyBondedToCarbonNitrogen($NitrogenAtom)) {
      $AtomType = 'HNC%';
      last ATOMTYPE;
    }

    # HPYL : H-N IN PYRROLE
    if ($This->_IsHydrogenAttachedToPyrroleNitrogen($NitrogenAtom)) {
      $AtomType = 'HPYL';
      last ATOMTYPE;
    }

    # HN5+ : H ON N5+, N5A+ OR N5B+
    if ($This->_IsHydrogenAttachedToFiveMemberedHetreoCyclicPostivelyChargedNitrogen($NitrogenAtom)) {
      $AtomType =  'HN5+';
      last ATOMTYPE;
    }

    # HGD+ : H ON GUANIDINIUM-TYPE NITROGEN
    if ($This->_IsHydrogenAttachedToGuanidiniumNitrogen($NitrogenAtom)) {
      $AtomType = 'HGD+';
      last ATOMTYPE;
    }

    # HNCS : THIOAMIDE HYDROGEN
    if ($This->_IsHydrogenAttachedToThioamideNitrogen($NitrogenAtom)) {
      $AtomType = 'HNCS';
      last ATOMTYPE;
    }

    # HNCO : AMIDE HYDROGEN
    if ($This->_IsHydrogenAttachedToAmideNitrogen($NitrogenAtom)) {
      $AtomType = 'HNCO';
      last ATOMTYPE;
    }

    # HNSO : H-N IN SULFONAMIDE
    if ($This->_IsHydrogenAttachedToSulfonamideNitrogen($NitrogenAtom)) {
      $AtomType = 'HNSO';
      last ATOMTYPE;
    }

    # HNPO : H-N IN PHOSPHONAMIDE
    if ($This->_IsHydrogenAttachedToPhosphonamideNitrogen($NitrogenAtom)) {
      $AtomType = 'HNPO';
      last ATOMTYPE;
    }

    # HNOX : H-N IN IN A N-OXIDE
    if ($This->_IsHydrogenAttachedToNOXideNitrogen($NitrogenAtom)) {
      $AtomType = 'HNOX';
      last ATOMTYPE;
    }

    # HNCC : H-N IN ENAMINES (H-N-C=C)
    if ($This->_IsHydrogenAttachedToEnamineNitrogen($NitrogenAtom)) {
      $AtomType = 'HNCC';
      last ATOMTYPE;
    }

    # HNCN : H-N IN H-N-C=N
    if ($This->_IsHydrogenAttachedToNCNNitrogen($NitrogenAtom)) {
      $AtomType = 'HNCN';
      last ATOMTYPE;
    }

    #  HNNC : H-N IN H-N-N=C
    if ($This->_IsHydrogenAttachedToNNCNitrogen($NitrogenAtom)) {
      $AtomType = 'HNNC';
      last ATOMTYPE;
    }

    # HNNN : H-N IN H-N-N=N
    if ($This->_IsHydrogenAttachedToNNNNitrogen($NitrogenAtom)) {
      $AtomType = 'HNNN';
      last ATOMTYPE;
    }

    # HNM : H ON DICOORD, NEGATIVELY CHARGED NITROGEN
    if ($This->_IsHydrogenAttachedToNegativelyChargedDicoordinatedNitrogen($NitrogenAtom)) {
      $AtomType = 'HNM';
      last ATOMTYPE;
    }

    # H3N : H-N(SP3), AMMONIA
    if ($This->_IsHydrogenAttachedToSP3AmmoniaNitrogen($NitrogenAtom)) {
      $AtomType = 'H3N';
      last ATOMTYPE;
    }

    # HNR : H-N(SP3)
    if ($This->_IsHydrogenAttachedToSP3Nitrogen($NitrogenAtom)) {
      $AtomType = 'HNR';
      last ATOMTYPE;
    }

    # HNR+ : H ON QUATERNARY NITROGEN
    if ($This->_IsHydrogenAttachedToQuaternaryNitrogen($NitrogenAtom)) {
      $AtomType = 'HNR+';
      last ATOMTYPE;
    }

    # HN : GENERAL H ON NITROGEN
    $AtomType = 'HN';
  }
  return $AtomType;
}

# Get atom type for Hydrogen attached to a Nitrogen with one pi bonds...
#
sub _GetAtomTypeForHydrogenAttachedToNitrogenWithOnePiBond {
  my($This, $NitrogenAtom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # HPYL : H-N IN PYRROLE
    if ($This->_IsHydrogenAttachedToPyrroleNitrogen($NitrogenAtom)) {
      $AtomType = 'HPYL';
      last ATOMTYPE;
    }

    # HIM+ : H ON IMIDAZOLIUM-TYPE NITROGEN
    if ($This->_IsHydrogenAttachedToImidazoliumNitrogen($NitrogenAtom)) {
      $AtomType =  'HIM+';
      last ATOMTYPE;
    }

    # HN5+ : H ON N5+, N5A+ OR N5B+
    if ($This->_IsHydrogenAttachedToFiveMemberedHetreoCyclicPostivelyChargedNitrogen($NitrogenAtom)) {
      $AtomType =  'HN5+';
      last ATOMTYPE;
    }

    # HPD+ : H ON PROTONATED PYRIDINE NITROGEN
    if ($This->_IsHydrogenAttachedToPositivelyChargedPyridineNitrogen($NitrogenAtom)) {
      $AtomType = 'HPD+';
      last ATOMTYPE;
    }

    # HNOX : H-N IN IN A N-OXIDE
    if ($This->_IsHydrogenAttachedToNOXideNitrogen($NitrogenAtom)) {
      $AtomType = 'HNOX';
      last ATOMTYPE;
    }

    # HGD+ : H ON GUANIDINIUM-TYPE NITROGEN
    if ($This->_IsHydrogenAttachedToGuanidiniumNitrogen($NitrogenAtom)) {
      $AtomType = 'HGD+';
      last ATOMTYPE;
    }

    # HNN+ : H ON AMIDINIUM-TYPE NITROGEN
    if ($This->_IsHydrogenAttachedToAmidiniumNitrogen($NitrogenAtom)) {
      $AtomType = 'HNN+';
      last ATOMTYPE;
    }

    # HNC+ : H ON PROTONATED IMINE NITROGEN
    if ($This->_IsHydrogenAttachedToPositivelyChargedImineNitrogen($NitrogenAtom)) {
      $AtomType = 'HNC+';
      last ATOMTYPE;
    }

    # HN=N : AZO HYDROGEN
    if ($This->_IsHydrogenAttachedToAzoNitrogen($NitrogenAtom)) {
      $AtomType = 'HN=N';
      last ATOMTYPE;
    }

    # HN=C : IMINE HYDROGEN
    if ($This->_IsHydrogenAttachedToImineNitrogen($NitrogenAtom)) {
      $AtomType = 'HN=C';
      last ATOMTYPE;
    }

    # HSP2: GENERAL H ON SP2 NITROGEN
    $AtomType = 'HSP2';
  }
  return $AtomType;
}

# Get atom type for Hydrogen attached to a Nitrogen with two pi bonds...
#
sub _GetAtomTypeForHydrogenAttachedToNitrogenWithTwoPiBonds {
  my($This, $NitrogenAtom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # HN : GENERAL H ON NITROGEN
    $AtomType = 'HN';
  }
  return $AtomType;
}

# HNR : H-N(SP3)
#
sub _IsHydrogenAttachedToSP3Nitrogen {
  my($This, $NitrogenAtom) = @_;

  return $NitrogenAtom->DoesAtomNeighborhoodMatch('N.T3.FC0', ['*', '*', '*'], ['-', '-', '-']) ? 1 : 0;
}

# H3N : H-N(SP3), AMMONIA
#
sub _IsHydrogenAttachedToSP3AmmoniaNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $NitrogenAtom->DoesAtomNeighborhoodMatch('N.T4.FC+1', ['H', 'H', 'H', 'H'], ['-', '-', '-', '-']) ? 1 : 0;
}

# HPYL : H-N IN PYRROLE
#
sub _IsHydrogenAttachedToPyrroleNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsPyrroleNitrogen($NitrogenAtom) ? 1 : 0;
}

# HNOX : H-N IN IN A N-OXIDE
#
sub _IsHydrogenAttachedToNOXideNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $NitrogenAtom->DoesAtomNeighborhoodMatch('N.FC+1', ['O.X1.FC-1'], ['-']) ? 1 : 0;
}

# HNM : H ON DICOORD, NEGATIVELY CHARGED NITROGEN
#
sub _IsHydrogenAttachedToNegativelyChargedDicoordinatedNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $NitrogenAtom->DoesAtomNeighborhoodMatch('N.T2.FC-1', ['*', '*'], ['-', '-']) ? 1 : 0;
}

# HN=N : AZO HYDROGEN
#
sub _IsHydrogenAttachedToAzoNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsAzoNitrogen($NitrogenAtom) ? 1 : 0;
}

# HN=C : IMINE HYDROGEN
#
sub _IsHydrogenAttachedToImineNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsImineNitrogen($NitrogenAtom) ? 1 : 0;
}

# HNCO : AMIDE HYDROGEN
#
sub _IsHydrogenAttachedToAmideNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsAmideNitrogen($NitrogenAtom) ? 1 : 0;
}

# HNCS : THIOAMIDE HYDROGEN
#
sub _IsHydrogenAttachedToThioamideNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsThioAmideNitrogen($NitrogenAtom) ? 1 : 0;
}

# HNCC : H-N IN ENAMINES (H-N-C=C)
#
sub _IsHydrogenAttachedToEnamineNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsNCCNitrogen($NitrogenAtom) ? 1 : 0;
}

# HNCN : H-N IN H-N-C=N
#
sub _IsHydrogenAttachedToNCNNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsNCNNitrogen($NitrogenAtom) ? 1 : 0;
}

#  HNNC : H-N IN H-N-N=C
#
sub _IsHydrogenAttachedToNNCNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsNNCNitrogen($NitrogenAtom) ? 1 : 0;
}

# HNNN : H-N IN H-N-N=N
#
sub _IsHydrogenAttachedToNNNNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsNNNNitrogen($NitrogenAtom) ? 1 : 0;
}

# HNSO : H-N IN SULFONAMIDE
#
sub _IsHydrogenAttachedToSulfonamideNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsNSO2SulfonamideNitrogen($NitrogenAtom) ? 1 : 0;
}

# HNPO : H-N IN PHOSPHONAMIDE
#
sub _IsHydrogenAttachedToPhosphonamideNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsNPO2PhosphonamideNitrogen($NitrogenAtom) ? 1 : 0;
}

# HNC% : HYDROGEN ON N ATTACHED TO TRIPLY BONDED CARBON
#
sub _IsHydrogenAttachedToTriplyBondedToCarbonNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsAttchedToCCTripleBondNitrogen($NitrogenAtom) ? 1 : 0;
}

# HSP2 : GENERAL H ON SP2 NITROGEN
#
sub _IsHydrogenAttachedToSP2Nitrogen {
  my($This, $NitrogenAtom) = @_;

  return $NitrogenAtom->DoesAtomNeighborhoodMatch('N.DB1') ? 1 : 0;
}

# HNR+ : H ON QUATERNARY NITROGEN
#
sub _IsHydrogenAttachedToQuaternaryNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $NitrogenAtom->DoesAtomNeighborhoodMatch('N.T4.FC+1', ['*', '*', '*', '*'], ['-', '-', '-', '-']) ? 1 : 0;
}

# HIM+ : H ON IMIDAZOLIUM-TYPE NITROGEN
#
sub _IsHydrogenAttachedToImidazoliumNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsImidazoliumNitrogen($NitrogenAtom) ? 1 : 0;
}

# HPD+ : H ON PROTONATED PYRIDINE NITROGEN
#
sub _IsHydrogenAttachedToPositivelyChargedPyridineNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsPyridiniumNitrogen($NitrogenAtom) ? 1 : 0;
}

# HNN+ : H ON AMIDINIUM-TYPE NITROGEN
#
sub _IsHydrogenAttachedToAmidiniumNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsPositivelyChargedAzoNitrogen($NitrogenAtom) ? 1 : 0;
}

# HNC+ : H ON PROTONATED IMINE NITROGEN
#
sub _IsHydrogenAttachedToPositivelyChargedImineNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsPositivelyChargedIminiumNitrogen($NitrogenAtom) ? 1 : 0;
}

# HGD+ : H ON GUANIDINIUM-TYPE NITROGEN
#
sub _IsHydrogenAttachedToGuanidiniumNitrogen {
  my($This, $NitrogenAtom) = @_;

  return $This->_IsGuanidiniumNitrogen($NitrogenAtom) ? 1 : 0;
}

# HN5+ : H ON N5+, N5A+ OR N5B+
#
sub _IsHydrogenAttachedToFiveMemberedHetreoCyclicPostivelyChargedNitrogen {
  my($This, $NitrogenAtom) = @_;

  if (!$NitrogenAtom->DoesAtomNeighborhoodMatch('N.RA5.FC+1')) {
    return 0;
  }

  return ($This->_IsPositivelyChargedFiveMemberedHeteroAromaticRingAlphaNitrogen($NitrogenAtom) ||
	 $This->_IsPositivelyChargedFiveMemberedHeteroAromaticRingBetaNitrogen($NitrogenAtom) ||
	 $This->_IsPositivelyChargedFiveMemberedHeteroCyclicRingNitrogen($NitrogenAtom)) ? 1 : 0;
}

# Get MMFF94 atom type for Hydrogen attached to Oxygen...
#
# 11 AtomTypeSymbols for element H:
#
# AtomTypeSymbol   AtomTypeNum   AtomTypeDefinition
#   HOR      21    HYDROGEN IN ALCOHOLS
#   HO       21    GENERAL H ON OXYGEN
#   HOM      21    HYDROGEN IN HYDROXIDE ANION
#   HOCO     24    H-O IN CARBOXYLIC ACIDS
#   HOP      24    HYDROGEN ON OXYGEN ATTACHED TO PHOSPHOROUS
#   HOCC     29    H-O IN ENOLS AND PHENOLS
#   HOCN     29    H-O IN HO-C=N
#   HOH      31    HYDROGEN IN H2O
#   HOS      33    H ON OXYGEN ATTACHED TO SULFUR
#   HO+      50    HYDROGEN ON O+ OXYGEN
#   HO=+     52    HYDROGEN ON OXENIUM OXYGEN
#
sub _GetAtomTypeForHydrogenAttachedToOxygen {
  my($This, $OxygenAtom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

    # HOP : HYDROGEN ON OXYGEN ATTACHED TO PHOSPHOROUS
    if ($This->_IsHydrogenAttachedToOPOxygen($OxygenAtom)) {
      $AtomType = 'HOP';
      last ATOMTYPE;
    }

    # HOS : H ON OXYGEN ATTACHED TO SULFUR
    if ($This->_IsHydrogenAttachedToOSOxygen($OxygenAtom)) {
      $AtomType = 'HOS';
      last ATOMTYPE;
    }

    # HOCO : H-O IN CARBOXYLIC ACIDS
    if ($This->_IsHydrogenAttachedToCarboxylicAcidOxygen($OxygenAtom)) {
      $AtomType = 'HOCO';
      last ATOMTYPE;
    }

    # HOCC : H-O IN ENOLS AND PHENOLS
    if ($This->_IsHydrogenAttachedToEnolOrPhenolOxygen($OxygenAtom)) {
      $AtomType = 'HOCC';
      last ATOMTYPE;
    }

    # HOCN : H-O IN HO-C=N
    if ($This->_IsHydrogenAttachedToOCNOxygen($OxygenAtom)) {
      $AtomType = 'HOCN';
      last ATOMTYPE;
    }

    # HOM : HYDROGEN IN HYDROXIDE ANION
    if ($This->_IsHydrogenAttachedToHydroxideAnionOxygen($OxygenAtom)) {
      $AtomType = 'HOM';
      last ATOMTYPE;
    }

    # HO+ : HYDROGEN ON O+ OXYGEN
    if ($This->_IsHydrogenAttachedToPositivelyChargedOxygen($OxygenAtom)) {
      $AtomType = 'HO+';
      last ATOMTYPE;
    }

    # HO=+ : HYDROGEN ON OXENIUM OXYGEN
    if ($This->_IsHydrogenAttachedToOxeniumOxygen($OxygenAtom)) {
      $AtomType = 'HO=+';
      last ATOMTYPE;
    }

    # HOR : HYDROGEN IN ALCOHOLS
    if ($This->_IsHydrogenAttachedToAlcoholOxygen($OxygenAtom)) {
      $AtomType = 'HOR';
      last ATOMTYPE;
    }

    # HOH : HYDROGEN IN H2O
    if ($This->_IsHydrogenAttachedToWaterOxygen($OxygenAtom)) {
      $AtomType = 'HOH';
      last ATOMTYPE;
    }

    # HO: GENERAL H ON OXYGEN
    $AtomType = 'HO';
  }
  return $AtomType;
}

# HOR : HYDROGEN IN ALCOHOLS
#
sub _IsHydrogenAttachedToAlcoholOxygen {
  my($This, $OxygenAtom) = @_;

  return $This->_IsAlcoholOxygen($OxygenAtom) ? 1 : 0;
}

# HO : GENERAL H ON OXYGEN
#
sub _IsHydrogenAttachedToOxygen {
  my($This, $OxygenAtom) = @_;

  return 1;
}

# HOM : HYDROGEN IN HYDROXIDE ANION
#
sub _IsHydrogenAttachedToHydroxideAnionOxygen {
  my($This, $OxygenAtom) = @_;

  # H-(O-)
  return $OxygenAtom->DoesAtomNeighborhoodMatch('O.FC-1.T1', ['H'], ['-']) ? 1 : 0;
}

# HOCO : H-O IN CARBOXYLIC ACIDS
#
sub _IsHydrogenAttachedToCarboxylicAcidOxygen {
  my($This, $OxygenAtom) = @_;

  return $This->_IsEsterOrCarboxylicAcidOxygen($OxygenAtom) ? 1 : 0;
}

# HOP : HYDROGEN ON OXYGEN ATTACHED TO PHOSPHOROUS
#
sub _IsHydrogenAttachedToOPOxygen {
  my($This, $OxygenAtom) = @_;

  return $OxygenAtom->DoesAtomNeighborhoodMatch('O', ['P'], ['*']) ? 1 : 0;
}

# HOCC : H-O IN ENOLS AND PHENOLS
#
sub _IsHydrogenAttachedToEnolOrPhenolOxygen {
  my($This, $OxygenAtom) = @_;

  return $This->_IsEnolicOrPhenolicOxygen($OxygenAtom) ? 1 : 0;
}

# HOCN : H-O IN HO-C=N
#
sub _IsHydrogenAttachedToOCNOxygen {
  my($This, $OxygenAtom) = @_;

  return $This->_IsOCNDivalentOxygen($OxygenAtom) ? 1 : 0;
}

# HOH : HYDROGEN IN H2O
#
sub _IsHydrogenAttachedToWaterOxygen {
  my($This, $OxygenAtom) = @_;

  return $This->_IsWaterOxygen($OxygenAtom) ? 1 : 0;
}

# HOS : H ON OXYGEN ATTACHED TO SULFUR
#
sub _IsHydrogenAttachedToOSOxygen {
  my($This, $OxygenAtom) = @_;

  return $OxygenAtom->DoesAtomNeighborhoodMatch('O', ['S'], ['*']) ? 1 : 0;
}

# HO+ : HYDROGEN ON O+ OXYGEN
#
sub _IsHydrogenAttachedToPositivelyChargedOxygen {
  my($This, $OxygenAtom) = @_;

  return $This->_IsPositivelyChargedOxoniumOxygen($OxygenAtom) ? 1 : 0;
}

# HO=+ : HYDROGEN ON OXENIUM OXYGEN
#
sub _IsHydrogenAttachedToOxeniumOxygen {
  my($This, $OxygenAtom) = @_;

  return $This->_IsPositivelyChargedOxeniumOxygen($OxygenAtom) ? 1 : 0;
}

# Get MMFF94 atom type for Hydrogen attached to Phosphorus...
#
sub _GetAtomTypeForHydrogenAttachedToPhosphorus {
  my($This, $PhosphorusAtom) = @_;
  my($AtomType);

  # HP : H ATTACHED TO TRI- OR TETRACOORDINATE PHOSPHORUS
  $AtomType = 'HP';

  return $AtomType;
}

# Get MMFF94 atom type for Hydrogen attached to Sulfur...
#
sub _GetAtomTypeForHydrogenAttachedToSulfur {
  my($This, $SulfurAtom) = @_;
  my($AtomType);

  $AtomType = 'None';

  ATOMTYPE: {

   # HS=N : H ATTACHED TO TETRAVALENT, TRICOODR S DBL BONDED TO N
    if ($This->_IsSNTricoordinatedSulfur($SulfurAtom)) {
      $AtomType = 'HS=N';
      last ATOMTYPE;
    }

    # HS : H ATTACHED TO DIVALENT, DICOORDINATE S
    if ($This->_IsDivalentDiCoordinatedSulfur($SulfurAtom)) {
      $AtomType = 'HS';
      last ATOMTYPE;
    }

    $AtomType = 'None';
    carp "Warning: ${ClassName}->_GetAtomTypeForHydrogenAttachedToSulfur: MMFF94 atom type for Sulfur cann't be assigned...";
  }
  return $AtomType;
}

# Get MMFF94 atom type for Hydrogen attached to Silicon...
#
sub _GetAtomTypeForHydrogenAttachedToSilicon {
  my($This, $SiliconAtom) = @_;
  my($AtomType);

  # HSI : H ATTACHED TO SI
  $AtomType = 'HSI';

  return $AtomType;
}

# Get information about number and types of hetero atoms present in ring...
#
# Note:
#   . Any atom other than Carbon and Hydrogen atom is considered a hetero atom.
#
sub _GetHeteroAtomsInformationInRing {
  my($This, $RingAtomsRef) = @_;
  my($RingAtom, $RingAtomSymbol, $RingIsAromatic, $NumOfAromaticAtoms, $NumOfHeteroAtoms, %HeteroAtomSymbolsMap);

  %HeteroAtomSymbolsMap = ();

  $NumOfAromaticAtoms = 0;
  $NumOfHeteroAtoms = 0;

  RINGATOM: for $RingAtom (@{$RingAtomsRef}) {
    if ($RingAtom->IsAromatic()) {
      $NumOfAromaticAtoms++;
    }

    if (!$This->_IsHeteroAtom($RingAtom)) {
      next RINGATOM;
    }
    $NumOfHeteroAtoms++;

    $RingAtomSymbol = $RingAtom->GetAtomSymbol();
    if (exists $HeteroAtomSymbolsMap{$RingAtomSymbol}) {
      $HeteroAtomSymbolsMap{$RingAtomSymbol} += 1;
    }
    else {
      $HeteroAtomSymbolsMap{$RingAtomSymbol} = 1;
    }
  }
  $RingIsAromatic = ($NumOfAromaticAtoms == scalar @{$RingAtomsRef}) ? 1 : 0;

  return ($RingIsAromatic, $NumOfHeteroAtoms, \%HeteroAtomSymbolsMap);
}

# Check whether specified atom has a hetero atom at alpha position in an aromatic ring...
#
sub _IsAtomPositionAlphaInHeteroAromaticRing {
  my($This, $Atom, $RingAtomsRef) = @_;
  my($CheckRingAromaticity);

  $CheckRingAromaticity = 1;

  return $This->_IsAtomPositionAlphaInHeteroRing($Atom, $RingAtomsRef, $CheckRingAromaticity);
}

# Check whether specified atom has a hetero atom at beta position in an aromatic ring...
#
sub _IsAtomPositionBetaInHeteroAromaticRing {
  my($This, $Atom, $RingAtomsRef) = @_;
  my($CheckRingAromaticity);

  $CheckRingAromaticity = 1;

  return $This->_IsAtomPositionBetaInHeteroRing($Atom, $RingAtomsRef, $CheckRingAromaticity);
}

# Check whether specified atom has a hetero atom at alpha position in an aromatic or non-aromatic
# ring...
#
sub _IsAtomPositionAlphaInHeteroRing {
  my($This, $Atom, $RingAtomsRef, $CheckRingAromaticity) = @_;
  my($RingIsAromatic, $NumOfHeteroAtoms, $NumOfAllowedHeteroAtoms, $HeteroAtomPositionsRef);

  $CheckRingAromaticity = defined($CheckRingAromaticity) && $CheckRingAromaticity ? 1 : 0;

  # Is it an aromatic ring containing appropriate number hertero atoms?
  ($RingIsAromatic, $NumOfHeteroAtoms) = $This->_GetHeteroAtomsInformationInRing($RingAtomsRef);
  $NumOfAllowedHeteroAtoms = $Atom->IsCarbon() ? 1 : 2;

  if ($CheckRingAromaticity && !$RingIsAromatic) {
    return 0;
  }

  if ($NumOfHeteroAtoms != $NumOfAllowedHeteroAtoms) {
    return 0;
  }

  # Does ring contain hetero atoms at alpha position?
  $HeteroAtomPositionsRef = $This->_GetHeteroAtomPositionsInRing($Atom, $RingAtomsRef);
  if (exists($HeteroAtomPositionsRef->{Alpha}) && !exists($HeteroAtomPositionsRef->{Beta})) {
    return 1;
  }

  return 0;
}

# Check whether specified atom has a hetero atom at alpha position in an aromatic or non-aromatic
# ring...
#
sub _IsAtomPositionBetaInHeteroRing {
  my($This, $Atom, $RingAtomsRef, $CheckRingAromaticity) = @_;
  my($RingIsAromatic, $NumOfHeteroAtoms, $NumOfAllowedHeteroAtoms, $HeteroAtomPositionsRef);

  $CheckRingAromaticity = defined($CheckRingAromaticity) && $CheckRingAromaticity ? 1 : 0;

  # Is it an aromatic ring containing hertero atoms?
  ($RingIsAromatic, $NumOfHeteroAtoms) = $This->_GetHeteroAtomsInformationInRing($RingAtomsRef);
  $NumOfAllowedHeteroAtoms = $Atom->IsCarbon() ? 1 : 2;

  if ($CheckRingAromaticity && !$RingIsAromatic) {
    return 0;
  }

  if ($NumOfHeteroAtoms != $NumOfAllowedHeteroAtoms) {
    return 0;
  }

  # Does ring contain hetero atoms at alpha position?
  $HeteroAtomPositionsRef = $This->_GetHeteroAtomPositionsInRing($Atom, $RingAtomsRef);
  if (exists($HeteroAtomPositionsRef->{Beta}) && !exists($HeteroAtomPositionsRef->{Alpha})) {
    return 1;
  }

  return 0;
}

# Get hetro atom positions relative to atom position...
#
# Notes:
#    . Any atom other than Carbon and Hydrogen atom is considered a hetero atom.
#
sub _GetHeteroAtomPositionsInRing {
  my($This, $Atom, $RingAtomsRef) = @_;
  my($RingAtom, $Index, $AtomIndex, $NumOfHeteroAtoms, %HeteroAtomPositionsMap);

  %HeteroAtomPositionsMap = ();

  $NumOfHeteroAtoms = 0;
  $AtomIndex = 0;
  $Index = 0;

  # Find position of specified atom in the ring and count hetero atoms...
  for $RingAtom (@{$RingAtomsRef}) {
    if ($This->_IsHeteroAtom($RingAtom)) {
      $NumOfHeteroAtoms++;
    }
    if ($RingAtom->GetID() == $Atom->GetID()) {
      $AtomIndex = $Index;
    }
    $Index++;
  }

  # Does ring contain any hetereo atoms?
  if (!$NumOfHeteroAtoms) {
    return \%HeteroAtomPositionsMap;
  }

  # Check hetero atoms around specified atom to determine their position using their
  # their distance from specified atom: 1 - Alpha, 2 - Beta, 3 - Gamma, 4 - Delta, 5 - Omega
  #
  my($RingSize, $MaxPositionNum, $PositionNum, $PositionName, $MaxAtomIndex, $BeforeAtomIndex, $AfterAtomIndex, $BeforeAtom, $AfterAtom, %PositionNumToNameMap);

  %PositionNumToNameMap = ('1' => 'Alpha', '2' => 'Beta', '3' => 'Gamma', '4' => 'Delta', '5' => 'Omega');

  $RingSize = scalar @{$RingAtomsRef};
  $MaxPositionNum = int $RingSize/2;
  $MaxAtomIndex = $RingSize - 1;

  POSITIONNUM: for $PositionNum (1 .. $MaxPositionNum) {
    # Get atom before atom at a specific position...
    $BeforeAtomIndex = $AtomIndex - $PositionNum;
    if ($BeforeAtomIndex < 0) {
      $BeforeAtomIndex = $MaxAtomIndex + $BeforeAtomIndex + 1;
    }
    $BeforeAtom = $RingAtomsRef->[$BeforeAtomIndex];

    $PositionName = exists $PositionNumToNameMap{$PositionNum} ? $PositionNumToNameMap{$PositionNum} : 'Unknown';

    # Is atom before atom at a specific position a hetero atom?
    if (!$BeforeAtom->IsCarbon()) {
      $This->_TrackHeteroAtomPositionInRing($BeforeAtom, $PositionName, \%HeteroAtomPositionsMap);
    }

    # Get atom after atom at a specific position...
    $AfterAtomIndex = $AtomIndex + $PositionNum;
    if ($AfterAtomIndex > $MaxAtomIndex) {
      $AfterAtomIndex = $AfterAtomIndex - $MaxAtomIndex - 1;
    }

    # Is it a different atom?
    if ($AfterAtomIndex == $BeforeAtomIndex) {
      next POSITIONNUM;
    }

    # Is atom after atom at a specific position a hetero atom?
    $AfterAtom = $RingAtomsRef->[$AfterAtomIndex];

    if (!$AfterAtom->IsCarbon()) {
      $This->_TrackHeteroAtomPositionInRing($AfterAtom, $PositionName, \%HeteroAtomPositionsMap);
    }
  }
  return \%HeteroAtomPositionsMap;
}

# Is it a hetero atom?
#
# Notes:
#    . Any atom other than Carbon and Hydrogen atom is considered a hetero atom.
#
sub _IsHeteroAtom {
  my($This, $Atom) = @_;

  return ($Atom->IsCarbon() || $Atom->IsHydrogen()) ? 0 : 1;
}

# Track hetero atom positions in ring by updating data in specified data hash reference...
#
sub _TrackHeteroAtomPositionInRing {
  my($This, $HeteroAtom, $PositionName, $HeteroAtomPositionsMapRef) = @_;
  my($HeteroAtomSymbol);

  # Is it a new hetero atom position?
  if (!exists $HeteroAtomPositionsMapRef->{$PositionName}) {
    %{$HeteroAtomPositionsMapRef->{$PositionName}} = ();
  }

  $HeteroAtomSymbol = $HeteroAtom->GetAtomSymbol();
  if (exists $HeteroAtomPositionsMapRef->{$PositionName}{$HeteroAtomSymbol}) {
    $HeteroAtomPositionsMapRef->{$PositionName}{$HeteroAtomSymbol} += 1;
  }
  else {
    $HeteroAtomPositionsMapRef->{$PositionName}{$HeteroAtomSymbol} += 1;
  }
  return $This;
}

# Return a string containg data for MMFF94AtomTypes object...
#
sub StringifyMMFF94AtomTypes {
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

# Is it a MMFF94AtomTypes object?
sub _IsMMFF94AtomTypes {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Check and load MMFF94 atom types data...
#
sub _CheckAndLoadMMFF94AtomTypesData {

  # Is it already loaded?
  if (exists $MMFF94AtomTypesDataMap{AtomTypes}) {
    return;
  }

  _LoadMMFF94AtomTypesData();
}

# Load MMFF94 atom types data from the file assuming first column to be atom type symbol..
#
# Format:
#
# "AtomTypeSymbol","AtomTypeNum","ElementSymbol","AtomTypeDefinition"
# "CR","1","C","ALKYL CARBON, SP3"
# "C=C","2","C","VINYLIC CARBON, SP2"
#
sub _LoadMMFF94AtomTypesData {
  my($AtomTypesDataFile, $MayaChemToolsLibDir);

  $MayaChemToolsLibDir = FileUtil::GetMayaChemToolsLibDirName();

  $AtomTypesDataFile =  "$MayaChemToolsLibDir" . "/data/MMFF94AtomTypes.csv";
  if (! -e "$AtomTypesDataFile") {
    croak "Error: MayaChemTools package file, $AtomTypesDataFile, is missing: Possible installation problems...";
  }

  %MMFF94AtomTypesDataMap = ();
  AtomTypes::AtomTypes::LoadAtomTypesData($AtomTypesDataFile, \%MMFF94AtomTypesDataMap);
}

1;

__END__

=head1 NAME

MMFF94AtomTypes

=head1 SYNOPSIS

use AtomTypes::MMFF94AtomTypes;

use AtomTypes::MMFF94AtomTypes qw(:all);

=head1 DESCRIPTION

B<MMFF94AtomTypes> class provides the following methods:

new, AssignAtomTypes, GetAllPossibleMMFF94AtomTypes,
GetAllPossibleMMFF94NonHydrogenAtomTypes, GetMMFF94AtomTypesData,
StringifyMMFF94AtomTypes

The following functions are available:

GetAllPossibleMMFF94AtomTypes,
GetAllPossibleMMFF94NonHydrogenAtomTypes, GetMMFF94AtomTypesData

B<MMFF94AtomTypes> is derived from B<AtomTypes> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<MMFF94AtomTypes>, B<AtomTypes> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

The data file MMFF94AtomTypes.csv distributed with MayaChemTools release contains
all possible MMFF94 [ Ref 83-87 ] atom types.

Examples of MMFF94 atom types:

    CR, C=C, C=N, C=S, NR, N=C, OR, OC=O and so on

=head2 METHODS

=over 4

=item B<new>

    $NewMMFF94AtomTypes = new AtomTypes::MMFF94AtomTypes(%NamesAndValues);

Using specified I<MMFF94AtomTypes> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<MMFF94AtomTypes>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'MMFF94'
    IgnoreHydrogens = 0

Examples:

    $MMFF94AtomTypes = new AtomTypes::MMFF94AtomTypes(
                              'Molecule' => $Molecule,
                              'IgnoreHydrogens' => 0);

=item B<AssignAtomTypes>

    $MMFF94AtomTypes->AssignAtomTypes();

Assigns MMFF94 atom types to all the atoms in a molecule and returns
I<MMFF94AtomTypes>.

=item B<GetAllPossibleMMFF94AtomTypes>

    $AllAtomTypesDataRef = $MMFF94AtomTypes->
                           GetAllPossibleMMFF94AtomTypes();
    $AllAtomTypesDataRef = AtomTypes::MMFF94AtomTypes::
                           GetAllPossibleMMFF94AtomTypes();

Returns all possible MMFF94 atom types corresponding to hydrogen and non-hydrogen
atoms as an array reference.

=item B<GetAllPossibleMMFF94NonHydrogenAtomTypes>

    $AtomTypesDataRef = $MMFF94AtomTypes->
                        GetAllPossibleMMFF94NonHydrogenAtomTypes();
    $AtomTypesDataRef = AtomTypes::MMFF94AtomTypes::
                        GetAllPossibleMMFF94NonHydrogenAtomTypes();

Returns all possible MMFF94 atom types corresponding to non-hydrogen atoms as
an array reference.

=item B<GetMMFF94AtomTypesData>

    $AtomTypesDataMapRef = $MMFF94AtomTypes->GetMMFF94AtomTypesData();
    $AtomTypesDataMapRef = AtomTypes::MMFF94AtomTypes::
                           GetMMFF94AtomTypesData();

Returns MMFF94 atom types and associated data loaded from MMFF94 data file as
a reference to hash with the following hash data format:

    @{$MMFF94AtomTypesDataMap{AtomTypes}} - Array of all possible atom
                              types for all atoms
    @{$MMFF94AtomTypesDataMap{NonHydrogenAtomTypes}} - Array of all
                              possible atom types for non-hydrogen atoms
    @{$MMFF94AtomTypesDataMap->{ColLabels}} - Array of column labels
    %{$MMFF94AtomTypesDataMap->{DataCol<Num>}} - Hash keys pair:
                                                 DataCol<Num>, AtomType


=item B<StringifyMMFF94AtomTypes>

    $String = $MMFF94AtomTypes->StringifyMMFF94AtomTypes();

Returns a string containing information about I<MMFF94AtomTypes> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AtomTypes.pm, AtomicInvariantsAtomTypes.pm, DREIDINGAtomTypes.pm,
EStateAtomTypes.pm, FunctionalClassAtomTypes.pm, SLogPAtomTypes.pm,
SYBYLAtomTypes.pm, TPSAAtomTypes.pm, UFFAtomTypes.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
