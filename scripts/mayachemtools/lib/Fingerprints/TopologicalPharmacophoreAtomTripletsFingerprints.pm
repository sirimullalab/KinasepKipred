package Fingerprints::TopologicalPharmacophoreAtomTripletsFingerprints;
#
# File: TopologicalPharmacophoreAtomTripletsFingerprints.pm
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
use Fingerprints::Fingerprints;
use TextUtil ();
use MathUtil ();
use Molecule;
use AtomTypes::FunctionalClassAtomTypes;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Fingerprints::Fingerprints Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyTopologicalPharmacophoreAtomTripletsFingerprints';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeTopologicalPharmacophoreAtomTripletsFingerprints();

  $This->_InitializeTopologicalPharmacophoreAtomTripletsFingerprintsProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeTopologicalPharmacophoreAtomTripletsFingerprints {
  my($This) = @_;

  # Type of fingerprint...
  $This->{Type} = 'TopologicalPharmacophoreAtomTriplets';

  # Type of vector...
  $This->{VectorType} = 'FingerprintsVector';

  # AtomTripletsSetSizeToUse...
  #
  # ArbitrarySize - Corrresponds to atom triplets with non-zero count
  # FixedSize - Corresponds to all atom triplets with zero and non-zero count
  #
  # Possible values: ArbitrarySize or FixedSize. Default: ArbitrarySize
  #
  $This->{AtomTripletsSetSizeToUse} = '';

  #
  # OrderedNumericalValues - For ArbitrarySize value of AtomTripletsSetSizeToUse
  # NumericalValues - For FixedSize value of AtomTripletsSetSizeToUse
  #
  # Possible values: OrderedNumericalValues or NumericalValues. Default: NumericalValues
  #
  $This->{FingerprintsVectorType} = '';

  # Minimum and maximum bond distance between pharmacophore atom pairs corresponding to
  # atom triplets and distance bin size used for binning distances.
  #
  # In order to distribute distance bins of equal size, the last bin is allowed to go past the
  # maximum distance specified by upto distance bin size.
  #
  # The default MinDistance and MaxDistance values of 1 and 10 with DistanceBinSize of
  # 2 [ Ref 70 ] generates the following 5 distance bins: [1, 2] [3, 4] [5, 6] [7, 8] [9 10]
  #
  $This->{MinDistance} = 1;
  $This->{MaxDistance} = 10;

  # Distance bin size used for binning distances...
  #
  $This->{DistanceBinSize} = 2;

  # Determines whether to apply triangle inequality to distances triplets during basis set generation...
  #
  $This->{UseTriangleInequality} = 1;

  # Initialize pharmacophore atom types information...
  $This->_InitializeToplogicalPharmacophoreAtomTypesInformation();

  # Pharmacophore types assigned to each heavy atom...
  #
  %{$This->{AssignedAtomTypes}} = ();

  # All pharmacophore atom triplets between minimum and maximum distance...
  #
  %{$This->{AtomTriplets}} = ();
  @{$This->{AtomTriplets}{IDs}} = ();
  %{$This->{AtomTriplets}{Count}} = ();
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize object properties....
sub _InitializeTopologicalPharmacophoreAtomTripletsFingerprintsProperties {
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
  $This->_InitializeTopologicalPharmacophoreAtomTripletsFingerprintsVector();

  return $This;
}

# Initialize fingerprints vector...
#
sub _InitializeTopologicalPharmacophoreAtomTripletsFingerprintsVector {
  my($This) = @_;

  if (!$This->{AtomTripletsSetSizeToUse}) {
    $This->{AtomTripletsSetSizeToUse} =  'ArbitrarySize';
  }

  # Vector type and type of values...
  $This->{VectorType} = 'FingerprintsVector';

  if ($This->{AtomTripletsSetSizeToUse} =~ /^FixedSize$/i) {
    $This->{FingerprintsVectorType} = 'OrderedNumericalValues';
  }
  else {
    $This->{FingerprintsVectorType} = 'NumericalValues';
  }

  $This->_InitializeFingerprintsVector();
}

# Set atom parits set size to use...
#
sub SetAtomTripletsSetSizeToUse {
  my($This, $Value) = @_;

  if ($This->{AtomTripletsSetSizeToUse}) {
    croak "Error: ${ClassName}->SetAtomTripletsSetSizeToUse: Can't change size:  It's already set...";
  }

  if ($Value !~ /^(ArbitrarySize|FixedSize)$/i) {
    croak "Error: ${ClassName}->SetAtomTripletsSetSizeToUse: Unknown AtomTripletsSetSizeToUse value: $Value; Supported values: ArbitrarySize or FixedSize";
  }

  $This->{AtomTripletsSetSizeToUse} = $Value;

  return $This;
}

# Initialize topological atom types, generated by AtomTypes::FunctionalClassAtomTypes
# class, to use for atom triplets fingerprint generation...
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
#   Functiononal class atom type specification for an atom corresponds to:
#
#     Ar.CA.H.HBA.HBD.Hal.NI.PI.RA
#
#   Default pharmacophore atom types [ Ref 71 ] to use for atom triplets fingerprint generation
#   are: HBD, HBA, PI, NI, H, Ar
#
#   FunctionalAtomTypes are assigned using the following definitions [ Ref 60-61, Ref 65-66 ]:
#
#     HydrogenBondDonor: NH, NH2, OH
#     HydrogenBondAcceptor: N[!H], O
#     PositivelyIonizable: +, NH2
#     NegativelyIonizable: -, C(=O)OH, S(=O)OH, P(=O)OH
#
sub _InitializeToplogicalPharmacophoreAtomTypesInformation {
  my($This) = @_;

  #   Default pharmacophore atom types to use for atom triplets fingerprint generation
  #   are: HBD, HBA, PI, NI, H, Ar
  #
  @{$This->{AtomTypesToUse}} = ();
  @{$This->{AtomTypesToUse}} = sort ('HBD', 'HBA', 'PI', 'NI', 'H', 'Ar');

  return $This;
}

# Set atom types to use for atom triplets...
#
sub SetAtomTypesToUse {
  my($This, @Values) = @_;
  my($FirstValue, $TypeOfFirstValue, $AtomType, $SpecifiedAtomType, @SpecifiedAtomTypes, @AtomTypesToUse);

  if (!@Values) {
    carp "Warning: ${ClassName}->SetAtomTypesToUse: No values specified...";
    return;
  }

  $FirstValue = $Values[0];
  $TypeOfFirstValue = ref $FirstValue;

  @SpecifiedAtomTypes = ();
  @AtomTypesToUse = ();

  if ($TypeOfFirstValue =~ /^ARRAY/) {
    push @SpecifiedAtomTypes, @{$FirstValue};
  }
  else {
    push @SpecifiedAtomTypes, @Values;
  }

  # Make sure specified AtomTypes are valid...
  for $SpecifiedAtomType (@SpecifiedAtomTypes) {
    if (!AtomTypes::FunctionalClassAtomTypes::IsFunctionalClassAvailable($SpecifiedAtomType)) {
      croak "Error: ${ClassName}->SetAtomTypesToUse: Specified atom type, $SpecifiedAtomType, is not supported...\n ";
    }
    $AtomType = $SpecifiedAtomType;
    push @AtomTypesToUse, $AtomType;
  }

  # Set atom types to use...
  @{$This->{AtomTypesToUse}} = ();
  push @{$This->{AtomTypesToUse}}, sort @AtomTypesToUse;

  return $This;
}

# Set minimum distance for pharmacophore atom pairs in atom triplets...
#
sub SetMinDistance {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetMinDistance: MinDistance value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{MinDistance} = $Value;

  return $This;
}

# Set maximum distance for pharmacophore atom pairs in atom triplets...
#
sub SetMaxDistance {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetMaxDistance: MaxDistance value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{MaxDistance} = $Value;

  return $This;
}

# Set distance bin size for binning pharmacophore atom pair distances in atom triplets...
#
sub SetDistanceBinSize {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetDistanceBinSize: DistanceBinSize value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{DistanceBinSize} = $Value;

  return $This;
}

# Generate fingerprints description...
#
sub GetDescription {
  my($This) = @_;

  # Is description explicity set?
  if (exists $This->{Description}) {
    return $This->{Description};
  }

  # Generate fingerprints description...

  return "$This->{Type}:$This->{AtomTripletsSetSizeToUse}:MinDistance$This->{MinDistance}:MaxDistance$This->{MaxDistance}";
}

# Generate topological pharmacophore atom triplets [ Ref 66, Ref 68-71 ]  fingerprints...
#
# Let:
#
#   P = Any of the supported pharmacophore atom types
#
#   Px = Pharmacophore atom x
#   Py = Pharmacophore atom y
#   Pz = Pharmacophore atom z
#
#   Dxy = Distance or lower bound of binned distance between Px and Py
#   Dxz = Distance or lower bound of binned distance between Px and Pz
#   Dyz = Distance or lower bound of binned distance between Py and Pz
#
# Then:
#   PxDyz-PyDxz-PzDxy = Pharmacophore atom triplet ID for atoms Px, Py and Pz
#
# For example: H1-H1-H1, H2-HBA-H2 and so on
#
# Methodology:
#   . Generate a distance matrix.
#   . Using specified minimum, maximum and distance bin size, generate a binned distance
#     matrix from distance matrix. The lower distance bound on the distance bin is used
#     in the binned distance matrix and atom triplet IDs.
#   . Assign pharmacophore atom types to all the atoms.
#   . Initialize pharmacophore atom triplets basis set for all unique triplets constituting
#     atom pairs binned distances between minimum and maximum distance.
#       . Optionally, trinagle inequality is also implied which means:
#         . Distance or binned distance between any two pairs in a triplet must be less than the
#            sum of distances or binned distances between other two pairs and greater than the
#            difference of distances between other pairs.
#   . Using binned distance matrix and pharmacophore atom types, count occurance of
#     unique atom triplets.
#
# Notes:
#   . Hydrogen atoms are ignored during the fingerprint generation.
#
sub GenerateFingerprints {
  my($This) = @_;

  if ($This->{MinDistance} > $This->{MaxDistance}) {
    croak "Error: ${ClassName}->GenerateTopologicalPharmacophoreAtomTripletsFingerprints: No fingerpritns generated: MinDistance, $This->{MinDistance}, must be <= MaxDistance, $This->{MaxDistance}...";
  }

  # Cache appropriate molecule data...
  $This->_SetupMoleculeDataCache();

  # Generate distance matrix...
  if (!$This->_SetupDistanceMatrix()) {
    carp "Warning: ${ClassName}->GenerateFingerprints: Fingerprints generation didn't succeed: Couldn't generate distance matrix...";
    return $This;
  }

  # Generate binned distance matrix...
  $This->_GenerateBinnedDistanceMatrix();

  # Assign pharmacohore atom types to all heavy atoms...
  $This->_AssignPharmacophoreAtomTypes();

  # Initialize values of all possible pharmacohore atom triplets...
  $This->_InitializePharmacophoreAtomTriplets();

  # Count atom triplets...
  $This->_CountPharmacohoreAtomTriplets();

  # Set final fingerprints...
  $This->_SetFinalFingerprints();

  # Clear cached molecule data...
  $This->_ClearMoleculeDataCache();

  return $This;
}

# Setup distance matrix...
#
sub _SetupDistanceMatrix {
  my($This) = @_;

  $This->{DistanceMatrix} = $This->GetMolecule()->GetDistanceMatrix();

  if (!$This->{DistanceMatrix}) {
    return undef;
  }

  return $This;
}

# Generate binned distance matrix for distances with in the specified distance ranges...
#
sub _GenerateBinnedDistanceMatrix {
  my($This) = @_;
  my($DistanceMatrix, $BinnedDistanceMatrix, $NumOfRows, $NumOfCols, $RowIndex, $ColIndex, $SkipIndexCheck);

  $DistanceMatrix = $This->{DistanceMatrix};
  ($NumOfRows, $NumOfCols) = $DistanceMatrix->GetSize();

  # Initialize binned distance matrix...
  $BinnedDistanceMatrix = new Matrix($NumOfRows, $NumOfCols);

  # Setup distance to binned distance map...
  my($BinnedDistance, $Distance, %DistanceToBinnedDistance);
  %DistanceToBinnedDistance = ();
  for ($BinnedDistance = $This->{MinDistance}; $BinnedDistance <= $This->{MaxDistance}; $BinnedDistance += $This->{DistanceBinSize}) {
    for $Distance ($BinnedDistance .. ($BinnedDistance + $This->{DistanceBinSize} - 1)) {
      $DistanceToBinnedDistance{$Distance} = $BinnedDistance;
    }
  }

  # Generate binned distance matrix...
  $SkipIndexCheck = 0;
  for $RowIndex (0 .. ($NumOfRows - 1) ) {
    COLINDEX: for $ColIndex (($RowIndex + 1) .. ($NumOfCols - 1) ) {
      $Distance = $DistanceMatrix->GetValue($RowIndex, $ColIndex, $SkipIndexCheck);
      if ($Distance < $This->{MinDistance} || $Distance > $This->{MaxDistance}) {
	next COLINDEX;
      }
      $BinnedDistance = $DistanceToBinnedDistance{$Distance};
      $BinnedDistanceMatrix->SetValue($RowIndex, $ColIndex, $BinnedDistance, $SkipIndexCheck);
      $BinnedDistanceMatrix->SetValue($ColIndex, $RowIndex, $BinnedDistance, $SkipIndexCheck);
    }
  }

  $This->{BinnedDistanceMatrix} = $BinnedDistanceMatrix;

  return $This;
}

# Assign pharmacohore atom types to all heavy atoms...
#
sub _AssignPharmacophoreAtomTypes {
  my($This) = @_;
  my($Atom, $AtomID, $AtomType, $FunctionalClassAtomTypes);

  # Assign topological pharmacophore atom types...
  $FunctionalClassAtomTypes = new AtomTypes::FunctionalClassAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => 1, 'FunctionalClassesToUse' => $This->{AtomTypesToUse});
  $FunctionalClassAtomTypes->AssignAtomTypes();

  %{$This->{AssignedAtomTypes}} = ();

  ATOM: for $Atom (@{$This->{Atoms}}) {
    if ($Atom->IsHydrogen()) {
      next ATOM;
    }
    $AtomID = $Atom->GetID();

    my(@AtomTypes);
    @AtomTypes = ();

    $AtomType = $FunctionalClassAtomTypes->GetAtomType($Atom);
    if ($AtomType && $AtomType !~ /^None$/i) {
      push @AtomTypes, split /\./, $AtomType;
    }
    # Assign phramacophore types list to atom...
    $This->{AssignedAtomTypes}{$AtomID} = \@AtomTypes;
  }
  return $This;
}

# Initialize pharmacophore atom triplets basis set for all unique triplets constituting atom pairs
# binned distances between minimum and maximum distance and optionally applying triangle
# inequality. The DistanceBinSize determines the size of the distance bins. The lower distance
# bound, along with specified pharmacophore types, is used during generation of atom triplet
# IDs.
#
#
sub _InitializePharmacophoreAtomTriplets {
  my($This) = @_;
  my($AtomType1, $AtomType2, $AtomType3, $BinnedDistance12, $BinnedDistance13, $BinnedDistance23, $AtomTripletID);

  # Initialize atom triplets information...
  for ($BinnedDistance12 = $This->{MinDistance}; $BinnedDistance12 <= $This->{MaxDistance}; $BinnedDistance12 += $This->{DistanceBinSize}) {
    for ($BinnedDistance13 = $This->{MinDistance}; $BinnedDistance13 <= $This->{MaxDistance}; $BinnedDistance13 += $This->{DistanceBinSize}) {
      DISTANCE23: for ($BinnedDistance23 = $BinnedDistance12; $BinnedDistance23 <= $This->{MaxDistance}; $BinnedDistance23 += $This->{DistanceBinSize}) {
	if ($This->{UseTriangleInequality} && !$This->_DoDistancesSatisfyTriangleInequality($BinnedDistance12, $BinnedDistance13, $BinnedDistance23)) {
	  next DISTANCE23;
	}
	for $AtomType1 (@{$This->{AtomTypesToUse}}) {
	  for $AtomType2 (@{$This->{AtomTypesToUse}}) {
	    ATOMTYPE3: for $AtomType3 (@{$This->{AtomTypesToUse}}) {
	      $AtomTripletID = $This->_GetAtomTripletID($AtomType1, $BinnedDistance23, $AtomType2, $BinnedDistance13, $AtomType3, $BinnedDistance12);
	      if (exists $This->{AtomTriplets}{Count}{$AtomTripletID}) {
	      	next ATOMTYPE3;
	      }
	      # Unique atom triplets information...
	      push @{$This->{AtomTriplets}{IDs}}, $AtomTripletID;
	      $This->{AtomTriplets}{Count}{$AtomTripletID} = 0;
	    }
	  }
	}
      }
    }
  }
  return $This;
}

# Check triangle inequality...
#
sub _DoDistancesSatisfyTriangleInequality {
  my($This, $Distance1, $Distance2, $Distance3) = @_;

  if ( !($Distance1 > abs($Distance2 - $Distance3) && $Distance1 < ($Distance2 + $Distance3)) ) {
    return 0;
  }
  if ( !($Distance2 > abs($Distance1 - $Distance3) && $Distance2 < ($Distance1 + $Distance3)) ) {
    return 0;
  }
  if ( !($Distance3 > abs($Distance1 - $Distance2) && $Distance3 < ($Distance1 + $Distance2)) ) {
    return 0;
  }
  return 1;
}

# Count pharmacophore atom triplets...
#
sub _CountPharmacohoreAtomTriplets {
  my($This) = @_;
  my($NumOfAtoms, $AtomIndex1, $AtomIndex2, $AtomIndex3, $AtomID1, $AtomID2, $AtomID3, $AtomType1, $AtomType2, $AtomType3, $BinnedDistance12, $BinnedDistance13, $BinnedDistance23, $SkipIndexCheck, $BinnedDistanceMatrix, $AtomTripletID);

  $NumOfAtoms = @{$This->{Atoms}};
  $BinnedDistanceMatrix = $This->{BinnedDistanceMatrix};
  $SkipIndexCheck = 0;

  ATOMINDEX1: for $AtomIndex1 (0 .. ($NumOfAtoms - 1)) {
    $AtomID1 = $This->{AtomIndexToID}{$AtomIndex1};
    if ( !((exists($This->{AssignedAtomTypes}{$AtomID1}) && @{$This->{AssignedAtomTypes}{$AtomID1}})) ) {
      next ATOMINDEX1;
    }

    ATOMINDEX2: for $AtomIndex2 (($AtomIndex1 + 1) .. ($NumOfAtoms - 1)) {
      $AtomID2 = $This->{AtomIndexToID}{$AtomIndex2};
      if ( !((exists($This->{AssignedAtomTypes}{$AtomID2}) && @{$This->{AssignedAtomTypes}{$AtomID2}})) ) {
	next ATOMINDEX2;
      }
      $BinnedDistance12 = $BinnedDistanceMatrix->GetValue($AtomIndex1, $AtomIndex2, $SkipIndexCheck);
      if ($BinnedDistance12 == 0) {
	next ATOMINDEX2;
      }

      ATOMINDEX3: for $AtomIndex3 (($AtomIndex2 + 1) .. ($NumOfAtoms - 1)) {
	$AtomID3 = $This->{AtomIndexToID}{$AtomIndex3};
	if ( !((exists($This->{AssignedAtomTypes}{$AtomID3}) && @{$This->{AssignedAtomTypes}{$AtomID3}})) ) {
	  next ATOMINDEX3;
	}
	$BinnedDistance13 = $BinnedDistanceMatrix->GetValue($AtomIndex1, $AtomIndex3, $SkipIndexCheck);
	$BinnedDistance23 = $BinnedDistanceMatrix->GetValue($AtomIndex2, $AtomIndex3, $SkipIndexCheck);
	if ($BinnedDistance13 == 0 || $BinnedDistance23 == 0) {
	  next ATOMINDEX3;
	}
	if ($This->{UseTriangleInequality} && !$This->_DoDistancesSatisfyTriangleInequality($BinnedDistance12, $BinnedDistance13, $BinnedDistance23)) {
	  next ATOMINDEX3;
	}

	# Go over possible pharmacohore triplets for the three pharmacophore atoms using the
	# binned distances...
	for $AtomType1 (@{$This->{AssignedAtomTypes}{$AtomID1}}) {
	  for $AtomType2 (@{$This->{AssignedAtomTypes}{$AtomID2}}) {
	    for $AtomType3 (@{$This->{AssignedAtomTypes}{$AtomID3}}) {
	      $AtomTripletID = $This->_GetAtomTripletID($AtomType1, $BinnedDistance23, $AtomType2, $BinnedDistance13, $AtomType3, $BinnedDistance12);
	      $This->{AtomTriplets}{Count}{$AtomTripletID} += 1;
	    }
	  }
	}
      }
    }
  }
  return $This;
}

# Set final fingerpritns vector...
#
sub _SetFinalFingerprints {
  my($This) = @_;
  my($UseArbitrarySetSize, $ID, $Value, @IDs, @Values);

  # Mark successful generation of fingerprints...
  $This->{FingerprintsGenerated} = 1;

  # Is it an ArbitraySize atom triplets set size?
  $UseArbitrarySetSize = $This->{AtomTripletsSetSizeToUse} =~ /^ArbitrarySize$/i ? 1 : 0;

  # Set atom triplet count values...
  @IDs = (); @Values = ();

  if ($UseArbitrarySetSize) {
    ID: for $ID (@{$This->{AtomTriplets}{IDs}}) {
      $Value = $This->{AtomTriplets}{Count}{$ID};
      if ($Value == 0) {
	next ID;
      }
      push @IDs, $ID;
      push @Values, $Value;
    }
  }
  else {
    @Values = map { $This->{AtomTriplets}{Count}{$_} } @{$This->{AtomTriplets}{IDs}};
  }

  # Set atom triplet IDs for fingerprint vector...
  if ($UseArbitrarySetSize) {
    $This->{FingerprintsVector}->AddValueIDs(\@IDs);
  }
  else {
    $This->{FingerprintsVector}->AddValueIDs(\@{$This->{AtomTriplets}{IDs}});
  }

  # Set atom triplets count values for fingerprint vector...
  $This->{FingerprintsVector}->AddValues(\@Values);

  return $This;
}

# Return an array or reference to an array containing atom triplet IDs...
#
sub GetAtomTripletIDs {
  my($This) = @_;

  return wantarray ? @{$This->{AtomTriplets}{IDs}} : \@{$This->{AtomTriplets}{IDs}};
}

# Get pharmacophore atom triplet ID corresponding to atom types and distances
# corresponding to atom triplet...
#
sub _GetAtomTripletID {
  my($This, $Px, $Dyz, $Py, $Dxz, $Pz, $Dxy) = @_;
  my($AtomTripletID, @AtomIDs);

  @AtomIDs = ();

  @AtomIDs = sort("${Px}${Dyz}", "${Py}${Dxz}", "${Pz}${Dxy}");
  $AtomTripletID = join "-", @AtomIDs;

  return $AtomTripletID;
}

# Cache  appropriate molecule data...
#
sub _SetupMoleculeDataCache {
  my($This) = @_;

  # Get all atoms including hydrogens to correctly map atom indices to atom IDs for
  # usage of distance matrix. The hydrogen atoms are ignored during processing...
  #
  @{$This->{Atoms}} = $This->GetMolecule()->GetAtoms();

  # Get all atom IDs...
  my(@AtomIDs);
  @AtomIDs = ();
  @AtomIDs =  map { $_->GetID() } @{$This->{Atoms}};

  # Set AtomIndex to AtomID hash...
  %{$This->{AtomIndexToID}} = ();
  @{$This->{AtomIndexToID}}{ (0 .. $#AtomIDs) } = @AtomIDs;

  return $This;
}

# Clear cached molecule data...
#
sub _ClearMoleculeDataCache {
  my($This) = @_;

  @{$This->{Atoms}} = ();

  return $This;
}


# Return a string containg data for TopologicalPharmacophoreAtomTripletsFingerprints object...
#
sub StringifyTopologicalPharmacophoreAtomTripletsFingerprints {
  my($This) = @_;
  my($FingerprintsString, $UseTriangleInequality);

  # Type of fingerprint...
  $FingerprintsString = "Fingerprint type: $This->{Type}; AtomTripletsSetSizeToUse: $This->{AtomTripletsSetSizeToUse}";

  # Distances information...
  $FingerprintsString .= "; MinDistance:  $This->{MinDistance}; MaxDistance: $This->{MaxDistance}; DistanceBinSize: $This->{DistanceBinSize}; UseTriangleInequality: " . ($This->{UseTriangleInequality} ? "Yes" : "No");

  # Pharmacophore atom type labels and description...
  my($AtomType, @AtomTypes, @AtomTypesOrder, %AvailableAtomTypes);

  @AtomTypesOrder = AtomTypes::FunctionalClassAtomTypes::GetFunctionalClassesOrder();
  %AvailableAtomTypes = AtomTypes::FunctionalClassAtomTypes::GetAvailableFunctionalClasses();

  @AtomTypes = ();
  for $AtomType (@AtomTypesOrder) {
    push @AtomTypes, "$AtomType: $AvailableAtomTypes{$AtomType}";
  }

  $FingerprintsString .= "; AtomTypesToUse: <" . TextUtil::JoinWords(\@{$This->{AtomTypesToUse}}, ", ", 0) . ">";
  $FingerprintsString .= "; AtomTypesOrder: <" . TextUtil::JoinWords(\@AtomTypesOrder, ", ", 0) . ">";
  $FingerprintsString .= "; AvailableAtomTypes: <" . TextUtil::JoinWords(\@AtomTypes, ", ", 0) . ">";

  # Total number of pharmacophore atom triplets...
  $FingerprintsString .= "; NumOfAtomTriplets: " . $This->{FingerprintsVector}->GetNumOfValues();

  # FingerprintsVector...
  $FingerprintsString .= "; FingerprintsVector: < $This->{FingerprintsVector} >";

  return $FingerprintsString;
}

1;

__END__

=head1 NAME

TopologicalPharmacophoreAtomTripletsFingerprints

=head1 SYNOPSIS

use Fingerprints::TopologicalPharmacophoreAtomTripletsFingerprints;

use Fingerprints::TopologicalPharmacophoreAtomTripletsFingerprints qw(:all);

=head1 DESCRIPTION

B<TopologicalPharmacophoreAtomTripletsFingerprints> [ Ref 66, Ref 68-71 ] class provides
the following methods:

new, GenerateFingerprints, , GetDescription, GetAtomTripletIDs,
SetAtomTypesToUse, SetDistanceBinSize, SetMaxDistance, SetMinDistance,
StringifyTopologicalPharmacophoreAtomTripletsFingerprints

B<TopologicalPharmacophoreAtomTripletsFingerprints> is derived from B<Fingerprints> class
which in turn is  derived from B<ObjectProperty> base class that provides methods not explicitly
defined in B<TopologicalPharmacophoreAtomTripletsFingerprints>, B<Fingerprints> or B<ObjectProperty>
classes using Perl's AUTOLOAD functionality. These methods are generated on-the-fly for a specified
object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

Based on the values specified for B<AtomTypesToUse>, pharmacophore atom types are
assigned to all non-hydrogen atoms in a molecule and a distance matrix is generated.
Using B<MinDistance>, B<MaxDistance>, and B<DistanceBinSize> values, a
binned distance matrix is generated with lower bound on the distance bin as the distance
in distance matrix; the lower bound on the distance bin is also used as the distance between
atom pairs for generation of atom triplet identifiers.

A pharmacophore atom triplets basis set is generated for all unique atom triplets constituting
atom pairs binned distances between B<--MinDistance> and B<--MaxDistance>. The value
of B<--UseTriangleInequality> determines whether the triangle inequality test is applied during
generation of atom triplets basis set. The lower distance bound, along with specified pharmacophore
types, is used during generation of atom triplet IDs.

    Let:

    P = Valid pharmacophore atom type

    Px = Pharmacophore atom x
    Py = Pharmacophore atom y
    Pz = Pharmacophore atom z

    Dmin = Minimum distance corresponding to number of bonds between two atoms
    Dmax = Maximum distance corresponding to number of bonds between two atoms
    D = Distance corresponding to number of bonds between two atom

    Bsize  = Distance bin size
    Nbins = Number of distance bins

    Dxy = Distance or lower bound of binned distance between Px and Py
    Dxz = Distance or lower bound of binned distance between Px and Pz
    Dyz = Distance or lower bound of binned distance between Py and Pz

    Then:

    PxDyz-PyDxz-PzDxy = Pharmacophore atom triplet IDs for atom types Px,
                        Py, and Pz

    For example: H1-H1-H1, H2-HBA-H2 and so on.

    For default values of Dmin = 1 , Dmax = 10 and Bsize = 2, the number of
    distance bins, Nbins = 5, are:

    [1, 2] [3, 4] [5, 6] [7, 8] [9 10]

    and atom triplet basis set size is 2692.

    Atom triplet basis set size for various values of Dmin, Dmax and Bsize in
    conjunction with usage of triangle inequality is:

    Dmin    Dmax   Bsize   UseTriangleInequality   TripletBasisSetSize
    1       10     2       No                      4960
    1       10     2       Yes                     2692 [ Default ]
    2       12     2       No                      8436
    2       12     2       Yes                     4494


Using binned distance matrix and pharmacohore atom types, occurrence of unique pharmacohore
atom triplets is counted.

The final pharmacophore atom triples count along with atom pair identifiers involving all non-hydrogen
atoms constitute pharmacophore topological atom triplets fingerprints of the molecule.

For I<ArbitrarySize> value of B<AtomTripletsSetSizeToUse>, the fingerprint vector correspond to
only those topological pharmacophore atom triplets which are present and have non-zero count. However,
for I<FixedSize> value of B<AtomTripletsSetSizeToUse>, the fingerprint vector contains all possible
valid topological pharmacophore atom triplets with both zero and non-zero count values.

The current release of MayaChemTools generates the following types of topological pharmacophore
atom triplets fingerprints vector strings:

    FingerprintsVector;TopologicalPharmacophoreAtomTriplets:ArbitrarySize:
    MinDistance1:MaxDistance10;696;NumericalValues;IDsAndValuesString;Ar1-
    Ar1-Ar1 Ar1-Ar1-H1 Ar1-Ar1-HBA1 Ar1-Ar1-HBD1 Ar1-H1-H1 Ar1-H1-HBA1 Ar1
    -H1-HBD1 Ar1-HBA1-HBD1 H1-H1-H1 H1-H1-HBA1 H1-H1-HBD1 H1-HBA1-HBA1 H1-
    HBA1-HBD1 H1-HBA1-NI1 H1-HBD1-NI1 HBA1-HBA1-NI1 HBA1-HBD1-NI1 Ar1-...;
    46 106 8 3 83 11 4 1 21 5 3 1 2 2 1 1 1 100 101 18 11 145 132 26 14 23
    28 3 3 5 4 61 45 10 4 16 20 7 5 1 3 4 5 3 1 1 1 1 5 4 2 1 2 2 2 1 1 1
    119 123 24 15 185 202 41 25 22 17 3 5 85 95 18 11 23 17 3 1 1 6 4 ...

    FingerprintsVector;TopologicalPharmacophoreAtomTriplets:FixedSize:MinD
    istance1:MaxDistance10;2692;OrderedNumericalValues;ValuesString;46 106
    8 3 0 0 83 11 4 0 0 0 1 0 0 0 0 0 0 0 0 21 5 3 0 0 1 2 2 0 0 1 0 0 0
    0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 101 18 11 0 0 145 132 26
    14 0 0 23 28 3 3 0 0 5 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 61 45 10 4 0
    0 16 20 7 5 1 0 3 4 5 3 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 5 ...

    FingerprintsVector;TopologicalPharmacophoreAtomTriplets:FixedSize:MinD
    istance1:MaxDistance10;2692;OrderedNumericalValues;IDsAndValuesString;
    Ar1-Ar1-Ar1 Ar1-Ar1-H1 Ar1-Ar1-HBA1 Ar1-Ar1-HBD1 Ar1-Ar1-NI1 Ar1-Ar1-P
    I1 Ar1-H1-H1 Ar1-H1-HBA1 Ar1-H1-HBD1 Ar1-H1-NI1 Ar1-H1-PI1 Ar1-HBA1-HB
    A1 Ar1-HBA1-HBD1 Ar1-HBA1-NI1 Ar1-HBA1-PI1 Ar1-HBD1-HBD1 Ar1-HBD1-...;
    46 106 8 3 0 0 83 11 4 0 0 0 1 0 0 0 0 0 0 0 0 21 5 3 0 0 1 2 2 0 0 1
    0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 101 18 11 0 0 145
    132 26 14 0 0 23 28 3 3 0 0 5 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 61 ...

=head2 METHODS

=over 4

=item B<new>

    $TPATFP = new TopologicalPharmacophoreAtomTripletsFingerprints(
                                                   %NamesAndValues);

Using specified I<TopologicalPharmacophoreAtomTripletsFingerprints> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<TopologicalPharmacophoreAtomTripletsFingerprints>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'TopologicalPharmacophoreAtomTriplets'
    MinDistance = 1
    MaxDistance = 10
    DistanceBinSize = 2
    UseTriangleInequality = 1
    AtomTypesToUse = ['HBD', 'HBA', 'PI', 'NI', 'H', 'Ar']

Examples:

    $TPATFP = new TopologicalPharmacophoreAtomTripletsFingerprints(
                              'Molecule' => $Molecule);

    $TPATFP = new TopologicalPharmacophoreAtomTripletsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomTripletsSetSizeToUse' => 'ArbitrarySize';
                              'MinDistance' => 1,
                              'MaxDistance' => 10,
                              'DistanceBinSize' => 2,
                              'AtomTypesToUse' => ['HBD', 'HBA', 'PI', 'NI', 'H', 'Ar'],
                              'UseTriangleInequality' => 1);

    $TPATFP = new TopologicalPharmacophoreAtomTripletsFingerprints(
                              'Molecule' => $Molecule,
                              'AtomTripletsSetSizeToUse' => 'FixedSize';
                              'MinDistance' => 1,
                              'MaxDistance' => 10,
                              'DistanceBinSize' => 2,
                              'AtomTypesToUse' => ['HBD', 'HBA', 'PI', 'NI', 'H', 'Ar'],
                              'UseTriangleInequality' => 1);

    $TPATFP->GenerateFingerprints();
    print "$TPATFP\n";

=item B<GetDescription>

    $Description = $TopologicalPharmacophoreAtomTripletsFP->GetDescription();

Returns a string containing description of topological pharmacophore atom triplets fingerprints.

=item B<GenerateFingerprints>

    $TopologicalPharmacophoreAtomTripletsFP->GenerateFingerprints();

Generates topological pharmacophore atom triplets fingerprints and returns
I<TopologicalPharmacophoreAtomTripletsFP>.

=item B<GetAtomTripletIDs>

    $AtomTripletsIDsRef = $TopologicalPharmacophoreATFP->GetAtomTripletIDs();
    @AtomTripletIDs = $TopologicalPharmacophoreATFP->GetAtomTripletIDs();

Returns atom triplet IDs corresponding to atom pairs count values in topological pharmacophore
atom triplet fingerprints vector as an array or reference to an array.

=item B<AtomTripletsSetSizeToUse>

    $TPAFP->AtomTripletsSetSizeToUse($Values);

Sets pharmacophore atom triplets set size to use for topological pharmacophore fingerprints
generation and returns I<TopologicalPharmacophoreAtomTripletsFingerprints>.

Possible values for pharmacophore atom triplets set size are: I<ArbitrarySize, FizedSize>.
Default value: I<ArbitrarySize>.

For I<ArbitrarySize> value of B<AtomTripletsSetSizeToUse>, the fingerprint vector correspond to
only those topological pharmacophore atom triplets which are present and have non-zero count. However,
for I<FixedSize> value of B<AtomTripletsSetSizeToUse>, the fingerprint vector contains all possible
valid topological pharmacophore atom triplets with both zero and non-zero count values.

=item B<SetAtomTypesToUse>

    $TopologicalPharmacophoreAtomTripletsFP->SetAtomTypesToUse($ValuesRef);
    $TopologicalPharmacophoreAtomTripletsFP->SetAtomTypesToUse(@Values);

Sets pharmacophore atom types to use for topological pharmacophore fingerprints
generation and returns I<TopologicalPharmacophoreAtomTripletsFingerprints>.

Possible values for pharmacophore atom types are: I<Ar, CA, H, HBA, HBD, Hal, NI, PI, RA>.
Default value [ Ref 71 ] : I<HBD,HBA,PI,NI,H,Ar>.

The pharmacophore atom types abbreviations correspond to:

    HBD: HydrogenBondDonor
    HBA: HydrogenBondAcceptor
    PI :  PositivelyIonizable
    NI : NegativelyIonizable
    Ar : Aromatic
    Hal : Halogen
    H : Hydrophobic
    RA : RingAtom
    CA : ChainAtom

I<AtomTypes::FunctionalClassAtomTypes> module is used to assign pharmacophore atom
types. It uses following definitions [ Ref 60-61, Ref 65-66 ]:

    HydrogenBondDonor: NH, NH2, OH
    HydrogenBondAcceptor: N[!H], O
    PositivelyIonizable: +, NH2
    NegativelyIonizable: -, C(=O)OH, S(=O)OH, P(=O)OH


=item B<SetDistanceBinSize>

    $TopologicalPharmacophoreAtomTripletsFP->SetDistanceBinSize($Value);

Sets distance bin size used to bin distances between atom pairs in atom triplets and returns
I<TopologicalPharmacophoreAtomTriplesFP>.

For default B<MinDistance> and B<MaxDistance> values of 1 and 10 with  B<DistanceBinSize>
of 2 [ Ref 70 ], the following 5 distance bins are generated:

    [1, 2] [3, 4] [5, 6] [7, 8] [9 10]

The lower distance bound on the distance bin is uses to bin the distance between atom pairs in
atom triplets. So in the previous example, atom pairs with distances 1 and 2 fall in first distance
bin, atom pairs with distances 3 and 4  fall in second distance bin and so on.

In order to distribute distance bins of equal size, the last bin is allowed to go past B<MaxDistance>
by up to distance bin size. For example, B<MinDistance> and B<MaxDistance> values of 2 and 10
with B<DistanceBinSize> of 2 generates the following 6 distance bins:

    [2, 3] [4, 5] [6, 7] [8, 9] [10 11]


=item B<SetMaxDistance>

    $TopologicalPharmacophoreAtomTriplesFP->SetMaxDistance($Value);

Sets maximum bond distance between atom pairs  corresponding to atom triplets for
generating topological pharmacophore atom triplets fingerprints and returns
I<TopologicalPharmacophoreAtomTriplesFP>.

=item B<SetMinDistance>

    $TopologicalPharmacophoreAtomTriplesFP->SetMinDistance($Value);

Sets minimum bond distance between atom pairs  corresponding to atom triplets for
generating topological pharmacophore atom triplets fingerprints and returns
I<TopologicalPharmacophoreAtomTriplesFP>.

=item B<StringifyTopologicalPharmacophoreAtomTripletsFingerprints>

    $String = $TopologicalPharmacophoreAtomTripletsFingerprints->
                   StringifyTopologicalPharmacophoreAtomTripletsFingerprints();

Returns a string containing information about I<TopologicalPharmacophoreAtomTripletsFingerprints> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Fingerprints.pm, FingerprintsStringUtil.pm, AtomNeighborhoodsFingerprints.pm,
AtomTypesFingerprints.pm, EStateIndiciesFingerprints.pm, ExtendedConnectivityFingerprints.pm,
MACCSKeys.pm, PathLengthFingerprints.pm, TopologicalAtomPairsFingerprints.pm,
TopologicalAtomTripletsFingerprints.pm, TopologicalAtomTorsionsFingerprints.pm,
TopologicalPharmacophoreAtomPairsFingerprints.pm,

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
