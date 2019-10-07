package AtomicDescriptors::EStateValuesDescriptors;
#
# File: EStateValuesDescriptors.pm
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
use Matrix;
use Constants;
use TextUtil ();
use MathUtil ();
use StatisticsUtil ();
use Atom;
use Molecule;
use AtomicDescriptors::AtomicDescriptors;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(AtomicDescriptors::AtomicDescriptors Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyEStateValuesDescriptors';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeEStateValuesDescriptors();

  $This->_InitializeEStateValuesDescriptorsProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}


# Initialize object data...
#
sub _InitializeEStateValuesDescriptors {
  my($This) = @_;

  # Type of AtomicDescriptor...
  $This->{Type} = 'EStateValue';

  # Intrinsic state values calculated for non-hydrogen atom...
  #
  %{$This->{IStateValues}} = ();

  # Calculatetion of E-state values for types for hydrogens is not supported...
  $This->{IgnoreHydrogens} = 1;

  # Perturbation to intrinsic state values of atoms corresponding to all non-hydrogen
  # atom pairs...
  #
  %{$This->{DeltaIStateMatrix}} = ();

  # E-state values calculated for non-hydrogen atoms...
  #
  %{$This->{EStateValues}} = ();

  return $This;
}

# Initialize object properties...
#
sub _InitializeEStateValuesDescriptorsProperties {
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

  # Intialize atomic descriptor values...
  $This->_InitializeDescriptorValues();

  return $This;
}

# Disable change of ignore hydrogens...
#
sub SetIgnoreHydrogens {
  my($This, $IgnoreHydrogens) = @_;

  carp "Warning: ${ClassName}->SetIgnoreHydrogens: Ignore hydrogens value can't be changed: It's not supported...";

  return $This;
}

# Generate electrotopological state (E-state) values [ Ref 75-78 ] for all atoms
# in the molecule...
#
# Calculation of E-state values for non-hydrogen atoms:
#
# Let:
#
#   N = Principal quantum number or period number corresponding to element symbol
#
#   Sigma = Number of sigma electrons involves in bonds to hydrogen and non-hydrogen atoms
#           attached to atom
#         = Number of sigma bonds to hydrogen and non-hydrogen atoms attached to atom
#   PI = Number of PI electrons involved in bonds to non-hydrogen atoms attached to atom
#       = Number of PI bonds to non-hydrogen atoms attached to atom
#
#   LP = Number of lone pair electrons on atom
#
#   Zv = Number of electrons in valence shell of atom
#
#   X = Number of non-hydrogen atom neighbors or heavy atoms attached to atom
#   H = Number of implicit and explicit hydrogens for atom
#
#   Delta = Number of sigma electrons involved to bonds to non-hydrogen atoms
#   DeltaV = ValenceDelta = Number of valence shell electrons not involved in bonding to hydrogen atoms
#
#   Ii = Intrinsic state value for atom i
#
#   DeltaIi = Sum of perturbations to intrinsic state value Ii of atom i by all other atoms besides atom i
#
#   DeltaIij = Perturbation to intrinsic state value Ii of atom i by atom j
#
#   Dij = Graph/bond distance between atom i and j
#   Rij = Dij + 1
#
#   Si = E-state value for atom i
#
#
# Then:
#
#   Delta = Sigma - H = X
#
#   DeltaV = Zv - H
#      = Sigma + PI + LP - H
#
#   Ii = ( ( ( 2 / N ) ** 2 ) * DeltaV + 1 ) / Delta
#
#   Si = Ii + DeltaIi
#
#   DeltaIij = (Ii - Ij) / (Rij ** 2) for j not equal to i
#
#   DeltaIji = - DeltaIij
#
#   DeltaIi = SUM ( (Ii - Ij) / (Rij ** 2) ) for j = 1 to num of atoms skipping atom i
#
# Methodology:
#   . Calculate intrinsic state values for atoms.
#   . Generate a distance matrix.
#   . Use distance matrix to calculate DeltaIij matrix with each row i containing
#     DeltaIij values corresponding to perturbation to intrinsic state value of atom
#     i by atom j.
#   . Calculate E-state values using DeltaIij matrix.
#   . Assign E-state values to atoms.
#
# Notes:
#   . The current release of MayaChemTools doesn't support calculation of Hydrogen
#     E-state values.
#
sub GenerateDescriptors {
  my($This) = @_;

  # Cache appropriate molecule data...
  $This->_SetupMoleculeDataCache();

  # Generate distance matrix...
  if (!$This->_SetupDistanceMatrix()) {
    carp "Warning: ${ClassName}->GenerateDescriptors: E-state values description generation didn't succeed: Couldn't generate a distance matrix...";
    return $This;
  }

  # Calculate EState values..
  if (!$This->_CalculateEStateValuesDescriptors()) {
    carp "Warning: ${ClassName}->GenerateDescriptors: E-state values description generation didn't succeed: Couldn't calculate IState values for all atoms...";
    return $This;
  }

  # Set final descriptor values...
  $This->_SetFinalValues();

  # Clear cached molecule data...
  $This->_ClearMoleculeDataCache();

  return $This;
}

# Calculate E-state values for non-hydrogen atoms...
#
sub _CalculateEStateValuesDescriptors {
  my($This) = @_;
  my($DeltaIStateMatrix, $NumOfRows, $NumOfCols, $RowIndex, $AtomID, $EStateValue, $IStateValue, $DeltaIStateValue, @DeltaIStateMatrixRowValues);

  %{$This->{EStateValues}} = ();

  # Calculate intrinsic state values for non-hydrogen atoms...
  if (!$This->_CalculateIStateValues()) {
    return undef;
  }

  # Calculate delta intrinsic state matrix for non-hydrogen atoms...
  $This->_CalculateDeltaIStateMatrix();

  # Get DeltaIState matrix information...
  $DeltaIStateMatrix = $This->{DeltaIStateMatrix};
  ($NumOfRows, $NumOfCols) = $DeltaIStateMatrix->GetSize();

  # Calculate EState values...
  ROWINDEX: for $RowIndex (0 .. ($NumOfRows - 1) ) {
    $AtomID = $This->{AtomIndexToID}{$RowIndex};
    if ( !(exists($This->{IStateValues}{$AtomID})) ) {
      next ROWINDEX;
    }
    $IStateValue = $This->{IStateValues}{$AtomID};

    @DeltaIStateMatrixRowValues = $DeltaIStateMatrix->GetRowValues($RowIndex);
    $DeltaIStateValue = StatisticsUtil::Sum(\@DeltaIStateMatrixRowValues);

    $EStateValue = $IStateValue + $DeltaIStateValue;

    $This->{EStateValues}{$AtomID} = $EStateValue;
  }
  return $This;
}

# Calculate intrinsic state values for non-hydrogen atoms...
#
sub _CalculateIStateValues {
  my($This) = @_;
  my($Atom, $AtomID);

  %{$This->{IStateValues}} = ();

  ATOM: for $Atom (@{$This->{Atoms}}) {
    # Irrespective of IgoreHydrogens value, just ignore hydrogens...
    if ($Atom->IsHydrogen()) {
      next ATOM;
    }
    $AtomID = $Atom->GetID();
    $This->{IStateValues}{$AtomID} = $This->_CalculateIStateValue($Atom);
    if (!defined($This->{IStateValues}{$AtomID})) {
      return undef;
    }
  }
  return $This;
}

# Calculation intrinsic state value for non-hydrogen...
#
sub _CalculateIStateValue {
  my($This, $Atom) = @_;
  my($IStateValue, $Delta, $DeltaV, $PeriodNumber);

  $PeriodNumber = $Atom->GetPeriodNumber();
  ($Delta, $DeltaV) = $This->_GetDeltaValues($Atom);

  if (!(defined($Delta) && defined($PeriodNumber) && defined($DeltaV))) {
    return undef;
  }

  $IStateValue = ($PeriodNumber && $Delta) ? (((2/$PeriodNumber)**2)*$DeltaV + 1)/$Delta : 0;

  return $IStateValue;
}

# Get Delta and DeltaV values for atom...
#
sub _GetDeltaValues {
  my($This, $Atom) = @_;
  my($Delta, $DeltaV, $ValenceElectrons, $NumOfHydrogens);

  ($Delta, $DeltaV) = (undef, undef);

  $ValenceElectrons = $Atom->GetValenceElectrons();
  $NumOfHydrogens = $Atom->GetAtomicInvariantValue('H');

  $Delta = $Atom->GetAtomicInvariantValue('X');
  if (defined($ValenceElectrons) && defined($NumOfHydrogens)) {
    $DeltaV = $ValenceElectrons - $NumOfHydrogens;
  }

  return ($Delta, $DeltaV);
}

# Calculate DeltaIState matrix for atoms with each row i containing DeltaIij values
# corresponding atom atoms i and j.
#
# Notes:
#   . Matrix elements corresponding to hydrogen atoms or unconnected
#     are assigned zero value.
#
sub _CalculateDeltaIStateMatrix {
  my($This) = @_;
  my($DistanceMatrix, $NumOfRows, $NumOfCols, $RowIndex, $ColIndex, $AtomID1, $AtomID2, $DeltaIStateMatrix, $IStateValue1, $IStateValue2, $GraphDistance, $DeltaIState12, $DeltaIState21, $SkipIndexCheck);

  # Get distance matrix information...
  $DistanceMatrix = $This->{DistanceMatrix};
  ($NumOfRows, $NumOfCols) = $DistanceMatrix->GetSize();

  # Initialize DeltaIState matrix...
  $This->{DeltaIStateMatrix} = new Matrix($NumOfRows, $NumOfCols);
  $DeltaIStateMatrix = $This->{DeltaIStateMatrix};

  $SkipIndexCheck = 1;

  # Calculate DeltaIState matrix values...
  ROWINDEX: for $RowIndex (0 .. ($NumOfRows - 1) ) {
    $AtomID1 = $This->{AtomIndexToID}{$RowIndex};
    if (!exists($This->{IStateValues}{$AtomID1})) {
      next ROWINDEX;
    }
    $IStateValue1 = $This->{IStateValues}{$AtomID1};

    COLINDEX: for $ColIndex (($RowIndex + 1) .. ($NumOfCols - 1) ) {
      $AtomID2 = $This->{AtomIndexToID}{$ColIndex};
      if (!exists($This->{IStateValues}{$AtomID2})) {
	next COLINDEX;
      }
      $IStateValue2 = $This->{IStateValues}{$AtomID2};

      # Make sure it's a connected atom...
      $GraphDistance = $DistanceMatrix->GetValue($RowIndex, $ColIndex, $SkipIndexCheck);
      if ($GraphDistance >= BigNumber) {
	next COLINDEX;
      }

      $DeltaIState12 = ($IStateValue1 - $IStateValue2)/(($GraphDistance + 1)**2);
      $DeltaIState21 = -$DeltaIState12;

      # Set DeltaIState values...
      $DeltaIStateMatrix->SetValue($RowIndex, $ColIndex, $DeltaIState12, $SkipIndexCheck);
      $DeltaIStateMatrix->SetValue($ColIndex, $RowIndex, $DeltaIState21, $SkipIndexCheck);
    }
  }
}

# Setup distance matrix...
#
sub _SetupDistanceMatrix {
  my($This) = @_;

  $This->{DistanceMatrix} = $This->GetMolecule()->GetDistanceMatrix();

  if (!defined($This->{DistanceMatrix})) {
    return undef;
  }

  return $This;
}

# Setup final descriptor values...
#
sub _SetFinalValues {
  my($This) = @_;
  my($Atom, $AtomID, $EStateValue);

  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomID = $Atom->GetID();
    if (!exists $This->{EStateValues}{$AtomID}) {
      next ATOM;
    }
    $EStateValue = $This->{EStateValues}{$AtomID};
    $This->SetDescriptorValue($Atom, $EStateValue);
  }

  return $This;
}

# Cache  appropriate molecule data...
#
sub _SetupMoleculeDataCache {
  my($This) = @_;

  # Get all atoms including hydrogens to correctly map atom indicies to atom IDs for
  # usage of distance matrix.
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

# Return a string containg data for EStateValuesDescriptors object...
#
sub StringifyEStateValuesDescriptors {
  my($This) = @_;
  my($EStateValuesDescriptorsString);

  # Type of AtomicValues...
  $EStateValuesDescriptorsString = "AtomicDescriptorType: $This->{Type}; IgnoreHydrogens: " . ($This->{IgnoreHydrogens} ? "Yes" : "No");

  # Setup atomic descriptor information...
  my($AtomID, $DescriptorValue, @DescriptorValuesInfo, %DescriptorValues);

  @DescriptorValuesInfo = ();
  %DescriptorValues = $This->GetDescriptorValues();

  for $AtomID (sort { $a <=> $b } keys %DescriptorValues) {
    $DescriptorValue = $DescriptorValues{$AtomID};
    $DescriptorValue = (TextUtil::IsEmpty($DescriptorValue) || $DescriptorValue =~ /^None$/i) ? 'None' : MathUtil::round($DescriptorValue, 3) + 0;
    push @DescriptorValuesInfo, "$AtomID:$DescriptorValue";
  }
  $EStateValuesDescriptorsString .= "; AtomIDs:EStateValuesDescriptors: <" . TextUtil::JoinWords(\@DescriptorValuesInfo, ", ", 0) . ">";

  return $EStateValuesDescriptorsString;
}

# Is it a EStateValuesDescriptors object?
sub _IsEStateValuesDescriptors {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

EStateValuesDescriptors

=head1 SYNOPSIS

use AtomicDescriptors::EStateValuesDescriptors;

use AtomicDescriptors::EStateValuesDescriptors qw(:all);

=head1 DESCRIPTION

B<EStateValuesDescriptors> class provides the following methods:

new, GenerateDescriptors, StringifyEStateValuesDescriptors

B<EStateValuesDescriptors> is derived from B<AtomicValues> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<EStateValuesDescriptors>, B<AtomicValues> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

For calculation of electrotopological state (E-state) values for non-hydrogen atoms:

Let:

    N = Principal quantum number or period number corresponding to
        element symbol

    Sigma = Number of sigma electrons involves in bonds to hydrogen and
            non-hydrogen atoms attached to atom
          = Number of sigma bonds to hydrogen and non-hydrogen atoms
            attached to atom
    PI = Number of PI electrons involved in bonds to non-hydrogen atoms
         attached to atom
       = Number of PI bonds to non-hydrogen atoms attached to atom

    LP = Number of lone pair electrons on atom

    Zv = Number of electrons in valence shell of atom

    X = Number of non-hydrogen atom neighbors or heavy atoms attached
        to atom
    H = Number of implicit and explicit hydrogens for atom

    Delta = Number of sigma electrons involved to bonds to non-hydrogen
            atoms
    DeltaV = ValenceDelta = Number of valence shell electrons not involved
             in bonding to hydrogen atoms

    Ii = Intrinsic state value for atom i

    DeltaIi = Sum of perturbations to intrinsic state value Ii of atom i
              by all other atoms besides atom i

    DeltaIij = Perturbation to intrinsic state value Ii of atom i by atom j

    Dij = Graph/bond distance between atom i and j
    Rij = Dij + 1

    Si = E-state value for atom i

Then:

    Delta = Sigma - H = X

    DeltaV = Zv - H
           = Sigma + PI + LP - H

    Ii = ( ( ( 2 / N ) ** 2 ) * DeltaV + 1 ) / Delta

    DeltaIi = SUM ( (Ii - Ij) / (Rij ** 2) ) for j = 1 to num of atoms skipping atom i

    Si = Ii + DeltaIi

The current release of MayaChemTools doesn't support calculation of E-state
values [ Ref 75-78 ] for hydrogens.

=head2 METHODS

=over 4

=item B<new>

    $NewEStateValuesDescriptors = new AtomicDescriptors::
                                  EStateValuesDescriptors(%NamesAndValues);

Using specified I<EStateValuesDescriptors> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<EStateValuesDescriptors>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'EState'
    IgnoreHydrogens = 1

Examples:

    $EStateValuesDescriptors = new AtomicDescriptors::EStateValuesDescriptors(
                              'Molecule' => $Molecule,
                              'IgnoreHydrogens' => 1);

=item B<GenerateDescriptors>

    $EStateValuesDescriptors->GenerateDescriptors();

Calculates E-state atomic descriptors for all the atoms in a molecule and returns
I<EStateValuesDescriptors>.

=item B<StringifyEStateValuesDescriptors>

    $String = $EStateValuesDescriptors->StringifyEStateValuesDescriptors();

Returns a string containing information about I<EStateValuesDescriptors> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AtomicDescriptors.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
