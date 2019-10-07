package MolecularDescriptors::RotatableBondsDescriptors;
#
# File: RotatableBondsDescriptors.pm
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
use TextUtil ();
use Atom;
use Molecule;
use MolecularDescriptors::MolecularDescriptors;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(MolecularDescriptors::MolecularDescriptors Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GetDescriptorNames);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, @DescriptorNames);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyRotatableBondsDescriptors';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeRotatableBondsDescriptors();

  $This->_InitializeRotatableBondsDescriptorsProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Descriptor names...
  @DescriptorNames = ('RotatableBonds');

}

# Get descriptor names as an array.
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetDescriptorNames {
  return @DescriptorNames;
}

# Initialize object data...
#
sub _InitializeRotatableBondsDescriptors {
  my($This) = @_;

  # Type of MolecularDescriptor...
  $This->{Type} = 'RotatableBonds';

  # MayaChemTools rotatable bonds default definition corresponds to modifed
  # version of rotatable bonds definition used by Veber et al. [ Ref 92 ]
  #
  $This->{IgnoreTerminalBonds} = 1;
  $This->{IgnoreBondsToTripleBonds} = 1;
  $This->{IgnoreAmideBonds} = 1;
  $This->{IgnoreThioamideBonds} = 1;
  $This->{IgnoreSulfonamideBonds} = 1;

  # Intialize descriptor names and values...
  $This->_InitializeDescriptorNamesAndValues(@DescriptorNames);

  return $This;
}

# Initialize object properties...
#
sub _InitializeRotatableBondsDescriptorsProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  return $This;
}

# Calculate number of rotatable bonds in a molecule...
#
# A rotatable bond is defined as any single bond which is not in a ring
# and involves only non-hydrogen atoms. By default, the following types
# of single bonds are not considered rotatable bonds:
#
#   . Terminal bonds
#   . Bonds attached to triple bonds
#   . Amide C-N bonds
#   . Thioamide C-N bond bonds
#   . Sulfonamide S-N bonds
#
# MayaChemTools rotatable bonds default definition corresponds to modifed
# version of rotatable bonds definition used by Veber et al. [ Ref 92 ]
#
sub GenerateDescriptors {
  my($This) = @_;

  # Initialize descriptor values...
  $This->_InitializeDescriptorValues();

  # Check availability of molecule...
  if (!$This->{Molecule}) {
    carp "Warning: ${ClassName}->GenerateDescriptors: $This->{Type} molecular descriptors generation didn't succeed: Molecule data is not available: Molecule object hasn't been set...";
    return undef;
  }

  # Calculate descriptor values...
  if (!$This->_CalculateDescriptorValues()) {
    carp "Warning: ${ClassName}->GenerateDescriptors: $This->{Type} molecular descriptors generation didn't succeed: Couldn't calculate RotatableBonds values...";
    return undef;
  }

  # Set final descriptor values...
  $This->_SetFinalDescriptorValues();

  return $This;
}

# Calculate RotatableBonds value...
#
sub _CalculateDescriptorValues {
  my($This) = @_;
  my($Bond, $RotatableBonds, $Atom1, $Atom2);

  $RotatableBonds = 0;

  BOND: for $Bond ($This->{Molecule}->GetBonds()) {
    # Is it a non-ring ring bond?
    if (!$Bond->IsSingle() || $Bond->IsInRing()) {
      next BOND;
    }

    ($Atom1, $Atom2) = $Bond->GetAtoms();

    # Does bond contain any Hydrogen atoms?
    if ($Atom1->IsHydrogen() || $Atom2->IsHydrogen()) {
      next BOND;
    }

    # Check for terminal bonds...
    if ($This->{IgnoreTerminalBonds} && $This->_IsTerminalBond($Atom1, $Atom2)) {
      next BOND;
    }

    # Check for bonds attached to triple bonds...
    if ($This->{IgnoreBondsToTripleBonds} && $This->_IsAttachedToTripleBond($Atom1, $Atom2)) {
      next BOND;
    }

    # Check for amide bonds...
    if ($This->{IgnoreAmideBonds} && $This->_IsAmideBond($Atom1, $Atom2)) {
      next BOND;
    }

    # Check for amide bonds...
    if ($This->{IgnoreThioamideBonds} && $This->_IsThioamideBond($Atom1, $Atom2)) {
      next BOND;
    }

    # Check for sulfonamide bonds...
    if ($This->{IgnoreSulfonamideBonds} && $This->_IsSulfonamideBond($Atom1, $Atom2)) {
      next BOND;
    }

    $RotatableBonds += 1;
  }

  # Track the calculated values...
  $This->{RotatableBonds} = $RotatableBonds;

  return $This;
}

# Is it a terminal bond?
#
sub _IsTerminalBond {
  my($This, $Atom1, $Atom2) = @_;

  return ($Atom1->GetAtomicInvariantValue('X') <= 1 || $Atom2->GetAtomicInvariantValue('X') <= 1 ) ? 1 : 0;
}

# Is it attached to a terminal bond?
#
sub _IsAttachedToTripleBond {
  my($This, $Atom1, $Atom2) = @_;

  return ($Atom1->GetAtomicInvariantValue('LBO') == 3 || $Atom2->GetAtomicInvariantValue('LBO') == 3) ? 1 : 0;
}

# Is it an amide bond?
#
# Amide: R-C(=O)-N(-R)(-R")
#
sub _IsAmideBond {
  my($This, $Atom1, $Atom2) = @_;
  my($CarbonAtom, $NitrogenAtom);

  ($CarbonAtom, $NitrogenAtom) = (undef, undef);

  if ($Atom1->IsCarbon() && $Atom2->IsNitrogen()) {
    ($CarbonAtom, $NitrogenAtom) = ($Atom1, $Atom2);
  }
  elsif ($Atom2->IsCarbon() && $Atom1->IsNitrogen()) {
    ($CarbonAtom, $NitrogenAtom) = ($Atom2, $Atom1);
  }

  if (!$CarbonAtom) {
    return 0;
  }

  return $CarbonAtom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['O', 'N', 'C,H'], ['=', '-', '-']) ? 1 : 0;
}

# Is it a thioamide bond?
#
# Thioamide: R-C(=S)-N(-R)(-R")
#
sub _IsThioamideBond {
  my($This, $Atom1, $Atom2) = @_;
  my($CarbonAtom, $NitrogenAtom);

  ($CarbonAtom, $NitrogenAtom) = (undef, undef);

  if ($Atom1->IsCarbon() && $Atom2->IsNitrogen()) {
    ($CarbonAtom, $NitrogenAtom) = ($Atom1, $Atom2);
  }
  elsif ($Atom2->IsCarbon() && $Atom1->IsNitrogen()) {
    ($CarbonAtom, $NitrogenAtom) = ($Atom2, $Atom1);
  }

  if (!$CarbonAtom) {
    return 0;
  }

  return $CarbonAtom->DoesAtomNeighborhoodMatch('C.T3.DB1', ['S', 'N', 'C,H'], ['=', '-', '-']) ? 1 : 0;
}

# Is it a sulfonamide bond?
#
# Sulfonamide: R-S(=O)(=O)-N(-R)(-R")
#
sub _IsSulfonamideBond {
  my($This, $Atom1, $Atom2) = @_;
  my($SulfurAtom, $NitrogenAtom);

  ($SulfurAtom, $NitrogenAtom) = (undef, undef);

  if ($Atom1->IsSulfur() && $Atom2->IsNitrogen()) {
    ($SulfurAtom, $NitrogenAtom) = ($Atom1, $Atom2);
  }
  elsif ($Atom2->IsSulfur() && $Atom1->IsNitrogen()) {
    ($SulfurAtom, $NitrogenAtom) = ($Atom2, $Atom1);
  }

  if (!$SulfurAtom) {
    return 0;
  }

  return $SulfurAtom->DoesAtomNeighborhoodMatch('S.T4.DB2', ['O', 'O', 'N', '!O'], ['=', '=', '-', '-']) ? 1 : 0;
}

# Setup final descriptor values...
#
sub _SetFinalDescriptorValues {
  my($This) = @_;

  $This->{DescriptorsGenerated} = 1;

  $This->SetDescriptorValues($This->{RotatableBonds});

  return $This;
}

# Return a string containg data for RotatableBondsDescriptors object...
#
sub StringifyRotatableBondsDescriptors {
  my($This) = @_;
  my($RotatableBondsDescriptorsString);

  # Type of MolecularDescriptors...
  $RotatableBondsDescriptorsString = "MolecularDescriptorType: $This->{Type}; IgnoreTerminalBonds: " . ($This->{IgnoreTerminalBonds} ? "Yes" : "No") . "; IgnoreBondsToTripleBonds: " .  ($This->{IgnoreBondsToTripleBonds} ? "Yes" : "No") . "; IgnoreAmideBonds: " .  ($This->{IgnoreAmideBonds} ? "Yes" : "No") . "; IgnoreThioamideBonds: " .  ($This->{IgnoreThioamideBonds} ? "Yes" : "No") . "; IgnoreSulfonamideBonds: " .  ($This->{IgnoreSulfonamideBonds} ? "Yes" : "No");

  # Setup molecular descriptor information...
  $RotatableBondsDescriptorsString .= "; " . $This->_StringifyDescriptorNamesAndValues();

  return $RotatableBondsDescriptorsString;
}

# Is it a RotatableBondsDescriptors object?
sub _IsRotatableBondsDescriptors {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

RotatableBondsDescriptors

=head1 SYNOPSIS

use MolecularDescriptors::RotatableBondsDescriptors;

use MolecularDescriptors::RotatableBondsDescriptors qw(:all);

=head1 DESCRIPTION

B<RotatableBondsDescriptors> class provides the following methods:

new, GenerateDescriptors, GetDescriptorNames,
StringifyRotatableBondsDescriptors

B<RotatableBondsDescriptors> is derived from B<MolecularDescriptors> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<RotatableBondsDescriptors>, B<MolecularDescriptors> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

A rotatable bond [ Ref 92 ] is defined as any single bond which is not in a ring and involves only non-hydrogen
atoms. By default, the following types of single bonds are not considered rotatable bonds:

    o Terminal bonds
    o Bonds attached to triple bonds
    o Amide C-N bonds
    o Thioamide C-N bond bonds
    o Sulfonamide S-N bonds

=head2 METHODS

=over 4

=item B<new>

    $RotatableBondsDescriptors = new MolecularDescriptors::
                                 RotatableBondsDescriptors(
                                 %NamesAndValues);

Using specified I<RotatableBondsDescriptors> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<RotatableBondsDescriptors>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'RotatableBonds'
    IgnoreTerminalBonds = 1
    IgnoreBondsToTripleBonds = 1
    IgnoreAmideBonds = 1
    IgnoreThioamideBonds = 1
    IgnoreSulfonamideBonds = 1
    @DescriptorNames = ('RotatableBonds')
    @DescriptorValues = ('None')

Examples:

    $RotatableBondsDescriptors = new MolecularDescriptors::
                                 RotatableBondsDescriptors();

    $RotatableBondsDescriptors = new MolecularDescriptors::
                                 RotatableBondsDescriptors(
                                 'IgnoreAmideBonds' => 0,
                                 'IgnoreThioamideBonds' => 0,
                                 'IgnoreSulfonamideBonds' => 0);

    $RotatableBondsDescriptors->SetMolecule($Molecule);
    $RotatableBondsDescriptors->GenerateDescriptors();
    print "RotatableBondsDescriptors: $RotatableBondsDescriptors\n";

=item B<GenerateDescriptors>

    $RotatableBondsDescriptors->GenerateDescriptors();

Calculates number of rotatable bonds descriptors in a molecule and returns
I<RotatableBondsDescriptors>.

=item B<GetDescriptorNames>

    @DescriptorNames = $RotatableBondsDescriptors->GetDescriptorNames();
    @DescriptorNames = MolecularDescriptors::RotatableBondsDescriptors::
                       GetDescriptorNames();

Returns all available descriptor names as an array.

=item B<StringifyRotatableBondsDescriptors>

    $String = $RotatableBondsDescriptors->StringifyRotatableBondsDescriptors();

Returns a string containing information about I<RotatableBondsDescriptors> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

MolecularDescriptors.pm, MolecularDescriptorsGenerator.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
