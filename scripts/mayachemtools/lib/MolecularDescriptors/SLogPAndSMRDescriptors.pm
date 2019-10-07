package MolecularDescriptors::SLogPAndSMRDescriptors;
#
# File: SLogPAndSMRDescriptors.pm
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
use MathUtil ();
use Atom;
use Molecule;
use MolecularDescriptors::MolecularDescriptors;
use AtomTypes::SLogPAtomTypes;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(MolecularDescriptors::MolecularDescriptors Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GetDescriptorNames);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, @DescriptorNames);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifySLogPAndSMRDescriptors';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeSLogPAndSMRDescriptors();

  $This->_InitializeSLogPAndSMRDescriptorsProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Descriptor names...
  @DescriptorNames = ('SLogP', 'SMR');
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
sub _InitializeSLogPAndSMRDescriptors {
  my($This) = @_;

  # Type of MolecularDescriptor...
  $This->{Type} = 'SLogPAndSMR';

  # SLogPAndSMR atom types assigned to hydrogen and non-hydrogen atoms...
  %{$This->{AtomTypes}} = ();

  # Intialize descriptor names and values...
  $This->_InitializeDescriptorNamesAndValues(@DescriptorNames);

  return $This;
}

# Initialize object properties...
#
sub _InitializeSLogPAndSMRDescriptorsProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  return $This;
}


# Calculate SLogPAndSMR value [ Ref 89 ] for a molecule...
#
# Methodology:
#   . Assign SLogP atom types to all atoms.
#   . Calculate SLogP and SMR value by adding contribution of each atom type.
#
# Caveats:
#   . All hydrogens must be added to molecule before calling GenerateDescriptors.
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

  # Cache appropriate molecule data...
  $This->_SetupMoleculeDataCache();

  # Assign SLogP atom types...
  if (!$This->_AssignAtomTypes()) {
    carp "Warning: ${ClassName}->GenerateDescriptors: $This->{Type} molecular descriptors generation didn't succeed: Couldn't assign valid SLogPAndSMR atom types to all atoms...";
    return undef;
  }

  # Calculate descriptor values...
  if (!$This->_CalculateDescriptorValues()) {
    carp "Warning: ${ClassName}->GenerateDescriptors: $This->{Type} molecular description generation didn't succeed: Couldn't calculate SLogPAndSMR values corresponding to assigned SLogP atom types...";
    return undef;
  }

  # Set final descriptor values...
  $This->_SetFinalDescriptorValues();

  # Clear cached molecule data...
  $This->_ClearMoleculeDataCache();

  return $This;
}

# Assign SLogPAndSMR atom types..
#
sub _AssignAtomTypes {
  my($This) = @_;
  my($SLogPAtomTypes, $Atom, $AtomID);

  %{$This->{AtomTypes}} = ();

  # Assign atom types...
  $SLogPAtomTypes = new AtomTypes::SLogPAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => 0);
  $SLogPAtomTypes->AssignAtomTypes();

  # Make sure SLogP atom types assignment is successful...
  if (!$SLogPAtomTypes->IsAtomTypesAssignmentSuccessful()) {
    return undef;
  }

  # Collect assigned atom types...
  for $Atom (@{$This->{Atoms}}) {
    $AtomID = $Atom->GetID();
    $This->{AtomTypes}{$AtomID} = $SLogPAtomTypes->GetAtomType($Atom);
  }

  return $This;
}

# Calculate SLogP and SMR values...
#
sub _CalculateDescriptorValues {
  my($This) = @_;
  my($SLogP, $SMR, $AtomID, $SLogPAndSMRDataRef, $AtomType);

  $SLogP = 0; $SMR = 0;

  # Get SLogP and SMR atom types data...
  $SLogPAndSMRDataRef = AtomTypes::SLogPAtomTypes::GetSLogPAtomTypesData();

  for $AtomID (keys %{$This->{AtomTypes}}) {
    $AtomType = $This->{AtomTypes}{$AtomID};

    # Makes sure data for SLogp and SMR contribution exists for atom type...
    if (!(exists($SLogPAndSMRDataRef->{DataCol4}{$AtomType}) && exists($SLogPAndSMRDataRef->{DataCol5}{$AtomType}))) {
      return undef;
    }

    # Data for SLogP contribution is in column number 4...
    $SLogP += $SLogPAndSMRDataRef->{DataCol4}{$AtomType};

    # Data for SMR contribution is in column number 5...
    $SMR += $SLogPAndSMRDataRef->{DataCol5}{$AtomType};
  }

  # Track the calculated values...
  $This->{SLogP} = MathUtil::round($SLogP, 2);
  $This->{SMR} = MathUtil::round($SMR, 2);

  return $This;
}

# Setup final descriptor values...
#
sub _SetFinalDescriptorValues {
  my($This) = @_;

  $This->{DescriptorsGenerated} = 1;

  $This->SetDescriptorValues($This->{SLogP}, $This->{SMR});

  return $This;
}

# Cache  appropriate molecule data...
#
sub _SetupMoleculeDataCache {
  my($This) = @_;

  @{$This->{Atoms}} = $This->GetMolecule()->GetAtoms();

  return $This;
}

# Clear cached molecule data...
#
sub _ClearMoleculeDataCache {
  my($This) = @_;

  @{$This->{Atoms}} = ();

  return $This;
}

# Return a string containg data for SLogPAndSMRDescriptors object...
#
sub StringifySLogPAndSMRDescriptors {
  my($This) = @_;
  my($SLogPAndSMRDescriptorsString);

  $SLogPAndSMRDescriptorsString = "MolecularDescriptorType: $This->{Type}; "  . $This->_StringifyDescriptorNamesAndValues();

  return $SLogPAndSMRDescriptorsString;
}

# Is it a SLogPAndSMRDescriptors object?
sub _IsSLogPAndSMRDescriptors {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

SLogPAndSMRDescriptors

=head1 SYNOPSIS

use MolecularDescriptors::SLogPAndSMRDescriptors;

use MolecularDescriptors::SLogPAndSMRDescriptors qw(:all);

=head1 DESCRIPTION

B<SLogPAndSMRDescriptors> class provides the following methods:

new, GenerateDescriptors, GetDescriptorNames,
StringifySLogPAndSMRDescriptors

B<SLogPAndSMRDescriptors> is derived from B<MolecularDescriptors> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<SLogPAndSMRDescriptors>, B<MolecularDescriptors> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

After SLogP atom types [ Ref 89 ] has been assigned to all atoms in a molecule using
AtomTypes::SLogPAndSMR.pm module,  SLogP (calculated logP) and SMR (calculated molar
refractivity) values are calculated by adding up LogP and MR contributions of each atom
type.

=head2 METHODS

=over 4

=item B<new>

    $NewSLogPAndSMRDescriptors = new MolecularDescriptors::
                                 SLogPAndSMRDescriptors(
                                 %NamesAndValues);

Using specified I<SLogPAndSMRDescriptors> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<SLogPAndSMRDescriptors>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'SLogPAndSMR'
    @DescriptorNames = ('SLogP', 'SMR')
    @DescriptorValues = ('None', 'None')

Examples:

    $SLogPAndSMRDescriptors = new MolecularDescriptors::
                              SLogPAndSMRDescriptors();

    $SLogPAndSMRDescriptors->SetMolecule($Molecule);
    $SLogPAndSMRDescriptors->GenerateDescriptors();
    print "SLogPAndSMRDescriptors: $SLogPAndSMRDescriptors\n";

=item B<GenerateDescriptors>

    $SLogPAndSMRDescriptors->GenerateDescriptors();

Calculate SLogP and SMR values for  a molecule and returns I<SLogPAndSMRDescriptors>.

=item B<GetDescriptorNames>

    @DescriptorNames = $SLogPAndSMRDescriptors->GetDescriptorNames();
    @DescriptorNames = MolecularDescriptors::SLogPAndSMRDescriptors::
                       GetDescriptorNames();

Returns all available descriptor names as an array.

=item B<StringifySLogPAndSMRDescriptors>

    $String = $SLogPAndSMRDescriptors->StringifySLogPAndSMRDescriptors();

Returns a string containing information about I<SLogPAndSMRDescriptors> object.

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
