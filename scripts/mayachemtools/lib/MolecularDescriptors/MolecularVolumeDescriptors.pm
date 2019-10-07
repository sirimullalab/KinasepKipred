package MolecularDescriptors::MolecularVolumeDescriptors;
#
# File: MolecularVolumeDescriptors.pm
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
use AtomTypes::AtomTypes;
use MolecularDescriptors::MolecularDescriptors;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(MolecularDescriptors::MolecularDescriptors Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GetDescriptorNames GetVDWAtomRadiiAndVolumesData);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, @DescriptorNames, %VDWAtomRadiiAndVolumesDataMap);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyMolecularVolumeDescriptors';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeMolecularVolumeDescriptors();

  $This->_InitializeMolecularVolumeDescriptorsProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Descriptor names...
  @DescriptorNames = ('MolecularVolume');

  # Initialize the data hash. It'll be loaded on demand later...
  %VDWAtomRadiiAndVolumesDataMap = ();

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
sub _InitializeMolecularVolumeDescriptors {
  my($This) = @_;

  # Type of MolecularDescriptor...
  $This->{Type} = 'MolecularVolume';

  # Intialize descriptor names and values...
  $This->_InitializeDescriptorNamesAndValues(@DescriptorNames);

  return $This;
}

# Initialize object properties...
#
sub _InitializeMolecularVolumeDescriptorsProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  return $This;
}

# Get VDW atom data loaded from VDW atom radii and and volumes data file as
# a reference to hash with the following hash data format:
#
# @{$VDWAtomRadiiAndVolumesDataMap{AtomTypes}} - Array of all possible atom type symbols for all atoms
# @{$VDWAtomRadiiAndVolumesDataMap->{ColLabels}} - Array of column labels
# %{$VDWAtomRadiiAndVolumesDataMap->{DataCol<Num>}} - Hash keys pair: <DataCol<Num>, AtomType>
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetVDWAtomRadiiAndVolumesData {

  # Make sure data is loaded...
  _CheckAndLoadVDWAtomRadiiAndVolumesData();

  return \%VDWAtomRadiiAndVolumesDataMap;
}

# Calculate van der Waals molecular volume [ Ref 93 ] of a molecule using
# atomic and bonds contributions...
#
# van der Waals molecular volume (A**3/molecule) is defined as:
#
# vdwMolecularVolume = SumOfAtomicVDWVolumeContributions - 5.92 * NumOfBonds
#                      - 14.7 * NumOfAromaticRings - 3.8 * NumOfNonAromaticRings
#
# Methodology:
#   . Add up van der Waals atom volumne of all atoms
#   . Calculate molecular volume by correcting sum of atom volumes for num of
#     bonds and rings
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

  # Calculate descriptor values...
  if (!$This->_CalculateDescriptorValues()) {
    carp "Warning: ${ClassName}->GenerateDescriptors: $This->{Type} molecular descriptors generation didn't succeed: Couldn't calculate MolecularVolume values: van der Waals atom volume data is not available for all atoms...";
    return undef;
  }

  # Set final descriptor values...
  $This->_SetFinalDescriptorValues();

  return $This;
}

# Calculate MolecularVolume value...
#
sub _CalculateDescriptorValues {
  my($This) = @_;
  my($Atom, $AtomID, $AtomSymbol, $SumOfVDWAtomVolumes, $Molecule, $MolecularVolume, $NumOfBonds, $NumOfAromaticRings, $NumOfNonAromaticRings, $VDWAtomRadiiAndVolumesDataMapRef);

  $MolecularVolume = 0;

  $VDWAtomRadiiAndVolumesDataMapRef = $This->GetVDWAtomRadiiAndVolumesData();
  $Molecule = $This->{Molecule};

  # Calculate atom volumes contribution to molecular volume...
  #
  $SumOfVDWAtomVolumes = 0;

  ATOM: for $Atom ($Molecule->GetAtoms()) {
    $AtomID = $Atom->GetID();
    $AtomSymbol = $Atom->GetAtomSymbol();

    # Make sure van der Waals atom volume is available...
    if (!exists $VDWAtomRadiiAndVolumesDataMap{DataCol3}{$AtomSymbol}) {
      return undef;
    }
    $SumOfVDWAtomVolumes += $VDWAtomRadiiAndVolumesDataMapRef->{DataCol3}{$AtomSymbol};
  }

  $NumOfBonds = $Molecule->GetNumOfBonds();
  $NumOfAromaticRings = $Molecule->GetNumOfAromaticRings();
  $NumOfNonAromaticRings = $Molecule->GetNumOfRings() - $NumOfAromaticRings;

  # Apply correction for bonds and rings...
  $MolecularVolume = $SumOfVDWAtomVolumes - 5.92 * $NumOfBonds - 14.7 * $NumOfAromaticRings - 3.8 * $NumOfNonAromaticRings;

  # Track the calculated values...
  $This->{MolecularVolume} = MathUtil::round($MolecularVolume, 2);

  return $This;
}

# Setup final descriptor values...
#
sub _SetFinalDescriptorValues {
  my($This) = @_;

  $This->{DescriptorsGenerated} = 1;

  $This->SetDescriptorValues($This->{MolecularVolume});

  return $This;
}

# Return a string containg data for MolecularVolumeDescriptors object...
#
sub StringifyMolecularVolumeDescriptors {
  my($This) = @_;
  my($MolecularVolumeDescriptorsString);

  $MolecularVolumeDescriptorsString = "MolecularDescriptorType: $This->{Type}; " . $This->_StringifyDescriptorNamesAndValues();

  return $MolecularVolumeDescriptorsString;
}

# Is it a MolecularVolumeDescriptors object?
sub _IsMolecularVolumeDescriptors {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Check and load van der Waals atom radii and volumes data...
#
sub _CheckAndLoadVDWAtomRadiiAndVolumesData {

  # Is it already loaded?
  if (exists $VDWAtomRadiiAndVolumesDataMap{AtomTypes}) {
    return;
  }

  _LoadVDWAtomRadiiAndVolumesData();
}

# Initialize van der Waals atom radii and volumes data from the file...
#
# Format:
#
# "AtomTypeSymbol","VDWAtomRadius(A)","VDWAtomVolume(A**3)/molecule"
# "H","1.20","7.24"
# "He","1.40","11.49"
#
sub  _LoadVDWAtomRadiiAndVolumesData {
  my($VDWAtomDataFile, $MayaChemToolsLibDir);

  $MayaChemToolsLibDir = FileUtil::GetMayaChemToolsLibDirName();

  $VDWAtomDataFile =  "$MayaChemToolsLibDir" . "/data/VDWAtomRadiiAndVolumes.csv";
  if (! -e "$VDWAtomDataFile") {
    croak "Error: MayaChemTools package file, $VDWAtomDataFile, is missing: Possible installation problems...";
  }

  %VDWAtomRadiiAndVolumesDataMap = ();
  AtomTypes::AtomTypes::LoadAtomTypesData($VDWAtomDataFile, \%VDWAtomRadiiAndVolumesDataMap);
};

1;

__END__

=head1 NAME

MolecularVolumeDescriptors

=head1 SYNOPSIS

use MolecularDescriptors::MolecularVolumeDescriptors;

use MolecularDescriptors::MolecularVolumeDescriptors qw(:all);

=head1 DESCRIPTION

B<MolecularVolumeDescriptors> class provides the following methods:

new, GenerateDescriptors, GetDescriptorNames,
GetVDWAtomRadiiAndVolumesData, StringifyMolecularVolumeDescriptors

B<MolecularVolumeDescriptors> is derived from B<MolecularDescriptors> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<MolecularVolumeDescriptors>, B<MolecularDescriptors> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

van der Waals molecular volume [ Ref 93 ] (A**3/molecule) of a molecule is
calculated using atomic and bonds contributions along with adjustments for
aromatic and non-aromatic rings using the following equation:

    vdwMolecularVolume = SumOfAtomicVDWVolumeContributions
                         - 5.92 * NumOfBonds
                         - 14.7 * NumOfAromaticRings
                         - 3.8 * NumOfNonAromaticRings

van der Waals atomic volume for atoms is taken from data file VDWAtomRadiiAndVolumes.csv
distributed with MayaChemTools. It contains van der Waals atom radii and atom and volumes
data for 38 elements; Table 2 [ Ref 93 ] contains data for only 15 elements. After converting
valid van der Waals atom radius data from pm (picometer) to A (Angstrom) available under column
name VanderWaalsRadius in PeriodicTableElementsData.csv data file, van der Waals atom volume
is calculated using: 4/3 * PI * (Radius ** 3). For elements specified in Table 2 [ Ref 93 ] -
H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, Te, I - the van der Waals atom radii and calculated
atom volumes match the values in the table.

=head2 METHODS

=over 4

=item B<new>

    $NewMolecularVolumeDescriptors = new MolecularDescriptors::
                                     MolecularVolumeDescriptors(
                                     %NamesAndValues);

Using specified I<MolecularVolumeDescriptors> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<MolecularVolumeDescriptors>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'MolecularVolume'
    @DescriptorNames = ('MolecularVolume')
    @DescriptorValues = ('None')

Examples:

    $MolecularVolumeDescriptors = new MolecularDescriptors::
                                  MolecularVolumeDescriptors();

    $MolecularVolumeDescriptors->SetMolecule($Molecule);
    $MolecularVolumeDescriptors->GenerateDescriptors();
    print "MolecularVolumeDescriptors: $MolecularVolumeDescriptors\n";

=item B<GenerateDescriptors>

    $MolecularVolumeDescriptors->GenerateDescriptors();

Calculate van der Waals molecular volume descriptor for a molecule and returns
I<MolecularVolumeDescriptors>.

=item B<GetDescriptorNames>

    @DescriptorNames = $MolecularVolumeDescriptors->GetDescriptorNames();
    @DescriptorNames = MolecularDescriptors::MolecularVolumeDescriptors::
                          GetDescriptorNames();

Returns all available descriptor names as an array.

=item B<GetVDWAtomRadiiAndVolumesData>

    $VDWVolumeDataMapRef = $MolecularVolumeDescriptors->
                              GetVDWAtomRadiiAndVolumesData();
    $VDWVolumeDataMapRef = MolecularDescriptors::MolecularVolumeDescriptors::
                              GetVDWAtomRadiiAndVolumesData();

Returns a hash reference to van der Waals atom symbols corresponding to atom types
and associated data loaded from VDWAtomRadiiAndVolumes.csv data file as a reference
to hash with the following hash data format:

    @{$VDWVolumeDataMap{AtomTypes}} - Array of all possible atom
                types for all atom symbols
    @{$VDWVolumeDataMap->{ColLabels}} - Array of column labels
    %{$VDWVolumeDataMap->{DataCol<Num>}} - Hash keys pair:
                                                   DataCol<Num>, AtomType

=item B<StringifyMolecularVolumeDescriptors>

    $String = $MolecularVolumeDescriptors->
                              StringifyMolecularVolumeDescriptors();

Returns a string containing information about I<MolecularVolumeDescriptors> object.

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
