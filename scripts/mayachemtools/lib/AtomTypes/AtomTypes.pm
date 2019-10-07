package AtomTypes::AtomTypes;
#
# File: AtomTypes.pm
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
use ObjectProperty;
use TextUtil ();

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(ObjectProperty Exporter);
@EXPORT = qw(LoadAtomTypesData);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializeAtomTypes();

  $This->_InitializeAtomTypesProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeAtomTypes {
  my($This) = @_;

  # Molecule object...
  $This->{Molecule} = '';

  # Type of AtomType...
  $This->{Type} = '';

  # By default, atom types are also assigned to hydrogens...
  $This->{IgnoreHydrogens} = 0;

}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}


# Initialize object properties....
sub _InitializeAtomTypesProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  return $This;
}

# Set molecule object and make sure it's not already set...
#
sub SetMolecule {
  my($This, $Molecule) = @_;

  if ($This->{Molecule}) {
    croak "Error: ${ClassName}->SetMolecule: Can't change molecule object:  It's already set...";
  }
  $This->{Molecule} = $Molecule;

  # Weaken the reference to disable increment of reference count...
  Scalar::Util::weaken($This->{Molecule});

  return $This;
}

# Set type and make sure it's not already set...
#
sub SetType {
  my($This, $Type) = @_;

  if ($This->{Type}) {
    croak "Error: ${ClassName}->SetType: Can't change AtomType type:  It's already set...";
  }
  $This->{Type} = $Type;

  return $This;
}

# Set specific atom type...
#
sub SetAtomType {
  my($This, $Atom, $AtomType) = @_;
  my($MethodName);

  # Assign AtomType to Atom...
  $MethodName = "Set" . $This->{Type} . "AtomType";
  $Atom->$MethodName($AtomType);

  return $This;
}

# Get specific atom type...
#
sub GetAtomType {
  my($This, $Atom) = @_;
  my($MethodName, $AtomType);

  $MethodName = "Get" . $This->{Type} . "AtomType";
  $AtomType = $Atom->$MethodName();

  return defined $AtomType ? $AtomType : 'None';
}

# Get atom types for all atoms as a hash with atom ID and atom types as
# key/value pairs.
#
# Notes:
#   . Irrespective of ignore hydrogens value, atom type for hydrogens are also
#     returned. Based on value of ignore hydrogens, atom type assignment methodology
#     might igonore hydrogens and value of None is returned for the hydrogens.
#
sub GetAtomTypes {
  my($This) = @_;
  my($Atom, $AtomID,  %AtomTypesMap);

  %AtomTypesMap = ();

  if (!$This->{Molecule}) {
    return %AtomTypesMap;
  }

  for $Atom ($This->{Molecule}->GetAtoms()) {
    $AtomID = $Atom->GetID();
    $AtomTypesMap{$AtomID} = $This->GetAtomType($Atom);
  }

  return %AtomTypesMap;
}

# Are all atoms types successfully assigned?
#
# Notes:
#   . Dynamic checking of atom types assignment for atoms eliminates the need
#     to check and synchronize valid atom types during SetAtomType.
#
sub IsAtomTypesAssignmentSuccessful {
  my($This) = @_;
  my($Atom, $AtomType);

  ATOM: for $Atom ($This->{Molecule}->GetAtoms()) {
    if ($Atom->IsHydrogen() && $This->{IgnoreHydrogens}) {
      next ATOM;
    }
    $AtomType = $This->GetAtomType($Atom);
    if ($AtomType =~ /^None$/i) {
      return 0;
    }
  }

  return 1;
}

# Load atom types data from the specified CSV atom type file into the specified
# hash reference.
#
# The lines starting with # are treated as comments and ignored. First line
# not starting with # must contain column labels and the number of columns in
# all other data rows must match the number of column labels.
#
# The first column is assumed to contain atom types; all other columns contain data
# as indicated in their column labels.
#
# In order to avoid dependence of data access on the specified column labels, the
# column data is loaded into hash with Column<Num> and AtomType as hash keys;
# however, the data for the first column which is treated as AtomTypes is also loaded
# into an array with AtomTypes as hash key. The format of the data structure loaded
# into a specified hash reference is:
#
# @{$AtomTypesDataMapRef->{AtomTypes}} - Array of all possible atom types for all atoms
# @{$AtomTypesDataMapRef->{NonHydrogenAtomTypes}} - Array of all possible atom types for non-hydrogen atoms
# @{$AtomTypesDataMapRef->{ColLabels}} - Array of column labels
# %{$AtomTypesDataMapRef->{DataCol<Num>}} - Hash keys pair: <DataCol<Num>, AtomType>
#
# Caveats:
#   . The column number start from 1.
#   . Column data for first column is not loaded into <Column<Num>, AtomType> hash keys pairs.
#
# AtomType file examples: SYBYLAtomTypes.csv, DREIDINGAtomTypes.csv,
# MMFF94AtomTypes.csv etc.
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub LoadAtomTypesData {
  my($FirstParameter, @OtherParamaters) = @_;
  my($AtomTypesDataFile, $AtomTypesDataMapRef, $InDelim, $Line, $NumOfCols, $ColIndex, $ColNum, $ColLabel, $ColValue, $AtomType, %AtomTypes, @LineWords, @ColLabels, @ColDataLabels);

  if (Scalar::Util::blessed($FirstParameter)) {
    ($AtomTypesDataFile, $AtomTypesDataMapRef) = @OtherParamaters;
  }
  else {
    ($AtomTypesDataFile, $AtomTypesDataMapRef) = @_;
  }

  $InDelim = "\,";
  open ATOMTYPESDATAFILE, "$AtomTypesDataFile" or croak "Couldn't open $AtomTypesDataFile: $! ...";

  # Skip lines up to column labels...
  LINE: while ($Line = TextUtil::GetTextLine(\*ATOMTYPESDATAFILE)) {
    if ($Line !~ /^#/) {
      last LINE;
    }
  }

  # Initialize data map...
  %{$AtomTypesDataMapRef} = ();
  @{$AtomTypesDataMapRef->{AtomTypes}} = ();
  @{$AtomTypesDataMapRef->{NonHydrogenAtomTypes}} = ();
  @{$AtomTypesDataMapRef->{ColLabels}} = ();

  %AtomTypes = ();

  # Process column labels...
  @ColLabels= quotewords($InDelim, 0, $Line);
  $NumOfCols = @ColLabels;
  push @{$AtomTypesDataMapRef->{ColLabels}}, @ColLabels;

  # Set up column data labels for storing the data...
  @ColDataLabels = ();
  for $ColNum (1 .. $NumOfCols) {
    $ColLabel = "DataCol${ColNum}";
    push @ColDataLabels, $ColLabel;
  }

  # Initialize column data hash...
  for $ColIndex (1 .. $#ColDataLabels) {
    $ColLabel = $ColDataLabels[$ColIndex];
    %{$AtomTypesDataMapRef->{$ColLabel}} = ();
  }

  # Process atom types data assuming first column to be atom type..
  LINE: while ($Line = TextUtil::GetTextLine(\*ATOMTYPESDATAFILE)) {
    if ($Line =~ /^#/) {
      next LINE;
    }
    @LineWords = quotewords($InDelim, 0, $Line);
    if (@LineWords != $NumOfCols) {
      croak "Error: The number of data fields, @LineWords, in $AtomTypesDataFile must be $NumOfCols.\nLine: $Line...";
    }
    $AtomType = $LineWords[0];
    if (exists $AtomTypes{$AtomType}) {
      carp "Warning: Ignoring data for atom type, $AtomType, in file $AtomTypesDataFile: It has already been loaded.\nLine: $Line....";
      next LINE;
    }

    $AtomTypes{$AtomType} = $AtomType;
    push @{$AtomTypesDataMapRef->{AtomTypes}}, $AtomType;

    # Is it a non-hydrogen atom type?
    if ($AtomType !~ /^H/i || $AtomType =~ /^(HAL|HET|HEV)$/i || $AtomType =~ /^(He4|Ho6|Hf3|Hg1)/) {
      # Non-hydrogen SYBYL atom types starting with H: HAL, HET, HEV
      # Non-hydrogen UFF atom types starting with H: He4+4, Ho6+3, Hf3+4, Hg1+2
      #
      push @{$AtomTypesDataMapRef->{NonHydrogenAtomTypes}}, $AtomType;
    }

    # Track column data values...
    for $ColIndex (1 .. $#LineWords) {
      $ColLabel = $ColDataLabels[$ColIndex]; $ColValue = $LineWords[$ColIndex];
      $AtomTypesDataMapRef->{$ColLabel}{$AtomType} = $ColValue;
    }
  }
  close ATOMTYPESDATAFILE;
}

1;

__END__

=head1 NAME

AtomTypes - AtomTypes class

=head1 SYNOPSIS

use AtomTypes::AtomTypes;

use AtomTypes::AtomTypes qw(:all);

=head1 DESCRIPTION

B<AtomTypes> base class used to derive all other atom types classes provides the following methods:

new, GetAtomType, GetAtomTypes, IsAtomTypesAssignmentSuccessful,
LoadAtomTypesData, SetAtomType, SetMolecule, SetType

B<AtomTypes> class is  derived from B<ObjectProperty> base class which provides methods not
explicitly defined in B<Fingerprints> or B<ObjectProperty> classes using Perl's AUTOLOAD functionality.
These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

=head2 METHODS

=over 4

=item B<new>

    $NewAtomTypes = new AtomTypes::AtomTypes(%NamesAndValues);

Using specified I<AtomTypes> property names and values hash, B<new> method creates a new object
and returns a reference to newly created B<AtomTypes> object. By default, following properties are
initialized:

    Molecule = '';
    Type = '';
    IgnoreHydrogens = 0;

=item B<GetAtomType>

    $AtomType = $AtomTypes->GetAtomType($Atom);

Returns B<AtomType> value string assigned to I<Atom> by I<AtomTypes> object.

=item B<GetAtomTypes>

    %AtomTypes = $AtomTypes->GetAtomTypes();

Returns atom types assigned to atoms by I<AtomTypes> object as a hash
with atom ID and atom types as key and value pairs.

=item B<IsAtomTypesAssignmentSuccessful>

    $Status = $AtomTypes->IsAtomTypesAssignmentSuccessful();

Returns 1 or 0 based on whether atom types assignment was successfully performed.
For a successful atom types assignment, all atoms must have an atom type other
than a string I<None>.

=item B<LoadAtomTypesData>

    $AtomTypes->LoadAtomTypesData($AtomTypesDataMapRef);
    AtomTypes::AtomTypes::LoadAtomTypesData($AtomTypesDataMapRef);

Loads atom types data from the specified CSV atom type file into the specified hash
reference.

The lines starting with # are treated as comments and ignored. First line not starting with
# must contain column labels and the number of columns in all other data rows must match
the number of column labels.

The first column is assumed to contain atom types; all other columns contain data
as indicated in their column labels.

In order to avoid dependence of data access on the specified column labels, the
column data is loaded into hash with I<DataColNum> and I<AtomType> as hash keys;
however, the data for the first column which is treated as AtomTypes is also loaded
into an array with AtomTypes as hash key. The format of the data structure loaded
into a specified hash reference is:

    @{$AtomTypesDataMapRef->{AtomTypes}} - Array of all possible atom
                                           types for all atoms
    @{$AtomTypesDataMapRef->{NonHydrogenAtomTypes}} - Array of all possible
                                           atom types for non-hydrogen atoms
    @{$AtomTypesDataMapRef->{ColLabels}} - Array of column labels
    %{$AtomTypesDataMapRef->{DataCol<ColNum>}} - Hash keys pair:
                                           <DataCol<ColNum>, AtomType>

I<ColNum> starts from 1. Column data for first column is not loaded into I<DataColNum>,
I<AtomType> hash keys pairs.

=item B<SetAtomType>

    $AtomTypes->SetAtomType($Atom, $AtomType);

Assigns specific I<AtomType> to I<Atom> and returns I<AtomTypes>.

=item B<SetMolecule>

    $AtomTypes->SetMolecule($Molecule);

Sets I<Molecule> object for I<AtomTypes> and retuens I<AtomTypes>.

=item B<SetType>

    $AtomTypes->SetType($Type);

Sets I<Type> for I<AtomTypes> object and retuens I<AtomTypes>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AtomicInvariantsAtomTypes.pm, DREIDINGAtomTypes.pm, EStateAtomTypes.pm,
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
