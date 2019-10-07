package FileIO::SDFileIO;
#
# File: SDFileIO.pm
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
use FileUtil ();
use SDFileUtil ();
use FileIO::FileIO;
use FileIO::MDLMolFileIO;
use Molecule;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(FileIO::FileIO Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(IsSDFile);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeSDFileIO();

  $This->_InitializeSDFileIOProperties(%NamesAndValues);

  return $This;
}

# Initialize any local object data...
#
sub _InitializeSDFileIO {
  my($This) = @_;

  # Sorting of MDL data fields during output: Keep the initial order or write 'em out alphabetically...
  $This->{SortDataFieldsDuringOutput} = 'No';

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

}

# Initialize object values...
sub _InitializeSDFileIOProperties {
  my($This, %NamesAndValues) = @_;

  # All other property names and values along with all Set/Get<PropertyName> methods
  # are implemented on-demand using ObjectProperty class.

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  if (!exists $NamesAndValues{Name}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying file name...";
  }

  # Make sure it's a SD file...
  $Name = $NamesAndValues{Name};
  if (!$This->IsSDFile($Name)) {
    croak "Error: ${ClassName}->New: Object can't be instantiated: File, $Name, doesn't appear to be SDF format...";
  }

  return $This;
}

# Is it a SD file?
sub IsSDFile ($;$) {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $FileName, $Status);

  if ((@_ == 2) && (_IsSDFileIO($FirstParameter))) {
    ($This, $FileName) = ($FirstParameter, $SecondParameter);
  }
  else {
    $FileName = $FirstParameter;
  }

  # Check file extension...
  $Status = FileUtil::CheckFileType($FileName, "sd sdf");

  return $Status;
}

# Read molecule from file and return molecule object...
sub ReadMolecule {
  my($This) = @_;
  my($FileHandle);

  $FileHandle = $This->GetFileHandle();
  return $This->ParseMoleculeString(SDFileUtil::ReadCmpdString($FileHandle));
}

# Write compound data along with any data field label and values using Molecule object...
sub WriteMolecule {
  my($This, $Molecule) = @_;

  if (!(defined($Molecule) && $Molecule->IsMolecule())) {
    carp "Warning: ${ClassName}->WriteMolecule: No data written: Molecule object is not specified...";
    return $This;
  }
  my($FileHandle);
  $FileHandle = $This->GetFileHandle();

  print $FileHandle $This->GenerateMoleculeString($Molecule) . "\n";

  return $This;
}

# Retrieve molecule string...
sub ReadMoleculeString {
  my($This) = @_;
  my($FileHandle);

  $FileHandle = $This->GetFileHandle();
  return SDFileUtil::ReadCmpdString($FileHandle);
}

# Parse molecule string and return molecule object. ParseMoleculeString supports two invocation methods: class
# method or a package function.
#
sub ParseMoleculeString {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $MoleculeString);

  if ((@_ == 2) && (_IsSDFileIO($FirstParameter))) {
    ($This, $MoleculeString) = ($FirstParameter, $SecondParameter);
  }
  else {
    $MoleculeString = $FirstParameter;
    $This = undef;
  }
  if (!$MoleculeString) {
    return undef;
  }
  # Parse molecule data...
  my($Molecule);
  $Molecule = FileIO::MDLMolFileIO::ParseMoleculeString($MoleculeString);

  # Process data label/value pairs...
  my(@MoleculeLines, @DataLabels, %DataLabelsAndValues);

  %DataLabelsAndValues = ();
  @MoleculeLines = split /\n/, $MoleculeString;
  @DataLabels = SDFileUtil::GetCmpdDataHeaderLabels(\@MoleculeLines);
  %DataLabelsAndValues = SDFileUtil::GetCmpdDataHeaderLabelsAndValues(\@MoleculeLines);

  # Store reference to data labels to keep track of their initial order in SD file...
  $Molecule->SetDataFieldLabels(\@DataLabels);

  # Store reference to SD data label/value pairs hash as a generic property of molecule...
  $Molecule->SetDataFieldLabelAndValues(\%DataLabelsAndValues);

  return $Molecule;
}

# Generate molecule string using molecule object...
sub GenerateMoleculeString {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $Molecule);

  if ((@_ == 2) && (_IsSDFileIO($FirstParameter))) {
    ($This, $Molecule) = ($FirstParameter, $SecondParameter);
  }
  else {
    $Molecule = $FirstParameter;
    $This = undef;
  }
  if (!defined($Molecule)) {
    return undef;
  }
  # Generate CTAB data...
  my($CmpdString);
  $CmpdString = FileIO::MDLMolFileIO::GenerateMoleculeString($Molecule);

  # Generate any data field labels and values...
  my($DataFieldLabelsAndValuesString);

  $DataFieldLabelsAndValuesString = '';
  if ($Molecule->HasProperty('DataFieldLabels')) {
    my($DataFieldLabelsRef, $DataFieldLabelAndValuesRef, $SortDataFields);

    $SortDataFields = (exists($This->{SortDataFieldsDuringOutput}) && $This->{SortDataFieldsDuringOutput} =~ /^Yes$/i) ? 1 : 0;

    $DataFieldLabelsRef = $Molecule->GetDataFieldLabels();
    $DataFieldLabelAndValuesRef = $Molecule->GetDataFieldLabelAndValues();
    $DataFieldLabelsAndValuesString = join "\n", SDFileUtil::GenerateCmpdDataHeaderLabelsAndValuesLines($DataFieldLabelsRef, $DataFieldLabelAndValuesRef, $SortDataFields);
  }

  return "${CmpdString }\n${DataFieldLabelsAndValuesString}\n\$\$\$\$";
}


# Is it a SDFileIO object?
sub _IsSDFileIO {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

SDFileIO

=head1 SYNOPSIS

use FileIO::SDFileIO;

use FileIO::SDFileIO qw(:all);

=head1 DESCRIPTION

B<SDFIleIO> class provides the following methods:

new, GenerateMoleculeString, IsSDFile, ParseMoleculeString, ReadMolecule,
ReadMoleculeString, WriteMolecule

The following methods can also be used as functions:

GenerateMoleculeString, IsSDFile, ParseMoleculeString

Data specific to B<SDFileIO> class not directly used by B<Molecule>, B<Atom> and
B<Bond> objects - data label/value pairs, atom SteroParity and so on - is associated to
and retrieved from approptiate objects using following methods:

    SetMDL<PropertyName>
    GetMDL<PropertyName>.

SD data label and values are attached to B<Molecule> object as a refernece to a hash
using SetDataFieldLabelAndValues and can be retrieved using GetDataFieldLabelAndValues
method.

B<SDFileIO> class is derived from I<FileIO> class and uses its methods to support
generic file related functionality.

=head2 METHODS

=over 4

=item B<new>

    $NewSDFileIO = new FileIO::SDFileIO(%NamesAndValues);

Using specified I<SDFileIO> property names and values hash, B<new> method creates a new object
and returns a reference to newly created B<SDFileIO> object.

=item B<GenerateMoleculeString>

    $MoleculeString = $SDFileIO->GenerateMoleculeString($Molecule);
    $MoleculeString = FileIO::SDFileIO::GenerateMoleculeString($Molecule);

Returns a B<MoleculeString> in SD format corresponding to I<Molecule>.

=item B<IsSDFile>

    $Status = $SDFileIO->IsSDFile($FileName);
    $Status = FileIO::SDFileIO::IsSDFile($FileName);

Returns 1 or 0 based on whether I<FileName> is a SD file.

=item B<ParseMoleculeString>

    $Molecule = $SDFileIO->ParseMoleculeString($MoleculeString);
    $Molecule = FileIO::SDFileIO::ParseMoleculeString($MoleculeString);

Parses I<MoleculeString> and returns a B<Molecule> object. SD data field label and value pairs
are associated to B<Molecule> object as a reference to a hash using:

    $Molecule->SetDataFieldLabelAndValues(\%DataLabelsAndValues)

The reference to hash can be retrieved by:

    $DataLabelsAndValues = $Molecule->GetDataFieldLabelAndValues();
    for $DataLabel (sort keys %{$DataLabelsAndValues}) {
        $DataValue = $DataLabelsAndValues->{$DataLabel};
    }

=item B<ReadMolecule>

    $Molecule = $SDFileIO->ReadMolecule($FileHandle);

Reads data for the next compound in a file using already opened I<FileHandle>, creates,
and returns a B<Molecule> object.

=item B<ReadMoleculeString>

    $MoleculeString = $SDFileIO->ReadMoleculeString($FileHandle);

Reads data for the next compound in a file using already opened I<FileHandle> and
returns a B<MoleculeString> corresponding to compound structure and other associated
data.

=item B<WriteMolecule>

    $SDFileIO->WriteMolecule($Molecule);

Writes I<Molecule> data to a file in MDLMol format and returns B<SDFileIO>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

MoleculeFileIO.pm, MDLMolFileIO.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
