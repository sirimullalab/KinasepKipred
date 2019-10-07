package MoleculeFileIO;
#
# File: MoleculeFileIO.pm
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
use FileIO::SDFileIO;
use FileIO::MDLMolFileIO;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(IsSupportedMoleculeFileFormat);

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
  $This->_InitializeMoleculeFileIO();

  $This->_InitializeMoleculeFileIOProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeMoleculeFileIO {
  my($This) = @_;

  # Reference to specific FileIO object...
  $This->{FileIORef} = '';

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

}

# Initialize object properties......
#
sub _InitializeMoleculeFileIOProperties {
  my($This, %NamesAndValues) = @_;

  if (!exists $NamesAndValues{Name}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying file name...";
  }

  if (!exists $NamesAndValues{Mode}) {
    $NamesAndValues{Mode} = 'Read';
  }

  # Make sure its a supported format and intialize FileIO object reference...
  $This->_SetFileIORef(%NamesAndValues);

  return $This;
}

# Setup FileIO object reference...
sub _SetFileIORef {
  my($This, %NamesAndValues) = @_;
  my($Name, $Status, $Format, $IOPackageName);

  $Name = $NamesAndValues{Name};

  ($Status, $Format, $IOPackageName) = $This->IsSupportedMoleculeFileFormat($Name);
  if (!$Status) {
    croak "Error: ${ClassName}->New: Object can't be instantiated: File format, $Name, is not supported: Currently supported file formats are: SDF, MDLMol...";
  }

  $This->{FileIORef} = ${IOPackageName}->new(%NamesAndValues);

  return $This;
}

# Is it a supported file format?
#
# In scalar context only status is returned; otherwise, file format and file IO package name is also
# returned.
#
# Note:
#   . To support additional file formats, this is the only method which needs to be changed.
#
#   . Currently supported file formats are:
#
#      SDF         .sdf, .sd
#      MDLMol   .mol
#
sub IsSupportedMoleculeFileFormat {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $Name);

  if ((@_ == 2) && (_IsMoleculeFileIO($FirstParameter))) {
    ($This, $Name) = ($FirstParameter, $SecondParameter);
  }
  else {
    ($Name) = ($FirstParameter);
  }
  my($Status, $Format, $IOPackageName);

  $Status = 0; $Format = 'NotSupported'; $IOPackageName = 'Unknown';

  FORMAT: {
    if (FileIO::SDFileIO::IsSDFile($Name)) { $Status = 1; $Format = 'SDF'; $IOPackageName = 'FileIO::SDFileIO'; last FORMAT; }
    if (FileIO::MDLMolFileIO::IsMDLMolFile($Name)) { $Status = 1; $Format = 'MDLMol'; $IOPackageName = 'FileIO::MDLMolFileIO'; last FORMAT; }
    $Status = 0; $Format = 'NotSupported'; $IOPackageName = 'Unknown';
  }

  return wantarray ? ($Status, $Format, $IOPackageName) : $Status;
}

# Prohibit file ref change...
#
sub SetFileIORef {
  my($This, $Value) = @_;

  carp "Warning: ${ClassName}->SetFileIORef: Explicit setting of file ref is not supported...";

  return $This;
}

# Prohibit file name change...
#
sub SetName {
  my($This, $Name) = @_;

  carp "Warning: ${ClassName}->SetName: Explicit setting of file name is not supported: It must be set during object instantiation...";

  return $This;
}

# Prohibit file mode change...
#
sub SetMode {
  my($This, $Mode) = @_;

  carp "Warning: ${ClassName}->SetMode: Explicit setting of file mode is not supported: It must be set during object instantiation...";

  return $This;
}

# Open file in a specific mode; default mode is Read only.
# Supported mode values are: Read, Write, Append, <, >, >>, r, w, a
#
sub Open {
  my($This, $Mode) = @_;

  return $This->{FileIORef}->Open($Mode);
}

# close file...
sub Close {
  my($This) = @_;

  return $This->{FileIORef}->Close();
}

# Read molecule string from file and return a molecule object...
sub ReadMolecule {
  my($This) = @_;

  return $This->{FileIORef}->ReadMolecule();
}

# Retrieve molecule string from file...
sub ReadMoleculeString {
  my($This) = @_;

  return $This->{FileIORef}->ReadMoleculeString();
}

# Write molecule using molecule object...
sub WriteMolecule {
  my($This, $Molecule) = @_;

  return $This->{FileIORef}->WriteMolecule($Molecule);
}

# Is it a MoleculeFileIO object?
sub _IsMoleculeFileIO {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

MoleculeFileIO

=head1 SYNOPSIS

use MoleculeFileIO;

use MoleculeFileIO qw(:all);

=head1 DESCRIPTION

B<MoleculeFileIO> class provides the following methods:

new, Close, IsSupportedMoleculeFileFormat, Open, ReadMolecule,
ReadMoleculeString, WriteMolecule

The following methods can also be used as functions:

IsSupportedMoleculeFileFormat

=head2 METHODS

=over 4

=item B<new>

    $NewMoleculeFileIO = new MoleculeFileIO([%PropertyNameAndValues]);

Using specified I<MoleculeFileIO> property names and values hash, B<new> method
creates a new object and returns a reference to newly created B<MoleculeFileIO> object.
By default, following properties are initialized:

    Name = ""
    Mode = ""
    FileIORef = ""

Based on extension of specified file I<Name>, an input class is automatically associated to
provide molecule read and write methods.

Examples:

    $Name = "Water.mol";
    $Mode = "Read";
    $MoleculeFileIO = new MoleculeFileIO('Name' => $Name,
                                         'Mode' => $Mode);
    $MoleculeFileIO->Open();
    $Molecule = $MoleculeFileIO->ReadMolecule();
    $Molecule->DetectRings();
    print "$Molecule\n";
    $MoleculeFileIO->Close();

    $MoleculeFileIO = new MoleculeFileIO('Name' => 'Sample1.sdf',
                                         'Mode' => 'Read');
    $MoleculeFileIO->Open();
    while ($Molecule = $MoleculeFileIO1->ReadMolecule()) {
        $Molecule->DetectRings();
        print "$Molecule\n";

        $DataLabelsAndValuesRef =
          $Molecule->GetDataFieldLabelAndValues();
        for $DataLabel (sort keys %{$DataLabelsAndValuesRef} ) {
            $DataValue = $DataLabelsAndValuesRef->{$DataLabel};
            print "<DataLabel: $DataLabel; DataValue: $DataValue>; ";
        }
        print "\n";
    }
    $MoleculeFileIO->Close();

=item B<Close>

    $MoleculeFileIO->Close();

Closes an open file

=item B<IsSupportedMoleculeFileFormat>

    $Status = $MoleculeFileIO->IsSupportedMoleculeFileFormat($Name);
    $Status = MoleculeFileIO::IsSupportedMoleculeFileFormat($Name);
    ($Status, $FormatType, $IOClassName) =
       $MoleculeFileIO::IsSupportedMoleculeFileFormat($Name);

Returns 1 or 0 based on whether input file I<Name> format is supported. In list context,
value of supported format type and name of associated IO class is also returned.

File extension is used to determine file format. Currently, following file extensions are
supported:

    FileExts - FormatType - AssociatedIOClassName

    .mol - MDLMOL - MDLMolFileIO
    .sdf, .sd - SDF - SDFileIO

=item B<Open>

    $MoleculeFileIO->Open([$Mode]);

Opens a file in a specified I<Mode>. Default mode value: I<Read>. Supported mode
values:

    Read, Write, Append, <, >, >>, r, w, or a

=item B<ReadMolecule>

    $Molecule = $MoleculeFileIO->ReadMolecule();

Reads molecule data from the file and returns a I<Molecule> object.

=item B<ReadMoleculeString>

    $MoleculeString = $MoleculeFileIO->ReadMoleculeString();

Reads molecule data from a file and returns a I<Molecule> string.

=item B<WriteMolecule>

    $MoleculeFileIO->WriteMolecule();

Write molecule data to a file for a I<Molecule>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

FileIO.pm, MDLMolFileIO.pm, SDFileIO.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
