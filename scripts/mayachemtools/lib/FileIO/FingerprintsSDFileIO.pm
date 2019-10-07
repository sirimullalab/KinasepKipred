package FileIO::FingerprintsSDFileIO;
#
# File: FingerprintsSDFileIO.pm
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
use Fingerprints::FingerprintsStringUtil ();
use FileIO::FileIO;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(FileIO::FileIO Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(IsFingerprintsSDFile);

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
  $This->_InitializeFingerprintsSDFileIO();

  $This->_InitializeFingerprintsSDFileIOProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeFingerprintsSDFileIO {
  my($This) = @_;

  # Fingerprints string data format during read/write...
  #
  # For file read:
  #
  # AutoDetect  - automatically detect format of fingerprints string
  # FingerprintsBitVectorString - Bit vector fingerprints string format
  # FingerprintsVectorString - Vector fingerprints string format
  #
  # Default value: AutoDetect
  #
  # For file write:
  #
  # FingerprintsBitVectorString - Bit vector fingerprints string format
  # FingerprintsVectorString - Vector fingerprints string format
  #
  # Default value: undef
  #
  $This->{FingerprintsStringMode} = undef;

  # For file read:
  #
  #   o Fingerprints bit-vector and vector object for current fingerprints string
  #
  # For file write:
  #
  #   o Fingerprints bit-vector and vector object for current fingerprints string
  #   o Any supported fingerprints object: PathLengthFingerprints, ExtendedConnectivity, and so on.
  #
  $This->{FingerprintsObject} = undef;

  # Fingerprints SD file data field label during read/write
  #
  # For file read:
  #
  # Value of AutoDetect implies use first data field containing the word Fingerprints in its
  # data field label to retrieve fingerprints string data. Othwewise, a valid data field name
  # must be specified.
  #
  # For file write:
  #
  # Data field label to use for writing fingerprints string. Default: Fingerprints
  #
  $This->{FingerprintsFieldLabel} = undef;

  # Fingepritns string for current line during read/write...
  $This->{FingerprintsString} = undef;

  # First compound data string read/write...
  $This->{FirstCompoundDataIO} = 1;

  # Current fingerprints string data compound number during read/write...
  $This->{CompoundNum} = 0;

  # Compound data string during read/write...
  $This->{CompoundString} = undef;

  # Initialize parameters for read...
  $This->_InitializeFingerprintsSDFileIORead();

  # Initialize parameters for write...
  $This->_InitializeFingerprintsSDFileIOWrite();

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

}

# Initialize object data for reading fingerprints SD file...
#
sub _InitializeFingerprintsSDFileIORead {
  my($This) = @_;

  # Compound ID mode to use for retrieving compound IDs for fingerprints...
  #
  # Specify how to generate compound IDs: use a SD file datafield value; use molname line from
  # SD file; generate a sequential ID with specific prefix; use combination of both MolName and
  # LabelPrefix with usage of LabelPrefix values for empty molname lines.
  #
  # Possible values: DataField, MolName, LabelPrefix, or MolNameOrLabelPrefix. Default: LabelPrefix.
  #
  # For MolNameAndLabelPrefix value of CompoundIDMode, molname line in SD file takes precedence over
  # sequential compound IDs generated using LabelPrefix and only empty molname values are replaced
  # with sequential compound IDs.
  #
  $This->{CompoundIDMode} = 'LabelPrefix';

  #
  # Compound ID data field label name whose value is used as compound ID during DatafField value of
  # CompoundIDMode
  #
  $This->{CompoundIDFieldLabel} = undef;

  # A prefix string used for generating compound IDs like LabelPrefixString<Number> during LabelPrefix
  # or MolNameOrLabelPrefix value of  CompoundIDMode. Default value, Cmpd, generates compound IDs
  # which look like Cmpd<Number>.
  #
  $This->{CompoundIDPrefix} = 'Cmpd';

  # By default, the fingerprints data corresponding to FingerprintsCol is assumed to
  # be valid and no validation is performed before generating fingerprints objects...
  #
  $This->{ValidateData} = 1;

  # Level of detail to print during validation of data for invalid or missing data...
  $This->{DetailLevel} = 1;

  # Number of missing and invalid fingerprints string data compound strings...
  $This->{NumOfCmpdsWithMissingData} = 0;
  $This->{NumOfCmpdsWithInvalidData} = 0;

  # Compound ID for current fingerprints string...
  $This->{CompoundID} = undef;

  # Compound data field labels and values map for current compound data...
  %{$This->{DataFieldLabelsAndValues}} = ();

  # Status of data in fingerprints SD file...
  $This->{ValidFileData} = 0;

  $This->{ValidCompoundIDField} = 0;
  $This->{ValidFingerprintsField} = 0;

  $This->{ValidFingerprintsStringMode} = 0;

  return $This;
}

# Initialize object data for writing fingerprints SD file...
#
sub _InitializeFingerprintsSDFileIOWrite {
  my($This) = @_;

  # Fingerprints bit vector string format...
  #
  # Possible values: BinaryString or HexadecimalString [Default]
  #
  # Default BitStringFormat is set during first write using Fingerprints::FingerprintsStringUtil::GetDefaultBitStringFormat.
  #
  $This->{BitStringFormat} = undef;

  # Bits order in fingerprints bit vector string...
  #
  # Ascending - First bit in each byte as the lowest bit [Default]
  # Descending - First bit in each byte as the highest bit
  #
  # Default BitsOrder is set during first write using Fingerprints::FingerprintsStringUtil::GetDefaultBitsOrder.
  #
  $This->{BitsOrder} = undef;

  # Fingerprints vector string format...
  #
  # Possible values: IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString, ValuesAndIDsPairsString, ValuesString
  #
  # Default VectorStringFormat is set during first write using Fingerprints::FingerprintsStringUtil::GetDefaultVectorStringFormat.
  # For fingerprints vector object containing vector NumericalValues, it corresponds to IDsAndValuesString; othwerwise,
  # it's set to ValuesString.
  #
  $This->{VectorStringFormat} = undef;

  # Overwriting existing file...
  $This->{Overwrite} = 0;

  return $This;
}

# Initialize object values...
sub _InitializeFingerprintsSDFileIOProperties {
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

  # Make sure it's a fingerprints file...
  $Name = $NamesAndValues{Name};
  if (!$This->IsFingerprintsSDFile($Name)) {
    croak "Error: ${ClassName}->New: Object can't be instantiated: File, $Name, doesn't appear to be fingerprints format...";
  }

  if ($This->GetMode() =~ /^Read$/i) {
    $This->_InitializeFingerprintsSDFileIOReadProperties(%NamesAndValues);
  }
  elsif ($This->GetMode() =~ /^(Write|Append)$/i) {
    $This->_InitializeFingerprintsSDFileIOWriteProperties(%NamesAndValues);
  }

  return $This;
}

# Initialize object properties for reading fingerprints SD file...
#
sub _InitializeFingerprintsSDFileIOReadProperties {
  my($This, %NamesAndValues) = @_;

  # Set default value for FingerprintsStringMode...
  if (!$This->{FingerprintsStringMode}) {
    $This->{FingerprintsStringMode} = 'AutoDetect';
  }

  # Set default value for FingerprintsFieldLabel...
  if (!$This->{FingerprintsFieldLabel}) {
    $This->{FingerprintsFieldLabel} = 'AutoDetect';
  }

  # Check compound ID data field...
  if (($This->{CompoundIDMode} =~ /^DataField$/i) && (!defined($This->{CompoundIDFieldLabel}))) {
    croak "Error: ${ClassName}->: Object can't be instantiated: Compound ID data field lable must be specifed using \"CompoundIDFieldLabel\" during \"DataField\" value of \"CompoundIDMode\"...";
  }

  $This->_PrepareForReadingFingerprintsSDFileData();

  return $This;
}

# Initialize object properties for writing fingerprints SD file...
#
sub _InitializeFingerprintsSDFileIOWriteProperties {
  my($This, %NamesAndValues) = @_;

  # Check FingerprintsStringMode value...
  if (!exists $NamesAndValues{FingerprintsStringMode}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying FingerprintsStringMode...";
  }

  if ($This->{FingerprintsStringMode} !~ /^(FingerprintsBitVectorString|FingerprintsVectorString)$/i) {
    croak "Error: ${ClassName}->: Object can't be instantiated: FingerprintsStringMode value, $This->{FingerprintsStringMode}, is not valid; Supported values for write/append: FingerprintsBitVectorString or FingerprintsVectorString...";
  }

  # Set default value for FingerprintsFieldLabel...
  if (!$This->{FingerprintsFieldLabel}) {
    $This->{FingerprintsFieldLabel} = 'Fingerprints';
  }

  $This->_PrepareForWritingFingerprintsSDFileData();

  return $This;
}

# Set FingerprintsStringMode...
#
sub SetFingerprintsStringMode {
  my($This, $Value) = @_;

  # AutoDetect - automatically detect format of fingerprints string
  # FingerprintsBitVectorString - Bit vector fingerprints string format
  # FingerprintsVectorString - Vector fingerprints string format

  if ($Value !~ /^(AutoDetect|FingerprintsBitVectorString|FingerprintsVectorString)$/i) {
    croak "Error: ${ClassName}->SetFingerprintsStringMode: FingerprintsStringMode value, $Value, is not valid; Supported values: AutoDetect, FingerprintsBitVectorString or FingerprintsVectorString...";
  }

  $This->{FingerprintsStringMode} = $Value;

  return $This;
}

# Set CompoundIDMode...
#
sub SetCompoundIDMode {
  my($This, $Value) = @_;

  if ($Value !~ /^(DataField|MolName|LabelPrefix|MolNameOrLabelPrefix)$/i) {
    croak "Error: ${ClassName}->SetFingerprintsStringMode: CompoundIDMode value, $Value, is not valid; Supported values: DataField, MolName, LabelPrefix, or MolNameOrLabelPrefix...";
  }

  $This->{CompoundIDMode} = $Value;

  return $This;
}

# Set DetailLevel...
#
sub SetDetailLevel {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetDetailLevel: DetailLevel value, $Value, is not valid; Supported values: > 0...";
  }

  $This->{DetailLevel} = $Value;

  return $This;
}

# Set BitStringFormat...
#
sub SetBitStringFormat {
  my($This, $Value) = @_;

  if ($Value !~ /^(BinaryString|HexadecimalString)$/i) {
    croak "Error: ${ClassName}->SetBitStringFormat: BitStringFormat value, $Value, is not valid; Supported values: BinaryString or HexadecimalString...";
  }

  $This->{BitStringFormat} = $Value;

  return $This;
}

# Set BitsOrder...
#
sub SetBitsOrder {
  my($This, $Value) = @_;

  # Ascending - First bit in each byte as the lowest bit
  # Descending - First bit in each byte as the highest bit
  #
  if ($Value !~ /^(Ascending|Descending)$/i) {
    croak "Error: ${ClassName}->SetBitsOrder: FingerprintsStringMode value, $Value, is not valid; Supported values: Ascending or Descending...";
  }

  $This->{BitsOrder} = $Value;

  return $This;
}

# Set VectorStringFormat...
#
sub SetVectorStringFormat {
  my($This, $Value) = @_;

  # Possible values: IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString, ValuesAndIDsPairsString, ValuesString

  if ($Value !~ /^(IDsAndValuesString|IDsAndValuesPairsString|ValuesAndIDsString|ValuesAndIDsPairsString|ValuesString)$/i) {
    croak "Error: ${ClassName}->SetVectorStringFormat: FingerprintsStringMode value, $Value, is not valid; Supported values: IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString, ValuesAndIDsPairsString, or ValuesString...";
  }

  $This->{VectorStringFormat} = $Value;

  return $This;
}

# Get compound string for current compound with optional removal of fingerprints data..
#
sub GetCompoundString {
  my($This, $RemoveFingerprintsData) = @_;

  $RemoveFingerprintsData = defined $RemoveFingerprintsData ? $RemoveFingerprintsData : 0;

  if ($RemoveFingerprintsData && $This->_IsFingerprintsDataPresentInCompoundString()) {
    return SDFileUtil::RemoveCmpdDataHeaderLabelAndValue($This->{CompoundString}, $This->{FingerprintsFieldLabel});
  }

  return $This->{CompoundString};
}

# Set compound string for current compound..
#
sub SetCompoundString {
  my($This, $CompoundString) = @_;

  $This->{CompoundString} = $CompoundString;

  return $This;
}

# Get fingerprints object for current compound using fingerprints, fingerprints bit-vector
# fingerprints vector object. Fingerprints object correspond to any of supported fingerprints
# objects such as PathLengthFingerprints, ExtendedConnectivity, and so on.
#
sub GetFingerprints {
  my($This) = @_;

  return $This->{FingerprintsObject};
}

# Set fingerprints object for current compound...
#
sub SetFingerprints {
  my($This, $FingerprintsObject) = @_;

  $This->{FingerprintsObject} = $FingerprintsObject;

  return $This;
}

# Get fingerprints string  for current compound...
#
sub GetFingerprintsString {
  my($This) = @_;

  return $This->{FingerprintsString} ? $This->{FingerprintsString} : 'None';
}

# Set fingerprints string for current compound...
#
sub SetFingerprintsString {
  my($This, $FingerprintsString) = @_;

  $This->{FingerprintsString} = $FingerprintsString;

  return $This;
}

# Does fingerprints SD file contain valid data?
#
sub IsFingerprintsFileDataValid {
  my($This) = @_;

  return $This->{ValidFileData} ? 1 : 0;
}

# Does current compound contains valid fingerprints object data?
#
sub IsFingerprintsDataValid {
  my($This) = @_;

  return defined $This->{FingerprintsObject} ? 1 : 0;
}

# Read next available compound data string,  process it and generate appropriate fingerprints
# objects...
#
sub Read {
  my($This) = @_;

  # Read compound data string...
  if (!$This->_ReadCompoundDataString()) {
    return undef;
  }

  # No need to process invalid SD file with invalid data...
  if (!$This->{ValidFileData}) {
    if ($This->{ValidateData}) {
      $This->{NumOfCmpdsWithMissingData} += 1;
    }
    return $This;
  }

  # Perform data validation...
  if ($This->{ValidateData}) {
    if (!$This->_ValidateReadCompoundDataString()) {
      return $This;
    }
  }

  # Setup fingerprints string after checking again to handle problematic data for
  # non-validated compound string data...
  #
  my($FingerprintsFieldLabel);
  $FingerprintsFieldLabel = $This->{FingerprintsFieldLabel};
  if (exists $This->{DataFieldLabelsAndValues}{$FingerprintsFieldLabel}) {
    $This->{FingerprintsString} = $This->{DataFieldLabelsAndValues}{$FingerprintsFieldLabel};
  }

  # Generate fingeprints object...
  $This->_GenerateFingerprintsObject();

  # Setup fingerprints compound ID for fingerprints string...
  $This->_GenerateCompoundID();

  return $This;
}

# Read next available compound data string,  process it and generate appropriate fingerprints
# objects...
#
sub Next {
  my($This) = @_;

  return $This->Read();
}

# Read compound data string...
#
sub _ReadCompoundDataString {
  my($This) = @_;
  my(@CmpdLines);

  if ($This->{FirstCompoundDataIO}) {
    $This->_ProcessFirstCompoundDataStringRead();
  }

  # Initialize data for current compound data string...
  $This->_InitializeReadCompoundDataString();

  # Get next compound data line...
  $This->{CompoundString} = SDFileUtil::ReadCmpdString($This->{FileHandle});
  if (!$This->{CompoundString}) {
    return 0;
  }

  $This->{CompoundNum} += 1;

  # Set up data field labels and values...
  @CmpdLines = split "\n", $This->{CompoundString};
  %{$This->{DataFieldLabelsAndValues}} = SDFileUtil::GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

  return 1;
}

# Initialize compound data string for reading...
#
sub _InitializeReadCompoundDataString {
  my($This) = @_;

  $This->{CompoundID} = undef;
  $This->{CompoundString} = undef;

  %{$This->{DataFieldLabelsAndValues}} = ();

  $This->{FingerprintsObject} = undef;
  $This->{FingerprintsString} = undef;

  return $This;
}

# Validate compound data string containing fingerprints data...
#
sub _ValidateReadCompoundDataString {
  my($This) = @_;
  my($FingerprintsFieldLabel);

  $FingerprintsFieldLabel = $This->{FingerprintsFieldLabel};

  # Check for missing data...
  if (!exists $This->{DataFieldLabelsAndValues}{$FingerprintsFieldLabel}) {
    # Missing data...
    $This->{NumOfCmpdsWithMissingData} += 1;
    if ($This->{DetailLevel} >= 3) {
      carp "Warning: ${ClassName}->_ValidateReadCompoundDataString: Compound number $This->{CompoundNum} contains no fingerprints data: $This->{CompoundString}...";
    }
    elsif ($This->{DetailLevel} >= 2) {
      carp "Warning: ${ClassName}->_ValidateReadCompoundDataString: Compound number $This->{CompoundNum} contains no fingerprints data...";
    }
    return 0;
  }

  # Check for invalid data...
  my($InvalidFingerprintsData, $FingerprintsType, $FingerprintsDescription);

  $InvalidFingerprintsData = 0;

  if (Fingerprints::FingerprintsStringUtil::AreFingerprintsStringValuesValid($This->{DataFieldLabelsAndValues}{$FingerprintsFieldLabel})) {
    ($FingerprintsType, $FingerprintsDescription) = Fingerprints::FingerprintsStringUtil::GetFingerprintsStringTypeAndDescription($This->{DataFieldLabelsAndValues}{$FingerprintsFieldLabel});
    if (defined($FingerprintsType) && defined($FingerprintsDescription)) {
      if ($This->{FirstFingerprintsStringType} !~ /^$FingerprintsType$/i || $This->{FirstFingerprintsStringDescription} !~ /^$FingerprintsDescription$/i) {
	$InvalidFingerprintsData = 1;
      }
    }
    else {
      $InvalidFingerprintsData = 1;
    }
  }
  else {
    $InvalidFingerprintsData = 1;
  }

  if ($InvalidFingerprintsData) {
    $This->{NumOfCmpdsWithInvalidData} += 1;
    if ($This->{DetailLevel} >= 3) {
      carp "Warning: ${ClassName}->_ValidateReadCompoundDataString: Compound number $This->{CompoundNum} contains invalid fingerprints data: $This->{DataLine}...";
    }
    elsif ($This->{DetailLevel} >= 2) {
      carp "Warning: ${ClassName}->_ValidateReadCompoundDataString: Compound number $This->{CompoundNum} contains invalid fingerprints data...";
    }
    return 0;
  }

  return 1;
}

# Setup fingerprints compound ID for fingerprints string...
sub _GenerateCompoundID {
  my($This) = @_;
  my($CompoundID, $MolName);

  $CompoundID = '';

  if ($This->{CompoundIDMode} =~ /^LabelPrefix$/i) {
    $CompoundID = "$This->{CompoundIDPrefix}$This->{CompoundNum}";
  }
  elsif ($This->{CompoundIDMode} =~ /^DataField$/i) {
    my($SpecifiedDataField);
    $SpecifiedDataField = $This->{CompoundIDFieldLabel};
    $CompoundID = exists $This->{DataFieldLabelsAndValues}{$SpecifiedDataField} ?  $This->{DataFieldLabelsAndValues}{$SpecifiedDataField} : '';
  }
  elsif ($This->{CompoundIDMode} =~ /^MolNameOrLabelPrefix$/i) {
    ($MolName) = split "\n", $This->{CompoundString};
    $CompoundID = TextUtil::IsNotEmpty($MolName) ? $MolName : "$This->{CompoundIDPrefix}$This->{CompoundNum}";
  }
  elsif ($This->{CompoundIDMode} =~ /^MolName$/i) {
    ($MolName) = split "\n", $This->{CompoundString};
    $CompoundID = $MolName;
  }

  $This->{CompoundID} = $CompoundID;

  return $This;
}

# Process first compound data string read...
#
sub _ProcessFirstCompoundDataStringRead {
  my($This) = @_;
  my($Line, $FileHandle);

  $This->{FirstCompoundDataIO} = 0;

  return $This;
}

# Get ready for reading fingerprints SD file...
#
sub _PrepareForReadingFingerprintsSDFileData {
  my($This) = @_;

  # Retrieve SD file data fields information....
  $This->_RetrieveSDFileDataFields();

  # Validate compound and fingerprints field information...
  $This->_ValidateReadCompoundIDField();
  $This->_ValidateReadFingerprintsField();

  # Validate fingeprints string mode information...
  if ($This->{ValidFingerprintsField}) {
    $This->_ValidateReadFingerprintsStringMode();
  }

  # Set status of SD file data...
  $This->{ValidFileData} = ($This->{ValidCompoundIDField} && $This->{ValidFingerprintsField} && $This->{ValidFingerprintsStringMode}) ? 1 : 0;

  return $This;
}

# Retrieve information data fields and fingerprints string...
#
sub _RetrieveSDFileDataFields {
  my($This) = @_;
  my($SDFile, $CmpdString, @CmpdLines);

  $SDFile = $This->{Name};

  if (!(-e $SDFile)) {
    croak "Error: ${ClassName}->New: Object can't be instantiated: File, $SDFile, doesn't exist...";
  }

  if (!open SDFILE, "$SDFile") {
    croak "Error: ${ClassName}->New: Object can't be instantiated: Couldn't open input SD file $SDFile: $! ...";
  }
  $CmpdString = SDFileUtil::ReadCmpdString(\*SDFILE);
  close SDFILE;

  # Set up data field labels and values for first compound string data...
  @CmpdLines = split "\n", $CmpdString;

  %{$This->{FirstDataFieldLabelsAndValues}} = ();
  %{$This->{FirstDataFieldLabelsAndValues}} = SDFileUtil::GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

  return $This;
}

# Validate compound ID field information...
#
sub _ValidateReadCompoundIDField {
  my($This) = @_;
  my($SpecifiedDataField);

  $This->{ValidCompoundIDField} = 0;

  if ($This->{CompoundIDMode} =~ /^DataField$/i) {
    $SpecifiedDataField = $This->{CompoundIDFieldLabel};
    if (! exists $This->{FirstDataFieldLabelsAndValues}{$SpecifiedDataField}) {
      carp "Warning: ${ClassName}->_ValidateReadCompoundIDField: Compound ID data field, $SpecifiedDataField, specified using \"CompoundIDField\" in \"DataField\" \"CompoundIDMode\" doesn't exist...";
      return 0;
    }
  }

  $This->{ValidCompoundIDField} = 1;

  return 1;
}

# Validate fingerprints string field information...
#
sub _ValidateReadFingerprintsField {
  my($This) = @_;
  my($FingerprintsFieldLabel);

  $This->{ValidFingerprintsField} = 0;

  $FingerprintsFieldLabel = '';

  if ($This->{FingerprintsFieldLabel} !~ /^AutoDetect$/i) {
    $FingerprintsFieldLabel = $This->{FingerprintsFieldLabel};
    if (! exists $This->{FirstDataFieldLabelsAndValues}{$FingerprintsFieldLabel}) {
      carp "Warning: ${ClassName}->_ValidateReadFingerprintsField: Fingerprints data field value, $FingerprintsFieldLabel, specified using \"FingerprintsFieldLabel\" doesn't exist...";
      return 0;
    }
  }
  else {
    # Make sure default fingerprints field does exist...
    my($FingerprintsFieldFound, $DataFieldLabel);
    $FingerprintsFieldFound = 0;

    DATAFIELDLABEL: for $DataFieldLabel (keys %{$This->{FirstDataFieldLabelsAndValues}}) {
      if ($DataFieldLabel =~ /Fingerprints/i) {
	$FingerprintsFieldFound = 1;
	$FingerprintsFieldLabel = $DataFieldLabel;
	last DATAFIELDLABEL;
      }
    }
    if (!$FingerprintsFieldFound) {
      carp "Warning: ${ClassName}->_ValidateReadFingerprintsField: Data field label containing \"Fingerprints\" string in its name doesn't exist...";
      return 0;
    }
  }

  $This->{ValidFingerprintsField} = 1;
  $This->{FingerprintsFieldLabel} = $FingerprintsFieldLabel;

  return 1;
}

# Validate fingerprints string mode information...
#
sub _ValidateReadFingerprintsStringMode {
  my($This) = @_;
  my($FingerprintsBitVectorStringMode, $FingerprintsVectorStringMode, $FirstFingerprintsStringType, $FirstFingerprintsStringDescription, $FingerprintsFieldLabel, $FingerprintsType, $FingerprintsDescription);

  $This->{ValidFingerprintsStringMode} = 0;

  $This->{FingerprintsBitVectorStringMode} = 0;
  $This->{FingerprintsVectorStringMode} = 0;

  $This->{FirstFingerprintsStringType} = '';
  $This->{FirstFingerprintsStringDescription} = '';

  $FingerprintsBitVectorStringMode = 0;
  $FingerprintsVectorStringMode = 0;

  $FirstFingerprintsStringType = '';
  $FirstFingerprintsStringDescription = '';

  $FingerprintsFieldLabel = $This->{FingerprintsFieldLabel};

  ($FingerprintsType, $FingerprintsDescription) = Fingerprints::FingerprintsStringUtil::GetFingerprintsStringTypeAndDescription($This->{FirstDataFieldLabelsAndValues}{$FingerprintsFieldLabel});

  if ($This->{FingerprintsStringMode} =~ /^FingerprintsBitVectorString$/i) {
    if ($FingerprintsType !~ /^FingerprintsBitVector$/i) {
      carp "Warning: ${ClassName}->_ValidateReadFingerprintsStringMode: First fingerprint string data type, $FingerprintsType, doesn't match value, FingerprintsBitVectorString, specified using \"FingerprintsStringMode\"...";
      return 0;
    }
    $FingerprintsBitVectorStringMode = 1;
    $FirstFingerprintsStringType = 'FingerprintsBitVector';
    $FirstFingerprintsStringDescription = $FingerprintsDescription;
  }
  elsif ($This->{FingerprintsStringMode} =~ /^FingerprintsVectorString$/i) {
    if ($FingerprintsType !~ /^FingerprintsVector$/i) {
      carp "Warning: ${ClassName}->_ValidateReadFingerprintsStringMode: First fingerprint string data type, $FingerprintsType, doesn't match value, FingerprintsVectorString, specified using \"FingerprintsStringMode\"...";
      return 0;
    }
    $FingerprintsVectorStringMode = 1;
    $FirstFingerprintsStringType = 'FingerprintsVector';
    $FirstFingerprintsStringDescription = $FingerprintsDescription;
  }
  else {
    # AutoDetect mode...
    if ($FingerprintsType =~ /^FingerprintsBitVector$/i) {
      $FingerprintsBitVectorStringMode = 1;
    }
    elsif ($FingerprintsType =~ /^FingerprintsVector$/i) {
      $FingerprintsVectorStringMode = 1;
    }
    else {
      carp "Warning: ${ClassName}->_ValidateReadFingerprintsStringMode: First fingerprint string data type, $FingerprintsType, identified during, AutoDetect, value of \"FingerprintsStringMode\" is not valid; Supported fingerprints types: FingerprintBitVector or FingerprintsVector...";
      return 0;
    }
    $FirstFingerprintsStringType = $FingerprintsType;
    $FirstFingerprintsStringDescription = $FingerprintsDescription;
  }

  $This->{ValidFingerprintsStringMode} = 1;

  $This->{FingerprintsBitVectorStringMode} = $FingerprintsBitVectorStringMode;
  $This->{FingerprintsVectorStringMode} = $FingerprintsVectorStringMode;

  $This->{FirstFingerprintsStringType} = $FirstFingerprintsStringType;
  $This->{FirstFingerprintsStringDescription} = $FirstFingerprintsStringDescription;

  return 1;
}

# Write fingerprints string generated from specified fingerprints, fingerprints-bit vector, or
# fingerprints vector object and other data to SD file...
#
sub WriteFingerprints {
  my($This, $FingerprintsObject, $CompoundString) = @_;

  # Initialize data for current compound...
  $This->_InitializeWriteCompoundDataString();

  # Set fingerprints object...
  $This->{FingerprintsObject} = $FingerprintsObject;

  # Generate fingerprints string...
  $This->_GenerateFingerprintsString();

  # Set and update compound string...
  $This->{CompoundString} = $CompoundString;
  $This->_AddFingerprintsDataToCompoundString();

  # Write it out...
  $This->_WriteCompoundDataString();

  return $This;
}

# Write fingerprints string and other data to SD file...
#
# Note:
#   o FingerprintsStringMode, BitStringFormat, BitsOrder, VectorStringFormat values
#     are ignored during writing of fingerprints and it's written to the file as it is.
#   o CompoundString contains no fingerprints data
#
sub WriteFingerprintsString {
  my($This, $FingerprintsString, $CompoundString) = @_;

  # Initialize data for current compound...
  $This->_InitializeWriteCompoundDataString();

  # Set fingerprints string...
  $This->{FingerprintsString} = $FingerprintsString;

  # Generate fingerprints object...
  $This->_GenerateFingerprintsObject();

  # Set and update compound string...
  $This->{CompoundString} = $CompoundString;
  $This->_AddFingerprintsDataToCompoundString();

  # Write it out...
  $This->_WriteCompoundDataString();

  return $This;
}

# Initialize compound data string for writing...
#
sub _InitializeWriteCompoundDataString {
  my($This) = @_;

  $This->{CompoundString} = undef;

  $This->{FingerprintsObject} = undef;
  $This->{FingerprintsString} = undef;

  return $This;
}

# Writi compound data string...
#
sub _WriteCompoundDataString {
  my($This) = @_;
  my($FileHandle);

  if ($This->{FirstCompoundDataIO}) {
    $This->_ProcessFirstCompoundDataStringWrite();
  }

  $This->{CompoundNum} += 1;
  $FileHandle = $This->{FileHandle};

  print $FileHandle "$This->{CompoundString}\n";

  return $This;
}

# Process first compound data string write...
#
sub _ProcessFirstCompoundDataStringWrite {
  my($This) = @_;
  my($Line, $FileHandle);

  $This->{FirstCompoundDataIO} = 0;

  return $This;
}

# Get ready for writing fingerprints SD file...
#
sub _PrepareForWritingFingerprintsSDFileData {
  my($This) = @_;
  my($SDFile);

  $SDFile = $This->{Name};
  if (!$This->{Overwrite}) {
    if (-e $SDFile) {
      croak "Error: ${ClassName}->New: Object can't be instantiated: File, $SDFile, already exist. Use overwrite option...";
    }
  }

  # Setup FingerprintsStringMode status...

  $This->{FingerprintsBitVectorStringMode} = 0;
  $This->{FingerprintsVectorStringMode} = 0;
  $This->{ValidFingerprintsStringMode} = 0;

  if ($This->{FingerprintsStringMode} =~ /^FingerprintsBitVectorString$/i) {
    $This->{FingerprintsBitVectorStringMode} = 1;
  }
  elsif ($This->{FingerprintsStringMode} =~ /^FingerprintsVectorString$/i) {
    $This->{FingerprintsVectorStringMode} = 1;
  }

  $This->{ValidFingerprintsStringMode} = ($This->{FingerprintsBitVectorStringMode} || $This->{FingerprintsVectorStringMode}) ? 1 : 0;

  if ($This->{FingerprintsBitVectorStringMode}) {
    $This->_SetDefaultBitStringFormat();
    $This->_SetDefaultBitsOrder();
  }
  elsif ($This->{FingerprintsVectorStringMode}) {
    $This->_SetDefaultVectorStringFormat();
  }

  return $This;
}

# Set default value for bit string format...
#
sub _SetDefaultBitStringFormat {
  my($This) = @_;

  if (!$This->{BitStringFormat}) {
    $This->{BitStringFormat} = Fingerprints::FingerprintsStringUtil::GetDefaultBitStringFormat();
  }

  return $This;
}

# Set default value for bit string format...
#
sub _SetDefaultBitsOrder {
  my($This) = @_;

  if (!$This->{BitsOrder}) {
    $This->{BitsOrder} = Fingerprints::FingerprintsStringUtil::GetDefaultBitsOrder();
  }

  return $This;
}

# Set default value for vector string format...
#
sub _SetDefaultVectorStringFormat {
  my($This) = @_;

  if (!$This->{VectorStringFormat} && $This->{FingerprintsObject}) {
    $This->{VectorStringFormat} = Fingerprints::FingerprintsStringUtil::GetDefaultVectorStringFormat($This->{FingerprintsObject});
  }

  return $This;
}

# Add fingerprints data to compound string...
#
sub _AddFingerprintsDataToCompoundString {
  my($This) = @_;
  my($CmpdString);

  # Check and remove existing fingerprints data...
  if ($This->_IsFingerprintsDataPresentInCompoundString()) {
    carp "Warning: ${ClassName}->_AddFingerprintsDataToCompoundString: The compound string already contains fingerprints data corresponding to fingerprints field label $This->{FingerprintsFieldLabel}; It has been replaced with new fingerprints data...";
    $This->{CompoundString} = SDFileUtil::RemoveCmpdDataHeaderLabelAndValue($This->{CompoundString}, $This->{FingerprintsFieldLabel});
  }

  $CmpdString = $This->{CompoundString};

  $CmpdString =~ s/\$\$\$\$$//;

  $This->{CompoundString} = "${CmpdString}>  <$This->{FingerprintsFieldLabel}>\n$This->{FingerprintsString}\n\n\$\$\$\$";

  return $This;
}

# Is fingerprints data already present in compound string?
#
sub _IsFingerprintsDataPresentInCompoundString {
  my($This) = @_;
  my($FingerprintsFieldLabel);

  if (TextUtil::IsEmpty($This->{CompoundString}) || TextUtil::IsEmpty($This->{FingerprintsFieldLabel})) {
    return 0;
  }

  $FingerprintsFieldLabel = $This->{FingerprintsFieldLabel};

  return ($This->{CompoundString} =~ /<$FingerprintsFieldLabel>/) ? 1 : 0;
}

# Generate fingerprints object using current fingerprints string...
#
sub _GenerateFingerprintsObject {
  my($This) = @_;

  $This->{FingerprintsObject} = undef;

  if (!$This->{FingerprintsString}) {
    return $This;
  }

  if ($This->{FingerprintsBitVectorStringMode}) {
    $This->{FingerprintsObject} = Fingerprints::FingerprintsStringUtil::ParseFingerprintsBitVectorString($This->{FingerprintsString});
  }
  elsif ($This->{FingerprintsVectorStringMode}) {
    $This->{FingerprintsObject} = Fingerprints::FingerprintsStringUtil::ParseFingerprintsVectorString($This->{FingerprintsString});
  }
  else {
    return undef;
  }

  return $This;
}

# Generate fingerprints string using current fingerprints object...
#
sub _GenerateFingerprintsString {
  my($This) = @_;

  $This->{FingerprintsString} = '';

  if (!$This->{FingerprintsObject}) {
    return $This;
  }

  if ($This->{FingerprintsBitVectorStringMode}) {
    $This->{FingerprintsString} = Fingerprints::FingerprintsStringUtil::GenerateFingerprintsString($This->{FingerprintsObject}, $This->{BitStringFormat}, $This->{BitsOrder});
  }
  elsif ($This->{FingerprintsVectorStringMode}) {
    $This->{FingerprintsString} = Fingerprints::FingerprintsStringUtil::GenerateFingerprintsString($This->{FingerprintsObject}, $This->{VectorStringFormat});
  }

  return $This;
}

# Is it a fingerprints file?
sub IsFingerprintsSDFile ($;$) {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $FileName, $Status);

  if ((@_ == 2) && (_IsFingerprintsSDFileIO($FirstParameter))) {
    ($This, $FileName) = ($FirstParameter, $SecondParameter);
  }
  else {
    $FileName = $FirstParameter;
  }

  # Check file extension...
  $Status = FileUtil::CheckFileType($FileName, "sdf sd");

  return $Status;
}

# Is it a FingerprintsSDFileIO object?
sub _IsFingerprintsSDFileIO {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

FingerprintsSDFileIO

=head1 SYNOPSIS

use FileIO::FingerprintsSDFileIO;

use FileIO::FingerprintsSDFileIO qw(:all);

=head1 DESCRIPTION

B<FingerprintsSDFileIO> class provides the following methods:

new, GetCompoundString, GetFingerprints, GetFingerprintsString,
IsFingerprintsDataValid, IsFingerprintsFileDataValid, IsFingerprintsSDFile, Next,
Read, SetBitStringFormat, SetBitsOrder, SetCompoundIDMode, SetCompoundString,
SetDetailLevel, SetFingerprints, SetFingerprintsString, SetFingerprintsStringMode,
SetVectorStringFormat, WriteFingerprints, WriteFingerprintsString

The following methods can also be used as functions:

IsFingerprintsSDFile

B<FingerprintsSDFileIO> class is derived from I<FileIO> class and uses its methods to support
generic file related functionality.

The fingerprints SD file format with B<.sdf> or B<.sd> file extensions supports two types of
fingerprints string data: fingerprints bit-vectors and fingerprints vector strings. The fingerprints
string data is treated as value of a fingerprints data field label in a SD file.

Example of SD file format containing fingerprints string data:

    ... ...
    ... ...
    $$$$
    ... ...
    ... ...
    ... ...
    41 44  0  0  0  0  0  0  0  0999 V2000
     -3.3652    1.4499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    ... ...
    2  3  1  0  0  0  0
    ... ...
    M  END
    >  <CmpdID>
    Test

    >  <PathLengthFingerprints>
    FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes:MinLengt
    h1:MaxLength8;1024;HexadecimalString;Ascending;9c8460989ec8a49913991a66
    03130b0a19e8051c89184414953800cc2151082844a201042800130860308e8204d4028
    00831048940e44281c00060449a5000ac80c894114e006321264401600846c050164462
    08190410805000304a10205b0100e04c0038ba0fad0209c0ca8b1200012268b61c0026a
    aa0660a11014a011d46

    $$$$
    ... ...
    ... ...

The current release of MayaChemTools supports the following types of fingerprint
bit-vector and vector strings:

    FingerprintsVector;AtomNeighborhoods:AtomicInvariantsAtomTypes:MinRadi
    us0:MaxRadius2;41;AlphaNumericalValues;ValuesString;NR0-C.X1.BO1.H3-AT
    C1:NR1-C.X3.BO3.H1-ATC1:NR2-C.X1.BO1.H3-ATC1:NR2-C.X3.BO4-ATC1 NR0-C.X
    1.BO1.H3-ATC1:NR1-C.X3.BO3.H1-ATC1:NR2-C.X1.BO1.H3-ATC1:NR2-C.X3.BO4-A
    TC1 NR0-C.X2.BO2.H2-ATC1:NR1-C.X2.BO2.H2-ATC1:NR1-C.X3.BO3.H1-ATC1:NR2
    -C.X2.BO2.H2-ATC1:NR2-N.X3.BO3-ATC1:NR2-O.X1.BO1.H1-ATC1 NR0-C.X2.B...

    FingerprintsVector;AtomTypesCount:AtomicInvariantsAtomTypes:ArbitraryS
    ize;10;NumericalValues;IDsAndValuesString;C.X1.BO1.H3 C.X2.BO2.H2 C.X2
    .BO3.H1 C.X3.BO3.H1 C.X3.BO4 F.X1.BO1 N.X2.BO2.H1 N.X3.BO3 O.X1.BO1.H1
    O.X1.BO2;2 4 14 3 10 1 1 1 3 2

    FingerprintsVector;AtomTypesCount:SLogPAtomTypes:ArbitrarySize;16;Nume
    ricalValues;IDsAndValuesString;C1 C10 C11 C14 C18 C20 C21 C22 C5 CS F
    N11 N4 O10 O2 O9;5 1 1 1 14 4 2 1 2 2 1 1 1 1 3 1

    FingerprintsVector;AtomTypesCount:SLogPAtomTypes:FixedSize;67;OrderedN
    umericalValues;IDsAndValuesString;C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C
    12 C13 C14 C15 C16 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C27 CS N1 N
    2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12 N13 N14 NS O1 O2 O3 O4 O5 O6 O7 O8
    O9 O10 O11 O12 OS F Cl Br I Hal P S1 S2 S3 Me1 Me2;5 0 0 0 2 0 0 0 0 1
    1 0 0 1 0 0 0 14 0 4 2 1 0 0 0 0 0 2 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0...

    FingerprintsVector;EStateIndicies:ArbitrarySize;11;NumericalValues;IDs
    AndValuesString;SaaCH SaasC SaasN SdO SdssC SsCH3 SsF SsOH SssCH2 SssN
    H SsssCH;24.778 4.387 1.993 25.023 -1.435 3.975 14.006 29.759 -0.073 3
    .024 -2.270

    FingerprintsVector;EStateIndicies:FixedSize;87;OrderedNumericalValues;
    ValuesString;0 0 0 0 0 0 0 3.975 0 -0.073 0 0 24.778 -2.270 0 0 -1.435
    4.387 0 0 0 0 0 0 3.024 0 0 0 0 0 0 0 1.993 0 29.759 25.023 0 0 0 0 1
    4.006 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0

    FingerprintsVector;ExtendedConnectivity:AtomicInvariantsAtomTypes:Radi
    us2;60;AlphaNumericalValues;ValuesString;73555770 333564680 352413391
    666191900 1001270906 1371674323 1481469939 1977749791 2006158649 21414
    08799 49532520 64643108 79385615 96062769 273726379 564565671 85514103
    5 906706094 988546669 1018231313 1032696425 1197507444 1331250018 1338
    532734 1455473691 1607485225 1609687129 1631614296 1670251330 17303...

    FingerprintsVector;ExtendedConnectivityCount:AtomicInvariantsAtomTypes
    :Radius2;60;NumericalValues;IDsAndValuesString;73555770 333564680 3524
    13391 666191900 1001270906 1371674323 1481469939 1977749791 2006158649
    2141408799 49532520 64643108 79385615 96062769 273726379 564565671...;
    3 2 1 1 14 1 2 10 4 3 1 1 1 1 2 1 2 1 1 1 2 3 1 1 2 1 3 3 8 2 2 2 6 2
    1 2 1 1 2 1 1 1 2 1 1 2 1 2 1 1 1 1 1 1 1 1 1 2 1 1

    FingerprintsBitVector;ExtendedConnectivityBits:AtomicInvariantsAtomTyp
    es:Radius2;1024;BinaryString;Ascending;0000000000000000000000000000100
    0000000001010000000110000011000000000000100000000000000000000000100001
    1000000110000000000000000000000000010011000000000000000000000000010000
    0000000000000000000000000010000000000000000001000000000000000000000000
    0000000000010000100001000000000000101000000000000000100000000000000...

    FingerprintsVector;ExtendedConnectivity:FunctionalClassAtomTypes:Radiu
    s2;57;AlphaNumericalValues;ValuesString;24769214 508787397 850393286 8
    62102353 981185303 1231636850 1649386610 1941540674 263599683 32920567
    1 571109041 639579325 683993318 723853089 810600886 885767127 90326012
    7 958841485 981022393 1126908698 1152248391 1317567065 1421489994 1455
    632544 1557272891 1826413669 1983319256 2015750777 2029559552 20404...

    FingerprintsVector;ExtendedConnectivity:EStateAtomTypes:Radius2;62;Alp
    haNumericalValues;ValuesString;25189973 528584866 662581668 671034184
    926543080 1347067490 1738510057 1759600920 2034425745 2097234755 21450
    44754 96779665 180364292 341712110 345278822 386540408 387387308 50430
    1706 617094135 771528807 957666640 997798220 1158349170 1291258082 134
    1138533 1395329837 1420277211 1479584608 1486476397 1487556246 1566...

    FingerprintsBitVector;MACCSKeyBits;166;BinaryString;Ascending;00000000
    0000000000000000000000000000000001001000010010000000010010000000011100
    0100101010111100011011000100110110000011011110100110111111111111011111
    11111111111110111000

    FingerprintsBitVector;MACCSKeyBits;322;BinaryString;Ascending;11101011
    1110011111100101111111000111101100110000000000000011100010000000000000
    0000000000000000000000000000000000000000000000101000000000000000000000
    0000000000000000000000000000000000000000000000000000000000000000000000
    0000000000000000000000000000000000000011000000000000000000000000000000
    0000000000000000000000000000000000000000

    FingerprintsVector;MACCSKeyCount;166;OrderedNumericalValues;ValuesStri
    ng;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 1 0 0 3 0 0 0 0 4 0 0 2 0 0 0 0 0 0 0 0 2 0 0 2 0 0 0 0
    0 0 0 0 1 1 8 0 0 0 1 0 0 1 0 1 0 1 0 3 1 3 1 0 0 0 1 2 0 11 1 0 0 0
    5 0 0 1 2 0 1 1 0 0 0 0 0 1 1 0 1 1 1 1 0 4 0 0 1 1 0 4 6 1 1 1 2 1 1
    3 5 2 2 0 5 3 5 1 1 2 5 1 2 1 2 4 8 3 5 5 2 2 0 3 5 4 1

    FingerprintsVector;MACCSKeyCount;322;OrderedNumericalValues;ValuesStri
    ng;14 8 2 0 2 0 4 4 2 1 4 0 0 2 5 10 5 2 1 0 0 2 0 5 13 3 28 5 5 3 0 0
    0 4 2 1 1 0 1 1 0 0 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 22 5 3 0 0 0 1 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 2 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...

    FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes:MinLeng
    th1:MaxLength8;1024;BinaryString;Ascending;001000010011010101011000110
    0100010101011000101001011100110001000010001001101000001001001001001000
    0010110100000111001001000001001010100100100000000011000000101001011100
    0010000001000101010100000100111100110111011011011000000010110111001101
    0101100011000000010001000011000010100011101100001000001000100000000...

    FingerprintsVector;PathLengthCount:AtomicInvariantsAtomTypes:MinLength
    1:MaxLength8;432;NumericalValues;IDsAndValuesPairsString;C.X1.BO1.H3 2
    C.X2.BO2.H2 4 C.X2.BO3.H1 14 C.X3.BO3.H1 3 C.X3.BO4 10 F.X1.BO1 1 N.X
    2.BO2.H1 1 N.X3.BO3 1 O.X1.BO1.H1 3 O.X1.BO2 2 C.X1.BO1.H3C.X3.BO3.H1
    2 C.X2.BO2.H2C.X2.BO2.H2 1 C.X2.BO2.H2C.X3.BO3.H1 4 C.X2.BO2.H2C.X3.BO
    4 1 C.X2.BO2.H2N.X3.BO3 1 C.X2.BO3.H1:C.X2.BO3.H1 10 C.X2.BO3.H1:C....

    FingerprintsVector;PathLengthCount:MMFF94AtomTypes:MinLength1:MaxLengt
    h8;463;NumericalValues;IDsAndValuesPairsString;C5A 2 C5B 2 C=ON 1 CB 1
    8 COO 1 CR 9 F 1 N5 1 NC=O 1 O=CN 1 O=CO 1 OC=O 1 OR 2 C5A:C5B 2 C5A:N
    5 2 C5ACB 1 C5ACR 1 C5B:C5B 1 C5BC=ON 1 C5BCB 1 C=ON=O=CN 1 C=ONNC=O 1
    CB:CB 18 CBF 1 CBNC=O 1 COO=O=CO 1 COOCR 1 COOOC=O 1 CRCR 7 CRN5 1 CR
    OR 2 C5A:C5B:C5B 2 C5A:C5BC=ON 1 C5A:C5BCB 1 C5A:N5:C5A 1 C5A:N5CR ...

    FingerprintsVector;TopologicalAtomPairs:AtomicInvariantsAtomTypes:MinD
    istance1:MaxDistance10;223;NumericalValues;IDsAndValuesString;C.X1.BO1
    .H3-D1-C.X3.BO3.H1 C.X2.BO2.H2-D1-C.X2.BO2.H2 C.X2.BO2.H2-D1-C.X3.BO3.
    H1 C.X2.BO2.H2-D1-C.X3.BO4 C.X2.BO2.H2-D1-N.X3.BO3 C.X2.BO3.H1-D1-...;
    2 1 4 1 1 10 8 1 2 6 1 2 2 1 2 1 2 2 1 2 1 5 1 10 12 2 2 1 2 1 9 1 3 1
    1 1 2 2 1 3 6 1 6 14 2 2 2 3 1 3 1 8 2 2 1 3 2 6 1 2 2 5 1 3 1 23 1...

    FingerprintsVector;TopologicalAtomPairs:FunctionalClassAtomTypes:MinDi
    stance1:MaxDistance10;144;NumericalValues;IDsAndValuesString;Ar-D1-Ar
    Ar-D1-Ar.HBA Ar-D1-HBD Ar-D1-Hal Ar-D1-None Ar.HBA-D1-None HBA-D1-NI H
    BA-D1-None HBA.HBD-D1-NI HBA.HBD-D1-None HBD-D1-None NI-D1-None No...;
    23 2 1 1 2 1 1 1 1 2 1 1 7 28 3 1 3 2 8 2 1 1 1 5 1 5 24 3 3 4 2 13 4
    1 1 4 1 5 22 4 4 3 1 19 1 1 1 1 1 2 2 3 1 1 8 25 4 5 2 3 1 26 1 4 1 ...

    FingerprintsVector;TopologicalAtomTorsions:AtomicInvariantsAtomTypes;3
    3;NumericalValues;IDsAndValuesString;C.X1.BO1.H3-C.X3.BO3.H1-C.X3.BO4-
    C.X3.BO4 C.X1.BO1.H3-C.X3.BO3.H1-C.X3.BO4-N.X3.BO3 C.X2.BO2.H2-C.X2.BO
    2.H2-C.X3.BO3.H1-C.X2.BO2.H2 C.X2.BO2.H2-C.X2.BO2.H2-C.X3.BO3.H1-O...;
    2 2 1 1 2 2 1 1 3 4 4 8 4 2 2 6 2 2 1 2 1 1 2 1 1 2 6 2 4 2 1 3 1

    FingerprintsVector;TopologicalAtomTorsions:EStateAtomTypes;36;Numerica
    lValues;IDsAndValuesString;aaCH-aaCH-aaCH-aaCH aaCH-aaCH-aaCH-aasC aaC
    H-aaCH-aasC-aaCH aaCH-aaCH-aasC-aasC aaCH-aaCH-aasC-sF aaCH-aaCH-aasC-
    ssNH aaCH-aasC-aasC-aasC aaCH-aasC-aasC-aasN aaCH-aasC-ssNH-dssC a...;
    4 4 8 4 2 2 6 2 2 2 4 3 2 1 3 3 2 2 2 1 2 1 1 1 2 1 1 1 1 1 1 1 2 1 1 2

    FingerprintsVector;TopologicalAtomTriplets:AtomicInvariantsAtomTypes:M
    inDistance1:MaxDistance10;3096;NumericalValues;IDsAndValuesString;C.X1
    .BO1.H3-D1-C.X1.BO1.H3-D1-C.X3.BO3.H1-D2 C.X1.BO1.H3-D1-C.X2.BO2.H2-D1
    0-C.X3.BO4-D9 C.X1.BO1.H3-D1-C.X2.BO2.H2-D3-N.X3.BO3-D4 C.X1.BO1.H3-D1
    -C.X2.BO2.H2-D4-C.X2.BO2.H2-D5 C.X1.BO1.H3-D1-C.X2.BO2.H2-D6-C.X3....;
    1 2 2 2 2 2 2 2 8 8 4 8 4 4 2 2 2 2 4 2 2 2 4 2 2 2 2 1 2 2 4 4 4 2 2
    2 4 4 4 8 4 4 2 4 4 4 2 4 4 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 8...

    FingerprintsVector;TopologicalAtomTriplets:SYBYLAtomTypes:MinDistance1
    :MaxDistance10;2332;NumericalValues;IDsAndValuesString;C.2-D1-C.2-D9-C
    .3-D10 C.2-D1-C.2-D9-C.ar-D10 C.2-D1-C.3-D1-C.3-D2 C.2-D1-C.3-D10-C.3-
    D9 C.2-D1-C.3-D2-C.3-D3 C.2-D1-C.3-D2-C.ar-D3 C.2-D1-C.3-D3-C.3-D4 C.2
    -D1-C.3-D3-N.ar-D4 C.2-D1-C.3-D3-O.3-D2 C.2-D1-C.3-D4-C.3-D5 C.2-D1-C.
    3-D5-C.3-D6 C.2-D1-C.3-D5-O.3-D4 C.2-D1-C.3-D6-C.3-D7 C.2-D1-C.3-D7...

    FingerprintsVector;TopologicalPharmacophoreAtomPairs:ArbitrarySize:Min
    Distance1:MaxDistance10;54;NumericalValues;IDsAndValuesString;H-D1-H H
    -D1-NI HBA-D1-NI HBD-D1-NI H-D2-H H-D2-HBA H-D2-HBD HBA-D2-HBA HBA-D2-
    HBD H-D3-H H-D3-HBA H-D3-HBD H-D3-NI HBA-D3-NI HBD-D3-NI H-D4-H H-D4-H
    BA H-D4-HBD HBA-D4-HBA HBA-D4-HBD HBD-D4-HBD H-D5-H H-D5-HBA H-D5-...;
    18 1 2 1 22 12 8 1 2 18 6 3 1 1 1 22 13 6 5 7 2 28 9 5 1 1 1 36 16 10
    3 4 1 37 10 8 1 35 10 9 3 3 1 28 7 7 4 18 16 12 5 1 2 1

    FingerprintsVector;TopologicalPharmacophoreAtomPairs:FixedSize:MinDist
    ance1:MaxDistance10;150;OrderedNumericalValues;ValuesString;18 0 0 1 0
    0 0 2 0 0 1 0 0 0 0 22 12 8 0 0 1 2 0 0 0 0 0 0 0 0 18 6 3 1 0 0 0 1
    0 0 1 0 0 0 0 22 13 6 0 0 5 7 0 0 2 0 0 0 0 0 28 9 5 1 0 0 0 1 0 0 1 0
    0 0 0 36 16 10 0 0 3 4 0 0 1 0 0 0 0 0 37 10 8 0 0 0 0 1 0 0 0 0 0 0
    0 35 10 9 0 0 3 3 0 0 1 0 0 0 0 0 28 7 7 4 0 0 0 0 0 0 0 0 0 0 0 18...

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

=head2 METHODS

=over 4

=item B<new>

    $NewFingerprintsSDFileIO = new FileIO::FingerprintsSDFileIO(%IOParameters);

Using specified I<IOParameters> names and values hash, B<new> method creates a new
object and returns a reference to a newly created B<FingerprintsSDFileIO> object. By default,
the following properties are initialized during I<Read> mode:

    Name = '';
    Mode = 'Read';
    Status = 0;
    FingerprintsStringMode = 'AutoDetect';
    FingerprintsFieldLabel = 'AutoDetect';
    CompoundIDMode = 'LabelPrefix';
    CompoundIDFieldLabel = undef;
    CompoundIDPrefix = 'Cmpd';
    ValidateData = 1;
    DetailLevel = 1;

During I<Write> mode, the following properties get initialize by default:

    FingerprintsStringMode = undef;

    BitStringFormat = HexadecimalString;
    BitsOrder = Ascending;

    VectorStringFormat = NumericalValuesString or ValuesString;

Examples:

    $NewFingerprintsSDFileIO = new FileIO::FingerprintsSDFileIO(
                               'Name' => 'Sample.sdf',
                               'Mode' => 'Read');

    $NewFingerprintsSDFileIO = new FileIO::FingerprintsSDFileIO(
                               'Name' => 'Sample.sdf',
                               'Mode' => 'Read',;
                               'FingerprintsStringMode' =>
                                       'AutoDetect',
                               'FingerprintsFieldLabel' =>
                                       'Fingerprints',
                               'CompoundIDMode' =>
                                       'DataField',
                               'CompoundIDFieldLabel' =>
                                       'CompoundID');

    $NewFingerprintsSDFileIO = new FileIO::FingerprintsSDFileIO(
                               'Name' => 'Sample.sdf',
                               'Mode' => 'Write',
                               'FingerprintsStringMode' =>
                                       'FingerprintsBitVectorString',
                               'Overwrite' => 1,
                               'BitStringFormat' => 'HexadecimalString',
                               'BitsOrder' => 'Ascending');

    $NewFingerprintsSDFileIO = new FileIO::FingerprintsSDFileIO(
                               'Name' => 'Sample.sd',
                               'Mode' => 'Write',
                               'FingerprintsStringMode' =>
                                       'FingerprintsVectorString',
                               'Overwrite' => 1,
                               'VectorStringFormat' => 'IDsAndValuesString',
                               'FingerprintsLabel' => 'Fingerprints');

=item B<GetCompoundString>

    $CompoundString = $FingerprintsSDFileIO->GetCompoundString();

Returns B<CompoundString> for current compound.

=item B<GetFingerprints>

    $FingerprintsObject = $FingerprintsSDFileIO->GetFingerprints();

Returns B<FingerprintsObject> generated for current compound using fingerprints bit-vector
or vector string data. The fingerprints object corresponds to any of the supported fingerprints
such as PathLengthFingerprints, ExtendedConnectivity, and so on.

=item B<GetFingerprintsString>

    $FingerprintsString = $FingerprintsSDFileIO->GetFingerprintsString();

Returns B<FingerprintsString> for current compound.

=item B<IsFingerprintsDataValid>

    $Status = $FingerprintsSDFileIO->IsFingerprintsDataValid();

Returns 1 or 0 based on whether B<FingerprintsObject> is valid.

=item B<IsFingerprintsFileDataValid>

    $Status = $FingerprintsSDFileIO->IsFingerprintsFileDataValid();

Returns 1 or 0 based on whether fingerprints file contains valid fingerprints data.

=item B<IsFingerprintsSDFile>

    $Status = $FingerprintsSDFileIO->IsFingerprintsSDFile($FileName);
    $Status = FileIO::FingerprintsSDFileIO::IsFingerprintsSDFile($FileName);

Returns 1 or 0 based on whether I<FileName> is a SD file.

=item B<Next or Read>

    $FingerprintsSDFileIO = $FingerprintsSDFileIO->Next();
    $FingerprintsSDFileIO = $FingerprintsSDFileIO->Read();

Reads next available compound fingerprints in SD file, processes the data, generates appropriate
fingerprints object, and returns B<FingerprintsSDFileIO>. The generated fingerprints object is available
using method B<GetFingerprints>.

=item B<SetBitStringFormat>

    $FingerprintsSDFileIO->SetBitStringFormat($Format);

Sets bit string I<Format> for fingerprints bit-vector string data in a SD file and returns B<FingerprintsSDFileIO>.
Possible values for B<BitStringFormat>: I<BinaryString or HexadecimalString>.

=item B<SetBitsOrder>

    $FingerprintsSDFileIO->SetBitsOrder($BitsOrder);

Sets I<BitsOrder> for fingerprints bit-vector string data in SD file and returns B<FingerprintsSDFileIO>.
Possible values for B<BitsOrder>: I<Ascending or Descending>.

=item B<SetCompoundIDMode>

    $FingerprintsSDFileIO->SetCompoundIDMode($Mode);

Sets compound ID I<Mode> for fingerprints bit-vector string data in a SD file and returns B<FingerprintsSDFileIO>.
Possible values for B<CompoundIDMode>: I<DataField, MolName, LabelPrefix, or MolNameOrLabelPrefix>.

=item B<SetCompoundString>

    $FingerprintsSDFileIO->SetCompoundString($CompoundString);

Sets I<CompoundString> and returns B<FingerprintsSDFileIO>.

=item B<SetDetailLevel>

    $FingerprintsSDFileIO->SetDetailLevel($Level);

Sets details I<Level> for generating diagnostics messages during SD file processing and returns
B<FingerprintsSDFileIO>. Possible values: I<Positive integers>.

=item B<SetFingerprints>

    $FingerprintsSDFileIO->SetFingerprints($FingerprintsObject);

Sets I<FingerprintsObject> for current data line and returns B<FingerprintsSDFileIO>.

=item B<SetFingerprintsString>

    $FingerprintsSDFileIO->SetFingerprintsString($FingerprintsString);

Sets I<FingerprintsString> for current data line and returns B<FingerprintsSDFileIO>.

=item B<SetFingerprintsStringMode>

    $FingerprintsSDFileIO->SetFingerprintsStringMode($Mode);

Sets I<FingerprintsStringMode> for SD file and returns B<FingerprintsFPFileIO>.
Possible values: I<AutoDetect, FingerprintsBitVectorString, or FingerprintsVectorString>

=item B<SetVectorStringFormat>

    $FingerprintsSDFileIO->SetVectorStringFormat($Format);

Sets I<VectorStringFormat> for SD file and returns B<FingerprintsFPFileIO>. Possible values:
I<IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString, ValuesAndIDsPairsString>.

=item B<WriteFingerprints>

    $FingerprintsFPFileIO->WriteFingerprints($FingerprintsObject,
                                            $CompoundID);

Writes fingerprints string generated from I<FingerprintsObject> object and other data including
I<CompoundID> to SD file and returns B<FingerprintsSDFileIO>.

=item B<WriteFingerprintsString>

    $FingerprintsSDFileIO->WriteFingerprints($FingerprintsString,
                                            $CompoundID);

Writes I<FingerprintsString> and other data including I<CompoundID> to SD file and returns
B<FingerprintsSDFileIO>.

Caveats:

    o FingerprintsStringMode, BitStringFormat, BitsOrder, VectorStringFormat
      values are ignored during writing of fingerprints and it's written to the file
      as it is.
    o CompoundString is not checked to remove any existing fingerprints data


=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

FingerprintsTextFileIO.pm, FingerprintsFPFileIO.pm, SDFileIO.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
