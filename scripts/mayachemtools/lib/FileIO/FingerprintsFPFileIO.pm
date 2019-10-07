package FileIO::FingerprintsFPFileIO;
#
# File: FingerprintsFPFileIO.pm
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
use TimeUtil ();
use Fingerprints::FingerprintsStringUtil ();
use PackageInfo ();
use FileIO::FileIO;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(FileIO::FileIO Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(IsFingerprintsFPFile);

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
  $This->_InitializeFingerprintsFPFileIO();

  $This->_InitializeFingerprintsFPFileIOProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeFingerprintsFPFileIO {
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

  # Fingeprints string for current line during read/write...
  $This->{FingerprintsString} = undef;

  # Partial fingeprints string corresponding to what's on the current line for current
  # line during read/write...
  $This->{PartialFingerprintsString} = undef;

  # Required header data keys and values during read/write...
  @{$This->{RequiredHeaderDataKeys}} = ();
  %{$This->{RequiredHeaderDataKeysAndValues}} = ();

  # First data line read/write...
  $This->{FirstDataLineIO} = 1;

  # Current fingerprints string data line number during read/write...
  $This->{LineNum} = 0;

  # FP line data during read/write...
  $This->{DataLine} = undef;

  # Initialize parameters for read...
  $This->_InitializeFingerprintsFPFileIORead();

  # Initialize parameters for write...
  $This->_InitializeFingerprintsFPFileIOWrite();

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

}

# Initialize object data for reading fingerprints FP file...
#
sub _InitializeFingerprintsFPFileIORead {
  my($This) = @_;

  # Header data keys and values...
  #
  @{$This->{HeaderDataKeys}} = ();
  %{$This->{HeaderDataKeysAndValues}} = ();
  %{$This->{CannonicalHeaderDataKeysAndValues}} = ();

  # By default, the fingerprints data is assumed to be valid and no validation is
  # performed before generating fingerprints objects...
  #
  $This->{ValidateData} = 1;

  # Level of detail to print during validation of data for invalid or missing data...
  $This->{DetailLevel} = 1;

  # Number of missing and invalid fingerprints string data lines...
  $This->{NumOfLinesWithMissingData} = 0;
  $This->{NumOfLinesWithInvalidData} = 0;

  # Compound ID for current fingerprints string...
  $This->{CompoundID} = undef;

  # Status of data in fingerprints FP file...
  $This->{ValidFileData} = 0;
  $This->{ValidRequiredHeaderDataKeys} = 0;
  $This->{ValidFingerprintsStringMode} = 0;

  return $This;
}

# Initialize object data for writing fingerprints FP file...
#
sub _InitializeFingerprintsFPFileIOWrite {
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
  # For fingerprints vector object containing vector NumericalValues, it corresponds to IDsAndValuesString; otherwise,
  # it's set to ValuesString.
  #
  $This->{VectorStringFormat} = undef;

  # Overwriting existing file...
  $This->{Overwrite} = 0;

  return $This;
}

# Initialize object values...
sub _InitializeFingerprintsFPFileIOProperties {
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
  if (!$This->IsFingerprintsFPFile($Name)) {
    croak "Error: ${ClassName}->New: Object can't be instantiated: File, $Name, doesn't appear to be fingerprints format...";
  }

  if ($This->GetMode() =~ /^Read$/i) {
    $This->_InitializeFingerprintsFPFileIOReadProperties(%NamesAndValues);
  }
  elsif ($This->GetMode() =~ /^(Write|Append)$/i) {
    $This->_InitializeFingerprintsFPFileIOWriteProperties(%NamesAndValues);
  }

  return $This;
}

# Initialize object properties for reading fingerprints FP file...
#
sub _InitializeFingerprintsFPFileIOReadProperties {
  my($This, %NamesAndValues) = @_;

  # Set default value for FingerprintsStringMode...
  if (!$This->{FingerprintsStringMode}) {
    $This->{FingerprintsStringMode} = 'AutoDetect';
  }

  $This->_PrepareForReadingFingerprintsFPFileData();

  return $This;
}

# Initialize object properties for writing fingerprints FP file...
#
sub _InitializeFingerprintsFPFileIOWriteProperties {
  my($This, %NamesAndValues) = @_;

  # Check FingerprintsStringMode value...
  if (!exists $NamesAndValues{FingerprintsStringMode}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying FingerprintsStringMode...";
  }

  if ($This->{FingerprintsStringMode} !~ /^(FingerprintsBitVectorString|FingerprintsVectorString)$/i) {
    croak "Error: ${ClassName}->New: Object can't be instantiated: FingerprintsStringMode value, $This->{FingerprintsStringMode}, is not valid; Supported values for write/append: FingerprintsBitVectorString or FingerprintsVectorString...";
  }

  $This->_PrepareForWritingFingerprintsFPFileData();

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

# Set compound ID...
#
sub SetCompoundID {
  my($This, $Value) = @_;

  if ($Value =~ / /) {
    $Value =~ s/ //g;
    carp "Warning: ${ClassName}->SetCompoundID: Spaces are not allowed in compound ID; They have been removed...";
  }

  $This->{CompoundID} = $Value;

  return $This;
}

# Set VectorStringFormat...
#
sub SetVectorStringFormat {
  my($This, $Value) = @_;

  # Possible values: IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString, ValuesAndIDsPairsString

  if ($Value !~ /^(IDsAndValuesString|IDsAndValuesPairsString|ValuesAndIDsString|ValuesAndIDsPairsString|ValuesString)$/i) {
    croak "Error: ${ClassName}->SetVectorStringFormat: FingerprintsStringMode value, $Value, is not valid; Supported values: IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString, ValuesAndIDsPairsString, or ValuesString...";
  }

  $This->{VectorStringFormat} = $Value;

  return $This;
}

# Get header data keys or number of header data keys in header data block...
#
sub GetHeaderDataKeys {
  my($This) = @_;

  return wantarray ? @{$This->{HeaderDataKeys}} : scalar @{$This->{HeaderDataKeys}};
}

# Set header data keys...
#
sub SetHeaderDataKeys {
  my($This, @Keys) = @_;

  croak "Error: ${ClassName}->SetHeaderDataKeys: Can't set HeaderDataKeys: Not allowed...";

  return $This;
}

# Get header data keys and values hash...
#
sub GetHeaderDataKeysAndValues {
  my($This) = @_;

  return %{$This->{HeaderDataKeysAndValues}};
}

# Set header data keys and values hash...
#
sub SetHeaderDataKeysAndValues {
  my($This, %KeysAndValues) = @_;

  croak "Error: ${ClassName}->SetHeaderDataKeysAndValues: Can't set HeaderDataKeysAndValues: Not allowed...";

  return $This;
}

# Get required header data keys or number of header data keys in header data block...
#
sub GetRequiredHeaderDataKeys {
  my($This) = @_;

  return wantarray ? @{$This->{RequiredHeaderDataKeys}} : scalar @{$This->{RequiredHeaderDataKeys}};
}

# Set required header data keys...
#
sub SetRequiredHeaderDataKeys {
  my($This, @Keys) = @_;

  croak "Error: ${ClassName}->SetRequiredHeaderDataKeys: Can't set RequiredHeaderDataKeys: Not allowed...";

  return $This;
}

# Get required header data keys and values hash...
#
sub GetRequiredHeaderDataKeysAndValues {
  my($This) = @_;

  return %{$This->{RequiredHeaderDataKeysAndValues}};
}

# Set required header data keys and values hash...
#
sub SetRequiredHeaderDataKeysAndValues {
  my($This, %KeysAndValues) = @_;

  croak "Error: ${ClassName}->SetRequiredHeaderDataKeysAndValues: Can't set RequiredHeaderDataKeysAndValues: Not allowed...";

  return $This;
}

# Get fingerprints object for current data line...
#
sub GetFingerprints {
  my($This) = @_;

  return $This->{FingerprintsObject};
}

# Set fingerprints object for current data line...
#
sub SetFingerprints {
  my($This, $FingerprintsObject) = @_;

  $This->{FingerprintsObject} = $FingerprintsObject;

  return $This;
}

# Get fingerprints string  for current data line...
#
sub GetFingerprintsString {
  my($This) = @_;

  return $This->{FingerprintsString} ? $This->{FingerprintsString} : 'None';
}

# Set fingerprints string for current data line...
#
sub SetFingerprintsString {
  my($This, $FingerprintsString) = @_;

  $This->{FingerprintsString} = $FingerprintsString;

  return $This;
}

# Get partial fingerprints string  for current data line...
#
sub GetPartialFingerprintsString {
  my($This) = @_;

  return $This->{PartialFingerprintsString} ? $This->{PartialFingerprintsString} : 'None';
}

# Set partial fingerprints string for current data line...
#
sub SetPartialFingerprintsString {
  my($This, $PartialFingerprintsString) = @_;

  $This->{PartialFingerprintsString} = $PartialFingerprintsString;

  return $This;
}

# Does fingerprints FP file contain valid data?
#
sub IsFingerprintsFileDataValid {
  my($This) = @_;

  return $This->{ValidFileData} ? 1 : 0;
}

# Does current data line contains valid fingerprints object data?
#
sub IsFingerprintsDataValid {
  my($This) = @_;

  return defined $This->{FingerprintsObject} ? 1 : 0;
}

# Check presence of a header data key...
#
sub IsHeaderDataKeyPresent {
  my($This, $Key) = @_;
  my($CannonicalKey);

  $CannonicalKey = lc $Key;

  return exists $This->{CannonicalHeaderDataKeysAndValues}{$CannonicalKey} ? 1 : 0;
}

# Get value of header data key...
#
sub GetHeaderDataKeyValue {
  my($This, $Key) = @_;
  my($CannonicalKey);

  $CannonicalKey = lc $Key;

  return exists $This->{CannonicalHeaderDataKeysAndValues}{$CannonicalKey} ?  $This->{CannonicalHeaderDataKeysAndValues}{$CannonicalKey} : undef;
}

#
# Read next available fingerprints line, process it and generate appropriate fingerprints
# objects...
#
sub Read {
  my($This) = @_;

  # Read data line...
  if (!$This->_ReadDataLine()) {
    return undef;
  }

  # No need to process invalid FP file with invalid data...
  if (!$This->{ValidFileData}) {
    if ($This->{ValidateData}) {
      $This->{NumOfLinesWithMissingData} += 1;
    }
    return $This;
  }

  # Perform data validation...
  if ($This->{ValidateData}) {
    if (!$This->_ValidateReadDataLine()) {
      return $This;
    }
  }

  # Check again to handle problematic data for non-validated data lines...
  if (!$This->{FingerprintsString}) {
    return $This;
  }

  # Generate fingeprints object...
  $This->_GenerateFingerprintsObject();

  # Setup fingerprints compound ID for fingerprints string...
  $This->_GenerateCompoundID();

  return $This;
}

# Read next available fingerprints line, process it and generate appropriate fingerprints
# objects...
#
sub Next {
  my($This) = @_;

  return $This->Read();
}

# Read fingerprints data line line...
#
sub _ReadDataLine {
  my($This) = @_;

  # Initialize data for current line...
  $This->_InitializeReadDataLine();

  if ($This->{FirstDataLineIO}) {
    # Get first data line...
    $This->_ProcessFirstDataLineRead();
  }
  else {
    # Get next data line...
    $This->{LineNum} += 1;
    $This->{DataLine} = TextUtil::GetTextLine($This->{FileHandle});
  }

  # Is it end of file?
  if (!$This->{DataLine}) {
    return 0;
  }

  # Process data line to retrieve compound ID and fingerprints string information...
  $This->_ProcessDataLineRead();

  return 1;
}

# Process data line to retrieve compound ID and fingerprints string information...
#
sub _ProcessDataLineRead {
  my($This) = @_;
  my($CompoundID, $PartialFingerprintsString);

  ($CompoundID, $PartialFingerprintsString) = $This->{DataLine} =~ /^(.*?)[ ]+(.*?)$/;

  if (!(defined($CompoundID) && defined($PartialFingerprintsString))) {
    return $This;
  }

  $This->{CompoundID} = $CompoundID;
  $This->{PartialFingerprintsString} = $PartialFingerprintsString;

  # Set up fingerprints string...
  $This->_GenerateFingerprintsStringFromPartialFingerprintsString();

  return $This;
}

# Initialize data line for reading...
#
sub _InitializeReadDataLine {
  my($This) = @_;

  $This->{CompoundID} = undef;
  $This->{DataLine} = undef;

  $This->{FingerprintsObject} = undef;

  $This->{FingerprintsString} = undef;
  $This->{PartialFingerprintsString} = undef;

  return $This;
}

# Validate fingerprints string data line...
#
sub _ValidateReadDataLine {
  my($This) = @_;

  # Check for missing data...
  if (!($This->{CompoundID} && $This->{PartialFingerprintsString})) {
    # Missing data...
    $This->{NumOfLinesWithMissingData} += 1;
    if ($This->{DetailLevel} >= 3) {
      carp "Warning: ${ClassName}->_ValidateReadDataLine: Data line number $This->{LineNum} contains no fingerprints data: $This->{DataLine}...";
    }
    elsif ($This->{DetailLevel} >= 2) {
      carp "Warning: ${ClassName}->_ValidateReadDataLine: Data line number $This->{LineNum} contains no fingerprints data...";
    }
    return 0;
  }

  # Check for invalid data...
  my($InvalidFingerprintsData);

  $InvalidFingerprintsData = 0;
  if ($This->{FingerprintsString}) {
    $InvalidFingerprintsData = Fingerprints::FingerprintsStringUtil::AreFingerprintsStringValuesValid($This->{FingerprintsString}) ? 0 : 1;
  }
  else {
    $InvalidFingerprintsData = 1;
  }

  if ($InvalidFingerprintsData) {
    $This->{NumOfLinesWithInvalidData} += 1;
    if ($This->{DetailLevel} >= 3) {
      carp "Warning: ${ClassName}->_ValidateReadDataLine: Data line number $This->{LineNum} contains invalid fingerprints data: $This->{DataLine}...";
    }
    elsif ($This->{DetailLevel} >= 2) {
      carp "Warning: ${ClassName}->_ValidateReadDataLine: Data line number $This->{LineNum} contains invalid fingerprints data...";
    }
    return 0;
  }

  return 1;
}

# Setup fingerprints compound ID for fingerprints string...
sub _GenerateCompoundID {
  my($This) = @_;

  # Set fingerprints ID...
  if ($This->{FingerprintsObject}) {
    $This->{FingerprintsObject}->SetID($This->{CompoundID});
  }

  return $This;
}

# Process first read...
#
sub _ProcessFirstDataLineRead {
  my($This) = @_;
  my($Line);

  $This->{FirstDataLineIO} = 0;

  # Skip over header data lines and collect first data line...

  LINE: while ($Line = TextUtil::GetTextLine($This->{FileHandle})) {
    $This->{LineNum} += 1;

    # Is it a header data line?
    if ($Line =~ /^#/) {
      next LINE;
    }
    $This->{DataLine} = $Line;
    last LINE;
  }

  return $This;
}

# Get ready for reading fingerprints FP file...
#
sub _PrepareForReadingFingerprintsFPFileData {
  my($This) = @_;

  # Retrieve FP file data headers information....
  $This->_RetrieveFPFileDataHeaders();

  # Validate header data keys and values information...
  $This->_ValidateReadHeaderDataKeysAndValues();

  # Validate fingeprints string mode information...
  if ($This->{ValidRequiredHeaderDataKeys}) {
    $This->_ValidateReadFingerprintsStringMode();
  }

  # Set status of FP file data...
  $This->{ValidFileData} = ($This->{ValidRequiredHeaderDataKeys} && $This->{ValidFingerprintsStringMode}) ? 1 : 0;

  return $This;
}

# Retrieve information about fingerprints date header in FP file...
#
sub _RetrieveFPFileDataHeaders {
  my($This) = @_;
  my($FPFile, $Line, $Index, $KeyValuePair, $Key, $Value, $KeyValueDelimiter, $KeyValuePairDelimiter, @LineKeyValuePairs);

  $FPFile = $This->{Name};

  if (!(-e $FPFile)) {
    croak "Error: ${ClassName}->_RetrieveFPFileDataHeaders: File, $FPFile, doesn't exist...";
  }

  if (!open FPFILE, "$FPFile") {
    croak "Error: ${ClassName}->_RetrieveFPFileDataHeaders: Couldn't open input FP file $FPFile: $! ...";
  }

  # Process header key/value pair data...
  #
  $KeyValueDelimiter = '=';
  $KeyValuePairDelimiter = ';';

  @{$This->{HeaderDataKeys}} = ();
  %{$This->{HeaderDataKeysAndValues}} = ();
  %{$This->{CannonicalHeaderDataKeysAndValues}} = ();

  LINE: while ($Line = TextUtil::GetTextLine(\*FPFILE)) {
    # Is it a key/value pairs line?
    if ($Line !~ /^#/) {
      last LINE;
    }

    # Take out starting hash mark before processing key/value pairs...
    $Line =~ s/^#//;
    if (TextUtil::IsEmpty($Line)) {
      next LINE;
    }

    @LineKeyValuePairs = ();

    for $KeyValuePair (split "$KeyValuePairDelimiter", $Line) {
      ($Key, $Value) = split "$KeyValueDelimiter", $KeyValuePair;

      $Key = defined($Key) ? TextUtil::RemoveLeadingAndTrailingWhiteSpaces($Key) : '';
      $Value = defined($Value) ? TextUtil::RemoveLeadingAndTrailingWhiteSpaces($Value) : '';

      if (TextUtil::IsEmpty($Key) || TextUtil::IsEmpty($Value)) {
	carp "Warning: ${ClassName}->_RetrieveFPFileDataHeaders: Data header line containing \"Key = Value\" pairs is not valid: It must contain even number of \"Key = Value\" pairs with valid values. Ignoring data header line: \"$Line\"...";
	next LINE;
      }
      push @{$This->{HeaderDataKeys}}, $Key;
      push @LineKeyValuePairs, ($Key, $Value);
    }

    for ($Index = 0; $Index < $#LineKeyValuePairs; $Index += 2) {
      $Key = $LineKeyValuePairs[$Index]; $Value = $LineKeyValuePairs[$Index + 1];

      $This->{HeaderDataKeysAndValues}{$Key} = $Value;
      $This->{CannonicalHeaderDataKeysAndValues}{lc($Key)} = $Value;
    }
  }
  close FPFILE;

  return $This;
}

# Validate header data and keys...
#
sub _ValidateReadHeaderDataKeysAndValues {
  my($This) = @_;
  my($FingerprintsStringType, $Key, $Value, @RequiredHeaderDataKeys);

  $This->{ValidRequiredHeaderDataKeys} = 0;
  @{$This->{RequiredHeaderDataKeys}} = ();

  # Is FingerprintsStringType key is present?
  if (!$This->IsHeaderDataKeyPresent('FingerprintsStringType')) {
    carp "carp: ${ClassName}->_ValidateReadHeaderDataKeysAndValues: FingerprintsStringType data header key is missing in fingerprints file...";
    return 0;
  }
  $FingerprintsStringType = $This->GetHeaderDataKeyValue('FingerprintsStringType');

  # Are all required data header keys present?
  #
  @RequiredHeaderDataKeys = ();

  if ($FingerprintsStringType =~ /^(FingerprintsBitVector|FingerprintsVector)$/i) {
    push @RequiredHeaderDataKeys, $This->_GetRequiredHeaderDataKeys($FingerprintsStringType);
  }
  else {
    carp "Warning: ${ClassName}->_ValidateReadHeaderDataKeysAndValues: FingerprintsStringType data header key value, $FingerprintsStringType, is not valid. SUpported values: FingerprintsBitVector or FingerprintsVector...";
    return 0;
  }

  for $Key (@RequiredHeaderDataKeys) {
    if (!$This->IsHeaderDataKeyPresent($Key)) {
      croak "Error: ${ClassName}->_ValidateReadHeaderDataKeysAndValues: Requires data header key, $Key, is missing in fingerprints file...";
    }
  }

  push @{$This->{RequiredHeaderDataKeys}}, @RequiredHeaderDataKeys;

  # Are all required data header key values valid?
  #
  if (!$This->_ValidateRequiredHeaderDataKeyValues()) {
    return 0;
  }

  # Process required header key values...
  #
  $This->_ProcessRequiredHeaderDataKeyValues();

  $This->{ValidRequiredHeaderDataKeys} = 1;

  return 1;
}

# Validate data header key values....
#
sub _ValidateRequiredHeaderDataKeyValues {
  my($This) = @_;
  my($Key, $Value);

  for $Key (@{$This->{RequiredHeaderDataKeys}}) {
    $Value = $This->GetHeaderDataKeyValue($Key);
    KEY: {
      if ($Key =~ /^FingerprintsStringType$/i) {
	if ($Value !~ /^(FingerprintsBitVector|FingerprintsVector)$/i) {
	  carp "Warning: ${ClassName}->_ValidateRequiredHeaderDataKeyValues: Required $Key data header key value, $Value, is not valid. Supported values: FingerprintsBitVector or FingerprintsVector...";
	  return 0;
	}
	last KEY;
      }
      if ($Key =~ /^Size$/i) {
	if (!TextUtil::IsPositiveInteger($Value)) {
	  carp "Warning: ${ClassName}->_ValidateRequiredHeaderDataKeyValues: Required $Key data header key value, $Value, is not valid. Supported values: > 0...";
	  return 0;
	}
	last KEY;
      }
      if ($Key =~ /^BitStringFormat$/i) {
	if ($Value !~ /^(BinaryString|HexadecimalString)$/i) {
	  carp "Warning: ${ClassName}->_ValidateRequiredHeaderDataKeyValues: Required $Key data header key value, $Value, is not valid. Supported values: BinaryString or HexadecimalString ...";
	  return 0;
	}
	last KEY;
      }
      if ($Key =~ /^BitsOrder$/i) {
	if ($Value !~ /^(Ascending|Descending)$/i) {
	  carp "Warning: ${ClassName}->_ValidateRequiredHeaderDataKeyValues: Required $Key data header key value, $Value, is not valid. Supported values: Ascending or Descending...";
	  return 0;
	}
	last KEY;
      }
      if ($Key =~ /^VectorStringFormat$/i) {
	if ($Value !~ /^(IDsAndValuesString|IDsAndValuesPairsString|ValuesAndIDsString|ValuesAndIDsPairsString|ValuesString)$/i) {
	  carp "Warning: ${ClassName}->_ValidateRequiredHeaderDataKeyValues: Required $Key data header key value, $Value, is not valid. Supported values: IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString, ValuesAndIDsPairsString, or ValuesString ...";
	  return 0;
	}
	last KEY;
      }
      if ($Key =~ /^VectorValuesType$/i) {
	if ($Value !~ /^(OrderedNumericalValues|NumericalValues|AlphaNumericalValues)$/i) {
	  carp "Warning: ${ClassName}->_ValidateRequiredHeaderDataKeyValues: Required $Key data header key value, $Value, is not valid. Supported values: OrderedNumericalValues, NumericalValues or AlphaNumericalValues...";
	  return 0;
	}
	last KEY;
      }
      if ($Key =~ /^Description$/i) {
	if (TextUtil::IsEmpty($Value)) {
	  carp "Warning: ${ClassName}->_ValidateRequiredHeaderDataKeyValues: Required $Key data header key value is not valid. Supported value: A no-empty text string...";
	  return 0;
	}
	last KEY;
      }
      carp "Warning: ${ClassName}->_ValidateRequiredHeaderDataKeyValues: Required $Key data header key is not not supported...";
      return 0;
    }
  }

  return 1;
}

# Process required header key valeues for access during complete fingerprints
# string generation from a partial fingerprints string specified on fingerprints
# line...
#
sub _ProcessRequiredHeaderDataKeyValues {
  my($This) = @_;
  my($Key, $Value, @Keys);

  %{$This->{RequiredHeaderDataKeysAndValues}} = ();

  for $Key (@{$This->{RequiredHeaderDataKeys}}) {
    $Value = $This->GetHeaderDataKeyValue($Key);
    $This->{RequiredHeaderDataKeysAndValues}{$Key} = $Value;
  }

  # Setup prefixes for generating fingerprints strings...
  $This->{FingerprintsBitVectorStringPrefix} = '';
  $This->{FingerprintsVectorStringPrefix1} = '';
  $This->{FingerprintsVectorStringPrefix2} = '';

  if ($This->{RequiredHeaderDataKeysAndValues}{FingerprintsStringType} =~ /^FingerprintsBitVector$/i) {
    @Keys = qw(FingerprintsStringType Description Size BitStringFormat BitsOrder);
    $This->{FingerprintsBitVectorStringPrefix} = $This->_GenerateFingerprintsPrefixUsingKeys(@Keys);
  }
  elsif ($This->{RequiredHeaderDataKeysAndValues}{FingerprintsStringType} =~ /^FingerprintsVector$/i) {
    @Keys = qw(FingerprintsStringType Description);
    $This->{FingerprintsVectorStringPrefix1} = $This->_GenerateFingerprintsPrefixUsingKeys(@Keys);

    @Keys = qw(VectorValuesType VectorStringFormat);
    $This->{FingerprintsVectorStringPrefix2} = $This->_GenerateFingerprintsPrefixUsingKeys(@Keys);
  }

  return $This;
}

# Generate fingerprints prefix using header keys data...
#
sub _GenerateFingerprintsPrefixUsingKeys {
  my($This, @Keys) = @_;
  my($Delimiter, $Key, @Values);

  $Delimiter = Fingerprints::FingerprintsStringUtil::GetFingeprintsStringDelimiter();

  @Values = ();
  for $Key (@Keys) {
    push @Values, $This->{RequiredHeaderDataKeysAndValues}{$Key};
  }

  return join($Delimiter, @Values)
}

# Get required header data keys...
#
sub _GetRequiredHeaderDataKeys {
  my($This, $FingerprintsStringType) = @_;
  my(@RequiredKeys);

  @RequiredKeys = ();

  if ($FingerprintsStringType =~ /FingerprintsBitVector$/i) {
    push @RequiredKeys, qw(FingerprintsStringType Description Size BitStringFormat BitsOrder);
  }
  elsif ($FingerprintsStringType =~ /^FingerprintsVector/i) {
    push @RequiredKeys, qw(FingerprintsStringType Description VectorStringFormat VectorValuesType);
  }
  else {
    carp "Warning: ${ClassName}->GetRequiredHeaderDataKeys: FingerprintsStringType value, $FingerprintsStringType, is not valid. Supported values: FingerprintsBitVector or FingerprintsVector...";
  }

  return @RequiredKeys;
}

# Validate fingerprints string mode information...
#
sub _ValidateReadFingerprintsStringMode {
  my($This) = @_;
  my($FingerprintsStringType, $FingerprintsStringDescription, $FingerprintsBitVectorStringMode, $FingerprintsVectorStringMode, $FirstFingerprintsStringType, $FirstFingerprintsStringDescription);

  $This->{ValidFingerprintsStringMode} = 0;
  $This->{FingerprintsBitVectorStringMode} = 0;
  $This->{FingerprintsVectorStringMode} = 0;

  $This->{FirstFingerprintsStringType} = '';
  $This->{FirstFingerprintsStringDescription} = '';

  $FingerprintsBitVectorStringMode = 0;
  $FingerprintsVectorStringMode = 0;

  $FirstFingerprintsStringType = '';
  $FirstFingerprintsStringDescription = '';

  $FingerprintsStringType = $This->GetHeaderDataKeyValue('FingerprintsStringType');
  $FingerprintsStringDescription = $This->GetHeaderDataKeyValue('Description');

  if ($This->{FingerprintsStringMode} =~ /^FingerprintsBitVectorString$/i) {
    if ($FingerprintsStringType !~ /^FingerprintsBitVector$/i) {
      carp "Warning: ${ClassName}->_ValidateReadFingerprintsStringMode: Fingerprints string data type, $FingerprintsStringType, doesn't correspond to, FingerprintsBitVectorString, specified using \"FingerprintsStringMode\"...";
      return 0;
    }
    $FingerprintsBitVectorStringMode = 1;
    $FirstFingerprintsStringType = 'FingerprintsBitVector';
    $FirstFingerprintsStringDescription = $FingerprintsStringDescription;
  }
  elsif ($This->{FingerprintsStringMode} =~ /^FingerprintsVectorString$/i) {
    if ($FingerprintsStringType !~ /^FingerprintsVector$/i) {
      carp "Warning: ${ClassName}->_ValidateReadFingerprintsStringMode: Fingerprints string data type, $FingerprintsStringType, doesn't correspond to, FingerprintsVectorString, specified using \"FingerprintsStringMode\"...";
      return 0;
    }
    $FingerprintsVectorStringMode = 1;
    $FirstFingerprintsStringType = 'FingerprintsVector';
    $FirstFingerprintsStringDescription = $FingerprintsStringDescription;
  }
  else {
    # AutoDetect mode...
    if ($FingerprintsStringType =~ /^FingerprintsBitVector$/i) {
      $FingerprintsBitVectorStringMode = 1;
    }
    elsif ($FingerprintsStringType =~ /^FingerprintsVector$/i) {
      $FingerprintsVectorStringMode = 1;
    }
    else {
      carp "Warning: ${ClassName}->_ValidateReadFingerprintsStringMode: Fingerprints string data type, $FingerprintsStringType, identified during, AutoDetect, value of \"FingerprintsStringMode\" is not valid; Supported fingerprints types: FingerprintBitVector or FingerprintsVector...";
      return 0;
    }
    $FirstFingerprintsStringType = $FingerprintsStringType;
    $FirstFingerprintsStringDescription = $FingerprintsStringDescription;
  }

  $This->{ValidFingerprintsStringMode} = 1;

  $This->{FingerprintsBitVectorStringMode} = $FingerprintsBitVectorStringMode;
  $This->{FingerprintsVectorStringMode} = $FingerprintsVectorStringMode;

  $This->{FirstFingerprintsStringType} = $FirstFingerprintsStringType;
  $This->{FirstFingerprintsStringDescription} = $FirstFingerprintsStringDescription;

  return 1;
}

# Write fingerprints string generated from specified fingerprints - fingerprints-bit vector or
# fingerprints vector - object and other data to FP file...
#
sub WriteFingerprints {
  my($This, $FingerprintsObject, $CompoundID) = @_;

  # Initialize data for current line...
  $This->_InitializeWriteDataLine();

  # Set fingerprints object and compound ID...
  $This->{FingerprintsObject} = $FingerprintsObject;
  $This->SetCompoundID($CompoundID);

  # Generate fingerprints string...
  $This->_GenerateFingerprintsString();

  # Generate partial fingerprints string...
  $This->_GeneratePartialFingerprintsStringFromFingerprintsString();

  # Write data line..
  $This->_WriteDataLine();

  return $This;
}

# Write fingerprints string and other data to FP file...
#
# Notes:
#   o FingerprintsStringMode, BitStringFormat, BitsOrder, VectorStringFormat values
#     are ignored during writing of fingerprints and it's written to the file as it is.
#   o FingerprintsString is a regular fingerprints string as oppose to a partial fingerprints
#     string.
#
sub WriteFingerprintsString {
  my($This, $FingerprintsString, $CompoundID) = @_;

  # Initialize data for current line...
  $This->_InitializeWriteDataLine();

  # Set fingerprints string and compound ID...
  $This->{FingerprintsString} = $FingerprintsString;
  $This->SetCompoundID($CompoundID);

  # Generate fingerprints object...
  $This->_GenerateFingerprintsObject();

  # Generate partial fingerprints string...
  $This->_GeneratePartialFingerprintsStringFromFingerprintsString();

  # Write data line..
  $This->_WriteDataLine();

  return $This;
}

# Initialize data line for reading...
#
sub _InitializeWriteDataLine {
  my($This) = @_;

  $This->{DataLine} = undef;
  $This->{CompoundID} = undef;

  $This->{FingerprintsObject} = undef;

  $This->{FingerprintsString} = undef;
  $This->{PartialFingerprintsString} = undef;

  return $This;
}

# Write fingerprints data line line...
#
sub _WriteDataLine {
  my($This) = @_;
  my($FileHandle, $Line);

  if ($This->{FirstDataLineIO}) {
    $This->_ProcessFirstDataLineWrite();
  }

  # Write data compound ID along with partial fingerprints string...
  $Line = $This->{CompoundID} . ' ' . $This->{PartialFingerprintsString};

  $This->{LineNum} += 1;
  $FileHandle = $This->{FileHandle};
  print $FileHandle "$Line\n";

  $This->{DataLine} = $Line;

  return $This;
}

# Process first write...
#
sub _ProcessFirstDataLineWrite {
  my($This) = @_;
  my($Line, $FileHandle);

  $This->{FirstDataLineIO} = 0;

  if ($This->GetMode() =~ /^Write$/i) {
    # Skip it for append mode...
    $This->_WritePackageAndTimeStampHeaderKeys();
    $This->_WriteRequiredHeaderDataKeys();
  }

  return $This;
}

# Write out package and time stamp information...
#
sub _WritePackageAndTimeStampHeaderKeys {
  my($This) = @_;
  my($FileHandle, $Key, $Value);

  $FileHandle = $This->{FileHandle};

  # Package information...
  $This->{LineNum} += 1;
  $Key = "Package"; $Value = PackageInfo::GetPackageName() . " " . PackageInfo::GetVersionNumber();
  print $FileHandle "# $Key = $Value\n";

  $This->{LineNum} += 1;
  $Key = "Release Date"; $Value = PackageInfo::GetReleaseDate();
  print $FileHandle "# $Key = $Value\n";

  # Timestamp information...
  $This->{LineNum} += 1;
  print $FileHandle "#\n";

  $This->{LineNum} += 1;
  $Key = "TimeStamp"; $Value = TimeUtil::FPFileTimeStamp();
  print $FileHandle "# $Key = $Value\n";

  return $This;
}

# Write out required header data keys...
#
sub _WriteRequiredHeaderDataKeys {
  my($This) = @_;
  my($FileHandle, $Key, $Value);

  $FileHandle = $This->{FileHandle};

  $This->_GenerateWriteRequiredHeaderDataKeys();

  $This->{LineNum} += 1;
  print $FileHandle "#\n";

  for $Key (@{$This->{RequiredHeaderDataKeys}}) {
    $Value = $This->{RequiredHeaderDataKeysAndValues}{$Key};

    $This->{LineNum} += 1;
    print $FileHandle "# $Key = $Value\n";

    if ($Key =~ /^FingerprintsStringType$/i) {
      $This->{LineNum} += 1;
      print $FileHandle "#\n";
    }
  }

  $This->{LineNum} += 1;
  print $FileHandle "#\n";

  return $This;
}

sub _GenerateWriteRequiredHeaderDataKeys {
  my($This) = @_;

  if ($This->{FingerprintsBitVectorStringMode} && ($This->{FingerprintsString} =~ /^FingerprintsBitVector/i)) {
    $This->_GenerateWriteRequiredHeaderDataKeysForBitVectorString();
  }
  elsif ($This->{FingerprintsVectorStringMode} && ($This->{FingerprintsString} =~ /^FingerprintsVector/i)) {
    $This->_GenerateWriteRequiredHeaderDataKeysForVectorString();
  }
  else {
    croak "Error: ${ClassName}->_GenerateWriteRequiredHeaderDataKeys: Required header data keys can't be generated: FingerprintsStringMode value, $This->{FingerprintsStringMode}, doesn't correspond to type of first FingerprintsString: $This->{FingerprintsString}...";
  }

  return $This;
}

# Generate required data header keys and values for writing fingerprints bit vector string...
#
sub _GenerateWriteRequiredHeaderDataKeysForBitVectorString {
  my($This) = @_;
  my($Key, $VectorType, $Description, $Size, $BitStringFormat, $BitsOrder);

  @{$This->{RequiredHeaderDataKeys}} = ();
  push @{$This->{RequiredHeaderDataKeys}}, $This->_GetRequiredHeaderDataKeys('FingerprintsBitVector');

  ($VectorType, $Description, $Size, $BitStringFormat, $BitsOrder) = Fingerprints::FingerprintsStringUtil::GetFingerprintsStringValues($This->{FingerprintsString});

  %{$This->{RequiredHeaderDataKeysAndValues}} = ();

  for $Key (@{$This->{RequiredHeaderDataKeys}}) {
    KEYTYPE: {
      if ($Key =~ /^FingerprintsStringType$/i) {
	$This->{RequiredHeaderDataKeysAndValues}{$Key} = $VectorType;
	last KEYTYPE;
      }
      if ($Key =~ /^Description$/i) {
	$This->{RequiredHeaderDataKeysAndValues}{$Key} = $Description;
	last KEYTYPE;
      }
      if ($Key =~ /^Size$/i) {
	$This->{RequiredHeaderDataKeysAndValues}{$Key} = $Size;
	last KEYTYPE;
      }
      if ($Key =~ /^BitStringFormat$/i) {
	$This->{RequiredHeaderDataKeysAndValues}{$Key} = $BitStringFormat;
	last KEYTYPE;
      }
      if ($Key =~ /^BitsOrder$/i) {
	$This->{RequiredHeaderDataKeysAndValues}{$Key} = $BitsOrder;
	last KEYTYPE;
      }
      croak "Error: ${ClassName}->_GenerateWriteRequiredHeaderDataKeysForBitVectorString: Required header data key, $Key, value can't be generated: It's not a known key ...";
    }
  }

  return $This;
}

# Generate required data header keys and values for writing fingerprints vector string...
#
sub _GenerateWriteRequiredHeaderDataKeysForVectorString {
  my($This) = @_;
  my($Key, $Value,  $VectorType, $Description, $NumOfValues, $VectorValuesType, $VectorStringFormat);

  @{$This->{RequiredHeaderDataKeys}} = ();
  push @{$This->{RequiredHeaderDataKeys}}, $This->_GetRequiredHeaderDataKeys('FingerprintsVector');

  ($VectorType, $Description, $NumOfValues, $VectorValuesType, $VectorStringFormat) = Fingerprints::FingerprintsStringUtil::GetFingerprintsStringValues($This->{FingerprintsString});

  %{$This->{RequiredHeaderDataKeysAndValues}} = ();

  for $Key (@{$This->{RequiredHeaderDataKeys}}) {
    KEYTYPE: {
      if ($Key =~ /^FingerprintsStringType$/i) {
	$This->{RequiredHeaderDataKeysAndValues}{$Key} = $VectorType;
	last KEYTYPE;
      }
      if ($Key =~ /^Description$/i) {
	$This->{RequiredHeaderDataKeysAndValues}{$Key} = $Description;
	last KEYTYPE;
      }
      if ($Key =~ /^VectorValuesType$/i) {
	$This->{RequiredHeaderDataKeysAndValues}{$Key} = $VectorValuesType;
	last KEYTYPE;
      }
      if ($Key =~ /^VectorStringFormat$/i) {
	$This->{RequiredHeaderDataKeysAndValues}{$Key} = $VectorStringFormat;
	last KEYTYPE;
      }
      croak "Error: ${ClassName}->_GenerateWriteRequiredHeaderDataKeysForVectorString: Required header data key, $Key, value can't be generated: It's not a known key ...";
    }
  }

  return $This;
}


# Get ready for writing fingerprints FP file...
#
sub _PrepareForWritingFingerprintsFPFileData {
  my($This) = @_;
  my($FPFile, $FileDir, $FileName, $FileExt, $OutDelim);

  $FPFile = $This->{Name};
  if (!$This->{Overwrite}) {
    if (-e $FPFile) {
      croak "Error: ${ClassName}->_PrepareForWritingFingerprintsFPFileData: File, $FPFile, already exist. Use overwrite option...";
    }
  }

  # Setup FingerprintsStringMode status...
  #
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

# Generate fingerprints string using partial fingerprints string and header keys data...
#
# Notes:
#   o FP file fingerprints data line only contain partial fingerprints data which
#     can't be used directly to create fingerprints bit-vector or vector objects
#     using functions available in FingerprintsStringUtil.pm module
#
sub _GenerateFingerprintsStringFromPartialFingerprintsString {
  my($This) = @_;
  my($FPStringDelim);

  $This->{FingerprintsString} = '';

  if (!$This->{PartialFingerprintsString}) {
    return $This;
  }

  $FPStringDelim = Fingerprints::FingerprintsStringUtil::GetFingeprintsStringDelimiter();

  if ($This->{FingerprintsBitVectorStringMode}) {
    $This->{FingerprintsString} = $This->{FingerprintsBitVectorStringPrefix} . $FPStringDelim . $This->{PartialFingerprintsString};
  }
  elsif ($This->{FingerprintsVectorStringMode}) {
    my($NumOfValues, $VectorStringData);

    ($NumOfValues, $VectorStringData) =  $This->{PartialFingerprintsString} =~ /^(.*?)$FPStringDelim(.*?)$/;
    if (!(defined($NumOfValues) && defined($VectorStringData) && $VectorStringData)) {
      return $This;
    }

    $This->{FingerprintsString} = $This->{FingerprintsVectorStringPrefix1} . $FPStringDelim . $NumOfValues . $FPStringDelim . $This->{FingerprintsVectorStringPrefix2} . $FPStringDelim . $VectorStringData;
  }

  return $This;
}

# Generate partial fingerprints string using fingerprints string and header keys data...
#
# Notes:
#   o FP file fingerprints data line only contain partial fingerprints data which
#     can't be used directly to create fingerprints bit-vector or vector objects
#     using functions available in FingerprintsStringUtil.pm module
#
sub _GeneratePartialFingerprintsStringFromFingerprintsString {
  my($This) = @_;

  $This->{PartialFingerprintsString} = '';

  if (!$This->{FingerprintsString}) {
    return $This;
  }

  if ($This->{FingerprintsBitVectorStringMode}) {
    my($VectorType, $Description, $Size, $BitStringFormat, $BitsOrder, $BitVectorString);

    ($VectorType, $Description, $Size, $BitStringFormat, $BitsOrder, $BitVectorString) = Fingerprints::FingerprintsStringUtil::GetFingerprintsStringValues($This->{FingerprintsString});
    $This->{PartialFingerprintsString} = $BitVectorString;
  }
  elsif ($This->{FingerprintsVectorStringMode}) {
    my($FPStringDelim, $VectorType, $Description, $NumOfValues, $VectorValuesType, $VectorStringFormat, $VectorString1, $VectorString2, $VectorString);

    $FPStringDelim = Fingerprints::FingerprintsStringUtil::GetFingeprintsStringDelimiter();

    ($VectorType, $Description, $NumOfValues, $VectorValuesType, $VectorStringFormat, $VectorString1, $VectorString2) = Fingerprints::FingerprintsStringUtil::GetFingerprintsStringValues($This->{FingerprintsString});
    $VectorString = TextUtil::IsEmpty($VectorString2) ? $VectorString1 : "${VectorString1}${FPStringDelim}${VectorString2}";

    $This->{PartialFingerprintsString} = $NumOfValues . $FPStringDelim . $VectorString;
  }

  return $This;
}

# Is it a fingerprints file?
sub IsFingerprintsFPFile ($;$) {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $FileName, $Status);

  if ((@_ == 2) && (_IsFingerprintsFPFileIO($FirstParameter))) {
    ($This, $FileName) = ($FirstParameter, $SecondParameter);
  }
  else {
    $FileName = $FirstParameter;
  }

  # Check file extension...
  $Status = FileUtil::CheckFileType($FileName, "fpf fp");

  return $Status;
}

# Is it a FingerprintsFPFileIO object?
sub _IsFingerprintsFPFileIO {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

FingerprintsFPFileIO

=head1 SYNOPSIS

use FileIO::FingerprintsFPFileIO;

use FileIO::FingerprintsFPFileIO qw(:all);

=head1 DESCRIPTION

B<FingerprintsFPFileIO> class provides the following methods:

new, GetFingerprints, GetFingerprintsString, GetHeaderDataKeyValue,
GetHeaderDataKeys, GetHeaderDataKeysAndValues, GetPartialFingerprintsString,
GetRequiredHeaderDataKeys, GetRequiredHeaderDataKeysAndValues,
IsFingerprintsDataValid, IsFingerprintsFPFile, IsFingerprintsFileDataValid,
IsHeaderDataKeyPresent, Next, Read, SetBitStringFormat, SetBitsOrder,
SetCompoundID, SetDetailLevel, SetFingerprints, SetFingerprintsString,
SetFingerprintsStringMode, SetPartialFingerprintsString, SetVectorStringFormat,
WriteFingerprints, WriteFingerprintsString

The following methods can also be used as functions:

IsFingerprintsFPFile

B<FingerprintsFPFileIO> class is derived from I<FileIO> class and uses its methods to support
generic file related functionality.

The MayaChemTools fingerprints file (FP) format with B<.fpf> or B<.fp> file extensions supports
two types of fingerprints data: fingerprints bit-vectors and fingerprints vectors.

Example of FP file format containing fingerprints bit-vector string data:

    #
    # Package = MayaChemTools 7.4
    # ReleaseDate = Oct 21, 2010
    #
    # TimeStamp =  Mon Mar 7 15:14:01 2011
    #
    # FingerprintsStringType = FingerprintsBitVector
    #
    # Description = PathLengthBits:AtomicInvariantsAtomTypes:MinLength1:...
    # Size = 1024
    # BitStringFormat = HexadecimalString
    # BitsOrder = Ascending
    #
    Cmpd1 9c8460989ec8a49913991a6603130b0a19e8051c89184414953800cc21510...
    Cmpd2 000000249400840040100042011001001980410c000000001010088001120...
    ... ...
    ... ..

Example of FP file format containing fingerprints vector string data:

    #
    # Package = MayaChemTools 7.4
    # ReleaseDate = Oct 21, 2010
    #
    # TimeStamp =  Mon Mar 7 15:14:01 2011
    #
    # FingerprintsStringType = FingerprintsVector
    #
    # Description = PathLengthBits:AtomicInvariantsAtomTypes:MinLength1:...
    # VectorStringFormat = IDsAndValuesString
    # VectorValuesType = NumericalValues
    #
    Cmpd1 338;C F N O C:C C:N C=O CC CF CN CO C:C:C C:C:N C:CC C:CF C:CN C:
    N:C C:NC CC:N CC=O CCC CCN CCO CNC NC=O O=CO C:C:C:C C:C:C:N C:C:CC...;
    33 1 2 5 21 2 2 12 1 3 3 20 2 10 2 2 1 2 2 2 8 2 5 1 1 1 19 2 8 2 2 2 2
    6 2 2 2 2 2 2 2 2 3 2 2 1 4 1 5 1 1 18 6 2 2 1 2 10 2 1 2 1 2 2 2 2 ...
    Cmpd2 103;C N O C=N C=O CC CN CO CC=O CCC CCN CCO CNC N=CN NC=O NCN O=C
    O C CC=O CCCC CCCN CCCO CCNC CNC=N CNC=O CNCN CCCC=O CCCCC CCCCN CC...;
    15 4 4 1 2 13 5 2 2 15 5 3 2 2 1 1 1 2 17 7 6 5 1 1 1 2 15 8 5 7 2 2 2 2
    1 2 1 1 3 15 7 6 8 3 4 4 3 2 2 1 2 3 14 2 4 7 4 4 4 4 1 1 1 2 1 1 1 ...
    ... ...
    ... ...

B<FP> file data format consists of two main sections: header section and fingerprints string
data section. The header section lines start with # and the first line not starting with # represents
the start of fingerprints string data section. The header section contains both the required and
optional information which is specified as key = value pairs. The required information
describes fingerprints bit-vector and vector strings and used to generate fingerprints objects;
the optional information is ignored during generation of fingerpints objects.

The key = value data specification in the header section and its processing follows these
rules:

    o Leading and trailing spaces for key = value pairs are ignored
    o Key and value strings may contain spaces
    o Multiple key = value pairs on a single are delimited by semicolon

The default optional header data section key = value pairs are:

    # Package = MayaChemTools 7.4
    # ReleaseDate = Oct 21, 2010

The B<FingerprintsStringType> key is required data header key for both fingerprints bit-vector
and vector strings. Possible key values: I<FingerprintsBitVector or FingerprintsVector>.
For example:

    # FingerprintsStringType = FingerprintsBitVector

The required data header keys for fingerprints bit-vector string are: B<Description, Size,
BitStringFormat, and BitsOrder>. Possible values for B<BitStringFormat>: I<HexadecimalString
or BinaryString>. Possible values for B<BitsOrder>: I<Ascending or Descending>. The B<Description>
key contains information about various parameters used to generate fingerprints bit-vector
string. The B<Size> corresponds to number of fingerprints bits and is always less than or equal
to number of bits in bit-vetor string which might contain extra bits at the end to round off the
size to make it multiple of 8. For example:

    # Description = PathLengthBits:AtomicInvariantsAtomTypes:MinLength1:...
    # Size = 1024
    # BitStringFormat = HexadecimalString
    # BitsOrder = Ascending

The required data header keys for fingerprints vector string are: B<Description, VectorStringFormat,
and VectorValuesType>.  Possible values for B<VectorStringFormat>: I<DsAndValuesString,
IDsAndValuesPairsString, ValuesAndIDsString, ValuesAndIDsPairsString or ValuesString>.
Possible values for B<VectorValuesType>: I<NumericalValues, OrderedNumericalValues or
AlphaNumericalValues>. The B<Description> keys contains information various parameters used
to generate fingerprints vector string. For example:

    # Description = PathLengthBits:AtomicInvariantsAtomTypes:MinLength1:...
    # VectorStringFormat = IDsAndValuesString
    # VectorValuesType = NumericalValues

The fingerprints data section for fingerprints bit-vector string contains data in the following
format:

    ... ...
    CmpdID FingerprintsPartialBitVectorString
    ... ...

For example:

    ... ...
    Cmpd1 9c8460989ec8a49913991a6603130b0a19e8051c89184414953800cc21510...
    ... ...

The fingerprints data section for fingerprints vector string contains data in the following
format:

    ... ...
    CmpdID Size;FingerprintsPartialVectorString
    ... ...

For example:

    ... ...
    Cmpd1 338;C F N O C:C C:N C=O CC CF CN CO C:C:C C:C:N C:CC C:CF C:CN C:
    N:C C:NC CC:N CC=O CCC CCN CCO CNC NC=O O=CO C:C:C:C C:C:C:N C:C:CC...;
    33 1 2 5 21 2 2 12 1 3 3 20 2 10 2 2 1 2 2 2 8 2 5 1 1 1 19 2 8 2 2 2 2
    6 2 2 2 2 2 2 2 2 3 2 2 1 4 1 5 1 1 18 6 2 2 1 2 10 2 1 2 1 2 2 2 2 ...
    ... ...

Unlike fingerprints bit-vector string, I<Size> is specified for each partial fingerprints vector string:
It may change from molecule to molecule for same type of fingerprints.

Values IDs are optional for fingerprints vector string containing I<OrderedNumericalValues or
AlphaNumericalValues>; however, they must be present for for I<NumericalValues>. Due to
various possible values for B<VectorStringFormat>, the fingerprints data section for fingerprints
vector string supports following type of data formats:

    CmpdID Size;ID1 ID2 ID3...;Value1 Value2 Value3...
    CmpdID Size;ID1 Value1 ID2 Value2 ID3 Value3... ...
    CmpdID Size;ValuesAndIDsString: Value1 Value2 Value3...;ID1 ID2 ID3...
    CmpdID Size;ValuesAndIDsPairsString: Value1 ID1 Value2 ID2 Value3 ID3... ...
    CmpdID Size;Value1 Value2 Value3 ...

However, all the fingerprints vector string data present in FP file must correspond to only
one of the formats shown above; multiple data formats in the same file are not allowed.

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

    $NewFingerprintsFPFileIO = new FileIO::FingerprintsFPFileIO(%IOParameters);

Using specified I<IOParameters> names and values hash, B<new> method creates a new
object and returns a reference to a newly created B<FingerprintsFPFileIO> object. By default,
the following properties are initialized during I<Read> mode:

    Name = '';
    Mode = 'Read';
    Status = 0;
    FingerprintsStringMode = 'AutoDetect';
    ValidateData = 1;
    DetailLevel = 1;

During I<Write> mode, the following properties get initialize by default:

    FingerprintsStringMode = undef;

    BitStringFormat = HexadecimalString;
    BitsOrder = Ascending;

    VectorStringFormat = NumericalValuesString or ValuesString;

Examples:

    $NewFingerprintsFPFileIO = new FileIO::FingerprintsFPFileIO(
                               'Name' => 'Sample.fpf',
                               'Mode' => 'Read',
                               'FingerprintsStringMode' =>
                                       'AutoDetect');

    $NewFingerprintsFPFileIO = new FileIO::FingerprintsFPFileIO(
                               'Name' => 'Sample.fpf',
                               'Mode' => 'Write',
                               'FingerprintsStringMode' =>
                                       'FingerprintsBitVectorString',
                               'Overwrite' => 1,
                               'BitStringFormat' => 'HexadecimalString',
                               'BitsOrder' => 'Ascending');

    $NewFingerprintsFPFileIO = new FileIO::FingerprintsFPFileIO(
                               'Name' => 'Sample.fp',
                               'Mode' => 'Write',
                               'FingerprintsStringMode' =>
                                       'FingerprintsVectorString',
                               'Overwrite' => 1,
                               'VectorStringFormat' => 'IDsAndValuesString');

=item B<GetFingerprints>

    $FingerprintsObject = $FingerprintsFPFileIO->GetFingerprints();

Returns B<FingerprintsObject> generated for current data line using fingerprints bit-vector
or vector string data. The fingerprints object corresponds to any of the supported fingerprints
such as PathLengthFingerprints, ExtendedConnectivity, and so on.

=item B<GetFingerprintsString>

    $FingerprintsString = $FingerprintsFPFileIO->GetFingerprintsString();

Returns B<FingerprintsString> for current data line.

=item B<GetHeaderDataKeyValue>

    $KeyValue = $FingerprintsFPFileIO->GetHeaderDataKeyValue($Key);

Returns B<KeyValue> of a data header I<Key>.

=item B<GetHeaderDataKeys>

    @Keys = $FingerprintsFPFileIO->GetHeaderDataKeys();
    $NumOfKeys = $FingerprintsFPFileIO->GetHeaderDataKeys();

Returns an array of data header B<Keys> retrieved from data header section of fingerprints
file. In scalar context, it returns number of keys.

=item B<GetHeaderDataKeysAndValues>

    %KeysAndValues = $FingerprintsFPFileIO->GetHeaderDataKeysAndValues();

Returns a hash of data header keys and values retrieved from data header section of fingerprints
file.

=item B<GetPartialFingerprintsString>

    $FingerprintsString = $FingerprintsFPFileIO->GetPartialFingerprintsString();

Returns partial B<FingerprintsString> for current data line. It corresponds to fingerprints string
specified present in a line.

=item B<GetRequiredHeaderDataKeys>

    @Keys = $FingerprintsFPFileIO->GetRequiredHeaderDataKeys();
    $NumOfKeys = $FingerprintsFPFileIO->GetRequiredHeaderDataKeys();

Returns an array of required data header B<Keys> for a fingerprints file containing bit-vector or
vector strings data. In scalar context, it returns number of keys.

=item B<GetRequiredHeaderDataKeysAndValues>

    %KeysAndValues = $FingerprintsFPFileIO->
                     GetRequiredHeaderDataKeysAndValues();

Returns a hash of required data header keys and values for a fingerprints file containing bit-vector or
vector strings data

=item B<IsFingerprintsDataValid>

    $Status = $FingerprintsFPFileIO->IsFingerprintsDataValid();

Returns 1 or 0 based on whether B<FingerprintsObject> is valid.

=item B<IsFingerprintsFPFile>

    $Status = $FingerprintsFPFileIO->IsFingerprintsFPFile($FileName);
    $Status = FileIO::FingerprintsFPFileIO::IsFingerprintsFPFile($FileName);

Returns 1 or 0 based on whether I<FileName> is a FP file.

=item B<IsFingerprintsFileDataValid>

    $Status = $FingerprintsFPFileIO->IsFingerprintsFileDataValid();

Returns 1 or 0 based on whether fingerprints file contains valid fingerprints data.

=item B<IsHeaderDataKeyPresent>

    $Status = $FingerprintsFPFileIO->IsHeaderDataKeyPresent($Key);

Returns 1 or 0 based on whether data header I<Key> is present in data header
section of a FP file.

=item B<Next or Read>

    $FingerprintsFPFileIO = $FingerprintsFPFileIO->Next();
    $FingerprintsFPFileIO = $FingerprintsFPFileIO->Read();

Reads next available fingerprints line in FP file, processes the data, generates appropriate fingerprints
object, and returns B<FingerprintsFPFileIO>. The generated fingerprints object is available using
method B<GetFingerprints>.

=item B<SetBitStringFormat>

    $FingerprintsFPFileIO->SetBitStringFormat($Format);

Sets bit string I<Format> for fingerprints bit-vector string data in a FP file and returns B<FingerprintsFPFileIO>.
Possible values for B<BitStringFormat>: I<BinaryString or HexadecimalString>.

=item B<SetBitsOrder>

    $FingerprintsFPFileIO->SetBitsOrder($BitsOrder);

Sets I<BitsOrder> for fingerprints bit-vector string data in a FP file and returns B<FingerprintsFPFileIO>.
Possible values for B<BitsOrder>: I<Ascending or Descending>.

=item B<SetCompoundID>

    $FingerprintsFPFileIO->SetCompoundID($ID);

Sets compound ID for current data line and returns B<FingerprintsFPFileIO>. Spaces are not allowed
in compound IDs.

=item B<SetDetailLevel>

    $FingerprintsFPFileIO->SetDetailLevel($Level);

Sets details I<Level> for generating diagnostics messages during FP file processing and returns
B<FingerprintsFPFileIO>. Possible values: I<Positive integers>.

=item B<SetFingerprints>

    $FingerprintsFPFileIO->SetFingerprints($FingerprintsObject);

Sets I<FingerprintsObject> for current data line and returns B<FingerprintsFPFileIO>.

=item B<SetFingerprintsString>

    $FingerprintsFPFileIO->SetFingerprintsString($FingerprintsString);

Sets I<FingerprintsString> for current data line and returns B<FingerprintsFPFileIO>.

=item B<SetFingerprintsStringMode>

    $FingerprintsFPFileIO->SetFingerprintsStringMode($Mode);

Sets I<FingerprintsStringMode> for FP file and returns B<FingerprintsFPFileIO>.
Possible values: I<AutoDetect, FingerprintsBitVectorString, or FingerprintsVectorString>

=item B<SetPartialFingerprintsString>

    $FingerprintsFPFileIO->SetPartialFingerprintsString($PartialString);

Sets I<PartialFingerprintsString> for current data line and returns B<FingerprintsFPFileIO>.

=item B<SetVectorStringFormat>

    $FingerprintsFPFileIO->SetVectorStringFormat($Format);

Sets I<VectorStringFormat> for FP file and returns B<FingerprintsFPFileIO>. Possible values:
I<IDsAndValuesString, IDsAndValuesPairsString, ValuesAndIDsString, ValuesAndIDsPairsString>.

=item B<WriteFingerprints>

    $FingerprintsFPFileIO->WriteFingerprints($FingerprintsObject,
                                            $CompoundID);

Writes fingerprints string generated from I<FingerprintsObject> object and other data including
I<CompoundID> to FP file and returns B<FingerprintsFPFileIO>.

=item B<WriteFingerprintsString>

    $FingerprintsFPFileIO->WriteFingerprints($FingerprintsString,
                                            $CompoundID);

Writes I<FingerprintsString> and other data including I<CompoundID> to FP file and returns
B<FingerprintsFPFileIO>.

Caveats:

    o FingerprintsStringMode, BitStringFormat, BitsOrder, VectorStringFormat
      values are ignored during writing of fingerprints and it's written to
      the file as it is.
    o FingerprintsString is a regular fingerprints string as oppose to a
      partial fingerprints string.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

FingerprintsSDFileIO.pm, FingerprintsTextFileIO.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
