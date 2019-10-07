package NucleicAcids;
#
# File: NucleicAcids.pm
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
use Text::ParseWords;
use TextUtil;
use FileUtil;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GetNucleicAcids GetNucleicAcidsByType GetNucleicAcidPropertiesData GetNucleicAcidPropertiesNames IsNucleicAcid IsNucleicAcidProperty IsNucleicAcidType);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

#
# Load nucleic acids data...
#
my(%NucleicAcidDataMap, %NucleicAcidCodeMap, %NucleicAcidOtherCodeMap, %NucleicAcidNameMap, @NucleicAcidCodes, @NucleicAcidPropertyNames, %NucleicAcidPropertyNamesMap, %NucleicAcidTypesMap);
_LoadNucleicAcidsData();

#
# Get a list of all known nucleic acids as one of these values:
# code or nucleic acid name...
#
sub GetNucleicAcids {
  my($NameType, $Code, $Name, @NucleicAcidNames);

  $NameType = 'Code';
  if (@_ >= 1) {
    ($NameType) = @_;
  }

  # Collect names...
  @NucleicAcidNames = ();
  for $Code (@NucleicAcidCodes) {
    NAME : {
      if ($NameType =~ /^Name$/i) {$Name = $NucleicAcidDataMap{$Code}{Name}; last NAME; }
      $Name = $Code;
    }
    push @NucleicAcidNames, $Name;
  }

  return (wantarray ? @NucleicAcidNames : \@NucleicAcidNames);
}

#
# Get a list of all known nucleic acids by one of these specified types:
# Nucleobase, Nucleoside, Deoxynucleoside, Nucleotide, Deoxynucleotide. Default: Nucleoside
#
sub GetNucleicAcidsByType {
  my($NameType, $Type, $Code, $Name, @NucleicAcidNames);

  $Type = 'Nucleoside';
  $NameType = 'Code';
  if (@_ == 2) {
    ($Type, $NameType) = @_;
  }
  elsif (@_ == 1) {
    ($Type) = @_;
  }

  # Collect names...
  @NucleicAcidNames = ();
  CODE: for $Code (@NucleicAcidCodes) {
    if ($NucleicAcidDataMap{$Code}{Type} !~ /^$Type$/i ) {
      next CODE;
    }
    NAME : {
      if ($NameType =~ /^Name$/i) {$Name = $NucleicAcidDataMap{$Code}{Name}; last NAME; }
      $Name = $Code;
    }
    push @NucleicAcidNames, $Name;
  }

  return (wantarray ? @NucleicAcidNames : \@NucleicAcidNames);
}

#
# Get all available properties data for an nucleic acid using any of these symbols:
# code, other code or name.
#
# A reference to a hash array is returned with keys and values representing property
# name and its values respectively.
#
sub GetNucleicAcidPropertiesData {
  my($NucleicAcidID) = @_;
  my($Code);

  if ($Code = _ValidateNucleicAcidID($NucleicAcidID)) {
    return \%{$NucleicAcidDataMap{$Code}};
  }
  else {
    return undef;
  }
}

#
# Get names of all available nucleic acid properties. A reference to  an array containing
# names of all available properties is returned.
#
sub GetNucleicAcidPropertiesNames {
  my($Mode);
  my($PropertyName, @PropertyNames);

  $Mode = 'ByGroup';
  if (@_ == 1) {
    ($Mode) = @_;
  }

  @PropertyNames = ();
  if ($Mode =~ /^Alphabetical$/i) {
    my($PropertyName);
    # Code, OtherCodes and Name are always listed first...
    push @PropertyNames, qw(Code OtherCodes Name);
    for $PropertyName (sort keys %NucleicAcidPropertyNamesMap) {
      if ($PropertyName !~ /^(Code|OtherCodes|Name)$/) {
	push @PropertyNames, $PropertyName;
      }
    }
  }
  else {
    push @PropertyNames, @NucleicAcidPropertyNames;
  }
  return (wantarray ? @PropertyNames : \@PropertyNames);
}

#
# Is it a known nucleic acid? Input is either a code or a name
#
sub IsNucleicAcid {
  my($NucleicAcidID) = @_;
  my($Status);

  $Status = (_ValidateNucleicAcidID($NucleicAcidID)) ? 1 : 0;

  return $Status;
}

#
# Is it an available nucleic acid property?
#
sub IsNucleicAcidProperty {
  my($PropertyName) = @_;
  my($Status);

  $Status = (exists($NucleicAcidPropertyNamesMap{$PropertyName})) ? 1 : 0;

  return $Status;
}

#
# Is it an available nucleic acid type?
#
sub IsNucleicAcidType {
  my($Type) = @_;
  my($Status);

  $Status = (exists($NucleicAcidTypesMap{lc($Type)})) ? 1 : 0;

  return $Status;
}

#
# Implents GetNucleicAcid<PropertyName> for a valid proprty name.
#
sub AUTOLOAD {
  my($NucleicAcidID) = @_;
  my($FunctionName, $PropertyName, $PropertyValue, $Code);

  $PropertyValue = undef;

  use vars qw($AUTOLOAD);
  $FunctionName = $AUTOLOAD;
  $FunctionName =~ s/.*:://;

  # Only Get<PropertyName> functions are supported...
  if ($FunctionName !~ /^Get/) {
    croak "Error: Function, NucleicAcid::$FunctionName, is not supported by AUTOLOAD in NucleicAcid module: Only Get<PropertyName> functions are implemented...";
  }

  $PropertyName = $FunctionName;
  $PropertyName =~  s/^GetNucleicAcid//;
  if (!exists $NucleicAcidPropertyNamesMap{$PropertyName}) {
    croak "Error: Function, NucleicAcid::$FunctionName, is not supported by AUTOLOAD in NucleicAcid module: Unknown nucleic acid property name, $PropertyName, specified...";
  }

  if (!($Code = _ValidateNucleicAcidID($NucleicAcidID))) {
    return undef;
  }
  $PropertyValue = $NucleicAcidDataMap{$Code}{$PropertyName};
  return $PropertyValue;
}

#
# Load NucleicAcidsData.csv files from <MayaChemTools>/lib directory...
#
sub _LoadNucleicAcidsData {
  my($NucleicAcidsDataFile, $MayaChemToolsLibDir);

  $MayaChemToolsLibDir = GetMayaChemToolsLibDirName();

  $NucleicAcidsDataFile =  "$MayaChemToolsLibDir" . "/data/NucleicAcidsData.csv";

  if (! -e "$NucleicAcidsDataFile") {
    croak "Error: MayaChemTools package file, $NucleicAcidsDataFile, is missing: Possible installation problems...";
  }

  _LoadData($NucleicAcidsDataFile);
}

#
# Load NucleicAcidsData.csv file from <MayaChemTools>/lib directory...
#
sub _LoadData {
  my($NucleicAcidsDataFile) = @_;

  %NucleicAcidDataMap = ();
  @NucleicAcidCodes = ();
  @NucleicAcidPropertyNames = ();
  %NucleicAcidPropertyNamesMap = ();
  %NucleicAcidCodeMap = ();
  %NucleicAcidOtherCodeMap = ();
  %NucleicAcidNameMap = ();
  %NucleicAcidTypesMap = ();

  # Load property data for all nucleic acids...
  #
  # File Format:
  # "Code","OtherCodes","BasePair","Name","Type","ChemicalFormula","ChemicalFormulaAtpH7.5","MolecularWeight","ExactMass","ElementalComposition"
  #
  my($Code, $OtherCodes, $NucleicAcidName, $Line, $NumOfCols, $InDelim, $Index, $Name, $Value, $Units, @LineWords, @ColLabels);

  $InDelim = "\,";
  open NUCLEICACIDSDATAFILE, "$NucleicAcidsDataFile" or croak "Couldn't open $NucleicAcidsDataFile: $! ...";

  # Skip lines up to column labels...
  LINE: while ($Line = GetTextLine(\*NUCLEICACIDSDATAFILE)) {
    if ($Line !~ /^#/) {
      last LINE;
    }
  }
  @ColLabels= quotewords($InDelim, 0, $Line);
  $NumOfCols = @ColLabels;

  # Extract property names from column labels...
  @NucleicAcidPropertyNames = ();
  for $Index (0 .. $#ColLabels) {
    $Name = $ColLabels[$Index];
    push @NucleicAcidPropertyNames, $Name;

    # Store property names...
    $NucleicAcidPropertyNamesMap{$Name} = $Name;
  }

  # Process nucleic acid data...
  LINE: while ($Line = GetTextLine(\*NUCLEICACIDSDATAFILE)) {
    if ($Line =~ /^#/) {
      next LINE;
    }
    @LineWords = ();
    @LineWords = quotewords($InDelim, 0, $Line);
    if (@LineWords != $NumOfCols) {
      croak "Error: The number of data fields, @LineWords, in $NucleicAcidsDataFile must be $NumOfCols.\nLine: $Line...";
    }
    $Code = $LineWords[0]; $OtherCodes = $LineWords[1]; $NucleicAcidName = $LineWords[3];
    if (exists $NucleicAcidDataMap{$Code}) {
      carp "Warning: Ignoring data for nucleic acid $Code: It has already been loaded.\nLine: $Line....";
      next LINE;
    }

    # Store all the values...
    push @NucleicAcidCodes, $Code;
    %{$NucleicAcidDataMap{$Code}} = ();
    for $Index (0 .. $#LineWords) {
      $Name = $NucleicAcidPropertyNames[$Index];
      $Value = $LineWords[$Index];
      $NucleicAcidDataMap{$Code}{$Name} = $Value;
    }
  }
  close NUCLEICACIDSDATAFILE;

  # Setup one letter and nucleic acid name maps...
  _SetupNucleicAcidIDMap();
}

#
# Setup lowercase other codes and name maps pointing
# to code as show in data file.
#
sub _SetupNucleicAcidIDMap {
  my($Code, @OtherCodes, $OtherCode, $NucleicAcidName, $NucleicAcidType);

  %NucleicAcidCodeMap = ();
  %NucleicAcidOtherCodeMap = ();
  %NucleicAcidNameMap = ();
  %NucleicAcidTypesMap = ();

  for $Code (keys %NucleicAcidDataMap) {
    $NucleicAcidCodeMap{lc($Code)} = $Code;

    $NucleicAcidName = $NucleicAcidDataMap{$Code}{Name};
    $NucleicAcidNameMap{lc($NucleicAcidName)} = $Code;

    $NucleicAcidType = $NucleicAcidDataMap{$Code}{Type};
    if (! exists $NucleicAcidTypesMap{$NucleicAcidType}) {
      $NucleicAcidTypesMap{lc($NucleicAcidType)} = $NucleicAcidType;
    }

    @OtherCodes = split /\,/, $NucleicAcidDataMap{$Code}{OtherCodes};
    OTHERCODE: for $OtherCode (@OtherCodes) {
      if (!$OtherCode) {
	next OTHERCODE;
      }
      $OtherCode = RemoveLeadingAndTrailingWhiteSpaces($OtherCode);
      $NucleicAcidOtherCodeMap{lc($OtherCode)} = $Code;
    }
  }
}

# Validate Nucleic acid ID...
sub _ValidateNucleicAcidID {
  my($NucleicAcidID) = @_;
  my($Code) = undef;

  if (exists $NucleicAcidCodeMap{lc($NucleicAcidID)}) {
    $Code = $NucleicAcidCodeMap{lc($NucleicAcidID)};
  }
  elsif (exists $NucleicAcidOtherCodeMap{lc($NucleicAcidID)}) {
    $Code = $NucleicAcidOtherCodeMap{lc($NucleicAcidID)};
  }
  elsif (exists $NucleicAcidNameMap{lc($NucleicAcidID)}) {
    $Code = $NucleicAcidNameMap{lc($NucleicAcidID)};
  }
  return $Code;
}


1;

__END__

=head1 NAME

NucleicAcids

=head1 SYNOPSIS

use NucleicAcids;

use NucleicAcids qw(:all);

=head1 DESCRIPTION

B<NucleicAcids> module the provides the following functions:

GetNucleicAcidPropertiesData, GetNucleicAcidPropertiesNames,
GetNucleicAcids, GetNucleicAcidsByType, IsNucleicAcid, IsNucleicAcidProperty,
IsNucleicAcidType

=head1 Functions

=over 4

=item B<GetNucleicAcids>

    (@Names) = GetNucleicAcids([$NameType]);
    $NamesRef = GetNucleicAcids([$NameType]);

Returns an array or a reference to an array containing names of nucleic acids
as a code or nucleic acid name controlled by optional parameter I<NameType>. By
default, nucleic acids names are returned as the code. Possible values for
I<NameType>: I<Code or Name>.

=item B<GetNucleicAcidsByType>

    (@Names) = GetNucleicAcidsByType([$Type, $NameType]);
    $NamesRef = GetNucleicAcidsByType([$Type, $NameType]);

Returns an array or a reference to an array containing names of nucleic acids
specified by parameter I<Type> as a code or name controlled by optional
parameter I<NameType>. Default values for I<Type>: I<Nucleoside>. Default value for
I<NameType>: I<Code>. Possible values for I<Type>: I<Nucleobase, Nucleoside, Deoxynucleoside,
Nucleotide, Deoxynucleotide>. Possible values for I<NameType>: I<Code or Name>.

=item B<GetNucleicAcidPropertiesData>

    $DataHashRef = GetNucleicAcidPropertiesData($NucleicAcidID);

Returns a reference to hash containing property names and values for a specified
I<NucleicAcidID>.

=item B<GetNucleicAcidPropertyName>

    $Value = GetNucleicAcid<PropertyName>($NucleicAcidID);

Returns nucleic acid property value for a specified I<NucleicAcidID>. This function is
implemented on-the-fly using Perl's AUTOLOAD functionality.

=item B<GetNucleicAcidPropertiesNames>

    @Names = GetNucleicAcidPropertiesNames([$Mode]);
    $NamesRef = GetNucleicAcidPropertiesNames([$Mode]);

Returns an array or a reference to an array containing names of properties for
nucleic acids. Order of nucleic acids properties is controlled by optional parameter
I<Mode>. Possible values for I<Mode>: I<Alphabetical or ByGroup>; Default: I<ByGroup>.

=item B<IsNucleicAcid>

    $Status = IsNucleicAcid($NucleicAcidID);

Returns 1 or 0 based on whether it's a known nucleic acid ID.

=item B<IsNucleicAcidProperty>

    $Status = IsNucleicAcid($PropertyName);

Returns 1 or 0 based on whether it's a known nucleic acid property name.

=item B<IsNucleicAcidType>

    $Status = IsNucleicAcidType();

Returns 1 or 0 based on whether it's a known nucleic acid type.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AminoAcids.pm, PeriodicTable.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
