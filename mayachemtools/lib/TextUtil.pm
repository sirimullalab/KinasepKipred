package TextUtil;
#
# File: TextUtil.pm
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
use Exporter;
use Text::ParseWords;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw(AddNumberSuffix ContainsWhiteSpaces GetTextLine GetTextFileDataByUniqueKey GetTextFileDataByNonUniqueKey HashCode IsEmpty IsNumberPowerOfNumber IsInteger IsPositiveInteger IsFloat IsNotEmpty IsNumerical JoinWords SplitWords  QuoteAWord RemoveLeadingWhiteSpaces RemoveTrailingWhiteSpaces RemoveLeadingAndTrailingWhiteSpaces WrapText);
@EXPORT_OK = qw();
%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Add number suffix...
sub AddNumberSuffix {
  my($Value) = @_;
  my($ValueWithSuffix, $Suffix);

  $ValueWithSuffix = $Value;
  if (!IsPositiveInteger($Value)) {
    return $ValueWithSuffix;
  }
  $Suffix = "th";
  if ($Value < 10 || $Value > 20) {
    my $Remainder = $Value % 10;
    $Suffix = ($Remainder == 1) ? "st" : (($Remainder == 2) ? "nd" : (($Remainder == 3) ? "rd" : "th"));
  }
  $ValueWithSuffix = "${ValueWithSuffix}${Suffix}";
  return $ValueWithSuffix;
}

# Check out the string: Doen it contain any white space characters?
sub ContainsWhiteSpaces {
  my($TheString) = @_;
  my($Status) = 0;

  if (defined($TheString) && length($TheString)) {
    $Status = ($TheString =~ /[ \t\r\n\f]/ ) ? 1 : 0;
  }
  return $Status;
}

# Read the line, change to UNIX new line char, and chop off new line char as well...
sub GetTextLine {
  my($TextFileRef) = @_;
  my($Line) = '';

  # Get the next non empty line...
  LINE: while (defined($_ = <$TextFileRef>)) {
    # Change Windows and Mac new line char to UNIX...
    s/(\r\n)|(\r)/\n/g;

    # Take out any new line char at the end by explicitly removing it instead of using
    # chomp, which might not always work correctly on files generated on a system
    # with a value of input line separator different from the current system...
    s/\n$//g;

    # Doesn't hurt to chomp...
    chomp;

    $Line = $_;
    if (length $Line) {
      last LINE;
    }
  }
  return $Line;
}

# Load data from a CSV file into the specified hash reference using a specific
# column for unique data key values.
#
# The lines starting with # are treated as comments and ignored. First line
# not starting with # must contain column labels and the number of columns in
# all other data rows must match the number of column labels.
#
# The first column is assumed to contain data key value by default; all other columns
# contain data as indicated in their column labels.
#
# In order to avoid dependence of data access on the specified column labels, the
# column data is loaded into hash with Column<Num> hash keys, where column number
# start from 1. The data key column is not available as Colnum<Num> hash key;
#
# The format of the data structure loaded into a specified hash reference is:
#
# @{$TextDataMapRef->{DataKeys}} - Array of unique data keys
# @{$TextDataMapRef->{ColLabels}} - Array of column labels
# @{$TextDataMapRef->{DataColIDs}} - Array of data column IDs
# $TextDataMapRef->{NumOfCols} - Number of columns
# %{$TextDataMapRef->{DataKey}} - Hash keys pair: <DataKey, DataKey>
# %{$TextDataMapRef->{DataCol<Num>}} - Hash keys pair: <DataCol<Num>, DataKey>
#
# Caveats:
#   . The column number start from 1.
#   . Column data for data key column column is not loaded into <Column<Num>, DataKey> hash keys pairs.
#
sub GetTextFileDataByUniqueKey {
  my($TextDataFile, $TextDataMapRef, $DataKeyColNum, $InDelim) = @_;

  return _GetTextFileData("UniqueKey", $TextDataFile, $TextDataMapRef, $DataKeyColNum, $InDelim);
}

# Load data from a CSV file into the specified hash reference using a specific
# column for non-unique data key values.
#
# The lines starting with # are treated as comments and ignored. First line
# not starting with # must contain column labels and the number of columns in
# all other data rows must match the number of column labels.
#
# The first column is assumed to contain data key value by default; all other columns
# contain data as indicated in their column labels.
#
# In order to avoid dependence of data access on the specified column labels, the
# column data is loaded into hash with Column<Num> hash keys, where column number
# start from 1. The data key column is not available as Colnum<Num> hash key;
#
# The format of the data structure loaded into a specified hash reference is:
#
# @{$TextDataMapRef->{DataKeys}} - Array of unique data keys
# @{$TextDataMapRef->{ColLabels}} - Array of column labels
# @{$TextDataMapRef->{DataColIDs}} - Array of data column IDs
# $TextDataMapRef->{NumOfCols} - Number of columns
# %{$TextDataMapRef->{DataKey}} - Hash keys pair: <DataKey, DataKey>
# @{$TextDataMapRef->{DataCol<Num>}} - Hash keys pair with data as an array: <DataCol<Num>, DataKey>
#
# Caveats:
#   . The column number start from 1.
#   . Column data for data key column column is not loaded into <Column<Num>, DataKey> hash keys pairs.
#
sub GetTextFileDataByNonUniqueKey {
  my($TextDataFile, $TextDataMapRef, $DataKeyColNum, $InDelim) = @_;

  return _GetTextFileData("NonUniqueKey", $TextDataFile, $TextDataMapRef, $DataKeyColNum, $InDelim);
}

# Loadtext file data using unique or non-uniqye data column key...
#
sub _GetTextFileData {
  my($DataKeyMode, $TextDataFile, $TextDataMapRef, $DataKeyColNum, $InDelim) = @_;
  my($DataKeyColIndex, $LineCount, $IgnoredLineCount, $UniqueDataKeyMode, $DataKey, $Line, $NumOfCols, $ColIndex, $ColNum, $ColID, $ColValue, @LineWords, @ColLabels, @DataColIDs, @DataColNums);

  print "\nProcessing text data file $TextDataFile...\n";

  $UniqueDataKeyMode = 0;
  if ($DataKeyMode =~ /^UniqueKey$/i) {
    $UniqueDataKeyMode = 1;
  }

  # Setup default values...
  $DataKeyColNum = defined  $DataKeyColNum ? $DataKeyColNum : 1;

  if (defined $InDelim) {
    if ($InDelim =~ /^tab$/i) {
      $InDelim = "\t";
    }
    elsif ($InDelim =~ /^semicolon$/i) {
      $InDelim = "\;";
    }
    elsif ($InDelim =~ /^comma$/i) {
      $InDelim = "\,";
    }
    else {
      warn "Warning: Ignoring specified input delimiter: $InDelim. Supported values: comma, semicolon or tab. Using default comma delimiter...";
      $InDelim = "\,";
    }
  }
  else {
    if ($TextDataFile =~ /\.tsv$/i) {
      $InDelim = "\t";
    }
    elsif ($TextDataFile =~ /\.csv$/i) {
      $InDelim = "\,";
    }
    else {
      warn "Warning: Unknown file extension. Using default comma delimiter...";
      $InDelim = "\,";
    }
  }

  ($LineCount, $IgnoredLineCount) = (0) x 2;

  open TEXTDATAFILE, "$TextDataFile" or die "Couldn't open $TextDataFile: $! ...";

  # Skip lines up to column labels...
  LINE: while ($Line = TextUtil::GetTextLine(\*TEXTDATAFILE)) {
    $LineCount++;
    if ($Line =~ /^#/) {
      $IgnoredLineCount++;
    }
    else {
      last LINE;
    }
  }

  # Initialize data map...
  %{$TextDataMapRef} = ();
  @{$TextDataMapRef->{DataKeys}} = ();
  @{$TextDataMapRef->{ColLabels}} = ();
  @{$TextDataMapRef->{DataColIDs}} = ();
  $TextDataMapRef->{NumOfCols} = undef;

  # Process column labels...
  @ColLabels= quotewords($InDelim, 0, $Line);
  $NumOfCols = @ColLabels;

  if ($DataKeyColNum < 1 || $DataKeyColNum > $NumOfCols) {
    warn "Warning: Ignoring text data file $TextDataFile: Invalid data key column number, $DataKeyColNum, specified. It must be > 0 or <= $NumOfCols, number of columns in the text file ...";
    return;
  }
  $DataKeyColIndex = $DataKeyColNum - 1;

  $TextDataMapRef->{NumOfCols} = $NumOfCols;
  push @{$TextDataMapRef->{ColLabels}}, @ColLabels;

  # Set up column data IDs for tracking the data...
  @DataColNums = ();
  @DataColIDs = ();
  COLNUM: for $ColNum (1 .. $NumOfCols) {
    if ($ColNum == $DataKeyColNum) {
      next COLNUM;
    }
    push @DataColNums, $ColNum;
    $ColID = "DataCol${ColNum}";
    push @DataColIDs, $ColID;
  }
  push @{$TextDataMapRef->{DataColIDs}}, @DataColIDs;

  # Initialize column data hash...
  %{$TextDataMapRef->{DataKey}} = ();
  for $ColIndex (0 .. $#DataColNums) {
    $ColNum = $DataColNums[$ColIndex];
    $ColID = $DataColIDs[$ColIndex];
    %{$TextDataMapRef->{$ColID}} = ();
  }

  LINE: while ($Line = TextUtil::GetTextLine(\*TEXTDATAFILE)) {
    $LineCount++;
    if ($Line =~ /^#/) {
      $IgnoredLineCount++;
      next LINE;
    }

    @LineWords = quotewords($InDelim, 0, $Line);
    if (@LineWords != $NumOfCols) {
      $IgnoredLineCount++;
      warn "Warning: The number of data fields, @LineWords, in $TextDataFile must be $NumOfCols.\nIgnoring line number $LineCount: $Line...\n";
      next LINE;
    }
    $DataKey = $LineWords[$DataKeyColIndex];

    if ($UniqueDataKeyMode) {
      if (exists $TextDataMapRef->{DataKey}{$DataKey}) {
	$IgnoredLineCount++;
	warn "Warning: The data key, $DataKey, in data column key number, $DataKeyColNum, is already present.\nIgnoring line number $LineCount: $Line...\n";
	next LINE;
      }
      push @{$TextDataMapRef->{DataKeys}}, $DataKey;
      $TextDataMapRef->{DataKey}{$DataKey} = $DataKey;
    }
    else {
      if (!exists $TextDataMapRef->{DataKey}{$DataKey}) {
	push @{$TextDataMapRef->{DataKeys}}, $DataKey;
	$TextDataMapRef->{DataKey}{$DataKey} = $DataKey;

	for $ColIndex (0 .. $#DataColNums) {
	  $ColNum = $DataColNums[$ColIndex];
	  $ColID = $DataColIDs[$ColIndex];
	  @{$TextDataMapRef->{$ColID}{$DataKey}} = ();
	}
      }
    }

    # Track column data values...
    for $ColIndex (0 .. $#DataColNums) {
      $ColID = $DataColIDs[$ColIndex];

      $ColNum = $DataColNums[$ColIndex];
      $ColValue = $LineWords[$ColNum - 1];

      if ($UniqueDataKeyMode) {
	$TextDataMapRef->{$ColID}{$DataKey} = $ColValue;
      }
      else {
	push @{$TextDataMapRef->{$ColID}{$DataKey}}, $ColValue;
      }
    }

  }

  print "\nTotal number of lines in file $TextDataFile: $LineCount\n";
  print "Total number of lines ignored: $IgnoredLineCount\n";

  close TEXTDATAFILE;
}

# Returns a 32 bit integer hash code using One-at-a-time algorithm By Bob Jenkins [Ref 38]. It's also implemented in
# Perl for internal hash keys in hv.h include file.
#
# It's not clear how to force Perl perform unsigned integer arithmetic irrespective of the OS/Platform and
# the value of use64bitint flag used during its compilation.
#
# In order to generate a consistent 32 bit has code across OS/platforms, the following methodology appear
# to work:
#
#    o Use MaxHashCodeMask to retrieve appropriate bits after left shifting by bit operators and additions
#    o Stay away from "use integer" to avoid signed integer arithmetic for bit operators
#
#
#   MaxHashCodeMask (2147483647) corresponds to the maximum value which can be stored in 31 bits
#
my($MaxHashCodeMask);
$MaxHashCodeMask = 2**31 - 1;

sub HashCode {
  my($String) = @_;
  my($HashCode, $Value, $ShiftedHashCode);

  $HashCode = 0;
  for $Value (unpack('C*', $String)) {
    $HashCode += $Value;

    $ShiftedHashCode = $HashCode << 10;
    if ($ShiftedHashCode > $MaxHashCodeMask) {
      $ShiftedHashCode = $ShiftedHashCode & $MaxHashCodeMask;
    }

    $HashCode += $ShiftedHashCode;
    if ($HashCode > $MaxHashCodeMask) {
      $HashCode = $HashCode & $MaxHashCodeMask;
    }

    $HashCode ^= ($HashCode >> 6);
  }

  $ShiftedHashCode = $HashCode << 3;
  if ($ShiftedHashCode > $MaxHashCodeMask) {
    $ShiftedHashCode = $ShiftedHashCode & $MaxHashCodeMask;
  }

  $HashCode += $ShiftedHashCode;
  if ($HashCode > $MaxHashCodeMask) {
    $HashCode = $HashCode & $MaxHashCodeMask;
  }
  $HashCode ^= ($HashCode >> 11);

  $ShiftedHashCode = $HashCode << 15;
  if ($ShiftedHashCode > $MaxHashCodeMask) {
    $ShiftedHashCode = $ShiftedHashCode & $MaxHashCodeMask;
  }

  $HashCode += $ShiftedHashCode;
  if ($HashCode > $MaxHashCodeMask) {
    $HashCode = $HashCode & $MaxHashCodeMask;
  }
  return $HashCode;
}

# Check out the string: Is it defined and has a non zero length?
sub IsEmpty {
  my($TheString) = @_;
  my($Status) = 1;

  $Status = (defined($TheString) && length($TheString)) ? 0 : 1;

  return $Status;
}

# Is first specified number power of second specified number...
sub IsNumberPowerOfNumber {
  my($FirstNum, $SecondNum) = @_;
  my($PowerValue);

  $PowerValue = log($FirstNum)/log($SecondNum);

  return IsInteger($PowerValue) ? 1 : 0;
}

# Check out the string: Is it an integer?
sub IsInteger {
  my($TheString) = @_;
  my($Status) = 0;

  if (defined($TheString) && length($TheString)) {
    $TheString = RemoveLeadingAndTrailingWhiteSpaces($TheString);
    $TheString =~ s/^[+-]//;
    $Status = ($TheString =~ /[^0-9]/) ? 0 : 1;
  }
  return $Status;
}

# Check out the string: Is it an integer with value > 0?
sub IsPositiveInteger {
  my($TheString) = @_;
  my($Status) = 0;

  $Status = IsInteger($TheString) ? ($TheString > 0 ? 1 : 0) : 0;

  return $Status;
}


# Check out the string: Is it a float?
sub IsFloat {
  my($TheString) = @_;
  my($Status) = 0;

  if (defined($TheString) && length($TheString)) {
    $TheString = RemoveLeadingAndTrailingWhiteSpaces($TheString);
    $TheString =~ s/^[+-]//;
    $Status = ($TheString =~ /[^0-9.eE]/) ? 0 : (((length($TheString) == 1) && ($TheString =~ /[.eE]/)) ? 0 : 1);
  }
  return $Status;
}

# Check out the string: Is it defined and has a non zero length?
sub IsNotEmpty {
  my($TheString) = @_;
  my($Status);

  $Status = IsEmpty($TheString) ? 0 : 1;

  return $Status;
}

# Check out the string: Does it only contain numerical data?
sub IsNumerical {
  my($TheString) = @_;
  my($Status) = 0;

  if (defined($TheString) && length($TheString)) {
    $TheString = RemoveLeadingAndTrailingWhiteSpaces($TheString);
    $TheString =~ s/^[+-]//;
    $Status = ($TheString =~ /[^0-9.eE]/) ? 0 : (((length($TheString) == 1) && ($TheString =~ /[.eE]/)) ? 0 : 1);
  }
  return $Status;
}

# Join different words using delimiter and quote parameters. And return as
# a string value.
sub JoinWords {
  my($Words, $Delim, $Quote) = @_;

  if (!@$Words) {
    return "";
  }

  $Quote = $Quote ? "\"" : "";
  my(@NewWords) = map { (defined($_) && length($_)) ? "${Quote}$_${Quote}" : "${Quote}${Quote}" } @$Words;

  return join $Delim, @NewWords;
}

# Split string value containing quoted or unquoted words in to an array containing
# unquoted words.
#
# This function is used to split strings generated by JoinWords.
#
sub SplitWords {
  my($Line, $Delim) = @_;

  if (!$Line) {
    return ();
  }

  # Is it a quoted string?
  if ($Line =~ /^\"/) {
    # Take out first and last quote...
    $Line =~ s/^\"//; $Line =~ s/\"$//;

    $Delim = "\"$Delim\"";
  }
  return split /$Delim/, $Line;
}

# Based on quote parameter, figure out what to do
sub QuoteAWord {
  my($Word, $Quote) = @_;
  my($QuotedWord);

  $QuotedWord = "";
  if ($Word) {
    $QuotedWord = $Word;
  }
  if ($Quote) {
    $QuotedWord = "\"$QuotedWord\"";
  }
  return ($QuotedWord);
}

# Remove leading white space characters from the string...
sub RemoveLeadingWhiteSpaces {
  my($InString) = @_;
  my($OutString, $TrailingString, $LeadingWhiteSpace);

  $OutString = $InString;
  if (length($InString) && ContainsWhiteSpaces($InString)) {
    $OutString =~ s/^([ \t\r\n\f]*)(.*?)$/$2/;
  }
  return $OutString;
}

# Remove Trailing white space characters from the string...
sub RemoveTrailingWhiteSpaces {
  my($InString) = @_;
  my($OutString, $LeadingString, $TrailingWhiteSpace);

  $OutString = $InString;
  if (length($InString) && ContainsWhiteSpaces($InString)) {
    $OutString =~ s/^(.*?)([ \t\r\n\f]*)$/$1/;
  }
  return $OutString;
}

# Remove both leading and trailing white space characters from the string...
sub RemoveLeadingAndTrailingWhiteSpaces {
  my($InString) = @_;
  my($OutString);

  $OutString = $InString;
  if (length($InString) && ContainsWhiteSpaces($InString)) {
    $OutString =~ s/^([ \t\r\n\f]*)(.*?)([ \t\r\n\f]*)$/$2/;
  }
  return $OutString;
}

# Wrap text string...
sub WrapText {
  my($InString, $WrapLength, $WrapDelimiter);
  my($OutString);

  $WrapLength = 40;
  $WrapDelimiter = "\n";
  if (@_ == 3) {
    ($InString, $WrapLength, $WrapDelimiter) = @_;
  }
  elsif (@_ == 2) {
    ($InString, $WrapLength) = @_;
  }
  else {
    ($InString, $WrapLength) = @_;
  }
  $OutString = $InString;
  if ($InString && (length($InString) > $WrapLength)) {
    $OutString = "";
    my($Index, $Length, $FirstPiece, $StringPiece);
    $Index = 0; $Length = length($InString);
    $FirstPiece = 1;
    for ($Index = 0; $Index < $Length; $Index += $WrapLength) {
      if (($Index + $WrapLength) < $Length) {
	$StringPiece = substr($InString, $Index, $WrapLength);
      }
      else {
	# Last piece of the string...
	$StringPiece = substr($InString, $Index, $WrapLength);
      }
      if ($FirstPiece) {
	$FirstPiece = 0;
	$OutString = $StringPiece;
      }
      else {
	$OutString .= "${WrapDelimiter}${StringPiece}";
      }
    }
  }
  return $OutString;
}

1;

__END__

=head1 NAME

TextUtil

=head1 SYNOPSIS

use TextUtil;

use TextUtil qw(:all);

=head1 DESCRIPTION

B<TextUtil> module provides the following functions:

AddNumberSuffix, ContainsWhiteSpaces, GetTextFileDataByNonUniqueKey,
GetTextFileDataByUniqueKey, GetTextLine, HashCode, IsEmpty, IsFloat, IsInteger,
IsNotEmpty, IsNumberPowerOfNumber, IsNumerical, IsPositiveInteger, JoinWords,
QuoteAWord, RemoveLeadingAndTrailingWhiteSpaces, RemoveLeadingWhiteSpaces,
RemoveTrailingWhiteSpaces, SplitWords, WrapText

=head1 FUNCTIONS

=over 4

=item B<AddNumberSuffix>

    $NumberWithSuffix = AddNumberSuffix($IntegerValue);

Returns number with appropriate suffix: 0, 1st, 2nd, 3rd, 4th, and so on.

=item B<ContainsWhiteSpaces>

    $Status = ContainsWhiteSpaces($TheString);

Returns 1 or 0 based on whether the string contains any white spaces.

=item B<GetTextLine>

    $Line = GetTextLine(\*TEXTFILE);

Reads next line from an already opened text file, takes out any carriage return,
and returns it as a string. NULL is returned for EOF.

=item B<GetTextFileDataByNonUniqueKey>

    GetTextFileDataByNonUniqueKey($TextDataFile, $TextDataMapRef,
                                  $DataKeyColNum, $InDelim);

Load data from a text file into the specified hash reference using a specific
column for non-unique data key values.

The lines starting with # are treated as comments and ignored. First line
not starting with # must contain column labels and the number of columns in
all other data rows must match the number of column labels.

The first column is assumed to contain data key value by default; all other columns
contain data as indicated in their column labels.

In order to avoid dependence of data access on the specified column labels, the
column data is loaded into hash with Column<Num> hash keys, where column number
start from 1. The data key column is not available as Colnum<Num> hash key;

The format of the data structure loaded into a specified hash reference is:

    @{$TextDataMapRef->{DataKeys}} - Array of unique data keys
    @{$TextDataMapRef->{ColLabels}} - Array of column labels
    @{$TextDataMapRef->{DataColIDs}} - Array of data column IDs
    $TextDataMapRef->{NumOfCols} - Number of columns
    %{$TextDataMapRef->{DataKey}} - Hash keys pair: <DataKey, DataKey>
    @{$TextDataMapRef->{DataCol<Num>}} - Hash keys pair with data as an array:
                                         <DataCol<Num>, DataKey>

=item B<GetTextFileDataByUniqueKey>

    GetTextFileDataByUniqueKey($TextDataFile, $TextDataMapRef, $DataKeyColNum,
                                $InDelim);

Load data from a text file into the specified hash reference using a a specific
column for unique data key values.

The lines starting with # are treated as comments and ignored. First line
not starting with # must contain column labels and the number of columns in
all other data rows must match the number of column labels.

The first column is assumed to contain data key value by default; all other columns
contain data as indicated in their column labels.

In order to avoid dependence of data access on the specified column labels, the
column data is loaded into hash with Column<Num> hash keys, where column number
start from 1. The data key column is not available as Colnum<Num> hash key;

The format of the data structure loaded into a specified hash reference is:

    @{$TextDataMapRef->{DataKeys}} - Array of unique data keys
    @{$TextDataMapRef->{ColLabels}} - Array of column labels
    @{$TextDataMapRef->{DataColIDs}} - Array of data column IDs
    $TextDataMapRef->{NumOfCols} - Number of columns
    %{$TextDataMapRef->{DataKey}} - Hash keys pair: <DataKey, DataKey>
    %{$TextDataMapRef->{DataCol<Num>}} - Hash keys pair: <DataCol<Num>, DataKey>

=item B<HashCode>

    $HashCode = HashCode($TheString);

Returns a 32 bit integer hash code using One-at-a-time algorithm By Bob Jenkins [Ref 38].
It's also implemented in Perl for internal hash keys in hv.h include file.

=item B<IsEmpty>

    $Status = IsEmpty($TheString);

Returns 1 or 0 based on whether the string is empty.

=item B<IsInteger>

    $Status = IsInteger($TheString);

Returns 1 or 0 based on whether the string is a positive integer.

=item B<IsPositiveInteger>

    $Status = IsPositiveInteger($TheString);

Returns 1 or 0 based on whether the string is an integer.

=item B<IsFloat>

    $Status = IsFloat($TheString);

Returns 1 or 0 based on whether the string is a float.

=item B<IsNotEmpty>

    $Status = IsNotEmpty($TheString);

Returns 0 or 1 based on whether the string is empty.

=item B<IsNumerical>

    $Status = IsNumerical($TheString);

Returns 1 or 0 based on whether the string is a number.

=item B<IsNumberPowerOfNumber>

    $Status = IsNumberPowerOfNumber($FirstNum, $SecondNum);

Returns 1 or 0 based on whether the first number is a power of second number.

=item B<JoinWords>

    $JoinedWords = JoinWords($Words, $Delim, $Quote);

Joins different words using delimiter and quote parameters, and returns it
as a string.

=item B<QuoteAWord>

    $QuotedWord = QuoteAWord($Word, $Quote);

Returns a quoted string based on I<Quote> value.

=item B<RemoveLeadingWhiteSpaces>

    $OutString = RemoveLeadingWhiteSpaces($InString);

Returns a string without any leading and traling white spaces.

=item B<RemoveTrailingWhiteSpaces>

    $OutString = RemoveTrailingWhiteSpaces($InString);

Returns a string without any trailing white spaces.

=item B<RemoveLeadingAndTrailingWhiteSpaces>

    $OutString = RemoveLeadingAndTrailingWhiteSpaces($InString);

Returns a string without any leading and traling white spaces.

=item B<SplitWords>

    @Words = SplitWords($Line, $Delimiter);

Returns an array I<Words> ontaining unquoted words generated after spliting
string value I<Line> containing quoted or unquoted words.

This function is used to split strings generated by JoinWords as replacement
for Perl's core module funtion Text::ParseWords::quotewords() which dumps core
on very long strings.

=item B<WrapText>

    $OutString = WrapText($InString, [$WrapLength, $WrapDelimiter]);

Returns a wrapped string. By default, I<WrapLenght> is I<40> and I<WrapDelimiter>
is Unix new line character.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

FileUtil.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
