package PackageInfo;
#
# File: PackageInfo.pm
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
use Carp;
use Text::ParseWords;
use TextUtil;
use FileUtil;

use vars qw($AUTOLOAD @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw(GetPackageKeyValue SetPackageKeyValue IsPackageKeyNameAvailable);
@EXPORT_OK = qw();
%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

#
# Load package data...
#
my(%PackageDataMap);
_LoadPackageData();

# Return value of a specific key...
#
sub GetPackageKeyValue {
  my($KeyName) = @_;

  return exists $PackageDataMap{$KeyName} ? $PackageDataMap{$KeyName} : 'Not Available';
}

# Set value of a specific key...
#
sub SetPackageKeyValue {
  my($KeyName, $KeyValue) = @_;

  $PackageDataMap{$KeyName} = $KeyValue;
}

# Check availability of a package key name...
#
sub IsPackageKeyNameAvailable {
  my($KeyName) = @_;

  return exists $PackageDataMap{$KeyName} ? 1 : 0;
}

# Implements Set<KeyName> and Get<KeyName> functions...
#
sub AUTOLOAD {
  my($KeyValue) = @_;
  my($PackageName, $FunctionName, $KeyName);

  ($PackageName, $FunctionName) = $AUTOLOAD =~ /^(.*?)::(.*?)$/;

  if (!($FunctionName =~ /^Get/ || $FunctionName =~ /^Set/)) {
    croak "Error: Can't locate function \"$FunctionName\" via package \"$PackageName\": This function is not automatically implemented by AUTOLOAD: Only Get<KeyName> and Set<KeyName> functions are implemented via AUTOLOAD...";
  }

  ($KeyName) = $FunctionName =~ /^[SG]et(.*?)$/;

  if ($FunctionName =~ /^Set/ && !defined($KeyValue)) {
    carp "Warning:  ${PackageName}::${FunctionName}: Didn't set value for key $KeyName: Key value for must be specified...\n";
    return undef;
  }

  if ($FunctionName =~ /^Get/) {
    return GetPackageKeyValue($KeyName);
  }
  elsif ($FunctionName =~ /^Set/) {
    return SetPackageKeyValue($KeyName, $KeyValue);
  }

}

# Load PackageInfo.csv files from <MayaChemTools>/lib directory...
#
sub _LoadPackageData {
  my($PackageDataFile, $MayaChemToolsLibDir, $Key, $Value, $Line, $InDelim, @LineWords);

  $MayaChemToolsLibDir = GetMayaChemToolsLibDirName();
  $PackageDataFile =  "$MayaChemToolsLibDir" . "/data/PackageInfo.csv";

  if (! -e "$PackageDataFile") {
    croak "Error: MayaChemTools package file, $PackageDataFile, is missing: Possible installation problems...";
  }

  #
  # Format:
  #
  # "Key","Value"
  # "PackageName","MayaChemTools"
  # ... ...
  #

  %PackageDataMap = ();
  $InDelim = "\,";

  open PACKAGEDATAFILE, "$PackageDataFile" or croak "Couldn't open $PackageDataFile: $! ...";

  # Skip lines up to column labels...
  LINE: while ($Line = GetTextLine(\*PACKAGEDATAFILE)) {
    if ($Line !~ /^#/) {
      last LINE;
    }
  }

  # Process key/value pairs...
  LINE: while ($Line = GetTextLine(\*PACKAGEDATAFILE)) {
    if ($Line =~ /^#/) {
      next LINE;
    }
    @LineWords = ();
    @LineWords = quotewords($InDelim, 0, $Line);
    if (@LineWords != 2) {
      croak "Error: The number of data fields, @LineWords, in $PackageDataFile must be 2.\nLine: $Line...";
    }
    ($Key, $Value) = @LineWords;

    if (exists $PackageDataMap{$Key}) {
      carp "Warning: Multiple entries for key, $Key, in $PackageDataFile. Ignoring current line.\nLine: $Line...";
      next LINE;
    }

    $PackageDataMap{$Key} = $Value;
  }
  close PACKAGEDATAFILE;
}

1;

__END__

=head1 NAME

PackageInfo

=head1 SYNOPSIS

use PackageInfo;

use PackageInfo qw(:all);

=head1 DESCRIPTION

B<PackageInfo> module provides the following functions:

GetPackageKeyValue, IsPackageKeyNameAvailable, SetPackageKeyValue

The functions to set and get package keyvalues not explicitly defined in this module
are implemented using Perl's AUTOLOAD functionality. These methods are generated
on-the-fly for a specified key:

    Set<KeyName>(<KeyValue>);
    $KeyValue = Get<KeyName>();

B<PackageInfo> module provides functionality to retrieve information about MayaChemTools
package from PackagaInfo.csv which contains the following types key name and values:

    "KeyName","KeyValue"
    "PackageName","MayaChemTools"
    "ReleaseDate","Oct 21, 2010"
    "VersionNumber","7.4"
    "DevSoftwareEnvironment","Cygwin on Windows XP"
    ... ...
    ... ...

=head1 FUNCTIONS

=over 4

=item B<GetPackageKeyValue>

    $KeyValue = GetPackageKeyValue($KeyName);

Returns B<KeyValue> for a specified I<KeyName>.

=item B<IsPackageKeyNameAvailable>

    $Status = IsPackageKeyNameAvailable($KeyName);

Returns 1 or 0 based on whether I<KeyName> is available in package info file.

=item B<SetPackageKeyValue>

    SetPackageKeyValue($KeyName, $KeyValue);

Sets I<KeyValue> for a I<KeyName>. No data is written to package info file.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
