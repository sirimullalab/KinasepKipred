package TimeUtil;
#
# File: TimeUtil.pm
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
use Time::localtime ();

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw(CTimeStamp FPFileTimeStamp ISO8601Date ISO8601Time ISO8601TimeStamp PDBFileTimeStamp SDFileTimeStamp TimeStamp MonthNameToNumber MonthNumberToFullName MonthNumberToAbbreviatedName WeekDayNameToNumber WeekDayNumberToFullName WeekDayNumberToAbbreviatedName);
@EXPORT_OK = qw();
%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

#
# Initialize package data...
#
my(%MonthNameToNumber, %MonthNumberToFullNameName, %MonthNumberToAbbreviatedName, %WeekDayNameToNumber, %WeekDayNumberToFullName, %WeekDayNumberToAbbreviatedName);
_InitializeData();

# Return CTime as default time stamp for MayaChemTools...
#
sub TimeStamp {
  return CTimeStamp();
}

# Generate ctime time stamp...
#
# Format: WDay Mon MDay HH:MM:SS YYYY
#
sub CTimeStamp {
  my($CTimeStamp);

  # Take out an extra space inserted between month name and day by ctime...
  $CTimeStamp = Time::localtime::ctime();
  $CTimeStamp =~ s/[ ]+/ /g;

  return $CTimeStamp;
}

# Generate ISO 8601 timestamp in extended format...
#
# Format: [YYYY]-[MM]-[DD]T[hh]:[mm]:[ss]
#
sub ISO8601TimeStamp {
  my($TimeStamp, $Sec, $Min, $Hour, $MDay, $Mon, $Year, $WDay, $YDay, $IsDst);

  ($Sec, $Min, $Hour, $MDay, $Mon, $Year, $WDay, $YDay, $IsDst) = _LocalTime();

  $TimeStamp = sprintf "%4i-%02i-%02iT%02i:%02i:%02i", $Year, $Mon, $MDay, $Hour, $Min, $Sec;

  return $TimeStamp;
}

# Generate ISO 8601 date...
#
# Format: [YYYY]-[MM]-[DD]
#
sub ISO8601Date {
  my($Date, $Sec, $Min, $Hour, $MDay, $Mon, $Year, $WDay, $YDay, $IsDst);

  ($Sec, $Min, $Hour, $MDay, $Mon, $Year, $WDay, $YDay, $IsDst) = _LocalTime();

  $Date = sprintf "%4i-%02i-%02i", $Year, $Mon, $MDay;

  return $Date;
}

# Generate ISO 8601 time in extended format...
#
# Format: [hh]:[mm]:[ss]
#
sub ISO8601Time {
  my($Time, $Sec, $Min, $Hour);

  ($Sec, $Min, $Hour) = _LocalTime();

  $Time = sprintf "%02i:%02i:%02i", $Hour, $Min, $Sec;

  return $Time;
}

# Generate MayaChemTools' FP file timestamp...
#
sub FPFileTimeStamp {
  return CTimeStamp();
}

# Generate PDB file timestamp...
#
sub PDBFileTimeStamp {
  my($TimeStamp, $Sec, $Min, $Hour, $MDay, $Mon, $Year, $WDay, $YDay, $IsDst, $MonthName);

  ($Sec, $Min, $Hour, $MDay, $Mon, $Year, $WDay, $YDay, $IsDst) = _LocalTime();

  $MonthName = uc MonthNumberToAbbreviatedName($Mon);
  $Year = substr($Year, -2, 2);

  $TimeStamp = sprintf "%02i-%3s-%2i", $MDay, $MonthName, $Year;

  return $TimeStamp;
}

# Generate SD file timestamp...
#
sub SDFileTimeStamp {
  my($TimeStamp, $Sec, $Min, $Hour, $MDay, $Mon, $Year, $WDay, $YDay, $IsDst);

  ($Sec, $Min, $Hour, $MDay, $Mon, $Year, $WDay, $YDay, $IsDst) = _LocalTime();

  $Year = substr($Year, -2, 2);

  $TimeStamp = sprintf "%02i%02i%02i%02i%02i", $Mon, $MDay, $Year, $Hour, $Min;

  return $TimeStamp;
}

# Get local time with modifications to data returned by native localtime function...
#
sub _LocalTime {
  my($Sec, $Min, $Hour, $MDay, $Mon, $Year, $WDay, $YDay, $IsDst);

  ($Sec, $Min, $Hour, $MDay, $Mon, $Year, $WDay, $YDay, $IsDst) = localtime;

  $Mon += 1;
  $Year += 1900;

  return ($Sec, $Min, $Hour, $MDay, $Mon, $Year, $WDay, $YDay, $IsDst);
}

# Return month number from full or three letter abbreviated month name...
sub MonthNameToNumber {
  my($Name) = @_;

  return (exists $MonthNameToNumber{lc $Name}) ? $MonthNameToNumber{lc $Name} : '';
}

# Return full month name from month number...
sub MonthNumberToFullName {
  my($Number) = @_;

  return (exists $MonthNumberToFullNameName{$Number}) ? $MonthNumberToFullNameName{$Number} : '';
}

# Return three letter abbreviated month name from month number...
sub MonthNumberToAbbreviatedName {
  my($Number) = @_;

  return (exists $MonthNumberToAbbreviatedName{$Number}) ? $MonthNumberToAbbreviatedName{$Number} : '';
}

# Return week daty number from full or three letter abbreviated week day name...
sub WeekDayNameToNumber {
  my($Name) = @_;

  return (exists $WeekDayNameToNumber{lc $Name}) ? $WeekDayNameToNumber{lc $Name} : '';
}

# Return full week day name from week day number...
sub WeekDayNumberToFullName {
  my($Number) = @_;

  return (exists $WeekDayNumberToFullName{$Number}) ? $WeekDayNumberToFullName{$Number} : '';
}

# Return three letter abbreviated week day name from week day number...
sub WeekDayNumberToAbbreviatedName {
  my($Number) = @_;

  return (exists $WeekDayNumberToAbbreviatedName{$Number}) ? $WeekDayNumberToAbbreviatedName{$Number} : '';
}

# Initialize week/month day/name data...
#
sub _InitializeData {

  %MonthNameToNumber = ('january' => 1, 'february' => 2, 'march' => 3, 'april' => 4,
			'may' => 5, 'june' => 6, 'july' => 7, 'august' => 8,
			'september' => 9, 'october' => 10, 'november' => 11, 'december' => 12,
			'jan' => 1, 'feb' => 2, 'mar' => 3, 'apr' => 4,
			'may' => 5, 'jun' => 6, 'jul' => 7, 'aug' => 8,
			'sep' => 9, 'oct' => 10, 'nov' => 11, 'dec' => 12);

  %MonthNumberToFullNameName = (1 => 'January', 2 => 'February', 3 => 'March', 4 => 'April',
				5 => 'May', 6 => 'June', 7 => 'July', 8 => 'August',
				9 => 'September', 10 => 'October', 11 => 'November', 12 => 'December');

  %MonthNumberToAbbreviatedName = (1 => 'Jan', 2 => 'Feb', 3 => 'Mar', 4 => 'Apr',
				   5 => 'May', 6 => 'Jun', 7 => 'Jul', 8 => 'Aug',
				   9 => 'Sep', 10 => 'Oct', 11 => 'Nov', 12 => 'Dec');

  %WeekDayNameToNumber = ('sunday' => 1, 'monday' => 2, 'tuesday' => 3, 'wednesday' => 4,
			  'thursday' => 5, 'friday' => 6, 'saturday' => 7,
			  'sun' => 1, 'mon' => 2, 'tue' => 3, 'wed' => 4,
			  'thu' => 5, 'fri' => 6, 'sat' => 7);

  %WeekDayNumberToFullName = (1 => 'Sunday', 2 => 'Monday', 3 => 'Tuesday',
			      4 => 'Wednesday', 5 => 'Thursday', 6 => 'Friday', 7 => 'Saturday');

  %WeekDayNumberToAbbreviatedName = (1 => 'Sun', 2 => 'Mon', 3 => 'Tue',
				     4 => 'Wed', 5 => 'Thu', 6 => 'Fri', 7 => 'Sat');
}

1;

__END__

=head1 NAME

TimeUtil

=head1 SYNOPSIS

use TimeUtil;

use TimeUtil qw(:all);

=head1 DESCRIPTION

B<TimeUtil> module provides the following functions:

CTimeStamp, FPFileTimeStamp, ISO8601Date, ISO8601Time, ISO8601TimeStamp,
MonthNameToNumber, MonthNumberToAbbreviatedName, MonthNumberToFullName,
PDBFileTimeStamp, SDFileTimeStamp, TimeStamp, WeekDayNameToNumber,
WeekDayNumberToAbbreviatedName, WeekDayNumberToFullName

=head1 FUNCTIONS

=over 4

=item B<CTimeStamp>

    $CTimeStamp = CTimeStamp();

Returns B<CTimeStamp> string using the following format: WDay Mon MDay HH:MM:SS YYYY

=item B<FPFileTimeStamp>

    $FPFileTimeStamp = FPFileTimeStamp();

Returns fingerints B<FP> file time stamp string for MayaChemTools package. It corresponds to
B<CTimeStamp>.

=item B<ISO8601Date>

    $Date = ISO8601Date();

Returns ISO8601 B<Date> string using the following format: [YYYY]-[MM]-[DD]

=item B<ISO8601Time>

    $Time = ISO8601Time();

Returns ISO8601 B<Time> string using the following extended format: [hh]:[mm]:[ss]

=item B<ISO8601TimeStamp>

    $TimeStamp = ISO8601TimeStamp();

Returns ISO8601 B<TimeStamp> string using the following extended format: [YYYY]-[MM]-[DD]T[hh]:[mm]:[ss]

=item B<MonthNameToNumber>

    $Number = MonthNameToNumber($Name);

Return month B<Number> for full month I<Name> or three letter abbreviated month I<Name>.

=item B<MonthNumberToAbbreviatedName>

    $AbbrevMonthName = MonthNumberToAbbreviatedName($Number);

Returns three letter B<AbbrevMonthName> for month I<Number>.

=item B<MonthNumberToFullName>

    $Name = MonthNumberToFullName($Number);

Returns full month B<Name> for month I<Number>.

=item B<PDBFileTimeStamp>

    $TimeStamp = PDBFileTimeStamp();

Returns PDB file B<TimeStamp> using the following format: DD-MMM-YY

=item B<SDFileTimeStamp>

    $TimeStamp = SDFileTimeStamp();

Returns SD file B<TimeStamp> using the following format: MMDDYYHHMM

=item B<TimeStamp>

    $TimeStamp = TimeStamp();

Returns deafult I<TimeStamp> for MayaChemTools. It corresponds to B<CTimeStamp>.

=item B<WeekDayNameToNumber>

    $Number = WeekDayNameToNumber($Name);

Returns week day B<Number> from full week day I<Name> or three letter abbreviated week
day I<Name>.

=item B<WeekDayNumberToAbbreviatedName>

    $Name = WeekDayNumberToAbbreviatedName($Number);

Returns three letter abbreviates week day B<Name> for week day I<Number>.

=item B<WeekDayNumberToFullName>

    $Name = WeekDayNumberToFullName($Number);

Returns full week day B<Name> for week day I<Number>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

FileUtil.pm, TextUtil.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
