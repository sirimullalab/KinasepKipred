#!/usr/bin/perl -w
#
# File: DBTablesToTextFiles.pl
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
use FindBin; use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use File::Basename;
use Text::ParseWords;
use Benchmark;
use FileUtil;
use DBUtil;

my($ScriptName, %Options, $StartTime, $EndTime, $TotalTime);

# Autoflush STDOUT
$| = 1;

# Starting message...
$ScriptName = basename($0);
print "\n$ScriptName: Starting...\n\n";
$StartTime = new Benchmark;

# Get the options and setup script...
SetupScriptUsage();
if ($Options{help} || @ARGV < 1) {
  die GetUsageFromPod("$FindBin::Bin/$ScriptName");
}

my($DBDriver, $DBHost, $DBName, $DBUser, $DBPassword, $DBMode, $ExportDataLabels, $ExportLOBs, $OutDelim, $OutQuote, $ReplaceNullStr);
ProcessOptions();

# Collect input parameters information...
print "Checking input parameter(s)...\n";
my(@DBSQLStatements, @DBTextFiles);
RetrieveDBInfo();

# Connect to database...
my($DBHandle);
print "Connecting to $DBDriver:database=$DBName as $DBUser...\n";
$DBHandle = DBConnect($DBDriver, $DBName, $DBHost, $DBUser, $DBPassword);

# Generate text files...
if (@DBTextFiles > 1) {
  print "Generating text files...\n";
}
my($Index, $TextFile, $SQL);
TEXTFILE: for $Index (0 .. $#DBTextFiles) {
  $TextFile = $DBTextFiles[$Index];
  $SQL = $DBSQLStatements[$Index];

  if (@DBTextFiles > 1) {
    print "\nGenerating text file $TextFile...\n";
  }
  else {
    print "Generating text file $TextFile...\n";
  }
  print "Processing SQL statement \"$SQL\"...\n";

  if (!open TEXTFILE, ">$TextFile") {
    warn "Warning: Abandoning $TextFile generation: Couldn't open it: $! \n";
    next TEXTFILE;
  }

  if (DBSQLToTextFile($DBHandle, $SQL, \*TEXTFILE, $OutDelim, $OutQuote, $ExportDataLabels, $ExportLOBs, $ReplaceNullStr)) {
    warn "Warning: Abandoning $TextFile generation...\n";
    next TEXTFILE;
  }
  close TEXTFILE;
}
print "\nDisconnecting from  $DBDriver:database=$DBName...\n";
DBDisconnect($DBHandle);

print "$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Collect input parameters information...
sub RetrieveDBInfo {
  my($FileExt, $UserFileName, $FileDBPrefix);

  # Setup out file ext...
  $FileExt = ($Options{outdelim} =~ /^tab$/i) ? "tsv" : "csv";

  # Get user specified information...
  $UserFileName = "";
  if ($Options{root} && (@ARGV == 1)) {
    my($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($Options{root});
    if ($RootFileName && $RootFileExt) {
      $UserFileName = $RootFileName;
    }
    else {
      $UserFileName = $Options{root};
    }
  }
  $FileDBPrefix = ($DBMode =~ /^exportdata$/i) ? "Export" : "Describe";

  my($TableName, $SQL, $FileName);
  # Go over all the input parameters...
  @DBSQLStatements = ();
  @DBTextFiles = ();
  for $TableName (@ARGV) {
    $TableName = $TableName;
    $SQL = ($DBMode =~ /^exportdata$/i) ? DBSetupSelectSQL($DBDriver, $TableName) : DBSetupDescribeSQL($DBDriver, $TableName);
    push @DBSQLStatements, $SQL;
    $FileName = $UserFileName ? $UserFileName : ("$FileDBPrefix" . "$TableName");
    $FileName .= ".$FileExt";
    if (!$Options{overwrite}) {
      if (-e $FileName) {
	die "Error: The file $FileName already exists.\n";
      }
    }
    push @DBTextFiles, $FileName;
  }
}

# Process option values...
sub ProcessOptions {

  $DBDriver = $Options{dbdriver} ? $Options{dbdriver} : (exists $ENV{DBI_DRIVER} ? $ENV{DBI_DRIVER} : "") ;
  if ($DBDriver) {
    if ($DBDriver =~ /^oracle$/i) {
      $DBDriver = "Oracle";
    }
    elsif ($DBDriver =~ /^mysql$/i) {
      $DBDriver = "mysql";
    }
    elsif ($DBDriver =~ /^(Pg|Postgres)$/i) {
      $DBDriver = "Pg";
    }
    else {
      if ($Options{dbdriver}) {
	die "Error: The value specified, $DBDriver, for option \"-d --dbdriver\" is not valid. Allowed values: MySQL, Oracle, Postgres or Pg\n";
      }
      else {
	die "Error: The value specified, $DBDriver, using environment variable DBI_DRIVER not valid. Allowed values: MySQL, Oracle, Postgres or Pg\n";
      }
    }
  }
  else {
    $DBDriver = "mysql";
  }
  $DBHost = $Options{dbhost} ? $Options{dbhost} : (exists $ENV{DBI_HOST} ? $ENV{DBI_HOST} : "127.0.0.1");
  $DBName = $Options{dbname} ? $Options{dbname} : (exists $ENV{DBI_NAME} ? $ENV{DBI_NAME} : "");
  if (!$DBName) {
    if ($DBDriver =~ /^mysql$/i) {
      $DBName = "mysql";
    }
    elsif ($DBDriver =~ /^pg|Postgres$/i) {
      $DBName = "postgres";
    }
  }
  $DBUser = $Options{dbusername} ? $Options{dbusername} : (exists $ENV{DBI_USER} ? $ENV{DBI_USER} : "") ;
  if (!$DBUser) {
    die "Error: No database username specified. Use \"--dbusername\" option or environment variable DBI_USER to enter a valid value.\n";
  }
  $DBPassword = $Options{dbpassword} ? $Options{dbpassword} : (exists $ENV{DBI_PASS} ? $ENV{DBI_PASS} : "") ;
  if (!$DBPassword) {
    die "Error: No database password specified. Use \"--dbpassword\" option or environment variable DBI_PASS to enter a valid value.\n";
  }
  $DBMode = $Options{mode};
  $ExportLOBs = ($Options{exportlobs} =~ /^yes$/i) ? 1 : 0;
  $ExportDataLabels = ($DBMode =~ /^describetable$/i) ? 1 : (($Options{exportdatalabels} =~ /^yes$/i) ? 1 : 0);

  $OutDelim = ($Options{outdelim} =~ /^tab$/i ) ? "\t" : (($Options{outdelim} =~ /^semicolon$/i) ? "\;" : "\,");
  $OutQuote = ($Options{quote} =~ /^yes$/i) ? 1 : 0;

  $ReplaceNullStr = (defined($Options{replacenullstr}) && length($Options{replacenullstr})) ? $Options{replacenullstr} : "";
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{mode} = "exportdata";
  $Options{exportlobs} = "no";
  $Options{exportdatalabels} = "yes";
  $Options{outdelim} = "comma";
  $Options{quote} = "yes";
  if (!GetOptions(\%Options, "dbdriver|d=s", "dbhost=s", "dbname=s", "dbpassword=s", "dbusername=s", "exportdatalabels=s", "exportlobs=s", "help|h",  "mode|m=s", "outdelim=s", "overwrite|o", "quote|q=s", "root|r=s", "replacenullstr=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{exportdatalabels} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{exportlobs}, for option \"--exportdatalabels\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{exportlobs} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{exportlobs}, for option \"--exportlobs\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{mode} !~ /^(exportdata|describetable)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: exportdata, or describetable\n";
  }
  if ($Options{outdelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{outdelim}, for option \"--outdelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if ($Options{quote} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not valid. Allowed values: yes or no\n";
  }
}

__END__

=head1 NAME

DBToTextFiles.pl - Export data from database TableName(s) into CSV/TSV text files

=head1 SYNOPSIS

DBTablesToTextFiles.pl TableName(s)...

DBTablesToTextFiles.pl [B<-d, --dbdriver> mysql | Oracle | Postgres or Pg] [B<--dbhost > hostname]
[B<--dbname> databasename] [B<--dbpassword> password] [B<--dbusername> username]
[B<--exportdatalabels> yes | no] [B<--exportlobs> yes | no] [B<-h, --help>]
[B<-m, --mode> exportdata | describetable] [B<-o, --overwrite>] [B<--outdelim> comma | tab | semicolon]
[B<-q, --quote> yes | no] [B<-r, --root> rootname] [B<--replacenullstr string>]
[B<-w --workingdir> dirname] TableName(s)...

=head1 DESCRIPTION

Export data from MySQL, Oracle or PostgreSQL database tables into CSV/TSV text files. Or perform
describe on all tables and store its output into CSV/TSV text files. A text file is generated
for each table.

=head1 OPTIONS

=over 4

=item B<-d, --dbdriver> I<mysql | Oracle | Postgres or Pg>

Database driver name. Possible values: I<mysql, Oracle, Postgres or Pg>. Default: I<MySQL> or value of
environment variable DBI_DRIVER. This script has only been tested with MySQL, Oracle
and PostgreSQL drivers.

=item B<--dbhost > I<hostname>

Database host name. Default: I<127.0.0.1> for both MySQL and Oracle. For remote
databases, specify complete remote host domain: I<dbhostname.org> or something
like it.

=item B<--dbname> I<databasename>

Database name. Default: mysql for MySQL, postgres for PostgreSQL and none for Oracle.
For connecting to local/remote Oracle databases, this value can be left undefined assuming
B<--dbhost> is correctly specified.

=item B<--dbpassword> I<password>

Database user password. Default: I<none> and value of environment variable DBI_PASS
is used for connecting to database.

=item B<--dbusername> I<username>

Database user name. Default: I<none> and value of environment variable DBI_USER is
used for connecting to database.

=item B<--exportdatalabels> I<yes | no>

This option is mode specific and controls exporting of column data labels during
exportdata mode. Possible values: I<yes or no>. Default: I<yes>.

=item B<--exportlobs> I<yes | no>

This option is mode specific and controls exporting of CLOB/BLOB or BYTEA data columns during
exportdata mode. Possible values: I<yes or no>. Default: I<no>.

=item B<-h, --help>

Print this help message.

=item B<-m, --mode> I<exportdata | describetable>

Data selection criterion from database. Possible values: I<exportdata or describetable>.
Default value: I<exportdata>.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<--outdelim> I<comma | tab | semicolon>

Output text file delimiter. Possible values: I<comma, tab, or semicolon>
Default value: I<comma>.

=item B<-q, --quote> I<yes | no>

Put quotes around column values in output text file. Possible values: I<yes or
no>. Default value: I<yes>.

=item B<-r, --root> I<rootname>

New file name is generated using the root:<Root>.<Ext>. Default new file
file names: <Mode><TableName>.<Ext>. The csv and tsv <Ext> values are used
for comma/semicolon, and tab delimited text files respectively.This option is
ignored for multiple input table names.

=item B<--replacenullstr> I<string>

Replace NULL or undefined row values with specified value. Default: I<none>

For importing output text files into MySQL database using "load data local infile '<tablename>.tsv'
into table <tablename>" command, use I<--raplacenullstr "NULL"> in conjunction with I<--exportdatalabels no>,
I<--quote no>, and I<--outdelim tab> options: it'll generate files for direct import into MySQL assuming
tables already exists.

=item B<-w --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To export all data in user and user_info tables from a MySQL server running on a local machine
using username/password from DBI_USER and DBI_PASS environmental variables, type:

    % DBTablesToTextFiles.pl -o user user_info

To describe user and user_info tables in a MySQL server running on a remote machine using explicit
username/password and capturing the output into a DescribeTables.csv file, type:

    % DBTablesToTextFiles.pl --dbdriver mysql --dbuser <name> --dbpassword
      <pasword> --dbname mysql --dbhost <mysqlhostname.org>
      -r DescribeTable -m describetable -o user user_info

To describe table all_tables in Oracle running on a remote machine using explicit
username/password and capturing the output into a DescribeAllTable.tsv file, type:

    % DBTablesToTextFiles.pl --dbdriver Oracle --dbuser <name> --dbpassword
      <pasword> --dbhost <oraclehostname.com> -r DescribeAllTable
      -m describetable --outdelim tab --quote no -o all_tables

To export all data in user and user_info tables from MySQL server running on a local machine
using explicit username/password and capturing the data in ExportTables.tsv file with empty
values substituted with NULL and no column labels, type:

    % DBTablesToTextFiles.pl --dbdriver Oracle --dbuser <name> --dbpassword
      <pasword> -r ExportTables --outdelim tab --quote no --replacenullstr
      "\N" -m exportdata --exportlobs no  --exportdatalabels no -o
      user user_info

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

DBSchemaTablesToTextFiles.pl, DBSQLToTextFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
