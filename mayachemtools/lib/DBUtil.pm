package DBUtil;
#
# File: DBUtil.pm
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
use DBI;
use TextUtil;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw(DBConnect DBDisconnect DBFetchSchemaTableNames DBSetupDescribeSQL DBSetupSelectSQL DBSQLToTextFile);
@EXPORT_OK = qw();
%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Connect to a specified database...
sub DBConnect {
  my($DBDriver, $DBName, $DBHost, $DBUser, $DBPassword) = @_;
  my($DBHandle, $DataSource);

  if ($DBDriver eq "Oracle") {
    $DataSource = qq(DBI:$DBDriver:$DBHost);
  }
  else {
    $DataSource = qq(DBI:$DBDriver:database=$DBName);
    if ($DBHost) {
      $DataSource .= qq(;host=$DBHost);
    }
  }

  # Don't raise the error; otherwise, DBI functions termiates on encountering an error.
  # All terminations decisions are made outside of DBI functions...
  $DBHandle = DBI->connect($DataSource, $DBUser, $DBPassword, { RaiseError => 0, AutoCommit => 0 }) or croak "Couldn't connect to database...";

  return $DBHandle;
}

# Disconnect from a database...
sub DBDisconnect {
  my($DBHandle) = @_;

  $DBHandle->disconnect or carp "Couldn't disconnect from a database...";
}

# Fetch all table name for a database schema...
sub DBFetchSchemaTableNames {
  my($DBDriver, $DBHandle, $SchemaName) = @_;
  my(@SchemaTableNames, $SQL, $SQLHandle);

  @SchemaTableNames = ();

  $SchemaName = (defined $SchemaName && length $SchemaName) ? $SchemaName : "";

  if ($DBDriver eq "mysql") {
    # Switch schemas...
    $SQL = qq(USE $SchemaName);
    $SQLHandle = $DBHandle->prepare($SQL) or return @SchemaTableNames;
    $SQLHandle->execute or return @SchemaTableNames;
    $SQLHandle->finish or return @SchemaTableNames;

    # Setup to fetch table names...
    $SQL = qq(SHOW TABLES);
  }
  elsif ($DBDriver eq "Oracle") {
    $SQL = qq(SELECT SEGMENT_NAME FROM DBA_SEGMENTS WHERE OWNER = '$SchemaName' AND SEGMENT_TYPE = 'TABLE' ORDER BY SEGMENT_NAME);
  }
  elsif ($DBDriver =~ /^(Pg|Postgres)$/i) {
    $SQL = qq(SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = '$SchemaName');
  }
  $SQLHandle = $DBHandle->prepare($SQL) or return @SchemaTableNames;
  $SQLHandle->execute or return @SchemaTableNames;

  my(@RowValues, $TableName);
  while (@RowValues = $SQLHandle->fetchrow_array) {
    $TableName = ($DBDriver =~ /^(mysql|Oracle)$/i) ? uc($RowValues[0]) : $RowValues[0];
    if (defined $TableName && length $TableName) {
      push @SchemaTableNames, $TableName;
    }
  }
  $SQLHandle->finish or return @SchemaTableNames;

  return @SchemaTableNames;
}

# Setup describe SQL statement...
sub DBSetupDescribeSQL {
  my($DBDriver, $TableName, $SchemaName);
  my($DescribeSQL);

  $DBDriver = ""; $TableName = ""; $SchemaName = "";
  if (@_ == 3) {
    ($DBDriver, $TableName, $SchemaName) = @_;
  }
  else {
    ($DBDriver, $TableName) = @_;
  }
  $TableName = (defined $TableName && length $TableName) ? $TableName : "";
  $SchemaName = (defined $SchemaName && length $SchemaName) ? $SchemaName : "";

  $DescribeSQL = ($SchemaName) ? ("DESCRIBE " . "$SchemaName" . ".$TableName") : "DESCRIBE $TableName";

  if ($DBDriver eq "Oracle") {
    $DescribeSQL = qq(SELECT COLUMN_NAME "Column_Name", DECODE(NULLABLE, 'N','Not Null','Y','Null') "Null", DATA_TYPE "Data_Type", DATA_LENGTH "Data_Length", DATA_PRECISION "Data_Precision" FROM DBA_TAB_COLUMNS WHERE TABLE_NAME = '$TableName');
    if ($SchemaName) {
      $DescribeSQL .= qq( AND OWNER = '$SchemaName');
    }
    $DescribeSQL .= qq( ORDER BY COLUMN_ID);
  }
  elsif ($DBDriver =~ /^(Pg|Postgres)$/i) {
    $DescribeSQL = qq(SELECT COLUMN_NAME "Column_Name", data_type "Data_Type" FROM information_schema.columns WHERE table_name ='$TableName');
    if ($SchemaName) {
      $DescribeSQL .= " and table_schema = '$SchemaName'";
    }
  }

  return $DescribeSQL;
}

# Setup describe SQL statement...
sub DBSetupSelectSQL {
  my($DBDriver, $TableName, $SchemaName);
  my($SelectSQL);

  $DBDriver = ""; $TableName = ""; $SchemaName = "";
  if (@_ == 3) {
    ($DBDriver, $TableName, $SchemaName) = @_;
  }
  else {
    ($DBDriver, $TableName) = @_;
  }
  $TableName = (defined $TableName && length $TableName) ? $TableName : "";
  $SchemaName = (defined $SchemaName && length $SchemaName) ? $SchemaName : "";

  $SelectSQL = ($SchemaName) ? ("SELECT * FROM " . "$SchemaName" . ".$TableName") : "SELECT * FROM $TableName";

  return $SelectSQL;
}

# Prepare and execute a SQL statement and write out results into
# a text file.
sub DBSQLToTextFile {
  my($DBHandle, $SQL, $TextFile, $OutDelim, $OutQuote, $ExportDataLabels, $ExportLOBs, $ReplaceNullStr);
  my($SQLHandle, $Status);

  $Status = 1;
  $ExportDataLabels = 1;
  $ExportLOBs = 0;
  $ReplaceNullStr = "";
  if (@_ == 8) {
    ($DBHandle, $SQL, $TextFile, $OutDelim, $OutQuote, $ExportDataLabels, $ExportLOBs, $ReplaceNullStr) = @_;
  }
  elsif (@_ == 7) {
    ($DBHandle, $SQL, $TextFile, $OutDelim, $OutQuote, $ExportDataLabels, $ExportLOBs) = @_;
  }
  elsif (@_ == 6) {
    ($DBHandle, $SQL, $TextFile, $OutDelim, $OutQuote, $ExportDataLabels) = @_;
  }
  else {
    ($DBHandle, $SQL, $TextFile, $OutDelim, $OutQuote) = @_;
  }

  # Execute SQL statement...
  $SQLHandle = $DBHandle->prepare($SQL) or return $Status;
  $SQLHandle->execute() or return $Status;

  my($FieldsNum, @FieldNames, @RowValues, @ColNumsToExport, @ColLabels, $ColNum, $ColLabelsLine, @Values, $Value, $ValuesLine);

  $Status = 0;
  # Figure out which column numbers need to be exported...
  $FieldsNum = $SQLHandle->{NUM_OF_FIELDS};
  @FieldNames = @{$SQLHandle->{NAME}};
  @ColNumsToExport = ();
  if ($ExportLOBs) {
    @ColNumsToExport = (0 .. $#FieldNames);
  }
  else {
    my(@FieldTypes, @FieldTypeNames, $Type, $TypeName);
    @FieldTypes = @{$SQLHandle->{TYPE}};
    @FieldTypeNames = map { scalar $DBHandle->type_info($_)->{TYPE_NAME} } @FieldTypes;
    for $ColNum (0 .. $#FieldNames) {
      if ($FieldTypeNames[$ColNum] !~ /lob|bytea/i ) {
	push @ColNumsToExport, $ColNum;
      }
    }
  }

  if ($ExportDataLabels) {
    # Print out column labels...
    @ColLabels = ();
    for $ColNum (@ColNumsToExport) {
      push @ColLabels, $FieldNames[$ColNum];
    }
    $ColLabelsLine = JoinWords(\@ColLabels, $OutDelim, $OutQuote);
    print $TextFile "$ColLabelsLine\n";
  }
  # Print out row values...
  while (@RowValues = $SQLHandle->fetchrow_array) {
    @Values = ();
    for $ColNum (@ColNumsToExport) {
      if (defined($RowValues[$ColNum]) && length($RowValues[$ColNum])) {
	$Value = $RowValues[$ColNum];
      }
      else {
	$Value = $ReplaceNullStr ? $ReplaceNullStr : "";
      }
      push @Values, $Value;
    }
    $ValuesLine = JoinWords(\@Values, $OutDelim, $OutQuote);
    print $TextFile "$ValuesLine\n";
  }
  $SQLHandle->finish or return $Status;
  $Status = 0;

  return $Status;
}

1;

__END__

=head1 NAME

DBUtil

=head1 SYNOPSIS

use DBUtil;

use DBUtil qw(:all);

=head1 DESCRIPTION

B<DBUtil> module provides the following functions:

DBConnect, DBDisconnect, DBFetchSchemaTableNames, DBSQLToTextFile,
DBSetupDescribeSQL, DBSetupSelectSQL

DBUtil package uses Perl DBI for interacting with MySQL Oracle, and PostgreSQL
databases.

=head1 FUNCTIONS

=over 4

=item B<DBConnect>

    $DBHandle = DBConnect($DBDriver, $DBName, $DBHost, $DBUser, $DBPassword);

Connects to a database using specified parameters and returns a B<DBHandle>.

=item B<DBDisconnect>

    DBDisconnect($DBHandle);

Disconnects from a database specified by I<DBHandle>.

=item B<DBFetchSchemaTableNames>

    @SchemaTableNames = DBFetchSchemaTableNames($DBDriver, $DBHandle,
                       $SchemaName);

Returns an array of all the table names in a database I<SchemaName>.

=item B<DBSetupDescribeSQL>

    $DescribeSQL = DBSetupDescribeSQL($DBDriver, $TableName, [$SchemaName]);

Sets up and returns a SQL statement to describe a table for MySQ, Oracle or PostgreSQL.

=item B<DBSetupSelectSQL>

    $SelectSQL = DBSetupSelectSQL($DBDriver, $TableName, $SchemaName);

Sets up and returns a SQL statement to retrieve all columns from a table for MySQL,
Oracle, or PostgreSQL.

=item B<DBSQLToTextFile>

    $Status = DBSQLToTextFile($DBHandle, $SQL, \*TEXTFILE, $OutDelim,
              $OutQuote, [$ExportDataLabels, $ExportLOBs,
              $ReplaceNullStr]);

Executes a I<SQL> statement and export all data into a text file.

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
