#!/usr/bin/perl -w
#
# File: SortSDFiles.pl
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
use SDFileUtil;
use TextUtil;

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

my(@SDFilesList);
@SDFilesList = ExpandFileNames(\@ARGV, "sdf sd");

print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

print "Checking input SD file(s)...\n";
my(%SDFilesInfo);
RetrieveSDFilesInfo();

# Generate output files...
my($FileIndex);
if (@SDFilesList > 1) {
  print "\nProcessing SD files...\n";
}
for $FileIndex (0 .. $#SDFilesList) {
  if ($SDFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $SDFilesList[$FileIndex]...\n";
    SortSDFile($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Sort it out...
sub SortSDFile {
  my($Index) = @_;
  my($SDFile, $NewSDFile, $KeyDataFieldName);

  $SDFile = $SDFilesList[$Index];
  $NewSDFile = $SDFilesInfo{OutFile}[$Index];
  $KeyDataFieldName = $SDFilesInfo{KeyDataFieldName}[$Index];

  print "Generating new SD file $NewSDFile...\n";
  open NEWSDFILE, ">$NewSDFile" or die "Error: Couldn't open $NewSDFile: $! \n";
  open SDFILE, "$SDFile" or die "Error: Can't open $SDFile: $! \n";

  # Go over all compound records and store 'em using key value as hash...
  my(%KeyToCompundRecordsMap, @InvalidCompoundRecords, $CmpdCount, $CmpdString, @CmpdLines, %DataFieldValues, $KeyDataFieldValue);
  %KeyToCompundRecordsMap = ();
  @InvalidCompoundRecords = ();
  $CmpdCount = 0;

  COMPOUND: while ($CmpdString = ReadCmpdString(\*SDFILE)) {
      $CmpdCount++;
      @CmpdLines = split "\n", $CmpdString;
      %DataFieldValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
      #Make sure data field value is okay...
      if (!(IsNotEmpty($DataFieldValues{$KeyDataFieldName}) && ($DataFieldValues{$KeyDataFieldName} !~ /\n/))) {
	push @InvalidCompoundRecords, $CmpdString;
	if ($OptionsInfo{DetailLevel} >= 3 ) {
	  print "Ignoring compound record $CmpdCount: Contains empty value for key data field $KeyDataFieldName :\n $CmpdString\n\n";
	}
	elsif ($OptionsInfo{DetailLevel} >= 2) {
	  print "Ignoring compound record $CmpdCount: Contains empty value for key data field $KeyDataFieldName...\n";
	}
	next COMPOUND;
      }
      $KeyDataFieldValue = $DataFieldValues{$KeyDataFieldName};
      if ($OptionsInfo{KeyData} =~ /^numeric$/i) {
	if (!IsFloat($KeyDataFieldValue)) {
	  push @InvalidCompoundRecords, $CmpdString;
	  if ($OptionsInfo{DetailLevel} >= 3 ) {
	    print "Ignoring compound record $CmpdCount: Contains non-numerical value for key data field $KeyDataFieldName :\n $CmpdString\n\n";
	  }
	  elsif ($OptionsInfo{DetailLevel} >= 2) {
	    print "Ignoring compound record $CmpdCount: Contains non-numerical value for key data field $KeyDataFieldName...\n";
	  }
	  next COMPOUND;
	}
      }
      if (exists($KeyToCompundRecordsMap{$KeyDataFieldValue})) {
	# Append to existing coompund data...
	$KeyToCompundRecordsMap{$KeyDataFieldValue} .= "\n" . $CmpdString;
      }
      else {
	$KeyToCompundRecordsMap{$KeyDataFieldValue} = $CmpdString;
      }
  }

  if ($OptionsInfo{Sort} =~ /^ascending$/i) {
    if ($OptionsInfo{KeyData} =~ /^alphanumeric$/i) {
      for $KeyDataFieldValue (sort { lc($a) cmp lc($b) } keys %KeyToCompundRecordsMap ) {
	print NEWSDFILE "$KeyToCompundRecordsMap{$KeyDataFieldValue}\n";
      }
    }
    else {
      for $KeyDataFieldValue (sort { $a <=> $b } keys %KeyToCompundRecordsMap ) {
	print NEWSDFILE "$KeyToCompundRecordsMap{$KeyDataFieldValue}\n";
      }
    }
  }
  else {
    if ($OptionsInfo{KeyData} =~ /^alphanumeric$/i) {
      for $KeyDataFieldValue (sort { lc($b) cmp lc($a) } keys %KeyToCompundRecordsMap ) {
	print NEWSDFILE "$KeyToCompundRecordsMap{$KeyDataFieldValue}\n";
      }
    }
    else {
      for $KeyDataFieldValue (sort { $b <=> $a } keys %KeyToCompundRecordsMap ) {
	print NEWSDFILE "$KeyToCompundRecordsMap{$KeyDataFieldValue}\n";
      }
    }
  }
  # Append the records containing data not appropriate for sorting...
  if (@InvalidCompoundRecords) {
    print "Placing ", scalar(@InvalidCompoundRecords)," compound record(s) with invalid data field key data the end...\n";
    for $CmpdString (@InvalidCompoundRecords) {
      print NEWSDFILE "$CmpdString\n";
    }
  }
  close NEWSDFILE;
  close SDFILE;
}

# Retrieve information about input SD files...
sub RetrieveSDFilesInfo {
  my($Index, $SDFile, $FileDir, $FileName, $FileExt, $OutFileRoot,  $OutFile, $DataFieldName);

  %SDFilesInfo = ();

  @{$SDFilesInfo{FileOkay}} = ();
  @{$SDFilesInfo{OutFile}} = ();
  @{$SDFilesInfo{KeyDataFieldName}} = ();

  FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFile = $SDFilesList[$Index];
    $SDFilesInfo{FileOkay}[$Index] = 0;
    $SDFilesInfo{OutFile}[$Index] = "";
    $SDFilesInfo{KeyDataFieldName}[$Index] = "";

    if (!(-e $SDFile)) {
      warn "Warning: Ignoring file $SDFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($SDFile, "sd sdf")) {
      warn "Warning: Ignoring file $SDFile: It's not a SD file\n";
      next FILELIST;
    }
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($SDFile);
    if ($Options{root} && (@SDFilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($Options{root});
      if ($RootFileName && $RootFileExt) {
	$FileName = $RootFileName;
      }
      else {
	$FileName = $Options{root};
      }
      $OutFileRoot = $FileName;
    }
    else {
      $OutFileRoot = $FileName . "SortedByDataField";
    }

    $OutFile = $OutFileRoot . ".$FileExt";
    if (lc($OutFile) eq lc($SDFile)) {
      warn "Warning: Ignoring file $SDFile:Output file name, $OutFile, is same as input SD file name, $SDFile\n";
      next FILELIST;
    }
    if (!$Options{overwrite}) {
      if (-e $OutFile) {
	warn "Warning: Ignoring file $SDFile: The file $OutFile already exists\n";
	next FILELIST;
      }
    }
    # Setup data field name...
    if ($OptionsInfo{SpecifiedDataFieldName}) {
      $DataFieldName = $OptionsInfo{SpecifiedDataFieldName};
    }
    else {
      my($CmpdString, @CmpdLines, @DataFieldNames);
      @DataFieldNames = ();
      if (!open(SDFILE, "$SDFile")) {
	warn "Warning: Ignoring file $SDFile: Couldn't open it: $! \n";
	next FILELIST;
      }
      $CmpdString = ReadCmpdString(\*SDFILE);
      close SDFILE;

      @CmpdLines = split "\n", $CmpdString;
      @DataFieldNames = GetCmpdDataHeaderLabels(\@CmpdLines);
      $DataFieldName = $DataFieldNames[0];
    }

    $SDFilesInfo{FileOkay}[$Index] = 1;
    $SDFilesInfo{OutFile}[$Index] = "$OutFile";
    $SDFilesInfo{KeyDataFieldName}[$Index] = $DataFieldName;
  }
}

# Process option values...
sub ProcessOptions {
  $OptionsInfo{DetailLevel} = $Options{detail};

  $OptionsInfo{Key} = defined $Options{key} ? $Options{key} : undef;
  $OptionsInfo{SpecifiedDataFieldName} = "";
  if (defined $Options{key}) {
    $OptionsInfo{SpecifiedDataFieldName} = $Options{key};
  }

  $OptionsInfo{KeyData} = $Options{keydata};
  $OptionsInfo{Sort} = $Options{sort};

  $OptionsInfo{Overwrite} = defined $Options{overwrite} ? $Options{overwrite} : undef;
  $OptionsInfo{Root} = defined $Options{root} ? $Options{root} : undef;
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{detail} = 1;
  $Options{sort} = "ascending";
  $Options{keydata} = "numeric";
  if (!GetOptions(\%Options, "detail|d=i", "help|h",  "key|k=s", "keydata=s", "overwrite|o", "root|r=s", "sort|s=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{keydata} !~ /^(numeric|alphanumeric)$/i) {
    die "Error: The value specified, $Options{keydata}, for option \"--keydata\" is not valid. Allowed values: numeric or alphanumeric\n";
  }
  if ($Options{sort} !~ /^(ascending|descending)$/i) {
    die "Error: The value specified, $Options{sort}, for option \"-s --sort\" is not valid. Allowed values: ascending or descending\n";
  }
  if (!IsPositiveInteger($Options{detail})) {
    die "Error: The value specified, $Options{detail}, for option \"-d --detail\" is not valid. Allowed values: > 0\n";
  }
}

__END__

=head1 NAME

SortSDFiles.pl - Sort SDFile(s) using values for a data field

=head1 SYNOPSIS

SortSDFiles.pl SDFile(s)...

SortSDFiles.pl [B<-d, --detail> infolevel] [B<-h, --help>] [B<-k, --key> I<SD data field name>]
[B<--keydata> numeric | alphanumeric] [B<-o, --overwrite>] [B<-r, --root> rootname]
[B<-s, --sort> ascending | descending] [B<-w, --workingdir> dirname] SDFile(s)...

=head1 DESCRIPTION

Sort I<SDFile(s)> using values for a specified data field name key. Only one SD
data field name key can be specified for sorting. In an event of conflict during sorting
process, two similar values for a SD data field name key are simply transferred to
output files in order of their presence in input files. Additionally, compound records
with no data field name, empty field values, or field values containing multiple lines
are simply placed at the end. The file names are separated by space.The valid file
extensions are I<.sdf> and I<.sd>. All other file names are ignored. All the SD files in a
current directory can be specified either by I<*.sdf> or the current directory name.

=head1 OPTIONS

=over 4

=item B<-d, --detail> I<infolevel>

Level of information to print about compound records being ignored. Default: I<1>. Possible
values: I<1, 2 or 3>.

=item B<-h, --help>

Print this help message.

=item B<-k, --key> I<SD data field name>

I<SDFile(s)> data field name used for sorting compound records. Default value: I<first
data field name>. Compound records with no I<sdfieldname>, empty field values, field
values containing multiple lines, or field values inappropriate for sorting are simply placed
at the end.

=item B<--keydata> I<numeric | alphanumeric>

Data type for I<sdfieldname> values. Possible values: I<numeric or alphanumeric>. Default
value: I<numeric>. For I<alphanumeric> data values, comparison is case insensitive.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

New SD file name is generated using the root: <Root>.<Ext>. Default new file
name: <InitialSDFileName>SortedByDataField.<Ext>. This option is ignored for multiple
input files.

=item B<-s, --sort> I<ascending | descending>

Sorting order for SD data field values. Possible values: I<ascending or descending>.
Default value: I<ascending>.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To perform numerical sort in ascending order using first data field values and
generate a new SD file NewSample1.sdf, type:

    % SortSDFiles.pl -o -r NewSample1 Sample1.sdf

To perform numerical sort in descending order using MolWeight data field and
generate a new SD text file NewSample1.sdf, type:

    % SortSDFiles.pl -k MolWeight --keydata numeric -s descending
      -r NewSample1 -o Sample1.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

JoinSDFiles.pl, MergeTextFilesWithSD.pl, SplitSDFiles.pl, SDFilesToHTML.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
