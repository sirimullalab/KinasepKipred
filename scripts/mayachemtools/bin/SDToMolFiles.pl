#!/usr/bin/perl -w
#
# File: SDToMolFiles.pl
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
use SDFileUtil;
use FileUtil;

my($ScriptName, %Options, $StartTime, $EndTime, $TotalTime);

# Autoflush STDOUT
$| = 1;

# Starting message...
$ScriptName = basename $0;
print "\n$ScriptName:Starting...\n\n";
$StartTime = new Benchmark;

# Get the options and setup script...
SetupScriptUsage();
if ($Options{help} || @ARGV < 1) {
  die GetUsageFromPod("$FindBin::Bin/$ScriptName");
}

my(@SDFilesList);
@SDFilesList = ExpandFileNames(\@ARGV, "sdf sd");

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Setup information about input files...
print "Checking input SD file(s)...\n";
my(%SDFilesInfo);
RetrieveSDFilesInfo();

# Process input files..
my($FileIndex);
if (@SDFilesList > 1) {
  print "\nProcessing SD files...\n";
}
for $FileIndex (0 .. $#SDFilesList) {
  if ($SDFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $SDFilesList[$FileIndex]...\n";
    GenerateMolFiles($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Generate MOL files for a SD file...
#
sub GenerateMolFiles {
  my($FileIndex) = @_;
  my($SDFile, $MOLFile, $MOLFileRoot, $OutFileRoot, $OverwriteFiles, $UseDataField, $DataFieldName, $UseMolName, $CmpdCount, $MolEndDelimiter, $CmpdString, @CmpdLines, %DataFieldValues);

  $SDFile = $SDFilesList[$FileIndex];

  if (!open SDFILE, "$SDFile") {
    warn "Warning: Ignoring file $SDFile: Couldn't open it: $! \n";
    return;
  }

  $CmpdCount = 0;
  $MolEndDelimiter = "M  END";

  $OutFileRoot = $SDFilesInfo{OutFileRoot}[$FileIndex];
  $OverwriteFiles = $OptionsInfo{OverwriteFiles};

  $UseDataField = ($OptionsInfo{Mode} =~ /^DataField$/i) ? 1 : 0;
  $DataFieldName = $OptionsInfo{DataField};

  $UseMolName = ($OptionsInfo{Mode} =~ /^MolName$/i) ? 1 : 0;

  CMPDSTRING: while ($CmpdString = ReadCmpdString(\*SDFILE)) {
    $CmpdCount++;

    # Setup MOL file name...
    $MOLFileRoot = '';
    if ($UseDataField) {
      @CmpdLines = split "\n", $CmpdString;
      %DataFieldValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
      if (exists $DataFieldValues{$DataFieldName}) {
	$MOLFileRoot = $DataFieldValues{$DataFieldName};
      }
    }
    elsif ($UseMolName) {
      @CmpdLines = split "\n", $CmpdString;
      $MOLFileRoot = $CmpdLines[0];
    }

    # Check for any invalid file name characters in data field or molname values...
    if ($MOLFileRoot && $MOLFileRoot =~ /[^a-zA-Z0-9_]/) {
      $MOLFileRoot =~ s/[^a-zA-Z0-9_]//g;
    }
    # Fall back plan for MOL file name...
    if (!$MOLFileRoot) {
      $MOLFileRoot = "${OutFileRoot}Cmpd${CmpdCount}";
    }

    $MOLFile = "${MOLFileRoot}.mol";

    if (!$OverwriteFiles) {
      if (-e $MOLFile) {
	warn "Warning: Ignoring compound number, $CmpdCount, in $SDFile: New MOL file, $MOLFile, already exists\n";
	next CMPDSTRING;
      }
    }

    if (!($CmpdString =~ /$MolEndDelimiter/)) {
      warn "Warning: Ignoring compound number, $CmpdCount, in $SDFile: Invalid compound data\n";
      next CMPDSTRING;
    }

    # Write out MOL file...

    print "Generating $MOLFile file...\n";
    open MOLFILE, ">$MOLFile" or die "Error: Can't open $MOLFile: $! \n";
    ($CmpdString) = split "$MolEndDelimiter", $CmpdString;
    print MOLFILE "$CmpdString";
    print MOLFILE "$MolEndDelimiter\n";
    close MOLFILE;

  }

  close SDFILE;
}

# Retrieve information about SD files...
#
sub RetrieveSDFilesInfo {
  my($SDFile, $Index, $FileDir, $FileName, $FileExt, $OutFileRoot);

  %SDFilesInfo = ();
  @{$SDFilesInfo{FileOkay}} = ();
  @{$SDFilesInfo{OutFileRoot}} = ();

  FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFile = $SDFilesList[$Index];

    $SDFilesInfo{FileOkay}[$Index] = 0;
    $SDFilesInfo{OutFileRoot}[$Index] = '';

    $SDFile = $SDFilesList[$Index];
    if (!(-e $SDFile)) {
      warn "Warning: Ignoring file $SDFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($SDFile, "sd sdf")) {
      warn "Warning: Ignoring file $SDFile: It's not a SD file\n";
      next FILELIST;
    }

    # Setup output file root...
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($SDFile);

    if ($OptionsInfo{OutFileRoot} && (@SDFilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($OptionsInfo{OutFileRoot});
      if ($RootFileName && $RootFileExt) {
	$FileName = $RootFileName;
      }
      else {
	$FileName = $OptionsInfo{OutFileRoot};
      }
      $OutFileRoot = $FileName;
    }
    else {
      $OutFileRoot = "$FileName";
    }

    $SDFilesInfo{FileOkay}[$Index] = 1;
    $SDFilesInfo{OutFileRoot}[$Index] = $OutFileRoot;
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{DataField} = '';
  if ($Options{mode} =~ /^DataField$/i) {
    if (!$Options{datafield}) {
      die "Error: You must specify a value for \"-d, --DataField\" option in \"DataField\" \"-m, --mode\". \n";
    }
    $OptionsInfo{DataField} = $Options{datafield};
  }

  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? 1 : 0;

  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : 0;
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{mode} = 'RootPrefix';

  if (!GetOptions(\%Options, "datafield|d=s", "help|h", "mode|m=s", "overwrite|o", "root|r=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }

  if ($Options{mode} !~ /^(DataField|MolName|RootPrefix)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: DataField, MolName, RootPrefix\n";
  }
}

__END__

=head1 NAME

SDToMolFiles.pl - Generate MDLMOL file(s) from SD file(s)

=head1 SYNOPSIS

SDToMolFiles.pl SDFile(s)...

SDToMolFiles.pl [B<-d, --DataField> DataFieldName]
[B<-m, --mode> DataField | MolName | RootPrefix] [B<-h, --help>]
[B<-o, --overwrite>] [B<-r, --root> rootname]
[B<-w, --workingdir> dirname] SDFile(s)...

=head1 DESCRIPTION

Generate MDLMOL file(s) from I<SDFile(s)>. All header data labels and values in
SDFile(s) are simply ignored; other appopriate data from SDFile(s) is transferred to MDLMOL
files. Multiple I<SDFile(s)> names are separated by spaces. The valid file extensions are
I<.sdf> and I<.sd>. All other file names are ignored. All the SD files in a current
directory can be specified either by I<*.sdf> or the current directory name.

=head1 OPTIONS

=over 4

=item B<-d, --DataField> I<DataFieldName>

Specify I<SDFile(s)> datafield label name whose value is used for generation of MDLMOL
file names. Default value: I<None>.

=item B<-h, --help>

Print this help message.

=item B<-m, --mode> I<DataField | MolName | RootPrefix>

Specify how to generate MDLMOL file names: use a I<SDFile(s)> datafield value; use
molname line from I<SDFile(s)>; generate a sequential ID using root prefix specified
by B<-r, --root> option.

Possible values: I<DataField | MolName | RootPrefix | RootPrefix>.
Default: I<RootPrefix>.

For empty I<MolName> and I<DataField> values during these specified modes, file
name is automatically generated using I<RootPrefix>.

For I<RootPrefix> value of B<-m, --mode> option, MDLMOL file names are generated
using by appending compound record number to value of B<-r, --root> option. For
example: I<RootName>Cmd<RecordNumber>.mol.

Allowed characters in file names are: a-zA-Z0-9_. All other characters in datafield
values, molname line, and root prefix are ignore during generation of file names.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

Specify root name to used during I<RootPrefix> B<-m, --mode> option value.
New MDLMOL file names are generated using the root: <Root>Cmpd<RecordNumber>.mol
Default for new file names: <InitialSDFileName>Cmpd<RecordNumber>.mol. This option
is ignored for multiple input files.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To generate MDLMOL files from Sample1*.sdf and Sample2*.sd files, type:

    % SDToMolFiles.pl -o Sample1*.sdf Sample2*.sd

To generate Sample*.mol files from Sample1.sdf, type:

    % SDToMolFiles.pl -r Sample -o Sample1.sdf

To generate MOL files from Sample1.sdf using molname line data for generating
MOL file names, type:

    % SDToMolFiles.pl -m MolName -r Sample -o Sample1.sdf

To generate MOL files from Sample1.sdf using a specific data field values for
generating MOL file names, type:

    % SDToMolFiles.pl -m DataField --DataField MolID -r Sample
      -o Sample1.sdf

=head1 AUTHOR


=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoSDFiles.pl, MolFilesToSD.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
