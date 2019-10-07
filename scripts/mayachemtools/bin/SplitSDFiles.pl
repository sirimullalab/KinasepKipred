#!/usr/bin/perl -w
#
# File: SplitSDFiles.pl
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
my(%SDFilesInfo);
print "Checking input SD file(s)...\n";
RetrieveSDFilesInfo();

# Process input files..
my($FileIndex);
if (@SDFilesList > 1) {
  print "\nProcessing SD files...\n";
}
for $FileIndex (0 .. $#SDFilesList) {
  if ($SDFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $SDFilesList[$FileIndex]...\n";
    SplitSDFile($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Split a SD file...
#
sub SplitSDFile {
  my($FileIndex) = @_;

  if ($OptionsInfo{Mode} =~ /^Files$/i) {
    SplitSDFileByNumOfFiles($FileIndex);
  }
  elsif ($OptionsInfo{Mode} =~ /^Cmpds$/i) {
    SplitSDFileByNumOfCmpds($FileIndex);
  }
}

# Split SD into specified number of files...
#
sub SplitSDFileByNumOfFiles {
  my($FileIndex) = @_;
  my($SDFile, $CmpdCount, $MaxCmpdsPerFile, $MaxNumOfFiles);

  $SDFile = $SDFilesList[$FileIndex];

  if (!open SDFILE, "$SDFile") {
    warn "Warning: Ignoring file $SDFile: Couldn't open it: $! \n";
    return;
  }

  $MaxNumOfFiles = $OptionsInfo{NumOfFiles};

  # Count number of compounds to figure out maximum number of compound per file...
  $CmpdCount = 0;
  while (<SDFILE>) {
    if (/^\$\$\$\$/) {
      $CmpdCount++;
    }
  }
  close SDFILE;

  if ($CmpdCount < $MaxNumOfFiles) {
    warn "Warning: Ignoring file $SDFile: Total number of compounds, $CmpdCount, is smaller than number of new files, $MaxNumOfFiles\n";
    return;
  }

  $MaxCmpdsPerFile = int $CmpdCount / $MaxNumOfFiles;

  SplitSDFileByNumOfFilesAndCmpds($FileIndex, $MaxNumOfFiles, $MaxCmpdsPerFile);
}

# Split SD into files containing specified number of compounds...
#
sub SplitSDFileByNumOfCmpds {
  my($FileIndex) = @_;

  if ($OptionsInfo{NumOfCmpds} == 1) {
    SplitSDFileByOneCmpdPerFile($FileIndex);
  }
  else {
    SplitSDFileByNumOfCmpdsPerFile($FileIndex);
  }
}

# Split SD into files containing one compound per file...
#
sub SplitSDFileByOneCmpdPerFile {
  my($FileIndex) = @_;
  my($SDFile, $NewSDFile, $NewSDFileRoot, $FileExt, $OutFileRoot, $OverwriteFiles, $UseDataField, $DataFieldName, $UseMolName, $CmpdCount, $CmpdString, @CmpdLines, %DataFieldValues);

  $SDFile = $SDFilesList[$FileIndex];

  if (!open SDFILE, "$SDFile") {
    warn "Warning: Ignoring file $SDFile: Couldn't open it: $! \n";
    return;
  }

  print "\n";

  $CmpdCount = 0;

  $FileExt = $SDFilesInfo{FileExt}[$FileIndex];

  $OutFileRoot = $SDFilesInfo{OutFileRoot}[$FileIndex];
  $OverwriteFiles = $OptionsInfo{OverwriteFiles};

  $UseDataField = ($OptionsInfo{CmpdsMode} =~ /^DataField$/i) ? 1 : 0;
  $DataFieldName = $OptionsInfo{DataField};

  $UseMolName = ($OptionsInfo{CmpdsMode} =~ /^MolName$/i) ? 1 : 0;

  CMPDSTRING: while ($CmpdString = ReadCmpdString(\*SDFILE)) {
    $CmpdCount++;

    # Setup SD file name...
    $NewSDFileRoot = '';
    if ($UseDataField) {
      @CmpdLines = split "\n", $CmpdString;
      %DataFieldValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
      if (exists $DataFieldValues{$DataFieldName}) {
	$NewSDFileRoot = $DataFieldValues{$DataFieldName};
      }
    }
    elsif ($UseMolName) {
      @CmpdLines = split "\n", $CmpdString;
      $NewSDFileRoot = $CmpdLines[0];
    }

    # Check for any invalid file name characters in data field or molname values...
    if ($NewSDFileRoot && $NewSDFileRoot =~ /[^a-zA-Z0-9_]/) {
      $NewSDFileRoot =~ s/[^a-zA-Z0-9_]//g;
    }

    # Fall back plan for SD file name...
    if (!$NewSDFileRoot) {
      $NewSDFileRoot = "${OutFileRoot}Cmpd${CmpdCount}";
    }

    $NewSDFile = "${NewSDFileRoot}.${FileExt}";

    if (!$OverwriteFiles) {
      if (-e $NewSDFile) {
	warn "Warning: Ignoring compound number, $CmpdCount, in $SDFile: New SD file, $NewSDFile, already exists\n";
	next CMPDSTRING;
      }
    }

    # Write out new SD file...

    print "Generating $NewSDFile file\n";
    open NEWSDFILE, ">$NewSDFile" or die "Error: Can't open $NewSDFile: $! \n";
    print NEWSDFILE "$CmpdString\n";
    close NEWSDFILE;

  }
  close SDFILE;
}

# Split SD into files containing specified number of compounds per file...
#
sub SplitSDFileByNumOfCmpdsPerFile {
  my($FileIndex) = @_;
  my($SDFile, $CmpdCount, $MaxCmpdsPerFile, $MaxNumOfFiles);

  $SDFile = $SDFilesList[$FileIndex];

  if (!open SDFILE, "$SDFile") {
    warn "Warning: Ignoring file $SDFile: Couldn't open it: $! \n";
    return;
  }

  $MaxCmpdsPerFile = $OptionsInfo{NumOfCmpds};

  # Count number of compounds to figure out maximum number of files...
  $CmpdCount = 0;
  while (<SDFILE>) {
    if (/^\$\$\$\$/) {
      $CmpdCount++;
    }
  }
  close SDFILE;

  $MaxNumOfFiles = int $CmpdCount / $MaxCmpdsPerFile;

  if (($MaxNumOfFiles * $MaxCmpdsPerFile) < $CmpdCount) {
    $MaxNumOfFiles++;
  }

  if ($CmpdCount <= $MaxCmpdsPerFile) {
    warn "Warning: Ignoring file $SDFile: Total number of compounds, $CmpdCount, is <= specified number of compunds per file, $MaxCmpdsPerFile\n";
    return;
  }

  SplitSDFileByNumOfFilesAndCmpds($FileIndex, $MaxNumOfFiles, $MaxCmpdsPerFile);
}

# Split SD files into specified number of files with specified number of compounds
# in each file...
#
sub SplitSDFileByNumOfFilesAndCmpds {
  my($FileIndex, $NumOfFiles, $NumOfCmpdsPerFile) = @_;
  my($SDFile, $CmpdCount, $NewFileIndex, $NewFileName, $MaxCmpdsCount, @NewSDFilesList);

  $SDFile = $SDFilesList[$FileIndex];

  if (!open SDFILE, "$SDFile") {
    warn "Warning: Ignoring file $SDFile: Couldn't open it: $! \n";
    return;
  }

  # Setup new file names list...
  @NewSDFilesList = ();
  for $NewFileIndex (1 .. $NumOfFiles) {
    $NewFileName = $SDFilesInfo{OutFileRoot}[$FileIndex] . "Part${NewFileIndex}." . $SDFilesInfo{FileExt}[$FileIndex];
    if (!$OptionsInfo{OverwriteFiles}) {
      if (-e $NewFileName) {
	warn "Warning: Ignoring file $SDFile: New SD file, $NewFileName, already exists\n";
	return;
      }
    }
    push @NewSDFilesList, $NewFileName;
  }

  $MaxCmpdsCount = $NumOfCmpdsPerFile;

  $CmpdCount = 0;
  $NewFileIndex = 1;

  open NEWSDFILE, ">$NewSDFilesList[$NewFileIndex - 1]" or die "Error: Can't open $NewSDFilesList[$NewFileIndex -1]: $! \n";
  print "\nGenerating $NewSDFilesList[$NewFileIndex - 1] file\n";

  open SDFILE, "$SDFile" or die "Error: Can't open $SDFile: $! \n";

  while (<SDFILE>) {
    s/(\r\n)|(\r)/\n/g;
    print NEWSDFILE;

    if ( /^\$\$\$\$/ ) {
      $CmpdCount++;
      if ($NewFileIndex <= $NumOfFiles) {
	if ($CmpdCount >= $MaxCmpdsCount) {
	  if ($NewFileIndex < $NumOfFiles) {
	    close NEWSDFILE;
	  }
	  $NewFileIndex++;
	  $MaxCmpdsCount = $NumOfCmpdsPerFile * $NewFileIndex;

	  if ($NewFileIndex <= $NumOfFiles) {
	    open NEWSDFILE, ">$NewSDFilesList[$NewFileIndex - 1]" or die "Error: Can't open $NewSDFilesList[$NewFileIndex - 1]: $! \n";
	    print "Generating $NewSDFilesList[$NewFileIndex - 1] file...\n";
	  }
	}
      }
    }
  }
  close NEWSDFILE;
}

# Retrieve information about SD files...
#
sub RetrieveSDFilesInfo {
  my($SDFile, $Index, $FileDir, $FileName, $FileExt, $OutFileRoot);

  %SDFilesInfo = ();
  @{$SDFilesInfo{FileOkay}} = ();
  @{$SDFilesInfo{FileExt}} = ();
  @{$SDFilesInfo{OutFileRoot}} = ();

  FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFile = $SDFilesList[$Index];

    $SDFilesInfo{FileOkay}[$Index] = 0;
    $SDFilesInfo{FileExt}[$Index] = '';
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
    $SDFilesInfo{FileExt}[$Index] = $FileExt;
    $SDFilesInfo{OutFileRoot}[$Index] = $OutFileRoot;
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{CmpdsMode} = $Options{cmpdsmode};

  $OptionsInfo{NumOfFiles} = $Options{numfiles};
  $OptionsInfo{NumOfCmpds} = $Options{numcmpds};

  $OptionsInfo{DataField} = '';
  if ($Options{mode} =~ /^Cmpds$/i && $Options{cmpdsmode} =~ /^DataField$/i) {
    if (!$Options{datafield}) {
      die "Error: You must specify a value for \"-d, --DataField\" option in \"DataField\" value of \"-c, --CmpdsMode\" during \"Cmpds\" \"-m, --mode\" value. \n";
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

  $Options{cmpdsmode} = 'RootPrefix';
  $Options{mode} = 'Files';

  $Options{numfiles} = 2;
  $Options{numcmpds} = 1;


  if (!GetOptions(\%Options, "cmpdsmode|c=s", "datafield|d=s", "help|h", "mode|m=s", "numfiles|n=i", "numcmpds=i", "overwrite|o", "root|r=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{cmpdsmode} !~ /^(DataField|MolName|RootPrefix)$/i) {
    die "Error: The value specified, $Options{cmpdsmode}, for option \"-c, --CmpdsMode\" is not valid. Allowed values: DataField, MolName, RootPrefix\n";
  }
  if ($Options{mode} !~ /^(Cmpds|Files)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m, --mode\" is not valid. Allowed values: Cmpds, Files\n";
  }
  if ($Options{numfiles} < 2) {
    die "Error: The value specified, $Options{numfiles}, for option \"-n --numfiles\" is not valid. Allowed values: >= 2 \n";
  }
  if ($Options{numcmpds} < 1) {
    die "Error: The value specified, $Options{numcmpds}, for option \"-n --numcmpds\" is not valid. Allowed values: >= 1 \n";
  }
}

__END__

=head1 NAME

SplitSDFiles.pl - Split SDFile(s) into multiple SD files

=head1 SYNOPSIS

SplitSDFiles.pl SDFile(s)...

SplitSDFiles.pl [B<-c, --CmpdsMode> DataField | MolName | RootPrefix]
[B<-d, --DataField> DataFieldName] [B<-h, --help>] [B<-m, --mode> Cmpds | Files]
[B<-n, --numfiles> number] [B<--numcmpds> number] [B<-o, --overwrite>]
[B<-r, --root> rootname] [B<-w,--workingdir> dirname] SDFile(s)...

=head1 DESCRIPTION

Split I<SDFile(s)> into multiple SD files. Each new SDFile contains a compound
subset of similar size from the initial file. Multiple I<SDFile(s)> names are separated
by space. The valid file extensions are I<.sdf> and I<.sd>. All other file names are
ignored. All the SD files in a current directory can be specified either by I<*.sdf>
or the current directory name.

=head1 OPTIONS

=over 4

=item B<-c, --CmpdsMode> I<DataField | MolName | RootPrefix>

This option is only used during I<Cmpds> value of <-m, --mode> option with
specified B<--numcmpds> value of 1.

Specify how to generate new file names during I<Cmpds> value of <-m, --mode>
option: use I<SDFile(s)> datafield value or molname line for a specific compound;
generate a sequential ID using root prefix specified by B<-r, --root> option.

Possible values: I<DataField | MolName | RootPrefix | RootPrefix>.
Default: I<RootPrefix>.

For empty I<MolName> and I<DataField> values during these specified modes, file
name is automatically generated using I<RootPrefix>.

For I<RootPrefix> value of B<-c, --CmpdsMode> option, new file names are
generated using by appending compound record number to value of B<-r, --root> option.
For example: I<RootName>Cmd<RecordNumber>.sdf.

Allowed characters in file names are: a-zA-Z0-9_. All other characters in datafield
values, molname line, and root prefix are ignore during generation of file names.

=item B<-d, --DataField> I<DataFieldName>

This option is only used during I<DataField> value of <-c, --CmpdsMode> option.

Specify I<SDFile(s)> datafield label name whose value is used for generation of new file
for a specific compound. Default value: I<None>.

=item B<-h, --help>

Print this help message.

=item B<-m, --mode> I<Cmpds | Files>

Specify how to split I<SDFile(s)>: split into files with each file containing specified
number of compounds or split into a specified number of files.

Possible values: I<Cmpds | Files>. Default: I<Files>.

For I<Cmpds> value of B<-m, --mode> option, value of B<--numcmpds> option
determines the number of new files. And value of B<-n, --numfiles> option is
used to figure out the number of new files for I<Files> value of B<-m, --mode> option.

=item B<-n, --numfiles> I<number>

Number of new files to generate for each I<SDFile(s)>. Default: I<2>.

This value is only used during I<Files> value of B<-m, --mode> option.

=item B<--numcmpds> I<number>

Number of compounds in each new file corresponding to each I<SDFile(s)>.
Default: I<1>.

This value is only used during I<Cmpds> value of B<-m, --mode> option.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

New SD file names are generated using the root: <Root>Part<Count>.sdf.
Default new file names: <InitialSDFileName> Part<Count>.sdf. This option
is ignored for multiple input files.

=item B<-w,--workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To split each SD file into 5 new SD files, type:

    % SplitSDFiles.pl -n 5 -o Sample1.sdf Sample2.sdf
    % SplitSDFiles.pl -n 5 -o *.sdf

To split Sample1.sdf into 10 new NewSample*.sdf files, type:

    % SplitSDFiles.pl -m Files -n 10 -r NewSample -o Sample1.sdf

To split Sample1.sdf into new NewSample*.sdf files containing maximum of 5 compounds
in each file, type:

    % SplitSDFiles.pl -m Cmpds --numcmpds 5 -r NewSample -o Sample1.sdf

To split Sample1.sdf into new SD files containing one compound each with new file
names corresponding to molname line, type:

    % SplitSDFiles.pl -m Cmpds --numcmpds 1 -c MolName -o Sample1.sdf

To split Sample1.sdf into new SD files containing one compound each with new file
names corresponding to value of datafield MolID, type:

    % SplitSDFiles.pl -m Cmpds --numcmpds 1 -c DataField -d MolID
      -o Sample1.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoSDFiles.pl, JoinSDFiles.pl, MolFilesToSD.pl, SDToMolFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
