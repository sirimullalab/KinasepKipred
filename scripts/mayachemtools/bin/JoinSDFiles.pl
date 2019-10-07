#!/usr/bin/perl -w
#
# File: JoinSDFiles.pl
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
if (@SDFilesList == 1) {
  die "Error: Specify more than one SD file.\n";
}

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Setup information about input files...
print "Checking input SD files...\n";
my(%SDFilesInfo);
RetrieveSDFilesInfo();

# Join files...
print "\nGenerating new SD file $OptionsInfo{NewSDFile}...\n";
JoinSDFiles();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Join all valid SD files...
sub JoinSDFiles {
  my($FileIndex, $SDFile, $NewSDFile);

  $NewSDFile = $OptionsInfo{NewSDFile};

  open NEWSDFILE, ">$NewSDFile" or die "Error: Couldn't open $NewSDFile: $! \n";
  FILELIST: for $FileIndex (0 .. $#SDFilesList) {
    if (!$SDFilesInfo{FileOkay}[$FileIndex]) {
      next FILELIST;
    }
    $SDFile = $SDFilesList[$FileIndex];
    print "\nProcessing file $SDFile...\n";

    open SDFILE, "$SDFile" or die "Error: Couldn't open $SDFile: $! \n";
    while (<SDFILE>) {
      s/(\r\n)|(\r)/\n/g;
      print NEWSDFILE;
    }
    close SDFILE;
  }

  close NEWSDFILE;
}

# Retrieve information about SD files...
sub RetrieveSDFilesInfo {
  my($Index, $SDFile);

  %SDFilesInfo = ();
  @{$SDFilesInfo{FileOkay}} = ();

  FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFilesInfo{FileOkay}[$Index] = 0;

    $SDFile = $SDFilesList[$Index];
    if (!(-e $SDFile)) {
      warn "Warning: Ignoring file $SDFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($SDFile, "sdf sd")) {
      warn "Warning: Ignoring file $SDFile: It's not a SD file\n";
      next FILELIST;
    }
    if (! open SDFILE, "$SDFile") {
      warn "Warning: Ignoring file $SDFile: Couldn't open it: $! \n";
      next FILELIST;
    }
    close SDFILE;

    $SDFilesInfo{FileOkay}[$Index] = 1;
  }
}

# Process option values...
sub ProcessOptions {
  my($FileDir, $FileName, $FileExt, $NewSDFile);

  %OptionsInfo = ();

  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : undef;
  $OptionsInfo{Overwrite} = $Options{overwrite} ? $Options{overwrite} : undef;

  if ($Options{root}) {
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($Options{root});
    if ($FileName && $FileExt) {
      $NewSDFile = $FileName . "." . $FileExt;
    }
    else {
      $NewSDFile =  $Options{root} . ".sdf";
    }
  }
  else {
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($SDFilesList[0]);
    $NewSDFile = $FileName . "1To" . @SDFilesList . "Joined.sdf";
  }

  if (!$Options{overwrite}) {
    if (-e $NewSDFile) {
      die "Error: The file $NewSDFile already exists.\n";
    }
  }
  if ($Options{root}) {
    my($FileIndex);
    for $FileIndex (0 .. $#SDFilesList) {
      if (lc($NewSDFile) eq lc($SDFilesList[$FileIndex])) {
	die "Error: Output filename, $NewSDFile, is similar to a input file name.\nSpecify a different name using \"-r --root\" option or use default name.\n";
      }
    }
  }
  $OptionsInfo{NewSDFile} = $NewSDFile;

}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  if (!GetOptions(\%Options, "help|h", "overwrite|o", "root|r=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
}

__END__

=head1 NAME

JoinSDFiles.pl - Join multiple SDFiles into a single SDFile

=head1 SYNOPSIS

JoinSDFiles.pl  SDFiles...

JoinSDFiles.pl [B<-h, --help>] [B<-o, --overwrite>] [B<-r, --root> rootname]
[B<-w, --workingdir> dirname] SDFiles...

=head1 DESCRIPTION

Multiple I<SDFiles> are joined to generate a single SDFile. The file names
are separated by spaces. The valid file extensions are I<.sdf> and I<.sd>.
All other file names are ignored. All the SD files in a current directory can be
specified either by I<*.sdf> or the current directory name.

=head1 OPTIONS

=over 4

=item B<-h, --help>

Print this help message.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

New SD file name is generated using the root: <Root>.sdf. Default file
name:<FirstSDFileName>1To<Count>Joined.sdf.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To join SD files, type:

    % JoinSDFiles.pl -o Sample1.sdf Sample2.sdf
    % JoinSDFiles.pl -o *.sdf

To join all Sample*.sdf files in a directory, SomeDir, and generate a new file NewSample.sdf, type:

    % JoinSDFiles.pl -r NewSample -w SomeDir -o *.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoSDFiles.pl, MolFilesToSD.pl, SDToMolFiles.pl, SplitSDFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
