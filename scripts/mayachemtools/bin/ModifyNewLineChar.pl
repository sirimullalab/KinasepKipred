#!/usr/bin/perl -w
#
# File: ModifyNewLineChar.pl
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
use TextUtil;

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

my( @FilesList);
@FilesList = ExpandFileNames(\@ARGV, "");

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

print "Checking input file(s)...\n";
my(%FilesInfo);
RetrieveFilesInfo();

# Generate output files...
my($FileIndex);
if (@FilesList > 1) {
  print "\nProcessing files...\n";
}
for $FileIndex (0 .. $#FilesList) {
  if ($FilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $FilesList[$FileIndex]...\n";
    ModifyNewLineChar($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Modify new line characters...
sub ModifyNewLineChar {
  my($Index) = @_;
  my($File, $NewFile, $Mode, $Nothing);

  $File = $FilesList[$Index];
  $NewFile = $FilesInfo{OutFile}[$Index];

  $Mode = $OptionsInfo{Mode};

  print "Generating new $NewFile file...\n";

  open NEWFILE, ">$NewFile" or die "Error: Can't open $NewFile: !$ \n";
  open FILE, "$File" or die "Error: Can't open $File: $! \n";

  while (<FILE>) {
    LINE: {
      if ($Mode =~ /^Unix$/i) { s/(\r\n)|(\r)|(\n)/\n/g; last LINE; }
      if ($Mode =~ /^Windows$/i) { s/(\r\n)|(\r)|(\n)/\r\n/g; last LINE; }
      if ($Mode =~ /^Mac$/i) { s/(\r\n)|(\r)|(\n)/\r/g; last LINE; }
      $Nothing = 1;
    }
    print NEWFILE;
  }

  close NEWFILE;
  close FILE;
}

# Retrieve input files info...
sub RetrieveFilesInfo {
  my($File, $Index, $FileDir, $FileName, $FileExt, $NewFileName);

  %FilesInfo = ();

  @{$FilesInfo{FileOkay}} = ();
  @{$FilesInfo{OutFile}} = ();

  FILELIST: for $Index (0 .. $#FilesList) {
    $File = $FilesList[$Index];

    $FilesInfo{FileOkay}[$Index] = 0;
    $FilesInfo{OutFile}[$Index] = "";

    if (!(-e $File)) {
      warn "Warning: Ignoring file $File: It doesn't exist\n";
      next FILELIST;
    }

    if (!open FILE, "$File") {
      warn "Warning: Ignoring file $File: Couldn't open it: $! \n";
      next FILELIST;
    }
    close FILE;

    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($File);
    $NewFileName = $FileName;
    if ($OptionsInfo{OutFileRoot} && (@FilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($OptionsInfo{OutFileRoot});
      if ($RootFileName && $RootFileExt) {
	$NewFileName = $RootFileName;
      }
      else {
	$NewFileName = $OptionsInfo{OutFileRoot};
      }
    }
    else {
      $NewFileName .= $OptionsInfo{Mode};
    }

    if ($FileExt) {
      $NewFileName .= ".$FileExt";
    }

    if (!$OptionsInfo{Overwrite}) {
      if (-e $NewFileName) {
	warn "Warning: Ignoring file $File: New Text file, $NewFileName, already exists\n";
	next FILELIST;
      }
    }
    $FilesInfo{FileOkay}[$Index] = 1;
    $FilesInfo{OutFile}[$Index] = $NewFileName;
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : undef;
  $OptionsInfo{Overwrite} = $Options{overwrite} ? $Options{overwrite} : undef;

}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{mode} = "Unix";

  if (!GetOptions(\%Options, "help|h", "mode|m=s", "overwrite|o", "root|r=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir},  for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{mode} !~ /^(Unix|Windows|Mac)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: Unix, Windows, or Mac\n";
  }
}

__END__

=head1 NAME

ModifyNewLineChar.pl - Modify new line char(s)

=head1 SYNOPSIS

ModifyNewLineChar.pl File(s)...

ModifyNewLineChar.pl [B<-h, --help>] [B<-m, --mode> Unix | Mac | Windows] [B<-o, --overwrite>]
[B<-r, --root> rootname] [B<-w, --workingdir> dirname] File(s)...

=head1 DESCRIPTION

Modify new line char(s) in ASCII files to interchange among Unix, Windows, and Mac
formats.

=head1 OPTIONS

=over 4

=item B<-h, --help>

Print this help message.

=item B<-m, --mode> I<Unix | Mac | Windows>

New line  char(s) mode. Possible values: I<Unix, Mac, or Windows>. Default: I<Unix>. Here
are default values for new line char(s): I<Unix - \n; Windows: \r\n; Mac - \r>

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

New text file name is generated using the root: <Root>.<Ext>. Default new file name:
<InitialFileName><Mode>.<InitialFileExt>. This option is ignored for multiple input files.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To use Unix new line char and generate NewSample1.csv file, type:

    % ModifyNewLineChar.pl -m Unix -r NewSample1 -o Sample1.csv

To use Mac new line char and generate NewSample1.sdf file, type:

    % ModifyNewLineChar.pl -m Mac -r NewSample1 -o Sample1.sdf

To use Windows new line chars and generate NewSample1.csv file, type:

    % ModifyNewLineChar.pl -m Windows -r NewSample1 -o Sample1.csv

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

ModifyTextFilesFormat.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
