#!/usr/bin/perl -w
#
# File: MolFilesToSD.pl
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

my(@MOLFilesList);
@MOLFilesList = ExpandFileNames(\@ARGV, "mol");

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

print "Generating SD file $OptionsInfo{SDFile}...\n";
GenerateSDFile();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Generate a SD file using all valid MDL MOL files...
sub GenerateSDFile {
  my($MOLFile, $Index, $FileCount, $FileOkayCount, $MolNameLine, $CmpdID, $FileDir, $FileName, $FileExt);

  open SDFILE, ">$OptionsInfo{SDFile}" or die "Error: Can't open $OptionsInfo{SDFile}: $! \n";
  $FileCount = 0;
  $FileOkayCount = 0;

  FILELIST: for $Index (0 .. $#MOLFilesList) {
    $MOLFile = $MOLFilesList[$Index];
    $FileCount++;

    print "Processing file $MOLFile...\n";

    if (!(-e $MOLFile)) {
      warn "Warning: Ignoring file $MOLFile: It doesn't exist\n";
      next FILELIST;
    }

    if (!CheckFileType($MOLFile, "mol")) {
      warn "Warning: Ignoring file $MOLFile: It's not a MDLMOL file\n";
      next FILELIST;
    }

    if (!open MOLFILE, "$MOLFile") {
      warn "Warning: Ignoring file $MOLFile: Couldn't open it: $! \n";
      next FILELIST;
    }

    $FileOkayCount++;

    if ($OptionsInfo{ModifyData}) {
      $MolNameLine = <MOLFILE>;
      if ($OptionsInfo{UseFilePrefix}) {
	($FileDir, $FileName, $FileExt) = ParseFileName($MOLFile);
	$CmpdID = $FileName;
      }
      else {
	$CmpdID = $OptionsInfo{CompoundID} . "$FileCount";
      }

      if ($OptionsInfo{AddMolNameLine}) {
	print SDFILE "$CmpdID\n";
      }
      else {
	$MolNameLine =~ s/(\r\n)|(\r)/\n/g;
	print SDFILE $MolNameLine;
      }

      while (<MOLFILE>) {
	s/(\r\n)|(\r)/\n/g;
	print SDFILE;
      }

      if ($OptionsInfo{AddDataField}) {
	print SDFILE ">  <$OptionsInfo{DataFieldLabel}>\n${CmpdID}\n";
      }
    }
    else {
      while (<MOLFILE>) {
	s/(\r\n)|(\r)/\n/g;
	print SDFILE;
      }
    }
    print SDFILE "\n\$\$\$\$\n";
    close MOLFILE;
  }
  close SDFILE;

  print "\nNumber of files: $FileCount\n";
  print "Number of files processed successfully: $FileOkayCount\n";
  print "Number of files ignored: " . ($FileCount - $FileOkayCount) . "\n";
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{CompoundID} = $Options{compoundid};
  $OptionsInfo{DataFieldLabel} = $Options{datafieldlabel};

  $OptionsInfo{Overwrite} = defined $Options{overwrite} ? $Options{overwrite} : undef;
  $OptionsInfo{OutFileRoot} = defined $Options{root} ? $Options{root} : undef;

  $OptionsInfo{AddMolNameLine} = ($Options{mode} =~ /^(molnameline|both)$/i) ? 1 : 0;
  $OptionsInfo{AddDataField} = ($Options{mode} =~ /^(datafield|both)$/i) ? 1 : 0;

  $OptionsInfo{AddMolNameLine} = ($Options{mode} =~ /^(molnameline|both)$/i) ? 1 : 0;
  $OptionsInfo{AddDataField} = ($Options{mode} =~ /^(datafield|both)$/i) ? 1 : 0;

  $OptionsInfo{ModifyData} = ($OptionsInfo{AddMolNameLine} || $OptionsInfo{AddDataField}) ? 1 : 0;

  $OptionsInfo{UseFilePrefix} = ($Options{compoundid} =~ /^usefileprefix$/i) ? 1 : 0;

  # Setup SD file name...
  my($FileDir, $FileName, $FileExt, $SDFile);
  if ($Options{root}) {
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($Options{root});
    if ($FileName && $FileExt) {
      $SDFile = $FileName;
    }
    else {
      $SDFile =  $Options{root};
    }
    $SDFile .=  ".sdf";
  }
  else {
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($MOLFilesList[0]);
    $SDFile = $FileName . "1To" . @MOLFilesList . ".sdf";
  }

  if (!$Options{overwrite}) {
    if (-e $SDFile) {
      die "Error: The file $SDFile already exists.\n";
    }
  }
  $OptionsInfo{SDFile} = $SDFile;

}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{compoundid} = "Cmpd";
  $Options{datafieldlabel} = "Cmpd_ID";
  $Options{mode} = "none";

  if (!GetOptions(\%Options, "compoundid|c=s", "datafieldlabel|d=s", "help|h", "mode|m=s", "overwrite|o", "root|r=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{mode} !~ /^(molnameline|datafield|both|none)$/i ) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: molnameline, datafield, both, or none\n";
  }
}

__END__

=head1 NAME

MolFilesToSD.pl - Generate a SD file from MDLMOL File(s)

=head1 SYNOPSIS

MolFilesToSD.pl  MDLMOLFile(s)...

MolFilesToSD.pl [B<-c, --compoundid> usefileprefix | idlabel] [B<-d, --datafieldlabel> fieldlabel]
[B<-h, --help>] [B<-m, --mode> molnameline | datafield | both | none] [B<-o, --overwrite>]
[B<-r, --root> rootname] [B<-w, --workingdir> dirname] MDLMOLFile(s)...

=head1 DESCRIPTION

Generate a SD file from I<MDLMOL File(s)>. Multiple file names are separated by spaces.
The valid file extension is I<.mol>. All other file names are ignored. All the files in a current
directory can be specified by I<*.mol>, or the current directory name.

=head1 OPTIONS

=over 4

=item B<-c, --compoundid> I<usefileprefix | idlabel>

Specify how to generate compound IDs: use MOL filename prefix or generate
a new compound ID by combining I<idlabel> with compound number. Possible
values: I<usefileprefix | idlabel>. By default, I<Cmd> is used as a I<idlabel> to generate
these types of compound IDs: Cmpd1, Cmpd2 and so on.

Example: To generate compound IDs like Mol_ID1, Mol_ID2 and so on, specify
"MolID" value for this option.

=item B<-d, --datafieldlabel> I<fieldlabel>

Specify data field label for adding compound ID field into SD file during I<datafield | both>
values of B<-m, --mode> option. Default: <Cmpd_ID>.

=item B<-h, --help>

Print this help message.

=item B<-m, --mode> I<molnameline | datafield | both | none>

Specify how to add compopund ID into SD file: relplace the molname line,
add a new data field, replace the molname line and add data field, or do
nothing. Possible values: I<molnameline | datafield | both | none>.
Default: I<nothing>.

Use B<-c, --compoundid> to specify compound ID generation process.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

New SD file name is generated using the root: <Root>.sdf. Default new file
name: <InitialMOLFileName>1To<Count>.sdf.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To generate NewSample.sdf file from Sample*.mol files, type:

    % MolFilesToSD.pl  -r NewSample -o Sample*.mol

To generate NewSample.sdf with Cmpd1, Cmpd2 and so on as compound ID in
MolName line and Cmpd_ID datafield  from Sample*.mol files, type:

    % MolFilesToSD.pl  -r NewSample -m both -o Sample*.mol

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoSDFiles.pl, SDToMolFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
