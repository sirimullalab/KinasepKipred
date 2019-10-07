#!/usr/bin/perl -w
#
# File: FilterSDFiles.pl
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

print "Checking input SD file(s)...\n";
my(%SDFilesInfo);
RetrieveSDFilesInfo();

# Generate output files...
my($FileIndex, %FilteredSDFileInfo);
if (@SDFilesList > 1) {
  print "\nProcessing SD files...\n";
}
for $FileIndex (0 .. $#SDFilesList) {
  if ($SDFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $SDFilesList[$FileIndex]...\n";
    FilterSDFile($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Filter SD file...
sub FilterSDFile {
  my($Index) = @_;
  my($SDFile, $NewSDFile, $NewKeepSDFile, $CtabLinesCount, $CmpdString, $PrintCmpdCounterHeader, @CmpdLines);

  $SDFile = $SDFilesList[$Index];
  $NewSDFile = $SDFilesInfo{OutFile}[$Index];
  $NewKeepSDFile = $SDFilesInfo{OutFileKeep}[$Index];

  open NEWSDFILE, ">$NewSDFile" or die "Error: Couldn't open $NewSDFile: $! \n";
  if ($OptionsInfo{Keep}) {
    open NEWKEEPSDFILE, ">$NewKeepSDFile" or die "Error: Couldn't open $NewKeepSDFile: $! \n";
  }
  open SDFILE, "$SDFile" or die "Error: Can't open $SDFile: $! \n";

  print "\nGenerating SD file $NewSDFile...\n";
  if ($OptionsInfo{Keep}) {
    print "Generating file $NewKeepSDFile...\n";
  }

  %FilteredSDFileInfo = ();

  $FilteredSDFileInfo{CmpdCount} = 0; $FilteredSDFileInfo{FilterCmpd} = 0;
  $FilteredSDFileInfo{FilteredCmpdCount} = 0; $FilteredSDFileInfo{KeepCmpdCount} = 0;

  $PrintCmpdCounterHeader = 1;

  CMPDSTRING: while ($CmpdString = ReadCmpdString(\*SDFILE)) {
    $FilteredSDFileInfo{CmpdCount} += 1;
    $FilteredSDFileInfo{FilterCmpd} = 0;
    if (($FilteredSDFileInfo{CmpdCount} % 5000) == 0) {
      if ($PrintCmpdCounterHeader) {
	$PrintCmpdCounterHeader = 0;
	print "\nProcessing compounds:";
      }
      print "$FilteredSDFileInfo{CmpdCount}...";
    }
    @CmpdLines = split "\n", $CmpdString;
    $CtabLinesCount = GetCtabLinesCount(\@CmpdLines);
    if ($CtabLinesCount <= 0) {
      $FilteredSDFileInfo{FilterCmpd} = 1;
      WriteOutCmpdString($CmpdString);
      next CMPDSTRING;
    }
    my ($AtomCount, $BondCount) = ParseCmpdCountsLine($CmpdLines[3]);
    if ($OptionsInfo{All} || $OptionsInfo{Mismatch}) {
      if ($CtabLinesCount != ($AtomCount + $BondCount)) {
	$FilteredSDFileInfo{FilterCmpd} = 1;
	WriteOutCmpdString($CmpdString);
	next CMPDSTRING;
      }
    }
    if ($CtabLinesCount == ($AtomCount + $BondCount)) {
      if ($OptionsInfo{All} || $OptionsInfo{UnknownAtoms}) {
	my($UnknownAtomCount, $UnknownAtoms, $UnknownAtomLines) = GetUnknownAtoms(\@CmpdLines);
	if ($UnknownAtomCount) {
	  $FilteredSDFileInfo{FilterCmpd} = 1;
	  WriteOutCmpdString($CmpdString);
	  next CMPDSTRING;
	}
      }
      if ($OptionsInfo{All} || $OptionsInfo{CleanSalts} || $OptionsInfo{Salts}) {
	my ($FragmentsCount, $Fragments, $WashedCmpdString) = WashCmpd(\@CmpdLines);
	if ($FragmentsCount > 1) {
	  if ($OptionsInfo{all} || $OptionsInfo{CleanSalts}) {
	    $CmpdString = $WashedCmpdString;
	  }
	  else {
	    $FilteredSDFileInfo{FilterCmpd} = 1;
	  }
	  WriteOutCmpdString($CmpdString);
	  next CMPDSTRING;
	}
      }
    }
    WriteOutCmpdString($CmpdString);
  }
  if (!$PrintCmpdCounterHeader) {
    print "\n";
  }

  close NEWSDFILE;
  if ($OptionsInfo{Keep}) {
    close NEWKEEPSDFILE;
  }
  close SDFILE;

  print "\nTotal Number of compounds: $FilteredSDFileInfo{CmpdCount}\n";
  print "Number of compounds left after filtering: $FilteredSDFileInfo{FilteredCmpdCount}\n";
  print "Number of compounds ignored: $FilteredSDFileInfo{KeepCmpdCount}\n";
}

# Write out the compound data...
sub WriteOutCmpdString {
  my($CmpdString) = @_;

  if ($FilteredSDFileInfo{FilterCmpd}) {
    $FilteredSDFileInfo{KeepCmpdCount} += 1;
    if ($OptionsInfo{Keep}) {
      print NEWKEEPSDFILE "$CmpdString\n";
    }
  }
  else {
    $FilteredSDFileInfo{FilteredCmpdCount} += 1;
    print NEWSDFILE "$CmpdString\n";
  }
}

# Retrieve information about input SD files...
sub RetrieveSDFilesInfo {
  my($Index, $SDFile, $FileDir, $FileName, $FileExt, $NewSDFile, $NewKeepSDFile);

  %SDFilesInfo = ();
  @{$SDFilesInfo{FileOkay}} = ();
  @{$SDFilesInfo{OutFile}} = ();
  @{$SDFilesInfo{OutFileKeep}} = ();

   FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFile = $SDFilesList[$Index];

    $SDFilesInfo{FileOkay}[$Index] = 0;
    $SDFilesInfo{OutFile}[$Index] = '';
    $SDFilesInfo{OutFileKeep}[$Index] = '';

    if (!(-e $SDFile)) {
      warn "Warning: Ignoring file $SDFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($SDFile, "sd sdf")) {
      warn "Warning: Ignoring file $SDFile: It's not a SD file\n";
      next FILELIST;
    }

    # Setup new file names...
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($SDFile);
    if ($Options{root} && (@SDFilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($Options{root});
      if ($RootFileName && $RootFileExt) {
	$NewSDFile = $RootFileName;
      }
      else {
	$NewSDFile = $Options{root};
      }
      $NewKeepSDFile = $NewSDFile;
    }
    else {
      $NewSDFile = $FileName . "Filtered";
      $NewKeepSDFile = $FileName;
    }
    $NewSDFile .= ".$FileExt";
    $NewKeepSDFile .= "Ignored" . ".$FileExt";
    if (!$Options{overwrite}) {
      if (-e $NewSDFile) {
	warn "Warning: Ignoring file $SDFile: New SD file, $NewSDFile, already exists\n";
	next FILELIST;
      }
      if ($Options{keep}) {
	if (-e $NewKeepSDFile) {
	  warn "Warning: Ignoring file $SDFile: New SD file, $NewKeepSDFile, already exists\n";
	  next FILELIST;
	}
      }
    }
    if (lc($NewSDFile) eq lc($SDFile)) {
      warn "Warning: Ignoring file $SDFile: Same output, $NewSDFile, and input file name\n";
      print "Specify a different name using \"-r --root\" option or use default name.\n";
      next FILELIST;
    }

    $SDFilesInfo{FileOkay}[$Index] = 1;
    $SDFilesInfo{OutFile}[$Index] = $NewSDFile;
    $SDFilesInfo{OutFileKeep}[$Index] = $NewKeepSDFile;
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{All} = $Options{all} ? $Options{all} : undef;
  $OptionsInfo{CleanSalts} = $Options{cleansalts} ? $Options{cleansalts} : undef;
  $OptionsInfo{Empty} = $Options{empty} ? $Options{empty} : undef;
  $OptionsInfo{Keep} = $Options{keep} ? $Options{keep} : undef;
  $OptionsInfo{Mismatch} = $Options{mismatch} ? $Options{mismatch} : undef;
  $OptionsInfo{Overwrite} = $Options{overwrite} ? $Options{overwrite} : undef;
  $OptionsInfo{Salts} = $Options{salts} ? $Options{salts} : undef;
  $OptionsInfo{UnknownAtoms} = $Options{unknownatoms} ? $Options{unknownatoms} : undef;

}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  if (!GetOptions(\%Options, "all|a", "cleansalts|c", "empty|e", "help|h", "keep|k", "mismatch|m", "overwrite|o", "root|r=s", "salts|s", "unknownatoms|u", "workingdir|w=s")) {
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

FilterSDFiles.pl - Filter compounds from SDFile(s)

=head1 SYNOPSIS

FilterSDFiles.pl SDFile(s)...

FilterSDFiles.pl [B<-a, --all>] [B<-e, --empty>] [B<-c, --cleansalts>] [B<-h, --help>]
[B<-k, --keep>] [B<-m, --mismatch>] [B<-o, --overwrite>] [B<-r, --root> I<rootname>]
[B<-s, --salts>] [B<-u, --unknownatoms>] [B<-w, --workingdir> I<dirname>] SDFile(s)...

=head1 DESCRIPTION

Filter specific compounds from I<SDFile(s)>. Available choices are: wash or
remove compounds with salts; take out compounds with no
structural data; remove compounds with mismatched atom/bond blocks data;
remove compounds which contain uknown atoms and so on. Multiple SDFile
names are separated by spaces. The valid file extensions are I<.sdf> and I<.sd>.
All other file names are ignored. All the SD files in a current directory can be
specified either by I<*.sdf> or the current directory name.

=head1 OPTIONS

=over 4

=item B<-a, --all>

Use all options to filter compounds.

=item B<-e, --empty>

Filter compounds with empty atom/bond blocks. This is B<default behavior>.

=item B<-c, --cleansalts>

Wash compounds which contain salts identified as disconnected structural
units. The largest fragment is kept.

=item B<-h, --help>

Print this help message.

=item B<-k, --keep>

Keep the compounds which were filtered in a separate file. Default: Just
ignore these compounds. Option B<-r --root> is used to generate the new file
name: <Root>Ignored.sdf. Default file name: <SDFileName>Ignored.sdf.

=item B<-m, --mismatch>

Remove compounds with mismatched atom/bond blocks and counts line
information specified by header block.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

New SD file name is generated using the root: <Root>.sdf. Default file
name:<SDFileName>Filtered.sdf. This option is ignored for multiple input files.

=item B<-s, --salts>

Remove compounds which contain salts identified as disconnected structural
units.

=item B<-u, --unknownatoms>

Remove compounds with atom blocks containing  special atom symbols such
as L, Q, * ,LP, X, R#, or any other non periodic table symbols.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To remove compounds from SD files which contain salts, unknown atoms, or
mismatched atom/bonds block data or no structural data, type:

    % FilterSDFiles.pl -a -o Sample.sdf
    % FilterSDFiles.pl -a -o *.sdf

And to generate a new NewSampleIgnored.sdf file for filtered compounds, type:

    % FilterSDFiles.pl -a -k -r NewSample -o Sample.sdf

To wash compounds in order to get rid of all disconnected fragments except for
the largest one, type:

    % FilterSDFiles.pl -c -o Sample.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

ExtractFromSDFiles.pl, InfoSDFiles.pl, MergeTextFilesWithSD.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
