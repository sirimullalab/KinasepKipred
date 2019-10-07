#!/usr/bin/perl -w
#
# File: InfoSDFiles.pl
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
use TextUtil;
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
my(%SDFilesInfo, %SDCmpdsInfo);
RetrieveSDFilesInfo();
InitializeSDCmpdsInfo();

# Process input files..
my($FileIndex);
if (@SDFilesList > 1) {
  print "\nProcessing SD files...\n";
}
for $FileIndex (0 .. $#SDFilesList) {
  if ($SDFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $SDFilesList[$FileIndex]...\n";
    ListSDFileInfo($FileIndex);
  }
}
ListTotalSizeOfFiles();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# List appropriate information...
sub ListSDFileInfo {
  my($Index) = @_;
  my($SDFile);

  $SDFile = $SDFilesList[$Index];

  if ($OptionsInfo{ProcessCmpdInfo}) {
    ListCompoundDetailsInfo($Index);
  }
  else {
    ListCompoundCountInfo($Index);
  }

  # File size and modification information...
  print "\nFile size: ", FormatFileSize($SDFilesInfo{FileSize}[$Index]), " \n";
  print "Last modified: ", $SDFilesInfo{FileLastModified}[$Index], " \n";
}

# List number of compounds in SD file...
sub ListCompoundCountInfo {
  my($Index) = @_;
  my($SDFile, $CmpdCount);

  $SDFile = $SDFilesList[$Index];

  $CmpdCount = 0;

  open SDFILE, "$SDFile" or die "Couldn't open $SDFile: $! \n";
  while (<SDFILE>) {
    if (/^\$\$\$\$/) {
      $CmpdCount++;
    }
  }
  close SDFILE;

  $SDCmpdsInfo{TotalCmpdCount} += $CmpdCount;

  print "\nNumber of compounds: $CmpdCount\n";
}

# List detailed compound information...
sub ListCompoundDetailsInfo {
  my($Index) = @_;
  my($SDFile, $CmpdCount, $EmptyCtabBlocksCount, $MismatchCtabBlockCount, $ChiralCtabBlockCount, $UnknownAtomsCtabBlockCount, $InvalidAtomNumbersCtabBlockCount, $SaltsCtabBlockCount, $CtabLinesCount, $PrintCmpdCounterHeader, $ProblematicCmpdData, $CmpdString, @CmpdLines);

  $SDFile = $SDFilesList[$Index];

  ($CmpdCount, $EmptyCtabBlocksCount, $MismatchCtabBlockCount, $ChiralCtabBlockCount, $UnknownAtomsCtabBlockCount, $InvalidAtomNumbersCtabBlockCount, $SaltsCtabBlockCount) = (0) x 7;

  InitializeSDCmpdsInfo();

  $PrintCmpdCounterHeader = 1;

  open SDFILE, "$SDFile" or die "Couldn't open $SDFile: $! \n";
  while ($CmpdString = ReadCmpdString(\*SDFILE)) {
    $CmpdCount++;
    $ProblematicCmpdData = 0;
    if ($OptionsInfo{Detail} <= 1) {
      if (($CmpdCount % 5000) == 0) {
	if ($PrintCmpdCounterHeader) {
	  $PrintCmpdCounterHeader = 0;
	  print "Processing compounds:";
	}
	print "$CmpdCount...";
      }
    }
    @CmpdLines = split "\n", $CmpdString;
    $CtabLinesCount = GetCtabLinesCount(\@CmpdLines);
    if ($OptionsInfo{All} || $OptionsInfo{Empty}) {
      if ($CtabLinesCount <= 0) {
	$EmptyCtabBlocksCount++;
	$ProblematicCmpdData = 1;
      }
    }
    if ($CtabLinesCount > 0) {
      my ($AtomCount, $BondCount, $ChiralFlag) = ParseCmpdCountsLine($CmpdLines[3]);
      if ($OptionsInfo{All} || $OptionsInfo{Mismatch}) {
	if ($CtabLinesCount != ($AtomCount + $BondCount)) {
	  $MismatchCtabBlockCount++;
	  $ProblematicCmpdData = 1;
	  if ($OptionsInfo{Detail} >= 2) {
	    print "\nMismatch found: Ctab lines count: $CtabLinesCount;  Atoms count: $AtomCount; Bond count: $BondCount\n";
	  }
	}
      }
      if ($OptionsInfo{All} || $OptionsInfo{Chiral}) {
	if ($ChiralFlag == 1) {
	  $ChiralCtabBlockCount++;
	}
      }
      if ($CtabLinesCount == ($AtomCount + $BondCount)) {
	if ($OptionsInfo{All} || $OptionsInfo{UnknownAtoms}) {
	  my($UnknownAtomCount, $UnknownAtoms, $UnknownAtomLines) = GetUnknownAtoms(\@CmpdLines);
	  if ($UnknownAtomCount) {
	    $UnknownAtomsCtabBlockCount++;
	    $ProblematicCmpdData = 1;
	    if ($OptionsInfo{Detail} >= 2) {
	      print "\nUnknown atom(s) found: $UnknownAtomCount\nUnknown atom(s) symbols:$UnknownAtoms\nUnknown atom(s) data lines:\n$UnknownAtomLines\n";
	    }
	  }
	}
	if ($OptionsInfo{All} || $OptionsInfo{InvalidAtomNumbers}) {
	  my($InvalidAtomNumbersCount, $InvalidAtomNumbers, $InvalidAtomNumberLines) = GetInvalidAtomNumbers(\@CmpdLines);
	  if ($InvalidAtomNumbersCount) {
	    $InvalidAtomNumbersCtabBlockCount++;
	    $ProblematicCmpdData = 1;
	    if ($OptionsInfo{Detail} >= 2) {
	      print "\nInvalid atom number(s) found: $InvalidAtomNumbersCount\nInvalid atom number(s):$InvalidAtomNumbers\nInvalid atom number(s) data lines:\n$InvalidAtomNumberLines\n";
	    }
	  }
	}
	if ($OptionsInfo{All} || $OptionsInfo{Salts}) {
	  my($FragmentsCount, $Fragments) = GetCmpdFragments(\@CmpdLines);
	  if ($FragmentsCount > 1) {
	    $SaltsCtabBlockCount++;
	    $ProblematicCmpdData = 1;
	    if ($OptionsInfo{Detail} >= 2) {
	      print "\nSalts found: $FragmentsCount\nSalts atom numbers:\n$Fragments\n";
	    }
	  }
	}
      }
    }
    if ($OptionsInfo{ProcessCmpdData}) {
      ProcessCmpdInfo(\@CmpdLines, $CmpdCount);
    }
    if ($OptionsInfo{Detail} >= 3) {
      if ($ProblematicCmpdData) {
	print "\nCompound data:\n$CmpdString\n\n";
      }
    }
  }
  if ($OptionsInfo{Detail} <= 1) {
    if (!$PrintCmpdCounterHeader) {
      print "\n";
    }
  }
  close SDFILE;

  $SDCmpdsInfo{TotalCmpdCount} += $CmpdCount;

  print "\nNumber of compounds: $CmpdCount\n";

  if ($OptionsInfo{All} || $OptionsInfo{Empty}) {
    print "Number of empty atom/bond blocks: $EmptyCtabBlocksCount\n";
  }
  if ($OptionsInfo{All} || $OptionsInfo{Mismatch}) {
    print "Number of mismatched atom/bond blocks: $MismatchCtabBlockCount\n";
  }
  if ($OptionsInfo{All} || $OptionsInfo{UnknownAtoms}) {
    print "Number of atom blocks with unknown atom labels: $UnknownAtomsCtabBlockCount\n";
  }
  if ($OptionsInfo{All} || $OptionsInfo{InvalidAtomNumbers}) {
    print "Number of bond blocks and atom property blocks with invalid atom numbers: $InvalidAtomNumbersCtabBlockCount\n";
  }
  if ($OptionsInfo{All} || $OptionsInfo{Salts}) {
    print "Number of atom blocks containing salts: $SaltsCtabBlockCount\n";
  }
  if ($OptionsInfo{All} || $OptionsInfo{Chiral}) {
    print "Number of chiral atom/bond blocks: $ChiralCtabBlockCount\n";
  }
  if ($OptionsInfo{ProcessCmpdData}) {
    PrintCmpdInfoSummary();
  }

}

# Initialize compound data information for a SD file...
sub InitializeSDCmpdsInfo {

  if (!exists $SDCmpdsInfo{TotalCmpdCount}) {
    $SDCmpdsInfo{TotalCmpdCount} = 0;
  }

  @{$SDCmpdsInfo{FieldLabels}} = ();
  %{$SDCmpdsInfo{FieldLabelsMap}} = ();
  %{$SDCmpdsInfo{NonEmptyFieldValuesCountMap}} = ();
  %{$SDCmpdsInfo{EmptyFieldValuesCountMap}} = ();
  %{$SDCmpdsInfo{NonNumericalFieldValuesCountMap}} = ();
  %{$SDCmpdsInfo{NumericalFieldValuesCountMap}} = ();
}

# Process compound data header labels and figure out which ones are present for
# all the compounds...
sub ProcessCmpdInfo {
  my($CmpdLinesRef, $CmpdCount) = @_;
  my($Label);

  if (@{$SDCmpdsInfo{FieldLabels}}) {
    my (@CmpdFieldLabels) = GetCmpdDataHeaderLabels($CmpdLinesRef);
    my(%CmpdFieldLabelsMap) = ();
    # Setup a map for the current labels...
    for $Label (@CmpdFieldLabels) {
      $CmpdFieldLabelsMap{$Label} = "PresentInSome";
    }
    # Check the presence old labels for this compound; otherwise, mark 'em new...
    for $Label (@{$SDCmpdsInfo{FieldLabels}}) {
      if (!$CmpdFieldLabelsMap{$Label}) {
	$SDCmpdsInfo{FieldLabelsMap}{$Label} = "PresentInSome";
      }
    }
    # Check the presence this compound in the old labels; otherwise, add 'em...
    for $Label (@CmpdFieldLabels ) {
      if (!$SDCmpdsInfo{FieldLabelsMap}{$Label}) {
	# It's a new label...
	push @{$SDCmpdsInfo{FieldLabels}}, $Label;
	$SDCmpdsInfo{FieldLabelsMap}{$Label} = "PresentInSome";
      }
    }
  }
  else {
    # Get the initial label set and set up a map...
    @{$SDCmpdsInfo{FieldLabels}} = GetCmpdDataHeaderLabels($CmpdLinesRef);
    for $Label (@{$SDCmpdsInfo{FieldLabels}}) {
      $SDCmpdsInfo{FieldLabelsMap}{$Label} = "PresentInAll";
    }
  }
  if ($OptionsInfo{CountEmptyData} || $OptionsInfo{CheckData}) {
    # Count empty data field values...
    my(%DataFieldAndValues, $Label, $Value);

    %DataFieldAndValues = GetCmpdDataHeaderLabelsAndValues($CmpdLinesRef);
    for $Label (keys %DataFieldAndValues) {
      $Value = $DataFieldAndValues{$Label};
      if ($OptionsInfo{CountEmptyData}) {
	if (IsNotEmpty($Value)) {
	  if (exists($SDCmpdsInfo{NonEmptyFieldValuesCountMap}{$Label})) {
	    $SDCmpdsInfo{NonEmptyFieldValuesCountMap}{$Label} += 1;
	  }
	  else {
	    $SDCmpdsInfo{NonEmptyFieldValuesCountMap}{$Label} = 1;
	  }
	}
	else {
	  if ($Options{detail} >= 2) {
	    print "Compound record $CmpdCount: Empty data field <$Label>\n";
	  }
	  if (exists($SDCmpdsInfo{EmptyFieldValuesCountMap}{$Label})) {
	    $SDCmpdsInfo{EmptyFieldValuesCountMap}{$Label} += 1;
	  }
	  else {
	    $SDCmpdsInfo{EmptyFieldValuesCountMap}{$Label} = 1;
	  }
	}
      }
      if ($OptionsInfo{CheckData}) {
	if (IsNumerical($Value)) {
	  if (exists($SDCmpdsInfo{NumericalFieldValuesCountMap}{$Label})) {
	    $SDCmpdsInfo{NumericalFieldValuesCountMap}{$Label} += 1;
	  }
	  else {
	    $SDCmpdsInfo{NumericalFieldValuesCountMap}{$Label} = 1;
	  }
	}
	else {
	  if (exists($SDCmpdsInfo{NonNumericalFieldValuesCountMap}{$Label})) {
	    $SDCmpdsInfo{NonNumericalFieldValuesCountMap}{$Label} += 1;
	  }
	  else {
	    $SDCmpdsInfo{NonNumericalFieldValuesCountMap}{$Label} = 1;
	  }
	}
      }
    }
  }
}

# Print compound summary...
sub PrintCmpdInfoSummary {
  if (@{$SDCmpdsInfo{FieldLabels}}) {
    my($PresentInAllCount, $Label, @FieldLabelsPresentInSome, @FieldLabelsPresentInAll);

    @FieldLabelsPresentInSome = ();
    @FieldLabelsPresentInAll = ();

    $PresentInAllCount = 0;
    print "\nNumber of data fields: ", scalar(@{$SDCmpdsInfo{FieldLabels}}), "\n";
    print "All data field labels: ";
    for $Label (sort keys %{$SDCmpdsInfo{FieldLabelsMap}}) {
      print "<$Label> ";
    }
    print "\n";
    for $Label (sort keys %{$SDCmpdsInfo{FieldLabelsMap}}) {
      if ($SDCmpdsInfo{FieldLabelsMap}{$Label} eq "PresentInAll") {
	$PresentInAllCount++;
	push @FieldLabelsPresentInAll, $Label;
      }
    }
    if ($PresentInAllCount != @{$SDCmpdsInfo{FieldLabels}}) {
      print "Data field labels present in all compounds: ";
      for $Label (sort keys %{$SDCmpdsInfo{FieldLabelsMap}}) {
	if ($SDCmpdsInfo{FieldLabelsMap}{$Label} eq "PresentInAll") {
	  print "<$Label> ";
	}
      }
      print "\n";
      print "Data field labels present in some compounds: ";
      for $Label (sort keys %{$SDCmpdsInfo{FieldLabelsMap}}) {
	if ($SDCmpdsInfo{FieldLabelsMap}{$Label} eq "PresentInSome") {
	  print "<$Label> ";
	  push @FieldLabelsPresentInSome, $Label;
	}
      }
      print "\n";
    }
    # List empty data field values count...
    if ($OptionsInfo{CountEmptyData}) {
      print "\n";
      if ($PresentInAllCount == @{$SDCmpdsInfo{FieldLabels}}) {
	PrintDataInformation("Number of non-empty values for data field(s)", \@{$SDCmpdsInfo{FieldLabels}}, \%{$SDCmpdsInfo{NonEmptyFieldValuesCountMap}});
	PrintDataInformation("Number of empty values for data field(s)", \@{$SDCmpdsInfo{FieldLabels}}, \%{$SDCmpdsInfo{EmptyFieldValuesCountMap}});
      }
      else {
	PrintDataInformation("Number of non-empty values for data field(s) present in all compounds", \@FieldLabelsPresentInAll, \%{$SDCmpdsInfo{NonEmptyFieldValuesCountMap}});
	PrintDataInformation("Number of empty values for data field(s) present in all compounds", \@FieldLabelsPresentInAll, \%{$SDCmpdsInfo{EmptyFieldValuesCountMap}});
	PrintDataInformation("Number of non-empty values for data field(s) present in some compounds", \@FieldLabelsPresentInSome, \%{$SDCmpdsInfo{NonEmptyFieldValuesCountMap}});
	PrintDataInformation("Number of empty values for data field(s) present in some compounds", \@FieldLabelsPresentInSome, \%{$SDCmpdsInfo{EmptyFieldValuesCountMap}});
      }
      print "\n";
    }
    # List numerical data values count...
    if ($OptionsInfo{CheckData}) {
      print "\n";
      if ($PresentInAllCount == @{$SDCmpdsInfo{FieldLabels}}) {
	PrintDataInformation("Number of non-numerical values for data field(s)", \@{$SDCmpdsInfo{FieldLabels}}, \%{$SDCmpdsInfo{NonNumericalFieldValuesCountMap}});
	PrintDataInformation("Number of numerical values for data field(s)", \@{$SDCmpdsInfo{FieldLabels}}, \%{$SDCmpdsInfo{NumericalFieldValuesCountMap}});
      }
      else {
	PrintDataInformation("Number of non-numerical values for data field(s) present in all compounds", \@FieldLabelsPresentInAll, \%{$SDCmpdsInfo{NonNumericalFieldValuesCountMap}});
	PrintDataInformation("Number of numerical values for data field(s) present in all compounds", \@FieldLabelsPresentInAll, \%{$SDCmpdsInfo{NumericalFieldValuesCountMap}});
	PrintDataInformation("Number of non-numerical values for data field(s) present in some compounds", \@FieldLabelsPresentInSome, \%{$SDCmpdsInfo{NonNumericalFieldValuesCountMap}});
	PrintDataInformation("Number of numerical values for data field(s) present in some compounds", \@FieldLabelsPresentInSome, \%{$SDCmpdsInfo{NumericalFieldValuesCountMap}});
      }
      print "\n";
    }
  }
  else {
    print "\nNumber of data fields: 0\n";
  }
}
# List data information...
sub PrintDataInformation {
  my($InfoLabel, $DataLabelRef, $DataLabelToValueMapRef) = @_;
  my($Line, $Label);

  $Line = "";
  for $Label (@{$DataLabelRef}) {
    $Line .= " <$Label> - " . (exists($DataLabelToValueMapRef->{$Label}) ? $DataLabelToValueMapRef->{$Label} : 0) . ",";
  }
  $Line =~ s/\,$//g;
  print "$InfoLabel: $Line\n";
}

# Total size of all the files...
sub ListTotalSizeOfFiles {
  my($FileOkayCount, $TotalSize, $Index);

  $FileOkayCount = 0;
  $TotalSize = 0;

  for $Index (0 .. $#SDFilesList) {
    if ($SDFilesInfo{FileOkay}[$Index]) {
      $FileOkayCount++;
      $TotalSize += $SDFilesInfo{FileSize}[$Index];
    }
  }
  if ($FileOkayCount > 1) {
    print "\nTotal number of compounds in  $FileOkayCount SD files: $SDCmpdsInfo{TotalCmpdCount}\n";
    print "\nTotal size of $FileOkayCount SD files: ", FormatFileSize($TotalSize), "\n";
  }

}

# Retrieve information about SD files...
sub RetrieveSDFilesInfo {
  my($Index, $SDFile, $ModifiedTimeString, $ModifiedDateString);

  %SDCmpdsInfo = ();

  %SDFilesInfo = ();
  @{$SDFilesInfo{FileOkay}} = ();
  @{$SDFilesInfo{FileSize}} = ();
  @{$SDFilesInfo{FileLastModified}} = ();

  FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFilesInfo{FileOkay}[$Index] = 0;
    $SDFilesInfo{FileSize}[$Index] = 0;
    $SDFilesInfo{FileLastModified}[$Index] = '';

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
    $SDFilesInfo{FileSize}[$Index] = FileSize($SDFile);
    ($ModifiedTimeString, $ModifiedDateString) = FormattedFileModificationTimeAndDate($SDFile);
    $SDFilesInfo{FileLastModified}[$Index] = "$ModifiedTimeString; $ModifiedDateString";
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{All} = $Options{all} ? $Options{all} : 0;
  $OptionsInfo{Chiral} = $Options{chiral} ? $Options{chiral} : 0;
  $OptionsInfo{Count} = $Options{count} ? $Options{count} : 0;
  $OptionsInfo{DataCheck} = $Options{datacheck} ? $Options{datacheck} : 0;
  $OptionsInfo{Empty} = $Options{empty} ? $Options{empty} : 0;
  $OptionsInfo{Fields} = $Options{fields} ? $Options{fields} : 0;
  $OptionsInfo{InvalidAtomNumbers} = $Options{invalidatomnumbers} ? $Options{invalidatomnumbers} : 0;
  $OptionsInfo{Mismatch} = $Options{mismatch} ? $Options{mismatch} : 0;
  $OptionsInfo{Salts} = $Options{salts} ? $Options{salts} : 0;
  $OptionsInfo{UnknownAtoms} = $Options{unknownatoms} ? $Options{unknownatoms} : 0;

  $OptionsInfo{Detail} = $Options{detail};

  $OptionsInfo{ProcessCmpdInfo} = ($Options{all} ||  $Options{chiral} || $Options{empty} || $Options{fields} || $Options{invalidatomnumbers}  || $Options{mismatch} || $Options{salts} || $Options{unknownatoms} || $Options{datacheck}) ? 1 : 0;

  $OptionsInfo{ProcessCmpdData} = ($Options{all} || $Options{fields} || $Options{empty} || $Options{datacheck}) ? 1 : 0;

  $OptionsInfo{CountEmptyData} = ($Options{all} || $Options{empty}) ? 1 : 0;
  $OptionsInfo{CheckData} = ($Options{all} || $Options{datacheck}) ? 1 : 0;
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Setup default and retrieve all the options...
  %Options = ();
  $Options{detail} = 1;
  if (!GetOptions(\%Options, "all|a", "count|c", "chiral", "datacheck", "detail|d:i", "empty|e", "fields|f", "help|h", "invalidatomnumbers|i", "mismatch|m", "salts|s", "unknownatoms|u", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{detail} <= 0 || $Options{detail} > 3) {
    die "Error: The value specified, $Options{detail}, for option \"-d --detail\" is not valid. Possible values: 1 to 3\n";
  }
}

__END__

=head1 NAME

InfoSDFiles.pl - List information about SDFile(s)

=head1 SYNOPSIS

InfoSDFile.pl SDFile(s)...

InfoSDFile.pl [B<-a --all>] [B<-c --count>] [B<--chiral>] [B<--datacheck>]
[B<-d --detail> infolevel] [B<-e --empty>] [B<-f, --fields>] [B<-h, --help>]
[B<-i, --invalidatomnumbers>] [B<-m, --mismatch>] [B<-s, --salts>] [B<-u, --unknownatoms>]
[B<-w, --workingdir> dirname] SDFile(s)...

=head1 DESCRIPTION

List information about I<SDFile(s)> contents: number of compounds, empty records
and so on. Multiple SDFile names are separated by spaces. The valid file extensions
are I<.sdf> and I<.sd>. All other file names are ignored. All the SD files in a current
directory can be specified either by I<*.sdf> or the current directory name.

=head1 OPTIONS

=over 4

=item B<-a, --all>

List all the available information.

=item B<-c, --count>

List number of compounds. This is B<default behavior>.

=item B<--chiral>

List number of empty atom/bond blocks for compounds with chiral flag set in
count line.

=item B<-d, --detail> I<infolevel>

Level of information to print. Default: 1. Possible values: I<1, 2, or 3>.

=item B<--datacheck>

List number of numerical and non-numerical values for each data field.

=item B<-e, --empty>

List number of empty atom/bond blocks and data fields for compounds.

=item B<-f, --fields>

List data field labels present for compounds.

=item B<-h, --help>

Print this help message.

=item B<-i, --invalidatomnumbers>

List number of bond blocks for compounds which contain invalid atom numbers.

=item B<-m, --mismatch>

List number of atom/bond blocks for compounds which don't match with counts
line information in header block.

=item B<-s, --salts>

List number of atom blocks for compounds which contain salts identified as
disconnected structural units.

=item B<-u, --unknownatoms>

List number of atom blocks for compounds which contain special atom symbols
such as L, Q, * ,LP, X, R#, or any other non periodic table symbols.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To count compounds in SD file(s), type:

    % InfoSDFiles.pl Sample1.sdf
    % InfoSDFiles.pl Sample1.sdf Sample2.sdf
    % InfoSDFiles.pl *.sdf

To list all available information for SD file(s), type:

    % InfoSDFiles.pl -a *.sdf

To list all data fields present in sample.sdf, type:

    % InfoSDFiles.pl -f Sample.sdf

To count number of compounds which contain salts and list associated structural
data, type:

    % InfoSDFiles.pl -s -d 3 Sample.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

ExtractFromSDFiles.pl, FilterSDFiles.pl, MergeTextFilesWithSD.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
