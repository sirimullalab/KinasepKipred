#!/usr/bin/perl -w
#
# File: ElementalAnalysisSDFiles.pl
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
use MolecularFormula;
use FileIO::SDFileIO;
use Molecule;

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
    PerformElementalAnalysis($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Perform analysis...
sub PerformElementalAnalysis {
  my($Index) = @_;
  my($SDFile, $NewSDFile, $KeyDataFieldName, $CmpdCount, $CurrentFormula, $FormulaFieldName, $CmpdString, $Value, $CalculationType, $CalculatedValue, $ErrorMsg, $Status, $ElementsRef, $ElementCompositionRef, $Molecule, @CalculatedValues, @CmpdLines, %DataFieldValuesMap);

  $SDFile = $SDFilesList[$Index];
  $NewSDFile = $SDFilesInfo{OutFile}[$Index];

  print "Generating new SD file $NewSDFile...\n";
  open NEWSDFILE, ">$NewSDFile" or die "Error: Couldn't open $NewSDFile: $! \n";
  open SDFILE, "$SDFile" or die "Error: Can't open $SDFile: $! \n";


  $CmpdCount = 0;
  $FormulaFieldName = $SDFilesInfo{FormulaFieldName}[$Index];

  COMPOUND: while ($CmpdString = ReadCmpdString(\*SDFILE)) {
    $CmpdCount++;
    @CmpdLines = split "\n", $CmpdString;
    %DataFieldValuesMap = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

    @CalculatedValues = ();
    for $Value (@{$OptionsInfo{SpecifiedCalculations}}) {
      push @CalculatedValues, '';
    }

    $CurrentFormula = undef;

    if ($OptionsInfo{UseStructureData}) {
      $Molecule = FileIO::SDFileIO::ParseMoleculeString($CmpdString);
      $CurrentFormula = $Molecule->GetMolecularFormula();
    }
    else {
      if (!exists $DataFieldValuesMap{$FormulaFieldName}) {
	$ErrorMsg = "Ignoring compound record $CmpdCount: Formula field $FormulaFieldName not found";
	PrintErrorMsg($CmpdString, $ErrorMsg);
	WriteNewCompoundRecord($Index, \*NEWSDFILE, \@CmpdLines, \@CalculatedValues);
	next COMPOUND;
      }

      # Make sure it's a valid molecular formula...
      $CurrentFormula = $DataFieldValuesMap{$FormulaFieldName};
      if ($OptionsInfo{CheckFormula}) {
	($Status, $ErrorMsg) = MolecularFormula::IsMolecularFormula($CurrentFormula);
	if (!$Status) {
	  $ErrorMsg = "Ignoring compound record $CmpdCount: Formula field value $CurrentFormula is not valid: $ErrorMsg";
	  PrintErrorMsg($CmpdString, $ErrorMsg);
	  WriteNewCompoundRecord($Index, \*NEWSDFILE, \@CmpdLines, \@CalculatedValues);
	  next COMPOUND;
	}
      }
    }

    # Calculate appropriate values and write 'em out...
    @CalculatedValues = ();
    for $CalculationType (@{$OptionsInfo{SpecifiedCalculations}}) {
      if ($CalculationType =~ /^ElementalAnalysis$/i) {
	($ElementsRef, $ElementCompositionRef) = MolecularFormula::CalculateElementalComposition($CurrentFormula);
	$CalculatedValue = (defined($ElementsRef) && defined($ElementCompositionRef)) ? MolecularFormula::FormatCompositionInfomation($ElementsRef, $ElementCompositionRef, $OptionsInfo{Precision}) : '';
      }
      elsif ($CalculationType =~ /^MolecularWeight$/i) {
	$CalculatedValue = MolecularFormula::CalculateMolecularWeight($CurrentFormula);
	$CalculatedValue = (defined($CalculatedValue) && length($CalculatedValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $CalculatedValue)) : "";
      }
      elsif ($CalculationType =~ /^ExactMass$/i) {
	$CalculatedValue = MolecularFormula::CalculateExactMass($CurrentFormula);
	$CalculatedValue = (defined($CalculatedValue) && length($CalculatedValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $CalculatedValue)) : "";
      }
      else {
	$CalculatedValue = '';
      }
      push @CalculatedValues, $CalculatedValue;
    }
    WriteNewCompoundRecord($Index, \*NEWSDFILE, \@CmpdLines, \@CalculatedValues, $CurrentFormula);
  }
  close NEWSDFILE;
  close SDFILE;
}

# Write out compound record with calculated values...
sub WriteNewCompoundRecord {
  my($Index, $SDFileRef, $CmpdLinesRef, $CalculatedValuesRef, $MolecularFormula) = @_;

  # Write out compound lines except the last line which contains $$$$...
  my($LineIndex);
  for $LineIndex (0 .. ($#{$CmpdLinesRef} - 1)) {
    print $SDFileRef "$CmpdLinesRef->[$LineIndex]\n";
  }

  # Write out calculated values...
  my($CalcIndex, $FieldName, $FieldValue);
  for $CalcIndex (0 .. $#{$OptionsInfo{SpecifiedCalculations}}) {
    $FieldName = $SDFilesInfo{ValueFieldNamesMap}[$Index]{$OptionsInfo{SpecifiedCalculations}[$CalcIndex]};
    $FieldValue = $CalculatedValuesRef->[$CalcIndex];
    print  $SDFileRef ">  <$FieldName>\n$FieldValue\n\n";
  }

  if ($OptionsInfo{UseStructureData} && $OptionsInfo{WriteOutFormula} && defined($MolecularFormula)) {
    $FieldName = $SDFilesInfo{ValueFieldNamesMap}[$Index]{MolecularFormula};
    $FieldValue = $MolecularFormula;
    print  $SDFileRef ">  <$FieldName>\n$FieldValue\n\n";
  }

  print $SDFileRef  "\$\$\$\$\n";
}

# Print out error message...
sub PrintErrorMsg {
  my($CmpdString, $ErrorMsg) = @_;

  if ($OptionsInfo{DetailLevel} >= 2 ) {
    print "$ErrorMsg:\n$CmpdString\n";
  }
  elsif ($OptionsInfo{DetailLevel} >= 1) {
    print "$ErrorMsg\n";
  }
}

# Retrieve information about input SD files...
sub RetrieveSDFilesInfo {
  my($Index, $SDFile, $FileDir, $FileName, $FileExt, $OutFileRoot,  $OutFile, $FormulaFieldName, $Value, $FieldName, $NewFieldName, $Count);

  my(%NewValueFieldNames) = (ElementalAnalysis => 'ElementalAnalysis', MolecularWeight => 'MolecularWeight', ExactMass => 'ExactMass', MolecularFormula => 'MolecularFormula');
  if (@{$OptionsInfo{SpecifiedValueFieldNames}}) {
    for ($Index = 0; $Index < @{$OptionsInfo{SpecifiedValueFieldNames}}; $Index +=2) {
      $Value = $OptionsInfo{SpecifiedValueFieldNames}[$Index];
      $FieldName = $OptionsInfo{SpecifiedValueFieldNames}[$Index + 1];
      if (exists $NewValueFieldNames{$Value}) {
	$NewValueFieldNames{$Value} = $FieldName;
      }
    }
  }

  %SDFilesInfo = ();

  @{$SDFilesInfo{FileOkay}} = ();
  @{$SDFilesInfo{OutFile}} = ();
  @{$SDFilesInfo{FormulaFieldName}} = ();
  @{$SDFilesInfo{ValueFieldNamesMap}} = ();

  FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFile = $SDFilesList[$Index];

    $SDFilesInfo{FileOkay}[$Index] = 0;
    $SDFilesInfo{OutFile}[$Index] = '';
    $SDFilesInfo{FormulaFieldName}[$Index] = '';

    %{$SDFilesInfo{ValueFieldNamesMap}[$Index]} = ();

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
      $OutFileRoot = $FileName . "ElementalAnalysis";
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
    # Get data field names and values...
    my($CmpdString, $FieldName, @CmpdLines, @DataFieldNames, %DataFieldNamesMap);
    @DataFieldNames = ();
    if (!open(SDFILE, "$SDFile")) {
      warn "Warning: Ignoring file $SDFile: Couldn't open it: $! \n";
      next FILELIST;
    }
    $CmpdString = ReadCmpdString(\*SDFILE);
    close SDFILE;

    @CmpdLines = split "\n", $CmpdString;
    @DataFieldNames = GetCmpdDataHeaderLabels(\@CmpdLines);
    %DataFieldNamesMap = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

    # Setup formula field name...
    $FormulaFieldName = '';
    if ($OptionsInfo{UseDataField}) {
      if ($OptionsInfo{SpecifiedFormulaFieldName}) {
	$FormulaFieldName = $OptionsInfo{SpecifiedFormulaFieldName};
      }
      else {
      FIELDNAME: for $FieldName (@DataFieldNames) {
	  if ($FieldName =~ /Formula/i) {
	    $FormulaFieldName = $FieldName;
	    last FIELDNAME;
	  }
	}
	if (!$FormulaFieldName) {
	  warn "Warning: Ignoring file $SDFile: Data field label containing the word Formula doesn't exist\n";
	  next FILELIST;
	}
      }
    }
    $SDFilesInfo{FileOkay}[$Index] = 1;
    $SDFilesInfo{OutFile}[$Index] = $OutFile;
    $SDFilesInfo{FormulaFieldName}[$Index] = $FormulaFieldName;

    # Setup value data field names for calculated values...
    for $Value (keys %NewValueFieldNames) {
      $FieldName = $NewValueFieldNames{$Value};

      # Make sure it doesn't already exists...
      $Count = 1;
      $NewFieldName = $FieldName;
      while (exists $DataFieldNamesMap{$NewFieldName}) {
	$Count++;
	$NewFieldName = $FieldName . $Count;
      }
      $SDFilesInfo{ValueFieldNamesMap}[$Index]{$Value} = $NewFieldName;
    }
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};
  $OptionsInfo{FormulaMode} = $Options{formulamode};

  $OptionsInfo{WriteOutFormula} = ($Options{formulaout} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{UseStructureData} = ($Options{formulamode} =~ /^StructureData$/i) ? 1 : 0;
  $OptionsInfo{UseDataField} = ($Options{formulamode} =~ /^DataField$/i) ? 1 : 0;

  $OptionsInfo{Fast} = defined $Options{fast} ? $Options{fast} : undef;

  $OptionsInfo{DetailLevel} = $Options{detail};
  $OptionsInfo{CheckFormula} = $Options{fast} ? 0 : 1;
  $OptionsInfo{Precision} = $Options{precision};

  $OptionsInfo{Overwrite} = defined $Options{overwrite} ? $Options{overwrite} : undef;

  $OptionsInfo{FormulaField} = defined $Options{formulafield} ? $Options{formulafield} : undef;
  $OptionsInfo{SpecifiedFormulaFieldName} = "";

  if (defined $Options{formulafield}) {
    $OptionsInfo{SpecifiedFormulaFieldName} = $Options{formulafield};
  }
  # Setup what to calculate...
  @{$OptionsInfo{SpecifiedCalculations}} = ();
  if ($Options{mode} =~ /^All$/i) {
    @{$OptionsInfo{SpecifiedCalculations}} = qw(ElementalAnalysis MolecularWeight ExactMass);
  }
  else {
    my($Mode, $ModeValue, @SpecifiedModeValues);
    $Mode = $Options{mode};
    $Mode =~ s/ //g;
    @SpecifiedModeValues = split /\,/, $Mode;
    for $ModeValue (@SpecifiedModeValues) {
      if ($ModeValue !~ /^(ElementalAnalysis|MolecularWeight|ExactMass)$/i) {
	if ($ModeValue =~ /^All$/i) {
	  die "Error: All value for option \"-m --mode\" is not allowed with other valid values.\n";
	}
	else {
	  die "Error: The value specified, $ModeValue, for option \"-m --mode\" is not valid. Allowed values: ElementalAnalysis, MolecularWeight, or ExactMass\n";
	}
      }
      push @{$OptionsInfo{SpecifiedCalculations}}, $ModeValue;
    }
  }

  $OptionsInfo{ValueFieldNames} = defined $Options{valuefieldnames} ? $Options{valuefieldnames} : undef;
  @{$OptionsInfo{SpecifiedValueFieldNames}} = ();

  if ($Options{valuefieldnames}) {
    my($Value, $Label, @ValueLabels);
    @ValueLabels = split /\,/, $Options{valuefieldnames};
    if (@ValueLabels % 2) {
      die "Error: The value specified, $Options{valuefieldnames}, for option \"-v --valuefieldnames\" is not valid: It must contain even number of comma delimited values\n";
    }
    my($Index);
    for ($Index = 0; $Index < @ValueLabels; $Index +=2) {
      $Value = $ValueLabels[$Index];
      $Value =~ s/ //g;
      $Label = $ValueLabels[$Index + 1];
      if ($Value !~ /^(ElementalAnalysis|MolecularWeight|ExactMass|MolecularFormula)$/i) {
	die "Error: The value specified, $Value, using option \"-v --valuefieldnames\" is not valid. Allowed values: ElementalAnalysis, MolecularWeight, ExactMass, or MolecularFormula\n";
      }
      push @{$OptionsInfo{SpecifiedValueFieldNames}}, ($Value, $Label);
    }
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{detail} = 1;
  $Options{formulamode} = "DataField";
  $Options{formulaout} = "No";
  $Options{mode} = "All";
  $Options{precision} = 2;

  if (!GetOptions(\%Options, "detail|d=i", "fast", "formulafield=s", "formulamode|f=s", "formulaout=s", "mode|m=s", "help|h", "overwrite|o", "precision|p=i", "root|r=s", "valuefieldnames|v=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if (!IsPositiveInteger($Options{detail})) {
    die "Error: The value specified, $Options{detail}, for option \"-d --detail\" is not valid. Allowed values: > 0\n";
  }
  if ($Options{formulamode} !~ /^(StructureData|DataField)$/i) {
    die "Error: The value specified, $Options{formulamode}, for option \"-f, --formulamode\" is not valid. Allowed values: StructureData or DataField \n";
  }
  if ($Options{formulaout} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{formulaout}, for option \"--formulaout\" is not valid. Allowed values: Yes or No \n";
  }
  if (!IsPositiveInteger($Options{precision})) {
    die "Error: The value specified, $Options{precision}, for option \"-p --precision\" is not valid. Allowed values: > 0 \n";
  }
}

__END__

=head1 NAME

ElementalAnalysisSDFiles.pl - Perform elemental analysis using formula data field in SDFile(s)

=head1 SYNOPSIS

ElementalAnalysisSDFiles.pl SDFile(s)...

ElementalAnalysisSDFiles.pl [B<-d, --detail> infolevel] [B<--fast>]
[B<--formulafield> SD data field name] [B<-f, --formulamode> I<DataField | StructureData>]
[B<--formulaout> yes or no] [B<-m, --mode> All | "ElementalAnalysis, [MolecularWeight, ExactMass]"]
[B<-h, --help>] [B<-o, --overwrite>] [B<-r, --root> rootname]
[B<-v --valuefieldnames> Name, Label, [Name, Label,...]] [B<-w, --workingdir> dirname] SDFile(s)...

=head1 DESCRIPTION

Perform elemental analysis using molecular formula specified by a data field name or generated
from structure data in I<SDFile(s)>.

In addition to straightforward molecular formulas - H2O, HCl, C3H7O2N -
other supported variations are: Ca3(PO4)2, [PCl4]+, [Fe(CN)6]4-, C37H42N2O6+2, Na2CO3.10H2O,
8H2S.46H2O, and so on. Charges are simply ignored. Isotope symbols in formulas specification, including
D and T, are not supported.

The file names are separated by space.The valid file extensions are I<.sdf> and I<.sd>.
All other file names are ignored. All the SD files in a current directory can be specified
either by I<*.sdf> or the current directory name.

=head1 OPTIONS

=over 4

=item B<-d, --detail> I<infolevel>

Level of information to print about compound records being ignored. Default: I<1>. Possible
values: I<1, 2 or 3>.

=item B<--fast>

In this mode, the formula data field specified using B<-f, --formulafield> option is assumed
to contain valid molecular formula data and initial formula validation check is skipped.

=item B<--formulafield> I<SD data field name>

I<SDFile(s)> data field name containing molecular formulas used for performing
elemental analysis during I<DataField> value of B<-f, --formulamode> option.
Default value: I<SD data field containing the word formula in its name>.

This option is ignore during I<StructureData> value of B<-f, --formulamode> option.

=item B<-f, --formulamode> I<DataField | StructureData>

Specify source of molecular formula used for performing elemental analysis: retrieve
formula using I<SDFile(s)> data field name or generate formula from structure. Possible
values: I<DataField or StructureData>. Default value: I<DataField>.

=item B<--formulaout> I<yes or no>

Specify whether to write out formula to SD file during I<StructureData> value of
B<-f, --formulamode> option. Possible values: I<Yes or No>. Default: I<No>.

=item B<-m, --mode> I<All | "ElementalAnalysis,[MolecularWeight,ExactMass]">

Specify what values to calculate using molecular formula data field or structure data from
I<SDFile(s)>: calculate all supported values or specify a comma delimited list of values.
Possible values: I<All | "ElementalAnalysis, [MolecularWeight, ExactMass]">. Default: I<All>

=item B<-h, --help>

Print this help message.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-p, --precision> I<number>

Precision of calculated values in the output file. Default: up to I<2> decimal places.
Valid values: positive integers.

=item B<-r, --root> I<rootname>

New SD file name is generated using the root: <Root>.<Ext>. Default new file
name: <InitialSDFileName>ElementalAnalysis.<Ext>. This option is ignored for multiple
input files.

=item B<-v --valuefieldnames> I<Name,Label,[Name,Label,...]>

Specify SD data field names to use for calculated values. In general, it's a comma delimited
list of value name and SD field name  pairs. Supported value names: I<ElementalAnalysis,
MolecularWeight, ExactMass, and MolecularFormula>. Default labels: I<ElementalAnalysis,
MolecularWeight, ExactMass, and MolecularFormula>.

I<MolecularFormula> label is only used during I<StructureData> value of
B<-f, --formulamode> option.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To perform elemental analysis, calculate molecular weight and exact mass using SD
field name value with the word Formula in its name and generate a new SD file
NewSample1.sdf, type:

    % ElementalAnalysisSDFiles.pl -o -r NewSample1 Sample1.sdf

To perform elemental analysis, calculate molecular weight and exact mass using
structure data in SD file and generate a new SD file NewSample1.sdf, type:

    % ElementalAnalysisSDFiles.pl --formulamode StructureData -o
      -r NewSample1 Sample1.sdf

To perform elemental analysis using formulas in SD field name Formula, use field name
Analysis for calculated data, and generate a new SD file NewSample1.sdf, type:

    % ElementalAnalysisSDFiles.pl --m ElementalAnalysis --formulafield
      Formula --valuefieldnames "ElementalAnalysis,Analysis" -o
      -r NewSample1 Sample1.sdf

To calculate molecular weight, using formulas in SD field name Formula, with four decimal
precision and generate a new SD file NewSample1.sdf, type

    % ElementalAnalysisSDFiles.pl --m MolecularWeight --formulafield
      Formula --precision 4 -o -r NewSample1 Sample1.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AnalyzeSDFilesData.pl, InfoSDFiles.pl, ExtractFromSDFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
