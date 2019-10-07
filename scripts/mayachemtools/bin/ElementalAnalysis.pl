#!/usr/bin/perl -w
#
# File: ElementalAnalysis.pl
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
use MolecularFormula;

my($ScriptName, %Options, $StartTime, $EndTime, $TotalTime);

# Autoflush STDOUT
$| = 1;

# Starting message...
$ScriptName = basename($0);
print "\n$ScriptName: Starting...\n\n";
$StartTime = new Benchmark;

# Get the options and setup script...
SetupScriptUsage();
if ($Options{help}) {
  die GetUsageFromPod("$FindBin::Bin/$ScriptName");
}

print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

PerformElementalAnalysis();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Elemental analysis...
sub PerformElementalAnalysis {
  my($Formula, $FormulaDataRef, $ValueLabel, $CalculationType, $CalculatedValue, $ElementsRef, $ElementCompositionRef, $Status, $ErrorMsg, @ValueLabels, @CalculatedValues);

  print "Performing elemental analysis...\n";

  if ($OptionsInfo{FileOutput}) {
    print "Generating file $OptionsInfo{OutFileName}...\n";
    open OUTFILE, ">$OptionsInfo{OutFileName}" or die "Couldn't open $OptionsInfo{OutFileName}: $!\n";
  }

  # Setup value labels...
  @ValueLabels = ();
  if ($OptionsInfo{RowsOutput}) {
    push @ValueLabels, "Formula";
  }
  for $CalculationType (@{$OptionsInfo{SpecifiedCalculations}}) {
    push @ValueLabels, $OptionsInfo{ValueLabelsMap}{$CalculationType};
  }

  if ($OptionsInfo{RowsOutput}) {
    ListHeaderRowData(\@ValueLabels);
  }

  # Go over specified properties...
  FORMULA: for $Formula (@{$OptionsInfo{SpecifiedFormulas}}) {
    if (!$OptionsInfo{RowsOutput}) {
      if ($OptionsInfo{FileOutput}) {
	print OUTFILE "\nPerforming elemental analysis using formula $Formula...\n\n";
      }
      else {
	print "\nPerforming elemental analysis using formula  $Formula...\n\n";
      }
    }
    # Calculate appropriate values and write 'em out...
    if ($OptionsInfo{CheckFormula}) {
      ($Status, $ErrorMsg) = MolecularFormula::IsMolecularFormula($Formula);
      if (!$Status) {
	warn("Warning: Ignoring formula $Formula: It's not a valid value: $ErrorMsg\n");
	next FORMULA;
      }
    }
    @CalculatedValues = ();
    if ($OptionsInfo{RowsOutput}) {
      push @CalculatedValues, $Formula;
    }
    for $CalculationType (@{$OptionsInfo{SpecifiedCalculations}}) {
      if ($CalculationType =~ /^ElementalAnalysis$/i) {
	($ElementsRef, $ElementCompositionRef) = MolecularFormula::CalculateElementalComposition($Formula);
	$CalculatedValue = (defined($ElementsRef) && defined($ElementCompositionRef)) ? MolecularFormula::FormatCompositionInfomation($ElementsRef, $ElementCompositionRef, $OptionsInfo{Precision}) : '';
      }
      elsif ($CalculationType =~ /^MolecularWeight$/i) {
	$CalculatedValue = MolecularFormula::CalculateMolecularWeight($Formula);
	$CalculatedValue = (defined($CalculatedValue) && length($CalculatedValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $CalculatedValue)) : "";
      }
      elsif ($CalculationType =~ /^ExactMass$/i) {
	$CalculatedValue = MolecularFormula::CalculateExactMass($Formula);
	$CalculatedValue = (defined($CalculatedValue) && length($CalculatedValue)) ? (sprintf("%.$OptionsInfo{Precision}f", $CalculatedValue)) : "";
      }
      else {
	$CalculatedValue = '';
      }
      push @CalculatedValues, $CalculatedValue;
    }
    # List data...
    ListFormulaData(\@ValueLabels, \@CalculatedValues);
  }
  if ($OptionsInfo{FileOutput}) {
    close OUTFILE;
  }
  print "\n";
}

# List calculated data values...
sub ListFormulaData {
  my($DataLabelRef, $DataValueRef) = @_;
  my($Index, $Line, $Value);

  if ($OptionsInfo{RowsOutput}) {
    $Line = '';
    # Format data...
    if ($OptionsInfo{OutQuote} || $Options{outdelim} !~ /^comma$/i) {
      $Line = JoinWords($DataValueRef, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    }
    else {
      # Always quote values containing commas...
      $Line = ($DataValueRef->[0] =~ /\,/) ? qq("$DataValueRef->[0]") : $DataValueRef->[0];
      for $Index (1 .. $#{$DataValueRef} ) {
	$Value = $DataValueRef->[$Index];
	if ($Value =~ /\,/) {
	  $Value = qq("$Value");
	}
	$Line .= $OptionsInfo{OutDelim} . $Value;
      }
    }
    if ($OptionsInfo{FileOutput}) {
      print OUTFILE "$Line\n";
    }
    else {
      print "$Line\n";
    }
  }
  else {
    # Format and list data...
    $Line = '';
    for $Index (0 .. $#{$DataLabelRef} ) {
      $Line = $DataLabelRef->[$Index] . ': ' . $DataValueRef->[$Index];
      if ($OptionsInfo{FileOutput}) {
	print OUTFILE "$Line\n";
      }
      else {
	print "$Line\n";
      }
    }
  }
}

# List calculated data for a formula...
sub ListHeaderRowData {
  my($DataLabelRef) = @_;
  my($Line);

  # Format data...
  $Line = JoinWords($DataLabelRef, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  $Line =~ s/\://g;
  # List data...
  if ($OptionsInfo{FileOutput}) {
    print OUTFILE "$Line\n";
  }
  else {
    print "$Line\n";
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /^tab$/i ) ? "\t" : (($Options{outdelim} =~ /^semicolon$/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^yes$/i) ? 1 : 0;

  $OptionsInfo{Output} = $Options{output};
  $OptionsInfo{OutputStyle} = $Options{outputstyle};

  $OptionsInfo{RowsOutput} = ($Options{outputstyle} =~ /^FormulaRows$/i) ? 1 : 0;
  $OptionsInfo{FileOutput} = ($Options{output} =~ /^File$/i) ? 1 : 0;
  $OptionsInfo{CheckFormula} = $Options{fast} ? 0 : 1;

  $OptionsInfo{Precision} = $Options{precision};

  $OptionsInfo{OutFileRoot} = defined $Options{root} ? $Options{root} : undef;
  $OptionsInfo{ValueLabels} = defined $Options{valuelabels} ? $Options{valuelabels} : undef;

  # Setup formulas...
  @{$OptionsInfo{SpecifiedFormulas}} = ();
  if (@ARGV >= 1) {
    push @{$OptionsInfo{SpecifiedFormulas}}, @ARGV;
  }
  else {
    # Setup mode specified default values...
    push @{$OptionsInfo{SpecifiedFormulas}}, 'H2O';
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
  my($Index, $Value, $Label);
  # Set up labels for calculated values...
  @{$OptionsInfo{SpecifiedValueLabels}} = ();
  if ($Options{valuelabels}) {
    my($Value, $Label, @ValueLabels);
    @ValueLabels = split /\,/, $Options{valuelabels};
    if (@ValueLabels % 2) {
      die "Error: The value specified, $Options{valuelabels}, for option \"-v --valuelabels\" is not valid: It must contain even number of comma delimited values\n";
    }
    for ($Index = 0; $Index < @ValueLabels; $Index +=2) {
      $Value = $ValueLabels[$Index];
      $Value =~ s/ //g;
      $Label = $ValueLabels[$Index + 1];
      if ($Value !~ /^(ElementalAnalysis|MolecularWeight|ExactMass)$/i) {
	die "Error: The value specified, $Value, using option \"-v --valuelabels\" is not valid. Allowed values: ElementalAnalysis, MolecularWeight, or ExactMass\n";
      }
      push @{$OptionsInfo{SpecifiedValueLabels}}, ($Value, $Label);
    }
  }

  # Set up calculation type to label map...
  %{$OptionsInfo{ValueLabelsMap}} = (ElementalAnalysis => 'ElementalAnalysis', MolecularWeight => 'MolecularWeight', ExactMass => 'ExactMass');
  if (@{$OptionsInfo{SpecifiedValueLabels}}) {
    for ($Index = 0; $Index < @{$OptionsInfo{SpecifiedValueLabels}}; $Index +=2) {
      $Value = $OptionsInfo{SpecifiedValueLabels}[$Index];
      $Label = $OptionsInfo{SpecifiedValueLabels}[$Index + 1];
      if (exists $OptionsInfo{ValueLabelsMap}{$Value}) {
	$OptionsInfo{ValueLabelsMap}{$Value} = $Label;
      }
    }
  }

  # Setup output file name...
  $OptionsInfo{OutFileName} = '';
  if ($OptionsInfo{FileOutput}) {
    my($OutFileRoot, $OutFileExt);

    $OutFileRoot = '';
    $OutFileExt = "csv";
    if ($Options{outdelim} =~ /^tab$/i) {
      $OutFileExt = "tsv";
    }
    if ($Options{root}) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($Options{root});
      if ($RootFileName && $RootFileExt) {
	$OutFileRoot = $RootFileName;
      }
      else {
	$OutFileRoot = $Options{root};
      }
    }
    else {
      $OutFileRoot = 'FormulasElementalAnalysis';
    }
    $OptionsInfo{OutFileName} = $OutFileRoot . '.' . $OutFileExt;
    if (!$Options{overwrite}) {
      if (-e $OptionsInfo{OutFileName}) {
	die "Error: Output file, $OptionsInfo{OutFileName}, already exists.\nUse \-o --overwrite\ option or specify a different name using \"-r --root\" option.\n";
      }
    }
  }
}


# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{mode} = "All";
  $Options{outdelim} = "comma";
  $Options{output} = "STDOUT";
  $Options{outputstyle} = "FormulaBlock";
  $Options{precision} = 2;
  $Options{quote} = "yes";

  if (!GetOptions(\%Options, "help|h", "fast", "mode|m=s", "outdelim=s", "output=s", "outputstyle=s", "overwrite|o", "precision=i", "quote|q=s", "root|r=s", "valuelabels|v=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{outdelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{outdelim}, for option \"--outdelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if ($Options{output} !~ /^(STDOUT|File)$/i) {
    die "Error: The value specified, $Options{output}, for option \"--output\" is not valid. Allowed values: STDOUT or File\n";
  }
  if ($Options{outputstyle} !~ /^(FormulaBlock|FormulaRows)$/i) {
    die "Error: The value specified, $Options{outputstyle}, for option \"--outputstyle\" is not valid. Allowed values: FormulaBlock or FormulaRows\n";
  }
  if (!IsPositiveInteger($Options{precision})) {
    die "Error: The value specified, $Options{precision}, for option \"-p --precision\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{quote} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not valid. Allowed values: yes or no\n";
  }
}

__END__

=head1 NAME

ElementalAnalysis.pl - Perform elemental analysis using specified formulas

=head1 SYNOPSIS

ElementalAnalysis.pl Formula(s)...

ElementalAnalysis.pl [B<-h, --help>]
[B<-m, --mode> All | "ElementalAnalysis, [MolecularWeight, ExactMass]"]
[B<--outdelim> comma | tab | semicolon]
[B<--output> STDOUT | File] [B<--outputstyle> FormulaBlock | FormulaRows]
[B<-o, --overwrite>] [B<--precision> number] [B<-q, --quote> yes | no] [B<-r, --root> rootname]
[B<-v --valuelabels> [Name, Label, [Name, Label,...]]
[B<-w, --workingdir> dirname] Formula(s)...

=head1 DESCRIPTION

Perform elemental analysis using molecular formula(s) specified on the command line.

In addition to straightforward molecular formulas - H2O, HCl, C3H7O2N -
other supported variations are: Ca3(PO4)2, [PCl4]+, [Fe(CN)6]4-, C37H42N2O6+2, Na2CO3.10H2O,
8H2S.46H2O, and so on. Charges are simply ignored. Isotope symbols in formulas specification, including
D and T, are not supported.

=head1 PARAMETERS

=over 4

=item B<Formulas> I<Formula1 [Formula2...]>

I<Formulas> is a space delimited list of molecular formulas to use for elemental analysis.

Input value format is: I<Formula1 [Formula2 Formula3...]>. Default: I<H2O>.
Examples:

    HCl
    HCl, C3H7O2N
    H2O2 Ca3(PO4)2 [PCl4]+

=back

=head1 OPTIONS

=over 4

=item B<-h, --help>

Print this help message.

=item B<--fast>

In this mode, the specified formulas are considered valid and initial formula
validation check is skipped.

=item B<-m, --mode> I<All | "ElementalAnalysis,[MolecularWeight,ExactMass]">

Specify what values to calculate using molecular formulas specified on command
line: calculate all supported values or specify a comma delimited list of values. Possible
values: I<All | "ElementalAnalysis, [MolecularWeight, ExactMass]">. Default: I<All>.

=item B<--outdelim> I<comma | tab | semicolon>

Output text file delimiter. Possible values: I<comma, tab, or semicolon>
Default value: I<comma>.

=item B<--output> I<STDOUT | File>

List information at STDOUT or write it to a file. Possible values: I<STDOUT or File>. Default:
I<STDOUT>. B<-r, --root> option is used to generate output file name.

=item B<--outputstyle> I<FormulaBlock | FormulaRows>

Specify how to list calculated values: add a new line for each property and present it as a block
for each formula; or include all properties in one line and show it as a single line.

Possible values: I<FormulaBlock | FormulaRows>. Default: I<FormulaBlock>

An example for I<FormulaBlock> output style:

    Formula: H2O
    ElementalAnalysis: H: H: 11.1898%; O: 88.8102%
    MolecularWeight: 18.0153
    ExactMass: 18.0106
    ... ...
    ... ...
    ... ...

    Formula: H2O2
    ElementalAnalysis: H: 5.9265%; O: 94.0735%
    MolecularWeight: 34.0147
    ExactMass: 34.0055
    ... ...
    ... ...
    ... ...

An example for I<FormulaRows> output style:

    Formula,ElementalAnalysis,MolecularWeight,ExactMass
    H2O,H: 11.1898%; O: 88.8102%,18.0153,18.0106
    H2O2,H: 5.9265%; O: 94.0735%,34.0147,34.0055

=item B<-o, --overwrite>

Overwrite existing files.

=item B<--precision> I<number>

Precision for listing numerical values. Default: up to I<4> decimal places.
Valid values: positive integers.

=item B<-r, --root> I<rootname>

New text file name is generated using the root: <Root>.<Ext>. File name is only
used during I<File> value of B<-o, --output> option.

Default file name: FormulsElementalAnalysis.<Ext>. The csv, and tsv
<Ext> values are used for comma/semicolon, and tab delimited text files respectively.

=item B<-v --valuelabels> I<Name,Label,[Name,Label,...]>

Specify labels to use for calculated values. In general, it's a comma delimited
list of value name and column label pairs. Supported value names: I<ElementalAnalysis,
MolecularWeight,  and ExactMass>. Default labels: I<ElementalAnalysis, MolecularWeight,
and ExactMass>.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To perform elemental analysis, calculate molecular weight and exact mass for H2O,
type:

    % ElementalAnalysis.pl

To perform elemental analysis, calculate molecular weight and exact mass for
Ca3(PO4)2 and [PCl4]+, type:

    % ElementalAnalysis.pl "Ca3(PO4)2" "[PCl4]+"

To perform elemental analysis, use label analysis for calculated data, and generate a
new CSV file ElementalAnalysis.csv for H2O and H2O2, type:

    % ElementalAnalysis.pl --m ElementalAnalysis --output File
      --valuelabels "ElementalAnalysis,Analysis" -o -r ElementalAnalysis.csv
      H2O H2O2

To calculate molecular weight and exact mass with four decimal precision and
generate a new CSV file WeightAndMass.csv with data rows for H2O and H2O2, type:

    % ElementalAnalysis.pl --m "MolecularWeight,ExactMass" --output File
      --outputstyle FormulaRows -o -r WeightAndMass.csv
      H2O H2O2

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

ElementalAnalysisSDFiles.pl, ElementalAnalysisTextFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
