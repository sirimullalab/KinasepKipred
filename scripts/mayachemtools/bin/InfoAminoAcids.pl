#!/usr/bin/perl -w
#
# File: InfoAminoAcids.pl
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
use AminoAcids;

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

ListAminoAcidProperties();
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# List data for an amino acid...
sub ListAminoAcidData {
  my($DataLabelRef, $DataValueRef) = @_;
  my($Index, $Line, $Value);

  if ($OptionsInfo{AminoAcidRowsOutput}) {
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
      $Line = $DataLabelRef->[$Index] . ' ' . $DataValueRef->[$Index];
      if ($OptionsInfo{FileOutput}) {
	print OUTFILE "$Line\n";
      }
      else {
	print "$Line\n";
      }
    }
  }
}

# List data for an amino acid...
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

# List properties for amino acids...
sub ListAminoAcidProperties {
  my($AminoAcidID, $AminoAcidDataRef, $PropertyName, $PropertyValue, @PropertyLabels, @PropertyValues);

  print "Listing information for amino acid(s)...\n";

  if ($OptionsInfo{FileOutput}) {
    print "Generating file $OptionsInfo{OutFileName}...\n";
    open OUTFILE, ">$OptionsInfo{OutFileName}" or die "Couldn't open $OptionsInfo{OutFileName}: $!\n";
  }

  # Setup property labels...
  @PropertyLabels = ();
  for $PropertyName (@{$OptionsInfo{SpecifiedProperies}}) {
    push @PropertyLabels, ("$PropertyName:");
  }

  if ($OptionsInfo{AminoAcidRowsOutput}) {
    ListHeaderRowData(\@PropertyLabels);
  }

  # Go over specified properties...
  for $AminoAcidID (@{$OptionsInfo{SpecifiedAminoAcidIDs}}) {
    $AminoAcidDataRef = AminoAcids::GetAminoAcidPropertiesData($AminoAcidID);

    if (!$OptionsInfo{AminoAcidRowsOutput}) {
      if ($OptionsInfo{FileOutput}) {
	print OUTFILE "\nListing properties for amino acid $AminoAcidID...\n\n";
      }
      else {
	print "\nListing properties for amino acid $AminoAcidID...\n\n";
      }
    }

    # Collect data..
    @PropertyValues = ();
    for $PropertyName (@{$OptionsInfo{SpecifiedProperies}}) {
      $PropertyValue = $AminoAcidDataRef->{$PropertyName};
      if (IsFloat($PropertyValue)) {
	$PropertyValue = sprintf("%.$OptionsInfo{Precision}f", $PropertyValue) + 0;
      }
      push @PropertyValues, $PropertyValue;
    }
    # List data...
    ListAminoAcidData(\@PropertyLabels, \@PropertyValues);
  }
  if ($OptionsInfo{FileOutput}) {
    close OUTFILE;
  }
  print "\n";
}

# Get propery names from categories...
sub GetPropertyNamesFromCategories {
  my($CategoryName) = @_;
  my(@PropertyNames);

  @PropertyNames = ();
  if ($CategoryName =~ /^Basic$/i) {
    @PropertyNames = ('ThreeLetterCode', 'OneLetterCode', 'AminoAcid', 'DNACodons', 'RNACodons', 'ChemicalFormula','MolecularWeight', 'LinearStructure', 'LinearStructureAtpH7.4');
  } elsif ($CategoryName =~ /^BasicPlus$/i) {
    @PropertyNames = ('ThreeLetterCode', 'OneLetterCode', 'AminoAcid', 'DNACodons', 'RNACodons', 'AcidicBasic', 'PolarNonpolar', 'Charged', 'Aromatic', 'HydrophobicHydophilic', 'IsoelectricPoint', 'pKCOOH', 'pKNH3+', 'ChemicalFormula', 'MolecularWeight', 'ExactMass', 'ChemicalFormulaMinusH2O', 'MolecularWeightMinusH2O(18.01524)', 'ExactMassMinusH2O(18.01056)','LinearStructure', 'LinearStructureAtpH7.4');
  } elsif ($CategoryName =~ /^BasicAndHydrophobicity$/i) {
    @PropertyNames = ('ThreeLetterCode', 'OneLetterCode', 'AminoAcid', 'DNACodons', 'RNACodons', 'ChemicalFormula', 'MolecularWeight', 'LinearStructure', 'LinearStructureAtpH7.4', 'HydrophobicityEisenbergAndOthers', 'HydrophobicityHoppAndWoods', 'HydrophobicityJanin', 'HydrophobicityKyteAndDoolittle', 'HydrophobicityRoseAndOthers', 'HydrophobicityWolfendenAndOthers');
  } elsif ($CategoryName =~ /^BasicAndHydrophobicityPlus$/i) {
    @PropertyNames = ('ThreeLetterCode', 'OneLetterCode', 'AminoAcid', 'DNACodons', 'RNACodons', 'ChemicalFormula', 'MolecularWeight', 'LinearStructure', 'LinearStructureAtpH7.4', 'HydrophobicityAbrahamAndLeo', 'HydrophobicityBlack', 'HydrophobicityBullAndBreese', 'HydrophobicityChothia', 'HydrophobicityEisenbergAndOthers', 'HydrophobicityFauchereAndOthers', 'HydrophobicityGuy', 'HydrophobicityHPLCAtpH3.4Cowan', 'HydrophobicityHPLCAtpH7.5Cowan', 'HydrophobicityHPLCParkerAndOthers', 'HydrophobicityHPLCWilsonAndOthers', 'HydrophobicityHoppAndWoods', 'HydrophobicityJanin', 'HydrophobicityKyteAndDoolittle', 'HydrophobicityManavalanAndOthers', 'HydrophobicityMiyazawaAndOthers', 'HydrophobicityOMHSweetAndOthers', 'HydrophobicityRaoAndArgos', 'HydrophobicityRfMobility', 'HydrophobicityRoseAndOthers', 'HydrophobicityRoseman', 'HydrophobicityWellingAndOthers', 'HydrophobicityWolfendenAndOthers');
  }

  return @PropertyNames;
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /^tab$/i ) ? "\t" : (($Options{outdelim} =~ /^semicolon$/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^yes$/i) ? 1 : 0;

  $OptionsInfo{Overwrite} = defined $Options{overwrite} ? $Options{overwrite} : undef;
  $OptionsInfo{OutFileRoot} = defined $Options{root} ? $Options{root} : undef;

  $OptionsInfo{OutputStyle} = $Options{outputstyle};

  $OptionsInfo{AminoAcidRowsOutput} = ($Options{outputstyle} =~ /^AminoAcidRows$/i) ? 1 : 0;
  $OptionsInfo{FileOutput} = ($Options{output} =~ /^File$/i) ? 1 : 0;

  $OptionsInfo{Precision} = $Options{precision};

  my($AminoAcidID, @AminoAcidIDs);

  @{$OptionsInfo{SpecifiedAminoAcidIDs}} = ();

  # Set up Amino Acids IDs except for All mode...
  @AminoAcidIDs = ();

  if (@ARGV >= 1) {
    push @AminoAcidIDs, @ARGV;
  }
  else {
    # Setup mode specified default values...
    push @AminoAcidIDs, 'Ala';
  }

  # Generate list of amino acids...
  if (@ARGV == 1 && $ARGV[0] =~ /^All$/i) {
    push @{$OptionsInfo{SpecifiedAminoAcidIDs}}, AminoAcids::GetAminoAcids();
  }
  else {
    ID: for $AminoAcidID (@AminoAcidIDs) {
      if (AminoAcids::IsAminoAcid($AminoAcidID)) {
	push @{$OptionsInfo{SpecifiedAminoAcidIDs}}, $AminoAcidID;
      }
      else {
	warn "Ignoring amino acid ID, $AminoAcidID, specified using command line parameter option: Unknown amino acid ID...\n";
	next ID;
      }
    }
  }
  SetupSpecifiedProperties();

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
      $OutFileRoot = 'AminoAcidsInfo';
    }
    $OptionsInfo{OutFileName} = $OutFileRoot . '.' . $OutFileExt;
    if (!$Options{overwrite}) {
      if (-e $OptionsInfo{OutFileName}) {
	die "Error: Output file, $OptionsInfo{OutFileName}, already exists.\nUse \-o --overwrite\ option or specify a different name using \"-r --root\" option.\n";
      }
    }
  }
}

# Setup properties to list...
sub SetupSpecifiedProperties {

  $OptionsInfo{Properties} = defined $Options{properties} ? $Options{properties} : undef;

  $OptionsInfo{PropertiesMode} = $Options{propertiesmode};
  $OptionsInfo{PropertiesListing} = $Options{propertieslisting};

  # Make sure appropriate properties/category names are specified...
  @{$OptionsInfo{SpecifiedProperies}} = ();
  if ($Options{properties} && ($Options{propertiesmode} =~ /^All$/i) ) {
    warn "Warning: Ignoring values specifed by \"-p --properties\" option: Not valid for All value of \"--propertiesmode\" option...\n";
  }
  if ($Options{propertiesmode} =~ /^All$/i) {
    if ($Options{propertieslisting} =~ /^Alphabetical$/i) {
      push @{$OptionsInfo{SpecifiedProperies}}, AminoAcids::GetAminoAcidPropertiesNames('Alphabetical');
    }
    else {
      push @{$OptionsInfo{SpecifiedProperies}}, AminoAcids::GetAminoAcidPropertiesNames();
    }
  }
  else {
    if ($Options{properties}) {
      if ($Options{propertiesmode} =~ /^Categories$/i) {
	# Check category name...
	if ($Options{properties} !~ /^(Basic|BasicPlus|BasicAndHydrophobicity|BasicAndHydrophobicityPlus)$/i) {
	  die "Error: The value specified, $Options{properties}, for option \"-p --properties\" in conjunction with \"Categories\" value for option \"--propertiesmode\" is not valid. Allowed values: Basic, BasicPlus, BasicAndHydrophobicity, and BasicAndHydrophobicityPlus\n";
	}
	# Set propertynames...
	push @{$OptionsInfo{SpecifiedProperies}}, GetPropertyNamesFromCategories($Options{properties});
      }
      else {
	# Check property names..
	my($Name, $PropertyName, @Names);
	@Names = split /\,/, $Options{properties};
	NAME: for $Name (@Names) {
	  $PropertyName = RemoveLeadingAndTrailingWhiteSpaces($Name);
	  if (AminoAcids::IsAminoAcidProperty($PropertyName)) {
	    push @{$OptionsInfo{SpecifiedProperies}}, $PropertyName;
	  }
	  else {
	    warn "Warning: Ignoring value, $Name, specifed by \"-p --properties\" option: Unknown property name...\n";
	  }
	}
	if ($Options{propertieslisting} =~ /^Alphabetical$/i) {
	  # ThreeLetterCode, OneLetterCode and AminoAcid are always listed first...
	  # NaturalIsotopeData in the end...
	  my($OneLetterCodePresent, $ThreeLetterCodePresent, $AminoAcidPresent,  @AlphabeticalProperties, %PropertiesMap);
	  %PropertiesMap = ();
	  @AlphabeticalProperties = ();
	  $OneLetterCodePresent = 0; $ThreeLetterCodePresent = 0; $AminoAcidPresent = 0;
	  NAME: for $Name (@{$OptionsInfo{SpecifiedProperies}}) {
	    if ($Name =~ /^OneLetterCode$/i) {
	      $OneLetterCodePresent = 1;
	      next NAME;
	    }
	    if ($Name =~ /^ThreeLetterCode$/i) {
	      $ThreeLetterCodePresent = 1;
	      next NAME;
	    }
	    if ($Name =~ /^AminoAcid$/i) {
	      $AminoAcidPresent = 1;
	      next NAME;
	    }
	    $PropertiesMap{$Name} = $Name;
	  }
	  # Setup the alphabetical list...
	  if ($ThreeLetterCodePresent) {
	    push @AlphabeticalProperties, 'ThreeLetterCode';
	  }
	  if ($OneLetterCodePresent) {
	    push @AlphabeticalProperties, 'OneLetterCode';
	  }
	  if ($AminoAcidPresent) {
	    push @AlphabeticalProperties, 'AminoAcid';
	  }
	  for $Name (sort keys %PropertiesMap) {
	    push @AlphabeticalProperties, $Name;
	  }
	  @{$OptionsInfo{SpecifiedProperies}} = ();
	  push @{$OptionsInfo{SpecifiedProperies}}, @AlphabeticalProperties;
	}
      }
    }
    else {
      # Set default value...
      push @{$OptionsInfo{SpecifiedProperies}}, GetPropertyNamesFromCategories('Basic');
    }
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{outdelim} = "comma";
  $Options{output} = "STDOUT";
  $Options{outputstyle} = "AminoAcidBlock";
  $Options{precision} = 4;
  $Options{propertiesmode} = "Categories";
  $Options{propertieslisting} = "ByGroup";
  $Options{quote} = "yes";

  if (!GetOptions(\%Options, "help|h", "outdelim=s", "output=s", "outputstyle=s", "overwrite|o", "precision=i", "properties|p=s", "propertieslisting=s", "propertiesmode=s", "quote|q=s", "root|r=s", "workingdir|w=s")) {
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
  if ($Options{outputstyle} !~ /^(AminoAcidBlock|AminoAcidRows)$/i) {
    die "Error: The value specified, $Options{outputstyle}, for option \"--outputstyle\" is not valid. Allowed values: AminoAcidBlock or AminoAcidRows\n";
  }
  if (!IsPositiveInteger($Options{precision})) {
    die "Error: The value specified, $Options{precision}, for option \"-p --precision\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{propertiesmode} !~ /^(Categories|Names|All)$/i) {
    die "Error: The value specified, $Options{propertiesmode}, for option \"--propertiesmode\" is not valid. Allowed values: Categories, Names, or All\n";
  }
  if ($Options{propertieslisting} !~ /^(ByGroup|Alphabetical)$/i) {
    die "Error: The value specified, $Options{propertieslisting}, for option \"--propertieslisting\" is not valid. Allowed values: ByGroup, or Alphabetical\n";
  }
  if ($Options{quote} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not valid. Allowed values: yes or no\n";
  }
}

__END__

=head1 NAME

InfoAminoAcids.pl - List properties of amino acids

=head1 SYNOPSIS

InfoAminoAcids.pl AminoAcidIDs...

InfoAminoAcids.pl [B<-h, --help>] [B<--outdelim> comma | tab | semicolon]
[B<--output> STDOUT | File] [B<--outputstyle> AminoAcidBlock | AminoAcidRows]
[B<-o, --overwrite>] [B<--precision> number] [B<--propertiesmode> Categories | Names | All]
[B<-p, --properties> CategoryName,[CategoryName,...] | PropertyName,[PropertyName,...]]
[B<--propertieslinting> ByGroup | Alphabetical] [B<-q, --quote> yes | no] [B<-r, --root> rootname]
[B<-w, --workingdir> dirname] AminoAcidIDs...

=head1 DESCRIPTION

List amino acid properties. Amino acids identification supports these three types of IDs: one letter
code, three letter code or name. Amino acid properties data, in addition to basic information about
amino acids - one and three letter codes, name, DNA and RNA codons, molecular weight - include
variety of other properties: polarity, acidity, hydrophobicity, and so on.

=head1 PARAMETERS

=over 4

=item B<AminoAcidIDs> I<ThreeLetterCode [OneLetterCode AminoAcidName...]>

I<AminoAcidIDs> is a space delimited list of values to identify amino acids.

Input value format is: I<ThreeLetterCode [OneLetterCode AminoAcidName...]>. Default: I<Ala>.
Examples:

    Ala
    Glu A
    Alanine Glu Y "Aspartic acid"

=back

=head1 OPTIONS

=over 4

=item B<-h, --help>

Print this help message.

=item B<--outdelim> I<comma | tab | semicolon>

Output text file delimiter. Possible values: I<comma, tab, or semicolon>
Default value: I<comma>.

=item B<--output> I<STDOUT | File>

List information at STDOUT or write it to a file. Possible values: I<STDOUT or File>. Default:
I<STDOUT>. B<-r, --root> option is used to generate output file name.

=item B<--outputstyle> I<AminoAcidBlock | AminoAcidRows>

Specify how to list amino acid information: add a new line for each property and present it as a block
for each amino acid; or include all properties in one line and show it as a single line.

Possible values: I<AminoAcidBlock | AminoAcidRows>. Default: I<AminoAcidBlock>

An example for I<AminoAcidBlock> output style:

    ThreeLetterCode: Ala
    OneLetterCode: A
    AminoAcid: Alanine
    MolecularWeight: 89.0941
    ... ...
    ... ...
    ... ...

    ThreeLetterCode: Glu
    OneLetterCode: E
    AminoAcid: Glutamic acid
    MolecularWeight: 147.1308
    ... ...
    ... ...
    ... ...

An example for I<AminoAcidRows> output style:

    ThreeLetterCode,OneLetterCode,AminoAcid,MolecularWeight
    Ala,A,Alanine,89.0941..
    Glu,E,Glutamic acid,147.1308..

=item B<-o, --overwrite>

Overwrite existing files.

=item B<--precision> I<number>

Precision for listing numerical values. Default: up to I<4> decimal places.
Valid values: positive integers.

=item B<--propertiesmode> I<Categories | Names | All>

Specify how property names are specified: use category names; explicit list of property names; or
use all available properties. Possible values: I<Categories, Names, or All>. Default: I<Categories>.

This option is used in conjunction with B<-p, --properties> option to specify properties of
interest.

=item B<-p, --properties> I<CategoryName,[CategoryName,...] | PropertyName,[PropertyName,...]>

This option is B<--propertiesmode> specific. In general, it's a list of comma separated category or
property names.

Specify which amino acid properties information to list for the amino acid IDs specified using command:
line parameters: list basic and/or hydrophobicity information; list all available information; or specify a comma
separated list of amino acid property names.

Possible values: I<Basic | BasicPlus | BasicAndHydrophobicity | BasicAndHydrophobicityPlus | PropertyName,[PropertyName,...]>.
Default: I<Basic>.

I<Basic> includes: I<ThreeLetterCode, OneLetterCode, AminoAcid, DNACodons, RNACodons, ChemicalFormula, MolecularWeight, LinearStructure, LinearStructureAtpH7.4>

I<BasicPlus> includes: I<ThreeLetterCode, OneLetterCode, AminoAcid, DNACodons, RNACodons, AcidicBasic, PolarNonpolar, Charged, Aromatic, HydrophobicHydophilic, IsoelectricPoint, pKCOOH, pKNH3+, ChemicalFormula, MolecularWeight, ExactMass, ChemicalFormulaMinusH2O, MolecularWeightMinusH2O(18.01524), ExactMassMinusH2O(18.01056), LinearStructure, LinearStructureAtpH7.4>

I<BasicAndHydrophobicity> includes: I<ThreeLetterCode, OneLetterCode, AminoAcid, DNACodons, RNACodons, ChemicalFormula, MolecularWeight, LinearStructure, LinearStructureAtpH7.4, HydrophobicityEisenbergAndOthers, HydrophobicityHoppAndWoods, HydrophobicityJanin, HydrophobicityKyteAndDoolittle, HydrophobicityRoseAndOthers, HydrophobicityWolfendenAndOthers>

I<BasicAndHydrophobicityPlus> includes: I<(ThreeLetterCode, OneLetterCode, AminoAcid, DNACodons, RNACodons, ChemicalFormula, MolecularWeight, LinearStructure, LinearStructureAtpH7.4, HydrophobicityAbrahamAndLeo, HydrophobicityBlack, HydrophobicityBullAndBreese, HydrophobicityChothia, HydrophobicityEisenbergAndOthers, HydrophobicityFauchereAndOthers, HydrophobicityGuy, HydrophobicityHPLCAtpH3.4Cowan, HydrophobicityHPLCAtpH7.5Cowan, HydrophobicityHPLCParkerAndOthers, HydrophobicityHPLCWilsonAndOthers, HydrophobicityHoppAndWoods, HydrophobicityJanin, HydrophobicityKyteAndDoolittle, HydrophobicityManavalanAndOthers, HydrophobicityMiyazawaAndOthers, HydrophobicityOMHSweetAndOthers, HydrophobicityRaoAndArgos, HydrophobicityRfMobility, HydrophobicityRoseAndOthers, HydrophobicityRoseman, HydrophobicityWellingAndOthers, HydrophobicityWolfendenAndOthers>

Here is a complete list of available properties: ThreeLetterCode, OneLetterCode, AminoAcid, DNACodons, RNACodons, AcidicBasic, PolarNonpolar, Charged, Aromatic, HydrophobicHydophilic, IsoelectricPoint, pKCOOH, pKNH3+, ChemicalFormula, MolecularWeight, ExactMass, ChemicalFormulaMinusH2O, MolecularWeightMinusH2O(18.01524), ExactMassMinusH2O(18.01056), vanderWaalsVolume, %AccessibleResidues, %BuriedResidues, AlphaHelixChouAndFasman, AlphaHelixDeleageAndRoux, AlphaHelixLevitt, AminoAcidsComposition, AminoAcidsCompositionInSwissProt, AntiparallelBetaStrand, AverageAreaBuried, AverageFlexibility, BetaSheetChouAndFasman, BetaSheetDeleageAndRoux, BetaSheetLevitt, BetaTurnChouAndFasman, BetaTurnDeleageAndRoux, BetaTurnLevitt, Bulkiness, CoilDeleageAndRoux, HPLCHFBARetention, HPLCRetentionAtpH2.1, HPLCRetentionAtpH7.4, HPLCTFARetention, HydrophobicityAbrahamAndLeo, HydrophobicityBlack, HydrophobicityBullAndBreese, HydrophobicityChothia, HydrophobicityEisenbergAndOthers, HydrophobicityFauchereAndOthers, HydrophobicityGuy, HydrophobicityHPLCAtpH3.4Cowan, HydrophobicityHPLCAtpH7.5Cowan, HydrophobicityHPLCParkerAndOthers, HydrophobicityHPLCWilsonAndOthers, HydrophobicityHoppAndWoods, HydrophobicityJanin, HydrophobicityKyteAndDoolittle, HydrophobicityManavalanAndOthers, HydrophobicityMiyazawaAndOthers, HydrophobicityOMHSweetAndOthers, HydrophobicityRaoAndArgos, HydrophobicityRfMobility, HydrophobicityRoseAndOthers, HydrophobicityRoseman, HydrophobicityWellingAndOthers, HydrophobicityWolfendenAndOthers, ParallelBetaStrand, PolarityGrantham, PolarityZimmerman, RatioHeteroEndToSide, RecognitionFactors, Refractivity, RelativeMutability, TotalBetaStrand, LinearStructure, LinearStructureAtpH7.4

=item B<--propertieslisting> I<ByGroup | Alphabetical>

Specify how to list properties for amino acids: group by category or an alphabetical by
property names. Possible values: I<ByGroup or Alphabetical>. Default: I<ByGroup>.

=item B<-q, --quote> I<yes | no>

Put quotes around column values in output text file. Possible values: I<yes or
no>. Default value: I<yes>.

=item B<-r, --root> I<rootname>

New text file name is generated using the root: <Root>.<Ext>. File name is only
used during I<File> value of B<-o, --output> option.

Default file name: AminoAcidInfo<mode>.<Ext>. The csv, and tsv
<Ext> values are used for comma/semicolon, and tab delimited text files respectively.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To list basic properties information for amino acid Ala, type:

    % InfoAminoAcids.pl

To list all available properties information for amino acid Ala, type:

    % InfoAminoAcids.pl --propertiesmode all

To list basic properties information for amino acids Ala, Arg, and Asp type:

    % InfoAminoAcids.pl Ala Arg Asp
    % InfoAminoAcids.pl A Arg "Aspartic acid"

To list all available properties information for amino acids Ala, Arg, and Asp type:

    % InfoAminoAcids.pl --propertiesmode all Ala Arg Asp

To list basic and hydrophobicty properties information for amino acids Ala, Arg, and Asp type:

    % InfoAminoAcids.pl --propertiesmode Categories
      --properties BasicAndHydrophobicity Ala Arg Asp

To list OneLetterCode, ThreeLetterCode, DNACodons, and MolecularWeight for amino
acids Ala, Arg, and Asp type:

    % InfoAminoAcids.pl --propertiesmode Names
      --properties OneLetterCode,ThreeLetterCode,DNACodons,MolecularWeight
      Ala Arg Asp

To alphabetically list basic and hydrophobicty properties information for amino acids Ala, Arg, and Asp
in rows insetad of amino acid blocks with quotes around the values, type:

    % InfoAminoAcids.pl --propertiesmode Categories
      --properties BasicAndHydrophobicity --propertieslisting alphabetical
      --outdelim comma --outputstyle AminoAcidRows --quote yes Ala Arg Asp

To alphabetically list basic and hydrophobicty properties information for amino acids Ala, Arg, and Asp
in rows insetad of amino acid blocks with quotes around the values and write them into a file
AminoAcidProperties.csv, type:

    % InfoAminoAcids.pl --propertiesmode Categories
      --properties BasicAndHydrophobicity --propertieslisting alphabetical
      --outdelim comma --outputstyle AminoAcidRows --quote yes
      --output File -r AminoAcidProperties -o Ala Arg Asp

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoNucleicAcids.pl InfoPeriodicTableElements.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
