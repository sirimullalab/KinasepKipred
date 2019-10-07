#!/usr/bin/perl -w
#
# File: InfoNucleicAcids.pl
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
use NucleicAcids;

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

ListNucleicAcidProperties();
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# List data for an nucleic acid...
sub ListNucleicAcidData {
  my($DataLabelRef, $DataValueRef) = @_;
  my($Index, $Line, $Value);

  if ($OptionsInfo{NucleicAcidRowsOutput}) {
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

# List data for an nucleic acid...
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

# List properties for nucleic acids...
sub ListNucleicAcidProperties {
  my($NucleicAcidID, $NucleicAcidDataRef, $PropertyName, $PropertyValue, @PropertyLabels, @PropertyValues);

  print "Listing information for nucleic acid(s)...\n";

  if ($OptionsInfo{FileOutput}) {
    print "Generating file $OptionsInfo{OutFileName}...\n";
    open OUTFILE, ">$OptionsInfo{OutFileName}" or die "Couldn't open $OptionsInfo{OutFileName}: $!\n";
  }

  # Setup property labels...
  @PropertyLabels = ();
  for $PropertyName (@{$OptionsInfo{SpecifiedProperies}}) {
    push @PropertyLabels, ("$PropertyName:");
  }

  if ($OptionsInfo{NucleicAcidRowsOutput}) {
    ListHeaderRowData(\@PropertyLabels);
  }

  # Go over specified properties...
  for $NucleicAcidID (@{$OptionsInfo{SpecifiedNucleicAcidIDs}}) {
    $NucleicAcidDataRef = NucleicAcids::GetNucleicAcidPropertiesData($NucleicAcidID);

    if (!$OptionsInfo{NucleicAcidRowsOutput}) {
      if ($OptionsInfo{FileOutput}) {
	print OUTFILE "\nListing properties for nucleic acid $NucleicAcidID...\n\n";
      }
      else {
	print "\nListing properties for nucleic acid $NucleicAcidID...\n\n";
      }
    }

    # Collect data..
    @PropertyValues = ();
    for $PropertyName (@{$OptionsInfo{SpecifiedProperies}}) {
      $PropertyValue = $NucleicAcidDataRef->{$PropertyName};
      if (IsFloat($PropertyValue)) {
	$PropertyValue = sprintf("%.$OptionsInfo{Precision}f", $PropertyValue) + 0;
      }
      push @PropertyValues, $PropertyValue;
    }
    # List data...
    ListNucleicAcidData(\@PropertyLabels, \@PropertyValues);
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
    @PropertyNames = ('Code', 'OtherCodes', 'Name', 'Type', 'MolecularFormula', 'MolecularWeight');
  } elsif ($CategoryName =~ /^BasicPlus$/i) {
    @PropertyNames = ('Code', 'OtherCodes', 'Name', 'Type', 'MolecularFormula', 'MolecularWeight', 'ExactMass', 'ElementalComposition');
  }

  return @PropertyNames;
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /^tab$/i ) ? "\t" : (($Options{outdelim} =~ /^semicolon$/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^yes$/i) ? 1 : 0;

  $OptionsInfo{Overwrite} = defined $Options{overwrite} ? $Options{overwrite} : undef;
  $OptionsInfo{OutFileRoot} = defined $Options{root} ? $Options{root} : undef;

  $OptionsInfo{Output} = $Options{output};
  $OptionsInfo{OutputStyle} = $Options{outputstyle};

  $OptionsInfo{NucleicAcidRowsOutput} = ($Options{outputstyle} =~ /^NucleicAcidRows$/i) ? 1 : 0;
  $OptionsInfo{FileOutput} = ($Options{output} =~ /^File$/i) ? 1 : 0;

  $OptionsInfo{Precision} = $Options{precision};

  my($NucleicAcidID, @NucleicAcidIDs);

  @{$OptionsInfo{SpecifiedNucleicAcidIDs}} = ();

  # Set up Nucleic Acids IDs except for All mode...
  @NucleicAcidIDs = ();

  if (@ARGV >= 1) {
    push @NucleicAcidIDs, @ARGV;
  }
  else {
    # Setup mode specified default values...
    if ($Options{mode} =~ /NucleicAcidID/i) {
      push @NucleicAcidIDs, 'A';
    }
    elsif ($Options{mode} =~ /NucleicAcidType/i) {
      push @NucleicAcidIDs, 'Nucleoside';
    }
    else {
      push @NucleicAcidIDs, 'A';
    }
  }

  # Generate list of nucleic acids...
  if (@ARGV == 1 && $ARGV[0] =~ /^All$/i) {
    push @{$OptionsInfo{SpecifiedNucleicAcidIDs}}, NucleicAcids::GetNucleicAcids();
  }
  else {
    if ($Options{mode} =~ /NucleicAcidID/i) {
      ID: for $NucleicAcidID (@NucleicAcidIDs) {
	if (NucleicAcids::IsNucleicAcid($NucleicAcidID)) {
	  push @{$OptionsInfo{SpecifiedNucleicAcidIDs}}, $NucleicAcidID;
	}
	else {
	  warn "Ignoring nucleic acid ID, $NucleicAcidID, specified using command line parameter option: Unknown nucleic acid ID...\n";
	  next ID;
	}
      }
    }
    elsif ($Options{mode} =~ /NucleicAcidType/i) {
      ID: for $NucleicAcidID (@NucleicAcidIDs) {
	  if (!NucleicAcids::IsNucleicAcidType($NucleicAcidID)) {
	    warn "Ignoring nucleic acid type, $NucleicAcidID, specified using command line parameter option: Unknown nucleic acid type...\n";
	    next ID;
	  }
	  push @{$OptionsInfo{SpecifiedNucleicAcidIDs}}, NucleicAcids::GetNucleicAcidsByType($NucleicAcidID);
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
      $OutFileRoot = 'NucleicAcidsInfo';
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
      push @{$OptionsInfo{SpecifiedProperies}}, NucleicAcids::GetNucleicAcidPropertiesNames('Alphabetical');
    }
    else {
      push @{$OptionsInfo{SpecifiedProperies}}, NucleicAcids::GetNucleicAcidPropertiesNames();
    }
  }
  else {
    if ($Options{properties}) {
      if ($Options{propertiesmode} =~ /^Categories$/i) {
	# Check category name...
	if ($Options{properties} !~ /^(Basic|BasicPlus)$/i) {
	  die "Error: The value specified, $Options{properties}, for option \"-p --properties\" in conjunction with \"Categories\" value for option \"--propertiesmode\" is not valid. Allowed values: Basic and BasicPlus\n";
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
	  if (NucleicAcids::IsNucleicAcidProperty($PropertyName)) {
	    push @{$OptionsInfo{SpecifiedProperies}}, $PropertyName;
	  }
	  else {
	    warn "Warning: Ignoring value, $Name, specifed by \"-p --properties\" option: Unknown property name...\n";
	  }
	}
	if ($Options{propertieslisting} =~ /^Alphabetical$/i) {
	  # Code, OtherCodes and Name are always listed first...
	  my($CodePresent, $OtherCodesPresent, $NamePresent,  @AlphabeticalProperties, %PropertiesMap);
	  %PropertiesMap = ();
	  @AlphabeticalProperties = ();
	  $CodePresent = 0; $OtherCodesPresent = 0; $NamePresent = 0;
	  NAME: for $Name (@{$OptionsInfo{SpecifiedProperies}}) {
	    if ($Name =~ /^Code$/i) {
	      $CodePresent = 1;
	      next NAME;
	    }
	    if ($Name =~ /^OtherCodes$/i) {
	      $OtherCodesPresent = 1;
	      next NAME;
	    }
	    if ($Name =~ /^Name$/i) {
	      $NamePresent = 1;
	      next NAME;
	    }
	    $PropertiesMap{$Name} = $Name;
	  }
	  # Setup the alphabetical list...
	  if ($CodePresent) {
	    push @AlphabeticalProperties, 'Code';
	  }
	  if ($OtherCodesPresent) {
	    push @AlphabeticalProperties, 'OtherCodesPresent';
	  }
	  if ($NamePresent) {
	    push @AlphabeticalProperties, 'Name';
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
  $Options{mode} = "NucleicAcidID";
  $Options{outdelim} = "comma";
  $Options{output} = "STDOUT";
  $Options{outputstyle} = "NucleicAcidBlock";
  $Options{precision} = 4;
  $Options{propertiesmode} = "Categories";
  $Options{propertieslisting} = "ByGroup";
  $Options{quote} = "yes";

  if (!GetOptions(\%Options, "help|h", "mode|m=s", "outdelim=s", "output=s", "outputstyle=s", "overwrite|o", "precision=i", "properties|p=s", "propertieslisting=s", "propertiesmode=s", "quote|q=s", "root|r=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{mode} !~ /^(NucleicAcidID|NucleicAcidType)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"--mode\" is not valid. Allowed values: NucleicAcidID or NucleicAcidType\n";
  }
  if ($Options{outdelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{outdelim}, for option \"--outdelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if ($Options{output} !~ /^(STDOUT|File)$/i) {
    die "Error: The value specified, $Options{output}, for option \"--output\" is not valid. Allowed values: STDOUT or File\n";
  }
  if ($Options{outputstyle} !~ /^(NucleicAcidBlock|NucleicAcidRows)$/i) {
    die "Error: The value specified, $Options{outputstyle}, for option \"--outputstyle\" is not valid. Allowed values: NucleicAcidBlock or NucleicAcidRows\n";
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

InfoNucleicAcids.pl - List properties of nucleic acids

=head1 SYNOPSIS

InfoNucleicAcids.pl NucleicAcidIDs...

InfoNucleicAcids.pl [B<-h, --help>] [B<-m, --mode> NucleicAcidID | NucleicAcidType]
[B<--OutDelim> comma | tab | semicolon]
[B<--output> STDOUT | File] [B<--OutputStyle> NucleicAcidBlock | NucleicAcidRows]
[B<-o, --overwrite>] [B<--precision> number] [B<--PropertiesMode> Categories | Names | All]
[B<-p, --properties> CategoryName, [CategoryName,...] | PropertyName, [PropertyName,...]]
[B<--PropertiesListing> ByGroup | Alphabetical] [B<-q, --quote> yes | no] [B<-r, --root> rootname]
[B<-w, --WorkingDir> dirname] NucleicAcidIDs...

=head1 DESCRIPTION

List nucleic acid properties. Nucleic acids identification supports two types of IDs: code
or name. Nucleic acid properties data, in addition to basic information about nucleic acids - code,
name, type, chemical formula and molecular weight - include information about exact mass and
elemental composition.

=head1 PARAMETERS

=over 4

=item B<NucleicAcidIDs> I<Code [NucleicAcidName...] | NucleicAcidType [NucleicAcidType...]>

I<NucleicAcidIDs> is a space delimited list of values to identify nucleic acids.

For I<NucleicAcidID> mode, input value format is: I<Code [NucleicAcidName...]>. Default: I<A>.
Examples:

    A
    dG AMP
    Cytidine T UDP dpppA "5'-dATP"

For I<NucleicAcidType> mode, input value format is: I<NucleicAcidType [NucleicAcidType...]>.
Default: I<A>. Possible values are: I<Nucleobase, Nucleoside, Deoxynucleoside, Nucleotide,
Deoxynucleotide>. Default: I<Nucleoside>.
Examples:

    Deoxynucleoside
    Nucleobase Nucleotide

=back

=head1 OPTIONS

=over 4

=item B<-h, --help>

Print this help message.

=item B<-m, --mode> I<NucleicAcidID | NucleicAcidType>

Specify nucleic acids for listing properties using one of these methods: nucleic acid
code and/or names or nucleic acid type.

Possible values: I<NucleicAcidID or NucleicAcidType>. Default: I<NucleicAcidID>

For I<NucleicAcidType>, command line parameters support these type: I<Nucleobase,
Nucleoside, Deoxynucleoside, Nucleotide, Deoxynucleotide>.

=item B<--OutDelim> I<comma | tab | semicolon>

Output text file delimiter. Possible values: I<comma, tab, or semicolon>
Default value: I<comma>.

=item B<--output> I<STDOUT | File>

List information at STDOUT or write it to a file. Possible values: I<STDOUT or File>. Default:
I<STDOUT>. B<-r, --root> option is used to generate output file name.

=item B<--OutputStyle> I<NucleicAcidBlock | NucleicAcidRows>

Specify how to list nucleic acid information: add a new line for each property and present it as a block
for each nucleic acid; or include all properties in one line and show it as a single line.

Possible values: I<NucleicAcidBlock | NucleicAcidRows>. Default: I<NucleicAcidBlock>

An example for I<NucleicAcidBlock> output style:

    Code: Ado
    OtherCodes: A
    Name: Adenosine
    Type: Nucleoside
    MolecularFormula: C10H13O4N5
    MolecularWeight: 267.2413
    ... ...

An example for I<NucleicAcidRows> output style:

    Code,OtherCodes,Name,Type,MolecularFormula,MolecularWeight

=item B<-o, --overwrite>

Overwrite existing files.

=item B<--precision> I<number>

Precision for listing numerical values. Default: up to I<4> decimal places.
Valid values: positive integers.

=item B<--PropertiesMode> I<Categories | Names | All>

Specify how property names are specified: use category names; explicit list of property names; or
use all available properties. Possible values: I<Categories, Names, or All>. Default: I<Categories>.

This option is used in conjunction with B<-p, --properties> option to specify properties of
interest.

=item B<-p, --properties> I<CategoryName,[CategoryName,...] | PropertyName,[PropertyName,...]>

This option is B<--propertiesmode> specific. In general, it's a list of comma separated category or
property names.

Specify which nucleic acid properties information to list for the nucleic acid IDs specified using
command line parameters: list basic information; list all available information; or specify a comma
separated list of nucleic acid property names.

Possible values: I<Basic | BasicPlus | PropertyName,[PropertyName,...]>.
Default: I<Basic>.

I<Basic> includes: I<Code, OtherCodes, Name, Type, MolecularFormula, MolecularWeight>

I<BasicPlus> includes: I<Code, OtherCodes, Name, Type, MolecularFormula, MolecularWeight, ExactMass,
ElementalComposition>

Here is a complete list of available properties: I<Code, OtherCodes, BasePair, Name, Type, MolecularFormula,
MolecularFormulaAtpH7.5, MolecularWeight, ExactMass, ElementalComposition>.

=item B<--PropertiesListing> I<ByGroup | Alphabetical>

Specify how to list properties for nucleic acids: group by category or an alphabetical by
property names. Possible values: I<ByGroup or Alphabetical>. Default: I<ByGroup>

=item B<-q, --quote> I<yes | no>

Put quotes around column values in output text file. Possible values: I<yes or
no>. Default value: I<yes>.

=item B<-r, --root> I<rootname>

New text file name is generated using the root: <Root>.<Ext>. File name is only
used during I<File> value of B<-o, --output> option.

Default file name: NucleicAcidInfo<mode>.<Ext>. The csv, and tsv
<Ext> values are used for comma/semicolon, and tab delimited text files respectively.

=item B<-w, --WorkingDir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To list basic properties information for nucleoside A, type:

    % InfoNucleicAcids.pl

To list all available properties information for nucleoside A, type:

    % InfoNucleicAcids.pl --propertiesmode all A

To list all available information for all available nucleic acids, type:

    % InfoNucleicAcids.pl --propertiesmode All All

To list basic properties information for all nucleobases, type:

    % InfoNucleicAcids.pl -m NucleicAcidType Nucleoside

To list basic properties information for all nucleotides and deoxynulceotides, type:

    % InfoNucleicAcids.pl -m NucleicAcidType Nucleotide Deoxynucleotide

To list basic properties information for variety of nucleic acids, type:

    % InfoNucleicAcids.pl A dG AMP Cytidine T UDP "5'-dATP"

To list code and molecular weights for nucleosides A, G, C and T, type:

    % InfoNucleicAcids.pl --PropertiesMode  Names --properties
      Code,MolecularWeight A G C T

To alphabetically list all the available properties for nucleotides dAMP, dGMP,
dCMP, and dTMP in rows instead of nucleic acid blocks with quotes around the values, type:

    % InfoNucleicAcids.pl --PropertiesMode All --PropertiesListing
      Alphabetical --OutputStyle NucleicAcidRows -q yes dAMP dGMP
      dCMP dTMP

To alphabetically list all the available properties for all available nucleic acids to
a file names NucleicAcidsProperties.csv with quotes around the values, type

    % InfoNucleicAcids.pl --PropertiesMode All --PropertiesListing
      Alphabetical --output File --OutputStyle NucleicAcidRows -r
      NucleicAcidsProperties -o -q Yes All

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoAminoAcids.pl, InfoPeriodicTableElements.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
