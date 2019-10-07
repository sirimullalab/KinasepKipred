#!/usr/bin/perl -w
#
# File: InfoPeriodicTableElements.pl
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
use PeriodicTable;

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

ListElementProperties();
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# List atomic properties for elements...
sub ListElementProperties {
  my($ElementID, $ElementDataRef, $PropertyName, $PropertyValue, $PropertyUnits, $PropertyUnitsRef, @PropertyLabels, @PropertyValues);

  print "Listing information for periodic table element(s)...\n";

  if ($OptionsInfo{FileOutput}) {
    print "Generating file $OptionsInfo{OutFileName}...\n";
    open OUTFILE, ">$OptionsInfo{OutFileName}" or die "Couldn't open $OptionsInfo{OutFileName}: $!\n";
  }

  # Setup property labels...
  @PropertyLabels = ();
  $PropertyUnitsRef = PeriodicTable::GetElementPropertiesNamesAndUnits();
  for $PropertyName (@{$OptionsInfo{SpecifiedProperies}}) {
    $PropertyUnits = (exists $PropertyUnitsRef->{$PropertyName}) ? $PropertyUnitsRef->{$PropertyName} : '';
    if ($PropertyName =~ /^NaturalIsotopeData$/i) {
      push @PropertyLabels, qw(MassNumber: RelativeAtomicMass: NaturalAbundance:);
    }
    else {
      push @PropertyLabels, ($PropertyUnits ? "$PropertyName ($PropertyUnits):" : "$PropertyName:");
    }
  }

  if ($OptionsInfo{ElementRowsOutput}) {
    ListHeaderRowData(\@PropertyLabels);
  }

  # Go over specified properties...
  for $ElementID (@{$OptionsInfo{SpecifiedElementIDs}}) {
    $ElementDataRef = PeriodicTable::GetElementPropertiesData($ElementID);

    if (!$OptionsInfo{ElementRowsOutput}) {
      if ($OptionsInfo{FileOutput}) {
	print OUTFILE "\nListing atomic properties for element $ElementID...\n\n";
      }
      else {
	print "\nListing atomic properties for element $ElementID...\n\n";
      }
    }

    # Collect data..
    @PropertyValues = ();
    for $PropertyName (@{$OptionsInfo{SpecifiedProperies}}) {
      if ($PropertyName =~ /^NaturalIsotopeData$/i) {
	push @PropertyValues, SetupIsotopeData($ElementID);
      }
      else {
	$PropertyValue = $ElementDataRef->{$PropertyName};
	if (IsFloat($PropertyValue)) {
	  $PropertyValue = sprintf("%.$OptionsInfo{Precision}f", $PropertyValue) + 0;
	}
	push @PropertyValues, $PropertyValue;
      }
    }
    # List data...
    ListElementData(\@PropertyLabels, \@PropertyValues);
  }
  if ($OptionsInfo{FileOutput}) {
    close OUTFILE;
  }
  print "\n";
}

# List data for an element...
sub ListElementData {
  my($DataLabelRef, $DataValueRef) = @_;
  my($Index, $Line, $Value);

  if ($OptionsInfo{ElementRowsOutput}) {
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

# List data for an element...
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

# Setup isotope data strings...
sub SetupIsotopeData {
  my($ElementID) = @_;
  my($MassNumber, $RelativeAtomicMass, $NaturalAbundance, $NaturalIsotopeDataRef, @MassNumbers, @RelativeAtomicMasses, @NaturalAbundances);

  # Get natural isotope data: MassNumber, RelativeAtomicMass and NaturalAbundance
  @MassNumbers = (); @RelativeAtomicMasses = (); @NaturalAbundances = ();
  $NaturalIsotopeDataRef = PeriodicTable::GetElementNaturalIsotopesData($ElementID);
  for $MassNumber (sort {$a <=> $b} keys %{$NaturalIsotopeDataRef}) {
    $RelativeAtomicMass = $NaturalIsotopeDataRef->{$MassNumber}{RelativeAtomicMass};
    $NaturalAbundance = $NaturalIsotopeDataRef->{$MassNumber}{NaturalAbundance};
    push @MassNumbers, $MassNumber;
    $RelativeAtomicMass = ($RelativeAtomicMass > 0) ? (sprintf("%.$OptionsInfo{Precision}f", $RelativeAtomicMass) + 0) : '';
    push @RelativeAtomicMasses, $RelativeAtomicMass;
    $NaturalAbundance = ($NaturalAbundance > 0) ? (sprintf("%.$OptionsInfo{Precision}f", $NaturalAbundance) + 0) : '';
    push @NaturalAbundances, $NaturalAbundance;
  }
  $MassNumber = JoinWords(\@MassNumbers, ",", 0);
  $RelativeAtomicMass = JoinWords(\@RelativeAtomicMasses, ",", 0);
  $NaturalAbundance = JoinWords(\@NaturalAbundances, ",", 0);
  return ($MassNumber, $RelativeAtomicMass, $NaturalAbundance);
}

# Get propery names from categories...
sub GetPropertyNamesFromCategories {
  my($CategoryName) = @_;
  my(@PropertyNames);

  @PropertyNames = ();
  if ($CategoryName =~ /^Basic$/i) {
    @PropertyNames = ('AtomicNumber', 'ElementSymbol', 'ElementName', 'AtomicWeight', 'GroundStateConfiguration', 'GroupNumber', 'PeriodNumber', 'FirstIonizationEnergy');
  } elsif ($CategoryName =~ /^BasicAndNaturalIsotope$/i) {
    # Natural isotope data includes: 'MassNumber', 'RelativeAtomicMass', 'NaturalAbundance'
    @PropertyNames = ('AtomicNumber', 'ElementSymbol', 'ElementName', 'AtomicWeight', 'GroundStateConfiguration', 'GroupNumber', 'PeriodNumber', 'FirstIonizationEnergy', 'NaturalIsotopeData');
  } elsif ($CategoryName =~ /^NaturalIsotope$/i) {
    @PropertyNames = ('AtomicNumber', 'ElementSymbol', 'ElementName', 'NaturalIsotopeData');
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

  $OptionsInfo{ElementRowsOutput} = ($Options{outputstyle} =~ /^ElementRows$/i) ? 1 : 0;
  $OptionsInfo{FileOutput} = ($Options{output} =~ /^File$/i) ? 1 : 0;

  $OptionsInfo{Precision} = $Options{precision};

  my($ElementID, @ElementIDs, @GroupElements, @PeriodElements, %GroupNamesMap);

  @{$OptionsInfo{SpecifiedElementIDs}} = ();
  if (@ARGV >=1 && ($Options{mode} =~ /^All$/i) ) {
    warn "Warning: Ignoring comman line element IDs: Not valid for All value of \"-m --mode\" option...\n";
  }

  # Set up element IDs except for All mode...
  @ElementIDs = ();
  %GroupNamesMap = ();

  if (@ARGV >=1 ) {
    if ($Options{mode} !~ /^All$/i) {
      push @ElementIDs, @ARGV;
    }
  }
  else {
    # Setup mode specified default values...
    my($Nothing);
    MODE: {
	if ($Options{mode} =~ /^ElementID$/i) { push @ElementIDs, 'H'; last MODE; };
	if ($Options{mode} =~ /^AmericanGroupLabel$/i) { push @ElementIDs, 'IA'; last MODE; };
	if ($Options{mode} =~ /^EuropeanGroupLabel$/i) { push @ElementIDs, 'IA'; last MODE; };
	if ($Options{mode} =~ /^GroupNumber$/i) { push @ElementIDs, '1'; last MODE; };
	if ($Options{mode} =~ /^GroupName$/i) { push @ElementIDs, 'AlkaliMetals'; last MODE; };
	if ($Options{mode} =~ /^PeriodNumber$/i) { push @ElementIDs, '1'; last MODE; };
	$Nothing = 1;
      }
  }
  if ($Options{mode} =~ /^GroupName$/i) {
    # Map group names to what's stored in Perioidic table data file...
    %GroupNamesMap = ('alkalimetals', 'Alkali metal', 'alkalineearthmetals', 'Alkaline earth metal', 'chalcogens', 'Chalcogen', 'coinagemetals', 'Coinage metal', 'halogens', 'Halogen', 'noblegases', 'Noble gas', 'pnictogens', 'Pnictogen', 'lanthanides', 'Lanthanoid', 'lanthanoids', 'Lanthanoid', 'actinides', 'Actinoid', 'actinoids', 'Actinoid' );
  }

  # Generate list of elements...
  if ($Options{mode} =~ /^All$/i) {
    push @{$OptionsInfo{SpecifiedElementIDs}}, PeriodicTable::GetElements();
  }
  else {
    ELEMENTID: for $ElementID (@ElementIDs) {
      if ($Options{mode} =~ /^ElementID$/i) {
	if (PeriodicTable::IsElement($ElementID)) {
	  push @{$OptionsInfo{SpecifiedElementIDs}}, $ElementID;
	}
	else {
	  warn "Ignoring element ID, $ElementID, specified using command line parameter: Unknown element ID...\n";
	  next ELEMENTID;
	}
      }
      elsif ($Options{mode} =~ /^AmericanGroupLabel$/i) {
	if (@GroupElements = PeriodicTable::GetElementsByAmericanStyleGroupLabel($ElementID)) {
	  push @{$OptionsInfo{SpecifiedElementIDs}}, @GroupElements;
	}
	else {
	  warn "Ignoring American style group label, $ElementID, specified using command line parameter: Unknown group label...\n";
	  next ELEMENTID;
	}
      }
      elsif ($Options{mode} =~ /^EuropeanGroupLabel$/i) {
	if (@GroupElements = PeriodicTable::GetElementsByEuropeanStyleGroupLabel($ElementID)) {
	  push @{$OptionsInfo{SpecifiedElementIDs}}, @GroupElements;
	}
	else {
	  warn "Ignoring American style group label, $ElementID, specified using command line parameter: Unknown group label...\n";
	  next ELEMENTID;
	}
      }
      elsif ($Options{mode} =~ /^GroupNumber$/i) {
	if (@GroupElements = PeriodicTable::GetElementsByGroupNumber($ElementID)) {
	  push @{$OptionsInfo{SpecifiedElementIDs}}, @GroupElements;
	}
	else {
	  warn "Ignoring group number, $ElementID, specified using command line parameter: Unknown group number...\n";
	  next ELEMENTID;
	}
      }
      elsif ($Options{mode} =~ /^GroupName$/i) {
	if (exists $GroupNamesMap{lc($ElementID)}) {
	  @GroupElements = PeriodicTable::GetElementsByGroupName($GroupNamesMap{lc($ElementID)});
	  push @{$OptionsInfo{SpecifiedElementIDs}}, @GroupElements;
	}
	else {
	  warn "Ignoring group name, $ElementID, specified using command line parameter: Unknown group name...\n";
	  next ELEMENTID;
	}
      }
      elsif ($Options{mode} =~ /^PeriodNumber$/i) {
	if (@GroupElements = PeriodicTable::GetElementsByPeriodNumber($ElementID)) {
	  push @{$OptionsInfo{SpecifiedElementIDs}}, @GroupElements;
	}
	else {
	  warn "Ignoring period number, $ElementID, specified using command line parameter: Unknown period number...\n";
	  next ELEMENTID;
	}
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
      $OutFileRoot = 'PeriodicTableElementsInfo' . $Options{mode};
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

  # Make sure atomic appropriate properties/category names are specified...
  @{$OptionsInfo{SpecifiedProperies}} = ();
  if ($Options{properties} && ($Options{propertiesmode} =~ /^All$/i) ) {
    warn "Warning: Ignoring values specifed by \"-p --properties\" option: Not valid for All value of \"--propertiesmode\" option...\n";
  }
  if ($Options{propertiesmode} =~ /^All$/i) {
    if ($Options{propertieslisting} =~ /^Alphabetical$/i) {
      push @{$OptionsInfo{SpecifiedProperies}}, PeriodicTable::GetElementPropertiesNames('Alphabetical');
    }
    else {
      push @{$OptionsInfo{SpecifiedProperies}}, PeriodicTable::GetElementPropertiesNames();
    }
    push @{$OptionsInfo{SpecifiedProperies}}, 'NaturalIsotopeData';
  }
  else {
    if ($Options{properties}) {
      if ($Options{propertiesmode} =~ /^Categories$/i) {
	# Check category name...
	if ($Options{properties} !~ /^(Basic|BasicAndNaturalIsotope|NaturalIsotope)$/i) {
	  die "Error: The value specified, $Options{properties}, for option \"-p --properties\" in conjunction with \"Categories\" value for option \"--propertiesmode\" is not valid. Allowed values: Basic, BasicAndNaturalIsotope, NaturalIsotope\n";
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
	  if ($PropertyName =~ /^NaturalIsotopeData$/i) {
	    push @{$OptionsInfo{SpecifiedProperies}}, $PropertyName;
	    next NAME;
	  }
	  if (PeriodicTable::IsElementProperty($PropertyName)) {
	    push @{$OptionsInfo{SpecifiedProperies}}, $PropertyName;
	  }
	  else {
	    warn "Warning: Ignoring value, $Name, specifed by \"-p --properties\" option: Unknown property name...\n";
	  }
	}
	if ($Options{propertieslisting} =~ /^Alphabetical$/i) {
	  # AtomicNumber, ElementSymbol and ElementName are always listed first and
	  # NaturalIsotopeData in the end...
	  my($AtomicNumberPresent, $ElementSymbolPresent, $ElementNamePresent, $NaturalIsotopeDataPresent, @AlphabeticalProperties, %PropertiesMap);
	  %PropertiesMap = ();
	  @AlphabeticalProperties = ();
	  $AtomicNumberPresent = 0; $ElementSymbolPresent = 0; $ElementNamePresent = 0; $NaturalIsotopeDataPresent = 0;
	  NAME: for $Name (@{$OptionsInfo{SpecifiedProperies}}) {
	    if ($Name =~ /^AtomicNumber$/i) {
	      $AtomicNumberPresent = 1;
	      next NAME;
	    }
	    if ($Name =~ /^ElementSymbol$/i) {
	      $ElementSymbolPresent = 1;
	      next NAME;
	    }
	    if ($Name =~ /^ElementName$/i) {
	      $ElementNamePresent = 1;
	      next NAME;
	    }
	    if ($Name =~ /^NaturalIsotopeData$/i) {
	      $NaturalIsotopeDataPresent = 1;
	      next NAME;
	    }
	    $PropertiesMap{$Name} = $Name;
	  }
	  # Setup the alphabetical list...
	  if ($AtomicNumberPresent) {
	    push @AlphabeticalProperties, 'AtomicNumber';
	  }
	  if ($ElementSymbolPresent) {
	    push @AlphabeticalProperties, 'ElementSymbol';
	  }
	  if ($ElementNamePresent) {
	    push @AlphabeticalProperties, 'ElementName';
	  }
	  for $Name (sort keys %PropertiesMap) {
	    push @AlphabeticalProperties, $Name;
	  }
	  if ($NaturalIsotopeDataPresent) {
	    push @AlphabeticalProperties, 'NaturalIsotopeData';
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
  $Options{mode} = "ElementID";
  $Options{outdelim} = "comma";
  $Options{output} = "STDOUT";
  $Options{outputstyle} = "ElementBlock";
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
  if ($Options{mode} !~ /^(ElementID|AmericanGroupLabel|EuropeanGroupLabel|GroupNumber|GroupName|PeriodNumber|All)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: ElementID, AmericanGroupLabel, EuropeanGroupLabel, GroupNumber, GroupName, PeriodNumber, or All\n";
  }
  if ($Options{outdelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{outdelim}, for option \"--outdelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if ($Options{output} !~ /^(STDOUT|File)$/i) {
    die "Error: The value specified, $Options{output}, for option \"--output\" is not valid. Allowed values: STDOUT or File\n";
  }
  if ($Options{outputstyle} !~ /^(ElementBlock|ElementRows)$/i) {
    die "Error: The value specified, $Options{outputstyle}, for option \"--outputstyle\" is not valid. Allowed values: ElementBlock or ElementRows\n";
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

InfoPeriodicTableElements.pl - List atomic properties of elements

=head1 SYNOPSIS

InfoPeriodicTableElements.pl ElementID(s)...

InfoPeriodicTableElements.pl [B<-h, --help>]
[B<-m, --mode> ElementID | AmericanGroupLabel | EuropeanGroupLabel | GroupNumber | GroupName | PeriodNumber | All]
[B<--outdelim> comma | tab | semicolon] [B<--output> STDOUT | File] [B<--outputstyle> ElementBlock | ElementRows]
[B<-o, --overwrite>] [B<--precision> number] [B<--propertiesmode> Categories | Names | All]
[B<-p, --properties> CategoryName,[CategoryName,...] | PropertyName,[PropertyName,...]]
[B<--propertieslinting> ByGroup | Alphabetical] [B<-q, --quote> yes | no] [B<-r, --root> rootname]
[B<-w, --workingdir> dirname] ElementID(s)...

=head1 DESCRIPTION

List atomic properties of elements in the periodic table. A variety of methods are available to
specify elements of interest: atomic numbers, element symbols, American or European style group
labels, IUPAC group numbers, period numbers, and group names.

Atomic properties data, in addition to basic information about the periodic table elements, is
also available for these categories: atomic radii, bulk properties, common valences, electronegativities,
electron affinities, historical data, ionization energies, natural isotopes, oxidation states,
and thermal properties.

Natural isotopes data include mass number, relative atomic mass and percent natural
abundance for each isotope of an element.

=head1 PARAMETERS

=over 4

=item B<ElementIDs> I<ElementSymbol [AtomicNumber...] | GroupLabel [GroupLabel...] | GroupNumbel [GroupNumber...] | PeriodNumber [PeriodNumbe...]>

Command line specification of elements is mode specific. In general, it's a space delimited list of values to identify
elements. All element IDs must correspond to a specific mode; mixed specifications is not supported.

For I<ElementID> mode, input value format is: I<AtomicNumber [ElementSymbol ...]>. Default: I<H>.
Examples:

    C
    6
    C N O P S Cl
    6 7 8 15 16 17
    C 7 8 15 S 17

For I<AmericanGroupLabel> mode, input value format is: I<GroupLabel [GroupLabel ...]>. Default: I<IA>. Possible
group label values are: I<IA IIA IIIB IVB VB VIB VIIB VIII or VIIIB IB IIB IIIA IVA VA,
VIA, VIIA, VIIA>. Examples:

    IA
    IA IVA IIB

For I<EuropeanGroupLabel> mode, input value format is: I<GroupLabel [GroupLabel ...]>. Default: I<IA>. Possible
group label values are: I<IA IIA IIIA IVA VA VIA VIIA VIII or VIIIA IB IIB IIIB IVB VB,
VIB VIIB VIIB>. Examples:

    IA
    IA IVB IIB

For IUPAC I<GroupNumber> mode, input value format is: I<GroupNumber [GroupNumber...]>. Default: I<1>. Possible
group label values are: I<1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18>. Examples:

    1
    1 14 12

For  I<GroupName> mode, input value format is: I<GroupName [GroupName...]>. Default: I<AlkaliMetals>. Possible
group name values are: I<AlkaliMetals AlkalineEarthMetals Chalcogens CoinageMetals Halogens
NobleGases Pnictogens Lanthanides or Lanthanoids, Actinides or Actinoids>. Examples:

    AlkaliMetals
    AlkaliMetals Halogens NobleGases

For I<PeriodNumber> mode, input value format is: I<PeriodNumber [PeriodNumber,...]>. Default: I<1>. Possible
group label values are: I<1 2 3 4 5 6 7>. Examples:

    1
    1 2 3

For I<All> mode, no input value is needed and atomic properties information is listed for all the
elements.

=back

=head1 OPTIONS

=over 4

=item B<-h, --help>

Print this help message.

=item B<-m, --mode> I<ElementID | AmericanGroupLabel | EuropeanGroupLabel | GroupNumber | GroupName | PeriodNumber | All>

Specify elements for listing atomic properties using one of these methods: atomic numbers
and/or element symbols list, American style group labels, European style group labels, IUPAC
group number, group names, period numbers, or all elements.

Possible values: I<ElementID, AmericanGroupLabel, EuropeanGroupLabel, GroupNumber,
GroupName, PeriodNumber, All>. Default: I<ElementID>.

=item B<--outdelim> I<comma | tab | semicolon>

Output text file delimiter. Possible values: I<comma, tab, or semicolon>
Default value: I<comma>.

=item B<--output> I<STDOUT | File>

List information at STDOUT or write it to a file. Possible values: I<STDOUT or File>. Default:
I<STDOUT>. B<-r, --root> option is used to generate output file name.

=item B<--outputstyle> I<ElementBlock | ElementRows>

Specify how to list element information: add a new line for each property and present it as a block
for each element; or include all properties in one line and show it as a single line.

Possible values: I<ElementBlock | ElementRows>. Default: I<ElementBlock>

An example for I<ElementBlock> output style:

    Atomic number: 1
    Element symbol: H
    Element name: Hydrogen
    Atomic weight: 1.00794
    ... ...
    ... ...

    Atomic number: 6
    Element symbol: C
    Element name: Carbon
    Atomic weight: 12.0107
    ... ...
    ... ...

An example for I<ElementRows> output style:

    Atomic number, Element symbol, Elemenet name, Atomic weight, ...
    1,H,Hydrogen,1.00794,..
    6,C,Carbon,12.0107,..

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

Specify which atomic properties information to list for the elements specified using command line
parameters: list basic and/or isotope information; list all available information; or specify a comma
separated list of atomic property names.

Possible values: I<Basic| BasicAndNaturalIsotope | NaturalIsotope | PropertyName,[PropertyName,...]>.
Default: I<Basic>.

I<Basic> includes: I<AtomicNumber, ElementSymbol, ElementName, AtomicWeight, GroundStateConfiguration,
GroupNumber, PeriodNumber, FirstIonizationEnergy>.

I<NaturalIsotope> includes: I<AtomicNumber, ElementSymbol, ElementName, MassNumber,
RelativeAtomicMass, NaturalAbundance>.

Here is a complete list of available properties: AllenElectronegativity, AllredRochowElectronegativity, AtomicNumber,
AtomicRadiusCalculated, AtomicRadiusEmpirical, AtomicWeight, Block, BoilingPoint, BondLength,
BrinellHardness, BulkModulus, Classification, CoefficientOfLinearExpansion, Color,
CommonValences, LowestCommonValence, HighestCommonValence,
CommonOxidationNumbers, LowestCommonOxidationNumber, HighestCommonOxidationNumber,
CovalentRadiusEmpirical, CriticalTemperature, DensityOfSolid, DiscoveredAt, DiscoveredBy,
DiscoveredWhen, ElectricalResistivity, ElectronAffinity, ElementName, ElementSymbol, EnthalpyOfAtmization,
EnthalpyOfFusion, EnthalpyOfVaporization, FirstIonizationEnergy, GroundStateConfiguration, GroundStateLevel,
GroupName, GroupNumber, NaturalIsotopeData, MeltingPoint, MineralHardness, MolarVolume,
MullikenJaffeElectronegativity, OriginOfName, PaulingElectronegativity, PeriodNumber, PoissonsRatio,
Reflectivity, RefractiveIndex, RigidityModulus, SandersonElectronegativity, StandardState,
SuperconductionTemperature, ThermalConductivity, VanderWaalsRadius, VelocityOfSound, VickersHardness,
YoungsModulus.

=item B<--propertieslisting> I<ByGroup | Alphabetical>

Specify how to list properties for elements: group by category or an alphabetical by
property names. Possible values: I<ByGroup or Alphabetical>. Default: I<ByGroup>.
During I<Alphabetical> listing, element identification data - I<AtomicNumber, ElementSymbol,
ElementName> - is show first,  and natural isotope data - I<MassNumber, RelativeAtomicMass,
NaturalAbundance> - is listed in the end.

=item B<-q, --quote> I<yes | no>

Put quotes around column values in output text file. Possible values: I<yes or
no>. Default value: I<yes>.

=item B<-r, --root> I<rootname>

New text file name is generated using the root: <Root>.<Ext>. File name is only
used during I<File> value of B<-o, --output> option.

Default file name: PeriodicTableElementsInfo<mode>.<Ext>. The csv, and tsv
<Ext> values are used for comma/semicolon, and tab delimited text files respectively.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To list basic atomic properties information for element H, type:

    % InfoPeriodicTableElements.pl

To list basic atomic properties information for elements C,N,O and F, type:

    % InfoPeriodicTableElements.pl C N O F

To list all available atomic properties information for elements C,N,O and F, type:

    % InfoPeriodicTableElements.pl --propertiesmode all 6 N O 9

To list basic and natural isotope information for elements C,N,O and F, type:

    % InfoPeriodicTableElements.pl --propertiesmode Categories
      --properties BasicAndNaturalIsotope  C N O F

To list AtomicNumber, ElementName, AtomicWeight and CommonValences information
for elements C,N,O and F, type:

    % InfoPeriodicTableElements.pl --propertiesmode Names
      --properties AtomicNumber,ElementName,AtomicWeight,CommonValences
      C N O F

To alphabetically list basic and natural isotope information for elements C,N,O and F in rows instead of
element blocks with quotes around the values, type:

    % InfoPeriodicTableElements.pl --propertiesmode Categories
      --properties BasicAndNaturalIsotope --propertieslisting alphabetical
      --outdelim comma --outputstyle ElementRows --quote yes C N O F

To alphabetically list all available atomic information for elements C,N,O and F in rows instead of
element blocks with quotes around the values and write them into a file ElementProperties.csv, type:

    % InfoPeriodicTableElements.pl --propertiesmode Categories
      --properties BasicAndNaturalIsotope --propertieslisting alphabetical
      --outdelim comma --outputstyle ElementRows --quote yes
      --output File -r ElementsProperties -o -m All

To list basic atomic properties information for elements in groups IA and VIA using American
style group labels, type:

    % InfoPeriodicTableElements.pl -m AmericanGroupLabel IA VIA

To list basic atomic properties information for elements in groups IA and VB using European
style group labels, type:

    % InfoPeriodicTableElements.pl -m AmericanGroupLabel IA VB

To list basic atomic properties information for elements in groups Halogens and NobleGases, type:

    % InfoPeriodicTableElements.pl -m GroupName Halogens NobleGases

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoAminoAcids.pl InfoNucleicAcids.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
