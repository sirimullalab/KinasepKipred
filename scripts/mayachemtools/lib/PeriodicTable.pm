package PeriodicTable;
#
# File: PeriodicTable.pm
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
use Carp;
use Text::ParseWords;
use TextUtil;
use FileUtil;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GetElements GetElementsByGroupName GetElementsByGroupNumber GetElementsByAmericanStyleGroupLabel GetElementsByEuropeanStyleGroupLabel GetElementsByPeriodNumber GetElementMostAbundantNaturalIsotopeData GetElementNaturalIsotopeCount GetElementNaturalIsotopesData GetElementNaturalIsotopeAbundance GetElementMostAbundantNaturalIsotopeMass GetElementMostAbundantNaturalIsotopeMassNumber GetElementNaturalIsotopeMass GetElementPropertiesData GetElementPropertiesNames GetElementPropertiesNamesAndUnits GetElementPropertyUnits GetIUPACGroupNumberFromAmericanStyleGroupLabel GetIUPACGroupNumberFromEuropeanStyleGroupLabel IsElement IsElementNaturalIsotopeMassNumber IsElementProperty);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

#
# Load atomic properties and isotope data for elements...
#
my(%ElementDataMap, %ElementIsotopeDataMap, %ElementSymbolMap, @ElementPropertyNames, %ElementPropertyNamesMap, %ElementIsotopeDerivedDataMap);
_LoadPeriodicTableElementData();

#
# Get a list of all known element symbols...
#
sub GetElements {
  my($AtomicNumber, @ElementSymbols);

  @ElementSymbols = ();
  for $AtomicNumber (sort {$a <=> $b} keys %ElementDataMap) {
      push @ElementSymbols, $ElementDataMap{$AtomicNumber}{ElementSymbol};
  }
  return (wantarray ? @ElementSymbols : \@ElementSymbols);
}

#
# Get element symbols of elements with a specific group name. Valid group
# names are: Alkali metals, Alkaline earth metals, Coinage metals, Pnictogens,
# Chalcogens, Halogens, Noble gases; Additionally, usage of Lanthanides (Lanthanoids)
# and Actinides (Actinoids) is also supported.
#
sub GetElementsByGroupName {
  my($SpecifiedGroupName) = @_;
  my($AtomicNumber, @ElementSymbols, $GroupName);

  if (IsEmpty($SpecifiedGroupName)) {
    return (wantarray ? () : undef);
  }
  if ($SpecifiedGroupName =~ /Lanthanide/i) {
    $SpecifiedGroupName = 'Lanthanoids';
  }
  elsif ($SpecifiedGroupName =~ /Actinide/i) {
    $SpecifiedGroupName = 'Actinoids';
  }
  @ElementSymbols = ();
  for $AtomicNumber (sort {$a <=> $b} keys %ElementDataMap) {
    $GroupName = $ElementDataMap{$AtomicNumber}{GroupName};
    if ($SpecifiedGroupName =~ /$GroupName/i ) {
      push @ElementSymbols, $ElementDataMap{$AtomicNumber}{ElementSymbol};
    }
  }
  return (wantarray ? @ElementSymbols : \@ElementSymbols);
}

#
# Get element symbols of elements in a specific IUPAC group number.
# A reference to an array containing element symbols is returned.
#
sub GetElementsByGroupNumber {
  my($GroupNumber) = @_;
  my($AtomicNumber, @ElementSymbols);

  if (!IsInteger($GroupNumber)) {
    return (wantarray ? () : undef);
  }

  @ElementSymbols = ();
  for $AtomicNumber (sort {$a <=> $b} keys %ElementDataMap) {
    if ($GroupNumber eq $ElementDataMap{$AtomicNumber}{GroupNumber}) {
      push @ElementSymbols, $ElementDataMap{$AtomicNumber}{ElementSymbol};
    }
  }
  return (wantarray ? @ElementSymbols : \@ElementSymbols);
}

#
# Get element symbols of elements in a specific American style group label.
# A reference to an array containing element symbols is returned.
#
sub GetElementsByAmericanStyleGroupLabel {
  my($GroupLabel) = @_;

  return _GetElementsByGroupLabel($GroupLabel, 'AmericanStyle');
}

#
# Get element symbols of elements in a specific European style group label.
# A reference to an array containing element symbols is returned.
#
sub GetElementsByEuropeanStyleGroupLabel {
  my($GroupLabel) = @_;

  return _GetElementsByGroupLabel($GroupLabel, 'EuropeanStyle');
}

#
# Get IUPAC group number from American style group label. A comma delimited
# string is returned for group VIII or VIIIB.
#
sub GetIUPACGroupNumberFromAmericanStyleGroupLabel {
  my($GroupLabel) = @_;
  my($GroupNumber);

  if (IsEmpty($GroupLabel)) {
    return undef;
  }
  $GroupNumber = "";
  SWITCH: {
      if ($GroupLabel =~ /^IA$/) { $GroupNumber = 1; last SWITCH;}
      if ($GroupLabel =~ /^IIA$/) { $GroupNumber = 2; last SWITCH;}
      if ($GroupLabel =~ /^IIIB$/) { $GroupNumber = 3; last SWITCH;}
      if ($GroupLabel =~ /^IVB$/) { $GroupNumber = 4; last SWITCH;}
      if ($GroupLabel =~ /^VB$/) { $GroupNumber = 5; last SWITCH;}
      if ($GroupLabel =~ /^VIB$/) { $GroupNumber = 6; last SWITCH;}
      if ($GroupLabel =~ /^VIIB$/) { $GroupNumber = 7; last SWITCH;}
      if ($GroupLabel =~ /^(VIII|VIIIB)$/) { $GroupNumber = '8,9,10'; last SWITCH;}
      if ($GroupLabel =~ /^IB$/) { $GroupNumber = 11; last SWITCH;}
      if ($GroupLabel =~ /^IIB$/) { $GroupNumber = 12; last SWITCH;}
      if ($GroupLabel =~ /^IIIA$/) { $GroupNumber = 13; last SWITCH;}
      if ($GroupLabel =~ /^IVA$/) { $GroupNumber = 14; last SWITCH;}
      if ($GroupLabel =~ /^VA$/) { $GroupNumber = 15; last SWITCH;}
      if ($GroupLabel =~ /^VIA$/) { $GroupNumber = 16; last SWITCH;}
      if ($GroupLabel =~ /^VIIA$/) { $GroupNumber = 17; last SWITCH;}
      if ($GroupLabel =~ /^VIIIA$/) { $GroupNumber = 18; last SWITCH;}
      $GroupNumber = "";
  }
  if (!$GroupNumber) {
    return undef;
  }
  return $GroupNumber;
}

#
# Get IUPAC group number from European style group label. A comma delimited
# string is returned for group VIII or VIIIA
#
sub GetIUPACGroupNumberFromEuropeanStyleGroupLabel {
  my($GroupLabel) = @_;
  my($GroupNumber);

  if (IsEmpty($GroupLabel)) {
    return undef;
  }
  $GroupNumber = "";
  SWITCH: {
      if ($GroupLabel =~ /^IA$/) { $GroupNumber = 1; last SWITCH;}
      if ($GroupLabel =~ /^IIA$/) { $GroupNumber = 2; last SWITCH;}
      if ($GroupLabel =~ /^IIIA$/) { $GroupNumber = 3; last SWITCH;}
      if ($GroupLabel =~ /^IVA$/) { $GroupNumber = 4; last SWITCH;}
      if ($GroupLabel =~ /^VA$/) { $GroupNumber = 5; last SWITCH;}
      if ($GroupLabel =~ /^VIA$/) { $GroupNumber = 6; last SWITCH;}
      if ($GroupLabel =~ /^VIIA$/) { $GroupNumber = 7; last SWITCH;}
      if ($GroupLabel =~ /^(VIII|VIIIA)$/) { $GroupNumber = '8,9,10'; last SWITCH;}
      if ($GroupLabel =~ /^IB$/) { $GroupNumber = 11; last SWITCH;}
      if ($GroupLabel =~ /^IIB$/) { $GroupNumber = 12; last SWITCH;}
      if ($GroupLabel =~ /^IIIB$/) { $GroupNumber = 13; last SWITCH;}
      if ($GroupLabel =~ /^IVB$/) { $GroupNumber = 14; last SWITCH;}
      if ($GroupLabel =~ /^VB$/) { $GroupNumber = 15; last SWITCH;}
      if ($GroupLabel =~ /^VIB$/) { $GroupNumber = 16; last SWITCH;}
      if ($GroupLabel =~ /^VIIB$/) { $GroupNumber = 17; last SWITCH;}
      if ($GroupLabel =~ /^VIIIB$/) { $GroupNumber = 18; last SWITCH;}
      $GroupNumber = "";
  }
  if (!$GroupNumber) {
    return undef;
  }
  return $GroupNumber;
}

#
# Get element symbols of elements in a specific period number.
# A reference to an array containing element symbols is returned.
#
sub GetElementsByPeriodNumber {
  my($SpecifiedPeriodNumber) = @_;
  my($AtomicNumber, $PeriodNumber, @ElementSymbols);

  if (!IsInteger($SpecifiedPeriodNumber)) {
    return (wantarray ? () : undef);
  }

  @ElementSymbols = ();
  for $AtomicNumber (sort {$a <=> $b} keys %ElementDataMap) {
    $PeriodNumber = $ElementDataMap{$AtomicNumber}{PeriodNumber};
    if ($PeriodNumber =~ /\(/) {
      # Lanthanides and Actinides...
      ($PeriodNumber) = split /\(/, $PeriodNumber;
    }
    if ($PeriodNumber == $SpecifiedPeriodNumber) {
      push @ElementSymbols, $ElementDataMap{$AtomicNumber}{ElementSymbol};
    }
  }
  return (wantarray ? @ElementSymbols : \@ElementSymbols);
}

#
# Get data for most abundant isotope of an element using element symbol or atomic number.
#
sub GetElementMostAbundantNaturalIsotopeData {
  my($ElementID) = @_;
  my($AtomicNumber);

  if (!($AtomicNumber = _ValidateElementID($ElementID))) {
    return (wantarray ? () : undef);
  }

  my(@IsotopeData, $IsotopeSymbol, $MassNumber, $RelativeAtomicMass, $NaturalAbundance);
  $MassNumber = $ElementIsotopeDerivedDataMap{$AtomicNumber}{MostAbundantMassNumber};
  $IsotopeSymbol = $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}{IsotopeSymbol};
  $RelativeAtomicMass = $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}{RelativeAtomicMass};
  $NaturalAbundance = $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}{NaturalAbundance};
  @IsotopeData = ();
  @IsotopeData = ($AtomicNumber, $IsotopeSymbol, $MassNumber, $RelativeAtomicMass, $NaturalAbundance);

  return (wantarray ? @IsotopeData : \@IsotopeData);

}
#
# Get natural isotope count for an element...
#
sub GetElementNaturalIsotopeCount {
  my($ElementID) = @_;
  my($AtomicNumber);

  if ($AtomicNumber = _ValidateElementID($ElementID)) {
    return $ElementIsotopeDerivedDataMap{$AtomicNumber}{IsotopeCount};
  }
  else {
    return undef;
  }
}

#
# Get all available isotope data for an element using element symbol or atomic number or
# data for a specific mass number using one of these two invocation methods:
#
# $HashRef = GetElementNaturalIsotopesData($ElementID);
#
# $HashRef = GetElementNaturalIsotopesData($ElementID, $MassNumber);
#
# In the first mode, a reference to a two-dimensional hash array is return where first
# and second dimension keys correspond to mass number and isotope data labels.
#
# And in the second mode, a refernce to one-dimensional hash array is returned with
# keys and values corresponding to isotope data label and values.
#
sub GetElementNaturalIsotopesData {
  my($ElementID, $MassNumber, $InvocationMode, $AtomicNumber);

  if (@_ == 2) {
    ($ElementID, $MassNumber) = @_;
    $InvocationMode = 2;
  }
  else {
    ($ElementID) = @_;
    $InvocationMode = 1;
  }
  if (!($AtomicNumber = _ValidateElementID($ElementID))) {
    return undef;
  }
  if ($InvocationMode == 1) {
    return \%{$ElementIsotopeDataMap{$AtomicNumber}};
  }
  elsif ($InvocationMode == 2) {
    if (exists $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}) {
      return \%{$ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}};
    }
    else {
      return undef;
    }
  }
  else {
    return undef;
  }
}

#
# Get relative atomic mass for an element with specfic mass number.
#
sub GetElementNaturalIsotopeMass {
  my($ElementID, $MassNumber) = @_;
  my($AtomicNumber);

  if (!($AtomicNumber = _ValidateElementID($ElementID))) {
    return undef;
  }
  if (exists $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}) {
    return $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}{RelativeAtomicMass};
  }
  else {
    return undef;
  }
}

#
# Get relative atomic mass of most abundant isotope for an element...
#
sub GetElementMostAbundantNaturalIsotopeMass {
  my($ElementID) = @_;
  my($AtomicNumber);

  if (!($AtomicNumber = _ValidateElementID($ElementID))) {
    return undef;
  }
  my($MassNumber, $RelativeAtomicMass);

  $MassNumber = $ElementIsotopeDerivedDataMap{$AtomicNumber}{MostAbundantMassNumber};
  $RelativeAtomicMass = $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}{RelativeAtomicMass};

  return $RelativeAtomicMass;
}

#
# Get mass number of most abundant isotope for an element...
#
sub GetElementMostAbundantNaturalIsotopeMassNumber {
  my($ElementID) = @_;
  my($AtomicNumber);

  if (!($AtomicNumber = _ValidateElementID($ElementID))) {
    return undef;
  }
  my($MassNumber);

  $MassNumber = $ElementIsotopeDerivedDataMap{$AtomicNumber}{MostAbundantMassNumber};

  return $MassNumber;
}
#
# Get % natural abundance of natural isotope for an element with specfic mass number.
#
sub GetElementNaturalIsotopeAbundance {
  my($ElementID, $MassNumber) = @_;
  my($AtomicNumber);

  if (!($AtomicNumber = _ValidateElementID($ElementID))) {
    return undef;
  }
  if (exists $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}) {
    return $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}{NaturalAbundance};
  }
  else {
    return undef;
  }
}

#
# Get all available properties data for an element using element symbol or atomic number.
# A reference to a hash array is returned with keys and values representing property
# name and its values respectively.
#
sub GetElementPropertiesData {
  my($ElementID) = @_;
  my($AtomicNumber);

  if ($AtomicNumber = _ValidateElementID($ElementID)) {
    return \%{$ElementDataMap{$AtomicNumber}};
  }
  else {
    return undef;
  }
}

#
# Get names of all available element properties. A reference to  an array containing
# names of all available properties is returned.
#
sub GetElementPropertiesNames {
  my($Mode);
  my($PropertyName, @PropertyNames);

  $Mode = 'ByGroup';
  if (@_ == 1) {
    ($Mode) = @_;
  }

  @PropertyNames = ();
  if ($Mode =~ /^Alphabetical$/i) {
    # AtomicNumber, ElementSymbol and ElementName are always listed first...
    push @PropertyNames, qw(AtomicNumber ElementSymbol ElementName);
    for $PropertyName (sort keys %ElementPropertyNamesMap) {
      if ($PropertyName !~ /^(AtomicNumber|ElementSymbol|ElementName)$/i) {
	push @PropertyNames, $PropertyName;
      }
    }
  }
  else {
    push @PropertyNames, @ElementPropertyNames;
  }
  return (wantarray ? @PropertyNames : \@PropertyNames);
}

#
# Get names and units of all available element properties...
# A reference to a hash array is returned with keys and values representing property
# name and its units respectively. Names with no units contains empty strings as hash
# values.
#
sub GetElementPropertiesNamesAndUnits {

   return \%ElementPropertyNamesMap;
}

#
# Get units for a specific element property. An empty string is returned for a property
# with no units.
#
sub GetElementPropertyUnits {
  my($PropertyName) = @_;
  my($PropertyUnits);

  $PropertyUnits = (exists($ElementPropertyNamesMap{$PropertyName})) ? $ElementPropertyNamesMap{$PropertyName} : undef;

  return $PropertyUnits;
}

#
# Is it a known element? Input is either an element symol or a atomic number.
#
sub IsElement {
  my($ElementID) = @_;
  my($Status);

  $Status = (_ValidateElementID($ElementID)) ? 1 : 0;

  return $Status;
}

#
# Is it a valid mass number for an element? Element ID is either an element symol or a atomic number.
#
sub IsElementNaturalIsotopeMassNumber {
  my($ElementID, $MassNumber) = @_;
  my($AtomicNumber, $Status);

  $Status = 0;
  if (!($AtomicNumber = _ValidateElementID($ElementID))) {
    return $Status;
  }
  if (exists $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}) {
    $Status = 1;
  }

  return $Status;
}

#
# Is it an available element property?
#
sub IsElementProperty {
  my($PropertyName) = @_;
  my($Status);

  $Status = (exists($ElementPropertyNamesMap{$PropertyName})) ? 1 : 0;

  return $Status;
}

#
# Implents GetElement<PropertyName> for a valid proprty name.
#
sub AUTOLOAD {
  my($ElementID) = @_;
  my($FunctionName, $PropertyName, $PropertyValue, $AtomicNumber);

  $PropertyValue = undef;

  use vars qw($AUTOLOAD);
  $FunctionName = $AUTOLOAD;
  $FunctionName =~ s/.*:://;

  # Only Get<PropertyName> functions are supported...
  if ($FunctionName !~ /^GetElement/) {
    croak "Error: Function, PeriodicTable::$FunctionName, is not supported by AUTOLOAD in PeriodicTable module: Only GetElement<PropertyName> functions are implemented...";
  }

  $PropertyName = $FunctionName;
  $PropertyName =~  s/^GetElement//;
  if (!exists $ElementPropertyNamesMap{$PropertyName}) {
    croak "Error: Function, PeriodicTable::$FunctionName, is not supported by AUTOLOAD in PeriodicTable module: Unknown element property name, $PropertyName, specified...";
  }

  if (!($AtomicNumber = _ValidateElementID($ElementID))) {
    return undef;
  }
  $PropertyValue = $ElementDataMap{$AtomicNumber}{$PropertyName};
  return $PropertyValue;
}

#
# Get elements labels for group name specified using American or European style...
#
sub _GetElementsByGroupLabel {
  my($GroupLabel, $LabelStyle) = @_;
  my($GroupNumber);

  if ($LabelStyle =~ /^AmericanStyle$/i) {
    $GroupNumber = GetIUPACGroupNumberFromAmericanStyleGroupLabel($GroupLabel);
  }
  elsif ($LabelStyle =~ /^EuropeanStyle$/i) {
    $GroupNumber = GetIUPACGroupNumberFromEuropeanStyleGroupLabel($GroupLabel);
  }

  if (IsEmpty($GroupNumber)) {
    return (wantarray ? () : undef);
  }

  my($AtomicNumber, @GroupElements, @ElementSymbols);
  @ElementSymbols = ();
  if ($GroupNumber =~ /\,/) {
    my(@GroupNumbers);

    @GroupNumbers = split /\,/, $GroupNumber;
    for $GroupNumber (@GroupNumbers) {
      @GroupElements =  GetElementsByGroupNumber($GroupNumber);
      push @ElementSymbols, @GroupElements;
    }
  }
  else {
    @GroupElements =  GetElementsByGroupNumber($GroupNumber);
    push @ElementSymbols, @GroupElements;
  }
  return (wantarray ? @ElementSymbols : \@ElementSymbols);
}

#
# Load PeriodicTableElementData.csv and PeriodicTableIsotopeData.csv files from
# <MayaChemTools>/lib directory...
#
sub _LoadPeriodicTableElementData {
  my($ElementDataFile, $ElementIsotopeDataFile, $MayaChemToolsLibDir);

  $MayaChemToolsLibDir = GetMayaChemToolsLibDirName();

  $ElementDataFile =  "$MayaChemToolsLibDir" . "/data/PeriodicTableElementData.csv";
  $ElementIsotopeDataFile = "$MayaChemToolsLibDir" . "/data/PeriodicTableIsotopeData.csv";

  if (! -e "$ElementDataFile") {
    croak "Error: MayaChemTools package file, $ElementDataFile, is missing: Possible installation problems...";
  }
  if (! -e "$ElementIsotopeDataFile") {
    croak "Error: MayaChemTools package file, $ElementIsotopeDataFile, is missing: Possible installation problems...";
  }

  _LoadElementData($ElementDataFile);
  _LoadElementIsotopeData($ElementIsotopeDataFile);
}

#
# Load PeriodicTableElementData.csv file from <MayaChemTools>/lib directory...
#
sub _LoadElementData {
  my($ElementDataFile) = @_;

  %ElementDataMap = ();
  @ElementPropertyNames = ();
  %ElementPropertyNamesMap = ();
  %ElementSymbolMap = ();

  # Load atomic properties data for all elements...
  #
  # File Format:
  #"AtomicNumber","ElementSymbol","ElementName","AtomicWeight","GroupNumber","GroupName","PeriodNumber","Block","GroundStateConfiguration","ValenceElectrons","GroundStateLevel","StandardState","CommonValences","LowestCommonValence","HighestCommonValence","CommonOxidationNumbers","LowestCommonOxidationNumber","HighestCommonOxidationNumber","BondLength(pm)","AtomicRadiusEmpirical(pm)","AtomicRadiusCalculated(pm)","CovalentRadiusEmpirical(pm)","VanderWaalsRadius(pm)","ElectronAffinity(kJ mol-1)","FirstIonizationEnergy(kJ mol-1)","PaulingElectronegativity(Pauling units)","SandersonElectronegativity(Pauling units)","AllredRochowElectronegativity(Pauling units)","MullikenJaffeElectronegativity(Pauling units)","AllenElectronegativity(Pauling units)","DensityOfSolid(kg m-3)","MolarVolume(cm3)","VelocityOfSound(m s-1)","YoungsModulus(GPa)","RigidityModulus(GPa)","BulkModulus(GPa)","PoissonsRatio(No units)","MineralHardness(No units)","BrinellHardness(MN m-2)","VickersHardness(MN m-2)","ElectricalResistivity(10-8 omega m)","Reflectivity(%)","RefractiveIndex(No units)","MeltingPoint(Celsius)","BoilingPoint(Celsius)","CriticalTemperature(Celsius)","SuperconductionTemperature(Celsius)","ThermalConductivity(W m-1 K-1)","CoefficientOfLinearExpansion(K-1 x 10^6)","EnthalpyOfFusion(kJ mol-1)","EnthalpyOfVaporization(kJ mol-1)","EnthalpyOfAtmization(kJ mol-1)","Color","Classification","DiscoveredBy","DiscoveredAt","DiscoveredWhen","OriginOfName"
  #
  #
  my($AtomicNumber, $ElementSymbol, $Line, $NumOfCols, $InDelim, $Index, $Name, $Value, $Units, @LineWords, @ColLabels);

  $InDelim = "\,";
  open ELEMENTDATAFILE, "$ElementDataFile" or croak "Couldn't open $ElementDataFile: $! ...";

  # Skip lines up to column labels...
  LINE: while ($Line = GetTextLine(\*ELEMENTDATAFILE)) {
    if ($Line !~ /^#/) {
      last LINE;
    }
  }
  @ColLabels= quotewords($InDelim, 0, $Line);
  $NumOfCols = @ColLabels;

  # Extract property names from column labels - and unit names where appropriate...
  @ElementPropertyNames = ();
  for $Index (0 .. $#ColLabels) {
    $Name = $ColLabels[$Index];
    $Units = "";
    if ($Name =~ /\(/) {
      ($Name, $Units) = split /\(/,  $Name;
      $Units =~ s/\)//g;
    }
    push @ElementPropertyNames, $Name;

    # Store element names and units...
    $ElementPropertyNamesMap{$Name} = $Units;
  }

  # Process element data...
  LINE: while ($Line = GetTextLine(\*ELEMENTDATAFILE)) {
    if ($Line =~ /^#/) {
      next LINE;
    }
    @LineWords = ();
    @LineWords = quotewords($InDelim, 0, $Line);
    if (@LineWords != $NumOfCols) {
      croak "Error: The number of data fields, @LineWords, in $ElementDataFile must be $NumOfCols.\nLine: $Line...";
    }
    $AtomicNumber = $LineWords[0]; $ElementSymbol = $LineWords[1];
    if (exists $ElementDataMap{$AtomicNumber}) {
      carp "Warning: Ignoring data for element $ElementSymbol: It has already been loaded.\nLine: $Line....";
      next LINE;
    }

    # Store all the values...
    %{$ElementDataMap{$AtomicNumber}} = ();
    for $Index (0 .. $#LineWords) {
      $Name = $ElementPropertyNames[$Index];
      $Value = $LineWords[$Index];
      $ElementDataMap{$AtomicNumber}{$Name} = $Value;
    }
  }
  close ELEMENTDATAFILE;

  # Setup the element symbol map as well...
  _SetupElementSymbolMap();
}

#
# Load PeriodicTableIsotopeData.csv files from <MayaChemTools>/lib directory...
#
sub _LoadElementIsotopeData {
  my($ElementIsotopeDataFile) = @_;

  %ElementIsotopeDataMap = ();
  %ElementIsotopeDerivedDataMap = ();

  # Load isotope data for all elements...
  #
  # File format:
  # "Atomic Number","Isotope Symbol","Mass Number","Relative Atomic Mass","% Natural Abundnace"
  #
  # Empty values for "Relative Atomic Mass" and "% Natural Abundnace" imply absence of any
  # naturally occuring isotopes for the element.
  #
  my($InDelim, $Line, $NumOfCols, @ColLabels, @LineWords);

  $InDelim = "\,";
  open ISOTOPEDATAFILE, "$ElementIsotopeDataFile" or croak "Couldn't open $ElementIsotopeDataFile: $! ...";

  # Skip lines up to column labels...
  LINE: while ($Line = GetTextLine(\*ISOTOPEDATAFILE)) {
    if ($Line !~ /^#/) {
      last LINE;
    }
  }
  @ColLabels= quotewords($InDelim, 0, $Line);
  $NumOfCols = @ColLabels;

  my($AtomicNumber, $IsotopeSymbol, $MassNumber, $RelativeAtomicMass, $NaturalAbundance, %ZeroNaturalAbundanceMap);
  %ZeroNaturalAbundanceMap = ();

  # Process element data...
  LINE: while ($Line = GetTextLine(\*ISOTOPEDATAFILE)) {
    if ($Line =~ /^#/) {
      next LINE;
    }
    @LineWords = ();
    @LineWords = quotewords($InDelim, 0, $Line);
    if (@LineWords != $NumOfCols) {
      croak "Error: The number of data fields, @LineWords, in $ElementIsotopeDataFile must be $NumOfCols.\nLine: $Line...";
    }
    ($AtomicNumber, $IsotopeSymbol, $MassNumber, $RelativeAtomicMass, $NaturalAbundance) = @LineWords;
    if (exists $ZeroNaturalAbundanceMap{$AtomicNumber}) {
      # Only one isotope data line allowed for elements with no natural isotopes...
      carp "Warning: Ignoring isotope data for element with atomic number $AtomicNumber: Only one data line allowed for an element with no natural isotopes.\nLine: $Line...";
      next LINE;
    }
    if (IsEmpty($NaturalAbundance)) {
      $RelativeAtomicMass = 0;
      $NaturalAbundance = 0;
      $ZeroNaturalAbundanceMap{$AtomicNumber} = 1;
    }
    if (exists $ElementIsotopeDataMap{$AtomicNumber}) {
      # Additional data for an existing element...
      if (exists $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}) {
	carp "Warning: Ignoring isotope data for element with atomic number $AtomicNumber: It has already been loaded.\nLine: $Line...";
	next LINE;
      }
    }
    else {
      # Data for a new element...
      %{$ElementIsotopeDataMap{$AtomicNumber}} = ();
    }
    %{$ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}} = ();
    $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}{IsotopeSymbol} = $IsotopeSymbol;
    $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}{RelativeAtomicMass} = $RelativeAtomicMass;
    $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}{NaturalAbundance} = $NaturalAbundance;
  }
  close ISOTOPEDATAFILE;

  _SetupElementIsotopeDerivedDataMap();
}

#
# Map mass number of most abundant isotope for each element; additionally,
# count number of isotopes as well.
#
sub _SetupElementIsotopeDerivedDataMap {
  my($AtomicNumber, $MassNumber, $NaturalAbundance, $MostNaturalAbundance, $MostAbundantMassNumber, $IsotopeCount);

  %ElementIsotopeDerivedDataMap = ();

  for $AtomicNumber (sort {$a <=> $b} keys %ElementIsotopeDataMap) {
    $IsotopeCount = 0;
    $MostAbundantMassNumber = 0;
    $MostNaturalAbundance = 0;
    MASSNUMBER: for $MassNumber (sort {$a <=> $b} keys %{$ElementIsotopeDataMap{$AtomicNumber}}) {
      $NaturalAbundance = $ElementIsotopeDataMap{$AtomicNumber}{$MassNumber}{NaturalAbundance};
      if (IsEmpty($NaturalAbundance)) {
	# No natural isotopes available...
	$MostAbundantMassNumber = $MassNumber;
	last MASSNUMBER;
      }
      if ($NaturalAbundance == 0) {
	# Not a natural isotope; Listed in periodic table data file to support non-natural
	# elements such as T. It's not included in natural isotope count...
	next MASSNUMBER;
      }
      $IsotopeCount++;
      if ($NaturalAbundance > $MostNaturalAbundance) {
	$MostAbundantMassNumber = $MassNumber;
	$MostNaturalAbundance = $NaturalAbundance;
      }
    }
    %{$ElementIsotopeDerivedDataMap{$AtomicNumber}} = ();
    $ElementIsotopeDerivedDataMap{$AtomicNumber}{IsotopeCount} = $IsotopeCount;
    $ElementIsotopeDerivedDataMap{$AtomicNumber}{MostAbundantMassNumber} = $MostAbundantMassNumber;
  }
}

#
# Setup element symbol map...
#
sub _SetupElementSymbolMap {
  my($AtomicNumber, $ElementSymbol);

  %ElementSymbolMap = ();

  for $AtomicNumber (keys %ElementDataMap) {
    $ElementSymbol = $ElementDataMap{$AtomicNumber}{ElementSymbol};
    $ElementSymbolMap{$ElementSymbol} = $AtomicNumber;
  }
}

# Validate element ID...
sub _ValidateElementID {
  my($ElementID) = @_;
  my($ElementSymbol, $AtomicNumber);

  if ($ElementID =~ /[^0-9]/) {
    # Assume atomic symbol...
    $ElementSymbol = $ElementID;
    if (exists $ElementSymbolMap{$ElementSymbol}) {
      $AtomicNumber = $ElementSymbolMap{$ElementSymbol};
    }
    else {
      return undef;
    }
  }
  else {
    # Assume atomic number...
    $AtomicNumber = $ElementID;
    if (!exists $ElementDataMap{$AtomicNumber}) {
      return undef;
    }
  }
  return $AtomicNumber;
}

1;

__END__

=head1 NAME

PeriodicTable

=head1 SYNOPSIS

use PeriodicTable;

use PeriodicTable qw(:all);

=head1 DESCRIPTION

B<PeriodicTable> module provides the following functions:

GetElementMostAbundantNaturalIsotopeData,
GetElementMostAbundantNaturalIsotopeMass,
GetElementMostAbundantNaturalIsotopeMassNumber, GetElementNaturalIsotopeAbundance,
GetElementNaturalIsotopeCount, GetElementNaturalIsotopeMass,
GetElementNaturalIsotopesData, GetElementPropertiesData,
GetElementPropertiesNames, GetElementPropertiesNamesAndUnits,
GetElementPropertyUnits, GetElements, GetElementsByAmericanStyleGroupLabel,
GetElementsByEuropeanStyleGroupLabel, GetElementsByGroupName,
GetElementsByGroupNumber, GetElementsByPeriodNumber,
GetIUPACGroupNumberFromAmericanStyleGroupLabel,
GetIUPACGroupNumberFromEuropeanStyleGroupLabel, IsElement,
IsElementNaturalIsotopeMassNumber, IsElementProperty

=head1 METHODS

=over 4

=item B<GetElements>

    @ElementSymbols = GetElements();
    $ElementSymbolsRef = GetElements();

Returns an array or a reference to an array of known element symbols

=item B<GetElementsByGroupName>

    @ElementSymbols = GetElementsByGroupName($GroupName);
    $ElementSymbolsRef = GetElementsByGroupName($GroupName);

Returns an array or a reference to an array of element symbols for a specified I<GroupName>.
Supported I<GroupName> values are: I<Alkali metals, Alkaline earth metals, Coinage metals, Pnictogens,
Chalcogens, Halogens, Noble gases>; Additionally, usage of I<Lanthanides> (Lanthanoids)
and I<Actinides> (Actinoids) is also supported.

=item B<GetElementsByGroupNumber>

    @ElementSymbols = GetElementsByGroupNumber($GroupNumber);
    $ElementSymbolsRef = GetElementsByGroupNumber($GroupNumber);

Returns an array or a reference to an array of element symbols for a specified I<GroupNumber>

=item B<GetElementsByAmericanStyleGroupLabel>

    @ElementSymbols = GetElementsByAmericanStyleGroupLabel($GroupLabel);
    $ElementSymbolsRef = GetElementsByAmericanStyleGroupLabel($GroupLabel);

Returns an array or a reference to an array of element symbols for a specified American
style I<GroupLabel>. Valid values for Amercian style group labels: I<IA to VIIIA, IB to VIIIB, VIII>.

=item B<GetElementsByEuropeanStyleGroupLabel>

    @ElementSymbols = GetElementsByEuropeanStyleGroupLabel($GroupLabel);
    $ElementSymbolsRef = GetElementsByEuropeanStyleGroupLabel($GroupLabel);

Returns an array or a reference to an array of element symbols for a specified European
style I<GroupLabel>. Valid values for European style group labels: I<IA to VIIIA, IB to VIIIB, VIII>.

=item B<GetElementsByPeriodNumber>

    @ElementSymbols = GetElementsByPeriodNumber($PeriodNumber);
    $ElementSymbolsRef = GetElementsByPeriodNumber($PeriodNumber);

Returns an array or a reference to an array of element symbols for a specified I<PeriodNumber>.

=item B<GetElementMostAbundantNaturalIsotopeData>

    @IsotopeData = GetElementMostAbundantNaturalIsotopeData(
                   $ElementID);
    $IsotopeDataRef = GetElementMostAbundantNaturalIsotopeData(
                   $ElementID);

Returns an array or reference to an array containing data for most abundant isotope of
an element specfied by element symbol or atomic number. Isotope data arrays contain these
values: I<AtomicNumber, IsotopeSymbol, MassNumber, RelativeAtomicMass, and NaturalAbundance>.

=item B<GetElementMostAbundantNaturalIsotopeMassNumber>

    $MassNumber = GetElementMostAbundantNaturalIsotopeMassNumber($ElementID);

Returns mass number of most abundant natural isotope of an element specfied by element
symbol or atomic number

=item B<GetElementNaturalIsotopeCount>

    $IsotopeCount = GetElementNaturalIsotopeCount($ElementID);

Returns natural isotope count for an element specfied by element symbol or
atomic number

=item B<GetElementNaturalIsotopesData>

    $DataHashRef = GetElementNaturalIsotopesData($ElementID,
                   [$MassNumber]);

Reurns a reference to a hash containingall available isotope data for an element specified
using element symbol or aromic number; an optional mass number indicates retrieve data
for a specific isotope

=item B<GetElementNaturalIsotopeAbundance>

    $Abundance = GetElementNaturalIsotopeAbundance($ElementID,
                 $MassNumber);

Returns percent abundance of natural isotope for an element with specfic mass
number.

=item B<GetElementMostAbundantNaturalIsotopeMass>

    $RelativeAtomicMass = GetElementMostAbundantNaturalIsotopeMass(
                          $ElementID);

Returns relative atomic mass of most abundant isotope for an element specified using
element symbol or aromic number.

=item B<GetElementNaturalIsotopeMass>

    $RelativeAtomicMass = GetElementNaturalIsotopeMass($ElementID,
                          $MassNumber);

Returns relative atomic mass of an element with specfic mass number.

=item B<GetElementPropertiesData>

    $PropertyDataHashRef = GetElementPropertiesData($ElementID);

Returns a reference to a hash containing all available properties data for an element
specified using element symbol or atomic number.

=item B<GetElementPropertyName>

    $PropertyValue = GetElement<PropertyName>($ElementID);

Returns value of an element for a element specified using element symbol or atomic number.

These functions are not defined in this modules; these are implemented on-the-fly using
Perl's AUTOLOAD funcionality.

Here is the list of known element I<property names>: AllenElectronegativity,
AllredRochowElectronegativity, AtomicNumber, AtomicRadiusCalculated,
AtomicRadiusEmpirical, AtomicWeight, Block, BoilingPoint, BondLength,
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

=item B<GetElementPropertiesNames>

    @PropertyNames = GetElementPropertiesNames([$Mode]);
    $PropertyNamesRef = GetElementPropertiesNames([$Mode]);

Returns names of all available element properties. Optional mode parameter controls
grouping of property names; Possible values: I<ByGroup or Alphabetical>. Default:
I<ByGroup>.

=item B<GetElementPropertiesNamesAndUnits>

    $NameUnitsHashRef = GetElementPropertiesNamesAndUnits();

Returns a reference to a hash of property names and units of all available element
properties. Names with no units contains empty strings.

=item B<GetElementPropertyUnits>

    $Units = GetElementPropertyUnits($PropertyName);

Returns units for a specific element property name. An empty string is returned for
a property with no units.

=item B<GetIUPACGroupNumberFromAmericanStyleGroupLabel>

    $GroupNumber = GetIUPACGroupNumberFromAmericanStyleGroupLabel(
                   $GroupLabel);

Returns IUPAC group numbers of a specific American style group label. A comma delimited
string is returned for group VIII or VIIIB.

=item B<GetIUPACGroupNumberFromEuropeanStyleGroupLabel>

    $GroupNumber = GetIUPACGroupNumberFromEuropeanStyleGroupLabel(
                   $GroupLabel);

Returns IUPAC group numbers of a specific European style group label. A comma delimited
string is returned for group VIII or VIIIA.

=item B<IsElement>

    $Status = IsElement($ElementID);

Returns 1 or 0 based on whether it's a known element symbol or atomic number.

=item B<IsElementNaturalIsotopeMassNumber>

    $Status = IsElementNaturalIsotopeMassNumber($ElementID, $MassNumber);

Returns 1 or 0 based on whether it's a valid mass number for an element symbol
or atomic number.

=item B<IsElementProperty>

    $Status = IsElementProperty($PropertyName);

Returns 1 or 0 based on whether it's a valid property name.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

AminoAcids.pm, NucleicAcids.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
