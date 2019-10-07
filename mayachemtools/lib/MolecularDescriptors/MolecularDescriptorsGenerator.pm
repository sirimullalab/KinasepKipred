package MolecularDescriptors::MolecularDescriptorsGenerator;
#
# File: MolecularDescriptorsGenerator.pm
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
use Exporter;
use Scalar::Util ();
use ObjectProperty;
use TextUtil ();
use FileUtil ();
use Molecule;
use MolecularDescriptors::MolecularDescriptors;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(ObjectProperty Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GetAvailableDescriptorClassNames GetAvailableClassAndDescriptorNames GetAvailableDescriptorNames GetAvailableDescriptorNamesForDescriptorClass GetAvailableClassNameForDescriptorName GetRuleOf5DescriptorNames GetRuleOf3DescriptorNames IsDescriptorClassNameAvailable IsDescriptorNameAvailable);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, %DescriptorsDataMap);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyMolecularDescriptorsGenerator';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializeMolecularDescriptorsGenerator();

  $This->_InitializeMolecularDescriptorsGeneratorProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Load available molecular descriptor classes...
  _LoadMolecularDescriptorsData();

}

# Initialize object data...
#
sub _InitializeMolecularDescriptorsGenerator {
  my($This) = @_;

  # Type of desciptors to generate...
  #
  # The current release of MayaChemTools supports generation of four sets of
  # descriptors: All available descriptors, rule of 5 or 3 descriptors or a specified
  # set of descriptors.
  #
  # Possible values: All, RuleOf5, RuleOf3 or Specify
  #
  # RuleOf5 [ Ref 91 ] descriptor names: MolecularWeight, HydrogenBondDonors, HydrogenBondAcceptors,
  # SLogP. RuleOf5 states: MolecularWeight <= 500, HydrogenBondDonors <= 5, HydrogenBondAcceptors <= 10,
  # and logP <= 5.
  #
  # RuleOf3 [ Ref 92 ] descriptor names: MolecularWeight, RotatableBonds, HydrogenBondDonors,
  # HydrogenBondAcceptors, SLogP, TPSA. RuleOf3 states: MolecularWeight <= 300, RotatableBonds <= 3,
  # HydrogenBondDonors <= 3, HydrogenBondAcceptors <= 3, logP <= 3, and TPSA <= 60.
  #
  # For Specify value of Mode, a set of descritor names must be specified using
  # DescriptorNames parameter.
  #
  # Default: All
  #
  $This->{Mode} = '';

  # Descriptor names used to generate descriptor values during a specified descriptor
  # generation mode...
  #
  @{$This->{DescriptorNames}} = ();

  # Descriptor calculation control parameters for specified descriptor class names...
  #
  # These parameters are passed on to appropriate descriptor classes during
  # instantiations of descriptor class objects.
  #
  %{$This->{DescriptorClassParameters}} = ();

  $This->{DescriptorClassesInstantiated} = 0;

  # Descriptor class names and objects corresponding to specified descriptor names...
  #
  @{$This->{DescriptorClassNames}} = ();
  %{$This->{DescriptorClassObjects}} = ();

  # Descriptor values generated for specified descriptor names...
  #
  @{$This->{DescriptorValues}} = ();

  return $This;
}

# Initialize object properties...
#
sub _InitializeMolecularDescriptorsGeneratorProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  # Set default value for Mode...
  if (!$This->{Mode}) {
    $This->{Mode} = 'All';
  }

  $This->_CheckAndInitializeDescriptorNames();

  return $This;
}

# Set descriptors generation mode......
#
sub SetMode {
  my($This, $Value) = @_;

  # All - all available descriptors
  # Specify - Specified set of descriptors

  if ($Value !~ /^(All|RuleOf5|RuleOf3|Specify)$/i) {
    croak "Error: ${ClassName}->SetMode: Mode value, $Value, is not valid; Supported values: All, RuleOf5, RuleOf3 or Specify...";
  }

  $This->{Mode} = $Value;

  return $This;
}

# Set descriptor names to use for generating descriptor values using an array
# or reference to an array...
#
sub SetDescriptorNames {
  my($This, @Values) = @_;

  if ($This->{Mode} =~ /^All$/i) {
    croak "Error: ${ClassName}->SetDescriptorNames: Descriptor names cann't be specified during \"All\" value of descsriptors generation \"Mode\"...";
  }

  if (!@Values) {
    return;
  }

  my($FirstValue, $TypeOfFirstValue);

  $FirstValue = $Values[0];
  $TypeOfFirstValue = ref $FirstValue;

  @{$This->{DescriptorNames}} = ();

  if ($TypeOfFirstValue =~ /^ARRAY/) {
    # Initialize using array refernce...
    push @{$This->{DescriptorNames}}, @{$FirstValue};
  }
  else {
    # It's a list of values...
    push @{$This->{DescriptorNames}}, @Values;
  }

  # Make sure specified descriptor names are valid...
  $This->_ValidateDescriptorNames();

  return $This;
}

# Get descriptor names as an array...
#
sub GetDescriptorNames {
  my($This) = @_;

  return wantarray ? @{$This->{DescriptorNames}} : scalar @{$This->{DescriptorNames}};
}

# Get all descriptor values as an array...
#
sub GetDescriptorValues {
  my($This) = @_;

  if ($This->{DescriptorsGenerated}) {
    return wantarray ? @{$This->{DescriptorValues}} : scalar @{$This->{DescriptorValues}};
  }
  else {
    my(@DescriptorValues);

    @DescriptorValues = ('None') x scalar @{$This->{DescriptorNames}};

    return wantarray ? @DescriptorValues : scalar @DescriptorValues;
  }
}

# Get descriptor value for a specified descriptor name...
#
sub GetDescriptorValueByName {
  my($This, $Name) = @_;
  my(%NamesAndValues);

  %NamesAndValues = $This->GetDescriptorNamesAndValues();

  return exists $NamesAndValues{$Name} ? $NamesAndValues{$Name} : 'None';

}

# Get calculated molecular descriptor names sand values as a  hash with names
# and values as key/value pairs...
#
sub GetDescriptorNamesAndValues {
  my($This) = @_;
  my(%NamesAndValues);

  %NamesAndValues = ();
  @NamesAndValues{ @{$This->{DescriptorNames}} } = $This->GetDescriptorValues();

  return %NamesAndValues;
}

# Set up descriptor calculation control parameters for a specified descriptor class name...
#
# The specified parameter names and values are simply passed on to specified descriptor
# class during instantiation of descriptor class object without any performing any validation
# of parameter names and associated values. It's up to the appropriate descriptor class methods
# to validate these parameters and values.
#
# In addition to specified parameter names and values, the parameter hash must also contain
# descriptor class name as key and value pair with DescriptorClassName as key with class
# name as value.
#
sub SetDescriptorClassParameters {
  my($This, %NamesAndValues) = @_;
  my($DescriptorClassName, $Name, $Value);

  if (!exists $NamesAndValues{DescriptorClassName}) {
    croak "Error: ${ClassName}->_SetDescriptorNameParameters: Can't set descriptor class name paramaters: DescriptorClassName is not specified...";
  }

  $DescriptorClassName = $NamesAndValues{DescriptorClassName};
  if (!IsDescriptorClassNameAvailable($DescriptorClassName)) {
    carp "Warning: ${ClassName}->_SetDescriptorClassParameters: Can't set descriptor class name paramaters: Specified descriptor class name, $DescriptorClassName, is not available...";
    return $This;
  }

  if (exists $This->{DescriptorClassParameters}{$DescriptorClassName}) {
    carp "Warning: ${ClassName}->SetDescriptorClassParameters: Class name parameters for $DescriptorClassName have already been specified: Replacing previous values...";
  }

  %{$This->{DescriptorClassParameters}{$DescriptorClassName}} = ();
  NAME: while (($Name, $Value) = each  %NamesAndValues) {
    if ($Name =~ /^DescriptorClassName$/) {
      next NAME;
    }
    $This->{DescriptorClassParameters}{$DescriptorClassName}{$Name} = $Value;
  }

  return $This;
}

# Get descriptor name parameters as a reference to hash of hashes with hash
# keys corresponding to class name and class parameter name with hash value
# as class parameter value...
#
sub GetDescriptorClassParameters {
  my($This) = @_;

  return \%{$This->{DescriptorClassParameters}};
}

# Get available descriptor class names as an array.
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAvailableDescriptorClassNames {

  return wantarray ? @{$DescriptorsDataMap{ClassNames}} : scalar @{$DescriptorsDataMap{ClassNames}};
}

# Get available descriptors class and descriptors names as a hash containing key/value
# pairs corresponding to class name and an array of descriptor names available for the
# class.
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAvailableClassAndDescriptorNames {
  my($DescriptorClassName, @DescriptorNames, %ClassAndDescriptorNames);

  %ClassAndDescriptorNames = ();
  for $DescriptorClassName (@{$DescriptorsDataMap{ClassNames}}) {
    @{$ClassAndDescriptorNames{$DescriptorClassName}} = ();
    push @{$ClassAndDescriptorNames{$DescriptorClassName}}, @{$DescriptorsDataMap{ClassToDescriptorNames}{$DescriptorClassName}};
  }

  return %ClassAndDescriptorNames;
}

# Get available descriptor names as an array.
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAvailableDescriptorNames {
  my(@DescriptorNames);

  @DescriptorNames = ();
  push @DescriptorNames, map { @{$DescriptorsDataMap{ClassToDescriptorNames}{$_}} } @{$DescriptorsDataMap{ClassNames}};

  return wantarray ? @DescriptorNames : scalar @DescriptorNames;
}

# Is it a valid descriptors class name?
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub IsDescriptorClassNameAvailable {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $DescriptorClassName);

  if ((@_ == 2) && (_IsMolecularDescriptorsGenerator($FirstParameter))) {
    ($This, $DescriptorClassName) = ($FirstParameter, $SecondParameter);
  }
  else {
    ($DescriptorClassName) = ($FirstParameter);
  }

  return (exists $DescriptorsDataMap{ClassToDescriptorNames}{$DescriptorClassName}) ? 1 : 0;
}

# Is it a valid descriptor name?
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub IsDescriptorNameAvailable {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $DescriptorName);

  if ((@_ == 2) && (_IsMolecularDescriptorsGenerator($FirstParameter))) {
    ($This, $DescriptorName) = ($FirstParameter, $SecondParameter);
  }
  else {
    ($DescriptorName) = ($FirstParameter);
  }

  return (exists $DescriptorsDataMap{DescriptorToClassName}{$DescriptorName}) ? 1 : 0;
}

# Get available descriptors names for a descriptor class as an array.
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAvailableDescriptorNamesForDescriptorClass {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $DescriptorClassName, @DescriptorNames);

  if ((@_ == 2) && (_IsMolecularDescriptorsGenerator($FirstParameter))) {
    ($This, $DescriptorClassName) = ($FirstParameter, $SecondParameter);
  }
  else {
    ($DescriptorClassName) = ($FirstParameter);
  }

  @DescriptorNames = ();
  if (exists $DescriptorsDataMap{ClassToDescriptorNames}{$DescriptorClassName}) {
    push @DescriptorNames, @{$DescriptorsDataMap{ClassToDescriptorNames}{$DescriptorClassName}};
  }

  return wantarray ? @DescriptorNames : scalar @DescriptorNames;
}

# Get available descriptors class name for a descriptor name.
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetAvailableClassNameForDescriptorName {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $DescriptorClassName, $DescriptorName);

  if ((@_ == 2) && (_IsMolecularDescriptorsGenerator($FirstParameter))) {
    ($This, $DescriptorName) = ($FirstParameter, $SecondParameter);
  }
  else {
    ($DescriptorName) = ($FirstParameter);
  }

  $DescriptorClassName = '';
  if (exists $DescriptorsDataMap{DescriptorToClassName}{$DescriptorName}) {
    $DescriptorClassName = $DescriptorsDataMap{DescriptorToClassName}{$DescriptorName};
  }

  return $DescriptorClassName;
}

# Get RuleOf5 descriptor names as an array.
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetRuleOf5DescriptorNames {
  my(@DescriptorNames);

  @DescriptorNames = qw(MolecularWeight HydrogenBondDonors HydrogenBondAcceptors SLogP);

  return wantarray ? @DescriptorNames : scalar @DescriptorNames;
}

# Get RuleOf3 descriptor names as an array.
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetRuleOf3DescriptorNames {
  my(@DescriptorNames);

  @DescriptorNames = qw(MolecularWeight RotatableBonds HydrogenBondDonors HydrogenBondAcceptors SLogP TPSA);

  return wantarray ? @DescriptorNames : scalar @DescriptorNames;
}


# Set molecule object...
#
sub SetMolecule {
  my($This, $Molecule) = @_;

  $This->{Molecule} = $Molecule;

  # Weaken the reference to disable increment of reference count...
  Scalar::Util::weaken($This->{Molecule});

  return $This;
}

# Generate specified molecular descriptors...
#
# After instantiating descriptor class objects at first invocation and  initialializing
# descriptor values during subsequent invocations, GenerateDescriptors method
# provided by each descriptor class is used to calculate descriptor values for
# specified descriptors.
#
sub GenerateDescriptors {
  my($This) = @_;
  my($DescriptorClassName, $DescriptorClassObject);

  # Initialize descriptor values...
  $This->_InitializeDescriptorValues();

  # Instantiate decriptor classed corresponding to specified descriptors...
  if (!$This->{DescriptorClassesInstantiated}) {
    $This->_InstantiateDescriptorClasses();
  }

  # Check availability of molecule...
  if (!$This->{Molecule}) {
    carp "Warning: ${ClassName}->GenerateDescriptors: $This->{Type} molecular descriptors generation didn't succeed: Molecule data is not available: Molecule object hasn't been set...";
    return undef;
  }

  # Calculate descriptor values...
  for $DescriptorClassName (@{$This->{DescriptorClassNames}}) {
    $DescriptorClassObject = $This->{DescriptorClassObjects}{$DescriptorClassName};

    $DescriptorClassObject->SetMolecule($This->{Molecule});
    $DescriptorClassObject->GenerateDescriptors();

    if (!$DescriptorClassObject->IsDescriptorsGenerationSuccessful()) {
      return undef;
    }
  }

  # Set final descriptor values...
  $This->_SetFinalDescriptorValues();

  return $This;
}

# Initialize descriptor values...
#
sub _InitializeDescriptorValues {
  my($This) = @_;

  $This->{DescriptorsGenerated} = 0;

  @{$This->{DescriptorValues}} = ();

  return $This;
}

# Setup final descriptor values...
#
sub _SetFinalDescriptorValues {
  my($This) = @_;
  my($DescriptorName, $DescriptorClassName, $DescriptorClassObject);

  $This->{DescriptorsGenerated} = 1;

  @{$This->{DescriptorValues}} = ();

  if ($This->{Mode} =~ /^All$/i) {
    # Set descriptor values for all available descriptors...
    for $DescriptorClassName (@{$This->{DescriptorClassNames}}) {
      $DescriptorClassObject = $This->{DescriptorClassObjects}{$DescriptorClassName};

      push @{$This->{DescriptorValues}}, $DescriptorClassObject->GetDescriptorValues();
    }
  }
  else {
    # Set descriptor values for a subset of available descriptors...
    for $DescriptorName (@{$This->{DescriptorNames}}) {
      $DescriptorClassName = $DescriptorsDataMap{DescriptorToClassName}{$DescriptorName};
      $DescriptorClassObject = $This->{DescriptorClassObjects}{$DescriptorClassName};

      push @{$This->{DescriptorValues}}, $DescriptorClassObject->GetDescriptorValueByName($DescriptorName);
    }
  }

  return $This;
}

# Is descriptors generation successful?
#
# Notes:
#   . After successful generation of descriptor values by each descriptor class
#     corresponding to specified descriptor names, DescriptorsCalculated  to 1;
#     otherwise, it's set to 0.
#
sub IsDescriptorsGenerationSuccessful {
  my($This) = @_;

  return $This->{DescriptorsGenerated} ? 1 : 0;
}

# Check and set default descriptor names for generating descriptor values...
#
sub _CheckAndInitializeDescriptorNames {
  my($This) = @_;

  if ($This->{Mode} =~ /^(All|RuleOf5|RuleOf3)$/i) {
    if (@{$This->{DescriptorNames}}) {
      croak "Error: ${ClassName}->_CheckAndInitializeDescriptorNames: Descriptor names can't be specified during \"All, RuleOf5 or RuleOf3\" values of descsriptors generation \"Mode\"...";
    }
  }

  if ($This->{Mode} =~ /^All$/i) {
    @{$This->{DescriptorNames}} = GetAvailableDescriptorNames();
  }
  elsif ($This->{Mode} =~ /^RuleOf5$/i) {
    @{$This->{DescriptorNames}} = GetRuleOf5DescriptorNames();
  }
  elsif ($This->{Mode} =~ /^RuleOf3$/i) {
    @{$This->{DescriptorNames}} = GetRuleOf3DescriptorNames();
  }
  elsif ($This->{Mode} =~ /^Specify$/i) {
    if (!@{$This->{DescriptorNames}}) {
      croak "Error: ${ClassName}->_CheckAndInitializeDescriptorNames: DescriptorNames must be specified during Specify value for Mode...";
    }
  }
  else {
    croak "Error: ${ClassName}->_CheckAndInitializeDescriptorNames: Mode value, $This->{Mode}, is not valid...";
  }
}

# Instantiate descriptor classes corresponding to specified descriptor names...
#
sub _InstantiateDescriptorClasses {
  my($This) = @_;
  my($DescriptorClassName, $DescriptorName, $DescriptorClassPath);

  $This->{DescriptorClassesInstantiated} = 1;

  @{$This->{DescriptorClassNames}} = ();
  %{$This->{DescriptorClassObjects}} = ();

  NAME: for $DescriptorName (@{$This->{DescriptorNames}}) {
    $DescriptorClassName = $DescriptorsDataMap{DescriptorToClassName}{$DescriptorName};

    if (exists $This->{DescriptorClassObjects}{$DescriptorClassName}) {
      next NAME;
    }
    push @{$This->{DescriptorClassNames}}, $DescriptorClassName;

    $DescriptorClassPath = $DescriptorsDataMap{ClassNameToClassPath}{$DescriptorClassName};

    if (exists $This->{DescriptorClassParameters}{$DescriptorClassName}) {
      $This->{DescriptorClassObjects}{$DescriptorClassName} = $DescriptorClassPath->new(%{$This->{DescriptorClassParameters}{$DescriptorClassName}});
    }
    else {
      $This->{DescriptorClassObjects}{$DescriptorClassName} = $DescriptorClassPath->new();
    }
  }

  return $This;
}

# Return a string containg data for MolecularDescriptorsGenerator object...
#
sub StringifyMolecularDescriptorsGenerator {
  my($This) = @_;
  my($TheString, $NamesAndValuesString, $Name, $Value, @NamesAndValuesInfo, %NamesAndValues);

  # Type of MolecularDescriptors...
  $TheString = "MolecularDescriptorsGenerator: Mode - $This->{Mode}; SpecifiedDescriptorNames - < @{$This->{DescriptorNames}} >; AvailableMolecularDescriptorClassNames - < @{$DescriptorsDataMap{ClassNames}} >";

  @NamesAndValuesInfo = ();
  %NamesAndValues = $This->GetDescriptorNamesAndValues();

  for $Name (@{$This->{DescriptorNames}}) {
    $Value = $NamesAndValues{$Name};
    $Value = (TextUtil::IsEmpty($Value) || $Value =~ /^None$/i) ? 'None' : $Value;
    push @NamesAndValuesInfo, "$Name - $Value";
  }
  if (@NamesAndValuesInfo) {
    $TheString .= "Names - Values: <" . TextUtil::JoinWords(\@NamesAndValuesInfo, ", ", 0) . ">";
  }
  else {
    $TheString .= "Names - Values: < None>";
  }

  return $TheString;
}

# Is it a MolecularDescriptorsGenerator object?
sub _IsMolecularDescriptorsGenerator {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Validate descriptor names for generating descriptor values...
#
sub _ValidateDescriptorNames {
  my($This) = @_;
  my($DescriptorName);

  for $DescriptorName (@{$This->{DescriptorNames}}) {
    if (!exists $DescriptorsDataMap{DescriptorToClassName}{$DescriptorName}) {
      croak "Error: ${ClassName}->_SetAndValidateDescriptorNames: Specified descriptor name, $DescriptorName, is not valid...";
    }
  }

  return $This;
}

#
# Load available molecular descriptors data...
#
# All available molecular descriptors classes are automatically detected in
# MolecularDescriptors directory under <MayaChemTools>/lib directory and
# information about available descriptor names is retrieved from each descriptor
# class using function GetDescriptorNames. The following %DescriptorsDataMap
# is setup containing all available molecular descriptors data:
#
#   @{$DescriptorsDataMap{ClassNames}}
#   %{$DescriptorsDataMap{ClassNameToPath}}
#   %{$DescriptorsDataMap{ClassToDescriptorNames}}
#   %{$DescriptorsDataMap{DescriptorToClassName}}
#
# GenerateDescriptors method is invoked fo each specified descriptor class
# object to calculate descriptor values for specified descriptors. After successful
# calculation of descriptors, GetDescriptorValues or GetDescriptorValueByName
# methods provided by descriptor objects are used to retrieve calculated
# descriptor values.
#
sub _LoadMolecularDescriptorsData {

  %DescriptorsDataMap = ();

  _RetrieveAndLoadDescriptorClasses();
  _SetupDescriptorsDataMap();
}

#
# Retrieve available molecular descriptors classes from MolecularDescriptors directory under
# <MayaChemTools>/lib directory...
#
sub _RetrieveAndLoadDescriptorClasses {
  my($DescriptorsDirName, $MayaChemToolsLibDir, $DescriptorsDirPath, $IncludeDirName, $DescriptorClassName, $DescriptorClassPath, $DescriptorsClassFileName, @FileNames, @DescriptorsClassFileNames);

  @{$DescriptorsDataMap{ClassNames}} = ();
  %{$DescriptorsDataMap{ClassNameToPath}} = ();

  $DescriptorsDirName = "MolecularDescriptors";
  $MayaChemToolsLibDir = FileUtil::GetMayaChemToolsLibDirName();

  $DescriptorsDirPath = "$MayaChemToolsLibDir/$DescriptorsDirName";

  if (! -d "$DescriptorsDirPath") {
    croak "Error: ${ClassName}::_RetrieveAndLoadDescriptorClasses: MayaChemTools package molecular descriptors directory, $DescriptorsDirPath, is missing: Possible installation problems...";
  }

  @FileNames = ("$DescriptorsDirPath/*");
  $IncludeDirName = 0;
  @DescriptorsClassFileNames = FileUtil::ExpandFileNames(\@FileNames, "pm", $IncludeDirName);

  if (!@DescriptorsClassFileNames) {
    croak "Error: ${ClassName}::_RetrieveAndLoadDescriptorClasses: MayaChemTools package molecular descriptors directory, $DescriptorsDirPath, doesn't contain any molecular descriptor class: Possible installation problems...";
  }

  FILENAME: for $DescriptorsClassFileName (sort @DescriptorsClassFileNames) {
    if ($DescriptorsClassFileName !~ /\.pm/) {
      croak "Error: ${ClassName}::_RetrieveAndLoadDescriptorClasses: MayaChemTools package molecular descriptors directory, $DescriptorsDirPath, contains invalid class file name $DescriptorsClassFileName: Possible installation problems...";
    }

    # Ignore base class and descriptors generator class...
    if ($DescriptorsClassFileName =~ /^(MolecularDescriptorsGenerator\.pm|MolecularDescriptors\.pm)$/) {
      next FILENAME;
    }

    ($DescriptorClassName) = split /\./, $DescriptorsClassFileName;
    $DescriptorClassPath = "${DescriptorsDirName}::${DescriptorClassName}";

    # Load descriptors class...
    eval "use $DescriptorClassPath";

    if ($@) {
      croak "Error: ${ClassName}::_RetrieveAndLoadDescriptorClasses: use $DescriptorClassPath failed: $@ ...";
    }

    push @{$DescriptorsDataMap{ClassNames}}, $DescriptorClassName;

    $DescriptorsDataMap{ClassNameToClassPath}{$DescriptorClassName} = $DescriptorClassPath;
  }
}

#
# Setup descriptors data map using loaded descriptor classes...
#
sub _SetupDescriptorsDataMap {
  my($DescriptorClassName, $DescriptorName, $DescriptorClassPath, @DescriptorNames);

  # Class to decriptor names map...
  %{$DescriptorsDataMap{ClassToDescriptorNames}} = ();

  # Descriptor to class name map...
  %{$DescriptorsDataMap{DescriptorToClassName}} = ();

  for $DescriptorClassName (@{$DescriptorsDataMap{ClassNames}}) {
    $DescriptorClassPath = $DescriptorsDataMap{ClassNameToClassPath}{$DescriptorClassName};

    @DescriptorNames = $DescriptorClassPath->GetDescriptorNames();

    if (!@DescriptorNames) {
      croak "Error: ${ClassName}::_SetupDescriptorsDataMap: Molecular descriptor class $DescriptorClassName doesn't provide any descriptor names...";
    }

    if (exists $DescriptorsDataMap{ClassToDescriptorNames}{$DescriptorClassName} ) {
      croak "Error: ${ClassName}::_SetupDescriptorsDataMap: Molecular descriptor class $DescriptorClassName has already been processed...";
    }

    @{$DescriptorsDataMap{ClassToDescriptorNames}{$DescriptorClassName}} = ();
    @{$DescriptorsDataMap{ClassToDescriptorNames}{$DescriptorClassName}} = @DescriptorNames;

    for $DescriptorName (@DescriptorNames) {
      if (exists $DescriptorsDataMap{DescriptorToClassName}{$DescriptorName}) {
	croak "Error: ${ClassName}::_SetupDescriptorsDataMap: Molecular descriptor name, $DescriptorName, in class name, $DescriptorClassName, has already been provided by class name $DescriptorsDataMap{DescriptorToClassName}{$DescriptorName}...";
      }

      $DescriptorsDataMap{DescriptorToClassName}{$DescriptorName} = $DescriptorClassName;
    }
  }
}

1;

__END__

=head1 NAME

MolecularDescriptorsGenerator

=head1 SYNOPSIS

use MolecularDescriptors::MolecularDescriptorsGenerator;

use MolecularDescriptors::MolecularDescriptorsGenerator qw(:all);

=head1 DESCRIPTION

B<MolecularDescriptorsGenerator> class provides the following methods:

new, GenerateDescriptors, GetAvailableClassAndDescriptorNames,
GetAvailableClassNameForDescriptorName, GetAvailableDescriptorClassNames,
GetAvailableDescriptorNames, GetAvailableDescriptorNamesForDescriptorClass,
GetDescriptorClassParameters, GetDescriptorNames, GetDescriptorNamesAndValues,
GetDescriptorValueByName, GetDescriptorValues, GetRuleOf3DescriptorNames,
GetRuleOf5DescriptorNames, IsDescriptorClassNameAvailable,
IsDescriptorNameAvailable, IsDescriptorsGenerationSuccessful,
SetDescriptorClassParameters, SetDescriptorNames, SetMode, SetMolecule,
StringifyMolecularDescriptorsGenerator

B<MolecularDescriptorsGenerator> is derived from is derived from B<ObjectProperty>
base class that provides methods not explicitly defined in B<MolecularDescriptorsGenerator>
or B<ObjectProperty> classes using Perl's AUTOLOAD functionality. These methods are
generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

B<MolecularDescriptorsGenerator> is designed to provide a plug-in environment for
molecular descriptors development. The molecular descriptor class modules available
in B<MolecularDescriptors> directory under B<MayaChemTools/lib> directory are
automatically detected and loaded into the system. The descriptor names provided
by each descriptor class module through its B<GetDescriptorNames> function are
retrieved and are made available for calculations of their values for a specified
molecule.

Any combination of available descriptor names can be specified during calculation
of descriptor values using B<GenerateDescriptors> method. The current release of
MayaChemTools supports generation of four sets of descriptors: All available
descriptors, rule of 5 or 3 descriptors, or a specified set of descriptor names.

RuleOf5 [ Ref 91 ] descriptor names are: MolecularWeight, HydrogenBondDonors,
HydrogenBondAcceptors, SLogP. RuleOf5 states: MolecularWeight <= 500,
HydrogenBondDonors <= 5, HydrogenBondAcceptors <= 10, and logP <= 5.

RuleOf3 [ Ref 92 ] descriptor names are: MolecularWeight, RotatableBonds,
HydrogenBondDonors, HydrogenBondAcceptors, SLogP, TPSA. RuleOf3 states:
MolecularWeight <= 300, RotatableBonds <= 3, HydrogenBondDonors <= 3,
HydrogenBondAcceptors <= 3, logP <= 3, and TPSA <= 60.

Before calculation of a specified set of descriptors by B<GenerateDescriptors>
method, a set of descriptor calculation control parameters for a specific descriptor
class name can be set using B<SetDescriptorClassParameters> method. The specified
control parameter names and values are simply passed on to specified descriptor
class during instantiation of descriptor class object without performing any validation
of parameter names and associated values. It's up to the appropriate descriptor class methods
to validate these parameters and values. In addition to specified parameter names and
values, the parameter hash must also contain descriptor class name as key and
value pair with DescriptorClassName as key with class name as value.

=head2 METHODS

=over 4

=item B<new>

    $NewMolecularDescriptorsGenerator = new MolecularDescriptors::
                                        MolecularDescriptorsGenerator(
                                        %NamesAndValues);

Using specified I<MolecularDescriptorsGenerator> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<MolecularDescriptorsGenerator>
object. By default, the following properties are initialized:

    Mode = 'All'
    @{$This->{DescriptorNames}} = ()
    %{$This->{DescriptorClassParameters}} = ()
    @{$This->{DescriptorClassNames}} = ()
    %{$This->{DescriptorClassObjects}} = ()
    @{$This->{DescriptorValues}} = ()

Examples:

    $MolecularDescriptorsGenerator = new MolecularDescriptors::
                                     MolecularDescriptorsGenerator(
                              'Molecule' => $Molecule);

    @DescriptorNames = qw(MolecularWeight HydrogenBondDonors Fsp3Carbons)
    $MolecularDescriptorsGenerator = new MolecularDescriptors::
                                     MolecularDescriptorsGenerator(
                              'Mode' => 'Specify',
                              'DescriptorNames' => \@DescriptorNames);

    $MolecularDescriptorsGenerator->SetDescriptorClassParameters(
                              'DescriptorClassName' => 'WeightAndMassDescriptors',
                              'WeightPrecision' => 2,
                              'MassPrecision' => 2);

    $MolecularDescriptorsGenerator->SetDescriptorClassParameters(
                              'DescriptorClassName' => 'HydrogenBondsDescriptors',
                              'HydrogenBondsType' => 'HBondsType1');

    $MolecularDescriptorsGenerator->SetMolecule($Molecule);
    $MolecularDescriptorsGenerator->GenerateDescriptors();
    print "MolecularDescriptorsGenerator: $MolecularDescriptorsGenerator\n";


=item B<GenerateDescriptors>

    $MolecularDescriptorsGenerator->GenerateDescriptors();

Calculates descriptor values for specified descriptors and returns I<MolecularDescriptorsGenerator>.

Descriptor class objects are instantiated only once at first invocation. During
subsequent calls to B<GenerateDescriptors> method, descriptor values are
initialized and B<GenerateDescriptors> method provided by descriptor class is
used to calculate descriptor values for specified descriptors.

=item B<GetAvailableClassAndDescriptorNames>

    %ClassAndDescriptorNames = $MolecularDescriptorsGenerator->
                              GetAvailableClassAndDescriptorNames();
    %ClassAndDescriptorNames = MolecularDescriptors::
                               MolecularDescriptorsGenerator::
                               GetAvailableClassAndDescriptorNames();

Returns available descriptors class and descriptors names as a hash containing key
and value pairs corresponding to class name and an array of descriptor names
available for the class.

=item B<GetAvailableClassNameForDescriptorName>

    $DescriptorClassName = $MolecularDescriptorsGenerator->
                      GetAvailableClassNameForDescriptorName($DescriptorName);

    $DescriptorClassName = MolecularDescriptors::MolecularDescriptorsGenerator::
                      GetAvailableClassNameForDescriptorName($DescriptorName);

Returns available descriptor class name for a descriptor name.

=item B<GetAvailableDescriptorClassNames>

    $Return = $MolecularDescriptorsGenerator->GetAvailableDescriptorClassNames();

    @DescriptorClassNames = $MolecularDescriptorsGenerator->
                              GetAvailableDescriptorClassNames();
    @DescriptorClassNames = MolecularDescriptors::
                              MolecularDescriptorsGenerator::
                              GetAvailableDescriptorClassNames();

Returns available descriptor class names as an array or number of available descriptor
class names in scalar context.

=item B<GetAvailableDescriptorNames>

    @DescriptorNames = $MolecularDescriptorsGenerator->
                              GetAvailableDescriptorNames();
    @DescriptorNames = MolecularDescriptors::
                              MolecularDescriptorsGenerator::
                              GetAvailableDescriptorNames();

Returns available descriptor names as an array or number of available descriptor
names in scalar context.

=item B<GetAvailableDescriptorNamesForDescriptorClass>

    @DescriptorNames = $MolecularDescriptorsGenerator->
          GetAvailableDescriptorNamesForDescriptorClass($DescriptorClassName);
    @DescriptorNames = MolecularDescriptors::
                       MolecularDescriptorsGenerator::
          GetAvailableDescriptorNamesForDescriptorClass($DescriptorClassName);

Returns available descriptors names for a descriptor class as an array or number
of available descriptor names in scalar context.

=item B<GetDescriptorClassParameters>

    $DescriptorClassParametersRef = $MolecularDescriptorsGenerator->
                              GetDescriptorClassParameters();
    $DescriptorClassParametersRef = MolecularDescriptors::
                                    MolecularDescriptorsGenerator::
                                    GetDescriptorClassParameters();

Returns descriptor name parameters as a reference to hash of hashes with hash
keys corresponding to class name and class parameter name with hash value
as class parameter value.

=item B<GetDescriptorNames>

    @DescriptorNames = $MolecularDescriptorsGenerator->GetDescriptorNames();
    @DescriptorNames = MolecularDescriptors::MolecularDescriptorsGenerator::
                       GetDescriptorNames();

Returns all available descriptor names as an array or number of available descriptors
in scalar context.

=item B<GetDescriptorNamesAndValues>

    %NamesAndValues = $MolecularDescriptorsGenerator->
                              GetDescriptorNamesAndValues();

Returns calculated molecular descriptor names and values as a hash with descriptor
names and values as hash key and value pairs.

=item B<GetDescriptorValueByName>

    $Value = $MolecularDescriptorsGenerator->
                              GetDescriptorValueByName($Name);

Returns calculated descriptor values for a specified descriptor name.

=item B<GetDescriptorValues>

    @DescriptorValues = $MolecularDescriptorsGenerator->GetDescriptorValues();

Returns all calculated descriptor values as an array corresponding to specified
descriptor names.

=item B<GetRuleOf3DescriptorNames>

    @DescriptorNames = $MolecularDescriptorsGenerator->
                       GetRuleOf3DescriptorNames();
    @DescriptorNames = MolecularDescriptors::
                       MolecularDescriptorsGenerator::
                       GetRuleOf3DescriptorNames();

Returns rule of 3  descriptor names as an array or number of rule of 3 descriptors in scalar
context.

RuleOf3 [ Ref 92 ] descriptor names are: MolecularWeight, RotatableBonds,
HydrogenBondDonors, HydrogenBondAcceptors, SLogP, TPSA. RuleOf3 states:
MolecularWeight <= 300, RotatableBonds <= 3, HydrogenBondDonors <= 3,
HydrogenBondAcceptors <= 3, logP <= 3, and TPSA <= 60.

=item B<GetRuleOf5DescriptorNames>

    @DescriptorNames = $MolecularDescriptorsGenerator->
                              GetRuleOf5DescriptorNames();
    @DescriptorNames = $MolecularDescriptorsGenerator::
                             GetRuleOf5DescriptorNames();

Returns rule of 5  descriptor names as an array or number of rule of 4 descriptors in scalar
context.

RuleOf5 [ Ref 91 ] descriptor names are: MolecularWeight, HydrogenBondDonors,
HydrogenBondAcceptors, SLogP. RuleOf5 states: MolecularWeight <= 500,
HydrogenBondDonors <= 5, HydrogenBondAcceptors <= 10, and logP <= 5.

=item B<IsDescriptorClassNameAvailable>

    $Status = $MolecularDescriptorsGenerator->
                              IsDescriptorClassNameAvailable($ClassName);
    $Status = MolecularDescriptors::
                              MolecularDescriptorsGenerator::
                              IsDescriptorClassNameAvailable($ClassName);

Returns 1 or 0 based on whether specified descriptor class name is available.

=item B<IsDescriptorNameAvailable>

    $Status = $MolecularDescriptorsGenerator->
                              IsDescriptorNameAvailable($DescriptorName);
    $Status = MolecularDescriptors::
                              MolecularDescriptorsGenerator::
                              IsDescriptorNameAvailable($DescriptorName);

Returns 1 or 0 based on whether specified descriptor name is available.

=item B<IsDescriptorsGenerationSuccessful>

    $Status = $MolecularDescriptorsGenerator->
                              IsDescriptorsGenerationSuccessful();

Returns 1 or 0 based on whether descriptors generation is successful.

=item B<SetDescriptorClassParameters>

    $MolecularDescriptorsGenerator->SetDescriptorClassParameters(
                              %NamesAndValues);

Sets descriptor calculation control parameters for a specified descriptor class name
and returns I<MolecularDescriptorsGenerator>.

The specified parameter names and values are simply passed on to specified descriptor
class during instantiation of descriptor class object without any performing any validation
of parameter names and associated values. It's up to the appropriate descriptor class methods
to validate these parameters and values.

In addition to specified parameter names and values, the parameter hash must also contain
descriptor class name as key and value pair with DescriptorClassName as key with class
name as value.

=item B<SetDescriptorNames>

    $MolecularDescriptorsGenerator->SetDescriptorNames(@Names);
    $MolecularDescriptorsGenerator->SetDescriptorNames(\@Names);

Sets descriptor names to use for generating descriptor values using an array
or reference to an array and returns I<MolecularDescriptorsGenerator>.

=item B<SetMode>

    $MolecularDescriptorsGenerator->SetMode($Mode);

Sets descriptors generation mode and returns I<MolecularDescriptorsGenerator>.
Possible I<Mode> values: I<All, RuleOf5, RuleOf3, Specify>.

=item B<SetMolecule>

    $MolecularDescriptorsGenerator->SetMolecule($Molecule);

Sets molecule to use during calculation of molecular descriptors and returns
I<MolecularDescriptorsGenerator>.

=item B<StringifyMolecularDescriptorsGenerator>

    $String = $MolecularDescriptorsGenerator->StringifyMolecularDescriptorsGenerator();

Returns a string containing information about I<MolecularDescriptorsGenerator> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

MolecularDescriptors.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
