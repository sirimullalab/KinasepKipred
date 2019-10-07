package Fingerprints::EStateIndiciesFingerprints;
#
# File: EStateIndiciesFingerprints.pm
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
use Text::ParseWords;
use TextUtil ();
use FileUtil ();
use MathUtil ();
use Fingerprints::Fingerprints;
use Molecule;
use AtomTypes::EStateAtomTypes;
use AtomicDescriptors::EStateValuesDescriptors;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Fingerprints::Fingerprints Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyEStateIndiciesFingerprints';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeEStateIndiciesFingerprints();

  $This->_InitializeEStateIndiciesFingerprintsProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeEStateIndiciesFingerprints {
  my($This) = @_;

  # EStateIndicies is a vector containing sum of E-state values for E-state atom types
  #
  $This->{Type} = 'EStateIndicies';

  # EStateAtomTypesSetToUse for EStateIndicies:
  #
  # ArbitrarySize - Corrresponds to only E-state atom types detected in molecule
  # FixedSize - Corresponds to fixed number of E-state atom types previously defined [ Ref 77 ]
  #
  # The default EStateAtomTypesSetToUse value for EStateIndicies fingerprints type: ArbitrarySize.
  # Possible values: ArbitrarySize or FixedSize.
  #
  $This->{EStateAtomTypesSetToUse} = '';

  # Assigned E-state atom types...
  %{$This->{EStateAtomTypes}} = ();

  # Vector values precision for real values during E-state indicies...
  $This->{ValuesPrecision} = 3;

  # Calculated E-state values and indicies for generating E-state indicies fingerprints...
  %{$This->{EStateValues}} = ();
  %{$This->{EStateIndicies}} = ();
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

}

# Initialize object properties....
sub _InitializeEStateIndiciesFingerprintsProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  # Make sure molecule object was specified...
  if (!exists $NamesAndValues{Molecule}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying molecule...";
  }

  $This->_InitializeEstateIndicies();

  return $This;
}

# Initialize E-state indicies...
#
sub _InitializeEstateIndicies {
  my($This) = @_;

  # Set default EStateAtomTypesSetToUse...
  if (!$This->{EStateAtomTypesSetToUse}) {
    $This->{EStateAtomTypesSetToUse} = 'ArbitrarySize';
  }

  # Vector type...
  $This->{VectorType} = 'FingerprintsVector';

  if ($This->{EStateAtomTypesSetToUse} =~ /^FixedSize$/i) {
    $This->{FingerprintsVectorType} = 'OrderedNumericalValues';
  }
  else {
    $This->{FingerprintsVectorType} = 'NumericalValues';
  }

  $This->_InitializeFingerprintsVector();

  return $This;
}

# Disable set size method...
#
sub SetSize {
  my($This, $Type) = @_;

  croak "Error: ${ClassName}->SetSize: Can't change size:  It's not allowed...";
}

# Set E-state atom types set to use...
#
sub SetEStateAtomTypesSetToUse {
  my($This, $Value) = @_;

  if ($This->{EStateAtomTypesSetToUse}) {
    croak "Error: ${ClassName}->SetEStateAtomTypesSetToUse: Can't change size:  It's already set...";
  }

  if ($Value !~ /^(ArbitrarySize|FixedSize)/i) {
    croak "Error: ${ClassName}->SetEStateAtomTypesSetToUse: Unknown EStateAtomTypesSetToUse value: $Value; Supported values: ArbitrarySize or FixedSize";
  }

  $This->{EStateAtomTypesSetToUse} = $Value;

  return $This;
}

# Set vector values precision for real values for E-state indicies...
#
sub SetValuesPrecision {
  my($This, $Value) = @_;

  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetValuesPrecision: ValuesPrecision value, $Value, is not valid:  It must be a positive integer...";
  }
  $This->{ValuesPrecision} = $Value;

  return $This;
}

# Generate fingerprints description...
#
sub GetDescription {
  my($This) = @_;

  # Is description explicity set?
  if (exists $This->{Description}) {
    return $This->{Description};
  }

  # Generate fingerprints description...

  return "$This->{Type}:$This->{EStateAtomTypesSetToUse}";
}

# Generate electrotopological state indicies (E-state) [ Ref 75-78 ] fingerprints for
# non-hydrogen atoms in a molecule...
#
# EStateIndicies fingerprints constitute a vector containing sum of E-state values
# for E-state atom types. Two types of E-state atom types set size are allowed:
#
# ArbitrarySize - Corrresponds to only E-state atom types detected in molecule
# FixedSize - Corresponds to fixed number of E-state atom types previously defined
#
# Module AtomTypes::EStateAtomTypes.pm is used to assign E-state atom types to
# non-hydrogen atoms in the molecule which is able to assign atom types to any valid
# atom group. However, for FixedSize value of EStateAtomTypesSetToUse, only a fixed
# set of E-state atom types corresponding to specific atom groups [ Appendix III in
# Ref 77 ] are used for fingerprints.
#
# The fixed size E-state atom type set size used during generation of fingerprints corresponding
# FixedSize value of EStateAtomTypesSetToUse contains 87 E-state non-hydrogen atom types
# in EStateAtomTypes.csv data file distributed with MayaChemTools.
#
# Combination of Type and EStateAtomTypesSetToUse allow generation of 2 different types of
# E-state indicies fingerprints:
#
# Type                        EStateAtomTypesSetToUse
#
# EStateIndicies               ArbitrarySize      [ default fingerprints ]
# EStateIndicies               FixedSize
#
# The default is generate EStateIndicies type fingeprints corresponding to ArbitrarySize as
# EStateAtomTypesSetToUse value.
#
#
sub GenerateFingerprints {
  my($This) = @_;

  # Cache appropriate molecule data...
  $This->_SetupMoleculeDataCache();

  # Assign E-state atom types...
  if (!$This->_AssignEStateAtomTypes()) {
    carp "Warning: ${ClassName}->GenerateFingerprints: $This->{Type} fingerprints generation didn't succeed: Couldn't assign valid E-state atom types to all atoms...";
    return $This;
  }

  # Calculate E-state indicies...
  if (!$This->_CalculateEStateIndicies()) {
    carp "Warning: ${ClassName}->GenerateFingerprints: $This->{Type} fingerprints generation didn't succeed: Couldn't calculate E-state values for all atoms...";
    return $This;
  }

  # Set final fingerprints...
  $This->_SetFinalFingerprints();

  # Clear cached molecule data...
  $This->_ClearMoleculeDataCache();

  return $This;
}

# Assign E-state atom types...
#
sub _AssignEStateAtomTypes {
  my($This) = @_;
  my($EStateAtomTypes, $Atom, $AtomID, $AtomType);

  %{$This->{EStateAtomTypes}} = ();

  # Assign E-state atom types...
  $EStateAtomTypes = new AtomTypes::EStateAtomTypes('Molecule' => $This->{Molecule}, 'IgnoreHydrogens' => 1);
  $EStateAtomTypes->AssignAtomTypes();

  # Make sure atom types assignment is successful...
  if (!$EStateAtomTypes->IsAtomTypesAssignmentSuccessful()) {
    return undef;
  }

  # Collect assigned atom types...
  for $Atom (@{$This->{Atoms}}) {
    $AtomID = $Atom->GetID();

    $AtomType = $EStateAtomTypes->GetAtomType($Atom);
    $This->{EStateAtomTypes}{$AtomID} = $AtomType;
  }
  return $This;
}

# Calculate E-state indicies by summing up E-state values for specific
# E-state atom types...
#
sub _CalculateEStateIndicies {
  my($This) = @_;
  my($Atom, $AtomID, $AtomType, $EStateValue);

  # Calculate E-state values to generate E-state indicies...
  if (!$This->_CalculateEStateValuesDescriptors()) {
    return undef;
  }

  # Calculate E-state indicies...
  for $Atom (@{$This->{Atoms}}) {
    $AtomID = $Atom->GetID();

    $AtomType = $This->{EStateAtomTypes}{$AtomID};
    $EStateValue = $This->{EStateValues}{$AtomID};

    if (!exists $This->{EStateIndicies}{$AtomType}) {
      $This->{EStateIndicies}{$AtomType} = 0;
    }

    $This->{EStateIndicies}{$AtomType} += $EStateValue;
  }
  return $This;
}

# Calculate E-state values for E-state indicies...
#
sub _CalculateEStateValuesDescriptors {
  my($This) = @_;
  my($EStateValuesDescriptors, $Atom, $AtomID, $EStateValue);

  %{$This->{EStateValues}} = ();

  # Calculate and assign E-state values...
  $EStateValuesDescriptors = new AtomicDescriptors::EStateValuesDescriptors('Molecule' => $This->{Molecule});
  $EStateValuesDescriptors->GenerateDescriptors();

  # Make sure E-state values calculation is successful...
  if (!$EStateValuesDescriptors->IsDescriptorsGenerationSuccessful()) {
    return undef;
  }

  # Collect assigned E-state values...
  for $Atom (@{$This->{Atoms}}) {
    $AtomID = $Atom->GetID();
    $EStateValue = $EStateValuesDescriptors->GetDescriptorValue($Atom);
    $This->{EStateValues}{$AtomID} = $EStateValue;
  }
  return $This;
}

# Set final final fingerpritns for E-state indicies...
#
sub _SetFinalFingerprints {
  my($This) = @_;
  my($AtomType, $ValuesPrecision, $EStateAtomTypesDataRef, @Values, @IDs);

  # Mark successful generation of fingerprints...
  $This->{FingerprintsGenerated} = 1;

  @Values = ();
  @IDs = ();

  $ValuesPrecision = $This->{ValuesPrecision};

  if ($This->{EStateAtomTypesSetToUse} =~ /^FixedSize$/i) {
    # Use fixed size E-state atom types set for non-hydrogen atoms...
    for $AtomType (@{AtomTypes::EStateAtomTypes::GetAllPossibleEStateNonHydrogenAtomTypes()}) {
      push @IDs, "S${AtomType}";
      push @Values, exists($This->{EStateIndicies}{$AtomType}) ? MathUtil::round($This->{EStateIndicies}{$AtomType}, $ValuesPrecision) : 0;
    }
  }
  else {
    for $AtomType (sort keys %{$This->{EStateIndicies}}) {
      push @IDs, "S${AtomType}";
      push @Values, MathUtil::round($This->{EStateIndicies}{$AtomType}, $ValuesPrecision);
    }
  }

  # Add IDs and values to fingerprint vector...
  if (@IDs) {
    $This->{FingerprintsVector}->AddValueIDs(\@IDs);
  }
  $This->{FingerprintsVector}->AddValues(\@Values);

  return $This;
}

# Cache  appropriate molecule data...
#
sub _SetupMoleculeDataCache {
  my($This) = @_;

  # Get all non-hydrogen atoms...
  my($NegateAtomCheckMethod);
  $NegateAtomCheckMethod = 1;
  @{$This->{Atoms}} = $This->GetMolecule()->GetAtoms("IsHydrogen", $NegateAtomCheckMethod);

  return $This;
}

# Clear cached molecule data...
#
sub _ClearMoleculeDataCache {
  my($This) = @_;

  @{$This->{Atoms}} = ();

  return $This;
}

# Return a string containg data for EStateIndiciesFingerprints object...
sub StringifyEStateIndiciesFingerprints {
  my($This) = @_;
  my($EStateIndiciesFingerprintsString);

  # Type of Keys...
  $EStateIndiciesFingerprintsString = "Type: $This->{Type}; EStateAtomTypesSetToUse: $This->{EStateAtomTypesSetToUse}";

  # Fingerprint vector...
  $EStateIndiciesFingerprintsString .= "; FingerprintsVector: < $This->{FingerprintsVector} >";

  return $EStateIndiciesFingerprintsString;
}

1;

__END__

=head1 NAME

EStateIndiciesFingerprints

=head1 SYNOPSIS

use Fingerprints::EStateIndiciesFingerprints;

use Fingerprints::EStateIndiciesFingerprints qw(:all);

=head1 DESCRIPTION

B<EStateIndiciesFingerprints> [ Ref 75-78 ] class provides the following methods:

new, GenerateFingerprints, GetDescription, SetEStateAtomTypesSetToUse,
SetValuesPrecision, StringifyEStateIndiciesFingerprints

B<EStateIndiciesFingerprints> is derived from B<Fingerprints> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<AtomNeighborhoodsFingerprints>, B<Fingerprints> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

E-state atom types are assigned to all non-hydrogen atoms in a molecule using module
AtomTypes::EStateAtomTypes.pm and E-state values are calculated using module
AtomicDescriptors::EStateValues.pm. Using E-state atom types and E-state values,
B<EStateIndiciesFingerprints> constituting sum of E-state values for E-sate atom types
are generated.

Two types of E-state atom types set size are allowed:

    ArbitrarySize - Corresponds to only E-state atom types detected
                    in molecule
    FixedSize - Corresponds to fixed number of E-state atom types previously
                defined

Module AtomTypes::EStateAtomTypes.pm, used to assign E-state atom types to
non-hydrogen atoms in the molecule, is able to assign atom types to any valid
atom group. However, for I<FixedSize> value of B<EStateAtomTypesSetToUse>, only a
fixed set of E-state atom types corresponding to specific atom groups [ Appendix III in
Ref 77 ] are used for fingerprints.

The fixed size E-state atom type set size used during generation of fingerprints contains
87 E-state non-hydrogen atom types in EStateAtomTypes.csv data file distributed with
MayaChemTools.

Combination of Type and EStateAtomTypesSetToUse allow generation of 2 different types of
E-state indicies fingerprints:

    Type                        EStateAtomTypesSetToUse

    EStateIndicies              ArbitrarySize      [ default fingerprints ]
    EStateIndicies              FixedSize

The current release of MayaChemTools generates the following types of E-state
fingerprints vector strings:

    FingerprintsVector;EStateIndicies:ArbitrarySize;11;NumericalValues;IDs
    AndValuesString;SaaCH SaasC SaasN SdO SdssC SsCH3 SsF SsOH SssCH2 SssN
    H SsssCH;24.778 4.387 1.993 25.023 -1.435 3.975 14.006 29.759 -0.073 3
    .024 -2.270

    FingerprintsVector;EStateIndicies:FixedSize;87;OrderedNumericalValues;
    ValuesString;0 0 0 0 0 0 0 3.975 0 -0.073 0 0 24.778 -2.270 0 0 -1.435
    4.387 0 0 0 0 0 0 3.024 0 0 0 0 0 0 0 1.993 0 29.759 25.023 0 0 0 0 1
    4.006 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0

    FingerprintsVector;EStateIndicies:FixedSize;87;OrderedNumericalValues;
    IDsAndValuesString;SsLi SssBe SssssBem SsBH2 SssBH SsssB SssssBm SsCH3
    SdCH2 SssCH2 StCH SdsCH SaaCH SsssCH SddC StsC SdssC SaasC SaaaC Sssss
    C SsNH3p SsNH2 SssNH2p SdNH SssNH SaaNH StN SsssNHp SdsN SaaN SsssN Sd
    0 0 0 0 0 0 0 3.975 0 -0.073 0 0 24.778 -2.270 0 0 -1.435 4.387 0 0 0
    0 0 0 3.024 0 0 0 0 0 0 0 1.993 0 29.759 25.023 0 0 0 0 14.006 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0...

=head2 METHODS

=over 4

=item B<new>

    $EStateIndiciesFingerprints = new EStateIndiciesFingerprints(%NamesAndValues);

Using specified I<EStateIndiciesFingerprints> property names and values hash, B<new> method creates a new object
and returns a reference to newly created B<PathLengthFingerprints> object. By default, the
following properties are initialized:

    Molecule = '';
    Type = 'EStateIndicies'
    EStateAtomTypesSetToUse = 'ArbitrarySize'
    ValuesPrecision = 3

Examples:

    $EStateIndiciesFingerprints = new AtomTypesFingerprints(
                              'Molecule' => $Molecule,
                              'EStateAtomTypesSetToUse' =>
                                              'ArbitrarySize');

    $EStateIndiciesFingerprints = new AtomTypesFingerprints(
                              'Molecule' => $Molecule,
                              'EStateAtomTypesSetToUse' =>
                                              'FixedSize');

    $EStateIndiciesFingerprints->GenerateFingerprints();
    print "$EStateIndiciesFingerprints\n";

=item B<GenerateFingerprints>

    $EStateIndiciesFingerprints = $EStateIndiciesFingerprints->
                                  GenerateEStateIndiciesFingerprints();

Generates EState keys fingerprints and returns I<EStateIndiciesFingerprints>.

=item B<GetDescription>

    $Description = $EStateIndiciesFingerprints->GetDescription();

Returns a string containing description of EState keys fingerprints.

=item B<SetEStateAtomTypesSetToUse>

    $EStateIndiciesFingerprints->SetEStateAtomTypesSetToUse($Value);

Sets I<Value> of I<EStateAtomTypesSetToUse> and returns I<EStateIndiciesFingerprints>.
Possible values: I<ArbitrarySize or FixedSize>. Default value: I<ArbitrarySize>.

=item B<SetValuesPrecision>

    $EStateIndiciesFingerprints->SetValuesPrecision($Precision);

Sets precesion of E-state values to use during generation of E-state indices fingerprints
and returns I<EStateIndiciesFingerprints>. Possible values: I<Positive integers>.
Default value: I<3>.

=item B<StringifyEStateIndiciesFingerprints>

    $String = $EStateIndiciesFingerprints->StringifyEStateIndiciesFingerprints();

Returns a string containing information about I<EStateIndiciesFingerprints> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Fingerprints.pm, FingerprintsStringUtil.pm, AtomNeighborhoodsFingerprints.pm, AtomTypesFingerprints.pm,
ExtendedConnectivityFingerprints.pm, MACCSKeys.pm, PathLengthFingerprints.pm,
TopologicalAtomPairsFingerprints.pm, TopologicalAtomTripletsFingerprints.pm,
TopologicalAtomTorsionsFingerprints.pm, TopologicalPharmacophoreAtomPairsFingerprints.pm,
TopologicalPharmacophoreAtomTripletsFingerprints.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
