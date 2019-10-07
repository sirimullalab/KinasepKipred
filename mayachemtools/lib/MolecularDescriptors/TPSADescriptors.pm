package MolecularDescriptors::TPSADescriptors;
#
# File: TPSADescriptors.pm
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
use TextUtil ();
use MathUtil ();
use Atom;
use Molecule;
use MolecularDescriptors::MolecularDescriptors;
use AtomTypes::TPSAAtomTypes;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(MolecularDescriptors::MolecularDescriptors Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GetDescriptorNames);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, @DescriptorNames);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyTPSADescriptors';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeTPSADescriptors();

  $This->_InitializeTPSADescriptorsProperties(%NamesAndValues);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Descriptor names...
  @DescriptorNames = ('TPSA');

}

# Get descriptor names as an array.
#
# This functionality can be either invoked as a class function or an
# object method.
#
sub GetDescriptorNames {
  return @DescriptorNames;
}

# Initialize object data...
#
sub _InitializeTPSADescriptors {
  my($This) = @_;

  # Type of MolecularDescriptor...
  $This->{Type} = 'TPSA';

  # By default, TPSA atom contributions from Phosphorus and Sulfur atoms
  # are not included during TPSA calculations. [ Ref 91 ]
  #
  $This->{IgnorePhosphorus} = 1;
  $This->{IgnoreSulfur} = 1;

  # TPSA atom types assigned to appropriate atoms...
  %{$This->{AtomTypes}} = ();

  # Intialize descriptor names and values...
  $This->_InitializeDescriptorNamesAndValues(@DescriptorNames);

  return $This;
}

# Initialize object properties...
#
sub _InitializeTPSADescriptorsProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  return $This;
}

# Calculate Topological Polar Surface Area (TPSA) value [ Ref 90-91 ] for molecule...
#
# Methodology:
#   . Assign TPSA atom types [ Ref 90-91 ] to Nitrogen and Oxygen
#     atoms with optional assignment to Phosphorus and Sulfur atoms.
#   . Calculate TPSA value adding contribution of appropriate atom types.
#
sub GenerateDescriptors {
  my($This) = @_;

  # Initialize descriptor values...
  $This->_InitializeDescriptorValues();

  # Check availability of molecule...
  if (!$This->{Molecule}) {
    carp "Warning: ${ClassName}->GenerateDescriptors: $This->{Type} molecular descriptors generation didn't succeed: Molecule data is not available: Molecule object hasn't been set...";
    return undef;
  }

  # Cache appropriate molecule data...
  $This->_SetupMoleculeDataCache();

  # Assign TPSA atom types...
  if (!$This->_AssignAtomTypes()) {
    carp "Warning: ${ClassName}->GenerateDescriptors: $This->{Type} molecular descriptors generation didn't succeed: Couldn't assign valid TPSA atom types to appropriate atoms...";
    return undef;
  }

  # Calculate descriptor values...
  if (!$This->_CalculateDescriptorValues()) {
    carp "Warning: ${ClassName}->GenerateDescriptors: $This->{Type} molecular descriptors generation didn't succeed: Couldn't calculate TPSA values corresponding to assigned TPSA atom types...";
    return undef;
  }

  # Set final descriptor values...
  $This->_SetFinalDescriptorValues();

  # Clear cached molecule data...
  $This->_ClearMoleculeDataCache();

  return $This;
}

# Assign TPSA atom types..
#
sub _AssignAtomTypes {
  my($This) = @_;
  my($TPSAAtomTypes, $Atom, $AtomID);

  %{$This->{AtomTypes}} = ();

  # Assign atom types...
  $TPSAAtomTypes = new AtomTypes::TPSAAtomTypes('Molecule' => $This->{Molecule}, 'IgnorePhosphorus' => $This->{IgnorePhosphorus}, 'IgnoreSulfur' => $This->{IgnoreSulfur});
  $TPSAAtomTypes->AssignAtomTypes();

  # Make sure TPSA atom types assignment is successful...
  if (!$TPSAAtomTypes->IsAtomTypesAssignmentSuccessful()) {
    return undef;
  }

  # Collect assigned atom types...
  for $Atom (@{$This->{Atoms}}) {
    $AtomID = $Atom->GetID();
    $This->{AtomTypes}{$AtomID} = $TPSAAtomTypes->GetAtomType($Atom);
  }

  return $This;
}

# Calculate TPSA value...
#
sub _CalculateDescriptorValues {
  my($This) = @_;
  my($Atom, $AtomID, $TPSA, $TPSAContribution, $TPSADataRef, $AtomType);

  $TPSA = 0;

  # Get TPSA atom types data...
  $TPSADataRef = AtomTypes::TPSAAtomTypes::GetTPSAAtomTypesData();

  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomID = $Atom->GetID();
    $AtomType = $This->{AtomTypes}{$AtomID};

    # Ignore inappropriate atoms...
    if ($AtomType =~ /^None$/i) {
      next ATOM;
    }

    $TPSAContribution = 0.0;

    if ($AtomType =~ /^(N|O)$/i) {
      # TPSA contributions for Nitrogen and Oxygen atoms not explicity defined using atom
      # environments in Table 1 [ Ref 90 ]
      if ($AtomType =~ /^N$/i) {
	# N = 30.5 - X*8.2 + H*1.5 or 0.0 for negative value
	$TPSAContribution = 30.5 - $Atom->GetAtomicInvariantValue('X') * 8.2 + $Atom->GetAtomicInvariantValue('H') * 1.5;
      }
      elsif ($AtomType =~ /^O$/i) {
	# O = 28.5 - X*8.6 + H*1.5 or 0.0 for negative value
	$TPSAContribution = 28.5 - $Atom->GetAtomicInvariantValue('X') * 8.6 + $Atom->GetAtomicInvariantValue('H') * 1.5;
      }
      if ($TPSAContribution < 0.0) {
	$TPSAContribution = 0.0;
      }
    }
    elsif (exists $TPSADataRef->{DataCol3}{$AtomType}) {
      # Data for TPSA contribution is in column number 3...
      $TPSAContribution = $TPSADataRef->{DataCol3}{$AtomType};
    }
    else {
      # No TPSA data for assigned atom type...
      return undef;
    }
    $TPSA += $TPSAContribution;
  }

  # Track the calculated values...
  $This->{TPSA} = MathUtil::round($TPSA, 2);

  return $This;
}

# Setup final descriptor values...
#
sub _SetFinalDescriptorValues {
  my($This) = @_;

  $This->{DescriptorsGenerated} = 1;

  $This->SetDescriptorValues($This->{TPSA});

  return $This;
}

# Cache  appropriate molecule data...
#
sub _SetupMoleculeDataCache {
  my($This) = @_;

  @{$This->{Atoms}} = $This->GetMolecule()->GetAtoms();

  return $This;
}

# Clear cached molecule data...
#
sub _ClearMoleculeDataCache {
  my($This) = @_;

  @{$This->{Atoms}} = ();

  return $This;
}

# Return a string containg data for TPSADescriptors object...
#
sub StringifyTPSADescriptors {
  my($This) = @_;
  my($TPSADescriptorsString);

  # Type of MolecularDescriptors...
  $TPSADescriptorsString = "MolecularDescriptorType: $This->{Type}; IgnorePhosphorus: " . ($This->{IgnorePhosphorus} ? "Yes" : "No") . "; IgnoreSulfur: " .  ($This->{IgnoreSulfur} ? "Yes" : "No");

  # Setup molecular descriptor information...
  $TPSADescriptorsString .= "; " . $This->_StringifyDescriptorNamesAndValues();

  return $TPSADescriptorsString;
}

# Is it a TPSADescriptors object?
sub _IsTPSADescriptors {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

TPSADescriptors

=head1 SYNOPSIS

use MolecularDescriptors::TPSADescriptors;

use MolecularDescriptors::TPSADescriptors qw(:all);

=head1 DESCRIPTION

B<TPSADescriptors> class provides the following methods:

new, GenerateDescriptors, GetDescriptorNames, StringifyTPSADescriptors

B<TPSADescriptors> is derived from B<MolecularDescriptors> class which in turn
is  derived from B<ObjectProperty> base class that provides methods not explicitly defined
in B<TPSADescriptors>, B<MolecularDescriptors> or B<ObjectProperty> classes using Perl's
AUTOLOAD functionality. These methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

After Topological Polar Surface Area (TPSA) atom types [ Ref 90-91 ] has been assigned
to appropriate atoms in a molecule using AtomTypes::TPSAAtomTypes.pm module,  TPSA value
is calculated by adding up contributions of each appropriate atom type.

By default, MayaChemTools only uses nitrogen and oxygen atoms during calculation
of TPSA and ignores phosphorous and sulfur atoms. [ Ref 90 - 91 ]

=head2 METHODS

=over 4

=item B<new>

    $NewTPSADescriptors = new MolecularDescriptors::TPSADescriptors(
                                                   %NamesAndValues);

Using specified I<TPSADescriptors> property names and values hash, B<new>
method creates a new object and returns a reference to newly created B<TPSADescriptors>
object. By default, the following properties are initialized:

    Molecule = ''
    Type = 'TPSA'
    IgnorePhosphorus = 1
    IgnoreSulfur = 1
    @DescriptorNames = ('TPSA')
    @DescriptorValues = ('None')

Examples:

    $TPSADescriptors = new MolecularDescriptors::TPSADescriptors(
                              'Molecule' => $Molecule);

    $TPSADescriptors = new MolecularDescriptors::TPSADescriptors();

    $TPSADescriptors = new MolecularDescriptors::TPSADescriptors(
                              'IgnorePhosphorus' => 0,
                              'IgnoreSulfur' => 0);

    $TPSADescriptors->SetMolecule($Molecule);
    $TPSADescriptors->GenerateDescriptors();
    print "TPSADescriptors: $TPSADescriptors\n";


=item B<GenerateDescriptors>

    $TPSADescriptors->GenerateDescriptors();

Calculate TPSA value for a molecule and returns I<TPSADescriptors>.

=item B<GetDescriptorNames>

    @DescriptorNames = $TPSADescriptors->GetDescriptorNames();
    @DescriptorNames = MolecularDescriptors::TPSADescriptors::
                       GetDescriptorNames();

Returns all available descriptor names as an array.

=item B<StringifyTPSADescriptors>

    $String = $TPSADescriptors->StringifyTPSADescriptors();

Returns a string containing information about I<TPSADescriptors> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

MolecularDescriptors.pm, MolecularDescriptorsGenerator.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
