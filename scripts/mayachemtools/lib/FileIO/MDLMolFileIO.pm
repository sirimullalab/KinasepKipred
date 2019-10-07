package FileIO::MDLMolFileIO;
#
# File: MDLMolFileIO.pm
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
use FileUtil ();
use SDFileUtil ();
use FileIO::FileIO;
use Molecule;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(FileIO::FileIO Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(IsMDLMolFile);

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeMDLMolFileIO();

  $This->_InitializeMDLMolFileIOProperties(%NamesAndValues);

  return $This;
}

# Initialize any local object data...
#
sub _InitializeMDLMolFileIO {
  my($This) = @_;

  # Nothing to do: Base class FileIO handles default class variables...

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

}

# Initialize object values...
sub _InitializeMDLMolFileIOProperties {
  my($This, %NamesAndValues) = @_;

  # All other property names and values along with all Set/Get<PropertyName> methods
  # are implemented on-demand using ObjectProperty class.

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  if (!exists $NamesAndValues{Name}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying file name...";
  }

  # Make sure it's a MDLMol file...
  $Name = $NamesAndValues{Name};
  if (!$This->IsMDLMolFile($Name)) {
    croak "Error: ${ClassName}->New: Object can't be instantiated: File, $Name, doesn't appear to be MDLMol format...";
  }

  return $This;
}

# Is it a MDLMol file?
sub IsMDLMolFile ($;$) {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $FileName, $Status);

  if ((@_ == 2) && (_IsMDLMolFileIO($FirstParameter))) {
    ($This, $FileName) = ($FirstParameter, $SecondParameter);
  }
  else {
    $FileName = $FirstParameter;
  }

  # Check file extension...
  $Status = FileUtil::CheckFileType($FileName, "mol");

  return $Status;
}

# Read molecule from file and return molecule object...
sub ReadMolecule {
  my($This) = @_;
  my($FileHandle);

  $FileHandle = $This->GetFileHandle();
  return $This->ParseMoleculeString(SDFileUtil::ReadCmpdString($FileHandle));
}

# Write compound data using Molecule object...
sub WriteMolecule {
  my($This, $Molecule) = @_;

  if (!(defined($Molecule) && $Molecule->IsMolecule())) {
    carp "Warning: ${ClassName}->WriteMolecule: No data written: Molecule object is not specified...";
    return $This;
  }
  my($FileHandle);
  $FileHandle = $This->GetFileHandle();

  print $FileHandle $This->GenerateMoleculeString($Molecule) . "\n";

  return $This;
}

# Retrieve molecule string...
sub ReadMoleculeString {
  my($This) = @_;
  my($FileHandle);

  $FileHandle = $This->GetFileHandle();
  return SDFileUtil::ReadCmpdString($FileHandle);
}

# Parse molecule string and return molecule object. ParseMoleculeString supports two invocation methods: class
# method or a package function.
#
sub ParseMoleculeString {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $MoleculeString);

  if ((@_ == 2) && (_IsMDLMolFileIO($FirstParameter))) {
    ($This, $MoleculeString) = ($FirstParameter, $SecondParameter);
  }
  else {
    $MoleculeString = $FirstParameter;
    $This = undef;
  }
  if (!$MoleculeString) {
    return undef;
  }
  my($LineIndex, @MoleculeLines);
  @MoleculeLines = split /\n/, $MoleculeString;

  # Create molecule object and set molecule level native and MDL properties...
  #
  my($Molecule);
  $Molecule = new Molecule();

  # Set valence model for calculating implicit hydrogens...
  $Molecule->SetValenceModel('MDLValenceModel');

  # Process headers data...
  $LineIndex = 0;
  my($MoleculeName) = SDFileUtil::ParseCmpdMolNameLine($MoleculeLines[$LineIndex]);
  $MoleculeName = TextUtil::RemoveTrailingWhiteSpaces($MoleculeName);
  $Molecule->SetName($MoleculeName);

  $LineIndex++;
  my($UserInitial, $ProgramName, $Date, $Code, $ScalingFactor1, $ScalingFactor2, $Energy, $RegistryNum) = SDFileUtil::ParseCmpdMiscInfoLine($MoleculeLines[$LineIndex]);
  $Molecule->SetProperties('MDLUserInitial' => $UserInitial, 'MDLProgramName' => $ProgramName, 'MDLDate' => $Date, 'MDLCode' => $Code, 'MDLScalingFactor1' => $ScalingFactor1, 'MDLScalingFactor2' => $ScalingFactor2, 'MDLEnergy' => $Energy, 'MDLRegistryNum' => $RegistryNum);

  $LineIndex++;
  my($Comments) = SDFileUtil::ParseCmpdCommentsLine($MoleculeLines[$LineIndex]);
  $Molecule->SetProperties('MDLComments' => $Comments);

  $LineIndex++;
  my($AtomCount, $BondCount, $ChiralFlag, $PropertyCount, $Version) = SDFileUtil::ParseCmpdCountsLine($MoleculeLines[$LineIndex]);

  $Molecule->SetProperties('MDLChiralFlag' => $ChiralFlag, 'MDLPropertyCount' => $PropertyCount, 'MDLVersion' => $Version);

  # Process atom data...
  my($FirstAtomLineIndex, $LastAtomLineIndex, $AtomNum, $AtomX, $AtomY, $AtomZ, $AtomSymbol, $MassDifference, $Charge, $StereoParity, $Atom, %AtomNumToAtomMap);

  $AtomNum = 0;
  %AtomNumToAtomMap = ();
  $FirstAtomLineIndex = 4; $LastAtomLineIndex = $FirstAtomLineIndex + $AtomCount - 1;

  for ($LineIndex = $FirstAtomLineIndex; $LineIndex <= $LastAtomLineIndex; $LineIndex++) {
    $AtomNum++;
    ($AtomSymbol, $AtomX, $AtomY, $AtomZ, $MassDifference, $Charge, $StereoParity) = SDFileUtil::ParseCmpdAtomLine($MoleculeLines[$LineIndex]);

    $Atom = new Atom('AtomSymbol' => $AtomSymbol, 'XYZ' => [$AtomX, $AtomY, $AtomZ]);

    if ($MassDifference && $MassDifference != 0) {
      _ProcessMassDifference($Atom, $MassDifference);
    }
    if ($Charge && $Charge != 0) {
      _ProcessCharge($Atom, $Charge);
    }
    if ($StereoParity && $StereoParity != 0) {
      _ProcessStereoParity($Atom, $StereoParity);
    }

    $AtomNumToAtomMap{$AtomNum} = $Atom;
    $Molecule->AddAtom($Atom);
  }

  # Process bond data...
  my($FirstBondLineIndex, $LastBondLineIndex, $FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo, $InternalBondOrder, $InternalBondType, $Bond, $Atom1, $Atom2);

  $FirstBondLineIndex = $FirstAtomLineIndex + $AtomCount;
  $LastBondLineIndex = $FirstAtomLineIndex + $AtomCount + $BondCount - 1;

  for ($LineIndex = $FirstBondLineIndex; $LineIndex <= $LastBondLineIndex; $LineIndex++) {
    ($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo) = SDFileUtil::ParseCmpdBondLine($MoleculeLines[$LineIndex]);

    $Atom1 = $AtomNumToAtomMap{$FirstAtomNum};
    $Atom2 = $AtomNumToAtomMap{$SecondAtomNum};

    ($InternalBondOrder, $InternalBondType) = SDFileUtil::MDLBondTypeToInternalBondOrder($BondType);
    $Bond = new Bond('Atoms' => [$Atom1, $Atom2], 'BondOrder' => $InternalBondOrder);
    $Bond->SetBondType($InternalBondType);

    if ($BondStereo && $BondStereo != 0) {
      _ProcessBondStereo($Bond, $BondStereo);
    }

    $Molecule->AddBond($Bond);
  }

  # Process available property block lines starting with A  aaa, M CHG, M ISO and M RAD. All other property blocks
  # lines are for query or specific display purposes and are ignored for now.
  #
  #
  my($PropertyLineIndex, $PropertyLine, $FirstChargeOrRadicalLine, @ValuePairs);

  $PropertyLineIndex = $FirstAtomLineIndex + $AtomCount + $BondCount;
  $PropertyLine = $MoleculeLines[$PropertyLineIndex];
  $FirstChargeOrRadicalLine = 1;

  PROPERTYLINE: while ($PropertyLine !~ /^M  END/i ) {
    if ($PropertyLine =~ /\$\$\$\$/) {
      last PROPERTYLINE;
    }
    if ($PropertyLine =~ /^(M  CHG|M  RAD)/i) {
      if ($FirstChargeOrRadicalLine) {
	$FirstChargeOrRadicalLine = 0;
	_ZeroOutAtomsChargeAndRadicalValues(\%AtomNumToAtomMap);
      }
      if ($PropertyLine =~ /^M  CHG/i) {
	@ValuePairs = SDFileUtil::ParseCmpdChargePropertyLine($PropertyLine);
	_ProcessChargeProperty(\@ValuePairs, \%AtomNumToAtomMap);
      }
      elsif ($PropertyLine =~ /^M  RAD/i) {
	@ValuePairs = SDFileUtil::ParseCmpdRadicalPropertyLine($PropertyLine);
	_ProcessRadicalProperty(\@ValuePairs, \%AtomNumToAtomMap);
      }
    }
    elsif ($PropertyLine =~ /^M  ISO/i) {
      @ValuePairs = SDFileUtil::ParseCmpdIsotopePropertyLine($PropertyLine);
      _ProcessIsotopeProperty(\@ValuePairs, \%AtomNumToAtomMap);
    }
    elsif ($PropertyLine =~ /^A  /i) {
      my($NextPropertyLine);
      $PropertyLineIndex++;
      $NextPropertyLine = $MoleculeLines[$PropertyLineIndex];
      @ValuePairs = SDFileUtil::ParseCmpdAtomAliasPropertyLine($PropertyLine, $NextPropertyLine);
      _ProcessAtomAliasProperty(\@ValuePairs, \%AtomNumToAtomMap);
    }
    $PropertyLineIndex++;
    $PropertyLine = $MoleculeLines[$PropertyLineIndex];
  }
  # Store input molecule string as generic property of molecule...
  $Molecule->SetInputMoleculeString($MoleculeString);

  return $Molecule;
}

# Generate molecule string using molecule object...
sub GenerateMoleculeString {
  my($FirstParameter, $SecondParameter) = @_;
  my($This, $Molecule);

  if ((@_ == 2) && (_IsMDLMolFileIO($FirstParameter))) {
    ($This, $Molecule) = ($FirstParameter, $SecondParameter);
  }
  else {
    $Molecule = $FirstParameter;
    $This = undef;
  }
  if (!defined($Molecule)) {
    return undef;
  }
  my(@MoleculeLines);
  @MoleculeLines = ();

  # First line: Molname line...
  push @MoleculeLines, SDFileUtil::GenerateCmpdMolNameLine($Molecule->GetName());

  # Second line: Misc info...
  my($ProgramName, $UserInitial, $Code);
  $ProgramName = ''; $UserInitial = ''; $Code = '';

  $Code = $Molecule->IsThreeDimensional() ? '3D' : '2D';

  push @MoleculeLines, SDFileUtil::GenerateCmpdMiscInfoLine($ProgramName, $UserInitial, $Code);

  # Third line: Comments line...
  my($Comments);
  $Comments = $Molecule->HasProperty('MDLComments') ? $Molecule->GetMDLComments() : ($Molecule->HasProperty('Comments') ? $Molecule->GetComments() : '');
  push @MoleculeLines, SDFileUtil::GenerateCmpdCommentsLine($Comments);

  # Fourth line: Counts line for V2000
  my($AtomCount, $BondCount, $ChiralFlag);
  $AtomCount = $Molecule->GetNumOfAtoms();
  $BondCount = $Molecule->GetNumOfBonds();
  $ChiralFlag = 0;
  push @MoleculeLines, SDFileUtil::GenerateCmpdCountsLine($AtomCount, $BondCount, $ChiralFlag);

  # Atom lines...
  my($Atom, $AtomSymbol, $AtomX, $AtomY, $AtomZ, $MassDifference, $Charge, $StereoParity, $AtomNum, $AtomID, @Atoms, %AtomIDToNum);
  my($ChargePropertyValue, $IsotopePropertyValue, $RadicalPropertyValue, $AtomAliasPropertyValue, @IsotopePropertyValuePairs, @ChargePropertyValuePairs, @RadicalPropertyValuePairs, @AtomAliasPropertyValuePairs);

  @ChargePropertyValuePairs = ();
  @IsotopePropertyValuePairs = ();
  @RadicalPropertyValuePairs = ();
  @AtomAliasPropertyValuePairs = ();

  @Atoms = $Molecule->GetAtoms();

  $AtomNum = 0;
  for $Atom (@Atoms) {
    $AtomNum++;
    $AtomID = $Atom->GetID();
    $AtomIDToNum{$AtomID} = $AtomNum;

    $AtomSymbol = $Atom->GetAtomSymbol();
    ($AtomX, $AtomY, $AtomZ) = $Atom->GetXYZ();

    # Setup mass difference...
    $MassDifference = _GetMassDifference($Atom);
    if ($MassDifference) {
      # Hold it for M  ISO property lines...
      $IsotopePropertyValue = _GetIsotopePropertyValue($Atom);
      if ($IsotopePropertyValue) {
	push @IsotopePropertyValuePairs, ($AtomNum, $IsotopePropertyValue);
      }
    }

    # Setup charge...
    $Charge = _GetCharge($Atom);
    if ($Charge) {
      # Hold it for M  CHG property lines...
      $ChargePropertyValue = _GetChargePropertyValue($Atom);
      if ($ChargePropertyValue) {
	push @ChargePropertyValuePairs, ($AtomNum, $ChargePropertyValue);
      }
    }

    # Hold any radical values for  for M  RAD property lines...
    $RadicalPropertyValue = _GetRadicalPropertyValue($Atom);
    if ($RadicalPropertyValue) {
      push @RadicalPropertyValuePairs, ($AtomNum, $RadicalPropertyValue);
    }

    # Hold any atom alias value for A  xxx property lines....
    $AtomAliasPropertyValue = _GetAtomAliasPropertyValue($Atom);
    if ($AtomAliasPropertyValue) {
      push @AtomAliasPropertyValuePairs, ($AtomNum, $AtomAliasPropertyValue);

      # Set AtomSymbol to carbon as atom alias would override its value during parsing...
      $AtomSymbol = "C";
    }

    # Setup stereo parity...
    $StereoParity = _GetStereoParity($Atom);

    push @MoleculeLines, SDFileUtil::GenerateCmpdAtomLine($AtomSymbol, $AtomX, $AtomY, $AtomZ, $MassDifference, $Charge, $StereoParity);
  }

  # Bond lines...
  my($FirstAtomID, $FirstAtom, $FirstAtomNum, $SecondAtomID, $SecondAtom, $SecondAtomNum, $MDLBondType, $BondOrder, $BondType, $MDLBondStereo, $Bond, @Bonds);
  for $FirstAtom (@Atoms) {
    $FirstAtomID = $FirstAtom->GetID();
    $FirstAtomNum = $AtomIDToNum{$FirstAtomID};

    @Bonds = ();
    @Bonds = $FirstAtom->GetBonds();
    BOND: for $Bond (@Bonds) {
      $SecondAtom = $Bond->GetBondedAtom($FirstAtom);
      $SecondAtomID = $SecondAtom->GetID();
      $SecondAtomNum = $AtomIDToNum{$SecondAtomID};
      if ($FirstAtomNum >= $SecondAtomNum) {
	next BOND;
      }
      # Setup BondType...
      $BondOrder = $Bond->GetBondOrder();
      $BondType = $Bond->GetBondType();
      $MDLBondType = SDFileUtil::InternalBondOrderToMDLBondType($BondOrder, $BondType);

      # Setup BondStereo...
      $MDLBondStereo = _GetBondStereo($Bond);

      push @MoleculeLines, SDFileUtil::GenerateCmpdBondLine($FirstAtomNum, $SecondAtomNum, $MDLBondType, $MDLBondStereo);
    }
  }
  # Property lines...
  if (@IsotopePropertyValuePairs) {
    push @MoleculeLines, SDFileUtil::GenerateCmpdIsotopePropertyLines(\@IsotopePropertyValuePairs);
  }
  if (@ChargePropertyValuePairs) {
    push @MoleculeLines, SDFileUtil::GenerateCmpdChargePropertyLines(\@ChargePropertyValuePairs);
  }
  if (@RadicalPropertyValuePairs) {
    push @MoleculeLines, SDFileUtil::GenerateCmpdRadicalPropertyLines(\@RadicalPropertyValuePairs);
  }
  if (@AtomAliasPropertyValuePairs) {
    push @MoleculeLines, SDFileUtil::GenerateCmpdAtomAliasPropertyLines(\@AtomAliasPropertyValuePairs);
  }

  push @MoleculeLines, "M  END";

  return join "\n", @MoleculeLines;
}

# Process MassDifference value and set atom's mass number...
#
sub _ProcessMassDifference {
  my($Atom, $MassDifference) = @_;
  my($MassNumber, $NewMassNumber, $AtomicNumber);

  $AtomicNumber = $Atom->GetAtomicNumber();

  if (!$AtomicNumber) {
    carp "Warning: ${ClassName}->_ProcessMassDifference: Ignoring specified mass difference value, $MassDifference, in SD file: Assigned to non standard element...";
    return;
  }
  $MassNumber = $Atom->GetMassNumber();
  if (!$MassDifference) {
    carp "Warning: ${ClassName}->_ProcessMassDifference: Ignoring specified mass difference value, $MassDifference, in SD file: Unknown MassNumber value...";
    return;
  }
  $NewMassNumber = $MassNumber + $MassDifference;
  if (!PeriodicTable::IsElementNaturalIsotopeMassNumber($AtomicNumber, $NewMassNumber)) {
    my($AtomSymbol) = $Atom->GetAtomSymbol();
    carp "Warning: ${ClassName}->_ProcessMassDifference: Unknown mass number, $MassNumber, corresponding to specified mass difference value, $MassDifference, in SD for atom with atomic number, $AtomicNumber, and atomic symbol, $AtomSymbol. The mass number value has been assigned. Don't forget to Set ExactMass property explicitly; otherwise, GetExactMass method would return mass of most abundant isotope...\n";
  }

  # Use SetProperty method instead of SetMassNumber to skip explicit checks on MassNumber value...
  $Atom->SetProperty('MassNumber', $NewMassNumber);
}

# Get mass difference value...
sub _GetMassDifference {
  my($Atom) = @_;
  my($MassDifference, $MassNumber, $MostAbundantMassNumber, $AtomicNumber);

  $MassDifference = 0;
  $MassNumber = $Atom->GetMassNumber();
  if (defined $MassNumber) {
    $AtomicNumber = $Atom->GetAtomicNumber();
    if (defined $AtomicNumber) {
      $MostAbundantMassNumber = PeriodicTable::GetElementMostAbundantNaturalIsotopeMassNumber($AtomicNumber);
      if (defined($MostAbundantMassNumber) && $MassNumber != $MostAbundantMassNumber) {
	$MassDifference = $MassNumber - $MostAbundantMassNumber;
      }
    }
  }
  return $MassDifference;
}

# Process formal charge value and assign it to atom as formal charge...
sub _ProcessCharge {
  my($Atom, $Charge) = @_;
  my($InternalCharge);

  $InternalCharge = SDFileUtil::MDLChargeToInternalCharge($Charge);
  $Atom->SetFormalCharge($InternalCharge);
}

# Get MDL formal charge value ...
sub _GetCharge {
  my($Atom) = @_;
  my($InternalCharge, $Charge);

  $Charge = 0;
  if ($Atom->HasProperty('FormalCharge')) {
    $InternalCharge = $Atom->GetFormalCharge();
    if ($InternalCharge) {
      $Charge = SDFileUtil::InternalChargeToMDLCharge($InternalCharge);
    }
  }
  return $Charge;
}

# Process stereo parity value and assign it to atom as MDL property...
#
# Notes:
#   . Mark atom as chiral center
#   . Assign any explicit Clockwise (parity 1), CounterClockwise (parity 2) or either value (parity 3) as property of atom.
#   . MDL values of Clockwise and CounterClockwise don't correspond to priority assigned to ligands around
#     stereo center using CIP scheme; consequently, these values can't be used to set internal Stereochemistry for
#     an atom.
#
sub _ProcessStereoParity {
  my($Atom, $StereoParity) = @_;

  $Atom->SetStereoCenter('1');
  $Atom->SetMDLStereoParity($StereoParity);
}

# Set stereo parity value to zero for now: The current release of MayaChemTools hasn't implemented
# functionality to determine chirality.
#
sub _GetStereoParity {
  my($Atom) = @_;
  my($StereoParity);

  $StereoParity = 0;

  return $StereoParity;
}

# Process bond stereo value...
sub _ProcessBondStereo {
  my($Bond, $BondStereo) = @_;
  my($InternalBondStereo);

  $InternalBondStereo = SDFileUtil::MDLBondStereoToInternalBondStereochemistry($BondStereo);
  if ($InternalBondStereo) {
    $Bond->SetBondStereochemistry($InternalBondStereo);
  }
}

# Get MDLBondStereo value...
sub _GetBondStereo {
  my($Bond) = @_;
  my($InternalBondStereo, $BondStereo);

  $BondStereo = 0;

  $InternalBondStereo = '';
  BONDSTEREO: {
    if ($Bond->IsUp()) {
      $InternalBondStereo = 'Up';
      last BONDSTEREO;
    }
    if ($Bond->IsDown()) {
      $InternalBondStereo = 'Down';
      last BONDSTEREO;
    }
    if ($Bond->IsUpOrDown()) {
      $InternalBondStereo = 'UpOrDown';
      last BONDSTEREO;
    }
    if ($Bond->IsCisOrTrans() || $Bond->IsCis() || $Bond->IsTrans()) {
      $InternalBondStereo = 'CisOrTrans';
      last BONDSTEREO;
    }
    $InternalBondStereo = '';
  }

  if ($InternalBondStereo) {
    $BondStereo = SDFileUtil::InternalBondStereochemistryToMDLBondStereo($InternalBondStereo);
  }

  return $BondStereo;
}

# Zero out charge and radical values specified for atoms...
sub _ZeroOutAtomsChargeAndRadicalValues {
  my($AtomNumToAtomMapRef) = @_;
  my($Atom);

  for $Atom (values %{$AtomNumToAtomMapRef}) {
    if ($Atom->HasProperty('FormalCharge')) {
      $Atom->DeleteProperty('FormalCharge');
    }
    elsif ($Atom->HasProperty('SpinMultiplicity')) {
      $Atom->DeleteProperty('SpinMultiplicity');
    }
  }
}

# Process charge property value pairs...
sub _ProcessChargeProperty {
  my($ValuePairsRef, $AtomNumToAtomMapRef) = @_;

  if (!(defined($ValuePairsRef) && @{$ValuePairsRef})) {
    return;
  }
  my($Index, $ValuePairsCount, $AtomNum, $Charge, $Atom);

  $ValuePairsCount = scalar @{$ValuePairsRef};
  VALUEPAIRS: for ($Index = 0; $Index < $ValuePairsCount; $Index +=2) {
    $AtomNum = $ValuePairsRef->[$Index]; $Charge = $ValuePairsRef->[$Index + 1];
    if (!$Charge) {
      next VALUEPAIRS;
    }
    if (!exists $AtomNumToAtomMapRef->{$AtomNum}) {
      next VALUEPAIRS;
    }
    $Atom = $AtomNumToAtomMapRef->{$AtomNum};
    if ($Atom->HasProperty('SpinMultiplicity')) {
      carp "Warning: ${ClassName}->_ProcessChargeProperty: Setting formal charge on atom number, $AtomNum,  with already assigned spin multiplicity value...";
    }
    $Atom->SetFormalCharge($Charge);
  }
}

# Get charge property value for an atom...
sub _GetChargePropertyValue {
  my($Atom) = @_;
  my($Charge);

  $Charge = 0;
  if ($Atom->HasProperty('FormalCharge')) {
    $Charge = $Atom->GetFormalCharge();
  }
  return $Charge;
}

# Process charge property value pairs...
sub _ProcessRadicalProperty {
  my($ValuePairsRef, $AtomNumToAtomMapRef) = @_;

  if (!(defined($ValuePairsRef) && @{$ValuePairsRef})) {
    return;
  }
  my($Index, $ValuePairsCount, $AtomNum, $Radical, $SpinMultiplicity, $Atom);

  $ValuePairsCount = scalar @{$ValuePairsRef};
  VALUEPAIRS: for ($Index = 0; $Index < $ValuePairsCount; $Index +=2) {
    $AtomNum = $ValuePairsRef->[$Index]; $Radical = $ValuePairsRef->[$Index + 1];
    if (!$Radical) {
      next VALUEPAIRS;
    }
    if (!exists $AtomNumToAtomMapRef->{$AtomNum}) {
      next VALUEPAIRS;
    }
    $Atom = $AtomNumToAtomMapRef->{$AtomNum};
    if ($Atom->HasProperty('FormalCharge')) {
      carp "Warning: ${ClassName}->_ProcessRadicalProperty: Setting spin multiplicity on atom number, $AtomNum,  with already assigned formal charge value...";
    }
    $SpinMultiplicity = SDFileUtil::MDLRadicalToInternalSpinMultiplicity($Radical);
    $Atom->SetSpinMultiplicity($SpinMultiplicity);
  }
}

# Get radical property value for an atom...
sub _GetRadicalPropertyValue {
  my($Atom) = @_;
  my($Radical, $SpinMultiplicity);

  $Radical = 0;
  if ($Atom->HasProperty('SpinMultiplicity')) {
    $SpinMultiplicity = $Atom->GetSpinMultiplicity();
    $Radical = SDFileUtil::InternalSpinMultiplicityToMDLRadical($SpinMultiplicity);
  }
  return $Radical;
}

# Process isotope property value pairs...
sub _ProcessIsotopeProperty {
  my($ValuePairsRef, $AtomNumToAtomMapRef) = @_;

  if (!(defined($ValuePairsRef) && @{$ValuePairsRef})) {
    return;
  }
  my($Index, $ValuePairsCount, $AtomNum, $MassNumber, $Atom, $AtomicNumber);

  $ValuePairsCount = scalar @{$ValuePairsRef};
  VALUEPAIRS: for ($Index = 0; $Index < $ValuePairsCount; $Index +=2) {
    $AtomNum = $ValuePairsRef->[$Index]; $MassNumber = $ValuePairsRef->[$Index + 1];
    if (!$MassNumber) {
      next VALUEPAIRS;
    }
    if (!exists $AtomNumToAtomMapRef->{$AtomNum}) {
      next VALUEPAIRS;
    }
    $Atom = $AtomNumToAtomMapRef->{$AtomNum};
    $AtomicNumber = $Atom->GetAtomicNumber();

    if (!PeriodicTable::IsElementNaturalIsotopeMassNumber($AtomicNumber, $MassNumber)) {
      my($AtomSymbol) = $Atom->GetAtomSymbol();
      carp "Warning: ${ClassName}->_ProcessProcessIsotopeProperty: Unknown mass number, $MassNumber, specified on M  ISO property line for atom number, $AtomNum,  in SD for atom with atomic number, $AtomicNumber, and atomic symbol, $AtomSymbol. The mass number value has been assigned. Don't forget to Set ExactMass property explicitly; otherwise, GetExactMass method would return mass of most abundant isotope...\n";
    }

    # Use SetProperty method instead of SetMassNumber to skip explicit checks on MassNumber value...
    $Atom->SetProperty('MassNumber', $MassNumber);
  }
}

# Get isotope property value for an atom...
sub _GetIsotopePropertyValue {
  my($Atom) = @_;
  my($MassNumber);

  $MassNumber = 0;
  if ($Atom->HasProperty('MassNumber')) {
    $MassNumber = $Atom->GetMassNumber();
  }
  return $MassNumber;
}

# Process atom alias property value pairs...
sub _ProcessAtomAliasProperty {
  my($ValuePairsRef, $AtomNumToAtomMapRef) = @_;

  if (!(defined($ValuePairsRef) && @{$ValuePairsRef})) {
    return;
  }
  my($Index, $ValuePairsCount, $AtomNum, $AtomAlias, $Atom);

  $ValuePairsCount = scalar @{$ValuePairsRef};
  VALUEPAIRS: for ($Index = 0; $Index < $ValuePairsCount; $Index +=2) {
    $AtomNum = $ValuePairsRef->[$Index]; $AtomAlias = $ValuePairsRef->[$Index + 1];
    if (!$AtomNum) {
      next VALUEPAIRS;
    }
    if (!exists $AtomNumToAtomMapRef->{$AtomNum}) {
      next VALUEPAIRS;
    }
    $AtomAlias = TextUtil::RemoveLeadingAndTrailingWhiteSpaces($AtomAlias);
    if (TextUtil::IsEmpty($AtomAlias)) {
      carp("Warning: ${ClassName}->_ProcessAtomAliasProperty: Ignoring atom alias property line: No Atom alias value specified...");
      next VALUEPAIRS;
    }

    # Set atom symbol to atom alias which sets atomic number automatically...
    $Atom = $AtomNumToAtomMapRef->{$AtomNum};
    $Atom->SetAtomSymbol($AtomAlias);

    $Atom->SetProperty('AtomAlias', $AtomAlias);
  }
}

# Get atom alias property value for an atom...
sub _GetAtomAliasPropertyValue {
  my($Atom) = @_;
  my($AtomAlias);

  $AtomAlias = undef;
  if ($Atom->HasProperty('AtomAlias')) {
    $AtomAlias = $Atom->GetAtomAlias();
  }
  return $AtomAlias;
}

# Is it a MDLMolFileIO object?
sub _IsMDLMolFileIO {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}


1;

__END__

=head1 NAME

MDLMolFileIO

=head1 SYNOPSIS

use FileIO::MDLMolFileIO;

use FileIO::MDLMolFileIO qw(:all);

=head1 DESCRIPTION

B<MDLMolFIleIO> class provides the following methods:

new, GenerateMoleculeString, IsMDLMolFile, ParseMoleculeString, ReadMolecule,
ReadMoleculeString, WriteMolecule

The following methods can also be used as functions:

GenerateMoleculeString, IsMDLMolFile, ParseMoleculeString

Data specific to B<MDLMolFileIO> class not directly used by B<Molecule>, B<Atom> and
B<Bond> objects - data label/value pairs, atom SteroParity and so on - is associated to
and retrieved from appropriate objects using following methods:

    SetMDL<PropertyName>
    GetMDL<PropertyName>.

B<MDLMolFileIO> class is derived from I<FileIO> class and uses its methods to support
generic file related functionality.

=head2 METHODS

=over 4

=item B<new>

    $NewMDLMolFileIO = new FileIO::MDLMolFileIO(%NamesAndValues);

Using specified I<MDLMolFileIO> property names and values hash, B<new> method creates a new object
and returns a reference to newly created B<MDLMolFileIO> object.

=item B<GenerateMoleculeString>

    $MoleculeString = $MDLMolFileIO->GenerateMoleculeString($Molecule);
    $MoleculeString = FileIO::MDLMolFileIO::GenerateMoleculeString($Molecule);

Returns a B<MoleculeString> in MDLMol format corresponding to I<Molecule>.

=item B<IsMDLMolFile>

    $Status = $MDLMolFileIO->IsMDLMolFile($FileName);
    $Status = FileIO::MDLMolFileIO::IsMDLMolFile($FileName);

Returns 1 or 0 based on whether I<FileName> is a MDLMol file.

=item B<ParseMoleculeString>

    $Molecule = $MDLMolFileIO->ParseMoleculeString($MoleculeString);
    $Molecule = FileIO::MDLMolFileIO::ParseMoleculeString($MoleculeString);

Parses I<MoleculeString> and returns a B<Molecule> object.

=item B<ReadMolecule>

    $Molecule = $MDLMolFileIO->ReadMolecule($FileHandle);

Reads data for the compound in a file using already opened I<FileHandle>, creates,
and returns a B<Molecule> object.

=item B<ReadMoleculeString>

    $MoleculeString = $MDLMolFileIO->ReadMoleculeString($FileHandle);

Reads data for the compound in a file using already opened I<FileHandle> and
returns a B<MoleculeString> corresponding to compound structure and other associated
data.

=item B<WriteMolecule>

    $MDLMolFileIO->WriteMolecule($Molecule);

Writes I<Molecule> data to a file in MDLMol format and returns B<MDLMolFileIO>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

MoleculeFileIO.pm, SDFileIO.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
