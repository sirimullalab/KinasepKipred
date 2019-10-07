#!/usr/bin/perl -w
#
# File: CalculatePhysicochemicalProperties.pl
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
use SDFileUtil;
use MoleculeFileIO;
use Molecule;
use AtomTypes::AtomicInvariantsAtomTypes;
use AtomTypes::FunctionalClassAtomTypes;
use MolecularDescriptors::MolecularDescriptorsGenerator;

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

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Setup information about input files...
print "Checking input SD file(s)...\n";
my(%SDFilesInfo);
RetrieveSDFilesInfo();

# Process input files..
my($FileIndex);
if (@SDFilesList > 1) {
  print "\nProcessing SD files...\n";
}
for $FileIndex (0 .. $#SDFilesList) {
  if ($SDFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $SDFilesList[$FileIndex]...\n";
    CalculatePhysicochemicalProperties($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Calculate physicochemical properties for a SD file...
#
sub CalculatePhysicochemicalProperties {
  my($FileIndex) = @_;
  my($CmpdCount, $IgnoredCmpdCount, $RuleOf5ViolationsCount, $RuleOf3ViolationsCount, $SDFile, $MoleculeFileIO, $Molecule, $MolecularDescriptorsGenerator, $PhysicochemicalPropertiesDataRef, $NewSDFileRef, $NewTextFileRef);

  $SDFile = $SDFilesList[$FileIndex];

  # Setup output files...
  $NewSDFileRef = ''; $NewTextFileRef = '';
  ($NewSDFileRef, $NewTextFileRef) = SetupAndOpenOutputFiles($FileIndex);

  # Setup molecular descriptor generator to calculate property values for specifed
  # property names...
  $MolecularDescriptorsGenerator = SetupMolecularDescriptorsGenerator();

  ($CmpdCount, $IgnoredCmpdCount, $RuleOf5ViolationsCount, $RuleOf3ViolationsCount) = ('0') x 4;

  $MoleculeFileIO = new MoleculeFileIO('Name' => $SDFile);
  $MoleculeFileIO->Open();

  COMPOUND: while ($Molecule = $MoleculeFileIO->ReadMolecule()) {
    $CmpdCount++;

    # Filter compound data before calculating physiochemical properties...
    if ($OptionsInfo{Filter}) {
      if (CheckAndFilterCompound($CmpdCount, $Molecule)) {
	$IgnoredCmpdCount++;
	next COMPOUND;
      }
    }

    # Calculate properties...
    $PhysicochemicalPropertiesDataRef = CalculateMoleculeProperties($MolecularDescriptorsGenerator, $Molecule);

    if (!defined($PhysicochemicalPropertiesDataRef)) {
      $IgnoredCmpdCount++;
      ProcessIgnoredCompound('PropertiesCalculationFailed', $CmpdCount, $Molecule);
      next COMPOUND;
    }

    # Calculate any rule violations...
    if ($OptionsInfo{RuleOf5Violations} && $PhysicochemicalPropertiesDataRef->{RuleOf5Violations}) {
      $RuleOf5ViolationsCount++;
    }

    if ($OptionsInfo{RuleOf3Violations} && $PhysicochemicalPropertiesDataRef->{RuleOf3Violations}) {
      $RuleOf3ViolationsCount++;
    }

    # Write out calculate properties...
    WriteDataToOutputFiles($FileIndex, $CmpdCount, $Molecule, $PhysicochemicalPropertiesDataRef, $NewSDFileRef, $NewTextFileRef);
  }
  $MoleculeFileIO->Close();

  if ($OptionsInfo{SDOutput} && $NewSDFileRef) {
    close $NewSDFileRef;
  }
  if ($OptionsInfo{TextOutput} && $NewTextFileRef) {
    close $NewTextFileRef;
  }

  WriteCalculationSummaryStatistics($CmpdCount, $IgnoredCmpdCount, $RuleOf5ViolationsCount, $RuleOf3ViolationsCount);
}

# Process compound being ignored due to problems in physicochemical properties calculation...
#
sub ProcessIgnoredCompound {
  my($Mode, $CmpdCount, $Molecule) = @_;
  my($CmpdID, $DataFieldLabelAndValuesRef);

  $DataFieldLabelAndValuesRef = $Molecule->GetDataFieldLabelAndValues();
  $CmpdID = SetupCmpdIDForTextFiles($CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);

  MODE: {
    if ($Mode =~ /^ContainsNonElementalData$/i) {
      warn "\nWarning: Ignoring compound record number $CmpdCount with ID $CmpdID: Compound contains atom data corresponding to non-elemental atom symbol(s)...\n\n";
      next MODE;
    }

    if ($Mode =~ /^ContainsNoElementalData$/i) {
      warn "\nWarning: Ignoring compound record number $CmpdCount with ID $CmpdID: Compound contains no atom data...\n\n";
      next MODE;
    }

    if ($Mode =~ /^PropertiesCalculationFailed$/i) {
      warn "\nWarning: Ignoring compound record number $CmpdCount with ID $CmpdID: Physicochemical properties calculation didn't succeed...\n\n";
      next MODE;
    }
    warn "\nWarning: Ignoring compound record number $CmpdCount with ID $CmpdID: Physicochemical properties calculation didn't succeed...\n\n";
  }
}

# Check and filter compounds....
#
sub CheckAndFilterCompound {
  my($CmpdCount, $Molecule) = @_;
  my($ElementCount, $NonElementCount);

  ($ElementCount, $NonElementCount) = $Molecule->GetNumOfElementsAndNonElements();

  if ($NonElementCount) {
    ProcessIgnoredCompound('ContainsNonElementalData', $CmpdCount, $Molecule);
    return 1;
  }

  if (!$ElementCount) {
    ProcessIgnoredCompound('ContainsNoElementalData', $CmpdCount, $Molecule);
    return 1;
  }

  return 0;
}

# Write out compounds physicochemical properties calculation summary statistics...
#
sub WriteCalculationSummaryStatistics {
  my($CmpdCount, $IgnoredCmpdCount, $RuleOf5ViolationsCount, $RuleOf3ViolationsCount) = @_;
  my($ProcessedCmpdCount);

  $ProcessedCmpdCount = $CmpdCount - $IgnoredCmpdCount;

  print "\nNumber of compounds: $CmpdCount\n";
  print "Number of compounds processed successfully during physicochemical properties calculation: $ProcessedCmpdCount\n";
  print "Number of compounds ignored during physicochemical properties calculation: $IgnoredCmpdCount\n";

  if ($OptionsInfo{RuleOf5Violations}) {
    print "Number of compounds with one or more RuleOf5 violations: $RuleOf5ViolationsCount\n";
  }

  if ($OptionsInfo{RuleOf3Violations}) {
    print "Number of compounds with one or more RuleOf3 violations: $RuleOf3ViolationsCount\n";
  }

}

# Open output files...
#
sub SetupAndOpenOutputFiles {
  my($FileIndex) = @_;
  my($NewSDFile, $NewTextFile, $NewSDFileRef, $NewTextFileRef);

  $NewSDFileRef = '';
  $NewTextFileRef = '';

  if ($OptionsInfo{SDOutput}) {
    $NewSDFile = $SDFilesInfo{SDOutFileNames}[$FileIndex];
    print "Generating SD file $NewSDFile...\n";
    open NEWSDFILE, ">$NewSDFile" or die "Error: Couldn't open $NewSDFile: $! \n";
    $NewSDFileRef = \*NEWSDFILE;
  }
  if ($OptionsInfo{TextOutput}) {
    $NewTextFile = $SDFilesInfo{TextOutFileNames}[$FileIndex];
    print "Generating text file $NewTextFile...\n";
    open NEWTEXTFILE, ">$NewTextFile" or die "Error: Couldn't open $NewTextFile: $! \n";
    WriteTextFileCoulmnLabels($FileIndex, \*NEWTEXTFILE);
    $NewTextFileRef = \*NEWTEXTFILE;
  }
  return ($NewSDFileRef, $NewTextFileRef);
}

# Write calculated physicochemical properties and other data to appropriate output files...
#
sub WriteDataToOutputFiles {
  my($FileIndex, $CmpdCount, $Molecule, $PhysicochemicalPropertiesDataRef, $NewSDFileRef, $NewTextFileRef) = @_;
  my($PropertyName, $PropertyValue);

  if ($OptionsInfo{SDOutput}) {
    # Retrieve input compound string used to create molecule and write it out
    # without last line containing a delimiter...
    my($CmpdString);
    $CmpdString = $Molecule->GetInputMoleculeString();
    $CmpdString =~ s/\$\$\$\$$//;
    print $NewSDFileRef "$CmpdString";

    # Write out calculated physicochemical properties data...
    for $PropertyName (@{$OptionsInfo{SpecifiedPropertyNames}}) {
      $PropertyValue = $PhysicochemicalPropertiesDataRef->{$PropertyName};
      print $NewSDFileRef  ">  <$PropertyName>\n$PropertyValue\n\n";
    }

    # Write out RuleOf5 violations for molecule....
    if ($OptionsInfo{RuleOf5Violations}) {
      $PropertyValue = $PhysicochemicalPropertiesDataRef->{RuleOf5Violations};
      print $NewSDFileRef  ">  <RuleOf5Violations>\n$PropertyValue\n\n";
    }

    # Write out RuleOf3 violations for molecule....
    if ($OptionsInfo{RuleOf3Violations}) {
      $PropertyValue = $PhysicochemicalPropertiesDataRef->{RuleOf3Violations};
      print $NewSDFileRef  ">  <RuleOf3Violations>\n$PropertyValue\n\n";
    }

    # Write out delimiter...
    print $NewSDFileRef "\$\$\$\$\n";
  }

  if ($OptionsInfo{TextOutput}) {
    my($Line, $DataFieldLabelAndValuesRef, $DataFieldLabel, $DataFieldValue, @LineWords,);

    $DataFieldLabelAndValuesRef = $Molecule->GetDataFieldLabelAndValues();
    @LineWords = ();
    if ($OptionsInfo{DataFieldsMode} =~ /^CompoundID$/i) {
      push @LineWords, SetupCmpdIDForTextFiles($CmpdCount, $Molecule, $DataFieldLabelAndValuesRef);
    }
    elsif ($OptionsInfo{DataFieldsMode} =~ /^All$/i) {
      @LineWords = map { exists $DataFieldLabelAndValuesRef->{$_} ? $DataFieldLabelAndValuesRef->{$_} : ''} @{$SDFilesInfo{AllDataFieldsRef}[$FileIndex]};
    }
    elsif ($OptionsInfo{DataFieldsMode} =~ /^Common$/i) {
      @LineWords = map { exists $DataFieldLabelAndValuesRef->{$_} ? $DataFieldLabelAndValuesRef->{$_} : ''} @{$SDFilesInfo{CommonDataFieldsRef}[$FileIndex]};
    }
    elsif ($OptionsInfo{DataFieldsMode} =~ /^Specify$/i) {
      @LineWords = map { exists $DataFieldLabelAndValuesRef->{$_} ? $DataFieldLabelAndValuesRef->{$_} : ''} @{$OptionsInfo{SpecifiedDataFields}};
    }

    # Append calculated physicochemical properties data...
    for $PropertyName (@{$OptionsInfo{SpecifiedPropertyNames}}) {
      $PropertyValue = $PhysicochemicalPropertiesDataRef->{$PropertyName};
      push @LineWords, $PropertyValue;
    }

    # Write out RuleOf5 violations for molecule....
    if ($OptionsInfo{RuleOf5Violations}) {
      $PropertyValue = $PhysicochemicalPropertiesDataRef->{RuleOf5Violations};
      push @LineWords, $PropertyValue;
    }

    # Write out RuleOf3 violations for molecule....
    if ($OptionsInfo{RuleOf3Violations}) {
      $PropertyValue = $PhysicochemicalPropertiesDataRef->{RuleOf3Violations};
      push @LineWords, $PropertyValue;
    }

    $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
    print $NewTextFileRef "$Line\n";
  }
}

# Write out approriate column labels to text file...
sub WriteTextFileCoulmnLabels {
  my($FileIndex, $NewTextFileRef) = @_;
  my($Line, @LineWords);

  @LineWords = ();
  if ($OptionsInfo{DataFieldsMode} =~ /^All$/i) {
    push @LineWords, @{$SDFilesInfo{AllDataFieldsRef}[$FileIndex]};
  }
  elsif ($OptionsInfo{DataFieldsMode} =~ /^Common$/i) {
    push @LineWords, @{$SDFilesInfo{CommonDataFieldsRef}[$FileIndex]};
  }
  elsif ($OptionsInfo{DataFieldsMode} =~ /^Specify$/i) {
    push @LineWords, @{$OptionsInfo{SpecifiedDataFields}};
  }
  elsif ($OptionsInfo{DataFieldsMode} =~ /^CompoundID$/i) {
    push @LineWords, $OptionsInfo{CompoundIDLabel};
  }
  my($SpecifiedPropertyName);

  # Append physicochemical properties column labels...
  push @LineWords,  @{$OptionsInfo{SpecifiedPropertyNames}};

  # Write out RuleOf5 violations label...
  if ($OptionsInfo{RuleOf5Violations}) {
    push @LineWords, 'RuleOf5Violations';
  }

  # Write out RuleOf3 violations label...
  if ($OptionsInfo{RuleOf3Violations}) {
    push @LineWords, 'RuleOf3Violations';
  }

  $Line = JoinWords(\@LineWords, $OptionsInfo{OutDelim}, $OptionsInfo{OutQuote});
  print $NewTextFileRef "$Line\n";
}

# Generate compound ID for text files..
#
sub SetupCmpdIDForTextFiles {
  my($CmpdCount, $Molecule, $DataFieldLabelAndValuesRef) = @_;
  my($CmpdID);

  $CmpdID = '';
  if ($OptionsInfo{CompoundIDMode} =~ /^MolNameOrLabelPrefix$/i) {
    my($MolName);
    $MolName = $Molecule->GetName();
    $CmpdID = $MolName ? $MolName : "$OptionsInfo{CompoundID}${CmpdCount}";
  }
  elsif ($OptionsInfo{CompoundIDMode} =~ /^LabelPrefix$/i) {
    $CmpdID = "$OptionsInfo{CompoundID}${CmpdCount}";
  }
  elsif ($OptionsInfo{CompoundIDMode} =~ /^DataField$/i) {
    my($SpecifiedDataField);
    $SpecifiedDataField = $OptionsInfo{CompoundID};
    $CmpdID = exists $DataFieldLabelAndValuesRef->{$SpecifiedDataField} ? $DataFieldLabelAndValuesRef->{$SpecifiedDataField} : '';
  }
  elsif ($OptionsInfo{CompoundIDMode} =~ /^MolName$/i) {
    $CmpdID = $Molecule->GetName();
  }
  return $CmpdID;
}

# Calculate physicochemical properties for molecule...
#
sub CalculateMoleculeProperties {
  my($MolecularDescriptorsGenerator, $Molecule) = @_;
  my($PropertyName, $PropertyValue, $MolecularDescriptorsObject, %CalculatedPhysicochemicalProperties);

  %CalculatedPhysicochemicalProperties = ();

  if ($OptionsInfo{KeepLargestComponent}) {
    $Molecule->KeepLargestComponent();
  }

  if (!$Molecule->DetectRings()) {
    return undef;
  }
  $Molecule->SetAromaticityModel($OptionsInfo{AromaticityModel});
  $Molecule->DetectAromaticity();

  if ($OptionsInfo{AddHydrogens}) {
    $Molecule->AddHydrogens();
  }

  # Calculate physicochemical properties...
  $MolecularDescriptorsGenerator->SetMolecule($Molecule);
  $MolecularDescriptorsGenerator->GenerateDescriptors();

  if (!$MolecularDescriptorsGenerator->IsDescriptorsGenerationSuccessful()) {
    return undef;
  }

  %CalculatedPhysicochemicalProperties = $MolecularDescriptorsGenerator->GetDescriptorNamesAndValues();

  # Count RuleOf3 violations...
  if ($OptionsInfo{RuleOf3Violations}) {
    CalculateRuleViolationsCount('RuleOf3Violations', \%CalculatedPhysicochemicalProperties);
  }

  # Count RuleOf5 violations...
  if ($OptionsInfo{RuleOf5Violations}) {
    CalculateRuleViolationsCount('RuleOf5Violations', \%CalculatedPhysicochemicalProperties);
  }

  return \%CalculatedPhysicochemicalProperties;
}

# Setup molecular descriptor generator to calculate property values for specifed
# property names...
#
sub SetupMolecularDescriptorsGenerator {
  my($PropertyName, $MolecularDescriptorsGenerator);

  $MolecularDescriptorsGenerator = new MolecularDescriptors::MolecularDescriptorsGenerator('Mode' => 'Specify', 'DescriptorNames' => \@{$OptionsInfo{SpecifiedPropertyNames}});

  # Setup molecular desciptor calculation parameters...
  if (exists($OptionsInfo{SpecifiedPropertyNamesMap}{lc('MolecularWeight')}) || exists($OptionsInfo{SpecifiedPropertyNamesMap}{lc('ExactMass')}) ) {
    $MolecularDescriptorsGenerator->SetDescriptorClassParameters('DescriptorClassName' => 'WeightAndMassDescriptors', %{$OptionsInfo{PrecisionParametersMap}});
  }

  if (exists($OptionsInfo{SpecifiedPropertyNamesMap}{lc('RotatableBonds')})) {
    $MolecularDescriptorsGenerator->SetDescriptorClassParameters('DescriptorClassName' => 'RotatableBondsDescriptors', %{$OptionsInfo{RotatableBondsParametersMap}});
  }

  if (exists($OptionsInfo{SpecifiedPropertyNamesMap}{lc('HydrogenBondDonors')}) || exists($OptionsInfo{SpecifiedPropertyNamesMap}{lc('HydrogenBondAcceptors')}) ) {
    $MolecularDescriptorsGenerator->SetDescriptorClassParameters('DescriptorClassName' => 'HydrogenBondsDescriptors', 'HydrogenBondsType' => $OptionsInfo{HydrogenBonds});
  }

  if (exists($OptionsInfo{SpecifiedPropertyNamesMap}{lc('TPSA')})) {
    $MolecularDescriptorsGenerator->SetDescriptorClassParameters('DescriptorClassName' => 'TPSADescriptors', %{$OptionsInfo{TPSAParametersMap}});
  }

  if (exists($OptionsInfo{SpecifiedPropertyNamesMap}{lc('MolecularComplexity')})) {
    $MolecularDescriptorsGenerator->SetDescriptorClassParameters('DescriptorClassName' => 'MolecularComplexityDescriptors', %{$OptionsInfo{MolecularComplexityParametersMap}});
  }

  return $MolecularDescriptorsGenerator;
}

# Calculate RuleOf3 or RuleOf5 violations count...
#
sub CalculateRuleViolationsCount {
  my($RuleViolationsType, $CalculatedPropertiesMapRef) = @_;
  my($RuleViolationsCount, $PropertyName);

  $RuleViolationsCount = 0;

  RULEVIOLATIONSTYPE: {
    if ($RuleViolationsType =~ /^RuleOf3Violations$/i) {
      for $PropertyName (@{$OptionsInfo{RuleOf3PropertyNames}}) {
	if ($CalculatedPropertiesMapRef->{$PropertyName} > $OptionsInfo{RuleOf3MaxPropertyValuesMap}{$PropertyName}) {
	  $RuleViolationsCount++;
	}
      }
      last RULEVIOLATIONSTYPE;
    }

    if ($RuleViolationsType =~ /^RuleOf5Violations$/i) {
      for $PropertyName (@{$OptionsInfo{RuleOf5PropertyNames}}) {
	if ($CalculatedPropertiesMapRef->{$PropertyName} > $OptionsInfo{RuleOf5MaxPropertyValuesMap}{$PropertyName}) {
	  $RuleViolationsCount++;
	}
      }
      last RULEVIOLATIONSTYPE;
    }

    die "Warning: Unknown rule violation type: $RuleViolationsType...";
  }

  # Set rule violation count...
  $CalculatedPropertiesMapRef->{$RuleViolationsType} = $RuleViolationsCount;

}

# Retrieve information about SD files...
#
sub RetrieveSDFilesInfo {
  my($SDFile, $Index, $FileDir, $FileExt, $FileName, $OutFileRoot, $TextOutFileExt, $SDOutFileExt, $NewSDFileName, $NewTextFileName, $CheckDataField, $CollectDataFields, $AllDataFieldsRef, $CommonDataFieldsRef);

  %SDFilesInfo = ();
  @{$SDFilesInfo{FileOkay}} = ();
  @{$SDFilesInfo{OutFileRoot}} = ();
  @{$SDFilesInfo{SDOutFileNames}} = ();
  @{$SDFilesInfo{TextOutFileNames}} = ();
  @{$SDFilesInfo{AllDataFieldsRef}} = ();
  @{$SDFilesInfo{CommonDataFieldsRef}} = ();

  $CheckDataField = ($OptionsInfo{TextOutput} && ($OptionsInfo{DataFieldsMode} =~ /^CompoundID$/i) && ($OptionsInfo{CompoundIDMode} =~ /^DataField$/i)) ? 1 : 0;
  $CollectDataFields = ($OptionsInfo{TextOutput} && ($OptionsInfo{DataFieldsMode} =~ /^(All|Common)$/i)) ? 1 : 0;

  FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFile = $SDFilesList[$Index];

    $SDFilesInfo{FileOkay}[$Index] = 0;
    $SDFilesInfo{OutFileRoot}[$Index] = '';
    $SDFilesInfo{SDOutFileNames}[$Index] = '';
    $SDFilesInfo{TextOutFileNames}[$Index] = '';

    $SDFile = $SDFilesList[$Index];
    if (!(-e $SDFile)) {
      warn "Warning: Ignoring file $SDFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($SDFile, "sd sdf")) {
      warn "Warning: Ignoring file $SDFile: It's not a SD file\n";
      next FILELIST;
    }

    if ($CheckDataField) {
      # Make sure data field exists in SD file..
      my($CmpdString, $SpecifiedDataField, @CmpdLines, %DataFieldValues);

      @CmpdLines = ();
      open SDFILE, "$SDFile" or die "Error: Couldn't open $SDFile: $! \n";
      $CmpdString = ReadCmpdString(\*SDFILE);
      close SDFILE;
      @CmpdLines = split "\n", $CmpdString;
      %DataFieldValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
      $SpecifiedDataField = $OptionsInfo{CompoundID};
      if (!exists $DataFieldValues{$SpecifiedDataField}) {
	warn "Warning: Ignoring file $SDFile: Data field value, $SpecifiedDataField, using  \"--CompoundID\" option in \"DataField\" \"--CompoundIDMode\" doesn't exist\n";
	next FILELIST;
      }
    }

    $AllDataFieldsRef = '';
    $CommonDataFieldsRef = '';
    if ($CollectDataFields) {
      my($CmpdCount);
      open SDFILE, "$SDFile" or die "Error: Couldn't open $SDFile: $! \n";
      ($CmpdCount, $AllDataFieldsRef, $CommonDataFieldsRef) = GetAllAndCommonCmpdDataHeaderLabels(\*SDFILE);
      close SDFILE;
    }

    # Setup output file names...
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($SDFile);

    $TextOutFileExt = "csv";
    if ($Options{outdelim} =~ /^tab$/i) {
      $TextOutFileExt = "tsv";
    }
    $SDOutFileExt = $FileExt;

    if ($OptionsInfo{OutFileRoot} && (@SDFilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($OptionsInfo{OutFileRoot});
      if ($RootFileName && $RootFileExt) {
	$FileName = $RootFileName;
      }
      else {
	$FileName = $OptionsInfo{OutFileRoot};
      }
      $OutFileRoot = $FileName;
    }
    else {
      $OutFileRoot = "${FileName}PhysicochemicalProperties";
    }

    $NewSDFileName = "${OutFileRoot}.${SDOutFileExt}";
    $NewTextFileName = "${OutFileRoot}.${TextOutFileExt}";

    if ($OptionsInfo{SDOutput}) {
      if ($SDFile =~ /$NewSDFileName/i) {
	warn "Warning: Ignoring input file $SDFile: Same output, $NewSDFileName, and input file names.\n";
	print "Specify a different name using \"-r --root\" option or use default name.\n";
	next FILELIST;
      }
    }

    if (!$OptionsInfo{OverwriteFiles}) {
      # Check SD and text outout files...
      if ($OptionsInfo{SDOutput}) {
	if (-e $NewSDFileName) {
	  warn "Warning: Ignoring file $SDFile: The file $NewSDFileName already exists\n";
	  next FILELIST;
	}
      }
      if ($OptionsInfo{TextOutput}) {
	if (-e $NewTextFileName) {
	  warn "Warning: Ignoring file $SDFile: The file $NewTextFileName already exists\n";
	  next FILELIST;
	}
      }
    }

    $SDFilesInfo{FileOkay}[$Index] = 1;

    $SDFilesInfo{OutFileRoot}[$Index] = $OutFileRoot;
    $SDFilesInfo{SDOutFileNames}[$Index] = $NewSDFileName;
    $SDFilesInfo{TextOutFileNames}[$Index] = $NewTextFileName;

    $SDFilesInfo{AllDataFieldsRef}[$Index] = $AllDataFieldsRef;
    $SDFilesInfo{CommonDataFieldsRef}[$Index] = $CommonDataFieldsRef;
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{AromaticityModel} = $Options{aromaticitymodel};

  # Process property name related options...
  ProcessPropertyNamesOption();

  # Setup RuleOf3 and RuleOf5 violation calculations...
  $OptionsInfo{RuleOf3Violations} = ($Options{ruleof3violations} =~ /^Yes$/i) ? 1 : 0;
  $OptionsInfo{RuleOf5Violations} = ($Options{ruleof5violations} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{CompoundIDMode} = $Options{compoundidmode};
  $OptionsInfo{CompoundIDLabel} = $Options{compoundidlabel};
  $OptionsInfo{DataFieldsMode} = $Options{datafieldsmode};

  my(@SpecifiedDataFields);
  @SpecifiedDataFields = ();

  @{$OptionsInfo{SpecifiedDataFields}} = ();
  $OptionsInfo{CompoundID} = '';

  if ($Options{datafieldsmode} =~ /^CompoundID$/i) {
    if ($Options{compoundidmode} =~ /^DataField$/i) {
      if (!$Options{compoundid}) {
	die "Error: You must specify a value for \"--CompoundID\" option in \"DataField\" \"--CompoundIDMode\". \n";
      }
      $OptionsInfo{CompoundID} = $Options{compoundid};
    }
    elsif ($Options{compoundidmode} =~ /^(LabelPrefix|MolNameOrLabelPrefix)$/i) {
      $OptionsInfo{CompoundID} = $Options{compoundid} ? $Options{compoundid} : 'Cmpd';
    }
  }
  elsif ($Options{datafieldsmode} =~ /^Specify$/i) {
    if (!$Options{datafields}) {
      die "Error: You must specify a value for \"--DataFields\" option in \"Specify\" \"-d, --DataFieldsMode\". \n";
    }
    @SpecifiedDataFields = split /\,/, $Options{datafields};
    push @{$OptionsInfo{SpecifiedDataFields}}, @SpecifiedDataFields;
  }

  # Types of hydrogen bonds...
  $OptionsInfo{HydrogenBonds} = $Options{hydrogenbonds};

  # Process precision value parameters...
  ProcessPrecisionOption();

  # Process rotatable bonds parameters...
  ProcessRotatableBondsOption();

  # Process TPSA parameters...
  ProcessTPSAOption();

  # Process molecular complexity parameters...
  ProcessMolecularComplexityOption();

  $OptionsInfo{Filter} = ($Options{filter} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{KeepLargestComponent} = ($Options{keeplargestcomponent} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{Output} = $Options{output};
  $OptionsInfo{SDOutput} = ($Options{output} =~ /^(SD|Both)$/i) ? 1 : 0;
  $OptionsInfo{TextOutput} = ($Options{output} =~ /^(Text|Both)$/i) ? 1 : 0;

  $OptionsInfo{OutDelim} = ($Options{outdelim} =~ /tab/i ) ? "\t" : (($Options{outdelim} =~ /semicolon/i) ? "\;" : "\,");
  $OptionsInfo{OutQuote} = ($Options{quote} =~ /^Yes$/i) ? 1 : 0;

  $OptionsInfo{OverwriteFiles} = $Options{overwrite} ? 1 : 0;
  $OptionsInfo{OutFileRoot} = $Options{root} ? $Options{root} : 0;
}

# Process property name related options...
#
sub ProcessPropertyNamesOption {

  # Setup supported physicochemical properties...
  my($SupportedProperty);

  @{$OptionsInfo{SupportedPropertyNames}} = ();
  %{$OptionsInfo{SupportedPropertyNamesMap}} = ();

  @{$OptionsInfo{RuleOf5PropertyNames}} = ();
  %{$OptionsInfo{RuleOf5MaxPropertyValuesMap}} = ();

  @{$OptionsInfo{RuleOf3PropertyNames}} = ();
  %{$OptionsInfo{RuleOf3MaxPropertyValuesMap}} = ();

  @{$OptionsInfo{DefaultPropertyNames}} = ();

  @{$OptionsInfo{SupportedPropertyNames}} = qw(MolecularWeight ExactMass HeavyAtoms Rings AromaticRings MolecularVolume RotatableBonds HydrogenBondDonors HydrogenBondAcceptors SLogP SMR TPSA Fsp3Carbons Sp3Carbons MolecularComplexity);

  @{$OptionsInfo{RuleOf5PropertyNames}} = qw(MolecularWeight HydrogenBondDonors HydrogenBondAcceptors SLogP);
  %{$OptionsInfo{RuleOf5MaxPropertyValuesMap}} = ('MolecularWeight' => 500, 'HydrogenBondDonors' => 5, 'HydrogenBondAcceptors' => 10,  'SLogP' => 5);

  @{$OptionsInfo{RuleOf3PropertyNames}} = qw(MolecularWeight RotatableBonds HydrogenBondDonors HydrogenBondAcceptors SLogP TPSA);
  %{$OptionsInfo{RuleOf3MaxPropertyValuesMap}} = ('MolecularWeight' => 300, 'RotatableBonds' => 3, 'HydrogenBondDonors' => 3, 'HydrogenBondAcceptors' => 3, 'SLogP' => 3, 'TPSA' => 60);

  @{$OptionsInfo{DefaultPropertyNames}} = qw(MolecularWeight HeavyAtoms MolecularVolume RotatableBonds HydrogenBondDonors HydrogenBondAcceptors SLogP TPSA);

  for $SupportedProperty (@{$OptionsInfo{SupportedPropertyNames}}) {
    $OptionsInfo{SupportedPropertyNamesMap}{lc($SupportedProperty)} = $SupportedProperty;
  }

  # Process specified properties....
  my($SpecifiedPropertyName, @SpecifiedPropertyNames, %SpecifiedPropertyNamesMap);

  @SpecifiedPropertyNames = ();
  %SpecifiedPropertyNamesMap = ();

  @{$OptionsInfo{SpecifiedPropertyNames}} = ();
  %{$OptionsInfo{SpecifiedPropertyNamesMap}} = ();

  if ($Options{mode} =~ /^All$/i) {
    @SpecifiedPropertyNames = @{$OptionsInfo{SupportedPropertyNames}};
  }
  elsif ($Options{mode} =~ /^RuleOf5$/i) {
    @SpecifiedPropertyNames = @{$OptionsInfo{RuleOf5PropertyNames}};
  }
  elsif ($Options{mode} =~ /^RuleOf3$/i) {
    @SpecifiedPropertyNames = @{$OptionsInfo{RuleOf3PropertyNames}};
  }
  elsif (IsEmpty($Options{mode})) {
    @SpecifiedPropertyNames = @{$OptionsInfo{DefaultPropertyNames}};
  }
  else {
    # Comma delimited lisr of specified property names...
    my($Mode, $PropertyName, @PropertyNames, @UnsupportedPropertyNames);

    $Mode = $Options{mode};
    $Mode =~ s/ //g;

    @PropertyNames = split ",", $Mode;
    @UnsupportedPropertyNames = ();

    for $PropertyName (@PropertyNames) {
      if (exists($OptionsInfo{SupportedPropertyNamesMap}{lc($PropertyName)})) {
	push @SpecifiedPropertyNames, $PropertyName;
      }
      else {
	push @UnsupportedPropertyNames, $PropertyName;
      }
    }
    if (@UnsupportedPropertyNames) {
      if (@UnsupportedPropertyNames > 1) {
	warn "Error: The physicochemical property names specified - ", JoinWords(\@UnsupportedPropertyNames, ", ", 0)," - for option \"-m --mode\" are not valid.\n";
      }
      else {
	warn "Error: The physicochemical property name specified, @UnsupportedPropertyNames , for option \"-m --mode\" is not valid.\n";
      }
      die "Allowed values:", JoinWords(\@{$OptionsInfo{SupportedPropertyNames}}, ", ", 0), "\n";
    }
    if (!@SpecifiedPropertyNames) {
      die "Error: No valid physicochemical property names specified for option \"-m --mode\".\n";
    }
  }

  # Set up specified property names map...
  PROPERTY: for $SpecifiedPropertyName (@SpecifiedPropertyNames) {
    if (exists $SpecifiedPropertyNamesMap{lc($SpecifiedPropertyName)}) {
      warn "Warning: The physicochemical property name, $SpecifiedPropertyName, is specified multiple times as value of option \"-m --mode\" .\n";
      next PROPERTY;
    }
    # Canonical specified property name...
    $SpecifiedPropertyNamesMap{lc($SpecifiedPropertyName)} = $OptionsInfo{SupportedPropertyNamesMap}{lc($SpecifiedPropertyName)};
  }

  # Make sure for calculation of  RuleOf3Violations, all appropriate property names are specified...
  if ($Options{ruleof3violations} =~ /^Yes$/i && $Options{mode} =~ /^RuleOf5$/i) {
    die "Error: The value specified, $Options{ruleof3violations}, for  \"--RuleOf3Violations\" option in \"RuleOf5\" \"-m --Mode\" is not valid. You must specify RuleOf3 value for \"-m --Mode\" to calculate RuleOf3 violations.\n";
  }

  if ($Options{ruleof3violations} =~ /^Yes$/i) {
    my($RuleOf3PropertyName, @MissingRuleOf3Names);

    @MissingRuleOf3Names = ();
    PROPERTY: for $RuleOf3PropertyName (@{$OptionsInfo{RuleOf3PropertyNames}}) {
      if (exists $SpecifiedPropertyNamesMap{lc($RuleOf3PropertyName)}) {
	next PROPERTY;
      }
      push @MissingRuleOf3Names, $RuleOf3PropertyName;

      # Add property name to specified properties names list and map...
      push @SpecifiedPropertyNames, $RuleOf3PropertyName;
      $SpecifiedPropertyNamesMap{lc($RuleOf3PropertyName)} = $OptionsInfo{SupportedPropertyNamesMap}{lc($RuleOf3PropertyName)};
    }
    if (@MissingRuleOf3Names) {
      warn "Warning: The following physicochemical property names not specified in \"-m --Mode\" option are required for calculating RuleOf3Violations and have been added to the list of property names: @MissingRuleOf3Names\n";
    }
  }

  # Make sure for calculation of  RuleOf5Violations, all appropriate property names are specified...
  if ($Options{ruleof5violations} =~ /^Yes$/i && $Options{mode} =~ /^RuleOf3$/i) {
    die "Error: The value specified, $Options{ruleof5violations}, for  \"--RuleOf5Violations\" option in \"RuleOf3\" \"-m --Mode\" is not valid. You must specify RuleOf5 value for \"-m --Mode\" to calculate RuleOf5 violations.\n";
  }

  if ($Options{ruleof5violations} =~ /^Yes$/i) {
    my($RuleOf5PropertyName, @MissingRuleOf5Names);

    @MissingRuleOf5Names = ();
    PROPERTY: for $RuleOf5PropertyName (@{$OptionsInfo{RuleOf5PropertyNames}}) {
      if (exists $SpecifiedPropertyNamesMap{lc($RuleOf5PropertyName)}) {
	next PROPERTY;
      }
      push @MissingRuleOf5Names, $RuleOf5PropertyName;

      # Add property name to specified properties names list and map...
      push @SpecifiedPropertyNames, $RuleOf5PropertyName;
      $SpecifiedPropertyNamesMap{lc($RuleOf5PropertyName)} = $OptionsInfo{SupportedPropertyNamesMap}{lc($RuleOf5PropertyName)};
    }
    if (@MissingRuleOf5Names) {
      warn "Warning: The following physicochemical property names not specified in \"-m --Mode\" option are required for calculating RuleOf5Violations and have been added to the list of property names: @MissingRuleOf5Names\n";
    }
  }
  $OptionsInfo{Mode} = $Options{mode};

  # Setup canonical specified property names corresponding to supported names in mixed case...
  my(@SpecifiedCanonicalPropertyNames);

  @SpecifiedCanonicalPropertyNames = ();
  for $SpecifiedPropertyName (@SpecifiedPropertyNames) {
    push @SpecifiedCanonicalPropertyNames, $SpecifiedPropertyNamesMap{lc($SpecifiedPropertyName)};
  }
  @{$OptionsInfo{SpecifiedPropertyNames}} = @SpecifiedCanonicalPropertyNames;
  %{$OptionsInfo{SpecifiedPropertyNamesMap}} = %SpecifiedPropertyNamesMap;

  # Based on specified property names, figure out whether hydrogens need to be added before
  # calculation of properties...
  #
  $OptionsInfo{AddHydrogens} = 0;
  if (exists($SpecifiedPropertyNamesMap{lc('MolecularVolume')}) || exists($SpecifiedPropertyNamesMap{lc('SLogP')}) || exists($SpecifiedPropertyNamesMap{lc('SMR')})) {
    $OptionsInfo{AddHydrogens} = 1;
  }
}

# Process precision option...
#
sub ProcessPrecisionOption {
  my($ParameterName, $ParameterValue, %PrecisionParametersMap, %PrecisionParameterNamesMap);

  %{$OptionsInfo{PrecisionParametersMap}} = ();

  %PrecisionParametersMap = ('WeightPrecision' => 2, 'MassPrecision' => 4);
  %PrecisionParameterNamesMap = ('molecularweight' => 'WeightPrecision', 'exactmass' => 'MassPrecision');

  if ($Options{precision}) {
    # Process specified values...
    my($Index, $SpecifiedPrecision, @SpecifiedPrecisionValuePairs);

    $SpecifiedPrecision = $Options{precision};
    $SpecifiedPrecision =~ s/ //g;
    @SpecifiedPrecisionValuePairs = split ",", $SpecifiedPrecision;
    if (@SpecifiedPrecisionValuePairs % 2) {
      die "Error: Invalid number of values specified using \"--Precision\" option: It must contain even number of values.\n";
    }
    for ($Index = 0; (($Index + 1) < @SpecifiedPrecisionValuePairs); $Index += 2 ) {
      $ParameterName = $SpecifiedPrecisionValuePairs[$Index];
      $ParameterValue = $SpecifiedPrecisionValuePairs[$Index + 1];
      if (!exists $PrecisionParameterNamesMap{lc($ParameterName)}) {
	die "Error: The precision parameter name specified, $ParameterName, for option \"--Precision\" is not valid.\n";
      }
      if (!IsPositiveInteger($ParameterValue)) {
	die "Error: The parameter value specified, $ParameterValue, for parameter name, $ParameterName in option \"--Precision\" is not valid. Allowed values: positive integer. \n";
      }
      $ParameterName = $PrecisionParameterNamesMap{lc($ParameterName)};
      $PrecisionParametersMap{$ParameterName} = $ParameterValue;
    }
  }
  $OptionsInfo{Precision} = $Options{precision};
  %{$OptionsInfo{PrecisionParametersMap}} = %PrecisionParametersMap;
}

# Process rotatable bonds option...
sub ProcessRotatableBondsOption {
  my($ParameterName, $ParameterValue, %RotatableBondsParametersMap, %RotatableBondsParameterNamesMap);

  %{$OptionsInfo{RotatableBondsParametersMap}} = ();
  %RotatableBondsParametersMap = ('IgnoreTerminalBonds' => 1, 'IgnoreBondsToTripleBonds' => 1, 'IgnoreAmideBonds' => 1, 'IgnoreThioamideBonds' => 1, 'IgnoreSulfonamideBonds' => 1);

  for $ParameterName (keys %RotatableBondsParametersMap) {
    $RotatableBondsParameterNamesMap{lc($ParameterName)} = $ParameterName;
  }

  if ($Options{rotatablebonds}) {
    # Process specified values...
    my($Index, $SpecifiedRotatableBonds, @SpecifiedRotatableBondsValuePairs);

    $SpecifiedRotatableBonds = $Options{rotatablebonds};
    $SpecifiedRotatableBonds =~ s/ //g;
    @SpecifiedRotatableBondsValuePairs = split ",", $SpecifiedRotatableBonds;
    if (@SpecifiedRotatableBondsValuePairs % 2) {
      die "Error: Invalid number of values specified using \"--RotatableBonds\" option: It must contain even number of values.\n";
    }
    for ($Index = 0; (($Index + 1) < @SpecifiedRotatableBondsValuePairs); $Index += 2 ) {
      $ParameterName = $SpecifiedRotatableBondsValuePairs[$Index];
      $ParameterValue = $SpecifiedRotatableBondsValuePairs[$Index + 1];
      if (!exists $RotatableBondsParameterNamesMap{lc($ParameterName)}) {
	die "Error: The rotatable bonds parameter name specified, $ParameterName, for option \"--RotatableBonds\" is not valid.\n";
      }
      if ($ParameterValue !~ /^(Yes|No)$/i) {
	die "Error: The parameter value specified, $ParameterValue, for parameter name, $ParameterName in option \"--RotatableBonds\" is not valid. Allowed values: Yes or No. \n";
      }
      $ParameterName = $RotatableBondsParameterNamesMap{lc($ParameterName)};
      $ParameterValue = ($ParameterValue =~ /^Yes$/i) ? 1 : 0;
      $RotatableBondsParametersMap{$ParameterName} = $ParameterValue;
    }
  }
  $OptionsInfo{RotatableBonds} = $Options{rotatablebonds};
  %{$OptionsInfo{RotatableBondsParametersMap}} = %RotatableBondsParametersMap;
}

# Process TPSA option...
#
sub ProcessTPSAOption {
  my($ParameterName, $ParameterValue, %TPSAParametersMap, %TPSAParameterNamesMap);

  %{$OptionsInfo{TPSAParametersMap}} = ();

  %TPSAParametersMap = ('IgnorePhosphorus' => 1, 'IgnoreSulfur' => 1);
  for $ParameterName (keys %TPSAParametersMap) {
    $TPSAParameterNamesMap{lc($ParameterName)} = $ParameterName;
  }

  if ($Options{tpsa}) {
    # Process specified values...
    my($Index, $SpecifiedTPSA, @SpecifiedTPSAValuePairs);

    $SpecifiedTPSA = $Options{tpsa};
    $SpecifiedTPSA =~ s/ //g;
    @SpecifiedTPSAValuePairs = split ",", $SpecifiedTPSA;
    if (@SpecifiedTPSAValuePairs % 2) {
      die "Error: Invalid number of values specified using \"--TPSA\" option: It must contain even number of values.\n";
    }
    for ($Index = 0; (($Index + 1) < @SpecifiedTPSAValuePairs); $Index += 2 ) {
      $ParameterName = $SpecifiedTPSAValuePairs[$Index];
      $ParameterValue = $SpecifiedTPSAValuePairs[$Index + 1];
      if (!exists $TPSAParameterNamesMap{lc($ParameterName)}) {
	die "Error: The TPSA parameter name specified, $ParameterName, for option \"--TPSA\" is not valid.\n";
      }
      if ($ParameterValue !~ /^(Yes|No)$/i) {
	die "Error: The parameter value specified, $ParameterValue, for parameter name, $ParameterName in option \"--TPSA\" is not valid. Allowed values: Yes or No. \n";
      }
      $ParameterName = $TPSAParameterNamesMap{lc($ParameterName)};
      $ParameterValue = ($ParameterValue =~ /^Yes$/i) ? 1 : 0;
      $TPSAParametersMap{$ParameterName} = $ParameterValue;
    }
  }
  $OptionsInfo{TPSA} = $Options{tpsa};
  %{$OptionsInfo{TPSAParametersMap}} = %TPSAParametersMap;
}

# Process molecular complexity parameters...
#
sub ProcessMolecularComplexityOption {
  my($MolecularComplexityType, $ParameterName, $ParameterValue, @ParameterNames, @ParameterValues, @AtomIdentifierTypeParameters, %ComplexityParametersMap, %ComplexityParameterNamesMap);

  %{$OptionsInfo{MolecularComplexityParametersMap}} = ();

  %ComplexityParametersMap = ('MolecularComplexityType' => '', 'AtomIdentifierType' => '',
			      'AtomicInvariantsToUse' => '', 'FunctionalClassesToUse' => '',
			      'MACCSKeysSize' => '166', 'NeighborhoodRadius' => '2',
			      'MinPathLength' => '1', 'MaxPathLength' => '8', 'UseBondSymbols' => '1',
			      'MinDistance' => '1', 'MaxDistance' => '10', 'UseTriangleInequality' => '',
			      'DistanceBinSize' => '2', 'NormalizationMethodology' => 'None');

  %ComplexityParameterNamesMap = ();
  for $ParameterName (keys %ComplexityParametersMap) {
    $ComplexityParameterNamesMap{lc($ParameterName)} = $ParameterName;
  }

  if ($Options{molecularcomplexity}) {
    # Process specified values...
    my($Index, $SpecifiedComplexity, @SpecifiedComplexityValuePairs);

    $SpecifiedComplexity = $Options{molecularcomplexity};

    @SpecifiedComplexityValuePairs = split ",", $SpecifiedComplexity;
    if (@SpecifiedComplexityValuePairs % 2) {
      die "Error: Invalid number of values specified using \"--MolecularComplexity\" option: It must contain even number of values.\n";
    }

    for ($Index = 0; (($Index + 1) < @SpecifiedComplexityValuePairs); $Index += 2 ) {
      $ParameterName = $SpecifiedComplexityValuePairs[$Index];
      $ParameterValue = $SpecifiedComplexityValuePairs[$Index + 1];

      $ParameterName = RemoveLeadingAndTrailingWhiteSpaces($ParameterName);
      $ParameterValue = RemoveLeadingAndTrailingWhiteSpaces($ParameterValue);

      if (!exists $ComplexityParameterNamesMap{lc($ParameterName)}) {
	die "Error: The molecular complexity parameter name specified, $ParameterName, for option \"--MolecularComplexity\" is not valid.\n";
      }
      $ParameterName = $ComplexityParameterNamesMap{lc($ParameterName)};

      if ($ParameterName =~ /^AtomicInvariantsToUse$/i) {
	my($AtomSymbolFound);

	$AtomSymbolFound = 0;
	@ParameterValues = split(' ', $ParameterValue);
	for $ParameterValue (@ParameterValues) {
	  if (!AtomTypes::AtomicInvariantsAtomTypes::IsAtomicInvariantAvailable($ParameterValue)) {
	    die "Error: The atomic invariant specified, $ParameterValue, for  AtomicInvariantsToUse in option \"--MolecularComplexity\" is not valid.\n";
	  }
	  if ($ParameterValue =~ /^(AS|AtomSymbol)$/i) {
	    $AtomSymbolFound = 1;
	  }
	}
	if (!$AtomSymbolFound) {
	  die "Error: The atomic invariants specified using AtomicInvariantsToUse in option \"--MolecularComplexity\" is not valid: AtomicInvariant atom symbol, AS or AtomSymbol, must be specified.\n";
	}
	$ParameterValue = JoinWords(\@ParameterValues, ",", 0);
      }
      elsif ($ParameterName =~ /^FunctionalClassesToUse$/i) {
	@ParameterValues = split(' ', $ParameterValue);
	for $ParameterValue (@ParameterValues) {
	  if (!AtomTypes::FunctionalClassAtomTypes::IsFunctionalClassAvailable($ParameterValue)) {
	    die "Error: The functional class specified, $ParameterValue, for  FunctionalClassesToUse in option \"--MolecularComplexity\" is not valid.\n";
	  }
	}
	$ParameterValue = JoinWords(\@ParameterValues, ",", 0);
      }
      else {
	if ($ParameterValue =~ / /) {
	  $ParameterValue =~ s/ //g;
	}
	if ($ParameterValue =~ /^(Yes|No)$/i) {
	  $ParameterValue = ($ParameterValue =~ /^Yes$/i) ? 1 : 0;
	}
      }

      if ($ParameterName =~ /^MolecularComplexityType$/i) {
	if ($ParameterValue !~ /^(AtomTypesFingerprints|ExtendedConnectivityFingerprints|MACCSKeys|PathLengthFingerprints|TopologicalAtomPairsFingerprints|TopologicalAtomTripletsFingerprints|TopologicalAtomTorsionsFingerprints|TopologicalPharmacophoreAtomPairsFingerprints|TopologicalPharmacophoreAtomTripletsFingerprints)$/i) {
	  die "Error: The parameter value specified, $ParameterValue, for parameter name, $ParameterName in option \"--MolecularComplexity\" is not valid. Allowed values: AtomTypesFingerprints, ExtendedConnectivityFingerprints, MACCSKeys, PathLengthFingerprints, TopologicalAtomPairsFingerprints, TopologicalAtomTripletsFingerprints, TopologicalAtomTorsionsFingerprints, TopologicalPharmacophoreAtomPairsFingerprints, or TopologicalPharmacophoreAtomTripletsFingerprints..\n";
	}
      }
      elsif ($ParameterName =~ /^AtomIdentifierType$/i) {
	if ($ParameterValue !~ /^(AtomicInvariantsAtomTypes|FunctionalClassAtomTypes|DREIDINGAtomTypes|EStateAtomTypes|MMFF94AtomTypes|SLogPAtomTypes|SYBYLAtomTypes|TPSAAtomTypes|UFFAtomTypes)$/i) {
	  die "Error: The parameter value specified, $ParameterValue, for parameter name, $ParameterName in option \"--MolecularComplexity\" is not valid. Supported atom identifier types in current release of MayaChemTools: AtomicInvariantsAtomTypes, FunctionalClassAtomTypes, DREIDINGAtomTypes, EStateAtomTypes, MMFF94AtomTypes, SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes and UFFAtomTypes.\n";
	}
      }
      elsif ($ParameterName =~ /^(MACCSKeysSize|MinPathLength|MaxPathLength|MinDistance|MaxDistance|DistanceBinSize)$/i) {
	if (!IsPositiveInteger($ParameterValue)) {
	  die "Error: The parameter value specified, $ParameterValue, for parameter name, $ParameterName in option \"--MolecularComplexity\" is not valid. Allowed values: positive integer. \n";
	}
      }
      elsif ($ParameterName =~ /^NeighborhoodRadius$/i) {
	if (!(IsInteger($ParameterValue) && $ParameterValue >=0)) {
	  die "Error: The parameter value specified, $ParameterValue, for parameter name, $ParameterName in option \"--MolecularComplexity\" is not valid. Allowed values: 0 or positive integer. \n";
	}
      }
      elsif ($ParameterName =~ /^NormalizationMethodology$/i) {
	if ($ParameterValue !~ /^(None|ByHeavyAtomsCount|ByPossibleKeysCount)$/i) {
	  die "Error: The parameter value specified, $ParameterValue, for parameter name, $ParameterName in option \"--MolecularComplexity\" is not valid. Allowed values: None, ByHeavyAtomsCount, or ByPossibleKeysCount\n";
	}
      }
      $ComplexityParametersMap{$ParameterName} = $ParameterValue;
    }

    if ($ComplexityParametersMap{MACCSKeysSize} !~ /^(166|322)$/i) {
      die "Error: The parameter value specified, $ComplexityParametersMap{MACCSKeysSize}, for parameter name, MACCSKeysSize in option \"--MolecularComplexity\" is not valid. Allowed values: 166 or 322\n";
    }
    if ($ComplexityParametersMap{MinPathLength} > $ComplexityParametersMap{MaxPathLength}) {
      die "Error: The parameter value specified for MinPathLength, $ComplexityParametersMap{MinPathLength}, must be <= MaxPathLength, $ComplexityParametersMap{MaxPathLength} ...\n";
    }
    if ($ComplexityParametersMap{MinDistance} > $ComplexityParametersMap{MaxDistance}) {
      die "Error: The parameter value specified for MinDistance, $ComplexityParametersMap{MinDistance}, must be <= MaxDistance, $ComplexityParametersMap{MaxDistance} ...\n";
    }
  }

  # Set default parameter values...

  if (IsEmpty($ComplexityParametersMap{MolecularComplexityType})) {
    $ComplexityParametersMap{MolecularComplexityType} = 'MACCSKeys';
  }
  $MolecularComplexityType = $ComplexityParametersMap{MolecularComplexityType};


  if (IsEmpty($ComplexityParametersMap{AtomIdentifierType})) {
    $ComplexityParametersMap{AtomIdentifierType} = ($MolecularComplexityType =~ /^(TopologicalPharmacophoreAtomPairsFingerprints|TopologicalPharmacophoreAtomTripletsFingerprints)$/i) ? "FunctionalClassAtomTypes" : "AtomicInvariantsAtomTypes";
  }

  if (IsEmpty($ComplexityParametersMap{AtomicInvariantsToUse})) {
    my($AtomicInvariantsToUse);

    if ($MolecularComplexityType =~ /^(AtomTypesFingerprints|TopologicalAtomPairsFingerprints|TopologicalAtomTripletsFingerprints|TopologicalAtomTorsionsFingerprints)$/i) {
      $AtomicInvariantsToUse = "AS,X,BO,H,FC";
    }
    elsif ($MolecularComplexityType =~ /^ExtendedConnectivityFingerprints$/i) {
      $AtomicInvariantsToUse = "AS,X,BO,H,FC,MN";
    }
    else {
      $AtomicInvariantsToUse = "AS";
    }
    $ComplexityParametersMap{AtomicInvariantsToUse} = $AtomicInvariantsToUse;
  }

  if (IsEmpty($ComplexityParametersMap{FunctionalClassesToUse})) {
    my($FunctionalClassesToUse);

    if ($MolecularComplexityType =~ /^TopologicalPharmacophoreAtomPairsFingerprints$/i) {
      $FunctionalClassesToUse = "HBD,HBA,PI,NI,H";
    }
    elsif ($MolecularComplexityType =~ /^TopologicalPharmacophoreAtomTripletsFingerprints$/i) {
      $FunctionalClassesToUse = "HBD,HBA,PI,NI,H,Ar";
    }
    else {
      $FunctionalClassesToUse = "HBD,HBA,PI,NI,H,Ar,Hal";
    }
    $ComplexityParametersMap{FunctionalClassesToUse} = $FunctionalClassesToUse;
  }

  my(@AtomicInvariantsToUse);
  @AtomicInvariantsToUse = split ',', $ComplexityParametersMap{AtomicInvariantsToUse};
  $ComplexityParametersMap{AtomicInvariantsToUse} = \@AtomicInvariantsToUse;

  my(@FunctionalClassesToUse);
  @FunctionalClassesToUse = split ',', $ComplexityParametersMap{FunctionalClassesToUse};
  $ComplexityParametersMap{FunctionalClassesToUse} = \@FunctionalClassesToUse;

  if (IsEmpty($ComplexityParametersMap{UseTriangleInequality})) {
    $ComplexityParametersMap{UseTriangleInequality} = 0;
    if ($MolecularComplexityType =~ /^TopologicalPharmacophoreAtomTripletsFingerprints$/i) {
      $ComplexityParametersMap{UseTriangleInequality} = 1;
    }
  }

  if ($MolecularComplexityType =~ /^(TopologicalPharmacophoreAtomPairsFingerprints|TopologicalPharmacophoreAtomTripletsFingerprints)$/i) {
    if ($ComplexityParametersMap{AtomIdentifierType} !~ /^FunctionalClassAtomTypes$/i) {
      die "Error: The parameter value specified for AtomIdentifierType, $ComplexityParametersMap{AtomIdentifierType}, in option \"--MolecularComplexity\" is not valid for MolecularComplexityType, $MolecularComplexityType: Allowed value: FunctionalClassAtomTypes...\n";
    }
  }

  # Set up approprate paremeter names for specified molecular complexity...

  @ParameterNames = ();
  push @ParameterNames, 'MolecularComplexityType';

  @AtomIdentifierTypeParameters = ();
  push @AtomIdentifierTypeParameters, 'AtomIdentifierType';
  if ($ComplexityParametersMap{AtomIdentifierType} =~ /^AtomicInvariantsAtomTypes$/i) {
    push @AtomIdentifierTypeParameters, 'AtomicInvariantsToUse';
  }
  elsif ($ComplexityParametersMap{AtomIdentifierType} =~ /^FunctionalClassAtomTypes$/i) {
    push @AtomIdentifierTypeParameters, 'FunctionalClassesToUse';
  }

  COMPLEXITYTYPE: {
    if ($MolecularComplexityType =~ /^AtomTypesFingerprints$/i) {
      push @ParameterNames, @AtomIdentifierTypeParameters;
      last COMPLEXITYTYPE;
    }
    if ($MolecularComplexityType =~ /^ExtendedConnectivityFingerprints$/i) {
      push @ParameterNames, @AtomIdentifierTypeParameters;
      push @ParameterNames, ('NeighborhoodRadius', 'NormalizationMethodology');
      last COMPLEXITYTYPE;
    }
    if ($MolecularComplexityType =~ /^MACCSKeys$/i) {
      push @ParameterNames, 'MACCSKeysSize';
      last COMPLEXITYTYPE;
    }
    if ($MolecularComplexityType =~ /^PathLengthFingerprints$/i) {
      push @ParameterNames, @AtomIdentifierTypeParameters;
      push @ParameterNames, ('MinPathLength', 'MaxPathLength', 'UseBondSymbols');
      last COMPLEXITYTYPE;
    }
    if ($MolecularComplexityType =~ /^TopologicalAtomPairsFingerprints$/i) {
      push @ParameterNames, @AtomIdentifierTypeParameters;
      push @ParameterNames, ('MinDistance', 'MaxDistance');
      last COMPLEXITYTYPE;
    }
    if ($MolecularComplexityType =~ /^TopologicalAtomTripletsFingerprints$/i) {
      push @ParameterNames, @AtomIdentifierTypeParameters;
      push @ParameterNames, ('MinDistance', 'MaxDistance', 'UseTriangleInequality');
      last COMPLEXITYTYPE;
    }
    if ($MolecularComplexityType =~ /^TopologicalAtomTorsionsFingerprints$/i) {
      push @ParameterNames, @AtomIdentifierTypeParameters;
      last COMPLEXITYTYPE;
    }
    if ($MolecularComplexityType =~ /^TopologicalPharmacophoreAtomPairsFingerprints$/i) {
      push @ParameterNames, ('AtomIdentifierType', 'FunctionalClassesToUse', 'MinDistance', 'MaxDistance', 'NormalizationMethodology');
      last COMPLEXITYTYPE;
    }
    if ($MolecularComplexityType =~ /^TopologicalPharmacophoreAtomTripletsFingerprints$/i) {
      push @ParameterNames, ('AtomIdentifierType', 'FunctionalClassesToUse', 'MinDistance', 'MaxDistance', 'UseTriangleInequality', 'NormalizationMethodology', 'DistanceBinSize');
      last COMPLEXITYTYPE;
    }
    die "Error: The parameter value specified, $ParameterValue, for parameter name MolecularComplexityType using \"--MolecularComplexity\" is not valid.\n";
  }

  $OptionsInfo{MolecularComplexity} = $Options{molecularcomplexity};

  %{$OptionsInfo{MolecularComplexityParametersMap}} = ();
  for $ParameterName (@ParameterNames) {
    $ParameterValue = $ComplexityParametersMap{$ParameterName};
    $OptionsInfo{MolecularComplexityParametersMap}{$ParameterName} = $ParameterValue;
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{aromaticitymodel} = 'MayaChemToolsAromaticityModel';

  $Options{compoundidmode} = 'LabelPrefix';
  $Options{compoundidlabel} = 'CompoundID';
  $Options{datafieldsmode} = 'CompoundID';

  $Options{filter} = 'Yes';

  $Options{hydrogenbonds} = 'HBondsType2';

  $Options{keeplargestcomponent} = 'Yes';

  # Default mode values are set later...
  $Options{mode} = '';

  # Default moelcular complexity values are set later...
  $Options{molecularcomplexity} = '';

  # Default precision values are set later...
  $Options{precision} = '';

  $Options{output} = 'text';
  $Options{outdelim} = 'comma';
  $Options{quote} = 'yes';

  # Default rotatable bond parameter values are set later...
  $Options{rotatablebonds} = '';

  $Options{ruleof3violations} = 'No';
  $Options{ruleof5violations} = 'No';

  # Default TPSA paramater values are set later...
  $Options{tpsa} = '';

  if (!GetOptions(\%Options, "aromaticitymodel=s", "compoundid=s", "compoundidlabel=s", "compoundidmode=s", "datafields=s", "datafieldsmode|d=s", "filter|f=s", "help|h", "hydrogenbonds=s", "keeplargestcomponent|k=s", "mode|m=s", "molecularcomplexity=s", "outdelim=s", "output=s", "overwrite|o", "precision=s", "rotatablebonds=s", "ruleof3violations=s", "ruleof5violations=s", "quote|q=s", "root|r=s", "tpsa=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if (!Molecule::IsSupportedAromaticityModel($Options{aromaticitymodel})) {
    my(@SupportedModels) = Molecule::GetSupportedAromaticityModels();
    die "Error: The value specified, $Options{aromaticitymodel}, for option \"--AromaticityModel\" is not valid. Supported aromaticity models in current release of MayaChemTools: @SupportedModels\n";
  }
  if ($Options{compoundidmode} !~ /^(DataField|MolName|LabelPrefix|MolNameOrLabelPrefix)$/i) {
    die "Error: The value specified, $Options{compoundidmode}, for option \"--CompoundIDMode\" is not valid. Allowed values: DataField, MolName, LabelPrefix or MolNameOrLabelPrefix\n";
  }
  if ($Options{datafieldsmode} !~ /^(All|Common|Specify|CompoundID)$/i) {
    die "Error: The value specified, $Options{datafieldsmode}, for option \"-d, --DataFieldsMode\" is not valid. Allowed values: All, Common, Specify or CompoundID\n";
  }
  if ($Options{filter} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{filter}, for option \"-f, --Filter\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{hydrogenbonds} !~ /^(HBondsType1|HydrogenBondsType1|HBondsType2|HydrogenBondsType2)$/i) {
    die "Error: The value specified, $Options{hydrogenbonds}, for option \"--HydrogenBonds\" is not valid. Allowed values: HBondsType1, HydrogenBondsType1, HBondsType2, HydrogenBondsType2\n";
  }
  if ($Options{keeplargestcomponent} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{keeplargestcomponent}, for option \"-k, --KeepLargestComponent\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{output} !~ /^(SD|text|both)$/i) {
    die "Error: The value specified, $Options{output}, for option \"--output\" is not valid. Allowed values: SD, text, or both\n";
  }
  if ($Options{outdelim} !~ /^(comma|semicolon|tab)$/i) {
    die "Error: The value specified, $Options{outdelim}, for option \"--outdelim\" is not valid. Allowed values: comma, tab, or semicolon\n";
  }
  if ($Options{quote} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{quote}, for option \"-q --quote\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{ruleof3violations} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{ruleof3violations}, for option \"--RuleOf3Violations\" is not valid. Allowed values: Yes or No\n";
  }
  if ($Options{ruleof5violations} !~ /^(Yes|No)$/i) {
    die "Error: The value specified, $Options{ruleof5violations}, for option \"--RuleOf5Violations\" is not valid. Allowed values: Yes or No\n";
  }
}

__END__

=head1 NAME

CalculatePhysicochemicalProperties.pl - Calculate physicochemical properties for SD files

=head1 SYNOPSIS

CalculatePhysicochemicalProperties.pl SDFile(s)...

PhysicochemicalProperties.pl  [B<--AromaticityModel> I<AromaticityModelType>]
[B<--CompoundID> DataFieldName or LabelPrefixString]
[B<--CompoundIDLabel> text] [B<--CompoundIDMode>] [B<--DataFields> "FieldLabel1, FieldLabel2,..."]
[B<-d, --DataFieldsMode> All | Common | Specify | CompoundID] [B<-f, --Filter> Yes | No] [B<-h, --help>]
[B<--HydrogenBonds> HBondsType1 | HBondsType2] [B<-k, --KeepLargestComponent> Yes | No]
[B<-m, --mode> All | RuleOf5 | RuleOf3 | "name1, [name2,...]"]
[B<--MolecularComplexity> I<Name,Value, [Name,Value,...]>]
[B<--OutDelim> comma | tab | semicolon] [B<--output> SD | text | both] [B<-o, --overwrite>]
[B<--Precision> Name,Number,[Name,Number,..]] [B<--RotatableBonds> Name,Value, [Name,Value,...]]
[B<--RuleOf3Violations> Yes | No] [B<--RuleOf5Violations> Yes | No]
[B<-q, --quote> Yes | No] [B<-r, --root> RootName]
[B<-w, --WorkingDir> dirname] SDFile(s)...

=head1 DESCRIPTION

Calculate physicochemical properties for I<SDFile(s)> and create appropriate SD or CSV/TSV
text file(s) containing calculated properties.

The current release of MayaChemTools supports the calculation of these physicochemical
properties:

    MolecularWeight, ExactMass, HeavyAtoms, Rings, AromaticRings,
    van der Waals MolecularVolume [ Ref 93 ], RotatableBonds,
    HydrogenBondDonors, HydrogenBondAcceptors, LogP and
    Molar Refractivity (SLogP and SMR) [ Ref 89 ], Topological Polar
    Surface Area (TPSA) [ Ref 90 ], Fraction of SP3 carbons (Fsp3Carbons)
    and SP3 carbons (Sp3Carbons) [ Ref 115-116, Ref 119 ],
    MolecularComplexity [ Ref 117-119 ]

Multiple SDFile names are separated by spaces. The valid file extensions are I<.sdf>
and I<.sd>. All other file names are ignored. All the SD files in a current directory
can be specified either by I<*.sdf> or the current directory name.

The calculation of molecular complexity using I<MolecularComplexityType> parameter
corresponds to the number of bits-set or unique keys [ Ref 117-119 ] in molecular  fingerprints.
Default value for I<MolecularComplexityType>: I<MACCSKeys> of size 166. The calculation
of MACCSKeys is relatively expensive and can take rather substantial amount of time.

=head1 OPTIONS

=over 4

=item B<--AromaticityModel> I<MDLAromaticityModel | TriposAromaticityModel | MMFFAromaticityModel | ChemAxonBasicAromaticityModel | ChemAxonGeneralAromaticityModel | DaylightAromaticityModel | MayaChemToolsAromaticityModel>

Specify aromaticity model to use during detection of aromaticity. Possible values in the current
release are: I<MDLAromaticityModel, TriposAromaticityModel, MMFFAromaticityModel,
ChemAxonBasicAromaticityModel, ChemAxonGeneralAromaticityModel, DaylightAromaticityModel
or MayaChemToolsAromaticityModel>. Default value: I<MayaChemToolsAromaticityModel>.

The supported aromaticity model names along with model specific control parameters
are defined in B<AromaticityModelsData.csv>, which is distributed with the current release
and is available under B<lib/data> directory. B<Molecule.pm> module retrieves data from
this file during class instantiation and makes it available to method B<DetectAromaticity>
for detecting aromaticity corresponding to a specific model.

=item B<--CompoundID> I<DataFieldName or LabelPrefixString>

This value is B<--CompoundIDMode> specific and indicates how compound ID is generated.

For I<DataField> value of B<--CompoundIDMode> option, it corresponds to datafield label name
whose value is used as compound ID; otherwise, it's a prefix string used for generating compound
IDs like LabelPrefixString<Number>. Default value, I<Cmpd>, generates compound IDs which
look like Cmpd<Number>.

Examples for I<DataField> value of B<--CompoundIDMode>:

    MolID
    ExtReg

Examples for I<LabelPrefix> or I<MolNameOrLabelPrefix> value of B<--CompoundIDMode>:

    Compound

The value specified above generates compound IDs which correspond to Compound<Number>
instead of default value of Cmpd<Number>.

=item B<--CompoundIDLabel> I<text>

Specify compound ID column label for CSV/TSV text file(s) used during I<CompoundID> value
of B<--DataFieldsMode> option. Default value: I<CompoundID>.

=item B<--CompoundIDMode> I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>

Specify how to generate compound IDs and write to CSV/TSV text file(s) along with calculated
physicochemical properties for I<text | both> values of B<--output> option: use a I<SDFile(s)>
datafield value; use molname line from I<SDFile(s)>; generate a sequential ID with specific prefix;
use combination of both MolName and LabelPrefix with usage of LabelPrefix values for empty
molname lines.

Possible values: I<DataField | MolName | LabelPrefix | MolNameOrLabelPrefix>.
Default value: I<LabelPrefix>.

For I<MolNameAndLabelPrefix> value of B<--CompoundIDMode>, molname line in I<SDFile(s)> takes
precedence over sequential compound IDs generated using I<LabelPrefix> and only empty molname
values are replaced with sequential compound IDs.

This is only used for I<CompoundID> value of B<--DataFieldsMode> option.

=item B<--DataFields> I<"FieldLabel1,FieldLabel2,...">

Comma delimited list of I<SDFiles(s)> data fields to extract and write to CSV/TSV text file(s) along
with calculated physicochemical properties for I<text | both> values of B<--output> option.

This is only used for I<Specify> value of B<--DataFieldsMode> option.

Examples:

    Extreg
    MolID,CompoundName

=item B<-d, --DataFieldsMode> I<All | Common | Specify | CompoundID>

Specify how data fields in I<SDFile(s)> are transferred to output CSV/TSV text file(s) along
with calculated physicochemical properties for I<text | both> values of B<--output> option:
transfer all SD data field; transfer SD data files common to all compounds; extract specified
data fields; generate a compound ID using molname line, a compound prefix, or a combination
of both. Possible values: I<All | Common | specify | CompoundID>. Default value: I<CompoundID>.

=item B<-f, --Filter> I<Yes | No>

Specify whether to check and filter compound data in SDFile(s). Possible values: I<Yes or No>.
Default value: I<Yes>.

By default, compound data is checked before calculating physiochemical properties and compounds
containing atom data corresponding to non-element symbols or no atom data are ignored.

=item B<-h, --help>

Print this help message.

=item B<--HydrogenBonds> I<HBondsType1 | HBondsType2>

Parameters to control calculation of hydrogen bond donors and acceptors. Possible values:
I<HBondsType1, HydrogenBondsType1, HBondsType2, HydrogenBondsType2>. Default value:
I<HBondsType2> which corresponds to B<RuleOf5> definition for number of hydrogen bond
donors and acceptors.

The current release of MayaChemTools supports identification of two types of hydrogen bond
donor and acceptor atoms with these names:

    HBondsType1 or HydrogenBondsType1
    HBondsType2 or HydrogenBondsType2

The names of these hydrogen bond types are rather arbitrary. However, their definitions have
specific meaning and are as follows:

    HydrogenBondsType1 [ Ref 60-61, Ref 65-66 ]:

        Donor: NH, NH2, OH - Any N and O with available H
        Acceptor: N[!H], O - Any N without available H and any O

    HydrogenBondsType2 [ Ref 91 ]:

        Donor: NH, NH2, OH - N and O with available H
        Acceptor: N, O - And N and O

=item B<-k, --KeepLargestComponent> I<Yes | No>

Calculate physicochemical properties for only the largest component in molecule. Possible values:
I<Yes or No>. Default value: I<Yes>.

For molecules containing multiple connected components, physicochemical properties can be
calculated in two different ways: use all connected components or just the largest connected
component. By default, all atoms except for the largest connected component are
deleted before calculation of physicochemical properties.

=item B<-m, --mode> I<All | RuleOf5 | RuleOf3 | "name1, [name2,...]">

Specify physicochemical properties to calculate for SDFile(s): calculate all available physical
chemical properties; calculate properties corresponding to Rule of 5; or use a comma delimited
list of supported physicochemical properties. Possible values: I<All | RuleOf5 | RuleOf3 |
"name1, [name2,...]">.

Default value: I<MolecularWeight, HeavyAtoms, MolecularVolume, RotatableBonds, HydrogenBondDonors,
HydrogenBondAcceptors, SLogP, TPSA>. These properties are calculated by default.

I<RuleOf5> [ Ref 91 ] includes these properties: I<MolecularWeight, HydrogenBondDonors, HydrogenBondAcceptors,
SLogP>. I<RuleOf5> states: MolecularWeight <= 500, HydrogenBondDonors <= 5, HydrogenBondAcceptors <= 10, and
logP <= 5.

I<RuleOf3> [ Ref 92 ] includes these properties: I<MolecularWeight, RotatableBonds, HydrogenBondDonors,
HydrogenBondAcceptors, SLogP, TPSA>. I<RuleOf3> states: MolecularWeight <= 300, RotatableBonds <= 3,
HydrogenBondDonors <= 3, HydrogenBondAcceptors <= 3, logP <= 3, and TPSA <= 60.

I<All> calculates all supported physicochemical properties: I<MolecularWeight, ExactMass,
HeavyAtoms, Rings, AromaticRings, MolecularVolume, RotatableBonds, HydrogenBondDonors,
HydrogenBondAcceptors, SLogP, SMR, TPSA, Fsp3Carbons, Sp3Carbons, MolecularComplexity>.

=item B<--MolecularComplexity> I<Name,Value, [Name,Value,...]>

Parameters to control calculation of molecular complexity: it's a comma delimited list of parameter
name and value pairs.

Possible parameter names: I<MolecularComplexityType, AtomIdentifierType,
AtomicInvariantsToUse, FunctionalClassesToUse, MACCSKeysSize, NeighborhoodRadius,
MinPathLength, MaxPathLength, UseBondSymbols, MinDistance, MaxDistance,
UseTriangleInequality, DistanceBinSize, NormalizationMethodology>.

The valid paramater valuse for each parameter name are described in the following sections.

The current release of MayaChemTools supports calculation of molecular complexity using
I<MolecularComplexityType> parameter corresponding to the number of bits-set or unique
keys [ Ref 117-119 ] in molecular  fingerprints. The valid values for I<MolecularComplexityType>
are:

    AtomTypesFingerprints
    ExtendedConnectivityFingerprints
    MACCSKeys
    PathLengthFingerprints
    TopologicalAtomPairsFingerprints
    TopologicalAtomTripletsFingerprints
    TopologicalAtomTorsionsFingerprints
    TopologicalPharmacophoreAtomPairsFingerprints
    TopologicalPharmacophoreAtomTripletsFingerprints

Default value for I<MolecularComplexityType>: I<MACCSKeys>.

I<AtomIdentifierType> parameter name correspods to atom types used during generation of
fingerprints. The valid values for I<AtomIdentifierType> are: I<AtomicInvariantsAtomTypes,
DREIDINGAtomTypes, EStateAtomTypes, FunctionalClassAtomTypes, MMFF94AtomTypes,
SLogPAtomTypes, SYBYLAtomTypes, TPSAAtomTypes, UFFAtomTypes>. I<AtomicInvariantsAtomTypes>
is not supported for during the following values of I<MolecularComplexityType>: I<MACCSKeys,
TopologicalPharmacophoreAtomPairsFingerprints, TopologicalPharmacophoreAtomTripletsFingerprints>.
I<FunctionalClassAtomTypes> is the only valid value for I<AtomIdentifierType> for topological
pharmacophore fingerprints.

Default value for I<AtomIdentifierType>: I<AtomicInvariantsAtomTypes>
for all except topological pharmacophore fingerprints where it is I<FunctionalClassAtomTypes>.

I<AtomicInvariantsToUse> parameter name and values are used during I<AtomicInvariantsAtomTypes>
value of parameter I<AtomIdentifierType>. It's a list of space separated valid atomic invariant atom types.

Possible values for atomic invariants are: I<AS, X, BO,  LBO, SB, DB, TB, H, Ar, RA, FC, MN, SM>.
Default value for I<AtomicInvariantsToUse> parameter are set differently for different fingerprints
using I<MolecularComplexityType> parameter as shown below:

    MolecularComplexityType              AtomicInvariantsToUse

    AtomTypesFingerprints                AS X BO H FC
    TopologicalAtomPairsFingerprints     AS X BO H FC
    TopologicalAtomTripletsFingerprints  AS X BO H FC
    TopologicalAtomTorsionsFingerprints  AS X BO H FC

    ExtendedConnectivityFingerprints     AS X  BO H FC MN
    PathLengthFingerprints               AS


The atomic invariants abbreviations correspond to:

    AS = Atom symbol corresponding to element symbol

    X<n>   = Number of non-hydrogen atom neighbors or heavy atoms
    BO<n> = Sum of bond orders to non-hydrogen atom neighbors or heavy atoms
    LBO<n> = Largest bond order of non-hydrogen atom neighbors or heavy atoms
    SB<n> = Number of single bonds to non-hydrogen atom neighbors or heavy atoms
    DB<n> = Number of double bonds to non-hydrogen atom neighbors or heavy atoms
    TB<n> = Number of triple bonds to non-hydrogen atom neighbors or heavy atoms
    H<n>   = Number of implicit and explicit hydrogens for atom
    Ar     = Aromatic annotation indicating whether atom is aromatic
    RA     = Ring atom annotation indicating whether atom is a ring
    FC<+n/-n> = Formal charge assigned to atom
    MN<n> = Mass number indicating isotope other than most abundant isotope
    SM<n> = Spin multiplicity of atom. Possible values: 1 (singlet), 2 (doublet) or
            3 (triplet)

Atom type generated by AtomTypes::AtomicInvariantsAtomTypes class corresponds to:

    AS.X<n>.BO<n>.LBO<n>.<SB><n>.<DB><n>.<TB><n>.H<n>.Ar.RA.FC<+n/-n>.MN<n>.SM<n>

Except for AS which is a required atomic invariant in atom types, all other atomic invariants are
optional. Atom type specification doesn't include atomic invariants with zero or undefined values.

In addition to usage of abbreviations for specifying atomic invariants, the following descriptive words
are also allowed:

    X : NumOfNonHydrogenAtomNeighbors or NumOfHeavyAtomNeighbors
    BO : SumOfBondOrdersToNonHydrogenAtoms or SumOfBondOrdersToHeavyAtoms
    LBO : LargestBondOrderToNonHydrogenAtoms or LargestBondOrderToHeavyAtoms
    SB :  NumOfSingleBondsToNonHydrogenAtoms or NumOfSingleBondsToHeavyAtoms
    DB : NumOfDoubleBondsToNonHydrogenAtoms or NumOfDoubleBondsToHeavyAtoms
    TB : NumOfTripleBondsToNonHydrogenAtoms or NumOfTripleBondsToHeavyAtoms
    H :  NumOfImplicitAndExplicitHydrogens
    Ar : Aromatic
    RA : RingAtom
    FC : FormalCharge
    MN : MassNumber
    SM : SpinMultiplicity

I<AtomTypes::AtomicInvariantsAtomTypes> module is used to assign atomic invariant
atom types.

I<FunctionalClassesToUse> parameter name and values are used during I<FunctionalClassAtomTypes>
value of parameter I<AtomIdentifierType>. It's a list of space separated valid atomic invariant atom types.

Possible values for atom functional classes are: I<Ar, CA, H, HBA, HBD, Hal, NI, PI, RA>.

Default value for I<FunctionalClassesToUse> parameter is set to:

    HBD HBA PI NI Ar Hal

for all fingerprints except for the following two I<MolecularComplexityType> fingerints:

    MolecularComplexityType                           FunctionalClassesToUse

    TopologicalPharmacophoreAtomPairsFingerprints     HBD HBA P, NI H
    TopologicalPharmacophoreAtomTripletsFingerprints  HBD HBA PI NI H Ar

The functional class abbreviations correspond to:

    HBD: HydrogenBondDonor
    HBA: HydrogenBondAcceptor
    PI :  PositivelyIonizable
    NI : NegativelyIonizable
    Ar : Aromatic
    Hal : Halogen
    H : Hydrophobic
    RA : RingAtom
    CA : ChainAtom

 Functional class atom type specification for an atom corresponds to:

    Ar.CA.H.HBA.HBD.Hal.NI.PI.RA

I<AtomTypes::FunctionalClassAtomTypes> module is used to assign functional class atom
types. It uses following definitions [ Ref 60-61, Ref 65-66 ]:

    HydrogenBondDonor: NH, NH2, OH
    HydrogenBondAcceptor: N[!H], O
    PositivelyIonizable: +, NH2
    NegativelyIonizable: -, C(=O)OH, S(=O)OH, P(=O)OH

I<MACCSKeysSize> parameter name is only used during I<MACCSKeys> value of
I<MolecularComplexityType> and corresponds to the size of MACCS key set. Possible
values: I<166 or 322>. Default value: I<166>.

I<NeighborhoodRadius> parameter name is only used during I<ExtendedConnectivityFingerprints>
value of I<MolecularComplexityType> and corresponds to atomic neighborhoods radius for
generating extended connectivity fingerprints. Possible values: positive integer. Default value:
I<2>.

I<MinPathLength> and I<MaxPathLength> parameters are only used during I<PathLengthFingerprints>
value of I<MolecularComplexityType> and correspond to minimum and maximum path lengths to use
for generating path length fingerprints. Possible values: positive integers. Default value: I<MinPathLength - 1>;
I<MaxPathLength - 8>.

I<UseBondSymbols> parameter is only used during I<PathLengthFingerprints> value of
I<MolecularComplexityType> and indicates whether bond symbols are included in atom path
strings used to generate path length fingerprints. Possible value: I<Yes or No>. Default value:
I<Yes>.

I<MinDistance> and I<MaxDistance> parameters are only used during I<TopologicalAtomPairsFingerprints>
and I<TopologicalAtomTripletsFingerprints> values of I<MolecularComplexityType> and correspond to
minimum and maximum bond distance between atom pairs during topological pharmacophore fingerprints.
Possible values: positive integers. Default value: I<MinDistance - 1>; I<MaxDistance - 10>.

I<UseTriangleInequality> parameter is used during these values for I<MolecularComplexityType>:
I<TopologicalAtomTripletsFingerprints> and I<TopologicalPharmacophoreAtomTripletsFingerprints>.
Possible values: I<Yes or No>. It determines wheter to apply triangle inequality to distance triplets.
Default value: I<TopologicalAtomTripletsFingerprints - No>;
I<TopologicalPharmacophoreAtomTripletsFingerprints - Yes>.

I<DistanceBinSize> parameter is used during I<TopologicalPharmacophoreAtomTripletsFingerprints>
value of I<MolecularComplexityType> and correspons to distance bin size used for binning
distances during generation of topological pharmacophore atom triplets fingerprints. Possible
value: positive integer. Default value: I<2>.

I<NormalizationMethodology> is only used for these values for I<MolecularComplexityType>:
I<ExtendedConnectivityFingerprints>, I<TopologicalPharmacophoreAtomPairsFingerprints>
and I<TopologicalPharmacophoreAtomTripletsFingerprints>. It corresponds to normalization
methodology to use for scaling the number of bits-set or unique keys during generation of
fingerprints. Possible values during I<ExtendedConnectivityFingerprints>: I<None or
ByHeavyAtomsCount>; Default value: I<None>. Possible values during topological
pharmacophore atom pairs and tripletes fingerprints: I<None or ByPossibleKeysCount>;
Default value: I<None>. I<ByPossibleKeysCount> corresponds to total number of
possible topological pharmacophore atom pairs or triplets in a molecule.

Examples of I<MolecularComplexity> name and value parameters:

    MolecularComplexityType,AtomTypesFingerprints,AtomIdentifierType,
    AtomicInvariantsAtomTypes,AtomicInvariantsToUse,AS X BO H FC

    MolecularComplexityType,ExtendedConnectivityFingerprints,
    AtomIdentifierType,AtomicInvariantsAtomTypes,
    AtomicInvariantsToUse,AS X BO H FC MN,NeighborhoodRadius,2,
    NormalizationMethodology,None

    MolecularComplexityType,MACCSKeys,MACCSKeysSize,166

    MolecularComplexityType,PathLengthFingerprints,AtomIdentifierType,
    AtomicInvariantsAtomTypes,AtomicInvariantsToUse,AS,MinPathLength,
    1,MaxPathLength,8,UseBondSymbols,Yes

    MolecularComplexityType,TopologicalAtomPairsFingerprints,
    AtomIdentifierType,AtomicInvariantsAtomTypes,AtomicInvariantsToUse,
    AS X BO H FC,MinDistance,1,MaxDistance,10

    MolecularComplexityType,TopologicalAtomTripletsFingerprints,
    AtomIdentifierType,AtomicInvariantsAtomTypes,AtomicInvariantsToUse,
    AS X BO H FC,MinDistance,1,MaxDistance,10,UseTriangleInequality,No

    MolecularComplexityType,TopologicalAtomTorsionsFingerprints,
    AtomIdentifierType,AtomicInvariantsAtomTypes,AtomicInvariantsToUse,
    AS X BO H FC

    MolecularComplexityType,TopologicalPharmacophoreAtomPairsFingerprints,
    AtomIdentifierType,FunctionalClassAtomTypes,FunctionalClassesToUse,
    HBD HBA PI NI H,MinDistance,1,MaxDistance,10,NormalizationMethodology,
    None

    MolecularComplexityType,TopologicalPharmacophoreAtomTripletsFingerprints,
    AtomIdentifierType,FunctionalClassAtomTypes,FunctionalClassesToUse,
    HBD HBA PI NI H Ar,MinDistance,1,MaxDistance,10,NormalizationMethodology,
    None,UseTriangleInequality,Yes,NormalizationMethodology,None,
    DistanceBinSize,2

=item B<--OutDelim> I<comma | tab | semicolon>

Delimiter for output CSV/TSV text file(s). Possible values: I<comma, tab, or semicolon>
Default value: I<comma>.

=item B<--output> I<SD | text | both>

Type of output files to generate. Possible values: I<SD, text, or both>. Default value: I<text>.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<--Precision> I<Name,Number,[Name,Number,..]>

Precision of calculated property values in the output file: it's a comma delimited list of
property name and precision value pairs. Possible property names: I<MolecularWeight,
ExactMass>. Possible values: positive intergers. Default value: I<MolecularWeight,2,
ExactMass,4>.

Examples:

    ExactMass,3
    MolecularWeight,1,ExactMass,2

=item B<-q, --quote> I<Yes | No>

Put quote around column values in output CSV/TSV text file(s). Possible values:
I<Yes or No>. Default value: I<Yes>.

=item B<-r, --root> I<RootName>

New file name is generated using the root: <Root>.<Ext>. Default for new file names:
<SDFileName><PhysicochemicalProperties>.<Ext>. The file type determines <Ext> value.
The sdf, csv, and tsv <Ext> values are used for SD, comma/semicolon, and tab
delimited text files, respectively.This option is ignored for multiple input files.

=item B<--RotatableBonds> I<Name,Value, [Name,Value,...]>

Parameters to control calculation of rotatable bonds [ Ref 92 ]: it's a comma delimited list of parameter
name and value pairs. Possible parameter names: I<IgnoreTerminalBonds, IgnoreBondsToTripleBonds,
IgnoreAmideBonds, IgnoreThioamideBonds, IgnoreSulfonamideBonds>. Possible parameter values:
I<Yes or No>. By default, value of all parameters is set to I<Yes>.

=item B<--RuleOf3Violations> I<Yes | No>

Specify whether to calculate B<RuleOf3Violations> for SDFile(s). Possible values: I<Yes or No>.
Default value: I<No>.

For I<Yes> value of B<RuleOf3Violations>, in addition to calculating total number of B<RuleOf3> violations,
individual violations for compounds are also written to output files.

B<RuleOf3> [ Ref 92 ] states: MolecularWeight <= 300, RotatableBonds <= 3, HydrogenBondDonors <= 3,
HydrogenBondAcceptors <= 3, logP <= 3, and TPSA <= 60.

=item B<--RuleOf5Violations> I<Yes | No>

Specify whether to calculate B<RuleOf5Violations> for SDFile(s). Possible values: I<Yes or No>.
Default value: I<No>.

For I<Yes> value of B<RuleOf5Violations>, in addition to calculating total number of B<RuleOf5> violations,
individual violations for compounds are also written to output files.

B<RuleOf5> [ Ref 91 ] states: MolecularWeight <= 500, HydrogenBondDonors <= 5, HydrogenBondAcceptors <= 10,
and logP <= 5.

=item B<--TPSA> I<Name,Value, [Name,Value,...]>

Parameters to control calculation of TPSA: it's a comma delimited list of parameter name and value
pairs. Possible parameter names: I<IgnorePhosphorus, IgnoreSulfur>. Possible parameter values:
I<Yes or No>. By default, value of all parameters is set to I<Yes>.

By default, TPSA atom contributions from Phosphorus and Sulfur atoms are not included during
TPSA calculations. [ Ref 91 ]

=item B<-w, --WorkingDir> I<DirName>

Location of working directory. Default value: current directory.

=back

=head1 EXAMPLES

To calculate default set of physicochemical properties - MolecularWeight, HeavyAtoms,
MolecularVolume, RotatableBonds, HydrogenBondDonor, HydrogenBondAcceptors, SLogP,
TPSA - and generate a SamplePhysicochemicalProperties.csv file containing sequential
compound IDs along with properties data, type:

    % CalculatePhysicochemicalProperties.pl -o Sample.sdf

To calculate all available physicochemical properties and generate both SampleAllProperties.csv
and SampleAllProperties.sdf files containing sequential compound IDs in CSV file along with
properties data, type:

    % CalculatePhysicochemicalProperties.pl -m All --output both
      -r SampleAllProperties -o Sample.sdf

To calculate RuleOf5 physicochemical properties and generate a SampleRuleOf5Properties.csv file
containing sequential compound IDs along with properties data, type:

    % CalculatePhysicochemicalProperties.pl -m RuleOf5
      -r SampleRuleOf5Properties -o Sample.sdf

To calculate RuleOf5 physicochemical properties along with counting RuleOf5 violations and generate
a SampleRuleOf5Properties.csv file containing sequential compound IDs along with properties data, type:

    % CalculatePhysicochemicalProperties.pl -m RuleOf5 --RuleOf5Violations Yes
      -r SampleRuleOf5Properties -o Sample.sdf

To calculate RuleOf3 physicochemical properties and generate a SampleRuleOf3Properties.csv file
containing sequential compound IDs along with properties data, type:

    % CalculatePhysicochemicalProperties.pl -m RuleOf3
      -r SampleRuleOf3Properties -o Sample.sdf

To calculate RuleOf3 physicochemical properties along with counting RuleOf3 violations and generate
a SampleRuleOf3Properties.csv file containing sequential compound IDs along with properties data, type:

    % CalculatePhysicochemicalProperties.pl -m RuleOf3 --RuleOf3Violations Yes
      -r SampleRuleOf3Properties -o Sample.sdf

To calculate a specific set of physicochemical properties and generate a SampleProperties.csv file
containing sequential compound IDs along with properties data, type:

    % CalculatePhysicochemicalProperties.pl -m "Rings,AromaticRings"
      -r SampleProperties -o Sample.sdf

To calculate HydrogenBondDonors and HydrogenBondAcceptors using HydrogenBondsType1 definition
and generate a SampleProperties.csv file containing sequential compound IDs along with properties
data, type:

    % CalculatePhysicochemicalProperties.pl -m "HydrogenBondDonors,HydrogenBondAcceptors"
      --HydrogenBonds HBondsType1 -r SampleProperties -o Sample.sdf

To calculate TPSA using sulfur and phosphorus atoms along with nitrogen and oxygen atoms and
generate a SampleProperties.csv file containing sequential compound IDs along with properties
data, type:

    % CalculatePhysicochemicalProperties.pl -m "TPSA" --TPSA "IgnorePhosphorus,No,
      IgnoreSulfur,No" -r SampleProperties -o Sample.sdf

To calculate MolecularComplexity using extendend connectivity fingerprints corresponding
to atom neighborhood radius of 2 with atomic invariant atom types without any scaling and
generate a SampleProperties.csv file containing sequential compound IDs along with properties
data, type:

    % CalculatePhysicochemicalProperties.pl -m MolecularComplexity --MolecularComplexity
      "MolecularComplexityType,ExtendedConnectivityFingerprints,NeighborhoodRadius,2,
      AtomIdentifierType, AtomicInvariantsAtomTypes,
      AtomicInvariantsToUse,AS X BO H FC MN,NormalizationMethodology,None"
      -r SampleProperties -o Sample.sdf

To calculate RuleOf5 physicochemical properties along with counting RuleOf5 violations and generate
a SampleRuleOf5Properties.csv file containing compound IDs from molecule name line along with
properties data, type:

    % CalculatePhysicochemicalProperties.pl -m RuleOf5 --RuleOf5Violations Yes
      --DataFieldsMode CompoundID --CompoundIDMode MolName
      -r SampleRuleOf5Properties -o Sample.sdf

To calculate all available physicochemical properties and generate a SampleAllProperties.csv
file containing compound ID using specified data field along with along with properties data,
type:

    % CalculatePhysicochemicalProperties.pl -m All
      --DataFieldsMode CompoundID --CompoundIDMode DataField --CompoundID Mol_ID
      -r SampleAllProperties -o Sample.sdf

To calculate all available physicochemical properties and generate a SampleAllProperties.csv
file containing compound ID using combination of molecule name line and an explicit compound
prefix along with properties data, type:

    % CalculatePhysicochemicalProperties.pl -m All
      --DataFieldsMode CompoundID --CompoundIDMode MolnameOrLabelPrefix
      --CompoundID Cmpd --CompoundIDLabel MolID  -r SampleAllProperties
      -o Sample.sdf

To calculate all available physicochemical properties and generate a SampleAllProperties.csv
file containing specific data fields columns along with with properties data, type:

    % CalculatePhysicochemicalProperties.pl -m All
      --DataFieldsMode Specify --DataFields Mol_ID -r SampleAllProperties
      -o Sample.sdf

To calculate all available physicochemical properties and generate a SampleAllProperties.csv
file containing common data fields columns along with with properties data, type:

    % CalculatePhysicochemicalProperties.pl -m All
      --DataFieldsMode Common -r SampleAllProperties -o Sample.sdf

To calculate all available physicochemical properties and generate both SampleAllProperties.csv
and CSV files containing all data fields columns in CSV files along with with properties data, type:

    % CalculatePhysicochemicalProperties.pl -m All
      --DataFieldsMode All  --output both -r SampleAllProperties
      -o Sample.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

ExtractFromSDtFiles.pl, ExtractFromTextFiles.pl, InfoSDFiles.pl, InfoTextFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
