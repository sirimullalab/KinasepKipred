#!/usr/bin/perl -w
#
# File: ModifySDFilesDataFields.pl
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
use SDFileUtil;
use TextUtil;

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

print "Checking input SD file(s)...\n";
my(%SDFilesInfo);
RetrieveSDFilesInfo();

# Generate output files...
my($FileIndex);
if (@SDFilesList > 1) {
  print "\nProcessing SD files...\n";
}
for $FileIndex (0 .. $#SDFilesList) {
  if ($SDFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $SDFilesList[$FileIndex]...\n";
    ModifySDFile($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Modify SD file data fields....
sub ModifySDFile {
  my($Index) = @_;
  my($SDFile, $NewSDFile);

  $SDFile = $SDFilesList[$Index];
  $NewSDFile = $SDFilesInfo{OutFile}[$Index];

  print "Generating new SD file $NewSDFile...\n";
  open NEWSDFILE, ">$NewSDFile" or die "Error: Couldn't open $NewSDFile: $! \n";
  open SDFILE, "$SDFile" or die "Error: Can't open $SDFile: $! \n";

  my($CmpdCount, $CmpdString, $CmpdData, $MolName, $OldSDField, $NewSDField, $CommonSDField, $Label, $Value, $FieldValues, $MolNameDataField, $URLCmpdIdFieldName, @CmpdLines, %DataFieldAndValues, @DataFieldLabels);
  $CmpdCount = 0;

  COMPOUND: while ($CmpdString = ReadCmpdString(\*SDFILE)) {
      $CmpdCount++;
      @CmpdLines = split "\n", $CmpdString;
      if ($OptionsInfo{UseDataFieldForMolName} || $OptionsInfo{ModifyDataFields}) {
	%DataFieldAndValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
      }
      if ($OptionsInfo{ModifyMolName}) {
	if ($OptionsInfo{AlwaysReplaceMolName} || !IsNotEmpty($CmpdLines[0])) {
	  $MolNameDataField = $OptionsInfo{MolNameDataField};
	  if ($OptionsInfo{UseDataFieldForMolName} && exists($DataFieldAndValues{$MolNameDataField})) {
	    $MolName = $DataFieldAndValues{$MolNameDataField};
	    if (length($MolName) > 80) {
	      $MolName = substr($MolName, 0, 80);
	    }
	  }
	  else {
	    $MolName = "$OptionsInfo{MolNamePrefix}${CmpdCount}";
	  }
	  $CmpdLines[0] = $MolName;
	  $CmpdString = join "\n", @CmpdLines;
	}
      }
      if (!$OptionsInfo{ModifyDataFields}) {
	# Just write the data and get the next compound...
	print NEWSDFILE "$CmpdString\n";
	next COMPOUND;
      }
      # Write out the structure data now and handle the old data fields later...
      ($CmpdData) = split /\n>/, $CmpdString;
      print NEWSDFILE "$CmpdData\n";

      # Modify specified data fields...
      for $NewSDField (sort keys %{$OptionsInfo{SpecifiedNewToOldSDFieldMap}}) {
	$FieldValues = "";
	for $OldSDField (@{$OptionsInfo{SpecifiedNewToOldSDFieldMap}{$NewSDField}}) {
	  if (exists($DataFieldAndValues{$OldSDField}) && length($DataFieldAndValues{$OldSDField})) {
	    $Value = $DataFieldAndValues{$OldSDField};
	    $FieldValues .= ($FieldValues) ? "\n$Value" : $Value;
	  }
	}
	print NEWSDFILE "> <$NewSDField>\n$FieldValues\n\n";
      }
      # Add specified common fields...
      for $CommonSDField (sort keys %{$OptionsInfo{SpecifiedCommonFieldMap}}) {
	$Value = $OptionsInfo{SpecifiedCommonFieldMap}{$CommonSDField};
	print NEWSDFILE "> <$CommonSDField>\n$Value\n\n";
      }
      if ($OptionsInfo{CreateDataFieldURL}) {
	$Value = "";
	$URLCmpdIdFieldName = $OptionsInfo{URLCmpdIdFieldName};
	if (exists($DataFieldAndValues{$URLCmpdIdFieldName}) && length($DataFieldAndValues{$URLCmpdIdFieldName})) {
	  $Value = $DataFieldAndValues{$URLCmpdIdFieldName};
	  $Value = "$OptionsInfo{URLCGIScriptName}?$OptionsInfo{URLParamName}=${Value}";
	}
	print NEWSDFILE "> <$OptionsInfo{URLDataFieldLabel}>\n$Value\n\n";
      }

      # Handle old data fields and write 'em in the same order as they appear in the input
      # files...
      if ($OptionsInfo{KeepAllOldDataFields} || $OptionsInfo{KeepUnMappedOldDataFields}) {
	my($KeepLabel);
	@DataFieldLabels = GetCmpdDataHeaderLabels(\@CmpdLines);
	LABEL: for $Label (@DataFieldLabels) {
	  $KeepLabel = $OptionsInfo{KeepAllOldDataFields} ? 1 : ( exists($OptionsInfo{SpecifiedOldToNewSDFieldMap}{$Label}) ? 0 : 1  );
	  if (!$KeepLabel) {
	    next LABEL;
	  }
	  $Value = $DataFieldAndValues{$Label};
	  print NEWSDFILE "> <$Label>\n$Value\n\n";
	}
      }

      print NEWSDFILE "\$\$\$\$\n";
  }
  close NEWSDFILE;
  close SDFILE;
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{Mode} = $Options{mode};

  $OptionsInfo{ModifyMolName} = 1; $OptionsInfo{ModifyDataFields} = 0;
  if ($Options{mode} =~ /^both$/i) {
    $OptionsInfo{ModifyMolName} = 1; $OptionsInfo{ModifyDataFields} = 1;
  }
  elsif ($Options{mode} =~ /^datafields$/i) {
    $OptionsInfo{ModifyMolName} = 0; $OptionsInfo{ModifyDataFields} = 1;
  }

  $OptionsInfo{KeepOldDataFields} = $Options{keepolddatafields};
  $OptionsInfo{KeepAllOldDataFields} = ($Options{keepolddatafields} =~ /^all$/i) ? 1 : 0;
  $OptionsInfo{KeepUnMappedOldDataFields} = ($Options{keepolddatafields} =~ /^unmappedonly$/i) ? 1 : 0;

  $OptionsInfo{MolNameMode} = $Options{molnamemode};
  $OptionsInfo{UseDataFieldForMolName} = ($Options{molnamemode} =~ /^datafield$/i) ? 1 : 0;

  $OptionsInfo{MolName} = $Options{molname};
  $OptionsInfo{MolNameDataField} = ""; $OptionsInfo{MolNamePrefix} = "Cmpd";
  if ($Options{molname}) {
    if ($OptionsInfo{UseDataFieldForMolName}) {
      $OptionsInfo{MolNameDataField} = $Options{molname};
    }
    else {
      $OptionsInfo{MolNamePrefix} = $Options{molname};
    }
  }

  $OptionsInfo{MolNameReplace} = $Options{molnamereplace};
  $OptionsInfo{AlwaysReplaceMolName} = ($Options{molnamereplace} =~ /^always$/i) ? 1 : 0;

  if ($Options{datafieldsmap} && $Options{datafieldsmapfile}) {
    die "Error: Both \"--datafieldsmap\" and  \"--datafieldsmapfile\" options specified: only one is allowed at a time\n";
  }

  $OptionsInfo{DataFieldsMap} = $Options{datafieldsmap} ? $Options{datafieldsmap} : '';
  $OptionsInfo{DataFieldsMapFile} = $Options{datafieldsmapfile} ? $Options{datafieldsmapfile} : '';

  my($SpecifiedDataFieldMap);

  %{$OptionsInfo{SpecifiedNewToOldSDFieldMap}} = ();
  %{$OptionsInfo{SpecifiedOldToNewSDFieldMap}} = ();

  $SpecifiedDataFieldMap = "";
  if ($Options{datafieldsmap}) {
    $SpecifiedDataFieldMap = $Options{datafieldsmap};
  }
  elsif ($Options{datafieldsmapfile}) {
    my($Line, @LineWords);
    open DATAFIELDSFILE, "$Options{datafieldsmapfile}" or die "Couldn't  open $Options{datafieldsmapfile}: $! \n";
    while ($Line = GetTextLine(\*DATAFIELDSFILE)) {
      @LineWords = quotewords(";", 0, $Line);
      $SpecifiedDataFieldMap .= JoinWords(\@LineWords, ";", 0);
    }
    close DATAFIELDSFILE;
  }

  if ($SpecifiedDataFieldMap) {
    my($DataFieldMap, $DataField, $NewSDField, @OldSDFields, @DataFieldMapSplit, @DataFieldsSplit, $FirstField);
    @DataFieldMapSplit = split ";", $SpecifiedDataFieldMap;
    for $DataFieldMap (@DataFieldMapSplit) {
      @DataFieldsSplit = split ",", $DataFieldMap;
      if (@DataFieldsSplit == 1) {
	die "Error: Invalid number of comma delimited values, ", scalar(@DataFieldsSplit), ", specified,  @DataFieldsSplit, using \"--datafieldsmap or --datafieldsmapfile\" option: it must contain more than one value.\n";
      }
      $FirstField = 1;
      @OldSDFields = ();
      for $DataField (@DataFieldsSplit) {
	if (!(defined($DataField) && length($DataField))) {
	  die "Error: One of the comma delimited values, \"", join(",", @DataFieldsSplit), "\", specified using \"--datafieldsmap or --datafieldsmapfile\" option is empty.\n";
	}
	if ($FirstField) {
	  $FirstField = 0;
	  $NewSDField = $DataField;
	}
	else {
	  push @OldSDFields, $DataField;
	}
      }
      # Make sure a datafield is only specified once...
      if (exists $OptionsInfo{SpecifiedNewToOldSDFieldMap}{$NewSDField}) {
	die "Error: New data field, $NewSDField, specified more than once using \"--datafieldsmap or --datafieldsmapfile\" option.\n";
      }
      @{$OptionsInfo{SpecifiedNewToOldSDFieldMap}{$NewSDField}} = ();
      push @{$OptionsInfo{SpecifiedNewToOldSDFieldMap}{$NewSDField}}, @OldSDFields;
      for $DataField (@OldSDFields) {
	if (exists $OptionsInfo{SpecifiedOldToNewSDFieldMap}{$DataField} ) {
	  die "Error: SD field, $DataField, specified more than once using \"--datafieldsmap or --datafieldsmapfile\" option.\n";
	}
	else {
	  $OptionsInfo{SpecifiedOldToNewSDFieldMap}{$DataField} = $NewSDField;
	}
      }

    }
  }

  $OptionsInfo{DataFieldsCommon} = $Options{datafieldscommon} ? $Options{datafieldscommon} : '';
  %{$OptionsInfo{SpecifiedCommonFieldMap}} = ();

  if ($Options{datafieldscommon}) {
    my($DataFieldName, $DataFieldValue, $Index, @CommonDataFieldsSplit);
    @CommonDataFieldsSplit = split ",", $Options{datafieldscommon};
    if (@CommonDataFieldsSplit % 2) {
	die "Error: Invalid number of comma delimited values, ", scalar(@CommonDataFieldsSplit), ", specified \"",  join(",", @CommonDataFieldsSplit), "\" using \"--datafieldscommon\" option: it must contain even number of values.\n";
    }
    for ($Index = 0; $Index < @CommonDataFieldsSplit; $Index += 2) {
      $DataFieldName = $CommonDataFieldsSplit[$Index];
      $DataFieldValue = $CommonDataFieldsSplit[$Index + 1];
      if (exists $OptionsInfo{SpecifiedCommonFieldMap}{$DataFieldName}) {
	die "Error: Common data field, $DataFieldName, specified more than once using \"--datafieldscommon\" option.\n";
      }
      if (exists($OptionsInfo{SpecifiedNewToOldSDFieldMap}{$DataFieldName}) || exists($OptionsInfo{SpecifiedOldToNewSDFieldMap}{$DataFieldName})) {
	die "Error: Common data field, $DataFieldName, specified using \"--datafieldscommon\" option cannot be specified in \"--datafieldsmap or --datafieldsmapfile\" option.\n";
      }
      $OptionsInfo{SpecifiedCommonFieldMap}{$DataFieldName} = $DataFieldValue;
    }
  }

  $OptionsInfo{DataFieldURL} = $Options{datafieldurl} ? $Options{datafieldurl} : '';
  $OptionsInfo{CreateDataFieldURL} = (exists($Options{datafieldurl}) && length($Options{datafieldurl}) ) ? 1 : 0;

  $OptionsInfo{URLDataFieldLabel} = ""; $OptionsInfo{URLCGIScriptName} = "";
  $OptionsInfo{URLParamName} = ""; $OptionsInfo{URLCmpdIdFieldName} = "";

  if ($OptionsInfo{CreateDataFieldURL}) {
    my(@DataFieldURLSplit, $Value);
    @DataFieldURLSplit = split ",", $Options{datafieldurl};
    if (@DataFieldURLSplit != 4) {
      die "Error: Invalid number of values, ", scalar(@DataFieldURLSplit), ", specified using \"--datafieldURL\" option: it must contain 4 values.\n";
    }
    for $Value (@DataFieldURLSplit) {
      if (!IsNotEmpty($Value)) {
	die "Error: One of the values, $Options{datafieldurl}, specified using \"--datafieldURL\" option is empty.\n";
      }
    }
    $OptionsInfo{URLDataFieldLabel} = $DataFieldURLSplit[0];
    $OptionsInfo{URLCGIScriptName} = $DataFieldURLSplit[1];
    $OptionsInfo{URLParamName}  = $DataFieldURLSplit[2];
    $OptionsInfo{URLCmpdIdFieldName} = $DataFieldURLSplit[3];
  }

}

# Retrieve information about input SD files...
sub RetrieveSDFilesInfo {
  my($Index, $SDFile, $FileDir, $FileName, $FileExt, $OutFileRoot,  $OutFile, $DataFieldName);

  %SDFilesInfo = ();
  @{$SDFilesInfo{FileOkay}} = ();
  @{$SDFilesInfo{OutFile}} = ();

   FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFile = $SDFilesList[$Index];

    $SDFilesInfo{FileOkay}[$Index] = 0;
    $SDFilesInfo{OutFile}[$Index] = '';

    if (!(-e $SDFile)) {
      warn "Warning: Ignoring file $SDFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($SDFile, "sd sdf")) {
      warn "Warning: Ignoring file $SDFile: It's not a SD file\n";
      next FILELIST;
    }
    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($SDFile);
    if ($Options{root} && (@SDFilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($Options{root});
      if ($RootFileName && $RootFileExt) {
	$FileName = $RootFileName;
      }
      else {
	$FileName = $Options{root};
      }
      $OutFileRoot = $FileName;
    }
    else {
      $OutFileRoot = $FileName . "ModifiedDataFields";
    }

    $OutFile = $OutFileRoot . ".$FileExt";
    if (lc($OutFile) eq lc($SDFile)) {
      warn "Warning: Ignoring file $SDFile:Output file name, $OutFile, is same as input SD file name, $SDFile\n";
      next FILELIST;
    }
    if (!$Options{overwrite}) {
      if (-e $OutFile) {
	warn "Warning: Ignoring file $SDFile: The file $OutFile already exists\n";
	next FILELIST;
      }
    }

    $SDFilesInfo{FileOkay}[$Index] = 1;
    $SDFilesInfo{OutFile}[$Index] = $OutFile;
  }
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{detail} = 1;
  $Options{keepolddatafields} = "none";
  $Options{mode} = "molname";
  $Options{molnamemode} = "labelprefix";
  $Options{molnamereplace} = "empty";

  if (!GetOptions(\%Options, "detail|d=i", "datafieldscommon=s", "datafieldsmap=s", "datafieldsmapfile=s", "datafieldurl=s", "help|h", "keepolddatafields|k=s", "mode|m=s", "molname=s", "molnamemode=s", "molnamereplace=s", "overwrite|o", "root|r=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{keepolddatafields} !~ /^(all|unmappedonly|none)$/i) {
    die "Error: The value specified, $Options{keepolddatafields}, for option \"-k --keepolddatafields\" is not valid. Allowed values: all, unmappedonly, or none\n";
  }
  if ($Options{mode} !~ /^(molname|datafields|both)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: molname, datafields, or both\n";
  }
  if ($Options{molnamemode} !~ /^(datafield|labelprefix)$/i) {
    die "Error: The value specified, $Options{molnamemode}, for option \"--molnamemode\" is not valid. Allowed values: datafield or labelprefix\n";
  }
  if ($Options{molnamereplace} !~ /^(always|empty)$/i) {
    die "Error: The value specified, $Options{molnamereplace}, for option \"--molnamereplace\" is not valid. Allowed values: always or empty\n";
  }
  if (!IsPositiveInteger($Options{detail})) {
    die "Error: The value specified, $Options{detail}, for option \"-d --detail\" is not valid. Allowed values: > 0\n";
  }
}

__END__

=head1 NAME

ModifySDFilesDataFields.pl - Modify data fields in SDFile(s)

=head1 SYNOPSIS

ModifySDFilesDataFields.pl SDFile(s)...

ModifySDFilesDataFields.pl [B<-d, --detail> infolevel]
[B<--datafieldscommon> newfieldlabel, newfieldvalue, [newfieldlabel, newfieldvalue,...]]
[B<--datafieldsmap> newfieldlabel, oldfieldlabel, [oldfieldlabel,...]; [newfieldlabel, oldfieldlabel, [oldfieldlabel,...]]]
[B<--datafieldsmapfile> filename] [B<--datafieldURL> URLDataFieldLabel, CGIScriptPath, CGIParamName, CmpdIDFieldLabel]
[B<-h, --help>] [B<-k, --keepolddatafields> all | unmappedonly | none] [B<-m, --mode> molname | datafields | both]
[B<--molnamemode> datafield | labelprefix] [B<--molname> datafieldname or prefixstring]
[B<--molnamereplace> always | empty] [B<-o, --overwrite>] [B<-r, --root> rootname]
[B<-w, --workingdir> dirname] SDFile(s)...

=head1 DESCRIPTION

Modify molname line and data fields in I<SDFile(s)>. Molname line can be replaced by a
data field value or assigned a sequential ID prefixed with a specific string. For data
fields and modification of their values, these types of options are supported: replace
data field labels by another set of labels; combine values of multiple data fields and
assign a new label; add specific set of data field labels and values to all compound
records; and others.

The file names are separated by space.The valid file extensions are I<.sdf> and I<.sd>.
All other file names are ignored. All the SD files in a current directory can be specified
either by I<*.sdf> or the current directory name.

=head1 OPTIONS

=over 4

=item B<-d, --detail> I<infolevel>

Level of information to print about compound records being ignored. Default: I<1>. Possible
values: I<1, 2 or 3>.

=item B<--datafieldscommon> I<newfieldlabel, newfieldvalue, [newfieldlabel, newfieldvalue,...]>

Specify data field labels and values for addition to each compound record. It's a comma delimited
list of data field label and values pair. Default: I<none>.

Examples:

    DepositionDate,YYYY-MM-DD
    Source,www.domainname.org,ReleaseData,YYYY-MM-DD

=item B<--datafieldsmap> I<newfieldlabel, oldfieldlabel, [oldfieldlabel,...]; [newfieldlabel, oldfieldlabel, [oldfieldlabel,...]]>

Specify how various data field labels and values are combined to generate a new data field
labels and their values. All the comma delimited data fields, with in a semicolon delimited set,
are mapped to the first new data field label along with the data field values joined via new
line character. Default: I<none>.

Examples:

    Synonym,Name,SystematicName,Synonym;CmpdID,Extreg
    HBondDonors,SumNHOH

=item B<--datafieldsmapfile> I<filename>

Filename containing mapping of data fields. Format of data fields line in this file corresponds
to B<--datafieldsmap> option. Example:

    Line 1: Synonym,Name,SystematicName,Synonym;CmpdID,Extreg
    Line 2: HBondDonors,SumNHOH


=item B<--datafieldURL> I<URLDataFieldLabel, CGIScriptPath, CGIParamName, CmpdIDFieldLabel>

Specify how to generate a URL for retrieving compound data from a web server and add it
to each compound record. I<URLDataFieldLabel> is used as the data field label for URL value
which is created by combining I<CGIScriptPath,CGIParamName,CmpdIDFieldLabel> values:
CGIScriptPath?CGIParamName=CmpdIDFieldLabelValue. Default: I<none>.

Example:

    Source,http://www.yourdomain.org/GetCmpd.pl,Reg_ID,Mol_ID

=item B<-h, --help>

Print this help message.

=item B<-k, --keepolddatafields> I<all | unmappedonly | none>

Specify how to transfer old data fields from input SDFile(s) to new SDFile(s) during
I<datafields | both> value of B<-m, --mode> option: keep all old data fields; write out the ones
not mapped to new fields as specified by B<--datafieldsmap> or <--datafieldsmapfile> options;
or ignore all old data field labels. For I<molname> B<-m --mode>, old datafields are always kept.
Possible values: I<all | unmappedonly | none>. Default: I<none>.

=item B<-m, --mode> I<molname | datafields | both>

Specify how to modify SDFile(s): I<molname> - change molname line by another datafield or value;
I<datafield> - modify data field labels and values by replacing one label by another, combining
multiple data field labels and values, adding specific set of data field labels and values to all compound, or
inserting an URL for compound retrieval to each record; I<both> - change molname line and datafields
simultaneously. Possible values: I<molname | datafields | both>. Default: I<molname>

=item B<--molnamemode> I<datafield | labelprefix>

Specify how to change molname line for B<-m --mode> option values of I<molname | both>: use
a datafield label value or assign a sequential ID prefixed with I<labelprefix>. Possible values:
I<datafield | labelprefix>. Default: I<labelprefix>.

=item B<--molname> I<datafieldname or prefixstring>

Molname generation method. For I<datafield> value of B<--molnamemode> option, it corresponds
to datafield label name whose value is used for molname; otherwise, it's a prefix string used for
generating compound IDs like labelprefixstring<Number>. Default value, I<Cmpd>, generates
compound IDs like Cmpd<Number> for molname.

=item B<--molnamereplace> I<always | empty>

Specify when to replace molname line for B<-m --mode> option values of I<molname | both>:
always replace the molname line using B<--molname> option or only when it's empty. Possible
values: I<always | empty>. Default: I<empty>.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

New SD file name is generated using the root: <Root>.<Ext>. Default new file
name: <InitialSDFileName>ModifiedDataFields.<Ext>. This option is ignored for multiple
input files.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To replace empty molname lines by Cmpd<CmpdNumber> and generate a new SD file
NewSample1.sdf, type:

    % ModifySDFilesDataFields.pl -o -r NewSample1 Sample1.sdf

To replace all molname lines by Mol_ID data field generate a new SD file
NewSample1.sdf, type:

    % ModifySDFilesDataFields.pl --molnamemode datafield
    --molnamereplace always -r NewSample1 -o Sample1.sdf

To replace all molname lines by Mol_ID data field, map Name and CompoundName to
a new datafield Synonym, and generate a new SD file NewSample1.sdf, type:

    % ModifySDFilesDataFields.pl --molnamemode datafield
      --molnamereplace always --molname Mol_ID --mode both
      --datafieldsmap "Synonym,Name,CompoundName" -r
      NewSample1 -o Sample1.sdf

To replace all molname lines by Mol_ID data field, map Name and CompoundName to
a new datafield Synonym, add common fields ReleaseDate and Source, and
generate a new SD file NewSample1.sdf without keeping any old SD data fields, type:

    % ModifySDFilesDataFields.pl --molnamemode datafield
      --molnamereplace always --molname Mol_ID --mode both
      --datafieldsmap "Synonym,Name,CompoundName"
      --datafieldscommon "ReleaseDate,yyyy-mm-dd,Source,
      www.mayachemtools.org" --keepolddatafields none -r
      NewSample1 -o Sample1.sdf

B<Preparing SD files PubChem deposition:>

Consider a SD file with these fields: Mol_ID, Name, Synonyms and Systematic_Name.
And Mol_ID data field uniquely identifies your compound.

To prepare a new SD file CmpdDataForPubChem.sdf containing only required
PUBCHEM_EXT_DATASOURCE_REGID field, type:

    % ModifySDFilesDataFields.pl --m datafields
      --datafieldsmap
      "PUBCHEM_EXT_DATASOURCE_REGID,Mol_ID"
      -r CmpdDataForPubChem -o Sample1.sdf

To prepare a new SD file CmpdDataForPubChem.sdf containing only required
PUBCHEM_EXT_DATASOURCE_REGID field and replace molname line with Mol_ID, type:

    % ModifySDFilesDataFields.pl --molnamemode datafield
      --molnamereplace always --molname Mol_ID --mode both
      --datafieldsmap
       "PUBCHEM_EXT_DATASOURCE_REGID,Mol_ID"
      -r CmpdDataForPubChem -o Sample1.sdf

In addition to required PubChem data field, you can also add optional PubChem data
fields.

To map your Name, Synonyms and Systematic_Name data fields to optional
PUBCHEM_SUBSTANCE_SYNONYM data field along with required ID field, type:

    % ModifySDFilesDataFields.pl --molnamemode datafield
      --molnamereplace always --molname Mol_ID --mode both
      --datafieldsmap
      "PUBCHEM_EXT_DATASOURCE_REGID,Mol_ID;
      PUBCHEM_SUBSTANCE_SYNONYM,Name,CompoundName"
      -r CmpdDataForPubChem -o Sample1.sdf

To add your <domain.org> as PUBCHEM_EXT_SUBSTANCE_URL and link substance
retrieval to your CGI script <http://www.yourdomain.org/GetCmpd.pl,Reg_ID,Mol_ID>
via PUBCHEM_EXT_DATASOURCE_REGID field along with optional and required
data fields, type:

    % ModifySDFilesDataFields.pl --molnamemode datafield
      --molnamereplace always --molname Mol_ID --mode both
      --datafieldsmap
      "PUBCHEM_EXT_DATASOURCE_REGID,Mol_ID;
      PUBCHEM_SUBSTANCE_SYNONYM,Name,CompoundName"
      --datafieldscommon
      "PUBCHEM_EXT_SUBSTANCE_URL,domain.org"
      --datafieldURL "PUBCHEM_EXT_DATASOURCE_URL,
      http://www.yourdomain.org/GetCmpd.pl,Reg_ID,Mol_ID"
      -r CmpdDataForPubChem -o Sample1.sdf

And to add a publication date and request a release data using
PUBCHEM_PUBLICATION_DATE and PUBCHEM_DEPOSITOR_RECORD_DATE data fields
along with all the data fields in earlier examples, type:
optional fields, type:

    % ModifySDFilesDataFields.pl --molnamemode datafield
      --molnamereplace always --molname Mol_ID --mode both
      --datafieldsmap
      "PUBCHEM_EXT_DATASOURCE_REGID,Mol_ID;
      PUBCHEM_SUBSTANCE_SYNONYM,Name,CompoundName"
      --datafieldURL "PUBCHEM_EXT_DATASOURCE_URL,
      http://www.yourdomain.org/GetCmpd.pl,Reg_ID,Mol_ID"
      --datafieldscommon
      "PUBCHEM_EXT_SUBSTANCE_URL,domain.org,
      PUBCHEM_PUBLICATION_DATE,YYY-MM-DD,
      PUBCHEM_DEPOSITOR_RECORD_DATE,YYYY-MM-DD"
      -r CmpdDataForPubChem -o Sample1.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

InfoSDFiles.pl, JoinSDFiles.pl, MergeTextFilesWithSD.pl, SplitSDFiles.pl, SDFilesToHTML.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
