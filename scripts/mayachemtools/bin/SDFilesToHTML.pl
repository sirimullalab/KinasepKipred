#!/usr/bin/perl -w
#
# File: SDFilesToHTML.pl
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
use File::Spec;
use Text::ParseWords;
use Benchmark;
use Cwd;
use FileUtil;
use SDFileUtil;
use TextUtil;
use HTMLUtil;

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

#Make sure appropriate mode specific option values are specified...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Collect information about SD files...
print "Checking input SD file(s)...\n";
my(%SDFilesInfo);
RetrieveSDFilesInfo();
SetupMultipleTablesAndMiscInfo();

# Generate output files...
my($FileIndex);
if (@SDFilesList > 1) {
  print "\nProcessing SD files...\n";
}
for $FileIndex (0 .. $#SDFilesList) {
  if ($SDFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $SDFilesList[$FileIndex]...\n";
    GenerateHTMLTable($FileIndex);
  }
}
print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Generate HTML table(s)...
sub GenerateHTMLTable {
  my($Index) = @_;

  if ($SDFilesInfo{MultipleHTMLTables}[$Index]) {
    GenerateMultipleHTMLTables($Index);
  }
  else {
    GenerateOneHTMLTable($Index);
  }
}

# Generate one HTML table...
sub GenerateOneHTMLTable {
  my($Index) = @_;
  my($SDFile, $TopHTMLDir, $HTMLFile, $StartCmpdNum, $EndCmpdNum, $CSSFile, $CSSRef, $CSSFilePath, $TableNum);

  $HTMLFile = $SDFilesInfo{HTMLRoot}[$Index] . ".html";
  $SDFile = $SDFilesList[$Index];

  # Setup data directories...
  ($TopHTMLDir) = SetupDataDirs($Index);

  # Setup stylesheet file...
  $CSSRef = "";
  if ($Options{stylesheet} =~ /^new$/i) {
    $CSSFile = $SDFilesInfo{HTMLRoot}[$Index] . ".css"; $CSSRef = ".\/" . "$CSSFile";
    $CSSFilePath = "$TopHTMLDir" . "\/" . $CSSFile;
    GenerateStyleSheetFile($CSSFilePath);
  }
  elsif ($Options{stylesheet} =~ /^old$/i) {
    $CSSRef = $Options{stylesheetname};
  }

  # Set HTML file location...
  $HTMLFile = "$TopHTMLDir" . "\/" . $HTMLFile;

  print "Generating HTML file $HTMLFile...\n";
  open HTMLFILE, ">$HTMLFile" or die "Error: Can't open $HTMLFile: $! \n";
  open SDFILE, "$SDFile" or die "Error: Can't open $SDFile: $! \n";

  # Write out HTML page header...
  print HTMLFILE SetupHTMLPageHeader($SDFilesInfo{HTMLTitle}[$Index], $CSSRef, $OptionsInfo{TopHTMLDirStrViewerJSFileRef});

  if ($OptionsInfo{StrViewerJSFileRef}) {
    print HTMLFILE SetupStrViewerJSInitCmd($OptionsInfo{StrViewerType}, $OptionsInfo{TopHTMLDirStrViewerCodeBase});
  }

  # Setup page title...
  if ($OptionsInfo{TitleDisplay}) {
    print HTMLFILE SetupHTMLPageTitle($SDFilesInfo{HTMLTitle}[$Index]);
  }
  else {
    print HTMLFILE SetupHTMLEmptyLines(1);
  }

  # Start the table...
  print HTMLFILE SetupHTMLAlignmentBegin("center");
  print HTMLFILE SetupHTMLTableHeader($OptionsInfo{TableBorder}, $OptionsInfo{TableCellPadding}, $OptionsInfo{TableCellSpacing});

  # Generate table rows...
  $StartCmpdNum = 1;
  $EndCmpdNum = $SDFilesInfo{CmpdCount}[$Index];
  $TableNum = 1;
  GenerateTableRows($Index, $TableNum, $StartCmpdNum, $EndCmpdNum, \*SDFILE, \*HTMLFILE);

  # Finish up the table...
  print HTMLFILE SetupHTMLTableEnd();
  print HTMLFILE SetupHTMLAlignmentEnd("center");

  # Write out HTML page end...
  print HTMLFILE SetupHTMLPageEnd($OptionsInfo{FooterMsg});

  close HTMLFILE;
  close SDFILE;
}

# Generate multiple tables...
sub GenerateMultipleHTMLTables {
  my($Index) = @_;
  my($TopHTMLDir, $SubHTMLDir, $SDFile, $HTMLFile, $TableNum, $TableCount, $TableIndex, $TableStartCmpdNum, $TableEndCmpdNum, $PrintMsg, $CSSFile, $CSSRef, $CSSFilePath, $NewStyleSheet, $StrViewerCodeBase, $StrViewerJSFileRef);

  # Open SD file...
  $SDFile = $SDFilesList[$Index];
  open SDFILE, "$SDFile" or die "Error: Can't open $SDFile: $! \n";

  # Set up data directories to hold various html files...
  ($TopHTMLDir, $SubHTMLDir) = SetupDataDirs($Index);

  # Create stylesheet file...
  $CSSRef = "";
  $NewStyleSheet = 0;
  if ($Options{stylesheet} =~ /^new$/i) {
    $NewStyleSheet = 1;
    $CSSFile = $SDFilesInfo{HTMLRoot}[$Index] . ".css";
    $CSSFilePath = "$TopHTMLDir" . "\/" . $CSSFile;
    GenerateStyleSheetFile($CSSFilePath);
  }
  elsif ($Options{stylesheet} =~ /^old$/i) {
    $CSSRef = $Options{stylesheetname};
  }

  $PrintMsg = 1;
  # Generate HTML files for all the tables...
  $TableCount = $SDFilesInfo{TableCount}[$Index];
  for $TableNum (1 .. $TableCount) {
    $TableIndex = $TableNum - 1;
    $HTMLFile = ${$SDFilesInfo{TableHTMLFiles}[$Index]}[$TableIndex];
    $TableStartCmpdNum = ${$SDFilesInfo{TableStartCmpdNum}[$Index]}[$TableIndex];
    $TableEndCmpdNum = ${$SDFilesInfo{TableEndCmpdNum}[$Index]}[$TableIndex];

    # Setup file name...
    if ($TableNum == 1) {
      $HTMLFile = "$TopHTMLDir" . "\/" . $HTMLFile;
      print "Generating HTML file $HTMLFile...\n";
    }
    else {
      $HTMLFile = "$SubHTMLDir" . "\/" . $HTMLFile;
      if ($PrintMsg) {
	$PrintMsg = 0;
	if ($TableCount == 2) {
	  print "Generating HTML file $HTMLFile...\n";
	}
	else {
	  print "Generating ", ($TableCount - 1), " other HTML files: $SubHTMLDir\/$SDFilesInfo{HTMLRoot}[$Index]\*.html...\n";
	}
      }
    }
    # Setup stylesheet reference...
    if ($NewStyleSheet) {
      $CSSRef = ($TableNum == 1) ? ".\/" : "..\/";
      $CSSRef .= $CSSFile;
    }

    open HTMLFILE, ">$HTMLFile" or die "Error: Can't open $HTMLFile: $! \n";
    # Write out HTML page header...
    $StrViewerJSFileRef = ($TableNum == 1) ? $OptionsInfo{TopHTMLDirStrViewerJSFileRef} : $OptionsInfo{SubHTMLDirStrViewerJSFileRef};
    print HTMLFILE SetupHTMLPageHeader($SDFilesInfo{HTMLTitle}[$Index], $CSSRef, $StrViewerJSFileRef);

    if ($OptionsInfo{StrViewerJSFileRef}) {
      $StrViewerCodeBase = ($TableNum == 1) ? $OptionsInfo{TopHTMLDirStrViewerCodeBase} : $OptionsInfo{SubHTMLDirStrViewerCodeBase};
      print HTMLFILE SetupStrViewerJSInitCmd($OptionsInfo{StrViewerType}, $StrViewerCodeBase);
    }

    # Set up the navigation links for this table...
    if ($OptionsInfo{NavLinksAtTop}) {
      WriteNavigationLinks($Index, $TableNum, \*HTMLFILE);
    }
    # Setup page title...
    if ($OptionsInfo{TitleDisplay}) {
      print HTMLFILE SetupHTMLPageTitle($SDFilesInfo{HTMLTitle}[$Index]);
    }
    else {
      print HTMLFILE SetupHTMLEmptyLines(1);
    }

    # Start the table...
    print HTMLFILE SetupHTMLAlignmentBegin("center");
    print HTMLFILE SetupHTMLTableHeader($OptionsInfo{TableBorder}, $OptionsInfo{TableCellPadding}, $OptionsInfo{TableCellSpacing});

    # Generate table content...
    GenerateTableRows($Index, $TableNum, $TableStartCmpdNum, $TableEndCmpdNum, \*SDFILE, \*HTMLFILE);

    # Finish up the table...
    print HTMLFILE SetupHTMLTableEnd();
    print HTMLFILE SetupHTMLAlignmentEnd("center");

    # Set up the navigation links for this table...
    if ($OptionsInfo{NavLinksAtBottom}) {
      print HTMLFILE SetupHTMLEmptyLines(1);
      WriteNavigationLinks($Index, $TableNum, \*HTMLFILE);
    }

    # Write out HTML page end...
    print HTMLFILE SetupHTMLPageEnd($OptionsInfo{FooterMsg});
    close HTMLFILE;
  }
  close SDFILE;

}

# Generate table content...
sub GenerateTableRows {
  my($Index, $TableNum, $StartCmpdNum, $EndCmpdNum, $SDFileRef, $HTMLFileRef) = @_;

  if ($OptionsInfo{StructuresOnlyMode}) {
    WriteRowStructures($Index, $TableNum, $StartCmpdNum, $EndCmpdNum, $SDFileRef, $HTMLFileRef);
  }
  else {
    WriteColLabels($Index, $SDFileRef, $HTMLFileRef);
    WriteRowValues($Index, $TableNum, $StartCmpdNum, $EndCmpdNum, $SDFileRef, $HTMLFileRef);
  }
}

# Create stylesheet file...
sub GenerateStyleSheetFile {
  my($CSSFile) = @_;
    print "Generating stylesheet file $CSSFile...\n";
    open CSSFILE, ">$CSSFile" or die "Error: Can't open $CSSFile: $! \n";
    print CSSFILE SetupHTMLStyleSheetTags();
    close CSSFILE;
}

# Write out table header using column labels...
sub WriteColLabels {
  my($Index, $SDFileRef, $HTMLFileRef) = @_;

  my(@ColLabels, $Label);
  print $HTMLFileRef $SDFilesInfo{TableRowHeaderTags};

  # Write out structure label...
  $Label = "Structure";
  print $HTMLFileRef SetupHTMLTableRowHeaderValue($Label);

  # Write out field values..
  @ColLabels = @{$SDFilesInfo{SpecifiedDataFieldLabels}[$Index]};
  for $Label (@ColLabels) {
    print $HTMLFileRef SetupHTMLTableRowHeaderValue($Label);
  }
  print $HTMLFileRef $SDFilesInfo{RowEndTags};
}

# Write out the rows value...
sub WriteRowValues {
  my($Index, $TableNum, $StartCmpdNum, $EndCmpdNum, $SDFileRef, $HTMLFileRef) = @_;
  my($BackgroundColor, $FontColor, $RowNum, $CmpdNum, $CmpdString, @CmpdLines, $Label, %DataFieldValues, $Value);

  $RowNum = 0;
  for $CmpdNum ($StartCmpdNum .. $EndCmpdNum) {
    $RowNum++;
    $CmpdString = ReadCmpdString($SDFileRef);
    if ($OptionsInfo{ShadeRowsStatus}) {
      print $HTMLFileRef ($RowNum % 2) ? $SDFilesInfo{BgFilledOddRowHeaderTags} : $SDFilesInfo{BgFilledEvenRowHeaderTags};
    }
    else {
      print $HTMLFileRef $SDFilesInfo{RowHeaderTags};
    }
    @CmpdLines = split "\n", $CmpdString;
    %DataFieldValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

    # Setup structure column...
    SetupStructureColumn($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef);
    # Write out field values..
    for $Label (@{$SDFilesInfo{SpecifiedDataFieldLabels}[$Index]}) {
      $Value = (IsNotEmpty($DataFieldValues{$Label})) ? $DataFieldValues{$Label} : "";
      $BackgroundColor = ""; $FontColor = "";
      if ($OptionsInfo{HighlightStatus}) {
	if (exists($OptionsInfo{SpecifiedHighlightDataFieldLabelsMap}{$Label})) {
	  ($BackgroundColor, $FontColor) = GetValueHighlightColors($Label, $Value);
	}
      }
      print $HTMLFileRef SetupHTMLTableRowDataValue($Value, $BackgroundColor, $FontColor);
    }
    print $HTMLFileRef $SDFilesInfo{RowEndTags};
  }
}

# Write only structures...
sub WriteRowStructures {
  my($Index, $TableNum, $StartCmpdNum, $EndCmpdNum, $SDFileRef, $HTMLFileRef) = @_;
  my($CmpdNum, $CmpdString, $StartRowFlag, $ColNum, $RowNum, $RowBgColor, $RowStartTags, $ColumnHeaderTags, $ColumnEndTags, $CmpdDataFieldValue, $CmpdHTMLFileRef, $Value);

  $StartRowFlag = 1; $ColNum = 0; $RowNum = 0;
  $ColumnHeaderTags = SetupHTMLTableColumnHeader();
  $ColumnEndTags = SetupHTMLTableColumnEnd();

  if ($OptionsInfo{StructuresOnlyMode} && !$OptionsInfo{TableBorder} && ($OptionsInfo{OddRowsShadeColor} =~ /^(#ffffff|white)$/i)) {
    print $HTMLFileRef SetupHTMLTableRowHeader($OptionsInfo{RowHAlignment}, $OptionsInfo{TableHeaderRowColor}, $OptionsInfo{RowVAlignment});
    $Value = SetupHTMLTableRowDataValue("");
    print $HTMLFileRef InsertHTMLTags($Value, "colspan", "$OptionsInfo{StrTableCols}");
    print $HTMLFileRef $SDFilesInfo{RowEndTags};
  }

  for $CmpdNum ($StartCmpdNum .. $EndCmpdNum) {
    $CmpdString = ReadCmpdString($SDFileRef);
    if ($StartRowFlag) {
      $StartRowFlag = 0;
      $RowNum++;
      if ($OptionsInfo{ShadeRowsStatus}) {
	print $HTMLFileRef ($RowNum % 2) ? $SDFilesInfo{BgFilledOddRowHeaderTags} : $SDFilesInfo{BgFilledEvenRowHeaderTags};
      }
      else {
	print $HTMLFileRef $SDFilesInfo{RowHeaderTags};
      }
    }
    $ColNum++;

    $CmpdDataFieldValue = "";
    if ($OptionsInfo{CmpdDataField}) {
      my($CmpdDataField, @CmpdLines, %DataFieldValues);
      $CmpdDataField = $OptionsInfo{CmpdDataField};
      @CmpdLines = split "\n", $CmpdString;
      %DataFieldValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);
      if (exists($DataFieldValues{$CmpdDataField}) && length($DataFieldValues{$CmpdDataField})) {
	$CmpdDataFieldValue = $DataFieldValues{$CmpdDataField};
	if ($OptionsInfo{CmpdDataFieldLabel} =~ /^yes$/i) {
	  $CmpdDataFieldValue = "${CmpdDataField}: ${CmpdDataFieldValue}";
	}
	# Make sure it's not to looong...
	if (length($CmpdDataFieldValue) > 30) {
	  $CmpdDataFieldValue = substr($CmpdDataFieldValue, 0, 30) . "...";
	}
      }
    }
    if ($CmpdDataFieldValue) {
      $RowBgColor = "";
      if ($OptionsInfo{ShadeRowsStatus}) {
	$RowBgColor = ($RowNum % 2) ? $OptionsInfo{OddRowsShadeColor} : $OptionsInfo{EvenRowsShadeColor};
      }
      $RowStartTags = SetupHTMLTableRowHeader($OptionsInfo{CmpdDataFieldAlignment}, $RowBgColor, "middle");
      # Start  a new table in current column...
      print $HTMLFileRef $ColumnHeaderTags;
      print $HTMLFileRef SetupHTMLAlignmentBegin("center");
      print $HTMLFileRef SetupHTMLTableHeader(0, 0, 0);

      if ($OptionsInfo{CmpdDataFieldPosition} =~ /^top$/i ) {
	# Add an empty row...
	print $HTMLFileRef $RowStartTags;
	print $HTMLFileRef SetupHTMLTableRowDataValue("");
	print $HTMLFileRef $SDFilesInfo{RowEndTags};

	# Display the label value...
	print $HTMLFileRef $RowStartTags;
	$CmpdHTMLFileRef = SetupCompoundSummaryFileAndLink($Index, $TableNum, $CmpdString, $CmpdNum);
	$Value = SetupHTMLHRef("$CmpdDataFieldValue", $CmpdHTMLFileRef, "Compound Summary");
	print $HTMLFileRef SetupHTMLTableRowDataValue($Value);
	print $HTMLFileRef $SDFilesInfo{RowEndTags};
      }
      # Display the structure...
      print $HTMLFileRef SetupHTMLTableRowHeader("center", $RowBgColor, "middle");
      SetupStructureDisplay($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef);
      print $HTMLFileRef $SDFilesInfo{RowEndTags};

      if ($OptionsInfo{CmpdDataFieldPosition} =~ /^bottom$/i ) {
	# Display the label value...
	print $HTMLFileRef $RowStartTags;
	$CmpdHTMLFileRef = SetupCompoundSummaryFileAndLink($Index, $TableNum, $CmpdString, $CmpdNum);
	$Value = SetupHTMLHRef("$CmpdDataFieldValue", $CmpdHTMLFileRef, "Compound Summary");
	print $HTMLFileRef SetupHTMLTableRowDataValue($Value);
	print $HTMLFileRef $SDFilesInfo{RowEndTags};

	# Add an empty row...
	print $HTMLFileRef $RowStartTags;
	print $HTMLFileRef SetupHTMLTableRowDataValue("");
	print $HTMLFileRef $SDFilesInfo{RowEndTags};
      }

      print $HTMLFileRef SetupHTMLTableEnd();
      print $HTMLFileRef SetupHTMLAlignmentEnd("center");
      print $HTMLFileRef $ColumnEndTags;
    }
    else {
      SetupStructureDisplay($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef);
    }

    if ($ColNum == $OptionsInfo{StrTableCols}) {
      # Finish up the current row and get ready for starting a new row...
      print $HTMLFileRef $SDFilesInfo{RowEndTags};
      $ColNum = 0;
      $StartRowFlag = 1;
    }
  }
  if (!$StartRowFlag) {
    # Finish up an existing row...
    my($ColIndex, $Value);
    $Value = "";
    for $ColIndex ($ColNum .. ($OptionsInfo{StrTableCols} - 1) ) {
      print $HTMLFileRef SetupHTMLTableRowDataValue($Value);
    }
    print $HTMLFileRef $SDFilesInfo{RowEndTags};
  }
}

# Setup structure column...
sub SetupStructureColumn {
  my($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef) = @_;

  if ($OptionsInfo{DisplayStructure}) {
    SetupStructureDisplay($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef);
  }
  else {
    SetupStructureLink($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef);
  }
}

# Setup structure display for compound summary page...
sub SetupStructureDisplayForCmpdSummaryPage {
  my($Index, $CmpdString, $CmpdNum, $HTMLFileRef) = @_;
  my($TableNum, $RowNum);

  # Use table num 0 to force usage of "../mol" prefix for all MOL file references. Row num
  # doesn't matter...
  $TableNum = 0;
  $RowNum = 1;

  $OptionsInfo{SettingUpCmpdSummaryPage} = 1;

  # Setup size and bgcolor parameters for linking structures...
  $OptionsInfo{StrViewerParams}{width} = $OptionsInfo{StrLinkWidth};
  $OptionsInfo{StrViewerParams}{height} = $OptionsInfo{StrLinkHeight};
  $OptionsInfo{StrViewerParams}{bgcolor} = $OptionsInfo{StrLinkBgColorSpecified};

  SetupStructureDisplay($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef);

  # Reset size and bgcolor parameters back to displaying structures in tables...
  $OptionsInfo{StrViewerParams}{width} = $OptionsInfo{StrWidth};
  $OptionsInfo{StrViewerParams}{height} = $OptionsInfo{StrHeight};
  $OptionsInfo{StrViewerParams}{bgcolor} = $OptionsInfo{StrBgColorSpecified} ? $OptionsInfo{StrBgColorSpecified} : "";

  $OptionsInfo{SettingUpCmpdSummaryPage} = 0;
}


# Setup structure column display...
sub SetupStructureDisplay {
  my($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef) = @_;
  my($Nothing);

 STRVIEWERTYPE: {
    if ($OptionsInfo{StrViewerType} =~ /^JME$/i) { SetupJMEDisplay($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef); last STRVIEWERTYPE; }
    if ($OptionsInfo{StrViewerType} =~ /^Jmol$/i) { SetupJmolDisplay($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef); last STRVIEWERTYPE; }
    if ($OptionsInfo{StrViewerType} =~ /^Chime$/i) { SetupChimeDisplay($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef); last STRVIEWERTYPE; }
    if ($OptionsInfo{StrViewerType} =~ /^(Chem3DActiveX|ChemDrawActiveX|ChemDrawPlugIn)$/i) { SetupCambridgeSoftDisplay($OptionsInfo{StrViewerType}, $Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef); last STRVIEWERTYPE; }
    if ($OptionsInfo{StrViewerType} =~ /^MarvinView$/i) { SetupMarvinDisplay($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef); last STRVIEWERTYPE; }
    if ($OptionsInfo{StrViewerType} =~ /^ViewerActiveX$/i) { SetupViewerAccelrysActiveXDisplay($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef); last STRVIEWERTYPE; }
    $Nothing = 1;
  }
}

# Setup JME display...
sub SetupJMEDisplay {
  my($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef) = @_;
  my($MolString, $AppletBGColor, $Value, $ValueTag, $AppletName, $StrViewerCodeBase);

  $Value = "";
  ($MolString) = split "$SDFilesInfo{MolEndTag}", $CmpdString;
  if (IsNotEmpty($MolString)) {
    $AppletBGColor = SetupStructureBGColor($RowNum);
    $MolString .= "$SDFilesInfo{MolEndTag}";

    # JME viewer doesn't appear to support "bgcolor" param. So, always use white background for
    # structure cell...
    $AppletName = "JME" . $CmpdNum;
    $OptionsInfo{StrViewerParams}{name} = $AppletName;
    if (!$OptionsInfo{SettingUpCmpdSummaryPage}) {
      if (!$OptionsInfo{StrBgColorSpecified}) {
	$OptionsInfo{StrViewerParams}{bgcolor} = $AppletBGColor;
      }
    }
    $StrViewerCodeBase = ($TableNum == 1) ? $OptionsInfo{TopHTMLDirStrViewerCodeBase} : $OptionsInfo{SubHTMLDirStrViewerCodeBase};
    $Value = SetupStrViewerJMEApplet($MolString, $StrViewerCodeBase, \%{$OptionsInfo{StrViewerParams}});
    $ValueTag = SetupHTMLTableRowDataValue($Value, $SDFilesInfo{White});
  }
  else {
    $ValueTag = SetupHTMLTableRowDataValue($Value);
  }
  if ($OptionsInfo{SettingUpCmpdSummaryPage}) {
    $ValueTag = InsertHTMLTags($ValueTag, ("class", "box"));
  }
  print $HTMLFileRef $ValueTag;
}

# Setup Marvin display...
sub SetupMarvinDisplay {
  my($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef) = @_;
  my($MolString, $AppletBGColor, $Value, $ValueTag, $AppletName, $StrViewerCodeBase);

  $Value = "";
  ($MolString) = split "$SDFilesInfo{MolEndTag}", $CmpdString;
  if (IsNotEmpty($MolString)) {
    $AppletBGColor = SetupStructureBGColor($RowNum);
    $MolString .= "$SDFilesInfo{MolEndTag}";

    $AppletName = "MView" . $CmpdNum;
    $OptionsInfo{StrViewerParams}{name} = $AppletName;
    if (!$OptionsInfo{SettingUpCmpdSummaryPage}) {
      if (!$OptionsInfo{StrBgColorSpecified}) {
	$OptionsInfo{StrViewerParams}{bgcolor} = $AppletBGColor;
      }
    }
    $StrViewerCodeBase = ($TableNum == 1) ? $OptionsInfo{TopHTMLDirStrViewerCodeBase} : $OptionsInfo{SubHTMLDirStrViewerCodeBase};
    $Value = SetupStrViewerMarvinViewApplet($MolString, $StrViewerCodeBase, \%{$OptionsInfo{StrViewerParams}});
    $ValueTag = SetupHTMLTableRowDataValue($Value);
  }
  else {
    $ValueTag = SetupHTMLTableRowDataValue($Value);
  }
  if ($OptionsInfo{SettingUpCmpdSummaryPage}) {
    $ValueTag = InsertHTMLTags($ValueTag, ("class", "box"));
  }
  print $HTMLFileRef $ValueTag;
}

# Setup Jmol display...
sub SetupJmolDisplay {
  my($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef) = @_;
  my($MolString, $AppletBGColor, $Value, $ValueTag, $AppletName, $StrViewerCodeBase);

  $Value = ""; $ValueTag = "";
  ($MolString) = split "$SDFilesInfo{MolEndTag}", $CmpdString;
  if (IsNotEmpty($MolString)) {
    $AppletBGColor = SetupStructureBGColor($RowNum);
    $MolString .= "$SDFilesInfo{MolEndTag}";

    # Make sure MolName line is not empty; otherwise, JMol doesn't display structure...
    my(@MolLines) = split "\n", $MolString;
    if (IsEmpty($MolLines[0])) {
      $MolLines[0] = "Cmpd${CmpdNum}";
      $MolString = join "\n", @MolLines;
    }

    # Setup the applet...
    $AppletName = "Jmol" . $CmpdNum;
    $OptionsInfo{StrViewerParams}{name} = $AppletName;
    if (!$OptionsInfo{SettingUpCmpdSummaryPage}) {
      if (!$OptionsInfo{StrBgColorSpecified}) {
	$OptionsInfo{StrViewerParams}{bgcolor} = $AppletBGColor;
      }
    }
    $StrViewerCodeBase = ($TableNum == 1) ? $OptionsInfo{TopHTMLDirStrViewerCodeBase} : $OptionsInfo{SubHTMLDirStrViewerCodeBase};
    $Value = SetupStrViewerJmolApplet($MolString, $StrViewerCodeBase, \%{$OptionsInfo{StrViewerParams}});
    $ValueTag = SetupHTMLTableRowDataValue($Value);
  }
  else {
    $ValueTag = SetupHTMLTableRowDataValue($Value);
  }
  if ($OptionsInfo{SettingUpCmpdSummaryPage}) {
    $ValueTag = InsertHTMLTags($ValueTag, ("class", "box"));
  }
  print $HTMLFileRef $ValueTag;
}

# Setup Chime display...
sub SetupChimeDisplay {
  my($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef) = @_;
  my($MolString, $BGColor, $Value, $ValueTag, $MolFileRef);

  $Value = "";
  ($MolString) = split "$SDFilesInfo{MolEndTag}", $CmpdString;
  if (IsNotEmpty($MolString)) {
    $BGColor = SetupStructureBGColor($RowNum);
    $MolString .= "$SDFilesInfo{MolEndTag}";
    # Write out MOL file...
    $MolFileRef = SetupMOLFile($Index, $TableNum, $MolString, $CmpdNum);
    # Setup the applet...
    if (!$OptionsInfo{SettingUpCmpdSummaryPage}) {
      if (!$OptionsInfo{StrBgColorSpecified}) {
	$OptionsInfo{StrViewerParams}{bgcolor} = $BGColor;
      }
    }
    $Value = SetupStrViewerChimePlugIn($MolFileRef, \%{$OptionsInfo{StrViewerParams}});
    $ValueTag = SetupHTMLTableRowDataValue($Value);
  }
  else {
    $ValueTag = SetupHTMLTableRowDataValue($Value);
  }
  if ($OptionsInfo{SettingUpCmpdSummaryPage}) {
    $ValueTag = InsertHTMLTags($ValueTag, ("class", "box"));
  }
  print $HTMLFileRef $ValueTag;
}

# Setup displays for various viewers available from CambridgeSoft...
sub SetupCambridgeSoftDisplay {
  my($ViewerType, $Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef) = @_;
  my($MolString, $BGColor, $Value, $ValueTag, $MolFileRef, $Name);

  $Value = "";
  ($MolString) = split "$SDFilesInfo{MolEndTag}", $CmpdString;
  if (IsNotEmpty($MolString)) {
    $BGColor = SetupStructureBGColor($RowNum);
    $MolString .= "$SDFilesInfo{MolEndTag}";
    # Write out MOL file...
    $MolFileRef = SetupMOLFile($Index, $TableNum, $MolString, $CmpdNum);
    # Setup the viewer...
    $Name = "CS" . $CmpdNum;
    if ($ViewerType =~ /^Chem3DActiveX$/) {
      # Use white background for Chem3D and cell; otherwise, it doesn't look good:
      # cell size is larger than Chem3D window size and different colors don't work
      $BGColor = $SDFilesInfo{White};
      $OptionsInfo{StrViewerParams}{name} = $Name;
      if (!$OptionsInfo{SettingUpCmpdSummaryPage}) {
	if (!$OptionsInfo{StrBgColorSpecified}) {
	  $OptionsInfo{StrViewerParams}{bgcolor} = $BGColor;
	}
      }
      $Value = SetupStrViewerChem3DActiveX($MolFileRef, \%{$OptionsInfo{StrViewerParams}});
      $ValueTag = SetupHTMLTableRowDataValue($Value, $BGColor);
    }
    elsif ($ViewerType =~ /^ChemDrawActiveX$/i) {
      # BGColor is not supported. So, make it all white...
      $BGColor = $SDFilesInfo{White};
      $OptionsInfo{StrViewerParams}{name} = $Name;
      if (!$OptionsInfo{SettingUpCmpdSummaryPage}) {
	if (!$OptionsInfo{StrBgColorSpecified}) {
	  $OptionsInfo{StrViewerParams}{bgcolor} = $BGColor;
	}
      }
      $Value = SetupStrViewerChemDrawActiveX($MolFileRef, \%{$OptionsInfo{StrViewerParams}});
      $ValueTag = SetupHTMLTableRowDataValue($Value, $BGColor);
    }
    elsif ($ViewerType =~ /^ChemDrawPlugIn$/i) {
      # BGColor is not supported. So, make it all white...
      $BGColor = $SDFilesInfo{White};
      if (!$OptionsInfo{SettingUpCmpdSummaryPage}) {
	if (!$OptionsInfo{StrBgColorSpecified}) {
	  $OptionsInfo{StrViewerParams}{bgcolor} = $BGColor;
	}
      }
      $Value = SetupStrViewerChemDrawPlugIn($MolFileRef, \%{$OptionsInfo{StrViewerParams}});
      $ValueTag = SetupHTMLTableRowDataValue($Value, $BGColor);
    }
  }
  else {
    $ValueTag = SetupHTMLTableRowDataValue($Value);
  }
  if ($OptionsInfo{SettingUpCmpdSummaryPage}) {
    $ValueTag = InsertHTMLTags($ValueTag, ("class", "box"));
  }
  print $HTMLFileRef $ValueTag;
}

# Setup Accelrys Viewer ActiveX display...
sub SetupViewerAccelrysActiveXDisplay {
  my($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef) = @_;
  my($MolString, $BGColor, $Value, $ValueTag, $Name, $MolFileRef);

  $Value = "";
  ($MolString) = split "$SDFilesInfo{MolEndTag}", $CmpdString;
  if (IsNotEmpty($MolString)) {
    $BGColor = SetupStructureBGColor($RowNum);
    $MolString .= "$SDFilesInfo{MolEndTag}";
    # Write out MOL file. Accelrys ActiveX viewer doesn't load mol files with relative path names.
    # So, set up a complete path name for now; however, it may lead to issues during web
    # deployment.
    my($CompletePath) = 1;
    $MolFileRef = SetupMOLFile($Index, $TableNum, $MolString, $CmpdNum, $CompletePath);
    # Setup the viewer...
    $Name = "ViewerActiveX" . $CmpdNum;
    if (!$OptionsInfo{SettingUpCmpdSummaryPage}) {
      if (!$OptionsInfo{StrBgColorSpecified}) {
	$OptionsInfo{StrViewerParams}{bgcolor} = $BGColor;
      }
    }
    $OptionsInfo{StrViewerParams}{name} = $Name;
    $Value = SetupStrViewerAccelrysActiveX($MolFileRef, \%{$OptionsInfo{StrViewerParams}});
    $ValueTag = SetupHTMLTableRowDataValue($Value);
  }
  else {
    $ValueTag = SetupHTMLTableRowDataValue($Value);
  }
  if ($OptionsInfo{SettingUpCmpdSummaryPage}) {
    $ValueTag = InsertHTMLTags($ValueTag, ("class", "box"));
  }
  print $HTMLFileRef $ValueTag;
}


# Setup structure background color...
sub SetupStructureBGColor {
  my($RowNum) = @_;
  my($BGColor);

  $BGColor = "";
  if ($OptionsInfo{ShadeRowsStatus}) {
    $BGColor =  ($RowNum % 2) ? $OptionsInfo{OddRowsShadeColor} : $OptionsInfo{EvenRowsShadeColor};
  }
  else {
    $BGColor = $SDFilesInfo{White};
  }
  return $BGColor;
}

# Setup  MDL MOL file...
sub SetupMOLFile {
  my($Index, $TableNum, $MolString, $CmpdNum, $CompletePath);
  my($SubMolDir, $MolFileName, $MolFile, $MolFileRef);

  $CompletePath = "";
  if (@_ == 5) {
    ($Index, $TableNum, $MolString, $CmpdNum, $CompletePath) = @_;
  }
  else {
    ($Index, $TableNum, $MolString, $CmpdNum) = @_;
  }

  $SubMolDir = $SDFilesInfo{SubMolDir}[$Index];
  $MolFileName = $SDFilesInfo{HTMLRoot}[$Index] . "Cmpd" . $CmpdNum . ".mol";
  $MolFile = $SubMolDir . "\/" . $MolFileName;

  open MOLFILE, ">$MolFile" or die "Error: Can't open $MolFile: $! \n";
  print MOLFILE "$MolString\n";
  close MOLFILE;

  if ($CompletePath) {
    my($CWD, $NewCWD);
    $CWD = cwd();
    $NewCWD = ConvertCygwinPath($CWD);
    $MolFileRef = $NewCWD . "\/" . $SDFilesInfo{TopHTMLDir}[$Index] .  "\/mol\/$MolFileName" ;
  }
  else {
    $MolFileRef = ($TableNum == 1) ? ".\/mol\/$MolFileName" : "..\/mol\/$MolFileName";
  }

  return $MolFileRef;
}

# Setup a link to structure and other available information...
sub SetupStructureLink {
  my($Index, $TableNum, $RowNum, $CmpdString, $CmpdNum, $HTMLFileRef) = @_;
  my($CmpdHTMLFileRef, $Value);

  $CmpdHTMLFileRef = SetupCompoundSummaryFileAndLink($Index, $TableNum, $CmpdString, $CmpdNum);

  if ($Options{strlinktype} =~ /^button$/i) {
    $Value = SetupHTMLButtonRef("View", $CmpdHTMLFileRef);
  }
  else {
    $Value = SetupHTMLHRef("View", $CmpdHTMLFileRef);
  }
  print $HTMLFileRef SetupHTMLTableRowDataValue($Value);
}

# Setup HTML compound summary file and link...
sub SetupCompoundSummaryFileAndLink {
  my($Index, $TableNum, $CmpdString, $CmpdNum) = @_;
  my($CmpdHTMLFile, $CmpdHTMLFileName, $CmpdHTMLFileRef, $CSSRef, @CmpdLines, $Label, @DataFieldLabels, %DataFieldValues, $Value, $Tag);

  # Setup compound file names...
  $CmpdHTMLFileName = $SDFilesInfo{HTMLRoot}[$Index] . "Cmpd" . $CmpdNum . ".html";
  $CmpdHTMLFile = $SDFilesInfo{SubHTMLDir}[$Index] . "\/" . $CmpdHTMLFileName;

  # Setup stylesheet reference....
  $CSSRef = "";
  if ($Options{stylesheet} =~ /^old$/i) {
    $CSSRef = $Options{stylesheetname};
  }
  else {
    $CSSRef = "..\/" . $SDFilesInfo{HTMLRoot}[$Index] . ".css";
  }

  # Write out compound data in a new HTML file. For summary page, usage of even and odd row shade color
  # is reversed: it causes structure background to be white by default...
  open CMPDHTMLFILE, ">$CmpdHTMLFile" or die "Error: Can't open $CmpdHTMLFile: $! \n";
  print CMPDHTMLFILE SetupHTMLPageHeader($OptionsInfo{StrLinkTitle}, $CSSRef, $OptionsInfo{SubHTMLDirStrViewerJSFileRef});

  if ($OptionsInfo{StrViewerJSFileRef}) {
    print CMPDHTMLFILE SetupStrViewerJSInitCmd($OptionsInfo{StrViewerType}, $OptionsInfo{SubHTMLDirStrViewerCodeBase});
  }

  if ($OptionsInfo{StrLinkTitleDisplay}) {
    print CMPDHTMLFILE SetupHTMLPageTitle($OptionsInfo{StrLinkTitle}, "center");
  }
  else {
    print CMPDHTMLFILE SetupHTMLEmptyLines(1);
  }
  print CMPDHTMLFILE SetupHTMLAlignmentBegin("center");

  # Setup structure display ...
  print CMPDHTMLFILE SetupHTMLTableHeader(0, 5, 2);

  print CMPDHTMLFILE SetupHTMLTableRowHeader("center", "#ffffff", "middle");

  SetupStructureDisplayForCmpdSummaryPage($Index, $CmpdString, $CmpdNum, \*CMPDHTMLFILE);
  print CMPDHTMLFILE $SDFilesInfo{RowEndTags};

  if ($Options{strlinkmode} =~ /^plain$/i) {
    print CMPDHTMLFILE SetupHTMLTableRowHeader("center", $OptionsInfo{StrLinkShadeColor});
    $Tag = SetupHTMLTableRowDataValue("");
    print CMPDHTMLFILE $Tag;
    print CMPDHTMLFILE $SDFilesInfo{RowEndTags};
  }

  print CMPDHTMLFILE SetupHTMLTableRowHeader("left", "", "middle");
  # Start a new table with two columns, one each for data field labels and values, in
  # current column...
  print CMPDHTMLFILE SetupHTMLTableColumnHeader();
  print CMPDHTMLFILE SetupHTMLAlignmentBegin("left");
  print CMPDHTMLFILE SetupHTMLTableHeader(0, 5, 2);

  # Setup table for other available data...
  my($CmpdRowHeaderTags);
  $CmpdRowHeaderTags = SetupHTMLTableRowHeader("left", "", "middle");

  @CmpdLines = split "\n", $CmpdString;

  @DataFieldLabels = GetCmpdDataHeaderLabels(\@CmpdLines);
  %DataFieldValues = GetCmpdDataHeaderLabelsAndValues(\@CmpdLines);

  my($LabelWrapLength, $ValueWrapLength, $LabelColWidth);
  $LabelWrapLength = 30; $ValueWrapLength = 60; $LabelColWidth = 40;

  for $Label (@DataFieldLabels) {
    $Value =  $DataFieldValues{$Label};
    $Label .= ":";
    if ($Label && (length($Label) > $LabelWrapLength)) {
      $Label = WrapText($Label,  $LabelWrapLength, "<br>");
    }
    print CMPDHTMLFILE $CmpdRowHeaderTags;
    if ($Options{strlinkmode} =~ /^plain$/i) {
      $Tag = SetupHTMLTableRowDataValue($Label, "", "", 1);
    }
    else {
      $Tag = SetupHTMLTableRowDataValue($Label, $OptionsInfo{StrLinkShadeColor});
    }
    $Tag = InsertHTMLTags($Tag, "width", "$LabelColWidth");
    print CMPDHTMLFILE $Tag;

    if ($Value && (length($Value) >=$ValueWrapLength) && $Value !~ /a href/i) {
      $Value =~ s/(\r\n)|(\r)|\n//g;
      $Value = WrapText($Value,  $ValueWrapLength, "<br>");
    }
    $Tag = SetupHTMLTableRowDataValue($Value);
    print CMPDHTMLFILE $Tag;
    print CMPDHTMLFILE $SDFilesInfo{RowEndTags};
  }

  # Finish up table holding numerical data...
  print CMPDHTMLFILE SetupHTMLTableEnd();
  print CMPDHTMLFILE SetupHTMLAlignmentEnd("left");
  print CMPDHTMLFILE SetupHTMLTableColumnEnd();
  print CMPDHTMLFILE $SDFilesInfo{RowEndTags};

  # Finish up main table...
  print CMPDHTMLFILE SetupHTMLTableEnd();
  print CMPDHTMLFILE SetupHTMLAlignmentEnd("center");

  if ($OptionsInfo{StrLinkNavigation} && ($SDFilesInfo{CmpdCount}[$Index] > 1) ) {
    print CMPDHTMLFILE SetupHTMLEmptyLines(1);
    WriteCompoundSummaryNavigationLinks($Index, $TableNum, $CmpdNum, \*CMPDHTMLFILE);
  }

  print CMPDHTMLFILE SetupHTMLPageEnd($OptionsInfo{FooterMsg});
  close CMPDHTMLFILE;

  # Add a link to compound's HTML file in table cell...
  $CmpdHTMLFileRef = ($TableNum == 1) ? ".\/html\/" : ".\/";
  $CmpdHTMLFileRef .= $CmpdHTMLFileName;

  return $CmpdHTMLFileRef;
}

# Write navigation link information for compound summary page...
sub WriteCompoundSummaryNavigationLinks {
  my($Index, $CurTableNum, $CurCmpdNum, $CmpdHTMLFileRef) = @_;
  my($FirstTableNum, $CurTableIndex, $FirstCmpdNum, $LastCmpdNum, $PreviousCmpdNum, $NextCmpdNum, $HTMLFile, $HTMLRefFile, $HTMLRefValue);

  $FirstTableNum = 1;
  $FirstCmpdNum = 1;

  $CurTableIndex = $CurTableNum - 1;

  if ($SDFilesInfo{MultipleHTMLTables}[$Index]) {
    my($FirstTableIndex, $LastTableNum, $LastTableIndex);
    $FirstTableIndex = $FirstTableNum - 1;
    $LastTableNum = $SDFilesInfo{TableCount}[$Index]; $LastTableIndex = $LastTableNum - 1;
    $LastCmpdNum = ${$SDFilesInfo{TableEndCmpdNum}[$Index]}[$LastTableIndex];
  }
  else {
    $LastCmpdNum = $SDFilesInfo{CmpdCount}[$Index];
  }

  $PreviousCmpdNum = ($CurCmpdNum == $FirstCmpdNum) ? 0 : ($CurCmpdNum - 1);
  $NextCmpdNum = ($CurCmpdNum == $LastCmpdNum) ? 0 : ($CurCmpdNum + 1);

  my($InactiveLinkNumColor, $InactiveLinkFontBold) = ("#8e2323", "1");
  my($LinkTextColor, $LinkBGColor, $LinkFontBold) = ("", "", "0");
  my($BorderWidth, $CellPadding, $CellSpacing) = (0, 2, 2);

  # Start link table...
  print $CmpdHTMLFileRef SetupHTMLAlignmentBegin("center");
  print $CmpdHTMLFileRef SetupHTMLDivBegin("tablenav");
  print $CmpdHTMLFileRef  SetupHTMLTableHeader($BorderWidth, $CellPadding, $CellSpacing);
  print $CmpdHTMLFileRef $SDFilesInfo{RowHeaderTags};

  print $CmpdHTMLFileRef SetupHTMLTableRowDataValue("Compounds: ");

  # Setup a link to first compound...
  if ($CurCmpdNum != $FirstCmpdNum) {
    $HTMLFile = $SDFilesInfo{HTMLRoot}[$Index] . "Cmpd" . $FirstCmpdNum . ".html";
    $HTMLRefFile = "./${HTMLFile}";
    $HTMLRefValue = SetupHTMLHRef("First", $HTMLRefFile, "First Compound");
    print $CmpdHTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue, $LinkBGColor, $LinkTextColor, $LinkFontBold);
  }

  # Setup a link to previous compund
  if ($PreviousCmpdNum) {
    $HTMLFile = $SDFilesInfo{HTMLRoot}[$Index] . "Cmpd" . $PreviousCmpdNum . ".html";
    $HTMLRefFile = "./${HTMLFile}";
    $HTMLRefValue = SetupHTMLHRef("Previous", $HTMLRefFile, "Previous Compound");
    print $CmpdHTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue, $LinkBGColor, $LinkTextColor, $LinkFontBold);
  }

  # Setup a link to compound table...
  if ($SDFilesInfo{MultipleHTMLTables}[$Index]) {
    $HTMLFile = ${$SDFilesInfo{TableHTMLFiles}[$Index]}[$CurTableIndex];
  }
  else {
    $HTMLFile = $SDFilesInfo{HTMLRoot}[$Index] . ".html";
  }
  $HTMLRefFile = (($CurTableNum == $FirstTableNum) ? "../" : "./") . $HTMLFile;
  $HTMLRefValue = SetupHTMLHRef("Table", $HTMLRefFile, "Table With This Compound");
  print $CmpdHTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue, $LinkBGColor, $LinkTextColor, $LinkFontBold);

  # Setup a link to next compound...
  if ($NextCmpdNum) {
    $HTMLFile = $SDFilesInfo{HTMLRoot}[$Index] . "Cmpd" . $NextCmpdNum . ".html";
    $HTMLRefFile = "./${HTMLFile}";
    $HTMLRefValue = SetupHTMLHRef("Next", $HTMLRefFile, "Next Compound");
    print $CmpdHTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue, $LinkBGColor, $LinkTextColor, $LinkFontBold);
  }

  # Setup a link to last compund
  if ($CurCmpdNum != $LastCmpdNum) {
    $HTMLFile = $SDFilesInfo{HTMLRoot}[$Index] . "Cmpd" . $LastCmpdNum . ".html";
    $HTMLRefFile = "./${HTMLFile}";
    $HTMLRefValue = SetupHTMLHRef("Last", $HTMLRefFile, "Last Compound");
    print $CmpdHTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue, $LinkBGColor, $LinkTextColor, $LinkFontBold);
  }

  # Setup current table info text....
  print $CmpdHTMLFileRef SetupHTMLTableRowDataValue("&nbsp");
  print $CmpdHTMLFileRef SetupHTMLTableRowDataValue("&nbsp");
  print $CmpdHTMLFileRef SetupHTMLTableRowDataValue("Showing $CurCmpdNum of $LastCmpdNum");

  print $CmpdHTMLFileRef $SDFilesInfo{RowEndTags};

  # End link table...
  print $CmpdHTMLFileRef SetupHTMLTableEnd();
  print $CmpdHTMLFileRef SetupHTMLDivEnd();
  print $CmpdHTMLFileRef SetupHTMLAlignmentEnd("center");
}

# Setup navigation link information for each table.
#
# All table sets besides first and last have these links: FirstTable, Previous, Current-1,Current,Current+1,  Next, and LastTable
# First set: Current, Next, and LastTable
# Last set: FirstTable, Previous and Current.
#
sub WriteNavigationLinks {
  my($Index, $CurTableNum, $HTMLFileRef) = @_;
  my($TableNum, $StartTableNum, $EndTableNum, $TableIndex, $BorderWidth, $CellPadding, $CellSpacing,$HTMLFile, $HTMLRefFile, $RelativeFileDir, $HTMLRefValue, $FirstTableNum, $FirstTableIndex, $LastTableNum, $LastTableIndex, $TableStartCmpdNum, $TableEndCmpdNum, $LastCmpdNum, $BGColor, $LinksOffSet);

  $LinksOffSet = 10;

  $FirstTableNum = 1; $FirstTableIndex = $FirstTableNum - 1;
  $LastTableNum = $SDFilesInfo{TableCount}[$Index]; $LastTableIndex = $LastTableNum - 1;
  $LastCmpdNum = ${$SDFilesInfo{TableEndCmpdNum}[$Index]}[$LastTableIndex];

  # Figure out which links to display for a particular table...
  $StartTableNum = $CurTableNum - $LinksOffSet + 1;
  $StartTableNum = ($StartTableNum < $FirstTableNum) ? $FirstTableNum : $StartTableNum;
  if ($CurTableNum < $LinksOffSet) {
    $EndTableNum = $LinksOffSet;
  }
  else {
    $EndTableNum = $CurTableNum + $LinksOffSet - 1;
  }
  $EndTableNum = ($EndTableNum > $LastTableNum) ? $LastTableNum : $EndTableNum;

  my($InactiveLinkNumColor, $InactiveLinkFontBold) = ("#8e2323", "1");
  my($LinkTextColor, $LinkBGColor, $LinkFontBold) = ("", "", "1");

  # Start link table...
  $BorderWidth = 0; $CellPadding = 2; $CellSpacing = 2;
  print $HTMLFileRef SetupHTMLAlignmentBegin("center");
  print $HTMLFileRef SetupHTMLDivBegin("tablenav");
  print $HTMLFileRef  SetupHTMLTableHeader($BorderWidth, $CellPadding, $CellSpacing);
  print $HTMLFileRef $SDFilesInfo{RowHeaderTags};

  if ($OptionsInfo{NavLinksTableInfo} && $OptionsInfo{NavLinksCmpdInfo}) {
    print $HTMLFileRef SetupHTMLTableRowDataValue("Showing table $CurTableNum of $LastTableNum");
    print $HTMLFileRef SetupHTMLTableRowDataValue("&nbsp");
    print $HTMLFileRef SetupHTMLTableRowDataValue("&nbsp");
  }

  print $HTMLFileRef SetupHTMLTableRowDataValue("Tables: ");
  # Setup a link to first table...
  if ($StartTableNum != $FirstTableNum) {
    $HTMLFile = ${$SDFilesInfo{TableHTMLFiles}[$Index]}[$FirstTableIndex];
    $HTMLRefFile = GetRelativeFileDir($CurTableNum, $FirstTableNum, $FirstTableNum) . $HTMLFile;
    $TableStartCmpdNum = ${$SDFilesInfo{TableStartCmpdNum}[$Index]}[$FirstTableIndex];
    $TableEndCmpdNum = ${$SDFilesInfo{TableEndCmpdNum}[$Index]}[$FirstTableIndex];
    $HTMLRefValue = SetupHTMLHRef("First", $HTMLRefFile, "First Table Containing Compounds $TableStartCmpdNum To $TableEndCmpdNum");
    print $HTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue, $LinkBGColor, $LinkTextColor, $LinkFontBold);
  }

  # Setup link to previous table...
  if ($CurTableNum != $FirstTableNum) {
    my($PreviousTableNum, $PreviousTableIndex);
    $PreviousTableNum = $CurTableNum - 1; $PreviousTableIndex = $PreviousTableNum - 1;
    $HTMLFile = ${$SDFilesInfo{TableHTMLFiles}[$Index]}[$PreviousTableIndex];
    $HTMLRefFile = GetRelativeFileDir($CurTableNum, $PreviousTableNum, $FirstTableNum) . $HTMLFile;
    $TableStartCmpdNum = ${$SDFilesInfo{TableStartCmpdNum}[$Index]}[$PreviousTableIndex];
    $TableEndCmpdNum = ${$SDFilesInfo{TableEndCmpdNum}[$Index]}[$PreviousTableIndex];
    $HTMLRefValue = SetupHTMLHRef("Previous", $HTMLRefFile, "Previous Table Containing Compounds $TableStartCmpdNum To $TableEndCmpdNum");
    print $HTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue, $LinkBGColor, $LinkTextColor, $LinkFontBold);
  }

  for $TableNum ($StartTableNum .. $EndTableNum) {
    $TableIndex = $TableNum - 1;
    $HTMLFile = ${$SDFilesInfo{TableHTMLFiles}[$Index]}[$TableIndex];
    if ($TableNum == $CurTableNum) {
      print $HTMLFileRef SetupHTMLTableRowDataValue($TableNum, $LinkBGColor, $InactiveLinkNumColor, $InactiveLinkFontBold);
    }
    else {
      # Setup the link...
      my($RefTitle);
      $TableStartCmpdNum = ${$SDFilesInfo{TableStartCmpdNum}[$Index]}[$TableIndex];
      $TableEndCmpdNum = ${$SDFilesInfo{TableEndCmpdNum}[$Index]}[$TableIndex];
      $RefTitle = AddNumberSuffix($TableNum) . " Table Containing Compounds $TableStartCmpdNum To $TableEndCmpdNum";
      $HTMLRefFile = GetRelativeFileDir($CurTableNum, $TableNum, $FirstTableNum) . $HTMLFile;
      $HTMLRefValue = SetupHTMLHRef($TableNum, $HTMLRefFile, $RefTitle);
      print $HTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue);
    }
  }

  # Setup link to next table...
  if ($CurTableNum != $LastTableNum) {
    my($NextTableNum, $NextTableIndex);
    $NextTableNum = $CurTableNum + 1; $NextTableIndex = $NextTableNum - 1;
    $HTMLFile = ${$SDFilesInfo{TableHTMLFiles}[$Index]}[$NextTableIndex];
    $HTMLRefFile = GetRelativeFileDir($CurTableNum, $NextTableNum, $FirstTableNum) . $HTMLFile;
    $TableStartCmpdNum = ${$SDFilesInfo{TableStartCmpdNum}[$Index]}[$NextTableIndex];
    $TableEndCmpdNum = ${$SDFilesInfo{TableEndCmpdNum}[$Index]}[$NextTableIndex];
    $HTMLRefValue = SetupHTMLHRef("Next", $HTMLRefFile, "Next Table Containing Compounds $TableStartCmpdNum To $TableEndCmpdNum");
    print $HTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue, $LinkBGColor, $LinkTextColor, $LinkFontBold);
  }

  # Setup link to last table...
  if ($EndTableNum != $LastTableNum) {
    $HTMLFile = ${$SDFilesInfo{TableHTMLFiles}[$Index]}[$LastTableIndex];
    $HTMLRefFile = GetRelativeFileDir($CurTableNum, $LastTableNum, $FirstTableNum) . $HTMLFile;
    $TableStartCmpdNum = ${$SDFilesInfo{TableStartCmpdNum}[$Index]}[$LastTableIndex];
    $TableEndCmpdNum = ${$SDFilesInfo{TableEndCmpdNum}[$Index]}[$LastTableIndex];
    $HTMLRefValue = SetupHTMLHRef("Last", $HTMLRefFile, "Last Table Containing Compounds $TableStartCmpdNum To $TableEndCmpdNum");
    print $HTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue, $LinkBGColor, $LinkTextColor, $LinkFontBold);
  }
  # Setup current table info text....
  print $HTMLFileRef SetupHTMLTableRowDataValue("&nbsp");
  print $HTMLFileRef SetupHTMLTableRowDataValue("&nbsp");
  $TableStartCmpdNum = ${$SDFilesInfo{TableStartCmpdNum}[$Index]}[$CurTableNum - 1];
  $TableEndCmpdNum = ${$SDFilesInfo{TableEndCmpdNum}[$Index]}[$CurTableNum - 1];
  if ($OptionsInfo{NavLinksCmpdInfo}) {
    print $HTMLFileRef SetupHTMLTableRowDataValue("Showing compounds $TableStartCmpdNum to $TableEndCmpdNum of $LastCmpdNum");
  }
  else {
    print $HTMLFileRef SetupHTMLTableRowDataValue("Showing table $CurTableNum of $LastTableNum");
  }

  print $HTMLFileRef $SDFilesInfo{RowEndTags};
  # End link table...
  print $HTMLFileRef SetupHTMLTableEnd();
  print $HTMLFileRef SetupHTMLDivEnd();
  print $HTMLFileRef SetupHTMLAlignmentEnd("center");
}

# Generate relative directory path...
sub GetRelativeFileDir {
  my($FromTableNum, $ToTableNum, $FirstTableNum) = @_;
  my($RelativeFileDir) = "";

  if ($FromTableNum == $FirstTableNum) {
    $RelativeFileDir = ($ToTableNum == $FirstTableNum) ? ".\/" : ".\/html\/";
  }
  else {
    $RelativeFileDir = ($ToTableNum == $FirstTableNum) ? "..\/" : ".\/";
  }
  return $RelativeFileDir;
}

# Based on hightlight stype, return appropriate colors for background or text...
sub GetValueHighlightColors {
  my($Label, $Value) = @_;
  my($DataType, $Criterion, $CriterionValue, $BgColor, $FontColor, $ValueOk, $Nothing);

  $BgColor = ""; $FontColor = "";
  $DataType = $OptionsInfo{SpecifiedHighlightDataFieldTypesMap}{$Label};
  $Criterion = $OptionsInfo{SpecifiedHighlightDataFieldCriteriaMap}{$Label};
  $CriterionValue = $OptionsInfo{SpecifiedHighlightDataFieldValueMap}{$Label};

  $ValueOk = 0;
  if ($DataType =~ /^numeric$/i) {
  NUMSWITCH: {
      if ($Criterion =~ /^ge$/i) { $ValueOk = ($Value >= $CriterionValue) ? 1 : 0; last NUMSWITCH; }
      if ($Criterion =~ /^le$/i) { $ValueOk = ($Value <= $CriterionValue) ? 1 : 0; last NUMSWITCH; }
      if ($Criterion =~ /^eq$/i) { $ValueOk = ($Value == $CriterionValue) ? 1 : 0; last NUMSWITCH; }
      $Nothing = 1;
    }
  }
  else {
  TEXTSWITCH: {
      if ($Criterion =~ /^ge$/i) { $ValueOk = ($Value ge $CriterionValue) ? 1 : 0; last TEXTSWITCH; }
      if ($Criterion =~ /^le$/i) { $ValueOk = ($Value le $CriterionValue) ? 1 : 0; last TEXTSWITCH; }
      if ($Criterion =~ /^eq$/i) { $ValueOk = ($Value eq $CriterionValue) ? 1 : 0; last TEXTSWITCH; }
      $Nothing = 1;
    }
  }
  $BgColor = $ValueOk ? $OptionsInfo{ValueOkColor} : $OptionsInfo{ValueNotOkColor};
  if ($Options{highlightstyle} =~ /^text$/i) {
    $BgColor = "";
    $FontColor = $ValueOk ? $OptionsInfo{ValueOkColor} : $OptionsInfo{ValueNotOkColor};
  }
  return ($BgColor, $FontColor);
}

#Make sure appropriate mode specific option values are specified...
sub ProcessOptions {

  %OptionsInfo = ();

  $OptionsInfo{TitleDisplay} = ($Options{titledisplay} =~ /^yes$/i) ? 1 : 0;

  $OptionsInfo{RowHAlignment} = "left"; $OptionsInfo{RowVAlignment} = "middle";
  if (exists($Options{align})) {
    my (@AlignValues) = split ",", $Options{align};
    if (@AlignValues == 2) {
      $OptionsInfo{RowHAlignment} = $AlignValues[0];
      $OptionsInfo{RowVAlignment} = $AlignValues[1];
    }
    elsif (@AlignValues == 1) {
      $OptionsInfo{RowHAlignment} = $AlignValues[0];
    }
    else {
      die "Error: Invalid number of values, ", scalar(@AlignValues) , ", specified by \"-a --align\" option.\nIt must contain only one or two values.\n";
    }
    if ($OptionsInfo{RowHAlignment} !~ /^(left|center|right)$/i) {
      die "Error: The horizontal alignment value specified, $Options{align}, for option \"-a --align\" is not valid. Allowed values: left, center, or right\n";
    }
    if ($OptionsInfo{RowVAlignment} !~ /^(top|middle|bottom)$/i) {
      die "Error: The horizontal alignment value specified, $Options{align}, for option \"-a --align\" is not valid. Allowed values: top, middle, or bottom\n";
    }
  }

  $OptionsInfo{TableHeaderRowHAlignment} = "center"; $OptionsInfo{TableHeaderRowVAlignment} = "middle";
  if (exists($Options{headeralign})) {
    my (@AlignValues) = split ",", $Options{headeralign};
    if (@AlignValues == 2) {
      $OptionsInfo{TableHeaderRowHAlignment} = $AlignValues[0];
      $OptionsInfo{TableHeaderRowVAlignment} = $AlignValues[1];
    }
    elsif (@AlignValues == 1) {
      $OptionsInfo{TableHeaderRowHAlignment} = $AlignValues[0];
    }
    else {
      die "Error: Invalid number of values, ", scalar(@AlignValues) , ", specified by \"--headeralign\" option.\nIt must contain only one or two value.\n";
    }
    if ($OptionsInfo{TableHeaderRowHAlignment} !~ /^(left|center|right)$/i) {
      die "Error: The horizontal alignment value specified, $Options{headeralign}, for option \"--headeralign\" is not valid. Allowed values: left, center, or right\n";
    }
    if ($OptionsInfo{TableHeaderRowVAlignment} !~ /^(top|middle|bottom)$/i) {
      die "Error: The horizontal alignment value specified, $Options{headeralign}, for option \"-a --headeralign\" is not valid. Allowed values: top, middle, or bottom\n";
    }
  }

  if (exists($Options{border})) {
    $OptionsInfo{TableBorder} = $Options{border};
  }
  else {
    $OptionsInfo{TableBorder} = ($Options{mode} =~ /^(plain|highlight)$/i) || $Options{mode} =~ /^structuresonly$/i ? 1 : 0;
  }
  $OptionsInfo{TableCellPadding} = $Options{cellpadding};
  $OptionsInfo{TableCellSpacing} = $Options{cellspacing};
  $OptionsInfo{FooterMsg} = $Options{footer} ? $Options{footer} : "";

  if ($Options{headercolor}) {
    $OptionsInfo{TableHeaderRowColor} = $Options{headercolor};
  }
  else {
    $OptionsInfo{TableHeaderRowColor} = ($Options{mode} =~ /^plain$/i) ? "" : "#e0e9eb";
  }

  $OptionsInfo{NavLinksAtBottom} = 1; $OptionsInfo{NavLinksAtTop} = 0;
  if ($Options{displaylinks} =~ /^(both|top)$/i) {
    $OptionsInfo{NavLinksAtTop} = 1;
  }
  $OptionsInfo{NavLinksTableInfo} = 1; $OptionsInfo{NavLinksCmpdInfo} = 0;
  if ($Options{displaylinksinfo} =~ /^both$/i) {
    $OptionsInfo{NavLinksCmpdInfo} = 1;
    $OptionsInfo{NavLinksTableInfo} = 1;
  }
  elsif ($Options{displaylinksinfo} =~ /^compound$/i) {
    $OptionsInfo{NavLinksCmpdInfo} = 1;
    $OptionsInfo{NavLinksTableInfo} = 0;
  }

  if ($Options{stylesheet} =~ /^old$/i ) {
    if (!$Options{stylesheetname}) {
      die "Error: No stylesheet name specified using \"--stylesheetname\" option: It is required for \"old\" value of \"-s --stylesheet\" option. \n";
    }
  }

  my(@ColorValues);
  $OptionsInfo{ShadeRowsStatus} = 0;
  $OptionsInfo{OddRowsShadeColor} = "#ffffff";
  $OptionsInfo{EvenRowsShadeColor} = "#e0e9eb";
  if ($Options{shadecolor}) {
    # Make sure only one value is specified...
    @ColorValues = split ",", $Options{shadecolor};
    if (@ColorValues == 2) {
      $OptionsInfo{OddRowsShadeColor} = $ColorValues[0];
      $OptionsInfo{EvenRowsShadeColor} = $ColorValues[1];
    }
    else {
      die "Error: Invalid number of values, ", scalar(@ColorValues) , ", specified by \"--shadecolor\" option.\nIt must contain only two value.\n";
    }
  }
  if ($Options{mode} =~ /^(shade|shadedhighlight|shadedstructuresonly)$/i) {
    $OptionsInfo{ShadeRowsStatus} = 1;
  }

  $OptionsInfo{SettingUpCmpdSummaryPage} = 0;
  $OptionsInfo{StrLinkShadeColor} = (exists $Options{strlinkshadecolor}) ? $Options{strlinkshadecolor} : "#e0e9eb";
  $OptionsInfo{DisplayStructure} = ($Options{structure} =~ /^display$/i) ? 1 : 0;
  $OptionsInfo{StrViewerType} = $Options{strviewertype};
  $OptionsInfo{StrLinkNavigation} = ($Options{strlinknavigation} =~ /^yes$/i) ? 1 : 0;
  $OptionsInfo{StrLinkTitleDisplay} = ($Options{strlinktitledisplay} =~ /^yes$/i) ? 1 : 0;
  $OptionsInfo{StrLinkTitle} = (exists($Options{strlinktitle}) && length($Options{strlinktitle})) ? "$Options{strlinktitle}" : "Compound Summary";

  my($StrViewerEmbedUsingJS) = (($Options{strviewerembed} =~ /^javascript$/i) && ($OptionsInfo{StrViewerType} =~ /^(Jmol|MarvinView|ChemDrawPlugIn|ChemDrawActiveX|Chem3DActiveX)$/i )) ? 1 : 0;

  $OptionsInfo{StrTableRows} = 6; $OptionsInfo{StrTableCols} = 4;
  if ($Options{strtablesize}) {
    my(@StrTableSizeValues) = split ",", $Options{strtablesize};
    if (@StrTableSizeValues == 2) {
      $OptionsInfo{StrTableRows} = $StrTableSizeValues[0];
      $OptionsInfo{StrTableCols} = $StrTableSizeValues[1];
      if (!IsPositiveInteger($OptionsInfo{StrTableRows})) {
	die "Error: The first value specified, $OptionsInfo{StrTableRows},  for option \"--strtablesize\" is not valid: Allowed integer values: > 0.\n";
      }
      if (!IsPositiveInteger($OptionsInfo{StrTableCols})) {
	die "Error: The first value specified, $OptionsInfo{StrTableCols},  for option \"--strtablesize\" is not valid: Allowed integer values: > 0.\n";
      }
    }
    else {
      die "Error: Invalid number of values, ", scalar(@StrTableSizeValues), ", specified by \"--strtablesize\" option.\nIt must contain only two value for structuresonly \"-m --mode\" option.\n";
    }
  }

  # Setup applet information...
  $OptionsInfo{StrViewerCodeBase} = GetMayaChemToolsLibDirName() . "/Jmol";
  $OptionsInfo{TopHTMLDirStrViewerCodeBase} = $OptionsInfo{StrViewerCodeBase};
  $OptionsInfo{SubHTMLDirStrViewerCodeBase} = $OptionsInfo{StrViewerCodeBase};

  my($StrViewerAppletArchive, $StrViewerAppletCode) = SetupDefaultAppletArchiveAndCode($OptionsInfo{StrViewerType});
  if ($Options{strviewerconfig}) {
    my(@StrViewerConfigParts) = split ",", $Options{strviewerconfig};
    if (@StrViewerConfigParts >=1 && @StrViewerConfigParts <= 3) {
      if (@StrViewerConfigParts == 3) {
	$OptionsInfo{StrViewerCodeBase} = $StrViewerConfigParts[0];
	$StrViewerAppletArchive = $StrViewerConfigParts[1];
	$StrViewerAppletCode = $StrViewerConfigParts[2];
      }
      elsif (@StrViewerConfigParts == 2) {
	$OptionsInfo{StrViewerCodeBase} = $StrViewerConfigParts[0];
	$StrViewerAppletArchive = $StrViewerConfigParts[1];
	my($AppletArchive, $AppletCode) = SetupDefaultAppletArchiveAndCode($OptionsInfo{StrViewerType});
	$StrViewerAppletCode = $AppletCode;
      }
      else {
	$OptionsInfo{StrViewerCodeBase} = $StrViewerConfigParts[0];
	($StrViewerAppletArchive, $StrViewerAppletCode) = SetupDefaultAppletArchiveAndCode($OptionsInfo{StrViewerType});
      }
    }
    else {
      die "Error: Invalid number of values, ", scalar(@StrViewerConfigParts), ", specified by \"--strviewerconfig\" option.\nNumver of allowed values:1 to 3 \n";
    }
  }
  else {
    if ($OptionsInfo{StrViewerType} =~ /^(JME|MarvinView)$/i ) {
      die "Error: No codebase specified using \"--strviewerconfig\" option for $OptionsInfo{StrViewerType} structure viewer\n";
    }
    if ($StrViewerEmbedUsingJS && $OptionsInfo{StrViewerType} !~ /^Jmol$/i) {
      die "Error: No codebase specified using \"--strviewerconfig\" option for javascript value of \"--strviewerembed\" option for $OptionsInfo{StrViewerType} structure viewer \n";
    }
  }

  if (-d $OptionsInfo{StrViewerCodeBase}) {
    # Change local code base direcrory name to a relative directory name based on the
    # current directory containing SD file; otherwise, Java applets and JavaScripts don't
    # get loaded into Firefox and Chrome browsers.
    #
    # For top and sub HTML directories, add prefix "../" and "../../" to relative path...
    $OptionsInfo{StrViewerCodeBase} = File::Spec->abs2rel($OptionsInfo{StrViewerCodeBase}, Cwd::cwd());

    $OptionsInfo{TopHTMLDirStrViewerCodeBase} = "../" . $OptionsInfo{StrViewerCodeBase};
    $OptionsInfo{SubHTMLDirStrViewerCodeBase} = "../../" . $OptionsInfo{StrViewerCodeBase};
  }

  # Setup structure viewer parameter information...
  %{$OptionsInfo{StrViewerParams}} = ();
  if ($Options{strviewerparams}) {
    my(@ParamsSplit, @ParamPairSplit, $ParamPair);
    #@ParamsSplit = split " ", $Options{strviewerparams};
    @ParamsSplit = quotewords(" ", 0, $Options{strviewerparams});
    for $ParamPair (@ParamsSplit) {
      @ParamPairSplit = split "=", $ParamPair;
      if (@ParamPairSplit == 2) {
	$OptionsInfo{StrViewerParams}{$ParamPairSplit[0]} = $ParamPairSplit[1];
      }
      else {
	die "Error: Invalid value, $ParamPair, specified by \"--strviewerparams\" option.\nValid values: name=value\n";
      }
    }
  }

  if ($OptionsInfo{StrViewerType} =~ /^(JME|Jmol|MarvinView)$/i ) {
    $OptionsInfo{StrViewerParams}{name} = $StrViewerAppletCode;
    $OptionsInfo{StrViewerParams}{archive} = $StrViewerAppletArchive;
    $OptionsInfo{StrViewerParams}{code} = $StrViewerAppletCode;
  }
  $OptionsInfo{StrWidth} = exists($OptionsInfo{StrViewerParams}{width}) ? $OptionsInfo{StrViewerParams}{width} : 250;
  $OptionsInfo{StrViewerParams}{width} = $OptionsInfo{StrWidth};
  $OptionsInfo{StrHeight} = exists($OptionsInfo{StrViewerParams}{height}) ? $OptionsInfo{StrViewerParams}{height} : 170;
  $OptionsInfo{StrViewerParams}{height} = $OptionsInfo{StrHeight};

  $OptionsInfo{StrLinkWidth} = 500;
  if (exists($OptionsInfo{StrViewerParams}{strlinkwidth})) {
    $OptionsInfo{StrLinkWidth} = $OptionsInfo{StrViewerParams}{strlinkwidth};
    $OptionsInfo{StrViewerParams}{strlinkwidth} = "";
  }
  $OptionsInfo{StrLinkHeight} = 295;
  if (exists($OptionsInfo{StrViewerParams}{strlinkheight})) {
    $OptionsInfo{StrLinkHeight} = $OptionsInfo{StrViewerParams}{strlinkheight};
    $OptionsInfo{StrViewerParams}{strlinkheight} = "";
  }

  $OptionsInfo{StrBgColorSpecified} = "";
  if (exists($OptionsInfo{StrViewerParams}{bgcolor})) {
    $OptionsInfo{StrBgColorSpecified} = $OptionsInfo{StrViewerParams}{bgcolor};
  }

  $OptionsInfo{StrLinkBgColorSpecified} = "#ffffff";
  if (exists($OptionsInfo{StrViewerParams}{strlinkbgcolor})) {
    $OptionsInfo{StrLinkBgColorSpecified} = $OptionsInfo{StrViewerParams}{strlinkbgcolor};
    $OptionsInfo{StrViewerParams}{strlinkbgcolor} = "";
  }

  # Setup Java Script usage...
  $OptionsInfo{StrViewerJSFileRef} = "";
  $OptionsInfo{TopHTMLDirStrViewerJSFileRef} = "";
  $OptionsInfo{SubHTMLDirStrViewerJSFileRef} = "";

  if ($StrViewerEmbedUsingJS) {
    my ($StrViewerJSFileName) = "";
    if ($Options{strviewerjsfile}) {
      $StrViewerJSFileName = $Options{strviewerjsfile};
    }
    else {
      if ($OptionsInfo{StrViewerType} =~ /^Jmol$/i) {
	$StrViewerJSFileName = "Jmol.js";
      }
      elsif ($OptionsInfo{StrViewerType} =~ /^MarvinView$/i) {
	$StrViewerJSFileName = "marvin.js";
      }
      elsif ($OptionsInfo{StrViewerType} =~ /^(ChemDrawPlugIn|ChemDrawActiveX)$/i) {
	$StrViewerJSFileName = "chemdraw.js";
      }
      elsif ($OptionsInfo{StrViewerType} =~ /^Chem3DActiveX$/i) {
	$StrViewerJSFileName = "chem3d.js";
      }
    }
    if ($StrViewerJSFileName) {
      $OptionsInfo{StrViewerParams}{usejavascript} = $StrViewerJSFileName;
      $OptionsInfo{StrViewerJSFileRef} = "$OptionsInfo{StrViewerCodeBase}" . "\/" . "$StrViewerJSFileName";
      $OptionsInfo{TopHTMLDirStrViewerJSFileRef} = "$OptionsInfo{TopHTMLDirStrViewerCodeBase}" . "\/" . "$StrViewerJSFileName";
      $OptionsInfo{SubHTMLDirStrViewerJSFileRef} = "$OptionsInfo{SubHTMLDirStrViewerCodeBase}" . "\/" . "$StrViewerJSFileName";
    }
  }

  # Check any other user specified parametrs applicable to all structure viewers...

  $OptionsInfo{StructuresOnlyMode} = 0;
  $OptionsInfo{MaxCmpdsPerTable} = ($Options{structure} =~ /^display$/i) ? 15 : 50;
  if (exists $Options{numcmpds}) {
    $OptionsInfo{MaxCmpdsPerTable} = $Options{numcmpds};
  }
  if ($Options{mode} =~ /^(structuresonly|shadedstructuresonly)$/i) {
    $OptionsInfo{MaxCmpdsPerTable} = ($OptionsInfo{MaxCmpdsPerTable} > 0) ? ($OptionsInfo{StrTableRows} * $OptionsInfo{StrTableCols}) : 0;
    $OptionsInfo{StructuresOnlyMode} = 1;
  }
  $OptionsInfo{CmpdDataField} = "";
  $OptionsInfo{CmpdDataFieldLabel} = "no";
  $OptionsInfo{CmpdDataFieldPosition} = "bottom";
  $OptionsInfo{CmpdDataFieldAlignment} = "center";
  if (exists($Options{cmpddatafield}) && length($Options{cmpddatafield})) {
    my (@CmpdDataFieldValues) = split ",", $Options{cmpddatafield};
    if (@CmpdDataFieldValues == 1) {
      $OptionsInfo{CmpdDataField} = $CmpdDataFieldValues[0];
    }
    elsif (@CmpdDataFieldValues == 2) {
      $OptionsInfo{CmpdDataField} = $CmpdDataFieldValues[0];
      $OptionsInfo{CmpdDataFieldLabel} = $CmpdDataFieldValues[1];
    }
    elsif (@CmpdDataFieldValues == 3) {
      $OptionsInfo{CmpdDataField} = $CmpdDataFieldValues[0];
      $OptionsInfo{CmpdDataFieldLabel} = $CmpdDataFieldValues[1];
      $OptionsInfo{CmpdDataFieldPosition} = $CmpdDataFieldValues[2];
    }
    elsif (@CmpdDataFieldValues == 4) {
      $OptionsInfo{CmpdDataField}  = $CmpdDataFieldValues[0];
      $OptionsInfo{CmpdDataFieldLabel} = $CmpdDataFieldValues[1];
      $OptionsInfo{CmpdDataFieldPosition} = $CmpdDataFieldValues[2];
      $OptionsInfo{CmpdDataFieldAlignment} = $CmpdDataFieldValues[3];
    }
    else {
      die "Error: Invalid number of values, ", scalar(@CmpdDataFieldValues) , ", specified by \"--cmpddatafield\" option.\nIt must contain only one, two, three, or four values.\n";
    }
    if ($OptionsInfo{CmpdDataFieldLabel} !~ /^(yes|no)$/ ) {
      die "Error: The label value specified, $Options{cmpddatafield}, for option \"--cmpddatafield\" is not valid. Allowed values: yes or no\n";
    }
    if ($OptionsInfo{CmpdDataFieldPosition} !~ /^(top|bottom)$/ ) {
      die "Error: The position value specified, $Options{cmpddatafield}, for option \"--cmpddatafield\" is not valid. Allowed values: top or bottom\n";
    }
    if ($OptionsInfo{CmpdDataFieldAlignment} !~ /^(left|center|right)$/ ) {
      die "Error: The alignment value specified, $Options{cmpddatafield}, for option \"--cmpddatafield\" is not valid. Allowed values: left, center, or right\n";
    }
  }

  # Process data fields to be displayed in tables...
  $OptionsInfo{SpecifiedDataFields} = exists($Options{datafields}) ? $Options{datafields} : "All";

  $OptionsInfo{ValueOkColor} = ""; $OptionsInfo{ValueNotOkColor} = ""; $OptionsInfo{HighlightStatus} = 0;
  if ($Options{mode} =~ /^(highlight|shadedhighlight)$/i) {
    my($HighlightMode, $HighlightBy);
    $HighlightMode = $Options{mode}; $HighlightBy = $Options{highlightby};

    $OptionsInfo{HighlightStatus} = 1;
    $OptionsInfo{ValueOkColor} = "#0fff0f";
    $OptionsInfo{ValueNotOkColor} = "#ff0f0f";
    if ($Options{highlightstyle} =~ /^text$/i) {
      $OptionsInfo{ValueOkColor} = "#0fbb0f";
      $OptionsInfo{ValueNotOkColor} = "#ff0f0f";
    }
    if ($Options{highlightcolor}) {
      # Make sure two values are specified...
      @ColorValues = split ",", $Options{highlightcolor};
      if (@ColorValues == 2) {
	$OptionsInfo{ValueOkColor} = $ColorValues[0];
	$OptionsInfo{ValueNotOkColor} = $ColorValues[1];
      }
      else {
	die "Error: Invalid number of values, ", scalar(@ColorValues), ", specified by \"--highlightcolor\" option.\nIt must contain only two value for $HighlightMode value specified using \"-m --mode\" option.\n";
      }
    }
    if (!$Options{highlight}) {
      die "Error: Specify columns to be highlighted using \"--hightlight\" option\n";
    }
    # Retrieve quartet values from "hightlight" option...
    my(@HighlightValueQuartets);

    @HighlightValueQuartets = ();
    @HighlightValueQuartets = split ",", $Options{highlight};
    if ((@HighlightValueQuartets % 4)) {
      die "Error: Quartets not found in values specified using \"--highlight\" option for $HighlightMode \"-m --mode\"\n";
    }
    # Process quartets...
    my($Index, $Label, $DataType, $Criterion, $Value);

    @{$OptionsInfo{SpecifiedHighlightDataFieldLabels}} = ();
    %{$OptionsInfo{SpecifiedHighlightDataFieldLabelsMap}} = ();
    %{$OptionsInfo{SpecifiedHighlightDataFieldTypesMap}} = ();
    %{$OptionsInfo{SpecifiedHighlightDataFieldCriteriaMap}} = ();
    %{$OptionsInfo{SpecifiedHighlightDataFieldValueMap}} = ();

    for ($Index = 0; $Index < @HighlightValueQuartets; $Index = $Index + 4) {
      $Label = $HighlightValueQuartets[$Index];
      $DataType = $HighlightValueQuartets[$Index + 1];
      $Criterion = $HighlightValueQuartets[$Index + 2];
      $Value = $HighlightValueQuartets[$Index + 3];
      if ($DataType !~ /^(numeric|text)$/i) {
	die "Error: Invalid column data type, $DataType, specified in quartet, \"$Label,$DataType,$Criterion,$Value\", using \"--hightlight\" option: Valid values: numeric or text\n";
      }
      if ($Criterion !~ /^(eq|le|ge)$/i) {
	die "Error: Invalid criterion value, $Criterion, specified in quartet, \"$Label,$DataType,$Criterion,$Value\", using \"--hightlight\" option: Valid values: le, ge, or eq\n";
      }
      if ($DataType =~ /^numeric$/i) {
	if (!IsFloat($Value)) {
	  die "Error: Invalid criterion value, $Value, specified in quartet, \"$Label,$DataType,$Criterion,$Value\", using \"--hightlight\" option: Numeric value required for numeric data type\n";
	}
      }
      if (exists($OptionsInfo{SpecifiedHighlightDataFieldLabelsMap}{$Label})) {
	die "Error: Invalid field label value, $Label, in quartet, \"$Label,$DataType,$Criterion,$Value\", using \"--hightlight\" option: Multiple occurences of label.  \n";
      }
      push @{$OptionsInfo{SpecifiedHighlightDataFieldLabels}}, $Label;
      $OptionsInfo{SpecifiedHighlightDataFieldLabelsMap}{$Label} = $Label;
      $OptionsInfo{SpecifiedHighlightDataFieldTypesMap}{$Label} = $DataType;
      $OptionsInfo{SpecifiedHighlightDataFieldCriteriaMap}{$Label} = $Criterion;
      $OptionsInfo{SpecifiedHighlightDataFieldValueMap}{$Label} = $Value;
    }
  }
}

# Set up default archive and code values for a specific applet...
sub SetupDefaultAppletArchiveAndCode {
  my($ViewerType) = @_;
  my($Archive, $Code, $Nothing);

 STRVIEWERTYPE: {
    if ($OptionsInfo{StrViewerType} =~ /^JME$/i) { $Archive = "JME.jar"; $Code = "JME"; last STRVIEWERTYPE; }
    if ($OptionsInfo{StrViewerType} =~ /^Jmol$/i) {$Archive = "JmolApplet.jar"; $Code = "JmolApplet";  last STRVIEWERTYPE; }
    if ($OptionsInfo{StrViewerType} =~ /^MarvinView$/i) { $Archive = "marvin.jar"; $Code = "MView"; last STRVIEWERTYPE; }
    $Nothing = 1;
  }
  return ($Archive, $Code);
}

# Retrieve information about input SD files...
sub RetrieveSDFilesInfo {
  my($SDFile, $FileDir, $FileName, $HTMLFile, $CSSFile, $HTMLRoot, $HTMLTitle, $FileExt, $Index, $TopHTMLDir);

  %SDFilesInfo = ();

  @{$SDFilesInfo{FileOkay}} = ();
  @{$SDFilesInfo{CmpdCount}} = ();
  @{$SDFilesInfo{SpecifiedDataFieldLabels}} = ();

  @{$SDFilesInfo{HTMLRoot}} = ();
  @{$SDFilesInfo{HTMLTitle}} = ();
  @{$SDFilesInfo{MultipleHTMLTables}} = ();

  @{$SDFilesInfo{TopHTMLDir}} = ();
  @{$SDFilesInfo{SubHTMLDir}} = ();
  @{$SDFilesInfo{SubMolDir}} = ();


  FILELIST: for $Index (0 .. $#SDFilesList) {
    $SDFile = $SDFilesList[$Index];

    $SDFilesInfo{FileOkay}[$Index] = 0;
    $SDFilesInfo{CmpdCount}[$Index] = 0;
    $SDFilesInfo{HTMLRoot}[$Index] = "";
    $SDFilesInfo{HTMLTitle}[$Index] = "";
    $SDFilesInfo{MultipleHTMLTables}[$Index] = 0;
    $SDFilesInfo{TopHTMLDir}[$Index] = "";
    $SDFilesInfo{SubHTMLDir}[$Index] = "";
    $SDFilesInfo{SubMolDir}[$Index] = "";

    @{$SDFilesInfo{SpecifiedDataFieldLabels}[$Index]} = ();

    if (!(-e $SDFile)) {
      warn "Warning: Ignoring file $SDFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($SDFile, "sd sdf")) {
      warn "Warning: Ignoring file $SDFile: It's not a SD file\n";
      next FILELIST;
    }
    ($FileDir, $FileName, $FileExt) = ParseFileName($SDFile);

    if (!open SDFILE, "$SDFile") {
      warn "Warning: Ignoring file $SDFile: Couldn't open it: $! \n";
      next FILELIST;
    }
    # Count number of compounds and collect all possible data field labels...
    my($CmpdCount, $CmpdString, @DataFieldLabels, @CommonDataFieldLabels);
    $CmpdCount = 0;
    @DataFieldLabels = ();
    @CommonDataFieldLabels = ();
    if ($OptionsInfo{SpecifiedDataFields} =~ /^(All|Common)$/i ) {
      my($DataFieldLabelsRef, $CommonDataFieldLabelsRef);
      ($CmpdCount, $DataFieldLabelsRef, $CommonDataFieldLabelsRef) = GetAllAndCommonCmpdDataHeaderLabels(\*SDFILE);
      push @DataFieldLabels, @{$DataFieldLabelsRef};
      push @CommonDataFieldLabels, @{$CommonDataFieldLabelsRef};
    }
    else {
      while ($CmpdString = ReadCmpdString(\*SDFILE)) {
	$CmpdCount++;
      }
    }
    close SDFILE;

    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($SDFile);
    $HTMLRoot = $FileName;
    if ($Options{root} && (@SDFilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($Options{root});
      if ($RootFileName && $RootFileExt) {
	$HTMLRoot = $RootFileName;
      }
      else {
	$HTMLRoot = $Options{root};
      }
    }
    $HTMLTitle = $HTMLRoot;
    if ($Options{title} && (@SDFilesList == 1)) {
      $HTMLTitle = $Options{title};
    }
    $HTMLFile = lc($HTMLRoot) . "-html";
    if (!$Options{overwrite}) {
      if (-d $HTMLFile) {
	warn "Warning: Ignoring file $SDFile: The directory $HTMLFile already exists\n";
	next FILELIST;
      }
    }
    $SDFilesInfo{FileOkay}[$Index] = 1;
    $SDFilesInfo{CmpdCount}[$Index] = $CmpdCount;
    $SDFilesInfo{HTMLRoot}[$Index] = "$HTMLRoot";
    $SDFilesInfo{HTMLTitle}[$Index] = "$HTMLTitle";
    if ($OptionsInfo{MaxCmpdsPerTable} == 0 || $CmpdCount <= $OptionsInfo{MaxCmpdsPerTable}) {
      $SDFilesInfo{MultipleHTMLTables}[$Index] = 0;
    }
    else {
      $SDFilesInfo{MultipleHTMLTables}[$Index] = 1;
    }
    if ($OptionsInfo{SpecifiedDataFields} =~ /^All$/i ) {
      push @{$SDFilesInfo{SpecifiedDataFieldLabels}[$Index]}, @DataFieldLabels;
    }
    elsif ($OptionsInfo{SpecifiedDataFields} =~ /^Common$/i) {
      push @{$SDFilesInfo{SpecifiedDataFieldLabels}[$Index]}, @CommonDataFieldLabels;
    }
    else {
      push @{$SDFilesInfo{SpecifiedDataFieldLabels}[$Index]}, split(",", $OptionsInfo{SpecifiedDataFields});
    }

    # Setup HTML data directories paths...
    $TopHTMLDir = lc($SDFilesInfo{HTMLRoot}[$Index]) . "-html";
    $SDFilesInfo{TopHTMLDir}[$Index] = "$TopHTMLDir";
    $SDFilesInfo{SubHTMLDir}[$Index] = "$TopHTMLDir\/html";
    $SDFilesInfo{SubMolDir}[$Index] = "$TopHTMLDir\/mol";
  }
}

# Setup information...
sub SetupMultipleTablesAndMiscInfo {
  SetupMultipleTablesInfo();
  SetupMiscInfo();
}

# Setup navigation link information for multiple tables...
sub SetupMultipleTablesInfo {
  my($Index, $LinesPerTable);

  $LinesPerTable = $OptionsInfo{MaxCmpdsPerTable};

  @{$SDFilesInfo{TableCount}} = ();
  @{$SDFilesInfo{TableHTMLFiles}} = ();
  @{$SDFilesInfo{TableStartCmpdNum}} = ();
  @{$SDFilesInfo{TableEndCmpdNum}} = ();

  for $Index (0 .. $#SDFilesList) {
    $SDFilesInfo{TableCount}[$Index] = 1;
    @{$SDFilesInfo{TableHTMLFiles}[$Index]} = ();
    @{$SDFilesInfo{TableStartCmpdNum}[$Index]} = ();
    @{$SDFilesInfo{TableEndCmpdNum}[$Index]} = ();

    if ($SDFilesInfo{FileOkay}[$Index]) {
      if ($SDFilesInfo{MultipleHTMLTables}[$Index]) {
	my($TableIndex, $TotalLines, $TableCount, $TableStartLineNum, $TableEndLineNum, $Name);

	$TotalLines = $SDFilesInfo{CmpdCount}[$Index];
	$TableCount = ($TotalLines % $LinesPerTable) ? (int($TotalLines/$LinesPerTable) + 1) : ($TotalLines/$LinesPerTable);
	$SDFilesInfo{TableCount}[$Index] = $TableCount;
	for $TableIndex (1 .. $TableCount) {
	  $TableStartLineNum = ($TableIndex - 1) * $LinesPerTable + 1;
	  $TableEndLineNum = ($TableIndex == $TableCount) ? $TotalLines : ($TableIndex * $LinesPerTable);
	  push @{$SDFilesInfo{TableStartCmpdNum}[$Index]}, $TableStartLineNum;
	  push @{$SDFilesInfo{TableEndCmpdNum}[$Index]}, $TableEndLineNum;

	  # Setup HTML file names for all the tables...
	  $Name = "Cmpd" . "$TableStartLineNum" . "To" . "$TableEndLineNum";
	  if ($TableIndex == 1) {
	    $Name = "";
	  }
	  $Name = $SDFilesInfo{HTMLRoot}[$Index] . $Name . ".html";
	  push @{$SDFilesInfo{TableHTMLFiles}[$Index]}, $Name;
	}
	#print "$SDFilesList[$Index]: $TableCount -  @{$SDFilesInfo{TableStartCmpdNum}[$Index]} - @{$SDFilesInfo{TableEndCmpdNum}[$Index]} -  @{$SDFilesInfo{TableHTMLFiles}[$Index]}\n";
      }
    }
  }
}

# Setup HTML tags and other information...
sub SetupMiscInfo {
  $SDFilesInfo{RowHeaderTags} = "";
  $SDFilesInfo{RowEndTags} = "";
  $SDFilesInfo{BgFilledOddRowHeaderTags} = "";
  $SDFilesInfo{BgFilledEvenRowHeaderTags} = "";
  $SDFilesInfo{TableRowHeaderTags} = "";

  $SDFilesInfo{RowHeaderTags} = SetupHTMLTableRowHeader($OptionsInfo{RowHAlignment}, "", $OptionsInfo{RowVAlignment});
  $SDFilesInfo{RowEndTags} = SetupHTMLTableRowEnd();

  if ($OptionsInfo{ShadeRowsStatus}) {
    $SDFilesInfo{BgFilledOddRowHeaderTags} = SetupHTMLTableRowHeader($OptionsInfo{RowHAlignment}, $OptionsInfo{OddRowsShadeColor}, $OptionsInfo{RowVAlignment});
    $SDFilesInfo{BgFilledEvenRowHeaderTags} = SetupHTMLTableRowHeader($OptionsInfo{RowHAlignment}, $OptionsInfo{EvenRowsShadeColor}, $OptionsInfo{RowVAlignment});
  }

  $SDFilesInfo{TableRowHeaderTags} = SetupHTMLTableRowHeader($OptionsInfo{TableHeaderRowHAlignment}, $OptionsInfo{TableHeaderRowColor}, $OptionsInfo{TableHeaderRowVAlignment});

  $SDFilesInfo{MolEndTag} = "M  END";
  $SDFilesInfo{White} = qq(#ffffff);
}

# Setup various data directories to hold HTML and other related files...
sub SetupDataDirs {
  my($Index) = @_;
  my($TopHTMLDir, $SubHTMLDir, $SubMolDir, $CreateTopHTMLDir, $CreateSubHTMLDir, $CreateSubMolDir);

  $TopHTMLDir = $SDFilesInfo{TopHTMLDir}[$Index];
  $SubHTMLDir = $SDFilesInfo{SubHTMLDir}[$Index];
  $SubMolDir = $SDFilesInfo{SubMolDir}[$Index];

  # Clean up existing directories...
  if (-d $TopHTMLDir) {
    unlink "<$TopHTMLDir/*.html>";
    unlink "<$TopHTMLDir/*.css>";
  }
  if (-d $SubHTMLDir) {
    unlink "<$SubHTMLDir/*.html>";
  }
  if (-d $SubMolDir) {
    unlink "<$SubMolDir/*.mol>";
  }

  # What directories need to be created...
  $CreateTopHTMLDir = (-d $TopHTMLDir) ? 0 : 1;
  $CreateSubHTMLDir = (-d $SubHTMLDir) ? 0 : 1;
  $CreateSubMolDir = 0;
  if ($OptionsInfo{StrViewerType} =~ /^(Jmol|Chime|Chem3DActiveX|ChemDrawActiveX|ChemDrawPlugIn|ViewerActiveX)$/i) {
    $CreateSubMolDir = (-d $SubMolDir) ? 0 : 1;
  }

  # Create appropriate directories...
  if ($CreateTopHTMLDir) {
    mkdir $TopHTMLDir or die "Couldn't mkdir $TopHTMLDir: $! \n";
  }
  if ($CreateSubHTMLDir) {
    mkdir $SubHTMLDir or die "Error: Couldn't mkdir $SubHTMLDir: $! \n";
  }
  else {
    unlink <$SubHTMLDir/*.html>;
  }
  if ($CreateSubMolDir) {
    mkdir $SubMolDir or die "Error: Couldn't mkdir $SubMolDir: $! \n";
  }
  return ($TopHTMLDir, $SubHTMLDir, $SubMolDir);
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{mode} = "shade";
  $Options{highlightstyle} = "background";

  $Options{cellpadding} = 2;
  $Options{cellspacing} = 1;

  $Options{displaylinks} = "both";
  $Options{displaylinksinfo} = "both";
  $Options{stylesheet} = "new";

  $Options{structure} = "display";
  $Options{strlinktype} = "href";
  $Options{strlinkmode} = "plain";
  $Options{strlinknavigation} = "yes";
  $Options{strlinktitledisplay} = "no";

  $Options{strviewertype} = "Jmol";
  $Options{strviewerembed} = "direct";

  $Options{titledisplay} = "yes";

  if (!GetOptions(\%Options, "align|a=s", "border|b=i", "cellpadding=i", "cellspacing=i", "cmpddatafield|c=s", "datafields=s", "footer=s", "displaylinks|d=s", "displaylinksinfo=s", "help|h", "headeralign=s", "headercolor=s", "highlight=s", "highlightcolor=s", "highlightstyle=s", "mode|m=s", "numcmpds|n=i", "overwrite|o", "root|r=s", "shadecolor=s", "stylesheet=s", "stylesheetname=s", "structure|s=s", "strlinkmode=s", "strlinknavigation=s", "strlinkshadecolor=s", "strlinktitle=s", "strlinktitledisplay=s", "strlinktype=s", "strviewertype=s", "strviewerconfig=s", "strviewerparams=s", "strviewerembed=s",  "strviewerjsfile=s", "strtablesize=s", "title|t=s", "titledisplay=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }

  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{displaylinks} !~ /^(top|bottom|both)$/i) {
    die "Error: The value specified, $Options{displaylinks}, for option \"-d --displaylinks\" is not valid. Allowed values: top, bottom, or both\n";
  }
  if ($Options{displaylinksinfo} !~ /^(compound|table|both)$/i) {
    die "Error: The value specified, $Options{displaylinksinfo}, for option \"--displaylinksinfo\" is not valid. Allowed values: compound, table, or both\n";
  }
  if ($Options{highlightstyle} !~ /^(background|text)$/i) {
    die "Error: The value specified, $Options{highlightstyle}, for option \"--highlightstyle\" is not valid. Allowed values: background or text\n";
  }
  if ($Options{mode} !~ /^(plain|shade|highlight|shadedhighlight|structuresonly|shadedstructuresonly)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: plain, shade, hightlight, shadedhighlight, structuresonly, or shadedstructuresonly\n";
  }
  if ($Options{stylesheet} !~ /^(old|new|none)$/i) {
    die "Error: The value specified, $Options{stylesheet}, for option \"-s --stylesheet\" is not valid. Allowed values: old, new, or none\n";
  }
  if ($Options{structure} !~ /^(display|link)$/i) {
    die "Error: The value specified, $Options{structure}, for option \"-s --structure\" is not valid. Allowed values: display or link\n";
  }
  if ($Options{strlinkmode} !~ /^(plain|shade)$/i) {
    die "Error: The value specified, $Options{strlinkmode}, for option \"--strlinkmode\" is not valid. Allowed values: plain or shade\n";
  }
  if ($Options{strlinktype} !~ /^(href|button)$/i) {
    die "Error: The value specified, $Options{strlinktype}, for option \"--strlinktype\" is not valid. Allowed values: href or button\n";
  }
  if ($Options{strlinknavigation} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{strlinknavigation}, for option \"--strlinknavigation\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{strlinktitledisplay} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{strlinktitledisplay}, for option \"--strlinktitledisplay\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{strviewertype} !~ /^(JME|Jmol|Chime|MarvinView|ChemDrawPlugIn|Chem3DActiveX|ChemDrawActiveX|ViewerActiveX)$/i) {
    die "Error: The value specified, $Options{strviewertype}, for option \"--strviewertype\" is not valid. Allowed values: Chem3DActiveX, ChemDrawActiveX, ChemDrawPlugIn, Chime, JME, Jmol, MarvinView, or ViewerActiveX.\n";
  }
  if ($Options{strviewerembed} !~ /^(direct|javascript)$/i) {
    die "Error: The value specified, $Options{strviewerembed},  for option \"--strviewerembed\" is not valid. Allowed values: direct or javascript \n";
  }
  if (exists $Options{numcmpds} && $Options{numcmpds} < 0) {
    die "Error: The value specified, $Options{numcmpds},  for option \"-n --numcmpds\" is not valid. Allowed values: >= 0 \n";
  }
  if ($Options{titledisplay} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{titledisplay}, for option \"--titledisplay\" is not valid. Allowed values: yes or no\n";
  }
  if (exists($Options{border})) {
    if ($Options{border} < 0) {
      die "Error: The value specified, $Options{border},  for option \"--border\" is not valid. Allowed values: >= 0 \n";
    }
  }
  if ($Options{cellpadding} < 0) {
    die "Error: The value specified, $Options{cellpadding},  for option \"--cellpadding\" is not valid. Allowed values: >= 0 \n";
  }
  if ($Options{cellspacing} < 0) {
    die "Error: The value specified, $Options{cellspacing},  for option \"--cellspacing\" is not valid. Allowed values: >= 0 \n";
  }
}

__END__

=head1 NAME

SDFilesToHTML.pl - Generate HTML table file(s) from SDFile(s)

=head1 SYNOPSIS

SDFilesToHTML.pl  SDFiles(s)...

SDFilesToHTML.pl [B<-a, --align> left | center | right,[top | middle | bottom]] [B<-b, --border> borderwidth] [B<--cellpadding> padding]
[B<--cellspacing> spacing] [B<--cmpddatafield> "fieldlabel,[label,position,alignment]"] [B<--datafields> "fieldlabel,[fieldlabel]..." | Common | All]
[B<--footer> string] [B<-d, --displaylinks> top | bottom | both] [B<--displaylinksinfo> compound | table | both] [B<-h, --help>]
[B<--headeralign> left | center | right,[top | middle | bottom]] [B<--headercolor> "#RRGGBB"]
[B<--highlight> "fieldlabel,datatype,criterion,value,[fieldlabel,datatype,criterion,value,...]"]
[B<--highlightcolor> "#RRGGBB,#RRGGBB"] [B<--highlightstyle> text | background]
[B<-m, --mode> plain | shade | highlight | shadedhighlight | structuresonly | shadedstructuresonly]
[B<-n, --numcmpds> number] [B<-o, --overwrite>] [B<-r, --root> rootname] [B<-s, --structure> display | link]
[B<--strlinkmode> plain | shaded] [B<--strlinknavigation> yes | no]
[B<--strlinkshadecolor> "#RRGGBB"] [B<--strlinktitle> string] [B<--strlinktitledisplay> yes | no] [B<--strlinktype> href | button]
[B<--strviewertype> Chem3DActiveX | ChemDrawActiveX | ChemDrawPlugIn | Chime | JME | Jmol | MarvinView | ViewerActiveX]
[B<--strviewerconfig> codebase[,archive,code]] [B<--strviewerparams> "name=value [name=value ...]"]
[B<--strviewerembed> direct | javascript] [B<--strviewerjsfile> javascriptfilename]
[B<--strtablesize> "numrows,numcols"] [B<--stylesheet> old | new | none]
[B<--stylesheetname> filename] [B<--shadecolor> "#RRGGBB,#RRGGBB"] [B<-t, --title> string] [B<--titledisplay> yes | no]
[B<-w, --workingdir> dirname] SDFiles(s)...

=head1 DESCRIPTION

Generate HTML file(s) from I<SDFile(s)>. The HTML file(s) contain data tables
and appropriate navigational links to view other tables; navigational links are also
provided on compound HTML pages. These files can be generated for local viewing or
deployment on a web server. A variety of options are provided to control style and
appearance of tables. And for viewing structures, options are available to use any one of
these viewers: Chem3DActiveX, ChemDrawActiveX, ChemDrawPlugIn, Chime, Jmol, JME,
MarvinView, or ViewerActiveX. Jmol is the default structure viewer and it is also distributed
along with this package; however, to use any other supported viewers, make sure it's available
in your environment.

Multiple I<SDFile(s)> names are separated by space. The valid file extensions are
I<.sdf> and I<.sd>. All other file names are ignored. All the SD files in a current directory can
be specified either by I<*.sdf> or the current directory name.

=head1 OPTIONS

=over 4

=item B<-a, --align> I<left | center | right,[top | middle | bottom]>

Horizontal and vertical alignment for table rows except for header row which is specified
using B<--headeralign> option. Possible horizontal alignment values: I<left, center, or right>.
Possible vertical alignment values: I<top, middle, or bottom>.

Default values: I<left,middle>

=item B<-b, --border> I<borderwidth>

Table border width. Default value: 1 for I<plain> and I<highlight> mode; 0 for I<shade>
and I<shadedhightlight> mode. Zero indicates no border.

=item B<--cellpadding> I<padding>

Table cell padding. Default value: I<2>.

=item B<--cellspacing> I<spacing>

Table cell spacing. Default value: I<1>.

=item B<--cmpddatafield> I<fieldlabel,[label,position,alignment]>

This value is mode specific. It indicates data field value to be displayed with the structure along
with its label, position and alignment during  I<structuresonly | shadedstructuresonly> value of B<-m, --mode>
option. Possible values: feldlabel - valid data field label; label - yes or no; position - I<top or bottom>; alignment
- I<left, center, or right>. Default: I<none,no,bottom,center>. Example:

    MolWt,no,bottom,middle

B<--cmpddatafield> option value is also linked to compound summary page.

=item B<--datafields> I<"fieldlabel,[fieldlabel]..." | Common | All>

Data fields to display in HTML table(s). Possible values: list of comma separated data field
labels, data fields common to all records, or all data fields. Default value: I<All>.
Examples:

    ALogP,MolWeight,EC50
    "MolWeight,PSA"

=item B<--footer> I<string>

Text string to be included at bottom of each HTML file. Default: none.

=item B<-d --displaylinks> I<top | bottom | both>

Specify where to display navigation links in each HTML file for accessing all other HTML
files. Possible values: I<top, bottom, or both>. Default: I<both>. This option is
only valid during multiple HTML files generation for an input file.

=item B<--displaylinksinfo> I<compound | table | both>

Control display of additional information along with navigational links: Showing compound
n of m is displyed for compound and showing table n of m for table. Possible values: I<compound
| table | both>. Default: I<both>. This option is only valid  during multiple HTML files generation.

=item B<-h, --help>

Print this help message.

=item B<--headeralign> I<left | center | right,[top | middle | bottom>

Horizontal and vertical alignment for table header rows. Possible horizontal alignment
values: I<left, center, or right>. Possible vertical alignment values: I<top, middle, or bottom>.

Default values: I<center,middle>

=item B<--headercolor> I<#RRGGBB>

Color used to fill background of table header row containing column labels
represented as a hexadecimal string. Default value: None for B<-m, --mode> option
value of I<plain> and I<#ccccff>, light blue, for others.

=item B<--highlight> I<"fieldlabel,datatype,criterion,value,[fieldlabel,datatype,criterion,value,...]">

Highlighting methodology used to highlight various SDFile(s) data field values in
HTML file(s). Same set of quartets values are applied to all SDFile(s).

Input text contains these quartets: I<fieldlabel,datatype,criterion,value,...>.
Possible datatype values: I<numeric or text>. Possible criterion values: I<le, ge, or eq>.
Examples:

    "MolWt,numeric,le,450"
    "MolWt,numeric,le,450,LogP,numeric,le,5"
    Name,text,eq,Aspirin

=item B<--highlightcolor> I<"#RRGGBB,#RRGGBB">

Colors used to highlight column values during I<highlight> and I<shadedhightlight>
mode represented as hexadecimal strings.

For B<--highlighstyle> option values of I<text> and I<background>, these colors represent
text or background colors respectively. For a specific column, first color string is used for
values which meet criterion indicated by B<--highlight> option; the second color is used for
 rest of the values.

Default values for I<background> B<--highlightstyle>: I<"#0fff0f,#ff0f0f">. And default values for
I<text> B<--highlightstyle>: I<"#0fbb0f,#ff0f0f">. Hexadecimal strings for both B<--highlightstyle>
colors correspond to I<reddish> and I<greenish>.

=item B<--highlightstyle> I<text | background>

This value is mode specific. It indicates highlight style used to differentiate column
values which meet a specified criterion in B<--highlight> option. Possible values: I<text or
background>. Default: I<background>.

=item B<-m, --mode> I<plain | shade | highlight | shadedhighlight | structuresonly | shadedstructuresonly>

Specify how to generate HTML table(s): plain tables with line borders, background of
alternate rows filled with a specified color, column values highlighted using a specified
criteria, combination of previous two styles,  tables containing only structures, or tables
containing only structures with filled background of alternate rows.

Possible values: I<plain, shade, highlight, shadedhighlight, structuresonly, or
shadedstructuresonly>. Default: I<shade>.

=item B<-n, --numcmpds> I<number>

Maximum number of compounds per table. Default value: I<15> for tables with structures and
I<50> for tables with links to structures. Use 0 to put all compounds into one table. For SDFile(s)
with more than maximum number of specified compounds, multiple HTML tables, with appropriate
navigation links, are created.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

New file or directory name is generated using the root: <root>.html or <root>-html.
Default new file name: <InitialSDFileName>.html. Default directory name:
<InitialSDFileName>-html.

For SDFile(s) with more than maximum number of specified compounds per table,
this directory tree is generated using <Name> where <Name> corresponds to <root>
or <InitialSDFileName>: Top dir - <Name>-html; Sub dirs - html and mols. <Top dir> contains
<Name>.html and <Name>.css files and <sub dir> html conatins various
<Name>Lines<Start>To<End>.html files; <sub dir> mols is created as needed and contains
MOL files.

This option is ignored for multiple input files.

=item B<-s, --structure> I<display | link>

Structure display control: display structures in a table column or set up a link for each
structure which opens up a new HTML page containing structure and other appropriate
information. Possible values: I<display or link>. Default value: I<display>

=item B<--strlinkmode> I<plain | shaded>

Specify how to display compound HTML page: plain or background of data field
field labels is filled with a specified color. Possible values: I<plain or shad>.
Default value: I<plane>.

Structure viewer background color is white. Use B<--strviewerparams> option to change
default behavior of structure viewers.

=item B<--strlinknavigation> I<yes | no>

Display navigation links to other compounds in compound HTML page. Possible values:
I<yes or no>. Default value: I<yes>.

=item B<--strlinkshadecolor> I<"#RRGGBB">

This value is B<--strlinkmode> specific. For I<shade> value of B<--strlinkmode> option, it
represents colors used to fill background of data field labels.

Default value: I<"#e0e9eb"> - it's a very light blue color.

=item B<--strlinktitle> I<string>

Title for compound HTML page. Default value: I<Compound Summary>.

=item B<--strlinktitledisplay> I<yes | no>

Display title for compound HTML page. Possible values: I<yes or no>. Default value: I<no>.

=item B<--strlinktype> I<href | button>

Type of structure link. Possible values: I<href or button>. Default: I<href>.

=item B<--strviewertype> I<Chem3DActiveX | ChemDrawActiveX | ChemDrawPlugIn | Chime | JME | Jmol | MarvinView | ViewerActiveX>

Structure viewer supported for viewing structures. Possible values: I<Chem3DActiveX,
ChemDrawActiveX, ChemDrawPlugIn, Chime, JME, Jmol, MarvinView, or ViewerActiveX>.
Default value: I<Jmol>.

Assuming you have access to one of these viewers on your machine, you are all set
to use this script. Otherwise, visit one of these web sites to download and install
your favorite viewer:

    accelrys.com: Viewer ActiveX 5.0
    cambridgesoft.com: Chem3DActiveX 8.0, ChemDrawActiveX 8.0,
                       ChemDrawPlugIn
    chemaxon.com: MarvinView applet
    mdli.com: Chime plug-in
    jmol.sourceforge.net: JmolApplet V10
    molinspiration.com: JME applet

The default viewer, JmolApplet V10, is distributed with MayaChemTools package.
Earlier versions of JmolApplet are not supported: due to applet security issues related to
reading files, this script uses in-line loading of MOL files and this option doesn't exist in
earlier version of JmolApplet.

=item B<--strviewerconfig> I<codebase[,archive,code]>

Configuration information for structure viewers. This option is only valid for structure
viewers which are applets: Jmol, JME and MarvinView. For other viewer types available via
B<--strviewertype> option  - MDL Chime, ChemDrawActiveX, ChemDrawPlugIn, and
Chem3DActiveX - this value is ignored.

Input text format: I<codebase[,archive,code]>. For an applet viewer, I<codebase> must be
specified; I<archive> and I<code> values are optional. Here are default I<archive> and
I<codebase> values for various applets: Jmol - JmolApplet, JmolApplet.jar; JME - JME, JME.jar;
 MarvinView: MView, marvin.jar

For local deployment of HTML files, I<codebase> must correspond to a complete path to
the local directory containing appropriate I<archive> file and the complete path is converted
into appropriate relative path during generation of HTML files.

By default, I<codebase> value of <this script dir>/../lib/Jmol is used for I<Jmol> applet viewer, and
HTML file(s) are generated for local deployment; however, you can specify any supported
applet viewer and generate HTML file(s) for deploying on a web server.

For deploying the HTML file(s) on a web server, specify a valid I<codebase> directory name
relative to <WWWRootDir>. Example when JME archive file, JME.jar, is available in
I</jme> directory on the web server:

    /jme

For local deployment of HTML file(s), specify a complete I<codebase> directory name.
Example when JmolApplet archive file, JmolApplet.jar, is present in <JMOLROOT> directory:

    <JMOLROOT>

In addition to I<codebase>, you can also specify I<archive> file name. Example for web
deployment:

    "/jme,JME.jar"
    "/jme"

Example for local deployment:

    "<JMEROOT>,JME.jar"
    "<JMEROOT>"

=item B<--strviewerparams> I<"name=value [name=value ...]">

Parameters name and value pairs for structure viewers. These name and value pairs
are used to control the appearance and behavior of structure viewers in tables and
compound HTML page during I<link> value for B<-s --structure> option.

The parameter names, along with their values,  are just passed to each structure viewer
in appropriate format without checking their validity. Check documentation of appropriate
structure viewers to figure out valid parameter names.

Input text format: I<name=value name=value ...> Example:

    "width=250 height=170"

Default for all structure viewers: I<width=250 height=170> for displaying structures in
tables, and I<strlinkwidth=500 strlinkheight=295> for compound HTML page during I<link> value
for B<-s --structure> option.

Default background color for all structure viewers: same as B<--shadecolor> value for
displaying structures in tables and I<strlinkbgcolor=#ffffff> for  compound HTML page;
however, explicit specification of background color in this option overrides default value.
To use black background for structures in tables and compound HTML page, specify I<bgcolor=#000000>
and I<strlinkbgcolor=#000000> respectively.  Keep this in mind: Some structure viewers
don't appear to support background color parameter.

Additional structure viewer specific default values:

    Chem3DActiveX: "displaytype=Ball&Stick rotationbars=false
                    moviecontroller=false"
    ChemDrawActiveX: "ViewOnly=1 ShrinkToFit=1 ShowToolsWhenVisible=1"
    ChemDrawPlugIn: "type=chemical/x-mdl-molfile ViewOnly=1
                     ShrinkToFit=1 ShowToolsWhenVisible=1"
    Chime: "display2d=true"
    JME: "options=depict"
    Jmol: "progressbar=true progresscolor=#0000ff boxbgcolor=#000000
           boxfgcolor=#ffffff script="select *; set frank off;
           wireframe on; spacefill off""
    MarvinView: "navmode=zoom"
    ViewerActiveX:"Mouse=4 Convert2Dto3D=0"

Try overriding default values or specify additional valid parameter/value pairs to get desired
results. Example for using CPK rendering scheme with Jmol viewer:

    "script="select *; set frank off; wireframe off; spacefill on""

=item B<--strviewerembed> I<direct | javascript>

Specify how to embed structure viewers in HTML pages. Possible values: I<direct> - use applet/object
tags to emded structure viewer; I<javascript> - use vendor supplied java scripts. Default value:
direct.

This option only applies to these vieweres: I<Chem3DActiveX, ChemDrawActiveX, ChemDrawPlugIn,
Jmol, and MarvinView>.

For marvin.js to work correctly on your browser, you may need to set I<marvin_jvm=builtin> or
I<marvin_jvm=plugin> using B<--strviewerparams> option. Additionally, MarvinView - at least
in my hands - also has problems during usage of JavaScript for local deployment; however, it
does work via web server.

As far as I can tell, Jmol.js supplied with Jmol10 release has these issues: jmolSetAppletColor
doesn't support background color; jmolInitialize disables relative specification of codebase
directroy which works okay. So, use Jmol.js supplied with MayaChemTools.

=item B<--strviewerjsfile> I<java script file name>

Name of vendor supplied java script file. Default values: Chem3DActiveX: I<chem3d.js>; ChemDrawActiveX,
and ChemDrawPlugIn: I<chemdraw.js>; Jmol: I<Jmol.js>, MarvinView: I<marvin.js>.

Directory location for these files is specified via I<codebase> value of B<--strviewerconfig> option.

=item B<--strtablesize> I<"numrows,numcols">

This option is only valid for I<structuresonly> and I<shadedstructuresonly> modes. And it indicates
maximum number of rows and columns per structure table. Default value:I<6,4>.

=item B<--stylesheet> I<old | new | none>

Controls usage of stylesheet for newly generated HTML file(s). Possible values: I<old,
new, or none>. Default value: I<new>.

Stylesheet file contains various properties which control appearance of HTML pages:
type, size, and color of fonts; background color; and so on.

For I<old> value, an existing stylesheet file specified by B<--stylesheetname> option is
used for each HTML file; no new stylesheet file is created. This option is quite handy
for deploying HTML file(s) on a web server: assuming you specify a valid stylesheet
file location relative to your WWWRoot, a reference to this stylesheet is added to each
HTML file. For local deployment of HTML file(s), a complete path to a local stylesheet
is fine as well.

For I<create> value, a new stylesheet is created and reference to this local stylesheet
is added to each HTML file. Use option B<--stylesheetname> to specify name.

For I<none> value, stylesheet usage is completely ignored.

=item B<--stylesheetname> I<filename>

Stylesheet file name to be used in conjunction with B<-s --stylesheet> option. It is only
valid for I<old> value of B<-s --stylesheet> option. Specify a valid stylesheet file location
relative to your WWWRoot and a reference to this stylesheet is added to each HTML
file. Example: I<"/stylesheets/MyStyleSheet.css">. Or a complete path name to a local
stylesheet file.

For I<create> value of B<-s --stylesheet> option, a new stylesheet file is created using
B<-r --root> option. And value of B<--stylesheetname> is simply ignored.

=item B<--shadecolor> I<"#RRGGBB,#RRGGBB">

Colors used to fill background of rows during I<shade> and I<shadedhightlight> mode
represented as a pair of hexadecimal string; the first and second color values
are used for odd and even number rows respectively.

Default value: I<"#ffffff,#e0e9eb"> - it's white and very light blue for odd and even number rows.

=item B<-t, --title> I<string>

Title for HTML table(s). Default value: I<SDFileName>. This option is ignored for
multiple input files. And B<-r --root> option is used to generate appropriate
titles.

=item B<--titledisplay> I<yes | no>

Display title for HTML table(s). Possible values: I<yes or no>. Default value: I<yes>.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

HTML table file(s), containing structures, can be used in two different ways: browsing on a
local machine or deployment via a web server. By default, HTML file(s) are created for viewing
on a local machine using Jmol viewer through a browser; however, you can specify any
supported applet viewer and  generate HTML file(s) for deploying on a web server.

First two sets of examples show generation of HTML file(s) using different applet viewers
and a variety of options for local browsing; last set deals with web deployment.

B<Local deployment: Usage of default JMol viewer distributed with MayaChemTools:>

To generate HTML tables with structure display using JMol viewer, rows background filled
with white and light blue colors, navigation links on top and botton of each page, type:

    % SDFilesToHTML.pl -o Sample1.sdf

To generate HTML tables with structure display using JMol viewer, rows background filled
with white and light blue colors, navigation links on top and botton of each page, and
only containing MolWeight and Mol_ID SD data fields, type:

    % SDFilesToHTML.pl --datafields "MolWeight,Mol_ID" -o Sample1.sdf

To generate HTML tables with CPK structure display using JMol viewer, rows
background filled with white and light blue colors, navigation links on top and botton of
each page, type:

    % SDFilesToHTML.pl --strviewerparams "script=\"select *; set frank off;
      wireframe off; spacefill on\"" -o Sample1.sdf

To generate HTML tables with structure display using JMol viewer and black background, rows
background filled with light golden and greyish colors, navigation links on top and botton of
each page, 10 rows in each table, greyish header row color, and cell spacing of 1, type:

    % SDFilesToHTML.pl -o -n 10 --headeralign "center" --headercolor
      "#a1a1a1" --shadecolor "#fafad2,#d1d1d1" --cellspacing 1
      --strviewerparams "bgcolor=#000000" Sample1.sdf

To highlight molecular weight values using specified highlight criteria and fill in default background
colors, type:

    % SDFilesToHTML.pl -n 10 --highlight "MolWeight,numeric,le,450"
      --highlightstyle background -m shadedhighlight -o Sample1.sdf

To highlight molecular weight values using specified highlight criteria, color the text using
default colors, and add a footer message in every page, type:

    % SDFilesToHTML.pl -n 4 --highlight "MolWeight,numeric,le,500"
      --highlightstyle text -m shadedhighlight -o
      --footer "Copyright (C) MayaChemTools" --cellspacing 1 Sample1.sdf

To generate tables containing only structures, type:

    % SDFilesToHTML.pl -d both -m shadedstructuresonly --strtablesize "6,4"
      --cellspacing 1 -b 1 -o Sample1.sdf

To generate tables containing only structures with molecular weight displayed above the
structure, type:

    % SDFilesToHTML.pl -d both -m shadedstructuresonly --strtablesize "6,4"
      --cmpddatafield "MolWeight,no,top,center"  --cellspacing 1 -b 1
      -o Sample1.sdf

To generate tables containing links to structures and highlight molecular weight data field values
using specified highlight criteria , type:

    % SDFilesToHTML.pl -n 4 --footer "Copyright (C) MayaChemTools"
      --highlight "MolWeight,numeric,le,450" --highlightstyle background
      -d both -m shadedhighlight  -s link --strlinktype button
      -o Sample1.sdf

B<Local deployment: Usage of other structure viewers:>

    % SDFilesToHTML.pl --strviewertype MarvinView --strviewerconfig
      "<Marvin dir path>" -o Sample1.sdf

    % SDFilesToHTML.pl -o -n 10 --headeralign "center" --headercolor
      "#a1a1a1" --shadecolor "#fafad2,#d1d1d1" --cellspacing 1
      --strviewerparams "bgcolor=#000000" --strviewertype Chime
      Sample1.sdf

    % SDFilesToHTML.pl -n 10 --highlight "MolWeight,numeric,le,450"
      --highlightstyle background -m shadedhighlight --strviewertype
      Chime -o Sample1.sdf

    % SDFilesToHTML.pl -d both -m shadedstructuresonly --strtablesize "6,4"
      --cellspacing 1 -b 1 -strviewertype JME -strviewerconfig "<JME dir
      path>" -o Sample1.sdf

B<Web deployment: Usage of different structure viewers and options:>

For deploying HTML file(s) on a web server, specify a valid I<codebase> directory name
relative to <WWWRootDir>. In addition to I<codebase>, you can also specify I<archive> file
name.

    % SDFilesToHTML.pl -m plain -s display --strviewertype Jmol
      -strviewerconfig "/jmol" -n 5 -d both -r PlainTable -t "Example
      using Jmol: Plain Table" -o Sample1.sdf

    % SDFilesToHTML.pl -n 5 -m shade  -s display -strviewertype JME
      -strviewerconfig "/jme,JME.jar" -r ShadeTable -t "Example using JME:
      Shaded Table" -o Sample.sdf

    % SDFilesToHTML.pl -n 5 --highlight "MolWeight,numeric,le,450"
      --highlightstyle background  -d both -m shadedhighlight  -s display
      -strviewertype MarvinView -strviewerconfig "/marvin" -r
      ShadedHightlightTable -t "Example using MarvinView: Shaded and
      Highlighted Table" -o Sample.sdf

    % SDFilesToHTML.pl -n 4 --highlight "MolWeight,numeric,le,450" -s link
      --strlinktype href --strviewertype ChemDrawPlugIn  --highlightstyle
      background -m shadedhighlight -t "Example using ChemDrawPlugIn:
      Shaded and Highlighted Table" -r ShadedHightlightTable -o Sample1.sdf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

FilterSDFiles.pl, InfoSDFiles.pl, SplitSDFiles.pl, MergeTextFilesWithSD.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
