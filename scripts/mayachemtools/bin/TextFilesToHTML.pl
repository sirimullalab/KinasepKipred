#!/usr/bin/perl -w
#
# File: TextFilesToHTML.pl
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

my(@TextFilesList);
@TextFilesList = ExpandFileNames(\@ARGV, "csv tsv");

print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

print "Checking input text file(s)...\n";
my(%TextFilesInfo);
RetrieveTextFilesInfo();
SetupCoulmnsTablesAndMiscInfo();

# Generate output files...
my($FileIndex);
if (@TextFilesList > 1) {
  print "\nProcessing text files...\n";
}
for $FileIndex (0 .. $#TextFilesList) {
  if ($TextFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $TextFilesList[$FileIndex]...\n";
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

  if ($TextFilesInfo{MultipleHTMLTables}[$Index]) {
    GenerateMultipleHTMLTable($Index);
  }
  else {
    GenerateOneHTMLTable($Index);
  }
}

# Generate one table...
sub GenerateOneHTMLTable {
  my($Index) = @_;
  my($TextFile, $TopHTMLDir, $HTMLFile, $Line, $StartRowNum, $EndRowNum, $CSSFile, $CSSFilePath, $CSSRef);

  $HTMLFile = $TextFilesInfo{HTMLRoot}[$Index] . ".html";
  $TextFile = $TextFilesList[$Index];

  # Setup data directories...
  ($TopHTMLDir) = SetupDataDirs($Index);

  # Setup stylesheet file...
  $CSSRef = "";
  if ($Options{stylesheet} =~ /^new$/i) {
    $CSSFile = $TextFilesInfo{HTMLRoot}[$Index] . ".css"; $CSSRef = ".\/" . "$CSSFile";
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
  open TEXTFILE, "$TextFile" or die "Error: Can't open $TextFile: $! \n";

  # Write out HTML page header...
  print HTMLFILE SetupHTMLPageHeader($TextFilesInfo{HTMLTitle}[$Index], $CSSRef);
  if ($OptionsInfo{TitleDisplay}) {
    print HTMLFILE SetupHTMLPageTitle($TextFilesInfo{HTMLTitle}[$Index]);
  }
  else {
    print HTMLFILE SetupHTMLEmptyLines(1);
  }

  # Start the table...
  print HTMLFILE SetupHTMLAlignmentBegin("center");
  print HTMLFILE SetupHTMLTableHeader($OptionsInfo{TableBorder}, $OptionsInfo{TableCellPadding}, $OptionsInfo{TableCellSpacing});

  WriteColLabels($Index, \*TEXTFILE, \*HTMLFILE);

  # Skip the labels and write out all the other rows...
  $Line = <TEXTFILE>;
  $StartRowNum = 1;
  $EndRowNum = $TextFilesInfo{LineCount}[$Index];
  WriteRowValues($Index, $StartRowNum, $EndRowNum, \*TEXTFILE, \*HTMLFILE);

  # Finish up the table...
  print HTMLFILE SetupHTMLTableEnd();
  print HTMLFILE SetupHTMLAlignmentEnd("center");

  # Write out HTML page end...
  print HTMLFILE SetupHTMLPageEnd($OptionsInfo{Footer});

  close HTMLFILE;
  close TEXTFILE;
}

# Generate multiple tables...
sub GenerateMultipleHTMLTable {
  my($Index) = @_;
  my($TopHTMLDir, $SubHTMLDir, $TextFile, $HTMLFile, $TableNum, $TableCount, $TableIndex, $TableStartLineNum, $TableEndLineNum, $Line, $InSubHTMLDir, $PrintMsg, $CSSFile, $CSSFilePath, $CSSRef, $NewStyleSheet);

  # Open text file and skip over label line...
  $TextFile = $TextFilesList[$Index];
  open TEXTFILE, "$TextFile" or die "Error: Can't open $TextFile: $! \n";
  $Line = <TEXTFILE>;

  # Set up data directories to hold various html files...
  ($TopHTMLDir, $SubHTMLDir) = SetupDataDirs($Index);

  # Create stylesheet file...
  $CSSRef = "";
  $NewStyleSheet = 0;
  if ($Options{stylesheet} =~ /^new$/i) {
    $NewStyleSheet = 1;
    $CSSFile = $TextFilesInfo{HTMLRoot}[$Index] . ".css";
    $CSSFilePath = "$TopHTMLDir" . "\/" . $CSSFile;
    GenerateStyleSheetFile($CSSFilePath);
  }
  elsif ($Options{stylesheet} =~ /^old$/i) {
    $CSSRef = $Options{stylesheetname};
  }

  $PrintMsg = 1;
  # Generate HTML files for all the tables...
  $TableCount = $TextFilesInfo{TableCount}[$Index];
  for $TableNum (1 .. $TableCount) {
    $TableIndex = $TableNum - 1;
    $HTMLFile = ${$TextFilesInfo{TableHTMLFiles}[$Index]}[$TableIndex];
    $TableStartLineNum = ${$TextFilesInfo{TableStartLineNum}[$Index]}[$TableIndex];
    $TableEndLineNum = ${$TextFilesInfo{TableEndLineNum}[$Index]}[$TableIndex];

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
	  print "Generating ", ($TableCount - 1), " other HTML files: $SubHTMLDir\/$TextFilesInfo{HTMLRoot}[$Index]\*.html...\n";
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
    print HTMLFILE SetupHTMLPageHeader($TextFilesInfo{HTMLTitle}[$Index], $CSSRef);

    # Set up the navigation links for this table...
    if ($OptionsInfo{NavLinksAtTop}) {
      WriteNavigationLinks($Index, $TableNum, \*HTMLFILE);
    }
    # Setup page title...
    if ($OptionsInfo{TitleDisplay}) {
      print HTMLFILE SetupHTMLPageTitle($TextFilesInfo{HTMLTitle}[$Index]);
    }
    else {
      print HTMLFILE SetupHTMLEmptyLines(1);
    }

    # Start the table...
    print HTMLFILE SetupHTMLAlignmentBegin("center");
    print HTMLFILE SetupHTMLTableHeader($OptionsInfo{TableBorder}, $OptionsInfo{TableCellPadding}, $OptionsInfo{TableCellSpacing});

    WriteColLabels($Index, \*TEXTFILE, \*HTMLFILE);

    # Write out appropriate row data for this table...
    WriteRowValues($Index, $TableStartLineNum, $TableEndLineNum, \*TEXTFILE, \*HTMLFILE);

    # Finish up the table...
    print HTMLFILE SetupHTMLTableEnd();
    print HTMLFILE SetupHTMLAlignmentEnd("center");

    # Set up the navigation links for this table...
    if ($OptionsInfo{NavLinksAtBottom}) {
      print HTMLFILE SetupHTMLEmptyLines(1);
      WriteNavigationLinks($Index, $TableNum, \*HTMLFILE);
    }

    # Write out HTML page end...
    print HTMLFILE SetupHTMLPageEnd($OptionsInfo{Footer});
    close HTMLFILE;
  }
  close TEXTFILE;

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
  my($Index, $TextFileRef, $HTMLFileRef) = @_;
  my(@ColLabels, $Label);

  print $HTMLFileRef $TextFilesInfo{TableRowHeaderTags};

  @ColLabels = @{$TextFilesInfo{ColLabels}[$Index]};
  for $Label (@ColLabels) {
    print $HTMLFileRef SetupHTMLTableRowHeaderValue($Label);
  }
  print $HTMLFileRef $TextFilesInfo{RowEndTags};
}

#Write out the rows value...
sub WriteRowValues {
  my($Index, $StartRowNum, $EndRowNum, $TextFileRef, $HTMLFileRef) = @_;
  my($ColNum, $BackgroundColor, $FontColor, $LineCount, $Line, @RowValues, $Value, $InDelim, $LastColNum);

  $InDelim = $TextFilesInfo{InDelim}[$Index];
  $LastColNum = @{$TextFilesInfo{ColLabels}[$Index]} - 1;

  for $LineCount ($StartRowNum .. $EndRowNum) {
    $Line = GetTextLine($TextFileRef);

    if ($OptionsInfo{ShadeRowsStatus}) {
      print $HTMLFileRef ($LineCount % 2) ? $TextFilesInfo{BgFilledOddRowHeaderTags} : $TextFilesInfo{BgFilledEvenRowHeaderTags};
    }
    else {
      print $HTMLFileRef $TextFilesInfo{RowHeaderTags};
    }
    @RowValues = quotewords($InDelim, 0, $Line);
    for $ColNum (0 .. $LastColNum) {
      $Value = ($ColNum <= $#RowValues) ? $RowValues[$ColNum] : "";
      $BackgroundColor = ""; $FontColor = "";
      if ($OptionsInfo{HighlightStatus}) {
	if (exists($TextFilesInfo{HightlightColNumMap}[$Index]{$ColNum})) {
	  ($BackgroundColor, $FontColor) = GetValueHighlightColors($Index, $ColNum, $Value);
	}
      }
      print $HTMLFileRef SetupHTMLTableRowDataValue($Value, $BackgroundColor, $FontColor);
    }
    print $HTMLFileRef $TextFilesInfo{RowEndTags};
  }
}

# Setup navigation link information for each table.
#
# All table sets besides first and last have these links: FirstTable, Previous, Current-1,Current,Current+1,  Next, and LastTable
# First set: Current, Next, and LastTable
# Last set: FirstTable, Previous and Current.
#
sub WriteNavigationLinks {
  my($Index, $CurTableNum, $HTMLFileRef) = @_;
  my($TableNum, $StartTableNum, $EndTableNum, $TableIndex, $BorderWidth, $CellPadding, $CellSpacing,$HTMLFile, $HTMLRefFile, $RelativeFileDir, $HTMLRefValue, $FirstTableNum, $FirstTableIndex, $LastTableNum, $LastTableIndex, $TableStartLineNum, $TableEndLineNum, $LastLineNum, $BGColor, $LinksOffSet);

  $LinksOffSet = 10;

  $FirstTableNum = 1; $FirstTableIndex = $FirstTableNum - 1;
  $LastTableNum = $TextFilesInfo{TableCount}[$Index]; $LastTableIndex = $LastTableNum - 1;
  $LastLineNum = ${$TextFilesInfo{TableEndLineNum}[$Index]}[$LastTableIndex];

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
  print $HTMLFileRef $TextFilesInfo{RowHeaderTags};

  if ($OptionsInfo{NavLinksTableInfo} && $OptionsInfo{NavLinksLineInfo}) {
    print $HTMLFileRef SetupHTMLTableRowDataValue("Showing table $CurTableNum of $LastTableNum");
    print $HTMLFileRef SetupHTMLTableRowDataValue("&nbsp");
    print $HTMLFileRef SetupHTMLTableRowDataValue("&nbsp");
  }

  print $HTMLFileRef SetupHTMLTableRowDataValue("Tables: ");
  # Setup a link to first table...
  if ($StartTableNum != $FirstTableNum) {
    $HTMLFile = ${$TextFilesInfo{TableHTMLFiles}[$Index]}[$FirstTableIndex];
    $HTMLRefFile = GetRelativeFileDir($CurTableNum, $FirstTableNum, $FirstTableNum) . $HTMLFile;
    $TableStartLineNum = ${$TextFilesInfo{TableStartLineNum}[$Index]}[$FirstTableIndex];
    $TableEndLineNum = ${$TextFilesInfo{TableEndLineNum}[$Index]}[$FirstTableIndex];
    $HTMLRefValue = SetupHTMLHRef("First", $HTMLRefFile, "First Table Containing Lines $TableStartLineNum To $TableEndLineNum");
    print $HTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue, $LinkBGColor, $LinkTextColor, $LinkFontBold);
  }

  # Setup link to previous table...
  if ($CurTableNum != $FirstTableNum) {
    my($PreviousTableNum, $PreviousTableIndex);
    $PreviousTableNum = $CurTableNum - 1; $PreviousTableIndex = $PreviousTableNum - 1;
    $HTMLFile = ${$TextFilesInfo{TableHTMLFiles}[$Index]}[$PreviousTableIndex];
    $HTMLRefFile = GetRelativeFileDir($CurTableNum, $PreviousTableNum, $FirstTableNum) . $HTMLFile;
    $TableStartLineNum = ${$TextFilesInfo{TableStartLineNum}[$Index]}[$PreviousTableIndex];
    $TableEndLineNum = ${$TextFilesInfo{TableEndLineNum}[$Index]}[$PreviousTableIndex];
    $HTMLRefValue = SetupHTMLHRef("Previous", $HTMLRefFile, "Previous Table Containing Lines $TableStartLineNum To $TableEndLineNum");
    print $HTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue, $LinkBGColor, $LinkTextColor, $LinkFontBold);
  }

  for $TableNum ($StartTableNum .. $EndTableNum) {
    $TableIndex = $TableNum - 1;
    $HTMLFile = ${$TextFilesInfo{TableHTMLFiles}[$Index]}[$TableIndex];
    if ($TableNum == $CurTableNum) {
      print $HTMLFileRef SetupHTMLTableRowDataValue($TableNum, $LinkBGColor, $InactiveLinkNumColor, $InactiveLinkFontBold);
    }
    else {
      # Setup the link...
      my($RefTitle);
      $TableStartLineNum = ${$TextFilesInfo{TableStartLineNum}[$Index]}[$TableIndex];
      $TableEndLineNum = ${$TextFilesInfo{TableEndLineNum}[$Index]}[$TableIndex];
      $RefTitle = AddNumberSuffix($TableNum) . " Table Containing Lines $TableStartLineNum To $TableEndLineNum";
      $HTMLRefFile = GetRelativeFileDir($CurTableNum, $TableNum, $FirstTableNum) . $HTMLFile;
      $HTMLRefValue = SetupHTMLHRef($TableNum, $HTMLRefFile, $RefTitle);
      print $HTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue);
    }
  }

  # Setup link to next table...
  if ($CurTableNum != $LastTableNum) {
    my($NextTableNum, $NextTableIndex);
    $NextTableNum = $CurTableNum + 1; $NextTableIndex = $NextTableNum - 1;
    $HTMLFile = ${$TextFilesInfo{TableHTMLFiles}[$Index]}[$NextTableIndex];
    $HTMLRefFile = GetRelativeFileDir($CurTableNum, $NextTableNum, $FirstTableNum) . $HTMLFile;
    $TableStartLineNum = ${$TextFilesInfo{TableStartLineNum}[$Index]}[$NextTableIndex];
    $TableEndLineNum = ${$TextFilesInfo{TableEndLineNum}[$Index]}[$NextTableIndex];
    $HTMLRefValue = SetupHTMLHRef("Next", $HTMLRefFile, "Next Table Containing Lines $TableStartLineNum To $TableEndLineNum");
    print $HTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue, $LinkBGColor, $LinkTextColor, $LinkFontBold);
  }

  # Setup link to last table...
  if ($EndTableNum != $LastTableNum) {
    $HTMLFile = ${$TextFilesInfo{TableHTMLFiles}[$Index]}[$LastTableIndex];
    $HTMLRefFile = GetRelativeFileDir($CurTableNum, $LastTableNum, $FirstTableNum) . $HTMLFile;
    $TableStartLineNum = ${$TextFilesInfo{TableStartLineNum}[$Index]}[$LastTableIndex];
    $TableEndLineNum = ${$TextFilesInfo{TableEndLineNum}[$Index]}[$LastTableIndex];
    $HTMLRefValue = SetupHTMLHRef("Last", $HTMLRefFile, "Last Table Containing Lines $TableStartLineNum To $TableEndLineNum");
    print $HTMLFileRef SetupHTMLTableRowDataValue($HTMLRefValue, $LinkBGColor, $LinkTextColor, $LinkFontBold);
  }
  # Setup current table info text....
  print $HTMLFileRef SetupHTMLTableRowDataValue("&nbsp");
  print $HTMLFileRef SetupHTMLTableRowDataValue("&nbsp");
  $TableStartLineNum = ${$TextFilesInfo{TableStartLineNum}[$Index]}[$CurTableNum - 1];
  $TableEndLineNum = ${$TextFilesInfo{TableEndLineNum}[$Index]}[$CurTableNum - 1];
  if ($OptionsInfo{NavLinksLineInfo}) {
    print $HTMLFileRef SetupHTMLTableRowDataValue("Showing lines $TableStartLineNum to $TableEndLineNum of $LastLineNum");
  }
  else {
    print $HTMLFileRef SetupHTMLTableRowDataValue("Showing table $CurTableNum of $LastTableNum");
  }

  print $HTMLFileRef $TextFilesInfo{RowEndTags};
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
  my($FileIndex, $ColNum, $Value) = @_;
  my($DataType, $Criterion, $CriterionValue, $BgColor, $FontColor, $ValueOk, $Nothing);

  $BgColor = ""; $FontColor = "";
  $DataType = ${$TextFilesInfo{HightlightDataMap}[$FileIndex]{$ColNum}}[0];
  $Criterion = ${$TextFilesInfo{HightlightDataMap}[$FileIndex]{$ColNum}}[1];
  $CriterionValue = ${$TextFilesInfo{HightlightDataMap}[$FileIndex]{$ColNum}}[2];

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

# Setup columns, tables and other information...
sub SetupCoulmnsTablesAndMiscInfo {
  SetupColumnsToHighlightInfo();
  SetupMultipleTablesInfo();
  SetupHTMLTagsInfo();
}

# Setup columns to highlight information...
sub SetupColumnsToHighlightInfo {
  my($ColID, $DataType, $Criterion, $Value, $Index, $ColNum, $ColLabel, $ColIndex);

  @{$TextFilesInfo{HightlightColNumMap}} = ();
  @{$TextFilesInfo{HightlightDataMap}} = ();

  for $Index (0 .. $#TextFilesList) {
    %{$TextFilesInfo{HightlightColNumMap}[$Index]} = ();
    %{$TextFilesInfo{HightlightDataMap}[$Index]} = ();
    if ($TextFilesInfo{FileOkay}[$Index]) {
      SPECIFIEDCOLS: for $ColIndex (0 .. $#{$OptionsInfo{SpecifiedColIds}}) {
	$ColID = $OptionsInfo{SpecifiedColIds}[$ColIndex];
	$DataType = $OptionsInfo{SpecifiedColDataTypes}[$ColIndex];
	$Criterion = $OptionsInfo{SpecifiedColCriteria}[$ColIndex];
	$Value = $OptionsInfo{SpecifiedColValues}[$ColIndex];
	if (!$OptionsInfo{HighlightStatus}) {
	  next SPECIFIEDCOLS;
	}
	if ($Options{highlightby} =~ /^colnum$/i) {
	  $ColNum = $ColID;
	  if ($ColNum > 0 && $ColNum <= $TextFilesInfo{ColCount}[$Index]) {
	    $ColNum -= 1;
	  }
	  else {
	    warn "Warning: Ignoring column number, $ColID, specifed in quartet, \"$ColID,$DataType,$Criterion,$Value\", using \"--highlight\" option for $TextFilesList[$Index]: it doesn't exists \n";
	    next SPECIFIEDCOLS;
	  }
	}
	else {
	  $ColLabel = $ColID;
	  if (exists($TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel})) {
	    $ColNum = $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel};
	  } else {
	    warn "Warning: Ignoring column label, $ColID, specifed in quartet, \"$ColID,$DataType,$Criterion,$Value\", using \"--highlight\" option for $TextFilesList[$Index]: it doesn't exists \n";
	    next SPECIFIEDCOLS;
	  }
	}
	$TextFilesInfo{HightlightColNumMap}[$Index]{$ColNum} = $ColNum;
	@{$TextFilesInfo{HightlightDataMap}[$Index]{$ColNum}} =();
	push @{$TextFilesInfo{HightlightDataMap}[$Index]{$ColNum}}, ($DataType, $Criterion, $Value);
      }
    }
  }
}

# Setup navigation link information for multiple tables...
sub SetupMultipleTablesInfo {
  my($Index, $LinesPerTable);

  $LinesPerTable = $Options{numrows};
  @{$TextFilesInfo{TableCount}} = ();
  @{$TextFilesInfo{TableHTMLFiles}} = ();
  @{$TextFilesInfo{TableStartLineNum}} = ();
  @{$TextFilesInfo{TableEndLineNum}} = ();

  for $Index (0 .. $#TextFilesList) {
    $TextFilesInfo{TableCount}[$Index] = 1;
    @{$TextFilesInfo{TableHTMLFiles}[$Index]} = ();
    @{$TextFilesInfo{TableStartLineNum}[$Index]} = ();
    @{$TextFilesInfo{TableEndLineNum}[$Index]} = ();

    if ($TextFilesInfo{FileOkay}[$Index]) {
      if ($TextFilesInfo{MultipleHTMLTables}[$Index]) {
	my($TableIndex, $TotalLines, $TableCount, $TableStartLineNum, $TableEndLineNum, $Name);

	$TotalLines = $TextFilesInfo{LineCount}[$Index];
	$TableCount = ($TotalLines % $LinesPerTable) ? (int($TotalLines/$LinesPerTable) + 1) : ($TotalLines/$LinesPerTable);
	$TextFilesInfo{TableCount}[$Index] = $TableCount;
	for $TableIndex (1 .. $TableCount) {
	  $TableStartLineNum = ($TableIndex - 1) * $LinesPerTable + 1;
	  $TableEndLineNum = ($TableIndex == $TableCount) ? $TotalLines : ($TableIndex * $LinesPerTable);
	  push @{$TextFilesInfo{TableStartLineNum}[$Index]}, $TableStartLineNum;
	  push @{$TextFilesInfo{TableEndLineNum}[$Index]}, $TableEndLineNum;

	  # Setup HTML file names for all the tables...
	  $Name = "Lines" . "$TableStartLineNum" . "To" . "$TableEndLineNum";
	  if ($TableIndex == 1) {
	    $Name = "";
	  }
	  $Name = $TextFilesInfo{HTMLRoot}[$Index] . $Name . ".html";
	  push @{$TextFilesInfo{TableHTMLFiles}[$Index]}, $Name;
	}
	#print "$TextFilesList[$Index]: $TableCount -  @{$TextFilesInfo{TableStartLineNum}[$Index]} - @{$TextFilesInfo{TableEndLineNum}[$Index]} -  @{$TextFilesInfo{TableHTMLFiles}[$Index]}\n";
      }
    }
  }
}

# Setup HTML tags information...
sub SetupHTMLTagsInfo {
  # Setup row tags...
  $TextFilesInfo{RowHeaderTags} = "";
  $TextFilesInfo{RowEndTags} = "";
  $TextFilesInfo{BgFilledOddRowHeaderTags} = "";
  $TextFilesInfo{BgFilledEvenRowHeaderTags} = "";
  $TextFilesInfo{TableRowHeaderTags} = "";

  $TextFilesInfo{RowHeaderTags} = SetupHTMLTableRowHeader($OptionsInfo{RowHAlignment}, "", $OptionsInfo{RowVAlignment});
  $TextFilesInfo{RowEndTags} = SetupHTMLTableRowEnd();

  if ($OptionsInfo{ShadeRowsStatus}) {
    $TextFilesInfo{BgFilledOddRowHeaderTags} = SetupHTMLTableRowHeader($OptionsInfo{RowHAlignment}, $OptionsInfo{OddRowsShadeColor}, $OptionsInfo{RowVAlignment});
    $TextFilesInfo{BgFilledEvenRowHeaderTags} = SetupHTMLTableRowHeader($OptionsInfo{RowHAlignment}, $OptionsInfo{EvenRowsShadeColor}, $OptionsInfo{RowVAlignment});
  }

  $TextFilesInfo{TableRowHeaderTags} = SetupHTMLTableRowHeader($OptionsInfo{TableHeaderRowHAlignment}, $OptionsInfo{TableHeaderRowColor}, $OptionsInfo{TableHeaderRowVAlignment});

}

#Make sure appropriate mode specific option values are specified...
sub ProcessOptions {

  %OptionsInfo = ();

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
      die "Error: Invalid number of values, ", scalar(@AlignValues) , ", specified by \"-a --align\" option.\nIt must contain only one or two value.\n";
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

  $OptionsInfo{TitleDisplay} = ($Options{titledisplay} =~ /^yes$/i) ? 1 : 0;

  if (exists($Options{border})) {
    $OptionsInfo{TableBorder} = $Options{border};
  }
  else {
    $OptionsInfo{TableBorder} = ($Options{mode} =~ /^(plain|highlight)$/i) ? 1 : 0;
  }
  $OptionsInfo{TableCellPadding} = $Options{cellpadding};
  $OptionsInfo{TableCellSpacing} = $Options{cellspacing};
  $OptionsInfo{Footer} = $Options{footer} ? $Options{footer} : "";

  if ($Options{headercolor}) {
    $OptionsInfo{TableHeaderRowColor} = $Options{headercolor};
  }
  else {
    $OptionsInfo{TableHeaderRowColor} = ($Options{mode} =~ /^plain$/i) ? "" : "#ccccff";
  }

  $OptionsInfo{NavLinksAtBottom} = 1; $OptionsInfo{NavLinksAtTop} = 0;
  if ($Options{displaylinks} =~ /^(both|top)$/i) {
    $OptionsInfo{NavLinksAtTop} = 1;
  }
  $OptionsInfo{NavLinksTableInfo} = 1; $OptionsInfo{NavLinksLineInfo} = 0;
  if ($Options{displaylinksinfo} =~ /^both$/i) {
    $OptionsInfo{NavLinksLineInfo} = 1;
    $OptionsInfo{NavLinksTableInfo} = 1;
  }
  elsif ($Options{displaylinksinfo} =~ /^line$/i) {
    $OptionsInfo{NavLinksLineInfo} = 1;
    $OptionsInfo{NavLinksTableInfo} = 0;
  }

  if ($Options{stylesheet} =~ /^old$/i ) {
    if (!$Options{stylesheetname}) {
      die "Error: No stylesheet name specified using \"--stylesheetname\" option: It is required for \"old\" value of \"-s --stylesheet\" option. \n";
    }
  }

  my(@ColorValues);
  $OptionsInfo{OddRowsShadeColor} = ""; $OptionsInfo{EvenRowsShadeColor} = ""; $OptionsInfo{ShadeRowsStatus} = 0;
  if ($Options{mode} =~ /^(shade|shadedhighlight)$/i) {
    $OptionsInfo{OddRowsShadeColor} = "#ffffff";
    $OptionsInfo{EvenRowsShadeColor} = "#e0e0eb";
    $OptionsInfo{ShadeRowsStatus} = 1;
    if ($Options{shadecolor}) {
      # Make sure only two value are specified...
      @ColorValues = split ",", $Options{shadecolor};
      if (@ColorValues == 2) {
	$OptionsInfo{OddRowsShadeColor} = $ColorValues[0];
	$OptionsInfo{EvenRowsShadeColor} = $ColorValues[1];
      }
      else {
	die "Error: Invalid number of values, ", scalar(@ColorValues) , ", specified by \"--shadecolor\" option.\nIt must contain only two values.\n";
      }
    }
  }
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
    my($Index, $Col, $DataType, $Criterion, $Value);

    @{$OptionsInfo{SpecifiedColIds}} = ();
    @{$OptionsInfo{SpecifiedColDataTypes}} = ();
    @{$OptionsInfo{SpecifiedColCriteria}} = ();
    @{$OptionsInfo{SpecifiedColValues}} = ();
    for ($Index = 0; $Index < @HighlightValueQuartets; $Index = $Index + 4) {
      $Col = $HighlightValueQuartets[$Index];
      $DataType = $HighlightValueQuartets[$Index + 1];
      $Criterion = $HighlightValueQuartets[$Index + 2];
      $Value = $HighlightValueQuartets[$Index + 3];
      if ($Options{highlightby} =~ /^colnum$/i ) {
	if (!IsPositiveInteger($Col)) {
	  die "Error: Invalid column id, $Col, specified in quartet, \"$Col,$DataType,$Criterion,$Value\", using \"--hightlight\" option: It must be an integer value > 0 for $HighlightMode \"-m --mode\" and $HighlightBy \"--highlightby\" option values.\n";
	}
      }
      if ($DataType !~ /^(numeric|text)$/i) {
	die "Error: Invalid column data type, $DataType, specified in quartet, \"$Col,$DataType,$Criterion,$Value\", using \"--hightlight\" option: Valid values: numeric or text\n";
      }
      if ($Criterion !~ /^(eq|le|ge)$/i) {
	die "Error: Invalid criterion value, $Criterion, specified in quartet, \"$Col,$DataType,$Criterion,$Value\", using \"--hightlight\" option: Valid values: le, ge, or eq\n";
      }
      if ($DataType =~ /^numeric$/i) {
	if (!IsFloat($Value)) {
	  die "Error: Invalid criterion value, $Value, specified in quartet, \"$Col,$DataType,$Criterion,$Value\", using \"--hightlight\" option: numeric value required for numeric data type\n";
	}
      }
      push @{$OptionsInfo{SpecifiedColIds}}, $Col;
      push @{$OptionsInfo{SpecifiedColDataTypes}}, $DataType;
      push @{$OptionsInfo{SpecifiedColCriteria}}, $Criterion;
      push @{$OptionsInfo{SpecifiedColValues}}, $Value;
    }
  }
}

# Retrieve information about input text files...
sub RetrieveTextFilesInfo {
  my($LineCount, $TextFile, $FileDir, $FileName, $HTMLFile, $CSSFile, $HTMLRoot, $HTMLTitle, $FileExt, $Index, $ColIndex, $ColNum, $ColLabel, $LinesCount, $InDelim, $Line, @LineWords, @ColLabels, $TopHTMLDir);

  %TextFilesInfo = ();

  @{$TextFilesInfo{FileOkay}} = ();
  @{$TextFilesInfo{ColCount}} = ();
  @{$TextFilesInfo{ColLabels}} = ();
  @{$TextFilesInfo{ColLabelToNumMap}} = ();
  @{$TextFilesInfo{LineCount}} = ();
  @{$TextFilesInfo{InDelim}} = ();

  @{$TextFilesInfo{HTMLRoot}} = ();
  @{$TextFilesInfo{HTMLTitle}} = ();
  @{$TextFilesInfo{MultipleHTMLTables}} = ();

  @{$TextFilesInfo{TopHTMLDir}} = ();
  @{$TextFilesInfo{SubHTMLDir}} = ();

  FILELIST: for $Index (0 .. $#TextFilesList) {
    $TextFile = $TextFilesList[$Index];

    $TextFilesInfo{FileOkay}[$Index] = 0;
    $TextFilesInfo{ColCount}[$Index] = 0;
    $TextFilesInfo{LineCount}[$Index] = 0;
    $TextFilesInfo{InDelim}[$Index] = "";
    $TextFilesInfo{HTMLRoot}[$Index] = "";
    $TextFilesInfo{HTMLTitle}[$Index] = "";
    $TextFilesInfo{MultipleHTMLTables}[$Index] = 0;

    @{$TextFilesInfo{ColLabels}[$Index]} = ();
    %{$TextFilesInfo{ColLabelToNumMap}[$Index]} = ();

    if (!(-e $TextFile)) {
      warn "Warning: Ignoring file $TextFile: It doesn't exist\n";
      next FILELIST;
    }
    if (!CheckFileType($TextFile, "csv tsv")) {
      warn "Warning: Ignoring file $TextFile: It's not a csv or tsv file\n";
      next FILELIST;
    }
    ($FileDir, $FileName, $FileExt) = ParseFileName($TextFile);
    if ($FileExt =~ /^tsv$/i) {
      $InDelim = "\t";
    }
    else {
      $InDelim = "\,";
      if ($Options{indelim} !~ /^(comma|semicolon)$/i) {
	warn "Warning: Ignoring file $TextFile: The value specified, $Options{indelim}, for option \"--indelim\" is not valid for csv files\n";
	next FILELIST;
      }
      if ($Options{indelim} =~ /^semicolon$/i) {
	$InDelim = "\;";
      }
    }

    if (!open TEXTFILE, "$TextFile") {
      warn "Warning: Ignoring file $TextFile: Couldn't open it: $! \n";
      next FILELIST;
    }

    $Line = GetTextLine(\*TEXTFILE);
    @ColLabels = quotewords($InDelim, 0, $Line);
    $LineCount = 0;
    while (<TEXTFILE>) {
      $LineCount++;
    }
    close TEXTFILE;

    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($TextFile);
    $HTMLRoot = $FileName;
    if ($Options{root} && (@TextFilesList == 1)) {
      my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($Options{root});
      if ($RootFileName && $RootFileExt) {
	$HTMLRoot = $RootFileName;
      }
      else {
	$HTMLRoot = $Options{root};
      }
    }
    $HTMLTitle = $HTMLRoot;
    if ($Options{title} && (@TextFilesList == 1)) {
      $HTMLTitle = $Options{title};
    }
    $HTMLFile = lc($HTMLRoot) . "-html";
    if (!$Options{overwrite}) {
      if (-d $HTMLFile) {
	warn "Warning: Ignoring file $TextFile: The directory $HTMLFile already exists\n";
	next FILELIST;
      }
    }

    $TextFilesInfo{FileOkay}[$Index] = 1;
    $TextFilesInfo{InDelim}[$Index] = $InDelim;
    $TextFilesInfo{HTMLRoot}[$Index] = "$HTMLRoot";
    $TextFilesInfo{HTMLTitle}[$Index] = "$HTMLTitle";

    $TextFilesInfo{ColCount}[$Index] = @ColLabels;
    push @{$TextFilesInfo{ColLabels}[$Index]}, @ColLabels;
    for $ColNum (0 .. $#ColLabels) {
      $ColLabel = $ColLabels[$ColNum];
      $TextFilesInfo{ColLabelToNumMap}[$Index]{$ColLabel} = $ColNum;
    }
    $TextFilesInfo{LineCount}[$Index] = $LineCount;

    if ($Options{numrows} == 0 || $LineCount <= $Options{numrows}) {
      $TextFilesInfo{MultipleHTMLTables}[$Index] = 0;
    }
    else {
      $TextFilesInfo{MultipleHTMLTables}[$Index] = 1;
    }
    # Setup HTML data directories paths...
    $TopHTMLDir = lc($TextFilesInfo{HTMLRoot}[$Index]) . "-html";
    $TextFilesInfo{TopHTMLDir}[$Index] = "$TopHTMLDir";
    $TextFilesInfo{SubHTMLDir}[$Index] = "$TopHTMLDir\/html";
  }
}

# Setup various data directories to hold HTML and other related files...
sub SetupDataDirs {
  my($Index) = @_;
  my($TopHTMLDir, $SubHTMLDir, $CreateTopHTMLDir, $CreateSubHTMLDir);

  $TopHTMLDir = $TextFilesInfo{TopHTMLDir}[$Index];
  $SubHTMLDir = $TextFilesInfo{SubHTMLDir}[$Index];

  # Clean up existing directories...
  if (-d $TopHTMLDir) {
    unlink "<$TopHTMLDir/*.html>";
    unlink "<$TopHTMLDir/*.css>";
  }
  if (-d $SubHTMLDir) {
    unlink "<$SubHTMLDir/*.html>";
  }
  # What directories need to be created...
  $CreateTopHTMLDir = (-d $TopHTMLDir) ? 0 : 1;

  $CreateSubHTMLDir = 0;
  if ($TextFilesInfo{MultipleHTMLTables}[$Index]) {
    $CreateSubHTMLDir = (-d $SubHTMLDir) ? 0 : 1;
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
  return ($TopHTMLDir, $SubHTMLDir);
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{indelim} = "comma";
  $Options{numrows} = 50;

  $Options{mode} = "shade";
  $Options{highlightby} = "colnum";
  $Options{highlightstyle} = "background";

  $Options{cellpadding} = 2;
  $Options{cellspacing} = 1;

  $Options{displaylinks} = "both";
  $Options{displaylinksinfo} = "both";
  $Options{stylesheet} = "new";

  $Options{titledisplay} = "yes";

  if (!GetOptions(\%Options, "align|a=s", "border|b=i", "cellpadding=i", "cellspacing=i", "color|c=s", "footer=s", "displaylinks|d=s", "displaylinksinfo=s", "help|h", "headeralign=s", "headercolor=s", "highlight=s", "highlightby=s", "highlightcolor=s", "highlightstyle=s", "indelim=s", "mode|m=s", "numrows|n=i", "overwrite|o", "root|r=s", "shadecolor=s", "stylesheet=s", "stylesheetname=s", "title|t=s", "titledisplay=s", "workingdir|w=s")) {
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
  if ($Options{displaylinksinfo} !~ /^(line|table|both)$/i) {
    die "Error: The value specified, $Options{displaylinksinfo}, for option \"--displaylinksinfo\" is not valid. Allowed values: line, table, or both\n";
  }
  if ($Options{indelim} !~ /^(comma|semicolon)$/i) {
    die "Error: The value specified, $Options{indelim}, for option \"--indelim\" is not valid. Allowed values: comma or semicolon\n";
  }
  if ($Options{highlightby} !~ /^(colnum|collabel)$/i) {
    die "Error: The value specified, $Options{highlightby}, for option \"--highlightby\" is not valid. Allowed values: colnum or collabel\n";
  }
  if ($Options{highlightstyle} !~ /^(background|text)$/i) {
    die "Error: The value specified, $Options{highlightstyle}, for option \"--highlightstyle\" is not valid. Allowed values: background or text\n";
  }
  if ($Options{mode} !~ /^(plain|shade|highlight|shadedhighlight)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: plain, shade, hightlight, or shadedhighlight\n";
  }
  if ($Options{stylesheet} !~ /^(old|new|none)$/i) {
    die "Error: The value specified, $Options{stylesheet}, for option \"-s --stylesheet\" is not valid. Allowed values: old, new, or none\n";
  }
  if ($Options{numrows} < 0) {
    die "Error: The value specified, $Options{numrows},  for option \"-n --numrows\" is not valid. Allowed values: >= 0 \n";
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

TextFilesToHTML.pl - Generate HTML table file(s) from TextFile(s)

=head1 SYNOPSIS

TextFilesToHTML.pl ... TextFile(s)...

TextFilesToHTML.pl [B<-a, --align> left | center | right,[top | middle | bottom]] [B<-b, --border> borderwidth] [B<--cellpadding> padding]
[B<--cellspacing> spacing] [B<--footer> string] [B<-d, --displaylinks> top | bottom | both]
[B<--displaylinksinfo> line | table | both] [B<-h, --help>]
[B<--headeralign> left | center | right,[top | middle | bottom]] [B<--headercolor> "#RRGGBB"]
[B<--highlight> "fieldlabel,datatype,criterion,value,[fieldlabel,datatype,criterion,value,]..."]
[B<--highlightby> colnum | collabel] [B<--highlightcolor> "#RRGGBB,#RRGGBB"]
[B<--highlightstyle> text | background] [B<--indelim> comma | semicolon] [B<-m, --mode> plain | shade | highlight | shadedhighlight]
[B<-n, --numrows> number] [B<-o, --overwrite>] [B<-r, --root> rootname]
[B<--stylesheet> old | new | none] [B<--stylesheetname> filename] [B< --shadecolor> "#RRGGBB,#RRGGBB"]
[B<-t, --title> string] [B<--titledisplay> yes | no] [B<-w, --workingdir> dirname] TextFile(s)...

=head1 DESCRIPTION

Generate HTML file(s) from I<TextFile(s)>. The HTML file(s) contain data tables and appropriate
navigational links to view other tables. These files can be generated for local viewing or
deployment on a web server. A variety of options are provided to control style and
appearence of tables.

Multiple I<TextFile(s)> names are separated by spaces. The valid file extensions are I<.csv> and
I<.tsv> for comma/semicolon and tab delimited text files respectively. All other file names
are ignored. All the text files in a current directory can be specified by I<*.csv>,
I<*.tsv>, or the current directory name. The B<--indelim> option determines the
format of I<TextFile(s)>. Any file which doesn't correspond to the format indicated
by B<--indelim> option is ignored.

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

=item B<--footer> I<string>

Text string to be included at bottom of each HTML file. Default: none.

=item B<-d, --displaylinks> I<top | bottom | both>

Specify where to display navigation links in each HTML file for accessing all other HTML
files. Possible values: I<top, bottom, or both>. Default value: I<both>. This option is
only valid during multiple HTML files generation for an input file.

=item B<--displaylinksinfo> I<line | table | both>

Control display of additional information along with navigational links: Showing line
n of m is displyed for line and showing table n of m for table. Possible values: I<line
| table | both>. Default: I<both>. This option is only valid  during multiple HTML files generation.

=item B<-h, --help>

Print this help message

=item B<--headeralign> I<left | center | right,[top | middle | bottom]>

Horizontal and vertical alignment for table header rows. Possible horizontal alignment
values: I<left, center, or right>. Possible vertical alignment values: I<top, middle, or bottom>.

Default values: I<center,middle>

=item B<--headercolor> I<"#RRGGBB">

Color used to fill background of table header row containing column labels
represented as a hexadecimal string. None for B<-m, --mode> option value
of I<plain> and I<#ccccff>, light blue, for others.

=item B<--highlight> I<"fieldlabel,datatype,criterion,value,[fieldlabel,datatype,criterion,value,]...">

This value is mode specific. It specifies how to highlight various column values
for each text file. Same set of quartets values are applied to all I<TextFile(s)>.

For I<highlightbycolnum> mode, input text format contains these quartets:
I<colnum,datatype,criterion,value,...>. Possible datatype values: I<numeric or text>.
Possible criterion values: I<le, ge, or eq>. Examples: "1,numeric,le,450>" or
"2,numeric,ge,150,6,numeric,le,10".

For I<highlightbycollabel> mode, input text format contains these quartets:
I<collabel,datatype,criterion,value,...>.

=item B<--highlightby> I<colnum | collabel>

This value is mode specific. It indicates how columns to be highlighted are specified
using B<--hightlight> option. Possible values: I<colnum or collabel>. Default value: I<colnum>.

=item B<--highlightcolor> I<"#RRGGBB,#RRGGBB">

Colors used to highlight column values during I<highlight> and I<shadedhightlight>
mode represented as hexadecimal strings.

For B<--highlighstyle> option values of I<text> and I<background>, these colors represent
text or background colors respectively. For a specific column, first color string is used for
values which meet criterion indicated by B<--highlight> option; the second color is used
for rest of the values.

Default values for I<background> B<--highlightstyle>: I<#0fff0f,#ff0f0f>. And default values for
I<text> B<--highlightstyle>: I<#0fbb0f,#ff0f0f>. Hexadecimal strings for both B<--highlightstyle>
colors correspond to I<reddish> and I<greenish>.

=item B<--highlightstyle> I<text | background>

This value is mode specific. It indicates highlight style used to differentiate column
values which pass a specified criterion from others. Possible values: I<text or
background>. Default: I<background>.

=item B<--indelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile(s)>. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a
delimiter.

=item B<-m, --mode> I<plain | shade | highlight | shadedhighlight>

Specify how to generate HTML table(s): plain tables with line borders, background of
alternate rows filled with a specified color, column values hightlighted using a specified
criteria, or combination of previous two styles.

Possible values: I<plain, shade, highlight, or shadedhighlight>. Default: I<shade>.

=item B<-n, --numrows> I<number>

Maximum number of rows per table. Default value: I<100>. Use 0 to put all rows into
one table. For I<TextFile(s)> with more than maximum number of specified rows,
multiple HTML tables, with appropriate navigation links, are created.

=item B<-o, --overwrite>

Overwrite existing files.

=item B<-r, --root> I<rootname>

New file or directory name is generated using the root: <root>.html or <root>-html.
Default new file name: <InitialTextFileName>.html. Default directory name:
<InitialTextFileName>-html.

For I<TextFile(s)> with more than maximum number of rows specified per table,
this directory tree is generated using <Name> where <Name> corresponds to <root>
or <InitialTextFileName>:Top dir - <Name>-html; Sub dirs - html and mols. <Top dir> contains
<Name>.html and <Name>.css files and <sub dir> html conatins various
<Name>Lines<Start>To<End>.html files; <sub dir> mols is created as needed and contains

This option is ignored for multiple input files.

=item B<--stylesheet> I<old | new | none>

Controls usage of stylesheet for newly generated HTML file(s). Possible values: I<old,
new, or none>. Default value: I<new>.

Stylesheet file contains various properties which control apperance of HTML pages:
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
file. Example: "/stylesheets/MyStyleSheet.css". Or a complete path name to a local
stylesheet file.

For I<create> value of B<-s --stylesheet> option, a new stylesheet file is created using
B<-r --root> option. And value of B<--stylesheetname> is simply ignored.

=item B< --shadecolor> I<"#RRGGBB,#RRGGBB">

Colors used to fill background of rows during I<shade> and I<shadedhightlight> mode
represented as a pair of hexadecimal string; the first and second color values
are used for odd and even number rows respectively.

Default value: I<"#ffffff,#e0e9eb"> - it's white and very light blue for odd and even number rows.

=item B<-t, --title> I<string>

Title for HTML table(s). Default value: <TextFileName>. For multiple input files,
B<-r --root> option is used to generate appropriate titles.

=item B<--titledisplay> I<yes | no>

Display title for HTML table(s). Possible values: I<yes or no>. Default value: I<yes>.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To generate HTML tables with rows background filled with white and greyish colors and
navigation links on top and botton of each page, type:

    % TextFilesToHTML.pl -o Sample1.csv

To generate HTML tables with rows background filled with golden and greyish colors,
navigation links on top and botton of each page, 10 rows in each table, greyish header
row color, and cell spacing of 1, type:

    % TextFilesToHTML.pl -o -n 10 --headeralign "center" --headercolor
      "#a1a1a1" --shadecolor "#ddd700,#d1d1d1" --cellspacing 1
      Sample1.csv

To generate plain HTML tables with 10 rows in each table and navigation links only at
the bottom, type:

    % TextFilesToHTML.pl -o -n 10 --displaylinks bottom -m plain
      Sample1.csv

To highlight values in column 3 using specified highlight criteria and fill in default background
colors, type:

    % TextFilesToHTML.pl -n 10 --highlight "3,numeric,le,450"
      --highlightby colnum --highlightstyle background -m
      shadedhighlight -o Sample1.csv

To highlight values in column MolWeight using specified highlight criteria, color the text using
default colors, and add a footer message in every page, type:

    % TextFilesToHTML.pl -n 4 --highlight "MolWeight,numeric,le,500"
      --highlightby collabel --highlightstyle text -m shadedhighlight -o
      --footer "Copyright (C) MayaChemTools" --cellspacing 1 Sample1.csv

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

JoinTextFiles.pl, MergeTextFilesWithSD.pl, ModifyTextFilesFormat.pl, SplitTextFiles.pl, SortTextFiles.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
