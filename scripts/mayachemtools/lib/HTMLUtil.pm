package HTMLUtil;
#
# File: HTMLUtil.pm
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
use Exporter;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw(InsertHTMLTags SetupHTMLAlignmentBegin SetupHTMLAlignmentEnd SetupHTMLButtonRef SetupHTMLDivBegin SetupHTMLDivEnd SetupHTMLEmptyLines SetupHTMLPageHeader SetupHTMLHRef SetupHTMLPageEnd SetupHTMLPageTitle SetupHTMLStyleSheetTags SetupHTMLTableHeader SetupHTMLTableEnd SetupHTMLTableColumnHeader SetupHTMLTableColumnEnd SetupHTMLTableRowHeader SetupHTMLTableRowEnd SetupHTMLTableRowHeaderValue SetupHTMLTableRowDataValue SetupJavaScriptCmds SetupStrViewerJSInitCmd SetupStrViewerJMEApplet SetupStrViewerJmolApplet SetupStrViewerChimePlugIn SetupStrViewerChem3DActiveX SetupStrViewerChemDrawActiveX SetupStrViewerChemDrawPlugIn SetupStrViewerMarvinViewApplet SetupStrViewerAccelrysActiveX);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Default window size for various supported structure viewers...
my($StrViewerWidth, $StrViewerHeight) = (250, 170);

# Insert specfied tags into existing tag string...
sub InsertHTMLTags {
  my($Tag, %TagsMap) = @_;
  my($NewTag, $TagName, $TagValue, $TagPart1, $TagPart2);

  $NewTag = $Tag; $TagPart1 = ""; $TagPart2 = "";
  ($TagPart1) = $Tag =~ /^(.*?)>/;

  if ($TagPart1 && length($TagPart1)) {
    $TagPart2 = $Tag;
    $TagPart2 =~ s/^(.*?)>//;
    if ($TagPart2 && length($TagPart2)) {
      for $TagName (keys %TagsMap) {
	$TagValue = $TagsMap{$TagName};
	$TagPart1 .= qq( $TagName="$TagValue" );
      }
      $NewTag = "${TagPart1}>${TagPart2}";
    }
  }

  return $NewTag;
}

sub SetupHTMLAlignmentBegin {
  my($AlignmentTag, $Alignment);

  $Alignment = (@_ == 1) ? $_[0] : "left";
  $AlignmentTag = qq(<$Alignment>\n);

  return $AlignmentTag;
}

sub SetupHTMLAlignmentEnd {
  my($AlignmentTag, $Alignment);

  $Alignment = (@_ == 1) ? $_[0] : "left";
  $AlignmentTag = qq(</$Alignment>\n);

  return $AlignmentTag;
}

# Setup a button reference...
sub SetupHTMLButtonRef {
  my($ButtonLabel, $RefFile, $ButtonTags);

  ($ButtonLabel, $RefFile) = @_;

  $ButtonTags = qq(<input type="button" value="$ButtonLabel" onClick="document.location='$RefFile'">);
  return $ButtonTags;
}

sub SetupHTMLDivBegin {
  my($Id) = @_;
  my($DivTag);

  $DivTag = qq(<div id="$Id">\n);

  return $DivTag;
}

sub SetupHTMLDivEnd {
  my($DivTag);

  $DivTag = qq(</div>\n);

  return $DivTag;
}
sub SetupHTMLEmptyLines {
  my($LineCount, $Index, $EmptyLineTags);

  $LineCount = 1;
  $EmptyLineTags = qq(<p>&nbsp</p>);
  ($LineCount) = @_;
  if ($LineCount > 1) {
    for $Index (2 .. $LineCount) {
      $EmptyLineTags .= qq(<p>&nbsp</p>);
    }
  }
  return $EmptyLineTags;
}

# Setup HTML page header...
sub SetupHTMLPageHeader {
  my($HeaderTitle, $Stylesheet, $JavaScript, $PageHeader);

  $HeaderTitle = "";  $Stylesheet = ""; $JavaScript = "";
  if (@_ == 3) {
    ($HeaderTitle, $Stylesheet, $JavaScript) = @_;
  }
  elsif (@_ == 2) {
    ($HeaderTitle, $Stylesheet) = @_;
  }
  else {
    ($HeaderTitle) = @_;
  }
  $PageHeader = qq(<html>\n);
  $PageHeader .= qq(<head>\n);
  $PageHeader .= qq(<title>$HeaderTitle</title>\n);
  $PageHeader .= qq(<meta http-equiv="content-type" content="text/html;charset=utf-8">\n);
  if ($Stylesheet) {
    $PageHeader .= qq(<link rel="stylesheet" type="text/css" href="$Stylesheet">\n);
  }
  if ($JavaScript) {
    $PageHeader .= qq(<script src="$JavaScript"></script>\n);
  }
  $PageHeader .= <<ENDPAGEHEADER;
</head>
<body>
<p>&nbsp</p>
ENDPAGEHEADER

  return $PageHeader;
}

# Setup page title...
sub SetupHTMLPageTitle {
  my($Title, $Alignment, $PageTitle);

  $Alignment = "center";
  if (@_ == 2) {
    ($Title, $Alignment) = @_;
  }
  else {
    ($Title) = @_;
  }

  $PageTitle=<<ENDPAGETITLE;
<$Alignment>
<h3>$Title</h3>
</$Alignment>
ENDPAGETITLE

  return $PageTitle;
}

# Setup HTML page end...
sub SetupHTMLPageEnd {
  my($PageEnd, $Footer);

  $Footer = "";
  if (@_ == 1) {
    ($Footer) = @_;
  }
  if ($Footer) {
    $Footer = qq(<span class="Footer">$Footer</span>);
  }
  $PageEnd=<<ENDPAGE;
<center>
<p>&nbsp</p>
$Footer
</center>
</body>
</html>
ENDPAGE

  return $PageEnd;
}

# Setup HTML link tags...
sub SetupHTMLHRef {
  my($Value, $RefFile, $HRef, $Title);

  $Title = "";
  if (@_ == 3) {
    ($Value, $RefFile, $Title) = @_;
  }
  else {
    ($Value, $RefFile) = @_;
  }

  $HRef = qq(<a href="$RefFile");
  if ($Title) {
    $HRef .= qq( title="$Title");
  }
  $HRef .= qq(>$Value</a>);
  return $HRef;
}

#
sub SetupHTMLStyleSheetTags {
  my($StyleSheetTags);

  $StyleSheetTags=<<ENDSTYLESHEET;
body
{
    background-color: #ffffff;
    font-family: Verdana, Arial, Helvetica, sans-serif;
    font-size: 11px;
}
p
{
    font-family: Verdana, Arial, Helvetica, sans-serif;
    font-size: 11px;
}
h1
{
    font-family: Verdana, Arial, Helvetica, sans-serif;
    font-size: 25px;
    font-weight: bold;
    color: #0054aa;
}
h2
{
    font-family: Verdana, Arial, Helvetica, sans-serif;
    font-size: 18px;
    font-weight: bold;
    color: #0054aa;
}
h3
{
    font-family: Verdana, Arial, Helvetica, sans-serif;
    font-size: 14px;
    font-weight: bold;
    color: #0054aa;
}
b
{
    font-weight: bold;
}
td
{
    font-family: Verdana, Arial, Helvetica, sans-serif;
    font-size: 11px;
}
th
{
    font-family: Verdana, Arial, Helvetica, sans-serif;
    font-size: 11px;
    color: #0054aa;
    font-weight: bold;
}
.box {
  border-color: #000000;
  border-style: solid;
  border-top-width: 1px;
  border-bottom-width: 1px;
  border-left-width: 1px;
  border-right-width: 1px;
}
a
{
    color: #0000bb;
    text-decoration: none;
    font-family: Verdana, Arial, Helvetica, sans-serif;
    font-size: 11px;
}
a:hover
{
    color: #ff0000;
}
#tablenav {
    font-family: Verdana, Arial, Helvetica, sans-serif;
    font-size: 11px;
}
#tablenav td
{
    font-family: Verdana, Arial, Helvetica, sans-serif;
    font-size: 11px;
}
#tablenav th
{
    font-family: Verdana, Arial, Helvetica, sans-serif;
    font-size: 11px;
    font-weight: bold;
}
#tablenav a
{
    color: #0000bb;
    text-decoration: none;
    font-family: Verdana, Arial, Helvetica, sans-serif;
    font-size: 11px;
}
#tablenav a:hover
{
    color: #ff0000;
}
.footer
{
    font-family: Arial, Verdana, Helvetica, sans-serif;
    font-size: 9px;
    color: #888888;
}
ENDSTYLESHEET

  return $StyleSheetTags;
}

# Setup HTML table header...
sub SetupHTMLTableHeader {
  my($TableHeader, $BorderWidth, $CellPadding, $CellSpacing, $Width, $Height);

  $BorderWidth = 1; $CellPadding = 2; $CellSpacing = 0; $Width = ""; $Height = "";
  if (@_ == 5) {
    ($BorderWidth, $CellPadding, $CellSpacing, $Width, $Height) = @_;
  }
  elsif (@_ == 4) {
    ($BorderWidth, $CellPadding, $CellSpacing, $Width) = @_;
  }
  elsif (@_ == 3) {
    ($BorderWidth, $CellPadding, $CellSpacing) = @_;
  }
  elsif (@_ == 2) {
    ($BorderWidth, $CellPadding) = @_;
  }
  elsif (@_ = 1) {
    ($BorderWidth) = @_;
  }
  $TableHeader = qq(<table border=$BorderWidth cellpadding=$CellPadding cellspacing=$CellSpacing);
  if ($Width) {
    $TableHeader .= qq( width=$Width);
  }
  if ($Height) {
    $TableHeader .= qq( height=$Height);
  }
  $TableHeader .= qq(>\n);

  return $TableHeader;
}

# Setup HTML table end...
sub SetupHTMLTableEnd {
  my($TableEnd);

  $TableEnd=<<ENDTABLE;
</table>
ENDTABLE

  return $TableEnd;
}

# Setup HTML table column header...
sub SetupHTMLTableColumnHeader {
  my($BgColor, $Width, $ColumnHeader);

  $BgColor = ""; $Width = "";
  if (@_ == 1) {
    ($BgColor) = @_;
  }
  elsif (@_ == 2) {
    ($BgColor, $Width) = @_;
  }
  $ColumnHeader = qq(<td);
  if ($BgColor) {
    $ColumnHeader .= qq( bgcolor="$BgColor")
  }
  if ($Width) {
    $ColumnHeader .= qq( width="$Width")
  }
  $ColumnHeader .= qq(>);
  return $ColumnHeader;
}

# Setup HTML table column end...
sub SetupHTMLTableColumnEnd {
  my($ColumnEnd);

  $ColumnEnd = qq(</td>);
  return $ColumnEnd;
}

# Setup HTML table row header...
sub SetupHTMLTableRowHeader {
  my($RowHeader, $HAlignment, $BgColor, $VAlignment);

  $HAlignment = "center"; $BgColor = ""; $VAlignment = "top";
  if (@_ == 3) {
    ($HAlignment, $BgColor, $VAlignment) = @_;
  }
  elsif (@_ == 2) {
    ($HAlignment, $BgColor) = @_;
  }
  elsif (@_ == 1) {
    ($HAlignment) = @_;
  }
  if ($BgColor) {
    $RowHeader = qq(<tr bgcolor="$BgColor" align="$HAlignment" valign="$VAlignment">);
  }
  else {
    $RowHeader = qq(<tr align="$HAlignment" valign="$VAlignment">);
  }

  return $RowHeader;
}

# Setup HTML table row end...
sub SetupHTMLTableRowEnd {
  my($RowEnd);

  $RowEnd = qq(</tr>\n);
  return $RowEnd;
}

# Setup HTML table header values...
sub SetupHTMLTableRowHeaderValue {
  my($HeaderValue, $Value);

  $Value = "";
  if (@_ >= 1) {
    ($Value) = @_;
  }
  if (defined $Value && length $Value) {
    $HeaderValue = qq(<th>$Value</th>);
  }
  else {
    $HeaderValue = qq(<th>&nbsp</th>);
  }
  return $HeaderValue;
}

# Setup HTML table row data values...
sub SetupHTMLTableRowDataValue {
  my($RowValues, $Value, $BgColor, $FontColor, $FontBold, $FontBoldTag1, $FontBoldTag2);

  $Value = ""; $BgColor = ""; $FontColor = ""; $FontBold = "";
  if (@_ == 1) {
    ($Value) = @_;
  }
  elsif (@_ == 2) {
    ($Value, $BgColor) = @_;
  }
  elsif (@_ == 3) {
    ($Value, $BgColor, $FontColor) = @_;
  }
  elsif (@_ == 4) {
    ($Value, $BgColor, $FontColor, $FontBold) = @_;
  }
  if (!(defined $Value && length $Value)) {
    $Value = qq(&nbsp);
  }
  $FontBoldTag1 = ""; $FontBoldTag2 = "";
  if ($FontBold) {
    $FontBoldTag1 = qq(<b>);
    $FontBoldTag2 = qq(</b>);
  }
  if ($BgColor || $FontColor) {
    my ($BgColorTag, $FontTag1, $FontTag2);

    $BgColorTag = "";
    if ($BgColor) {
      $BgColorTag = qq( bgcolor="$BgColor");
    }
    $FontTag1 = ""; $FontTag2 = "";
    if ($FontColor) {
      $FontTag1 = qq(<font color="$FontColor">);
      $FontTag2 = qq(</font>);
    }
    if ($FontBold) {
      $RowValues = "<td" . $BgColorTag . ">" . $FontBoldTag1 . $FontTag1 . "$Value" . $FontTag2 . $FontBoldTag2 .  "</td>";
    }
    else {
      $RowValues = "<td" . $BgColorTag . ">" . $FontTag1 . "$Value" . $FontTag2 .  "</td>";
    }
  }
  elsif ($FontBold) {
    $RowValues = "<td>" . $FontBoldTag1 . "$Value" . $FontBoldTag2 . "</td>";
  }
  else {
    $RowValues = qq(<td>$Value</td>);
  }
  return $RowValues;
}

# Setup Java scripts command...
sub SetupJavaScriptCmds {
  my(@JSCmdList) = @_;
  my($JSTags, $JSCmd);

  $JSTags = qq(<script>\n);
  for $JSCmd (@JSCmdList) {
    $JSTags .= qq($JSCmd\n);
  }
  $JSTags .= qq(</script>\n);

  return $JSTags;
}

# Setup Java script initialize command...
sub SetupStrViewerJSInitCmd {
  my($StrViewerType, $CodeBase) = @_;
  my($JSTag);

  $JSTag = "";
  if ($StrViewerType eq "Jmol") {
    $JSTag = qq(<script>jmolInitialize("$CodeBase", "JmolApplet.jar");</script>\n);
  }
  elsif ($StrViewerType eq "ChemDrawPlugIn" || $StrViewerType eq "ChemDrawActiveX") {
    $JSTag = qq(<script>cd_includeWrapperFile("$CodeBase/");</script>\n);
  }
  elsif ($StrViewerType eq "Chem3DActiveX") {
  }
  return $JSTag;
}


# Setup Jmol applet...
sub SetupStrViewerJmolApplet {
  my($MolString, $CodeBase, $ParamsMapRef, %ParamsMap, $AppletTags, $JavaScriptTags, $ReturnTags, $Name, $Code, $Archive, $Width, $Height, $ParamName, $ParamValue, $JSFileName, $UseJavaScript);
  my($ProgressBar, $ProgressColor, $BoxMessage, $BoxFgColor, $BoxBgColor, $BgColor, $JMolScript);

  $AppletTags = "";  $JavaScriptTags = ""; $ReturnTags = "";
  $ParamsMapRef = ""; %ParamsMap = ();
  $Name = "Jmol"; $Code = "JmolApplet"; $Archive = "JmolApplet.jar"; $Width = $StrViewerWidth; $Height = $StrViewerHeight;
  $ProgressBar = "true"; $ProgressColor = "#0000ff"; $BgColor = "#000000";
  $BoxMessage = "Setting up JmolApplet..."; $BoxFgColor = "#000000"; $BoxBgColor = "#ffffff";
  $UseJavaScript = 0; $JSFileName = "";
  $JMolScript = "select *; set frank off; wireframe on; spacefill off";

 PARAMS: {
    if (@_ == 3) { ($MolString, $CodeBase, $ParamsMapRef) = @_; last PARAMS; }
    ($MolString, $CodeBase) = @_;
  }

  if ($ParamsMapRef) {
    %ParamsMap = %$ParamsMapRef;
    if (exists $ParamsMap{usejavascript} ) { $JSFileName = $ParamsMap{usejavascript}; $UseJavaScript = 1; $ParamsMap{usejavascript} = ""; }
    if (exists $ParamsMap{name} ) { $Name = $ParamsMap{name}; $ParamsMap{name} = ""; }
    if (exists $ParamsMap{code} ) { $Code = $ParamsMap{code}; $ParamsMap{code} = ""; }
    # if (exists $ParamsMap{archive} ) { $Archive = $ParamsMap{archive}; $ParamsMap{archive} = ""; }
    if (exists $ParamsMap{archive} ) { $ParamsMap{archive} = ""; }
    if (exists $ParamsMap{height} ) { $Height = $ParamsMap{height}; $ParamsMap{height} = ""; }
    if (exists $ParamsMap{width} ) { $Width = $ParamsMap{width}; $ParamsMap{width} = ""; }
    if (exists $ParamsMap{progressbar} ) { $ProgressBar = $ParamsMap{progressbar}; $ParamsMap{progressbar} = ""; }
    if (exists $ParamsMap{progresscolor} ) { $ProgressColor = $ParamsMap{progresscolor}; $ParamsMap{progresscolor} = ""; }
    if (exists $ParamsMap{boxmessage} ) { $BoxMessage = $ParamsMap{boxmessage}; $ParamsMap{boxmessage} = ""; }
    if (exists $ParamsMap{script} ) { $JMolScript = $ParamsMap{script}; $ParamsMap{script} = ""; }
    if (exists $ParamsMap{bgcolor}) {
      $BgColor = $ParamsMap{bgcolor};
      if (length($BgColor)) {
	if ($BgColor =~ /black/i || $BgColor =~ /#000000/ ) {
	  $BoxFgColor = "#ffffff";
	  $BoxBgColor = "#000000";
	}
      }
    }
    if (exists $ParamsMap{boxbgcolor} ) { $BoxBgColor = $ParamsMap{boxbgcolor}; $ParamsMap{boxbgcolor} = ""; }
    if (exists $ParamsMap{boxfgcolor} ) { $BoxFgColor = $ParamsMap{boxfgcolor}; $ParamsMap{boxfgcolor} = ""; }
  }

  $MolString =~ s/(\r\n)|(\r)|(\n)/|/g;
  if ($UseJavaScript) {
    $JavaScriptTags = qq(\n<script>\n);
    my($Size) = ($Width > $Height ) ? $Width : $Height;
    $JavaScriptTags .= qq(var $Name = \n);
    my(@MolLines) = split /\|/, $MolString;
    my($LineIndex);
    $JavaScriptTags .= qq("$MolLines[0]\\n");
    for $LineIndex (1 .. $#MolLines) {
      $JavaScriptTags .= qq( + \n"$MolLines[$LineIndex]\\n");
    }
    $JavaScriptTags .= qq(;\n);
    $JavaScriptTags .= qq(jmolSetAppletColor("$BgColor", "$BoxBgColor", "$BoxFgColor", "$ProgressColor");\n);
    # "set frank off turns" off JMol logo. For wireframe display; use wireframe on; spacefill off...
    # $JavaScriptTags .= qq(jmolAppletInline($Size, $Name);\n);
    $JavaScriptTags .= qq(jmolAppletInline([$Width,$Height], $Name, \"$JMolScript\");\n);
    $JavaScriptTags .= qq(</script>\n);
    $ReturnTags = $JavaScriptTags;
  }
  else {
    # Setup applet header...
    $AppletTags = qq(\n<applet name="$Name" id="$Name" code="$Code" archive="$Archive" codebase="$CodeBase" width="$Width" height="$Height">\n);

    # Setup molecule data...
    $AppletTags .= qq(<param name="loadInline" value="$MolString">\n);

    # Setup prograss bar...
    $AppletTags .= qq(<param name="progressbar" value="$ProgressBar">\n<param name="progresscolor" value="$ProgressColor">\n<param name="boxmessage" value="$BoxMessage">\n<param name="boxbgcolor" value="$BoxBgColor">\n<param name="boxfgcolor" value="$BoxFgColor">\n);

    # "set frank off turns" off JMol logo. For wireframe display; use wireframe on; spacefill off...
    $AppletTags .= qq(<param name="script" value="$JMolScript">);

    #Setup other parameters...
    for $ParamName (sort keys %ParamsMap) {
      $ParamValue = $ParamsMap{$ParamName};
      if (length $ParamValue) {
	$AppletTags .= qq(<param name="$ParamName" value="$ParamValue">\n);
      }
    }
    #Finish it up...
    $AppletTags .= qq(</applet>\n);
    $ReturnTags = $AppletTags;
  }
  return $ReturnTags;
}

# Setup JME applet...
sub SetupStrViewerJMEApplet {
  my($MolString, $CodeBase, $ParamsMapRef, %ParamsMap, $AppletTags, $Name, $Code, $Archive, $Width, $Height, $ParamName, $ParamValue);
  my($Options);

  $AppletTags = "";  $ParamsMapRef = ""; %ParamsMap = ();
  $Name = "JME"; $Code = "JME"; $Archive = "JME.jar"; $Width = $StrViewerWidth; $Height = $StrViewerHeight;
  $Options = "depict";

  if (@_ == 3) {
    ($MolString, $CodeBase, $ParamsMapRef) = @_;
  }
  else {
    ($MolString, $CodeBase) = @_;
  }
  $MolString =~ s/(\r\n)|(\r)|(\n)/|/g;

  if ($ParamsMapRef) {
    %ParamsMap = %$ParamsMapRef;
    if (exists $ParamsMap{name} ) { $Name = $ParamsMap{name}; $ParamsMap{name} = ""; }
    if (exists $ParamsMap{code} ) { $Code = $ParamsMap{code}; $ParamsMap{code} = ""; }
    if (exists $ParamsMap{archive} ) { $Archive = $ParamsMap{archive}; $ParamsMap{archive} = ""; }
    if (exists $ParamsMap{height} ) { $Height = $ParamsMap{height}; $ParamsMap{height} = ""; }
    if (exists $ParamsMap{width} ) { $Width = $ParamsMap{width}; $ParamsMap{width} = ""; }
    if (exists $ParamsMap{options} ) { $Options = $ParamsMap{options}; $ParamsMap{options} = ""; }
  }

  # Setup applet header...
  $AppletTags = qq(\n<applet name="$Name" id="$Name" code="$Code" archive="$Archive" codebase="$CodeBase" width="$Width" height="$Height">\n);

  # Setup molecule data...
  $AppletTags .= qq(<param name="mol" value="$MolString">\n<param name="options" value="$Options">\n);

  #Setup other parameters...
  for $ParamName (sort keys %ParamsMap) {
    $ParamValue = $ParamsMap{$ParamName};
    if (length $ParamValue) {
      $AppletTags .= qq(<param name="$ParamName" value="$ParamValue">\n);
    }
  }

  #Finish it up...
  $AppletTags .= qq(</applet>\n);

  return $AppletTags;
}

# Setup MarvinView applet...
sub SetupStrViewerMarvinViewApplet {
  my($MolString, $CodeBase, $ParamsMapRef, %ParamsMap, $AppletTags, $JavaScriptTags, $ReturnTags, $Name, $Code, $Archive, $Width, $Height, $ParamName, $NavMode, $ParamValue, $JSFileName, $UseJavaScript, $MarvinJVM);

  $AppletTags = "";  $JavaScriptTags = ""; $ReturnTags = "";
  $ParamsMapRef = ""; %ParamsMap = ();
  $Name = "MView"; $Code = "MView"; $Archive = "marvin.jar"; $Width = $StrViewerWidth; $Height = $StrViewerHeight;
  $NavMode = "zoom";
  $UseJavaScript = 0; $JSFileName = ""; $MarvinJVM = "";

  if (@_ == 3) {
    ($MolString, $CodeBase, $ParamsMapRef) = @_;
  }
  else {
    ($MolString, $CodeBase) = @_;
  }

  if ($ParamsMapRef) {
    %ParamsMap = %$ParamsMapRef;
    if (exists $ParamsMap{usejavascript} ) { $JSFileName = $ParamsMap{usejavascript}; $UseJavaScript = 1; $ParamsMap{usejavascript} = ""; }
    if (exists $ParamsMap{marvin_jvm} ) { $MarvinJVM = $ParamsMap{marvin_jvm}; $ParamsMap{marvin_jvm} = ""; }
    if (exists $ParamsMap{name} ) { $Name = $ParamsMap{name}; $ParamsMap{name} = ""; }
    if (exists $ParamsMap{code} ) { $Code = $ParamsMap{code}; $ParamsMap{code} = ""; }
    if (exists $ParamsMap{archive} ) { $Archive = $ParamsMap{archive}; $ParamsMap{archive} = ""; }
    if (exists $ParamsMap{height} ) { $Height = $ParamsMap{height}; $ParamsMap{height} = ""; }
    if (exists $ParamsMap{width} ) { $Width = $ParamsMap{width}; $ParamsMap{width} = ""; }
    if (exists $ParamsMap{bgcolor}) {
      $ParamsMap{background} = "$ParamsMap{bgcolor}";
      $ParamsMap{molbg} = "$ParamsMap{bgcolor}";
      $ParamsMap{bgcolor} = "";
    }
    if (exists $ParamsMap{navmode}) {
      $NavMode = $ParamsMap{navmode};
      $ParamsMap{navmode} = "";
    }
  }
  $MolString =~ s/(\r\n)|(\r)|(\n)/\\/g;
  if ($UseJavaScript) {
    $JavaScriptTags = qq(\n<script>\n);
    $JavaScriptTags .= qq(var marvin_name = "$Name";\n);
    $JavaScriptTags .= qq(var marvin_jvm = "$MarvinJVM";\n);

    $JavaScriptTags .= qq(mview_begin("$CodeBase", $Width, $Height);\n);

    $JavaScriptTags .= qq(var $Name = \n);
    my(@MolLines) = split /\\/, $MolString;
    my($LineIndex);
    $JavaScriptTags .= qq("$MolLines[0]\\n");
    for $LineIndex (1 .. $#MolLines) {
      $JavaScriptTags .= qq( + \n"$MolLines[$LineIndex]\\n");
    }
    $JavaScriptTags .= qq(;\n);
    $JavaScriptTags .= qq(mview_param("mol", $Name);\n);

    $JavaScriptTags .= qq(mview_param("navmode", "$NavMode");\n);

    for $ParamName (sort keys %ParamsMap) {
      $ParamValue = $ParamsMap{$ParamName};
      if (length $ParamValue) {
	$JavaScriptTags .= qq(mview_param("$ParamName", "$ParamValue");\n);
      }
    }
    $JavaScriptTags .= qq(mview_end();\n);
    $JavaScriptTags .= qq(</script>\n);
    $ReturnTags = $JavaScriptTags;
  }
  else {
    # Setup applet header...
    $AppletTags = qq(\n<applet name="$Name" id="$Name" code="$Code" archive="$Archive" codebase="$CodeBase" width="$Width" height="$Height">\n);

    # Setup molecule data...
    $AppletTags .= qq(<param name="mol" value="$MolString">\n);

    $AppletTags .= qq(<param name="navmode" value="$NavMode">\n);

    #Setup other parameters...
    for $ParamName (sort keys %ParamsMap) {
      $ParamValue = $ParamsMap{$ParamName};
      if (length $ParamValue) {
	$AppletTags .= qq(<param name="$ParamName" value="$ParamValue">\n);
      }
    }
    $AppletTags .= qq(</applet>\n);
    $ReturnTags = $AppletTags;
  }
  return $ReturnTags;
}

# Setup MDL chime plug-in...
sub SetupStrViewerChimePlugIn {
  my($MolFile, $ParamsMapRef, %ParamsMap, $Width, $Height, $ParamName, $ParamValue, $PlugInTags);
  my($Display2D);

  $PlugInTags = ""; $ParamsMapRef = ""; %ParamsMap = ();
  $Width = $StrViewerWidth; $Height = $StrViewerHeight;
  $Display2D = "true";

  if (@_ == 2) {
    ($MolFile, $ParamsMapRef) = @_;
  }
  else {
    ($MolFile) = @_;
  }

  if ($ParamsMapRef) {
    %ParamsMap = %$ParamsMapRef;
    if (exists $ParamsMap{height} ) { $Height = $ParamsMap{height}; $ParamsMap{height} = ""; }
    if (exists $ParamsMap{width} ) { $Width = $ParamsMap{width}; $ParamsMap{width} = ""; }
    if (exists $ParamsMap{display2d} ) { $Display2D = $ParamsMap{display2d}; $ParamsMap{display2d} = ""; }
  }
  # Start plug-in tag...
  $PlugInTags = qq(<embed src="$MolFile" width="$Width" height="$Height" display2d="$Display2D");

  #Setup other parameters...
  for $ParamName (sort keys %ParamsMap) {
    $ParamValue = $ParamsMap{$ParamName};
    if (length $ParamValue) {
      $PlugInTags .= qq( $ParamName="$ParamValue");
    }
  }

  # Finish it off...
  $PlugInTags .= qq( >);

  return $PlugInTags;
}

# Setup Accelrys ViewerActiveX controls...
sub SetupStrViewerAccelrysActiveX {
  my($MolFile, $ParamsMapRef, %ParamsMap, $ActiveXTags, $Name, $Width, $Height, $ParamName, $ParamValue);
  my($ClassId, $Convert2DTo3D, $Style, $Mouse);

  $ActiveXTags = "";  $ParamsMapRef = ""; %ParamsMap = ();
  $Name = "ViewerActiveX"; $Width = $StrViewerWidth; $Height = $StrViewerHeight;
  $ClassId = "clsid:90690CB6-BC07-11D4-AEF7-0050DA948176";
  $Convert2DTo3D = "0";
  $Mouse = 4;

  if (@_ == 2) {
    ($MolFile, $ParamsMapRef) = @_;
  }
  else {
    ($MolFile) = @_;
  }

  if ($ParamsMapRef) {
    %ParamsMap = %$ParamsMapRef;
    if (exists $ParamsMap{classid} ) { $ClassId = $ParamsMap{classid}; $ParamsMap{classid} = ""; }
    if (exists $ParamsMap{name} ) { $Name = $ParamsMap{name}; $ParamsMap{name} = ""; }
    if (exists $ParamsMap{height} ) { $Height = $ParamsMap{height}; $ParamsMap{height} = ""; }
    if (exists $ParamsMap{width} ) { $Width = $ParamsMap{width}; $ParamsMap{width} = ""; }
    if (exists $ParamsMap{Convert2Dto3D} ) { $Convert2DTo3D = $ParamsMap{Convert2Dto3D}; $ParamsMap{Convert2Dto3D} = ""; }
    if (exists $ParamsMap{Mouse} ) { $Mouse = $ParamsMap{Mouse}; $ParamsMap{Mouse} = ""; }
    if (exists $ParamsMap{bgcolor} ) {
      my($BgColor) = $ParamsMap{bgcolor};
      $ParamsMap{bgcolor} = "";
      # Get OLE color value: &aabbggrr&
      # Set it to white for now...
      $BgColor = "16777215";
      $ParamsMap{BackColor} = "$BgColor";
    }
  }
  $Style = qq(style="height: ) . $Height . qq(px; width: ) . $Width . qq(px");

  # Setup object header...
  $ActiveXTags = qq(\n<object id="$Name" classid="$ClassId" $Style>\n);

  # Setup molecule data...
  $ActiveXTags .= qq(<param name="Source" value="$MolFile">\n<param name="Mouse" value="$Mouse">\n<param name="Convert2Dto3D" value="$Convert2DTo3D">\n);

  #Setup other parameters...
  for $ParamName (sort keys %ParamsMap) {
    $ParamValue = $ParamsMap{$ParamName};
    if (length $ParamValue) {
      $ActiveXTags .= qq(<param name="$ParamName" value="$ParamValue">\n);
    }
  }

  # Finish it off...
  $ActiveXTags .= qq(</object>\n);

  return $ActiveXTags;
}

# Setup Chem3D ActiveX 8.0 control...
sub SetupStrViewerChem3DActiveX {
  my($MolFile, $ParamsMapRef, %ParamsMap, $ActiveXTags, $JavaScriptTags, $ReturnTags, $Name, $Width, $Height, $ParamName, $ParamValue);
  my($ClassId, $Style, $DisplayType, $RotationBars, $MovieController, $JSFileName, $UseJavaScript);

  $ActiveXTags = ""; $JavaScriptTags = ""; $ReturnTags = "";
  $ParamsMapRef = ""; %ParamsMap = ();
  $Name = "Chem3D"; $Width = $StrViewerWidth; $Height = $StrViewerHeight;
  $ClassId = "clsid:B7A6B8E4-3E8B-4D18-8F8F-B4057EFC784B";
  $DisplayType = "Ball&Stick";
  $RotationBars = "false";
  $MovieController = "false";

  if (@_ == 2) {
    ($MolFile, $ParamsMapRef) = @_;
  }
  else {
    ($MolFile) = @_;
  }

  if ($ParamsMapRef) {
    %ParamsMap = %$ParamsMapRef;
    if (exists $ParamsMap{usejavascript} ) { $JSFileName = $ParamsMap{usejavascript}; $UseJavaScript = 1; $ParamsMap{usejavascript} = ""; }
    if (exists $ParamsMap{classid} ) { $ClassId = $ParamsMap{classid}; $ParamsMap{classid} = ""; }
    if (exists $ParamsMap{name} ) { $Name = $ParamsMap{name}; $ParamsMap{name} = ""; }
    if (exists $ParamsMap{height} ) { $Height = $ParamsMap{height}; $ParamsMap{height} = ""; }
    if (exists $ParamsMap{width} ) { $Width = $ParamsMap{width}; $ParamsMap{width} = ""; }
    if (exists $ParamsMap{displaytype} ) { $DisplayType = $ParamsMap{displaytype}; $ParamsMap{displaytype} = ""; }
    if (exists $ParamsMap{rotationbars} ) { $RotationBars = $ParamsMap{rotationbars}; $ParamsMap{rotationbars} = ""; }
    if (exists $ParamsMap{moviecontroller} ) { $MovieController = $ParamsMap{moviecontroller}; $ParamsMap{moviecontroller} = ""; }
  }
  $Style = qq(style="height: ) . $Height . qq(px; width: ) . $Width . qq(px");

  if ($UseJavaScript) {
    #Setup parameters...
    my($Params) = "";
    for $ParamName (sort keys %ParamsMap) {
      $ParamValue = $ParamsMap{$ParamName};
      if (length $ParamValue) {
	$Params .= qq( $ParamName='$ParamValue');
      }
    }
    $JavaScriptTags = qq(\n<script>\n);
    $JavaScriptTags .= qq(c3d_insert3dStr("name='$Name' src='$MolFile' width='$Width' height='$Height' displaytype='$DisplayType' rotation_bars_visible='$RotationBars' movie_controller_visible='$MovieController' $Params");\n);
    $JavaScriptTags .= qq(</script>\n);
    $ReturnTags = $JavaScriptTags;
  }
  else {
    # Setup object header...
    $ActiveXTags = qq(\n<object id="$Name" classid="$ClassId" $Style>\n);

    # Setup molecule data...
    $ActiveXTags .= qq(<param name="src" value="$MolFile">\n<param name="displaytype" value="$DisplayType">\n<param name="rotationbars" value="$RotationBars">\n<param name="moviecontroller" value="$MovieController">\n);

    #Setup other parameters...
    for $ParamName (sort keys %ParamsMap) {
      $ParamValue = $ParamsMap{$ParamName};
      if (length $ParamValue) {
	$ActiveXTags .= qq(<param name="$ParamName" value="$ParamValue">\n);
      }
    }
    $ActiveXTags .= qq(</object>\n);
    $ReturnTags = $ActiveXTags;
  }
  return $ReturnTags;
}

# Setup ChemDraw ActiveX 8.0 control...
# Problems: "bgcolor" parameter doesn't work.
sub SetupStrViewerChemDrawActiveX {
  my($MolFile, $ParamsMapRef, %ParamsMap, $ActiveXTags, $JavaScriptTags, $ReturnTags, $Name, $Width, $Height, $ParamName, $ParamValue);
  my($ClassId, $Style, $ViewOnly, $ShrinkToFit, $ShowToolsWhenVisible, $JSFileName, $UseJavaScript);

  $ActiveXTags = "";  $JavaScriptTags = ""; $ReturnTags = "";
  $ParamsMapRef = ""; %ParamsMap = ();
  $Name = "ChemDraw"; $Width = $StrViewerWidth; $Height = $StrViewerHeight;
  $ClassId = "clsid:51A649C4-3E3D-4557-9BD8-B14C0AD44B0C";
  $ViewOnly = "1"; $JavaScriptTags = "";
  $ShrinkToFit = "1";
  $ShowToolsWhenVisible = "1";

  if (@_ == 2) {
    ($MolFile, $ParamsMapRef) = @_;
  }
  else {
    ($MolFile) = @_;
  }

  if ($ParamsMapRef) {
    %ParamsMap = %$ParamsMapRef;
    if (exists $ParamsMap{usejavascript} ) { $JSFileName = $ParamsMap{usejavascript}; $UseJavaScript = 1; $ParamsMap{usejavascript} = ""; }
    if (exists $ParamsMap{classid} ) { $ClassId = $ParamsMap{classid}; $ParamsMap{classid} = ""; }
    if (exists $ParamsMap{name} ) { $Name = $ParamsMap{name}; $ParamsMap{name} = ""; }
    if (exists $ParamsMap{height} ) { $Height = $ParamsMap{height}; $ParamsMap{height} = ""; }
    if (exists $ParamsMap{width} ) { $Width = $ParamsMap{width}; $ParamsMap{width} = ""; }
    if (exists $ParamsMap{ViewOnly} ) { $ViewOnly = $ParamsMap{ViewOnly}; $ParamsMap{ViewOnly} = ""; }
    if (exists $ParamsMap{ShrinkToFit} ) { $ShrinkToFit = $ParamsMap{ShrinkToFit}; $ParamsMap{ShrinkToFit} = ""; }
    if (exists $ParamsMap{ShowToolsWhenVisible} ) { $ShowToolsWhenVisible = $ParamsMap{ShowToolsWhenVisible}; $ParamsMap{ShowToolsWhenVisible} = ""; }
  }
  if ($UseJavaScript) {
    #Setup parameter...
    my($Params) = "";
    for $ParamName (sort keys %ParamsMap) {
      $ParamValue = $ParamsMap{$ParamName};
      if (length $ParamValue) {
	$Params .= qq( $ParamName='$ParamValue');
      }
    }
    $JavaScriptTags = qq(\n<script>\n);
    $JavaScriptTags .= qq(cd_insertObjectStr("name='$Name' src='$MolFile' width='$Width' height='$Height' shrinktofit='$ShrinkToFit' viewonly='$ViewOnly' $Params");\n);
    $JavaScriptTags .= qq(</script>\n);
    $ReturnTags = $JavaScriptTags;
  }
  else {
    $Style = qq(style="height: ) . $Height . qq(px; width: ) . $Width . qq(px");

    # Setup object header...
    $ActiveXTags = qq(\n<object id="$Name" classid="$ClassId" $Style>\n);

    # Setup molecule data...
    $ActiveXTags .= qq(<param name="SourceURL" value="$MolFile">\n<param name="ShrinkToFit" value="$ShrinkToFit">\n<param name="ViewOnly" value="$ViewOnly">\n<param name="ShowToolsWhenVisible" value="$ShowToolsWhenVisible">\n);

    #Setup other parameters...
    for $ParamName (sort keys %ParamsMap) {
      $ParamValue = $ParamsMap{$ParamName};
      if (length $ParamValue) {
	$ActiveXTags .= qq(<param name="$ParamName" value="$ParamValue">\n);
      }
    }
    $ActiveXTags .= qq(</object>\n);
    $ReturnTags = $ActiveXTags;
  }
  return $ReturnTags;
}

# Setup ChemDraw plug-in used for Netscape browsers...
# Problems: "bgcolor" parameter doesn't work.
sub SetupStrViewerChemDrawPlugIn {
  my($MolFile, $Name, $ParamsMapRef, %ParamsMap, $Width, $Height, $ParamName, $ParamValue, $PlugInTags, $JavaScriptTags, $ReturnTags,);
  my($MimeType, $ViewOnly, $ShrinkToFit, $ShowToolsWhenVisible, $JSFileName, $UseJavaScript);

  $Name = "ChemDraw"; $PlugInTags = ""; $ParamsMapRef = ""; %ParamsMap = ();
  $Width = $StrViewerWidth; $Height = $StrViewerHeight;
  $MimeType = "chemical/x-mdl-molfile";
  $ViewOnly = "1";
  $ShrinkToFit = "1";
  $ShowToolsWhenVisible = "1"; $JavaScriptTags = "";

  if (@_ == 2) {
    ($MolFile, $ParamsMapRef) = @_;
  }
  else {
    ($MolFile) = @_;
  }

  if ($ParamsMapRef) {
    %ParamsMap = %$ParamsMapRef;
    if (exists $ParamsMap{usejavascript} ) { $JSFileName = $ParamsMap{usejavascript}; $UseJavaScript = 1; $ParamsMap{usejavascript} = ""; }
    if (exists $ParamsMap{height} ) { $Height = $ParamsMap{height}; $ParamsMap{height} = ""; }
    if (exists $ParamsMap{width} ) { $Width = $ParamsMap{width}; $ParamsMap{width} = ""; }
    if (exists $ParamsMap{type} ) { $MimeType = $ParamsMap{type}; $ParamsMap{type} = ""; }
    if (exists $ParamsMap{viewonly} ) { $ViewOnly = $ParamsMap{viewonly}; $ParamsMap{viewonly} = ""; }
    if (exists $ParamsMap{shrinktofit} ) { $ShrinkToFit = $ParamsMap{shrinktofit}; $ParamsMap{shrinktofit} = ""; }
    if (exists $ParamsMap{showtoolswhenvisible} ) { $ShowToolsWhenVisible = $ParamsMap{showtoolswhenvisible}; $ParamsMap{showtoolswhenvisible} = ""; }
  }
  if ($UseJavaScript) {
    $JavaScriptTags = qq(\n<script>\n);
    $JavaScriptTags .= qq(cd_insertObjectStr("name='$Name' src='$MolFile' type='$MimeType' width='$Width' height='$Height' shrinktofit='$ShrinkToFit' viewonly='$ViewOnly'");\n);
    $JavaScriptTags .= qq(</script>\n);
    $ReturnTags = $JavaScriptTags;
  }
  else {
    # Start plug-in tag...
    $PlugInTags = qq(<embed src="$MolFile" width="$Width" height="$Height" type="$MimeType" viewonly="$ViewOnly" shrinktofit="$ShrinkToFit" showtoolswhenvisible=''$ShowToolsWhenVisible");

    #Setup other parameters...
    for $ParamName (sort keys %ParamsMap) {
      $ParamValue = $ParamsMap{$ParamName};
      if (length $ParamValue) {
	$PlugInTags .= qq(" $ParamName"="$ParamValue");
      }
    }
    # Finish it off...
    $PlugInTags .= qq( >);
    $ReturnTags = $PlugInTags;
  }

  return $ReturnTags;
}


1;

__END__

=head1 NAME

HTMLUtil

=head1 SYNOPSIS

use HTMLUtil;

use HTMLUtil qw(:all);

=head1 DESCRIPTION

B<HTMLUtil> module provides the following functions:

InsertHTMLTags, SetupHTMLAlignmentBegin, SetupHTMLAlignmentEnd,
SetupHTMLButtonRef, SetupHTMLDivBegin, SetupHTMLDivEnd, SetupHTMLEmptyLines,
SetupHTMLHRef, SetupHTMLPageEnd, SetupHTMLPageHeader, SetupHTMLPageTitle,
SetupHTMLStyleSheetTags, SetupHTMLTableColumnEnd, SetupHTMLTableColumnHeader,
SetupHTMLTableEnd, SetupHTMLTableHeader, SetupHTMLTableRowDataValue,
SetupHTMLTableRowEnd, SetupHTMLTableRowHeader, SetupHTMLTableRowHeaderValue,
SetupJavaScriptCmds, SetupStrViewerAccelrysActiveX, SetupStrViewerChem3DActiveX,
SetupStrViewerChemDrawActiveX, SetupStrViewerChemDrawPlugIn,
SetupStrViewerChimePlugIn, SetupStrViewerJMEApplet, SetupStrViewerJSInitCmd,
SetupStrViewerJmolApplet, SetupStrViewerMarvinViewApplet

=head2 FUNCTIONS

=over 4

=item B<InsertHTMLTags>

    $NewTag = InsertHTMLTags($Tag, @TagsNameValue);

Inserts tag name and value pair from I<TagsNameValue> into a exisiting I<Tag> as I<TagName = "TagValue">
and returns B<NewTag> string.

=item B<SetupHTMLAlignmentBegin>

    $AlignmentTag = SetupHTMLAlignmentBegin([$Alignment]);

Returns an alignment begin tag string. Possible I<Alignment> values: I<left, center, or right>.
Default: I<left>.

=item B<SetupHTMLAlignmentEnd>

    $AlignmentTag = SetupHTMLAlignmentBegin([$Alignment]);

Returns an alignment end tag string.

=item B<SetupHTMLButtonRef>

    $ButtonTag = SetupHTMLButtonRef($ButtonLabel, $FileName);

Returns a button tag string for associating B<onClick> button event of a button with label I<ButtonLabel>
to open a file I<FileName>.

=item B<SetupHTMLDivBegin>

    $DivTag = SetupHTMLDivBegin($ID);

Returns a div begin tag string for div I<ID>.

=item B<SetupHTMLDivEnd>

    $DivTag = SetupHTMLDivEnd();

Returns a div end tag string.

=item B<SetupHTMLTableEnd>

    $TableEndTag = SetupHTMLTableEnd();

Returns a table end tag string.

=item B<SetupHTMLEmptyLines>

    $EmptyLineTags = SetupHTMLEmptyLines([$LineCount]);

Returns an empty lines tag string for empty I<LineCount>. Default line count: I<1>.

=item B<SetupHTMLPageHeader>

    $PageHeaderTag = SetupHTMLPageHeader($HeaderTitle, [$Stylesheet,
                     $JavaScript]);

Returns a page header tag string using I<HeaderTitle> and using optionally specifed
values for I<Stylesheet> and I<JavaScript>.

=item B<SetupHTMLHRef>

    $HRef = SetupHTMLHRef($Label, $URL, [$Title]);

Returns a HRef tag string for setting up a URL with I<Label> and I<URL> with optional I<Title>.

=item B<SetupHTMLPageEnd>

    $PageEndTag = SetupHTMLPageEnd([$FooterMsg]);

Returns a page end tag string conating optional I<FooterMsg>.

=item B<SetupHTMLPageTitle>

    $PageTitleTag = SetupHTMLPageTitle($Title, [$Alignment]);

Returns a page title tag string with optional alignment. Valid alignment value: I<left, center, right>
Default alignment: I<center>.

=item B<SetupHTMLStyleSheetTags>

    $StyleSheetTags = SetupHTMLStyleSheetTags();

Returns a default style sheet tag string to be used for HTML files generated by MayaChemTools.

=item B<SetupHTMLTableHeader>

    $TableHeaderTags = SetupHTMLTableHeader([$BorderWidth,
                       $CellPadding, $CellSpacing, $Width, $Height]);

Returns a table header tag string containing specified values for I<BorderWidth, CellPadding, CellSpacing,
Width, and Height>. Default values: I<BorderWidth = 1; CellPadding = 2; CellSpacing = 0; Width = NotUsed;
Height = NotUsed>.

=item <SetupHTMLTableEnd>

    $TableEndTag = SetupHTMLTableEnd();

Returns a table end tag string.

=item B<SetupHTMLTableColumnHeader>

    $ColumnHeaderTag = SetupHTMLTableColumnHeader([$BgColor, $Width]);

Returns a table column header tag string containing specified values for I<BgColor, Width>. Default
values: I<BgColor = NotUsed; Width = NotUsed>.

=item B<SetupHTMLTableColumnEnd>

    $ColumnEndTag = SetupHTMLTableColumnEnd();

Returns a table column end tag string.

=item B<SetupHTMLTableRowHeader>

    $RowHeaderTag = SetupHTMLTableRowHeader([$HAlignment, $BgColor,
                    $VAlignment]);

Returns a table row header tag string containing specified values for I<HAlignment, BgColor, and VAlignment>.
Default values: I<HAlignment = center; $BgColor = NotUsed; $VAlignment = top>.

=item B<SetupHTMLTableRowEnd>

    $RowEndTag = SetupHTMLTableRowEnd();

Returns a table row end tag string.

=item B<SetupHTMLTableRowHeaderValue>

    $HeaderValueTag = SetupHTMLTableRowHeaderValue([$Value]);

Returns a table header row tag string using specifed I<Value>. Default value: I<EmptySpace>.

=item B<SetupHTMLTableRowDataValue>

    $RowValueTag = SetupHTMLTableRowDataValue([$Value, $BgColor,
                   $FontColor, $FontBold]);

Returns a table row column value tag string using specified values for I<Value, BgColor,
FontColor, and FontBold>. Default values: I<Value = EmptySpace; BgColor = NotUsed;
FontColor = NotUsed; $FontBold = NotUsed>.

=item B<SetupJavaScriptCmds>

   $JSTag = SetupJavaScriptCmds(@JSCmdList);

Returns a Java script tag string using java script commands specified in I<JSCmdList>.

=item B<SetupStrViewerJSInitCmd>

   $JSTag = SetupStrViewerJSInitCmd($StrViewerType, $CodeBase);

Returns a Java script command tag string for intializing structure viewers with specified I<CodeBase>
location for viewers to be invoked as Java Applets. Supported values for I<StrViewerType>: I<Jmol,
ChemDrawPlugIn, ChemDrawActiveX, Chem3DActiveX>.

=item B<SetupStrViewerJMEApplet>

    $JMEAppletTag = SetupStrViewerJMEApplet($MolString, $CodeBase,
                    [{param => "value"}]);

Returns a JME tag string for displaying molecule using I<MolString> along with valid optional applet
parameters specified as name and value pairs. Defaul JME parameter values: I<name = JME; id = JME;
width = 250; height = 170>.

=item B<SetupStrViewerJmolApplet>

    $JmolAppletTag = SetupStrViewerJmolApplet($MolString, $CodeBase,
                     [{param => "value"}]);

Returns a JMol tag string for displaying molecule using I<MolString> along with valid optional applet
parameters specified as name and value pairs. Defaul JMol parameter values: I<name = Jmol; id = Jmol;
width = 250; height = 170; progressbar = true; progresscolor = 0000ff; bgcolor = 000000; JMolScript =
select *; set frank off; wireframe on; spacefill off>.

=item B<SetupStrViewerMarvinViewApplet>

    $MarvinAppletTag = SetupStrViewerMarvinViewApplet($MolString,
                       $CodeBase, [{param => "value"}]);

Returns a MarvinView tag string for displaying molecule using I<MolString> along with valid optional applet
parameters specified as name and value pairs. Defaul MarvinView parameter values: I<name = MView; id = MView;
width = 250; height = 170; navmode = zoom>.

=item B<SetupStrViewerChimePlugIn>

    $ChimePlugInTag = SetupStrViewerChimePlugIn($MolFile,
                      [{param => "value"}]);

Returns a MDL Chime tag string for displaying molecule using I<MolFile> along with valid optional
parameters specified as name and value pairs. Defaul Chime parameter values: I<width = 250; height = 170;
display2d = true>.

=item B<SetupStrViewerChem3DActiveX>

    $ChemDraw3DActiveXTags = SetupStrViewerChemDrawActiveX($MolFile,
                             [{param => "value"}]);

Returns a CambridgeSoft Chem3D tag string for displaying molecule using I<MolFile> along with valid optional
parameters specified as name and value pairs. Defaul Chime parameter values: I<width = 250; height = 170;
displaytype = BallAndStick; rotationbars = false; moviecontroller = false>.

=item B<SetupStrViewerChemDrawActiveX>

    $ChemDrawActiveXTags = SetupStrViewerChem3DActiveX($MolFile,
                           [{param => "value"}]);

Returns a CambridgeSoft ChemDraw ActiveX tag string for displaying molecule using I<MolFile> along with valid optional
parameters specified as name and value pairs. Defaul ChemDraw ActiveX parameter values: I<width = 250; height = 170;
ViewOnly = 1; ShrinkToFit = 1; ShowToolsWhenVisible = 1>.

=item B<SetupStrViewerChemDrawPlugIn>

    $ChemDrawPlugInTag = SetupStrViewerChemDrawPlugIn($MolFile,
                         [{param => "value"}]);

Returns a CambridgeSoft ChemDraw PlugIn tag string for displaying molecule using I<MolFile> along with valid optional
parameters specified as name and value pairs. Defaul ChemDraw PlugIn parameter values: I<width = 250; height = 170;
ViewOnly = 1; ShrinkToFit = 1; ShowToolsWhenVisible = 1>.

=item B<SetupStrViewerAccelrysActiveX>

    $AccelrysActiveXTags = SetupStrViewerAccelrysActiveX($MolFile,
                           [{param => "value"}]);

Returns a Accelrys ViewerActiveX tag string for displaying molecule using I<MolFile> along with valid optional
parameters specified as name and value pairs. Defaul ViewerActiveX parameter values: I<width = 250; height = 170;
Convert2Dto3D = 0; Mouse = 4>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
