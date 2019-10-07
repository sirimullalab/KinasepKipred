#!/bin/env python
#
# File: RDKitDrawMoleculesAndDataTable.py
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2018 Manish Sud. All rights reserved.
#
# The functionality available in this script is implemented using RDKit, an
# open source toolkit for cheminformatics developed by Greg Landrum.
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

from __future__ import print_function

# Add local python path to the global path and import standard library modules...
import os
import sys;  sys.path.insert(0, os.path.join(os.path.dirname(sys.argv[0]), "..", "lib", "Python"))
import time
import re
import random

# RDKit imports...
try:
    from rdkit import rdBase
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError as ErrMsg:
    sys.stderr.write("\nFailed to import RDKit module/package: %s\n" % ErrMsg)
    sys.stderr.write("Check/update your RDKit environment and try again.\n\n")
    sys.exit(1)

# MayaChemTools imports...
try:
    from docopt import docopt
    import MiscUtil
    import RDKitUtil
except ImportError as ErrMsg:
    sys.stderr.write("\nFailed to import MayaChemTools module/package: %s\n" % ErrMsg)
    sys.stderr.write("Check/update your MayaChemTools environment and try again.\n\n")
    sys.exit(1)

ScriptName = os.path.basename(sys.argv[0])
Options = {}
OptionsInfo = {}

def main():
    """Start execution of the script"""
    
    MiscUtil.PrintInfo("\n%s (RDK v%s; %s): Starting...\n" % (ScriptName, rdBase.rdkitVersion, time.asctime()))
    
    (WallClockTime, ProcessorTime) = MiscUtil.GetWallClockAndProcessorTime()
    
    # Retrieve command line arguments and options...
    RetrieveOptions()
    
    # Process and validate command line arguments and options...
    ProcessOptions()
    
    # Perform actions required by the script...
    GenerateMoleculesAndDataTable()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def GenerateMoleculesAndDataTable():
    """Generate a HTML table containing molecules and alphanumerical data."""
    
    # Retrieve data...
    ValidMols = RetrieveMoleculesAndData()

    # Setup data type map...
    DataMap = IdentifyStructureAndNumericalData(ValidMols)

    # Validate data labels used to specify highlighting data...
    ValidateSpecifiedDataLabels(DataMap)

    # Validate show molecule name option...
    ValidateShowMolNameOption(DataMap)
    
    # Compute 2D coordinates before alignment...
    if OptionsInfo["Compute2DCoords"]:
        MiscUtil.PrintInfo("\nComputing 2D coordinates for primary structure data...")
        for Mol in ValidMols:
            AllChem.Compute2DCoords(Mol)
            
    # Perform alignment to a common template for primary molecular structure data...
    PerformAlignment(ValidMols)

    # Write out a HTML file...
    WriteHTMLTableFile(ValidMols, DataMap)

def WriteHTMLTableFile(ValidMols, DataMap):
    """Write out a HTML table file."""
    
    Outfile = OptionsInfo["Outfile"]
    
    Writer = open(Outfile, "w")
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
    
    MiscUtil.PrintInfo("\nGenerating file %s..." % Outfile)

    WriteHTMLPageHeader(Writer, DataMap)
    WriteHTMLPageTitle(Writer)
    
    WriteHTMLTableHeader(Writer)
    WriteHTMLTableRows(Writer, ValidMols, DataMap)
    WriteHTMLTableEnd(Writer)
    
    WriteHTMLPageFooter(Writer)
    WriteHTMLPageEnd(Writer)
    
    if Writer is not None:
        Writer.close()

def WriteHTMLTableRows(Writer, ValidMols, DataMap):
    """Write out HTML table rows."""

    WriteTableHeaderRow(Writer, ValidMols, DataMap)
    WriteTableDataRows(Writer, ValidMols, DataMap)
    WriteTableFooterRow(Writer, ValidMols, DataMap)

def WriteTableDataRows(Writer, ValidMols, DataMap):
    """Write out table data row."""

    Writer.write("""        <tbody>\n""")
    
    MolCount = 0
    for Mol in ValidMols:
        MolCount += 1
        Writer.write("""          <tr>\n""")

        if OptionsInfo["CounterCol"]:
            Writer.write("""            <td></td>\n""")
            
        SetupPrimaryStructureTableData(Writer, Mol)
        
        if OptionsInfo["ShowMolName"]:
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            WrappedMolName = MiscUtil.WrapText(MolName, "<br/>", OptionsInfo["WrapTextWidth"])
            Writer.write("""            <td>%s</td>\n""" % WrappedMolName)

        # Set up rest of the data..
        AvailableDataLabelsMap = Mol.GetPropsAsDict(includePrivate = False, includeComputed = False)
        for DataLabel in DataMap["DataLabels"]:
            if not DataLabel in AvailableDataLabelsMap:
                Writer.write("""            <td></td>\n""")
                continue

            # Check for empty value...
            DataValue = "%s" % AvailableDataLabelsMap[DataLabel]
            DataValue = DataValue.strip()
            if MiscUtil.IsEmpty(DataValue):
                Writer.write("""            <td></td>\n""")
                continue
            
            if DataMap["StructureDataMap"][DataLabel]:
                SetupNonPrimaryStructureTableData(Writer, DataLabel, DataValue, DataMap)
            else:
                SetupAlphanumericTableData(Writer, DataLabel, DataValue, DataMap)
        
        Writer.write("""          </tr>\n""")
        
    Writer.write("""        </tbody>\n""")

def SetupPrimaryStructureTableData(Writer, Mol):
    """Set up an inline SVG image for primary structure data for a table cell."""
    
    HightlightAtomList = SetupAtomListToHighlight(Mol, "Structure")
    SVGImageTag = SetupMolInLineSVGImageTag(Mol, HightlightAtomList)
    
    Writer.write("""            <td bgcolor="white"><%s></td>\n""" % SVGImageTag)

def SetupNonPrimaryStructureTableData(Writer, DataLabel, DataValue, DataMap):
    """Set up an inline SVG image for non primary structure data cell."""

    WrappedDataValue = DataValue
    if OptionsInfo["WrapText"]:
        WrappedDataValue = MiscUtil.WrapText(DataValue, "<br/>", OptionsInfo["WrapTextWidth"])

    if DataMap["SMILESDataMap"][DataLabel]:
        Mol = Chem.MolFromSmiles(DataValue, sanitize = False)
        Mol.UpdatePropertyCache(strict = False)
    else:
        MiscUtil.PrintWarning("\nIgnoring uknown structure data column type with column label %s: %s\n" % (DataLabel, DataValue))
        Writer.write("""            <td>%s</td>\n""" % WrappedDataValue)
        return

    if Mol is None:
        MiscUtil.PrintWarning("\nSMILES parsing failed for data label %s: %s\n" % (DataLabel, DataValue))
        Writer.write("""            <td>%s</td>\n""" % WrappedDataValue)
        return
    elif  not Mol.GetNumHeavyAtoms():
        Writer.write("""            <td>%s</td>\n""" % WrappedDataValue)
        return
    elif AllChem.Compute2DCoords(Mol) < 0:
        Writer.write("""            <td>%s</td>\n""" % WrappedDataValue)
        return
    
    HightlightAtomList = SetupAtomListToHighlight(Mol, DataLabel)
    SVGImageTag = SetupMolInLineSVGImageTag(Mol, HightlightAtomList)
    
    Writer.write("""            <td bgcolor="white"><%s></td>\n""" % SVGImageTag)

def SetupAlphanumericTableData(Writer, DataLabel, DataValue, DataMap):
    """Set up alphanumeric data."""

    BackgroundColor, BackgroundColorType = GetAlphanumeircValueHighlightBackgroundColor(DataLabel, DataValue, DataMap)
    SetupAlphanumericTableDataValue(Writer, DataValue, BackgroundColor, BackgroundColorType)

def WriteTableHeaderRow(Writer, ValidMols, DataMap):
    """Write out table header row."""

    TableHeaderStyle = OptionsInfo["TableHeaderStyle"]
    if TableHeaderStyle is None:
        Writer.write("""      <thead>\n""")
        Writer.write("""        <tr>\n""")
    elif re.match("^(thead|table)", TableHeaderStyle):
        Writer.write("""      <thead class="%s">\n""" % TableHeaderStyle)
        Writer.write("""        <tr>\n""")
    else:
        Writer.write("""      <thead>\n""")
        Writer.write("""        <tr bgcolor="%s"\n""" % TableHeaderStyle)

    if OptionsInfo["CounterCol"]:
        Writer.write("""          <th></th>\n""")
    Writer.write("""          <th>Structure</th>\n""")
    if OptionsInfo["ShowMolName"]:
        Writer.write("""          <th>%s</th>\n""" % OptionsInfo["ShowMolNameDataLabel"])
    
    # Write out rest of the column headers...
    for DataLabel in DataMap["DataLabels"]:
        Writer.write("""          <th>%s</th>\n""" % DataLabel)
        
    Writer.write("""        </tr>\n""")
    Writer.write("""      </thead>\n""")

def WriteTableFooterRow(Writer, ValidMols, DataMap):
    """Write out table footer row."""

    if not OptionsInfo["TableFooter"]:
        return
    
    Writer.write("""      <tfoot>\n""")
    Writer.write("""        <tr>\n""")

    if OptionsInfo["CounterCol"]:
        Writer.write("""          <td></td>\n""")
    Writer.write("""          <td>Structure</td>\n""")
    if OptionsInfo["ShowMolName"]:
        Writer.write("""          <td>%s</td>\n""" % OptionsInfo["ShowMolNameDataLabel"])

    # Write out rest of the column headers...
    for DataLabel in DataMap["DataLabels"]:
        Writer.write("""          <td>%s</td>\n""" % DataLabel)
        
    Writer.write("""        </tr>\n""")
    Writer.write("""      </tfoot>\n""")

def WriteHTMLPageHeader(Writer, DataMap):
    """Write out HTML page header."""

    # Collect column indices containing counter and structure data to disable
    # sorting and searching. In addition, set up a list to exclude counter and
    # primary structure columns from column visibility pulldown along with
    # any other columns...
    #
    if OptionsInfo["CounterCol"]:
        StrColIndicesList = ["0", "1"]
        ColVisibilityExcludeColIndicesList = ["0", "1"]
        ColIndexOffset = 2
        FreezeLeftColumns = "2"
    else:
        StrColIndicesList = ["0"]
        ColVisibilityExcludeColIndicesList = ["0"]
        ColIndexOffset = 1
        FreezeLeftColumns = "1"

    if OptionsInfo["ShowMolName"]:
        ColIndexOffset += 1
        
    MaxColVisColCount = OptionsInfo["ColVisibilityCtrlMax"]
    MaxDataColVisColCount = MaxColVisColCount - len(ColVisibilityExcludeColIndicesList)
    MaxDataColVisColCount = MaxColVisColCount

    DataColVisibilityExclude = False
    ColCount = len(DataMap["DataLabels"])
    if OptionsInfo["ColVisibility"]:
        if ColCount > MaxDataColVisColCount:
            DataColVisibilityExclude = True
            MiscUtil.PrintWarning("The number of data columns, %d, is more than %d. Only first %d data columns will be available in column visibility pulldown." % (ColCount, MaxColVisColCount, MaxColVisColCount))
    
    DisplayButtons = False
    if OptionsInfo["ColVisibility"]:
        if ColCount > 0 or OptionsInfo["ShowMolName"]:
            DisplayButtons = True
    
    FreezeCols = False
    if (OptionsInfo["FreezeCols"] and OptionsInfo["ScrollX"]):
        FreezeCols = True
    
    for Index, DataLabel in enumerate(DataMap["DataLabels"]):
        if DataMap["StructureDataMap"][DataLabel]:
            StrColIndex = Index + ColIndexOffset
            StrColIndicesList.append("%s" % StrColIndex)
        
        if OptionsInfo["ColVisibility"]:
            if Index >= MaxDataColVisColCount:
                ColIndex = Index + ColIndexOffset
                ColVisibilityExcludeColIndicesList.append("%s" %ColIndex)
            
    StrColIndices = MiscUtil.JoinWords(StrColIndicesList, ", ")
    ColVisibilityExcludeColIndices = MiscUtil.JoinWords(ColVisibilityExcludeColIndicesList, ", ")

    Paging = "true" if OptionsInfo["Paging"] else "false"
    PageLength = "%d" % OptionsInfo["PageLength"]
    PagingType = "\"%s\"" % OptionsInfo["PagingType"]

    ScrollX = "true" if OptionsInfo["ScrollX"] else "false"
    
    ScrollY = ""
    if OptionsInfo["ScrollY"]:
        if re.search("vh$", OptionsInfo["ScrollYSize"]):
            ScrollY = "\"%s\"" % OptionsInfo["ScrollYSize"]
        else:
            ScrollY = "%s" % OptionsInfo["ScrollYSize"]
    
    RegexSearch = "true" if  OptionsInfo["RegexSearch"] else "false"

    # Start HTML header...
    Title = "Molecules and data table" if OptionsInfo["Header"] is None else OptionsInfo["Header"]

    Writer.write("""\
<!doctype html>
<html lang="en">
<head>
    <title>%s</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css">
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.16/css/dataTables.bootstrap4.min.css">
  
""" % (Title))

    if (FreezeCols):
        Writer.write("""\
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/fixedcolumns/3.2.4/css/fixedColumns.bootstrap4.min.css">
""")
    
    if (OptionsInfo["KeysNavigation"]):
        Writer.write("""\
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/keytable/2.3.2/css/keyTable.bootstrap4.min.css">
""")
    
    Writer.write("""\

    <script type="text/javascript" language="javascript" src="https://code.jquery.com/jquery-1.12.4.js"></script>
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"></script>
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/1.10.16/js/dataTables.bootstrap4.min.js"></script>

""")

    if DisplayButtons:
        Writer.write("""\
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/buttons/1.5.1/js/dataTables.buttons.min.js"></script>
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/buttons/1.5.1/js/buttons.bootstrap4.min.js"></script>
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/buttons/1.5.1/js/buttons.colVis.min.js"></script>

""")
    
    if (FreezeCols):
        Writer.write("""\
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/fixedcolumns/3.2.4/js/dataTables.fixedColumns.min.js"></script>
""")
    
    if (OptionsInfo["KeysNavigation"]):
        Writer.write("""\
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/keytable/2.3.2/js/dataTables.keyTable.min.js"></script>
""")
    
    # Intialize table using Bootstrap, DataTables and JQuery frameworks...
    Writer.write("""\

    <script type="text/javascript" class="init">

$(document).ready(function() {
    var MolsAndDataTable = $('#MolsAndDataTable').DataTable( {
        "columnDefs": [
            {
                "orderable": false,
                "searchable": false,
                "targets": [%s]
            },
""" % (StrColIndices))

    if OptionsInfo["ColVisibility"]:
        Writer.write("""\
            {
                "className": "noColVisCtrl",
                "targets": [%s]
            }
""" % (ColVisibilityExcludeColIndices))

    Writer.write("""\
        ],
""")
    
    # Setup column visibility control pulldown by excluding counter column
    # and primary structure column from the list...
    #
    if OptionsInfo["ColVisibility"]:
        # Set up dom for button display...
        if OptionsInfo["Paging"]:
            Writer.write("""\
        "dom":  "<'row'<'col'l><'col'B><'col'f>>" +
            "<'row'<'col-sm-12'tr>>" +
            "<'row'<'col-sm-5'i><'col-sm-7'p>>",
""")
        else:
            Writer.write("""\
        "dom":  "<'row'<'col-sm-6'B><'col-sm-6'f>>" +
            "<'row'<'col-sm-12'tr>>" +
            "<'row'<'col-sm-5'i><'col-sm-7'p>>",
""")
        # Set up buttons...
        Writer.write("""\
        "buttons": [
            {
                "extend": "colvis",
                "text": "Column visibility",
                "className": "btn btn-outline-light text-dark",
                "columns": ":not(.noColVisCtrl)",
""")
        if not DataColVisibilityExclude:
            Writer.write("""\
                "prefixButtons": [ "colvisRestore" ],
""")
        
        Writer.write("""\
                "columnText": function ( dt, colIndex, colLabel ) {
                    return (colIndex + 1) + ": " + colLabel;
                },
            }
        ],
""")

    # Write out rest of the variables for DataTables...
    if FreezeCols:
        Writer.write("""\
        "fixedColumns": {
            "leftColumns": %s
        },
""" % (FreezeLeftColumns))
    
    if (OptionsInfo["KeysNavigation"]):
        Writer.write("""\
        "keys": true,
""")
    
    Writer.write("""\
        "pageLength": %s,
        "lengthMenu": [ [10, 15, 25, 50, 100, 500, 1000, -1], [10, 15, 25, 50, 100, 500, 1000, "All"] ],
        "paging": %s,
        "pagingType": %s,
        "scrollX": %s,
        "scrollY": %s,
        "scrollCollapse": true,
        "order": [],
        "search" : {"regex" : %s},
    } );
""" % (PageLength, Paging, PagingType, ScrollX, ScrollY, RegexSearch))
    
    if OptionsInfo["CounterCol"]:
        Writer.write("""\
    MolsAndDataTable.on( 'order.dt search.dt', function () {
        MolsAndDataTable.column(0, {search:'applied', order:'applied'}).nodes().each( function (cell, rowIndex) {
            cell.innerHTML = rowIndex + 1;
        } );
    } ).draw();
""")
    
    # End of Javacscript code...
    Writer.write("""\
} );

    </script>
""")

    # Finish up HTML header...
    Writer.write("""\
  
</head>
<body>
  <div class="container-fluid">
    <br/>
""")

def WriteHTMLPageEnd(Writer):
    """Write out HTML page end."""

    Writer.write("""\
  </div>
</body>
</html>
""")

def WriteHTMLPageTitle(Writer):
    """Write out HTML page title."""

    if OptionsInfo["Header"] is None:
        return

    Writer.write("""    <%s class="text-center">%s</%s>\n""" % (OptionsInfo["HeaderStyle"], OptionsInfo["Header"], OptionsInfo["HeaderStyle"]))

def WriteHTMLPageFooter(Writer):
    """Write out HTML page footer."""

    if OptionsInfo["Footer"] is None:
        return

    Writer.write("""    <br/>\n    <p class="%s">%s</p>\n""" % (OptionsInfo["FooterClass"], OptionsInfo["Footer"]))

def WriteHTMLTableHeader(Writer):
    """Write out HTML table header."""

    if OptionsInfo["TableStyle"] is None:
        Writer.write("""\n    <table id="MolsAndDataTable" cellspacing="0" width="100%">\n""")
    else:
        Writer.write("""    <table id="MolsAndDataTable" class="%s" cellspacing="0" width="100%s">\n""" % (OptionsInfo["TableStyle"], "%"))
        
def WriteHTMLTableEnd(Writer):
    """Write out HTML table end."""

    Writer.write("""    </table>\n\n""")

def RetrieveMoleculesAndData():
    """Retrieve molecules and data from input file."""

    MiscUtil.PrintInfo("\nReading file %s..." % OptionsInfo["Infile"])

    if MiscUtil.CheckFileExt(OptionsInfo["Infile"] ,"smi csv tsv txt"):
        # Check for the presence of SMILES column name in title line...
        Infile = open(OptionsInfo["Infile"], "r")
        if Infile is None:
            MiscUtil.PrintError("Couldn't open file %s..." % OptionsInfo["Infile"])
        Line = Infile.readline()
        Infile.close()
        
        if not re.search("SMILES", Line, re.I):
            MiscUtil.PrintError("The input file, %s, must contain a title line containing a column name with SMILES in its name." % OptionsInfo["Infile"])
    
    if MiscUtil.CheckFileExt(OptionsInfo["Infile"],"sdf sd smi"):
        ValidMols, MolCount, ValidMolCount  = RDKitUtil.ReadAndValidateMolecules(OptionsInfo["Infile"], **OptionsInfo["InfileParams"])
    else:
        ValidMols, MolCount, ValidMolCount  = RetrieveMoleculesFromTextFile(OptionsInfo["Infile"])
    
    MiscUtil.PrintInfo("Total number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))

    return ValidMols

def RetrieveMoleculesFromTextFile(Infile):
    """Retrieve molecules from a CSV/TSV text file."""

    # Read and parse text lines...
    Delimiter = "," if MiscUtil.CheckFileExt(Infile ,"csv") else "\t"
    QuoteChar = '"'
    IgnoreHeaderLine = False
    TextLinesWords = MiscUtil.GetTextLinesWords(Infile, Delimiter, QuoteChar, IgnoreHeaderLine)

    # Process column names...
    ColNames = TextLinesWords[0]
    ColCount = len(ColNames)
    
    MolColIndex = None
    MolDataColIndices = []

    FirstSMILES = True
    for ColIndex in range(0, ColCount):
        if re.search("SMILES", ColNames[ColIndex], re.I) and FirstSMILES:
            MolColIndex = ColIndex
            FirstSMILES = False
            continue
        
        MolDataColIndices.append(ColIndex)
    
    if MolColIndex is None:
            MiscUtil.PrintError("The input file, %s, must contain a title line containing a column name with SMILES in its name." % Infile)

    ValidMols = []
    MolCount = 0

    Sanitize = OptionsInfo["InfileParams"]["Sanitize"]
    
    # Process data lines...
    for LineIndex in range(1, len(TextLinesWords)):
        MolCount += 1
        LineWords = TextLinesWords[LineIndex]
        if len(LineWords) != ColCount:
            MiscUtil.PrintWarning("Ignoring text line number %d: Number of columns, %d, must match number of columns, %d, in title line.\nLine: %s" % (MolCount, len(LineWords), ColCount, Delimiter.join(LineWords)))
            continue

        # Process molecule column...
        MolSMILES = LineWords[MolColIndex]
        Mol = Chem.MolFromSmiles(MolSMILES, sanitize = Sanitize)
        if Mol is None:
            MiscUtil.PrintWarning("Ignoring text line number %d: SMILES parsing failed\nLine: %s" % (MolCount, Delimiter.join(LineWords)))
            continue

        # Process molecule data columns...
        for ColIndex in MolDataColIndices:
            Name = ColNames[ColIndex]
            Value = LineWords[ColIndex]
            Mol.SetProp(Name, Value)
        
        ValidMols.append(Mol)
    
    ValidMolCount = len(ValidMols)
    
    return (ValidMols, MolCount, ValidMolCount)

def IdentifyStructureAndNumericalData(ValidMols):
    """Identify structure and alphanumerical data."""

    DataMap = {}
    DataMap["DataLabels"] = []
    DataMap["DataLabelsMap"] = {}
    DataMap["CanonicalDataLabelsMap"] = {}
    
    DataMap["StructureDataMap"] = {}
    DataMap["SMILESDataMap"] = {}

    # Retrieve all possible data labels...
    if MiscUtil.CheckFileExt(OptionsInfo["Infile"] ,"smi csv tsv txt"):
        # First molecule contains all possible data fields...
        Mol = ValidMols[0]
        ProcessMolDataLabels(ValidMols[0], DataMap)
    else:
        # Go over all molecules to identify unique data labels...
        MiscUtil.PrintInfo("\nRetrieving unique data labels for data in file %s..." % OptionsInfo["Infile"])
        for Mol in ValidMols:
            ProcessMolDataLabels(Mol, DataMap)

    return DataMap

def ProcessMolDataLabels(Mol, DataMap):
    """Process data label to identify and track its type"""
    
    for DataLabel in Mol.GetPropNames(includePrivate = False, includeComputed = False):
        if DataLabel in DataMap["DataLabelsMap"]:
            continue

        # Track labels...
        DataMap["DataLabels"].append(DataLabel)
        DataMap["DataLabelsMap"][DataLabel] = DataLabel
        DataMap["CanonicalDataLabelsMap"][DataLabel.lower()] = DataLabel
        
        DataMap["StructureDataMap"][DataLabel] = False
        DataMap["SMILESDataMap"][DataLabel] = False
        
        if re.search("SMILES", DataLabel, re.I):
            DataMap["StructureDataMap"][DataLabel] = True
            DataMap["SMILESDataMap"][DataLabel] = True

def ValidateShowMolNameOption(DataMap):
    """Validate show molecule name option. """
    
    if not OptionsInfo["ShowMolName"]:
        return

    if not MiscUtil.CheckFileExt(OptionsInfo["Infile"],"sdf sd smi"):
        OptionsInfo["ShowMolName"] = False
        return

    CanonicalDataLabel = OptionsInfo["ShowMolNameDataLabel"].lower()
    if CanonicalDataLabel in DataMap["CanonicalDataLabelsMap"]:
        OptionsInfo["ShowMolName"] = False
        if not OptionsInfo["ShowMolNameAuto"]:
            MiscUtil.PrintWarning("Ignoring \"--showMolName\" option: Data label \"Name\" corresponding to molecule name is already present in input file.")

def ValidateSpecifiedDataLabels(DataMap):
    """Validate data labels used to specify highlighting data. """
    
    ValidateSpecifiedDataLabelsForHighlightSMARTS(DataMap)
    
    ValidateSpecifiedDataLabelsForHighlightValues(DataMap)
    ValidateSpecifiedDataLabelsForHighlightRanges(DataMap)
    ValidateSpecifiedDataLabelsForHighlightClasses(DataMap)
    
def ValidateSpecifiedDataLabelsForHighlightSMARTS(DataMap):
    """Validate data labels used to specify highlighting SMARTS option. """

    if OptionsInfo["HighlightSMARTSAllMode"]:
        return

    for DataLabel in OptionsInfo["HighlightSMARTSDataLabels"]:
        if re.match("^Structure$", DataLabel, re.I):
            continue
        
        CanonicalDataLabel = DataLabel.lower()
        if not CanonicalDataLabel in DataMap["CanonicalDataLabelsMap"]:
            MiscUtil.PrintError("The data label specified, %s, using option \"--highlightSMARTS\" doesn't exist in input file." % DataLabel)
        
        Label = DataMap["CanonicalDataLabelsMap"][CanonicalDataLabel]
        if not DataMap["StructureDataMap"][Label]:
            MiscUtil.PrintError("The data label specified, %s, using option \"--highlightSMARTS\" doesn't correspond to structure data: Valid structure data labels: SMILES in data label." % DataLabel)
            
def ValidateSpecifiedDataLabelsForHighlightValues(DataMap):
    """Validate data labels used to specify highlighting values option. """

    ValidateDataLabels("--highlightValues", DataMap, OptionsInfo["HighlightValuesLabels"])

def ValidateSpecifiedDataLabelsForHighlightRanges(DataMap):
    """Validate data labels used to specify highlighting ranges option. """
    
    ValidateDataLabels("--highlightRanges", DataMap, OptionsInfo["HighlightRangesLabels"])

def ValidateSpecifiedDataLabelsForHighlightClasses(DataMap):
    """Validate data labels used to specify highlighting classes option. """

    if OptionsInfo["HighlightClassesRules"] is None:
        return
    
    ValidDataLabelsList = []
    NotValidDataLabelsList = []
    for Label in OptionsInfo["HighlightClassesLabels"]:
        ValidCanonicalLabel = None
        
        for LabelSynonym in OptionsInfo["HighlightClassesSynonymsMap"][Label]:
            CanonicalLabel = LabelSynonym.lower()
            
            # Is this label already in use...
            if CanonicalLabel in OptionsInfo["HighlightValuesCanonicalLabelsMap"]:
                MiscUtil.PrintInfo("")
                MiscUtil.PrintWarning("The data label, %s, for class, %s , in option \"--highlightValuesClasses\" has already been used in \"--highlightValues\" option. It'll be ignored during highlighting." % (LabelSynonym, OptionsInfo["HighlightClasses"]))
                continue

            if CanonicalLabel in OptionsInfo["HighlightRangesCanonicalLabelsMap"]:
                MiscUtil.PrintInfo("")
                MiscUtil.PrintWarning("The data label, %s, for class, %s , in option \"--highlightValuesClasses\" has already been used in \"--highlightValuesRanges\" option. It'll be ignored during highlighting." % (LabelSynonym, OptionsInfo["HighlightClasses"]))
                continue
            
            # Is this label present in input file...
            if  CanonicalLabel in DataMap["CanonicalDataLabelsMap"]:
                ValidCanonicalLabel = CanonicalLabel
                break
        
        if ValidCanonicalLabel is None:
            MiscUtil.PrintWarning("The data label or its synonyms - %s - for class, %s , in option \"--highlightValuesClasses\" either don't exist in input file or have already been used for highlighting in option \"--highlightValuesClasses\" or \"--highlightValuesRanges\". It'll be ignored during highlighting." % (MiscUtil.JoinWords(OptionsInfo["HighlightClassesSynonymsMap"][Label], ", "), OptionsInfo["HighlightClasses"]))
            NotValidDataLabelsList.append(Label)
            continue
        
        # Track label...
        OptionsInfo["HighlightClassesCanonicalLabelsMap"][ValidCanonicalLabel] = Label
        ValidDataLabelsList.append(DataMap["CanonicalDataLabelsMap"][ValidCanonicalLabel])

    ValidDataLabelsCount = len(ValidDataLabelsList)
    NotValidDataLabelsCount = len(NotValidDataLabelsList)
    DataLabelsCount = len(OptionsInfo["HighlightClassesLabels"])
    
    if ValidDataLabelsCount == 0:
        MiscUtil.PrintInfo("")
        MiscUtil.PrintWarning("The data labels and their synonyms for class, %s , in option \"--highlightValuesClasses\" either don't exists in input file or have already been used for highlighting in option \"--highlightValuesClasses\" or \"--highlightValuesRanges\". No class highlighting will be performed. Missing data labels:  %s" % (OptionsInfo["HighlightClasses"], MiscUtil.JoinWords(OptionsInfo["HighlightClassesLabels"], ", ")))
    elif  ValidDataLabelsCount < DataLabelsCount:
        MiscUtil.PrintInfo("")
        MiscUtil.PrintWarning("The class, %s, based highlighting specified using \"--highlightValuesClasses\" option will be performed using only, %d, out of, %d, data labels: %s\nThe rest of the data label(s) - %s - either don't exist in the input file or have aready been used for highlighting in option \"--highlightValuesClasses\" or \"--highlightValuesRanges\"." % (OptionsInfo["HighlightClasses"], ValidDataLabelsCount,  DataLabelsCount,  MiscUtil.JoinWords(ValidDataLabelsList, ", "),  MiscUtil.JoinWords(NotValidDataLabelsList, ", ") ))
        
def ValidateDataLabels(OptionName, DataMap, DataLabels):
    """Validate data labels."""
    
    for DataLabel in DataLabels:
        if re.match("^Structure$", DataLabel, re.I):
            MiscUtil.PrintError("The data label specified, %s, using option \"-%s\" must not correspond to structure data. Structure label is not allowed." % (DataLabel, OptionName))
            
        CanonicalDataLabel = DataLabel.lower()
        if not CanonicalDataLabel in DataMap["CanonicalDataLabelsMap"]:
            MiscUtil.PrintError("The data label specified, %s, using option \"%s\" doesn't exist in input file." % (DataLabel, OptionName))
        
        Label = DataMap["CanonicalDataLabelsMap"][CanonicalDataLabel]
        if DataMap["StructureDataMap"][Label]:
            MiscUtil.PrintError("The data label specified, %s, using option \"%s\" must not correspond to structure data: Valid structure data labels contain \"SMILES\" in their name.." % (DataLabel, OptionName))

def SetupMolInLineSVGImageTag(Mol, HightlightAtomList):
    """Setup a inline SVG image tag for molecule."""
    
    SVGText = RDKitUtil.GetInlineSVGForMolecule(Mol, OptionsInfo["MolImageWidth"], OptionsInfo["MolImageHeight"], AtomListToHighlight = HightlightAtomList)
    SVGInlineImageTag = "img src=\"data:image/svg+xml;charset=UTF-8,\n%s\"" % SVGText
    
    return SVGInlineImageTag

def SetupAtomListToHighlight(Mol, DataLabel):
    """Set up atom list to highlight using specified SMARTS patterns."""

    HighlightAtomList = None
    if OptionsInfo["HighlightSMARTS"] is None:
        return  HighlightAtomList

    if OptionsInfo["HighlightSMARTSAllMode"]:
        PatternMol = OptionsInfo["HighlightSMARTSPatternMol"]
    else:
        CanonicalDataLabel = DataLabel.lower()
        if not CanonicalDataLabel in OptionsInfo["HighlightSMARTSCanonicalDataLabelsMap"]:
            return  HighlightAtomList
        
        Label = OptionsInfo["HighlightSMARTSCanonicalDataLabelsMap"][CanonicalDataLabel]
        PatternMol = OptionsInfo["HighlightSMARTSPatternMolsMap"][Label]
        
    # Get matched atom lists and flatten it...
    MatchedAtomsLists = Mol.GetSubstructMatches(PatternMol)
    MatchedAtoms = [ Atom for AtomsList in MatchedAtomsLists for Atom in AtomsList]

    if len(MatchedAtoms):
        HighlightAtomList = MatchedAtoms
    
    return  HighlightAtomList

def GetAlphanumeircValueHighlightBackgroundColor(DataLabel, DataValue, DataMap):
    """Get background highlight color for a value."""

    BackgroundColor = None
    BackgroundColorType = None
    
    CanonicalDataLabel =DataLabel.lower()
    if CanonicalDataLabel in OptionsInfo["HighlightValuesCanonicalLabelsMap"]:
        return GetBackgroundColorUsingHighlightValuesMode(DataLabel, DataValue, DataMap)
    elif CanonicalDataLabel in OptionsInfo["HighlightRangesCanonicalLabelsMap"]:
        return GetBackgroundColorUsingHighlightRangesMode(DataLabel, DataValue, DataMap)
    elif CanonicalDataLabel in OptionsInfo["HighlightClassesCanonicalLabelsMap"]:
        return GetBackgroundColorUsingHighlightClassesMode(DataLabel, DataValue, DataMap)
    elif OptionsInfo["HighlightClassesRandom"]:
        return GetBackgroundColorUsingRandomMode(DataLabel, DataValue, DataMap)
    
    return (BackgroundColor, BackgroundColorType)
    
def GetBackgroundColorUsingHighlightValuesMode(DataLabel, DataValue, DataMap):
    """Get background highlight color for a value."""

    BackgroundColor = None
    BackgroundColorType = None

    CanonicalDataLabel =DataLabel.lower()
    if not CanonicalDataLabel in OptionsInfo["HighlightValuesCanonicalLabelsMap"]:
        return (BackgroundColor, BackgroundColorType)

    Label = OptionsInfo["HighlightValuesCanonicalLabelsMap"][CanonicalDataLabel]
    DataType = OptionsInfo["HighlightValuesTypesMap"][Label]
    Criterion = OptionsInfo["HighlightValuesCriteriaMap"][Label]
    CriterionValue = OptionsInfo["HighlightValuesCriteriaValuesMap"][Label]

    return GetBackgroundColorForHighlightingValue(DataLabel, DataValue, DataType, Criterion, CriterionValue)
    
def GetBackgroundColorUsingHighlightClassesMode(DataLabel, DataValue, DataMap):
    """Get background highlight color for a value."""

    BackgroundColor = None
    BackgroundColorType = None

    CanonicalDataLabel =DataLabel.lower()
    if not CanonicalDataLabel in OptionsInfo["HighlightClassesCanonicalLabelsMap"]:
        return (BackgroundColor, BackgroundColorType)

    Label = OptionsInfo["HighlightClassesCanonicalLabelsMap"][CanonicalDataLabel]
    DataType = OptionsInfo["HighlightClassesTypesMap"][Label]
    Criterion = OptionsInfo["HighlightClassesCriteriaMap"][Label]
    CriterionValue = OptionsInfo["HighlightClassesCriteriaValuesMap"][Label]

    return GetBackgroundColorForHighlightingValue(DataLabel, DataValue, DataType, Criterion, CriterionValue)

def GetBackgroundColorForHighlightingValue(DataLabel, DataValue, DataType, Criterion, CriterionValue):
    """Get background color for highlighting a value."""
    
    ValueOkay = False
    if re.match("^numeric$", DataType, re.I):
        if not MiscUtil.IsNumber(DataValue):
            MiscUtil.PrintWarning("Ignoring data value, %s, for data label, %s, during numeric highlighting: It must be a number" % (DataValue, DataLabel))
            return (BackgroundColor, BackgroundColorType)
        
        DataValue = float(DataValue)
        if re.match("^gt$", Criterion, re.I):
            ValueOkay = True if DataValue > CriterionValue else False
        elif re.match("^lt$", Criterion, re.I):
            ValueOkay = True if DataValue < CriterionValue else False
        elif re.match("^ge$", Criterion, re.I):
            ValueOkay = True if DataValue >= CriterionValue else False
        elif re.match("^le$", Criterion, re.I):
            ValueOkay = True if DataValue <= CriterionValue else False
        elif re.match("^eq$", Criterion, re.I):
            ValueOkay = True if DataValue == CriterionValue else False
        elif re.match("^ne$", Criterion, re.I):
            ValueOkay = True if DataValue != CriterionValue else False
        else:
            return (BackgroundColor, BackgroundColorType)
    elif re.match("^text$", DataType, re.I):
        DataValue = "%s" % DataValue
        if re.match("^gt$", Criterion, re.I):
            ValueOkay = True if DataValue > CriterionValue else False
        elif re.match("^lt$", Criterion, re.I):
            ValueOkay = True if DataValue < CriterionValue else False
        elif re.match("^ge$", Criterion, re.I):
            ValueOkay = True if DataValue >= CriterionValue else False
        elif re.match("^le$", Criterion, re.I):
            ValueOkay = True if DataValue <= CriterionValue else False
        elif re.match("^eq$", Criterion, re.I):
            ValueOkay = True if DataValue == CriterionValue else False
        elif re.match("^ne$", Criterion, re.I):
            ValueOkay = True if DataValue != CriterionValue else False
        else:
            return (BackgroundColor, BackgroundColorType)
    elif re.match("^regex$", DataType, re.I):
        DataValue = "%s" % DataValue
        if re.match("^eq$", Criterion, re.I):
            ValueOkay = True if re.search("%s" % CriterionValue, DataValue, re.I) else False
        elif re.match("^ne$", Criterion, re.I):
            ValueOkay = False if re.search("%s" % CriterionValue, DataValue, re.I) else True
        else:
            return (BackgroundColor, BackgroundColorType)

    BackgroundColor = OptionsInfo["HighlightColorsList"][0] if ValueOkay else OptionsInfo["HighlightColorsList"][1]
    BackgroundColorType = OptionsInfo["HighlightColorsType"]
    
    return (BackgroundColor, BackgroundColorType)

def GetBackgroundColorUsingHighlightRangesMode(DataLabel, DataValue, DataMap):
    """Get background highlight color for value range."""

    BackgroundColor = None
    BackgroundColorType = None

    CanonicalDataLabel =DataLabel.lower()
    if not CanonicalDataLabel in OptionsInfo["HighlightRangesCanonicalLabelsMap"]:
        return (BackgroundColor, BackgroundColorType)

    Label = OptionsInfo["HighlightRangesCanonicalLabelsMap"][CanonicalDataLabel]
    DataType = OptionsInfo["HighlightRangesTypesMap"][Label]
    CriterionLower = OptionsInfo["HighlightRangesCriteriaLowerMap"][Label]
    CriterionLowerValue = OptionsInfo["HighlightRangesCriteriaLowerValuesMap"][Label]
    CriterionUpper = OptionsInfo["HighlightRangesCriteriaUpperMap"][Label]
    CriterionUpperValue = OptionsInfo["HighlightRangesCriteriaUpperValuesMap"][Label]
    
    if re.match("^numeric$", DataType, re.I):
        if not MiscUtil.IsNumber(DataValue):
            MiscUtil.PrintWarning("Ignoring data value, %s, for data label, %s, during numeric highlighting: It must be a number" % (DataValue, DataLabel))
            return (BackgroundColor, BackgroundColorType)
        
        DataValue = float(DataValue)
        ColorIndex = 1
        
        if DataValue < CriterionLowerValue and re.match("^lt$", CriterionLower, re.I):
            ColorIndex = 0
        elif DataValue <= CriterionLowerValue and re.match("^le$", CriterionLower, re.I):
            ColorIndex = 0
        elif DataValue > CriterionUpperValue and re.match("^gt$", CriterionUpper, re.I):
            ColorIndex = 2
        elif DataValue >= CriterionUpperValue and re.match("^ge$", CriterionUpper, re.I):
            ColorIndex = 2
    elif re.match("^text$", DataType, re.I):
        DataValue = "%s" % DataValue
        ColorIndex = 1
        
        if DataValue < CriterionLowerValue and re.match("^lt$", CriterionLower, re.I):
            ColorIndex = 0
        elif DataValue <= CriterionLowerValue and re.match("^le$", CriterionLower, re.I):
            ColorIndex = 0
        elif DataValue > CriterionUpperValue and re.match("^gt$", CriterionUpper, re.I):
            ColorIndex = 2
        elif DataValue >= CriterionUpperValue and re.match("^ge$", CriterionUpper, re.I):
            ColorIndex = 2
    else:
        return (BackgroundColor, BackgroundColorType)

    BackgroundColor = OptionsInfo["HighlightColorsRangesList"][ColorIndex]
    BackgroundColorType = OptionsInfo["HighlightColorsRangesType"]
    
    return (BackgroundColor, BackgroundColorType)

def GetBackgroundColorUsingRandomMode(DataLabel, DataValue, DataMap):
    """Get a random background highlight color for a value."""

    BackgroundColor = random.choice(OptionsInfo["HighlightColorsRandomList"])
    BackgroundColorType = OptionsInfo["HighlightColorsRandomType"]
    
    return (BackgroundColor, BackgroundColorType)
    
def SetupAlphanumericTableDataValue(Writer, DataValue, BackgroundColor, BackgroundColorType):
    """Set up alphanumeric data value for a table cell."""
    
    WrappedDataValue = "%s" % DataValue

    # Look for new lines...
    Delim = "<br/>"
    if re.search("(\r\n|\r|\n)", WrappedDataValue):
        WrappedDataValue = re.sub("(\r\n|\r|\n)", "<br/>", DataValue)

    # Wrap text...
    if OptionsInfo["WrapText"] and len(WrappedDataValue) > OptionsInfo["WrapTextWidth"]:
        WrappedDataLines = []
        for DataLine in WrappedDataValue.split("<br/>"):
            WrappedDataLine = MiscUtil.WrapText(DataLine, "<br/>", OptionsInfo["WrapTextWidth"])
            WrappedDataLines.append(WrappedDataLine)
        
        WrappedDataValue = "<br/>".join(WrappedDataLines)
    
    # Highlight value...
    if BackgroundColor is not None:
        ColorTypeTag = GetBackgroundColorTypeTagForTableValue(BackgroundColor, BackgroundColorType)
        Writer.write("""            <td %s = "%s">%s</td>\n""" % (ColorTypeTag, BackgroundColor, WrappedDataValue))
    else:
        Writer.write("""            <td>%s</td>\n""" % WrappedDataValue)

def GetBackgroundColorTypeTagForTableValue(Color, ColorType):
    """Setup color type tage for setting background of a table value."""

    ColorTypeTag = "class" if re.match("^colorclass", ColorType, re.I) else "bgcolor"
    
    return ColorTypeTag
    
def PerformAlignment(ValidMols):
    """Perform alignment to a common template specified by a SMARTS pattern."""
    
    if OptionsInfo["AlignmentSMARTSPattern"] is None:
        return
    
    MiscUtil.PrintInfo("\nPerforming alignment for primary structure data...")
    
    PatternMol = Chem.MolFromSmarts(OptionsInfo["AlignmentSMARTSPattern"])
    AllChem.Compute2DCoords(PatternMol)
        
    MatchedValidMols = [ValidMol for ValidMol in ValidMols if ValidMol.HasSubstructMatch(PatternMol)]
    for ValidMol in MatchedValidMols:
        AllChem.GenerateDepictionMatching2DStructure(ValidMol, PatternMol)

def ProcessHighlightSMARTSOption():
    """Process highlight SMARTS option"""

    OptionsInfo["HighlightSMARTS"]  = None
    OptionsInfo["HighlightSMARTSAllMode"]  = False
    OptionsInfo["HighlightSMARTSPatternMol"]  = None
    
    OptionsInfo["HighlightSMARTSDataLabels"]  = []
    OptionsInfo["HighlightSMARTSDataLabelsMap"]  = {}
    
    OptionsInfo["HighlightSMARTSCanonicalDataLabelsMap"]  = {}
    OptionsInfo["HighlightSMARTSPatternsMap"]  = {}
    OptionsInfo["HighlightSMARTSPatternMolsMap"]  = {}

    if re.match("^None$", Options["--highlightSMARTS"], re.I):
        # Nothing to proecess...
        return
    
    HighlightSMARTS = Options["--highlightSMARTS"].strip()
    if not HighlightSMARTS:
        MiscUtil.PrintError("No valid values specified using \"--highlightSMARTS\" option.")

    OptionsInfo["HighlightSMARTS"]  = HighlightSMARTS
    HighlightSMARTSWords = HighlightSMARTS.split(",")

    if len(HighlightSMARTSWords) == 1:
        PatternMol = Chem.MolFromSmarts(HighlightSMARTS)
        if PatternMol is None:
            MiscUtil.PrintError("The value specified, %s, using option \"--highlightSMARTS\" is not a valid SMARTS: Failed to create pattern molecule" % Options["--highlightSMARTS"])
        OptionsInfo["HighlightSMARTSAllMode"]  = True
        OptionsInfo["HighlightSMARTSPatternMol"]  = PatternMol
        return
    
    if len(HighlightSMARTSWords) % 2:
        MiscUtil.PrintError("The number of comma delimited paramater names and values, %d, specified using \"--highlightSMARTS\" option must be an even number." % (len(HighlightSMARTSWords)))
    
    HighlightSMARTSAllMode = False
            
    for Index in range(0, len(HighlightSMARTSWords), 2):
        DataLabel = HighlightSMARTSWords[Index].strip()
        SMARTSPattern = HighlightSMARTSWords[Index + 1].strip()
                
        PatternMol = Chem.MolFromSmarts(SMARTSPattern)
        if PatternMol is None:
            MiscUtil.PrintError("The value specified, %s, using option \"--highlightSMARTS\" is not a valid SMARTS: Failed to create pattern molecule" % Options["--highlightSMARTS"])

        if DataLabel in OptionsInfo["HighlightSMARTSDataLabelsMap"]:
            MiscUtil.PrintError("The datalabel, %s, specified in pair, \"%s, %s\", using option \"--highlightSMARTS\" is not a valid: Multiple occurences of data label" % (DataLabel, DataLabel, SMARTSPattern))
            
        OptionsInfo["HighlightSMARTSDataLabels"].append(DataLabel)
        OptionsInfo["HighlightSMARTSDataLabelsMap"][DataLabel] = DataLabel
        OptionsInfo["HighlightSMARTSCanonicalDataLabelsMap"][DataLabel.lower()] = DataLabel
        OptionsInfo["HighlightSMARTSPatternsMap"][DataLabel] = SMARTSPattern
        OptionsInfo["HighlightSMARTSPatternMolsMap"][DataLabel] = PatternMol

def ProcessHighlightDataOptions():
    """Process highlight values and colors option"""
    
    ProcessHighlightValuesOption()
    ProcessHighlightValuesRangesOption()
    ProcessHighlightValuesClassesOption()
    
    ProcessHighlightColorsOption()
    ProcessHighlightColorsRangesOption()
    ProcessHighlightColorsRandomOption()

def ProcessHighlightValuesOption():
    """Process highlight values option"""

    OptionsInfo["HighlightValues"]  = None
    OptionsInfo["HighlightValuesLabels"]  = []
    
    OptionsInfo["HighlightValuesLabelsMap"]  = {}
    OptionsInfo["HighlightValuesCanonicalLabelsMap"]  = {}
    
    OptionsInfo["HighlightValuesTypesMap"]  = {}
    OptionsInfo["HighlightValuesCriteriaMap"]  = {}
    OptionsInfo["HighlightValuesCriteriaValuesMap"]  = {}

    HighlightValues = Options["--highlightValues"].strip()
    if re.match("^None$", HighlightValues, re.I):
        return
    
    OptionsInfo["HighlightValues"]  = HighlightValues
    HighlightValuesWords = HighlightValues.split(",")
    
    if len(HighlightValuesWords) % 4:
        MiscUtil.PrintError("The number of comma delimited paramater names and values, %d, specified using \"--highlightValues\" option must be a multiple of 4." % (len(HighlightValuesWords)))

    for Index in range(0, len(HighlightValuesWords), 4):
        DataLabel = HighlightValuesWords[Index].strip()
        DataType = HighlightValuesWords[Index + 1].strip()
        DataCriterion = HighlightValuesWords[Index + 2].strip()
        DataValue = HighlightValuesWords[Index + 3].strip()
        
        if not re.match("^(numeric|text|regex)$", DataType, re.I):
            MiscUtil.PrintError("The data type, %s, specified in quratet \"%s,%s,%s,%s\", using \"--highlightValues\" option is not valid. Supported values: numeric, regex or text." % (DataType, DataLabel, DataType, DataCriterion, DataValue))
        
        if re.match("^regex$", DataType, re.I):
            if not re.match("^(eq|ne)$", DataCriterion, re.I):
                MiscUtil.PrintError("The data criterion, %s, specified in quratet \"%s,%s,%s,%s\", using \"--highlightValues\" option is not valid. Supported values: eq or ne" % (DataType, DataLabel, DataType, DataCriterion, DataValue))
        else:
            if not re.match("^(gt|lt|ge|le|eq|ne)$", DataCriterion, re.I):
                MiscUtil.PrintError("The data criterion, %s, specified in quratet \"%s,%s,%s,%s\", using \"--highlightValues\" option is not valid. Supported values: gt, lt, ge, le, eq, or  ne." % (DataType, DataLabel, DataType, DataCriterion, DataValue))

        # Check criterion value...
        if re.match("^numeric$", DataType, re.I):
            if not MiscUtil.IsNumber(DataValue):
                MiscUtil.PrintError("The data value, %s, specified in quratet \"%s,%s,%s,%s\", using \"--highlightValues\" option is not valid. It must be a number for data type, %s" % (DataType, DataLabel, DataType, DataCriterion, DataValue, DataType))
            DataValue = float(DataValue)
        
        # Track values...
        if DataLabel in OptionsInfo["HighlightValuesLabelsMap"]:
            MiscUtil.PrintError("The data label, %s, specified in quratet \"%s,%s,%s,%s\", using \"--highlightValues\" option is not valid: Multiple occurences of data label" % (DataLabel, DataLabel, DataType, DataCriterion, DataValue))

        OptionsInfo["HighlightValuesLabels"].append(DataLabel)
        OptionsInfo["HighlightValuesLabelsMap"][DataLabel]  = DataLabel
        OptionsInfo["HighlightValuesCanonicalLabelsMap"][DataLabel.lower()]  = DataLabel
        
        OptionsInfo["HighlightValuesTypesMap"][DataLabel]  = DataType
        OptionsInfo["HighlightValuesCriteriaMap"][DataLabel]  = DataCriterion 
        OptionsInfo["HighlightValuesCriteriaValuesMap"][DataLabel]  = DataValue

def ProcessHighlightValuesRangesOption():
    """Process highlight values ranges option"""
    
    OptionsInfo["HighlightRanges"]  = None
    OptionsInfo["HighlightRangesLabels"]  = []
    
    OptionsInfo["HighlightRangesLabelsMap"]  = {}
    OptionsInfo["HighlightRangesCanonicalLabelsMap"]  = {}
    
    OptionsInfo["HighlightRangesTypesMap"]  = {}
    OptionsInfo["HighlightRangesCriteriaLowerMap"]  = {}
    OptionsInfo["HighlightRangesCriteriaLowerValuesMap"]  = {}
    OptionsInfo["HighlightRangesCriteriaUpperMap"]  = {}
    OptionsInfo["HighlightRangesCriteriaUpperValuesMap"]  = {}

    HighlightRanges = Options["--highlightValuesRanges"].strip()
    if re.match("^None$", HighlightRanges, re.I):
        return
    
    OptionsInfo["HighlightRanges"]  = HighlightRanges
    HighlightRangesWords = HighlightRanges.split(",")
    
    if len(HighlightRangesWords) % 6:
        MiscUtil.PrintError("The number of comma delimited paramater names and values, %d, specified in sextet \"%s\" using \"--highlightValuesRanges\" option must be a multiple of 6." % (len(HighlightRangesWords), HighlightRanges))

    for Index in range(0, len(HighlightRangesWords), 6):
        DataLabel = HighlightRangesWords[Index].strip()
        DataType = HighlightRangesWords[Index + 1].strip()
        LowerBoundDataCriterion = HighlightRangesWords[Index + 2].strip()
        LowerBoundDataValue = HighlightRangesWords[Index + 3].strip()
        UpperBoundDataCriterion = HighlightRangesWords[Index + 4].strip()
        UpperBoundDataValue = HighlightRangesWords[Index + 5].strip()

        SpecifiedSextet = "%s,%s,%s,%s,%s,%s" % (DataLabel, DataType, LowerBoundDataCriterion, LowerBoundDataValue, UpperBoundDataCriterion, UpperBoundDataValue)

        CanonicalDataLabel = DataLabel.lower()
        if CanonicalDataLabel in OptionsInfo["HighlightValuesCanonicalLabelsMap"]:
            MiscUtil.PrintError("The data label specified, %s, using option \"--highlightRanges\" has already been used in \"--highlightValues\" option" % DataLabel)
        
        if not re.match("^(numeric|text)$", DataType, re.I):
            MiscUtil.PrintError("The data type, %s, specified in sextet \"%s\" using \"--highlightValuesRanges\" option is not valid. Supported values: numeric text." % (DataType, SpecifiedSextet))
        
        if not re.match("^(lt|le)$", LowerBoundDataCriterion, re.I):
            MiscUtil.PrintError("The lower bound criterion, %s, specified in sextet \"%s\" using \"--highlightValuesRanges\" option is not valid. Supported values: lt or le." % (LowerBoundDataCriterion, SpecifiedSextet))
        
        if not re.match("^(gt|ge)$", UpperBoundDataCriterion, re.I):
            MiscUtil.PrintError("The upper bound criterion, %s, specified in sextet \"%s\" using \"--highlightValuesRanges\" option is not valid. Supported values: gt or ge." % (UpperBoundDataCriterion, SpecifiedSextet))

        if re.match("^numeric$", DataType, re.I):
            if not MiscUtil.IsNumber(LowerBoundDataValue):
                MiscUtil.PrintError("The lower bound data value, %s, specified in sextet \"%s\", using \"--highlightValuesRanges\" option is not valid. It must be a number for \"%s\" data type." % (LowerBoundDataValue, SpecifiedSextet, DataType))
            
            if not MiscUtil.IsNumber(UpperBoundDataValue):
                MiscUtil.PrintError("The upper bound data value, %s, specified in sextet \"%s\", using \"--highlightValuesRanges\" option is not valid. It must be a number for \"%s\"data type." % (UpperBoundDataValue, SpecifiedSextet, DataType))
            
            if float(LowerBoundDataValue) >= float(UpperBoundDataValue):
                MiscUtil.PrintError("The lower bound data value, %s, must be less than upper bound value, %s, specified in sextet \"%s\" using \"--highlightValuesRanges\" option." % (LowerBoundDataValue, UpperBoundDataValue,  SpecifiedSextet))
            
            LowerBoundDataValue = float(LowerBoundDataValue)
            UpperBoundDataValue = float(UpperBoundDataValue)
        else:
            if LowerBoundDataValue >= UpperBoundDataValue:
                MiscUtil.PrintError("The lower bound data value, %s, must be less than upper bound value, %s, specified in sextet \"%s\", using \"--highlightValuesRanges\" option is not valid. It must be a number for data type, %s" % (LowerBoundDataValue, UpperBoundDataValue,  SpecifiedSextet, DataType))

        # Track values...
        if DataLabel in OptionsInfo["HighlightRangesLabelsMap"]:
            MiscUtil.PrintError("The data label, %s, specified in sextet \"%s\", using \"--highlightValuesRanges\" option is not valid. Multiple occurences of data label" % (DataLabel, SpecifiedSextet))

        OptionsInfo["HighlightRangesLabels"].append(DataLabel)
        OptionsInfo["HighlightRangesLabelsMap"][DataLabel]  = DataLabel
        OptionsInfo["HighlightRangesCanonicalLabelsMap"][CanonicalDataLabel]  = DataLabel
        
        OptionsInfo["HighlightRangesTypesMap"][DataLabel]  = DataType
        
        OptionsInfo["HighlightRangesCriteriaLowerMap"][DataLabel]  = LowerBoundDataCriterion 
        OptionsInfo["HighlightRangesCriteriaLowerValuesMap"][DataLabel]  = LowerBoundDataValue
        OptionsInfo["HighlightRangesCriteriaUpperMap"][DataLabel]  = UpperBoundDataCriterion 
        OptionsInfo["HighlightRangesCriteriaUpperValuesMap"][DataLabel]  = UpperBoundDataValue

def ProcessHighlightValuesClassesOption():
    """Process highlight values classes option"""
    
    OptionsInfo["HighlightClasses"]  = None
    OptionsInfo["HighlightClassesRules"]  = None
    OptionsInfo["HighlightClassesSynonymsMap"]  = None
    OptionsInfo["HighlightClassesRandom"]  = False

    OptionsInfo["HighlightClassesLabels"]  = []
    OptionsInfo["HighlightClassesLabelsMap"]  = {}
    OptionsInfo["HighlightClassesCanonicalLabelsMap"]  = {}
    
    OptionsInfo["HighlightClassesTypesMap"]  = {}
    OptionsInfo["HighlightClassesCriteriaMap"]  = {}
    OptionsInfo["HighlightClassesCriteriaValuesMap"]  = {}
    
    HighlightClasses = Options["--highlightValuesClasses"].strip()
    if re.match("^None$", HighlightClasses, re.I):
        return

    OptionsInfo["HighlightClasses"]  = HighlightClasses
    
    if re.match("^RuleOf5$", HighlightClasses, re.I):
        HighlightClassessRules = "MolecularWeight,numeric,le,500,HydrogenBondDonors,numeric,le,5,HydrogenBondAcceptors,numeric,le,10,LogP,numeric,le,5"
    elif re.match("^RuleOf3$", HighlightClasses, re.I):
        HighlightClassessRules = "MolecularWeight,numeric,le,300,HydrogenBondDonors,numeric,le,3,HydrogenBondAcceptors,numeric,le,3,LogP,numeric,le,3,RotatableBonds,numeric,le,3,TPSA,numeric,le,60"
    elif re.match("^DrugLike$", HighlightClasses, re.I):
        HighlightClassessRules = "MolecularWeight,numeric,le,500,HydrogenBondDonors,numeric,le,5,HydrogenBondAcceptors,numeric,le,10,LogP,numeric,le,5,RotatableBonds,numeric,le,10,TPSA,numeric,le,140"
    elif re.match("^Random$", HighlightClasses, re.I):
        if OptionsInfo["HighlightValues"] is not None:
            MiscUtil.PrintError("The value specified, %s, using option \"--highlightValuesClasses\" is not allowed in conjunction with \"--highlightValues\" option." % HighlightClasses)
        if OptionsInfo["HighlightRanges"] is not None:
            MiscUtil.PrintError("The value specified, %s, using option \"--highlightValuesClasses\" is not allowed in conjunction with \"--highlightRanges\" option ." % HighlightClasses)
            
        OptionsInfo["HighlightClassesRandom"]  = True
        return
    else:
        MiscUtil.PrintError("The value specified, %d, using option \"--highlightValuesClasses\" is not supported." % HighlightClasses)
        return
        
    OptionsInfo["HighlightClassesRules"]  = HighlightClassessRules

    # Process rules for highlighting values...
    HighlightClassesWords = HighlightClassessRules.split(",")
    for Index in range(0, len(HighlightClassesWords), 4):
        DataLabel = HighlightClassesWords[Index].strip()
        DataType = HighlightClassesWords[Index + 1].strip()
        DataCriterion = HighlightClassesWords[Index + 2].strip()
        DataValue = HighlightClassesWords[Index + 3].strip()
        
        DataValue = float(DataValue)

        if DataLabel in OptionsInfo["HighlightClassesLabelsMap"]:
            MiscUtil.PrintWarning("Ignoring duplicate datalabel, %s, specified in highlighting values rule for class, %s, in \"--highlightClassesValue\" option..." % (DataLabel, HighlightClasses))
            continue
            
        OptionsInfo["HighlightClassesLabels"].append(DataLabel)
        OptionsInfo["HighlightClassesLabelsMap"][DataLabel] = DataLabel
        
        OptionsInfo["HighlightClassesTypesMap"][DataLabel] = DataType
        OptionsInfo["HighlightClassesCriteriaMap"][DataLabel] = DataCriterion
        OptionsInfo["HighlightClassesCriteriaValuesMap"][DataLabel] = DataValue
    
    # Set up synonyms for data labels corresponding to physicochemical properties
    # calculated by MayaChemTools and RDKit...
    OptionsInfo["HighlightClassesSynonymsMap"]  = {}
    OptionsInfo["HighlightClassesSynonymsMap"]["MolecularWeight"]  = ["MolecularWeight", "MolWt"]
    OptionsInfo["HighlightClassesSynonymsMap"]["HydrogenBondDonors"]  = ["HydrogenBondDonors", "NHOHCount"]
    OptionsInfo["HighlightClassesSynonymsMap"]["HydrogenBondAcceptors"]  = ["HydrogenBondAcceptors", "NOCount"]
    OptionsInfo["HighlightClassesSynonymsMap"]["LogP"]  = ["SLogP", "MolLogP"]
    OptionsInfo["HighlightClassesSynonymsMap"]["RotatableBonds"]  = ["RotatableBonds", "NumRotatableBonds"]
    OptionsInfo["HighlightClassesSynonymsMap"]["TPSA"]  = ["TPSA", "TPSA"]

def ProcessHighlightColorsOption():
    """Process highlight colors option"""

    OptionsInfo["HighlightColors"] = None
    OptionsInfo["HighlightColorsType"] = None
    OptionsInfo["HighlightColorsList"] = None

    HighlightColors = "colorclass,table-success, table-danger"
    if not re.match("^auto$", Options["--highlightColors"], re.I):
        HighlightColors = Options["--highlightColors"].strip()
        if MiscUtil.IsEmpty(HighlightColors):
            MiscUtil.PrintError("The value specified using \"--highlightColors\" is empty.")

    OptionsInfo["HighlightColors"] = re.sub(" ", "", HighlightColors)
    HighlightColorsList = [Color.lower() for Color in OptionsInfo["HighlightColors"].split(",")]
        
    if len(HighlightColorsList) != 3:
        MiscUtil.PrintError("The number of comma delimited paramater names and values, %d, specified using \"--highlightColors\" option must be 3." % (len(HighlightColorsList)))
        
    ColorsType, Color1, Color2 = HighlightColorsList
    if not re.match("^(colorclass|colorspec)$", ColorsType, re.I):
        MiscUtil.PrintError("The color type, %s, specified using \"--highlightColors\" option is not valid. Supported values: colorclass or colorspec." % ColorsType)

    ColorsList = [Color1, Color2]
    if re.match("^colorclass$", ColorsType, re.I):
        CheckOptionTableClassColorValues("--highlightColors", ColorsList)
        
    OptionsInfo["HighlightColorsList"] = ColorsList
    OptionsInfo["HighlightColorsType"] = ColorsType

def ProcessHighlightColorsRangesOption():
    """Process highlight colors ranges option"""

    OptionsInfo["HighlightColorsRanges"] = None
    OptionsInfo["HighlightColorsRangesType"] = None
    OptionsInfo["HighlightColorsRangesList"] = None

    HighlightColors = "colorclass,table-success, table-warning, table-danger"
    if not re.match("^auto$", Options["--highlightColorsRanges"], re.I):
        HighlightColors = Options["--highlightColorsRanges"].strip()
        if MiscUtil.IsEmpty(HighlightColors):
            MiscUtil.PrintError("The value specified using \"--highlightColorsRanges\" is empty.")

    OptionsInfo["HighlightColorsRanges"] = re.sub(" ", "", HighlightColors)
    HighlightColorsList = [Color.lower() for Color in OptionsInfo["HighlightColorsRanges"].split(",")]
        
    if len(HighlightColorsList) != 4:
        MiscUtil.PrintError("The number of comma delimited paramater names and values, %d, specified using \"--highlightColorsRanges\" option must be 4." % (len(HighlightColorsList)))
    
    ColorsType, Color1, Color2, Color3 = HighlightColorsList
    if not re.match("^(colorclass|colorspec)$", ColorsType, re.I):
        MiscUtil.PrintError("The color type, %s, specified using \"--highlightColorsRanges\" option is not valid. Supported values: colorclass or colorspec." % ColorsType)
    
    ColorsList = [Color1, Color2, Color3]
    if re.match("^colorclass$", ColorsType, re.I):
        CheckOptionTableClassColorValues("--highlightColorsRanges", ColorsList)
    
    OptionsInfo["HighlightColorsRangesList"] = ColorsList
    OptionsInfo["HighlightColorsRangesType"] = ColorsType

def ProcessHighlightColorsRandomOption():
    """Process highlight colors random option"""

    OptionsInfo["HighlightColorsRandom"] = None
    OptionsInfo["HighlightColorsRandomType"] = None
    OptionsInfo["HighlightColorsRandomList"] = None

    HighlightColors = "colorclass,table-primary,table-success,table-danger,table-info,table-warning,table-secondary"
    if not re.match("^auto$", Options["--highlightColorsRandom"], re.I):
        HighlightColors = Options["--highlightColorsRandom"].strip()
        if MiscUtil.IsEmpty(HighlightColors):
            MiscUtil.PrintError("The value specified using \"--highlightColorsRandom\" is empty.")

    OptionsInfo["HighlightColorsRandom"] = re.sub(" ", "", HighlightColors)
    HighlightColorsList = [Color.lower() for Color in OptionsInfo["HighlightColorsRandom"].split(",")]
        
    if len(HighlightColorsList) <= 1:
        MiscUtil.PrintError("The number of comma delimited paramater names and values, %d, specified using \"--highlightColorsRandom\" option must be > 1." % (len(HighlightColorsList)))
        
    ColorsType = HighlightColorsList[0]
    ColorsList = HighlightColorsList[1:]
    
    if not re.match("^(colorclass|colorspec)$", ColorsType, re.I):
        MiscUtil.PrintError("The color type, %s, specified using \"--highlightColorsRandim\" option is not valid. Supported values: colorclass or colorspec." % ColorsType)

    if re.match("^colorclass$", ColorsType, re.I):
        CheckOptionTableClassColorValues("--highlightColorsRandom", ColorsList)
        
    OptionsInfo["HighlightColorsRandomList"] = ColorsList
    OptionsInfo["HighlightColorsRandomType"] = ColorsType

def CheckOptionTableClassColorValues(OptionName, ColorsList):
    """Check names of table color classes and issue a warning for unknown names."""

    TableClassColors = ["thead-dark", "thead-light", "table-primary", "table-success", "table-danger", "table-info", "table-warning", "table-active", "table-secondary", "table-light", "table-dark", "bg-primary", "bg-success", "bg-danger",  "bg-info", "bg-warning", "bg-secondary", "bg-dark", "bg-light"]

    for Color in ColorsList:
        if not Color in TableClassColors:
            MiscUtil.PrintWarning("The color class name, %s, specified using option \"%s\" appears to be a unknown name..." % (Color, OptionName))
        
def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    
    # No need for any RDKit specific --outfileParams....
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], OptionsInfo["Infile"])

    AlignmentSMARTSPattern = None
    if not re.match("^None$", Options["--alignmentSMARTS"], re.I):
        AlignmentSMARTSPattern = Options["--alignmentSMARTS"]
    OptionsInfo["AlignmentSMARTSPattern"]  = AlignmentSMARTSPattern
    
    Compute2DCoords = True
    if re.match("^no$", Options["--compute2DCoords"], re.I):
        Compute2DCoords = False
    OptionsInfo["Compute2DCoords"]  = Compute2DCoords

    CounterCol = True
    if re.match("^no$", Options["--counterCol"], re.I):
        CounterCol = False
    OptionsInfo["CounterCol"]  = CounterCol
    
    ColVisibility = True
    if re.match("^no$", Options["--colVisibility"], re.I):
        ColVisibility = False
    OptionsInfo["ColVisibility"]  = ColVisibility
    
    OptionsInfo["ColVisibilityCtrlMax"]  = int(Options["--colVisibilityCtrlMax"])
    
    Footer = None
    if not re.match("^None$", Options["--footer"], re.I):
        Footer = Options["--footer"]
    OptionsInfo["Footer"]  = Footer

    FooterClass = Options["--footerClass"].strip()
    if MiscUtil.IsEmpty(FooterClass):
        MiscUtil.PrintError("The value specified using option \"--footerClass\" is empty.")
    OptionsInfo["FooterClass"]  = FooterClass
    
    FreezeCols = True
    if re.match("^no$", Options["--freezeCols"], re.I):
        FreezeCols = False
    OptionsInfo["FreezeCols"]  = FreezeCols
    
    Header = None
    if not re.match("^None$", Options["--header"], re.I):
        Header = Options["--header"]
    OptionsInfo["Header"]  = Header
    
    HeaderStyle = Options["--headerStyle"].strip()
    if MiscUtil.IsEmpty(HeaderStyle):
        MiscUtil.PrintError("The value specified using option \"--headerStyle\" is empty.")
    OptionsInfo["HeaderStyle"]  = HeaderStyle

    ProcessHighlightSMARTSOption()
    ProcessHighlightDataOptions()
    
    OptionsInfo["KeysNavigation"] = True
    if re.match("^no$", Options["--keysNavigation"], re.I):
        OptionsInfo["KeysNavigation"] = False
    
    SizeValues = Options["--molImageSize"].split(",")
    OptionsInfo["MolImageWidth"] = int(SizeValues[0])
    OptionsInfo["MolImageHeight"] = int(SizeValues[1])
    
    OptionsInfo["Paging"] = True
    if re.match("^no$", Options["--paging"], re.I):
        OptionsInfo["Paging"] = False
    
    PagingType = Options["--pagingType"]
    if not re.match("^(numbers|simple|simple_numbers|full|full_numbers|simple_number)$", Options["--pagingType"], re.I):
        MiscUtil.PrintWarning("The paging type name, %s, specified using option \"--pagingType\" appears to be a unknown type..." % (PagingType))
    OptionsInfo["PagingType"] = PagingType.lower()
    
    OptionsInfo["PageLength"] = int(Options["--pageLength"])
    
    OptionsInfo["RegexSearch"] = True
    if re.match("^no$", Options["--regexSearch"], re.I):
        OptionsInfo["RegexSearch"] = False
    
    OptionsInfo["ShowMolName"] = True
    OptionsInfo["ShowMolNameDataLabel"] = "Name"
    if re.match("^no$", Options["--showMolName"], re.I):
        OptionsInfo["ShowMolName"] = False
    
    OptionsInfo["ShowMolNameAuto"] = True if re.match("^auto$", Options["--showMolName"], re.I) else False

    OptionsInfo["ScrollX"] = True
    if re.match("^no$", Options["--scrollX"], re.I):
        OptionsInfo["ScrollX"] = False
        
    OptionsInfo["ScrollY"] = True
    if re.match("^no$", Options["--scrollY"], re.I):
        OptionsInfo["ScrollY"] = False

    OptionsInfo["ScrollYSize"] = Options["--scrollYSize"]
    if re.match("vh$", Options["--scrollYSize"], re.I):
        ScrollYSize = int(re.sub("vh$", "", Options["--scrollYSize"]))
        if ScrollYSize <= 0:
            MiscUtil.PrintError("The value specified, %s, for option \"--scrollYSize\" is not valid. Supported value: > 0 followed by \"vh\"" % Options["--scrollYSize"])
    
    TableStyle = None
    if not re.match("^None$", Options["--tableStyle"], re.I):
        if re.match("^All$", Options["--tableStyle"], re.I):
            TableStyle = "table table-striped table-bordered table-hover table-dark"
        else:
            TableStyle = re.sub(" ", "", Options["--tableStyle"])
            for Style in [Style for Style in TableStyle.split(",")]:
                if not re.match("^(table|table-striped|table-bordered|table-hover|table-dark|table-sm)$", Style, re.I):
                    MiscUtil.PrintWarning("The table style name, %s, specified using option \"-t, --tableStyle\" appears to be a unknown style..." % (Style))
            TableStyle = re.sub(",", " ", TableStyle.lower())
    OptionsInfo["TableStyle"]  = TableStyle

    TableHeaderStyle = None
    if not re.match("^None$", Options["--tableHeaderStyle"], re.I):
        TableHeaderStyle = Options["--tableHeaderStyle"]
        TableHeaderStyle = TableHeaderStyle.lower()
        CheckOptionTableClassColorValues("--tableHeaderStyle", [TableHeaderStyle])
    OptionsInfo["TableHeaderStyle"]  = TableHeaderStyle
    
    OptionsInfo["TableFooter"] = True
    if re.match("^no$", Options["--tableFooter"], re.I):
        OptionsInfo["TableFooter"] = False

    OptionsInfo["WrapText"] = True
    if re.match("^no$", Options["--wrapText"], re.I):
        OptionsInfo["WrapText"] = False

    OptionsInfo["WrapTextWidth"] = int(Options["--wrapTextWidth"])

def RetrieveOptions():
    """Retrieve command line arguments and options"""
    
    # Get options...
    global Options
    Options = docopt(_docoptUsage_)
    
    # Set current working directory to the specified directory...
    WorkingDir = Options["--workingdir"]
    if WorkingDir:
        os.chdir(WorkingDir)
    
    # Handle examples option...
    if "--examples" in Options and Options["--examples"]:
        MiscUtil.PrintInfo(MiscUtil.GetExamplesTextFromDocOptText(_docoptUsage_))
        sys.exit(0)

def ValidateOptions():
    """Validate option values"""
    
    if not re.match("^None$", Options["--alignmentSMARTS"], re.I):
        PatternMol = Chem.MolFromSmarts(Options["--alignmentSMARTS"])
        if PatternMol is None:
            MiscUtil.PrintError("The value specified, %s, using option \"--alignmentSMARTS\" is not a valid SMARTS: Failed to create pattern molecule" % Options["--alignmentSMARTS"])
    
    MiscUtil.ValidateOptionTextValue("-c, --compute2DCoords", Options["--compute2DCoords"], "yes no auto")
    
    MiscUtil.ValidateOptionTextValue("--counterCol", Options["--counterCol"], "yes no")
    MiscUtil.ValidateOptionTextValue("--colVisibility", Options["--colVisibility"], "yes no")
    MiscUtil.ValidateOptionIntegerValue("--colVisibilityCtrlMax", Options["--colVisibilityCtrlMax"], {">": 0})
    
    MiscUtil.ValidateOptionTextValue("--freezeCols", Options["--freezeCols"], "yes no")
    MiscUtil.ValidateOptionTextValue("--highlightValuesClasses", Options["--highlightValuesClasses"], "RuleOf5 RuleOf3 DrugLike Random None")
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi csv tsv txt")

    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "html")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    
    MiscUtil.ValidateOptionTextValue("-k, --keysNavigation", Options["--keysNavigation"], "yes no")
    
    MiscUtil.ValidateOptionNumberValues("-m, --molImageSize", Options["--molImageSize"], 2, ",", "integer", {">": 0})
    
    MiscUtil.ValidateOptionTextValue("-p, --paging", Options["--paging"], "yes no")
    MiscUtil.ValidateOptionIntegerValue("--pageLength", Options["--pageLength"], {">": 0})
    MiscUtil.ValidateOptionTextValue("-r, --regexSearch", Options["--regexSearch"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--showMolName", Options["--showMolName"], "yes no auto")
    
    MiscUtil.ValidateOptionTextValue("--scrollX", Options["--scrollX"], "yes no")
    MiscUtil.ValidateOptionTextValue("--scrollY", Options["--scrollY"], "yes no")
    if not re.search("vh$", Options["--scrollYSize"], re.I):
        MiscUtil.ValidateOptionIntegerValue("--scrollYSize", Options["--scrollYSize"], {">": 0})
    
    MiscUtil.ValidateOptionTextValue("--tableFooter", Options["--tableFooter"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--wrapText", Options["--wrapText"], "yes no")
    MiscUtil.ValidateOptionIntegerValue("--wrapTextWidth", Options["--wrapTextWidth"], {">": 0})

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitDrawMoleculesAndDataTable.py - Generate a HTML data table

Usage:
    RDKitDrawMoleculesAndDataTable.py [--alignmentSMARTS <SMARTS>]
                             [--compute2DCoords <yes or  no>] [--counterCol <yes or no>]
                             [--colVisibility <yes or no>] [--colVisibilityCtrlMax <number>] [--footer <text>]
                             [--footerClass <text>] [--freezeCols <yes or no>] [--header <text>]
                             [--headerStyle <text>] [--highlightSMARTS <SMARTS,...>]
                             [--highlightValues <datalabel,datatype,criterion,value,...>]
                             [--highlightValuesRanges <datalabel,datatype,criterion1,vaue1,criterion2,value2...>]
                             [--highlightValuesClasses <RuleOf5,RuleOf3,...>]
                             [--highlightColors <colortype,color1,color2>]
                             [--highlightColorsRanges <colortype,color1,color2,color3>]
                             [--highlightColorsRandom <colottype,color1,color2,...>]
                             [--infileParams <Name,Value,...>] [--keysNavigation <yes or no>]
                             [--molImageSize <width,height>] [--overwrite] [--paging <yes or no>]
                             [--pagingType <numbers,simple, ...>] [--pageLength <number>]
                             [--regexSearch <yes or no>] [--showMolName <yes or no>]
                             [--scrollX <yes or no>] [--scrollY <yes or no>] [--scrollYSize <number>]
                             [--tableStyle <table,table-striped,...>] [--tableFooter <yes or no>]
                             [--tableHeaderStyle <thead-dark,thead-light,...>] [--wrapText <yes or no>] 
                             [--wrapTextWidth <number>] [-w <dir>] -i <infile> -o <outfile>
    RDKitDrawMoleculesAndDataTable.py -h | --help | -e | --examples

Description:
    Generate an interactive HTML table with columns corresponding to molecules
    and available alphanumerical data in an input file. The drawing of molecules are
    embedded in the columns as in line SVG images.

    The interactive HTML table may contain multiple columns with drawing of
    molecules. These columns are automatically generated for each data field in SD
    file or a column name in SMILES and CSV/TSV file containing SMILES
    string in their names. The first molecular drawing column in the HTML table
    represents primary molecular structure data available in an input file. It
    corresponds to MOL block is SD file or a first column containing SMILES string
    in its name in SMILES and CSV/TSV files.
 
    The interactive table requires internet access for viewing in a browser and
    employs the following frameworks: JQuery, Bootstrap, and DataTable. It provides
    the following functionality: sorting by columns, page length control, page 
    navigation, searching data with regular expressions, and horizontal/vertical
    scrolling, row highlighting during hovering, a counter column, freezing of primary
    structure and counter columns, and column visibility control.

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd), SMILES (.smi),
    CSV/TSV (.csv, .tsv, .txt)

    The supported output file format is HTML (.html).

Options:
    -a, --alignmentSMARTS <SMARTS>  [default: none]
        SMARTS pattern for aligning molecules to a common template. This option is
        only used for primary molecular data in SD, SMILES and CSV/TSV files. It is 
        ignored for all other molecular coordinates corresponding to data fields in SD
        file or columns in SMILES and CSV/TSV files containing SMILES string in their
        names.
    -c, --compute2DCoords <yes or no>  [default: yes]
        Compute 2D coordinates of molecules before drawing. Default: yes for SMILES
        strings in SMILES, CSV/TSV, and SD file data fields. In addition, 2D coordinated are
        always calculated for molecules corresponding to data fields in SD file or columns
        in SMILES and CSV/TSV files containing SMILES string in their names.
    --counterCol <yes or no>  [default: yes]
        Show a counter column as the first column in the table. It contains the position
        for each row in the table.
    --colVisibility <yes or no>  [default: yes]
        Show a dropdown button to toggle visibility of columns in the table. The counter
        and primary structure columns are excluded from the list.
    --colVisibilityCtrlMax <number>  [default: 25]
        Maximum number of columns to show in column visibility dropdown button. The
        rest of the data columns are not listed in the dropdown and are shown in the table.
        A word to the wise: The display of too many columns appear to hang interactive
        Javascript framework for Bootstrap and DataTables.
    --freezeCols <yes or no>  [default: yes]
        Lock counter and primary structure columns in place during horizontal scrolling.
    --footer <text>  [default: none]
        Footer text to insert at the bottom of the HTML page after the table.
    --footerClass <text>  [default: small text-center text-muted]
        Footer class style to use with <p> tag.
    -h, --help
        Print this help message.
    --header <text>  [default: none]
        Header text to insert at the top of the HTML page before the table.
    --headerStyle <text>  [default: h5]
        Header style to use. Possible values: h1 to h6.
    --highlightSMARTS <SMARTS,...>  [default: none]
        SMARTS pattern for highlighting atoms and bonds in molecules. All matched
        substructures are highlighted.
        
        The SMARTS string is used to highlight atoms and bonds in drawing of
        molecules present in a HTML table across multiple columns. These columns
        correspond to data field labels in SD file or a column name in SMILES and
        CSV/TSV file containing SMILES string in their names. The first molecular
        drawing column in HTML table corresponds to primary molecular structure
        data available in an input file. It is identified by a label 'Structure' across
        all input formats.
        
        A single SMARTS string is used to highlight a common substructure across
        all columns containing drawing of molecules in HTML table.
        
        Format:
            
            SMARTS
            Structure,SMARTS1,DataLabel,SMARTS2,...
            Structure,SMARTS1,Collabel,SMARTS2,...
            
        Example:
            
            c1ccccc1
            Structure,c1ccccc1,SMILESR1,c1ccccc1,SMILESR2,c1ccccc1
            
    --highlightValues <datalabel,datatype,criterion,value,...>  [default: none]
        Highlighting methodology to use for highlighting  alphanumerical data
        corresponding to data fields in SD file or column names in SMILES and
        CSV/TSV text files.
        
        Input text contains these quartets: DataLabel, DataType, Criterion, Value.
        Possible datatype values: numeric, text. Possible criterion values for numeric
        and text: gt, lt, ge, le.
        
        The 'datalabel' corresponds to either data field label in SD file or column name
        in SMILES and CSV/TSV text files.
        
        Examples:
            
            MolecularWeight,numeric,le,500
            MolecularWeight,numeric,le,450,SLogP,numeric,le,5
            Name,text,eq,Aspirin
            Name,regex,eq,acid|amine
            
    --highlightValuesRanges <datalabel,datatype,...>  [default: none]
        Highlighting methodology to use for highlighting ranges of alphanumerical
        data corresponding to data fields in SD file or column names in SMILES and
        CSV/TSV text files.
        
        Input text contains these sextets: DataLabel, DataType, CriterionLowerBound,
        LowerBoundValue, CriterionUpperBound, UpperBoundValue.
        
        Possible datatype values: numeric or text. Possible criterion values: Lower
        bound value - lt, le; Upper bound value: gt, ge.
        
        The 'datalabel' corresponds to either data field label in SD file or column name
        in SMILES and CSV/TSV text files.
        
        Examples:
            
            MolecularWeight,numeric,lt,450,gt,1000
            MolecularWeight,numeric,lt,450,gt,1000,SLogP,numeric,lt,0,gt,5
            
    --highlightValuesClasses <RuleOf5,RuleOf3,...>  [default: none]
        Highlighting methodology to use for highlighting ranges of numerical data
        data corresponding to specific set of data fields in SD file or column names in
        SMILES and CSV/TSV text files. Possible values: RuleOf5, RuleOf3, DrugLike,
        Random.
        
        The following value classes are supported: RuleOf5, RuleOf3, LeadLike, DrugLike.
        LeadLike is equivalent to RuleOf3.
        
        Each supported class encompasses a specific set of data labels along with
        appropriate criteria to compare and highlight column values, except for
        'Random' class. The data labels in these classes are automatically associated
        with appropriate data fields in SD file or column names in SMILES and CSV/TSV
        text files.
        
        No data labels are associated with 'Random' class. It is used to highlight
        available alphanumeric data by randomly selecting a highlight color from the
        list of colors specified using '--highlightColorsRandom' option. The 'Random'
        class value is not allowed in conjunction with '--highlightValues' or
        '--highlightValuesRanges'.
        
        The rules to highlight values for the supported classes are as follows.
        
        RuleOf5 [ Ref 91 ]:
         
            MolecularWeight,numeric,le,500 (MolecularWeight <= 500)
            HydrogenBondDonors,numeric,le,5 (HydrogenBondDonors <= 5)
            HydrogenBondAcceptors,numeric,le,10 (HydrogenBondAcceptors <= 10)
            LogP,numeric,le,5 (LogP <= 5)
         
        RuleOf3 or LeadLike [ Ref 92 ]:
         
            MolecularWeight,numeric,le,300 (MolecularWeight <= 300)
            HydrogenBondDonors,numeric,le,3 (HydrogenBondDonors <= 3)
            HydrogenBondAcceptors,numeric,le,3 (HydrogenBondAcceptors <= 3)
            LogP,numeric,le,3 (LogP <= 3)
            RotatableBonds,numeric,le,3 (RotatableBonds <= 3)
            TPSA,numeric,le,60 (TPSA <= 60)
         
        DrugLike:
         
            MolecularWeight,numeric,le,500 (MolecularWeight <= 500)
            HydrogenBondDonors,numeric,le,5 (HydrogenBondDonors <= 5)
            HydrogenBondAcceptors,numeric,le,10 (HydrogenBondAcceptors <= 10)
            LogP,numeric,le,5 (LogP <= 5)
            RotatableBonds,numeric,le,10 (RotatableBonds <= 10)
            TPSA,numeric,le,140 (TPSA <= 140)
            
        The following synonyms are automatically detected for data labels used
        by MayaChemTools and RDKit packages during the calculation of
        physicochemical properties.
        
        MayaChemTools: MolecularWeight, HydrogenBondDonors, HydrogenBondAcceptors,
        SLogP, RotatableBonds, TPSA.
            
        RDKit: MolWt,  NHOHCount, NOCount, MolLogP, NumRotatableBonds, TPSA
        
    --highlightColors <colortype,color1,color2>  [default: auto]
        Background colors used to highlight column values based on criterion
        specified by '--highlightValues' and '--highlightColorsClasses' option. Default
        value: colorclass,table-success, table-danger.
        
        The first color is used to highlight column values that satisfy the specified
        criterion for the column. The second color highlights the rest of the values
        in the column. 
        
        Possible values for colortype: colorclass or colorspec.
        
        Any valid bootstrap contextual color class is supported for 'colorclass'
        color type. For example: table-primary (Blue), table-success (Green),
        table-danger (Red), table-info (Light blue), table-warning (Orange),
        table-secondary (Grey), table-light (Light grey), and  table-dark (Dark grey).
        
        The following bootstrap color classes may also used: bg-primary bg-success,
        bg-danger bg-info, bg-warning, bg-secondary.
        
        Any valid color name or hexadecimal color specification is supported for
        'colorspec' color type: For example: red, green, blue, #ff000, #00ff00, #0000ff.
    --highlightColorsRanges <colortype,color1,color2,color3>  [default: auto]
        Background colors used to highlight column values using criteria specified
        by '--highlightValuesRanges' option. Default value:  colorclass, table-success,
        table-warning, table-danger.
        
        The first and third color are used to highlight column values lower and higher
        than the specified values for the lower and upper bound. The middle color highlights
        the rest of the values in the column.
        
        The supported color type and values are explained in the section for '--highlightColors'.
    --highlightColorsRandom <colortype,color1,color2,...>  [default: auto]
        Background color list to use for randomly selecting a color  to highlight
        column values during 'Random" value of '--highlightValuesClasses' option.
        
        Default value:  colorclass,table-primary,table-success,table-danger,table-info,
        table-warning,table-secondary.
        
        The supported color type and values are explained in the section for '--highlightColors'.
    -i, --infile <infile>
        Input file name.
    --infileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for reading
        molecules from files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD, MOL: removeHydrogens,yes,sanitize,yes,strictParsing,yes
            SMILES: smilesColumn,1,smilesNameColumn,2,smilesDelimiter,space,
                sanitize,yes
            
        Possible values for smilesDelimiter: space, comma or tab.
    -k, --keysNavigation <yes or no>  [default: yes]
        Provide Excel like keyboard cell navigation for the table.
    -m, --molImageSize <width,height>  [default: 200,150]
        Image size of a molecule in pixels.
    -o, --outfile <outfile>
        Output file name.
    --overwrite
        Overwrite existing files.
    -p, --paging <yes or no>  [default: yes]
        Provide page navigation for browsing data in the table.
    --pagingType <numbers, simple, ...>  [default: full_numbers]
        Type of page navigation. Possible values: numbers, simple, simple_numbers,
        full, full_numbers, or first_last_numbers.
            
            numbers - Page number buttons only
            simple - 'Previous' and 'Next' buttons only
            simple_numbers - 'Previous' and 'Next' buttons, plus page numbers
            full - 'First', 'Previous', 'Next' and 'Last' buttons
            full_numbers - 'First', 'Previous', 'Next' and 'Last' buttons, plus
                page numbers
            first_last_numbers - 'First' and 'Last' buttons, plus page numbers
            
    --pageLength <number>  [default: 15]
        Number of rows to show per page.
    -r, --regexSearch <yes or no>  [default: yes]
        Allow regular expression search through alphanumerical data in the table.
    -s, --showMolName <yes or no>  [default: auto]
        Show molecule names in a column next to the column corresponding to primary
        structure data in SD and SMILES file. The default value is yes for SD and SMILES file.
        This option is ignored for CSV/TSV text files.
    --scrollX <yes or no>  [default: yes]
        Provide horizontal scroll bar in the table as needed.
    --scrollY <yes or no>  [default: yes]
        Provide vertical scroll bar in the table as needed.
    --scrollYSize <number>  [default: 75vh]
        Maximum height of table viewport either in pixels or percentage of the browser
        window height before providing a vertical scroll bar. Default: 75% of the height of
        browser window.
    -t, --tableStyle <table,table-striped,...>  [default: table,table-hover,table-sm]
        Style of table. Possible values: table, table-striped, table-bordered,
        table-hover, table-dark, table-sm, none, or All. Default: 'table,table-hover'. A
        comma delimited list of any valid Bootstrap table styles is also supported.
    --tableFooter <yes or no>  [default: yes]
        Show column headers at the end of the table.
    --tableHeaderStyle <thead-dark,thead-light,...>  [default: thead-dark]
        Style of table header. Possible values: thead-dark, thead-light, or none.
        The names of the following contextual color classes are also supported:
        table-primary (Blue), table-success (Green), table-danger (Red), table-info
        (Light blue), table-warning (Orange), table-active (Grey), table-light (Light
        grey), and  table-dark (Dark grey).
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.
    --wrapText <yes or no>  [default: yes]
        Wrap alphanumeric text using <br/> delimiter for display in a HTML table.
    --wrapTextWidth <number>  [default: 40]
        Maximum width in characters before wraping alphanumeric text for display
        in a HTML table.

Examples:
    To generate a HTML table containing structure and alphanumeric data for
    molecules in a SD file along with all the bells and whistles to interact with
    the table, type:

        % RDKitDrawMoleculesAndDataTable.py -i Sample.sdf -o SampleOut.html

    To generate a HTML table containing structure and alphanumeric data for
    molecules in a SMILES file along with all the bells and whistles to interact
    with the table, type:

        % RDKitDrawMoleculesAndDataTable.py  -i Sample.smi -o SampleOut.html

    To generate a HTML table containing multiple structure columns for molecules
    in a CSV file along with all the bells and whistles to interact with the table, type:

        % RDKitDrawMoleculesAndDataTable.py -i SampleSeriesRGroupsD3R.csv
          -o SampleSeriesRGroupsD3ROut.html

    To generate a HTML table containing structure and alphanumeric data for
    molecules in a SD file along without any bells and whistles to interact with
    the table, type:

        % RDKitDrawMoleculesAndDataTable.py --colVisibility no --freezeCols no
          --keysNavigation no --paging no --regexSearch no --scrollX no
          --scrollY no -i Sample.sdf -o SampleOut.html

    To generate a HTML table containing structure and alphanumeric data for
    molecules in a SD file along with highlighting molecular weight values
    using a specified criterion, type:

        % RDKitDrawMoleculesAndDataTable.py  --highlightValues
          "MolecularWeight,numeric,le,500" -i Sample.sdf -o SampleOut.html

    To generate a HTML table containing structure and alphanumeric data for
    molecules in a SD file along with highlighting range of molecular weight values
    using a specified criterion, type:

        % RDKitDrawMoleculesAndDataTable.py  --highlightValuesRanges
          "MolecularWeight,numeric,lt,400,gt,500" -i Sample.sdf -o SampleOut.html

    To generate a HTML table containing structure and alphanumeric data for
    molecules in a SD file along with highlighting molecular weight values and
    ranges of SLogP values using a specified criterion and color schemes, type:

        % RDKitDrawMoleculesAndDataTable.py  --highlightValues
          "MolecularWeight,numeric,le,500" --highlightValuesRanges
          "SLogP,numeric,lt,0,gt,5" --highlightColors "colorclass,table-success,
          table-danger" --highlightColorsRanges "colorclass,table-danger,
          table-success,table-warning" -i Sample.sdf -o SampleOut.html

    To generate a HTML table containing structure and alphanumeric data for
    molecules in a SD file along with highlighting RuleOf5 physicochemical
    properties using a pre-defined set of criteria, type:

        % RDKitDrawMoleculesAndDataTable.py  --highlightValuesClasses RuleOf5
          -i Sample.sdf -o SampleOut.html

    To generate a HTML table containing structure and alphanumeric data for
    molecules in a SD file along with all the bells and whistles to interact
    with the table and highlight a specific SMARTS pattern in molecules, type:

        % RDKitDrawMoleculesAndDataTable.py  --highlightSMARTS "c1ccccc1"
          -i Sample.sdf -o SampleOut.html

    To generate a HTML table containing structure and alphanumeric data for
    molecules in a SD file along with highlighting of values using random colors
    from a default list of colors, type:

        % RDKitDrawMoleculesAndDataTable.py --highlightValuesClasses Random
          -i Sample.sdf -o SampleOut.html

    To generate a HTML table containing structure and alphanumeric data for
    molecules in a SD file along with highlighting of values using random colors
    from a specified list of colors, type:

        % RDKitDrawMoleculesAndDataTable.py --highlightValuesClasses Random
          --highlightColorsRandom "colorspec,Lavendar,MediumPurple,SkyBlue,
          CornflowerBlue,LightGreen,MediumSeaGreen,Orange,Coral,Khaki,Gold,
          Salmon,LightPink,Aquamarine,MediumTurquoise,LightGray" 
          -i Sample.sdf -o SampleOut.html

    To generate a HTML table containing structure and alphanumeric data for
    molecules in a SMILES file specific columns, type:

        % RDKitDrawMoleculesAndDataTable.py --infileParams "smilesDelimiter,
          comma, smilesColumn,1,smilesNameColumn,2"
          -i SampleSMILES.csv -o SampleOut.html

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitDrawMolecules.py, RDKitRemoveDuplicateMolecules.py,
    RDKitSearchFunctionalGroups.py, RDKitSearchSMARTS.py

Copyright:
    Copyright (C) 2018 Manish Sud. All rights reserved.

    The functionality available in this script is implemented using RDKit, an
    open source toolkit for cheminformatics developed by Greg Landrum.

    This file is part of MayaChemTools.

    MayaChemTools is free software; you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your option) any
    later version.

"""

if __name__ == "__main__":
    main()
