#!/bin/env python
#
# File: RDKitDrawMolecules.py
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

# RDKit imports...
try:
    from rdkit import rdBase
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw.MolDrawing import DrawingOptions
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
    DrawMolecules()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def DrawMolecules():
    """Draw molecules"""
    
    Infile = OptionsInfo["Infile"]
    Outfile = OptionsInfo["Outfile"]
    
    # Read molecules...
    MiscUtil.PrintInfo("\nReading file %s..." % Infile)

    ValidMols, MolCount, ValidMolCount  = RDKitUtil.ReadAndValidateMolecules(Infile, **OptionsInfo["InfileParams"])
    
    MiscUtil.PrintInfo("Total number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))

    # Compute 2D coordinates...
    if OptionsInfo["Compute2DCoords"]:
        MiscUtil.PrintInfo("\nComputing 2D coordinates...")
        for Mol in ValidMols:
            AllChem.Compute2DCoords(Mol)
    
    MiscUtil.PrintInfo("Generating image grid...")
    
    # Setup atoms lists for highlighting atoms and bonds...
    AtomLists = SetupAtomListsToHighlight(ValidMols)
    BondLists = None
        
    # Set up legends...
    MolNames = None
    if OptionsInfo["ShowMolName"]:
        MolNames = []
        MolCount = 0
        for Mol in ValidMols:
            MolCount += 1
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            MolNames.append(MolName)
    
    # Perform alignment to a common template...
    PerformAlignment(ValidMols)

    # Generate appropriate output files...
    if MiscUtil.CheckFileExt(Outfile, "svg"):
        GenerateSVGImageFile(ValidMols, MolNames, AtomLists, BondLists)
    elif MiscUtil.CheckFileExt(Outfile, "html htm"):
        GenerateHTMLTableFile(ValidMols, MolNames, AtomLists, BondLists)
    else:
        GenerateImageFile(ValidMols, MolNames, AtomLists, BondLists)
    
def GenerateSVGImageFile(ValidMols, MolNames, AtomLists, BondLists):
    """Generate a SVG image file."""
    
    MolsSVGText = RDKitUtil.GetSVGForMolecules(ValidMols, OptionsInfo["NumOfMolsPerRow"], OptionsInfo["MolImageWidth"], OptionsInfo["MolImageHeight"], Legends = MolNames, AtomListsToHighlight = AtomLists, BondListsToHighlight = BondLists, BoldText =  OptionsInfo["FontBold"])
    
    MiscUtil.PrintInfo("\nGenerating SVG image file %s..." % OptionsInfo["Outfile"])
    
    OutFH = open(OptionsInfo["Outfile"], "w")
    OutFH.write(MolsSVGText)
    OutFH.close()
    
def GenerateImageFile(ValidMols, MolNames, AtomLists, BondLists):
    """Generate a non SVG image file."""
    
    Outfile = OptionsInfo["Outfile"]
    
    NumOfMolsPerRow = OptionsInfo["NumOfMolsPerRow"]
    Width = OptionsInfo["MolImageWidth"]
    Height = OptionsInfo["MolImageHeight"]
    
    # Setup drawing options...
    UpdatedDrawingOptions = DrawingOptions()
    UpdatedDrawingOptions.atomLabelFontSize = int(OptionsInfo["AtomLabelFontSize"])
    UpdatedDrawingOptions.bondLineWidth = float(OptionsInfo["BondLineWidth"])
    
    MolsImage = Draw.MolsToGridImage(ValidMols, molsPerRow = NumOfMolsPerRow,  subImgSize = (Width,Height), legends = MolNames, highlightAtomLists = AtomLists, highlightBondLists = BondLists, useSVG = False, kekulize = OptionsInfo["Kekulize"],  options = UpdatedDrawingOptions)
    
    MiscUtil.PrintInfo("\nGenerating image file %s..." % Outfile)
    
    if MiscUtil.CheckFileExt(Outfile, "pdf"):
        if MolsImage.mode == 'RGBA':
            MolsImage = MolsImage.convert('RGB')
    
    MolsImage.save(Outfile)
    
def GenerateHTMLTableFile(ValidMols, MolNames, HighlightAtomLists, HighlightBondLists):
    """Generate a HTML table file."""

    Outfile = OptionsInfo["Outfile"]
    
    Writer = open(Outfile, "w")
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
    
    MiscUtil.PrintInfo("\nGenerating HTML table file %s..." % Outfile)

    WriteHTMLPageHeader(Writer, len(ValidMols))
    WriteHTMLPageTitle(Writer)
    
    WriteHTMLTableHeader(Writer)
    WriteHTMLTableRows(Writer, ValidMols, MolNames, HighlightAtomLists, HighlightBondLists)
    WriteHTMLTableEnd(Writer)
    
    WriteHTMLPageFooter(Writer)
    WriteHTMLPageEnd(Writer)
    
    if Writer is not None:
        Writer.close()

def WriteHTMLTableRows(Writer, ValidMols, MolNames, HighlightAtomLists, HighlightBondLists):
    """Write out HTML table rows."""

    WriteTableHeaderRow(Writer, ValidMols)
    WriteTableDataRows(Writer, ValidMols, MolNames, HighlightAtomLists, HighlightBondLists)
    WriteTableFooterRow(Writer, ValidMols)

def WriteTableDataRows(Writer, ValidMols, MolNames, HighlightAtomLists, HighlightBondLists):
    """Write out table data row."""

    Writer.write("""        <tbody>\n""")

    MolCount = len(ValidMols)
    ColCount = GetColCount(MolCount)

    for Index in range(0, MolCount, ColCount):
        Writer.write("""          <tr>\n""")
    
        if OptionsInfo["CounterCol"]:
            Writer.write("""            <td></td>\n""")

        for MolIndex in range(Index, (Index + ColCount)):
            SetupStructureDataDrawing(Writer, MolIndex, ValidMols, MolNames, HighlightAtomLists, HighlightBondLists)
        
        Writer.write("""          </tr>\n""")
        
    Writer.write("""        </tbody>\n""")

def SetupStructureDataDrawing(Writer, MolIndex, Mols, MolNames, HighlightAtomLists, HighlightBondLists):
    """Setup structure data drawing for a tabel cell."""
    
    if MolIndex >= len(Mols):
        Writer.write("""            <td></td>\n""")
        return

    Mol = Mols[MolIndex]
    MolName = None if MolNames is None else MolNames[MolIndex]
    HighlightAtomList = None if HighlightAtomLists is None else HighlightAtomLists[MolIndex]
    HighlightBondList = None if HighlightBondLists is None else HighlightBondLists[MolIndex]
    
    SVGText = RDKitUtil.GetInlineSVGForMolecule(Mol, OptionsInfo["MolImageWidth"], OptionsInfo["MolImageHeight"], Legend = MolName, AtomListToHighlight = HighlightAtomList, BondListToHighlight = HighlightBondList, BoldText = OptionsInfo["FontBold"])

    PopoverTag = GetMolPopoverTag(Mol)
    ImageTag = "img" if PopoverTag is None else "img %s" % PopoverTag
    SVGInlineImageTag = "%s src=\"data:image/svg+xml;charset=UTF-8,\n%s\"" % (ImageTag, SVGText)
    
    Writer.write("""            <td bgcolor="white"><%s></td>\n""" % SVGInlineImageTag)

def WriteTableHeaderRow(Writer, ValidMols):
    """Write out table header row."""

    if not OptionsInfo["TableHeader"]:
        return
    
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
    
    # Write out rest of the column headers...
    MolCount = len(ValidMols)
    ColCount = GetColCount(MolCount)
    for ColIndex in range(0, ColCount):
        ColLabel = MiscUtil.GetExcelStyleColumnLabel(ColIndex + 1)
        Writer.write("""          <th>%s</th>\n""" % ColLabel)
        
    Writer.write("""        </tr>\n""")
    Writer.write("""      </thead>\n""")

def WriteTableFooterRow(Writer, ValidMols):
    """Write out table footer row."""

    if not OptionsInfo["TableFooter"]:
        return
    
    Writer.write("""      <tfoot>\n""")
    Writer.write("""        <tr>\n""")

    if OptionsInfo["CounterCol"]:
        Writer.write("""          <td></td>\n""")

    # Write out rest of the column headers...
    MolCount = len(ValidMols)
    ColCount = GetColCount(MolCount)
    for ColIndex in range(0, ColCount):
        ColLabel = MiscUtil.GetExcelStyleColumnLabel(ColIndex + 1)
        Writer.write("""          <td>%s</td>\n""" % ColLabel)
        
    Writer.write("""        </tr>\n""")
    Writer.write("""      </tfoot>\n""")

def WriteHTMLPageHeader(Writer, MolCount):
    """Write out HTML page header."""

    ColCount = GetColCount(MolCount)
    
    # Exclude counter and structure columns from sorting and searching...
    if OptionsInfo["CounterCol"]:
        ColIndicesList = ["0"]
        ColVisibilityExcludeColIndicesList = ["0"]
        ColIndexOffset = 1
        FreezeLeftColumns = "1"
    else:
        ColIndicesList = []
        ColVisibilityExcludeColIndicesList = []
        ColIndexOffset = 0

    MaxDataColVisColCount = 25
    for Index in range(0, ColCount):
        ColIndex = Index + ColIndexOffset
        ColIndicesList.append("%s" % ColIndex)
        
        if OptionsInfo["ColVisibility"]:
            if Index >= MaxDataColVisColCount:
                ColVisibilityExcludeColIndicesList.append("%s" %ColIndex)
    
    ColIndices = MiscUtil.JoinWords(ColIndicesList, ", ") if  len(ColIndicesList) else ""
    ColVisibilityExcludeColIndices = MiscUtil.JoinWords(ColVisibilityExcludeColIndicesList, ", ") if len(ColVisibilityExcludeColIndicesList) else ""
        
    DataColVisibilityExclude = False
    if OptionsInfo["ColVisibility"]:
        if ColCount > MaxDataColVisColCount:
            DataColVisibilityExclude = True
            MiscUtil.PrintWarning("The number of data columns, %d, is quite large. Only first %d data columns will be available in column visibility pulldown." % (ColCount, MaxDataColVisColCount))
            
    DisplayButtons = False
    if OptionsInfo["ColVisibility"]:
        if ColCount > 0:
            DisplayButtons = True
    
    FreezeCols = False
    if (OptionsInfo["CounterCol"] and OptionsInfo["ScrollX"]):
        FreezeCols = True
    
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
    
    # Start HTML header...
    Title = "Molecules table" if OptionsInfo["Header"] is None else OptionsInfo["Header"]
    
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

    if (OptionsInfo["Popover"]):
        Writer.write("""\
    <script type="text/javascript" language="javascript" src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js"></script>
    <script type="text/javascript" language="javascript" src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js"></script>

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
""")

    if (OptionsInfo["Popover"]):
        Writer.write("""\
    $('.MolPopover').popover();

""")
        
    Writer.write("""\
    var MolsTable = $('#MolsTable').DataTable( {
        "columnDefs": [
            {
                "orderable": false,
                "searchable": false,
                "targets": [%s]
            },
""" % (ColIndices))

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

    # Set up dom for displaying button and other options...
    if OptionsInfo["ColVisibility"]:
        if OptionsInfo["Paging"]:
            Writer.write("""\
        "dom":  "<'row'<'col-sm-6'l><'col-sm-6'<'float-right'B>>>" +
            "<'row'<'col-sm-12'tr>>" +
            "<'row'<'col-sm-5'i><'col-sm-7'p>>",
""")
        else:
            Writer.write("""\
        "dom":  "<'row'<'col'<'float-right'B>>>" +
            "<'row'<'col-sm-12'tr>>" +
            "<'row'<'col-sm-5'i><'col-sm-7'p>>",
""")
    else:
            Writer.write("""\
        "dom":  "<'row'<'col'l>>" +
            "<'row'<'col-sm-12'tr>>" +
            "<'row'<'col-sm-5'i><'col-sm-7'p>>",
""")
    
    #
    if OptionsInfo["ColVisibility"]:
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
                    return "Column " + (colIndex + 1);
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
        "lengthMenu": [ [5, 10, 15, 25, 50, 100, 500, 1000, -1], [5, 10, 15, 25, 50, 100, 500, 1000, "All"] ],
        "paging": %s,
        "pagingType": %s,
        "scrollX": %s,
        "scrollY": %s,
        "scrollCollapse": true,
        "order": [],
    } );
""" % (PageLength, Paging, PagingType, ScrollX, ScrollY))
    
    if OptionsInfo["CounterCol"]:
        Writer.write("""\
    MolsTable.on( 'order.dt search.dt', function () {
        MolsTable.column(0, {search:'applied', order:'applied'}).nodes().each( function (cell, rowIndex) {
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
        Writer.write("""\n    <table id="MolsTable" cellspacing="0" width="100%">\n""")
    else:
        Writer.write("""    <table id="MolsTable" class="%s" cellspacing="0" width="100%s">\n""" % (OptionsInfo["TableStyle"], "%"))
        
def WriteHTMLTableEnd(Writer):
    """Write out HTML table end."""

    Writer.write("""    </table>\n\n""")

def GetColCount(MolCount):
    """Get tabke column count."""
    
    ColCount = OptionsInfo["NumOfMolsPerRow"] if OptionsInfo["NumOfMolsPerRow"] <= MolCount else MolCount
    
    return ColCount

def SetupAtomListsToHighlight(ValidMols):
    """Set up atom lists to highlight using specified SMARTS pattern."""

    AtomListsToHighlight = None
    if OptionsInfo["HighlightSMARTSPattern"] is None:
        return  AtomListsToHighlight
    
    PatternMol = Chem.MolFromSmarts(OptionsInfo["HighlightSMARTSPattern"])
    AtomListsToHighlight = []
    for ValidMol in ValidMols:
        # Get matched atom lists and flatten it...
        MatchedAtomsLists = ValidMol.GetSubstructMatches(PatternMol)
        MatchedAtoms = [ Atom for AtomsList in MatchedAtomsLists for Atom in AtomsList] 
        AtomListsToHighlight.append(MatchedAtoms)
    
    return  AtomListsToHighlight

def PerformAlignment(ValidMols):
    """Perform alignment to a common template specified by a SMARTS pattern."""
    
    if OptionsInfo["AlignmentSMARTSPattern"] is None:
        return
    
    PatternMol = Chem.MolFromSmarts(OptionsInfo["AlignmentSMARTSPattern"])
    AllChem.Compute2DCoords(PatternMol)
        
    MatchedValidMols = [ValidMol for ValidMol in ValidMols if ValidMol.HasSubstructMatch(PatternMol)]
    for ValidMol in MatchedValidMols:
        AllChem.GenerateDepictionMatching2DStructure(ValidMol, PatternMol)

def GetMolPopoverTag(Mol):
    """Set up a popover window containing any additional information about molecule."""
    
    if not OptionsInfo["Popover"]:
        return None
    
    # Set up data label and values...
    AvailableDataLabels = Mol.GetPropNames(includePrivate = False, includeComputed = False)
    
    DataContentLines = []
    MaxDataCharWidth = OptionsInfo["PopoverTextWidth"]
    MaxDataDisplayCount = OptionsInfo["PopoverDataCount"]

    DataDisplayCount = 0
    SkippedDataDisplay = False
    for DataLabel in AvailableDataLabels:
        DataDisplayCount += 1
        if DataDisplayCount > MaxDataDisplayCount:
            SkippedDataDisplay = True
            break
        
        DataValue = "%s" % Mol.GetProp(DataLabel)
        DataValue = DataValue.strip()
        if MiscUtil.IsEmpty(DataValue):
            continue

        # Change any new lines to ;
        if re.search("(\r\n|\r|\n)", DataValue):
            DataValue = re.sub("(\r\n|\r|\n)", "; ", DataValue)
        
        DataValue = MiscUtil.TruncateText(DataValue, MaxDataCharWidth, "...")
        DataValue = MiscUtil.ReplaceHTMLEntitiesInText(DataValue)
        
        DataContent = "<b>%s</b>: %s" % (DataLabel, DataValue)
        DataContentLines.append(DataContent)

    if not len(DataContentLines):
        return None

    if SkippedDataDisplay:
        DataContent = "<b>... ... ...</b>"
        DataContentLines.append(DataContent)
        
        DataContent = "Showing 1 to %s of %s" % (MaxDataDisplayCount, len(AvailableDataLabels))
        DataContentLines.append(DataContent)
    else:
        DataContent = "Showing 1 to %s of %s" % (DataDisplayCount, len(AvailableDataLabels))
        DataContentLines.append(DataContent)
        
    DataContent = MiscUtil.JoinWords(DataContentLines, "<br/>")
    PopoverTag = """class="MolPopover" data-toggle="popover" data-html="true" data-trigger="click" data-placement="right" title="<span class='small'><b>Additional Information</b></span>" data-content="<span class='small'>%s</span>" """ % DataContent
    
    return PopoverTag

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
    
    OptionsInfo["AtomLabelFontSize"] = Options["--atomLabelFontSize"]
    OptionsInfo["BondLineWidth"] = Options["--bondLineWidth"]
    
    Compute2DCoords = True
    if re.match("^yes$", Options["--compute2DCoords"], re.I):
        Compute2DCoords = True
    elif re.match("^no$", Options["--compute2DCoords"], re.I):
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
    
    OptionsInfo["FontBold"] = True
    if re.match("^no$", Options["--fontBold"], re.I):
        OptionsInfo["FontBold"] = False
        
    Footer = None
    if not re.match("^None$", Options["--footer"], re.I):
        Footer = Options["--footer"]
    OptionsInfo["Footer"]  = Footer

    FooterClass = Options["--footerClass"].strip()
    if MiscUtil.IsEmpty(FooterClass):
        MiscUtil.PrintError("The value specified using option \"--footerClass\" is empty.")
    OptionsInfo["FooterClass"]  = FooterClass
    
    Header = None
    if not re.match("^None$", Options["--header"], re.I):
        Header = Options["--header"]
    OptionsInfo["Header"]  = Header
    
    HeaderStyle = Options["--headerStyle"].strip()
    if MiscUtil.IsEmpty(HeaderStyle):
        MiscUtil.PrintError("The value specified using option \"--headerStyle\" is empty.")
    OptionsInfo["HeaderStyle"]  = HeaderStyle
    
    HighlightSMARTSPattern = None
    if not re.match("^None$", Options["--highlightSMARTS"], re.I):
        HighlightSMARTSPattern = Options["--highlightSMARTS"]
    OptionsInfo["HighlightSMARTSPattern"]  = HighlightSMARTSPattern
    
    OptionsInfo["Kekulize"] = True
    if re.match("^no$", Options["--kekulize"], re.I):
        OptionsInfo["Kekulize"] = False
        
    OptionsInfo["KeysNavigation"] = True
    if re.match("^no$", Options["--keysNavigation"], re.I):
        OptionsInfo["KeysNavigation"] = False
    
    SizeValues = Options["--molImageSize"].split(",")
    OptionsInfo["MolImageWidth"] = int(SizeValues[0])
    OptionsInfo["MolImageHeight"] = int(SizeValues[1])

    OptionsInfo["NumOfMolsPerRow"] = int(Options["--numOfMolsPerRow"])

    OptionsInfo["Paging"] = True
    if re.match("^no$", Options["--paging"], re.I):
        OptionsInfo["Paging"] = False
    
    PagingType = Options["--pagingType"]
    if not re.match("^(numbers|simple|simple_numbers|full|full_numbers|simple_number)$", Options["--pagingType"], re.I):
        MiscUtil.PrintWarning("The paging type name, %s, specified using option \"--pagingType\" appears to be a unknown type..." % (PagingType))
    OptionsInfo["PagingType"] = PagingType.lower()
    
    OptionsInfo["PageLength"] = int(Options["--pageLength"])
    
    OptionsInfo["Popover"] = True
    if re.match("^no$", Options["--popover"], re.I):
        OptionsInfo["Popover"] = False
    OptionsInfo["PopoverDataCount"] = int(Options["--popoverDataCount"])
    OptionsInfo["PopoverTextWidth"] = int(Options["--popoverTextWidth"])
    
    OptionsInfo["ShowMolName"] = True
    if re.match("^no$", Options["--showMolName"], re.I):
        OptionsInfo["ShowMolName"] = False
    
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

    OptionsInfo["TableFooter"] = True
    if re.match("^no$", Options["--tableFooter"], re.I):
        OptionsInfo["TableFooter"] = False
    
    OptionsInfo["TableHeader"] = True
    if re.match("^no$", Options["--tableHeader"], re.I):
        OptionsInfo["TableHeader"] = False
    
    TableHeaderStyle = None
    if not re.match("^None$", Options["--tableHeaderStyle"], re.I):
        TableHeaderStyle = Options["--tableHeaderStyle"]
        TableHeaderStyle = TableHeaderStyle.lower()
        CheckOptionTableClassColorValues("--tableHeaderStyle", [TableHeaderStyle])
    OptionsInfo["TableHeaderStyle"]  = TableHeaderStyle

def CheckOptionTableClassColorValues(OptionName, ColorsList):
    """Check names of table color classes and issue a warning for unknown names."""

    TableClassColors = ["thead-dark", "thead-light", "table-primary", "table-success", "table-danger", "table-info", "table-warning", "table-active", "table-secondary", "table-light", "table-dark", "bg-primary", "bg-success", "bg-danger",  "bg-info", "bg-warning", "bg-secondary", "bg-dark", "bg-light"]

    for Color in ColorsList:
        if not Color in TableClassColors:
            MiscUtil.PrintWarning("The color class name, %s, specified using option \"%s\" appears to be a unknown name..." % (Color, OptionName))
        
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
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi csv tsv txt")

    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    
    if not re.match("^None$", Options["--alignmentSMARTS"], re.I):
        PatternMol = Chem.MolFromSmarts(Options["--alignmentSMARTS"])
        if PatternMol is None:
            MiscUtil.PrintError("The value specified, %s, using option \"--alignmentSMARTS\" is not a valid SMARTS: Failed to create pattern molecule" % Options["--alignmentSMARTS"])
    
    MiscUtil.ValidateOptionIntegerValue("--atomLabelFontSize", Options["--atomLabelFontSize"], {">": 0})
    MiscUtil.ValidateOptionFloatValue("-b, --bondLineWidth", Options["--bondLineWidth"], {">": 0.0})
    
    MiscUtil.ValidateOptionTextValue("--compute2DCoords", Options["--compute2DCoords"], "yes no auto")
    
    MiscUtil.ValidateOptionTextValue("--counterCol", Options["--counterCol"], "yes no")
    MiscUtil.ValidateOptionTextValue("--colVisibility", Options["--colVisibility"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--f, -fontBold", Options["--fontBold"], "yes no")
    
    if not re.match("^None$", Options["--highlightSMARTS"], re.I):
        PatternMol = Chem.MolFromSmarts(Options["--highlightSMARTS"])
        if PatternMol is None:
            MiscUtil.PrintError("The value specified, %s, using option \"--highlightSMARTS\" is not a valid SMARTS: Failed to create pattern molecule" % Options["--highlightSMARTS"])
    
    MiscUtil.ValidateOptionTextValue("--kekulize", Options["--kekulize"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("-k, --keysNavigation", Options["--keysNavigation"], "yes no")
    
    MiscUtil.ValidateOptionNumberValues("-m, --molImageSize", Options["--molImageSize"], 2, ",", "integer", {">": 0})
    MiscUtil.ValidateOptionIntegerValue("--numOfMolsPerRow", Options["--numOfMolsPerRow"], {">": 0})
    
    MiscUtil.ValidateOptionTextValue("-p, --paging", Options["--paging"], "yes no")
    MiscUtil.ValidateOptionIntegerValue("--pageLength", Options["--pageLength"], {">": 0})
    
    MiscUtil.ValidateOptionTextValue("--popover", Options["--popover"], "yes no")
    MiscUtil.ValidateOptionIntegerValue("--popoverDataCount", Options["--popoverDataCount"], {">": 0})
    MiscUtil.ValidateOptionIntegerValue("--popoverTextWidth", Options["--popoverTextWidth"], {">": 0})
    
    MiscUtil.ValidateOptionTextValue("--showMolName", Options["--showMolName"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--scrollX", Options["--scrollX"], "yes no")
    MiscUtil.ValidateOptionTextValue("--scrollY", Options["--scrollY"], "yes no")
    if not re.search("vh$", Options["--scrollYSize"], re.I):
        MiscUtil.ValidateOptionIntegerValue("--scrollYSize", Options["--scrollYSize"], {">": 0})

    MiscUtil.ValidateOptionTextValue("--tableFooter", Options["--tableFooter"], "yes no")
    MiscUtil.ValidateOptionTextValue("--tableHeader", Options["--tableHeader"], "yes no")
    
# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitDrawMolecules.py - Draw molecules and generate an image or HTML file

Usage:
    RDKitDrawMolecules.py [--alignmentSMARTS <SMARTS>] [--atomLabelFontSize <number>]
                             [--bondLineWidth <number>] [--compute2DCoords <yes | no>] [--counterCol <yes or no>]
                             [--colVisibility <yes or no>] [--fontBold <yes or no>] [--footer <text>] [--footerClass <text>] 
                             [--header <text>] [--headerStyle <text>] [--highlightSMARTS <SMARTS>]
                             [--infileParams <Name,Value,...>] [--kekulize <yes or no>] [--keysNavigation <yes or no>]
                             [--molImageSize <width,height>] [--numOfMolsPerRow <number>] [--overwrite]
                             [--paging <yes or no>] [--pagingType <numbers, simple, ...>] [--pageLength <number>]
                             [--popover <yes or no>] [--popoverDataCount <number>] [--popoverTextWidth <number>]
                             [--showMolName <yes or no>] [--scrollX <yes or no>] [--scrollY <yes or no>]
                             [--scrollYSize <number>] [--tableFooter <yes or no>] [--tableHeader <yes or no>]
                             [--tableHeaderStyle <thead-dark,thead-light,...>]
                             [--tableStyle <table,table-striped,...>] [-w <dir>] -i <infile> -o <outfile>
    RDKitDrawMolecules.py -h | --help | -e | --examples

Description:
    Draw molecules in a grid and write them out as an image file or a HTML table file. The
    SVG image or HTML table file appears to be the best among all the available image file
    options, as rendered in a browser. The Python modules aggdraw/cairo are required to
    generate high quality PNG images.
    
    The drawing of the molecules are embedded in HTML table columns as in line SVG
    images. The HTML table is an interactive table and requires internet access for viewing
    in a browser. It employs he following frameworks: JQuery, Bootstrap, and DataTable.
    
    The options '--atomLabelFontSize' and '--bondLineWidth' don't appear
    to work during the generation of a SVG image.

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd), SMILES (.smi,
    .txt, .csv, .tsv)

    The output image file can be saved in any format supported by the Python Image
    Library (PIL). The image format is automatically detected from the output file extension.

    Some of the most common output image file formats are: GIF (.gif), JPEG (.jpg),
    PNG (.png), SVG (.svg), TIFF (.tif). In addition, a HTML (.html) file format
    containing a table is supported.

Options:
    -a, --alignmentSMARTS <SMARTS>  [default: none]
        SMARTS pattern for aligning molecules to a common template.
    --atomLabelFontSize <number>  [default: 12]
        Font size for drawing atom labels. This option is ignored during the generation of
        a SVG and HTML output file.
    -b, --bondLineWidth <number>  [default: 1.2]
        Line width for drawing bonds. This option is ignored during the generation of a SVG
        and HTML output file.
    -c, --compute2DCoords <yes or no>  [default: auto]
        Compute 2D coordinates of molecules before drawing. Default: yes for all file
        formats.
    --counterCol <yes or no>  [default: yes]
        Show a counter column as the first column in the table. It contains the position
        for each row in the HTML table. This option is only used during the generation of
        a HTML table file.
    --colVisibility <yes or no>  [default: yes]
        Show a dropdown button to toggle visibility of columns in the table. This option is
        only used during the generation of a HTML table file.
    -e, --examples
        Print examples.
    -f --fontBold <yes or no>  [default: yes]
        Make all text fonts bold during the generation of  a SVG and HTML output file. This
        option is ignored for all other output files.
    --footer <text>  [default: none]
        Footer text to insert at the bottom of the HTML page after the table. This option is
        only used during the generation of a HTML table file.
    --footerClass <text>  [default: small text-center text-muted]
        Footer class style to use with <p> tag. This option is only used during the
        generation of a HTML table file.
    -h, --help
        Print this help message.
    --header <text>  [default: none]
        Header text to insert at the top of the HTML page before the table. This option is
        only used during the generation of a HTML table file.
    --headerStyle <text>  [default: h5]
        Header style to use. Possible values: h1 to h6. This option is only used during the
        generation of a HTML table file.
    --highlightSMARTS <SMARTS>  [default: none]
        SMARTS pattern for highlighting atoms and bonds in molecules. All matched
        substructures are highlighted.
    -i, --infile <infile>
        Input file name.
    --infileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for reading
        molecules from files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD, MOL: removeHydrogens,yes,sanitize,yes,strictParsing,yes
            SMILES: smilesColumn,1,smilesNameColumn,2,smilesDelimiter,space,
                smilesTitleLine,auto,sanitize,yes
            
        Possible values for smilesDelimiter: space, comma or tab.
    -k, --kekulize <yes or no>  [default: yes]
        Perform kekulization on molecules. This option is ignored during the generation of
        a SVG and HTML output file.
    --keysNavigation <yes or no>  [default: yes]
        Provide Excel like keyboard cell navigation for the table. This option is only used
        during the generation of a HTML table file.
    -m, --molImageSize <width,height>  [default: 250,200]
        Image size of a molecule in pixels.
    -n, --numOfMolsPerRow <number>  [default: 4]
        Number of molecules to draw in a row.
    -o, --outfile <outfile>
        Output file name.
    --overwrite
        Overwrite existing files.
    -p, --paging <yes or no>  [default: yes]
        Provide page navigation for browsing data in the table. This option is only used
        during the generation of a HTML table file.
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
            
        This option is only used during the generation of a HTML table file.
    --pageLength <number>  [default: 5]
        Number of rows to show per page. This option is only used during the
        generation of a HTML table file.
    --popover <yes or no>  [default: yes]
        Display a popover window containing additional information about the
        molecule. The popover is opened after a click on the drawing of a
        molecule. A subsequent click on the same drawing closes the popover.
        This option is only used during the generation of a HTML table file.
    --popoverDataCount <number>  [default: 25]
        Maximum number of data fields to show in a popover window. This option is
        only used during the generation of a HTML table file.
    --popoverTextWidth <number>  [default: 50]
        Maximum width in characters for text display in a popover window before
        truncating the text. This option is only used during the generation of a HTML
        table file.
    -s, --showMolName <yes or no>  [default: yes]
        Show molecule names under the images.This option is only used during the
        generation of a HTML table file.
    --scrollX <yes or no>  [default: yes]
        Provide horizontal scroll bar in the table as needed.This option is only used
        during the generation of a HTML table file.
    --scrollY <yes or no>  [default: yes]
        Provide vertical scroll bar in the table as needed.This option is only used during
        the generation of a HTML table file.
    --scrollYSize <number>  [default: 75vh]
        Maximum height of table viewport either in pixels or percentage of the browser
        window height before providing a vertical scroll bar. Default: 75% of the height of
        browser window.This option is only used during the generation of a HTML table file.
    -t, --tableStyle <table,table-striped,...>  [default: table,table-hover,table-sm]
        Style of table. Possible values: table, table-striped, table-bordered,
        table-hover, table-dark, table-sm, none, or All. Default: 'table,table-hover'. A
        comma delimited list of any valid Bootstrap table styles is also supported
        
        This option is only used during the generation of a HTML table file.
    --tableFooter <yes or no>  [default: yes]
        Show Excel style column headers at the end of  the table. This option is only
        used during the generation of a HTML table file.
    --tableHeader <yes or no>  [default: yes]
        Show Excel style column headers in the table. This option is only used
        during the generation of a HTML table file.
    --tableHeaderStyle <thead-dark,thead-light,...>  [default: thead-dark]
        Style of table header. Possible values: thead-dark, thead-light, or none.
        The names of the following contextual color classes are also supported:
        table-primary (Blue), table-success (Green), table-danger (Red), table-info
        (Light blue), table-warning (Orange), table-active (Grey), table-light (Light
        grey), and  table-dark (Dark grey).
        
        This option is only used during the generation of a HTML table file.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To automatically compute 2D coordinates for molecules in a SMILES file and
    generate a SVG image file containing 4 molecules per row in a grid with cell
    size of 250 x 200 pixels, type:

        % RDKitDrawMolecules.py -i Sample.smi -o SampleOut.svg

    To automatically compute 2D coordinates for molecules in a SMILES file and
    generate a SVG image file containing 2 molecules per row in a grid with cell
    size of 400 x 300 pixels and without any keulization along with highlighting
    a specific set of atoms and bonds indicated by a SMARTS pattern, type:

        % RDKitDrawMolecules.py -n 2 -m "400,300" -k no --fontBold no
          --highlightSMARTS  'c1ccccc1' -i Sample.smi -o SampleOut.svg

    To generate a PNG image file for molecules in a SD file using existing 2D
    coordinates, type

        % RDKitDrawMolecules.py --compute2DCoords no -i Sample.sdf
          -o SampleOut.png

    To automatically compute 2D coordinates for molecules in a SD file and
    generate a HTML file containing 4 molecules per row in a table, along with
    all the bells and whistles to interact with the table, type:

        % RDKitDrawMolecules.py -i Sample.sdf -o SampleOut.html

    To automatically compute 2D coordinates for molecules in a SD file and
    generate a HTML file containing 4 molecules per row in a table without
    any bells and whistles to interact with the table, type:

        % RDKitDrawMolecules.py --counterCol no --colVisibility no
          --keysNavigation no --paging  no --popover no --scrollX no
          --scroll no --tableFooter no --tableHeader  no -i Sample.sdf
          -o SampleOut.html

    To automatically compute 2D coordinates for molecules in a CSV SMILES file
    with column headers, SMILES strings in column 1, and name in column 2 and
    generate a PDF image file, type:

        % RDKitDrawMolecules.py --infileParams "smilesDelimiter,comma,
          smilesTitleLine,yes,smilesColumn,1,smilesNameColumn,2"
          -i SampleSMILES.csv -o SampleOut.pdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitDrawMoleculesAndDataTable.py, RDKitRemoveDuplicateMolecules.py,
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
