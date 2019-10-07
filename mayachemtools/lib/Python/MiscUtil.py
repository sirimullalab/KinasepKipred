#
# File: MiscUtil.py
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

from __future__ import print_function

import os
import sys
import time
import re
import csv
import textwrap

__all__ = ["CheckFileExt", "CheckTextValue", "DoesSMILESFileContainTitleLine", "GetExamplesTextFromDocOptText",  "GetExcelStyleColumnLabel", "GetMayaChemToolsLibDataPath", "GetTextLinesWords", "GetWallClockAndProcessorTime", "GetFormattedElapsedTime", "IsEmpty", "IsFloat", "IsInteger", "IsNumber", "JoinWords", "ParseFileName", "PrintError", "PrintInfo", "PrintWarning",  "ProcessOptionInfileParameters", "ProcessOptionOutfileParameters", "ReplaceHTMLEntitiesInText", "ValidateOptionsDistinctFileNames", "ValidateOptionFileExt", "ValidateOptionFilePath", "ValidateOptionFloatValue", "ValidateOptionIntegerValue", "ValidateOptionNumberValue", "ValidateOptionNumberValues", "ValidateOptionsOutputFileOverwrite", "ValidateOptionTextValue",  "TruncateText", "WrapText"]

def CheckFileExt(FileName, FileExts):
    """Check file type based on the specified file extensions delimited by spaces.
    
    Arguments:
        FileName (str): Name of a file.
        FileExts (str): Space delimited string containing valid file extensions.

    Returns:
        bool : True, FileName contains a valid file extension; Otherwise, False.

    """
    
    for FileExt in FileExts.split():
        if re.search(r"\.%s$" % FileExt, FileName, re.IGNORECASE):
            return True
    
    return False

def CheckTextValue(Value, ValidValues):
    """Check text value based on the specified valid values delimited by spaces.

    Arguments:
        Value (str): Text value
        ValidValues (str): Space delimited string containing valid values.

    Returns:
        bool : True, Value is valid; Otherwise, False.

    """
    
    ValidValues = re.sub(' ', '|', ValidValues)
    if re.match("^(%s)$" % ValidValues, Value, re.IGNORECASE):
        return True
    
    return False

def GetTextLinesWords(TextFilePath, Delimiter, QuoteChar, IgnoreHeaderLine):
    """Parse lines in the specified text file into words in a line and return a list containing
    list of parsed line words.

    Arguments:
        TextFilePath (str): Text file name including file path.
        Delimiter (str): Delimiter for parsing text lines.
        QuoteChar (str): Quote character for line words.
        IgnoreHeaderLine (bool): A flag indicating whether to ignore first
            valid data line corresponding to header line.

    Returns:
        list : A list of lists containing parsed words for lines.

    Notes:
        The lines starting with # or // are considered comment lines and are
        ignored during parsing along with any empty lines.

    """
    TextFile = open(TextFilePath, "r")
    if TextFile is None:
        PrintError("Couldn't open text file: %s.\n" % (TextFilePath))

    # Collect text lines...
    TextLines = []
    FirstValidLine = True
    for Line in TextFile:
        Line = Line.strip()
        
        # Ignore empty lines...
        if not len(Line):
            continue
        
        # Ignore comments...
        if re.match("^(#|\/\/)", Line, re.I):
            continue

        # Ignore header line...
        if FirstValidLine:
            FirstValidLine = False
            if IgnoreHeaderLine:
                continue
        
        TextLines.append(Line)
        
    TextFile.close()

    # Parse text lines...
    TextLinesWords = []
    
    TextLinesReader = csv.reader(TextLines, delimiter = Delimiter, quotechar = QuoteChar)
    for LineWords in TextLinesReader:
        TextLinesWords.append(LineWords)
    
    return TextLinesWords
    
def DoesSMILESFileContainTitleLine(FileName):
    """Determine whether the SMILES file contain a title line based on the  presence
    of a string SMILES, Name or ID in the first line.

    Arguments:
        FileName (str): Name of a file.

    Returns:
        bool : True, File contains title line; Otherwise, False.

    """
    
    Infile = open(FileName, "r")
    if Infile is None:
        return False

    Line = Infile.readline()
    Infile.close()

    if re.search("(SMILES|Name|ID)", Line, re.I):
        return True
        
    return False
    
def GetExamplesTextFromDocOptText(DocOptText):
    """Get script usage example lines from a docopt doc string. The example text
    line start from a line containing `Examples:`  keyword at the beginning of the line.
    
    Arguments:
        DocOptText (str): Doc string containing script usage examples lines starting with
            a line marked by `Examples:` keyword at the beginning of a line.

    Returns:
        str : A string containing text lines retrieved from the examples section of
            DocOptText parameter.

    """
    
    ExamplesStart = re.compile("^Examples:", re.IGNORECASE)
    ExamplesEnd = re.compile("^(Author:|See also:|Copyright:)", re.IGNORECASE)
    
    ExamplesText = 'Examples text is not available';
    ExamplesTextFound = False
    
    for Line in DocOptText.splitlines():
        if ExamplesStart.match(Line):
            ExamplesText = 'Examples:';
            ExamplesTextFound = True
            continue
        
        if ExamplesEnd.match(Line):
            break
        
        if ExamplesTextFound:
            ExamplesText += "\n" + Line;
    
    return ExamplesText

def GetExcelStyleColumnLabel(ColNum):
    """Return Excel style column label for a colum number.
    
    Arguments:
        ColNum (int): Column number

    Returns:
        str : Excel style column label.

    """
    Letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    ColLabelList = []
    while ColNum:
        ColNum, SubColNum = divmod(ColNum - 1, 26)
        ColLabelList[:0] = Letters[SubColNum]
    
    return ''.join(ColLabelList)
    
def GetWallClockAndProcessorTime():
    """Get wallclock and processor times in seconds.
    
    Returns:
        float : Wallclock time.
        float : Processor time.

    """
    return (time.time(), time.clock())

def GetMayaChemToolsLibDataPath():
    """Get location of MayaChemTools lib data directory.
    
    Returns:
        str : Location of MayaChemTools lib data directory.

    Notes:
        The location of MayaChemTools lib data directory is determined relative to
        MayaChemTools python lib directory name available through sys.path.

    """
    MayaChemToolsDataPath = ""
    
    for PathEntry in sys.path:
        if re.search("MayaChemTools", PathEntry, re.I) and re.search("Python", PathEntry, re.I):
            MayaChemToolsDataPath = os.path.join( PathEntry, "..",  "data")
            break
        else:
            PrintInfo("PathEntry didn't match")
    
    if not len(MayaChemToolsDataPath):
        PrintWarning("MayaChemTools lib directory location doesn't appear to exist in system search path specified by sys.path...")
        
    return MayaChemToolsDataPath
    
def GetFormattedElapsedTime(StartingWallClockTime, StartingProcessorTime):
    """Get elapsed wallclock and processor times  as a string in the following
    format: %d wallclock secs ( %.2f process secs).
    
    Arguments:
        StartingWallClockTime (float): Starting wallclock time in seconds.
        StartingProcessorTime (float): Starting processor time in seconds.

    Returns:
        str : Elapsed time formatted as: %d wallclock secs ( %.2f process secs)

    """
    
    ElapsedWallClockTime = time.time() - StartingWallClockTime
    ElapsedProcessorTime = time.clock() - StartingProcessorTime
    
    ElapsedTime = "%d wallclock secs ( %.2f process secs)" % (ElapsedWallClockTime, ElapsedProcessorTime)
    
    return ElapsedTime

def IsEmpty(Value):
    """Determine whether the specified value is empty after converting
    it in to a string and removing all leading and trailing white spaces. A  value
    of type None is considered empty.
    
    Arguments:
        Value (str, int or float): Text or a value

    Returns:
        bool : True, Text string is empty; Otherwsie, False.

    """

    if Value is None:
        return True

    TextValue = "%s" % Value
    TextValue = TextValue.strip()

    return False if len(TextValue) else True

def IsFloat(Value):
    """Determine whether the specified value is a float by converting it
    into a float.
    
    Arguments:
        Value (str, int or float): Text

    Returns:
        bool : True, Value is a float; Otherwsie, False.

    """

    return IsNumber(Value)

def IsInteger(Value):
    """Determine whether the specified value is an integer by converting it
    into an int.
    
    Arguments:
        Value (str, int or float): Text

    Returns:
        bool : True, Value is an integer; Otherwsie, False.

    """

    Status = True
    
    if Value is None:
        return False

    try:
        Value = int(Value)
        Status = True
    except ValueError:
        Status = False
    
    return Status

def IsNumber(Value):
    """Determine whether the specified value is a number by converting it
    into a float.
    
    Arguments:
        Value (str, int or float): Text

    Returns:
        bool : True, Value is a number; Otherwsie, False.

    """

    Status = True
    
    if Value is None:
        return Status

    try:
        Value = float(Value)
        Status = True
    except ValueError:
        Status = False
    
    return Status

def JoinWords(Words, Delimiter, Quote = False):
    """Join words in a list using specified delimiter with optional quotes around words.
    
    Arguments:
        Words (list): List containing words to join.
        Delimiter (string): Delimiter for joining words.
        Quote (boolean): Put quotes around words.

    Returns:
        str : String containing joined words.

    """
    
    if Quote:
        JoinedWords = Delimiter.join('"{0}"'.format(Word) for Word in Words)
    else:
        JoinedWords = Delimiter.join(Words)
        
    return JoinedWords
    
def ParseFileName(FilePath):
    """Parse specified file path and return file dir, file name, and file extension.
    
    Arguments:
        FilePath (str): Name of a file with complete file path.

    Returns:
        str : File directory.
        str : File name without file extension.
        str : File extension.

    """
    FileDir, FileBaseName = os.path.split(FilePath)
    FileName, FileExt = os.path.splitext(FileBaseName)
    
    if re.match("^\.", FileExt):
        FileExt = re.sub("^\.", "", FileExt)
        
    return (FileDir, FileName, FileExt)
    
def PrintError(Msg, Status=2):
    """Print message to stderr along with flushing stderr and exit with a specified
    status. An `Error` prefix is placed before the message.
    
    Arguments:
        Msg (str): Text message.
        Status (int): Exit status.

    """
    
    PrintInfo("Error: %s" % Msg)
    sys.exit(Status)

def PrintInfo(Msg=''):
    """Print message to stderr along with flushing stderr.
    
    Arguments:
        Msg (str): Text message.

    """
    
    print(Msg, sep=' ', end='\n', file=sys.stderr)
    sys.stderr.flush()

def PrintWarning(msg):
    """Print message to stderr along with flushing stderr. An `Warning` prefix
    is placed before the message.
    
    Arguments:
        Msg (str): Text message.

    """
    
    PrintInfo("Warning: %s" % msg)

def ValidateOptionFileExt(OptionName, FileName, FileExts):
    """Validate file type based on the specified file extensions delimited by spaces.
    
    Arguments:
        OptionName (str): Command line option name.
        FileName (str): Name of a file.
        FileExts (str): Space delimited string containing valid file extensions.

    Notes:
        The function exits with an error message for a file name containing
        invalid file extension.

    """
    
    if not CheckFileExt(FileName, FileExts):
        PrintError("The file name specified , %s, for option \"%s\" is not valid. Supported file formats: %s\n" % (FileName, OptionName, FileExts))

def ValidateOptionFilePath(OptionName, FilePath):
    """Validate presence of the file.
    
    Arguments:
        OptionName (str): Command line option name.
        FilePath (str): Name of a file with complete path.

    Notes:
        The function exits with an error message for a file path that doesn't exist.

    """
    
    if not os.path.exists(FilePath):
        PrintError("The file specified, %s, for option \"%s\" doesn't exist.\n" % (FilePath, OptionName))

def ValidateOptionFloatValue(OptionName, OptionValue, CmpOpValueMap):
    """Validate option value using comparison operater and value pairs in specified in
    a map.
    
    Arguments:
        OptionName (str): Command line option name.
        OptionValue (float or str): Command line option value.
        CmpOpValueMap (dictionary): Comparison operator key and value pairs to
            validate values specified in OptionValue.

    Notes:
        The function exits with an error message for an invalid option values specified
        in OptionValue.

    Examples:

        ValidateOptionNumberValue("-b, --butinaSimilarityCutoff", 
            Options["--butinaSimilarityCutoff"],
            {">": 0.0, "<=" : 1.0})

    """

    if not IsFloat(OptionValue):
        PrintError("The value specified, %s, for option \"%s\" must be a float." % (OptionValue, OptionName))
    
    return ValidateOptionNumberValue(OptionName, float(OptionValue), CmpOpValueMap)

def ValidateOptionIntegerValue(OptionName, OptionValue, CmpOpValueMap):
    """Validate option value using comparison operater and value pairs in specified in
    a map.
    
    Arguments:
        OptionName (str): Command line option name.
        OptionValue (int or str): Command line option value.
        CmpOpValueMap (dictionary): Comparison operator key and value pairs to
            validate values specified in OptionValue.

    Notes:
        The function exits with an error message for an invalid option values specified
        in OptionValue.

    Examples:

        ValidateOptionIntegerValue("--maxConfs", Options["--maxConfs"],
            {">": 0})

    """

    if not IsInteger(OptionValue):
        PrintError("The value specified, %s, for option \"%s\" must be an integer." % (OptionValue, OptionName))
    
    return ValidateOptionNumberValue(OptionName, int(OptionValue), CmpOpValueMap)

def ValidateOptionNumberValue(OptionName, OptionValue, CmpOpValueMap):
    """Validate option value using comparison operater and value pairs in specified in
    a map.
    
    Arguments:
        OptionName (str): Command line option name.
        OptionValue (int or float): Command line option value.
        CmpOpValueMap (dictionary): Comparison operator key and value pairs to
            validate values specified in OptionValue.

    Notes:
        The function exits with an error message for an invalid option values specified
        in OptionValue.

    Examples:

        ValidateOptionNumberValue("--maxConfs", int(Options["--maxConfs"]),
            {">": 0})
        ValidateOptionNumberValue("-b, --butinaSimilarityCutoff", 
            float(Options["--butinaSimilarityCutoff"]),
            {">": 0.0, "<=" : 1.0})

    """
    
    Status = True
    for CmpOp in CmpOpValueMap:
        Value = CmpOpValueMap[CmpOp]
        if re.match("^>$", CmpOp, re.I):
            if OptionValue <= Value:
                Status = False
                break
        elif re.match("^>=$", CmpOp, re.I):
            if OptionValue < Value:
                Status = False
                break
        elif re.match("^<$", CmpOp, re.I):
            if OptionValue >= Value:
                Status = False
                break
        elif re.match("^<=$", CmpOp, re.I):
            if OptionValue > Value:
                Status = False
                break
        else:
            PrintError("The specified comparison operator, %s, for function MiscUtil.ValidateOptionNumberValue is not supported\n" % (CmpOp))
    
    if not Status:
        FirstValue = True
        SupportedValues = ""
        for CmpOp in CmpOpValueMap:
            Value = CmpOpValueMap[CmpOp]
            if FirstValue:
                FirstValue = False
                SupportedValues = "%s %s" % (CmpOp, Value)
            else:
                SupportedValues = "%s and %s %s" % (SupportedValues, CmpOp, Value)
        
        PrintError("The value specified, %s, for option \"%s\" is not valid. Supported value(s): %s " % (OptionValue, OptionName, SupportedValues))

def ValidateOptionNumberValues(OptionName, OptionValueString, OptionValueCount, OptionValueDelimiter, OptionValueType, CmpOpValueMap):
    """Validate numerical option values using option value string, delimiter, value type,
    and a specified map containing comparison operator and value pairs.
    
    Arguments:
        OptionName (str): Command line option name.
        OptionValueString (str): Command line option value.
        OptionValueCount (int): Number of values in OptionValueString.
        OptionValueDelimiter (str): Delimiter used for values in OptionValueString.
        OptionValueType (str): Valid number types (integer or float)
        CmpOpValueMap (dictionary): Comparison operator key and value pairs to
            validate values specified in OptionValueString.

    Notes:
        The function exits with an error message for invalid option values specified
        in OptionValueString

    Examples:

        ValidateOptionNumberValues("-m, --molImageSize",
            Options["--molImageSize"], 2, ",", "integer", {">": 0})

    """
    if not CheckTextValue(OptionValueType, "integer float"):
        PrintError("The option value type specified, %s, for function MiscUtil.ValidateOptionNumberValues  is not valid. Supported value: integer float " % (OptionValueType))
        
    Values = OptionValueString.split(OptionValueDelimiter)
    if len(Values) != OptionValueCount:
        PrintError("The value specified, %s, for option \"%s\" is not valid. It must contain %d %s values separated by \"%s\"" % (OptionValueString, OptionName, OptionValueCount, OptionValueType, OptionValueDelimiter))

    IsIntergerValue = True
    if re.match("^float$", OptionValueType, re.I):
        IsIntergerValue = False
    
    for Value in Values:
        if IsIntergerValue:
            if not IsInteger(Value):
                PrintError("The value specified, %s, for option \"%s\" in string \"%s\" must be an integer." % (Value, OptionName, OptionValueString))
            Value = int(Value)
        else:
            if not IsFloat(Value):
                PrintError("The value specified, %s, for option \"%s\" in string \"%s\" must be a float." % (Value, OptionName, OptionValueString))
            Value = float(Value)
        ValidateOptionNumberValue(OptionName, Value, CmpOpValueMap)
    
def ValidateOptionTextValue(OptionName, OptionValue, ValidValues):
    """Validate option value based on the valid specified values separated by spaces.
    
    Arguments:
        OptionName (str): Command line option name.
        OptionValue (str): Command line option value.
        ValidValues (str): Space delimited string containing valid values.

    Notes:
        The function exits with an error message for an invalid option value.

    """
    
    if not CheckTextValue(OptionValue, ValidValues):
        PrintError("The value specified, %s, for option \"%s\" is not valid. Supported value(s): %s " % (OptionValue, OptionName, ValidValues))

def ValidateOptionsOutputFileOverwrite(OptionName, FilePath, OverwriteOptionName, OverwriteStatus):
    """Validate overwriting of output file.
    
    Arguments:
        OptionName (str): Command line option name.
        FilePath (str): Name of a file with complete file path.
        OverwriteOptionName (str): Overwrite command line option name.
        OverwriteStatus (bool): True, overwrite

    Notes:
        The function exits with an error message for a file that is present and is not allowed
        to be written as indicated by value of OverwriteStatus.

    """
    
    if os.path.exists(FilePath):
        if not OverwriteStatus:
            if len(OverwriteOptionName) > 4:
                ShortOverwriteOptionName = OverwriteOptionName[:4]
            else:
                ShortOverwriteOptionName = OverwriteOptionName
            
            PrintError("The file specified, %s, for option \"%s\" already exist. Use option \"%s\" or \"%s\" and try again.\n" % (FilePath, OptionName, ShortOverwriteOptionName, OverwriteOptionName))

def ValidateOptionsDistinctFileNames(OptionName1, FilePath1, OptionName2, FilePath2):
    """Validate two distinct file names.

    Arguments:
        OptionName1 (str): Command line option name.
        FilePath1 (str): Name of a file with complete file path.
        OptionName2 (str): Command line option name.
        FilePath2 (str): Name of a file with complete file path.

    Notes:
        The function exits with an error message for two non distinct file names.
    
    """
    
    FilePath1Pattern = r"^" + re.escape(FilePath1) + r"$";
    if re.match(FilePath1Pattern, FilePath2, re.I):
        PrintError("The file name specified, %s, for options \"%s\" and \"%s\" must be different.\n" % (FilePath1, OptionName1, OptionName2))

def ProcessOptionInfileParameters(ParamsOptionName, ParamsOptionValue, InfileName = None, OutfileName = None):
    """Process parameters for reading input files and return a map containing
    processed parameter names and values.
    
    Arguments:
        ParamsOptionName (str): Command line input parameters option name.
        ParamsOptionValues (str): Comma delimited list of parameter name and value pairs.
        InfileName (str): Name of input file.
        OutfileName (str): Name of output file.

    Returns:
        dictionary: Processed parameter name and value pairs.

    Notes:
        The parameter name and values specified in ParamsOptionValues are validated before
        returning them in a dictionary.

    """
    
    ParamsInfo = {'RemoveHydrogens': True, 'Sanitize' : True, 'StrictParsing': True,  'SMILESColumn': 1, 'SMILESNameColumn' : 2, 'SMILESDelimiter': ' ', 'SMILESTitleLine': 'auto'}
    _ProcessInfileAndOutfileParameters('Infile', ParamsInfo, ParamsOptionName, ParamsOptionValue, InfileName, OutfileName,)
    
    return ParamsInfo

def ProcessOptionOutfileParameters(ParamsOptionName, ParamsOptionValue, InfileName = None, OutfileName = None):
    """Process parameters for writing output files and return a map containing
    processed parameter names and values.
    
    Arguments:
        ParamsOptionName (str): Command line input parameters option name.
        ParamsOptionValues (str): Comma delimited list of parameter name and value pairs.
        InfileName (str): Name of input file.
        OutfileName (str): Name of output file.

    Returns:
        dictionary: Processed parameter name and value pairs.

    Notes:
        The parameter name and values specified in ParamsOptionValues are validated before
        returning them in a dictionary.

        The default value of some parameters may depend on type of input file. Consequently,
        the input file name is also needed.

    """
    ParamsInfo = {'Compute2DCoords' : 'auto', 'Kekulize': False, 'SMILESDelimiter': ' ', 'SMILESIsomeric': True,  'SMILESTitleLine': True}
    _ProcessInfileAndOutfileParameters('Outfile', ParamsInfo, ParamsOptionName, ParamsOptionValue, InfileName, OutfileName)
    
    return ParamsInfo
    
def _ProcessInfileAndOutfileParameters(Mode, ParamsInfo, ParamsOptionName, ParamsOptionValue, InfileName, OutfileName):
    """Process specified infile and outfile paramaters.
    
    """
    if re.match("^auto$", ParamsOptionValue, re.I):
        # No specific parameters to process except for parameters with possible auto value...
        _ProcessInfileAndOutfileAutoParameters(Mode, ParamsInfo, ParamsOptionName, ParamsOptionValue, InfileName, OutfileName)
        return
    
    ParamsOptionValue = re.sub(" ", "", ParamsOptionValue)
    if not ParamsOptionValue:
        PrintError("No valid parameter name and value pairs specified using \"%s\" option" % ParamsOptionName)
    
    ParamsOptionValueWords = ParamsOptionValue.split(",")
    if len(ParamsOptionValueWords) % 2:
        PrintError("The number of comma delimited paramater names and values, %d, specified using \"%s\" option must be an even number." % (len(ParamsOptionValueWords), ParamsOptionName))
        
    # Setup a canonical paramater names...
    ValidParamNames = []
    CanonicalParamNamesMap = {}
    for ParamName in sorted(ParamsInfo):
        ValidParamNames.append(ParamName)
        CanonicalParamNamesMap[ParamName.lower()] = ParamName
    
    # Validate paramater name and value pairs...
    for Index in range(0, len(ParamsOptionValueWords), 2):
        Name = ParamsOptionValueWords[Index]
        Value = ParamsOptionValueWords[Index + 1]

        CanonicalName = Name.lower()
        if  not CanonicalName in CanonicalParamNamesMap:
            PrintError("The parameter name, %s, specified using \"%s\" is not a valid name. Supported parameter names: %s" % (Name, ParamsOptionName, " ".join(ValidParamNames)))

        ParamName = CanonicalParamNamesMap[CanonicalName]
        ParamValue = Value
        
        if re.match("^(Sanitize|StrictParsing|RemoveHydrogens|Kekulize|SMILESIsomeric)$", ParamName, re.I):
            if not re.match("^(Yes|No|True|False)$", Value, re.I):
                PrintError("The parameter value, %s, specified for parameter name, %s, using \"%s\" option is not a valid value. Supported values: Yes No True False" % (Value, Name, ParamsOptionName))
            ParamValue = True
            if re.match("^(No|False)$", Value, re.I):
                ParamValue = False
        elif re.match("^SMILESTitleLine$", ParamName, re.I):
            if re.match("^Infile$", Mode, re.I):
                if not re.match("^(Yes|No|True|False|Auto)$", Value, re.I):
                    PrintError("The parameter value, %s, specified for paramater name, %s, using \"%s\" option is not a valid value. Supported values: Yes No True False Auto" % (Value, Name, ParamsOptionName))
            elif re.match("^Outfile$", Mode, re.I):
                if not re.match("^(Yes|No|True|False)$", Value, re.I):
                    PrintError("The parameter value, %s, specified for parameter name, %s, using \"%s\" option is not a valid value. Supported values: Yes No True False" % (Value, Name, ParamsOptionName))
                ParamValue = True
                if re.match("^(No|False)$", Value, re.I):
                    ParamValue = False
        elif re.match("^SMILESDelimiter$", ParamName, re.I):
            if not re.match("^(space|tab|comma)$", Value, re.I):
                PrintError("The parameter value, %s, specified for parameter name, %s, using \"%s\" option is not a valid value. Supported values: space tab comma" % (Value, Name, ParamsOptionName))
            ParamValue = " "
            if re.match("^tab$", Value, re.I):
                ParamValue = "\t"
            elif re.match("^comma$", Value, re.I):
                ParamValue = ","
        elif re.match("^Compute2DCoords$", ParamName, re.I):
            # No need to set the value. It would be processed later to handle "auto" value...
            if not re.match("^(Yes|No|True|False|Auto)$", Value, re.I):
                PrintError("The parameter value, %s, specified for paramater name, %s, using \"%s\" option is not a valid value. Supported values: Yes No True False Auto" % (Value, Name, ParamsOptionName))
        else:
            ParamValue = int(Value)
            if ParamValue <= 0:
                PrintError("The parameter value, %s, specified for parameter name, %s, using \"%s\" option is not a valid value. Supported values: > 0" % (Value, Name, ParamsOptionName))
        
        # Set value...
        ParamsInfo[ParamName] = ParamValue
        
    # Handle paramaters with possible auto values...
    _ProcessInfileAndOutfileAutoParameters(Mode, ParamsInfo, ParamsOptionName, ParamsOptionValue, InfileName, OutfileName)

def _ProcessInfileAndOutfileAutoParameters(Mode, ParamsInfo, ParamsOptionName, ParamsOptionValue, InfileName, OutfileName):
    """Process parameters with possible auto values.
    
    """
    if re.match("^Infile$", Mode, re.I):
        # SMILESTitleLine parameter...
        Value = ParamsInfo["SMILESTitleLine"]
        ParamValue = False
        if re.match("^auto$", Value, re.I):
            if InfileName is not None:
                if CheckFileExt(InfileName, "smi csv tsv txt"):
                    ParamValue = DoesSMILESFileContainTitleLine(InfileName)
        elif re.match("^(Yes|True)$", Value, re.I):
            ParamValue = True
        ParamsInfo["SMILESTitleLine"] = ParamValue
    elif re.match("^Outfile$", Mode, re.I):
        # Compute2DCoords parameter...
        Value = ParamsInfo["Compute2DCoords"]
        ParamValue = False
        if re.match("^auto$", Value, re.I):
            if InfileName is not None:
                if CheckFileExt(InfileName, "smi csv tsv txt"):
                    ParamValue = True
            if OutfileName is not None:
                if CheckFileExt(OutfileName, "smi csv tsv txt"):
                    # No need to compute 2D coords for SMILES file...
                    ParamValue = False
        elif re.match("^(Yes|True)$", Value, re.I):
            ParamValue = True
        ParamsInfo["Compute2DCoords"] = ParamValue

def ReplaceHTMLEntitiesInText(Text):
    """Check and replace the followng HTML entities to their respective code
    for display in a browser: < (less than), > (greater than), & (ampersand),
    " (double quote),  and ' (single quote).

    Arguments:
        Text (str): Text value.

    Returns:
        str : Modifed text value.

    """

    if re.search("""(<|>|&|"|')""", Text):
        return Text.replace("<", "&lt;").replace(">", "&gt;").replace("&", "&amp;").replace('"', "&quot;").replace("'","&apos;")
    else:
        return Text

def TruncateText(Text, Width, TrailingChars = "..."):
    """Truncate text using specified width along with appending any trailing
    characters.
    
    Arguments:
        Text (string): Input text.
        Width (int): Max number of characters before truncating text.
        Delimiter (string): Trailing characters to append or None.

    Returns:
        str : Truncated text

    """
    
    if len(Text) < Width:
        return Text
    
    TruncatedText = (Text[:Width] + TrailingChars) if not IsEmpty(TrailingChars) else Text[:Width] 
    
    return TruncatedText

def WrapText(Text, Delimiter, Width):
    """Wrap text using specified delimiter and width.
    
    Arguments:
        Text (string): Input text
        Delimiter (string): Delimiter for wrapping text
        Width (int): Max number of characters before wrapping text

    Returns:
        str : Wrapped text

    """
    WrappedText = Delimiter.join(textwrap.wrap(Text, width = Width))
    
    return WrappedText
