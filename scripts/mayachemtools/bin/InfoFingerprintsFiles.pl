#!/usr/bin/perl -w
#
# File: InfoFingerprintsFiles.pl
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
use Fingerprints::FingerprintsFileUtil;
use Fingerprints::FingerprintsStringUtil;

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

my(@FingerprintsFilesList);
@FingerprintsFilesList = ExpandFileNames(\@ARGV, "sdf sd fpf fp csv tsv");

# Process options...
print "Processing options...\n";
my(%OptionsInfo);
ProcessOptions();

# Setup information about input files...
print "Checking input fingerprints file(s)...\n";
my(%FingerprintsFilesInfo);
RetrieveFingerprintsFilesInfo();

# Process input files..
my($FileIndex);
if (@FingerprintsFilesList > 1) {
  print "\nProcessing fingerprints files...\n";
}
for $FileIndex (0 .. $#FingerprintsFilesList) {
  if ($FingerprintsFilesInfo{FileOkay}[$FileIndex]) {
    print "\nProcessing file $FingerprintsFilesList[$FileIndex]...\n";
    ListFingerprintsFileInfo($FileIndex);
  }
}
ListTotalSizeOfFiles();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# List approptiate information...
#
sub ListFingerprintsFileInfo {
  my($FileIndex) = @_;
  my($FileName, $FingerprintsFileIO, $InvalidFingerprintsFileData, $InvalidFingerprintsData, $DataEntryCount, $ValidDataEntryCount, $InvalidDataEntryCount, $MissingDataEntryCount, $BitVectorDataEntryCount, $VectorDataEntryCount, $FingerprintsObject, $FingerprintsType, $TotalBitDensity, $FileType, $DataEntryLabel);

  $FileType = $FingerprintsFilesInfo{FileType}[$FileIndex];
  $DataEntryLabel = ($FileType =~ /^SD$/i) ? 'compounds' : 'lines';

  ($DataEntryCount, $ValidDataEntryCount, $InvalidDataEntryCount, $MissingDataEntryCount, $BitVectorDataEntryCount, $VectorDataEntryCount, $TotalBitDensity) = (0) x 7;

  $FingerprintsFileIO = Fingerprints::FingerprintsFileUtil::NewFingerprintsFileIO(%{$FingerprintsFilesInfo{FingerprintsFileIOParameters}[$FileIndex]});
  $FingerprintsFileIO->Open();

  $InvalidFingerprintsFileData = $FingerprintsFileIO->IsFingerprintsFileDataValid() ? 0 : 1;

  FINGERPRINTS: while ($FingerprintsFileIO->Read()) {
    $DataEntryCount++;

    # Missing data...
    if ($InvalidFingerprintsFileData) {
      $MissingDataEntryCount++;
      if ($OptionsInfo{ValidateData} || $OptionsInfo{CountEmptyFingerprints}) {
	ListEmptyOrInvalidFingerprintsDataInfo('EmptyData', $FingerprintsFileIO, $FileType);
      }
      next FINGERPRINTS;
    }
    $InvalidFingerprintsData = $FingerprintsFileIO->IsFingerprintsDataValid() ? 0 : 1;

    # Invalid data...
    if ($InvalidFingerprintsData) {
      $InvalidDataEntryCount++;
      if ($OptionsInfo{ValidateData}) {
	ListEmptyOrInvalidFingerprintsDataInfo('InvalidData', $FingerprintsFileIO, $FileType);
      }
      next FINGERPRINTS;
    }
    $ValidDataEntryCount++;

    $FingerprintsObject = $FingerprintsFileIO->GetFingerprints();
    $FingerprintsType = $FingerprintsObject->GetVectorType();

    if ($FingerprintsType =~ /^FingerprintsBitVector$/i) {
      $BitVectorDataEntryCount++;
      if ($OptionsInfo{ListAverageBitDensity}) {
	$TotalBitDensity += $FingerprintsObject->GetFingerprintsBitDensity();
      }
    }
    elsif ($FingerprintsType =~ /^FingerprintsVector$/i) {
      $VectorDataEntryCount++;
    }

    if ($OptionsInfo{ListFingerprintsDataEntryInfo}) {
      ListFingerprintsDataEntryInfo($FingerprintsFileIO, $FileType);
    }

  }
  $FingerprintsFileIO->Close();

  print "\nFingerprints file type: $FileType\n";
  if ($FileType =~ /^SD$/i) {
    print "Number of compounds: $DataEntryCount\n";
  }
  else {
    print "Number of data lines: $DataEntryCount\n";
  }

  ListFileTypeHeaderInfo($FingerprintsFileIO, $FileType);

  print "\nNumber of $DataEntryLabel with valid fingerprints string data: $ValidDataEntryCount\n";
  print "Number of $DataEntryLabel with bit-vector fingerprints string data: $BitVectorDataEntryCount\n";
  print "Number of $DataEntryLabel with vector fingerprints string data: $VectorDataEntryCount\n";

  if ($OptionsInfo{CountEmptyFingerprints}) {
    print "Number of $DataEntryLabel with missing fingerprints data: $MissingDataEntryCount\n";
    print "Number of $DataEntryLabel with invalid fingerprints data: $InvalidDataEntryCount\n";
  }

  if ($OptionsInfo{ListAverageBitDensity} && $BitVectorDataEntryCount) {
    my($AverageBitDensity);
    $AverageBitDensity = $TotalBitDensity/$BitVectorDataEntryCount;
    $AverageBitDensity = sprintf("%.2f", $AverageBitDensity) + 0;
    print "\nAverage bit density: $AverageBitDensity\n";
  }


  # File size and modification information...
  print "\nFile size: ", FormatFileSize($FingerprintsFilesInfo{FileSize}[$FileIndex]), " \n";
  print "Last modified: ", $FingerprintsFilesInfo{FileLastModified}[$FileIndex], " \n";
}

# List empty or invalid fingerprints file data information...
#
sub ListEmptyOrInvalidFingerprintsDataInfo {
  my($Mode, $FingerprintsFileIO, $FileType) = @_;
  my($ModeInfo);

  $ModeInfo = ($Mode =~ /^EmptyData$/i) ? "no" : "invalid";

  if ($FileType =~ /^SD$/i) {
    my($CmpdNum, $CmpdString);

    $CmpdNum = $FingerprintsFileIO->GetCompoundNum();
    if ($OptionsInfo{DetailLevel} >= 3 ) {
      $CmpdString = $FingerprintsFileIO->GetCompoundString();
      print "Compound number $CmpdNum contains $ModeInfo fingerprints data: $CmpdString \n";
    }
    elsif ($OptionsInfo{DetailLevel} >= 1 ) {
      print "Compound number $CmpdNum contains $ModeInfo fingerprints data...\n";
    }
  }
  else {
    my($LineNum, $DataLine);

    $LineNum = $FingerprintsFileIO->GetLineNum();
    if ($OptionsInfo{DetailLevel} >= 3 ) {
      $DataLine = $FingerprintsFileIO->GetDataLine();
      print "Data line number $LineNum contains $ModeInfo fingerprints data: $DataLine \n";
    }
    elsif ($OptionsInfo{DetailLevel} >= 1 ) {
      print "Data line number $LineNum contains $ModeInfo fingerprints data...\n";
    }
  }
}

# List detailed information about fingerprints data entry...
#
sub ListFingerprintsDataEntryInfo {
  my($FingerprintsFileIO, $FileType) = @_;
  my($FingerprintsObject, $FingerprintsString, $FingerprintsType, $FingerprintsDescription, $FingerprintsSize, $FingerprintsBitStringFormat, $FingerprintsBitOrder, $BitDensity, $NumOfOnBits, $FingerprintsVectorValuesType, $FingerprintsVectorValuesFormat, $NumOfNonZeroValues);

  $FingerprintsObject = $FingerprintsFileIO->GetFingerprints();
  $FingerprintsString = $FingerprintsFileIO->GetFingerprintsString();

  $FingerprintsType = $FingerprintsObject->GetVectorType();

  if ($FingerprintsType =~ /^FingerprintsBitVector$/i) {
    $BitDensity = '';
    $NumOfOnBits = '';

    ($FingerprintsType, $FingerprintsDescription, $FingerprintsSize, $FingerprintsBitStringFormat, $FingerprintsBitOrder) = Fingerprints::FingerprintsStringUtil::GetFingerprintsStringValues($FingerprintsString);

    if ($OptionsInfo{ListBitDensity} || $OptionsInfo{ListNumOfOnBits}) {
      if ($OptionsInfo{ListBitDensity}) {
	$BitDensity = $FingerprintsObject->GetFingerprintsBitDensity();
      }
      if ($OptionsInfo{ListNumOfOnBits}) {
	$NumOfOnBits = $FingerprintsObject->GetNumOfSetBits();
      }
    }
  }
  elsif ($FingerprintsType =~ /^FingerprintsVector$/i) {
    $NumOfNonZeroValues = '';

    ($FingerprintsType, $FingerprintsDescription, $FingerprintsSize, $FingerprintsVectorValuesType, $FingerprintsVectorValuesFormat) = Fingerprints::FingerprintsStringUtil::GetFingerprintsStringValues($FingerprintsString);

    if ($OptionsInfo{ListNumOfNonZeroValues}) {
      if ($FingerprintsVectorValuesType =~ /^AlphaNumericalValues$/i) {
	$NumOfNonZeroValues = 'NA';
      }
      else {
	$NumOfNonZeroValues = $FingerprintsObject->GetNumOfNonZeroValues();
      }
    }
  }

  if ($FileType =~ /^SD$/i) {
    print "Compound number: " . $FingerprintsFileIO->GetCompoundNum();
  }
  else {
    print "Data line number: " . $FingerprintsFileIO->GetLineNum();
  }

  if ($OptionsInfo{ListFingerprintsType}) {
    print "; FPType: $FingerprintsType";
  }
  if ($OptionsInfo{ListFingerprintsDescription}) {
    print "; FPDescription: $FingerprintsDescription";
  }
  if ($OptionsInfo{ListFingerprintsSize}) {
    print "; FPSize: $FingerprintsSize";
  }

  if ($FingerprintsType =~ /^FingerprintsBitVector$/i) {
    if ($OptionsInfo{ListFingerprintsBitStringFormat}) {
      print "; FPBitStringFormat: $FingerprintsBitStringFormat";
    }
    if ($OptionsInfo{ListFingerprintsBitOrder}) {
      print "; FPBitOrder: $FingerprintsBitOrder";
    }
    if ($OptionsInfo{ListBitDensity}) {
      print "; BitDensity: $BitDensity";
    }
    if ($OptionsInfo{ListNumOfOnBits}) {
      print "; NumOfOnBits: $NumOfOnBits";
    }
  }
  elsif ($FingerprintsType =~ /^FingerprintsVector$/i) {
    if ($OptionsInfo{ListFingerprintsVectorValuesType}) {
      print "; FPVectorValuesType: $FingerprintsVectorValuesType";
    }
    if ($OptionsInfo{ListFingerprintsVectorValuesFormat}) {
      print "; FPVectorValuesFormat: $FingerprintsVectorValuesFormat";
    }
    if ($OptionsInfo{ListNumOfNonZeroValues}) {
      print "; NumOfNonZeroValues: $NumOfNonZeroValues";
    }
  }
  print "\n";
}

# List file type header information...
#
sub ListFileTypeHeaderInfo {
  my($FingerprintsFileIO, $FileType) = @_;
  my($Key, $Value, @DataColLabels, %HeaderDataKeysAndValues);

  if ($FileType =~ /^Text$/i) {
    @DataColLabels = $FingerprintsFileIO->GetDataColLabels();
    print "Number of columns: " . scalar @DataColLabels . "\n";
    print "Column labels: ", JoinWords(\@DataColLabels, ", ", 1), "\n";
  }
  elsif ($FileType =~ /^FP$/i) {
    %HeaderDataKeysAndValues = $FingerprintsFileIO->GetHeaderDataKeysAndValues();

    print "\nFP file header data keys and values: \n#\n";
    for $Key ($FingerprintsFileIO->GetHeaderDataKeys()) {
      $Value = $HeaderDataKeysAndValues{$Key};
      print "# $Key = $Value\n";
    }
    print "#\n";
  }
}

# Total size of all the fiels...
sub ListTotalSizeOfFiles {
  my($FileOkayCount, $TotalSize, $Index);

  $FileOkayCount = 0;
  $TotalSize = 0;

  for $Index (0 .. $#FingerprintsFilesList) {
    if ($FingerprintsFilesList[$Index]) {
      $FileOkayCount++;
      $TotalSize += $FingerprintsFilesInfo{FileSize}[$Index];
    }
  }
  if ($FileOkayCount > 1) {
    print "\nTotal size of $FileOkayCount files: ", FormatFileSize($TotalSize), "\n";
  }
}

# Retrieve information about fingerprints files...
#
sub RetrieveFingerprintsFilesInfo {
  my($FingerprintsFile, $Index, $FileDir, $FileExt, $FileName, $FileType, $InDelim, $ModifiedTimeString, $ModifiedDateString, %FingerprintsFileIOParameters);

  %FingerprintsFilesInfo = ();
  @{$FingerprintsFilesInfo{FileOkay}} = ();
  @{$FingerprintsFilesInfo{FileType}} = ();
  @{$FingerprintsFilesInfo{FileSize}} = ();
  @{$FingerprintsFilesInfo{FileLastModified}} = ();
  @{$FingerprintsFilesInfo{InDelim}} = ();

  @{$FingerprintsFilesInfo{FingerprintsFileIOParameters}} = ();

  FILELIST: for $Index (0 .. $#FingerprintsFilesList) {
    $FingerprintsFile = $FingerprintsFilesList[$Index];

    $FingerprintsFilesInfo{FileOkay}[$Index] = 0;
    $FingerprintsFilesInfo{FileType}[$Index] = '';
    $FingerprintsFilesInfo{FileSize}[$Index] = 0;
    $FingerprintsFilesInfo{FileLastModified}[$Index] = '';
    $FingerprintsFilesInfo{InDelim}[$Index] = "";

    %{$FingerprintsFilesInfo{FingerprintsFileIOParameters}[$Index]} = ();

    $FingerprintsFile = $FingerprintsFilesList[$Index];
    if (!(-e $FingerprintsFile)) {
      warn "Warning: Ignoring file $FingerprintsFile: It doesn't exist\n";
      next FILELIST;
    }

    $FileType = Fingerprints::FingerprintsFileUtil::GetFingerprintsFileType($FingerprintsFile);
    if (IsEmpty($FileType)) {
      warn "Warning: Ignoring file $FingerprintsFile: It's not a fingerprints file\n";
      next FILELIST;
    }

    $FileDir = ""; $FileName = ""; $FileExt = "";
    ($FileDir, $FileName, $FileExt) = ParseFileName($FingerprintsFile);

    $InDelim = ($FileExt =~ /^tsv$/i) ? 'Tab' : $OptionsInfo{InDelim};

    # Setup FingerprintsFileIO parameters...
    %FingerprintsFileIOParameters = ();
    FILEIOPARAMETERS: {
      if ($FileType =~ /^SD$/i) {
	%FingerprintsFileIOParameters = ('Name' => $FingerprintsFile, 'Mode' => 'Read', 'ValidateData' => $OptionsInfo{ValidateData}, 'DetailLevel' => 1, 'FingerprintsFieldLabel' => $OptionsInfo{FingerprintsFieldLabel});
	last FILEIOPARAMETERS;
      }
      if ($FileType =~ /^FP$/i) {
	%FingerprintsFileIOParameters = ('Name' => $FingerprintsFile, 'Mode' => 'Read', 'ValidateData' => $OptionsInfo{ValidateData}, 'DetailLevel' => 1);
	last FILEIOPARAMETERS;
      }
      if ($FileType =~ /^Text$/i) {
	%FingerprintsFileIOParameters = ('Name' => $FingerprintsFile, 'Mode' => 'Read', 'ValidateData' => $OptionsInfo{ValidateData}, 'DetailLevel' => 1, 'FingerprintsCol' => $OptionsInfo{FingerprintsCol}, 'ColMode' => $OptionsInfo{ColMode}, 'InDelim' => $OptionsInfo{InDelim});
	last FILEIOPARAMETERS;
      }
      warn "Warning: File type for fingerprints file, $FingerprintsFile, is not valid. Supported file types: SD, FP or Text\n";
      next FILELIST;
    }

    $FingerprintsFilesInfo{FileOkay}[$Index] = 1;
    $FingerprintsFilesInfo{FileType}[$Index] = $FileType;

    $FingerprintsFilesInfo{FileSize}[$Index] = FileSize($FingerprintsFile);
    ($ModifiedTimeString, $ModifiedDateString) = FormattedFileModificationTimeAndDate($FingerprintsFile);
    $FingerprintsFilesInfo{FileLastModified}[$Index] = "$ModifiedTimeString; $ModifiedDateString";

    $FingerprintsFilesInfo{InDelim}[$Index] = $InDelim;

    %{$FingerprintsFilesInfo{FingerprintsFileIOParameters}[$Index]} = %FingerprintsFileIOParameters;
  }
}

# Process option values...
sub ProcessOptions {
  %OptionsInfo = ();

  $OptionsInfo{ListAverageBitDensity} = ($Options{all} || $Options{averagebitdensity}) ? 1 :0;
  $OptionsInfo{ListBitDensity} = ($Options{all} || $Options{bitdensity}) ? 1 :0;

  if ($OptionsInfo{ListAverageBitDensity}) {
    # List bit density as well...
    $OptionsInfo{ListBitDensity} = 1;
  }

  # By default, count number of rows containing fingerprints data...
  $OptionsInfo{CountFingerprints} = 1;
  $OptionsInfo{CountEmptyFingerprints} = ($Options{all} || $Options{empty}) ? 1 :0;

  $OptionsInfo{ColMode} = $Options{colmode};
  if (IsNotEmpty($Options{fingerprintscol})) {
    if ($Options{colmode} =~ /^ColNum$/i) {
      if (!IsPositiveInteger($Options{fingerprintscol})) {
	die "Error: Column value, $Options{fingerprintscol}, specified using \"--FingerprintsCol\" is not valid: Allowed integer values: > 0.\n";
      }
    }
    $OptionsInfo{FingerprintsCol} = $Options{fingerprintscol};
  }
  else {
    $OptionsInfo{FingerprintsCol} = 'AutoDetect';
  }

  if (IsNotEmpty($Options{fingerprintsfield})) {
    $OptionsInfo{FingerprintsFieldLabel} = $Options{fingerprintsfield};
  }
  else {
    $OptionsInfo{FingerprintsFieldLabel} = 'AutoDetect';
  }

  $OptionsInfo{ValidateData} = ($Options{all} || $Options{datacheck}) ? 1 :0;
  $OptionsInfo{DetailLevel} = $Options{detail};

  $OptionsInfo{ListFingerprintsType} = ($Options{all} || $Options{fingerprintstype}) ? 1 :0;
  $OptionsInfo{ListFingerprintsDescription} = ($Options{all} || $Options{fingerprintsdescription}) ? 1 :0;
  $OptionsInfo{ListFingerprintsSize} = ($Options{all} || $Options{fingerprintssize}) ? 1 :0;

  $OptionsInfo{ListFingerprintsBitStringFormat} = ($Options{all} || $Options{fingerprintsbitstringformat}) ? 1 :0;
  $OptionsInfo{ListFingerprintsBitOrder} = ($Options{all} || $Options{fingerprintsbitorder}) ? 1 :0;

  $OptionsInfo{ListFingerprintsVectorValuesType} = ($Options{all} || $Options{fingerprintsvectorvaluestype}) ? 1 :0;
  $OptionsInfo{ListFingerprintsVectorValuesFormat} = ($Options{all} || $Options{fingerprintsvectorvaluesformat}) ? 1 :0;

  $OptionsInfo{InDelim} = $Options{indelim};

  $OptionsInfo{ListNumOfOnBits} = ($Options{all} || $Options{numofonbits}) ? 1 :0;
  $OptionsInfo{ListNumOfNonZeroValues} = ($Options{all} || $Options{numofnonzerovalues}) ? 1 :0;

  $OptionsInfo{ListFingerprintsDataEntryInfo} = ($OptionsInfo{ListFingerprintsType} || $OptionsInfo{ListFingerprintsDescription} || $OptionsInfo{ListFingerprintsSize} || $OptionsInfo{ListFingerprintsBitStringFormat} || $OptionsInfo{ListFingerprintsBitOrder} || $OptionsInfo{ListFingerprintsVectorValuesType} || $OptionsInfo{ListFingerprintsVectorValuesFormat} || $OptionsInfo{ListBitDensity} || $OptionsInfo{ListAverageBitDensity} || $OptionsInfo{ListNumOfOnBits} ||  $OptionsInfo{ListNumOfNonZeroValues}) ? 1 : 0;

}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{colmode} = 'colnum';
  $Options{detail} = 1;
  $Options{indelim} = 'comma';

  if (!GetOptions(\%Options, "all|a", "averagebitdensity", "bitdensity", "count", "colmode|c=s", "detail|d=i", "datacheck", "empty|e", "fingerprintsfield=s", "fingerprintscol=s", "fingerprintstype", "fingerprintsdescription", "fingerprintssize", "fingerprintsbitstringformat", "fingerprintsbitorder", "fingerprintsvectorvaluestype", "fingerprintsvectorvaluesformat", "help|h", "indelim=s", "numofonbits", "numofnonzerovalues", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{colmode} !~ /^(ColNum|ColLabel)$/i) {
    die "Error: The value specified, $Options{colmode}, for option \"-c, --ColMode\" is not valid. Allowed values: ColNum, or ColLabel\n";
  }
  if (!IsPositiveInteger($Options{detail})) {
    die "Error: The value specified, $Options{detail}, for option \"-d, --detail\" is not valid. Allowed values: > 0 \n";
  }
  if ($Options{indelim} !~ /^(comma|semicolon)$/i) {
    die "Error: The value specified, $Options{indelim}, for option \"--InDelim\" is not valid. Allowed values: comma, or semicolon\n";
  }
}

__END__

=head1 NAME

InfoFingerprintsFiles.pl - List information about fingerprints data in SD, FP and CSV/TSV text file(s)

=head1 SYNOPSIS

InfoFingerprintsFiles.pl SDFile(s) FPFile(s) TextFile(s)...

InfoFingerprintsFiles.pl [B<-a, --all>] [B<--AverageBitDensity>] [B<--BitDensity>]
[B<-c, --count>] [B<-c, --ColMode> I<ColNum | ColLabel>] [B<--DataCheck>]
[B<-d, --detail> I<InfoLevel>] [B<-e, --empty>] [B<--FingerprintsCol> I<col number | col name>]
[B<--FingerprintsField> I<FieldLabel>] [B<--FingerprintsType>] [B<--FingerprintsDescription>]
[B<--FingerprintsSize>] [B<--FingerprintsBitStringFormat>] [B<--FingerprintsBitOrder>]
[B<--FingerprintsVectorValuesType>] [B<--FingerprintsVectorValuesFormat>]
[B<-h, --help>] [B<--InDelim> I<comma | semicolon>]
[B<--NumOfOnBits>] [B<--NumOfNonZeroValues>]
[B<-w, --WorkingDir> dirname] SDFile(s) FPFile(s) TextFile(s)...

=head1 DESCRIPTION

List information about fingerprints data in I<SD, FP and CSV/TSV> text file(s): number of
rows containing fingerprints data, type of fingerprints vector, description and size of fingerprints,
bit density and average bit density for bit-vector fingerprints strings, and so on.

The scripts InfoFingerprintsSDFiles.pl and InfoFingerprintsTextFiles.pl have been removed from the
current release of MayaChemTools and their functionality merged with this script.

The valid I<SDFile> extensions are I<.sdf> and I<.sd>. All SD files in a current directory
can be specified either by I<*.sdf> or the current directory name.

The valid I<FPFile> extensions are I<.fpf> and I<.fp>. All FP files in a current directory
can be specified either by I<*.fpf> or the current directory name.

The valid I<TextFile> extensions are I<.csv> and I<.tsv> for comma/semicolon and tab
delimited text files respectively. All other file names are ignored. All text files in a
current directory can be specified by I<*.csv>, I<*.tsv>, or the current directory
name. The B<--indelim> option determines the format of I<TextFile(s)>. Any file
which doesn't correspond to the format indicated by B<--indelim> option is ignored.

Format of fingerprint strings data in I<SDFile(s), FPFile(s) and TextFile(s)> is automatically
detected.

Example of I<FP> file containing fingerprints bit-vector string data:

    #
    # Package = MayaChemTools 7.4
    # ReleaseDate = Oct 21, 2010
    #
    # TimeStamp =  Mon Mar 7 15:14:01 2011
    #
    # FingerprintsStringType = FingerprintsBitVector
    #
    # Description = PathLengthBits:AtomicInvariantsAtomTypes:MinLength1:...
    # Size = 1024
    # BitStringFormat = HexadecimalString
    # BitsOrder = Ascending
    #
    Cmpd1 9c8460989ec8a49913991a6603130b0a19e8051c89184414953800cc21510...
    Cmpd2 000000249400840040100042011001001980410c000000001010088001120...
    ... ...
    ... ..

Example of I<FP> file containing fingerprints vector string data:

    #
    # Package = MayaChemTools 7.4
    # ReleaseDate = Oct 21, 2010
    #
    # TimeStamp =  Mon Mar 7 15:14:01 2011
    #
    # FingerprintsStringType = FingerprintsVector
    #
    # Description = PathLengthBits:AtomicInvariantsAtomTypes:MinLength1:...
    # VectorStringFormat = IDsAndValuesString
    # VectorValuesType = NumericalValues
    #
    Cmpd1 338;C F N O C:C C:N C=O CC CF CN CO C:C:C C:C:N C:CC C:CF C:CN C:
    N:C C:NC CC:N CC=O CCC CCN CCO CNC NC=O O=CO C:C:C:C C:C:C:N C:C:CC...;
    33 1 2 5 21 2 2 12 1 3 3 20 2 10 2 2 1 2 2 2 8 2 5 1 1 1 19 2 8 2 2 2 2
    6 2 2 2 2 2 2 2 2 3 2 2 1 4 1 5 1 1 18 6 2 2 1 2 10 2 1 2 1 2 2 2 2 ...
    Cmpd2 103;C N O C=N C=O CC CN CO CC=O CCC CCN CCO CNC N=CN NC=O NCN O=C
    O C CC=O CCCC CCCN CCCO CCNC CNC=N CNC=O CNCN CCCC=O CCCCC CCCCN CC...;
    15 4 4 1 2 13 5 2 2 15 5 3 2 2 1 1 1 2 17 7 6 5 1 1 1 2 15 8 5 7 2 2 2 2
    1 2 1 1 3 15 7 6 8 3 4 4 3 2 2 1 2 3 14 2 4 7 4 4 4 4 1 1 1 2 1 1 1 ...
    ... ...
    ... ...

Example of I<SD> file containing fingerprints bit-vector string data:

    ... ...
    ... ...
    $$$$
    ... ...
    ... ...
    ... ...
    41 44  0  0  0  0  0  0  0  0999 V2000
     -3.3652    1.4499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    ... ...
    2  3  1  0  0  0  0
    ... ...
    M  END
    >  <CmpdID>
    Cmpd1

    >  <PathLengthFingerprints>
    FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes:MinLengt
    h1:MaxLength8;1024;HexadecimalString;Ascending;9c8460989ec8a49913991a66
    03130b0a19e8051c89184414953800cc2151082844a201042800130860308e8204d4028
    00831048940e44281c00060449a5000ac80c894114e006321264401600846c050164462
    08190410805000304a10205b0100e04c0038ba0fad0209c0ca8b1200012268b61c0026a
    aa0660a11014a011d46

    $$$$
    ... ...
    ... ...

Example of CSV I<Text> file containing fingerprints bit-vector string data:

    "CompoundID","PathLengthFingerprints"
    "Cmpd1","FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes
    :MinLength1:MaxLength8;1024;HexadecimalString;Ascending;9c8460989ec8a4
    9913991a6603130b0a19e8051c89184414953800cc2151082844a20104280013086030
    8e8204d402800831048940e44281c00060449a5000ac80c894114e006321264401..."
    ... ...
    ... ...

The current release of MayaChemTools supports the following types of fingerprint
bit-vector and vector strings:

    FingerprintsVector;AtomNeighborhoods:AtomicInvariantsAtomTypes:MinRadi
    us0:MaxRadius2;41;AlphaNumericalValues;ValuesString;NR0-C.X1.BO1.H3-AT
    C1:NR1-C.X3.BO3.H1-ATC1:NR2-C.X1.BO1.H3-ATC1:NR2-C.X3.BO4-ATC1 NR0-C.X
    1.BO1.H3-ATC1:NR1-C.X3.BO3.H1-ATC1:NR2-C.X1.BO1.H3-ATC1:NR2-C.X3.BO4-A
    TC1 NR0-C.X2.BO2.H2-ATC1:NR1-C.X2.BO2.H2-ATC1:NR1-C.X3.BO3.H1-ATC1:NR2
    -C.X2.BO2.H2-ATC1:NR2-N.X3.BO3-ATC1:NR2-O.X1.BO1.H1-ATC1 NR0-C.X2.B...

    FingerprintsVector;AtomTypesCount:AtomicInvariantsAtomTypes:ArbitraryS
    ize;10;NumericalValues;IDsAndValuesString;C.X1.BO1.H3 C.X2.BO2.H2 C.X2
    .BO3.H1 C.X3.BO3.H1 C.X3.BO4 F.X1.BO1 N.X2.BO2.H1 N.X3.BO3 O.X1.BO1.H1
    O.X1.BO2;2 4 14 3 10 1 1 1 3 2

    FingerprintsVector;AtomTypesCount:SLogPAtomTypes:ArbitrarySize;16;Nume
    ricalValues;IDsAndValuesString;C1 C10 C11 C14 C18 C20 C21 C22 C5 CS F
    N11 N4 O10 O2 O9;5 1 1 1 14 4 2 1 2 2 1 1 1 1 3 1

    FingerprintsVector;AtomTypesCount:SLogPAtomTypes:FixedSize;67;OrderedN
    umericalValues;IDsAndValuesString;C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C
    12 C13 C14 C15 C16 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C27 CS N1 N
    2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12 N13 N14 NS O1 O2 O3 O4 O5 O6 O7 O8
    O9 O10 O11 O12 OS F Cl Br I Hal P S1 S2 S3 Me1 Me2;5 0 0 0 2 0 0 0 0 1
    1 0 0 1 0 0 0 14 0 4 2 1 0 0 0 0 0 2 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0...

    FingerprintsVector;EStateIndicies:ArbitrarySize;11;NumericalValues;IDs
    AndValuesString;SaaCH SaasC SaasN SdO SdssC SsCH3 SsF SsOH SssCH2 SssN
    H SsssCH;24.778 4.387 1.993 25.023 -1.435 3.975 14.006 29.759 -0.073 3
    .024 -2.270

    FingerprintsVector;EStateIndicies:FixedSize;87;OrderedNumericalValues;
    ValuesString;0 0 0 0 0 0 0 3.975 0 -0.073 0 0 24.778 -2.270 0 0 -1.435
    4.387 0 0 0 0 0 0 3.024 0 0 0 0 0 0 0 1.993 0 29.759 25.023 0 0 0 0 1
    4.006 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0

    FingerprintsVector;ExtendedConnectivity:AtomicInvariantsAtomTypes:Radi
    us2;60;AlphaNumericalValues;ValuesString;73555770 333564680 352413391
    666191900 1001270906 1371674323 1481469939 1977749791 2006158649 21414
    08799 49532520 64643108 79385615 96062769 273726379 564565671 85514103
    5 906706094 988546669 1018231313 1032696425 1197507444 1331250018 1338
    532734 1455473691 1607485225 1609687129 1631614296 1670251330 17303...

    FingerprintsVector;ExtendedConnectivityCount:AtomicInvariantsAtomTypes
    :Radius2;60;NumericalValues;IDsAndValuesString;73555770 333564680 3524
    13391 666191900 1001270906 1371674323 1481469939 1977749791 2006158649
    2141408799 49532520 64643108 79385615 96062769 273726379 564565671...;
    3 2 1 1 14 1 2 10 4 3 1 1 1 1 2 1 2 1 1 1 2 3 1 1 2 1 3 3 8 2 2 2 6 2
    1 2 1 1 2 1 1 1 2 1 1 2 1 2 1 1 1 1 1 1 1 1 1 2 1 1

    FingerprintsBitVector;ExtendedConnectivityBits:AtomicInvariantsAtomTyp
    es:Radius2;1024;BinaryString;Ascending;0000000000000000000000000000100
    0000000001010000000110000011000000000000100000000000000000000000100001
    1000000110000000000000000000000000010011000000000000000000000000010000
    0000000000000000000000000010000000000000000001000000000000000000000000
    0000000000010000100001000000000000101000000000000000100000000000000...

    FingerprintsVector;ExtendedConnectivity:FunctionalClassAtomTypes:Radiu
    s2;57;AlphaNumericalValues;ValuesString;24769214 508787397 850393286 8
    62102353 981185303 1231636850 1649386610 1941540674 263599683 32920567
    1 571109041 639579325 683993318 723853089 810600886 885767127 90326012
    7 958841485 981022393 1126908698 1152248391 1317567065 1421489994 1455
    632544 1557272891 1826413669 1983319256 2015750777 2029559552 20404...

    FingerprintsVector;ExtendedConnectivity:EStateAtomTypes:Radius2;62;Alp
    haNumericalValues;ValuesString;25189973 528584866 662581668 671034184
    926543080 1347067490 1738510057 1759600920 2034425745 2097234755 21450
    44754 96779665 180364292 341712110 345278822 386540408 387387308 50430
    1706 617094135 771528807 957666640 997798220 1158349170 1291258082 134
    1138533 1395329837 1420277211 1479584608 1486476397 1487556246 1566...

    FingerprintsBitVector;MACCSKeyBits;166;BinaryString;Ascending;00000000
    0000000000000000000000000000000001001000010010000000010010000000011100
    0100101010111100011011000100110110000011011110100110111111111111011111
    11111111111110111000

    FingerprintsBitVector;MACCSKeyBits;322;BinaryString;Ascending;11101011
    1110011111100101111111000111101100110000000000000011100010000000000000
    0000000000000000000000000000000000000000000000101000000000000000000000
    0000000000000000000000000000000000000000000000000000000000000000000000
    0000000000000000000000000000000000000011000000000000000000000000000000
    0000000000000000000000000000000000000000

    FingerprintsVector;MACCSKeyCount;166;OrderedNumericalValues;ValuesStri
    ng;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 1 0 0 3 0 0 0 0 4 0 0 2 0 0 0 0 0 0 0 0 2 0 0 2 0 0 0 0
    0 0 0 0 1 1 8 0 0 0 1 0 0 1 0 1 0 1 0 3 1 3 1 0 0 0 1 2 0 11 1 0 0 0
    5 0 0 1 2 0 1 1 0 0 0 0 0 1 1 0 1 1 1 1 0 4 0 0 1 1 0 4 6 1 1 1 2 1 1
    3 5 2 2 0 5 3 5 1 1 2 5 1 2 1 2 4 8 3 5 5 2 2 0 3 5 4 1

    FingerprintsVector;MACCSKeyCount;322;OrderedNumericalValues;ValuesStri
    ng;14 8 2 0 2 0 4 4 2 1 4 0 0 2 5 10 5 2 1 0 0 2 0 5 13 3 28 5 5 3 0 0
    0 4 2 1 1 0 1 1 0 0 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 22 5 3 0 0 0 1 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 2 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...

    FingerprintsBitVector;PathLengthBits:AtomicInvariantsAtomTypes:MinLeng
    th1:MaxLength8;1024;BinaryString;Ascending;001000010011010101011000110
    0100010101011000101001011100110001000010001001101000001001001001001000
    0010110100000111001001000001001010100100100000000011000000101001011100
    0010000001000101010100000100111100110111011011011000000010110111001101
    0101100011000000010001000011000010100011101100001000001000100000000...

    FingerprintsVector;PathLengthCount:AtomicInvariantsAtomTypes:MinLength
    1:MaxLength8;432;NumericalValues;IDsAndValuesPairsString;C.X1.BO1.H3 2
    C.X2.BO2.H2 4 C.X2.BO3.H1 14 C.X3.BO3.H1 3 C.X3.BO4 10 F.X1.BO1 1 N.X
    2.BO2.H1 1 N.X3.BO3 1 O.X1.BO1.H1 3 O.X1.BO2 2 C.X1.BO1.H3C.X3.BO3.H1
    2 C.X2.BO2.H2C.X2.BO2.H2 1 C.X2.BO2.H2C.X3.BO3.H1 4 C.X2.BO2.H2C.X3.BO
    4 1 C.X2.BO2.H2N.X3.BO3 1 C.X2.BO3.H1:C.X2.BO3.H1 10 C.X2.BO3.H1:C....

    FingerprintsVector;PathLengthCount:MMFF94AtomTypes:MinLength1:MaxLengt
    h8;463;NumericalValues;IDsAndValuesPairsString;C5A 2 C5B 2 C=ON 1 CB 1
    8 COO 1 CR 9 F 1 N5 1 NC=O 1 O=CN 1 O=CO 1 OC=O 1 OR 2 C5A:C5B 2 C5A:N
    5 2 C5ACB 1 C5ACR 1 C5B:C5B 1 C5BC=ON 1 C5BCB 1 C=ON=O=CN 1 C=ONNC=O 1
    CB:CB 18 CBF 1 CBNC=O 1 COO=O=CO 1 COOCR 1 COOOC=O 1 CRCR 7 CRN5 1 CR
    OR 2 C5A:C5B:C5B 2 C5A:C5BC=ON 1 C5A:C5BCB 1 C5A:N5:C5A 1 C5A:N5CR ...

    FingerprintsVector;TopologicalAtomPairs:AtomicInvariantsAtomTypes:MinD
    istance1:MaxDistance10;223;NumericalValues;IDsAndValuesString;C.X1.BO1
    .H3-D1-C.X3.BO3.H1 C.X2.BO2.H2-D1-C.X2.BO2.H2 C.X2.BO2.H2-D1-C.X3.BO3.
    H1 C.X2.BO2.H2-D1-C.X3.BO4 C.X2.BO2.H2-D1-N.X3.BO3 C.X2.BO3.H1-D1-...;
    2 1 4 1 1 10 8 1 2 6 1 2 2 1 2 1 2 2 1 2 1 5 1 10 12 2 2 1 2 1 9 1 3 1
    1 1 2 2 1 3 6 1 6 14 2 2 2 3 1 3 1 8 2 2 1 3 2 6 1 2 2 5 1 3 1 23 1...

    FingerprintsVector;TopologicalAtomPairs:FunctionalClassAtomTypes:MinDi
    stance1:MaxDistance10;144;NumericalValues;IDsAndValuesString;Ar-D1-Ar
    Ar-D1-Ar.HBA Ar-D1-HBD Ar-D1-Hal Ar-D1-None Ar.HBA-D1-None HBA-D1-NI H
    BA-D1-None HBA.HBD-D1-NI HBA.HBD-D1-None HBD-D1-None NI-D1-None No...;
    23 2 1 1 2 1 1 1 1 2 1 1 7 28 3 1 3 2 8 2 1 1 1 5 1 5 24 3 3 4 2 13 4
    1 1 4 1 5 22 4 4 3 1 19 1 1 1 1 1 2 2 3 1 1 8 25 4 5 2 3 1 26 1 4 1 ...

    FingerprintsVector;TopologicalAtomTorsions:AtomicInvariantsAtomTypes;3
    3;NumericalValues;IDsAndValuesString;C.X1.BO1.H3-C.X3.BO3.H1-C.X3.BO4-
    C.X3.BO4 C.X1.BO1.H3-C.X3.BO3.H1-C.X3.BO4-N.X3.BO3 C.X2.BO2.H2-C.X2.BO
    2.H2-C.X3.BO3.H1-C.X2.BO2.H2 C.X2.BO2.H2-C.X2.BO2.H2-C.X3.BO3.H1-O...;
    2 2 1 1 2 2 1 1 3 4 4 8 4 2 2 6 2 2 1 2 1 1 2 1 1 2 6 2 4 2 1 3 1

    FingerprintsVector;TopologicalAtomTorsions:EStateAtomTypes;36;Numerica
    lValues;IDsAndValuesString;aaCH-aaCH-aaCH-aaCH aaCH-aaCH-aaCH-aasC aaC
    H-aaCH-aasC-aaCH aaCH-aaCH-aasC-aasC aaCH-aaCH-aasC-sF aaCH-aaCH-aasC-
    ssNH aaCH-aasC-aasC-aasC aaCH-aasC-aasC-aasN aaCH-aasC-ssNH-dssC a...;
    4 4 8 4 2 2 6 2 2 2 4 3 2 1 3 3 2 2 2 1 2 1 1 1 2 1 1 1 1 1 1 1 2 1 1 2

    FingerprintsVector;TopologicalAtomTriplets:AtomicInvariantsAtomTypes:M
    inDistance1:MaxDistance10;3096;NumericalValues;IDsAndValuesString;C.X1
    .BO1.H3-D1-C.X1.BO1.H3-D1-C.X3.BO3.H1-D2 C.X1.BO1.H3-D1-C.X2.BO2.H2-D1
    0-C.X3.BO4-D9 C.X1.BO1.H3-D1-C.X2.BO2.H2-D3-N.X3.BO3-D4 C.X1.BO1.H3-D1
    -C.X2.BO2.H2-D4-C.X2.BO2.H2-D5 C.X1.BO1.H3-D1-C.X2.BO2.H2-D6-C.X3....;
    1 2 2 2 2 2 2 2 8 8 4 8 4 4 2 2 2 2 4 2 2 2 4 2 2 2 2 1 2 2 4 4 4 2 2
    2 4 4 4 8 4 4 2 4 4 4 2 4 4 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 8...

    FingerprintsVector;TopologicalAtomTriplets:SYBYLAtomTypes:MinDistance1
    :MaxDistance10;2332;NumericalValues;IDsAndValuesString;C.2-D1-C.2-D9-C
    .3-D10 C.2-D1-C.2-D9-C.ar-D10 C.2-D1-C.3-D1-C.3-D2 C.2-D1-C.3-D10-C.3-
    D9 C.2-D1-C.3-D2-C.3-D3 C.2-D1-C.3-D2-C.ar-D3 C.2-D1-C.3-D3-C.3-D4 C.2
    -D1-C.3-D3-N.ar-D4 C.2-D1-C.3-D3-O.3-D2 C.2-D1-C.3-D4-C.3-D5 C.2-D1-C.
    3-D5-C.3-D6 C.2-D1-C.3-D5-O.3-D4 C.2-D1-C.3-D6-C.3-D7 C.2-D1-C.3-D7...

    FingerprintsVector;TopologicalPharmacophoreAtomPairs:ArbitrarySize:Min
    Distance1:MaxDistance10;54;NumericalValues;IDsAndValuesString;H-D1-H H
    -D1-NI HBA-D1-NI HBD-D1-NI H-D2-H H-D2-HBA H-D2-HBD HBA-D2-HBA HBA-D2-
    HBD H-D3-H H-D3-HBA H-D3-HBD H-D3-NI HBA-D3-NI HBD-D3-NI H-D4-H H-D4-H
    BA H-D4-HBD HBA-D4-HBA HBA-D4-HBD HBD-D4-HBD H-D5-H H-D5-HBA H-D5-...;
    18 1 2 1 22 12 8 1 2 18 6 3 1 1 1 22 13 6 5 7 2 28 9 5 1 1 1 36 16 10
    3 4 1 37 10 8 1 35 10 9 3 3 1 28 7 7 4 18 16 12 5 1 2 1

    FingerprintsVector;TopologicalPharmacophoreAtomPairs:FixedSize:MinDist
    ance1:MaxDistance10;150;OrderedNumericalValues;ValuesString;18 0 0 1 0
    0 0 2 0 0 1 0 0 0 0 22 12 8 0 0 1 2 0 0 0 0 0 0 0 0 18 6 3 1 0 0 0 1
    0 0 1 0 0 0 0 22 13 6 0 0 5 7 0 0 2 0 0 0 0 0 28 9 5 1 0 0 0 1 0 0 1 0
    0 0 0 36 16 10 0 0 3 4 0 0 1 0 0 0 0 0 37 10 8 0 0 0 0 1 0 0 0 0 0 0
    0 35 10 9 0 0 3 3 0 0 1 0 0 0 0 0 28 7 7 4 0 0 0 0 0 0 0 0 0 0 0 18...

    FingerprintsVector;TopologicalPharmacophoreAtomTriplets:ArbitrarySize:
    MinDistance1:MaxDistance10;696;NumericalValues;IDsAndValuesString;Ar1-
    Ar1-Ar1 Ar1-Ar1-H1 Ar1-Ar1-HBA1 Ar1-Ar1-HBD1 Ar1-H1-H1 Ar1-H1-HBA1 Ar1
    -H1-HBD1 Ar1-HBA1-HBD1 H1-H1-H1 H1-H1-HBA1 H1-H1-HBD1 H1-HBA1-HBA1 H1-
    HBA1-HBD1 H1-HBA1-NI1 H1-HBD1-NI1 HBA1-HBA1-NI1 HBA1-HBD1-NI1 Ar1-...;
    46 106 8 3 83 11 4 1 21 5 3 1 2 2 1 1 1 100 101 18 11 145 132 26 14 23
    28 3 3 5 4 61 45 10 4 16 20 7 5 1 3 4 5 3 1 1 1 1 5 4 2 1 2 2 2 1 1 1
    119 123 24 15 185 202 41 25 22 17 3 5 85 95 18 11 23 17 3 1 1 6 4 ...

    FingerprintsVector;TopologicalPharmacophoreAtomTriplets:FixedSize:MinD
    istance1:MaxDistance10;2692;OrderedNumericalValues;ValuesString;46 106
    8 3 0 0 83 11 4 0 0 0 1 0 0 0 0 0 0 0 0 21 5 3 0 0 1 2 2 0 0 1 0 0 0
    0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 101 18 11 0 0 145 132 26
    14 0 0 23 28 3 3 0 0 5 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 61 45 10 4 0
    0 16 20 7 5 1 0 3 4 5 3 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 5 ...

=head1 OPTIONS

=over 4

=item B<-a, --all>

List all the available information.

=item B<--AverageBitDensity>

List average bit density of fingerprint bit-vector strings.

=item B<--BitDensity>

List bit density of fingerprints bit-vector strings data in each row.

=item B<--count>

List number of data entries containing fingerprints bit-vector or vector strings data. This
is B<default behavior>.

=item B<-c, --ColMode> I<ColNum | ColLabel>

Specify how columns are identified in CSV/TSV I<TextFile(s)>: using column number or column
label. Possible values: I<ColNum or ColLabel>. Default value: I<ColNum>

=item B<-d, --detail> I<InfoLevel>

Level of information to print about lines being ignored. Default: I<1>. Possible values:
I<1, 2 or 3>.

=item B<--DataCheck>

Validate fingerprints data specified using B<--FingerprintsCol> and list information
about missing and invalid data.

=item B<-e, --empty>

List number of rows containing no fingerprints data.

=item B<--FingerprintsCol> I<col number | col name>

This value is B<-c, --colmode> specific. It corresponds to column in CSV/TSV I<TextFile(s)>
containing fingerprints data. Possible values: I<col number or col label>.
Default value: I<first column containing the word Fingerprints in its column label>.

=item B<--FingerprintsField> I<FieldLabel>

Fingerprints field label to use during listing of fingerprints information for I<SDFile(s)>.
Default value: I<first data field label containing the word Fingerprints in its label>.

=item B<--FingerprintsType>

List types of fingerprint strings: FingerprintsBitVector or FingerprintsVector.

=item B<--FingerprintsDescription>

List types of fingerprints: PathLengthBits, PathLengthCount, MACCSKeyCount,
ExtendedConnectivity and so on.

=item B<--FingerprintsSize>

List size of fingerprints.

=item B<--FingerprintsBitStringFormat>

List format of fingerprint bit-vector strings: BinaryString or HexadecimalString.

=item B<--FingerprintsBitOrder>

List order of bits data in fingerprint bit-vector bit strings: Ascending or Descending.

=item B<--FingerprintsVectorValuesType>

List type of values in fingerprint vector strings: OrderedNumericalValues, NumericalValues or
AlphaNumericalValues.

=item B<--FingerprintsVectorValuesFormat>

List format of values in fingerprint vector strings: ValuesString, IDsAndValuesString,
IDsAndValuesPairsString, ValuesAndIDsString or ValuesAndIDsPairsString.

=item B<-h, --help>

Print this help message.

=item B<--InDelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile(s)>. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a
delimiter.

=item B<--NumOfOnBits>

List number of on bits in fingerprints bit-vector strings data in each row.

=item B<--NumOfNonZeroValues>

List number of non-zero values in fingerprints vector strings data in each row.

=item B<-w, --WorkingDir> I<DirName>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To count number of lines containing fingerprints bit-vector or vector strings data present
in FP file, in a column name containing Fingerprint substring in text file, and in a data
field with Fingerprint substring in its label, type:

    % InfoFingerprintsFiles.pl SampleFPBin.csv

    % InfoFingerprintsFiles.pl SampleFPBin.sdf SampleFPBin.fpf
      SampleFPBin.csv

    % InfoFingerprintsFiles.pl SampleFPHex.sdf SampleFPHex.fpf
      SampleFPHex.csv

    % InfoFingerprintsFiles.pl SampleFPcount.sdf SampleFPcount.fpf
      SampleFPcount.csv

To list all available information about fingerprints bit-vector or vector strings data present
in FP file, in a column name containing Fingerprint substring in text file, and in a data
field with Fingerprint substring in its label, type:

    % InfoFingerprintsFiles.pl -a SampleFPHex.sdf SampleFPHex.fpf
      SampleFPHex.csv

    % InfoFingerprintsFiles.pl -a SampleFPcount.sdf SampleFPcount.fpf
      SampleFPcount.csv

To list all available information about fingerprints bit-vector or vector strings data present in a
column named Fingerprints in text file, type:

    % InfoFingerprintsFiles.pl -a --ColMode ColLabel --FingerprintsCol
      Fingerprints SampleFPHex.sdf

    % InfoFingerprintsFiles.pl -a --ColMode ColLabel --FingerprintsCol
      Fingerprints SampleFPcount.csv

To list all available information about fingerprints bit-vector or vector strings data present in a
data field names Fingerprints in SD file, type:

    % InfoFingerprintsFiles.pl -a --FingerprintsField Fingerprints
      SampleFPHex.sdf

    % InfoFingerprintsFiles.pl -a --FingerprintsField Fingerprints
      SampleFPcount.sdf

To list bit density, average bit density, and number of on bits for fingerprints bit-vector strings data
present in FP file, in a column name containing Fingerprint substring in text file, and in a data
field with Fingerprint substring in its label, type:

    % InfoFingerprintsFiles.pl --BitDensity --AverageBitDensity
      --NumOfOnBits SampleFPBin.csv SampleFPBin.sdf SampleFPBin.fpf

To list vector values type, format and number of non-zero values for fingerprints vector strings
data present in FP file, in a column name containing Fingerprint substring in text file, and in a data
field with Fingerprint substring in its label along with fingerprints type and description, type:

    % InfoFingerprintsFiles.pl --FingerprintsType --FingerprintsDescription
      --FingerprintsVectorValuesType --FingerprintsVectorValuesFormat
      --NumOfNonZeroValues SampleFPcount.csv SampleFPcount.sdf
      SampleFPcount.fpf

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

SimilarityMatricesFingerprints.pl, SimilaritySearchingFingerprints.pl, AtomNeighborhoodsFingerprints.pl,
AtomNeighborhoodsFingerprints.pl, ExtendedConnectivityFingerprints.pl, MACCSKeysFingerprints.pl,
PathLengthFingerprints.pl, TopologicalAtomPairsFingerprints.pl, TopologicalAtomTorsionsFingerprints.pl,
TopologicalPharmacophoreAtomPairsFingerprints.pl, TopologicalPharmacophoreAtomTripletsFingerprints.pl

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
