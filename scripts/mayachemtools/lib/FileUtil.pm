package FileUtil;
#
# File: FileUtil.pm
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
use Carp;
use File::stat;
use Time::localtime ();
use TextUtil ();

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw(CheckFileType ConvertCygwinPath ExpandFileNames FileModificationTimeAndDate FormattedFileModificationTimeAndDate FileSize FormatFileSize GetMayaChemToolsLibDirName GetUsageFromPod ParseFileName);
@EXPORT_OK = qw();
%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup package variables...
my($MayaChemToolsLibDir);

# Check to see path contains cygdrive and convert it into windows path...
sub ConvertCygwinPath {
  my($Path) = @_;
  my($NewPath, $OSName);

  $NewPath = $Path; $OSName = $^O;
  if ($OSName =~ /cygwin/i && $Path =~ /cygdrive/i ) {
    my(@PathParts) = split "\/", $Path;
    my($Drive) = $PathParts[2];
    shift @PathParts; shift @PathParts; shift @PathParts;
    $NewPath = join "\/", @PathParts;
    $NewPath = $Drive. ":\/" . $NewPath;
  }
  return $NewPath;
}

# Based on the file name extension, figure out its type.
sub CheckFileType {
  my($FileName, $FileExts) = @_;
  my($Status, @FileExtsList, $Index, $Ext);

  $Status = 0;
  @FileExtsList = split " ", $FileExts;
  for $Index (0 .. $#FileExtsList) {
    $Ext = $FileExtsList[$Index];
    if ($FileName =~ /(\.$Ext)$/i) {
      $Status = 1;
    }
  }
  return ($Status);
}

# Expand file names using specified directory and/or file names along with any
# file extensions containing one or more wild cards. And return the expanded
# list.
#
# IncludeDirName controls whether directory prefixes are included in expanded
# file names. Default is to always append directory name before expanded file
# name.
#
# Notes:
#   . Multiple file extensions are delimited by spaces.
#   . Wild card, *, is supported in directory and file names along with file
#     extensions.
#   . For a specified directory name in the files list, all the files in the
#     directory are retrieved using Perl opendir function and files filtered using file
#     extensions. The file names "." and ".." returned by opendir are ignored.
#   . For file names containing wild cards with and without any explicit file
#     extension specification in the file name, all the files in the directory are retrieved
#     using Perl opendir function and files filtered using file name along with any
#     file extension. The file names "." and ".." returned by opendir are ignored.
#
sub ExpandFileNames {
  my($Files, $FileExts, $IncludeDirName) = @_;
  my($FileName, $Index, $Delimiter, $FileExtsPattern, @FilesList, @DirFileNames);

  # Check whether to include directory name in expanded file names...
  $IncludeDirName = defined $IncludeDirName ? $IncludeDirName : 1;

  # Setup file externsions...
  $FileExtsPattern = "";
  if ($FileExts) {
    $FileExtsPattern = join "|", split " ", $FileExts;
    if ($FileExtsPattern =~ /\*/) {
      # Replace * by .*? for greedy match...
      $FileExtsPattern =~ s/\*/\.\*\?/g;
    }
  }

  @FilesList = ();

  FILEINDEX: for ($Index = 0; $Index < @$Files; $Index++) {
    $FileName = @$Files[$Index];
    $Delimiter = "\/";
    if ($FileName =~ /\\/ ) {
      $Delimiter = "\\";
    }

    if (-d $FileName) {
      my($DirName, $DirNamePrefix);

      $DirName = $FileName;
      $DirNamePrefix = $IncludeDirName ? "$DirName$Delimiter" : "";

      # glob doesn't appear to work during command line invocation from Windows.
      # So, use opendir to make it work...
      #
      # push @FilesList,  map {glob("$DirName/*.$_")} split " ", $FileExts;
      #
      @DirFileNames = ();
      if (!opendir DIRNAME, $DirName) {
	carp "Warning: Ignoring directory $DirName: Couldn't open it: $! ...";
	next FILEINDEX;
      }

      # Collect file names without '.' and '..' as readdir function places them on the list...
      #
      @DirFileNames = map { "$DirNamePrefix$_"  } grep { !/^(\.|\.\.)$/ } readdir DIRNAME;
      closedir DIRNAME;

      # Collect files with any specified file extensions...
      if ($FileExtsPattern) {
	@DirFileNames = grep { /\.$FileExtsPattern$/ } @DirFileNames;
      }

      push @FilesList, @DirFileNames;
    }
    elsif ($FileName =~ /\*/) {
      my($FileDir, $Name, $FileExt, $DirNamePrefix);

      # Filenames are not expanded during command line invocation from Windows...
      ($FileDir, $Name, $FileExt) = ParseFileName($FileName);

      $DirNamePrefix = $IncludeDirName ? "$FileDir$Delimiter" : "";

      @DirFileNames = ();
      if (!opendir FILEDIR, $FileDir) {
	carp "Warning: Ignoring files $FileName: Couldn't open directory $FileDir: $! ...";
	next FILEINDEX;
      }

      # Collect file names without '.' and '..' as readdir function places them on the list...
      #
      @DirFileNames = map { "$DirNamePrefix$_"  } grep { !/^(\.|\.\.)$/ } readdir FILEDIR;
      closedir FILEDIR;

      if (length($Name) > 1) {
	# Replace * by .*? for greedy match...
	$Name =~ s/\*/\.\*\?/g;
	@DirFileNames =  grep { /$Name/ } @DirFileNames;
      }

      if ($FileExt) {
	$FileExt =~ s/\*/\.\*\?/g;
	@DirFileNames =  grep { /\.$FileExt$/ } @DirFileNames;
      }
      elsif ($FileExtsPattern) {
	@DirFileNames = grep { /\.$FileExtsPattern$/ } @DirFileNames;
      }

      push @FilesList, @DirFileNames;
    }
    else {
      push @FilesList, $FileName;
    }
  }
  return @FilesList;
}

# Formatted file modification time...
sub FormattedFileModificationTimeAndDate {
  my($FileName) = @_;
  my($TimeString, $DateString) = ('') x 2;

  if (! -e $FileName) {
    return ($TimeString, $DateString);
  }
  my($Hours, $Mins, $Secs, $DayName, $MonthName, $Month, $Year) = FileModificationTimeAndDate($FileName);

  # Setup time suffix...
  my($TimeSuffix) = '';
  if ($Hours < 12) {
    $TimeSuffix = 'AM';
  }
  elsif ($Hours > 12) {
    $TimeSuffix = 'PM';
    $Hours = $Hours - 12;
  }
  elsif ($Hours == 12 && ($Mins > 0 || $Secs > 0)) {
    $TimeSuffix = 'PM';
  }
  elsif ($Hours == 12 && $Mins == 0 && $Secs == 0) {
    $TimeSuffix = 'Noon';
  }

  $Month = TextUtil::AddNumberSuffix($Month);

  $TimeString = "${DayName} ${Hours}:${Mins}:${Secs} ${TimeSuffix}";
  $DateString = "${MonthName} ${Month}, ${Year}";

  return ($TimeString, $DateString);
}

# File modifcation time and date...
sub FileModificationTimeAndDate {
  my($FileName) = @_;
  my($Hours, $Mins, $Secs, $DayName, $MonthName, $Month, $Year) = ('') x 7;

  if (! -e $FileName) {
    return ($Hours, $Mins, $Secs, $DayName, $MonthName, $Month, $Year);
  }

  my($CTimeString, $FileStatRef, $TimeStamp);
  $FileStatRef = stat($FileName);

  $CTimeString = Time::localtime::ctime($FileStatRef->mtime);

  # ctime returns: Thu Aug 3 10:13:53 2006
  ($DayName, $MonthName, $Month, $TimeStamp, $Year) = split /[ ]+/, $CTimeString;
  ($Hours, $Mins, $Secs) = split /\:/, $TimeStamp;

  return ($Hours, $Mins, $Secs, $DayName, $MonthName, $Month, $Year);
}

# Format file size...
sub FormatFileSize {
  my($Precision, $Size);

  $Precision = 1;
  if (@_ == 2) {
    ($Size, $Precision) = @_;
  }
  else {
    ($Size) = @_;
  }
  my($SizeDenominator, $SizeSuffix);
  FORMAT: {
      if ($Size < 1024) { $SizeDenominator = 1; $SizeSuffix = 'bytes'; last FORMAT;}
      if ($Size < (1024*1024)) { $SizeDenominator = 1024; $SizeSuffix = 'KB'; last FORMAT;}
      if ($Size < (1024*1024*1024)) { $SizeDenominator = 1024*1024; $SizeSuffix = 'MB'; last FORMAT;}
      if ($Size < (1024*1024*1024*1024)) { $SizeDenominator = 1024*1024*1024; $SizeSuffix = 'GB'; last FORMAT;}
      $SizeDenominator = 1; $SizeSuffix = 'bytes';
    }
  $Size /= $SizeDenominator;
  $Size = sprintf("%.${Precision}f", $Size) + 0;
  $Size = "$Size $SizeSuffix";

  return $Size;
}

# Get file size in bytes...
sub FileSize {
  my($File) = @_;

  if (! -e $File) {
    return undef;
  }
  return (-s $File)
}

# Get MayaChemTool's lib directory name using @INC to extract
# <MAYACHEMTOOLS>/lib directory location: first entry in @INC path should contain
# MayaChemTools modules location
#
sub GetMayaChemToolsLibDirName {

  if (defined $MayaChemToolsLibDir) {
    return $MayaChemToolsLibDir;
  }

  $MayaChemToolsLibDir = "";
  if ($INC[0] =~ /MayaChemTools/i) {
    $MayaChemToolsLibDir = $INC[0];
  }
  else {
    # Go through rest of the entries...
    my($Index);
    INDEX: for $Index (1 .. $#INC) {
      if ($INC[$Index] =~ /MayaChemTools/i) {
	$MayaChemToolsLibDir = $INC[$Index];
	last INDEX;
      }
    }
    if (!$MayaChemToolsLibDir) {
      carp "Warning: MayaChemTools lib directory location doesn't appear to exist in library search path specified by \@INC ...";
    }
  }
  return $MayaChemToolsLibDir;
}

# Get Usage from Pod...
sub GetUsageFromPod {
  my($Usage, $ScriptPath);

  ($ScriptPath) = @_;
  $Usage = "Script usage not available: pod2text or pod2text.bat doesn't exist in your Perl installation and direct invocation of Pod::Text also failed\n";

  # Get pod documentation: try pod2text first followed by perdoc.bat in case it fails to
  # to handle ActiveState Perl...
  my($PodStatus);
  $PodStatus = (open CMD, "pod2text $ScriptPath|") ? 1 : ((open CMD, "pod2text.bat $ScriptPath|") ? 1 : 0);
  if (!$PodStatus) {
    # Try direct invocation of Pod::Text before giving up...
    my($PodTextCmd);
    $PodTextCmd = "perl -e \'use Pod::Text (); \$TextFormatter = Pod::Text->new(); \$TextFormatter->parse_from_file(\"$ScriptPath\");\'";
    $PodStatus = (open CMD, "$PodTextCmd|") ? 1 : 0;
    if (!$PodStatus) {
      return $Usage;
    }
  }
  my($ProcessingSection, $InParametersSection, $InOptionsSection, @LineWords);
  $ProcessingSection = 0; $InParametersSection = 0; $InOptionsSection = 0;
  PODLINE: while (<CMD>) {
    if (/^SYNOPSIS/) {
      $_ = <CMD>; chomp; s/^ +//g;
      (@LineWords) = split / /;
      $Usage = qq(Usage: $LineWords[0] [-options]... );
      shift @LineWords;
      $Usage .= join(" ", @LineWords) . "\n";
    }
    elsif (/^(DESCRIPTION|PARAMETERS|OPTIONS|EXAMPLES|AUTHOR|SEE ALSO|COPYRIGHT)/i) {
      # Various sections...
      chomp;
      $Usage .= ucfirst(lc($_)) . ":\n";
      $ProcessingSection = 1;
      $InOptionsSection = /^OPTIONS/ ? 1 : 0;
      $InParametersSection = /^PARAMETERS/ ? 1 : 0;
    }
    elsif ($InParametersSection|$InOptionsSection) {
      if (/^[ ]+\-/ || /^[ ]{4,4}/) {
	# Start of option line: any number of spaces followed by - sign.
	# Put back in <> which pod2text replaced to **
	my($OptionLine) = qq($_);
           OPTIONLINE: while (<CMD>) {
	  if (/^(    )/) {
	    $OptionLine .= qq($_);
	  }
	  else {
	    $OptionLine =~ s/\*(([a-zA-Z0-9])|(\[)|(\#)|(\"))/"\<" . substr($&, -1, 1)/e;
	    $OptionLine =~ s/(([a-zA-Z0-9])|(\])|(\#)|(\"))\*/substr($&, 0, 1) . "\>"/e;
	    $Usage .= qq($OptionLine$_);
	    last OPTIONLINE;
	  }
	}
      }
    }
    else {
      if ($ProcessingSection) { $Usage .= qq($_); }
    }
  }
  close CMD;

  # Take out **which pod2text puts in for <>
  $Usage =~ s/\*(([a-zA-Z0-9;#-])|(\")|(\()|(\[)|(\.))/substr($&, -1, 1)/eg;
  $Usage =~ s/(([a-zA-Z0-9;#-])|(\")|(\))|(\])|(\.))\*/substr($&, 0, 1)/eg;

  return $Usage;
}

# Split full file name into directory path, file name, and the extension.
sub ParseFileName {
  my($FullName) = @_;
  my($FileDir, $FileName, $FileExt, @FullFileNameParts, @FileNameParts, $Delimiter);

  $Delimiter = "\/";
  if ($FullName =~ /\\/ ) {
    $Delimiter = "\\";
    $FullName =~ s/\\/\//g;
  }
  $FileDir = ""; $FileName = ""; $FileExt = "";
  @FullFileNameParts = (); @FileNameParts = ();

  @FullFileNameParts = split "\/", $FullName;
  @FileNameParts = split /\./, $FullFileNameParts[$#FullFileNameParts];

  # Setup file dir...
  if (@FullFileNameParts == 1) {
    $FileDir = "\.";
  }
  else {
    pop @FullFileNameParts;
    $FileDir = join $Delimiter, @FullFileNameParts;
  }

  # Setup file name and ext...
  if (@FileNameParts == 1) {
    $FileName = $FileNameParts[0];
    $FileExt = "";
  }
  elsif (@FileNameParts == 2) {
    $FileName = $FileNameParts[0];
    $FileExt = $FileNameParts[1];
  }
  elsif (@FileNameParts > 2) {
    # Use the last entry as file extension and the rest for file name...
    $FileExt = $FileNameParts[$#FileNameParts];
    pop @FileNameParts;
    $FileName = join '.', @FileNameParts;
  }
  return ($FileDir, $FileName, $FileExt);
}

1;

__END__

=head1 NAME

FileUtil

=head1 SYNOPSIS

use FileUtil;

use FileUtil qw(:all);

=head1 DESCRIPTION

B<FileUtil> module provides the following functions:

CheckFileType, ConvertCygwinPath, ExpandFileNames, FileModificationTimeAndDate,
FileSize, FormatFileSize, FormattedFileModificationTimeAndDate,
GetMayaChemToolsLibDirName, GetUsageFromPod, ParseFileName

=head1 FUNCTIONS

=over 4

=item B<CheckFileType>

    $Status = CheckFileType($FileName, $FileExts);

Based on I<FileExts>, decides type of I<FileName> and return 1 or 0.

=item B<ConvertCygwinPath>

    $NewPath = ConvertCygwinPath($Path);

Check to see whether I<Path> contains any Cygwin drive specification and convert
it into Windows path.

=item B<ExpandFileNames>

    @FilesList = ExpandFileNames(\@Files, $FileExts);
    @FilesList = ExpandFileNames(\@Files, $FileExts, $IncludeDirName);

For each directory name or wild card file name in I<Files>, generate all file names which
correspond to the specification along with match to any extensions in I<FileExts> and return an
array B<FileList> containing these file names and other names. I<IncludeDirName> controls
controls whether directory prefixes are included in expanded file names. Default is to always
append directory name before expanded file name.

Notes:

    . Multiple file extensions are delimited by spaces.
    . Wild card, *, is supported in directory and file names along with file
      extensions.
    . For a specified directory name in the files list, all the files in the
      directory are retrieved using Perl opendir function and files filtered using file
      extensions. The file names "." and ".." returned by opendir are ignored.
    . For file names containing wild cards with and without any explicit file
      extension specification in the file name, all the files in the directory are retrieved
      using Perl opendir function and files filtered using file name along with any
      file extension. The file names "." and ".." returned by opendir are ignored.

=item B<FormattedFileModificationTimeAndDate>

    ($TimeString, $DateString) =
         FormattedFileModificationTimeAndDate($FileName);

Returns a formatted time and date string corresponding to I<FileName> modification time.

=item B<FileModificationTimeAndDate>

    ($Hours, $Mins, $Secs, $DayName, $MonthName, $Month, $Year) =
         FileModificationTimeAndDate($FileName);

Returns file modification time and date values for specified I<FileName>.

=item B<FormatFileSize>

    $FormattedSize= FormatFileSize($Size, [$Precision]);

Formats the file size in bytes to human readable value and returns a formatted file
size string.

=item B<FileSize>

    $Size= FileSize($FileName);

Returns size of I<FileName> in bytes

=item B<GetMayaChemToolsLibDirName>

    $MayaChemToolsLibDir = GetMayaChemToolsLibDirName();

Returns MayaChemTools lib directory name by parsing B<INC> values to extract
B<MAYACHEMTOOLS/lib> directory location: first entry in B<INC> path should contain
MayaChemTools lib location.

=item B<GetUsageFromPod>

    $ScriptUsage = GetUsageFromPod($AbsoluteScriptPath);

Generates a B<ScriptUsage> string from pod documentation in the script file using
pod2text or perdoc.bat Perl utitities.

=item B<ParseFileName>

    ($FileDir, $FileName, $FileExt) = ParseFileName($FullFileName);

Splits I<FullFileName> into directory name, file name, and extension. B<FileDir> is
set to current directory for absent directory name in I<FullFileName>. And I<FileExt>
is set to NULL string for I<FullFileName> without any extension.

This function doesn't perform checking ragarding the presence of the directory I<FileDir>
and I<FullFileName> and the I<FullFileName> without any extension is assumed to be
a file instead of a directory.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

TextUtil.pm, TimeUtil.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
