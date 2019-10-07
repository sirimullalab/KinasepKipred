package SequenceFileUtil;
#
# File: SequenceFileUtil.pm
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
use Text::ParseWords;
use TextUtil;
use FileUtil;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw(AreSequenceLengthsIdentical CalcuatePercentSequenceIdentity CalculatePercentSequenceIdentityMatrix GetLongestSequence GetShortestSequence GetSequenceLength IsGapResidue IsSupportedSequenceFile IsClustalWSequenceFile IsPearsonFastaSequenceFile IsMSFSequenceFile ReadSequenceFile RemoveSequenceGaps RemoveSequenceAlignmentGapColumns WritePearsonFastaSequenceFile);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Compare lengths of all sequences...
sub AreSequenceLengthsIdentical {
  my($SequencesDataRef) = @_;
  my($Status, $ID, $FirstID, $FirstSeqLen, $FirstDifferentLenID, $SeqLen);

  $Status = 1;
  $FirstID = '';
  $FirstDifferentLenID = '';

  ID: for $ID (@{$SequencesDataRef->{IDs}}) {
    if (!$FirstID) {
      $FirstID = $ID;
      $FirstSeqLen = length($SequencesDataRef->{Sequence}{$ID});
      next ID;
    }
    $SeqLen = length($SequencesDataRef->{Sequence}{$ID});
    if ($SeqLen != $FirstSeqLen) {
      $Status = 0;
      $FirstDifferentLenID = $ID;
      last ID;
    }
  }
  return ($Status);
}

# Calculate percent identity between two sequences. By default, gaps are ignored.
sub CalcuatePercentSequenceIdentity {
  my($Sequence1, $Sequence2, $PercentIdentity, $IgnoreGaps, $Precision);

  $PercentIdentity = '';
  $Precision = 1;
  $IgnoreGaps = 1;
  if (@_ == 4) {
    ($Sequence1, $Sequence2, $IgnoreGaps, $Precision) = @_;
  }
  elsif (@_ == 3) {
    ($Sequence1, $Sequence2, $IgnoreGaps) = @_;
  }
  elsif (@_ == 2) {
    ($Sequence1, $Sequence2) = @_;
  }
  else {
    return $PercentIdentity;
  }
  if (!(IsNotEmpty($Sequence1) && IsNotEmpty($Sequence2))) {
    return $PercentIdentity;
  }
  my($Index, $Identity, $Sequence1Len, $Sequence2Len, $Residue1, $Residue2, $ResMatchCount, $ResCount);

  $Sequence1Len = length($Sequence1);
  $Sequence2Len = length($Sequence2);

  $ResMatchCount = 0;
  $ResCount = 0;
  RESIDUE: for $Index (0 .. ($Sequence1Len - 1)) {
    $Residue1 = substr($Sequence1, $Index, 1);
    $Residue2 = ($Index < $Sequence2Len) ? substr($Sequence2, $Index, 1) : '';
    if ($IgnoreGaps) {
      if ($Residue1 !~ /[A-Z]/i || $Residue2 !~ /[A-Z]/i) {
	next RESIDUE;
      }
    }
    if ($Residue1 eq $Residue2) {
      $ResMatchCount++;
    }
    $ResCount++;
  }
  $Identity = $ResCount ? ($ResMatchCount/$ResCount) : 0.0;
  $PercentIdentity = sprintf("%.${Precision}f", ($Identity * 100));

  return $PercentIdentity;
}

# Calculate pairwise identify matrix for all the sequences and return a reference
# to a hash with the following keys:
#
# {IDs} - Sequence IDs
# {Count} - Number of IDs
# {PercentIdentity}{$RowID}{$ColID} - Percent identify for a pair of sequences
#
sub CalculatePercentSequenceIdentityMatrix {
  my($SequencesDataRef, $IgnoreGaps, , $Precision, $ID, $RowID, $ColID, $RowIDSeq, $ColIDSeq, $PercentIdentity, %IdentityMatrixData);

  $IgnoreGaps = 1;
  $Precision = 1;
  if (@_ == 3) {
    ($SequencesDataRef, $IgnoreGaps, $Precision) = @_;
  }
  elsif (@_ == 2) {
    ($SequencesDataRef, $IgnoreGaps) = @_;
  }
  else {
    ($SequencesDataRef) = @_;
  }

  %IdentityMatrixData = ();
  @{$IdentityMatrixData{IDs}} = ();
  %{$IdentityMatrixData{PercentIdentity}} = ();
  $IdentityMatrixData{Count} = 0;

  for $ID (@{$SequencesDataRef->{IDs}}) {
    push @{$IdentityMatrixData{IDs}}, $ID;
    $IdentityMatrixData{Count} += 1;
  }
  # Initialize and calculate percent identity data values...
  for $RowID (@{$SequencesDataRef->{IDs}}) {
    %{$IdentityMatrixData{PercentIdentity}{$RowID}} = ();
    $RowIDSeq = $SequencesDataRef->{Sequence}{$RowID};
    for $ColID (@{$SequencesDataRef->{IDs}}) {
      $IdentityMatrixData{$RowID}{$ColID} = '';
      $ColIDSeq = $SequencesDataRef->{Sequence}{$ColID};
      $PercentIdentity = CalcuatePercentSequenceIdentity($RowIDSeq, $ColIDSeq, $IgnoreGaps, $Precision);
      $IdentityMatrixData{PercentIdentity}{$RowID}{$ColID} = $PercentIdentity;
    }
  }
  return \%IdentityMatrixData;
}

# Retrieve information about shortest sequence...
sub GetShortestSequence {
  my($SequencesDataRef, $IgnoreGaps, $ID, $Sequence, $SeqLen, $Description);

  $IgnoreGaps = 1;
  if (@_ == 2) {
    ($SequencesDataRef, $IgnoreGaps) = @_;
  }
  else {
    ($SequencesDataRef) = @_;
  }

  ($ID, $Sequence, $SeqLen, $Description) =  _GetShortestOrLongestSequence($SequencesDataRef, 'Shortest', $IgnoreGaps);
  return ($ID, $Sequence, $SeqLen, $Description);
}

# Retrieve information about longest sequence..
sub GetLongestSequence {
  my($SequencesDataRef, $IgnoreGaps, $ID, $Sequence, $SeqLen, $Description);

  $IgnoreGaps = 1;
  if (@_ == 2) {
    ($SequencesDataRef, $IgnoreGaps) = @_;
  }
  else {
    ($SequencesDataRef) = @_;
  }

  ($ID, $Sequence, $SeqLen, $Description) =  _GetShortestOrLongestSequence($SequencesDataRef, 'Longest', $IgnoreGaps);
  return ($ID, $Sequence, $SeqLen, $Description);
}

# Get sequence length...
sub GetSequenceLength {
  my($Seq, $SeqLen, $IgnoreGaps);

  $SeqLen = ''; $IgnoreGaps = 1;
  if (@_ == 2) {
    ($Seq, $IgnoreGaps) = @_;
  }
  else {
    ($Seq) = @_;
  }
  if ($IgnoreGaps) {
    my($Index, $Residue);
    $SeqLen = 0;
    for $Index (0 .. (length($Seq) - 1)) {
      $Residue = substr($Seq, $Index, 1);
      if ($Residue =~ /[A-Z]/i) {
	$SeqLen++;
      }
    }
  }
  else {
    $SeqLen = length($Seq);
  }

  return $SeqLen;
}

# Is it a gap residue...
sub IsGapResidue {
  my($Residue) = @_;
  my($Status);

  $Status = ($Residue !~ /[A-Z]/i ) ? 1 : 0;

  return $Status;
}

# Is it a supported sequence file?
#
# Supported seqence formats are:
#
# ALN/ClustalW   .aln
# GCG/MSF         .msf
# PILEUP/MSF     .msf
# Fasts(Pearson) .fasta, .fta
# NBRF/PIR         .pir
#
sub IsSupportedSequenceFile {
  my($SequenceFile) = @_;
  my($Status, $SequenceFormat);
  $Status = 0; $SequenceFormat = 'NotSupported';

  SEQFORMAT: {
      if (IsClustalWSequenceFile($SequenceFile)) {$Status = 1; $SequenceFormat = 'ClustalW'; last SEQFORMAT}
      if (IsPearsonFastaSequenceFile($SequenceFile)) {$Status = 1; $SequenceFormat = 'Pearson'; last SEQFORMAT}
      if (IsPIRFastaSequenceFile($SequenceFile)) {$Status = 1; $SequenceFormat = 'PIR'; last SEQFORMAT}
      if (IsMSFSequenceFile($SequenceFile)) {$Status = 1; $SequenceFormat = 'MSF'; last SEQFORMAT}
      $Status = 0; $SequenceFormat = 'NotSupported';
  }
  return ($Status, $SequenceFormat);
}

# Is it a ClustalW multiple sequence sequence file...
sub IsClustalWSequenceFile {
  my($SequenceFile) = @_;
  my($Status, $Line);

  $Status = 0;

  open SEQUENCEFILE, "$SequenceFile" or die "Couldn't open $SequenceFile: $!\n";
  $Line = GetTextLine(\*SEQUENCEFILE);
  $Status = ($Line =~ /(ClustalW|Clustal W|Clustal)/i ) ? 1 : 0;
  close SEQUENCEFILE;

  return $Status;
}

# Is it a valid Pearson fasta sequence or alignment file?
#
sub IsPearsonFastaSequenceFile {
  my($FastaFile, $Line, $Status);

  ($FastaFile) = @_;
  $Status = 0;

  open FASTAFILE, "$FastaFile" or die "Couldn't open $FastaFile: $!\n";
  $Line = GetTextLine(\*FASTAFILE);

  # First line starts with > and the fourth character is not ';'; otherwise, it's
  # PIR FASTA format.
  if ($Line =~ /^>/) {
    my($FourthChar);
    $FourthChar = substr($Line, 3, 1);
    $Status = ($FourthChar !~ /\;/) ? 1 : 0;
  }
  close FASTAFILE;

  return $Status;
}

# Is it a valid NBRF/PIR fasta sequence or alignment file?
#
sub IsPIRFastaSequenceFile {
  my($FastaFile, $Line, $Status);

  ($FastaFile) = @_;
  $Status = 0;

  open FASTAFILE, "$FastaFile" or die "Couldn't open $FastaFile: $!\n";
  $Line = GetTextLine(\*FASTAFILE);

  # First line starts with > and the fourth character is ';'; otherwise, it's
  # a Pearson FASTA format.
  if ($Line =~ /^>/) {
    my($FourthChar);
    $FourthChar = substr($Line, 3, 1);
    $Status = ($FourthChar =~ /\;/) ? 1 : 0;
  }
  close FASTAFILE;

  return $Status;
}

# Is it a valid MSF sequence or alignment file?
#
sub IsMSFSequenceFile {
  my($MSFFile) = @_;

  open MSFFILE, "$MSFFile" or die "Couldn't open $MSFFile: $!\n";

  my($Line, $Status);

  $Status = 0;
  # Find a line that contains MSF: keyword and ends with '..'
  LINE: while ($Line = GetTextLine(\*MSFFILE)) {
    $Line = RemoveLeadingWhiteSpaces($Line);
    if ($Line =~ /MSF:/i && $Line =~ /\.\.[ ]*$/) {
      $Status = 1;
      last LINE;
    }
    elsif ($Line =~ /(!!AA_MULTIPLE_ALIGNMENT|!!NA_MULTIPLE_ALIGNMENT|PILEUP)/i) {
      # Pileup MSF...
      $Status = 1;
      last LINE;
    }
  }
  close MSFFILE;

  return $Status;
}

# Read sequence or sequence alignment file...
sub ReadSequenceFile {
  my($SequenceFile) = @_;

  if (IsPearsonFastaSequenceFile($SequenceFile)) {
    return ReadPearsonFastaSequenceFile($SequenceFile);
  }
  elsif (IsPIRFastaSequenceFile($SequenceFile)) {
    return ReadPIRFastaSequenceFile($SequenceFile);
  }
  elsif (IsMSFSequenceFile($SequenceFile)) {
    return ReadMSFSequenceFile($SequenceFile);
  }
  elsif (IsClustalWSequenceFile($SequenceFile)) {
    return ReadClustalWSequenceFile($SequenceFile);
  }
  else {
    return undef;
  }
}

# Read file and setup alignment data...
sub ReadClustalWSequenceFile {
  my($SequenceFile) = @_;

  return _ReadFileAndSetupSequencesData($SequenceFile, 'ClustalW');
}

# Read file and setup alignment data...
sub ReadPearsonFastaSequenceFile {
  my($SequenceFile) = @_;

  return _ReadFileAndSetupSequencesData($SequenceFile, 'Pearson');
}

# Read file and setup alignment data...
sub ReadPIRFastaSequenceFile {
  my($SequenceFile) = @_;

  return _ReadFileAndSetupSequencesData($SequenceFile, 'PIR');
}


# Read file and setup sequence data...
sub ReadMSFSequenceFile {
  my($SequenceFile) = @_;

  return _ReadFileAndSetupSequencesData($SequenceFile, 'MSF');
}

# Write out a Pearson FASTA file...
sub WritePearsonFastaSequenceFile {
  my($SequenceFileName, $SequenceDataRef, $MaxLength, $ID, $Description, $Sequence, $WrappedSequence);

  $MaxLength = 80;
  if (@_ == 3) {
    ($SequenceFileName, $SequenceDataRef, $MaxLength) = @_;
  }
  elsif (@_ == 2) {
    ($SequenceFileName, $SequenceDataRef) = @_;
  }
  open SEQUENCEFILE, ">$SequenceFileName" or die "Can't open $SequenceFileName: $!\n";
  for $ID (@{$SequenceDataRef->{IDs}}) {
    $Description = $SequenceDataRef->{Description}{$ID};
    $Sequence = $SequenceDataRef->{Sequence}{$ID};
    $WrappedSequence = WrapText($Sequence, $MaxLength, "\n");

    # Description also contains ID...
    print SEQUENCEFILE ">$Description\n";
    print SEQUENCEFILE "$WrappedSequence\n";
  }
  close SEQUENCEFILE;
}

# Get ID, Sequence and Length for smallest or longest sequence
sub _GetShortestOrLongestSequence {
  my($SequencesDataRef, $SequenceType, $IgnoreGaps) = @_;
  my($ID, $Seq, $SeqLen, $Description, $FirstID, $FirstSeqLen, $CurrentID, $CurrentSeq, $CurrentSeqLen, $CurrentDescription);

  ($ID, $Seq, $SeqLen) = ('', '', '');
  $FirstID = '';

  ID: for $CurrentID (@{$SequencesDataRef->{IDs}}) {
    $CurrentSeq = $IgnoreGaps ? RemoveSequenceGaps($SequencesDataRef->{Sequence}{$CurrentID}) : $SequencesDataRef->{Sequence}{$CurrentID};
    $CurrentSeqLen = GetSequenceLength($CurrentSeq, $IgnoreGaps);
    $CurrentDescription = $SequencesDataRef->{Description}{$CurrentID};
    if (!$FirstID) {
      $FirstID = $ID; $FirstSeqLen = $CurrentSeqLen;
      ($ID, $Seq, $SeqLen, $Description) = ($CurrentID, $CurrentSeq, $CurrentSeqLen, $CurrentDescription);
      next ID;
    }
    if ($CurrentSeqLen != $SeqLen) {
      if (($SequenceType =~ /Shortest/i) && ($CurrentSeqLen < $SeqLen)) {
	($ID, $Seq, $SeqLen, $Description) = ($CurrentID, $CurrentSeq, $CurrentSeqLen, $CurrentDescription);
      }
      elsif (($SequenceType =~ /Longest/i) && ($CurrentSeqLen > $SeqLen) ) {
	($ID, $Seq, $SeqLen, $Description) = ($CurrentID, $CurrentSeq, $CurrentSeqLen, $CurrentDescription);
      }
    }
  }
  return ($ID, $Seq, $SeqLen, $Description);
}

# Remove gaps in the sequence and return new sequence...
sub RemoveSequenceGaps {
  my($Seq) = @_;
  my($SeqWithoutGaps, $SeqLen, $Index, $Residue);

  $SeqWithoutGaps = '';
  $SeqLen = length($Seq);
  for $Index (0 .. ($SeqLen - 1)) {
    $Residue = substr($Seq, $Index, 1);
    if ($Residue =~ /[A-Z]/i) {
      $SeqWithoutGaps .= $Residue;
    }
  }

  return $SeqWithoutGaps;
}

# Using input alignment data map ref containing following keys, generate
# a new hash with same set of keys after residue columns containg only
# gaps have been removed:
#
# {IDs} : Array of IDs in order as they appear in file
# {Count}: ID count...
# {Description}{$ID} : Description data...
# {Sequence}{$ID} : Sequence data...
#
sub RemoveSequenceAlignmentGapColumns {
  my($ID, $AlignmentDataMapRef, %NewAlignmentDataMap);

  ($AlignmentDataMapRef) = @_;

  %NewAlignmentDataMap = ();
  @{$NewAlignmentDataMap{IDs}} =();
  %{$NewAlignmentDataMap{Description}} =();
  %{$NewAlignmentDataMap{Sequence}} =();
  $NewAlignmentDataMap{Count} = 0;

  # Transfer ID and count information...
  for $ID (@{$AlignmentDataMapRef->{IDs}}) {
    push @{$NewAlignmentDataMap{IDs}}, $ID;
    $NewAlignmentDataMap{Description}{$ID} = $AlignmentDataMapRef->{Description}{$ID};
    $NewAlignmentDataMap{Sequence}{$ID} = '';
    $NewAlignmentDataMap{Count} += 1;
  }

  # Go over residue columns and transfer the data...
  my($FirstID, $FirstSeq, $FirstSeqLen, $Index, $Res, $GapColumn);

  $FirstID = $AlignmentDataMapRef->{IDs}[0];
  $FirstSeq = $AlignmentDataMapRef->{Sequence}{$FirstID};
  $FirstSeqLen = length($FirstSeq);

  RES: for $Index (0 .. ($FirstSeqLen - 1)) {
    # Is this a gap column?
    $GapColumn = 1;
    ID: for $ID (@{$AlignmentDataMapRef->{IDs}}) {
      $Res = substr($AlignmentDataMapRef->{Sequence}{$ID}, $Index, 1);
      if ($Res =~ /[A-Z]/i) {
	$GapColumn = 0;
	last ID;
      }
    }
    if ($GapColumn) {
      next RES;
    }
    # Transfer this residue...
    for $ID (@{$AlignmentDataMapRef->{IDs}}) {
      $Res = substr($AlignmentDataMapRef->{Sequence}{$ID}, $Index, 1);
      $NewAlignmentDataMap{Sequence}{$ID} .= $Res;
    }
  }

  return (\%NewAlignmentDataMap);
}

#
# Read sequences file and return a reference to hash with the following keys:
#
# {IDs} - Array of sequence IDs
# {Count} - Number of sequences
# {Description}{$ID} - Sequence description
# {Sequence}{$ID} - Sequence for a specific ID
# {InputFileType} - Sequence file format
# {ConservedAnnotation} - Conserved residue annonation
#
# Note:
#   . Conserved residue annotation either exist in the input sequence alignment file or set
#     for a file containing same number of residues for all the sequence using the following
#     notation: * - Residue conserved; ' ' - Residue not conserved.
#
sub _ReadFileAndSetupSequencesData {
  my($SequenceFile, $SequenceType) = @_;
  my($SequenceDataMapRef);

  $SequenceDataMapRef = undef;

  # Read sequence file...
  $SequenceDataMapRef = '';
  if ($SequenceType =~ /^ClustalW$/i) {
    $SequenceDataMapRef = _ReadClustalWFile($SequenceFile);
  }
  elsif ($SequenceType =~ /^Pearson$/i) {
    $SequenceDataMapRef = _ReadPearsonFastaFile($SequenceFile);
  }
  elsif ($SequenceType =~ /^PIR$/i) {
    $SequenceDataMapRef = _ReadPIRFastaFile($SequenceFile);
  }
  elsif ($SequenceType =~ /^MSF$/i) {
    $SequenceDataMapRef = _ReadMSFFile($SequenceFile);
  }
  else {
    return $SequenceDataMapRef;
  }

  if (exists $SequenceDataMapRef->{ConservedAnnotation}) {
    return ($SequenceDataMapRef);
  }
  if (!(($SequenceDataMapRef->{Count} > 1) && (AreSequenceLengthsIdentical($SequenceDataMapRef)))) {
    return ($SequenceDataMapRef);
  }

  # Use the first sequence to setup an empty ConservedAnnotation key...
  # And mark fully conserved residues...
  #
  my($ID, $Sequence, $FirstSequence, $FirstSeqLen, $Res, $FirstRes, $ResConserved, $Index);
  $ID = $SequenceDataMapRef->{IDs}[0];
  $FirstSequence = $SequenceDataMapRef->{Sequence}{$ID};
  $FirstSeqLen = length($FirstSequence);
  $SequenceDataMapRef->{ConservedAnnotation} = '';
  for $Index (0 .. ($FirstSeqLen - 1)) {
    $FirstRes = '';
    $ResConserved = 1;
    ID: for $ID (@{$SequenceDataMapRef->{IDs}}) {
      $Sequence = $SequenceDataMapRef->{Sequence}{$ID};
      $Res = substr($Sequence, $Index, 1);
      if (!$FirstRes) {
	$FirstRes = $Res;
	next ID;
      }
      if (($Res !~ /[A-Z]/i) || ($Res ne $FirstRes)) {
	$ResConserved = 0;
	last ID;
      }
    }
    if ($ResConserved) {
      $SequenceDataMapRef->{ConservedAnnotation} .= '*';
    }
    else {
      $SequenceDataMapRef->{ConservedAnnotation} .= ' ';
    }
  }

  return ($SequenceDataMapRef);
}

# Read sequence data in ClustalW multiple sequence alignment file and
# return a reference to hash with these keys and values:
#
# {IDs} - Array of sequence IDs
# {Count} - Number of sequences
# {Description}{$ID} - Sequence description
# {Sequence}{$ID} - Sequence for a specific ID
# {InputFileType} - Sequence file format
# {ConservedAnnotation} - Conserved residue annonations: space, *, : , .
#
#
#
# And based on ClustalW/X manual, here is what the ConservedAnnonations mean:
#
# '*' indicates positions which have a single, fully conserved residue
#
# ':' indicates that one of the following 'strong' groups is fully conserved: STA
#    NEQK NHQK NDEQ QHRK MILV MILF HY FYW

# '.' indicates that one of the following 'weaker' groups is fully conserved:
#     CSA ATV SAG STNK STPA SGND SNDEQK NDEQHK NEQHRK FVLIM HFY
#
# These are all the positively scoring groups that occur in the Gonnet Pam250
# matrix. The strong and weak groups are defined as strong score >0.5 and weak
# score =<0.5 respectively.
#
sub _ReadClustalWFile {
  my($SequenceFile) = @_;
  my(%SequencesDataMap);

  # Initialize data...
  %SequencesDataMap = ();
  @{$SequencesDataMap{IDs}} = ();
  %{$SequencesDataMap{Description}} = ();
  %{$SequencesDataMap{Sequence}} = ();
  $SequencesDataMap{Count} = 0;
  $SequencesDataMap{ConservedAnnotation} = '';
  $SequencesDataMap{InputFileType} = 'ClustalW';

  open SEQUENCEFILE, "$SequenceFile" or die "Couldn't open $SequenceFile: $!\n";

  my($Line, $LineLength, $AnnotationStart, $AnnotationLength, $Annotation, $Sequence, $SequenceLength, $ID, $IDIndex);

  # Ignore the header line...
  $Line = <SEQUENCEFILE>;

  LINE: while ($Line = GetTextLine(\*SEQUENCEFILE)) {
    if (($Line =~ /^[ \*\:\.]/) && ($Line !~ /[A-Z]/i)) {
      # Annotation for sequences: fully conserverd, weaker or stronger group conserverd.
      # Extract it and save...
      $LineLength = length($Line);
      $AnnotationStart = $LineLength - $SequenceLength;
      $AnnotationLength = $SequenceLength;
      $Annotation = substr($Line, $AnnotationStart, $AnnotationLength);
      $SequencesDataMap{ConservedAnnotation} .= $Annotation;
    }
    else {
      # Extract ID and sequences...
      ($ID, $Sequence)= $Line =~ /^[ ]*(.*?)[ ]+(.*?)[ 01-9]*$/;
      $Sequence =~ s/ //g;
      if (!($ID && $Sequence)) {
	next LINE;
      }

      if (exists $SequencesDataMap{Sequence}{$ID}) {
	# Append to existing alignment value...
	$SequenceLength = length($Sequence);
	$SequencesDataMap{Sequence}{$ID} .= $Sequence;
      }
      else {
	# New alignment data...
	$SequencesDataMap{Count} += 1;
	push @{$SequencesDataMap{IDs}}, $ID;
	$SequencesDataMap{Description}{$ID} = $ID;
	$SequencesDataMap{Sequence}{$ID} = $Sequence;
	$SequenceLength = length($Sequence);
      }
    }
  }
  close SEQUENCEFILE;
  return (\%SequencesDataMap);
}

# Read Pearson fasta file and return a reference to hash with these keys:
#
# {IDs} - Array of sequence IDs
# {Count} - Number of sequences
# {Description}{$ID} - Sequence description
# {Sequence}{$ID} - Sequence for a specific ID
# {InputFileType} - Sequence file format
# {ConservedAnnotation} - Conserved residue annonation
#
sub _ReadPearsonFastaFile {
  my($FastaFileName, $ID, $Description, $Line, $IgnoreID, @LineWords, %FastaDataMap);

  ($FastaFileName) = @_;

  %FastaDataMap = ();
  @{$FastaDataMap{IDs}} =();
  %{$FastaDataMap{Description}} =();
  %{$FastaDataMap{Sequence}} =();
  $FastaDataMap{Count} = 0;
  $FastaDataMap{InputFileType} = 'Pearson';

  open FASTAFILE, "$FastaFileName" or die "Couldn't open $FastaFileName: $!\n";
  $ID = '';
  $IgnoreID = 0;
  LINE: while ($Line = GetTextLine(\*FASTAFILE)) {
    if ($Line =~ /^\>/) {
      # Start of a new ID...
      $Line =~ s/^\>//;
      $Line = RemoveLeadingWhiteSpaces($Line);
      @LineWords = ();
      @LineWords = split / /, $Line;

      $ID = $LineWords[0];
      $ID =~ s/ //g;
      $Description = $Line;

      $IgnoreID = 0;
      if (exists $FastaDataMap{Sequence}{$ID}) {
	$IgnoreID = 1;
	warn "Warning: ID, $ID, in Fasta file already exists. Ignoring ID and sequence data...\n";
	next LINE;
      }
      push @{$FastaDataMap{IDs}}, $ID;
      $FastaDataMap{Description}{$ID} = $Description;
      $FastaDataMap{Count} += 1;
      next LINE;
    }
    if ($IgnoreID) { next LINE; }

    # Remove any spaces in the sequence...
    $Line =~ s/ //g;
    # Sequence data for active ID...
    if (exists $FastaDataMap{Sequence}{$ID}) {
      $FastaDataMap{Sequence}{$ID} .= $Line;
    }
    else {
      $FastaDataMap{Sequence}{$ID} = $Line;
    }
  }
  close FASTAFILE;
  return \%FastaDataMap;
}

# Read PIR fasta file and return a reference to hash with these keys:
#
# {IDs} - Array of sequence IDs
# {Count} - Number of sequences
# {Description}{$ID} - Sequence description
# {Sequence}{$ID} - Sequence for a specific ID
# {InputFileType} - Sequence file format
# {ConservedAnnotation} - Conserved residue annonation
#
# Format:
# A sequence in PIR format consists of:
# One line starting with
#   a ">" (greater-than) sign, followed by
#   a two-letter code describing the sequence type code (P1, F1, DL, DC, RL, RC, N3, N1 or XX), followed by
#   a semicolon, followed by
#   the sequence identification code (the database ID-code).
# One line containing a textual description of the sequence.
# One or more lines containing the sequence itself. The end of the
# sequence is marked by a "*" (asterisk) character.
#
# A file in PIR format may comprise more than one sequence.
#
# The PIR format is also often referred to as the NBRF format.
#
# Code SequenceType
# P1    Protein (complete)
# F1    Protein (fragment)
# DL    DNA (linear)
# DC    DNA (circular)
# RL    RNA (linear)
# RC   RNA (circular)
# N3    tRNA
# N1    Other functional RNA
#

sub _ReadPIRFastaFile {
  my($FastaFileName, $ID, $Description, $Line, $SequenceTypeCode, $ReadingSequenceData, %FastaDataMap);

  ($FastaFileName) = @_;

  %FastaDataMap = ();
  @{$FastaDataMap{IDs}} =();
  %{$FastaDataMap{Description}} =();
  %{$FastaDataMap{Sequence}} =();
  %{$FastaDataMap{SequenceTypeCode}} =();
  $FastaDataMap{Count} = 0;
  $FastaDataMap{InputFileType} = 'PIR';

  open FASTAFILE, "$FastaFileName" or die "Couldn't open $FastaFileName: $!\n";
  $ID = '';
  $ReadingSequenceData = 0;
  LINE: while ($Line = GetTextLine(\*FASTAFILE)) {
    if ($Line =~ /^\>/) {
      # Start of a new ID...
      $Line =~ s/^\>//;
      $Line = RemoveLeadingWhiteSpaces($Line);
      ($SequenceTypeCode, $ID) = /^\>(.*?)\;(.*?)$/;

      # Use next line to retrieve sequence description...
      $Line = GetTextLine(\*FASTAFILE);
      $Line = RemoveLeadingWhiteSpaces($Line);
      $Description = $Line;

      if (exists $FastaDataMap{Sequence}{$ID}) {
	warn "Warning: ID, $ID, in Fasta file already exists. Ignoring ID and sequence data...\n";
	next LINE;
      }
      $ReadingSequenceData = 1;
      push @{$FastaDataMap{IDs}}, $ID;
      $FastaDataMap{SequenceTypeCode}{$ID} = $SequenceTypeCode;
      $FastaDataMap{Description}{$ID} = $Description;
      $FastaDataMap{Count} += 1;
      next LINE;
    }
    if (!$ReadingSequenceData) { next LINE; }

    # Remove any spaces in the sequence...
    $Line =~ s/ //g;
    if ($Line =~ /[\*]$/) {
      # End of sequence...
      $ReadingSequenceData = 0;
      $Line =~ s/[\*]$//;
    }
    # Sequence data for active ID...
    if (exists $FastaDataMap{Sequence}{$ID}) {
      $FastaDataMap{Sequence}{$ID} .= $Line;
    }
    else {
      $FastaDataMap{Sequence}{$ID} = $Line;
    }
  }
  close FASTAFILE;
  return \%FastaDataMap;
}

# Read MSF file and return a reference to hash with these keys:
#
# {IDs} : Array of IDs in order as they appear in file
# {Count}: ID count...
# {Description}{$ID} : Description data...
# {Sequence}{$ID} : Sequence data...
#
sub _ReadMSFFile {
  my($MSFFileName, $Line, @LineWords, %MSFDataMap);

  ($MSFFileName) = @_;

  %MSFDataMap = ();
  @{$MSFDataMap{IDs}} =();
  %{$MSFDataMap{Description}} =();
  %{$MSFDataMap{Sequence}} =();
  $MSFDataMap{Count} = 0;
  $MSFDataMap{InputFileType} = 'MSF';

  open MSFFILE, "$MSFFileName" or die "Couldn't open $MSFFileName: $!\n";

  # Collect sequences and IDs...
  #
  # '//' after the name fields indicates end of header list and start of sequence data.
  #
  my($ID, $Len, $Check, $Weight, $Sequence, $NameFieldsFound, %MSFIDsMap);
  %MSFIDsMap = ();
  $NameFieldsFound = 0;
  LINE: while ($Line = GetTextLine(\*MSFFILE)) {
    if ($Line =~ /Name:/) {
      $NameFieldsFound++;
      ($ID, $Len, $Check, $Weight) = $Line =~ /^[ ]*Name:[ ]+(.*?)[ ]+Len:[ ]+(.*?)[ ]+Check:[ ]+(.*?)[ ]+Weight:[ ]+(.*?)[ ]*$/;
      if ($ID =~ / /) {
	($ID) = $ID =~ /^(.*?)[ ]+/
      }
      if (exists $MSFIDsMap{$ID}) {
	warn "Warning: ID, $ID, in MSF file already exists. Ignoring ID and sequence data...\n";
	next LINE;
      }
      $MSFIDsMap{$ID} = $ID;
      push @{$MSFDataMap{IDs}}, $ID;
      $MSFDataMap{Description}{$ID} = $ID;
      $MSFDataMap{Count} += 1;
    }
    elsif ( /\/\// && $NameFieldsFound) {
      # End of header list...
      last LINE;
    }
  }
  # Collect all sequences...
  #
  my($FirstField, $SecondField);
  while ($Line = GetTextLine(\*MSFFILE)) {
    ($FirstField, $SecondField) = $Line =~ /^[ ]*(.*?)[ ]+(.*?)$/;
    if (exists $MSFIDsMap{$FirstField}) {
      # It's ID and sequence data...
      $ID = $FirstField;
      $Sequence = $SecondField;
      # Take out spaces and leave the gap characters...
      $Sequence =~ s/ //g;
      if ($MSFDataMap{Sequence}{$ID}) {
	$MSFDataMap{Sequence}{$ID} .= $Sequence;
      }
      else {
	$MSFDataMap{Sequence}{$ID} = $Sequence;
      }
    }
  }

  close MSFFILE;
  return \%MSFDataMap;
}


1;

__END__

=head1 NAME

SequenceFileUtil

=head1 SYNOPSIS

use SequenceFileUtil ;

use SequenceFileUtil qw(:all);

=head1 DESCRIPTION

B<SequenceFileUtil> module provides the following functions:

AreSequenceLengthsIdentical, CalcuatePercentSequenceIdentity,
CalculatePercentSequenceIdentityMatrix, GetLongestSequence, GetSequenceLength,
GetShortestSequence, IsClustalWSequenceFile, IsGapResidue, IsMSFSequenceFile,
IsPIRFastaSequenceFile, IsPearsonFastaSequenceFile, IsSupportedSequenceFile,
ReadClustalWSequenceFile, ReadMSFSequenceFile, ReadPIRFastaSequenceFile,
ReadPearsonFastaSequenceFile, ReadSequenceFile, RemoveSequenceAlignmentGapColumns,
RemoveSequenceGaps, WritePearsonFastaSequenceFile
SequenceFileUtil module provides various methods to process sequence
files and retreive appropriate information.

=head1 FUNCTIONS

=over 4

=item B<AreSequenceLengthsIdentical>

    $Status = AreSequenceLengthsIdentical($SequencesDataRef);

Checks the lengths of all the sequences available in I<SequencesDataRef> and returns 1
or 0 based whether lengths of all the sequence is same.

=item B<CalcuatePercentSequenceIdentity>

    $PercentIdentity =
       AreSequenceLengthsIdenticalAreSequenceLengthsIdentical(
          $Sequence1, $Sequence2, [$IgnoreGaps, $Precision]);

Returns percent identity between I<Sequence1> and I<Sequence2>. Optional arguments
I<IgnoreGaps> and I<Precision> control handling of gaps in sequences and precision of the
returned value. By default, gaps are ignored and precision is set up to 1 decimal.

=item B<CalculatePercentSequenceIdentityMatrix>

    $IdentityMatrixDataRef = CalculatePercentSequenceIdentityMatrix(
                             $SequencesDataRef, [$IgnoreGaps,
                             $Precision]);

Calculate pairwise percent identity between all the sequences available in I<SequencesDataRef>
and returns a reference to identity matrix hash. Optional arguments I<IgnoreGaps> and
I<Precision> control handling of gaps in sequences and precision of the returned value. By default, gaps
are ignored and precision is set up to 1 decimal.

=item B<GetSequenceLength>

    $SeqquenceLength = GetSequenceLength($Sequence, [$IgnoreGaps]);

Returns length of the specified sequence. Optional argument I<IgnoreGaps> controls handling
of gaps. By default, gaps are ignored.

=item B<GetShortestSequence>

   ($ID, $Sequence, $SeqLen, $Description) = GetShortestSequence(
       $SequencesDataRef, [$IgnoreGaps]);

Checks the lengths of all the sequences available in $SequencesDataRef and returns $ID,
$Sequence, $SeqLen, and $Description values for the shortest sequence. Optional arguments $IgnoreGaps
controls handling of gaps in sequences. By default, gaps are ignored.

=item B<GetLongestSequence>

   ($ID, $Sequence, $SeqLen, $Description) = GetLongestSequence(
       $SequencesDataRef, [$IgnoreGaps]);

Checks the lengths of all the sequences available in I<SequencesDataRef> and returns B<ID>,
B<Sequence>, B<SeqLen>, and B<Description> values for the longest sequence. Optional argument
$I<IgnoreGaps> controls handling of gaps in sequences. By default, gaps are ignored.

=item B<IsGapResidue>

    $Status = AreSequenceLengthsIdentical($Residue);

Returns 1 or 0 based on whether I<Residue> corresponds to a gap. Any character other than A to Z is
considered a gap residue.

=item B<IsSupportedSequenceFile>

    $Status = IsSupportedSequenceFile($SequenceFile);

Returns 1 or 0 based on whether I<SequenceFile> corresponds to a supported sequence
format.

=item B<IsClustalWSequenceFile>

    $Status = IsClustalWSequenceFile($SequenceFile);

Returns 1 or 0 based on whether I<SequenceFile> corresponds to Clustal sequence alignment
format.

=item B<IsPearsonFastaSequenceFile>

    $Status = IsPearsonFastaSequenceFile($SequenceFile);

Returns 1 or 0 based on whether I<SequenceFile> corresponds to Pearson FASTA sequence
format.

=item B<IsPIRFastaSequenceFile>

    $Status = IsPIRFastaSequenceFile($SequenceFile);

Returns 1 or 0 based on whether I<SequenceFile> corresponds to PIR FASTA sequence
format.

=item B<IsMSFSequenceFile>

    $Status = IsClustalWSequenceFile($SequenceFile);

Returns 1 or 0 based on whether I<SequenceFile> corresponds to MSF sequence alignment
format.

=item B<ReadSequenceFile>

    $SequenceDataMapRef = ReadSequenceFile($SequenceFile);

Reads I<SequenceFile> and returns reference to a hash containing following key/value
pairs:

    $SequenceDataMapRef->{IDs} - Array of sequence IDs
    $SequenceDataMapRef->{Count} - Number of sequences
    $SequenceDataMapRef->{Description}{$ID} - Sequence description
    $SequenceDataMapRef->{Sequence}{$ID} - Sequence for a specific ID
    $SequenceDataMapRef->{Sequence}{InputFileType} - File format

=item B<ReadClustalWSequenceFile>

    $SequenceDataMapRef = ReadClustalWSequenceFile($SequenceFile);

Reads ClustalW I<SequenceFile> and returns reference to a hash containing following key/value
pairs as describes in B<ReadSequenceFile> method.

=item B<ReadMSFSequenceFile>

    $SequenceDataMapRef = ReadMSFSequenceFile($SequenceFile);

Reads MSF I<SequenceFile> and returns reference to a hash containing following key/value
pairs as describes in B<ReadSequenceFile> method.

=item B<ReadPIRFastaSequenceFile>

    $SequenceDataMapRef = ReadPIRFastaSequenceFile($SequenceFile);

Reads PIR FASTA I<SequenceFile> and returns reference to a hash containing following key/value
pairs as describes in B<ReadSequenceFile> method.

=item B<ReadPearsonFastaSequenceFile>

    $SequenceDataMapRef = ReadPearsonFastaSequenceFile($SequenceFile);

Reads Pearson FASTA I<SequenceFile> and returns reference to a hash containing following key/value
pairs as describes in B<ReadSequenceFile> method.

=item B<RemoveSequenceGaps>

    $SeqWithoutGaps = RemoveSequenceGaps($Sequence);

Removes gaps from I<Sequence> and return a sequence without any gaps.

=item B<RemoveSequenceAlignmentGapColumns>

    $NewAlignmentDataMapRef = RemoveSequenceAlignmentGapColumns(
                              $AlignmentDataMapRef);

Using input alignment data map ref containing following keys, generate a new hash with
same set of keys after residue columns containg only gaps have been removed:

    {IDs} : Array of IDs in order as they appear in file
    {Count}: ID count
    {Description}{$ID} : Description data
    {Sequence}{$ID} : Sequence data

=item B<WritePearsonFastaSequenceFile>

    WritePearsonFastaSequenceFile($SequenceFileName, $SequenceDataRef,
                                  [$MaxLength]);

Using sequence data specified via I<SequenceDataRef>, write out a Pearson FASTA sequence
file. Optional argument I<MaxLength> controls maximum length sequence in each line; default is
80.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

PDBFileUtil.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
