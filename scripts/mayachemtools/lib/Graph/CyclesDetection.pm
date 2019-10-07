package Graph::CyclesDetection;
#
# File: CyclesDetection.pm
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
use Carp;
use Exporter;
use Graph;
use Graph::Path;
use Graph::PathGraph;
use BitVector;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyCyclesDetection';

# Class constructor...
sub new {
  my($Class, $Graph) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializeCyclesDetection($Graph);

  return $This;
}

# Initialize object data...
sub _InitializeCyclesDetection {
  my($This, $Graph) = @_;

  $This->{Graph} = $Graph;

  # Cycles list...
  @{$This->{AllCyclicPaths}} = ();

  # Cyclic paths which are not part of any other cycle...
  @{$This->{IndependentCyclicPaths}} = ();

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Detect all cycles in graph...
#
sub DetectCycles {
  my($This) = @_;

  return $This->DetectCyclesUsingCollapsingPathGraphMethodology();
}

# Detect all cycles in the graph using collapsing path graph [Ref 31] methodology...
#
# Note:
#   . For topologically complex graphs containing large number of cycles,
#     CollapseVertexAndCollectCyclicPathsDetectCycles method implemented in
#     PathGraph can time out in which case no cycles are detected or
#     assigned.
#
sub DetectCyclesUsingCollapsingPathGraphMethodology {
  my($This) = @_;
  my($PathGraph);

  # Create a path graph...
  $PathGraph = new Graph::PathGraph($This->{Graph});

  # For a vertex to be in a cycle, its degree must be >=2. So delete vertices recursively
  # till all vertices with degree less than 2 have been deleted...
  $PathGraph->DeleteVerticesWithDegreeLessThan(2);

  # Setup a VertexID and EdgeID to index map usage during retrieval of independent cycles to
  # avoid going over all vertices in all cylic paths later...
  #
  my($VertexIDsToIndicesRef, $LargestVertexIndex, $EdgeIDsToIndicesRef, $LargestEdgeIDIndex);
  ($VertexIDsToIndicesRef, $LargestVertexIndex) = $This->_SetupVertexIDsToIndicesMap($PathGraph);
  ($EdgeIDsToIndicesRef, $LargestEdgeIDIndex) = $This->_SetupEdgeIDsToIndicesMap($PathGraph);

  # Recursively collapse vertices with lowest degree...
  my($VertexID, $CycleVertexID);
  while ($VertexID = $PathGraph->GetVertexWithSmallestDegree()) {
      if (!$PathGraph->CollapseVertexAndCollectCyclicPaths($VertexID)) {
	# Cycles detection didn't finish...
	return undef;
      }
  }

  # Get detected cycles and save 'em sorted by size...
  push @{$This->{AllCyclicPaths}}, sort { $a->GetLength() <=> $b->GetLength() } $PathGraph->GetCyclicPaths();

  # Retrieve independent cyclic paths...
  return $This->_RetrieveIndependentCycles($VertexIDsToIndicesRef, $LargestVertexIndex, $EdgeIDsToIndicesRef, $LargestEdgeIDIndex);
}

# Retrieve and save independent cyclic paths...
#
# Set of independent cycles identified using this method doesn't correspond to basis set of
# rings or smallest set of smallest rings (SSSR) [ Refs 29-30 ]; instead, set of cycles identified
# as independent cycles simply correspond to cycles which contain no other cycle as their
# subcycles and can't be described as linear combination of smaller cycles. And it also happen
# to contain all the rings in basis set of rings and SSSR. In other words, it's a superset of basis set
# of cycles and SSSR. For example, six four membered cycles are identified for cubane which is one
# more than the basis set of cycles.
#
sub _RetrieveIndependentCycles {
  my($This, $VertexIDsToIndicesRef, $LargestVertexIndex, $EdgeIDsToIndicesRef, $LargestEdgeIDIndex) = @_;

  # Only 1 or 0 cyclic paths...
  if (@{$This->{AllCyclicPaths}} <= 1) {
    push @{$This->{IndependentCyclicPaths}}, @{$This->{AllCyclicPaths}};
    return $This;
  }

  # Setup bit vectors for each cyclic path with size of each bit vector corresponding to
  # maximum vertex index plus one...
  my($CyclicPath, $BitVector, @BitNums, @CyclicPathBitVectors, @CyclicPathEdgeIDsBitVectors);

  @CyclicPathBitVectors = (); @CyclicPathEdgeIDsBitVectors = ();

  # Set bits corresponding to VertexIDs EdgeIDs using their indices...
  for $CyclicPath (@{$This->{AllCyclicPaths}}) {
    $BitVector = new BitVector($LargestVertexIndex);
    @BitNums = map { $VertexIDsToIndicesRef->{$_} } $CyclicPath->GetVertices();
    $BitVector->SetBits(@BitNums);
    push @CyclicPathBitVectors, $BitVector;

    $BitVector = new BitVector($LargestEdgeIDIndex);
    @BitNums = map { $EdgeIDsToIndicesRef->{$_} } $This->_GetPathEdgeIDs($CyclicPath);
    $BitVector->SetBits(@BitNums);
    push @CyclicPathEdgeIDsBitVectors, $BitVector;
  }

  # First smallest cyclic path always ends up as an independent cyclic path...
  push @{$This->{IndependentCyclicPaths}}, $This->{AllCyclicPaths}[0];

  # Retrieve other independent cyclic paths...
  my($CurrentIndex, $PreviousIndex, $CurrentCyclicPath, $PreviousCyclicPath, $CurrentPathLength, $PreviousPathLength, $CurrentBitVector, $PreviousBitVector, $CurrentAndPreviousBitVectpor, $AllPreviousSmallerPathsBitVector, $AllPreviousSmallerPathsEdgeIDsBitVector, $CurrentEdgeIDsBitVector, $AndBitVector, %SmallerPathAlreadyAdded, %SkipPath);

  %SkipPath = ();
  %SmallerPathAlreadyAdded = ();
  $AllPreviousSmallerPathsBitVector = new BitVector($LargestVertexIndex);
  $AllPreviousSmallerPathsEdgeIDsBitVector = new BitVector($LargestEdgeIDIndex);

  CURRENT: for $CurrentIndex (1 .. $#{$This->{AllCyclicPaths}}) {
    if (exists $SkipPath{$CurrentIndex}) {
      next CURRENT;
    }
    $CurrentCyclicPath = $This->{AllCyclicPaths}[$CurrentIndex];
    $CurrentBitVector = $CyclicPathBitVectors[$CurrentIndex];
    $CurrentPathLength = $CurrentCyclicPath->GetLength();

    PREVIOUS: for $PreviousIndex (0 .. ($CurrentIndex - 1)) {
      if (exists $SkipPath{$PreviousIndex}) {
	next PREVIOUS;
      }
      $PreviousCyclicPath = $This->{AllCyclicPaths}[$PreviousIndex];
      $PreviousBitVector = $CyclicPathBitVectors[$PreviousIndex];

      # Is previous path a subset of current path?
      $CurrentAndPreviousBitVectpor = $PreviousBitVector &  $CurrentBitVector;
      if ($PreviousBitVector->GetNumOfSetBits() == $CurrentAndPreviousBitVectpor->GetNumOfSetBits()) {
	# Current path doesn't qualify an independent path...
	$SkipPath{$CurrentIndex} = 1;
	next CURRENT;
      }

      $PreviousPathLength = $PreviousCyclicPath->GetLength();
      if ($PreviousPathLength < $CurrentPathLength) {
	# Mark cycle vertices seen in cyclic paths with length smaller than current path...
	if (! exists $SmallerPathAlreadyAdded{$PreviousIndex}) {
	  $SmallerPathAlreadyAdded{$PreviousIndex} = 1;
	  $AllPreviousSmallerPathsBitVector = $AllPreviousSmallerPathsBitVector | $PreviousBitVector;
	  $AllPreviousSmallerPathsEdgeIDsBitVector = $AllPreviousSmallerPathsEdgeIDsBitVector | $CyclicPathEdgeIDsBitVectors[$PreviousIndex];
	}
      }
    }
    if ($AllPreviousSmallerPathsBitVector->GetNumOfSetBits()) {
      # Is current path a linear combination of smaller paths?
      $AndBitVector = $AllPreviousSmallerPathsBitVector &  $CurrentBitVector;
      if ($CurrentBitVector->GetNumOfSetBits() == $AndBitVector->GetNumOfSetBits()) {
	# Are all edges in the current path already present in smaller paths?
	$CurrentEdgeIDsBitVector = $CyclicPathEdgeIDsBitVectors[$CurrentIndex];
	$AndBitVector = $AllPreviousSmallerPathsEdgeIDsBitVector &  $CurrentEdgeIDsBitVector;

	if ($CurrentEdgeIDsBitVector->GetNumOfSetBits() == $AndBitVector->GetNumOfSetBits()) {
	  # Current path doesn't qualify an independent path...
	  $SkipPath{$CurrentIndex} = 1;
	  next CURRENT;
	}
      }
    }
    # It's an independent cyclic path...
    push @{$This->{IndependentCyclicPaths}}, $CurrentCyclicPath;
  }
  return $This;
}

# Setup a VertexID to index map...
#
sub _SetupVertexIDsToIndicesMap {
  my($This, $PathGraph) = @_;
  my($LargestVertexIndex, @VertexIDs, %VertexIDsMap);

  %VertexIDsMap = (); @VertexIDs = (); $LargestVertexIndex = 0;

  @VertexIDs = $PathGraph->GetVertices();
  if (! @VertexIDs) {
    return (\%VertexIDsMap, $LargestVertexIndex);
  }
  @VertexIDsMap{ @VertexIDs } = (0 .. $#VertexIDs);
  $LargestVertexIndex = scalar @VertexIDs;

  return (\%VertexIDsMap, $LargestVertexIndex);
}

# Setup a Edge to index map using paths associated to egdes in an intial
# path graph...
#
sub _SetupEdgeIDsToIndicesMap {
  my($This, $PathGraph) = @_;
  my($Path, $LargestEdgeIndex, @EdgeIDs, %EdgeIDsMap);

  %EdgeIDsMap = (); @EdgeIDs = (); $LargestEdgeIndex = 0;

  @EdgeIDs = ();
  for $Path ($PathGraph->GetPaths()) {
    push @EdgeIDs, $This->_GetPathEdgeIDs($Path);
  }

  if (! @EdgeIDs) {
    return (\%EdgeIDsMap, $LargestEdgeIndex);
  }

  @EdgeIDsMap{ @EdgeIDs } = (0 .. $#EdgeIDs);
  $LargestEdgeIndex = scalar @EdgeIDs;

  return (\%EdgeIDsMap, $LargestEdgeIndex);
}

# Get path edge IDs or number of edges. The edge IDs are generated from
# edge vertices and correpond to VertexID1-VertexID2 where VertexID1 <=
# VertexID2.
#
sub _GetPathEdgeIDs {
  my($This, $Path) = @_;
  my(@EdgesVertexIDs, @EdgeIDs);

  @EdgesVertexIDs = (); @EdgeIDs = ();
  @EdgesVertexIDs = $Path->GetEdges();
  if (!@EdgesVertexIDs) {
    return wantarray ? @EdgeIDs : (scalar @EdgeIDs);
  }

  # Set up edge IDs...
  my($Index, $VertexID1, $VertexID2, $EdgeID);

  for ($Index = 0; $Index < $#EdgesVertexIDs; $Index += 2) {
    $VertexID1 = $EdgesVertexIDs[$Index]; $VertexID2 = $EdgesVertexIDs[$Index + 1];
    $EdgeID = ($VertexID1 <= $VertexID2) ? "$VertexID1-$VertexID2" : "$VertexID2-$VertexID1";
    push @EdgeIDs, $EdgeID;
  }

  return wantarray ? @EdgeIDs : (scalar @EdgeIDs);
}

# Return an array containing references to cyclic paths or number of cylic paths...
#
sub GetAllCyclicPaths {
  my($This) = @_;

  return wantarray ? @{$This->{AllCyclicPaths}} : scalar @{$This->{AllCyclicPaths}};
}

# Get cyclic paths which are independent. In otherwords, cycles which don't contain any other
# cycle as their subset...
#
sub GetIndependentCyclicPaths {
  my($This) = @_;

  return wantarray ? @{$This->{IndependentCyclicPaths}} : scalar @{$This->{IndependentCyclicPaths}};
}

# Return a string containg data for CyclesDetection object...
sub StringifyCyclesDetection {
  my($This) = @_;
  my($CyclesDetectionString, $CyclesCount, $CyclicPath);

  $CyclesCount = @{$This->{AllCyclicPaths}};
  $CyclesDetectionString = "AllCycles: Count - $CyclesCount";
  for $CyclicPath (@{$This->{AllCyclicPaths}}) {
    $CyclesDetectionString .= "; Cycle: " . join('-', $CyclicPath->GetVertices());
  }

  $CyclesCount = @{$This->{IndependentCyclicPaths}};
  $CyclesDetectionString .= "\nIndependentCycles: Count - $CyclesCount";
  for $CyclicPath (@{$This->{IndependentCyclicPaths}}) {
    $CyclesDetectionString .= "; Cycle: " . join('-', $CyclicPath->GetVertices());
  }

  return $CyclesDetectionString;
}

# Return a reference to new cycle detection object...
sub Copy {
  my($This) = @_;
  my($NewCyclesDetection);

  $NewCyclesDetection = Storable::dclone($This);

  return $NewCyclesDetection;
}

1;

__END__

=head1 NAME

CyclesDetection

=head1 SYNOPSIS

use Graph::CyclesDetection;

use Graph::CyclesDetection qw(:all);

=head1 DESCRIPTION

B<CyclesDetection> class provides the following methods:

new, Copy, DetectCycles, DetectCyclesUsingCollapsingPathGraphMethodology,
GetAllCyclicPaths, GetIndependentCyclicPaths, StringifyCyclesDetection

Cycles in a B<Graph> are detected using collapsing path graph [Ref 31]
methodology.

=head2 METHODS

=over 4

=item B<new>

    $NewCyclesDetection = new Graph::CyclesDetection($Graph);

Using specified I<Graph>, B<new> method creates a new B<CyclesDetection> object and returns
newly created B<CyclesDetection> object.

=item B<Copy>

    $NewCyclesDetection = $CyclesDetection->Copy();

Copies I<CyclesDetection> and its associated data using B<Storable::dclone> and returns a new
B<CyclesDetection> object.

=item B<DetectCycles>

    $CyclesDetection->DetectCycles();

Detects all cycles in a graph and returns I<CyclesDetection>.

=item B<DetectCyclesUsingCollapsingPathGraphMethodology>

    $CyclesDetection->DetectCyclesUsingCollapsingPathGraphMethodology();

Detects all cycles in a graph using collapsing path graph [Ref 31] methodology
and returns I<CyclesDetection>.

=item B<GetAllCyclicPaths>

    @AllCyclicPaths = $CyclesDetection->GetAllCyclicPaths();
    $NumOfAllCyclicPaths = $CyclesDetection->GetAllCyclicPaths();

Returns an array containing references to all cyclic paths identified during cycles
detection. In scalar text, number of cycles is returned.

=item B<GetIndependentCyclicPaths>

    @IndependentCyclicPaths = $CyclesDetection->GetAllCyclicPaths();
    $NumOfIndependentCyclicPaths = $CyclesDetection->GetAllCyclicPaths();

Returns an array containing references to independent cyclic paths identified during cycles
detection. In scalar text, number of cycles is returned.

A set of independent cycles identified during cycles detection doesn't correspond to the basis set of
rings or smallest set of smallest rings (SSSR) [ Refs 29-30 ]; instead, set of cycles indentified
as independent cycles simply correpond to cycles which contain no other cycle as their
subcycles and can't be described as a linear combination of smaller cycles. And it also happens
to contain all the rings in basis set of rings and SSSR. In other words, it's a superset of a basis set
of cycles and SSSR. For example, six four membered cycles are indentified for cubane, which is one
more than the basis set of cycles.

=item B<StringifyCyclesDetection>

    $String = $CyclesDetection->StringifyCyclesDetection();

Returns a string containing information about I<CyclesDetection> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Graph.pm, Path.pm, PathGraph.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
