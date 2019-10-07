package Graph::PathsTraversal;
#
# File: PathsTraversal.pm
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

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyPathsTraversal';

# Class constructor...
sub new {
  my($Class, $Graph) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializePathsTraversal($Graph);

  return $This;
}

# Initialize object data...
sub _InitializePathsTraversal {
  my($This, $Graph) = @_;

  # Graph object...
  $This->{Graph} = $Graph;

  # Traversal mode: Vertex or Path
  $This->{TraversalMode} = '';

  # Traversal type: DFS, DFSWithLimit, BFS, BFSWithLimit...
  $This->{TraversalType} = '';

  # For finding root vertex and controlling search...
  my(@VertexIDs);
  @VertexIDs = $This->{Graph}->GetVertices();
  %{$This->{VerticesToVisit}} = ();
  @{$This->{VerticesToVisit}}{ @VertexIDs } = @VertexIDs;

  # Root vertex of all visited vertices...
  %{$This->{VerticesRoots}} = ();

  # Visited vertices...
  %{$This->{VisitedVertices}} = ();

  # Finished vertices...
  %{$This->{FinishedVertices}} = ();

  # List of active vertices during DFS/BFS search...
  @{$This->{ActiveVertices}} = ();

  # List of ordered vertices traversed during DFS/BFS search...
  @{$This->{Vertices}} = ();

  # Vertex neighbors during traversal...
  %{$This->{VerticesNeighbors}} = ();

  # Vertices depth from root...
  %{$This->{VerticesDepth}} = ();

  # Predecessor of each vertex during vertex traversal. For root vertex, it's root itself...
  %{$This->{VerticesPredecessors}} = ();

  # Successors of each vertex during vertex traversal...
  %{$This->{VerticesSuccessors}} = ();

  # Vertices at different neighborhood levels during vertex traversal...
  @{$This->{VerticesNeighborhoods}} = ();

  # Vertices, along with their successors, at different neighborhood levels during vertex traversal...
  @{$This->{VerticesNeighborhoodsWithSuccessors}} = ();

  # Visited edges during Path TraversalMode...
  %{$This->{VisitedEdges}} = ();
  %{$This->{VisitedEdges}->{From}} = ();
  %{$This->{VisitedEdges}->{To}} = ();

  # Vertex path during Path TraversalMode...
  %{$This->{VerticesPaths}} = ();

  # Allow cycles in paths during VertexNeighborhood TraversalMode. By default, cycles are not allowed
  # during vertex traversal: a vertex is only visited once during BFS search for neighborhoods. For
  # neighborhood vertices search during successors identification, vertex cycles are explicity allowed
  # to indentify all successors.
  $This->{AllowVertexCycles} = 0;

  # Allow cycles in paths during Path TraversalMode...
  $This->{AllowPathCycles} = 1;

  # Cycle closure vertices during Path TraversalMode...
  %{$This->{CycleClosureVertices}} = ();

  # Paths traversed during search...
  @{$This->{Paths}} = ();

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Perform a depth first search (DFS)...
#
sub PerformDepthFirstSearch {
  my($This, $RootVertexID) = @_;

  if (defined $RootVertexID) {
    if (!$This->{Graph}->HasVertex($RootVertexID)) {
      carp "Warning: ${ClassName}->PerformDepthFirstSearch: No search performed: Vertex $RootVertexID doesn't exist...";
      return undef;
    }
  }
  return $This->_PerformVertexSearch("DFS", $RootVertexID);
}

# Perform a depth first search (DFS) with limit on depth...
#
sub PerformDepthFirstSearchWithLimit {
  my($This, $DepthLimit, $RootVertexID) = @_;

  if (!defined $DepthLimit) {
      carp "Warning: ${ClassName}->PerformDepthFirstSearchWithLimit: No search performed: Depth limit must be specified...";
      return undef;
  }
  if ($DepthLimit < 0) {
    carp "Warning: ${ClassName}->PerformDepthFirstSearchWithLimit: No search performed: Specified depth limit, $DepthLimit, must be a positive integer...";
    return undef;
  }
  if (defined $RootVertexID) {
    if (!$This->{Graph}->HasVertex($RootVertexID)) {
      carp "Warning: ${ClassName}->PerformDepthFirstSearchWithLimit: No search performed: Vertex $RootVertexID doesn't exist...";
      return undef;
    }
  }
  return $This->_PerformVertexSearch("DFSWithLimit", $RootVertexID, $DepthLimit);
}

# Perform a breadth first search (BFS)...
#
sub PerformBreadthFirstSearch {
  my($This, $RootVertexID) = @_;

  if (defined $RootVertexID) {
    if (!$This->{Graph}->HasVertex($RootVertexID)) {
      carp "Warning: ${ClassName}->PerformBreadthFirstSearch: No search performed: Vertex $RootVertexID doesn't exist...";
      return undef;
    }
  }
  return $This->_PerformVertexSearch("BFS", $RootVertexID);
}

# Perform a breadth first search (BFS) with limit...
#
sub PerformBreadthFirstSearchWithLimit {
  my($This, $DepthLimit, $RootVertexID) = @_;

  if (!defined $DepthLimit) {
      carp "Warning: ${ClassName}->PerformBreadthFirstSearchWithLimit: No search performed: Depth limit must be specified...";
      return undef;
  }
  if ($DepthLimit < 0) {
    carp "Warning: ${ClassName}->PerformBreadthFirstSearchWithLimit: No search performed: Specified depth limit, $DepthLimit, must be a positive integer...";
    return undef;
  }
  if (defined $RootVertexID) {
    if (!$This->{Graph}->HasVertex($RootVertexID)) {
      carp "Warning: ${ClassName}->PerformDepthFirstSearchWithLimit: No search performed: Vertex $RootVertexID doesn't exist...";
      return undef;
    }
  }
  return $This->_PerformVertexSearch("BFSWithLimit", $RootVertexID, $DepthLimit);
}

# Perform appropriate vertex search...
#
sub _PerformVertexSearch {
  my($This, $SearchType, $RootVertexID, $DepthLimit, $TargetVertexID) = @_;

  # Setup search...
  $This->{TraversalMode} = 'Vertex';
  $This->{TraversalType} = $SearchType;

  if (defined $RootVertexID) {
    $This->{RootVertex} = $RootVertexID;
  }
  if (defined $DepthLimit) {
    $This->{DepthLimit} = $DepthLimit;
  }
  if (defined $TargetVertexID) {
    $This->{TargetVertex} = $TargetVertexID;
  }

  # Perform search...
  return $This->_TraverseGraph();
}

# Perform DFS or BFS traversal with or without any limits...
#
sub _TraverseGraph {
  my($This) = @_;
  my($ProcessingVertices, $CurrentVertexID, $NeighborVertexID, $VertexID);

  if ($This->{TraversalMode} !~ /^(Vertex|Path|VertexNeighborhood)$/i) {
    return $This;
  }

  $ProcessingVertices = 1;

  VERTICES: while ($ProcessingVertices) {
    # Set root vertex...
    if (!@{$This->{ActiveVertices}}) {
      my($RootVertexID);

      $RootVertexID = $This->_GetRootVertex();
      if (!defined $RootVertexID) {
	$ProcessingVertices = 0; next VERTICES;
      }
      $This->_ProcessVisitedVertex($RootVertexID, $RootVertexID);
    }

    # Get current active vertex...
    $CurrentVertexID = $This->_GetActiveVertex();
    if (!defined $CurrentVertexID) {
      $ProcessingVertices = 0; next VERTICES;
    }

    # Get next available neighbor of current vertex...
    #
    $NeighborVertexID = $This->_GetNeighborVertex($CurrentVertexID);

    # Process neighbor or current vertex...
    if (defined $NeighborVertexID) {
      $This->_ProcessVisitedVertex($NeighborVertexID, $CurrentVertexID);
    }
    else {
      # Finished with all neighbors for current vertex...
      $This->_ProcessFinishedVertex($CurrentVertexID);
    }
  }
  return $This;
}

# Get root vertex to start the search...
#
# Notes:
#   . User specification of root vertex forces traversal in a specific connected component
#     of graph; To traverse find all connected components, perform traversal without specification of
#     a root vertex.
#
sub _GetRootVertex {
  my($This) = @_;
  my($RootVertexID);

  # Check for specified root vertex and constrain traversal to specific connected
  # component by setting root limit...
  if (exists $This->{RootVertex}) {
    $RootVertexID = $This->{RootVertex};
    delete $This->{RootVertex};
    $This->{RootVertexSpecified} = 1;

    return $RootVertexID;
  }
  # Traversal limited to connected component containing specified root vertex...
  if (exists $This->{RootVertexSpecified}) {
    return undef;
  }

  # Use first vertex in sorted available vertices list to get root vertex. Vertex
  # with largest degree could also be used as root vertex. However, for all
  # practical purposes, any arbitrary vertex can be used as root vertex to
  # start search for another disconnected component of the graph.
  #
  my(@VerticesToVisit);

  $RootVertexID = undef; @VerticesToVisit = ();
  @VerticesToVisit = sort { $a <=> $b } keys %{$This->{VerticesToVisit}};
  if (@VerticesToVisit) {
    $RootVertexID = $VerticesToVisit[0];
  }
  return $RootVertexID;
}

# Get current or new active vertex for DFS/BFS traversals...
#
sub _GetActiveVertex {
  my($This) = @_;
  my($ActiveVertexID);

  $ActiveVertexID = undef;
  if ($This->{TraversalType} =~ /^(DFS|DFSWithLimit)$/i) {
    # For DFS, it's last vertex in  LIFO queue...
    $ActiveVertexID = $This->{ActiveVertices}[-1];
  }
  elsif ($This->{TraversalType} =~ /^(BFS|BFSWithLimit)$/i) {
    # For BFS, it's first vertex in FIFO queue...
    $ActiveVertexID = $This->{ActiveVertices}[0];
  }
  return $ActiveVertexID;
}

# Get available neigbor of specified vertex...
#
sub _GetNeighborVertex {
  my($This, $VertexID) = @_;

  # Retrieve neighbors for vertex...
  if (!exists $This->{VerticesNeighbors}{$VertexID}) {
    @{$This->{VerticesNeighbors}{$VertexID}}  = ();

    if (exists $This->{DepthLimit}) {
      # Only collect neighbors to visit below specified depth limit...
      if ($This->{VerticesDepth}{$VertexID} < $This->{DepthLimit}) {
	push @{$This->{VerticesNeighbors}{$VertexID}}, $This->{Graph}->GetNeighbors($VertexID);
      }
      else {
	if (!exists $This->{RootVertexSpecified}) {
	  # Mark all other downstream neighbor vertices to be ignored from any further
	  # processing and avoid selection of a new root...
	  $This->_IgnoreDownstreamNeighbors($VertexID);
	}
      }
    }
    elsif (exists $This->{TargetVertex}) {
      if ($VertexID != $This->{TargetVertex}) {
	push @{$This->{VerticesNeighbors}{$VertexID}}, $This->{Graph}->GetNeighbors($VertexID);
      }
    }
    else {
      push @{$This->{VerticesNeighbors}{$VertexID}}, $This->{Graph}->GetNeighbors($VertexID);
    }
  }

  if ($This->{TraversalMode} =~ /^Path$/i) {
    # Get available neighbor for path search...
    return $This->_GetNeighborVertexDuringPathTraversal($VertexID);
  }
  elsif ($This->{TraversalMode} =~ /^Vertex$/i) {
    # Get unvisited neighbor for vertex search...
    return $This->_GetNeighborVertexDuringVertexTraversal($VertexID);
  }
  elsif ($This->{TraversalMode} =~ /^VertexNeighborhood$/i) {
    # Get available neighbor during vertex neighborhood search...
    return $This->_GetNeighborVertexDuringVertexNeighborhoodTraversal($VertexID);
  }
  return undef;
}

# Get unvisited neighbor of specified vertex during vertex traversal...
#
sub _GetNeighborVertexDuringVertexTraversal {
  my($This, $VertexID) = @_;
  my($NeighborVertexID, $UnvisitedNeighborVertexID);

  # Get unvisited neighbor...
  $UnvisitedNeighborVertexID = undef;
  NEIGHBOR: for $NeighborVertexID (@{$This->{VerticesNeighbors}{$VertexID}}) {
    if (!exists $This->{VisitedVertices}{$NeighborVertexID}) {
      $UnvisitedNeighborVertexID = $NeighborVertexID;
      last NEIGHBOR;
    }
  }
  return $UnvisitedNeighborVertexID;
}

# Get available neighbor of specified vertex during vertex neighborhood traversal...
#
sub _GetNeighborVertexDuringVertexNeighborhoodTraversal {
  my($This, $VertexID) = @_;
  my($NeighborVertexID, $UnvisitedNeighborVertexID);

  # Get available neighbor...
  $UnvisitedNeighborVertexID = undef;
  NEIGHBOR: for $NeighborVertexID (@{$This->{VerticesNeighbors}{$VertexID}}) {
    if (!exists $This->{VisitedVertices}{$NeighborVertexID}) {
      $UnvisitedNeighborVertexID = $NeighborVertexID;
      last NEIGHBOR;
    }
    # Look for any unvisited edge back to visited vertex...
    if ($This->_IsVisitedEdge($VertexID, $NeighborVertexID) || $This->_IsVisitedEdge($NeighborVertexID, $VertexID)) {
      next NEIGHBOR;
    }
    # Check its depth...
    if (exists $This->{DepthLimit}) {
      if (($This->{VerticesDepth}{$VertexID} + 1) > $This->{DepthLimit}) {
    	next NEIGHBOR;
      }
    }
    # Its an edge that makes a cycle during BFS search...
    if ($This->{AllowVertexCycles}) {
      $This->{CycleClosureVertices}{$NeighborVertexID} = 1;
      $UnvisitedNeighborVertexID = $NeighborVertexID;
      last NEIGHBOR;
    }
  }
  return $UnvisitedNeighborVertexID;
}

# Get available neighbor of specified vertex during path traversal...
#
sub _GetNeighborVertexDuringPathTraversal {
  my($This, $VertexID) = @_;
  my($NeighborVertexID, $UnvisitedNeighborVertexID);

  # Get unvisited neighbor...
  $UnvisitedNeighborVertexID = undef;
  NEIGHBOR: for $NeighborVertexID (@{$This->{VerticesNeighbors}{$VertexID}}) {
    if (!exists $This->{VisitedVertices}{$NeighborVertexID}) {
      # An unvisited vertex...
      $UnvisitedNeighborVertexID = $NeighborVertexID;
      last NEIGHBOR;
    }
    # Look for any unvisited edge back to visited vertex...
    if ($This->_IsVisitedEdge($VertexID, $NeighborVertexID) || $This->_IsVisitedEdge($NeighborVertexID, $VertexID)) {
      next NEIGHBOR;
    }
    # Check its depth...
    if (exists $This->{DepthLimit}) {
      if (($This->{VerticesDepth}{$VertexID} + 1) >= $This->{DepthLimit}) {
    	next NEIGHBOR;
      }
    }

    # It's the edge final edge of a cycle in case $NeighborVertexID is already in the path; otherwise, it's
    # part of the path from a different direction in a cycle or a left over vertex during Limit search.
    #
    if ($This->_IsCycleClosureEdge($VertexID, $NeighborVertexID)) {
      if ($This->{AllowPathCycles}) {
	$This->{CycleClosureVertices}{$NeighborVertexID} = 1;
	$UnvisitedNeighborVertexID = $NeighborVertexID;
	last NEIGHBOR;
      }
    }
    else {
      $UnvisitedNeighborVertexID = $NeighborVertexID;
      last NEIGHBOR;
    }
  }
  return $UnvisitedNeighborVertexID;
}

# Process visited vertex...
#
sub _ProcessVisitedVertex {
  my($This, $VertexID, $PredecessorVertexID) = @_;

  if (!exists $This->{VisitedVertices}{$VertexID}) {
    # Add it to active vertices list...
    push @{$This->{ActiveVertices}}, $VertexID;

    # Mark vertex as visited vertex and take it out from the list of vertices to visit...
    $This->{VisitedVertices}{$VertexID} = 1;
    delete $This->{VerticesToVisit}{$VertexID};
  }

  # Set up root vertex, predecessor vertex and distance from root...
  if ($VertexID == $PredecessorVertexID) {
    $This->{VerticesRoots}{$VertexID} = $VertexID;

    $This->{VerticesPredecessors}{$VertexID} = $VertexID;
    if (!exists $This->{VerticesSuccessors}{$VertexID}) {
      @{$This->{VerticesSuccessors}{$VertexID}} = ();
    }

    $This->{VerticesDepth}{$VertexID} = 0;

    if ($This->{TraversalMode} =~ /^Path$/i) {
      $This->_ProcessVisitedPath($VertexID, $PredecessorVertexID);
    }
  }
  else {
    $This->{VerticesRoots}{$VertexID} = $This->{VerticesRoots}{$PredecessorVertexID};

    $This->{VerticesPredecessors}{$VertexID} = $PredecessorVertexID;
    if (!exists $This->{VerticesSuccessors}{$PredecessorVertexID}) {
      @{$This->{VerticesSuccessors}{$PredecessorVertexID}} = ();
    }
    push @{$This->{VerticesSuccessors}{$PredecessorVertexID}}, $VertexID;

    if (!exists $This->{VerticesDepth}{$VertexID}) {
      $This->{VerticesDepth}{$VertexID} = $This->{VerticesDepth}{$PredecessorVertexID} + 1;
    }

    if ($This->{TraversalMode} =~ /^Path$/i) {
      $This->_ProcessVisitedPath($VertexID, $PredecessorVertexID);
      $This->_ProcessVisitedEdge($PredecessorVertexID, $VertexID);
    }
    elsif ($This->{TraversalMode} =~ /^VertexNeighborhood$/i) {
      $This->_ProcessVisitedEdge($PredecessorVertexID, $VertexID);
    }
  }
  return $This;
}

# Process visited path...
#
sub _ProcessVisitedPath {
  my($This, $VertexID, $PredecessorVertexID) = @_;

  # Initialize VerticesPath...
  if (!exists $This->{VerticesPaths}{$VertexID}) {
    @{$This->{VerticesPaths}{$VertexID}} = ();
  }

  if ($VertexID == $PredecessorVertexID) {
    # Starting of a path from root...
    push @{$This->{VerticesPaths}{$VertexID}}, $VertexID;
  }
  else {
    # Setup path for a vertex using path information from predecessor vertex...
    if (exists $This->{CycleClosureVertices}{$PredecessorVertexID}) {
      # Start of a new path from predecessor vertex...
      push @{$This->{VerticesPaths}{$VertexID}}, "${PredecessorVertexID}-${VertexID}";
    }
    else {
      my($PredecessorVertexPath);
      for $PredecessorVertexPath (@{$This->{VerticesPaths}{$PredecessorVertexID}}) {
	push @{$This->{VerticesPaths}{$VertexID}}, "${PredecessorVertexPath}-${VertexID}";
      }
    }
  }
  return $This;
}

# Process visited edge...
#
sub _ProcessVisitedEdge {
  my($This, $VertexID1, $VertexID2) = @_;

  if (!exists $This->{VisitedEdges}->{From}->{$VertexID1}) {
    %{$This->{VisitedEdges}->{From}->{$VertexID1}} = ();
  }
  $This->{VisitedEdges}->{From}->{$VertexID1}->{$VertexID2} = $VertexID2;

  if (!exists $This->{VisitedEdges}->{To}->{$VertexID2}) {
    %{$This->{VisitedEdges}->{To}->{$VertexID2}} = ();
  }
  $This->{VisitedEdges}->{To}->{$VertexID2}->{$VertexID1} = $VertexID1;

  return $This;
}

# Finished processing active vertex...
#
sub _ProcessFinishedVertex {
  my($This, $VertexID) = @_;

  if (!exists $This->{FinishedVertices}{$VertexID}) {
    $This->{FinishedVertices}{$VertexID} = $VertexID;
    # Add vertex to list of vertices found by traversal...
    push @{$This->{Vertices}}, $VertexID;
  }

  # Any active vertices left...
  if (!@{$This->{ActiveVertices}}) {
    return $This;
  }

  # Take it off active vertices list...
  if ($This->{TraversalType} =~ /^(DFS|DFSWithLimit)$/i) {
    # For DFS, it's last vertex in LIFO queue...
    pop @{$This->{ActiveVertices}};
  }
  elsif ($This->{TraversalType} =~ /^(BFS|BFSWithLimit)$/i) {
    # For BFS, it's first vertex in  FIFO queue...
    shift @{$This->{ActiveVertices}};
  }
  return $This;
}

# Mark all other downstream neighbor vertices to be ignored from any further
# processing...
#
sub _IgnoreDownstreamNeighbors {
  my($This, $VertexID, $PredecessorVertexID) = @_;

  if (exists $This->{VerticesToVisit}{$VertexID}) {
    # Mark vertex as visited vertex and take it out from the list of vertices to visit...
    $This->{VisitedVertices}{$VertexID} = 1;
    delete $This->{VerticesToVisit}{$VertexID};

    if (defined($PredecessorVertexID) && $This->{TraversalMode} =~ /^(Path|VertexNeighborhood)$/i) {
      $This->_ProcessVisitedEdge($VertexID, $PredecessorVertexID);
    }
  }
  my($NeighborVertexID, @NeighborsVertexIDs);

  @NeighborsVertexIDs = ();
  @NeighborsVertexIDs = $This->{Graph}->GetNeighbors($VertexID);
  NEIGHBOR: for $NeighborVertexID (@NeighborsVertexIDs) {
    if (!exists $This->{VerticesToVisit}{$NeighborVertexID}) {
      # Avoid going back to predecessor vertex which has already been ignored...
      next NEIGHBOR;
    }
    $This->_IgnoreDownstreamNeighbors($NeighborVertexID, $VertexID);
  }
  return $This;
}

# Is it a visited edge?
#
sub _IsVisitedEdge {
  my($This, $VertexID1, $VertexID2) = @_;

  if (exists $This->{VisitedEdges}->{From}->{$VertexID1}) {
    if (exists $This->{VisitedEdges}->{From}->{$VertexID1}->{$VertexID2}) {
      return 1;
    }
  }
  elsif (exists $This->{VisitedEdges}->{To}->{$VertexID2}) {
    if (exists $This->{VisitedEdges}->{To}->{$VertexID2}->{$VertexID1}) {
      return 1;
    }
  }
  return 0;
}

# Is it a cycle closure edge?
#
# Notes:
#   . Presence of VertexID2 in DFS path traversed for VertexID1 make it a cycle
#     closure edge...
#
sub _IsCycleClosureEdge {
  my($This, $VertexID1, $VertexID2) = @_;

  if (!exists $This->{VerticesPaths}{$VertexID1}) {
    return 0;
  }
  my($Path);
  for $Path (@{$This->{VerticesPaths}{$VertexID1}}) {
    if (($Path =~ /-$VertexID2-/ || $Path =~ /^$VertexID2-/ || $Path =~ /-$VertexID2$/)) {
      return 1;
    }
  }
  return 0;
}

# Search paths starting from a specified vertex with no sharing of edges in paths traversed.
# By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
# completing the cycle.
#
sub PerformPathsSearch {
  my($This, $StartVertexID, $AllowCycles) = @_;

  # Make sure start vertex is defined...
  if (!defined $StartVertexID) {
    carp "Warning: ${ClassName}->PerformPathsSearch: No paths search performed: Start vertex  must be specified...";
    return undef;
  }

  # Make sure start vertex is valid...
  if (!$This->{Graph}->HasVertex($StartVertexID)) {
    carp "Warning: ${ClassName}->PerformPathsSearch: No paths search performed: Vertex $StartVertexID doesn't exist...";
    return undef;
  }

  if (!defined $AllowCycles) {
    $AllowCycles = 1;
  }

  # Perform paths search...
  return $This->_PerformPathsSearch("AllLengths", $StartVertexID, $AllowCycles);
}

# Search paths starting from a specified vertex with length upto a specified length
# with no sharing of edges in paths traversed...
#
# By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
# completing the cycle.
#
sub PerformPathsSearchWithLengthUpto {
  my($This, $StartVertexID, $Length, $AllowCycles) = @_;

  return $This->_PerformPathsSearchWithLength("LengthUpto", $StartVertexID, $Length, $AllowCycles);
}

# Search paths starting from a specified vertex with specified length
# with no sharing of edges in paths traversed...
#
# By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
# completing the cycle.
#
sub PerformPathsSearchWithLength {
  my($This, $StartVertexID, $Length, $AllowCycles) = @_;

  return $This->_PerformPathsSearchWithLength("Length", $StartVertexID, $Length, $AllowCycles);
}


# Search paths starting from a specified vertex with length upto a specified length
# with no sharing of edges in paths traversed...
#
# By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
# completing the cycle.
#
sub _PerformPathsSearchWithLength {
  my($This, $Mode, $StartVertexID, $Length, $AllowCycles) = @_;

  # Make sure both start vertex and length are defined...
  if (!defined $StartVertexID) {
    carp "Warning: ${ClassName}->_PerformPathsSearchWithLength: No paths search performed: Start vertex  must be specified...";
    return undef;
  }
  if (!defined $Length) {
    carp "Warning: ${ClassName}->_PerformPathsSearchWithLength: No paths search performed: Length must be specified...";
    return undef;
  }

  if (!defined $AllowCycles) {
    $AllowCycles = 1;
  }

  # Make sure both start vertex and length are valid...
  if (!$This->{Graph}->HasVertex($StartVertexID)) {
    carp "Warning: ${ClassName}->_PerformPathsSearchWithLength: No paths search performed: Vertex $StartVertexID doesn't exist...";
    return undef;
  }

  if ($Length < 1) {
    carp "Warning: ${ClassName}->_PerformPathsSearchWithLength: No paths search performed: Specified length, $Length, must be a positive integer with value greater than 1...";
    return undef;
  }

  # Perform paths search...
  return $This->_PerformPathsSearch($Mode, $StartVertexID, $AllowCycles, $Length);
}

# Search all paths starting from a specified vertex with sharing of edges in paths traversed...
# By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
# completing the cycle.
#
sub PerformAllPathsSearch {
  my($This, $StartVertexID, $AllowCycles) = @_;

  # Make sure start vertex is defined...
  if (!defined $StartVertexID) {
    carp "Warning: ${ClassName}->PerformAllPathsSearch: No paths search performed: Start vertex  must be specified...";
    return undef;
  }

  # Make sure start vertex is valid...
  if (!$This->{Graph}->HasVertex($StartVertexID)) {
    carp "Warning: ${ClassName}->PerformAllPathsSearch: No paths search performed: Vertex $StartVertexID doesn't exist...";
    return undef;
  }

  if (!defined $AllowCycles) {
    $AllowCycles = 1;
  }

  # Perform paths search...
  return $This->_PerformAllPathsSearch("AllLengths", $StartVertexID, $AllowCycles);
}

# Search all paths starting from a specified vertex with length upto a specified length with sharing of
# edges in paths traversed.
#
# By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
# completing the cycle.
#
sub PerformAllPathsSearchWithLengthUpto {
  my($This, $StartVertexID, $Length, $AllowCycles) = @_;

  return $This->_PerformAllPathsSearchWithLength("LengthUpto", $StartVertexID, $Length, $AllowCycles);
}

# Search all paths starting from a specified vertex with  specified length with sharing of
# edges in paths traversed.
#
# By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
# completing the cycle.
#
sub PerformAllPathsSearchWithLength {
  my($This, $StartVertexID, $Length, $AllowCycles) = @_;

  return $This->_PerformAllPathsSearchWithLength("Length", $StartVertexID, $Length, $AllowCycles);
}

# Search all paths starting from a specified vertex with length upto a specified length with sharing of
# edges in paths traversed.
#
# By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
# completing the cycle.
#
sub _PerformAllPathsSearchWithLength {
  my($This, $Mode, $StartVertexID, $Length, $AllowCycles) = @_;

  # Make sure both start vertex and length are defined...
  if (!defined $StartVertexID) {
    carp "Warning: ${ClassName}->_PerformAllPathsSearchWithLength: No paths search performed: Start vertex  must be specified...";
    return undef;
  }
  if (!defined $Length) {
    carp "Warning: ${ClassName}->_PerformAllPathsSearchWithLength: No paths search performed: Length must be specified...";
    return undef;
  }

  if (!defined $AllowCycles) {
    $AllowCycles = 1;
  }

  # Make sure both start vertex and length are valid...
  if (!$This->{Graph}->HasVertex($StartVertexID)) {
    carp "Warning: ${ClassName}->_PerformAllPathsSearchWithLength: No paths search performed: Vertex $StartVertexID doesn't exist...";
    return undef;
  }

  if ($Length < 1) {
    carp "Warning: ${ClassName}->_PerformAllPathsSearchWithLength: No paths search performed: Specified length, $Length, must be a positive integer with value greater than 1...";
    return undef;
  }

  # Perform paths search...
  return $This->_PerformAllPathsSearch($Mode, $StartVertexID, $AllowCycles, $Length);
}

# Search paths between two vertices...
#
sub PerformPathsSearchBetween {
  my($This, $StartVertexID, $EndVertexID) = @_;

  # Make sure start and end vertices are defined...
  if (!defined $StartVertexID) {
    carp "Warning: ${ClassName}->PerformPathsSearchBetweeb: No paths search performed: Start vertex  must be specified...";
    return undef;
  }
  if (!defined $EndVertexID) {
    carp "Warning: ${ClassName}->PerformPathsSearchBetweeb: No paths search performed: EndVertex vertex  must be specified...";
    return undef;
  }
  # Make sure start and end vertices are valid...
  if (!$This->{Graph}->HasVertex($StartVertexID)) {
    carp "Warning: ${ClassName}->PerformPathsSearchBetween: No paths search performed: Vertex $StartVertexID doesn't exist...";
    return undef;
  }
  if (!$This->{Graph}->HasVertex($EndVertexID)) {
    carp "Warning: ${ClassName}->PerformPathsSearchBetween: No paths search performed: Vertex $EndVertexID doesn't exist...";
    return undef;
  }

  # Perform paths search...
  return $This->_PerformPathsSearchBetween($StartVertexID, $EndVertexID);
}

# Search paths starting from root vertex with no sharing of edges...
#
# Notes:
#   . Possible paths searche modes are: DFSPathsWithLimit, DFSPaths. And each
#     of these modes supports any combination of two options: CommonEdges, Cycles.
#     Default for CommonEdges - No; Cycles - No.
#
sub _PerformPathsSearch {
  my($This, $Mode, $RootVertexID, $AllowCycles, $Length) = @_;

  # Perform DFS path search...

  $This->{TraversalMode} = 'Path';

  if ($Mode =~ /^(LengthUpto|Length)$/i) {
    my($DepthLimit);

    $DepthLimit = $Length - 1;
    $This->{TraversalType} = 'DFSWithLimit';
    $This->{DepthLimit} = $DepthLimit;
  }
  else {
    $This->{TraversalType} = 'DFS';
  }
  if (defined $RootVertexID) {
    $This->{RootVertex} = $RootVertexID;
  }

  $This->{AllowPathCycles} = $AllowCycles;

  # Perform search...
  $This->_TraverseGraph();

  # Make sure traversal did get the root vertex...
  if (!exists $This->{VerticesDepth}{$RootVertexID}) {
    return $This;
  }
  if ($Mode =~ /^Length$/i) {
    push @{$This->{Paths}}, $This->_CollectPathsTraversedDuringPathsSearchWithLength($Length);
  }
  else {
    push @{$This->{Paths}}, $This->_CollectPathsTraversedDuringPathsSearch();
  }

  return $This;
}

# Search all paths starting from root vertex with sharing of edges...
#
sub _PerformAllPathsSearch {
  my($This, $Mode, $RootVertexID, $AllowCycles, $Length) = @_;

  # Perform DFS path search...

  $This->{TraversalMode} = 'AllPaths';

  if ($Mode =~ /^(LengthUpto|Length)$/i) {
    my($DepthLimit);

    $DepthLimit = $Length - 1;
    $This->{TraversalType} = 'DFSWithLimit';
    $This->{DepthLimit} = $DepthLimit;
  }
  else {
    $This->{TraversalType} = 'DFS';
  }
  $This->{RootVertex} = $RootVertexID;
  $This->{AllowPathCycles} = $AllowCycles;

  # Traverse all paths search using DFS search...
  $This->_TraverseAllPathsInGraph($Mode, $Length);

  return $This;
}

# Travese all paths in graph starting from a specified root vertex...
#
sub _TraverseAllPathsInGraph {
  my($This, $Mode, $Length) = @_;

  if ($This->{TraversalMode} !~ /^AllPaths$/i) {
    return $This;
  }
  my($CurrentVertexID, $PredecessorVertexID, $CurrentDepth, $CurrentPath);

  $CurrentVertexID = $This->{RootVertex};
  $PredecessorVertexID = $CurrentVertexID;
  $CurrentDepth = 0;
  $CurrentPath = "$CurrentVertexID";

  $This->_TraverseAllPaths($CurrentVertexID, $PredecessorVertexID, $CurrentDepth, $CurrentPath);

  if ($Mode =~ /^Length$/i) {
    push @{$This->{Paths}}, $This->_CollectPathsTraversedDuringPathsSearchWithLength($Length);
  }
  else {
    push @{$This->{Paths}}, $This->_CollectPathsTraversedDuringPathsSearch();
  }

  return $This;
}

# Traverse and collect all paths recuresively..
#
sub _TraverseAllPaths {
  my($This, $CurrentVertexID, $PredecessorVertexID, $CurrentDepth, $CurrentPath) = @_;

  # Save path traversed for current vertex..
  if (!exists $This->{VerticesPaths}{$CurrentVertexID}) {
    @{$This->{VerticesPaths}{$CurrentVertexID}} = ();
    $This->{VerticesDepth}{$CurrentVertexID} = 0;
  }
  push @{$This->{VerticesPaths}{$CurrentVertexID}}, $CurrentPath;
  $This->{VerticesDepth}{$CurrentVertexID} = $CurrentDepth;

  $CurrentDepth++;
  if (exists $This->{DepthLimit}) {
    if ($CurrentDepth > $This->{DepthLimit}) {
      # Nothing more to do...
      return $This;
    }
  }
  my($NeighborVertexID, $NewPath);

  NEIGHBOR: for $NeighborVertexID ($This->{Graph}->GetNeighbors($CurrentVertexID)) {
    if ($NeighborVertexID == $PredecessorVertexID) {
      next NEIGHBOR;
    }
    if ($This->_IsVertexInTraversedPath($NeighborVertexID, $CurrentPath)) {
      # It's a cycle...
      if ($This->{AllowPathCycles}) {
	$NewPath = "${CurrentPath}-${NeighborVertexID}";
	if (!exists $This->{VerticesPaths}{$NeighborVertexID}) {
	  @{$This->{VerticesPaths}{$NeighborVertexID}} = ();
	}
	push @{$This->{VerticesPaths}{$NeighborVertexID}}, $NewPath;
      }
      next NEIGHBOR;
    }
    $NewPath = "${CurrentPath}-${NeighborVertexID}";
    $This->_TraverseAllPaths($NeighborVertexID, $CurrentVertexID, $CurrentDepth, $NewPath);
  }
  return $This;
}

# Is vertex already in traversed path?
#
sub _IsVertexInTraversedPath {
  my($This, $VertexID, $Path) = @_;

  return ($Path =~ /-$VertexID-/ || $Path =~ /^$VertexID-/ || $Path =~ /-$VertexID$/) ? 1 : 0;
}

# Collect all paths traversed during Path TraversalMode and sort 'em in
# ascending order of lengths
#
sub _CollectPathsTraversedDuringPathsSearch {
  my($This) = @_;
  my($VertexID, @Paths, @SortedPaths);

  @Paths = (); @SortedPaths = ();

  # Create path objects from path vertex strings...
  for $VertexID (keys %{$This->{VerticesPaths}}) {
    push @Paths, map { new Graph::Path(split /-/, $_) } @{$This->{VerticesPaths}{$VertexID}};
  }

  # Sort paths in ascending order of lengths...
  push @SortedPaths, sort { $a->GetLength() <=> $b->GetLength() } @Paths;

  return @SortedPaths;
}

# Collect paths traversed during Path TraversalMode with specific length...
#
sub _CollectPathsTraversedDuringPathsSearchWithLength {
  my($This, $Length) = @_;
  my($VertexID, $Depth, $PathString, @VertexIDs, @Paths);

  @Paths = ();
  $Depth = $Length - 1;

  # Create path objects from path vertex strings...
  VERTEXID: for $VertexID (keys %{$This->{VerticesPaths}}) {
    if ($This->{VerticesDepth}{$VertexID} != $Depth) {
      next VERTEXID;
    }
    # For vertices involved in cycles, the path might also contain some shorter paths. So check
    # the lengths before its collection...
    PATHSTRING: for $PathString (@{$This->{VerticesPaths}{$VertexID}}) {
      @VertexIDs = split /-/, $PathString;
      if ($Length != @VertexIDs) {
	next PATHSTRING;
      }
      push @Paths, new Graph::Path(@VertexIDs);
    }
  }
  return @Paths;
}

# Collect paths traversed during Vertex TraversalMode...
#
sub _CollectPathsTraversedDuringVertexSearch {
  my($This, $RootVertexID) = @_;
  my($Depth, @Paths, @VerticesAtDepth);
  @Paths = ();

  # Get vertices at specific depths...
  @VerticesAtDepth = ();
  @VerticesAtDepth = $This->_CollectVerticesAtSpecificDepths();
  if (!@VerticesAtDepth) {
    return @Paths;
  }

  # Make sure search found only one root vertex and it corresponds to
  # what was specified...
  $Depth = 0;
  if ((@{$VerticesAtDepth[$Depth]} > 1) || ($VerticesAtDepth[$Depth][0] != $RootVertexID)) {
    carp "Warning: ${ClassName}->_PerformPathsSearch: No paths found: Root vertex, $VerticesAtDepth[$Depth][0], identified by paths traversal doen't match specified root vertex $RootVertexID...";
    return @Paths;
  }

  # Setup root vertex at depth 0. And set its path...
  my($Path, $VertexID, $SuccessorVertexID, @VertexIDs, %PathAtVertex);
  %PathAtVertex = ();
  $PathAtVertex{$RootVertexID} = new Graph::Path($RootVertexID);

  for $Depth (0 .. $#VerticesAtDepth) {
    # Go over all vertices at current depth...
    VERTEX: for $VertexID (@{$VerticesAtDepth[$Depth]}) {
      if (!exists $This->{VerticesSuccessors}{$VertexID}) {
	next VERTEX;
      }
      # Get vertices for current path...
      @VertexIDs = ();
      push @VertexIDs, $PathAtVertex{$VertexID}->GetVertices;

      # Expand path to successor vertex found during traversal...
      for $SuccessorVertexID (@{$This->{VerticesSuccessors}{$VertexID}}) {
	$Path = new Graph::Path(@VertexIDs);
	$Path->AddVertex($SuccessorVertexID);
	$PathAtVertex{$SuccessorVertexID} = $Path;
      }
    }
  }
  # Sort paths in ascending order of lengths...
  push @Paths, sort { $a->GetLength() <=> $b->GetLength() } values %PathAtVertex;

  return @Paths;
}

# Collect vertices at specific depths. Depth values start from 0...
#
sub _CollectVerticesAtSpecificDepths {
  my($This) = @_;
  my($VertexID, $Depth, @VerticesAtDepth);

  @VerticesAtDepth = ();
  while (($VertexID, $Depth) = each %{$This->{VerticesDepth}}) {
    push @{$VerticesAtDepth[$Depth]}, $VertexID;
  }
  return @VerticesAtDepth;
}

# Collect vertices, along with their successors, at specific depths and return  a list containing references to
# lists with first value corresponding to vertex ID and second value a reference to a list containing
# its successors.
#
# Depth values start from 0...
#
sub _CollectVerticesWithSuccessorsAtSpecificDepths {
  my($This) = @_;
  my($VertexID, $Depth, @VerticesWithSuccessorsAtDepth);

  @VerticesWithSuccessorsAtDepth = ();
  while (($VertexID, $Depth) = each %{$This->{VerticesDepth}}) {
    my(@VertexWithSuccessors, @VertexSuccessors);

    @VertexWithSuccessors = (); @VertexSuccessors = ();
    if (exists $This->{VerticesSuccessors}{$VertexID}) {
      push @VertexSuccessors, @{$This->{VerticesSuccessors}{$VertexID}};
    }
    push @VertexWithSuccessors, ($VertexID, \@VertexSuccessors);
    # Multiple entries for a vertex and its successors could be present at a specific depth...
    push @{$VerticesWithSuccessorsAtDepth[$Depth]}, \@VertexWithSuccessors;
  }
  return @VerticesWithSuccessorsAtDepth;
}

# Search paths between two vertices...
#
sub _PerformPathsSearchBetween {
  my($This, $RootVertexID, $TargetVertexID) = @_;
  my($DepthLimit);

  # Perform a targeted DFS search...
  $DepthLimit = undef;
  $This->_PerformVertexSearch("DFS", $RootVertexID, $DepthLimit, $TargetVertexID);

  my($Path);
  $Path = $This->_CollectPathBetween($RootVertexID, $TargetVertexID);

  if (defined $Path) {
    push  @{$This->{Paths}}, $Path;
  }
  return $This;
}

# Collect path between root and target vertex after the search...
#
sub _CollectPathBetween {
  my($This, $RootVertexID, $TargetVertexID) = @_;

  # Does a path from root to target vertex exist?
  if (!(exists($This->{VerticesRoots}{$TargetVertexID}) && ($This->{VerticesRoots}{$TargetVertexID} == $RootVertexID))) {
    return undef;
  }

  # Add target vertex ID path vertices...
  my($VertexID, $Path, @VertexIDs);
  @VertexIDs = ();
  $VertexID = $TargetVertexID;
  push @VertexIDs, $VertexID;

  # Backtrack to root vertex ID...
  while ($This->{VerticesPredecessors}{$VertexID} != $VertexID) {
    $VertexID = $This->{VerticesPredecessors}{$VertexID};
    push @VertexIDs, $VertexID;
  }

  # Create path from target to root and reverse it...
  $Path = new Graph::Path(@VertexIDs);
  $Path->Reverse();

  return $Path;
}

# Search vertices around specified root vertex with in specific neighborhood radius...
#
sub PerformNeighborhoodVerticesSearchWithRadiusUpto {
  my($This, $StartVertexID, $Radius) = @_;

  # Make sure both start vertex and radius are defined...
  if (!defined $StartVertexID) {
    carp "Warning: ${ClassName}->PerformNeighborhoodVerticesSearchWithRadiusUpto: No vertices search performed: Start vertex  must be specified...";
    return undef;
  }
  if (!defined $Radius) {
    carp "Warning: ${ClassName}->PerformNeighborhoodVerticesSearchWithRadiusUpto: No vertices search performed: Radius must be specified...";
    return undef;
  }

  # Make sure both start vertex and length are valid...
  if (!$This->{Graph}->HasVertex($StartVertexID)) {
    carp "Warning: ${ClassName}->PerformNeighborhoodVerticesSearchWithRadiusUpto: No vertices search performed: Vertex $StartVertexID doesn't exist...";
    return undef;
  }
  if ($Radius < 0) {
    carp "Warning: ${ClassName}->PerformNeighborhoodVerticesSearchWithRadiusUpto: No vertices search performed: Specified radius, $Radius, must be a positive integer...";
    return undef;
  }

  # Perform vertices search...
  return $This->_PerformNeighborhoodVerticesSearch("RadiusUpto", $StartVertexID, $Radius);
}

# Search vertices around specified root vertex...
#
sub PerformNeighborhoodVerticesSearch {
  my($This, $StartVertexID) = @_;

  # Make sure start vertex is defined...
  if (!defined $StartVertexID) {
    carp "Warning: ${ClassName}->PerformNeighborhoodVerticesSearch: No vertices search performed: Start vertex  must be specified...";
    return undef;
  }

  # Make sure start vertex is valid...
  if (!$This->{Graph}->HasVertex($StartVertexID)) {
    carp "Warning: ${ClassName}->PerformNeighborhoodVerticesSearch: No vertices search performed: Vertex $StartVertexID doesn't exist...";
    return undef;
  }
  # Perform vertices search...
  return $This->_PerformNeighborhoodVerticesSearch("AllRadii", $StartVertexID);
}

# Search vertices around specified root vertex with in specific neighborhood radius along with
# identification of successors of each vertex found during the search...
#
sub PerformNeighborhoodVerticesSearchWithSuccessorsAndRadiusUpto {
  my($This, $StartVertexID, $Radius) = @_;

  # Make sure both start vertex and radius are defined...
  if (!defined $StartVertexID) {
    carp "Warning: ${ClassName}->PerformNeighborhoodVerticesSearchWithSuccessorsAndRadiusUpto: No vertices search performed: Start vertex  must be specified...";
    return undef;
  }
  if (!defined $Radius) {
    carp "Warning: ${ClassName}->PerformNeighborhoodVerticesSearchWithSuccessorsAndRadiusUpto: No vertices search performed: Radius must be specified...";
    return undef;
  }

  # Make sure both start vertex and length are valid...
  if (!$This->{Graph}->HasVertex($StartVertexID)) {
    carp "Warning: ${ClassName}->PerformNeighborhoodVerticesSearchWithSuccessorsAndRadiusUpto: No vertices search performed: Vertex $StartVertexID doesn't exist...";
    return undef;
  }
  if ($Radius < 0) {
    carp "Warning: ${ClassName}->PerformNeighborhoodVerticesSearchWithSuccessorsAndRadiusUpto: No vertices search performed: Specified radius, $Radius, must be a positive integer...";
    return undef;
  }

  # Perform vertices search...
  return $This->_PerformNeighborhoodVerticesSearch("WithSuccessorsAndRadiusUpto", $StartVertexID, $Radius);
}

# Search vertices around specified root vertex along with identification of
# successors of each vertex found during the search...
#
sub PerformNeighborhoodVerticesSearchWithSuccessors {
  my($This, $StartVertexID) = @_;

  # Make sure start vertex is defined...
  if (!defined $StartVertexID) {
    carp "Warning: ${ClassName}->PerformNeighborhoodVerticesSearchWithSuccessors: No vertices search performed: Start vertex  must be specified...";
    return undef;
  }

  # Make sure start vertex is valid...
  if (!$This->{Graph}->HasVertex($StartVertexID)) {
    carp "Warning: ${ClassName}->PerformNeighborhoodVerticesSearchWithSuccessors: No vertices search performed: Vertex $StartVertexID doesn't exist...";
    return undef;
  }
  # Perform vertices search...
  return $This->_PerformNeighborhoodVerticesSearch("WithSuccessorsAndAllRadii", $StartVertexID);
}

# Search vertices at successive neighborhood radii levels...
#
sub _PerformNeighborhoodVerticesSearch {
  my($This, $Mode, $RootVertexID, $Radius) = @_;
  my($DepthLimit, $AllowCycles);

  $DepthLimit = defined $Radius ? $Radius :  undef;
  $AllowCycles = undef;

  # Perform BFS search...
  if ($Mode =~ /^RadiusUpto$/i) {
    $This->_PerformVertexNeighborhoodSearch("BFSWithLimit", $RootVertexID, $DepthLimit);
  }
  elsif ($Mode =~ /^(AllRadii)$/i) {
    $This->_PerformVertexNeighborhoodSearch("BFS", $RootVertexID);
  }
  elsif ($Mode =~ /^WithSuccessorsAndRadiusUpto$/i) {
    $AllowCycles = 1;
    $This->_PerformVertexNeighborhoodSearch("BFSWithLimit", $RootVertexID, $DepthLimit, $AllowCycles);
  }
  elsif ($Mode =~ /^WithSuccessorsAndAllRadii$/i) {
    $AllowCycles = 1;
    $This->_PerformVertexNeighborhoodSearch("BFSWithLimit", $RootVertexID, $DepthLimit, $AllowCycles);
  }

  # Make sure traversal did get the root vertex...
  if (!exists $This->{VerticesDepth}{$RootVertexID}) {
    return $This;
  }

  if ($Mode =~ /^(RadiusUpto|AllRadii)$/i) {
    push @{$This->{VerticesNeighborhoods}}, $This->_CollectVerticesAtSpecificDepths();
  }
  elsif ($Mode =~ /^(WithSuccessorsAndRadiusUpto|WithSuccessorsAndAllRadii)$/i) {
    push @{$This->{VerticesNeighborhoodsWithSuccessors}}, $This->_CollectVerticesWithSuccessorsAtSpecificDepths();
  }

  return $This;
}

# Perform appropriate vertex search...
#
sub _PerformVertexNeighborhoodSearch {
  my($This, $SearchType, $RootVertexID, $DepthLimit, $AllowCycles) = @_;

  # Setup search...
  $This->{TraversalMode} = 'VertexNeighborhood';
  $This->{TraversalType} = $SearchType;

  if (defined $RootVertexID) {
    $This->{RootVertex} = $RootVertexID;
  }
  if (defined $DepthLimit) {
    $This->{DepthLimit} = $DepthLimit;
  }
  if (defined $AllowCycles) {
    $This->{AllowVertexCycles} = $AllowCycles;
  }

  # Perform search...
  return $This->_TraverseGraph();
}

# Get orderded list of vertices after DFS/BFS search...
#
sub GetVertices {
  my($This) = @_;

  return wantarray ? @{$This->{Vertices}} : scalar @{$This->{Vertices}};
}

# Get a hash list containing vertex and root vertex as key/value pair for all vertices
# ordered using DFS/BFS search available via GetVertices method...
#
sub GetVerticesRoots {
  my($This) = @_;

  return %{$This->{VerticesRoots}};
}

# Get a list containing lists of vertices in connected components of graph after DFS/BFS
# search...
#
# Note:
#   . List is sorted in descending order of number of vertices in each connected component.
#
sub GetConnectedComponentsVertices {
  my($This) = @_;
  my($VertexID, $VertexRoot, @ConnectedVertices, %VerticesAtRoot);

  @ConnectedVertices = ();
  %VerticesAtRoot = ();
  for $VertexID (@{$This->{Vertices}}) {
    $VertexRoot = $This->{VerticesRoots}{$VertexID};
    if (!exists $VerticesAtRoot{$VertexRoot}) {
      @{$VerticesAtRoot{$VertexRoot}} = ();
    }
    push @{$VerticesAtRoot{$VertexRoot}}, $VertexID;
  }
  push @ConnectedVertices, sort { @{$b} <=> @{$a} }  values %VerticesAtRoot;

  return wantarray ? @ConnectedVertices : scalar @ConnectedVertices;
}

# Get predecessor vertices...
#
sub GetVerticesPredecessors {
  my($This) = @_;

  return %{$This->{VerticesPredecessors}};
}

# Get a hash list containing vertex and depth from root vertex as key/value pair for all vertices
# ordered using DFS/BFS search available via GetVertices method...
#
sub GetVerticesDepth {
  my($This) = @_;

  return %{$This->{VerticesDepth}};
}

# Get paths found during paths search...
#
sub GetPaths {
  my($This) = @_;

  return wantarray ? @{$This->{Paths}} : scalar @{$This->{Paths}};
}

# Get vertices collected at various neighborhood radii...
#
sub GetVerticesNeighborhoods {
  my($This) = @_;

  return wantarray ? @{$This->{VerticesNeighborhoods}} : scalar @{$This->{VerticesNeighborhoods}};
}

# Get vertices, along with their successor vertices, collected at various neighborhood radii as
# a list containing references to lists with first value corresponding to vertex ID and second value
# a reference to a list containing its successors.
#
sub GetVerticesNeighborhoodsWithSuccessors {
  my($This) = @_;

  return wantarray ? @{$This->{VerticesNeighborhoodsWithSuccessors}} : scalar @{$This->{VerticesNeighborhoodsWithSuccessors}};
}

# Return a string containg data for PathsTraversal object...
sub StringifyPathsTraversal {
  my($This) = @_;
  my($PathsTraversalString);

  $PathsTraversalString = "PathsTraversalMode: " . $This->{TraversalMode};
  $PathsTraversalString .= "; PathsTraversalType: " . $This->{TraversalType};

  # Vertices ordered by traversal...
  $PathsTraversalString .= "; Vertices: " . join(' ', @{$This->{Vertices}});

  # Stringify depths of vertices...
  $PathsTraversalString .= "; " . $This->StringifyVerticesDepth();

  # Stringify roots of vertices...
  $PathsTraversalString .= "; " . $This->StringifyVerticesRoots();

  # Stringify predecessor of vertices...
  $PathsTraversalString .= "; " . $This->StringifyVerticesPredecessors();

  # Stringify successor vertices...
  $PathsTraversalString .= "; " . $This->StringifyVerticesSuccessors();

  # Stringify paths...
  $PathsTraversalString .= "; " . $This->StringifyPaths();

  # Stringify vertices neighborhoods...
  $PathsTraversalString .= "; " . $This->StringifyVerticesNeighborhoods();

  # Stringify vertices neighborhoods with successors...
  $PathsTraversalString .= "; " . $This->StringifyVerticesNeighborhoodsWithSuccessors();

  return $PathsTraversalString;
}

# Stringify vertices depth...
#
sub StringifyVerticesDepth {
  my($This) = @_;
  my($VertexID, $VertexDepth, $DepthString);

  if (!@{$This->{Vertices}}) {
    $DepthString = "<Vertex-Depth>: None";
    return $DepthString;
  }

  $DepthString = "<Vertex-Depth>: ";
  for $VertexID (@{$This->{Vertices}}) {
    $VertexDepth = $This->{VerticesDepth}{$VertexID};
    $DepthString .= " <$VertexID-$VertexDepth>";
  }
  return $DepthString;
}

# Stringify roots of vertices...
#
sub StringifyVerticesRoots {
  my($This) = @_;
  my($VertexID, $RootVertexID, $RootsString);

  if (!@{$This->{Vertices}}) {
    $RootsString = "<Vertex-RootVertex>: None";
    return $RootsString;
  }

  $RootsString = "<Vertex-RootVertex>: ";
  for $VertexID (@{$This->{Vertices}}) {
    $RootVertexID = $This->{VerticesRoots}{$VertexID};
    $RootsString .= " <$VertexID-$RootVertexID>";
  }
  return $RootsString;
}

# Stringify predecessor of vertices...
#
sub StringifyVerticesPredecessors {
  my($This) = @_;
  my($VertexID, $PredecessorVertexID, $PredecessorString);

  if (!@{$This->{Vertices}}) {
    $PredecessorString = "<Vertex-PredecessorVertex>: None";
    return $PredecessorString;
  }

  $PredecessorString = "<Vertex-PredecessorVertex>: ";
  for $VertexID (@{$This->{Vertices}}) {
    $PredecessorVertexID = $This->{VerticesPredecessors}{$VertexID};
    $PredecessorString .= " <$VertexID-$PredecessorVertexID>";
  }
  return $PredecessorString;
}

# Stringify successor vertices...
#
sub StringifyVerticesSuccessors {
  my($This) = @_;
  my($VertexID, $SuccessorString, $VerticesSuccessorsString);

  if (!@{$This->{Vertices}}) {
    $SuccessorString = "<Vertex-VerticesSuccessorsList>: None";
    return $SuccessorString;
  }

  $SuccessorString = "<Vertex-VerticesSuccessorsList>: ";
  for $VertexID (@{$This->{Vertices}}) {
    if (exists($This->{VerticesSuccessors}{$VertexID}) && @{$This->{VerticesSuccessors}{$VertexID}}) {
      $VerticesSuccessorsString = join(',', @{$This->{VerticesSuccessors}{$VertexID}});
    }
    else {
      $VerticesSuccessorsString = "None";
    }
    $SuccessorString .= " <$VertexID-$VerticesSuccessorsString>";
  }
  return $SuccessorString;
}

# Strinigify paths...
#
sub StringifyPaths {
  my($This) = @_;
  my($PathsString, $Path);

  if (!@{$This->{Paths}}) {
    $PathsString = "Paths: None";
    return $PathsString;
  }

  my($FirstPath);
  $PathsString = "Paths: ";
  $FirstPath = 1;
  for $Path (@{$This->{Paths}}) {
    if ($FirstPath) {
      $FirstPath = 0;
    }
    else {
      $PathsString .= " ";
    }
    $PathsString .= "<" . join('-', $Path->GetVertices()) . ">";
  }
  return $PathsString;
}

# Strinigify vertices neighborhoods...
#
sub StringifyVerticesNeighborhoods {
  my($This) = @_;
  my($NeighborhoodsString, $NeighborhoodVerticesString, $Radius);

  if (!@{$This->{VerticesNeighborhoods}}) {
    $NeighborhoodsString = "<NeighborhoodRadius-NeighborhoodVerticesList>: None";
    return $NeighborhoodsString;
  }
  $NeighborhoodsString = "<NeighborhoodRadius-NeighborhoodVerticesList>:";
  for $Radius (0 .. $#{$This->{VerticesNeighborhoods}}) {
    $NeighborhoodVerticesString = join(',', @{$This->{VerticesNeighborhoods}[$Radius]});
    $NeighborhoodsString .= " <$Radius-$NeighborhoodVerticesString>";
  }

  return $NeighborhoodsString;
}

# Strinigify vertices neighborhoods...
#
sub StringifyVerticesNeighborhoodsWithSuccessors {
  my($This) = @_;
  my($NeighborhoodsString, $NeighborhoodVertexSuccessorsString, $Radius, $NeighborhoodVertericesWithSuccessorsRef, $NeighborhoodVertexWithSuccessorsRef, $VertexID, $NeighborhoodVertexSuccessorsRef);

  if (!@{$This->{VerticesNeighborhoodsWithSuccessors}}) {
    $NeighborhoodsString = "<NeighborhoodRadius-NeighborhoodVertex-NeighborhoodVerticeSuccessorsList>: None";
    return $NeighborhoodsString;
  }
  $NeighborhoodsString = "<NeighborhoodRadius-NeighborhoodVertex-NeighborhoodVerticeSuccessorsList>: None";

  $Radius = 0;
  for $NeighborhoodVertericesWithSuccessorsRef (@{$This->{VerticesNeighborhoodsWithSuccessors}}) {
    for $NeighborhoodVertexWithSuccessorsRef (@{$NeighborhoodVertericesWithSuccessorsRef}) {
      ($VertexID, $NeighborhoodVertexSuccessorsRef) = @{$NeighborhoodVertexWithSuccessorsRef};
      $NeighborhoodVertexSuccessorsString = 'None';
      if (@{$NeighborhoodVertexSuccessorsRef}) {
	$NeighborhoodVertexSuccessorsString = join(',', @{$NeighborhoodVertexSuccessorsRef});
      }
      $NeighborhoodsString .= " <$Radius-$VertexID-$NeighborhoodVertexSuccessorsString>";
    }
    $Radius++;
  }
  return $NeighborhoodsString;
}

# Return a reference to new paths traversal object...
sub Copy {
  my($This) = @_;
  my($NewPathsTraversal);

  $NewPathsTraversal = Storable::dclone($This);

  return $NewPathsTraversal;
}

1;

__END__

=head1 NAME

PathsTraversal

=head1 SYNOPSIS

use Graph::PathsTraversal;

use Graph::PathsTraversal qw(:all);

=head1 DESCRIPTION

B<PathsTraversal> class provides the following methods:

new, Copy, GetConnectedComponentsVertices, GetPaths, GetVertices,
GetVerticesDepth, GetVerticesNeighborhoods,
GetVerticesNeighborhoodsWithSuccessors, GetVerticesPredecessors, GetVerticesRoots,
PerformAllPathsSearch, PerformAllPathsSearchWithLength,
PerformAllPathsSearchWithLengthUpto, PerformBreadthFirstSearch,
PerformBreadthFirstSearchWithLimit, PerformDepthFirstSearch,
PerformDepthFirstSearchWithLimit, PerformNeighborhoodVerticesSearch,
PerformNeighborhoodVerticesSearchWithRadiusUpto,
PerformNeighborhoodVerticesSearchWithSuccessors,
PerformNeighborhoodVerticesSearchWithSuccessorsAndRadiusUpto, PerformPathsSearch,
PerformPathsSearchBetween, PerformPathsSearchWithLength,
PerformPathsSearchWithLengthUpto, StringifyPaths, StringifyPathsTraversal,
StringifyVerticesDepth, StringifyVerticesNeighborhoods,
StringifyVerticesNeighborhoodsWithSuccessors, StringifyVerticesPredecessors,
StringifyVerticesRoots, StringifyVerticesSuccessors

=head2 METHODS

=over 4

=item B<new>

    $PathsTraversal = new Graph::PathsTraversal($Graph);

Using specified I<Graph>, B<new> method creates a new B<PathsTraversal> object and returns
newly created B<PathsTraversal> object.

=item B<Copy>

    $PathsTraversal = $PathsTraversal->Copy();

Copies I<PathsTraversal> and its associated data using B<Storable::dclone> and returns a new
B<PathsTraversal> object.

=item B<GetConnectedComponentsVertices>

    @Components = $PathsTraversal->GetConnectedComponentsVertices();
    $NumOfComponents = $PathsTraversal->GetConnectedComponentsVertices();

Returns an array of B<Components> containing references to arrays of vertex IDs corresponding
to connected components of graph after a search. In scalar context, the number of connected
components is returned.

Connected B<Components> is sorted in descending order of number of vertices in each
connected component.

=item B<GetPaths>

    @Paths = $PathsTraversal->GetPaths();
    $NumOfPaths = $PathsTraversal->GetPaths();

Returns an array of B<Paths> containing references to arrays of vertex IDs corresponding to
to paths traversed in a graph after a search. In scalar context, number of paths is returned.

B<Paths> array is sorted in ascending order of path lengths.

=item B<GetVertices>

    @Vertices = $PathsTraversal->GetVertices();
    $NumOfVertices = $PathsTraversal->GetVertices();

Returns an array containing an ordered list of vertex IDs traversed during a search. In
scalar context, the number of vertices is returned.

=item B<GetVerticesDepth>

    %VerticesDepth = $PathsTraversal->GetVerticesDepth();

Returns a hash I<VerticesDepth> containing vertex ID and depth from root vertex as a key and
value pair for all vertices traversed during a search.

=item B<GetVerticesNeighborhoods>

    @VerticesNeighborhoods =
       $PathsTraversal->GetVerticesNeighborhoods();
    $NumOfVerticesNeighborhoods =
      $PathsTraversal->GetVerticesNeighborhoods();

Returns an array I<VerticesNeighborhoods> containing references to arrays corresponding
to vertices collected at various neighborhood radii around a specified vertex during a vertex
neighborhood search. In scalar context, the number of neighborhoods is returned.

=item B<GetVerticesNeighborhoodsWithSuccessors>

    @VerticesNeighborhoodsWithSucceessors =
       $PathsTraversal->GetVerticesNeighborhoodsWithSuccessors();
    $NumOfVerticesNeighborhoodsWithSucceessors =
      $PathsTraversal->GetVerticesNeighborhoodsWithSuccessors();

Returns an array I<VerticesNeighborhoodsWithSucceessors> containing references to arrays
with first value corresponding to vertex IDs corresponding to a vertex at a specific neighborhood
radius level and second value a reference to an arraty containing its successors.

=item B<GetVerticesPredecessors>

    %VerticesPredecessors = $PathsTraversal->GetVerticesPredecessors();

Returns a hash I<VerticesPredecessors> containing vertex ID and predecessor vertex ID as key
and value pair for all vertices traversed during a search.

=item B<GetVerticesRoots>

    %VerticesRoots = $PathsTraversal->GetVerticesRoots();

Returns a hash I<VerticesPredecessors> containing vertex ID and root vertex ID as a key
and value pair for all vertices traversed during a search.

=item B<PerformAllPathsSearch>

    $PathsTraversal->PerformAllPathsSearch($StartVertexID, [$AllowCycles]);

Searches all paths starting from a I<StartVertexID> with sharing of edges in paths traversed and
returns I<PathsTraversal>.

By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
completing the cycle.

=item B<PerformAllPathsSearchWithLength>

    $PathsTraversal->PerformAllPathsSearchWithLength($StartVertexID,
                     $Length, [$AllowCycles]);

Searches all paths starting from I<StartVertexID> of specific I<Length> with sharing of
edges in paths traversed and returns I<PathsTraversal>.

By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
completing the cycle.

=item B<PerformAllPathsSearchWithLengthUpto>

    $PathsTraversal->PerformAllPathsSearchWithLengthUpto($StartVertexID,
                     $Length, [$AllowCycles]);

Searches all paths starting from I<StartVertexID> of length upto a I<Length> with sharing of
edges in paths traversed and returns I<PathsTraversal>.

By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
completing the cycle.

=item B<PerformBreadthFirstSearch>

    $PathsTraversal->PerformBreadthFirstSearch();

Performs Breadth First Search (BFS) and returns I<PathsTraversal>.

=item B<PerformBreadthFirstSearchWithLimit>

    $PathsTraversal->PerformBreadthFirstSearchWithLimit($DepthLimit,
                     [$RootVertexID]);

Performs BFS with depth up to I<DepthLimit> starting at I<RootVertexID> and returns
I<PathsTraversal>. By default, root vertex ID corresponds to an arbitrary vertex.

=item B<PerformDepthFirstSearch>

    $Return = $PathsTraversal->PerformDepthFirstSearch();

Performs Depth First Search (DFS) and returns I<PathsTraversal>.

=item B<PerformDepthFirstSearchWithLimit>

    $PathsTraversal->PerformDepthFirstSearchWithLimit($DepthLimit,
                     [$RootVertexID]);

Performs DFS with depth up to I<DepthLimit> starting at I<RootVertexID> and returns
I<PathsTraversal>. By default, root vertex ID corresponds to an arbitrary vertex.

=item B<PerformNeighborhoodVerticesSearch>

    $PathsTraversal->PerformNeighborhoodVerticesSearch($StartVertexID);

Searches vertices around I<StartVertexID> at all neighborhood radii and returns
I<PathsTraversal> object.

=item B<PerformNeighborhoodVerticesSearchWithRadiusUpto>

    $PathsTraversal->PerformNeighborhoodVerticesSearchWithRadiusUpto(
                     $StartVertexID, $Radius);

Searches vertices around I<StartVertexID> with neighborhood radius up to I<Radius> and returns
I<PathsTraversal> object.

=item B<PerformNeighborhoodVerticesSearchWithSuccessors>

    $PathsTraversal->PerformNeighborhoodVerticesSearchWithSuccessors(
                     $StartVertexID);

Searches vertices around I<StartVertexID> at all neighborhood radii along with identification of
successor vertices for each vertex found during the traversal and returns I<PathsTraversal>.

=item B<PerformNeighborhoodVerticesSearchWithSuccessorsAndRadiusUpto>

    $PathsTraversal->
                PerformNeighborhoodVerticesSearchWithSuccessorsAndRadiusUpto(
                     $StartVertexID, $Radius);

Searches vertices around I<StartVertexID> with neighborhood radius upto I<Radius> along with
identification of successor vertices for each vertex found during the traversal and returns
I<PathsTraversal>.

=item B<PerformPathsSearch>

    $PathsTraversal->PerformPathsSearch($StartVertexID, [$AllowCycles]);

Searches paths starting from I<StartVertexID> with no sharing of edges in paths traversed and
returns I<PathsTraversal>.

By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
completing the cycle.

=item B<PerformPathsSearchBetween>

    $PathsTraversal->PerformPathsSearchBetween($StartVertexID, $EndVertexID);

Searches paths between I<StartVertexID> and I<EndVertexID> and returns I<PathsTraversal>

=item B<PerformPathsSearchWithLength>

    $PathsTraversal->PerformPathsSearchWithLength($StartVertexID, $Length,
                     [$AllowCycles]);

Searches paths starting from I<StartVertexID>  with length I<Length> with no sharing of
edges in paths traversed and returns I<PathsTraversal>.

By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
completing the cycle.

=item B<PerformPathsSearchWithLengthUpto>

    $PathsTraversal->PerformPathsSearchWithLengthUpto($StartVertexID, $Length,
                     [$AllowCycles]);

Searches paths starting from I<StartVertexID>  with length upto I<Length> with no sharing of
edges in paths traversed and returns I<PathsTraversal>.

By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
completing the cycle.

=item B<StringifyPaths>

    $String = $PathsTraversal->StringifyPaths();

Returns a string containing information about traversed paths in I<PathsTraversal> object

=item B<StringifyPathsTraversal>

    $String = $PathsTraversal->StringifyPathsTraversal();

Returns a string containing information about I<PathsTraversal> object.

=item B<StringifyVerticesDepth>

    $String = $PathsTraversal->StringifyVerticesDepth();

Returns a string containing information about depth of vertices found during search by
I<PathsTraversal> object.

=item B<StringifyVerticesNeighborhoods>

    $String = $PathsTraversal->StringifyVerticesNeighborhoods();

Returns a string containing information about neighborhoods of vertices found during search by
I<PathsTraversal> object.

=item B<StringifyVerticesNeighborhoodsWithSuccessors>

    $String = $PathsTraversal->StringifyVerticesNeighborhoodsWithSuccessors();

Returns a string containing information about neighborhoods of vertices along with their successors
found during search by I<PathsTraversal> object.

=item B<StringifyVerticesPredecessors>

    $String = $PathsTraversal->StringifyVerticesPredecessors();

Returns a string containing information about predecessors of vertices found during search by
I<PathsTraversal> object.

=item B<StringifyVerticesRoots>

    $String = $PathsTraversal->StringifyVerticesRoots();

Returns a string containing information about roots of vertices found during search by
I<PathsTraversal> object.

=item B<StringifyVerticesSuccessors>

    $String = $PathsTraversal->StringifyVerticesSuccessors();

Returns a string containing information about successors of vertices found during search by
I<PathsTraversal> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Graph.pm, Path.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
