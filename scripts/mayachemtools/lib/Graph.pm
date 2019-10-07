package Graph;
#
# File: Graph.pm
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
use Storable ();
use Scalar::Util ();
use Graph::CyclesDetection;
use Graph::PathsTraversal;
use Graph::GraphMatrix;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw(IsGraph);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyGraph';

# Class constructor...
sub new {
  my($Class, @VertexIDs) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializeGraph();

  if (@VertexIDs) { $This->AddVertices(@VertexIDs); }

  return $This;
}

# Initialize object data...
sub _InitializeGraph {
  my($This) = @_;

  %{$This->{Vertices}} = ();

  %{$This->{Edges}} = ();
  %{$This->{Edges}->{From}} = ();
  %{$This->{Edges}->{To}} = ();

  %{$This->{Properties}} = ();
  %{$This->{Properties}->{Graph}} = ();
  %{$This->{Properties}->{Vertices}} = ();
  %{$This->{Properties}->{Edges}} = ();
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Add a vertex...
sub AddVertex {
  my($This, $VertexID) = @_;

  if (!defined $VertexID ) {
    carp "Warning: ${ClassName}->AddVertex: No vertex added: Vertex ID must be specified...";
    return undef;
  }
  if (exists $This->{Vertices}->{$VertexID}) {
    carp "Warning: ${ClassName}->AddVertex: Didn't add vertex $VertexID: Already exists in the graph...";
    return undef;
  }

  $This->{Vertices}->{$VertexID} = $VertexID;

  return $This;
}

# Add vertices to the graph and return graph...
sub AddVertices {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    carp "Warning: ${ClassName}->AddVertices: No vertices added: Vertices list is empty...";
    return undef;
  }

  my($VertexID);
  for $VertexID (@VertexIDs) {
    $This->AddVertex($VertexID);
  }

  return $This;
}

# Delete a vertex...
sub DeleteVertex {
  my($This, $VertexID) = @_;

  if (!defined $VertexID ) {
    carp "Warning: ${ClassName}->DeleteVertex: No vertex deleted: Vertex ID must be specified...";
    return undef;
  }
  if (!$This->HasVertex($VertexID)) {
    carp "Warning: ${ClassName}->DeleteVertex: Didn't delete vertex $VertexID: Vertex $VertexID doesn't exist...";
    return undef;
  }
  $This->_DeleteVertex($VertexID);

  return $This;
}

# Delete vertex...
sub _DeleteVertex {
  my($This, $VertexID) = @_;

  # Delete corresponding edges; the corresponding edge properties are deleted during
  # edges deletetion...
  my(@VertexIDs);
  @VertexIDs = $This->GetEdges($VertexID);
  if (@VertexIDs) {
    $This->DeleteEdges(@VertexIDs);
  }

  # Delete the vertex and any properties associated with vertex...
  $This->DeleteVertexProperties($VertexID);
  delete $This->{Vertices}->{$VertexID};
}

# Delete vertices...
sub DeleteVertices {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    carp "Warning: ${ClassName}->DeleteVertices: No vertices deleted: Vertices list is empty...";
    return undef;
  }
  my($VertexID);
  for $VertexID (@VertexIDs) {
    $This->DeleteVertex($VertexID);
  }

  return $This;
}

# Get vertex data...
sub GetVertex {
  my($This, $VertexID) = @_;

  if (!defined $VertexID) {
    return undef;
  }

  return (exists $This->{Vertices}->{$VertexID}) ? $This->{Vertices}->{$VertexID} : undef;
}

# Get data for all vertices or those specifed in the list. In scalar context, returned
# the number of vertices found.
#
sub GetVertices {
  my($This, @VertexIDs) = @_;
  my($ValuesCount, @VertexValues);

  @VertexValues = ();
  if (@VertexIDs) {
    @VertexValues = map { $This->GetVertex($_) } @VertexIDs;
    $ValuesCount = grep { 1 } @VertexValues;
  }
  else {
    @VertexValues = sort { $a <=> $b } keys %{$This->{Vertices}};
    $ValuesCount = @VertexValues;
  }

  return wantarray ? @VertexValues : $ValuesCount;
}

# Is this vertex present?
sub HasVertex {
  my($This, $VertexID) = @_;

  if (!defined $VertexID) {
    return 0;
  }
  return (exists $This->{Vertices}->{$VertexID}) ? 1 : 0;
}

# Are these vertices present? Return an array containing 1 or 0  for each vertex.
# In scalar context, return number of vertices found.
sub HasVertices {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    return undef;
  }
  my($VerticesCount, @VerticesStatus);

  @VerticesStatus = map { $This->HasVertex($_) } @VertexIDs;
  $VerticesCount = grep { 1 }  @VerticesStatus;

  return wantarray ? @VerticesStatus : $VerticesCount;
}

# Add an edge...
sub AddEdge {
  my($This, $VertexID1, $VertexID2) = @_;

  if (!(defined($VertexID1) && defined($VertexID2))) {
    carp "Warning: ${ClassName}->AddEdge: No edge added: Both vertices must be defined...";
    return undef;
  }
  if (!$This->HasVertex($VertexID1)) {
    carp "Warning: ${ClassName}->AddEdge: Didn't add edge between vertices $VertexID1 and $VertexID2: Vertex $VertexID1 doesn's exist...";
    return undef;
  }
  if (!$This->HasVertex($VertexID2)) {
    carp "Warning: ${ClassName}->AddEdge: Didn't add edge between vertices $VertexID1 and $VertexID2: Vertex $VertexID2 doesn's exist...";
    return undef;
  }
  if ($VertexID1 == $VertexID2) {
    carp "Warning: ${ClassName}->AddEdge: Didn't add edge between vertices $VertexID1 and $VertexID2: Vertices must be different...";
    return undef;
  }
  if ($This->HasEdge($VertexID1, $VertexID2)) {
    carp "Warning: ${ClassName}->AddEdge: Didn't add edge between vertices $VertexID1 and $VertexID2: Edge already exists...";
    return undef;
  }

  if (!exists $This->{Edges}->{From}->{$VertexID1}) {
    %{$This->{Edges}->{From}->{$VertexID1}} = ();
  }
  $This->{Edges}->{From}->{$VertexID1}->{$VertexID2} = $VertexID2;

  if (!exists $This->{Edges}->{To}->{$VertexID2}) {
    %{$This->{Edges}->{To}->{$VertexID2}} = ();
  }
  $This->{Edges}->{To}->{$VertexID2}->{$VertexID1} = $VertexID1;

  return $This;
}

# Add edges...
sub AddEdges {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    carp "Warning: ${ClassName}->AddEdges: No edges added: Vertices list is empty...";
    return undef;
  }
  if (@VertexIDs % 2) {
    carp "Warning: ${ClassName}->AddEdges: No edges added: Invalid vertices data: Input list must contain even number of vertex IDs...";
    return undef;
  }
  my($VertexID1, $VertexID2, $Index);
  for ($Index = 0; $Index < $#VertexIDs; $Index += 2) {
    $VertexID1 = $VertexIDs[$Index]; $VertexID2 = $VertexIDs[$Index + 1];
    $This->AddEdge($VertexID1, $VertexID2);
  }

  return $This;
}

# Delete an edge...
sub DeleteEdge {
  my($This, $VertexID1, $VertexID2) = @_;

  if (!(defined($VertexID1) && defined($VertexID2))) {
    carp "Warning: ${ClassName}->Delete: No edge deleted: Both vertices must be defined...";
    return undef;
  }
  if (!$This->HasVertex($VertexID1)) {
    carp "Warning: ${ClassName}->DeleteEdge: Didn't delete edge between vertices $VertexID1 and $VertexID2: Vertex $VertexID1 doesn's exist...";
    return undef;
  }
  if (!$This->HasVertex($VertexID2)) {
    carp "Warning: ${ClassName}->DeleteEdge: Didn't delete edge between vertices $VertexID1 and $VertexID2: Vertex $VertexID2 doesn's exist...";
    return undef;
  }
  if (!$This->HasEdge($VertexID1, $VertexID2)) {
    carp "Warning: ${ClassName}->DeleteEdge: Didn't delete edge between vertices $VertexID1 and $VertexID2: Edge doesn't exist...";
    return undef;
  }
  $This->_DeleteEdge($VertexID1, $VertexID2);
  $This->_DeleteEdge($VertexID2, $VertexID1);
}

# Delete edge...
sub _DeleteEdge {
  my($This, $VertexID1, $VertexID2) = @_;

  # Delete the edge...
  if (exists $This->{Edges}->{From}->{$VertexID1}) {
    if (exists $This->{Edges}->{From}->{$VertexID1}->{$VertexID2}) {
      delete $This->{Edges}->{From}->{$VertexID1}->{$VertexID2};
    }
    if (! keys %{$This->{Edges}->{From}->{$VertexID1}}) {
      delete $This->{Edges}->{From}->{$VertexID1};
    }
  }

  if (exists $This->{Edges}->{To}->{$VertexID2}) {
    if (exists $This->{Edges}->{To}->{$VertexID2}->{$VertexID1}) {
      delete $This->{Edges}->{To}->{$VertexID2}->{$VertexID1};
    }
    if (! keys %{$This->{Edges}->{To}->{$VertexID2}}) {
      delete $This->{Edges}->{To}->{$VertexID2};
    }
  }

  # Delete properties associated with the edge...
  $This->DeleteEdgeProperties($VertexID1, $VertexID2);
}

# Delete edges...
sub DeleteEdges {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    carp "Warning: ${ClassName}->DeleteEdges: No edges deleted: Vertices list is empty...";
    return undef;
  }
  if (@VertexIDs % 2) {
    carp "Warning: ${ClassName}->DeleteEdges: No edges deleted: Invalid vertices data: Input list must contain even number of vertex IDs...";
    return undef;
  }
  my($VertexID1, $VertexID2, $Index);
  for ($Index = 0; $Index < $#VertexIDs; $Index += 2) {
    $VertexID1 = $VertexIDs[$Index]; $VertexID2 = $VertexIDs[$Index + 1];
    $This->DeleteEdge($VertexID1, $VertexID2);
  }

  return $This;
}

# Does the edge defiend by a vertex pair exists? Edges defined from VertexID1 to VertecID2
# and VertexID2 to VertexID1 are considered equivalent...
sub HasEdge {
  my($This, $VertexID1, $VertexID2) = @_;

  if (!(defined($VertexID1) && defined($VertexID2))) {
    return 0;
  }

  return ($This->_HasEdge($VertexID1, $VertexID2) || $This->_HasEdge($VertexID2, $VertexID1)) ? 1 : 0;
}

# Does edge exists?
sub _HasEdge {
  my($This, $VertexID1, $VertexID2) = @_;

  if (exists $This->{Edges}->{From}->{$VertexID1}) {
    if (exists $This->{Edges}->{From}->{$VertexID1}->{$VertexID2}) {
      return 1;
    }
  }
  elsif (exists $This->{Edges}->{To}->{$VertexID2}) {
    if (exists $This->{Edges}->{To}->{$VertexID2}->{$VertexID1}) {
      return 1;
    }
  }
  return 0;
}

# Do the edges defiend by vertex pairs exist? In scalar context, return the number
# of edges found...
sub HasEdges {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    return 0;
  }
  if (@VertexIDs % 2) {
    return 0;
  }
  my($VertexID1, $VertexID2, $Index, $Status, $EdgesCount, @EdgesStatus);
  @EdgesStatus = ();
  $EdgesCount = 0;
  for ($Index = 0; $Index < $#VertexIDs; $Index += 2) {
    $VertexID1 = $VertexIDs[$Index]; $VertexID2 = $VertexIDs[$Index + 1];
    $Status = $This->HasEdge($VertexID1, $VertexID2);
    push @EdgesStatus, ($Status);
    if (defined($Status) && $Status) {
      $EdgesCount++;
    }
  }
  return wantarray ? @EdgesStatus : $EdgesCount;
}

# Get edges for a vertex ID or retrieve all the edges. In scalar context,
#  return the number of edges.
#
sub GetEdges {
  my($This, $VertexID) = @_;
  my(@VertexIDs);

  @VertexIDs = ();
  if (defined $VertexID) {
    push @VertexIDs, ($This->_GetEdgesFrom($VertexID), $This->_GetEdgesTo($VertexID))
  }
  else {
    push @VertexIDs, $This->_GetEdges();
  }
  return (wantarray ? @VertexIDs : @VertexIDs/2);
}

# Get edge starting from the vertex to its successor vertices...
sub _GetEdgesFrom {
  my($This, $VertexID1) = @_;
  my($VertexID2) = undef;

  return $This->_GetEdges($VertexID1, $VertexID2);
}

# Get edge starting from predecessors to the vertex...
sub _GetEdgesTo {
  my($This, $VertexID2) = @_;
  my($VertexID1) = undef;

  return $This->_GetEdges($VertexID1, $VertexID2);
}

# Get edges as pair of vertex IDs. Edges data can be retrieved in three
# different ways:
#
# Both vertex IDs are defined: Returns existing edge between the vertices
# Only first vertex ID defined: Returns all edges at the vertex
# Only second vertex defined: Returns all edges at the vertex
# No vertex IDs defined: Returns all edges
#
sub _GetEdges {
  my($This, $VertexID1, $VertexID2) = @_;
  my($VertexID, @VertexIDs);

  @VertexIDs = ();

  if (defined($VertexID1) && defined($VertexID2)) {
    if ($This->HasEdge($VertexID1, $VertexID2)) {
      push @VertexIDs, ($VertexID1, $VertexID2);
    }
  }
  elsif (defined($VertexID1)) {
    for $VertexID ($This->_GetNeighborsFrom($VertexID1)) {
      push @VertexIDs, $This->_GetEdges($VertexID1, $VertexID);
    }
  }
  elsif (defined($VertexID2)) {
    for $VertexID ($This->_GetNeighborsTo($VertexID2)) {
      push @VertexIDs, $This->_GetEdges($VertexID, $VertexID2);
    }
  }
  else {
    for $VertexID ($This->GetVertices()) {
      push @VertexIDs, $This->_GetEdges($VertexID);
    }
  }

  return @VertexIDs;
}

# Add edges between successive pair of vertex IDs......
sub AddPath {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    carp "Warning: ${ClassName}->AddPath: No path added: Vertices list is empty...";
    return undef;
  }
  if (@VertexIDs == 1) {
    carp "Warning: ${ClassName}->AddPath: No path added: Invalid vertices data: Input list must contain more than on vertex ID...";
    return undef;
  }
  if (!$This->HasVertices(@VertexIDs)) {
    carp "Warning: ${ClassName}->AddPath: No path added: Some of the vertex IDs don't exist in the graph...";
    return undef;
  }
  if ($This->HasPath(@VertexIDs)) {
    carp "Warning: ${ClassName}->AddPath: No path added: Path already exist in the graph...";
    return undef;
  }
  my(@PathVertexIDs);
  @PathVertexIDs =$This-> _SetupPathVertices(@VertexIDs);

  return $This->AddEdges(@PathVertexIDs);
}


# Delete edges between successive pair of vertex IDs......
sub DeletePath {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    carp "Warning: ${ClassName}->DeletePath: No path deleted: Vertices list is empty...";
    return undef;
  }
  if (@VertexIDs == 1) {
    carp "Warning: ${ClassName}->DeletePath: No path deleted: Invalid vertices data: Input list must contain more than on vertex ID...";
    return undef;
  }
  if (!$This->HasVertices(@VertexIDs)) {
    carp "Warning: ${ClassName}->DeletePath: No path deleted: Some of the vertex IDs don't exist in the graph...";
    return undef;
  }
  if (!$This->HasPath(@VertexIDs)) {
    carp "Warning: ${ClassName}->DeletePath: No path deleted: Path doesn't exist in the graph...";
    return undef;
  }
  my(@PathVertexIDs);
  @PathVertexIDs = $This->_SetupPathVertices(@VertexIDs);

  return $This->DeleteEdges(@PathVertexIDs);
}

# Does the path defiend by edges between successive pairs of vertex IDs exist?
sub HasPath {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    return 0;
  }
  if (@VertexIDs == 1) {
    return 0;
  }
  if (!$This->HasVertices(@VertexIDs)) {
    return 0;
  }
  my($Status, @PathVertexIDs);
  @PathVertexIDs = $This->_SetupPathVertices(@VertexIDs);
  $Status = ($This->HasEdges(@PathVertexIDs) == (@PathVertexIDs/2)) ? 1 : 0;

  return $Status;
}

# Setup vertices for the path to define edges between successive pair of vertex IDs...
sub _SetupPathVertices {
  my($This, @VertexIDs) = @_;
  my($VertexID1, $VertexID2, $Index, @PathVertexIDs);

  @PathVertexIDs = ();
  for $Index (0 .. ($#VertexIDs - 1)) {
    $VertexID1 = $VertexIDs[$Index];
    $VertexID2 = $VertexIDs[$Index + 1];
    push @PathVertexIDs, ($VertexID1, $VertexID2);
  }

  return @PathVertexIDs;
}

# Add edges between successive pair of vertex IDs and an additional edge from the last to
# the first ID to complete the cycle......
sub AddCycle {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    carp "Warning: ${ClassName}->AddCycle: No cycle added: Vertices list is empty...";
    return undef;
  }
  if (@VertexIDs == 1) {
    carp "Warning: ${ClassName}->AddCycle: No cycle added: Invalid vertices data: Input list must contain more than on vertex ID...";
    return undef;
  }
  if (!$This->HasVertices(@VertexIDs)) {
    carp "Warning: ${ClassName}->AddCycle: No cycle added: Some of the vertex IDs don't exist in the graph...";
    return undef;
  }
  my($FirstVertextID) = $VertexIDs[0];
  push @VertexIDs, ($FirstVertextID);

  if ($This->HasCycle(@VertexIDs)) {
    carp "Warning: ${ClassName}->AddCycle: No cycle added: Cycle already exist in the graph...";
    return undef;
  }

  return $This->AddPath(@VertexIDs);
}

# Delete edges between successive pair of vertex IDs and an additional edge from the last to
# the first ID to complete the cycle......
sub DeleteCycle {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    carp "Warning: ${ClassName}->DeleteCycle: No cycle deleted: Vertices list is empty...";
    return undef;
  }
  if (@VertexIDs == 1) {
    carp "Warning: ${ClassName}->DeleteCycle: No cycle deleted: Invalid vertices data: Input list must contain more than on vertex ID...";
    return undef;
  }
  if (!$This->HasVertices(@VertexIDs)) {
    carp "Warning: ${ClassName}->DeleteCycle: No cycle deleted: Some of the vertex IDs don't exist in the graph...";
    return undef;
  }
  if (!$This->HasCycle(@VertexIDs)) {
    carp "Warning: ${ClassName}->DeleteCycle: No cycle deleted: Cycle doesn't exist in the graph...";
    return undef;
  }

  my($FirstVertextID) = $VertexIDs[0];
  push @VertexIDs, ($FirstVertextID);

  return $This->DeletePath(@VertexIDs);
}

# Does the cycle defiend by edges between successive pairs of vertex IDs along with an additional
# edge from last to first vertex ID exist?
sub HasCycle {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    return 0;
  }
  if (@VertexIDs == 1) {
    return 0;
  }
  my($FirstVertextID) = $VertexIDs[0];
  push @VertexIDs, ($FirstVertextID);

  return $This->HasPath(@VertexIDs);
}

# Get neighbors...
sub GetNeighbors {
  my($This, $VertexID) = @_;

  if (!defined $VertexID) {
    return undef;
  }
  if (! exists $This->{Vertices}->{$VertexID}) {
    return undef;
  }

  # Get a list of unsorted vertices and sort 'em once before returning...
  #
  my($VerticesCount, $SortVertices, @VertexIDs);

  $SortVertices = 0;
  @VertexIDs = ();

  push @VertexIDs, $This->_GetNeighborsFrom($VertexID, $SortVertices);
  push @VertexIDs, $This->_GetNeighborsTo($VertexID, $SortVertices);
  $VerticesCount = @VertexIDs;

  return wantarray ? sort {$a <=> $b} @VertexIDs : $VerticesCount;
}

# Get neighbors added by defining edges from the vertex. For undirected graph, it has no
# strict meaning...
sub _GetNeighborsFrom {
  my($This, $VertexID, $SortVertices) = @_;
  my(@VertexIDs);

  $SortVertices = defined $SortVertices ? $SortVertices : 1;
  @VertexIDs = ();

  if (exists $This->{Edges}->{From}->{$VertexID}) {
    push @VertexIDs, map { $This->{Edges}->{From}->{$VertexID}->{$_} }  keys %{$This->{Edges}->{From}->{$VertexID}};
  }
  return $SortVertices ? sort {$a <=> $b} @VertexIDs : @VertexIDs;
}

# Get neighbors added by defining edges to the vertex. For undirected graphs, it has no
# strict meaning.
sub _GetNeighborsTo {
  my($This, $VertexID, $SortVertices) = @_;
  my(@VertexIDs);

  $SortVertices = defined $SortVertices ? $SortVertices : 1;
  @VertexIDs = ();

  if (exists $This->{Edges}->{To}->{$VertexID}) {
    push @VertexIDs, map { $This->{Edges}->{To}->{$VertexID}->{$_} } keys %{$This->{Edges}->{To}->{$VertexID}};
  }
  return $SortVertices ? sort {$a <=> $b} @VertexIDs : @VertexIDs;
}

# Get vertex degree...
#
sub GetDegree {
  my($This, $VertexID) = @_;

  if (!defined $VertexID) {
    return undef;
  }
  if (! exists $This->{Vertices}->{$VertexID}) {
    return undef;
  }
  my($Degree);
  $Degree = $This->_GetInDegree($VertexID) + $This->_GetOutDegree($VertexID);

  return $Degree;
}

# Get in degree added by defining edges to the vertex. For undirected graphs, it has no
# strict meaning.
#
sub _GetInDegree {
  my($This, $VertexID) = @_;
  my($Degree);

  $Degree = 0;
  if (exists $This->{Edges}->{To}->{$VertexID}) {
    $Degree = keys %{$This->{Edges}->{To}->{$VertexID}};
  }
  return $Degree;
}

# Get out degree added by defining edges from the vertex. For undirected graphs, it has no
# strict meaning.
#
sub _GetOutDegree {
  my($This, $VertexID) = @_;
  my($Degree);

  $Degree = 0;
  if (exists $This->{Edges}->{From}->{$VertexID}) {
    $Degree = keys %{$This->{Edges}->{From}->{$VertexID}};
  }
  return $Degree;
}

# Get vertex with smallest degree...
#
sub GetVertexWithSmallestDegree {
  my($This) = @_;
  my($Degree, $SmallestDegree, $SmallestDegreeVertexID, $VertexID, @VertexIDs);

  @VertexIDs = ();
  @VertexIDs = $This->GetVertices();
  if (! @VertexIDs) {
    return undef;
  }
  $SmallestDegree = 99999; $SmallestDegreeVertexID = undef;

  for $VertexID (@VertexIDs) {
    $Degree = $This->_GetInDegree($VertexID) + $This->_GetOutDegree($VertexID);
    if ($Degree < $SmallestDegree) {
      $SmallestDegree = $Degree;
      $SmallestDegreeVertexID = $VertexID;
    }
  }
  return $SmallestDegreeVertexID;
}

# Get vertex with largest degree...
#
sub GetVertexWithLargestDegree {
  my($This) = @_;
  my($Degree, $LargestDegree, $LargestDegreeVertexID, $VertexID, @VertexIDs);

  @VertexIDs = ();
  @VertexIDs = $This->GetVertices();
  if (! @VertexIDs) {
    return undef;
  }
  $LargestDegree = 0; $LargestDegreeVertexID = undef;

  for $VertexID (@VertexIDs) {
    $Degree = $This->_GetInDegree($VertexID) + $This->_GetOutDegree($VertexID);
    if ($Degree > $LargestDegree) {
      $LargestDegree = $Degree;
      $LargestDegreeVertexID = $VertexID;
    }
  }
  return $LargestDegreeVertexID;
}

# Get maximum degree in the graph...
#
sub GetMaximumDegree {
  my($This) = @_;
  my($Degree, $MaximumDegree, $VertexID, @VertexIDs);

  @VertexIDs = ();
  @VertexIDs = $This->GetVertices();
  if (! @VertexIDs) {
    return undef;
  }
  $MaximumDegree = 0;

  for $VertexID (@VertexIDs) {
    $Degree = $This->GetDegree($VertexID);
    if ($Degree > $MaximumDegree) {
      $MaximumDegree = $Degree;
    }
  }
  return $MaximumDegree;
}

# Get minimum degree in the graph...
#
sub GetMininumDegree {
  my($This) = @_;
  my($Degree, $MininumDegree, $VertexID, @VertexIDs);

  @VertexIDs = ();
  @VertexIDs = $This->GetVertices();
  if (! @VertexIDs) {
    return undef;
  }
  $MininumDegree = 99999;

  for $VertexID (@VertexIDs) {
    $Degree = $This->GetDegree($VertexID);
    if ($Degree < $MininumDegree) {
      $MininumDegree = $Degree;
    }
  }
  return $MininumDegree;
}

# Is it a isolated vertex?
#
sub IsIsolatedVertex {
  my($This, $VertexID) = @_;

  if (!defined $VertexID) {
    return undef;
  }
  if (! exists $This->{Vertices}->{$VertexID}) {
    return undef;
  }
  return ($This->GetDegree() == 0) ? 1 : 0;
}

# Get all isolated vertices...
#
sub GetIsolatedVertices {
  my($This) = @_;

  return $This->GetVerticesWithDegreeLessThan(1);
}

# Is it a leaf vertex?
#
sub IsLeafVertex {
  my($This, $VertexID) = @_;

  if (!defined $VertexID) {
    return undef;
  }
  if (! exists $This->{Vertices}->{$VertexID}) {
    return undef;
  }
  return ($This->GetDegree() == 1) ? 1 : 0;
}

# Get all leaf vertices...
#
sub GetLeafVertices {
  my($This) = @_;

  return $This->GetVerticesWithDegreeLessThan(2);
}

# Get vertices  with degree less than a specified value...
#
sub GetVerticesWithDegreeLessThan {
  my($This, $SpecifiedDegree) = @_;
  my($Degree, $VertexID, @VertexIDs, @FilteredVertexIDs);

  @VertexIDs = (); @FilteredVertexIDs = ();

  @VertexIDs = $This->GetVertices();
  if (! @VertexIDs) {
    return @FilteredVertexIDs;
  }

  for $VertexID (@VertexIDs) {
    $Degree = $This->_GetInDegree($VertexID) + $This->_GetOutDegree($VertexID);
    if ($Degree < $SpecifiedDegree) {
      push @FilteredVertexIDs, $VertexID;
    }
  }
  return @FilteredVertexIDs;
}

# Set a property for graph...
sub SetGraphProperty {
  my($This, $Name, $Value) = @_;

  if (!(defined($Name) && defined($Value))) {
    carp "Warning: ${ClassName}->SetGraphProperty: Didn't set property: Both property name and value should be specified...";
    return undef;
  }
  if (exists $This->{Properties}->{Graph}->{$Name}) {
    carp "Warning: ${ClassName}->SetGraphProperty: Didn't set property $Name: Already exists in the graph...";
    return undef;
  }

  $This->{Properties}->{Graph}->{$Name} = $Value;

  return $This;
}

# Set a properties for graph...
sub SetGraphProperties {
  my($This, %NamesAndValues) = @_;

  if (!(keys %NamesAndValues)) {
    carp "Warning: ${ClassName}->SetGraphProperties: Didn't set properties: Names and values list is empty...";
    return undef;
  }

  my($Name, $Value);
  while (($Name, $Value) = each  %NamesAndValues) {
    $This->SetGraphProperty($Name, $Value);
  }

  return $This;
}

# Get a property value for graph...
sub GetGraphProperty {
  my($This, $Name) = @_;

  if (!$This->HasGraphProperty($Name)) {
    return undef;
  }

  return $This->{Properties}->{Graph}->{$Name};
}

# Get all poperty name/value pairs for graph...
sub GetGraphProperties {
  my($This) = @_;

  return %{$This->{Properties}->{Graph}};
}

# Delete a property for graph...
sub DeleteGraphProperty {
  my($This, $Name) = @_;

  if (!defined $Name) {
    carp "Warning: ${ClassName}->DeleteGraphProperty: Didn't delete graph property: Name must be specified...";
    return undef;
  }
  if (!$This->HasGraphProperty($Name)) {
    carp "Warning: ${ClassName}-> DeleteGraphProperty: Didn't delete graph property $Name: Property doesn't exist...";
    return undef;
  }
  delete $This->{Properties}->{Graph}->{$Name};

  return $This;
}

# Delete graph properites...
sub DeleteGraphProperties {
  my($This) = @_;

  %{$This->{Properties}->{Graph}} = ();

  return $This;
}


# Is this property associated with graph?
sub HasGraphProperty {
  my($This, $Name) = @_;

  if (!defined $Name) {
    return 0;
  }
  return (exists $This->{Properties}->{Graph}->{$Name}) ? 1 : 0;
}

# Set a property for vertex...
sub SetVertexProperty {
  my($This, $Name, $Value, $VertexID) = @_;

  if (!(defined($Name) && defined($Value) && defined($VertexID))) {
    carp "Warning: ${ClassName}->SetVertexProperty: Didn't set property: Property name, value and vertexID should be specified...";
    return undef;
  }
  if (!$This->HasVertex($VertexID)) {
    carp "Warning: ${ClassName}->SetVertexProperty: Didn't set property $Name for vertex $VertexID: Vertex doesn't exist...";
    return undef;
  }
  if ($This->HasVertexProperty($Name, $VertexID)) {
    carp "Warning: ${ClassName}->SetVertexProperty: Didn't set property $Name for vertex $VertexID: Property already exists...";
    return undef;
  }

  $This->_SetVertexProperty($Name, $Value, $VertexID);

  return $This;
}

# Update a property for vertex...
sub UpdateVertexProperty {
  my($This, $Name, $Value, $VertexID) = @_;

  if (!(defined($Name) && defined($Value) && defined($VertexID))) {
    carp "Warning: ${ClassName}->UpdateVextexProperty: Didn't update property: Property name, value and vertexID should be specified...";
    return undef;
  }
  if (!$This->HasVertex($VertexID)) {
    carp "Warning: ${ClassName}->UpdateVextexProperty: Didn't update property $Name for vertex $VertexID: Vertex doesn't exist...";
    return undef;
  }
  if (!$This->HasVertexProperty($Name, $VertexID)) {
    carp "Warning: ${ClassName}->UpdateVextexProperty: Didn't update property $Name for vertex $VertexID: Property doesn't exists...";
    return undef;
  }
  $This->_SetVertexProperty($Name, $Value, $VertexID);

  return $This;
}

# Set a vextex property...
sub _SetVertexProperty {
  my($This, $Name, $Value, $VertexID) = @_;

  if (!exists $This->{Properties}->{Vertices}->{$VertexID}) {
    %{$This->{Properties}->{Vertices}->{$VertexID}} = ();
  }
  $This->{Properties}->{Vertices}->{$VertexID}->{$Name} = $Value;

  return $This;
}

# Set a properties for vertex..
sub SetVertexProperties {
  my($This, $VertexID, @NamesAndValues) = @_;

  if (!defined $VertexID) {
    carp "Warning: ${ClassName}->SetVertexProperties: Didn't set property: VertexID should be specified...";
    return undef;
  }
  if (!@NamesAndValues) {
    carp "Warning: ${ClassName}->SetVertexProperties: Didn't set property: Names and values list is empty...";
    return undef;
  }
  if (@NamesAndValues % 2) {
    carp "Warning: ${ClassName}->SetVertexProperties: Didn't set property: Invalid property name and values IDs data: Input list must contain even number of values...";
    return undef;
  }

  my($Name, $Value, $Index);
  for ($Index = 0; $Index < $#NamesAndValues; $Index += 2) {
    $Name = $NamesAndValues[$Index]; $Value = $NamesAndValues[$Index + 1];
    $This->SetVertexProperty($Name, $Value, $VertexID);
  }

  return $This;
}


# Set a property for vertices...
sub SetVerticesProperty {
  my($This, $Name, @ValuesAndVertexIDs) = @_;

  if (!defined $Name) {
    carp "Warning: ${ClassName}->SetVerticesProperty: Didn't set property: Property name should be specified...";
    return undef;
  }
  if (!@ValuesAndVertexIDs) {
    carp "Warning: ${ClassName}->SetVerticesProperty: Didn't set property: Values and vertex IDs list is empty...";
    return undef;
  }
  if (@ValuesAndVertexIDs % 2) {
    carp "Warning: ${ClassName}->SetVerticesProperty: Didn't set property: Invalid property values and vertex IDs data: Input list must contain even number of values...";
    return undef;
  }

  my($Value, $VertexID, $Index);
  for ($Index = 0; $Index < $#ValuesAndVertexIDs; $Index += 2) {
    $Value = $ValuesAndVertexIDs[$Index]; $VertexID = $ValuesAndVertexIDs[$Index + 1];
    $This->SetVertexProperty($Name, $Value, $VertexID);
  }

  return $This;
}

# Get a property value for vertex...
sub GetVertexProperty {
  my($This, $Name, $VertexID) = @_;

  if (!$This->HasVertexProperty($Name, $VertexID)) {
    return undef;
  }

  return $This->{Properties}->{Vertices}->{$VertexID}->{$Name};
}

# Get a property values for vertices...
sub GetVerticesProperty {
  my($This, $Name, @VertexIDs) = @_;
  my($ValuesCount, @Values);

  if (!@VertexIDs) {
    @VertexIDs = ();
    @VertexIDs = $This->GetVertices();
  }
  @Values = ();
  @Values = map { $This->GetVertexProperty($Name, $_ ) } @VertexIDs;
  $ValuesCount = grep { defined($_) } @Values;

  return wantarray ? @Values : $ValuesCount;
}

# Get all property name/values pairs for vertex...
sub GetVertexProperties {
  my($This, $VertexID) = @_;

  if (!$This->HasVertex($VertexID)) {
    return ();
  }

  if (exists $This->{Properties}->{Vertices}->{$VertexID}) {
    return %{$This->{Properties}->{Vertices}->{$VertexID}};
  }
  else {
    return ();
  }
}


# Delete a property for vertex...
sub DeleteVertexProperty {
  my($This, $Name, $VertexID) = @_;

  if (!(defined($Name) && defined($VertexID))) {
    carp "Warning: ${ClassName}->DeleteVertexProperty: Didn't delete property: Property name and vertex ID must be defined...";
    return undef;
  }
  if (!$This->HasVertexProperty($Name, $VertexID)) {
    carp "Warning: ${ClassName}->DeleteVertexProperty: Didn't delete property $Name for vertex $VertexID: Vertex property doesn't exist...";
    return undef;
  }
  $This->_DeleteVertexProperty($VertexID, $Name);

  return $This;
}

# Delete vertex property or all properties...
sub _DeleteVertexProperty {
  my($This, $VertexID, $Name) = @_;

  if (exists $This->{Properties}->{Vertices}->{$VertexID}) {
    if (defined $Name) {
      if (exists $This->{Properties}->{Vertices}->{$VertexID}->{$Name}) {
	delete $This->{Properties}->{Vertices}->{$VertexID}->{$Name};
      }
    }
    else {
      %{$This->{Properties}->{Vertices}->{$VertexID}} = ();
    }
    if (! keys %{$This->{Properties}->{Vertices}->{$VertexID}}) {
      delete $This->{Properties}->{Vertices}->{$VertexID};
    }
  }
  return $This;
}

# Delete a property for specified vertices or all the vertices...
sub DeleteVerticesProperty {
  my($This, $Name, @VertexIDs) = @_;

  if (!defined $Name) {
    carp "Warning: ${ClassName}->DeleteVerticesProperty: Didn't delete property: Property name should be specified...";
    return undef;
  }
  if (!@VertexIDs) {
    @VertexIDs = ();
    @VertexIDs = $This->GetVertices();
  }
  my($VertexID);
  for $VertexID (@VertexIDs) {
    $This->DeleteVertexProperty($Name, $VertexID);
  }

  return $This;
}

# Delete all properities for vertex...
sub DeleteVertexProperties {
  my($This, $VertexID) = @_;

  if (!defined $VertexID) {
    carp "Warning: ${ClassName}->DeleteVertexProperties: Didn't delete properties: Vertex ID must be defined...";
    return undef;
  }
  $This->_DeleteVertexProperty($VertexID);

  return $This;
}

# Is this property associated with vertex?
sub HasVertexProperty {
  my($This, $Name, $VertexID) = @_;

  if (!(defined($Name) && defined($VertexID))) {
    return 0;
  }

  if (exists $This->{Properties}->{Vertices}->{$VertexID}) {
    if (exists $This->{Properties}->{Vertices}->{$VertexID}->{$Name}) {
      return 1;
    }
  }
  return 0;
}

# Set a property for edge...
sub SetEdgeProperty {
  my($This, $Name, $Value, $VertexID1, $VertexID2) = @_;

  if (!(defined($Name) && defined($Value) && defined($VertexID1) && defined($VertexID2))) {
    carp "Warning: ${ClassName}->SetEdgeProperty: Didn't set property: Property name, value, vertexID1 and vertexID2 should be specified...";
    return undef;
  }
  if (!$This->HasEdge($VertexID1, $VertexID2)) {
    carp "Warning: ${ClassName}->SetEdgeProperty: Didn't set property $Name for edge between vertices $VertexID1 and $VertexID2: Edge doesn't exist...";
    return undef;
  }
  if ($This->HasEdgeProperty($Name, $VertexID1, $VertexID2)) {
    carp "Warning: ${ClassName}->SetEdgeProperty: Didn't set property $Name for edge between vertices $VertexID1 and $VertexID2: Edge property already exists...";
    return undef;
  }
  $This->_SetEdgeProperty($Name, $Value, $VertexID1, $VertexID2);
  $This->_SetEdgeProperty($Name, $Value, $VertexID2, $VertexID1);

  return $This;
}

# Update a property for edge...
sub UpdateEdgeProperty {
  my($This, $Name, $Value, $VertexID1, $VertexID2) = @_;

  if (!(defined($Name) && defined($Value) && defined($VertexID1) && defined($VertexID2))) {
    carp "Warning: ${ClassName}->UpdateEdgeProperty: Didn't update property: Property name, value, vertexID1 and vertexID2 should be specified...";
    return undef;
  }
  if (!$This->HasEdge($VertexID1, $VertexID2)) {
    carp "Warning: ${ClassName}->UpdateEdgeProperty: Didn't update property $Name for edge between vertices $VertexID1 and $VertexID2: Edge doesn't exist...";
    return undef;
  }
  if (!$This->HasEdgeProperty($Name, $VertexID1, $VertexID2)) {
    carp "Warning: ${ClassName}->UpdateEdgeProperty: Didn't update property $Name for edge between vertices $VertexID1 and $VertexID2: Edge property doesn't exists...";
    return undef;
  }
  $This->_SetEdgeProperty($Name, $Value, $VertexID1, $VertexID2);
  $This->_SetEdgeProperty($Name, $Value, $VertexID2, $VertexID1);

  return $This;
}
# Set a property for edges...
sub SetEdgesProperty {
  my($This, $Name, @ValuesAndVertexIDs) = @_;

  if (!defined $Name) {
    carp "Warning: ${ClassName}->SetEdgesProperty: Didn't set property: Property name should be specified...";
    return undef;
  }
  if (!@ValuesAndVertexIDs) {
    carp "Warning: ${ClassName}->SetEdgesProperty: Didn't set property: Values and vertex IDs list is empty...";
    return undef;
  }
  if (@ValuesAndVertexIDs % 3) {
    carp "Warning: ${ClassName}->SetEdgesProperty: Didn't set property: Invalid property values and vertex IDs data: Input list must contain triplets of values...";
    return undef;
  }

  my($Value, $VertexID1, $VertexID2, $Index);
  for ($Index = 0; $Index < $#ValuesAndVertexIDs; $Index += 3) {
    $Value = $ValuesAndVertexIDs[$Index];
    $VertexID1 = $ValuesAndVertexIDs[$Index + 1]; $VertexID2 = $ValuesAndVertexIDs[$Index + 2];
    $This->SetEdgeProperty($Name, $Value, $VertexID1, $VertexID2);
  }

  return $This;
}

# Set a properties for a edge...
sub SetEdgeProperties {
  my($This, $VertexID1, $VertexID2, @NamesAndValues) = @_;

  if (!(defined($VertexID1) && defined($VertexID2))) {
    carp "Warning: ${ClassName}->SetEdgeProperties: Didn't set property: Both vertexID1 and vertexID2 should be specified...";
    return undef;
  }
  if (!@NamesAndValues) {
    carp "Warning: ${ClassName}->SetEdgeProperties: Didn't set property: Property name and values ist is empty...";
    return undef;
  }
  if (@NamesAndValues % 2) {
    carp "Warning: ${ClassName}->SetEdgeProperties: Didn't set property: Invalid property name and values data: Input list must contain triplets of values...";
    return undef;
  }

  my($Name, $Value, $Index);
  for ($Index = 0; $Index < $#NamesAndValues; $Index += 2) {
    $Name = $NamesAndValues[$Index];
    $Value = $NamesAndValues[$Index + 1];
    $This->SetEdgeProperty($Name, $Value, $VertexID1, $VertexID2);
  }

  return $This;
}

# Set edge property...
sub _SetEdgeProperty {
  my($This, $Name, $Value, $VertexID1, $VertexID2) = @_;

  if (!exists $This->{Properties}->{Edges}->{$VertexID1}) {
    %{$This->{Properties}->{Edges}->{$VertexID1}} = ();
  }
  if (!exists $This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2}) {
    %{$This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2}} = ();
  }
  $This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2}->{$Name} = $Value;

  return $This;
}

# Get a property value for edge...
sub GetEdgeProperty {
  my($This, $Name, $VertexID1, $VertexID2) = @_;

  if (!$This->HasEdgeProperty($Name, $VertexID1, $VertexID2)) {
    return undef;
  }
  return $This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2}->{$Name};
}

# Get a property values for edges...
sub GetEdgesProperty {
  my($This, $Name, @VertexIDs) = @_;

  if (!@VertexIDs) {
    @VertexIDs = ();
    @VertexIDs = $This->GetEdges();
  }
  if (@VertexIDs % 2) {
    return undef;
  }

  my($VertexID1, $VertexID2, $Index, $ValuesCount, @Values);
  @Values = ();
  for ($Index = 0; $Index < $#VertexIDs; $Index += 2) {
    $VertexID1 = $VertexIDs[$Index]; $VertexID2 = $VertexIDs[$Index + 1];
    push @Values, $This->GetEdgeProperty($Name, $VertexID1, $VertexID2);
  }
  $ValuesCount = grep { defined($_) } @Values;

  return wantarray ? @Values : $ValuesCount;
}

# Get all property name/value paries for edge...
sub GetEdgeProperties {
  my($This, $VertexID1, $VertexID2) = @_;

  if (!(defined($VertexID1) && defined($VertexID2))) {
    return ();
  }
  if (!exists $This->{Properties}->{Edges}->{$VertexID1}) {
    return ();
  }
  if (!exists $This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2}) {
    return ();
  }
  return %{$This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2}};
}

# Delete a property for edge...
sub DeleteEdgeProperty {
  my($This, $Name, $VertexID1, $VertexID2) = @_;

  if (!(defined($Name) && defined($VertexID1) && defined($VertexID2))) {
    carp "Warning: ${ClassName}->DeleteEdgeProperty: Didn't delete property: Property name, vertexID1 and vertexID2 should be specified...";
    return undef;
  }
  if (!$This->HasEdgeProperty($Name, $VertexID1, $VertexID2)) {
    carp "Warning: ${ClassName}->DeleteEdgeProperty: Didn't delete property $Name for edge between vertices $VertexID1 and $VertexID2: Edge property doesn't exist...";
    return undef;
  }
  $This->_DeleteEdgeProperty($VertexID1, $VertexID2, $Name);
  $This->_DeleteEdgeProperty($VertexID2, $VertexID1, $Name);

  return $This;
}

# Delete a property for edges...
sub DeleteEdgesProperty {
  my($This, $Name, @VertexIDs) = @_;

  if (!defined $Name) {
    carp "Warning: ${ClassName}->DeleteEdgesProperty: Didn't delete property: Property name should be specified...";
    return undef;
  }
  if (!@VertexIDs) {
    @VertexIDs = ();
    @VertexIDs = $This->GetEdges();
  }
  if (@VertexIDs % 2) {
    carp "Warning: ${ClassName}->DeleteEdgesProperty: Didn't set property: Invalid property values and vertex IDs data: Input list must contain even number of values...";
    return undef;
  }
  my($Index, $VertexID1, $VertexID2);
  for ($Index = 0; $Index < $#VertexIDs; $Index += 2) {
    $VertexID1 = $VertexIDs[$Index]; $VertexID2 = $VertexIDs[$Index + 1];
    $This->DeleteEdgeProperty($Name, $VertexID1, $VertexID2);
  }

  return $This;
}

# Delete all properties for edge...
sub DeleteEdgeProperties {
  my($This, $VertexID1, $VertexID2) = @_;

  if (!(defined($VertexID1) && defined($VertexID2))) {
    carp "Warning: ${ClassName}->DeleteEdgeProperties: Didn't delete property: VertexID1 and vertexID2 should be specified...";
    return undef;
  }
  $This->_DeleteEdgeProperty($VertexID1, $VertexID2);
  $This->_DeleteEdgeProperty($VertexID2, $VertexID1);

  return $This;
}

# Delete properties for edges...
sub DeleteEdgesProperties {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    @VertexIDs = ();
    @VertexIDs = $This->GetEdges();
  }
  if (@VertexIDs % 2) {
    carp "Warning: ${ClassName}->DeleteEdgesProperty: Didn't set property: Invalid property values and vertex IDs data: Input list must contain even number of values...";
    return undef;
  }
  my($Index, $VertexID1, $VertexID2);
  for ($Index = 0; $Index < $#VertexIDs; $Index += 2) {
    $VertexID1 = $VertexIDs[$Index]; $VertexID2 = $VertexIDs[$Index + 1];
    $This->DeleteEdgeProperties($VertexID1, $VertexID2);
  }
  return $This;
}


# Delete a specific edge property or all edge properties...
sub _DeleteEdgeProperty {
  my($This, $VertexID1, $VertexID2, $Name) = @_;

  if (exists $This->{Properties}->{Edges}->{$VertexID1}) {
    if (exists $This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2}) {
      if (defined $Name) {
	if (exists $This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2}->{$Name}) {
	  delete $This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2}->{$Name};
	}
      }
      else {
	%{$This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2}} = ();
      }
      if (! keys %{$This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2}}) {
	delete $This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2};
      }
    }
    if (! keys %{$This->{Properties}->{Edges}->{$VertexID1}}) {
      delete $This->{Properties}->{Edges}->{$VertexID1};
    }
  }

  return $This;
}

# Is this property associated with edge?
sub HasEdgeProperty {
  my($This, $Name, $VertexID1, $VertexID2) = @_;

  if (!(defined($Name) && defined($VertexID1) && defined($VertexID2))) {
    return 0;
  }
  my($Status);

  $Status = ($This->_HasEdgeProperty($Name, $VertexID1, $VertexID2) || $This->_HasEdgeProperty($Name, $VertexID2, $VertexID1)) ? 1 : 0;

  return $Status;
}

# Does edge proprty exists?
sub _HasEdgeProperty {
  my($This, $Name, $VertexID1, $VertexID2) = @_;

  if (exists $This->{Properties}->{Edges}->{$VertexID1}) {
    if (exists $This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2}) {
      if (exists $This->{Properties}->{Edges}->{$VertexID1}->{$VertexID2}->{$Name}) {
	return 1;
      }
    }
  }
  return 0;
}

# Is it a graph object?
sub IsGraph ($) {
  my($Object) = @_;

  return _IsGraph($Object);
}

# Return a string containg vertices, edges and other properties...
sub StringifyGraph {
  my($This) = @_;
  my($GraphString);

  $GraphString = 'Graph:' . $This->StringifyVerticesAndEdges() . '; ' . $This->StringifyProperties();

  return $GraphString;
}

# Return a string containg vertices, edges and other properties...
sub StringifyProperties {
  my($This) = @_;
  my($PropertiesString);

  $PropertiesString = $This->StringifyGraphProperties() . '; ' . $This->StringifyVerticesProperties(). '; ' . $This->StringifyEdgesProperties();

  return $PropertiesString;
}

# Return a string containg vertices and edges...
sub StringifyVerticesAndEdges {
  my($This) = @_;
  my($GraphString, $Index, $VertexID, $VertexID1, $VertexID2, $Count, @EdgeVertexIDs, @VertexIDs);

  # Get vertices and edges...
  $GraphString = '';
  @VertexIDs = ();
  @VertexIDs = $This->GetVertices();
  $Count = 0;
  for $VertexID (@VertexIDs) {
    $Count++;
    @EdgeVertexIDs = ();
    @EdgeVertexIDs = $This->_GetEdgesFrom($VertexID);
    if (@EdgeVertexIDs) {
      for ($Index = 0; $Index < $#EdgeVertexIDs; $Index += 2) {
	$VertexID1 = $EdgeVertexIDs[$Index]; $VertexID2 = $EdgeVertexIDs[$Index + 1];
	$GraphString .= " ${VertexID1}-${VertexID2}";
      }
    }
    else {
      $GraphString .= " ${VertexID}";
    }
  }
  if (!$Count) {
    $GraphString = 'Graph: None';
  }
  return $GraphString;
}

# Return a string containg graph properties...
sub StringifyGraphProperties {
  my($This) = @_;
  my($Name, $Value, $Count, $GraphPropertyString, %GraphProperties);

  $GraphPropertyString = "GraphProperties: ";
  %GraphProperties = ();
  %GraphProperties = $This->GetGraphProperties();
  if (keys %GraphProperties) {
    for $Name (sort keys %GraphProperties) {
      $Value = $GraphProperties{$Name};
      if (ref($Value) =~ /^ARRAY/) {
	if (@{$Value}) {
	  $GraphPropertyString .= " ${Name}=(Count: " . scalar @{$Value} . "; " . join(', ', @{$Value}) .  ")";
	}
	else {
	  $GraphPropertyString .= " ${Name}=None";
	}
      }
      else {
	$GraphPropertyString .= " ${Name}=${Value}";
      }
    }
  }
  else {
    $GraphPropertyString .= " None";
  }
  return $GraphPropertyString;
}

# Return a string containg vertices  properties...
sub StringifyVerticesProperties {
  my($This) = @_;
  my($Name, $Value, $Count, $VertexPropertyString, $VertexID, @VertexIDs, %VertexProperties);

  @VertexIDs = ();
  @VertexIDs = $This->GetVertices();
  $Count = 0;
  $VertexPropertyString = "VertexProperties:";
  for $VertexID (@VertexIDs) {
    %VertexProperties = ();
    %VertexProperties = $This->GetVertexProperties($VertexID);
    if (keys %VertexProperties) {
      $Count++;
      $VertexPropertyString .= " <Vertex ${VertexID}: ";
      for $Name (sort keys %VertexProperties) {
	$Value = $VertexProperties{$Name};
	if (ref($Value) =~ /^ARRAY/) {
	  if (@{$Value}) {
	    $VertexPropertyString .= " ${Name}=(" . join(', ', @{$Value}) .  ")";
	  }
	  else {
	    $VertexPropertyString .= " ${Name}=None";
	  }
	}
	else {
	  $VertexPropertyString .= " ${Name}=${Value}";
	}
      }
      $VertexPropertyString .= ">";
    }
  }
  if (!$Count) {
    $VertexPropertyString = "VertexProperties: None";
  }
  return $VertexPropertyString;
}

# Return a string containg edges properties...
sub StringifyEdgesProperties {
  my($This) = @_;
  my($Name, $Value, $Index, $EdgePropertyString, $Count, $VertexID, $VertexID1, $VertexID2, @EdgesVertexIDs, %EdgeProperties);

  @EdgesVertexIDs = ();
  @EdgesVertexIDs = $This->GetEdges();
  $Count = 0;
  $EdgePropertyString = "EdgeProperties:";
  for ($Index = 0; $Index < $#EdgesVertexIDs; $Index += 2) {
    $VertexID1 = $EdgesVertexIDs[$Index]; $VertexID2 = $EdgesVertexIDs[$Index + 1];
    %EdgeProperties = ();
    %EdgeProperties = $This->GetEdgeProperties($VertexID1, $VertexID2);
    if (keys %EdgeProperties) {
      $Count++;
      $EdgePropertyString .= " <Edge ${VertexID1}-${VertexID2}:";
      for $Name (sort keys %EdgeProperties) {
	$Value = $EdgeProperties{$Name};
	if (ref($Value) =~ /^ARRAY/) {
	  if (@{$Value}) {
	    $EdgePropertyString .= " ${Name}=(" . join(', ', @{$Value}) .  ")";
	  }
	  else {
	    $EdgePropertyString .= " ${Name}=None";
	  }
	}
	else {
	  $EdgePropertyString .= " ${Name}=${Value}";
	}
      }
      $EdgePropertyString .= ">";
    }
  }
  if (!$Count) {
    $EdgePropertyString = "EdgeProperties: None";
  }

  return $EdgePropertyString;
}

# Is it a graph object?
sub _IsGraph {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Copy graph and its associated data using Storable::dclone and return a new graph...
#
sub Copy {
  my($This) = @_;
  my($NewGraph);

  $NewGraph = Storable::dclone($This);

  return $NewGraph;
}

# Copy vertrices and edges from This graph to NewGraph specified...
#
sub CopyVerticesAndEdges {
  my($This, $NewGraph) = @_;

  # Copy vertices and edges...
  my(@Vertices, @Edges);
  @Vertices = $This->GetVertices();
  if (@Vertices) {
    $NewGraph->AddVertices(@Vertices);
  }
  @Edges = $This->GetEdges();
  if (@Edges) {
    $NewGraph->AddEdges(@Edges);
  }

  return $NewGraph;
}

# Copy properties of vertices from This graph to NewGraph specified...
#
sub CopyVerticesProperties {
  my($This, $NewGraph) = @_;

  my($VertexID, @VertexIDs, %VertexProperties);
  @VertexIDs = ();
  @VertexIDs = $This->GetVertices();
  for $VertexID (@VertexIDs) {
    %VertexProperties = (); %VertexProperties = $This->GetVertexProperties($VertexID);
    if (keys %VertexProperties) {
      $NewGraph->SetVertexProperties($VertexID, %VertexProperties);
    }
  }
  return $NewGraph;
}

# Copy properties of edges from This graph to NewGraph specified...
#
sub CopyEdgesProperties {
  my($This, $NewGraph) = @_;

  my($Index, $VertexID1, $VertexID2, @EdgesVertexIDs, %EdgeProperties);
  @EdgesVertexIDs = ();
  @EdgesVertexIDs = $This->GetEdges();
  for ($Index = 0; $Index < $#EdgesVertexIDs; $Index += 2) {
    $VertexID1 = $EdgesVertexIDs[$Index]; $VertexID2 = $EdgesVertexIDs[$Index + 1];
    %EdgeProperties = (); %EdgeProperties = $This->GetEdgeProperties($VertexID1, $VertexID2);
    if (keys %EdgeProperties) {
      $NewGraph->SetEdgeProperties($VertexID1, $VertexID2, %EdgeProperties);
    }
  }
  return $NewGraph;
}

# Detect cycles and associate 'em to graph as graph property...
#
# Note:
#   . CyclesDetection class detects all cycles in the graph and filters 'em to find
#     independent cycles.
#   . All cycles related methods in the graph operate on active cyclic paths. By default,
#     active cyclic paths correspond to independent cycles. This behavior can be changed
#     using SetActiveCyclicPaths method.
#   . For topologically complex graphs containing large number of cycles, DetectCycles method
#     implemented in CyclesDetection can time out in which case no cycles are detected or
#     assigned.
#
sub DetectCycles {
  my($This) = @_;
  my($CyclesDetection);

  # Delete existing graph cycles...
  $This->_DeleteCyclesAssociatedWithGraph();

  # Detect cycles...
  $CyclesDetection = new Graph::CyclesDetection($This);
  if (!$CyclesDetection->DetectCycles()) {
    # Cycles detection didn't finish...
    return undef;
  }

  # Get cycles and associate 'em to graph as properties...
  my(@AllCyclicPaths, @IndependentCyclicPaths);
  @AllCyclicPaths = $CyclesDetection->GetAllCyclicPaths();
  @IndependentCyclicPaths = $CyclesDetection->GetIndependentCyclicPaths();

  $This->SetGraphProperty('ActiveCyclicPaths', \@IndependentCyclicPaths);
  $This->SetGraphProperty('AllCyclicPaths', \@AllCyclicPaths);
  $This->SetGraphProperty('IndependentCyclicPaths', \@IndependentCyclicPaths);

  # Map cycles information to vertices and edges; identify fused cycles as well...
  return $This->_ProcessDetectedCycles();
}

# Delete any cycle properties assigned to graph, vertices and edges during detect cycles operation...
#
sub ClearCycles {
  my($This) = @_;

  # Delete cycle properties associated with graph...
  $This->_DeleteCyclesAssociatedWithGraph();
  $This->_DeleteFusedCyclesAssociatedWithGraph();

  # Delete cycle properties associated with vertices and edges...
  $This->_DeleteCyclesAssociatedWithVertices();
  $This->_DeleteCyclesAssociatedWithEdges();

  return $This;
}

# Setup cyclic paths to use during all cycle related methods. Possible values:
# Independent or All. Default is to use Independent cyclic paths.
#
sub SetActiveCyclicPaths {
  my($This, $CyclicPathsType) = @_;

  if (!defined $CyclicPathsType) {
    carp "Warning: ${ClassName}->SetActiveCyclicPaths: Didn't set active cyclic path: Cyclic path must be specified...";
    return undef;
  }
  if ($CyclicPathsType !~ /^(Independent|All)$/i) {
    carp "Warning: ${ClassName}->SetActiveCyclicPaths: Didn't set active cyclic path: Specified path type, $CyclicPathsType, is not valid. Supported valeus: Independent or All...";
    return undef;
  }
  if (!$This->HasGraphProperty('ActiveCyclicPaths')) {
    carp "Warning: ${ClassName}->SetActiveCyclicPaths: Didn't set active cyclic path: Cycles haven't been detected yet...";
    return undef;
  }
  $This->DeleteGraphProperty('ActiveCyclicPaths');

  my($ActiveCyclicPathsRef);
  if ($CyclicPathsType =~ /^Independent$/i) {
    $ActiveCyclicPathsRef = $This->GetGraphProperty('IndependentCyclicPaths');
  }
  elsif ($CyclicPathsType =~ /^All$/i) {
    $ActiveCyclicPathsRef = $This->GetGraphProperty('AllCyclicPaths');
  }
  $This->SetGraphProperty('ActiveCyclicPaths', $ActiveCyclicPathsRef);

  # Map cycles information to vertices and edges; identify fused cycles as well...
  $This->_ProcessDetectedCycles();

  return $This;
}

# Assign cycles information on to vertices and edges as vertex edge properties properties;
# identify fused cycles as well...
#
sub _ProcessDetectedCycles {
  my($This) = @_;

  $This->_AssociateCyclesWithVertices();
  $This->_AssociateCyclesWithEdgesAndIdentifyFusedCycles();

  return $This;
}

# Associate cycles information to vertices as vertex properties...
#
sub _AssociateCyclesWithVertices {
  my($This) = @_;

  # Clear up any exisiting properties...
  $This->_DeleteCyclesAssociatedWithVertices();

  # Collects CyclicPaths for each vertex...
  my($VertexID, $ActiveCyclicPath, $ActiveCyclicPathsRef, @CyclicPathVertexIDs, %VertexIDToCylicPaths);

  %VertexIDToCylicPaths = ();
  $ActiveCyclicPathsRef = $This->GetGraphProperty('ActiveCyclicPaths');

  if (!@{$ActiveCyclicPathsRef}) {
    # No cycles out there...
    return $This;
  }

  for $ActiveCyclicPath (@{$ActiveCyclicPathsRef}) {
    @CyclicPathVertexIDs = ();
    @CyclicPathVertexIDs = $ActiveCyclicPath->GetVertices();
    # Take out end vertex: It's same as start vertex for cyclic path...
    pop @CyclicPathVertexIDs;
    for $VertexID (@CyclicPathVertexIDs) {
      if (!exists $VertexIDToCylicPaths{$VertexID}) {
	@{$VertexIDToCylicPaths{$VertexID}} = ();
      }
      push @{$VertexIDToCylicPaths{$VertexID}}, $ActiveCyclicPath;
    }
  }

  # Associate CyclicPaths to vertices...
  for $VertexID (keys %VertexIDToCylicPaths) {
    $This->SetVertexProperty('ActiveCyclicPaths', \@{$VertexIDToCylicPaths{$VertexID}}, $VertexID);
  }
  return $This;
}

# Associate cycles information to edges as edge properties and identify fused
# cycles...
#
sub _AssociateCyclesWithEdgesAndIdentifyFusedCycles {
  my($This) = @_;

  # Delete existing cycles...
  $This->_DeleteCyclesAssociatedWithEdges();
  $This->_DeleteFusedCyclesAssociatedWithGraph();

  # Collect cyclic paths for each edge...
  my($Index, $VertexID1, $VertexID2, $ActiveCyclicPath, $ActiveCyclicPathsRef, $EdgeID, $EdgeIDDelimiter, $CyclesWithCommonEdgesPresent, @CyclicPathEdgeVertexIDs, %EdgeIDToCylicPaths);

  %EdgeIDToCylicPaths = ();
  $EdgeIDDelimiter = "~";
  $ActiveCyclicPathsRef = $This->GetGraphProperty('ActiveCyclicPaths');

  if (!@{$ActiveCyclicPathsRef}) {
    # No cycles out there...
    return $This;
  }

  $CyclesWithCommonEdgesPresent = 0;
  for $ActiveCyclicPath (@{$ActiveCyclicPathsRef}) {
    @CyclicPathEdgeVertexIDs = ();
    @CyclicPathEdgeVertexIDs = $ActiveCyclicPath->GetEdges();
    for ($Index = 0; $Index < $#CyclicPathEdgeVertexIDs; $Index += 2) {
      $VertexID1 = $CyclicPathEdgeVertexIDs[$Index]; $VertexID2 = $CyclicPathEdgeVertexIDs[$Index + 1];
      $EdgeID = ($VertexID1 < $VertexID2) ? "${VertexID1}${EdgeIDDelimiter}${VertexID2}" : "${VertexID2}${EdgeIDDelimiter}${VertexID1}";
      if (exists $EdgeIDToCylicPaths{$EdgeID}) {
	# A common edge between two cycles indicates a potential fused cycle...
	$CyclesWithCommonEdgesPresent = 1;
      }
      else {
	@{$EdgeIDToCylicPaths{$EdgeID}} = ();
      }
      push @{$EdgeIDToCylicPaths{$EdgeID}}, $ActiveCyclicPath;
    }
  }

  # Associate CyclicPaths with edges...
  for $EdgeID (keys %EdgeIDToCylicPaths) {
    ($VertexID1, $VertexID2) = split($EdgeIDDelimiter, $EdgeID);
    $This->SetEdgeProperty('ActiveCyclicPaths', \@{$EdgeIDToCylicPaths{$EdgeID}}, $VertexID1, $VertexID2);
  }

  if ($CyclesWithCommonEdgesPresent) {
    # Identify fused cycles...
    $This->_IdentifyAndAssociateFusedCyclesWithGraph();
  }

  return $This;
}

# Identify fused cycles and associate them to graph as graph property after cycles
# have been associated with edges...
#
# Note:
#   . During aromaticity detection, fused cycles are treated as one set for counting
#     number of available pi electrons to check against Huckel's rule.
#   . Fused cylce sets contain cycles with at least one common edge between pair
#     of cycles. A specific pair of cycles might not have a direct common edge, but
#     ends up in the same set due to a common edge with another cycle.
#   . Fused cycles are attached to graph as 'FusedActiveCyclicPaths' property with
#     its value as a reference to list of reference where each refernece corresponds
#     to a list of cyclic path objects in a fused set.
#   . For graphs containing fused cycles, non-fused cycles are separeted from fused
#     cycles and attached to the graph as 'NonFusedActiveCyclicPaths'. It's a reference
#     to list containing cylic path objects.
#
sub _IdentifyAndAssociateFusedCyclesWithGraph {
  my($This) = @_;

  # Delete exisiting fused and non-fused cycles...
  $This->_DeleteFusedCyclesAssociatedWithGraph();

  my($ActiveCyclicPathsRef);
  $ActiveCyclicPathsRef = $This->GetGraphProperty('ActiveCyclicPaths');
  if (!@{$ActiveCyclicPathsRef}) {
    # No cycles out there...
    return $This;
  }

  # Get fused cycle pairs...
  my($FusedCyclePairsRef, $FusedCyclesRef, $InValidFusedCycleRef);
  ($FusedCyclePairsRef, $FusedCyclesRef, $InValidFusedCycleRef) = $This->_GetFusedCyclePairs($ActiveCyclicPathsRef);

  # Get fused cycle set indices...
  my($FusedCycleSetsIndicesRef, $FusedCycleSetsCommonEdgesRef);
  $FusedCycleSetsIndicesRef = $This->_GetFusedCycleSetsIndices($FusedCyclePairsRef, $FusedCyclesRef);
  if (!@{$FusedCycleSetsIndicesRef}) {
    # No fused cycles out there...
    return $This;
  }

  # Get fused and non-fused cycles...
  my($FusedCycleSetsRef, $NonFusedCyclesRef);
  ($FusedCycleSetsRef, $NonFusedCyclesRef) = $This->_GetFusedAndNonFusedCycles($ActiveCyclicPathsRef, $FusedCycleSetsIndicesRef, $InValidFusedCycleRef);
  if (!@{$FusedCycleSetsRef}) {
    # No fused cycles out there...
    return $This;
  }

  # Associate fused and non fused cycles to graph....
  $This->SetGraphProperty('FusedActiveCyclicPaths', $FusedCycleSetsRef);
  $This->SetGraphProperty('NonFusedActiveCyclicPaths', $NonFusedCyclesRef);

  return $This;
}

# Collect fused cycle pairs...
#
sub _GetFusedCyclePairs {
  my($This, $ActiveCyclicPathsRef) = @_;

  # Setup a CyclicPathID to CyclicPathIndex map...
  my($CyclicPathIndex, $CyclicPathID, $ActiveCyclicPath, %CyclicPathIDToIndex);

  %CyclicPathIDToIndex = ();
  for $CyclicPathIndex (0 .. $#{$ActiveCyclicPathsRef}) {
    $ActiveCyclicPath = $ActiveCyclicPathsRef->[$CyclicPathIndex];
    $CyclicPathID = "$ActiveCyclicPath";
    $CyclicPathIDToIndex{$CyclicPathID} = $CyclicPathIndex;
  }
  # Go over cycle edges and collect fused cycle pairs...
  my($Index, $VertexID1, $VertexID2, $EdgeCyclicPathsRef, $EdgeID, $CyclicPath1, $CyclicPath2, $CyclicPathID1, $CyclicPathID2, $FusedCyclePairID, $FusedCyclicPath1, $FusedCyclicPath2, $FusedCyclicPathID1, $FusedCyclicPathID2, $FusedCyclicPathIndex1, $FusedCyclicPathIndex2, $FusedCyclePairEdgeCount, @CyclicPathEdgeVertexIDs, %FusedCyclePairs, %CommonEdgeVisited, %CommonEdgesCount, %FusedCycles, %InValidFusedCycles);

  %FusedCyclePairs = (); %CommonEdgeVisited = ();
  %CommonEdgesCount = ();
  %InValidFusedCycles = (); %FusedCycles = ();

  for $ActiveCyclicPath (@{$ActiveCyclicPathsRef}) {
    @CyclicPathEdgeVertexIDs = ();
    @CyclicPathEdgeVertexIDs = $ActiveCyclicPath->GetEdges();
    EDGE: for ($Index = 0; $Index < $#CyclicPathEdgeVertexIDs; $Index += 2) {
      $VertexID1 = $CyclicPathEdgeVertexIDs[$Index]; $VertexID2 = $CyclicPathEdgeVertexIDs[$Index + 1];
      $EdgeCyclicPathsRef = $This->GetEdgeProperty('ActiveCyclicPaths', $VertexID1, $VertexID2);
      if (@{$EdgeCyclicPathsRef} != 2) {
	# Not considered a fused edge...
	next EDGE;
      }
      # Set up a fused cycle pair...
      ($FusedCyclicPath1, $FusedCyclicPath2) = @{$EdgeCyclicPathsRef};
      ($FusedCyclicPathID1, $FusedCyclicPathID2) = ("${FusedCyclicPath1}", "${FusedCyclicPath2}");
      ($FusedCyclicPathIndex1, $FusedCyclicPathIndex2) = ($CyclicPathIDToIndex{$FusedCyclicPathID1}, $CyclicPathIDToIndex{$FusedCyclicPathID2});
      $FusedCyclePairID = ($FusedCyclicPathIndex1 < $FusedCyclicPathIndex2) ? "${FusedCyclicPathIndex1}-${FusedCyclicPathIndex2}" : "${FusedCyclicPathIndex2}-${FusedCyclicPathIndex1}";
      $EdgeID = ($VertexID1 < $VertexID2) ? "${VertexID1}-${VertexID2}" : "${VertexID2}-${VertexID1}";

      if (exists $FusedCyclePairs{$FusedCyclePairID}) {
	if (exists $CommonEdgeVisited{$FusedCyclePairID}{$EdgeID}) {
	  # Edge already processed...
	  next EDGE;
	}
	$CommonEdgeVisited{$FusedCyclePairID}{$EdgeID} = $EdgeID;

	$CommonEdgesCount{$FusedCyclePairID} += 1;
	push @{$FusedCyclePairs{$FusedCyclePairID}}, $EdgeID;
      }
      else {
	@{$FusedCyclePairs{$FusedCyclePairID}} = ();
	push @{$FusedCyclePairs{$FusedCyclePairID}}, $EdgeID;

	%{$CommonEdgeVisited{$FusedCyclePairID}} = ();
	$CommonEdgeVisited{$FusedCyclePairID}{$EdgeID} = $EdgeID;

	$CommonEdgesCount{$FusedCyclePairID} = 1;
      }
    }
  }
  # Valid fused cyle in fused cycle pairs must have only one common egde...
  for $FusedCyclePairID (keys %FusedCyclePairs) {
    ($FusedCyclicPathIndex1, $FusedCyclicPathIndex2) = split /-/, $FusedCyclePairID;
    $FusedCycles{$FusedCyclicPathIndex1} = $FusedCyclicPathIndex1;
    $FusedCycles{$FusedCyclicPathIndex2} = $FusedCyclicPathIndex2;
    if (@{$FusedCyclePairs{$FusedCyclePairID}} != 1) {
      # Mark the cycles involved as invalid fused cycles...
      $InValidFusedCycles{$FusedCyclicPathIndex1} = $FusedCyclicPathIndex1;
      $InValidFusedCycles{$FusedCyclicPathIndex2} = $FusedCyclicPathIndex2;
    }
  }
  return (\%FusedCyclePairs, \%FusedCycles, \%InValidFusedCycles);
}

# Go over fused cycles and set up a graph to collect fused cycle sets. Graph vertices
# correspond to cylce indices; edges correspond to pair of fused cylcles; fused cycle
# sets correspond to connected components. Addionally set up common edges for
# fused cycle sets.
#
sub _GetFusedCycleSetsIndices {
  my($This, $FusedCyclePairsRef, $FusedCyclesRef) = @_;
  my($FusedCyclesGraph, @FusedCycleIndices, @FusedCyclePairIndices, @FusedCycleSetsIndices);

  @FusedCycleIndices = (); @FusedCyclePairIndices = ();
  @FusedCycleSetsIndices = ();

  @FusedCycleIndices = sort { $a <=> $b } keys %{$FusedCyclesRef};
  @FusedCyclePairIndices = map { split /-/, $_; } keys %{$FusedCyclePairsRef};
  if (!@FusedCycleIndices) {
    # No fused cycles out there...
    return \@FusedCycleSetsIndices;
  }
  $FusedCyclesGraph = new Graph(@FusedCycleIndices);
  $FusedCyclesGraph->AddEdges(@FusedCyclePairIndices);

  @FusedCycleSetsIndices = $FusedCyclesGraph->GetConnectedComponentsVertices();

  return \@FusedCycleSetsIndices;
}

# Go over indices of fused cycle sets and map cyclic path indices to cyclic path objects.
# For fused sets containing a cycle with more than one common edge, the whole set is treated
# as non-fused set...
#
sub _GetFusedAndNonFusedCycles {
  my($This, $ActiveCyclicPathsRef, $FusedCycleSetsIndicesRef, $InValidFusedCycleRef) = @_;
  my($CycleSetIndicesRef, $CyclicPathIndex, $ValidFusedCycleSet, @FusedCycleSets, @UnsortedNonFusedCycles, @NonFusedCycles, %CycleIndexVisited);

  @FusedCycleSets = (); @NonFusedCycles = (); @UnsortedNonFusedCycles = ();
  %CycleIndexVisited = ();
  for $CycleSetIndicesRef (@{$FusedCycleSetsIndicesRef}) {
    # Is it a valid fused cycle set? Fused cycle set containing any cycle with more than one common
    # edge is considered invalid and all its cycles are treated as non-fused cycles.
    $ValidFusedCycleSet = 1;
    for $CyclicPathIndex (@{$CycleSetIndicesRef}) {
      $CycleIndexVisited{$CyclicPathIndex} = $CyclicPathIndex;
      if (exists $InValidFusedCycleRef->{$CyclicPathIndex}) {
	$ValidFusedCycleSet = 0;
      }
    }
    if ($ValidFusedCycleSet) {
      my(@FusedCycleSet);
      @FusedCycleSet = ();
      push @FusedCycleSet, sort { $a->GetLength() <=> $b->GetLength() } map { $ActiveCyclicPathsRef->[$_] } @{$CycleSetIndicesRef};
      push @FusedCycleSets, \@FusedCycleSet;
    }
    else {
      push @UnsortedNonFusedCycles, map { $ActiveCyclicPathsRef->[$_] } @{$CycleSetIndicesRef};
    }
  }

  # Add any leftover cycles to non-fused cycles list...
  CYCLICPATH: for $CyclicPathIndex (0 .. $#{$ActiveCyclicPathsRef}) {
    if (exists $CycleIndexVisited{$CyclicPathIndex}) {
      next CYCLICPATH;
    }
    push @UnsortedNonFusedCycles, $ActiveCyclicPathsRef->[$CyclicPathIndex];
  }
  @NonFusedCycles = sort { $a->GetLength() <=> $b->GetLength() } @UnsortedNonFusedCycles;

  return (\@FusedCycleSets, \@NonFusedCycles);
}

# Delete cycles associated with graph...
#
sub _DeleteCyclesAssociatedWithGraph {
  my($This) = @_;

  if ($This->HasGraphProperty('ActiveCyclicPaths')) {
    $This->DeleteGraphProperty('ActiveCyclicPaths');
    $This->DeleteGraphProperty('AllCyclicPaths');
    $This->DeleteGraphProperty('IndependentCyclicPaths');
  }
  return $This;
}

# Delete cycles associated with vertices...
#
sub _DeleteCyclesAssociatedWithVertices {
  my($This) = @_;
  my($VertexID, @VertexIDs);

  @VertexIDs = ();
  @VertexIDs = $This->GetVertices();
  for $VertexID (@VertexIDs) {
    if ($This->HasVertexProperty('ActiveCyclicPaths', $VertexID)) {
      $This->DeleteVertexProperty('ActiveCyclicPaths', $VertexID);
    }
  }
  return $This;
}

# Delete cycles associated with edges...
#
sub _DeleteCyclesAssociatedWithEdges {
  my($This) = @_;
  my($Index, $VertexID1, $VertexID2, @EdgeVertexIDs);

  @EdgeVertexIDs = ();
  @EdgeVertexIDs = $This->GetEdges();
  for ($Index = 0; $Index < $#EdgeVertexIDs; $Index += 2) {
    $VertexID1 = $EdgeVertexIDs[$Index]; $VertexID2 = $EdgeVertexIDs[$Index + 1];
    if ($This->HasEdgeProperty('ActiveCyclicPaths', $VertexID1, $VertexID2)) {
      $This->DeleteEdgeProperty('ActiveCyclicPaths', $VertexID1, $VertexID2);
    }
  }
  return $This;
}

# Delete fused cycles associated with edges...
#
sub _DeleteFusedCyclesAssociatedWithGraph {
  my($This) = @_;

  # Delete exisiting cycles...
  if ($This->HasGraphProperty('FusedActiveCyclicPaths')) {
    $This->DeleteGraphProperty('FusedActiveCyclicPaths');
    $This->DeleteGraphProperty('NonFusedActiveCyclicPaths');
  }
  return $This;
}

# Does graph contains any cycles?
#
sub IsAcyclic {
  my($This) = @_;

  return $This->GetNumOfCycles() ? 0 : 1;
}

# Does graph contains cycles?
#
sub IsCyclic {
  my($This) = @_;

  return $This->GetNumOfCycles() ? 1 : 0;
}

# Does graph contains only any cycle?
#
sub IsUnicyclic {
  my($This) = @_;

  return ($This->GetNumOfCycles() == 1) ? 1 : 0;
}

# Get size of smallest cycle in graph...
#
sub GetGirth {
  my($This) = @_;

  return $This->GetSizeOfSmallestCycle();
}

# Get size of smallest cycle in graph...
#
sub GetSizeOfSmallestCycle {
  my($This) = @_;

  return $This->_GetCycleSize('GraphCycle', 'SmallestCycle');
}

# Get size of largest cycle in graph...
#
sub GetCircumference {
  my($This) = @_;

  return $This->GetSizeOfLargestCycle();
}

# Get size of largest cycle in graph...
#
sub GetSizeOfLargestCycle {
  my($This) = @_;

  return $This->_GetCycleSize('GraphCycle', 'LargestCycle');
}

# Get number of cycles in graph...
#
sub GetNumOfCycles {
  my($This) = @_;

  return $This->_GetNumOfCyclesWithSize('GraphCycle', 'AllSizes');
}

# Get number of cycles with odd size in graph...
#
sub GetNumOfCyclesWithOddSize {
  my($This) = @_;

  return $This->_GetNumOfCyclesWithSize('GraphCycle', 'OddSize');
}

# Get number of cycles with even size in graph...
#
sub GetNumOfCyclesWithEvenSize {
  my($This) = @_;

  return $This->_GetNumOfCyclesWithSize('GraphCycle', 'EvenSize');
}

# Get number of cycles with specific size in graph...
#
sub GetNumOfCyclesWithSize {
  my($This, $CycleSize) = @_;

  return $This->_GetNumOfCyclesWithSize('GraphCycle', 'SpecifiedSize', $CycleSize);
}

# Get number of cycles with size less than a specific size in graph...
#
sub GetNumOfCyclesWithSizeLessThan {
  my($This, $CycleSize) = @_;

  return $This->_GetNumOfCyclesWithSize('GraphCycle', 'SizeLessThan', $CycleSize);
}

# Get number of cycles with size greater than a specific size in graph...
#
sub GetNumOfCyclesWithSizeGreaterThan {
  my($This, $CycleSize) = @_;

  return $This->_GetNumOfCyclesWithSize('GraphCycle', 'SizeGreaterThan', $CycleSize);
}

# Get largest cyclic path in graph...
#
sub GetLargestCycle {
  my($This) = @_;

  return $This->_GetCycle('GraphCycle', 'LargestCycle');
}

# Get smallest cyclic path in graph...
#
sub GetSmallestCycle {
  my($This) = @_;

  return $This->_GetCycle('GraphCycle', 'SmallestCycle');
}

# Get all cycles in graph...
#
sub GetCycles {
  my($This) = @_;

  return $This->_GetCyclesWithSize('GraphCycle', 'AllSizes');
}

# Get cycles with odd size in graph...
#
sub GetCyclesWithOddSize {
  my($This) = @_;

  return $This->_GetCyclesWithSize('GraphCycle', 'OddSize');
}

# Get cycles with even size in graph...
#
sub GetCyclesWithEvenSize {
  my($This) = @_;

  return $This->_GetCyclesWithSize('GraphCycle', 'EvenSize');
}

# Get cycles with specific size in graph...
#
sub GetCyclesWithSize {
  my($This, $CycleSize) = @_;

  return $This->_GetCyclesWithSize('GraphCycle', 'SpecifiedSize', $CycleSize);
}

# Get cycles with size less than a specific size in graph...
#
sub GetCyclesWithSizeLessThan {
  my($This, $CycleSize) = @_;

  return $This->_GetCyclesWithSize('GraphCycle', 'SizeLessThan', $CycleSize);
}

# Get cycles with size greater than a specific size in graph...
#
sub GetCyclesWithSizeGreaterThan {
  my($This, $CycleSize) = @_;

  return $This->_GetCyclesWithSize('GraphCycle', 'SizeGreaterThan', $CycleSize);
}

# Is vertex in a cycle?
#
sub IsCyclicVertex {
  my($This, $VertexID) = @_;

  return $This->GetNumOfVertexCycles($VertexID) ? 1 : 0;
}

# Is vertex in a only one cycle?
#
sub IsUnicyclicVertex {
  my($This, $VertexID) = @_;

  return ($This->GetNumOfVertexCycles($VertexID)  == 1) ? 1 : 0;
}

# Is vertex not in a cycle?
#
sub IsAcyclicVertex {
  my($This, $VertexID) = @_;

  return $This->GetNumOfVertexCycles($VertexID) ? 0 : 1;
}

# Get size of smallest cycle containing specified vertex...
#
sub GetSizeOfSmallestVertexCycle {
  my($This, $VertexID) = @_;

  return $This->_GetCycleSize('VertexCycle', 'SmallestCycle', $VertexID);
}

# Get size of largest cycle containing specified vertex...
#
sub GetSizeOfLargestVertexCycle {
  my($This, $VertexID) = @_;

  return $This->_GetCycleSize('VertexCycle', 'LargestCycle', $VertexID);
}

# Get number of cycles containing specified vertex...
#
sub GetNumOfVertexCycles {
  my($This, $VertexID) = @_;

  return $This->_GetNumOfCyclesWithSize('VertexCycle', 'AllSizes', 0, $VertexID);
}

# Get number of cycles with odd size containing specified vertex...
#
sub GetNumOfVertexCyclesWithOddSize {
  my($This, $VertexID) = @_;

  return $This->_GetNumOfCyclesWithSize('VertexCycle', 'OddSize', 0, $VertexID);
}

# Get number of cycles with even size containing specified vertex...
#
sub GetNumOfVertexCyclesWithEvenSize {
  my($This, $VertexID) = @_;

  return $This->_GetNumOfCyclesWithSize('VertexCycle', 'EvenSize', 0, $VertexID);
}

# Get number of cycles with specified size containing specified vertex...
#
sub GetNumOfVertexCyclesWithSize {
  my($This, $VertexID, $CycleSize) = @_;

  return $This->_GetNumOfCyclesWithSize('VertexCycle', 'SpecifiedSize', $CycleSize, $VertexID);
}

# Get number of cycles with size less than specified size containing specified vertex...
#
sub GetNumOfVertexCyclesWithSizeLessThan {
  my($This, $VertexID, $CycleSize) = @_;

  return $This->_GetNumOfCyclesWithSize('VertexCycle', 'SizeLessThan', $CycleSize, $VertexID);
}

# Get number of cycles with size greater than specified size containing specified vertex...
#
sub GetNumOfVertexCyclesWithSizeGreaterThan {
  my($This, $VertexID, $CycleSize) = @_;

  return $This->_GetNumOfCyclesWithSize('VertexCycle', 'SizeGreaterThan', $CycleSize, $VertexID);
}

# Get smallest cycle containing specified vertex...
#
sub GetSmallestVertexCycle {
  my($This, $VertexID) = @_;

  return $This->_GetCycle('VertexCycle', 'SmallestCycle', $VertexID);
}

# Get largest cycle containing specified vertex...
#
sub GetLargestVertexCycle {
  my($This, $VertexID) = @_;

  return $This->_GetCycle('VertexCycle', 'LargestCycle', $VertexID);
}

# Get cycles containing specified vertex...
#
sub GetVertexCycles {
  my($This, $VertexID) = @_;

  return $This->_GetCyclesWithSize('VertexCycle', 'AllSizes', 0, $VertexID);
}

# Get cycles with odd size containing specified vertex...
#
sub GetVertexCyclesWithOddSize {
  my($This, $VertexID) = @_;

  return $This->_GetCyclesWithSize('VertexCycle', 'OddSize', 0, $VertexID);
}

# Get cycles with even size containing specified vertex...
#
sub GetVertexCyclesWithEvenSize {
  my($This, $VertexID) = @_;

  return $This->_GetCyclesWithSize('VertexCycle', 'EvenSize', 0, $VertexID);
}

# Get cycles with specified size containing specified vertex...
#
sub GetVertexCyclesWithSize {
  my($This, $VertexID, $CycleSize) = @_;

  return $This->_GetCyclesWithSize('VertexCycle', 'SpecifiedSize', $CycleSize, $VertexID);
}

# Get cycles with size less than specified size containing specified vertex...
#
sub GetVertexCyclesWithSizeLessThan {
  my($This, $VertexID, $CycleSize) = @_;

  return $This->_GetCyclesWithSize('VertexCycle', 'SizeLessThan', $CycleSize, $VertexID);
}

# Get cycles with size greater than specified size containing specified vertex...
#
sub GetVertexCyclesWithSizeGreaterThan {
  my($This, $VertexID, $CycleSize) = @_;

  return $This->_GetCyclesWithSize('VertexCycle', 'SizeGreaterThan', $CycleSize, $VertexID);
}

# Is edge in a cycle?
#
sub IsCyclicEdge {
  my($This, $VertexID1, $VertexID2) = @_;

  return $This->GetNumOfEdgeCycles($VertexID1, $VertexID2) ? 1 : 0;
}

# Is edge in a only one cycle?
#
sub IsUnicyclicEdge {
  my($This, $VertexID1, $VertexID2) = @_;

  return ($This->GetNumOfEdgeCycles($VertexID1, $VertexID2)  == 1) ? 1 : 0;
}

# Is Edge not in a cycle?
#
sub IsAcyclicEdge {
  my($This, $VertexID1, $VertexID2) = @_;

  return $This->GetNumOfEdgeCycles($VertexID1, $VertexID2) ? 0 : 1;
}

# Get size of smallest cycle containing specified edge...
#
sub GetSizeOfSmallestEdgeCycle {
  my($This, $VertexID1, $VertexID2) = @_;

  return $This->_GetCycleSize('EdgeCycle', 'SmallestCycle', $VertexID1, $VertexID2);
}

# Get size of largest cycle containing specified edge...
#
sub GetSizeOfLargestEdgeCycle {
  my($This, $VertexID1, $VertexID2) = @_;

  return $This->_GetCycleSize('EdgeCycle', 'LargestCycle', $VertexID1, $VertexID2);
}

# Get number of cycles containing specified edge...
#
sub GetNumOfEdgeCycles {
  my($This, $VertexID1, $VertexID2) = @_;

  return $This->_GetNumOfCyclesWithSize('EdgeCycle', 'AllSizes', 0, $VertexID1, $VertexID2);
}

# Get number of cycles with odd size containing specified edge...
#
sub GetNumOfEdgeCyclesWithOddSize {
  my($This, $VertexID1, $VertexID2) = @_;

  return $This->_GetNumOfCyclesWithSize('EdgeCycle', 'OddSize', 0, $VertexID1, $VertexID2);
}

# Get number of cycles with even size containing specified edge...
#
sub GetNumOfEdgeCyclesWithEvenSize {
  my($This, $VertexID1, $VertexID2) = @_;

  return $This->_GetNumOfCyclesWithSize('EdgeCycle', 'EvenSize', 0, $VertexID1, $VertexID2);
}

# Get number of cycles with specified size containing specified edge...
#
sub GetNumOfEdgeCyclesWithSize {
  my($This, $VertexID1, $VertexID2, $CycleSize) = @_;

  return $This->_GetNumOfCyclesWithSize('EdgeCycle', 'SpecifiedSize', $CycleSize, $VertexID1, $VertexID2);
}

# Get number of cycles with size less than specified size containing specified edge...
#
sub GetNumOfEdgeCyclesWithSizeLessThan {
  my($This, $VertexID1, $VertexID2, $CycleSize) = @_;

  return $This->_GetNumOfCyclesWithSize('EdgeCycle', 'SizeLessThan', $CycleSize, $VertexID1, $VertexID2);
}

# Get number of cycles with size greater than specified size containing specified edge...
#
sub GetNumOfEdgeCyclesWithSizeGreaterThan {
  my($This, $VertexID1, $VertexID2, $CycleSize) = @_;

  return $This->_GetNumOfCyclesWithSize('EdgeCycle', 'SizeGreaterThan', $CycleSize, $VertexID1, $VertexID2);
}

# Get smallest cycle containing specified edge...
#
sub GetSmallestEdgeCycle {
  my($This, $VertexID1, $VertexID2) = @_;

  return $This->_GetCycle('EdgeCycle', 'SmallestCycle', $VertexID1, $VertexID2);
}

# Get largest cycle containing specified edge...
#
sub GetLargestEdgeCycle {
  my($This, $VertexID1, $VertexID2) = @_;

  return $This->_GetCycle('EdgeCycle', 'LargestCycle', $VertexID1, $VertexID2);
}

# Get cycles containing specified edge...
#
sub GetEdgeCycles {
  my($This, $VertexID1, $VertexID2) = @_;

  return $This->_GetCyclesWithSize('EdgeCycle', 'AllSizes', 0, $VertexID1, $VertexID2);
}

# Get cycles with odd size containing specified edge...
#
sub GetEdgeCyclesWithOddSize {
  my($This, $VertexID1, $VertexID2) = @_;

  return $This->_GetCyclesWithSize('EdgeCycle', 'OddSize', 0, $VertexID1, $VertexID2);
}

# Get cycles with even size containing specified edge...
#
sub GetEdgeCyclesWithEvenSize {
  my($This, $VertexID1, $VertexID2) = @_;

  return $This->_GetCyclesWithSize('EdgeCycle', 'EvenSize', 0, $VertexID1, $VertexID2);
}

# Get cycles with specified size containing specified edge...
#
sub GetEdgeCyclesWithSize {
  my($This, $VertexID1, $VertexID2, $CycleSize) = @_;

  return $This->_GetCyclesWithSize('EdgeCycle', 'SpecifiedSize', $CycleSize, $VertexID1, $VertexID2);
}

# Get cycles with size less than specified size containing specified edge...
#
sub GetEdgeCyclesWithSizeLessThan {
  my($This, $VertexID1, $VertexID2, $CycleSize) = @_;

  return $This->_GetCyclesWithSize('EdgeCycle', 'SizeLessThan', $CycleSize, $VertexID1, $VertexID2);
}

# Get cycles with size greater than specified size containing specified edge...
#
sub GetEdgeCyclesWithSizeGreaterThan {
  my($This, $VertexID1, $VertexID2, $CycleSize) = @_;

  return $This->_GetCyclesWithSize('EdgeCycle', 'SizeGreaterThan', $CycleSize, $VertexID1, $VertexID2);
}

# Get size of specified cycle type...
#
sub _GetCycleSize {
  my($This, $Mode, $CycleSize, $VertexID1, $VertexID2) = @_;
  my($ActiveCyclicPathsRef, $CyclicPath, $Size);

  if (!$This->HasGraphProperty('ActiveCyclicPaths')) {
    return 0;
  }
  if ($Mode =~ /^VertexCycle$/i) {
    if (!$This->HasVertexProperty('ActiveCyclicPaths', $VertexID1)) {
      return 0;
    }
  }
  elsif ($Mode =~ /^EdgeCycle$/i) {
    if (!$This->HasEdgeProperty('ActiveCyclicPaths', $VertexID1, $VertexID2)) {
      return 0;
    }
  }

  MODE: {
      if ($Mode =~ /^GraphCycle$/i) { $ActiveCyclicPathsRef = $This->GetGraphProperty('ActiveCyclicPaths'); last MODE; }
      if ($Mode =~ /^VertexCycle$/i) { $ActiveCyclicPathsRef = $This->GetVertexProperty('ActiveCyclicPaths', $VertexID1); last MODE; }
      if ($Mode =~ /^EdgeCycle$/i) { $ActiveCyclicPathsRef = $This->GetEdgeProperty('ActiveCyclicPaths', $VertexID1, $VertexID2); last MODE; }
      return 0;
  }

  if (!@{$ActiveCyclicPathsRef}) {
    return 0;
  }

  CYCLESIZE: {
      if ($CycleSize =~ /^SmallestCycle$/i) { $CyclicPath = $ActiveCyclicPathsRef->[0]; last CYCLESIZE; }
      if ($CycleSize =~ /^LargestCycle$/i) { $CyclicPath = $ActiveCyclicPathsRef->[$#{$ActiveCyclicPathsRef}]; last CYCLESIZE; }
      return 0;
  }
  $Size = $CyclicPath->GetLength() - 1;

  return $Size;
}

# Get of specified cycle size...
#
sub _GetCycle {
  my($This, $Mode, $CycleSize, $VertexID1, $VertexID2) = @_;
  my($ActiveCyclicPathsRef, $CyclicPath, $Size);

  if (!$This->HasGraphProperty('ActiveCyclicPaths')) {
    return $This->_GetEmptyCycles();
  }
  if ($Mode =~ /^VertexCycle$/i) {
    if (!$This->HasVertexProperty('ActiveCyclicPaths', $VertexID1)) {
      return $This->_GetEmptyCycles();
    }
  }
  elsif ($Mode =~ /^EdgeCycle$/i) {
    if (!$This->HasEdgeProperty('ActiveCyclicPaths', $VertexID1, $VertexID2)) {
      return $This->_GetEmptyCycles();
    }
  }

  MODE: {
      if ($Mode =~ /^GraphCycle$/i) { $ActiveCyclicPathsRef = $This->GetGraphProperty('ActiveCyclicPaths'); last MODE; }
      if ($Mode =~ /^VertexCycle$/i) { $ActiveCyclicPathsRef = $This->GetVertexProperty('ActiveCyclicPaths', $VertexID1); last MODE; }
      if ($Mode =~ /^EdgeCycle$/i) { $ActiveCyclicPathsRef = $This->GetEdgeProperty('ActiveCyclicPaths', $VertexID1, $VertexID2); last MODE; }
      return $This->_GetEmptyCycles();
  }

  if (!@{$ActiveCyclicPathsRef}) {
    return $This->_GetEmptyCycles();
  }

  CYCLESIZE: {
      if ($CycleSize =~ /^SmallestCycle$/i) { $CyclicPath = $ActiveCyclicPathsRef->[0]; last CYCLESIZE; }
      if ($CycleSize =~ /^LargestCycle$/i) { $CyclicPath = $ActiveCyclicPathsRef->[$#{$ActiveCyclicPathsRef}]; last CYCLESIZE; }
      return $This->_GetEmptyCycles();
  }
  return $CyclicPath;
}

# Get num of cycles in graph...
#
sub _GetNumOfCyclesWithSize {
  my($This, $Mode, $SizeMode, $SpecifiedSize, $VertexID1, $VertexID2) = @_;
  my($ActiveCyclicPathsRef);

  if (!$This->HasGraphProperty('ActiveCyclicPaths')) {
    return 0;
  }
  if ($Mode =~ /^VertexCycle$/i) {
    if (!$This->HasVertexProperty('ActiveCyclicPaths', $VertexID1)) {
      return 0;
    }
  }
  elsif ($Mode =~ /^EdgeCycle$/i) {
    if (!$This->HasEdgeProperty('ActiveCyclicPaths', $VertexID1, $VertexID2)) {
      return 0;
    }
  }

  if ($SizeMode =~ /^(SizeLessThan|SizeGreaterThan|SpecifiedSize)$/i) {
    if (!defined $SpecifiedSize) {
      carp "Warning: ${ClassName}->_GetNumOfCyclesWithSize: Cycle size muse be defined...";
      return 0;
    }
    if ($SpecifiedSize < 0) {
      carp "Warning: ${ClassName}->_GetNumOfCyclesWithSize: Specified cycle size, $SpecifiedSize, must be > 0 ...";
      return 0;
    }
  }

  MODE: {
      if ($Mode =~ /^GraphCycle$/i) { $ActiveCyclicPathsRef = $This->GetGraphProperty('ActiveCyclicPaths'); last MODE; }
      if ($Mode =~ /^VertexCycle$/i) { $ActiveCyclicPathsRef = $This->GetVertexProperty('ActiveCyclicPaths', $VertexID1); last MODE; }
      if ($Mode =~ /^EdgeCycle$/i) { $ActiveCyclicPathsRef = $This->GetEdgeProperty('ActiveCyclicPaths', $VertexID1, $VertexID2); last MODE; }
      return 0;
  }

  if (!@{$ActiveCyclicPathsRef}) {
    return 0;
  }
  my($NumOfCycles);

  $NumOfCycles = $This->_GetCycles($Mode, $ActiveCyclicPathsRef, $SizeMode, $SpecifiedSize);

  return $NumOfCycles;
}

# Get cycles in graph...
#
sub _GetCyclesWithSize {
  my($This, $Mode, $SizeMode, $SpecifiedSize, $VertexID1, $VertexID2) = @_;
  my($ActiveCyclicPathsRef);

  if (!$This->HasGraphProperty('ActiveCyclicPaths')) {
    return $This->_GetEmptyCycles();
  }
  if ($Mode =~ /^VertexCycle$/i) {
    if (!$This->HasVertexProperty('ActiveCyclicPaths', $VertexID1)) {
      return $This->_GetEmptyCycles();
    }
  }
  elsif ($Mode =~ /^EdgeCycle$/i) {
    if (!$This->HasEdgeProperty('ActiveCyclicPaths', $VertexID1, $VertexID2)) {
      return $This->_GetEmptyCycles();
    }
  }

  if ($SizeMode =~ /^(SizeLessThan|SizeGreaterThan|SpecifiedSize)$/i) {
    if (!defined $SpecifiedSize) {
      carp "Warning: ${ClassName}->_GetCyclesWithSize: Cycle size must be defined...";
      return $This->_GetEmptyCycles();
    }
    if ($SpecifiedSize < 0) {
      carp "Warning: ${ClassName}->_GetCyclesWithSize: Specified cycle size, $SpecifiedSize, must be > 0 ...";
      return $This->_GetEmptyCycles();
    }
  }

  MODE: {
      if ($Mode =~ /^GraphCycle$/i) { $ActiveCyclicPathsRef = $This->GetGraphProperty('ActiveCyclicPaths'); last MODE; }
      if ($Mode =~ /^VertexCycle$/i) { $ActiveCyclicPathsRef = $This->GetVertexProperty('ActiveCyclicPaths', $VertexID1); last MODE; }
      if ($Mode =~ /^EdgeCycle$/i) { $ActiveCyclicPathsRef = $This->GetEdgeProperty('ActiveCyclicPaths', $VertexID1, $VertexID2); last MODE; }
      return $This->_GetEmptyCycles();
    }

  if (!@{$ActiveCyclicPathsRef}) {
    return $This->_GetEmptyCycles();
  }
  return $This->_GetCycles($Mode, $ActiveCyclicPathsRef, $SizeMode, $SpecifiedSize);
}

# Get cycles information...
#
sub _GetCycles {
  my($This, $Mode, $ActiveCyclicPathsRef, $SizeMode, $SpecifiedSize) = @_;

  if (!@{$ActiveCyclicPathsRef}) {
    return $This->_GetEmptyCycles();
  }

  if ($SizeMode =~ /^AllSizes$/i) {
    return wantarray ? @{$ActiveCyclicPathsRef} : scalar @{$ActiveCyclicPathsRef};
  }

  # Get appropriate cycles...
  my($Size, $CyclicPath, @FilteredCyclicPaths);
  @FilteredCyclicPaths = ();

  for $CyclicPath (@{$ActiveCyclicPathsRef}) {
    $Size = $CyclicPath->GetLength() - 1;
    SIZEMODE: {
      if ($SizeMode =~ /^OddSize$/i) { if ($Size % 2) { push @FilteredCyclicPaths, $CyclicPath; } last SIZEMODE; }
      if ($SizeMode =~ /^EvenSize$/i) { if (!($Size % 2)) { push @FilteredCyclicPaths, $CyclicPath; } last SIZEMODE; }
      if ($SizeMode =~ /^SizeLessThan$/i) { if ($Size < $SpecifiedSize) { push @FilteredCyclicPaths, $CyclicPath; } last SIZEMODE; }
      if ($SizeMode =~ /^SizeGreaterThan$/i) { if ($Size > $SpecifiedSize) { push @FilteredCyclicPaths, $CyclicPath; } last SIZEMODE; }
      if ($SizeMode =~ /^SpecifiedSize$/i) { if ($Size == $SpecifiedSize) { push @FilteredCyclicPaths, $CyclicPath; } last SIZEMODE; }
      return undef;
    }
  }
  return wantarray ? @FilteredCyclicPaths : scalar @FilteredCyclicPaths;
}

# Return empty cyles array...
#
sub _GetEmptyCycles {
  my($This) = @_;
  my(@CyclicPaths);

  @CyclicPaths = ();

  return wantarray ? @CyclicPaths : scalar @CyclicPaths;
}

# Does graph contains fused cycles?
sub HasFusedCycles {
  my($This) = @_;

  return ($This->HasGraphProperty('FusedActiveCyclicPaths')) ? 1 : 0;
}

# Return a reference to fused cycle sets lists containing references to lists of cyclic path objects
# in each fused cycle set and a reference to a list containing non-fused cyclic paths...
#
sub GetFusedAndNonFusedCycles {
  my($This) = @_;
  my($FusedCycleSetsRef, $NonFusedCyclesRef);

  $FusedCycleSetsRef = $This->HasGraphProperty('FusedActiveCyclicPaths') ? $This->GetGraphProperty('FusedActiveCyclicPaths') : undef;
  $NonFusedCyclesRef = $This->HasGraphProperty('NonFusedActiveCyclicPaths') ? $This->GetGraphProperty('NonFusedActiveCyclicPaths') : undef;

  return ($FusedCycleSetsRef, $NonFusedCyclesRef);
}

# Get vertices of connected components as a list containing references to
# lists of vertices for each component  sorted in order of its decreasing size...
#
sub GetConnectedComponentsVertices {
  my($This) = @_;
  my($PathsTraversal);

  $PathsTraversal = new Graph::PathsTraversal($This);
  $PathsTraversal->PerformDepthFirstSearch();

  return $PathsTraversal->GetConnectedComponentsVertices();
}

# Get a list of topologically sorted vertrices starting from a specified vertex or
# an arbitrary vertex in the graph...
#
sub GetTopologicallySortedVertices {
  my($This, $RootVertexID) = @_;
  my($PathsTraversal);

  $PathsTraversal = new Graph::PathsTraversal($This);
  $PathsTraversal->PerformBreadthFirstSearch($RootVertexID);

  return $PathsTraversal->GetVertices();
}

# Get a list of paths starting from a specified vertex with length upto specified length
# and no sharing of edges in paths traversed. By default, cycles are included in paths.
# A path containing a cycle is terminated at a vertex completing the cycle.
#
sub GetPathsStartingAtWithLengthUpto {
  my($This, $StartVertexID, $Length, $AllowCycles) = @_;
  my($PathsTraversal);

  $PathsTraversal = new Graph::PathsTraversal($This);
  $PathsTraversal->PerformPathsSearchWithLengthUpto($StartVertexID, $Length, $AllowCycles);

  return $PathsTraversal->GetPaths();
}

# Get a list of paths starting from a specified vertex with specified length
# and no sharing of edges in paths traversed. By default, cycles are included in paths.
# A path containing a cycle is terminated at a vertex completing the cycle.
#
sub GetPathsStartingAtWithLength {
  my($This, $StartVertexID, $Length, $AllowCycles) = @_;
  my($PathsTraversal);

  $PathsTraversal = new Graph::PathsTraversal($This);
  $PathsTraversal->PerformPathsSearchWithLength($StartVertexID, $Length, $AllowCycles);

  return $PathsTraversal->GetPaths();
}

# Get a list of paths with all possible lengths starting from a specified vertex
# with no sharing of edges in paths traversed. By default, cycles are included in paths.
# A path containing a cycle is terminated at a vertex completing the cycle.
#
sub GetPathsStartingAt {
  my($This, $StartVertexID, $AllowCycles) = @_;
  my($PathsTraversal);

  $PathsTraversal = new Graph::PathsTraversal($This);
  $PathsTraversal->PerformPathsSearch($StartVertexID, $AllowCycles);

  return $PathsTraversal->GetPaths();
}

# Get a list of all paths starting from a specified vertex with length upto a specified length
# with sharing of edges in paths traversed. By default, cycles are included in paths.
# A path containing a cycle is terminated at a vertex completing the cycle.
#
sub GetAllPathsStartingAtWithLengthUpto {
  my($This, $StartVertexID, $Length, $AllowCycles) = @_;
  my($PathsTraversal);

  $PathsTraversal = new Graph::PathsTraversal($This);
  $PathsTraversal->PerformAllPathsSearchWithLengthUpto($StartVertexID, $Length, $AllowCycles);

  return $PathsTraversal->GetPaths();
}

# Get a list of all paths starting from a specified vertex with specified length
# with sharing of edges in paths traversed. By default, cycles are included in paths.
# A path containing a cycle is terminated at a vertex completing the cycle.
#
sub GetAllPathsStartingAtWithLength {
  my($This, $StartVertexID, $Length, $AllowCycles) = @_;
  my($PathsTraversal);

  $PathsTraversal = new Graph::PathsTraversal($This);
  $PathsTraversal->PerformAllPathsSearchWithLength($StartVertexID, $Length, $AllowCycles);

  return $PathsTraversal->GetPaths();
}


# Get a list of all paths with all possible lengths starting from a specified vertex
# with sharing of edges in paths traversed. By default, cycles are included in paths.
# A path containing a cycle is terminated at a vertex completing the cycle.
#
sub GetAllPathsStartingAt {
  my($This, $StartVertexID, $AllowCycles) = @_;
  my($PathsTraversal);

  $PathsTraversal = new Graph::PathsTraversal($This);
  $PathsTraversal->PerformAllPathsSearch($StartVertexID, $AllowCycles);

  return $PathsTraversal->GetPaths();
}

# Get a reference to list of paths starting from each vertex in graph with length upto specified
# length and no sharing of edges in paths traversed. By default, cycles are included in paths.
# A path containing a cycle is terminated at a vertex completing the cycle.
#
sub GetPathsWithLengthUpto {
  my($This, $Length, $AllowCycles) = @_;

  $AllowCycles = (defined $AllowCycles) ? $AllowCycles : 1;

  return $This->_GetPaths('PathsWithLengthUpto', $Length, $AllowCycles);
}

# Get a reference to list of paths starting from each vertex in graph with specified
# length and no sharing of edges in paths traversed. By default, cycles are included in paths.
# A path containing a cycle is terminated at a vertex completing the cycle.
#
sub GetPathsWithLength {
  my($This, $Length, $AllowCycles) = @_;

  $AllowCycles = (defined $AllowCycles) ? $AllowCycles : 1;

  return $This->_GetPaths('PathsWithLength', $Length, $AllowCycles);
}

# Get a reference to list of paths with all possible lengths starting from each vertex
# with no sharing of edges in paths traversed. By default, cycles are included in paths.
# A path containing a cycle is terminated at a vertex completing the cycle.
#
#
sub GetPaths {
  my($This, $AllowCycles) = @_;

  $AllowCycles = (defined $AllowCycles) ? $AllowCycles : 1;

  return $This->_GetPaths('PathsWithAllLengths', undef, $AllowCycles);
}

# Get a reference to list of all paths starting from each vertex in graph with length upto a specified
# length with sharing of edges in paths traversed. By default, cycles are included in paths. A path
# containing a cycle is terminated at a vertex completing the cycle.
#
# Note:
#   . Duplicate paths are not removed.
#
sub GetAllPathsWithLengthUpto {
  my($This, $Length, $AllowCycles) = @_;

  $AllowCycles = (defined $AllowCycles) ? $AllowCycles : 1;

  return $This->_GetPaths('AllPathsWithLengthUpto', $Length, $AllowCycles);
}

# Get a reference to list of all paths starting from each vertex in graph with specified
# length with sharing of edges in paths traversed. By default, cycles are included in paths. A path
# containing a cycle is terminated at a vertex completing the cycle.
#
# Note:
#   . Duplicate paths are not removed.
#
sub GetAllPathsWithLength {
  my($This, $Length, $AllowCycles) = @_;

  $AllowCycles = (defined $AllowCycles) ? $AllowCycles : 1;

  return $This->_GetPaths('AllPathsWithLength', $Length, $AllowCycles);
}

# Get a reference to list of all paths with all possible lengths starting from each vertex in graph
# with sharing of edges in paths traversed. By default, cycles are included in paths. A path
# containing a cycle is terminated at a vertex completing the cycle.
#
# Note:
#   . Duplicate paths are not removed.
#
sub GetAllPaths {
  my($This, $AllowCycles) = @_;

  $AllowCycles = (defined $AllowCycles) ? $AllowCycles : 1;

  return $This->_GetPaths('AllPathsWithAllLengths', undef, $AllowCycles);
}


# Retrieve appropriate paths for each vertex in graph and return a referernce to list
# containing path objects...
#
sub _GetPaths {
  my($This, $Mode, $Length, $AllowCycles) = @_;
  my($VertexID, @EmptyPaths, @Paths);

  @Paths = (); @EmptyPaths = ();

  for $VertexID ($This->GetVertices()) {
    my($Status, $PathsTraversal);

    $PathsTraversal = new Graph::PathsTraversal($This);
    MODE: {
      if ($Mode =~ /^PathsWithLengthUpto$/i) { $Status = $PathsTraversal->PerformPathsSearchWithLengthUpto($VertexID, $Length, $AllowCycles); last MODE; }
      if ($Mode =~ /^PathsWithLength$/i) { $Status = $PathsTraversal->PerformPathsSearchWithLength($VertexID, $Length, $AllowCycles); last MODE; }
      if ($Mode =~ /^PathsWithAllLengths$/i) { $Status = $PathsTraversal->PerformPathsSearch($VertexID, $AllowCycles); last MODE; }

      if ($Mode =~ /^AllPathsWithLengthUpto$/i) { $Status = $PathsTraversal->PerformAllPathsSearchWithLengthUpto($VertexID, $Length, $AllowCycles); last MODE; }
      if ($Mode =~ /^AllPathsWithLength$/i) { $Status = $PathsTraversal->PerformAllPathsSearchWithLength($VertexID, $Length, $AllowCycles); last MODE; }
      if ($Mode =~ /^AllPathsWithAllLengths$/i) { $Status = $PathsTraversal->PerformAllPathsSearch($VertexID, $AllowCycles); last MODE; }

      return \@EmptyPaths;
    }
    if (!defined $Status) {
      return \@EmptyPaths;
    }
    push @Paths, $PathsTraversal->GetPaths();
  }
  return \@Paths;
}

# Get a list of paths between two vertices. For cyclic graphs, the list contains
# may contain two paths...
#
sub GetPathsBetween {
  my($This, $StartVertexID, $EndVertexID) = @_;
  my($Path, $ReversePath, @Paths);

  @Paths = ();

  $Path = $This->_GetPathBetween($StartVertexID, $EndVertexID);
  if (!defined $Path) {
    return \@Paths;
  }

  $ReversePath = $This->_GetPathBetween($EndVertexID, $StartVertexID);
  if (!defined $ReversePath) {
    return \@Paths;
  }
  if ($Path eq $ReversePath) {
    push @Paths, $Path;
  }
  else {
    # Make sure first vertex in reverse path corresponds to specified start vertex ID...
    $ReversePath->Reverse();
    push @Paths, ($Path->GetLength <= $ReversePath->GetLength()) ? ($Path, $ReversePath) : ($ReversePath, $Path);
  }
  return @Paths;
}

# Get a path beween two vertices...
#
sub _GetPathBetween {
  my($This, $StartVertexID, $EndVertexID) = @_;
  my($PathsTraversal,  @Paths);

  $PathsTraversal = new Graph::PathsTraversal($This);
  $PathsTraversal->PerformPathsSearchBetween($StartVertexID, $EndVertexID);

  @Paths = $PathsTraversal->GetPaths();

  return (@Paths) ? $Paths[0] : undef;
}

# Get a list containing lists of neighborhood vertices around a specified vertex with in a
# specified radius...
#
sub GetNeighborhoodVerticesWithRadiusUpto {
  my($This, $StartVertexID, $Radius) = @_;
  my($PathsTraversal);

  $PathsTraversal = new Graph::PathsTraversal($This);
  $PathsTraversal->PerformNeighborhoodVerticesSearchWithRadiusUpto($StartVertexID, $Radius);

  return $PathsTraversal->GetVerticesNeighborhoods();
}

# Get a list containing lists of neighborhood vertices around a specified vertex at all
# radii levels...
#
sub GetNeighborhoodVertices {
  my($This, $StartVertexID) = @_;
  my($PathsTraversal);

  $PathsTraversal = new Graph::PathsTraversal($This);
  $PathsTraversal->PerformNeighborhoodVerticesSearch($StartVertexID);

  return $PathsTraversal->GetVerticesNeighborhoods();
}

# Get neighborhood vertices around a specified vertex, along with their successor connected vertices, collected
# with in a specified radius as a list containing references to lists with first value corresponding to vertex
# ID and second value as reference to a list containing its successor connected vertices.
#
# For a neighborhood vertex at each radius level, the successor connected vertices correspond to the
# neighborhood vertices at the next radius level. Consequently, the neighborhood vertices at the last
# radius level don't contain any successor vertices which fall outside the range of specified radius.
#
sub GetNeighborhoodVerticesWithSuccessorsAndRadiusUpto {
  my($This, $StartVertexID, $Radius) = @_;
  my($PathsTraversal);

  $PathsTraversal = new Graph::PathsTraversal($This);
  $PathsTraversal->PerformNeighborhoodVerticesSearchWithSuccessorsAndRadiusUpto($StartVertexID, $Radius);

  return $PathsTraversal->GetVerticesNeighborhoodsWithSuccessors();
}

# Get neighborhood vertices around a specified vertex, along with their successor connected vertices, collected
# at all neighborhood radii as a list containing references to lists with first value corresponding to vertex
# ID and second value as reference to a list containing its successor connected vertices.
#
# For a neighborhood vertex at each radius level, the successor connected vertices correspond to the
# neighborhood vertices at the next radius level. Consequently, the neighborhood vertices at the last
# radius level don't contain any successor vertices which fall outside the range of specified radius.
#
sub GetNeighborhoodVerticesWithSuccessors {
  my($This, $StartVertexID) = @_;
  my($PathsTraversal);

  $PathsTraversal = new Graph::PathsTraversal($This);
  $PathsTraversal->PerformNeighborhoodVerticesSearchWithSuccessors($StartVertexID);

  return $PathsTraversal->GetVerticesNeighborhoodsWithSuccessors();
}

# Get adjacency matrix for the graph as a Matrix object with row and column indices
# corresponding to graph vertices returned by GetVertices method.
#
# For a simple graph G with n vertices, the adjacency matrix for G is a n x n square matrix and
# its elements Mij are:
#
#   . 0    if i == j
#   . 1    if i != j and vertex Vi is adjacent to vertex Vj
#   . 0    if i != j and vertex Vi is not adjacent to vertex Vj
#
sub GetAdjacencyMatrix {
  my($This) = @_;
  my($GraphMatrix);

  $GraphMatrix = new Graph::GraphMatrix($This);
  $GraphMatrix->GenerateAdjacencyMatrix();

  return $GraphMatrix->GetMatrix();
}

# Get Siedel adjacency matrix for the graph as a Matrix object with row and column indices
# corresponding to graph vertices returned by GetVertices method.
#
# For a simple graph G with n vertices, the Siedal adjacency matrix for G is a n x n square matrix and
# its elements Mij are:
#
#   . 0    if i == j
#   . -1   if i != j and vertex Vi is adjacent to vertex Vj
#   . 1    if i != j and vertex Vi is not adjacent to vertex Vj
#
sub GetSiedelAdjacencyMatrix {
  my($This) = @_;
  my($GraphMatrix);

  $GraphMatrix = new Graph::GraphMatrix($This);
  $GraphMatrix->GenerateSiedelAdjacencyMatrix();

  return $GraphMatrix->GetMatrix();
}

# Get distance matrix for the graph as a Matrix object with row and column indices
# corresponding to graph vertices returned by GetVertices method.
#
# For a simple graph G with n vertices, the distance matrix for G is a n x n square matrix and
# its elements Mij are:
#
#   . 0    if i == j
#   . d    if i != j and d is the shortest distance between vertex Vi and vertex Vj
#
# Note:
#   . In the final matrix, BigNumber values correspond to vertices with no edges.
#
sub GetDistanceMatrix {
  my($This) = @_;
  my($GraphMatrix);

  $GraphMatrix = new Graph::GraphMatrix($This);
  $GraphMatrix->GenerateDistanceMatrix();

  return $GraphMatrix->GetMatrix();
}

# Get incidence matrix for the graph as a Matrix object with row and column indices
# corresponding to graph vertices and edges returned by GetVertices and GetEdges
# methods respectively.
#
# For a simple graph G with n vertices and e edges, the incidence matrix for G is a n x e matrix
# its elements Mij are:
#
#   . 1    if vertex Vi and the edge Ej are incident; in other words, Vi and Ej are related
#   . 0    otherwise
#
sub GetIncidenceMatrix {
  my($This) = @_;
  my($GraphMatrix);

  $GraphMatrix = new Graph::GraphMatrix($This);
  $GraphMatrix->GenerateIncidenceMatrix();

  return $GraphMatrix->GetMatrix();
}

# Get degree matrix for the graph as a Matrix object with row and column indices
# corresponding to graph vertices returned by GetVertices method.
#
# For a simple graph G with n vertices, the degree matrix for G is a n x n square matrix and
# its elements Mij are:
#
#   . deg(Vi)   if i == j and deg(Vi) is the degree of vertex Vi
#   . 0         otherwise
#
sub GetDegreeMatrix {
  my($This) = @_;
  my($GraphMatrix);

  $GraphMatrix = new Graph::GraphMatrix($This);
  $GraphMatrix->GenerateDegreeMatrix();

  return $GraphMatrix->GetMatrix();
}

# Get Laplacian matrix for the graph as a Matrix object with row and column indices
# corresponding to graph vertices returned by GetVertices method.
#
# For a simple graph G with n vertices, the Laplacian matrix for G is a n x n square matrix and
# its elements Mij are:
#
#   . deg(Vi)   if i == j and deg(Vi) is the degree of vertex Vi
#   . -1        if i != j and vertex Vi is adjacent to vertex Vj
#   . 0         otherwise
#
# Note: The Laplacian matrix is the difference between the degree matrix and adjacency matrix.
#
sub GetLaplacianMatrix {
  my($This) = @_;
  my($GraphMatrix);

  $GraphMatrix = new Graph::GraphMatrix($This);
  $GraphMatrix->GenerateLaplacianMatrix();

  return $GraphMatrix->GetMatrix();
}

# Get normalized Laplacian matrix for the graph as a Matrix object with row and column indices
# corresponding to graph vertices returned by GetVertices method.
#
# For a simple graph G with n vertices, the normalized Laplacian matrix L for G is a n x n square matrix and
# its elements Lij are:
#
#   . 1                           if i == j and deg(Vi) != 0
#   . -1/SQRT(deg(Vi) * deg(Vj))  if i != j and vertex Vi is adjacent to vertex Vj
#   . 0                           otherwise
#
#
sub GetNormalizedLaplacianMatrix {
  my($This) = @_;
  my($GraphMatrix);

  $GraphMatrix = new Graph::GraphMatrix($This);
  $GraphMatrix->GenerateNormalizedLaplacianMatrix();

  return $GraphMatrix->GetMatrix();
}

# Get admittance matrix for the graph as a Matrix object with row and column indices
# corresponding to graph vertices returned by GetVertices method.
#
sub GetAdmittanceMatrix {
  my($This) = @_;

  return $This->GetLaplacianMatrix();
}

# Get Kirchhoff matrix for the graph as a Matrix object with row and column indices
# corresponding to graph vertices returned by GetVertices method.
#
sub GetKirchhoffMatrix {
  my($This) = @_;

  return $This->GetLaplacianMatrix();
}

1;

__END__

=head1 NAME

Graph

=head1 SYNOPSIS

use Graph;

use Graph qw(:all);

=head1 DESCRIPTION

B<Graph> class provides the following methods:

new, AddCycle, AddEdge, AddEdges, AddPath, AddVertex, AddVertices, ClearCycles,
Copy, CopyEdgesProperties, CopyVerticesAndEdges, CopyVerticesProperties,
DeleteCycle, DeleteEdge, DeleteEdgeProperties, DeleteEdgeProperty, DeleteEdges,
DeleteEdgesProperties, DeleteEdgesProperty, DeleteGraphProperties,
DeleteGraphProperty, DeletePath, DeleteVertex, DeleteVertexProperties,
DeleteVertexProperty, DeleteVertices, DeleteVerticesProperty, DetectCycles,
GetAdjacencyMatrix, GetAdmittanceMatrix, GetAllPaths, GetAllPathsStartingAt,
GetAllPathsStartingAtWithLength, GetAllPathsStartingAtWithLengthUpto,
GetAllPathsWithLength, GetAllPathsWithLengthUpto, GetCircumference,
GetConnectedComponentsVertices, GetCycles, GetCyclesWithEvenSize,
GetCyclesWithOddSize, GetCyclesWithSize, GetCyclesWithSizeGreaterThan,
GetCyclesWithSizeLessThan, GetDegree, GetDegreeMatrix, GetDistanceMatrix,
GetEdgeCycles, GetEdgeCyclesWithEvenSize, GetEdgeCyclesWithOddSize,
GetEdgeCyclesWithSize, GetEdgeCyclesWithSizeGreaterThan,
GetEdgeCyclesWithSizeLessThan, GetEdgeProperties, GetEdgeProperty, GetEdges,
GetEdgesProperty, GetFusedAndNonFusedCycles, GetGirth, GetGraphProperties,
GetGraphProperty, GetIncidenceMatrix, GetIsolatedVertices, GetKirchhoffMatrix,
GetLaplacianMatrix, GetLargestCycle, GetLargestEdgeCycle, GetLargestVertexCycle,
GetLeafVertices, GetMaximumDegree, GetMininumDegree, GetNeighborhoodVertices,
GetNeighborhoodVerticesWithRadiusUpto, GetNeighborhoodVerticesWithSuccessors,
GetNeighborhoodVerticesWithSuccessorsAndRadiusUpto, GetNeighbors,
GetNormalizedLaplacianMatrix, GetNumOfCycles, GetNumOfCyclesWithEvenSize,
GetNumOfCyclesWithOddSize, GetNumOfCyclesWithSize,
GetNumOfCyclesWithSizeGreaterThan, GetNumOfCyclesWithSizeLessThan,
GetNumOfEdgeCycles, GetNumOfEdgeCyclesWithEvenSize, GetNumOfEdgeCyclesWithOddSize,
GetNumOfEdgeCyclesWithSize, GetNumOfEdgeCyclesWithSizeGreaterThan,
GetNumOfEdgeCyclesWithSizeLessThan, GetNumOfVertexCycles,
GetNumOfVertexCyclesWithEvenSize, GetNumOfVertexCyclesWithOddSize,
GetNumOfVertexCyclesWithSize, GetNumOfVertexCyclesWithSizeGreaterThan,
GetNumOfVertexCyclesWithSizeLessThan, GetPaths, GetPathsBetween,
GetPathsStartingAt, GetPathsStartingAtWithLength,
GetPathsStartingAtWithLengthUpto, GetPathsWithLength, GetPathsWithLengthUpto,
GetSiedelAdjacencyMatrix, GetSizeOfLargestCycle, GetSizeOfLargestEdgeCycle,
GetSizeOfLargestVertexCycle, GetSizeOfSmallestCycle, GetSizeOfSmallestEdgeCycle,
GetSizeOfSmallestVertexCycle, GetSmallestCycle, GetSmallestEdgeCycle,
GetSmallestVertexCycle, GetTopologicallySortedVertices, GetVertex,
GetVertexCycles, GetVertexCyclesWithEvenSize, GetVertexCyclesWithOddSize,
GetVertexCyclesWithSize, GetVertexCyclesWithSizeGreaterThan,
GetVertexCyclesWithSizeLessThan, GetVertexProperties, GetVertexProperty,
GetVertexWithLargestDegree, GetVertexWithSmallestDegree, GetVertices,
GetVerticesProperty, GetVerticesWithDegreeLessThan, HasCycle, HasEdge,
HasEdgeProperty, HasEdges, HasFusedCycles, HasGraphProperty, HasPath, HasVertex,
HasVertexProperty, HasVertices, IsAcyclic, IsAcyclicEdge, IsAcyclicVertex,
IsCyclic, IsCyclicEdge, IsCyclicVertex, IsGraph, IsIsolatedVertex, IsLeafVertex,
IsUnicyclic, IsUnicyclicEdge, IsUnicyclicVertex, SetActiveCyclicPaths,
SetEdgeProperties, SetEdgeProperty, SetEdgesProperty, SetGraphProperties,
SetGraphProperty, SetVertexProperties, SetVertexProperty, SetVerticesProperty,
StringifyEdgesProperties, StringifyGraph, StringifyGraphProperties,
StringifyProperties, StringifyVerticesAndEdges, StringifyVerticesProperties,
UpdateEdgeProperty, UpdateVertexProperty

=head2 METHODS

=over 4

=item B<new>

    $NewGraph = new Graph([@VertexIDs]);

Using specified I<Graph> I<VertexIDs>, B<new> method creates a new B<Graph> object and returns
newly created B<Graph> object.

Examples:

    $Graph = new Graph();
    $Graph = new Graph(@VertexIDs);

=item B<AddCycle>

    $Graph->AddCycle(@VertexIDs);

Adds edges between successive pair of I<VertexIDs> including an additional edge from the last
to first vertex ID to complete the cycle to I<Graph> and returns I<Graph>.

=item B<AddEdge>

    $Graph->AddEdge($VertexID1, $VertexID2);

Adds an edge between I<VertexID1> and I<VertexID2> in a I<Graph> and returns I<Graph>.

=item B<AddEdges>

    $Graph->AddEdges(@VertexIDs);

Adds edges between successive pair of I<VertexIDs> in a I<Graph> and returns I<Graph>.

=item B<AddPath>

    $Graph->AddPath(@VertexIDs);

Adds edges between successive pair of I<VertexIDs> in a I<Graph> and returns I<Graph>.

=item B<AddVertex>

    $Graph->AddVertex($VertexID);

Adds I<VertexID> to I<Graph> and returns I<Graph>.

=item B<AddVertices>

    $Graph->AddVertices(@VertexIDs);

Adds vertices using I<VertexIDs> to I<Graph> and returns I<Graph>.

=item B<ClearCycles>

    $Graph->ClearCycles();

Delete all cycle properties assigned to graph, vertices, and edges by I<DetectCycles> method.

=item B<Copy>

    $NewGraph = $Graph->Copy();

Copies I<Graph> and its associated data using B<Storable::dclone> and returns a new
B<Graph> object.

=item B<CopyEdgesProperties>

    $OtherGraph = $Graph->CopyEdgesProperties($OtherGraph);

Copies all properties associated with edges from I<Graph> to I<$OtherGraph> and
returns I<OtherGraph>.

=item B<CopyVerticesAndEdges>

    $OtherGraph = $Graph->CopyVerticesAndEdges($OtherGraph);

Copies all vertices and edges from I<Graph> to I<$OtherGraph> and returns I<OtherGraph>.

=item B<CopyVerticesProperties>

    $OtherGraph = $Graph->CopyVerticesProperties($OtherGraph);

Copies all properties associated with vertices from I<Graph> to I<$OtherGraph> and
returns I<OtherGraph>.

=item B<DeleteCycle>

    $Graph->DeleteCycle(@VertexIDs);

Deletes edges between successive pair of I<VertexIDs> including an additional edge from the last
to first vertex ID to complete the cycle to I<Graph> and returns I<Graph>.

=item B<DeleteEdge>

    $Graph->DeleteEdge($VertexID1, $VertexID2);

Deletes an edge between I<VertexID1> and I<VertexID2> in a I<Graph> and returns I<Graph>.

=item B<DeleteEdgeProperties>

    $Graph->DeleteEdgeProperties($VertexID1, $VertexID2);

Deletes all properties associated with edge between I<VertexID1> and I<VertexID2> in a I<Graph>
and returns I<Graph>.

=item B<DeleteEdgeProperty>

    $Graph->DeleteEdgeProperty($PropertyName, $VertexID1, $VertexID2);

Deletes I<PropertyName> associated with edge between I<VertexID1> and I<VertexID2> in a I<Graph>
and returns I<Graph>.

=item B<DeleteEdges>

    $Graph->DeleteEdges(@VertexIDs);

Deletes edges between successive pair of I<VertexIDs> and returns I<Graph>.

=item B<DeleteEdgesProperties>

    $Graph->DeleteEdgesProperties(@VertexIDs);

Deletes all properties associated with edges between successive pair of I<VertexIDs> and
returns I<Graph>.

=item B<DeleteEdgesProperty>

    $Graph->DeleteEdgesProperty($PropertyName, @VertexIDs);

Deletes I<PropertyName> associated with edges between successive pair of I<VertexIDs>
and returns I<Graph>.

=item B<DeleteGraphProperties>

    $Graph->DeleteGraphProperties();

Deletes all properties associated as graph not including properties associated to vertices
or edges and returns I<Graph>.

=item B<DeleteGraphProperty>

    $Graph->DeleteGraphProperty($PropertyName);

Deletes a I<PropertyName> associated as graph property and returns I<Graph>.

=item B<DeletePath>

    $Graph->DeletePath(@VertexIDs);

Deletes edges between successive pair of I<VertexIDs> in a I<Graph> and returns I<Graph>.

=item B<DeleteVertex>

    $Graph->DeleteVertex($VertexID);

Deletes I<VertexID> to I<Graph> and returns I<Graph>.

=item B<DeleteVertexProperties>

    $Graph->DeleteVertexProperties($VertexID);

Deletes all properties associated with I<VertexID> and returns I<Graph>.

=item B<DeleteVertexProperty>

    $Graph->DeleteVertexProperty($PropertyName, $VertexID);

Deletes a I<PropertyName> associated with I<VertexID> and returns I<Graph>.

=item B<DeleteVertices>

    $Graph->DeleteVertices(@VertexIDs);

Deletes vertices specified in I<VertexIDs> and returns I<Graph>.

=item B<DeleteVerticesProperty>

    $Graph->DeleteVerticesProperty($PropertyName, @VertexIDs);

Deletes a I<PropertyName> associated with I<VertexIDs> and returns I<Graph>.

=item B<DetectCycles>

    $Graph->DetectCycles();

Detect cycles using B<CyclesDetection> class and associate found cycles to I<Graph>
object as graph properties: I<ActiveCyclicPaths, AllCyclicPaths, IndependentCyclicPaths>.

Notes:

    . CyclesDetection class detects all cycles in the graph and filters
      them to find independent cycles.
    . All cycles related methods in the graph operate on
      ActiveCyclicPaths. By default, active cyclic paths correspond
      to IndependentCyclicPaths. This behavior can be changed
      using SetActiveCyclicPaths method.

=item B<GetAdjacencyMatrix>

    $GraphMatrix = $Graph->GetAdjacencyMatrix();

Returns adjacency matrix for I<Graph> as a I<GraphMatrix> object with row and column indices
corresponding to graph vertices returned by GetVertices method.

For a simple graph G with n vertices, the adjacency matrix for G is a n x n square matrix and
its elements Mij are:

    . 0    if i == j
    . 1    if i != j and vertex Vi is adjacent to vertex Vj
    . 0    if i != j and vertex Vi is not adjacent to vertex Vj

=item B<GetAdmittanceMatrix>

    $GraphMatrix = $Graph->GetAdmittanceMatrix();

Returns admittance matrix for I<Graph> as a I<GraphMatrix> object with row and column indices
corresponding to graph vertices returned by GetVertices method.

For a simple graph G with n vertices, the adjacency matrix for G is a n x n square matrix and
its elements Mij are:

    . 0    if i == j
    . 1    if i != j and vertex Vi is adjacent to vertex Vj
    . 0    if i != j and vertex Vi is not adjacent to vertex Vj

=item B<GetAllPaths>

    $PathsRef = $Graph->GetAllPaths([$AllowCycles]);

Returns a reference to an array containing B<Path> objects corresponding to all possible
lengths starting from each vertex in graph with sharing of edges in paths traversed.
By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
completing the cycle. Duplicate paths are not removed.

=item B<GetAllPathsStartingAt>

    @Paths = $Graph->GetAllPathsStartingAt($StartVertexID,
             [$AllowCycles]);

Returns an array of I<Path> objects starting from a I<StartVertexID> of any length
with sharing of edges in paths traversed. By default, cycles are included in paths.
A path containing a cycle is terminated at a vertex completing the cycle.

=item B<GetAllPathsStartingAtWithLength>

    @Paths = $Graph->GetAllPathsStartingAtWithLength($StartVertexID,
             $Length, [$AllowCycles]);

Returns an array of I<Path> objects starting from a I<StartVertexID> of specified I<Length>
with sharing of edges in paths traversed. By default, cycles are included in paths.
A path containing a cycle is terminated at a vertex completing the cycle.

=item B<GetAllPathsStartingAtWithLengthUpto>

    @Paths = $Graph->GetAllPathsStartingAtWithLengthUpto($StartVertexID,
             $Length, [$AllowCycles]);

Returns an array of I<Path> objects starting from a I<StartVertexID> with length upto a
I<Length> with sharing of edges in paths traversed. By default, cycles are included in paths.
A path containing a cycle is terminated at a vertex completing the cycle.

=item B<GetAllPathsWithLength>

    $PathsRef = $Graph->GetAllPathsWithLength($Length,
                [$AllowCycles]);

Returns a reference to an array containing B<Path> objects corresponding to paths with
I<Length> starting from each vertex in graph with sharing of edges in paths traversed.
By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
completing the cycle. Duplicate paths are not removed.

=item B<GetAllPathsWithLengthUpto>

    $PathsRef = $Graph->GetAllPathsWithLengthUpto($Length,
                [$AllowCycles]);

Returns a reference to an array containing B<Path> objects corresponding to paths up to
specified I<Length> starting from each vertex in graph with sharing of edges in paths traversed.
By default, cycles are included in paths. A path containing a cycle is terminated at a vertex
completing the cycle. Duplicate paths are not removed.

=item B<GetCircumference>

    $Circumference = $Graph->GetCircumference();

Returns size of largest cycle in a I<Graph>

=item B<GetConnectedComponentsVertices>

    @ConnectedComponents = $Graph->GetConnectedComponentsVertices();

Returns an array I<ConnectedComponents> containing referecens to arrays with vertex
IDs for each component sorted in order of their decreasing size.

=item B<GetCycles>

    @CyclicPaths = $Graphs->GetCycles();

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles
in a I<Graph>.

=item B<GetCyclesWithEvenSize>

    @CyclicPaths = $Graph->GetCyclesWithEvenSize();

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
even size in a I<Graph>.

=item B<GetCyclesWithOddSize>

    @CyclicPaths = $Graph->GetCyclesWithOddSize();

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
odd size in a I<Graph>.

=item B<GetCyclesWithSize>

    @CyclicPaths = $Graph->GetCyclesWithSize($CycleSize);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
I<CycleSize> in a I<Graph>.

=item B<GetCyclesWithSizeGreaterThan>

    @CyclicPaths = $Graph->GetCyclesWithSizeGreaterThan($CycleSize);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
size greater than I<CycleSize> in a I<Graph>.

=item B<GetCyclesWithSizeLessThan>

    @CyclicPaths = $Graph->GetCyclesWithSizeGreaterThan($CycleSize);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
size less than I<CycleSize> in a I<Graph>.

=item B<GetDegree>

    $Degree = $Graph->GetDegree($VertexID);

Returns B<Degree> for I<VertexID> in a I<Graph> corresponding to sum of in and out vertex
degree values.

=item B<GetDegreeMatrix>

    $GraphMatrix = $Graph->GetDegreeMatrix();

Returns degree matrix for I<Graph> as a I<GraphMatrix> object with row and column indices
corresponding to graph vertices returned by GetVertices method.

For a simple graph G with n vertices, the degree matrix for G is a n x n square matrix and
its elements Mij are:

    . deg(Vi)   if i == j and deg(Vi) is the degree of vertex Vi
    . 0         otherwise

=item B<GetDistanceMatrix>

    $GraphMatrix = $Graph->GetDistanceMatrix();

Returns distance matrix for I<Graph> as a I<GraphMatrix> object with row and column indices
corresponding to graph vertices returned by GetVertices method.

For a simple graph G with n vertices, the distance matrix for G is a n x n square matrix and
its elements Mij are:

    . 0    if i == j
    . d    if i != j and d is the shortest distance between vertex Vi and vertex Vj

In the final matrix, value of constant B<BigNumber> defined in B<Constants.pm> module
corresponds to vertices with no edges.

=item B<GetEdgeCycles>

    @CyclicPaths = $Graph->GetEdgeCycles($VertexID1, $VertexID2);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to all cycles containing
edge between I<VertexID1> and I<VertexID2> in a I<Graph>.

=item B<GetEdgeCyclesWithEvenSize>

    @CyclicPaths = $Graph->GetEdgeCyclesWithEvenSize($VertexID1,
                   $VertexID2);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
even size containing edge between I<VertexID1> and I<VertexID2> in a I<Graph>.

=item B<GetEdgeCyclesWithOddSize>

    @CyclicPaths = $Graph->GetEdgeCyclesWithOddSize($VertexID1,
                   $VertexID2);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
odd size containing edge between I<VertexID1> and I<VertexID2> in a I<Graph>.

=item B<GetEdgeCyclesWithSize>

    @CyclicPaths = $Graph->GetEdgeCyclesWithSize($VertexID1, $VertexID2,
                   $CycleSize);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
size I<CycleSize> containing edge between I<VertexID1> and I<VertexID2> in a I<Graph>.

=item B<GetEdgeCyclesWithSizeGreaterThan>

    @CyclicPaths = $Graph->GetEdgeCyclesWithSizeGreaterThan($VertexID1,
                   $VertexID2, $CycleSize);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
size greater than I<CycleSize> containing edge between I<VertexID1> and I<VertexID2>
in a I<Graph>.

=item B<GetEdgeCyclesWithSizeLessThan>

    @CyclicPaths = $Graph->GetEdgeCyclesWithSizeLessThan($VertexID1,
                   $VertexID2, $CycleSize);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
size less than I<CycleSize> containing edge between I<VertexID1> and I<VertexID2>.

=item B<GetEdgeProperties>

    %EdgeProperties = $Graph->GetEdgeProperties($VertexID1, $VertexID2);

Returns a hash B<EdgeProperties> containing all B<PropertyName> and B<PropertyValue>
pairs associated with an edge between I<VertexID1> and I<VertexID2> in a I<Graph>.

=item B<GetEdgeProperty>

    $Value = $Graph->GetEdgeProperty($PropertyName, $VertexID1, $VertexID2);

Returns value of I<PropertyName> associated with an edge between I<VertexID1>
and I<VertexID2> in a I<Graph>.

=item B<GetEdges>

    @EdgeVertexIDs = $Graph->GetEdges($VertexID);
    $NumOfEdges = $Graph->GetEdges($VertexID);

Returns an array I<EdgeVertexIDs> with successive pair of IDs corresponding to edges involving
I<VertexID> or number of edges for I<VertexID> in a I<Graph>.

=item B<GetEdgesProperty>

    @PropertyValues = $Graph->GetEdgesProperty($PropertyName, @VertexIDs);

Returns an array I<PropertyValues> containing property values corresponding to
I<PropertyName> associated with edges between successive pair of I<VertexIDs>.

=item B<GetFusedAndNonFusedCycles>

    ($FusedCycleSetsRef, $NonFusedCyclesRef) =
       $Graph->GetFusedAndNonFusedCycles();

Returns references to arrays I<FusedCycleSetsRef> and I<NonFusedCyclesRef>
containing references to arrays of cyclic I<Path> objects corresponding to fuses and
non-fused cyclic paths.

=item B<GetGirth>

    $Girth = $Graph->GetGirth();

Returns size of smallest cycle in a I<Graph>.

=item B<GetGraphProperties>

    %GraphProperties = $Graph->GetGraphProperties();

Returns a hash B<EdgeProperties> containing all B<PropertyName> and B<PropertyValue>
pairs associated with graph in a I<Graph>.

=item B<GetGraphProperty>

    $PropertyValue = $Graph->GetGraphProperty($PropertyName);

Returns value of I<PropertyName> associated with graph in a I<Graph>.

=item B<GetIncidenceMatrix>

    $GraphMatrix = $Graph->GetIncidenceMatrix();

Returns incidence matrix for I<Graph> as a I<GraphMatrix> object with row and column indices
corresponding to graph vertices returned by GetVertices method.

For a simple graph G with n vertices and e edges, the incidence matrix for G is a n x e matrix
its elements Mij are:

    . 1    if vertex Vi and the edge Ej are incident; in other words, Vi and Ej are related
    . 0    otherwise

=item B<GetIsolatedVertices>

    @VertexIDs = $Graph->GetIsolatedVertices();

Returns an array I<VertexIDs> containing vertices without any edges in I<Graph>.

=item B<GetKirchhoffMatrix>

    $GraphMatrix = $Graph->GetGetKirchhoffMatrix();

Returns Kirchhoff matrix for I<Graph> as a I<GraphMatrix> object with row and column indices
corresponding to graph vertices returned by GetVertices method.

B<KirchhoffMatrix> is another name for B<LaplacianMatrix>.

=item B<GetLaplacianMatrix>

    $GraphMatrix = $Graph->GetLaplacianMatrix();

Returns Laplacian matrix for I<Graph> as a I<GraphMatrix> object with row and column indices
corresponding to graph vertices returned by GetVertices method.

For a simple graph G with n vertices, the Laplacian matrix for G is a n x n square matrix and
its elements Mij are:

    . deg(Vi)   if i == j and deg(Vi) is the degree of vertex Vi
    . -1        if i != j and vertex Vi is adjacent to vertex Vj
    . 0         otherwise

=item B<GetLargestCycle>

    $CyclicPath = $Graph->GetLargestCycle();

Returns a cyclic I<Path> object corresponding to largest cycle in a I<Graph>.

=item B<GetLargestEdgeCycle>

    $CyclicPath = $Graph->GetLargestEdgeCycle($VertexID1, $VertexID2);

Returns a cyclic I<Path> object corresponding to largest cycle containing edge between
I<VertexID1> and I<VertexID2> in a I<Graph>.

=item B<GetLargestVertexCycle>

    $CyclicPath = $Graph->GetLargestVertexCycle($VertexID);

Returns a cyclic I<Path> object corresponding to largest cycle containing I<VertexID>
in a I<Graph>.

=item B<GetLeafVertices>

    @VertexIDs = $Graph->GetLeafVertices();

Returns an array I<VertexIDs> containing vertices with degree of 1 in a I<Graph>.

=item B<GetMaximumDegree>

    $Degree = $Graph->GetMaximumDegree();

Returns value of maximum vertex degree in a I<Graph>.

=item B<GetMininumDegree>

    $Degree = $Graph->GetMininumDegree();

Returns value of minimum vertex degree in a I<Graph>.

=item B<GetNeighborhoodVertices>

    @VertexNeighborhoods = GetNeighborhoodVertices($StartVertexID);

Returns an array I<VertexNeighborhoods> containing references to arrays corresponding to
neighborhood vertices around a specified I<StartVertexID> at all possible radii levels.

=item B<GetNeighborhoodVerticesWithRadiusUpto>

    @VertexNeighborhoods = GetNeighborhoodVerticesWithRadiusUpto(
                           $StartVertexID, $Radius);

Returns an array I<VertexNeighborhoods> containing references to arrays corresponding to
neighborhood vertices around a specified I<StartVertexID> upto specified I<Radius> levels.

=item B<GetNeighborhoodVerticesWithSuccessors>

    @VertexNeighborhoods = GetNeighborhoodVerticesWithSuccessors(
                           $StartVertexID);

Returns vertex neighborhoods around a specified I<StartVertexID>, along with their successor
connected vertices, collected at all neighborhood radii as an array I<VertexNeighborhoods>
containing references to arrays with first value corresponding to vertex ID and second
value as reference to an array  containing its successor connected vertices.

For a neighborhood vertex at each radius level, the successor connected vertices correspond to the
neighborhood vertices at the next radius level. Consequently, the neighborhood vertices at the last
radius level don't contain any successor vertices which fall outside the range of specified radius.

=item B<GetNeighborhoodVerticesWithSuccessorsAndRadiusUpto>

    @VertexNeighborhoods = GetNeighborhoodVerticesWithSuccessors(
                           $StartVertexID, $Radius);

Returns vertex neighborhoods around a specified I<StartVertexID>, along with their successor
connected vertices, collected with in a specified I<Radius> as an array I<VertexNeighborhoods>
containing references to arrays with first value corresponding to vertex ID and second value
as reference to a list containing its successor connected vertices.

For a neighborhood vertex at each radius level, the successor connected vertices correspond to the
neighborhood vertices at the next radius level. Consequently, the neighborhood vertices at the last
radius level don't contain any successor vertices which fall outside the range of specified radius.

=item B<GetNeighbors>

    @VertexIDs = $Graph->GetNeighbors($VertexID);
    $NumOfNeighbors = $Graph->GetNeighbors($VertexID);

Returns an array I<VertexIDs> containing vertices connected to I<VertexID> of number of
neighbors of a I<VertextID> in a I<Graph>.

=item B<GetNormalizedLaplacianMatrix>

    $GraphMatrix = $Graph->GetNormalizedLaplacianMatrix();

Returns normalized Laplacian matrix for I<Graph> as a I<GraphMatrix> object with row and column indices
corresponding to graph vertices returned by GetVertices method.

For a simple graph G with n vertices, the normalized Laplacian matrix L for G is a n x n square
matrix and its elements Lij are:

    .  1                           if i == j and deg(Vi) != 0
    .  -1/SQRT(deg(Vi) * deg(Vj))  if i != j and vertex Vi is adjacent to vertex Vj
    .  0                           otherwise

=item B<GetNumOfCycles>

    $NumOfCycles = $Graph->GetNumOfCycles();

Returns number of cycles in a I<Graph>.

=item B<GetNumOfCyclesWithEvenSize>

    $NumOfCycles = $Graph->GetNumOfCyclesWithEvenSize();

Returns number of cycles with even size in a I<Graph>.

=item B<GetNumOfCyclesWithOddSize>

    $NumOfCycles = $Graph->GetNumOfCyclesWithOddSize();

Returns number of cycles with odd size in a I<Graph>.

=item B<GetNumOfCyclesWithSize>

    $NumOfCycles = $Graph->GetNumOfCyclesWithSize($CycleSize);

Returns number of cycles with I<CyclesSize> in a I<Graph>.

=item B<GetNumOfCyclesWithSizeGreaterThan>

    $NumOfCycles = $Graph->GetNumOfCyclesWithSizeGreaterThan(
                   $CycleSize);

Returns number of cycles with size greater than I<CyclesSize> in a I<Graph>.

=item B<GetNumOfCyclesWithSizeLessThan>

    $NumOfCycles = $Graph->GetNumOfCyclesWithSizeLessThan($CycleSize);

Returns number of cycles with size less than I<CyclesSize> in a I<Graph>.

=item B<GetNumOfEdgeCycles>

    $NumOfCycles = $Graph->GetNumOfEdgeCycles($VertexID1, $VertexID2);

Returns number of cycles containing edge between I<VertexID1> and I<VertexID2>
in a I<Graph>.

=item B<GetNumOfEdgeCyclesWithEvenSize>

    $NumOfCycles = $Graph->GetNumOfEdgeCyclesWithEvenSize($VertexID1,
                   $VertexID2);

Returns number of cycles containing edge between I<VertexID1> and I<VertexID2> with even
size in a I<Graph>.

=item B<GetNumOfEdgeCyclesWithOddSize>

    $NumOfCycles = $Graph->GetNumOfEdgeCyclesWithOddSize($VertexID1,
                   $VertexID2);

Returns number of cycles containing edge between I<VertexID1> and I<VertexID2> with odd
size in a I<Graph>.

=item B<GetNumOfEdgeCyclesWithSize>

    $NumOfCycles = $Graph->GetNumOfEdgeCyclesWithSize($VertexID1,
                   $VertexID2, $CycleSize);

Returns number of cycles containing edge between I<VertexID1> and I<VertexID2> with
I<CycleSize> size in a I<Graph>.

=item B<GetNumOfEdgeCyclesWithSizeGreaterThan>

    $NumOfCycles = $Graph->GetNumOfEdgeCyclesWithSizeGreaterThan(
                   $VertexID1, $VertexID2, $CycleSize);

Returns number of cycles containing edge between I<VertexID1> and I<VertexID2> with
size greater than I<CycleSize> size in a I<Graph>.

=item B<GetNumOfEdgeCyclesWithSizeLessThan>

    $NumOfCycles = $Graph->GetNumOfEdgeCyclesWithSizeLessThan(
                   $VertexID1, $VertexID2, $CycleSize);

Returns number of cycles containing edge between I<VertexID1> and I<VertexID2> with
size less than I<CycleSize> size in a I<Graph>.

=item B<GetNumOfVertexCycles>

    $NumOfCycles = $Graph->GetNumOfVertexCycles($VertexID);

Returns number of cycles containing I<VertexID> in a I<Graph>.

=item B<GetNumOfVertexCyclesWithEvenSize>

    $NumOfCycles = $Graph->GetNumOfVertexCyclesWithEvenSize($VertexID);

Returns number of cycles containing I<VertexID> with even size in a I<Graph>.

=item B<GetNumOfVertexCyclesWithOddSize>

    $NumOfCycles = $Graph->GetNumOfVertexCyclesWithOddSize($VertexID);

Returns number of cycles containing I<VertexID> with odd size in a I<Graph>.

=item B<GetNumOfVertexCyclesWithSize>

    $NumOfCycles = $Graph->GetNumOfVertexCyclesWithSize($VertexID);

Returns number of cycles containing I<VertexID> with even size in a I<Graph>.

=item B<GetNumOfVertexCyclesWithSizeGreaterThan>

    $NumOfCycles = $Graph->GetNumOfVertexCyclesWithSizeGreaterThan(
                   $VertexID, $CycleSize);

Returns number of cycles containing I<VertexID> with size greater than I<CycleSize>
in a I<Graph>.

=item B<GetNumOfVertexCyclesWithSizeLessThan>

    $NumOfCycles = $Graph->GetNumOfVertexCyclesWithSizeLessThan(
                   $VertexID, $CycleSize);

Returns number of cycles containing I<VertexID> with size less than I<CycleSize>
in a I<Graph>.

=item B<GetPaths>

    $PathsRefs = $Graph->GetPaths([$AllowCycles]);

Returns a reference to an array of I<Path> objects corresponding to paths of all possible
lengths starting from each vertex with no sharing of edges in paths traversed. By default,
cycles are included in paths. A path containing a cycle is terminated at a vertex completing
the cycle.

=item B<GetPathsBetween>

    @Paths = $Graph->GetPathsBetween($StartVertexID, $EndVertexID);

Returns an arrays of I<Path> objects list of paths between I<StartVertexID> and I<EndVertexID>.
For cyclic graphs, the list contains may contain more than one I<Path> object.

=item B<GetPathsStartingAt>

    @Paths = $Graph->GetPathsStartingAt($StartVertexID, [$AllowCycles]);

Returns an array of I<Path> objects corresponding to all possible lengths starting from a
specified I<StartVertexID> with no sharing of edges in paths traversed. By default, cycles
are included in paths. A path containing a cycle is terminated at a vertex completing the cycle.

=item B<GetPathsStartingAtWithLength>

    @Paths = $Graph->StartingAtWithLength($StartVertexID, $Length,
             $AllowCycles);

Returns an array of I<Path> objects corresponding to all paths starting from a specified I<StartVertexID>
with length I<Length> and no sharing of edges in paths traversed. By default, cycles are included in paths.
A path containing a cycle is terminated at a vertex completing the cycle.

=item B<GetPathsStartingAtWithLengthUpto>

    @Paths = $Graph->StartingAtWithLengthUpto($StartVertexID, $Length,
             $AllowCycles);

Returns an array of I<Path> objects corresponding to all paths starting from a specified I<StartVertexID>
with length upto I<Length> and no sharing of edges in paths traversed. By default, cycles are included in paths.
A path containing a cycle is terminated at a vertex completing the cycle.

=item B<GetPathsWithLength>

    @Paths = $Graph->GetPathsWithLength($Length, $AllowCycles);

Returns an array of I<Path> objects corresponding to to paths starting from each vertex in graph
with specified <Length> and no sharing of edges in paths traversed. By default, cycles are included
in paths. A path containing a cycle is terminated at a vertex completing the cycle.

=item B<GetPathsWithLengthUpto>

    @Paths = $Graph->GetPathsWithLengthUpto($Length, $AllowCycles);

Returns an array of I<Path> objects corresponding to to paths starting from each vertex in graph
with length upto specified I<Length> and no sharing of edges in paths traversed. By default,
cycles are included in paths. A path containing a cycle is terminated at a vertex completing the cycle.

=item B<GetSiedelAdjacencyMatrix>

    $GraphMatrix = $Graph->GetSiedelAdjacencyMatrix();

Returns Siedel admittance matrix for I<Graph> as a I<GraphMatrix> object with row and column indices
corresponding to graph vertices returned by GetVertices method.

For a simple graph G with n vertices, the Siedal adjacency matrix for G is a n x n square matrix and
its elements Mij are:

    . 0    if i == j
    . -1   if i != j and vertex Vi is adjacent to vertex Vj
    . 1    if i != j and vertex Vi is not adjacent to vertex Vj

=item B<GetSizeOfLargestCycle>

    $Size = $Graph->GetSizeOfLargestCycle();

Returns size of the largest cycle in a I<Graph>.

=item B<GetSizeOfLargestEdgeCycle>

    $Size = $Graph->GetSizeOfLargestEdgeCycle($VertexID1, $VertexID2);

Returns size of the largest cycle containing egde between I<VertextID1> and  I<VertexID2>
in a I<Graph>.

=item B<GetSizeOfLargestVertexCycle>

    $Size = $Graph->GetSizeOfLargestVertexCycle($VertexID);

Returns size of the largest cycle containing  I<VertextID> in a I<Graph>.

=item B<GetSizeOfSmallestCycle>

    $Size = $Graph->GetSizeOfSmallestCycle();

Returns size of the smallest cycle in a I<Graph>.

=item B<GetSizeOfSmallestEdgeCycle>

    $Size = $Graph->GetSizeOfSmallestEdgeCycle($VertexID1, $VertexID2);

Returns size of the smallest cycle containing egde between I<VertextID1> and  I<VertexID2>
in a I<Graph>.

=item B<GetSizeOfSmallestVertexCycle>

    $Size = $Graph->GetSizeOfSmallestVertexCycle($VertexID);

Returns size of the smallest cycle containing  I<VertextID> in a I<Graph>.

=item B<GetSmallestCycle>

    $CyclicPath = $Graph->GetSmallestCycle();

Returns a cyclic I<Path> object corresponding to smallest cycle in a I<Graph>.

=item B<GetSmallestEdgeCycle>

    $CyclicPath = $Graph->GetSmallestEdgeCycle($VertexID1, $VertexID2);

Returns a cyclic I<Path> object corresponding to smallest cycle containing edge between
I<VertexID1> and I<VertexID2> in a I<Graph>.

=item B<GetSmallestVertexCycle>

    $CyclicPath = $Graph->GetSmallestVertexCycle($VertexID);

Returns a cyclic I<Path> object corresponding to smallest cycle containing I<VertexID> in a I<Graph>.

=item B<GetTopologicallySortedVertices>

    @VertexIDs = $Graph->GetTopologicallySortedVertices(
                 [$RootVertexID]);

Returns an array of I<VertexIDs> sorted topologically starting from a specified I<RootVertexID> or
from an arbitrary vertex ID.

=item B<GetVertex>

    $VertexValue = $Graph->GetVertex($VertexID);

Returns vartex value for I<VertexID> in a I<Graph>. Vartex IDs and values are equivalent
in the current implementation of B<Graph>.

=item B<GetVertexCycles>

    @CyclicPaths = $Graph->GetVertexCycles($VertexID);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to all cycles containing
I<VertexID> in a I<Graph>.

=item B<GetVertexCyclesWithEvenSize>

    @CyclicPaths = $Graph->GetVertexCyclesWithEvenSize($VertexID);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
even size containing I<VertexID> in a I<Graph>.

=item B<GetVertexCyclesWithOddSize>

    @CyclicPaths = $Graph->GetVertexCyclesWithOddSize($VertexID);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
odd size containing I<VertexID> in a I<Graph>.

=item B<GetVertexCyclesWithSize>

    @CyclicPaths = $Graph->GetVertexCyclesWithSize($VertexID,
                   $CycleSize);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
size I<CycleSize> containing I<VertexID> in a I<Graph>.

=item B<GetVertexCyclesWithSizeGreaterThan>

    @CyclicPaths = $Graph->GetVertexCyclesWithSizeGreaterThan($VertexID,
                   $CycleSize);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
size greater than I<CycleSize> containing I<VertexID> in a I<Graph>.

=item B<GetVertexCyclesWithSizeLessThan>

    @CyclicPaths = $Graph->GetVertexCyclesWithSizeLessThan($VertexID,
                   $CycleSize);

Returns an array I<CyclicPaths> containing I<Path> objects corresponding to cycles with
size less than I<CycleSize> containing I<VertexID> in a I<Graph>.

=item B<GetVertexProperties>

    %VertexProperties = $Graph->GetVertexProperties($VertexID);

Returns a hash B<VertexProperties> containing all B<PropertyName> and B<PropertyValue>
pairs associated with a I<VertexID> in a I<Graph>.

=item B<GetVertexProperty>

    $Value = $Graph->GetVertexProperty($PropertyName, $VertexID);

Returns value of I<PropertyName> associated with a I<VertexID> in a I<Graph>.

=item B<GetVertexWithLargestDegree>

    $VertexID = $Graph->GetVertexWithLargestDegree();

Returns B<VertexID> with largest degree in a I<Graph>.

=item B<GetVertexWithSmallestDegree>

    $VertexID = $Graph->GetVertexWithSmallestDegree();

Returns B<VertexID> with smallest degree in a I<Graph>.

=item B<GetVertices>

    @VertexIDs = $Graph->GetVertices();
    $VertexCount = $Graph->GetVertices();

Returns an array of I<VertexIDs> corresponding to all vertices in a I<Graph>; in a scalar context,
number of vertices is returned.

=item B<GetVerticesProperty>

    @PropertyValues = $Graph->GetVerticesProperty($PropertyName, @VertexIDs);

Returns an array I<PropertyValues> containing property values corresponding to
I<PropertyName> associated with with I<VertexIDs> in a I<Graph>.

=item B<GetVerticesWithDegreeLessThan>

    @VertexIDs  = $Graph->GetVerticesWithDegreeLessThan($Degree);

Returns an array of I<VertexIDs> containing vertices with degree less than I<Degree> in
a I<Graph>.

=item B<HasCycle>

    $Status = $Graph->HasCycle(@VertexIDs);

Returns 1 or 0 based on whether edges between successive pair of I<VertexIDs> including
an additional edge from the last to first vertex ID exists in a I<Graph>.

=item B<HasEdge>

    $Status = $Graph->HasEdge($VertexID1, $VertexID2);

Returns 1 or 0 based on whether an edge between I<VertexID1> and I<VertexID2> exist in
a I<Graph>.

=item B<HasEdgeProperty>

    $Status = $Graph->HasEdgeProperty($PropertyName, $VertexID1,
              $VertexID2);

Returns 1 or 0 based on whether I<PropertyName> has already been associated with an edge
between I<VertexID1> and I<VertexID2> in a I<Graph>.

=item B<HasEdges>

    @EdgesStatus = $Graph->HasEdges(@VertexIDs);
    $FoundEdgesCount = $Graph->HasEdges(@VertexIDs);

Returns an array I<EdgesStatus> containing 1s and 0s corresponding to whether edges between
successive pairs of I<VertexIDs> exist in a I<Graph>. In a scalar context, number of edges found
is returned.

=item B<HasFusedCycles>

    $Status = $Graph->HasFusedCycles();

Returns 1 or 0 based on whether any fused cycles exist in a I<Graph>.

=item B<HasGraphProperty>

    $Status = $Graph->HasGraphProperty($PropertyName);

Returns 1 or 0 based on whether I<PropertyName> has already been associated as a
graph property as opposed to vertex or edge property in a I<Graph>.

=item B<HasPath>

    $Status = $Graph->HasPath(@VertexIDs));

Returns 1 or 0 based on whether edges between all successive pairs of I<VertexIDs> exist in
a I<Graph>.

=item B<HasVertex>

    $Status = $Graph->HasVertex($VertexID);

Returns 1 or 0 based on whether I<VertexID> exists in a I<Graph>.

=item B<HasVertexProperty>

    $Status = $Graph->HasGraphProperty($HasVertexProperty, $VertexID);

Returns 1 or 0 based on whether I<PropertyName> has already been associated with
I<VertexID> in a I<Graph>.

=item B<HasVertices>

    @VerticesStatus = $Graph->HasVertices(@VertexIDs);
    $VerticesFoundCount = $Graph->HasVertices(@VertexIDs);

Returns an array I<> containing 1s and 0s corresponding to whether I<VertexIDs> exist in a 
I<Graph>. In a scalar context, number of vertices found is returned.

=item B<IsAcyclic>

    $Status = $Graph->IsAcyclic();

Returns 0 or 1 based on whether a cycle exist in a I<Graph>.

=item B<IsAcyclicEdge>

    $Status = $Graph->IsAcyclicEdge($VertexID1, $VertexID2);

Returns 0 or 1 based on whether a cycle containing an edge between I<VertexID1> and
I<VertexID2> exists in a I<Graph>.

=item B<IsAcyclicVertex>

    $Status = $Graph->IsAcyclicVertex($VertexID1);

Returns 0 or 1 based on whether a cycle containing a I<VertexID> exists in a I<Graph>.

=item B<IsCyclic>

    $Status = $Graph->IsCyclic();

Returns 1 or 0 based on whether a cycle exist in a I<Graph>.

=item B<IsCyclicEdge>

    $Status = $Graph->IsCyclicEdge($VertexID1, $VertexID2);

Returns 1 or 0 based on whether a cycle containing an edge between I<VertexID1> and
I<VertexID2> exists in a I<Graph>.

=item B<IsCyclicVertex>

    $Status = $Graph->IsCyclicVertex($VertexID1);

Returns 1 or 0 based on whether a cycle containing a I<VertexID> exists in a I<Graph>.

=item B<IsGraph>

    $Status = Graph::IsGraph($Object);

Returns 1 or 0 based on whether I<Object> is a B<Graph> object.

=item B<IsIsolatedVertex>

    $Status = $Graph->IsIsolatedVertex($VertexID);

Returns 1 or 0 based on whether I<VertexID> is an isolated vertex in a I<Graph>. A vertex
with zero as its degree value is considered an isolated vertex.

=item B<IsLeafVertex>

    $Status = $Graph->IsLeafVertex($VertexID);

Returns 1 or 0 based on whether I<VertexID> is an isolated vertex in a I<Graph>. A vertex
with one as its degree value is considered an isolated vertex.

=item B<IsUnicyclic>

    $Status = $Graph->IsUnicyclic();

Returns 1 or 0 based on whether only one cycle is present in a I<Graph>.

=item B<IsUnicyclicEdge>

    $Status = $Graph->IsUnicyclicEdge($VertexID1, $VertexID2);

Returns 1 or 0 based on whether only one cycle contains the edge between I<VertexID1> and
I<VertexID2> in a I<Graph>.

=item B<IsUnicyclicVertex>

    $Status = $Graph->IsUnicyclicVertex($VertexID);

Returns 1 or 0 based on whether only one cycle contains I<VertexID> in a I<Graph>.

=item B<SetActiveCyclicPaths>

    $Graph->SetActiveCyclicPaths($CyclicPathsType);

Sets the type of cyclic paths to use during all methods related to cycles and returns I<Graph>.
Possible values for cyclic paths: I<Independent or All>.

=item B<SetEdgeProperties>

    $Graph->SetEdgeProperties($VertexID1, $VertexID2, @NamesAndValues);

Associates property names and values corresponding to successive pairs of values in
I<NamesAndValues> to an edge between I<VertexID1> and I<VertexID2> in a I<Graph>
and returns I<Graph>.

=item B<SetEdgeProperty>

    $Graph->SetEdgeProperty($Name, $Value, $VertexID1, $VertexID2);

Associates property I<Name> and I<Value> to an edge between I<VertexID1> and I<VertexID2>
in a I<Graph> and returns I<Graph>.

=item B<SetEdgesProperty>

    $Graph->SetEdgesProperty($Name, @ValuesAndVertexIDs);

Associates a same property I<Name> but different I<Values> for different edges specified using
triplets of I<PropertyValue, $VertexID1, $VertexID2> via I<ValuesAndVertexIDs> in a I<graph>.

=item B<SetGraphProperties>

    $Graph->SetGraphProperties(%NamesAndValues);

Associates property names and values I<NamesAndValues> hash to graph as opposed to vertex
or edge and returns I<Graph>.

=item B<SetGraphProperty>

    $Graph->SetGraphProperty($Name, $Value);

Associates property I<Name> and I<Value> to graph as opposed to vertex
or edge and returns I<Graph>.

=item B<SetVertexProperties>

    $Graph->SetVertexProperties($VertexID, @NamesAndValues);

Associates property names and values corresponding to successive pairs of values in
I<NamesAndValues> to I<VertexID> in a I<Graph> and returns I<Graph>.

=item B<SetVertexProperty>

    $Graph->SetVertexProperty($Name, $Value, $VertexID);

Associates property I<Name> and I<Value> to I<VertexID> in a I<Graph> and returns I<Graph>.

=item B<SetVerticesProperty>

    $Graph->SetVerticesProperty($Name, @ValuesAndVertexIDs));

Associates a same property I<Name> but different I<Values> for different vertices specified using
doublets of I<PropertyValue, $VertexID> via I<ValuesAndVertexIDs> in a I<graph>.

=item B<StringifyEdgesProperties>

    $String = $Graph->StringifyEdgesProperties();

Returns a string containing information about properties associated with all edges in a I<Graph> object.

=item B<StringifyGraph>

    $String = $Graph->StringifyGraph();

Returns a string containing information about I<Graph> object.

=item B<StringifyGraphProperties>

    $String = $Graph->StringifyGraphProperties();

Returns a string containing information about properties associated with graph as opposed to vertex.
or an edge in a I<Graph> object

=item B<StringifyProperties>

    $String = $Graph->StringifyProperties();

Returns a string containing information about properties associated with graph, vertices, and edges in
a I<Graph> object.

=item B<StringifyVerticesAndEdges>

    $String = $Graph->StringifyVerticesAndEdges();

Returns a string containing information about vertices and edges in a I<Graph> object.

=item B<StringifyVerticesProperties>

    $String = $Graph->StringifyVerticesProperties();

Returns a string containing information about properties associated with vertices a I<Graph> object.

=item B<UpdateEdgeProperty>

    $Graph->UpdateEdgeProperty($Name, $Value, $VertexID1, $VertexID2);

Updates property I<Value> for I<Name> associated with an edge between I<VertexID1> and
I<VertexID1> and returns I<Graph>.

=item B<UpdateVertexProperty>

    $Graph->UpdateVertexProperty($Name, $Value, $VertexID);

Updates property I<Value> for I<Name> associated with I<VertexID> and returns I<Graph>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

CyclesDetection.pm, Path.pm, PathGraph.pm, PathsTraversal.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
