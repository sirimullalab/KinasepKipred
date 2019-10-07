package Graph::PathGraph;
#
# File: PathGraph.pm
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
use Graph;
use Graph::Path;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Graph Exporter);
@EXPORT = qw(IsPathGraph);
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, $PathsPropertyName, $CyclicPathsPropertyName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyPathGraph';

# Class constructor...
sub new {
  my($Class, $Graph) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializePathGraph($Graph);

  $This->_ConvertGraphIntoPathGraph($Graph);

  return $This;
}

# Initialize object data...
sub _InitializePathGraph {
  my($This, $Graph) = @_;

  if (!(defined($Graph) && Graph::IsGraph($Graph))) {
    croak "Error: ${ClassName}->new: PathGraph object can't be instantiated without a Graph object...";
  }

  $This->{Graph} = $Graph;

  # Maximum time allowed for cycles detection during collapse vertex cycles detection
  # methodology in seconds...
  $This->{MaxAllowedTime} = 30;

  # Starting time for cycles detection during collapse vertex cycles detection
  # methodology...
  $This->{StartTime} = time;

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Path edge property name...
  $PathsPropertyName = 'Paths';

  # Cyclic path vertex property name...
  $CyclicPathsPropertyName = 'CyclicPaths';
}

# Convert graph into a path graph...
#
sub _ConvertGraphIntoPathGraph {
  my($This, $Graph) = @_;

  # Copy graph vertices and edges without any associated properties data
  # from Graph to This: Graph properties data is available using Graph object reference
  # store in This object...
  #
  $Graph->CopyVerticesAndEdges($This);

  # . Attach Path property to each edge...
  #
  my($Index, $VertexID1, $VertexID2, $Path, @EdgesVertexIDs);

  @EdgesVertexIDs = ();
  @EdgesVertexIDs = $This->GetEdges();
  for ($Index = 0; $Index < $#EdgesVertexIDs; $Index += 2) {
    $VertexID1 = $EdgesVertexIDs[$Index]; $VertexID2 = $EdgesVertexIDs[$Index + 1];
    $Path = new Graph::Path($VertexID1, $VertexID2);
    my(@Paths) = ();
    push @Paths, $Path;
    $This->SetEdgeProperty($PathsPropertyName, \@Paths, $VertexID1, $VertexID2);
  }
  return $This;
}

# Collapse paths around a specified vertex by updating paths around the vertex
# and adding any resulting cyclic paths to vertices attached to specified vertex.
#
# Notes:
#   . Path object references are stored as a list attached to Paths property on edges.
#     Usage of list allows multiple paths attached to the egde between a pair of vertices;
#     Graph doesn't support multiple egdes between a pair of vertices.
#
#   . Cyclic path object references are stored as list on vertices as CyclicPaths graph property.
#     List allows multiple Loop properties attached to a vertex.
#
#   . For topologically complex graphs containing large number of cycles, cycles detection algorithm
#     [ Ref 31 ] as implemented implemented in CollapseVertexAndCollectCyclicPathsDetectCycles
#     might not be able to find all the cycles in a reasonable amount of time and is designed to
#     abandon cycles detection after MaxAllowedTime. Consequently, no cycles are detected
#     or assigned.
#
sub CollapseVertexAndCollectCyclicPaths {
  my($This, $VertexID) = @_;

  if (!$This->HasVertex($VertexID)) {
    carp "Warning: ${ClassName}->CollapseVertexAndCollectCyclicPaths: Didn't collapse vertex $VertexID: Vertex $VertexID doesn't exist...";
    return undef;
  }
  # Collect all paths around specified VertexID by going over paths associated with its edges...
  my($Index, $EdgePathsRef, $EdgeVertexID1, $EdgeVertexID2, @Paths, @EdgesVertexIDs);

  @EdgesVertexIDs = ();
  @EdgesVertexIDs = $This->GetEdges($VertexID);

  @Paths = ();
  for ($Index = 0; $Index < $#EdgesVertexIDs; $Index += 2) {
    ($EdgeVertexID1, $EdgeVertexID2) = ($EdgesVertexIDs[$Index], $EdgesVertexIDs[$Index + 1]);
    $EdgePathsRef = $This->GetEdgeProperty($PathsPropertyName, $EdgeVertexID1, $EdgeVertexID2);
    push @Paths, @{$EdgePathsRef};
  }

  # Go over each pair of paths around the specified vertex, join paths and associate
  # joined path to appropriate edge...
  my($Index1, $Index2, $Path1, $Path2, $JoinedPath, $JoinedPathStartVertexID, $JoinedPathEndVertexID, @CommonVertices);

  for ($Index1 = 0; $Index1 < $#Paths; $Index1 +=1 ) {
    $Path1 = $Paths[$Index1];

    PATH2: for ($Index2 = $Index1 + 1; $Index2 <= $#Paths; $Index2 +=1 ) {
      $Path2 = $Paths[$Index2];

      # For JoinedPath to be valid cycle, Path1 and Path2 must have exactly two vertices in common.
      # Otherwise, joined path contains duplicate vertices besides the terminal vertices and
      # indicates a path from a different direction.
      #
      # For paths leading to cycles, it only makes sense to join paths with only one common vertex;
      # otherwise, it wouldn't lead to a cycle and can be ignored.
      #
      @CommonVertices = $Path1->GetCommonVertices($Path2);
      if (!(@CommonVertices <= 2 && ($CommonVertices[0] == $VertexID || $CommonVertices[1] == $VertexID))) {
	next PATH2;
      }

      $JoinedPath = $Path1->JoinAtVertex($Path2, $VertexID);
      ($JoinedPathStartVertexID, $JoinedPathEndVertexID) = $JoinedPath->GetTerminalVertices();

      if (!$JoinedPath->IsIndependentPath()) {
	next PATH2;
      }

      # Decide whether to give up or keep going...
      if ($This->_IsTimeToGiveUpCyclesDetection()) {
	warn "Warning: ${ClassName}->CollapseVertexAndCollectCyclicPaths: Cycles detection algorithm [ Ref 31 ] as implemented in the current release of MayaChemTools didn't finish with in the maximum allowed time of $This->{MaxAllowedTime} seconds; Cycles detection has been abandoned...";
	return undef;
      }

      if ($JoinedPathStartVertexID == $JoinedPathEndVertexID) {
	# It's a cycle. Attach it to the graph as CylicPaths property...
	if ($This->HasGraphProperty($CyclicPathsPropertyName)) {
	  my($ExistingCyclicPathsRef);
	  $ExistingCyclicPathsRef = $This->GetGraphProperty($CyclicPathsPropertyName);
	  push @{$ExistingCyclicPathsRef}, $JoinedPath;
	}
	else {
	  my(@NewCyclicPaths) = ();
	  push @NewCyclicPaths, $JoinedPath;
	  $This->SetGraphProperty($CyclicPathsPropertyName, \@NewCyclicPaths, $JoinedPathStartVertexID);
	}
      }
      else {
	if ($This->HasEdge($JoinedPathStartVertexID, $JoinedPathEndVertexID)) {
	  # Append to the list of exisiting paths property of the edge...
	  my($ExistingPathsRef);
	  $ExistingPathsRef = $This->GetEdgeProperty($PathsPropertyName, $JoinedPathStartVertexID, $JoinedPathEndVertexID);
	  push @{$ExistingPathsRef}, $JoinedPath;
	}
	else {
	  # Create a new edge and associate path property...
	  my(@NewPaths) = ();
	  push @NewPaths, $JoinedPath;
	  $This->AddEdge($JoinedPathStartVertexID, $JoinedPathEndVertexID);
	  $This->SetEdgeProperty($PathsPropertyName, \@NewPaths, $JoinedPathStartVertexID, $JoinedPathEndVertexID);
	}
      }
    }
  }
  $This->DeleteVertex($VertexID);

  return $This;
}

# Decide whether to give up cycles detection using collapse vertex methodology...
#
sub _IsTimeToGiveUpCyclesDetection {
  my($This) = @_;

  return ((time - $This->{StartTime}) > $This->{MaxAllowedTime}) ? 1 : 0;
}

# Delete vertices with degree less than a specifed degree...
#
sub DeleteVerticesWithDegreeLessThan {
  my($This, $Degree) = @_;
  my($VertexID, @VertexIDs);

  while (@VertexIDs = $This->GetVerticesWithDegreeLessThan($Degree)) {
    for $VertexID (@VertexIDs) {
      $This->DeleteVertex($VertexID);
    }
  }
  return $This;
}

# Get paths associated with edges...
#
sub GetPaths {
  my($This) = @_;
  my($PathsRef, @Paths, @PathsList);

  @Paths = (); @PathsList = ();
  @PathsList = $This->GetEdgesProperty($PathsPropertyName);
  for $PathsRef (@PathsList) {
    push @Paths, @{$PathsRef};
  }
  return wantarray ? @Paths : scalar @Paths;
}

# Get paths associated with edges which make a cylce...
#
sub GetCyclicPaths {
  my($This) = @_;
  my($PathsRef, @Paths, @PathsList);

  @Paths = (); @PathsList = ();
  @PathsList = $This->GetGraphProperty($CyclicPathsPropertyName);
  PATHS: for $PathsRef (@PathsList) {
    if (!(defined($PathsRef) && @{$PathsRef})) {
      next PATHS;
    }
    push @Paths, @{$PathsRef};
  }
  return wantarray ? @Paths : scalar @Paths;
}

# Is it a path graph object?
sub IsPathGraph ($) {
  my($Object) = @_;

  return _IsPathGraph($Object);
}

# Return a string containg data for PathGraph object...
sub StringifyPathGraph {
  my($This) = @_;
  my($PathGraphString);

  $PathGraphString = 'PathGraph:' . $This->StringifyVerticesAndEdges() . '; ' . $This->StringifyProperties();

  return $PathGraphString;
}

# Is it a PathGraph object?
sub _IsPathGraph {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

PathGraph

=head1 SYNOPSIS

use Graph::PathGraph;

use Graph::PathGraph qw(:all);

=head1 DESCRIPTION

B<PathGraph> class provides the following methods:

new, CollapseVertexAndCollectCyclicPaths, DeleteVerticesWithDegreeLessThan,
GetCyclicPaths, GetPaths, IsPathGraph, StringifyPathGraph

B<PathGraph> class is derived from I<Graph> class.

=head2 METHODS

=over 4

=item B<new>

    $NewPathGraph = new Graph::PathGraph($Graph);

Using specified I<Graph>, B<new> method creates a new B<PathGraph> object and returns
newly created B<PathGraph> object.

I<Graph> is converted into a B<PathGraph> by copying all its vertices and edges without any
associated properties data and associating a I<Path> object to each edge containing edge
vertex IDs as intial path.

=item B<CollapseVertexAndCollectCyclicPaths>

    $PathGraph->CollapseVertexAndCollectCyclicPaths($VertexID);

Collapses paths around a I<VertexID> by updating paths around the vertex [Ref 31] and associating any
resulting cyclic paths to graph as B<CyclicPaths> property name. And returns I<PathGraph>.

=item B<DeleteVerticesWithDegreeLessThan>

    $Return = $PathGraph->DeleteVerticesWithDegreeLessThan($Degree);

Deletes vertices with degree less than I<Degree> from I<PathGraph> and returns I<PathGraph>.

=item B<GetCyclicPaths>

    @CyclicPaths = $PathGraph->GetCyclicPaths();
    $NumOfPaths = $PathGraph->GetCyclicPaths();

Returns an array of cyclic I<Paths> associated with edges in I<PathGraph>. In scalar context, number
of cyclic paths is returned.

=item B<GetPaths>

    @Paths = $PathGraph->GetPaths();
    $NumOfPaths = $PathGraph->GetPaths();

Returns an array of I<Paths> associated with edges in I<PathGraph>. In scalar context, number
of paths is returned.

=item B<IsPathGraph>

    $Status = Graph::PathGraph::IsPathGraph($Object);

Returns 1 or 0 based on whether I<Object> is a B<PathGraph> object.

=item B<StringifyPathGraph>

    $String = $PathGraph->StringifyPathGraph();

Returns a string containing information about traversed paths in I<PathGraph> object.

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
