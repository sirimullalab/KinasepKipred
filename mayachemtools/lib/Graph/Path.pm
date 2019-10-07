package Graph::Path;
#
# File: Path.pm
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

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, $ObjectID);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyPath',

  '==' => '_PathEqualOperator',
  'eq' => '_PathEqualOperator',

  'fallback' => undef;

# Class constructor...
sub new {
  my($Class, @VertexIDs) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializePath();

  if (@VertexIDs) { $This->AddVertices(@VertexIDs); }

  return $This;
}

# Initialize object data...
#
sub _InitializePath {
  my($This) = @_;

  @{$This->{Vertices}} = ();
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Add a vertex to path after the end vertex...
#
sub AddVertex {
  my($This, $VertexID) = @_;

  if (!defined $VertexID ) {
    carp "Warning: ${ClassName}->AddVertex: No vertex added: Vertex ID must be specified...";
    return undef;
  }
  push @{$This->{Vertices}}, $VertexID;

  return $This;
}

# Add vertices to the path after the end vertex...
#
sub AddVertices {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    carp "Warning: ${ClassName}->AddVertices: No vertices added: Vertices list is empty...";
    return undef;
  }
  push @{$This->{Vertices}}, @VertexIDs;

  return $This;
}

# Add a vertex to path after the end vertex...
#
sub PushVertex {
  my($This, $VertexID) = @_;

  return $This->AddVertex($VertexID);
}

# Add vertices to the path after the end vertex...
#
sub PushVertices {
  my($This, @VertexIDs) = @_;

  return $This->AddVertices(@VertexIDs);
}

# Remove end vertex from path...
#
sub PopVertex {
  my($This) = @_;

  if (!@{$This->{Vertices}}) {
    carp "Warning: ${ClassName}->PopVertex: No vertex removed: Path is empty...";
    return undef;
  }
  pop @{$This->{Vertices}};

  return $This;
}

# Remove start vertex from path...
#
sub ShiftVertex {
  my($This) = @_;

  if (!@{$This->{Vertices}}) {
    carp "Warning: ${ClassName}->ShiftVertex: No vertex removed: Path is empty...";
    return undef;
  }
  shift @{$This->{Vertices}};

  return $This;
}

# Add a vertex to path before the start vertex...
#
sub UnshiftVertex {
  my($This, $VertexID) = @_;

  if (!defined $VertexID ) {
    carp "Warning: ${ClassName}->UnshiftVertex: No vertex added: Vertex ID must be specified...";
    return undef;
  }
  unshift @{$This->{Vertices}}, $VertexID;

  return $This;
}

# Add vertices to the path before the start vertex...
#
sub UnshiftVertices {
  my($This, @VertexIDs) = @_;

  if (!@VertexIDs) {
    carp "Warning: ${ClassName}->UnshiftVertices: No vertices added: Vertices list is empty...";
    return undef;
  }
  unshift @{$This->{Vertices}}, @VertexIDs;

  return $This;
}

# Get length...
#
sub GetLength {
  my($This) = @_;

  return scalar @{$This->{Vertices}};
}

# Get start vertex...
#
sub GetStartVertex {
  my($This) = @_;

  if (!$This->GetLength()) {
    return undef;
  }
  my($Index) = 0;
  return $This->_GetVertex($Index);
}

# Get end vertex...
#
sub GetEndVertex {
  my($This) = @_;

  if (!$This->GetLength()) {
    return undef;
  }
  my($Index);

  $Index = $This->GetLength() - 1;
  return $This->_GetVertex($Index);
}

# Get start and end vertices...
#
sub GetTerminalVertices {
  my($This) = @_;

  return ( $This->GetStartVertex(), $This->GetEndVertex() ),
}

# Get path vertices...
#
sub GetVertices {
  my($This) = @_;

  return wantarray ? @{$This->{Vertices}} : scalar @{$This->{Vertices}};
}

# Get a specific vertex from path with indicies starting from 0...
#
sub GetVertex {
  my($This, $Index) = @_;

  if ($Index < 0) {
    croak "Error: ${ClassName}->GetValue: Index value must be a positive number...";
  }
  if ($Index >= $This->GetLength()) {
    croak "Error: ${ClassName}->GetValue: Index vaue must be less than length of path...";
  }
  if (!$This->GetLength()) {
    return undef;
  }
  return $This->_GetVertex($Index);
}

# Get a vertex...
#
sub _GetVertex {
  my($This, $Index) = @_;

  return $This->{Vertices}[$Index];
}

# Get path edges as pair of vertices or number of edges...
#
sub GetEdges {
  my($This) = @_;

  if ($This->GetLength < 1) {
    return undef;
  }
  # Set up edges...
  my($Index, $VertexID1, $VertexID2, @Vertices, @Edges);

  @Edges = ();
  for $Index (0 .. ($#{$This->{Vertices}} - 1) ) {
    $VertexID1 = $This->{Vertices}[$Index];
    $VertexID2 = $This->{Vertices}[$Index + 1];
    push @Edges, ($VertexID1, $VertexID2);
  }

  return wantarray ? @Edges : ((scalar @Edges)/2);
}

# Is it a cycle?
#
sub IsCycle {
  my($This) = @_;
  my($StartVertex, $EndVertex);

  ($StartVertex, $EndVertex) = $This->GetTerminalVertices();

  return ($StartVertex == $EndVertex) ? 1 : 0;
}

# For a path to be an independent path, it must meet the following conditions:
#   . All other vertices are unique.
#
sub IsIndependentPath {
  my($This) = @_;

  # Make sure it has at least two vertices...
  if ($This->GetLength() < 2) {
    return 0;
  }

  # Check frequency of occurence for non-terminal vertices...
  my($VertexID, $IndependenceStatus, @Vertices, %VerticesMap);

  @Vertices = $This->GetVertices();
  shift @Vertices; pop @Vertices;

  %VerticesMap = ();
  $IndependenceStatus = 1;

  VERTEXID: for $VertexID (@Vertices) {
    if (exists $VerticesMap{$VertexID} ) {
      $IndependenceStatus = 0;
      last VERTEXID;
    }
    $VerticesMap{$VertexID} = $VertexID;
  }
  return $IndependenceStatus;
}

# For a path to be an independent cyclic path, it must meet the following conditions:
#   . Termimal vertices are the same
#   . All other vertices are unique.
#
sub IsIndependentCyclicPath {
  my($This) = @_;

  # Make sure it's a cycle...
  if (!($This->GetLength() >= 3 && $This->IsCycle())) {
    return 0;
  }
  return $This->IsIndependentPath();
}

# Is it a path object?
sub IsPath ($) {
  my($Object) = @_;

  return _IsPath($Object);
}

# Copy path...
#
sub Copy {
  my($This) = @_;
  my($NewPath);

  $NewPath = Storable::dclone($This);

  return $NewPath;
}

# Reverse order of vertices in path...
#
sub Reverse {
  my($This) = @_;
  my(@VertexIDs);

  @VertexIDs = (); push @VertexIDs, @{$This->{Vertices}};

  @{$This->{Vertices}} = (); push @{$This->{Vertices}}, reverse @VertexIDs;

  return $This;
}

# Get vertices common between two paths...
#
sub GetCommonVertices {
  my($This, $Other) = @_;
  my($VertexID, @CommonVertices, %OtherVerticesMap);

  # Setup a vertices hash for a quick look up...
  %OtherVerticesMap = ();
  for $VertexID ($Other->GetVertices()) {
    $OtherVerticesMap{$VertexID} = $VertexID;
  }

  @CommonVertices = ();
  for $VertexID ($This->GetVertices()) {
    if ($OtherVerticesMap{$VertexID}) {
      push @CommonVertices, $VertexID
    }
  }
  return wantarray ? @CommonVertices : scalar @CommonVertices;
}

# Join the existing path with a new path specifed using a path object of a list of
# verticies.
#
sub Join {
  my($This, @Values) = @_;

  return $This->_Join(@Values);
}

# Join the existing path with a new path specifed using a path object at a specified
# vertex.
#
sub JoinAtVertex {
  my($This, $Other, $CenterVertexID) = @_;

  # Make sure CenterVertexID is end vertex in This and start vertex in Other before
  # joining them...
  if ($This->GetEndVertex() != $CenterVertexID) {
    $This->Reverse();
  }
  if ($Other->GetStartVertex() != $CenterVertexID) {
    $Other->Reverse();
  }
  return $This->_Join($Other);
}

# Join the existing path with a new path specifed using a path object of a list of
# verticies.
#
# Notes:
#  . Paths must have a common terminal vertex.
#  . Based on the common terminal vertex found, new path vertices are added to the
#    current path in one of the four ways:
#    . New path at end of current path with same vertices order : EndVertex = NewStartVertex
#    . New path at end of current path with reversed vertices order: EndVertex = NewEndVertex
#    . New path at front of current path with same vertices order: StartVertex = NewEndVertex
#    . New path at front of current path with reversed vertices order: StartVertex = NewStartVertex
#
sub _Join {
  my($This, @Values) = @_;

  if (!@Values) {
    return;
  }

  # Get a list of new vertex IDs..
  my($NewPath, $FirstValue, $TypeOfFirstValue, @NewVertexIDs);

  $NewPath = $This->Copy();

  @NewVertexIDs = ();
  $FirstValue = $Values[0];
  $TypeOfFirstValue = ref $FirstValue;
  if ($TypeOfFirstValue =~ /^(SCALAR|HASH|CODE|REF|GLOB)/) {
    croak "Error: ${ClassName}->JoinPath: Trying to add vertices to path object with a reference to unsupported value format...";
  }

  if (_IsPath($FirstValue)) {
    # It's another path object...
    push @NewVertexIDs,  @{$FirstValue->{Vertices}};
  }
  elsif ($TypeOfFirstValue =~ /^ARRAY/) {
    # It's array reference...
    push @NewVertexIDs,  @{$FirstValue};
  }
  else {
    # It's a list of values...
    push @NewVertexIDs,  @Values;
  }
  my($StartVertex, $EndVertex, $NewStartVertex, $NewEndVertex);

  ($StartVertex, $EndVertex) = $NewPath->GetTerminalVertices();
  ($NewStartVertex, $NewEndVertex) = ($NewVertexIDs[0], $NewVertexIDs[$#NewVertexIDs]);

  if (!($EndVertex == $NewStartVertex || $EndVertex == $NewEndVertex || $StartVertex == $NewEndVertex || $StartVertex == $NewStartVertex)) {
    carp "Warning: ${ClassName}->JoinPath: Paths can't be joined: No common terminal vertex found...";
    return undef;
  }

  if ($EndVertex == $NewStartVertex) {
    # Take out EndVertex and add new path at the end...
    pop @{$NewPath->{Vertices}};
    push @{$NewPath->{Vertices}}, @NewVertexIDs;
  }
  elsif ($EndVertex == $NewEndVertex) {
    # Take out EndVertex and add new path at the end with reversed vertex order...
    pop @{$NewPath->{Vertices}};
    push @{$NewPath->{Vertices}}, reverse @NewVertexIDs;
  }
  elsif ($StartVertex == $NewEndVertex) {
    # Take out NewEndVertex and add new path at the front...
    pop @NewVertexIDs;
    unshift @{$NewPath->{Vertices}}, @NewVertexIDs;
  }
  elsif ($StartVertex == $NewStartVertex) {
    # Take out NewStartVertex and add new path at the front...
    shift @NewVertexIDs;
    unshift @{$NewPath->{Vertices}}, reverse @NewVertexIDs;
  }

  return $NewPath,
}

# Compare two paths...
#
sub _PathEqualOperator {
  my($This, $Other) = @_;

  if (!(defined($This) && _IsPath($This) && defined($Other) && _IsPath($Other))) {
    croak "Error: ${ClassName}->_PathEqualOperator: Path equal comparison failed: Both object must be paths...";
  }

  if ($This->GetLength() != $Other->GetLength()) {
    return 0;
  }
  my($ThisID, $OtherID, $ReverseOtherID);

  $ThisID = join('-', @{$This->{Vertices}});
  $OtherID = join('-', @{$Other->{Vertices}});
  $ReverseOtherID = join('-', reverse(@{$Other->{Vertices}}));

  return ($ThisID =~ /^($OtherID|$ReverseOtherID)$/) ? 1 : 0;
}

# Return a string containing vertices in the path...
sub StringifyPath {
  my($This) = @_;
  my($PathString);

  $PathString = "Path: " . join('-', @{$This->{Vertices}});

  return $PathString;
}

# Is it a path object?
sub _IsPath {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

1;

__END__

=head1 NAME

Path - Path class

=head1 SYNOPSIS

use Graph::Path;

use Graph::Path qw(:all);

=head1 DESCRIPTION

B<Path> class provides the following methods:

new, AddVertex, AddVertices, Copy, GetCommonVertices, GetEdges, GetEndVertex,
GetLength, GetStartVertex, GetTerminalVertices, GetVertex, GetVertices, IsCycle,
IsIndependentCyclicPath, IsIndependentPath, IsPath, Join, JoinAtVertex, PopVertex,
PushVertex, PushVertices, Reverse, ShiftVertex, StringifyPath, UnshiftVertex,
UnshiftVertices

Path is a sequential list of vertices with an edge between two successive vertices. The path
becomes a cycle when start vertex and end vertex are the same.

The following operators are overloaded:

    "" == eq

=head2 METHODS

=over 4

=item B<new>

    $NewPath = new Path();
    $NewPath = new Path(@VertexIDs);

Using specified I<VertexIDs>, B<new> method creates a new B<Path> object and returns
newly created B<Path> object.

=item B<AddVertex>

    $Path->AddVertex($VertexID);

Adds I<VertexID> to I<Path> and returns I<Path>.

=item B<AddVertices>

    $Path->AddVertices(@VertexIDs);

Adds vertices using I<VertexIDs> to I<Path> and returns I<Graph>.

=item B<Copy>

    $Return = $Path->Copy();

Copies I<Path> and its associated data using B<Storable::dclone> and returns a new
B<Path> object.

=item B<GetCommonVertices>

    @CommonVertices = $Path->GetCommonVertices($OtherPath);
    $NumOfCommonVertices = $Path->GetCommonVertices($OtherPath);

Returns an array containing common vertex IDs between two paths. In scalar context, number
of common vertices is returned.

=item B<GetEdges>

    @EdgesVertexIDs = $Path->GetEdges();
    $NumOfEdges = $Path->GetEdges();

Returns an array containg successive paris of vertex IDs corresponding to all edges in I<Path>.
In scalar context, the number of edges is returned.

=item B<GetEndVertex>

    $VertexID = $Path->GetEndVertex();

Returns B<VertexID> of end vertex in I<Path>.

=item B<GetLength>

    $Length = $Path->GetLength();

Returns B<Length> of I<Path> corresponding to number of vertices in I<Path>.

=item B<GetStartVertex>

    $VertexID = $Path->GetStartVertex();

Returns B<VertexID> of start vertex in I<Path>.

=item B<GetTerminalVertices>

    ($StartVertexID, $EndVertexID) = $Path->GetTerminalVertices();

Returns vertex IDs of start and end vertices in I<Path>.

=item B<GetVertex>

    $VertexID = $Path->GetVertex($Index);

Returns specific vertex ID from I<Path> corresponding to I<Index> with indicies starting from 0.

=item B<GetVertices>

    @Vertices = $Path->GetVertices();
    $NumOfVertices = $Path->GetVertices();

Returns an array containing all vertex IDs in I<Path>. In scalar context, number of vertices
is returned.

=item B<IsCycle>

    $Status = $Path->IsCycle();

Returns 1 or 0 based on whether I<Path> is a B<CyclicPath> which has the same start and
end vertex IDs.

=item B<IsIndependentCyclicPath>

    $Status = $Path->IsIndependentCyclicPath();

Returns 1 or 0 based on whether I<Path> is an independent B<CyclicPath>. For a I<Path> to be
an independent cyclic path, it must be a cyclic path and have unique vertices.

=item B<IsIndependentPath>

    $Status = $Path->IsIndependentPath();

Returns 1 or 0 based on whether I<Path> is an independent B<Path>. For a I<Path> to be
an independent path, it must have unique vertices.

=item B<IsPath>

    $Status = Graph::Path::IsPath();

Returns 1 or 0 based on whether I<Object> is a B<Path> object

=item B<Join>

    $NewPath = $Path->Join($OtherPath);
    $NewPath = $Path->Join(@VertexIDs);

Joins existing I<Path> with a new path specified as a I<OtherPath> object or an array of I<VertexIDs>
and returns I<NewPath>.

In order to successfully join two paths, terminal vertices must have a common vertex. Based on the
common terminal vertex found, additional path vertices are added to the current I<Path> in one of
the following four ways:

    . EndVertex = NewStartVertex: New path at end of current path with
      same vertices order

    . EndVertex = NewEndVertex: New path at end of current path with
      reversed vertices order

    . StartVertex = NewEndVertex: New path at front of current path
      with same vertices order

    . StartVertex = NewStartVertex: New path at front of current path
      with reversed vertices order

=item B<JoinAtVertex>

    $NewPath = $Path->JoinAtVertex($OtherPath, $CenterVertexID);

Joins existing I<Path> with I<OtherPath> at a specified I<CeterVertexID> and returns a I<NewPath>.

=item B<PopVertex>

    $Path->PopVertex();

Removes end vertex from I<Path> and returns I<Path>.

=item B<PushVertex>

    $Path->PushVertex($VertexID);

Adds I<VertexID> to I<Path> after end vertex and returns I<Path>.

=item B<PushVertices>

    $Path->PushVertices(@VertexIDs);

Adds I<VertexIDs> to I<Path> after end vertex and returns I<Path>.

=item B<Reverse>

    $Path->Reverse();

Reverses order of vertices in I<Path> and returns I<Path>.

=item B<ShiftVertex>

    $Path->ShiftVertex();

Removes start vertex from I<Path> and returns I<Path>.

=item B<StringifyPath>

    $String = $Path->StringifyPath();

Returns a string containing information about I<Path> object.

=item B<UnshiftVertex>

    $Path->UnshiftVertex($VertexID);

Adds I<VertexID> to I<Path> before start vertex and returns I<Path>.

=item B<UnshiftVertices>

    $Path->UnshiftVertices(@VertexIDs);

Adds I<VertexIDs> to I<Path> before start vertex and returns I<Path>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

PathGraph.pm, PathsTraversal.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
