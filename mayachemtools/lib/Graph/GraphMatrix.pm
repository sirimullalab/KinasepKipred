package Graph::GraphMatrix;
#
# File: GraphMatrix.pm
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
use Matrix;
use Constants;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyGraphMatrix';

# Class constructor...
sub new {
  my($Class, $Graph) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializeGraphMatrix($Graph);

  return $This;
}

# Initialize object data...
sub _InitializeGraphMatrix {
  my($This, $Graph) = @_;

  # Specified graph object...
  $This->{Graph} = $Graph;

  # Generated matrix...
  $This->{Matrix} = undef;
  $This->{MatrixType} = undef;

  # Row and columns IDs...
  @{$This->{RowIDs}} = ();
  @{$This->{ColumnIDs}} = ();

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Generate the adjacency matrix for a simple graph.
#
# For a simple graph G with n vertices, the adjacency matrix for G is a n x n square matrix and
# its elements Mij are:
#
#   . 0    if i == j
#   . 1    if i != j and vertex Vi is adjacent to vertex Vj
#   . 0    if i != j and vertex Vi is not adjacent to vertex Vj
#
sub GenerateAdjacencyMatrix {
  my($This) = @_;
  my($Graph, $NumOfVertices, @VertexIDs);

  # Get graph vertices...
  $Graph = $This->{Graph};
  @VertexIDs = $Graph->GetVertices();
  $NumOfVertices = scalar @VertexIDs;

  if ($NumOfVertices == 0) {
    carp "Warning: ${ClassName}->GenerateAdjacencyMatrix: Specified graph doesn't contain any vertices: No matrix generated...";
    return $This;
  }

  # Create adjacency matrix...
  my($Matrix, $RowIndex, $ColIndex, $RowVertexID, $ColVertexID, $Value, $SkipIndexCheck);

  $Matrix = new Matrix($NumOfVertices, $NumOfVertices);
  $SkipIndexCheck = 1;

  for $RowIndex (0 .. ($NumOfVertices - 1)) {
    $RowVertexID = $VertexIDs[$RowIndex];
    for $ColIndex (0 .. ($NumOfVertices - 1)) {
      $ColVertexID = $VertexIDs[$ColIndex];
      $Value = ($RowIndex == $ColIndex) ? 0 : ($Graph->HasEdge($RowVertexID, $ColVertexID) ? 1 : 0);
      $Matrix->SetValue($RowIndex, $ColIndex, $Value, $SkipIndexCheck);
    }
  }
  $This->_SetMatrixAndAssociatedInformation($Matrix, "AdjacencyMatrix", \@VertexIDs);

  return $This;
}

# Generate the Siedel adjacency matrix for a simple graph.
#
# For a simple graph G with n vertices, the Siedal adjacency matrix for G is a n x n square matrix and
# its elements Mij are:
#
#   . 0    if i == j
#   . -1   if i != j and vertex Vi is adjacent to vertex Vj
#   . 1    if i != j and vertex Vi is not adjacent to vertex Vj
#
sub GenerateSiedelAdjacencyMatrix {
  my($This) = @_;
  my($Graph, $NumOfVertices, @VertexIDs);

  # Get graph vertices...
  $Graph = $This->{Graph};
  @VertexIDs = $Graph->GetVertices();
  $NumOfVertices = scalar @VertexIDs;

  if ($NumOfVertices == 0) {
    carp "Warning: ${ClassName}->GenerateSiedelAdjacencyMatrix: Specified graph doesn't contain any vertices: No matrix generated...";
    return $This;
  }

  # Create Siedel adjacency matrix...
  my($Matrix, $RowIndex, $ColIndex, $RowVertexID, $ColVertexID, $Value, $SkipIndexCheck);

  $Matrix = new Matrix($NumOfVertices, $NumOfVertices);
  $SkipIndexCheck = 1;

  for $RowIndex (0 .. ($NumOfVertices - 1)) {
    $RowVertexID = $VertexIDs[$RowIndex];
    for $ColIndex (0 .. ($NumOfVertices - 1)) {
      $ColVertexID = $VertexIDs[$ColIndex];
      $Value = ($RowIndex == $ColIndex) ? 0 : ($Graph->HasEdge($RowVertexID, $ColVertexID) ? -1 : 1);
      $Matrix->SetValue($RowIndex, $ColIndex, $Value, $SkipIndexCheck);
    }
  }
  $This->_SetMatrixAndAssociatedInformation($Matrix, "SiedelAdjacencyMatrix", \@VertexIDs);

  return $This;
}

# Generate distance matrix for a simple graph using Floyd-Marshall algorithm [Ref 67].
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
sub GenerateDistanceMatrix {
  my($This) = @_;
  my($Graph, $NumOfVertices, @VertexIDs);

  # Get graph vertices...
  $Graph = $This->{Graph};
  @VertexIDs = $Graph->GetVertices();
  $NumOfVertices = scalar @VertexIDs;

  if ($NumOfVertices == 0) {
    carp "Warning: ${ClassName}->GenerateDistanceMatrix: Specified graph doesn't contain any vertices: No matrix generated...";
    return $This;
  }

  # Initialize matrix...
  my($Matrix, $MatrixValuesRef, $RowIndex, $ColIndex, $RowVertexID, $ColVertexID, $Value);

  $Matrix = new Matrix($NumOfVertices, $NumOfVertices);
  $MatrixValuesRef = $Matrix->GetMatrixValuesReference();

  for $RowIndex (0 .. ($NumOfVertices - 1)) {
    $RowVertexID = $VertexIDs[$RowIndex];
    for $ColIndex (0 .. ($NumOfVertices - 1)) {
      $ColVertexID = $VertexIDs[$ColIndex];
      $Value = ($RowIndex == $ColIndex) ? 0 : ($Graph->HasEdge($RowVertexID, $ColVertexID) ? 1 : BigNumber);
      $MatrixValuesRef->[$RowIndex][$ColIndex] = $Value;
    }
  }

  # Create distance matrix...
  my($i, $j, $k, $Valuejk);

  for $i (0 .. ($NumOfVertices - 1)) {
    for $j (0 .. ($NumOfVertices - 1)) {
      for $k (0 .. ($NumOfVertices - 1)) {
	$Valuejk = $MatrixValuesRef->[$j][$i] + $MatrixValuesRef->[$i][$k];
	if ($Valuejk < $MatrixValuesRef->[$j][$k]) {
	  $MatrixValuesRef->[$j][$k] = $Valuejk;
	}
      }
    }
  }
  $This->_SetMatrixAndAssociatedInformation($Matrix, "DistanceMatrix", \@VertexIDs);

  return $This;
}

# Generate the incidence matrix for a simple graph.
#
# For a simple graph G with n vertices and e edges, the incidence matrix for G is a n x e matrix
# its elements Mij are:
#
#   . 1    if vertex Vi and the edge Ej are incident; in other words, Vi and Ej are related
#   . 0    otherwise
#
sub GenerateIncidenceMatrix {
  my($This) = @_;
  my($Graph, $NumOfVertices, $NumOfEdges, @VertexIDs, @EdgeVertexIDs);

  # Get graph vertices and edges...
  $Graph = $This->{Graph};
  @VertexIDs = $Graph->GetVertices();
  $NumOfVertices = scalar @VertexIDs;

  @EdgeVertexIDs = $Graph->GetEdges();
  $NumOfEdges = @EdgeVertexIDs/2;

  if ($NumOfVertices == 0) {
    carp "Warning: ${ClassName}->GenerateIncidenceMatrix: Specified graph doesn't contain any vertices: No matrix generated...";
    return $This;
  }

  # Create incidence matrix...
  my($Matrix, $RowIndex, $ColIndex, $EdgeVertexIndex, $VertexID, $EdgeVertexID1, $EdgeVertexID2, $Value, $SkipIndexCheck);

  $Matrix = new Matrix($NumOfVertices, $NumOfEdges);

  $SkipIndexCheck = 1;
  for $RowIndex (0 .. ($NumOfVertices - 1)) {
    $VertexID = $VertexIDs[$RowIndex];
    $EdgeVertexIndex = 0;
    for $ColIndex (0 .. ($NumOfEdges - 1)) {
      $EdgeVertexID1 = $EdgeVertexIDs[$EdgeVertexIndex]; $EdgeVertexIndex++;
      $EdgeVertexID2 = $EdgeVertexIDs[$EdgeVertexIndex]; $EdgeVertexIndex++;

      $Value = ($VertexID == $EdgeVertexID1 || $VertexID == $EdgeVertexID2) ? 1 : 0;
      $Matrix->SetValue($RowIndex, $ColIndex, $Value, $SkipIndexCheck);
    }
  }
  $This->_SetMatrixAndAssociatedInformation($Matrix, "IncidenceMatrix", \@VertexIDs, \@EdgeVertexIDs);

  return $This;
}

# Generate the degree matrix for a simple graph.
#
# For a simple graph G with n vertices, the degree matrix for G is a n x n square matrix and
# its elements Mij are:
#
#   . deg(Vi)   if i == j and deg(Vi) is the degree of vertex Vi
#   . 0         otherwise
#
sub GenerateDegreeMatrix {
  my($This) = @_;
  my($Graph, $NumOfVertices, @VertexIDs);

  # Get graph vertices...
  $Graph = $This->{Graph};
  @VertexIDs = $Graph->GetVertices();
  $NumOfVertices = scalar @VertexIDs;

  if ($NumOfVertices == 0) {
    carp "Warning: ${ClassName}->GenerateDegreeMatrix: Specified graph doesn't contain any vertices: No matrix generated...";
    return $This;
  }

  # Create degree matrix...
  my($Matrix, $Index, $VertexID, $Degree, $SkipIndexCheck);

  $Matrix = new Matrix($NumOfVertices, $NumOfVertices);
  $SkipIndexCheck = 1;

  for $Index (0 .. ($NumOfVertices - 1)) {
    $VertexID = $VertexIDs[$Index];
    $Degree = $Graph->GetDegree($VertexID);
    $Matrix->SetValue($Index, $Index, $Degree, $SkipIndexCheck);
  }
  $This->_SetMatrixAndAssociatedInformation($Matrix, "DegreeMatrix", \@VertexIDs);

  return $This;
}

# Generate the Laplacian matrix for a simple graph.
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
sub GenerateLaplacianMatrix {
  my($This) = @_;
  my($Graph, $NumOfVertices, @VertexIDs);

  # Get graph vertices...
  $Graph = $This->{Graph};
  @VertexIDs = $Graph->GetVertices();
  $NumOfVertices = scalar @VertexIDs;

  if ($NumOfVertices == 0) {
    carp "Warning: ${ClassName}->GenerateLaplacianMatrix: Specified graph doesn't contain any vertices: No matrix generated...";
    return $This;
  }

  # Create adjacency matrix...
  my($Matrix, $RowIndex, $ColIndex, $RowVertexID, $ColVertexID, $Value, $SkipIndexCheck);

  $Matrix = new Matrix($NumOfVertices, $NumOfVertices);
  $SkipIndexCheck = 1;

  for $RowIndex (0 .. ($NumOfVertices - 1)) {
    $RowVertexID = $VertexIDs[$RowIndex];
    for $ColIndex (0 .. ($NumOfVertices - 1)) {
      $ColVertexID = $VertexIDs[$ColIndex];
      $Value = ($RowIndex == $ColIndex) ? $Graph->GetDegree($RowVertexID) : ($Graph->HasEdge($RowVertexID, $ColVertexID) ? -1 : 0);
      $Matrix->SetValue($RowIndex, $ColIndex, $Value, $SkipIndexCheck);
    }
  }
  $This->_SetMatrixAndAssociatedInformation($Matrix, "LaplacianMatrix", \@VertexIDs);

  return $This;
}

# Generate the normalized Laplacian matrix for a simple graph.
#
# For a simple graph G with n vertices, the normalized Laplacian matrix L for G is a n x n square matrix and
# its elements Lij are:
#
#   . 1                           if i == j and deg(Vi) != 0
#   . -1/SQRT(deg(Vi) * deg(Vj))  if i != j and vertex Vi is adjacent to vertex Vj
#   . 0                           otherwise
#
#
sub GenerateNormalizedLaplacianMatrix {
  my($This) = @_;
  my($Graph, $NumOfVertices, @VertexIDs);

  # Get graph vertices...
  $Graph = $This->{Graph};
  @VertexIDs = $Graph->GetVertices();
  $NumOfVertices = scalar @VertexIDs;

  if ($NumOfVertices == 0) {
    carp "Warning: ${ClassName}->GenerateNormalizedLaplacianMatrix: Specified graph doesn't contain any vertices: No matrix generated...";
    return $This;
  }

  # Create adjacency matrix...
  my($Matrix, $RowIndex, $ColIndex, $RowVertexID, $ColVertexID, $RowVertexDegree, $ColVertexDegree, $Value, $SkipIndexCheck);

  $Matrix = new Matrix($NumOfVertices, $NumOfVertices);
  $SkipIndexCheck = 1;

  for $RowIndex (0 .. ($NumOfVertices - 1)) {
    $RowVertexID = $VertexIDs[$RowIndex];
    $RowVertexDegree = $Graph->GetDegree($RowVertexID);
    for $ColIndex (0 .. ($NumOfVertices - 1)) {
      $ColVertexID = $VertexIDs[$ColIndex];
      $Value = 0;
      if ($RowIndex == $ColIndex) {
	$Value = ($RowVertexDegree == 0) ? 0 : 1;
      }
      else {
	$ColVertexDegree = $Graph->GetDegree($ColVertexID);
	$Value = $Graph->HasEdge($RowVertexID, $ColVertexID) ? (-1/sqrt($RowVertexDegree * $ColVertexDegree)) : 0;
      }
      $Matrix->SetValue($RowIndex, $ColIndex, $Value, $SkipIndexCheck);
    }
  }
  $This->_SetMatrixAndAssociatedInformation($Matrix, "NormalizedLaplacianMatrix", \@VertexIDs);

  return $This;
}

# Generate the admittance matrix for a simple graph.
#
sub GenerateAdmittanceMatrix {
  my($This) = @_;

  $This->GenerateLaplacianMatrix();
  $This->{MatrixType} = "AdmittanceMatrix";

  return $This;
}

# Generate the Kirchhoff matrix for a simple graph.
#
sub GenerateKirchhoffMatrix {
  my($This) = @_;

  $This->GenerateLaplacianMatrix();
  $This->{MatrixType} = "KirchhoffMatrix";

  return $This;
}

# Get generated matrix...
#
sub GetMatrix {
  my($This) = @_;

  return $This->{Matrix};
}

# Get matrix type...
#
sub GetMatrixType {
  my($This) = @_;

  return $This->{MatrixType};
}

# Get row IDs...
#
sub GetRowIDs {
  my($This) = @_;

  return @{$This->{RowIDs}};
}

# Get column IDs...
#
sub GetColumnIDs {
  my($This) = @_;

  return @{$This->{ColumnIDs}};
}


# Setup matrix and other associated information...
#
sub _SetMatrixAndAssociatedInformation {
  my($This, $Matrix, $MatrixType, $VertexIDsRef, $EdgeVertexIDsRef) = @_;

  $This->{Matrix} = $Matrix;
  $This->{MatrixType} = $MatrixType;

  @{$This->{RowIDs}} = ();
  @{$This->{ColumnIDs}} = ();

  @{$This->{RowIDs}} = @{$VertexIDsRef};

  if ($MatrixType =~ /^IncidenceMatrix$/i) {
    # Setup column IDs using edge vertex IDs...
    my($NumOfEdges, $EdgeVertexIndex, $ColIndex, $EdgeVertexID1, $EdgeVertexID2, $EdgeID);

    $NumOfEdges = (scalar @{$EdgeVertexIDsRef})/2;
    $EdgeVertexIndex = 0;

    for $ColIndex (0 .. ($NumOfEdges - 1)) {
      $EdgeVertexID1 = $EdgeVertexIDsRef->[$EdgeVertexIndex]; $EdgeVertexIndex++;
      $EdgeVertexID2 = $EdgeVertexIDsRef->[$EdgeVertexIndex]; $EdgeVertexIndex++;
      $EdgeID = ($EdgeVertexID1 < $EdgeVertexID2) ? "${EdgeVertexID1}-${EdgeVertexID2}" : "${EdgeVertexID2}-${EdgeVertexID1}";

      push @{$This->{ColumnIDs}}, $EdgeID;
    }
  }
  else {
    push @{$This->{ColumnIDs}}, @{$VertexIDsRef};
  }
  return $This;
}

# Return a string containg data for GraphMatrix object...
sub StringifyGraphMatrix {
  my($This) = @_;
  my($GraphMatrixString, $RowIDs);

  $GraphMatrixString = "GraphMatrix:";

  $GraphMatrixString .= (defined $This->{MatrixType}) ? " MatrixType: $This->{MatrixType}" : "MatrixType: <Undefined>";
  $GraphMatrixString .= (scalar @{$This->{RowIDs}}) ? "; RowIDs: @{$This->{RowIDs}}" : "; RowIDs: <Undefined>";
  $GraphMatrixString .= (scalar @{$This->{ColumnIDs}}) ? "; ColumnIDs: @{$This->{ColumnIDs}}" : "; ColumnIDs: <Undefined>";

  $GraphMatrixString .= (defined $This->{Matrix}) ? "; Matrix: $This->{Matrix}" : "; Matrix: <Undefined>";

  return $GraphMatrixString;
}

1;

__END__

=head1 NAME

GraphMatrix

=head1 SYNOPSIS

use Graph::GraphMatrix;

use Graph::GraphMatrix qw(:all);

=head1 DESCRIPTION

B<GraphMatrix> class provides the following methods:

new, GenerateAdjacencyMatrix, GenerateAdmittanceMatrix, GenerateDegreeMatrix,
GenerateDistanceMatrix, GenerateIncidenceMatrix, GenerateKirchhoffMatrix,
GenerateLaplacianMatrix, GenerateNormalizedLaplacianMatrix,
GenerateSiedelAdjacencyMatrix, GetColumnIDs, GetMatrix, GetMatrixType, GetRowIDs,
StringifyGraphMatrix

=head2 METHODS

=over 4

=item B<new>

    $NewGraphMatrix = new Graph::GraphMatrix($Graph);

Using specified I<Graph>, B<new> method creates a new B<GraphMatrix> and returns
newly created B<GraphMatrix>.

=item B<GenerateAdjacencyMatrix>

    $AdjacencyGraphMatrix = $GraphMatrix->GenerateAdjacencyMatrix();

Generates a new I<AdjacencyGraphMatrix> for specified B<Graph> and returns
I<AdjacencyGraphMatrix>.

For a simple graph G with n vertices, the adjacency matrix for G is a n x n square matrix and
its elements Mij are:

    . 0    if i == j
    . 1    if i != j and vertex Vi is adjacent to vertex Vj
    . 0    if i != j and vertex Vi is not adjacent to vertex Vj

=item B<GenerateAdmittanceMatrix>

    $AdmittanceGraphMatrix = $GraphMatrix->GenerateAdmittanceMatrix();

Generates a new I<AdmittanceGraphMatrix> for specified B<Graph> and returns
I<AdmittanceGraphMatrix>.

B<AdmittanceMatrix> is another name for B<LaplacianMatrix>.

=item B<GenerateDegreeMatrix>

    $DegreeGraphMatrix = $GraphMatrix->GenerateDegreeMatrix();

Generates a new I<DegreeGraphMatrix> for specified B<Graph> and returns
I<DegreeGraphMatrix>.

For a simple graph G with n vertices, the degree matrix for G is a n x n square matrix and
its elements Mij are:

    . deg(Vi)   if i == j and deg(Vi) is the degree of vertex Vi
    . 0         otherwise

=item B<GenerateDistanceMatrix>

    $DistanceGraphMatrix = $GraphMatrix->GenerateDistanceMatrix();

Generates a new I<DistanceGraphMatrix> for specified B<Graph> using Floyd-Marshall
algorithm [Ref 67] and returns I<DistanceGraphMatrix>.

For a simple graph G with n vertices, the distance matrix for G is a n x n square matrix and
its elements Mij are:

    . 0    if i == j
    . d    if i != j and d is the shortest distance between vertex Vi and vertex Vj

In the final matrix, value of constant B<BigNumber> defined in B<Constants.pm> module
corresponds to vertices with no edges.

=item B<GenerateIncidenceMatrix>

    $IncidenceGraphMatrix = $GraphMatrix->GenerateIncidenceMatrix();

Generates a new I<IncidenceGraphMatrix> for specified B<Graph> and returns
I<IncidenceGraphMatrix>.

For a simple graph G with n vertices and e edges, the incidence matrix for G is a n x e matrix
its elements Mij are:

    . 1    if vertex Vi and the edge Ej are incident; in other words, Vi and Ej are related
    . 0    otherwise

=item B<GenerateKirchhoffMatrix>

    $KirchhoffGraphMatrix = $GraphMatrix->GenerateKirchhoffMatrix();

Generates a new I<KirchhoffGraphMatrix> for specified B<Graph> and returns
I<KirchhoffGraphMatrix>.

B<KirchhoffMatrix> is another name for B<LaplacianMatrix>.

=item B<GenerateLaplacianMatrix>

    $LaplacianGraphMatrix = $GraphMatrix->GenerateLaplacianMatrix();

Generates a new I<LaplacianGraphMatrix> for specified B<Graph> and returns
I<LaplacianGraphMatrix>.

For a simple graph G with n vertices, the Laplacian matrix for G is a n x n square matrix and
its elements Mij are:

    . deg(Vi)   if i == j and deg(Vi) is the degree of vertex Vi
    . -1        if i != j and vertex Vi is adjacent to vertex Vj
    . 0         otherwise

The Laplacian matrix is the difference between the degree matrix and adjacency matrix.

=item B<GenerateNormalizedLaplacianMatrix>

    $NormalizedLaplacianGraphMatrix = $GraphMatrix->GenerateNormalizedLaplacianMatrix();

Generates a new I<NormalizedLaplacianGraphMatrix> for specified B<Graph> and returns
I<NormalizedLaplacianGraphMatrix>.

For a simple graph G with n vertices, the normalized Laplacian matrix L for G is a n x n square
matrix and its elements Lij are:

    .  1                           if i == j and deg(Vi) != 0
    .  -1/SQRT(deg(Vi) * deg(Vj))  if i != j and vertex Vi is adjacent to vertex Vj
    .  0                           otherwise

=item B<GenerateSiedelAdjacencyMatrix>

    $SiedelAdjacencyGraphMatrix = $GraphMatrix->GenerateSiedelAdjacencyMatrix();

Generates a new I<SiedelAdjacencyGraphMatrix> for specified B<Graph> and returns
I<SiedelAdjacencyGraphMatrix>.

For a simple graph G with n vertices, the Siedal adjacency matrix for G is a n x n square matrix and
its elements Mij are:

    . 0    if i == j
    . -1   if i != j and vertex Vi is adjacent to vertex Vj
    . 1    if i != j and vertex Vi is not adjacent to vertex Vj

=item B<GetColumnIDs>

    @ColumnIDs = $GraphMatrix->GetColumnIDs();

Returns an array containing any specified column IDs for I<GraphMatrix>.

=item B<GetMatrix>

    $Matrix = $GraphMatrix->GetMatrix();

Returns I<Matrix> object corresponding to I<GraphMatrix> object.

=item B<GetMatrixType>

    $MatrixType = $GraphMatrix->GetMatrixType();

Returns B<MatrixType> of I<GraphMatrix>.

=item B<GetRowIDs>

    @RowIDs = $GraphMatrix->GetRowIDs();

Returns an array containing any specified rowIDs IDs for I<GraphMatrix>.

=item B<StringifyGraphMatrix>

    $String = $GraphMatrix->StringifyGraphMatrix();

Returns a string containing information about I<GraphMatrix> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Constants.pm, Graph.pm, Matrix.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
