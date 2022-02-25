"""
Helper functions for searching the Pachner graph for counterexamples.
"""
from regina import *


def pinchGivesHandlebody( edge, genus ):
    """
    Returns True if and only if pinching the given edge yields an ideal
    triangulation of the orientable handlebody with the given genus.

    Note that this routine relies on normal surface theory, so it may be very
    slow if the given edge belongs to a large triangulation.
    """
    tri = Triangulation3( edge.triangulation() )
    tri.pinchEdge( tri.edge( edge.index() ) )
    return ( tri.isHandlebody() == genus )


def isTunnel( edge ):
    """
    Returns True if and only if pinching the given edge yields an ideal
    triangulation of the orientable genus-2 handlebody; if edge belongs to a
    triangulation of the complement of a knot with tunnel number 1, then this
    corresponds to edge forming a tunnel.
    """
    return pinchGivesHandlebody( edge, 2 )


def mulDefect(e):
    """
    Returns the multiplicity defect of the given edge e.

    The multiplicity defect of e is defined to be:
        abs( e.degree() - len(tets) ), where
        tets = { emb.simplex.index() for emb in e.embeddings() }.
    """
    tets = { emb.simplex().index() for emb in e.embeddings() }
    return abs( e.degree() - len(tets) )


def degDefect(e):
    """
    Returns the degree defect of the given edge e.

    The degree defect of e is defined to be:
        abs( e.degree() - 3 ).
    """
    return abs( e.degree() - 3 )


def relabelEdges( tri, indices, isom ):
    """
    Given a list of edge indices, returns the list of new indices obtained by
    applying the given isomorphism to the given triangulation.
    """
    newTri = isom.apply(tri)
    newIndices = []
    for i in indices:
        emb = tri.edge(i).embedding(0)
        simp = emb.simplex().index()
        vert = [ emb.vertices()[j] for j in range(2) ]

        # What happens after applying isom?
        newSimp = isom.simpImage(simp)
        newVert = [ isom.facetPerm(simp)[j] for j in vert ]
        newIndices.append(
                newTri.tetrahedron(newSimp).edge(*newVert).index() )
    return newIndices

