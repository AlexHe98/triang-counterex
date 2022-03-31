"""
Helper functions for searching the Pachner graph for counterexamples.
"""
from regina import *
from seifert import Fibre, isFibre


def isCore( edge ):
    """
    Returns True if and only if pinching the given edge yields an ideal
    triangulation of the solid torus; if edge belongs to a triangulation of
    a lens space, then this corresponds to edge forming a core loop.

    Note that this routine relies on normal surface theory, so it may be very
    slow if the given edge belongs to a large triangulation.
    """
    tri = Triangulation3( edge.triangulation() )
    tri.pinchEdge( tri.edge( edge.index() ) )
    return tri.isSolidTorus()


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

    Note that this routine relies on normal surface theory, so it may be very
    slow if the given edge belongs to a large triangulation.
    """
    return pinchGivesHandlebody( edge, 2 )


def possibleFibre( edge ):
    """
    Returns True if we cannot rule out the possibility that the given edge is
    isotopic to a fibre in a small Seifert fibre space.

    Note that this routine relies on normal surface theory, so it may be very
    slow if the given edge belongs to a large triangulation.
    """
    RandomEngine.reseedWithHardware()
    fibreType = isFibre(edge)
    if fibreType is Fibre.UNKNOWN:
        RandomEngine.reseedWithHardware()
        fibreType = isFibre(edge)
    return ( fibreType is not Fibre.NONFIBRE )


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

