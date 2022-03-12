"""
Functions for determining whether or not a given edge e is isotopic to a
Seifert fibre, where we know in advance that e belongs to a one-vertex
triangulation of a small Seifert fibre space (i.e., a Seifert fibre space
over S^2 with three exceptional fibres).
"""
from enum import Enum, auto
from regina import *
from pachner import twoThree


class Fibre(Enum):
    UNKNOWN = auto()
    NONFIBRE = auto()
    EXCEPTIONAL = auto()
    REGULAR = auto()


def isBlockedSFSOverDisc_(tri):
    """
    Is the given triangulation a blocked Seifert fibre space over the disc?

    If this routine is able to determine that the answer is yes, then it
    returns the number of exceptional fibres. If this routine is able to
    determine that the answer is no, then it returns -1. This routine does
    not guarantee that it will always come to one of these two conclusions,
    and will raise a ValueError whenever it encounters an unsolved case.
    """
    blocked = BlockedSFS.recognise(tri)
    if blocked is None:
        return False
    sfs = blocked.manifold()
    if sfs.name() == "M/n2 x~ S1":
        # Homeomorphic to SFS [D: (2,1) (2,-1)].
        return 2
    if ( sfs.baseOrientable() and sfs.baseGenus() == 0 and
            sfs.punctures() == 1 ):
        # Base is a disc.
        return sfs.fibreCount()
    else:
        # We have a blocked Seifert fibre space whose base is not a disc, but
        # since we do not want to do any more work at this point, we simply
        # mark this case as unsolved.
        raise ValueError(
                "Unsolved case: {}".format( sfs.detail().rstrip() ) )


def isFibre(edge):
    """
    Assuming that the given edge belongs to a one-vertex triangulation of a
    small Seifert fibre space, determines whether this edge is isotopic to a
    Seifert fibre.
    """
    # Pinch the given edge. If the resulting ideal triangulation admits a
    # strict angle structure, then the edge cannot be a Seifert fibre.
    drilled = Triangulation3( edge.triangulation() )
    drilled.pinchEdge( drilled.edge( edge.index() ) )
    if drilled.hasStrictAngleStructure():
        return Fibre.NONFIBRE

    # Now truncate, and see whether we can prove that the edge is isotopic to
    # a Seifert fibre using combinatorial recognition.
    drilled.idealToFinite()
    drilled.intelligentSimplify()
    drilled.intelligentSimplify()
    counts = []
    def recognise( sig, tri ):
        # The following may raise a ValueError in an unsolved case.
        fibreCount = isBlockedSFSOverDisc_(tri)
        if fibreCount in {2,3}:
            counts.append(fibreCount)
            return True
        else:
            return False
    height = 0
    nThreads = 1
    if drilled.retriangulate( height, nThreads, recognise ):
        if len(counts) != 1:
            raise ValueError( "Wrong number of counts: {}".format(counts) )
        fibreCount = counts[0]
        if fibreCount == 2:
            return Fibre.EXCEPTIONAL
        elif fibreCount == 3:
            return Fibre.REGULAR
        else:
            raise ValueError(
                    "Should never get {} exceptional fibres.".format(
                        fibreCount ) )

    # We need to try more powerful techniques: search for essential annuli
    # and essential tori.
    # TODO
    return Fibre.UNKNOWN


# Test code.
if __name__ == "__main__":
    # Generate test triangulations.
    tests = [
            ( "cPcbbbqxh", "SFS [S2: (2,1) (2,1) (2,-1)]" ),
            ( "dLQbccchhrw", "SFS [S2: (2,1) (2,1) (3,-2)]" ) ]
    def genEdges(sig):
        tri = Triangulation3.fromIsoSig(sig)
        # What are the edges currently in the triangulation?
        for i in range( tri.countEdges() ):
            yield tri, i, "e{}".format(i)
        # What new edges can we introduce using 2-3 moves?
        for i in range( tri.countTriangles() ):
            result = twoThree( tri.triangle(i) )
            if result is None:
                continue
            yield result[0], result[1][-1], "f{}".format(i)

    # Perform tests.
    for sig, name in tests:
        print( sig, name )
        for tri, i, name in genEdges(sig):
            msg = "    " + name + ": {}"
            try:
                fibreType = isFibre( tri.edge(i) )
            except ValueError as err:
                print( msg.format( err ) )
            else:
                print( msg.format( fibreType.name ) )
                if fibreType is Fibre.UNKNOWN:
                    # Can we find annuli among the vertex normal surfaces?
                    bounded = Triangulation3(tri)
                    bounded.pinchEdge( bounded.edge(i) )
                    bounded.idealToFinite()
                    bounded.intelligentSimplify()
                    bounded.intelligentSimplify()
                    surfs = NormalSurfaces.enumerate( bounded, NS_STANDARD )
                    for s in surfs:
                        if ( s.isOrientable() and s.eulerChar() == 0 and
                                s.hasRealBoundary() ):
                            print( "        {} ...".format(
                                s.detail()[:60] ) )
        print()

