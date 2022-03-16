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
        try:
            # May raise a ValueError in an unsolved case.
            fibreCount = isBlockedSFSOverDisc_(tri)
        except ValueError:
            return False
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

    # TODO Use compressing discs to rule out "degenerate" cases.
    # We need to try more powerful techniques: search for essential annuli.
    # Fast tests have failed us, so we turn to something conclusive: search
    # for essential annuli.
    # If we started with an edge isotopic to a fibre in a small Seifert fibre
    # space, then after drilling we should have a Seifert fibre space over
    # the disc with either two or three exceptional fibres.
    #   --> If there are two exceptional fibres, then there will be an
    #       (essential) annulus that cuts the space into two solid tori.
    #   --> If there are three exceptional fibres, then there will be an
    #       (essential) annulus that cuts the space into a solid torus and
    #       a Seifert fibre space over the disc with two exceptional fibres.
    # Since the annuli of interest are essential, such annuli will appear as
    # normal surfaces whenever they exist.
    surfs = NormalSurfaces.enumerate( drilled, NS_STANDARD )
    annuli = []
    for s in surfs:
        if not s.isOrientable() or not s.hasRealBoundary():
            continue
        euler = s.eulerChar()
        if euler == 1:
            # The surface s is a disc. Check that it isn't a compressing disc
            # by cutting along s and verifying that the boundary is still
            # just a torus.
            compress = s.cutAlong()
            if ( compress.countBoundaryComponents() != 1 or
                    compress.boundaryComponent(0).eulerChar() != 0 ):
                return Fibre.NONFIBRE
        elif euler == 0:
            annuli.append(s)
    if not annuli:
        # The annulus that we are looking for doesn't exist.
        return Fibre.NONFIBRE
    for a in annuli:
        # Try cutting along the annulus. If this is the annulus we are
        # looking for, then we must get two pieces, at least one of which is
        # a solid torus.
        cut = a.cutAlong()
        if cut.countComponents() != 2:
            # Cutting along a only yields one piece, so it isn't the annulus
            # we're looking for. Move along to the next annulus.
            continue
        notSolidTorus = []
        for comp in cut.triangulateComponents():
            comp.intelligentSimplify()
            comp.intelligentSimplify()
            if not comp.isSolidTorus():
                notSolidTorus.append(comp)
        nonSolidTorusCount = len(notSolidTorus)
        if nonSolidTorusCount == 0:
            return Fibre.EXCEPTIONAL
        elif nonSolidTorusCount == 2:
            # Cutting along a doesn't yield any solid torus pieces, so it
            # isn't the annulus we're looking for. Move along.
            continue

        # Getting to this point means that cutting along a yields two pieces,
        # exactly one of which is a solid torus. If a is the annulus we're
        # looking for, then the other piece should be a Seifert fibre space
        # over the disc with two exceptional fibres. First try to verify this
        # using combinatorial recognition.
        other = notSolidTorus[0]
        def recognise2( sig, tri ):
            try:
                # May raise a ValueError in an unsolved case.
                fibreCount = isBlockedSFSOverDisc_(tri)
            except ValueError:
                return False
            return ( fibreCount == 2 )
        height = 0
        nThreads = 1
        if other.retriangulate( height, nThreads, recognise2 ):
            return Fibre.REGULAR

        # We need to resort to searching for an annulus again. This time, we
        # are looking specifically for an annulus that cuts the other piece
        # into two solid tori.
        cutSurfs = NormalSurfaces.enumerate( notSolidTorus[0], NS_STANDARD )
        cutAnnuli = []
        hasCompress = False
        for s in cutSurfs:
            if not s.isOrientable() or not s.hasRealBoundary():
                continue
            euler = s.eulerChar()
            if euler == 1:
                # The surface s is a disc. Check that it isn't a compressing
                # disc by cutting along s and verifying that the boundary is
                # still just a torus.
                compress = s.cutAlong()
                if ( compress.countBoundaryComponents() != 1 or
                        compress.boundaryComponent(0).eulerChar != 0 ):
                    hasCompress = True
                    break
            elif euler == 0:
                cutAnnuli.append(s)
        if hasCompress or not cutAnnuli:
            # We have a compressing disc, or we don't have any annuli. In
            # either case, the annulus a isn't the one we're looking for.
            # Move along.
            continue
        for aa in cutAnnuli:
            # Try cutting along aa.
            c = aa.cutAlong()
            if c.countComponents() != 2:
                # Cutting along aa only yields one piece, so it isn't the
                # annulus we're looking for. Move along.
                continue
            onlySolidTori = True
            for cc in c.triangulateComponents():
                cc.intelligentSimplify()
                cc.intelligentSimplify()
                if not cc.isSolidTorus():
                    # Cutting along aa doesn't yield two solid tori, so it
                    # isn't the annulus we're looking for. Move along.
                    onlySolidTori = False
                    break
            if onlySolidTori:
                return Fibre.REGULAR

    # Surviving to this point means that the annulus we were looking for
    # doesn't exist.
    return Fibre.NONFIBRE


# Test code.
if __name__ == "__main__":
    # Generate test triangulations.
    tests = [
            ( "cPcbbbqxh", "SFS [S2: (2,1) (2,1) (2,-1)]" ),
            ( "dLQbccchhrw", "SFS [S2: (2,1) (2,1) (3,-2)]" ),
            ( "eLAkbccddemken", "SFS [S2: (2,1) (2,1) (2,1)]" ),
            ( "eLAkbccddemkij", "SFS [S2: (2,1) (2,1) (3,-1)] : #1" ),
            ( "eLPkbcddddrwos", "SFS [S2: (2,1) (2,1) (3,-1)] : #2" ),
            ( "eLMkbcdddhhhqx", "SFS [S2: (2,1) (2,1) (4,-3)]" ),
            ( "eLPkbcdddhrrnk", "SFS [S2: (2,1) (3,1) (3,-2)]" ) ]
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
        for i in range(10):
            RandomEngine.reseedWithHardware()
            ii = RandomEngine.rand( tri.countTriangles() )
            result = twoThree( tri.triangle(ii) )
            if result is None:
                break
            tri = result[0]
            yield tri, result[1][-1], "r{}".format(i)

    # Perform tests.
    for sig, name in tests:
        print()
        print( sig, name )
        for tri, i, name in genEdges(sig):
            msg = "    " + name + ": {}"
            try:
                fibreType = isFibre( tri.edge(i) )
            except ValueError as err:
                print( msg.format( err ) )
            else:
                if fibreType is Fibre.NONFIBRE:
                    print( msg.format( fibreType.name + "   <--" ) )
                else:
                    print( msg.format( fibreType.name ) )

