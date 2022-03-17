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
    UNKNOWN = auto()        # Curve of entirely unknown type.
    NONFIBRE = auto()       # Curve that is not a fibre.
    FIBRE = auto()          # Curve that is a fibre of unknown type.
    EXCEPTIONAL = auto()    # Curve that is an exceptional fibre.
    REGULAR = auto()        # Curve that is a regular fibre.


class SFS(Enum):
    UNKNOWN = auto()    # Unknown whether SFS or not.
    NOTSFS = auto()     # Not an SFS.
    NONDISC = auto()    # SFS, but not over disc.
    DISCTWO = auto()    # SFS over disc; 2 exceptional fibres.
    DISCTHREE = auto()  # SFS over disc; 3 exceptional fibres.
    DISCOTHER = auto()  # SFS over disc; neither 2 nor 3 exceptional fibres.


def recogniseSFSOverDisc(tri):
    """
    Is the given triangulation a Seifert fibre space over the disc?

    If tri is recognised as a Seifert fibre space, then this routine will
    conclusively determine whether it has a Seifert fibration over the disc.

    This routine relies exclusively on combinatorial recognition, so it will
    never be able to prove that tri is not a Seifert fibre space (i.e., it
    will never return SFS.NOTSFS).
    """
    blocked = BlockedSFS.recognise(tri)
    if blocked is None:
        return SFS.UNKNOWN
    sfs = blocked.manifold()
    if sfs.name() == "M/n2 x~ S1":
        # Homeomorphic to SFS [D: (2,1) (2,-1)]. Apart from fibred solid
        # tori, this is the only bounded Seifert fibre space with non-unique
        # Seifert fibration.
        return SFS.DISCTWO
    elif ( sfs.baseOrientable() and sfs.baseGenus() == 0 and
            sfs.punctures() == 1 ):
        # Base is a disc.
        fc = sfs.fibreCount()
        if fc == 2:
            return SFS.DISCTWO
        elif fc == 3:
            return SFS.DISCTHREE
        else:
            return SFS.DISCOTHER
    else:
        # We have a blocked Seifert fibre space whose base is not a disc.
        return SFS.NONDISC


def isFibre(edge):
    """
    Assuming that the given edge belongs to a one-vertex triangulation of a
    small Seifert fibre space that is not also fibred over the projective
    plane, determines whether this edge is isotopic to a Seifert fibre.
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
    classifications = []
    def recognise( sig, tri ):
        c = recogniseSFSOverDisc(tri)
        if c is SFS.UNKNOWN:
            # We don't know if we have a Seifert fibre space.
            return False
        # Since recogniseSFSOverDisc() can never prove that something is not
        # a Seifert fibre space, we now know that this tri is definitely a
        # Seifert fibre space.
        classifications.append(c)
        return True
    height = 0
    nThreads = 1
    if drilled.retriangulate( height, nThreads, recognise ):
        # The drilled triangulation is definitely a Seifert fibre space.
        # Moreover, we know that if the drilled triangulation has a Seifert
        # fibration over the disc, then we will have recognised this.
        if len(classifications) != 1:
            raise ValueError( "Wrong number of classifications: {}".format(
                classifications ) )
        c = classifications[0]
        if c is SFS.DISCTWO:
            # SFS over disc with 2 exceptional fibres. We must have drilled
            # out an exceptional fibre.
            return Fibre.EXCEPTIONAL
        elif c is SFS.DISCTHREE:
            # SFS over disc with 3 exceptional fibres. We must have drilled
            # out a regular fibre.
            return Fibre.REGULAR
        elif c is SFS.NONDISC:
            # SFS, but not over the disc. We must have drilled out a curve
            # that is not isotopic to a Seifert fibre.
            return Fibre.NONFIBRE
        else:
            raise ValueError(
                    "Should never recognise something as {}.".format(c) )

    # Fast tests have failed us. Try looking at normal surfaces.
    surfs = NormalSurfaces.enumerate( drilled, NS_STANDARD )
    annuli = []
    for s in surfs:
        if not s.isOrientable() or not s.hasRealBoundary():
            continue
        euler = s.eulerChar()
        if euler == 1:
            # The surface s is a disc. Determine whether it is a compressing
            # disc by cutting along s and checking whether the torus boundary
            # gets compressed down to a single sphere boundary.
            compress = s.cutAlong()
            compress.intelligentSimplify()
            compress.intelligentSimplify()
            if ( compress.countBoundaryComponents() == 1 and
                    compress.boundaryComponent(0).eulerChar() == 2 ):
                return Fibre.NONFIBRE
        elif euler == 0:
            annuli.append(s)
    if not annuli:
        # The drilled triangulation has no essential annuli, so we cannot
        # have drilled out a Seifert fibre.
        return Fibre.NONFIBRE
    # We are looking for an annulus that either:
    #   --> cuts the drilled triangulation into two solid tori (in which case
    #       we know that we drilled out an exceptional fibre); or
    #   --> cuts the drilled triangulation into a solid torus and a Seifert
    #       fibre space over the disc with 2 exceptional fibres (in which
    #       case we possibly cannot determine whether we drilled out a
    #       Seifert fibre).
    # If we never find such an annulus, then we definitely did not drill out
    # a Seifert fibre.
    for a in annuli:
        # Try cutting along the annulus a. If this causes the drilled
        # triangulation to fall apart into two solid tori, then we know that
        # we drilled out an exceptional fibre.
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
        # exactly one of which is a solid torus. Is the other piece a Seifert
        # fibre space over the disc with 2 exceptional fibres? First try to
        # answer this using combinatorial recognition.
        other = notSolidTorus[0]
        classifications = []
        height = 0
        nThreads = 1
        if other.retriangulate( height, nThreads, recognise ):
            # The other piece is definitely a Seifert fibre space. Moreover,
            # we know that if the other piece has a Seifert fibration over
            # the disc, then we will have recognised this.
            if len(classifications) != 1:
                raise ValueError(
                        "Wrong number of classifications: {}".format(
                            classifications ) )
            c = classifications[0]
            if c is SFS.DISCTWO:
                # SFS over disc with 2 exceptional fibres. We possibly cannot
                # conclude anything.
                return Fibre.UNKNOWN
            elif c is SFS.DISCTHREE:
                # SFS over disc with 3 exceptional fibres.
                # TODO I'm not sure what to do here.
                pass
            elif c is SFS.NONDISC:
                # SFS, but not over the disc.
                # TODO I'm not sure what to do here.
                pass
            else:
                raise ValueError(
                        "Should never recognise something as {} {}.".format(
                            c, "after cutting" ) )

        # Now try looking at normal surfaces in the other triangulation.
        cutSurfs = NormalSurfaces.enumerate( other, NS_STANDARD )
        cutAnnuli = []
        hasCompress = False
        for s in cutSurfs:
            if not s.isOrientable() or not s.hasRealBoundary():
                continue
            euler = s.eulerChar()
            if euler == 1:
                # The surface s is a disc. Determine whether it is a
                # compressing disc by cutting along s and checking whether
                # the torus boundary gets compressed down to a single sphere
                # boundary.
                compress = s.cutAlong()
                compress.intelligentSimplify()
                compress.intelligentSimplify()
                if ( compress.countBoundaryComponents() == 1 and
                        compress.boundaryComponent(0).eulerChar() == 2 ):
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
                # TODO We possibly cannot say anything more definitive.
                return Fibre.UNKNOWN

    # Surviving to this point means that the annulus we were looking for
    # doesn't exist.
    return Fibre.NONFIBRE


# Test code.
if __name__ == "__main__":
    # Generate test triangulations.
    tests = [
#            ( "cPcbbbqxh", "SFS [S2: (2,1) (2,1) (2,-1)]" ),
#            ( "dLQbccchhrw", "SFS [S2: (2,1) (2,1) (3,-2)]" ),
#            ( "eLAkbccddemken", "SFS [S2: (2,1) (2,1) (2,1)]" ),
#            ( "eLAkbccddemkij", "SFS [S2: (2,1) (2,1) (3,-1)] : #1" ),
#            ( "eLPkbcddddrwos", "SFS [S2: (2,1) (2,1) (3,-1)] : #2" ),
#            ( "eLMkbcdddhhhqx", "SFS [S2: (2,1) (2,1) (4,-3)]" ),
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

    # How common are edges that are *not* isotopic to Seifert fibres?
    def findNonFibres(sigs):
        nonFibreSigs = set()
        fibreSigs = set()
        for s in sigs:
            t = Triangulation3.fromIsoSig(s)
            for f in t.triangles():
                result = twoThree(f)
                if result is None:
                    continue
                newt = result[0]
                if isFibre( newt.edge( result[1][-1] ) ) is Fibre.NONFIBRE:
                    nonFibreSigs.add( newt.isoSig() )
                elif not nonFibreSigs:
                    fibreSigs.add( newt.isoSig() )
        if nonFibreSigs:
            return ( True, nonFibreSigs )
        else:
            return ( False, fibreSigs )

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
        print()
        height = 0
        sigSet = {sig}
        newNonFibres = 0
        maxHeight = 6
        maxSize = 15
        while height < maxHeight and len(sigSet) < maxSize:
            height += 1
            foundNew, sigSet = findNonFibres(sigSet)
            if foundNew:
                newNonFibres += 1
            print( "Height {}: Found {} sigs with {} new non-fibres.".format(
                height, len(sigSet), newNonFibres ) )

