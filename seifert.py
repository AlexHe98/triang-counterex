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
    OTHERSFS = auto()   # SFS, but not one of the types listed below.
    MOBIUSONE = auto()  # SFS over Mobius band, class bn2; 1 exceptional fibre.
    DISCMOBIUS = auto() # SFS over disc that also fibres over Mobius band.
    DISCTWO = auto()    # SFS over disc, class bo1; 2 exceptional fibres.
    DISCTHREE = auto()  # SFS over disc, class bo1; 3 exceptional fibres.


def recogniseSFS(tri):
    """
    Uses combinatorial recognition to determine whether the given
    triangulation is a Seifert fibre space, and if so classify it into one of
    a number of "interesting" types of Seifert fibre space.

    See the enumeration SFS for the list of all possible classifications. We
    currently consider the following types of Seifert fibre space to be
    "interesting":
    --> Seifert fibre spaces of class bo1 over the disc, with two
        exceptional fibres;
    --> Seifert fibre spaces of class bo1 over the disc with three
        exceptional fibres; and
    --> Seifert fibre spaces of class bn2 over the Mobius band with one
        exceptional fibre.
    Note that there is exactly one Seifert fibre space over the disc that
    also fibres over the Mobius band; this has class bo1 and two exceptional
    fibres over the disc, and has class bn2 and no exceptional fibres over
    the Mobius band.

    This routine relies exclusively on combinatorial recognition, so it will
    never be able to prove that tri is not a Seifert fibre space (i.e., it
    will never return SFS.NOTSFS).
    """
    blocked = BlockedSFS.recognise(tri)
    if blocked is None:
        return SFS.UNKNOWN
    sfs = blocked.manifold()
    if sfs.punctures() != 1:
        # Definitely not fibred over the disc or the Mobius band.
        return SFS.OTHERSFS
    if sfs.baseClass() == SFSpace.bo1 and sfs.baseGenus() == 0:
        # Fibred over disc, class bo1.
        if sfs.fibreCount() == 2:
            # When we have two exceptional fibres, we could also be fibred
            # over the Mobius band.
            if ( sfs.fibre(0).alpha == 2 and sfs.fibre(0).beta == 1 and
                    sfs.fibre(1).alpha == 2 and sfs.fibre(1).beta == 1 ):
                return SFS.DISCMOBIUS
            else:
                return SFS.DISCTWO
        elif sfs.fibreCount() == 3:
            return SFS.DISCTHREE
        else:
            return SFS.OTHERSFS
    elif sfs.baseClass() == SFSpace.bn2 and sfs.baseGenus() == 1:
        # Fibred over Mobius band, class bn2.
        if sfs.fibreCount() == 0:
            return SFS.DISCMOBIUS
        elif sfs.fibreCount() == 1:
            return SFS.MOBIUSONE
        else:
            return SFS.OTHERSFS
    else:
        return SFS.OTHERSFS


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
    # Retriangulate with height 0 (single-threaded), and see if we can ever
    # recognise the drilled triangulation as a Seifert fibre space.
    classifications = []
    def recognise( sig, tri ):
        c = recogniseSFS(tri)
        if c is SFS.UNKNOWN:
            # We don't know if we have a Seifert fibre space.
            return False
        elif c is SFS.NOTSFS:
            # Combinatorial recognition should never be able to prove that
            # tri is not a Seifert fibre space.
            raise ValueError( "Recognised {}".format(c) )
        # We definitely have a Seifert fibre space.
        classifications.append(c)
        return True
    if drilled.retriangulate( 0, 1, recognise ):
        # The drilled triangulation is definitely a Seifert fibre space.
        if len(classifications) != 1:
            raise ValueError( "Wrong number of classifications: {}".format(
                classifications ) )
        c = classifications[0]
        if ( c is SFS.DISCMOBIUS ) or ( c is SFS.DISCTWO ):
            return Fibre.EXCEPTIONAL
        elif ( c is SFS.DISCTHREE ) or ( c is SFS.MOBIUSONE ):
            return Fibre.REGULAR
        elif c is SFS.OTHERSFS:
            return Fibre.NONFIBRE
        else:
            raise ValueError(
                    "Should never recognise something as {}".format(c) )

    # Fast tests have failed us. Try looking at normal surfaces.
    # Look for essential annuli. By Corollary 6.8 of "Algorithms for the
    # complete decomposition of a closed 3-manifold" (Jaco and Tollefson,
    # 1995), such annuli are guaranteed to appear either as a two-sided
    # vertex normal surface or as the double of a one-sided vertex normal
    # surface (the latter case is necessary because Jaco and Tollefson do not
    # consider one-sided surfaces to be vertex surfaces).
    surfs = NormalSurfaces.enumerate( drilled, NS_STANDARD )
    annuli = []
    for s in surfs:
        if not s.hasRealBoundary():
            if not s.isOrientable() or s.eulerChar() != 0:
                continue
            # TODO Experiment with vertex normal tori.
            # We have a torus.
#            print( "            Torus: {}".format( s.detail()[:60] ) )
            continue
        euler = s.eulerChar()
        if not s.isOrientable():
            # The surface s is bounded and non-orientable. If it is a Mobius
            # band, then we need to double it.
            if euler != 0:
                continue
            doubleSurf = s.doubleSurface()
            # We know for sure that doubleSurf has real boundary and has
            # Euler characteristic 0. We just need to check that it is
            # connected and orientable.
            if doubleSurf.isOrientable() and doubleSurf.isConnected():
                annuli.append(doubleSurf)
            else:
                raise ValueError( "Unexpected double of Mobius band." )
            continue
        if euler == 1:
            # The surface s is a disc. Determine whether it is a compressing
            # disc by cutting along s and checking whether the torus boundary
            # gets compressed down to a single sphere boundary.
            compress = s.cutAlong()
            compress.intelligentSimplify()
            compress.intelligentSimplify()
            if ( compress.countBoundaryComponents() == 1 and
                    compress.boundaryComponent(0).eulerChar() == 2 ):
                # TODO Test.
#                print( "        Compressing disc!" )
                return Fibre.NONFIBRE
        elif euler == 0:
            annuli.append(s)
    if not annuli:
        # The drilled triangulation has no essential annuli, so we cannot
        # have drilled out a Seifert fibre.
        # TODO Test.
#        print( "        No annuli!" )
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
    otherPieces = []
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
        # exactly one of which is a solid torus. Deal with the other piece
        # later.
        otherPieces.append( notSolidTorus[0] )

    # We definitely didn't drill out an exceptional fibre. Did we perhaps
    # drill out a regular fibre instead?
    for other in otherPieces:
        # Retriangulate with height 0 (single-threaded), and see if we can
        # recognise this other piece as a Seifert fibre space.
        classifications = []
        if other.retriangulate( 0, 1, recognise ):
            # This other piece is definitely a Seifert fibre space.
            if len(classifications) != 1:
                raise ValueError(
                        "Wrong number of classifications: {}".format(
                            classifications ) )
            c = classifications[0]
            if ( c is SFS.DISCMOBIUS ) or ( c is SFS.DISCTWO ):
                # We cannot say anything definitive without doing more work.
                return Fibre.UNKNOWN
            elif ( c is SFS.DISCTHREE ) or ( c is SFS.MOBIUSONE ):
                # We must have drilled out a regular fibre and then cut along
                # a boundary-parallel annulus.
                # TODO Check this.
                return Fibre.REGULAR
            elif c is SFS.OTHERSFS:
                # We cannot get this if we drilled out a Seifert fibre.
                # TODO Check this.
                return Fibre.NONFIBRE
            else:
                raise ValueError(
                        "Should never recognise something as {} {}".format(
                            c, "after cutting" ) )

        # Now try looking at normal surfaces in this other piece.
        cutSurfs = NormalSurfaces.enumerate( other, NS_STANDARD )
        cutAnnuli = []
        hasCompress = False
        for s in cutSurfs:
            if not s.hasRealBoundary():
                continue
            euler = s.eulerChar()
            if not s.isOrientable():
                # The surface is bounded and non-orientable. If it is a
                # Mobius band, then we need to double it.
                if euler != 0:
                    continue
                doubleSurf = s.doubleSurface()
                # We know for sure that doubleSurf has real boundary and has
                # Euler characteristic 0. We just need to check that it is
                # connected and orientable.
                if doubleSurf.isOrientable() and doubleSurf.isConnected():
                    cutAnnuli.append(doubleSurf)
                else:
                    raise ValueError( "Unexpected double of Mobius band." )
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
                # We cannot say anything definitive without doing more work.
                return Fibre.UNKNOWN

    # Surviving to this point means that the annulus we were looking for
    # doesn't exist.
#    print( "        Survived cutting!" )
    return Fibre.NONFIBRE


def isExceptionalFibre(edge):
    """
    Assuming that the given edge belongs to a one-vertex triangulation of a
    small Seifert fibre space, determines whether this edge is isotopic to an
    exceptional fibre.
    """
    # Drill out the given edge.
    drilled = Triangulation3( edge.triangulation() )
    drilled.pinchEdge( drilled.edge( edge.index() ) )
    drilled.idealToFinite()
    drilled.intelligentSimplify()
    drilled.intelligentSimplify()

    # Retriangulate with height 0 (single-threaded), and see if we can ever
    # recognise the drilled triangulation as a Seifert fibre space.
    classifications = []
    def recognise( sig, tri ):
        c = recogniseSFS(tri)
        if c is SFS.UNKNOWN:
            # We don't know if we have a Seifert fibre space.
        elif c is SFS.NOTSFS:
            # Combinatorial recognition should never be able to prove that
            # tri is not a Seifert fibre space.
            raise ValueError( "Recognises {}".format(c) )
        # We definitely have a Seifert fibre space.
        classifications.append(c)
        return True
    if drilled.retriangulate( 0, 1, recognise ):
        # The drilled triangulation is definitely a Seifert fibre space.
        if len(classifications) != 1:
            raise ValueError( "Wrong number of classifications: {}".format(
                classifications ) )
        c = classifications[0]
        return ( ( c is SFS.DISCMOBIUS ) or ( c is SFS.DISCTWO ) )

    # Look for essential annuli. By Corollary 6.8 of "Algorithms for the
    # complete decomposition of a closed 3-manifold" (Jaco and Tollefson,
    # 1995), such annuli are guaranteed to appear either as a two-sided
    # vertex normal surface or as the double of a one-sided vertex normal
    # surface (the latter case is necessary because Jaco and Tollefson do not
    # consider one-sided surfaces to be vertex surfaces).
    surfs = NormalSurfaces.enumerate( drilled, NS_STANDARD )
    annuli = []
    for s in surfs:
        if not s.hasRealBoundary():
            continue
        euler = s.eulerChar()
        if not s.isOrientable():
            # The surface is bounded and non-orientable. If it is a Mobius
            # band, then we need to double it.
            if euler != 0:
                continue
            doubleSurf = s.doubleSurface()
            # We know for sure that doubleSurf has real boundary and has
            # Euler characteristic 0. For peace of mind, we double check that
            # it is connected and orientable.
            if doubleSurf.isOrientable() and doubleSurf.isConnected():
                annuli.append(doubleSurf)
            else:
                raise ValueError( "Unexpected double of Mobius band." )
            continue
        if euler == 1:
            # The surface s is a disc. Determine whether it is a compressing
            # disc by cutting along s and checking whether the torus boundary
            # gets compressed down to a single sphere boundary.
            compress = s.cutAlong()
            compress.intelligentSimplify()
            compress.intelligentSimplify()
            if ( compress.countBoundaryComponents() == 1 and
                    compress.boundaryComponent(0).eulerChar() == 2 ):
                return False
        elif euler == 0:
            annuli.append(s)
    if not annuli:
        # The drilled triangulation has no essential annuli, so we cannot
        # have drilled out an exceptional fibre.
        return False
    # We are looking for an annulus that cuts the drilled triangulation into
    # two solid tori.
    for a in annuli:
        # Try cutting along the annulus a. If this causes the drilled
        # triangulation to fall apart into two solid tori, then we know that
        # we drilled out an exceptional fibre.
        cut = a.cutAlong()
        if cut.countComponents() != 2:
            # Cutting along a only yields one piece, so it isn't the annulus
            # we're looking for. Move along to the next annulus.
            continue
        onlySolidTori = True
        for comp in cut.triangulateComponents():
            comp.intelligentSimplify()
            comp.intelligentSimplify()
            if not comp.isSolidTorus():
                onlySolidTori = False
                break
        if onlySolidTori:
            return True
    # Surviving to this point means that the annulus we were looking doesn't
    # exist, which implies that we did not drill out an exceptional fibre.
    return False


# Test code.
if __name__ == "__main__":
    # Generate test triangulations.
    tests = [
#            ( "cPcbbbqxh", "SFS [S2: (2,1) (2,1) (2,-1)]" ) ]
            ( "cPcbbbqxh", "SFS [S2: (2,1) (2,1) (2,-1)]" ),
            ( "dLQbccchhrw", "SFS [S2: (2,1) (2,1) (3,-2)]" ),
            ( "eLAkbccddemken", "SFS [S2: (2,1) (2,1) (2,1)]" ),
            ( "eLAkbccddemkij", "SFS [S2: (2,1) (2,1) (3,-1)] : #1" ),
            ( "eLPkbcddddrwos", "SFS [S2: (2,1) (2,1) (3,-1)] : #2" ),
            ( "eLMkbcdddhhhqx", "SFS [S2: (2,1) (2,1) (4,-3)]" ),
#            ( "eLPkbcdddhrrnk", "SFS [S2: (2,1) (3,1) (3,-2)]" ) ]
            ( "eLPkbcdddhrrnk", "SFS [S2: (2,1) (3,1) (3,-2)]" ),
            ( "fLLQcaceeedjkuxkn", "SFS [S2: (2,1) (3,1) (3,-1)] : #1" ),
            ( "fLLQcbeddeehhokum", "SFS [S2: (2,1) (3,1) (3,-1)] : #2" ),
            ( "fLLQcbeddeehhnkxx", "SFS [S2: (2,1) (3,1) (4,-3)] : #1" ),
            ( "fvPQcceddeerrnskr", "SFS [S2: (2,1) (3,1) (4,-3)] : #2" ),
            ( "fvPQcdecedekrsnrs", "SFS [S2: (2,1) (3,1) (5,-4)]" ),
            ( "fLLQcaceeedjkuxkj", "SFS [S2: (2,1) (3,2) (3,-1)]" ) ]
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
        maxHeight = 4
        maxSize = 15
        while height < maxHeight and len(sigSet) < maxSize:
            height += 1
            foundNew, sigSet = findNonFibres(sigSet)
            if foundNew:
                newNonFibres += 1
            print( "Height {}: Found {} sigs with {} new non-fibres.".format(
                height, len(sigSet), newNonFibres ) )

