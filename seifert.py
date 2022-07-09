"""
Functions for determining whether or not a given edge e is isotopic to a
Seifert fibre, where we know in advance that e belongs to a one-vertex
triangulation of a small Seifert fibre space (i.e., a Seifert fibre space
over S^2 with three exceptional fibres).
"""
from sys import stdout
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
    ANNULUSONE = auto() # SFS over annulus, class bo1; 1 exceptional fibre.


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
        exceptional fibres;
    --> Seifert fibre spaces of class bn2 over the Mobius band with one
        exceptional fibre; and
    --> Seifert fibre spaces of class bo1 over the annulus, with one
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
    punctures = sfs.punctures()
    baseClass = sfs.baseClass()
    baseGenus = sfs.baseGenus()
    fibreCount = sfs.fibreCount()
    if punctures == 2:
        # Are we fibred over the annulus, class bo1, with one exceptional
        # fibre?
        if baseClass == SFSpace.bo1 and baseGenus == 0 and fibreCount == 1:
            return SFS.ANNULUSONE
        else:
            return SFS.OTHERSFS
    elif punctures != 1:
        # Definitely not fibred over the disc, annulus or Mobius band.
        return SFS.OTHERSFS
    elif baseClass == SFSpace.bo1 and baseGenus == 0:
        # Fibred over disc, class bo1.
        if fibreCount == 2:
            # When we have two exceptional fibres, we could also be fibred
            # over the Mobius band.
            if ( sfs.fibre(0).alpha == 2 and sfs.fibre(0).beta == 1 and
                    sfs.fibre(1).alpha == 2 and sfs.fibre(1).beta == 1 ):
                return SFS.DISCMOBIUS
            else:
                return SFS.DISCTWO
        elif fibreCount == 3:
            return SFS.DISCTHREE
        else:
            return SFS.OTHERSFS
    elif baseClass == SFSpace.bn2 and baseGenus == 1:
        # Fibred over Mobius band, class bn2.
        if fibreCount == 0:
            return SFS.DISCMOBIUS
        elif fibreCount == 1:
            return SFS.MOBIUSONE
        else:
            return SFS.OTHERSFS
    else:
        return SFS.OTHERSFS


def cutOneSolidTorus(tri):
    """
    Does the given triangulation represent a boundary-irreducible 3-manifold
    that can be decomposed into a single solid torus by cutting along a
    (non-separating) essential annulus?

    Pre-condition:
    --> The given triangulation is an orientable and irreducible 3-manifold
        with exactly two real torus boundary components (and no other
        boundary components).
    """
    surfs = NormalSurfaces( tri, NS_STANDARD )
    annuli = []
    for s in surfs:
        if not s.hasRealBoundary():
            continue
        euler = s.eulerChar()
        if not s.isOrientable():
            # Is the double of s an annulus?
            if euler != 0:
                continue
            doubleSurf = s.doubleSurface()
            if doubleSurf.isOrientable() and doubleSurf.isConnected():
                annuli.append(doubleSurf)
            else:
                raise ValueError( "Unexpected double of Mobius band" )
        elif euler == 1:
            # Is s a compressing disc?
            compress = s.cutAlong()
            compress.intelligentSimplify()
            compress.intelligentSimplify()
            if compress.countBoundaryComponents() != 2:
                continue
            for b in compress.boundaryComponents():
                if b.eulerChar() == 2:
                    return False
        elif euler == 0:
            if s.isThinEdgeLink()[0] is None:
                annuli.append(s)
    if not annuli:
        return False
    for a in annuli:
        c = a.cutAlong()
        if c.isSolidTorus():
            return True
    return False


def cutTwoSolidTori(tri):
    """
    Does the given triangulation represent a boundary-irreducible 3-manifold
    that can be decomposed into two solid tori by cutting along a
    (separating) essential annulus?

    Pre-condition:
    --> The given triangulation is an orientable and irreducible 3-manifold
        with exactly one real torus boundary component (and no other boundary
        components).
    """
    surfs = NormalSurfaces( tri, NS_STANDARD )
    annuli = []
    for s in surfs:
        if not s.hasRealBoundary():
            continue
        euler = s.eulerChar()
        if not s.isOrientable():
            # Is the double of s an annulus?
            if euler != 0:
                continue
            doubleSurf = s.doubleSurface()
            if doubleSurf.isOrientable() and doubleSurf.isConnected():
                annuli.append(doubleSurf)
            else:
                raise ValueError( "Unexpected double of Mobius band" )
        elif euler == 1:
            # Is s a compressing disc?
            compress = s.cutAlong()
            compress.intelligentSimplify()
            compress.intelligentSimplify()
            if ( compress.countBoundaryComponents() == 1 and
                    compress.boundaryComponent(0).eulerChar() == 2 ):
                return False
        elif euler == 0:
            if s.isThinEdgeLink()[0] is None:
                annuli.append(s)
    if not annuli:
        return False
    for a in annuli:
        c = a.cutAlong()
        if c.countComponents() != 2:
            continue
        onlySolidTori = True
        for cc in c.triangulateComponents():
            cc.intelligentSimplify()
            cc.intelligentSimplify()
            if not cc.isSolidTorus():
                onlySolidTori = False
                break
        if onlySolidTori:
            return True
    return False


def isFibre(edge):
    """
    Assuming that the given edge belongs to a one-vertex triangulation of a
    small Seifert fibre space (i.e., a Seifert fibre space over the sphere
    three exceptional fibres), determines whether this edge is isotopic to a
    Seifert fibre in the fibring over the disc (note that the space could
    hava a different Seifert fibring over the projective plane, but this
    routine will *not* consider the regular fibres of such a fibring to be a
    fibre).
    """
    # Pinch the given edge. If the resulting ideal triangulation admits a
    # strict angle structure, then the edge cannot be a Seifert fibre.
    drilled = Triangulation3( edge.triangulation() )
    drilled.pinchEdge( drilled.edge( edge.index() ) )
    drilled.intelligentSimplify()
    drilled.intelligentSimplify()
    if drilled.hasStrictAngleStructure():
        return Fibre.NONFIBRE

    # Now truncate, which corresponds to drilling out the edge that we
    # pinched. If this edge was isotopic to a Seifert fibre, then the
    # drilled triangulation will be a Seifert fibre space over the disc with
    # either two or three exceptional fibres. Our goal is to show that the
    # drilled triangulation is not a Seifert fibre space of the correct type,
    # and hence that we did not drill out a Seifert fibre.
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
        if c in { SFS.DISCMOBIUS, SFS.DISCTWO, SFS.DISCTHREE }:
            # The drilled triangulation is a Seifert fibre space of the
            # correct type, so we cannot rule out the possibility that we
            # drilled out a Seifert fibre.
            return Fibre.UNKNOWN
        elif c in { SFS.OTHERSFS, SFS.MOBIUSONE }:
            # The drilled triangulation is a Seifert fibre space, but not of
            # the correct type, so we definitely didn't drill out a Seifert
            # fibre.
            return Fibre.NONFIBRE
        else:
            raise ValueError(
                    "Unexpected {}".format(c) )

    # Fast tests have failed us. Try looking at embedded surfaces:
    #   --> If the drilled triangulation is a Seifert fibre space of the
    #       correct type, then it cannot contain a compressing disc. Thus,
    #       finding such a disc tells us that we didn't drill out a Seifert
    #       fibre.
    #   --> If the drilled triangulation is fibred over the disc with two
    #       exceptional fibres, then there is an essential annulus that cuts
    #       the triangulation into two solid tori. Thus, showing that no such
    #       annulus exists would tell us that we didn't drill out an
    #       exceptional fibre.
    #   --> If the drilled triangulation is fibred over the disc with three
    #       exceptional fibres, then there exist:
    #       (a) three essential annuli, each of which cuts the triangulation
    #           into a solid torus and a Seifert fibre space over the disc
    #           with two exceptional fibres; and
    #       (b) three essential tori, each of which cuts the triangulation
    #           a Seifert fibre space over the disc with two exceptional
    #           fibres and a Seifert fibre space over the disc with one
    #           exceptional fibre.
    # We use the following facts, which are consequences of various results
    # from Section 6 of the paper "Algorithms for the complete decomposition
    # of a closed 3-manifold" (Jaco and Tollefson, 1995):
    #   --> Since our drilled triangulation is guaranteed to represent a
    #       compact and irreducible 3-manifold, we know that if it has a
    #       compressing disc then at least one such disc must appear among
    #       the vertex normal surfaces.
    #   --> Since the Seifert fibre spaces of interest are orientable,
    #       compact, irreducible and boundary-irreducible, we know that at
    #       least one of the essential annuli of interest must appear as
    #       either: a vertex normal annulus, or the double of a vertex normal
    #       Mobius band.
    #   --> Similarly, in the Seifert fibre spaces of interest, we know that
    #       at least two of the essential tori of interest must appear as
    #       either: a vertex normal torus, or the double of a vertex normal
    #       Klein bottle.
    surfs = NormalSurfaces( drilled, NS_STANDARD )
    annuli = []
    tori = []
    for s in surfs:
        euler = s.eulerChar()
        if not s.hasRealBoundary():
            # Is s a torus? What about its double?
            if euler != 0:
                continue
            if s.isOrientable():
                # We have a torus (but ignore it if it is obviously trivial).
                if s.isThinEdgeLink()[0] is None:
                    tori.append(s)
            else:
                # We have a Klein bottle. Double it (we should get a torus,
                # but do some sanity checks).
                doubleSurf = s.doubleSurface()
                if doubleSurf.isOrientable() and doubleSurf.isConnected():
                    tori.append(doubleSurf)
                else:
                    raise ValueError( "Unexpected double of Klein bottle" )
        elif not s.isOrientable():
            # Is the double of s an annulus?
            if euler != 0:
                continue
            doubleSurf = s.doubleSurface()
            if doubleSurf.isOrientable() and doubleSurf.isConnected():
                annuli.append(doubleSurf)
            else:
                raise ValueError( "Unexpected double of Mobius band" )
        elif euler == 1:
            # Is s a compressing disc?
            compress = s.cutAlong()
            compress.intelligentSimplify()
            compress.intelligentSimplify()
            if ( compress.countBoundaryComponents() == 1 and
                    compress.boundaryComponent(0).eulerChar() == 2 ):
                # Found a compressing disc in the drilled triangulation. This
                # means that we didn't drill out a Seifert fibre.
                return Fibre.NONFIBRE
        elif euler == 0:
            if s.isThinEdgeLink()[0] is None:
                annuli.append(s)
    if not annuli:
        # The drilled triangulation has no essential annuli, so we cannot
        # have drilled out a Seifert fibre.
        return Fibre.NONFIBRE

    # If none of the annuli cut the drilled triangulation into two solid
    # tori, then we definitely didn't drill out an exceptional fibre. While
    # we check this, also pay attention to annuli that cut the drilled
    # triangulation into two pieces of which only one is a solid torus, since
    # this could tell us whether we drilled out a regular fibre.
    otherPieces = []
    for a in annuli:
        cut = a.cutAlong()
        if cut.countComponents() != 2:
            continue

        # Cutting along a decomposed the drilled triangulation into two
        # pieces. How many of these pieces are solid tori?
        notSolidTorus = []
        for comp in cut.triangulateComponents():
            comp.intelligentSimplify()
            comp.intelligentSimplify()
            if not comp.isSolidTorus():
                notSolidTorus.append(comp)
        nonSolidTorusCount = len(notSolidTorus)
        if nonSolidTorusCount == 0:
            # The drilled triangulation decomposes into two solid tori, so we
            # cannot rule out the possibility that we drilled out an
            # exceptional fibre.
            return Fibre.UNKNOWN
        elif nonSolidTorusCount == 1:
            # We cut off a single solid torus from the drilled triangulation.
            # The topology of the other piece could tell us whether we
            # drilled out a regular fibre.
            otherPieces.append( notSolidTorus[0] )

    # We didn't drill out an exceptional fibre, but did we perhaps drill out
    # a regular fibre instead?
    if len(tori) < 2:
        # As mentioned above, we must find at least two essential tori among
        # the vertex normal surfaces.
        return Fibre.NONFIBRE

    # If we drilled out a regular fibre, at least one of the other pieces
    # must be fibred over the disc with two exceptional fibres.
    foundOther = False
    for other in otherPieces:
        # Start with combinatorial recognition.
        classifications = []
        if other.retriangulate( 0, 1, recognise ):
            # This other piece is definitely a Seifert fibre space.
            if len(classifications) != 1:
                raise ValueError(
                        "Wrong number of classifications: {}".format(
                            classifications ) )
            c = classifications[0]
            if c in { SFS.DISCMOBIUS, SFS.DISCTWO }:
                # Found a piece of the required type.
                foundOther = True
                break
            elif c in { SFS.DISCTHREE, SFS.OTHERSFS, SFS.MOBIUSONE }:
                continue
            else:
                raise ValueError(
                        "Unexpected {}".format(c) )

        # Now try looking at normal surfaces in this other piece.
        if cutTwoSolidTori(other):
            # We found a boundary-irreducible piece that can be decomposed
            # into two solid tori by cutting along an essential annulus.
            foundOther = True
            break
    if not foundOther:
        # We did not find a piece of the required type, so we could not have
        # drilled out a regular fibre.
        return Fibre.NONFIBRE

    # We are still not sure whether we drilled out a regular fibre. We use
    # one final test: if we did drill out a regular fibre, then we have at
    # least two vertex normal essential tori that cut the drilled
    # triangulation into:
    #   --> a Seifert fibre space over the disc with two exceptional fibres;
    #       and
    #   --> a Seifert fibre space over the annulus with one exceptional
    #       fibre.
    torusCount = 0
    for t in tori:
        cut = t.cutAlong()
        if cut.countComponents() != 2:
            continue
        foundDiscPiece = False
        foundAnnulusPiece = False
        for piece in cut.triangulateComponents():
            piece.intelligentSimplify()
            piece.intelligentSimplify()
            bdries = piece.countBoundaryComponents()
            if bdries == 1:
                # Is this piece the required Seifert fibre space over the
                # disc? Start with combinatorial recognition.
                classifications = []
                if piece.retriangulate( 0, 1, recognise ):
                    # This piece is definitely a Seifert fibre space.
                    if len(classifications) != 1:
                        raise ValueError(
                                "Wrong number of classifications: {}".format(
                                    classifications ) )
                    c = classifications[0]
                    if c in { SFS.DISCMOBIUS, SFS.DISCTWO }:
                        foundDiscPiece = True
                        continue
                    elif c in { SFS.DISCTHREE, SFS.OTHERSFS, SFS.MOBIUSONE }:
                        break
                    else:
                        raise ValueError(
                                "Unexpected {}".format(c) )

                # Now try looking at normal surfaces in this piece.
                if cutTwoSolidTori(piece):
                    foundDiscPiece = True
            elif bdries == 2:
                # Is this piece the required Seifert fibre space over the
                # annulus? Start with combinatorial recognition.
                classifications = []
                if piece.retriangulate( 0, 1, recognise ):
                    # This piece is definitely a Seifert fibre space.
                    if len(classifications) != 1:
                        raise ValueError(
                                "Wrong number of classifications: {}".format(
                                    classifications ) )
                    c = classifications[0]
                    if c is SFS.ANNULUSONE:
                        foundAnnulusPiece = True
                        continue
                    elif c is SFS.OTHERSFS:
                        break
                    else:
                        raise ValueError(
                                "Unexpected {}".format(c) )

                # Now try looking at normal surfaces in this piece.
                if cutOneSolidTorus(piece):
                    foundAnnulusPiece = True
            else:
                raise ValueError(
                        "Piece with {} boundary components".format(bdries) )

        # Did cutting along t give two pieces of the required types? If so,
        # then t is one of the two tori of the required type.
        if foundDiscPiece and foundAnnulusPiece:
            torusCount += 1
        if torusCount == 2:
            break
    if torusCount < 2:
        # We did not find two tori of the required type, so we could not have
        # drilled out a regular fibre.
        return Fibre.NONFIBRE

    # If we survived to this point, then nothing we did could rule out the
    # possibility that we drilled out a regular fibre.
    return Fibre.UNKNOWN


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
            return False
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
    surfs = NormalSurfaces( drilled, NS_STANDARD )
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

    # How common are edges that are *not* isotopic to exceptional fibres?
    def findNonExcFibres(sigs):
        nonFibreSigs = set()
        fibreSigs = set()
        for s in sigs:
            t = Triangulation3.fromIsoSig(s)
            for f in t.triangles():
                result = twoThree(f)
                if result is None:
                    continue
                newt = result[0]
                if isExceptionalFibre( newt.edge( result[1][-1] ) ):
                    fibreSigs.add( newt.isoSig() )
                else:
                    nonFibreSigs.add( newt.isoSig() )
        if nonFibreSigs:
            return ( True, nonFibreSigs )
        else:
            return ( False, fibreSigs )

    # Perform tests.
#    for sig, name in tests:
#        print()
#        print( sig, name )
#        for tri, i, name in genEdges(sig):
#            msg = "    " + name + ": {}"
#            fibreType = isFibre( tri.edge(i) )
#            print( msg.format( fibreType.name ) )
#            stdout.flush()

#    for sig, name in tests:
#        print()
#        print( sig, name )
#        for tri, i, name in genEdges(sig):
#            msg = "    " + name + ": {}"
#            try:
#                fibreType = isFibre( tri.edge(i) )
#            except ValueError as err:
#                print( msg.format( err ) )
#            else:
#                if fibreType is Fibre.NONFIBRE:
#                    print( msg.format( fibreType.name + "   <--" ) )
#                else:
#                    print( msg.format( fibreType.name ) )
#            stdout.flush()
#        print()
#        height = 0
#        sigSet = {sig}
#        newNonFibres = 0
#        maxHeight = 3
#        maxSize = 10
#        while height < maxHeight and len(sigSet) < maxSize:
#            height += 1
#            foundNew, sigSet = findNonFibres(sigSet)
#            if foundNew:
#                newNonFibres += 1
#            print( "Height {}: Found {} sigs with {} new non-fibres.".format(
#                height, len(sigSet), newNonFibres ) )
#            stdout.flush()
#        print()
#        height = 0
#        sigSet = {sig}
#        newNonFibres = 0
#        maxHeight = 3
#        maxSize = 10
#        while height < maxHeight and len(sigSet) < maxSize:
#            height += 1
#            foundNew, sigSet = findNonExcFibres(sigSet)
#            if foundNew:
#                newNonFibres += 1
#            print( "Height {}: Found {} sigs with {} new {}.".format(
#                height, len(sigSet), newNonFibres, "non-exc-fibres" ) )
#            stdout.flush()

    # Tests for SFS with no edges isotopic to exceptional fibres.
    noExcFibres = [
            # No exceptional fibres.
#            "oLvPvLQLQQcccgkhjlkmknlmnnhrauchalaahjggf",
#            "oLvPvLMPQQcccgljlmknnjkkmnqjaqlaaiqaidkgb",
            # No fibres at all.
            "mLLvAwAQQcdfehijjlklklhsaagatthofwj",
            "nLLvPPvQQkcdegijlmjlkmmlhsagvaahhluovn",
#            "oLLLvMMLQQcbcgijkilmlnlmnnlsmjaxftatvfrcv",
            "nLvPwLzQQkccgfiikjmklmlmhnahlupmtrsvgb",
            "nLLLvLQQPkccfghililkmklmlnacnbdwathjsn",
            "nLLLvLQQAkccfghliiljlkmmlnawnpjsqqjjxr",
#            "oLALLvALQQccbcegjkjlnnmnmmudbsaausjjckrkw",
            "lLLvPAPQccdeghiihkkkjhsaggvgndqkj",
#            "mLLvPzPQQcdegiililjlkkhsagvgxgqhjsj" ]
            "mLLvPzPQQcdegiililjlkkhsagvgxgqhjsj",
#            "oLLLMwMMLQccdgfghjikknmnmnhshagcqiacinggn",
#            "oLLLMwMMwQccdgfghjikklnmnnhshagcqiacecrrf",
#            "oLLLMwwLQQccdgfghjlmnmnlmnhshagcqrclruobm" ]
#            "oLLLMwwLQQccdgfghjlmnmnlmnhshagcqrclruobm",
#            "oLLLwwAPPQccdgfhkilmkmlnnnhshacqaacnccxhu" ]
#            "oLLLwwAPPQccdgfhkilmkmlnnnhshacqaacnccxhu",
#            "oLLLMwLPAQccdgfghkjkmllmnnhshagcabrcoobej",
#            "oLLLMwLQwQccdgfghilklkmmnnhshagciaccwaarj" ]
#            "oLLvMvMQQMccdfgilkljkjmmnnhsegspmimddggpw",
#            "oLLvvLQPQQccdimjkhnmnllmlnhsvasantlauuusr",
#            "oLLvvLQPQQccdimkkhnmlllnmnhsvaaantlpuuuur",
#            "oLLLwvAPQQccdgfhijklnmnmnmhshesadccjcvcnk",
#            "oLLvLLAAQQccdgiljhmnkklmnnhsgvasndhuuubnn" ]
            "nLLvMwzQQkcdfgikkmjljmlmhsagctthtaegwj",
            "nLvvAAwQQkcgfhjjlmkmkjlmhahgooatptigof" ]
#            "oLLvLLLQQQccdgimlmlkkjnknnhsgvtptiipolpof",
#            "oLLvLMLPQQccdgjimmnkljknnmhsgcviihttoeggo" ]
#            "oLLvMvAQQPccdfgilkmjmjllnnhsegstmhmqdffdk" ]

    for sig in noExcFibres:
        print()
        tri = Triangulation3.fromIsoSig(sig)
        simplified = Triangulation3(tri)
        simplified.intelligentSimplify()
        simplified.intelligentSimplify()
        blocked = BlockedSFS.recognise(simplified)
        print(sig)
        if blocked is not None:
            print( blocked.manifold().detail().rstrip() )
#        census = Census.lookup(simplified)
#        print(census)
        for i in range( tri.countEdges() ):
            msg = "    e{}".format(i) + ": {}"
            RandomEngine.reseedWithHardware()
            try:
                fibreType = isFibre( tri.edge(i) )
            except ValueError as err:
                print( msg.format(err) )
            else:
                output = True
                if fibreType is Fibre.UNKNOWN:
#                    print( "            Try again!" )
                    RandomEngine.reseedWithHardware()
                    try:
                        fibreType = isFibre( tri.edge(i) )
                    except ValueError as err:
                        print( msg.format(err) )
                        output = False
                if output:
                    if fibreType is Fibre.NONFIBRE:
                        print( msg.format(
                            "--> " + fibreType.name + " <--" ) )
                    else:
                        print( msg.format( fibreType.name ) )
            stdout.flush()

