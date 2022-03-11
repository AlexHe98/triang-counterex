"""
Functions for determining whether or not a given edge e is isotopic to a
Seifert fibre, where we know in advance that e belongs to a one-vertex
triangulation of a small Seifert fibre space (i.e., a Seifert fibre space
over S^2 with three exceptional fibres).
"""
from regina import *
from pachner import twoThree


def fastNotFibre(ideal):
    """
    Given an ideal triangulation that results from pinching an edge e in a
    one-vertex triangulation of a small Seifert fibre space, uses fast tests
    to see whether we can verify that e is not isotopic to a Seifert fibre.

    Currently, this routine attempts to find a strict angle structure in the
    given ideal triangulation. If such an angle structure exists, then the
    given manifold is hyperbolic, which means it could not have resulted from
    drilling out a Seifert fibre; thus, this routine returns True in this
    case. Otherwise, this routine returns False because it could not make a
    conclusive decision.
    """
    return ideal.hasStrictAngleStructure()


def fastFibre(bounded):
    """
    Given a bounded triangulation that results from drilling out an edge e in
    a one-vertex triangulation of a small Seifert fibre space, uses fast
    tests to see whether we can verify that e is isotopic to a Seifert fibre.

    Currently, this routine uses combinatorial recognition to try to prove
    that the given bounded triangulation represents a Seifert fibre space
    over the disc with either two or three exceptional fibres, which means
    that filling the edge e back into the torus boundary would recover the
    original Seifert fibre space with e as a fibre; thus, this routine
    returns True in this case. Otherwise, this routine returns False because
    it could not make a conclusive decision.
    """
    # Retriangulate with height 0, and check whether combinatorial
    # recognition works for any of the resulting triangulations.
    def recognise( sig, tri ):
        blocked = BlockedSFS.recognise(tri)
        if blocked is None:
            return False
        sfs = blocked.manifold()
        # Verify that we have a Seifert fibre space over the disc with either
        # two or three exceptional fibres.
        differentReps = [
                "M/n2 x~ S1" ] # Homeomorphic to SFS [D: (2,1) (2,-1)].
        if sfs.name() in differentReps:
            return True
        elif ( sfs.baseOrientable() and sfs.baseGenus() == 0 and
                sfs.punctures() == 1 and ( sfs.fibreCount() in {2,3} ) ):
            return True
        else:
            # TODO I'm not sure what to do in this case, so for now I'm just
            #   going to raise an Exception.
            raise ValueError(
                    "Unexpected SFS: {}".format( sfs.detail().rstrip() ) )
    height = 0
    nThreads = 1
    return bounded.retriangulate( height, nThreads, recognise ) 


def annulusNotFibre( bounded ):
    """
    Given a bounded triangulation that results from drilling out an edge e in
    a one-vertex triangulation of a small Seifert fibre space, uses essential
    annuli to see whether we can verify that e is not isotopic to a Seifert
    fibre.
    """
    # Search for annuli among the vertex normal surfaces.
    # TODO
    pass


# Test code.
if __name__ == "__main__":
    # Generate test triangulations.
    tests = [
            ( "cPcbbbqxh", "SFS [S2: (2,1) (2,1) (2,-1)]" ),
            ( "dLQbccchhrw", "SFS [S2: (2,1) (2,1) (3,-2)]" ) ]
    def genDrilled(sig):
        tri = Triangulation3.fromIsoSig(sig)
        # What are the edges currently in the triangulation?
        for i in range( tri.countEdges() ):
            yield Triangulation3(tri), i, "e{}".format(i)
        # What new edges can we introduce using 2-3 moves?
        for i in range( tri.countTriangles() ):
            result = twoThree( tri.triangle(i) )
            if result is None:
                continue
            yield result[0], result[1][-1], "f{}".format(i)

    # Perform tests.
    for sig, name in tests:
        print( sig, name )
        for drilled, i, name in genDrilled(sig):
            drilled.pinchEdge( drilled.edge(i) )
            msg = "    " + name + ": {}"
            if fastNotFibre(drilled):
                print( msg.format("Not a fibre!" ) )
                continue
            drilled.idealToFinite()
            drilled.intelligentSimplify()
            drilled.intelligentSimplify()
            try:
                verifyFibre = fastFibre(drilled)
            except ValueError as err:
                print( msg.format(err) )
            else:
                if verifyFibre:
                    print( msg.format("Fibre!") )
                else:
                    print( msg.format("Inconclusive...") )

                    # Can we find annuli among the vertex normal surfaces?
                    surfs = NormalSurfaces.enumerate( drilled, NS_STANDARD )
                    for s in surfs:
                        if not s.isOrientable():
                            continue
                        if s.countBoundaries() != 2:
                            continue
                        if s.eulerChar() == 0:
                            print( "        {} ...".format(
                                s.detail()[:60] ) )
        print()

