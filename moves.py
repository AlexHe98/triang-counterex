"""
Perform 2-3, 3-2 and 2-0 moves while tracking how edges are relabelled.
"""
from sys import argv, stdout
from regina import *


def twoThree(triangle):
    """
    Performs a 2-3 move about the given triangle, and returns a dictionary r
    that describes how the edges were renumbered.

    This routine directly modifies the triangulation that contains the given
    triangle. If an edge is currently numbered i in the triangulation, then
    it will be numbered r[i] in the triangulation after the requested 2-3
    move has been performed (provided the move is actually legal); also, a
    2-3 move always creates a new edge, and the index of this new edge will
    be given by r[-1]. If the move is not legal, the triangulation is left
    untouched and this routine returns None.

    If the triangulation containing the given triangle is currently oriented,
    then this orientation will be preserved by the requested 2-3 move.
    """
    tri = triangle.triangulation()

    # Is the requested 2-3 move legal?
    if not tri.pachner(triangle, True, False):
        return None

    # We need to work out the gluings that we need to perform before we make
    # any changes.
    doomed = []
    verts = []
    for emb in triangle.embeddings():
        doomed.append( emb.simplex() )
        verts.append( emb.vertices() )
    newIndex = [] # Value of newIndex[ doomed[i].index() ] is meaningless.
    doomedIndices = sorted( [ d.index() for d in doomed ] )
    for k in range( tri.size() ):
        if k < doomedIndices[0]:
            newIndex.append(k)
        elif k < doomedIndices[1]:
            newIndex.append(k-1)
        else:
            newIndex.append(k-2)
    survive = tri.size() - 2
    gluings = []
    triGlu = doomed[0].adjacentGluing( verts[0][3] )
    for j in range(3):
        jj = (j+1)%3

        # Glue the new tetrahedron that will eventually be numbered survive+j
        # to the new tetrahedron that will eventually be numbered survive+jj.
        gluings.append(
                ( survive+j, verts[0][jj], survive+jj,
                    Perm4( verts[0][j], verts[0][jj] ) ) )

        # Glue triangle verts[0][j] of survive+j.
        oldAdj = doomed[0].adjacentSimplex( verts[0][j] )
        if oldAdj is not None:
            oldGlu = doomed[0].adjacentGluing( verts[0][j] )
            try:
                i = doomed.index(oldAdj)
            except ValueError:
                # We are gluing to a tetrahedron that survives.
                newAdj = newIndex[ oldAdj.index() ]
                gluings.append(
                        ( survive+j, verts[0][j], newAdj, oldGlu ) )
            else:
                # We are gluing to a doomed tetrahedron, which means that we
                # are actually supposed to glue to one of the new tetrahedra.
                jjj = verts[i].inverse()[ oldGlu[ verts[0][j] ] ]

                # Do the gluing if we haven't already done it from the other
                # side.
                if jjj >= j:
                    if i == 0:
                        gluings.append(
                                ( survive+j, verts[0][j],
                                    survive+jjj, oldGlu ) )
                    else:
                        newGlu = ( triGlu.inverse() *
                                Perm4( verts[1][3], verts[1][jjj] ) *
                                oldGlu )
                        gluings.append(
                                ( survive+j, verts[0][j],
                                    survive+jjj, newGlu ) )

        # Glue triangle verts[0][3] of survive+j.
        oldAdj = doomed[1].adjacentSimplex( verts[1][j] )
        if oldAdj is not None:
            oldGlu = doomed[1].adjacentGluing( verts[1][j] )
            try:
                i = doomed.index(oldAdj)
            except ValueError:
                # We are gluing to a tetrahedron that survives.
                newAdj = newIndex[ oldAdj.index() ]
                newGlu = oldGlu * Perm4( verts[1][3], verts[1][j] ) * triGlu
                gluings.append(
                        ( survive+j, verts[0][3], newAdj, newGlu ) )
            else:
                # We are gluing to a doomed tetrahedron, which means that we
                # are actually supposed to glue to one of the new tetrahedra.
                jjj = verts[i].inverse()[ oldGlu[ verts[1][j] ] ]

                # Do the gluing if we haven't already done it from the other
                # side.
                if jjj > j:
                    if i == 0:
                        newGlu = ( oldGlu *
                                Perm4( verts[1][3], verts[1][j] ) *
                                triGlu )
                        gluings.append(
                                ( survive+j, verts[0][3],
                                    survive+jjj, newGlu ) )
                    else:
                        newGlu = ( triGlu.inverse() *
                                Perm4( verts[1][3], verts[1][jjj] ) *
                                oldGlu *
                                Perm4( verts[1][3], verts[1][j] ) *
                                triGlu )
                        gluings.append(
                                ( survive+j, verts[0][3],
                                    survive+jjj, newGlu ) )

    # For each edge e in tri, find one tetrahedron that will meet e after we
    # have performed the requested 2-3 move.
    edgeLocations = []
    for e in tri.edges():
        oldInd = e.index()
        emb = e.embedding(0)
        oldTet = emb.simplex()
        oldNum = emb.face()
        try:
            i = doomed.index(oldTet)
        except ValueError:
            # The tetrahedron oldTet survives.
            edgeLocations.append(
                    ( oldInd, newIndex[ oldTet.index() ], oldNum ) )
        else:
            # The tetrahedron oldTet is doomed, but this means that e will
            # meet one of the new tetrahedra.
            p = verts[i].inverse() * emb.vertices()
            # By construction, verts[i][ p[j] ] == emb.vertices()[j].

            apex = p.inverse()[3]
            if apex in {0,1}:
                j = p[ 1 - apex ]
                # We have e == oldTet.edge( verts[i][j], verts[i][3] ).
                ii = (j+1) % 3

                # Construct a permutation pp so that e will meet tetrahedron
                # survive+ii in the edge numbered Face3_1.faceNumber(pp).
                if i == 0:
                    pp = verts[0] * p
                else:
                    pp = verts[0] * Perm4( ii, 3 ) * p
            else:
                # We have e == oldTet.edge( verts[i][p[0]], verts[i][p[1]] ).
                ii = p[ 5 - apex ]
                
                # Construct a permutation pp so that e will meet tetrahedron
                # survive+ii in the edge numbered Face3_1.faceNumber(pp).
                pp = verts[0] * p
            edgeLocations.append(
                    ( oldInd, survive+ii, Face3_1.faceNumber(pp) ) )

    # Remove the two doomed tetrahedra, add three new tetrahedra, and then
    # perform the gluings that we just computed.
    for d in doomed:
        tri.removeTetrahedron(d)
    tri.newTetrahedra(3)
    for me, face, you, glu in gluings:
        tri.tetrahedron(me).join( face, tri.tetrahedron(you), glu )

    # How did the edges get renumbered?
    renum = {}
    for k, m, n in edgeLocations:
        renum[k] = tri.tetrahedron(m).edge(n).index()

    # Don't forget that we created a new edge!
    renum[-1] = tri.tetrahedron(survive).edge(
            verts[0][0], verts[0][3] ).index()
    return renum


def threeTwo(edge):
    """
    Performs a 3-2 move about the given edge, and returns a dictionary r that
    describes how the edges were renumbered.

    This routine directly modifies the triangulation that contains the given
    edge. If an edge is currently numbered i in the triangulation, then it
    will be numbered r[i] in the triangulation after the requested 3-2 move
    has been performed (provided the move is actually legal); also, a 3-2
    move removes the given edge, and this will be indicated by the fact that
    r[i] == -1, where i is the index of the given edge. If the move is not
    legal, the triangulation is left untouched and this routine returns None.

    If the triangulation containing the given edge is currently oriented,
    then this orientation will be preserved by the requested 3-2 move.
    """
    tri = edge.triangulation()

    # Is the requested 3-2 move legal?
    if not tri.pachner(edge, True, False):
        return None

    # We need to work out the gluings that we need to perform before we make
    # any changes.
    doomed = []
    verts = []
    for emb in edge.embeddings():
        doomed.append( emb.simplex() )
        verts.append( emb.vertices() )
    newIndex = [] # Value of newIndex[ doomed[i].index() ] is meaningless.
    doomedIndices = sorted( [ d.index() for d in doomed ] )
    for k in range( tri.size() ):
        if k < doomedIndices[0]:
            newIndex.append(k)
        elif k < doomedIndices[1]:
            newIndex.append(k-1)
        elif k < doomedIndices[2]:
            newIndex.append(k-2)
        else:
            newIndex.append(k-3)
    survive = tri.size() - 3
    gluings = [ ( survive, verts[1][1], survive+1,
        Perm4( verts[1][0], verts[1][1] ) ) ]
    # The two new tetrahedra will end up being numbered survive+j, for j in
    # {0,1}. These tetrahedra will be labelled in such a way that triangle
    # verts[1][j] of survive+j corresponds precisely to triangle verts[1][j]
    # of doomed[1], with precisely the same vertex numberings. Since this
    # doesn't tell us how doomed[1] meets the other two doomed tetrahedra,
    # we record this information now in the variable triGlu.
    triGlu = [
            doomed[1].adjacentGluing( verts[1][3] ),  # Gluing to doomed[0]
            None,
            doomed[1].adjacentGluing( verts[1][2] ) ] # Gluing to doomed[2]
    for j in range(2):
        # We have already glued triangle verts[1][1-j] of survive+j (to
        # triangle verts[1][j] of survive+(1-j)). We now glue the other
        # triangles of survive+j.

        # Glue triangle verts[1][j] of survive+j.
        oldAdj = doomed[1].adjacentSimplex( verts[1][j] )
        if oldAdj is not None:
            oldGlu = doomed[1].adjacentGluing( verts[1][j] )
            try:
                i = doomed.index(oldAdj)
            except ValueError:
                # We are gluing to a tetrahedron that survives.
                newAdj = newIndex[ oldAdj.index() ]
                gluings.append(
                        ( survive+j, verts[1][j], newAdj, oldGlu ) )
            else:
                # We are gluing to a doomed tetrahedron, which means that we
                # are actually supposed to glue to one of the new tetrahedra.
                jj = verts[i].inverse()[ oldGlu[ verts[1][j] ] ]

                # Do the gluing if we haven't already done it from the other
                # side.
                if jj >= j:
                    if i == 1:
                        # In this case, we know that jj == 1 and j == 0.
                        gluings.append(
                                ( survive, verts[1][0], survive+1, oldGlu ) )
                    else:
                        newGlu = ( triGlu[i].inverse() *
                                Perm4( verts[i][2+i//2], verts[i][jj] ) *
                                oldGlu )
                        gluings.append(
                                ( survive+j, verts[1][j],
                                    survive+jj, newGlu ) )

        # Glue triangle verts[1][2] of survive+j.
        oldAdj = doomed[2].adjacentSimplex( verts[2][j] )
        if oldAdj is not None:
            oldGlu = doomed[2].adjacentGluing( verts[2][j] )
            try:
                i = doomed.index(oldAdj)
            except ValueError:
                # We are gluing to a tetrahedron that survives.
                newAdj = newIndex[ oldAdj.index() ]
                newGlu = ( oldGlu *
                        Perm4( verts[2][3], verts[2][j] ) *
                        triGlu[2] )
                gluings.append(
                        ( survive+j, verts[1][2], newAdj, newGlu ) )
            else:
                # We are gluing to a doomed tetrahedron, which means that we
                # are actually supposed to glue to one of the new tetrahedra.
                jj = verts[i].inverse()[ oldGlu[ verts[2][j] ] ]

                # Do the gluing if we haven't already done it from the other
                # side.
                if jj >= j:
                    if i == 1:
                        # In this case, if jj == j then we would already have
                        # glued from the other side.
                        if jj > j:
                            # In this case, we know that jj == 1 and j == 0.
                            newGlu = ( oldGlu *
                                    Perm4( verts[2][3], verts[2][0] ) *
                                    triGlu[2] )
                            gluings.append(
                                    ( survive, verts[1][2],
                                        survive+1, newGlu ) )
                    else:
                        newGlu = ( triGlu[i].inverse() *
                                Perm4( verts[i][2+i//2], verts[i][jj] ) *
                                oldGlu *
                                Perm4( verts[2][3], verts[2][j] ) *
                                triGlu[2] )
                        gluings.append(
                                ( survive+j, verts[1][2],
                                    survive+jj, newGlu ) )

        # Glue triangle verts[1][3] of survive+j.
        oldAdj = doomed[0].adjacentSimplex( verts[0][j] )
        if oldAdj is not None:
            oldGlu = doomed[0].adjacentGluing( verts[0][j] )
            try:
                i = doomed.index(oldAdj)
            except ValueError:
                # We are gluing to a tetrahedron that survives.
                newAdj = newIndex[ oldAdj.index() ]
                newGlu = ( oldGlu *
                        Perm4( verts[0][2], verts[0][j] ) *
                        triGlu[0] )
                gluings.append(
                        ( survive+j, verts[1][3], newAdj, newGlu ) )
            else:
                # We are gluing to a doomed tetrahedron, which means that we
                # are actually supposed to glue to one of the new tetrahedra.
                jj = verts[i].inverse()[ oldGlu[ verts[0][j] ] ]

                # Do the gluing if we haven't already done it from the other
                # side.
                if jj > j:
                    # We know that jj == 1 and j == 0.
                    if i == 1:
                        newGlu = ( oldGlu *
                                Perm4( verts[0][2], verts[0][0] ) *
                                triGlu[0] )
                        gluings.append(
                                ( survive, verts[1][3], survive+1, newGlu ) )
                    else:
                        newGlu = ( triGlu[i].inverse() *
                                Perm4( verts[i][2+i//2], verts[i][1] ) *
                                oldGlu *
                                Perm4( verts[0][2], verts[0][0] ) *
                                triGlu[0] )
                        gluings.append(
                                ( survive, verts[1][3], survive+1, newGlu ) )

    # For each edge e in tri, find one tetrahedron that will meet e after we
    # have performed the requested 3-2 move.
    edgeLocations = []
    for e in tri.edges():
        # The edge about which we perform the 3-2 move gets removed entirely.
        if e == edge:
            continue
        oldInd = e.index()
        emb = e.embedding(0)
        oldTet = emb.simplex()
        oldNum = emb.face()
        try:
            i = doomed.index(oldTet)
        except ValueError:
            # The tetrahedron oldTet survives.
            edgeLocations.append(
                    ( oldInd, newIndex[ oldTet.index() ], oldNum ) )
        else:
            # The tetrahedron oldTet is doomed, but this means that e will
            # meet one of the new tetrahedra.
            p = verts[i].inverse() * emb.vertices()
            # By construction, verts[i][ p[j] ] == emb.vertices()[j].

            top = p.inverse()[0]
            bot = p.inverse()[1]
            # Note that at most one of top or bot is in {0,1} (otherwise we
            # would have e == edge).
            if top in {0,1}:
                j = p[ 1 - top ]
                # We have e == oldTet.edge( verts[i][j], verts[i][0] ).
                ii = 1

                # Construct a permutation pp so that e will meet tetrahedron
                # survive+ii in the edge numbered Face3_1.faceNumber(pp).
                swap = { 2: 1, 0: 2, 1: 3 }
                pp = verts[1] * Perm4( j, swap[(i+j)%3] ) * p
            elif bot in {0,1}:
                j = p[ 1 - bot ]
                # We have e == oldTet.edge( verts[i][j], verts[i][1] ).
                ii = 0

                # Construct a permutation pp so that e will meet tetrahedron
                # survive+ii in the edge numbered Face3_1.faceNumber(pp).
                swap = { 2: 0, 0: 2, 1: 3 }
                pp = verts[1] * Perm4( j, swap[(i+j)%3] ) * p
            else:
                # We have e == oldTet.edge( verts[i][2], verts[i][3] ).
                ii = 1

                # Construct a permutation pp so that e will meet tetrahedron
                # survive+ii in the edge numbered Face3_1.faceNumber(pp).
                r = Perm3.rot(i)
                pp = verts[1] * Perm4( r[0]+1, r[1]+1, r[2]+1, 0 )
            edgeLocations.append(
                    ( oldInd, survive+ii, Face3_1.faceNumber(pp) ) )

    # Remove the three doomed tetrahedra, add two new tetrahedra, and then
    # perform the gluings that we just computed.
    removeIndex = edge.index() # Need to remember this for later.
    for d in doomed:
        tri.removeTetrahedron(d)
    tri.newTetrahedra(2)
    for me, face, you, glu in gluings:
        tri.tetrahedron(me).join( face, tri.tetrahedron(you), glu )

    # How did the edges get renumbered?
    renum = {}
    for k, m, n in edgeLocations:
        renum[k] = tri.tetrahedron(m).edge(n).index()

    # Don't forget that we removed the given edge!
    renum[ removeIndex ] = -1
    return renum


def twoZero(edge):
    """
    Performs a 2-0 move about the given edge, and returns a dictionary r that
    describes how the edges were renumbered.

    This routine directly modifies the triangulation that contains the given
    edge. If an edge is currently numbered i in the triangulation, then it
    will be numbered r[i] in the triangulation after the requested 2-0 move
    has been performed (provided the move is actually legal); also, a 2-0
    move always removes the given edge, and this will be indicated by the
    fact that r[i] == -1, where i is the index of the given edge. If the move
    is not legal, the triangulation is left untouched and this routine
    returns None.

    If the triangulation containing the given edge is currently oriented,
    then this orientation will be preserved by the requested 2-0 move.
    """
    tri = edge.triangulation()

    # Is the requested 2-0 move legal?
    if not tri.twoZeroMove(edge, True, False):
        return None

    # How will the tetrahedra in tri get renumbered after we perform the
    # requested 2-0 move?
    doomed = sorted(
            [ emb.simplex().index() for emb in edge.embeddings() ] )
    newIndex = [] # Value of newIndex[ doomed[i] ] is meaningless.
    for k in range( tri.size() ):
        if k < doomed[0]:
            newIndex.append(k)
        elif k < doomed[1]:
            newIndex.append(k-1)
        else:
            newIndex.append(k-2)

    # For each edge e in tri, find one tetrahedron that will meet e after we
    # have performed the requested 2-0 move.
    # Note that the two edges "opposite" e become identified, so we need to
    # treat these edges a bit differently from the others.
    opp = { emb.simplex().edge( 5 - emb.face() )
            for emb in edge.embeddings() }
    oppDest = None
    edgeLocations = []
    for e in tri.edges():
        # The edge about which we perform the 3-2 move gets removed entirely.
        if e == edge:
            continue
        # If e is one of the opposite edges and we have already processed the
        # other opposite edge, then we already know where e will end up.
        if ( e in opp ) and ( oppDest is not None ):
            edgeLocations.append(
                    ( e.index(), oppDest[0], oppDest[1] ) )
            continue
        for emb in e.embeddings():
            oldTet = emb.simplex().index()
            if oldTet not in doomed:
                edgeLocations.append(
                        ( e.index(), newIndex[oldTet], emb.face() ) )
                # If e is one of the opposite edges (and we haven't yet
                # processed the other opposite edge), then remember where e
                # ends up (because the other opposite edge will end up in the
                # same place).
                if e in opp:
                    oppDest = ( newIndex[oldTet], emb.face() )
                # Move on to the next edge.
                break
        else:
            # We did not break out of this for loop early. The only way this
            # can happen is if e is one of the opposite edges, in which case
            # what we should do is find a tetrahedron that will meet the
            # *other* opposite edge.
            for otherOpp in opp:
                if otherOpp != e:
                    break
            for otherEmb in otherOpp.embeddings():
                otherTet = otherEmb.simplex().index()
                if otherTet not in doomed:
                    oppDest = ( newIndex[otherTet], otherEmb.face() )
                    edgeLocations.append(
                            ( e.index(), oppDest[0], oppDest[1] ) )
                    # Move on to the next edge.
                    break

    # Perform the 2-0 move, and work out how the edges were renumbered.
    removeIndex = edge.index() # Need to remember this for later.
    tri.twoZeroMove( edge, False, True )
    renum = {}
    for k, m, n in edgeLocations:
        renum[k] = tri.tetrahedron(m).edge(n).index()

    # Don't forget that we removed the given edge!
    renum[ removeIndex ] = -1
    return renum


def fourFour( edge, newAxis ):
    """
    Performs a 4-4 move about the given edge, and returns a dictionary r that
    describes how the edges were renumbered.

    If the 4-4 move is legal, then the given edge forms the axis of a
    four-tetrahedron octahedron. There are two other choices of axis for this
    octahedron, and the requested 4-4 move retriangulates this octahedron to
    use one of these other choices of axis edge. Specifically, number the
    original four tetrahedra 0, 1, 2 and 3 according to the order described
    by e.embeddings(). If newAxis is 0, then the new axis will separate
    tetrahedra 0 and 1 from tetrahedra 2 and 3. If newAxis is 1, then the new
    axis will separate tetrahedra 1 and 2 from tetrahedra 3 and 0.

    This routine directly modifies the triangulation that contains the given
    edge. If an edge is currently numbered i in the triangulation, then it
    will be numbered r[i] in the triangulation after the requested 4-4 move
    has been performed (provided the move is actually legal); also, we will
    have r[i] == -1, where i is the index of the given edge (since this edge
    gets removed), and we will have r[-1] == j, where j is the index of the
    newly-created axis edge. If the move is not legal, the triangulation is
    left untouched and this routine returns None.

    If the triangulation containing the given edge is currently oriented,
    then this orientation will be preserved by the requested 4-4 move.
    """
    tri = edge.triangulation()

    # Is the requested 4-4 move legal?
    if not tri.fourFourMove(edge, newAxis, True, False):
        return None

    # Find the doomed tetrahedra.
    doomed = []
    verts = []
    for emb in edge.embeddings():
        doomed.append( emb.simplex() )
        verts.append( emb.vertices() )

    # Perform the 4-4 move as a 2-3 move followed by a 3-2 move.
    if newAxis == 0:
        f23 = doomed[0].triangle( verts[0][2] )
    else:
        f23 = doomed[1].triangle( verts[1][2] )
    e32 = edge.embedding(3).edge()
    rInc = twoThree(f23)
    rDec = threeTwo( doomed[3].edge(e32) )

    # Work out how the edges were renumbered. Don't forget that we have
    # removed an edge, and also created an edge.
    renum = {}
    for e in range( -1, tri.countEdges() ):
        renum[e] = rDec[ rInc[e] ]
    return renum


def twoOne( edge, edgeEnd ):
    #TODO Write documentation.
    """
    Performs a 2-1 move about the given edge, and returns a dictionary r that
    describes how the edges were renumbered.

    If the 2-1 move is legal, then the given edge forms the internal edge of
    a snapped ball B, and the move involves the following two tetrahedra:
    --> the tetrahedron T that forms B; and
    --> the tetrahedron that is glued to T along the triangle opposite the
        vertex number edgeEnd of the given edge.

    This routine directly modifies the triangulation that contains the given
    edge. If an edge is currently numbered i in the triangulation, then it
    will be numbered r[i] in the triangulation after the requested 2-1 move
    has been performed (provided the move is actually legal); also, we will
    have r[i] == -1, where i is the index of the given edge (since this edge
    gets removed), and we will have r[-1] == j, where j is the index of the
    newly-created edge. If the move is not legal, the triangulation is left
    untouched and this routine returns None.

    If the triangulation containing the given edge is currently oriented,
    then this orientation will be preserved by the requested 2-1 move.
    """
    tri = edge.triangulation()

    # Is the requested 2-1 move legal?
    if not tri.twoOneMove( edge, edgeEnd, True, False ):
        return None

    # Find the doomed tetrahedra.

    # Perform the 2-1 move as a 2-3 move followed by a 2-0 move.
    emb = edge.embedding(0)
    f23 = emb.tetrahedron().triangle( emb.vertices()[edgeEnd] )
    e20 = edge.index()
    rInc = twoThree(f23)
    rDec = twoZero( tri.edge( rInc[e20] ) )

    # Work out how the edges were renumbered. Don't forget that we have
    # removed an edge, and also created an edge.
    renum = {}
    for e in range( -1, tri.countEdges() + 1 ):
        renum[e] = rDec[ rInc[e] ]
    return renum


# Tests.
if __name__ == "__main__":
    print()
    print( "TESTS FOR ELEMENTARY MOVES" )
    print( "==========================" )
    print()

    t = Triangulation3("gLLPQcdefeffpvauppb")
    t.orient()
    RandomEngine.reseedWithHardware()
    iso = Isomorphism3.random(t.size(),True)
    iso.applyInPlace(t)

    # Test 2-3 and 3-2 moves.
    print( "2-3 and 3-2 moves on \"gLLPQcdefeffpvauppb\"" )
    print()
    for f in t.triangles():
        # Test that twoThree gives correct isomorphism type, and that it
        # preserves orientation.
        pach = Triangulation3(t)
        t23 = Triangulation3(t)
        pach.pachner( pach.triangle( f.index() ) )
        trenum = twoThree( t23.triangle( f.index() ) )
        if not t23.isOriented():
            print("{}: Not oriented!".format(f.index()))
            print(t23)
            break
        if not pach.isIsomorphicTo(t23):
            print("{}: Oriented, but not isomorphic!".format(f.index()))
            print(t)
            print(pach)
            print(t23)
            break

        # Test that threeTwo inverts twoThree correctly, and that it
        # preserves orientation.
        # This test makes assumptions about the implementation of twoThree.
        tinv = Triangulation3(t23)
        tinnum = threeTwo( tinv.tetrahedron( tinv.size()-1 ).edge(
            f.embedding(0).vertices()[2], f.embedding(0).vertices()[3] ) )
        if not tinv.isOriented():
            print("{}: Inverse not oriented!".format(f.index()))
            print(tinv)
            break
        if not t.isIsomorphicTo(tinv):
            print("{}: Inverse oriented, but not isomorphic!".format(
                f.index() ))
            print(t)
            print(tinv)
            print(tinv.detail())
            break

        # Now check that twoThree gives correct isomorphism type after a
        # random relabelling.
        RandomEngine.reseedWithHardware()
        r = Triangulation3(t)
        iso = Isomorphism3.random(r.size())
        iso.applyInPlace(r)
        source = f.embedding(0).simplex().index()
        num = f.embedding(0).face()
        simpImage = iso.simpImage(source)
        faceImage = iso.facetPerm(source)[num]
        r23 = Triangulation3(r)
        rrenum = twoThree(
                r23.tetrahedron(simpImage).triangle(faceImage) )
        if not pach.isIsomorphicTo(r23):
            print("{}: Not isomorphic!".format(f.index()))
            print(r)
            print(pach)
            print(r23)
            break

        # Finally, check that the renumberings are sensible by comparing edge
        # degrees in the isomorphic triangulations t and tinv.
        error = ""
        for i in range( t.countEdges() ):
            deg = t.edge(i).degree()
            comDeg = tinv.edge( tinnum[ trenum[i] ] ).degree()
            if deg != comDeg:
                error = "<--"
                break
        print( "{}: {}, {}; {} {}".format(
            f.index(), trenum, rrenum, tinnum, error ) )
    else:
        print()
        print( "2-3 and 3-2 moves: Success!" )

    print()
    print("========================")
    print()
    print( "0-2 and 2-0 moves on \"gLLPQcdefeffpvauppb\"" )
    print()

    # Test 2-0 moves.
    for e in t.edges():
        deg = e.degree()
        for i in range(deg):
            RandomEngine.reseedWithHardware()
            ii = i + RandomEngine.rand( deg - i )
            tt = Triangulation3(t)
            if not tt.zeroTwoMove( tt.edge( e.index() ), i, ii ):
                continue

            # Test the inverse 2-0 move using twoZero().
            RandomEngine.reseedWithHardware()
            iso = Isomorphism3.random(tt.size())
            iso.applyInPlace(tt)
            simpImage = iso.simpImage( tt.size() - 1 )
            facetPerm = iso.facetPerm( tt.size() - 1 )
            testInd = tt.tetrahedron( simpImage ).edge(
                    facetPerm[2], facetPerm[3] ).index()
            renum = twoZero( tt.edge(testInd) )
            if not t.isIsomorphicTo(tt):
                print( "{}/{}/{}".format( e.index(), i, ii ) +
                        ": 0-2 followed by 2-0 not isomorphic!" )

            # Very basic sanity check for the renumberings.
            bef = len( set( renum.keys() ) )
            aft = len( set( renum.values() ) )
            if ( renum[testInd] != -1 ) or ( bef != 1 + aft ):
                error = "<--"
            else:
                error = ""
            print( "{}/{}/{}: {} {}".format(
                e.index(), i, ii, renum, error ) )
    else:
        print()
        print( "0-2 and 2-0 moves: Success!" )

    print()
    print("========================")
    print()
    print( "4-4 moves on \"gLLPQcdefeffpvauppb\"" )
    print()

    # Test 4-4 moves.
    for e in t.edges():
        for newAxis in range(2):
            reg44 = Triangulation3(t)
            new44 = Triangulation3(t)
            if not reg44.fourFourMove( reg44.edge( e.index() ), newAxis ):
                continue
            renum = fourFour( new44.edge( e.index() ), newAxis )

            # Test that fourFour gives the right isomorphism type, and that
            # it outputs sensible renumberings.
            if not reg44.isIsomorphicTo(new44):
                print("{}/{}: 4-4 not isomorphic!".format(
                    e.index(), newAxis ))
                print(t)
                print(reg44)
                print(new44)
                break
            if set( renum.keys() ) != set( renum.values() ):
                error = "<--"
            else:
                error = ""
            print( "{}/{}: {} {}".format(
                e.index(), newAxis, renum, error ) )
    else:
        print()
        print( "4-4 moves: Success!" )

    # Perform a bunch of moves to check that we don't get any exceptions.
    print()
    print("========================")
    print()
    print( "Perform moves on \"cMcabbgqs\" up to height 3." )
    print()
    initSig = "cMcabbgqs"
    stack = [ initSig ]
    seen = { initSig }
    initSize = 2
    maxSize = 5
    counts = { "2-1": 0, "2-0": 0, "3-2": 0, "4-4": 0, "2-3": 0 }
    while stack:
        sig = stack.pop()
        tri = Triangulation3.fromIsoSig(sig)
        if tri.size() <= initSize:
            print(sig)
            stdout.flush()

        # Moves on edges.
        for e in tri.edges():
            if e.degree() == 1:
                for edgeEnd in {0,1}:
                    newTri = Triangulation3(tri)
                    renum = twoOne( newTri.edge( e.index() ), edgeEnd )
                    if renum is None:
                        continue
                    counts["2-1"] += 1

                    # If we haven't seen the new triangulation before (up to
                    # combinatorial isomorphism), then add it to the stack.
                    newSig = newTri.isoSig()
                    if newSig not in seen:
                        seen.add(newSig)
                        stack.append(newSig)
            elif e.degree() == 2:
                newTri = Triangulation3(tri)
                renum = twoZero( newTri.edge( e.index() ) )
                if renum is None:
                    continue
                counts["2-0"] += 1

                # If we haven't seen the new triangulation before (up to
                # combinatorial isomorphism), then add it to the stack.
                newSig = newTri.isoSig()
                if newSig not in seen:
                    seen.add(newSig)
                    stack.append(newSig)
            elif e.degree() == 3:
                newTri = Triangulation3(tri)
                renum = threeTwo( newTri.edge( e.index() ) )
                if renum is None:
                    continue
                counts["3-2"] += 1

                # If we haven't seen the new triangulation before (up to
                # combinatorial isomorphism), then add it to the stack.
                newSig = newTri.isoSig()
                if newSig not in seen:
                    seen.add(newSig)
                    stack.append(newSig)
            elif e.degree() == 4:
                for newAxis in {0,1}:
                    newTri = Triangulation3(tri)
                    renum = fourFour( newTri.edge( e.index() ), newAxis )
                    if renum is None:
                        continue
                    counts["4-4"] += 1

                    # If we haven't seen the new triangulation before (up to
                    # combinatorial isomorphism), then add it to the stack.
                    newSig = newTri.isoSig()
                    if newSig not in seen:
                        seen.add(newSig)
                        stack.append(newSig)

        # 2-3 moves (only if such moves would not exceed the given height).
        if tri.size() == maxSize:
            continue
        for f in tri.triangles():
            newTri = Triangulation3(tri)
            renum = twoThree( newTri.triangle( f.index() ) )
            if renum is None:
                continue
            counts["2-3"] += 1

            # If we haven't seen the new triangulation before (up to
            # combinatorial isomorphism), then add it to the stack.
            newSig = newTri.isoSig()
            if newSig not in seen:
                seen.add(newSig)
                stack.append(newSig)

    # Print move counts.
    print()
    print( "Tested: {} {}; {} {}; {} {}; {} {}; {} {}.".format(
        counts["2-1"], "2-1 moves",
        counts["2-0"], "2-0 moves",
        counts["3-2"], "3-2 moves",
        counts["4-4"], "4-4 moves",
        counts["2-3"], "2-3 moves" ) )
