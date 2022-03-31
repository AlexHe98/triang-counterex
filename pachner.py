"""
Perform 2-3 and 3-2 moves while tracking how edges are relabelled.
"""
from regina import *


def twoThree(triangle):
    """
    Performs a 2-3 move about the given triangle, and returns both: the
    triangulation obtained from the move, and a dictionary that describes
    how the edges were renumbered.

    If the requested 2-3 move may be performed without changing the
    topology of the original triangulation that contains the given
    triangle, then returns a pair (t, d), where t is the new triangulation
    that is obtained after performing the move and d is the dictionary
    mentioned above (and described in more detail below). Otherwise, simply
    returns None. In either case, this routine does *not* modify the
    original triangulation.

    The returned dictionary d has keys corresponding to edge indices in the
    original triangulation and values corresponding to edge indices in the
    new triangulation, with the exception that the newly-created edge has
    key -1 (since this new edge has no corresponding edge in the original
    triangulation).

    Note that the implementation of this routine is deliberately simple,
    and does nothing to ensure that tetrahedra are labelled in a way that
    preserves orientation.

    We do, however, ensure the following:
    --> The last three tetrahedra in the returned triangulation are exactly
        the three tetrahedra that are created by the 2-3 move.
    --> In each of these three tetrahedra, the newly-created edge appears
        as edge(0,1).
    --> The order of the indices of these three tetrahedra respects the
        vertex numbering of the given triangle in the original
        triangulation.
    --> Moreover, walking from vertex 2 to vertex 3 of each of the new
        tetrahedra corresponds to walking from vertex i of the given
        triangle to vertex i+1 (mod 3).
    """
    oldTri = triangle.triangulation()

    # First, use Regina to check eligibility of the 2-3 move (but do *not*
    # perform the move, as this routine will use its own implementation of
    # performing a 2-3 move).
    if not oldTri.pachner(triangle, True, False):
        return None

    # Now we know that the 2-3 move is legal. We perform the 2-3 move on a
    # *copy* of the original triangulation.
    newTri = Triangulation3(oldTri)

    # Create the three new tetrahedra and glue their "internal faces".
    newTet = [newTri.newTetrahedron() for _ in range(3)]
    for i in range(3):
        # Join newTet[i](013) to newTet[i+1](012).
        newTet[i].join(
            2,
            newTet[(i+1) % 3],
            Perm4(*[0,1,3,2]) )

    # Calculate the "external gluings" for the three new tetrahedra.
    # If the gluing is between two faces of these new tetrahedra, then just
    # perform the gluing immediately.
    # Otherwise, if the gluing is between a face of a new tetrahedron and a
    # face of an old tetrahedron, then save the gluing to perform later.
    oldDestroy = [e.simplex().index() for e in triangle.embeddings()]
    newGluings = []
    for i, emb in enumerate(triangle.embeddings()):
        oldTet = emb.simplex()
        for j in range(3):
            # Continue to the next face if face (i,2,3) of newTet[j] has
            # already been glued from the "other side".
            newFaceNum = 1 - i
            if newTet[j].adjacentSimplex(newFaceNum) is not None:
                continue

            # The face in oldTri of oldTet with vertices
            # ( emb.face(),
            #   emb.vertices()[j],
            #   emb.vertices()[j+1] )
            # corresponds to the face in newTri of newTet[j] with vertices
            # (i,2,3).
            oldVertices = [
                emb.face(),
                emb.vertices()[j],
                emb.vertices()[(j+1) % 3] ]
            oldFaceNum = [k for k in range(4)
                          if k not in oldVertices].pop()
            oldGluingPerm = oldTet.adjacentGluing(oldFaceNum)
            newGluingList = []
            for k in range(4):
                if k == newFaceNum:
                    image = oldGluingPerm[oldFaceNum]
                elif k == i:
                    image = oldGluingPerm[emb.face()]
                else:
                    image = oldGluingPerm[
                        emb.vertices()[(j+k-2) % 3]]
                newGluingList.append(image)

            # The easy case is if the gluing of oldTet is to some
            # tetrahedron that will *not* be destroyed.
            oldGluingTetIndex = oldTet.adjacentSimplex(oldFaceNum).index()
            if oldGluingTetIndex not in oldDestroy:
                newGluings.append((
                    newTet[j],
                    newFaceNum,
                    newTri.tetrahedron(oldGluingTetIndex),
                    Perm4(*newGluingList)
                    ))
                continue

            # If we reach this stage, then the gluing of oldTet is to some
            # tetrahedron that *will* be destroyed.
            # In other words, face (i,2,3) of newTet[j] needs to be glued
            # (in some orientation determined by triangle.embedding(ii)) to
            # some face (ii,2,3) of some newTet[jj].
            # oldGluingTetIndex tells us what ii should be, and we can work
            # out jj from triangle.embedding(ii), together with
            # newGluingList[newFaceNum].
            ii = oldDestroy.index(oldGluingTetIndex)
            workingEmb = triangle.embedding(ii)
            for k in range(3):
                if workingEmb.vertices()[k] == newGluingList[newFaceNum]:
                    jj = (k+1) % 3
                    break

            # We now need to adjust the gluing permutation newGluingList.
            for k in range(4):
                oldImage = newGluingList[k]
                if k == newFaceNum:
                    newImage = 1 - ii
                elif oldImage == workingEmb.vertices()[jj]:
                    newImage = 2
                elif oldImage == workingEmb.vertices()[(jj+1) % 3]:
                    newImage = 3
                else:
                    newImage = ii
                newGluingList[k] = newImage

            # Immediately glue!
            newTet[j].join(
                newFaceNum,
                newTet[jj],
                Perm4(*newGluingList) )

    # Destroy the tetrahedra that need to be destroyed, and then perform
    # the remaining new gluings.
    newDestroy = [newTri.tetrahedron(i) for i in oldDestroy]
    for tet in newDestroy:
        newTri.removeTetrahedron(tet)
    for myTet, myFace, yourTet, gluing in newGluings:
        myTet.join(myFace, yourTet, gluing)

    # We have finished performing the 2-3 move, but we still need to
    # calculate the renumbering dictionary.
    renumbering = dict()
    for oldEdgeIndex in range(oldTri.countEdges()):
        oldEmb = oldTri.edge(oldEdgeIndex).embedding(0)
        oldTetIndex = oldEmb.simplex().index()

        # The easy case is if oldTetIndex does *not* correspond to one of
        # the tetrahedra that got destroyed as a result of the 2-3 move.
        if oldTetIndex not in oldDestroy:
            # Adjust the index to account for the removed tetrahedra.
            newTetIndex = oldTetIndex - len([
                i for i in oldDestroy if (i < oldTetIndex) ])

            # It is now straightforward to calculate the renumbering.
            renumbering[oldEdgeIndex] = newTri.tetrahedron(
                newTetIndex).edge(oldEmb.face()).index()
            continue

        # We are now in the slightly harder case: oldTetIndex got
        # destroyed, so we have to convert oldTetIndex and oldEmb.face()
        # into newTet numberings.
        # Thinking of the two destroyed tetrahedra as forming a bipyramid
        # based around the given triangle, we have two cases to consider:
        # - oldEmb.face() is an edge of the given triangle;
        # - oldEmb.face() is a "slope" of the bipyramid.
        for i, e in enumerate(triangle.embeddings()):
            if e.simplex().index() == oldTetIndex:
                triEmbIndex = i
                triEmb = e
                break
        oldEmbVerts = [oldEmb.vertices()[i] for i in range(2)]

        # Handle the "slope" case first.
        if triEmb.face() in oldEmbVerts:
            oldEmbVerts.remove(triEmb.face())
            oldOtherVert = oldEmbVerts.pop()

            # The "index" i of otherVert in triEmb.vertices() tells us that
            # oldOtherVert corresponds to vertex 2 of newTet[i].
            for i in range(3):
                if triEmb.vertices()[i] == oldOtherVert:
                    newTetIndex = i

            # We also know that triEmb.face() corresponds to vertex
            # triEmbIndex of newTet[i]. This suffices to define the
            # renumbering so we can continue to the next edge.
            renumbering[oldEdgeIndex] = newTet[newTetIndex].edge(
                2, triEmbIndex).index()
            continue

        # We now know that oldEmb.face() is an edge of the given triangle.
        # The vertex numbers are 2 and 3, so we just need to work out which
        # newTet[j] we have.
        for i in range(3):
            if triEmb.vertices()[i] not in oldEmbVerts:
                newTetIndex = (i+1) % 3
                break
        renumbering[oldEdgeIndex] = newTet[newTetIndex].edge(
            2, 3).index()

    # Don't forget to include the newly-created edge in the renumbering.
    renumbering[-1] = newTet[0].edge(0, 1).index()

    # All done!
    return (newTri, renumbering)


def threeTwo(edge):
    """
    Performs a 3-2 move about the given edge, and returns both: the
    triangulation obtained from the move, and a dictionary that describes
    how the edges were renumbered.

    If the requested 3-2 move may be performed without changing the
    topology of the original triangulation that contains the given edge,
    then returns a pair (t, d), where t is the new triangulation that is
    obtained after performing the move and d is the dictionary mentioned
    above (and described in more detail below). Otherwise, simply
    returns None. In either case, this routine does *not* modify the
    original triangulation.

    The returned dictionary d has keys corresponding to edge indices in the
    original triangulation and values corresponding to edge indices in the
    new triangulation, with the exception that the given edge maps to the
    value -1 (since the given edge will be removed, and hence has no
    corresponding edge in the new triangulation).

    Note that the implementation of this routine is deliberately simple,
    and does nothing to ensure that tetrahedra are labelled in a way that
    preserves orientation.

    We do, however, ensure the following:
    --> The last two tetrahedra in the returned triangulation are exactly
        the two tetrahedra that are created by the 3-2 move.
    --> In each of these two tetrahedra, the newly-created triangle is the
        the one with vertices numbered (0,1,2).
    --> The order of the indices of these two tetrahedra respects the
        vertex numbering of the given edge in the original triangulation.
    --> Moreover, for i in {0,1,2}, walking from vertex i to vertex
        i+1 (mod 3) of either of the newly-created tetrahedra corresponds
        to walking through edge.embedding(i).simplex() in the direction
        from edge.embedding(i-1 (mod 3)).simplex() to
        edge.embedding(i+1 (mod 3)).simplex().
    """
    oldTri = edge.triangulation()

    # First, use Regina to check eligibility of the 3-2 move (but do *not*
    # perform the move, as this routine will use its own implementation of
    # performs a 3-2 move).
    if not oldTri.pachner(edge, True, False):
        return None

    # Now we know that the 3-2 move is legal. We perform the 3-2 move on a
    # *copy* of the original triangulation.
    newTri = Triangulation3(oldTri)

    # Create the two new tetrahedra and glue their "internal face" (i.e.,
    # join newTet[0](012) to newTet[1](012)).
    newTet = [newTri.newTetrahedron() for _ in range(2)]
    newTet[0].join(
        3,
        newTet[1],
        Perm4(*[0,1,2,3]) )

    # Before we can calculate "external gluings" for the two new
    # tetrahedra, we need to associate each vertex number of the new
    # triangle with "apex numbers" of the old triangles emanating from the
    # given edge.
    embTet = [emb.simplex() for emb in edge.embeddings()]
    embVerts = []
    oppVerts = []
    apexFaceIndices = []
    for i, emb in enumerate(edge.embeddings()):
        embVerts.append( [emb.vertices()[j] for j in range(2)] )
        oppVerts.append(
            [j for j in range(4) if (j not in embVerts[i])] )
        apexFaceIndices.append(
            [embTet[i].triangle(j).index() for j in oppVerts[i][::-1]] )
    numbering = [
        set(apexFaceIndices[0]).difference(apexFaceIndices[1]).pop(),
        set(apexFaceIndices[0]).intersection(apexFaceIndices[1]).pop(),
        set(apexFaceIndices[1]).difference(apexFaceIndices[0]).pop() ]

    # Calculate the "external gluings" for the two new tetrahedra.
    # If the gluing is between two faces of these new tetrahedra, then just
    # perform the gluing immediately.
    # Otherwise, if the gluing is between a face of a new tetrahedron and a
    # face of an old tetrahedron, then save the gluing to perform later.
    oldDestroy = [tet.index() for tet in embTet]
    newGluings = []
    for i, emb in enumerate(edge.embeddings()):
        oldTet = emb.simplex()
        for j in range(2):
            # Continue to the next face if face (i,i+1,3) of newTet[j] has
            # already been glued from the "other side".
            newFaceNum = (i-1) % 3
            if newTet[j].adjacentSimplex(newFaceNum) is not None:
                continue

            # The face in oldTri of oldTet with vertices
            # ( emb.vertices()[j],
            #   oppVerts[i][0],
            #   oppVerts[i][1] )
            # corresponds to the face in newTri of newTet[j] with vertices
            # ( 3,
            #   numbering.index(apexFaceIndices[i][0]),
            #   numbering.index(apexFaceIndices[i][1]) ).
            oldVertices = [
                emb.vertices()[j],
                oppVerts[i][0],
                oppVerts[i][1] ]
            oldFaceNum = [k for k in range(4)
                          if k not in oldVertices].pop()
            oldGluingPerm = oldTet.adjacentGluing(oldFaceNum)
            newGluingList = []
            for k in range(4):
                if k == newFaceNum:
                    image = oldGluingPerm[oldFaceNum]
                elif k == 3:
                    image = oldGluingPerm[emb.vertices()[j]]
                elif k == numbering.index(apexFaceIndices[i][0]):
                    image = oldGluingPerm[oppVerts[i][0]]
                else:
                    image = oldGluingPerm[oppVerts[i][1]]
                newGluingList.append(image)

            # The easy case is if the gluing of oldTet is to some
            # tetrahedron that will *not* be destroyed.
            oldGluingTetIndex = oldTet.adjacentSimplex(oldFaceNum).index()
            if oldGluingTetIndex not in oldDestroy:
                newGluings.append((
                    newTet[j],
                    newFaceNum,
                    newTri.tetrahedron(oldGluingTetIndex),
                    Perm4(*newGluingList)
                    ))
                continue

            # If we reach this stage, then the gluing of oldTet is to some
            # tetrahedron that *will* be destroyed.
            # In other words, face (i,i+1,3) of newTet[j] needs to be glued
            # (in some orientation determined by edge.embedding(ii)) to
            # some face (ii,ii+1,3) of some newTet[jj].
            # oldGluingTetIndex tells us what ii should be, and we can work
            # out jj from edge.embedding(ii), together with
            # newGluingList[newFaceNum].
            ii = oldDestroy.index(oldGluingTetIndex)
            workingEmb = edge.embedding(ii)
            for k in range(2):
                if workingEmb.vertices()[k] == newGluingList[newFaceNum]:
                    jj = 1 - k
                    break

            # We now need to "adjust" the gluing permutation newGluingList.
            for k in range(4):
                oldImage = newGluingList[k]
                if k == newFaceNum:
                    newImage = (ii-1) % 3
                elif oldImage == oppVerts[ii][0]:
                    newImage = numbering.index(apexFaceIndices[ii][0])
                elif oldImage == oppVerts[ii][1]:
                    newImage = numbering.index(apexFaceIndices[ii][1])
                else:
                    newImage = 3
                newGluingList[k] = newImage

            # Immediately glue!
            newTet[j].join(
                newFaceNum,
                newTet[jj],
                Perm4(*newGluingList) )

    # Destroy the tetrahedra that need to be destroyed, and then perform
    # the remaining new gluings.
    newDestroy = [newTri.tetrahedron(i) for i in oldDestroy]
    for tet in newDestroy:
        newTri.removeTetrahedron(tet)
    for myTet, myFace, yourTet, gluing in newGluings:
        myTet.join(myFace, yourTet, gluing)

    # We have finished performing the 3-2 move, but we still need to
    # calculate the renumbering dictionary.
    renumbering = dict()
    for oldEdgeIndex in range(oldTri.countEdges()):
        # Get the trivial "given edge" case out of the way.
        if oldEdgeIndex == edge.index():
            renumbering[oldEdgeIndex] = -1
            continue

        # We now know that oldEdgeIndex does *not* correspond to the given
        # edge.
        oldEmb = oldTri.edge(oldEdgeIndex).embedding(0)
        oldTetIndex = oldEmb.simplex().index()

        # The easy case is if oldTetIndex does *not* correspond to one of
        # the tetrahedra that got destroyed as a result of the 3-2 move.
        if oldTetIndex not in oldDestroy:
            # Adjust the index to account for the removed tetrahedra.
            newTetIndex = oldTetIndex - len([
                i for i in oldDestroy if (i < oldTetIndex) ])

            # It is now straightforward to calculate the renumbering.
            renumbering[oldEdgeIndex] = newTri.tetrahedron(
                newTetIndex).edge(oldEmb.face()).index()
            continue

        # We are now in the slightly harder case: oldTetIndex got
        # destroyed, so we have to convert oldTetIndex and oldEmb.face()
        # into newTet numberings.
        # Thinking of the three destroyed tetrahedra as forming a bipyramid
        # based around the newly-created triangle, we have two cases to
        # consider:
        # - oldEmb.face() is an edge of the newly-created triangle;
        # - oldEmb.face() is a "slope" of the bipyramid.
        for i, e in enumerate(edge.embeddings()):
            if e.simplex().index() == oldTetIndex:
                edgeEmbIndex = i
                edgeEmb = e
                break
        oldEmbVerts = [oldEmb.vertices()[i] for i in range(2)]

        # Handle the "slope" case first.
        if edgeEmb.vertices()[0] in oldEmbVerts:
            oldEmbVerts.remove(edgeEmb.vertices()[0])
            oldOtherVert = oldEmbVerts.pop()

            # The presence of edgeEmb.vertices()[0] tells us that we are in
            # edge (j, 3) of newTet[0], for some j in {0,1,2}. To work out
            # j, examine the oldOtherVert.
            for i, o in enumerate(oppVerts[edgeEmbIndex]):
                if o == oldOtherVert:
                    apexFaceIndex = apexFaceIndices[edgeEmbIndex][i]
                    break
            newOtherVert = numbering.index(apexFaceIndex)

            # This suffices to define the renumbering so we can continue to
            # the next edge.
            renumbering[oldEdgeIndex] = newTet[0].edge(
                newOtherVert, 3).index()
            continue
        elif edgeEmb.vertices()[1] in oldEmbVerts:
            oldEmbVerts.remove(edgeEmb.vertices()[1])
            oldOtherVert = oldEmbVerts.pop()

            # The presence of edgeEmb.vertices()[1] tells us that we are in
            # edge (j, 3) of newTet[1], for some j in {0,1,2}. To work out
            # j, examine the oldOtherVert.
            for i, o in enumerate(oppVerts[edgeEmbIndex]):
                if o == oldOtherVert:
                    apexFaceIndex = apexFaceIndices[edgeEmbIndex][i]
                    break
            newOtherVert = numbering.index(apexFaceIndex)

            # This suffices to define the renumbering so we can continue to
            # the next edge.
            renumbering[oldEdgeIndex] = newTet[1].edge(
                newOtherVert, 3).index()
            continue

        # We now know that oldEmb.face() is an edge of the newly-created
        # triangle. We know that is edge is exactly the edge of newTet[0]
        # with vertices (edgeEmbIndex, edgeEmbIndex + 1).
        renumbering[oldEdgeIndex] = newTet[0].edge(
            edgeEmbIndex, (edgeEmbIndex+1) % 3).index()

    # All done!
    return (newTri, renumbering)
