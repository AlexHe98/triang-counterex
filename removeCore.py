"""
Greedily search for a 1-vertex triangulation of a lens space that has no
"bad" edges. Here, we consider an edge e to be bad if it forms a "core loop",
meaning that pinching e yields an ideal triangulation of the solid torus.

To do this, we start with a triangulation that has at least one bad edge; for
example, the triangulation given by the isomorphism signature "cMcabbgqs" is
a 1-vertex triangulation of the 3-sphere with three edges, exactly two of
which form core loops. We then search the Pachner graph for a triangulation
with one less bad edge.

The search prioritises triangulations in which:
(1) the maximum "multiplicity defect" among bad edges is small (see the
    mulDefect() routine);
(2) the maximum "degree defect" among bad edges is small (see the degDefect()
    routine); and
(3) the number of tetrahedra is small.
"""
from sys import argv, stdout
from triangCounterexHelpers import isCore
from removeBadEdge import removeBadEdge


if __name__ == "__main__":
    sig = argv[1]
    nProcesses = int( argv[2] )
    interval = int( argv[3] )
    print( "Greedily removing core edges. Sig: {}. #Processes: {}.".format(
            sig, nProcesses ) )
    stdout.flush()
    print( "Results: {}\n".format( removeBadEdge(
        sig, nProcesses, interval, isCore, "Core edges" ) ) )

