"""
Naively search for a 1-vertex triangulation of a lens space that has no "bad"
edges. Here, we consider an edge e to be bad if it forms a "core loop",
meaning that pinching e yields an ideal triangulation of the solid torus.

To do this, we start with a triangulation that has at least one bad edge; for
example, the triangulation given by the isomorphism signature "cMcabbgqs" is
a 1-vertex triangulation of the 3-sphere with three edges, exactly two of
which form core loops. We then search the Pachner graph for a triangulation
with one less bad edge.

The search naively prioritises triangulations with a small number of
tetrahedra. The main purpose is to demonstrate the effectiveness of the more
sophisticated heuristic used in removeCore.py.
"""
from sys import argv, stdout
from triangCounterexHelpers import isCore
from naiveBadEdge import naiveBadEdge


if __name__ == "__main__":
    sig = argv[1]
    nProcesses = int( argv[2] )
    interval = int( argv[3] )
    print( "Naively removing core edges. Sig: {}. #Processes: {}.".format(
            sig, nProcesses ) )
    stdout.flush()
    print( "Results: {}\n".format( naiveBadEdge(
        sig, nProcesses, interval, isCore, "Core edges" ) ) )

