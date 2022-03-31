"""
Search for a 1-vertex triangulation of a lens space in which no edge forms a
"core loop"; in other words, for each edge e of the triangulation, pinching
e yields an ideal triangulation that does *not* represent the solid torus.

To do this, we start with a triangulation that has at least one core edge. We
then search the Pachner graph for a 3-2 move that removes a core edge.
"""
from sys import argv, stdout
from triangCounterexHelpers import isCore
from removeBadEdge import removeBadEdge


if __name__ == "__main__":
    sig = argv[1]
    nProcesses = int( argv[2] )
    interval = int( argv[3] )
    print( "Removing core edges. Sig: {}. #Processes: {}.".format(
            sig, nProcesses ) )
    stdout.flush()
    print( "Results: {}\n".format( removeBadEdge(
        sig, nProcesses, interval, isCore, "Core edges" ) ) )

