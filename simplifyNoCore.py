"""
Simplify a 1-vertex triangulation of a lens space without ever performing a
2-3 move that introduces an edge that forms a "core loop".
"""
from sys import argv, stdout
from triangCounterexHelpers import isCore
from simplify import simplifyNoBadEdge


if __name__ == "__main__":
    sig = argv[1]
    nProcesses = int( argv[2] )
    interval = int( argv[3] )
    height = int( argv[4] )
    print( "Simplifying without creating core edges." +
            " Sig: {}. #Processes: {}. Height: {}.".format(
                sig, nProcesses, height ) )
    stdout.flush()
    print( "Results: {}\n".format( simplifyNoBadEdge(
        sig, nProcesses, interval, isCore, height ) ) )
