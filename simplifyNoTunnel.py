"""
Simplify a 1-vertex ideal triangulation of a knot complement without ever
performing a 2-3 move that introduces an edge that forms a "tunnel".
"""
from sys import argv, stdout
from triangCounterexHelpers import isTunnel
from simplify import simplifyNoBadEdge


if __name__ == "__main__":
    sig = argv[1]
    nProcesses = int( argv[2] )
    interval = int( argv[3] )
    height = int( argv[4] )
    print( "Simplifying without creating tunnel edges." +
            " Sig: {}. #Processes: {}. Height: {}.".format(
                sig, nProcesses, height ) )
    stdout.flush()
    print( "Results: {}\n".format( simplifyNoBadEdge(
        sig, nProcesses, interval, isTunnel, height ) ) )
