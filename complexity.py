"""
Convenience script for computing complexities of triangulations.
"""
import sys
from triangCounterexHelpers import isCore, isTunnel, possibleFibre
from removeBadEdge import complexity


if __name__ == "__main__":
    badEdge = sys.argv[1]
    if badEdge == "core":
        isBad = isCore
    elif badEdge == "tunnel":
        isBad = isTunnel
    elif badEdge == "fibre":
        isBad = possibleFibre
    else:
        sys.exit( "Unknown bad edge type: {}".format(badEdge) )

    tri = Triangulation3.fromIsoSig( sys.argv[2] )
    badEdges = [ e.index() for e in tri.edges() if isBad(e) ]
    print( complexity( tri, badEdges ) )
