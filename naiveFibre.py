"""
Naively search for a 1-vertex triangulation of a small Seifert fibre space
that has no "bad" edges. Here, we consider an edge to be bad if we cannot
rule out the possibility that it is isotopic to a Seifert fibre.

To do this, we start with a triangulation that has at least one bad edge. We
then search the Pachner graph for a triangulation with one less bad edge.

The search naively prioritises triangulations with a small number of
tetrahedra. The main purpose is to demonstrate the effectiveness of the more
sophisticated heuristic used in removeFibre.py.
"""
import sys
from triangCounterexHelpers import possibleFibre
from naiveBadEdge import naiveBadEdge


if __name__ == "__main__":
    sig = sys.argv[1]
    nProcesses = int( sys.argv[2] )
    interval = int( sys.argv[3] )
    try:
        badEdgeString = sys.argv[4]
    except IndexError:
        badEdges = None
    else:
        # Parse badEdgeString as a list of edge indices.
        badEdges = [ int(e) for e in badEdgeString.split(",") ]
    print( "Naively removing {}. Sig: {}, #Processes: {}.".format(
        "Seifert fibres", sig, nProcesses ) )
    sys.stdout.flush()
    print( "Results: {}\n".format( naiveBadEdge(
        sig, nProcesses, interval, possibleFibre,
        "Seifert fibres", badEdges ) ) )

