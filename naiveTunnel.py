"""
Naively search for a 1-vertex ideal triangulation of a knot complement that
has no "bad" edges. Here, we consider an edge e to be bad if it forms a
"tunnel", meaning that pinching e yields an ideal triangulation of the
orientable genus-2 handlebody.

To do this, we start with a triangulation that has at least one bad edge;
for example, the triangulation given by the isomorphism signature "cPcbbbadu"
is an ideal triangulation of the trefoil knot with two edges, exactly one of
which is a tunnel. We then search the Pachner graph for a triangulation with
one less bad edge.

The search naively prioritises triangulations with a small number of
tetrahedra. The main purpose is to demonstrate the effectiveness of the more
sophisticated heuristic used in removeTunnel.py.
"""
from sys import argv, stdout
from triangCounterexHelpers import isTunnel
from naiveBadEdge import naiveBadEdge


if __name__ == "__main__":
    sig = argv[1]
    nProcesses = int( argv[2] )
    interval = int( argv[3] )
    print( "Naively removing tunnels. Sig: {}. #Processes: {}.".format(
            sig, nProcesses ) )
    stdout.flush()
    print( "Results: {}\n".format( naiveBadEdge(
        sig, nProcesses, interval, isTunnel, "Tunnels" ) ) )

