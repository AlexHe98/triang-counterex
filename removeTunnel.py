"""
Search for a 1-vertex ideal triangulation of a knot complement in which no
edge forms a "tunnel"; in other words, for each edge e of the triangulation,
pinching e yields an ideal triangulation that does *not* represent an
orientable genus-2 handlebody.

To do this, we start with a triangulation that has at least one tunnel edge;
for example, the triangulation given by the isomorphism signature "cPcbbbadu"
is an ideal triangulation of the trefoil knot with two edges, exactly one of
which is a tunnel. We then search the Pachner graph for a 3-2 move that
removes a tunnel edge, prioritising triangulations in which:
(1) the maximum "multiplicity defect" among tunnel edges is small (see the
    mulDefect() routine);
(2) the maximum "degree defect" among tunnel edges is small (see the
    degDefect() routine); and
(3) the number of tetrahedra is small.

This is written to run with either Python 2 or Python 3, without needing to
rewrite any code.
"""
from __future__ import print_function
from sys import argv, stdout
from triangCounterexHelpers import isTunnel
from removeBadEdge import removeBadEdge


if __name__ == "__main__":
    sig = argv[1]
    nProcesses = int( argv[2] )
    interval = int( argv[3] )
    print( "Removing tunnels. Sig: {}. #Processes: {}.".format(
            sig, nProcesses ) )
    stdout.flush()
    print( "Results: {}\n".format( removeBadEdge(
        sig, nProcesses, interval, isTunnel, "Tunnels" ) ) )

