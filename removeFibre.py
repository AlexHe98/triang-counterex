"""
Search for a 1-vertex triangulation of a small Seifert fibre space that has
no "bad" edges, where we consider an edge to be bad if we cannot rule out the
possibility that it is isotopic to a Seifert fibre.

To do this, we start with a triangulation that has at least one bad edge. We
then search the Pachner graph for a triangulation with one less bad edge. The
search prioritises triangulations in which:
(1) the maximum "multiplicity defect" among bad edges is small (see the
    mulDefect() routine);
(2) the maximum "degree defect" among bad edges is small (see the degDefect()
    routine); and
(3) the number of tetrahedra is small.
"""
from sys import argv, stdout
from triangCounterexHelpers import possibleFibre
from removeBadEdge import removeBadEdge


if __name__ == "__main__":
    sig = argv[1]
    nProcesses = int( argv[2] )
    interval = int( argv[3] )
    print( "Removing Seifert fibres. Sig: {}, #Processes: {}.".format(
        sig, nProcesses ) )
    stdout.flush()
    print( "Results: {}\n".format( removeBadEdge(
        sig, nProcesses, interval, possibleFibre, "Seifert fibres", 1 ) ) )

