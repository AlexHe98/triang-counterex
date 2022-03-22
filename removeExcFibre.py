"""
Search for a 1-vertex triangulation of a small Seifert fibre space that has
no "bad" edges, where we consider an edge to be bad if it is isotopic to an
exceptional fibre.

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
from seifert import isExceptionalFibre
from removeBadEdge import removeBadEdge


if __name__ == "__main__":
    sig = argv[1]
    nProcesses = int( argv[2] )
    interval = int( argv[3] )
    try:
        minArg = argv[4]
    except IndexError:
        useMin = False
    else:
        if minArg == "min":
            useMin = True
        else:
            print( "    Unknown argument \"{}\".".format(minArg) )
            useMin = False
    print( "Removing {}. Sig: {}, #Processes: {}. UseMin?: {}".format(
        "exceptional fibres", sig, nProcesses, useMin ) )
    stdout.flush()
    print( "Results: {}\n".format( removeBadEdge(
        sig, nProcesses, interval, isExceptionalFibre,
        "Exceptional fibres", 1, useMin ) ) )

