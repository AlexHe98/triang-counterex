"""
Simplify a 1-vertex triangulation of a small Seifert fibre space without
ever performing a 2-3 move that introduces a "bad" edge; we consider an edge
to be bad if we cannot rule out the possibility that it is isotopic to a
Seifert fibre.
"""
from sys import argv, stdout
from triangCounterexHelpers import possibleFibre
from simplify import simplifyNoBadEdge


if __name__ == "__main__":
    sig = argv[1]
    nProcesses = int( argv[2] )
    interval = int( argv[3] )
    height = int( argv[4] )
    print( "Simplifying without creating edges isotopic to Seifert fibres." +
            " Sig: {}. #Processes: {}. Height: {}.".format(
                sig, nProcesses, height ) )
    stdout.flush()
    print( "Results: {}\n".format( simplifyNoBadEdge(
        sig, nProcesses, interval, possibleFibre, height ) ) )
