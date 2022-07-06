"""
Given a collection of isomorphism signatures, searches for a signature whose
corresponding triangulation is 1-vertex and has no edge that forms a "core
loop"; in other words, for each edge e of the triangulation, pinching e
yields an ideal triangulation that does *not* represent the solid torus.
"""
import fileinput
from sys import stdout
from timeit import default_timer
from regina import *
from triangCounterexHelpers import isCore

if __name__ == "__main__":
    start = default_timer()
    prev = start
    interval = 10800 # Print update every 3 hours.
    tested = 0
    found = 0
    msg = "Time: {:.6f}. Tested: {}. Found: {}.{}"
    for line in fileinput.input():
        time = default_timer()
        if time - prev > interval:
            prev = time
            print( msg.format( time - start, tested, found, "" ) )
        tested += 1
        sig = line.rstrip()
        tri = Triangulation3.fromIsoSig(sig)
        noCore = True
        for edge in tri.edges():
            if isCore(edge):
                noCore = False
                break
        if noCore:
            found += 1
            prev = default_timer()
            print( msg.format( prev - start, tested, found, " " + sig ) )
    print( msg.format( default_timer() - start, tested, found, "" ) )
    print()
