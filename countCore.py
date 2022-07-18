"""
Given a collection of isomorphism signatures, finds the minimum number of
core edges that occurs in a triangulation that arises from an isomorphism
signature in this list; also counts how many triangulations achieve this
minimum number.

By "core edge", we mean an edge e of the triangulation such that pinching e
yields an ideal triangulation of the solid torus.
"""
import fileinput
from sys import stdout
from timeit import default_timer
from regina import *
from triangCounterexHelpers import isCore

if __name__ == "__main__":
    start = default_timer()
    prev = start
    interval = 21600 # Print update every 6 hours.
    tested = 0
    minCore = None
    found = 0
    msg = "Time: {:.6f}. Tested: {}. Min: {}. Found: {}.{}"
    for line in fileinput.input():
        time = default_timer()
        if time - prev > interval:
            prev = time
            print( msg.format( time - start, tested, minCore, found, "" ) )
            stdout.flush()
        tested += 1
        sig = line.rstrip()
        tri = Triangulation3.fromIsoSig(sig)
        coreCount = 0
        for edge in tri.edges():
            if isCore(edge):
                coreCount += 1
                # Break out early wherever possible.
                if ( minCore is not None ) and ( coreCount > minCore ):
                    break
        if ( minCore is None ) or ( coreCount < minCore ):
            found = 1
            minCore = coreCount
            prev = default_timer()
            print( msg.format(
                prev - start, tested, minCore, found, " " + sig ) )
            stdout.flush()
        elif coreCount == minCore:
            found += 1
    print( msg.format(
        default_timer() - start, tested, minCore, found, "" ) )
    print()
