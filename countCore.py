"""
Given a series of isomorphism signatures via standard input, finds the
minimum number of core edges that occurs in a triangulation that arises from
one of these isomorphism signatures. Also counts how many triangulations
achieve this minimum number.

By "core edge", we mean an edge e of the triangulation such that pinching e
yields an ideal triangulation of the solid torus.
"""
from sys import argv, stdin, stdout
from timeit import default_timer
from multiprocessing import Pool, Value, Array
from regina import *
from triangCounterexHelpers import isCore

if __name__ == "__main__":
    start = default_timer()
    doubles = Array( "d", [
        start, # Start time.
        start ] ) # Previous update time (initially equal to start time).
    ints = Array( "i", [
        21600, # Print update every 6 hours.
        0, # Number of iso sigs tested.
        -1, # Minimum number of core edges (negative indicates no bound).
        0 ] ) # Number of iso sigs that achieve the current minimum. 

    def testSig(line):
        time = default_timer()
        prev = 1
        interval = 0
        tested = 1
        minCore = 2
        found = 3
        msg = "Time: {:.6f}. Tested: {}. Min: {}. Found: {}.{}"
        with ints.get_lock():
            # When ints and doubles are both locked, always lock ints first.
            with doubles.get_lock():
                if time - doubles[prev] > ints[interval]:
                    doubles[prev] = time
                    print( msg.format( time - doubles[0],
                        ints[tested], ints[minCore], ints[found], "" ) )
                    stdout.flush()
            ints[tested] += 1
            currentMin = ints[minCore]
        sig = line.rstrip()
        tri = Triangulation3.fromIsoSig(sig)
        coreCount = 0
        for edge in tri.edges():
            if isCore(edge):
                coreCount += 1
                # Break out early wherever possible.
                if ( currentMin >= 0 ) and ( coreCount > currentMin ):
                    break
        with ints.get_lock():
            if ( ints[minCore] < 0 ) or ( coreCount < ints[minCore] ):
                ints[found] = 1
                ints[minCore] = coreCount
                # When ints and doubles are both locked, always lock ints
                # first.
                with doubles.get_lock():
                    doubles[prev] = default_timer()
                    print( msg.format( doubles[prev] - doubles[0],
                        ints[tested], ints[minCore], ints[found],
                        " " + sig ) )
                    stdout.flush()
            elif coreCount == ints[minCore]:
                ints[found] += 1

    # Read iso sigs from standard input, and distribute the testing workload
    # across the given number of processes.
    with Pool( processes=int( argv[1] ) ) as pool:
        pool.imap( testSig, stdin )
    print( "Time: {:.6f}. Tested: {}. Min: {}. Found: {}.{}".format(
        default_timer() - start, ints[1], ints[2], ints[3], "" ) )
    print()
