"""
Helper routines for running multiprocessed brute-force searches.
"""
from sys import stdout
from timeit import default_timer
from regina import *
from triangCounterexHelpers import isCore


def setup(d, i):
    global doubles, ints
    doubles = d
    ints = i


def countCoreTest(line):
    start = 0 # doubles[start] == start time.
    prev = 1 # doubles[prev] == previous update time.
    interval = 0 # ints[interval] == update interval.
    tested = 1 # ints[tested] == total number of iso sigs processed.
    minCore = 2 # ints[minCore] == current minimum number of core edges.
    found = 3 # ints[found] == number of sigs that achieve current minimum.
    msg = "Time: {:.6f}. Tested: {}. Min: {}. Found: {}.{}"
    with ints.get_lock():
        # When ints and doubles are both locked, always lock ints first.
        with doubles.get_lock():
            time = default_timer()
            if time - doubles[prev] > ints[interval]:
                doubles[prev] = time
                print( msg.format( time - doubles[start],
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
                print( msg.format( doubles[prev] - doubles[start],
                    ints[tested], ints[minCore], ints[found],
                    " " + sig ) )
                stdout.flush()
        elif coreCount == ints[minCore]:
            ints[found] += 1
