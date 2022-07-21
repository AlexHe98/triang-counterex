from sys import stdout
from timeit import default_timer
from regina import *
from triangCounterexHelpers import isCore


def setup(d, i):
    global doubles, ints
    doubles = d
    ints = i


def testSig(line):
    prev = 1
    interval = 0
    tested = 1
    minCore = 2
    found = 3
    msg = "Time: {:.6f}. Tested: {}. Min: {}. Found: {}.{}"
    with ints.get_lock():
        # When ints and doubles are both locked, always lock ints first.
        with doubles.get_lock():
            time = default_timer()
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
