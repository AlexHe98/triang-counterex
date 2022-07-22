"""
Given a series of isomorphism via standard input, searches for isomorphism
signatures whose triangulations have no core edges.

By "core edge", we mean an edge e of the triangulation such that pinching e
yields an ideal triangulation of the solid torus.
"""
from sys import argv, stdin
from timeit import default_timer
from multiprocessing import Pool, Array
from bruteForceHelpers import setup, bruteCoreTest


if __name__ == "__main__":
    start = default_timer()
    doubles = Array( "d", [
        start, # Start time.
        start ] ) # Previous update time (initially equal to start time).
    ints = Array( "i", [
        21600, # Print update every 6 hours.
        0, # Number of iso sigs tested.
        0 ] ) # Number of iso sigs that have no core edges.

    # For each iso sig that we read from standard input, check whether the
    # corresponding triangulation has any core edges.
    with Pool(
            processes=int( argv[1] ),
            initializer=setup,
            initargs=[doubles, ints] ) as pool:
        for _ in pool.imap( bruteCoreTest, stdin ):
            pass
    print( "Time: {:.6f}. Tested: {}. Found: {}.".format(
        default_timer() - start, ints[1], ints[2] ) )
    print()
