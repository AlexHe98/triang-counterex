"""
Given a series of isomorphism signatures via standard input, finds the
minimum number of core edges that occurs in a triangulation that arises from
one of these isomorphism signatures. Also counts how many triangulations
achieve this minimum number.

By "core edge", we mean an edge e of the triangulation such that pinching e
yields an ideal triangulation of the solid torus.
"""
from sys import argv, stdin
from timeit import default_timer
from multiprocessing import Pool, Array
from bruteForceHelpers import setup, countCoreTest


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

    # For each iso sig that we read from standard input, count the number of
    # core edges in the corresponding triangulation (but terminate early if
    # we find too many core edges).
    with Pool(
            processes=int( argv[1] ),
            initializer=setup,
            initargs=[doubles, ints] ) as pool:
        for _ in pool.imap( countCoreTest, stdin ):
            pass
    print( "Time: {:.6f}. Tested: {}. Min: {}. Found: {}.".format(
        default_timer() - start, ints[1], ints[2], ints[3] ) )
    print()
