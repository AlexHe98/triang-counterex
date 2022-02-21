"""
Search for a 1-vertex ideal triangulation of a knot complement in which no
edge forms a "tunnel"; in other words, for each edge e of the triangulation,
pinching e yields an ideal triangulation that does *not* represent an
orientable genus-2 handlebody.

To do this, we start with a triangulation that has at least one tunnel edge;
for example, the triangulation given by the isomorphism signature "cPcbbbadu"
is an ideal triangulation of the trefoil knot with two edges, exactly one of
which is a tunnel. We then search the Pachner graph for a 3-2 move that
removes a tunnel edge, prioritising triangulations in which:
(1) the maximum "multiplicity defect" among tunnel edges is small (see the
    mulDefect() routine);
(2) the maximum "degree defect" among tunnel edges is small (see the
    degDefect() routine); and
(3) the number of tetrahedra is small.

This is written to run with either Python 2 or Python 3, without needing to
rewrite any code.
"""
from sys import argv
from multiprocessing import Process
from removeTunnel import BADNAME
from removeTunnel import RemoveBadEdgeManager
from removeTunnel import removeBadEdgeLoop


def removeBadEdge( s, nProcesses, interval ):
    """
    Starting at the given isomorphism signature s, and using the given number
    of parallel processes, searches for a 3-2 move that removes a bad edge
    from s, where we consider an edge e to be bad if and only if isBad(e)
    returns True.

    The search abides by the following constraints:
    (1) Any 2-3 moves that introduce bad edges are ignored.
    (2) Priority is given to triangulations in which the maximum
        "multiplicity defect" (see the mulDefect() routine) among bad edges
        is small, then to triangulations in which the maximum "degree defect"
        (see the degDefect() routine) among bad edges is small, and then to
        triangulations with fewer tetrahedra.
    (3) Triangulations with equal priority are processed in order of
        insertion.

    If the search successfully finds one or more 3-2 moves that remove a bad
    edge, then returns a list containing the isomorphism signatures of the
    triangulations that result from performing these moves (it is possible
    for multiple parallel processes to each find such a move at roughly the
    same time, in which case the returned list will contain more than one
    isomorphism signature).

    It is also possible for the search to "get stuck" in a situation where
    the only moves available are 2-3 moves that introduced bad edges (which
    we ignore). In this case, returns an empty list.

    Note: When run with multiple processes (i.e., nProcesses > 1), this
    routine is not completely deterministic because insertion order may vary
    across different run-throughs.

    The given interval determines the number of seconds that we should wait
    before printing some info to standard output.

    Pre-condition:
    --- The triangulation represented by s has at least one bad edge.
    """
    manager = RemoveBadEdgeManager()
    manager.start()
    lock = manager.Lock()
    cond = manager.Condition(lock)
    shared = manager.RemoveBadEdgeData( s, nProcesses, interval )
    processes = [
            Process( target=removeBadEdgeLoop, args=( shared, lock, cond ) )
            for _ in range(nProcesses) ]
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    # Print some info, and then return results.
    print( shared.info() )
    results = shared.getResults()
    for r, comp in results:
        rem = shared.badEdges(r)
        print( "Result: {}. Complexity: {}. {}: {}".format(
            r, comp, BADNAME, rem ) )
        path = shared.backtrack(r)
        for s in path:
            print( "    {}".format(s) )
    return results


if __name__ == "__main__":
    sig = argv[1]
    nProcesses = int( argv[2] )
    interval = int( argv[3] )
    print( "Results: {}\n".format( removeBadEdge(
        sig, nProcesses, interval ) ) )
