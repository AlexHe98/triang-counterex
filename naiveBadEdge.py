"""
Given a triangulation with at least one "bad" edge (where "bad" stands for
any computable topological property of an edge), search for a 3-2 move that
removes a bad edge.

This script naively prioritises triangulations with fewer tetrahedra. The
main purpose is to demonstrate the effectiveness of the more sophisticated
heuristic used in removeBadEdge.py.
"""
from sys import stdout
from timeit import default_timer
from multiprocessing import Condition, Process
from multiprocessing.managers import SyncManager
from queue import Empty
from regina import *
from search import SearchData
from pachner import twoThree, threeTwo
from triangCounterexHelpers import relabelEdges


class NaiveBadEdgeData(SearchData):
    """
    Stores all the data that needs to be shared between processes when
    running the naiveBadEdge() routine.
    """
    def __init__( self, s, nProcesses, interval, isBad, badEdges ):
        """
        Initialises shared search data.

        Among other things, enqueues the initial isomorphism signature s, and
        provides the isBad() routine for processes to use to test whether an
        edge is bad.
        """
        super().__init__( s, nProcesses, interval, isBad )

        # Initial triangulation.
        t = Triangulation3.fromIsoSig(s)
        if badEdges is None:
            badEdges = [ e.index() for e in t.edges() if isBad(e) ]
        self._height = t.size()
        badEdgeCount = len(badEdges)

        # Our target is to get below the current number of bad edges.
        self._target = badEdgeCount

        # Enqueue the initial isomorphism signature.
        self.enqueue( s, self._height, badEdges,
                None ) # Initial iso sig has no source.

    def info(self):
        """
        Returns a string describing the time elapsed, the number of
        triangulations discovered, the number of triangulations explored, and
        the maximum height reached in the Pachner graph.
        """
        return super().info() + " Height: {}.".format( self._height )

    def _printInfo( self, check ):
        time = default_timer()
        if check and time - self._previousInfo < self._interval:
            return
        self._previousInfo = time
        print( self.info() )
        stdout.flush()

    def _updateHeight( self, height ):
        if height > self._height:
            self._height = height
            self._printInfo(False)

    def target(self):
        """
        Returns the target number of bad edges that we need to get below.
        """
        return self._target

    def newResult( self, sig, height, badEdges, source ):
        """
        Records the given sig as a new result, and instructs all processes to
        stop as soon as possible.
        """
        self._results.append( ( sig, height ) )
        self._details[sig] = ( badEdges, source, height )
        self._updateHeight(height)
        self._stop = True

    def badEdges( self, sig ):
        """
        Returns the tuple of bad edge indices in sig.

        Pre-condition:
        --> The given sig was enqueued at some point in time.
        """
        return self._details[sig][0]

    def enqueue( self, sig, height, badEdges, source ):
        """
        Pushes the given sig onto the queue, assigns badEdges as the tuple
        returned by self.badEdges(sig), and records the given source
        information.
        """
        self._queue.put(
                ( height, self._total, sig ),
                False )
        self._total += 1
        self._details[sig] = ( badEdges, source, height )
        self._updateHeight(height)

    def dequeue(self):
        """
        Prints some info if a lot of time has passed, and then pops the
        isomorphism signature with the highest priority from the queue.

        Exceptions:
        --> Raises Empty if the queue is empty.
        """
        self._printInfo(True)

        # Try to pop from queue immediately, since we don't want to change
        # anything if it turns out that the queue is currently empty.
        height, _, sig = self._queue.get(False)
        self._explored += 1
        source = self._details[sig][1]
        return ( sig, height, source )

    def backtrack( self, sig ):
        """
        Returns a list that describes the sequence of 2-3 and 3-2 moves that
        produced the given sig.
        """
        back = []
        source = self._details[sig][1]
        while source is not None:
            back.append(source)
            source = self._details[ source[3] ][1]
        return back[::-1]


class NaiveBadEdgeManager(SyncManager):
    pass
NaiveBadEdgeManager.register( "NaiveBadEdgeData", NaiveBadEdgeData )


def naiveBadEdgeGreedy_(
        tri, sig, height, source, badEdges, target, shared, lock, cond ):
    if source is None:
        steps = 1
    else:
        steps = source[0] + 1

    # Try to remove bad edges using 3-2 moves.
    newHeight = height - 1
    atTarget = ( len(badEdges) == target )
    for i in badEdges:
        result = threeTwo( tri.edge(i) )
        if result is None:
            continue

        # We successfully removed a bad edge!
        newTri, newTrans = result
        newSig, newIsom = newTri.isoSigDetail()
        with lock:
            if shared.alreadySeen(newSig):
                continue

        # We haven't seen newSig before. Did we get below the target?
        newBadEdges = [ newTrans[j] for j in badEdges if j != i ]
        newSource = ( steps,
                ( "3-2", i ),
                height,
                sig )
        relBadEdges = tuple( relabelEdges( newTri, newBadEdges, newIsom ) )
        if atTarget:
            # We started at the target, so removing a bad edge gets us below!
            with lock:
                shared.newResult(
                        newSig, newHeight, relBadEdges, newSource )
            return True

        # We didn't get below the target, so we should continue trying to
        # remove bad edges.
        with lock:
            if shared.queueEmpty():
                shared.enqueue(
                        newSig, newHeight, relBadEdges, newSource )
                cond.notify_all()
            else:
                shared.enqueue(
                        newSig, newHeight, relBadEdges, newSource )
        if naiveBadEdgeGreedy_( newIsom(newTri), newSig, newHeight,
                newSource, relBadEdges, target, shared, lock, cond ):
            return True
    return False


def naiveBadEdgeExplore_( sig, height, source, shared, lock, cond ):
    # NOTES
    """
    --> Many of the computations can be done independently of other
        processes, but be careful to acquire the given lock every time we
        need to access shared.
    """
    tri = Triangulation3.fromIsoSig(sig)
    with lock:
        badEdges = shared.badEdges(sig)
        target = shared.target()

    # If sig was generated by the search (i.e., source is not None), then we
    # would already have greedily attempted to remove bad edges.
    if source is None:
        if naiveBadEdgeGreedy_( tri, sig, height, source, badEdges, target,
                shared, lock, cond ):
            # Terminate immediately if we found a result.
            return
        steps = 1
    else:
        steps = source[0] + 1

    # Try all possible 3-2 moves on *good* edges.
    newHeight = height - 1
    for edge in tri.edges():
        if edge.index() in badEdges:
            continue
        result = threeTwo(edge)
        if result is None:
            continue
        newTri, newTrans = result
        newSig, newIsom = newTri.isoSigDetail()

        # If newSig is genuinely new, then we should enqueue it.
        with lock:
            if shared.alreadySeen(newSig):
                continue
        newBadEdges = [ newTrans[i] for i in badEdges ]
        newSource = ( steps,
                ( "3-2", edge.index() ),
                height,
                sig )
        relBadEdges = tuple( relabelEdges(
            newTri, newBadEdges, newIsom ) )
        with lock:
            if shared.queueEmpty():
                shared.enqueue(
                        newSig, newHeight, relBadEdges, newSource )
                cond.notify_all()
            else:
                shared.enqueue(
                        newSig, newHeight, relBadEdges, newSource )

        # Can we remove bad edges?
        if naiveBadEdgeGreedy_( newIsom(newTri), newSig, newHeight,
                newSource, relBadEdges, target, shared, lock, cond ):
            # Terminate immediately if we found a result.
            return

    # Now try all 2-3 moves.
    with lock:
        isBad = shared.getIsBadRoutine()
    newHeight = height + 1
    for triangle in tri.triangles():
        result = twoThree(triangle)
        if result is None:
            continue
        newTri, newTrans = result
        newSig, newIsom = newTri.isoSigDetail()

        # Have we already seen newSig?
        with lock:
            if shared.alreadySeen(newSig):
                continue

        # Did we introduce a bad edge?
        newBadEdges = [ newTrans[i] for i in badEdges ]
        if isBad( newTri.edge( newTrans[-1] ) ):
            newBadEdges.append( newTrans[-1] )
        
        # We should enqueue the newSig.
        newSource = ( steps,
                ( "2-3", triangle.index() ),
                height,
                sig )
        relBadEdges = tuple( relabelEdges(
            newTri, newBadEdges, newIsom ) )
        with lock:
            if shared.queueEmpty():
                shared.enqueue(
                        newSig, newHeight, relBadEdges, newSource )
                cond.notify_all()
            else:
                shared.enqueue(
                        newSig, newHeight, relBadEdges, newSource )

        # Can we remove bad edges?
        if naiveBadEdgeGreedy_( newIsom(newTri), newSig, newHeight,
                newSource, relBadEdges, target, shared, lock, cond ):
            # Terminate immediately if we found a result.
            return


def naiveBadEdgeLoop_( shared, lock, cond ):
    # NOTES
    """
    --> To minimise testing whether edges are bad (which could potentially be
        very expensive), take advantage of transitions to keep track of which
        edges are bad.
    --> *Wait* (but don't stop) if the shared queue is empty but there is
        still at least one other process running (since new sigs could still
        be added to the queue).
    --> Two reasons to stop:
        (1) We have been instructed to stop (because we have successfully
            reduced the number of bad edges).
        (2) The shared queue is empty and no other processes are running.
    """
    while True:
        # Keep exploring sigs until we are instructed to stop.
        with lock:
            if shared.stop():
                # We have been instructed to stop. Wake up any waiting
                # processes so that everyone else can stop too.
                cond.notify_all()
                return

            # Try popping the next sig.
            try:
                sig, height, source = shared.dequeue()
            except Empty:
                # There is currently nothing to explore in the shared queue,
                # so tell shared that we are about to pause this process.
                shared.tellPaused()
                if shared.nRunning() == 0:
                    # Every other process is already paused. Wake them up so
                    # that everyone can stop.
                    cond.notify_all()
                    return

                # There is still at least one other process running, so wait
                # for other processes to either:
                # - stop; or
                # - push more sigs onto the shared queue.
                # Note: cond.wait() should automatically release the given
                # lock, and should automatically re-acquire this lock when
                # awakened.
                cond.wait()

                # We woke up! There are two possibilities from here:
                # - shared.stop() is True or shared.nRunning() == 0, in which
                #   case we stop immediately.
                # - shared.nRunning() > 0, in which case there should now be
                #   new sigs in the shared queue to explore, so we should
                #   resume.
                if shared.stop() or shared.nRunning() == 0:
                    return

                # Tell shared that we are resuming this process, and then
                # continue exploring.
                shared.tellResumed()
                continue

        # Explore sig. We do most of this independently of other processes.
        naiveBadEdgeExplore_( sig, height, source, shared, lock, cond )


def naiveBadEdge( s, nProcesses, interval, isBad, badName, badEdges=None ):
    """
    Starting at the given isomorphism signature s, and using the given number
    of parallel processes, naively searches the Pachner graph for a
    triangulation with fewer bad edges than s.

    We consider an edge e to be bad if and only if isBad(e) returns True.

    The search prioritises triangulations with fewer tetrahedra. When the
    number of tetrahedra is equal, the triangulations are processed in order
    of insertion into a queue.

    If the search successfully finds one or more triangulations with fewer
    bad edges than s, then this routine returns a list containing all of the
    corresponding isomorphism signatures (it is possible for multiple
    parallel processes to each find such a triangulation at roughly the same
    time, which is why the returned list may contain more than one
    isomorphism signature).

    Note: When run with multiple processes (i.e., nProcesses > 1), this
    routine is not completely deterministic because insertion order may vary
    across different run-throughs.

    The given interval determines the number of seconds that we should wait
    before printing some info to standard output.

    There is an option to provide a pre-computed list of bad edges for s, by
    specifying badEdges to be a list of edge indices; in this case, this
    routine assumes that this list correctly identifies the bad edges in s.
    Otherwise, if badEdges is None (the default), then this routine will do
    the work of computing the bad edges.

    Pre-condition:
    --> The triangulation represented by s has at least one bad edge.
    --> The function isBad takes a single edge as input, and
        deterministically outputs either True or False (corresponding to
        whether the edge should be considered bad).
    """
    manager = NaiveBadEdgeManager()
    manager.start()
    lock = manager.Lock()
    cond = manager.Condition(lock)
    shared = manager.NaiveBadEdgeData(
            s, nProcesses, interval, isBad, badEdges )
    processes = [
            Process( target=naiveBadEdgeLoop_,
                args=( shared, lock, cond ) )
            for _ in range(nProcesses) ]
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    # Print some info, and then return results.
    print( shared.info() )
    results = shared.getResults()
    for r, nTet in results:
        rem = shared.badEdges(r)
        print( "Result: {}. Size: {}. {}: {}".format(
            r, nTet, badName, rem ) )
        path = shared.backtrack(r)
        for s in path:
            print( "    {}".format(s) )
    return results
