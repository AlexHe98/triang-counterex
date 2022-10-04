"""
Given a triangulation with at least one "bad" edge (where "bad" stands for
any computable topological property of an edge), search for a 3-2 move that
removes a bad edge, prioritising triangulations in which:
(1) the maximum "multiplicity defect" among bad edges is small (see the
    mulDefect() routine);
(2) the maximum "degree defect" among bad edges is small (see the degDefect()
    routine); and
(3) the number of tetrahedra is small.
"""
from sys import stdout
from timeit import default_timer
from multiprocessing import Condition, Process
from multiprocessing.managers import SyncManager
from queue import Empty
from regina import *
from search import SearchData
from pachner import twoThree, threeTwo
from triangCounterexHelpers import mulDefect
from triangCounterexHelpers import degDefect
from triangCounterexHelpers import relabelEdges


def complexity( tri, badEdges, excess=None ):
    """
    Returns a measure for how far away we are from removing the given
    collection of bad edges from the given triangulation.

    In general, our ultimate goal is to remove all of the bad edges in a
    triangulation, which means that we should try to simultaneously minimise
    the complexity of all the bad edges. To take all of the bad edges into
    account, the value of excess should be set to None (this is the default
    value). In this case, this routine returns a tuple consisting of:
    (0) The number of bad edges (i.e., the length of the given collection of
        bad edge indices).
    (1) The maximum "multiplicity defect" among the given bad edges (see the
        mulDefect() routine).
    (2) The maximum "degree defect" among the given bad edges (see the
        degDefect() routine).
    (3) The number of tetrahedra in the given triangulation.

    Sometimes, it may be helpful to focus on an intermediate goal: reducing
    the number of bad edges by one. So, if we started with a triangulation
    that had n bad edges, and if the given triangulation tri now has n+x bad
    edges, for some non-negative integer x, we might like to focus only on
    minimising the complexity of x+1 bad edges. To do this, the value of
    excess should be set to the integer x. In this case, this routine returns
    a tuple consisting of:
    (0) The number of bad edges.
    (1) The maximum "multiplicity defect" among the n+1 bad edges with the
        smallest such defect.
    (2) The maximum "degree defect" among the n+1 bad edges with the smallest
        such defect.
    (3) The number of tetrahedra in the given triangulation.
    """
    if badEdges:
        mulDefs = [ mulDefect( tri.edge(i) ) for i in badEdges ]
        degDefs = [ degDefect( tri.edge(i) ) for i in badEdges ]
    else:
        mulDefs = [0]
        degDefs = [0]
    if excess is None:
        return ( len(badEdges),
                max(mulDefs),
                max(degDefs),
                tri.size() )
    else:
        # Use minimums instead of maximums.
        trunc = max( 0, excess )
        return ( len(badEdges),
                sorted(mulDefs)[trunc],
                sorted(degDefs)[trunc],
                tri.size() )


class RemoveBadEdgeData(SearchData):
    """
    Stores all the data that needs to be shared between processes when
    running the removeBadEdge() routine.
    """
    def __init__( self, s, nProcesses, interval, isBad, excess, inter,
            badEdges ):
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
        if inter:
            c = complexity( t, badEdges, 0 )
        else:
            c = complexity( t, badEdges )
        self._minSig = { c[0]: s }
        self._current = s

        # The current number of bad edges is given by c[0]. Our target is to
        # get below this number. Moreover, during our search, we never allow
        # the number of bad edges to exceed c[0]+excess.
        self._target = c[0]
        self._maxBadEdges = self._target + excess

        # Enqueue the initial isomorphism signature.
        self.enqueue( s, c, tuple(badEdges),
                None ) # Initial iso sig has no source.

    def _minComplexity( self, badEdgeCount ):
        return self._details[ self._minSig[badEdgeCount] ][2]

    def info(self):
        """
        Returns a string describing the time elapsed, the number of
        triangulations discovered, the number of triangulations explored, and
        details of the triangulation with the smallest complexity that has
        been discovered so far.
        """
        minFormat = " Min{}: {}, {}."
        msg = super().info()
        for badEdgeCount in self._minSig:
            msg += minFormat.format(
                    badEdgeCount,
                    self._minComplexity(badEdgeCount),
                    self._minSig[badEdgeCount] )
        msg += " Current: {}, {}.".format(
                self._details[ self._current ][2],
                self._current )
        return msg

    def _updateMinSig( self, sig, priority ):
        badEdgeCount = priority[0]
        if ( ( badEdgeCount not in self._minSig ) or
                ( priority < self._minComplexity(badEdgeCount) ) ):
            self._minSig[badEdgeCount] = sig

    def target(self):
        """
        Returns the target number of bad edges that we need to get below.
        """
        return self._target

    def maxBadEdges(self):
        """
        Returns the maximum allowed number of bad edges.
        """
        return self._maxBadEdges

    def newResult( self, sig, priority, badEdges, source ):
        """
        Records the given sig as a new result, and instructs all processes to
        stop as soon as possible.
        """
        self._results.append( ( sig, priority ) )
        self._details[sig] = ( badEdges, source, priority )
        self._updateMinSig( sig, priority )
        self._stop = True

    def badEdges( self, sig ):
        """
        Returns the tuple of bad edge indices in sig.

        Pre-condition:
        --> The given sig was enqueued at some point in time.
        """
        return self._details[sig][0]

    def enqueue( self, sig, priority, badEdges, source ):
        """
        Pushes the given sig onto the queue with the given priority, assigns
        badEdges as the tuple returned by self.badEdges(sig), and records the
        given source information.
        """
        self._queue.put(
                ( priority, self._total, sig ),
                False )
        self._total += 1
        self._details[sig] = ( badEdges, source, priority )
        self._updateMinSig( sig, priority )

    def dequeue(self):
        """
        Prints some info if a lot of time has passed, and then pops the
        isomorphism signature with the highest priority from the queue.

        Exceptions:
        --> Raises Empty if the queue is empty.
        """
        time = default_timer()
        if time - self._previousInfo > self._interval:
            self._previousInfo = time
            print( self.info() )
            stdout.flush()

        # Try to pop from queue immediately, since we don't want to change
        # anything if it turns out that the queue is currently empty.
        priority, _, sig = self._queue.get(False)
        self._explored += 1
        source = self._details[sig][1]
        self._current = sig
        return ( sig, priority, source )

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


class RemoveBadEdgeManager(SyncManager):
    pass
RemoveBadEdgeManager.register( "RemoveBadEdgeData", RemoveBadEdgeData )


def removeBadEdgeGreedy_(
        tri, sig, priority, source, badEdges, target, inter,
        shared, lock, cond ):
    if source is None:
        steps = 1
    else:
        steps = source[0] + 1

    # Try to remove bad edges using 3-2 moves.
    atTarget = ( priority[0] == target )
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
        if inter:
            newPriority = complexity(
                    newTri, newBadEdges, priority[0] - 1 - target )
        else:
            newPriority = complexity( newTri, newBadEdges )
        newSource = ( steps,
                ( "3-2", i ),
                priority,
                sig )
        relBadEdges = tuple( relabelEdges( newTri, newBadEdges, newIsom ) )
        if atTarget:
            # We started at the target, so removing a bad edge gets us below!
            with lock:
                shared.newResult(
                        newSig, newPriority, relBadEdges, newSource )
            return True

        # We didn't get below the target, so we should continue trying to
        # remove bad edges.
        with lock:
            if shared.queueEmpty():
                shared.enqueue(
                        newSig, newPriority, relBadEdges, newSource )
                cond.notify_all()
            else:
                shared.enqueue(
                        newSig, newPriority, relBadEdges, newSource )
        if removeBadEdgeGreedy_( newIsom(newTri), newSig, newPriority,
                newSource, relBadEdges, target, inter, shared, lock, cond ):
            return True
    return False


def removeBadEdgeExplore_(
        sig, priority, source, inter, shared, lock, cond ):
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
        if removeBadEdgeGreedy_( tri, sig, priority, source,
                badEdges, target, inter, shared, lock, cond ):
            # Terminate immediately if we found a result.
            return
        steps = 1
    else:
        steps = source[0] + 1

    # Try all possible 3-2 moves on *good* edges.
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
        if inter:
            newPriority = complexity(
                    newTri, newBadEdges, priority[0] - target )
        else:
            newPriority = complexity( newTri, newBadEdges )
        newSource = ( steps,
                ( "3-2", edge.index() ),
                priority,
                sig )
        relBadEdges = tuple( relabelEdges(
            newTri, newBadEdges, newIsom ) )
        with lock:
            if shared.queueEmpty():
                shared.enqueue(
                        newSig, newPriority, relBadEdges, newSource )
                cond.notify_all()
            else:
                shared.enqueue(
                        newSig, newPriority, relBadEdges, newSource )

        # Can we remove bad edges?
        if removeBadEdgeGreedy_( newIsom(newTri), newSig, newPriority,
                newSource, relBadEdges, target, inter, shared, lock, cond ):
            # Terminate immediately if we found a result.
            return

    # Now try all 2-3 moves, as long as they don't take us beyond the maximum
    # allowed number of bad edges.
    with lock:
        maxBadEdges = shared.maxBadEdges()
        isBad = shared.getIsBadRoutine()
    ignore = ( priority[0] == maxBadEdges )
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
        newEdge = newTri.edge( newTrans[-1] )
        newEdgeIsBad = isBad(newEdge)
        if ignore and newEdgeIsBad:
            continue

        # We should enqueue the newSig.
        newBadEdges = [ newTrans[i] for i in badEdges ]
        if newEdgeIsBad:
            newBadEdges.append( newTrans[-1] )
        if inter:
            if newEdgeIsBad:
                newPriority = complexity(
                        newTri, newBadEdges, priority[0] + 1 - target )
            else:
                newPriority = complexity(
                        newTri, newBadEdges, priority[0] - target )
        else:
            newPriority = complexity( newTri, newBadEdges )
        newSource = ( steps,
                ( "2-3", triangle.index() ),
                priority,
                sig )
        relBadEdges = tuple( relabelEdges(
            newTri, newBadEdges, newIsom ) )
        with lock:
            if shared.queueEmpty():
                shared.enqueue(
                        newSig, newPriority, relBadEdges, newSource )
                cond.notify_all()
            else:
                shared.enqueue(
                        newSig, newPriority, relBadEdges, newSource )

        # Can we remove bad edges?
        if removeBadEdgeGreedy_( newIsom(newTri), newSig, newPriority,
                newSource, relBadEdges, target, inter, shared, lock, cond ):
            # Terminate immediately if we found a result.
            return


def removeBadEdgeLoop_( shared, lock, cond, inter ):
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
            removed a bad edge).
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
                sig, priority, source = shared.dequeue()
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
        removeBadEdgeExplore_(
                sig, priority, source, inter, shared, lock, cond )


def removeBadEdge(
        s, nProcesses, interval, isBad, badName, excess=0, inter=False,
        badEdges=None ):
    """
    Starting at the given isomorphism signature s, and using the given number
    of parallel processes, greedily searches the Pachner graph for a
    triangulation with fewer bad edges than s.

    We consider an edge e to be bad if and only if isBad(e) returns True.

    The search abides by the following constraints:
    (1) If the initial isomorphism signature s corresponds to a triangulation
        with b bad edges, then we ignore any 2-3 moves that cause the number
        of bad edges to exceed b by more than the given excess.
    (2) Priority is given to triangulations in which the maximum
        "multiplicity defect" (see the mulDefect() routine) among bad edges
        is small, then to triangulations in which the maximum "degree defect"
        (see the degDefect() routine) among bad edges is small, and then to
        triangulations with fewer tetrahedra.
    (3) Triangulations with equal priority are processed in order of
        insertion into a queue.

    If the search successfully finds one or more triangulations with fewer
    bad edges than s, then this routine returns a list containing all of the
    corresponding isomorphism signatures (it is possible for multiple
    parallel processes to each find such a triangulation at roughly the same
    time, which is why the returned list may contain more than one
    isomorphism signature).

    It is also possible for the search to "get stuck" in a situation where
    the only moves available are 2-3 moves that increase the number of bad
    edges beyond the allowed excess. Since we ignore all such moves, this
    routine returns an empty list whenever we get stuck in this way.

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
    --> The given excess is a non-negative integer.
    """
    manager = RemoveBadEdgeManager()
    manager.start()
    lock = manager.Lock()
    cond = manager.Condition(lock)
    shared = manager.RemoveBadEdgeData(
            s, nProcesses, interval, isBad, excess, inter, badEdges )
    processes = [
            Process( target=removeBadEdgeLoop_,
                args=( shared, lock, cond, inter ) )
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
            r, comp, badName, rem ) )
        path = shared.backtrack(r)
        for s in path:
            print( "    {}".format(s) )
    return results
