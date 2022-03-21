"""
Given a triangulation with at least one "bad" edge (where "bad" stands for
any property of an edge that can be computed algorithmically), search for a
3-2 move that removes a bad edge, prioritising triangulations in which:
(1) the maximum "multiplicity defect" among bad edges is small (see the
    mulDefect() routine);
(2) the maximum "degree defect" among bad edges is small (see the degDefect()
    routine); and
(3) the number of tetrahedra is small.

This is written to run with either Python 2 or Python 3, without needing to
rewrite any code.
"""
from __future__ import print_function
from sys import stdout
from timeit import default_timer
from multiprocessing import Condition, Process
from multiprocessing.managers import SyncManager
try:
    from queue import Empty, PriorityQueue
except ImportError:
    from Queue import Empty, PriorityQueue
from regina import *
from pachner import twoThree, threeTwo
from triangCounterexHelpers import mulDefect
from triangCounterexHelpers import degDefect
from triangCounterexHelpers import relabelEdges


def complexity( tri, badEdges, excess=None ):
    """
    Returns a measure for how far away we are from removing the given
    collection of bad edges from the given triangulation.

    If excess is None (the default), then the returned complexity is a tuple
    consisting of the following items:
    (0) The number of bad edges (i.e., the length of the given collection of
        bad edge indices).
    (1) The maximum "multiplicity defect" among the given bad edges (see the
        mulDefect() routine).
    (2) The maximum "degree defect" among the given bad edges (see the
        degDefect() routine).
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


class RemoveBadEdgeData:
    """
    Stores all the data that needs to be shared between processes when
    running the removeBadEdge() routine.
    """
    def __init__( self, s, nProcesses, interval, isBad, excess, useMin ):
        """
        Initialises shared search data.

        Among other things, enqueues the initial isomorphism signature s.
        """
        # Progress tracking.
        self._start = default_timer()
        self._previousInfo = self._start
        self._interval = interval

        # Number of processes currently running.
        self._nRunning = nProcesses

        # Details of all discovered isomorphism signatures.
        self._total = 0
        self._explored = 0
        self._details = dict()

        # Initial triangulation.
        t = Triangulation3.fromIsoSig(s)
        b = [ e.index() for e in t.edges() if isBad(e) ]
        if useMin:
            c = complexity( t, b, 0 )
        else:
            c = complexity( t, b )
        self._minSig = { c[0]: s }
        self._current = s

        # We can stop the computation as soon as we find a result.
        self._stop = False
        self._target = c[0]
        self._maxBadEdges = self._target + excess
        self._results = []

        # Share isBad routine with all processes.
        self._isBad = isBad

        # Priority queue of isomorphism signatures to explore.
        # Note that everything else needs to have been initialised before
        # calling enqueue().
        self._queue = PriorityQueue()
        self.enqueue( s, c, tuple(b),
                None ) # Initial sig has no source.

    def totalTime(self):
        """
        Returns the time elapsed since self was initialised.
        """
        return default_timer() - self._start

    def getIsBadRoutine(self):
        """
        Returns a routine that takes a single edge as input and determines
        whether this edge is bad.
        """
        return self._isBad

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
        msg = ""
        for badEdgeCount in self._minSig:
            msg += minFormat.format(
                    badEdgeCount,
                    self._minComplexity(badEdgeCount),
                    self._minSig[badEdgeCount] )
        msg += " Current: {}, {}.".format(
                self._details[ self._current ][2],
                self._current )
        return "Time: {:.6f}. Seen: {}. Visited: {}.{}".format(
                self.totalTime(), self._total, self._explored, msg )

    def nRunning(self):
        """
        Returns the number of processes that are currently running (i.e.,
        not paused).
        """
        return self._nRunning

    def tellPaused(self):
        """
        Tells self that a process has paused.
        """
        self._nRunning -= 1

    def tellResumed(self):
        """
        Tells self that a process has resumed.
        """
        self._nRunning += 1

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

    def getResults(self):
        """
        Returns the (possibly empty) list of results.
        """
        return self._results

    def newResult( self, sig, priority, badEdges, source ):
        """
        Records the given sig as a new result, and instructs all processes
        to stop as soon as possible.
        """
        self._results.append( ( sig, priority ) )
        self._details[sig] = ( badEdges, source, priority )
        self._updateMinSig( sig, priority )
        self._stop = True

    def stop(self):
        """
        Checks whether the computation should be stopped because
        self.getResults() is no longer empty.
        """
        return self._stop

    def badEdges( self, sig ):
        """
        Returns the tuple of bad edge indices in sig.

        Pre-condition:
        --- The given sig was enqueued at some point in time.
        """
        return self._details[sig][0]

    def alreadySeen( self, sig ):
        """
        Returns True if and only if the given sig has already been seen.
        """
        return sig in self._details

    def queueEmpty(self):
        """
        Returns True if and only if the queue is currently empty.
        """
        return self._queue.empty()

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
        Prints some info is a lot of time has passed, and then pops the
        isomorphism signature with the highest priority from the queue.

        Exceptions:
        --- Raises Empty if the queue is empty.
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


def removeBadEdgeGreedy(
        tri, sig, priority, source, badEdges, target, useMin,
        shared, lock, cond ):
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
        if useMin:
            newPriority = complexity(
                    newTri, newBadEdges, priority[0] - 1 - target )
        else:
            newPriority = complexity( newTri, newBadEdges )
        newSource = ( source[0] + 1,
                ( "3-2", i ),
                priority,
                sig )
        relBadEdges = tuple( relabelEdges( newTri, newBadEdges, newIsom ) )
        if priority[0] == target:
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
        newIsom.applyInPlace(newTri)
        if removeBadEdgeGreedy( newTri, newSig, newPriority, newSource,
                relBadEdges, target, useMin, shared, lock, cond ):
            return True
    return False


def removeBadEdgeExplore(
        sig, priority, source, useMin, shared, lock, cond ):
    # NOTES
    """
    --- Many of the computations can be done independently of other
        processes, but be careful to acquire the given lock every time we
        need to access shared.
    """
    tri = Triangulation3.fromIsoSig(sig)
    tetCount = tri.size()
    with lock:
        badEdges = shared.badEdges(sig)
        target = shared.target()

    # If one of the bad edges could be removed, then we would already have
    # done so. Thus, we only need to try to remove *good* edges.
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
        if useMin:
            newPriority = complexity(
                    newTri, newBadEdges, priority[0] - target )
        else:
            newPriority = complexity( newTri, newBadEdges )
        if source is None:
            steps = 1
        else:
            steps = source[0] + 1
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
        newIsom.applyInPlace(newTri)
        if removeBadEdgeGreedy( newTri, newSig, newPriority, newSource,
                relBadEdges, target, useMin, shared, lock, cond ):
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
        if useMin:
            if newEdgeIsBad:
                newPriority = complexity(
                        newTri, newBadEdges, priority[0] + 1 - target )
            else:
                newPriority = complexity(
                        newTri, newBadEdges, priority[0] - target )
        else:
            newPriority = complexity( newTri, newBadEdges )
        if source is None:
            steps = 1
        else:
            steps = source[0] + 1
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
        newIsom.applyInPlace(newTri)
        if removeBadEdgeGreedy( newTri, newSig, newPriority, newSource,
                relBadEdges, target, useMin, shared, lock, cond ):
            # Terminate immediately if we found a result:
            return


def removeBadEdgeLoop( shared, lock, cond, useMin ):
    # NOTES
    """
    --- To minimise testing whether edges are bad (which could potentially be
        very expensive), take advantage of transitions to keep track of which
        edges are bad.
    --- *Wait* (but don't stop) if the shared queue is empty but there is
        still at least one other process running (since new sigs could still
        be added to the queue).
    --- Two reasons to stop:
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
        removeBadEdgeExplore(
                sig, priority, source, useMin, shared, lock, cond )


def removeBadEdge(
        s, nProcesses, interval, isBad, badName, excess=0, useMin=False ):
    """
    Starting at the given isomorphism signature s, and using the given number
    of parallel processes, searches for a 3-2 move that reduces the number of
    bad edges in s, where we consider an edge e to be bad if and only if
    isBad(e) returns True.

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
        insertion.

    If the search successfully finds one or more 3-2 moves that reduce the
    number of bad edges, then this routine returns a list containing the
    isomorphism signatures of the triangulations that result from performing
    these moves (it is possible for multiple parallel processes to each find
    such a move at roughly the same time, in which case the returned list
    will contain more than one isomorphism signature).

    It is also possible for the search to "get stuck" in a situation where
    the only moves available are 2-3 moves increase the number of bad edges
    beyond the allowed excess. Since we ignore all such moves, this routine
    returns an empty list whenever we get stuck in this way.

    Note: When run with multiple processes (i.e., nProcesses > 1), this
    routine is not completely deterministic because insertion order may vary
    across different run-throughs.

    The given interval determines the number of seconds that we should wait
    before printing some info to standard output.

    Pre-condition:
    --- The triangulation represented by s has at least one bad edge.
    --- The function isBad takes a single edge e as input, and outputs either
        True or False (corresponding to whether the edge e should be
        considered bad).
    --- The given excess is a non-negative integer.
    """
    manager = RemoveBadEdgeManager()
    manager.start()
    lock = manager.Lock()
    cond = manager.Condition(lock)
    shared = manager.RemoveBadEdgeData(
            s, nProcesses, interval, isBad, excess, useMin )
    processes = [
            Process( target=removeBadEdgeLoop,
                args=( shared, lock, cond, useMin ) )
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

