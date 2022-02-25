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
from __future__ import print_function
from sys import stdout
from timeit import default_timer
from multiprocessing import Condition
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
from triangCounterexHelpers import pinchGivesHandlebody


BADNAME = "Tunnels"


def isBad( edge ):
    return pinchGivesHandlebody(edge, 2)


def complexity( tri, badEdges ):
    mulDefs = [ mulDefect( tri.edge(i) ) for i in badEdges ]
    degDefs = [ degDefect( tri.edge(i) ) for i in badEdges ]
    mulDefs.append(0)
    degDefs.append(0)
    return ( len(badEdges),
             max(mulDefs),
             max(degDefs),
             tri.size() )


class RemoveBadEdgeData:
    """
    Stores all the data that needs to be shared between processes when
    running the removeBadEdge() routine.
    """
    def __init__( self, s, nProcesses, interval ):
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

        # We can stop the computation as soon as we find a result.
        self._stop = False
        self._results = []

        # Details of all discovered isomorphism signatures.
        self._total = 0
        self._explored = 0
        self._details = dict()

        # Initial triangulation.
        self._minSig = s
        t = Triangulation3.fromIsoSig(s)
        b = [ e.index() for e in t.edges() if isBad(e) ]
        c = complexity( t, b )

        # Priority queue of isomorphism signatures to explore.
        # Note that everything else needs to have been initialised before
        # calling enqueue()
        self._queue = PriorityQueue()
        self.enqueue( s, c, tuple(b),
                None ) # Initial sig has no source.

    def totalTime(self):
        """
        Returns the time elapsed since self was initialised.
        """
        return default_timer() - self._start

    def _minComplexity(self):
        return self._details[ self._minSig ][2]

    def info(self):
        """
        Returns a string describing the time elapsed, the number of
        triangulations discovered, the number of triangulations explored, and
        details of the triangulation with the smallest complexity that has
        been discovered so far.
        """
        return "Time: {:.6f}. Seen: {}. Visited: {}. Min: {}, {}.".format(
                self.totalTime(), self._total, self._explored,
                self._minComplexity(), self._minSig )

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
        if priority < self._minComplexity():
            self._minSig = sig

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


def removeBadEdgeExplore(
        sig, oldPriority, oldSource, shared, lock, cond ):
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

    # If one of the bad edges could be removed, then we would already have
    # done so. Thus, we only need to try to remove *good* (i.e., non-bad)
    # edges using 3-2 moves.
    for edge in tri.edges():
        if edge.index() in badEdges:
            continue
        result = threeTwo(edge)
        if result is None:
            continue
        newTri, transition = result
        newSig, isom = newTri.isoSigDetail()

        # Have we already seen newSig?
        with lock:
            if shared.alreadySeen(newSig):
                continue

        # Priority and source.
        newBadEdges = [ transition[i] for i in badEdges ]
        priority = complexity( newTri, newBadEdges )
        if oldSource is None:
            steps = 1
        else:
            steps = oldSource[0] + 1
        newSource = ( steps,
                ( "3-2", edge.index() ),
                oldPriority,
                sig )

        # We now know that newSig is genuinely new. Check whether we can
        # remove any bad edges using a 3-2 moves.
        relBadEdges = tuple( relabelEdges(
            newTri, newBadEdges, isom ) )
        isom.applyInPlace(newTri)
        foundResult = False
        for i in relBadEdges:
            result = threeTwo( newTri.edge(i) )
            if result is None:
                continue

            # We have found a result!
            outTri, outTrans = result
            outSig, outIsom = outTri.isoSigDetail()
            outBadEdges = [ outTrans[j] for j in relBadEdges ]
            outBadEdges.remove(-1)
            outPriority = complexity( outTri, outBadEdges )
            outSource = ( newSource[0] + 1,
                    ( "3-2", i ),
                    priority,
                    newSig )
            relOutBadEdges = tuple( relabelEdges(
                outTri, outBadEdges, outIsom ) )
            with lock:
                if shared.alreadySeen(outSig):
                    continue
                shared.newResult(
                        outSig, outPriority, relOutBadEdges, outSource )
            foundResult = True
        if foundResult:
            return

        # No new result yet, so we need to enqueue newSig.
        # If shared.queueEmpty(), then don't forget to notify any paused
        # processes to resume.
        with lock:
            if shared.queueEmpty():
                shared.enqueue(
                        newSig, priority, relBadEdges, newSource )
                cond.notify_all()
            else:
                shared.enqueue(
                        newSig, priority, relBadEdges, newSource )

    # Now try all 2-3 moves that introduce a good edge.
    for triangle in tri.triangles():
        result = twoThree(triangle)
        if result is None:
            continue
        newTri, transition = result
        newSig, isom = newTri.isoSigDetail()

        # Have we already seen newSig?
        with lock:
            if shared.alreadySeen(newSig):
                continue

        # Did we introduce a bad edge?
        # This computation could potentially be very expensive, so we should
        # make sure to parallelise.
        newEdge = newTri.edge( transition[-1] )
        if isBad(newEdge):
            continue

        # Priority and source.
        newBadEdges = [ transition[i] for i in badEdges ]
        priority = complexity( newTri, newBadEdges )
        if oldSource is None:
            steps = 1
        else:
            steps = oldSource[0] + 1
        newSource = ( steps,
                ( "2-3", triangle.index() ),
                oldPriority,
                sig )

        # We now know that newSig is genuinely new, and that the newly
        # introduced edge is good. Check whether we can remove any bad edges
        # using a 3-2 move.
        relBadEdges = tuple( relabelEdges(
            newTri, newBadEdges, isom ) )
        isom.applyInPlace(newTri)
        foundResult = False
        for i in relBadEdges:
            result = threeTwo( newTri.edge(i) )
            if result is None:
                continue

            # We have found a new result!
            outTri, outTrans = result
            outSig, outIsom = outTri.isoSigDetail()
            outBadEdges = [ outTrans[j] for j in relBadEdges ]
            outBadEdges.remove(-1)
            outPriority = complexity( outTri, outBadEdges )
            outSource = ( newSource[0] + 1,
                    ( "3-2", i ),
                    priority,
                    newSig )
            relOutBadEdges = tuple( relabelEdges(
                outTri, outBadEdges, outIsom ) )
            with lock:
                if shared.alreadySeen(outSig):
                    continue
                shared.newResult(
                        outSig, outPriority, relOutBadEdges, outSource )
        if foundResult:
            return

        # No new result yet, so we need to enqueue newSig.
        # If shared.queueEmpty(), then don't forget to notify any paused
        # processes to resume.
        with lock:
            if shared.queueEmpty():
                shared.enqueue(
                        newSig, priority, relBadEdges, newSource )
                cond.notify_all()
            else:
                shared.enqueue(
                        newSig, priority, relBadEdges, newSource )


def removeBadEdgeLoop( shared, lock, cond ):
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
        removeBadEdgeExplore( sig, priority, source, shared, lock, cond )

