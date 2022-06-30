"""
Given a triangulation with no "bad" edges (where "bad" stands for any
property of an edge that can be computed algorithmically), search the Pachner
graph for a smaller example of such a triangulation. Since 2-3 moves create
new edges in the triangulation, we only consider such moves when the new edge
is not bad.
"""
from sys import stdout
from timeit import default_timer
from multiprocessing import Condition, Process
from multiprocessing.managers import SyncManager
from queue import Empty
from regina import *
from search import SearchData
from pachner import twoThree, threeTwo


class SimplifyData(SearchData):
    """
    Stores all the data that needs to be shared between processes when
    running the simplifyNoBadEdge() routine.
    """
    def __init__( self, s, nProcesses, interval, isBad, height ):
        """
        Initialises shared search data.
        
        Among other things, enqueues the initial isomorphism signature s, and
        provides the isBad() routine for processes to use to test whether an
        edge is bad.
        """
        super().__init__( s, nProcesses, interval, isBad )

        # Initial triangulation.
        t = Triangulation3.fromIsoSig(s)
        self._current = s

        # Our target is to reduce the number of tetrahedra without introducing
        # bad edges.
        self._target = t.size()
        self._maxSize = self._target + height

        # Enqueue the initial isomorphism signature.
        self.enqueue( s,
                self._target, # Prioritise by size.
                None ) # Initial iso sig has no source.

    def info(self):
        """
        Returns a string describing the time elapsed, the number of
        triangulations discovered, the number of triangulations explored, and
        the latest isomorphism signature to be dequeued.
        """
        return super().info() + " Current: {}.".format( self._current )

    def target(self):
        """
        Returns the target number of tetrahedra that we need to get below.
        """
        return self._target

    def maxSize(self):
        """
        Returns the maximum allowed number of tetrahedra.
        """
        return self._maxSize

    def newResult( self, sig, source ):
        """
        Records the given sig as a new result, and instructs all processes to
        stop as soon as possible.
        """
        self._results.append(sig)
        self._details[sig] = source
        self._stop = True

    def enqueue( self, sig, priority, source ):
        """
        Pushes the given sig onto the queue with the given priority, and
        records the given source information.
        """
        #TODO
        # We prioritise by size, and 3-2 moves will be performed greedily
        # without ever needing to first dequeue the given sig, so if the
        # given priority is already equal to the maxSize then there is no
        # need to actually enqueue.
        if ( source is None ) or ( priority < self._maxSize ):
            self._queue.put(
                    ( priority, self._total, sig ),
                    False )
        self._total += 1
        self._details[sig] = source

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
        source = self._details[sig]
        self._current = sig
        return ( sig, priority, source )

    def backtrack( self, sig ):
        """
        Returns a list that describes the sequence of 2-3 and 3-2 moves that
        produces the given sig.
        """
        back = []
        source = self._details[sig]
        while source is not None:
            back.append(source)
            source = self._details[ source[2] ]
        return back[::-1]


class SimplifyManager(SyncManager):
    pass
SimplifyManager.register( "SimplifyData", SimplifyData )


def simplifyGreedy_( sig, priority, source, shared, lock, cond ):
    tri = Triangulation3.fromIsoSig(sig)
    tetCount = tri.size()
    with lock:
        atTarget = ( shared.target() == tetCount )
    if source is None:
        steps = 1
    else:
        steps = source[0] + 1

    # Try all possible 3-2 moves.
    newPriority = tetCount - 1 # Prioritise by size.
    for edge in tri.edges():
        result = threeTwo(edge)
        if result is None:
            continue
        newSig = result[0].isoSig()
        with lock:
            if shared.alreadySeen(newSig):
                continue
        newSource = ( steps, ( "3-2", edge.index() ), sig )

        # Did we get below the target? If we did, then terminate immediately.
        if atTarget:
            # We started at the target, so a 3-2 move gets us below!
            with lock:
                shared.newResult( newSig, newSource )
            return True

        # We didn't get below the target. Enqueue newSig so that we can
        # explore further later on.
        with lock:
            if shared.queueEmpty():
                shared.enqueue( newSig, newPriority, newSource )
                cond.notify_all()
            else:
                shared.enqueue( newSig, newPriority, newSource )

        # Greedily attempt 3-2 moves.
        if simplifyGreedy_(
                newSig, newPriority, newSource, shared, lock, cond ):
            return True
    return False


def simplifyExplore_( sig, priority, source, shared, lock, cond ):
    # NOTES
    """
    --> Many of the computations can be done independently of other
        processes, but be careful to acquire the given lock every time we
        need to access shared.
    """
    # If sig was generated by the search (i.e., source is not None), then
    # we would already have greedily attempted all possible 3-2 moves.
    if source is None:
        if simplifyGreedy_( sig, priority, source, shared, lock, cond ):
            # Terminate immediately if we found a result.
            return
        steps = 1
    else:
        steps = source[0] + 1

    # Try all possible 2-3 moves, as long as they don't introduce any new bad
    # edges.
    # NOTE: We rely on the fact that shared.enqueue() will take care of
    #   ensuring that we will never exceed the maxSize.
    tri = Triangulation3.fromIsoSig(sig)
    tetCount = tri.size()
    with lock:
        if tetCount == shared.maxSize():
            return
        isBad = shared.getIsBadRoutine()
    newPriority = tetCount + 1 # Prioritise by size.
    for triangle in tri.triangles():
        result = twoThree(triangle)
        if result is None:
            continue
        newTri, newTrans = result
        newSig = newTri.isoSig()

        # Have we already seen newSig?
        with lock:
            if shared.alreadySeen(newSig):
                continue

        # Did we introduce a bad edge?
        if isBad( newTri.edge( newTrans[-1] ) ):
            continue

        # Enqueue newSig so that we can explore further later on.
        newSource = ( steps,
                ( "2-3", triangle.index() ),
                sig )
        with lock:
            if shared.queueEmpty():
                shared.enqueue( newSig, newPriority, newSource )
                cond.notify_all()
            else:
                shared.enqueue( newSig, newPriority, newSource )

        # Greedily attempt 3-2 moves.
        if simplifyGreedy_(
                newSig, newPriority, newSource, shared, lock, cond ):
            # Terminate immediately if we found a result.
            return


def simplifyLoop_( shared, lock, cond ):
    # NOTES
    """
    --> *Wait* (but don't stop) if the shared queue is empty but there is
        still at least one other process running (since new sigs could still
        be added to the queue).
    --> Two reasons to stop:
        (1) We have been instructed to stop (because we have successfully
            reduced the number of tetrahedra).
        (2) The shared queue is empty and no other processes are running.
    """
    while True:
        # Keep explorings sigs until we are instructed to stop.
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
        simplifyExplore_(
                sig, priority, source, shared, lock, cond )


def simplifyNoBadEdge(
        s, nProcesses, interval, isBad, height=0 ):
    """
    Starting at the given isomorphism signature s, and using the given number
    of parallel processes, attempts to reduce the number of tetrahedra
    without ever introducing a new bad edge, where we consider an edge e to
    be bad if and only if isBad(e) returns True.

    If the search successfully reduces the number of tetrahedra, then this
    routine returns a list containing the isomorphism signatures of the
    simplified triangulations (it is possible for multiple parallel processes
    to each find such a triangulation at roughly the same time, in which case
    the returned list will contain more than one isomorphism signature).

    It is also possible for the search to "get stuck" in a situation where
    the only moves available are 2-3 moves that either introduce a new bad
    edge or increase the number of additional tetrahedra above the given
    height. Since we ignore all such moves, this routine returns an empty
    list whenever we get stuck in this way.

    Note: When run with multiple processes (i.e., nProcesses > 1), this
    routine is not completely deterministic because isomorphism signatures
    are partly prioritised by insertion order, and this order may vary across
    different run-throughs.

    The given interval determines the number of seconds that we should wait
    before printing some info to standard output.

    Pre-condition:
    --> The function isBad takes a single edge e as input, and outputs either
        True or False (corresponding to whether the edge e should be
        considered bad).
    --> The given height is a non-negative integer.
    """
    manager = SimplifyManager()
    manager.start()
    lock = manager.Lock()
    cond = manager.Condition(lock)
    shared = manager.SimplifyData(
            s, nProcesses, interval, isBad, height )
    processes = [
            Process( target=simplifyLoop_,
                args=( shared, lock, cond ) )
            for _ in range(nProcesses) ]
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    # Print some info, and then return results.
    print( shared.info() )
    results = shared.getResults()
    for r in results:
        print( "Result: {}.".format(r) )
        path = shared.backtrack(r)
        for s in path:
            print( "    {}".format(s) )
    return results
