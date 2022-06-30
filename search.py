"""
Given a one-vertex 3-manifold triangulation, perform a parallelised search
through the Pachner graph.
"""
from timeit import default_timer
from queue import PriorityQueue


class SearchData:
    """
    A base class for storing the data that needs to be shared between
    processes when performing a parallelised search through a Pachner graph.

    Subclasses must provide the following routines, as they are not
    implemented in this class:
    --> newResult()
    --> enqueue()
    --> dequeue()
    --> getSource()
    """
    def __init__( self, s, nProcesses, interval, isBad ):
        """
        Initialises the shared search data.

        Among other things, enqueues the initial isomorphism signature s, and
        provides the isBad() routine for processes to use to test whether an
        edge is bad.

        Subclasses must call enqueue() as the very last step of their
        __init__() routines.
        """
        # Progress tracking.
        self._start = default_timer()
        self._previousInfo = self._start
        self._interval = interval

        # Number of processes currently running.
        self._nRunning = nProcesses

        # Keep track of relevant details of all isomorphism signatures that
        # we have encountered.
        self._total = 0
        self._explored = 0
        self._details = dict()

        # Store any results we find, and track whether we should stop the
        # computation early.
        self._stop = False
        self._results = []

        # Share isBad routine with all processes.
        self._isBad = isBad

        # Priority queue of isomorphism signatures to explore.
        self._queue = PriorityQueue()

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

    def info(self):
        """
        Returns a string describing the time elapsed, the number of
        triangulations encountered, and the number of triangulations
        explored.

        Subclasses should probably override this routine to provide more
        specialised information.
        """
        return "Time: {:.6f}. Seen: {}. Visited: {}.".format(
                self.totalTime(), self._total, self._explored )

    def nRunning(self):
        """
        Returns the number of processes that are currently running (i.e., not
        paused).
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

    def getResults(self):
        """
        Returns the (possibly empty) list of results.
        """
        return self._results

    def newResult( self, sig ):
        """
        Records the given sig a new result.

        Subclasses must provide their own implementation for this routine.
        This routine may be a good place to stop the computation early.
        """
        raise NotImplementedError()

    def stop(self):
        """
        Checks whether the computation should be stopped early.
        """
        return self._stop

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

    def enqueue( self, sig, priority, source ):
        """
        Pushes the given sig onto the queue with the given priority, and
        records the given source information.

        Subclasses must provide their own implementation for this routine.
        """
        raise NotImplementedError()

    def dequeue(self):
        """
        Pops the isomorphism signature with the highest priority from the
        queue.

        Subclasses must provide their own implementation for this routine,
        and must raise Empty if the queue is empty. This routine may be a
        good place to print info if a lot of time has passed.
        """
        raise NotImplementedError()

    def getSource( self, sig ):
        """
        Returns a tuple describing the move that produced the given sig.

        Subclasses must provide their own implementation for this routine.
        """
        raise NotImplementedError()

    def backtrack( self, sig ):
        """
        Returns a list that describes the sequence of 2-3 and 3-2 moves that
        produced the given sig.
        """
        back = []
        source = self.getSource(sig)
