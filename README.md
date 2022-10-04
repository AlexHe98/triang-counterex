Counterexamples for triangulations
----------------------------------

This repository contains source code for finding counterexamples for
one-vertex 3-manifold triangulations by using Pachner moves to remove
distinguished edges.

The code is written to run using Regina's regina-python command-line program.
For details on how to install Regina, visit https://regina-normal.github.io/.
See below for usage notes.

We currently recognise three types of distinguished edges:
* **Core edges**: Deleting a solid torus neighbourhood of such an edge leaves
    behind another solid torus. The existence of a core edge implies that the
    triangulation represents a lens space.
* **Tunnel edges**: In an ideal triangulation of a knot complement, deleting
    a small regular neighbourhood of such an edge turns the knot complement
    into a genus-2 orientable handlebody. The existence of a tunnel edge
    implies that the knot is either trivial or has tunnel number one.
* **Edges isotopic to Seifert fibres**

The main scripts in this repository are:
* ``removeCore.py``: Attempts to remove a core edge.
* ``removeTunnel.py``: Attempts to remove a tunnel edge.
* ``removeFibre.py``: Attempts to remove an edge isotopic to a Seifert fibre.
* ``simplifyNoCore.py``: Attempts to simplify a triangulation without
    introducing new core edges.
* ``simplifyNoTunnel.py``: Attempts to simplify a triangulation without
    introducing new tunnel edges.
* ``simplifyNoFibre.py``: Attempts to simplify a triangulation without
    introducing new edges isotopic to Seifert fibres.

The bulk of the supporting code can be found in ``search.py``,
``removeBadEdge.py`` and ``simplify.py``.

The script ``removeCore.py`` takes three command-line arguments:
* An isomorphism signature for the triangulation from which we would like to
    remove a core edge.
* The number of concurrent processes to use.
* The number of seconds to wait before printing a brief progress update.

In the following example, assume that the regina-python command-line program
is in the directory ``dir``. To use 4 concurrent processes to remove a core
edge from the triangulation with isomorphism signature "cMcabbgqs", and to
request a progress update every 10 seconds, run the following:
> ``dir/regina-python removeCore.py cMcabbgqs 4 10``

The scripts ``removeTunnel.py`` and ``removeFibre.py`` are similar, except
that ``removeFibre.py`` takes one additional mandatory command-line argument
and one additional optional command-line argument. For details, see the
documentation in ``removeFibre.py``.

The scripts ``simplifyNoCore.py``, ``simplifyNoTunnel.py`` and
``simplifyNoFibre.py`` each take four command-line arguments:
* An isomorphism signature.
* The number of concurrent processes to use.
* The number of seconds to wait before printing a brief progress update.
* The allowed excess height (i.e., the amount by which we are allowed to
    increase the size of the input triangulation).

— *Alex He (a.he@uqconnect.edu.au)*
