PolyKnot
========

A C-library to heuristically identify and analyse 'knots' in open chains such as
polymers. On receipt of the positions of a beads along the chain or string the
function finds: crossing number, knot size and knot position.

Usage
-----

* Create knot struct.
* Vector inputted needs to describe the positions of the 'joints' or 'beads' in the
chain with length 3N like so: **r** = ( x1, y1, z1, x2, y2, z2 ... xN, yN, zN )
* Run analysis on knot object and chosen vector.
* Clear knot struct.

Theory
------

PolyKnot 'unties' an arbitrary length chain via Reidemeister moves rather than
attempting to classify the knot through symbolic representation of the
crossings and subsequent simplification and recognition of a polynomial. The
benefit of this, aside from being more a efficient routine in terms of memory and
CPU expense, is the possibility of identifying a knot's 'position' and 'size'
along a open chain without forcing the ends together.

To Do
-----

* More efficient exploration of Reidemeister 3 moves.
* Find best combination of projections.
* Memory optimise.
* categorise knots with cross count (via recognition)
* Two knots?
