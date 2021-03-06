# SpatialIndex

## Spatial Indexing Library for Astronomy

One of the most common database searches in astronomy involves finding all the 
objects in a region of the sky.  Basic DBMS engines provide efficient indexed 
searches only for single column (i.e. something that can be sorted linearly).
Some servers are starting to support true multidimensional indexing (e.g.,
R-Trees) but this is neither universal nor uniform.

It is, however, possible to leverage the standard internal database indexing
(B-Tree) to support impressively efficient 2D indexing.  If one tesselates the
sky into a hierarchical, Z-ordered space (using a scheme like a Heirarchical
Triangular Mesh (HTM)  or Hierarchical Equal Area isoLatitude Pixelization
(HEALPix)), then any "object" (sky coordinate) belongs to a specific 
tesselation cell and each cell has a unique ID (number).  This number can
be stored in an integer column in the database table.

Then when one wants to perform a spatial search, the region can be converted
to a range (actually set of ranges) of spatial index cell numbers and
(because of the Z-ordering) these ranges tend to be co-located and can be
accessed efficiently.  For large tables, this can speed up queries by orders
of magnitude.

To perform a spatial search, the region can be converted to a range of spatial
index cell numbers and, because of the Z-ordering,  these ranges tend to be 
co-located and can therefore accessed efficiently.  For large tables, this can
speed up queries by orders of magnitude. The resultant records are actually a
superset of those desired, but a simple geometric filtering, performed as part
of the database query, removes the extraneous records.

The spatal searches require that the DBMS table be augmented with four extra
columns: a spatial index column at some pre-chosen HTM/HEALPix level and the
three-vector (x,y,z) coordinate of the point on the sky.

Our library converts a geometric constraint into a pair of constraints that
can be added to  SQL to turn it into a spatially-indexed region query.  For
instance, if the DBMS table has been indexed at HTM level 7, a cone on the
sky at latitude 43.7, longitude 129.4 with radius 0.5 degrees gets turned
into the following SQL constraints 

>   WHERE (   (htm7 = 245093) 
>          OR (htm7 = 245098) 
>          OR (htm7 = 245100)
>          OR (htm7 = 245105) 
>          OR (htm7 = 245110)
>          OR (htm7 = 245118))
>
>     AND (-0.45888930755155893*x)+(0.55866098617988125*y)+(0.69088241107685844*z)
>         >=9.99975630705394747e-01

These constraints are inserted into the SQL statement  submitted to the DBMS.

Building the Python library:

The C spatial index code needs to be wrapped for Python use.  This is taken care
of in the Makefile and consists of using Cython to compile the spatial-index.pyx
code (Python with Cython directives) into "spatial_index.c", then building 
a LINUX library from the result plus our C libraries (our C code plus tinyhtm).
This results in a library file "spatial_index.cpython-37m-x86_64-linux-gnu.so"
which has the right content to be loadable by the Python runtime and accessed by 
Python calls.

The last step is to turn this into a Wheel file that can be pip-install in
our Python distribution and/or uploaded to PyPI.  This file is
(dist/spatial_index-0.9-cp37-cp37m-linux_x86_64.whl).  There is bookkeepping
as well, but the primary files in this sequence are again

  spatial_index.c
  spatial_index.cpython-37m-x86_64-linux-gnu.so
  dist/spatial_index-0.9-cp37-cp37m-linux_x86_64.whl

*The NASA Exoplanet Science Institute is operated by the California Institute
of Technology, under contract with the National Aeronautics and Space Administration
under the Exoplanet Exploration Program.*
