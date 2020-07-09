# SpatialIndex: Spatial Indexing Library for Astronomy
One of the most common database searches in astronomy involves finding all the  objects in a region of the sky.  DBMS engines provide efficient indexed  searches only for single column (i.e. something that can be sorted linearly). Some DMMSs are starting to support true multidimensional indexing (e.g.,R-Trees) but at present this is neither universal nor uniform.

It is, however, possible to leverage the standard internal, B-tree database indexing to support  efficient 2D indexing.  Tesselating the
sky into a heirarchical, Z-ordered space,  using a scheme like a Heirarchical
Triangular Mesh (HTM)  or Hierarchical Equal Area isoLatitude Pixelization
(HEALPix), places any sky coordinate) in a specific tesselation cell and each cell is identified by a unique ID number)  This number is stored in an integer column in the database table.

To perform a spatial search, the region can be converted to a range of spatial index cell numbers and, because of the Z-ordering,  these ranges tend to be co-located and can therefoe accessed efficiently.  For large tables, this can speed up queries by orders of magnitude. The resultant records are actually a superset of those desired, but a simple geometric filtering, performed as part of the database query, removes the extraneous records.

The spatal searches require that the DBMS table be augmented with four extra columns: a spatial index column at some pre-chosen HTM/HEALPix level and the three-vector (x,y,z) coordinate of the point on the sky.

Our library converts a geometric constraint into a pair of constraints that
can be added to  SQL to turn it into a spatially-indexed region query.  For instance, if the DBMS table has been indexed at HTM level
7, a cone on the sky at latitude 43.7, longitude 129.4 with radius 0.5 degrees
gets turned into the following SQL constraints 

   WHERE (   (htm7 = 245093) 
          OR (htm7 = 245098) 
          OR (htm7 = 245100)
          OR (htm7 = 245105) 
          OR (htm7 = 245110)
          OR (htm7 = 245118))

     AND (-0.45888930755155893*x)+(0.55866098617988125*y)+(0.69088241107685844*z)
         >=9.99975630705394747e-01

These constarints are inserted into the SQL statement that is submitted to the DBMS.

Building the Python library:

The spatial index code need to be wrapped for Python use.  This is taken care of in the Makefile. It consists of using Cython to compile the spatial-index.pyxcode (Python with Cython directives) into "spatial_index.c", then building a LINUX library from this. This results in "spatial_index.cpython-37m-x86_64-linux-gnu.so". which can be loaded by Python at runtime. Finally, we tuen this  into a Wheel file that can be pip-installed within
our Python distribution, or uploaded to PyPI.  This file is
(dist/spatial_index-0.9-cp37-cp37m-linux_x86_64.whl).   The primary files created in this sequence are therefore:

  spatial_index.c
  spatial_index.cpython-37m-x86_64-linux-gnu.so
  dist/spatial_index-0.9-cp37-cp37m-linux_x86_64.whl
