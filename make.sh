#!/bin/sh

mkdir -p SpatialIndex
mkdir -p lib


# Get and build SpatialIndex

git clone https://github.com/Caltech-IPAC/SpatialIndex.git 

(cd SpatialIndex && make libs)
find . -name \*.o -print


# Copy the files from the SpatialIndex build that we need to build the wheel

cp -r SpatialIndex/lib/* lib 

cp    SpatialIndex/pyproject.toml .
cp    SpatialIndex/README.md      .   
cp    SpatialIndex/LICENSE        .   
