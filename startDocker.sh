#!/bin/bash

DOCKCROSS_IMAGE=dockcross/manylinux2010-x64

docker run -i -t \
   -v $PWD:/build \
   $DOCKCROSS_IMAGE bash
