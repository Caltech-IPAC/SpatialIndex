#!/bin/sh

# The cibuildwheel utility can handle both linux and macos builds but in different
# ways.  For Linux it uses specially constructed Docker linux containers that have all
# the needed Python versions installed (Docker itself must already be configured installed
# on the system) and for MacOS it uses the local Mac compiler (which supports explicit
# deployment targets) but the versions of Python have to be manually installed.
#

# The cibuildwheel utility uses environment variables to control what it actually
# builds.  cibuildwheel in turn uses setup.py to define much of the detailed 
# configuration.

make clean
rm -rf build/*
rm -f wheelhouse/*


# LINUX

export OS='linux'
export CIBW_BUILD='*'
export CIBW_BEFORE_ALL='sh make.sh'
export CIBW_BUILD_FRONTEND='build'

echo "CIBW_BUILD>         " "$CIBW_BUILD"
echo "CIBW_ARCHS>         " "$CIBW_ARCHS"
echo "CIBW_BEFORE_ALL>    " "$CIBW_BEFORE_ALL"
echo "CIBW_BUILD_FRONTEND>" "$CIBW_BUILD_FRONTEND"

pip install pipx

pipx run cibuildwheel --platform $OS


# MACOS

# export OS='macos'
# export CIBW_BUILD='pp311-**'
# export CIBW_SKIP='cp36-* cp37-*'
# export CIBW_ARCHS='x86_64 universal2 arm64'
# export CIBW_BUILD_FRONTEND='build'
# export MACOSX_DEPLOYMENT_TARGET='11.1'
# 
# echo "OS>                      " "$OS"
# echo "CIBW_BUILD>              " "$CIBW_BUILD"
# echo "CIBW_SKIP>               " "$CIBW_SKIP"
# echo "CIBW_ARCHS>              " "$CIBW_ARCHS"
# echo "CIBW_BUILD_FRONTEND>     " "$CIBW_BUILD_FRONTEND"
# echo "MACOSX_DEPLOYMENT_TARGET>" "$MACOSX_DEPLOYMENT_TARGET"
# 
# pip install pipx
# 
# pipx run cibuildwheel --platform $OS
