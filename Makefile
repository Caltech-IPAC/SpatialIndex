LIB_DIR = lib

all:				libs spatial_index pgm

libs:
					(cd lib; make; make install)

spatial_index:	setup.py spatial_index.pyx
					python setup.py bdist_wheel
					pip install -t ../lib --upgrade dist/*.whl

pgm:				
					(cd src; make; make install)

clean:
					(cd lib; make clean)
					(cd src; make clean)
					rm -rf build
					rm -rf spatial_index.egg-info
					rm spatial_index.c
					rm -f make.out
