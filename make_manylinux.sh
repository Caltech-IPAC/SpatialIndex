# /bin/sh

# This script is specific to building the 'manylinux' wheels
# on a Centos06 platform (in our case under Docker).  Things
# like python2/python3 and auditwheel are already installed
# in the Docker image we ae using.
#
# In the Docker container, Python 2.7 is installed in
#
#    /opt/python/cp27-cp27m/bin
#
# Python 3.6, 3.7 and 3.8 are installed in
#
#    /opt/python/cp36-cp36m/bin
#    /opt/python/cp37-cp37m/bin
#    /opt/python/cp38-cp38/bin
#
# All have 'python' and 'pip' but only 3.7 has 'auditwheel'.
# To avoid confusion, we'll use full paths.


# ------------------------------------------------------------
# IMPORTANT:  This only works in you are in the right Docker 
# container so first start Docker (using startDocker.sh) then 
# you will have to navigate back to this directory from inside
# the Docker container ("cd /build").
# ------------------------------------------------------------


rm -rf final_dist
rm -rf final_wheel

mkdir final_dist
mkdir final_wheel


# To be safe, we will recompile the C code for this
# (Centos06) platform

make clean
make





# Python 2.7

/opt/python/cp27-cp27m/bin/pip install Cython
/opt/python/cp27-cp27m/bin/pip install jinja2

rm -rf wheelhouse dist

/opt/python/cp27-cp27m/bin/python setup.py build bdist_wheel

/opt/python/cp37-cp37m/bin/auditwheel repair dist/*.whl

cp wheelhouse/* final_wheel
cp dist/* final_dist





# Python 3.6

/opt/python/cp36-cp36m/bin/pip install Cython
/opt/python/cp36-cp36m/bin/pip install jinja2

rm -rf wheelhouse dist

/opt/python/cp36-cp36m/bin/python setup.py build bdist_wheel

/opt/python/cp37-cp37m/bin/auditwheel repair dist/*.whl

cp wheelhouse/* final_wheel
cp dist/* final_dist





# Python 3.7

/opt/python/cp37-cp37m/bin/pip install Cython
/opt/python/cp37-cp37m/bin/pip install jinja2

rm -rf wheelhouse dist

/opt/python/cp37-cp37m/bin/python setup.py build bdist_wheel

/opt/python/cp37-cp37m/bin/auditwheel repair dist/*.whl

cp wheelhouse/* final_wheel
cp dist/* final_dist



# Python 3.8

/opt/python/cp38-cp38/bin/pip install Cython
/opt/python/cp38-cp38/bin/pip install jinja2

rm -rf wheelhouse dist

/opt/python/cp38-cp38/bin/python setup.py build bdist_wheel

/opt/python/cp37-cp37m/bin/auditwheel repair dist/*.whl

cp wheelhouse/* final_wheel
cp dist/* final_dist



rm -rf wheelhouse dist build
rm -rf spatial_index.egg-info
