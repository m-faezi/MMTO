MDMTO
==============

Multi-dimensional faint object detection


Installation
------------

**Requires a C++ 14 compiler and cmake**

 - `python -m venv mmto`
 - `source mmto/bin/activate`
 - `./recompile.sh (make sure to have gsl)`
 - `pip install ./MDMTO`

Build a binary wheel
--------------------
 
A binary wheel ease the redistribution of your project and can be installed with *pip* on a client machine without a compiler.

**Create wheel**

 - `cd MDMTO`
 - `python setup.py bdist_wheel`
 - `pip install ./MDMTO`
 
 The wheel is created in the directory `MDMTO/dist`, it will be named `mdmto-XXXXX.whl` where `XXXXXX` are name tags identifying the current platform and Python version. 
 
**Install wheel**
 
A wheel can be installed with *pip*:
 
 - `pip install wheel_name.whl`
 
 Note that a binary wheel is specific to a platform and to a python version (a wheel built on Windows with Python 3.5 can only be installed on Windows with Python 3.5).

Tests
-----

Tests are run automatically at the end of a build: the build will fail if tests are not successful. 

Known Issues
------------

Clang on Linux may not work due to ABI compatibilty issues.
