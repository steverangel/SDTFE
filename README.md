# SDTFE
This code implements the method for computing a surface density field from particle data, as described in the paper: Parallel DTFE Surface Density Field Reconstruction, see author section for details. It was developed for computing gravitational lensing effects using the flat-sky approximation from N-body cosmological simulations.

## Installation

These instructions will get a copy of the project up and running on your local machine for development and testing purposes. Tested on macOSX 10.13.3 High Sierra.

Download SDTFE
```bash
git clone git@github.com:steverangel/SDTFE.git
```

### Building Using C-Make

CMake will automatically download the Qhull and libtiff dependencies.

Prerequisites:

- CMake v3.11 or later
- compiler supporting C++14 (for the python extension)

**Python library**

The python library can be installed with ``pip`` from within the SDTFE folder.
The ``setup.py`` script will internally call CMake to compile the DTFE library

```bash
cd SDTFE
pip install .
```


**Executables**

To build the executables, create a build directory and run ``cmake``

```bash
cd SDTFE
mkdir build && cd build
cmake ..
make -j4
```

### Building Using Make

Prerequisites:

- This software uses [Qhull](http://www.qhull.org) for the Delaunay triangulation library.
- Optional, but recomended, is [libtiff](http://www.libtiff.org) for visualization of the resulting field.


Download and install Qhull.

```bash
git clone git@github.com:qhull/qhull.git
cd qhull
make
```

Download and install libtiff.

```bash
curl -O http://download.osgeo.org/libtiff/tiff-4.0.9.tar.gz
tar xvfz tiff-4.0.9.tar.gz
cd tiff-4.0.9/
mkdir install
./configure --prefix=/your/local/install
make
make install
```

Edit the SDTFE Makefile to link Qhull and libtiff.

```
cd SDTFE
vim Makefile
...
QHULLLIBDIR = ../qhull/lib
QHULLINCDIR = -I../qhull/src/libqhull

TIFFLIBDIR = ../tiff-4.0.9/install/lib
TIFFINCDIR = -I../tiff-4.0.9/libtiff
...
```

Compile SDTFE using make.

```bash
make
```

--------------------------------------------------------------------------------
## Authors

* **Esteban Rangel** - *Parallel DTFE Surface Density Field Reconstruction* - [pdf](http://cucis.ece.northwestern.edu/publications/pdf/RLH16.pdf)

## License

This project is licensed under the GNU License - see the [LICENSE](https://github.com/steverangel/SDTFE/blob/master/LICENSE) file for details.

## Acknowledgments

* The Delaunay Tessellation Field Estimator - [pdf](https://arxiv.org/pdf/astro-ph/0011007.pdf)
* Fast Ray–Tetrahedron Intersection using Plücker Coordinates - [pdf](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.565.3129&rep=rep1&type=pdf)
