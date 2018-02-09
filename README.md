# SDTFE
This code computes the surface density field from a particle set using the DTFE method described in the paper, Parallel DTFE Surface Density Field Reconstruction.  

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

This software uses Qhull for the Delaunay triangulation library. http://www.qhull.org

Optionally, but recomended, is libtiff for visualization of the resulting field. http://www.libtiff.org

### Installing

Download and install Qhull.

```
git clone git@github.com:qhull/qhull.git
cd qhull
make
```

Download and install libtiff.

```
curl -O http://download.osgeo.org/libtiff/tiff-4.0.9.tar.gz
./configure
make
```

Download SDTFE.

```
git clone git@github.com:steverangel/SDTFE.git
```

Edit the Makefile for SDTFE for linking to Qhull and libtiff.

```
QHULLLIBDIR = ../qhull/lib 
QHULLINCDIR = -I../qhull/src/libqhull
...
TIFFLIBDIR = ../tiff-4.0.9/install/lib
TIFFINCDIR = -I../tiff-4.0.9/libtiff
```

Compile SDTFE using make.

```
make
```

## Running the examples

Explain how to run the automated tests for this system


## Authors

* **Esteban Rangel** - *Parallel DTFE Surface Density Field Reconstruction* - [paper](http://ieeexplore.ieee.org/document/7776476/)

## License

This project is licensed under the GNU License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* The Delaunay Tessellation Field Estimator - [paper](https://arxiv.org/pdf/astro-ph/0011007.pdf)
* Fast Ray–Tetrahedron Intersection using Plücker Coordinates - [paper](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.565.3129&rep=rep1&type=pdf)
