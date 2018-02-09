# SDTFE
This code computes the surface density field from a particle set using the DTFE method described in the paper, Parallel DTFE Surface Density Field Reconstruction.  

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

This software uses Qhull for the Delaunay triangulation library. http://www.qhull.org

Optional, but recomended, is libtiff for visualization of the resulting field. http://www.libtiff.org

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
tar xvfz tiff-4.0.9.tar.gz
cd tiff-4.0.9/
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

Download and extract example data

```
curl -O http://users.eecs.northwestern.edu/~emr126/sdtfe_examples.tar.gz
tar xvfz sdtfe_examples.tar.gz
```

Run the examples with all the parameters

```
./bin/dtfe ../data/913571938961.bin 216683 768 2556.9 1510.4 1986.6 6.0 4.0 1 0.01 5 0.25
```

or by prompts

```
./bin/dtfe
```

## Authors

* **Esteban Rangel** - *Parallel DTFE Surface Density Field Reconstruction* - [paper](http://ieeexplore.ieee.org/document/7776476/)

## License

This project is licensed under the GNU License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* The Delaunay Tessellation Field Estimator - [paper](https://arxiv.org/pdf/astro-ph/0011007.pdf)
* Fast Ray–Tetrahedron Intersection using Plücker Coordinates - [paper](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.565.3129&rep=rep1&type=pdf)
