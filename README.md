# SDTFE
This code implements the method for computing a surface density field from particle data, as described in the paper: Parallel DTFE Surface Density Field Reconstruction, see author section for details. It was developed for computing gravitational lensing effects using the flat-sky approximation from N-body cosmological simulations.  

## Getting Started

These instructions will get a copy of the project up and running on your local machine for development and testing purposes. Tested on macOSX 10.13.3 High Sierra. 

### Prerequisites

This software uses [Qhull](http://www.qhull.org) for the Delaunay triangulation library. 

Optional, but recomended, is [libtiff](http://www.libtiff.org) for visualization of the resulting field. 

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
mkdir install
./configure --prefix=/your/local/install
make
make install
```

Download SDTFE.

```
git clone git@github.com:steverangel/SDTFE.git
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

```
make
```

## Running the examples

Download and extract example data.

```
curl -O http://users.eecs.northwestern.edu/~emr126/sdtfe_examples.tar.gz
tar xvfz sdtfe_examples.tar.gz
```

Run the examples with all the parameters.

Usage: dtfe [ path\_to\_file n\_particles grid\_dim center\_x center\_y center\_z field\_width field\_depth particle\_mass mc\_box\_width n\_mc\_samples sample\_factor ]

```
./bin/dtfe ../data/913571938961.bin 216683 768 2556.9 1510.4 1986.6 6.0 4.0 1 0.01 5 0.25
```

Or by terminal prompts.

```
./bin/dtfe
```

## Authors

* **Esteban Rangel** - *Parallel DTFE Surface Density Field Reconstruction* - [pdf](http://cucis.ece.northwestern.edu/publications/pdf/RLH16.pdf)

## License

This project is licensed under the GNU License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

* The Delaunay Tessellation Field Estimator - [pdf](https://arxiv.org/pdf/astro-ph/0011007.pdf)
* Fast Ray–Tetrahedron Intersection using Plücker Coordinates - [pdf](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.565.3129&rep=rep1&type=pdf)
