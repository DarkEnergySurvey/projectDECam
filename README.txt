Creates full focal plane fits files and grayscale PNGs
of DECam exposure images. This version uses a filelist on input CCDs
and catalogs. These codes have been migrated from qatoolkit

Felipe Menanteau, Apr 2014. 

Requirements:
-------------

- despyastro (loads numpy, scipy, matplotlib and scipy)
- PIL
- fitsio
- stiff
- swarp
- pyfits
- drawDECam

To set up for testing:
----------------------
 setup -r ~/DESDM-Code/devel/projectDECam/trunk
   or
 setup -r ~/DESDM-Code/devel/projectDECam/tags/XX.YY.ZZ

To install:
-----------
 
 python setup install --home=$HOME/Python
  or
 python setup.py install --prefix=$PRODUCT_DIR --install-lib=$PRODUCT_DIR/python 


Scripts: (inside bin/)
----------------------

 - projectDECamPNG : Projects a DECam exposure using SWarp and creates
 grayscale PNGs using Python's native PIL/Image module. It reads the
 input files from lists provided command-line. This script is aimed for the
 Refacted system.

