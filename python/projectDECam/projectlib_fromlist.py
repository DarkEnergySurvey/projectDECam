#!/usr/bin/env python
#

"""

 Set of functions to call SWarp in order to project a set or a
 single exposure into the focal plane and to draw ellipses from the
 SEx catalogs using matplotlib and PIL

 The steps involved are:
  1. Reading the files from lists
  2. SWarp the 62 chips into a single fits file
  3. Create a png [optional]
  4. Create a png with detections [optional]

 The code, will create the following files, based on the input BASENAME
 provided by the user where BASENAME can be a path+name
 such as /somedir/anotherdir/myname

   $BASENAME.png
   $BASENAME_TN.png
   $BASENAME_ell.png
   $BASENAME_ell_TN.png
   $BASENAME_proj.fits

  /somedir/anotherdir/myname.png
  /somedir/anotherdir/myname_TN.png
  /somedir/anotherdir/myname_ell.png
  /somedir/anotherdir/myname_ell_TN.png
  /somedir/anotherdir/myname_proj.fits

 Author:
  Felipe Menanteau, NCSA, Sept-Oct 2013.
  into a class library file May/June 2014
  into a new class library that read input lists

"""

import os
import sys
import time
import math

# Python external packages
import fitsio
import numpy
# ------------------------------------------
# trick to avoid X11 crash when no display
# Needs to be done before calling pylab
import matplotlib
import pylab
from PIL import Image

from despyastro import wcsutil
from despyastro import tableio
from drawDECam import drawDECam as draw  # to use ellipses
matplotlib.use('Agg')

sout = sys.stdout


class project_DECam_fromlist:

    """
    A class to project DECam exposures using SWarp and create PNGs
    using PIL/Image
    """

    def __init__(self, **kwargs):

        """ Run by path exposure """

        # Fill up the dictionaries
        for key in list(kwargs.keys()):
            self.__dict__[key] = kwargs[key]
            # print key, kwargs[key]

        # Estimate the "scaled" weight-threshold due the new pixscale
        # assumed DECam 0.263 arcsec/pix as static value
        self.INPUT_PIXEL_SCALE = 0.263

        if self.weight_thresh:
            self.weight_thresh = float(self.weight_thresh)*(self.INPUT_PIXEL_SCALE/self.pixscale)**2

        # Read in the imagalist and catlist (if present)
        self.read_filelists()

        # Now make sure that the output path exists
        self.outpath = os.path.split(self.basename)[0]
        if self.outpath == '':
            self.outpath = os.getcwd()

        print(f"# Will write files to: {self.outpath}")
        if not os.path.exists(self.outpath):
            print(f"# Making {self.outpath}")
            os.makedirs(self.outpath)

        # SWarp the exposure
        self.swarp_exposure(noSWarp=self.noSWarp,
                            noBack=self.noBack,
                            keep=self.keepfiles)
        # Create PNGs
        if not self.noPNG:
            print("# Running PNG Creation")
            self.stiff_exposure()
            self.make_png_thumbnail()
        else:
            print("# Skipping PNG Creation")

        # Draw detection
        if not self.noEll:
            print("# Running Ell Creation")
            self.read_exposure_catalogs_files()  # This one is diferent for SQL
            self.make_ell_thumbnail()
        else:
            print("# Skipping Ell on PNG")

        return

    def read_filelists(self):

        """
        Read image list (self.imglist) and optional
        catalog list (self.cataloglist)
        """

        # initialize lists
        self.scilist = []
        self.wgtlist = []

        # Read in the filelist
        print(f"# Reading image list from {self.imglist}")
        self.imgfiles = tableio.get_str(self.imglist, cols=0)

        # If we want catalogs
        if self.cataloglist:
            print(f"# Reading catalogs list from: {self.cataloglist}")
            self.catlist = tableio.get_str(self.cataloglist, cols=0)
            self.catlist.sort()

        # Sort the files
        self.imgfiles.sort()

        # Create list of images and weight
        for fname in self.imgfiles:
            scifile = "%s[0]" % (fname)
            wgtfile = "%s[2]" % (fname)
            self.scilist.append(scifile)
            self.wgtlist.append(wgtfile)
        return

    def swarp_exposure(self, noSWarp, noBack=False, keep=False):

        """ Project using Swarp, the files that make and exposure"""

        # Search swarp in the path
        swarp_exe = "swarp"
        if not inpath(swarp_exe, verb='yes'):
            sys.exit("No SWarp executable on path")

        # Get ready for SWArp, create strings
        self.scinames = ",".join(self.scilist)
        self.wgtnames = ",".join(self.wgtlist)
        self.swarp_outname = os.path.join("%s_proj.fits" % self.basename)

        if os.path.exists(self.swarp_outname) and not self.force:
            print("# SWarped file already exists")
            print("# Skipping SWArped image creation")
            return

        # Now let swarp them
        t1 = time.time()

        # Call SWarp with the tweaked options from the production pipeline
        opts = ''
        # NOTE: This only work if we use eups,might want more general solution
        try:
            print("# Loading configuration via PROJECTDECAM_DIR evironment variable")
            opts = opts + f" -c  {os.environ['PROJECTDECAM_DIR']}/etc/default.swarp"
        except KeyError:
            sys.exit("# ERROR: could not define files via PROJECTDECAM_DIR evironment variable")

        #######################################################################
        # Tries to fix the backgr calculation around large bright objects/stars
        opts = opts + ' -BACK_SIZE 128'
        opts = opts + ' -BACK_FILTERSIZE 7'
        #######################################################################
        opts = opts + ' -WEIGHT_TYPE MAP_WEIGHT'
        # Only if we want it
        if self.weight_thresh:
            opts = opts + ' -WEIGHT_THRESH %s' % self.weight_thresh
        opts = opts + ' -WEIGHT_IMAGE %s' % self.wgtnames
        opts = opts + ' -BLANK_BADPIXELS Y'
        opts = opts + ' -FSCALASTRO_TYPE VARIABLE'
        if noBack:
            opts = opts + ' -SUBTRACT_BACK N'
        else:
            opts = opts + ' -SUBTRACT_BACK Y'
        opts = opts + ' -RESAMPLING_TYPE NEAREST'  # Much faster than LANCZOS3!
        opts = opts + ' -COMBINE Y'
        opts = opts + ' -COMBINE_TYPE WEIGHTED'
        if keep:
            opts = opts + ' -DELETE_TMPFILES N'
            opts = opts + ' -RESAMPLE_DIR %s' % self.outdir
        else:
            opts = opts + ' -DELETE_TMPFILES Y'
        opts = opts + ' -WRITE_XML   N'
        opts = opts + ' -HEADER_ONLY N'
        opts = opts + ' -VERBOSE_TYPE FULL'
        opts = opts + ' -PIXELSCALE_TYPE MANUAL'
        opts = opts + ' -PIXEL_SCALE %s' % self.pixscale
        opts = opts + ' -NTHREADS %d' % self.NTHREADS_swarp
        opts = opts + ' -IMAGEOUT_NAME %s' % self.swarp_outname

        #######################################################################
        # Extra option in case we cross RA=0 and need to perform
        # manual centering
        # Check if cross RA=0 -- return self.crossRA [True/False]
        #######################################################################
        self.cross_RA_zero_center()
        if self.crossRA:
            print("# Manualy centering exposure...")
            self.get_exposure_imsize_center()
            opts = opts + ' -CENTER_TYPE MANUAL'
            opts = opts + ' -CENTER  \"%s, %s\"' % (self.RA0, self.DEC0)
            opts = opts + ' -IMAGE_SIZE %s,%s' % (self.NX, self.NY)

        swarp_cmd = "%s %s %s" % (swarp_exe, self.scinames, opts)

        if noSWarp:
            print("noSWarp invoked -- Skipping SWArp")
        else:
            os.system(swarp_cmd)
        print(f"SWarp time {elapsed_time(t1)}")

        # Clean up
        self.clean_up_weight()

        return

    def stiff_exposure(self):

        """ Stiff a DECam exposure and create a png of it """

        stiff_exe = 'stiff'
        if not inpath(stiff_exe, verb='yes'):
            sys.exit("try:\nsetup stiff 2.1.3+0\n")

        self.pngfile = f"{self.basename}.png"
        self.tiffile = f"{self.basename}.tif"

        if os.path.exists(self.pngfile) and not self.force:
            print(f"# PNG for file %s already exists {self.basename}")
            print(f"# Skipping PNG creation")
            return

        # Very Explicit Call to stiff
        opts = ''

        # NOTE: This only work if we use eups might want more general solution
        try:
            print("# Loading configuration via PROJECTDECAM_DIR evironment variable")
            opts = opts + f" -c  {os.environ['PROJECTDECAM_DIR']}/etc/default.stiff"
        except KeyError:
            sys.exit("# ERROR: could not define files via PROJECTDECAM_DIR evironment variable")

        opts = opts + " -IMAGE_TYPE   TIFF"        # Output image format
        opts = opts + " -COMPRESSION_TYPE JPEG"     # Compression type
        opts = opts + " -BINNING      1"           # Binning factor for the data
        opts = opts + " -GAMMA        2.2"         # Display gamma
        opts = opts + " -GAMMA_FAC    1.0"         # Luminance gamma correction factor
        opts = opts + " -COLOUR_SAT   1.0"         # Colour saturation (0.0 = B&W)
        opts = opts + " -NEGATIVE     N"           # Make negative of the image
        opts = opts + " -SKY_TYPE     AUTO"        # Sky-level: "AUTO" or "MANUAL"
        opts = opts + " -SKY_LEVEL    0.0"         # Background level for each image
        opts = opts + " -MIN_TYPE     GREYLEVEL"   # Min-level: "QUANTILE", "MANUAL" or "GREYLEVEL"
        opts = opts + " -MIN_LEVEL    0.005"       # Minimum value or quantile
        opts = opts + " -MAX_TYPE     QUANTILE"    # Max-level: "QUANTILE" or "MANUAL"
        opts = opts + " -MAX_LEVEL    {}".format(self.grayscale)  # Maximum value or quantile
        opts = opts + " -SATUR_LEVEL  40000.0"     # FITS data saturation level(s)
        opts = opts + " -WRITE_XML    N"           # Write XML file (Y/N)?
        opts = opts + " -COPYRIGHT    DES/NCSA"    # Copyright
        opts = opts + " -COPY_HEADER  Y"           # Copy FITS header to description field?
        opts = opts + ' -NTHREADS {}'.format(self.NTHREADS_stiff)  # Number of simultaneous threads
        opts = opts + " -OUTFILE_NAME %s" % self.tiffile
        stiff_cmd = "%s %s %s" % (stiff_exe, self.swarp_outname, opts)

        t0 = time.time()
        print(stiff_cmd)
        os.system(stiff_cmd)
        print(f"# stiff time: {elapsed_time(t0)}")

        # Create PNG using PIL/Image
        t1 = time.time()
        im = Image.open(self.tiffile)
        im.save(self.pngfile, "png", options='optimize')
        print(f"# PIL time: {elapsed_time(t1)}")

        # Clean up the tiff file
        print(f"# Cleaning tiff: {self.tiffile}")
        os.remove(self.tiffile)
        return

    def read_exposure_catalogs_files(self):

        """ Read the 62 exposure SEx catalogs"""

        # Define the output name
        self.pngfile_ell = f"{self.basename}_ell.png"

        if os.path.exists(self.pngfile_ell) and not self.force:
            print("# PNG/Ell file already exists")
            print("# Skipping PNG/Ell creation")
            return

        t0 = time.time()
        print(f"# Reading {self.pngfile}")

        image = Image.open(self.pngfile).convert("L")
        self.png_array = numpy.asarray(image)
        self.png_array = self.png_array[::-1, :]
        print(f"# Shape (ny,nx): {self.png_array.shape}")
        print(f"# Done in {time.time()-t0} sec.")

        # Figure out the size
        dpi = 90.
        (self.ny, self.nx) = self.png_array.shape
        x_size = float(self.nx)/float(dpi)
        y_size = float(self.ny)/float(dpi)

        pylab.figure(1, figsize=(x_size, y_size))
        pylab.axes([0, 0, 1, 1], frameon=False)
        pylab.imshow(self.png_array, origin='lower', cmap='gray', interpolation='none')
        ec_ima1 = 'red'
        ec_ima2 = 'blue'
        pylab.axis('off')

        i = 0
        t0 = time.time()

        for catfile in self.catlist:
            print(f"# Reading {catfile}")
            tbdata = fitsio.read(catfile, ext=2)
            # Store the Relevant information to draw ellipse
            if i == 0:
                ra = tbdata['ALPHA_J2000']
                dec = tbdata['DELTA_J2000']
                a_image = tbdata['A_IMAGE']*tbdata['KRON_RADIUS']
                b_image = tbdata['B_IMAGE']*tbdata['KRON_RADIUS']
                theta = tbdata['THETA_IMAGE']
                imaflag_iso = tbdata['IMAFLAGS_ISO']

            else:
                ra = numpy.append(ra, tbdata['ALPHA_J2000'])
                dec = numpy.append(dec, tbdata['DELTA_J2000'])
                a_image = numpy.append(a_image, tbdata['A_IMAGE']*tbdata['KRON_RADIUS'])
                b_image = numpy.append(b_image, tbdata['B_IMAGE']*tbdata['KRON_RADIUS'])
                theta = numpy.append(theta, tbdata['THETA_IMAGE'])
                imaflag_iso = numpy.append(imaflag_iso, tbdata['IMAFLAGS_ISO'])

            i = i+1

        print(f"# Read {i} SEx catalogs in time: {elapsed_time(t0)}")
        # Now let's put the positions on the projected image
        hdr = fitsio.read_header(self.swarp_outname)
        wcs = wcsutil.WCS(hdr)
        x, y = wcs.sky2image(ra, dec)

        # Draw all at once -- faster
        t1 = time.time()
        print(f"# Drawing ellipses for {len(x)} objects")

        # Drawing imaflags > 0 blue and the rest 'red'
        idx1 = numpy.where(imaflag_iso == 0)
        idx2 = numpy.where(imaflag_iso > 0)
        draw.PEllipse_multi((x[idx1], y[idx1]),
                            (a_image[idx1], b_image[idx1]),
                            resolution=60,
                            angle=theta,
                            facecolor='none',
                            edgecolor=ec_ima1,
                            linewidth=0.5)
        draw.PEllipse_multi((x[idx2], y[idx2]),
                            (a_image[idx2], b_image[idx2]),
                            resolution=60,
                            angle=theta,
                            facecolor='none',
                            edgecolor=ec_ima2,
                            linewidth=0.5)
        print(f"# Ellipses draw time: {elapsed_time(t1)}")
        print("# Saving PNG file with ellipses")
        pylab.savefig(self.pngfile_ell, dpi=dpi)
        print("# Done")
        pylab.close()
        return

    def make_png_thumbnail(self):

        """ Make the thumbnail of the PNG for the webpages"""

        self.TN_png = "%s_TN.png" % self.basename
        if os.path.exists(self.TN_png) and not self.force:
            print(f"# File exists -- Skipping creation of {self.TN_png}")
        else:
            print(f"# Creating TN {self.TN_png}")
            im = Image.open(self.pngfile)
            (w, h) = im.size  # Get width (w) and height (h)
            TNshape = self.TNsize, int(self.TNsize*h/w)
            im.thumbnail(TNshape, Image.ANTIALIAS)
            im.save(self.TN_png, "png", options='optimize')
        return

    def make_ell_thumbnail(self):

        """ Make TN for the PNG with ellipses"""

        self.TN_ell = '%s_ell_TN.png' % self.basename
        if os.path.exists(self.TN_ell) and not self.force:
            print(f"# File exists -- Skipping creation of {self.TN_ell}")
        else:
            print(f"# Creating TN {self.TN_ell}")
            im = Image.open(self.pngfile_ell)
            (w, h) = im.size  # Get width (w) and height (h)
            TNshape = self.TNsize, int(self.TNsize*h/w)
            im.thumbnail(TNshape, Image.ANTIALIAS)
            im.save(self.TN_ell, "png", options='optimize')
        return

    def cross_RA_zero_center(self):

        """
        Check if crossing the RA=0.0 by the min and max
        values for the centers of all CCDs. We loop over all of the
        wcs of the images... slower (a few secs) but safer in case
        some CCDs are missing.
        """

        d2r = math.pi/180.  # degrees to radians shorthand
        DECam_width = 2.3   # approx width in degrees

        # Collect the centers of all CCDS first and store to compare later
        # Slowers than asking for particular CCDs, but safer
        print("# Figuring out edges of CCDs")
        t0 = time.time()
        ra0 = []
        dec0 = []
        for filename in self.imgfiles:
            hdr = fitsio.read_header(filename)
            wcs = wcsutil.WCS(hdr)
            nx = hdr['NAXIS1']
            ny = hdr['NAXIS2']
            ra, dec = wcs.image2sky(nx/2.0, ny/2.0)
            ra0.append(ra)
            dec0.append(dec)

        self.ra0 = numpy.array(ra0)
        self.dec0 = numpy.array(dec0)

        # Get min and max values
        dec1 = self.dec0.min()
        dec2 = self.dec0.max()
        ra1 = self.ra0.min()
        ra2 = self.ra0.max()

        print(f"# Cross RA examination: {elapsed_time(t0)}")
        dec0 = (dec1+dec2)/2.0
        D_ra = abs(ra1 - ra2)
        D_dec = abs(dec1 - dec2)

        # in case ras near 360 have negative values
        # so, we need to check that both ras have the same sign
        if math.copysign(1, ra1) == math.copysign(1, ra2):
            same_sign = True
        else:
            same_sign = False

        # Define the tolerances for ra/dec
        tol_ra = 2.5*DECam_width/math.cos(dec0*d2r)
        tol_dec = 2.5*DECam_width

        # Check WCS solution is sane, by making sure that the distance
        # between ra.min and ra.max or dec.min/dec.max in D > tol and
        # D < 360-tol (for RA)
        if (D_ra > tol_ra and D_ra < 360 - tol_ra) or D_dec > tol_dec:
            print("# ***************************************************************")
            print("# **  WARNING: Distance between CCDs greated that DECam FOV    **")
            print("# **  WARNING: Projection will not be performed -- bye         **")
            print("# ***************************************************************")
            sys.exit()

        # The tolerace for deciding when we crossed RA=0 is ~2x a
        # DECam FOV at the declination
        if D_ra > tol_ra or same_sign is False:
            self.crossRA = True
            print("# *****************************************")
            print("# **  WARNING: exposure crosses RA=0.0   **")
            print("# *****************************************")

            # Bring to zero the values near 360
            idx = numpy.where(self.ra0 > tol_ra)
            self.ra0[idx] = self.ra0[idx] - 360
        else:
            self.crossRA = False

    def get_exposure_imsize_center(self):

        """
        Get the center and size of the output projected
        exposure image if required
        """

        # Get the right size, the default size of a projected image in the
        # native pixscale of DECam is 31123x28149 pixels
        NX = 31123  #
        NY = 28149  # at 0.263 arcsec/pixel
        self.NX = int(NX*self.INPUT_PIXEL_SCALE/self.pixscale)
        self.NY = int(NY*self.INPUT_PIXEL_SCALE/self.pixscale)

        self.DEC0 = (self.dec0.min() + self.dec0.max())/2.0
        self.RA0 = (self.ra0.min() + self.ra0.max())/2.0

        # in case we fall bellow zero
        if self.RA0 < 0:
            self.RA0 = self.RA0+360

        print("# RE-CENTER is:")
        print(f"# RA : {self.RA0}")
        print(f"# DEC: {self.DEC0}")
        print(f"# ra_min,  ra_max: {self.ra0.min()},{self.ra0.max()} [degrees]")
        print(f"# dec_min, dec_max: {self.dec0.min()},{self.dec0.max()} [degrees]")

    def clean_up_weight(self):

        """ Clean up the weight.fits image created by swarp"""
        wgt = "coadd.weight.fits"
        if os.path.exists(wgt):
            print(f"# Cleaning up: {wgt}")
            os.remove(wgt)
        return


def elapsed_time(t1, verb=False):
    """
    Returns the time between t1 and the current time now
    I can can also print the formatted elapsed time.
    ----------
    t1: float
        The initial time (in seconds)
    verb: bool, optional
        Optionally print the formatted elapsed time
    returns
    -------
    stime: float
        The elapsed time in seconds since t1
    """
    t2 = time.time()
    stime = "%dm %2.2fs" % (int((t2-t1)/60.), (t2-t1) - 60*int((t2-t1)/60.))
    if verb:
        print(f"Elapsed time: {stime}")
    return stime


# Check if executable is in path of user
def inpath(program, verb=None):
    """ Checks if program is in the user's path """
    import os.path
    for path in os.environ['PATH'].split(':'):
        if os.path.exists(os.path.join(path, program)):
            if verb:
                print(f"# program: {program} found in: {os.path.join(path, program)}")
            return 1
    if verb:
        print(f"# program: {program} NOT found in user's path")
    return 0


def cmdline():

    """ Parse the command line arguments and options using argparse"""

    import argparse

    USAGE = "\n"
    USAGE = USAGE + "  %(prog)s <imglist> <basename> [options] \n"
    USAGE = USAGE + "  i.e.: \n"
    USAGE = USAGE + "  %(prog)s file-with-list-of-images /someplace/anotherplace/somename\n"

    epilog = "Author: Felipe Menanteau, NCSA/University of Illinois (felipe@illinois.edu)"
    description = "Projects a DECam exposure using SWarp and creates grayscale PNGs" \
                  "using Python's native PIL/Image module"
    parser = argparse.ArgumentParser(usage=USAGE,
                                     epilog=epilog,
                                     description=description,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The positional arguments
    parser.add_argument("imglist", default=None,
                        help="Image list to project")
    parser.add_argument("basename", action="store",
                        help="Output Directory w/BASENAME")
    parser.add_argument("--catlist", dest='cataloglist',
                        help="List of catalogs")
    parser.add_argument("--noPNG", action="store_true", default=False,
                        help="Skip creation of PNG files")
    parser.add_argument("--noEll", action="store_true", default=False,
                        help="Skip creation of PNG file with overlaid ellipse")
    parser.add_argument("--TNsize", type=int, default=800,
                        help="Size of thumbnail [pixels]")
    parser.add_argument("--grayscale", type=float, default=0.98,
                        help="grayscale for png creation [i.e. 0.98]")
    parser.add_argument("--pixscale", type=float, action="store", default=1.0,
                        help="pixel-scale in arcsec/pix")
    parser.add_argument("--noSWarp", action="store_true", default=False,
                        help="No SWarp -- dry run")
    parser.add_argument("--weight_thresh", action="store", default=None,
                        help="SWarps WEIGHT_THRESH value")
    parser.add_argument("--force", action="store_true", default=False,
                        help="Forces the re-creation of existing files")
    parser.add_argument("--dryrun", action="store_true", default=False,
                        help="Dry Run -- only build the lists")
    parser.add_argument("--keepfiles", action="store_true", default=False,
                        help="Keep each CCD projected file")
    parser.add_argument("--noBack", action="store_true", default=False,
                        help="Avoids Background substraction on SWarp call")
    parser.add_argument("--NTHREADS_swarp", type=int, default=0,
                        help="NTHREADS for SWarp [0=auto]")
    parser.add_argument("--NTHREADS_stiff", type=int, default=0,
                        help="NTHREADS for stiff [0=auto]")

    args = parser.parse_args()

    # noPNG turns off noEll
    if args.noPNG:
        args.noEll = True

    # Dry run turns off everything... almost
    if args.dryrun:
        args.noPNG = True
        args.noEll = True
        args.noSWarp = True

    print("# Will run:")
    print(f"# {parser.prog}")
    for key, val in sorted(vars(args).items()):
        print("# \t--%-10s\t%s" % (key, val))
    return args
