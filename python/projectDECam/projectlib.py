#!/usr/bin/env python
#

"""

 Set of function class to call SWarp in order to project a set or a
 single exposure into the focal plane.

 The steps involved are:
  0. Figure out which files will be used.
  1. funpack the files into .fits MEF files
  2. SWarp the 62 chips into a single fits file
  3. Create a png [optional]
  4. Create a png with detections [optional]

 Author:
  Felipe Menanteau, NCSA, Sept-Oct 2013.
  into a class library file May/June 2014
  
"""

import os,sys
import time
import glob
import shlex
import multiprocessing
import subprocess
import math

# Python external packages
import pyfits
import numpy
# ------------------------------------------
# trick to avoid X11 crash when no display
# Needs to be done before calling pylab
import matplotlib
matplotlib.use('Agg')
import pylab
from PIL import Image
import despydb.desdbi
from despyastro import wcsutil
from drawDECam import drawDECam as draw # to use ellipses

sout = sys.stdout

class project_DECam:

    """
    A class to project DECam exposures using SWarp and create PNGs using PIL/Image
    proj = project_DECam(runID,projectID,**kwargs)
    proj.run(SQL=True,**kwargs)
     or
    proj.run(SQL=False,**kwargs)
    """

    def __init__(self,runID,exposureID):
        self.runID      = runID
        self.exposureID = exposureID

    def run(self,SQL=False,**kwargs):

        """ Run by path exposure """

        t0 = time.time()
        # In case we can to use a query to figure out the location of the files
        self.SQL = SQL
        # Fill up the dictionaries
        for key in kwargs.keys():
            self.__dict__[key] = kwargs[key]
            
        # Estimate the "scaled" weight-threshold due the new pixscale
        # assumed DECam 0.263 arcsec/pix as static value
        self.INPUT_PIXEL_SCALE = 0.263
        self.weight_thresh = self.weight_thresh*(self.INPUT_PIXEL_SCALE/self.pixscale)**2
        
        # Get the information on the project for the run via SQL
        if self.SQL:
            self.connectDB()
            self.query_desar_project()

        # Build the list of files using the path and glob
        self.build_filelist()
        
        # Make sure that the output directory exists
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
            
        # Funapck the files
        self.funpack_directory() # This sets the names

        # SWarp the exposure
        self.swarp_exposure(noSWarp=self.noSWarp,
                            noBack=self.noBack,
                            keep=self.keepfiles)
        # Create PNGs
        if not self.noPNG:
            self.stiff_exposure()
            self.make_png_thumbnail()
        else:
            print "# Skipping PNG Creation"
    
        # Draw detection
        if not self.noEll:
            if self.SQL:
                self.read_exposure_catalogs_SQL() # If SQL
            else:
                self.read_exposure_catalogs_files() # This one is diferent for SQL
                
            self.make_ell_thumbnail()
        else:
            print "# Skipping Ell on PNG"
            print "# Exposure %s Total time: %s" % (self.exposure_name,elapsed_time(t0))

        return

    def connectDB(self):

        """ Setup desar connection for queries """

        section = "db-desoper"
        try:
            desdmfile = os.environ["des_services"]
        except KeyError:
            desdmfile = None
        self.dbh = despydb.desdbi.DesDbi(desdmfile,section)

    def query_desar_project(self):

        """
        Find out the project for the runID, in order to build the path
        and search for the files.
        """
        
        cur    = self.dbh.cursor()

        #####################################################################
        # PART 1 -- Get the archive root -- usually /archive_data/Archive/
        #####################################################################
        query = "select archive_root from archive_sites where location_name='desardata'"
        print "# Getting the archive root name"
        print "# Will execute the SQL query:\n********\n** %s\n********" % query
        cur.execute(query)
        archive_root = cur.fetchone()[0]

        #####################################################################
        # PART 2 --  Get the location of the filenames
        #####################################################################
        query = """select path from filepath_desar, location
           where
            filepath_desar.id = location.id and
            location.run = '%s' and
            location.filetype='red' """ % self.runID

        print "# Will execute the SQL query:\n********\n** %s\n********" % query
        cur.execute(query)
        item = cur.fetchone()[0] # Relative path
        rel_path = os.path.split(os.path.dirname(item))[0]
        
        self.runID_path = os.path.join(archive_root,rel_path)
        print "# Will search for files in:"
        print "#  %s" % self.runID_path
        return

    def build_filelist(self):

        """ Build file list """

        print "# Building the list of files to be used"
        self.exposure_path = os.path.join(self.runID,self.exposureID)
        self.exposure_name = self.exposureID
        print "# Will project single exposure %s" % self.exposureID
        return

    def funpack_directory(self): 

        """ funpack/copy or link a DECam exposure directory at a time"""

        # initialize lists
        self.catlist   = []
        self.outlist   = []
        self.scilist   = []
        self.wgtlist   = []
        self.fpack_cmd = []

        # 1. Search for fz fits files,
        # 2. otherwise search for normal fist files
        pattern = "%s/%s_[0-9][0-9].fits.fz" % (self.exposure_path,self.exposure_name)
        #pattern = "%s/%s_[0][0-9].fits.fz" % (self.exposure_path,self.exposure_name) # for just a few to run faster while debugging.
        files = glob.glob(pattern)

        # Use this opportunity to get the cat.fits files
        cat_pattern = "%s/%s_[0-9][0-9]_cat.fits" % (self.exposure_path,self.exposure_name)
        self.catlist =  glob.glob(cat_pattern)
        self.catlist.sort()
        
        if len(files) > 0:
            print "# Found .fz files -- will funzip now"
            do_funpack = True
        else:
            pattern = "%s/%s_[0-9][0-9].fits" % (self.exposure_path,self.exposure_name)
            files = glob.glob(pattern)
            do_funpack = False

        # Sort the files
        files.sort()

        # Make sure we can execute funpack
        if not inpath("funpack",verb='yes'):
            sys.exit("No funpack executable found\n\tTry setup cfitsio:")

        # Now make sure that the output path exists
        self.outpath = os.path.join(self.outdir,self.exposure_name)
        print "# Will store file in: %s" % self.outpath
        if not os.path.exists(self.outpath):
            print "# Making %s" % self.outpath
            os.mkdir(self.outpath)

        # Go trough all 62 images
        t0 = time.time()
        print "# funziping files in: %s" % self.exposure_path
        for fname in files:

            # Check if dealing with .fz or .fits files
            if os.path.basename(os.path.splitext(fname)[-1]) == '.fz':
                base = os.path.basename(os.path.splitext(fname)[0])
            elif os.path.basename(os.path.splitext(fname)[-1]) == '.fits':
                base = os.path.basename(fname)
            else:
                print "ERROR: No .fz or .fits files found"
                
            outfile = "%s/%s"    % (self.outpath,base)
            scifile = "%s/%s[0]" % (self.outpath,base)
            wgtfile = "%s/%s[2]" % (self.outpath,base)
         
	    if do_funpack:
		# funpack Full Image
   	        cmd = "funpack  -O %s %s" % (outfile,fname)
	    else:
                if self.force: # -n for noclobber
                    cmd = "cp -pv  %s %s" % (fname,outfile)
                elif self.noCopy:
                    cmd = 'echo %s' % fname
                    outfile = fname
                    scifile = "%s[0]" % (fname)
                    wgtfile = "%s[2]" % (fname)
                else:
                    cmd = "cp -pvn %s %s" % (fname,outfile)

            # MP-mode
            if self.MP:
                self.fpack_cmd.append(shlex.split(cmd))
            # Non-MP mode
            else:
                sout.write("\rfunpacking %s " % fname);sout.flush()
                os.system(cmd)
            
            # Create list of images and weight
            self.scilist.append(scifile)
            self.wgtlist.append(wgtfile)
            self.outlist.append(outfile)
            
        # Clean terminal
        sout.write("\n")

        if self.MP:
            count = multiprocessing.cpu_count()
            print "# Will Use %s threads" % count
            pool = multiprocessing.Pool(processes=count)
            pool.map(work, self.fpack_cmd)

        # Time in fpack
	if do_funpack:
           print "fpack time: %s" % elapsed_time(t0)
	else: 	
	   print "copy time: %s" % elapsed_time(t0)

        return

    def swarp_exposure(self,noSWarp,noBack=False,keep=False):
        
        # Search swarp in the path
        swarp_exe = "swarp"
        if not inpath(swarp_exe,verb='yes'):
            sys.exit("No SWarp executable on path")

        # Get ready for SWArp, create strings
        self.scinames =  ",".join(self.scilist)
        self.wgtnames =  ",".join(self.wgtlist)
        self.base_outname  = os.path.join(self.outpath,"%s_proj"      % self.exposure_name)
        self.swarp_outname = os.path.join(self.outpath,"%s_proj.fits" % self.exposure_name)

        if os.path.exists(self.swarp_outname) and not self.force:
            print "# SWarped file already exists"
            print "# Skipping SWArped image creation"
            return
        
        # Now let swarp them
        t1 = time.time()

        # Call SWarp with the tweaked options from the production pipeline
        opts = ''
        # CHANGE: This only work if we use eups, might want more general solution
        try:
            print "# Loading configuration via QATOOLKIT_DIR evironment variable"
            opts = opts + ' -c  %s/etc/default.swarp' % os.environ['QATOOLKIT_DIR']
        except:
            sys.exit("# ERROR: could not define files via QATOOLKIT_DIR evironment variable")
        ##############################################################################
        # Tries to fix the background calculation around large bright objects/stars
        opts = opts + ' -BACK_SIZE 128'
        opts = opts + ' -BACK_FILTERSIZE 7'
        ##############################################################################
        opts = opts + ' -WEIGHT_TYPE MAP_WEIGHT'
        opts = opts + ' -WEIGHT_THRESH %s' % self.weight_thresh
        opts = opts + ' -WEIGHT_IMAGE %s' % self.wgtnames
        opts = opts + ' -BLANK_BADPIXELS Y'
        opts = opts + ' -FSCALASTRO_TYPE VARIABLE'
        if noBack:
            opts = opts + ' -SUBTRACT_BACK N'
        else:
            opts = opts + ' -SUBTRACT_BACK Y'
        opts = opts + ' -RESAMPLING_TYPE NEAREST' # Much faster than LANCZOS3!
        opts = opts + ' -COMBINE Y'
        opts = opts + ' -COMBINE_TYPE WEIGHTED'
        if keep:
            opts = opts + ' -DELETE_TMPFILES N'
            opts = opts + ' -RESAMPLE_DIR %s' % self.outdir
        else:
            opts = opts + ' -DELETE_TMPFILES Y'
        opts = opts + ' -HEADER_ONLY N'
        opts = opts + ' -VERBOSE_TYPE FULL'
        opts = opts + ' -PIXELSCALE_TYPE MANUAL'
        opts = opts + ' -PIXEL_SCALE %s' % self.pixscale
        opts = opts + ' -NTHREADS 0'
        opts = opts + ' -IMAGEOUT_NAME %s' % self.swarp_outname

        ############################################################################
        # Extra option in case we cross RA=0 and need to perform manual centering
        # Check if cross RA=0 -- return self.crossRA [True/False]
        ############################################################################
        self.cross_RA_zero()
        if self.crossRA:
            print "# Manualy centering exposure..."
            self.get_exposure_imsize_center()
            opts = opts + ' -CENTER_TYPE MANUAL'
            opts = opts + ' -CENTER  \"%s, %s\"' % (self.RA0,self.DEC0)
            opts = opts + ' -IMAGE_SIZE %s,%s' % (self.NX,self.NY)

        swarp_cmd = "%s %s %s" % (swarp_exe,self.scinames,opts)

        if noSWarp:
            print "noSWarp invoked -- Skipping SWArp"
        else:
            os.system(swarp_cmd)
        print "SWarp time %s" % elapsed_time(t1)
        return


    def stiff_exposure(self):

        """ Stiff a DECam exposure and create a png of it """

        stiff_exe = 'stiff'
        if not inpath(stiff_exe,verb='yes'):
            sys.exit("try:\nsetup stiff 2.1.3+0\n")

        self.pngfile = "%s.png" % self.base_outname
        self.tiffile = "%s.tif" % self.base_outname
        
        if os.path.exists(self.pngfile) and not self.force:
            print "# PNG for exposure %s file already exists" % self.exposureID
            print "# Skipping PNG creation"
            return

        # Very Explicit Call to stiff
        opts = ''
        opts = opts + " -IMAGE_TYPE   TIFF"        # Output image format
        opts = opts + " -COMPRESSION_TYPE LZW"     # Compression type
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
        opts = opts + " -MAX_LEVEL    %s" % self.grayscale # Maximum value or quantile 
        opts = opts + " -SATUR_LEVEL  40000.0"     # FITS data saturation level(s)
        opts = opts + " -WRITE_XML    N"           # Write XML file (Y/N)?
        opts = opts + " -COPYRIGHT    DES/NCSA"    # Copyright
        opts = opts + " -COPY_HEADER  Y"           # Copy FITS header to description field?
        opts = opts + " -OUTFILE_NAME %s" % self.tiffile
        stiff_cmd = "%s %s %s" % (stiff_exe,self.swarp_outname,opts)

        t0 = time.time()
        print stiff_cmd
        os.system(stiff_cmd)
        print "# stiff time: %s" % elapsed_time(t0)
        # Create PNG using PIL/Image
        t1 = time.time()
        im = Image.open(self.tiffile)
        im.save(self.pngfile, "png",options='optimize')
        print "# PIL time: %s" % elapsed_time(t1)

        # Clean up the tiff file
        print "# Cleaning tiff: %s" % self.tiffile
        os.remove(self.tiffile)
        return

    def read_exposure_catalogs_files(self):

        """ Read the 62 exposure catalogs"""

        # Define the output name
        self.pngfile_ell = '%s_ell.png' % (self.base_outname)
        
        if os.path.exists(self.pngfile_ell) and not self.force:
            print "# PNG/Ell file already exists"
            print "# Skipping PNG/Ell creation"
            return

        t0 = time.time()
        print >>sys.stderr,"# Reading %s" % self.pngfile

        image = Image.open(self.pngfile).convert("L")
        self.png_array = numpy.asarray(image)
        self.png_array = self.png_array[::-1,:]
        print "# Shape (ny,nx): ", self.png_array.shape
        print "# Done in %.3f sec." % (time.time()-t0)

        # Figure out the size
        dpi = 90.
        (self.ny,self.nx) = self.png_array.shape
        x_size = float(self.nx)/float(dpi)
        y_size = float(self.ny)/float(dpi)

        pylab.figure(1,figsize=(x_size,y_size))
        pylab.axes([0,0,1,1], frameon=False)
        pylab.imshow(self.png_array,origin='lower',cmap='gray',interpolation='none')
        ax = pylab.gca()
        ec = 'red'
        pylab.axis('off')

        i = 0
        t0 = time.time()
        for catfile in self.catlist:

            print "# Reading %s" % catfile
            hdulist = pyfits.open(catfile)
            tbdata = hdulist[2].data
            # Store the Relevant information to draw ellipse
            if i==0:
                ra      = tbdata['ALPHA_J2000']
                dec     = tbdata['DELTA_J2000']
                a_image = tbdata['A_IMAGE']*tbdata['KRON_RADIUS']
                b_image = tbdata['B_IMAGE']*tbdata['KRON_RADIUS']
                theta   = tbdata['THETA_IMAGE']
            else:
                ra      = numpy.append(ra,tbdata['ALPHA_J2000'])
                dec     = numpy.append(dec,tbdata['DELTA_J2000'])
                a_image = numpy.append(a_image,tbdata['A_IMAGE']*tbdata['KRON_RADIUS'])
                b_image = numpy.append(b_image,tbdata['B_IMAGE']*tbdata['KRON_RADIUS'])
                theta   = numpy.append(theta,tbdata['THETA_IMAGE'])

            hdulist.close()
            i = i+1

        print "# Read %s SEx catalogs in time: %s" % (i,elapsed_time(t0))
        # Now let's put the positions on the projected image
        hdr = pyfits.getheader(self.swarp_outname)
        wcs = wcsutil.WCS(hdr)
        x,y = wcs.sky2image(ra,dec)

        ###########################################################################################
        # Too slow
        # Draw each at a time -- required for ellipse call
        #for k in range(len(x)):
        #    E = PEllipse((x[k],y[k]),(a_image[k],b_image[k]),resolution=60,angle=theta[k],fill=0, edgecolor=ec,linewidth=0.5)
        #    ax.add_patch(E)
        ###########################################################################################
        
        # Draw all at once -- faster
        t1 = time.time()
        print "# Drawing ellipses for %s objects" % len(x)
        draw.PEllipse_multi((x,y),(a_image,b_image),resolution=60,angle=theta,facecolor='none',edgecolor=ec,linewidth=0.5)
        print "# Ellipses draw time: %s" % elapsed_time(t1)

        print "# Saving png with ellipses"
        pylab.savefig(self.pngfile_ell,dpi=dpi)
        print "# Done"
        pylab.close()
        return


    def read_exposure_catalogs_SQL(self):

        """ Read the 62 exposure catalogs via SQL"""

        # Define the output name
        self.pngfile_ell = '%s_ell.png' % (self.base_outname)
        
        if os.path.exists(self.pngfile_ell) and not self.force:
            print "# PNG/Ell file already exists"
            print "# Skipping PNG/Ell creation"
            return

        t0 = time.time()
        print >>sys.stderr,"# Reading %s" % self.pngfile

        image = Image.open(self.pngfile).convert("L")
        self.png_array = numpy.asarray(image)
        self.png_array = self.png_array[::-1,:]
        print "# Shape (ny,nx): ", self.png_array.shape
        print "# Done in %.3f sec." % (time.time()-t0)

        # Figure out the size
        dpi = 90.
        (self.ny,self.nx) = self.png_array.shape
        x_size = float(self.nx)/float(dpi)
        y_size = float(self.ny)/float(dpi)

        pylab.figure(1,figsize=(x_size,y_size))
        pylab.axes([0,0,1,1], frameon=False)
        pylab.imshow(self.png_array,origin='lower',cmap='gray',interpolation='none')
        ax = pylab.gca()
        ec = 'red'
        pylab.axis('off')

        # Here we read in the relevant information to draw the ellipse
        self.query_cat_ellipses()

        # Now let's put it on the projected image
        hdr = pyfits.getheader(self.swarp_outname)
        wcs = wcsutil.WCS(hdr)
        x,y = wcs.sky2image(self.ra,self.dec)

        ###########################################################################################
        # Too slow - deprecated
        # Draw each at a time -- required for ellipse call
        #for k in range(len(x)):
        #    E = draw.PEllipse((x[k],y[k]),
        #                 (self.a_image[k],self.b_image[k]),
        #                 resolution=60,angle=self.theta[k],fill=0, edgecolor=ec,linewidth=0.5)
        #    ax.add_patch(E)
        ###########################################################################################
        
        # Draw all at once -- faster
        t1 = time.time()
        print "# Drawing ellipses for %s objects" % len(x)
        draw.PEllipse_multi((x,y),(self.a_image,self.b_image),resolution=60,angle=self.theta,facecolor='none',edgecolor=ec,linewidth=0.5)
        print "# Ellipses draw time: %s" % elapsed_time(t1)
        
        print "# Saving png with ellipses"
        pylab.savefig(self.pngfile_ell,dpi=dpi)
        print "# Done"
        pylab.close()
        return

    # Get the parameters for the ellipses from the catalogs using and SQL query in the database
    def query_cat_ellipses(self):

        t0 = time.time()

        # Do some magic to extract the integer describing the expnum
        # 'DECam_00220340' -->  220340
        expnum = int(self.exposureID.split('_')[1].lstrip("0"))
        run    = self.runID
        
        #run    = '20130902212714_20130831'
        #expnum = 229342
        
        cur = self.dbh.cursor()

        # ------------------------------------------------------------- 
        # 1.- Collect the catalog id corresponding to the run/exposure/ccd
        # -------------------------------------------------------------                 
        t1 = time.time()
        query1 = """select c.id from image i, exposure e, catalog c 
           where
            i.run         = '%s' and
            e.expnum      = %s and
            i.exposureid  = e.id and
            c.exposureid  = e.id and
            c.parentid    = i.id and
            c.run         = i.run and
            i.imagetype   = 'red' and
            c.catalogtype = 'red_cat'
            order by i.ccd """ % ( run, expnum )

        print "# Will execute the SQL query:\n********\n** %s\n********" % query1
        # Execute the query
        cur.execute(query1)
        print "# Query run/exposure/ccd time  : %s" % elapsed_time(t1)
        # Make a comma-separated list of the catalog's ID that correspond to the exposure/images/run we need
        cat_list = numpy.array(cur.fetchall()).flatten()
        qcatlist = ','.join(map(str, cat_list))

        # ------------------------------------------------------------- 
        # 2.- Gather the SEx entries from the catalogs we need to get
        # ------------------------------------------------------------- 
        t2 = time.time()
        querylist_ell  = "c.ra, c.dec, c.a_image, c.b_image, c.kron_radius, c.theta_image"
        query2 = """select %s
                from objects_current c where c.catalogid in
                (%s) """ % ( querylist_ell, qcatlist )
        print "# Will execute the SQL query:\n********\n** %s\n********" % query2

        cur.execute(query2)
        qout = numpy.array(cur.fetchall())
        self.ra          = qout[:,0]
        self.dec         = qout[:,1]
        self.kron_radius = qout[:,4]
        self.a_image     = qout[:,2]*self.kron_radius
        self.b_image     = qout[:,3]*self.kron_radius
        self.theta       = qout[:,5]
        print "# Catalog query time  : %s" % elapsed_time(t2)
        print "# Total SQL query time: %s" % elapsed_time(t0)
        
        return


    def make_png_thumbnail(self):

        """ Make the thumbnail of the PNG for the webpages"""        

        self.TN_png = "%s_TN.png"     % self.base_outname
        if os.path.exists(self.TN_png) and not self.force:
            print "# File exists -- Skipping creation of %s" % self.TN_png
        else:
            print "# Creating TN %s" % self.TN_png
            im = Image.open(self.pngfile)
            (w,h) = im.size # Get width (w) and height (h)
            TNshape = self.TNsize,int(self.TNsize*h/w)
            im.thumbnail(TNshape, Image.ANTIALIAS)
            im.save(self.TN_png,"png",options='optimize')
        return


    def make_ell_thumbnail(self):

        """ Make TN for the PNG with ellipses"""

        self.TN_ell = '%s_ell_TN.png' % self.base_outname
        if os.path.exists(self.TN_ell) and not self.force:
            print "# File exists -- Skipping creation of %s" % self.TN_ell
        else:
            print "# Creating TN %s" % self.TN_ell
            im = Image.open(self.pngfile_ell)
            (w,h) = im.size # Get width (w) and height (h)
            TNshape = self.TNsize,int(self.TNsize*h/w)
            im.thumbnail(TNshape, Image.ANTIALIAS)
            im.save(self.TN_ell,"png",options='optimize')
        return


    def cross_RA_zero(self):

        """
        Check that we are not crossing the RA=0.0 by comparing the
        distance in RA of oppostive CCDs making sure is not
        greated than the DEcam width.

        Alernative solutions to scanning just a few CCD images:

        1) use the function drawDECam.createDECam_TANheader() to
        generate a fake TAN for the image. It's quick, but if the WCS
        of the image has changed a lot, might not be accurate

        2) loop over all of the wcs of the images... slower (a few secs) but safer.

        """

        d2r = math.pi/180. # degrees to radians shorthand
        DECam_width = 2.3  # approx width in degrees

        # CHANGE: Hard-coded for now, should fix to more general code
        # Center of chips 2048x4096
        xo = 1024
        yo = 4096

        # Figure out the min and max RA of the image, looping trough all of them
        print "# Figuring out edges of CCDs"
        foundRA1 = False
        foundRA2 = False
        foundCenter = False
        t0 = time.time()
        for filename in self.outlist:
            hdr = pyfits.getheader(filename)

            if hdr['CCDNUM'] == 25 or hdr['CCDNUM'] == 32 or not foundRA1:
                wcs = wcsutil.WCS(hdr)
                ra1,dec1 = wcs.image2sky(xo,yo)
                foundRA1 = True 
                
            if hdr['CCDNUM'] == 31 or hdr['CCDNUM'] == 38 or not foundRA2:
                wcs = wcsutil.WCS(hdr)
                ra2,dec2 = wcs.image2sky(xo,yo)
                foundRA2 = True 

            # And the center image
            if hdr['CCDNUM'] == 28:
                self.centerCCD = filename
                foundCenter = True

        if not foundRA1 or not foundRA1:
            print "# Warning: could not find exposure edges"
            
        print "# Cross RA examination: %s" % elapsed_time(t0)
        dec = (dec1+dec2)/2.0
        D = abs(ra1-ra2)

        # in case ras near 360 have negative values
        # so, we need to check that both ras have the same sign
        if math.copysign(1,ra1) == math.copysign(1,ra2):
            same_sign = True
        else:
            same_sign = False

        if D > 2*DECam_width/math.cos(dec*d2r) or same_sign is False:

            self.crossRA = True
            print "# *****************************************"
            print "# **  WARNING: exposure crosses RA=0.0   **"
            print "# *****************************************"
        else:
            self.crossRA = False
        return 
     
    def get_exposure_imsize_center(self):

        """ Get the center and size of the output projected exposure image if required """

        # 1. Get the right size
        # The default size of a projected image in the native pixscale of DECam is 31123x28149 pixels
        NX = 31123 #
        NY = 28149 # at 0.263 arcsec/pixel
        self.NX = int(NX*self.INPUT_PIXEL_SCALE/self.pixscale)
        self.NY = int(NY*self.INPUT_PIXEL_SCALE/self.pixscale)
        # 2. Get the center in RA,Dec
        xo = 2048+104 # 208 is the space between chips in pixels
        yo = 2048
        print "# Will recompute the center using %s" % self.centerCCD
        hdr = pyfits.getheader(self.centerCCD)
        wcs = wcsutil.WCS(hdr)         # use wcsutils
        self.RA0,self.DEC0 = wcs.image2sky(xo,yo)
        return

def elapsed_time(t1,verb=False):
    import time
    t2    = time.time()
    stime = "%dm %2.2fs" % ( int( (t2-t1)/60.), (t2-t1) - 60*int((t2-t1)/60.))
    if verb:
        print >>sys.stderr,"Elapsed time: %s" % stime
        
    return stime

# Dummy function to call in multiprocess
def work(cmd):
    return subprocess.call(cmd, shell=False)

# Check if executable is in path of user
def inpath(program,verb=None):
    """ Checks if program is in the user's path """
    import os.path
    for path in os.environ['PATH'].split(':'):
        if os.path.exists( os.path.join(path,program) ):
            if verb: print "# program: %s found in: %s" % (program , os.path.join(path,program))
            return 1
    if verb: print "# program: %s NOT found in user's path " % program
    return 0


def cmdline():

    """ Parse the command line arguments and options using argparse"""

    import argparse

    USAGE = "\n"
    USAGE = USAGE + "  %(prog)s <path-to-run> <exposure-name> [options] \n" 
    USAGE = USAGE + "  i.e.: \n"
    USAGE = USAGE + "  %(prog)s /archive_data/Archive/OPS/red/20130908135040_20130907/red/ DECam_00231577\n"

    parser = argparse.ArgumentParser(usage=USAGE,
                                     epilog="Author: Felipe Menanteau, NCSA/University of Illinois (felipe@illinois.edu)",
                                     description="Projects a DECam exposure using SWarp and creates grayscale PNGs using Python's native PIL/Image module",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The positional arguments
    parser.add_argument("runID", 
                        help="RunID path (path-mode) or name (SQL-mode).")
    parser.add_argument("exposureID",   
                        help="DECam exposureID.")
    parser.add_argument("--outdir", action="store", default='projected',
                        help="Output Directory to put files in")
    parser.add_argument("--noPNG", action="store_true", default=False,
                        help="Skip creation of PNG files")
    parser.add_argument("--noEll", action="store_true", default=False,
                        help="Skip creation of PNG file with overlaid ellipse")
    parser.add_argument("--MP", action="store_true", default=True,
                        help="Use multiprocess")
    parser.add_argument("--noMP", action="store_true", default=False,
                        help="Don't use multiprocess")
    parser.add_argument("--TNsize", type=int,default=800,
                        help="Size of thumbnail [pixels]")
    parser.add_argument("--grayscale", type=float,default=0.98,
                        help="grayscale for png creation [i.e. 0.98]")
    parser.add_argument("--pixscale",  type=float, action="store", default=1.0,
                        help="pixel-scale in arcsec/pix")
    parser.add_argument("--noSWarp", action="store_true",default=False,
                        help="No SWarp -- dry run")
    parser.add_argument("--weight_thresh", type=float, action="store",default=1e-4,
                        help="SWarps WEIGHT_THRESH value")
    parser.add_argument("--force", action="store_true", default=False,
                        help="Forces the re-creation of existing files")
    parser.add_argument("--dryrun", action="store_true", default=False,
                        help="Dry Run -- only build the lists")
    parser.add_argument("--keepfiles", action="store_true", default=False,
                        help="Keep each CCD projected file")
    parser.add_argument("--noBack", action="store_true", default=False,
                        help="Avoids Background substraction on SWarp call")
    parser.add_argument("--noCopy", action="store_true", default=False,
                        help="Avoids copying input files")
    
    args = parser.parse_args()

    # noPNG turns of noEll
    if args.noPNG:
        args.noEll = 1
    # One turns off the other
    if args.noMP:
        args.MP = False
    # Dry run turns off everything... almost
    if args.dryrun:
        args.noPNG   = True
        args.noEll   = True
        args.noSWarp = True

    print "# Will run:"
    print "# %s " % parser.prog
    for key,val in sorted(vars(args).items()):
       print "# \t--%-10s\t%s" % (key,val)
       
    return args

