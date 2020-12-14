from distutils.core import setup

# The main call
setup(name='projectDECam',
      version='3.1.0',
      license="GPL",
      description="Created full focal plane fits files and grayscale PNGs of DECam exposures",
      author="Felipe Menanteau",
      author_email="felipe@illinois.edu",
      packages=['projectDECam'],
      package_dir={'': 'python'},
      scripts=['bin/projectDECamPNG', ],
      data_files=[('ups', ['ups/projectDECam.table']),
                  ('etc', ['etc/default.stiff', 'etc/default.swarp'])]
      )
