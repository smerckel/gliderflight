Installing GliderFlight
=======================


Download
--------

The software's repository is hosted on `github
<https://github.com/smerckel/gliderflight>`_, from where you can
download the source using git or download it as a zip file.

Installing
----------
After having obtained the source code, you can build it using the
standard way of installing python code. On linux this would be ::
  
  $ python3 setup.py build
  $ python3 setup.py install

Depending on your system setup, the ``install`` command may require root privileges.

Dependencies
^^^^^^^^^^^^
The gliderflight module depends on numpy.

The glidertrim script -- useful to check and adjust the ballast trim
of glider during deployment -- additionally depends on scipy,
matplotlib, gsw, and dbdreader. All these packages are available from
PyPi and will be downloaded automatically when not present. However,
you may prefer to install the complex packages numpy, scipy and
matplotlib using your distribution's package manager (when on linux).

See also the file requirements.txt in the root directory.

PyPi
----
Glider flight is available from pypi ::

  pip install gliderflight


