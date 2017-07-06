
cosmopy is python package that computes the many cosmological
distances for Lambda non-zero cosmology, such as the luminosity,
angular and proper distances at a given redshift. It also computes the
age of the universe, the look-back time and the age of galaxy for a
given redshift of formation. cosmopy follows the definitions from Hogg
(1998) astro-ph/9905116. One of the main advantages of cosmopy is that
it works with arrays, so it can compute dl(z) where z(n) is very big
array. It numerically integrates in E(z) integral using a fast Fortran
77 library.


Requirements:
-------------

1) First of all, you need python 2.2 or greater
   http://www.python.org. If you don;t use python, I recommend to read
   one the tutorials from their website.

2) A Fortran 77 compiler ( i.e. g77/gcc will do the trick)

3) Numpy. Get it from http://numpy.scipy.org/


Installation:
-------------

a) Build and install system-wide


%> ./setup.py build install 

or 

%> ./setup.py build install --prefix=/usr/local


b) Build and install in the user's python repository


%> ./setup.py build install --home=/home/felipe/Python
(where /home/felipe/Python is the repository)

In case you need to pass special flags to the compiler use the 
config_fc option

%> ./setup.py config_fc "--f77exec=/usr/bin/f77" install --home=/home/felipe/Python

Examples:
---------

# Fits import the necesary packages
import numpy
import cosmopy

Ho = 70.0
Om = 0.3
OL = 0.7
   
dz = 0.2
z  = numarray.arange(0.0,2.1,dz)

# Set the morphological parameters
flat = (Om, OL, Ho/100.)	
c    = cosmopy.set(flat)
	
# The available functions are:

# Luminosity distance
c.dlum(z)

# Proper distace
c.proper_distance(z)

# Angular distance
c.dang(z)

# Volumen element at z
c.dvol_comov(z)

# Volumen at redshift z
c.vol_comov(z)

# Volumen between z1 and z2 (z2 > z1!!)
c.vol_comov12(z1,z2)

# Age of the universe at z
c.age(z)

# Lookback time at z
c.lookback(z)

# Finds the redshift of a given age
c.get_redshift(ages,units='yr')

# The age of galaxy for a give z and z_formation
c.age_galaxy(z,zf=5.0)

#
# If we now want to try a diferent cosmology
open = (1.,0.,Ho/100.)
o    = cosmopy.set(open)

# The distances for the new cosmology are
o.dlum(z)
o.dang(z)
o.age(z)
o.lookback(z)
# ... and so on...
