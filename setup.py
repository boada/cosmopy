#!/usr/bin/env python

from numpy.distutils.core import Extension

#
# Construct the ext for fcosmo.f
#
ext = Extension(name='cosmopy.fcosmo',
                sources=['cosmopy/fcosmo.f'], )

ldesc = "cosmopy is a python package that computes the many cosmological\n\
         distances for Lambda non-zero cosmology, such as the\n\
         luminosity, angular and comoving distances at a given\n\
         redshift. It also computes the age of the universe, the\n\
         look-back time and the age of galaxy for a given redshift of\n\
         formation. cosmopy follows the definitions from Hogg (1998)\n\
         astro-ph/9905116."

if __name__ == "__main__":

    from numpy.distutils.core import setup
    setup(
        name='cosmopy',
        description="A Python class for Cosmological Distances",
        author="Felipe Menanteau",
        author_email="felipe@physics.rutgers.edu",
        version="0.2",
        url="http://felipe.menanteau.org",
        license="GNU General Public License, http://www.gnu.org/copyleft/gpl.html",
        platforms=["Linux", "Solaris"],
        long_description=ldesc,
        packages=['cosmopy'],
        ext_modules=[ext])
