from __future__ import print_function
from __future__ import division
# Set of routines to compute a number of
# useful cosmological distances

# Changed fron numarry to numpy -- Aug 2009

from builtins import object
from past.utils import old_div
import numpy
import sys
from numpy import sqrt as sqrt
from numpy import sin as sin
from numpy import sinh as sinh
from numpy import arcsin as asin
from numpy import arcsinh as asinh
from . import fcosmo


class set(object):
    def __init__(self, cosmology=(0.3, 0.7, 0.7)):

        c = 299792.458  # speed of light in km/s
        Mpc2km = 3.08568025e19  # kms in 1 Mpc
        yr2sec = 31536000.  # seconds in 1 yr
        self.cosmo = cosmology
        self.Om = float(cosmology[0])
        self.OL = float(cosmology[1])
        self.Ok = 1.0 - self.Om - self.OL
        self.h0 = float(cosmology[2])
        self.dh = old_div(c, (100. * self.h0))
        self.th = 100. * self.h0 * yr2sec / Mpc2km
        self.th = old_div(1.0, self.th)  # in Hubble time in yr^-1
        # Set the fortran function params
        self.set_cosmo()
        return

    def set_cosmo(self):
        fcosmo.pars.ol = self.OL
        fcosmo.pars.om = self.Om
        fcosmo.pars.ok = self.Ok
        return

    def dlum(self, z):
        ''' Return the luminosity distance and redshift z
        for a give cosmology'''
        # We use the transverse comoving distance
        self.comoving_distance_tranverse(z)
        self.dl = (1. + z) * self.DM
        return self.dl

    # A short hand for the comoving distance
    def dc(self, z):
        return self.comoving_distance_tranverse(z)

    def comoving_distance_tranverse(self, z):
        self.set_cosmo()
        ''' Returns the transverse comoving distance (DM) as in Hogg
        (1998) equations 16 and 17'''
        z = numpy.asarray(z)

        # Simple case OL = 0
        if (self.OL == 0 and self.Om != 0.):
            self.DM = self.dh * 2. * (
                2. - self.Om * (1. - z) -
                (2. - self.Om) * sqrt(1. + self.Om * z)) / (self.Om**2 *
                                                            (1. + z))
            return self.DM

        # Flat cosmo Om + OL = 1
        if (self.Ok == 0):
            self.DM = self.dh * fcosmo.dist_comov(z)
            return self.DM
        # Ok not zero
        else:
            DC = fcosmo.dist_comov(z)
            kappa = sqrt(abs(self.Ok))

            # Open
        if (self.Ok > 0):
            self.DM = self.dh * (old_div(1., kappa)) * sinh(kappa * DC)
            return self.DM
        # Closed
        if (self.Ok < 0):
            self.DM = self.dh * (old_div(1., kappa)) * sin(kappa * DC)
            return self.DM

        print("# ERROR:Could not compute comoving distance", file=sys.sdterr)
        sys.exit()
        return

    def dang(self, z):
        ''' Angular distance'''
        return old_div(self.dlum(z), (1. + z)**2)

    def dvol_comov(self, z):
        ''' Comoving volume element dv/dz/dAngle Mpc^3/dz/dangle'''
        self.dlum(z)
        Ez = self.ez(z)
        self.dvol = self.dh * self.dl**2 / ((1. + z)**2 * Ez)
        return self.dvol

    def vol_comov(self, z):
        ''' Comoving volumen element at z in Mpc^3 as using Hogg
        (1998) eq. 29 for Ok=0 and Ok not 0'''

        pi = 3.14159265
        # Get the comoving distance tranverse first
        DM = old_div(self.comoving_distance_tranverse(z), self.dh)

        # Flat cosmology
        if (self.Ok == 0):
            return (4. * pi / 3.) * self.DM**3

        # Ok not zero
        else:
            kappa = sqrt(abs(self.Ok))
            f = 4.0 * pi * self.dh**3 / (2. * self.Ok)

        # Open
        if (self.Ok > 0):
            return f * DM * sqrt(1. + self.Ok * DM**2) - (old_div(1., kappa)
                                                          ) * asinh(kappa * DM)

        # Close
        if (self.Ok < 0):
            return f * DM * sqrt(1. + self.Ok * DM**2) - (old_div(1., kappa)
                                                          ) * asin(kappa * DM)

        return

    def ez(self, z):
        return sqrt(self.Om * (1. + z)**3 + self.Ok * (1. + z)**2 + self.OL)

    def vol_comov12(self, z1, z2):
        z1 = numpy.asarray(z1)
        z2 = numpy.asarray(z2)
        '''Returns the comoving volumen between z1 and z2 in Mpc^3 '''

        self.DM = None
        vol1 = self.vol_comov(z1)
        self.DM = None
        vol2 = self.vol_comov(z2)

        return vol2 - vol1

    def age(self, z):
        ''' return the sge of universe at redshift z '''
        z = numpy.asarray(z)
        if (self.Ok == 0.0 and self.OL != 0.0):
            A = sqrt(self.Om * (1. + z)**3 + self.OL) + sqrt(self.OL)
            B = sqrt(self.Om * (1. + z)**3)
            C = self.th * 2. / 3. / sqrt(self.OL)
            return C * numpy.log(old_div(A, B))

        self.set_cosmo()
        return self.th * fcosmo.time(z)

    def lookback(self, z):
        ''' return the look backtime at redshift z'''
        z = numpy.asarray(z)
        if (self.Ok == 0.0 and self.OL != 0.0):
            self.to = self.age(0.0)
            A = sqrt(self.Om * (1. + z)**3 + self.OL) + sqrt(self.OL)
            B = sqrt(self.Om * (1. + z)**3)
            C = self.th * 2. / 3. / sqrt(self.OL)
            return self.to - C * numpy.log(old_div(A, B))

        self.set_cosmo()
        return self.th * fcosmo.lkbt(z)

    def get_redshift(self, ages, units='yr'):
        ''' Finds the redshift for a given age'''
        ages = numpy.asarray(ages)
        self.set_cosmo()
        if units is "Gyr" or units is "gyr":
            ages = ages / 1.e9 / self.th
        elif units is "yr":
            ages = old_div(ages, self.th)
        else:
            print("# error: get_redshift() specify units of age [Gyr] or [yr]",
                  file=sys.stderr)
            sys.exit()
        zx = fcosmo.get_z(ages)
        return zx

    def age_galaxy(self, z, zf=5.0):
        z = numpy.asarray(z)
        return self.age(z) - self.age(zf)
