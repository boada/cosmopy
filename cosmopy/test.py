#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import range
from past.utils import old_div
import numpy
import cosmopy
import sys

Ho = 70.0
Om = 0.3
OL = 0.7

dz = 0.2
z = numpy.arange(0.0, 2.1, dz)  # .astype('Float')

# Set the morphological parameters
flat = (Om, OL, old_div(Ho, 100.))
c = cosmopy.set(flat)

dl = c.dlum(z)
da = c.dang(z)
age = old_div(c.age(z), 1.e9)
lbt = c.lookback(z)

print(" --------------------------------------------------------------------")
print(" cosmology: Om=%.2f, OL=%.2f,  h0=%.2f" % flat)
print("%6s %12s %12s %12s %12s" % ('z', 'dl(z)', 'da(z)', 'age(z)[Gyr]',
                                   'lookbacktime(z)'))
print(" --------------------------------------------------------------------")
for i in range(len(z)):
    print("%6.2f %12.5f %12.5f %.6e %.6e" % (z[i], dl[i], da[i], age[i],
                                             lbt[i]))

open = (1., 0., old_div(Ho, 100.))
o = cosmopy.set(open)

dl = o.dlum(z)
da = o.dang(z)
age = old_div(o.age(z), 1.e9)
lbt = o.lookback(z)

print(
    " ----------------------------------------------------------------------")
print(" cosmology: Om=%.2f, OL=%.2f,  h0=%.2f" % open)
print("%6s %12s %12s %12s %12s" % ('z', 'dl(z)', 'da(z)', 'age(z)[Gyr]',
                                   'lookbacktime(z)'))
print(
    " ----------------------------------------------------------------------")
for i in range(len(z)):
    print("%6.2f %12.5f %12.5f %.6e %.6e" % (z[i], dl[i], da[i], age[i],
                                             lbt[i]))
