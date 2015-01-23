#!/usr/bin/python

import numpy
import sys
import matplotlib.pylab as plt

if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = "result.txt"
data = numpy.loadtxt(filename)
Nsq = len(data)
N = numpy.sqrt(Nsq)
print N
print data.max(), data.sum()
# plt.imshow(data.reshape((N, N)), cmap=plt.cm.gray, interpolation='none')
plt.imshow(data.reshape((N, N)), interpolation='none')
plt.colorbar()
plt.show()
