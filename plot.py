#!/usr/bin/python

import numpy
import sys
import matplotlib.pylab as plt

if len(sys.argv) > 1:
    filenames = sys.argv[1: ]
else:
    filenames = ["result.txt"]

for filename in filenames:
    data = numpy.loadtxt(filename)
    Nsq = len(data)
    N = numpy.sqrt(Nsq)
    print N
    print data.max(), data.sum()
    # plt.imshow(data.reshape((N, N)), cmap=plt.cm.gray, interpolation='none')
    plt.imshow(data.reshape((N, N)), interpolation='none')
    # plt.xlim(270, 330)
    # plt.ylim(270, 330)
    plt.colorbar()
    plt.savefig("%s.png" % filename)
    plt.show()
    plt.clf()
