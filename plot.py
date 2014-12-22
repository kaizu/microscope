import numpy
import sys
import matplotlib.pylab as plt

filename = sys.argv[1]
data = numpy.loadtxt(filename)
Nsq = len(data)
N = numpy.sqrt(Nsq)
plt.imshow(data.reshape((N, N)))
plt.colorbar()
plt.show()
