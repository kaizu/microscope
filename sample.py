import numpy as np
import microscope

N_pixel = 600
k = 2 * np.pi / 508.0
N_A = 1.4
pixel_length = 6500.0 / 100
focal_point = np.zeros(3, dtype=np.float64)
cutoff = 2000

N_point = 2620

L = pixel_length * N_pixel
L_2 = L * 0.5

points = np.zeros((N_point, 3), dtype=np.float64)

for i in range(200):
    points[i][0] = np.random.uniform(-L_2, +L_2)
    points[i][1] = np.random.uniform(-L_2, +L_2)
    # points[i][2] = 0.0

for i in range(200, N_point):
    points[i][0] = np.random.uniform(-L_2, +L_2)
    points[i][1] = np.random.uniform(-L_2, +L_2)
    points[i][2] = np.random.uniform(0, 3000)

intensity = np.zeros(N_point, dtype=np.float64)
microscope.emission(points, intensity)

data = np.zeros(N_pixel * N_pixel, dtype=np.float64)

for i in range(N_point):
    microscope.overlay_psf(
        data, N_pixel, pixel_length, points[i], intensity[i],
        focal_point, k, N_A, cutoff)

microscope.detection(data, data, N_pixel * N_pixel)

import matplotlib.pylab as plt
plt.imshow(data.reshape((N_pixel, N_pixel)), interpolation='none')
# plt.xlim(270, 330)
# plt.ylim(270, 330)
plt.colorbar()
# plt.savefig("%s.png" % filename)
plt.show()
