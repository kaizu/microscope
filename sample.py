import csv
import math

import numpy as np
import microscope


def generate_random_points(points, N_point, L):
    assert N_point > 200

    for i in range(200):
        points[i][0] = np.random.uniform(-L_2, +L_2)
        points[i][1] = np.random.uniform(-L_2, +L_2)
        # points[i][2] = 0.0

    for i in range(200, N_point):
        points[i][0] = np.random.uniform(-L_2, +L_2)
        points[i][1] = np.random.uniform(-L_2, +L_2)
        points[i][2] = np.random.uniform(0, 3000)

def read_input(filename, points, data_size, shift, scale):
    i = 0
    with open(filename, 'r') as fin:
        header = fin.readline()
        for row in csv.reader(fin):
            if i >= data_size:
                break

            x, y, z, r, sid = (
                float(row[0]), float(row[1]), float(row[2]), float(row[3]), int(row[4]))

            # if sid != 0:
            #     continue

            points[i][0] = (x - shift[0]) * scale
            points[i][1] = (y - shift[1]) * scale
            points[i][2] = math.fabs(z - shift[2]) * scale
            i += 1

def plot_data(data, N_pixel):
    import matplotlib.pylab as plt

    plt.imshow(data.reshape((N_pixel, N_pixel)), interpolation='none')
    # plt.xlim(270, 330)
    # plt.ylim(270, 330)
    plt.colorbar()
    # plt.savefig("%s.png" % filename)
    plt.show()


if __name__ == "__main__":
    import sys

    N_pixel = 600
    k = 2 * np.pi / 508.0
    N_A = 1.4
    pixel_length = 6500.0 / 100
    focal_point = np.zeros(3, dtype=np.float64)
    cutoff = 2000

    L = pixel_length * N_pixel
    L_2 = L * 0.5

    data = np.zeros(N_pixel * N_pixel, dtype=np.float64)

    if len(sys.argv) == 1:
        N_point = 2620
        points = np.zeros((N_point, 3), dtype=np.float64)
        generate_random_points(points, N_point, L)

        intensity = np.zeros(N_point, dtype=np.float64)
        microscope.emission(points, intensity)

        for i in range(N_point):
            microscope.overlay_psf(
                data, N_pixel, pixel_length, points[i], intensity[i],
                focal_point, k, N_A, cutoff)
    else:
        shift = np.array([15, 15, 1.5], dtype=np.float64)
        scale = 1000.0

        N_point = 2620

        points = np.zeros((N_point, 3), dtype=np.float64)

        filenames = sys.argv[1: ]
        for filename in filenames:
            read_input(filename, points, N_point, shift, scale)

            intensity = np.zeros(N_point, dtype=np.float64)
            microscope.emission(points, intensity)

            for i in range(N_point):
                I = intensity[i] * len(filenames)
                microscope.overlay_psf(
                    data, N_pixel, pixel_length, points[i], intensity[i],
                    focal_point, k, N_A, cutoff)

    microscope.detection(data, data, N_pixel * N_pixel)

    plot_data(data, N_pixel)
