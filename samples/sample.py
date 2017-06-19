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

            points[i][0] = (y - shift[0]) * scale
            points[i][1] = (z - shift[1]) * scale
            points[i][2] = math.fabs(x - shift[2]) * scale + 650
            # points[i][0] = (x - shift[0]) * scale
            # points[i][1] = (y - shift[1]) * scale
            # points[i][2] = math.fabs(z - shift[2]) * scale + 650
            i += 1

def plot_data(data, N_pixel):
    import matplotlib.pylab as plt

    #plt.imshow(data.reshape((N_pixel, N_pixel)), interpolation='none', cmap=plt.cm.gray)
    plt.imshow(data.reshape((N_pixel, N_pixel)), interpolation='none', cmap=plt.cm.viridis)
    # plt.xlim(270, 330)
    # plt.ylim(270, 330)
    plt.colorbar()
    plt.savefig("result.txt.png")
    plt.show()

def overlay_psf_from_files(
        filenames, N_pixel, N_point, pixel_length, focal_point, k, N_A, cutoff):
    data = np.zeros(N_pixel * N_pixel, dtype=np.float64)

    shift = np.array([15, 15, 0.5], dtype=np.float64)
    scale = 1000.0

    points = np.zeros((N_point, 3), dtype=np.float64)
    intensity = np.zeros(N_point, dtype=np.float64)

    for filename in filenames:
        print('=> {}'.format(filename))

        read_input(filename, points, N_point, shift, scale)
        microscope.emission(points, intensity)

        for i in range(N_point):
            I = intensity[i]
            microscope.overlay_psf(
                data, N_pixel, pixel_length, points[i], intensity[i],
                focal_point, k, N_A, cutoff)
    return data

def pool_function(*args):
    try:
        return overlay_psf_from_files(args[-1], *args[: -1])
    except KeyboardInterrupt:
        return None


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--np", help="the number of processes", type=int)
    parser.add_argument(
        'filenames', metavar='FILENAME', type=str, nargs='*', help='an input filename')
    args = parser.parse_args()
    proc = 1 if args.np is None else args.np
    filenames = args.filenames

    N_pixel = 600
    k = 2 * np.pi / 508.0
    N_A = 1.4
    pixel_length = 6500.0 / 100
    focal_point = np.zeros(3, dtype=np.float64)
    cutoff = 2000

    L = pixel_length * N_pixel
    L_2 = L * 0.5

    data = np.zeros(N_pixel * N_pixel, dtype=np.float64)

    if len(filenames) == 0:
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
        import multiprocessing as mp
        import functools

        # N_point = 17200
        # N_point = 2620
        # N_point = 1000
        N_point = 1000

        pool = mp.Pool(proc)
        m = int(math.ceil(len(filenames) / float(proc)))

        try:
            result = pool.map_async(
                functools.partial(
                    pool_function,
                    N_pixel, N_point, pixel_length, focal_point, k, N_A, cutoff),
                [filenames[i * m: (i + 1) * m] for i in range(proc)])
        except KeyboardInterrupt as e:
            pool.terminate()
            raise e
        else:
            pool.close()
            pool.join()

        data = sum(result.get()) / len(filenames)

    microscope.detection(data, data, N_pixel * N_pixel)
    np.savetxt("result.txt", data)

    plot_data(data, N_pixel)
