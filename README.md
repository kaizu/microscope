microscope
==========

`libboost-dev` and `libgsl-dev` are required.

```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

```
$ OMP_NUM_THREADS=5 ./main test*.csv
```

```
$ g++ -Icubature-1.0.2 -DCUBATURE main.cpp cubature-1.0.2/hcubature.c -lgsl -lcblas
```

python
------

```
$ mkdir build
$ cd build
$ cmake ..
$ make psf_tables
$ cd python
$ python setup.py build_ext --inplace
$ PYTHONPATH=. python ../../samples/sample.py
```

gpu
---

`cuda` and `cuda-drivers` are also required.

```
$ mkdir build
$ cd build
$ cmake -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-9.1 ..
$ make
$ ./microscope/main_gpu
```

result
------

```
$ python ../samples/plot.py result.txt
```

![doc/result.txt.png](doc/result.txt.png)
