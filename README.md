microscope
==========

```
$ cmake .
$ make
$ ./main
$ python plot.py
```

```
$ g++ generate_table.cpp -lgsl -lcblas
$ ./a.out
$ g++ main.cpp -lgsl -lcblas
$ ./a.out
$ python plot.py
```

```
$ g++ -Icubature-1.0.2 -DCUBATURE main.cpp cubature-1.0.2/hcubature.c -lgsl -lcblas
```

```
$ OMP_NUM_THREADS=5 ./main test*.csv
```

```
$ cmake .
$ make psf_tables
$ python setup.py build_ext --inplace
$ python sample.py
```

gpu
---

`libboost-dev`, `libgsl-dev`, `cuda` and `cuda-drivers` are required.

```
$ mkdir build
$ cd build
$ cmake -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-9.1 ..
$ make
$ ./microscope/main
$ python ../samples/plot.py result.txt
```

![doc/result.txt.png](doc/result.txt.png)
