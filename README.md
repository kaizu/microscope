microscope
==========

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/42c73eb47bf240c38a730689a3b73188)](https://app.codacy.com/app/kwaizu/microscope?utm_source=github.com&utm_medium=referral&utm_content=kaizu/microscope&utm_campaign=Badge_Grade_Dashboard)

`libboost-dev` and `libgsl-dev` are required. `cuda` and `cuda-drivers` are also required for GPU acceleration.

```
$ mkdir build
$ cd build
$ cmake -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-9.1 -DCMAKE_INSTALL_PREFIX=./local ..
$ make
$ make install
```

```
$ ./local/bin/main_gpu
```

```
$ OMP_NUM_THREADS=5 ./local/bin/main test*.csv
```

```
$ PYTHONPATH=./local/lib/python3.5/site-packages python ../samples/sample.py
```

result
------

```
$ python ../samples/plot.py result.txt
```

![doc/result.txt.png](doc/result.txt.png)

more
----

```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

```
$ g++ -Icubature-1.0.2 -DCUBATURE main.cpp cubature-1.0.2/hcubature.c -lgsl -lcblas
```

```
$ mkdir build
$ cd build
$ cmake ..
$ make psf_tables
$ cd python
$ python setup.py build_ext --inplace
$ PYTHONPATH=. python ../../samples/sample.py
```

```
$ mkdir build
$ cd build
$ cmake -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-9.1 ..
$ make
$ ./microscope/main_gpu
```
