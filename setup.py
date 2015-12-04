from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

setup(
    name="microscope",
    cmdclass={"build_ext": build_ext},
    ext_modules=[
        Extension("microscope", ["python/microscope.pyx"], language="c++",
                  depends=[
                      "born_wolf_psf.hpp", "born_wolf_psf_table.hpp",
                      "detection.hpp", "emission.hpp", "integration.hpp",
                      "microscope.hpp", "point_spreading_functions.hpp"],
                  libraries=["gsl", "gslcblas"],
                  include_dirs=[".", numpy.get_include()])
        ]
)
