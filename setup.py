from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

setup(
    name="microscope",
    cmdclass={"build_ext": build_ext},
    ext_modules=[
        Extension("microscope", ["microscope.pyx"], language="c++",
                  libraries=["gsl", "gslcblas"],
                  include_dirs=[numpy.get_include()])
        ]
)
