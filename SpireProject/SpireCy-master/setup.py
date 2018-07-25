from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

ext_modules = cythonize([Extension("client", ["cy_files/client.pyx", 'c_files/utils.c', 'c_files/factorizations.c']),
                         Extension("server", ["cy_files/server.pyx"]),
                         Extension("seq", ["cy_files/seq.pyx", 'c_files/utils.c', 'c_files/factorizations.c'])])

setup(
  name='SpireCy',
  cmdclass={'build_ext': build_ext},
  ext_modules=ext_modules
)