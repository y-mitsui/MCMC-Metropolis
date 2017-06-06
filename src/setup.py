from distutils.core import setup, Extension
from Cython.Build import cythonize
import os

ext1 = Extension("mixture_bernoulli_wrap",
                sources=["mixture_bernoulli_wrap.pyx", "mixture_bernoulli.c", "common.c", "metropolis.c", "random.c"],
                extra_compile_args=["-I/usr/local/include/"],
                extra_link_args=["-L/usr/local/lib", "-lgsl", "-lgslcblas", "-lm"],
                )
                
setup(name = 'mixture_bernoulli_wrap', ext_modules = cythonize([ext1]))