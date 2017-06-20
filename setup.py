#!/usr/bin/env python

"""
setup.py file for pyTFMPval
"""

from distutils.core import setup, Extension


pytfmpval_module = Extension('_pytfmpval',
                             sources=['src/Matrix.cpp', 'pytfmpval_wrap.cxx'],
                             swig_opts=['-c++'],
                             language='c++'
                             )

setup(name='pytfmpval',
      version='0.1',
      author="Jared Andrews",
      author_email='jared.andrews07@gmail.com',
      description="""Python binding for the TFMPvalue program.""",
      ext_modules=[pytfmpval_module],
      py_modules=["pytfmpval"]
      )
