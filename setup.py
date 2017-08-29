#!/usr/bin/env python

"""
setup.py file for pytfmpval
"""

from distutils.core import setup, Extension


pytfmpval_module = Extension('_pytfmpval',
                             sources=['src/Matrix.cpp', 'pytfmpval/pytfmpval_wrap.cxx'],
                             swig_opts=['-c++'],
                             language='c++'
                             )

setup(name='pytfmpval',
      version='0.2.0',
      author="Jared Andrews",
      author_email='jared.andrews07@gmail.com',
      url='https://github.com/j-andrews7/pytfmpval',
      description="""Python bindings for the TFM-Pvalue program.""",
      license='GPL-3.0',
      keywords='bioinformatics tfmpvalue motifs transcription factor genomics science',
      install_requires=["psutil"],
      ext_modules=[pytfmpval_module],
      py_modules=["pytfmpval"],
      classifiers=['Development Status :: 4 - Beta',
                   'Intended Audience :: Science/Research',
                   'Intended Audience :: Developers',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3.3',
                   'Programming Language :: Python :: 3.4',
                   'Programming Language :: Python :: 3.5',
                   'Programming Language :: Python :: 3.6',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Mathematics',
                   'Topic :: Software Development :: Libraries :: Python Modules']
      )
