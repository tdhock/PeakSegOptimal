#!/usr/bin/env python
import os
import sys
from distutils.core import setup, Extension

is_darwin = sys.platform=='darwin'
is_linux = 'linux' in sys.platform

extra_compile_args = ['-std=c++11']

if is_darwin:
  os.environ["CC"] = 'clang++'
  extra_compile_args.append('-stdlib=libc++')
elif is_linux:
  os.environ["CC"] = 'g++'

setup(name='FastLZeroSpikeInference',
      version='2018.09.29',
      description='Python wrapper for FastLZeroSpikeInference',
      author='Sean Jewell',
      author_email='swjewell@uw.edu',
      url='https://github.com/jewellsean/FastLZeroSpikeInference',
      packages=['FastLZeroSpikeInference'],
      ext_modules=[
          Extension(
              name = 'FastLZeroSpikeInference', 
              sources = ['src/funPieceListLog.cpp', 
                         'src/ARFPOP.cpp', 
                         'src/FitSegmentModel.cpp', 
                         'src/python_interface.cpp'
                         ],
              include_dirs = ['src/'], 
              language = 'c++',
              extra_compile_args = extra_compile_args
              )
          ]
     )
