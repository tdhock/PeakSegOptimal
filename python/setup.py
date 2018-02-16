#!/usr/bin/env python
import os
import sys
from distutils.core import setup, Extension

os.environ["CC"] = "clang++"

is_darwin = sys.platform=='darwin'

setup(name='FastLZeroSpikeInference',
      version='1.0',
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
                         'src/IsotonicFPOP.cpp', 
                         'src/interface.cpp'
                         ],
              include_dirs = ['src/'], 
              language = 'c++',
              extra_compile_args = [
                  '-std=c++11', 
                  '-stdlib={}'.format('libc++' if is_darwin else 'libstdc++')
                  ]
              )
          ]
     )

