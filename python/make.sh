#!/bin/bash
mkdir src
cp -r ../src/ src/
python setup.py install
rm -r src  
