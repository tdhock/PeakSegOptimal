#!/bin/bash
cp -r ../src ./
python setup.py install
rm -r src  
