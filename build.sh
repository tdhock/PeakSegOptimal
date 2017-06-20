#!/bin/bash
cd ..
rm -rf coseg-release
cp -r coseg coseg-release
PKG_TGZ=$(R CMD build coseg-release|grep coseg_|sed 's/.*‘//'|sed 's/’.*//')
R CMD check --as-cran $PKG_TGZ
