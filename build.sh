#!/bin/bash
cd ..
rm -rf coseg-release
cp -r coseg coseg-release
grep -v cosegData coseg/DESCRIPTION > coseg-release/DESCRIPTION
rm coseg-release/tests/testthat/test-cosegData.R
PKG_TGZ=$(R CMD build coseg-release|grep building|sed 's/.*‘//'|sed 's/’.*//')
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
