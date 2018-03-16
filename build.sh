#!/bin/bash
cd ..
rm -rf PeakSegOptimal-release
cp -r PeakSegOptimal PeakSegOptimal-release
grep -v cosegData PeakSegOptimal/DESCRIPTION | grep -v Remotes > PeakSegOptimal-release/DESCRIPTION
rm PeakSegOptimal-release/tests/testthat/test-PeakSegOptimalData.R
PKG_TGZ=$(R CMD build PeakSegOptimal-release|grep building|sed 's/.*‘//'|sed 's/’.*//')
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
