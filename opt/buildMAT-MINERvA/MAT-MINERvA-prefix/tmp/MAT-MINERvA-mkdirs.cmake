# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/exp/minerva/app/users/qvuong/MAT_AL9/MAT-MINERvA/bootstrap/../../MAT-MINERvA"
  "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-MINERvA-prefix/src/MAT-MINERvA-build"
  "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-MINERvA-prefix"
  "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-MINERvA-prefix/tmp"
  "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-MINERvA-prefix/src/MAT-MINERvA-stamp"
  "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-MINERvA-prefix/src"
  "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-MINERvA-prefix/src/MAT-MINERvA-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-MINERvA-prefix/src/MAT-MINERvA-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-MINERvA-prefix/src/MAT-MINERvA-stamp${cfgdir}") # cfgdir has leading slash
endif()
