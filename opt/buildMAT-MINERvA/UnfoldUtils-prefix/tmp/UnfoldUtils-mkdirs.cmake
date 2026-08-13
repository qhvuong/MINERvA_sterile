# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/exp/minerva/app/users/qvuong/MAT_AL9/MAT-MINERvA/bootstrap/../../UnfoldUtils"
  "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/UnfoldUtils-prefix/src/UnfoldUtils-build"
  "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/UnfoldUtils-prefix"
  "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/UnfoldUtils-prefix/tmp"
  "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/UnfoldUtils-prefix/src/UnfoldUtils-stamp"
  "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/UnfoldUtils-prefix/src"
  "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/UnfoldUtils-prefix/src/UnfoldUtils-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/UnfoldUtils-prefix/src/UnfoldUtils-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/UnfoldUtils-prefix/src/UnfoldUtils-stamp${cfgdir}") # cfgdir has leading slash
endif()
