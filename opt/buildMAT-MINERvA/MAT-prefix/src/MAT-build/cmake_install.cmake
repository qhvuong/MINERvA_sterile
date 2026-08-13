# Install script for directory: /exp/minerva/app/users/qvuong/MAT_AL9/MAT

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/exp/minerva/app/users/qvuong/MAT_AL9/opt")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/opt/spack/linux-almalinux9-x86_64_v3/gcc-11.4.1/binutils-2.42-exh76ctar6d7u2dbwoflh74fapgrrkra/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE FILE FILES "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/setup_MAT.sh")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/MAT" TYPE FILE FILES
    "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/MATConfig.cmake"
    "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/MATConfigVersion.cmake"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/PlotUtils/cmake_install.cmake")
  include("/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/macros/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
