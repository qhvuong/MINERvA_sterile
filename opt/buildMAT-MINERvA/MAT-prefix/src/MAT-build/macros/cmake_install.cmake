# Install script for directory: /exp/minerva/app/users/qvuong/MAT_AL9/MAT/macros

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/madd" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/madd")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/madd"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/macros/madd")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/madd" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/madd")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/madd"
         OLD_RPATH "/cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/opt/spack/linux-almalinux9-x86_64_v3/gcc-12.2.0/root-6.28.12-hljl7gyomotoqht2uzvhnf73337jq67q/lib/root:/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/PlotUtils:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/opt/spack/linux-almalinux9-x86_64_v3/gcc-11.4.1/binutils-2.42-exh76ctar6d7u2dbwoflh74fapgrrkra/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/madd")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE FILE FILES
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/macros/ProcessMCSampleSizeScan.py"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/macros/PlotMCSampleSizeScan.py"
    )
endif()

