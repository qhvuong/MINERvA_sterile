# Install script for directory: /exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE FILES "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/PlotUtils/MATDict.rootmap")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE FILES "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/PlotUtils/MAT_rdict.pcm")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMAT.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMAT.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMAT.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/PlotUtils/libMAT.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMAT.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMAT.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMAT.so"
         OLD_RPATH "/cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/opt/spack/linux-almalinux9-x86_64_v3/gcc-12.2.0/root-6.28.12-hljl7gyomotoqht2uzvhnf73337jq67q/lib/root:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/opt/spack/linux-almalinux9-x86_64_v3/gcc-11.4.1/binutils-2.42-exh76ctar6d7u2dbwoflh74fapgrrkra/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMAT.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  EXECUTE_PROCESS( COMMAND /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/opt/spack/linux-almalinux9-x86_64_v3/gcc-11.4.1/cmake-3.27.9-dxlizspondhgno5v6vg5wpmvtkapxkyb/bin/cmake -E create_symlink libMAT.so /exp/minerva/app/users/qvuong/MAT_AL9/opt/lib/libMATDict.so )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/PlotUtils" TYPE FILE FILES
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/BaseUniverse.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/ChainWrapper.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/CrashOnROOTMessage.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/Cut.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/Cutter.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/ErrorHandler.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/Exceptions.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/FluxReweighter.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/FluxSystematics.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/GridCanvas.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/Hist2DWrapper.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/HistFolio.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/HistWrapper.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/HistogramUtils.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/Linkdef.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MacroUtil.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvApplication.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvColors.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvFluxConstraint.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvH1D.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvH1DToCSV.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvH2D.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvH2DLog.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvH2DToCSV.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvH3D.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvLatErrorBand.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvLatErrorBand2D.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvLatErrorBand3D.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvNuclearModelWeight.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvPlotter.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvVertErrorBand.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvVertErrorBand2D.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvVertErrorBand3D.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/NamedCategory.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/NuclModUtils.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/ROOTglob.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/Table.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/TreeWrapper.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/Variable2DBase.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/VariableBase.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/genie_particle.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/makeChainWrapper.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/modifier_Eb.h"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/BaseUniverse.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/ChainWrapper.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/Cutter.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/ErrorHandler.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/FluxReweighter.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/FluxSystematics.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/GridCanvas.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/Hist2DWrapper.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/HistFolio.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/HistWrapper.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/HistogramUtils.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MacroUtil.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvApplication.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvColors.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvFluxConstraint.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvH1D.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvH1DToCSV.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvH2D.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvH2DLog.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvH2DToCSV.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvH3D.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvLatErrorBand.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvLatErrorBand2D.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvLatErrorBand3D.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvNuclearModelWeight.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvPlotter.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvVertErrorBand.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvVertErrorBand2D.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/MnvVertErrorBand3D.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/NuclModUtils.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/ROOTglob.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/Table.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/TreeWrapper.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/Variable2DBase.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/VariableBase.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/genie_particle.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/makeChainWrapper.cxx"
    "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/modifier_Eb.cxx"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/MAT/MATTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/MAT/MATTargets.cmake"
         "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/PlotUtils/CMakeFiles/Export/31467188b90ce0aecf35f9e1df0f7389/MATTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/MAT/MATTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/MAT/MATTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/MAT" TYPE FILE FILES "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/PlotUtils/CMakeFiles/Export/31467188b90ce0aecf35f9e1df0f7389/MATTargets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/MAT" TYPE FILE FILES "/exp/minerva/app/users/qvuong/MAT_AL9/opt/buildMAT-MINERvA/MAT-prefix/src/MAT-build/PlotUtils/CMakeFiles/Export/31467188b90ce0aecf35f9e1df0f7389/MATTargets-release.cmake")
  endif()
endif()

