# Install script for directory: /home/gulsum/galaxy/tools/racon/vendor/ClaraGenomicsAnalysis

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/gulsum/galaxy/tools/racon/install")
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
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/docs" TYPE DIRECTORY FILES "/home/gulsum/galaxy/tools/racon/ClaraGenomicsAnalysis/html/")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/gulsum/galaxy/tools/racon/ClaraGenomicsAnalysis/common/logging/cmake_install.cmake")
  include("/home/gulsum/galaxy/tools/racon/ClaraGenomicsAnalysis/common/utils/cmake_install.cmake")
  include("/home/gulsum/galaxy/tools/racon/ClaraGenomicsAnalysis/common/io/cmake_install.cmake")
  include("/home/gulsum/galaxy/tools/racon/ClaraGenomicsAnalysis/cudapoa/cmake_install.cmake")
  include("/home/gulsum/galaxy/tools/racon/ClaraGenomicsAnalysis/cudamapper/cmake_install.cmake")
  include("/home/gulsum/galaxy/tools/racon/ClaraGenomicsAnalysis/cudaaligner/cmake_install.cmake")

endif()

