cmake_minimum_required(VERSION 3.0)

project(bart_gadget)

cmake_policy(SET CMP0042 NEW)  # CMake 3.0 ``MACOSX_RPATH`` is enabled by default.

# ==============================================================================

include(CheckCXXCompilerFlag)
set(CXX_STANDARDS 17 14) # could maybe add 11 with some changes to the code

set(CXX_STANDARD_FOUND 0)
foreach(year ${CXX_STANDARDS})
  check_cxx_compiler_flag(-std=c++${year} COMPILER_SUPPORTS_CXX${year})
  if(COMPILER_SUPPORTS_CXX${year})
    set(CXX_STANDARD_FOUND 1)
    message(STATUS "Compiling using C++${year} standard")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${year}")
    break()
  endif(COMPILER_SUPPORTS_CXX${year})
endforeach()

if(NOT CXX_STANDARD_FOUND)
  message(FATAL_ERROR "The compiler does not support the most recent C++ standards. Please change compiler")
endif(NOT CXX_STANDARD_FOUND)

# ------------------------------------------------------------------------------

if (WIN32)
  add_definitions(-DWIN32 -D_WIN32 -D_WINDOWS)
  add_definitions(-DUNICODE -D_UNICODE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /W3")
  add_definitions(-D__BUILD_GADGETRON_bartgadget__)
endif (WIN32)

# ==============================================================================
# Bootstrap search for libraries 
# (We need to find cmake modules in Gadgetron)

find_path(GADGETRON_CMAKE_MODULES FindGadgetron.cmake
  HINTS
  $ENV{GADGETRON_HOME}/share/gadgetron/cmake
  /usr/local/gadgetron)

if(NOT GADGETRON_CMAKE_MODULES)
  message(FATAL_ERROR "GADGETRON_CMAKE_MODULES cannot be found. 
   Try to set GADGETRON_HOME environment variable.")
endif(NOT GADGETRON_CMAKE_MODULES)

set(CMAKE_MODULE_PATH ${GADGETRON_CMAKE_MODULES})

# ==============================================================================

find_package(ISMRMRD REQUIRED)
find_package(Gadgetron REQUIRED)

# We actually only use system and thread explicitly, but they require linking in date_time and chrono
if (WIN32)
  find_package(Boost COMPONENTS system thread date_time chrono program_options filesystem timer REQUIRED)
else()
  find_package(Boost COMPONENTS system thread program_options filesystem timer REQUIRED)
endif()

find_package(ACE REQUIRED)
find_package(Armadillo 4.600 REQUIRED)

set(CMAKE_INSTALL_PREFIX ${GADGETRON_HOME})

# ------------------------------------------------------------------------------

include_directories(
  ${GADGETRON_INCLUDE_DIR}
  ${GADGETRON_HOME}/gadgets/mri_core
  ${GADGETRON_HOME}/toolboxes/mri_core
  ${GADGETRON_HOME}/toolboxes/fft/cpu
  ${GADGETRON_HOME}/toolboxes/klt/cpu 
  ${GADGETRON_HOME}/toolboxes/image_io   
  ${GADGETRON_HOME}/toolboxes/core/cpu/math 
  ${Boost_INCLUDE_DIR}
  ${ISMRMRD_INCLUDE_DIR}
  ${ISMRMRD_SCHEMA_DIR}
  ${ISMRMRD_XSD_INCLUDE_DIR}
  ${ACE_INCLUDE_DIR})

link_directories(${GADGETRON_LIB_DIR})

# ------------------------------------------------------------------------------

add_library(gadgetron_baselbart SHARED 
  bartgadget.h
  bartgadget.cpp
  BART_Recon.xml 
  BART_Recon_cloud.xml
  BART_Recon_cloud_Standard.xml
)

set_target_properties(gadgetron_baselbart
  PROPERTIES
  VERSION ${GADGETRON_VERSION_STRING}
  SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_baselbart
  gadgetron_gadgetbase
  gadgetron_mricore 
  gadgetron_toolbox_mri_core 
  gadgetron_toolbox_log
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpufft
  gadgetron_toolbox_cpuklt  
  ${ISMRMRD_LIBRARIES}
  optimized ${ACE_LIBRARIES}
  debug ${ACE_DEBUG_LIBRARY}   
  ${Boost_LIBRARIES}
  )
  
if(ARMADILLO_FOUND)
  target_link_libraries(gadgetron_baselbart gadgetron_toolbox_cpucore_math )
endif()
  
# ------------------------------------------------------------------------------

set(GADGETRON_INSTALL_CMAKE_PATH share/gadgetron/cmake)
set(GADGETRON_INSTALL_CONFIG_PATH share/gadgetron/config)
set(GADGETRON_INSTALL_INCLUDE_PATH include/gadgetron)

install(FILES bartgadget.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH})

install(TARGETS gadgetron_baselbart DESTINATION lib)

install(FILES Sample_Grappa_Recon.sh Sample_Grappa_Recon_Standard.sh
  DESTINATION share/gadgetron/bart)
install(FILES BART_Recon.xml BART_Recon_cloud.xml BART_Recon_cloud_Standard.xml
  DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH})

# ------------------------------------------------------------------------------

# Add uninstall target
add_custom_target ( uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/uninstall.cmake)
