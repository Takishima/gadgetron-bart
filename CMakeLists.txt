cmake_minimum_required(VERSION 3.0)

project(bart_gadget CXX)

cmake_policy(SET CMP0042 NEW)  # CMake 3.0 ``MACOSX_RPATH`` is enabled by default.

set(BART_GIT_REPOSITORY "git@github.com:Takishima/bart.git")
set(BART_GIT_TAG bart-single-library)

# ==============================================================================

option(USE_CUDA "Build using support for CUDA (should match the library BART(main) compile option)" OFF)
if(USE_CUDA)
  find_package(CUDA)
endif(USE_CUDA)

option(USE_OPENMP "Build using support for OpenMP" OFF)
if(USE_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif(OPENMP_FOUND)
endif(USE_OPENMP)

# ==============================================================================

macro(setup_default_bart_options)
  set(BARTMAIN_AUTOCOPY      OFF CACHE STRING "" ${ARGN})
  set(BUILD_TESTING          OFF CACHE STRING "" ${ARGN})
  set(CREATE_PYTHON_MODULE   OFF CACHE STRING "" ${ARGN})
  set(DISABLE_PNG            ON  CACHE STRING "" ${ARGN})
  set(ENABLE_LONGJUMP        ON  CACHE STRING "" ${ARGN})
  set(ENABLE_MEM_CFL         ON  CACHE STRING "" ${ARGN})
  set(FORCE_BUILTIN_COMMANDS ON  CACHE STRING "" ${ARGN})
  set(GENERATE_DOC           OFF CACHE STRING "" ${ARGN})
  set(USE_FFTWTHREADS        ON  CACHE STRING "" ${ARGN})
  set(USE_FPIC               ON  CACHE STRING "" ${ARGN})
  set(USE_LOCAL_FFTW         OFF CACHE STRING "" ${ARGN})
  set(USE_LOG_BACKEND        OFF CACHE STRING "" ${ARGN})
  set(USE_MEMONLY_CFL        ON  CACHE STRING "" ${ARGN})

  set(LINALG_VENDOR ${BART_LINALG_VENDOR} CACHE STRING "" ${ARGN})
endmacro()

# ------------------------------------------------------------------------------

if(NOT DEFINED BART_LINALG_VENDOR)
  set(BART_LINALG_VENDOR "NOT-DEFINED" CACHE FILEPATH "Corresponds to LINALG_VENDOR for BART compilation")
endif(NOT DEFINED BART_LINALG_VENDOR)
if(NOT DEFINED BART_LIBRARIES)
  set(BART_LIBRARIES "NOT-DEFINED" CACHE STRING "List of libraries required to link to BART(main)")
endif(NOT DEFINED BART_LIBRARIES)

option(BARTMAIN_DOWNLOAD_AND_BUILD "Automatically download and build the bartmain library as part of the compilation process" OFF)
if(BARTMAIN_DOWNLOAD_AND_BUILD)
  if(BART_LINALG_VENDOR STREQUAL "NOT-DEFINED")
    message(FATAL_ERROR "You must specify the BART_LINALG_VENDOR variable in order to properly build BART (will be passed on to the BART build process as LINALG_VENDOR)")
  endif(BART_LINALG_VENDOR STREQUAL "NOT-DEFINED")

  set(DOWNLOAD_BUILD_GOT_FOLDER OFF)
  mark_as_advanced(DOWNLOAD_BUILD_GOT_FOLDER)

  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.7)
    # If CMake >= 3.7 we can use the FetchContent module
    set(DOWNLOAD_BUILD_GOT_FOLDER ON)
    include(FetchContent)
    FetchContent_Declare(bartmain_git
      GIT_REPOSITORY ${BART_GIT_REPOSITORY}
      GIT_TAG        ${BART_GIT_TAG})

    FetchContent_GetProperties(bartmain_git)
    if(NOT bartmain_git_POPULATED)
      setup_default_bart_options(FORCE)
      FetchContent_Populate(bartmain_git)
      add_subdirectory(${bartmain_git_SOURCE_DIR} ${bartmain_git_BINARY_DIR}
	EXCLUDE_FROM_ALL)
    endif(NOT bartmain_git_POPULATED)
  else(CMAKE_VERSION VERSION_GREATER_EQUAL 3.7)
    # If CMake < 3.7 we cannot use the FetchContent module, but we should be
    # able to use the ExternalProject module
    # We therefore force the user to download and configure BART instead of
    # starting the building process when he/she runs "make" the first time
    # A rebuilding of the cache after downloading the data is automatically
    # started in order to properly setup the downloaded BART as a
    # sub-directory
    if(EXISTS ${CMAKE_SOURCE_DIR}/bart/CMakeLists.txt)
      set(DOWNLOAD_BUILD_GOT_FOLDER ON)
      setup_default_bart_options()
      add_subdirectory(${CMAKE_SOURCE_DIR}/bart ${CMAKE_BINARY_DIR}/bart
	EXCLUDE_FROM_ALL)
    else(EXISTS ${CMAKE_SOURCE_DIR}/bart/CMakeLists.txt)
      setup_default_bart_options(FORCE)
      message(WARNING
	"                   ********** IMPORTANT **********\n"
	"Since you're using CMake < 3.7, I cannot use the FetchContent module\n"
	"After this configuration ends, running \"make all\" will automatically\n"
	"download and (default) configure BART for compilation into the BartGadget\n"
	"A second call to \"make all\" will then be required in order to actually\n"
	"build the project.\n"
	"If you're ok with the defaults, simply run \"make make_2nd_stage\"\n"
	"                   ********** IMPORTANT **********")

      include(ExternalProject)
      ExternalProject_Add(bartmain_git
	GIT_REPOSITORY ${BART_GIT_REPOSITORY}
	GIT_TAG        ${BART_GIT_TAG}

	CONFIGURE_COMMAND ""
	UPDATE_COMMAND ""
	PATCH_COMMAND ""
	BUILD_COMMAND ""
	INSTALL_COMMAND ""

	SOURCE_DIR "${CMAKE_SOURCE_DIR}/bart"
	BINARY_DIR "${CMAKE_BINARY_DIR}/bart"
	)
    endif(EXISTS ${CMAKE_SOURCE_DIR}/bart/CMakeLists.txt)
  endif(CMAKE_VERSION VERSION_GREATER_EQUAL 3.7)

else(BARTMAIN_DOWNLOAD_AND_BUILD)
  set(BART_LINALG_VENDOR "" CHACHE FILEPATH "Corresponds to LINALG_VENDOR for BART compilation")
  find_file(BART_LIBMAIN
    NAMES libbartmain.a libbartmain.so
    PATH /usr/local/bart /opt/bart /opt/local/ /usr/local /usr
    PATH_SUFFIXES lib lib/static lib64)

  if(NOT BART_LIBMAIN)
    message(FATAL_ERROR "Unable to locate the BART(main) library")
  else(NOT BART_LIBMAIN)
    message(STATUS "Using BART(main) library: ${BART_LIBMAIN}")
  endif(NOT BART_LIBMAIN)
endif(BARTMAIN_DOWNLOAD_AND_BUILD)

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

list(APPEND CMAKE_MODULE_PATH ${GADGETRON_CMAKE_MODULES})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

# ==============================================================================

find_package(FFTW REQUIRED)
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

# ==============================================================================

set(GADGET_FILES
  bartgadget.h
  bartgadget.cpp
  BART_Recon.xml
  BART_Recon_cloud.xml
  BART_Recon_cloud_Standard.xml
)

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

# ==============================================================================

if(BARTMAIN_DOWNLOAD_AND_BUILD)
  if(DOWNLOAD_BUILD_GOT_FOLDER)
    add_library(gadgetron_baselbart SHARED
      ${GADGET_FILES}
      $<TARGET_OBJECTS:bartsupport_objs>
      $<TARGET_OBJECTS:bartmain_objs>
      ${BARTSUPPORT_PYTHON_SEPARATE_COMPILATION_FILES}
      ${BARTMAIN_PYTHON_SEPARATE_COMPILATION_FILES}
      ${bartsupport_cuda_objs})
    get_property(BARTMAIN_LIBS TARGET bartmain PROPERTY LINK_LIBRARIES)
    get_property(BARTSUPPORT_LIBS TARGET bartsupport PROPERTY LINK_LIBRARIES)
    set(BART_LIBRARIES "${BARTSUPPORT_LIBS};${BARTMAIN_LIBS};${LINALG_LIBRARIES}" CACHE STRING "" FORCE)
    target_link_libraries(gadgetron_baselbart ${BART_LIBRARIES})
  else(DOWNLOAD_BUILD_GOT_FOLDER)
    # In this case, the BART source folder does not exist yet, so populate it
    # the next time the build command is issued.
    # Also automatically reconfigure CMake after
    add_custom_target(make_1st_stage
      COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target bartmain_git --config $<CONFIG>)
    if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
      add_custom_target(reconfigure_for_bart
	COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target rebuild_cache --config $<CONFIG>)
    elseif(CMAKE_GENERATOR MATCHES "Visual Studio.*")
      add_custom_target(reconfigure_for_bart
	COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ZERO_CHECK --config $<CONFIG>)
    endif(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    add_custom_target(print_info ALL
      COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "Downloading and setup of BART done in ${CMAKE_SOURCE_DIR}/bart. Please check CMake cache properties and run make again")
    add_dependencies(reconfigure_for_bart make_1st_stage)
    add_dependencies(print_info reconfigure_for_bart make_1st_stage)

    # For convenience
    add_custom_target(make_2nd_stage
      COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --config $<CONFIG>)
    add_dependencies(make_2nd_stage reconfigure_for_bart make_1st_stage)
  endif(DOWNLOAD_BUILD_GOT_FOLDER)
else(BARTMAIN_DOWNLOAD_AND_BUILD)
  add_library(gadgetron_baselbart SHARED ${GADGET_FILES})
  target_link_libraries(main ${BART_LIBMAIN})
  target_link_libraries(main ${BART_LIBRARIES})
endif(BARTMAIN_DOWNLOAD_AND_BUILD)

# ------------------------------------------------------------------------------

if(TARGET gadgetron_baselbart)
  set_target_properties(gadgetron_baselbart
    PROPERTIES
    VERSION   ${GADGETRON_VERSION_STRING}
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

  # if(USE_CUDA)
  #   CUDA_ADD_CUFFT_TO_TARGET(gadgetron_baselbart)
  #   CUDA_ADD_CUBLAS_TO_TARGET(gadgetron_baselbart)
  #   target_link_libraries(gadgetron_baselbart ${CUDA_LIBRARIES})
  # endif(USE_CUDA)

  if(ARMADILLO_FOUND)
    target_link_libraries(gadgetron_baselbart gadgetron_toolbox_cpucore_math )
  endif(ARMADILLO_FOUND)
  
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
  add_custom_target (uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/uninstall.cmake)
endif(TARGET gadgetron_baselbart)

