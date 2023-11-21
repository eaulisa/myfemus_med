cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(myfemus_med
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

include(FetchContent)
FetchContent_Declare(fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG master
)
FetchContent_MakeAvailable(fmt)

set(MPIEXEC_EXECUTABLE "$ENV{PLAT_THIRD_PARTY_DIR}/openmpi/bin/mpiexec")
find_package(MPI REQUIRED COMPONENTS CXX)

set(hdf5_DIR "$ENV{PLAT_THIRD_PARTY_DIR}/hdf5/share/cmake/hdf5")
find_package(hdf5 REQUIRED NO_MODULE)

set(MEDCoupling_DIR "$ENV{PLAT_THIRD_PARTY_DIR}/medcoupling/cmake_files")
find_package(MEDCoupling REQUIRED PATHS $ENV{SPACK_ENV}/.spack-env/view/cmake_files)

include(FeatureSummary)

if (DEFINED MYFEMUS_ROOT)
  include(${MYFEMUS_ROOT}/UseFemus.CMake)
elseif(DEFINED ENV{MYFEMUS_ROOT})
  include($ENV{MYFEMUS_ROOT}/UseFemus.CMake)
elseif(EXISTS $ENV{PLAT_CODES_DIR}/MyFEMuS/build)
  include($ENV{PLAT_CODES_DIR}/MyFEMuS/build/UseFemus.CMake)
elseif(EXISTS $ENV{PLAT_CODES_DIR}/myfemus_build)
  include($ENV{PLAT_CODES_DIR}/myfemus_build/UseFemus.CMake)
else()
  message(FATAL_ERROR "MYFEMUS build directory not found. please set MYFEMUS_ROOT.")
endif()

add_executable(${PROJECT_NAME})

set(SRC_FILES
  main.cpp
  fem_ttu.cpp
  FemusMedCoupling.cpp
)
set(HEADER_FILES
  fem_ttu.hpp
  FemusMedCoupling.hpp
)

target_sources(${PROJECT_NAME}
  PRIVATE ${SRC_FILES}
  PUBLIC ${HEADER_FILES}
)
target_include_directories(${PROJECT_NAME} PUBLIC ${FEMUS_INCLUDES})
target_link_libraries(${PROJECT_NAME}
  PUBLIC
    MPI::MPI_CXX
    ${FEMUS_LIBRARIES}
    interpkernel
    medcouplingcpp
    medcouplingremapper
    medloader
    hdf5::hdf5-shared
    fmt::fmt
)

feature_summary(WHAT ALL)

find_program(CLANGFORMAT clang-format)
if (CLANGFORMAT)
  add_custom_target(format
    COMMAND ${CLANGFORMAT}
    --verbose -i
    ${SRC_FILES} ${HEADER_FILES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
else()
  add_custom_target(format
    COMMAND echo "clang-format could not be found"
  )
endif()

