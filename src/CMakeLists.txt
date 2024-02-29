cmake_minimum_required(VERSION 3.24)
project(ar_gibbs)

set(CMAKE_CXX_STANDARD 17)

## External libraries
set(EIGEN_DIR "/users/daniel/Desktop/libraries/Eigen3")
set(BOOST_DIR "/users/daniel/Desktop/libraries/boost")
set(H_FILES "/users/daniel/Desktop/libraries/h_files")
set(STAN_DIR "/users/daniel/Desktop/libraries/math")
if(NOT EXISTS ${EIGEN_DIR})
    message(FATAL_ERROR "Please check that the set Eigen directory is valid!")
endif()
if(NOT EXISTS ${BOOST_DIR})
    message(FATAL_ERROR "Please check that the set boost directory is valid!")
endif()
if(NOT EXISTS ${STAN_DIR})
    message(FATAL_ERROR "Please check that the set stan directory is valid!")
endif()
if(NOT EXISTS ${H_FILES})
    message(FATAL_ERROR "Please check that the set sampler directory is valid!")
endif()

include_directories(${EIGEN_DIR})
include_directories(${BOOST_DIR})
include_directories(${STAN_DIR})
include_directories(${H_FILES})


find_package(Protobuf REQUIRED)

ADD_SUBDIRECTORY(proto)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${Protobuf_INCLUDE_DIRS})

add_executable(ar_gibbs
        main.cpp
        ${H_FILES}/eigenmvn.h
        matern.h
        coordinates/coordinates.h
        coordinates/coordinates.cpp
        calc_posterior/posterior.h
        calc_posterior/posterior.cpp
        ar_model/ar_class.cpp
        ar_model/ar_class.h
        debug_functions/debug.h
        protocpp/serialize.h
        protocpp/serialize.cpp)

TARGET_LINK_LIBRARIES(ar_gibbs proto ${Protobuf_LIBRARIES})