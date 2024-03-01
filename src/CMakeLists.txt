cmake_minimum_required(VERSION 3.24)
project(starr)

set(CMAKE_CXX_STANDARD 17)

## External libraries
#set(EIGEN_DIR "/users/daniel/Desktop/libraries/Eigen3")
#set(BOOST_DIR "/users/daniel/Desktop/libraries/boost")
set(H_FILES "${CMAKE_CURRENT_LIST_DIR}/h_files")
#if(NOT EXISTS ${EIGEN_DIR})
#    message(FATAL_ERROR "Please check that the set Eigen directory is valid!")
#endif()
#if(NOT EXISTS ${BOOST_DIR})
 #   message(FATAL_ERROR "Please check that the set boost directory is valid!")
#endif()

if(NOT EXISTS ${H_FILES})
    message(FATAL_ERROR "Please check that the set sampler directory is valid!")
endif()

#include_directories(${EIGEN_DIR})
##include_directories(${BOOST_DIR})
include_directories(${H_FILES})

#find_package(Protobuf REQUIRED)

##Eigen
find_package(Eigen3 3.1.2 QUIET)
if(NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found. Fetching and installing...")
    include(FetchContent)
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0 
    )
    FetchContent_MakeAvailable(eigen)
    message(STATUS "Eigen3 installed.")
endif()

find_package(Boost QUIET)
if(NOT Boost_FOUND)
    message(STATUS "Boost not found. Fetching and installing...")
    include(FetchContent)
    FetchContent_Declare(
        boost
        URL https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.tar.bz2
        URL_HASH SHA256=a6e1ab9b0860e6a2881dd7b21fe9f737a095e5f33a3a874afc6a345228597ee6

    )
    FetchContent_MakeAvailable(boost)
    message(STATUS "Boost installed.")
    set(BOOST_DIR "${CMAKE_CURRENT_LIST_DIR}/build/_deps/boost-src")
    if(NOT EXISTS ${BOOST_DIR})
        message(FATAL_ERROR "Please check that the set boost directory is valid!")
    endif()
endif()

include_directories(${BOOST_DIR})

find_package(Protobuf QUIET)
if(NOT Protobuf_FOUND)
    message(STATUS "Boost not found. Fetching and installing...")
    include(FetchContent)
    FetchContent_Declare(
        protobuf
        URL https://github.com/protocolbuffers/protobuf/releases/download/v25.3/protobuf-25.3.tar.gz
    )
    FetchContent_MakeAvailable(protobuf)
    message(STATUS "protobuf installed.")
endif()

ADD_SUBDIRECTORY(proto)


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${Protobuf_INCLUDE_DIRS})

add_executable(main
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

TARGET_LINK_LIBRARIES(main proto protobuf::libprotobuf Eigen3::Eigen)
