cmake_minimum_required(VERSION 3.16)
project(CurlDownload)

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

include_directories(SYSTEM /usr/include/opencv4 )

add_library(CurlDownload STATIC
        SafeQueue.h
        SafeQueueTask.h
        ThreadPool.h
        Task.h
        download.cpp
        download.h
        cv.cpp
        file.cpp
        npy.h
)


FIND_PACKAGE(CURL)

if (CURL_FOUND)
    INCLUDE_DIRECTORIES(${CURL_INCLUDE_DIR})
    TARGET_LINK_LIBRARIES(CurlDownload ${CURL_LIBRARY})
ELSE (CURL_FOUND)
    MESSAGE(FATAL_ERROR "Failed to find CURL")
ENDIF(CURL_FOUND)