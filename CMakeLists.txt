cmake_minimum_required(VERSION 3.16)
project(CUDAHelloWorld CUDA)

set(CMAKE_CUDA_STANDARD 14)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

link_directories(${CMAKE_SOURCE_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR}/include)

set(SOURCE_FILES pil.cu)
add_executable(CUDAHelloWorld main.cu
        pil.cu
        ipl.h
        Imaging.h
        imageOps.cpp
        jpeg.cu)

set_target_properties(CUDAHelloWorld
        PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(${PROJECT_NAME} libjpeg.a)
