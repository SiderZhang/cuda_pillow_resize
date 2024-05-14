cmake_minimum_required(VERSION 3.16)

project(CUDAHelloWorld CUDA)

set(CMAKE_CUDA_STANDARD 11)

#link_directories(${CMAKE_SOURCE_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(/usr/include/opencv4 )
#set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
#add_library(CUDAHelloWorld SHARED main.cu)
add_executable(CUDAHelloWorld main.cu
        file.cu)
add_subdirectory(cv)

set_target_properties(CUDAHelloWorld
        PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

target_include_directories(CUDAHelloWorld PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} CVLoad)
target_link_libraries(${PROJECT_NAME} opencv_imgcodecs opencv_core opencv_imgproc)