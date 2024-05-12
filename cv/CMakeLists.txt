cmake_minimum_required(VERSION 3.16)

project(CVLoad)

#set(CMAKE_CUDA_STANDARD 11)

include_directories(SYSTEM /usr/include/opencv4 )
#add_definitions(-D _GLIBCXX_USE_CXX11_ABI=0)
add_library(CVLoad STATIC
        cv.cpp)

set_target_properties(CUDAHelloWorld
        PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

#target_link_libraries(${PROJECT_NAME} opencv_imgcodecs opencv_core opencv_imgproc)
#target_link_libraries(${PROJECT_NAME} `pkg-config --libs opencv4`)