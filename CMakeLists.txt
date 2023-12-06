cmake_minimum_required(VERSION 3.27)
project(DiffInvertedPendulum)

#CHAI-3D
FIND_PACKAGE(CHAI3D REQUIRED CONFIG)
ADD_DEFINITIONS(${CHAI3D_DEFINITIONS})
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
    ADD_DEFINITIONS(-D_MACOSX -DMACOSX)
elseif (CMAKE_SYSTEM_NAME MATCHES Linux)
    ADD_DEFINITIONS(-D_LINUX -DLINUX)
endif () #NOTE: Windows is currently unsupported
LINK_DIRECTORIES(${CHAI3D_LIBRARY_DIRS})
INCLUDE_DIRECTORIES(${CHAI3D_INCLUDE_DIRS})#CHAI-3D

#EIGEN3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

#link to glfw
find_package(glfw3 REQUIRED)

set(CMAKE_CXX_STANDARD 17)

add_executable(DiffInvertedPendulum main.cpp)

target_include_directories(DiffInvertedPendulum PRIVATE ${GLFW3_INCLUDE_DIRS})

target_link_libraries(DiffInvertedPendulum
                PUBLIC
                Eigen3::Eigen
                glfw
                ${CHAI3D_LIBRARIES})
