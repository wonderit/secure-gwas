cmake_minimum_required(VERSION 3.13)  # CMake version check
project(mpc_multiply)               # Create project "simple_example"
set(CMAKE_CXX_STANDARD 11)            # Enable c++14 standard
#set(INCPATHS -I/Users/wonsuk/projects/local/include)
#set(LDLIBS -lssl -lcrypto -lntl -lgmp -lm -ldl)
#set(LDPATH -L/usr/local/lib)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -O1 -march=native -std=c++11 -pthread")
message(STATUS "start running cmake...")

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES
        crypto.cpp
        aesstream.cpp
        param.cpp
        DataSharingClient.cpp
        mpc_multiply.cpp
        )

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(mpc_multiply ${SOURCE_FILES})

include_directories(/usr/local/include/)