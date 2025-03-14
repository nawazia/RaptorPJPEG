cmake_minimum_required(VERSION 3.30)
project(jpeg)

set(CMAKE_CXX_STANDARD 20)

add_executable(jpeg
        decoder.cpp     # fragmentation decoder
        jpg.h
#        decode2.cpp    # og decoder
)
