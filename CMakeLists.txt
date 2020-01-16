cmake_minimum_required(VERSION 3.5)
project(DQBF-BDD VERSION 0.1 LANGUAGES C CXX)

# Enforce C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# enable logging for HQSpre
add_definitions("-DLOGGING=true")

add_subdirectory("libs")
add_subdirectory("src")
