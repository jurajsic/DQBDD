cmake_minimum_required(VERSION 3.10)
project(DQBF-BDD VERSION 0.1 LANGUAGES C CXX)

# Enforce C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANARD_REQUIRED ON)

# enable logging for HQSpre
add_definitions("-DLOGGING=true")

add_subdirectory("libs")
add_subdirectory("src")