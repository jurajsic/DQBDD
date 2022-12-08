cmake_minimum_required(VERSION 3.5)
project(dqbdd VERSION 1.3 LANGUAGES C CXX)

# Enforce C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckIncludeFiles)
include(CheckFunctionExists)
include(CheckLibraryExists)
include(CheckTypeSize)
include(TestBigEndian)

# Static linking is enabled by default...
option(STATIC "Use static linking" ON)
# ...except for mac, where static linking does not work
if (STATIC AND NOT APPLE)
    if (WIN32)
        # this is only for completness, compiling DQBDD for windows will almost definitely not work
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES}")
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    endif()
    message("-- Enabled static linking.")
endif()

configure_file(src/version.hpp.in src/version.hpp)

########################################################################################
################################### Stuff for HQSpre ###################################
########################################################################################

# enable logging for HQSpre
add_definitions("-DLOGGING=true")

# Check if valgrind is installed. This is used by the macro
# val_assert(...) to print a stack trace when the condition
# is violated and run through valgrind. (used in HQSpre)
check_include_files("valgrind/valgrind.h" HAVE_VALGRIND)
if(HAVE_VALGRIND)
  add_definitions("-DHAVE_VALGRIND")
endif()

########################################################################################
########################################################################################
########################################################################################

add_subdirectory("libs")
add_subdirectory("src")
