cmake_minimum_required(VERSION 3.5)
project(HQSpre VERSION 1.4 LANGUAGES C CXX)

# Enforce C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANARD_REQUIRED ON)
#set(CMAKE_BUILD_TYPE "Release")
#set(CMAKE_BUILD_TYPE "Debug")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Enable warnings
#set(CMAKE_EXE_LINKER_FLAGS "-Wall -Wextra -Wpedantic -g ${CMAKE_EXE_LINKER_FLAGS}" )
#set(CMAKE_C_FLAGS          "-Wall -Wextra -Wpedantic -g ${CMAKE_C_FLAGS}")
#set(CMAKE_CXX_FLAGS        "-Wall -Wextra -Wpedantic -g ${CMAKE_CXX_FLAGS}")

#if(CMAKE_BUILD_TYPE STREQUAL "Release")
#    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")
#    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto-partition=none")
#        set(CMAKE_AR  "gcc-ar")
#        set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
#        set(CMAKE_CXX_ARCHIVE_FINISH true)
#    endif()
#    message("-- Enabled LTO.")
#endif()


option(STATIC:BOOL "Use static linking" OFF)
if (STATIC)
    set(CMAKE_EXE_LINKER_FLAGS "-static -static-libstdc++ ${CMAKE_EXE_LINKER_FLAGS}" )
#    set(CMAKE_C_FLAGS          "-static ${CMAKE_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS        "-static -static-libstdc++ ${CMAKE_CXX_FLAGS}")
    if (WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES}")
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    endif()
    message("-- Enabled static linking.")
endif()


# Enable the selection of different SAT solvers
option(SAT_SOLVER:String "Select the SAT solver to be used (Antom|Picosat)" "Antom")
if (SAT_SOLVER STREQUAL "Antom")
    add_definitions("-DSatSolver=Antom")
    message("-- Using Antom as a SAT solver")
elseif(SAT_SOLVER STREQUAL "Picosat")
    add_definitions("-DSatSolver=Picosat")
    message("-- Using Picosat as a SAT solver")
endif()


# Enable the computation of Skolem functions
option(SKOLEM:BOOL "Enable the computation of Skolem functions" OFF)
if (SKOLEM)
    add_definitions("-DSKOLEM")
    message("-- Enabled Skolem function computation")
else()
    message("-- Disabled Skolem function computation")
endif()

#option(LOGGING:BOOL "Enable logging" ON)
#if (NOT LOGGING)
#    add_definitions("-DELPP_DISABLE_LOGS -Wno-unused-parameter")
#    message("-- Logging disabled")
#endif()



if (NOT TARGET zlib)
    # Find libz (used by boost-iostreams)
    find_package(ZLIB REQUIRED)
    add_library(zlib IMPORTED INTERFACE)
    set_property(TARGET zlib PROPERTY INTERFACE_LINK_LIBRARIES  "${ZLIB_LIBRARIES}")
    set_property(TARGET zlib PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${ZLIB_INCLUDES}")
endif()


# Find librt (used for timeout handling in hqspre)
if (NOT TARGET rt)
    find_package(RT REQUIRED)
    add_library(rt IMPORTED INTERFACE)
    set_property(TARGET rt PROPERTY INTERFACE_LINK_LIBRARIES  "${RT_LIBRARIES}")
    set_property(TARGET rt PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${RT_INCLUDES}")
endif()


# Find Boost libraries
if (NOT TARGET boost)
    set(Boost_USE_STATIC_LIBS    OFF)
    set(Boost_USE_MULTITHREADED  OFF)
    set(Boost_USE_STATIC_RUNTIME OFF)
    find_package(Boost 1.40 REQUIRED COMPONENTS iostreams program_options)
    add_library(boost IMPORTED INTERFACE)
    set_property(TARGET boost PROPERTY INTERFACE_LINK_LIBRARIES  "${Boost_LIBRARIES}" zlib)
    set_property(TARGET boost PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}")
endif()


# Try to find GLPK (used as MILP solver)
if (NOT TARGET glpk)
    find_package(GLPK)
    if (GLPK_FOUND)
        add_library(glpk IMPORTED INTERFACE)
        set_property(TARGET glpk PROPERTY INTERFACE_COMPILE_OPTIONS     "-DHAVE_GLPK")
        set_property(TARGET glpk PROPERTY INTERFACE_LINK_LIBRARIES      "${GLPK_LIBRARIES}")
        set_property(TARGET glpk PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${GLPK_INCLUDE_DIRS}")
    endif()
endif()


# Enable OpenMP (used by the Antom SAT solver)
if (NOT TARGET openmp)
    find_package(OpenMP REQUIRED)
    add_library(openmp IMPORTED INTERFACE)
    set_property(TARGET openmp PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    set_property(TARGET openmp PROPERTY INTERFACE_LINK_LIBRARIES  ${OpenMP_CXX_FLAGS})
endif()

# Check if Doxygen is available; if yes, add the build target 'doc'
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)


# Check if valgrind is installed. This is used by the macro
# val_assert(...) to print a stack trace when the condition
# is violated and run through valgrind.
include(CheckIncludeFiles)
check_include_files("valgrind/valgrind.h" HAVE_VALGRIND)


add_subdirectory("libs/antom")
#add_subdirectory("libs/pacose")
add_subdirectory("libs/picosat")
add_subdirectory("libs/easylogging")
add_subdirectory("src")
