cmake_minimum_required(VERSION 3.5)
project(DQBDD VERSION 0.2 LANGUAGES C CXX)

# Enforce C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# use only static libraries
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

########################################################################################
################################### Stuff for HQSpre ###################################
########################################################################################

include(CheckIncludeFiles)
include(CheckFunctionExists)
include(CheckLibraryExists)
include(CheckTypeSize)
include(TestBigEndian)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")


# enable logging for HQSpre
add_definitions("-DLOGGING=true")

# Find librt (used for timeout handling in HQSpre)
find_package(RT REQUIRED)
add_library(rt IMPORTED INTERFACE)
set_property(TARGET rt PROPERTY INTERFACE_LINK_LIBRARIES  "${RT_LIBRARIES}")
set_property(TARGET rt PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${RT_INCLUDES}")

# Try to find GLPK (used as MILP solver in HQSpre)
find_package(GLPK)
if (GLPK_FOUND)
    add_library(glpk IMPORTED INTERFACE)
    set_property(TARGET glpk PROPERTY INTERFACE_COMPILE_OPTIONS     "-DHAVE_GLPK")
    set_property(TARGET glpk PROPERTY INTERFACE_LINK_LIBRARIES      "${GLPK_LIBRARIES};${LTDL_LIBRARIES}")
    set_property(TARGET glpk PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${GLPK_INCLUDE_DIRS}")
endif()

find_package(LTDL REQUIRED)
add_library(ltdl IMPORTED INTERFACE)
set_property(TARGET ltdl PROPERTY INTERFACE_LINK_LIBRARIES  "${LTDL_LIBRARIES}")
set_property(TARGET ltdl PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${LTDL_INCLUDES}")

# Enable OpenMP (used by the Antom SAT solver which is used in HQSpre)
find_package(OpenMP REQUIRED)
add_library(openmp IMPORTED INTERFACE)
set_property(TARGET openmp PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
set_property(TARGET openmp PROPERTY INTERFACE_LINK_LIBRARIES  ${OpenMP_CXX_FLAGS})

# Check if valgrind is installed. This is used by the macro
# val_assert(...) to print a stack trace when the condition
# is violated and run through valgrind.
check_include_files("valgrind/valgrind.h" HAVE_VALGRIND)
if(HAVE_VALGRIND)
  add_definitions("-DHAVE_VALGRIND")
endif()

########################################################################################
########################################################################################
########################################################################################

########################################################################################
################################# Stuff for CUDD #######################################
########################################################################################

check_include_files("assert.h" HAVE_ASSERT_H)
check_include_files("float.h" HAVE_FLOAT_H)
check_include_files("dlfcn.h" HAVE_DLFCN_H)
check_include_files("inttypes.h" HAVE_INTTYPES_H)
check_include_files("limits.h" HAVE_LIMITS_H)
check_include_files("math.h" HAVE_MATH_H)
check_include_files("memory.h" HAVE_MEMORY_H)
check_include_files("stddef.h" HAVE_STDDEF_H)
check_include_files("stdint.h" HAVE_STDINT_H)
check_include_files("stdlib.h" HAVE_STDLIB_H)
check_include_files("strings.h" HAVE_STRINGS_H)
check_include_files("string.h" HAVE_STRING_H)
check_include_files("unistd.h" HAVE_UNISTD_H)
check_include_files("sys/resource.h" HAVE_SYS_RESOURCE_H)
check_include_files("sys/stat.h" HAVE_SYS_STAT_H)
check_include_files("sys/times.h" HAVE_SYS_TIMES_H)
check_include_files("sys/time.h" HAVE_SYS_TIME_H)
check_include_files("sys/types.h" HAVE_SYS_TYPES_H)
check_include_files("sys/wait.h" HAVE_SYS_WAIT_H)
check_function_exists("gethostname" HAVE_GETHOSTNAME)
check_function_exists("getrlimit" HAVE_GETRLIMIT)
check_function_exists("getrusage" HAVE_GETRUSAGE)
check_library_exists("m" "pow" "math.h" HAVE_POW)
check_library_exists("m" "powl" "math.h" HAVE_POWL)
check_library_exists("m" "sqrt" "math.h" HAVE_SQRT)
check_function_exists("strchr" HAVE_STRCHR)
check_function_exists("strstr" HAVE_STRSTR)
check_function_exists("sysconf" HAVE_SYSCONF)
check_type_size("int" SIZEOF_INT)
check_type_size("long" SIZEOF_LONG)
check_type_size("long double" SIZEOF_LONG_DOUBLE)
check_type_size("void*" SIZEOF_VOID_P)
set(CMAKE_EXTRA_INCLUDE_FILES "stddef.h")
check_type_size("ptrdiff_t" PTRDIFF_T)
set(CMAKE_EXTRA_INCLUDE_FILES)
test_big_endian(WORDS_BIGENDIAN)
set(CUDD_VERSION "3.0.0")

configure_file(libs/cudd-3.0.0/cudd_config.h.in libs/cudd-3.0.0/cudd_config.h @ONLY)

########################################################################################
########################################################################################
########################################################################################

add_subdirectory("libs")
add_subdirectory("src")
