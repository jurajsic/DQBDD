cmake_minimum_required(VERSION 3.5)
project(dqbdd VERSION 1.2 LANGUAGES C CXX)

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
