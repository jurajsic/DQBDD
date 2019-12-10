# Create the PicoSat SAT solver library
add_library(picosat STATIC picosat.c version.c)
target_compile_definitions(picosat PRIVATE "-DTRACE")
target_include_directories(picosat PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
