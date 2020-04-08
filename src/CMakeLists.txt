# Create the DQBF-BDD solver binary
add_executable(DQBDD
                DQBDDexceptions.hpp
                DQBDDformula.cpp
                DQBDDformula.hpp
                HQSpreinterface.cpp
                HQSpreinterface.hpp
                main.cpp
                parser.hpp
                quantifiedvariablesmanipulator.cpp
                quantifiedvariablesmanipulator.hpp
                quantifiertree.cpp
                quantifiertree.hpp
                DQBDDvariable.cpp
                DQBDDvariable.hpp
                DQDIMACSparser.cpp
                DQDIMACSparser.hpp
                )
target_link_libraries(DQBDD PUBLIC cudd
                            PUBLIC hqspre
                            PUBLIC cxxopts)

# old way to link to cudd
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../cudd/)
#target_link_libraries(solver ${CMAKE_CURRENT_SOURCE_DIR}/../cudd/libs/libobj.a)
#target_link_libraries(solver ${CMAKE_CURRENT_SOURCE_DIR}/../cudd/libs/libcudd.a)

# statically link to stuff
target_link_libraries(DQBDD PRIVATE -static)
target_link_libraries(DQBDD PRIVATE ltdl)
set_target_properties(DQBDD PROPERTIES LINK_SEARCH_END_STATIC 1)

target_compile_options(DQBDD PRIVATE -Wall -Wextra -Wpedantic)

install(TARGETS DQBDD
        RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin
        LIBRARY DESTINATION ${CMAKE_BINARY_DIR}/lib
        ARCHIVE DESTINATION ${CMAKE_BINARY_DIR}/lib
        INCLUDES DESTINATION ${CMAKE_BINARY_DIR}/include
       )

