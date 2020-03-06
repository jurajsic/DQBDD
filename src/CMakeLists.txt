# Create the DQBF-BDD solver binary
add_executable(solver
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
                #simplesolver.cpp
                #simplesolver.hpp
                #solver.cpp
                #solver.hpp
                #treesolver.cpp
                #treesolver.hpp
                DQBDDvariable.cpp
                DQBDDvariable.hpp
                )
target_link_libraries(solver hqspre)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../cudd/)
target_link_libraries(solver ${CMAKE_CURRENT_SOURCE_DIR}/../cudd/libs/libobj.a)
target_link_libraries(solver ${CMAKE_CURRENT_SOURCE_DIR}/../cudd/libs/libcudd.a)
# statically link to standard library
target_link_libraries(solver -static-libgcc -static-libstdc++)

install(TARGETS solver
        RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin
        LIBRARY DESTINATION ${CMAKE_BINARY_DIR}/lib
        ARCHIVE DESTINATION ${CMAKE_BINARY_DIR}/lib
        INCLUDES DESTINATION ${CMAKE_BINARY_DIR}/include
       )
