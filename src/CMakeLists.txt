# Create the DQBF-BDD solver binary
add_executable(dqbdd
                dqbddexceptions.hpp
                dqbddformula.cpp
                dqbddformula.hpp
                hqspreinterface.cpp
                hqspreinterface.hpp
                main.cpp
                parser.hpp
                quantifiedvariablesmanipulator.cpp
                quantifiedvariablesmanipulator.hpp
                quantifiertree.cpp
                quantifiertree.hpp
                dqbddvariable.cpp
                dqbddvariable.hpp
                dqdimacsparser.cpp
                dqdimacsparser.hpp
                prenexcleansedqcirparser.cpp
                prenexcleansedqcirparser.hpp
                )
target_link_libraries(dqbdd PUBLIC cudd
                            PUBLIC hqspre
                            PUBLIC cxxopts)

# if we want, we can link dqbdd statically to standard library (except on mac, where it is not possible)
if (STATIC AND NOT APPLE)
    target_link_libraries(dqbdd PUBLIC "-static")
    set_target_properties(dqbdd PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()


target_compile_options(dqbdd PRIVATE -Wall -Wextra -Wpedantic)
