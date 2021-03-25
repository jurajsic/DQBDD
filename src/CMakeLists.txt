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
                prenexcleansedqcirparser.cpp
                prenexcleansedqcirparser.hpp
                )
target_link_libraries(DQBDD PUBLIC cudd
                            PUBLIC hqspre
                            PUBLIC cxxopts)

# if we want, we can link DQBDD statically to standard library (except on mac, where it is not possible)
if (STATIC AND NOT APPLE)
    target_link_libraries(DQBDD PUBLIC "-static")
    set_target_properties(DQBDD PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()


target_compile_options(DQBDD PRIVATE -Wall -Wextra -Wpedantic)
