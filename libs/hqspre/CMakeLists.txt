# Create the HQSpre library

add_library(hqspre STATIC
              aux.hpp
              bool_vector.hpp
              bool_vector.ipp
              bool_vector.cpp
              clause.hpp
              clause.ipp
              clause.cpp
              exceptions.hpp
              formula.hpp
              formula.ipp
              formula.cpp
              formula_blocked_et_al.cpp
              formula_blocked_et_al.ipp
              formula_debug.cpp
              formula_dep_elim.cpp
              formula_dep_schemes.cpp
              formula_elim_set.cpp
              formula_equiv_contra.cpp
              formula_experimental.cpp
              formula_gate_interface.cpp
              formula_gates.cpp
              formula_inout.cpp
              formula_resolution.cpp
              formula_sat.cpp
              formula_substitute.cpp
              formula_subsumption.cpp
              formula_subsumption.ipp
              formula_unit_pure.cpp
              formula_univ_expand.cpp
              formula_univ_red.cpp
              formula_upla.cpp
              formula_vivification.cpp
              gate.hpp
              gate.cpp
              literal.hpp
              milp_solver.hpp
              milp_solver.cpp
              prefix.hpp
              prefix.ipp
              prefix.cpp
              process_limits.hpp
              propagator.hpp
              propagator.cpp
              sat_solver.hpp
              settings.hpp
              settings.cpp
              skolem.hpp
              timer.hpp
              timer.cpp
              varheap.hpp
           )
target_link_libraries(hqspre PUBLIC antom
                             PUBLIC picosat
                             PUBLIC easylogging
                             PUBLIC rt
                             PUBLIC glpk)
target_compile_options(hqspre PRIVATE "-Wconversion" "-Wuninitialized" "-Wpedantic")

if(HAVE_VALGRIND)
    target_compile_definitions(hqspre PUBLIC "-DHAVE_VALGRIND")
endif()
target_include_directories(hqspre PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")



# Create the HQSpre binary
# add_executable(preprocess main.cpp)
# if (STATIC)
#     target_link_libraries(preprocess PRIVATE "-static")
#     target_link_libraries(preprocess PRIVATE ltdl)
#     set_target_properties(preprocess PROPERTIES LINK_SEARCH_END_STATIC 1)
# endif()
# set_target_properties(preprocess PROPERTIES OUTPUT_NAME "hqspre")
# target_link_libraries(preprocess PUBLIC hqspre)
# target_link_libraries(preprocess PUBLIC easylogging)
# target_link_libraries(preprocess PUBLIC boost)
# target_compile_options(preprocess PRIVATE "-Wconversion")

# install(TARGETS preprocess
#         RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin
#         LIBRARY DESTINATION ${CMAKE_BINARY_DIR}/lib
#         ARCHIVE DESTINATION ${CMAKE_BINARY_DIR}/lib
#         INCLUDES DESTINATION ${CMAKE_BINARY_DIR}/include
#        )
