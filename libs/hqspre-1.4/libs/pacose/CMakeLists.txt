# Set the include paths

# Create the antom library
add_library(pacose STATIC
    mtl/Alloc.h
    mtl/Alg.h
    mtl/XAlloc.h
    mtl/Vec.h
    mtl/Map.h
    mtl/Queue.h
    mtl/Sort.h
    mtl/IntTypes.h
    mtl/VecThreads.h
    mtl/Clone.h
    mtl/Heap.h
    simp/SimpSolver.h
    utils/System.h
    utils/ParseUtils.h
    utils/Options.h
    pacose/Pacose.h
    pacose/Encodings.h
    pacose/Softclause.h
    pacose/Settings.h
    pacose/MaxSATDimacs.h
    core/Solver.h
    core/Constants.h
    core/Dimacs.h
    core/SolverTypes.h
    core/SolverStats.h
    core/BoundedQueue.h
    simp/SimpSolver.cc
    utils/Options.cc
    utils/System.cc
    pacose/Settings.cc
    pacose/Encodings.cc
    pacose/Pacose.cc
    core/Solver.cc
)

#target_compile_definitions(pacose PUBLIC "-DINCREMENTAL")
target_include_directories(pacose PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
