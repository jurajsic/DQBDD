# Set the include paths

# Create the antom library
add_library(antom STATIC
  Core/allocator.h
  Core/boundedQueue.h
  Core/clause.h
  Core/control.h
  Core/core.h
  Core/helper.h
  Core/reason.h
  Core/solverstate.cpp
  Core/statistics.h
  Core/varheap.h
  Core/watcher.h
  prepro/bce.h
  prepro/bce.cpp
  prepro/bva.h
  prepro/bva.cpp
  prepro/hte.h
  prepro/hte.cpp
  prepro/modelrebuilder.h
  prepro/modelrebuilder.cpp
  prepro/preprocessor.h
  prepro/preprocessor.cpp
  prepro/upla.h
  prepro/upla.cpp
  prepro/varcandidate.h
  solver/antombase.h
  solver/antombase.cpp
  solver/antom.h
  solver/antom.cpp
  solver/bucket.h
  solver/bucket.cpp
  solver/cascade.h
  solver/cascade.cpp
  solver/debug.cpp
  solver/multiplecascade.h
  solver/multiplecascade.cpp
  solver/parser.h
  solver/parser.cpp
  solver/settings.cpp
  solver/settings.h
  solver/softclause.h
  solver/softclausenodes.h
  solver/solverproxy.h
  solver/solverproxy.cpp
  solver/sorter.h
  solver/sorter.cpp
  solver/timemeasurement.h
  solver/timevariables.h
  solver/totalizerencodetree.h
  )

target_include_directories(antom PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(antom PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Core")
target_include_directories(antom PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/prepro")
target_include_directories(antom PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/solver")
#target_link_libraries(antom PUBLIC openmp)
