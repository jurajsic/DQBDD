# Create the easylogging++ library
add_library(easylogging STATIC easylogging++.hpp easylogging++.cpp)
target_include_directories(easylogging PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(easylogging PUBLIC "-DELPP_NO_DEFAULT_LOG_FILE -Wno-unused-variable")
