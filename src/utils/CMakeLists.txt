add_library(utils stl_to_string.cpp clock.cpp config.cpp logging1.cpp)
target_link_libraries(utils ${Boost_PROGRAM_OPTIONS_LIBRARY})
