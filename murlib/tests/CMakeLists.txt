project(tests)

add_executable(tests "tests.cpp")

set_property(TARGET tests PROPERTY CXX_STANDARD 20)

target_link_libraries(tests murlib doctests)