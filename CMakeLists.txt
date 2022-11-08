cmake_minimum_required (VERSION  3.14)

project ("kushtau")

set(BUILD_SHARED_LIBS OFF)

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
endif(CCACHE_FOUND)


add_subdirectory(source)
add_subdirectory(dependencies)
add_subdirectory(tests)
add_subdirectory(murlib)

find_package(Python COMPONENTS Interpreter)