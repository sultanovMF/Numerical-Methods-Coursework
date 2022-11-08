project (murlib LANGUAGES CXX)
add_library(murlib INTERFACE)

target_include_directories(murlib INTERFACE ${PROJECT_SOURCE_DIR}/include)