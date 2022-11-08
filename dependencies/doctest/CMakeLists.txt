project (doctests LANGUAGES CXX)
add_library(doctests INTERFACE)

target_include_directories(doctests INTERFACE ${PROJECT_SOURCE_DIR}/include)