project ("tests")

add_executable(${PROJECT_NAME} "tests.cpp")

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 20)