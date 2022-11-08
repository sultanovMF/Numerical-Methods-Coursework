add_executable(${PROJECT_NAME} "main.cpp")

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 20)
target_link_libraries(${PROJECT_NAME} murlib indicators)
if (MSVC)
    add_compile_options(/Ox)
    add_compile_options(/QxHost)
endif()

