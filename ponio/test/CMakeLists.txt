find_package(doctest)

add_executable(ponio_tests main.cpp)
target_link_libraries(ponio_tests ponio doctest::doctest)
