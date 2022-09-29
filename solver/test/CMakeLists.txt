find_package(doctest)


if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

add_executable(ponio_tests main.cpp)
target_link_libraries(ponio_tests ponio doctest::doctest)