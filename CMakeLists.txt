cmake_minimum_required(VERSION 3.15)

project(ponio CXX)

set(CMAKE_CXX_STANDARD 20)

set(PONIO_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/solver/include)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(BUTCHER_METHODS "${CMAKE_CURRENT_SOURCE_DIR}/solver/include/solver/butcher_methods.hpp")
set(TEST_ORDER "${CMAKE_CURRENT_SOURCE_DIR}/solver/test/test_order.hxx")

set(Ndigit 36)

file(GLOB json_files
    "database/*.json")

add_custom_command(
    OUTPUT ${BUTCHER_METHODS} ${TEST_ORDER}
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/solver/code_generator.py ${json_files}
                -o ${BUTCHER_METHODS}
                -t ${TEST_ORDER}
                --Ndigit=${Ndigit}
    VERBATIM
)

add_library(ponio INTERFACE ${BUTCHER_METHODS})
target_include_directories(
        ponio
        INTERFACE $<BUILD_INTERFACE:${PONIO_INCLUDE_DIRS}>
                  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

option(BUILD_TESTS "build the tests" OFF)

if (BUILD_TESTS)
    add_subdirectory(solver/test)
endif()

add_subdirectory(solver/demos)
