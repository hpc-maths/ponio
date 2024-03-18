cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24:
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

# set version
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.release-please-manifest.json")
  file(READ "${CMAKE_CURRENT_SOURCE_DIR}/.release-please-manifest.json" PONIO_VERSION_JSON)
  string(JSON PONIO_VERSION GET "${PONIO_VERSION_JSON}" ".")
else()
  message(FATAL_ERROR "File ${CMAKE_CURRENT_SOURCE_DIR}/.release-please-manifest.json not found")
endif()

# add project_options v0.33.0
# https://github.com/aminya/project_options
include(FetchContent)
set(PROJECT_OPTIONS_VERSION "v0.33.0")
FetchContent_Declare(_project_options URL https://github.com/aminya/project_options/archive/refs/tags/${PROJECT_OPTIONS_VERSION}.zip)
FetchContent_MakeAvailable(_project_options)
include(${_project_options_SOURCE_DIR}/Index.cmake)

# install vcpkg dependencies: - should be called before defining project()
option(ENABLE_VCPKG "Use vcpkg to install dependencies" OFF)
option(ENABLE_CONAN_OPTION "Use Conan to install dependencies" OFF)

if(${ENABLE_VCPKG})
  run_vcpkg()
  set(VCPKG_BUILD_TYPE release)
endif()

if(${ENABLE_CONAN_OPTION})
  set(ENABLE_CONAN "ENABLE_CONAN")
endif()

# enable sanitizers and static analyzers when running the tests
option(CLANG_TIDY "Activate clang-tidy" OFF)
option(CPPCHECK "Activate cppcheck" OFF)
option(IWYU "Activate include-what-you-use" OFF)
option(SANITIZERS "Activate sanitizers" OFF)
option(ENABLE_COVERAGE "Activate coverage" OFF)

project(ponio VERSION ${PONIO_VERSION} LANGUAGES CXX)

SET(FEATURES)

if(${CLANG_TIDY})
  LIST(APPEND FEATURES ENABLE_CLANG_TIDY)
endif()
if(${CPPCHECK})
  LIST(APPEND FEATURES ENABLE_CPPCHECK)
endif()
if(${IWYU})
  LIST(APPEND FEATURES ENABLE_INCLUDE_WHAT_YOU_USE)
endif()
if(${SANITIZERS})
  LIST(APPEND FEATURES ENABLE_SANITIZER_ADDRESS)
  LIST(APPEND FEATURES ENABLE_SANITIZER_UNDEFINED_BEHAVIOR)
endif()
if(${COVERAGE})
  LIST(APPEND FEATURES ENABLE_COVERAGE)
endif()

message(STATUS "Available FEATURES: ${FEATURES}")
project_options(
  ${FEATURES}
  ENABLE_VS_ANALYSIS
  ${ENABLE_CONAN}
)

# includes
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ponio/include)

# generate butcher_methods.hpp
find_package(Python COMPONENTS Interpreter REQUIRED)
set(BUTCHER_METHODS "${INCLUDE_DIR}/ponio/runge_kutta/butcher_methods.hpp")

set(Ndigit 36)

file(GLOB json_files
    "database/*.json")

add_custom_command(
    OUTPUT ${BUTCHER_METHODS}
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ponio/code_generator.py ${json_files}
                -o ${BUTCHER_METHODS}
                --Ndigit=${Ndigit}
    COMMENT "Generating code..."
    VERBATIM
)

add_library(ponio INTERFACE ${BUTCHER_METHODS})
target_include_directories(
        ponio
        INTERFACE $<BUILD_INTERFACE:${INCLUDE_DIR}>
                  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ponio INTERFACE project_options project_warnings)
set_target_properties(ponio PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO)
target_compile_features(ponio INTERFACE cxx_std_20)

target_include_directories(ponio INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

# find dependencies:
set(DEPENDENCIES_CONFIGURED "")

if(${ENABLE_VCPKG})
  list(APPEND DEPENDENCIES_CONFIGURED hdf5)
endif()

foreach(DEPENDENCY ${DEPENDENCIES_CONFIGURED})
  find_package(${DEPENDENCY} CONFIG REQUIRED)
endforeach()

set(PONIO_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/ponio/include)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

option(BUILD_TESTS "build the tests" OFF)
option(BUILD_DEMOS "build the demos" OFF)
option(BUILD_SAMURAI_DEMOS "build demos working with samurai" OFF)

if (BUILD_TESTS)
    add_subdirectory(ponio/test)
endif()

if (BUILD_DEMOS)
  add_subdirectory(ponio/demos)
endif()

# package the project
package_project(
  TARGETS ponio project_options project_warnings
  INTERFACE_DEPENDENCIES_CONFIGURED ${DEPENDENCIES_CONFIGURED}
  INTERFACE_INCLUDES ${INCLUDE_DIR}
)
