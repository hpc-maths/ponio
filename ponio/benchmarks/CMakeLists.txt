set(benchmarks heat_bench)

find_package(CLI11 REQUIRED)

message(STATUS "list of benchmarks:")
foreach(bench IN LISTS benchmarks)
    message(STATUS "+ ${bench}")
    add_executable(${bench} ${bench}.cpp)

    # -g gives profilers source lines and inline frames without changing the generated code,
    # frame pointers give them cheap call stacks
    if(NOT MSVC)
        target_compile_options(${bench} PRIVATE -g -fno-omit-frame-pointer)
    endif()

    target_link_libraries(${bench} ponio CLI11::CLI11)
endforeach()
