set(demos arenstorf brownian brusselator brusselator_dirk exp heat pendulum lorenz lotka_volterra)

foreach(demo IN LISTS demos)
    message(STATUS "demo=${demo}")
    add_executable(${demo} ${demo}.cpp)
    target_link_libraries(${demo} ponio)
endforeach()

foreach(demo IN LISTS demos)
    add_custom_target(visu_${demo} python ${CMAKE_SOURCE_DIR}/solver/demos/${demo}.py)
    add_custom_target(${demo}_visu make visu_${demo})
endforeach()
