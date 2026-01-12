set(examples arenstorf brusselator exp exp_splitting heat heat_rock pendulum lorenz curtiss_hirschfelder curtiss_hirschfelder_exprk nagumo)

set(examples_eigen brusselator_dirk lorenz_tuto lorenz_all)
set(examples_cli11 brownian lotka_volterra)
set(examples_samurai heat_samurai nagumo_pirock belousov_zhabotinsky_pirock bz_2d_pirock)

if(BUILD_ALL_EXAMPLES)
    set(BUILD_EIGEN_EXAMPLES ON)
    set(BUILD_CLI11_EXAMPLES ON)
    set(BUILD_SAMURAI_EXAMPLES ON)
endif()

if(BUILD_EIGEN_EXAMPLES)
    list(APPEND examples ${examples_eigen})
    find_package (Eigen3 3.4 REQUIRED NO_MODULE)
endif()

if(BUILD_CLI11_EXAMPLES)
    list(APPEND examples ${examples_cli11})
    find_package(CLI11)
endif()

if(BUILD_SAMURAI_EXAMPLES)
    list(APPEND examples ${examples_samurai})
    set(SAMURAI_WITH_PETSC ON)
    find_package(samurai CONFIG REQUIRED)
endif()

message(STATUS "list of examples:")

foreach(demo IN LISTS demos)
    message(STATUS "+ ${demo}")
    add_executable(${demo} ${demo}.cpp)

    if(${demo} IN_LIST demos_eigen)
        message(STATUS "    - eigen")
        target_link_libraries(${demo} Eigen3::Eigen)
    endif()

    if(${demo} IN_LIST demos_cli11)
        message(STATUS "    - cli11")
        target_link_libraries(${demo} CLI11::CLI11)
    endif()

    if(${demo} IN_LIST examples_samurai)
        message(STATUS "    - samurai")
        target_link_libraries(${demo} samurai::samurai)
    endif()

    target_link_libraries(${demo} ponio)
endforeach()

set(visu_list "")
set(save_visu_list "")
foreach(demo IN LISTS demos)
    add_custom_target(visu_${demo} python ${CMAKE_SOURCE_DIR}/ponio/examples/${demo}.py)
    add_custom_target(${demo}_visu make visu_${demo})
    list(APPEND visu_list visu_${demo})

    add_custom_target(save_visu_${demo} python ${CMAKE_SOURCE_DIR}/ponio/examples/${demo}.py --only-save)
    add_custom_target(${demo}_save_visu make save_visu_${demo})
    list(APPEND save_visu_list save_visu_${demo})
endforeach()

add_custom_target(visu DEPENDS ${visu_list})
add_custom_target(save_visu DEPENDS ${save_visu_list})
