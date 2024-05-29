set(demos arenstorf brownian brusselator brusselator_dirk exp heat heat_rock pendulum lorenz lorenz_tuto lotka_volterra curtiss_hirschfelder curtiss_hirschfelder_exprk nagumo lorenz_all)
set(demos_eigen brusselator_dirk lorenz_tuto lorenz_all)
set(demos_cli11 brownian lotka_volterra)

set(demos_samurai heat_samurai nagumo_pirock belousov_zhabotinsky_pirock)

if(BUILD_SAMURAI_DEMOS)
    list(APPEND demos ${demos_samurai})
endif()

message(STATUS "list of examples:")
foreach(demo IN LISTS demos)
    message(STATUS "+ ${demo}")
    add_executable(${demo} ${demo}.cpp)

    if(${demo} IN_LIST demos_eigen)
        message(STATUS "    - eigen")
        find_package (Eigen3 3.4 REQUIRED NO_MODULE)
        target_link_libraries(${demo} Eigen3::Eigen)
    endif()

    if(${demo} IN_LIST demos_cli11)
        message(STATUS "    - cli11")
        find_package(CLI11)
        target_link_libraries(${demo} CLI11::CLI11)
    endif()

    if(${demo} IN_LIST demos_samurai)
        message(STATUS "    - samurai")
        find_package(samurai)
        target_link_libraries(${demo} samurai::samurai)

        find_package(PkgConfig)
        pkg_check_modules(PETSC PETSc)

        if(PETSC_FOUND)
            find_package(MPI)
            target_link_libraries(${demo} ${PETSC_LIBRARIES} ${MPI_LIBRARIES})
        endif()
    endif()

    target_link_libraries(${demo} ponio)
endforeach()

foreach(demo IN LISTS demos)
    add_custom_target(visu_${demo} python ${CMAKE_SOURCE_DIR}/ponio/examples/${demo}.py)
    add_custom_target(${demo}_visu make visu_${demo})
endforeach()
