find_package(doctest)

add_executable(ponio_tests main.cpp)
target_link_libraries(ponio_tests ponio doctest::doctest)

if (BUILD_SAMURAI_DEMOS)
    find_package(samurai)
    find_package(PkgConfig)
    pkg_check_modules(PETSC PETSc)

    if(PETSC_FOUND)
        find_package(MPI)
        add_definitions(-DBUILD_SAMURAI_DEMOS)

        target_link_libraries(ponio_tests samurai::samurai ${PETSC_LIBRARIES} ${MPI_LIBRARIES})
    endif()
endif()
