cmake_minimum_required(VERSION 2.6)
project(microscope)

find_package(OpenMP)
if(OPENMP_FOUND)
    add_definitions(${OpenMP_CXX_FLAGS})
    set(CMAKE_EXE_LINKER_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_MODULE_LINKER_FLAGS}")
endif()

find_library(GSL_LIBRARIES gsl)
find_library(GSL_CBLAS_LIBRARIES gslcblas)

add_executable(generate_table generate_table.cpp)
target_link_libraries(generate_table ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})
add_custom_command(
    OUTPUT ${PROJECT_BINARY_DIR}/born_wolf_psf_table.hpp
    COMMAND generate_table ${PROJECT_BINARY_DIR}/born_wolf_psf_table.hpp
    DEPENDS generate_table)
add_custom_target(psf_tables ALL DEPENDS ${PROJECT_BINARY_DIR}/born_wolf_psf_table.hpp)

add_executable(microscope main.cpp)
target_link_libraries(microscope ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})
add_dependencies(microscope psf_tables)

add_executable(main main.cpp)
target_link_libraries(main ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})
add_dependencies(main psf_tables)

add_executable(test test.cpp)
target_link_libraries(test ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})
add_dependencies(test psf_tables)
