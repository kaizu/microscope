find_package(OpenMP)
if(OPENMP_FOUND)
    add_definitions(${OpenMP_CXX_FLAGS})
    set(CMAKE_EXE_LINKER_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_MODULE_LINKER_FLAGS}")
endif()

find_library(GSL_LIBRARIES gsl)
find_library(GSL_CBLAS_LIBRARIES gslcblas)

find_library(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

add_executable(generate_table generate_table.cpp)
target_link_libraries(generate_table ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/born_wolf_psf_table.hpp
    COMMAND generate_table ${CMAKE_CURRENT_BINARY_DIR}/born_wolf_psf_table.hpp
    DEPENDS generate_table)
add_custom_target(psf_tables ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/born_wolf_psf_table.hpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_executable(main main.cpp)
target_link_libraries(main ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})
add_dependencies(main psf_tables)
install(TARGETS main DESTINATION bin)

# add_executable(test_microscope test.cpp)
# target_link_libraries(test_microscope ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})
# add_dependencies(test_microscope psf_tables)
# 
# add_executable(psf_property psf_property.cpp)
# target_link_libraries(psf_property ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})
# add_dependencies(psf_property psf_tables)

# cmake --help-module FindCUDA
find_package(CUDA)
if(CUDA_FOUND)
    cuda_add_library(microscope_gpu microscope_gpu.cu)
    add_dependencies(microscope_gpu psf_tables)
    add_executable(main_gpu main_gpu.cpp)
    target_link_libraries(main_gpu microscope_gpu ${CUDA_LIBRARIES} ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})
    install(TARGETS main_gpu DESTINATION bin)
endif()
