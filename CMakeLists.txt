cmake_minimum_required(VERSION 3.26)

project(cloverLeaf3d LANGUAGES C CXX Fortran)

#set(CMAKE_MODULE_PATH /Users/dpn/sw/adios2/master/install/lib/cmake/adios2)
#set(ADIOS_DIR /Users/dpn/sw/adios2/master/install/lib/cmake/adios2)
find_package(MPI REQUIRED)
find_package(ADIOS2 REQUIRED)

set(clover_sources
    PdV.f90
    PdV_kernel.f90
    accelerate.f90
    accelerate_kernel.f90
    advec_cell_driver.f90
    advec_cell_kernel.f90
    advec_mom_driver.f90
    advec_mom_kernel.f90
    advection.f90
    build_field.f90
    calc_dt.f90
    calc_dt_kernel.f90
    clover.F90
    clover_leaf.f90
    data.f90
    definitions.f90
    field_summary.f90
    field_summary_kernel.f90
    flux_calc.f90
    flux_calc_kernel.f90
    generate_chunk.f90
    generate_chunk_kernel.f90
    hydro.f90
    ideal_gas.f90
    ideal_gas_kernel.f90
    initialise.f90
    initialise_chunk.f90
    initialise_chunk_kernel.f90
    pack_kernel.f90
    parse.f90
    read_input.f90
    report.f90
    reset_field.f90
    reset_field_kernel.f90
    revert.f90
    revert_kernel.f90
    start.f90
    timer.f90
    timestep.f90
    update_halo.f90
    update_halo_kernel.f90
    viscosity.f90
    viscosity_kernel.f90
    visit.F90
    timer_c.c
    clover_main.cpp
    adiosio.f90)

set(CMAKE_Fortran_FLAGS "-DUSE_MOD")

add_executable(cloverLeaf3d ${clover_sources})


target_include_directories(cloverLeaf3d PRIVATE ${ADIOS_INCLUDE_DIRS})
target_link_libraries(cloverLeaf3d PRIVATE ${ADIOS_LIBRARIES})

target_compile_definitions(cloverLeaf3d PRIVATE ADIOS2)
target_link_libraries(cloverLeaf3d  adios2::fortran_mpi MPI::MPI_Fortran)
