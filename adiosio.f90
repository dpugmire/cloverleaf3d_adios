module adiosio
    use adios2
    USE MPI

    implicit none

    !! Global variables
    type(adios2_adios) :: adios2obj
    type(adios2_io) :: adios2IO
    type(adios2_engine) :: adios2Engine

    !! variables
    type(adios2_variable) :: step_var
    type(adios2_variable) :: coordsX_var, coordsY_var, coordsZ_var
    type(adios2_variable) :: density_var, energy_var, pressure_var, ghostzone_var
    type(adios2_variable) :: velocityX_var, velocityY_var, velocityZ_var

contains

subroutine init_adiosio()
    !use mpivars
    use adios2
    use MPI

    implicit none
    integer*4 :: err

    !write(*,*) "----- DRP: init_adiosio"
    !MPI_COMM_WORLD
    !call adios2_init(adios2Obj, "xmlstring", ...)

    call adios2_init(adios2Obj, "adios2.xml", MPI_COMM_WORLD, err);

    !call adios2_init(adios2obj, app_comm, err)
    !if (.not.adios2obj%valid) then
    !    write(*,*) "Error when creating ADIOS2 object"
    !endif

    call adios2_declare_io (adios2IO, adios2obj, 'Output', err)
    call adios2_open(adios2Engine, adios2IO, "output.bp", adios2_mode_write, err)
    call adios2_define_variable(step_var, adios2IO, "step", adios2_type_integer4, err)
    !XYZ coords
    call adios2_define_variable(coordsX_var, adios2IO, "coordsX", adios2_type_real4, 1, &
         (/ 20_8 /), (/ 16_8 /), (/4_8 /), .FALSE., err)
    call adios2_define_variable(coordsY_var, adios2IO, "coordsY", adios2_type_real4, 1, &
         (/ 20_8 /), (/ 16_8 /), (/4_8 /), .FALSE., err)
    call adios2_define_variable(coordsZ_var, adios2IO, "coordsZ", adios2_type_real4, 1, &
         (/ 20_8 /), (/ 16_8 /), (/4_8 /), .FALSE., err)

    !scalar values
    call adios2_define_variable(density_var, adios2IO, "density", adios2_type_real, 1, &
         (/ 20_8 /), (/ 16_8 /), (/4_8 /), .FALSE., err)
    call adios2_define_variable(energy_var, adios2IO, "energy", adios2_type_real, 1, &
         (/ 20_8 /), (/ 16_8 /), (/4_8 /), .FALSE., err)
    call adios2_define_variable(pressure_var, adios2IO, "pressure", adios2_type_real, 1, &
         (/ 20_8 /), (/ 16_8 /), (/4_8 /), .FALSE., err)
    call adios2_define_variable(ghostzone_var, adios2IO, "ghost_zones", adios2_type_real, 1, &
         (/ 20_8 /), (/ 16_8 /), (/4_8 /), .FALSE., err)

    !velocity
    call adios2_define_variable(velocityX_var, adios2IO, "velocityX", adios2_type_real, 1, &
         (/ 20_8 /), (/ 16_8 /), (/4_8 /), .FALSE., err)
    call adios2_define_variable(velocityY_var, adios2IO, "velocityY", adios2_type_real, 1, &
         (/ 20_8 /), (/ 16_8 /), (/4_8 /), .FALSE., err)
    call adios2_define_variable(velocityZ_var, adios2IO, "velocityZ", adios2_type_real, 1, &
         (/ 20_8 /), (/ 16_8 /), (/4_8 /), .FALSE., err)


end subroutine init_adiosio

subroutine finalize_adiosio()
    use adios2
    implicit none
    integer*4 :: err
    ! close the engine now if it was not closed yet
    !write(*,*) "----- DRP: finalize_adiosio"

    if (adios2Engine%valid) then
        call adios2_close(adios2Engine, err)
    endif
    if (adios2obj%valid) then
        call adios2_finalize(adios2obj, err)
    endif
end subroutine finalize_adiosio

subroutine put_adiosio1D(adiosVar, rank, numRanks, Nvals, adiosData)
    use adios2
    implicit none
    INTEGER*8, intent(in) :: Nvals
    integer*4, intent(in) :: rank, numRanks
    type(adios2_variable), intent(in) :: adiosVar
    real(kind=8), dimension(:), intent(in) :: adiosData  ! Define the type and dimensions for adiosData
    integer*4 :: err
    integer(kind=8), dimension(1) :: shape_dims, start_dims, count_dims
    integer, parameter :: dp = kind(1.0d0)    ! Double precision kind
    integer, parameter :: sp = kind(1.0)      ! Single precision kind
    real(kind=4), dimension(:), allocatable :: adiosData_single

    shape_dims(1) = Nvals * numRanks
    start_dims(1) = rank * Nvals
    count_dims(1) = Nvals

    CALL adios2_set_shape(adiosVar, 1, shape_dims, err)
    CALL adios2_set_selection(adiosVar, 1, start_dims, count_dims, err)

    !convert to single precision
    adiosData_single = real(adiosData, sp)
    CALL adios2_put(adios2Engine, adiosVar, adiosData_single, err)

end subroutine put_adiosio1D

subroutine put_adiosio3D(adiosVar, rank, numRanks, Nx, Ny, Nz, adiosData)
    use adios2
    implicit none
    INTEGER*8, intent(in) :: Nx, Ny, Nz
    integer*4, intent(in) :: rank, numRanks
    type(adios2_variable), intent(in) :: adiosVar
    real(kind=8), dimension(:,:,:), intent(in) :: adiosData  ! Define the type and dimensions for adiosData
    integer*4 :: err
    integer(kind=8), dimension(3) :: shape_dims, start_dims, count_dims
    real(kind=8), pointer :: array1d(:)
    integer, parameter :: dp = kind(1.0d0)    ! Double precision kind
    integer, parameter :: sp = kind(1.0)      ! Single precision kind
    real(kind=4), dimension(:,:,:), allocatable :: adiosData_single

    shape_dims(1) = Nx * numRanks
    shape_dims(2) = Ny * numRanks
    shape_dims(3) = Nz * numRanks
    start_dims(1) = rank * Nx
    start_dims(2) = rank * Ny
    start_dims(3) = rank * Nz
    count_dims(1) = Nx
    count_dims(2) = Ny
    count_dims(3) = Nz

    CALL adios2_set_shape(adiosVar, 3, shape_dims, err)
    CALL adios2_set_selection(adiosVar, 3, start_dims, count_dims, err)

    !convert to single precision
    adiosData_single = real(adiosData, sp)
    CALL adios2_put(adios2Engine, adiosVar, adiosData_single, err)

end subroutine put_adiosio3D

end module adiosio
