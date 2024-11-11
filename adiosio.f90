module adiosio
    use adios2
    USE MPI

    implicit none

    !! Global variables
    type(adios2_adios) :: adios2obj
    type(adios2_io) :: adios2IO
    type(adios2_engine) :: adios2Engine

    type(adios2_variable) :: stepVar, arrVar

contains

subroutine init_adiosio()
    !use mpivars
    use adios2
    use MPI

    implicit none
    integer*4 :: err

    write(*,*) "----- DRP: init_adiosio"
    !MPI_COMM_WORLD
    !call adios2_init(adios2Obj, "xmlstring", ...)

    call adios2_init(adios2Obj, MPI_COMM_WORLD, err);


    !call adios2_init(adios2obj, app_comm, err)
    !if (.not.adios2obj%valid) then
    !    write(*,*) "Error when creating ADIOS2 object"
    !endif

    call adios2_declare_io (adios2IO, adios2obj, 'Output', err)
    call adios2_open(adios2Engine, adios2IO, "output.bp", adios2_mode_write, err)
    call adios2_define_variable(stepVar, adios2IO, "step", adios2_type_integer4, err)
    call adios2_define_variable(arrVar, adios2IO, "var", adios2_type_real, 1, (/ 20_8 /), (/ 16_8 /), (/4_8 /), .FALSE., err)

end subroutine init_adiosio

subroutine finalize_adiosio()
    use adios2
    implicit none
    integer*4 :: err
    ! close the engine now if it was not closed yet
    write(*,*) "----- DRP: finalize_adiosio"

    if (adios2Engine%valid) then
        call adios2_close(adios2Engine, err)
    endif
    if (adios2obj%valid) then
        call adios2_finalize(adios2obj, err)
    endif
end subroutine finalize_adiosio

end module adiosio
