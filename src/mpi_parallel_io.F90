!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Stefan Mauerberger <mauerberger@geophysik.uni-muenchen.de>
!!
!! This program is free software: you can redistribute it and/or modify
!! under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!! $Date: 2013-11-16 10:55:05 +0100 (Sat, 16 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 785 $
!! @copyright GNU General Public License version 3 or later



!> MPI Parallel I/O
!!
!! Procedures included are ought to preform asynchronous parallel file in- and 
!! output. 
MODULE mpi_parallel_io_mod
    USE parameters_mod, ONLY : real_kind, isl, is_mpi_master_rank, &
                               my_mpi_real, init_my_mpi_real, &
                               my_mpi_integer, init_my_mpi_integer
    USE configuration_mod, ONLY : conf => configuration
! line control through pre-processor directives
#ifndef MPI_INCLUDE
    !! RECOMMENDED !!
    ! in case MPI_INCLUDE is NOT defined the mpi module will be used 
    USE mpi
#endif

    IMPLICIT NONE
    
#ifdef MPI_INCLUDE
        !! DEPRECATED !! 
        ! when MPI_INCLUDE is defined, mpi header-files will be included 
        INCLUDE 'mpif.h'
#endif

    PRIVATE

    !> MPI data type distributed array (handle)
    INTEGER, ALLOCATABLE :: my_mpi_darray

    PUBLIC :: mpi_write_parallel, mpi_write_parallel_begin, mpi_write_parallel_end, &
              mpi_read_parallel, mpi_read_parallel_begin, mpi_read_parallel_end

CONTAINS

    !> Creates and commits a MPI distributed array datatype
    SUBROUTINE init_my_mpi_darray()
        ! Distribution of array in each dimension (array of state)
        INTEGER :: distrib(6) = [ MPI_DISTRIBUTE_BLOCK, &
                                  MPI_DISTRIBUTE_BLOCK, &
                                  MPI_DISTRIBUTE_BLOCK, &
                                  MPI_DISTRIBUTE_NONE, & 
                                  MPI_DISTRIBUTE_NONE, &
                                  MPI_DISTRIBUTE_NONE ]
        ! Distribution argument in each dimension (array of positive integers)
        INTEGER :: drags(6) = [ MPI_DISTRIBUTE_DFLT_DARG, & 
                                MPI_DISTRIBUTE_DFLT_DARG, &
                                MPI_DISTRIBUTE_DFLT_DARG, &
                                MPI_DISTRIBUTE_DFLT_DARG, &
                                MPI_DISTRIBUTE_DFLT_DARG, &
                                MPI_DISTRIBUTE_DFLT_DARG ]
        INTEGER :: mpi_err

        ! In case, initialize my_mpi_real
        IF ( my_mpi_real == -1 ) &
            CALL init_my_mpi_real()

        ! Allocate my_mpi_darray
        IF ( .NOT. ALLOCATED(my_mpi_darray) ) &
            ALLOCATE( my_mpi_darray )

        ASSOCIATE( px => conf%px(), & 
                   py => conf%py(), & 
                   pz => conf%pz()  )

        ! Creates a 6d distributed array 
        CALL MPI_TYPE_CREATE_DARRAY( conf%mpi_size(), & ! n_cpus
                                     conf%mpi_rank(), & ! MPI-rank
                                     6, & ! array-rank
                                     conf%shape_glb(), & ! global shape
                                     distrib, & 
                                     drags, & 
                                     [ px, py, pz, 1, 1, 1 ], &
                                     MPI_ORDER_FORTRAN, &
                                     my_mpi_real, & ! kind type
                                     my_mpi_darray, & ! handle
                                     mpi_err ) 

        ! Initialize distributed array
        ! XXX very important question:
        !   may my_mpi_darray be used simultaneously for IO? 
        CALL MPI_TYPE_COMMIT( my_mpi_darray, mpi_err ) 

        END ASSOCIATE

    END SUBROUTINE init_my_mpi_darray



    !> Beginning part of a split collective routine (nonblocking) writing 
    !! the passed field to a file. It returns a file handle.
    !!
    !! @param file_name Filename 
    !! @param filed The field must no be altered while writing 
    !! @param override If true existing files are overwritten 
    !! @result A MPI file handle 
    FUNCTION mpi_write_parallel_begin( file_name, field, override ) &
        RESULT( mpi_file_handle )
        USE error_mod, ONLY : warn
        CHARACTER(LEN=*), INTENT(IN):: file_name
        REAL(real_kind), INTENT(IN) :: field(:,:,:,:,:,:)
        LOGICAL, INTENT(IN), OPTIONAL :: override
        INTEGER :: mpi_file_handle, amode, mpi_err, ierr2, str_len
        CHARACTER(LEN=isl) :: mpi_err_string

        ! In case initialize my_mpi_darray 
        IF ( .NOT. ALLOCATED(my_mpi_darray) ) &
            CALL init_my_mpi_darray()

        ! Set the file access mode to non-override 
        amode = MPI_MODE_CREATE+MPI_MODE_WRONLY+MPI_MODE_EXCL
        ! In case override == true 
        IF ( PRESENT(override) ) THEN
            IF ( override .EQV. .TRUE. ) &
                amode = MPI_MODE_CREATE+MPI_MODE_WRONLY
        END IF

        ! Opens a file (collective)
        CALL MPI_FILE_OPEN( conf%mpi_comm_cart(), &
                            file_name, &
                            amode, &
                            MPI_INFO_NULL, &
                            mpi_file_handle, &
                            mpi_err )

        IF ( mpi_err /= MPI_SUCCESS ) THEN
            CALL MPI_ERROR_STRING( mpi_err, mpi_err_string, str_len, ierr2 )
            mpi_err_string = "Error while writing raw output " // & 
                TRIM(file_name) // ' - ' // mpi_err_string(1:str_len)
            IF ( conf%log_file_true() ) &
                CALL warn( TRIM(mpi_err_string), conf%log_file_unit )
            IF ( is_mpi_master_rank() ) &
                CALL warn( TRIM(mpi_err_string) )
            RETURN
        END IF

        ! Set all MPI_FILE_x errors being treated fatal
        CALL MPI_FILE_SET_ERRHANDLER( mpi_file_handle, &
                                      MPI_ERRORS_ARE_FATAL, &
                                      mpi_err )

        ! Changes process's view of data in file (collective)
        CALL MPI_FILE_SET_VIEW( mpi_file_handle, &
                                0_MPI_OFFSET_KIND, &
                                my_mpi_real, &
                                my_mpi_darray, &
                                "native", &
                                MPI_INFO_NULL, &
                                mpi_err )

        ! Starts writing field to file_name (nonblocking)
        CALL MPI_FILE_WRITE_ALL_BEGIN( mpi_file_handle, &
                                       field, &
                                       PRODUCT( SHAPE( field ) ), & 
                                       my_mpi_real, &
                                       mpi_err ) 

    END FUNCTION mpi_write_parallel_begin

    !> Ending part of the split collective write routine (blocking). 
    !!
    !! @param mpi_file_handle A MPI file handle 
    !! @param filed The corresponding field (not altered)
    SUBROUTINE mpi_write_parallel_end( mpi_file_handle, field )
        REAL(real_kind), INTENT(IN) :: field(:,:,:,:,:,:)
        INTEGER :: mpi_file_handle 
        INTEGER :: mpi_status(MPI_STATUS_SIZE), mpi_err

        ! Waits till file has been written (blocking)
        CALL MPI_FILE_WRITE_ALL_END( mpi_file_handle, &
                                     field, &
                                     mpi_status, &
                                     mpi_err )

        ! Closes the file (collective)
        CALL MPI_FILE_CLOSE( mpi_file_handle, mpi_err )

    END SUBROUTINE mpi_write_parallel_end

    !> Convenience procedure to put mpi_write_parallel_begin() and 
    !! mpi_write_parallel_end() into one single call.
    !!
    !! @param file_name Filename
    !! @param field Field to write 
    !! @param override If true existing files are overwritten 
    SUBROUTINE mpi_write_parallel( file_name, field, override )
        CHARACTER(LEN=*), INTENT(IN):: file_name
        REAL(real_kind), INTENT(IN) :: field(:,:,:,:,:,:)
        LOGICAL, INTENT(IN), OPTIONAL :: override
        INTEGER :: fh

        ! Initiate writing field to file file_name 
        fh = mpi_write_parallel_begin( file_name, field, override )

        ! Wait till file has been written
        CALL mpi_write_parallel_end( fh, field )

    END SUBROUTINE mpi_write_parallel



    !> Beginning part of a split collective routine (nonblocking) reading 
    !! from a file into the passed field. 
    !!
    !! @param file_name Filename 
    !! @param field Field
    !! @result A MPI file handle 
    FUNCTION mpi_read_parallel_begin(file_name, field) RESULT(mpi_file_handle)
        CHARACTER(LEN=*), INTENT(IN):: file_name
        REAL(real_kind), INTENT(INOUT) :: field(:,:,:,:,:,:)
        INTEGER :: mpi_file_handle
        INTEGER :: mpi_err

        IF ( .NOT. ALLOCATED(my_mpi_darray) ) &
            CALL init_my_mpi_darray()

        CALL MPI_FILE_OPEN( conf%mpi_comm_cart(), &
                            file_name, &
                            MPI_MODE_RDONLY, &
                            MPI_INFO_NULL, &
                            mpi_file_handle, &
                            mpi_err )

        ! Set all MPI_FILE_x errors being treated fatal
        CALL MPI_FILE_SET_ERRHANDLER( mpi_file_handle, &
                                      MPI_ERRORS_ARE_FATAL, &
                                      mpi_err )

        CALL MPI_FILE_SET_VIEW( mpi_file_handle, &
                                0_MPI_OFFSET_KIND, &
                                my_mpi_real, &
                                my_mpi_darray, &
                                "native", &
                                MPI_INFO_NULL, &
                                mpi_err )

        CALL MPI_FILE_READ_ALL_BEGIN( mpi_file_handle, &
                                      field, &
                                      PRODUCT( SHAPE( field ) ), & 
                                      my_mpi_real, &
                                      mpi_err ) 

    END FUNCTION mpi_read_parallel_begin

    !> Ending part of the split collective read routine (blocking). 
    !!
    !! @param mpi_file_handle A MPI file handle 
    !! @param filed The corresponding field 
    SUBROUTINE mpi_read_parallel_end( mpi_file_handle, field )
        USE error_mod, ONLY : abort
        REAL(real_kind), INTENT(INOUT) :: field(:,:,:,:,:,:)
        INTEGER, INTENT(INOUT) :: mpi_file_handle 
        INTEGER :: mpi_status(MPI_STATUS_SIZE), mpi_err, &
                   local_counts, global_counts, my_mpi_real_size
        INTEGER(MPI_OFFSET_KIND) :: file_size

        ! Waits till file has been read (blocking)
        CALL MPI_FILE_READ_ALL_END( mpi_file_handle, &
                                    field, &
                                    mpi_status, &
                                    mpi_err )

        ! Returns the number of entries received
        CALL MPI_GET_COUNT( mpi_status, &
                            my_mpi_real, &
                            local_counts, &
                            mpi_err )

        ! Sums up values from all processes; Total number of entities received
        CALL MPI_ALLREDUCE( local_counts, &
                            global_counts, &
                            1, &
                            MPI_INTEGER, &
                            MPI_SUM, &
                            conf%mpi_comm_cart(), &
                            mpi_err ) 
        ! Returns the number of bytes occupied my_mpi_real_size
        CALL MPI_TYPE_SIZE( my_mpi_real, my_mpi_real_size, mpi_err )

        ! Returns the size of the file
        CALL MPI_FILE_GET_SIZE( mpi_file_handle, file_size, mpi_err )

        IF ( file_size /= global_counts*my_mpi_real_size ) THEN
            CALL abort( "Could not read file, wrong file size!" )
        END IF

        ! Closes the file (collective)
        CALL MPI_FILE_CLOSE( mpi_file_handle, mpi_err )
        
    END SUBROUTINE mpi_read_parallel_end

    !> Convenience procedure to put mpi_read_parallel_begin() and 
    !! mpi_read_parallel_end() into one single call.
    !!
    !! @bug GCC bug: field has to allocatable 
    !! (see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=57966)
    !!
    !! @param file_name Filename 
    !! @param field Field to write 
    !! @param override If true existing files are overwritten 
    FUNCTION mpi_read_parallel( file_name ) RESULT( field )
        CHARACTER(LEN=*), INTENT(IN) :: file_name
        REAL(real_kind), ALLOCATABLE :: field (:,:,:,:,:,:)
        INTEGER :: fh
       
        ! Unfortunately an allocatable needs to be used 
        ! placing conf%... at initialization does not work
        ALLOCATE( field(0:conf%nx_loc(),0:conf%ny_loc(),0:conf%nz_loc(),&
                        0:conf%lpd(),   0:conf%lpd(),   0:conf%lpd()  ) )

        ! Initiate writing field to file file_name 
        fh = mpi_read_parallel_begin( file_name, field )

        ! Wait till file has been written
        CALL mpi_read_parallel_end( fh, field )

    END FUNCTION mpi_read_parallel

END MODULE mpi_parallel_io_mod
