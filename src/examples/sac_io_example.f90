!> @example sac_io_example.f90
!!
!! SAC I/O Example Program. It describes how the modules sac_io_mod 
!! and sac_header_mod can be used. Invoking the sac_cls structure constructors
!! we can instantiate and read SAC-objects.
!!
!! You can use Python to visualize the generated SAC file
!! @code{.py} 
!!import obspy
!!st = obspy.read('./*.sac')
!!st.plot( )
!! @endcode
!!
!! Compile by simply using the Makefile:
!!
!!     make sac_io_example
!!
!! ...
PROGRAM sac_io_example
    USE sac_io_mod, ONLY : sac_cls ! Make the sac_cls class available
    USE sac_header_mod, ONLY : rk, & ! Real kind
                               enum ! The SAC Enumerated Header Variables

    IMPLICIT NONE

    INTEGER :: i, npts
    REAL(rk) :: dt
    REAL(rk), ALLOCATABLE :: data1(:)

    ! Definition of a SAC-object
    TYPE(sac_cls) :: sac_obj


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Creating an writing SAC-objects !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Just a very simple data section 
    data1 = [ (REAL(i, rk), i=1, 100) ] 

    ! Instantiate a SAC-object by invoking its structure constructor
    ! sac_cls(...) is not only a type specifier but also
    ! a generic (function) interface witch returns a sac_cls object
    sac_obj = sac_cls(lat      =1.0_rk, & ! Station latitude
                      lon      =1.0_rk, & ! Station longitude
                      depth    =1.0_rk, & ! Station depth 
                      data     =data1,  & ! Actual data; array of rank 1
                      dt       =0.1_rk, & ! Time increment
                      network  ='a', & ! Network name; up to 8 chars
                      station  ='b', & ! Station name; up to 8 chars
                      location ='c', & ! Location; up to 8 chars
                      channel  ='d', & ! Channel; up to 8 chars
                      event    ='e', & ! Event name; up to 16 chars
                      override =.TRUE., & ! Whether to override existing files
                      date_time=[1,1,1,1,1,1,1,1]) ! yyyy,mm,dd,gmt,h,m,s,ms

    ! Writs the SAC-trac (arguments iostat, iomsg are optional)
    ! Default file name is "./[network].[station].[location].[channel].sac"
    CALL sac_obj%write()
    ! Passing a file-name is optional
    CALL sac_obj%write('./test.sac')



    !!!!!!!!!!!!!!!!!!!!!!!
    ! Reading SAC-objects !
    !!!!!!!!!!!!!!!!!!!!!!!

    ! Read and instantiate sac-object using one of its structure constructors
    sac_obj = sac_cls( './test.sac' )



    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Accessing components !
    !!!!!!!!!!!!!!!!!!!!!!!!

    ! Prints header to stdout; a unit identifier is optional
    CALL sac_obj%print()
    ! Prints data to std-out
    PRINT *, sac_obj%get_data1()
    ! Ge the number of sampling-points
    npts = sac_obj%get_header('npts') ! Case-insensitive
    dt = sac_obj%get_header('delta') ! Case-insensitive
    ! Print npts and dt to stdout
    PRINT *, 'npts=', npts, 'dt=', dt



    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Modifying components !
    !!!!!!!!!!!!!!!!!!!!!!!!

    ! Make data=sin(time) where time is the time vector
    dt = 0.12; npts = 154
    data1 = SIN( [ (REAL(i, rk), i=1, npts) ]*dt )

    ! Assign data 
    CALL sac_obj%set_header('NPTS', npts)
    CALL sac_obj%set_header('DELTA', dt)
    CALL sac_obj%set_data1(data1)
    ! Mark override flag as true
    CALL sac_obj%set_header('LOVROK', .TRUE.)
    ! Set network, location, station, channel and event
    CALL sac_obj%set_header('KNETWK', 'network')
    CALL sac_obj%set_header('KHOLE', 'location')
    CALL sac_obj%set_header('KSTNM', 'station')
    CALL sac_obj%set_header('KCMPNM', 'channel')
    CALL sac_obj%set_header('KEVNM', '1234567891011131') 
    ! Mark override flag as true
    CALL sac_obj%set_header('LEVEN', .TRUE.)
    ! Write SAC-object to file "./sin.sac"
    CALL sac_obj%write( './sin.sac' )


END PROGRAM sac_io_example


!> @file
!! $Date: 2013-10-30 16:25:25 +0100 (Wed, 30 Oct 2013) $
!! $Author: mauerberger $
!! $Revision: 743 $
!! @copyright GNU General Public License version 3 or later 
