! compile with:
!   gfortran ../../src/examples/dummy_parameters.f90 ../../src/misc/coordinate_utilities.F90 ../../src/misc/string_utilities.f90 cart2sphere_HSP1a.f90 -o cart2sphere_HSPa1
PROGRAM cart2sphere
    USE parameters_mod, ONLY : real_kind, earth_radius
    USE string_utilities, ONLY : i2c
    USE coordinate_utilities, ONLY : xyz2tpr, rad2deg, colat2lat

    IMPLICIT NONE

    INTEGER :: i

    REAL(real_kind), PARAMETER, DIMENSION(12) :: &
        x = [0.0, 0.0, 0.0, 490.0, 3919.0, 7348.0, 400.0, 3200.0, 6000.0, 555.0, 4443.0, 8331.0], &
        y = [693.0, 5543.0, 10392.0, 490.0, 3919.0, 7348.0, 400.0, 3200.0, 6000.0, 370.0, 2952.0, 5554.0], &
        z = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 400.0, 3200.0, 6000.0, 184.0, 1481.0, 2777.0]

    REAL(real_kind), DIMENSION(12) :: t, p, r

    REAL(real_kind) lat, lon, depth
    CHARACTER(LEN=8) :: network='n', station='', location='n', attributes(3) = ['vx', 'vy', 'vz'], prefix='./'
    NAMELIST /receiver/ network, station, location, lat, lon, depth, attributes, prefix

    CALL xyz2tpr(y=x, z=y, x=earth_radius-20000.0-z, &
                 theta=t, phi=p, r=r)

    DO i=1, 12
        station = i2c(i)
        lat = colat2lat(rad2deg(t(i)))
        lon = rad2deg(p(i))
        depth = earth_radius-r(i)
        WRITE(*,NML=receiver)
    END DO

END PROGRAM cart2sphere
