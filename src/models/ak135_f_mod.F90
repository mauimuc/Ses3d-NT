!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Stefan Mauerberger
!!    and Maksym Melnyk
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
!! $Date: 2013-11-24 14:42:01 +0100 (Sun, 24 Nov 2013) $
!! $Author: mmelnyk $
!! $Revision: 806 $
!> @copyright GNU General Public License version 3 or later



!> This module uses the ak135 velocity model to define rho, vp, vs, qmu,
!! qkappa model parameters
!!
!! For further details please see the model websites:
!!  * http://rses.anu.edu.au/seismology/ak135/ak135f.html
!!  * http://rses.anu.edu.au/seismology/ak135/ak135t.html
!!
!! @todo continental structure
!! @todo ak135_t
!!
!! To specify the precision of floats this module expects to find a parameter
!! parameters_mod::real_kind in a module named parameters_mod.
MODULE ak135_f_mod
    USE parameters_mod, ONLY : real_kind

    IMPLICIT NONE

    PRIVATE

    !> Size of arrays holding values for ak135 full
    INTEGER, PARAMETER :: f = 145

    !> Size of arrays holding values for ak135 full
    INTEGER, PARAMETER :: c = 139

    !> Size of arrays holding values for ak135 full
    INTEGER, PARAMETER :: t = 145

    !> Depth
    REAL(real_kind), PARAMETER :: ak135_r(f) = [ &
        6371.00, 6368.00, 6368.00, 6367.70, 6367.70, 6361.00, &
        6361.00, 6353.00, 6353.00, 6328.00, 6291.00, 6291.00, &
        6251.00, 6251.00, 6206.00, 6161.00, 6161.00, 6111.00, &
        6061.00, 6011.00, 5961.00, 5961.00, 5911.00, 5861.00, &
        5811.00, 5761.00, 5711.00, 5711.00, 5661.00, 5611.00, &
        5561.50, 5512.00, 5462.50, 5413.00, 5363.50, 5314.00, &
        5264.50, 5215.00, 5165.50, 5116.00, 5066.50, 5017.00, &
        4967.50, 4918.00, 4868.50, 4819.00, 4769.50, 4720.00, &
        4670.50, 4621.00, 4571.50, 4522.00, 4472.50, 4423.00, &
        4373.50, 4324.00, 4274.50, 4225.00, 4175.50, 4126.00, &
        4076.50, 4027.00, 3977.50, 3928.00, 3878.50, 3829.00, &
        3779.50, 3731.00, 3681.00, 3631.00, 3631.00, 3581.33, &
        3531.67, 3479.50, 3479.50, 3431.67, 3381.34, 3331.01, &
        3280.68, 3230.34, 3180.01, 3129.68, 3079.35, 3029.02, &
        2978.69, 2928.36, 2878.03, 2827.70, 2777.36, 2727.03, &
        2676.70, 2626.37, 2576.04, 2525.71, 2475.38, 2425.05, &
        2374.72, 2324.38, 2274.05, 2223.72, 2173.39, 2123.06, &
        2072.73, 2022.40, 1972.07, 1921.74, 1871.40, 1821.07, &
        1770.74, 1720.41, 1670.08, 1619.75, 1569.42, 1519.09, &
        1468.76, 1418.42, 1368.09, 1317.76, 1267.43, 1217.50, &
        1217.50, 1166.39, 1115.68, 1064.96, 1014.25,  963.54, &
         912.83,  862.11,  811.40,  760.69,  709.98,  659.26, &
         608.55,  557.84,  507.13,  456.41,  405.70,  354.99, &
         304.28,  253.56,  202.85,  152.14,  101.43,   50.71, 0.00 ]

    !> Depth for continental structures
    REAL(real_kind), PARAMETER :: ak135_c_r(c) = [ &
        6371.0 - [0.0, 20.0, 20.0, 35.0, 35.0, 71.0], &
        ak135_r(13:)]

    !> Rho
    REAL(real_kind), PARAMETER :: ak135_rho(f) = [ &
         1.0200,  1.0200,  2.0000,  2.0000,  2.6000,  2.6000, &
         2.9200,  2.9200,  3.6410,  3.5801,  3.5020,  3.5020, &
         3.4268,  3.4268,  3.3711,  3.3243,  3.3243,  3.3663, &
         3.4110,  3.4577,  3.5068,  3.9317,  3.9273,  3.9233, &
         3.9218,  3.9206,  3.9201,  4.2387,  4.2986,  4.3565, &
         4.4118,  4.4650,  4.5162,  4.5654,  4.5926,  4.6198, &
         4.6467,  4.6735,  4.7001,  4.7266,  4.7528,  4.7790, &
         4.8050,  4.8307,  4.8562,  4.8817,  4.9069,  4.9321, &
         4.9570,  4.9817,  5.0062,  5.0306,  5.0548,  5.0789, &
         5.1027,  5.1264,  5.1499,  5.1732,  5.1963,  5.2192, &
         5.2420,  5.2646,  5.2870,  5.3092,  5.3313,  5.3531, &
         5.3748,  5.3962,  5.4176,  5.4387,  5.6934,  5.7196, &
         5.7458,  5.7721,  9.9145,  9.9942, 10.0722, 10.1485, &
        10.2233, 10.2964, 10.3679, 10.4378, 10.5062, 10.5731, &
        10.6385, 10.7023, 10.7647, 10.8257, 10.8852, 10.9434, &
        11.0001, 11.0555, 11.1095, 11.1623, 11.2137, 11.2639, &
        11.3127, 11.3604, 11.4069, 11.4521, 11.4962, 11.5391, &
        11.5809, 11.6216, 11.6612, 11.6998, 11.7373, 11.7737, &
        11.8092, 11.8437, 11.8772, 11.9098, 11.9414, 11.9722, &
        12.0001, 12.0311, 12.0593, 12.0867, 12.1133, 12.1391, &
        12.7037, 12.7289, 12.7530, 12.7760, 12.7980, 12.8188, &
        12.8387, 12.8574, 12.8751, 12.8917, 12.9072, 12.9217, &
        12.9351, 12.9474, 12.9586, 12.9688, 12.9779, 12.9859, &
        12.9929, 12.9988, 13.0036, 13.0074, 13.0100, 13.0117, 13.0122 ]

    !> V_s
    REAL(real_kind), PARAMETER :: ak135_vs(f) = [ &
        0.0000, 0.0000, 1.0000, 1.0000, 3.2000, 3.2000, &
        3.9000, 3.9000, 4.4839, 4.4856, 4.4800, 4.4900, &
        4.5000, 4.5000, 4.5090, 4.5184, 4.5184, 4.6094, &
        4.6964, 4.7832, 4.8702, 5.0806, 5.1864, 5.2922, &
        5.3989, 5.5047, 5.6104, 5.9607, 6.0898, 6.2100, &
        6.2424, 6.2799, 6.3164, 6.3519, 6.3860, 6.4182, &
        6.4514, 6.4822, 6.5131, 6.5431, 6.5728, 6.6009, &
        6.6285, 6.6554, 6.6813, 6.7070, 6.7323, 6.7579, &
        6.7820, 6.8056, 6.8289, 6.8517, 6.8743, 6.8972, &
        6.9194, 6.9416, 6.9625, 6.9852, 7.0069, 7.0286, &
        7.0504, 7.0722, 7.0932, 7.1144, 7.1368, 7.1584, &
        7.1804, 7.2031, 7.2253, 7.2485, 7.2485, 7.2593, &
        7.2700, 7.2817, 0.0000, 0.0000, 0.0000, 0.0000, &
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
        3.5043, 3.5187, 3.5314, 3.5435, 3.5551, 3.5661, &
        3.5765, 3.5864, 3.5957, 3.6044, 3.6126, 3.6202, &
        3.6272, 3.6337, 3.6396, 3.6450, 3.6498, 3.6540, &
        3.6577, 3.6608, 3.6633, 3.6653, 3.6667, 3.6675, 3.6678 ]

    !> V_p
    REAL(real_kind), PARAMETER :: ak135_vp(f) = [ &
         1.4500,  1.4500,  1.6500,  1.6500,  5.8000,  5.8000, &
         6.8000,  6.8000,  8.0355,  8.0379,  8.0400,  8.0450, &
         8.0505,  8.0505,  8.1750,  8.3007,  8.3007,  8.4822, &
         8.6650,  8.8476,  9.0302,  9.3601,  9.5280,  9.6962, &
         9.8640, 10.0320, 10.2000, 10.7909, 10.9222, 11.0553, &
        11.1355, 11.2228, 11.3068, 11.3897, 11.4704, 11.5493, &
        11.6265, 11.7020, 11.7768, 11.8491, 11.9208, 11.9891, &
        12.0571, 12.1247, 12.1912, 12.2558, 12.3181, 12.3813, &
        12.4427, 12.5030, 12.5638, 12.6226, 12.6807, 12.7384, &
        12.7956, 12.8524, 12.9093, 12.9663, 13.0226, 13.0786, &
        13.1337, 13.1895, 13.2465, 13.3017, 13.3584, 13.4156, &
        13.4741, 13.5311, 13.5899, 13.6498, 13.6498, 13.6533, &
        13.6570, 13.6601,  8.0000,  8.0382,  8.1283,  8.2213, &
         8.3122,  8.4001,  8.4861,  8.5692,  8.6496,  8.7283, &
         8.8036,  8.8761,  8.9461,  9.0138,  9.0792,  9.1426, &
         9.2042,  9.2634,  9.3205,  9.3760,  9.4297,  9.4814, &
         9.5306,  9.5777,  9.6232,  9.6673,  9.7100,  9.7513, &
         9.7914,  9.8304,  9.8682,  9.9051,  9.9410,  9.9761, &
        10.0103, 10.0439, 10.0768, 10.1095, 10.1415, 10.1739, &
        10.2049, 10.2329, 10.2565, 10.2745, 10.2854, 10.2890, &
        11.0427, 11.0585, 11.0718, 11.0850, 11.0983, 11.1166, &
        11.1316, 11.1457, 11.1590, 11.1715, 11.1832, 11.1941, &
        11.2041, 11.2134, 11.2219, 11.2295, 11.2364, 11.2424, &
        11.2477, 11.2521, 11.2557, 11.2586, 11.2606, 11.2618, 11.2622 ]

    !> Q_kappa
    REAL(real_kind), PARAMETER :: ak135_qkappa(f) = [ &
        57822.00, 57822.00,   163.35,   163.35,  1478.30,  1478.30,&
         1368.02,  1368.02,   950.50,   972.77,  1008.71,   182.03, &
          182.57,   182.57,   188.72,   200.97,   338.47,   346.37, &
          355.85,   366.34,   377.93,   413.66,   417.32,   419.94, &
          422.55,   425.51,   428.69,  1350.54,  1311.17,  1277.93, &
         1269.44,  1260.68,  1251.69,  1243.02,  1234.54,  1226.52, &
         1217.91,  1210.02,  1202.04,  1193.99,  1186.06,  1178.19, &
         1170.53,  1163.16,  1156.04,  1148.76,  1141.32,  1134.01, &
         1127.02,  1120.09,  1108.58,  1097.16,  1085.97,  1070.38, &
         1064.23,  1058.03,  1048.09,  1042.07,  1032.14,  1018.38, &
         1008.79,   999.44,   990.77,   985.63,   976.81,   968.46, &
          960.36,   952.00,   940.88,   933.21,   722.73,   726.87, &
          725.11,   723.12, 57822.00, 57822.00, 57822.00, 57822.00, &
        57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, &
        57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, &
        57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, &
        57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, &
        57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, &
        57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, &
        57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, &
          633.26,   629.89,   626.87,   624.08,   621.50,   619.71, &
          617.78,   615.93,   614.21,   612.62,   611.12,   609.74, &
          608.48,   607.31,   606.26,   605.28,   604.44,   603.69, &
          603.04,   602.49,   602.05,   601.70,   601.46,   601.32, 601.27 ]

    !> Q_mu
    REAL(real_kind), PARAMETER :: ak135_qmu(f) = [ &
          0.00,   0.00,  80.00,  80.00, 599.99, 599.99, &
        599.99, 599.99, 394.62, 403.93, 417.59,  75.60, &
         76.06,  76.06,  76.55,  79.40, 133.72, 136.38, &
        139.38, 142.76, 146.57, 162.50, 164.87, 166.80, &
        168.78, 170.82, 172.93, 549.45, 543.48, 537.63, &
        531.91, 526.32, 520.83, 515.46, 510.20, 505.05, &
        500.00, 495.05, 490.20, 485.44, 480.77, 476.19, &
        471.70, 467.29, 462.96, 458.72, 454.55, 450.45, &
        446.43, 442.48, 436.68, 431.03, 425.53, 418.41, &
        414.94, 411.52, 406.50, 403.23, 398.41, 392.16, &
        387.60, 383.14, 378.79, 375.94, 371.75, 367.65, &
        363.64, 359.71, 354.61, 350.88, 271.74, 273.97, &
        273.97, 273.97,   0.00,   0.00,   0.00,   0.00, &
          0.00,   0.00,   0.00,   0.00,   0.00,   0.00, &
          0.00,   0.00,   0.00,   0.00,   0.00,   0.00, &
          0.00,   0.00,   0.00,   0.00,   0.00,   0.00, &
          0.00,   0.00,   0.00,   0.00,   0.00,   0.00, &
          0.00,   0.00,   0.00,   0.00,   0.00,   0.00, &
          0.00,   0.00,   0.00,   0.00,   0.00,   0.00, &
          0.00,   0.00,   0.00,   0.00,   0.00,   0.00, &
         85.03,  85.03,  85.03,  85.03,  85.03,  85.03, &
         85.03,  85.03,  85.03,  85.03,  85.03,  85.03, &
         85.03,  85.03,  85.03,  85.03,  85.03,  85.03, &
         85.03,  85.03,  85.03,  85.03,  85.03,  85.03, 85.03 ]

    PUBLIC :: ak135_f

CONTAINS

    !> Returns the values of rho, vp, vs, qmu, qkappa corresponding to the
    !! passed `x` value using ak135_f velocity model
    !!
    !! @param x Radius corresponding to the depth [km]
    !! @param rho Density [Mg/km3]
    !! @param vp Velocity of the P-waves [km/s]
    !! @param vs Velocity of the S-waves [km/s]
    !! @param qmu Shear quality factor
    !! @param qkappa Bulk quality factor
    ELEMENTAL SUBROUTINE ak135_f( x, rho, vp, vs, qmu, qkappa )
        ! Passed dummy variables
        REAL(real_kind), INTENT(IN) :: x
        REAL(real_kind), INTENT(OUT), OPTIONAL :: rho, vp, vs, qmu, qkappa
        ! Local variables
        REAL(real_kind) :: delta
        INTEGER :: n, l, r

        IF ( ABS(x) > 6371.0_real_kind ) THEN
            IF ( PRESENT( rho ) ) rho = 0.0
            IF ( PRESENT( vp ) ) vp = 0.0
            IF ( PRESENT( vs ) ) vs = 0.0
            IF ( PRESENT( qkappa ) ) qkappa = HUGE(0.0_real_kind)
            IF ( PRESENT( qmu ) ) qmu = HUGE(0.0_real_kind)
        END IF

        ! Determine nearest index
        n = MINLOC( ABS( ak135_r(:) - ABS(x) ), DIM=1 )

        ! Find neighbouring indices
        IF ( ak135_r(n) < ABS(x) ) THEN
            l = n
            r = n - 1
        ELSE
            l = n + 1
            r = n
        END IF

        ! Calculate radial gradient for interpolation
        IF ( ak135_r(l) == ak135_r(r) ) THEN
            ! If x is located at a discontinuity, returned values will be
            ! totally arbitrary. This however is necessary to suppress
            ! singularities in the latter linear interpolation.
            delta = 0.0
        ELSE
            delta = ( x - ak135_r(l) ) / &
                    ( ak135_r(r) - ak135_r(l) )
        END IF


        ! Interpolate linearly and assign return values
        IF ( PRESENT( rho ) ) &
            rho = ak135_rho(l) + ( ak135_rho(r) - ak135_rho(l) ) * delta
        IF ( PRESENT( vp ) ) &
            vp = ak135_vp(l) + ( ak135_vp(r) - ak135_vp(l) ) * delta
        IF ( PRESENT( vs ) ) &
            vs = ak135_vs(l) + ( ak135_vs(r) - ak135_vs(l) ) * delta
        IF ( PRESENT( qmu ) ) &
            qmu = ak135_qmu(l) + ( ak135_qmu(r) - ak135_qmu(l) ) * delta
        IF ( PRESENT( qkappa ) ) &
            qkappa = ak135_qkappa(l) + &
                ( ak135_qkappa(r) - ak135_qkappa(l) ) * delta

    END SUBROUTINE ak135_f

END MODULE ak135_f_mod

