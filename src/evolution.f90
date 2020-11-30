! SES3D - simulation of elastic wave propagation in spherical sections

! (c) by Andreas Fichtner
!        Tarje Nissen-Meyer
!        Heiner Igel
!        Stefan Mauerberger

! This program is free software: you can redistribute it and/or modify
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Last changed: $Date: 2013-09-07 19:27:29 +0200 (Sat, 07 Sep 2013) $
! By: $Author: mauerberger $
! Revision: $Revision: 669 $

!------------------------------------------------------------------------------
! Time evolution
!------------------------------------------------------------------------------

MODULE evolution_mod
    USE parameters_mod
    USE configuration_mod, ONLY : config => configuration

    IMPLICIT NONE

    PRIVATE
    REAL(real_kind), ALLOCATABLE :: dl(:,:)
    REAL(real_kind) :: dt, dx_elm, dy_elm, dz_elm, ispml
    INTEGER :: lpd, nx, ny, nz
    LOGICAL :: anisotropy, visco_elasticity

    PUBLIC :: ses3d_evolution, init_temp

CONTAINS

!======================================================================
! initialisations this is just a temporal hack
!======================================================================
SUBROUTINE init_temp
    USE gll_mod, only : get_dlgll
    USE model_paras_mod, ONLY : a, b, c, qq
    INTEGER :: i,j

    ! Just a few temporal hacks
    ispml = 0.0
    dt = config%dt()
    dx_elm = config%dx_elm()
    dy_elm = config%dy_elm()
    dz_elm = config%dz_elm()
    lpd = config%lpd()
    nx = config%nx_loc()
    ny = config%ny_loc()
    nz = config%nz_loc()

    IF ( ANY( [a,b,c] /= 0.0_real_kind ) ) THEN
        anisotropy = .TRUE.
    ELSE
        anisotropy = .FALSE.
    END IF

    IF ( ANY( qq /= HUGE(0.0_real_kind ) ) ) THEN
        visco_elasticity = .TRUE.
    ELSE
        visco_elasticity = .FALSE.
    END IF

    ! Initialisation of LAGRANGE polynomial derivatives
    ALLOCATE( dl(0:config%lpd(),0:config%lpd()) )
    DO i=0,config%lpd()
        DO j=0,config%lpd()
            dl(i,j) = get_dlgll( config%lpd(), i, j )
        ENDDO
    ENDDO

END SUBROUTINE init_temp


SUBROUTINE ses3d_evolution
    USE elastic_vars_mod, ONLY : vx, vy, vz, sx, sy, sz, &
                                 exx, eyy, ezz, exy, eyx, &
                                 exz, ezx, ezy, eyz, &
                                 sxx_pml, sxy_pml, sxz_pml, &
                                 syy_pml, syz_pml, szz_pml, &
                                 dxux, dxuy, dxuz, dyux, dyuy, dyuz, &
                                 dzux, dzuy, dzuz, &
                                 !Mxx, Myy, Mzz, Mxz, Mxy, Myz, &
                                 src_x, src_y, src_z, src_xx, src_xy, src_xz, &
                                 src_yy, src_yz, src_zz
    USE geometric_paras_mod, ONLY : cot_theta, r, r_sin_theta, rr_sin_theta, &
                                    !tau_p, prof_x, prof_y, prof_z, prof, &
                                    taper, jac, wx_wy_wz, nrdiss
    USE model_paras_mod, ONLY : lambda, two_mu, a, b, c, mm!, kappa, tau, mu_tau
    USE communicate_fields_mod, ONLY : communicate_local_boundaries, &
                                       communicate_global_boundaries
    INTEGER :: i, j, k, q
    REAL(real_kind) :: dummy_1, dummy_2, dummy_3
    REAL(real_kind) :: dummy_x(0:nx,0:ny,0:nz), &
                       dummy_y(0:nx,0:ny,0:nz), &
                       dummy_z(0:nx,0:ny,0:nz)
    REAL(real_kind) :: L(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)!, &
                       !M_dummy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
                       !Mdx(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
                       !Mdy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
                       !Mdz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)
    REAL(real_kind) :: sxx(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
                       sxy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
                       sxz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
                       syy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
                       syz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
                       szz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)


    !----------------------------------------------------------------------
    !- derivative proxies of velocity field
    !----------------------------------------------------------------------

    exx(:,:,:,:,:,:)=0.0; exy(:,:,:,:,:,:)=0.0; exz(:,:,:,:,:,:)=0.0
    eyx(:,:,:,:,:,:)=0.0; eyy(:,:,:,:,:,:)=0.0; eyz(:,:,:,:,:,:)=0.0
    ezx(:,:,:,:,:,:)=0.0; ezy(:,:,:,:,:,:)=0.0; ezz(:,:,:,:,:,:)=0.0

    DO i=0,lpd
        DO q=0,lpd

            dummy_3 = 2.0*dl(q,i)/dz_elm
            dummy_2 = 2.0*dl(q,i)/dy_elm
            dummy_1 = 2.0*dl(q,i)/dx_elm

            ! (grad v)_(r r)
            ezz(:,:,:,:,:,i)=ezz(:,:,:,:,:,i)+vz(:,:,:,:,:,q)*dummy_3
            ! (grad v)_(r phi)
            ezy(:,:,:,:,:,i)=ezy(:,:,:,:,:,i)+vy(:,:,:,:,:,q)*dummy_3
            ! (grad v)_(r theta)
            ezx(:,:,:,:,:,i)=ezx(:,:,:,:,:,i)+vx(:,:,:,:,:,q)*dummy_3

            ! (grad v)_(phi phi)
            eyy(:,:,:,:,i,:)=eyy(:,:,:,:,i,:)+vy(:,:,:,:,q,:)*dummy_2
            ! (grad v)_(phi theta)
            eyx(:,:,:,:,i,:)=eyx(:,:,:,:,i,:)+vx(:,:,:,:,q,:)*dummy_2
            ! (grad v)_(phi r)
            eyz(:,:,:,:,i,:)=eyz(:,:,:,:,i,:)+vz(:,:,:,:,q,:)*dummy_2

            ! (grad v)_(theta theta)
            exx(:,:,:,i,:,:)=exx(:,:,:,i,:,:)+vx(:,:,:,q,:,:)*dummy_1
            ! (grad v)_(theta phi)
            exy(:,:,:,i,:,:)=exy(:,:,:,i,:,:)+vy(:,:,:,q,:,:)*dummy_1
            ! (grad v)_(theta r)
            exz(:,:,:,i,:,:)=exz(:,:,:,i,:,:)+vz(:,:,:,q,:,:)*dummy_1

        ENDDO
    ENDDO

    !-------------------------------------------------------------------------
    !- complete velocity gradients in spherical coordinates ------------------
    !-------------------------------------------------------------------------

    ezz(:,:,:,:,:,:) = -ezz                                  ! (grad v)_(r r)
    eyy(:,:,:,:,:,:) = eyy/(r_sin_theta)+(vz+vx*cot_theta)/r ! (grad v)_(phi phi)
    exx(:,:,:,:,:,:) = (vz+exx)/r                            ! (grad v)_(theta theta)
    ezx(:,:,:,:,:,:) = -ezx                                  ! (grad v)_(r theta)
    exz(:,:,:,:,:,:) = (exz-vx)/r                            ! (grad v)_(theta r)
    ezy(:,:,:,:,:,:) = -ezy                                  ! (grad v)_(r phi)
    eyz(:,:,:,:,:,:) = -vy/r+eyz/(r_sin_theta)               ! (grad v)_(phi r)
    exy(:,:,:,:,:,:) = exy/r                                 ! (grad v)_(theta phi)
    eyx(:,:,:,:,:,:) = -vy*cot_theta/r+eyx/(r_sin_theta)     ! (grad v)_(phi theta)

    !--------------------------------------------------------------------------
    !- integrate velocity gradients to displacement gradients -----------------
    !--------------------------------------------------------------------------

    dzuz(:,:,:,:,:,:)=dzuz+dt*ezz
    dyuy(:,:,:,:,:,:)=dyuy+dt*eyy
    dxux(:,:,:,:,:,:)=dxux+dt*exx
    dzux(:,:,:,:,:,:)=dzux+dt*ezx
    dxuz(:,:,:,:,:,:)=dxuz+dt*exz
    dzuy(:,:,:,:,:,:)=dzuy+dt*ezy
    dyuz(:,:,:,:,:,:)=dyuz+dt*eyz
    dxuy(:,:,:,:,:,:)=dxuy+dt*exy
    dyux(:,:,:,:,:,:)=dyux+dt*eyx

    !--------------------------------------------------------------------------
    !- build apml strain rates ------------------------------------------------
    !--------------------------------------------------------------------------

    !exx(:,:,:,:,:,:)=exx+(prof_z+prof_y)*dxux*ispml
    !eyy(:,:,:,:,:,:)=eyy+(prof_x+prof_z)*dyuy*ispml
    !ezz(:,:,:,:,:,:)=ezz+(prof_x+prof_y)*dzuz*ispml

    !exy(:,:,:,:,:,:)=exy+(prof_y+prof_z)*dxuy*ispml
    !eyx(:,:,:,:,:,:)=eyx+(prof_x+prof_z)*dyux*ispml

    !exz(:,:,:,:,:,:)=exz+(prof_y+prof_z)*dxuz*ispml
    !ezx(:,:,:,:,:,:)=ezx+(prof_x+prof_y)*dzux*ispml

    !eyz(:,:,:,:,:,:)=eyz+(prof_x+prof_z)*dyuz*ispml
    !ezy(:,:,:,:,:,:)=ezy+(prof_x+prof_y)*dzuy*ispml

    exy(:,:,:,:,:,:)=(exy+eyx)/2.0
    exz(:,:,:,:,:,:)=(exz+ezx)/2.0
    eyz(:,:,:,:,:,:)=(eyz+ezy)/2.0

    !======================================================================
    ! make strong form stress rates
    !======================================================================

    ! - diagonal stress rates ---------------------------------------------

    !IF (is_diss==1) THEN    !- dissipation on

    !    !Mdx(:,:,:,:,:,:) =0.0
    !    !Mdy(:,:,:,:,:,:) =0.0
    !    !Mdz(:,:,:,:,:,:) =0.0

    !    !DO k=1,nrdiss

    !    !    Mdx(:,:,:,:,:,:) =Mdx+Mxx(k,:,:,:,:,:,:)
    !    !    Mdy(:,:,:,:,:,:) =Mdy+Myy(k,:,:,:,:,:,:)
    !    !    Mdz(:,:,:,:,:,:) =Mdz+Mzz(k,:,:,:,:,:,:)

    !    !ENDDO

    !    Mdx = SUM( Mxx, dim=1 )
    !    Mdy = SUM( Myy, dim=1 )
    !    Mdz = SUM( Mzz, dim=1 )

    !    M_dummy(:,:,:,:,:,:)=Mdx+Mdy+Mdz

    !    L(:,:,:,:,:,:) = (kappa-2.0*mu_tau/3.0)*(exx+eyy+ezz)-two_mu*M_dummy/3.0

    !    !IF (is_aniso==1) THEN    !- anisotropy on

    !    sxx(:,:,:,:,:,:) = L + 2.0*mu_tau*exx + two_mu*Mdx + C*ezz + A*(eyy+exx)
    !    syy(:,:,:,:,:,:) = L + 2.0*mu_tau*eyy + two_mu*Mdy + C*ezz + A*(eyy+exx)
    !    szz(:,:,:,:,:,:) = L + 2.0*mu_tau*ezz + two_mu*Mdz + C*(eyy+exx)

    !    !ELSE          !- anisotropy off

    !    !    sxx(:,:,:,:,:,:) =L+2.0*mu_tau*exx+2.0*mu*Mdx
    !    !    syy(:,:,:,:,:,:) =L+2.0*mu_tau*eyy+2.0*mu*Mdy
    !    !    szz(:,:,:,:,:,:) =L+2.0*mu_tau*ezz+2.0*mu*Mdz

    !    !ENDIF

    !!- off-diagonal stress rates -------------------------------------------------

    !    !Mdx=0.0
    !    !Mdy=0.0
    !    !Mdz=0.0
    !
    !    !DO k=1,nrdiss

    !    !    Mdx=Mdx+Mxy(k,:,:,:,:,:,:)
    !    !    Mdy=Mdy+Mxz(k,:,:,:,:,:,:)
    !    !    Mdz=Mdz+Myz(k,:,:,:,:,:,:)

    !    !ENDDO

    !    Mdx(:,:,:,:,:,:) = SUM( Mxx, dim=1 )
    !    Mdy(:,:,:,:,:,:) = SUM( Myy, dim=1 )
    !    Mdz(:,:,:,:,:,:) = SUM( Mzz, dim=1 )

    !    !IF (is_aniso==1) THEN        !- anisotropy on

    !    sxy(:,:,:,:,:,:) = 2.0*mu_tau*exy + two_mu*Mdx
    !    sxz(:,:,:,:,:,:) = 2.0*mu_tau*exz + two_mu*Mdy + 2.0*B*exz
    !    syz(:,:,:,:,:,:) = 2.0*mu_tau*eyz + two_mu*Mdz + 2.0*B*eyz

    !    !ELSE                !- anisotropy off

    !    !    sxy(:,:,:,:,:,:)=2.0*mu_tau*exy+2.0*mu*Mdx
    !    !    sxz(:,:,:,:,:,:)=2.0*mu_tau*exz+2.0*mu*Mdy
    !    !    syz(:,:,:,:,:,:)=2.0*mu_tau*eyz+2.0*mu*Mdz

    !    !ENDIF

    !ELSE            !- dissipation off

        ! Just for optimization
        L(:,:,:,:,:,:) = lambda * ( exx + eyy + ezz )

        ! Diagonal stress rates ----------------------------------------------
        sxx(:,:,:,:,:,:) = L + two_mu*exx + C*ezz + A*(eyy+exx)
        syy(:,:,:,:,:,:) = L + two_mu*eyy + C*ezz + A*(eyy+exx)
        szz(:,:,:,:,:,:) = L + two_mu*ezz + C*(eyy+exx)
        ! Off-diagonal stress rates ------------------------------------------
        sxy(:,:,:,:,:,:) = two_mu*exy
        sxz(:,:,:,:,:,:) = two_mu*exz + 2.0*B*exz
        syz(:,:,:,:,:,:) = two_mu*eyz + 2.0*B*eyz


    !ENDIF


    !======================================================================
    ! march pml stresses (integrate stress rates)
    !======================================================================

    sxx_pml(:,:,:,:,:,:)=sxx_pml+dt*(sxx)!-ispml*prof_x*sxx_pml) src_xx+
    sxy_pml(:,:,:,:,:,:)=sxy_pml+dt*(sxy)!-ispml*prof_x*sxy_pml) src_xy+
    sxz_pml(:,:,:,:,:,:)=sxz_pml+dt*(sxz)!-ispml*prof_x*sxz_pml) src_xz+
    syy_pml(:,:,:,:,:,:)=syy_pml+dt*(syy)!-ispml*prof_y*syy_pml) src_yy+
    syz_pml(:,:,:,:,:,:)=syz_pml+dt*(syz)!-ispml*prof_y*syz_pml) src_yz+
    szz_pml(:,:,:,:,:,:)=szz_pml+dt*(szz)!-ispml*prof_z*szz_pml) src_zz+


    ! add moment-tensor force
    sxx(:,:,:,:,:,:) = sxx_pml + src_xx
    sxy(:,:,:,:,:,:) = sxy_pml + src_xy
    sxz(:,:,:,:,:,:) = sxz_pml + src_xz
    syy(:,:,:,:,:,:) = syy_pml + src_yy
    syz(:,:,:,:,:,:) = syz_pml + src_yz
    szz(:,:,:,:,:,:) = szz_pml + src_zz

    ! add single force
    sx(:,:,:,:,:,:) = src_x
    sy(:,:,:,:,:,:) = src_y
    sz(:,:,:,:,:,:) = src_z

    !======================================================================
    ! make weak form of stress divergence
    !======================================================================

    !- compute stress divergence (with moment tensor added to the stress tensor)
    DO i=0,lpd
        DO j=0,lpd
            DO k=0,lpd

                ! XXX constant for all timesteps
                dummy_1 =  (2.0/dx_elm)*Jac
                dummy_2 =  (2.0/dy_elm)*Jac
                dummy_3 = -(2.0/dz_elm)*Jac

                DO q=0,lpd

                    ! XXX constant for all timesteps
                    dummy_x(:,:,:)=wx_wy_wz(:,:,:,q,j,k)*dl(i,q)*r_sin_theta(:,:,:,q,j,k)*dummy_1
                    dummy_y(:,:,:)=wx_wy_wz(:,:,:,i,q,k)*dl(j,q)*r(:,:,:,i,q,k)*dummy_2
                    dummy_z(:,:,:)=wx_wy_wz(:,:,:,i,j,q)*dl(k,q)*rr_sin_theta(:,:,:,i,j,q)*dummy_3

                    sx(:,:,:,i,j,k)=sx(:,:,:,i,j,k) &
                        +(sxz(:,:,:,i,j,q))*dummy_z &
                        +(sxy(:,:,:,i,q,k))*dummy_y &
                        +(sxx(:,:,:,q,j,k))*dummy_x

                    sy(:,:,:,i,j,k)=sy(:,:,:,i,j,k) &
                        +(syz(:,:,:,i,j,q))*dummy_z &
                        +(syy(:,:,:,i,q,k))*dummy_y &
                        +(sxy(:,:,:,q,j,k))*dummy_x

                    sz(:,:,:,i,j,k)=sz(:,:,:,i,j,k) &
                        +(szz(:,:,:,i,j,q))*dummy_z &
                        +(syz(:,:,:,i,q,k))*dummy_y &
                        +(sxz(:,:,:,q,j,k))*dummy_x

                ENDDO

            ENDDO
        ENDDO
    ENDDO

    !======================================================================
    ! make accelerations / divide weak stresses by mass matrix
    !======================================================================
    sx(:,:,:,:,:,:) = sx / MM
    sy(:,:,:,:,:,:) = sy / MM
    sz(:,:,:,:,:,:) = sz / MM

    !======================================================================
    ! communicate boundaries inside local process rank
    !======================================================================
    call communicate_local_boundaries( sx )
    call communicate_local_boundaries( sy )
    call communicate_local_boundaries( sz )

    !======================================================================
    ! communication boundaries amongst processors
    !======================================================================
    call communicate_global_boundaries( sx )
    call communicate_global_boundaries( sy )
    call communicate_global_boundaries( sz )

    !======================================================================
    ! time extrapolation of the displacement field
    !======================================================================
    vx(:,:,:,:,:,:) = vx + dt * (sx)!-ispml*prof*vx)
    vy(:,:,:,:,:,:) = vy + dt * (sy)!-ispml*prof*vy)
    vz(:,:,:,:,:,:) = vz + dt * (sz)!-ispml*prof*vz)

    !=======================================================================
    ! taper boundaries when pmls are turned off
    !=======================================================================
    vx(:,:,:,:,:,:) = vx * taper
    vy(:,:,:,:,:,:) = vy * taper
    vz(:,:,:,:,:,:) = vz * taper

    sxx_pml(:,:,:,:,:,:) = sxx_pml * taper
    sxy_pml(:,:,:,:,:,:) = sxy_pml * taper
    sxz_pml(:,:,:,:,:,:) = sxz_pml * taper
    syy_pml(:,:,:,:,:,:) = syy_pml * taper
    syz_pml(:,:,:,:,:,:) = syz_pml * taper
    szz_pml(:,:,:,:,:,:) = szz_pml * taper

    !======================================================================
    ! time extrapolation of memory variables
    !======================================================================

    !IF (is_diss==1) THEN

    !    DO k=1,nrdiss
    !
    !        ! this is constant for all timesteps
    !        L(:,:,:,:,:,:) = dt * tau / ( nrdiss * tau_p(k) )

    !        Mxx(k,:,:,:,:,:,:)=Mxx(k,:,:,:,:,:,:)-dt*Mxx(k,:,:,:,:,:,:)/tau_p(k) &
    !            -L(:,:,:,:,:,:)*exx(:,:,:,:,:,:)

    !        Myy(k,:,:,:,:,:,:)=Myy(k,:,:,:,:,:,:)-dt*Myy(k,:,:,:,:,:,:)/tau_p(k) &
    !            -L(:,:,:,:,:,:)*eyy(:,:,:,:,:,:)

    !        Mzz(k,:,:,:,:,:,:)=Mzz(k,:,:,:,:,:,:)-dt*Mzz(k,:,:,:,:,:,:)/tau_p(k) &
    !            -L(:,:,:,:,:,:)*ezz(:,:,:,:,:,:)
    !
    !        Mxy(k,:,:,:,:,:,:)=Mxy(k,:,:,:,:,:,:)-dt*Mxy(k,:,:,:,:,:,:)/tau_p(k) &
    !            -L(:,:,:,:,:,:)*exy(:,:,:,:,:,:)

    !        Mxz(k,:,:,:,:,:,:)=Mxz(k,:,:,:,:,:,:)-dt*Mxz(k,:,:,:,:,:,:)/tau_p(k) &
    !            -L(:,:,:,:,:,:)*exz(:,:,:,:,:,:)
    !
    !        Myz(k,:,:,:,:,:,:)=Myz(k,:,:,:,:,:,:)-dt*Myz(k,:,:,:,:,:,:)/tau_p(k) &
    !            -L(:,:,:,:,:,:)*eyz(:,:,:,:,:,:)

    !    ENDDO

    !ENDIF

END SUBROUTINE ses3d_evolution

END MODULE evolution_mod
