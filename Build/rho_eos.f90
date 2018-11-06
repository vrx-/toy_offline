      MODULE rho_eos_mod
!
!svn $Id: rho_eos.F 854 2017-07-18 23:28:45Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes  "in situ" density and other associated       !
!  quantitites as a function of potential temperature,  salinity,      !
!  and pressure from a polynomial expression (Jackett & McDougall,     !
!  1992). The polynomial expression was found from fitting to 248      !
!  values  in the  oceanographic  ranges of  salinity,  potential      !
!  temperature,  and pressure.  It  assumes no pressure variation      !
!  along geopotential surfaces, that is, depth (meters; negative)      !
!  and pressure (dbar; assumed negative here) are interchangeable.     !
!                                                                      !
!  Check Values: (T=3 C, S=35.5 PSU, Z=-5000 m)                        !
!                                                                      !
!     alpha = 2.1014611551470d-04 (1/Celsius)                          !
!     beta  = 7.2575037309946d-04 (1/PSU)                              !
!     gamma = 3.9684764511766d-06 (1/Pa)                               !
!     den   = 1050.3639165364     (kg/m3)                              !
!     den1  = 1028.2845117925     (kg/m3)                              !
!     sound = 1548.8815240223     (m/s)                                !
!     bulk  = 23786.056026320     (Pa)                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of     !
!    Hydrostatic Profiles to Achieve Static Stability, J. of Atmos.    !
!    and Oceanic Techn., vol. 12, pp. 381-389.                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: rho_eos
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rho_eos (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL wclock_on (ng, model, 14, 66, "ROMS/Nonlinear/rho_eos.F")
      CALL rho_eos_tile (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   OCEAN(ng) % t,                                 &
     &                   COUPLING(ng) % rhoA,                           &
     &                   COUPLING(ng) % rhoS,                           &
     &                   OCEAN(ng) % pden,                              &
     &                   OCEAN(ng) % rho)
      CALL wclock_off (ng, model, 14, 99, "ROMS/Nonlinear/rho_eos.F")
      RETURN
      END SUBROUTINE rho_eos
!
!***********************************************************************
      SUBROUTINE rho_eos_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
     &                         Hz,                                      &
     &                         z_r, z_w, t,                             &
     &                         rhoA, rhoS,                              &
     &                         pden,                                    &
     &                         rho)
!***********************************************************************
!
      USE mod_param
      USE mod_eoscoef
      USE mod_scalars
!
      USE exchange_2d_mod
      USE exchange_3d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: rhoA(LBi:,LBj:)
      real(r8), intent(out) :: rhoS(LBi:,LBj:)
      real(r8), intent(out) :: pden(LBi:,LBj:,:)
      real(r8), intent(out) :: rho(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, ised, itrc, j, k
      real(r8) :: SedDen, Tp, Tpr10, Ts, Tt, sqrtTs
      real(r8) :: cff, cff1, cff2
      real(r8), dimension(0:9) :: C
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk0
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den1
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!=======================================================================
!  Nonlinear equation of state.  Notice that this equation of state
!  is only valid for potential temperature range of -2C to 40C and
!  a salinity range of 0 PSU to 42 PSU.
!  Because of the pressure effect on freezing and the iceshelf, we've
!  extended temp limit to -2.5 (msd 4/25/06)
!=======================================================================
!
      DO j=JstrT,JendT
        DO k=1,N(ng)
          DO i=IstrT,IendT
!
!  Check temperature and salinity lower values. Assign depth to the
!  pressure.
!
            Tt=MAX(-2.5_r8,t(i,j,k,nrhs,itemp))
            Tt=MIN(40.0_r8,Tt)
            Ts=0.0_r8
            sqrtTs=0.0_r8
            Tp=z_r(i,j,k)
            Tpr10=0.1_r8*Tp
!
!-----------------------------------------------------------------------
!  Compute density (kg/m3) at standard one atmosphere pressure.
!-----------------------------------------------------------------------
!
            C(0)=Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
            C(1)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
            C(2)=V00+Tt*(V01+Tt*V02)
!
            den1(i,k)=C(0)+Ts*(C(1)+sqrtTs*C(2)+Ts*W00)
!
!-----------------------------------------------------------------------
!  Compute secant bulk modulus.
!-----------------------------------------------------------------------
!
            C(3)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
            C(4)=B00+Tt*(B01+Tt*(B02+Tt*B03))
            C(5)=D00+Tt*(D01+Tt*D02)
            C(6)=E00+Tt*(E01+Tt*(E02+Tt*E03))
            C(7)=F00+Tt*(F01+Tt*F02)
            C(8)=G01+Tt*(G02+Tt*G03)
            C(9)=H00+Tt*(H01+Tt*H02)
!
            bulk0(i,k)=C(3)+Ts*(C(4)+sqrtTs*C(5))
            bulk1(i,k)=C(6)+Ts*(C(7)+sqrtTs*G00)
            bulk2(i,k)=C(8)+Ts*C(9)
            bulk (i,k)=bulk0(i,k)-Tp*(bulk1(i,k)-Tp*bulk2(i,k))
!
!-----------------------------------------------------------------------
!  Compute local "in situ" density anomaly (kg/m3 - 1000).
!-----------------------------------------------------------------------
!
            cff=1.0_r8/(bulk(i,k)+Tpr10)
            den(i,k)=den1(i,k)*bulk(i,k)*cff
            den(i,k)=den(i,k)-1000.0_r8
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Compute vertical averaged density (rhoA) and density perturbation
!  (rhoS) used in barotropic pressure gradient.
!-----------------------------------------------------------------------
!
        DO i=IstrT,IendT
          cff1=den(i,N(ng))*Hz(i,j,N(ng))
          rhoS(i,j)=0.5_r8*cff1*Hz(i,j,N(ng))
          rhoA(i,j)=cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrT,IendT
            cff1=den(i,k)*Hz(i,j,k)
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5_r8*cff1)
            rhoA(i,j)=rhoA(i,j)+cff1
          END DO
        END DO
        cff2=1.0_r8/rho0
        DO i=IstrT,IendT
          cff1=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          rhoA(i,j)=cff2*cff1*rhoA(i,j)
          rhoS(i,j)=2.0_r8*cff1*cff1*cff2*rhoS(i,j)
        END DO
!
!-----------------------------------------------------------------------
!  Load "in situ" density anomaly (kg/m3 - 1000) and potential
!  density anomaly (kg/m3 - 1000) referenced to the surface into global
!  arrays. Notice that this is done in a separate (i,k) DO-loops to
!  facilitate the adjoint.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)
          DO i=IstrT,IendT
            rho(i,j,k)=den(i,k)
            pden(i,j,k)=(den1(i,k)-1000.0_r8)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rho)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          pden)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoA)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoS)
      END IF
      RETURN
      END SUBROUTINE rho_eos_tile
      END MODULE rho_eos_mod
