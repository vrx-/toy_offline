      MODULE rhs3d_mod
!
!svn $Id: rhs3d.F 854 2017-07-18 23:28:45Z arango $
!=======================================================================
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine evaluates right-hand-side terms for 3D momentum     !
!  and tracers equations.                                              !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: rhs3d
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rhs3d (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
      USE pre_step3d_mod, ONLY : pre_step3d
      USE prsgrd_mod, ONLY : prsgrd
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
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
!-----------------------------------------------------------------------
!  Initialize computations for new time step of the 3D primitive
!  variables.
!-----------------------------------------------------------------------
!
      CALL pre_step3d (ng, tile)
!
!-----------------------------------------------------------------------
!  Compute baroclinic pressure gradient.
!-----------------------------------------------------------------------
!
      CALL prsgrd (ng, tile)
!
!
      RETURN
      END SUBROUTINE rhs3d
!
!***********************************************************************
      SUBROUTINE rhs3d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nrhs,                                      &
     &                       Hz, Huon, Hvom,                            &
     &                       fomn,                                      &
     &                       om_u, om_v, on_u, on_v, pm, pn,            &
     &                       bustr, bvstr,                              &
     &                       sustr, svstr,                              &
     &                       u, v, W,                                   &
     &                       rufrc, rvfrc,                              &
     &                       ru, rv)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: fomn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
      real(r8), intent(out) :: rufrc(LBi:,LBj:)
      real(r8), intent(out) :: rvfrc(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8), parameter :: Gadv = -0.25_r8
      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: cff5, cff6, cff7, cff8
      real(r8) :: fac, fac1, fac2
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
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
!
      K_LOOP : DO k=1,N(ng)
!
!-----------------------------------------------------------------------
!  Add in nudging of 3D momentum climatology.
!-----------------------------------------------------------------------
!
        IF (LnudgeM3CLM(ng)) THEN
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=0.25_r8*(CLIMA(ng)%M3nudgcof(i-1,j,k)+                &
     &                     CLIMA(ng)%M3nudgcof(i  ,j,k))*               &
     &            om_u(i,j)*on_u(i,j)
              ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+                            &
     &                       cff*(Hz(i-1,j,k)+Hz(i,j,k))*               &
     &                       (CLIMA(ng)%uclm(i,j,k)-                    &
     &                        u(i,j,k,nrhs))
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=0.25_r8*(CLIMA(ng)%M3nudgcof(i,j-1,k)+                &
     &                     CLIMA(ng)%M3nudgcof(i,j  ,k))*               &
     &            om_v(i,j)*on_v(i,j)
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)+                            &
     &                       cff*(Hz(i,j-1,k)+Hz(i,j,k))*               &
     &                       (CLIMA(ng)%vclm(i,j,k)-                    &
     &                        v(i,j,k,nrhs))
            END DO
          END DO
        END IF
      END DO K_LOOP
!
      J_LOOP : DO j=Jstr,Jend
!
!-----------------------------------------------------------------------
!  Compute forcing term for the 2D momentum equations.
!-----------------------------------------------------------------------
!
!  Vertically integrate baroclinic right-hand-side terms. If not
!  body force stresses, add in the difference between surface and
!  bottom stresses.
!
        DO i=IstrU,Iend
          rufrc(i,j)=ru(i,j,1,nrhs)
        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend
            rufrc(i,j)=rufrc(i,j)+ru(i,j,k,nrhs)
          END DO
        END DO
        DO i=IstrU,Iend
          cff=om_u(i,j)*on_u(i,j)
          cff1= sustr(i,j)*cff
          cff2=-bustr(i,j)*cff
          rufrc(i,j)=rufrc(i,j)+cff1+cff2
        END DO
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            rvfrc(i,j)=rv(i,j,1,nrhs)
          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend
              rvfrc(i,j)=rvfrc(i,j)+rv(i,j,k,nrhs)
            END DO
          END DO
          DO i=Istr,Iend
            cff=om_v(i,j)*on_v(i,j)
            cff1= svstr(i,j)*cff
            cff2=-bvstr(i,j)*cff
            rvfrc(i,j)=rvfrc(i,j)+cff1+cff2
          END DO
        END IF
      END DO J_LOOP
      RETURN
      END SUBROUTINE rhs3d_tile
      END MODULE rhs3d_mod
