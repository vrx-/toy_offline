      MODULE pt3dbc_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for the ITRC-th    !
!  tracer field.                                                       !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: pt3dbc_tile
      CONTAINS
!
!***********************************************************************
      SUBROUTINE pt3dbc (ng, tile, nout, itrc)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nout, itrc
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
      CALL pt3dbc_tile (ng, tile, itrc,                                 &
     &                  LBi, UBi, LBj, UBj, N(ng), NT(ng),              &
     &                  IminS, ImaxS, JminS, JmaxS,                      &
     &                  nstp(ng), nout,                                 &
     &                  OCEAN(ng)% t)
      RETURN
      END SUBROUTINE pt3dbc
!
!***********************************************************************
      SUBROUTINE pt3dbc_tile (ng, tile, itrc,                           &
     &                        LBi, UBi, LBj, UBj, UBk, UBt,             &
     &                        IminS, ImaxS, JminS, JmaxS,                &
     &                        nstp, nout,                               &
     &                       t)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, itrc
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nout
!
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8), parameter :: eps =1.0E-20_r8
      real(r8) :: Ce, Cx, cff, dTde, dTdt, dTdx, tau
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
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
      IF (.not.EWperiodic(ng)) THEN
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, closed boundary condition.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend
              t(Istr-1,j,k,nout,itrc)=t(Istr,j,k,nout,itrc)
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, closed boundary condition.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend
              t(Iend+1,j,k,nout,itrc)=t(Iend,j,k,nout,itrc)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, closed boundary condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,Jstr-1,k,nout,itrc)=t(i,Jstr,k,nout,itrc)
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, closed boundary condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,Jend+1,k,nout,itrc)=t(i,Jend,k,nout,itrc)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng).and. .not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          DO k=1,N(ng)
            t(Istr-1,Jstr-1,k,nout,itrc)=0.5_r8*                        &
     &                                   (t(Istr  ,Jstr-1,k,nout,itrc)+ &
     &                                    t(Istr-1,Jstr  ,k,nout,itrc))
          END DO
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          DO k=1,N(ng)
            t(Iend+1,Jstr-1,k,nout,itrc)=0.5_r8*                        &
     &                                   (t(Iend  ,Jstr-1,k,nout,itrc)+ &
     &                                    t(Iend+1,Jstr  ,k,nout,itrc))
          END DO
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          DO k=1,N(ng)
            t(Istr-1,Jend+1,k,nout,itrc)=0.5_r8*                        &
     &                                   (t(Istr-1,Jend  ,k,nout,itrc)+ &
     &                                    t(Istr  ,Jend+1,k,nout,itrc))
          END DO
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          DO k=1,N(ng)
            t(Iend+1,Jend+1,k,nout,itrc)=0.5_r8*                        &
     &                                   (t(Iend+1,Jend  ,k,nout,itrc)+ &
     &                                    t(Iend  ,Jend+1,k,nout,itrc))
          END DO
        END IF
      END IF
      RETURN
      END SUBROUTINE pt3dbc_tile
      END MODULE pt3dbc_mod
