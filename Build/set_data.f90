      SUBROUTINE set_data (ng, tile)
!
!svn $Id: set_data.F 858 2017-07-31 23:02:30Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine processes forcing, boundary, climatology, and       !
!  other input data. It time-interpolates between snapshots.           !
!                                                                      !
!=======================================================================
!
      USE mod_param
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
      CALL wclock_on (ng, iNLM, 4, 28, "ROMS/Nonlinear/set_data.F")
      CALL set_data_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS)
      CALL wclock_off (ng, iNLM, 4, 34, "ROMS/Nonlinear/set_data.F")
      RETURN
      END SUBROUTINE set_data
!
!***********************************************************************
      SUBROUTINE set_data_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
      USE mod_scalars
      USE mod_sources
!
      USE analytical_mod
      USE exchange_2d_mod
      USE set_2dfld_mod
      USE set_3dfld_mod
      USE strings_mod,     ONLY : FoundError
      USE mod_stepping
      USE exchange_2d_mod
      USE exchange_3d_mod
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      logical :: SetBC
      logical :: update = .FALSE.
      integer :: ILB, IUB, JLB, JUB
      integer :: i, ic, itrc, j, k, my_tile
      real(r8) :: cff, cff1, cff2
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
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      my_tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(my_tile)
      IUB=BOUNDS(ng)%UBi(my_tile)
      JLB=BOUNDS(ng)%LBj(my_tile)
      JUB=BOUNDS(ng)%UBj(my_tile)
!
!-----------------------------------------------------------------------
!  Set kinematic surface solar shortwave radiation flux (degC m/s).
!-----------------------------------------------------------------------
!
      CALL set_2dfld_tile (ng, tile, iNLM, idSrad,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%srflxG,                           &
     &                     FORCES(ng)%srflx,                            &
     &                     update)
      IF (FoundError(exit_flag, NoError, 378,                           &
     &               "ROMS/Nonlinear/set_data.F")) RETURN
!
!  Modulate the averaged shortwave radiation flux by the local diurnal
!  cycle.
!
      CALL ana_srflux (ng, tile, iNLM)
!
!-----------------------------------------------------------------------
!  Set kinematic surface net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
      CALL set_2dfld_tile (ng, tile, iNLM, idTsur(itemp),               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%stflxG(:,:,:,itemp),              &
     &                     FORCES(ng)%stflx (:,:,itemp),                &
     &                     update)
      IF (FoundError(exit_flag, NoError, 569,                           &
     &               "ROMS/Nonlinear/set_data.F")) RETURN
!
!-----------------------------------------------------------------------
!  Set kinematic bottom net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
      CALL set_2dfld_tile (ng, tile, iNLM, idTbot(itemp),               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%btflxG(:,:,:,itemp),              &
     &                     FORCES(ng)%btflx (:,:,itemp),                &
     &                     update)
      IF (FoundError(exit_flag, NoError, 636,                           &
     &               "ROMS/Nonlinear/set_data.F")) RETURN
!
!-----------------------------------------------------------------------
!  Set kinematic surface and bottom pasive tracer fluxes (T m/s).
!-----------------------------------------------------------------------
!
      DO itrc=NAT+1,NT(ng)
        CALL ana_stflux (ng, tile, iNLM, itrc)
        CALL ana_btflux (ng, tile, iNLM, itrc)
      END DO
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (m2/s2).
!-----------------------------------------------------------------------
!
      CALL ana_smflux (ng, tile, iNLM)
!
!-----------------------------------------------------------------------
!  Set point Sources/Sinks (river runoff).
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
        IF (LuvSrc(ng).or.LwSrc(ng)) THEN
          CALL set_ngfld (ng, iNLM, idRtra, 1, Nsrc(ng), 1,             &
     &                    1, Nsrc(ng), 1,                               &
     &                    SOURCES(ng) % QbarG,                          &
     &                    SOURCES(ng) % Qbar,                           &
     &                    update)
          IF (FoundError(exit_flag, NoError, 1057,                      &
     &                   "ROMS/Nonlinear/set_data.F")) RETURN
          DO k=1,N(ng)
            DO i=1,Nsrc(ng)
              SOURCES(ng)%Qsrc(i,k)=SOURCES(ng)%Qbar(i)*                &
     &                              SOURCES(ng)%Qshape(i,k)
            END DO
          END DO
        END IF
        DO itrc=1,NT(ng)
          IF (LtracerSrc(itrc,ng)) THEN
            CALL set_ngfld (ng, iNLM, idRtrc(itrc), 1, Nsrc(ng), N(ng), &
     &                      1, Nsrc(ng), N(ng),                         &
     &                      SOURCES(ng) % TsrcG(:,:,:,itrc),            &
     &                      SOURCES(ng) % Tsrc(:,:,itrc),               &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1092,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Set open boundary conditions fields.
!-----------------------------------------------------------------------
!
!  Free-surface
!
      IF (LprocessOBC(ng)) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isFsur,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idZbry(iwest), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % zetaG_west,                  &
     &                      BOUNDARY(ng) % zeta_west,                   &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1151,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(ieast,isFsur,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idZbry(ieast), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % zetaG_east,                  &
     &                      BOUNDARY(ng) % zeta_east,                   &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1161,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(isouth,isFsur,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idZbry(isouth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1 ,1,                             &
     &                      BOUNDARY(ng) % zetaG_south,                 &
     &                      BOUNDARY(ng) % zeta_south,                  &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1171,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(inorth,isFsur,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idZbry(inorth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % zetaG_north,                 &
     &                      BOUNDARY(ng) % zeta_north,                  &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1181,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
        END IF
      END IF
!
!  2D momentum.
!
      IF (LprocessOBC(ng)) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isUbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU2bc(iwest), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % ubarG_west,                  &
     &                      BOUNDARY(ng) % ubar_west,                   &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1278,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(iwest,isVbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV2bc(iwest), JLB, JUB, 1,       &
     &                      1, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % vbarG_west,                  &
     &                      BOUNDARY(ng) % vbar_west,                   &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1288,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(ieast,isUbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU2bc(ieast), JLB, JUB, 1,       &
     &                      0, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % ubarG_east,                  &
     &                      BOUNDARY(ng) % ubar_east,                   &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1298,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(ieast,isVbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV2bc(ieast), JLB, JUB, 1,       &
     &                      1, Mm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % vbarG_east,                  &
     &                      BOUNDARY(ng) % vbar_east,                   &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1308,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(isouth,isUbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU2bc(isouth), ILB, IUB, 1,      &
     &                      1, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % ubarG_south,                 &
     &                      BOUNDARY(ng) % ubar_south,                  &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1318,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(isouth,isVbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV2bc(isouth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % vbarG_south,                 &
     &                      BOUNDARY(ng) % vbar_south,                  &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1328,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(inorth,isUbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU2bc(inorth), ILB, IUB, 1,      &
     &                      1, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % ubarG_north,                 &
     &                      BOUNDARY(ng) % ubar_north,                  &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1338,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(inorth,isVbar,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV2bc(inorth), ILB, IUB, 1,      &
     &                      0, Lm(ng)+1, 1,                             &
     &                      BOUNDARY(ng) % vbarG_north,                 &
     &                      BOUNDARY(ng) % vbar_north,                  &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1348,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
        END IF
      END IF
!
!  3D momentum.
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (LBC(iwest,isUvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU3bc(iwest), JLB, JUB, N(ng),   &
     &                      0, Mm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % uG_west,                     &
     &                      BOUNDARY(ng) % u_west,                      &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1862,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(iwest,isVvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV3bc(iwest), JLB, JUB, N(ng),   &
     &                      1, Mm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % vG_west,                     &
     &                      BOUNDARY(ng) % v_west,                      &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1872,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(ieast,isUvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU3bc(ieast), JLB, JUB, N(ng),   &
     &                      0, Mm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % uG_east,                     &
     &                      BOUNDARY(ng) % u_east,                      &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1882,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(ieast,isVvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV3bc(ieast), JLB, JUB, N(ng),   &
     &                      1, Mm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % vG_east,                     &
     &                      BOUNDARY(ng) % v_east,                      &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1892,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(isouth,isUvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU3bc(isouth), ILB, IUB, N(ng),  &
     &                      1, Lm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % uG_south,                    &
     &                      BOUNDARY(ng) % u_south,                     &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1902,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(isouth,isVvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV3bc(isouth), ILB, IUB, N(ng),  &
     &                      0, Lm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % vG_south,                    &
     &                      BOUNDARY(ng) % v_south,                     &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1912,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(inorth,isUvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idU3bc(inorth), ILB, IUB, N(ng),  &
     &                      1, Lm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % uG_north,                    &
     &                      BOUNDARY(ng) % u_north,                     &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1922,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
          IF (LBC(inorth,isVvel,ng)%acquire) THEN
            CALL set_ngfld (ng, iNLM, idV3bc(inorth), ILB, IUB, N(ng),  &
     &                      0, Lm(ng)+1, N(ng),                         &
     &                      BOUNDARY(ng) % vG_north,                    &
     &                      BOUNDARY(ng) % v_north,                     &
     &                      update)
            IF (FoundError(exit_flag, NoError, 1932,                    &
     &                     "ROMS/Nonlinear/set_data.F")) RETURN
          END IF
        END IF
      END IF
!
!  Tracers.
!
      IF (LprocessOBC(ng)) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
           DO itrc=1,NT(ng)
            IF (LBC(iwest,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfld (ng, iNLM, idTbry(iwest,itrc),             &
     &                        JLB, JUB, N(ng), 0, Mm(ng)+1, N(ng),      &
     &                        BOUNDARY(ng) % tG_west(:,:,:,itrc),       &
     &                        BOUNDARY(ng) % t_west(:,:,itrc),          &
     &                        update)
              IF (FoundError(exit_flag, NoError, 1953,                  &
     &                       "ROMS/Nonlinear/set_data.F")) RETURN
            END IF
            IF (LBC(ieast,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfld (ng, iNLM, idTbry(ieast,itrc),             &
     &                        JLB, JUB, N(ng), 0, Mm(ng)+1, N(ng),      &
     &                        BOUNDARY(ng) % tG_east(:,:,:,itrc),       &
     &                        BOUNDARY(ng) % t_east(:,:,itrc),          &
     &                        update)
              IF (FoundError(exit_flag, NoError, 1963,                  &
     &                       "ROMS/Nonlinear/set_data.F")) RETURN
            END IF
            IF (LBC(isouth,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfld (ng, iNLM, idTbry(isouth,itrc),            &
     &                        ILB, IUB, N(ng), 0, Lm(ng)+1, N(ng),      &
     &                        BOUNDARY(ng) % tG_south(:,:,:,itrc),      &
     &                        BOUNDARY(ng) % t_south(:,:,itrc),         &
     &                        update)
              IF (FoundError(exit_flag, NoError, 1973,                  &
     &                       "ROMS/Nonlinear/set_data.F")) RETURN
            END IF
            IF (LBC(inorth,isTvar(itrc),ng)%acquire) THEN
              CALL set_ngfld (ng, iNLM, idTbry(inorth,itrc),            &
     &                        ILB, IUB, N(ng), 0, Lm(ng)+1, N(ng),      &
     &                        BOUNDARY(ng) % tG_north(:,:,:,itrc),      &
     &                        BOUNDARY(ng) % t_north(:,:,itrc),         &
     &                        update)
              IF (FoundError(exit_flag, NoError, 1983,                  &
     &                       "ROMS/Nonlinear/set_data.F")) RETURN
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Set sea surface height climatology (m).
!-----------------------------------------------------------------------
!
      IF (LsshCLM(ng)) THEN
        CALL set_2dfld_tile (ng, tile, iNLM, idSSHc,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%sshG,                            &
     &                       CLIMA(ng)%ssh,                             &
     &                       update)
        IF (FoundError(exit_flag, NoError, 2005,                        &
     &                 "ROMS/Nonlinear/set_data.F")) RETURN
        DO j=Jstr,Jend
          DO i=Istr,Iend
            OCEAN(ng)%zeta(i,j,1) = CLIMA(ng)%ssh(i,j)
            OCEAN(ng)%zeta(i,j,2) = CLIMA(ng)%ssh(i,j)
          END DO
        END DO
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          DO i=1,2
            CALL exchange_r2d_tile (ng, tile,                           &
     &                      LBi, UBi, LBj, UBj, OCEAN(ng)%zeta(:,:,i))
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Set 2D momentum climatology (m/s).
!-----------------------------------------------------------------------
!
      IF (Lm2CLM(ng)) THEN
        CALL set_2dfld_tile (ng, tile, iNLM, idUbcl,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%ubarclmG,                        &
     &                       CLIMA(ng)%ubarclm,                         &
     &                       update)
        IF (FoundError(exit_flag, NoError, 2044,                        &
     &                 "ROMS/Nonlinear/set_data.F")) RETURN
!
        CALL set_2dfld_tile (ng, tile, iNLM, idVbcl,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%vbarclmG,                        &
     &                       CLIMA(ng)%vbarclm,                         &
     &                       update)
        IF (FoundError(exit_flag, NoError, 2052,                        &
     &                 "ROMS/Nonlinear/set_data.F")) RETURN
        DO j=Jstr,Jend
          DO i=Istr,Iend
            OCEAN(ng)%ubar(i,j,1) = CLIMA(ng)%ubarclm(i,j)
            OCEAN(ng)%ubar(i,j,2) = CLIMA(ng)%ubarclm(i,j)
            OCEAN(ng)%vbar(i,j,1) = CLIMA(ng)%vbarclm(i,j)
            OCEAN(ng)%vbar(i,j,2) = CLIMA(ng)%vbarclm(i,j)
          ENDDO
        ENDDO
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          DO i=1,2
            CALL exchange_u2d_tile (ng, tile,                           &
     &                       LBi, UBi, LBj, UBj, OCEAN(ng)%ubar(:,:,i))
            CALL exchange_v2d_tile (ng, tile,                           &
     &                       LBi, UBi, LBj, UBj, OCEAN(ng)%vbar(:,:,i))
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Set 3D momentum climatology (m/s).
!-----------------------------------------------------------------------
!
      IF (Lm3CLM(ng)) THEN
        CALL set_3dfld_tile (ng, tile, iNLM, idUclm,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       CLIMA(ng)%uclmG,                           &
     &                       CLIMA(ng)%uclm,                            &
     &                       update)
        IF (FoundError(exit_flag, NoError, 2097,                        &
     &                 "ROMS/Nonlinear/set_data.F")) RETURN
!
        CALL set_3dfld_tile (ng, tile, iNLM, idVclm,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       CLIMA(ng)%vclmG,                           &
     &                       CLIMA(ng)%vclm,                            &
     &                       update)
        IF (FoundError(exit_flag, NoError, 2105,                        &
     &                 "ROMS/Nonlinear/set_data.F")) RETURN
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend
              OCEAN(ng)%u(i,j,k,nnew(ng)) = CLIMA(ng)%uclm(i,j,k)
              OCEAN(ng)%v(i,j,k,nnew(ng)) = CLIMA(ng)%vclm(i,j,k)
            END DO
          END DO
        END DO
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          DO i=1,2
            CALL exchange_u3d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            OCEAN(ng)%u(:,:,:,i))
            CALL exchange_v3d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            OCEAN(ng)%v(:,:,:,i))
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Set tracers climatology.
!-----------------------------------------------------------------------
!
      ic=0
      DO itrc=1,NT(ng)
        IF (LtracerCLM(itrc,ng)) THEN
          ic=ic+1
          CALL set_3dfld_tile (ng, tile, iNLM, idTclm(itrc),            &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         CLIMA(ng)%tclmG(:,:,:,:,ic),             &
     &                         CLIMA(ng)%tclm (:,:,:,ic),               &
     &                         update)
          IF (FoundError(exit_flag, NoError, 2158,                      &
     &                   "ROMS/Nonlinear/set_data.F")) RETURN
        END IF
      END DO
      DO itrc=1,NAT
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend
              OCEAN(ng)%t(i,j,k,nnew(ng),itrc) =                        &
     &                  CLIMA(ng)%tclm(i,j,k,itrc)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NAT
          DO i=1,2
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              OCEAN(ng)%t(:,:,:,i,itrc))
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE set_data_tile
