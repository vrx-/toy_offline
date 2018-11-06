      SUBROUTINE wrt_info (ng, model, ncid, ncname)
!
!svn $Id: wrt_info.F 857 2017-07-29 04:05:27Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes out information variables into requested        !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer)                            !
!     model    Calling model identifier (integer)                      !
!     ncid     NetCDF file ID (integer)                                !
!     ncname   NetCDF file name (string)                               !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_grid
      Use mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_sources
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE strings_mod,     ONLY : FoundError, find_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      logical :: Cgrid = .TRUE.
      integer :: LBi, UBi, LBj, UBj
      integer :: i, j, k, ibry, ilev, itrc, status, varid
      integer :: ifield = 0
      real(r8) :: scale
      real(r8), dimension(NT(ng)) :: nudg
      real(r8), dimension(NT(ng),4) :: Tobc
!
      SourceFile="ROMS/Utility/wrt_info.F"
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out running parameters.
!-----------------------------------------------------------------------
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname, ncid)
      IF (FoundError(exit_flag, NoError, 114,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  Time stepping parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'ntimes',                &
     &                      ntimes(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 122,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ndtfast',               &
     &                      ndtfast(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 128,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'dt',                    &
     &                      dt(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 134,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'dtfast',                &
     &                      dtfast(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 140,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'dstart',                &
     &                      dstart, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 146,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'nHIS',                  &
     &                      nHIS(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 172,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ndefHIS',               &
     &                      ndefHIS(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 178,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'nRST',                  &
     &                      nRST(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 184,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ntsAVG',                &
     &                      ntsAVG(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 194,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'nAVG',                  &
     &                      nAVG(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 200,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ndefAVG',               &
     &                      ndefAVG(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 206,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ntsDIA',                &
     &                      ntsDIA(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 334,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'nDIA',                  &
     &                      nDIA(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 340,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ndefDIA',               &
     &                      ndefDIA(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 346,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Falpha',                &
     &                      Falpha, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 380,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Fbeta',                 &
     &                      Fbeta, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 386,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Fgamma',                &
     &                      Fgamma, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 392,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  Horizontal mixing coefficients.
!
!
!  Background vertical mixing coefficients.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Akt_bak',               &
     &                      Akt_bak(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 545,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Akv_bak',               &
     &                      Akv_bak(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 551,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  Drag coefficients.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'rdrg',                  &
     &                      rdrg(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 611,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'rdrg2',                 &
     &                      rdrg2(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 617,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Zob',                   &
     &                      Zob(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 624,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Zos',                   &
     &                      Zos(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 630,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  Nudging inverse time scales used in various tasks.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Znudg',                 &
     &                      Znudg(ng)/sec2day, (/0/), (/0/),            &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 747,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'M2nudg',                &
     &                      M2nudg(ng)/sec2day, (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 753,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'M3nudg',                &
     &                      M3nudg(ng)/sec2day, (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 760,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      DO itrc=1,NT(ng)
        nudg(itrc)=Tnudg(itrc,ng)/sec2day
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'Tnudg',                 &
     &                      nudg, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 769,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  Open boundary nudging, inverse time scales.
!
      IF (NudgingCoeff(ng)) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'FSobc_in',            &
     &                        FSobc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 781,                         &
     &                 "ROMS/Utility/wrt_info.F")) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'FSobc_out',           &
     &                        FSobc_out(ng,:), (/1/), (/4/),            &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 787,                         &
     &                 "ROMS/Utility/wrt_info.F")) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M2obc_in',            &
     &                        M2obc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 793,                         &
     &                 "ROMS/Utility/wrt_info.F")) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M2obc_out',           &
     &                        M2obc_out(ng,:), (/1/), (/4/),            &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 799,                         &
     &                 "ROMS/Utility/wrt_info.F")) RETURN
        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_in(itrc,ng,ibry)
          END DO
        END DO
        CALL netcdf_put_fvar (ng, model, ncname, 'Tobc_in',             &
     &                        Tobc, (/1,1/), (/NT(ng),4/),              &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 811,                         &
     &                 "ROMS/Utility/wrt_info.F")) RETURN
        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_out(itrc,ng,ibry)
          END DO
        END DO
        CALL netcdf_put_fvar (ng, model, ncname, 'Tobc_out',            &
     &                        Tobc, (/1,1/), (/NT(ng),4/),              &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 822,                         &
     &                 "ROMS/Utility/wrt_info.F")) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M3obc_in',            &
     &                        M3obc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 828,                         &
     &                 "ROMS/Utility/wrt_info.F")) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M3obc_out',           &
     &                        M3obc_out(ng,:), (/1/), (/4/),            &
     &                      ncid = ncid)
        IF (FoundError(exit_flag, NoError, 834,                         &
     &                 "ROMS/Utility/wrt_info.F")) RETURN
      END IF
!
!  Equation of State parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'rho0',                  &
     &                      rho0, (/0/), (/0/),                         &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 845,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  Slipperiness parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'gamma2',                &
     &                      gamma2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 900,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
! Logical switches to activate horizontal momentum transport
! point Sources/Sinks (like river runoff transport) and mass point
! Sources/Sinks (like volume vertical influx).
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LuvSrc',                &
     &                      LuvSrc(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 910,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LwSrc',                 &
     &                      LwSrc(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 916,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  Logical switches to activate tracer point Sources/Sinks.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LtracerSrc',            &
     &                      LtracerSrc(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 926,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  Logical switches to process climatology fields.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LsshCLM',               &
     &                      LsshCLM(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 935,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'Lm2CLM',                &
     &                      Lm2CLM(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 941,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'Lm3CLM',                &
     &                      Lm3CLM(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 948,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LtracerCLM',            &
     &                      LtracerCLM(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 954,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  Logical switches for nudging climatology fields.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LnudgeM2CLM',           &
     &                      LnudgeM2CLM(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 963,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LnudgeM3CLM',           &
     &                      LnudgeM3CLM(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 970,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LnudgeTCLM',            &
     &                      LnudgeTCLM(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 976,                           &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  Write out Fennel et al. (2006) ecosystem model parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'BioIter',               &
     &                      BioIter(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 21,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'AttSW',                 &
     &                      AttSW(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 27,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'AttChl',                &
     &                      AttChl(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 33,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'PARfrac',               &
     &                      PARfrac(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 39,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Vp0',                   &
     &                      Vp0(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 45,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'I_thNH4',               &
     &                      I_thNH4(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 51,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'D_p5NH4',               &
     &                      D_p5NH4(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 57,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'NitriR',                &
     &                      NitriR(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 63,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'K_NO3',                 &
     &                      K_NO3(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 69,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'K_NH4',                 &
     &                      K_NH4(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 75,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'K_Phy',                 &
     &                      K_Phy(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 81,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Chl2C_m',               &
     &                      Chl2C_m(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 87,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ChlMin',                &
     &                      ChlMin(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 93,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'PhyCN',                 &
     &                      PhyCN(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 99,                            &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'PhyIP',                 &
     &                      PhyIP(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 105,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'PhyIS',                 &
     &                      PhyIS(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 111,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'PhyMin',                &
     &                      PhyMin(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 117,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'PhyMR',                 &
     &                      PhyMR(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 123,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooAE_N',               &
     &                      ZooAE_N(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 129,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooBM',                 &
     &                      ZooBM(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 135,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooCN',                 &
     &                      ZooCN(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 141,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooER',                 &
     &                      ZooER(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 147,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooGR',                 &
     &                      ZooGR(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 153,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooMin',                &
     &                      ZooMin(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 159,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooMR',                 &
     &                      ZooMR(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 165,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'LDeRRN',                &
     &                      LDeRRN(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 171,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'LDeRRC',                &
     &                      LDeRRC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 177,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'CoagR',                 &
     &                      CoagR(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 183,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'SDeRRN',                &
     &                      SDeRRN(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 189,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'SDeRRC',                &
     &                      SDeRRC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 195,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'wPhy',                  &
     &                      wPhy(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 201,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'wLDet',                 &
     &                      wLDet(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 207,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'wSDet',                 &
     &                      wSDet(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 213,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'pCO2air',               &
     &                      pCO2air(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 219,                           &
     &               "ROMS/Nonlinear/Biology/fennel_wrt.h")) RETURN
!
!-----------------------------------------------------------------------
!  Write out grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch. Writing characters in parallel I/O is extremely
!  inefficient.  It is better to write this as an integer switch:
!  0=Cartesian, 1=spherical.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'spherical',             &
     &                      spherical, (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1514,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  Domain Length.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'xl',                    &
     &                      xl(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1522,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'el',                    &
     &                      el(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1528,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  S-coordinate parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'Vtransform',            &
     &                      Vtransform(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1538,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'Vstretching',           &
     &                      Vstretching(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1544,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'theta_s',               &
     &                      theta_s(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1550,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'theta_b',               &
     &                      theta_b(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1556,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Tcline',                &
     &                      Tcline(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1562,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'hc',                    &
     &                      hc(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1568,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  SGRID conventions for staggered data on structured grids. The value
!  is arbitrary but is set to unity so it can be used as logical during
!  post-processing.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'grid',                  &
     &                      (/1/), (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1578,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  S-coordinate non-dimensional independent variables.
!
      CALL netcdf_put_fvar (ng, model, ncname, 's_rho',                 &
     &                      SCALARS(ng)%sc_r(:), (/1/), (/N(ng)/),      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1586,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 's_w',                   &
     &                      SCALARS(ng)%sc_w(0:), (/1/), (/N(ng)+1/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1592,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  S-coordinate non-dimensional stretching curves.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Cs_r',                  &
     &                      SCALARS(ng)%Cs_r(:), (/1/), (/N(ng)/),      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1600,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Cs_w',                  &
     &                      SCALARS(ng)%Cs_w(0:), (/1/), (/N(ng)+1/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1606,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'user',                &
     &                        user, (/1/), (/Nuser/),                   &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 1616,                        &
     &                 "ROMS/Utility/wrt_info.F")) RETURN
      END IF
      GRID_VARS : IF (ncid.ne.FLT(ng)%ncid) THEN
!
!  Bathymetry.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_r8
          IF (ncid.ne.STA(ng)%ncid) THEN
            IF (find_string(var_name, n_var, 'h', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % h,                          &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 1659,                  &
     &                       "ROMS/Utility/wrt_info.F")) THEN
                IF (Master) WRITE (stdout,10) 'h', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'h', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Coriolis parameter.
!
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'f', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % f,                          &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 1729,                  &
     &                       "ROMS/Utility/wrt_info.F")) THEN
                IF (Master) WRITE (stdout,10) 'f', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'f', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Curvilinear transformation metrics.
!
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'pm', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % pm,                         &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 1756,                  &
     &                       "ROMS/Utility/wrt_info.F")) THEN
                IF (Master) WRITE (stdout,10) 'pm', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'pm', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'pn', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % pn,                         &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 1781,                  &
     &                       "ROMS/Utility/wrt_info.F")) THEN
                IF (Master) WRITE (stdout,10) 'pn', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'pn', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_r8
            IF (ncid.ne.STA(ng)%ncid) THEN
              IF (find_string(var_name, n_var, 'lon_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % lonr,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 1809,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'lon_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lon_rho', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_r8
            IF (ncid.ne.STA(ng)%ncid) THEN
              IF (find_string(var_name, n_var, 'lat_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % latr,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 1847,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'lat_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_rho', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_r8
            IF (ncid.ne.STA(ng)%ncid) THEN
              IF (find_string(var_name, n_var, 'x_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % xr,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 1887,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'x_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_rho', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_r8
            IF (ncid.ne.STA(ng)%ncid) THEN
              IF (find_string(var_name, n_var, 'y_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % yr,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 1925,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'y_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_rho', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lon_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % lonu,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 1967,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'lon_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lon_u', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lat_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % latu,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 1992,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'lat_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_u', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'x_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % xu,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 2019,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'x_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_u', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'y_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % yu,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 2044,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'y_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_u', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lon_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % lonv,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 2073,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'lon_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) 'lon_v', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lat_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % latv,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 2098,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'lat_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_v', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'x_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % xv,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 2125,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'x_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) 'x_v', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'y_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % yv,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 2150,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'y_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_v', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lon_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % lonp,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 2179,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'lon_p', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lon_p', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lat_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % latp,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 2204,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'lat_p', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_p', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'x_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % xp,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 2231,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'x_p', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_p', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'y_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % yp,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr, 2256,                &
     &                         "ROMS/Utility/wrt_info.F")) THEN
                  IF (Master) WRITE (stdout,10) 'y_p', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_p', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
      END IF GRID_VARS
!
!-----------------------------------------------------------------------
!  Synchronize NetCDF file to disk to allow other processes to access
!  data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, ncname, ncid)
      IF (FoundError(exit_flag, NoError, 2582,                          &
     &               "ROMS/Utility/wrt_info.F")) RETURN
!
  10  FORMAT (/,' WRT_INFO - error while writing variable: ',a,/,       &
     &        12x,'into file: ',a)
  20  FORMAT (/,' WRT_INFO - error while inquiring ID for variable: ',  &
     &        a,/,12x,'in file: ',a)
  30  FORMAT (/,' WRT_INFO - unable to synchronize to disk file: ',     &
     &        /,12x,a)
      RETURN
      END SUBROUTINE wrt_info
