      MODULE biology_mod
!
!svn $Id: biology.F 851 2017-06-23 23:22:54Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the source and sink terms for selected        !
!  biology model.                                                      !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: biology
      CONTAINS
      SUBROUTINE biology (ng,tile)
!
!svn $Id: fennel.h 864 2017-08-10 04:11:10Z arango $
!***********************************************************************
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                               Katja Fennel   !
!****************************************** Alexander F. Shchepetkin ***
!                                                                      !
!  This routine computes the  biological sources and sinks for the     !
!  Fennel et at. (2006) ecosystem model. Then, it adds those terms     !
!  to the global biological fields.                                    !
!                                                                      !
!  This model is loosely based on the model by Fasham et al. (1990)    !
!  but it differs in many respects.  The detailed equations of the     !
!  nitrogen cycling component  are given in  Fennel et al. (2006).     !
!  Nitrogen is the  fundamental elemental  currency in this model.     !
!  This model was adapted from a code written originally  by  John     !
!  Moisan and Emanule DiLorenzo.                                       !
!                                                                      !
!  It is recommended to activate always the  "BIO_SEDIMENT" option     !
!  to ensure conservation of mass by converting the organic matter     !
!  that is sinking out of the bottom most grid cell into inorganic     !
!  nutrients (i.e.,  instantanaous remineralization  at the water-     !
!  sediment interface). Additionally, the "DENITRIFICATION" option     !
!  can be activated.  Hence, a fraction of the instantenous bottom     !
!  remineralization is  assumed to  occur  through  the  anearobic     !
!  (denitrification)  pathway  and  thus  lost  from the  pool  of     !
!  biologically availalbe fixed nitrogen. See Fennel et al. (2006)     !
!  for details.                                                        !
!                                                                      !
!  Additional  options can be  activated to  enable  simulation of     !
!  inorganic carbon and dissolved oxygen.  Accounting of inorganic     !
!  carbon is activated by the "CARBON" option,  and results in two     !
!  additional  biological  tracer  variables:  DIC and alkalinity.     !
!  See Fennel et al. (2008) for details.                               !
!                                                                      !
!  If the "pCO2_RZ" options is activated, in addition to "CARBON",     !
!  the carbonate system  routines by Zeebe and Wolf-Gladrow (2001)     !
!  are used,  while the  OCMIP  standard routines are the default.     !
!  There are two different ways of treating alkalinity.  It can be     !
!  treated diagnostically (default),  in this case alkalinity acts     !
!  like a passive tracer  that is  not affected  by changes in the     !
!  concentration of  nitrate or ammonium.  However,  if the option     !
!  "TALK_NONCONSERV" is used,  the alkalinity  will be affected by     !
!  sources and sinks in nitrate. See Fennel et al. (2008) for more     !
!  details.                                                            !
!                                                                      !
!  If the "OXYGEN" option is activated,  one additional biological     !
!  tracer variable for dissolved oxygen. "OXYGEN" can be activated     !
!  independently of the  "CARBON"  option. If "OCMIP_OXYGEN_SC" is     !
!  used, in addition to "OXYGEN",  the Schmidt number of oxygen in     !
!  seawater will be  computed  using the  formulation  proposed by     !
!  Keeling et al. (1998, Global Biogeochem. Cycles,  12, 141-163).     !
!  Otherwise, the Wanninkhof^s (1992) formula will be used.            !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Fennel, K., Wilkin, J., Levin, J., Moisan, J., O^Reilly, J.,      !
!      Haidvogel, D., 2006: Nitrogen cycling in the Mid Atlantic       !
!      Bight and implications for the North Atlantic nitrogen          !
!      budget: Results from a three-dimensional model.  Global         !
!      Biogeochemical Cycles 20, GB3007, doi:10.1029/2005GB002456.     !
!                                                                      !
!    Fennel, K., Wilkin, J., Previdi, M., Najjar, R. 2008:             !
!      Denitrification effects on air-sea CO2 flux in the coastal      !
!      ocean: Simulations for the Northwest North Atlantic.            !
!      Geophys. Res. Letters 35, L24608, doi:10.1029/2008GL036147.     !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_diags
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
      implicit none
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
!  Set header file name.
!
      IF (Lbiofile(iNLM).and.(tile.eq.0)) THEN
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)="ROMS/Nonlinear/Biology/fennel.h"
      END IF
!
      CALL wclock_on (ng, iNLM, 15, 104, "ROMS/Nonlinear/Biology/fennel.h")
      CALL biology_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   FORCES(ng) % srflx,                            &
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
     &                   OCEAN(ng) % pH,                                &
     &                   DIAGS(ng) % DiaBio2d,                          &
     &                   DIAGS(ng) % DiaBio3d,                          &
     &                   OCEAN(ng) % t)
      CALL wclock_off (ng, iNLM, 15, 139, "ROMS/Nonlinear/Biology/fennel.h")
      RETURN
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
     &                         Hz, z_r, z_w, srflx,                     &
     &                         sustr, svstr,                            &
     &                         pH,                                      &
     &                         DiaBio2d, DiaBio3d,                      &
     &                         t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
      USE dateclock_mod, ONLY : caldate
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(inout) :: pH(LBi:,LBj:)
      real(r8), intent(inout) :: DiaBio2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaBio3d(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
!
!  Local variable declarations.
!
      integer, parameter :: Nsink = 6
      integer :: Iter, i, ibio, isink, itrc, ivar, j, k, ks
      integer, dimension(Nsink) :: idsink
      real(r8), parameter :: eps = 1.0e-20_r8
      real(r8) :: u10squ
      integer, parameter :: DoNewton = 0            ! pCO2 solver
      real(r8), parameter :: Acoef = 2073.1_r8      ! Schmidt
      real(r8), parameter :: Bcoef = 125.62_r8      ! number
      real(r8), parameter :: Ccoef = 3.6276_r8      ! transfer
      real(r8), parameter :: Dcoef = 0.043219_r8    ! coefficients
      real(r8), parameter :: A1 = -60.2409_r8       ! surface
      real(r8), parameter :: A2 = 93.4517_r8        ! CO2
      real(r8), parameter :: A3 = 23.3585_r8        ! solubility
      real(r8), parameter :: B1 = 0.023517_r8       ! coefficients
      real(r8), parameter :: B2 = -0.023656_r8
      real(r8), parameter :: B3 = 0.0047036_r8
      real(r8) :: pmonth                         ! months since Jan 1951
      real(r8) :: pCO2air_secular
      real(r8) :: yday
      real(r8), parameter :: pi2 = 6.2831853071796_r8
      real(r8), parameter :: D0 = 282.6_r8          ! coefficients
      real(r8), parameter :: D1 = 0.125_r8          ! to calculate
      real(r8), parameter :: D2 =-7.18_r8           ! secular trend in
      real(r8), parameter :: D3 = 0.86_r8           ! atmospheric pCO2
      real(r8), parameter :: D4 =-0.99_r8
      real(r8), parameter :: D5 = 0.28_r8
      real(r8), parameter :: D6 =-0.80_r8
      real(r8), parameter :: D7 = 0.06_r8
      real(r8) :: Att, AttFac, ExpAtt, Itop, PAR
      real(r8) :: Epp, L_NH4, L_NO3, LTOT, Vp
      real(r8) :: Chl2C, dtdays, t_PPmax, inhNH4
      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5
      real(r8) :: fac1, fac2, fac3
      real(r8) :: cffL, cffR, cu, dltL, dltR
      real(r8) :: total_N
      real(r8) :: fiter
      real(r8) :: C_Flux_RemineL, C_Flux_RemineS
      real(r8) :: CO2_Flux, CO2_sol, SchmidtN, TempK
      real(r8) :: N_Flux_Assim
      real(r8) :: N_Flux_CoagD, N_Flux_CoagP
      real(r8) :: N_Flux_Egest
      real(r8) :: N_Flux_NewProd, N_Flux_RegProd
      real(r8) :: N_Flux_Nitrifi
      real(r8) :: N_Flux_Pmortal, N_Flux_Zmortal
      real(r8) :: N_Flux_RemineL, N_Flux_RemineS
      real(r8) :: N_Flux_Zexcret, N_Flux_Zmetabo
      real(r8), dimension(Nsink) :: Wbio
      integer, dimension(IminS:ImaxS,N(ng)) :: ksource
      real(r8), dimension(IminS:ImaxS) :: PARsur
      real(r8), dimension(IminS:ImaxS) :: pCO2
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc
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
!-----------------------------------------------------------------------
! If appropriate, initialize time-averaged diagnostic arrays.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsDIA(ng)).and.                                 &
     &     (MOD(iic(ng),nDIA(ng)).eq.1)).or.                            &
     &    ((iic(ng).ge.ntsDIA(ng)).and.(nDIA(ng).eq.1)).or.             &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
        DO ivar=1,NDbio2d
          DO j=Jstr,Jend
            DO i=Istr,Iend
              DiaBio2d(i,j,ivar)=0.0_r8
            END DO
          END DO
        END DO
        DO ivar=1,NDbio3d
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                DiaBio3d(i,j,k,ivar)=0.0_r8
              END DO
            END DO
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Avoid computing source/sink terms if no biological iterations.
!
      IF (BioIter(ng).le.0) RETURN
!
!  Set time-stepping according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)
!
!  A factor to account for the number of iterations in accumulating
!  diagnostic rate variables.
!
      fiter=1.0_r8/REAL(BioIter(ng),r8)
!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iPhyt
      idsink(2)=iChlo
      idsink(3)=iSDeN
      idsink(4)=iLDeN
      idsink(5)=iSDeC
      idsink(6)=iLDeC
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wPhy(ng)                ! phytoplankton
      Wbio(2)=wPhy(ng)                ! chlorophyll
      Wbio(3)=wSDet(ng)               ! small Nitrogen-detritus
      Wbio(4)=wLDet(ng)               ! large Nitrogen-detritus
      Wbio(5)=wSDet(ng)               ! small Carbon-detritus
      Wbio(6)=wLDet(ng)               ! large Carbon-detritus
!
!  Compute inverse thickness to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_old(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
              Bio(i,k,ibio)=Bio_old(i,k,ibio)
            END DO
          END DO
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio_old(i,k,iTIC_)=MIN(Bio_old(i,k,iTIC_),3000.0_r8)
            Bio_old(i,k,iTIC_)=MAX(Bio_old(i,k,iTIC_),400.0_r8)
            Bio(i,k,iTIC_)=Bio_old(i,k,iTIC_)
          END DO
        END DO
!
!  Extract potential temperature and salinity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),35.0_r8)
            Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
!  During the iterative procedure a series of fractional time steps are
!  performed in a chained mode (splitting by different biological
!  conversion processes) in sequence of the main food chain.  In all
!  stages the concentration of the component being consumed is treated
!  in fully implicit manner, so the algorithm guarantees non-negative
!  values, no matter how strong s the concentration of active consuming
!  component (Phytoplankton or Zooplankton).  The overall algorithm,
!  as well as any stage of it, is formulated in conservative form
!  (except explicit sinking) in sense that the sum of concentration of
!  all components is conserved.
!
!
!  In the implicit algorithm, we have for example (N: nitrate,
!                                                  P: phytoplankton),
!
!     N(new) = N(old) - uptake * P(old)     uptake = mu * N / (Kn + N)
!                                                    {Michaelis-Menten}
!  below, we set
!                                           The N in the numerator of
!     cff = mu * P(old) / (Kn + N(old))     uptake is treated implicitly
!                                           as N(new)
!
!  so the time-stepping of the equations becomes:
!
!     N(new) = N(old) / (1 + cff)     (1) when substracting a sink term,
!                                         consuming, divide by (1 + cff)
!  and
!
!     P(new) = P(old) + cff * N(new)  (2) when adding a source term,
!                                         growing, add (cff * source)
!
!  Notice that if you substitute (1) in (2), you will get:
!
!     P(new) = P(old) + cff * N(old) / (1 + cff)    (3)
!
!  If you add (1) and (3), you get
!
!     N(new) + P(new) = N(old) + P(old)
!
!  implying conservation regardless how "cff" is computed. Therefore,
!  this scheme is unconditionally stable regardless of the conversion
!  rate. It does not generate negative values since the constituent
!  to be consumed is always treated implicitly. It is also biased
!  toward damping oscillations.
!
!  The iterative loop below is to iterate toward an universal Backward-
!  Euler treatment of all terms. So if there are oscillations in the
!  system, they are only physical oscillations. These iterations,
!  however, do not improve the accuaracy of the solution.
!
        ITER_LOOP: DO Iter=1,BioIter(ng)
!
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
!
!  Compute attenuation coefficient based on the concentration of
!  chlorophyll-a within each grid box.  Then, attenuate surface
!  photosynthetically available radiation (PARsur) down inot the
!  water column.  Thus, PAR at certain depth depends on the whole
!  distribution of chlorophyll-a above.
!  To compute rate of maximum primary productivity (t_PPmax), one needs
!  PAR somewhat in the middle of the gridbox, so that attenuation "Att"
!  corresponds to half of the grid box height, while PAR is multiplied
!  by it twice: once to get it in the middle of grid-box and once the
!  compute on the lower grid-box interface.
!
          DO i=Istr,Iend
            PAR=PARsur(i)
            AttFac=0.0_r8
            IF (PARsur(i).gt.0.0_r8) THEN
              DO k=N(ng),1,-1
!
!  Compute average light attenuation for each grid cell. To include
!  other attenuation contributions like suspended sediment or CDOM
!  modify AttFac.
!
                Att=(AttSW(ng)+                                         &
     &               AttChl(ng)*Bio(i,k,iChlo)+                         &
     &               AttFac)*                                           &
     &               (z_w(i,j,k)-z_w(i,j,k-1))
                ExpAtt=EXP(-Att)
                Itop=PAR
                PAR=Itop*(1.0_r8-ExpAtt)/Att    ! average at cell center
!
!  Compute Chlorophyll-a phytoplankton ratio, [mg Chla / (mg C)].
!
                cff=PhyCN(ng)*12.0_r8
                Chl2C=MIN(Bio(i,k,iChlo)/(Bio(i,k,iPhyt)*cff+eps),      &
     &                    Chl2C_m(ng))
!
!  Temperature-limited and light-limited growth rate (Eppley, R.W.,
!  1972, Fishery Bulletin, 70: 1063-1085; here 0.59=ln(2)*0.851).
!  Check value for Vp is 2.9124317 at 19.25 degC.
!
                Vp=Vp0(ng)*0.59_r8*(1.066_r8**Bio(i,k,itemp))
                fac1=PAR*PhyIS(ng)
                Epp=Vp/SQRT(Vp*Vp+fac1*fac1)
                t_PPmax=Epp*fac1
!
!  Nutrient-limitation terms (Parker 1993 Ecol Mod., 66, 113-120).
!
                cff1=Bio(i,k,iNH4_)*K_NH4(ng)
                cff2=Bio(i,k,iNO3_)*K_NO3(ng)
                inhNH4=1.0_r8/(1.0_r8+cff1)
                L_NH4=cff1/(1.0_r8+cff1)
                L_NO3=cff2*inhNH4/(1.0_r8+cff2)
                LTOT=L_NO3+L_NH4
!
!  Nitrate and ammonium uptake by Phytoplankton.
!
                fac1=dtdays*t_PPmax
                cff4=fac1*K_NO3(ng)*inhNH4/(1.0_r8+cff2)*Bio(i,k,iPhyt)
                cff5=fac1*K_NH4(ng)/(1.0_r8+cff1)*Bio(i,k,iPhyt)
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                N_Flux_NewProd=Bio(i,k,iNO3_)*cff4
                N_Flux_RegProd=Bio(i,k,iNH4_)*cff5
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                          &
     &                         N_Flux_NewProd+N_Flux_RegProd
!
                Bio(i,k,iChlo)=Bio(i,k,iChlo)+                          &
     &                         (dtdays*t_PPmax*t_PPmax*LTOT*LTOT*       &
     &                          Chl2C_m(ng)*Bio(i,k,iChlo))/            &
     &                         (PhyIS(ng)*MAX(Chl2C,eps)*PAR+eps)
                DiaBio3d(i,j,k,iPPro)=DiaBio3d(i,j,k,iPPro)+            &
     &                                (N_Flux_NewProd+N_Flux_RegProd)*  &
     &                                fiter
                DiaBio3d(i,j,k,iNO3u)=DiaBio3d(i,j,k,iNO3u)+            &
     &                                N_Flux_NewProd*fiter
!
!  Total inorganic carbon (CO2) uptake during phytoplankton growth.
!
                cff1=PhyCN(ng)*(N_Flux_NewProd+N_Flux_RegProd)
                Bio(i,k,iTIC_)=Bio(i,k,iTIC_)-cff1
!
! The Nitrification of NH4 ==> NO3 is thought to occur only in dark and
! only in aerobic water (see Olson, R. J., 1981, JMR: (39), 227-238.).
!
!         NH4+ + 3/2 O2  ==> NO2- + H2O;  via Nitrosomonas bacteria
!         NO2-  + 1/2 O2 ==> NO3-      ;  via Nitrobacter  bacteria
!
! Note that the entire process has a total loss of two moles of O2 per
! mole of NH4. If we were to resolve NO2 profiles, this is where we
! would change the code to split out the differential effects of the
! two different bacteria types. If OXYGEN is defined, nitrification is
! inhibited at low oxygen concentrations using a Michaelis-Menten term.
!
                fac1=dtdays*NitriR(ng)
                cff1=(PAR-I_thNH4(ng))/                                 &
     &               (D_p5NH4(ng)+PAR-2.0_r8*I_thNH4(ng))
                cff2=1.0_r8-MAX(0.0_r8,cff1)
                cff3=fac1*cff2
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                PAR=Itop*ExpAtt
              END DO
!
!  If PARsur=0, nitrification occurs at the maximum rate (NitriR).
!
            ELSE
              cff3=dtdays*NitriR(ng)
              DO k=N(ng),1,-1
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
              END DO
            END IF
          END DO
!
!-----------------------------------------------------------------------
!  Phytoplankton grazing by zooplankton (rate: ZooGR), phytoplankton
!  assimilated to zooplankton (fraction: ZooAE_N) and egested to small
!  detritus, and phytoplankton mortality (rate: PhyMR) to small
!  detritus. [Landry 1993 L&O 38:468-472]
!-----------------------------------------------------------------------
!
          fac1=dtdays*ZooGR(ng)
          cff2=dtdays*PhyMR(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
!
! Phytoplankton grazing by zooplankton.
!
              cff1=fac1*Bio(i,k,iZoop)*Bio(i,k,iPhyt)/                  &
     &             (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
              cff3=1.0_r8/(1.0_r8+cff1)
              Bio(i,k,iPhyt)=cff3*Bio(i,k,iPhyt)
              Bio(i,k,iChlo)=cff3*Bio(i,k,iChlo)
!
! Phytoplankton assimilated to zooplankton and egested to small
! detritus.
!
              N_Flux_Assim=cff1*Bio(i,k,iPhyt)*ZooAE_N(ng)
              N_Flux_Egest=Bio(i,k,iPhyt)*cff1*(1.0_r8-ZooAE_N(ng))
              Bio(i,k,iZoop)=Bio(i,k,iZoop)+                            &
     &                       N_Flux_Assim
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+                            &
     &                       N_Flux_Egest
!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
              N_Flux_Pmortal=cff2*MAX(Bio(i,k,iPhyt)-PhyMin(ng),0.0_r8)
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)-N_Flux_Pmortal
              Bio(i,k,iChlo)=Bio(i,k,iChlo)-                            &
     &                       cff2*MAX(Bio(i,k,iChlo)-ChlMin(ng),0.0_r8)
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+                            &
     &                       N_Flux_Pmortal
              Bio(i,k,iSDeC)=Bio(i,k,iSDeC)+                            &
     &                       PhyCN(ng)*(N_Flux_Egest+N_Flux_Pmortal)+   &
     &                       (PhyCN(ng)-ZooCN(ng))*N_Flux_Assim
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Zooplankton basal metabolism to NH4  (rate: ZooBM), zooplankton
!  mortality to small detritus (rate: ZooMR), zooplankton ingestion
!  related excretion (rate: ZooER).
!-----------------------------------------------------------------------
!
          cff1=dtdays*ZooBM(ng)
          fac2=dtdays*ZooMR(ng)
          fac3=dtdays*ZooER(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              fac1=fac3*Bio(i,k,iPhyt)*Bio(i,k,iPhyt)/                  &
     &             (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
              cff2=fac2*Bio(i,k,iZoop)
              cff3=fac1*ZooAE_N(ng)
              Bio(i,k,iZoop)=Bio(i,k,iZoop)/                            &
     &                       (1.0_r8+cff2+cff3)
!
!  Zooplankton mortality and excretion.
!
              N_Flux_Zmortal=cff2*Bio(i,k,iZoop)
              N_Flux_Zexcret=cff3*Bio(i,k,iZoop)
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zexcret
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Zmortal
!
!  Zooplankton basal metabolism (limited by a zooplankton minimum).
!
              N_Flux_Zmetabo=cff1*MAX(Bio(i,k,iZoop)-ZooMin(ng),0.0_r8)
              Bio(i,k,iZoop)=Bio(i,k,iZoop)-N_Flux_Zmetabo
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zmetabo
              Bio(i,k,iSDeC)=Bio(i,k,iSDeC)+                            &
     &                       ZooCN(ng)*N_Flux_Zmortal
              Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+                            &
     &                       ZooCN(ng)*(N_Flux_Zmetabo+N_Flux_Zexcret)
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Coagulation of phytoplankton and small detritus to large detritus.
!-----------------------------------------------------------------------
!
          fac1=dtdays*CoagR(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=fac1*(Bio(i,k,iSDeN)+Bio(i,k,iPhyt))
              cff2=1.0_r8/(1.0_r8+cff1)
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff2
              Bio(i,k,iChlo)=Bio(i,k,iChlo)*cff2
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
              N_Flux_CoagP=Bio(i,k,iPhyt)*cff1
              N_Flux_CoagD=Bio(i,k,iSDeN)*cff1
              Bio(i,k,iLDeN)=Bio(i,k,iLDeN)+                            &
     &                       N_Flux_CoagP+N_Flux_CoagD
              Bio(i,k,iSDeC)=Bio(i,k,iSDeC)-PhyCN(ng)*N_Flux_CoagD
              Bio(i,k,iLDeC)=Bio(i,k,iLDeC)+                            &
     &                       PhyCN(ng)*(N_Flux_CoagP+N_Flux_CoagD)
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Detritus recycling to NH4, remineralization.
!-----------------------------------------------------------------------
!
          cff1=dtdays*SDeRRN(ng)
          cff2=1.0_r8/(1.0_r8+cff1)
          cff3=dtdays*LDeRRN(ng)
          cff4=1.0_r8/(1.0_r8+cff3)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
              Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
              N_Flux_RemineS=Bio(i,k,iSDeN)*cff1
              N_Flux_RemineL=Bio(i,k,iLDeN)*cff3
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+                            &
     &                       N_Flux_RemineS+N_Flux_RemineL
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Allow different remineralization rates for detrital C and detrital N.
!-----------------------------------------------------------------------
!
          cff1=dtdays*SDeRRC(ng)
          cff2=1.0_r8/(1.0_r8+cff1)
          cff3=dtdays*LDeRRC(ng)
          cff4=1.0_r8/(1.0_r8+cff3)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,iSDeC)=Bio(i,k,iSDeC)*cff2
              Bio(i,k,iLDeC)=Bio(i,k,iLDeC)*cff4
              C_Flux_RemineS=Bio(i,k,iSDeC)*cff1
              C_Flux_RemineL=Bio(i,k,iLDeC)*cff3
              Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+                            &
     &                       C_Flux_RemineS+C_Flux_RemineL
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Surface CO2 gas exchange.
!-----------------------------------------------------------------------
!
!  Compute equilibrium partial pressure inorganic carbon (ppmv) at the
!  surface.
!
          k=N(ng)
          CALL pCO2_water (Istr, Iend, LBi, UBi, LBj, UBj,              &
     &                     IminS, ImaxS, j, DoNewton,                   &
     &                     Bio(IminS:,k,itemp), Bio(IminS:,k,isalt),    &
     &                     Bio(IminS:,k,iTIC_), Bio(IminS:,k,iTAlk),    &
     &                     0.0_r8, 0.0_r8, pH, pCO2)
!
!  Compute surface CO2 gas exchange.
!
          cff1=rho0*550.0_r8
          cff2=dtdays*0.31_r8*24.0_r8/100.0_r8
          DO i=Istr,Iend
!
!  Compute CO2 transfer velocity : u10squared (u10 in m/s)
!
            u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+     &
     &                       (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
            SchmidtN=Acoef-                                             &
     &               Bio(i,k,itemp)*(Bcoef-                             &
     &                               Bio(i,k,itemp)*(Ccoef-             &
     &                               Bio(i,k,itemp)*Dcoef))
            cff3=cff2*u10squ*SQRT(660.0_r8/SchmidtN)
!
!  Calculate CO2 solubility [mol/(kg.atm)] using Weiss (1974) formula.
!
            TempK=0.01_r8*(Bio(i,k,itemp)+273.15_r8)
            CO2_sol=EXP(A1+                                             &
     &                  A2/TempK+                                       &
     &                  A3*LOG(TempK)+                                  &
     &                  Bio(i,k,isalt)*(B1+TempK*(B2+B3*TempK)))
!
!  Add in CO2 gas exchange.
!
            CALL caldate (tdays(ng), yd_r8=yday)
            pmonth=2003.0_r8-1951.0_r8+yday/365.0_r8
            CO2_Flux=cff3*CO2_sol*(pCO2air(ng)-pCO2(i))
            Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+                              &
     &                     CO2_Flux*Hz_inv(i,k)
            DiaBio2d(i,j,iCOfx)=DiaBio2d(i,j,iCOfx)+                    &
     &                          CO2_Flux*fiter
            DiaBio2d(i,j,ipCO2)=pCO2(i)
          END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
          SINK_LOOP: DO isink=1,Nsink
            ibio=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=Bio(i,k,ibio)
              END DO
            END DO
!
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                cff=(dltR-dltL)*Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                WL(i,k)=(dltR-2.0_r8*dltL)**2
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                dltR=MAX(cff,WR(i,k+1))
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
            END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                dltR=bR(i,k)-qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                cffR=2.0_r8*dltR
                cffL=2.0_r8*dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
              END DO
            END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
            cff=dtdays*ABS(Wbio(isink))
            DO k=1,N(ng)
              DO i=Istr,Iend
                FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                  END IF
                END DO
              END DO
            END DO
!
!  Finalize computation of flux: add fractional part.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,ibio)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO
!
!  Particulate flux reaching the seafloor is remineralized and returned
!  to the dissolved nitrate pool. Without this conversion, particulate
!  material falls out of the system. This is a temporary fix to restore
!  total nitrogen conservation. It will be replaced later by a
!  parameterization that includes the time delay of remineralization
!  and dissolved oxygen.
!
            cff2=4.0_r8/16.0_r8
            IF ((ibio.eq.iPhyt).or.                                     &
     &          (ibio.eq.iSDeN).or.                                     &
     &          (ibio.eq.iLDeN)) THEN
              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1*cff2
                DiaBio2d(i,j,iDNIT)=DiaBio2d(i,j,iDNIT)+                &
     &                              (1.0_r8-cff2)*cff1*Hz(i,j,1)*fiter
              END DO
            END IF
            cff3=12.0_r8
            cff4=0.74_r8
            IF ((ibio.eq.iSDeC).or.                                     &
     &          (ibio.eq.iLDeC))THEN
              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iTIC_)=Bio(i,1,iTIC_)+cff1
              END DO
            END IF
            IF (ibio.eq.iPhyt)THEN
              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iTIC_)=Bio(i,1,iTIC_)+cff1*PhyCN(ng)
              END DO
            END IF
          END DO SINK_LOOP
        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables: Add increment due to BGC processes
!  to tracer array in time index "nnew". Index "nnew" is solution after
!  advection and mixing and has transport units (m Tunits) hence the
!  increment is multiplied by Hz.  Notice that we need to subtract
!  original values "Bio_old" at the top of the routine to just account
!  for the concentractions affected by BGC processes. This also takes
!  into account any constraints (non-negative concentrations, carbon
!  concentration range) specified before entering BGC kernel. If "Bio"
!  were unchanged by BGC processes, the increment would be exactly
!  zero. Notice that final tracer values, t(:,:,:,nnew,:) are not
!  bounded >=0 so that we can preserve total inventory of N and
!  C even when advection causes tracer concentration to go negative.
!  (J. Wilkin and H. Arango, Apr 27, 2012)
!-----------------------------------------------------------------------
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio(i,k,ibio)-Bio_old(i,k,ibio)
              t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
            END DO
          END DO
        END DO
      END DO J_LOOP
      RETURN
      END SUBROUTINE biology_tile
      SUBROUTINE pCO2_water (Istr, Iend,                                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, j, DoNewton,                 &
     &                       T, S, TIC, TAlk, PO4, SiO3, pH, pCO2)
!
!***********************************************************************
!                                                                      !
!  This routine computes equilibrium partial pressure of CO2 (pCO2)    !
!  in the surface seawater.                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     LBj        J-dimension lower bound.                              !
!     UBj        J-dimension upper bound.                              !
!     IminS      I-dimension lower bound for private arrays.           !
!     ImaxS      I-dimension upper bound for private arrays.           !
!     j          j-pipelined index.                                    !
!     DoNewton   Iteration solver:                                     !
!                  [0] Bracket and bisection.                          !
!                  [1] Newton-Raphson method.                          !
!     rmask      Land/Sea masking.                                     !
!     T          Surface temperature (Celsius).                        !
!     S          Surface salinity (PSS).                               !
!     TIC        Total inorganic carbon (millimol/m3).                 !
!     TAlk       Total alkalinity (milli-equivalents/m3).              !
!     PO4        Inorganic phosphate (millimol/m3).                    !
!     SiO3       Inorganic silicate (millimol/m3).                     !
!     pH         Best pH guess.                                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     pCO2       partial pressure of CO2 (ppmv).                       !
!                                                                      !
!  Check Value:  (T=24, S=36.6, TIC=2040, TAlk=2390, PO4=0,            !
!                 SiO3=0, pH=8)                                        !
!                                                                      !
!                pcO2=0.35074945E+03 ppmv  (DoNewton=0)                !
!                pCO2=0.35073560E+03 ppmv  (DoNewton=1)                !
!                                                                      !
!  This subroutine was adapted by Mick Follows (Oct 1999) from OCMIP2  !
!  code CO2CALC. Modified for ROMS by Hernan Arango (Nov 2003).        !
!                                                                      !
!***********************************************************************
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j, DoNewton
!
      real(r8), intent(in) :: T(IminS:)
      real(r8), intent(in) :: S(IminS:)
      real(r8), intent(in) :: TIC(IminS:)
      real(r8), intent(in) :: TAlk(IminS:)
      real(r8), intent(inout) :: pH(LBi:,LBj:)
      real(r8), intent(in) :: PO4
      real(r8), intent(in) :: SiO3
      real(r8), intent(out) :: pCO2(IminS:ImaxS)
!
!  Local variable declarations.
!
      integer, parameter :: InewtonMax = 10
      integer, parameter :: IbrackMax = 30
      integer :: Hstep, Ibrack, Inewton, i
      real(r8) :: Tk, centiTk, invTk, logTk
      real(r8) :: SO4, scl, sqrtS, sqrtSO4
      real(r8) :: alk, dic, phos, sili
      real(r8) :: borate, sulfate, fluoride
      real(r8) :: ff, K1, K2, K1p, K2p, K3p, Kb, Kf, Ks, Ksi, Kw
      real(r8) :: K12, K12p, K123p, invKb, invKs, invKsi
      real(r8) :: A, A2, B, B2, C, dA, dB
      real(r8) :: df, fn, fni(3), ftest
      real(r8) :: deltaX, invX, invX2, X, X2, X3
      real(r8) :: pH_guess, pH_hi, pH_lo
      real(r8) :: X_guess, X_hi, X_lo, X_mid
      real(r8) :: CO2star, Htotal, Htotal2
!
!=======================================================================
!  Determine coefficients for surface carbon chemisty.  If land/sea
!  masking, compute only on water points.
!=======================================================================
!
      I_LOOP: DO i=Istr,Iend
        Tk=T(i)+273.15_r8
        centiTk=0.01_r8*Tk
        invTk=1.0_r8/Tk
        logTk=LOG(Tk)
        sqrtS=SQRT(S(i))
        SO4=19.924_r8*S(i)/(1000.0_r8-1.005_r8*S(i))
        sqrtSO4=SQRT(SO4)
        scl=S(i)/1.80655_r8
        alk=TAlk(i)*0.000001_r8
        dic=TIC(i)*0.000001_r8
        phos=PO4*0.000001_r8
        sili=SiO3*0.000001_r8
!
!-----------------------------------------------------------------------
!  Correction term for non-ideality, ff=k0*(1-pH2O). Equation 13 with
!  table 6 values from Weiss and Price (1980, Mar. Chem., 8, 347-359).
!-----------------------------------------------------------------------
!
        ff=EXP(-162.8301_r8+                                            &
     &         218.2968_r8/centiTk+                                     &
     &         LOG(centiTk)*90.9241_r8-                                 &
     &         centiTk*centiTk*1.47696_r8+                              &
     &         S(i)*(0.025695_r8-                                       &
     &               centiTk*(0.025225_r8-                              &
     &                        centiTk*0.0049867_r8)))
!
!-----------------------------------------------------------------------
!  Compute first (K1) and second (K2) dissociation constant of carboinic
!  acid:
!
!           K1 = [H][HCO3]/[H2CO3]
!           K2 = [H][CO3]/[HCO3]
!
!  From Millero (1995; page 664) using Mehrbach et al. (1973) data on
!  seawater scale.
!-----------------------------------------------------------------------
!
        K1=10.0_r8**(62.008_r8-                                         &
     &               invTk*3670.7_r8-                                   &
     &               logTk*9.7944_r8+                                   &
     &               S(i)*(0.0118_r8-                                   &
     &                     S(i)*0.000116_r8))
        K2=10.0_r8**(-4.777_r8-                                         &
     &               invTk*1394.7_r8+                                   &
     &               S(i)*(0.0184_r8-                                   &
     &                     S(i)*0.000118_r8))
!
!-----------------------------------------------------------------------
!  Compute dissociation constant of boric acid, Kb=[H][BO2]/[HBO2].
!  From Millero (1995; page 669) using data from Dickson (1990).
!-----------------------------------------------------------------------
!
        Kb=EXP(-invTk*(8966.90_r8+                                      &
     &                 sqrtS*(2890.53_r8+                               &
     &                        sqrtS*(77.942_r8-                         &
     &                               sqrtS*(1.728_r8-                   &
     &                                      sqrtS*0.0996_r8))))-        &
     &         logTk*(24.4344_r8+                                       &
     &                sqrtS*(25.085_r8+                                 &
     &                       sqrtS*0.2474_r8))+                         &
     &         Tk*(sqrtS*0.053105_r8)+                                  &
     &         148.0248_r8+                                             &
     &         sqrtS*(137.1942_r8+                                      &
     &                sqrtS*1.62142_r8))
!
!-----------------------------------------------------------------------
!  Compute first (K1p), second (K2p), and third (K3p) dissociation
!  constant of phosphoric acid:
!
!           K1p = [H][H2PO4]/[H3PO4]
!           K2p = [H][HPO4]/[H2PO4]
!           K3p = [H][PO4]/[HPO4]
!
!  From DOE (1994) equations 7.2.20, 7.2.23, and 7.2.26, respectively.
!  With footnote using data from Millero (1974).
!-----------------------------------------------------------------------
!
        K1p=EXP(115.525_r8-                                             &
     &          invTk*4576.752_r8-                                      &
     &          logTk*18.453_r8+                                        &
     &          sqrtS*(0.69171_r8-invTk*106.736_r8)-                    &
     &          S(i)*(0.01844_r8+invTk*0.65643_r8))
        K2p=EXP(172.0883_r8-                                            &
     &          invTk*8814.715_r8-                                      &
     &          logTk*27.927_r8+                                        &
     &          sqrtS*(1.3566_r8-invTk*160.340_r8)-                     &
     &          S(i)*(0.05778_r8-invTk*0.37335_r8))
        K3p=EXP(-18.141_r8-                                             &
     &          invTk*3070.75_r8+                                       &
     &          sqrtS*(2.81197_r8+invTk*17.27039_r8)-                   &
     &          S(i)*(0.09984_r8+invTk*44.99486_r8))
!
!-----------------------------------------------------------------------
!  Compute dissociation constant of silica, Ksi=[H][SiO(OH)3]/[Si(OH)4].
!  From Millero (1995; page 671) using data from Yao and Millero (1995).
!-----------------------------------------------------------------------
!
        Ksi=EXP(117.385_r8-                                             &
     &          invTk*8904.2_r8-                                        &
     &          logTk*19.334_r8+                                        &
     &          sqrtSO4*(3.5913_r8-invTk*458.79_r8)-                    &
     &          SO4*(1.5998_r8-invTk*188.74_r8-                         &
     &               SO4*(0.07871_r8-invTk*12.1652_r8))+                &
     &          LOG(1.0_r8-0.001005_r8*S(i)))
!
!-----------------------------------------------------------------------
!  Compute ion product of whater, Kw = [H][OH].
!  From Millero (1995; page 670) using composite data.
!-----------------------------------------------------------------------
!
        Kw=EXP(148.9652_r8-                                             &
     &         invTk*13847.26_r8-                                       &
     &         logTk*23.6521_r8-                                        &
     &         sqrtS*(5.977_r8-                                         &
     &                invTk*118.67_r8-                                  &
     &                logTk*1.0495_r8)-                                 &
     &         S(i)*0.01615_r8)
!
!------------------------------------------------------------------------
!  Compute salinity constant of hydrogen sulfate, Ks = [H][SO4]/[HSO4].
!  From Dickson (1990, J. chem. Thermodynamics 22, 113)
!------------------------------------------------------------------------
!
        Ks=EXP(141.328_r8-                                              &
     &         invTk*4276.1_r8-                                         &
     &         logTk*23.093_r8+                                         &
     &         sqrtSO4*(324.57_r8-invTk*13856.0_r8-logTk*47.986_r8-     &
     &                  SO4*invTk*2698.0_r8)-                           &
     &         SO4*(771.54_r8-invTk*35474.0_r8-logTk*114.723_r8-        &
     &              SO4*invTk*1776.0_r8)+                               &
     &         LOG(1.0_r8-0.001005_r8*S(i)))
!
!-----------------------------------------------------------------------
!  Compute stability constant of hydrogen fluorid, Kf = [H][F]/[HF].
!  From Dickson and Riley (1979) -- change pH scale to total.
!-----------------------------------------------------------------------
!
        Kf=EXP(-12.641_r8+                                              &
     &         invTk*1590.2_r8+                                         &
     &         sqrtSO4*1.525_r8+                                        &
     &         LOG(1.0_r8-0.001005_r8*S(i))+                            &
     &         LOG(1.0_r8+0.1400_r8*scl/(96.062_r8*Ks)))
!
!-----------------------------------------------------------------------
! Calculate concentrations for borate (Uppstrom, 1974), sulfate (Morris
! and Riley, 1966), and fluoride (Riley, 1965).
!-----------------------------------------------------------------------
!
        borate=0.000232_r8*scl/10.811_r8
        sulfate=0.14_r8*scl/96.062_r8
        fluoride=0.000067_r8*scl/18.9984_r8
!
!=======================================================================
!  Iteratively solver for computing hydrogen ions [H+] using either:
!
!    (1) Newton-Raphson method with fixed number of iterations,
!        use previous [H+] as first guess, or
!    (2) bracket and bisection
!=======================================================================
!
!  Set first guess and brackets for [H+] solvers.
!
        pH_guess=pH(i,j)         ! Newton-Raphson
        pH_hi=10.0_r8            ! high bracket/bisection
        pH_lo=5.0_r8             ! low bracket/bisection
!
!  Convert to [H+].
!
        X_guess=10.0_r8**(-pH_guess)
        X_lo=10.0_r8**(-pH_hi)
        X_hi=10.0_r8**(-pH_lo)
        X_mid=0.5_r8*(X_lo+X_hi)
!
!-----------------------------------------------------------------------
!  Newton-Raphson method.
!-----------------------------------------------------------------------
!
        IF (DoNewton.eq.1) THEN
          X=X_guess
          K12=K1*K2
          K12p=K1p*K2p
          K123p=K12p*K3p
          invKb=1.0_r8/Kb
          invKs=1.0_r8/Ks
          invKsi=1.0_r8/Ksi
!
          DO Inewton=1,InewtonMax
!
!  Set some common combinations of parameters used in the iterative [H+]
!  solver.
!
            X2=X*X
            X3=X2*X
            invX=1.0_r8/X
            invX2=1.0_r8/X2
            A=X*(K12p+X*(K1p+X))
            B=X*(K1+X)+K12
            C=1.0_r8/(1.0_r8+sulfate*invKs)
            A2=A*A
            B2=B*B
            dA=X*(2.0_r8*K1p+3.0_r8*X)+K12p
            dB=2.0_r8*X+K1
!
!  Evaluate f([H+]):
!
!     fn=HCO3+CO3+borate+OH+HPO4+2*PO4+H3PO4+silicate+Hfree+HSO4+HF-TALK
!
            fn=dic*K1*(X+2.0_r8*K2)/B+                                  &
     &         borate/(1.0_r8+X*invKb)+                                 &
     &         Kw*invX+                                                 &
     &         phos*(K12p*X+2.0_r8*K123p-X3)/A+                         &
     &         sili/(1.0_r8+X*invKsi)-                                  &
     &         X*C-                                                     &
     &         sulfate/(1.0_r8+Ks*invX*C)-                              &
     &         fluoride/(1.0_r8+Kf*invX)-                               &
     &         alk
!
!  Evaluate derivative, f(prime)([H+]):
!
!     df= d(fn)/d(X)
!
            df=dic*K1*(B-dB*(X+2.0_r8*K2))/B2-                          &
     &         borate/(invKb*(1.0+X*invKb)**2)-                         &
     &         Kw*invX2+                                                &
     &         phos*(A*(K12p-3.0_r8*X2)-dA*(K12p*X+2.0_r8*K123p-X3))/A2-&
     &         sili/(invKsi*(1.0_r8+X*invKsi)**2)+                      &
     &         C+                                                       &
     &         sulfate*Ks*C*invX2/((1.0_r8+Ks*invX*C)**2)+              &
     &         fluoride*Kf*invX2/((1.0_r8+Kf*invX)**2)
!
!  Evaluate increment in [H+].
!
            deltaX=-fn/df
!
!  Update estimate of [H+].
!
            X=X+deltaX
          END DO
!
!-----------------------------------------------------------------------
!  Bracket and bisection method.
!-----------------------------------------------------------------------
!
        ELSE
!
!  If first step, use Bracket and Bisection method with fixed, large
!  number of iterations
!
          K12=K1*K2
          K12p=K1p*K2p
          K123p=K12p*K3p
          invKb=1.0_r8/Kb
          invKs=1.0_r8/Ks
          invKsi=1.0_r8/Ksi
!
          BRACK_IT: DO Ibrack=1,IbrackMax
            DO Hstep=1,3
              IF (Hstep.eq.1) X=X_hi
              IF (Hstep.eq.2) X=X_lo
              IF (Hstep.eq.3) X=X_mid
!
!  Set some common combinations of parameters used in the iterative [H+]
!  solver.
!
              X2=X*X
              X3=X2*X
              invX=1.0_r8/X
              A=X*(K12p+X*(K1p+X))+K123p
              B=X*(K1+X)+K12
              C=1.0_r8/(1.0_r8+sulfate*invKs)
              A2=A*A
              B2=B*B
              dA=X*(K1p*2.0_r8+3.0_r8*X2)+K12p
              dB=2.0_r8*X+K1
!
!  Evaluate f([H+]) for bracketing and mid-value cases.
!
              fni(Hstep)=dic*(K1*X+2.0_r8*K12)/B+                       &
     &                   borate/(1.0_r8+X*invKb)+                       &
     &                   Kw*invX+                                       &
     &                   phos*(K12p*X+2.0_r8*K123p-X3)/A+               &
     &                   sili/(1.0_r8+X*invKsi)-                        &
     &                   X*C-                                           &
     &                   sulfate/(1.0_r8+Ks*invX*C)-                    &
     &                   fluoride/(1.0_r8+Kf*invX)-                     &
     &                   alk
            END DO
!
!  Now, bracket solution within two of three.
!
            IF (fni(3).eq.0.0_r8) THEN
              EXIT BRACK_IT
            ELSE
              ftest=fni(1)/fni(3)
              IF (ftest.gt.0.0) THEN
                X_hi=X_mid
              ELSE
                X_lo=X_mid
              END IF
              X_mid=0.5_r8*(X_lo+X_hi)
            END IF
          END DO BRACK_IT
!
! Last iteration gives value.
!
          X=X_mid
        END IF
!
!-----------------------------------------------------------------------
!  Determine pCO2.
!-----------------------------------------------------------------------
!
!  Total Hydrogen ion concentration, Htotal = [H+].
!
        Htotal=X
        Htotal2=Htotal*Htotal
!
!  Calculate [CO2*] (mole/m3) as defined in DOE Methods Handbook 1994
!  Version 2, ORNL/CDIAC-74, Dickson and Goyet, Eds. (Chapter 2,
!  page 10, Eq A.49).
!
        CO2star=dic*Htotal2/(Htotal2+K1*Htotal+K1*K2)
!
!  Save pH is used again outside this routine.
!
        pH(i,j)=-LOG10(Htotal)
!
!  Add two output arguments for storing pCO2surf.
!
        pCO2(i)=CO2star*1000000.0_r8/ff
      END DO I_LOOP
      RETURN
      END SUBROUTINE pCO2_water
      END MODULE biology_mod
