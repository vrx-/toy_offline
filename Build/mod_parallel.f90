      MODULE mod_parallel
!
!svn $Id: mod_parallel.F 875 2017-11-03 01:10:02Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains all variables used for parallelization         !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Switch to identify master processor. In serial and shared-memory
!  applications it is always true.
!
      logical :: Master
!
!  Switch to identify which thread is processing input/output files.
!  In distributed-memory applications, this thread can be the master
!  thread or all threads in case of parallel output. In serial and
!  shared-memory applications it is always true.
!
      logical :: InpThread
      logical :: OutThread
!$OMP THREADPRIVATE (InpThread, OutThread)
!
!  Number of shared-memory parallel threads or distributed-memory
!  parallel nodes.
!
      integer :: numthreads
!
!  First and last tile to process in a tiled application.
!
      integer, allocatable :: first_tile(:)
      integer, allocatable :: last_tile(:)
!$OMP THREADPRIVATE (first_tile, last_tile)
!
!  Parallel nodes assined to the ocean model.
!
      integer :: peOCN_frst            ! first ocean parallel node
      integer :: peOCN_last            ! last  ocean parallel node
!
!  Parallel threads/nodes counters used in critical parallel regions.
!
      integer :: tile_count = 0
      integer :: block_count = 0
      integer :: thread_count = 0
!
!  Profiling variables as function of parallel thread:
!
!    proc          Parallel process ID.
!    Cstr          Starting time for program region.
!    Cend          Ending time for program region.
!    Csum          Accumulated time for progam region.
!    Ctotal        Total time profiled
!    total_cpu     Total elapsed CPU
!
      integer, allocatable :: proc(:,:,:)
      real(r8) :: Ctotal, total_cpu, total_model(4)
      real(r8), allocatable :: Cstr(:,:,:)
      real(r8), allocatable :: Cend(:,:,:)
      real(r8), allocatable :: Csum(:,:,:)
!$OMP THREADPRIVATE (proc)
!$OMP THREADPRIVATE (Cstr, Cend)
!
!  Distributed-memory master process.
!
      integer, parameter :: MyMaster = 0
!
!  Rank of the parallel local process.
!
      integer :: MyRank = 0
      integer :: MyThread = 0
!$OMP THREADPRIVATE (MyThread)
!
      CONTAINS
!
      SUBROUTINE allocate_parallel
!
!=======================================================================
!                                                                      !
!  This routine allocates several variables in the module that depend  !
!  on the number of nested grids.                                      !
!                                                                      !
!=======================================================================
!
      USE mod_strings, ONLY: Nregion
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL
!  First and last tile to process in a tiled application.
!
      allocate ( first_tile(Ngrids) )
      allocate ( last_tile (Ngrids) )
!
!  Time profiling variables.
!
      allocate ( proc(0:1,4,Ngrids) )
      proc(0:1,1:4,1:Ngrids)=0
      allocate ( Cstr(0:Nregion,4,Ngrids) )
      Cstr(0:Nregion,1:4,1:Ngrids)=0.0_r8
      allocate ( Cend(0:Nregion,4,Ngrids) )
      Cend(0:Nregion,1:4,1:Ngrids)=0.0_r8
!$OMP END PARALLEL
      allocate ( Csum(0:Nregion,4,Ngrids) )
      Csum(0:Nregion,1:4,1:Ngrids)=0.0_r8
!
! Initialize other profiling variables.
!
      Ctotal=0.0_r8
      total_cpu=0.0_r8
      total_model=0.0_r8
      RETURN
      END SUBROUTINE allocate_parallel
      SUBROUTINE initialize_parallel
!
!=======================================================================
!                                                                      !
!  This routine initializes and spawn distribute-memory nodes.         !
!                                                                      !
!=======================================================================
!
      USE mod_iounits
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: i
      integer :: my_numthreads, my_threadnum
!
!-----------------------------------------------------------------------
!  Initialize serial configuration.
!-----------------------------------------------------------------------
!
      numthreads=my_numthreads()
      Master=.TRUE.
      InpThread=.TRUE.
      OutThread=.TRUE.
      RETURN
      END SUBROUTINE initialize_parallel
      END MODULE mod_parallel
