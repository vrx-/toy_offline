/*
 * ** svn $Id$
 * *******************************************************************************
 * ** Copyright (c) 2002-2018 The ROMS/TOMS Group                               **
 * **   Licensed under a MIT/X style license                                    **
 * **   See License_ROMS.txt                                                    **
 * *******************************************************************************
 * **
 * ** Options for Upwelling Test.
 * **
 * ** Application flag:   UPWELLING
 * ** Input script:       ocean_upwelling.in
 * */

#define ROMS_MODEL

#define UV_ADV
#define UV_COR
#define UV_LDRAG
#define UV_VIS2
#define MIX_S_UV
#define SPLINES_VVISC
#define SPLINES_VDIFF
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define TS_DIF2
#define DJ_GRADPS
#define MIX_S_TS

#define SALINITY
#define SOLVE3D
#define AVERAGES

#define ANA_GRID
#define ANA_INITIAL
#define ANA_VMIX
#define ANA_SMFLUX
#define ANA_STFLUX

#define NONLIN_EOS
#define SOLAR_SOURCE
#define DIURNAL_SRFLUX

#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX

#undef BIOLOGY

#if defined BIOLOGY
# define BIO_FENNEL
# undef ANA_BIOLOGY
# define CARBON
# define DENITRIFICATION
# define BIO_SEDIMENT
#endif

#define ANA_SPFLUX
#define ANA_BPFLUX
#define ANA_SRFLUX
#define ANA_CLOUD
