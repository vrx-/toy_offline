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
 * ** Application flag:   UPWELLING_OFFLINE
 * ** Input script:       ocean_upwelling_offline.in
 * */

#define ROMS_MODEL

#define UV_LDRAG
#define SOLVE3D
#define AVERAGES

#define UPWELLING
#define ANA_GRID
#undef ANA_INITIAL
#undef ANA_SMFLUX
#define ANA_SPFLUX
#define ANA_BPFLUX

#define OFFLINE
#define BIOLOGY
#define OFFLINE_BIOLOGY
#define BIO_FENNEL
#undef ANA_BIOLOGY

#define NONLIN_EOS
#define SOLAR_SOURCE
#define DIURNAL_SRFLUX

#define CARBON
#define DENITRIFICATION
#define BIO_SEDIMENT
#define DIAGNOSTICS_BIO

#define OCLIMATOLOGY
#define USE_DEBUG
