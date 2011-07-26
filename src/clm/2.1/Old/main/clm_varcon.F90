#include <misc.h>
#include <preproc.h>

module clm_varcon

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clm_varcon
!
! !DESCRIPTION: 
! Module containing various model constants 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_G,SHR_CONST_STEBOL,SHR_CONST_KARMAN,     &
       SHR_CONST_RWV,SHR_CONST_RDAIR,SHR_CONST_CPFW,      &
       SHR_CONST_CPICE,SHR_CONST_CPDAIR,SHR_CONST_LATVAP, &
       SHR_CONST_LATSUB,SHR_CONST_LATICE,SHR_CONST_RHOFW, &
       SHR_CONST_RHOICE,SHR_CONST_TKFRZ,SHR_CONST_REARTH
  use clm_varpar, only : numcol, numrad, nlevsoi, nlevlak
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!----------------------------------------------------------------------- 

  !------------------------------------------------------------------
  ! Initialize physical constants
  !------------------------------------------------------------------

  real(r8) :: grav   = SHR_CONST_G      !gravity constant [m/s2]
  real(r8) :: sb     = SHR_CONST_STEBOL !stefan-boltzmann constant  [W/m2/K4]
  real(r8) :: vkc    = SHR_CONST_KARMAN !von Karman constant [-]
  real(r8) :: rwat   = SHR_CONST_RWV    !gas constant for water vapor [J/(kg K)]
  real(r8) :: rair   = SHR_CONST_RDAIR  !gas constant for dry air [J/kg/K]
  real(r8) :: roverg = SHR_CONST_RWV/SHR_CONST_G*1000. !Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
  real(r8) :: cpliq  = SHR_CONST_CPFW   !Specific heat of water [J/kg-K]
  real(r8) :: cpice  = SHR_CONST_CPICE  !Specific heat of ice [J/kg-K]
  real(r8) :: cpair  = SHR_CONST_CPDAIR !specific heat of dry air [J/kg/K]
  real(r8) :: hvap   = SHR_CONST_LATVAP !Latent heat of evap for water [J/kg]
  real(r8) :: hsub   = SHR_CONST_LATSUB !Latent heat of sublimation    [J/kg]
  real(r8) :: hfus   = SHR_CONST_LATICE !Latent heat of fusion for ice [J/kg]
  real(r8) :: denh2o = SHR_CONST_RHOFW  !density of liquid water [kg/m3]
  real(r8) :: denice = SHR_CONST_RHOICE !density of ice [kg/m3]
  real(r8) :: tkair  = 0.023     !thermal conductivity of air   [W/m/k]
  real(r8) :: tkice  = 2.290     !thermal conductivity of ice   [W/m/k]
  real(r8) :: tkwat  = 0.6       !thermal conductivity of water [W/m/k]
  real(r8) :: tfrz   = SHR_CONST_TKFRZ  !freezing temperature [K]
  real(r8) :: tcrit  = 2.5       !critical temperature to determine rain or snow
  real(r8) :: po2    = 0.209     !constant atmospheric partial pressure  O2 (mol/mol)
  real(r8) :: pco2   = 355.e-06  !constant atmospheric partial pressure CO2 (mol/mol)

  real(r8) :: bdsno = 250.       !bulk density snow (kg/m**3)

  real(r8) :: re = SHR_CONST_REARTH*0.001 !radius of earth (km)

  real(r8), public, parameter :: spval = 1.e36  !special value for missing data (ocean)

  ! These are tunable constants from clm2_3

  real(r8) :: zlnd = 0.01      !Roughness length for soil [m]
  real(r8) :: zsno = 0.0024    !Roughness length for snow [m]
  real(r8) :: csoilc = 0.004   !Drag coefficient for soil under canopy [-]
  real(r8) :: capr   = 0.34    !Tuning factor to turn first layer T into surface T  
  real(r8) :: cnfac  = 0.5     !Crank Nicholson factor between 0 and 1
  real(r8) :: ssi    = 0.033   !Irreducible water saturation of snow
  real(r8) :: wimp   = 0.05    !Water impremeable if porosity less than wimp
  real(r8) :: pondmx = 10.0    !Ponding depth (mm)


  !------------------------------------------------------------------
  ! Initialize water type constants
  !------------------------------------------------------------------

  ! "water" types 
  !   1     soil
  !   2     land ice (glacier)
  !   3     deep lake
  !   4     shallow lake
  !   5     wetland: swamp, marsh, etc

  integer :: istsoil = 1  !soil         "water" type
  integer :: istice  = 2  !land ice     "water" type
  integer :: istdlak = 3  !deep lake    "water" type
  integer :: istslak = 4  !shallow lake "water" type
  integer :: istwet  = 5  !wetland      "water" type

  !------------------------------------------------------------------
  ! Initialize miscellaneous radiation constants
  !------------------------------------------------------------------

  integer, private :: i  ! loop index

  ! saturated soil albedos for 8 color classes: 1=vis, 2=nir

  real(r8) :: albsat(numcol,numrad) !wet soil albedo by color class and waveband
  data(albsat(i,1),i=1,8)/0.12,0.11,0.10,0.09,0.08,0.07,0.06,0.05/
  data(albsat(i,2),i=1,8)/0.24,0.22,0.20,0.18,0.16,0.14,0.12,0.10/

  ! dry soil albedos for 8 color classes: 1=vis, 2=nir 

  real(r8) :: albdry(numcol,numrad) !dry soil albedo by color class and waveband
  data(albdry(i,1),i=1,8)/0.24,0.22,0.20,0.18,0.16,0.14,0.12,0.10/
  data(albdry(i,2),i=1,8)/0.48,0.44,0.40,0.36,0.32,0.28,0.24,0.20/

  ! albedo land ice: 1=vis, 2=nir

  real(r8) :: albice(numrad)        !albedo land ice by waveband
  data (albice(i),i=1,numrad) /0.80, 0.55/

  ! albedo frozen lakes: 1=vis, 2=nir 

  real(r8) :: alblak(numrad)        !albedo frozen lakes by waveband
  data (alblak(i),i=1,numrad) /0.60, 0.40/

  ! omega,betad,betai for snow 

  real(r8) :: betads  = 0.5       !two-stream parameter betad for snow
  real(r8) :: betais  = 0.5       !two-stream parameter betai for snow
  real(r8) :: omegas(numrad)      !two-stream parameter omega for snow by band
  data (omegas(i),i=1,numrad) /0.8, 0.4/

  !------------------------------------------------------------------
  ! Soil and Lake depths are constants for now
  ! The values for the following arrays are set in routine iniTimeConst
  !------------------------------------------------------------------

  real(r8) :: zlak(1:nlevlak)     !lake z  (layers) 
  real(r8) :: dzlak(1:nlevlak)    !lake dz (thickness)
  real(r8) :: zsoi(1:nlevsoi)     !soil z  (layers)
  real(r8) :: dzsoi(1:nlevsoi)    !soil dz (thickness)
  real(r8) :: zisoi(0:nlevsoi)    !soil zi (interfaces)  

end module clm_varcon
