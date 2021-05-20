!##############################################################################
! MSIS® (NRL-SOF-014-1) SOFTWARE
!
! MSIS® is a registered trademark of the Government of the United States of 
! America, as represented by the Secretary of the Navy. Unauthorized use of 
! the trademark is prohibited. 
!
! The MSIS® Software (hereinafter Software) is property of the United States 
! Government, as represented by the Secretary of the Navy. Methods performed
! by this software are covered by U.S. Patent Number 10,641,925. The Government
! of the United States of America, as represented by the Secretary of the Navy, 
! herein grants a non-exclusive, non-transferable license to the Software for 
! academic, non-commercial, purposes only. A user of the Software shall not: 
! (i) use the Software for any non-academic, commercial purposes, (ii) make 
! any modification or improvement to the Software, (iii) disseminate the 
! Software or any supporting data to any other person or entity who will use 
! the Software for any non-academic, commercial purposes, or (iv) copy the 
! Software or any documentation related thereto except for (a) distribution 
! among the user’s personal computer systems, archival, or emergency repair 
! purposes, or (b) distribution for non-commercial, academic purposes, without 
! first obtaining the written consent of IP Counsel for the Naval Research 
! Laboratory. 
!
! As the owner of MSIS®, the United States, the United States Department of 
! Defense, and their employees: (1) Disclaim any warranties, express, or 
! implied, including but not limited to any implied warranties of 
! merchantability, fitness for a particular purpose, title or non-infringement, 
! (2) Do not assume any legal liability or responsibility for the accuracy, 
! completeness, or usefulness of the software, (3) Do not represent that use of 
! the software would not infringe privately owned rights, (4) Do not warrant 
! that the software will function uninterrupted, that is error-free or that any 
! errors will be corrected.
!
! BY USING THIS SOFTWARE YOU ARE AGREEING TO THE ABOVE TERMS AND CONDITIONS.  
!##############################################################################

!!! ===========================================================================
!!! NRLMSIS 2.0:
!!! Neutral atmosphere empirical model from the surface to lower exosphere
!!! John Emmert (john.emmert@nrl.navy.mil)
!!! Doug Drob (douglas.drob@nrl.navy.mil)
!!! ===========================================================================

!**************************************************************************************************
! MSIS_CONSTANTS Module: Contains constants and hardwired parameters
!**************************************************************************************************
module msis_constants

  implicit none

  ! Floating Point Precision
#ifdef DBLE
  integer, parameter         :: rp = 8
#else
  integer, parameter         :: rp = 4
#endif

  ! Missing density value
  real(kind=rp),parameter    :: dmissing = 9.999e-38_rp

  ! Trigonometric constants
  real(kind=rp), parameter   :: pi = 3.1415926535897932384626433832795_rp
  real(kind=rp), parameter   :: deg2rad = pi / 180.0_rp
  real(kind=rp), parameter   :: doy2rad = 2.0_rp*pi / 365.0_rp
  real(kind=rp), parameter   :: lst2rad = pi / 12.0_rp
  !real(kind=rp), parameter   :: tanh1 = 0.761594155955765485_rp  ! tanh(1.0)
  real(kind=rp), parameter   :: tanh1 = tanh(1.0_rp)

  ! Thermodynamic constants
  ! Boltzmann constant (CODATA 2018) (J/kg)
  real(kind=rp), parameter   :: kB = 1.380649e-23_rp
  ! Avogadro constant (CODATA 2018)
  real(kind=rp), parameter   :: NA = 6.02214076e23_rp
  ! Reference gravity (CIMO Guide 2014) (m/s^2) (specified separately in alt2gph.f90)
  real(kind=rp), parameter   :: g0 = 9.80665_rp
  ! Species molecular masses (kg/molecule) (CIPM 2007)
  real(kind=rp), parameter   :: specmass(1:10) = (/  0.0_rp,                          & ! Mass density (dummy value)
                                                    28.0134_rp,                       & ! N2
                                                    31.9988_rp,                       & ! O2
                                                    31.9988_rp/2.0_rp,                & ! O
                                                     4.0_rp,                          & ! He
                                                     1.0_rp,                          & ! H
                                                    39.948_rp,                        & ! Ar
                                                    28.0134_rp/2.0_rp,                & ! N
                                                    31.9988_rp/2.0_rp,                & ! Anomalous O
                                                    (28.0134_rp+31.9988_rp)/2.0_rp /) & ! NO
                                                    / (1.0e3_rp * NA)                   ! Convert from g/mol to kg/molecule
  ! Dry air mean mass in fully mixed atmosphere (CIPM 2007) (includes CO2 and other trace species that are not yet in MSIS)
  real(kind=rp), parameter   :: Mbar = 28.96546_rp / (1.0e3_rp * NA)   ! kg/molecule
  ! Dry air log volume mixing ratios (CIPM 2007)
  real(kind=rp), parameter   :: lnvmr(1:10) = log( (/ 1.0_rp,        & ! Mass density (dummy value)
                                                      0.780848_rp,   & ! N2
                                                      0.209390_rp,   & ! O2
                                                      1.0_rp,        & ! O (dummy value)
                                                      0.0000052_rp,  & ! He
                                                      1.0_rp,        & ! H (dummy value)
                                                      0.009332_rp,   & ! Ar
                                                      1.0_rp,        & ! N (dummy value)
                                                      1.0_rp,        & ! Anomalous O (dummy value)
                                                      1.0_rp /) )      ! NO (dummy value)
  ! Natural log of global average surface pressure (Pa)
  !real(kind=rp), parameter   :: lnP0 = 11.5080482 !+ 0.00759597 After calibration with MERRA2
  real(kind=rp), parameter   :: lnP0 = 11.515614
  ! Derived constants
  real(kind=rp), parameter   :: g0divkB = g0/kB * 1.0e3_rp  ! K/(kg km)
  real(kind=rp), parameter   :: Mbarg0divkB = Mbar*g0/kB * 1.0e3_rp   ! K/km
  ! References:
  ! CODATA Internationally recommended 2018 values of the fundamental physical constants.
  !   https://pml.nist.gov/cuu/Constants/; https://pml.nist.gov/cuu/pdf/wallet_2018.pdf
  ! Picard, A., Davis, R. S., Glaeser, M., and Fujii, K. (2007). Revised formula for the density of
  !   air (CIPM 2007). Metrologia 45, 149–155. doi:10.1088/0026-1394/45/2/004
  ! World Meteorological Organization (2014). WMO guide to meteorological instruments and methods of observation
  !   (the CIMO Guide). Part I, Chapter 12. https://www.wmo.int/pages/prog/www/IMOP/CIMO-Guide.html

  ! Vertical profile parameters
  integer, parameter         :: nspec = 11  !Number of species including temperature
  integer, parameter         :: nd = 27     !Number of temperature profile nodes
  integer, parameter         :: p = 4       !Spline order
  integer, parameter         :: nl = nd - p !Last temperature profile level index
  integer, parameter         :: nls = 9     !Last parameter index for each species (excluding O, NO splines)
  real(kind=rp), parameter   :: bwalt = 122.5_rp ! Reference geopotential height for Bates Profile
  real(kind=rp), parameter   :: zetaF = 70.0_rp  ! Fully mixed below this, uses constant mixing ratios
  real(kind=rp), parameter   :: zetaB = bwalt    ! Bates Profile above this altitude
  real(kind=rp), parameter   :: zetaA = 85.0_rp  ! Default reference height for active minor species
  real(kind=rp), parameter   :: zetagamma = 100.0_rp  ! Reference height of tanh taper of chemical/dynamical correction scale height
  real(kind=rp), parameter   :: Hgamma = 1.0_rp/30.0_rp  ! Inverse scale height of tanh taper of chemical/dynamical correction scale height
  real(kind=rp), parameter   :: nodesTN(0:nd+2) = &  !Nodes for temperature profile splines
      (/ -15., -10.,  -5.,   0.,   5., 10., 15., 20.,  25.,  30.,  35.,  40., 45., 50., &
          55.,  60.,  65.,  70.,  75., 80., 85., 92.5, 102.5, 112.5, 122.5, 132.5, 142.5, &
          152.5, 162.5, 172.5/)
  integer, parameter         :: izfmx = 13       ! fully mixed below this spline index
  integer, parameter         :: izfx = 14        ! Spline index at zetaF
  integer, parameter         :: izax = 17        ! Spline index at zetaA
  integer, parameter         :: itex = nl        ! Index of Bates exospheric temperature
  integer, parameter         :: itgb0 = nl - 1   ! Index of Bates temperature gradient at lower boundary
  integer, parameter         :: itb0 = nl - 2    ! Index of Bates temperature at lower boundary
  ! O1 Spline parameters
  integer, parameter         :: ndO1 = 13
  integer, parameter         :: nsplO1 = ndO1-5     !Number of unconstrained spline parameters for O1 (there are 2 additional C1-constrained splines)
  real(kind=rp), parameter   :: nodesO1(0:ndO1) = & !Nodes for O1 splines (Domain 50-85 km)
      (/ 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 92.5, 102.5, 112.5/)
  real(kind=rp), parameter   :: zetarefO1 = zetaA   !Joining height for O1 splines, and reference height for O1 density
  ! NO Spline parameters
  integer, parameter         :: ndNO = 13
  integer, parameter         :: nsplNO = ndNO-5     !Number of unconstrained spline parameters for NO (there are 2 additional C1-constrained splines)
  real(kind=rp), parameter   :: nodesNO(0:ndNO) = & !Nodes for NO splines (Domain 70-122.5 km)
      (/ 47.5, 55., 62.5, 70., 77.5, 85., 92.5, 100., 107.5, 115., 122.5, 130., 137.5, 145./)
  real(kind=rp), parameter   :: zetarefNO = zetaB   !Joining height for NO splines, and reference height for NO density
  !C2 Continuity matrix for temperature; Last 3 splines are constrained (must be recomputed if nodes change)
  real(kind=rp), parameter   :: c2tn(3,3) = reshape((/1.0_rp, -10.0_rp,  33.333333333333336_rp, &
                                                      1.0_rp,   0.0_rp, -16.666666666666668_rp, &
                                                      1.0_rp,  10.0_rp,  33.333333333333336_rp/), &
                                                    (/3,3/))
  !C1 Continuity for O1; Last 2 splines are constrained (must be recomputed if nodes change)
  real(kind=rp), parameter   :: c1o1(2,2) = reshape((/ 1.75_rp,               -2.916666573405061_rp, &
                                                      -1.624999900076852_rp,  21.458332647194382_rp /), &
                                                    (/2,2/))
  real(kind=rp), parameter   :: c1o1adj(2) = (/0.257142857142857_rp, -0.102857142686844_rp/) !Weights for coefficents on 3rd to last spline; product to be subtracted from RHS of continuity equation
  !C1 Continuity for NO; Last 2 splines are constrained (must be recomputed if nodes change)
  real(kind=rp), parameter   :: c1NO(2,2) = reshape((/ 1.5_rp, -3.75_rp, &
                                                       0.0_rp,  15.0_rp /), &
                                                    (/2,2/))
  real(kind=rp), parameter   :: c1NOadj(2) = (/0.166666666666667_rp, -0.066666666666667_rp/) !Weights for coefficents on 3rd to last spline; product to be subtracted from RHS of continuity equation
  ! Anomalous Oxygen parameters (legacy profile from NRLMSISE-00)
  real(kind=rp),parameter    :: zetarefOA = zetaB   !Reference height for anomalous oxygen density
  real(kind=rp),parameter    :: TOA = 4000.         !Temperature of anomalous oxygen density (K)
  real(kind=rp),parameter    :: HOA = (kB * TOA) / ( (16.0_rp/(1.0e3_rp*NA)) * g0 ) * 1.0e-3_rp  !Hydrostatic scale height of anomalous oxygen density (km)
    
  ! Horizontal and time-dependent basis function (gfn) parameters
  integer, parameter      :: maxnbf = 512   ! Number of basis functions to be allocated
  integer, parameter      :: maxn = 6       ! Maximum latitude (Legendre) spectral degree
  integer, parameter      :: maxl = 3       ! Maximum local time (tidal) spectral order
  integer, parameter      :: maxm = 2       ! Maximum longitude (stationary planetary wave) order
  integer, parameter      :: maxs = 2       ! Maximimum day of year (intra-annual) Fourier order
  integer, parameter      :: amaxn = 6      ! Maximum Legendre degree used in time independent and intra-annual zonal mean terms
  integer, parameter      :: amaxs = 2      ! Maximum intra-annual order used in zonal mean terms
  integer, parameter      :: tmaxl = 3      ! Maximum tidal order used
  integer, parameter      :: tmaxn = 6      ! Maximum Legendre degree coupled with tides
  integer, parameter      :: tmaxs = 2      ! Maximum intra-annual order coupled with tides
  integer, parameter      :: pmaxm = 2      ! Maximum stationary planetary wave order used
  integer, parameter      :: pmaxn = 6      ! Maximum Legendre degree coupled with SPW
  integer, parameter      :: pmaxs = 2      ! Maximum intra-annual order coupled with SPW
  integer, parameter      :: nsfx = 5       ! Number of linear solar flux terms
  integer, parameter      :: nsfxmod = 5    ! Number of nonlinear modulating solar flux terms (legacy NRLMSISE-00 terms)
  integer, parameter      :: nmag = 54      ! Number of terms in NRLMSISE-00 legacy geomagnetic parameterization
  integer, parameter      :: nut = 12       ! Number of terms in NRLMSISE-00 legacy UT parameterization
  integer, parameter      :: ctimeind = 0             ! Starting index of time-independent terms
  integer, parameter      :: cintann = ctimeind + (amaxn+1)   ! Starting index of zonal mean intra-annual terms
  integer, parameter      :: ctide = cintann + ((amaxn+1)*2*amaxs)   ! Starting index of zonal mean intra-annual terms
  integer, parameter      :: cspw = ctide + (4*tmaxs+2)*(tmaxl*(tmaxn+1)-(tmaxl*(tmaxl+1))/2) ! Starting index of SPW terms
  integer, parameter      :: csfx = cspw + (4*pmaxs+2)*(pmaxm*(pmaxn+1)-(pmaxm*(pmaxm+1))/2)   ! Starting index of linear solar flux terms
  integer, parameter      :: cextra = csfx + nsfx     ! Starting index of time-independent terms
  integer, parameter      :: mbf = 383                ! Last index of linear terms
  integer, parameter      :: cnonlin = mbf + 1        ! Starting index of nonlinear terms
  integer, parameter      :: csfxmod = cnonlin        ! Starting index of modulating solar flux terms
  integer, parameter      :: cmag = csfxmod + nsfxmod ! Starting index of daily geomagnetic terms
  integer, parameter      :: cut = cmag + nmag        ! Starting index of UT terms
    
  ! Weights for calculation log pressure spline coefficients from temperature coefficients (must be recalcuated if nodes change)
  real(kind=rp), parameter   :: gwht(0:3) =  (/ 5.0_rp/24.0_rp, 55.0_rp/24.0_rp, 55.0_rp/24.0_rp, 5.0_rp/24.0_rp /)

  ! Constants needed for analytical integration by parts of hydrostatic piecewise effective mass profile
  real(kind=rp), parameter   :: wbeta(0:nl) =  (nodesTN(4:nd)  - nodesTN(0:nl)) / 4.0_rp !Weights for 1st spline integration
  real(kind=rp), parameter   :: wgamma(0:nl) = (nodesTN(5:nd+1)- nodesTN(0:nl)) / 5.0_rp !Weights for 2nd spline integration
  ! Non-zero bspline values at zetaB (5th and 6th order) (must be recalcuated if nodes change)
  real(kind=rp), parameter   :: S5zetaB(0:3) = (/0.041666666666667_rp, 0.458333333333333_rp, 0.458333333333333_rp, &
                                                 0.041666666666667_rp/)
  real(kind=rp), parameter   :: S6zetaB(0:4) = (/0.008771929824561_rp, 0.216228070175439_rp, 0.550000000000000_rp, &
                                                 0.216666666666667_rp, 0.008333333333333_rp/)
  !Weights for calculating temperature gradient at zetaA (must be recalcuated if nodes change)
  real(kind=rp), parameter   :: wghtAxdz(0:2) = (/-0.102857142857_rp, 0.0495238095238_rp, 0.053333333333_rp/)
  !Non-zero bspline values at zetaA (4th, 5th and 6th order) (must be recalcuated if nodes change)
  real(kind=rp), parameter   :: S4zetaA(0:2) = (/0.257142857142857_rp, 0.653968253968254_rp, 0.088888888888889_rp/)
  real(kind=rp), parameter   :: S5zetaA(0:3) = (/0.085714285714286_rp, 0.587590187590188_rp, 0.313020313020313_rp, &
                                                 0.013675213675214_rp/)
  real(kind=rp), parameter   :: S6zetaA(0:4) = (/0.023376623376623_rp, 0.378732378732379_rp, 0.500743700743701_rp, &
                                                 0.095538448479625_rp, 0.001608848667672_rp/)
  !Non-zero bspline values at zetaF (4th and 5th order) (must be recalcuated if nodes change)
  real(kind=rp), parameter   :: S4zetaF(0:2) = (/0.166666666666667_rp, 0.666666666666667_rp, 0.166666666666667_rp/)
  real(kind=rp), parameter   :: S5zetaF(0:3) = (/0.041666666666667_rp, 0.458333333333333_rp, 0.458333333333333_rp, &
                                                 0.041666666666667_rp/)
  !Non-zero bspline values at zeta=0 (5th order) (must be recalcuated if nodes change)
  real(kind=rp), parameter   :: S5zeta0(0:2) = (/0.458333333333333_rp, 0.458333333333333_rp, 0.041666666666667_rp/)

end module msis_constants
