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

!==================================================================================================
! ALT2GPH: Altitude to Geopotential Height
! References:
!   DMA Technical Report TR8350.2 (1987),
!     http://earth-info.nga.mil/GandG/publications/historic/historic.html
!   Featherstone, W. E., and S. J. Claessens (2008), Closed-form transformation between
!     geodetic and ellipsoidal coordinates, Studia Geophysica et Geodaetica, 52, 1-18
!   Jekeli, C. (2009), Potential theory and static gravity field of the Earth, in 
!     Treatise on Geophysics, ed. T. Herring, vol 3, 11-42
!   NIMA Technical Report TR8350.2 (2000, 3rd edition, Amendment1), 
!     http://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350_2.html
!==================================================================================================
real(8) function alt2gph(lat,alt)

  implicit none

  ! Input variables
  real(8), intent(in) :: lat    !Geodetic latitude (deg)
  real(8), intent(in) :: alt    !Geodetic altitude (km)

  real(8), parameter  :: deg2rad = 0.017453292519943295d0

  ! WGS84 Defining parameters
  real(8), parameter  :: a = 6378.1370d0 * 1d3 !Semi-major axis of reference ellipsoid (m)
  real(8), parameter  :: finv = 298.257223563d0 ! 1/f = Reciprocal of flattening
  real(8), parameter  :: w = 7292115d-11 !Angular velocity of Earth rotation (rad/s)
  real(8), parameter  :: GM = 398600.4418 * 1d9 !Gravitational constant x Earth mass (m^3/s^2)

  ! WGS84 Derived parameters
  real(8), parameter  :: asq = a*a
  real(8), parameter  :: wsq = w*w
  real(8), parameter  :: f = 1.0d0 / finv
  real(8), parameter  :: esq = 2*f - f*f
  real(8), parameter  :: e = sqrt(esq)  !Ellipsoid eccentricity
  real(8), parameter  :: Elin = a*e     !Linear eccentricity of ellipsoid
  real(8), parameter  :: Elinsq = Elin*Elin
  real(8), parameter  :: epr = e / (1-f)  !Second eccentricity
  real(8), parameter  :: q0 = ((1.0d0 + 3.0d0/(epr*epr))*atan(epr) - 3.0d0/epr)/2.0d0  !DMA Technical Report tr8350.2, Eq. 3-25
  real(8), parameter  :: U0 = -GM*atan(epr)/Elin - wsq*asq/3d0 !Theoretical potential of reference ellipsoid (m^2/s^2), DMA Technical Report tr8350.2, Eq. 3-51
  real(8), parameter  :: g0 = 9.80665d0 !Standard gravity (m/s^2), CGPM 1901; WMO
  real(8), parameter  :: GMdivElin = GM / Elin
  
  ! Parameters for centrifugal potential taper
  real(8), parameter  :: x0sq = 2d7**2   !Axial distance squared at which tapering begins (m^2)
  real(8), parameter  :: Hsq = 1.2d7**2  !Relaxation scale length of taper (m^2)

  ! Working variables
  real(8)             :: altm, sinsqlat, v, xsq, zsq
  real(8)             :: rsqminElinsq, usq, cossqdelta, epru, atanepru, q, U, Vc

  ! Compute Cartesian and ellipsoidal coordinates
  altm = alt * 1000.0d0
  sinsqlat = sin(lat*deg2rad)**2
  v = a / sqrt(1-esq*sinsqlat)           !Radius of curvature of the reference ellipsoid, Featherstone eq. 4
  xsq = (v + altm)**2 * (1 - sinsqlat)   !Squared x-coordinate of geocentric system, Featherstone eq. 1
  zsq = (v*(1-esq) + altm)**2 * sinsqlat !Squared z-coordinate of geocentric system, Featherstone eq. 3
  rsqminElinsq = xsq + zsq - Elinsq
  usq = rsqminElinsq/2.0d0 + sqrt(rsqminElinsq**2 / 4.0d0 + Elinsq*zsq)  !Ellipsoidal distance coordinate, Featherstone eq. 19 
  cossqdelta = zsq / usq                 !Ellipsoidal polar angle, Featherstone eq. 21

  ! Compute gravitational potential
  epru = Elin / sqrt(usq)                !Second eccentricity at ellipsoidal coordinate u
  atanepru = atan(epru)
  q = ((1+3.0d0/(epru*epru))*atanepru - 3.0d0/epru)/2.0d0   !Jekeli, eq. 114
  U = -GMdivElin * atanepru - wsq * ( asq * q * (cossqdelta - 1/3.0d0) / q0 ) / 2.0d0   !Jekeli, eq. 113

  ! Compute centrifugal potential and adjust total potential
  if (xsq .le. x0sq) then
    Vc = (wsq/2.0d0) * xsq
  else
    Vc = (wsq/2.0d0) * (Hsq*tanh((xsq-x0sq)/Hsq) + x0sq) !Centrifugal potential taper
  endif
  U = U - Vc
  
  ! Compute geopotential height
  alt2gph = (U - U0) / g0 / 1000.0d0

  return

end function alt2gph

!==================================================================================================
! GPH2ALT: Geopotential Height to Altitude
!==================================================================================================
real(8) function gph2alt(theta,gph)

  implicit none

  real(8), external    :: alt2gph

  real(8), intent(in)  :: theta
  real(8), intent(in)  :: gph

  integer, parameter   :: maxn = 10
  real(8), parameter   :: epsilon = 0.0005

  real(8)              :: x,dx,y,dydz
  integer              :: n

  x = gph
  n = 0
  dx = epsilon + epsilon
  do while ((abs(dx) .gt. epsilon) .and. (n .lt. 10))
    y = alt2gph(theta,x)
    dydz = (alt2gph(theta,x+dx) - y)/dx
    dx = (gph - y)/dydz
    x = x + dx
    n = n + 1
  end do

  gph2alt = x

end function gph2alt
    
