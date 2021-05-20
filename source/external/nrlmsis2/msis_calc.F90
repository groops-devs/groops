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
!!!
!!! MSISCALC: New interface with re-ordered input arguments and output arrays.
!
!     PREREQUISITES:
!       Must first run MSISINIT to load parameters and set switches. The 
!       MSISCALC subroutine checks for initialization and does a default
!       initialization if necessary. This self-initialization will be removed
!       in future versions.
!
!     CALLING SEQUENCE:
!       CALL MSISCALC(DAY, UTSEC, Z, LAT, LON, SFLUXAVG, SFLUX, AP, TN, DN, [TEX])
!  
!     INPUT VARIABLES:
!       DAY       Day of year (1.0 to 365.0 or 366.0)
!       UTSEC     Universal time (seconds)
!       Z         Geodetic altitude (km) (default) or Geopotential height (km)
!       LAT       Geodetic latitude (deg)
!       LON       Geodetic longitude (deg)
!       SFLUXAVG  81 day average, centered on input time, of F10.7 solar
!                 activity index
!       SFLUX     Daily F10.7 for previous day
!       AP        Geomagnetic activity index array:
!                   (1) Daily Ap
!                   (2) 3 hr ap index for current time
!                   (3) 3 hr ap index for 3 hrs before current time
!                   (4) 3 hr ap index for 6 hrs before current time
!                   (5) 3 hr ap index for 9 hrs before current time
!                   (6) Average of eight 3 hr ap indices from 12 to 33 hrs
!                       prior to current time
!                   (7) Average of eight 3 hr ap indices from 36 to 57 hrs
!                       prior to current time
!                 AP(2:7) are only used when switch_legacy(9) = -1.0 in MSISINIT
!
!     NOTES ON INPUT VARIABLES: 
!       - The day-of-year dependence of the model only uses the DAY argument. If
!         a continuous day-of-year dependence is desired, this argument should
!         include the fractional day (e.g., DAY = <day of year> + UTSEC/86400.0
!       - If lzalt_type = .true. (default) in the MSISINIT call, then Z is
!         treated as geodetic altitude.
!         If lzalt_type = .false., then Z is treated as geopotential height.
!       - F107 and F107A values are the 10.7 cm radio flux at the Sun-Earth
!         distance, not the radio flux at 1 AU. 
!
!     OUTPUT VARIABLES:
!       TN     Temperature at altitude (K)
!       DN(1)  Total mass density (kg/m3)
!       DN(2)  N2 number density (m-3)
!       DN(3)  O2 number density (m-3)
!       DN(4)  O number density (m-3)
!       DN(5)  He number density (m-3)
!       DN(6)  H number density (m-3)
!       DN(7)  Ar number density (m-3)
!       DN(8)  N number density (m-3)
!       DN(9)  Anomalous oxygen number density (m-3)
!       DN(10) Not used in NRLMSIS 2.0 (will contain NO in future release)
!       TEX    Exospheric temperature (K) (optional argument)
!
!     NOTES ON OUTPUT VARIABLES: 
!       - Missing density values are returned as 9.999e-38
!       - Species included in mass density calculation are set in MSISINIT
!
!!! =========================================================================

!**************************************************************************************************
! MSIS_CALC Module: Contains main MSIS entry point
!**************************************************************************************************
module msis_calc

contains

  !==================================================================================================
  ! MSISCALC: The main MSIS subroutine entry point
  !==================================================================================================
  subroutine msiscalc(day,utsec,z,lat,lon,sfluxavg,sflux,ap,tn,dn,tex)

    use msis_constants, only    : rp, dmissing, lnp0, Mbarg0divkB, kB, nspec, nodesTN, nd, zetaF, zetaB, &
                                  Hgamma, zetagamma, maxnbf
    use msis_init, only         : msisinit, initflag, zaltflag, specflag, massflag, masswgt, etaTN
    use msis_gfn, only          : globe
    use msis_tfn, only          : tnparm, tfnparm, tfnx
    use msis_dfn, only          : dnparm, dfnparm, dfnx

    implicit none

    real(8), external          :: alt2gph
    real(kind=rp), external    :: dilog

    real(kind=rp), intent(in)  :: day
    real(kind=rp), intent(in)  :: utsec
    real(kind=rp), intent(in)  :: z
    real(kind=rp), intent(in)  :: lat
    real(kind=rp), intent(in)  :: lon
    real(kind=rp), intent(in)  :: sfluxavg,sflux,ap(1:7)
    real(kind=rp), intent(out) :: tn, dn(1:10)
    real(kind=rp), intent(out), optional :: tex
  
    real(kind=rp), save        :: lastday = -9999.0
    real(kind=rp), save        :: lastutsec = -9999.0
    real(kind=rp), save        :: lastlat = -9999.0
    real(kind=rp), save        :: lastlon = -9999.0
    real(kind=rp), save        :: lastz = -9999.0
    real(kind=rp), save        :: lastsflux = -9999.0
    real(kind=rp), save        :: lastsfluxavg = -9999.0
    real(kind=rp), save        :: lastap(1:7) = -9999.0
    real(kind=rp), save        :: gf(0:maxnbf-1)
    real(kind=rp), save        :: Sz(-5:0,2:6)
    integer, save              :: iz
    type(tnparm), save         :: tpro
    type(dnparm), save         :: dpro(1:nspec-1)

    real(8)                    :: zaltd, latd
    real(kind=rp)              :: zeta, lndtotz, Vz, Wz, HRfact, lnPz, delz
    integer                    :: i, j, kmax, ispec

    ! Check if model has been initialized; if not, perform default initialization
    if (.not. initflag) call msisinit()

    ! Calculate geopotential height, if necessary
    if(zaltflag) then
      zaltd = dble(z)
      latd = dble(lat)
      zeta = alt2gph(latd,zaltd)
    else
      zeta = z
    endif

    ! If only altitude changes then update the local spline weights
    if (zeta .lt. zetaB) then
      if (zeta .ne. lastz) then
        if (zeta .lt. zetaF) then
          kmax = 5
        else
          kmax = 6
        endif
        call bspline(zeta,nodesTN,nd+2,kmax,etaTN,Sz,iz)
        lastz = zeta
      endif
    endif

    ! If location, time, or solar/geomagnetic conditions change then recompute the profile parameters
    if ((day .ne. lastday)     .or. (utsec .ne. lastutsec)       .or. &
        (lat .ne. lastlat)     .or. (lon .ne. lastlon)           .or. &
        (sflux .ne. lastsflux) .or. (sfluxavg .ne. lastsfluxavg) .or. &
        any(ap .ne. lastap)) then
      call globe(day,utsec,lat,lon,sfluxavg,sflux,ap,gf)
      call tfnparm(gf,tpro)
      do ispec = 2, nspec-1
        if (specflag(ispec)) call dfnparm(ispec,gf,tpro,dpro(ispec))
      enddo
      lastday = day
      lastutsec = utsec
      lastlat = lat
      lastlon = lon
      lastsflux = sflux
      lastsfluxavg = sfluxavg
      lastap = ap
    endif

    ! Exospheric temperature
    if (present(tex)) then
      tex = tpro%tex
    endif

    ! Temperature at altitude
    tn = tfnx(zeta,iz,Sz(-3:0,4),tpro)

    ! Temperature integration terms at altitude, total number density
    delz = zeta - zetaB
    if (zeta .lt. zetaF) then
      i = max(iz-4,0)
      if (iz .lt. 4) then
        j = -iz
      else
        j = -4
      endif
      Vz = dot_product(tpro%beta(i:iz),Sz(j:0,5)) + tpro%cVS
      Wz = 0.0_rp
      lnPz = lnP0 - Mbarg0divkB*(Vz - tpro%Vzeta0)
      lndtotz = lnPz - log(kB*tn)
    else
      if (zeta .lt. zetaB) then
        Vz = dot_product(tpro%beta(iz-4:iz),Sz(-4:0,5)) + tpro%cVS
        Wz = dot_product(tpro%gamma(iz-5:iz),Sz(-5:0,6)) + tpro%cVS*delz + tpro%cWS
      else
        Vz = (delz + log(tn/tpro%tex)/tpro%sigma)/tpro%tex + tpro%cVB
        Wz = (0.5_rp*delz*delz + dilog(tpro%b*exp(-tpro%sigma*delz))/tpro%sigmasq)/tpro%tex &
              + tpro%cVB*delz + tpro%cWB
      endif
    endif
        
    ! Species number densities at altitude
    HRfact = 0.5_rp * (1.0_rp + tanh(Hgamma*(zeta - zetagamma)))  !Reduction factor for chemical/dynamical correction scale height below zetagamma
    do ispec = 2, nspec-1
      if (specflag(ispec)) then
        dn(ispec) = dfnx(zeta,tn,lndtotz,Vz,Wz,HRfact,tpro,dpro(ispec))
      else
        dn(ispec) = dmissing
      endif
    enddo

    ! Mass density
    if (specflag(1)) then
      dn(1) = dot_product(dn,masswgt)
    else
      dn(1) = dmissing
    endif

    return

  end subroutine msiscalc

end module msis_calc

!==================================================================================================
! BSPLINE: Returns array of nonzero b-spline values, for all orders up to specified order (max 6)
!==================================================================================================
subroutine bspline(x,nodes,nd,kmax,eta,S,i)

  use msis_constants, only:  rp

  implicit none

  ! Input variables
  real(kind=rp), intent(in)  :: x             !Location at which splines are to be evaluated
  real(kind=rp), intent(in)  :: nodes(0:30)   !Spline node locations
  integer, intent(in)        :: nd            !Number of spline nodes minus one (0:nd)
  integer, intent(in)        :: kmax          !Maximum order (up to 6 allowed) of evaluated splines
  real(kind=rp), intent(in)  :: eta(0:30,2:6) !Array of precomputed weights for recursion (reciprocals of node differences)
  ! Ouput variables
  real(kind=rp), intent(out) :: S(-5:0,2:6)   !Array of b-spline values (spline index relative to i (-5:0), spline order (2:6))
  integer, intent(out)       :: i             !Index of last nonzero b-spline

  ! Working variables
  integer                    :: j, k, l
  integer                    :: low, high
  real(kind=rp)              :: w(-4:0) !Weights for recursion relation
 
  ! Initialize to zero
  S(:,:) = 0.0_rp

  ! Find index of last (rightmost) nonzero spline
  if (x .ge. nodes(nd)) then
    i = nd
    return
  endif
  if (x .le. nodes(0)) then
    i = -1
    return
  endif
  low = 0
  high = nd
  i = (low + high)/2
  do while (x .lt. nodes(i) .or. x .ge. nodes(i + 1))
      if (x .lt. nodes(i)) then
        high = i
      else
        low = i
      endif
      i = (low + high)/2
  end do

  ! Initialize with linear splines
  S(0,2) = (x - nodes(i)) * eta(i,2)
  if (i .gt. 0) S(-1,2) = 1 - S(0,2)
  if (i .ge. nd-1) S(0,2) = 0.0_rp   !Reset out-of-bounds spline to zero

  ! k = 3 (quadratic splines)
  w(:) = 0.0_rp
  w(0) = (x - nodes(i)) * eta(i,3)
  if (i .ne. 0) w(-1) = (x - nodes(i-1)) * eta(i-1,3)
  if (i .lt. (nd-2)) S(0,3) = w(0)*S(0,2)
  if ( ((i-1) .ge. 0) .and. ((i-1) .lt. (nd-2)) ) &
      S(-1,3) = w(-1) * S(-1,2) + (1.0_rp - w(0))*S(0,2)
  if ((i-2) .ge. 0) S(-2,3) = (1.0_rp - w(-1))*S(-1,2)
    
  ! k = 4 (cubic splines)
  do l = 0, -2, -1
    j = i + l
    if (j .lt. 0) exit  !Skip out-of-bounds splines
    w(l) = (x - nodes(j)) * eta(j,4)
  enddo
  if (i .lt. (nd-3)) S(0,4) = w(0)*S(0,3)
  do l = -1, -2, -1
      if ( ((i+l) .ge. 0) .and. ((i+l) .lt. (nd-3)) ) &
          S(l,4) = w(l)*S(l,3) + (1.0_rp - w(l+1))*S(l+1,3)
  enddo
  if ((i-3) .ge. 0) S(-3,4) = (1.0_rp - w(-2))*S(-2,3)
  
  ! k = 5
  do l = 0, -3, -1
    j = i + l
    if (j .lt. 0) exit  !Skip out-of-bounds splines
    w(l) = (x - nodes(j)) * eta(j,5)
  enddo
  if (i .lt. (nd-4)) S(0,5) = w(0)*S(0,4)
  do l = -1, -3, -1
      if ( ((i+l) .ge. 0) .and. ((i+l) .lt. (nd-4)) ) &
          S(l,5) = w(l)*S(l,4) + (1.0_rp - w(l+1))*S(l+1,4)
  enddo
  if ((i-4) .ge. 0) S(-4,5) = (1.0_rp - w(-3))*S(-3,4)
  if (kmax .eq. 5) return  !Exit if only 5th order spline is needed

  ! k = 6
  do l = 0, -4, -1
    j = i + l
    if (j .lt. 0) exit  !Skip out-of-bounds splines
    w(l) = (x - nodes(j)) * eta(j,6)
  enddo
  if (i .lt. (nd-5)) S(0,6) = w(0)*S(0,5)
  do l = -1, -4, -1
    if ( ((i+l) .ge. 0) .and. ((i+l) .lt. (nd-5)) ) &
        S(l,6) = w(l)*S(l,5) + (1.0_rp - w(l+1))*S(l+1,5)
  enddo
  if ((i-5) .ge. 0) S(-5,6) = (1.0_rp - w(-4))*S(-4,5)

  return

end subroutine bspline

!==================================================================================================
! DILOG: Calculate dilogarithm in the domain [0,1)
! Retains terms up to order 3 in the expansion, which results in relative errors less than 1E-5.
! Reference: 
!   Ginsberg, E. S., and D. Zaborowski (1975), The Dilogarithm function of a real argument, 
!   Commun. ACM, 18, 200–202.
!==================================================================================================
real(kind=rp) function dilog(x0)

  use msis_constants, only     : rp, pi

  implicit none

  real(kind=rp), intent(in)   :: x0
  real(kind=rp), parameter    :: pi2_6 = pi*pi / 6.0_rp
  real(kind=rp)               :: x, xx, x4, lnx

  x = x0
  if (x .gt. 0.5_rp) then
    lnx = log(x)
    x = 1.0_rp - x          !Reflect argument into [0,0.5] range
    xx = x*x
    x4 = 4.0_rp*x
    dilog = pi2_6 - lnx*log(x) &
            - (4.0_rp*xx*(23.0_rp/16.0_rp + x/36.0_rp + xx/576.0_rp + xx*x/3600.0_rp) &
                + x4 + 3.0_rp*(1.0_rp - xx)*lnx) / (1.0_rp + x4 + xx)
  else
    xx = x*x
    x4 = 4.0_rp*x
    dilog = (4.0_rp*xx*(23.0_rp/16.0_rp + x/36.0_rp + xx/576.0_rp + xx*x/3600.0_rp) &
              + x4 + 3.0_rp*(1.0_rp - xx)*log(1.0_rp - x)) / (1.0_rp + x4 + xx)
  endif

  return

end function dilog
