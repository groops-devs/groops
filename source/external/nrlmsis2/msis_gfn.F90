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
! MSIS_GFN Module: Contains subroutines to calculate global (horizontal and time-dependent) model 
!                  basis functions
!**************************************************************************************************
module msis_gfn

  use msis_constants, only : rp, maxn
  use msis_init, only      : TN,PR,N2,O2,O1,HE,H1,AR,N1,OA,NO, swg

  implicit none
  
  real(kind=rp)                :: plg(0:maxn,0:maxn)
  real(kind=rp)                :: cdoy(2), sdoy(2)
  real(kind=rp)                :: clst(3), slst(3)
  real(kind=rp)                :: clon(2), slon(2)
  real(kind=rp)                :: sfluxavgref = 150.0 ! Reference F10.7 value (=150 in NRLMSISE-00)
  real(kind=rp)                :: lastlat = -999.9
  real(kind=rp)                :: lastdoy = -999.9
  real(kind=rp)                :: lastlst = -999.9
  real(kind=rp)                :: lastlon = -999.9
 
contains

  !==================================================================================================
  ! GLOBE: Calculate horizontal and time-dependent basis functions
  !        (Same purpose as NRLMSISE-00 "GLOBE7" subroutine)
  !==================================================================================================
  subroutine globe(doy,utsec,lat,lon,sfluxavg,sflux,ap,bf)

    use msis_constants, only    : deg2rad, doy2rad, lst2rad, &
                                  maxnbf, mbf, maxn, amaxn, amaxs, tmaxl, tmaxn, tmaxs, pmaxm, pmaxn, pmaxs, &
                                  nsfx, nsfxmod, ctimeind, cintann, ctide, cspw, csfx, cextra, cnonlin, csfxmod, cmag, cut
    implicit none

    real(kind=rp), intent(in)  :: doy                       ! Day of year
    real(kind=rp), intent(in)  :: utsec                     ! Universal time in seconds
    real(kind=rp), intent(in)  :: lat                       ! Latitude
    real(kind=rp), intent(in)  :: lon                       ! Longitdue
    real(kind=rp), intent(in)  :: sfluxavg                  ! 81-day average F10.7
    real(kind=rp), intent(in)  :: sflux                     ! Daily F10.7
    real(kind=rp), intent(in)  :: ap(1:7)                   ! Ap geomagnetic activity index history array
    real(kind=rp), intent(out) :: bf(0:maxnbf-1)            ! Output array of basis function terms

    real(kind=rp)              :: lst
    real(kind=rp)              :: slat, clat, clat2, clat4, slat2
    real(kind=rp)              :: cosdoy, sindoy
    real(kind=rp)              :: coslon, sinlon
    real(kind=rp)              :: pl
    real(kind=rp)              :: coslst, sinlst
    real(kind=rp)              :: dfa, df
    real(kind=rp)              :: theta
    real(kind=rp)              :: sza
    integer                    :: n, m, l, s, c

    ! Associated Legendre polynomials
    if (lat .ne. lastlat) then
      clat = sin(lat*deg2rad)  ! clat <=> sin, Legendre polyomial defined in colat
      slat = cos(lat*deg2rad)  ! slat <=> cos, Legendre polyomial defined in colat
      clat2 = clat*clat
      clat4 = clat2*clat2
      slat2 = slat*slat

      plg(0,0) = 1.0_rp
      plg(1,0) = clat
      plg(2,0) = 0.5_rp * (3.0_rp * clat2 - 1.0_rp)
      plg(3,0) = 0.5_rp * (5.0_rp * clat * clat2 - 3.0_rp * clat)
      plg(4,0) = (35.0_rp * clat4 - 30.0_rp * clat2 + 3.0_rp)/8.0_rp
      plg(5,0) = (63.0_rp * clat2 * clat2 * clat - 70.0_rp * clat2 * clat + 15.0_rp * clat)/8.0_rp
      plg(6,0) = (11.0_rp * clat * plg(5, 0) - 5.0_rp * plg(4, 0))/6.0_rp

      plg(1,1) = slat
      plg(2,1) = 3.0_rp * clat * slat
      plg(3,1) = 1.5_rp * (5.0_rp * clat2 - 1.0_rp) * slat
      plg(4,1) = 2.5_rp * (7.0_rp * clat2 * clat - 3.0_rp * clat) * slat
      plg(5,1) = 1.875_rp * (21.0_rp * clat4 - 14.0_rp * clat2 + 1.0_rp) * slat
      plg(6,1) = (11.0_rp * clat * plg(5, 1) - 6.0_rp * plg(4, 1))/5.0_rp

      plg(2,2) = 3.0_rp * slat2
      plg(3,2) = 15.0_rp * slat2 * clat
      plg(4,2) = 7.5_rp * (7.0_rp * clat2 - 1.0_rp) * slat2
      plg(5,2) = 3.0_rp * clat * plg(4, 2) - 2.0_rp * plg(3, 2)
      plg(6,2) = (11.0_rp * clat * plg(5, 2) - 7.0_rp * plg(4, 2))/4.0_rp

      plg(3,3) = 15.0_rp * slat2 * slat
      plg(4,3) = 105.0_rp * slat2 * slat * clat
      plg(5,3) = (9.0_rp * clat * plg(4, 3) - 7.0_rp * plg(3, 3))/2.0_rp
      plg(6,3) = (11.0_rp * clat * plg(5, 3) - 8.0_rp * plg(4, 3))/3.0_rp

      lastlat = lat
    endif

    ! Fourier harmonics of day of year
    if (doy .ne. lastdoy) then
      cdoy(1) = cos(doy2rad*doy)
      sdoy(1) = sin(doy2rad*doy)
      cdoy(2) = cos(doy2rad*doy*2.0_rp)
      sdoy(2) = sin(doy2rad*doy*2.0_rp)
      lastdoy = doy
    endif

    ! Fourier harmonics of local time
    lst = mod(utsec/3600.0_rp + lon/15.0_rp + 24.0_rp, 24.0_rp)
    if (lst .ne. lastlst) then
      clst(1) = cos(lst2rad*lst)
      slst(1) = sin(lst2rad*lst)
      clst(2) = cos(lst2rad*lst*2.0_rp)
      slst(2) = sin(lst2rad*lst*2.0_rp)
      clst(3) = cos(lst2rad*lst*3.0_rp)
      slst(3) = sin(lst2rad*lst*3.0_rp)
      lastlst = lst
    endif

    ! Fourier harmonics of longitude
    if (lon .ne. lastlon) then
      clon(1) = cos(deg2rad*lon)
      slon(1) = sin(deg2rad*lon)
      clon(2) = cos(deg2rad*lon*2.0_rp)
      slon(2) = sin(deg2rad*lon*2.0_rp)
      lastlon = lon
    endif

    !---------------------------------------------
    ! Coupled Linear Terms
    !---------------------------------------------

    ! Reset basis functions
    bf(:) = 0.0_rp

    ! Time-independent (pure latitude dependence)
    c = ctimeind
    do n = 0, amaxn
      bf(c) = plg(n,0)
      c = c + 1
    enddo

    ! Intra-annual (annual and semiannual)
    if (c .ne. cintann) stop 'problem with basis definitions'
    do s = 1, amaxs
      cosdoy = cdoy(s)
      sindoy = sdoy(s)
      do n = 0, amaxn
        pl = plg(n,0)
        bf(c) = pl*cosdoy
        bf(c+1) = pl*sindoy
        c = c + 2
      enddo
    enddo

    ! Migrating Tides (local time dependence)
    if (c .ne. ctide) stop 'problem with basis definitions'
    do l = 1, tmaxl
      coslst = clst(l)
      sinlst = slst(l)
      do n = l, tmaxn
        pl = plg(n,l)
        bf(c) = pl*coslst
        bf(c+1) = pl*sinlst
        c = c + 2
      enddo
      ! Intra-annual modulation of tides
      do s = 1, tmaxs
        cosdoy = cdoy(s)
        sindoy = sdoy(s)
        do n = l, tmaxn
          pl = plg(n,l)
          bf(c) = pl*coslst*cosdoy
          bf(c+1) = pl*sinlst*cosdoy
          bf(c+2) = pl*coslst*sindoy
          bf(c+3) = pl*sinlst*sindoy
          c = c + 4
        enddo
      enddo
    enddo

    ! Stationary Planetary Waves (longitude dependence)
    if (c .ne. cspw) stop 'problem with basis definitions'
    do m = 1, pmaxm
      coslon = clon(m)
      sinlon = slon(m)
      do n = m, pmaxn
        pl = plg(n,m)
        bf(c) = pl*coslon
        bf(c+1) = pl*sinlon
        c = c + 2
      enddo
      ! Intra-annual modulation of SPWs
      do s = 1, pmaxs
        cosdoy = cdoy(s)
        sindoy = sdoy(s)
        do n = m, pmaxn
          pl = plg(n,m)
          bf(c) = pl*coslon*cosdoy
          bf(c+1) = pl*sinlon*cosdoy
          bf(c+2) = pl*coslon*sindoy
          bf(c+3) = pl*sinlon*sindoy
          c = c + 4
        enddo
      enddo
    enddo
    
    ! Linear solar flux terms
    if (c .ne. csfx) stop 'problem with basis definitions'
    dfa = sfluxavg - sfluxavgref
    df = sflux - sfluxavg
    bf(c) = dfa
    bf(c+1) = dfa*dfa
    bf(c+2) = df
    bf(c+3) = df*df
    bf(c+4) = df*dfa
    c = c + nsfx

    ! Additional linear terms
    if (c .ne. cextra) stop 'problem with basis definitions'
    sza = solzen(doy,lst,lat,lon)
    bf(c)    = -0.5_rp*tanh((sza-98.0_rp)/6.0_rp)  !Solar zenith angle logistic function for O, H (transition width 3 deg, transition sza for horizon at ~65 km altitude)
    bf(c+1)  = -0.5_rp*tanh((sza-101.5_rp)/20.0_rp) !Solar zenith angle logistic function for NO (transition width 10 deg, transition sza for horizon at ~130 km altitude)
    bf(c+2)  = dfa*bf(c)                        !Solar flux modulation of logistic sza term
    bf(c+3)  = dfa*bf(c+1)                      !Solar flux modulation of logistic sza term
    bf(c+4)  = dfa*plg(2,0)                     !Solar flux modulation of P(2,0) term
    bf(c+5)  = dfa*plg(4,0)                     !Solar flux modulation of P(4,0) term
    bf(c+6)  = dfa*plg(0,0)*cdoy(1)             !Solar flux modulation of global AO
    bf(c+7)  = dfa*plg(0,0)*sdoy(1)             !Solar flux modulation of global AO
    bf(c+8) = dfa*plg(0,0)*cdoy(2)              !Solar flux modulation of global SAO
    bf(c+9) = dfa*plg(0,0)*sdoy(2)              !Solar flux modulation of global SAO
    
    !---------------------------------------------
    ! Nonlinear Terms
    !---------------------------------------------

    c = cnonlin

    ! Solar flux modulation terms
    if (c .ne. csfxmod) stop 'problem with basis definitions'
    bf(c) = dfa  
    bf(c+1) = dfa*dfa
    bf(c+2) = df 
    bf(c+3) = df*df
    bf(c+4) = df*dfa
    c = c + nsfxmod

    ! Terms needed for legacy geomagnetic activity dependence
    if (c .ne. cmag) stop 'problem with basis set'
    bf(c:c+6) = ap - 4.0
    bf(c+8) =   doy2rad*doy
    bf(c+9) =   lst2rad*lst
    bf(c+10) =  deg2rad*lon
    bf(c+11) =  lst2rad*utsec/3600.0
    bf(c+12) =  abs(lat)
    c = c + 13
    do m = 0,1
      do n = 0,amaxn
        bf(c) = plg(n,m)
        c = c + 1
      enddo
    enddo

    ! Terms needed for legacy UT dependence
    c = cut
    bf(c) =   lst2rad*utsec/3600.0
    bf(c+1) = doy2rad*doy
    bf(c+2) = dfa
    bf(c+3) = deg2rad*lon
    bf(c+4) = plg(1,0)
    bf(c+5) = plg(3,0)
    bf(c+6) = plg(5,0)
    bf(c+7) = plg(3,2)
    bf(c+8) = plg(5,2)

    !---------------------------------------------
    ! Apply Switches
    !---------------------------------------------
    where(.not. swg(0:mbf)) bf(0:mbf) = 0.0_rp
    
    return

  end subroutine globe

  !==================================================================================================
  ! SOLZEN: Calculate solar zenith angle (adapted from IRI subroutine)
  !==================================================================================================
  real(kind=rp) function solzen(ddd,lst,lat,lon)

    use msis_constants, only    : pi, deg2rad

    implicit none

    real(kind=rp), intent(in)  :: ddd
    real(kind=rp), intent(in)  :: lst
    real(kind=rp), intent(in)  :: lat
    real(kind=rp), intent(in)  :: lon

    real(kind=rp)              :: wlon,dec
    real(kind=rp)              :: teqnx,tf,teqt
    real(kind=rp)              :: rlat,phi,cosx
    real(kind=rp), parameter   :: humr = pi/12.0_rp
    real(kind=rp), parameter   :: dumr = pi/182.5_rp
    real(kind=rp), parameter   :: p(5) = (/0.017203534,0.034407068,0.051610602,0.068814136,0.103221204/)

    wlon = 360.0 - lon
    teqnx = ddd + (lst + wlon / 15.0_rp) / 24.0_rp + 0.9369_rp
    teqnx = ddd + 0.9369_rp

    ! Solar declination
    dec = 23.256_rp * sin(p(1) * (teqnx - 82.242_rp)) + 0.381_rp * sin(p(2)*(teqnx - 44.855_rp))  &
         + 0.167_rp * sin(p(3) * (teqnx - 23.355_rp)) - 0.013_rp * sin(p(4)*(teqnx + 11.97_rp)) &
         + 0.011_rp * sin(p(5) * (teqnx - 10.410_rp)) + 0.339137_rp
    dec = dec * deg2rad

    ! Equation of time
    tf = teqnx - 0.5_rp
    teqt = -7.38_rp * sin(p(1) * (tf -  4.0_rp)) - 9.87_rp * sin(p(2) * (tf +  9.0_rp)) &
          + 0.27_rp * sin(p(3) * (tf - 53.0_rp)) -  0.2_rp * cos(p(4) * (tf - 17.0_rp))

    phi = humr * (lst - 12.0_rp) + teqt * deg2rad / 4.0_rp
    rlat = lat * deg2rad

    ! Cosine of solar zenith angle
    cosx = sin(rlat) * sin(dec) + cos(rlat) * cos(dec) * cos(phi)
    if (abs(cosx) .gt. 1.0_rp) cosx = sign(1.0_rp,cosx)

    solzen = acos(cosx) / deg2rad

    return

  end function solzen

  !==================================================================================================
  ! SFLUXMOD: Legacy nonlinear modulation of intra-annual, tide, and SPW terms
  !==================================================================================================
  real(kind=rp) function sfluxmod(iz,gf,parmset,dffact)

    use msis_constants, only       : maxnbf, mbf, csfx, csfxmod
    use msis_init, only            : basissubset, zsfx, tsfx, psfx

    implicit none

    integer, intent(in)           :: iz
    real(kind=rp), intent(in)     :: gf(0:maxnbf-1)
    type(basissubset), intent(in) :: parmset
    real(kind=rp), intent(in)     :: dffact  !Turns on or adjusts the delta-F terms added to F1 and F2 (eqns. A22b and A22c in Hedin (1987)).

    real(kind=rp)                 :: f1, f2, f3, sum
    integer                       :: j

    ! Intra-annual modulation factor
    if (swg(csfxmod)) then
      f1 = parmset%beta(csfxmod,iz) * gf(csfxmod) &
           + (parmset%beta(csfx+2,iz) * gf(csfxmod+2) + parmset%beta(csfx+3,iz) * gf(csfxmod+3) ) * dffact
    else
      f1 = 0.0_rp
    endif

    ! Migrating tide (local time) modulation factor
    if (swg(csfxmod+1)) then
      f2 = parmset%beta(csfxmod+1,iz) * gf(csfxmod) &
           + (parmset%beta(csfx+2,iz) * gf(csfxmod+2) + parmset%beta(csfx+3,iz) * gf(csfxmod+3) ) * dffact
    else
      f2 = 0.0_rp
    endif

    ! SPW (longitude) modulation factor
    if (swg(csfxmod+2)) then
      f3 = parmset%beta(csfxmod+2,iz) * gf(csfxmod)
    else
      f3 = 0.0_rp
    endif

    sum = 0.0
    do j = 0, mbf
      ! Apply intra-annual modulation
      if (zsfx(j)) then
        sum = sum + parmset%beta(j,iz)*gf(j)*f1
        cycle
      endif
      ! Apply migrating tide modulation
      if (tsfx(j)) then
        sum = sum + parmset%beta(j,iz)*gf(j)*f2
        cycle
      endif
      ! Apply SPW modulation
      if (psfx(j)) then
        sum = sum + parmset%beta(j,iz)*gf(j)*f3
        cycle
      endif
    enddo

    sfluxmod = sum

    return

  end function sfluxmod

  !==================================================================================================
  ! GEOMAG: Legacy nonlinear ap dependence (daily ap mode and ap history mode), including mixed 
  !         ap/UT/Longitude terms.
  ! Master switch control is as follows:
  !   swg(cmag) .nor. swg(cmag+1)   Do nothing: Return zero
  !   swg(cmag) .and. swg(cmag+1)   Daily Ap mode
  !   swg(cmag) .neqv. swg(cmag+1)  3-hour ap history mode
  !==================================================================================================
  real(kind=rp) function geomag(p0,bf,plg)

    use msis_constants, only    : nmag, cmag

    implicit none

    real(kind=rp), intent(in)  :: p0(0:nmag-1)
    real(kind=rp), intent(in)  :: bf(0:12)
    real(kind=rp), intent(in)  :: plg(0:6,0:1)

    logical                    :: swg1(0:nmag-1) !Copy of switches
    real(kind=rp)              :: p(0:nmag-1)    !Copy of parameters used to apply switches
    real(kind=rp)              :: delA, gbeta, ex, sumex, G(1:6)
    integer(4)                 :: i

    ! Return zero if both master switches are off    
    if (.not. (swg(cmag) .or. swg(cmag+1))) then
      geomag = 0.0_rp
      return
    endif

    ! Copy parameters
    p = p0
    swg1 = swg(cmag:cmag+nmag-1)

    ! Calculate function
    if (swg1(0) .eqv. swg1(1)) then
      ! Daily Ap mode
      if (p(1) .eq. 0) then  !If k00s is zero, then cannot compute function
        geomag = 0.0_rp
        return
      endif
      where(.not. swg1(2:25)) p(2:25) = 0.0_rp !Apply switches
      p(8) = p0(8) !Need doy phase term
      delA = G0fn(bf(0),p(0),p(1))
      geomag = ( p(2)*plg(0,0) +  p(3)*plg(2,0) +  p(4)*plg(4,0)                     &  ! time independent
        + (p(5)*plg(1,0) + p(6)*plg(3,0) + p(7)*plg(5,0)) * cos(bf(8) - p(8))        &  ! doy modulation
        + (p(9)*plg(1,1) + p(10)*plg(3,1) + p(11)*plg(5,1)) * cos(bf(9) - p(12))     &  ! local time modulation
        + (1.0_rp + p(13)*plg(1,0)) *                                                &
          (p(14)*plg(2,1) + p(15)*plg(4,1) + p(16)*plg(6,1)) * cos(bf(10) - p(17))   &  ! longitude effect
        + (p(18)*plg(1,1) + p(19)*plg(3,1) + p(20)*plg(5,1)) * cos(bf(10) - p(21)) * &
          cos(bf(8) - p(8))                                                          &  ! longitude with doy modulaiton
        + (p(22)*plg(1,0) + p(23)*plg(3,0) + p(24)*plg(5,0)) * cos(bf(11) - p(25)) ) &  ! universal time
        *delA
    else
      ! 3-hour ap history mode
      if (p(28) .eq. 0) then  !If beta00 is zero, then cannot compute function
        geomag = 0.0
        return
      endif
      where(.not. swg1(30:)) p(30:) = 0.0  !Apply switches
      p(36) = p0(36) !Need doy phase term
      gbeta = p(28)/(1 + p(29)*(45.0_rp - bf(12)))
      ex = exp(-10800.0_rp*gbeta)
      sumex = 1 + (1 - ex**19.0_rp) * ex**(0.5_rp) / (1 - ex)
      do i = 1, 6
        G(i) = G0fn(bf(i),p(26),p(27))
      enddo
      delA = ( G(1)                                                                 &
                    + ( G(2)*ex + G(3)*ex*ex + G(4)*ex**3.0_rp                      &
                       +(G(5)*ex**4.0_rp + G(6)*ex**12.0_rp)*(1-ex**8.0_rp)/(1-ex) ) ) / sumex
      geomag = ( p(30)*plg(0,0) +  p(31)*plg(2,0) +  p(32)*plg(4,0)                  &  ! time independent
        + (p(33)*plg(1,0) + p(34)*plg(3,0) + p(35)*plg(5,0)) * cos(bf(8) - p(36))    &  ! doy modulation
        + (p(37)*plg(1,1) + p(38)*plg(3,1) + p(39)*plg(5,1)) * cos(bf(9) - p(40))    &  ! local time modulation
        + (1.0_rp + p(41)*plg(1,0)) *                                                &
          (p(42)*plg(2,1) + p(43)*plg(4,1) + p(44)*plg(6,1)) * cos(bf(10) - p(45))   &  ! longitude effect
        + (p(46)*plg(1,1) + p(47)*plg(3,1) + p(48)*plg(5,1)) * cos(bf(10) - p(49)) * &
          cos(bf(8) - p(36))                                                         &  ! longitude with doy modulaiton
        + (p(50)*plg(1,0) + p(51)*plg(3,0) + p(52)*plg(5,0)) * cos(bf(11) - p(53)) ) &  ! universal time
        *delA
    endif

    return

    contains

      real(kind=rp) function G0fn(a,k00r,k00s)
          real(kind=rp),intent(in)  :: a, k00r, k00s
          G0fn = a + (k00r - 1.0_rp) * (a + (exp(-a*k00s) - 1.0_rp)/k00s)
          return
      end function G0fn

  end function geomag

  !==================================================================================================
  ! UTDEP: Legacy nonlinear UT dependence
  !==================================================================================================
  real(kind=rp) function utdep(p0,bf)

    use msis_constants, only    : nut, cut

    implicit none

    real(kind=rp), intent(in)  :: p0(0:nut-1)
    real(kind=rp), intent(in)  :: bf(0:8)

    real(kind=rp)              :: p(0:nut-1)    !Copy of parameters used to apply switches
    logical                    :: swg1(0:nut-1) !Copy of switches

    !Copy parameters
    p = p0
    swg1 = swg(cut:cut+nut-1)
    where(.not. swg1(3:nut-1)) p(3:nut-1) = 0.0  !Apply switches

    ! Calculate function
    utdep = cos(bf(0)-p(0)) *                          &
            (1 + p(3)*bf(4)*cos(bf(1)-p(1)))  *        &
            (1 + p(4)*bf(2)) * (1 + p(5)*bf(4)) *      &
            (p(6)*bf(4) + p(7)*bf(5) + p(8)*bf(6)) +   &
            cos(bf(0)-p(2)+2*bf(3)) * (p(9)*bf(7) + p(10)*bf(8)) * (1 + p(11)*bf(2))

    return

  end function utdep

end module msis_gfn
