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
! MSIS_DFN Module: Contains vertical species density profile parameters and subroutines
!**************************************************************************************************
module msis_dfn

  use msis_constants, only : rp, nl, nsplO1, nsplNO

  type dnparm
    sequence
    real(kind=rp)         :: lnPhiF            ! (Except O, H) Natural log of mixing ratio at zetaF (70 km), before chemical and dynamical corrections are applied (ln m^-3) (global term only)
    real(kind=rp)         :: lndref            ! Natural log of number density at reference height
    real(kind=rp)         :: zetaM             ! "Turbopause Height": Height of midpoint of effective mass transition (km)
    real(kind=rp)         :: HML               ! Scale height of lower portion of effective mass profile (km)
    real(kind=rp)         :: HMU               ! Scale height of upper portion of effective mass profile (km)
    real(kind=rp)         :: C                 ! Chapman term coefficient
    real(kind=rp)         :: zetaC             ! Chapman term reference height (km)
    real(kind=rp)         :: HC                ! Chapman term scale height (km)
    real(kind=rp)         :: R                 ! Chemical/dynamical term coefficient
    real(kind=rp)         :: zetaR             ! Chemical/dynamical term reference height (km)
    real(kind=rp)         :: HR                ! Chemical/dynamical term scale height (km)
    real(kind=rp)         :: cf(0:nsplO1+1)    ! Merged spline coefficients (for chemistry-dominated region of O1, NO, and (eventually), H, N)
    real(kind=rp)         :: zref              ! Reference height for hydrostatic integral and ideal gas terms
    real(kind=rp)         :: Mi(0:4)           ! Effective mass at nodes of piecewise mass profile (derived from zetaM, HML, HMU)
    real(kind=rp)         :: zetaMi(0:4)       ! Height of nodes of piecewise mass profile (derived from zetaM, HML, HMU)
    real(kind=rp)         :: aMi(0:4) = 0.0_rp ! Slopes of piecewise mass profile segments (derived from zetaM, HML, HMU)
    real(kind=rp)         :: WMi(0:4) = 0.0_rp ! 2nd indefinite integral of 1/T at mass profile nodes
    real(kind=rp)         :: XMi(0:4) = 0.0_rp ! Cumulative adjustment to M/T integral due to changing effective mass
    real(kind=rp)         :: Izref             ! Indefinite hydrostatic integral at reference height
    real(kind=rp)         :: Tref              ! Temperature at reference height (for ideal gas law term)
    real(kind=rp)         :: zmin              ! Minimum height of profile (missing values below)
    real(kind=rp)         :: zhyd              ! Hydrostatic terms needed above this height
    integer(kind=rp)      :: ispec             ! Species index
  end type dnparm

  contains

  !==================================================================================================
  ! DFNPARM: Compute the species density profile parameters
  !==================================================================================================
  subroutine dfnparm(ispec,gf,tpro,dpro)

    use msis_constants, only   : tanh1, specmass, lnvmr, Mbar, g0divkB, &
                                 nd, zetaF, zetaB, zetaA, nodesTN, &
                                 nodesO1, zetarefO1, c1o1, c1o1adj, &
                                 nodesNO, zetarefNO, c1NO, c1NOadj, &
                                 zetarefOA, &
                                 maxnbf, mbf, nmag, nut, cmag, cut
    use msis_init, only        : etaTN, TN,PR,N2,O2,O1,HE,H1,AR,N1,OA,NO, N2Rflag, &
                                 HRfactO1ref, dHRfactO1ref, HRfactNOref, dHRfactNOref
    use msis_gfn, only         : sfluxmod, geomag, utdep
    use msis_tfn, only         : tnparm

    implicit none

    real(kind=rp), external  :: dilog

    integer, intent(in)       :: ispec          ! Species index
    real(kind=rp), intent(in) :: gf(0:maxnbf-1) ! Array of horizontal and temporal basis function terms   
    type(tnparm), intent(in)  :: tpro           ! Structure containing temperature vertical profile parameters
    type(dnparm), intent(out) :: dpro           ! Output structure containing density vertical profile parameters
        
    integer                   :: izf, i, i1, iz
    real(kind=rp)             :: Cterm, Rterm0, Rterm
    real(kind=rp)             :: bc(2)
    real(kind=rp)             :: hbetaL,hbetaU
    real(kind=rp)             :: delM, delz
    real(kind=rp)             :: Wi           ! 2nd indefinite integral at a piecewise mass profile node
    real(kind=rp)             :: Si(-5:0,2:6) ! Array of b-spline values at a mass profile node
    real(kind=rp)             :: Mzref        ! Effective mass at reference altitude

    dpro%ispec = ispec

    select case(ispec)

    ! Molecular Nitrogen ----------------------
    case(2)
      ! Mixing ratio and reference number density
      dpro%lnPhiF = lnvmr(ispec)
      dpro%lndref = tpro%lndtotF + dpro%lnPhiF
      dpro%zref = zetaF
      dpro%zmin = -1.0_rp
      dpro%zhyd = zetaF
      ! Effective mass
      dpro%zetaM = dot_product(N2%beta(0:mbf,1),gf(0:mbf))
      dpro%HML   = N2%beta(0,2)
      dpro%HMU   = N2%beta(0,3)
      ! Photochemical correction
      dpro%R     = 0.0_rp
      if (N2Rflag) dpro%R = dot_product(N2%beta(0:mbf,7),gf(0:mbf))
      dpro%zetaR = N2%beta(0,8)
      dpro%HR    = N2%beta(0,9)

    ! Molecular Oxygen ------------------------
    case(3)
      ! Mixing ratio and reference number density
      dpro%lnPhiF = lnvmr(ispec)
      dpro%lndref = tpro%lndtotF + dpro%lnPhiF
      dpro%zref = zetaF
      dpro%zmin = -1.0_rp
      dpro%zhyd = zetaF
      ! Effective mass
      dpro%zetaM = O2%beta(0,1)
      dpro%HML   = O2%beta(0,2)
      dpro%HMU   = O2%beta(0,3)
      ! Photochemical correction
      dpro%R     = dot_product(O2%beta(0:mbf,7),gf(0:mbf))
      dpro%R     = dpro%R + geomag(O2%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%zetaR = O2%beta(0,8)
      dpro%HR    = O2%beta(0,9)

    ! Atomic Oxygen --------------------------
    case(4)
      ! Reference number density
      dpro%lnPhiF = 0.0_rp
      dpro%lndref = dot_product(O1%beta(0:mbf,0),gf(0:mbf))
      dpro%zref = zetarefO1
      dpro%zmin = nodesO1(3)
      dpro%zhyd = zetarefO1
      ! Effective mass
      dpro%zetaM = O1%beta(0,1)
      dpro%HML   = O1%beta(0,2)
      dpro%HMU   = O1%beta(0,3)
      ! Chapman correction
      dpro%C     = dot_product(O1%beta(0:mbf,4),gf(0:mbf))
      dpro%zetaC = O1%beta(0,5)
      dpro%HC    = O1%beta(0,6)
      ! Dynamical correction
      dpro%R     = dot_product(O1%beta(0:mbf,7),gf(0:mbf))
      dpro%R     = dpro%R + sfluxmod(7,gf,O1,0.0_rp)         
      dpro%R     = dpro%R + geomag(O1%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%R     = dpro%R + utdep(O1%beta(cut:cut+nut-1,7),gf(cut:cut+8))
      dpro%zetaR = O1%beta(0,8)
      dpro%HR    = O1%beta(0,9)
      ! Unconstrained splines
      do izf = 0, nsplO1-1
        dpro%cf(izf) = dot_product(O1%beta(0:mbf,izf+10),gf(0:mbf))
      enddo
      ! Constrained splines calculated after case statement

    ! Helium ----------------------
    case(5)
      ! Mixing ratio and reference number density
      dpro%lnPhiF = lnvmr(ispec)
      dpro%lndref = tpro%lndtotF + dpro%lnPhiF
      dpro%zref = zetaF
      dpro%zmin = -1.0_rp
      dpro%zhyd = zetaF
      ! Effective mass
      dpro%zetaM = HE%beta(0,1)
      dpro%HML   = HE%beta(0,2)
      dpro%HMU   = HE%beta(0,3)
      ! Dynamical correction
      dpro%R     = dot_product(HE%beta(0:mbf,7),gf(0:mbf))
      dpro%R     = dpro%R + sfluxmod(7,gf,HE,1.0_rp)         
      dpro%R     = dpro%R + geomag(HE%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%R     = dpro%R + utdep(HE%beta(cut:cut+nut-1,7),gf(cut:cut+8))
      dpro%zetaR = HE%beta(0,8)
      dpro%HR    = HE%beta(0,9)

    ! Atomic Hydrogen ----------------------
    case(6)
      ! Reference number density
      dpro%lnPhiF = 0.0_rp
      dpro%lndref = dot_product(H1%beta(0:mbf,0),gf(0:mbf))
      dpro%zref = zetaA
      dpro%zmin = 75.0_rp
      dpro%zhyd = zetaF
      ! Effective mass
      dpro%zetaM = H1%beta(0,1)
      dpro%HML   = H1%beta(0,2)
      dpro%HMU   = H1%beta(0,3)
      ! Chapman correction
      dpro%C     = dot_product(H1%beta(0:mbf,4),gf(0:mbf))
      dpro%zetaC = dot_product(H1%beta(0:mbf,5),gf(0:mbf))
      dpro%HC    = H1%beta(0,6)
      ! Dynamical correction
      dpro%R     = dot_product(H1%beta(0:mbf,7),gf(0:mbf))
      dpro%R     = dpro%R + sfluxmod(7,gf,H1,0.0_rp)        
      dpro%R     = dpro%R + geomag(H1%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%R     = dpro%R + utdep(H1%beta(cut:cut+nut-1,7),gf(cut:cut+8))
      dpro%zetaR = H1%beta(0,8)
      dpro%HR    = H1%beta(0,9)

    ! Argon ----------------------
    case(7)
      ! Mixing ratio and reference number density
      dpro%lnPhiF = lnvmr(ispec)
      dpro%lndref = tpro%lndtotF + dpro%lnPhiF
      dpro%zref = zetaF
      dpro%zmin = -1.0_rp
      dpro%zhyd = zetaF
      ! Effective mass
      dpro%zetaM = AR%beta(0,1)
      dpro%HML   = AR%beta(0,2)
      dpro%HMU   = AR%beta(0,3)
      ! Dynamical correction
      dpro%R     = dot_product(AR%beta(0:mbf,7),gf(0:mbf))
      dpro%R     = dpro%R + geomag(AR%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%R     = dpro%R + utdep(AR%beta(cut:cut+nut-1,7),gf(cut:cut+8))
      dpro%zetaR = AR%beta(0,8)
      dpro%HR    = AR%beta(0,9)

    ! Atomic Nitrogen ----------------------
    case(8)
      ! Reference number density
      dpro%lnPhiF = 0.0_rp
      dpro%lndref = dot_product(N1%beta(0:mbf,0),gf(0:mbf))
      dpro%lndref = dpro%lndref + sfluxmod(0,gf,N1,0.0_rp)         
      dpro%lndref = dpro%lndref + geomag(N1%beta(cmag:cmag+nmag-1,0),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%lndref = dpro%lndref + utdep(N1%beta(cut:cut+nut-1,0),gf(cut:cut+8))
      dpro%zref = zetaB
      dpro%zmin = 90.0_rp
      dpro%zhyd = zetaF
      ! Effective mass
      dpro%zetaM = N1%beta(0,1)
      dpro%HML   = N1%beta(0,2)
      dpro%HMU   = N1%beta(0,3)
      ! Chapman correction
      dpro%C     = N1%beta(0,4)
      dpro%zetaC = N1%beta(0,5)
      dpro%HC    = N1%beta(0,6)
      ! Dynamical correction
      dpro%R     = dot_product(N1%beta(0:mbf,7),gf(0:mbf))
      dpro%zetaR = N1%beta(0,8)
      dpro%HR    = N1%beta(0,9)

    ! Anomalous Oxygen ----------------------
    case(9)
      dpro%lndref = dot_product(OA%beta(0:mbf,0),gf(0:mbf))
      dpro%lndref = dpro%lndref + geomag(OA%beta(cmag:cmag+nmag-1,0),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%zref = zetarefOA
      dpro%zmin = 120.0_rp
      dpro%zhyd = 0.0_rp
      dpro%C     = OA%beta(0,4)
      dpro%zetaC = OA%beta(0,5)
      dpro%HC    = OA%beta(0,6)
      return !No further parameters needed for legacy anomalous oxygen profile

    ! Nitic Oxide ----------------------
    case(10)
      ! Skip if parameters are not defined
      if (NO%beta(0,0) .eq. 0.0_rp) then
          dpro%lndref = 0.0_rp
          return
      endif
      ! Reference number density
      dpro%lnPhiF = 0.0_rp
      dpro%lndref = dot_product(NO%beta(0:mbf,0),gf(0:mbf))
      dpro%zref = zetarefNO
      dpro%zmin = nodesNO(3)
      dpro%zhyd = zetarefNO
      ! Effective mass
      dpro%zetaM = dot_product(NO%beta(0:mbf,1),gf(0:mbf))
      dpro%HML   = dot_product(NO%beta(0:mbf,2),gf(0:mbf))
      dpro%HMU   = dot_product(NO%beta(0:mbf,3),gf(0:mbf))
      ! Chapman correction
      dpro%C     = dot_product(NO%beta(0:mbf,4),gf(0:mbf))
      dpro%zetaC = dot_product(NO%beta(0:mbf,5),gf(0:mbf))
      dpro%HC    = dot_product(NO%beta(0:mbf,6),gf(0:mbf))
      ! Dynamical correction
      dpro%R     = dot_product(NO%beta(0:mbf,7),gf(0:mbf))
      dpro%zetaR = dot_product(NO%beta(0:mbf,8),gf(0:mbf))
      dpro%HR    = dot_product(NO%beta(0:mbf,9),gf(0:mbf))
      ! Unconstrained splines
      do izf = 0,nsplNO-1
          dpro%cf(izf) = dot_product(NO%beta(0:mbf,izf+10),gf(0:mbf))
      enddo
      ! Constrained splines calculated after case statement

! Failsafe -----   ---------------------------
    case default
      stop 'Species not yet implemented'

    endselect
        
    ! Compute piecewise mass profile values and integration terms
    dpro%zetaMi(0) = dpro%zetaM - 2.0_rp*dpro%HML
    dpro%zetaMi(1) = dpro%zetaM - dpro%HML
    dpro%zetaMi(2) = dpro%zetaM
    dpro%zetaMi(3) = dpro%zetaM + dpro%HMU
    dpro%zetaMi(4) = dpro%zetaM + 2.0_rp*dpro%HMU
    dpro%Mi(0) = Mbar
    dpro%Mi(4) = specmass(ispec)
    dpro%Mi(2) = (dpro%Mi(0) + dpro%Mi(4)) / 2.0_rp
    delM = tanh1 * (dpro%Mi(4) - dpro%Mi(0)) / 2.0_rp
    dpro%Mi(1) = dpro%Mi(2) - delM
    dpro%Mi(3) = dpro%Mi(2) + delM
    !do i = 0, 4
    !  i1 = i + 1
    !  if (i .lt. 4) dpro%aMi(i) = (dpro%Mi(i1) - dpro%Mi(i)) / (dpro%zetaMi(i1) - dpro%zetaMi(i))
    !  delz = dpro%zetaMi(i) - zetaB
    !  if (dpro%zetaMi(i) .lt. zetaB) then
    !    call bspline(dpro%zetaMi(i),nodesTN,nd+2,6,etaTN,Si,iz)
    !    dpro%WMi(i) = dot_product(tpro%gamma(iz-5:iz),Si(:,6)) + tpro%cVS*delz + tpro%cWS
    !  else
    !    dpro%WMi(i) = (0.5_rp*delz*delz + dilog(tpro%b*exp(-tpro%sigma*delz))/tpro%sigmasq)/tpro%tex &
    !                  + tpro%cVB*delz + tpro%cWB
    !  endif
    !end do
    do i = 0, 3
      dpro%aMi(i) = (dpro%Mi(i+1) - dpro%Mi(i)) / (dpro%zetaMi(i+1) - dpro%zetaMi(i))
    enddo
    do i = 0, 4
      delz = dpro%zetaMi(i) - zetaB
      if (dpro%zetaMi(i) .lt. zetaB) then
        call bspline(dpro%zetaMi(i),nodesTN,nd+2,6,etaTN,Si,iz)
        dpro%WMi(i) = dot_product(tpro%gamma(iz-5:iz),Si(:,6)) + tpro%cVS*delz + tpro%cWS
      else
        dpro%WMi(i) = (0.5_rp*delz*delz + dilog(tpro%b*exp(-tpro%sigma*delz))/tpro%sigmasq)/tpro%tex &
                      + tpro%cVB*delz + tpro%cWB
      endif
    end do
    dpro%XMi(0) = -dpro%aMi(0) * dpro%WMi(0)
    do i = 1, 3
      dpro%XMi(i) = dpro%XMi(i-1) - dpro%WMi(i) * (dpro%aMi(i) - dpro%aMi(i-1))
    end do
    dpro%XMi(4) = dpro%XMi(3) + dpro%WMi(4) * dpro%aMi(3)

    ! Calculate hydrostatic integral at reference height, and copy temperature
    if (dpro%zref .eq. zetaF) then
      Mzref = Mbar
      dpro%Tref = tpro%TzetaF
      dpro%Izref = Mbar * tpro%VzetaF
    else if (dpro%zref .eq. zetaB) then
      Mzref = pwmp(dpro%zref,dpro%zetaMi,dpro%Mi,dpro%aMi)
      dpro%Tref = tpro%Tb0
      dpro%Izref = 0.0_rp
      if ((zetaB .gt. dpro%zetaMi(0)) .and. (zetaB .lt. dpro%zetaMi(4))) then
        i = 0
        do i1 = 1, 3
          if (zetaB .lt. dpro%zetaMi(i1)) then
            exit
          else
            i = i1
          endif
        enddo
        dpro%Izref = dpro%Izref -  dpro%XMi(i)
      else
        dpro%Izref = dpro%Izref - dpro%XMi(4)                
      endif
    else if (dpro%zref .eq. zetaA) then
      Mzref = pwmp(dpro%zref,dpro%zetaMi,dpro%Mi,dpro%aMi)
      dpro%Tref = tpro%TzetaA
      dpro%Izref = Mzref * tpro%VzetaA
      if ((zetaA .gt. dpro%zetaMi(0)) .and. (zetaA .lt. dpro%zetaMi(4))) then
        i = 0
        do i1 = 1, 3
          if (zetaA .lt. dpro%zetaMi(i1)) then
            exit
          else
            i = i1
          endif
        enddo
        dpro%Izref = dpro%Izref - (dpro%aMi(i)*tpro%WzetaA + dpro%XMi(i))
      else
        dpro%Izref = dpro%Izref - dpro%XMi(4)                
      endif
    else 
      stop 'Integrals at reference height not available'
    endif

    ! C1 constraint for O1 at 85 km
    if (ispec .eq. 4) then
      Cterm = dpro%C*exp(-(dpro%zref-dpro%zetaC)/dpro%HC)
      Rterm0 = tanh((dpro%zref-dpro%zetaR)/(HRfactO1ref*dpro%HR))
      Rterm = dpro%R*(1+Rterm0)
      bc(1) = dpro%lndref - Cterm + Rterm - dpro%cf(7)*c1o1adj(1)      !Reference density, Chapman term, logistic term, and subtraction of last unconstrained spline(7)
      bc(2) = -Mzref*g0divkB/tpro%tzetaA &   !Gradient of hydrostatic term
              -tpro%dlntdzA &  !Gradient of ideal gas law term
              +Cterm/dpro%HC & !Gradient of Chapman term
              +Rterm*(1-Rterm0)/dpro%HR*dHrfactO1ref  & !Gradient of tapered logistic term
              -dpro%cf(7)*c1o1adj(2)  !Subtraction of gradient of last unconstrained spline(7)
      ! Compute coefficients for constrained splines
      dpro%cf(8:9) = matmul(bc,c1o1)
    endif
        
    ! C1 constraint for NO at 122.5 km
    if (ispec .eq. 10) then
      Cterm = dpro%C*exp(-(dpro%zref - dpro%zetaC)/dpro%HC)
      Rterm0 = tanh((dpro%zref-dpro%zetaR)/(HRfactNOref*dpro%HR))
      Rterm = dpro%R*(1+Rterm0)
      bc(1) = dpro%lndref - Cterm + Rterm - dpro%cf(7)*c1noadj(1)      !Reference density, Chapman term, logistic term, and subtraction of last unconstrained spline(7)
      bc(2) = -Mzref*g0divkB/tpro%tb0 &   !Gradient of hydrostatic term
              -tpro%tgb0/tpro%tb0 &  !Gradient of ideal gas law term
              +Cterm/dpro%HC & !Gradient of Chapman term
              +Rterm*(1-Rterm0)/dpro%HR*dHrfactNOref  & !Gradient of tapered logistic term
              -dpro%cf(7)*c1noadj(2)  !Subtraction of gradient of last unconstrained spline(7)
      ! Compute coefficients for constrained splines
      dpro%cf(8:9) = matmul(bc,c1no)
    endif
        
  return

  end subroutine dfnparm

  !==================================================================================================
  ! DFNX: Compute a species density at specified geopotential height
  !==================================================================================================
  real(kind=rp) function dfnx(z,tnz,lndtotz,Vz,Wz,HRfact,tpro,dpro)

    use msis_constants, only   : dmissing, g0divkB, ndO1, nodesO1, ndNO, nodesNO, HOA
    use msis_init, only        : etaO1, etaNO
    use msis_tfn, only         : tnparm

    implicit none

    real(kind=rp), intent(in) :: z            ! Geopotential height
    real(kind=rp), intent(in) :: tnz, lndtotz ! Temperature, total number density at input z 
    real(kind=rp), intent(in) :: Vz, Wz       ! First and second indefinite integrals of 1/T at z
    real(kind=rp), intent(in) :: HRfact       ! Reduction factor for chemical/dynamical correction scale height below zetaF
    type(tnparm), intent(in)  :: tpro         ! Structure containing temperature vertical profile parameters
    type(dnparm), intent(in)  :: dpro         ! Structure containing density vertical profile parameters

    integer(4)                :: i, i1, iz
    real(kind=rp)             :: Mz
    real(kind=rp)             :: Sz(-5:0,2:6)
    real(kind=rp)             :: Ihyd         ! Hydrostatic definite integral
    real(kind=rp)             :: ccor         ! Chapman and logistical corrections

    ! Below minimum height of profile
    if (z .lt. dpro%zmin) then
      dfnx = dmissing
      return
    endif
        
    ! Anomalous Oxygen (legacy MSISE-00 formulation)
    if (dpro%ispec .eq. 9) then
      dfnx = dpro%lndref - (z - dpro%zref)/HOA - dpro%C*exp(-(z-dpro%zetaC)/dpro%HC)
      dfnx = exp(dfnx)
      return               !No further calculation needed for anomalous oxygen
    endif

    ! Nitric Oxide: Skip if parameters are not defined
    if (dpro%ispec .eq. 10) then
      if (dpro%lndref .eq. 0.0_rp) then
        dfnx = dmissing
        return
      endif
    endif
        
    ! Chapman and logistic corrections
    select case(dpro%ispec)
    case(2,3,5,7)             !For N2, O2, He, and Ar: logistic correction only
      ccor =   dpro%R*(1+tanh((z-dpro%zetaR)/(HRfact*dpro%HR)))
    case(4,6,8,10)            !For O, H, N, and NO: Chapman and logistic corrections
      ccor = - dpro%C*exp(-(z-dpro%zetaC)/dpro%HC) &
             + dpro%R*(1+tanh((z-dpro%zetaR)/(HRfact*dpro%HR)))
    endselect

    ! Below height where hydrostatic terms are needed
    if (z .lt. dpro%zhyd) then
      select case(dpro%ispec)
      case(2,3,5,7)           !For N2, O2, He, and Ar, apply mixing ratios and exit
        dfnx = exp(lndtotz + dpro%lnPhiF + ccor)
        return
      case(4)                 !For O, evaluate splines
        call bspline(z,nodesO1,ndO1,4,etaO1,Sz,iz)
        dfnx = exp(dot_product(dpro%cf(iz-3:iz),Sz(-3:0,4)))
        return
      case(10)                !For NO, evaluate splines
        call bspline(z,nodesNO,ndNO,4,etaNO,Sz,iz)
        dfnx = exp(dot_product(dpro%cf(iz-3:iz),Sz(-3:0,4)))
        return
      endselect
    endif
        
    ! Calculate hydrostatic term and apply to reference density
    Mz = pwmp(z,dpro%zetaMi,dpro%Mi,dpro%aMi)
    Ihyd = Mz * Vz - dpro%Izref
    if ((z .gt. dpro%zetaMi(0)) .and. (z .lt. dpro%zetaMi(4))) then
      i = 0
      do i1 = 1, 3
        if (z .lt. dpro%zetaMi(i1)) then
          exit
        else
          i = i1
        endif
      enddo
      Ihyd = Ihyd - (dpro%aMi(i)*Wz + dpro%XMi(i))
    else if (z .ge. dpro%zetaMi(4)) then
      Ihyd = Ihyd - dpro%XMi(4)                
    endif
    dfnx = dpro%lndref - Ihyd * g0divkB  + ccor
        
    ! Apply ideal gas law                
    dfnx = exp(dfnx) * dpro%Tref/tnz

    return

  end function dfnx

  !==================================================================================================
  ! PWMP: Piecewise effective mass profile interpolation
  !==================================================================================================
  real(kind=rp) function pwmp(z,zm,m,dmdz)

    use msis_constants, only   : rp

    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: zm(0:4)
    real(kind=rp), intent(in) :: m(0:4)
    real(kind=rp), intent(in) :: dmdz(0:3)

    integer                   :: irng  !Index of piecwise interval
    integer                   :: inode

    ! Most probable case
    if (z .ge. zm(4)) then
      pwmp = m(4)
      return
    endif

    ! Second most probable case
    if (z .le. zm(0)) then
      pwmp = m(0)
      return
    endif

    ! None of the above
    do inode = 0,3
      if (z .lt. zm(inode+1)) then
        pwmp = m(inode) + dmdz(inode)*(z - zm(inode))
        return
      endif
    enddo

    ! If we are here this is a problem
    stop 'Error in pwmp'

  end function pwmp

end module msis_dfn

