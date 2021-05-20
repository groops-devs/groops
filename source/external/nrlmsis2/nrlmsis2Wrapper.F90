!***********************************************/
!**
!* @file nrlmsis2wrapper.F90
!*
!* @brief Fortran Wrapper.
!*
!* @author Sandro Krauss
!* @date   2021-03-18
!*
!/***********************************************/

SUBROUTINE msisinitWrapper(inputPath, inputFile) bind(C, name="msisinitWrapper")

  use msis_init, only                          : msisinit
  use, intrinsic                              :: ISO_C_BINDING

  IMPLICIT NONE
  real(4)                                     :: switches(1:25)=1.0
  character(len=1,kind=C_char), intent(in)    :: inputPath(*)
  character(len=1,kind=C_char), intent(in)    :: inputFile(*)
  character(len=:), allocatable               :: str1, str2
  integer                                     :: i, nchars


  ! DEFINE: Options for initializing NRLMSIS 2.0 model
  ! --------------------------------------------------
  !  switch_legacy   Floating point array (1:25) of legacy switches that
  !  control groups of terms:
  ! 1 - F10.7
  ! 2 - Time independent
  ! 3 - Symmetrical annual
  ! 4 - Symmetrical semiannual
  ! 5 - Asymmetrical annual
  ! 6 - Asymmetrical semiannual
  ! 7 - Diurnal
  ! 8 - Semidiurnal
  ! 9 - Geomagnetic activity:
  !     1.0 = Daily Ap mode
  !    -1.0 = Storm-time ap mode
  ! 10 - All UT/long effects
  ! 11 - Longitudinal
  ! 12 - UT and mixed UT/long
  ! 13 - Mixed Ap/UT/long
  ! 14 - Terdiurnal
  ! 15-25 - Not used in NRLMSIS 2.0
  switches(9) = -1.0

  ! CONVERT: c_char -> char necessary to call the msisinit subroutine
  i = 1
  do
     if (inputPath(i) == c_null_char) exit
     i = i + 1
  end do
  nchars = i - 1  ! Exclude null character from Fortran string
  allocate(character(len=nchars) :: str1)
  str1 = transfer(inputPath(1:nchars), str1)

  i = 1
  do
     if (inputFile(i) == c_null_char) exit
     i = i + 1
  end do
  nchars = i - 1  ! Exclude null character from Fortran string
  allocate(character(len=nchars) :: str2)
  str2 = transfer(inputFile(1:nchars), str2)

  ! initialise the model
  CALL MSISINIT(str1, str2, switch_legacy=switches)

  deallocate(str1)
  deallocate(str2)

END SUBROUTINE msisinitWrapper

!**********************************************

SUBROUTINE msiscalcWrapper(day,utsec,z,lat,lon,sfluxavg,sflux,ap,tn,dn,tex) bind(C, name="msiscalcWrapper")

  use msis_calc, only   : msiscalc
  use, intrinsic       :: ISO_C_BINDING

  IMPLICIT NONE
  real(4), intent(in)  :: day
  real(4), intent(in)  :: utsec
  real(4), intent(in)  :: z
  real(4), intent(in)  :: lat
  real(4), intent(in)  :: lon
  real(4), intent(in)  :: sfluxavg,sflux,ap(1:7)
  real(4), intent(out) :: tn, dn(1:10)
  real(4), intent(out) :: tex

  CALL msiscalc(day,utsec,z,lat,lon,sfluxavg,sflux,ap,tn,dn,tex)

END SUBROUTINE msiscalcWrapper

!**********************************************
