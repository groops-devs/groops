! Minimal kind definitions used by the NRLMSISE-00 Fortran implementation.
module groops_msise00_kinds
  use, intrinsic :: iso_c_binding
  implicit none
  integer, parameter :: RL = c_double
  integer, parameter :: IT = c_int
end module groops_msise00_kinds
