Horizontal Wind Model 2014 (HWM14)
Version HWM14.123114

VERSION HISTORY
  31 Dec 2014 (Release version)

PACKAGE CONTENTS

  README.txt           This file
  hwm14.f90            Source code for HWM14 and DWM07
  hwm123114.bin        Binary Data file containing quiet-time HWM parameters
  dwm07b104i.dat       Data file containing DWM parameters
  gd2qd.dat            Data file containing parameters for QD coordinate
                       conversion
  makefile             Very simple example Makefile
  checkhwm14.f90       Test driver and test output data
  Check/ifort.opt.txt  Test driver results, Intel Fortran (v14.0.0) default flags
  Check/ifort.nopt.txt Test driver results, Intel Fortran (v14.0.0) no optimization
  Check/gfortan.txt    Test driver results, gfortran (v4.8.1) default flags

POINT OF CONTACT
  Douglas P. Drob and John Emmert
  Space Science Division
  Naval Research Laboratory
  4555 Overlook Ave.
  Washington, DC 20375
  douglas.drob@nrl.navy.mil


REFERENCES
  Drob, D. P., et al. (2015), An Update to the Horizontal Wind Model (HWM): The
    Quiet Time Thermosphere, Earth and Space Science, submitted.
  Drob, D. P., et al. (2008), An Empirical Model of the Earth's Horizontal
    Wind Fields: HWM07, J. Geophys Res., 113, doi:10.1029/2008JA013668.
  Emmert, J. T., et al. (2008), DWM07 global empirical model of upper
    thermospheric storm-induced disturbance winds, J. Geophys Res., 113,
    doi:10.1029/2008JA013541.

COMPILER NOTES

  - The model package has been compiled succesfully with the following
    FORTRAN compiler system combinations: Intel (Windows, Linux, Mac),
    pgi (Linux), and gfortran (Mac, Linux).

  - No special compiler flags or options should be needed; but as the result
    of rounding, optimization, and floating point instruction reordering
    various compiler result in slightly different numerical results on the order
    of 1 x 10^-4 (i.e. +/- 1 mm/s) in certain regions of space and time. These
    differences are negligable with respect to the HWM14 empirical model's other
    climatological and observational uncertainties.  The results from gfortran
    (v4.8.1) and the ifort (intel) compiler without optimizations (-00) are
    identical. The default optimization level for the intel fortran compiler is
    '-O2', so example reaults from 'ifort -O0' and 'ifort' differ slightly.
    See included files Check/ifort.nopt.txt (ifort -O0), Check/gfortran.txt
    (gfortran), and Check/ifort.opt.txt (ifort -O2).

RELEASE NOTES

  HWM14 contains several updates to HWM07 and HWM93
  Improvements to the code include:

      - A more efficient routine (ALFBASIS) for computing associated
        Legendre functions and the latitude portion of vector spherical
        harmonic functions.

      - A new routine for the conversion of geodetic to Quasi-Dipole
        coordinates for the disturbance wind component. This is now done using
        a spherical harmonic expansion of QD coordinates, rather than gridded
        interpolation. This modification results in slightly different
        disturbance winds than in the previous version, because the
        representation is not exact. The difference is less than 1 m/s.

      - The code has been streamlined and combined into a single file.

      - The FORTRAN90 module layout has been reorganized to be compatable
        with the python NUMPY f2py feature.

  HWM14 requires three binary data files: hwm123114.bin,
  dwm07b104i.dat, and gd2qd.dat. These file should be located in one of the
  following locations; 1) directory from which the program is run, 2) a
  subdirectory "../Meta/" with respect to the
  directory from which the program is run; and/or 3) in a directory
  defined by the shell variable HWMPATH. See the subroutine findandopen()
  for additional details.

MODEL FORMULATION
  This is an update to the HWM07 empirical model of horizontal winds in the
  troposphere, stratosphere, mesosphere, and thermosphere, [Drob, D. P., et al.
  (2008), An Empirical Model of the Earth's Horizontal Wind Fields: HWM07, J.
  Geophys Res., 113, doi:10.1029/2008JA013668.] which also supersedes  HWM93
  [Hedin et al., J. Atmos. Terr. Phys., vol. 58, 1421-1447, 1996]. In addition
  to the date used in HWM93, the HWM07 used extensive new ground-based and space-
  based wind measurements, including height profiles from NASA-UARS/WINDII,
  NASA-UARS/HRDI, measurements from ground-based optical and radar instruments
  obtained from the NSF-CEDAR database, and lower atmospheric NCEP data. This
  model (HWM14) now includes additional ground-based 630-nm FPI wind
  measurements at both equatorial and high latitudes (~250 km altitude), in
  addition to cross track accelerometer wind measurements from the GOCE
  satellite (~250 km altitude).
     As in HWM07, in the thermosphere, the model consists of two parts: a quiet-time
  portion, and a geomagnetically disturbed portion. The quiet-time part
  represents average wind conditions when ap <= 12. The disturbed part
  represents average perturbation winds for the specified ap input.
      The quiet part is represented by vector spherical harmonics in geodetic
  latitude, geodetic longitude, and solar local time, up to wave number 8 in
  latitude, 2 in longitude, and 3 in local time. The seasonal dependence is
  represented by harmonic terms up to semiannual. The vertical structure is
  represented below 250 km by cubic B-splines with node spacing of 5 km below
  110 km and higher nodes at 110, 117, 130, 140, 155, 200, and 250 km (modified
  from HWM07). Above 250 km, an exponential decay function with a scale height
  of 60 km is used; continuity up to the second derivative is imposed at 250 km.
      The disturbance winds depend on magnetic latitude, magnetic local time,
  and Kp. The Quasi-Dipole magnetic coordinates described by Richmond
  [J Geomagn. Geoelectr., vol. 47, 191-212, 1995] are used for the magnetic
  coordinates; the code was obtained from the NSF-CEDAR database, and the
  interpolation grid was computed from IGRF at 250 km and epoch 1994.0 (a
  modified FORTRAN-90 version of the code, apexcord.f90, is included in this
  package; it contains only those subroutines needed for reading the
  interpolation grid and computing QD coordinates). The magnetic latitude and
  magnetic local time dependence of the disturbance winds is represented by
  vector spherical harmonics up to wave number 10 in magnetic latitude and wave
  number 3 in magnetic local time. At mid and low latitudes, only latitudinal
  terms up to wave number 4 are used; the transition from low resolution at low
  latitudes to high resolution at high latitudes occurs at a pre-determined
  latitude that depends on local time and Kp. The transition is made with an
  exponential function with a width of 4 degrees. The Kp dependence is
  represented by cubic splines with nodes at 0, 2, 5, and 8. The Kp dependence
  is constrained to have zero slope at Kp=0 and Kp=8, and is constant above
  Kp=8.


HWM07 ARGUMENT CHANGES FROM HWM93
  The input and output argument list remains unchanged from that of HWM93,
  however some of specifics of the arguments have been modified as follows:
   - The STL (solar local time) input argument is ignored. Local time is
     computed from the SEC (UT) and GLON (geographic longitude) input arguments.
   - Only the second element of the AP argument (the 3-hour ap) is used to
     compute disturbance winds.
   - The F107 and F107A arguments are ignored in this version.
   - Module component on/off switches are not yet fully implemented in this
     version. There is no longer a TSELEC subroutine. For now, for example if
     longitudinally averaged winds are desired, the model must be evaluated on
     a longitude grid and the output averaged by the user.

HWM14 ARGUMENT CHANGES FROM HWM07
   - None

MODEL LIMITATIONS AND FUTURE IMPROVEMENTS
   - The model currently contains no solar activity dependence, and the F107 and
     F107A arguments are ignored. During the day, solar activity has a
     relatively insignificant effect on thermospheric winds. At night, however,
     the effect is signficant. Solar activity dependencies will be included in
     a future version.
   - The disturbed part depends only on magnetic latitude, magnetic local time,
     and Kp (via the ap argument), and represents average disturbance winds in
     the upper thermosphere (above 225 km). The disturbance winds are assumed
     to be constant with height, with a smooth artificial cutoff below 125 km.
     Height, seasonal, and solar activity dependences will be included in a
     future version.


TO GET TOTAL WINDS:
  Call HWM14 with the following input arguments:
        IYD - YEAR AND DAY AS YYDDD
        SEC - UT(SEC)
        ALT - ALTITUDE(KM)
        GLAT - GEODETIC LATITUDE(DEG)
        GLON - GEODETIC LONGITUDE(DEG)
        STL - Not used
        F107A - Not used
        F107 - Not used
        AP - Two element array with
             AP(1) = Not used
             AP(2) = CURRENT 3HR ap INDEX
  The output argument is
        W(1) = MERIDIONAL WIND (m/sec + Northward)
        W(2) = ZONAL WIND (m/sec + Eastward)


TO GET QUIET TIME WINDS:
  Call HWM14 with a negative value for AP(2)


TO GET DISTURBANCE WINDS IN GEOGRAPHIC COORDINATES:
  Call DWM07 with the following input arguments:
        IYD - YEAR AND DAY AS YYDDD
        SEC - UT(SEC)
        ALT - ALTITUDE(KM)
        GLAT - GEODETIC LATITUDE(DEG)
        GLON - GEODETIC LONGITUDE(DEG)
        AP - Two element array with
             AP(1) = Not used
             AP(2) = CURRENT 3HR ap INDEX
  The output argument is
        DW(1) = MERIDIONAL DISTURBANCE WIND (m/sec + Geo. Northward)
        DW(2) = ZONAL DISTURBANCE WIND (m/sec + Geo. Eastward)


TO GET DISTURBANCE WINDS IN MAGNETIC COORDINATES:
  Call DWM07b with the following input arguments:
        MLT - MAGNETIC LOCAL TIME (HRS)
        MLAT - MAGNETIC LATITUDE(DEG)
        KP - CURRENT 3HR Kp INDEX
  The output arguments are
        MMPWIND = MERIDIONAL DISTURBANCE WIND (m/sec + Magnetic Northward)
        MZPWIND = ZONAL DISTURBANCE WIND (m/sec + Magnetic Eastward)
