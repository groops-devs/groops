===============================================================================
MSIS® (NRL-SOF-014-1) SOFTWARE

MSIS® is a registered trademark of the Government of the United States of 
America, as represented by the Secretary of the Navy. Unauthorized use of 
the trademark is prohibited. 

The MSIS® Software (hereinafter Software) is property of the United States 
Government, as represented by the Secretary of the Navy. Methods performed
by this software are covered by U.S. Patent Number 10,641,925. The Government
of the United States of America, as represented by the Secretary of the Navy, 
herein grants a non-exclusive, non-transferable license to the Software for 
academic, non-commercial, purposes only. A user of the Software shall not: 
(i) use the Software for any non-academic, commercial purposes, (ii) make 
any modification or improvement to the Software, (iii) disseminate the 
Software or any supporting data to any other person or entity who will use 
the Software for any non-academic, commercial purposes, or (iv) copy the 
Software or any documentation related thereto except for (a) distribution 
among the user’s personal computer systems, archival, or emergency repair 
purposes, or (b) distribution for non-commercial, academic purposes, without 
first obtaining the written consent of IP Counsel for the Naval Research 
Laboratory. 

As the owner of MSIS®, the United States, the United States Department of 
Defense, and their employees: (1) Disclaim any warranties, express, or 
implied, including but not limited to any implied warranties of 
merchantability, fitness for a particular purpose, title or non-infringement, 
(2) Do not assume any legal liability or responsibility for the accuracy, 
completeness, or usefulness of the software, (3) Do not represent that use of 
the software would not infringe privately owned rights, (4) Do not warrant 
that the software will function uninterrupted, that is error-free or that any 
errors will be corrected.

BY USING THIS SOFTWARE YOU ARE AGREEING TO THE ABOVE TERMS AND CONDITIONS.  
===============================================================================

NRLMSIS 2.0 Whole-Atmosphere Empirical Model of Temperature and Neutral Species
  Densities

This software package is a major upgrade to the NRLMSISE-00 model.

VERSION HISTORY
  08 MAR 19 Version 1.97 (Beta version)
  27 APR 20 Version 1.99 (Beta version)
  26 MAY 20 Version 2.0 (Release version)

AUTHORS
  Douglas Drob (douglas.drob@nrl.navy.mil)
  John Emmert (john.emmert@nrl.navy.mil)
  
REFERENCE
  Emmert, J.T., Drob, D. P., Picone, J. M., Siskind, D. E., Jones Jr., M.,
  et al. (2020). NRLMSIS 2.0: A whole-atmosphere empirical model of temperature
  and neutral species densities. Manuscript in preparation, to be submitted to
  Earth and Space Science.

PACKAGE CONTENTS
  readme.txt               This file
  msis2.0_test.F90         Test program
  msis_init.F90            Subroutines to initialize the model, set switches and
                             options, and load parameter file
  msis_gtd8d.F90           Subroutine to evaluate the model using the legacy
                             interface
  msis_calc.F90            Subroutine for evaluating the model using the new
                             interface
  msis_constants.F90       Module containing model constants
  msis_gfn.F90             Subroutines to calculate horizontal expansion
                             functions
  msis_tfn.F90             Subroutines to calculate the vertical temperature
                             profile
  msis_dfn.F90             Subroutines to calculate vertical density profiles
  alt2gph.F90              Subroutines to convert between geodetic height and
                             geopotential height
  subroutine_dir.jpg       Graphical directory of modules and subroutines
                             contained in the source code files
  msis2.0.parm             Binary data file containing model parameters
  msis2.0_test_in.txt      ASCII file containing input for test program.
  msis2.0_test_ref_dp.txt  ASCII file containing expected output of test program
                             (double-precision internally)
  msis2.0_test_ref_sp.txt  ASCII file containing expected output of test program
                             (single-precision internally)
 
RELEASE NOTES: MODEL FORMULATION
  Major changes to the NRLMSISE-00 formulation include:
  - The transition from a fully mixed atmosphere to diffusive separation is now
    represented via height-dependent effective mass for each species.
  - The temperature profile is now C2 continuous.
  - A global geopotential height function is now used internally.
  - Atomic oxygen now extends down to 50 km; below 85 km, the O profile is
    represented by cubic B splines decoupled from temperature.
  - Thermal diffusion is no longer applied to any species.
 
RELEASE NOTES: PARAMETER ESTIMATION
  - The parameters of the model were tuned to extensive new data in the
    troposphere, stratosphere, and mesosphere, and to NRLMSISE-00 in the
    thermosphere. In the thermosphere, the model atomic oxygen was additionally
    tuned to global average orbit-derived density data.
  - Thermospheric N2 is now controlled entirely by the temperature profile and
    the constant mixing ratio of N2 in the lower atmosphere. See paper for
    details. To instead retrieve the NRLMSISE-00 thermospheric N2, set 
    lN2_msis00 = .true. in the call to MSISINIT.

COMPILING THE MODEL CODE
  The model package was tested on Windows, Linux, and Mac systems using the
    following Fortran compilers and compile statements:
    gfortran 4.8.5, 6.3.0, 8.1.0
      gfortran -O3 -cpp -o msis2.0_test.exe alt2gph.F90 msis_constants.F90
      msis_init.F90 msis_gfn.F90 msis_tfn.F90 msis_dfn.F90 msis_calc.F90
      msis_gtd8d.F90 msis2.0_test.F90
        NOTES:
        - For double precision, add the flag -DDBLE
        - The following optimization flags may improve performance:
            -march=native -ffast-math
    Intel 17.0.2.163, 18.0.1.156
      ifort -O2 -fpp -o msis2.0_test.exe alt2gph.F90 msis_constants.F90
      msis_init.F90 msis_gfn.F90 msis_tfn.F90 msis_dfn.F90 msis_calc.F90
      msis_gtd8d.F90 msis2.0_test.F90
        NOTES:
        - For double precision, add the flag -DDBLE
        - The following optimization flags may improve performance:
            Windows:     -Qipo -QxHost
            Linux/macOS: -ipo -xHost

INITIALIZING AND RUNNING THE MODEL
  - The model must be initialized using the MSISINIT subroutine, which sets
    switches and options and loads the model parameter values from a file.
  - The switch_legacy optional argument to MSISINIT performs the same function
    as TSELEC(SW) in NRLSMSISE-00, except that switches 15-25 are not used in
    NRLMSIS 2.0. The change in the switch-setting call is illustrated as
    follows, where SW is the 25-element array of switches:
        NRLMSISE-00: CALL TSELEC(SW)
        NRLMSIS 2.0: call msisinit(switch_legacy=SW)
  - The MSISCALC subroutine checks for initialization and does a default
    initialization if necessary. This self-initialization will be removed in
    future versions.
  - The model can be called using either the legacy interface (subroutine
    GTD8D) or the new interface (subroutine MSISCALC).
  - Details of the input and output arguments of MSISINIT, GTD8D, and MSISCALC
    are provided in the headers of the respective source code files.
  
COMPUTATIONAL PERFORMANCE
  - The model should be initialized (by calling MSISINIT) only before first use
    or when switches and options need to be changed. Re-initializing before each
    model call will significantly slow down the code.
  - The code is most efficient when computing vertical profiles at the same
    horizontal location and time, in which case the computed vertical profile
    parameters are re-used from call to call. Therefore, altitude should be in
    the inner loop of multiple calls, when possible.
  - For "column" calls (inner loop on altitude) with 100 vertical levels, 
    NRLMSIS 2.0 is 4-5 times faster than NRLMSISE-00.
  - For applications in which the horizontal location or time changes with every
    call (e.g., satellite ephemerides), users may find that NRLMSIS 2.0 is 1/3
    to 1/2 as fast as NRLMSISE-00. The reduction in speed is due to the greater
    number of model parameters in NRLMSIS 2.0, as well its stronger coupling
    with the lower atmosphere, compared to NRLMSISE-00. However, the model is
    still very fast: a typical desktop system can process at least 100,000
    serial calls per second.
  - Additional improvement in speed can be obtained by computing only those
    species that are needed, via the lspec_select input argument of MSISINIT.

ACKNOWLEDGEMENTS
  This work was supported by the Chief of Naval Research, the Office of Naval
  Research (BSION program), and NASA (Grants NNH16ZDA001N-HSR/ITM16_2-0013, 
  NNH14ZDA001N-GIODDE14_2/NNH15AZ72I)