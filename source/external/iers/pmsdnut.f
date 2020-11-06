      subroutine PMsdnut (rmjd, pm)
c
c*******************************************************************
c   Fortran subroutine to evaluate the model of polar motion for
c   a nonrigid Earth due to tidal gravitation. This polar motion
c   is equivalent to the so-called "subdiurnal nutation". The model
c   is a sum of a first order polynomial and 25 trigonometric terms
c   (15 long periodic and 10 quasi diurnal) with coefficients given
c   in Table 5.1 of the IERS Conventions (2003).
c
c********IMPORTANT*********
c   In the present version this subroutine neglects the linear trend
c   and the long periodic terms of the expansion, for the reasons
c   explained in Section 5.4.2 of the IERS Conventions (2003), last
c   paragraph before Table 5.1. If the full expansion is needed set
c   the parameter iband to 0 instead of 1, that is replace the statement
c      parameter ( iband = 1 )
c   below by
c      parameter ( iband = 0 )
c*************************
c
c   Written by Aleksander Brzezinski,
c              Space Research Centre PAS, Warsaw, March 2005
c
c INPUT:
c   rmjd    r*8  time expressed as modified Julian date
c
c OUTPUT:
c   pm      r*8  vector of length 2 with the incremental polar motion
c                coordinates (dx,dy) expressed in microarcseconds
c
c  External calls: PMargs
c*******************************************************************
c
c Declarations
      real*8 rmjd, pm(2)
c
c  iband  - parameter defining the range of periods for the terms which
c           are included in comutations; if equal to 1 only the quasi
c           diurnal terms are computed, otherwise the full model
      integer iband
      parameter ( iband = 1 )
c
c  iarg   - array defining for each of the 25 trigonometric terms a set
c           of 6 integer multipliers of the fundamental angular arguments
c  arg    - vector of the following 6 fundamental arguments used to
c           compute the angular argument of the trigonometric functions
c           arg(1:6) = [ GMST+pi, el, elp, f, d, om ]; this vector is
c           evaluated by the subroutine PMargs which is enclosed below
c  period - array of periods of the trigonometric terms of expansion, in
c           mean solar days; only for a check - not used in computations
c  xs, xc - sine and cosine coefficients of the x coordinate of the pole,
c           in microarcseconds
c  ys, yc - sine and cosine coefficients of the y coordinate of the pole,
c           in microarcseconds
c  angle  - angular argument of the trigonometric functions
c           angle = Sum(i=1:6) iarg(i,j)*arg(i), for j=1,25
c
      integer iarg(6,25)
      real*8 arg(6)
      real*8 per(25), xs(25), xc(25), ys(25), yc(25)
      real*8 angle
c
c Set constants
c
c  rmjd0   - modified Julian date of J2000
c  twopi   - 2*pi
c
      real*8 rmjd0, twopi
      parameter ( rmjd0   = 51544.5d0                )
      parameter ( twopi   = 6.28318530717958647692d0 )
c
c Coefficients of the long periodic terms in polar motion
c Source: IERS Conventions 2003, Table 5.1
c
      data
     * (  (iarg(i,j),i=1,6),    per(j),   xs(j),  xc(j),  ys(j),  yc(j),
     *                                                           j=1,15)
     */ 0,  0, 0,  0,  0, -1, 6798.3837,  -0.03,   0.63,  -0.05,  -0.55,
     *  0, -1, 0,  1,  0,  2, 6159.1355,   1.46,   0.00,  -0.18,   0.11,
     *  0, -1, 0,  1,  0,  1, 3231.4956, -28.53,  -0.23,   3.42,  -3.86,
     *  0, -1, 0,  1,  0,  0, 2190.3501,  -4.65,  -0.08,   0.55,  -0.92,
     *  0,  1, 1, -1,  0,  0, 438.35990,  -0.69,   0.15,  -0.15,  -0.68,
     *  0,  1, 1, -1,  0, -1, 411.80661,   0.99,   0.26,  -0.25,   1.04,
     *  0,  0, 0,  1, -1,  1, 365.24219,   1.19,   0.21,  -0.19,   1.40,
     *  0,  1, 0,  1, -2,  1, 193.55971,   1.30,   0.37,  -0.17,   2.91,
     *  0,  0, 0,  1,  0,  2, 27.431826,  -0.05,  -0.21,   0.01,  -1.68,
     *  0,  0, 0,  1,  0,  1, 27.321582,   0.89,   3.97,  -0.11,  32.39,
     *  0,  0, 0,  1,  0,  0, 27.212221,   0.14,   0.62,  -0.02,   5.09,
     *  0, -1, 0,  1,  2,  1, 14.698136,  -0.02,   0.07,   0.00,   0.56,
     *  0,  1, 0,  1,  0,  1, 13.718786,  -0.11,   0.33,   0.01,   2.66,
     *  0,  0, 0,  3,  0,  3, 9.1071941,  -0.08,   0.11,   0.01,   0.88,
     *  0,  0, 0,  3,  0,  2, 9.0950103,  -0.05,   0.07,   0.01,   0.55/
c
c Coefficients of the quasi diurnal terms in polar motion
c Source: IERS Conventions 2003, Table 5.1
c
      data
     *(  (iarg(i,j),i=1,6),     per(j),   xs(j),  xc(j),  ys(j),  yc(j),
     *                                                          j=16,25)
     */ 1, -1, 0, -2,  0, -1, 1.1196992,  -0.44,   0.25,  -0.25,  -0.44,
     *  1, -1, 0, -2,  0, -2, 1.1195149,  -2.31,   1.32,  -1.32,  -2.31,
     *  1,  1, 0, -2, -2, -2, 1.1134606,  -0.44,   0.25,  -0.25,  -0.44,
     *  1,  0, 0, -2,  0, -1, 1.0759762,  -2.14,   1.23,  -1.23,  -2.14,
     *  1,  0, 0, -2,  0, -2, 1.0758059, -11.36,   6.52,  -6.52, -11.36,
     *  1, -1, 0,  0,  0,  0, 1.0347187,   0.84,  -0.48,   0.48,   0.84,
     *  1,  0, 0, -2,  2, -2, 1.0027454,  -4.76,   2.73,  -2.73,  -4.76,
     *  1,  0, 0,  0,  0,  0, 0.9972696,  14.27,  -8.19,   8.19,  14.27,
     *  1,  0, 0,  0,  0, -1, 0.9971233,   1.93,  -1.11,   1.11,   1.93,
     *  1,  1, 0,  0,  0,  0, 0.9624365,   0.76,  -0.43,   0.43,   0.76/
c
c Rate of secular polar motion, in microarcseconds per year
c Source: IERS Conventions 2003, Table 5.1
c
      data xrate, yrate / -3.80, -4.31/
c
c Compute the periodical part of the model
c Coordinates of the pole are set to zero first
      pm(1) = 0.d0
      pm(2) = 0.d0
c
c Evaluate the vector of the fundamental arguments
c arg(1:6) = [ GMST+pi, el, elp, f, d, om ] at t = rmjd
      call PMargs (rmjd,arg)
c
      if (iband.eq.1) then
        jstart = 16
      else
        jstart = 1
      endif
      do 20 j=jstart,25
c For the j-th term of the trigonometric expansion, compute the angular
c argument angle of sine and cosine functions as a linear integer
c combination of the 6 fundamental arguments
        angle = 0.d0
        do 10 i=1,6
          angle = angle + iarg(i,j) * arg(i)
   10   continue
        angle = dmod( angle, twopi)
c
c Compute contribution from the j-th term to the polar motion coordinates
        pm(1) = pm(1) + xs(j)*dsin(angle) + xc(j)*dcos(angle)
        pm(2) = pm(2) + ys(j)*dsin(angle) + yc(j)*dcos(angle)
   20 continue
      if (iband.eq.1) return
c
c Add the secular term of the model
      pm(1) = pm(1) + xrate * (rmjd-rmjd0) / 365.25d0
      pm(2) = pm(2) + yrate * (rmjd-rmjd0) / 365.25d0
c
      return
      end
c
c
c
      subroutine PMargs (rmjd, pmarg)
c
c*******************************************************************
c   Fortran subroutine to compute the angular arguments used in
c   the trigonometric expansion of the model of polar motion for
c   a nonrigid Earth due to tidal gravitation, as describe by
c   Table 5.1 of the IERS Conventions (2003).
c
c   Written by Aleksander Brzezinski,
c              Space Research Centre PAS, Warsaw, March 2005
c
c INPUT:
c   rmjd    r*8  time expressed as modified Julian date
c
c OUTPUT:
c   pmarg   r*8  vector of length 6 with the following angular arguments
c                expressed in radians:
c   1.           GMST (Greenwich mean sidereal time) + pi
c   2.           el  - mean anomaly of the Moon
c   3.           elp - mean anomaly of the Sun
c   4.           f   - L minus om, where L is the mean longitude
c                      of the Moon and om is defined below
c   5.           d   - mean elongation of the Moon from the Sun
c   6.           om  - mean longitude of the ascending node of the Moon
c
c  There is no external calls in this subroutine
c*******************************************************************
c
c Declarations
      real*8 rmjd, pmarg(6), GMST, el, elp, f, d, om
c
c  GMSTc, elc, elpc, fc, dc, omc - vectors of length 5 containing
c                        coefficients of the polynomial expansion
c                        of GMST, el, elp, f, d, om
c  cent                - epoch in Julian centuries of TDB since J2000
c  cent2, cent3, cent4 - powers of cent
      real*8 GMSTc(5), elc(5), elpc(5), fc(5), dc(5), omc(5)
      real*8 cent, cent2, cent3, cent4
c Constants
c
c   pi, twopi   - pi, 2*pi
c   rmjd0       - modified Julian date of J2000
c   sec360      - number of arcseconds in 360 degrees
c   rad2deg     - conversion from radians to degrees
c   rad2sec     - conversion from radians to time seconds
c
      real*8 pi, twopi, rmjd0, sec360, rad2deg, rad2sec
c
      parameter ( pi      = 3.14159265358979323846 d0 )
      parameter ( twopi   = 6.28318530717958647692 d0 )
      parameter ( rmjd0   = 51544.5d0                 )
      parameter ( sec360  = 1296000.d0                )
      parameter ( rad2deg = 180.d0/pi                 )
      parameter ( rad2sec = 86400.d0/twopi            )
c
c  Coefficients of the polynomial expansion of GMST (in seconds) and of
c  the fundamental arguments of nutation (in arcseconds). Source:
c  IERS Conventions (1996), Chapter 5 for GMST
c  IERS Conventions (2003), eq.(40) for the nutation arguments
c
      data GMSTc  / -0.0000062d0,  0.093104d0,  8640184.812866d0,
     *         3155760000.d0,   67310.54841d0 /
      data elc    /    -0.00024470d0,     0.051635d0,  31.8792d0,
     *         1717915923.2178d0,   485868.249036d0 /
      data elpc   /    -0.00001149d0,    +0.000136d0,  -0.5532d0,
     *          129596581.0481d0,   1287104.79305d0 /
      data fc     /     0.00000417d0,    -0.001037d0,  -12.7512d0,
     *         1739527262.8478d0,   335779.526232d0 /
      data dc     /    -0.00003169d0,     0.006593d0,   -6.3706d0,
     *         1602961601.2090d0,   1072260.70369d0 /
      data omc    /    -0.00005939d0,     0.007702d0,    7.4722d0,
     *           -6962890.5431d0,   450160.398036d0 /
c
c  Convert the input epoch to Julian centuries of TDB since J2000,
c  and compute its powers
c
      cent = (rmjd-rmjd0)/36525.d0
      cent2 = cent*cent
      cent3 = cent2*cent
      cent4 = cent3*cent
c
c  Evaluate the polynomial expansions of GMST and of
c  the fundamental arguments of nutation
c
      GMST = GMSTc(1)*cent3 + GMSTc(2)*cent2
     *      + (GMSTc(3) + GMSTc(4))*cent + GMSTc(5)
      GMST = dmod( GMST, 86400.d0 )
c
      el = elc(1)*cent4 + elc(2)*cent3 + elc(3)*cent2
     *    + elc(4)*cent + elc(5)
      el = dmod( el, sec360 )
c
      elp = elpc(1)*cent4 + elpc(2)*cent3 + elpc(3)*cent2
     *     + elpc(4)*cent + elpc(5)
      elp = dmod( elp, sec360 )
c
      f = fc(1)*cent4 + fc(2)*cent3 + fc(3)*cent2
     *   + fc(4)*cent + fc(5)
      f = dmod( f, sec360 )
c
      d = dc(1)*cent4 + dc(2)*cent3 + dc(3)*cent2
     *   + dc(4)*cent + dc(5)
      d = dmod( d, sec360 )
c
      om = omc(1)*cent4 + omc(2)*cent3 + omc(3)*cent2
     *    + omc(4)*cent + omc(5)
      om = dmod( om, sec360 )
c
c  Convert the arguments to radians
c
      pmarg(1) = GMST / rad2sec + pi
      pmarg(1) = dmod( pmarg(1), twopi )
      pmarg(2) = el   / (3600.d0*rad2deg)
      pmarg(3) = elp  / (3600.d0*rad2deg)
      pmarg(4) = f    / (3600.d0*rad2deg)
      pmarg(5) = d    / (3600.d0*rad2deg)
      pmarg(6) = om   / (3600.d0*rad2deg)
c
      return
      end

