C
      SUBROUTINE JB2008 (AMJD,SUN,SAT,F10,F10B,S10,S10B,XM10,XM10B,
     *                   Y10,Y10B,DSTDTC,TEMP,RHO)
C
C     Jacchia-Bowman 2008 Model Atmosphere
C
C
C***********************************************************************
C
C     This is the CIRA "Integration Form" of a Jacchia Model.
C     There are no tabular values of density.  Instead, the barometric
C     equation and diffusion equation are integrated numerically using
C     the Newton-Coates method to produce the density profile up to the
C     input position.
C
C     INPUT:
C
C           AMJD   : Date and Time, in modified Julian Days
C                    and Fraction (MJD = JD-2400000.5)
C           SUN(1) : Right Ascension of Sun (radians)
C           SUN(2) : Declination of Sun (radians)
C           SAT(1) : Right Ascension of Position (radians)
C           SAT(2) : Geocentric Latitude of Position (radians)
C           SAT(3) : Height of Position (km)
C           F10    : 10.7-cm Solar Flux (1.0E-22*Watt/(M**2*Hertz))
C                    (Tabular time 1.0 day earlier)
C           F10B   : 10.7-cm Solar Flux, ave.
C                    81-day centered on the input time
C                    (Tabular time 1.0 day earlier)
C           S10    : EUV index (26-34 nm) scaled to F10
C                    (Tabular time 1.0 day earlier)
C           S10B   : EUV 81-day ave. centered index
C                    (Tabular time 1.0 day earlier)
C           XM10   : MG2 index scaled to F10
C                    (Tabular time 2.0 days earlier)
C           XM10B  : MG2 81-day ave. centered index
C                    (Tabular time 2.0 days earlier)
C           Y10    : Solar X-Ray & Lya index scaled to F10
C                    (Tabular time 5.0 days earlier)
C           Y10B   : Solar X-Ray & Lya 81-day ave. centered index
C                    (Tabular time 5.0 days earlier)
C           DSTDTC : Temperature change computed from Dst index
C
C     OUTPUT:
C
C           TEMP(1): Exospheric Temperature above Input Position (deg K)
C           TEMP(2): Temperature at Input Position (deg K)
C           RHO    : Total Mass-Desnity at Input Position (kg/m**3)
C
C
C     JB2008 Model Development: (Ref. 7)
C
C
C     A. Development of the JB2006 model:
C
C       1. Started with the CIRA72 model (Jacchia 71).
C
C       2. Converted to CIRA70 model replacing equations from Jacchia 70
C          model (Ref. 5)
C
C       3. Replaced Tc equation using new solar indices (Ref. 1 and 2)
C
C       4. Replaced semiannual equation with new global model based
C          on F10B (Ref. 1 and 3)
C
C       5. Added correction for local solar time and latitude errors
C          (Ref. 1)
C          Added smooth transition between altitude bands
C
C       6. Added high altitude ( z > 1500 km ) correction
C          (Ref. 1 and 4)
C
C       7. REV A of JB2006 - Oct 2006
C                Smoothing of density corrections and scale height
C                through different altitude bands in the latitude-
C                local time correction subroutine DTSUB
C                dTx correction replaced with dTc correction
C
C     B. Modification to develop JB2008 model:
C
C       1. Replaced Tc equation in JB2006 using new solar indices
C          (Ref. 7)
C
C       2. Replaced semiannual equation with new global model based
C          on F10B and S10B (Ref. 6)
C
C       3. Use dTc value based on Dst geomagnetic storm index
C          (This replaces ap use) (Ref. 7)
C
C
C         All equation references below refer to the original
C         Jacchia 1971 (CIRA 1972) model papers.
C
C
C     References:
C
C      1. Bowman, Bruce R., etc. : "A New Empirical Thermospheric
C         Density Model JB2006 Using New Solar Indices",
C         AIAA/AAS Astrodynamics Specialists Conference, Keystone, CO,
C         21-24 Aug 2006, (Paper AIAA 2006-6166).
C
C      2. Bowman, Bruce R., etc. : "Improvements in Modeling
C         Thermospheric Densities Using New EUV and FUV Solar Indices",
C         AAS/AIAA Space Flight Mechanics Meeting, Tampa, FL,
C         23-26 Jan 2006, (Paper AAS 06-237).
C
C      3. Bowman, Bruce R.: "The Semiannual Thermospheric Density
C         Variation From 1970 to 2002 Between 200-1100 km",
C         AAS/AIAA Space Flight Mechanics Meeting, Maui, HI,
C         8-12 Feb 2004, (Paper AAS 04-174).
C
C      4. Bowman, Bruce R.; "Atmospheric Density Variations at
C         1500 km to 4000 km Height Determined from Long Term
C         Orbit Perturbation Analysis", AAS/AIAA Space Flight
C         Mechanics Meeting, Santa Barbara, CA, 11-14 Feb 2001,
C         (Paper AAS 01-132).
C
C      5. Jacchia, Luigi G.; "New Static Models of the
C         Thermosphere and Exosphere with Empirical Temperature
C         Profiles", (Smithsonian Astrophysical Observatory
C         Special Report 313), 6 May 1970.
C
C      6. Bowman, Bruce R., etc. : "The Thermospheric Semiannual Density
C         Response to Solar EUV Heating," JASTP, 2008
C
C      7. Bowman, Bruce R., etc. : "A New Empirical Thermospheric
C         Density Model JB2008 Using New Solar and Geomagnetic Indices",
C         AIAA/AAS 2008, COSPAR CIRA 2008 Model
C
C
C
C     Written by: Bruce R Bowman (HQ AFSPC, Space Analysis Division),
C                 2008
C
C
C
C***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)

      DIMENSION SUN(2), SAT(3), TEMP(2), AL10N(6), CHT(4)
      DIMENSION ALN(6), ALPHA(5), AMW(6), FRAC(4), TC(4), WT(5)

C     The alpha are the thermal diffusion coefficients in Eq. (6)

      DATA ALPHA(1) /0.D0/
      DATA ALPHA(2) /0.D0/
      DATA ALPHA(3) /0.D0/
      DATA ALPHA(4) /0.D0/
      DATA ALPHA(5) /-0.38D0/

C     AL10 is DLOG(10.0)

      DATA AL10 /2.3025851D0/

C     The AMW are the molecular weights in order: N2, O2, O, Ar, He & H

      DATA AMW(1) /28.0134D0/
      DATA AMW(2) /31.9988D0/
      DATA AMW(3) /15.9994D0/
      DATA AMW(4) /39.9480D0/
      DATA AMW(5) /4.0026D0/
      DATA AMW(6) /1.00797D0/

C     AVOGAD is Avogadro's number in mks units (molecules/kmol)

      DATA AVOGAD /6.02257D26/

      DATA FOURPI /12.566371D0/
      DATA TWOPI  /6.2831853D0/
      DATA PI     /3.1415927D0/
      DATA PIOV2  /1.5707963D0/
      DATA PIOV4  /0.78539816D0/

C     The FRAC are the assumed sea-level volume fractions in order:
C     N2, O2, Ar, and He

      DATA FRAC(1) /0.78110D0/
      DATA FRAC(2) /0.20955D0/
      DATA FRAC(3) /9.3400D-3/
      DATA FRAC(4) /1.2890D-5/

C     RSTAR is the universal gas-constant in mks units (joules/K/kmol)

      DATA RSTAR /8314.32D0/

C     The R# are values used to establish height step sizes in
C     the regimes 90km to 105km, 105km to 500km and 500km upward.

      DATA R1 /0.010D0/
      DATA R2 /0.025D0/
      DATA R3 /0.075D0/

C     The WT are weights for the Newton-Cotes Five-Point Quad. formula

      DATA WT(1) /0.311111111111111D0/
      DATA WT(2) /1.422222222222222D0/
      DATA WT(3) /0.533333333333333D0/
      DATA WT(4) /1.422222222222222D0/
      DATA WT(5) /0.311111111111111D0/

C     The CHT are coefficients for high altitude density correction

      DATA CHT /0.22D0,-0.20D-02,0.115D-02,-0.211D-05/
C
      DEGRAD  =   PI / 180.D0
C
      IDEBUG = 0
C
      IF (IDEBUG.EQ.1) THEN
        WRITE(8,101) (SAT(I),I=1,3),F10,F10B,XM10B
 101    FORMAT(' JB2008 - SAT,F10,F10B ',3F20.6,3F8.1)
      ENDIF
C
C     Equation (14)
C
      FN = (F10B/240.)**(1./4.)
      IF (FN.GT.1.) FN = 1.0
      FSB = F10B*FN + S10B*(1.-FN)
      TSUBC = 392.4D0 + 3.227D0*FSB + 0.298D0*(F10-F10B)
     *                + 2.259D0*(S10-S10B) + 0.312D0*(XM10-XM10B)
     *                + 0.178D0*(Y10-Y10B)
C
      IF (IDEBUG.EQ.1) THEN
        WRITE(8,130) TSUBC,F10B,S10B,XM10B
 130    FORMAT(' JB2008 - TSUBC',4F15.2)
      ENDIF
C
C     Equation (15)

      ETA =   0.5D0 * DABS(SAT(2) - SUN(2))
      THETA = 0.5D0 * DABS(SAT(2) + SUN(2))

C     Equation (16)

      H = SAT(1) - SUN(1)
      TAU = H - 0.64577182D0 + 0.10471976D0 * DSIN(H + 0.75049158D0)
C
      GLAT  = SAT(2)
      ZHT   = SAT(3)
      GLST  = H + PI
      GLSTHR = (GLST/DEGRAD)*(24.D0/360.D0)
      IF (GLSTHR.GE.24.D0) GLSTHR = GLSTHR - 24.D0
      IF (GLSTHR.LT. 0.D0) GLSTHR = GLSTHR + 24.D0
C
      IF (IDEBUG.EQ.1) THEN
        WRITE(8,131) ETA,THETA,TAU
 131    FORMAT(' JB2008 - ETA,THETA,TAU',3D25.12)
      ENDIF
C
C     Equation (17)

      C = DCOS(ETA)**2.5
      S = DSIN(THETA)**2.5
C
      IF (IDEBUG.EQ.1) THEN
        WRITE(8,132) C,S
 132    FORMAT(' JB2008 - C,S',2D25.12)
      ENDIF
C
      DF = S + (C - S) * DABS(DCOS(0.5 * TAU))**3
      TSUBL = TSUBC * (1.D0 + 0.31D0 * DF)

C     Compute correction to dTc for local solar time and lat correction
C
      CALL DTSUB (F10,GLSTHR,GLAT,ZHT,DTCLST)
C
      IF (IDEBUG.EQ.2) THEN
        WRITE(8,133) SAT(3),DTCLST
 133    FORMAT(' HT,DTCLST ',2F10.1)
      ENDIF
C
C     Compute the local exospheric temperature.
C     Add geomagnetic storm effect from input dTc value

      TEMP(1) = TSUBL + DSTDTC
      TINF    = TSUBL + DSTDTC + DTCLST

C     Equation (9)

      TSUBX = 444.3807D0 + 0.02385D0 * TINF - 392.8292D0 *
     >        DEXP(-0.0021357D0 * TINF)
C
      IF (IDEBUG.EQ.1) THEN
        WRITE(8,142) TINF,TSUBX,CAPPHI
 142    FORMAT(' JB2008 - TINF,TSUBX,CAPPHI',3D20.8)
      ENDIF
C
C     Equation (11)

      GSUBX = 0.054285714D0 * (TSUBX - 183.D0)

C     The TC array will be an argument in the call to
C     XLOCAL, which evaluates Equation (10) or Equation (13)

      TC(1) = TSUBX
      TC(2) = GSUBX

C     A AND GSUBX/A OF Equation (13)

      TC(3) = (TINF - TSUBX)/PIOV2
      TC(4) = GSUBX/TC(3)
C
      IF (IDEBUG.EQ.1) THEN
        WRITE(8,123) TC
 123    FORMAT(' JB2008 - TC',4F20.6)
      ENDIF
C
C     Equation (5)

      Z1 = 90.D0
      Z2 = DMIN1(SAT(3),105.D0)
      AL = DLOG(Z2/Z1)
      N = IDINT(AL/R1) + 1
      ZR = DEXP(AL/DFLOAT(N))
      AMBAR1 = XAMBAR(Z1)
      TLOC1 = XLOCAL(Z1,TC)
      ZEND = Z1
      SUM2 = 0.D0
      AIN = AMBAR1 * XGRAV(Z1)/TLOC1
C
      IF (IDEBUG.EQ.1) THEN
        WRITE(8,163) ZR,AMBAR1,TLOC1,AIN
 163    FORMAT(' JB2008 - ZR-AIN',4D20.8)
      ENDIF
C
      DO 2 I = 1,N
        Z = ZEND
        ZEND = ZR * Z
        DZ = 0.25D0 * (ZEND-Z)
        SUM1 = WT(1)*AIN
        DO 1 J = 2,5
          Z = Z + DZ
          AMBAR2 = XAMBAR(Z)
          TLOC2 = XLOCAL(Z,TC)
          GRAVL = XGRAV(Z)
          AIN = AMBAR2 * GRAVL/TLOC2
C
C     IF (IDEBUG.EQ.1) THEN
C       WRITE(8,168) I,J,SUM1,TLOC2
C168    FORMAT(' JB2008 - SUM1,TLOC2',2I4,2D20.8)
C     ENDIF
C
    1 SUM1 = SUM1 + WT(J) * AIN
    2 SUM2 = SUM2 + DZ * SUM1

      FACT1 = 1000.0/RSTAR
      RHO = 3.46D-6 * AMBAR2 * TLOC1 * DEXP(-FACT1*SUM2) /AMBAR1 /TLOC2

C     Equation (2)

      ANM = AVOGAD * RHO
      AN  = ANM/AMBAR2

C     Equation (3)

      FACT2  = ANM/28.960D0
      ALN(1) = DLOG(FRAC(1)*FACT2)
      ALN(4) = DLOG(FRAC(3)*FACT2)
      ALN(5) = DLOG(FRAC(4)*FACT2)

C     Equation (4)

      ALN(2) = DLOG(FACT2 * (1.D0 + FRAC(2)) - AN)
      ALN(3) = DLOG(2.D0 * (AN - FACT2))
C
      IF (IDEBUG.EQ.1) THEN
        WRITE(8,143) ALN
 143    FORMAT(' JB2008 - ALN',6D20.8)
      ENDIF
C
      IF (SAT(3) .GT. 105.D0) GO TO 3
      TEMP(2) = TLOC2

C     Put in negligible hydrogen for use in DO-LOOP 13

      ALN(6) = ALN(5) - 25.D0
      GO TO 11

C     Equation (6)

    3 Z3 = DMIN1(SAT(3),500.D0)
      AL = DLOG(Z3/Z)
      N = IDINT(AL/R2) + 1
      ZR = DEXP(AL/DFLOAT(N))
      SUM2 = 0.D0
      AIN = GRAVL/TLOC2

      DO 5 I = 1,N
        Z = ZEND
        ZEND = ZR * Z
        DZ = 0.25D0 * (ZEND - Z)
        SUM1 = WT(1) * AIN
        DO 4 J = 2,5
          Z = Z + DZ
          TLOC3 = XLOCAL(Z,TC)
          GRAVL = XGRAV(Z)
          AIN = GRAVL/TLOC3
    4   SUM1 = SUM1 + WT(J) * AIN
    5 SUM2 = SUM2 + DZ * SUM1
C
      IF (IDEBUG.EQ.1) THEN
        WRITE(8,124) SUM2,TLOC1,TLOC2,TLOC3
 124    FORMAT(' JB2008 - SUM2,TLOC..',4D20.6)
      ENDIF
C
      Z4 = DMAX1(SAT(3),500.D0)
      AL = DLOG(Z4/Z)
      R = R2
      IF (SAT(3) .GT. 500.D0) R = R3
      N = IDINT(AL/R) + 1
      ZR = DEXP(AL/DFLOAT(N))
      SUM3 = 0.D0

      DO 7 I=1,N
        Z = ZEND
        ZEND = ZR * Z
        DZ = 0.25D0 * (ZEND - Z)
        SUM1 = WT(1) * AIN
        DO 6 J = 2,5
          Z = Z + DZ
          TLOC4 = XLOCAL(Z,TC)
          GRAVL = XGRAV(Z)
          AIN = GRAVL/TLOC4
    6   SUM1 = SUM1 + WT(J) * AIN
    7 SUM3 = SUM3 + DZ * SUM1
C
      IF(IDEBUG.EQ.1) THEN
        WRITE(8,106) SUM3,TLOC4
  106   FORMAT(' SUM3,TLOC4 ',D20.9,F16.5)
      ENDIF
C
      IF (SAT(3) .GT. 500.D0) GO TO 8
      T500 = TLOC4
      TEMP(2) = TLOC3
      ALTR = DLOG(TLOC3/TLOC2)
      FACT2 = FACT1 * SUM2
      HSIGN = 1.D0
      GO TO 9
    8 T500 = TLOC3
      TEMP(2) = TLOC4
      ALTR = DLOG(TLOC4/TLOC2)
      FACT2 = FACT1 * (SUM2 + SUM3)
      HSIGN = -1.D0

    9 DO 10 I = 1,5
   10 ALN(I) = ALN(I) - (1.0 + ALPHA(I)) * ALTR - FACT2 * AMW(I)
C     Equation (7) - Note that in CIRA72, AL10T5 = DLOG10(T500)
C
      IF(IDEBUG.EQ.1) THEN
        WRITE(8,149) ALTR,FACT2
  149   FORMAT(' JB2008 ALTR,FACT2 ',2D20.8)
      ENDIF
C
C
      AL10T5 = DLOG10(TINF)
      ALNH5 = (5.5D0 * AL10T5 - 39.40D0) * AL10T5 + 73.13D0
      ALN(6) = AL10 * (ALNH5 + 6.D0) + HSIGN * (DLOG(TLOC4/TLOC3)
     >       + FACT1 * SUM3 * AMW(6))
C
      IF(IDEBUG.EQ.1) THEN
        WRITE(8,146) ALN
  146   FORMAT(' JB2008 ALN ',6D20.8)
      ENDIF
C
   11 CONTINUE

C     Equation (24)  - J70 Seasonal-Latitudinal Variation
C
      TRASH = (AMJD - 36204.D0) / 365.2422D0
      CAPPHI = DMOD(TRASH,1.D0)
C
      DLRSL = 0.02D0 * (SAT(3) - 90.D0)
     >      * DEXP(-0.045D0 * (SAT(3) - 90.D0))
     >      * DSIGN(1.D0,SAT(2)) * DSIN(TWOPI * CAPPHI+ 1.72D0)
     >      * DSIN(SAT(2))**2
C
C     Equation (23) - Computes the semiannual variation
C
      DLRSA = 0.D0
      IF (Z.LT.2000.) THEN
        D1950 = AMJD - 33281.D0
        CALL TMOUTD (D1950,IYR,YRDAY)
C       Use new semiannual model
        CALL SEMIAN08 (YRDAY,ZHT,F10B,S10B,XM10B,FZZ,GTZ,DLRSA)
C
        IF (FZZ.LT.0.0D0) DLRSA   = 0.D0
      ENDIF
C
C
      IF(IDEBUG.EQ.1) THEN
        WRITE(8,147) DLRSL, DLRSA
  147   FORMAT(' JB2008 DLRSL',2D20.8)
      ENDIF
C
C
C     Sum the delta-log-rhos and apply to the number densities.
C     In CIRA72 the following equation contains an actual sum,
C     namely DLR = AL10 * (DLRGM + DLRSA + DLRSL)
C     However, for Jacchia 70, there is no DLRGM or DLRSA.

      DLR = AL10 * (DLRSL + DLRSA)

      DO 12 I = 1,6
   12 ALN(I) = ALN(I) + DLR

C     Compute mass-density and mean-molecular-weight and
C     convert number density logs from natural to common.

      SUMN = 0.D0
      SUMNM = 0.D0

      DO 13 I = 1,6
        AN = DEXP(ALN(I))
        SUMN = SUMN + AN
        SUMNM = SUMNM + AN*AMW(I)
        AL10N(I) = ALN(I)/AL10
        IF(IDEBUG.EQ.1) THEN
          XVAL = 10.D0**AL10N(I)
          WRITE (8,100) I,SUMNM,XVAL,ALN(I)
 100      FORMAT(' JB2008 I,SUMN,XVAL,ALN ',I5,3D20.9)
        ENDIF
   13 CONTINUE

      RHO = SUMNM/AVOGAD

C     Compute the high altitude exospheric density correction factor

      FEX = 1.D0
C
      IF ((ZHT.GE.1000.D0).AND.(ZHT.LT.1500.D0)) THEN
        ZETA   = (ZHT - 1000.D0) * 0.002D0
        ZETA2  =  ZETA * ZETA
        ZETA3  =  ZETA * ZETA2
        F15C   = CHT(1) + CHT(2)*F10B + CHT(3)*1500.D0
     *                  + CHT(4)*F10B*1500.D0
        F15C_ZETA = (CHT(3) + CHT(4)*F10B) * 500.D0
        FEX2   = 3.D0 * F15C - F15C_ZETA - 3.D0
        FEX3   = F15C_ZETA - 2.D0 * F15C + 2.D0
        FEX    = 1.D0 + FEX2 * ZETA2 + FEX3 * ZETA3
      ENDIF
      IF (ZHT .GE. 1500.D0) THEN
        FEX    = CHT(1) + CHT(2)*F10B + CHT(3)*ZHT + CHT(4)*F10B*ZHT
      ENDIF

C     Apply the exospheric density correction factor.

      RHO    = FEX * RHO

C
      IF(IDEBUG.EQ.1) THEN
        WRITE (8,160) RHO
 160    FORMAT(' JB2008 RHO ',D20.8,F10.4)
      ENDIF
C
      RETURN
      END

C***********************************************************************

      DOUBLE PRECISION FUNCTION XAMBAR(Z)

C     Evaluates Equation (1)

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION C(7)

      DATA C(1) /28.15204D0/
      DATA C(2) /-8.5586D-2/
      DATA C(3) /+1.2840D-4/
      DATA C(4) /-1.0056D-5/
      DATA C(5) /-1.0210D-5/
      DATA C(6) /+1.5044D-6/
      DATA C(7) /+9.9826D-8/

      DZ = Z - 100.D0
      AMB = C(7)
      DO 1 I = 1,6
      J = 7-I
    1 AMB = DZ * AMB + C(J)
      XAMBAR = AMB
      RETURN
      END

C***********************************************************************

      DOUBLE PRECISION FUNCTION XGRAV(Z)

C     Evaluates Equation (8)

      IMPLICIT REAL*8(A-H,O-Z)

      XGRAV = 9.80665D0/(1.D0 + Z/6356.766D0)**2
      RETURN
      END

C***********************************************************************

      DOUBLE PRECISION FUNCTION XLOCAL(Z,TC)

C     Evaluates Equation (10) or Equation (13), depending on Z

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION TC(4)

      DZ = Z - 125.D0
      IF (DZ .GT. 0.D0) GO TO 1

      XLOCAL = ((-9.8204695D-6 * DZ - 7.3039742D-4) * DZ**2 + 1.0)
     >       * DZ * TC(2) + TC(1)
      RETURN
 1    XLOCAL = TC(1) + TC(3) * DATAN(TC(4)*DZ*(1.D0 + 4.5D-6*DZ**2.5))
      RETURN
      END

C***********************************************************************

      SUBROUTINE DTSUB (F10,XLST,XLAT,ZHT,DTC)
C
C     COMPUTE dTc correction for Jacchia-Bowman model
C
C        Calling Args:
C        ------------
C        F10       = (I)   F10 FLUX
C        XLST      = (I)   LOCAL SOLAR TIME (HOURS 0-23.999)
C        XLAT      = (I)   XLAT = SAT LAT (RAD)
C        ZHT       = (I)   ZHT = HEIGHT (KM)
C        DTC       = (O)   dTc correction
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION B(19),C(23)
C
      DATA B /
C
     1       -0.457512297D+01, -0.512114909D+01, -0.693003609D+02,
     2        0.203716701D+03,  0.703316291D+03, -0.194349234D+04,
     3        0.110651308D+04, -0.174378996D+03,  0.188594601D+04,
     4       -0.709371517D+04,  0.922454523D+04, -0.384508073D+04,
     5       -0.645841789D+01,  0.409703319D+02, -0.482006560D+03,
     6        0.181870931D+04, -0.237389204D+04,  0.996703815D+03,
     7        0.361416936D+02/
C
      DATA C /
C
     1       -0.155986211D+02, -0.512114909D+01, -0.693003609D+02,
     2        0.203716701D+03,  0.703316291D+03, -0.194349234D+04,
     3        0.110651308D+04, -0.220835117D+03,  0.143256989D+04,
     4       -0.318481844D+04,  0.328981513D+04, -0.135332119D+04,
     5        0.199956489D+02, -0.127093998D+02,  0.212825156D+02,
     6       -0.275555432D+01,  0.110234982D+02,  0.148881951D+03,
     7       -0.751640284D+03,  0.637876542D+03,  0.127093998D+02,
     8       -0.212825156D+02,  0.275555432D+01/
C
C
      DTC = 0.D0
C
      tx  = XLST/24.D0
      ycs = DCOS(XLAT)
      F   = (F10 - 100.D0)/100.D0
C
C         calculates dTc
C
      IF (ZHT.GE.120.0.AND.ZHT.LE.200.0) THEN
        H = (ZHT - 200.D0)/50.D0
        DTC200 =
     1    + C(17)             + C(18)*tx*ycs      + C(19)*tx**2*ycs
     2    + C(20)*tx**3*ycs   + C(21)*F*ycs       + C(22)*tx*F*ycs
     3    + C(23)*tx**2*F*ycs
        sum = C(1) + B(2)*F + C(3)*tx*F     + C(4)*tx**2*F
     1    + C(5)*tx**3*F    + C(6)*tx**4*F    + C(7)*tx**5*F
     2    + C(8)*tx*ycs     + C(9)*tx**2*ycs  + C(10)*tx**3*ycs
     3    + C(11)*tx**4*ycs + C(12)*tx**5*ycs + C(13)*ycs
     4    + C(14)*F*ycs     + C(15)*tx*F*ycs  + C(16)*tx**2*F*ycs
        DTC200DZ = SUM
        CC  = 3.*DTC200 - DTC200DZ
        DD  = DTC200 - CC
        ZP  = (ZHT-120.)/80.D0
        DTC = CC*ZP*ZP + DD*ZP*ZP*ZP
      ENDIF
C
C
      IF (ZHT.GT.200.0.AND.ZHT.LE.240.0) THEN
        H = (ZHT - 200.D0)/50.D0
        sum = C(1)*H + B(2)*F*H + C(3)*tx*F*H     + C(4)*tx**2*F*H
     1    + C(5)*tx**3*F*H    + C(6)*tx**4*F*H    + C(7)*tx**5*F*H
     2    + C(8)*tx*ycs*H     + C(9)*tx**2*ycs*H  + C(10)*tx**3*ycs*H
     3    + C(11)*tx**4*ycs*H + C(12)*tx**5*ycs*H + C(13)*ycs*H
     4    + C(14)*F*ycs*H     + C(15)*tx*F*ycs*H  + C(16)*tx**2*F*ycs*H
     5    + C(17)             + C(18)*tx*ycs      + C(19)*tx**2*ycs
     6    + C(20)*tx**3*ycs   + C(21)*F*ycs       + C(22)*tx*F*ycs
     7    + C(23)*tx**2*F*ycs
        DTC = sum
      ENDIF
C
C
      IF (ZHT.GT.240.0.AND.ZHT.LE.300.0) THEN
        H = (40.D0)/50.D0
        sum = C(1)*H + B(2)*F*H + C(3)*tx*F*H     + C(4)*tx**2*F*H
     1    + C(5)*tx**3*F*H    + C(6)*tx**4*F*H    + C(7)*tx**5*F*H
     2    + C(8)*tx*ycs*H     + C(9)*tx**2*ycs*H  + C(10)*tx**3*ycs*H
     3    + C(11)*tx**4*ycs*H + C(12)*tx**5*ycs*H + C(13)*ycs*H
     4    + C(14)*F*ycs*H     + C(15)*tx*F*ycs*H  + C(16)*tx**2*F*ycs*H
     5    + C(17)             + C(18)*tx*ycs      + C(19)*tx**2*ycs
     6    + C(20)*tx**3*ycs   + C(21)*F*ycs       + C(22)*tx*F*ycs
     7    + C(23)*tx**2*F*ycs
        AA = SUM
        BB = C(1) + B(2)*F  + C(3)*tx*F       + C(4)*tx**2*F
     1    + C(5)*tx**3*F    + C(6)*tx**4*F    + C(7)*tx**5*F
     2    + C(8)*tx*ycs     + C(9)*tx**2*ycs  + C(10)*tx**3*ycs
     3    + C(11)*tx**4*ycs + C(12)*tx**5*ycs + C(13)*ycs
     4    + C(14)*F*ycs     + C(15)*tx*F*ycs  + C(16)*tx**2*F*ycs
        H   = 300./100.D0
        sum = B(1)    + B(2)*F  + B(3)*tx*F         + B(4)*tx**2*F
     1      + B(5)*tx**3*F      + B(6)*tx**4*F      + B(7)*tx**5*F
     2      + B(8)*tx*ycs       + B(9)*tx**2*ycs    + B(10)*tx**3*ycs
     3      + B(11)*tx**4*ycs   + B(12)*tx**5*ycs   + B(13)*H*ycs
     4      + B(14)*tx*H*ycs    + B(15)*tx**2*H*ycs + B(16)*tx**3*H*ycs
     5      + B(17)*tx**4*H*ycs + B(18)*tx**5*H*ycs + B(19)*ycs
        DTC300 = SUM
        sum = B(13)*ycs
     1      + B(14)*tx*ycs    + B(15)*tx**2*ycs + B(16)*tx**3*ycs
     2      + B(17)*tx**4*ycs + B(18)*tx**5*ycs
        DTC300DZ = SUM
        CC = 3.*DTC300 - DTC300DZ - 3.*AA - 2.*BB
        DD = DTC300 - AA - BB - CC
        ZP  = (ZHT-240.)/60.D0
        DTC = AA + BB*ZP + CC*ZP*ZP + DD*ZP*ZP*ZP
      ENDIF
C
C
      IF (ZHT.GT.300.0.AND.ZHT.LE.600.0) THEN
        H   = ZHT/100.D0
        sum = B(1)    + B(2)*F  + B(3)*tx*F         + B(4)*tx**2*F
     1      + B(5)*tx**3*F      + B(6)*tx**4*F      + B(7)*tx**5*F
     2      + B(8)*tx*ycs       + B(9)*tx**2*ycs    + B(10)*tx**3*ycs
     3      + B(11)*tx**4*ycs   + B(12)*tx**5*ycs   + B(13)*H*ycs
     4      + B(14)*tx*H*ycs    + B(15)*tx**2*H*ycs + B(16)*tx**3*H*ycs
     5      + B(17)*tx**4*H*ycs + B(18)*tx**5*H*ycs + B(19)*ycs
        DTC = sum
      ENDIF
C
C
      IF (ZHT.GT.600.0.AND.ZHT.LE.800.0) THEN
        ZP = (ZHT - 600.D0)/100.
        HP = 600.D0/100.D0
        AA  = B(1)    + B(2)*F  + B(3)*tx*F         + B(4)*tx**2*F
     1      + B(5)*tx**3*F      + B(6)*tx**4*F      + B(7)*tx**5*F
     2      + B(8)*tx*ycs       + B(9)*tx**2*ycs    + B(10)*tx**3*ycs
     3      + B(11)*tx**4*ycs   + B(12)*tx**5*ycs   + B(13)*HP*ycs
     4      + B(14)*tx*HP*ycs   + B(15)*tx**2*HP*ycs+ B(16)*tx**3*HP*ycs
     5      + B(17)*tx**4*HP*ycs + B(18)*tx**5*HP*ycs + B(19)*ycs
        BB  = B(13)*ycs
     1      + B(14)*tx*ycs    + B(15)*tx**2*ycs + B(16)*tx**3*ycs
     2      + B(17)*tx**4*ycs + B(18)*tx**5*ycs
        CC  = -(3.*AA+4.*BB)/4.
        DD  = (AA+BB)/4.
        DTC = AA + BB*ZP + CC*ZP*ZP + DD*ZP*ZP*ZP
      ENDIF
C
      RETURN
      END SUBROUTINE DTSUB
C
C******************************************************************
C
      SUBROUTINE SEMIAN08 (DAY,HT,F10B,S10B,XM10B,FZZ,GTZ,DRLOG)
C
C     COMPUTE SEMIANNUAL VARIATION (DELTA LOG RHO)
C     INPUT DAY, HEIGHT, F10B, S10B, M10B FSMB
C           025.  650.   150.  148.  147. 151.
C     OUTPUT FUNCTIONS FZ, GT, AND DEL LOG RHO VALUE
C
C     DAY     (I)   DAY OF YEAR
C     HT      (I)   HEIGHT (KM)
C     F10B    (I)   AVE 81-DAY CENTERED F10
C     S10B    (I)   AVE 81-DAY CENTERED S10
C     XM10B   (I)   AVE 81-DAY CENTERED M10
C     FZZ     (O)   SEMIANNUAL AMPLITUDE
C     GTZ     (O)   SEMIANNUAL PHASE FUNCTION
C     DRLOG   (O)   DELTA LOG RHO
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      REAL*8 FZM(5),GTM(10)
C
      DATA TWOPI/6.2831853072D0/
C
C     FZ GLOBAL MODEL VALUES
C     1997-2006 FIT:
      DATA FZM /
     1   0.2689D+00,-0.1176D-01, 0.2782D-01,
     2  -0.2782D-01, 0.3470D-03/
C
C     GT GLOBAL MODEL VALUES
C     1997-2006 FIT:
      DATA GTM /
     1  -0.3633D+00, 0.8506D-01, 0.2401D+00,-0.1897D+00,
     2  -0.2554D+00,-0.1790D-01, 0.5650D-03,-0.6407D-03,
     3  -0.3418D-02,-0.1252D-02/
C
C
C     COMPUTE NEW 81-DAY CENTERED SOLAR INDEX FOR FZ
      FSMB  = 1.00*F10B - 0.70*S10B - 0.04*XM10B
C
      HTZ = HT/1000.D0
C
      FZZ = FZM(1) + FZM(2)*FSMB  + FZM(3)*FSMB*HTZ
     *    + FZM(4)*FSMB*HTZ**2    + FZM(5)*FSMB**2*HTZ
C
C
C     COMPUTE DAILY 81-DAY CENTERED SOLAR INDEX FOR GT
      FSMB  = 1.00*F10B - 0.75*S10B - 0.37*XM10B
C
      TAU   = (DAY-1.D0)/365.D0
      SIN1P = DSIN(TWOPI*TAU)
      COS1P = DCOS(TWOPI*TAU)
      SIN2P = DSIN(2.*TWOPI*TAU)
      COS2P = DCOS(2.*TWOPI*TAU)
C
      GTZ = GTM(1) + GTM(2)*SIN1P + GTM(3)*COS1P
     1             + GTM(4)*SIN2P + GTM(5)*COS2P
     2             + GTM(6)*FSMB
     3             + GTM(7)*FSMB*SIN1P + GTM( 8)*FSMB*COS1P
     4             + GTM(9)*FSMB*SIN2P + GTM(10)*FSMB*COS2P
C
      IF (FZZ.LT.1.D-6) FZZ = 1.D-6
C
      DRLOG = FZZ*GTZ
C
      RETURN
      END SUBROUTINE SEMIAN08
C
C******************************************************************
C
      SUBROUTINE TMOUTD(D1950,IYR,DAY)
C
C     COMPUTE DAY AND YEAR FROM TIME D1950 (DAYS SINCE 1950)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      IYDAY = D1950
      FRAC = D1950 - IYDAY
      IYDAY = IYDAY + 364
      ITEMP = IYDAY/1461
      IYDAY = IYDAY - ITEMP*1461
      IYR = 1949 + 4*ITEMP
      ITEMP = IYDAY/365
      IF (ITEMP.GE.3) ITEMP = 3
      IYR = IYR + ITEMP
      IYDAY = IYDAY - 365*ITEMP + 1
      IYR = IYR - 1900
      DAY = IYDAY + FRAC
      IF (IYR.GE.100) IYR = IYR - 100
C
      RETURN
      END SUBROUTINE TMOUTD
