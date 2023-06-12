#ifndef __sgp4header__
#define __sgp4header__

/**
 * This class implements the elsetrec data type from Vallado's SGP4 code.
 *
 * From SGP4.h
 * #define SGP4Version  "SGP4 Version 2016-03-09"
 *
 * @author aholinch
 *
 */
class ElsetRec
{
  public:
    int whichconst;
    int epochtynumrev;
    int error;
    char operationmode;
    char init;
    char method;
    double a;
    double altp;
    double alta;
    double jdsatepoch;
    double jdsatepochF;
    double nddot;
    double ndot;
    double bstar;
    double rcse;
    double inclo;
    double nodeo;
    double ecco;
    double argpo;
    double mo;
    double no_kozai;

    // sgp4fix add unkozai'd variable
    double no_unkozai;

    // sgp4fix add singly averaged variables
    double am;
    double em;
    double im;
    double Om;
    double om;
    double mm;
    double nm;
    double t;

    // sgp4fix add constant parameters to eliminate mutliple calls during execution
    double tumin;
    double mu;
    double radiusearthkm;
    double xke;
    double j2;
    double j3;
    double j4;
    double j3oj2;

    // Additional elements to capture relevant TLE and object information:
    long dia_mm; // RSO dia in mm
    double period_sec; // Period in seconds
    char active; // "Active S/C" flag (0=n, 1=y)
    char not_orbital; // "Orbiting S/C" flag (0=n, 1=y)
    double rcs_m2; // "RCS (m^2)" storage

    // temporary variables because the original authors call the same method with different variables
    double ep;
    double inclp;
    double nodep;
    double argpp;
    double mp;

    int isimp;
    double aycof;
    double con41;
    double cc1;
    double cc4;
    double cc5;
    double d2;
    double d3;
    double d4;
    double delmo;
    double eta;
    double argpdot;
    double omgcof;
    double sinmao;
    double t2cof;
    double t3cof;
    double t4cof;
    double t5cof;
    double x1mth2;
    double x7thm1;
    double mdot;
    double nodedot;
    double xlcof;
    double xmcof;
    double nodecf;

    // deep space
    int irez;
    double d2201;
    double d2211;
    double d3210;
    double d3222;
    double d4410;
    double d4422;
    double d5220;
    double d5232;
    double d5421;
    double d5433;
    double dedt;
    double del1;
    double del2;
    double del3;
    double didt;
    double dmdt;
    double dnodt;
    double domdt;
    double e3;
    double ee2;
    double peo;
    double pgho;
    double pho;
    double pinco;
    double plo;
    double se2;
    double se3;
    double sgh2;
    double sgh3;
    double sgh4;
    double sh2;
    double sh3;
    double si2;
    double si3;
    double sl2;
    double sl3;
    double sl4;
    double gsto;
    double xfact;
    double xgh2;
    double xgh3;
    double xgh4;
    double xh2;
    double xh3;
    double xi2;
    double xi3;
    double xl2;
    double xl3;
    double xl4;
    double xlamo;
    double zmol;
    double zmos;
    double atime;
    double xli;
    double xni;
    double snodm;
    double cnodm;
    double sinim;
    double cosim;
    double sinomm;
    double cosomm;
    double day;
    double emsq;
    double gam;
    double rtemsq;
    double s1;
    double s2;
    double s3;
    double s4;
    double s5;
    double s6;
    double s7;
    double ss1;
    double ss2;
    double ss3;
    double ss4;
    double ss5;
    double ss6;
    double ss7;
    double sz1;
    double sz2;
    double sz3;
    double sz11;
    double sz12;
    double sz13;
    double sz21;
    double sz22;
    double sz23;
    double sz31;
    double sz32;
    double sz33;
    double z1;
    double z2;
    double z3;
    double z11;
    double z12;
    double z13;
    double z21;
    double z22;
    double z23;
    double z31;
    double z32;
    double z33;
    double argpm;
    double inclm;
    double nodem;
    double dndt;
    double eccsq;

    // for initl
    double ainv;
    double ao;
    double con42;
    double cosio;
    double cosio2;
    double omeosq;
    double posq;
    double rp;
    double rteosq;
    double sinio;
};  // end struct


bool sgp4init(char opsmode, ElsetRec &rec);

bool sgp4(ElsetRec &rec, double tsince, double r[3], double v[3]);

#endif
