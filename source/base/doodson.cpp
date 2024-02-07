/***********************************************/
/**
* @file doodson.cpp
*
* @brief Doodson arguments and multipliers.
*
* @author Torsten Mayer-Guerr
* @author Daniel Rieser
* @date 2005-07-15
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/constants.h"
#include "base/matrix.h"
#include "base/time.h"
#include "base/planets.h"
#include "base/doodson.h"

/***********************************************/

struct DoodsonName
{
  const char *code;
  const char *name;
};

static const DoodsonName doodsonName[] = {
      // long periodic
      {"055.565", "om1"},     // PI
      {"055.575", "om2"},     // 0
      {"056.554", "sa"},      // 0
      {"056.555", "sa"},      // 0
      {"057.555", "ssa"},     // 0
      {"058.554", "sta"},     // 0
      {"063.655", "msm"},     // 0
      {"065.455", "mm"},      // 0
      {"073.555", "msf"},     // 0
      {"075.555", "mf"},      // 0
      {"083.655", "mstm"},    // 0
      {"085.455", "mtm"},     // 0
      {"093.555", "msq"},     // 0
      {"093.555", "msqm"},    // 0 // alternative name
      // diurnal
      {"125.755", "2q1"},     //  PI/2
      {"127.555", "sig1"},    //  PI/2
      {"127.555", "sigma1"},  //  PI/2 // alternative name
      {"135.655", "q1"},      //  PI/2
      {"137.455", "ro1"},     //  PI/2
      {"137.455", "rho1"},    //  PI/2 // alternative name
      {"145.555", "o1"},      //  PI/2
      {"147.555", "tau1"},    // -PI/2
      {"155.655", "m1"},      // -PI/2
      {"157.455", "chi1"},    // -PI/2
      {"162.556", "pi1"},     //  PI/2
      {"163.555", "p1"},      //  PI/2
      {"164.555", "s1"},      //  PI   // !!! sometimes defined as 164.556
      {"164.556", "s1"},      //  PI   // !!! sometimes defined as 164.556
      {"165.555", "k1"},      // -PI/2
      {"166.554", "psi1"},    // -PI/2
      {"167.555", "fi1"},     // -PI/2
      {"167.555", "phi1"},    // -PI/2 // alternative name
      {"173.655", "the1"},    // -PI/2
      {"173.655", "theta1"},  // -PI/2 // alternative name
      {"175.455", "j1"},      // -PI/2
      {"183.555", "so1"},     // -PI/2
      {"185.555", "oo1"},     // -PI/2
      {"195.455", "v1"},      // -PI/2
      // semidiurnal
      {"225.855", "3n2"},     // 0
      {"227.655", "eps2"},    // 0
      {"235.755", "2n2"},     // 0
      {"237.555", "mu2"},     // 0
      {"237.555", "mi2"},     // 0   // alternative name
      {"245.655", "n2"},      // 0
      {"247.455", "nu2"},     // 0
      {"247.455", "ni2"},     // 0   // alternative name
      {"253.755", "gam2"},    // PI
      {"254.556", "alf2"},    // PI
      {"255.555", "m2"},      // 0
      {"256.554", "bet2"},    // 0
      {"257.555", "dlt2"},    // 0
      {"263.655", "la2"},     // PI
      {"263.655", "lmb2"},    // PI  // alternative name
      {"263.655", "lambda2"}, // PI  // alternative name
      {"265.455", "l2"},      // PI
      {"271.557", "2t2"},     // 0
      {"272.556", "t2"},      // 0
      {"273.555", "s2"},      // 0
      {"274.554", "r2"},      // PI
      {"275.555", "k2"},      // 0
      {"283.655", "ksi2"},    // 0
      {"285.455", "eta2"},    // 0
      // nonlinear
      {"355.555", "m3"},      // PI
      {"381.555", "t3"},      // PI
      {"382.555", "s3"},      // PI
      {"383.555", "r3"},      // PI
      {"435.755", "n4"},      // 0
      {"445.655", "mn4"},     // 0
      {"455.555", "m4"},      // 0
      {"473.555", "ms4"},     // 0
      {"491.555", "s4"},      // 0
      {"5a0.555", "s5"},      // 0
      {"655.555", "m6"},      // 0
      {"6bz.555", "s6"},      // 0
      {"855.555", "m8"},      // 0

      {nullptr,   nullptr}};

/***** CLASS ***********************************/

Doodson::Doodson(const std::vector<Int> &v)
{
  try
  {
    for(UInt i=0; i<6; i++)
      d[i] = v.at(i);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Doodson::Doodson(const std::string &str)
{
  try
  {
    // convert to upper case
    std::string str2 = str;
    std::transform(str2.begin(), str2.end(), str2.begin(), ::tolower);
    // find a name?
    for(UInt i=0; doodsonName[i].name!=nullptr; i++)
      if(str2 == doodsonName[i].name)
      {
        str2 = doodsonName[i].code;
        break;
      }

    auto doodNumber = [](char c)
    {
      if(std::isdigit(c))          return static_cast<Int>(c-'0');
      if(('a' <= c) && (c <= 'm')) return static_cast<Int>(c-'a')+10;
      if(('n' <= c) && (c <= 'z')) return static_cast<Int>(c-'n')-13;
      throw(Exception("unknown character: "s+c));
    };

    d[0] = doodNumber(str2.at(0));
    d[1] = doodNumber(str2.at(1));
    d[2] = doodNumber(str2.at(2));
    if(str2.at(3) != '.') throw(Exception("dot expected"));
    d[3] = doodNumber(str2.at(4));
    d[4] = doodNumber(str2.at(5));
    d[5] = doodNumber(str2.at(6));
    for(UInt i=1; i<6; i++)
      d[i] -= 5;

    if(str2 != code())
      throw(Exception("something strange: '"+str2+"' != "+code()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("not a Doodson name or number: '"+str+"'", e)
  }
}

/***********************************************/

std::string Doodson::code() const
{
  std::stringstream ss;

  auto doodNumber = [](Int d)
  {
    if(d < 0) return static_cast<char>('n'+(d+13));
    if(d > 9) return static_cast<char>('a'+(d-10));
    return static_cast<char>('0'+d);
  };

  ss<<doodNumber(d[0]);
  ss<<doodNumber(d[1]+5);
  ss<<doodNumber(d[2]+5);
  ss<<".";
  ss<<doodNumber(d[3]+5);
  ss<<doodNumber(d[4]+5);
  ss<<doodNumber(d[5]+5);
  return ss.str();
}

/***********************************************/

std::string Doodson::name() const
{
  std::string c = code();
  for(UInt i=0; doodsonName[i].name!=nullptr; i++)
    if(c == doodsonName[i].code)
      return doodsonName[i].name;
  return code();
}

/***********************************************/

Double Doodson::thetaf(const Time &timeGPS) const
{
  Vector a = arguments(timeGPS);
  return a(0)*d[0]+a(1)*d[1]+a(2)*d[2]+a(3)*d[3]+a(4)*d[4]+a(5)*d[5];
}

/***********************************************/

Double Doodson::frequency(Time timeGPS) const
{
  Vector f(6);
  Double t = timeGPS2JC(timeGPS)/10.; // centuries

  f(0) = (127037328.88553056 + (2* 0.17696111 + (3*-0.00183140 + 4* 0.00008824*t)*t)*t) * DEG2RAD;
  f(1) = (  4812678.81195750 + (2*-0.14663889 + (3* 0.00185140 + 4*-0.00015355*t)*t)*t) * DEG2RAD;
  f(2) = (   360007.69748806 + (2* 0.03032222 + (3* 0.00002000 + 4*-0.00006532*t)*t)*t) * DEG2RAD;
  f(3) = (    40690.13635250 + (2*-1.03217222 + (3*-0.01249168 + 4* 0.00052655*t)*t)*t) * DEG2RAD;
  f(4) = (    19341.36261972 + (2*-0.20756111 + (3*-0.00213942 + 4* 0.00016501*t)*t)*t) * DEG2RAD;
  f(5) = (       17.19457667 + (2* 0.04568889 + (3*-0.00001776 + 4*-0.00003323*t)*t)*t) * DEG2RAD;

  return (f(0)*d[0]+f(1)*d[1]+f(2)*d[2]+f(3)*d[3]+f(4)*d[4]+f(5)*d[5])/365250.;
}

/***********************************************/

Vector Doodson::arguments(const Time &timeGPS)
{
  Vector a(6);
  const Vector f = Planets::fundamentals(timeGPS);
  a(1) = f(2)+f(4);
  a(2) = f(2)+f(4)-f(3);
  a(3) = f(2)+f(4)-f(0);
  a(4) = -f(4);
  a(5) = f(2)+f(4)-f(3)-f(1);
  a(0) = Planets::gmst(timeGPS2UTC(timeGPS)) + PI - a(1);
  return a;
}

/***********************************************/

Matrix Doodson::matrix(const std::vector<Doodson> &doodson)
{
  Matrix A(doodson.size(),6);
  for(UInt i=0; i<doodson.size(); i++)
  {
    A(i,0) = doodson.at(i).d[0];
    A(i,1) = doodson.at(i).d[1];
    A(i,2) = doodson.at(i).d[2];
    A(i,3) = doodson.at(i).d[3];
    A(i,4) = doodson.at(i).d[4];
    A(i,5) = doodson.at(i).d[5];
  }
  return A;
}

/***********************************************/

Matrix Doodson::nodeCorr(const std::vector<Doodson> &doodson, const Time &timeGPS, const UInt &nCorr)
{
  try
  {
    Vector f = arguments(timeGPS);    //Doodson arguments
    Matrix fu(doodson.size(),2);      //nodal factor and phase corrections

    for(UInt i=0; i<6;i++) // bring f between 0-2PI
    {
      f[i] = fmod(f[i], 2.*PI);
      if(f[i] < 0)
        f[i] += 2.*PI;
    };

    Double N = -f(4)+2.0*PI;  // longitude of moon's ascending node in the range of 0 to +2pi
    Double p = f(3);          // longitude of moon's perigee

    //according to IHO
    if (nCorr==1)
    {
      for(UInt i=0; i<doodson.size(); i++)
      {
        //compute nodal factors and phases for constituents given in 'doodson'
        if (doodson.at(i).code() == "065.455" || doodson.at(i).code() == "085.455") //Mm or Mtm
        {
          fu(i,0) = 1.0 - 0.1311*cos(N) + 0.0538*cos(2.0*p) + 0.0505*cos(2.0*p-N); //factor
          fu(i,1) = 0;                                                       //phase
        }
        else if (doodson.at(i).code() == "075.555") //Mf
        {
          fu(i,0) = 1.084 + 0.415*cos(N) + 0.039*cos(2.0*N);
          fu(i,1) = (-23.7*sin(N) + 2.7*sin(2*N) - 0.4*sin(3.0*N))*DEG2RAD;
        }
        else if (doodson.at(i).code() == "093.555") //MSQ
        {
          fu(i,0) = 1.0007 - 0.0373*cos(N) + 0.0002*cos(2.0*N);
          fu(i,1) = (2.14*sin(N))*DEG2RAD;
        }
        else if (doodson.at(i).code() == "145.555" || doodson.at(i).code() == "135.655" ) //O1 or Q1
        {
          fu(i,0) = 1.0176 + 0.1871*cos(N) - 0.0147*cos(2*N);
          fu(i,1) = (10.80*sin(N) - 1.34*sin(2.0*N) + 0.19*sin(3.0*N))*DEG2RAD;
        }
        else if (doodson.at(i).code() == "165.555")  //K1
        {
          fu(i,0) = 1.0060 + 0.1150*cos(N) - 0.0088*cos(2.0*N) + 0.0006*cos(3.0*N);
          fu(i,1) = (-8.86*sin(N)+0.68*sin(2.0*N)-0.07*sin(3.0*N))*DEG2RAD;
        }
        else if (doodson.at(i).code() == "255.555" || doodson.at(i).code() == "245.655" || doodson.at(i).code() == "235.755") //M2 or N2 or 2N2
        {
          fu(i,0) = 1.0007 - 0.0373*cos(N) + 0.0002*cos(2.0*N);
          fu(i,1) = (-2.14*sin(N))*DEG2RAD;
        }
        else if (doodson.at(i).code() == "275.555") //K2
        {
          fu(i,0) = 1.0246 + 0.2863*cos(N) + 0.0083*cos(2.0*N) - 0.0015*cos(3.0*N);
          fu(i,1) = (-17.74*sin(N) + 0.68*sin(2.0*N) - 0.04*sin(3.0*N))*DEG2RAD;
        }
        else if (doodson.at(i).code() == "455.555") //M4 (=2 x M2)
        {
          fu(i,0) = pow(1.0007 - 0.0373*cos(N) + 0.0002*cos(2.0*N),2);
          fu(i,1) = 2.0*(-2.14*sin(N))*DEG2RAD;
        }
        else  //for all other constituents
        {
          fu(i,0) = 1.0;
          fu(i,1) = 0.0;
        }
      } //end for

      return fu;
    } // if(nCorr==1)

    // according to Schureman
    if(nCorr==2)
    {
      Double i         = (5.0 + 8.0/60.0 + 43.3546/3600.0)*DEG2RAD;              // inclination of moon's orbit to ecliptic
      Double omega     = (23.0 + 27.0/60.0 + 8.26/3600.0)*DEG2RAD;               // obliquity of the ecliptic at epoch 1/1/1900
      Double I         = acos(cos(omega)*cos(i) - sin(omega)*sin(i)*cos(N));     // I... inclination of moon's orbit to celestial equator
      Double nu        = asin(sin(i) * sin(N) / sin(I));                         // nu
      Double nu_s      = atan(sin(2.0*I)*sin(nu)/(sin(2*I)*cos(nu)+0.3347));     // nu', Schureman f 224
      Double nu_2s     = atan(pow(sin(I),2)*sin(2.0*nu)/(pow(sin(I),2)*cos(2.0*nu)+0.0727))/2.0; // nü'', Schureman f 232
      Double omega_cap = acos(cos(N)*cos(nu)+sin(N)*sin(nu)*cos(omega));         // capital omega
      Double xi = 0.0;
      if(N>PI)
        xi = N + omega_cap -2.0*PI;
      else
        xi = N - omega_cap;

      for(UInt k=0; k<doodson.size(); k++)
      {
        //compute nodal factors and phases for constituents given in 'doodson'
        if (doodson.at(k).code() == "065.455" || doodson.at(k).code() == "085.455") //Mm or Mtm
        {
          fu(k,0) = (2.0/3.0 - pow(sin(I),2))/0.5021; //factor
          fu(k,1) = 0.0;                                                       //phase
        }
        else if (doodson.at(k).code() == "075.555") //Mf
        {
          fu(k,0) = pow(sin(I),2)/0.1578;
          fu(k,1) = -2.0*xi;
        }
        else if (doodson.at(k).code() == "093.555") //MSQ
        {
          fu(k,0) = pow(sin(I),2)/0.1578;
          fu(k,1) = -2.0*xi;
        }
        else if (doodson.at(k).code() == "145.555" || doodson.at(k).code() == "135.655" ) //O1 or Q1
        {
          fu(k,0) = sin(I)*pow(cos(I/2.0),2)/0.38;
          fu(k,1) = 2.0*xi-nu;
        }
        else if (doodson.at(k).code() == "165.555")  //K1
        {
          fu(k,0) = pow(0.8965*pow(sin(2.0*I),2)+0.6001*sin(2.0*I)*cos(nu)+0.1006,0.5);
          fu(k,1) = -nu_s;
        }
        else if (doodson.at(k).code() == "255.555" || doodson.at(k).code() == "245.655" || doodson.at(k).code() == "235.755") //M2 or N2 or 2N2
        {
          fu(k,0) = pow(cos(I/2.0),4)/0.9154;
          fu(k,1) = 2.0*xi-2.0*nu;
        }
        else if (doodson.at(k).code() == "275.555") //K2
        {
          fu(k,0) = pow(19.0444*pow(sin(I),4)+2.7702*pow(sin(I),2)*cos(2.0*nu)+0.0981,0.5);
          fu(k,1) = -2.0*nu_2s;
        }
        else if (doodson.at(k).code() == "455.555") //M4 (=2 x M2)
        {
          fu(k,0) = pow(cos(I/2.0),4)/0.9154 * pow(cos(I/2.0),4)/0.9154;
          fu(k,1) = 2.0*(2.0*xi-2.0*nu);
        }
        else  // for all other constituents
        {
          fu(k,0)=1.0;
          fu(k,1)=0.0;
        }
      } //end for

      return fu;
    } // if(nCorr==2)

    // no nodal corrections
    for(UInt i=0; i<doodson.size(); i++)
    {
      fu(i,0)=1.0;
      fu(i,1)=0.0;
    }
    return fu;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/
