/***********************************************/
/**
* @file fileEarthTide.cpp
*
* @brief Read/write EarthTide.
*
* @author Torsten Mayer-Guerr
* @date 2017-12-09
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_EarthTide

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileEarthTide.h"

GROOPS_REGISTER_FILEFORMAT(EarthTide, FILE_EARTHTIDE_TYPE)

/***********************************************/

void readFileEarthTide(const FileName &fileName,
                       Matrix &kReal, Matrix &kImag, Matrix &kPlus,
                       Matrix &doodson20, Matrix &doodson21, Matrix &doodson22,
                       Vector &ampIp20, Vector &ampOp20, Vector &ampIp21, Vector &ampOp21, Vector &amp22,
                       Double &h2_0, Double &h2_2,
                       Double &l2_0, Double &l2_2, Double &l21_1, Double &l22_1,
                       Double &h21_imag, Double &l21_imag, Double &h22_imag, Double &l22_imag,
                       Double &h3, Double &l3,
                       Matrix &deformationArg21, Matrix &deformationArg20,
                       Vector &dR21_ip, Vector &dR21_op, Vector &dR20_ip, Vector &dR20_op,
                       Vector &dT21_ip, Vector &dT21_op, Vector &dT20_ip, Vector &dT20_op)
{
  try
  {
    UInt count;
    std::vector<Doodson> doodson;

    InFileArchive file(fileName, FILE_EARTHTIDE_TYPE, FILE_EARTHTIDE_VERSION);

    file>>beginGroup("loveNumber");
    file>>nameValue("kReal", kReal);
    file>>nameValue("kImag", kImag);
    file>>nameValue("kPlus", kPlus);
    file>>endGroup("loveNumber");

    file>>beginGroup("longPeriodic");
    file>>nameValue("count", count);
    doodson.resize(count);
    ampIp20   = Vector(count);
    ampOp20   = Vector(count);
    for(UInt i=0; i<count; i++)
    {
      file>>beginGroup("constituent");
      file>>nameValue("doodson", doodson.at(i));
      file>>nameValue("ampInPhase",  ampIp20(i));
      file>>nameValue("ampOutPhase", ampOp20(i));
      file>>endGroup("constituent");
    }
    doodson20 = Doodson::matrix(doodson);
    file>>endGroup("longPeriodic");

    file>>beginGroup("diurnal");
    file>>nameValue("count", count);
    doodson.resize(count);
    ampIp21   = Vector(count);
    ampOp21   = Vector(count);
    for(UInt i=0; i<count; i++)
    {
      file>>beginGroup("constituent");
      file>>nameValue("doodson", doodson.at(i));
      file>>nameValue("ampInPhase",  ampIp21(i));
      file>>nameValue("ampOutPhase", ampOp21(i));
      file>>endGroup("constituent");
    }
    doodson21 = Doodson::matrix(doodson);
    file>>endGroup("diurnal");

    file>>beginGroup("semidiurnal");
    file>>nameValue("count", count);
    doodson.resize(count);
    amp22 = Vector(count);
    for(UInt i=0; i<count; i++)
    {
      file>>beginGroup("constituent");
      file>>nameValue("doodson", doodson.at(i));
      file>>nameValue("ampInPhase",  amp22(i));
      file>>endGroup("constituent");
    }
    doodson22 = Doodson::matrix(doodson);
    file>>endGroup("semidiurnal");

    // read displacement values
    file>>beginGroup("displacement");
    file>>nameValue("h2_0",     h2_0);
    file>>nameValue("h2_2",     h2_2);
    file>>nameValue("l2_0",     l2_0);
    file>>nameValue("l2_2",     l2_2);
    file>>nameValue("h3",       h3);
    file>>nameValue("l3",       l3);
    file>>nameValue("h21_imag", h21_imag);
    file>>nameValue("l21_imag", l21_imag);
    file>>nameValue("h22_imag", h22_imag);
    file>>nameValue("l22_imag", l22_imag);
    file>>nameValue("l21_1",    l21_1);
    file>>nameValue("l22_1",    l22_1);

    file>>beginGroup("longPeriodic");
    file>>nameValue("count", count);
    doodson.resize(count);
    dR20_ip = Vector(count);
    dR20_op = Vector(count);
    dT20_ip = Vector(count);
    dT20_op = Vector(count);
    for(UInt i=0; i<count; i++)
    {
      file>>beginGroup("constituent");
      file>>nameValue("doodson", doodson.at(i));
      file>>nameValue("dRInPhase",  dR20_ip(i));
      file>>nameValue("dROutPhase", dR20_op(i));
      file>>nameValue("dTInPhase",  dT20_ip(i));
      file>>nameValue("dTOutPhase", dT20_op(i));
      file>>endGroup("constituent");
    }
    deformationArg20 = Doodson::matrix(doodson);
    file>>endGroup("longPeriodic");

    file>>beginGroup("diurnal");
    file>>nameValue("count", count);
    doodson.resize(count);
    dR21_ip = Vector(count);
    dR21_op = Vector(count);
    dT21_ip = Vector(count);
    dT21_op = Vector(count);
    for(UInt i=0; i<count; i++)
    {
      file>>beginGroup("constituent");
      file>>nameValue("doodson", doodson.at(i));
      file>>nameValue("dRInPhase",  dR21_ip(i));
      file>>nameValue("dROutPhase", dR21_op(i));
      file>>nameValue("dTInPhase",  dT21_ip(i));
      file>>nameValue("dTOutPhase", dT21_op(i));
      file>>endGroup("constituent");
    }
    deformationArg21 = Doodson::matrix(doodson);
    file>>endGroup("diurnal");

    file>>endGroup("displacement");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFileEarthTide(const FileName &fileName,
                        const_MatrixSliceRef kReal, const_MatrixSliceRef kImag, const_MatrixSliceRef kPlus,
                        const_MatrixSliceRef doodson20, const_MatrixSliceRef doodson21, const_MatrixSliceRef doodson22,
                        const_MatrixSliceRef ampIp20, const_MatrixSliceRef ampOp20, const_MatrixSliceRef ampIp21, const_MatrixSliceRef ampOp21, const_MatrixSliceRef amp22,
                        Double h2_0, Double h2_2,
                        Double l2_0, Double l2_2, Double l21_1, Double l22_1,
                        Double h21_imag, Double l21_imag, Double h22_imag, Double l22_imag,
                        Double h3, Double l3,
                        const_MatrixSliceRef deformationArg21, const_MatrixSliceRef deformationArg20,
                        const_MatrixSliceRef dR21_ip, const_MatrixSliceRef dR21_op, const_MatrixSliceRef dR20_ip, const_MatrixSliceRef dR20_op,
                        const_MatrixSliceRef dT21_ip, const_MatrixSliceRef dT21_op, const_MatrixSliceRef dT20_ip, const_MatrixSliceRef dT20_op)
{
  try
  {
    OutFileArchive file(fileName, FILE_EARTHTIDE_TYPE, FILE_EARTHTIDE_VERSION);

    file<<beginGroup("loveNumber");
    file<<nameValue("kReal", kReal);
    file<<nameValue("kImag", kImag);
    file<<nameValue("kPlus", kPlus);
    file<<endGroup("loveNumber");

    file<<beginGroup("longPeriodic");
    file<<nameValue("count", doodson20.rows());
    for(UInt i=0; i<doodson20.rows(); i++)
    {
      Doodson doodson;
      for(UInt k=0; k<6; k++)
        doodson.d[k] = static_cast<Int>(doodson20(i, k));
      file<<beginGroup("constituent");
      file<<nameValue("doodson",     doodson);
      file<<nameValue("ampInPhase",  ampIp20(i,0));
      file<<nameValue("ampOutPhase", ampOp20(i,0));
      file<<endGroup("constituent");
    }
    file<<endGroup("longPeriodic");

    file<<beginGroup("diurnal");
    file<<nameValue("count", doodson21.rows());
    for(UInt i=0; i<doodson21.rows(); i++)
    {
      Doodson doodson;
      for(UInt k=0; k<6; k++)
        doodson.d[k] = static_cast<Int>(doodson21(i, k));
      file<<beginGroup("constituent");
      file<<nameValue("doodson",     doodson);
      file<<nameValue("ampInPhase",  ampIp21(i,0));
      file<<nameValue("ampOutPhase", ampOp21(i,0));
      file<<endGroup("constituent");
    }
    file<<endGroup("diurnal");

    file<<beginGroup("semidiurnal");
    file<<nameValue("count", doodson22.rows());
    for(UInt i=0; i<doodson22.rows(); i++)
    {
      Doodson doodson;
      for(UInt k=0; k<6; k++)
        doodson.d[k] = static_cast<Int>(doodson22(i, k));
      file<<beginGroup("constituent");
      file<<nameValue("doodson",    doodson);
      file<<nameValue("ampInPhase", amp22(i,0));
      file<<endGroup("constituent");
    }
    file<<endGroup("semidiurnal");

    // read displacement values
    file<<beginGroup("displacement");
    file<<nameValue("h2_0",     h2_0);
    file<<nameValue("h2_2",     h2_2);
    file<<nameValue("l2_0",     l2_0);
    file<<nameValue("l2_2",     l2_2);
    file<<nameValue("h3",       h3);
    file<<nameValue("l3",       l3);
    file<<nameValue("h21_imag", h21_imag);
    file<<nameValue("l21_imag", l21_imag);
    file<<nameValue("h22_imag", h22_imag);
    file<<nameValue("l22_imag", l22_imag);
    file<<nameValue("l21_1",    l21_1);
    file<<nameValue("l22_1",    l22_1);

    file<<beginGroup("longPeriodic");
    file<<nameValue("count", deformationArg20.rows());
    for(UInt i=0; i<deformationArg20.rows(); i++)
    {
      Doodson doodson;
      for(UInt k=0; k<6; k++)
        doodson.d[k] = static_cast<Int>(deformationArg20(i, k));
      file<<beginGroup("constituent");
      file<<nameValue("doodson",    doodson);
      file<<nameValue("dRInPhase",  dR20_ip(i,0));
      file<<nameValue("dROutPhase", dR20_op(i,0));
      file<<nameValue("dTInPhase",  dT20_ip(i,0));
      file<<nameValue("dTOutPhase", dT20_op(i,0));
      file<<endGroup("constituent");
    }
    file<<endGroup("longPeriodic");

    file<<beginGroup("diurnal");
    file<<nameValue("count", deformationArg21.rows());
    for(UInt i=0; i<deformationArg21.rows(); i++)
    {
      Doodson doodson;
      for(UInt k=0; k<6; k++)
        doodson.d[k] = static_cast<Int>(deformationArg21(i, k));
      file<<beginGroup("constituent");
      file<<nameValue("doodson",    doodson);
      file<<nameValue("dRInPhase",  dR21_ip(i,0));
      file<<nameValue("dROutPhase", dR21_op(i,0));
      file<<nameValue("dTInPhase",  dT21_ip(i,0));
      file<<nameValue("dTOutPhase", dT21_op(i,0));
      file<<endGroup("constituent");
    }
    file<<endGroup("diurnal");

    file<<endGroup("displacement");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
