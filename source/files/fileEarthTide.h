/***********************************************/
/**
* @file fileEarthTide.h
*
* @brief Read/write EarthTide.
*
* @author Torsten Mayer-Guerr
* @date 2017-12-09
*
*/
/***********************************************/

#ifndef __GROOPS_FILEMEARTHTIDE__
#define __GROOPS_FILEMEARTHTIDE__

#include "base/doodson.h"
#include "inputOutput/fileArchive.h"

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_EarthTide
static const char *docstringEarthTide = R"(
Containing the Love numbers together with frequency corrections to compute
the gravitational potential and the geometric displacements due to solid Earth tides.
It is used by \configClass{tides}{tidesType}.
)";
#endif

/***********************************************/

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_EARTHTIDE_TYPE    = "earthTide";
constexpr UInt    FILE_EARTHTIDE_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** FUNCTIONS *******************************/

/** @brief Read from a EarthTide file. */
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
                       Vector &dT21_ip, Vector &dT21_op, Vector &dT20_ip, Vector &dT20_op);

/** @brief Write an EarthTide file. */
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
                        const_MatrixSliceRef dT21_ip, const_MatrixSliceRef dT21_op, const_MatrixSliceRef dT20_ip, const_MatrixSliceRef dT20_op);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
